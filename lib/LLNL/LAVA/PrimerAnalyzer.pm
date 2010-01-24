################################################################################
#
# Copyright (c) 2010, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory
# Written by Clinton Torres <clinton.torres@llnl.gov>.
# CODE-42036.
# All rights reserved.
#
# This file is part of LAVA (LAMP Assay Versatile Analysis). For details, 
# see http://code.google.com/p/lava-dna/ . 
# Please also read the Additional BSD Notice.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# . Redistributions of source code must retain the above copyright notice, 
#   this list of conditions and the disclaimer below.
# . Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the disclaimer (as noted below) in the 
#   documentation and/or other materials provided with the distribution.
# . Neither the name of the LLNS/LLNL nor the names of its contributors may be 
#   used to endorse or promote products derived from this software without 
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, 
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Additional BSD Notice
# 1. This notice is required to be provided under our contract with the 
#    U.S. Department of Energy (DOE). This work was produced at Lawrence 
#    Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with 
#    the DOE.
# 2. Neither the United States Government nor Lawrence Livermore National 
#    Security, LLC nor any of their employees, makes any warranty, express or 
#    implied, or assumes any liability or responsibility for the accuracy, 
#    completeness, or usefulness of any information, apparatus, product, or 
#    process disclosed, or represents that its use would not infringe 
#    privately-owned rights.
# 3. Also, reference herein to any specific commercial products, process, or 
#    services by trade name, trademark, manufacturer or otherwise does not 
#    necessarily constitute or imply its endorsement, recommendation, or 
#    favoring by the United States Government or Lawrence Livermore National 
#    Security, LLC. The views and opinions of authors expressed herein do not 
#    necessarily state or reflect those of the United States Government or 
#    Lawrence Livermore National Security, LLC, and shall not be used for 
#    advertising or product endorsement purposes.
################################################################################

package LLNL::LAVA::PrimerAnalyzer;

use strict;
use vars qw(@ISA);
use Carp;

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::Oligo; # Subclasses recieve these as parameters
use LLNL::LAVA::PrimerInfo; # Subclasses are esponsible for returning these

use LLNL::LAVA::TagHolder; # is-a
@ISA = ("LLNL::LAVA::TagHolder");

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::PrimerAnalyzer - Abstract base class that defines the interface
different PrimerAnalyzers must implement

=head1 SYNOPSIS

do not use PrimerAnalyzer directly, use a concrete sub-class.

  # Get the analysis results for an oligo.
  $primerInfo = $analyzer->analyze($oligo);

  # Set parameter values.  Causes an error if a value is set for
  # a parameter name that doesn't already exist.  This means any parameters
  # for this PrimerAnalyzer need to be declared during instantiation.
  $analyzer->setParameters(
    {
      "parameter_name" => $value,
      ...
    });

  # Retrive a parameter value.  Causes an error if the value for an unknown
  # parameter is requested.
  $value = $analyzer->getParameterValue("parameter_name");

  # Tag use
  $analyzer->setTag("useful_tag_name", $value);
  $exists = $analyzer->getTagExists("useful_tag_name");
  $value = $analyzer->getTag("useful_tag_name");

=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

PrimerAnalyzer is an abstract base class which defines the interface for
concrete PrimerAnalyzer implementations (such as PCRPrimer).

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $analyzer = LLNL::LAVA::PrimerAnalyzer->new();
 Function  : Creates a new LLNL::LAVA::PrimerAnalyzer, only sub-classes 
             should use this
 Arguments : <n/a>
 Example   : See Usage
 Returns   : A new LLNL::LAVA::PrimerAnalyzer

=cut

sub new
{
  my ($classType) = @_;

  # Initialize via TagHolder
  my $this = $classType->SUPER::new();

  # An empty default set of parameters, which should probably be
  # replaced by sub-classes
  $this->{"d_analyzerParameters"} = {};

  return $this;
}

#-------------------------------------------------------------------------------

=head2 analyze

 Usage     : Do not use this directly 
 Function  : Sub-class implementation must return a LLNL::LAVA::PrimerInfo
 Arguments : L<LLNL::LAVA::Oligo> - subclass should accept an individual 
                                    oligo to analyze as a primer
 Example   : See Usage
 Returns   : L<LLNL::LAVA::PrimerInfo> - the result of analyzing the primer

=cut

sub analyze
{
  my ($this) = @_;

  confess("failure - function analyze() not implemented for class type \"" . 
    ref($this) . "\"");

}

#-------------------------------------------------------------------------------

=head2 setParameters

 Usage     : $analyzer->setParameters({param_name => $param_value...});
 Function  : Sets the analyzer parameters for analysis.  Will cause an error
             if a parameter is used that this module doesn't have a default 
             value for.
 Arguments : Hash ref of name->value pairs
 Example   : See Usage
 Returns   : <n/a>

=cut

sub setParameters
{
  my ($this, $paramHash_r) = @_;

  if(!defined $paramHash_r)
  {
    confess("programming error - first parameter is a required hash ref");
  }
  if(ref($paramHash_r) ne "HASH")
  {
    confess("programming error - first parameter must be a hash reference");
  } 

  my $analyzerParam_r = $this->{"d_analyzerParameters"};

  # for each param passed in
  foreach my $currTarget(keys(%{$paramHash_r}))
  {
    # oops if it doesn't exist
    if(! exists $analyzerParam_r->{$currTarget})
    {
      confess("programming error - failing to set analysis parameter " .
        "\"$currTarget\" because no default value for that parameter was " .
        "created during PrimerAnalyzer's instantiation");
    }

    # Since it does exist, set the correct tag in Targets to be the new value
    $analyzerParam_r->{$currTarget} = $paramHash_r->{$currTarget};
  }
}

#-------------------------------------------------------------------------------

=head2 getParameterValue

 Usage     : $paramValue = $analyzer->getParameterValue("parameter_name");
 Function  : Returns the value of a named PrimerAnalyzer parameter.
             Will cause an error (using the Options module) if a parameter 
             name is requested that there is no value for.
 Arguments : String - parameter name to get the value of
 Example   : See Usage
 Returns   : String - value of the named parameter

=cut

sub getParameterValue
{
  my ($this, $paramName) = @_;

  if(! defined $paramName)
  {
    confess("programming error - first parameter is a required parameter name");
  }

  # Could build this to be able to return arrays with a wantarray() check
  # also, but it doesn't look strictly necessary on this first attempt
  # at writing

  my $analyzerParam_r = $this->{"d_analyzerParameters"};

  # Using the Options module to handle strict enforcement of param existence
  my $errorMessage = "programming error - failing to retrieve value of " .
    "analysis parameter \"$paramName\" because that parameter was not " .
    "created during PrimerAnalyzer instantiation";

  my $paramValue = optionRequired($analyzerParam_r, $paramName, $errorMessage);

  # Since it does exist, spit it back out
  return $paramValue
}


1; # Got one?

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
