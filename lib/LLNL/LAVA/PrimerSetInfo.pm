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

package LLNL::LAVA::PrimerSetInfo;

use strict;
use vars qw(@ISA);
use Carp;

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::TagHolder; # is-a
@ISA = ("LLNL::LAVA::TagHolder");

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::PrimerSetInfo - The basic functionality of PrimerSetInfo, which 
represents an analysis result of a PrimerSet within a specific context.

=head1 SYNOPSIS

  use LLNL::LAVA::PrimerSetInfo;

  # Instantiation.  Usually, this is performed by a PrimerAnalyzer so you 
  # won't need to instantiate your own.
  # Instantiation will also usually be of a more specific sub-class such
  # as LLNL::LAVA::PrimerSetInfo::PCRPair
  my $primerInfo = LLNL::LAVA::PrimerSetInfo->new(
    {
      "penalty" => $penalty,
    });"

  # All non-tag functions are read-only because a PrimerSetInfo represents
  # an analysis result within a context, and I want to prevent accidents.
  $penatly = $primerInfo->getPenalty();

  # Tag use
  $primerInfo->setTag("useful_tag_name", $value);
  $exists = $primerInfo->getTagExists("useful_tag_name");
  $value = $primerInfo->getTag("useful_tag_name");
   
=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

PrimerSetInfo is the vehicle for primer set analysis results.  

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $primerInfo = LLNL::LAVA::PrimerSetInfo->new(
               {
                 "penalty" => $penalty,
               });"
 Function  : Creates a new LLNL::LAVA::PrimerSetInfo
 Arguments : Hash ref containing:
               penalty - the PrimerSet's penalty score from an analysis context
 Example   : See Usage
 Returns   : A new LLNL::LAVA::PrimerSetInfo

=cut

sub new
{
  my ($classType, $paramHash_r) = @_;

  if(!defined $paramHash_r)
  {
    confess("programming error - first parameter is a required hash ref");
  }
  if(ref($paramHash_r) ne "HASH")
  {
    confess("programming error - first parameter must be a hash reference");
  } 

  # Read the properties
  my $penalty = optionRequired($paramHash_r, "penalty");

  # Initialize via TagHolder
  my $this = $classType->SUPER::new({});

  # Handle component properties
  $this->{"d_penalty"} = $penalty;

  return $this;
}

#-------------------------------------------------------------------------------

=head2 getPenalty

 Usage     : $penalty = $primerInfo->getPenalty();
 Function  : Gets the penalty for this PrimerSetInfo, built from a specific
             analysis context.
 Arguments : <n/a> 
 Example   : See Usage
 Returns   : Float - in the simplest case, this can be the Primer3 penalty

=cut

sub getPenalty
{
  my ($this) = @_;

  return $this->{"d_penalty"};
}

1; # Got one?

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
