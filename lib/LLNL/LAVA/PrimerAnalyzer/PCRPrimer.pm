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

package LLNL::LAVA::PrimerAnalyzer::PCRPrimer;

use strict;
use vars qw(@ISA);
use Carp;

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::Oligo; # Recieves these as parameters
use LLNL::LAVA::PrimerInfo; # Returns these as results

use LLNL::LAVA::PrimerAnalyzer; # is-a
@ISA = ("LLNL::LAVA::PrimerAnalyzer");

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::PrimerAnalyzer::PRCPrimer - Builds PrimerInfo results that contain
the results of an oligo being analyzed as a PCR primer.

=head1 SYNOPSIS

  use LLNL::LAVA::PrimerAnalyzer::PCRPrimer;

  # Instantiation
  $analyzer = LLNL::LAVA::PrimerAnalyzer::PCRPrimer->new();

  # Parameters aren't relevant for PCRPrimer right now.  Getting the penalty
  # of an ordinary Oligo as a PCRPrimer is equivalent to asking for the 
  # Primer3 penalty of the Oligo.
  
  # Eventually, the primer penalty will also include a C-to-G ratio penalty.


  # Get the analysis results for an oligo.
  $primerInfo = $analyzer->analyze($oligo);

  # Tag use
  $analyzer->setTag("useful_tag_name", $value);
  $exists = $analyzer->getTagExists("useful_tag_name");
  $value = $analyzer->getTag("useful_tag_name");

=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

PrimerAnalyzer::PRCPrimer is a PrimerAnalyzer that builds and returns PrimerInfo
objects that contain the results of analyzing an LLNL::LAVA::Oligo 

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $analyzer = LLNL::LAVA::PrimerAnalyzer::PRCPrimer->new();
 Function  : Creates a new LLNL::LAVA::PrimerAnalyzer::PRCPrimer, only sub-classes 
             should use this
 Arguments : <n/a>
 Example   : See Usage
 Returns   : A new LLNL::LAVA::PrimerAnalyzer::PRCPrimer

=cut

sub new
{
  my ($classType) = @_;

  # Initialize via PrimerAnalyzer (and TagHolder)
  my $this = $classType->SUPER::new();

  return $this;
}

#-------------------------------------------------------------------------------

=head2 analyze

 Usage     : $primerInfo = $analyzer->analyze($oligo);
 Function  : Build a PrimerInfo object that represents the results of analyzing
             an oligo within this analyzer's context. 
 Arguments : L<LLNL::LAVA::Oligo> - individual oligo to analyze as a primer
 Example   : See Usage
 Returns   : L<LLNL::LAVA::PrimerInfo> - the result of analyzing the primer

=cut

sub analyze
{
  my ($this, $oligo) = @_;

  if(! defined($oligo))
  {
    confess("programming error - first parameter is a required " .
      "\"LLNL::LAVA::Oligo\" object");
  }
  if(! $oligo->isa("LLNL::LAVA::Oligo"))
  {
    confess("programming error - first parameter must be of type " .
      "\"LLNL::LAVA::Oligo\"");
  }

  # This is the simplest case.  There's only a single sequence for this type 
  # of oligo, so we build a PrimerInfo where the penalty is the primer3 
  # penalty of the Oligo, and the rest of the PrimerInfo values are copied
  # straignt from the Oligo
  my $penalty = $oligo->getTag("primer3_penalty");

  my $sequence = $oligo->sequence();
  my $location = $oligo->location();
  my $length = $oligo->length();
  my $strand =  $oligo->getTag("strand");
  my $tm = $oligo->getTag("primer3_tm");

  # TODO: Put the C-to-G ratio penalty back in as a factor

  my $primerInfo = LLNL::LAVA::PrimerInfo->new(
    {
      "penalty" => $penalty,
      "sequence" => $sequence,
      "location" => $location,
      "length" => $length,
      "analyzed_primer" => $oligo,
    });

  # Think we still want to keep the notion of "strand" kinda fuzzy, since we
  # will be mixing the "strand" of different parts of a single oligo for LAMP
  $primerInfo->setTag("strand", $strand);

  # TM is a tag because TM could be assessed in several contexts, and might
  # be particularly difficult to represent with a data member when paired 
  # with different specific anti-sense sequences
  $primerInfo->setTag("melting_temperature", $tm);

  #print "\"$penalty\"";

  return $primerInfo;
}

1; # Lame!

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
