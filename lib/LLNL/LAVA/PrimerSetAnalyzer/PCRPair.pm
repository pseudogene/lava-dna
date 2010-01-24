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

package LLNL::LAVA::PrimerSetAnalyzer::PCRPair;

use strict;
use warnings;
use Carp;

use vars qw(@ISA);

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::PrimerSet::PCRPair; # Recieves these as parameters
use LLNL::LAVA::PrimerSetInfo; # Creates and returns these as results

use LLNL::LAVA::PrimerAnalyzer; # is-a  - NOT named as a direct sub-class!
@ISA = ("LLNL::LAVA::PrimerAnalyzer");

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::PrimerSetAnalyzer::PCRPair - Builds PrimerSetInfo results that
contain the results of a PCRPair being analyzed.

=head1 SYNOPSIS

  use LLNL::LAVA::PrimerSetAnalyzer::PCRPair;

  # Instantiation
  $analyzer = LLNL::LAVA::PrimerSetAnalyzer::PCRPair->new();

  # For now...
  # The penalty of the PCRPair is calculated as:
  #   Sum of Primer3 penalties for both primers
  # + Target amplicon size penalty
  # + Difference in melting temperature penalty
  
  # Eventually...
  # + Dimerization penalty?

  # Get the analysis results for a LLNL::LAVA::PCRPair
  $primerSetInfo = $analyzer->analyze($pcrPair);

  # Tag use
  $analyzer->setTag("useful_tag_name", $value);
  $exists = $analyzer->getTagExists("useful_tag_name");
  $value = $analyzer->getTag("useful_tag_name");

=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

PrimerSetAnalyzer::PCRPair is a PrimerSetAnalyzer that builds and returns PrimerInfo
objects that contain the results of analyzing an LLNL::LAVA::Oligo 

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $analyzer = LLNL::LAVA::PrimerSetAnalyzer::PCRPair->new();
 Function  : Creates a new LLNL::LAVA::PrimerSetAnalyzer::PCRPair
 Arguments : <n/a>
 Example   : See Usage
 Returns   : A new LLNL::LAVA::PrimerAnalyzer::PCRPair

=cut

sub new
{
  my ($classType, $paramHash_r) = @_;

  # Initializing as a PrimerAnalyzer (and TagHolder)
  my $this = $classType->SUPER::new();

  # Calculate the internal slope equations given targets
  # The weights of each parameter can be set dynamically for the Analyzer,
  # but the penalty equations will be set on Analyzer instantiation
  my $targetAmpliconLength = 
    optionWithDefault($paramHash_r, "target_amplicon_length", 80);

  # These are the min and max UN-WEIGHTED length penalties
  my $minLengthPenalty = 0;
  my $maxLengthPenalty = 10;

  # The current goal is to have a linear penalty increase
  # from zero at target length, increasing to the maximum penalty at 4x the
  # target length.  Currently no under-size penalty is assessed.

  # Slope of length penalty line
  # y = mx + b
  my $lengthCap = 4 * $targetAmpliconLength;
  my $lengthM = ($maxLengthPenalty - $minLengthPenalty) / 
    ($lengthCap - $targetAmpliconLength);
  my $lengthB = $minLengthPenalty - ($lengthM * $targetAmpliconLength); 
 
  # These are the min and max UN-WEIGHTED tm difference penalties
  my $minTMDiffPenalty = 0;
  my $maxTMDiffPenalty = 10; 

  # Linear penalty increas from zero at no TM difference to maximum penalty
  # at 15 degrees C

  # Slope of TM penalty line (0,0 is a guaranteed point on the line)
  # y = mx 
  my $tmDifferenceCap = 15;
  my $tmDiffM = ($maxTMDiffPenalty - $minTMDiffPenalty) / $tmDifferenceCap;
  
  # This sets the dictionary of available parameters, and sets the default
  # weights.
  $this->{"d_analyzerParameters"} =
    {
      # Objective function weights
      "primer3_weight" => "1.0", # Scale is approx. 0-15 for most single oligos
      "amplicon_length_weight" => "1.0",
      "tm_difference_weight" => "1.0",
    };

  # Default targets and penalty parameters
  $this->{"d_targetAmpliconLength"} = $targetAmpliconLength;
  $this->{"d_ampliconPenaltyLengthCap"} = $lengthCap;
  $this->{"d_lengthM"} = $lengthM;
  $this->{"d_lengthB"} = $lengthB;

  $this->{"d_tmPenaltyDifferenceCap"} = $tmDifferenceCap;
  $this->{"d_tmDifferenceM"} = $tmDiffM;

  return $this;
}

#-------------------------------------------------------------------------------

=head2 analyze

 Usage     : $primerSetInfo = $analyzer->analyze($pcrPair);
 Function  : Build a PrimerSetInfo object that represents the results of 
             analyzing an LLNL::LAVA::PCRPair within this analyzer's context. 
 Arguments : L<LLNL::LAVA::PrimerSet::PCRPair> - Primer pair to analyze 
 Example   : See Usage
 Returns   : L<LLNL::LAVA::PrimerSetInfo::PCRPair> - result of analyzing the primer

=cut

sub analyze
{
  my ($this, $pair) = @_;

  if(! defined($pair))
  {
    confess("programming error - first parameter is a required " .
      "\"LLNL::LAVA::PrimerSet::PCRPair\" object");
  }
  if(! $pair->isa("LLNL::LAVA::PrimerSet::PCRPair"))
  {
    confess("programming error - first parameter must be of type " .
      "\"LLNL::LAVA::PrimerSet::PCRPair\"");
  }

  # For now...
  # The penalty of the PCRPair is calculated as:
  #   Sum of Primer3 penalties for both primers
  # + Target amplicon size penalty
  # + Difference in melting temperature penalty
  
  # Eventually...
  # + Dimerization penalty

  my $forwardInfo = $pair->getForwardInfo();
  my $reverseInfo = $pair->getReverseInfo();

  my $primer3Penalty = $forwardInfo->getPenalty() + $reverseInfo->getPenalty();

  # Target amplicon size penalty
  # Not calculating an under-size penalty for now
  my $ampliconLength = $pair->getLength();

  # Bound the lengths the amplicon penalty is assessed on (clip the ends)
  my $ampLengthForPenalty = $ampliconLength;
  if($ampLengthForPenalty < $this->{"d_targetAmpliconLength"})
  {
    $ampLengthForPenalty = $this->{"d_targetAmpliconLength"};
  }
  elsif($ampLengthForPenalty > $this->{"d_ampliconPenaltyLengthCap"})
  {
    $ampLengthForPenalty = $this->{"d_ampliconPenaltyLengthCap"};
  }

  # Line equation is y = mx + b
  my $ampliconLengthPenalty = 
    ($this->{"d_lengthM"} * $ampLengthForPenalty) + $this->{"d_lengthB"};

  # Melting temperature difference penalty
  # min/max "likely" temperature differences for primers to help set penalty
  # Cap out at 20 degrees off? Maybe 10? (remember - half of the range)
  my $forwardTM = $forwardInfo->getTag("melting_temperature");
  my $reverseTM = $reverseInfo->getTag("melting_temperature");

  # Bound the range the TM difference penalty is assessed (clip the end)
  my $tmDiff = abs($forwardTM - $reverseTM);
  my $tmDiffForPenalty = $tmDiff;
  if($tmDiffForPenalty > $this->{"d_tmPenaltyDifferenceCap"})
  {
    $tmDiffForPenalty = $this->{"d_tmPenaltyDifferenceCap"};
  }

  # Line equation is y = mx
  my $tmDiffPenalty = $this->{"d_tmDifferenceM"} * $tmDiffForPenalty; 

  # PrimerAnalyzer has setParameters() and getParameterValue()
  my $primer3Weight = $this->getParameterValue("primer3_weight");
  my $ampliconLengthWeight = $this->getParameterValue("amplicon_length_weight");
  my $tmWeight = $this->getParameterValue("tm_difference_weight");

  my $penalty = ($primer3Weight * $primer3Penalty) +
    ($ampliconLengthWeight * $ampliconLengthPenalty) +
    ($tmWeight * $tmDiffPenalty);

  my $primerSetInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new(
    {
      "penalty" => $penalty,
      "analyzed_pair" => $pair,
    });

  $primerSetInfo->setTag("analyzed_pair", $pair);

  # Temporarily keeping these for debugging and analytics
  $primerSetInfo->setTag("primer3_penalty", $primer3Penalty);
  $primerSetInfo->setTag("amplicon_length_penalty", $ampliconLengthPenalty);
  $primerSetInfo->setTag("tm_difference_penalty", $tmDiffPenalty);

  return $primerSetInfo;
}

1; # Lame!

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
