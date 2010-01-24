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

package LLNL::LAVA::PrimerSet::LAMP;

use strict;
use warnings;
use Carp;

use vars qw(@ISA);

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::PrimerSetInfo::PCRPair; # Recieves as parameter and holds

use LLNL::LAVA::PrimerSet; # is-a
@ISA = ("LLNL::LAVA::PrimerSet");

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::PrimerSet::LAMP - Three LLNL::LAVA::PrimerSetInfo::PCRPair objects

=head1 SYNOPSIS

  # Instantiation
  $lampSignature = LLNL::LAVA::PrimerSet::LAMP->new(
    {
      # Strand orientation gets pretty weird here, so we're going to require
      # that all PCRPair-Infos have plus strand forward primers and minus
      # strand reverse primers.
      "inner_info" => $innerInfo,
      "middle_info" => $middleInfo,
      "outer_info" => $outerInfo,
    });

  $innerInfo = $lampSignature->getInnerInfo();
  $middleInfo = $lampSignature->getMiddleInfo();
  $outerInfo = $lampSignature->getOuterInfo();
  
  # But since we're going to be strict about how our internal pairs are, we can
  # also do cool things like:
  $outerForward = $lampSignature->getF3();
  $outerReverse = $lampSignature->getB3();
  $fip = $lampSignature->getBIP();
  $bip = $lampSignature->getFIP();
 
  # Helper for debugging and status printout
  $locationSummary = $lampSignature->getLocationSummary();
   
  # Get/set pair for the BIP and FIP linker
  $linker = $lampSignature->linker();
  $lampSignature->linker($newLinker);
 
  
  # PrimerSet inherited functions
  $location = $pair->getStartLocation();
  $location = $pair->getEndLocation();
  $length = $pair->getLength();
  $rangeString = $pair->getRangeAsString(); # String in the form of "start-end"
   
  # Tag use (inherited from TagHolder)
  $oligo->setTag("useful_statistic", $usefulStatistic);
  $exists = $oligo->getTagExists("useful_statistic");
  $tm = $oligo->getTag("useful_statistic");

=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

PrimerSet::LAMP is the object model of a pair of oligos from opposite strands that
can be used as PCR primers.

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $sigComponent = LLNL::LAVA::PrimerSet::LAMP->new(
               {
                 # We're requiring that all PCRPair-Infos have plus strand 
                 # forward primers and minus strand reverse primers.
                 "inner_info" => $innerInfo,
                 "middle_info" => $middleInfo,
                 "outer_info" => $outerInfo,
               });
 Function  : Creates a new LLNL::LAVA::PrimerSet::LAMP 
 Arguments : Hash ref - Parameters for initialization including:
               inner_info - LLNL::LAVA::PrimerSetInfo::PCRPair
               middle_info - LLNL::LAVA::PrimerSetInfo::PCRPair
               outer_info - LLNL::LAVA::PrimerSetInfo::PCRPair
 Example   : See Usage
 Returns   : A new LLNL::LAVA::PrimerSet::LAMP

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

  # Read component properties
  my $innerInfo = optionRequired($paramHash_r, "inner_info");
  my $middleInfo = optionRequired($paramHash_r, "middle_info");
  my $outerInfo = optionRequired($paramHash_r, "outer_info");

  my $linker = optionWithDefault($paramHash_r, "linker", "");

  my $innerPair = $innerInfo->getAnalyzedPair();
  my $middlePair = $middleInfo->getAnalyzedPair();
  my $outerPair = $outerInfo->getAnalyzedPair();
  
  # Let's enforce all 6 strand orientations...
  if($innerPair->getForwardStrand() ne "plus" ||
     $middlePair->getForwardStrand() ne "plus" ||
     $outerPair->getForwardStrand() ne "plus")
  {
    confess("data error - all forward primers need to be on the plus strand so " .
      "we can guarantee proper rendering of the LAMP signature.");
  }
  if($innerPair->getReverseStrand() ne "minus" ||
     $middlePair->getReverseStrand() ne "minus" ||
     $outerPair->getReverseStrand() ne "minus")
  {
    confess("data error - all reverse primers need to be on the minus strand so " .
      "we can guarantee proper rendering of the LAMP signature.");
  }

  # Initialize with locations and length based on the outer primers - this has
  # a subtle dependence on the guaranteed strand-orientation of the primers
  # because of how the ranges are calculated
  my $outerStart = $outerPair->getStartLocation();
  my $outerEnd = $outerPair->getEndLocation();
  my $outerLength = $outerPair->getLength(); 
   
  # Initialize via PrimerSet (and TagHolder)
  my $this = $classType->SUPER::new(
    {
      "start_location" => $outerStart,
      "end_location" => $outerEnd,
      "length" => $outerLength,
    });

  # Handle component properties specialized for LAMP PrimerSets
  $this->{"d_innerInfo"} = $innerInfo;
  $this->{"d_middleInfo"} = $middleInfo;
  $this->{"d_outerInfo"} = $outerInfo;
  $this->{"d_linker"} = $linker;

  return $this;
}

#-------------------------------------------------------------------------------

=head2 getInnerInfo

 Usage     : $pcrPairInfo = $lampSignature->getInnerInfo();
 Function  : Get a reference to this signature's inner pair's info
 Arguments : <n/a>
 Example   : See Usage
 Returns   : LLNL::LAVA::PrimerSetInfo::PCRPair - inner pair's info

=cut

sub getInnerInfo
{
  my ($this) = @_;

  return $this->{"d_innerInfo"};
}

#-------------------------------------------------------------------------------

=head2 getMiddleInfo

 Usage     : $pcrPairInfo = $lampSignature->getMiddleInfo();
 Function  : Get a reference to this signature's middle pair's info
 Arguments : <n/a>
 Example   : See Usage
 Returns   : LLNL::LAVA::PrimerSetInfo::PCRPair - middle pair's info

=cut

sub getMiddleInfo
{
  my ($this) = @_;

  return $this->{"d_middleInfo"};
}

#-------------------------------------------------------------------------------

=head2 getOuterInfo

 Usage     : $pcrPairInfo = $lampSignature->getOuterInfo();
 Function  : Get a reference to this signature's outer pair's info
 Arguments : <n/a>
 Example   : See Usage
 Returns   : LLNL::LAVA::PrimerSetInfo::PCRPair - outer pair's info

=cut

sub getOuterInfo
{
  my ($this) = @_;

  return $this->{"d_outerInfo"};
}

#-------------------------------------------------------------------------------

=head2 getF3

 Usage     : $f3Sequence = $lampSignature->getF3();
 Function  : Get the sequence for this signature's outer forward primer
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - outer forward primer's sequence

=cut

sub getF3
{
  my ($this) = @_;

  my $pair = $this->getOuterInfo()->getAnalyzedPair();
  my $primer = $pair->getForwardInfo()->getSequence();
    
  return $primer;
}

#-------------------------------------------------------------------------------

=head2 getB3

 Usage     : $b3Sequence = $lampSignature->getB3();
 Function  : Get the sequence for this signature's outer reverse primer
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - outer reverse primer's sequence

=cut

sub getB3
{
  my ($this) = @_;

  my $pair = $this->getOuterInfo()->getAnalyzedPair();
  my $primer = $pair->getReverseInfo()->getSequence();
    
  return $primer;
}

#-------------------------------------------------------------------------------

=head2 getFIP

 Usage     : $fipSequence = $lampSignature->getFIP();
 Function  : Get the sequence for this signature's FIP oligo complex 
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - Complete FIP oligo sequence 

=cut

sub getFIP
{
  my ($this) = @_;

  my $innerPair = $this->getInnerInfo()->getAnalyzedPair();
  my $middlePair = $this->getMiddleInfo()->getAnalyzedPair();

  # FIP is inner-forward-revcomp -> linker -> middle-forward
  my $innerForward = $innerPair->getForwardInfo->getSequence();
  my $middleForward = $middlePair->getForwardInfo->getSequence();
  my $linker = $this->linker();

  my $revCompInnerForward = $this->getRevcomp($innerForward);

  my $sequence = $revCompInnerForward . $linker . $middleForward;

  return $sequence;
}

#-------------------------------------------------------------------------------

=head2 getBIP

 Usage     : $bipSequence = $lampSignature->getBIP();
 Function  : Get the sequence for this signature's BIP oligo complex 
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - Complete BIP oligo sequence 

=cut

sub getBIP
{
  my ($this) = @_;

  my $innerPair = $this->getInnerInfo()->getAnalyzedPair();
  my $middlePair = $this->getMiddleInfo()->getAnalyzedPair();

  # BIP is inner-reverse-revcomp -> linker -> middle-reverse
  my $innerReverse = $innerPair->getReverseInfo->getSequence();
  my $middleReverse = $middlePair->getReverseInfo->getSequence();
  my $linker = $this->linker();

  my $revCompInnerReverse = $this->getRevcomp($innerReverse);

  my $sequence = $revCompInnerReverse . $linker . $middleReverse;

  return $sequence;
}

#-------------------------------------------------------------------------------

=head2 getF2

 Usage     : $f2Sequence = $lampSignature->getF2();
 Function  : Get the sequence for this signature's middle forward primer
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - middle forward primer's sequence

=cut

sub getF2
{
  my ($this) = @_;

  my $pair = $this->getMiddleInfo()->getAnalyzedPair();
  my $primer = $pair->getForwardInfo()->getSequence();
    
  return $primer;
}

#-------------------------------------------------------------------------------

=head2 getB2

 Usage     : $b2Sequence = $lampSignature->getB2();
 Function  : Get the sequence for this signature's middle reverse primer
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - middle reverse primer's sequence

=cut

sub getB2
{
  my ($this) = @_;

  my $pair = $this->getMiddleInfo()->getAnalyzedPair();
  my $primer = $pair->getReverseInfo()->getSequence();
    
  return $primer;
}

#-------------------------------------------------------------------------------

=head2 getF1

 Usage     : $f1Sequence = $lampSignature->getF1();
 Function  : Get the sequence for this signature's inner forward primer
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - inner forward primer's sequence

=cut

sub getF1
{
  my ($this) = @_;

  my $pair = $this->getInnerInfo()->getAnalyzedPair();
  my $primer = $pair->getForwardInfo()->getSequence();
    
  return $primer;
}

#-------------------------------------------------------------------------------

=head2 getB1

 Usage     : $b1Sequence = $lampSignature->getB1();
 Function  : Get the sequence for this signature's inner reverse primer
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - inner reverse primer's sequence

=cut

sub getB1
{
  my ($this) = @_;

  my $pair = $this->getInnerInfo()->getAnalyzedPair();
  my $primer = $pair->getReverseInfo()->getSequence();
    
  return $primer;
}

#-------------------------------------------------------------------------------

=head2 linker

 Usage     : $linker = $lampSignature->linker();
             $lampSignature->linker($newLinker);
 Function  : Get/Set pair for the FIP and BIP linker
 Arguments : (optional) String - new linker
 Example   : See Usage
 Returns   : String - the FIP and BIP linker to use

=cut

sub linker
{
  my ($this, $newLinker) = @_;

  if(defined $newLinker)
  {
    $this->{"d_linker"} = $newLinker;
  }

  return $this->{"d_linker"};
}

#-------------------------------------------------------------------------------

=head2 getRevcomp

 Usage     : $reverseComplementSequence = $lampSignature->getRevcomp($sequence);
 Function  : More general helper function (TODO: belongs somewhere else).  
             This revcomps the parameter sequence, and does nothing to this object
 Arguments : String - Sequence to find the reverse complement of
 Example   : See Usage
 Returns   : String - The reverse-complemented sequence

=cut

sub getRevcomp
{
  my ($this, $sequence) = @_;

  # Invert the sequence alphabet, and reverse the sequence in 2 steps, then
  # set this object's "sequence" to the newly calculated reverse complement
  $sequence =~ 
    tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  my $revComp = reverse($sequence);

  return $revComp;
}

#-------------------------------------------------------------------------------

=head2 getLocationSummary

 Usage     : $locationSummary = $lampSignature->getLocationSummary();
 Function  : Get a string describing the primer locations for this signature -
             This is primarily for debugging
             Kinda relies on forward primers being plus strand and reverse being minus
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - locations and spans of each primer locus

=cut

sub getLocationSummary
{
  my ($this) = @_;

  my $outerInfo = $this->getOuterInfo();
  my $middleInfo = $this->getMiddleInfo();
  my $innerInfo = $this->getInnerInfo();

  my $summary = "";
  my @pairInfos = ($outerInfo, $middleInfo, $innerInfo);
  foreach my $pairInfo(@pairInfos)
  {
    my $pair = $pairInfo->getAnalyzedPair();

    my $forwardInfo = $pair->getForwardInfo();
    my $forwardLocation = $forwardInfo->getLocation();
    my $forwardLength = $forwardInfo->getLength();

    my $reverseInfo = $pair->getReverseInfo();
    my $reverseLocation = $reverseInfo->getLocation();
    my $reverseLength = $reverseInfo->getLength();

    $summary .= "($forwardLocation-" .
      ($forwardLocation + $forwardLength - 1) .
      ", " .
      ($reverseLocation - $reverseLength + 1) .
      "-$reverseLocation)"; 
  } 

  return $summary;
}

#-------------------------------------------------------------------------------

=head2 getLoopLocationSummary

 Usage     : $locationSummary = $lampSignature->getLoopLocationSummary();
 Function  : Get a string describing the loop primer locations for this signature -
             This is primarily for debugging
             Kinda relies on left loop primers being on the minus strand, and right
	       loop primers being on the plus strand
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - locations and spans of each primer locus

=cut

sub getLoopLocationSummary
{
  my ($this) = @_;

  # TODO: this might belong in a sub-class or a composed class for
  # signatures with loop primers
  if($this->getTagExists("has_loop_primers") == $FALSE)
  {
    confess("Cannot retrieve loop primer location summary since this signature " . 
      "doesn't appear to have loop primers.");
  }

  my $floopInfo = $this->getTag("floop_info");
  my $bloopInfo = $this->getTag("bloop_info");

  my $floopLocation = $floopInfo->getLocation();
  my $floopLength = $floopInfo->getLength();
  my $bloopLocation = $bloopInfo->getLocation();
  my $bloopLength = $bloopInfo->getLength();

  my $summary = "(" .
    ($floopLocation - $floopLength + 1) .
    "-$floopLocation, $bloopLocation-" .
    ($bloopLocation + $bloopLength - 1) .
    ")";

  return $summary;
}

#-------------------------------------------------------------------------------

=head2 getPenaltySummary

 Usage     : $penaltySummary = $lampSignature->getPenaltySummary();
 Function  : Get a string describing the primer penalties for this signature -
             This is primarily for debugging
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - penalties of each primer and pair

=cut

sub getPenaltySummary
{
  my ($this) = @_;

  my $outerInfo = $this->getOuterInfo();
  my $middleInfo = $this->getMiddleInfo();
  my $innerInfo = $this->getInnerInfo();

  my $summary = "";
  my @pairInfos = ($outerInfo, $middleInfo, $innerInfo);
  foreach my $pairInfo(@pairInfos)
  {
    my $pairPenalty = $pairInfo->getPenalty();
    my $pair = $pairInfo->getAnalyzedPair();

    my $forwardInfo = $pair->getForwardInfo();
    my $forwardPenalty = $forwardInfo->getPenalty();

    my $reverseInfo = $pair->getReverseInfo();
    my $reversePenalty = $reverseInfo->getPenalty();

    $summary .= "($pairPenalty" . # Newline to avoid array context...
      "[$forwardPenalty, $reversePenalty])";
  } 

  return $summary;
}

#-------------------------------------------------------------------------------

=head2 getTMSummary

 Usage     : $tmSummary = $lampSignature->getTMSummary();
 Function  : Get a string describing the primer tms for this signature -
             This is primarily for debugging
             Relies on primerInfo's having a "melting_temperature" tag
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - penalties of each primer and pair

=cut

sub getTMSummary
{
  my ($this) = @_;

  my $outerInfo = $this->getOuterInfo();
  my $middleInfo = $this->getMiddleInfo();
  my $innerInfo = $this->getInnerInfo();

  my $summary = "";
  my @pairInfos = ($outerInfo, $middleInfo, $innerInfo);
  foreach my $pairInfo(@pairInfos)
  {
    my $pair = $pairInfo->getAnalyzedPair();

    my $forwardInfo = $pair->getForwardInfo();
    my $forwardTM = $forwardInfo->getTag("melting_temperature");

    my $reverseInfo = $pair->getReverseInfo();
    my $reverseTM = $reverseInfo->getTag("melting_temperature");

    $summary .= "($forwardTM, $reverseTM)";
  } 

  return $summary;
}

1; # Lame!

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
