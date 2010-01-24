#!/usr/local/bin/perl -w

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

use strict;
use warnings;
use Carp;

use Getopt::Long;
use Data::Dumper; # Probably just for debugging

use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::LocatableSeq;

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::OligoEnumerator::Primer3Conserved;

use LLNL::LAVA::PrimerAnalyzer::PCRPrimer;
use LLNL::LAVA::PrimerInfo;
use LLNL::LAVA::PrimerSet::PCRPair;

use LLNL::LAVA::PrimerSetAnalyzer::PCRPair;
use LLNL::LAVA::PrimerSetInfo::PCRPair;

use LLNL::LAVA::PrimerSet::LAMP;


$| = 0;

{ # Fake main() to enforce scope
  my %options;
  my %optionMap =
    (
      "alignment_fasta=s" => \$options{"alignment_fasta"},
      "output_file=s" => \$options{"output_file"}, 
      "signature_max_length=i" => \$options{"signature_max_length"},

      "outer_primer_target_length=i" => \$options{"outer_primer_target_length"},
      "outer_primer_min_length=i" => \$options{"outer_primer_min_length"},
      "outer_primer_max_length=i" => \$options{"outer_primer_max_length"},
      "outer_primer_target_tm=f" => \$options{"outer_primer_target_tm"},
      "outer_primer_min_tm=f" => \$options{"outer_primer_min_tm"},
      "outer_primer_max_tm=f" => \$options{"outer_primer_max_tm"},

      "loop_primer_target_length=i" => \$options{"loop_primer_target_length"},
      "loop_primer_min_length=i" => \$options{"loop_primer_min_length"},
      "loop_primer_max_length=i" => \$options{"loop_primer_max_length"},
      "loop_primer_target_tm=f" => \$options{"loop_primer_target_tm"},
      "loop_primer_min_tm=f" => \$options{"loop_primer_min_tm"},
      "loop_primer_max_tm=f" => \$options{"loop_primer_max_tm"},

      "middle_primer_target_length=i" => \$options{"middle_primer_target_length"},
      "middle_primer_min_length=i" => \$options{"middle_primer_min_length"},
      "middle_primer_max_length=i" => \$options{"middle_primer_max_length"},
      "middle_primer_target_tm=f" => \$options{"middle_primer_target_tm"},
      "middle_primer_min_tm=f" => \$options{"middle_primer_min_tm"},
      "middle_primer_max_tm=f" => \$options{"middle_primer_max_tm"},

      "inner_primer_target_length=i" => \$options{"inner_primer_target_length"},
      "inner_primer_min_length=i" => \$options{"inner_primer_min_length"},
      "inner_primer_max_length=i" => \$options{"inner_primer_max_length"},
      "inner_primer_target_tm=f" => \$options{"inner_primer_target_tm"},
      "inner_primer_min_tm=f" => \$options{"inner_primer_min_tm"},
      "inner_primer_max_tm=f" => \$options{"inner_primer_max_tm"},
    
      "max_poly_bases=i" => \$options{"max_poly_bases"}, 

      "outer_pair_target_length=i" => \$options{"outer_pair_target_length"},
      "middle_pair_target_length=i" => \$options{"middle_pair_target_length"},
      "inner_pair_target_length=i" => \$options{"inner_pair_target_length"},

      "include_loop_primers!" => \$options{"include_loop_primers"},
      "loop_min_gap=i" => \$options{"loop_min_gap"},
      "min_signatures_for_success=i" => \$options{"min_signatures_for_success"},
      "min_primer_spacing=i" => \$options{"min_primer_spacing"},
      "min_inner_pair_spacing=i" => \$options{"min_inner_pair_spacing"},

      "primer3_executable=s" => \$options{"primer3_executable"},
      "alignment_format=s" => \$options{"alignment_format"},

      # TODO: Not sure if the pair target lengths should be exposed to the 
      # user, or adjusted based on other parameters
      #"outer_pair_target_length=i" => \$options{"outer_pair_target_length"}, 
      #"middle_pair_target_length=i" => \$options{"middle_pair_target_length"}, 
      #"inner_pair_target_length=i" => \$options{"inner_pair_target_length"}, 

      "option_file|options_file=s" => \$options{"option_file"},
    );

  my %optionDefaults =
    (
      "signature_max_length" => 320,
      "outer_primer_target_length" => 20,
      "outer_primer_min_length" => 18,
      "outer_primer_max_length" => 23,
      "outer_primer_target_tm" => "60.0",
      "loop_primer_target_length" => 20,
      "loop_primer_min_length" => 18,
      "loop_primer_max_length" => 23,
      "loop_primer_target_tm" => "60.0",
      "middle_primer_target_length" => 20,
      "middle_primer_min_length" => 18,
      "middle_primer_max_length" => 23,
      "middle_primer_target_tm" => "60.0",
      "inner_primer_target_length" => 23,
      "inner_primer_min_length" => 20,
      "inner_primer_max_length" => 26,
      "inner_primer_target_tm" => "62.0",
      #"inner_primer_target_length" => 19,
      #"inner_primer_min_length" => 17,
      #"inner_primer_max_length" => 24,
      #"inner_primer_target_tm" => "57.0",
      "max_poly_bases" => 5,
      "include_loop_primers" => $TRUE,
      "loop_min_gap" => 25,
      "min_signatures_for_success" => 1, # Should probably never go lower
      "min_primer_spacing" => 1,
      "min_inner_pair_spacing" => 1,
      # Some LAMP-specific approximate targets for a "minimum sized" signature
      # Currently, no penalty is assessed for lengths under the target size, so
      # these sizes are a little larger than they need to be.
      "outer_pair_target_length" => 200, 
      "middle_pair_target_length" => 160, 
      "inner_pair_target_length" => 50, 
      "primer3_executable" => "/usr/bin/primer3_core",
      "alignment_format" => "fasta",
    );

  my $usageString = "Usage:\n" .
    "./lava.pl \n" .
      "    --alignment_fasta <fasta_file>\n" .
      "    --output_file <output_file>\n" .
      "    [--signature_max_length <length, deafult=" .
        $optionDefaults{"signature_max_length"} .
	">]\n" .
      # Outer primer options
      "    [--outer_primer_target_length <length, default=" .
        $optionDefaults{"outer_primer_target_length"} .
	">]\n" .
      "    [--outer_primer_min_length <length, default=" .
        $optionDefaults{"outer_primer_min_length"} .
	">]\n" .
      "    [--outer_primer_max_length <length, default=" .
        $optionDefaults{"outer_primer_max_length"} .
	">]\n" .
      "    [--outer_primer_target_tm <tm, default=" .
        $optionDefaults{"outer_primer_target_tm"} .
        "C>]\n" .
      "    [--outer_primer_min_tm <tm, default=outer_primer_target_tm - 1.0>]\n" .
      "    [--outer_primer_max_tm <tm, default=outer_primer_target_tm + 1.0>]\n" .
      # Loop primer options
      "    [--loop_primer_target_length <length, default=" . 
        $optionDefaults{"loop_primer_target_length"} .
        ">]\n" .
      "    [--loop_primer_min_length <length, default=" .
        $optionDefaults{"loop_primer_min_length"} .
	">]\n" .
      "    [--loop_primer_max_length <length, default=" .
        $optionDefaults{"loop_primer_max_length"} .
	">]\n" .
      "    [--loop_primer_target_tm <tm, default=" .
        $optionDefaults{"loop_primer_target_tm"} .
	"C>]\n" .
      "    [--loop_primer_min_tm <tm, default=loop_primer_target_tm - 1.0>]\n" .
      "    [--loop_primer_max_tm <tm, default=loop_primer_target_tm + 1.0>]\n" .
      # Middle primer options
      "    [--middle_primer_target_length <length, default=" .
        $optionDefaults{"middle_primer_target_length"} .
	">]\n" .
      "    [--middle_primer_min_length <length, default=" .
        $optionDefaults{"middle_primer_min_length"} .
	">]\n" .
      "    [--middle_primer_max_length <length, default=" .
        $optionDefaults{"middle_primer_max_length"} .
	">]\n" .
      "    [--middle_primer_target_tm <tm, default=" .
        $optionDefaults{"middle_primer_target_tm"} .
	"C>]\n" .
      "    [--middle_primer_min_tm <tm, default=middle_primer_target_tm - 1.0>]\n" .
      "    [--middle_primer_max_tm <tm, default=middle_primer_target_tm + 1.0>]\n" .
      # Outer primer options
      "    [--inner_primer_target_length <length, default=" .
        $optionDefaults{"inner_primer_target_length"} .
	">]\n" .
      "    [--inner_primer_min_length <length, default=" .
        $optionDefaults{"inner_primer_min_length"} .
	">]\n" .
      "    [--inner_primer_max_length <length, default=" .
        $optionDefaults{"inner_primer_max_length"} .
	">]\n" .
      "    [--inner_primer_target_tm <tm, default=" .
        $optionDefaults{"inner_primer_target_tm"} .
	"C>]\n" .
      "    [--inner_primer_min_tm <tm, default=inner_primer_target_tm - 1.0>]\n" .
      "    [--inner_primer_max_tm <tm, default=inner_primer_target_tm + 1.0>]\n" .
      # Other kinds of options
      "    [--max_poly_bases <max, default=" .
        $optionDefaults{"max_poly_bases"} .
	">]\n" .
      "    [--min_primer_spacing <max, default=" .
        $optionDefaults{"min_primer_spacing"} .  
        ">]\n" .
      "    [--min_inner_pair_spacing <max, default=" .
        $optionDefaults{"min_inner_pair_spacing"} .
        ">]\n" .

      # TODO: Not sure if the pair target lengths should be exposed to the 
      # user, or adjusted based on other parameters
      #      "    [--outer_pair_target_length <length, default=" .
      #        $optionDefaults{"outer_pair_target_length"} .
      #        ">]\n" .
      #      "    [--middle_pair_target_length <length, default=" .
      #      $optionDefaults{"middle_pair_target_length"} .
      #      ">]\n" .
      #      "    [--inner_pair_target_length <length, default=" .
      #      $optionDefaults{"inner_pair_target_length"} .
      #      ">]\n" .

      "    [--include_loop_primers <length, default=" .
        $optionDefaults{"include_loop_primers"} .
	">]\n" .
      # Loop gap is the distance between the middle and inner primers
      "    [--loop_min_gap <length, default=" .
        $optionDefaults{"loop_min_gap"} .
	">]\n" .
      "    [--min_signatures_for_success <length, default=" .
        $optionDefaults{"min_signatures_for_success"} .
	">]\n" .
      "    [--primer3_executable <path_to_primer3, default=" .
        $optionDefaults{"primer3_executable"} .
	">]\n" .
      "    [--alignment_format <file format of alignment, default=\"" .
        $optionDefaults{"alignment_format"} .
	"\">]\n" .
      "    [--option_file <options_xml> (cmd line options take precedence)]\n";

  # TODO: Probably want to be able to use multiple files for parameter
  # definition, so we can have the thermo parameters set, and separately have
  # the file IO parameters.
  GetOptions(%optionMap);
  loadOptionsFromFile(\%options);
  my $options_r = \%options;

  # TODO: perl-check for file existence cause BioPerl dump isn't useful
  my $alignmentFastaName = optionRequired($options_r, "alignment_fasta", $usageString);
  my $outputFileName = optionRequired($options_r, "output_file", $usageString);

  my $signatureMaxLength = 
    optionWithDefault($options_r, "signature_max_length", 
      $optionDefaults{"signature_max_length"});

  my $outerPrimerTargetLength =
    optionWithDefault($options_r, "outer_primer_target_length", 
      $optionDefaults{"outer_primer_target_length"});
  my $outerPrimerMinLength =
    optionWithDefault($options_r, "outer_primer_min_length", 
      $optionDefaults{"outer_primer_min_length"});
  if($outerPrimerMinLength > $outerPrimerTargetLength)
  {
    $outerPrimerMinLength = $outerPrimerTargetLength;
  }
  my $outerPrimerMaxLength =
    optionWithDefault($options_r, "outer_primer_max_length", 
      $optionDefaults{"outer_primer_max_length"});
  if($outerPrimerMaxLength < $outerPrimerTargetLength)
  {
    $outerPrimerMaxLength = $outerPrimerTargetLength;
  }

  my $outerPrimerTargetTM =
    optionWithDefault($options_r, "outer_primer_target_tm", 
      $optionDefaults{"outer_primer_target_tm"});
  my $outerPrimerMinTM =
    optionWithDefault($options_r, "outer_primer_min_tm", 
      ($outerPrimerTargetTM - 1.0));
  if($outerPrimerMinTM > $outerPrimerTargetTM)
  {
    $outerPrimerMinTM = $outerPrimerTargetTM;
  }
  my $outerPrimerMaxTM =
    optionWithDefault($options_r, "outer_primer_max_tm", 
      ($outerPrimerTargetTM + 1.0));
  if($outerPrimerMaxTM < $outerPrimerTargetTM)
  {
    $outerPrimerMaxTM = $outerPrimerTargetTM;
  }

  my $loopPrimerTargetLength =
    optionWithDefault($options_r, "loop_primer_target_length", 
      $optionDefaults{"loop_primer_target_length"});
  my $loopPrimerMinLength =
    optionWithDefault($options_r, "loop_primer_min_length", 
      $optionDefaults{"loop_primer_min_length"});
  if($loopPrimerMinLength > $loopPrimerTargetLength)
  {
    $loopPrimerMinLength = $loopPrimerTargetLength;
  }
  my $loopPrimerMaxLength =
    optionWithDefault($options_r, "loop_primer_max_length", 
      $optionDefaults{"loop_primer_max_length"});
  if($loopPrimerMaxLength < $loopPrimerTargetLength)
  {
    $loopPrimerMaxLength = $loopPrimerTargetLength;
  }

  my $loopPrimerTargetTM =
    optionWithDefault($options_r, "loop_primer_target_tm", 
      $optionDefaults{"loop_primer_target_tm"});
  my $loopPrimerMinTM =
    optionWithDefault($options_r, "loop_primer_min_tm", 
      ($loopPrimerTargetTM - 1.0));
  if($loopPrimerMinTM > $loopPrimerTargetTM)
  {
    $loopPrimerMinTM = $loopPrimerTargetTM;
  }
  my $loopPrimerMaxTM =
    optionWithDefault($options_r, "loop_primer_max_tm", 
      ($loopPrimerTargetTM + 1.0));
  if($loopPrimerMaxTM < $loopPrimerTargetTM)
  {
    $loopPrimerMaxTM = $loopPrimerTargetTM;
  }


  my $middlePrimerTargetLength =
    optionWithDefault($options_r, "middle_primer_target_length", 
      $optionDefaults{"middle_primer_target_length"});
  my $middlePrimerMinLength =
    optionWithDefault($options_r, "middle_primer_min_length", 
      $optionDefaults{"middle_primer_min_length"});
  if($middlePrimerMinLength > $middlePrimerTargetLength)
  {
    $middlePrimerMinLength = $middlePrimerTargetLength;
  }
  my $middlePrimerMaxLength =
    optionWithDefault($options_r, "middle_primer_max_length", 
      $optionDefaults{"middle_primer_max_length"});
  if($middlePrimerMaxLength < $middlePrimerTargetLength)
  {
    $middlePrimerMaxLength = $middlePrimerTargetLength;
  }

  my $middlePrimerTargetTM =
    optionWithDefault($options_r, "middle_primer_target_tm", 
      $optionDefaults{"middle_primer_target_tm"});
  my $middlePrimerMinTM =
    optionWithDefault($options_r, "middle_primer_min_tm", 
      ($middlePrimerTargetTM - 1.0));
  if($middlePrimerMinTM > $middlePrimerTargetTM)
  {
    $middlePrimerMinTM = $middlePrimerTargetTM;
  }
  my $middlePrimerMaxTM =
    optionWithDefault($options_r, "middle_primer_max_tm", 
      ($middlePrimerTargetTM + 1.0));
  if($middlePrimerMaxTM < $middlePrimerTargetTM)
  {
    $middlePrimerMaxTM = $middlePrimerTargetTM;
  }

  my $innerPrimerTargetLength =
    optionWithDefault($options_r, "inner_primer_target_length", 
      $optionDefaults{"inner_primer_target_length"});
  my $innerPrimerMinLength =
    optionWithDefault($options_r, "inner_primer_min_length", 
      $optionDefaults{"inner_primer_min_length"});
  if($innerPrimerMinLength > $innerPrimerTargetLength)
  {
    $innerPrimerMinLength = $innerPrimerTargetLength;
  }
  my $innerPrimerMaxLength =
    optionWithDefault($options_r, "inner_primer_max_length", 
      $optionDefaults{"inner_primer_max_length"});
  if($innerPrimerMaxLength < $innerPrimerTargetLength)
  {
    $innerPrimerMaxLength = $innerPrimerTargetLength;
  }

  my $innerPrimerTargetTM =
    optionWithDefault($options_r, "inner_primer_target_tm", 
      $optionDefaults{"inner_primer_target_tm"});
  my $innerPrimerMinTM =
    optionWithDefault($options_r, "inner_primer_min_tm", 
      ($innerPrimerTargetTM - 1.0));
  if($innerPrimerMinTM > $innerPrimerTargetTM)
  {
    $innerPrimerMinTM = $innerPrimerTargetTM;
  }
  my $innerPrimerMaxTM =
    optionWithDefault($options_r, "inner_primer_max_tm", 
      ($innerPrimerTargetTM + 1.0));
  if($innerPrimerMaxTM < $innerPrimerTargetTM)
  {
    $innerPrimerMaxTM = $innerPrimerTargetTM;
  }

  my $maxPolyBases = 
    optionWithDefault($options_r, "max_poly_bases", 
      $optionDefaults{"max_poly_bases"});
  
  my $includeLoopPrimers = optionWithDefault($options_r, "include_loop_primers",
    $optionDefaults{"include_loop_primers"});
  my $loopMinGap = 
    optionWithDefault($options_r, "loop_min_gap", 
      $optionDefaults{"loop_min_gap"});
  my $minSignaturesForSuccess =
    optionWithDefault($options_r, "min_signatures_for_success",
      $optionDefaults{"min_signatures_for_success"});
  my $minPrimerSpacing = 
    optionWithDefault($options_r, "min_primer_spacing", 
      $optionDefaults{"min_primer_spacing"});
  my $minInnerPairSpacing =
    optionWithDefault($options_r, "min_inner_pair_spacing", 
      $optionDefaults{"min_inner_pair_spacing"});
  #print "Loop min gap: $loopMinGap\n";
  #print "Max poly: $maxPolyBases\n";

  my $outerPairTargetLength = 
    optionWithDefault($options_r, "outer_pair_target_length", 
      $optionDefaults{"outer_pair_target_length"});
  my $middlePairTargetLength = 
    optionWithDefault($options_r, "middle_pair_target_length", 
      $optionDefaults{"middle_pair_target_length"});
  my $innerPairTargetLength =
    optionWithDefault($options_r, "inner_pair_target_length", 
      $optionDefaults{"inner_pair_target_length"});

  # Eventually want to let the user specify which penalty method
  # is used to calculate the spacing penalty, making the objective function
  # more customizable

  my $primer3ExecutablePath = optionWithDefault($options_r, "primer3_executable",
    $optionDefaults{"primer3_executable"});
  my $alignmehtFormat = optionWithDefault($options_r, "alignment_format",
    $optionDefaults{"alignment_format"});

  # In theory, the overall score logic belongs in a PrimerSetAnalyzer, 
  # but I hope this helps me optimize the inner loop implementing it
  # here, and only instantiating LAMP signatures for the best combinations
  my $innerPenaltyWeight = "1.2";
  my $loopPenaltyWeight = ".7";
  my $middlePenaltyWeight = "1.1";
  my $outerPenaltyWeight = "1.0";

  my $innerToLoopPenaltyWeight = 1.0;
  my $loopToMiddlePenaltyWeight = 1.0;
  my $innerToMiddlePenaltyWeight = 1.0;
  my $middleToOuterPenaltyWeight = 1.0;
  my $innerForwardToReversePenaltyWeight = 1.0;

  # Let the games begin...

  # Load the input alignment, could be a single sequence
  # TODO: # Make sure the alignment format option suggestion is working
  my $alignIN = Bio::AlignIO->new(-file => "< $alignmentFastaName");
  my $inputMSA = $alignIN->next_aln();
  my $sequenceLength = $inputMSA->length;

  # Ideally we would  have separate forward and reverse primer generation,
  # But since Primer3 doesn't accept "PRIMER_INTERNAL_OLIGO_MAX_STABILITY", 
  # we're going to have to filter that out ourselves, but it does mean that we can
  # cheat and just reverse complement the forward primers to get the reverse primers.
  my $maxEnumeratedPrimers = 20001; # Off-by-one in primer3?

  # Enumerate outer primers
  my $outerEnumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => $primer3ExecutablePath,
    });
  $outerEnumerator->setPrimer3Targets(
    {
      "target_length" => $outerPrimerTargetLength,
      "min_length" => $outerPrimerMinLength,
      "max_length" => $outerPrimerMaxLength,
      "target_tm" => $outerPrimerTargetTM,
      "min_tm" => $outerPrimerMinTM,
      "max_tm" => $outerPrimerMaxTM,
      "max_poly_bases" => $maxPolyBases,
      "most_to_return" => $maxEnumeratedPrimers,
    });

  print "Enumerating outer forward primers\n";
  my @outerForwardPrimers = $outerEnumerator->getOligos($inputMSA);

  print "  Generated \"" .
    scalar(@outerForwardPrimers) .
    "\" outer forward primers\n";


  print "Building outer reverse primers from outer forward primers\n";
  my @outerReversePrimers = buildReversePrimers(\@outerForwardPrimers);


  # Enumerate loop primers, since the loop primers extend in the opposite 
  # direction of the other LAMP primers, the back-loop primers are 
  # generated on the as-is sequence, and the forward-loop primers are 
  # built in the opposite orientation
  my $loopEnumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => $primer3ExecutablePath,
    });
  $loopEnumerator->setPrimer3Targets(
    {
      "target_length" => $loopPrimerTargetLength,
      "min_length" => $loopPrimerMinLength,
      "max_length" => $loopPrimerMaxLength,
      "target_tm" => $loopPrimerTargetTM,
      "min_tm" => $loopPrimerMinTM,
      "max_tm" => $loopPrimerMaxTM,
      "max_poly_bases" => $maxPolyBases,
      "most_to_return" => $maxEnumeratedPrimers,
    });

  # This difference in naming is intentional for now (loopBackPrimers instead of 
  # loopReversePrimers), to serve as a reminder that
  # loop primers extend the other direction, and that their locations need to be 
  # with the opposite orientation
  print "Enumerating loop \"back\" primers\n";
  my @loopBackPrimers = $loopEnumerator->getOligos($inputMSA);

  print "  Generated \"" .
    scalar(@loopBackPrimers) .
    "\" loop \"back\" primers\n";

  print "Building loop \"forward\" primers from loop \"back\" primers\n";
  my @loopForwardPrimers = buildReversePrimers(\@loopBackPrimers);

  # Enumerate middle primers
  my $middleEnumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => $primer3ExecutablePath,
    });
  $middleEnumerator->setPrimer3Targets(
    {
      "target_length" => $middlePrimerTargetLength,
      "min_length" => $middlePrimerMinLength,
      "max_length" => $middlePrimerMaxLength,
      "target_tm" => $middlePrimerTargetTM,
      "min_tm" => $middlePrimerMinTM,
      "max_tm" => $middlePrimerMaxTM,
      "max_poly_bases" => $maxPolyBases,
      "most_to_return" => $maxEnumeratedPrimers,
    });

  print "Enumerating middle forward primers\n";
  my @middleForwardPrimers = $middleEnumerator->getOligos($inputMSA);

  print "  Generated \"" .
    scalar(@middleForwardPrimers) .
    "\" middle primers\n";

  # Keeping minus strand primers in the same sorted order, which means
  # they're sorted by the 3' end of the primers, although that's the
  # 5' end of the primer span on the positive strand
  print "Building middle reverse primers from middle forward primers\n";
  my @middleReversePrimers = buildReversePrimers(\@middleForwardPrimers);

  # Enumerate inner primers 
  my $innerEnumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => $primer3ExecutablePath,
    });
  $innerEnumerator->setPrimer3Targets(
    {
      "target_length" => $innerPrimerTargetLength,
      "min_length" => $innerPrimerMinLength,
      "max_length" => $innerPrimerMaxLength,
      "target_tm" => $innerPrimerTargetTM,
      "min_tm" => $innerPrimerMinTM,
      "max_tm" => $innerPrimerMaxTM,
      "max_poly_bases" => $maxPolyBases,
      "most_to_return" => $maxEnumeratedPrimers,
    });

  print "Enumerating inner forward primers\n";
  my @innerForwardPrimers = $innerEnumerator->getOligos($inputMSA);

  print "  Generated \"" .
    scalar(@innerForwardPrimers) .
    "\" inner primers\n";

  # Keeping minus strand primers in the same sorted order, which means
  # they're sorted by the 3' end of the primers, although that's the
  # 5' end of the primer span on the positive strand
  print "Building inner reverse primers from inner forward primers\n";
  my @innerReversePrimers = buildReversePrimers(\@innerForwardPrimers);

  # TODO: want to flip any primer locations to reflect the standard
  # positive strand 5' location notation if they were generated
  # on an anti-sense strand, so all the length-based calculations
  # are handled only here, and the locations are standardized for the
  # rest of the process.

  # Analyze every oligo to get oligo penalty scores
  # Currently sharing one default analyzer for all the primers
  my $outerPrimerAnalyzer = LLNL::LAVA::PrimerAnalyzer::PCRPrimer->new();
  my $middlePrimerAnalyzer = $outerPrimerAnalyzer;
  my $innerPrimerAnalyzer = $outerPrimerAnalyzer;
  my $loopPrimerAnalyzer = $outerPrimerAnalyzer;

  print "Analyzing outer forward primers\n";
  my $outerForwardPrimerMeasurements_r =
    analyzeAll(\@outerForwardPrimers, $outerPrimerAnalyzer);
  print "Analyzing outer reverse primers\n";
  my $outerReversePrimerMeasurements_r =
    analyzeAll(\@outerReversePrimers, $outerPrimerAnalyzer);

  print "Analyzing loop \"forward\" primers\n";
  my $loopForwardPrimerMeasurements_r =
    analyzeAll(\@loopForwardPrimers, $loopPrimerAnalyzer);
  print "Analyzing loop \"back\" primers\n";
  my $loopBackPrimerMeasurements_r =
    analyzeAll(\@loopBackPrimers, $loopPrimerAnalyzer);

  print "Analyzing middle forward primers\n";
  my $middleForwardPrimerMeasurements_r = 
    analyzeAll(\@middleForwardPrimers, $middlePrimerAnalyzer);
  print "Analyzing middle reverse primers\n";
  my $middleReversePrimerMeasurements_r = 
    analyzeAll(\@middleReversePrimers, $middlePrimerAnalyzer);

  print "Analyzing inner forward primers\n";
  my $innerForwardPrimerMeasurements_r = 
    analyzeAll(\@innerForwardPrimers, $innerPrimerAnalyzer);
  print "Analyzing inner reverse primers\n";
  my $innerReversePrimerMeasurements_r = 
    analyzeAll(\@innerReversePrimers, $innerPrimerAnalyzer);

  # Sort all primers by 5' start location, and separately by score
  # It's tempting to rely on their current order, but I want to make
  # double sure we get increasing penalty sorting, so I'll do it explicitly.
  print "Sorting primer sets\n";

  # Not using an identifier to cross-reference between the sets, because
  # each location+length pair should be unique
  
  # Outer primers sorted 2 ways
  my @outerForwardInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$outerForwardPrimerMeasurements_r};
  my @outerReverseInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$outerReversePrimerMeasurements_r};

  my @outerForwardInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$outerForwardPrimerMeasurements_r};
  my @outerReverseInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$outerReversePrimerMeasurements_r};

  # Loop primers sorted 2 ways
  my @loopForwardInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$loopForwardPrimerMeasurements_r};
  my @loopBackInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$loopBackPrimerMeasurements_r};

  my @loopForwardInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$loopForwardPrimerMeasurements_r};
  my @loopBackInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$loopBackPrimerMeasurements_r};

  # Middle primers sorted 2 ways
  my @middleForwardInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$middleForwardPrimerMeasurements_r};
  my @middleReverseInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$middleReversePrimerMeasurements_r};

  my @middleForwardInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$middleForwardPrimerMeasurements_r};
  my @middleReverseInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$middleReversePrimerMeasurements_r};

  # Inner primers sorted 2 ways
  my @innerForwardInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$innerForwardPrimerMeasurements_r};
  my @innerReverseInfoByLocation =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @{$innerReversePrimerMeasurements_r};

  my @innerForwardInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$innerForwardPrimerMeasurements_r};
  my @innerReverseInfoByPenalty =
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getPenalty()] } 
    @{$innerReversePrimerMeasurements_r};

  #  print "Primers available:\n" .
  #    "  " . scalar(@{$innerForwardPrimerMeasurements_r}) . " Inner Forward\n" .
  #    "  " . scalar(@{$innerReversePrimerMeasurements_r}) . " Inner Reverse\n" .
  #    "  " . scalar(@{$loopForwardPrimerMeasurements_r}) . " Loop Forward\n" .
  #    "  " . scalar(@{$loopBackPrimerMeasurements_r}) . " Loop Reverse\n" .
  #    "  " . scalar(@{$middleForwardPrimerMeasurements_r}) . " Middle Forward\n" .
  #    "  " . scalar(@{$middleReversePrimerMeasurements_r}) . " Middle Reverse\n" .
  #    "  " . scalar(@{$outerForwardPrimerMeasurements_r}) . " Outer Forward\n" .
  #    "  " . scalar(@{$outerReversePrimerMeasurements_r}) . " Outer Reverse\n\n";

  print "Enumerating signatures\n";

  # Attempts will be made for combinations of different reduced primer sets.
  # The order attempts are made depends on the plan below.
  # 
  # Subgroups are built from the possible primer sets,  based on the overlap percent
  # specified in the subgroup schedule

  # Subgroup Schedule
  # Lookup of subgroup ID to array of cutoff percentages to use for the subgroup
  # Columns = Outer, Middle, Loop, Inner sets
  # Data cells are percentage overlap permitted among primers
  # (lower number  = fewer primers survive filtering)

  # This schedule works and is closely matching to Eiken
  my %subgroupSchedule =
    (
      1 => [50, 40, 50, 60], # Subgroup 1 cutoffs 
      2 => [70, 60, 70, 80], # Subgroup 1 cutoffs 
      3 => [80, 75, 80, 90], # Subgroup 2 cutoffs...
      4 => [90, 90, 90, 94], 
      #4 => [100, 100, 100, 100], 
    );
    #(
    #  1 => [70, 60, 70, 80], # Subgroup 1 cutoffs 
    #  2 => [80, 75, 80, 90], # Subgroup 2 cutoffs...
    #  3 => [90, 90, 90, 94], 
    #  4 => [100, 100, 100, 100], 
    #);

  # Columns for this plan definition are Outer, Middle, Loop, Inner subgroup IDs
  # Rows are sets of subgroup IDs to use for an iteration.  
  # Currently subgroup IDs are restricted to only numbers 1-4.
  # Plan is executed one row at a time, until a signature is found.

  # Looks like middle primer is the real floodgate control, although
  # loop would be preferrably...
  my @combinationPlan = 
    (
      [1, 1, 1, 1],
      [1, 1, 1, 2],
      [1, 1, 2, 2],
      [1, 2, 1, 1],
      [1, 2, 1, 2],
      [1, 2, 2, 2],
      [2, 2, 2, 3],
      [2, 2, 3, 3],
      [2, 2, 3, 4],
      [3, 2, 3, 3],
      [3, 3, 3, 3],
      [3, 3, 3, 4],
#      [4, 3, 4, 4],
#      [4, 4, 4, 4],
    ); 
#  my @combinationPlan = 
#    (
#      [1, 1, 1, 1],
#      [2, 2, 2, 1],
#      [2, 2, 2, 2],
#      [3, 3, 3, 1],
#      [3, 3, 3, 2],
#      [3, 3, 3, 3],
#      [4, 4, 4, 1],
#      [4, 4, 4, 2],
#      [4, 4, 4, 3],
#      [4, 4, 4, 4],
#    ); 

  # Plan execution, first steps are to build the subgroups needed. 
  my $possibleSignatures_r = [];
  my %cachedSubgroups = ();
  my %cachedSubgroupData = (); # Matching structure for %cachedSubgroups, holds array data
  foreach my $comboPlanStep_r (@combinationPlan)
  {
    my ($outerGroupID, $middleGroupID, $loopGroupID, $innerGroupID) = 
      @{$comboPlanStep_r};

    # Make sure each step of the combo plan gets a clean slate of signatures
    $possibleSignatures_r = [];

    # Create the subgroups if they don't already exist
    my $innerFName = "InnerF$innerGroupID";
    my $innerRName = "InnerR$innerGroupID";
    my $loopFName = "LoopF$loopGroupID";
    my $loopRName = "LoopR$loopGroupID";
    my $middleFName = "MiddleF$middleGroupID";
    my $middleRName = "MiddleR$middleGroupID";
    my $outerFName = "OuterF$outerGroupID";
    my $outerRName = "OuterR$outerGroupID";

    print "Iterating with groups: $innerFName, $innerRName, $loopFName, $loopRName, " .
      "$middleFName, $middleRName, $outerFName, $outerRName\n";

    if(! exists($cachedSubgroups{$innerFName}))
    {
      # Look up the percentage to filter by.
      my $maxOverlapPercent = $subgroupSchedule{$innerGroupID}->[3]; # 3 = inner index

      # Send the two separately sorted lists to reduce with
      print "Building $innerFName (at $maxOverlapPercent, ";
      $cachedSubgroups{$innerFName} = 
        reducePrimersByOverlap(
	  {
	    "max_overlap_percent" => $maxOverlapPercent,
	    "info_sorted_by_location"=> \@innerForwardInfoByLocation, 
	    "info_sorted_by_penalty" => \@innerForwardInfoByPenalty,
	  });
      $cachedSubgroupData{$innerFName} =
        flattenInfoData(
	  {
	    "info_set_ref" => $cachedSubgroups{$innerFName},
	  }); 
    }
    if(! exists($cachedSubgroups{$innerRName}))
    {
      # Look up the percentage to filter by.
      my $maxOverlapPercent = $subgroupSchedule{$innerGroupID}->[3]; # 3 = inner index

      # Send the two separately sorted lists to reduce with
      print "Building $innerRName (at $maxOverlapPercent, ";
      $cachedSubgroups{$innerRName} = 
        reducePrimersByOverlap(
	  {
	    "max_overlap_percent" => $maxOverlapPercent,
	    "info_sorted_by_location" => \@innerReverseInfoByLocation, 
	    "info_sorted_by_penalty" => \@innerReverseInfoByPenalty,
	  });
      $cachedSubgroupData{$innerRName} =
        flattenInfoData(
	  {
	    "info_set_ref" => $cachedSubgroups{$innerRName},
	  }); 
    }
    if(! exists($cachedSubgroups{$loopFName}))
    {
      # Look up the percentage to filter by.
      my $maxOverlapPercent = $subgroupSchedule{$loopGroupID}->[2]; # 2 = loop index

      # Send the two separately sorted lists to reduce with
      print "Building $loopFName (at $maxOverlapPercent";
      $cachedSubgroups{$loopFName} = 
        reducePrimersByOverlap(
	  {
	    "max_overlap_percent" => $maxOverlapPercent, 
	    "info_sorted_by_location" => \@loopForwardInfoByLocation, 
	    "info_sorted_by_penalty" => \@loopForwardInfoByPenalty,
	  });
      $cachedSubgroupData{$loopFName} =
        flattenInfoData(
	  {
	    "info_set_ref" => $cachedSubgroups{$loopFName},
	  }); 
    }
    if(! exists($cachedSubgroups{$loopRName}))
    {
      # Look up the percentage to filter by.
      my $maxOverlapPercent = $subgroupSchedule{$loopGroupID}->[2]; # 2 = loop index

      # Send the two separately sorted lists to reduce with
      print "Building $loopRName (at $maxOverlapPercent";
      $cachedSubgroups{$loopRName} = 
        reducePrimersByOverlap(
	  {
	    "max_overlap_percent" => $maxOverlapPercent, 
	    "info_sorted_by_location" => \@loopBackInfoByLocation, 
	    "info_sorted_by_penalty" => \@loopBackInfoByPenalty,
	  });
      $cachedSubgroupData{$loopRName} =
        flattenInfoData(
	  {
	    "info_set_ref" => $cachedSubgroups{$loopRName},
	  }); 
    }
    if(! exists($cachedSubgroups{$middleFName}))
    {
      # Look up the percentage to filter by.
      my $maxOverlapPercent = $subgroupSchedule{$middleGroupID}->[1]; # 1 = middle index

      # Send the two separately sorted lists to reduce with
      print "Building $middleFName (at $maxOverlapPercent";
      $cachedSubgroups{$middleFName} = 
        reducePrimersByOverlap(
	  {
	    "max_overlap_percent" => $maxOverlapPercent, 
	    "info_sorted_by_location" => \@middleForwardInfoByLocation, 
	    "info_sorted_by_penalty" => \@middleForwardInfoByPenalty,
	  });
      $cachedSubgroupData{$middleFName} =
        flattenInfoData(
	  {
	    "info_set_ref" => $cachedSubgroups{$middleFName},
	  }); 
    }
    if(! exists($cachedSubgroups{$middleRName}))
    {
      # Look up the percentage to filter by.
      my $maxOverlapPercent = $subgroupSchedule{$middleGroupID}->[1]; # 1 = middle index

      # Send the two separately sorted lists to reduce with
      print "Building $middleRName (at $maxOverlapPercent";
      $cachedSubgroups{$middleRName} = 
        reducePrimersByOverlap(
	  {
	    "max_overlap_percent" => $maxOverlapPercent, 
	    "info_sorted_by_location" => \@middleReverseInfoByLocation, 
	    "info_sorted_by_penalty" => \@middleReverseInfoByPenalty,
	  });
      $cachedSubgroupData{$middleRName} =
        flattenInfoData(
	  {
	    "info_set_ref" => $cachedSubgroups{$middleRName},
	  }); 
    }
    if(! exists($cachedSubgroups{$outerFName}))
    {
      # Look up the percentage to filter by.
      my $maxOverlapPercent = $subgroupSchedule{$outerGroupID}->[0]; # 0 = outer

      # Send the two separately sorted lists to reduce with
      print "Building $outerFName (at $maxOverlapPercent";
      $cachedSubgroups{$outerFName} = 
        reducePrimersByOverlap(
	  {
	    "max_overlap_percent" => $maxOverlapPercent, 
	    "info_sorted_by_location" => \@outerForwardInfoByLocation, 
	    "info_sorted_by_penalty" => \@outerForwardInfoByPenalty,
	  });
      $cachedSubgroupData{$outerFName} =
        flattenInfoData(
	  {
	    "info_set_ref" => $cachedSubgroups{$outerFName},
	  }); 
    }
    if(! exists($cachedSubgroups{$outerRName}))
    {
      # Look up the percentage to filter by.
      my $maxOverlapPercent = $subgroupSchedule{$outerGroupID}->[0]; # 0 = outer

      # Send the two separately sorted lists to reduce with
      print "Building $outerRName (at $maxOverlapPercent";
      $cachedSubgroups{$outerRName} = 
        reducePrimersByOverlap(
	  {
	    "max_overlap_percent" => $maxOverlapPercent, 
	    "info_sorted_by_location" => \@outerReverseInfoByLocation, 
	    "info_sorted_by_penalty" => \@outerReverseInfoByPenalty,
	  });
      $cachedSubgroupData{$outerRName} =
        flattenInfoData(
	  {
	    "info_set_ref" => $cachedSubgroups{$outerRName},
	  }); 
    }

    my $innerForwardSubset_r = $cachedSubgroups{$innerFName};
    my $innerForwardSubsetData_r = $cachedSubgroupData{$innerFName};
    my $innerReverseSubset_r = $cachedSubgroups{$innerRName};
    my $innerReverseSubsetData_r = $cachedSubgroupData{$innerRName};
    my $loopForwardSubset_r = $cachedSubgroups{$loopFName};
    my $loopForwardSubsetData_r = $cachedSubgroupData{$loopFName};
    my $loopReverseSubset_r = $cachedSubgroups{$loopRName};
    my $loopReverseSubsetData_r = $cachedSubgroupData{$loopRName};
    my $middleForwardSubset_r = $cachedSubgroups{$middleFName};
    my $middleForwardSubsetData_r = $cachedSubgroupData{$middleFName};
    my $middleReverseSubset_r = $cachedSubgroups{$middleRName};
    my $middleReverseSubsetData_r = $cachedSubgroupData{$middleRName};
    my $outerForwardSubset_r = $cachedSubgroups{$outerFName};
    my $outerForwardSubsetData_r = $cachedSubgroupData{$outerFName};
    my $outerReverseSubset_r = $cachedSubgroups{$outerRName};
    my $outerReverseSubsetData_r = $cachedSubgroupData{$outerRName};

################################################################################
# Now that the subgroups are picked, try to create signatures for this set of 
# combination of subgroups
################################################################################

    print "Primer counts used for this plan iteration:\n" .
      "  " . scalar(@{$innerForwardSubset_r}) . " Inner Forward\n" .
      "  " . scalar(@{$innerReverseSubset_r}) . " Inner Reverse\n" .
      "  " . scalar(@{$loopForwardSubset_r}) . " Loop Forward\n" .
      "  " . scalar(@{$loopReverseSubset_r}) . " Loop Reverse\n" .
      "  " . scalar(@{$middleForwardSubset_r}) . " Middle Forward\n" .
      "  " . scalar(@{$middleReverseSubset_r}) . " Middle Reverse\n" .
      "  " . scalar(@{$outerForwardSubset_r}) . " Outer Forward\n" .
      "  " . scalar(@{$outerReverseSubset_r}) . " Outer Reverse\n\n";

    
    my $innerForwardCount = scalar(@{$innerForwardSubset_r});
    my $innerReverseCount = scalar(@{$innerReverseSubset_r});
    my $loopForwardCount = scalar(@{$loopForwardSubset_r});
    my $loopReverseCount = scalar(@{$loopReverseSubset_r});
    my $middleForwardCount = scalar(@{$middleForwardSubset_r});
    my $middleReverseCount = scalar(@{$middleReverseSubset_r});
    my $outerForwardCount = scalar(@{$outerForwardSubset_r});
    my $outerReverseCount = scalar(@{$outerReverseSubset_r});


    # TODO: HERE! (should)
    # Pre-compute top Middle->Outer pairings, so we don't keep iterating over them
    #
    # Can't pre-compute others though, because loop primer may or may not be needed?
    # (Maybe conditionally compute other sets...?)

    #my $innerPairCount = scalar(@{$innerPairs_r});
    my @bestSignatureForInnerForward = ();
    #my $bestSignatureCount = 0;

    # Pre-compute a set of distance penalties for faster use
    # Let's set up maxSigLength bases worth of penalty, should be more than enough
    my $distancePenalties_r = generateDistancePenalties($signatureMaxLength); 

    # To remember the optimum combination with
    # 3 columns: loop, middle, outer
    my @bestForwardInfos = (); 
    # 2 columns: spacing_penalty, primer3_penalty
    my @bestForwardPenalties = ();

    # To help short-cut when no possibilities are found
    my $forwardSetCount = 0;

    for(my $innerIndex = 0; $innerIndex < $innerForwardCount; $innerIndex++)
    {
      my $innerInfo = $innerForwardSubset_r->[$innerIndex];
      my ($innerLocation, $innerLength, $innerPenalty) = 
        @{$innerForwardSubsetData_r->[$innerIndex]};

      my $bestSetPenalty = 1000000; # Riduculously large starting value

      # Calculate the first and last base locations to consider for the 
      # loop forward primer (on the plus strand darn it... need an inversion)
      my $searchStartAt = $innerLocation - $signatureMaxLength +
        $innerLength + 20; # 20 represents 2 other primer min lengths...?
      if($searchStartAt < 0)
      {
	$searchStartAt = 0;
      }

      # Calculates opposite strand position off by one to not reuse the last base	
      # Loop Start At and End At are indexes where the loop SEARCH starts and ends
      my $loopStartAt = $searchStartAt;
      my $loopEndAt = $innerLocation - 1 - $minPrimerSpacing;
      if($loopEndAt < 0)
      {
	$loopEndAt = 0;
      }

      print "\nNew inner\n";
      #print "(PA $sequenceLength long, start at $searchStartAt, maxLen $signatureMaxLength, $innerLocation->$loopStartAt-$loopEndAt)";
      #print "(PA* $loopStartAt-$loopEndAt)";

      # If no loop primers sought, then overwrite the loop primer list with the
      # single placeholder, one-length (but ideally zero-length), zero-penalty 
      # loop primer, placed at the end of the inner primer, to make sure it 
      # appears to fit within the acceptable locations
      if($includeLoopPrimers == $FALSE)
      {
	#print "\n\n\nNO LOOPS?!\n\n\n";
        my $placeHolderPrimer = LLNL::LAVA::Oligo->new(
	  {
            "sequence" => "N",
	    "location" => $loopEndAt + 1, # Extra position to un-do length of 1
	    "strand" => "minus",
	  });
        $placeHolderPrimer->setTag("primer3_penalty", 0);
        $placeHolderPrimer->setTag("primer3_tm", 0);

        my $placeHolderInfo = LLNL::LAVA::PrimerInfo->new(
	  {
            "penalty" => 0,
	    "sequence" => $placeHolderPrimer->sequence(),
	    "location" =>$placeHolderPrimer->location(),
	    "length" => $placeHolderPrimer->length(),
	    "analyzed_primer" => $placeHolderPrimer,
	  })

        $loopForwardSubset_r = [$placeHolderInfo];
	$loopForwardCount = 1; 
      }

      # Start of the 3-level nested loop for forward primers.
      # Should  exhaustively iterate over loop, middle, outer 
      # combinations based on the inner pair
      for(my $i = 0; $i < $loopForwardCount; $i++)
      {
	my $loopInfo = $loopForwardSubset_r->[$i];
	my ($loopLocation, $loopLength, $loopPenalty) = 
	  @{$loopForwardSubsetData_r->[$i]};
        #my $loopLocation = $loopInfo->getLocation();
        #my $loopLength = $loopInfo->getLength();

	# Special inversion for loop primer, because forward loop primer was
	# designed on the minus strand.
	
        # Seek to the first loop primer within range
	# but, accept placeholder loop primer
        if($loopLocation < $loopStartAt &&
	   $loopLength != 1)
	{
	  next;
	}

	# Stop when loop primer goes out of range	
	if($loopLocation > $loopEndAt &&
	   $loopLength != 1)
	{
	  last;
	}

        print "\nL";

        # No check for loop->inner overlap, becausel loopEndAt is the location limiter

	# First legitimate middle forward location is EITHER the first
	# that fits past the loop primer, OR the first that fits past
	# the minimum loop length
        my $middleStartAt = $searchStartAt;
        my $middleEndAt = $loopLocation - $loopLength - $minPrimerSpacing;
	my $altMiddleEndAt = $innerLocation - ($loopMinGap + 1);
	if($altMiddleEndAt < $middleEndAt)
	{
	  $middleEndAt = $altMiddleEndAt;
	}  
	if($middleEndAt < 0)
	{
	  $middleEndAt = 0;
	}
        #print "(LA $middleStartAt-$middleEndAt)";
        # Inter-primer distance used later for calculating spacing penalty
        my $innerToLoopDistance = $innerLocation - ($loopLocation + 1);
        # Can cause a negative 1 distance

        for(my $j = 0; $j < $middleForwardCount; $j++)
	{
	  my $middleInfo = $middleForwardSubset_r->[$j];
	  my ($middleLocation, $middleLength, $middlePenalty) = 
	    @{$middleForwardSubsetData_r->[$j]};

          #my $middleLocation = $middleInfo->getLocation();
          #my $middleLength = $middleInfo->getLength();

	  # Seek to the first middle primer within range
	  if($middleLocation < $middleStartAt)
	  {
	    next;
	  }

	  # Stop when middle primer goes out of range
	  if($middleLocation > $middleEndAt)
	  {
	    last;
	  }

          print "M";

	  # Next primer if this middle doesn't leave enough spacing to the loop primer,
	  # or doesn't leave enough of the loop minimum gap to the inner primer
	  if(($middleLocation + $middleLength + $minPrimerSpacing > 
	      $loopLocation - $loopLength + 1) ||
	     ($middleLocation + $middleLength + $loopMinGap >
	      $innerLocation) )
	  {
	    next;
	  }

	  my $outerStartAt = $searchStartAt;
          my $outerEndAt = $middleLocation - 1 - $minPrimerSpacing;

         #print "(MA $outerStartAt-$outerEndAt)";

         # Inter-primer distance used later for calculating spacing penalty 
          my $loopToMiddleDistance = ($loopLocation - $loopLength + 1) - 
	    ($middleLocation + $middleLength);
	  
          for(my $k = 0; $k < $outerForwardCount; $k++)
	  {
	    my $outerInfo = $outerForwardSubset_r->[$k];
	    my ($outerLocation, $outerLength, $outerPenalty) = 
	      @{$outerForwardSubsetData_r->[$k]};
	    #my $outerLocation = $outerInfo->getLocation();
            #my $outerLength = $outerInfo->getLength();

            # Seek to first outer primer within range
	    if($outerLocation < $outerStartAt)
	    {
	      next;
	    }

	    # Stop when outer primer goes out of range
	    if($outerLocation > $outerEndAt)
	    {
              last;
	    }

            #print "O";
            # Next primer if this outer doesn't leave enough spacing to the middle primer
            if($outerLocation + $outerLength + $minPrimerSpacing >
	      $middleLocation)
	    {
	      next;
	    }
            #print "(OA)";

	    # Inter-primer distance used for calculating spacing penalty
            my $middleToOuterDistance = $middleLocation - 
	      ($outerLocation + $outerLength);
           
	    # The inner-to-middle distance is only used for non-loop primer
	    # signatures
	    my $innerToMiddleDistance = $innerLocation - 
	      ($middleLocation + $middleLength);
	 
	    my $spacingPenalty = 0;
	    my $primer3Penalty = 0;
	    if($includeLoopPrimers == $TRUE)
	    {
              $spacingPenalty = 
	        ($distancePenalties_r->[$innerToLoopDistance] * 
	         $innerToLoopPenaltyWeight) +
                ($distancePenalties_r->[$loopToMiddleDistance] *
		 $loopToMiddlePenaltyWeight) +
                ($distancePenalties_r->[$middleToOuterDistance] *
	         $middleToOuterPenaltyWeight);
	      $primer3Penalty = 
	        $innerPenalty * $innerPenaltyWeight +
	        $loopPenalty * $loopPenaltyWeight +
	        $middlePenalty * $middlePenaltyWeight +
	        $outerPenalty * $outerPenaltyWeight;
	    }
	    else
	    {
	      $spacingPenalty = 
	        ($distancePenalties_r->[$innerToMiddleDistance] * 
	         $innerToMiddlePenaltyWeight) +
                ($distancePenalties_r->[$middleToOuterDistance] *
	          $middleToOuterPenaltyWeight);
	      $primer3Penalty = 
	        $innerPenalty * $innerPenaltyWeight +
	        $middlePenalty * $middlePenaltyWeight +
	        $outerPenalty * $outerPenaltyWeight;
	    }
	   
	    my $forwardSetPenalty = $spacingPenalty + $primer3Penalty;
	    if($forwardSetPenalty < $bestSetPenalty)
	    {
	      $bestForwardInfos[$innerIndex] = [$loopInfo, $middleInfo, $outerInfo];
	      $bestSetPenalty = $forwardSetPenalty;
	      $forwardSetCount++;
              #print "*";             
	      # Not strictly necessary, but helpful for fine-tuning
              $bestForwardPenalties[$innerIndex] = [$spacingPenalty, $primer3Penalty];
	    }
	  } # End forward outer iteration
	} # End forward middle iteration
      } # End forward loop iteration
    } # End forward inner iteration

    # Stop trying if no forward primer sets were found
    if($forwardSetCount == 0)
    {
      #print "F";
      next;
    }

    my @bestSigantureForInnerReverse = ();

    # To remember the optimum combination with
    # 3 columns: loop, middle, outer
    my @bestReverseInfos = (); 
    # 2 columns: spacing_penalty, primer3_penalty
    my @bestReversePenalties = ();

    ## To help short-cut when no possibilities are found
    #my $reverseSetCount = 0;

    # Yes, this is the exact same logic as the above block, for some reason it was
    # easier for me to implement as two separate blocks, probably because
    # of the number of locations that had to be inverted to do it in a single block
    for(my $innerIndex = 0; $innerIndex < $innerReverseCount; $innerIndex++)
    {
      my $innerInfo = $innerReverseSubset_r->[$innerIndex];
      my ($innerLocation, $innerLength, $innerPenalty) = 
        @{$innerReverseSubsetData_r->[$innerIndex]};

      #my $innerLocation = $innerInfo->getLocation();
      #my $innerLength = $innerInfo->getLength();

      my $bestSetPenalty = 1000000; # Riduculously large starting value

      # Calculate the first and last base locations to consider for the 
      # loop reverse primer (on the plus strand, so need inversion again?)
      my $searchEndAt = $innerLocation + $signatureMaxLength - 
	$innerLength - 20; # -20 represents 2 other primer min lengths.
	
      # All the -1s are to deal with zero-indexing with counts that start at 1 
      my $loopStartAt = $innerLocation + 1 + $minPrimerSpacing;
      my $loopEndAt = $searchEndAt;

      # If no loop primers sought, then overwrite the loop primer list with the
      # single placeholder, one-length (but ideally zero-length), zero-penalty
      # loop primer, placed at the end of the inner primer, to make sure it 
      # appears to fit within the acceptable locations
      if($includeLoopPrimers == $FALSE)
      {
        my $placeHolderPrimer = LLNL::LAVA::Oligo->new(
	  {
            "sequence" => "N",
	    "location" => $loopStartAt - 1, # Extra position to un-do length of 1
	    "strand" => "plus",
	  });
        $placeHolderPrimer->setTag("primer3_penalty", 0);
        $placeHolderPrimer->setTag("primer3_tm", 0);

        my $placeHolderInfo = LLNL::LAVA::PrimerInfo->new(
	  { 
            "penalty" => 0,
	    "sequence" => $placeHolderPrimer->sequence(),
	    "location" =>$placeHolderPrimer->location(),
	    "length" => $placeHolderPrimer->length(),
	    "analyzed_primer" => $placeHolderPrimer,
	  });

        $loopReverseSubset_r = [$placeHolderInfo];
	$loopReverseCount = 1; 
      }

      # Start of the 3-level nested loop for reverse primers.
      # Should  exhaustively iterate over loop, middle, outer 
      # combinations based on the inner pair
      for(my $i = 0; $i < $loopReverseCount; $i++)
      {
	my $loopInfo = $loopReverseSubset_r->[$i];
	my ($loopLocation, $loopLength, $loopPenalty) = 
	  @{$loopReverseSubsetData_r->[$i]};

        #my $loopLocation = $loopInfo->getLocation();
        #my $loopLength = $loopInfo->getLength();

	# Special inversion for loop primer, because reverse loop primer was
	# designed on the minus strand.
	
        # Seek to the first loop primer within range 
	# but, accept placeholder loop primer
        if($loopLocation < $loopStartAt &&
	   $loopLength != 1)
	{
	  next;
	}

	# Stop when loop primer goes out of range	
	if($loopLocation > $loopEndAt &&
	   $loopLength != 1)
	{
	  last;
	}

        #print "L";

        # No check for loop->inner overlap because loopStartAt is the location limiter

	# First legitimate middle forward location is EITHER the first
	# that fits past the loop primer, OR the first that fits past
	# the minimum loop length
        my $middleStartAt = $loopLocation + $loopLength + $minPrimerSpacing;
	my $altMiddleStartAt = $innerLocation + ($loopMinGap + 1);
	if($altMiddleStartAt > $middleStartAt)
	{
	  $middleStartAt = $altMiddleStartAt;
	}  
        my $middleEndAt = $searchEndAt;

        # Inter-primer distance used later for calculating spacing penalty
	my $innerToLoopDistance = $loopLocation - ($innerLocation + 1);
        # Can cause a negative 1 distance

        for(my $j = 0; $j < $middleReverseCount; $j++)
	{
	  my $middleInfo = $middleReverseSubset_r->[$j];
	  my ($middleLocation, $middleLength, $middlePenalty) = 
	    @{$middleReverseSubsetData_r->[$j]};

          #my $middleLocation = $middleInfo->getLocation();
          #my $middleLength = $middleInfo->getLength();

	  # Seek to the first middle primer within range
	  if($middleLocation < $middleStartAt)
	  {
	    next;
	  }

	  # Stop when middle primer goes out of range
	  if($middleLocation > $middleEndAt)
	  {
	    last;
	  }

          #print "M";

          # Next primer if this middle doesn't leave enough spacing to the loop primer,
	  # or doesn't leave enough of the loop minimum gap to the inner primer
          if(($middleLocation - $middleLength - $minPrimerSpacing <
	      $loopLocation + $loopLength - 1) ||
             ($middleLocation - $middleLength - $loopMinGap <
	      $innerLocation) )
          {
            next;
	  }

          my $outerStartAt = $middleLocation + 1 + $minPrimerSpacing;
          my $outerEndAt = $searchEndAt;

          # Inter-primer distance used later for calculating spacing penalty 
          my $loopToMiddleDistance = ($middleLocation - $middleLength + 1) -
	    ($loopLocation + $loopLength);

          for(my $k = 0; $k < $outerReverseCount; $k++)
	  {
	    my $outerInfo = $outerReverseSubset_r->[$k];
	    my ($outerLocation, $outerLength, $outerPenalty) = 
	      @{$outerReverseSubsetData_r->[$k]};

	    #my $outerLocation = $outerInfo->getLocation();
            #my $outerLength = $outerInfo->getLength();

            # Seek to first outer primer within range
	    if($outerLocation < $outerStartAt)
	    {
	      next;
	    }

	    # Stop when outer primer goes out of range
	    if($outerLocation > $outerEndAt)
	    {
              last;
	    }

            #print "O";
            # Next primer if this outer doesn't leave enough spacing to the middle primer
	    if($outerLocation - $outerLength - $minPrimerSpacing <
	       $middleLocation)
	    {
	      next;
	    }

	    # Inter-primer distance used for calculating spacing penalty
            my $middleToOuterDistance = ($outerLocation - $outerLength) -
	      $middleLocation;
	       
	    # The inner-to-middle distance is only used for non-loop primer
	    # signatures
	    my $innerToMiddleDistance = ($middleLocation - $middleLength) -
	      $innerLocation;
	 
	    my $spacingPenalty = 0;
	    my $primer3Penalty = 0;
	    if($includeLoopPrimers == $TRUE)
	    {
              $spacingPenalty = 
	        ($distancePenalties_r->[$innerToLoopDistance] * 
	         $innerToLoopPenaltyWeight) +
                ($distancePenalties_r->[$loopToMiddleDistance] *
		 $loopToMiddlePenaltyWeight) +
                ($distancePenalties_r->[$middleToOuterDistance] *
	         $middleToOuterPenaltyWeight);
	      $primer3Penalty = 
	        $innerPenalty * $innerPenaltyWeight +
	        $loopPenalty * $loopPenaltyWeight +
	        $middlePenalty * $middlePenaltyWeight +
	        $outerPenalty * $outerPenaltyWeight;
	    }
	    else
	    {
	      $spacingPenalty = 
	        ($distancePenalties_r->[$innerToMiddleDistance] * 
	         $innerToMiddlePenaltyWeight) +
                ($distancePenalties_r->[$middleToOuterDistance] *
	          $middleToOuterPenaltyWeight);
	      $primer3Penalty = 
	        $innerPenalty * $innerPenaltyWeight +
	        $middlePenalty * $middlePenaltyWeight +
	        $outerPenalty * $outerPenaltyWeight;
	    }
 
	    my $reverseSetPenalty = $spacingPenalty + $primer3Penalty;
	    if($reverseSetPenalty < $bestSetPenalty)
	    {
	      $bestReverseInfos[$innerIndex] = [$loopInfo, $middleInfo, $outerInfo];
	      $bestSetPenalty = $reverseSetPenalty;
	      #$reverseSetCount++;

	      # Not strictly necessary, but helpful for fine-tuning
              $bestReversePenalties[$innerIndex] = [$spacingPenalty, $primer3Penalty];
	    }
	  } # End reverse outer iteration
	} # End reverse middle iteration
      } # End reverse loop iteration
    } # End reverse inner iteration

    ## Stop trying if no reverse primer sets were found (probably an un-needed optimization)
    #if($reverseSetCount == 0)
    #{
    #  print "R";
    #  next;
    #}

    # Now, try to combine forward and reverse primer sets into full signatures
    my $previousFirstCompatibleIndex = 0; # Bound the lower end of the inner loop
    for(my $i = 0; $i < $innerForwardCount; $i++)
    {
      # Skip inner primers without primer sets
      if(! exists($bestForwardInfos[$i]))
      {
	next;
      }

      my $finnerInfo = $innerForwardSubset_r->[$i];
      my ($floopInfo, $fmiddleInfo, $fouterInfo) = @{$bestForwardInfos[$i]};
      my ($forwardSpacingPenalty, $forwardPrimer3Penalty) = 
        @{$bestForwardPenalties[$i]};

      my $forwardStart = $fouterInfo->getLocation();
      my $forwardEnd = $finnerInfo->getLocation() + $finnerInfo->getLength() - 1;

      #print "\nInner $forwardStart -> $forwardEnd";

      # Used to bound the upper end of the inner loop search
      my $maxReverseLocation = $forwardStart + $signatureMaxLength - 1;
    
      # Used to help bound the lower end of the inner loop search
      my $previousCompatibleIndexFound = $FALSE;
      
      for(my $j = $previousFirstCompatibleIndex; $j < $innerReverseCount; $j++)
      {
        # Skip reverse primers without primer sets
        if(! exists($bestReverseInfos[$j]))
        {
	  next;
        }
	my $binnerInfo = $innerReverseSubset_r->[$j];
	my ($bloopInfo, $bmiddleInfo, $bouterInfo) = @{$bestReverseInfos[$j]};
        my ($reverseSpacingPenalty, $reversePrimer3Penalty) = 
	  @{$bestReversePenalties[$j]};

        my $reverseEnd = $bouterInfo->getLocation();
        my $reverseStart = $binnerInfo->getLocation() - $binnerInfo->getLength() + 1;
        #print "\n  Outer $reverseStart -> $reverseEnd";

        # Advance to the next compatible reverse primer by skipping all the
        # primers located too far 5' with respect to the forward primer
        if($previousCompatibleIndexFound == $FALSE)
        {
          if($reverseStart <= $forwardEnd)
          {
            next;
          }
          else
          {
            $previousFirstCompatibleIndex = $j;
            $previousCompatibleIndexFound = $TRUE;
          }
        }

	# Stop searching if the inner loop bounds are exceeded
        if($reverseStart > $maxReverseLocation)
        {
          last;
        }
 
        # Enforce minimum inner spacing distance
	my $innerSpacing = $reverseStart - ($forwardEnd + 1);
        if($innerSpacing < $minInnerPairSpacing)
	{
	  next;
	}
       
        # Enforce max signature length
        if($reverseEnd - ($forwardStart + 1) > $signatureMaxLength)
	{
	  next;
	}

        #print "*";
    
        # TODO: Spacing penalty should probably exclude minimum required spacings?
        my $innerSpacingPenalty = 
          ($distancePenalties_r->[$innerSpacing] *
	   $innerForwardToReversePenaltyWeight);

        my $totalPenalty = $forwardSpacingPenalty +
	  $innerSpacingPenalty +
	  $reverseSpacingPenalty +
	  $forwardPrimer3Penalty +
	  $reversePrimer3Penalty;

        my $innerPair = LLNL::LAVA::PrimerSet::PCRPair->new(
	  {
	    "forward_info" => $finnerInfo,
	    "reverse_info" => $binnerInfo,
	  });
        my $middlePair = LLNL::LAVA::PrimerSet::PCRPair->new(
	  {
	    "forward_info" => $fmiddleInfo,
	    "reverse_info" => $bmiddleInfo,
	  });
        my $outerPair = LLNL::LAVA::PrimerSet::PCRPair->new(
          {
            "forward_info" => $fouterInfo,
            "reverse_info" => $bouterInfo,
	  });
  
        my $innerPairInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new(
          {
            "penalty" => 0,
            "analyzed_pair" => $innerPair,
          });
        my $middlePairInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new(
          {
            "penalty" => 0,
            "analyzed_pair" => $middlePair,
          });
        my $outerPairInfo = LLNL::LAVA::PrimerSetInfo::PCRPair->new(
          {
            "penalty" => 0,
            "analyzed_pair" => $outerPair,
          });
  
        my $signature = LLNL::LAVA::PrimerSet::LAMP->new(
          {
            "inner_info" => $innerPairInfo,
            "middle_info" => $middlePairInfo,
            "outer_info" => $outerPairInfo,
          });
  
        $signature->setTag("lamp_penalty", $totalPenalty);
  
        # Just for fine-tuning
        $signature->setTag("penalty_notes", 
          "(spacing: $forwardSpacingPenalty, $innerSpacingPenalty, $reverseSpacingPenalty)" .
	  "(p3: $forwardPrimer3Penalty, $reversePrimer3Penalty)");
  
        if($includeLoopPrimers == $TRUE)
        {
          $signature->setTag("floop_info", $floopInfo);
          $signature->setTag("bloop_info", $bloopInfo);
          $signature->setTag("has_loop_primers", $TRUE);
        }
         
        push(@{$possibleSignatures_r}, $signature);
      } # End forward sets iteration
    } # End reverse sets iteration

    $possibleSignatures_r = reduceSignaturesByOverlap(
      {
	"signatures" => $possibleSignatures_r,
	"max_overlap_percent" => 0,
      });

    # Stop iterating over plans if we have enough signatures to be done with
    # TODO: Toss a useful warning if a high sig count is requested with only low
    # overlap allowed, or on a "small" region...
    my $signatureCount = scalar(@{$possibleSignatures_r});
    if($signatureCount > 0 && $signatureCount >= $minSignaturesForSuccess)
    {
      #print "G";
      last;
    }  
  } # End combination plan iteration

  # TODO: do something more reaosnable here
  if(scalar(@{$possibleSignatures_r}) == 0)
  {
    print "Faled to identify signatures - exiting normally\n";
    exit(0);
  }
  
  print "Found " .
    scalar(@{$possibleSignatures_r}) .
    " possible signatures\n";

  # Sort signatures by score
  my @possibleSignatures = 
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1] }
    map {[$_, $_->getTag("lamp_penalty")]}
    @{$possibleSignatures_r};

  # Write the output fasta
  #TODO: watch for stompping?
  open(OUTANSWER, "> $outputFileName") || 
    confess("file error - failed to open output file \"$outputFileName\" " .
      "for writing: $!");

  my $possibleSignatureCount = scalar(@possibleSignatures);
  for(my $i = 0; $i < $possibleSignatureCount; $i++)
  {
    my $signature = $possibleSignatures[$i];
    my $signatureName = "$i";
 
    # Set the linker to the dash, but restore the exsting linker afterwards
    my $originalLinker = $signature->linker();  
    $signature->linker("");

    my $penalty = $signature->getTag("lamp_penalty");
    my $locationSummary = $signature->getLocationSummary();
    #my $penaltySummary = $signature->getPenaltySummary(); 
    #my $tmSummary = $signature->getTMSummary();
     
    #print OUTANSWER ">$signatureName F3 (penatly: $penalty) (locations: $locationSummary) " .
    #  "(sub-penalties: $penaltySummary) (tms: $tmSummary)\n"; 
    # TODO: update the sig reader to load this data too! (need something more flexible!)

    #print OUTANSWER ">$signatureName F3 (penalty: $penalty) (locations: $locationSummary)\n";
    my $penaltyNotes = $signature->getTag("penalty_notes");
    print OUTANSWER ">$signatureName F3 (penalty: $penalty) $penaltyNotes (locations: $locationSummary)\n";
    print OUTANSWER $signature->getF3() . "\n";
    print OUTANSWER ">$signatureName B3\n"; 
    print OUTANSWER $signature->getB3() . "\n";
    print OUTANSWER ">$signatureName FIP\n"; 
    print OUTANSWER $signature->getFIP() . "\n";
    print OUTANSWER ">$signatureName BIP\n"; 
    print OUTANSWER $signature->getBIP() . "\n";
    if($includeLoopPrimers == $TRUE)
    {
      my $floopSequence = ($signature->getTag("floop_info"))->getSequence();
      my $bloopSequence = ($signature->getTag("bloop_info"))->getSequence();
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      print OUTANSWER ">$signatureName FLOOP (locations: $loopLocationSummary)\n";
      print OUTANSWER $floopSequence . "\n";
      print OUTANSWER ">$signatureName BLOOP\n";
      print OUTANSWER $bloopSequence . "\n";
    }
      
    # Return the linker back to its original state
    $signature->linker($originalLinker);
  }

  close(OUTANSWER) ||
    confess("file error - failed to cose output file \"$outputFileName\": $!");

  # Write the linked-marker file, using the dash-linker in context
  #TODO: watch for stompping?
  my $dashFileName = "$outputFileName.dash";
  open(OUTDASH, "> $dashFileName") || 
    confess("file error - failed to open output file \"$dashFileName\" " .
      "for writing: $!");

  for(my $i = 0; $i < $possibleSignatureCount; $i++)
  {
    my $signature = $possibleSignatures[$i];
    my $signatureName = "$i";
 
    # Set the linker to the dash, but restore the exsting linker afterwards
    my $originalLinker = $signature->linker();  
    $signature->linker("-");

    my $penalty = $signature->getTag("lamp_penalty");
    my $locationSummary = $signature->getLocationSummary();
    #my $penaltySummary = $signature->getPenaltySummary(); 
    #my $tmSummary = $signature->getTMSummary();
     
    #print OUTDASH ">$signatureName F3 (penatly: $penalty) (locations: $locationSummary) " .
    #  "(sub-penalties: $penaltySummary) (tms: $tmSummary)\n"; 
    # TODO: update the sig reader to load this data too! (need something more flexible!)

    #print OUTDASH ">$signatureName F3 (penalty: $penalty) (locations: $locationSummary)\n";
    my $penaltyNotes = $signature->getTag("penalty_notes");
    print OUTDASH ">$signatureName F3 (penalty: $penalty) $penaltyNotes (locations: $locationSummary)\n";
    print OUTDASH $signature->getF3() . "\n";
    print OUTDASH ">$signatureName B3\n"; 
    print OUTDASH $signature->getB3() . "\n";
    print OUTDASH ">$signatureName FIP\n"; 
    print OUTDASH $signature->getFIP() . "\n";
    print OUTDASH ">$signatureName BIP\n"; 
    print OUTDASH $signature->getBIP() . "\n";
    if($includeLoopPrimers == $TRUE)
    {
      my $floopSequence = ($signature->getTag("floop_info"))->getSequence();
      my $bloopSequence = ($signature->getTag("bloop_info"))->getSequence();
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      print OUTDASH ">$signatureName FLOOP (locations: $loopLocationSummary)\n";
      print OUTDASH $floopSequence . "\n";
      print OUTDASH ">$signatureName BLOOP\n";
      print OUTDASH $bloopSequence . "\n";
    }
      
    # Return the linker back to its original state
    $signature->linker($originalLinker);
  }

  close(OUTDASH) ||
    confess("file error - failed to cose output file \"$dashFileName\": $!");

  # Write the individual primers (in extension orientation) as an answer file
  #TODO: watch for stompping?
  my $primersFileName = "$outputFileName.primers";
  open(OUTPRIMERS, "> $primersFileName") || 
    confess("file error - failed to open output file \"$primersFileName\" " .
      "for writing: $!");

  for(my $i = 0; $i < $possibleSignatureCount; $i++)
  {
    my $signature = $possibleSignatures[$i];
 
    my $penalty = $signature->getTag("lamp_penalty");
    my $locationSummary = $signature->getLocationSummary();
    #my $penaltySummary = $signature->getPenaltySummary(); 
    #my $tmSummary = $signature->getTMSummary();
     
    #print OUTPRIMERS ">$signatureName F3 (penatly: $penalty) (locations: $locationSummary) " .
    #  "(sub-penalties: $penaltySummary) (tms: $tmSummary)\n"; 
    # TODO: update the sig reader to load this data too! (need something more flexible!)

    #print OUTPRIMERS ">$signatureName F3 (penalty: $penalty) (locations: $locationSummary)\n";
    my $penaltyNotes = $signature->getTag("penalty_notes");
    print OUTPRIMERS ">$i F3 (penalty: $penalty) $penaltyNotes (locations: $locationSummary)\n";
    print OUTPRIMERS $signature->getF3() . "\n";
    print OUTPRIMERS ">$i B3\n"; 
    print OUTPRIMERS $signature->getB3() . "\n";
    print OUTPRIMERS ">$i F2\n";
    print OUTPRIMERS $signature->getF2() . "\n";
    print OUTPRIMERS ">$i B2\n";
    print OUTPRIMERS $signature->getB2() . "\n";
    print OUTPRIMERS ">$i F1\n";
    print OUTPRIMERS $signature->getF1() . "\n";
    print OUTPRIMERS ">$i B1\n";
    print OUTPRIMERS $signature->getB1() . "\n";
    if($includeLoopPrimers == $TRUE)
    {
      my $floopSequence = ($signature->getTag("floop_info"))->getSequence();
      my $bloopSequence = ($signature->getTag("bloop_info"))->getSequence();
      my $loopLocationSummary = $signature->getLoopLocationSummary();
      print OUTPRIMERS ">$i FLOOP (locations: $loopLocationSummary)\n";
      print OUTPRIMERS $floopSequence . "\n";
      print OUTPRIMERS ">$i BLOOP\n";
      print OUTPRIMERS $bloopSequence . "\n";
    }
  }

  close(OUTPRIMERS) ||
    confess("file error - failed to cose output file \"$primersFileName\": $!");



##  # Dump the signature structures for debugging
##  open(DUMP, "> test_signature_dump.perl") ||
##    confess("failed to open debug output dump file: $!");
##  print DUMP Data::Dumper->Dump(\@signatures);
##  close(DUMP);
#
#  # Write the results file (default linker is the empty string - no linker)
#  # TODO: watch for stompping?
#  open(OUTRESULTS, "> $outputFileName") || 
#    confess("file error - failed to open output file \"$outputFileName\" " .
#      "for writing: $!");
#  for(my $i = 0; $i < $signatureCount; $i++)
#  {
#    my $signature = $signatures[$i];
#
#    my $penalty = $signature->getTag("lamp_penalty");
#    my $locationSummary = $signature->getLocationSummary();
#    my $penaltySummary = $signature->getPenaltySummary(); 
#    my $tmSummary = $signature->getTMSummary();
#    
#    print OUTRESULTS ">$i F3 (penatly: $penalty) (locations: $locationSummary) " .
#      "(sub-penalties: $penaltySummary) (tms: $tmSummary)\n"; 
#    print OUTRESULTS $signature->getF3() . "\n";
#    print OUTRESULTS ">$i B3\n"; 
#    print OUTRESULTS $signature->getB3() . "\n";
#    print OUTRESULTS ">$i FIP\n"; 
#    print OUTRESULTS $signature->getFIP() . "\n";
#    print OUTRESULTS ">$i BIP\n"; 
#    print OUTRESULTS $signature->getBIP() . "\n";
#  }
#
#  close(OUTRESULTS) ||
#    confess("file error - failed to cose output file \"$outputFileName\": $!");
#
#  # Write the linked-marker file, using the dash-linker in context
#  #TODO: watch for stompping?
#  my $dashFileName = "$outputFileName.dash";
#  open(OUTDASH, "> $dashFileName") || 
#    confess("file error - failed to open output file \"$dashFileName\" " .
#      "for writing: $!");
#  for(my $i = 0; $i < $signatureCount; $i++)
#  {
#    my $signature = $signatures[$i];
#
#    # Set the linker to the dash, but restore the exsting linker afterwards
#    my $originalLinker = $signature->linker();  
#    $signature->linker("-");
#
#    my $penalty = $signature->getTag("lamp_penalty");
#    my $locationSummary = $signature->getLocationSummary();
#    my $penaltySummary = $signature->getPenaltySummary(); 
#    my $tmSummary = $signature->getTMSummary();
#    
#    print OUTDASH ">$i F3 (penatly: $penalty) (locations: $locationSummary) " .
#      "(sub-penalties: $penaltySummary) (tms: $tmSummary)\n"; 
#    print OUTDASH $signature->getF3() . "\n";
#    print OUTDASH ">$i B3\n"; 
#    print OUTDASH $signature->getB3() . "\n";
#    print OUTDASH ">$i FIP\n"; 
#    print OUTDASH $signature->getFIP() . "\n";
#    print OUTDASH ">$i BIP\n"; 
#    print OUTDASH $signature->getBIP() . "\n";
#
#    $signature->linker($originalLinker);
#  }
#
#  close(OUTDASH) ||
#    confess("file error - failed to cose output file \"$dashFileName\": $!");
#
#  # Write the individual primers as an answer file - 6 per signature 
#  #TODO: watch for stompping?
#  my $primerFileName = "$outputFileName.primers";
#  open(OUTPRIMER, "> $primerFileName") || 
#    confess("file error - failed to open output file \"$primerFileName\" " .
#      "for writing: $!");
#  for(my $i = 0; $i < $signatureCount; $i++)
#  {
#    my $signature = $signatures[$i];
#
#    my $penalty = $signature->getTag("lamp_penalty");
#    my $locationSummary = $signature->getLocationSummary();
#    my $penaltySummary = $signature->getPenaltySummary(); 
#    my $tmSummary = $signature->getTMSummary();
#    
#    print OUTPRIMER ">$i F3 (penatly: $penalty) (locations: $locationSummary) " .
#      "(sub-penalties: $penaltySummary) (tms: $tmSummary)\n"; 
#    print OUTPRIMER $signature->getF3() . "\n";
#    print OUTPRIMER ">$i B3\n"; 
#    print OUTPRIMER $signature->getB3() . "\n";
#    print OUTPRIMER ">$i F2\n";
#    print OUTPRIMER $signature->getF2() . "\n";
#    print OUTPRIMER ">$i B2\n";
#    print OUTPRIMER $signature->getB2() . "\n";
#    print OUTPRIMER ">$i F1\n";
#    print OUTPRIMER $signature->getF1() . "\n";
#    print OUTPRIMER ">$i B1\n";
#    print OUTPRIMER $signature->getB1() . "\n";
#  }
#
#  close(OUTPRIMER) ||
#    confess("file error - failed to cose output file \"$primerFileName\": $!");

  print "Exiting normally\n";
}

# This doesn't seem needed if we don't generate primers on the opposite strand.
# If we DO start using it, make sure locations get inverted on conversion back
# to plus-strand locations, because they will be backwards!
## Helper function, accepts a Bio::Align::AlignI (from Bio::AlignIO), and
## creates a reverse-complement version of the alignment by converting
## each individual sequence in the alignment
#sub reverseAlignmentStrand
#{
#  my ($alignment) = @_;
#
#  my $sequenceCount = $alignment->no_sequences();
#  if($sequenceCount <= 0)
#  {
#    confess("data error - MSA must contain at least one sequence");
#  }
#
#  my $newAlignment = Bio::SimpleAlign->new();
#  for(my $i = 1; $i < $sequenceCount + 1 ; $i++) # 1- indexed for BioPerl seq pos!
#  {
#    my $newSequence = ($alignment->get_seq_by_pos($i))->revcom();
#
#  print "New sequence\n";
#
#    $newAlignment->add_seq($newSequence);
#  }
#
#  return $newAlignment; 
#}


# Helper function, accepts a reference to an array of primers, and
# returns a parallel-constructed array containing primers that are 
# reverse-complements of the input primers
sub buildReversePrimers
{
  my ($primers_r) = @_;

  my @reversePrimers = ();
  foreach my $primer(@{$primers_r})
  {
    my $minusStrandPrimer = $primer->clone();
    $minusStrandPrimer->reverseComplement();
    push(@reversePrimers, $minusStrandPrimer);
  }

  return @reversePrimers;
}

# Helper function, accepts a reference to an array of primers, and
# an analyzer, and returns an array of PrimerInfo objects, which are
# the analysis results
sub analyzeAll
{
  my ($primers_r, $analyzer) = @_;

  # TODO: probably want to sanity check these parameters
  my $measurements_r = [];
  foreach my $primer(@{$primers_r})
  {
    #print "_";
    my $primerInfo = $analyzer->analyze($primer);
    push(@{$measurements_r}, $primerInfo);
  }
  #print "\n";

  return $measurements_r;
}

# Helper function, build all the possible combinations of primers, and
# return the complete set of pairs
sub enumeratePairs
{
  my ($paramHash_r) = @_;

  my $forwardInfos_r = optionRequired($paramHash_r, "sorted_forward_infos");
  my $reverseInfos_r = optionRequired($paramHash_r, "sorted_reverse_infos");
  my $maxLength = optionRequired($paramHash_r, "max_length");
  my $minLength = optionWithDefault($paramHash_r, "min_length", 1);

  my $pcrPairs_r = [];
  my $forwardInfoCount = scalar(@{$forwardInfos_r});
  my $reverseInfoCount = scalar(@{$reverseInfos_r});
  my $previousFirstCompatibleIndex = 0; # Bound the lower end of the inner loop
  for(my $i = 0; $i < $forwardInfoCount; $i++)
  {
    my $forwardInfo = $forwardInfos_r->[$i];
    my $forwardStart = $forwardInfo->getLocation();
    my $forwardEnd = $forwardStart + $forwardInfo->getLength() - 1;

    # Used to bound the upper end of the inner loop search
    my $maxReverseLocation = $forwardStart + $maxLength - 1;
    
    # Used to help bound the lower end of the inner loop search
    my $previousCompatibleIndexFound = $FALSE;

    for(my $j = $previousFirstCompatibleIndex; $j < $reverseInfoCount; $j++)
    {
      # Careful - these reverse primer locations are expressed as their
      # positive strand starts and ends.  This means the reverse primer's
      # "start" location IS NOT the 5' start of the primer.
      my $reverseInfo = $reverseInfos_r->[$j];
      my $reverseEnd = $reverseInfo->getLocation();
      my $reverseStart = $reverseEnd - $reverseInfo->getLength() + 1;
      
      #print "    Trying pair $i ($forwardStart-$forwardEnd), $j ($reverseStart-$reverseEnd)\n";

      # Advance to the next compatible reverse primer by skipping all the
      # primers located too far 5' with respect to the forward primer
      if($previousCompatibleIndexFound == $FALSE)
      {
        if($reverseStart <= $forwardEnd)
        {
          next;
        }
        else
        {
          $previousFirstCompatibleIndex = $j;
          $previousCompatibleIndexFound = $TRUE;
        }
      }

      # Enforce lower length boundary
      if(($reverseEnd - $forwardStart + 1) < $minLength)
      {
        next;
      }

      # Enforce upper length boundary by using the calculated max location
      if($reverseEnd > $maxReverseLocation)
      {
        last;
      }
       
      #print "New pair adding to list of \"" .
      #  scalar(@oligoPairs) .
      #  "\" pairs.\n";
         
      my $newPair = LLNL::LAVA::PrimerSet::PCRPair->new(
        {
          "forward_info" => $forwardInfo,
          "reverse_info" => $reverseInfo,
        });
      push(@{$pcrPairs_r}, $newPair);
    }
  }

  return $pcrPairs_r;
}

sub buildMetricsArray
{
  my ($pairInfos_r) = @_;

  my $pairCount = scalar(@{$pairInfos_r});

  my $metricsArray_r = [];
  for(my $i = 0; $i < $pairCount; $i++)
  { 
    my $pairInfo = $pairInfos_r->[$i];

    my $penalty = $pairInfo->getPenalty();
    my $pair = $pairInfo->getAnalyzedPair();

    my $start = $pair->getStartLocation();
    my $end = $pair->getEndLocation();
    my $length = $pair->getLength();

    my $forwardInfo = $pair->getForwardInfo();
    my $reverseInfo = $pair->getReverseInfo();

    my $forwardLength = $forwardInfo->getLength();
    my $reverseLength = $reverseInfo->getLength();

    my $clearAt = $start + $forwardLength;
    my $clearThrough = $end - $reverseLength;

    # Guess we should embed the pair in the metric array too, so we
    # can handle it independently of the source array
    $metricsArray_r->[$i] = [$penalty, $start, $end, $length, 
      $forwardLength, $reverseLength, $clearAt, $clearThrough];
  }

  return $metricsArray_r;
}

sub reducePairInfosByPenalty
{
  my ($paramHash_r) = @_;

  my $pairInfos_r = optionRequired($paramHash_r, "pair_infos");
  my $maxPairs = optionRequired($paramHash_r, "max_pairs");

  # Sort the pairs by penalty
  my @sortedInfos =
    map { $_->[0] }
    sort { $a->[1] <=> $b->[1] }
    map {[$_, $_->getPenalty()] }
    @{$pairInfos_r}; 

  # Extract the X top pairs
  my @bestInfos = ();
  my $infoCount = scalar(@sortedInfos);
  for(my $i = 0; ($i < $infoCount) && ($i < $maxPairs); $i++)
  {
    $bestInfos[$i] = $sortedInfos[$i];
  }

  # Re-sort the pairs by start location
  @bestInfos =
    map { $_->[0] }
    sort { $a->[1] <=> $b->[1] }
    map {[$_, $_->getAnalyzedPair()->getStartLocation()]}
    @bestInfos;

  return \@bestInfos;
}

# Primers are identified by location:length, the combination of which is currently
# assumed to be unique, so be careful!
sub reducePrimersByOverlap
{
  my ($paramHash_r) = @_;

  my $maxOverlapPercent = optionRequired($paramHash_r, "max_overlap_percent");
  my $sortedByLocation_r = optionRequired($paramHash_r, "info_sorted_by_location");
  my $sortedByPenalty_r = optionRequired($paramHash_r, "info_sorted_by_penalty");

  # Quick dumb check for list length
  if(scalar(@{$sortedByLocation_r}) != scalar(@{$sortedByPenalty_r}))
  {
    confess("data error - lists don't have identical length! " .
      scalar(@{$sortedByLocation_r}) .
      " for location sorting, and " .
      scalar(@{$sortedByPenalty_r}) .
      " for penalty sorting.");
  }

  my $primerCount = scalar(@{$sortedByLocation_r});
  #print "Starting with $primerCount primers\n";
  print " $primerCount->";
   
  # Short-cut if we're at 100% overlap
  if($maxOverlapPercent == 100)
  {
    my @primerList = ();
    foreach my $primerInfo(@{$sortedByPenalty_r})
    {
      push(@primerList, $primerInfo);
    }

    # Go ahead and sort primers by location for their return :)
    @primerList = 
      map {$_->[0]}
      sort {$a->[1] <=> $b->[1]}
      map {[$_, $_->getLocation()] } 
      @primerList;
  
    print "" .
      scalar(@primerList) .
      ")\n";
    #print "Ending with " .
    #  scalar(@primerList) .
    #  " primers\n";
  
    return \@primerList;
  } # End shortcut for 100% overlap

  # The combination of location/length should be unique for each primer 
  my %unavailablePrimers = ();
  my @primerList = ();

  my @byPenaltyInfoLookup = ();
  for(my $infoIndex = 0; $infoIndex < $primerCount; $infoIndex++)
  {
    my $primerInfo = $sortedByPenalty_r->[$infoIndex];
    my $location = $primerInfo->getLocation();
    my $length = $primerInfo->getLength();
    my $primerTitle = "$location:$length"; 
    $byPenaltyInfoLookup[$infoIndex] = [$location, $length, $primerTitle];
  }
  my @byLocationInfoLookup = ();
  for(my $infoIndex = 0; $infoIndex < $primerCount; $infoIndex++)
  {
    my $primerInfo = $sortedByLocation_r->[$infoIndex];
    my $location = $primerInfo->getLocation();
    my $length = $primerInfo->getLength();
    my $primerTitle = "$location:$length"; 
    $byLocationInfoLookup[$infoIndex] = [$location, $length, $primerTitle];
  }

  # Priority to lowest penalty primers, knock out overlaps for each accepted primer. 
  for(my $outerIndex = 0; $outerIndex < $primerCount; $outerIndex++)
  {
    my $primerInfo = $sortedByPenalty_r->[$outerIndex];
    my ($location, $length, $primerTitle) = @{$byPenaltyInfoLookup[$outerIndex]};
    #my $location = $primerInfo->getLocation();
    #my $length = $primerInfo->getLength();
    #my $primerTitle = "$location:$length"; 

    #print "\n[outer $primerTitle]";
    # Skip this primer if its unavailable
    if(exists $unavailablePrimers{$primerTitle})
    {
      #print "X";
      next;
    }

    # Mark this new primer as unavailable now that it's chosen
    $unavailablePrimers{$primerTitle} = $TRUE;
    push(@primerList, $primerInfo);

    # Set upstream and downstream bounds for overlap checks
    # These will have to be adjusted if primer lengths are ever drasically different
    my $upstreamStart = $location - (2 * $length);
    if($upstreamStart < 0)
    {
      $upstreamStart = 0;
    }
    my $downstreamEnd = $location + $length - 1;

    # Iterate over primers within bounds and check overlap percent
    # Mark all primers with more than the cutoff percent of overlap as unavailable
    for(my $innerIndex = 0; $innerIndex < $primerCount; $innerIndex++)
    {
      my $innerInfo = $sortedByLocation_r->[$innerIndex];
      my ($innerLocation, $innerLength, $innerPrimerTitle) = @{$byLocationInfoLookup[$innerIndex]};

      #my $innerLocation = $innerInfo->getLocation();

      if($innerLocation < $upstreamStart)
      {
	#print "U";
	next;
      }

      if($innerLocation > $downstreamEnd)
      {
	#print "D";
	last;
      }

      #my $innerLength = $innerInfo->getLength();
      #my $innerPrimerTitle = "$innerLocation:$innerLength";
      #print "[inner $innerPrimerTitle]";
       
      # Assume inner is longer, swap if not true 
      my $shorterStart = $location;
      my $shorterLength = $length;
      my $shorterEnd = $shorterStart + $shorterLength - 1;
      
      my $longerStart = $innerLocation;
      my $longerLength = $innerLength;
      my $longerEnd = $longerStart + $longerLength - 1;
     
      if($innerLength > $length)
      {
        $shorterStart = $innerLocation;
	$shorterLength = $innerLength;
	$shorterEnd = $shorterStart + $shorterLength - 1;

	$longerStart = $location;
	$longerLength = $length;
	$longerEnd = $longerStart + $longerLength - 1;
      }
    
      my $earlierStart = $shorterStart;
      if($earlierStart > $longerStart)
      {
	$earlierStart = $longerStart;
      }
      my $laterEnd = $shorterEnd;
      if($laterEnd < $longerEnd)
      {
	$laterEnd = $longerEnd;
      }

      
      my $overlapCount = ($shorterLength + $longerLength) - ($laterEnd - $earlierStart + 1);
      # Negative value is distance between oligos, but not using that right now.
      if($overlapCount < 0)
      {
	$overlapCount = 0;
      }
      my $overlapPercent = ($overlapCount / $shorterLength) * 100;

      if($overlapPercent > $maxOverlapPercent)
      {
	#print "O";
        $unavailablePrimers{$innerPrimerTitle} = $TRUE;
      }
    } # End foreach primer by location
  } # End foreach primer by penalty

  # Go ahead and sort primers by location for their return :)
  @primerList = 
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->getLocation()] } 
    @primerList;

  print "" .
    scalar(@primerList) .
    ")\n";
  #print "Ending with " .
  #  scalar(@primerList) .
  #  " primers\n";

  return \@primerList;
}

# Final reduction step to keep from reporting heavily overlapping signatures
# Initially implemented to allow no overlap at all between signatures
sub reduceSignaturesByOverlap
{
  my ($paramHash_r) = @_;

  my $signatures_r = optionRequired($paramHash_r, "signatures");
  my $maxOverlapPercent = optionWithDefault($paramHash_r, "max_overlap_percent", 0);
  my $sortByScore = optionWithDefault($paramHash_r, "sort_by_score", $FALSE);
  my $sortByLocation = optionWithDefault($paramHash_r, "sort_by_location", $FALSE);
  if($sortByScore == $FALSE && $sortByLocation == $FALSE)
  {
    $sortByScore = $TRUE;
  }
  if($sortByScore == $TRUE && $sortByLocation == $TRUE)
  {
    $sortByLocation = $FALSE;
  }

  my $signatureCount= scalar(@{$signatures_r});
  #print "Starting with $primerCount primers\n";
  print " Reducing signatures $signatureCount->";
  
  # Short-cut if we're at 100% overlap
  if($maxOverlapPercent == 100)
  {
    my @signatureList = ();
    foreach my $signature(@{$signatures_r})
    {
      push(@signatureList, $signature);
    }

    # Go ahead and sort primers by location for their return :)
    @signatureList = 
      map {$_->[0]}
      sort {$a->[1] <=> $b->[1]}
      map {[$_, $_->getStartLocation()] } 
      @signatureList;
  
    print "" .
      scalar(@signatureList) .
      ")\n";
  
    return \@signatureList;
  } # End shortcut for 100% overlap

  my @byPenaltyLookup = ();
  my @byLocationLookup = ();
  for(my $sigIndex = 0; $sigIndex < $signatureCount; $sigIndex++)
  {
    my $signature = $signatures_r->[$sigIndex];

    my $penalty = $signature->getTag("lamp_penalty");
    my $startLocation = $signature->getStartLocation();
    my $length = $signature->getLength();

    $byPenaltyLookup[$sigIndex] = [$sigIndex, $penalty, $startLocation,
      $length];
    $byLocationLookup[$sigIndex] = [$sigIndex, $penalty, $startLocation,
      $length];
  }
  # Sort the data lookup by penalty
  # (doesn't need to be a map-sort-map now that the function calls
  # have been removed)
  @byPenaltyLookup = 
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->[1]]} 
    @byPenaltyLookup;

  # Sort the other data lookup by location
  @byLocationLookup = 
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->[2]]} # Here's the subtle difference.
    @byLocationLookup;

  # The combination of location/length should be unique for each primer 
  my %unavailableSignatures = (); # Hashed by original index
  my @signatureInfoList = (); # Accepted signature lookup data - not signatures!

  # Priority to lowest penalty primers, knock out overlaps for each accepted primer. 
  for(my $outerIndex = 0; $outerIndex < $signatureCount; $outerIndex++)
  {
    my ($originalIndex, $penalty, $location, $length) =
      @{$byPenaltyLookup[$outerIndex]};

    # Skip this signature if its unavailable
    if(exists $unavailableSignatures{$originalIndex})
    {
      #print "X";
      next;
    }

    # Mark this new primer as unavailable now that it's chosen
    $unavailableSignatures{$originalIndex} = $TRUE;
    push(@signatureInfoList, $byPenaltyLookup[$outerIndex]);

    # Set upstream and downstream bounds for overlap checks
    # These will have to be adjusted if signature lengths are ever drasically different
    my $upstreamStart = $location - (5 * $length);
    if($upstreamStart < 0)
    {
      $upstreamStart = 0;
    }
    my $downstreamEnd = $location + $length - 1;

    # Iterate over signatures within bounds and check overlap percent
    # Mark all primers with more than the cutoff percent of overlap as unavailable
    for(my $innerIndex = 0; $innerIndex < $signatureCount; $innerIndex++)
    {
      my ($innerOriginalIndex, $innerPenalty, $innerLocation, $innerLength) =
        @{$byLocationLookup[$innerIndex]};

      if($innerLocation < $upstreamStart)
      {
	#print "U";
	next;
      }

      if($innerLocation > $downstreamEnd)
      {
	#print "D";
	last;
      }

      # Assume inner is longer, swap if not true 
      my $shorterStart = $location;
      my $shorterLength = $length;
      my $shorterEnd = $shorterStart + $shorterLength - 1;
      
      my $longerStart = $innerLocation;
      my $longerLength = $innerLength;
      my $longerEnd = $longerStart + $longerLength - 1;
     
      if($innerLength > $length)
      {
        $shorterStart = $innerLocation;
	$shorterLength = $innerLength;
	$shorterEnd = $shorterStart + $shorterLength - 1;

	$longerStart = $location;
	$longerLength = $length;
	$longerEnd = $longerStart + $longerLength - 1;
      }
    
      my $earlierStart = $shorterStart;
      if($earlierStart > $longerStart)
      {
	$earlierStart = $longerStart;
      }
      my $laterEnd = $shorterEnd;
      if($laterEnd < $longerEnd)
      {
	$laterEnd = $longerEnd;
      }

      my $overlapCount = ($shorterLength + $longerLength) - ($laterEnd - $earlierStart + 1);
      # Negative value is distance between oligos, but not using that right now.
      if($overlapCount < 0)
      {
	$overlapCount = 0;
      }
      my $overlapPercent = ($overlapCount / $shorterLength) * 100;

      if($overlapPercent > $maxOverlapPercent)
      {
	#print "O";
        $unavailableSignatures{$innerOriginalIndex} = $TRUE;
      }
    } # End foreach primer by location
  } # End foreach primer by penalty

  # Go ahead and sort primers by for their return
  my $sortingIndex = 2; # Assume sorting by location, switch if not so
  if($sortByScore == $TRUE)
  {
    $sortingIndex = 1;
  }
  @signatureInfoList = 
    map {$_->[0]}
    sort {$a->[1] <=> $b->[1]}
    map {[$_, $_->[$sortingIndex]]} 
    @signatureInfoList;

  # Build the answer array from the original signature set indexes
  my @signatureList = ();
  foreach my $info_r(@signatureInfoList)
  {
    my ($originalIndex, $penalty, $location, $length) = @{$info_r};
    push(@signatureList, $signatures_r->[$originalIndex]);
  }

  print "" .
    scalar(@signatureList) .
    ")\n";
  #print "Ending with " .
  #  scalar(@primerList) .
  #  " primers\n";

  return \@signatureList;
}

sub flattenInfoData
{
  my ($paramHash_r) = @_;

  my $infoSet_r = optionRequired($paramHash_r, "info_set_ref");

  my $flattenedData_r = [];
  my $infoSetLength = scalar(@{$infoSet_r});
  for(my $infoIndex = 0; $infoIndex < $infoSetLength; $infoIndex++)
  {
    my $info = $infoSet_r->[$infoIndex];
    $flattenedData_r->[$infoIndex] = 
      [$info->getLocation(), $info->getLength(), $info->getPenalty()];
  }

  return $flattenedData_r;
}

sub generateDistancePenalties
{
  my ($maxDistance) = @_;

  # Going with a quadratic curve through the points (0,0), (35, 10), (60, 20)
  # Using y = ax^2 + bx + c form
  my ($aTerm, $bTerm, $cTerm) = solveCoefficients(0, 30, 10, 50, 20);
  #  print "aTerm $aTerm\n" .
  #    "bTerm $bTerm\n" .
  #    "cTerm $cTerm\n";

  my @penalties = ();
  for(my $i = 0; $i < $maxDistance; $i++)
  {
    $penalties[$i] = $aTerm * ($i ** 2) + $bTerm * $i + $cTerm;
  }

  return \@penalties;
}

# Cramer's rule in two easy steps to get the equations of the parabola satisfying
# the three given points, expanding the determinant the old fashioned way
sub solveCoefficients
{
  my ($zeroAt, $distanceA, $penaltyA, $distanceB, $penaltyB) = @_;

  # To simplify the terms
  my $a1 = $zeroAt ** 2;
  my $b1 = $zeroAt;
  my $c1 = 1;
  my $d1 = 0;
  my $a2 = $distanceA ** 2;
  my $b2 = $distanceA;
  my $c2 = 1;
  my $d2 = $penaltyA;

  my $a3 = $distanceB ** 2;
  my $b3 = $distanceB;
  my $c3 = 1;
  my $d3 = $penaltyB;

  # P = determinant of left-side terms
  # Px = determinant of left-side with respect to x
  # ...
  # So,
  # x = Px/P
  # y = Py/P
  # z = Pz/P

  # Solve Determinants
  my $p = findDeterminant($a1, $b1, $c1, $a2, $b2, $c2, $a3, $b3, $c3);

  my $px = findDeterminant($d1, $b1, $c1, $d2, $b2, $c2, $d3, $b3, $c3); 
  my $py = findDeterminant($a1, $d1, $c1, $a2, $d2, $c2, $a3, $d3, $c3); 
  my $pz = findDeterminant($a1, $b1, $d1, $a2, $b2, $d2, $a3, $b3, $d3); 

  my $solvedX = $px / $p;
  my $solvedY = $py / $p;
  my $solvedZ = $pz / $p;

  return ($solvedX, $solvedY, $solvedZ);
}

sub findDeterminant
{
  my ($a1, $b1, $c1, $a2, $b2, $c2, $a3, $b3, $c3) = @_;

  my $determinant = ($a1 * $b2 * $c3) + ($b1 * $c2 * $a3) + ($c1 * $a2 * $b3) - ($c1 * $b2 * $a3) - ($a1 * $c2 * $b3) - ($b1 * $a2 * $c3);

  #print "Det: $determinant\n";
  return $determinant;
}

