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

use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::LocatableSeq;

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

$| = 0;

{ # Fake main() to enforce scope
  my %options;
  my %optionMap =
    (
      "slava_alignment_fasta=s" => \$options{"slava_alignment_fasta"},
      "slava_segment_max_length=i" => \$options{"slava_segment_max_length"},
      "slava_segment_overlap_length=i" => \$options{"slava_segment_overlap_length"},
      "slava_output_file=s" => \$options{"slava_output_file"},
      "slava_temp_file_path=s" => \$options{"slava_temp_file_path"},
      "slava_path_to_lava=s" => \$options{"slava_path_to_lava"},
      "option_file|options_file=s" => \$options{"option_file"},
    );

  # Take the option file as parameter, so I don't have to process 
  # parameters here yet.

  # Given an input alignment and a length, prepare however many LAVA runs are
  # needed (try middle-hard, so need flag for disabling deep parts of the
  # solution plan)
 
  my %optionDefaults =
    (
      "slava_output_file" => "./lamp_signatures.fasta",
      "slava_temp_file_path" => "./",
      "slava_path_to_lava" => "./lava.pl",

      "slava_segment_max_length" => 5000,
      "slava_segment_overlap_length" => 400,
    );

  my $usageString = "Usage:\n" .
    "./slava.pl  \n" .
      "    --slava_alignment_fasta <fasta_file>\n" .
      "    --option_file <options_xml>\n" .
      "    [--slava_output_file <filename, default=" .
        $optionDefaults{"slava_output_file"} . ">]" .
      "    [--slava_temp_file_path <path, default=" .
        $optionDefaults{"slava_temp_file_path"} . ">]" .
      "    [--slava_path_to_lava <path, default=" .
        $optionDefaults{"slava_path_to_lava"} . ">]" .
      "    [--slava_segment_max_length <length, default=" .
        $optionDefaults{"slava_segment_max_length"} . ">]" .
      "    [--slava_segment_overlap_length <length, default=" .
        $optionDefaults{"slava_segment_overlap_length"} .
	">]\n";

  GetOptions(%optionMap);
  loadOptionsFromFile(\%options);
  my $options_r = \%options;


  my $alignmentFastaName = optionRequired($options_r, "slava_alignment_fasta", 
    $usageString);

  my $optionFile = optionRequired($options_r, "option_file", $usageString);

  my $outputFileName = optionWithDefault($options_r, "slava_output_file", 
    $optionDefaults{"slava_output_file"});
  my $tempFilePath = optionWithDefault($options_r, "slava_temp_file_path",
    $optionDefaults{"slava_temp_file_path"});
  if(-e $tempFilePath)
  {
    if(! -d $tempFilePath)
    {
      confess("Cannot create directory \"$tempFilePath\" for slava_temp_file_path " .
        "because a file already exists with that name");
    }
  }
  else
  {
    mkdir($tempFilePath);
  }
  my $pathToLava = optionWithDefault($options_r, "slava_path_to_lava",
    $optionDefaults{"slava_path_to_lava"});

  my $segmentMaxLength = optionWithDefault($options_r, "slava_segment_max_length",
    $optionDefaults{"slava_segment_max_length"});
  my $segmentOverlapLength = optionWithDefault($options_r, 
    "slava_segment_overlap_length", 
    $optionDefaults{"slava_segment_overlap_length"});

  if($segmentMaxLength <= $segmentOverlapLength)
  {
    confess("Parameter error - \"slava_segment_max_length\" must be larger " .
      "than \"slava_segment_overlap_length\"");
  }

  # Load the input alignment, could be a single sequence
  my $alignIN = Bio::AlignIO->new(-file => "< $alignmentFastaName");
  my $inputMSA = $alignIN->next_aln();
  my $inputLength = $inputMSA->length();

  # Slice the segments up.
    # create a log file of all the individual segments, to help with checkpointing

  my @alignmentSegments = ();
  my $offset = 0;
  while($offset < $inputLength)
  {
    # If this segment is smaller than the overlap length, then we're at
    # the end of the alignment, so stop, and don't add the under-length segment

    # CAUTION - slice() requires 1-indexed positions (inclusive), so 
    # coordinates are shifted by 1 for the actual function call (but not 
    # before because it felt awkward)
    my $startPosition = $offset;
    my $endPosition = $startPosition + $segmentMaxLength - 1;
    if($endPosition >= $inputLength)
    {
      $endPosition = $inputLength - 1;
    }
    
    # Don't accept a final tiny fragment
    my $sliceLength = $endPosition - $startPosition + 1;
    if($sliceLength < $segmentOverlapLength)
    {
      last;
    }

    my $subAlignment = $inputMSA->slice($startPosition + 1, $endPosition + 1);

    push(@alignmentSegments, $subAlignment);

    # Housekeeping
    $offset = $offset + $segmentMaxLength - $segmentOverlapLength;
  }

  my $numberOfSegments = scalar(@alignmentSegments); 
  print "Split into \"$numberOfSegments\" segments\n" .
    "  $inputLength bp alignment\n" .
    "  $segmentMaxLength bp per segment\n" .
    "  $segmentOverlapLength bp overlap\n"; 
    
  # Prepare commands to execute each segment
  # (log file tracks command status for each segment?)
  # Write each sub-alignment out to disk, and prepare the command line for each
  my @commandSet = ();
  for(my $i = 0; $i < $numberOfSegments; $i++)
  {
    my $subAlignment = $alignmentSegments[$i];

    # Pad the identifier with zeroes
    my $printableNumber = sprintf("%05d", $i); 
    my $subAlignmentFile = $tempFilePath . "/sub_alignment_$printableNumber.fasta";

    my $subAlignOut = Bio::AlignIO->new(-file => "> $subAlignmentFile", -format => "fasta");
    $subAlignOut->write_aln($subAlignment);

    my $subAlignmentResultsFile = $tempFilePath . "/segment_results_$printableNumber.lamp";
    my $command = "$pathToLava --alignment_fasta $subAlignmentFile" .
      " --output_file $subAlignmentResultsFile" .
      " --option_file $optionFile";
    push(@commandSet, $command);
    # TODO: Keep a separate "output file set" list, which we check for existence and load results from...
  }
  
  # Execute all prepared commands
  foreach my $command(@commandSet)
  {
    print "Executing: $command\n";
    my $commandOutput = `$command`;
    print "$commandOutput\n";
  }

  # Combine results from all the individual segments
    # filter down to max overlap per signature
  # Make this an external script, so users can run it on their own if
  # they checkpointed back in
      
  # Write out answers, really need some more input and output options

  
} # End fake main()


