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

package LLNL::LAVA::Oligo;

use strict;
use warnings;
use Carp;

use vars qw(@ISA);

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::TagHolder; # is-a
@ISA = ("LLNL::LAVA::TagHolder");

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::Oligo - Individual oligo sequence representation

=head1 SYNOPSIS

  use LLNL::LAVA::Oligo;

  # Instantiation
  $oligo = LLNL::LAVA::Oligo->new(
    {
      "sequence" => $sequence,
      "strand" => $strand, # becomes a tag, not a member, so retrieve the 
                           # value with $oligo->getTag("strand");
      "location" => $location, # 5' location, should be 0-indexed
    });

  # Can also instantiate with results from Bio::Tools::Run::Primer3
  $oligo = LLNL::LAVA::Oligo->newFromPrimer3($result);

  # Get the current sequence
  $sequence = $oligo->sequence();
  # Set to new sequence
  $oligo->sequence($newSequence);

  # Get the current location
  $location = $oligo->location();
  # Set to new location
  $oligo->location($newLocation);

  # Get the length
  $length = $oligo->length();

  # This special mutator does NOT return a reverse complement Oligo, it 
  # transforms this oligo into the reverse complement of itself.
  $oligo->reverseComplement(); # ACTUALLY CHANGES THE OLIGO!

  # Awkwardly duplicate - should probably move to using Storable's dclone()
  $newOligo = $oligo->clone();

  # Tag use
  $oligo->setTag("primer3_tm", $primer3TM);
  $exists = $oligo->getTagExists("primer3_tm");
  $tm = $oligo->getTag("primer3_tm");

=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

Oligo is the object model of a DNA oligo that is used as a 
component in a signature.

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $oligo = LLNL::LAVA::Oligo->new(
               {
                 "sequence" => $sequence,                # required
                 "strand" => $strand, #"plus" or "minus" # optional 
                     # default strand is "plus", but strand becomes a tag, not
                     # a member, so retrieve with $oligo->getTag("strand");
                 "location" => $location,                # optional (0 = default)
               });
 Function  : Creates a new LLNL::LAVA::Oligo 
 Arguments : Hash ref - Parameters for initialization include:
               sequence - sequence for the component 
               strand - "plus" or "minus" (stored as a tag!)
               location - 5-prime sequence index of this component
 Example   : See Usage
 Returns   : A new LLNL::LAVA::Oligo

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
  my $sequence = optionRequired($paramHash_r, "sequence");
  my $location = optionWithDefault($paramHash_r, "location", 0);
  # Strand becomes a tag, not a member
  my $strand = optionWithDefault($paramHash_r, "strand", "plus"); 

  # Initialize via TagHolder, not expecting default information 
  my $this = $classType->SUPER::new({});

  # Handle component properties
  $this->{"d_sequence"} = $sequence;
  $this->{"d_location"} = $location;

  # Strand of an oligo isn't always straight-forward, because oligos can
  # be built based on both strands (like the FIP and BIP of a LAMP signature).
  # Because it's not binary and required, we're tracking strand with a tag.
  $this->setTag("strand", $strand);

  return $this;
}

#-------------------------------------------------------------------------------

=head2 newFromPrimer3

 Usage     : $oligo = LLNL::LAVA::Oligo->newFromPrimer3($primer3Answer);
 Function  : Factory function for creating a new LLNL::LAVA::Oligo 
 Arguments : Hash ref - oligo returned from Bio::Tools::Run::Primer3
 Example   : See Usage
 Returns   : A new LLNL::LAVA::Oligo

=cut

sub newFromPrimer3
{
  my ($classType, $result) = @_;

  if(!defined $result)
  {
    confess("programming error - first parameter is a required result hash " .
    "ref as output from Bio::Tools::Run::Primer3");
  }
  if(ref($result) ne "HASH")
  {
    confess("programming error - first parameter must be a result hash " .
    "ref as output from Bio::Tools::Run::Primer3");
  } 

  #foreach my $key(keys(%{$result}))
  #{
  #  print "  $key->" .
  #    $result->{$key} .
  #    "\n";
  #}

  my $locationAndLength = $result->{"PRIMER_INTERNAL_OLIGO"};
  my ($location, $length) = split(",", $locationAndLength);

  my $oligo = $classType->new(
      {
      "sequence" => $result->{"PRIMER_INTERNAL_OLIGO_SEQUENCE"},
      "location" => $location,
      "strand" => "plus",
      });

  # Everything else is a tag
  $oligo->setTag("gc_percent",
      $result->{"PRIMER_INTERNAL_OLIGO_GC_PERCENT"});
  $oligo->setTag("primer3_tm",
      $result->{"PRIMER_INTERNAL_OLIGO_TM"});
  $oligo->setTag("primer3_penalty",
      $result->{"PRIMER_INTERNAL_OLIGO_PENALTY"});
  # Not actually using these right now
  #$oligo->setTag("self_end_score", 
  #    $result->{"PRIMER_INTERNAL_OLIGO_SELF_END"});
  #$oligo->setTag("self_any_score",
  #    $result->{"PRIMER_INTERNAL_OLIGO_SELF_ANY"});

  # Commented out because redundant with the instantiation above
  ## Oligo results are always on the plus strand 
  #$oligo->setTag("strand", "plus");

  return $oligo;
}

#-------------------------------------------------------------------------------

=head2 sequence

 Usage     : $sequence = $oligo->sequence();
             $oligo->sequence($newSequence);
 Function  : Get/Set pair for sequence
 Arguments : (optional) String - new sequence for the oligo
 Example   : See Usage
 Returns   : String - sequence of the oligo

=cut

sub sequence
{
  my ($this, $newSequence) = @_;

  if(defined $newSequence)
  {
    $this->{"d_sequence"} = $newSequence;
  }

  return $this->{"d_sequence"};
}

#-------------------------------------------------------------------------------

=head2 location

 Usage     : $location = $oligo->location();
             $oligo->location($newLocation);
 Function  : Get/Set pair for location
 Arguments : (optional) Integer - new location for the oligo
 Example   : See Usage
 Returns   : Integer - location of the oligo;

=cut

sub location
{
  my ($this, $newLocation) = @_;

  if(defined $newLocation)
  {
    $this->{"d_location"} = $newLocation;
  }

  return $this->{"d_location"};
}

#-------------------------------------------------------------------------------

=head2 length

 Usage     : $length = $oligo->length();
 Function  : Get the length of this oligo, calculated from the sequence content
 Arguments : <n/a>
 Example   : See Usage
 Returns   : Integer - length of the oligo

=cut

sub length
{
  my ($this) = @_;

  # A hair more expensive to calculate than memoizing the sequence length,
  # but it cleans up the code.
  return length($this->{"d_sequence"});
}

#-------------------------------------------------------------------------------

=head2 reverseComplement

 Usage     : $oligo->reverseComplement();
 Function  : This special mutator does NOT return a reverse complement of this
             Oligo, it transforms this oligo into the reverse complement
             of itself. 
 Arguments : <n/a>
 Example   : See Usage
 Returns   : <n/a>

=cut

sub reverseComplement
{
  my ($this) = @_;

  my $sequence = $this->sequence();
  my $origLocation = $this->location();
  my $strand = $this->getTag("strand");
  if($strand ne "plus" &&
     $strand ne "minus")
  {
    confess("data error - unrecognized strand type \"$strand\"");
  }

  if($sequence =~ /[^ACGTRYMKSWHBVDNX]/i)
  {
    confess("data error - only DNA vocabulary (acgtrymkswhbvdnxACGTRYMKSWHBVDNX) " .
      "is valid for an oligo - unrecognized base in the oligo \"$sequence\"");
  }

  # Invert the sequence alphabet, and reverse the sequence in 2 steps, then
  # set this object's "sequence" to the newly calculated reverse complement
  # "i" and "I" are left untouched
  $sequence =~ 
    tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  my $revComp = reverse($sequence);
  $this->sequence($revComp);

  # IF this oligo's strand is simply "plus" or "minus", invert the location, 
  # which tracks the outer 5' base, and invert the strand.
  # When we add more complicated strand descriptors, this needs to be re-visited. 
  if($strand eq "plus")
  {
    $this->location($origLocation + $this->length() - 1);
    $this->setTag("strand", "minus");
  }
  elsif($strand eq "minus")
  {
    $this->location($origLocation - $this->length() + 1);
    $this->setTag("strand", "plus");
  }

  # What's a logical return value for a mutator like this?
  return $TRUE;
}

#-------------------------------------------------------------------------------

=head2 clone

 Usage     : $newOligo = $oligo->clone();
 Function  : Duplicates an oligo, but DOES NOT perform a deep copy of tag values
 Arguments : <n/a>
 Example   : See Usage
 Returns   : A new LLNL::LAVA::Oligo

=cut


sub clone
{
  my ($this) = @_;

  my $newOligo = LLNL::LAVA::Oligo->new(
    {
      "sequence" => $this->sequence(),
      "location" => $this->location()
    });

  # Ideally we'd lean on the TagHolder to help us with cloning
  # We could also just use Storable's dclone()?
  my @tagNames = $this->getAllTagNames();
  foreach my $tagName(@tagNames)
  {
    $newOligo->setTag($tagName, $this->getTag($tagName));
  }

  return $newOligo;
}

1; # Lame!

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
