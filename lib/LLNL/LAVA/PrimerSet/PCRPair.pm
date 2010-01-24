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

package LLNL::LAVA::PrimerSet::PCRPair;

use strict;
use warnings;
use Carp;

use vars qw(@ISA);

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::PrimerInfo; # Recieves as parameter and holds

use LLNL::LAVA::PrimerSet; # is-a
@ISA = ("LLNL::LAVA::PrimerSet");

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::PrimerSet::PCRPair - A pair of LLNL::LAVA::PrimerInfo objects

=head1 SYNOPSIS

  # Instantiation
  $pair = LLNL::LAVA::PrimerSet::PCRPair->new(
    {
      # PrimerInfo for PCR pairs need to be from opposite strands
      "forward_info" => $forwardInfo,
      "reverse_info" => $reverseInfo,
    });

  # Member functions are named with "get" prefixes as an indicator that
  # they can't be used to modify the PCR pair, just to get info about it.
  $location = $pair->getStartLocation();
  $location = $pair->getEndLocation();
  $length = $pair->getLength();
  
  $forwardInfo = $pair->getForwardInfo();
  $reverseInfo = $pair->getReverseInfo();
 
  $forwardStrand = $pair->getForwardStrand();
  $reverseStrand = $pair->getReverseStrand();
   
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

PrimerSet::PCRPair is the object model of a pair of oligos from opposite strands that
can be used as PCR primers.

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $sigComponent = LLNL::LAVA::PrimerSet::PCRPair->new(
               {
                 # Forward primer needs to be opposite strand of the reverse
                 # primer.
                 "forward_info" => $forwardInfo,
                 "reverse_info" => $reverseInfo,
               });
 Function  : Creates a new LLNL::LAVA::PrimerSet::PCRPair 
 Arguments : Hash ref - Parameters for initialization including:
               forward_info - LLNL::LAVA::PrimerInfo as a forward primer 
               reverse_info - LLNL::LAVA::PrimerInfo as a reverse primer
 Example   : See Usage
 Returns   : A new LLNL::LAVA::PrimerSet::PCRPair

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
  my $forwardInfo = optionRequired($paramHash_r, "forward_info");
  my $reverseInfo = optionRequired($paramHash_r, "reverse_info");
  # Might want to check that each oligo conforms to an interface, but it feels
  # like that is getting in the way right now and causing strict typing where
  # it isn't needed.

  my $forwardStrand = $forwardInfo->getTag("strand");
  my $reverseStrand = $reverseInfo->getTag("strand");
  # Kinda a strange check, but we're not trying to keep users from reversing
  # the primers, just to make sure the strands aren't identical.
  # Also, phrasing it this way makes it possible for primers to have other
  # strand designations without completely destroying things. 
  if($forwardStrand eq "plus" && $reverseStrand eq "plus")
  {
    confess("data error - forward and reverse primers were both \"plus\" strand, " .
      "and one of them needs to be different");
  }
  if($forwardStrand eq "minus" && $reverseStrand eq "minus")
  {
    confess("data error - forward and reverse primers were both \"minus\" strand, " .
      "and one of them needs to be different");
  }

  # With plus strand forward primers, reverse primers are minus stranded, so their
  # location must be greater than the forward primer's location 
  # Vice-versa for minus strand forward primers

  # Initialize as if we have a plus strand forward primer
  my $pairStart = $forwardInfo->getLocation();
  my $pairEnd = $reverseInfo->getLocation();

  # Again, being very specific with the strand designations because we don't
  # want to have a chance at working if the primer strands aren't this simple
  if($forwardStrand eq "minus" &&
     $reverseStrand eq "plus")
  {
    $pairStart = $reverseInfo->getLocation();
    $pairEnd = $forwardInfo->getLocation();
  }

  # Start should be a smaller number than end at this point
  my $length = $pairEnd - $pairStart + 1;

  # Initialize via PrimerSet (and TagHolder)
  my $this = $classType->SUPER::new(
    {
      "start_location" => $pairStart,
      "end_location" => $pairEnd,
      "length" => $length,
    });

  # Handle component properties specialized for PCRPair PrimerSets
  $this->{"d_forwardInfo"} = $forwardInfo;
  $this->{"d_forwardStrand"} = $forwardStrand;
  $this->{"d_reverseInfo"} = $reverseInfo;
  $this->{"d_reverseStrand"} = $reverseStrand;

  return $this;
}

#-------------------------------------------------------------------------------

=head2 getForwardInfo

 Usage     : $oligo = $pair->getForwardInfo();
 Function  : Get a reference to this pair's forward primer info
 Arguments : <n/a>
 Example   : See Usage
 Returns   : LLNL::LAVA::PrimerInfo - the forward primer info

=cut

sub getForwardInfo
{
  my ($this) = @_;

  return $this->{"d_forwardInfo"};
}

#-------------------------------------------------------------------------------

=head2 getForwardStrand

 Usage     : $oligo = $pair->getForwardStrand();
 Function  : Get the strand that the forward primer is on
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - "plus" or "minus" (or perhaps something else eventually)

=cut

sub getForwardStrand
{
  my ($this) = @_;

  return $this->{"d_forwardStrand"};
}

#-------------------------------------------------------------------------------

=head2 getReverseInfo

 Usage     : $oligo = $pair->getReverseInfo();
 Function  : Get a reference to this pair's reverse primer info
 Arguments : <n/a>
 Example   : See Usage
 Returns   : LLNL::LAVA::PrimerInfo - the reverse primer info

=cut

sub getReverseInfo
{
  my ($this) = @_;

  return $this->{"d_reverseInfo"};
}

#-------------------------------------------------------------------------------

=head2 getReverseStrand

 Usage     : $oligo = $pair->getReverseStrand();
 Function  : Get the strand that the reverse primer is on
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - "plus" or "minus" (or perhaps something else eventually)

=cut

sub getReverseStrand
{
  my ($this) = @_;

  return $this->{"d_reverseStrand"};
}

1; # Lame!

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
