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

package LLNL::LAVA::PrimerSet;

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

LLNL::LAVA::PrimerSet - A pair of LLNL::LAVA::PrimerInfo objects

=head1 SYNOPSIS

  # Instantiation
  $pair = LLNL::LAVA::PrimerSet->new(
    {
      "start_location" => $startLocation,
      "end_location" => $endLocation,
      "length" => $length,
    });

  # Member functions are named with "get" prefixes as an indicator that
  # they can't be used to modify the PCR pair, just to get info about it.
  $location = $pair->getStartLocation();
  $location = $pair->getEndLocation();
  $length = $pair->getLength();
  
   
  # String in the form of "start-end".  Depending on how you created the pair,
  # the end location might be a smaller number than the start location.  It's
  # possible to have this behavior, we just don't explicitly support it yet.
  $rangeString = $pair->getRangeAsString();
      
  # Tag use (inherited from TagHolder)
  $oligo->setTag("useful_statistic", $usefulStatistic);
  $exists = $oligo->getTagExists("useful_statistic");
  $tm = $oligo->getTag("useful_statistic");

=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

PrimerSet is the object model of a group of oligos.  Sub-classes implement
behavior specific to different types of oligo combinations such as PCRPairs.

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $sigComponent = LLNL::LAVA::PrimerSet->new(
               {
                 "start_location" => $startLocation,
                 "end_location" => $endLocation,
                 "length" => $length, # Redundant, but keeping for now.
               });
 Function  : Creates a new LLNL::LAVA::PrimerSet 
 Arguments : Hash ref - Parameters for initialization including:
               start_location - Integer start location for the entire set
               end_location - Integer end location for the entire set
               length - Integer length/span of the entire set
 Example   : See Usage
 Returns   : A new LLNL::LAVA::PrimerSet

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

  # Start should be a smaller number than end at this point
  my $startLocation = optionRequired($paramHash_r, "start_location");
  my $endLocation = optionRequired($paramHash_r, "end_location");
  my $length = optionRequired($paramHash_r, "length");

  # Initialize via TagHolder
  my $this = $classType->SUPER::new({});

  # Start location is the up-strand 5' location
  $this->{"d_startLocation"} = $startLocation;
  $this->{"d_endLocation"} = $endLocation;
  $this->{"d_length"} = $length;

  return $this;
}

#-------------------------------------------------------------------------------

=head2 getStartLocation

 Usage     : $location = $pair->getStartLocation();
 Function  : Get the 5' start location of the forward primer
 Arguments : <n/a>
 Example   : See Usage
 Returns   : Integer - start location of the pair

=cut

sub getStartLocation
{
  my ($this) = @_;

  return $this->{"d_startLocation"};
}

#-------------------------------------------------------------------------------

=head2 getEndLocation

 Usage     : $location = $pair->getEndLocation();
 Function  : Get the 5' start location of the reverse primer, which is the
             3' end location of the PCR pair!
 Arguments : <n/a>
 Example   : See Usage
 Returns   : Integer - end location of the pair

=cut

sub getEndLocation
{
  my ($this) = @_;

  return $this->{"d_endLocation"};
}

#-------------------------------------------------------------------------------

=head2 getLength

 Usage     : $length = $pair->getLength();
 Function  : Get the lengh of this PCR pair, calculated from primer locations
 Arguments : <n/a>
 Example   : See Usage
 Returns   : Integer - length of the expected PCR amplicon

=cut

sub getLength
{
  my ($this) = @_;

  return $this->{"d_length"};
}

#-------------------------------------------------------------------------------

=head2 getRangeAsString

 Usage     : $rangeString = $pair->getRangeAsString();
 Function  : Get a string in the form of "start-end" (e.g. "123-456"), where
             the two numbers depend on primer locations
 Arguments : <n/a>
 Example   : See Usage
 Returns   : String - representation of the range of the primer pair

=cut

sub getRangeAsString
{
  my ($this) = @_;

  my $rangeString = $this->{"d_startLocation"} .
    "-" .
    $this->{"d_endLocation"};

  return $rangeString;
}

1; # Lame!

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
