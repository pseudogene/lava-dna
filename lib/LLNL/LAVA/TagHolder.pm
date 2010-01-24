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

package LLNL::LAVA::TagHolder;

use strict;
use Carp;

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::TagHolder - Generic name/value pair functionality with get/set 

=head1 SYNOPSIS

  # This is a trimmed down version of a tag library that was designed to
  # inter-operate with tag tables in a relational database.

  # Since we don't have to interact with a DB, the primary benefit of
  # using the TagHolder library is that it can cause a fatal run-time error
  # when a tag lookup fails. 

  use LLNL::LAVA::TagHolder;

  # Empty initialization
  my $tagHolder = LLNL::LAVA::TagHolder->new();

  # Setting a tag value
  $tagHolder->setTag("foo", $fooValue);
  # Setting a tag value for non-scalars doesn't work, so use references
  $tagHolder->setTag("foo", \@fooValues);

  # Get a tag's value 
  my $value = getTag("foo");

  # Returns $TRUE or $FALSE
  my $tagExists = $tagHolder->getTagExists("foo");

  # Get all tag names
  my @allTagNames = $tagHolder->getAllTagNames();

=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

TagHolder is a container for the tags of an object, held in
name/value pairs 

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : my $tagHolder = LLNL::LAVA::TagHolder->new();
 Function  : Creates a new LLNL::LAVA::TagHolder object 
 Arguments : <n/a>
 Example   : See Usage
 Returns   : A new LLNL::LAVA::TagHolder

=cut

sub new
{
  my ($classType, $paramHashRef) = @_;
  # Hmm, used to accept defaults as parameters, not that this is a simplified
  # version of the module, not sure we need to.

  my $this = {};
  bless($this, $classType);

  $this->{"d_th_tagsRef"} = {}; # Actual tag info
  
  return $this;
} # End new

#-------------------------------------------------------------------------------

=head2 getTagExists

 Usage     : $existenceBool = $tagHolder->getTagExists("foo");
 Function  : Probes the tagHolder for whether the requested name/value
             property exists 
 Arguments : String - name of the property to test for
 Example   : See Usage
 Returns   : Bool ($TRUE or $FALSE) - whether the requested property exists

=cut

sub getTagExists
{
  my ($this, $tagName) = @_;

  if(!defined $tagName)
  {
    confess("programming error - getTagExists() requires a tag name " .
      "as first parameter");
  }

  if(exists $this->{"d_th_tagsRef"}->{$tagName})
  {
    return $TRUE;
  }
  else
  {
    return $FALSE;
  }
}

#-------------------------------------------------------------------------------

=head2 getTag

 Usage     : $fooValue = $tagHolder->getTag("foo"); 
 Function  : Get the value of the requested property
 Arguments : String - name of the property to find the value of
 Example   : See Usage
 Returns   : Scalar - value of the property

=cut

sub getTag
{
  my ($this, $tagName) = @_;

  if(!defined $tagName)
  {
    confess("programming error - getTag() requires a tag name " .
      "as first parameter");
  }

  # Gonna be mean and cause error if it doesn't exist
  if(!exists $this->{"d_th_tagsRef"}->{$tagName})
  {
    confess("data error - did not retrieve tag \"$tagName\" since " .
      "it didn't exist");
  }

  return $this->{"d_th_tagsRef"}->{$tagName};
} # End getTag

#-------------------------------------------------------------------------------

=head2 setTag

 Usage     : $storage->setTag($tagName, $tagValue); 
             $storage->setTag($tagName, \@tagValues); 
 Function  : Adds the tag name/value pair to this result
 Arguments : String - name of the property to set
             Scalar - value of the property 
 Example   : See Usage
 Returns   : Scalar - the new tag value

=cut

sub setTag
{
  my ($this, $tagName, $tagValue) = @_;

  if(!defined $tagName)
  {
    confess("programming error - setTag() requires a tag name " .
      "as first parameter");
  }
  if(!defined $tagValue)
  {
    confess("programming error - setTag() requires a tag value " .
      "as second parameter");
  }

  $this->{"d_th_tagsRef"}->{$tagName} = $tagValue;

  return $this->{"d_th_tagsRef"}->{$tagName};
} # End setTag

#-------------------------------------------------------------------------------

=head2 getAllTagNames

 Usage     : @tagNames = $tagHolder->getAllTagNames(); 
 Function  : Get an array of all the tag names 
 Arguments : <n/a>
 Example   : See Usage
 Returns   : Array of String - all the tag names for this storage

=cut

sub getAllTagNames
{
  my ($this) = @_;

  my @tagNames = keys(%{$this->{"d_th_tagsRef"}});

  return @tagNames;
}

1; # Got one?

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
