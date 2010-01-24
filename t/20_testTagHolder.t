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

use Carp;

# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'

######################### We start with some black magic to print on failure.

# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)
BEGIN {print "1..7\n";}

my $TRUE = 1;
my $FALSE = 0;

my $loaded = 0;

use LLNL::LAVA::TagHolder ":standard";

END { if($loaded == 0)
      {
        print "NOT ok 1\n";
      }
    }
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code)

# Turn off buffering of output
$| = 1;

{ # fake main()
  my $testSuccess = $TRUE; # Status of the test


  my $testTagName = "foo_test";
  my $testTagValue = "bar_value";

  my $secondTagName = "another_test";
  my $secondTagValue = "another_value";
   
  # Test taht instantiation works 
  my $tagHolder;
  eval
  {
    $tagHolder = LLNL::LAVA::TagHolder->new();
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Failed to instantiate a new TagHolder object";
  }
  printTestResult(2, $testSuccess);

  # Test that a tag value doesn't exist before setting any.
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $impossibleValue = $tagHolder->getTag($testTagName);  
  };
  if(! $@)
  {
    $testSuccess = $FALSE;
    print "Tag retrieval should have failed for uninitialized tag name.";
  }
  printTestResult(3, $testSuccess);

  # Test that a tag value doesn't appear to exist yet
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $tagExists = $tagHolder->getTagExists($testTagName);
    if($tagExists == $TRUE)
    {
      confess("getTagExists() returned TRUE for a non-existent tag");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error while testing tag existence: $@";
  }
  printTestResult(4, $testSuccess);

  # Test setting and getting a tag
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    $tagHolder->setTag($testTagName, $testTagValue);

    my $retrievedValue = $tagHolder->getTag($testTagName);
    if($retrievedValue ne $testTagValue)
    {
      confess("getTag() returned a different value than was set in setTag()");
    }

    $tagHolder->setTag($secondTagName, $secondTagValue);
    $retrievedValue = $tagHolder->getTag($testTagName);
    my $secondRetrievedValue = $tagHolder->getTag($secondTagName);
    if($retrievedValue ne $testTagValue ||
        $secondRetrievedValue ne $secondTagValue)
    {
      confess("when multiple tags exist, getTag() returned a different value " .
        "than was set in setTag()");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error retrieving a tag that was set: $@";
  }
  printTestResult(5, $testSuccess);

  # Test that an existing tag appears to exist with getTagExists()
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $tagExists = $tagHolder->getTagExists($testTagName); 
    if($tagExists == $FALSE)
    {
      confess("getTagExists() returned FALSE for a tag that exists");
    }
    $tagExists = $tagHolder->getTagExists($secondTagName); 
    if($tagExists == $FALSE)
    {
      confess("getTagExists() returned FALSE for a tag that exists");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when testing tag existence option: $@";
  }
  printTestResult(6, $testSuccess);

  # Test that tag names return properly trough getAllTagNames()
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my @allTagNames = $tagHolder->getAllTagNames();
    if(scalar(@allTagNames) == 0)
    {
      confess("retrieval of a single tag with getAllTagNames() failed");
    }
    if(scalar(@allTagNames) != 2)
    {
      confess("getAllTagNames() failed to return the exact set of tags " .
        "that should exist");
    }

    # Get 2 tag names back, and test that they match the 2 we set
    my $firstWasFirst = $FALSE;
    my $foundFirst = $FALSE;
    my $foundSecond = $FALSE;
    if($allTagNames[0] eq $testTagName)
    {
      $firstWasFirst = $TRUE;
      $foundFirst = $TRUE;
    }
    elsif($allTagNames[1] eq $testTagName)
    {
      $foundFirst = $TRUE;
    }

    if($firstWasFirst == $TRUE &&
       $allTagNames[1] eq $secondTagName)
    {
      $foundSecond = $TRUE;
    }
    elsif($firstWasFirst == $FALSE &&
       $allTagNames[0] eq $secondTagName)
    {
      $foundSecond = $TRUE;
    }

    if($foundFirst == $FALSE ||
       $foundSecond == $FALSE)
    {
      confess("getAllTagNames() returned scrambled values when 2 tags existed");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when testing getAllTagNames(): $@";
  }
  printTestResult(7, $testSuccess);

  print "END\n";
} # end fake main()

sub printTestResult
{
  my ($testCount, $testSuccess) = @_;
  
  if($testSuccess == $TRUE)
  {
    print "ok $testCount\n";
  }
  else
  {
    print "NOT ok $testCount\n";
  }
}

