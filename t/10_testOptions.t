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
BEGIN {print "1..10\n";}

my $TRUE = 1;
my $FALSE = 0;

my $loaded = 0;

use LLNL::LAVA::Options ":standard";

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

  # Set up option hash
  my %options;
  my $testBazDefault = "foobie";
  $options{"baz"} = $testBazDefault;
  my $listValue01 = "1_$testBazDefault";
  my $listValue02 = "2_$testBazDefault";
  $options{"baz_list"} = [$listValue01, $listValue02];

  # Test that optionRequired() causes an error when called with non-existent
  # values - the error here is intentional
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $nonexistentValue = optionRequired(\%options, "foo");
  };
  if(! $@)
  {
    $testSuccess = $FALSE;
    print "Should have attempted to exit when option \"foo\" didn't exist\n";
  }
  printTestResult(2, $testSuccess);

  # Test that optionRequired does not cause error when option exists
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $bazValue = optionRequired(\%options, "baz");
    if($bazValue ne $testBazDefault)
    {
      confess("optionRequired appears to corrupt the value!");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when retrieving option: $@";
  }
  printTestResult(3, $testSuccess);

  # Test that optionRequired properly returns list-style values
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my @listValue = optionRequired(\%options, "baz_list");
    if(scalar(@listValue) != 2 ||
       $listValue[0] ne $listValue01 ||
       $listValue[1] ne $listValue02)
    {
      confess("optionRequired appears to corrupt the list value!");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when retrieving list-style option: $@";
  }
  printTestResult(4, $testSuccess);

  # Test that optionWithDefault returns existing values  
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $bazValue = optionWithDefault(\%options, "baz", "bum");
    if($bazValue ne $testBazDefault)
    {
      confess("optionWithDefault appears to corrupt the value!");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when retrieving option: $@";
  }
  printTestResult(5, $testSuccess);

  # Test that optionWithDefault returns default values  
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $localTestDefault = "buzz";
    my $bazValue = optionWithDefault(\%options, "bee", $localTestDefault);
    if($bazValue ne $localTestDefault)
    {
      confess("optionWithDefault appears to corrupt the default!"); 
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when retrieving option: $@";
  }
  printTestResult(6, $testSuccess);

  # Test that optionWithDefault properly returns individual list-style values
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $localFirstDefault = "01_buzz";
    my @localListDefault = ($localFirstDefault);

    my @bazValue = optionWithDefault(\%options, "bee_list", @localListDefault);
    if(scalar(@bazValue) != 1 ||
       $bazValue[0] ne $localFirstDefault)
    {
      confess("optionWithDefault appears to corrupt the single-item-list value");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when retrieving list-style option: $@";
  }
  printTestResult(7, $testSuccess);

  # Test that optionWithDefault properly returns multiple list-style values
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $localFirstDefault = "01_buzz";
    my $localSecondDefault = "02_buzz";
    my @multipleDefaultList = ($localFirstDefault, $localSecondDefault);

    my @listValue = 
      optionWithDefault(\%options, "bee_list", @multipleDefaultList);
    if(scalar(@listValue) != 2 ||
       $listValue[0] ne $localFirstDefault ||
       $listValue[1] ne $localSecondDefault)
    {
      confess("optionRequired appears to corrupt the multi-item-list value!");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when retrieving list-style option: $@";
  }
  printTestResult(8, $testSuccess);
    
  # Rough test that multi-valued options are being read from the option file 
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    my $optionsFileName = "t_data/test_options_file.xml";
    my %options = ();
    loadOptionsFromFile(\%options, $optionsFileName);

    my $localFirstValue = "search_me.fasta";
    my $localSecondValue = "search_me_too.fasta";
    my @knownDatabaseFileValues = ($localFirstValue, $localSecondValue);
    my @readValues = optionRequired(\%options, "database_file");

    if($readValues[0] ne $localFirstValue ||
       $readValues[1] ne $localSecondValue)
    {
      confess("loadOptionsFromFile failed to load multi-valued item");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when loading option file: $@";
  }
  printTestResult(9, $testSuccess);

  # Rough test that multi-valued options are being written to the option file 
  $testSuccess = $TRUE; # Reset for this test
  eval
  {
    # This should load correctly, right?!
    my $optionsFileName = "t_data/test_options_file.xml";
    my %firstOptions = ();
    loadOptionsFromFile(\%firstOptions, $optionsFileName);

    my $tempFileName = "/tmp/options_test$$.xml";
    writeOptionsToFile(\%firstOptions, $tempFileName);
     
    my %options = ();
    loadOptionsFromFile(\%options, $tempFileName);

    unlink($tempFileName);

    my $localFirstValue = "search_me.fasta";
    my $localSecondValue = "search_me_too.fasta";
    my @knownDatabaseFileValues = ($localFirstValue, $localSecondValue);
    my @readValues = optionRequired(\%options, "database_file");

    if($readValues[0] ne $localFirstValue ||
       $readValues[1] ne $localSecondValue)
    {
      confess("loadOptionsFromFile failed to write multi-valued item");
    }
  };
  if($@)
  {
    $testSuccess = $FALSE;
    print "Error when loading option file: $@";
  }
  printTestResult(10, $testSuccess);
     
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

