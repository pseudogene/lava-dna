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

package LLNL::LAVA::Options;

use strict;
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);

use Carp ":DEFAULT", "cluck"; # For informative error messages

# Requires XML::LibXML for file operations

require Exporter;

@ISA = qw(Exporter);

@EXPORT_OK = ("optionRequired",
              "optionWithDefault",
              "loadOptionsFromFile",
              "writeOptionsToFile");
%EXPORT_TAGS = (standard => ["optionRequired",
                             "optionWithDefault",
                             "loadOptionsFromFile",
                             "writeOptionsToFile"]);

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::Options - A way to enforce options and parameters for functions 

=head1 SYNOPSIS

  use LLNL::LAVA::Options ":options"; 
    # optionRequired, 
    # optionWithDefault, 
    # loadOptionsFromFile, 
    # writeOptionsToFile

  # A couple easy examples

  # Parameter hash is a reference to name->value pairs
  my $optionValue = optionRequired(\%options, "option_name");
  my $optionValue = optionRequired(\%options, "option_name", 
    $informativeErrorMessage);
  my @listValue = optionRequired(\%options, "same_named_options");
  my @listValue = optionRequired(\%options, "same_named_options",
    $informativeErrorMessage);

  my $optionValue = optionWithDefault(\%options, "option_name", $defaultValue);
  my @listValue = optionWithDefault(\%options, "same_named_options", 
    @myDefaultValues);

  # Option file functions
  #
  # Pre-existing option values take precedence when there is a conflict.
  # This is intended to make it simple for command-line arguments to override
  # file-defined option values
  # Yeah, this might not do FULL command-line representation, but it should
  # be good enough to start with :)
  #
  # App::Options is promising as an alternative to the file IO part of this,
  # but its documentation is a little opaque for now.

  # Pass in a hash ref which will be modified
  # If you already know the file name to read from, you can specify it
  loadOptionsFromFile(\%options, $optionsFileName); 

  # If the file name is one of the options, you can use this style where the
  # hash ref MUST CONTAIN the key/value "option_file" which is the file 
  # name to read.
  loadOptionsFromFile(\%options); 

  # Write an options ML file
  writeOptionsToFile(\%options, $outputFileName);  

  # You can also specify write mode if you want, valid modes are:
  #   0 - No stomping of file (default if not supplied)
  #   1 - [NOT IMPLEMENTED] append to file
  #   2 - Stop existing file if it's there
  writeOptionsToFile(\%options, $outputFileName, $writeMode); 

  # Sample option file
  <option_set>
    <sequence_file>virus_sequence.fasta</sequence_file>
    <database_file>search_me.fasta</database_file>
    <word_size>7</word_size>
    <output_file>answer.out</output_file>
  </option_set>

=head1 EXAMPLES

See Synopsis

=head1 DESCRIPTION

LLNL::LAVA::Options defines functions that many perl scripts need

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 optionRequired 

 Usage     : my $optionValue = optionRequired(\%params, "option_name");
             my $optionValue = optionRequired(\%params, "option_name", 
               "informative error message");
 Function  : Will retrieve the value of the element from the passed-in hash
             reference.  If there is no element of that name, or if the value 
             is blank, then will Carp out. 
 Arguments : String - the element name to retrieve the value of
             String - error message to show if parameter loading fails 
 Example   : See Usage
 Returns   : String - value of the argument 

=cut

sub optionRequired
{
  my ($optionHashRef, $optionName, $usageString) = @_;

  if(!defined $usageString)
  {
    $usageString = " ";
  }

  # Crap-out if the option isn't available
  # (Can't die on an empty-list here, because an empty array-ref
  # is a perfectly valid option value in some contexts)
  if(! exists $optionHashRef->{$optionName} ||
     ! defined $optionHashRef->{$optionName} ||
     $optionHashRef->{$optionName} =~ /^\s*$/)
  {
    confess("option \"$optionName\" is required\n$usageString");
  }

  # We made it this far, so we know the value is available
  if(wantarray())
  {
    # List context (return array)
    my @paramValues = ();

    # Read from the option hash, could be an
    # array ref, or a single value (assume single if not array ref, so 
    # objects are scalar values!), either way they need to be returned
    # as an array
    if(ref($optionHashRef->{$optionName}) eq "ARRAY")
    {
      @paramValues = @{$optionHashRef->{$optionName}};
      #print "Array option assigned " .
      #  scalar(@paramValues) .
      #  " values recognized\n";
    }
    else
    {
      # Scalar type, return the scalar value of the element (might actually
      # be a reference) 
      push(@paramValues, $optionHashRef->{$optionName});
      #print "Scalar option assigned\n";
    }

    return @paramValues;
  }
  else
  {
    # Scalar context (return single value)

    return $optionHashRef->{$optionName};
  }
}

#-------------------------------------------------------------------------------

=head2 optionWithDefault 

 Usage     : my $optionValue = optionWithDefault(\%params, "option_name", 
               "default value");
 Function  : Will retrieve the value of the element from the passed-in hash
             reference.  If there is no element of that name, or if the value 
             is blank, then the default parameter value will be returned. 
 Arguments : String - element name to retrieve the value of
             String - default to return in absence of a value
 Example   : See Usage
 Returns   : String - value of the element, or the default value 

=cut

sub optionWithDefault
{
  my ($optionHashRef, $optionName, $firstDefault, @otherDefaults) = @_;

  if(wantarray())
  {
    # List context (return array)
    my @paramValues = ();

    # Return the defaults if not available
    if(! exists $optionHashRef->{$optionName} ||
       ! defined $optionHashRef->{$optionName} ||
       $optionHashRef->{$optionName} =~ /^\s*$/)
    {
      push(@paramValues, $firstDefault);
      if(scalar(@otherDefaults) != 0)
      {
        push(@paramValues, @otherDefaults);
      }

      return @paramValues;
    }

    # Value is avalable, so read from the option hash, could be an array ref,
    # or a single value (always assuming single if not array ref, so 
    # objects are scalar values!), either way, they need to be returned
    # as an array
    if(ref($optionHashRef->{$optionName}) eq "ARRAY")
    {
      @paramValues = @{$optionHashRef->{$optionName}}; 
    }
    else
    {
      # Scalar version - return the "scalar" value of the hash element, which
      # might actually be a reference 
      push(@paramValues, $optionHashRef->{$optionName});
    }

    return @paramValues;
  }
  else
  {
    # Scalar context (return single value)

    # Return default if not available
    if(! exists $optionHashRef->{$optionName} ||
       ! defined $optionHashRef->{$optionName} ||
       $optionHashRef->{$optionName} =~ /^\s*$/)
    {
      return $firstDefault;
    }

    return $optionHashRef->{$optionName};
  }
}

#-------------------------------------------------------------------------------

=head2 loadOptionsFromFile 

 Usage     : loadOptionsFromFile(\%options, $fileName); 
             loadOptionsFromFile(\%options); # This way is OK when %options 
               # has an "option_file" key with the filename as value
 Function  : Loads an option file, filling keys and values into the
             option hash whenever the options weren't previously defined.
             This function will silently do nothing if an option file is not
             specified properly - this is so you don't need to include more
             logic to only check file options when available.  Sorry if it
             causes you heartburn
 Arguments : Hash ref - Hash of name/value pairs 
             String (optional with caveat) - File name to read
 Example   : See Usage
 Returns   : <n/a>

=cut

sub loadOptionsFromFile
{
  my ($optionRef, $optionalFileName) = @_;

  # This use statement SHOULD probably be at the global scope, but I put
  # it at the start of both the option file functions to keep it from being
  # loaded unless actually used.  I hope this means that pre 5.6.1 scripts
  # don't crash when using Options!
  use XML::LibXML;

  # Catch bad parameter error 
  if(! defined $optionRef ||
     ref($optionRef) ne "HASH")
  {
    confess("parameter must be a hash reference of name/value options");
  }

  # Load options from file if available, otherwise toss a warning and
  # return in a giving-up-but-not-quitting fashion
  my $optionFile;
  if(defined $optionalFileName)
  {
    $optionFile = $optionalFileName;
  }
  else
  {
    # Since optional filename argument not provided, a "option_file"
    # element must exist in the options ref in order to do anything useful
    my $optionTitle = "option_file";
    if(! exists($optionRef->{$optionTitle}))
    {
#      cluck("No option file visible: option \"$optionTitle\" must exist in " .
#        "the option set, or a file name must be provided as second argument to " .
#        "loadOptionsFromFile()");

      return;
    }

    $optionFile = $optionRef->{$optionTitle};
  }

  # If file parameter doesn't have a value, toss hands (and warning) into
  # the air and return
  if(! defined $optionFile ||
     $optionFile =~ /^\s*$/)
  {
#    cluck("No option file visible: the supplied option file was not " .
#      "defined, or is the empty string");
    return;
  }

  # Catch non-existent file errors early
  if(! -e $optionFile)
  {
    confess("option file \"$optionFile\" not found");
  }

  # Load the option file
  my $xmlParser = XML::LibXML->new();
  my $xmlDoc = ($xmlParser->parse_file($optionFile))->getDocumentElement();

  # Take the first option set (someday might be able to specify a name or 
  # ID for which option set to use)
  my @setList = $xmlDoc->findnodes("//option_set");
  if(scalar(@setList) < 1)
  {
    confess("failure - options file \"$optionFile\" had no valid \"option_set\" " .
      "in it");
  }

  # Just take the first "option_set" 
  my $optionsNode = $setList[0];
  #my @nameValueNodes = $optionsNode->childNodes();
  my @nameValueNodes = $optionsNode->nonBlankChildNodes();

  # Populate a temporary hash with the options in the file
  my %optionsFromFile;
  foreach my $currNameValueNode(@nameValueNodes)
  {
    my $currOptionName = $currNameValueNode->nodeName();

    # Skip orphan text nodes, only look at sub-nodes
    if($currOptionName eq "text")
    {
      next;
    }

    my $currOptionValue = $currNameValueNode->textContent();
    #print "option from file \"$currOptionName\" has value \"$currOptionValue\"\n";

    # Advance multiply-listed values into array-space by representing
    # the value as an array-ref
    if(! exists $optionsFromFile{$currOptionName})
    {
      # -Doesn't exist - get the value assigned directly
      $optionsFromFile{$currOptionName} = $currOptionValue;      
      #print "File load - \"$currOptionName\" didn't exist, assigning value \"$currOptionValue\"\n";
    }
    else
    {
      if(ref($optionsFromFile{$currOptionName}) eq "ARRAY")
      {
        # -Exists as array - get the value added to the array
        push(@{$optionsFromFile{$currOptionName}}, $currOptionValue);
        #print "File load - \"$currOptionName\" existed, pushing \"$currOptionValue\"\n";
      }
      else
      {
        # -Exists as scalar - change value into array ref with both values
        my $firstOptionValue = $optionsFromFile{$currOptionName};
        $optionsFromFile{$currOptionName} = [$firstOptionValue, $currOptionValue];
        #print "File load - \"$currOptionName\" converted to array to add value \"$currOptionValue\"\n";
      }
    }
  } # End foreach currNameValueNode(@nameValueNodes)

  # Prefer pre-existing option values when there is a conflict - this
  # is intented to make it simple for command-line options to take
  # precedence over file-defined options.
  # This won't work if command-line option specified NO VALUES, because
  # we're using the empty array as a signal to override with file-based values
  foreach my $currOptionName(keys %optionsFromFile)
  {
    if(! exists $optionRef->{$currOptionName} || 
       ! defined $optionRef->{$currOptionName} ||
       (ref($optionRef->{$currOptionName}) && # 3-step check for empty array ref
        ref($optionRef->{$currOptionName}) eq "ARRAY" &&
        scalar(@{$optionRef->{$currOptionName}}) == 0) )
    {
      $optionRef->{$currOptionName} = $optionsFromFile{$currOptionName};
      #print "Allowing \"$currOptionName\" to be read from file\n";
    }
  }

  return;
}

#-------------------------------------------------------------------------------

=head2 writeOptionsToFile 

 Usage     : writeOptionsToFile(\%options, $fileName, $writeMode); 
 Function  : Writes an option file, printing name/value pairs as XML elements 
 Arguments : Hash ref - Hash of name/value pairs 
             String - File name to write
             Integer (optional) - Write mode, valid values are:
                 0 - No stomping of file (default if not supplied)
                 1 - [NOT IMPLEMENTED] append to file
                 2 - Stop existing file if it's there
 Example   : See Usage
 Returns   : <n/a>

=cut

sub writeOptionsToFile
{
  my ($optionRef, $optionFileName, $writeMode) = @_;

  # This use statement SHOULD probably be at the global scope, but I put
  # it at the start of both the option file functions to keep it from being
  # loaded unless actually used.  I hope this means that pre 5.6.1 scripts
  # don't crash when using Options!
  use XML::LibXML;

  # Validate parameters
  if(! defined $optionRef ||
     ref($optionRef) ne "HASH")
  {
    confess("Hash reference of options expected as first argument");
  }
  if(!defined $optionFileName ||
     $optionFileName =~ /^\s*$/)
  {
    confess("File name expected as second argument");
  }

  # Write modes are:
  # 0 - no stomping of files
  # 1[NOT IMPLEMENTED] - append to existing file (or create if nonexistent)
  # 2 - stomp existing file if it's there
  if(! defined $writeMode ||
     ($writeMode != 0 &&
      $writeMode != 1 &&
      $writeMode != 2))
  {
    $writeMode = 0; # Default conservatively
  }

  # Practice file safety kids
  if(-e $optionFileName)
  {
    if($writeMode == 0 ||
       $writeMode == 1)
    {
      confess("Chickening out because asked to stomp the pre-existing " .
        "file \"$optionFileName\"");
    }
  }

  # Build the XML document to write
  my $xmlDoc = XML::LibXML::Document->createDocument();
  my $optionRoot = $xmlDoc->createElement("option_set");
  $xmlDoc->setDocumentElement($optionRoot); # Set up the node as the root

  foreach my $currOption(keys(%{$optionRef}))
  {
    # Blank options weren't defined on the command line, and weren't 
    # specified in the file (unless improperly so)
    if(! defined $optionRef->{$currOption})
    {
      next;
    }

    # Could just use getOptionRequired to return array or scalar, but I
    # wanted descriptive error messages, so I'm doing it manually here

    # Catching non-scalar refs here, only scalar values, and array
    # refs of scalar values are allowed
    my $currOptionValue = $optionRef->{$currOption};
    if( ! ref($currOptionValue))
    {
      # Scalar value, append text child
      $optionRoot->appendTextChild($currOption, $currOptionValue);
    }
    else
    {
      # Is a reference, don't accept non-array references for writing to file
      if(ref($currOptionValue) ne "ARRAY")
      {
        confess("programming error - options ref sent to writeOptionsToFile() " .
          "contained non-scalar non-arrayref value for \"$currOption\"");
      }
      
      # Append all text children - don't accept non-scalar values
      my @allValues = @{$currOptionValue};
      foreach my $nestedValue(@allValues)
      {
        if(ref($nestedValue))
        {
          confess("programming error - options ref sent to writeOptionsToFile() " .
            "contained non-scalar value in arrayref of \"$currOption\"");
        }

        $optionRoot->appendTextChild($currOption, $nestedValue);
      } 
    } 

    # TODO: catch arrayrefs here! write each value out
    # Crap-out if asked to write any type of reference out, only allowing
    # scalar values to be output into option files.

  }

  # The second arg of toFile() is 0 or undef for normal, 1 for human readable
  $xmlDoc->toFile($optionFileName, 0);

  return;
}

1; # Got one?

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
