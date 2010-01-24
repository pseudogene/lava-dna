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

package LLNL::LAVA::Constants;

use strict;
use vars qw($VERSION @ISA @EXPORT_OK %EXPORT_TAGS $TRUE $FALSE $debug %options %optionMap);

require Exporter;

@ISA = qw(Exporter);

@EXPORT_OK = qw(*TRUE 
                *FALSE 
                *debug
		*options
		*optionMap
                *VERSION);
%EXPORT_TAGS = 
  (standard => [qw(*TRUE
                   *FALSE
                   *debug
		   *options
		   *optionMap
                   *VERSION)] );

# Reference to a scalar is basicallya Perl constant
*VERSION = \"0.1";

# POD-formatted documentation

=head1 NAME

LLNL::LAVA::Constants - Exports symbols that every Perl script seems to need.

=head1 SYNOPSIS

  # Constants categories available are:
  #  ":standard"

  use LLNL::LAVA::Constants ":standard";
  $debug = $TRUE; # Both symbols were imported to the current namespace

=head1 EXAMPLES 

  # Constants categories available are:
  #  ":standard"

  use LLNL::LAVA::Constants ":standard";
  $debug = $TRUE; Both symbols were imported to the current namespace

=head1 DESCRIPTION

This module was desiged to make available the things that every Perl script
seems to need. 

All symbols imported with a tag group are currently available for 
explicit individual import as well.

Symbols imported with ":standard"

=over 4 

=item * $TRUE

integer value 1 

=item * $FALSE

integer value 0 

=item * $debug

initialized to $FALSE

=back

=cut

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

# Export these values.
# References to scalars can't be modified in Perl, so they're kinda constants
*TRUE  = \1; # Declared in "use vars"
*FALSE = \0; # Declared in "use vars"

# Export debug, default to not debugging 
$debug = $FALSE; # Declared in use vars

#-----------------------------------------------------------------------

=head2 import()

 Usage     : use LLNL::LAVA::Constants ":standard"; 
             Should only be used by Exporter while executing 
             use directives.
 Function  : This method overrides Exporter's import() function.  By
             taking control of the use statement, we can pass 
             arbitrary information to the module to be used in
             static initialization. 
             This sends all of the information passed in directly
             to Exporter::import().
 Example   : use LLNL::LAVA::Constants ":standard";
 Arguments : String  - Tag set for symbol export (currently available 
             set is ":standard") 
 Returns   : <n/a> 

=cut

sub import
{
  my $moduleName = shift; # This appears to be a simple string
  my @argumentList = @_;

  # Build the list of symbols that Exporter should know about
  # Initialize symbols with the module name
  my @symbolsToExporter;
  push(@symbolsToExporter, $moduleName);
  
  foreach my $currArgument(@argumentList)
  {
    push(@symbolsToExporter, $currArgument);
  }

  # Now call Exporter's import
  # We have to adjust the "level" since we don't want the symbols exported 
  # to this scope, but to the caller script 
  my $originalExportLevel = $Exporter::ExportLevel;
  $Exporter::ExportLevel = 1;
  Exporter::import(@symbolsToExporter);
  $Exporter::ExportLevel = $originalExportLevel; 
}

1; # Got one?

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
