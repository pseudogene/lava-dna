package LLNL::LAVA::OligoEnumerator::Primer3Conserved;

use strict;
use vars qw(@ISA);
use Carp;

use Bio::SimpleAlign; # Recieves as a parameter
use Bio::Tools::Run::Primer3; # Uses to enumerate oligos over sequences

use LLNL::LAVA::Constants ":standard";
use LLNL::LAVA::Options ":standard";

use LLNL::LAVA::Oligo; # Builds and returns

use LLNL::LAVA::OligoEnumerator; # is-a
@ISA = ("LLNL::LAVA::OligoEnumerator");

# POD-formatted documentation
#-------------------------------------------------------------------------------

=head1 NAME

LLNL::LAVA::OligoEnumerator::Primer3Conserved - Finds all oligos that could be
used as primers based on perfect conservation across the input BioPerl MSA, 
using Primer3 for enumeration.
   
=head1 SYNOPSIS

  use LLNL::LAVA::OligoEnumerator::Primer3Conserved;

  # Instantiation
  $enumerator = LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
    {
      "primer3_executable" => "/path/to/primer3_core",
    });
  
  # Set primer3 operating parameters 
  # ONLY these tags are recognized for now.
  # These values aren't sanity checked, so be very careful what you ask for!
  $enumerator->setPrimer3Targets(
    {
      "target_length" => 20,
      "min_length" => 17,
      "max_length" => 27,
      "target_tm" => 62,
      "min_tm" => 61,
      "max_tm" => 63,
      "most_to_return" => 20001,
    });

  # Get oligo results of enumerating over an MSA
  @oligos = $enumerator->getOligos($bioPerlMSA);

=head1 EXAMPLES

See synopsis.

=head1 DESCRIPTION

Primer3Conserved is a concrete OligoEnumerator.  A BioPerl MSA is accepted
as input to getOligos().  Oligos are enumerated across the first sequence
using BioPerl's Bio::Tools::Run::Primer3.  Oligos that have exactly conserved
matches in all the remaining MSA sequences are returned.

=head1 APPENDIX

The following documentation describes the functions in this package

=cut

#-------------------------------------------------------------------------------

=head2 new

 Usage     : $oligoEnumerator = 
               LLNL::LAVA::OligoEnumerator::Primer3Conserved->new(
                 {"primer3_executable" => "/usr/bin/primer3_core"} );
 Function  : Creates a new OligoEnumerator for Primer3 and perfect consensus
 Arguments : Hash Ref - options for OligoEnumerators, including:
               primer3_executable - path to primer3 (usually primer3_core)
 Example   : See Usage
 Returns   : A new LLNL::LAVA::OligoEnumerator

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

  # Path to the primer3 executable
  my $primer3Executable = optionRequired($paramHash_r, "primer3_executable");

  my $this = $classType->SUPER::new();

  # Map of parameter name to official Primer3 target name
  # No values can be set through setPrimer3Targets() that aren't listed here 
  my $p3Names_r = {
    "target_length" => "PRIMER_INTERNAL_OLIGO_OPT_SIZE",
    "min_length" => "PRIMER_INTERNAL_OLIGO_MIN_SIZE",
    "max_length" => "PRIMER_INTERNAL_OLIGO_MAX_SIZE",
    "target_tm" => "PRIMER_INTERNAL_OLIGO_OPT_TM",
    "min_tm" => "PRIMER_INTERNAL_OLIGO_MIN_TM",
    "max_tm" => "PRIMER_INTERNAL_OLIGO_MAX_TM",
    "max_poly_bases" => "PRIMER_INTERNAL_OLIGO_MAX_POLY_X",
    "most_to_return" => "PRIMER_NUM_RETURN", 
  };

  # Set of default primer3 targets (primer3 target name => value)
  my $p3Targets_r = {
    "PRIMER_TASK" => "pick_hyb_probe_only",
    #"PRIMER_INTERNAL_OLIGO_SELF_ANY" => "12.00",
    "PRIMER_INTERNAL_OLIGO_SELF_ANY" => "8.00",

    "PRIMER_SALT_CONC" => 50, # Milli-molar
    "PRIMER_DNA_CONC" => 50, # Nano-molar
    "PRIMER_INTERNAL_OLIGO_MAX_POLY_X" => 4,
    #"PRIMER_NUM_NS_ACCEPTED" => 0, # Already the default?

    # Not sure if we should make these adjustable or not
    "PRIMER_INTERNAL_OLIGO_MIN_GC" => 20,
    "PRIMER_INTERNAL_OLIGO_MAX_GC" => 80,

    # Default suggested by Primer3 documentation
    #"PRIMER_INTERNAL_OLIGO_MAX_END_STABILITY" => 9.0, 
    #"PRIMER_MAX_END_STABILITY" => 9.0, 
  
    $p3Names_r->{"target_length"} => 20,
    $p3Names_r->{"min_length"} => 18,
    $p3Names_r->{"max_length"} => 27,
    $p3Names_r->{"target_tm"} => 62,
    $p3Names_r->{"min_tm"} => 61,
    $p3Names_r->{"max_tm"} => 63,
    $p3Names_r->{"most_to_return"} => 20001, # Off-by-one error in primer3?
  };

  $this->{"d_primer3Executable"} = $primer3Executable;
  $this->{"d_primer3NameConversion"} = $p3Names_r; 
  $this->{"d_primer3Targets"} = $p3Targets_r;
 
  return $this;
}

#-------------------------------------------------------------------------------

=head2 getOligos

 Usage     : $enumerator->getOligos($bioPerlMSA);
 Function  : Creates an array of LLNL::LAVA::Oligo objects based on the
             input BioPerl MSA.  Oligos are enumerated across the first 
             sequence using BioPerl's Bio::Tools::Run::Primer3.  Oligos that 
             have exactly conserved matches in all the remaining MSA 
             sequences are returned.
 Arguments : Bio::SimpleAlign - MSA of the targets
 Example   : See Usage
 Returns   : Array of L<LLNL::LAVA::Oligo> - the oligos found for the MSA

=cut

sub getOligos
{
  my ($this, $alignment) = @_;

  my $sequenceCount = $alignment->no_sequences();
  if($sequenceCount <= 0)
  {
    confess("data error - MSA must contain at least one sequence");
  }

  # Make sure identifiers are unique, and that lengths are identical
  my %sequenceIDs = (); 
  my $observedLength = -1;
  foreach my $sequence($alignment->each_seq())
  {
    # Make any non-ATCGN an N, should handle this more gracefully, but 
    # primer3 doesn't really like non-ATCGN's
    my $seqContent = $sequence->seq();
    $seqContent =~ s/[^ATCG]/N/g;
    $sequence->seq($seqContent);

    my $id = $sequence->id();
    my $headerPiece = $sequence->desc();
    if(defined $headerPiece &&
       $headerPiece ne "")
    {
      $id .= " " . $headerPiece;
    }

    if(exists $sequenceIDs{$id})
    {
      confess("data error - alignment has 2 entries with the ID \"$id\", but " .
        "we can only handle one sequence of each ID here");
    }
    
    # Check for matched lengths
    my $currLength = $sequence->length();
    if($observedLength == -1) # Record first length no matter what
    {
      $observedLength = $currLength;
    }
    if($observedLength != $currLength)
    {
      confess("data error - the first sequence was $observedLength long, but " .
        "the sequence \"$id\" was $currLength long");
    }

    $sequenceIDs{$id} = $TRUE;
  }

  # TODO: Going to have to scale the number of primers returned along
  # with the sequence length, or enforce a max length for this locus?

  my %oligosBySequenceID = ();
  my $sequenceIndex = 1;
  
  # Running primer3 on the first sequence
  my $firstSequence = $alignment->get_seq_by_pos(1); # 1-indexed!
  my $primer3 = Bio::Tools::Run::Primer3->new(
    -seq => $firstSequence,
    -path => $this->{"d_primer3Executable"});
  my $p3Targets_r = $this->{"d_primer3Targets"};
  #foreach my $targetName (keys(%{$p3Targets_r}))
  #{
  #  print "  $targetName\t" .
  #    $p3Targets_r->{$targetName} .
  #    "\n"; 
  #}

  $primer3->add_targets(%{$p3Targets_r});

  my $primer3Results_r = $primer3->run();
  my $oligoCount = $primer3Results_r->number_of_results;
  print "Primer3Conserved getOligos had $oligoCount oligos\n";

  # Transition the results into an array, can't use the Primer3 object
  # interface, because it only undestands primer pairs, not individual
  # oligos
  my $resultsList_r = [];
  for(my $currIndex = 0; $currIndex < $oligoCount; $currIndex++)
  {
    #print "  Loading $currIndex\n";

    my $result = $primer3Results_r->primer_results($currIndex);

    #print "  Result: \"" .
    #  ref($result) .
    #  "\"\n";

    my $oligo = LLNL::LAVA::Oligo->newFromPrimer3($result);

    #print "  Oligo: \"" .
    #  ref($oligo) .
    #  "\"\n";

    push(@{$resultsList_r}, $oligo);
  }

  # TODO: Add a special primer3 oligo filtering step to remove
  # primers that would violate "PRIMER_INTERNAL_OLIGO_MAX_STABILITY", which
  # can't be specified to primer3 right now.
  # This limits the 5' stability of the last 5 bases of the primer. 
     
  my @survivingOligos = @{$resultsList_r};
  # Starting with the second alignment in the sequence, make sure the
  # oligo is at the same location 
  for(my $i = 2; $i <= $sequenceCount; $i++) # 1-indexed again!
  {
    # Grab the data for this aligned sequence
    my $currSequence = ($alignment->get_seq_by_pos($i))->seq();

    # Check for equivalence at the correct location of all potential primers
    # Store the list of still-valid primers at every step since we expect
    # the number to fall off quickly.
    my @stillSurvivingOligos = ();
    foreach my $oligo(@survivingOligos)
    {
      my $oligoLocation = $oligo->location();
      my $oligoLength = $oligo->length();

      if($oligo->sequence() eq substr($currSequence, $oligoLocation, $oligoLength))
      {
        push(@stillSurvivingOligos, $oligo);
      }
    }

    # Update the survivor list
    @survivingOligos = @stillSurvivingOligos;   
  }

  return @survivingOligos;
}

#-------------------------------------------------------------------------------

=head2 setPrimer3Targets

 Usage     : $enumerator->setPrimer3Targets({target_name => $target_value...});
 Function  : Sets the oligo properties for executing Primer3.  Will cause an
             error if a target is used that this module doesn't have a
             default value for.
 Arguments : Hash ref of target->value pairs
 Example   : See Usage
 Returns   : <n/a>

=cut

sub setPrimer3Targets
{
  my ($this, $paramHash_r) = @_;

  if(!defined $paramHash_r)
  {
    confess("programming error - first parameter is a required hash ref");
  }
  if(ref($paramHash_r) ne "HASH")
  {
    confess("programming error - first parameter must be a hash reference");
  } 

  my $p3Targets_r = $this->{"d_primer3Targets"};
  my $p3Names_r = $this->{"d_primer3NameConversion"}; 

  # for each param passed in
  foreach my $currTarget(keys(%{$paramHash_r}))
  {
    # oops if it doesn't exist
    if(! exists $p3Names_r->{$currTarget})
    {
      confess("programming error - failing to set primer3 target " .
        "\"$currTarget\" because no default value for that target was " .
        "created during Primer3Conserved's instantiation");
    }

    # Since it does exist, set the correct tag in Targets to be the new value
    $p3Targets_r->{$p3Names_r->{$currTarget}} = $paramHash_r->{$currTarget};
  }
}

1; # Lame

__END__

=head1 AUTHOR

Clinton Torres (clinton.torres@llnl.gov)

=head1 SEE ALSO

=cut
