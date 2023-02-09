#!/usr/bin/env perl
# Converted from parse_blastx.pl to parse_blast.pl 
# 
# program written by Alejandro Schaffer to convert a blastx 
# output file to a structured file that can be more easily
# parsed. This program is inspired by a previous program
# blast2flat, written by Yuri Wolf for a similar purpose
# of parsing blastpgp output files. Some of the code is
# reused from my parse_vecscreen.pl program in the
# VecScreen_plus_taxonomy package
#
# Modifications in vadr package by Eric Nawrocki
# https://github.com/ncbi/vadr
# 
# Usage: parse_blast.pl --program x --input <blastx output file>
# OR     parse_blast.pl --program n --input <blastn output file>
#        
# vadr 1.5.1 Feb 2023

use strict;
use warnings;
use Getopt::Long qw(:config no_auto_abbrev);

require "sqp_opts.pm";

# different possible states for the finite automaton used to help parse the input file
my $State_Naive                    = 0;
my $State_SeenDefline              = 1;
my $State_SeenSubjectLength        = 2;
my $State_SeenSubjectScore         = 3;
my $State_SeenSubjectIndentities   = 4;
my $State_FoundQuery               = 5;
my $State_InAlignment              = 6;

# different line types in the input file
my $Linetype_BLASTX        = 0;
my $Linetype_match         = 1;
my $Linetype_position      = 2;
my $Linetype_query_name    = 3;
my $Linetype_length        = 4;
my $Linetype_subject_name  = 5;
my $Linetype_score         = 6;
my $Linetype_identities    = 7;
my $Linetype_query_aln     = 8;
my $Linetype_subject_aln   = 9;
my $Linetype_frame         = 10;
my $Linetype_other         = 11;
my $Linetype_searchspace   = 12;
my $Linetype_strand        = 13;

# variables related to command line options
my $program = undef;             #either "x" or "n", argument for --program
my $input_file = undef;          #file with blastx output to parse
my $allow_subject_minus = undef; #usually 1, set to 0 if --splus used

# variables used while parsing input file
my $state = $State_Naive;   #what is the current state of scanning the input
my $linetype = $Linetype_other; #what is the type of this line
my $qdef;                   #query line
my $qacc;                   #accession of the query
my $qlen;                   #length of query
my $nextline;               #one line of file
my $hacc;                   #accession of matching subject
my $slen;                   #length of subject 
my $hlen;                   #number of positions in the alignment
my $nhsp;                   #number of alignments for this match
my $evalue;                 #E-value for one alignment
my $bitscore;               #bit score for one alignment
my $rawscore;               #raw score for one alignment
my $hdef;                   #defline for matching sequence
my $nmatches;               #number of matching sequences
my $subject_start = 0;      # first aligned position in one block of the subject sequence
my $subject_end = 0;        # last aligned position in one block of the subject sequence
my $old_subject_end = 0;    # last aligned position in previous block of the subject sequence
my $query_start = 0;        # first aligned position in one block of the query sequence
my $query_end = 0;          # last aligned position in one block of the query  sequence
my $old_query_end = 0;      # last aligned position in previous block of the query  sequence
my $query_alignment_part;   # string representing the query part of the alignment       
my $subject_alignment_part; # string representing the subject part of the alignment       
my $processing_alignment;   # are we possibly in the middle of processing an alignment
my $overall_subject_start;  # subject start for an alignment that can have multiple blocks
my $overall_subject_end;    # subject end for an alignment that can have multiple blocks
my $overall_query_start;    # query start for an alignment that can have multiple blocks
my $overall_query_end;      # query end for an alignment that can have multiple blocks
my $identities;             # number of identities in the alignment
my $frame;                  # frame of the nucleotide to protein alignment
my $gaps;                   # number of gaps in the alignment
my $first_line_query;       # are we looking for the first line of the query?
my $first_line_subject;     # are we looking for the first line of the subject?
my $query_strand;           # query strand (+ or -) 
my $subject_strand;         # subject strand (+ or -)

my $overall_query_stop_representation = "";     #representation of stop codons in one row of query alignment
my $one_row_query_stop_representation = "";     #representation of stop codons in one row of query alignment
my $overall_subject_stop_representation = "";   #representation of stop codons in one row of subject alignment
my $one_row_subject_stop_representation = "";   #representation of stop codons in one row of subject alignment
my $overall_query_gap_representation = "";      #representation of gaps in one row of query alignment
my $one_row_query_gap_representation = "";      #representation of gaps in one row of query alignment
my $old_query_gap_overhang;                     #how many gaps in the query are at the end of the current line
my $new_query_gap_overhang;                     #how many gaps in the query are at the end of the current line
my $overall_subject_gap_representation = "";    #representation of gaps in one row of subject alignment
my $one_row_subject_gap_representation = "";    #representation of gaps in one row of subject alignment
my $old_subject_gap_overhang = 0;               #how many gaps in the subject are at the end of the current line
my $new_subject_gap_overhang = 0;               #how many gaps in the subject are at the end of the current line
my $maximum_query_gap       = 0;                #largest gap in the query
my $maximum_subject_gap     = 0;                #largest gap in the subject
my $maximum_query_gap_str   = "";               #string representing largest gap in the query, format Q<d1>:S<d2>:-<d3>
my $maximum_subject_gap_str = "";               #string representing largest gap in the subject, format Q<d1>:S<d2>:+<d3>
my $one_line_query_gap = 0;                     #largest gap in query in one alignment line
my $one_line_subject_gap = 0;                   #largest gap in query in one alignment line
my $one_line_query_gap_str = "";                #string representing largest gap in query in one alignment line, format Q<d1>:S<d2>:-<d3>
my $one_line_subject_gap_str = "";              #string representing largest gap in subject in one alignment line, format Q<d1>:S<d2>:+<d3>
my $DEBUG = 0;                                  #are we in debugging mode?



# variables related to command line options, see sqp_opts.pm
my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
# The opt_Add() function is the way we add options to %opt_HH.
# It takes values of for each of the 2nd dim keys listed above.
#       option                  type       default group   requires incompat  preamble-outfile                                         help-outfile
opt_Add("-h",                   "boolean", 0,          0,    undef, undef,    undef,                                                   "display this help",  \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"1"} = "required options";
opt_Add("--program",            "string",  undef,      1,    undef, undef,   "blast program used is blast<s>, <s> must be 'x' or 'n'", "REQUIRED: blast program used is blast<s>, <s> must be 'x' or 'n'",           \%opt_HH, \@opt_order_A);
opt_Add("--input",              "string",  undef,      1,    undef, undef,   "input fasta file",                                       "REQUIRED: input file name <s> with vecscreen output",     \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"2"} = "other options";
opt_Add("--splus",             "boolean",  0,          1,    undef, undef,   "only parse alignments if subject is + strand",           "only parse alignments if subject is plus strand",         \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $all_options_recognized =
    &GetOptions(
  'h'                  => \$GetOptions_H{"-h"},
  'program=s'          => \$GetOptions_H{"--program"},
  'input=s'            => \$GetOptions_H{"--input"},
  'splus'              => \$GetOptions_H{"--splus"});

my $synopsis = "parse_blast.pl :: convert blastx or blastn output file to a more structured intermediate representation\n";
my $usage    = "Usage: parse_blast.pl ";    

# set options in %opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# get option values
$program             = opt_Get("--program", \%opt_HH);
$input_file          = opt_Get("--input",   \%opt_HH);
$allow_subject_minus = opt_Get("--splus", \%opt_HH) ? 0 : 1;

# Die if any of:
# - non-existent option is used
# - any of the required options are not used.
# - -h is used
my $reqopts_errmsg = "";
if(! defined $program)              { $reqopts_errmsg .= "ERROR, --program not used. It is required.\n"; }
if(! defined $input_file)           { $reqopts_errmsg .= "ERROR, --input not used. It is required.\n"; }
if($reqopts_errmsg ne "") { die $reqopts_errmsg; }

# die if program is not x or n
if(($program ne "x") && ($program ne "n")) { 
  die "ERROR with --program <s>, <s> must be \"x\" or \"n\"\n";
}

# open the input and output files
open(INPUT,        "<", $input_file)           or die "Cannot open $input_file for reading\n";

my $keep_going = 1; # set to '0' when we get to the end of the file
my $new_query  = 0; # set to '1' when we read a new query
$qdef = undef;
while($keep_going) { 
  $new_query = 0;
#-----------------------------------------------------------
#	obtain query data (for first query)
#-----------------------------------------------------------
  if(! defined $qdef) { 
    while(<INPUT>){    #read in a line
      chomp;           #remove new line character
      ($qdef) = m/^Query= (.*)$/;   #assign to qdef if line starts with Query
      if($qdef) { 
        # believe it or not, blast sometimes splits your query sequence name over multiple lines, example:
        #Query= HQ693446.1:dnaorg-
        #duplicated:HQ693446.1:1:2690:+:HQ693446.1:1:2690:+/2379-1738,
        #1645-1204
        # for the query: HQ693446.1:dnaorg-duplicated:HQ693446.1:1:2690:+:HQ693446.1:1:2690:+/2379-1738,1645-1204
        my $append_to_query = <INPUT>; 
        while($append_to_query =~ m/\S/) { # read until we hit a blank line
          chomp $append_to_query; 
          $qdef .= $append_to_query;
          $append_to_query = <INPUT>;
        }
      }
      last if($qdef);               #if qdef got defined, exit the loop
    }
  }
  # else $qdef is defined
  $state = $State_FoundQuery;
  
  ($qacc) = ($qdef =~ m/^(\S+)/); #get the accession for this query 
  
  while(<INPUT>){   #read a line of input
    chomp;          #remove newline character
    ($qlen) = m/Length=(\d+)/;  #find length designation
    last if($qlen);                  #if qlength defined, exit loop
    my ($qrest) = m/^\s*(\S.*\S)\s*$/; #otherwise read the query
    $qdef .= " ".$qrest if($qrest);    #keep concatenating $qrest to $qdef
  }
  
#-----------------------------------------------------------
#	output query data
#-----------------------------------------------------------
  print "QACC\t".$qacc."\n";       #a blank line except for QACC on it
  print "QDEF\t".$qdef."\n";     #QDEF followed by the def line
  print "QLEN\t".$qlen."\n";     #QLEN followed by the length of the query
  
  $nmatches = 0;        #number of matches initialized to 0
#-----------------------------------------------------------
#	skip the junk before hits
#-----------------------------------------------------------
  
  $processing_alignment = 0;
  while((! $new_query) && ($keep_going)) { 
    $nextline = <INPUT>; #read a line of input
    chomp($nextline);
    $linetype = determineLineType($nextline,
                                  $Linetype_BLASTX,
                                  $Linetype_match,
                                  $Linetype_position,
                                  $Linetype_query_name,
                                  $Linetype_length,
                                  $Linetype_subject_name,
                                  $Linetype_score,
                                  $Linetype_identities,				  
                                  $Linetype_query_aln,
                                  $Linetype_subject_aln,
                                  $Linetype_searchspace,
                                  $Linetype_other); 
    if ($DEBUG) {
      print "saw $nextline and classified as type $linetype\n";
    }
    if($Linetype_subject_name == $linetype){        #if the line starts with >
      $state = $State_SeenDefline;
      if ($DEBUG) {
        print "set state to $State_SeenDefline $nextline\n";
      }
#      ($hdef) = ($nextline =~ m/>\s*([A-Z]\S+.*)/);
      ($hdef) = ($nextline =~ m/>\s*(\S+.*)/);
      ($hacc) = ($hdef =~ m/^(\S+)/);
      $nhsp = 0;        #number of alignments for one match initialized to 0
      if ($processing_alignment) {
        # finished with previous query/target pair: 
        # output info from final HSP from previous query/target pair
        if ($overall_query_stop_representation) {
          print "QSTOP\t" .$overall_query_stop_representation . "\n";
          print "HSTOP\t" .$overall_subject_stop_representation . "\n";
        }
        if ($overall_query_gap_representation) {
          print "DEL\t" .$overall_query_gap_representation . "\n";
          print "MAXDE\t" .$maximum_query_gap_str ."\n";
        }
        if ($overall_subject_gap_representation) {
          print "INS\t" .$overall_subject_gap_representation . "\n";
          print "MAXIN\t" .$maximum_subject_gap_str ."\n"; 		
        }	    	    
        print "SRANGE\t".$overall_subject_start."..".$overall_subject_end."\n"; #print SRANGE and the match range
        print "QRANGE\t".$overall_query_start."..".$overall_query_end."\n"; #print QRANGE and the query range
        # finished with previous query/target pair: re-initialize per-HSP variables
        $processing_alignment = 0;
        $query_strand = "";
        $subject_strand = "";
        $overall_query_stop_representation = "";
        $one_row_query_stop_representation = "";	    
        $overall_subject_stop_representation = "";
        $one_row_subject_stop_representation = "";	    
        $overall_query_gap_representation = "";
        $one_row_query_gap_representation = "";
        $overall_subject_gap_representation = "";
        $one_row_subject_gap_representation = "";
        $maximum_query_gap = 0;	    
        $maximum_subject_gap = 0;
        $maximum_query_gap_str = "";	    
        $maximum_subject_gap_str = "";
        $query_start     = 0;
        $query_end       = 0;
        $old_query_end   = 0;
        $subject_start   = 0;
        $subject_end     = 0;
        $old_subject_end = 0;
      }	
    }
    if (($state == $State_SeenDefline) && ($Linetype_length == $linetype)) {
      ($slen) =  ($nextline =~ m/Length=(\d+)/);  #find length designation
      $nmatches++;            #increment the number of matches
      if ($nmatches > 1) {
        print "END_MATCH\n";  #print END_MATCH
      }
    }
    #    if (($state == $State_InAlignment) && ($linetype == $Linetype_score) && ($nhsp == 0)) {
    if (($linetype == $Linetype_score) && ($nhsp == 0)) {    
      print "MATCH\t".$nmatches."\n";  #print HIT and the number of the hit
      print "HACC\t".$hacc."\n";   #print HACC and the accession
      print "HDEF\t".$hdef."\n"; #print HDEF and the defline
      print "SLEN\t".$slen."\n"; #print SLEN and the subject length
    }
    if ($linetype == $Linetype_score) {
      ($bitscore, $rawscore, $evalue) = getScoreExpect($nextline);
      $nhsp++;            #increment the number of matches
      if ($processing_alignment) {
        if ($DEBUG) {
          print "Going to print alignment summary with overal_subject_representation as $overall_subject_gap_representation\n";
        }
        # finished with previous HSP 
        # output info from previous HSP
        if ($overall_query_stop_representation) {
          print "QSTOP\t" .$overall_query_stop_representation . "\n";
          print "HSTOP\t" .$overall_subject_stop_representation . "\n";
        }
        if ($overall_query_gap_representation) {
          print "DEL\t" .$overall_query_gap_representation . "\n";
          print "MAXDE\t" .$maximum_query_gap_str ."\n"; 
        }
        if ($overall_subject_gap_representation) {
          print "INS\t" .$overall_subject_gap_representation . "\n";
          print "MAXIN\t" .$maximum_subject_gap_str ."\n"; 				
        }
        print "SRANGE\t".$overall_subject_start."..".$overall_subject_end."\n"; #print SRANGE and the match range
        print "QRANGE\t".$overall_query_start."..".$overall_query_end."\n"; #print QRANGE and the query range
        # finished with previous HSP: re-initialize per-HSP variables
        $processing_alignment = 0;
        $query_strand = "";
        $subject_strand = "";
        $overall_query_stop_representation = "";
        $one_row_query_stop_representation = "";	    
        $overall_subject_stop_representation = "";
        $one_row_subject_stop_representation = "";	    
        $overall_query_gap_representation = "";
        $one_row_query_gap_representation = "";
        $overall_subject_gap_representation = "";
        $one_row_subject_gap_representation = "";	    
        $maximum_query_gap = 0;
        $maximum_subject_gap = 0;
        $maximum_query_gap_str = "";
        $maximum_subject_gap_str = "";
        $query_start     = 0;
        $query_end       = 0;
        $old_query_end   = 0;
        $subject_start   = 0;
        $subject_end     = 0;
        $old_subject_end = 0;
      }
      print "HSP\t".$nhsp."\n";          #print HSP and the number of the hsp
      print "BITSCORE\t".$bitscore."\n"; #print SCORE and the bit score
      print "RAWSCORE\t".$rawscore."\n"; #print SCORE and the raw score
      print "EVALUE\t".$evalue."\n";     #print EVALUE and the e-value
    }
    if ($linetype == $Linetype_identities) {
      ($hlen, $identities, $gaps) = getIdentitiesGaps($nextline);
      $first_line_query = $first_line_subject = 1;
    }
    if ($linetype == $Linetype_frame) {
      if($program eq "n") { die "ERROR read Frame line: $nextline but program is blastn (--program n), this line is only expected in blastx output"; }
      $frame = getFrame($nextline);	
      $query_strand   = ($frame =~ m/^\+/) ? "+" : "-";
      $subject_strand = "+"; # for blastx, subject is a protein so always "+"
      print "HLEN\t".$hlen."\n";   #print HLEN and the alignment length
      print "IDENT\t".$identities."\n";   #print IDENTITIES and the number of identities
      print "GAPS\t".$gaps."\n";   #print GAPS and the number of gaps
      print "FRAME\t".$frame."\n";   #print FRAME and the frame of the nucleotide to protein alignment
      $new_query_gap_overhang = $old_query_gap_overhang = 0;
      $one_line_query_gap = $one_line_subject_gap = 0;
      $one_line_query_gap_str = $one_line_subject_gap_str = "";
    }
    if ($linetype == $Linetype_strand) {
      if($program eq "x") { die "ERROR read Strand line: $nextline but program is blastx (--program x), this line is only expected in blastn output"; }
      ($query_strand, $subject_strand) = getStrand($nextline);	
      print "HLEN\t".$hlen."\n";   #print HLEN and the alignment length
      print "IDENT\t".$identities."\n";   #print IDENTITIES and the number of identities
      print "GAPS\t".$gaps."\n";   #print GAPS and the number of gaps
      print "QSTRAND\t".$query_strand."\n";   #print QSTRAND and the query strand of the HSP
      print "SSTRAND\t".$subject_strand."\n"; #print SSTRAND and the subject strand of the HSP
      $new_query_gap_overhang = $old_query_gap_overhang = 0;
      $one_line_query_gap = $one_line_subject_gap = 0;
      $one_line_query_gap_str = $one_line_subject_gap_str = "";
    }
    if ($linetype == $Linetype_query_aln) {
      if (! $processing_alignment) {
        $processing_alignment = 1;
      }
      # Query alignment line:line that begins with "Query "including a segment of the query in an alignment
      ($query_start, $query_end, $query_alignment_part) = getAlignmentPositions($nextline, $old_query_end);
      $old_query_end = $query_end;
      if ($first_line_query) {
        $overall_query_start = $query_start;
        $overall_query_end = $query_end;
        $first_line_query = 0;
      }
      else {
        updateOverallPositions($processing_alignment, $query_start, $query_end, \$overall_query_start, \$overall_query_end);
      }
    }
    if($linetype == $Linetype_subject_aln) {
      # Subject alignment line: line that begins with "Sbjct "  including a segment of the vector (subject) in an alignment
      ($subject_start, $subject_end, $subject_alignment_part) = getAlignmentPositions($nextline, $old_subject_end);
      $old_subject_end = $subject_end;
      if ($first_line_subject) {
        $overall_subject_start = $subject_start;
        $overall_subject_end = $subject_end;
        $first_line_subject = 0;
      }
      else {
        updateOverallPositions($processing_alignment, $subject_start, $subject_end, \$overall_subject_start, \$overall_subject_end);
      }
      ($one_row_query_stop_representation, $one_row_subject_stop_representation) = 
          ($program eq "x") ? findQueryStops($query_alignment_part, $subject_alignment_part, $query_start, $query_end, $subject_start, $subject_end, $query_strand) : "";

      if ($one_row_query_stop_representation) {
        $overall_query_stop_representation   .= $one_row_query_stop_representation;
        $overall_subject_stop_representation .= $one_row_subject_stop_representation;
      }
      if ($gaps > 0) {
        ($new_query_gap_overhang, $one_row_query_gap_representation, $one_line_query_gap_str, $one_line_query_gap) = findQueryGaps($query_alignment_part, $subject_alignment_part, $query_start, $query_end, $subject_start, $subject_end, $old_query_gap_overhang, $program, $query_strand);
        if ($DEBUG) {
          print "From findQueryGaps returned $one_row_query_gap_representation $new_query_gap_overhang\n";
        }
        if ($one_line_query_gap > $maximum_query_gap) {
          $maximum_query_gap     = $one_line_query_gap;
          $maximum_query_gap_str = $one_line_query_gap_str;
        }
        if (($overall_query_gap_representation) && ($one_row_query_gap_representation)) {
          if ($subject_strand eq "+") { 
            $overall_query_gap_representation = $overall_query_gap_representation . ";" . $one_row_query_gap_representation;
          }
          else {
            $overall_query_gap_representation = $one_row_query_gap_representation . ";" . $overall_query_gap_representation;
          }
        }
        else {
          if ($subject_strand eq "+") { 
            $overall_query_gap_representation = $overall_query_gap_representation . $one_row_query_gap_representation;
          }
          else {
            $overall_query_gap_representation = $one_row_query_gap_representation . $overall_query_gap_representation;
          }
        }
        $old_query_gap_overhang = $new_query_gap_overhang;
        if(($subject_strand eq "+") || ($allow_subject_minus)) { 
          ($new_subject_gap_overhang, $one_row_subject_gap_representation, $one_line_subject_gap_str, $one_line_subject_gap) = findSubjectGaps($query_alignment_part, $subject_alignment_part, $query_start, $query_end, $subject_start, $subject_end, $old_subject_gap_overhang, $program, $query_strand, $subject_strand);
          if ($one_line_subject_gap > $maximum_subject_gap) {
            $maximum_subject_gap     = $one_line_subject_gap;
            $maximum_subject_gap_str = $one_line_subject_gap_str;
          }
          if (($overall_subject_gap_representation) && ($one_row_subject_gap_representation)) {
            $overall_subject_gap_representation = $overall_subject_gap_representation . ";" . $one_row_subject_gap_representation;
          }
          else {
            $overall_subject_gap_representation = $overall_subject_gap_representation . $one_row_subject_gap_representation;
          }
        }
        if ($DEBUG) {
          print "Updated subject gap representation is $overall_subject_gap_representation\n";
        }
        $old_subject_gap_overhang = $new_subject_gap_overhang;	    
      }
    }
    if($linetype == $Linetype_searchspace) { 
      # end of query, look for next query
#-----------------------------------------------------------
#	obtain query data (for query #2 or greater)
#-----------------------------------------------------------
      $qdef = undef;
      while(<INPUT>){    #read in a line
        chomp;           #remove new line character
        ($qdef) = m/^Query= (.*)$/;   #assign to qdef if line starts with Query
        if($qdef) { 
          $new_query = 1; 
          # believe it or not, blast sometimes splits your query sequence name over multiple lines, example:
          #Query= HQ693446.1:dnaorg-
          #duplicated:HQ693446.1:1:2690:+:HQ693446.1:1:2690:+/2379-1738,
          #1645-1204
          # for the query: HQ693446.1:dnaorg-duplicated:HQ693446.1:1:2690:+:HQ693446.1:1:2690:+/2379-1738,1645-1204
          my $append_to_query = <INPUT>; 
          while($append_to_query =~ m/\S/) { # read until we hit a blank line
            chomp $append_to_query; 
            $qdef .= $append_to_query;
            $append_to_query = <INPUT>;
          }
        }
        last if($qdef);               #if qdef got defined, exit the loop
      }
      if(! defined $qdef) { # reached end of file
        $keep_going = 0;
      }
    }
  }
  # finished with previous query
  # output info from final HSP from final query/target pair for previous query
  if ($overall_query_stop_representation) {
    print "QSTOP\t" .$overall_query_stop_representation . "\n";
    print "HSTOP\t" .$overall_subject_stop_representation . "\n";
  }
  if ($overall_query_gap_representation) {
    print "DEL\t" .$overall_query_gap_representation . "\n";
    print "MAXDE\t" .$maximum_query_gap_str ."\n";     
  }
  if ($overall_subject_gap_representation) {
    print "INS\t" .$overall_subject_gap_representation . "\n";
    print "MAXIN\t" .$maximum_subject_gap_str ."\n"; 		    
  }
  if(! defined $overall_query_start)   { $overall_query_start = ""; }
  if(! defined $overall_query_end)     { $overall_query_end   = ""; }
  if(! defined $overall_subject_start) { $overall_subject_start = ""; }
  if(! defined $overall_subject_end)   { $overall_subject_end   = ""; }
  print "SRANGE\t".$overall_subject_start."..".$overall_subject_end."\n"; #print SRANGE and the match range
  print "QRANGE\t".$overall_query_start."..".$overall_query_end."\n"; #print QRANGE and the query range
  print "END_MATCH\n";  #print END_MATCH
  # finished with previous query: re-initialize per-HSP variables
  $overall_query_stop_representation   = "";
  $one_row_query_stop_representation   = "";	    
  $overall_subject_stop_representation = "";
  $one_row_subject_stop_representation = "";	    
  $overall_query_gap_representation    = "";
  $one_row_query_gap_representation    = "";
  $overall_subject_gap_representation  = "";
  $one_row_subject_gap_representation  = "";
  $maximum_query_gap = 0;	    
  $maximum_subject_gap = 0;
  $maximum_query_gap_str = "";
  $maximum_subject_gap_str = "";
  $query_start = 0;
  $query_end   = 0;
  $subject_start = 0;
  $subject_end   = 0;
}
################################################
# List of subroutines:
#
# determineLineType();       Returns type of line that $line is
# getScoreExpect();          Extracts the bit and raw scores and the expect value from a score line
# getIdentitiesGaps();       Extracts the numbere of identities and gaps from a suitable line
# getAlignmentPositions();   Extracts two numerical positions from a query or subject alignment line.
# getFrame();                Extracts the frame of a nucleotide to protein alignment from a suitable line
# getStrand();               Extracts the query and subject strands of a nucleotide to nucleotide alignment from a suitable line
# updateOverallPositions();  Update the $$overall_start_R and $$overall_end_R in an alignment.
# findQueryGaps();           Find the positions of any gaps in the query, which represent deletions (in-frame ones if program eq "x")
# findSubjectGaps();         Find the positions of any gaps in the subject, which represent insertions (in-frame ones if program eq "x")
# findQueryStops();          Find the positions of any predicted stop codons in the query. 
#
################################################



##########################################################################################
# Subroutine: determineLineType()
#
# Synopsis: Returns type of line that $line is.
#
# Args: $line:                  the line of input file
#       $Linetype_BLASTX:       '0' (hard-coded 'line type' value for BLASTN lines
#       $Linetype_match:        '1' (hard-coded 'line type' value for match lines
#       $Linetype_position:     '2' (hard-coded 'line type' value for position lines
#       $Linetype_query_name:   '3' (hard-coded 'line type' value for query name lines
#       $Linetype_length:       '4' (hard-coded 'line type' value for length lines
#       $Linetype_subject_name:  '5' (hard-coded 'line type' value for vector name lines
#       $Linetype_score:        '6' (hard-coded 'line type' value for score lines
#       $Linetype_identities:   '7' (hard-coded 'line type' value for score lines
#       $Linetype_query_aln:    '8' (hard-coded 'line type' value for query alignment lines
#       $Linetype_subject_aln:   '9' (hard-coded 'line type' value for vector alignment lines
#       $Linetype_other:        '10' (hard-coded 'line type' value for other (non-parsed) lines
#
# Returns: line type of $line
#
# Dies: Never
#
##########################################################################################

##########################################################################################
sub determineLineType {
    my $sub_name = "determineLineType()";
    my $nargs_exp = 13;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($line, $Linetype_BLASTX, $Linetype_match, $Linetype_position, $Linetype_query_name,
	$Linetype_length, $Linetype_subject_name, $Linetype_score, $linetypeidentitis, $Linetype_query_aln,
	$Linetype_vector_aln, $Linetype_searchspace, $Linetype_other) = @_;

    if($line =~ m/^BLASTX/)                              { return $Linetype_BLASTX;      }
    if($line =~ m/^Strong/)                              { return $Linetype_match;       }
    if($line =~ m/^Moderate/)                            { return $Linetype_match;       }
    if($line =~ m/^Weak/)                                { return $Linetype_match;       }
    if($line =~ m/^Suspect/)                             { return $Linetype_match;       }
    if($line =~ m/^\d/)                                  { return $Linetype_position;    }
    if($line =~ m/^Query=/)                              { return $Linetype_query_name;  }
    if($line =~ m/^Length=\d+/)                          { return $Linetype_length;      }
#    if($line =~ m/^>\s*[A-Z]/)                           { return $Linetype_subject_name; }
    if($line =~ m/^>\s*\S+/)                           { return $Linetype_subject_name; }
    if(($line =~ m/^\s*Score\s*=\s+\S+\s+bits\s+\(\d+\)/) && ($line =~m/Expect/))
      { return $Linetype_score;       }
    if($line =~ m/^\s*Identities\s*=\s+\S+/)             { return $Linetype_identities;  }    
    if($line =~ m/^Query\s+\d+\s+\D+\d+/)                { return $Linetype_query_aln;   }
    if($line =~ m/^Query\s+\-+\s*$/)                     { return $Linetype_query_aln;   } # all gaps in query
    if($line =~ m/^Sbjct\s+\d+\s+\D+\d+/)                { return $Linetype_subject_aln; } 
    if($line =~ m/^Sbjct\s+\-+\s*$/)                     { return $Linetype_subject_aln; } # all gaps in subject
    if($line =~ m/^\s*Frame\s*=\s+\S+/)                  { return $Linetype_frame;       }        
    if($line =~ m/^\s*Strand\s*=\S+/)                    { return $Linetype_strand;      }        
    if($line =~ m/^Effective search space used/)         { return $Linetype_searchspace; }        
    if ($DEBUG) {
	print "In determineLineType; about to return other, which is $Linetype_other \n";
    }
    return $Linetype_other;
}
##########################################################################################
# Subroutine: getScoreExpect()
# Synopsis:   Extracts the bit and raw scores and the expect value from a
#             score line.
#
# Args: $line: score line of input file
#
# Returns: score and expect value
#
# Dies: if $line is not of the type that has the Score and Expect value;
#
##########################################################################################
sub getScoreExpect {
    my $sub_name = "getScoreExpect()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_line) = @_;
    my $local_bit_score = "BLASTNULL"; #stores the bit score
    my $local_raw_score = "BLASTNULL"; #stores the raw score
    my $local_expect    = "BLASTNULL"; #stores the expect value
     
    if($local_line =~ m/^\s*Score\s*=\s+(\S+)\s+bits\s+\((\d+)\)/) {
      ($local_bit_score, $local_raw_score) = ($1, $2);
      ($local_expect) = ($local_line =~ m/Expect\S*\s+=\s+(\S+)\,?/);
      return ($local_bit_score, $local_raw_score, $local_expect);
    }
    else {
      die "ERROR in $sub_name, unexpected format of line: $local_line\n";
    }
}

##########################################################################################
# Subroutine: getIdentitiesGaps()
# Synopsis:   Extracts the numbers of identities and gaps from a suitable line
#             
#
# Args: $line: identities and gaps line of input file
#
# Returns: number of positions in alignment, number of identities and number of gaps
#
# Dies: if $line is not of the type that has the Identities and Gaps;
#
##########################################################################################
sub getIdentitiesGaps {
    my $sub_name = "getIdentitiesGaps()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_line) = @_;
    my $local_length;     #stores the alignment length
    my $local_identities; #stores the number of identities
    my $local_gaps;       #stores the number of gaps
     

    if($local_line =~ m/^\s*Identities\s*=\s*(\d+)/) {
	$local_identities = $1;
	($local_length) = ($local_line =~ m/^\s*Identities\s*=\s*\d+\/(\d+)/);	
	($local_gaps) = ($local_line =~ m/Gaps\s*=\s*(\d+)/);
	if ($DEBUG) {
	    print "extracted local_gaps as $local_gaps\n";
	}
	return ($local_length,$local_identities,$local_gaps) ;
    }
    else {
	die "ERROR in $sub_name, unexpected format of line: $local_line\n";
    }
}

##########################################################################################
# Subroutine: getFrame()
# Synopsis:   Extracts the frame of a nucleotide to protein alignment from a suitable line
#             
#
# Args: $line: frame line of input file
#
# Returns: frame as a string
#
# Dies: if $line is not of the type that has the Frame;
#
##########################################################################################
sub getFrame {
    my $sub_name = "getFrame()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_line) = @_;
    my $local_frame;      #stores the frame
     

    if($local_line =~ m/^\s*Frame\s*=\s*(\S+)/) {
	$local_frame = $1;
	return ($local_frame) ;
    }
    else {
	die "ERROR in $sub_name, unexpected format of line: $local_line\n";
    }
}

##########################################################################################
# Subroutine: getStrand()
# Synopsis:   Extracts the query and subject strands of a nucleotide to nucleotide alignment from a suitable line
#             
# Args: $line: strand line of input file
#
# Returns: Two values:
#          query strand: "+" or "-"
#          subject strand: "+" or "-"
#
# Dies: if $line is not of the type that has the Strand;
#
##########################################################################################
sub getStrand {
    my $sub_name = "getStrand()";
    my $nargs_exp = 1;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

    my ($local_line) = @_;
    my $local_query_strand;      #stores the query strand
    my $local_subject_strand;    #stores the subject strand

    if($local_line =~ m/^\s*Strand\s*=\s*Plus\/Plus\s*/) {
      return("+", "+"); 
    }
    elsif($local_line =~ m/^\s*Strand\s*=\s*Plus\/Minus\s*/) {
      return("+", "-"); 
    }
    elsif($local_line =~ m/^\s*Strand\s*=\s*Minus\/Plus\s*/) {
      return("-", "+"); 
    }
    elsif($local_line =~ m/^\s*Strand\s*=\s*Minus\/Minus\s*/) {
      return("-", "-"); 
    }
    else {
      die "ERROR in $sub_name, unexpected format of line: $local_line\n";
    }
}

##########################################################################################
# Subroutine: getAlignmentPositions()
# Synopsis:   Extracts two numerical positions and
#             the alignment string from
#             a query or subject alignment line.
#             In rare cases, the alignment string can
#             be all gaps in which case no numerical positions
#             are included. This function can deal with that
#             by returning $prv_end as return values for both 
#             $start and $end.
#
# Args: $line: query alignment or vector alignment line of
#              input file
#       $prv_end: end position of previous alignment string
#      
# Returns: Three values:
#          $start: the start position
#          $end:   the end position
#          $aln:   alignment string
#
# Dies: if $line is not a query or vector alignment line;
#       does not match either of
#       /Query\s+(\d+)\s+(\S+)\s+(\d+)/ or
#       /Sbjct\s+(\d+)\s+(\S+)\s+(\d+)/ or
#       /Query\s+(\-+)\s* or (all gaps in query)
#       /Sbjct\s+(\-+)\s* or (all gaps in subject)
#
##########################################################################################
sub getAlignmentPositions {
    my $sub_name = "getAlignmentPositions()";
    my $nargs_exp = 2;
    if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
    
    my ($line, $prv_end) = @_;

    if($DEBUG) { 
      printf("in getAlignmentPositions(), line: $line, prv_end: $prv_end\n");
    }

    if($line =~ m/^Query\s+(\d+)\s+(\S+)\s+(\d+)/) {
      return ($1, $3, $2);
    }
    elsif($line =~ m/^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)/) {
      return ($1, $3, $2);
    }
    elsif($line =~ m/^Query\s+(\-+)\s*/) {
      return ($prv_end, $prv_end, $1);
    }
    elsif($line =~ m/^Sbjct\s+(\-+)\s*/) {
      return ($prv_end, $prv_end, $1);
    }
    else {
	die "in $sub_name(), unexpected format for input line: $line\n";
    }
}

##########################################################################################
# Subroutine:  updateOverallPositions()
#
# Purpose:     Update the $$overall_start_R and $$overall_end_R
#              which keep track of overall alignment positions
#              (start of full alignment and end of full alignment)
#              based on $start and $end (positions of current
#              alignment block). If ($processing_alignment == 1)
#              then this is not the first block, else it is
#              the first block.
#
# Arguments:
#   $processing_alignment: '1' if not first alignment block
#   $start:                start for current block we are processing
#   $end:                  end  for current block we are processing
#   $overall_start_R:      reference, overall start position
#   $overall_end_R:        reference, overall end position
#
# Returns:    void, but updates $$overall_start_R and $$overall_end_R
#
# Dies:       Never
#
##########################################################################################
sub updateOverallPositions {
    my $sub_name = "updateOverallPositions()";
    my $nargs_expected = 5;
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

    my ($processing_alignment, $start, $end, $overall_start_R, $overall_end_R) = @_;

    if ($processing_alignment) {
	# not the first Query Alignment line we've seen for this match,
	# update $overall_end, but not $overall_start
	if (1 == abs($start - $$overall_end_R)) {
	    $$overall_end_R = $end;
	}
    }
    else {
	# the first alignment line we've seen for this match
	$$overall_start_R = $start;
	$$overall_end_R   = $end;
    }
    return;
}

##########################################################################################
# Subroutine:  findQueryGaps()
#
# Purpose:     Find the positions of any gaps in the query, which represent in-frame deletions
#              in the query. Return any completed gap intervals as gap start and the number of
#               of positions  and also return the number of gap characters at the end of the line
# Arguments:
#   $query_alignment_string:     a subsequence of the query, possibly containing gaps
#   $subject_alignment_string:   a subsequence of the query, possibly containing gaps
#   $query_start:                query start position for current block we are processing
#   $query_end:                  query end  position for current block we are processing
#   $subject_start:              subject start position for current block we are processing
#   $subject_end:                subject end  position for current block we are processing
#   $old_gap_overhang:           number of gap positions overhanging from previous row
#   $program:                    from --program, "x" if blastx, "n" if blastn
#   $query_strand:               "+" if query is positive strand, "-" if minus strand
#
# Returns:    Four values:
#             1) the number of gap characters at the end of the row 
#             2) a string of gap-length pairs in query or "" in >= 1 ';' delimited tokens
#                each token in format format: "Q<d1>:S<d2>+<d3>" 
#                <d1> is query position, <d2> is subject position, <d3> is length in nt
#             3) string reporting largest gap on row in same format as 2, but only 1 token
#             4) largest gap length in nt (integer)
#
# Dies:       Never
#
##########################################################################################
sub findQueryGaps {
  my $sub_name = "findQueryGaps()";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($query_alignment_string, $subject_alignment_string, $query_start, $query_end, $subject_start, $subject_end, $old_gap_overhang, $program, $query_strand) = @_;

  if($DEBUG) { 
    print("in $sub_name,\n\tquery_alignment_string: $query_alignment_string\n\tsubject_alignment_string: $subject_alignment_string\n");
  }

  if((! defined $query_strand) || 
     (($query_strand ne "+") && ($query_strand ne "-"))) { 
    die "ERROR read in $sub_name, strand is not + or - (or is undefed"; 
  }

  my $local_residue_qlen = ($program eq "x") ? 3 : 1; # length of one query position in nt in query (3 if blastx, else 1 (blast))
  my $local_length = length($query_alignment_string);  #length of this alignment block
  my @local_query_array = split("",$query_alignment_string);
  my @local_subject_array = split("", $subject_alignment_string);
  my $local_i;                      #loop index
  my $local_state;                  #state of a two-state automaton
  my $local_StateNotGap = 0;        #one state
  my $local_StateInGap = 1;         #the other satet
  my $local_QueryPositionIndex;     #where are we in absolute query positions
  my $local_SubjectPositionIndex;   #where are we in absolute subject positions
  my $local_QueryStartDeletion;     #absolute position in query starting a gap in the query
  my $local_QueryEndDeletion;       #absolute position in query ending a gap in the query
  my $local_SubjectStartDeletion;   #absolute position in subject starting a gap in the subject
  my $local_OneDeletionString;      #string to represent one gap
  my $local_AllDeletionString = ""; #string to represent all gaps        
  my $local_new_gap_overhang = $old_gap_overhang;  #number of gap characters at the end of this algnment row
  my $local_longest_gap = 0;
  my $local_longest_gap_str = "0"; # string representing longest gap "Q<d1>:S<d2>+<d3>" <d1> is query position, <d2> is subject position, <d3> is length in nt
  
  $local_QueryPositionIndex = $query_start;
  $local_SubjectPositionIndex = $subject_start;
  if ($DEBUG) {
    print "In $sub_name: Subject start is $subject_start, old_gap_overhang is: $old_gap_overhang query_start: $query_start query_end: $query_end\n";
  }
  if ($query_strand eq "+") { 
    if ($old_gap_overhang) {
      $local_QueryStartDeletion = $query_start - 1;
      $local_SubjectStartDeletion = $subject_start - $old_gap_overhang - 1;
      $local_state = $local_StateInGap;
    }
    else {
      $local_state = $local_StateNotGap;	    
    }

    # go through each position and keep track of deletions (gaps in query aligned to residues in sbjct)) 
    for ($local_i = 0; $local_i < $local_length; $local_i++) {
      my $query_is_gap = ($local_query_array[$local_i]   eq "-") ? 1 : 0;
      my $sbjct_is_gap = ($local_subject_array[$local_i] eq "-") ? 1 : 0;

      if (! $query_is_gap) { 
        # query is a residue, subject is either a gap or residue
        # this is either a continutation of non-gap in query or end of gap in query
        if ($local_StateNotGap == $local_state) {
          # previous position was a residue in the query, do nothing
          ;
        }
        else {
          # previous position was a gap in query, so this is the end of a deletion
          $local_QueryEndDeletion = $local_QueryPositionIndex - 1;
          $local_OneDeletionString = "Q" . $local_QueryStartDeletion . ":" . "S" . $local_SubjectStartDeletion . "-" . ($local_residue_qlen * $local_new_gap_overhang);

          # if this is the max length deletion we've seen, update the max length deletion
          if (($local_residue_qlen * $local_new_gap_overhang) > $local_longest_gap) {
            $local_longest_gap = $local_residue_qlen * $local_new_gap_overhang;
            $local_longest_gap_str = $local_OneDeletionString;
          }
          if ($local_AllDeletionString) {
            $local_AllDeletionString = $local_AllDeletionString . ";" . $local_OneDeletionString;
          }
          else {
            $local_AllDeletionString = $local_OneDeletionString;			
          }
          $local_state = $local_StateNotGap;
          $local_new_gap_overhang = 0;
        }
        if(! $sbjct_is_gap) { 
          $local_SubjectPositionIndex++;
        }
        $local_QueryPositionIndex += $local_residue_qlen; # 1 (blastn) or 3 (blastx)
      }
      elsif($query_is_gap && (! $sbjct_is_gap)) { 
        # query is a gap, sbjct is a residue
        if ($local_StateNotGap == $local_state) {
          # previous position was a residue in the query, so this is the start of a query gap (deletion)
          if ($DEBUG) {
            print "In findQueryGaps: start of a gap in the query at subject $local_SubjectPositionIndex\n";
          }
          $local_new_gap_overhang = 1;
          $local_QueryStartDeletion = $local_QueryPositionIndex-1; 
          $local_SubjectStartDeletion = $local_SubjectPositionIndex-1; 
          $local_state = $local_StateInGap;
        }
        else {
          # previous position was a gap in the query, so this is a continutation of a query gap, extend it
          $local_new_gap_overhang++;
        }
        $local_SubjectPositionIndex += 1;
      }
    }
    return($local_new_gap_overhang,$local_AllDeletionString, $local_longest_gap_str, $local_longest_gap);
  } # end of if($query_start eq "+"
  else { 
    # $query_strand is "-", this block is analogous to "+" strand block
    if ($old_gap_overhang) {
      $local_QueryStartDeletion = $query_start + 1;
      $local_SubjectStartDeletion = $subject_start - $old_gap_overhang - 1;
      $local_state = $local_StateInGap;
    }
    else {
      $local_state = $local_StateNotGap;	    
    }
    for ($local_i = 0; $local_i < $local_length; $local_i++) {
      my $query_is_gap = ($local_query_array[$local_i]   eq "-") ? 1 : 0;
      my $sbjct_is_gap = ($local_subject_array[$local_i] eq "-") ? 1 : 0;

      if (! $query_is_gap) { 
        #query is a residue, either continutation of non-gap in query or end of gap in query
# NEW
        # query is a residue, sbjct is a residue
        if ($local_StateNotGap == $local_state) {
          # previous position was a residue in the query too, so do nothing
          ;
        }
        else {
          # previous position was a gap in query, so this is the end of a deletion
          $local_QueryEndDeletion = $local_QueryPositionIndex + 1;
          $local_OneDeletionString = "Q" . $local_QueryEndDeletion . ":" . "S" . $local_SubjectStartDeletion . "-" . ($local_residue_qlen * $local_new_gap_overhang);
          # if this is the max length deletion we've seen, update the max length deletion
          if (($local_residue_qlen * $local_new_gap_overhang) > $local_longest_gap) {
            $local_longest_gap = $local_residue_qlen * $local_new_gap_overhang;
            $local_longest_gap_str = $local_OneDeletionString;
          }
          if ($local_AllDeletionString) {
            $local_AllDeletionString = $local_OneDeletionString . ";" . $local_AllDeletionString;
          }
          else {
            $local_AllDeletionString = $local_OneDeletionString;			
          }
          $local_state = $local_StateNotGap;
          $local_new_gap_overhang = 0;
        }
        if(! $sbjct_is_gap) { 
          $local_SubjectPositionIndex += 1;
        }
        $local_QueryPositionIndex -= $local_residue_qlen;
      }
      elsif($query_is_gap && (! $sbjct_is_gap)) { 
        if ($local_StateNotGap == $local_state) {   #this is the start of an insertion in the subject
          # previous position was a residue in the query, so this is the start of a query gap (deletion)
          if ($DEBUG) {
            print "In findQueryGaps: start of a gap in the query at subject $local_SubjectPositionIndex\n";
          }
          $local_new_gap_overhang = 1;
          $local_QueryStartDeletion = $local_QueryPositionIndex+1; 
          $local_SubjectStartDeletion = $local_SubjectPositionIndex-1; 
          $local_state = $local_StateInGap;
        }
        else {
          # previous position was a gap in the query, so this is a continutation of a query gap, extend it
          $local_new_gap_overhang++;
        }
        $local_SubjectPositionIndex += 1;
      }
    }
    return($local_new_gap_overhang,$local_AllDeletionString, $local_longest_gap_str, $local_longest_gap);
  }
}

##########################################################################################
# Subroutine:  findSubjectGaps()
#
# Purpose:     Find the positions of any gaps in the subject, which represent in-frame insertions 
#              in the query. Return any completed gap starts and gap lengths and the number of gap characters
#              at the end of the line
# Arguments:
#   $query_alignment_string:     a subsequence of the query, possibly containing gaps
#   $subject_alignment_string:   a subsequence of the query, possibly containing gaps
#   $query_start:                query start position for current block we are processing
#   $query_end:                  query end  position for current block we are processing
#   $subject_start:              subject start position for current block we are processing
#   $subject_end:                subject end  position for current block we are processing
#   $old_gap_overhang:           number of gap positions overhanging from previous row
#   $program:                    from --program, "x" if blastx, "n" if blastn
#   $query_strand:               "+" if query is positive strand, "-" if minus strand
#   $subject_strand:             "+" if subject is positive strand, "-" if minus strand
#
#
#  Returns:   4 values:
#             1) the number of gap characters at the end of the row 
#             2) a string of gap-length pairs in subject or "" in >= 1 ';' delimited tokens
#                each token in format format: "Q<d1>:S<d2>+<d3>" 
#                <d1> is query position, <d2> is subject position, <d3> is length in nt
#             3) string reporting largest gap on row in same format as 2, but only 1 token
#             4) largest gap length in nt (integer)
#
# Dies:       Never
#
##########################################################################################
sub findSubjectGaps {
  my $sub_name = "findSubjectGaps()";
  my $nargs_expected = 10;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($query_alignment_string, $subject_alignment_string, $query_start, $query_end, $subject_start, $subject_end, $old_gap_overhang, $program, $query_strand, $subject_strand) = @_;

  my $local_residue_qlen = ($program eq "x") ? 3 : 1; # length of one query position in nt in query (3 if blastx, else 1 (blast))
  my $local_length = length($query_alignment_string);  #length of this alignment block
  my @local_query_array = split("",$query_alignment_string);
  my @local_subject_array = split("", $subject_alignment_string);
  my $local_i;                    #loop index
  my $local_state;                #state of a two-state automaton
  my $local_StateNotGap = 0;      #one state
  my $local_StateInGap = 1;       #the other state
  my $local_QueryPositionIndex;   #where are we in absolute query positions
  my $local_SubjectPositionIndex; #where are we in absolute subject positions
  my $local_QueryStartDeletion;   #absolute position in query starting a gap in the subject
  my $local_QueryEndDeletion;     #absolute position in query ending a gap in the subject
  my $local_SubjectStartDeletion; #absolute position in subject starting a gap in the subject
  my $local_SubjectEndDeletion;   #absolute position in subject ending a gap in the subject
  my $local_OneInsertionString;   #string to represent one gap
  my $local_AllInsertionString = ""; #string to represent all gaps        
  my $local_new_gap_overhang = $old_gap_overhang;   #number of gap characters at the end of the alignment row
  my $local_longest_gap     = 0;
  my $local_longest_gap_str = "0"; # string representing longest gap "Q<d1>:S<d2>+<d3>" <d1> is query position, <d2> is subject position, <d3> is length in nt
  
  if((! defined $query_strand) || 
     (($query_strand ne "+") && ($query_strand ne "-"))) { 
    die "ERROR read in $sub_name, query strand is not + or - (or is undefed"; 
  }
  if((! defined $subject_strand) || 
     (($subject_strand ne "+") && ($subject_strand ne "-"))) { 
    die "ERROR read in $sub_name, subject strand is not + or - (or is undefed"; 
  }

  # It is always true that $subject_start < $subject_end
  if ($subject_strand eq "-") { 
    die "ERROR in $sub_name, subject is minus strand, NOT YET IMPLEMENTED (use --splus to skip alignment parsing for such hits)";
  }
  # following block assumes $subject_strand is "+", which we just checked for above
  $local_QueryPositionIndex = $query_start;
  $local_SubjectPositionIndex = $subject_start;
  if ($old_gap_overhang) {
    $local_QueryStartDeletion = $query_start - $old_gap_overhang - 1;
    $local_SubjectStartDeletion = $subject_start - 1;
    $local_state = $local_StateInGap;
  }
  else {
    $local_state = $local_StateNotGap;	    
  }
  for ($local_i = 0; $local_i < $local_length; $local_i++) {
    my $query_is_gap = ($local_query_array[$local_i]   eq "-") ? 1 : 0;
    my $sbjct_is_gap = ($local_subject_array[$local_i] eq "-") ? 1 : 0;

    if ((! $sbjct_is_gap)) { 
      # sbjct is a residue, query is either a residue or a gap
      if ($local_StateNotGap == $local_state) {
        # previous position was a residue in the sbjct too, so do nothing
        ;
      }
      else {
        # previous position was a gap in subject, so this is the end of a deletion in sbjct (insertion in query)
        $local_QueryEndDeletion = $local_QueryPositionIndex;
        $local_SubjectEndDeletion = $local_SubjectPositionIndex;
        if ($DEBUG) {
          print "In findSubjectGaps: Setting local_QueryEndDeletion to $local_QueryEndDeletion, when local_QueryStartDeletion is $local_QueryStartDeletion\n";
        }
        if ($query_strand eq "+") { 
          $local_OneInsertionString = "Q" . $local_QueryStartDeletion . ":" . "S" . $local_SubjectStartDeletion . "+" . ($local_residue_qlen * $local_new_gap_overhang);
        }
        else { # $query_strand is "-"
          $local_OneInsertionString = "Q" . $local_QueryEndDeletion . ":" . "S" . $local_SubjectEndDeletion . "+" . ($local_residue_qlen * $local_new_gap_overhang);		    
        }
        # if this is the max length sbjct deletion/query insertion we've seen, update the max length
        if (($local_residue_qlen * $local_new_gap_overhang) > $local_longest_gap) {
          $local_longest_gap = $local_residue_qlen * $local_new_gap_overhang;
          $local_longest_gap_str = $local_OneInsertionString;
        }
        if ($local_AllInsertionString) {
          if ($query_strand eq "+") { 
            $local_AllInsertionString = $local_AllInsertionString . ";" . $local_OneInsertionString;
          }
          else { # query_strand is -
            $local_AllInsertionString = $local_OneInsertionString . ";" . $local_AllInsertionString;			
          }
        }
        else {
          $local_AllInsertionString = $local_OneInsertionString;			
        }
        $local_state = $local_StateNotGap;
        $local_new_gap_overhang = 0;
      }
      $local_SubjectPositionIndex += 1;
      if(! $query_is_gap) { 
        if ($query_strand eq "+") { 
          $local_QueryPositionIndex += $local_residue_qlen;
        }
        else { # $query_strand is -
          $local_QueryPositionIndex -= $local_residue_qlen;		
        }
      }
    }
    elsif((! $query_is_gap) && ($sbjct_is_gap)) { 
      if ($local_StateNotGap == $local_state) {   #this is the start of an insertion in the query
        # previous position was a residue in the subject, so this is the start of an subject gap (insertion)
        $local_new_gap_overhang = 1;
        $local_SubjectStartDeletion = $local_SubjectPositionIndex-1;
        if ($query_strand eq "+") { 
          $local_QueryStartDeletion = $local_QueryPositionIndex-1;
        }
        else { # query_strand is -
          $local_QueryStartDeletion = $local_QueryPositionIndex+1;			
        }
        $local_state = $local_StateInGap;
      }
      else {
        # previous position was a gap in the subject, so this is a continutation of a subject gap, extend it
        $local_new_gap_overhang++;
      }
      if ($query_strand eq "+") { 
        $local_QueryPositionIndex += $local_residue_qlen;
      }
      else { # query_strand is -
        $local_QueryPositionIndex -= $local_residue_qlen;		
      }		    
    }
  }
  if ($DEBUG) {
    print "In findSubjectGaps: about to return $local_new_gap_overhang  $local_AllInsertionString\n";
  }    
  return($local_new_gap_overhang,$local_AllInsertionString,$local_longest_gap_str, $local_longest_gap);
}

##########################################################################################
# Subroutine:  findQueryStops()
#
# Purpose:     Find the positions of any predicted stop codons in the query. 
#              Return a list of the positions.
# Arguments:
#   $query_alignment_string:     a subsequence of the query, possibly containing gaps
#   $subject_alignment_string:   a subsequence of the subject, possibly containing gaps
#   $query_start:                query start position for current block we are processing
#   $query_end:                  query end  position for current block we are processing
#   $subject_start:              subject start position for current block we are processing
#   $subject_end:                subject end  position for current block we are processing
#   $query_strand:               "+" if query is positive strand, "-" if minus strand
#
# Returns:    1) a string of query positions (nt coords) with stop codons or NULL. 
#             2) a string of subject positions (AA coords) that align to stop codons or NULL
# Dies:       Never
#
##########################################################################################
sub findQueryStops {
    my $sub_name = "findQueryStops()";
    my $nargs_expected = 7;
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

    my ($query_alignment_string, $subject_alignment_string, $query_start, $query_end, $subject_start, $subject_end, $query_strand) = @_;

    my $local_length = length($query_alignment_string);          #length of this alignment block
    my @local_query_array   = split("",$query_alignment_string);   #query line in the alignment split into an array of characters
    my @local_subject_array = split("",$subject_alignment_string); #subject line in the alignment split into an array of characters
    my $local_i;                                                 #loop index
    my $local_QueryPositionIndex;                                #where are we in absolute query positions
    my $local_SubjectPositionIndex;                              #where are we in absolute query positions
    my $local_QueryAllStopString = "";                           #string to represent all query stop codons        
    my $local_SubjectAllStopString = "";                         #string to represent all subject stop codons        

    if ($DEBUG) {
      print "In findQueryStops: Query start is $query_start\n$query_alignment_string\n$subject_alignment_string";
    }

    if((! defined $query_strand) || 
       (($query_strand ne "+") && ($query_strand ne "-"))) { 
      die "ERROR read in $sub_name, strand is not + or - (or is undefed"; 
    }
    
    my $query_idx   = $query_start;
    my $subject_idx = $subject_start;
    if ($query_strand eq "+") { 
      for ($local_i = 0; $local_i < $local_length; $local_i++) {
        if ($local_query_array[$local_i] eq "*") {
          $local_QueryAllStopString   .= sprintf("%d..%d:+;", $query_idx, $query_idx + 2); # nt coords (3 positions)
          $local_SubjectAllStopString .= sprintf("%d..%d:+;", $subject_idx, $subject_idx); # AA coords (1 position)
        }
        if ($local_query_array[$local_i] ne "-") { 
          $query_idx += 3;
        }
        if ($local_subject_array[$local_i] ne "-") { 
          $subject_idx++;
        }
      }
    }
    else { # query_strand is "-"
      for ($local_i = 0; $local_i < $local_length; $local_i++) {
        if ($local_query_array[$local_i] eq "*") {
          $local_QueryAllStopString   .= sprintf("%d..%d:-;", $query_idx, $query_idx - 2); # nt coords (3 positions)
          $local_SubjectAllStopString .= sprintf("%d..%d:+;", $subject_idx, $subject_idx); # AA coords (1 position)
        }
        if ($local_query_array[$local_i] ne "-") { 
          $query_idx -= 3;
        }
        if ($local_subject_array[$local_i] ne "-") { 
          $subject_idx++;
        }
      }
    }
    return($local_QueryAllStopString, $local_SubjectAllStopString);
}
