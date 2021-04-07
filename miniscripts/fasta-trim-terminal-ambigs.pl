#!/usr/bin/env perl

### comments that start with ### pertain to code that could be used
### to try and reproduce how submission portal trims seqs for sars-cov-2
### - trim in 10nt chunks if > 5 ambigs
### - trim in 50nt chunks if > 15 ambigs
### - trim terminal nt
### - min length 50
### - max length 30000
### - max fraction Ns is 0.5

use strict;
use warnings;
use Bio::Easel::MSA;
use Bio::Easel::SqFile;
require "sqp_seq.pm";
require "sqp_utils.pm";

my $usage;
###$usage  = "fasta-sarscov2-trim-ambigs-genbank.pl\n\n";
$usage  = "fasta-trim-terminal-ambigs.pl\n\n";
$usage .= "Usage:\n\n";
$usage .= "perl fasta-trim-terminal-ambigs.pl [OPTIONS] <fasta file";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--minlen <n> : min allowed sequence length after trimming [df: 1]\n";
$usage .= "\t\t--maxlen <n> : max allowed sequence length after trimming [df: 1Gb]\n";
###$usage .= "\t\t--maxfrac <f>: max allowed fraction of sequence that can be Ns after trimming [0.5]\n";
###$usage .= "\t\t--ten <n>    : max number of ambiguous nucleotides allowed in first/final 10 [5]\n";
###$usage .= "\t\t--fifty <n>  : max number of ambiguous nucleotides allowed in first/final 10 [15]\n";
$usage .= "\t\t--sfx <s>    : suffix to add to each sequence name is <s> [default: do not change names]\n";
$usage .= "\t\t--strict     : die (instead of skipping a seq) if any seq is not within minlen..maxlen range after trimming\n";

# set defaults
###my $minlen          = 50;
###my $maxlen          = 30000;
my $minlen          = 1;
my $maxlen          = 1000000000; # 1Gb

### my $maxfrac_Ns      = 0.5;
### my $ten_max_ambig   = 5;
### my $fifty_max_ambig = 15;
my $sfx             = undef;
my $do_strict       = 0;

if(scalar(@ARGV) != 1) { 
  die $usage;

&GetOptions( "minlen=s" => \$minlen,
             "maxlen=s" => \$maxlen,
###             "maxfrac=s"=> \$maxfrac_Ns,
###             "ten=s"    => \$ten_max_ambig,
###             "fifty=s"  => \$fifty_max_ambig,
             "sfx=s"    => \$sfx, 
             "strict"   => \$do_strict);
}

if(scalar(@ARGV) != 1) { die $usage; }
my ($fasta_file) = @ARGV;

# enforce ranges that make sense for options
if($minlen < 0) { 
  die "ERROR with --minlen <n>, <n> must be >= 0";
}
if($maxlen < 0) { 
  die "ERROR with --maxlen <n>, <n> must be >= 0";
}
###if(($ten_max_ambig < 0) || ($ten_max_ambig >= 10)) { 
###  die "ERROR with --ten <n>, <n> must in range [0..9]";
###}
###if(($fifty_max_ambig < 0) || ($fifty_max_ambig >= 49)) { 
###  die "ERROR with --ten <n>, <n> must in range [0..49]";
###}
if($minlen > $maxlen) { 
  die "ERROR with --minlen <n1> and --maxlen <n2>, <n1> must not be greater than <n2>";
}
###if($ten_max_ambig > $fifty_max_ambig) { 
###  die "ERROR with --ten <n1> and --fifty <n2>, <n1> must not be greater than <n2>";
###}
###if(($maxfrac_Ns < -0.00001) || ($maxfrac_Ns > 1.00001)) { 
###  die "ERROR with --maxfrac <f>, <f> must be between 0 and 1";
###}

if(! -s $fasta_file) { 
  die "ERROR fasta file $fasta_file does not exist or is empty";
}

# determine number of sequences in file
my $nseq = `grep ^\\> $fasta_file | wc -l`;
chomp $nseq;
if($nseq !~ /^\d+$/) { 
  die "ERROR could not determine number of sequences in file using grep and wc";
}

# open Bio::Easel SqFile object
my $sqfile  = Bio::Easel::SqFile->new({ fileLocation => $fasta_file }); # the sequence file object

my $fasta_seq = undef;
my @fasta_seq_A = ();
my ($header, $sqstring, $out_header) = (undef, undef, undef);

# for each sequence, remove Ns and output to STDOUT
for(my $i = 0; $i < $nseq; $i++) { 
  $fasta_seq = $sqfile->fetch_consecutive_seqs(1, "", -1); # -1: unlimited line length
  @fasta_seq_A = split(/\n/, $fasta_seq);
  if(scalar(@fasta_seq_A) != 2) { 
    die "ERROR fetching sequence, could not split into header and sequence lines:\n$fasta_seq\n"; 
  }
  ($header, $sqstring) = (@fasta_seq_A);

  # add sfx to sequence name, if --sfx used
  if(defined $sfx) { 
    if($header =~ /^\>(\S+)(\s*.*)$/) { 
      $out_header = ">" . $1 . $sfx . $2;
    }
    else { 
      die "ERROR unable to parse header line:\n$header\n";
    }      
  }        
  else { # --sfx not used
    $out_header = $header;
  }

#####################################################
###  
###    # trim the sequence at 5' end, then reverse it and trim at 3' end (this makes it so we can reuse same code for both ends)
###    # then reverse it back and output it (if there's any sequence left after trimming)
###    $sqstring = trim_5p_end_using_three_rules($sqstring, $ten_max_ambig, $fifty_max_ambig);
###    #printf("5' trimmed length: %d\n", length($sqstring));
###    if($sqstring ne "") { 
###      $sqstring = reverse($sqstring);
###      $sqstring = trim_5p_end_using_three_rules($sqstring, $ten_max_ambig, $fifty_max_ambig);
###      #printf("3' trimmed length: %d\n", length($sqstring));
###      if($sqstring ne "") { # reverse it back to original forward direction
###        $sqstring = reverse($sqstring);
###      }
###    }
######################################################

  # the following two lines should be removed if above block is commented out and 
  # trim_5p_end_using_three_rules() gets called 
  $sqstring =~ s/^[^ACGTUacgtu]+//; # remove any 5'-terminal ambiguous nts
  $sqstring =~ s/[^ACGTUacgtu]+$//; # remove any 3'-terminal ambiguous nts

  my $seqlen = length($sqstring);
###    my $n_N    = () = $sqstring =~ /[Nn]/g; # count Ns
###    my $frac_N = $n_N / $seqlen;

  if(($seqlen < $minlen) ||     # too short (after trimming)
     ($seqlen > $maxlen) {      # too long  (after trimming)
###       ($frac_N > $maxfrac_Ns)) { # too many Ns (after trimming)
    # do not output sequence
    if(! $do_strict) { 
      ; # --strict not used, it's ok to skip seqs, do nothing
    }
    else { 
      # --strict used, it's not ok to skip seqs, die
      if($seqlen < $minlen) { 
        die "ERROR sequence with the header line below is too short (after trimming) length is $seqlen minlen is $minlen.\n$header\n";
      }
      if($seqlen > $maxlen) { 
        die "ERROR sequence with the header line below is too short (after trimming) length is $seqlen minlen is $minlen.\n$header\n";
      }
      else { 
        die "ERROR sequence with the header line below has too many Ns (after trimming) fraction Ns is $frac_N, max allowed is $maxfrac_Ns\n$header\n";
      }
    }
  }
  else { 
    # $sqstring is not empty, output it
    print $out_header . "\n";
    print seq_SqstringAddNewlines($sqstring, 60);
  }
}

exit 0;

#################################################################
# THIS SUBROUTINE IS EXPERIMENTAL AND NOT USED
# COULD BE USED TO TRY AND REPRODUCE WHAT 
# SUBMISSION PORTAL DOES TO SARS-COV-2 SEQS
# 
# Subroutine:  trim_5p_end_using_three_rules()
# Incept:      EPN, Wed Apr  7 07:22:43 2021
#
# 
# Purpose:   Given a sequence string <$sqstring>, remove ambiguous
#            nts from the 5' end the same way that GenBank processing 
#            does, using 2 rules:
#            while any rule results in trimming:
#            - rule 1: remove the 10 5'-most nt if > $ten_max_ambig 
#                      ambiguous nt exist in 10 5'-most nt [5 by default]
#            - rule 2: remove the 50 5'-most nt if > $fifty_max_ambig 
#                      ambiguous nt exist in 50 5'-most nt [15 by default]
#            Then after that trimming is finished, trim terminal ambiguous
#            nucleotides before returning.
#
# Arguments:
#   $sqstring:        the nucleotide sequence string
#   $ten_max_ambig:   max number of ambiguous nucleotides allowed in first 10
#   $fifty_max_ambig: max number of ambiguous nucleotides allowed in first 50
#
# Returns:    the input <$sqstring> trimmed as above, or "" if none is left after trimming
#
#################################################################
sub trim_5p_end_using_three_rules { 
  my $sub_name = "trim_5p_end_using_three_rules";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($sqstring, $ten_max_ambig, $fifty_max_ambig) = @_;

  my $keep_going = 1;
  my ($next_10, $next_50, $nambig, $trim_offset); 
  #my $ntrimmed = 0;
  while(($keep_going) && ($sqstring ne "")) { 
    $keep_going = 0; # set to 1 if we use any rules below
    # rule 1: remove the 10 5'-most nt if > $ten_max_ambig ambiguous nt exist in 10 5'-most nt
    $next_10 = substr($sqstring, 0, 10);
    $nambig = () = $next_10 =~ /[^ACGTUacgtu]/g;
    if($nambig > $ten_max_ambig) { 
      # determine how many to trim, stop trimming at final ambig char in region
      # remove all non-ambiguous nt from the end to determine length to trim
      $next_10 =~ s/[ACGTUacgtu]+$//;
      $trim_offset = length($next_10);

      #$ntrimmed += $trim_offset;
      #printf("nambig: $nambig, trimming $trim_offset: %s (ntrimmed: $ntrimmed)\n", substr($sqstring, 0, $trim_offset));

      $sqstring = substr($sqstring, $trim_offset);
      $keep_going = 1;
    }
    # rule 2: remove the 50 5'-most nt if > $fifty_max_ambig ambiguous nt exist in 50 5'-most nt
    if(! $keep_going) { # only enforce this rule if we didn't use rule 1
      $next_50 = substr($sqstring, 0, 50);
      $nambig = () = $next_50 =~ /[^ACGTUacgtu]/g;
      if($nambig > $fifty_max_ambig) { 
        # determine how many to trim, stop trimming at final ambig char in region
        # remove all non-ambiguous nt from the end to determine length to trim
        $next_50 =~ s/[ACGTUacgtu]+$//;
        $trim_offset = length($next_50);
        $sqstring = substr($sqstring, $trim_offset);
        $keep_going = 1;
      }
    }
  }

  # finally, trim any leading (terminal) ambiguities
  $sqstring =~ s/^[^ACGTUacgtu]+//;

  return $sqstring;
}

