#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Easel::MSA;
use Bio::Easel::SqFile;
require "sqp_seq.pm";
require "sqp_utils.pm";

my $usage;
$usage  = "fasta-nucleotide-trim-ambigs-genbank.pl\n\n";
$usage .= "Usage:\n\n";
$usage .= "perl fasta-nucleotide-trim-ambigs-genbank.pl [OPTIONS] <fasta file";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--ten <n>    : max number of ambiguous nucleotides allowed in first/final 10 [5]\n";
$usage .= "\t\t--fifty <n>  : max number of ambiguous nucleotides allowed in first/final 10 [15]\n";
$usage .= "\t\t--sfx <s>    : suffix to add to each sequence name is <s> [default: do not change names]\n";
$usage .= "\t\t--skipzero   : silently skip seqs of length 0 after trimming [default: exit in error if seq becomes length 0]\n";

# set defaults
my $ten_max_ambig   = 5;
my $fifty_max_ambig = 15;
my $sfx             = undef;
my $do_skipzero     = 0;

if(scalar(@ARGV) != 1) { 
  die $usage;

&GetOptions( "ten=s"    => \$ten_max_ambig,
             "fifty=s"  => \$fifty_max_ambig,
             "sfx=s"    => \$sfx, 
             "skipzero" => \$do_skipzero);
}

if(scalar(@ARGV) != 1) { die $usage; }
my ($fasta_file) = @ARGV;

if(($ten_max_ambig < 0) || ($ten_max_ambig >= 10)) { 
  die "ERROR with --ten <n>, <n> must in range [0..9]";
}
if(($fifty_max_ambig < 0) || ($fifty_max_ambig >= 49)) { 
  die "ERROR with --ten <n>, <n> must in range [0..49]";
}

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

  my $orig_sqstring = $sqstring;

  #print $out_header . "\n";
  #printf("original length: %d\n", length($sqstring));

  # trim the sequence at 5' end, then reverse it and trim at 3' end (this makes it so we can reuse same code for both ends)
  # then reverse it back and output it (if there's any sequence left after trimming)
  $sqstring = trim_5p_end_using_three_rules($sqstring, $ten_max_ambig, $fifty_max_ambig);
  #printf("5' trimmed length: %d\n", length($sqstring));
  if($sqstring ne "") { 
    $sqstring = reverse($sqstring);
    $sqstring = trim_5p_end_using_three_rules($sqstring, $ten_max_ambig, $fifty_max_ambig);
    #printf("3' trimmed length: %d\n", length($sqstring));
    if($sqstring ne "") { # reverse it back to original forward direction
      $sqstring = reverse($sqstring);
    }
  }
  if($sqstring eq "") { 
    # zero length sequence
    if($do_skipzero) { 
      ; # --skipzero used, do nothing
    }
    else { 
      # --skipzero not used, die
      die "ERROR sequence with the header line below is trimmed to length 0.\nTo silently skip these seqs (and not output them) use --skipzero.\n$header\n";
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
# Subroutine:  trim_5p_end_using_three_rules()
# Incept:      EPN, Wed Apr  7 07:22:43 2021
#
# 
# Purpose:   Given a sequence string <$sqstring>, remove ambiguous
#            nts from the 5' end the same way that GenBank processing 
#            does, using 3 rules:
#            while any rule results in trimming 1 or more nt:
#            - rule 1: remove 5' terminal non-ACGTU nts
#            - rule 2: remove the 10 5'-most nt if > $ten_max_ambig 
#                      ambiguous nt exist in 10 5'-most nt [5 by default]
#            - rule 3: remove the 50 5'-most nt if > $fifty_max_ambig 
#                      ambiguous nt exist in 50 5'-most nt [15 by default]
#           
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
  while(($keep_going) && ($sqstring ne "")) { 
    $keep_going = 0; # set to 1 if we use any rules below
    # rule 1: remove 5' terminal non-ACGTU nts
    if($sqstring =~ m/[^ACGTUacgtu]/) { 
      $sqstring =~ s/[^ACGTUacgtu]+//;
      $keep_going = 1;
    }
    # rule 2: remove the 10 5'-most nt if > $ten_max_ambig ambiguous nt exist in 10 5'-most nt
    if(! $keep_going) { # only enforce next rule if we didn't use any previous rules
      my $next_10 = substr($sqstring, 0, 10);
      my $nambig = () = $next_10 =~ /[^ACGTUacgtu]/g;
      if($nambig > $ten_max_ambig) { 
        # trim first 10 away
        $sqstring = substr($sqstring, 10);
        $keep_going = 1;
      }
    }
    if(! $keep_going) { # only enforce next rule if we didn't use any previous rules
      my $next_50 = substr($sqstring, 0, 50);
      my $nambig = () = $next_50 =~ /[^ACGTUacgtu]/g;
      if($nambig > $fifty_max_ambig) { 
        # trim first 50 away
        $sqstring = substr($sqstring, 50);
        $keep_going = 1;
      }
    }
  }
  return $sqstring;
}
