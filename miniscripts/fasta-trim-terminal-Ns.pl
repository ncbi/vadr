#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Easel::MSA;
use Bio::Easel::SqFile;
require "sqp_seq.pm";
require "sqp_utils.pm";

my $usage = "perl fasta-trim-terminal-Ns.pl\n\t<fasta file>\n\t<suffix to append to sequence names, 'NONE' to leave sequence names unchanged>\n\n";

if(scalar(@ARGV) != 2) { 
  die $usage;
}

my ($fasta_file, $sfx) = (@ARGV);

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

  # add sfx to sequence name unless NONE used
  if($sfx ne "NONE") { 
    if($header =~ /^\>(\S+)(\s*.*)$/) { 
      $out_header = ">" . $1 . $sfx . $2;
    }
    else { 
      die "ERROR unable to parse header line:\n$header\n";
    }      
  }        
  else { # sfx set as "NONE" on command-line
    $out_header = $header;
  }
  $sqstring =~ s/^[Nn]+//; # remove leading sequences
  $sqstring =~ s/[Nn]+$//; # remove trailing sequences
  print $out_header . "\n";
  print seq_SqstringAddNewlines($sqstring, 60);
}

exit 0;

