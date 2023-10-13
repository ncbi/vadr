#!/usr/bin/env perl
$usage = "perl esl-alipid-per-seq-stats.pl <alipid output file>";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t-i <f>: ignore the sequences listed in file f\n\n";

use Getopt::Long;

my $ignore_file = undef;
&GetOptions( "i=s" => \$ignore_file);

if(scalar(@ARGV) != 1) { die $usage; };

my ($alipid_file) = (@ARGV);

my %pid_H   = ();
my %denom_H = ();
my @seq_A = ();
my %seq_H = ();
my %max_H = ();
my %min_H = ();

my %ignore_H = ();
if(defined $ignore_file) { 
  open(IGNORE, $ignore_file) || die "ERROR unable to open $ignore_file"; 
  while($line = <IGNORE>) { 
    chomp $line;
    $ignore_H{$line} = 1;
  }
  close(IGNORE);
}

my $width = length("minpidseq");
open(ALIPID, $alipid_file) || die "ERROR unable to open $alipid_file";

while($line = <ALIPID>) { 
  if($line !~ m/^\#/) { 
    ## seqname1 seqname2 %id nid denomid %match nmatch denommatch
    #AY743607.1 AY387239.1  94.61    527    557  90.15    540    599
    my($seq1, $seq2, $id1, $match1, $len1, $id2, $match2, $match2) = split(/\s+/, $line);
    $pid_H{$seq1} += $id1;
    $pid_H{$seq2} += $id1;
    $denom_H{$seq1}++;
    $denom_H{$seq2}++;

    if((! exists $ignore_H{$seq1}) && 
       (! exists $ignore_H{$seq2})) { 

      if(! exists $seq_H{$seq1}) { 
        push(@seq_A, $seq1);
        $seq_H{$seq1} = 1;
        $min_H{$seq1} = 100.;
        $max_H{$seq1} = 0.;
        $argmin_H{$seq1} = undef;
        $argmax_H{$seq1} = undef;
        if(length($seq1) > $width) { $width = length($seq1); }
      }
      if(! exists $seq_H{$seq2}) { 
        push(@seq_A, $seq2);
        $seq_H{$seq2} = 1;
        $min_H{$seq2} = 100.;
        $max_H{$seq2} = 0.;
        $argmin_H{$seq2} = undef;
        $argmax_H{$seq2} = undef;
        if(length($seq2) > $width) { $width = length($seq2); }
      }
      if($id1 > $max_H{$seq1}) { 
        $max_H{$seq1}    = $id1;
        $argmax_H{$seq1} = $seq2; 
      }
      if($id1 < $min_H{$seq1}) { 
        $min_H{$seq1}    = $id1;
        $argmin_H{$seq1} = $seq2; 
      }
      if($id1 > $max_H{$seq2}) { 
        $max_H{$seq2}    = $id1;
        $argmax_H{$seq2} = $seq1; 
      }
      if($id1 < $min_H{$seq2}) { 
        $min_H{$seq2}    = $id1;
        $argmin_H{$seq2} = $seq1; 
      }
    }
  }
}

printf("%-*s  %6s  %-*s  %6s  %-*s  %6s\n", $width, "#seq", "avgpid", $width, "minpidseq", "minpid", $width, "maxpidseq", "maxpid");
foreach $seq (@seq_A) { 
  if(! defined $argmax_H{$seq}) { 
    die "ERROR did not read any lines with $seq\n";
  }
  printf("%-*s  %.3f  %-*s  %.3f  %-*s  %.3f\n", 
         $width, $seq, $pid_H{$seq}/$denom_H{$seq}, 
         $width, $argmin_H{$seq}, $min_H{$seq},
         $width, $argmax_H{$seq}, $max_H{$seq});
}


