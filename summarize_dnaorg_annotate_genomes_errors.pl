#!/usr/bin/env perl
# EPN, Mon Dec 21 09:21:19 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

my $usage = "\ncat <.errors output file(s) created by dnaorg_annotate_genomes.pl> | summarize_dnaorg_annotate_genomes_errors.pl\n";
$usage .= "\n"; 
$usage .= " This script summarizes an errors output file created by dnaorg_annotate_genomes.pl\n";
$usage .= " by listing counts of each error.\n";
$usage .= " BASIC OPTIONS:\n";
$usage .= "  -n <n>    : in normal output mode, also report fraction of sequences with each error setting total as <n>\n";
$usage .= "  -pairwise : output pairwise co-occurence matrix, not standard output\n";
$usage .= "  -pwlist   : with -pairwise, list pairwise counts instead of printing a table\n";
$usage .= "\n";

# general options:
my $do_fract = 0;     # set to '1' if -n enabled
my $n        = undef; # defined if -n used.
my $do_pairwise = 0;  # set to '1' if -pw enabled
my $do_pwlist   = 0;  # set to '1' if -pwlist enabled
&GetOptions("n=s"      => \$n, 
            "pairwise" => \$do_pairwise,
            "pwlist"   => \$do_pwlist) || 
    die "Unknown option";

if(defined $n) { 
  $do_fract = 1;
}

# check for required option combinations
if(($do_pwlist) && (! $do_pairwise)) { die "ERROR -pwlist requires -pairwise also used"; }

# check for incompatible option combinations
if($do_pairwise && $do_fract) { 
  die "ERROR -pairwise and -n cannot be used in combination";
}

my %tot_code_ct_H    = (); # key: error code, e.g. trc, value total number of times code occurs
my %has_code_ct_H    = (); # key: error code, e.g. trc, value total number of accessions that have this error code >= 1 time
my %accn_code_HH = (); # key 1D: accession, 2D: error code, value 1 if accession has this error code
my %pairwise_ct_HH = (); # key 1D: error code 1, 2D: error code 2, value, number of times error code 1 and error code 2 occur >= 1 time for the same accession
my @code_A = ("ori", "nop", "nm3", "bd5", "bd3", "olp", "str", "stp", "ajb", "aja", "trc", "ext", "ntr", "nst", "aji", "int", "inp");
my $any_code_ct = 0;  # number of accession with >= 1 error code reported >= 1 times
# initialize
foreach my $code (@code_A) { 
  $tot_code_ct_H{$code} = 0;
  $has_code_ct_H{$code} = 0;
  foreach my $code2 (@code_A) { 
    $pairwise_ct_HH{$code}{$code2} = 0;
  }
}

while(my $line = <>) { 
  if($line !~ m/^\#/) { 
#     #accn       idx  desc   code  error-message
#    FJ882131      3  CDS#3   trc  in-frame stop codon exists 5' of stop position predicted by homology to reference [homology search predicted 1801..1360 exon 2 of 2 revised to 1801..1376 (stop shifted 16 nt)]
    chomp $line;
    if($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\w\w\w)\s+(.+)$/) { 
      my ($accn, $idx, $desc, $code, $errmsg) = ($1, $2, $3, $4, $5);
      $tot_code_ct_H{$code}++;
      if(! exists $accn_code_HH{$accn}) { 
        $any_code_ct++;
        %{$accn_code_HH{$accn}} = ();
      }
      if(! exists $accn_code_HH{$accn}{$code}) { 
        $has_code_ct_H{$code}++;
        foreach my $code2 (keys (%{$accn_code_HH{$accn}})) { 
          if($code2 ne $code) { 
            $pairwise_ct_HH{$code}{$code2}++;
            $pairwise_ct_HH{$code2}{$code}++;
          }
        }
        $accn_code_HH{$accn}{$code} = 1;
      }
    }
    else { 
      die "ERROR unable to parse line: $line\n";
    }
  }
}

# print standard counts table
if(! $do_pairwise) { 
  printf("# Table below includes counts of error codes\n");
  printf("#\n");
  printf("# Explanation of columns:\n");
  printf("# \"code\"%s: the error code\n", ($do_fract) ? "       " : " ");
  printf("# \"\#tot\"%s: total number of occurences of code , possibly > 1 for some accessions\n", ($do_fract) ? "       " : " ");
  printf("# \"\#accn\"%s: number of accessions with at least 1 occurence of code\n", ($do_fract) ? "      " : "");
  if($do_fract) { 
    printf("# \"fraction...\": fraction of $n accessions with at least 1 occurence of code\n");
  }    
  printf("#\n");
  printf("# Explanation of final two rows beginning with \"total\" and \"any\":\n");
  printf("#   \"total\":\"#tot\"   column is total number of error codes reported\n");
  printf("#   \"any\":\"#tot\"     column is number of accessions with >= 1 error code\n");
  if($do_fract) { 
    printf("#   \"any\":\"fraction\" column is fraction of accessions with >= 1 error code");
  }
  print("\n");

  print("#\n");

  printf("#code  \#tot  \#accn");
  if($do_fract) { 
    printf("  fraction-of-all-$n-accn");
  }
  print("\n");

  # print divider row
  printf("#---- -----  -----");
  if($do_fract) { 
    printf("  ------");
  }
  print("\n");

  my $sum_tot_code = 0;
  my $sum_has_code = 0;
  my $sum_fract    = 0;
  foreach my $code (@code_A) { 
    printf("$code   %5d  %5d", $tot_code_ct_H{$code}, $has_code_ct_H{$code});
    $sum_tot_code += $tot_code_ct_H{$code};
    $sum_has_code += $has_code_ct_H{$code};
    if($do_fract) { 
      printf("  %6.4f", $has_code_ct_H{$code} / $n); 
      $sum_fract += $has_code_ct_H{$code} / $n; 
    }
    print("\n");
  }

  # print divider row
  printf("#---- -----  -----");
  if($do_fract) { 
    printf("  ------");
  }
  print("\n");

  # print any row
  printf("total %5d  %5s", $sum_tot_code, "-");
  if($do_fract) { 
    printf("  %6s", "-");
  }
  print("\n");

  printf("any   %5d  %5s", $any_code_ct, "-");
  if($do_fract) { 
    printf("  %6.4f", $any_code_ct / $n); 
  }
  print("\n");
}

# print pairwise co-occurence table:
if($do_pairwise) { 
  if($do_pwlist) { 
    printf("# Error code pairwise combinations that co-occur for the same accession are listed below.\n");
    printf("#\n");
    printf("# First  token is first  error code.\n");
    printf("# Second token is second error code.\n");
    printf("# Third  token is number of accessions for which the two error codes are reported >= 1 time each.\n");
    printf("#\n");
  }
  else { 
    printf("# Table below includes counts of co-occurences of all pairs of codes.\n");
    printf("#\n");
    printf("# Each number indicates the number of accessions for which the two error codes\n");
    printf("# are reported >= 1 time.\n");
    printf("# A \"-\" indicates a '0' count.\n");
    printf("#\n");
  }
  my $ncodes = scalar(@code_A);
  if(! $do_pwlist) { printf("   "); }
  for(my $c1 = 0; $c1 < $ncodes; $c1++) { 
    my $code1 = $code_A[$c1];
    if(! $do_pwlist) { printf("  %5s", $code1); }
  }
  if(! $do_pwlist) { printf("\n"); }
         
  for(my $c1 = 0; $c1 < $ncodes; $c1++) { 
    my $code1 = $code_A[$c1];
    if(! $do_pwlist) { 
      printf("$code1");
    }
    for(my $c2 = 0; $c2 < $ncodes; $c2++) { 
      my $code2 = $code_A[$c2];
      if(! $do_pwlist) { 
        if($c1 == $c2) { printf("  %5s", "-"); }
        else           { printf("  %5s",  $pairwise_ct_HH{$code1}{$code2} == 0 ? "-" : sprintf("%5d", $pairwise_ct_HH{$code1}{$code2})); }
      }
      if($do_pwlist) { 
        if(($c1 < $c2) && ($pairwise_ct_HH{$code1}{$code2} > 0)) { 
          printf("$code1 $code2 $pairwise_ct_HH{$code1}{$code2}\n");
        }
      }
    }
    if(! $do_pwlist) { 
      print "\n";
    }
  }
}
