#!/usr/bin/env perl
# EPN, Thu Oct  1 14:41:00 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

my $usage = "\nsummarize_dnaorg_annotate_genomes.pl\n";
$usage .= "\t<output file created by dnaorg_annotate_genomes.pl>\n";
$usage .= "\n"; 
$usage .= " This script summarizes an output file created by dnaorg_annotate_genomes.pl\n";
$usage .= " by listing counts of each value observed for each category.\n";
$usage .= " BASIC OPTIONS:\n";
$usage .= "  -passfail : only print counts of PF and result values\n";
$usage .= "\n";

# general options:
my $do_passfail = 0; # set to '1' if -passfail enabled, only print counts for 'PF' and 'result' values

&GetOptions("passfail" => \$do_passfail) ||
    die "Unknown option";

if(scalar(@ARGV) != 1) { die $usage; }
my ($infile) = (@ARGV);

my $do_seqcol = 1;

my %key_ct_HH      = (); # 2D hash, first  dim key is a row heading in seqcol mode, or a column heading in seq-row mode,
                         #          second dim key is an observed value for the first key in the table
                         #          value is the number of occurences of the 2nd dim key
my @key_order_A    = (); # array of order of 1st dim keys in key_ct_HH, as they're observed in the file
my %key_exists_H   = (); # hash of 1st dim keys in key_ct_HH, used to help populate key_order_A only
my %value_order_HA = (); # hash of arrays: key: a 1st dim key in $key_ct_HH, value array of all 2nd dim keys in key_ct_HH for that 1st key, in order
my %wvalue_H       = (); # hash, key: 1st dim key from key_ct_HH, value is the max width of all 2nd dim keys in key_ct_HH for that 1st key
my $key;
my $value;

open(IN, $infile) || die "ERROR unable to open $infile for writing";

my %skip_value_H = ();
$skip_value_H{"accession"} = 1;
$skip_value_H{"totlen"} = 1;

if($do_seqcol) { 
  my $naccn = 0; # number of elements we expect in each line, final 
  while(my $line = <IN>) { 
    if($line !~ m/^\#/) { 
      chomp $line;
      ## example:
      # idx                                                                                        1                     2
      if($line =~ m/^idx/) { 
        my @el_A = split(/\s+/, $line);
        $naccn = scalar(@el_A) - 1; # number of accessions we expect in each subsequent line
      }
      elsif($line =~ m/^accession/) { 
        ; # do nothing
      }
      else { # normal line
        ## example: 
        # MP #14 [single exon; +]:RNA-dependent RNA polymerase NS5:start1                         7565                [7193]
        # determine number of space-delimited tokens, final $naccn will be values
        my @el_A = split(/\s+/, $line);
        my $nel = scalar(@el_A);
        
        # determine key, this is the value in the first column, "MP #14 [single exon; +]:RNA-dependent RNA polymerase NS5:start1" in the example above
        my $key = $line;
        my $ntok_per_accn = 1; # number of tokens we expect per accession, usually 1
        if($line =~ m/^overlaps/) { $ntok_per_accn = 3; } # 3 tokens per accession for overlaps columns
        if($line =~ m/^result/)   { $ntok_per_accn = 2; } # 2 tokens per accession for result columns
        for(my $i = 0; $i < ($naccn * $ntok_per_accn); $i++) { 
          $key =~ s/\S+\s*$//; # remove final token, ($naccn * $ntok_per_accn) times to get key
        }
        $key =~ s/\s+$//;
        # initialize count hash for this key if it doesn't already exist
        if(! exists $key_ct_HH{$key}) {
          %{$key_ct_HH{$key}} = ();
          @{$value_order_HA{$key}} = ();
          $wvalue_H{$key} = 0;
        }
        # keep track of what order this key is, in the @key_order_A array
        if(! exists $key_exists_H{$key}) { 
          push(@key_order_A, $key); 
          $key_exists_H{$key} = 1;
        }
        
        # determine the first token that is an accession value
        my $first_accn = $nel - ($naccn * $ntok_per_accn);
        for(my $i = $first_accn; $i < ($first_accn + ($naccn * $ntok_per_accn)); $i++) { 
          my $value = "";
          my $j;
          for($j = $i; $j < ($i + ($ntok_per_accn)); $j++) { # we have to add all the tokens for multi-token cases ($ntok_per_accn > 1)
            if($j > $i) { $value .= " "; }
            $value .= $el_A[$j];
          }
          $i = $j-1; 
          
          if(! exists $key_ct_HH{$key}{$value}) { 
            $key_ct_HH{$key}{$value} = 1;
            push(@{$value_order_HA{$key}}, $value);
            if(length($value) > $wvalue_H{$key}) { 
              $wvalue_H{$key} = length($value); 
            }
          }
          else { 
            $key_ct_HH{$key}{$value}++;
          }
        }
      }
    }
  }
}

# output summary
foreach $key (@key_order_A) { 
  if((! $do_passfail) || ($key =~ m/PF$/) || ($key eq "result")) { 
    printf("$key\n");
    foreach $value (@{$value_order_HA{$key}}) { 
      printf("%-*s %d\n", $wvalue_H{$key}, $value, $key_ct_HH{$key}{$value});
    }
    printf("\n");
  }
}

