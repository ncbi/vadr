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
$usage .= "  -passfail    : only print counts of PF and result values\n";
$usage .= "  -outkey <s>  : save file with values for each accession for keys matching <s>\n";
$usage .= "  -outname <s> : name output file <s>\n";
$usage .= "\n";

# general options:
my $do_passfail = 0;     # set to '1' if -passfail enabled, only print counts for 'PF' and 'result' values
my $out_key     = undef; # defined if -outkey used
my $out_name    = undef; # defined if -outname used

&GetOptions("passfail"  => \$do_passfail,
            "outkey=s"  => \$out_key,
            "outname=s" => \$out_name) ||
    die "Unknown option";

if(scalar(@ARGV) != 1) { die $usage; }
my ($infile) = (@ARGV);

if((  defined $out_key) && (! defined $out_name)) { die "ERROR -outname requires -outkey also used"; }
if((! defined $out_key) && (  defined $out_name)) { die "ERROR -outkey requires -outname also used"; }

if(defined $out_name) { 
  open(OUT, ">", $out_name) || die "ERROR unable to open $out_name for writing"; 
}

my $do_seqcol = 1;       # output file is in sequences-as-columns format, currently hard-coded to be true

my $out_match      = 0;  # set to 1 if $out_key is defined and current key matches $out_key, we'll print this value to output file
my %key_ct_HH      = (); # 2D hash, first  dim key is a row heading in seqcol mode, or a column heading in seq-row mode,
                         #          second dim key is an observed value for the first key in the table
                         #          value is the number of occurences of the 2nd dim key
my @key_order_A    = (); # array of order of 1st dim keys in key_ct_HH, as they're observed in the file
my %key_exists_H   = (); # hash of 1st dim keys in key_ct_HH, used to help populate key_order_A only
my %value_order_HA = (); # hash of arrays: key: a 1st dim key in $key_ct_HH, value array of all 2nd dim keys in key_ct_HH for that 1st key, in order
my %wvalue_H       = (); # hash, key: 1st dim key from key_ct_HH, value is the max width of all 2nd dim keys in key_ct_HH for that 1st key
my $key;
my $value;
my @cur_accn_A = ();       # current accessions for this page, cleared when we see a 'result' line
my %processed_accn_H = (); # key is an accession, value is '1' if we've already processed this accession,
                           # this allows us to not multiple count the reference even though it appears many times in the output
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
      else { # normal line
        ## example: 
        # MP #14 [single exon; +]:RNA-dependent RNA polymerase NS5:start1                         7565                [7193]
        # determine number of space-delimited tokens, final $naccn will be values
        my @el_A = split(/\s+/, $line);
        my $nel = scalar(@el_A);

        my $is_accession_line = ($line =~ m/^accession/) ? 1 : 0;
        my $is_result_line    = ($line =~ m/^result/)    ? 1 : 0;
        
        # determine key, this is the value in the first column, "MP #14 [single exon; +]:RNA-dependent RNA polymerase NS5:start1" in the example above
        my $key = $line;
        my $ntok_per_accn = 1; # number of tokens we expect per accession, usually 1
        if($line =~ m/^overlaps/) { $ntok_per_accn = 3; } # 3 tokens per accession for overlaps columns
        if($line =~ m/^result/)   { $ntok_per_accn = 2; } # 2 tokens per accession for result columns
        for(my $i = 0; $i < ($naccn * $ntok_per_accn); $i++) { 
          $key =~ s/\S+\s*$//; # remove final token, ($naccn * $ntok_per_accn) times to get key
        }
        $key =~ s/\s+$//;
        $out_match = ((defined $out_key) && $key =~ m/$out_key/) ? 1 : 0;

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
        my $column = 0;
        for(my $i = $first_accn; $i < ($first_accn + ($naccn * $ntok_per_accn)); $i++) { 
          my $value = "";
          my $j;
          for($j = $i; $j < ($i + ($ntok_per_accn)); $j++) { # we have to add all the tokens for multi-token cases ($ntok_per_accn > 1)
            if($j > $i) { $value .= " "; }
            $value .= $el_A[$j];
          }
          $i = $j-1; 

          if($is_accession_line) { 
            push(@cur_accn_A, $value); 
            # printf("pushing $value to cur_accn_A size %d\n", scalar(@cur_accn_A));
          }

          my $cur_accn = $cur_accn_A[$column];
          
          if(! exists $processed_accn_H{$cur_accn}) { 
             # only store data for accessions we haven't processed already, this prevents us from multiple counting of reference data
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
            if($out_match) { 
              print OUT $cur_accn . " " . $value . "\n"; 
            }
          }
          $column++;
        } # end of 'for(my $i = $first_accn; $i < ($first_accn + ($naccn * $ntok_per_accn)); $i++) {'
        if($is_result_line) { # we're done with this page, store all accessions as having been processed
          foreach my $accn (@cur_accn_A) { 
            $processed_accn_H{$accn} = 1; 
          }
          @cur_accn_A = (); # move onto next page with blank cur_accn_A
        }
      }
    }
  }
}

# output summary
foreach $key (@key_order_A) { 
  if((! $do_passfail) || ($key =~ m/PF$/) || ($key =~ m/overlaps/) || ($key eq "result")) { 
    printf("$key\n");
    foreach $value (@{$value_order_HA{$key}}) { 
      printf("%-*s %d\n", $wvalue_H{$key}, $value, $key_ct_HH{$key}{$value});
    }
    printf("\n");
  }
}

if(defined $out_name) { 
  close(OUT); 
  printf("# Saved accessions and values for headings matching $out_key to file $out_name.\n");
}

