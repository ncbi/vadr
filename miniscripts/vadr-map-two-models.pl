#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Easel::MSA;
use Bio::Easel::SqFile;
require "sqp_seq.pm";
require "sqp_utils.pm";
require "sqp_ofile.pm";
require "vadr.pm";

my $usage;
$usage  = "vadr-map-two-models.pl\n\n";
$usage .= "Usage:\n\n";
$usage .= "perl vadr-map-two-models.pl [OPTIONS]\n\t<mmap file path>\n\t<model to map positions from>\n\t<model to map positions to>\n\n";

if(scalar(@ARGV) != 3) { die $usage; }
my ($mmap_file, $in_mfrom, $in_mto) = @ARGV;

if(! -s $mmap_file) { 
  die "ERROR model map file $mmap_file does not exist";
}

# read in map file and store arrays of positions from each model to $in_model_to
my %tofrom_HA = ();
open(MMAP, $mmap_file) || die "ERROR unable to open $mmap_file for reading";
my ($mfrom, $mfrom_len, $mto, $mto_len, $cigar) = (undef, undef, undef, undef, undef);
while(my $line = <MMAP>) { 
  chomp $line;
  if(($line !~ m/^#/) && ($line =~ m/\w/)) { 
    if($line =~ /^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)$/) { 
    #NC_045512           29903  NC_045512-MW422255  29884  11287M9I10469M6I222M3I6277M1I1629M
      ($mfrom, $mfrom_len, $mto, $mto_len, $cigar) = ($1, $2, $3, $4, $5);
      if(($mfrom eq $in_mfrom) && ($mto eq $in_mto)) { 
        # create the map and store this information
        if(defined $tofrom_HA{$mfrom}) { 
          die "ERROR read two lines mapping model $mfrom to $mto in $mmap_file";
        }
        my @tofrom_A = ();
        vdr_CigarToPositionMap(\@tofrom_A, $cigar, $mfrom_len, $mto_len, undef);

        # print map
        printf("#%s\t%s\t%s\t%s\n", "model-from", "model-from-position", "model-to", "model-to-position (if negative:gap in model-to)");
        for(my $i = 1; $i <= $mfrom_len; $i++) { 
          printf("%s\t%s\t%s\t%s\n", $mfrom, $i, $mto, $tofrom_A[$i]);
        }
        close(MMAP);
        exit 0;
      }
    }
    else { 
      die "ERROR unable to parse non-comment and non-blank line in map file: $line\n";
    }
  }
}
close(MMAP);
die "ERROR did not find map info from model $in_mfrom to model $in_mto in mmap file $mmap_file";

exit 0;
