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
$usage  = "vadr-map-model-coords.pl\n\n";
$usage .= "Usage:\n\n";
$usage .= "perl vadr-map-model-coords.pl [OPTIONS]\n\t<detailed-error-report.tsv or .alt.list file path>\n\t<mmap file path>\n\t<model to map positions to>\n\n";
$usage .= "\tOPTIONS:\n";
$usage .= "\t\t--alt : file listed in first argument is an .alt file (not a detailed-error-report.tsv or .alt.list file)\n\n";

my $do_alt = 0; 
&GetOptions( "alt" => \$do_alt);

if(scalar(@ARGV) != 3) { die $usage; }
my ($in_file, $mmap_file, $in_mto) = @ARGV;

if(! -s $mmap_file) { 
  die "ERROR model map file $mmap_file does not exist";
}
if(! -e $in_file) { 
  die "ERROR input file $in_file does not exist";
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
      if($mto eq $in_mto) { 
        # create the map and store this information
        if(defined $tofrom_HA{$mfrom}) { 
          die "ERROR read two lines mapping model $mfrom to $mto in $mmap_file";
        }
        @{$tofrom_HA{$mfrom}} = ();
        vdr_CigarToPositionMap(\@{$tofrom_HA{$mfrom}}, $cigar, $mfrom_len, $mto_len, undef);

        # for debugging: print map
        #for(my $i = 1; $i <= $mfrom_len; $i++) { 
        #  printf("MAP from %s to %s [%4d]: %4d\n", $mfrom, $mto, $i, $tofrom_HA{$mfrom}[$i]);
        #}
      }
    }
    else { 
      die "ERROR unable to parse non-comment and non-blank line in map file: $line\n";
    }
  }
}

my $new_mcoords_header = $in_mto . "-mdl-coords(mapped)";

# variables used only if $do_alt
my @orig_line_A  = (); # only used if $do_alt, because we need to know width of first column before outputting
my @new_mcoords_A = (); # only used if $do_alt, because we need to know width of first column before outputting
my $w_new_mcoords = length("#" . $new_mcoords_header);

# open .alt.list or .alt file and output 
open(IN, $in_file) || die "ERROR, unable to open $in_file for readin";
while(my $line = <IN>) { 
  chomp $line;
  my $orig_line = $line;
  my @el_A = ();
  my ($mfrom, $orig_mcoords);
  if(($line =~ m/^\#/) || ($line =~ m/^sequence\tmodel\tfeature-type/)) { 
    # header line, make a new one
    if(! $do_alt) { # default: alt.list
      print $orig_line . "\t" . $new_mcoords_header . "\n";
    }
    else { # $do_alt == 1, .alt file
      $orig_line =~ s/^\#/ /;
      push(@orig_line_A, $orig_line);
      if(scalar(@orig_line_A) == 2) { # only print header once
        push(@new_mcoords_A, "#" . $new_mcoords_header); # the new_mcoords header
      }
      else { 
        push(@new_mcoords_A, "#");
      }
    }
  }
  else { 
    if(! $do_alt) { 
      @el_A = split(/\t+/, $line);
      
      ##sequence	model	feature-type	feature-name	error	seq-coords	mdl-coords	error-description
      #MW691153	NC_045512-del28254	-	*sequence*	DELETION_OF_FEATURE	-	27756..27887:+	internal deletion of a complete feature [gene feature number 8: ORF7b]

      if(scalar(@el_A) != 8) { 
        die "ERROR: expected 8 tab-delimited tokens on each line of $in_file, but got " . scalar(@el_A) . " on line:\n$line\n";
      }
      ($mfrom, $orig_mcoords) = ($el_A[1], $el_A[6]);
    }
    else { # $do_alt == 1, .alt file
      @el_A = split(/\s+/, $line);
      ##        seq                           ftr   ftr            ftr  alert           alert                                                              seq    seq                                mdl    mdl  alert 
      ##idx     name      model               type  name           idx  code      fail  description                                                     coords    len                             coords    len  detail
      ##------  --------  ------------------  ----  -------------  ---  --------  ----  ----------------------------  ----------------------------------------  -----  ---------------------------------  -----  ------
      #1.1.1    MW286742  NC_045512           CDS   ORF8_protein    19  cdsstopn  no    CDS_HAS_STOP_CODON                                      27907..27909:+      3                     27945..27947:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TAA, shifted S:312,M:312]
      if(scalar(@el_A) < 14) { 
        die "ERROR, expected at least 14 whitespace-delimited tokens on each line of $in_file, but got " . scalar(@el_A) . " on line:\n$line\n";
      }
      ($mfrom, $orig_mcoords) = ($el_A[2], $el_A[11]);
    }

    my $new_mcoords = "";
    if($mfrom ne $in_mto) { 
      # model other than $mto, map coords
      if(! defined $tofrom_HA{$mfrom}) { 
        die "ERROR: read model $mfrom for which no map existed to $mto in input model map file on line:\n$line\n";
      }
      if($orig_mcoords eq "-") { 
        $new_mcoords = "-";
      }
      else { 
        my @orig_start_A = ();
        my @orig_stop_A  = ();
        my @strand_A     = ();
        vdr_FeatureStartStopStrandArrays($orig_mcoords, \@orig_start_A, \@orig_stop_A, \@strand_A, undef);
        my $nsgm = scalar(@orig_start_A);
        for(my $i = 0; $i < $nsgm; $i++) { 
          my $new_start = abs($tofrom_HA{$mfrom}[($orig_start_A[$i])]); 
          my $new_stop  = abs($tofrom_HA{$mfrom}[($orig_stop_A[$i])]);
          # abs() gets rid of negative value, user won't know it's a gap in $mto, I guess this is ok
          
          if($new_mcoords ne "") { $new_mcoords .= ","; }
          $new_mcoords .= vdr_CoordsSegmentCreate($new_start, $new_stop, $strand_A[$i], undef);
        }
        #printf("orig_mcoords: $orig_mcoords, new_mcoords: $new_mcoords\n");
      }
    }
    else { # same model
      $new_mcoords = $orig_mcoords;
    }
    if(! $do_alt) { 
      print $orig_line . "\t" . $new_mcoords . "\n";
    }
    else { # $do_alt == 1
      push(@orig_line_A, $orig_line);
      push(@new_mcoords_A, $new_mcoords); # the new_mcoords header
      if($w_new_mcoords < length($new_mcoords)) { 
        $w_new_mcoords = length($new_mcoords);
      }
    }
  }
}
if($do_alt) { # output all lines, now that we know max width
  for(my $i = 0; $i < scalar(@orig_line_A); $i++) { 
    if($orig_line_A[$i] ne "#") { 
      if(($new_mcoords_A[$i] =~ m/^\#/) && ($i == 2)) { # special case, make monocharacter string
        printf("%-*s  %s\n", $w_new_mcoords, "#" . utl_StringMonoChar($w_new_mcoords-1, "-", undef), $orig_line_A[$i]);
      }
      else { 
        printf("%-*s  %s\n", $w_new_mcoords, $new_mcoords_A[$i], $orig_line_A[$i]);
      }
    }
    else { 
      print $orig_line_A[$i] . "\n";
    }
  }
}
  
