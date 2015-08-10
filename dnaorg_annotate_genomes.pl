#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

# hard-coded-paths:
my $idfetch       = "/netopt/ncbi_tools64/bin/idfetch";
my $esl_fetch_cds = "/panfs/pan1/dnaorg/programs/esl-fetch-cds.pl";
#my $esl_fetch_cds = "/home/nawrocke/notebook/15_0518_dnaorg_virus_compare_script/wd-esl-fetch-cds/esl-fetch-cds.pl";
my $usearch       = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/dnaorg/2015.05.29/usearch";

my $do_uc    = 0; 
my $uc_id    = undef;
my $df_uc_id = 0.9;
my $nuc      = 0;   # number of uclust commands written to uclust script

# The definition of $usage explains the script and usage:
my $usage = "\ndnaorg_compare_genomes.pl\n";
$usage .= "\t<directory created by dnaorg_fetch_dna_wrapper>\n";
$usage .= "\t<list file with all accessions>\n";
$usage .= "\n"; 
$usage .= " This script compares genomes from the same species based\n";
$usage .= " mostly on files containing parsed information from a\n";
$usage .= " 'feature table' file which must already exist. That file is\n";
$usage .= " created by running 'dnaorg_fetch_dna_wrapper.pl -ftable' and\n";
$usage .= " subsequently 'parse-ftable.pl'.\n";
$usage .= "\n";
$usage .= " BASIC OPTIONS:\n";
$usage .= "  -t <f>      : fractional length difference threshold for mismatch [default: 0.1]\n";
$usage .= "  -s          : use short names: all CDS seqs will be named as sequence accession\n";
$usage .= "  -product    : add CDS 'product' qualifier value to output sequence files deflines\n";
$usage .= "  -protid     : add CDS 'protein_id' qualifier value to output sequence files deflines\n";
$usage .= "  -codonstart : add CDS 'protein_id' qualifier value to output sequence files deflines\n";
$usage .= "  -uc         : create a shell script to run uclust jobs for all fasta files\n";
$usage .= "  -uc_id <f>  : with -uc, set sequence identity cutoff to <f> [default: $df_uc_id]\n";
$usage .= "\n";

my ($seconds, $microseconds) = gettimeofday();
my $start_secs    = ($seconds + ($microseconds / 1000000.));
my $executable    = $0;
my $be_verbose    = 1;
my $df_fraclen    = 0.1;
my $fraclen       = undef;
my $do_shortnames = 0; # changed to '1' if -s enabled
my $do_product    = 0; # changed to '1' if -product enabled
my $do_protid     = 0; # changed to '1' if -protid enabled
my $do_codonstart = 0; # changed to '1' if -codonstart enabled

&GetOptions( "t"          => \$fraclen, 
             "s"          => \$do_shortnames, 
             "product"    => \$do_product,
             "protid"     => \$do_protid,
             "codonstart" => \$do_codonstart,
             "uc"         => \$do_uc, 
             "uc_id"      => \$uc_id) || 
    die "Unknown option";

if(scalar(@ARGV) != 2) { die $usage; }
my ($dir, $listfile) = (@ARGV);

$dir =~ s/\/*$//; # remove trailing '/' if there is one

# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
if(defined $fraclen) { 
  $opts_used_short .= "-t $fraclen";
  $opts_used_long  .= "# option:  setting fractional length threshold to $fraclen [-t]\n";
}
if($do_shortnames) { 
  $opts_used_short .= "-s ";
  $opts_used_long  .= "# option:  outputting CDS names as accessions [-s]\n";
}
if($do_product) { 
  $opts_used_short .= "-product ";
  $opts_used_long  .= "# option:  adding \'product\' qualifier values to output sequence file deflines [-product]\n";
}
if($do_protid) { 
  $opts_used_short .= "-protid ";
  $opts_used_long  .= "# option:  adding \'protein_id\' qualifier values to output sequence file deflines [-protid]\n";
}
if($do_codonstart) { 
  $opts_used_short .= "-codonstart ";
  $opts_used_long  .= "# option:  adding \'codon_start\' qualifier values to output sequence file deflines [-codonstart]\n";
}
if($do_uc) { 
  $opts_used_short .= "-uc ";
  $opts_used_long  .= "# option:  outputting uclust shell script [-uc]\n";
}
if(defined $uc_id) { 
  $opts_used_short .= "-uc_id ";
  $opts_used_long  .= "# option:  uclust cluster identity set to $uc_id [-uc_id]\n";
}

# check for incompatible option values/combinations:
if(defined $fraclen && ($fraclen < 0 || $fraclen > 1)) { 
  die "ERROR with -t <f>, <f> must be a number between 0 and 1."; 
}
if((! $do_uc) && (defined $uc_id)) { 
  die "ERROR -uc_id requires -uc";
}

# set default values if user didn't specify otherwise on the command line
if(! defined $fraclen) { 
  $fraclen = $df_fraclen; 
}
if(! defined $uc_id) { 
  $uc_id = $df_uc_id;
}

###############
# Preliminaries
###############
# check if the $dir exists, and that it contains a .gene.tbl file, and a .length file
if(! -d $dir)      { die "ERROR directory $dir does not exist"; }
if(! -s $listfile) { die "ERROR list file $listfile does not exist, or is empty"; }
my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
#my $gene_tbl_file  = $dir . "/" . $dir_tail . ".gene.tbl";
my $cds_tbl_file   = $dir . "/" . $dir_tail . ".CDS.tbl";
my $length_file    = $dir . "/" . $dir_tail . ".length";
my $out_root = $dir . "/" . $dir_tail;
#if(! -s $gene_tbl_file) { die "ERROR $gene_tbl_file does not exist."; }
if(! -s $cds_tbl_file)  { die "ERROR $cds_tbl_file does not exist."; }
if(! -s $length_file)   { die "ERROR $length_file does not exist."; }

# output banner
my $script_name = "dnaorg_compare_genomes.pl";
my $script_desc = "Compare GenBank annotation of genomes";
print ("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
print ("# $script_name: $script_desc\n");
print ("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
print ("# command: $executable $opts_used_short $dir $listfile\n");
printf("# date:    %s\n", scalar localtime());
if($opts_used_long ne "") { 
  print $opts_used_long;
}

#####################
# parse the list file
#####################
my @accn_A = (); # array of accessions
open(IN, $listfile) || die "ERROR unable to open $listfile for reading"; 
my $waccn = 0; # max length of all accessions
while(my $accn = <IN>) { 
  chomp $accn;
  stripVersion(\$accn); # remove version
  push(@accn_A, $accn);
  if(length($accn) > $waccn) { $waccn = length($accn); }
}
close(IN); 

my $head_accn = $accn_A[0];


##################################
# parse the table and length files
##################################
my %gene_tbl_HHA = ();  # Data from .gene.tbl file
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %cds_tbl_HHA = ();   # Data from .cds.tbl file
                        # hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %totlen_H = (); # key: accession, value length read from length file

parseLength($length_file, \%totlen_H);

#parseTable($gene_tbl_file, \%gene_tbl_HHA);
parseTable($cds_tbl_file, \%cds_tbl_HHA);

########
# output 
########
# now create the output
# the verbose output
my $wstrand_str = 0;
my %class_strand_str_H = (); # key: strand string, class number for this strand string 
my %ct_strand_str_H = ();    # key: strand string, value: number of accessions in the class defined by this strand string 
my %idx_strand_str_H = ();   # key: class number,  value: strand string
my %fa_strand_str_H = ();    # key: strand string, value: name of output fasta file for the class defined by this strand string 
my %out_strand_str_HA = ();  # key: strand string, value: array of output strings for the class defined by this strand string 
my @out_fetch_gnm_A = ();    # out_fetch_gnm_A[$c]: for class $c+1, input for idfetch to fetch full genome seqs
my @ct_fetch_gnm_A = ();     # ct_fetch_gnm_A[$c]:  for class $c+1, number of genomes to fetch
my @out_fetch_cds_AA = ();   # out_fetch_cds_AA[$c][$i]: for class $c+1, gene $i+1, input for esl-fetch-cds.pl
my @ct_fetch_cds_AA = ();    # ct_fetch_cds_AA[$c][$i]:  for class $c+1, gene $i+1, number of sequences to fetch 

my $class = undef;
my $nclasses = 0;
my @ngenes_per_class_A = ();
my $max_ngenes = 0;

my $uc_script = $out_root . ".uclust.sh";

for(my $a = 0; $a < scalar(@accn_A); $a++) { 
  my $accn = $accn_A[$a];

  # sanity checks
  if($a == 0 && (! exists $cds_tbl_HHA{$accn})) { die "ERROR didn't read any CDS table information for first accession in $listfile: $accn\n"; } 
  if(! exists $totlen_H{$accn}) { die "ERROR accession $accn does not exist in the length file $length_file"; }

  # set defaults that will stay if we don't have any CDS information
  my $ncds = 0; 
  my $npos = 0;
  my $nneg = 0;
  my $nunc = 0;
  my $nbth = 0; 
  my $strand_str = "";
  my @cds_len_A = ();
  my @cds_coords_A = ();
  my @cds_product_A = ();    # will remain empty unless $do_product is 1 (-product enabled at cmdline)
  my @cds_protid_A = ();     # will remain empty unless $do_protid is 1 (-protid enabled at cmdline)
  my @cds_codonstart_A = (); # will remain empty unless $do_codonstart is 1 (-codonstart enabled at cmdline)
  my $do_desc = ($do_product || $do_protid || $do_codonstart) ? 1 : 0;

  if(exists ($cds_tbl_HHA{$accn})) { 
    ($ncds, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats(\%cds_tbl_HHA, $accn);
    getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $accn, \@cds_len_A, \@cds_coords_A);
    if($do_product) { 
      getQualifierValues(\%cds_tbl_HHA, $accn, "product", \@cds_product_A);
    }
    if($do_protid) { 
      getQualifierValues(\%cds_tbl_HHA, $accn, "protein_id", \@cds_protid_A);
    }
    if($do_codonstart) { 
      getQualifierValues(\%cds_tbl_HHA, $accn, "codon_start", \@cds_codonstart_A);
    }
  }
  if($a == 0) { 
    $wstrand_str = length($strand_str) + 2; 
    if($wstrand_str < length("strand-string")) { $wstrand_str = length("strand_string"); }
  }
  
  if(! exists $class_strand_str_H{$strand_str}) { 
    $nclasses++;
    $class_strand_str_H{$strand_str} = $nclasses;
    $idx_strand_str_H{$nclasses}   = $strand_str;
    $ct_strand_str_H{$strand_str} = 0;
    $fa_strand_str_H{$strand_str} = $dir . "/" . $dir . "." . $nclasses . ".fa";
    @{$out_strand_str_HA{$strand_str}} = ();
  }
  $class = $class_strand_str_H{$strand_str};
  my $c = $class - 1;
  $ct_strand_str_H{$strand_str}++;
  my $ngenes = scalar(@cds_len_A);
  $ngenes_per_class_A[($class-1)] = $ngenes;
  if($ngenes > $max_ngenes) { $max_ngenes = $ngenes; }
  if(! exists $out_fetch_cds_AA[$c]) { @{$out_fetch_cds_AA[$c]} = (); }
  if(! exists $ct_fetch_cds_AA[$c])  { @{$ct_fetch_cds_AA[$c]} = (); }
  
  my $outline = sprintf("%-*s  %5d  %5d  %5d  %5d  %5d  %-*s  %3d  %7d  ", $waccn, $accn, $ncds, $npos, $nneg, $nbth, $nunc, $wstrand_str, $strand_str, $class, $totlen_H{$accn});
  $out_fetch_gnm_A[$c] .= sprintf("%s\n", $accn);
  $ct_fetch_gnm_A[$c]++;

  for(my $i = 0; $i < $ngenes; $i++) { 
    my $desc = "";
    # determine DESCRIPTION: string, esl-fetch-cds.pl will parse this and add it as a description in the defline after fetching the sequence
    if($do_desc) { 
      $desc = "\tDESCRIPTION:";
      if($do_product) { 
        $desc .= "product:";
        if(scalar(@cds_product_A) > 0) { 
          if($i >= scalar(@cds_product_A)) { die "ERROR ran out of products too early for $accn\n"; }
          $desc .= $cds_product_A[$i];
        }
        else { 
          $desc .= "none-annotated";
        }
        if($do_protid || $do_codonstart) { $desc .= " "; }
      }
      if($do_protid) { 
        $desc .= "protein_id:";
        if(scalar(@cds_protid_A) > 0) { 
          if($i >= scalar(@cds_protid_A)) { die "ERROR ran out of protein_ids too early for $accn\n"; }
          $desc .= $cds_protid_A[$i];
        }
        else { 
          $desc .= "none-annotated";
        }
        if($do_codonstart) { $desc .= " "; }
      }
      if($do_codonstart) { 
        $desc .= "codon_start:";
        if(scalar(@cds_codonstart_A) > 0) { 
          if($i >= scalar(@cds_codonstart_A)) { die "ERROR ran out of codon_starts too early for $accn\n"; }
          $desc .= $cds_codonstart_A[$i];
        }
        else { 
          $desc .= "none-annotated";
        }
      }
    }
    $outline .= sprintf("  %5d", $cds_len_A[$i]);
    # create line of input for esl-fetch-cds.pl for fetching the genes of this genome
    if($do_shortnames) { 
      $out_fetch_cds_AA[$c][$i] .= "$accn\t$cds_coords_A[$i]";
    }
    else { 
      $out_fetch_cds_AA[$c][$i] .= sprintf("%s:%s%d:%s%d\t$cds_coords_A[$i]", $head_accn, "class", $class, "gene", ($i+1));
    }
    $out_fetch_cds_AA[$c][$i] .= $desc . "\n"; # $desc will be "" unless $do_product is true
    $ct_fetch_cds_AA[$c][$i]++;
  }
  $outline .= "\n";

  push(@{$out_strand_str_HA{$strand_str}}, $outline);
}
# print header line
printf("#\n");
printf("#%-*s  %5s  %5s  %5s  %5s  %5s  %-*s  %3s  %7s  ", $waccn-1, "accn", "#cds", "#pos", "#neg", "#both", "#unkn", $wstrand_str, "strand-string", "cls", "tot-len");
for(my $i = 0; $i < $max_ngenes; $i++) { 
  printf("  %5s", sprintf("g%d", ($i+1)));
}
printf("\n");

# output stats
for(my $c = 0; $c < $nclasses; $c++) { 
  my $strand_str = $idx_strand_str_H{($c+1)};
  foreach my $outline (@{$out_strand_str_HA{$strand_str}}) { 
    print $outline;
  }
  print "\n";
  @{$out_strand_str_HA{$strand_str}} = (); # clear it for consise output
}

# summarize
#
printf("\n");
printf("# Number-of-classes: $nclasses\n");
printf("# class  #accn  #genes  strand-string\n");
printf("# -----  -----  ------  -------------\n");
my $tot_ngenes = 0;
my $tot_ct = 0;
for(my $c = 0; $c < $nclasses; $c++) { 
  my $strand_str = $idx_strand_str_H{($c+1)};
  printf("%7d  %5d  %6d  %s\n", ($c+1), $ct_strand_str_H{$strand_str}, $ngenes_per_class_A[$c], $strand_str);
  $tot_ngenes += $ngenes_per_class_A[$c];
  $tot_ct     += $ct_strand_str_H{$strand_str};
}
printf("# -----  -----  ------  -------------\n");
printf("%7s  %5d  %6d  %s\n", "total",     $tot_ct,             $tot_ngenes,             "N/A");
printf("%7s  %5.1f  %6.1f  %s\n", "avg",   $tot_ct / $nclasses, $tot_ngenes / $nclasses, "N/A");
printf("\n");

if($do_uc) { 
  open(UC, ">" . $uc_script);
}

# output esl-fetch-cds input, and run esl-fetch-cds.pl for each:
for(my $c = 0; $c < $nclasses; $c++) { 
  # fetch the full genomes
  my $out_fetch_gnm_file   = $out_root . ".c" . ($c+1) . ".fg.idfetch.in";
  my $tmp_out_fetch_gnm_fa = $out_root . ".c" . ($c+1) . ".fg.fa.tmp";
  my $out_fetch_gnm_fa     = $out_root . ".c" . ($c+1) . ".fg.fa";
  open(OUT, ">" . $out_fetch_gnm_file) || die "ERROR unable to open $out_fetch_gnm_file for writing";
  print OUT $out_fetch_gnm_A[$c];
  close OUT;
  sleep(0.1);
  printf("# Fetching %3d full genome sequences for class %2d ... ", $ct_fetch_gnm_A[$c], $c+1);
  my $cmd = "$idfetch -t 5 -c 1 -G $out_fetch_gnm_file > $tmp_out_fetch_gnm_fa";
  runCommand($cmd, 0);
  # now open up the file and change the sequence names manually
  open(IN, $tmp_out_fetch_gnm_fa)    || die "ERROR unable to open $tmp_out_fetch_gnm_fa for reading";
  open(OUT, ">" . $out_fetch_gnm_fa) || die "ERROR unable to open $out_fetch_gnm_fa for writing";
  while(my $line = <IN>) { 
    if($line =~ m/^>/) { 
      chomp $line;
      if($line =~ /^>.*\|(\S+)\|\S*\s+(.+)$/) { 
        print OUT (">$1 $2\n");
      }
      else { die "ERROR unable to parse defline $line in file $tmp_out_fetch_gnm_fa"; }
    }
    else { print OUT $line; }
  }
  close(IN);
  close(OUT);
  unlink $tmp_out_fetch_gnm_fa;
  printf("done. [$out_fetch_gnm_fa]\n");

  # fetch the cds'
  for(my $i = 0; $i < scalar(@{$out_fetch_cds_AA[$c]}); $i++) { 
    my $cg_substr = ".c" . ($c+1) . ".g" . ($i+1);
    my $out_fetch_cds_file = $out_root . $cg_substr . ".esl-fetch-cds.in";
    my $out_fetch_cds_fa   = $out_root . $cg_substr . ".fa";
    my $np_out_fetch_cds_fa    = stripPath($out_fetch_cds_fa);
    my $np_out_centroids_fa    = stripPath($out_root . $cg_substr . ".centroids.fa");
    my $np_out_uclust_sum_file = stripPath($out_root . $cg_substr . ".cluster_summary.txt");
    my $np_out_msa_root        = stripPath($out_root . $cg_substr . ".msa_cluster_");
    open(OUT, ">" . $out_fetch_cds_file) || die "ERROR unable to open $out_fetch_cds_file for writing";
    print OUT $out_fetch_cds_AA[$c][$i];
    close OUT;
    sleep(0.1);
    printf("# Fetching %3d CDS sequences for class %2d gene %2d ... ", $ct_fetch_cds_AA[$c][$i], $c+1, $i+1);
    my $cmd = "";
    if($do_shortnames) { 
      $cmd = "perl $esl_fetch_cds -onlyaccn $out_fetch_cds_file > $out_fetch_cds_fa";
      #printf("$cmd\n");
    }
    else { 
      $cmd = "perl $esl_fetch_cds -nocodon $out_fetch_cds_file > $out_fetch_cds_fa";
    }
    runCommand($cmd, 0);

    printf("done. [$out_fetch_cds_fa]\n");

    if($do_uc){ 
      print UC ("$usearch -cluster_fast $np_out_fetch_cds_fa -id $uc_id -centroids $np_out_centroids_fa -uc $np_out_uclust_sum_file -msaout $np_out_msa_root\n");
      $nuc++;
    }
  }
}

if($do_uc) { 
  close(UC);
  printf("#\n# Shell script for running $nuc usearch commands saved to $uc_script.\n");
}
printf("#[ok]\n");

########################################################
# CURRENTLY NO CONCISE OUTPUT IS PRINTED, BUT I'VE
# LEFT THE BELOW CODE BLOCK FOR REFERENCE AND/OR FUTURE USE
#
# the concise output
my ($ncds0, $npos0, $nneg0, $nunc0, $nbth0, $strand_str0) = getStrandStats(\%cds_tbl_HHA, $head_accn);
my @cds_len0_A = (); 
my @cds_coords0_A = (); 
getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $head_accn, \@cds_len0_A, \@cds_coords0_A);

my $mintotlen = $totlen_H{$head_accn} - ($fraclen * $totlen_H{$head_accn});
my $maxtotlen = $totlen_H{$head_accn} + ($fraclen * $totlen_H{$head_accn});

my @minlen_A = (); # [0..$i..scalar(@cds_len0_A)-1] minimum length for a length match for gene $i
my @maxlen_A = (); # [0..$i..scalar(@cds_len0_A)-1] maximum length for a length match for gene $i
for(my $i = 0; $i < scalar(@cds_len0_A); $i++) { 
  $minlen_A[$i] = $cds_len0_A[$i] - ($fraclen * $cds_len0_A[$i]);
  $maxlen_A[$i] = $cds_len0_A[$i] + ($fraclen * $cds_len0_A[$i]);
}

for(my $a = 0; $a < scalar(@accn_A); $a++) { 
  my $accn = $accn_A[$a];
  
  # set defaults that will stay if we don't have any CDS information
  my $ncds = 0; 
  my $npos = 0;
  my $nneg = 0;
  my $nunc = 0;
  my $nbth = 0; 
  my $strand_str = "";
  my @cds_len_A = ();
  my @cds_coords_A = ();
  if(exists ($cds_tbl_HHA{$accn})) { 
    ($ncds, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats(\%cds_tbl_HHA, $accn);
    getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $accn, \@cds_len_A, \@cds_coords_A);
  }    

  my $output_line = sprintf("%-*s  %s%s%s%s%s%s ", $waccn, $accn,
                            $ncds       == $ncds0   ? "*" : "!",
                            $npos       == $npos0   ? "*" : "!",
                            $nneg       == $nneg0   ? "*" : "!",
                            $nbth       == $nbth0   ? "*" : "!",
                            $nunc       == $nunc0   ? "*" : "!",
                            $strand_str eq $strand_str0 ? "*" : "!");

  if($totlen_H{$accn} >= $mintotlen && $totlen_H{$accn} <= $maxtotlen) { 
    $output_line .= "*";
  }
  else { 
    $output_line .= "!";
  }
  $output_line .= " ";

  for(my $i = 0; $i < scalar(@cds_len0_A); $i++) { 
    if($i < scalar(@cds_len_A)) { 
      if($cds_len_A[$i] >= $minlen_A[$i] && $cds_len_A[$i] <= $maxlen_A[$i]) { 
        $output_line .= "*";
      }
      else { 
        $output_line .= "!"; 
      }
    }
    else { $output_line .= " "; }
  }
  my $pass_or_fail = ($output_line =~ m/\!/) ? "FAIL" : "PASS";

  # print $output_line . " " . $pass_or_fail . "\n";
}


#############
# SUBROUTINES
#############
# Subroutine: runCommand()
# Args:       $cmd:            command to run, with a "system" command;
#             $be_verbose:     '1' to output command to stdout before we run it, '0' not to
#
# Returns:    amount of time the command took, in seconds
# Dies:       if $cmd fails

sub runCommand {
  my $sub_name = "runCommand()";
  my $nargs_exp = 2;

  my ($cmd, $be_verbose) = @_;

  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));
  system($cmd);
  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if($? != 0) { die "ERROR command failed:\n$cmd\n"; }

  return ($stop_time - $start_time);
}

# Subroutine: parseLength()
# Synopsis:   Parses a length file and stores the lengths read
#             into %{$len_HR}.
# Args:       $lenfile: full path to a length file
#             $len_HR:  ref to hash of lengths, key is accession
#
# Returns:    void; fills %{$len_HR}
#
# Dies:       if problem parsing $lenfile

sub parseLength {
  my $sub_name = "parseLength()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($lenfile, $len_HR) = @_;

#HM448898.1	2751

  open(LEN, $lenfile) || die "ERROR unable to open $lenfile for reading";

  while(my $line = <LEN>) { 
    chomp $line;
    my ($accn, $length) = split(/\s+/, $line);
    if($length !~ m/^\d+$/) { die "ERROR couldn't parse length file line: $line\n"; } 

    stripVersion(\$accn);
    $len_HR->{$accn} = $length;
  }
  close(LEN);

  return;
}

# Subroutine: parseTable()
# Synopsis:   Parses a table file and stores the relevant info in it 
#             into $values_HAR.
# Args:       $tblfile:      full path to a table file
#             $values_HHAR:  ref to hash of hash of arrays
#
# Returns:    void; fills @{$values_HHAR}
#
# Dies:       if problem parsing $tblfile

sub parseTable {
  my $sub_name = "parseTable()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($tblfile, $values_HHAR) = @_;

##full-accession	accession	coords	strand	min-coord	gene
#gb|HM448898.1|	HM448898.1	129..476	+	129	AV2

  open(TBL, $tblfile) || die "ERROR unable to open $tblfile for reading";

  # get column header line:
  my $line_ctr = 0;
  my @colnames_A = ();
  my $line = <TBL>;
  my $ncols = undef;
  $line_ctr++;
  if(! defined $line) { die "ERROR did not read any lines from file $tblfile"; }
  chomp $line;
  if($line =~ s/^\#//) { 
    @colnames_A = split(/\t/, $line);
    $ncols = scalar(@colnames_A);
  }
  else { 
    die "ERROR first line of $tblfile did not start with \"#\"";
  }
  if($colnames_A[0] ne "full-accession") { die "ERROR first column name is not full-accession"; }
  if($colnames_A[1] ne "accession")      { die "ERROR second column name is not accession"; }
  if($colnames_A[2] ne "coords")         { die "ERROR third column name is not coords"; }

  # read remaining lines
  while($line = <TBL>) { 
    chomp $line;
    $line_ctr++;
    if($line =~ m/^\#/) { die "ERROR, line $line_ctr of $tblfile begins with \"#\""; }
    my @el_A = split(/\t/, $line);
    if(scalar(@el_A) != $ncols) { 
      die "ERROR, read wrong number of columns in line $line_ctr of file $tblfile";
    }
    my $prv_min_coord = 0;
    # get accession
    my $accn = $el_A[1]; 
    stripVersion(\$accn);
    if(! exists $values_HHAR->{$accn}) { 
      %{$values_HHAR->{$accn}} = (); 
    }

    for(my $i = 0; $i < $ncols; $i++) { 
      my $colname = $colnames_A[$i];
      my $value   = $el_A[$i];
      if($colname eq "min-coord") { 
        if($value < $prv_min_coord) { 
          die "ERROR, minimum coordinates out of order at line $line_ctr and previous line of file $tblfile"; 
        }
        $prv_min_coord = $value; 
        # printf("prv_min_coord: $prv_min_coord\n");
      }

      if(! exists $values_HHAR->{$accn}{$colname}) { 
        @{$values_HHAR->{$accn}{$colname}} = ();
      }
      push(@{$values_HHAR->{$accn}{$colname}}, $el_A[$i]);
      #printf("pushed $accn $colname $el_A[$i]\n");
    }
  }
  close(TBL);
  return;
}

# Subroutine: getStrandStats()
# Synopsis:   Retreive strand stats.
# Args:       $tbl_HHAR:  ref to hash of hash of arrays
#             $accn:      1D key to print for
#
# Returns:    6 values:
#             $nfeatures:  number of features
#             $npos:       number of genes with all segments on positive strand
#             $nneg:       number of genes with all segmenst on negative strand
#             $nunc:       number of genes with all segments on unknown strand 
#             $nbth:       number of genes with that don't fit above 3 categories
#             $strand_str: strand string, summarizing strand of all genes, in order
#
sub getStrandStats {
  my $sub_name = "getStrandStats()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($tbl_HHAR, $accn) = @_;

  my $nfeatures; # number of genes in this genome
  my $npos = 0;  # number of genes on positive strand 
  my $nneg = 0;  # number of genes on negative strand 
  my $nbth = 0;  # number of genes with >= 1 segment on both strands (usually 0)
  my $nunc = 0;  # number of genes with >= 1 segments that are uncertain (usually 0)
  my $strand_str = "";

  if(! exists $tbl_HHAR->{$accn}{"strand"}) { die("ERROR didn't read strand information for accn: $accn\n"); }

  $nfeatures = scalar(@{$tbl_HHAR->{$accn}{"accession"}});
  if ($nfeatures > 0) { 
    for(my $i = 0; $i < $nfeatures; $i++) { 

      # sanity check
      my $accn2 = $tbl_HHAR->{$accn}{"accession"}[$i];
      stripVersion(\$accn2);
      if($accn ne $accn2) { die "ERROR accession mismatch in gene ftable file ($accn ne $accn2)"; }

      if   ($tbl_HHAR->{$accn}{"strand"}[$i] eq "+") { $npos++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "-") { $nneg++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "!") { $nbth++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "?") { $nunc++; }
      else { die("ERROR unable to parse strand for feature %d for $accn\n", $i+1); }
      $strand_str .= $tbl_HHAR->{$accn}{"strand"}[$i];
    }
  }

  return ($nfeatures, $npos, $nneg, $nunc, $nbth, $strand_str);
}


# Subroutine: getLengthStatsAndCoordStrings()
# Synopsis:   Retreive length stats for an accession
#             the length of all annotated genes.
# Args:       $tbl_HHAR:  ref to hash of hash of arrays
#             $accn:      accession we're interested in
#             $len_AR:    ref to array to fill with lengths of features in %{$tbl_HAR}
#             $coords_AR: ref to array to fill with coordinates for each gene
# Returns:    void; fills @{$len_AR} and @{$coords_AR}
#
sub getLengthStatsAndCoordStrings { 
  my $sub_name = "getLengthStatsAndCoordStrings()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($tbl_HHAR, $accn, $len_AR, $coords_AR) = @_;

  if(! exists $tbl_HHAR->{$accn}) { die "ERROR in $sub_name, no data for accession: $accn"; }
  if(! exists $tbl_HHAR->{$accn}{"coords"}) { die "ERROR in $sub_name, no coords data for accession: $accn"; }

  my $ngenes = scalar(@{$tbl_HHAR->{$accn}{"coords"}});

  if ($ngenes > 0) { 
    for(my $i = 0; $i < $ngenes; $i++) { 
      push(@{$len_AR},    lengthFromCoords($tbl_HHAR->{$accn}{"coords"}[$i]));
      push(@{$coords_AR}, addAccnToCoords($tbl_HHAR->{$accn}{"coords"}[$i], $accn));
    }
  }

  return;
}

# Subroutine: getQualifierValues()
# Synopsis:   Retreive values for the qualifier $qualifier in the given %tbl_HHAR
#             and return the values in $values_AR.
#             the length of all annotated genes.
# Args:       $tbl_HHAR:  ref to hash of hash of arrays
#             $accn:      accession we're interested in
#             $qualifier: qualifier we're interested in (e.g. 'Product')
#             $values_AR: ref to array to fill with values of $qualifier
# Returns:    void; fills @{$values_AR}
#
sub getQualifierValues {
  my $sub_name = "getQualifierValues()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($tbl_HHAR, $accn, $qualifier, $values_AR) = @_;

  if(! exists $tbl_HHAR->{$accn}) { die "ERROR in $sub_name, no data for accession: $accn"; }

  if(! exists $tbl_HHAR->{$accn}{$qualifier}) { return; } # no annotation for $qualifier, do not update arrays

  my $nvalues = scalar(@{$tbl_HHAR->{$accn}{$qualifier}});

  if ($nvalues > 0) { 
    for(my $i = 0; $i < $nvalues; $i++) { 
      push(@{$values_AR},  $tbl_HHAR->{$accn}{$qualifier}[$i]);
    }
  }

  return;
}


# Subroutine: lengthFromCoords()
# Synopsis:   Determine the length of a region give its coords in NCBI format.
#
# Args:       $coords:  the coords string
#
# Returns:    length in nucleotides implied by $coords  
#
sub lengthFromCoords { 
  my $sub_name = "lengthFromCoords()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($coords) = @_;

  my $orig_coords = $coords;
  # Examples:
  # complement(2173412..2176090)
  # complement(join(226623..226774, 226854..229725))

  # remove 'complement('  ')'
  $coords =~ s/^complement\(//;
  $coords =~ s/\)$//;

  # remove 'join('  ')'
  $coords =~ s/^join\(//;
  $coords =~ s/\)$//;

  my @el_A = split(/\s*\,\s*/, $coords);

  my $length = 0;
  foreach my $el (@el_A) { 
    # rare case: remove 'complement(' ')' that still exists:
    $el =~ s/^complement\(//;
    $el =~ s/\)$//;
    $el =~ s/\<//; # remove '<'
    $el =~ s/\>//; # remove '>'
    if($el =~ m/^(\d+)\.\.(\d+)$/) { 
      my ($start, $stop) = ($1, $2);
      $length += abs($start - $stop) + 1;
    }
    else { 
      die "ERROR unable to parse $orig_coords in $sub_name"; 
    }
  }

  # printf("in lengthFromCoords(): orig_coords: $orig_coords returning length: $length\n");
  return $length;
}

# Subroutine: addAccnToCoords()
# Synopsis:   Add accession Determine the length of a region give its coords in NCBI format.
#
# Args:       $coords:  the coords string
#             $accn:    accession to add
# Returns:    The accession to add to the coords string.
#
sub addAccnToCoords { 
  my $sub_name = "addAccnToCoords()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($coords, $accn) = @_;

  my $ret_coords = $coords;
  # deal with simple case of \d+..\d+
  if($ret_coords =~ /^\<?\d+\.\.\>?\d+/) { 
    $ret_coords = $accn . ":" . $ret_coords;
  }
  # replace 'complement(\d' with 'complement($accn:\d+'
  while($ret_coords =~ /complement\(\<?\d+/) { 
    $ret_coords =~ s/complement\((\<?\d+)/complement\($accn:$1/;
  }
  # replace 'join(\d' with 'join($accn:\d+'
  while($ret_coords =~ /join\(\<?\d+/) { 
    $ret_coords =~ s/join\((\<?\d+)/join\($accn:$1/;
  }
  # replace ',\d+' with ',$accn:\d+'
  while($ret_coords =~ /\,\s*\<?\d+/) { 
    $ret_coords =~ s/\,\s*(\<?\d+)/\,$accn:$1/;
  }

  #print("addAccnToCoords(), input $coords, returning $ret_coords\n");
  return $ret_coords;
}

# Subroutine: stripVersion()
# Purpose:    Given a ref to an accession.version string, remove the version.
# Args:       $accver_R: ref to accession version string
# Returns:    Nothing, $$accver_R has version removed
sub stripVersion {
  my $sub_name  = "stripVersion()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($accver_R) = (@_);

  $$accver_R =~ s/\.[0-9]*$//; # strip version

  return;
}

# Subroutine: stripPath()
# Purpose:    Given a file path, remove the all directories and leave only the file name.
# Args:       $filename: full path to file
# Returns:    only the file name, without any directory structure
sub stripPath {
  my $sub_name  = "stripPath()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($filename) = (@_);

  $filename =~ s/^.+\///;

  return $filename;
}
