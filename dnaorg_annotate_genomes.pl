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

# The definition of $usage explains the script and usage:
my $usage = "\ndnaorg_annotate_genomes.pl\n";
$usage .= "\t<directory created by dnaorg_fetch_dna_wrapper>\n";
$usage .= "\t<list file with all accessions>\n";
$usage .= "\n"; 
$usage .= " This script annotates genomes from the same species based\n";
$usage .= " on reference annotation.\n";
$usage .= "\n";

my ($seconds, $microseconds) = gettimeofday();
my $start_secs    = ($seconds + ($microseconds / 1000000.));
my $executable    = $0;

#&GetOptions("noexp"      => \$do_noexp) || 
#    die "Unknown option";

if(scalar(@ARGV) != 2) { die $usage; }
my ($dir, $listfile) = (@ARGV);

#$dir =~ s/\/*$//; # remove trailing '/' if there is one
#my $outdir     = $dir;
#my $outdirroot = $outdir;
#$outdirroot =~ s/^.+\///;

# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
# NONE YET: one left commented out below for future reference
#if(defined $do_noexp) { 
#  $opts_used_short .= "-noexp ";
#  $opts_used_long  .= "# option:  verbose explanation of class definition not being output to stdout [-noexp]\n";
#}
# 
# check for incompatible option values/combinations:
# NONE YET

###############
# Preliminaries
###############
# check if the $dir exists, and that it contains a .gene.tbl file, and a .length file
if(! -d $dir)      { die "ERROR directory $dir does not exist"; }
if(! -s $listfile) { die "ERROR list file $listfile does not exist, or is empty"; }
my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $gene_tbl_file  = $dir . "/" . $dir_tail . ".gene.tbl";
my $cds_tbl_file   = $dir . "/" . $dir_tail . ".CDS.tbl";
my $length_file    = $dir . "/" . $dir_tail . ".length";
my $out_root = $dir . "/" . $dir_tail;
#if(! -s $gene_tbl_file) { die "ERROR $gene_tbl_file does not exist."; }
if(! -s $cds_tbl_file)  { die "ERROR $cds_tbl_file does not exist."; }
if(! -s $length_file)   { die "ERROR $length_file does not exist."; }

# output banner
my $script_name = "dnaorg_annotate_genomes.pl";
my $script_desc = "Annotate genomes based on a reference and homology search";
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

parseTable($gene_tbl_file, \%gene_tbl_HHA);
parseTable($cds_tbl_file, \%cds_tbl_HHA);

###########################################################################
# for each non-reference accession: 
#   for each reference CDS: 
#     look for homologous CDS in genome
###########################################################################

###########################################################################
#my $wstrand_str = 0;         # for output formatting, width of max expected length strand string
#my $wlabel_str = 0;          # for output formatting, width of max expected length label string
#my %label_str2idx_H = ();    # key: label string, class number for this label string 
#my %idx2label_str_H = ();    # key: class number,  value: label string
#my %ct_label_str_H = ();     # key: label string, value: number of accessions in the class defined by this label string 
#my %fa_label_str_H = ();     # key: strand string, value: name of output fasta file for the class defined by this label string 
#my %out_label_str_HA = ();   # key: strand string, value: array of output strings for the class defined by this label string 
#my @out_fetch_gnm_A = ();    # out_fetch_gnm_A[$c]: for class $c+1, input for idfetch to fetch full genome seqs
#my @ct_fetch_gnm_A = ();     # ct_fetch_gnm_A[$c]:  for class $c+1, number of genomes to fetch
#my @out_fetch_cds_AA = ();   # out_fetch_cds_AA[$c][$i]: for class $c+1, gene $i+1, input for esl-fetch-cds.pl
#my @ct_fetch_cds_AA = ();    # ct_fetch_cds_AA[$c][$i]:  for class $c+1, gene $i+1, number of sequences to fetch 
#my %accn2label_str_H = ();   # key: $accn, value: $label_str for that accn
#my $label_str = "";          # a label string, e.g. 'ABCD:===='
my $strand_str;              # +/- string for all CDS for an accession: e.g. '+-+': 1st and 3rd CDS are + strand, 2nd is -
#my $class_idx = undef;       # class index for current accession
#my $nclasses = 0;            # total number of classes
#my @ncds_per_class_A = ();   # number of CDS per class
#my $max_ncds = 0;            # max number of CDS for any class

# variables related to a reference accession
my $ref_accn      = undef; # changed to <s> with -ref <s>
my $ref_label_str = undef; # label string for reference accn

#my $uc_script = $out_root . ".uclust.sh";

# reference information on reference accession, first accession read in ntlist file
my $ref_ncds          = 0;  # number of CDS in reference
my $ref_strand_str    = ""; # strand string for reference 
my @ref_cds_len_A     = (); # [0..$i..$ref_ncds-1]: length of each reference CDS
#my @ref_cds_len_tol_A = (); # [0..$i..$ref_ncds-1]: length tolerance, any gene that is within this fraction of the lenght of the ref gene is a match
my @ref_cds_coords_A  = (); # [0..$i..$ref_ncds-1]: CDS coords for reference
my @ref_cds_product_A = (); # CDS:product qualifier data for reference 

my $ncds = 0;              # number of CDS
my $npos = 0;              # number of CDS on positive strand
my $nneg = 0;              # number of CDS on negative strand
my $nunc = 0;              # number of CDS on uncertain strand
my $nbth = 0;              # number of CDS on both strands
my @cds_len_A = ();        # [0..$i..$ncds-1] length of CDS $i
my @cds_coords_A = ();     # [0..$i..$ncds-1] coords of CDS $i
my @cds_product_A = ();    # [0..$i..$ncds-1] CDS:product annotation for CDS $i
my @cds_protid_A = ();     # will remain empty unless $do_protid is 1 (-protid enabled at cmdline)
my @cds_codonstart_A = (); # will remain empty unless $do_codonstart is 1 (-codonstart enabled at cmdline)
#my $do_desc = ($do_product || $do_protid || $do_codonstart) ? 1 : 0; # '1' to create a description to add to defline of fetched sequences, '0' not to

# Get information on the reference accession
my $naccn = scalar(@accn_A);
$ref_accn = $accn_A[0];
if(! exists ($cds_tbl_HHA{$ref_accn})) { die "ERROR no CDS information stored for reference accession"; }
(undef, undef, undef, undef, undef, $ref_strand_str) = getStrandStats(\%cds_tbl_HHA, $ref_accn);
getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $ref_accn, \@ref_cds_len_A, \@ref_cds_coords_A);
getQualifierValues(\%cds_tbl_HHA, $ref_accn, "product", \@ref_cds_product_A);
$ref_ncds = scalar(@ref_cds_len_A);

# Create input query files for each of the reference CDS:
for(my $i = 0; $i < $ref_ncds; $i++) { 
  printf("# Fetching reference CDS sequence %d ... ", ($i+1));

  my $strand = substr($ref_strand_str, $i, 1);
  my $coords_with_accn = addAccnToCoords($ref_cds_coords_A[$i], $ref_accn);

  my $fetch_input = $out_root . ".ref.cds." . ($i+1) . ".esl-fetch.in";
  my $fetch_fasta = $out_root . ".ref.cds." . ($i+1) . ".fa";
  open(FOUT, ">" . $fetch_input) || die "ERROR unable to open $fetch_input for writing";
  printf FOUT ("$ref_accn:cds%d\t$coords_with_accn\n", ($i+1));
  close FOUT;

  my $cmd = "perl $esl_fetch_cds -nocodon $fetch_input > $fetch_fasta";
  runCommand($cmd, 0);
      
  printf("done. [$fetch_fasta]\n");

  # if CDS has multiple exons, fetch each independently into it's own file
  my @starts_A = ();
  my @stops_A  = ();
  my $nexons   = 0;
  startStopsFromCoords($ref_cds_coords_A[$i], \@starts_A, \@stops_A, \$nexons);
  if($nexons > 1) { 
    for(my $e = 0; $e < $nexons; $e++) { 
      printf("# Fetching reference CDS %d sequence, exon %d ... ", ($i+1), ($e+1));
      $fetch_input = $out_root . ".ref.cds." . ($i+1) . ".exon." . ($e+1) . ".esl-fetch.in";
      $fetch_fasta = $out_root . ".ref.cds." . ($i+1) . ".exon." . ($e+1) . ".fa";
      my $fetch_string = "";
      $fetch_string = $ref_accn . ":" . $starts_A[$e] . ".." . $stops_A[$e];
      if($strand eq "-") { # reverse strand
        $fetch_string = "complement(" . $fetch_string . ")";
      }
      elsif($strand ne "+") { 
        die "ERROR reference CDS $i has strand $strand, which we can't deal with";
      }
      open(FOUT, ">" . $fetch_input) || die "ERROR unable to open $fetch_input for writing";
      printf FOUT ("$ref_accn:cds%d:exon%d\t$fetch_string\n", ($i+1), ($e+1));
      close FOUT;
      my $cmd = "perl $esl_fetch_cds -nocodon $fetch_input > $fetch_fasta";
      runCommand($cmd, 0);
      printf("done. [$fetch_fasta]\n");
    }
  }      
}

# validate all reference CDS have unique CDS:product annotation
validateRefCDSAreUnique($ref_ncds, \@ref_cds_product_A);

print("#\n");
printf("# Reference accession: $ref_accn\n");
#print("#\n");
#for(my $rc = 0; $rc < $ref_ncds; $rc++) { 
#$ref_cds_len_tol_A[$rc] = $fraclen * $ref_cds_len_A[$rc];
#printf("# Reference CDS product %d: %-30s reference-length: %6d   length-range-for-match: %8.1f to %8.1f\n", $rc+1, $ref_cds_product_A[$rc], $ref_cds_len_A[$rc], ($ref_cds_len_A[$rc] - $ref_cds_len_tol_A[$rc]), ($ref_cds_len_A[$rc] + $ref_cds_len_tol_A[$rc]));
#}

#$wlabel_str  = (2*$ref_ncds+1) * 2;
#$wstrand_str = $ref_ncds * 2;
#if($wlabel_str  < length("label"))         { $wlabel_str  = length("label"); }
#if($wstrand_str < length("strand-string")) { $wstrand_str = length("strand_string"); }
# }    

# Fetch all genome sequences, including the reference
my $gnm_fetch_file = $out_root . ".fg.idfetch.in";
my $gnm_fasta_file = $out_root . ".fg.fa";
open(OUT, ">" . $gnm_fetch_file) || die "ERROR unable to open $gnm_fetch_file";
for(my $a = 0; $a < $naccn; $a++) { 
#  print OUT $accn_A[$a] . "\n";
  my $fetch_string = "join(" . $accn_A[$a] . ":1.." . $totlen_H{$accn_A[$a]} . "," . $accn_A[$a] . ":1.." . $totlen_H{$accn_A[$a]} . ")\n";
  print OUT $accn_A[$a] . ":" . "genome-duplicated" . "\t" . "join($accn_A[$a]:1..$totlen_H{$accn_A[$a]},$accn_A[$a]:1..$totlen_H{$accn_A[$a]})\n";
}
close(OUT);

printf("# Fetching $naccn full (duplicated) genome sequences... ");
# my $cmd = "$idfetch -t 5 -c 1 -G $gnm_fetch_file > $gnm_fasta_file";
my $cmd = "perl $esl_fetch_cds -nocodon $gnm_fetch_file > $gnm_fasta_file";
runCommand($cmd, 0);
printf("done. [$gnm_fasta_file]\n");

# Pass through all accessions, and annotate them
for(my $a = 0; $a < $naccn; $a++) { 
  my $accn = $accn_A[$a];
  # sanity checks
  if(! exists $totlen_H{$accn}) { die "ERROR accession $accn does not exist in the length file $length_file"; }

  # set defaults that will stay if we don't have any CDS information
  $ncds = 0; 
  $npos = 0;
  $nneg = 0;
  $nunc = 0;
  $nbth = 0; 
  $strand_str = "";
  @cds_len_A = ();
  @cds_coords_A = ();
  @cds_product_A = ();    # will remain empty unless $do_product is 1 (-product enabled at cmdline)
  @cds_protid_A = ();     # will remain empty unless $do_protid is 1 (-protid enabled at cmdline)
  @cds_codonstart_A = (); # will remain empty unless $do_codonstart is 1 (-codonstart enabled at cmdline)

  if(exists ($cds_tbl_HHA{$accn})) { 
    ($ncds, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats(\%cds_tbl_HHA, $accn);
    getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $accn, \@cds_len_A, \@cds_coords_A);
    getQualifierValues(\%cds_tbl_HHA, $accn, "product", \@cds_product_A);
  }
} # end of first pass through all accessions 'for(my $a = 0; $a < $naccn; $a++)'

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
      push(@{$coords_AR}, $tbl_HHAR->{$accn}{"coords"}[$i]);
      #push(@{$coords_AR}, addAccnToCoords($tbl_HHAR->{$accn}{"coords"}[$i], $accn));
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

# Subroutine: readInfile()
# Purpose:    Read an input file specified with -i <s> and store the information therein.
# Args:       $infile: name of the infile to read
#             $in_cds_product_HAR: reference to a hash of arrays to store CDS:product info
#             $in_cds_len_HAR:     reference to a hash of arrays to store CDS length info
# Returns:    void
sub readInfile {
  my $sub_name  = "readInfile()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($infile, $in_cds_product_HAR, $in_cds_len_HAR) = (@_);

  open(IN, $infile) || die "ERROR unable to open $infile for reading.";
  my $ctr = 0;
  while(my $line = <IN>) { 
    $ctr++;
    chomp $line;
    $line =~ s/^\s+//; # remove leading  whitespace
    $line =~ s/\s+$//; # remove trailing whitespace
    if($line !~ m/^\#/ && $line =~ m/\w/) { 
      if($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$/) { 
        my ($type, $idx, $qual, $qval, $value) = ($1, $2, $3, $4, $5);
        if($type ne "ref" && $type ne "alw" && $type ne "kma") { 
          die "ERROR first token not 'ref' nor 'alw' nor 'kma' in infile: $infile on line $ctr:\n$line\n"; 
        }
        if($qual ne "CDS") { 
          die "ERROR read non-CDS qualifier in infile: $infile on line $ctr:\n$line\n"; 
        }
        if($qval eq "product") { 
          if(! exists $in_cds_product_HAR->{$type}) { 
            @{$in_cds_product_HAR->{$type}} = ();
          }
          $in_cds_product_HAR->{$type}[$idx-1] = $value;
        }
        elsif($qval eq "length") { 
          if(! exists $in_cds_len_HAR->{$type}) { 
            @{$in_cds_len_HAR->{$type}} = ();
          }
          $in_cds_len_HAR->{$type}[$idx-1] = $value;
        }
        else { 
          die "ERROR unexpected qualifier value $qval in infile: $infile on line $ctr:\n$line\n"; 
        }
      }
      else { 
        die "ERROR read less than 5 whitespace-delimited tokens in infile: $infile on line $ctr:\n$line\n"; 
      }
    }
  } # end of 'while ($line = <IN>)'
  return;
}

# Subroutine: validateInfileGivenRefInfo()
# Purpose:    Given information read from the parsed ftable about the reference accession, 
#             validate the information we read from the infile with the -i option.
# Args:       $ref_ncds:           number of reference CDS
#             $ref_cds_product_AR: ref to array of CDS:product annotations for the $ref_ncds reference CDS 
#             $ref_cds_len_AR:     ref to array of lengths for the $ref_ncds reference CDS 
#             $in_cds_product_HAR: ref to hash of arrays of CDS:product annotations read from the infile (-i)
#             $in_cds_len_HAR:     ref to hash of arrays of length annotations read from the infile (-i)
# Returns:    void
# Dies:       if we didn't read product and length values for all reference accession from the infile
#             if product and length values read from the infile don't match reference accession info
sub validateInfileGivenRefInfo {
  my $sub_name  = "validateInfileGivenRefInfo()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ref_ncds, $infile, $ref_cds_product_AR, $ref_cds_len_AR, $in_cds_product_HAR, $in_cds_len_HAR) = @_;

  if(! defined $in_cds_product_HAR)           { die "ERROR didn't read any CDS:product information from infile $infile\n"; }
  if(! defined $in_cds_len_HAR)               { die "ERROR didn't read any CDS length information from infile $infile\n"; }
  if(! defined $in_cds_product_HAR->{"ref"})  { die "ERROR didn't read any reference CDS:product information from infile $infile\n"; }
  if(! defined $in_cds_len_HAR->{"ref"})      { die "ERROR didn't read any reference CDS length information from infile $infile\n"; }
  if(scalar(@{$in_cds_product_HAR->{"ref"}}) ne $ref_ncds) { die "ERROR didn't read correct number of CDS:product values from infile $infile\n"; }
  if(scalar(@{$in_cds_len_HAR->{"ref"}})     ne $ref_ncds) { die "ERROR didn't read correct number of CDS length values from infile $infile\n"; }
  for(my $i = 1; $i <= $ref_ncds; $i++) { 
    printf("$sub_name validating CDS $i ... ");
    if(! defined $in_cds_product_HAR->{"ref"}[($i-1)]) { 
      die "ERROR CDS:product information not read for CDS $i in infile: $infile"; 
    }
    if(! defined $in_cds_len_HAR->{"ref"}[($i-1)]) {
      die "ERROR CDS length information not read for CDS $i in infile: $infile"; 
    } 
    if($in_cds_product_HAR->{"ref"}[($i-1)] ne $ref_cds_product_AR->[($i-1)]) { 
      die sprintf("ERROR CDS:product information for CDS reference accession does not match feature table for CDS $i in infile: $infile\ninfile: %s\nftable: %s\n", $in_cds_product_HAR->{"ref"}[($i-1)], $ref_cds_product_AR->[($i-1)]);
    }
    if($in_cds_len_HAR->{"ref"}[($i-1)] ne $ref_cds_len_AR->[($i-1)]) { 
      die sprintf("ERROR CDS length information for CDS reference accession does not match feature table for CDS $i in infile: $infile\ninfile: %d\nftable: %d", $in_cds_len_HAR->{"ref"}[($i-1)], $ref_cds_len_AR->[($i-1)]);
                  
    }
    printf(" done.\n");
  }
  return;
}

# Subroutine: validateRefCDSAreUnique()
# Purpose:    Validate that all CDS:product annotation for all reference CDS 
#             are unique, i.e. there are no two CDS that have the same value
#             in their CDS:product annotation.
# Args:       $ref_ncds:           number of reference CDS
#             $ref_cds_product_AR: ref to array of CDS:product annotations for the $ref_ncds reference CDS 
# Returns:    void
# Dies:       if more than one ref CDS have same CDS:product annotation.
sub validateRefCDSAreUnique {
  my $sub_name  = "validateRefCDSAreUnique()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ref_ncds, $ref_cds_product_AR) = @_;

  my %exists_H = ();
  for(my $i = 0; $i < $ref_ncds; $i++) { 
    if(exists $exists_H{$ref_cds_product_AR->[$i]}) { die sprintf("ERROR %s is CDS:product value for more than one reference CDS!", $ref_cds_product_AR->[$i]); }
  }

  return;
}

# Subroutine: startStopsFromCoords()
# Synopsis:   Extract the starts and stops from a coords string.
#
# Args:       $coords:  the coords string
#             $starts_AR: ref to array to fill with start positions
#             $stops_AR:  ref to array to fill with stop positions
#             $nexons_R:  ref to scalar that fill with the number of exons
#
# Returns:    void; but fills
#
sub startStopsFromCoords { 
  my $sub_name = "startStopsFromCoords()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($coords, $starts_AR, $stops_AR, $nexons_R) = @_;

  @{$starts_AR} = ();
  @{$stops_AR}  = ();
  $$nexons_R    = 0;
  
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
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $2);
      $$nexons_R++;
    }
    else { 
      die "ERROR unable to parse $orig_coords in $sub_name"; 
    }
  }

  # printf("in startStopsFromCoords(): orig_coords: $orig_coords returning length: $length\n");
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

  my @starts_A = ();
  my @stops_A  = ();
  my $nexons   = 0;

  startStopsFromCoords($coords, \@starts_A, \@stops_A, \$nexons);

  my $length = 0;
  for(my $i = 0; $i < $nexons; $i++) { 
    $length += abs($starts_A[$i] - $stops_A[$i]) + 1;
  }

  return $length;
}
