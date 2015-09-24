#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

# hard-coded-paths:
my $idfetch        = "/netopt/ncbi_tools64/bin/idfetch";
my $esl_fetch_cds  = "/panfs/pan1/dnaorg/programs/esl-fetch-cds.pl";
my $esl_ssplit     = "/panfs/pan1/dnaorg/programs/Bio-Easel/scripts/esl-ssplit.pl";
my $esl_translate  = "/home/nawrocke/src/wd-easel/miniapps/esl-translate";
my $hmmer_exec_dir = "/home/nawrocke/bin/";
my $inf_exec_dir   = "/usr/local/infernal/1.1.1/bin/";

# The definition of $usage explains the script and usage:
my $usage = "\ndnaorg_annotate_genomes.pl\n";
$usage .= "\t<directory created by dnaorg_virus_wrapper.pl>\n";
$usage .= "\t<list file with all accessions>\n";
$usage .= "\n"; 
$usage .= " This script annotates genomes from the same species based\n";
$usage .= " on reference annotation. The reference accession is the\n";
$usage .= " first accession listed in <list file with all accessions>\n";
$usage .= "\n";
$usage .= " BASIC OPTIONS:\n";
$usage .= "  -nocorrect : do not correct annotations based on internal start/stop codons in predicted exons/CDS\n";
$usage .= "  -matpept   : use mat_peptide info instead of CDS info\n";
$usage .= "  -oseq <s>  : identify origin sequence <s> in genomes, put | at site of origin, e.g. TAATATT|AC\n";
$usage .= "  -strict    : require matching annotations to match CDS/exon index\n";
$usage .= "  -nodup     : do not duplicate each genome to allow identification of features that span stop..start\n";
$usage .= "  -notexon   : do not use exon-specific models\n";
$usage .= "  -onlybuild : exit after building reference models\n";
$usage .= "  -model <s> : use model file <s>, instead of building one\n";
$usage .= "\n OPTIONS CONTROLLING OUTPUT TABLE:\n";
$usage .= "  -c         : concise output mode (enables -nomdlb -noexist -nobrack and -nostop)\n";
$usage .= "  -nomdlb    : do not add model boundary annotation to output\n";
$usage .= "  -noexist   : do not include information on existing annotation\n";
$usage .= "  -nobrack   : do not include brackets around predicted annotations that do not match existing\n";
$usage .= "  -nostop    : do not output stop codon for each predicted CDS\n";
$usage .= "  -nofid     : do not output fractional identity relative to reference for each CDS/exon\n";
$usage .= "  -noss3     : do not output results of start codon, stop codon and multiple of 3 tests\n";
$usage .= "  -noolap    : do not output overlap information\n";
$usage .= "  -noexp     : do not output explanation of column headings\n";
$usage .= "\n OPTIONS FOR SELECTING HOMOLOGY SEARCH ALGORITHM:\n";
$usage .= "  -hmmer     : use HMMER for predicting annotations, default: use Infernal\n";
$usage .= "\n OPTIONS SPECIFIC TO HMMER3:\n";
$usage .= "  -hmmenv  : use HMM envelope boundaries for predicted annotations, default: use window boundaries\n";
$usage .= "\n OPTIONS SPECIFIC TO INFERNAL:\n";
$usage .= "  -iglocal      : use the -g option with cmsearch for glocal searches\n";
$usage .= "  -cslow        : use default cmcalibrate parameters, not parameters optimized for speed\n";
$usage .= "  -ccluster     : submit calibration to cluster and exit (requires --onlybuild)\n";
$usage .= "  -scluster <n> : split genome file into <n> pieces, submit <n> cmscan jobs and wait 3 minutes\n";
$usage .= "                  (changeable with -swait) before concatenating all the output files and continuing\n";
$usage .= "  -swait <n>    : with -scluster, set number of minutes to wait for cmscan jobs to finish to <n> [df: 3]\n";
$usage .= "\n OPTIONS USEFUL FOR DEVELOPMENT/DEBUGGING:\n";
$usage .= "  -skipfetch : use existing cmscan/hmmscan results, don't actually run it\n";
$usage .= "  -skipscan  : use existing cmscan/hmmscan results, don't actually run it\n";
$usage .= "  -skipaln   : use existing cmscan/hmmscan and alignment results, don't actually run it\n";

$usage .= "\n";

my ($seconds, $microseconds) = gettimeofday();
my $start_secs      = ($seconds + ($microseconds / 1000000.));
my $executable      = $0;
my $hmmbuild        = $hmmer_exec_dir  . "hmmbuild";
my $hmmpress        = $hmmer_exec_dir  . "hmmpress";
my $hmmalign        = $hmmer_exec_dir  . "hmmalign";
my $hmmfetch        = $hmmer_exec_dir  . "hmmfetch";
my $nhmmscan        = $hmmer_exec_dir  . "nhmmscan";
my $cmbuild         = $inf_exec_dir . "cmbuild";
my $cmcalibrate     = $inf_exec_dir . "cmcalibrate";
my $cmpress         = $inf_exec_dir . "cmpress";
my $cmscan          = $inf_exec_dir . "cmscan";
my $cmalign         = $inf_exec_dir . "cmalign";
my $cmfetch         = $inf_exec_dir . "cmfetch";

foreach my $x ($hmmbuild, $hmmpress, $nhmmscan, $cmbuild, $cmcalibrate, $cmpress, $cmscan, $esl_translate, $esl_ssplit) { 
  if(! -x $x) { die "ERROR executable file $x does not exist (or is not executable)"; }
}

# general options:
my $do_nocorrect = 0; # set to '1' if -nocorrect  enabled, do not look for internal start/stop codons and update annotation if found
my $do_matpept   = 0; # set to '1' if -matpept    enabled, genome has a single polyprotein, use mat_peptide info, not CDS
my $origin_seq   = undef; # defined if -oseq      enabled
my $do_strict    = 0; # set to '1' if -strict     enabled, matching annotations must be same index CDS+exon, else any will do
my $do_nodup     = 0; # set to '1' if -nodup      enabled, do not duplicate each genome, else do 
my $do_notexon   = 0; # set to '1' if -noexon     enabled, do not use exon-specific models, else do
my $do_onlybuild = 0; # set to '1' if -onlybuild  enabled, exit after building the model
my $in_model_db  = undef; # defined if -model <s> enabled, use <s> as the model file instead of building one
# options for controlling output table
my $do_concise   = 0; # set to '1' if -c       enabled, invoke concise output mode, set's all $do_no* variables below to '1'
my $do_nomdlb    = 0; # set to '1' if -nomdlb  or -c enabled, do not print HMM boundary info for annotations, else do
my $do_noexist   = 0; # set to '1' if -noexist or -c enabled, do not output information on existing annotations
my $do_nobrack   = 0; # set to '1' if -nobrack or -c enabled, do not output brackets around predicted annotations that do not match any existing annotation
my $do_nostop    = 0; # set to '1' if -nostop  or -c enabled, do not output stop codon for predicted annotations
my $do_nofid     = 0; # set to '1' if -nofid   or -c enabled, do not output fractional identities relative to the reference
my $do_noss3     = 0; # set to '1' if -noss3   or -c enabled, do not output SS3 columns: 'S'tart codon check, 'S'top codon check and multiple of '3' check
my $do_noolap    = 0; # set to '1' if -noolap  or -c enabled, do not output information on overlapping CDS/exons
my $do_noexp     = 0; # set to '1' if -noexp   or -c enabled, do not output explanatory information about column headings
# options for controlling homology search method
my $do_hmmer     = 0; # set to '1' if -hmmer      enabled, use HMMER3's nhmmscan, not Infernal 1.1
# options specific to HMMER3
my $do_hmmenv    = 0; # set to '1' if -hmmenv     enabled, use HMM envelope boundaries as predicted annotations, else use window boundaries
# options specific to Infernal
my $do_iglocal       = 0; # set to '1' if -iglocal    enabled, use -g with cmsearch
my $do_cslow         = 0; # set to '1' if -cslow      enabled, use default, slow, cmcalibrate parameters instead of speed optimized ones
my $do_ccluster      = 0; # set to '1' if -ccluster   enabled, submit cmcalibrate job to cluster
my $do_scluster      = 0; # set to '1' if -scluster   enabled, submit cmscan jobs to cluster
my $scluster_njobs   = undef; # set to a value if -scluster is used
my $df_scluster_wait = 3; # default value for number of minutes to wait for cmscan jobs to finish
my $scluster_wait    = $df_scluster_wait; # changeable to <n> with -swait <n>
# options for development/debugging
my $do_skipfetch = 0; # set to '1' if -skipfetch  enabled, skip genome fetch step, use pre-existing output from a previous run
my $do_skipscan  = 0; # set to '1' if -skipscan   enabled, skip cmscan/hmmscan step, use pre-existing output from a previous run
my $do_skipaln   = 0; # set to '1' if -skipaln    enabled, skip cmscan/hmmscan and alignment step, use pre-existing output from a previous run

&GetOptions("nocorrect" => \$do_nocorrect,
            "matpept"   => \$do_matpept,
            "oseq=s"    => \$origin_seq,
            "strict"    => \$do_strict,
            "nodup"     => \$do_nodup,
            "notexon"   => \$do_notexon,
            "onlybuild" => \$do_onlybuild,
            "model=s"   => \$in_model_db,
            "c"         => \$do_concise,
            "nomdlb"    => \$do_nomdlb,
            "noexist"   => \$do_noexist,
            "nobrack"   => \$do_nobrack,
            "nostop"    => \$do_nostop,
            "nofid"     => \$do_nofid,
            "noss3"     => \$do_noss3,
            "noolap"    => \$do_noolap,
            "noexp"     => \$do_noexp,
            "hmmer"     => \$do_hmmer,
            "hmmenv"    => \$do_hmmenv,
            "iglocal"   => \$do_iglocal,
            "cslow"     => \$do_cslow, 
            "ccluster"  => \$do_ccluster,
            "scluster=s"=> \$scluster_njobs,
            "swait=s"   => \$scluster_wait,
            "skipfetch" => \$do_skipfetch,
            "skipscan"  => \$do_skipscan,
            "skipaln"   => \$do_skipaln) ||
    die "Unknown option";

if(scalar(@ARGV) != 2) { die $usage; }
my ($dir, $listfile) = (@ARGV);
my $cmd; # a command to run with runCommand()

#$dir =~ s/\/*$//; # remove trailing '/' if there is one
#my $outdir     = $dir;
#my $outdirroot = $outdir;
#$outdirroot =~ s/^.+\///;

# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
if($do_nocorrect) { 
  $opts_used_short .= "-nocorrect ";
  $opts_used_long  .= "# option:  do not correct after identifying internal start/stops in predicted CDS [-nocorrect]\n";
}
if($do_matpept) { 
  $opts_used_short .= "-matpept ";
  $opts_used_long  .= "# option:  using mat_peptide info instead of CDS info [-matpept]\n";
}
if(defined $origin_seq) { 
  $opts_used_short .= "-oseq ";
  $opts_used_long  .= "# option:  searching for origin sequence of $origin_seq [-oseq]\n";
}
if($do_strict) { 
  $opts_used_short .= "-strict ";
  $opts_used_long  .= "# option:  demand matching annotations are same indexed CDS/exon [-strict]\n";
}
if($do_nodup) { 
  $opts_used_short .= "-nodup ";
  $opts_used_long  .= "# option:  not duplicating genomes, features that span the end..start will be undetectable [-nodup]\n";
}
if($do_notexon) { 
  $opts_used_short .= "-notexon ";
  $opts_used_long  .= "# option:  using full CDS, and not exon-specific models, for CDS with multiple exons [-noexon]\n";
}
if($do_onlybuild) { 
  $opts_used_short .= "-onlybuild ";
  $opts_used_long  .= "# option:  exit after model construction step [-onlybuild]\n";
}
if(defined $in_model_db) { 
  $opts_used_short .= "-model $in_model_db ";
  $opts_used_long  .= "# option:  use model in $in_model_db instead of building one here [-model]\n";
}
if($do_concise) { 
  $opts_used_short .= "-c ";
  $opts_used_long  .= "# option:  concise output mode [-c]\n";
}
if($do_nomdlb) { 
  $opts_used_short .= "-nomdlb ";
  $opts_used_long  .= "# option:  do not output HMM boundaries of predicted annotations [-nomdlb]\n";
}
if($do_noexist) { 
  $opts_used_short .= "-noexist";
  $opts_used_long  .= "# option:  not outputting info on existing annotations [-noexist]\n";
}
if($do_nobrack) { 
  $opts_used_short .= "-nobrack";
  $opts_used_long  .= "# option:  not putting brackets around predicted start/stop positions [-nobrack]\n";
}
if($do_nostop) { 
  $opts_used_short .= "-nostop";
  $opts_used_long  .= "# option:  do not output stop codons [-nostop]\n";
}
if($do_nofid) { 
  $opts_used_short .= "-nofid";
  $opts_used_long  .= "# option:  do not output fractional identities [-nofid]\n";
}
if($do_noss3) { 
  $opts_used_short .= "-noss3";
  $opts_used_long  .= "# option:  do not output start codon, stop codon or multiple of 3 test results [-noss3]\n";
}
if($do_noolap) { 
  $opts_used_short .= "-noolap";
  $opts_used_long  .= "# option:  do not output information on overlaps [-noolap]\n";
}
if($do_noexp) { 
  $opts_used_short .= "-noexp";
  $opts_used_long  .= "# option:  do not output information on column headings [-noexp]\n";
}
if($do_hmmer) { 
  $opts_used_short .= "-hmmer";
  $opts_used_long  .= "# option:  using HMMER3 for predicting annotation [-hmmer]\n";
}
if($do_hmmenv) { 
  $opts_used_short .= "-hmmenv ";
  $opts_used_long  .= "# option:  use HMM envelope boundaries as predicted annotations, not window boundaries [-hmmenv]\n";
}
if($do_iglocal) { 
  $opts_used_short .= "-iglocal ";
  $opts_used_long  .= "# option:  use glocal search option with Infernal [-iglocal]\n";
}
if($do_cslow) { 
  $opts_used_short .= "-cslow ";
  $opts_used_long  .= "# option:  run cmcalibrate in default (slow) mode [-cslow]\n";
}
if($do_ccluster) { 
  $opts_used_short .= "-ccluster ";
  $opts_used_long  .= "# option:  submit calibration job to cluster [-ccluster]\n";
}
if(defined $scluster_njobs) { 
  $do_scluster = 1;
  $opts_used_short .= "-scluster ";
  $opts_used_long  .= "# option:  submit cmscan jobs to cluster and wait for them to finish [-scluster]\n";
}
if($scluster_wait != $df_scluster_wait) { 
  $opts_used_short .= "-swait ";
  $opts_used_long  .= "# option:  submit cmscan jobs to cluster and wait for them to finish [-scluster]\n";
}
if($do_skipfetch) { 
  $opts_used_short .= "-skipfetch ";
  $opts_used_long  .= "# option:  use existing fetched genome file, don't create a new one [-skipfetch]\n";
}
if($do_skipscan) { 
  $opts_used_short .= "-skipscan ";
  $opts_used_long  .= "# option:  use existing cmscan/hmmscan results, don't actually run it [-skipscan]\n";
}
if($do_skipaln) { 
  $opts_used_short .= "-skipaln ";
  $opts_used_long  .= "# option:  use existing cmscan/hmmscan and alignment results, don't actually run either [-skipaln]\n";
}
# 
# check for incompatible option values/combinations:
if((! $do_hmmer) && $do_hmmenv) { 
  die "ERROR -hmmenv requires -hmmer"; 
}
if($do_ccluster && $do_hmmer) { 
  die "ERROR -ccluster is incompatible with -hmmer"; 
}
if($do_cslow && $do_hmmer) { 
  die "ERROR -cslow is incompatible with -hmmer"; 
}
if($do_skipscan && ($do_ccluster || $do_cslow || $do_onlybuild)) { 
  die "ERROR -skipscan is incompatible with -ccluster, -cslow, and -onlybuild";
}
if($do_scluster && ($do_skipaln || $do_skipscan)) { 
  die "ERROR -skipscan and -skipaln are incompatible with -scluster";
}
if($do_scluster && $do_hmmer) { 
  die "ERROR -scluster and -hmmer are incompatible";
}

# check that options that must occur in combination, do
if($do_ccluster && (! $do_onlybuild)) { 
  die "ERROR -ccluster must be used in combination with -onlybuild"; 
}
if(($scluster_wait != $df_scluster_wait) && (! $do_scluster)) { 
  die "ERROR -swait must be used in combination with -scluster"; 
}

# check that input files related to options actually exist
if(defined $in_model_db) { 
  if(! -s $in_model_db) { die "ERROR: $in_model_db file does not exist"; }
}
# verify origin sequence if necessary
my $origin_offset = undef;
if(defined $origin_seq) { 
  $origin_offset = validateOriginSeq($origin_seq);
  $origin_seq =~ s/\|//;
}

# if in $concise output mode, turn on other affected options:
if($do_concise) { 
  $do_nomdlb =  1;
  $do_noexist = 1;
  $do_nobrack = 1;
  $do_nostop  = 1;
  $do_nofid   = 1;
  $do_noss3   = 1;
  $do_noolap  = 1;
}

###############
# Preliminaries
###############
# check if the $dir exists, and that it contains the files we need
if(! -d $dir)      { die "ERROR directory $dir does not exist"; }
if(! -s $listfile) { die "ERROR list file $listfile does not exist, or is empty"; }
my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $gene_tbl_file    = $dir . "/" . $dir_tail . ".gene.tbl";
my $cds_tbl_file     = $dir . "/" . $dir_tail . ".CDS.tbl";
my $matpept_tbl_file = $dir . "/" . $dir_tail . ".mat_peptide.tbl";
my $length_file      = $dir . "/" . $dir_tail . ".length";
my $out_root = $dir . "/" . $dir_tail;
#if(! -s $gene_tbl_file) { die "ERROR $gene_tbl_file does not exist."; }
if($do_matpept) { if(! -s $matpept_tbl_file)  { die "ERROR $matpept_tbl_file does not exist."; } }
else            { if(! -s $cds_tbl_file)      { die "ERROR $cds_tbl_file does not exist."; } }
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
printf("#\n");

#####################
# parse the list file
#####################
my @accn_A = (); # array of accessions
open(IN, $listfile) || die "ERROR unable to open $listfile for reading"; 
my $waccn = 0; # max length of all accessions
while(my $accn = <IN>) { 
  if($accn =~ m/\w/) { 
    chomp $accn;
    stripVersion(\$accn); # remove version
    push(@accn_A, $accn);
    if(length($accn) > $waccn) { $waccn = length($accn); }
  }
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

# parseTable($gene_tbl_file, \%gene_tbl_HHA);
if($do_matpept) { 
  parseTable($matpept_tbl_file, \%cds_tbl_HHA);
}
else {
  parseTable($cds_tbl_file, \%cds_tbl_HHA);
}

#######################
# variable declarations
#######################
my $strand_str;              # +/- string for all CDS for an accession: e.g. '+-+': 1st and 3rd CDS are + strand, 2nd is -

# reference information on reference accession, first accession read in ntlist file
my $ref_accn          = undef; # changed to <s> with -ref <s>
my $ref_label_str     = undef; # label string for reference accn
my $ref_ncds          = 0;     # number of CDS in reference
my $ref_strand_str    = "";    # strand string for reference 
my @ref_cds_len_A     = ();    # [0..$i..$ref_ncds-1]: length of each reference CDS
#my @ref_cds_len_tol_A = ();   # [0..$i..$ref_ncds-1]: length tolerance, any gene that is within this fraction of the lenght of the ref gene is a match
my @ref_cds_coords_A  = ();    # [0..$i..$ref_ncds-1]: CDS coords for reference
my @ref_cds_product_A = ();    # CDS:product qualifier data for reference 

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

#####################################################
# Fetch all genome sequences, including the reference
#####################################################
my $naccn = scalar(@accn_A);
my $gnm_fetch_file = $out_root . ".fg.idfetch.in";
my $gnm_fasta_file = $out_root . ".fg.fa";
my @seq_accn_A = (); # [0..$naccn-1] name of genome fasta sequence for each accn
my $seq_accn;     # temp fasta sequence name
my $fetch_string = undef;
my $ref_seq_accn; # name of fasta sequence for reference

# another error check
if($do_scluster && ($scluster_njobs > $naccn)) { 
  die "ERROR with -scluster <n>, <n> must be <= number of genomes ($naccn), but you used $scluster_njobs"; 
}

open(OUT, ">" . $gnm_fetch_file) || die "ERROR unable to open $gnm_fetch_file";
for(my $a = 0; $a < $naccn; $a++) { 
#  print OUT $accn_A[$a] . "\n";
  my $accn = $accn_A[$a];
  if(! exists $totlen_H{$accn}) { die "ERROR no total length read for accession $accn"; } 
  if($do_nodup) { 
    $fetch_string = $accn . ":1.." . $totlen_H{$accn} . "\n";
    print OUT $accn . ":" . "genome" . "\t" . $fetch_string;
    $seq_accn = $accn . ":genome:" . $accn . ":1:" . $totlen_H{$accn} . ":+:";
  }
  else { 
    $fetch_string = "join(" . $accn . ":1.." . $totlen_H{$accn} . "," . $accn . ":1.." . $totlen_H{$accn} . ")\n";
    print OUT $accn . ":" . "genome-duplicated" . "\t" . $fetch_string;
    $seq_accn = $accn . ":genome-duplicated:" . $accn . ":1:" . $totlen_H{$accn} . ":+:" . $accn . ":1:" . $totlen_H{$accn} . ":+:";
  }
  push(@seq_accn_A, $seq_accn);
  if($a == 0) { $ref_seq_accn = $seq_accn; }
}
close(OUT);

if($do_skipfetch) { 
  # trying to skip fetch step, make sure we have the fasta file
  printf("%-65s ... ", "# Skipping genome fetch step");
  if(! -s $gnm_fasta_file) { die "ERROR, with -skipfetch genome fasta file is expected to exist, but $gnm_fasta_file does not or is empty"; }
  printf("done. [-skipfetch]\n");
}
else { 
  # remove the file we're about to create ($gnm_fasta_file), and any .ssi index that may exist with it
  if(-e $gnm_fasta_file)          { unlink $gnm_fasta_file; }
  if(-e $gnm_fasta_file . ".ssi") { unlink $gnm_fasta_file . ".ssi"; }

  printf("%-65s ... ", sprintf("# Fetching $naccn full%s genome sequences", $do_nodup ? "" : " (duplicated)"));
# my $cmd = "$idfetch -t 5 -c 1 -G $gnm_fetch_file > $gnm_fasta_file";
  $cmd = "perl $esl_fetch_cds -nocodon $gnm_fetch_file > $gnm_fasta_file";
  my $secs_elapsed = runCommand($cmd, 0);
  printf("done. [%.1f seconds]\n", $secs_elapsed);
}

# and open the sequence file using BioEasel
my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $gnm_fasta_file });

##################################################
# If we're looking for an origin sequence, do that
##################################################
my %origin_coords_HA = ();
if(defined $origin_seq) {
  findSeqInFile($sqfile, $origin_seq, $do_nodup, \%origin_coords_HA);
}

########################################################
# Gather information and sequence data on the reference.
# Use each reference CDS and reference CDS exon as a 
# homology search model against all the genomes.
#######################################################
$ref_accn = $accn_A[0];
if(! exists ($cds_tbl_HHA{$ref_accn})) { die "ERROR no CDS information stored for reference accession"; }
(undef, undef, undef, undef, undef, $ref_strand_str) = getStrandStats(\%cds_tbl_HHA, $ref_accn);
getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $ref_accn, \@ref_cds_len_A, \@ref_cds_coords_A);
getQualifierValues(\%cds_tbl_HHA, $ref_accn, "product", \@ref_cds_product_A);
$ref_ncds = scalar(@ref_cds_len_A);

my $all_stk_file = $out_root . ".ref.all.stk";

($seconds, $microseconds) = gettimeofday();
my $start_time = ($seconds + ($microseconds / 1000000.));
printf("%-65s ... ", "# Fetching reference CDS sequences");
my $cur_out_root;
my $cur_name_root;
my $fetch_input;
my $fetch_output;
my $nhmm = 0;               # number of HMMs (and alignments used to build those HMMs)
my @hmm2cds_map_A = ();     # [0..h..$nhmm-1]: $i: HMM ($h+1) maps to reference cds ($i+1)
my @hmm2exon_map_A = ();    # [0..h..$nhmm-1]: $e: HMM ($h+1) maps to exon ($e+1) of reference cds $hmm2cds_map_A[$h]+1
my @hmm_is_first_A = ();    # [0..h..$nhmm-1]: '1' if HMM ($h+1) is the first one for cds $hmm2cds_map_A[$h], else 0
my @hmm_is_final_A = ();    # [0..h..$nhmm-1]: '1' if HMM ($h+1) is the final one for cds $hmm2cds_map_A[$h], else 0
my @cds2first_hmm_A = ();   # [0..$c..ncds-1]: $h, first exon of CDS $c+1 is modeled by HMM $h+1
my @cds2final_hmm_A = ();   # [0..$c..ncds-1]: $h, final exon of CDS $c+1 is modeled by HMM $h+1
my @model_A = ();           # [0..$nhmm-1]: array of model HMM names, also name of stockholm alignments used to build those HMMs
my @cds_out_short_A   = (); # [0..$ref_ncds-1]: array of abbreviated model CDS names to print
my @cds_out_product_A = (); # [0..$ref_ncds-1]: array of 'CDS:product' qualifier (protein names)
my %mdllen_H          = (); # key: model name from @model_A, value is model length
my @ref_nexons_A      = (); # [0..$c..$ref_ncds-1]: number of exons in CDS $c+1
my $ref_tot_nexons    = 0;  # total number of exons in all CDS

# for each reference CDS, fetch each exon (or the full CDS if -notexon enabled)
for(my $i = 0; $i < $ref_ncds; $i++) { 
  # printf("REF CDS $i $ref_cds_product_A[$i]\n");
  
  # determine start and stop positions of all exons
  my @starts_A = ();
  my @stops_A  = ();
  my $nexons   = 0;
  startStopsFromCoords($ref_cds_coords_A[$i], \@starts_A, \@stops_A, \$nexons);
  push(@ref_nexons_A, $nexons);
  $ref_tot_nexons += $nexons;

  # if we're on the negative strand, reverse the arrays, they'll be in the incorrect order
  my $strand = substr($ref_strand_str, $i, 1);
  if($strand eq "-") { 
    @starts_A = reverse @starts_A;
    @stops_A  = reverse @stops_A;
  }

  # are we going to fetch multiple exons?
  my $cur_multi_exon = ($nexons == 1 || $do_notexon) ? 0 : 1;
  my $act_nexons = $nexons;
  if(! $cur_multi_exon) { $nexons = 1; }

  # for each exon, note that if $do_notexon is true, $nexons was redefined as 1 above
  for(my $e = 0; $e < $nexons; $e++) { 
    if($cur_multi_exon) { 
      $cur_out_root  = $out_root . ".ref.cds." . ($i+1) . ".exon." . ($e+1);
      $cur_name_root = $dir_tail . ".ref.cds." . ($i+1) . ".exon." . ($e+1);
    }
    else { 
      $cur_out_root  = $out_root . ".ref.cds." . ($i+1);
      $cur_name_root = $dir_tail . ".ref.cds." . ($i+1);
    }
    
    # determine start and stop of the region we are going to fetch
    my $start = ($cur_multi_exon) ? $starts_A[$e] : $starts_A[0];
    my $stop  = ($cur_multi_exon) ? $stops_A[$e]  : $stops_A[$nexons-1]; 
    if($strand eq "-") { # swap start and stop
      my $tmp = $start;
      $start = $stop;
      $stop  = $tmp;
    }
    my @fetch_AA = ();
    push(@fetch_AA, [$cur_name_root, $start, $stop, $ref_seq_accn]);
    
    # fetch the sequence
    my $cur_fafile = $cur_out_root . ".fa";
    $sqfile->fetch_subseqs(\@fetch_AA, undef, $cur_fafile);
    
    # reformat to stockholm
    my $cur_stkfile = $cur_out_root . ".stk";
    my $cmd = "esl-reformat --informat afa stockholm $cur_fafile > $cur_stkfile";
    runCommand($cmd, 0);
    
    # annotate the stockholm file with a blank SS and with a name
    my $do_blank_ss = (! $do_hmmer); # add a blank SS_cons line if we're using Infernal
    my $cur_named_stkfile = $cur_out_root . ".named.stk";
    my $mdllen = annotateStockholmAlignment($cur_name_root, $do_blank_ss, $cur_stkfile, $cur_named_stkfile);

    # store information on this model's name for output purposes
    if($e == ($nexons-1)) { 
      my $short = sprintf("CDS #%d", ($i+1));
      if($act_nexons > 1) { $short .= " [$act_nexons exons; $strand]"; }
      else                { $short .= " [single exon; $strand]"; }
      push(@cds_out_short_A,   $short);
      push(@cds_out_product_A, $ref_cds_product_A[$i]);
    }
    push(@model_A, $cur_name_root);

    $mdllen_H{$cur_name_root} = $mdllen;

    # now append the named alignment to the growing stockholm alignment database $all-stk_file
    $cmd = "cat $cur_named_stkfile";
    if($nhmm == 0) { $cmd .= " >  $all_stk_file"; }
    else           { $cmd .= " >> $all_stk_file"; }
    runCommand($cmd, 0);
    push(@hmm2cds_map_A,  $i);
    push(@hmm2exon_map_A, $e);
    push(@hmm_is_first_A, ($e == 0)           ? 1 : 0);
    push(@hmm_is_final_A, ($e == ($nexons-1)) ? 1 : 0);
    if($e == 0)           { $cds2first_hmm_A[$i] = $nhmm; }
    if($e == ($nexons-1)) { $cds2final_hmm_A[$i] = $nhmm; }
    $nhmm++;
  }
}
($seconds, $microseconds) = gettimeofday();
my $stop_time = ($seconds + ($microseconds / 1000000.));
printf("done. [%.1f seconds]\n", ($stop_time - $start_time));

######################
# HOMOLOGY SEARCH STEP
######################
my $model_db; # model database file, either HMMs or CMs

# first, create the model database, unless it was passed in:
if(defined $in_model_db) { 
  $model_db = $in_model_db;
}
else { 
  if(! $do_hmmer) { 
    createCmDb($cmbuild, $cmcalibrate, $cmpress, $nhmm, $do_cslow, $do_ccluster, $all_stk_file, $out_root . ".ref");
    if($do_onlybuild) { 
      printf("#\n# Model calibration %s. Exiting.\n", ($do_ccluster) ? "job submitted" : "complete");
      exit 0;
    }
    $model_db = $out_root . ".ref.cm";
  }
  else { # use HMMER3's nhmmscan
    createHmmDb($hmmbuild, $hmmpress, $nhmm, $all_stk_file, $out_root . ".ref");
    if($do_onlybuild) { 
      printf("#\n# Model construction complete. Exiting.\n");
      exit 0;
    }
    $model_db = $out_root . ".ref.hmm";
  }
}
  
# now we know we have a model database, perform the search and parse the results
# output files from homology searches
my $tblout = $out_root . ".tblout"; # tabular output file, created by either nhmmscan or cmsearch
my $stdout = $out_root . ".stdout"; # standard output file from either nhmmscan or cmsearch

# 2D hashes for storing the search results
# For all of these: 1D key: model name, 2D key: sequence name,
# and the value is for the top-scoring hit only
my %p_start_HH    = (); # start positions of hits
my %p_stop_HH     = (); # stop positions of hits
my %p_strand_HH   = (); # strands of hits
my %p_score_HH    = (); # bit score of hits 
my %p_hangover_HH = (); # "<a>:<b>" where <a> is number of 5' model positions not in the alignment
                        # and <b> is number of 3' model positions not in the alignment
my %p_fid2ref_HH  = (); # fractional identity to reference 
my %p_refdel_HHA  = (); # array of reference positions that are deleted for the alignment of this sequence
my %p_refins_HHA  = (); # array of strings ("<rfpos>:<count>") reference positions that are deleted for the alignment of this sequence

my %pred_fafile_H = (); # hash of names of fasta files for predicted exon sequences, keys: model name from @model_A, value: name of file

if($do_hmmer) { 
  runNhmmscan($nhmmscan, ($do_skipscan || $do_skipaln), $model_db, $gnm_fasta_file, $tblout, $stdout);
  parseNhmmscanTblout($tblout, $do_hmmenv, \%totlen_H, \%p_start_HH, \%p_stop_HH, \%p_strand_HH, \%p_score_HH, \%p_hangover_HH);
  fetchHits($sqfile, $do_skipaln, "predicted", \@model_A, \@seq_accn_A, \%totlen_H, \%p_start_HH, \%p_stop_HH, \%p_strand_HH, $out_root, \%pred_fafile_H);
}
else { 
  if(! $do_scluster) { 
    # default method, run cmscan on full genome file
    runCmscan($cmscan, $do_iglocal, ($do_skipscan || $do_skipaln), 0, $model_db, $gnm_fasta_file, $tblout, $stdout);
  }
  else { # we need to split up the genome fasta file and submit a different cmscan job for each fasta file
    splitFastaFile($esl_ssplit, $gnm_fasta_file, $scluster_njobs);
    # now submit a job for each
    printf("%-65s ... ", sprintf("# Submitting %d cmscan jobs", $scluster_njobs));
    for(my $z = 1; $z <= $scluster_njobs; $z++) { 
      my $cur_gnm_fasta_file = $gnm_fasta_file . "." . $z;
      my $cur_tblout = $tblout . "." . $z;
      my $cur_stdout = $stdout . "." . $z;
      runCmscan($cmscan, $do_iglocal, 0, 1, $model_db, $cur_gnm_fasta_file, $cur_tblout, $cur_stdout);
    }
    printf("done. [waiting max of $scluster_wait minutes for them to finish]\n");

    # now check that they're all finished every 15 seconds until we hit our max wait time
    my $all_finished = 0;
    my $nfinished = 0;
    for(my $zz = 0; $zz < ($scluster_wait*4); $zz++) { 
      # check to see if jobs are finished, every 30 seconds
      sleep(15);
      $nfinished = 0;
      for(my $z = 1; $z <= $scluster_njobs; $z++) { 
        my $cur_stdout = $stdout . "." . $z;
        if(-s $cur_stdout) { 
          my $final_line = `tail -n 1 $cur_stdout`;
          chomp $final_line;
          if($final_line eq "[ok]") { 
            $nfinished++;
          }
        }
      }
      if($nfinished == $scluster_njobs) { 
        # we're done break out of it
        $all_finished = 1;
        $zz = ($scluster_wait * 4);
      }
    }
    if($all_finished) { 
      # concatenate all outputs into one main one
      for(my $z = 1; $z <= $scluster_njobs; $z++) { 
        my $cur_gnm_fasta_file = $gnm_fasta_file . "." . $z;
        my $cur_stdout = $stdout . "." . $z;
        my $cur_tblout = $tblout . "." . $z;
        my $output_char = ($z == 1) ? ">" : ">>";
        $cmd = "cat $cur_stdout $output_char $stdout; rm $cur_stdout; rm $cur_gnm_fasta_file";
        runCommand($cmd, 0);
        $cmd = "cat $cur_tblout $output_char $tblout; rm $cur_tblout";
        runCommand($cmd, 0);
      }
    }
    else { 
      die "ERROR only $nfinished of the $scluster_njobs are finished after $scluster_wait minutes. Increase limit with -swait";
    }
  } # end of 'else' entered if we're submitting jobs to the cluster

  parseCmscanTblout($tblout, \%totlen_H, \%mdllen_H, \%p_start_HH, \%p_stop_HH, \%p_strand_HH, \%p_score_HH, \%p_hangover_HH);
  fetchHits($sqfile, $do_skipaln, "predicted", \@model_A, \@seq_accn_A, \%totlen_H, \%p_start_HH, \%p_stop_HH, \%p_strand_HH, $out_root, \%pred_fafile_H);
}

########################################################
# COMBINE MULTI-EXON SEQUENCES INTO SINGLE CDS SEQUENCES
########################################################
my @pred_cds_fafile_A = (); # array of predicted CDS sequence files, filled below
wrapperCombineExonsIntoCDS($nhmm, $dir, "predicted", \@model_A, \@hmm2cds_map_A, \@hmm_is_first_A, \@hmm_is_final_A, \@pred_cds_fafile_A);

###########################################################################
# CORRECT PREDICTIONS
# - look for internal starts and stops and make corrections based on those
# - combine new exons into new CDS
###########################################################################
($seconds, $microseconds) = gettimeofday();
my $start_time = ($seconds + ($microseconds / 1000000.));
printf("%-65s ... ", "# Translating predicted CDS to identify internal starts/stops");

my $source_accn;   # name of CDS sequence that was translated
my $source_length; # length of CDS sequence that was translated
my $coords_from;   # start nt coordinate of current translation
my $coords_to;     # stop nt coordinate of current translation
my @corr_cds_start_AH = ();  # [0..$i..ncds-1], each element is a hash with keys $key as sequence accessions and values 
                         # are number of nucleotides that the prediction of the start coordinate should be corrected
                         # based on an esl-translate translation of the predicted CDS sequence, values can be negative
                         # or positive
my @corr_cds_stop_AH = ();   # [0..$i..ncds-1], each element is a hash with keys $key as sequence accessions and values 
                         # are number of nucleotides that the prediction of the stop coordinate should be corrected
                         # based on an esl-translate translation of the predicted CDS sequence, values can be negative
                         # or positive
my @coords_len_AH = ();  # [0..$i..ncds-1], each element is a hash with keys $key as sequence accessions and values 
                         # are lengths of translated protein sequences (in nucleotides corresponding to values in $corr_cds_start_AH
                         # and $corr_cds_stop_AH
my @did_corr_exon_start_AH = ();  # [0..$i..ref_nexons-1], each element is a hash with keys $key as sequence accessions and values 
                                  # are '1' if this exon's start position was corrected
my @did_corr_exon_stop_AH = ();   # [0..$i..ref_nexons-1], each element is a hash with keys $key as sequence accessions and values 
                                  # are '1' if this exon's stop position was corrected
my %c_start_HH = (); # corrected start positions of hits, start with a copy of p_start_HH
my %c_stop_HH  = ();  # corrected stop positions of hits,  start with a copy of p_stop_HH
my %corr_fafile_H = ();     # hash of names of fasta files for corrected exon sequences, keys: model name from @model_A, value: name of file
my @corr_cds_fafile_A = (); # array of corrected CDS sequence files, filled below

if(! $do_nocorrect) { 
  # translate predicted CDS sequences using esl-translate to see if any corrections to predictions are necessary
  for(my $c = 0; $c < $ref_ncds; $c++) { 
    # determine first and final HMM for this 

    my $cur_fafile = $pred_cds_fafile_A[$c];
    # translate into AA sequences
    my $tmp_esl_translate_output = $cur_fafile;
    $tmp_esl_translate_output    =~ s/\.cds/\.esl-translate/;
    # --watson specifies only translate top strand (not reverse strand)
    # -m specifies only AUG allowed as start
    $cmd = $esl_translate . " --watson -m $cur_fafile | grep ^\\> | grep \"frame=1\" > $tmp_esl_translate_output";
    runCommand($cmd, 0);

    # now based on the esl-translate deflines, determine if there are any internal starts or stops
    # note that we've grep'ed for only frame=1 translations, and we've specified only top strand with
    # --watson, so we should only see at most 1 translation of each sequence
    open(IN, $tmp_esl_translate_output) || die "ERROR unable to open $tmp_esl_translate_output for reading";
    while(my $line = <IN>) { 
      parseEslTranslateDefline($line, \$coords_from, \$coords_to, \$source_accn, undef, \$source_length);
      # keep longest translated stretch
      my $coords_len = ($coords_to - $coords_from + 1);
      if((! exists $corr_cds_start_AH[$c]{$source_accn}) || ($coords_len > $coords_len_AH[$c]{$source_accn})) { 
        $corr_cds_start_AH[$c]{$source_accn} = $coords_from - 1;
        $corr_cds_stop_AH[$c]{$source_accn}  = -1 * (($source_length - 3) - $coords_to); # source length includes stop codon
        $coords_len_AH[$c]{$source_accn} = $coords_len;
      }
    }
    close(IN);
  }

  # now we have all corrections, go back and create c_start, c_stop data structures based 
  # on p_start, p_stop data structures and these corrections, this is tricky for multi-exon
  # CDS because we have to make sure we correct the proper exon
  for(my $a = 0; $a < $naccn; $a++) { 
    my $accn     = $accn_A[$a];
    my $seq_accn = $seq_accn_A[$a];
    my $mdl; # name of a model
    for(my $c = 0; $c < $ref_ncds; $c++) { 
      if(! exists $corr_cds_start_AH[$c]{$accn}) { 
        die sprintf("ERROR no corrected start position value for CDS %d sequence $accn", $c+1); 
      }
      if(! exists $corr_cds_stop_AH[$c]{$accn}) { 
        die sprintf("ERROR no corrected stop position value for CDS %d sequence $accn", $c+1); 
      }
      my $corr_start = $corr_cds_start_AH[$c]{$accn};
      my $corr_stop  = $corr_cds_stop_AH[$c]{$accn};
      # sanity check
      if($corr_start < 0) { die "ERROR corrected start less than 0, can't deal with that yet";  }
      if($corr_stop  > 0) { die "ERROR corrected stop greather than 0, can't deal with that yet"; }

      my $first_hmm = $cds2first_hmm_A[$c];
      my $final_hmm = $cds2final_hmm_A[$c];
      my $cur_nhmm = $final_hmm - $first_hmm + 1;

      if($first_hmm == $final_hmm) {
        # easy case: single exon CDS:
        my $mdl = $model_A[$first_hmm];
        if($p_strand_HH{$mdl}{$seq_accn} eq "+") { 
          $c_start_HH{$mdl}{$seq_accn} = $p_start_HH{$mdl}{$seq_accn} + $corr_start; # note that $corr_start may be 0
          $c_stop_HH{$mdl}{$seq_accn}  = $p_stop_HH{$mdl}{$seq_accn}  + $corr_stop;  # note that $corr_stop  may be 0
        }
        else { 
          $c_start_HH{$mdl}{$seq_accn} = $p_start_HH{$mdl}{$seq_accn} - $corr_start; # note that $corr_start may be 0
          $c_stop_HH{$mdl}{$seq_accn}  = $p_stop_HH{$mdl}{$seq_accn}  - $corr_stop;  # note that $corr_stop  may be 0
        }
      }
      else { # multi exon CDS
        # get temporary array of exon lengths and determine full CDS length
        my $cds_len    = 0;
        my @exon_len_A = ();
        for(my $h = $first_hmm; $h <= $final_hmm; $h++) { 
          my $hp = $h - $first_hmm;
          $exon_len_A[$hp] = ($p_strand_HH{$model_A[$h]}{$seq_accn} eq "+") ? 
              ($p_stop_HH{$model_A[$h]}{$seq_accn}  - $p_start_HH{$model_A[$h]}{$seq_accn} + 1) : 
              ($p_start_HH{$model_A[$h]}{$seq_accn} - $p_stop_HH{$model_A[$h]}{$seq_accn}  + 1);
          $cds_len += $exon_len_A[$hp];
          
          # IMPORTANT: copy predicted to corrected, we may 'correct' this below, but if we don't the prediction stays the same.
          #            In case when we correct a multi-exon gene by putting the stop in the non-final exon or the start in
          #            non-first exon, note that this means we will only 'correct' the exons that have the modified stop
          #            or start, and NOT any others. For example if we have a 2-exon CDS and our correction puts the stop
          #            in the first exon, we will not change the annotation of the second exon.
          $mdl = $model_A[$h];
          $c_start_HH{$mdl}{$seq_accn} = $p_start_HH{$mdl}{$seq_accn};
          $c_stop_HH{$mdl}{$seq_accn}  = $p_stop_HH{$mdl}{$seq_accn};
        }        

        my $len_so_far; # total length of all exons so far
        my $found_exon; # set to '1' when we find the correct exon to put the corrected start in
        if($corr_start > 0) { # determine which exon the corrected start is in, and correct it:
          $len_so_far = 0; 
          $found_exon = 0; 
          for(my $h = $first_hmm; $h <= $final_hmm; $h++) {
            if(! $found_exon) { 
              my $hp = $h - $first_hmm;
              $mdl = $model_A[$h];
              $len_so_far += $exon_len_A[$hp];
              if($corr_start <= $len_so_far) { 
                if($p_strand_HH{$mdl}{$seq_accn} eq "+") { 
                  $c_start_HH{$mdl}{$seq_accn} += ($corr_start - ($len_so_far - $exon_len_A[$hp]));
                }
                else { 
                  $c_start_HH{$mdl}{$seq_accn} -= ($corr_start - ($len_so_far - $exon_len_A[$hp]));
                }
                $found_exon = 1;
                $did_corr_exon_start_AH[$h]{$accn} = 1;
                #printf("corr_start: $corr_start\n");
                #printf("len_so_far: $len_so_far\n");
                #printf("exon_len_A[$hp]: $exon_len_A[$hp]\n");
                #printf("p_start_HH{$mdl}{$seq_accn} exon %d: $p_start_HH{$mdl}{$seq_accn}\n", ($h-$first_hmm+1));
                #printf("c_start_HH{$mdl}{$seq_accn} exon %d: $c_start_HH{$mdl}{$seq_accn}\n", ($h-$first_hmm+1));
              }
            }
          }
          if(! $found_exon) { die sprintf("ERROR unable to find proper exon for corrected start for CDS: %d, accn: $accn\n", $c+1); }
        } # end of 'if($corr_start > 0)'
        
        if($corr_stop < 0) { # determine which exon the corrected stop is in, and correct it:
          $len_so_far = 0; 
          $found_exon = 0; 
          for(my $h = $final_hmm; $h >= $first_hmm; $h--) { 
            if(! $found_exon) {
              my $hp = $h - $first_hmm;
              $mdl = $model_A[$h];
              $len_so_far += $exon_len_A[$hp];
              if((-1 * $corr_stop) <= $len_so_far) { 
                if($p_strand_HH{$mdl}{$seq_accn} eq "+") { 
                  $c_stop_HH{$mdl}{$seq_accn}  -= ((-1 * $corr_stop) - ($len_so_far - $exon_len_A[$hp]));
                }
                else { 
                  $c_stop_HH{$mdl}{$seq_accn}  += ((-1 * $corr_stop) - ($len_so_far - $exon_len_A[$hp]));
                }
                $found_exon = 1;
                $did_corr_exon_stop_AH[$h]{$accn} = 1;
                #printf("corr_stop: $corr_stop\n");
                #printf("len_so_far: $len_so_far\n");
                #printf("exon_len_A[$hp]: $exon_len_A[$hp]\n");
                #printf("p_stop_HH{$mdl}{$seq_accn} exon %d: $p_stop_HH{$mdl}{$seq_accn}\n", ($h-$first_hmm+1));
                #printf("c_stop_HH{$mdl}{$seq_accn} exon %d: $c_stop_HH{$mdl}{$seq_accn}\n", ($h-$first_hmm+1));
              }
            }
          }
          if(! $found_exon) { die sprintf("ERROR unable to find proper exon for corrected stop for CDS: %d, accn: $accn\n", $c+1); }
        } # end of 'if($corr_stop < 0)'
      } # end of 'else' entered if we're a multi-exon CDS
    } # end of loop over CDS
  } # end of loop over sequences
  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));
  printf("done. [%.1f seconds]\n", ($stop_time - $start_time));

  # fetch corrected hits into new files
  fetchHits($sqfile, $do_skipaln, "corrected", \@model_A, \@seq_accn_A, \%totlen_H, \%c_start_HH, \%c_stop_HH, \%p_strand_HH, $out_root, \%corr_fafile_H);

  # combine multi-exon sequences into CDS:
  wrapperCombineExonsIntoCDS($nhmm, $dir, "corrected", \@model_A, \@hmm2cds_map_A, \@hmm_is_first_A, \@hmm_is_final_A, \@corr_cds_fafile_A);
} # end of if(! $do_nocorrect)

#########################################
# TRANSLATE PREDICTIONS INTO PROTEIN SEQS
#########################################
my @nfullprot_A   = ();  # [0..$i..ncds-1], number of accessions we have a full protein for, for CDS $i
my @fullprot_AH   = ();  # [0..$i..ncds-1], each element is a hash with keys $key as sequence accessions and values 
                         # of number of full length protein sequences we have for CDS $i for key $key. Values should
                         # always be '1', more than '1' is an error, and we never create a value of '0'.
my @ntruncprot_A   = (); # [0..$i..ncds-1], number of accessions we do not have a full protein for, for CDS $i
my @truncprot_AH   = (); # [0..$i..ncds-1], each element is a hash with keys $key as sequence accessions and values 
                           # of number of non-full length protein sequences we have for CDS $i for key $key. Values should
                           # always be '1', more than '1' is an error, and we never create a value of '0'.
for(my $c = 0; $c < $ref_ncds; $c++) { 
  my $cur_fafile = ($do_nocorrect) ? $pred_cds_fafile_A[$c] : $corr_cds_fafile_A[$c];
  # translate into AA sequences
  my $tmp_aa_fafile   = $cur_fafile;
  my $aa_full_fafile  = $cur_fafile;
  my $aa_trunc_fafile = $cur_fafile;
  $tmp_aa_fafile      =~ s/\.cds/\.aa.tmp/;
  $aa_full_fafile     =~ s/\.cds/\.aa.full/;
  # $aa_trunc_fafile    =~ s/\.cds/\.aa.trunc/;
  
  # --watson specifies only translate top strand (not reverse strand)
  # -m specifies only AUG allowed as start
  $cmd = $esl_translate . " --watson -m $cur_fafile > $tmp_aa_fafile";
  
  runCommand($cmd, 0);

  # now we have to parse that file to only keep the full length protein seqs
  # we also keep track of which sequence accessions we have full length proteins for
  my $prot_must_start_at_posn_1 = 1; # we only want to fetch sequences that start at position 1
  my $prot_must_stop_at_posn_L  = 1; # we only want to fetch sequences that stop  at position L (final position, that is a valid stop exists there)
  parseEslTranslateOutput($tmp_aa_fafile, $aa_full_fafile, $prot_must_start_at_posn_1, $prot_must_stop_at_posn_L, \%{$fullprot_AH[$c]}, \$nfullprot_A[$c]);
  #printf("CDS: %d nfullprot: %d\n", ($c+1), $nfullprot_A[$c]);
  
  ## now fetch any protein sequences that start at position 1 but do not end at the final predicted position
  #$prot_must_start_at_posn_1 = 1; # we only want to fetch sequences that start at position 1
  #$prot_must_stop_at_posn_L  = 0; # we only want to fetch sequences that stop  at position L (final position, that is a valid stop exists there)
  #parseEslTranslateOutput($tmp_aa_fafile, $aa_trunc_fafile, $prot_must_start_at_posn_1, $prot_must_stop_at_posn_L, \%{$truncprot_AH[$c]}, \$ntruncprot_A[$c]);
  #printf("CDS: %d ntruncprot: %d\n", ($c+1), $ntruncprot_A[$c]);
}

#############################
# CREATE MULTIPLE ALIGNMENTS
#############################
if($do_hmmer) { 
  alignHits($hmmalign, $hmmfetch, $model_db, $do_skipaln, \@model_A, \@seq_accn_A, 
            (($do_nocorrect) ? \%pred_fafile_H : \%corr_fafile_H), 
            \%p_start_HH, \%p_fid2ref_HH, \%p_refdel_HHA, \%p_refins_HHA);
}
else { 
  alignHits($cmalign, $cmfetch, $model_db, $do_skipaln, \@model_A, \@seq_accn_A, 
            (($do_nocorrect) ? \%pred_fafile_H : \%corr_fafile_H), 
            \%p_start_HH, \%p_fid2ref_HH, \%p_refdel_HHA, \%p_refins_HHA);
}

#########################
# OUTPUT ANNOTATION TABLE
#########################

if(1) { # output sequences as rows 
  outputSeqRowHeadings($do_nofid, $do_nomdlb, $do_noss3, $do_nostop, $origin_seq, $ref_tot_nexons, $nhmm, \@hmm2cds_map_A, \@hmm2exon_map_A, \@hmm_is_first_A, \@hmm_is_final_A, \@cds_out_short_A, \@cds_out_product_A);
}
# eventually we may have additional output modes here, not just 1 seq per row
#
#######################################################################
# Pass through all accessions, and output predicted annotation for each
#######################################################################
my @ref_ol_AA = (); # 2D array that describes the overlaps in the reference, $ref_ol_AA[$i][$j] is '1' if the exons modeled by model $i and $j overlap
my $width;          # width of a field

for(my $a = 0; $a < $naccn; $a++) { 
  my $accn = $accn_A[$a];
  my $seq_accn = $seq_accn_A[$a];
  # sanity checks
  if(! exists $totlen_H{$accn}) { die "ERROR accession $accn does not exist in the length file $length_file"; }
  
  # Create the initial portion of the output line, the accession and length
  printf("%-20s  %6d ", $accn, $totlen_H{$accn});

  #########################################################
  # Get information on the actual annotation of this genome
  #########################################################
  my @act_exon_starts_AA = (); # [0..$ncds-1][0..$nexons-1] start positions of actual annotations of exons for this accn, $nexons is CDS specific
  my @act_exon_stops_AA  = (); # [0..$ncds-1][0..$nexons-1] stop  positions of actual annotations of exons for this accn, $nexons is CDS specific
  my $tot_nexons = 0;
  if(exists ($cds_tbl_HHA{$accn})) { 
    ($ncds, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats(\%cds_tbl_HHA, $accn);
    my @cds_len_A = ();
    my @cds_coords_A = ();
    my @cds_product_A = ();
    getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $accn, \@cds_len_A, \@cds_coords_A);
    getQualifierValues(\%cds_tbl_HHA, $accn, "product", \@cds_product_A);
    for(my $i = 0; $i < $ncds; $i++) { 
      # determine start and stop positions of all exons
      my @starts_A = ();
      my @stops_A  = ();
      my $nexons   = 0;
      @{$act_exon_starts_AA[$i]} = ();
      @{$act_exon_stops_AA[$i]}  = ();
      startStopsFromCoords($cds_coords_A[$i], \@starts_A, \@stops_A, \$nexons);

      my $strand = substr($strand_str, $i, 1);
      if($strand eq "-") { # switch order of starts and stops, because 1st exon is really last and vice versa
        @starts_A = reverse @starts_A;          # exons will be in reverse order, b/c we're on the negative strand
        @stops_A  = reverse @stops_A;           # exons will be in reverse order, b/c we're on the negative strand
        @{$act_exon_starts_AA[$i]} = @stops_A;  # save stops  to starts array b/c we're on the negative strand
        @{$act_exon_stops_AA[$i]}  = @starts_A; # save starts to stops  array b/c we're on the negative strand
      }
      else { 
        @{$act_exon_starts_AA[$i]} = @starts_A;
        @{$act_exon_stops_AA[$i]}  = @stops_A;
      }
      $tot_nexons += $nexons;
    }

    # Debugging print statements:
    # printf("\n");
    # printf("$accn\n");
    # for(my $zz = 0; $zz < scalar(@act_exon_starts_AA); $zz++) { 
    #  for(my $zzz = 0; $zzz < scalar(@{$act_exon_starts_AA[$zz]}); $zzz++) { 
    #    printf("act_exon_AA[$zz][$zzz]: $act_exon_starts_AA[$zz][$zzz]  $act_exon_stops_AA[$zz][$zzz]\n");
    #  }
    #  printf("\n");
    #}
    #printf("\n");
  }
  else { 
    $ncds       = 0;
    $tot_nexons = 0;
  }

  ############################################################
  # Create the predicted annotation portion of the output line
  my $predicted_string = "";
  my $nmatch_boundaries = 0;
  my $start_codon_posn;
  my $stop_codon_posn;
  my $start_codon;
  my $stop_codon;
  my $start_codon_char;
  my $stop_codon_char;
  my $multiple_of_3_char;
  my $ss3_yes_char = ".";
  my $ss3_no_char  = "!";
  my $hit_length;
  my $at_least_one_fail; # set to '1' for each CDS if any of the 'tests' for that CDS fail
  my $pass_fail_char; # "P" or "F"
  my $pass_fail_str;  # string of pass_fail_chars

  # data structures we use for checking for overlapping annotation
  my @ol_name_A   = ();  # [0..$nhmm-1]: name of CDS/exons to print if/when outputting information on overlaps
  my @ol_start_A  = ();  # [0..$nhmm-1]: start  position of CDS/exon for use when checking for overlaps
  my @ol_stop_A   = ();  # [0..$nhmm-1]: stop   position of CDS/exon for use when checking for overlaps
  my @ol_strand_A = ();  # [0..$nhmm-1]: strand position of CDS/exon for use when checking for overlaps
 
  ###############################################################
  # create the origin sequence portion of the output line, if nec
  my $oseq_string = "";
  if(defined $origin_seq) { 
    my $norigin = (exists $origin_coords_HA{$accn}) ? scalar(@{$origin_coords_HA{$accn}}) : 0;;
    if($norigin == 1) { 
      my ($ostart, $ostop) = split(":", $origin_coords_HA{$accn}[0]);
      my $predicted_offset = ($ostart < 0) ? ($ostart + $origin_offset) : ($ostart + $origin_offset - 1);
      # $predicted_offset is now number of nts to shift origin in counterclockwise direction
      if($predicted_offset > ($totlen_H{$accn} / 2)) { # simpler (shorter distance) to move origin clockwise
        $predicted_offset = ($totlen_H{$accn} - $predicted_offset + 1);
      }
      else { # simpler to shift origin in counterclockwise direction, we denote this as a negative offset
        $predicted_offset *= -1;
      }
      $pass_fail_char = "P";
      $oseq_string .= sprintf("%2d %5d %5d %5d  %s  ", 1, $ostart, $ostop, $predicted_offset, "P");
    }
    else { 
      $pass_fail_char = "F";
      $oseq_string .= sprintf("%2d %5s %5s %5s  %s  ", $norigin, "-", "-", "-", "F");
    }
    $pass_fail_str .= $pass_fail_char;
  }
  print $oseq_string;
  ###############################################################

  # now the per-exon predictions:
  my $tot_fid = 0.; # all fractional identities added together
  my $n_fid = 0;    # number of fractional identities

  # arrays of data structures that pertain to all exons of current CDS
  my @cur_model_A = (); # array of all model names of exons in current CDS
  my $cur_nexons  = 0;  # number of exons in current CDS

  for(my $h = 0; $h < $nhmm; $h++) { 
    my $model   = $model_A[$h];
    my $cds_i   = $hmm2cds_map_A[$h];
    my $exon_i  = $hmm2exon_map_A[$h];

    if($hmm_is_first_A[$h]) {
      # reset these
      $hit_length = 0; 
      $at_least_one_fail = 0;
      @cur_model_A = ();
      $cur_nexons = 0;
    }
    push(@cur_model_A, $model);
    $cur_nexons++;

    if($predicted_string ne "") { $predicted_string .= "  "; }
    if(exists $p_start_HH{$model}{$seq_accn}) { 
      my $start_corrected_exon = (exists $did_corr_exon_start_AH[$h]{$accn}) ? 1 : 0;
      my $stop_corrected_exon  = (exists $did_corr_exon_stop_AH[$h]{$accn})  ? 1 : 0;
      my $start_corrected_cds  = ($corr_cds_start_AH[$hmm2cds_map_A[$h]]{$accn} != 0) ? 1 : 0;
      my $stop_corrected_cds   = ($corr_cds_stop_AH[$hmm2cds_map_A[$h]]{$accn}  != 0) ? 1 : 0;

      my $start    = ($do_nocorrect) ? $p_start_HH{$model}{$seq_accn} : $c_start_HH{$model}{$seq_accn};
      my $stop     = ($do_nocorrect) ? $p_stop_HH{$model}{$seq_accn}  : $c_stop_HH{$model}{$seq_accn};
      my $hangover = $p_hangover_HH{$model}{$seq_accn};

      my ($hang5, $hang3) = split(":", $hangover);
      if($hang5    >  9) { $hang5 = "+"; $at_least_one_fail = 1; }
      elsif($hang5 == 0) { $hang5 = "."; }
      if($start_corrected_exon) { $hang5 = "c"; }

      if($hang3       >  9) { $hang3 = "+"; $at_least_one_fail = 1; }
      elsif($hang3    == 0) { $hang3 = "."; }
      if($stop_corrected_exon) { $hang3 = "c"; }

      my ($start_match, $stop_match);
      ($start_match, $stop_match) = ($do_strict) ? 
          checkStrictBoundaryMatch   (\@act_exon_starts_AA, \@act_exon_stops_AA, $cds_i, $exon_i, $start, $stop) :
          checkNonStrictBoundaryMatch(\@act_exon_starts_AA, \@act_exon_stops_AA, $start, $stop);
      if($start_match) { $nmatch_boundaries++; }
      if($stop_match)  { $nmatch_boundaries++; }
 
      if($do_nobrack) { # set to '1' so brackets are never printed
        $start_match = 1;
        $stop_match  = 1; 
      }

      $hit_length += abs($stop-$start) + 1;
      if(($stop < 0 && $start > 0) || 
         ($stop > 0 && $start < 0)) { 
        # correct for off-by-one induced by the way we use negative indices distance from -1..1 is 1 nt, not 2
        $hit_length -= 1;
      }

      # TODO: MODIFY ANNOTATION FOR EXONS WITHOUT CORRECTED STARTS OR STOPS IN CDS WITH
      #       OTHER EXONS THAT HAVE CORRECTED STARTS OR STOPS
      $predicted_string .= sprintf("%8s %8s",
                                   ($start_match ? " " . $start . " " : "[" . $start . "]"), 
                                   ($stop_match  ? " " . $stop .  " " : "[" . $stop . "]"));
      if(! $do_nofid) { 
        $predicted_string .= sprintf(" %5.3f",
                                   $p_fid2ref_HH{$model}{$seq_accn});
      }
      $tot_fid += $p_fid2ref_HH{$model}{$seq_accn};
      $n_fid++;

      if(! $do_nomdlb) { 
        $predicted_string .= "  " . $hang5 . $hang3;
      }        
                                   
      if($hmm_is_first_A[$h]) { # determine $start_codon_char
        if($p_strand_HH{$model}{$seq_accn} eq "-") { 
          $start_codon_posn = (($start-2) < 0) ? $start + $totlen_H{$accn} + 1 : $start;
        }
        else { 
          $start_codon_posn = ($start < 0) ? $start + $totlen_H{$accn} + 1 : $start;
        }
        $start_codon = fetchCodon($sqfile, $seq_accn, $start_codon_posn, $p_strand_HH{$model}{$seq_accn});
        if($start_codon eq "ATG") { 
          $start_codon_char = $ss3_yes_char;
        }
        else { 
          $start_codon_char = $ss3_no_char;
          $at_least_one_fail = 1;
        }
      }
      
      if($hmm_is_final_A[$h]) { 
        if($p_strand_HH{$model}{$seq_accn} eq "-") { 
          $stop_codon_posn    = ($stop < 0) ? ($stop + $totlen_H{$accn}) + 1 + 2 : $stop + 2;
        }
        else { 
          $stop_codon_posn    = (($stop-2) < 0) ? ($stop + $totlen_H{$accn}) + 1 - 2 : $stop - 2;
        }
        $stop_codon         = fetchCodon($sqfile, $seq_accn, $stop_codon_posn, $p_strand_HH{$model}{$seq_accn});

        if($stop_codon eq "TAG" || $stop_codon eq "TAA" || $stop_codon eq "TGA") { 
          $stop_codon_char = $ss3_yes_char;
        }
        else { 
          $stop_codon_char = $ss3_no_char;
          $at_least_one_fail = 1;
        }
        if(($hit_length % 3) == 0) { 
          $multiple_of_3_char = $ss3_yes_char;
        }
        else { 
          $multiple_of_3_char = $ss3_no_char;
          $at_least_one_fail = 1;
        }
        # append the ss3 (start/stop/multiple of 3 info)
        $predicted_string .= sprintf(" %6d", $hit_length);
        if(! $do_noss3) { 
          $predicted_string .= sprintf(" %s%s%s", $start_codon_char, $stop_codon_char, $multiple_of_3_char);
        }
        if(! $do_nostop) { 
          $predicted_string .= sprintf(" %3s", $stop_codon);
        }

        if($at_least_one_fail) { 
          $pass_fail_char = "F"; 
          #### ensure we didn't fetch a full length protein for this CDS
          ###if(exists $fullprot_AH[$cds_i]{$accn}) { 
          ###  die sprintf("ERROR, incorrectly translated full length protein for CDS: %d for seq accn: $accn, but at least one test failed...", $cds_i+1); 
          ###}
          ###if(! exists $truncprot_AH[$cds_i]{$accn}) { 
          ###  die sprintf("ERROR, failed to translate truncated protein for CDS: %d for seq accn: $accn, but at least one test failed...", $cds_i+1); 
          ###}
        }
        else { 
          $pass_fail_char = "P"; 
          #### verify that we have a translated protein sequence for this CDS
          ###if(! exists $fullprot_AH[$cds_i]{$accn}) { 
          ###  die sprintf("ERROR, failed to translate full length protein for CDS: %d for seq accn: $accn, but all tests passed...", $cds_i+1); 
          ###}
          ###if(exists $truncprot_AH[$cds_i]{$accn}) { 
          ###  die sprintf("ERROR, incorrectly translated truncated protein for CDS: %d for seq accn: $%accn, but all tests passed...", $cds_i+1);
          ###}
        }
        $pass_fail_char = ($at_least_one_fail) ? "F" : "P";
        $predicted_string .= sprintf(" %2s", $pass_fail_char);
        $pass_fail_str .= $pass_fail_char;
      }

      # save information for overlap check
      push(@ol_name_A, sprintf "%d.%d", $hmm2cds_map_A[$h]+1, $hmm2exon_map_A[$h]+1);
      push(@ol_start_A, $start);
      push(@ol_stop_A,  $stop);
      push(@ol_strand_A, $p_strand_HH{$model}{$seq_accn});
    }
    else { 
      # printf("no hits for $model $seq_accn\n");
      if($do_nomdlb) { 
        $width = ($hmm_is_final_A[$h]) ? 34 : 23;
        $predicted_string .= sprintf("%*s", $width, "NO PREDICTION");
      }
      else { 
        $width = ($hmm_is_final_A[$h]) ? 38 : 27;
        $predicted_string .= sprintf("%*s", $width, "NO PREDICTION");
      }
    }
  }
  print $predicted_string;

  printf("  %6d", $totlen_H{$accn});
  printf("  %5.3f", $tot_fid / $n_fid);

  # output number of actually annotated CDS and summed total of exons in those CDS, if nec
  if(! $do_noexist) { 
    printf("  %5d  %5d  %5d", $ncds, $tot_nexons, $nmatch_boundaries);
  }

  # check for overlaps
  my $overlap_notes;
  if($a == 0) { 
    # the reference, determine which overlaps are allowed, we'll fill @allowed_ol_AA in checkForOverlaps
    ($pass_fail_char, $overlap_notes) = checkForOverlaps(\@ol_name_A, \@ol_start_A, \@ol_stop_A, \@ol_strand_A, undef, \@ref_ol_AA);
  }
  else { 
    # not the reference, we'll determine if this accession 'passes' the overlap test based on whether it's observed
    # overlaps match those in \@allowed_ol_AA exactly or not
    ($pass_fail_char, $overlap_notes) = checkForOverlaps(\@ol_name_A, \@ol_start_A, \@ol_stop_A, \@ol_strand_A, \@ref_ol_AA, undef);
  }
  $pass_fail_str .= $pass_fail_char;

  # output overlap info, if nec
  if((! $do_noolap) && ($overlap_notes ne "")) { 
    printf("  %20s", $overlap_notes);
  }

  my $result_str = ($pass_fail_str =~ m/F/) ? "FAIL" : "PASS";
  $result_str .= " " . $pass_fail_str;
  printf("  %s", $result_str);

  print "\n";
}

##########################
# OUTPUT EXPLANATORY TEXT 
##########################
if(! $do_noexp) { 
  printColumnHeaderExplanations((defined $origin_seq), $do_nomdlb, $do_noexist, $do_nobrack, $do_nostop, $do_nofid, $do_noss3, $do_noolap);
}

###########################
# OUTPUT THE GAP INFO FILES
###########################
# first, the gap file that lists all gaps, and all gaps that are not multiples of 3
my $gap_perseq_all_file     = $out_root . ".gap.perseq-all.txt";
my $gap_perseq_not3_file    = $out_root . ".gap.perseq-not3.txt";
my $gap_perseq_special_file = $out_root . ".gap.perseq-special.txt";

my $gap_pergap_all_file     = $out_root . ".gap.pergap-all.txt";
my $gap_pergap_not3_file    = $out_root . ".gap.pergap-not3.txt";
my $gap_pergap_special_file = $out_root . ".gap.pergap-special.txt";

my $gap_perseq_all_FH     = undef;
my $gap_perseq_not3_FH    = undef;
my $gap_perseq_special_FH = undef;

my $gap_pergap_all_FH     = undef;
my $gap_pergap_not3_FH    = undef;
my $gap_pergap_special_FH = undef;

open($gap_perseq_all_FH,     ">", $gap_perseq_all_file)     || die "ERROR unable to open $gap_perseq_all_file"; 
open($gap_perseq_not3_FH,    ">", $gap_perseq_not3_file)    || die "ERROR unable to open $gap_perseq_not3_file"; 
open($gap_perseq_special_FH, ">", $gap_perseq_special_file) || die "ERROR unable to open $gap_perseq_special_file"; 
open($gap_pergap_all_FH,     ">", $gap_pergap_all_file)     || die "ERROR unable to open $gap_pergap_all_file"; 
open($gap_pergap_not3_FH,    ">", $gap_pergap_not3_file)    || die "ERROR unable to open $gap_pergap_not3_file"; 
open($gap_pergap_special_FH, ">", $gap_pergap_special_file) || die "ERROR unable to open $gap_pergap_special_file"; 

# 3 calls to outputGapInfo to create all the files we need:
outputGapInfo($gap_perseq_all_FH,     $gap_pergap_all_FH,     0, 1, 0, 0, \@model_A, \@seq_accn_A, \%mdllen_H, \@hmm2cds_map_A, \%p_refdel_HHA, \%p_refins_HHA);
outputGapInfo($gap_perseq_not3_FH,    $gap_pergap_not3_FH,    0, 0, 1, 0, \@model_A, \@seq_accn_A, \%mdllen_H, \@hmm2cds_map_A, \%p_refdel_HHA, \%p_refins_HHA);
outputGapInfo($gap_perseq_special_FH, $gap_pergap_special_FH, 0, 0, 0, 1, \@model_A, \@seq_accn_A, \%mdllen_H, \@hmm2cds_map_A, \%p_refdel_HHA, \%p_refins_HHA);

close($gap_perseq_all_FH);
close($gap_perseq_not3_FH);
close($gap_perseq_special_FH);
close($gap_pergap_all_FH);
close($gap_pergap_not3_FH);
close($gap_pergap_special_FH);

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


# Subroutine: validateRefCDSAreUnique()
# Purpose:    Validate that all CDS:product annotation for all reference CDS 
#             are unique, i.e. there are no two CDS that have the same value
#             in their CDS:product annotation.
#
# Args:       $ref_ncds:           number of reference CDS
#             $ref_cds_product_AR: ref to array of CDS:product annotations for the $ref_ncds reference CDS 
#
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
# Returns:    void; but fills @{$starts_AR}, @{$stops_AR}, and $$nexons_R.
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

# Subroutine: createHmmDb()
# Synopsis:   Create an HMM Db from a stockholm database file.
#
# Args:       $hmmbuild:   path to 'hmmbuild' executable
#             $hmmpress:   path to 'hmmpress' executable
#             $nmodel:     number of models we're creating
#             $stk_file:   stockholm DB file
#             $out_root:   string for naming output files
#
# Returns:    void
#
sub createHmmDb { 
  my $sub_name = "createHmmDb()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($hmmbuild, $hmmpress, $nmodel, $stk_file, $out_root) = @_;

  if(! -s $stk_file)  { die "ERROR in $sub_name, $stk_file file does not exist or is empty"; }

  # remove the binary files, possibly from an earlier hmmbuild/hmmpress:
  for my $suffix ("h3f", "h3i", "h3m", "h3p") { 
    my $file = $out_root . ".hmm." . $suffix;
    if(-e $file) { unlink $file; }
  }

  # first build the models

  printf("%-65s ... ", sprintf("# Running hmmbuild to build %d HMMs", $nmodel));
  my $cmd = "$hmmbuild --dna $out_root.hmm $stk_file > $out_root.hmmbuild";
  my $secs_elapsed = runCommand($cmd, 0);
#  printf("done. [$out_root.nhmmer and $out_root.tblout]\n");
  printf("done. [%.1f seconds]\n", $secs_elapsed);

  # next, press the HMM DB we just created
  printf("%-65s ... ", "# Running hmmpress");
  $cmd = "$hmmpress $out_root.hmm > $out_root.hmmpress";
  $secs_elapsed = runCommand($cmd, 0);
#  printf("done. [$out_root.nhmmer and $out_root.tblout]\n");
  printf("done. [%.1f seconds]\n", $secs_elapsed);

  return;
}

# Subroutine: createCmDb()
# Synopsis:   Create an CM database from a stockholm database file
#             for use with Infernal 1.1.
#             the $cmbuild executable. If $cmcalibrate is defined
#             also run cmcalibrate. 
#
# Args:       $cmbuild:          path to 'cmbuild' executable
#             $cmcalibrate:      path to 'cmcalibrate' executable
#             $cmpress:          path to 'cmpress' executable
#             $nmodel:           number of models we're creating/calibrating
#             $do_calib_slow:    '1' to calibrate using default parameters instead of
#                                options to make it go much faster
#             $do_calib_cluster: '1' to submit calibration job to cluster, '0' to do it locally
#             $stk_file:         stockholm DB file
#             $out_root:         string for naming output files
#
# Returns:    void
#
sub createCmDb { 
  my $sub_name = "createCmDb()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cmbuild, $cmcalibrate, $cmpress, $nmodel, $do_calib_slow, $do_calib_cluster, $stk_file, $out_root) = @_;

  if(! -s $stk_file)  { die "ERROR in $sub_name, $stk_file file does not exist or is empty"; }

  # remove the binary files, possibly from an earlier cmbuild/cmpress:
  for my $suffix ("i1m", "i1i", "i1f", "i1p") { 
    my $file = $out_root . ".cm." . $suffix;
    if(-e $file) { unlink $file; }
  }

  my ($cmbuild_opts,     $cmbuild_cmd);
  my ($cmcalibrate_opts, $cmcalibrate_cmd);
  my ($cmpress_opts,     $cmpress_cmd);

  $cmbuild_opts = "-F";
  $cmbuild_cmd  = "$cmbuild $cmbuild_opts $out_root.cm $stk_file > $out_root.cmbuild";

  $cmcalibrate_opts = " --cpu 4 ";
  if(! $do_calib_slow) { $cmcalibrate_opts .= " -L 0.04 "; }
  $cmcalibrate_cmd  = "$cmcalibrate $cmcalibrate_opts $out_root.cm > $out_root.cmcalibrate";
  
  $cmpress_cmd = "$cmpress $out_root.cm > $out_root.cmpress";

  # first build the models
  printf("%-65s ... ", sprintf("# Running cmbuild to build %d CMs", $nmodel));
  my $secs_elapsed = runCommand($cmbuild_cmd, 0);
  printf("done. [%.1f seconds]\n", $secs_elapsed);

  if($do_calib_cluster) { 
    # submit a job to the cluster and exit. 
    my $out_tail = $out_root;
    $out_tail =~ s/^.+\///;
    my $jobname = "cp." . $out_tail;
    my $errfile = $out_root . ".err";
    my $cluster_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n -l h_rt=288000,h_vmem=8G,mem_free=8G -pe multicore 4 -R y " . "\"" . $cmcalibrate_cmd . ";" . $cmpress_cmd . ";\"\n";
    # print("$cluster_cmd\n");
    runCommand($cluster_cmd, 0);
  }
  else { 
    # calibrate the model
    printf("%-65s ... ", "# Running cmcalibrate");
    $secs_elapsed = runCommand($cmcalibrate_cmd, 0);
    #printf("\n$cmcalibrate_cmd\n");
    printf("done. [%.1f seconds]\n", $secs_elapsed);

    # press the model
    printf("%-65s ... ", "# Running cmpress");
    $secs_elapsed = runCommand($cmpress_cmd, 0);
    #printf("\n$cmpress_cmd\n");
    printf("done [%.1f seconds]\n", $secs_elapsed);
  } # end of 'else' entered if $do_calib_cluster is false

  return;
}

# Subroutine: runNhmmscan()
# Synopsis:   Perform a homology search using nhmmscan.
#
# Args:       $nhmmscan:     path to nhmmscan executable
#             $do_skip:      '1' to not actually run hmmscan, but verify expected output already exists
#             $model_db:     path to model HMM database
#             $seq_fasta:    path to seq fasta file
#             $tblout_file:  path to --tblout output file to create, undef to not create one
#             $stdout_file:  path to output file to create with standard output from nhmmscan, undef 
#                            to pipe to /dev/null
#
# Returns:    void
#
sub runNhmmscan { 
  my $sub_name = "runNhmmscan()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($nhmmscan, $do_skip, $model_db, $seq_fasta, $tblout_file, $stdout_file) = @_;

  if($do_skip) { 
    printf("%-65s ... ", "# Skipping nhmmscan step");
    if(! -s $tblout_file) { die "ERROR, with -skipscan or -skipaln nhmmscan output is expected to already exist, but $tblout_file does not or is empty"; }
    if(! -s $stdout_file) { die "ERROR, with -skipscan or -skipaln nhmmscan output is expected to already exist, but $stdout_file does not or is empty"; }
    printf("done. [-skipscan or -skipaln]\n");
    return;
  }

  # my $opts = " --noali --tblout $tblout_file ";
  my $opts = " --tblout $tblout_file ";

  if(! defined $stdout_file) { $stdout_file = "/dev/null"; }

  if(! -s $model_db)  { die "ERROR in $sub_name, $model_db file does not exist or is empty"; }
  if(! -s $seq_fasta) { die "ERROR in $sub_name, $seq_fasta file does not exist or is empty"; }

  printf("%-65s ... ", "# Running nhmmscan");
  my $cmd = "$nhmmscan $opts $model_db $seq_fasta > $stdout_file";
  my $secs_elapsed = runCommand($cmd, 0);
  printf("done. [%.1f seconds]\n", $secs_elapsed);

  return;
}

# Subroutine: runCmscan()
# Synopsis:   Run Infernal 1.1's cmscan.
#
# Args:       $cmscan:      path to cmscan executable
#             $do_glocal:   '1' to use the -g option, '0' not to
#             $do_skip:     '1' to not actually run cmscan, but verify expected output already exists
#             $do_cluster:  '1' to submit job to cluster, instead of running it locally
#             $model_db:    path to model CM database
#             $seq_fasta:   path to seq fasta file
#             $tblout_file: path to --tblout output file to create, undef to not create one
#             $stdout_file: path to output file to create with standard output from cmsearch, undef 
#                           to pipe to /dev/null
#
# Returns:    void
#
sub runCmscan { 
  my $sub_name = "runCmscan()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cmscan, $do_glocal, $do_skip, $do_cluster, $model_db, $seq_fasta, $tblout_file, $stdout_file) = @_;

  if($do_skip) { 
    printf("%-65s ... ", "# Skipping cmscan step");
    if(! -s $tblout_file) { die "ERROR, with -skipscan or -skipaln cmscan output is expected to already exist, but $tblout_file does not or is empty"; }
    if(! -s $stdout_file) { die "ERROR, with -skipscan or -skipaln cmscan output is expected to already exist, but $stdout_file does not or is empty"; }
    printf("done. [-skipscan or -skipaln]\n");
    return;
  }

  my $opts = "";
  if($do_iglocal) { $opts .= "-g "; }
  $opts .= " --cpu 0 --rfam --tblout $tblout_file --verbose --nohmmonly ";
  if(! defined $stdout_file) { $stdout_file = "/dev/null"; }

  if(! -s $model_db)   { die "ERROR in $sub_name, $model_db file does not exist or is empty"; }
  if(! -s $seq_fasta) { die "ERROR in $sub_name, $seq_fasta file does not exist or is empty"; }

  my $cmd = "$cmscan $opts $model_db $seq_fasta > $stdout_file";
  if(! $do_cluster) { 
    # default mode, run job locally
    printf("%-65s ... ", "# Running cmscan");
    my $secs_elapsed = runCommand($cmd, 0);
    printf("done. [%.1f seconds]\n", $secs_elapsed);
  }
  else { 
    # submit job to cluster and return
    my $jobname = $seq_fasta;
    my $errfile = $stdout . ".err";
    $jobname =~ s/^.+\///; # remove everything up until final '/'
    my $cluster_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n -l h_rt=288000,h_vmem=8G,mem_free=8G " . "\"" . $cmd . "\" > /dev/null\n";
    runCommand($cluster_cmd, 0);
  }
  return;
}

# Subroutine: annotateStockholmAlignment
#
# Synopsis:   Read in a stockholm alignment ($in_file), and 
#             add a name ($name) to it, then optionally add
#             a blank SS (#=GC SS_cons) annotation to it
#             and output a new file ($out_file) that is
#             identical to it but with the name annotation
#             and possibly a blank SS.
#
# Args:       $name:          name to add to alignment
#             $do_blank_ss:   '1' to add a blank SS, else '0'
#             $in_file:       input stockholm alignment
#             $out_file:      output stockholm alignment to create
#
# Returns:    alignment length
#
sub annotateStockholmAlignment {
  my $sub_name = "annotateStockholmAlignment";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($name, $do_blank_ss, $in_file, $out_file) = @_;

  # printf("# Naming alignment in $in_file to $name ... "); 

  # open and validate file
  my $msa = Bio::Easel::MSA->new({
    fileLocation => $in_file,
                                 });  
  $msa->set_name($name);
  if($do_blank_ss) { 
    $msa->set_blank_ss_cons;
  }
  $msa->write_msa($out_file);

  # printf("done. [$out_file]\n");
  return $msa->alen;
}

# Subroutine: parseNhmmscanTblout
#
# Synopsis:   Parse nhmmscan tblout output into 5 2D hashes.
#             For each 2D hash first key is seq name, second key
#             is model name, value is either start, stop, strand,
#             score or hangover (number of model positions not included
#             on 5' and 3' end). Information for the lowest E-value hit
#             for each seq/model pair is stored. This will be the
#             first hit encountered in the file for each seq/model
#             pair.
#
# Args:       $tblout_file:   tblout file to parse
#             $do_hmmenv:     '1' to use envelope boundaries, else use window boundaries
#             $totlen_HR:     ref to hash, key is accession, value is length, pre-filled
#             $start_HHR:     ref to 2D hash of start values, to fill here
#             $stop_HHR:      ref to 2D hash of stop values, to fill here
#             $strand_HHR:    ref to 2D hash of strand value, to fill here
#             $score_HHR:     ref to 2D hash of score values, to fill here
#             $hangoverHHR:   ref to 2D hash of model hangover values, to fill here
#
# Returns:    void
#
sub parseNhmmscanTblout { 
  my $sub_name = "parseNhmmscanTblout";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($tblout_file, $do_hmmenv, $totlen_HR, $start_HHR, $stop_HHR, $strand_HHR, $score_HHR, $hangover_HHR) = @_;
  
  open(IN, $tblout_file) || die "ERROR unable to open $tblout_file for reading";

  my $line_ctr = 0;
  while(my $line = <IN>) { 
    $line_ctr++;
    if($line =~ m/^\# \S/ && $line_ctr == 1) { 
      # sanity check, make sure the fields are what we expect
      chomp $line;
      if($line !~ m/#\s+target\s+name\s+accession\s+query name\s+accession\s+hmmfrom\s+hmm to\s+alifrom\s+ali to\s+envfrom  env to\s+modlen\s+strand\s+E-value\s+score\s+bias\s+description of target/) { 
        die "ERROR unexpected field names in $tblout\n$line\n";
      }
    }
    elsif($line !~ m/^\#/) { 
      chomp $line;
      #Maize-streak_r23.NC_001346/Maize-streak_r23.NC_001346.ref.cds.4        -          NC_001346:genome:NC_001346:1:2689:+: -                1     819    2 527    1709    2527    1709     819    -    9.8e-261  856.5  12.1  -
      my @elA = split(/\s+/, $line);
      my ($mdl, $seq, $hmmfrom, $hmmto, $alifrom, $alito, $envfrom, $envto, $mdllen, $strand, $score) = 
          ($elA[0], $elA[2], $elA[4], $elA[5], $elA[6], $elA[7], $elA[8], $elA[9], $elA[10], $elA[11], $elA[13]);

      my $from = ($do_hmmenv) ? $envfrom : $alifrom;
      my $to   = ($do_hmmenv) ? $envto   : $alito;

      my $accn = $seq;
      $accn =~ s/\:.+$//;
      if(! exists $totlen_HR->{$accn}) { die "ERROR unable to determine accession with stored length from fasta sequence $mdl (seq: $seq)"; }
      my $L = $totlen_HR->{$accn};

      storeHit($mdl, $seq, $mdllen, $L, $hmmfrom, $hmmto, $from, $to, $strand, $score, $start_HHR, $stop_HHR, $strand_HHR, $score_HHR, $hangover_HHR);
    }
  }
  close(IN);
  
  return;
}

# Subroutine: parseCmscanTblout
#
# Synopsis:   Parse Infernal 1.1 cmscan --tblout output.
#             For each 2D hash first key is seq name, second key
#             is model name, value is either start, stop, strand,
#             score or hangover (number of model positions not included
#             on 5' and 3' end). Information for the lowest E-value hit
#             for each seq/model pair is stored. This will be the
#             first hit encountered in the file for each seq/model
#             pair.
#
# Args:       $tblout_file:   tblout file to parse
#             $totlen_HR:     ref to hash, key is accession, value is length, pre-filled
#             $mdllen_HR:     ref to hash, key is model name, value is model length, pre-filled
#             $start_HHR:     ref to 2D hash of start values, to fill here
#             $stop_HHR:      ref to 2D hash of stop values, to fill here
#             $strand_HHR:    ref to 2D hash of strand value, to fill here
#             $score_HHR:     ref to 2D hash of score values, to fill here
#             $hangoverHHR:   ref to 2D hash of model hangover values, to fill here
#
# Returns:    void
#
sub parseCmscanTblout { 
  my $sub_name = "parseCmscanTblout";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($tblout_file, $totlen_HR, $mdllen_HR, $start_HHR, $stop_HHR, $strand_HHR, $score_HHR, $hangover_HHR) = @_;
  
  open(IN, $tblout_file) || die "ERROR unable to open $tblout_file for reading";

  my $line_ctr = 0;
  while(my $line = <IN>) { 
    $line_ctr++;
    if($line =~ m/^\# \S/ && $line_ctr == 1) { 
      # sanity check, make sure the fields are what we expect
      chomp $line;
      if($line !~ m/#target name\s+accession\s+query name\s+accession\s+mdl\s+mdl\s+from\s+mdl to\s+seq from\s+seq to\s+strand\s+trunc\s+pass\s+gc\s+bias\s+score\s+E-value inc description of target/) { 
        die "ERROR unexpected field names in $tblout\n$line\n";
      }

    }
    elsif($line !~ m/^\#/) { 
      chomp $line;
      #Maize-streak_r23.NC_001346.ref.cds.4        -         NC_001346:genome-duplicated:NC_001346:1:2689:+:NC_001346:1:2689:+: -          cm        1      819     2527     1709      -    no    1 0.44   0.2  892.0         0 !   -
      my @elA = split(/\s+/, $line);
      my ($mdl, $seq, $mod, $mdlfrom, $mdlto, $from, $to, $strand, $score) = 
          ($elA[0], $elA[2], $elA[4], $elA[5], $elA[6], $elA[7], $elA[8], $elA[9], $elA[14]);

      my $accn = $seq;
      $accn =~ s/\:.+$//;
      if(! exists $totlen_HR->{$accn})  { die "ERROR unable to determine accession with stored length from fasta sequence $mdl"; }
      if(! exists $mdllen_HR->{$mdl})   { die "ERROR do not have model length information for model $mdl"; }
      my $L      = $totlen_HR->{$accn};
      my $mdllen = $mdllen_HR->{$mdl};

      storeHit($mdl, $seq, $mdllen, $L, $mdlfrom, $mdlto, $from, $to, $strand, $score, $start_HHR, $stop_HHR, $strand_HHR, $score_HHR, $hangover_HHR);
    }
  }
  close(IN);
  
  return;
}

# Subroutine: storeHit
#
# Synopsis:   Helper function for parseNhmmscanTblout and parseCmscanTblout.
#             Given info on a hit and refs to hashes to store info on it in,
#             store it.
#
# Args:       $mdl:           model name
#             $seq:           sequence name
#             $mdllen:        model length
#             $L:             target sequence length
#             $mdlfrom:       start position of hit
#             $mdlto:         stop position of hit
#             $seqfrom:       start position of hit
#             $seqto:         stop position of hit
#             $strand:        strand of hit
#             $score:         bit score of hit
#             $start_HHR:     ref to 2D hash of start values, to fill here
#             $stop_HHR:      ref to 2D hash of stop values, to fill here
#             $strand_HHR:    ref to 2D hash of strand value, to fill here
#             $score_HHR:     ref to 2D hash of score values, to fill here
#             $hangover_HHR:  ref to 2D hash of model hangover values, to fill here
#                             start..stop coordinates for $qseq.
# Returns:    void
#
sub storeHit { 
  my $sub_name = "storeHit";
  my $nargs_exp = 15;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl, $seq, $mdllen, $L, $mdlfrom, $mdlto, $seqfrom, $seqto, $strand, $score, 
      $start_HHR, $stop_HHR, $strand_HHR, $score_HHR, $hangover_HHR) = @_;

  # only consider hits where either the start or end are less than the total length
  # of the genome. Since we typically duplicate all genomes, this avoids storing 
  # duplicate hits at different positions.
  if(($seqfrom <= $L) || ($seqto <= $L)) { 
    
    # deal with case where one but not both of from to is > L:
    if($seqfrom > $L || $seqto > $L) { 
      $seqfrom -= $L; 
      $seqto   -= $L; 
      if($seqfrom < 0)  { $seqfrom--; }
      if($seqto   < 0)  { $seqto--; }
    }
    
    if(! exists $start_HHR->{$mdl}) { # initialize
      %{$start_HHR->{$mdl}}    = ();
      %{$stop_HHR->{$mdl}}     = ();
      %{$strand_HHR->{$mdl}}   = ();
      %{$score_HHR->{$mdl}}    = ();
      %{$hangover_HHR->{$mdl}} = ();
    }
    if(! exists $start_HHR->{$mdl}{$seq})    { $start_HHR->{$mdl}{$seq}    = $seqfrom; }
    if(! exists $stop_HHR->{$mdl}{$seq})     { $stop_HHR->{$mdl}{$seq}     = $seqto; }
    if(! exists $strand_HHR->{$mdl}{$seq})   { $strand_HHR->{$mdl}{$seq}   = $strand; }
    if(! exists $score_HHR->{$mdl}{$seq})    { $score_HHR->{$mdl}{$seq}    = $score; }
    if(! exists $hangover_HHR->{$mdl}{$seq}) { $hangover_HHR->{$mdl}{$seq} = ($mdlfrom - 1) . ":" . ($mdllen - $mdlto); }
  }

  return;
}

# Subroutine: findSeqInFile
#
# Synopsis:   Identify all exact occurences of a sequence in a file
#             of sequences, and store the coordinates of the
#             matches in %{$coords_HAR}.
#
# Args:       $sqfile:        Bio::Easel::SqFile object, the sequence file to search in
#             $qseq:          query sequence we're looking for
#             $do_nodup:      '1' if -nodup was used at command line, else '0'
#             $coords_HAR:    ref to hash of arrays to store coords in
#                             key is accession, value is array of 
#                             start..stop coordinates for $qseq.
# Returns:    void
#
sub findSeqInFile { 
  my $sub_name = "findSeqInFile";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $qseq, $do_nodup, $coords_HAR) = @_;

  # printf("# Naming alignment in $in_file to $name ... "); 

  # fetch each sequence and look for $qseq in it
  # (could (and probably should) make this more efficient...)
  my $nseq = $sqfile->nseq_ssi();
  for(my $i = 0; $i < $nseq; $i++) { 
    my $seqname   = $sqfile->fetch_seq_name_given_ssi_number($i);
    my $accn      = $seqname;
    # example: AJ224507:genome-duplicated:AJ224507:1:2685:+:AJ224507:1:2685:+:
    $accn =~ s/\:.+$//;
    my $fasta_seq = $sqfile->fetch_seq_to_fasta_string($seqname, -1);
    my ($header, $seq) = split(/\n/, $fasta_seq);
    chomp $seq;
    my $L = ($do_nodup) ? length($seq) : (length($seq) / 2);
    # now use Perl's index() function to find all occurrences of $qseq
    my $qseq_posn = index($seq, $qseq);
    while($qseq_posn != -1) { 
      $qseq_posn++;
      if($qseq_posn <= $L) { 
        my $qseq_start = $qseq_posn;
        my $qseq_stop  = $qseq_posn + length($qseq) - 1;;
        if($qseq_stop > $L) { 
          $qseq_start -= $L;
          $qseq_start -= 1; # off-by-one issue with negative indexing
          $qseq_stop  -= $L;
        }
        if(! exists $coords_HAR->{$accn}) { 
          @{$coords_HAR->{$accn}} = ();
        }
        push(@{$coords_HAR->{$accn}}, $qseq_start . ":" . $qseq_stop);
        # printf("Found $qseq in $accn at position %d..%d\n", $qseq_start, $qseq_stop);
      }
      $qseq_posn = index($seq, $qseq, $qseq_posn);
    }
  }
  
  return;
}

# Subroutine: checkStrictBoundaryMatch
#
# Synopsis:   Check if a given start..stop boundary set matches the 
#             actual annotation in $act_AAR->[$cds_i][$exon_i]
#             (if that array element even exists).
#
# Args:       $act_start_AAR: ref to 2D array [0..i..$ncds-1][0..e..$nexon-1], start for cds $i+1 exon $e+1
#             $act_stop_AAR:  ref to 2D array [0..i..$ncds-1][0..e..$nexon-1], stop for cds $i+1 exon $e+1
#             $cds_i:         CDS index we want to check against
#             $exon_i:        exon index we want to check against
#             $pstart:        predicted start boundary
#             $pstop:         predicted stop boundary
# Returns:    Two values:
#             '1' if $pstart == $act_start_AAR->[$cds_i][$exon_i], else '0'
#             '1' if $pstop  == $act_stop_AAR->[$cds_i][$exon_i], else '0'
#
sub checkStrictBoundaryMatch {
  my $sub_name = "checkStrictBoundaryMatch";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($act_start_AAR, $act_stop_AAR, $cds_i, $exon_i, $pstart, $pstop) = @_;

  my $retval1 = 
      ((exists ($act_start_AAR->[$cds_i])) && 
       (exists ($act_start_AAR->[$cds_i][$exon_i])) && 
       ($pstart == $act_start_AAR->[$cds_i][$exon_i])) ? 
       1 : 0;

  my $retval2 = 
      ((exists ($act_stop_AAR->[$cds_i])) && 
       (exists ($act_stop_AAR->[$cds_i][$exon_i])) && 
       ($pstop == $act_stop_AAR->[$cds_i][$exon_i])) ? 
       1 : 0;

  return ($retval1, $retval2);
}

# Subroutine: checkNonStrictBoundaryMatch
#
# Synopsis:   Check if a given boundary matches any
#             annotation in the 2D array referred to
#             by $act_AAR.
#
# Args:       $act_start_AAR: ref to 2D array [0..i..$ncds-1][0..e..$nexon-1], start for cds $i+1 exon $e+1
#     :       $act_stop_AAR:  ref to 2D array [0..i..$ncds-1][0..e..$nexon-1], stop for cds $i+1 exon $e+1
#             $pstart:        predicted start position
#             $pstop:         predicted stop position
# Returns:    Two values:
#             '1' if $pstart == $act_start_AAR->[$i][$e], else '0'
#             '1' if $pstop  == $act_stop_AAR->[$i][$e], else '0'
#             For any possible $i and $e values, as long as they are the same for start and stop
#             If both ('1', '0') and ('0', '1') are possible sets of return values, 
#             ('1', '0') is returned.
#
sub checkNonStrictBoundaryMatch {
  my $sub_name = "checkNonStrictBoundaryMatch";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($act_start_AAR, $act_stop_AAR, $pstart, $pstop) = @_;

  my $found_start_match = 0;
  my $found_stop_match = 0;
  my $found_both_match = 0;
  my $start_match = 0;
  my $stop_match = 0;
  my $ncds = scalar(@{$act_start_AAR});
  for(my $i = 0; $i < $ncds; $i++) { 
    my $nexons = scalar(@{$act_start_AAR->[$i]});
    if(! exists $act_stop_AAR->[$i]) { die "ERROR in checkNonStrictBoundaryMatch() $i exists in first dimension of start coords, but not stop coords"; }
    for(my $e = 0; $e < $nexons; $e++) { 
      if(! exists $act_stop_AAR->[$i][$e]) { die "ERROR in checkNonStrictBoundaryMatch() $i $e exists in start coords, but not stop coords"; }
      $start_match = ($pstart == $act_start_AAR->[$i][$e]) ? 1 : 0;
      $stop_match  = ($pstop  == $act_stop_AAR->[$i][$e])  ? 1 : 0;
      if($start_match && $stop_match) { $found_both_match  = 1; }
      elsif($start_match)             { $found_start_match = 1; }
      elsif($stop_match)              { $found_stop_match = 1; }
    }
  }

  if   ($found_both_match)  { return (1, 1); }
  elsif($found_start_match) { return (1, 0); }
  elsif($found_stop_match)  { return (0, 1); }
  else                      { return (0, 0); }
}

# Subroutine: validateOriginSeq
#
# Synopsis:   Validate an origin sequence passed in
#             as <s> with --oseq <s>. It should have 
#             a single '|' in it, which occurs 
#             just before what should be the first nt
#             of the genome. Return the origin offset:
#             the number of nts before the "|".
#
#             For example: "TAATATT|AC"
#             indicates that the final 7 nts of each
#             genome should be "TAATAAT" and the first
#             two should be "AC". In this case the origin
#             offset is 7.
#
# Args:       $origin_seq: the origin sequence
#
# Returns:    Origin offset, as explained in synopsis, above.
#
sub validateOriginSeq {
  my $sub_name = "validateOriginSeq";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($origin_seq) = @_;

  my $origin_offset = index($origin_seq, "|");
  if($origin_offset == -1) { 
    die "ERROR with --oseq <s>, <s> must contain a single | character immediately before the nucleotide that should be the first nt of the genome";
  }
  my $second_offset = index($origin_seq, "|", $origin_offset+1);
  if($second_offset != -1) { 

    die "ERROR with --oseq <s>, <s> must contain a single | character, $origin_seq has more than one";
  }

  #printf("in $sub_name, $origin_seq returning $origin_offset\n");

  return $origin_offset;
}

# Subroutine: fetchCodon()
#
# Synopsis:   Fetch a codon given it's first position
#             and the strand and a Bio::Easel::SqFile object
#             that is the open sequence file with the desired
#             sequence.
#
# Args:       $sqfile:  Bio::Easel::SqFile object, open sequence
#                       file containing $seqname;
#             $seqname: name of sequence to fetch part of
#             $start:   start position of the codon
#             $strand:  strand we want ("+" or "-")
#
# Returns:    The codon as a string
#
sub fetchCodon {
  my $sub_name = "fetchCodon";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seqname, $start, $strand) = @_;

  my $codon_start = $start;
  my $codon_stop  = ($strand eq "-") ? $start - 2 : $start + 2; 

  my $newname = $seqname . "/" . $codon_start . "-" . $codon_stop;

  my @fetch_AA = ();
  push(@fetch_AA, [$newname, $codon_start, $codon_stop, $seqname]);

  my $faseq = $sqfile->fetch_subseqs(\@fetch_AA, -1);

  my ($header, $seq) = split("\n", $faseq);

  # printf("$faseq");
  
  return $seq;
}

# Subroutine: monocharacterString()
#
# Synopsis:   Return a string of length $len of repeated
#             instances of the character $char.
#
# Args:       $len:   desired length of the string to return
#             $char:  desired character
#
# Returns:    A string of $char repeated $len times.
#
sub monocharacterString {
  my $sub_name = "monocharacterString";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($len, $char) = @_;

  my $ret_str = "";
  for(my $i = 0; $i < $len; $i++) { 
    $ret_str .= $char;
  }

  return $ret_str;
}

# Subroutine: fetchHits()
#
# Synopsis:   Given 2D hashes that describe all hits, fetch
#             the hits for each CDS/exon to files.
#
# Args:       $sqfile:        Bio::Easel::SqFile object, open sequence
#                             file containing $seqname;
#             $do_skip:       '1' to skip alignment step after verifying alignments exist
#             $key:           string for naming output files (e.g.: "predicted", "corrected";
#             $mdl_order_AR:  ref to array of model names in order
#             $seq_order_AR:  ref to array of sequence names in order
#             $seqlen_HR:     ref to hash of total lengths
#             $start_HHR:     ref to 2D hash of start values, pre-filled
#             $stop_HHR:      ref to 2D hash of stop values, pre-filled
#             $strand_HHR:    ref to 2D hash of strand values, pre-filled
#             $out_aln_root:  root name for output files
#             $fafile_HR:     ref to hash: key: model name, value name of fasta file
#                             with fetched sequences, FILLED HERE
# 
# Returns:    void
#
sub fetchHits {
  my $sub_name = "fetchHits";
  my $nargs_exp = 11;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $do_skip, $key, $mdl_order_AR, $seq_order_AR, $seqlen_HR, $start_HHR, $stop_HHR, $strand_HHR, $out_aln_root, $fafile_HR) = @_;

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));

  $out_aln_root =~ s/\/.+$/\//; # remove everything after the first '/'

  printf("%-65s ... ", ($do_skip) ? "# Skipping fetch of $key exons" : "# Fetching $key exon sequences");
    
  foreach my $mdl (@{$mdl_order_AR}) { 
    my $out_key = $mdl;
    $out_key =~ s/ref/$key/;
    my @fetch_AA = ();
    my $nseq = 0;
    for(my $seq_i = 0; $seq_i < scalar(@{$seq_order_AR}); $seq_i++) { 
      my $seq = $seq_order_AR->[$seq_i];
      if($seq_i == 0 && (! exists $start_HHR->{$mdl}{$seq})) {
        die "ERROR in $sub_name(), no hit from model $mdl to the reference sequence $seq"; 
      }
      if(exists $start_HHR->{$mdl}{$seq}) { 
        my $accn = $seq;
        $accn =~ s/\:.+$//;
        my $newname .= $accn . "/" . $start_HHR->{$mdl}{$seq} . "-" . $stop_HHR->{$mdl}{$seq};
        my $start = $start_HHR->{$mdl}{$seq};
        my $stop  = $stop_HHR->{$mdl}{$seq};
        if($start < 0 || $stop < 0) { 
          # first, take care of off-by-one we introduced for coordinates that wrap around start (e.g. -2..3 in a length 
          # 10 genome is really 9..10..11..12..13 in duplicated genome, not 8..13)
          if($start < 0) { $start += 1; }
          if($stop  < 0) { $stop  += 1; }
          $start += $seqlen_HR->{$accn};
          $stop  += $seqlen_HR->{$accn};
        }
        push(@fetch_AA, [$newname, $start, $stop, $seq]);
        $nseq++;
      }
    }
    if($nseq > 0) { 
      my $cur_fafile = $out_aln_root . $out_key . ".fa";
      $sqfile->fetch_subseqs(\@fetch_AA, undef, $cur_fafile);
      # printf("Saved $nseq sequences to $cur_fafile.\n");
      $fafile_HR->{$mdl} = $cur_fafile;
    }
  }
  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));
  printf("done. [%s]\n", ($do_skip) ? "-skipaln" : sprintf("%.1f seconds", ($stop_time - $start_time)));

  return;
}

# Subroutine: alignHits()
#
# Synopsis:   Given a hash of sequence file names (already created)
#             to align for each model, align each file to the appropriate 
#             model to create a multiple alignment. Then parse the 
#             alignments to get info on all gaps in the alignment.
#
# Args:       $align:         path to hmmalign or cmalign executable
#             $fetch:         path to hmmfetch or cmfetch executable
#             $model_db:      model database file to fetch the models from
#             $do_skip:       '1' to skip alignment step after verifying alignments exist
#             $mdl_order_AR:  ref to array of model names in order
#             $seq_order_AR:  ref to array of sequence names in order
#             $fafile_HR:     ref to hash: key: model name, value: name of fasta file with 
#                             all fetched seqs to align to that model
#             $start_HHR:     ref to 2D hash of start values, pre-filled, 
#                             used only so we know what seqs should exist in alignments
#             $fid2ref_HHR:   ref to 2D hash of fractional identity values 
#                             of aligned sequences to the reference, FILLED HERE
#             $refdel_HHAR:   ref to 2D hash where value is an array, each element is
#                             an rf position that is deleted in the alignment of the $seq_accn
#                             FILLED HERE
#             $refins_HHAR:   ref to 2D hash where value is an array, each element is
#                             a string <rfpos>:<ct> where <rfpos> is a rf pos in the alignment
#                             after which $seq_accn has an insert and <ct> is the number of inserted
#                             positions.
#                             FILLED HERE
# 
# Returns:    void
#
sub alignHits {
  my $sub_name = "alignHits";
  my $nargs_exp = 11;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($align, $fetch, $model_db, $do_skip, $mdl_order_AR, $seq_order_AR, $fafile_HR, $start_HHR, $fid2ref_HHR, $refdel_HHAR, $refins_HHAR) = @_;

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));

  printf("%-65s ... ", ($do_skip) ? "# Skipping multiple alignment creation" : "# Creating multiple alignments");
    
  foreach my $mdl (@{$mdl_order_AR}) { 
    if(exists($fafile_HR->{$mdl})) { 
      my $cur_fafile  = $fafile_HR->{$mdl};
      my $cur_stkfile = $cur_fafile;
      $cur_stkfile =~ s/\.fa/\.stk/;
      my $cmd = "$fetch $model_db $mdl | $align - $cur_fafile > $cur_stkfile";
      #print $cmd . "\n";
      if($do_skip) { 
        if(! -s $cur_stkfile) { die "ERROR, with -skipaln, alignments are expected to already exist, but $cur_stkfile does not or is empty"; }
      }
      else { 
        runCommand($cmd, 0);
      }
      # printf("Saved $nseq aligned sequences to $cur_stkfile.\n");

      # store the fractional identities between each sequence and the reference
      # first we need to read in the MSA we just created 
      my $msa = Bio::Easel::MSA->new({
        fileLocation => $cur_stkfile,
                                     });  

      # determine which positions are RF positions
      my $rfseq = $msa->get_rf();
      my @rfseq_A = split("", $rfseq);
      my $alen  = $msa->alen();
      my @i_am_rfpos_A = (); # 0..$apos..$alen-1: '1' if RF is a nongap at position $apos+1, else '0'
      for(my $apos = 0; $apos < $alen; $apos++) { 
        $i_am_rfpos_A[$apos] = ($rfseq_A[$apos] eq ".") ? 0 : 1;
      }          

      my $i = 0; # this will remain '0', which is the reference sequence
      my $j = 0; # we'll increment this from 0..$nseq-1
      foreach my $seq (@{$seq_order_AR}) { 
        if(exists $start_HHR->{$mdl}{$seq}) { 
          $fid2ref_HHR->{$mdl}{$seq} = $msa->pairwise_identity($i, $j);
          # printf("storing percent id of $fid2ref_HHR->{$mdl}{$seq} for $mdl $seq\n"); 

          # determine the RF positions that are gaps in this sequence
          # and the positions of the inserted residues in this sequence
          # and store them
          my $aseqstring  = $msa->get_sqstring_aligned($j);
          my @aseq_A = split("", $aseqstring);
          my $rfpos = 0;
          for(my $apos = 0; $apos < $alen; $apos++) { 
            if($i_am_rfpos_A[$apos]) { 
              $rfpos++; 
            }
            if($aseq_A[$apos] =~ m/[\.\-]/) { # a gap in the sequence
              if($i_am_rfpos_A[$apos]) { # not a gap in the RF sequence
                # deletion (gap) relative to the reference sequence
                if(! exists $refdel_HHAR->{$mdl}{$seq}) { 
                  @{$refdel_HHAR->{$mdl}{$seq}} = (); 
                }
                updateGapArray(\@{$refdel_HHAR->{$mdl}{$seq}}, $rfpos, 1); # 1 informs the subroutine that this is a delete array
              }
            }
            else { # nongap in the sequence
              if(! $i_am_rfpos_A[$apos]) { # gap in the RF sequence
                # insertion in sequence relative to the reference sequence
                if(! exists $refins_HHAR->{$mdl}{$seq}) { 
                  @{$refins_HHAR->{$mdl}{$seq}} = (); 
                }
                updateGapArray(\@{$refins_HHAR->{$mdl}{$seq}}, $rfpos, 0); # 1 informs the subroutine that this is a delete array
              }
            }
          }
          $j++; # increment sequence index in msa
          # printf("printing insert info for $mdl $seq\n");
          # debugPrintGapArray(\@{$refins_HHAR->{$mdl}{$seq}});
          # printf("printing delete info for $mdl $seq\n");
          # debugPrintGapArray(\@{$refdel_HHAR->{$mdl}{$seq}});
        }
      }
    }
  }
  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));
  printf("done. [%s]\n", ($do_skip) ? "-skipaln" : sprintf("%.1f seconds", ($stop_time - $start_time)));

  return;
}

# Subroutine: checkForOverlaps()
#
# Synopsis:   Given refs to three arrays that describe 
#             a list of hits, check if any of them overlap.
#
# Args:       $name_AR:         ref to array of short names for each annotation
#             $start_AR:        ref to array of start positions
#             $stop_AR:         ref to array of stop positions
#             $strand_AR:       ref to array of strands
#             $expected_ol_AAR: ref to 2D array of expected overlaps $expected_ol_AAR->[$i][$j] is '1' if 
#                               those two exons are expected to overlap, PRE-FILLED, 
#                               can be undefined
#             $return_ol_AAR:   ref to 2D array of observed overlaps $observed_ol_AAR->[$i][$j] is '1' if 
#                               those two exons are observed to overlap here, FILLED HERE,
#                               can be undefined
# 
# Returns:    Two values:
#             $pass_fail_char: "P" if overlaps match those in $expected_ol_AAR, else "F"
#                              if $expected_ol_AAR is undef, always return "P"
#             $overlap_notes:  string describing the overlaps, empty string ("") if no overlaps.
#
sub checkForOverlaps {
  my $sub_name = "checkForOverlaps";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($name_AR, $start_AR, $stop_AR, $strand_AR, $expected_ol_AAR, $return_ol_AAR) = @_;

  my $ret_str = "";

  my @observed_ol_AA = (); # [$i][$j] is '1' if we observe an overlap between $i and $j (or $j and $i), else 0

  my $nhits = scalar(@{$name_AR});

  # initialize
  for(my $i = 0; $i < $nhits; $i++) { 
    for(my $j = 0; $j < $nhits; $j++) { 
      $observed_ol_AA[$i][$j] = 0;
    }
  }

  my $noverlaps = 0;
  for(my $i = 0; $i < $nhits; $i++) { 
    my $start_i  = $start_AR->[$i];
    my $stop_i   = $stop_AR->[$i];
    if($start_i > $stop_i) { 
      my $tmp  = $start_i;
      $start_i = $stop_i;
      $stop_i  = $tmp;
    }
    for(my $j = $i+1; $j < $nhits; $j++) { 
      my $start_j  = $start_AR->[$j];
      my $stop_j   = $stop_AR->[$j];
      if($start_j > $stop_j) { 
        my $tmp  = $start_j;
        $start_j = $stop_j;
        $stop_j  = $tmp;
      }
      my $nres_overlap = get_nres_overlap($start_i, $stop_i, $start_j, $stop_j);
      if($nres_overlap > 0) { 
        $noverlaps++;
        $observed_ol_AA[$i][$j] = 1;
        $observed_ol_AA[$j][$i] = 1;
        # $ret_str .= sprintf("%s overlaps with %s (%s);", $name_AR->[$i], $name_AR->[$j], ($strand_AR->[$i] eq $strand_AR->[$j]) ? "same strand" : "opposite strands");
        if($ret_str ne "") { 
          $ret_str .= " ";
        }
        $ret_str .= sprintf("%s/%s", $name_AR->[$i], $name_AR->[$j]); 
      }
      else { # no overlap
        $observed_ol_AA[$i][$j] = 0;
        $observed_ol_AA[$j][$i] = 0;
      }
    }
  }

  # check @observed_ol_AA against @{$expected_ol_AAR} and/or
  # copy @observed_ol_AA to @{$return_ol_AAR}
  my $pass_fail_char = "P";
  for(my $i = 0; $i < $nhits; $i++) { 
    for(my $j = 0; $j < $nhits; $j++) { 
      if(defined $expected_ol_AAR) { 
        if($observed_ol_AA[$i][$j] ne $expected_ol_AAR->[$i][$j]) { 
          $pass_fail_char = "F";
        }
      }
      if(defined $return_ol_AAR) {
        $return_ol_AAR->[$i][$j] = $observed_ol_AA[$i][$j];
      }
    }
  }

  if($ret_str ne "") { 
    $ret_str = $pass_fail_char . " " . $noverlaps . " " . $ret_str;
  }

  return ($pass_fail_char, $ret_str);
}

# Subroutine: get_nres_overlap()
# Args:       $start1: start position of hit 1 (must be <= $end1)
#             $end1:   end   position of hit 1 (must be >= $end1)
#             $start2: start position of hit 2 (must be <= $end2)
#             $end2:   end   position of hit 2 (must be >= $end2)
#
# Returns:    Number of residues of overlap between hit1 and hit2,
#             0 if none
# Dies:       if $end1 < $start1 or $end2 < $start2.

sub get_nres_overlap {
  if(scalar(@_) != 4) { die "ERROR get_nres_overlap() entered with wrong number of input args"; }

  my ($start1, $end1, $start2, $end2) = @_; 

  #printf("in get_nres_overlap $start1..$end1 $start2..$end2\n");

  if($start1 > $end1) { die "ERROR start1 > end1 ($start1 > $end1) in get_nres_overlap()"; }
  if($start2 > $end2) { die "ERROR start2 > end2 ($start2 > $end2) in get_nres_overlap()"; }

  # Given: $start1 <= $end1 and $start2 <= $end2.
  
  # Swap if nec so that $start1 <= $start2.
  if($start1 > $start2) { 
    my $tmp;
    $tmp   = $start1; $start1 = $start2; $start2 = $tmp;
    $tmp   =   $end1;   $end1 =   $end2;   $end2 = $tmp;
  }
  
  # 3 possible cases:
  # Case 1. $start1 <=   $end1 <  $start2 <=   $end2  Overlap is 0
  # Case 2. $start1 <= $start2 <=   $end1 <    $end2  
  # Case 3. $start1 <= $start2 <=   $end2 <=   $end1
  if($end1 < $start2) { return 0; }                      # case 1
  if($end1 <   $end2) { return ($end1 - $start2 + 1); }  # case 2
  if($end2 <=  $end1) { return ($end2 - $start2 + 1); }  # case 3
  die "Unforeseen case in get_nres_overlap $start1..$end1 and $start2..$end2";

  return; # NOT REACHED
}

# Subroutine: updateGapArray()
#
# Synopsis:   Given an rfpos that we have a gap in, update a 
#             array that stores all the gap information for a given
#             sequence. This can be called for an array that holds
#             gaps in the reference (deletions, refdel* data 
#             structures) or in the other sequence (insertions,
#             refins* data structures).
#
# Args:       $AR:            reference to the array to update
#             $rfpos:         reference position the gap occurs at or after
#             $is_delete:     '1' if we're updating a delete array, else we're
#                             updating an insert array (we do the update slightly
#                             differently for each type)
#
# Returns:    void
#
sub updateGapArray {
  my $sub_name = "updateGapArray";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($AR, $rfpos, $is_delete) = @_;

  my $nel = scalar(@{$AR});
  my $same_as_prv = 0;
  if($nel > 0) { 
    # need to check latest element to see if it pertains to the same rfpos
    my ($prv_rfpos, $prv_cnt) = split(":", $AR->[$nel-1]);
    # printf("\tsplit %s into %s and %d\n", $AR->[$nel-1], $prv_rfpos, $prv_cnt);
    if((  $is_delete) && (($prv_rfpos + $prv_cnt) == $rfpos) ||
       (! $is_delete) && ($prv_rfpos == $rfpos)) { 
      # this is not the first insert/delete at this position, update previous value
      $same_as_prv = 1;
      $AR->[($nel-1)] = $prv_rfpos . ":" . ($prv_cnt+1);
    }
  }
  if(! $same_as_prv) { 
    # either this is the first insert/delete or not the same as the previous one
    push(@{$AR}, $rfpos . ":" . "1");
  }

  return;
}


# Subroutine: findSpecialGap()
#
# Synopsis:   Given a gap string that lists all gaps in a CDS, determine if the gaps cause
#             the predicted length of the CDS to be non-modulo 3. If so, determine if there's
#             exactly one gap that if we remove it, the remaining gaps will make the predicted
#             length of the CDS to be modulo 3. 
#
#             I think there can only be 0 or 1 such gap, but I'm not sure.
# 
#             Return two values: 
#             - the length of the speical gap if there is one, or 0 if there is not
#               (this value will be negative if it is a deletion relative to the 
#               reference, and positive if it is an insertion relative to the reference)
#             - the string that describes the special gap, 
#               or '-' if the predicted length of the CDS is modulo 3 if we include all gaps
#               or '?' if the predicted length of the CDS is not modulo 3 but there is no
#               special gap that can explain the non-modulo-3ness.
#
#             If the net gap length is zero modulo 3, return (0, "-").
#             If the net gap length is nonzero modulo 3, and there's exactly 1 gap that equals the net length modulo 3, return its length and the string..
#             If the net gap length is nonzero modulo 3, and there's not exactly 1 gap that equals the net length modulo 3, return (0, "?");
#            
#             Examples: 
#
#               Input:    $gapstr = "D377:2,I388:2,I1066:1 net gap is 1
#               Returns:  (1, "I1066:1")
#
#               Input:    $gapstr = "D377:2,I388:2", , net gap is 0;
#               Returns:  (0, "-")
#
#               Input:    $gapstr = "D377:5,I388:2", net gap is -3;
#               Returns:  (0, "-")
#
#               Input:    $gapstr = "D377:2,I388:5", net gap is 3;
#               Returns:  (0, "-")
#
#               Input:    $gapstr = "D377:2,I388:3,I1066:1", net gap is 2;
#               Returns:  (0, "?")
#
#               Input:    $gapstr = "D277:2,D377:2,I388:2", $net_gap = 2;
#               Returns:  (0, "?")
#
sub findSpecialGap {
  my $sub_name = "findSpecialGap";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($gapstr) = (@_);
  # printf("in $sub_name gapstr: $gapstr\n");

  if($gapstr eq "-") { # no gaps
    return ("-"); 
  }

  # determine net number of gaps if we include them all
  my @el_A = split(",", $gapstr);
  my $ngaps = scalar(@el_A);
  my @gapsize_A = (); # [0..$ngaps-1], length of each gap
  my $tot_gaps = 0;

  for(my $e = 0; $e < $ngaps; $e++) { 
    my $el = $el_A[$e];
    my ($type_loc, $len) = split(":", $el);

    if(($len % 3) == 0) { 
      $gapsize_A[$e] = 0;
    }
    else { 
      if($type_loc =~ /^I(\d+)/) { # insert
        $gapsize_A[$e] = $len;
        $tot_gaps += $len;
      }
      elsif($type_loc =~ /^D(\d+)/) { # delete
        $gapsize_A[$e] = -1 * $len;
        $tot_gaps -= $len;
      }
      else { 
        die "ERROR: in $sub_name, unable to parse gap string: $gapstr element $el"; 
      }
    }
    if(($tot_gaps % 3) == 0) { 
      return ("-"); 
    }
  }

  # find gap that gives minimal net gap length
  
  my $min_len = undef;
  for(my $e = 0; $e < $ngaps; $e++) { 
    my $net_gap_len = $tot_gaps - $gapsize_A[$e];
    if(($net_gap_len % 3) == 0) { 
      if(! defined $min_len || ($net_gap_len < $min_len)) { 
        $min_len = $net_gap_len;
      }
    }
  }

  my $nspecial = 0;
  my $special_str = "";
  my $special_len = 0;

  if(! defined $min_len) { # no gaps can be removed and give modulo 3
    return ("?");
  }
  else { # at least 1 gap can be removed and give modulo 3
    for(my $e = 0; $e < $ngaps; $e++) { 
      my $net_gap_len = $tot_gaps - $gapsize_A[$e];
      if(($net_gap_len % 3) == 0) { 
        if($net_gap_len == $min_len) { 
          if($special_str ne "") { $special_str .= ","; }
          $special_str .= $el_A[$e];
        }
      }
    }
    # printf("returning special_str: $special_str\n");
    return $special_str;
  }
}

# Subroutine: splitFastaFile()
#
# Synopsis:   Split up a fasta file into <n> smaller files using 
#             the esl-ssplit script.
#
# Args:       $esl_ssplit:   path to perl esl-sscript.pl script to use
#             $fafile:       fasta file to split up
#             $nfiles:       number of files to split $fafile into
#
# Returns:    void
# Dies:       if command fails or there is some other problem
#
sub splitFastaFile {
  my $sub_name = "splitFastaFile";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($esl_ssplit, $fafile, $nfiles) = @_;

  if(! -s $esl_ssplit) { die "ERROR the $esl_ssplit file does not exist or is empty"; }
  if(! -x $esl_ssplit) { die "ERROR the $esl_ssplit file is not executable"; }

  my $cmd = "$esl_ssplit -n $fafile $nfiles > /dev/null";
  runCommand($cmd, 0);

  return;
}

# Subroutine: wrapperCombineExonsIntoCDS()
#
# Synopsis:   For all CDS, combine all exons into CDS. A wrapper function
#             for combineExonsIntoCDS().
#
# Args:       $nhmm:            total number of exons/models
#             $dir:             directory for output files
#             $key:             string for naming output files (e.g.: "predicted" or "corrected")
#             $model_AR:        ref to array of model names 
#             $hmm2cds_map_AR:  ref to array mapping models to CDS
#             $hmm_is_first_AR: ref to array signifying if a model is first in a CDS or not
#             $hmm_is_final_AR: ref to array signifying if a model is final in a CDS or not
#             $outfile_AR:      ref to array of output files, filled here
#
# Returns:    void
# Dies:       if something unexpected happens when reading the exon fasta files
#
sub wrapperCombineExonsIntoCDS {
  my $sub_name = "combineExonsIntoCDS";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($nhmm, $dir, $key, $model_AR, $hmm2cds_map_AR, $hmm_is_first_AR, $hmm_is_final_AR, $outfile_AR) = @_;

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));
  printf("%-65s ... ", "# Combining multi-exon $key CDS ");

  my @tmp_exon_fafile_A = (); # temporary array of exon fafiles for all exons in current CDS
  
  for(my $h = 0; $h < $nhmm; $h++) { 
    my $cds_i      = $hmm2cds_map_AR->[$h];
    my $cur_fafile = $dir . "/" . $model_AR->[$h] . ".fa";
    $cur_fafile    =~ s/ref/$key/;
    push(@tmp_exon_fafile_A, $cur_fafile);
    if($hmm_is_final_A[$h]) { 
      if(! $hmm_is_first_A[$h]) { 
        $cur_fafile =~ s/\.exon\.\d+//; # remove exon.<d> part of file name
        combineExonsIntoCDS(\@tmp_exon_fafile_A, $cur_fafile);
      }
      else { # a single exon gene, we should already have the sequence from alignHits
        if(! -s $cur_fafile) { die sprintf("ERROR, expected output fasta file for CDS %s does not exist: $cur_fafile", $cds_i+1); }
      }
      push(@{$outfile_AR}, $cur_fafile);
      @tmp_exon_fafile_A = (); # reset to empty, for next CDS
    }
  }

  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));
  printf("done. [%.1f seconds]\n", ($stop_time - $start_time));

  return;
}
# Subroutine: combineExonsIntoCDS()
#
# Synopsis:   Given an array of fasta files each with a different exon 
#             of a CDS, create a single new fasta file that has the complete
#             CDS by combining all the exon sequences together.
#
# Args:       $exon_fafile_AR: ref to array of all exon files for this CDS
#             $cds_fafile:     name of CDS fasta file to create
#
# Returns:    void
# Dies:       if something unexpected happens when reading the exon fasta files
#
sub combineExonsIntoCDS {
  my $sub_name = "combineExonsIntoCDS";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($exon_fafile_AR, $cds_fafile) = @_;

  # we do this manually, by opening each exon file into a separate file handle and then 
  # outputting each exon as we read it
  my $nexons = scalar(@{$exon_fafile_AR});
  my @FH_A = ();
  my $FH;
  for(my $i = 0; $i < $nexons; $i++) { 
    open($FH_A[$i], "<", $exon_fafile_AR->[$i]) || die "ERROR unable to open $exon_fafile_AR->[$i] for reading";
  }

  open(OUT, ">", $cds_fafile) || die "ERROR unable to open file $cds_fafile for writing";

  my $nseqs_read      = 0;  # number of sequences read so far
  my @read_header_A   = (); # current header line read from file [0..$nexons-1]
  my @read_sequence_A = (); # current sequence read from file [0..$nexons-1]
  my $header          = ""; # name of sequence to output

  # read all exons of each sequence and then output it
  while(($nseqs_read == 0) || (defined $read_header_A[0])) { 
    # first read the header line for each exon, if this is the first sequence, else we already read them at the end of this loop
    if($nseqs_read == 0) { 
      for(my $i = 0; $i < $nexons; $i++) { 
        my $FH = $FH_A[$i];
        $read_header_A[$i] = <$FH>;
        chomp $read_header_A[$i];
      }
    }
    # split up header lines to get new header name
    my @read_seqname_A = ();
    my @read_coords_A = ();
    for(my $i = 0; $i < $nexons; $i++) { 
      ($read_seqname_A[$i], $read_coords_A[$i]) = split(/\//, $read_header_A[$i]);
      if(($i > 0) && ($read_seqname_A[$i] ne $read_seqname_A[0])) { 
        die "ERROR, expected same sequence names for each exon, but read $read_seqname_A[0] and $read_seqname_A[$i]"; 
      }
      # add to growing header
      $header .= ($i == 0) ? 
          ($read_seqname_A[$i] . "/" . $read_coords_A[$i]) : 
          ("," . $read_coords_A[$i]);
    }
    # printf("nseqs_read: $nseqs_read\nheader: $header\n");

    # we've just read the header lines for each sequence, now read all sequence lines
    for(my $i = 0; $i < $nexons; $i++) { 
      if(! defined $read_header_A[$i]) { die "ERROR in $sub_name, ran out of header lines too soon for file $i"; }
      my $FH = $FH_A[$i];
      my $read_line = <$FH>;
      while((defined $read_line) && ($read_line !~ m/^\>/)) { 
        $read_sequence_A[$i] .= $read_line;
        $read_line = <$FH>;
      }
      # if(defined $read_line) { printf("setting read_header_A[$i] to $read_line\n"); }
      $read_header_A[$i] = $read_line;
      if(defined $read_header_A[$i]) { chomp($read_header_A[$i]); }
    }

    # output full CDS for the current accesion
    print OUT $header . "\n";
    for(my $i = 0; $i < $nexons; $i++) { 
      print OUT $read_sequence_A[$i];
      $read_sequence_A[$i] = "";
    }
    $nseqs_read++;
    $header = "";
  }

  close(OUT);

  return;
}

# Subroutine: parseEslTranslateOutput()
#
# Synopsis:   Given an output file from EslTranslateOutput(), keep only
#             full length translations. 
#
# Args:       $esl_translate_output: output fasta file from esl-translate
#             $new_output:           new output file to create with subset
#                                    of sequences from $esl_translate_output
#             $prot_must_start_at_1: '1' if translated protein sequence must start 
#                                    at position 1 of the CDS sequence
#             $prot_must_stop_at_L:  'L' if translated protein sequence must stop
#                                    at final position of the CDS sequence
#             $prot_HR:              ref to hash to fill here, keys are seq accessions
#                                    values are number of protein sequences fetch for that
#                                    accession, any value >1 is an error, and we'll die
#             $nprot_R:              ref to scalar to fill with number of accessions
#                                    for which we fetched a full length protein
#
# Returns:    void
# Dies:       if something unexpected happens when reading the exon fasta files
#
sub parseEslTranslateOutput {
  my $sub_name = "parseEslTranslateOutput";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my $nprot = 0; # number of full length proteins kept

  my ($esl_translate_output, $new_output, $prot_must_start_at_1, $prot_must_stop_at_L, $prot_HR, $nprot_R) = @_;

  if(! $prot_must_start_at_1) { die "ERROR in $sub_name, prot_must_start_at_1 is FALSE, you haven't implemented this case yet..."; }

  open(IN,  "<", $esl_translate_output) || die "ERROR unable to open $esl_translate_output for reading";
  open(OUT, ">", $new_output)           || die "ERROR unable to open $new_output for writing";

  my ($source_name, $source_coords, $coords_from, $coords_to, $length, $frame);
  my $source_length = 0; # length of original CDS sequence, including stop
  my $coords_length = 0; # length of original CDS sequence that was translated, not included stop codon
  my $print_flag = 0;
  my $line = <IN>;
  while(defined $line) { 
    if($line =~ m/^\>/) { 
      parseEslTranslateDefline($line, \$coords_from, \$coords_to, \$source_name, \$source_coords, \$source_length);

      $coords_length = ($coords_to - $coords_from) + 1;
      # now determine if this protein passes or not
      $print_flag = 0; # by default it fails, but we redefine this below if it passes
      if($prot_must_start_at_1 && $prot_must_stop_at_L && 
         ($coords_from == 1) && ($coords_to == ($source_length - 3))) { 
        # we're looking only for full length proteins, and this is one
        $print_flag = 1;
      }
      if($prot_must_start_at_1 && (! $prot_must_stop_at_L) && 
         ($coords_from == 1) && $coords_to != ($source_length - 3)) { 
        # we're NOT looking for full length proteins, but rather for proteins
        # that start at the beginning of the CDS but do NOT end at the end of the 
        # CDS, they may be shorter or longer
        $print_flag = 1;
      }
      
      if($print_flag) { # print out new sequence name
        # print OUT $line; 
        printf OUT (">%s/%s-translated\n", $source_name, $source_coords);
        $nprot++;
        $prot_HR->{$source_name}++;
        if($prot_HR->{$source_name} > 1) { 
          die "ERROR in $sub_name, found more than 1 full length protein translation for $source_name"; 
        }
      }
      $line = <IN>;
      while(defined $line && $line !~ m/^\>/) { 
        if($print_flag) { print OUT $line; }
        $line = <IN>;
      }
    }
    else { 
      die "ERROR in $sub_name, expected header line but read line: $line";
    }
  }
  close(OUT);

  $$nprot_R = $nprot;
  return;
}


# Subroutine: parseEslTranslateDefline()
#
# Synopsis:   Given a defline output from esl-translate, parse it and
#             return some relevant info.
#
# Args:       $line:                 the line to parse
#             $coords_from_R:        ref to coordinate of start nt of translation, filled here, can be undef
#             $coords_to_R:          ref to coordinate of stop  nt of translation, filled here, can be undef
#             $source_name_R:        ref to name of source CDS, filled here, can be undef
#             $source_coords_R:      ref to coordinates part of name of source CDS, filled here, can be undef
#             $source_len_R:         ref to coordinate of length of source CDS, filled here, can be undef
#
# Returns:    void
# Dies:       if format of defline is unexpected
#
sub parseEslTranslateDefline {
  my $sub_name = "parseEslTranslateDefline";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($line, $coords_from_R, $coords_to_R, $source_name_R, $source_coords_R, $source_len_R) = @_;

  my ($source_name, $source_coords, $coords_from, $coords_to, $length, $frame, $source_len);

  # example $line's:
  #>orf1 source=NC_001346/150-455 coords=1..303 length=101 frame=1  
  # OR
  #>orf1 source=NC_001346/2527-1886,1793-1353 coords=87..161 length=25 frame=3  
  if($line =~ /^\>orf\d+\s+source\=(\S+)\/(\S+)\s+coords\=(\d+)\.\.(\d+)\s+length\=(\d+)\s+frame=(\d+)/) { 
    ($source_name, $source_coords, $coords_from, $coords_to, $length, $frame) = ($1, $2, $3, $4, $5, $6);
    # $source_coords = '150-455' or
    # $source_coords = '2527-1886,1793-1353'
    my @source_coords_A = split(",", $source_coords);
    $source_len = 0;
    foreach my $cur_coords (@source_coords_A) { 
      if($cur_coords =~ m/^(\-?\d+)\-(\-?\d+)$/) { 
        my($cur_from, $cur_to) = ($1, $2);
        $source_len += ($cur_from < $cur_to) ? ($cur_to - $cur_from + 1) : ($cur_from - $cur_to + 1);
        if(($cur_from < 0) && ($cur_to > 0)) { $source_len--; } # fix off-by-one introduced by negative indexing in sequence position
        if(($cur_from > 0) && ($cur_to < 0)) { $source_len--; } # fix off-by-one introduced by negative indexing in sequence position
      }
      else { 
        die "ERROR in $sub_name, unable to parse source coordinates in header line: $line";
      }
    }
  }
  else { 
    die "ERROR in $sub_name, unable to parse header line: $line";
  }

  if(defined $coords_from_R)   { $$coords_from_R   = $coords_from; }
  if(defined $coords_to_R)     { $$coords_to_R     = $coords_to;   }
  if(defined $source_name_R)   { $$source_name_R   = $source_name; }
  if(defined $source_coords_R) { $$source_coords_R = $source_coords; }
  if(defined $source_len_R)    { $$source_len_R    = $source_len;  }
  return;
}

######################
# OUTPUT subroutines #
######################
#
# Subroutine: outputSeqRowHeadings()
#
# Synopsis:   For 'sequences are rows' tabular output, output the headings.
#
# Args:       $do_nofid:           '1' if we're not printing fractional ids, else '0'
#             $do_mdlb:            '1' if we're not printing model boundaries, else '0'
#             $do_noss3:           '1' if we're not printing SS3 columns, else '0'
#             $do_nostop:          '1' if we're not printing Stop codons, else '0'
#             $origin_seq:         origin sequence, or undef if ! defined
#             $ref_tot_nexons:     number of total exons in reference
#             $nhmm:               number of total HMMs
#             $hmm2cds_map_AR:     ref to @hmm2cds_map array
#             $hmm2exon_map_AR:    ref to @hmm2exon_map array
#             $hmm_is_final_AR:    ref to @hmm_is_final_A
#             $hmm_is_first_AR:    ref to @hmm_is_first_A
#             $cds_out_short_AR:   ref to @cds_out_short_A
#             $cds_out_product_AR: ref to @cds_out_short_A
#
sub outputSeqRowHeadings {
  my $sub_name = "outputSeqRowHeadings";
  my $nargs_exp = 13;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($do_nofid, $do_mdlb, $do_noss3, $do_nostop, $origin_seq, $ref_tot_nexons, $nhmm, $hmm2cds_map_AR, $hmm2exon_map_AR, $hmm_is_first_AR, $hmm_is_final_AR, $cds_out_short_AR, $cds_out_product_AR) = @_;

  my $width;  # width of a field
  my $pad;    # string of all spaces used for pretty formatting
  my $width_result = 5 + $ref_tot_nexons + 2;

  # line 1 of column headers
  printf("%-20s  %6s", "#", "");
  if(defined $origin_seq) { 
    printf("  %22s", "");
  }
  # for each CDS, output the topmost column header
  $width = 0;
  for(my $h = 0; $h < $nhmm; $h++) { 
    $width += 18;
    my $cds_i = $hmm2cds_map_AR->[$h];
    if(! $do_nofid)  { $width += 6; }
    if(! $do_nomdlb) { $width += 4; }
    if($hmm_is_final_AR->[$h]) { 
      $width += 7;
      if(! $do_noss3)  { $width += 4; }
      if(! $do_nostop) { $width += 4; }
      printf("    %*s", $width, $cds_out_short_AR->[$cds_i] . monocharacterString(($width-length($cds_out_short_AR->[$cds_i]))/2, " "));
      $width = 0;
    }
  }
  printf("  %6s", "");
  printf("  %5s", "");
  if(! $do_noexist) { 
    printf("    %19s", "");
  }
  printf("  %*s", $width_result, "");
  printf("\n");
  
  # line 2 of column headers
  printf("%-20s  %6s", "#", "");
  if(defined $origin_seq) { 
    printf("  %22s", "   origin sequence");
  }
  # for each CDS, output the second column header
  $width = 0;
  for(my $h = 0; $h < $nhmm; $h++) { 
    $pad = "";
    $width += 18;
    my $cds_i = $hmm2cds_map_AR->[$h];
    if(! $do_nofid)  { $width += 6; }
    if(! $do_nomdlb) { $width += 4; }
    if($hmm_is_final_AR->[$h]) { 
      $width += 7;
      if(! $do_noss3)  { $width += 4; }
      if(! $do_nostop) { $width += 4; }
      printf("    %*s", $width, substr($cds_out_product_AR->[$cds_i], 0, $width) . monocharacterString(($width-length($cds_out_product_AR->[$cds_i]))/2, " "));
      $width = 0;
    }
    else { 
      $width += 2;
    }
  }
  printf("  %6s", "");
  printf("  %5s", "");
  if(! $do_noexist) { 
    printf(" %19s", "existing annotation");
  }
  printf("  %*s", $width_result, "");
  printf("\n");
  
  # line 3 of column headers 
  printf("%-20s  %6s", "#", "");
  if(defined $origin_seq) { 
    printf("  %22s", "----------------------");
  }
  $width = 0;
  for(my $h = 0; $h < $nhmm; $h++) { 
    $width += 18;
    if(! $do_nofid)  { $width += 6; }
    if(! $do_nomdlb) { $width += 4; }
    if($hmm_is_final_AR->[$h]) { 
      $width += 9;
      if(! $do_noss3)  { $width += 4; }
      if(! $do_nostop) { $width += 4; }
      printf("  %s", monocharacterString($width, "-"));
      $width = 0;
    }
    else { 
      $width += 1;
    }
  }
#    printf("  %6s  %-*s", "", $width_result, "");
  printf("  %6s", "");
  printf("  %5s", "");
  if(! $do_noexist) { 
    printf("  %19s", "-------------------");
  }
  printf("  %-*s", $width_result, "");
  printf("\n");
  
  # line 4 of column headers
  printf("%-20s  %6s", "# accession", "totlen");
  if(defined $origin_seq) {
    printf(" %2s %5s %5s %5s %2s", " #", "start", "stop", "offst", "PF");
  }
  for(my $h = 0; $h < $nhmm; $h++) { 
    printf("  %8s %8s", 
           sprintf("%s%s", "start", $hmm2exon_map_AR->[$h]+1), 
           sprintf("%s%s", "stop",  $hmm2exon_map_AR->[$h]+1));
    if(! $do_nofid) { 
      printf(" %5s", sprintf("%s%s", "fid", $hmm2exon_map_AR->[$h]+1));
    }
    if(! $do_nomdlb) { 
      printf(" %3s", sprintf("%s%s", "md", $hmm2exon_map_AR->[$h]+1));
    }
    if($hmm_is_final_AR->[$h]) { 
      printf(" %6s", "length");
      if(! $do_noss3) { 
        printf(" %3s", "SS3");
      }
      if(! $do_nostop) { 
        printf(" %3s", "stp");
      }
      printf(" %2s", "PF");
    }
  }
  printf("  %6s", "totlen");
  printf("  %5s", "avgid");    
  if(! $do_noexist) { 
    printf("  %5s  %5s  %5s", "cds", "exons", "match");
  }
  if(! $do_noolap) { 
    printf("  %20s", " overlaps?");
  }
  
  printf("  %-*s", $width_result, "result");
  print "\n";
  
  # line 5 of column headers
  printf("%-20s  %6s", "#-------------------", "------");
  if(defined $origin_seq) {
    printf(" %2s %5s %5s %5s %2s", "--", "-----", "-----", "-----", "--");
  }
  for(my $h = 0; $h < $nhmm; $h++) { 
    printf("  %8s %8s", "--------", "--------");
    if(! $do_nofid) { 
      printf(" %5s", "-----");
    }
    if(! $do_nomdlb) { 
      printf(" %3s", "---");
    }
    if($hmm_is_final_AR->[$h]) { 
      printf(" %6s", "------");
      if(! $do_noss3) { 
        printf(" %3s", "---");
      }
      if(! $do_nostop) { 
        printf(" %3s", "---");
      }
      printf(" --");
    }
  }
  printf("  %6s", "------");
  printf("  %5s", "-----");
  
  if(! $do_noexist) { 
    printf("  %5s  %5s  %5s", "-----", "-----", "-----");
  }
  if(! $do_noolap) { 
    printf("  %20s", monocharacterString(20, "-"));
  }
  printf("  %-*s", $width_result, monocharacterString($width_result, "-"));
  
  print "\n";
}

# Subroutine: printColumnHeaderExplanations()
# Args:       $do_oseq: '1' if -oseq was enabled
#
# Returns:    void

sub printColumnHeaderExplanations {
  my $sub_name = "printColumnHeaderExplanations";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($do_oseq, $do_nomdlb, $do_noexist, $do_nobrack, $do_nostop, $do_nofid, $do_noss3, $do_noolap) = @_; 
  
  print("#\n");
  print("# Explanations of column headings (in left to right order):\n");

  my $width = 35;

  printf("# %-*s %s\n", $width, "\"accession\":", "GenBank accession for genomic sequence");
  printf("# %-*s %s\n", $width, "\"totlen\":",    "total length (nt) for accession");

  if($do_oseq) {
    printf("#\n");
    printf("# %-*s %s\n", $width, "\"origin sequence: #\":",      "number of occurences of origin sequence (input with -oseq) in genome");
    printf("# %-*s %s\n", $width, "\"origin sequence: start\":",  "start position of lone occurence of origin sequence (if only 1 exists)");
    printf("# %-*s %s\n", $width, "\"origin sequence: stop\":",   "stop  position of lone occurence of origin sequence (if only 1 exists)");
    printf("# %-*s %s\n", $width, "\"origin sequence: offst\":",  "predicted offset of genome, number of nucleotides to shift start (>0: clockwise; <0: counterclockwise)");
    printf("# %-*s %s\n", $width, "\"origin sequence: PF\":",     "'P' (for PASS) if there is exactly 1 occurence of the offset, else 'F' for FAIL");
  }

  printf("#\n");
  printf("# %-*s %s%s\n", $width, "\"CDS #<i>: start<j>\":", "start position of exon #<j> of CDS #<i>", ($do_nobrack) ? "" : "enclosed in brackets \"\[\]\" if different from all exon starts in existing annotation");
  printf("# %-*s %s%s\n", $width, "\"CDS #<i>: stop<j>\":",  "stop  position of exon #<j> of CDS #<i>", ($do_nobrack) ? "" : "enclosed in brackets \"\[\]\" if different from all exon starts in existing annotation");

  if(! $do_nofid) { 
    printf("# %-*s %s\n", $width, "\"CDS #<i>: fid<j>\":",  "fractional identity between exon #<j> of CDS #<i> and reference genome");
  }

  if(! $do_nomdlb) { 
    printf("# %-*s %s\n", $width, "\"CDS #<i>: md<j>\":",  "annotation indicating if alignment to reference extends to 5' and 3' end of reference exon.");
    printf("# %-*s %s\n", $width, "",                      "first character pertains to 5' end and second character pertains to 3' end.");
    printf("# %-*s %s\n", $width, "",                      "possible values for each of the two characters:");
    printf("# %-*s %s\n", $width, "",                      "  \".\":   alignment extends to boundary of reference");
    printf("# %-*s %s\n", $width, "",                      "  \"<d>\": alignment truncates <d> nucleotides short of boundary of reference (1 <= <d> <= 9)");
    printf("# %-*s %s\n", $width, "",                      "  \"+\":   alignment truncates >= 10 nucleotides short of boundary of reference");
  }

  printf("# %-*s %s\n", $width, "\"CDS #<i>: length\":",   "length of CDS #<i> (all exons summed)");

  if(! $do_noss3) { 
    print("#\n");
    printf("# %-*s %s\n", $width, "\"CDS #<i>: SS3\":",   "annotation indicating if predicted CDS has a valid start codon, stop codon and is a multiple of 3");
    printf("# %-*s %s\n", $width, "",                      "first  character: '.' if predicted CDS has a valid start codon, else '!'");
    printf("# %-*s %s\n", $width, "",                      "second character: '.' if predicted CDS has a valid stop  codon, else '!'");
    printf("# %-*s %s\n", $width, "",                      "third  character: '.' if predicted CDS has a length which is a multiple of three, else '!'");
  }

  if(! $do_nostop) { 
    printf("# %-*s %s\n", $width, "\"CDS #<i>: stp\":",   "the predicted stop codon for this CDS");
  }

  printf("# %-*s %s\n", $width, "\"CDS #<i>: PF\":",      "annotation indicating if this exon PASSED ('P') or FAILED ('F')");
  printf("# %-*s %s\n", $width, "",                       "  a CDS PASSES ('P') if and only if ALL of its exons have a valid start codon,");
  printf("# %-*s %s\n", $width, "",                       "  a valid stop codon, are lengths that are a multiple of 3, and have an");
  printf("# %-*s %s\n", $width, "",                       "  alignment to the corresponding reference exon that extends to the 5'");
  printf("# %-*s %s\n", $width, "",                       "  and 3' boundary of the reference annotation.");
  printf("# %-*s %s\n", $width, "",                       "  If >= 1 of these conditions is not met then the CDS FAILS ('F').");

  print("#\n");
  printf("# %-*s %s\n", $width, "\"totlen\":",            "total length (nt) for accession (repeated for convenience)"); 
  
  if(! $do_noexist) { 
    printf("#\n");
    printf("# %-*s %s\n", $width, "\"existing annotation: cds\"",   "number of CDS in the existing NCBI annotation for this accession");
    printf("# %-*s %s\n", $width, "\"existing annotation: exons\"", "total number of exons in the existing NCBI annotation for this accession");
    printf("# %-*s %s\n", $width, "\"existing annotation: match\"", "number of exons in existing NCBI annotation for which existing and predicted annotation agree exactly");
  }

  if(! $do_noolap) { 
    printf("#\n");
    printf("# %-*s %s\n", $width, "\"overlaps\?\"",   "text describing which (if any) of the predicted exons overlap with each other");
    printf("# %-*s %s\n", $width, "",                 "first character:   'P' for PASS if predicted annotation for this accession has same overlaps as the reference");
    printf("# %-*s %s\n", $width, "",                 "                   'F' for FAIL if it does not");
    printf("# %-*s %s\n", $width, "",                 "second character:  number of overlaps between any two exons");
    printf("# %-*s %s\n", $width, "",                 "remainder of line: text explaining which exons overlap");
    printf("# %-*s %s\n", $width, "",                 "  e.g.: \"3.2/4.1\" indicates exon #2 of CDS #3 overlaps with exon #1 of CDS #4 on either strand");
  }  

  print("#\n");
  printf("# %-*s %s\n", $width, "\"result\":",            "\"PASS\" or \"FAIL\". \"PASS\" if and only if all tests for this accession PASSED ('P')");
  printf("# %-*s %s\n", $width, "",                       "as indicated in the \"PF\" columns. Followed by the individual P/F results in order.");
  if($do_noolap) { 
    printf("# %-*s %s\n", $width, "",                       "Final P/F in the results pertains to the overlap check: 'P' if this accession has the same");
    printf("# %-*s %s\n", $width, "",                       "set of overlaps as the reference accession, and 'F' if not.");
  }    

  return; 
}

# Subroutine: debugPrintGapArray()
#
# Synopsis:   Print a gap array, for debugging purposes.
#
# Args:       $AR:            reference to the array to update
#             $rfpos:         reference position the gap occurs at or after
#             $is_delete:     '1' if we're updating a delete array, else we're
#                             updating an insert array (we do the update slightly
#                             differently for each type)
#
# Returns:    void
#
sub debugPrintGapArray {
  my $sub_name = "debugPrintGapArray";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($AR) = @_;

  if(defined $AR) { 
    my $nel = scalar(@{$AR});
    if($nel > 0) { 
      for(my $zz = 0; $zz < scalar(@{$AR}); $zz++) { 
        printf("$AR->[$zz]\n");
      }
    }
  }
  return;
}

# Subroutine: outputGapInfo()
#
# Synopsis:   Output gap information for all predicted CDS in a single sequence 
#             to a file handle.
#
# Args:       $perseq_FH:          output file handle to print per-sequence gap info to, undef to not print perseq info
#             $pergap_FH:          output file handle to print per-gap info to
#             $do_perseq_tbl:      '1' to output per sequence gaps as a table, '0' to print a list
#             $do_gap_all:         '1' to output per gap info for all gaps
#             $do_gap_not3:        '1' to output per gap info for gaps that are not a multiple of 3, not for all gaps
#             $do_gap_special:     '1' to output per gap info for special gaps that are possibly causative of a frameshift
#                                  Only 1 of $pergap_all, $pergap_not3, and $pergap_special can be '1'.
#             $mdl_AR:             reference to array of all model names
#             $seq_AR:             reference to array of all sequence names
#             $mdllen_HR:          ref to hash of model lengths
#             $hmm2cds_map_AR:     ref to array that maps each hmm to a CDS
#             $refdel_HHAR:        ref to 2D hash where value is an array, each element is
#                                  an rf position that is deleted in the alignment of the $seq_accn
#                                  pre-filled
#             $refins_HHAR:        ref to 2D hash where value is an array, each element is
#                                  a string <rfpos>:<ct> where <rfpos> is a rf pos in the alignment
#                                  after which $seq_accn has an insert and <ct> is the number of inserted
#                                  positions.
#
# Returns:    void
#
sub outputGapInfo {
  my $sub_name = "outputGapInfo";
  my $nargs_exp = 12;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($perseq_FH, $pergap_FH, $do_perseq_tbl, $do_gap_all, $do_gap_not3, $do_gap_special, $mdl_AR, $seq_AR, $mdllen_HR, $hmm2cds_map_AR, $refdel_HHAR, $refins_HHAR) = @_;

  if($do_gap_all) {
    if($do_gap_not3 || $do_gap_special) { die "ERROR in $sub_name, exactly one of $do_gap_all, $do_gap_not3, and $do_gap_special must be true"; }
  }
  elsif($do_gap_not3) { 
    if($do_gap_all || $do_gap_special) { die "ERROR in $sub_name, exactly one of $do_gap_all, $do_gap_not3, and $do_gap_special must be true"; }
  }
  elsif($do_gap_special) { 
    if($do_gap_all || $do_gap_not3) { die "ERROR in $sub_name, exactly one of $do_gap_all, $do_gap_not3, and $do_gap_special must be true"; }
  }
  else { 
    die "ERROR in $sub_name, exactly one of $do_gap_all, $do_gap_not3, and $do_gap_special must be true"; 
  }
 
  my @gapstr_AA = ();           # [0..$nseq-1][0..$ncds-1]: string describing all gaps for this sequence and this CDS
  my @w_gapstr_A = ();          # [0..$i..$ncds-1] max width of $gapstr_AA[$0..nseq-1][$i] for cds $i over all sequences

  my @tot_gap_length_AA = ();   # [0..$nseq-1][0..$ncds-1]: total number of gap positions for this sequence and this CDS
  my @w_tot_gap_length_A = ();  # [0..$i..$ncds-1] max width of $tot_gap_length_AA[$0..nseq-1][$i] for cds $i over all sequences

  my @net_gap_length_AA = ();   # [0..$nseq-1][0..$ncds-1]: net number of gap positions for this sequence and this CDS
  my @w_net_gap_length_A = ();  # [0..$i..$ncds-1] max width of $net_gap_length_AA[$0..nseq-1][$i] for cds $i over all sequences

  my @cds_gapstr_AH      = ();  # [0..$i..$ncds-1]: hash w/key: gapstring for a single position, value: number of times that gapstring occurs in any of @gapstr_AA for CDS $i
  
  my $ch_gapstr         = "string";
  my $ch_tot_gap_length = "tot";
  my $ch_net_gap_length = "net";

  $w_gapstr_A[0]         = length($ch_gapstr);
  $w_tot_gap_length_A[0] = length($ch_tot_gap_length);
  $w_net_gap_length_A[0] = length($ch_net_gap_length);
  %{$cds_gapstr_AH[0]}   = ();

  my $width_seq = length("#accession");
  my $cds_idx = 0;

  for(my $i = 0; $i < scalar(@{$seq_AR}); $i++) { 
    my $seq = $seq_AR->[$i];
    my $seq2print = $seq;
    $seq2print =~ s/\:.+$//;
    if(length($seq2print) > $width_seq) { 
      $width_seq = length($seq2print);
    }
    @{$gapstr_AA[$i]} = ();
    my $offset = 0;
    my $ins_idx = 0;
    my $del_idx = 0;
    my $next_ins_str   = undef;
    my $next_del_str   = undef;
    my $next_ins_rfpos = undef;
    my $next_del_rfpos = undef;
    my $next_ins_count = undef;
    my $next_del_count = undef;
    my $nhmm = scalar(@{$hmm2cds_map_AR});
    my $gapstr = "";
    my $substr;
    my $tot_gap_length = 0;
    my $net_gap_length = 0;
    $cds_idx = 0;

    for(my $h = 0; $h < $nhmm; $h++) { 
      my $mdl = $mdl_AR->[$h];
      my $ndel = (exists $refdel_HHAR->{$mdl}{$seq}) ? scalar(@{$refdel_HHAR->{$mdl}{$seq}}) : 0;
      my $nins = (exists $refins_HHAR->{$mdl}{$seq}) ? scalar(@{$refins_HHAR->{$mdl}{$seq}}) : 0;
      my $hmm_is_final = (($h == ($nhmm-1)) || ($hmm2cds_map_AR->[($h+1)] != $hmm2cds_map_AR->[$h])) ? 1 : 0; 
      my $hmm_is_first = (($h == 0)         || ($hmm2cds_map_AR->[($h-1)] != $hmm2cds_map_AR->[$h])) ? 1 : 0; 
      $ins_idx = 0;
      $del_idx = 0;
      $next_del_str = ($del_idx < $ndel) ? $refdel_HHAR->{$mdl}{$seq}[$del_idx] : undef;
      $next_ins_str = ($ins_idx < $nins) ? $refins_HHAR->{$mdl}{$seq}[$ins_idx] : undef;
      
      while(defined $next_ins_str || defined $next_del_str) { 
        # printf("next_ins_str: %s\n", (defined $next_ins_str) ? $next_ins_str : "undefined");
        # printf("next_del_str: %s\n", (defined $next_del_str) ? $next_del_str : "undefined");
        ($next_del_rfpos, $next_del_count) = (defined $next_del_str) ? split(":", $refdel_HHAR->{$mdl}{$seq}[$del_idx]) : (undef, undef);
        ($next_ins_rfpos, $next_ins_count) = (defined $next_ins_str) ? split(":", $refins_HHAR->{$mdl}{$seq}[$ins_idx]) : (undef, undef);
        if(defined $next_del_rfpos) { $next_del_rfpos += $offset; }
        if(defined $next_ins_rfpos) { $next_ins_rfpos += $offset; }
        
        if($gapstr ne "") { $gapstr .= ","; }
        if(defined $next_ins_str && defined $next_del_str) { 
          if($next_del_rfpos <= $next_ins_rfpos) { # delete comes first, print it
            $substr = "D" . $next_del_rfpos . ":" . $next_del_count;
            if($do_gap_all || (($next_del_count % 3) != 0)) { 
              $gapstr .= $substr;
              $tot_gap_length += $next_del_count;
              $net_gap_length -= $next_del_count;
              if(! $do_gap_special) { 
                $cds_gapstr_AH[$cds_idx]{$substr}++;
              }
            }
            $del_idx++;
            $next_del_str = ($del_idx < $ndel) ? $refdel_HHAR->{$mdl}{$seq}[$del_idx] : undef;
          }
          elsif($next_ins_rfpos < $next_del_rfpos) { # insert comes first, print it
            $substr = "I" . $next_ins_rfpos . ":" . $next_ins_count;
            if($do_gap_all || (($next_ins_count % 3) != 0)) { 
              $gapstr .= $substr;
              $tot_gap_length += $next_ins_count;
              $net_gap_length += $next_ins_count;
              if(! $do_gap_special) { 
                $cds_gapstr_AH[$cds_idx]{$substr}++;
              }
            }
            $ins_idx++;
            $next_ins_str = ($ins_idx < $nins) ? $refins_HHAR->{$mdl}{$seq}[$ins_idx] : undef;
          }
        }
        elsif(defined $next_del_str) { # $next_ins is undefined
          $substr = "D" . $next_del_rfpos . ":" . $next_del_count;
          if($do_gap_all || (($next_del_count % 3) != 0)) { 
            $gapstr .= $substr;
            $tot_gap_length += $next_del_count;
            $net_gap_length -= $next_del_count;
            if(! $do_gap_special) {
              $cds_gapstr_AH[$cds_idx]{$substr}++;
            }
          }
          $del_idx++;
          $next_del_str = ($del_idx < $ndel) ? $refdel_HHAR->{$mdl}{$seq}[$del_idx] : undef;
        }
        elsif(defined $next_ins_str) { # $next_del is undefined
          $substr = "I" . $next_ins_rfpos . ":" . $next_ins_count;
          if($do_gap_all | (($next_ins_count % 3) != 0))  { 
            $gapstr .= $substr;
            $tot_gap_length += $next_ins_count;
            $net_gap_length += $next_ins_count;
            if(! $do_gap_special) { 
              $cds_gapstr_AH[$cds_idx]{$substr}++;
            }
          }
          $ins_idx++;
          $next_ins_str = ($ins_idx < $nins) ? $refins_HHAR->{$mdl}{$seq}[$ins_idx] : undef;
        }        
      } # end of 'while(defined $next_ins_str || defined $next_del_str) { 
      
      $offset += $mdllen_HR->{$mdl};
      
      if($hmm_is_final) {
        if($gapstr eq "") { $gapstr = "-"; }

        # important to update $gapstr here, before we store it if $do_gap_special is true
        if($do_gap_special) { 
          # printf("calling findSpecialGap with $gapstr\n");
          $gapstr = findSpecialGap($gapstr);
          $net_gap_length = $tot_gap_length;
          if($gapstr ne "?" && $gapstr ne "-") { 
            my @el_A = split(",", $gapstr); 
            foreach my $el (@el_A) { 
              $cds_gapstr_AH[$cds_idx]{$el}++;
            }
          }
        }

        push(@{$gapstr_AA[$i]}, $gapstr);
        if(length($gapstr) > $w_gapstr_A[$cds_idx]) { $w_gapstr_A[$cds_idx] = length($gapstr); }

        push(@{$tot_gap_length_AA[$i]}, $tot_gap_length);
        if(length($tot_gap_length) > $w_tot_gap_length_A[$cds_idx]) { $w_tot_gap_length_A[$cds_idx] = length($tot_gap_length); }
        $tot_gap_length = 0;

        push(@{$net_gap_length_AA[$i]}, $net_gap_length);
        if(length($net_gap_length) > $w_net_gap_length_A[$cds_idx]) { $w_net_gap_length_A[$cds_idx] = length($net_gap_length); }
        $net_gap_length = 0;

        $gapstr = "";
        $offset = 0;
        $cds_idx++;

        if(scalar(@w_gapstr_A)         <= $cds_idx) { $w_gapstr_A[$cds_idx]         = length($ch_gapstr); }
        if(scalar(@w_tot_gap_length_A) <= $cds_idx) { $w_tot_gap_length_A[$cds_idx] = length($ch_tot_gap_length); }
        if(scalar(@w_net_gap_length_A) <= $cds_idx) { $w_net_gap_length_A[$cds_idx] = length($ch_net_gap_length); }
        if(scalar(@cds_gapstr_AH)      <= $cds_idx) { %{$cds_gapstr_AH[$cds_idx]}   = (); }
      }
    }
  }
  my $ncds = $cds_idx;

  #################################################################
  # Output table of all gap info
  # output line 1 of the column headers:
  if($do_gap_all) { 
    printf $perseq_FH ("# List of all gaps in alignment of each sequence:\n#\n");
  }
  elsif($do_gap_not3) { 
    printf $perseq_FH ("# List of all gaps with length that is not a multiple of 3 in alignment of each sequence:\n#\n");
  }
  else { # $do_gap_special 
    printf $perseq_FH ("# List of all gaps that may solely explain a CDS' length not being a multiple of 3, for each sequence:\n#\n");
  }

  if(! $do_gap_special) { 
    printf $perseq_FH  ("%-*s  ", $width_seq, "#");
    for(my $c = 0; $c < $ncds; $c++) { 
      my $w_cur = $w_tot_gap_length_A[$c] + 2 + $w_net_gap_length_A[$c] + 2 + $w_gapstr_A[$c];
      if($c > 0) { print $perseq_FH "  "; }
      printf $perseq_FH ("%-*s", $w_cur, "CDS#" . ($c+1));
    }
    print $perseq_FH "\n";
    
    # output line 2 (dashes under line 1 of column headers)
    printf $perseq_FH ("%-*s  ", $width_seq, "#");
    for(my $c = 0; $c < $ncds; $c++) { 
      my $w_cur = $w_tot_gap_length_A[$c] + 2 + $w_net_gap_length_A[$c] + 2 + $w_gapstr_A[$c];
      if($c > 0) { print $perseq_FH "  "; }
      printf $perseq_FH ("%-*s", $w_cur, monocharacterString($w_cur, "="));
    }    
    print $perseq_FH "\n";
  }

  # output line 3 of the column headers:
  printf $perseq_FH ("%-*s  ", $width_seq, "#accession");
  for(my $c = 0; $c < $ncds; $c++) { 
    if($c > 0) { print $perseq_FH "  "; }
    if(! $do_gap_special) { 
      printf $perseq_FH ("%-*s  ", $w_tot_gap_length_A[$c], $ch_tot_gap_length);
      printf $perseq_FH ("%-*s  ", $w_net_gap_length_A[$c], $ch_net_gap_length);
      printf $perseq_FH ("%-*s", $w_gapstr_A[$c], $ch_gapstr);
    }
    else { 
      printf $perseq_FH ("%-*s", $w_gapstr_A[$c], "CDS#" . ($c+1));
    }
  }
  print $perseq_FH "\n";

  # output line 4 (dashes under line 3 of column headers)
  printf $perseq_FH ("%-*s  ", $width_seq, "#" . monocharacterString($width_seq-1, "-"));
  for(my $c = 0; $c < $ncds; $c++) { 
    if($c > 0) { print $perseq_FH "  "; }
    if(!  $do_gap_special) { 
      printf $perseq_FH ("%-*s  ", $w_tot_gap_length_A[$c], monocharacterString($w_tot_gap_length_A[$c], "-"));
      printf $perseq_FH ("%-*s  ", $w_net_gap_length_A[$c], monocharacterString($w_net_gap_length_A[$c], "-"));
    }
    printf $perseq_FH ("%-*s", $w_gapstr_A[$c], monocharacterString($w_gapstr_A[$c], "-"));
  }
  print $perseq_FH "\n";

  # output actual data, for each sequence
  for(my $i = 0; $i < scalar(@{$seq_AR}); $i++) { 
    my $seq2print = $seq_AR->[$i];
    $seq2print =~ s/\:.+$//;
    printf $perseq_FH ("%-*s  ", $width_seq, $seq2print);
    for(my $c = 0; $c < $ncds; $c++) { 
      if($c > 0) { print $perseq_FH "  "; }
      if(! $do_gap_special) { 
        printf $perseq_FH ("%-*s  ", $w_tot_gap_length_A[$c], $tot_gap_length_AA[$i][$c]);
        printf $perseq_FH ("%-*s  ", $w_net_gap_length_A[$c], $net_gap_length_AA[$i][$c]);
      }
      printf $perseq_FH ("%-*s", $w_gapstr_A[$c], $gapstr_AA[$i][$c]);
    }
    print $perseq_FH ("\n");
  }

  # output explanatory text
  print $perseq_FH "#\n";
  print $perseq_FH ("# Explanation of the above table:\n");
  if($do_gap_all) { 
    print $perseq_FH ("# The table includes information on all gaps that exist between all pairwise alignments of\n");
    print $perseq_FH ("# the reference CDS and the predicted homologous CDS for each sequence.\n");
  }
  elsif($do_gap_not3) { 
    print $perseq_FH ("# The table includes information on all gaps of lengths that are not multiples of 3 that exist\n");
    print $perseq_FH ("# between all pairwise alignments of the reference CDS and the predicted homologous CDS for each sequence.\n");
  }
  else { 
    print $perseq_FH ("# The table includes information on some gaps that can solely explain a CDS not being a multiple of length 3.\n");
    print $perseq_FH ("# This is (probably) not an exhaustive list of all such gaps.\n");
    print $perseq_FH ("# Specifically it is only gaps X in a CDS Y, such that the following criteria are met:\n");
    print $perseq_FH ("#   - length of CDS Y is not a multiple of 3\n");
    print $perseq_FH ("#   - if you remove only X from list of all gaps, total CDS length of Y becomes a multiple of 3\n");
    print $perseq_FH ("#       with length difference of D with reference CDS\n");
    print $perseq_FH ("#    - there are no other gaps Z such that if you remove only Z then length of Y becomes a multiple\n");
    print $perseq_FH ("#       of 3 with length difference D2 from reference where D2 < D.\n");
  }
  print $perseq_FH ("#\n");
  if($do_gap_all || $do_gap_not3) { 
    printf $perseq_FH ("# There are 3 columns under each header \"CDS#<n> (%s)\" named \"tot\", \"net\",\n", ($do_gap_all) ? "all gaps" : "gaps %3 != 0");
    print $perseq_FH ("# and \"string\".\n");
    print $perseq_FH ("# The \"tot\" columns list the total number of gap positions in either sequence in the pairwise alignment.\n");
    print $perseq_FH ("# The \"net\" columns list the net number of the listed gaps in the pairwise alignment; this is the number\n");
    print $perseq_FH ("#   of gaps in the reference sequence minus the number of gaps in the current sequence (inserts minus deletes)\n");
    print $perseq_FH ("# The \"string\" columns include a list of <n> tokens, each of which describes a gap of length >= 1 nucleotide.\n");
  }
  print $perseq_FH ("#\n");
  print $perseq_FH ("# Tokens are in the form: <char><position><length>\n");
  print $perseq_FH ("#   <char>     is 'I' for an insertion relative to the reference CDS (gap in reference sequence)\n");
  print $perseq_FH ("#              or 'D' for a  deletion  relative to the reference CDS (gap in current sequence)\n");
  print $perseq_FH ("#   <position> is the nucleotide position of the gap in reference coordinates.\n");
  print $perseq_FH ("#              For insertions this is the reference position after which the insertion occurs.\n");
  print $perseq_FH ("#              For deletions  this is the first reference position for this deletion.\n");
  print $perseq_FH ("#   <length>   length of the gap in nucleotides.\n");
  print $perseq_FH ("#              For insertions this is the number of nucleotides inserted relative to the reference\n");
  print $perseq_FH ("#              For deletions  this is the number of reference positions deleted.\n");
  print $perseq_FH ("#\n");
  if($do_gap_special) { 
    print $perseq_FH ("#\n");
    print $perseq_FH ("# \"-\" tokens indicate the CDS is a multiple of length 3\n");
    print $perseq_FH ("# \"?\" tokens indicate the CDS is not a multiple of length 3, but that no gaps that satisfy our criteria exist.\n");
    print $perseq_FH ("#\n");
  }

  # Now the per gap information
  if($do_gap_all) { 
    printf $pergap_FH ("# Counts of all gaps:\n#\n");
  }
  elsif($do_gap_not3) { 
    printf $pergap_FH ("# Counts of gaps that have length that is not a multiple of 3:\n#\n");
  }
  else { # $do_gap_special == 1
    printf $pergap_FH ("# Counts of gaps that are special (responsible for net gap length that is not a multiple of 3):\n#\n");
    print $perseq_FH ("# Specifically, these are counts of gaps X in a CDS Y, such that the following criteria are met:\n");
    print $perseq_FH ("#   - length of CDS Y is not a multiple of 3\n");
    print $perseq_FH ("#   - if you remove only X from list of all gaps, total CDS length of Y becomes a multiple of 3\n");
    print $perseq_FH ("#       with length difference of D with reference CDS\n");
    print $perseq_FH ("#    - there are no other gaps Z such that if you remove only Z then length of Y becomes a multiple\n");
    print $perseq_FH ("#       of 3 with length difference D2 from reference where D2 < D.\n");
    print $perseq_FH ("#\n");
  }
  my $nprinted = 0;
  for(my $c = 0; $c < $ncds; $c++) { 
    if((scalar(keys %{$cds_gapstr_AH[$c]})) > 0) { 
      foreach my $key (sort keys %{$cds_gapstr_AH[$c]}) { 
        printf $pergap_FH ("CDS#" . ($c+1) . " " . $key . " " . $cds_gapstr_AH[$c]{$key} . "\n");
        $nprinted++;
      }
      printf $pergap_FH ("#\n");
    }
  }
  if($nprinted == 0) { 
    printf $pergap_FH ("# NONE\n");
  }

  return;
}
###########################################################

