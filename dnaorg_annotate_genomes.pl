#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began]
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
$usage .= "  -oseq <s>  : identify origin seq <s> in genomes, put \"|\" at site of origin, e.g. \"TAATATT\\|AC\"\n";
$usage .= "               <s> must be a string consisting of only A,C,G,T and | characters. No regular expressions allowed.\n"; 
$usage .= "               Note that \"|\" must be escaped, i.e. \"\\|\"; the \"|\" can be any char, incl. first or last\n";
$usage .= "               This option is relevant only for circular genomes, e.g. Maize streak virus\n";
$usage .= "  -strict    : require matching annotations to match CDS/exon index\n";
$usage .= "  -nodup     : do not duplicate genome seqs to identify features that span stop..start (for circular genomes)\n";
$usage .= "  -notexon   : do not use exon-specific models\n";
$usage .= "  -onlybuild : exit after building reference models\n";
$usage .= "  -skipbuild : skip the build and calibration step because you already did an -onlybuild run\n";
$usage .= "  -model <s> : use model file <s>, instead of building one\n";
$usage .= "\n OPTIONS CONTROLLING OUTPUT TABLE:\n";
$usage .= "  -c        : concise output mode (enables -nomdlb -noexist -nobrack and -nostop)\n";
$usage .= "  -seqrow   : force sequence-as-rows    output mode (default for <= 5 models)\n";
$usage .= "  -seqcol   : force sequence-as-columns output mode (default for >  5 models)\n";
$usage .= "  -nseqcol  : with -seqcol, number of sequences per page of output (default: 10)\n";
$usage .= "  -nomdlb   : do not add model boundary annotation to output\n";
$usage .= "  -noexist  : do not include information on existing GenBank annotation\n";
$usage .= "  -nobrack  : do not include brackets around predicted annotations that do not match existing\n";
$usage .= "  -nostop   : do not output stop codon for each predicted CDS/mat_peptide\n";
$usage .= "  -nofid    : do not output fractional identity relative to reference for each CDS/exon\n";
$usage .= "  -noss3    : do not output results of start codon, stop codon and multiple of 3 tests\n";
$usage .= "  -noolap   : do not output overlap information about which predicted features overlap\n";
$usage .= "  -noexp    : do not output explanation of column headings\n";
$usage .= "  -fullolap : do output full overlap   string (default: do not)\n";
$usage .= "  -fulladj  : do output full adjacency string (default: do not)\n";
$usage .= "\n OPTIONS FOR SELECTING HOMOLOGY SEARCH ALGORITHM:\n";
$usage .= "  -hmmer     : use HMMER for predicting annotations, default: use Infernal\n";
$usage .= "\n OPTIONS SPECIFIC TO HMMER3:\n";
$usage .= "  -hmmenv  : use HMM envelope boundaries for predicted annotations, default: use window boundaries\n";
$usage .= "\n OPTIONS SPECIFIC TO INFERNAL:\n";
$usage .= "  -iglocal   : use the -g option with cmsearch for glocal searches\n";
$usage .= "  -cslow     : use default cmcalibrate parameters, not parameters optimized for speed\n";
$usage .= "  -cfarm     : submit calibration jobs for each CM to compute farm and exit (requires --onlybuild)\n";
$usage .= "  -sfarm <n> : split genome file into <n> pieces, submit <n> cmscan jobs to farm and wait 3 minutes\n";
$usage .= "               (changeable with -swait) before concatenating all the output files and continuing\n";
$usage .= "  -swait <n>  : with -sfarm, set number of minutes to wait for cmscan jobs to finish to <n> [df: 3]\n";
$usage .= "\n OPTIONS USEFUL FOR DEVELOPMENT/DEBUGGING:\n";
$usage .= "  -skipfetch : don't fetch whole genome sequences, we already have them from a previous run\n";
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
my $cmfetch         = $inf_exec_dir . "cmfetch";
my $cmpress         = $inf_exec_dir . "cmpress";
my $cmscan          = $inf_exec_dir . "cmscan";
my $cmalign         = $inf_exec_dir . "cmalign";

foreach my $x ($hmmbuild, $hmmpress, $nhmmscan, $cmbuild, $cmcalibrate, $cmfetch, $cmpress, $cmscan, $esl_translate, $esl_ssplit) { 
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
my $do_skipbuild = 0; # set to '1' if -skipbuild  enabled, skip the build step
my $in_model_db  = undef; # defined if -model <s> enabled, use <s> as the model file instead of building one
# options for controlling output table
my $do_concise   = 0; # set to '1' if -c       enabled, invoke concise output mode, set's all $do_no* variables below to '1'
my $do_seqrow    = 0; # set to '1' if -seqrow  enabled, force sequence-as-rows output model
my $do_seqcol    = 0; # set to '1' if -seqcol  enabled, force sequence-as-cols output model
my $df_nseqcol   = 10;          # default number of sequences to a page in sequence-as-cols mode
my $nseqcol      = $df_nseqcol; # changed to <n> if -nseqcol <n> enabled, set number of sequences per page in sequence-as-cols mode to $nseqcol
my $do_nomdlb    = 0; # set to '1' if -nomdlb  or -c enabled, do not print model boundary info for annotations, else do
my $do_noexist   = 0; # set to '1' if -noexist or -c enabled, do not output information on existing GenBank annotations
my $do_nobrack   = 0; # set to '1' if -nobrack or -c enabled, do not output brackets around predicted annotations that do not match any existing GenBank annotation
my $do_nostop    = 0; # set to '1' if -nostop  or -c enabled, do not output stop codon for predicted annotations
my $do_nofid     = 0; # set to '1' if -nofid   or -c enabled, do not output fractional identities relative to the reference
my $do_noss3     = 0; # set to '1' if -noss3   or -c enabled, do not output SS3 columns: 'S'tart codon check, 'S'top codon check and multiple of '3' check
my $do_noolap    = 0; # set to '1' if -noolap  or -c enabled, do not output information on overlapping features
my $do_noexp     = 0; # set to '1' if -noexp   or -c enabled, do not output explanatory information about column headings
my $do_fullolap  = 0; # set to '1' if -fullolap enabled, do output full overlap   string for all features (default: do not)
my $do_fulladj   = 0; # set to '1' if -fullolap enabled, do output full adjacency string for all features (default: do not)
# options for controlling homology search method
my $do_hmmer     = 0; # set to '1' if -hmmer      enabled, use HMMER3's nhmmscan, not Infernal 1.1
# options specific to HMMER3
my $do_hmmenv    = 0; # set to '1' if -hmmenv     enabled, use HMM envelope boundaries as predicted annotations, else use window boundaries
# options specific to Infernal
my $do_iglocal    = 0; # set to '1' if -iglocal enabled, use -g with cmsearch
my $do_cslow      = 0; # set to '1' if -cslow   enabled, use default, slow, cmcalibrate parameters instead of speed optimized ones
my $do_cfarm      = 0; # set to '1' if -cfarm   enabled, submit cmcalibrate job to farm
my $do_sfarm      = 0; # set to '1' if -sfarm   enabled, submit cmscan jobs to farm
my $sfarm_njobs   = undef; # set to a value if -sfarm is used
my $df_sfarm_wait = 3; # default value for number of minutes to wait for cmscan jobs to finish
my $sfarm_wait    = $df_sfarm_wait; # changeable to <n> with -swait <n>
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
            "skipbuild" => \$do_skipbuild,
            "model=s"   => \$in_model_db,
            "c"         => \$do_concise,
            "seqrow"    => \$do_seqrow,
            "seqcol"    => \$do_seqcol,
            "nseqcol=s" => \$nseqcol,
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
            "cfarm"     => \$do_cfarm,
            "sfarm=s"   => \$sfarm_njobs,
            "swait=s"   => \$sfarm_wait,
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
  $opts_used_long  .= "# option:  do not correct after identifying internal start/stops in predicted features [-nocorrect]\n";
}
if($do_matpept) { 
  $opts_used_short .= "-matpept ";
  $opts_used_long  .= "# option:  using mat_peptide info instead of CDS info [-matpept]\n";
}
if(defined $origin_seq) { 
  $opts_used_short .= "-oseq $origin_seq ";
  $opts_used_long  .= "# option:  searching for origin sequence of $origin_seq [-oseq]\n";
}
if($do_strict) { 
  $opts_used_short .= "-strict ";
  $opts_used_long  .= "# option:  demand matching annotations are same indexed CDS/exon or mat_peptide [-strict]\n";
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
if($do_skipbuild) { 
  $opts_used_short .= "-skipbuild ";
  $opts_used_long  .= "# option:  skipping the build step [-skipbuild]\n";
}
if(defined $in_model_db) { 
  $opts_used_short .= "-model $in_model_db ";
  $opts_used_long  .= "# option:  use model in $in_model_db instead of building one here [-model]\n";
}
if($do_concise) { 
  $opts_used_short .= "-c ";
  $opts_used_long  .= "# option:  concise output mode [-c]\n";
}
if($do_seqrow) { 
  $opts_used_short .= "-seqrow ";
  $opts_used_long  .= "# option:  forcing seqs-as-rows output mode [-seqrow]\n";
}
if($do_seqcol) { 
  $opts_used_short .= "-seqcol ";
  $opts_used_long  .= "# option:  forcing seqs-as-columns output mode [-seqcol]\n";
}
if($nseqcol ne $df_nseqcol) { 
  $opts_used_short .= "-nseqcol $nseqcol";
  $opts_used_long  .= "# option:  setting number of sequences per output page to $nseqcol [-nseqcol]\n";
}
if($do_nomdlb) { 
  $opts_used_short .= "-nomdlb ";
  $opts_used_long  .= "# option:  do not output model boundaries of predicted annotations [-nomdlb]\n";
}
if($do_noexist) { 
  $opts_used_short .= "-noexist";
  $opts_used_long  .= "# option:  not outputting info on existing GenBank annotations [-noexist]\n";
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
if($do_fullolap) { 
  $opts_used_short .= "-fullolap";
  $opts_used_long  .= "# option:  do output full overlap string for all genes [-fullolap]\n";
}
if($do_fulladj) { 
  $opts_used_short .= "-fulladj";
  $opts_used_long  .= "# option:  do output full adjacency string for all genes [-fulladj]\n";
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
if($do_cfarm) { 
  $opts_used_short .= "-cfarm ";
  $opts_used_long  .= "# option:  submitting calibration jobs to farm [-cfarm]\n";
}
if(defined $sfarm_njobs) { 
  $do_sfarm = 1;
  $opts_used_short .= "-sfarm $sfarm_njobs ";
  $opts_used_long  .= "# option:  submit cmscan jobs to farm and wait for them to finish [-sfarm]\n";
}
if($sfarm_wait != $df_sfarm_wait) { 
  $opts_used_short .= "-swait $sfarm_wait ";
  $opts_used_long  .= "# option:  submit cmscan jobs to farm and wait for them to finish [-sfarm]\n";
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
if($do_cfarm && $do_hmmer) { 
  die "ERROR -cfarm is incompatible with -hmmer"; 
}
if($do_cslow && $do_hmmer) { 
  die "ERROR -cslow is incompatible with -hmmer"; 
}
if($do_skipscan && ($do_cfarm || $do_cslow || $do_onlybuild)) { 
  die "ERROR -skipscan is incompatible with -cfarm, -cslow, and -onlybuild";
}
if($do_sfarm && ($do_skipaln || $do_skipscan)) { 
  die "ERROR -skipscan and -skipaln are incompatible with -sfarm";
}
if($do_sfarm && $do_hmmer) { 
  die "ERROR -sfarm and -hmmer are incompatible";
}
if($do_onlybuild && $do_skipbuild) { 
  die "ERROR -onlybuild and -skipbuild are incompatible";
}
if($do_onlybuild && (defined $in_model_db)) { 
  die "ERROR -onlybuild and -model are incompatible";
}
if($do_notexon && $do_matpept) { 
  die "ERROR -noexon is incompatible with -matpept";
}

# check that options that must occur in combination, do
if($do_cfarm && (! $do_onlybuild)) { 
  die "ERROR -cfarm must be used in combination with -onlybuild"; 
}
if(($sfarm_wait != $df_sfarm_wait) && (! $do_sfarm)) { 
  die "ERROR -swait must be used in combination with -sfarm"; 
}
if((defined $in_model_db) && (! $do_skipbuild)) {
  die "ERROR -skipbuild must be used in combination with -model"; 
}
# check that options that must occur in combination, do
if($nseqcol ne $df_nseqcol && (! $do_seqcol)) { 
  die "ERROR -nseqcol must be used in combination with -seqcol";
}

# check that input files related to options actually exist
if(defined $in_model_db) { 
  if(! -s $in_model_db) { die "ERROR: $in_model_db file does not exist"; }
}
# verify origin sequence if necessary
my $origin_offset = undef;
if(defined $origin_seq) { 
  $origin_seq =~ tr/a-z/A-Z/; # capitalize origin seq
  $origin_offset = validateOriginSeq($origin_seq);
  # printf("origin_offset: $origin_offset\n");
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
# NOTE: the table and length files get created external to this program and have prescribed names
#####################
my @accn_A        = (); # array of accessions
my %accn_exists_H = (); # hash of accessions, key is accession, value is always '1', used only to check for duplicates
open(IN, $listfile) || die "ERROR unable to open $listfile for reading"; 
my $waccn = 0; # max length of all accessions
while(my $accn = <IN>) { 
  if($accn =~ m/\w/) { 
    chomp $accn;
    stripVersion(\$accn); # remove version
    if(exists $accn_exists_H{$accn}) { 
      die "ERROR, the accession $accn appears twice in the input list file $listfile"; 
    }
    push(@accn_A, $accn);
    $accn_exists_H{$accn} = 1;
    if(length($accn) > $waccn) { $waccn = length($accn); }
  }
}
close(IN); 
%accn_exists_H = (); # free the memory, we don't need this anymore

##################################
# parse the table and length files
##################################
my %mft_tbl_HHA = ();   # Main feature (CDS or mat_peptide) data 
                        # from .cds.tbl or .mat_peptide.tbl file
                        # hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %totlen_H = (); # key: accession, value length read from length file

parseLength($length_file, \%totlen_H);

if($do_matpept) { 
  parseTable($matpept_tbl_file, \%mft_tbl_HHA);
}
else {
  parseTable($cds_tbl_file, \%mft_tbl_HHA);
}

# if we're in matpept mode, determine which peptides should be adjacent to each other at 
# the nucleotide level

#######################
# variable declarations
#######################
my $strand_str;                # +/- string for all main features (CDS or mat_peptide) for an accession: e.g. '+-+': 1st and 3rd CDS are + strand, 2nd is -

# reference information on reference accession, first accession read in ntlist file
my $ref_accn          = undef; # changed to <s> with -ref <s>
my $ref_label_str     = undef; # label string for reference accn
my $ref_nmft          = 0;     # number of main features (mft) (CDS or mat_peptide) in reference
my $ref_strand_str    = "";    # strand string for reference 
my @ref_mft_len_A     = ();    # [0..$i..$ref_nmft-1]: length of each reference main feature (CDS or mat_peptide)
#my @ref_mft_len_tol_A = ();   # [0..$i..$ref_nmft-1]: length tolerance, any gene that is within this fraction of the lenght of the ref gene is a match
my @ref_mft_coords_A  = ();    # [0..$i..$ref_nmft-1]: main feature coords (CDS or mat_peptide) for reference
my @ref_mft_product_A = ();    # product qualifier data for reference CDS or mat_peptide

my $nmft = 0;              # number of main feature (mft, CDS or mat_peptide)
my $npos = 0;              # number of main features on positive strand
my $nneg = 0;              # number of main features on negative strand
my $nunc = 0;              # number of main features on uncertain strand
my $nbth = 0;              # number of main features on both strands
my @mft_len_A = ();        # [0..$i..$nmft-1] length of main feature $i
my @mft_coords_A = ();     # [0..$i..$nmft-1] coords of main feature $i
my @mft_product_A = ();    # [0..$i..$nmft-1] product annotation for main feature $i

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
if($do_sfarm && ($sfarm_njobs > $naccn)) { 
  die "ERROR with -sfarm <n>, <n> must be <= number of genomes ($naccn), but you used $sfarm_njobs"; 
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

# make sure we have our reference accession in the sequence file
my $test_len = $sqfile->fetch_seq_length_given_name($ref_seq_accn);
if((! defined $test_len) || ($test_len == -1)) { 
  my $errmsg;
  if($do_nodup) { 
    $errmsg = "ERROR, unable to fetch reference sequence $ref_seq_accn from $gnm_fasta_file.\nDid you use -nodup when you did the first run of the script with -onlybuild?\nIf not, do not use it here. If so, and you want to use it, redo this run without -skipfetch.\n"; 
  }
  else { 
    $errmsg = "ERROR, unable to fetch reference sequence $ref_seq_accn from $gnm_fasta_file.\nDid you use -nodup when you did the first run of the script with -onlybuild?\nIf so, use it here too.\n";
  }
  die $errmsg;
}

##################################################
# If we're looking for an origin sequence, do that
##################################################
my %origin_coords_HA = ();
if(defined $origin_seq) {
  findSeqInFile($sqfile, $origin_seq, $do_nodup, \%origin_coords_HA);
}

########################################################
# Gather information and sequence data on the reference.
# Use each reference CDS exon or mat_peptide (if $do_matpept)
# as a homology search model against all the genomes.
#######################################################
$ref_accn = $accn_A[0];
if(! exists ($mft_tbl_HHA{$ref_accn})) { die sprintf("ERROR no %s information stored for reference accession", ($do_matpept) ? "mat_peptide" : "CDS"); }
(undef, undef, undef, undef, undef, $ref_strand_str) = getStrandStats(\%mft_tbl_HHA, $ref_accn);
getLengthStatsAndCoordStrings(\%mft_tbl_HHA, $ref_accn, \@ref_mft_len_A, \@ref_mft_coords_A);
getQualifierValues(\%mft_tbl_HHA, $ref_accn, "product", \@ref_mft_product_A);
$ref_nmft = scalar(@ref_mft_len_A);

## if we're in matpept mode, determine which peptides should be adjacent to each other
##my %adj2mft_5p_H = (); # key is CDS/peptide index in ref_mft* data structures, value is index of CDS/peptide that is immediately adjacent 5', if any
##my %adj2mft_3p_H = (); # key is CDS/peptide index in ref_mft* data structures, value is index of CDS/peptide that is immediately adjacent 3', if any
##matpeptFindAdjacentPeptides(\@ref_mft_coords_A, \%adj2mft_5p_H, \%adj2mft_3p_H);

my $all_stk_file = $out_root . ".ref.all.stk";

($seconds, $microseconds) = gettimeofday();
my $start_time = ($seconds + ($microseconds / 1000000.));
printf("%-65s ... ", sprintf("# Fetching reference %s sequences", ($do_matpept) ? "mat_peptide" : "CDS"));
my $cur_out_root;
my $cur_name_root;
my $fetch_input;
my $fetch_output;
my $nmdl = 0;               # number of HMMs (and alignments used to build those HMMs)
my @mdl2mft_map_A = ();     # [0..h..$nmdl-1]: $i: model ($m+1) maps to reference feature ($i+1)
my @mft2exon_map_A = ();    # [0..h..$nmdl-1]: $e: model ($m+1) maps to exon ($e+1) of reference cds $mdl2mft_map_A[$h]+1
my @mdl_is_first_A = ();    # [0..h..$nmdl-1]: '1' if model ($m+1) is the first one for cds $mdl2mft_map_A[$h], else 0
my @mdl_is_final_A = ();    # [0..h..$nmdl-1]: '1' if model ($m+1) is the final one for cds $mdl2mft_map_A[$h], else 0
my @mdl_is_primary_A = ();  # [0..h..$nmdl-1]: '1' if model ($m+1) is for a primary matpept, else 0, only created/used if $do_matpept
my @mft2first_mdl_A = ();   # [0..$c..ncds-1]: $h, first exon of feature $f+1 is modeled by model $m+1
my @mft2final_mdl_A = ();   # [0..$c..ncds-1]: $h, final exon of feature $f+1 is modeled by model $m+1
my @mdl_A = ();             # [0..$nmdl-1]: array of model names, also name of stockholm alignments used to build those models
my @mft_out_short_A   = (); # [0..$ref_nmft-1]: array of abbreviated model main feature names to print
my @mft_out_product_A = (); # [0..$ref_nmft-1]: array of 'product' qualifiers for a main feature (protein names)
my %mdllen_H          = (); # key: model name from @mdl_A, value is model length
my @ref_nexons_A      = (); # [0..$c..$ref_nmft-1]: number of exons in CDS or mat_peptide $c+1
my $ref_tot_nexons    = 0;  # total number of exons in all CDS or mat_peptides
my @indi_ref_name_A   = (); # [0..$nmdl-1]: name of individual stockholm alignments and models
my @indi_cksum_stk_A  = (); # [0..$nmdl-1]: checksum's of each named individual stockholm alignment

my $cds_or_mp = ($do_matpept) ? "mp" : "cds";
# for each reference main feature (CDS or mat_peptides), fetch each exon (or the full CDS if -notexon enabled)
for(my $i = 0; $i < $ref_nmft; $i++) { 
  # determine start and stop positions of all exons
  my @starts_A = ();
  my @stops_A  = ();
  my $nexons   = 0;
  startStopsFromCoords($ref_mft_coords_A[$i], \@starts_A, \@stops_A, \$nexons);
  if($do_matpept && $nexons > 1) { die "ERROR multi-exon CDS in matpept mode"; }
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
      $cur_out_root  = $out_root . ".ref." . $cds_or_mp . "." . ($i+1) . ".exon." . ($e+1);
      $cur_name_root = $dir_tail . ".ref." . $cds_or_mp . "." . ($i+1) . ".exon." . ($e+1);
    }
    else { 
      $cur_out_root  = $out_root . ".ref." . $cds_or_mp . "." . ($i+1);
      $cur_name_root = $dir_tail . ".ref." . $cds_or_mp . "." . ($i+1);
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
    my ($mdllen, $cksum) = annotateStockholmAlignment($cur_name_root, $do_blank_ss, $cur_stkfile, $cur_named_stkfile);
    push(@indi_ref_name_A, $cur_name_root);
    push(@indi_cksum_stk_A, $cksum);

    # store information on this model's name for output purposes
    if($e == ($nexons-1)) { 
      my $short = ($do_matpept) ? sprintf("MP #%d", ($i+1)) : sprintf("CDS #%d", ($i+1));
      if($act_nexons > 1) { $short .= " [$act_nexons exons; $strand]"; }
      else                { $short .= " [single exon; $strand]"; }
      push(@mft_out_short_A,   $short);
      push(@mft_out_product_A, $ref_mft_product_A[$i]);
    }
    push(@mdl_A, $cur_name_root);

    $mdllen_H{$cur_name_root} = $mdllen;

    # now append the named alignment to the growing stockholm alignment database $all-stk_file
    $cmd = "cat $cur_named_stkfile";
    if($nmdl == 0) { $cmd .= " >  $all_stk_file"; }
    else           { $cmd .= " >> $all_stk_file"; }
    runCommand($cmd, 0);
    push(@mdl2mft_map_A,  $i);
    push(@mft2exon_map_A, $e);
    push(@mdl_is_first_A, ($e == 0)           ? 1 : 0);
    push(@mdl_is_final_A, ($e == ($nexons-1)) ? 1 : 0);
    if($e == 0)           { $mft2first_mdl_A[$i] = $nmdl; }
    if($e == ($nexons-1)) { $mft2final_mdl_A[$i] = $nmdl; }
    $nmdl++;
  }
}
($seconds, $microseconds) = gettimeofday();
my $stop_time = ($seconds + ($microseconds / 1000000.));
printf("done. [%.1f seconds]\n", ($stop_time - $start_time));

# Now that we know the number of models, 
# determine the output mode if it was not specified on the cmdline:
my $df_max_nseqrow = 5;
if((! $do_seqrow) && (! $do_seqcol)) { 
  if($nmdl <= $df_max_nseqrow) { 
    $do_seqrow = 1; 
  }
  else { 
    $do_seqcol = 1;
  }
}

######################
# HOMOLOGY SEARCH STEP
######################
my $model_db; # model database file, either HMMs or CMs

# first, create the model database, unless we're skipping that step
if($do_skipbuild) { 
  if(defined $in_model_db) { 
    $model_db = $in_model_db;
    if(! -s $model_db) { die "ERROR: -model used but $in_model_db file does not exist"; }
  }
  else { 
    # -model not used, but we're trying to skip the build step
    # make sure we already have the CM DB file already (this will likely be true
    # if we've run this script with -skipbuild for this dataset >= 1 times already)
    # or if not (this will likely be true if we've run this script 0 times with 
    # -skipbuild already), then create that file by concatenating the individually
    # calibrated model files we created on the previous run of this script with
    # -onlybuild -sfarm.
    $model_db = $out_root . ".ref.cm";
    if(! -s $model_db) { 
      # case 2
      # concatenate the files into a CM DB and run cmpress on it
      for(my $i = 0; $i < $nmdl; $i++) { 
        my $indi_model = $out_root . ".ref.$i.cm";
        if(! -s ($indi_model)) { 
          die "ERROR: -skipbuild used but model db file $model_db does not exist, nor does individual model $indi_model;\nMake sure you've already run this script with -onlybuild -cfarm on this dataset -- it doesn't seem like you have.";
        }
      }
      my $cat_cmd = "cat $out_root.ref.0.cm > $model_db";
      runCommand($cat_cmd, 0);
      for(my $i = 1; $i < $nmdl; $i++) { 
        $cat_cmd = "cat $out_root.ref.$i.cm >> $model_db";
        runCommand($cat_cmd, 0);
      }
      # remove the binary files, possibly from an earlier cmbuild/cmpress:
      for my $suffix ("i1m", "i1i", "i1f", "i1p") { 
        my $file = $out_root . ".cm." . $suffix;
        if(-e $file) { unlink $file; }
      }
      my $cmpress_cmd = "$cmpress $model_db > $model_db.cmpress";
      # we do have all the models if we get here, so concatenate them and press it
      printf("%-65s ... ", "# Running cmpress");
      my $secs_elapsed = runCommand($cmpress_cmd, 0);
      #printf("\n$cmpress_cmd\n");
      printf("done. [%.1f seconds]\n", $secs_elapsed);
    }
  }
  # now we know we have the model file, extract the checksum values
  # and check them against the individual stockholm alignment checksum values

  ($seconds, $microseconds) = gettimeofday();
  my $cur_start_time = ($seconds + ($microseconds / 1000000.));
  printf("%-65s ... ", "# Validating CM file was built from correct alignment files");
  my $cksum_file = $model_db . ".cksum";
  $cmd = "grep ^CKSUM $model_db | awk '{ print \$2 '} > $cksum_file";
  runCommand($cmd, 0);
  open(CKSUM, $cksum_file) || die "ERROR unable to open $cksum_file for reading";
  my $i = 0;
  while(my $cksum = <CKSUM>) { 
    chomp $cksum;
    if($cksum != $indi_cksum_stk_A[$i]) { 
      die sprintf("ERROR checksum mismatch for CM %d (CM: %d != alignment: %d)", $i+1, $cksum, $indi_cksum_stk_A[$i]); 
    }
    $i++;
  }
  close(CKSUM);
  ($seconds, $microseconds) = gettimeofday();
  my $cur_stop_time = ($seconds + ($microseconds / 1000000.));
  printf("done. [%.1f seconds]\n", ($cur_stop_time - $cur_start_time));
}
else { 
  if(! $do_hmmer) { # use Infernal (default)
    createCmDb($cmbuild, $cmcalibrate, $cmpress, $cmfetch, $nmdl, $do_cslow, $do_cfarm, $all_stk_file, $out_root . ".ref", \@indi_ref_name_A);
    if($do_onlybuild) { 
      printf("#\n# Model calibration %s. Exiting.\n", ($do_cfarm) ? "job(s) submitted" : "complete");
      exit 0;
    }
    $model_db = $out_root . ".ref.cm";
  }
  else { # use HMMER3's nhmmscan
    createHmmDb($hmmbuild, $hmmpress, $nmdl, $all_stk_file, $out_root . ".ref");
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

my %pred_fafile_H = (); # hash of names of fasta files for predicted exon sequences, keys: model name from @mdl_A, value: name of file

if($do_hmmer) { 
  runNhmmscan($nhmmscan, ($do_skipscan || $do_skipaln), $model_db, $gnm_fasta_file, $tblout, $stdout);
  parseNhmmscanTblout($tblout, $do_hmmenv, \%totlen_H, \%p_start_HH, \%p_stop_HH, \%p_strand_HH, \%p_score_HH, \%p_hangover_HH);
  fetchHits($sqfile, $do_skipaln, "predicted", \@mdl_A, \@seq_accn_A, \%totlen_H, \%p_start_HH, \%p_stop_HH, \%p_strand_HH, $out_root, \%pred_fafile_H);
}
else { 
  if(! $do_sfarm) { 
    # default method, run cmscan on full genome file
    runCmscan($cmscan, $do_iglocal, ($do_skipscan || $do_skipaln), 0, $model_db, $gnm_fasta_file, $tblout, $stdout);
  }
  else { # we need to split up the genome fasta file and submit a different cmscan job for each fasta file
    # note we may redefined $sfarm_njobs here, splitFastaFile will return actual number of fasta files created, 
    # which can differ from the requested amount (which is $snfarm_njobs) that we pass in
    my $nfasta_created = splitFastaFile($esl_ssplit, $gnm_fasta_file, $sfarm_njobs);
    $sfarm_njobs = $nfasta_created;
    # now submit a job for each
    printf("%-65s ... ", sprintf("# Submitting %d cmscan jobs", $sfarm_njobs));
    for(my $z = 1; $z <= $sfarm_njobs; $z++) { 
      my $cur_gnm_fasta_file = $gnm_fasta_file . "." . $z;
      my $cur_tblout = $tblout . "." . $z;
      my $cur_stdout = $stdout . "." . $z;
      runCmscan($cmscan, $do_iglocal, 0, 1, $model_db, $cur_gnm_fasta_file, $cur_tblout, $cur_stdout);
    }
    printf("done. [waiting max of $sfarm_wait minutes for them to finish]\n");

    # now check that they're all finished every 15 seconds until we hit our max wait time
    my $all_finished = 0;
    my $nfinished = 0;
    for(my $zz = 0; $zz < ($sfarm_wait*4); $zz++) { 
      # check to see if jobs are finished, every 30 seconds
      sleep(15);
      $nfinished = 0;
      for(my $z = 1; $z <= $sfarm_njobs; $z++) { 
        my $cur_stdout = $stdout . "." . $z;
        if(-s $cur_stdout) { 
          my $final_line = `tail -n 1 $cur_stdout`;
          chomp $final_line;
          if($final_line eq "[ok]") { 
            $nfinished++;
          }
        }
      }
      if($nfinished == $sfarm_njobs) { 
        # we're done break out of it
        $all_finished = 1;
        $zz = ($sfarm_wait * 4);
      }
    }
    if($all_finished) { 
      # concatenate all outputs into one main one
      for(my $z = 1; $z <= $sfarm_njobs; $z++) { 
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
      die "ERROR only $nfinished of the $sfarm_njobs are finished after $sfarm_wait minutes. Increase limit with -swait";
    }
  } # end of 'else' entered if we're submitting jobs to the farm

  parseCmscanTblout($tblout, \%totlen_H, \%mdllen_H, \%p_start_HH, \%p_stop_HH, \%p_strand_HH, \%p_score_HH, \%p_hangover_HH);
  fetchHits($sqfile, $do_skipaln, "predicted", \@mdl_A, \@seq_accn_A, \%totlen_H, \%p_start_HH, \%p_stop_HH, \%p_strand_HH, $out_root, \%pred_fafile_H);
}

####################################################################
# COMBINE MULTI-EXON SEQUENCES INTO SINGLE CDS/MAT_PEPTIDE SEQUENCES
####################################################################
my @pred_mft_fafile_A = (); # array of predicted CDS sequence files, filled below
wrapperCombineExonsIntoCDS($nmdl, $dir, "predicted", \@mdl_A, \@mdl2mft_map_A, \@mdl_is_first_A, \@mdl_is_final_A, \@pred_mft_fafile_A);

###########################################################################
# CORRECT PREDICTIONS 
# - look for internal starts and stops and make corrections based on those
# - combine new exons into new CDS
###########################################################################
my $source_accn;   # name of CDS/mat_peptide sequence that was translated
my $source_length; # length of CDS/mat_peptide sequence that was translated
my $coords_from;   # start nt coordinate of current translation
my $coords_to;     # stop nt coordinate of current translation
my @corr_mft_start_AH = ();  # [0..$i..nmft-1], each element is a hash with keys $key as sequence accessions and values 
                             # are number of nucleotides that the prediction of the start coordinate should be corrected
                             # based on an esl-translate translation of the predicted CDS/mat_peptide sequence, values can be negative
                             # or positive
my @corr_mft_stop_AH = ();   # [0..$i..nmft-1], each element is a hash with keys $key as sequence accessions and values 
                             # are number of nucleotides that the prediction of the stop coordinate should be corrected
                             # based on an esl-translate translation of the predicted CDS/mat_peptide sequence, values can be negative
                             # or positive
my @coords_len_AH = ();  # [0..$i..nmft-1], each element is a hash with keys $key as sequence accessions and values 
                         # are lengths of translated protein sequences (in nucleotides corresponding to values in $corr_mft_start_AH
                         # and $corr_mft_stop_AH
my @did_corr_exon_start_AH = ();  # [0..$i..ref_nexons-1], each element is a hash with keys $key as sequence accessions and values 
                                  # are '1' if this exon's start position was corrected
my @did_corr_exon_stop_AH = ();   # [0..$i..ref_nexons-1], each element is a hash with keys $key as sequence accessions and values 
                                  # are '1' if this exon's stop position was corrected
my %c_start_HH = ();        # corrected start positions of hits, start with a copy of p_start_HH
my %c_stop_HH  = ();        # corrected stop positions of hits,  start with a copy of p_stop_HH
my %corr_fafile_H = ();     # hash of names of fasta files for corrected exon sequences, keys: model name from @mdl_A, value: name of file
my @corr_mft_fafile_A = (); # array of corrected CDS/mat_peptide sequence files, filled below

if((! $do_nocorrect) && (! $do_matpept)) { 
  ($seconds, $microseconds) = gettimeofday();
  $start_time = ($seconds + ($microseconds / 1000000.));
  printf("%-65s ... ", sprintf("# Translating predicted %s to identify internal starts/stops", ($do_matpept) ? "mat_peptides" : "CDS"));

  # translate predicted CDS/mat_peptide sequences using esl-translate to see if any corrections to predictions are necessary
  for(my $c = 0; $c < $ref_nmft; $c++) { 
    # determine first and final model for this feature

    my $cur_fafile = $pred_mft_fafile_A[$c];
    # translate into AA sequences
    my $tmp_esl_translate_output = $cur_fafile;
    if($do_matpept) { 
      $tmp_esl_translate_output    =~ s/\.mp/\.esl-translate/;
    }
    else { 
      $tmp_esl_translate_output    =~ s/\.cds/\.esl-translate/;
    }
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
      if((! exists $corr_mft_start_AH[$c]{$source_accn}) || ($coords_len > $coords_len_AH[$c]{$source_accn})) { 
        $corr_mft_start_AH[$c]{$source_accn} = $coords_from - 1;
        $corr_mft_stop_AH[$c]{$source_accn}  = -1 * (($source_length - 3) - $coords_to); # source length includes stop codon
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
    for(my $c = 0; $c < $ref_nmft; $c++) { 
      if(! exists $corr_mft_start_AH[$c]{$accn}) { 
        die sprintf("ERROR no corrected start position value for main feature %d sequence $accn", $c+1); 
      }
      if(! exists $corr_mft_stop_AH[$c]{$accn}) { 
        die sprintf("ERROR no corrected stop position value for main feature %d sequence $accn", $c+1); 
      }
      my $corr_start = $corr_mft_start_AH[$c]{$accn};
      my $corr_stop  = $corr_mft_stop_AH[$c]{$accn};
      # sanity check
      if($corr_start < 0) { die "ERROR corrected start less than 0, can't deal with that yet";  }
      if($corr_stop  > 0) { die "ERROR corrected stop greather than 0, can't deal with that yet"; }

      my $first_mdl = $mft2first_mdl_A[$c];
      my $final_mdl = $mft2final_mdl_A[$c];
      my $cur_nmdl = $final_mdl - $first_mdl + 1;

      if($first_mdl == $final_mdl) {
        # easy case: single exon CDS/mat_peptide:
        my $mdl = $mdl_A[$first_mdl];
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
        for(my $h = $first_mdl; $h <= $final_mdl; $h++) { 
          my $hp = $h - $first_mdl;
          $exon_len_A[$hp] = ($p_strand_HH{$mdl_A[$h]}{$seq_accn} eq "+") ? 
              ($p_stop_HH{$mdl_A[$h]}{$seq_accn}  - $p_start_HH{$mdl_A[$h]}{$seq_accn} + 1) : 
              ($p_start_HH{$mdl_A[$h]}{$seq_accn} - $p_stop_HH{$mdl_A[$h]}{$seq_accn}  + 1);
          $cds_len += $exon_len_A[$hp];
          
          # IMPORTANT: copy predicted to corrected, we may 'correct' this below, but if we don't the prediction stays the same.
          #            In case when we correct a multi-exon gene by putting the stop in the non-final exon or the start in
          #            non-first exon, note that this means we will only 'correct' the exons that have the modified stop
          #            or start, and NOT any others. For example if we have a 2-exon CDS and our correction puts the stop
          #            in the first exon, we will not change the annotation of the second exon.
          $mdl = $mdl_A[$h];
          $c_start_HH{$mdl}{$seq_accn} = $p_start_HH{$mdl}{$seq_accn};
          $c_stop_HH{$mdl}{$seq_accn}  = $p_stop_HH{$mdl}{$seq_accn};
        }        

        my $len_so_far; # total length of all exons so far
        my $found_exon; # set to '1' when we find the correct exon to put the corrected start in
        if($corr_start > 0) { # determine which exon the corrected start is in, and correct it:
          $len_so_far = 0; 
          $found_exon = 0; 
          for(my $h = $first_mdl; $h <= $final_mdl; $h++) {
            if(! $found_exon) { 
              my $hp = $h - $first_mdl;
              $mdl = $mdl_A[$h];
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
                #printf("p_start_HH{$mdl}{$seq_accn} exon %d: $p_start_HH{$mdl}{$seq_accn}\n", ($h-$first_mdl+1));
                #printf("c_start_HH{$mdl}{$seq_accn} exon %d: $c_start_HH{$mdl}{$seq_accn}\n", ($h-$first_mdl+1));
              }
            }
          }
          if(! $found_exon) { die sprintf("ERROR unable to find proper exon for corrected start for %s: %d, accn: $accn\n", ($do_matpept) ? "mat_peptide" : "CDS", $c+1); }
        } # end of 'if($corr_start > 0)'
        
        if($corr_stop < 0) { # determine which exon the corrected stop is in, and correct it:
          $len_so_far = 0; 
          $found_exon = 0; 
          for(my $h = $final_mdl; $h >= $first_mdl; $h--) { 
            if(! $found_exon) {
              my $hp = $h - $first_mdl;
              $mdl = $mdl_A[$h];
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
                #printf("p_stop_HH{$mdl}{$seq_accn} exon %d: $p_stop_HH{$mdl}{$seq_accn}\n", ($h-$first_mdl+1));
                #printf("c_stop_HH{$mdl}{$seq_accn} exon %d: $c_stop_HH{$mdl}{$seq_accn}\n", ($h-$first_mdl+1));
              }
            }
          }
          if(! $found_exon) { die sprintf("ERROR unable to find proper exon for corrected stop for %s: %d, accn: $accn\n", ($do_matpept) ? "mat_peptide" : "CDS", $c+1); }
        } # end of 'if($corr_stop < 0)'
      } # end of 'else' entered if we're a multi-exon CDS
    } # end of loop over CDS
  } # end of loop over sequences
  ($seconds, $microseconds) = gettimeofday();
  $stop_time = ($seconds + ($microseconds / 1000000.));
  printf("done. [%.1f seconds]\n", ($stop_time - $start_time));

  # fetch corrected hits into new files
  fetchHits($sqfile, $do_skipaln, "corrected", \@mdl_A, \@seq_accn_A, \%totlen_H, \%c_start_HH, \%c_stop_HH, \%p_strand_HH, $out_root, \%corr_fafile_H);

  # combine multi-exon sequences into CDS:
  wrapperCombineExonsIntoCDS($nmdl, $dir, "corrected", \@mdl_A, \@mdl2mft_map_A, \@mdl_is_first_A, \@mdl_is_final_A, \@corr_mft_fafile_A);
} # end of if(! $do_nocorrect)

#########################################
# TRANSLATE PREDICTIONS INTO PROTEIN SEQS
#########################################
my @nfullprot_A   = ();  # [0..$i..nmft-1], number of accessions we have a full protein for, for CDS $i
my @fullprot_AH   = ();  # [0..$i..nmft-1], each element is a hash with keys $key as sequence accessions and values 
                         # of number of full length protein sequences we have for CDS/mat_peptide $i for key $key. Values should
                         # always be '1', more than '1' is an error, and we never create a value of '0'.
my @ntruncprot_A   = (); # [0..$i..nmft-1], number of accessions we do not have a full protein for, for CDS $i
my @truncprot_AH   = (); # [0..$i..nmft-1], each element is a hash with keys $key as sequence accessions and values 
                          # of number of non-full length protein sequences we have for CDS/mat_peptide $i for key $key. Values should
                          # always be '1', more than '1' is an error, and we never create a value of '0'.
my @aa_full_files_A = (); # array of protein sequence files we are about to create
for(my $c = 0; $c < $ref_nmft; $c++) { 
  my $cur_fafile = ($do_nocorrect || $do_matpept) ? $pred_mft_fafile_A[$c] : $corr_mft_fafile_A[$c];
  # translate into AA sequences
  my $tmp_aa_fafile   = $cur_fafile;
  my $aa_full_fafile  = $cur_fafile;
  my $aa_trunc_fafile = $cur_fafile;
  if($do_matpept) { 
    $tmp_aa_fafile   =~ s/\.mp/\.aa.tmp/;
    $aa_full_fafile  =~ s/\.mp/\.aa.full/;
  }
  else { 
    $tmp_aa_fafile   =~ s/\.cds/\.aa.tmp/;
    $aa_full_fafile  =~ s/\.cds/\.aa.full/;
  }
  # $aa_trunc_fafile    =~ s/\.mft/\.aa.trunc/;
  
  # --watson specifies only translate top strand (not reverse strand)
  # -m specifies only AUG allowed as start
  $cmd = $esl_translate . " --watson -m $cur_fafile > $tmp_aa_fafile";

  runCommand($cmd, 0);

  # now we have to parse that file to only keep the full length protein seqs
  # we also keep track of which sequence accessions we have full length proteins for
  my $prot_must_start_at_posn_1 = 1; # we only want to fetch sequences that start at position 1
  my $prot_must_stop_at_posn_L  = 1; # we only want to fetch sequences that stop  at position L (final position, that is a valid stop exists there)
  parseEslTranslateOutput($tmp_aa_fafile, $aa_full_fafile, $prot_must_start_at_posn_1, $prot_must_stop_at_posn_L, \%{$fullprot_AH[$c]}, \$nfullprot_A[$c]);

  push(@aa_full_files_A, $aa_full_fafile);

  #printf("CDS: %d nfullprot: %d\n", ($c+1), $nfullprot_A[$c]);
  
  ## now fetch any protein sequences that start at position 1 but do not end at the final predicted position
  #$prot_must_start_at_posn_1 = 1; # we only want to fetch sequences that start at position 1
  #$prot_must_stop_at_posn_L  = 0; # we only want to fetch sequences that stop  at position L (final position, that is a valid stop exists there)
  #parseEslTranslateOutput($tmp_aa_fafile, $aa_trunc_fafile, $prot_must_start_at_posn_1, $prot_must_stop_at_posn_L, \%{$truncprot_AH[$c]}, \$ntruncprot_A[$c]);
  #printf("CDS: %d ntruncprot: %d\n", ($c+1), $ntruncprot_A[$c]);
}

#################################
# CREATE MULTIPLE DNA ALIGNMENTS
#################################
if($do_hmmer) { 
  alignHits($hmmalign, $hmmfetch, $model_db, $do_skipaln, \@mdl_A, \@seq_accn_A, 
            (($do_nocorrect || $do_matpept) ? \%pred_fafile_H : \%corr_fafile_H), 
            \%p_start_HH, \%p_fid2ref_HH, \%p_refdel_HHA, \%p_refins_HHA);
}
else { 
  alignHits($cmalign, $cmfetch, $model_db, $do_skipaln, \@mdl_A, \@seq_accn_A, 
            (($do_nocorrect || $do_matpept) ? \%pred_fafile_H : \%corr_fafile_H), 
            \%p_start_HH, \%p_fid2ref_HH, \%p_refdel_HHA, \%p_refins_HHA);
}

#####################################
# CREATE MULTIPLE PROTEIN ALIGNMENTS
#####################################
# Make sure that the predictions in the reference of each feature are the same as the GenBank annotation
# If they're not we can't use the current method which fetches the CDS/mat_peptide sequence based on the
# predicion
if(! $do_skipaln) {
  ($seconds, $microseconds) = gettimeofday();
  $start_time = ($seconds + ($microseconds / 1000000.));
  printf("%-65s ... ", sprintf("# Creating multiple alignments of protein sequences"));
  my @tmp_ref_act_exon_starts_AA = (); # [0..$nmft-1][0..$nexons-1] start positions of actual annotations of exons for ref accn, $nexons is main-feature (CDS or mat_peptide) specific
  my @tmp_ref_act_exon_stops_AA  = (); # [0..$nmft-1][0..$nexons-1] stop  positions of actual annotations of exons for ref accn, $nexons is main-feature (CDS or mat_peptide) specific
  my $tmp_ref_tot_nexons = 0;
  my $tmp_ref_nmft       = 0;
  
  getActualAnnotations($accn_A[0], \%mft_tbl_HHA, \@tmp_ref_act_exon_starts_AA, \@tmp_ref_act_exon_stops_AA, \$tmp_ref_tot_nexons, \$tmp_ref_nmft);
  for(my $h = 0; $h < $nmdl; $h++) { 
    my $model   = $mdl_A[$h];
    my $mft_i   = $mdl2mft_map_A[$h];
    my $exon_i  = $mft2exon_map_A[$h];
    if(! exists $p_start_HH{$model}{$ref_seq_accn}) { die "ERROR no prediction in reference for feature $mft_i exon $exon_i.\n"; }
    my $start  = ($do_nocorrect || $do_matpept) ? $p_start_HH{$model}{$ref_seq_accn} : $c_start_HH{$model}{$ref_seq_accn};
    my $stop   = ($do_nocorrect || $do_matpept) ? $p_stop_HH{$model}{$ref_seq_accn}  : $c_stop_HH{$model}{$ref_seq_accn};
    
    my ($ref_start_match, $ref_stop_match) = checkStrictBoundaryMatch(\@tmp_ref_act_exon_starts_AA, \@tmp_ref_act_exon_stops_AA, $mft_i, $exon_i, $start, $stop);
    if(! $ref_start_match) { die "ERROR, predicted reference feature $mft_i exon $exon_i does not match GenBank annotation (start position mismatch)."; }
    if(! $ref_stop_match)  { die "ERROR, predicted reference feature $mft_i exon $exon_i does not match GenBank annotation (stop position mismatch)."; }
    
    if($mdl_is_final_A[$h]) { 
      # fetch the protein sequence
      my $prot_fafile       = $aa_full_files_A[$mft_i];
      my $prot_stkfile      = $prot_fafile; 
      my $prot_hmmfile      = $prot_fafile; 
      my $prot_hmmbuildfile = $prot_fafile; 
      my $prot_alnfile      = $prot_fafile; 
      
      $prot_stkfile         =~ s/\.fa$/\.ref.stk/;
      $prot_hmmfile         =~ s/\.fa$/\.hmm/;
      $prot_hmmbuildfile    =~ s/\.fa$/\.hmmbuild/;
      $prot_alnfile         =~  s/\.fa$/\.stk/;
      
      my $prot_sqfile  = Bio::Easel::SqFile->new({ fileLocation => $aa_full_files_A[$mft_i] });
      my $ref_prot_str = $prot_sqfile->fetch_consecutive_seqs(1, "", -1, undef);
      my ($ref_prot_name, $ref_prot_seq) = split(/\n/, $ref_prot_str);
      $ref_prot_name =~ s/^\>//; # remove fasta header line '>'
      
      # write it out to a new stockholm alignment file
      open(OUT, ">", $prot_stkfile) || die "ERROR, unable to open $prot_stkfile for writing.";
      print OUT ("# STOCKHOLM 1.0\n");
      print OUT ("$ref_prot_name $ref_prot_seq\n");
      print OUT ("//\n");
      close OUT;
      
      # build an HMM from this single sequence alignment:
      my $cmd = "$hmmbuild $prot_hmmfile $prot_stkfile > $prot_hmmbuildfile";
      runCommand($cmd, 0);
      
      # align all sequences to this HMM
      $cmd = "$hmmalign $prot_hmmfile $prot_fafile > $prot_alnfile";
      runCommand($cmd, 0);
    }
  }
  ($seconds, $microseconds) = gettimeofday();
  $stop_time = ($seconds + ($microseconds / 1000000.));
  printf("done. [%.1f seconds]\n", ($stop_time - $start_time));
}

#########################
# OUTPUT ANNOTATION TABLE
#########################
print "#\n";

my @out_col_header_AA = (); # used only if $do_seqrow
                            # ref to 2D array of output tokens for column headers
                            # first dim: [0..4], 0, 1, and 3 are arrays of header tokens, 2 and 4 are dashed lines that are
my @out_row_header_A  = (); # used only if $do_seqcol
                            # ref to array of output tokens for column or row headers
my @out_header_exp_A  = (); # same size of 1st dim of @out_col_header_AA and only dim of @out_row_header_A
                            # explanations of each header
getHeadings($do_seqrow, $do_seqcol, $do_matpept, $do_nofid, $do_nomdlb, $do_noss3, $do_nostop, $do_fullolap, $do_fulladj, 
            $origin_seq, $ref_tot_nexons, $nmdl, \@mdl2mft_map_A, \@mft2exon_map_A, \@mdl_is_first_A, \@mdl_is_final_A, 
            \@mft_out_short_A, \@mft_out_product_A, ($do_seqrow) ? (\@out_col_header_AA) : undef,
            ($do_seqcol) ? (\@out_row_header_A)  : undef,
            (\@out_header_exp_A));


if($do_seqrow) { # output sequences as rows 
  for(my $i = 0; $i < 5; $i++) { 
    print "#";
    my $ncols = scalar(@{$out_col_header_AA[$i]});
    for(my $j = 0; $j < $ncols; $j++) { 
      print $out_col_header_AA[$i][$j];
    }
    print "\n";
  }
}

########################################################################
# Pass through all accessions, and gather and output annotation for each
########################################################################
my @ref_ol_AA     = (); # 2D array that describes the overlaps in the reference, $ref_ol_AA[$i][$j] is '1' if the exons modeled by model $i and $j overlap
my @ref_adj_AA    = (); # 2D array that describes the adjacencies in the reference, $ref_ol_AA[$i][$j] is '1' if the exons modeled by model $i and $j are adjacent,
                        # only used if $do_matpept
my @ref_prv_adj_A = (); # [0..$h..$nmdl-1]: $i, mdl $i is previous and adjacent to $h, $i and $h are both primary peptides
my @ref_nxt_adj_A = (); # [0..$h..$nmdl-1]: $i, mdl $i is after    and adjacent to $h, $i and $h are both primary peptides
my @ref_start_A   = (); # [0..$h..$nmdl-1]: predicted start  for reference for model $h
my @ref_stop_A    = (); # [0..$h..$nmdl-1]: predicted stop   for reference for model $h
my @ref_strand_A  = (); # [0..$h..$nmdl-1]: predicted strand for reference for model $h
my @ref_len_A     = (); # [0..$h..$nmdl-1]: predicted length for reference for model $h
my $width;              # width of a field

# data structures necessary only for seqs-as-columns output mode (if $do_seqcol) is TRUE
my @ref_out_A   = (); # array of output fields for the reference accession, used if $do_seqcol
my @page_accn_A = (); # [0..$cur_pagesize]2D array, each element is an array of output tokens for one accession
my @page_out_AA = (); # 2D array, [0..$a..$cur_pagesize-1][0..(ntoks-1)] first dimension is of size $cur_pagesize, 
                      # each element is an array of output tokens for accession $page_accn_A[$a]
my $cur_pagesize = 0; # current number of accessions we have info for in page_out_AA (size of first dimension in page_out_AA)
                      # when this hits $nseqcol, we dump the output
my $npages = 0;       # number of pages output, only used if $do_seqcol

for(my $a = 0; $a < $naccn; $a++) { 
  my $accn      = $accn_A[$a];
  my $seq_accn  = $seq_accn_A[$a];
  my @cur_out_A = (); # array of current tokens to print

  # sanity checks
  if(! exists $totlen_H{$accn}) { die "ERROR accession $accn does not exist in the length file $length_file"; }
  
  # Create the initial portion of the output line, the accession and length
  push(@cur_out_A, sprintf("%-5d  ", ($a+1)));
  push(@cur_out_A, sprintf("%-20s  ", $accn)); 
  push(@cur_out_A, sprintf("%6d ", $totlen_H{$accn}));

  #########################################################
  # Get information on the actual annotation of this genome
  #########################################################
  my @act_exon_starts_AA = (); # [0..$nmft-1][0..$nexons-1] start positions of actual annotations of exons for this accn, $nexons is main-feature (CDS or mat_peptide) specific
  my @act_exon_stops_AA  = (); # [0..$nmft-1][0..$nexons-1] stop  positions of actual annotations of exons for this accn, $nexons is main-feature (CDS or mat_peptide) specific
  my $tot_nexons = 0;

  getActualAnnotations($accn, \%mft_tbl_HHA, \@act_exon_starts_AA, \@act_exon_stops_AA, \$tot_nexons, \$nmft);

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
  my $ss3_yes_char    = ".";
  my $ss3_unsure_char = "?";
  my $ss3_no_char     = "!";
  my $hit_length;
  my $at_least_one_fail; # set to '1' for each main feature (CDS or mat_pept) if any of the 'tests' for that feature fail
  my $pass_fail_char; # "P" or "F"
  my $pass_fail_str = "";  # string of pass_fail_chars
  my @cur_name_A   = (); # [0..$h..$nmdl-1]: short name model $h
  my @cur_start_A  = (); # [0..$h..$nmdl-1]: predicted start  for model $h
  my @cur_stop_A   = (); # [0..$h..$nmdl-1]: predicted stop   for model $h
  my @cur_strand_A = (); # [0..$h..$nmdl-1]: predicted strand for model $h
  my @cur_len_A    = (); # [0..$h..$nmdl-1]: predicted length for model $h
  my @cur_ol_AA    = (); # [0..$h..$nmdl-1][0..$hp..$nmdl-1]: '1' if $h and $hp predictions overlap
                         # overlap means that the two predictions are on the same strand and overlap by at least 1 nt
  my @cur_adj_AA   = (); # [0..$h..$nmdl-1][0..$hp..$nmdl-1]: '1' if $h and $hp predictions are 'adjacent' 
                         # adjacent means that min distance between any two nt in each prediction is 1

  ###############################################################
  if(defined $origin_seq) { 
    my ($oseq_ct, $oseq_start, $oseq_stop, $oseq_offset, $oseq_passfail) = getOseqOutput($accn, $origin_seq, $origin_offset, \%origin_coords_HA);
    push(@cur_out_A, sprintf("%2d ", $oseq_ct));
    push(@cur_out_A, sprintf("%5s ", $oseq_start));
    push(@cur_out_A, sprintf("%5s ", $oseq_stop));
    push(@cur_out_A, sprintf("%5s ", $oseq_offset));
    push(@cur_out_A, sprintf(" %s", $oseq_passfail));
    $pass_fail_str .= $oseq_passfail;
  }
  ###############################################################

  # now the per-exon predictions:
  my $tot_fid = 0.; # all fractional identities added together
  my $n_fid = 0;    # number of fractional identities

  # arrays of data structures that pertain to all exons of current feature
  my @cur_mdl_A = (); # array of all model names of exons in current feature
  my $cur_nexons  = 0;  # number of exons in current feature

  # first, get starts, stops and strands for all models, so we can determine
  # adjacencies if necessary (if $do_matpept)
  for(my $h = 0; $h < $nmdl; $h++) { 
    my $model  = $mdl_A[$h];
    my ($start, $stop, $strand);
    if(exists $p_start_HH{$model}{$seq_accn}) { 
      $start  = ($do_nocorrect || $do_matpept) ? $p_start_HH{$model}{$seq_accn} : $c_start_HH{$model}{$seq_accn};
      $stop   = ($do_nocorrect || $do_matpept) ? $p_stop_HH{$model}{$seq_accn}  : $c_stop_HH{$model}{$seq_accn};
      $strand = $p_strand_HH{$model}{$seq_accn};
    }
    else { 
      $start  = 0;
      $stop   = 0;
      $strand = "?";
    }
    push(@cur_name_A,   sprintf "%d.%d", $mdl2mft_map_A[$h]+1, $mft2exon_map_A[$h]+1);
    push(@cur_start_A,  $start);
    push(@cur_stop_A,   $stop);
    push(@cur_strand_A, $strand);
    push(@cur_len_A,    (abs($stop-$start) + 1));
    if($a == 0) { 
      @ref_start_A  = @cur_start_A;
      @ref_stop_A   = @cur_stop_A;
      @ref_strand_A = @cur_strand_A;
      @ref_len_A    = @cur_len_A;
    }
  }

  # check for overlaps
  checkForOverlapsOrAdjacencies(0, \@cur_start_A, \@cur_stop_A, \@cur_strand_A, \@cur_ol_AA); # '0': do overlap, not adjacency
  if($a == 0) { 
    @ref_ol_AA = @cur_ol_AA; 
    if($do_matpept) { # determine which models are for primary mature peptides 
      matpeptFindPrimaryPeptides(\@ref_ol_AA, \@ref_len_A, \@mdl_is_primary_A);
    }      
  }
  # check for adjacencies if necessary
  if($do_matpept) {
    # get adjacency information
    checkForOverlapsOrAdjacencies(1, \@cur_start_A, \@cur_stop_A, \@cur_strand_A, \@cur_adj_AA);
    if($a == 0) { 
      @ref_adj_AA = @cur_adj_AA;
      # for all primary peptides: we assume there is at most 1 peptide before it that is adjacent 
      # and at most 1 other peptide after it that is adjacent, determine what those are and 
      # store them
      matpeptFindPrimaryAdjacencies(\@ref_adj_AA, \@mdl_is_primary_A, \@ref_prv_adj_A, \@ref_nxt_adj_A);
    }
  }

  # handle 5' UTR, if $do_matpept
  if($do_matpept) { 
    if($cur_start_A[0] == 0) { # no annotated first matpept, 
      push(@cur_out_A, sprintf("  %6s", "?")); # start
      push(@cur_out_A, sprintf(" %6s", "?"));  # stop
      push(@cur_out_A, sprintf(" %6s", "?"));  # length
    }
    elsif($cur_start_A[0] == 1) { # 1st matpept starts at nt 1!
      push(@cur_out_A, sprintf("  %6d", 0)); # start 
      push(@cur_out_A, sprintf("  %6d", 0)); # stop 
      push(@cur_out_A, sprintf("  %6d", 0)); # length
    }
    else { # 1st matpept does not start at nt 1!
      push(@cur_out_A, sprintf("  %6d", 1)); # start 
      push(@cur_out_A, sprintf("  %6d", $cur_start_A[0]-1)); # stop 
      push(@cur_out_A, sprintf("  %6d", $cur_start_A[0]-1)); # length
    }
  }    

  # now go through each model and collect output tokens
  for(my $h = 0; $h < $nmdl; $h++) { 
    my $model   = $mdl_A[$h];
    my $mft_i   = $mdl2mft_map_A[$h];
    my $exon_i  = $mft2exon_map_A[$h];

    if($mdl_is_first_A[$h]) {
      # reset these
      $hit_length = 0; 
      $at_least_one_fail = 0;
      @cur_mdl_A = ();
      $cur_nexons = 0;
    }
    push(@cur_mdl_A, $model);
    $cur_nexons++;

    if($predicted_string ne "") { $predicted_string .= "  "; }
    if(exists $p_start_HH{$model}{$seq_accn}) { 
      my $start_corrected_exon = (exists $did_corr_exon_start_AH[$h]{$accn}) ? 1 : 0;
      my $stop_corrected_exon  = (exists $did_corr_exon_stop_AH[$h]{$accn})  ? 1 : 0;
      my $start_corrected_mft = 0; # possibly redefined below
      my $stop_corrected_mft  = 0; # possibly redefined below
      if((! $do_nocorrect) && (! $do_matpept)) { 
        $start_corrected_mft  = ($corr_mft_start_AH[$mdl2mft_map_A[$h]]{$accn} != 0) ? 1 : 0;
        $stop_corrected_mft   = ($corr_mft_stop_AH[$mdl2mft_map_A[$h]]{$accn}  != 0) ? 1 : 0;
      }
      my $start    = ($do_nocorrect || $do_matpept) ? $p_start_HH{$model}{$seq_accn} : $c_start_HH{$model}{$seq_accn};
      my $stop     = ($do_nocorrect || $do_matpept) ? $p_stop_HH{$model}{$seq_accn}  : $c_stop_HH{$model}{$seq_accn};
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
          checkStrictBoundaryMatch   (\@act_exon_starts_AA, \@act_exon_stops_AA, $mft_i, $exon_i, $start, $stop) :
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

      # TODO: MODIFY ANNOTATION FOR EXONS WITHOUT CORRECTED STARTS OR STOPS IN FEATURES WITH
      #       OTHER EXONS THAT HAVE CORRECTED STARTS OR STOPS
      push(@cur_out_A, sprintf("  %8s ", ($start_match ? " " . $start . " " : "[" . $start . "]")));
      push(@cur_out_A, sprintf("%8s",  ($stop_match  ? " " . $stop .  " " : "[" . $stop . "]")));
      if(! $do_nofid) { 
        push(@cur_out_A, sprintf(" %5.3f", $p_fid2ref_HH{$model}{$seq_accn}));
      }
      $tot_fid += $p_fid2ref_HH{$model}{$seq_accn};
      $n_fid++;

      if(! $do_nomdlb) { 
        push(@cur_out_A, "  " . $hang5 . $hang3);
      }        

      # overlap string for this feature
      if(! $do_noolap) { 
        my $ol_pf_char = undef;
        my $ol_str     = undef;
        ($ol_pf_char, $ol_str) = compareOverlapsOrAdjacencies(\@cur_name_A, "/", $h, \@ref_ol_AA, \@cur_ol_AA);
        push(@cur_out_A, $ol_str);
        if($ol_pf_char eq "F") { 
          $at_least_one_fail = 1;
        }
      }

      # adjacency string for this feature
      if($do_matpept) { 
        my $adj_pf_char = undef;
        my $adj_str     = undef;
        ($adj_pf_char, $adj_str) = compareOverlapsOrAdjacencies(\@cur_name_A, "|", $h, \@ref_adj_AA, \@cur_adj_AA);
        push(@cur_out_A, $adj_str);
        if($adj_pf_char eq "F") { 
          $at_least_one_fail = 1;
        }
      }

      if($mdl_is_first_A[$h]) { # determine $start_codon_char
        if($p_strand_HH{$model}{$seq_accn} eq "-") { 
          $start_codon_posn = (($start-2) < 0) ? $start + $totlen_H{$accn} + 1 : $start;
        }
        else { 
          $start_codon_posn = ($start < 0) ? $start + $totlen_H{$accn} + 1 : $start;
        }
        $start_codon = fetchCodon($sqfile, $seq_accn, $start_codon_posn, $p_strand_HH{$model}{$seq_accn});
        if($do_matpept) { 
          if($ref_prv_adj_A[$h] != -1) { # in reference, primary mat_peptide $prv_adj_A[$h] is adjacent to this one
            $start_codon_char = ($cur_adj_AA[$ref_prv_adj_A[$h]][$h]) ? $ss3_yes_char : $ss3_no_char;
          }
          elsif($h == 0) { # first mat_peptide, look for start codon
            $start_codon_char = ($start_codon eq "ATG") ? $ss3_yes_char : $ss3_no_char;
          }
          else { 
            $start_codon_char = $ss3_unsure_char; # not adjacent, and not supposed to have an ATG that we know of
          }
        }
        else { # ! $do_matpept
          $start_codon_char = ($start_codon eq "ATG") ? $ss3_yes_char : $ss3_no_char;
        }
        if($start_codon_char eq $ss3_no_char) { 
          $at_least_one_fail = 1;
        }
      }
      
      if($mdl_is_final_A[$h]) { 
        if($p_strand_HH{$model}{$seq_accn} eq "-") { 
          $stop_codon_posn    = ($stop < 0) ? ($stop + $totlen_H{$accn}) + 1 + 2 : $stop + 2;
        }
        else { 
          $stop_codon_posn = (($stop-2) < 0) ? ($stop + $totlen_H{$accn}) + 1 - 2 : $stop - 2;
        }
        $stop_codon = fetchCodon($sqfile, $seq_accn, $stop_codon_posn, $p_strand_HH{$model}{$seq_accn});

        if($do_matpept) { 
          if($ref_nxt_adj_A[$h] != -1) { # in reference, primary mat_peptide $nxt_adj_A[$h] is adjacent to this one
            $stop_codon_char = ($cur_adj_AA[$h][$ref_nxt_adj_A[$h]]) ? $ss3_yes_char : $ss3_no_char;
          }
          # TODO: reexamine this, for Dengue, final matpept encoding seq DOES NOT end with a valid stop, so don't enforce it
          #elsif($h == ($nmdl-1)) { # final mat_peptide, look for stop codon
          #$stop_codon_char = ($stop_codon eq "TAG" || $stop_codon eq "TAA" || $stop_codon eq "TGA") ? $ss3_yes_char : $ss3_no_char;
        #}
          else { 
            $stop_codon_char = $ss3_unsure_char; # not adjacent, and not supposed to have a stop that we know of
          }
        }
        else { # ! $do_matpept
          $stop_codon_char = ($stop_codon eq "TAG" || $stop_codon eq "TAA" || $stop_codon eq "TGA") ? $ss3_yes_char : $ss3_no_char;
        }
        if($stop_codon_char eq $ss3_no_char) { 
          $at_least_one_fail = 1;
        }

        $multiple_of_3_char = (($hit_length % 3) == 0) ? $ss3_yes_char : $ss3_no_char;
        if($multiple_of_3_char eq $ss3_no_char) { 
          $at_least_one_fail = 1;
        }

        # append the ss3 (start/stop/multiple of 3 info)
        push(@cur_out_A, sprintf(" %6d", $hit_length));
        if(! $do_noss3) { 
          push(@cur_out_A,  sprintf(" %s%s%s", $start_codon_char, $stop_codon_char, $multiple_of_3_char));
        }
        if(! $do_nostop) { 
          push(@cur_out_A, sprintf(" %3s", $stop_codon));
        }

        if($at_least_one_fail) { 
          $pass_fail_char = "F"; 
          #### ensure we didn't fetch a full length protein for this feature
          ###if(exists $fullprot_AH[$mft_i]{$accn}) { 
          ###  die sprintf("ERROR, incorrectly translated full length protein for feature: %d for seq accn: $accn, but at least one test failed...", $mft_i+1); 
          ###}
          ###if(! exists $truncprot_AH[$mft_i]{$accn}) { 
          ###  die sprintf("ERROR, failed to translate truncated protein for feature: %d for seq accn: $accn, but at least one test failed...", $mft_i+1); 
          ###}
        }
        else { 
          $pass_fail_char = "P"; 
          #### verify that we have a translated protein sequence for this feature
          ###if(! exists $fullprot_AH[$mft_i]{$accn}) { 
          ###  die sprintf("ERROR, failed to translate full length protein for feature: %d for seq accn: $accn, but all tests passed...", $mft_i+1); 
          ###}
          ###if(exists $truncprot_AH[$mft_i]{$accn}) { 
          ###  die sprintf("ERROR, incorrectly translated truncated protein for feature: %d for seq accn: $%accn, but all tests passed...", $mft_i+1);
          ###}
        }
        $pass_fail_char = ($at_least_one_fail) ? "F" : "P";
        push(@cur_out_A, sprintf(" %2s", $pass_fail_char));
        $pass_fail_str .= $pass_fail_char;
      }
    }
    else { 
      # printf("no hits for $model $seq_accn\n");
      $at_least_one_fail = 1;
      push(@cur_out_A, sprintf("  %8s ", "NP")); # start position
      push(@cur_out_A, sprintf("%8s",  "NP"));   # stop position
      if(! $do_nofid) { 
        push(@cur_out_A, sprintf(" %5s", "NP")); # fid
      }        
      if(! $do_nomdlb) { 
        push(@cur_out_A, "  NP"); # model boundaries
      }
      if($mdl_is_final_A[$h]) { 
        push(@cur_out_A, sprintf(" %6s", "NP")); # length
        if(! $do_noss3) { 
          push(@cur_out_A, "  NP"); # ss3
        }
        if(! $do_nostop) { 
          push(@cur_out_A, sprintf(" %3s", "NP")); # stop
        }
        $pass_fail_char = "F";
        push(@cur_out_A, sprintf(" %2s", $pass_fail_char));
        $pass_fail_str .= $pass_fail_char;
      }
    }
  }

  # handle 3' UTR, if $do_matpept
  if($do_matpept) { 
    if($cur_stop_A[($nmdl-1)] == 0) { # no annotated final matpept
      push(@cur_out_A, sprintf("  %6s", "?")); # start
      push(@cur_out_A, sprintf(" %6s", "?"));  # stop
      push(@cur_out_A, sprintf(" %6s", "?"));  # length
    }
    elsif($cur_stop_A[($nmdl-1)] == $totlen_H{$accn}) { # final matpept stops at final nt!
      push(@cur_out_A, sprintf("  %6d", 0)); # start 
      push(@cur_out_A, sprintf("  %6d", 0)); # stop 
      push(@cur_out_A, sprintf("  %6d", 0)); # length
    }
    else { # final matpept does not stop at final nt
      push(@cur_out_A, sprintf("  %6d", $cur_stop_A[($nmdl-1)]+1)); # start 
      push(@cur_out_A, sprintf("  %6d", $totlen_H{$accn})); # stop 
      push(@cur_out_A, sprintf("  %6d", ($totlen_H{$accn} - ($cur_stop_A[($nmdl-1)]+1) + 1))); # length
    }
  }    

  # total length
  push(@cur_out_A, sprintf("  %6d", $totlen_H{$accn}));

  # average fid
  push(@cur_out_A, sprintf("  %5.3f", $tot_fid / $n_fid));

  # output number of actually annotated features and summed total of exons in those features, if nec
  if(! $do_noexist) { 
    push(@cur_out_A, sprintf("  %5d", $nmft));
    push(@cur_out_A, sprintf("  %5d", $tot_nexons));
    push(@cur_out_A, sprintf("  %5d", $nmatch_boundaries));
  }

  # output overlap info, if nec
  if($do_fullolap) { 
    my $ol_pass_fail_char = undef;
    my $overlap_notes = undef;
    ($ol_pass_fail_char, $overlap_notes) = compareOverlapsOrAdjacencies(\@cur_name_A, "/", undef, \@ref_ol_AA, \@cur_ol_AA);
    push(@cur_out_A, sprintf("  %20s", $overlap_notes));
    $pass_fail_str .= $ol_pass_fail_char;
  }
  
  # output adjacency info, if nec
  if($do_fulladj) { 
    my $adj_pass_fail_char;
    my $adjacency_notes = undef;
    ($adj_pass_fail_char, $adjacency_notes) = compareOverlapsOrAdjacencies(\@cur_name_A, "|", undef, \@ref_adj_AA, \@cur_adj_AA);
    $pass_fail_str .= $adj_pass_fail_char;
    push(@cur_out_A, sprintf("  %20s", $adjacency_notes));
  }

  my $result_str = ($pass_fail_str =~ m/F/) ? "FAIL" : "PASS";
  $result_str .= " " . $pass_fail_str;
  push(@cur_out_A, sprintf("  %s", $result_str));

  if($do_seqrow) { 
    foreach my $el (@cur_out_A) { 
      print $el;
    }
    print "\n";
  }
  elsif($do_seqcol) { 
    if($a == 0) { 
      # copy reference info if this is the reference
      @ref_out_A = @cur_out_A; 
    } 
    else { 
      push(@page_out_AA, [@cur_out_A]);
      $cur_pagesize++;
    }
    if(($cur_pagesize+1) == $nseqcol) { 
      $npages++;
      outputSeqAsColumnsPage(\@out_row_header_A, \@page_out_AA, \@ref_out_A, $npages);
      @page_out_AA = ();
      $cur_pagesize = 0;
    }
  }
}
# print final page (if non-empty)
if(($do_seqcol) && ($cur_pagesize > 0)) { 
  $npages++;
  outputSeqAsColumnsPage(\@out_row_header_A, \@page_out_AA, \@ref_out_A, $npages);
}

##########################
# OUTPUT EXPLANATORY TEXT 
##########################
if(! $do_noexp) { 
  outputColumnHeaderExplanations(\@out_header_exp_A);
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
outputGapInfo($gap_perseq_all_FH,     $gap_pergap_all_FH,     $do_matpept, 0, 1, 0, 0, \@mdl_A, \@seq_accn_A, \%mdllen_H, \@mdl2mft_map_A, \%p_refdel_HHA, \%p_refins_HHA);
outputGapInfo($gap_perseq_not3_FH,    $gap_pergap_not3_FH,    $do_matpept, 0, 0, 1, 0, \@mdl_A, \@seq_accn_A, \%mdllen_H, \@mdl2mft_map_A, \%p_refdel_HHA, \%p_refins_HHA);
outputGapInfo($gap_perseq_special_FH, $gap_pergap_special_FH, $do_matpept, 0, 0, 0, 1, \@mdl_A, \@seq_accn_A, \%mdllen_H, \@mdl2mft_map_A, \%p_refdel_HHA, \%p_refins_HHA);

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


# Subroutine: validateRefFeaturesAreUnique()
# Purpose:    Validate that all product annotation for all reference features
#             are unique, i.e. there are no two features that have the same value
#             in their CDS:product or mat_peptide:product annotation.
#
# Args:       $ref_nmft:           number of reference features
#             $ref_mft_product_AR: ref to array of product annotations for the $ref_nmft reference feature
#
# Returns:    void
# Dies:       if more than one ref feature have same product annotation.
sub validateRefFeaturesAreUnique {
  my $sub_name  = "validateRefFeaturesAreUnique()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ref_nmft, $ref_mft_product_AR) = @_;

  my %exists_H = ();
  for(my $i = 0; $i < $ref_nmft; $i++) { 
    if(exists $exists_H{$ref_mft_product_AR->[$i]}) { die sprintf("ERROR %s is product value for more than one reference features!", $ref_mft_product_AR->[$i]); }
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
#  printf("done. [$out_root.nmdler and $out_root.tblout]\n");
  printf("done. [%.1f seconds]\n", $secs_elapsed);

  # next, press the HMM DB we just created
  printf("%-65s ... ", "# Running hmmpress");
  $cmd = "$hmmpress $out_root.hmm > $out_root.hmmpress";
  $secs_elapsed = runCommand($cmd, 0);
#  printf("done. [$out_root.nmdler and $out_root.tblout]\n");
  printf("done. [%.1f seconds]\n", $secs_elapsed);

  return;
}

# Subroutine: createCmDb()
# Synopsis:   Create an CM database from a stockholm database file
#             for use with Infernal 1.1.
#             the $cmbuild executable. If $cmcalibrate is defined
#             also run cmcalibrate. 
#
# Args:       $cmbuild:       path to 'cmbuild' executable
#             $cmcalibrate:   path to 'cmcalibrate' executable
#             $cmpress:       path to 'cmpress' executable
#             $cmfetch:       path to 'cmfetch' executable
#             $nmodel:        number of models we're creating/calibrating
#             $do_calib_slow: '1' to calibrate using default parameters instead of
#                             options to make it go much faster
#             $do_calib_farm: '1' to submit calibration job to farm, '0' to do it locally
#             $stk_file:      stockholm DB file
#             $out_root:      string for naming output files
#             $indi_name_AR:  ref to array of individual model names, we only use this if 
#                             $do_calib_farm is true.
#
# Returns:    void
#
sub createCmDb { 
  my $sub_name = "createCmDb()";
  my $nargs_exp = 10;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cmbuild, $cmcalibrate, $cmpress, $cmfetch, $nmodel, $do_calib_slow, $do_calib_farm, $stk_file, $out_root, $indi_name_AR) = @_;

  if(! -s $stk_file)  { die "ERROR in $sub_name, $stk_file file does not exist or is empty"; }

  # remove the binary files, possibly from an earlier cmbuild/cmpress:
  for my $suffix ("i1m", "i1i", "i1f", "i1p") { 
    my $file = $out_root . ".cm." . $suffix;
    if(-e $file) { unlink $file; }
  }

  my ($cmbuild_opts,     $cmbuild_cmd);
  my ($cmcalibrate_opts, $cmcalibrate_cmd);
  my ($cmpress_opts,     $cmpress_cmd);

  # build step:
  $cmbuild_opts = "-F";
  $cmbuild_cmd  = "$cmbuild $cmbuild_opts $out_root.cm $stk_file > $out_root.cmbuild";
  printf("%-65s ... ", sprintf("# Running cmbuild to build %d CMs", $nmodel));
  my $secs_elapsed = runCommand($cmbuild_cmd, 0);
  printf("done. [%.1f seconds]\n", $secs_elapsed);

  # calibration step:
  $cmcalibrate_opts = " --cpu 4 ";
  if(! $do_calib_slow) { $cmcalibrate_opts .= " -L 0.04 "; }
  $cmcalibrate_cmd  = "$cmcalibrate $cmcalibrate_opts $out_root.cm > $out_root.cmcalibrate";
  
  if($do_calib_farm) { 
    # split up model file into individual CM files, then submit a job to calibrate each one, and exit. 
    for(my $i = 0; $i < $nmodel; $i++) { 
      my $cmfetch_cmd = "$cmfetch $out_root.cm $indi_name_AR->[$i] > $out_root.$i.cm";
      runCommand($cmfetch_cmd, 0);
      my $out_tail    = $out_root;
      $out_tail       =~ s/^.+\///;
      my $jobname     = "c." . $out_tail . $i;
      my $errfile     = $out_root . "." . $i . ".err";
      $cmcalibrate_cmd  = "$cmcalibrate $cmcalibrate_opts $out_root.$i.cm > $out_root.$i.cmcalibrate";
      my $farm_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n -l h_rt=288000,h_vmem=8G,mem_free=8G -pe multicore 4 -R y " . "\"" . $cmcalibrate_cmd . "\"\n";
      # print("$farm_cmd\n");
      runCommand($farm_cmd, 0);
    }
    # final step, remove the master CM file, so we can create a new one
    $cmd = "rm $out_root.cm";
    runCommand($cmd, 0);
  }
  else { 
    # calibrate the model
    printf("%-65s ... ", "# Running cmcalibrate");
    $secs_elapsed = runCommand($cmcalibrate_cmd, 0);
    #printf("\n$cmcalibrate_cmd\n");
    printf("done. [%.1f seconds]\n", $secs_elapsed);

    # press the model
    $cmpress_cmd = "$cmpress $out_root.cm > $out_root.cmpress";
    printf("%-65s ... ", "# Running cmpress");
    $secs_elapsed = runCommand($cmpress_cmd, 0);
    #printf("\n$cmpress_cmd\n");
    printf("done [%.1f seconds]\n", $secs_elapsed);
  } # end of 'else' entered if $do_calib_farm is false

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
#             $do_farm:     '1' to submit job to farm, instead of running it locally
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

  my ($cmscan, $do_glocal, $do_skip, $do_farm, $model_db, $seq_fasta, $tblout_file, $stdout_file) = @_;

  if($do_skip) { 
    printf("%-65s ... ", "# Skipping cmscan step");
    if(! -s $tblout_file) { die "ERROR, with -skipscan or -skipaln cmscan output is expected to already exist, but $tblout_file does not or is empty"; }
    if(! -s $stdout_file) { die "ERROR, with -skipscan or -skipaln cmscan output is expected to already exist, but $stdout_file does not or is empty"; }
    printf("done. [-skipscan or -skipaln]\n");
    return;
  }

  my $opts = "";
  if($do_iglocal) { $opts .= "-g "; }
  #$opts .= " --cpu 0 --rfam --tblout $tblout_file --verbose --nohmmonly ";
  # opts: nmdler F1, F2, F2b, F3 and F3b; Infernal --rfam F4, F4b, F5, and F6
  $opts .= " --cpu 0 --F1 0.02 --F2 0.001 --F2b 0.001 --F3 0.00001 --F3b 0.00001 --F4 0.0002 --F4b 0.0002 --F5 0.0002 --F6 0.0001 --tblout $tblout_file --verbose --nohmmonly ";
  if(! defined $stdout_file) { $stdout_file = "/dev/null"; }

  if(! -s $model_db)   { die "ERROR in $sub_name, $model_db file does not exist or is empty"; }
  if(! -s $seq_fasta) { die "ERROR in $sub_name, $seq_fasta file does not exist or is empty"; }

  my $cmd = "$cmscan $opts $model_db $seq_fasta > $stdout_file";
  if(! $do_farm) { 
    # default mode, run job locally
    printf("%-65s ... ", "# Running cmscan");
    my $secs_elapsed = runCommand($cmd, 0);
    printf("done. [%.1f seconds]\n", $secs_elapsed);
  }
  else { 
    # submit job to farm and return
    my $jobname = $seq_fasta;
    my $errfile = $stdout . ".err";
    $jobname =~ s/^.+\///; # remove everything up until final '/'
    my $farm_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n -l h_rt=288000,h_vmem=8G,mem_free=8G " . "\"" . $cmd . "\" > /dev/null\n";
    runCommand($farm_cmd, 0);
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
# Returns:    Two values:
#             alignment length 
#             checksum of alignment
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
  return ($msa->alen, $msa->checksum);
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
      #Maize-streak_r23.NC_001346/Maize-streak_r23.NC_001346.ref.mft.4        -          NC_001346:genome:NC_001346:1:2689:+: -                1     819    2 527    1709    2527    1709     819    -    9.8e-261  856.5  12.1  -
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
      #Maize-streak_r23.NC_001346.ref.mft.4        -         NC_001346:genome-duplicated:NC_001346:1:2689:+:NC_001346:1:2689:+: -          cm        1      819     2527     1709      -    no    1 0.44   0.2  892.0         0 !   -
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
    # if $qseq_posn == -1, then no occurences were found. In this case we don't store 
    # any entry in coords_HAR for this $accn. The caller needs to know what to do
    # if it sees no entry for an $accn

    while($qseq_posn != -1) { 
      $qseq_posn++;
#      print("$accn qseq_posn 0: $qseq_posn\n");
      if($qseq_posn <= $L) { # note: we've just incremented qseq_posn by 1 in prv line so now it is in 1..length($seq) coords, not 0..length($seq)-1
        my $qseq_start = $qseq_posn;
        my $qseq_stop  = $qseq_posn + length($qseq) - 1;
        if($qseq_stop > $L) { 
          $qseq_start -= $L;
          $qseq_start -= 1; # off-by-one issue with negative indexing
          $qseq_stop  -= $L;
        }
        if(! exists $coords_HAR->{$accn}) { # initialize
          @{$coords_HAR->{$accn}} = ();
        }
        push(@{$coords_HAR->{$accn}}, $qseq_start . ":" . $qseq_stop);
        # printf("Found $qseq in $accn at position %d..%d\n", $qseq_start, $qseq_stop);
      }
      if($qseq_posn > $L) { 
        $qseq_posn = -1; 
        # this breaks the while loop because we're searching in a duplicated genome, 
        # and we're into the second copy, no need to keep searching the same seq
      }
      else { 
        $qseq_posn = index($seq, $qseq, $qseq_posn);
      }
    }
  }
  
  return;
}

# Subroutine: checkStrictBoundaryMatch
#
# Synopsis:   Check if a given start..stop boundary set matches the 
#             actual annotation in $act_AAR->[$mft_i][$exon_i]
#             (if that array element even exists).
#
# Args:       $act_start_AAR: ref to 2D array [0..i..$nmft-1][0..e..$nexon-1], start for feature $i+1 exon $e+1
#             $act_stop_AAR:  ref to 2D array [0..i..$nmft-1][0..e..$nexon-1], stop for feature $i+1 exon $e+1
#             $mft_i:         feature index we want to check against
#             $exon_i:        exon index we want to check against
#             $pstart:        predicted start boundary
#             $pstop:         predicted stop boundary
# Returns:    Two values:
#             '1' if $pstart == $act_start_AAR->[$mft_i][$exon_i], else '0'
#             '1' if $pstop  == $act_stop_AAR->[$mft_i][$exon_i], else '0'
#
sub checkStrictBoundaryMatch {
  my $sub_name = "checkStrictBoundaryMatch";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($act_start_AAR, $act_stop_AAR, $mft_i, $exon_i, $pstart, $pstop) = @_;

  my $retval1 = 
      ((exists ($act_start_AAR->[$mft_i])) && 
       (exists ($act_start_AAR->[$mft_i][$exon_i])) && 
       ($pstart == $act_start_AAR->[$mft_i][$exon_i])) ? 
       1 : 0;

  my $retval2 = 
      ((exists ($act_stop_AAR->[$mft_i])) && 
       (exists ($act_stop_AAR->[$mft_i][$exon_i])) && 
       ($pstop == $act_stop_AAR->[$mft_i][$exon_i])) ? 
       1 : 0;

  return ($retval1, $retval2);
}

# Subroutine: checkNonStrictBoundaryMatch
#
# Synopsis:   Check if a given boundary matches any
#             annotation in the 2D array referred to
#             by $act_AAR.
#
# Args:       $act_start_AAR: ref to 2D array [0..i..$nmft-1][0..e..$nexon-1], start for feature $i+1 exon $e+1
#     :       $act_stop_AAR:  ref to 2D array [0..i..$nmft-1][0..e..$nexon-1], stop for feature $i+1 exon $e+1
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
  my $nmft = scalar(@{$act_start_AAR});
  for(my $i = 0; $i < $nmft; $i++) { 
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

  if($origin_seq =~ m/[^ACGT\|\\]/) { die "ERROR with -oseq <s>, <s> can only contain characters A, C, G, T or \|"; }

  my $origin_offset = index($origin_seq, "|");
  if($origin_offset == -1) { 
    die "ERROR with --oseq <s>, <s> must contain a single | character immediately before the nucleotide that should be the first nt of the genome";
  }
  elsif($origin_offset < (length($origin_seq)-1)) { # if this isn't the final character of the string, look for another one
    my $second_offset = index($origin_seq, "|", $origin_offset+1);
    if($second_offset != -1) { 
      die "ERROR with --oseq <s>, <s> must contain a single | character, $origin_seq has more than one";
    }
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
#             the hits for each feaute to files.
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

  printf("%-65s ... ", ($do_skip) ? "# Skipping multiple alignment creation" : "# Creating multiple alignments of DNA sequences");
    
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

# Subroutine: checkForOverlapsOrAdjacencies()
#
# Synopsis:   Given refs to three arrays that describe 
#             a list of hits, check if any of them overlap or
#             are adjacent to each other.
#
# Args:       $do_adj:         '1' to not check for overlaps but rather check for adjacency
#                              $j is adjacent to $i if the minimum distance between any nt
#                              in $i and any nt in $j is exactly 1 nt.
#             $enforce_strand: '1' to enforce that overlaps or adjacencies are only valid if on same strand
#             $start_AR:       ref to array of start positions
#             $stop_AR:        ref to array of stop positions
#             $strand_AR:      ref to array of strands
#             $observed_AAR:   ref to 2D array of observed overlaps or adjacencies
#                              $observed_AAR->[$i][$j] is '1' if 
#                              those two exons are observed to overlap or be adjacecent
#                              FILLED HERE
# 
# Returns:    void, fills $observed_AAR;
#
sub checkForOverlapsOrAdjacencies {
  my $sub_name = "checkForOverlapsOrAdjacencies";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($do_adj, $start_AR, $stop_AR, $strand_AR, $observed_AAR) = @_;

  my $nhits = scalar(@{$start_AR});

  # initialize
  for(my $i = 0; $i < $nhits; $i++) { 
    for(my $j = 0; $j < $nhits; $j++) { 
      $observed_AAR->[$i][$j] = 0;
    }
  }

  my $nres_overlap;  # number of nucleotides two exons overlap by
  my $nres_distance; # distance between two closest nucleotides in two exons
  for(my $i = 0; $i < $nhits; $i++) { 
    my $start_i  = $start_AR->[$i];
    my $stop_i   = $stop_AR->[$i];
    my $strand_i = "+";
    if($start_i > $stop_i) { # swap them
      my $tmp   = $start_i;
      $start_i  = $stop_i;
      $stop_i   = $tmp;
      $strand_i = "-";
    }
    for(my $j = $i+1; $j < $nhits; $j++) { 
      my $start_j  = $start_AR->[$j];
      my $stop_j   = $stop_AR->[$j];
      my $strand_j = "+";
      if($start_j > $stop_j) { # swap them
        my $tmp   = $start_j;
        $start_j  = $stop_j;
        $stop_j   = $tmp;
        $strand_j = "-";
      }
      my $found_ol_or_adj = 0;
      if($strand_i eq $strand_j) { 
        $nres_overlap  = getOverlap ($start_i, $stop_i, $start_j, $stop_j);
        $nres_distance = getDistance($start_i, $stop_i, $start_j, $stop_j);
        if(((! $do_adj) && ($nres_overlap > 0)) ||
           ((  $do_adj) && ($nres_distance == 1))) { 
          $observed_AAR->[$i][$j] = 1;
          $observed_AAR->[$j][$i] = 1;
          # printf("found %s between $i  and $j ($start_i..$stop_i and $start_j..$stop_j)\n", ($do_adj) ? "adjacency" : "overlap");
          $found_ol_or_adj = 1;
        }
      }
      if(! $found_ol_or_adj) { # no overlap or adjacency
        $observed_AAR->[$i][$j] = 0;
        $observed_AAR->[$j][$i] = 0;
      }
    }
  }

  return;
}

# Subroutine: compareOverlapsOrAdjacencies()
#
# Synopsis:   Compare one 2D array describing overlaps or adjacencies
#             to another.
#
# Args:       $name_AR:      ref to array of short names for each annotation
#             $div_char:     divider character
#             $index:        if   defined, the only feature $index we check overlaps/adjacencies for
#                            if ! defined, check overlaps/adjacencies for all features
#             $expected_AAR: ref to 2D array of expected overlaps $expected_ol_AAR->[$i][$j] is '1' if 
#                            those two exons are expected to overlap or be adjacent, PRE-FILLED
#             $observed_AAR: ref to 2D array of test overlaps $observed_ol_AAR->[$i][$j] is '1' if 
#                            those two exons are observed to overlap or be adjacent, PRE-FILLED
#             
# Returns:    Two values:
#             $pass_fail_char: "P" if overlaps/adjacencies match b/t observed and expected, else "F"
#             $notes:          string describing the overlaps/adjacencies, empty string ("") if no overlaps.
#
sub compareOverlapsOrAdjacencies {
  my $sub_name = "compareOverlapsOrAdjacencies";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($name_AR, $div_char, $index, $expected_AAR, $observed_AAR) = @_;

  my $size = scalar(@{$name_AR});
  my $ret_str = "";
  my $nfound = 0;

  # check @observed_AA against @{$expected_AAR}
  my $pass_fail_char = "P";
  for(my $i = 0; $i < $size; $i++) { 
    if((! defined $index) || ($index == $i)) { 
      for(my $j = 0; $j < $size; $j++) { 
        if($observed_AAR->[$i][$j] ne $expected_AAR->[$i][$j]) { 
          $pass_fail_char = "F";
        }
        if($observed_AAR->[$i][$j] && ($i < $j)) { 
          if($ret_str ne "") { 
            $ret_str .= ",";
          }
          if(defined $index) { 
            $ret_str .= sprintf("%s", $name_AR->[$j]); 
          }
          else { 
            $ret_str .= sprintf("%s%s%s", $name_AR->[$i], $div_char, $name_AR->[$j]); 
          }
          $nfound++;
        }
      }
    }
  }
  if($ret_str eq "") { $ret_str = "NONE"; }

  if(defined $index) { 
    $ret_str = $pass_fail_char . ":" . $ret_str;
  }
  else { 
    $ret_str = $pass_fail_char . ":" . $nfound . ":" . $ret_str;
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

# Subroutine: getOverlap()
# Args:       $start1: start position of hit 1 (must be <= $end1)
#             $end1:   end   position of hit 1 (must be >= $end1)
#             $start2: start position of hit 2 (must be <= $end2)
#             $end2:   end   position of hit 2 (must be >= $end2)
#
# Returns:    Number of nucleotides of overlap between hit1 and hit2,
#             0 if none
# Dies:       if $end1 < $start1 or $end2 < $start2.

sub getOverlap {
  my $sub_name = "getOverlap";
  my $nargs_exp = 4;
  if(scalar(@_) != 4) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($start1, $end1, $start2, $end2) = @_; 

  #printf("in $sub_name $start1..$end1 $start2..$end2\n");

  if($start1 > $end1) { die "ERROR in $sub_name start1 > end1 ($start1 > $end1)"; }
  if($start2 > $end2) { die "ERROR in $sub_name start2 > end2 ($start2 > $end2)"; }

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
  die "ERROR in $sub_name, unforeseen case in $start1..$end1 and $start2..$end2";

  return; # NOT REACHED
}

# Subroutine: getDistance()
# Args:       $start1: start position of hit 1 (must be <= $end1)
#             $end1:   end   position of hit 1 (must be >= $end1)
#             $start2: start position of hit 2 (must be <= $end2)
#             $end2:   end   position of hit 2 (must be >= $end2)
#
# Returns:    Distance between closest nucleotides in hit1 and hit2,
#             where hit1 = start1..end1 and hit2 = start2..end2.
#             If they overlap, return 0. 
#
# Dies:       if $end1 < $start1 or $end2 < $start2.

sub getDistance {
  my $sub_name = "getDistance";
  my $nargs_exp = 4;
  if(scalar(@_) != 4) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($start1, $end1, $start2, $end2) = @_; 

  #printf("in $sub_name $start1..$end1 $start2..$end2\n");

  if($start1 > $end1) { die "ERROR in $sub_name start1 > end1 ($start1 > $end1)"; }
  if($start2 > $end2) { die "ERROR in $sub_name start2 > end2 ($start2 > $end2)"; }

  # Given: $start1 <= $end1 and $start2 <= $end2.
  
  # Swap if nec so that $start1 <= $start2.
  if($start1 > $start2) { 
    my $tmp;
    $tmp   = $start1; $start1 = $start2; $start2 = $tmp;
    $tmp   =   $end1;   $end1 =   $end2;   $end2 = $tmp;
  }
  
  # 3 possible cases:
  # Case 1. $start1 <=   $end1 <  $start2 <=   $end2  distance is $start2-$end1
  # Case 2. $start1 <= $start2 <=   $end1 <    $end2  overlap: return 0
  # Case 3. $start1 <= $start2 <=   $end2 <=   $end1  overlap: return 0
  if($end1 < $start2) { return $start2-$end1; } # case 1
  if($end1 <   $end2) { return 0; }             # case 2
  if($end2 <=  $end1) { return 0; }             # case 3
  die "ERROR in $sub_name, unforeseen case in $start1..$end1 and $start2..$end2";

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
# Synopsis: Given a gap string that lists all gaps of a feature in an
#             alignment, determine if the gaps cause the predicted
#             length of the feature to be non-modulo 3. If so,
#             determine if there's exactly one gap that if we remove
#             it, the remaining gaps will make the predicted length of
#             the CDS to be modulo 3.
#
#             I think there can only be 0 or 1 such gap, but I'm not sure.
# 
#             Return two values: 
#             - the length of the speical gap if there is one, or 0 if there is not
#               (this value will be negative if it is a deletion relative to the 
#               reference, and positive if it is an insertion relative to the reference)
#             - the string that describes the special gap, 
#               or '-' if the predicted length of the feature is modulo 3 if we include all gaps
#               or '?' if the predicted length of the feature is not modulo 3 but there is no
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
# Returns:    Number of files actually created (can differ from requested
#             amount (which is $nfiles)).
# Dies:       if command fails or there is some other problem
#
sub splitFastaFile {
  my $sub_name = "splitFastaFile";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($esl_ssplit, $fafile, $nfiles) = @_;

  if(! -s $esl_ssplit) { die "ERROR the $esl_ssplit file does not exist or is empty"; }
  if(! -x $esl_ssplit) { die "ERROR the $esl_ssplit file is not executable"; }

  my $outfile = $fafile . ".esl-ssplit";
  my $cmd = "$esl_ssplit -v -n $fafile $nfiles > $outfile";
  runCommand($cmd, 0);
  # parse output to determine exactly how many files were created:
  open(IN, $outfile) || die "ERROR unable to open $outfile for reading";
  my $nfiles_created = 0;
  while(my $line = <IN>) { # we'll have one line per file created
    $nfiles_created++;
  }
  close(IN);
  unlink $outfile;

  return $nfiles_created;
}

# Subroutine: wrapperCombineExonsIntoCDS()
#
# Synopsis:   For all CDS, combine all exons into CDS. A wrapper function
#             for combineExonsIntoCDS().
#
# Args:       $nmdl:            total number of exons/models
#             $dir:             directory for output files
#             $key:             string for naming output files (e.g.: "predicted" or "corrected")
#             $mdl_AR:        ref to array of model names 
#             $mdl2mft_map_AR:  ref to array mapping models to CDS
#             $mdl_is_first_AR: ref to array signifying if a model is first in a CDS or not
#             $mdl_is_final_AR: ref to array signifying if a model is final in a CDS or not
#             $outfile_AR:      ref to array of output files, filled here
#
# Returns:    void
# Dies:       if something unexpected happens when reading the exon fasta files
#
sub wrapperCombineExonsIntoCDS {
  my $sub_name = "combineExonsIntoCDS";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($nmdl, $dir, $key, $mdl_AR, $mdl2mft_map_AR, $mdl_is_first_AR, $mdl_is_final_AR, $outfile_AR) = @_;

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));
  printf("%-65s ... ", "# Combining multi-exon $key CDS ");

  my @tmp_exon_fafile_A = (); # temporary array of exon fafiles for all exons in current CDS
  
  for(my $h = 0; $h < $nmdl; $h++) { 
    my $mft_i      = $mdl2mft_map_AR->[$h];
    my $cur_fafile = $dir . "/" . $mdl_AR->[$h] . ".fa";
    $cur_fafile    =~ s/ref/$key/;
    push(@tmp_exon_fafile_A, $cur_fafile);
    if($mdl_is_final_A[$h]) { 
      if(! $mdl_is_first_A[$h]) { 
        $cur_fafile =~ s/\.exon\.\d+//; # remove exon.<d> part of file name
        combineExonsIntoCDS(\@tmp_exon_fafile_A, $cur_fafile);
      }
      else { # a single exon gene, we should already have the sequence from alignHits
        if(! -s $cur_fafile) { die sprintf("ERROR, expected output fasta file for CDS %s does not exist: $cur_fafile", $mft_i+1); }
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
#                                    at position 1 of the feature sequence
#             $prot_must_stop_at_L:  'L' if translated protein sequence must stop
#                                    at final position of the feature sequence
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
  my $source_length = 0; # length of original feature sequence, including stop
  my $coords_length = 0; # length of original feature sequence that was translated, not included stop codon
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
        # that start at the beginning of the feature but do NOT end at the end of the 
        # feature, they may be shorter or longer
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
#             $source_name_R:        ref to name of source feature, filled here, can be undef
#             $source_coords_R:      ref to coordinates part of name of source feature, filled here, can be undef
#             $source_len_R:         ref to coordinate of length of source feature, filled here, can be undef
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

# Subroutine: matpeptFindAdjacentPeptides()
#
# Synopsis:   Find which peptides are adjacent to each other in the
#             reference and update %{$adj2mft_5p_H} and %{$adj2mft_3p_H}
#             with that information.
#
# Args:       $ref_mft_coords_AR:    reference to reference coordinates for each feature
#             $adj2mft_5p_HR:        ref to hash: key is index $i of feature in @{$ref_mft_coords_AR} and value
#                                    is $j: index of feature that it is immediately adjacent to on 5' end,
#                                    that is final nt of $j is immediately before first nt of $i
#                                    FILLED HERE
#             $adj2mft_3p_HR:        ref to hash: key is index $i of feature in @{$ref_mft_coords_AR} and value
#                                    is $j: index of feature that it is immediately adjacent to on 3' end,
#                                    that is first nt of $j is immediately after final nt of $i
#                                    FILLED HERE
#
#             If more than one feature is immediately adjacent to another, we store the one that is longest
#
# Returns:    void
# Dies:       if more than one feature is immediately adjacent to any other
#
sub matpeptFindAdjacentPeptides {
  my $sub_name = "matpeptFindAdjacentPeptides";
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

# Subroutine: matpeptFindPrimaryPeptides()
#
# Synopsis:   Find which peptides are 'primary'. A peptide is primary
#             if no other longer peptide overlaps with it.
#
# Args:       $ref_ol_AAR:         [0..$i..$nmdl-1][0..$j..$nmdl-1]: 1 if $i and $j overlap, else 0, PRE-FILLED
#             $ref_len_AR:         [0..$i..$nmdl-1]: length of predicted peptide $i (in nt), PRE-FILLED
#             $pept_is_primary_AR: [0..$i..$nmdl-1]: '1' if $i is a primary peptide, else 0, FILLED HERE
#
# Returns:    void, fills @{$mdl_is_primary_AR}
#
sub matpeptFindPrimaryPeptides {
  my $sub_name = "matpeptFindPrimaryPeptides";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ref_ol_AAR, $ref_len_AR, $pept_is_primary_AR) = @_;

  my $npept = scalar(@{$ref_len_AR});
  for(my $i = 0; $i < $npept; $i++) {
    my $is_primary = 1; # until proven otherwise
    my $ilen = $ref_len_AR->[$i];
    for(my $j = 0; $j < $npept; $j++) {
      if($ref_ol_AAR->[$i][$j] && $ref_len_AR->[$j] > $ilen) { 
        $is_primary = 0;
      }
    }
    $pept_is_primary_AR->[$i] = ($is_primary) ? 1 : 0;
    # if($pept_is_primary_AR->[$i]) { printf("pept $i is primary!\n"); }
  }

  return;
}

# Subroutine: matpeptFindPrimaryAdjacencies()
#
# Synopsis:   For each 'primary' peptide $i, find the other primary peptide that is adjacent before it
#             (and store it in $prv_adj_A[$i]) if any, and find the other primary peptide that is 
#             adjacent after it (and store it in $nxt_adj_A[$i]), if any. Here $i is before $j if $i < $j
#             and $i is after $j if $i > $j.
#
# Args:       $ref_adj_AAR:        [0..$i..$nmdl-1][0..$j..$nmdl-1]: 1 if $i and $j are adjacent, else 0, PRE-FILLED
#             $pept_is_primary_AR: [0..$i..$nmdl-1]: '1' if $i is a primary peptide, else 0, PRE-FILLED
#             $prv_adj_AR:         [0..$i..$nmdl-1]: $j if $j < $i and $i and $j are adjacent and $i and $j are primary peptides, FILLED HERE
#             $nxt_adj_AR:         [0..$i..$nmdl-1]: $j if $j > $i and $i and $j are adjacent and $i and $j are primary peptides, FILLED HERE
#
# Returns:    void, fills @{$prv_adj_AR} and @{$nxt_adj_AR}
# Dies:       if a primary peptide $i has more than 1 other primary adjacent peptides $j1 and $j2 such that $j1 < $i and $j2 < $i
#             if a primary peptide $i has more than 1 other primary adjacent peptides $j1 and $j2 such that $j1 > $i and $j2 > $i
#             both of these cases are expected to not happen
#
sub matpeptFindPrimaryAdjacencies {
  my $sub_name = "matpeptFindPrimaryAdjacencies";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ref_adj_AAR, $pept_is_primary_AR, $prv_adj_AR, $nxt_adj_AR) = @_;

  my $npept = scalar(@{$pept_is_primary_AR});
  my ($i, $j);
  for($i = 0; $i < $npept; $i++) {
    $prv_adj_AR->[$i] = -1;
    $nxt_adj_AR->[$i] = -1;
    if($pept_is_primary_AR->[$i]) { 
      # look for previous peptide
      for($j = 0; $j < $i; $j++) { 
        if($pept_is_primary_AR->[$j] && $ref_adj_AAR->[$i][$j]) { 
          if($prv_adj_AR->[$i] != -1) { 
            die sprintf("ERROR in $sub_name, found two primary peptides %d and %d which are prior to and adjacent to primary peptide %d\n", 
                        $prv_adj_AR->[$i], $j, $i);
          }
          $prv_adj_AR->[$i] = $j;
          # printf("set prv_adj_AR->[$i] to $j\n");
        }
      }
      # look for next peptide
      for($j = $i+1; $j < $npept; $j++) { 
        if($pept_is_primary_AR->[$j] && $ref_adj_AAR->[$i][$j]) { 
          if($nxt_adj_AR->[$i] != -1) { 
            die sprintf("ERROR in $sub_name, found two primary peptides %d and %d which are after than and adjacent to primary peptide %d\n", 
                        $nxt_adj_AR->[$i], $j, $i);
          }
          $nxt_adj_AR->[$i] = $j;
          # printf("set nxt_adj_AR->[$i] to $j\n");
        }
      }
    }
  }

  return;
}


######################
# OUTPUT subroutines #
######################
#
# Subroutine: getHeadings()
#
# Synopsis:   For 'sequences are rows' tabular output, output the headings.
#
# IMPORTANT:  This function must stay in sync with the long block of code
#             in the main script entitled 'Pass through all accessions, and 
#             gather and output annotation for each'. Here we define the
#             headers of the output, in the main script we add output for
#             each of those headers, so they must stay in sync.
# 
# Args:       $do_seqrow:          '1' if we're outputting in sequence-as-rows mode
#             $do_seqcol:          '1' if we're outputting in sequence-as-columns mode
#             $do_matpept:         '1' if we're in matpept mode
#             $do_nofid:           '1' if we're not printing fractional ids, else '0'
#             $do_mdlb:            '1' if we're not printing model boundaries, else '0'
#             $do_noss3:           '1' if we're not printing SS3 columns, else '0'
#             $do_nostop:          '1' if we're not printing Stop codons, else '0'
#             $do_fullolap         '1' if we're printing full overlap string
#             $do_fulladj          '1' if we're printing full adjacency string
#             $origin_seq:         origin sequence, or undef if ! defined
#             $ref_tot_nexons:     number of total exons in reference
#             $nmdl:               number of total models
#             $mdl2mft_map_AR:     ref to @mdl2mft_map array, already filled
#             $mft2exon_map_AR:    ref to @mft2exon_map array, already filled
#             $mdl_is_final_AR:    ref to @mdl_is_final_A, already filled
#             $mdl_is_first_AR:    ref to @mdl_is_first_A, already filled
#             $mft_out_short_AR:   ref to @mft_out_short_A, already filled
#             $mft_out_product_AR: ref to @mft_out_product_A, already filled
#             $out_col_header_AAR: ref to 2D array of column headers, filled here
#                                  undef unless $do_seqrow is '1'
#             $out_row_header_AR:  ref to 1D array of row headers, filled here
#                                  undef unless $do_seqcol is '1'
#             $out_header_exp_AR:  ref to 1D array of header explanations, each
#                                  element is a line to be printed in explanatory
#                                  section of the output; filled here
sub getHeadings {
  my $sub_name = "getHeadings";
  my $nargs_exp = 21;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($do_seqrow, $do_seqcol, $do_matpept, $do_nofid, $do_mdlb, $do_noss3, $do_nostop, $do_fullolap, $do_fulladj, $origin_seq, $ref_tot_nexons, $nmdl, $mdl2mft_map_AR, $mft2exon_map_AR, $mdl_is_first_AR, $mdl_is_final_AR, $mft_out_short_AR, $mft_out_product_AR, $out_col_header_AAR, $out_row_header_AR, $out_header_exp_AR) = @_;

  # contract checks
  if($do_seqrow     &&    $do_seqcol)  { die "ERROR in $sub_name, both $do_seqrow and $do_seqcol are '1'"; }
  if((! $do_seqrow) && (! $do_seqcol)) { die "ERROR in $sub_name, neither $do_seqrow nor $do_seqcol are '1'"; }
  if($do_seqrow && (! defined $out_col_header_AAR)) { die "ERROR in $sub_name, $do_seqrow is '1' but out_col_header_AAR is undefined"; }
  if($do_seqrow && defined $out_row_header_AR)      { die "ERROR in $sub_name, $do_seqrow is '1' but out_row_header_AR is defined"; }
  if($do_seqcol && (! defined $out_row_header_AR))  { die "ERROR in $sub_name, $do_seqcol is '1' but out_row_header_AR is undefined"; }
  if($do_seqcol && defined $out_col_header_AAR)     { die "ERROR in $sub_name, $do_seqcol is '1' but out_col_header_AAR is defined"; }

  my $width;  # width of a field
  my $pad;    # string of all spaces used for pretty formatting
  my $width_result = 5 + $ref_tot_nexons + 2;
  my $row_div_char = ":";

## If $do_seqrow: 
## We store the column headers in a 2D array @{$out_col_header_AAR}
## The first dimension is of size 5, one array for each 'line' of
## output in the sequence-per-row output format. The 2nd dimension
## arrays are different sizes because at higher levels single
## columns can pertain to multiple columns at lower levels, for
## example, "origin sequence" at level 2 pertains to all of "#", 
## "start", "stop", "offst", and "PF" at level 4. When we output
## the headings we don't have to know the relationship between
## the levels, we simply format them to the proper width and
## print them.
##
## Example of sequence-per-row output:
#                                                                      CDS #1 [single exon; +]                          CDS #2 [2 exons; -]                                        
#                                          origin sequence                 movement protein                        replication-associated protein                                         GenBank annotation
#                                   ----------------------  ---------------------------------------------  --------------------------------------------------------------------------    -------------------                   
# idx accession           totlen  # start  stop offst PF    start1    stop1  fid1 md1 length SS3 stp PF    start1    stop1  fid1 md1   start2    stop2  fid2 md2 length SS3 stp PF     cds   exons  match             overlaps?  result      
#---- -------------------  ------ -- ----- ----- ----- --  -------- -------- ----- --- ------ --- --- --  -------- -------- ----- ---   ------    ------ ---- --- ------ --- --- --     ----- -----  -----   -------------------  ------------
#
#
## If $do_seqrow: 
##  We store the row headers in a 1D array @{$out_row_header_AR}.
##  We have the same values as in the column headers, but each
##  level is concatenated together per row. Here's the row header
##  information that pertains to the example column header example
##  above.
#
#  idx
#  accession
#  totlen
#  origin sequence:#
#  origin sequence:start
#  origin sequence:stop
#  origin sequence:offst
#  origin sequence:PF
#  CDS #1 [single exon; +]:movement protein:start1
#  CDS #1 [single exon; +]:movement protein:stop1
#  CDS #1 [single exon; +]:movement protein:fid1
#  CDS #1 [single exon; +]:movement protein:md1
#  CDS #1 [single exon; +]:movement protein:length
#  CDS #1 [single exon; +]:movement protein:SS3
#  CDS #1 [single exon; +]:movement protein:stp
#  CDS #1 [single exon; +]:movement protein:PF
#  CDS #2 [2 exons; -]:replication associated protein:start1
#  CDS #2 [2 exons; -]:replication associated protein:stop1
#  CDS #2 [2 exons; -]:replication associated protein:fid1
#  CDS #2 [2 exons; -]:replication associated protein:md1
#  CDS #2 [2 exons; -]:replication associated protein:start2
#  CDS #2 [2 exons; -]:replication associated protein:stop2
#  CDS #2 [2 exons; -]:replication associated protein:fid2
#  CDS #2 [2 exons; -]:replication associated protein:md2
#  CDS #2 [2 exons; -]:replication associated protein:length
#  CDS #2 [2 exons; -]:replication associated protein:SS3
#  CDS #2 [2 exons; -]:replication associated protein:stp
#  CDS #2 [2 exons; -]:replication associated protein:PF
#  GenBank annotation:cds
#  GenBank annotation:exons
#  GenBank annotation:match
#  overlaps?
#  result

  my $tok1; # first  level token (line 1 of column headers) 
  my $tok2; # second level token (line 2 of column headers) 
  my $tok3; # third  level token (line 3 of column headers) 
  my $tok4; # fourth level token (line 4 of column headers) 
  my $tok5; # fifth  level token (line 5 of column headers) 

  my $exp_tok1; # first  level explanation token, only used if we want this to be different from $tok1
  my $exp_tok4; # fourth level explanation token, only used if we want this to be different from $tok4

  # first, initialize the @{$out_header_exp_AR} with the first two lines:
  push(@{$out_header_exp_AR}, "#\n");
  if($do_seqrow) { 
    push(@{$out_header_exp_AR}, "# Explanations of column headings (in left to right order):\n");
  }
  elsif($do_seqcol) { 
    push(@{$out_header_exp_AR}, "# Explanations of row headings on each page:\n");
  }
  push(@{$out_header_exp_AR}, "#\n");

  # column/row #2: 'idx'
  $tok1 = sprintf("%-4s  ", "");
  $tok2 = $tok1;
  $tok3 = $tok1;
  $tok4 = sprintf("%-4s  ", " idx");
  $tok5 = sprintf("%-4s  ", "----");
  if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
  elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); }
  getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "index of genome in list");

  # column/row #2: 'accession'
  $tok1 = sprintf("%-19s  ", "");
  $tok2 = $tok1;
  $tok3 = $tok1;
  $tok4 = sprintf("%-19s  ", " accession");
  $tok5 = sprintf("%-19s  ", "-------------------");
  if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
  elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); }
  getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "GenBank accession for genomic sequence");

  # column/row #3: 'totlen'
  $tok1 = sprintf("%-6s", "");
  $tok2 = $tok1;
  $tok3 = $tok1;
  $tok4 = sprintf("%-6s", "totlen");
  $tok5 = sprintf("%-6s", "------");
  if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
  elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); }
  getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "total length (nt) for accession");
  getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef); # adds a blank line

  if(defined $origin_seq) { 
    # column/row #4: 'origin sequence:#'
    $tok1 = sprintf("  %22s", "");
    $tok2 = sprintf("  %22s", "   origin sequence");
    $tok3 = sprintf("  %22s", "----------------------");
    $tok4 = sprintf(" %2s", " #");
    $tok5 = sprintf(" %2s", "--");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "number of occurences of origin sequence (input with -oseq) in genome");

    # column/row #5: 'origin sequence:start'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %5s", "start");
    $tok5 = sprintf(" %5s", "-----");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "start position of lone occurence of origin sequence (if only 1 exists)");

    # column/row #6: 'origin sequence:stop'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %5s", "stop");
    $tok5 = sprintf(" %5s", "-----");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "stop  position of lone occurence of origin sequence (if only 1 exists)");

    # column/row #7: 'origin sequence:offst'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %5s", "offst");
    $tok5 = sprintf(" %5s", "-----");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "predicted offset of genome, number of nucleotides to shift start (>0: clockwise; <0: counterclockwise)");

    # column/row #7: 'origin sequence:PF'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %2s", "PF");
    $tok5 = sprintf(" %2s", "--");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "'P' (for PASS) if there is exactly 1 occurence of the offset, else 'F' for FAIL");
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef); # adds a blank line
  } # end of 'if(defined $orig_seq)'

  # create columns for 5'UTR, if $do_matpept:
  if($do_matpept) { 
    $width = 6 + 1 + 6 + 1 + 6; #20
    $tok1 = sprintf("  %*s", $width, "");
    $tok2 = sprintf("         %*s", $width, "5' UTR");
    $tok3 = sprintf("  %*s", $width, monocharacterString($width, "-"));
    $tok4 = sprintf("  %6s", "start");
    $tok5 = sprintf("  ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "start position of 5' UTR (inferred from other predictions)");

    $tok4 = sprintf(" %6s", "stop");
    $tok5 = sprintf(" ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "stop  position of 5' UTR (inferred from other predictions)");

    $tok4 = sprintf(" %6s", "length");
    $tok5 = sprintf(" ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "length of 5' UTR (inferred from other predictions)");
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef); # adds a blank line
  }

  # create columns for main features (CDS/mat_peptide)
  $width = 0;
  my $do_explanation = 1;
  for(my $h = 0; $h < $nmdl; $h++) { 
    $width += 18;
    my $mft_i = $mdl2mft_map_AR->[$h];
    if(! $do_nofid)  { $width += 6;  }
    if(! $do_nomdlb) { $width += 4;  }
    if(! $do_noolap) { $width += 11; }
    if($do_matpept)  { $width += 11; }
    if($mdl_is_final_AR->[$h]) { 
      $width += 9;
      if(! $do_noss3)  { $width += 4; }
      if(! $do_nostop) { $width += 4; }
      $tok1     = sprintf("  %*s", $width, $mft_out_short_AR->[$mft_i] . monocharacterString(($width-length($mft_out_short_AR->[$mft_i]))/2, " "));
      $exp_tok1 = ($do_matpept) ? "MP #<i>" : "CDS #<i>";
      $tok2 = sprintf("  %*s", $width, substr($mft_out_product_AR->[$mft_i], 0, $width) . monocharacterString(($width-length($mft_out_product_AR->[$mft_i]))/2, " "));
      $tok3 = sprintf("  %s", monocharacterString($width, "-"));
      $tok4 = sprintf("  %8s", sprintf("%s%s", "start", $mft2exon_map_AR->[$h]+1));
      $exp_tok4 = "start<j>";
      $tok5 = sprintf("  %8s", "--------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      $width = 0; # reset width, this is impt
    }
    else { # not the final exon, still need start coordinate
      $width += 1;
      $tok1 = sprintf("    %s", $mft_out_short_AR->[$mft_i]);   # used only for getHeadingsSeqColHelper
      $tok2 = sprintf("    %s", $mft_out_product_AR->[$mft_i]); # used only for getHeadingsSeqColHelper
      $tok4 = sprintf("  %8s", sprintf("%s%s", "start", $mft2exon_map_AR->[$h]+1));
      $exp_tok4 = "start<j>";
      $tok5 = sprintf("  %8s", "--------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
    }
    my $exp_substr = ($do_matpept) ? "coding sequence part <j> of mat_peptide" : "exon #<j> of CDS #<i>";
    if($do_explanation) { 
      getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "start position of $exp_substr");
      getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "enclosed in brackets \"\[e\]\" if different from all exon starts in existing GenBank annotation");
    }
    
    # stop, fid, and md rows take place for all exons
    # only token 4 changes
    $tok4 = sprintf(" %8s", sprintf("%s%s", "stop", $mft2exon_map_AR->[$h]+1));
    $exp_tok4 = "stop<j>";
    $tok5 = sprintf(" %8s", "--------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
    if($do_explanation) { 
      getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "stop  position of $exp_substr");
      getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "enclosed in brackets \"\[e\]\" if different from all exon stops in existing GenBank annotation");
    }
    
    if(! $do_nofid) { 
      $tok4 = sprintf(" %5s", sprintf("%s%s", "fid", $mft2exon_map_AR->[$h]+1));
      $exp_tok4 = "fid<j>";
      $tok5 = sprintf(" %5s", "-----");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "fractional identity between $exp_substr and reference genome");
      }
    }

    $exp_substr = $do_matpept ? "mat_peptide coding sequence" : "exon coding sequence";
    if(! $do_nomdlb) { 
      $tok4 = sprintf(" %3s", sprintf("%s%s", "md", $mft2exon_map_AR->[$h]+1));
      $exp_tok4 = "md<j>";
      $tok5 = sprintf(" %3s", "---");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, sprintf("annotation indicating if alignment to reference extends to 5' and 3' end of reference $exp_substr."));
      }
    }
    
    if(! $do_noolap) { 
      $tok4 = sprintf(" %10s", sprintf("%s%s", "overlaps", $mft2exon_map_AR->[$h]+1));
      $exp_tok4 = "overlaps<j>";
      $tok5 = sprintf(" %10s", "----------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, sprintf("'P' or 'F' followed by list of %s this %s overlaps with", ($do_matpept) ? "mat_peptide" : "exon", ($do_matpept) ? "mat_peptide" : "exon"));
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "first letter is 'P' if agrees exactly with reference, else 'F'"); # adds a second line to explanation
      }
}

    if($do_matpept) { 
      $tok4 = sprintf(" %10s", sprintf("%s%s", "adjcnces", $mft2exon_map_AR->[$h]+1));
      $exp_tok4 = "adjcnces<j>";
      $tok5 = sprintf(" %10s", "----------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, sprintf("'P' or 'F' followed by list of %s this %s is adjacent with", ($do_matpept) ? "mat_peptide" : "exon", ($do_matpept) ? "mat_peptide" : "exon"));
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "first letter is 'P' if agrees exactly with reference, else 'F'"); # adds a second line to explanation
      }
    }

    $exp_substr = $do_matpept ? "mat_peptide coding sequence" : "CDS";
    if($mdl_is_final_AR->[$h]) { 
      $tok4 = sprintf(" %6s", "length");
      $tok5 = sprintf(" %6s", "------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $tok4, undef, "length of $exp_substr #<i> (all exons summed)");
      }      

      if(! $do_noss3) { 
        $tok4 = sprintf(" %3s", "ss3");
        $tok5 = sprintf(" %3s", "---");
        if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
        elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      }
      if($do_explanation) { 
        if($do_matpept) { 
          #TODO FIX ME! This is inaccurate for mat_peptides
          getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if predicted mat_peptide has a valid start codon, stop codon and is a multiple of 3");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "first  character: '.' if predicted mat_peptide has a valid start codon, else '!'");          
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "second character: '.' if predicted mat_peptide has a valid stop  codon, else '!'");      
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "third  character: '.' if predicted mat_peptide has a length which is a multiple of three, else '!'");
        }      
        else {
          getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if predicted CDS has a valid start codon, stop codon and is a multiple of 3");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "first  character: '.' if predicted CDS has a valid start codon, else '!'");          
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "second character: '.' if predicted CDS has a valid stop  codon, else '!'");      
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "third  character: '.' if predicted CDS has a length which is a multiple of three, else '!'");
        }
      }
      if(! $do_nostop) { 
        $tok4 = sprintf(" %3s", "stp");
        $tok5 = sprintf(" %3s", "---");
        if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
        elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
        if($do_explanation) { 
          getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, ($do_matpept) ? "the final codon for this mat_peptide sequence" : "the predicted stop codon for this CDS");
        }
      }
      
      $tok4 = sprintf(" %2s", "PF");
      $tok5 = sprintf(" %2s", "--");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      if($do_explanation) { 
        if($do_matpept) { 
          getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if this mat_peptide PASSED ('P') or FAILED ('F')");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  a mat_peptide coding sequence PASSES ('P') if and only if");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  conditions are met (else it FAILS):");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (1) it has a valid start codon or homologous reference mat_peptide");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "      does not");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (2) it has a valid stop  codon immediately after its predicted");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "      end or reference mat_peptide does not");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (3) its length is a multiple of 3");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (4) has a pairwise alignment to the homologous reference met_peptide that");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "      extends to the 5' and 3' boundary of the reference annotation");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (5) overlaps with exact same set of other mat_peptides as the homologous");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "      reference mat_peptide");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (6) is adjacent to the exact same set of other mat_peptides as the");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "      homologous reference mat_peptide");

        }
        else { 
          getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if this CDS PASSED ('P') or FAILED ('F')");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  a CDS sequence PASSES ('P') if and only if all of the following");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  conditions are met (else it FAILS):");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (1) it has a valid start codon at beginning of its first exon");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (2) it has a valid stop  codon at end of its final exon");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (3) its length is a multiple of 3");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (4) all of its exons have a pairwise alignment to the homologous");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "      reference exon that extends to the 5' and 3' boundary of the");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "      reference annotation.");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (5) all of its exons overlap with exact same set of other exons as the");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "      homologous reference CDS");
        }
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
      }
    }
    $do_explanation = 0; # once we see the final exon of the first CDS, we don't need to print CDS explanations anymore
  } # end of 'if($mdl_is_final_AR->[$h])'

  # create columns for 3'UTR, if $do_matpept:
  if($do_matpept) { 
    $width = 6 + 1 + 6 + 1 + 6; #20
    $tok1 = sprintf("  %*s", $width, "");
    $tok2 = sprintf("         %*s", $width, "3' UTR");
    $tok3 = sprintf("  %*s", $width, monocharacterString($width, "-"));
    $tok4 = sprintf("  %6s", "start");
    $tok5 = sprintf("  ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "start position of 3' UTR (inferred from other predictions)");

    $tok4 = sprintf(" %6s", "stop");
    $tok5 = sprintf(" ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "stop  position of 3' UTR (inferred from other predictions)");

    $tok4 = sprintf(" %6s", "length");
    $tok5 = sprintf(" ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "length of 3' UTR (inferred from other predictions)");

    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
  }

  # "totlen"
  $tok1 = sprintf("  %6s", "");
  $tok2 = sprintf("  %6s", "");
  $tok3 = sprintf("  %6s", "");
  $tok4 = sprintf("  %6s", "totlen");
  $tok5 = sprintf("  %6s", "------");
  if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
  elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); }
  getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "total length (nt) for accession (repeated for convenience)");

  # "avgid"
  $tok1 = sprintf("  %5s", "");
  $tok2 = sprintf("  %5s", "");
  $tok3 = sprintf("  %5s", "");
  $tok4 = sprintf("  %5s", "avgid");
  $tok5 = sprintf("  %5s", "-----");
  if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
  elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); }
  getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "average fractional identity of all pairwise alignments for this accession");
  getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);

  # existing GenBank annotation
  if(! $do_noexist) { 
    $tok1 = sprintf("  %19s", "");
    $tok2 = sprintf("  %19s", "GenBank annotation");
    $tok3 = sprintf("  %19s", "-------------------");
    $tok4 = sprintf("  %5s", ($do_matpept) ? "mp" : "cds");
    $tok5 = sprintf("  %5s", "-----");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, sprintf("number of %s in the existing GenBank annotation for this accession", ($do_matpept) ? "mat_peptides" : "CDS"));
    
    $tok4 = sprintf("  %5s", "exons");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "total number of exons in the existing GenBank annotation for this accession");
    
    $tok4 = sprintf("  %5s", "match");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "number of exons in the existing GenBank annotation for which existing and predicted annotation agree exactly");
  }
  getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);

  # overlaps?
  if($do_fullolap) { 
    $tok1 = sprintf("  %20s", "");
    $tok2 = sprintf("  %20s", "");
    $tok3 = sprintf("  %20s", "");
    $tok4 = sprintf("  %20s", " overlaps?");
    $tok5 = sprintf("  %20s", "--------------------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("text describing which (if any) of the predicted %s overlap with each other", $do_matpept ? "mat_peptides" : "exons"));
    getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "first character:   'P' for PASS if predicted annotation for this accession has same overlaps as the reference");
    getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "                   'F' for FAIL if it does not");
    getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("second character:  number of overlaps between any two %s", $do_matpept ? "mat_peptides" : "exons"));
    getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("remainder of line: text explaining which exons overlap", $do_matpept ? "mat_peptides" : "exons"));
    if($do_matpept) { 
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("  e.g.: \"3.1|4.1\" indicates exon #1 of mat_peptide #3 overlaps with exon #1 of mat_peptide #4 on either strand"));
    }
    else { 
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("  e.g.: \"3.2/4.1\" indicates exon #2 of CDS #3 overlaps with exon #1 of CDS #4 on either strand"));
    }
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
  }

  # adjacencies?
  if($do_fulladj) { 
      $tok1 = sprintf("  %122s", "");
      $tok2 = sprintf("  %122s", "");
      $tok3 = sprintf("  %122s", "");
      $tok4 = sprintf("  %122s", " adjacencies?");
      $tok5 = sprintf("  %122s", "--------------------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); }
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("text describing which (if any) of the predicted %s are adjacent to each other", $do_matpept ? "mat_peptides" : "exons"));
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("two %s i and j are adjacent if i < j and final nt of i is 1 less than first nt of j", $do_matpept ? "mat_peptides" : "exons"));
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "first character:   'P' for PASS if predicted annotation for this accession has same adjacencies as the reference");
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "                   'F' for FAIL if it does not");
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("second character:  number of overlaps between any two %s", $do_matpept ? "mat_peptides" : "exons"));
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("remainder of line: text explaining which exons overlap", $do_matpept ? "mat_peptides" : "exons"));
      if($do_matpept) { 
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("  e.g.: \"3.1|4.1\" indicates exon #1 of mat_peptide #3 overlaps with exon #1 of mat_peptide #4 on either strand"));
    }
    else { 
      getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, sprintf("  e.g.: \"3.2/4.1\" indicates exon #2 of CDS #3 overlaps with exon #1 of CDS #4 on either strand"));
    }
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
  }

  # result
  $tok1 = sprintf("  %*s",  $width_result, "");
  $tok2 = sprintf("  %*s",  $width_result, "");
  $tok3 = sprintf("  %*s",  $width_result, "");
  $tok4 = sprintf("  %-*s", $width_result, "result");
  $tok5 = sprintf("  %-*s", $width_result, monocharacterString($width_result, "-"));
  if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
  elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); }

  getHeadingsExplanationHelper($out_header_exp_AR, $tok4, undef, undef, "\"PASS\" or \"FAIL\". \"PASS\" if and only if all tests for this accession PASSED ('P')");
  getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "as indicated in the \"PF\" columns. Followed by the individual P/F results in order.");
  getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "ADD TEXT ABOUT OVERLAP?");
  return;
}

#
# Subroutine: getHeadingsSeqRowHelper()
#
# Synopsis:   Helper function for getHeadings() when used in sequences-as-rows mode.
#             Given up to 5 tokens, add them to the appropriate place in @{$out_col_header_AAR}.
#
# Args:       $out_col_header_AAR:  ref to output column header 2D array
#             $tok1:                token 1, can be undef
#             $tok2:                token 2, can be undef
#             $tok3:                token 3, can be undef
#             $tok4:                token 4, can be undef
#             $tok5:                token 5, can be undef
#
sub getHeadingsSeqRowHelper { 
  my $sub_name = "getHeadingsSeqRowHelper";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_col_header_AAR, $tok1, $tok2, $tok3, $tok4, $tok5) = @_;

  if(defined $tok1) { push(@{$out_col_header_AAR->[0]}, $tok1); }
  if(defined $tok2) { push(@{$out_col_header_AAR->[1]}, $tok2); }
  if(defined $tok3) { push(@{$out_col_header_AAR->[2]}, $tok3); }
  if(defined $tok4) { push(@{$out_col_header_AAR->[3]}, $tok4); }
  if(defined $tok5) { push(@{$out_col_header_AAR->[4]}, $tok5); }
  
  return;
}

#
# Subroutine: getHeadingsSeqColHelper()
#
# Synopsis:   Helper function for getHeadings() when used in sequences-as-columns mode.
#             Given up to 5 tokens, add them to the appropriate place in @{$out_row_header_AR}.
#
# Args:       $out_row_header_AR:  ref to output column header 2D array
#             $div_char:           divider character to put between tokens
#             $tok1:               token 1, can be undef
#             $tok2:               token 2, can be undef
#             $tok3:               token 3, can be undef
#             
sub getHeadingsSeqColHelper { 
  my $sub_name = "getHeadingsSeqColHelper";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_row_header_AR, $div_char, $tok1, $tok2, $tok3) = @_;

  # remove whitespace at beginning and end of tokens
  if(defined $tok1) { $tok1 =~ s/^\s+//; $tok1 =~ s/\s+$//; }
  if(defined $tok2) { $tok2 =~ s/^\s+//; $tok2 =~ s/\s+$//; }
  if(defined $tok3) { $tok3 =~ s/^\s+//; $tok3 =~ s/\s+$//; }

  my $toadd = $tok1;
  if(defined $tok2) { $toadd .= $div_char . $tok2; }
  if(defined $tok3) { $toadd .= $div_char . $tok3; }

  push(@{$out_row_header_AR}, $toadd); 

  return;
}

#
# Subroutine: getHeadingsExplanationHelper()
#
# Synopsis:   Helper function for getHeadings() for adding explanatory text to the 
#             @{$out_header_exp_AR} array.
#             Given up to 3 tokens that define the header, and one that is the 
#             explanatory text. Can be used in 3 modes:
#
#             Mode 1: at least one of $tok1, $tok2, $tok3 is defined
#                     and $desc is defined 
#                     In this mode, determine header by concatenating all of
#                     $tok1, $tok2, and $tok3 that are defined and use
#                     $desc as the description.
#
#             Mode 2: none of $tok1, $tok2, $tok3 is defined and $desc is
#                     defined.
#                     In this mode, header is blank, and use $desc as 
#                     the description.
#
#             Mode 3: none of $tok1, $tok2, $tok3 is defined and $desc is
#                     not defined either
#                     In this mode, add a blank line to @{$out_row_header_AR}.
#
# Args:       $out_header_exp_AR:  ref to output column header 2D array
#             $tok1:               token 1, can be undef
#             $tok2:               token 2, can be undef
#             $tok3:               token 3, can be undef
#             $desc:               description text, can be undef
#             
sub getHeadingsExplanationHelper { 
  my $sub_name = "getHeadingsExplanationHelper";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_header_exp_AR, $tok1, $tok2, $tok3, $desc) = @_;

  my $width = 35;
  if(defined $tok1) { $tok1 =~ s/^\s+//; $tok1 =~ s/\s+$//; }
  if(defined $tok2) { $tok2 =~ s/^\s+//; $tok2 =~ s/\s+$//; }
  if(defined $tok3) { $tok3 =~ s/^\s+//; $tok3 =~ s/\s+$//; }

  my $header = "";
  if(defined $tok1) { 
    $header .= $tok1;
  }
  if(defined $tok2) { 
    if($header ne "") { $header .= ":"; }
    $header .= $tok2;
  }
  if(defined $tok3) { 
    if($header ne "") { $header .= ":"; }
    $header .= $tok3;
  }
  if($header ne "") { 
    $header = "\"" . $header . "\":";
  }

  if(defined $desc) { 
    push(@{$out_header_exp_AR}, sprintf("# %-*s %s\n", $width, $header, $desc)); 
  }
  else { 
    if($header ne "") { 
      die "ERROR in $sub_name, desc is not defined but one of the header tokens is"; 
    }
    push(@{$out_header_exp_AR}, "#\n");
  }

  return;
}

#
# Subroutine: getActualAnnotations()
#
# Synopsis:   Store actual annotations for an accession in preparation for outputting
#             a summary for this accession.
#
# Args:       $accn:                accession
#             $mft_tbl_HHAR:        ref to mft_tbl_HHA
#             $act_exon_starts_AAR: FILLED HERE [0..$nmft-1][0..$nexons-1] start positions of actual annotations of exons for this accn, $nexons is CDS specific
#             $act_exon_stops_AAR:  FILLED HERE [0..$nmft-1][0..$nexons-1] stop  positions of actual annotations of exons for this accn, $nexons is CDS specific
#             $tot_nexons_R:        FILLED HERE, total number of annotated exons for this accession
#             $nmft_R:              FILLED HERE, number of CDS for this accession
#
sub getActualAnnotations { 
  my $sub_name = "getActualAnnotations";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($accn, $mft_tbl_HHAR, $act_exon_starts_AAR, $act_exon_stops_AAR, $tot_nexons_R, $nmft_R) = @_;

  my ($nmft, $npos, $nneg, $nunc, $nbth, $strand_str, $tot_nexons); 

  if(exists ($mft_tbl_HHAR->{$accn})) { 
    ($nmft, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats($mft_tbl_HHAR, $accn);
    my @mft_len_A = ();
    my @mft_coords_A = ();
    my @mft_product_A = ();
    getLengthStatsAndCoordStrings($mft_tbl_HHAR, $accn, \@mft_len_A, \@mft_coords_A);
    getQualifierValues($mft_tbl_HHAR, $accn, "product", \@mft_product_A);
    for(my $i = 0; $i < $nmft; $i++) { 
      # determine start and stop positions of all exons
      my @starts_A = ();
      my @stops_A  = ();
      my $nexons   = 0;
      @{$act_exon_starts_AAR->[$i]} = ();
      @{$act_exon_stops_AAR->[$i]}  = ();
      startStopsFromCoords($mft_coords_A[$i], \@starts_A, \@stops_A, \$nexons);
      
      my $strand = substr($strand_str, $i, 1);
      if($strand eq "-") { # switch order of starts and stops, because 1st exon is really last and vice versa
        @starts_A = reverse @starts_A;             # exons will be in reverse order, b/c we're on the negative strand
        @stops_A  = reverse @stops_A;              # exons will be in reverse order, b/c we're on the negative strand
        @{$act_exon_starts_AAR->[$i]} = @stops_A;  # save stops  to starts array b/c we're on the negative strand
        @{$act_exon_stops_AAR->[$i]}  = @starts_A; # save starts to stops  array b/c we're on the negative strand
      }
      else { 
        @{$act_exon_starts_AAR->[$i]} = @starts_A;
        @{$act_exon_stops_AAR->[$i]}  = @stops_A;
      }
      $tot_nexons += $nexons;
    }

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
    $nmft       = 0;
    $tot_nexons = 0;
  }

  $$nmft_R       = $nmft;
  $$tot_nexons_R = $tot_nexons;

  return;
}

#
# Subroutine: getOseqOutput()
#
# Synopsis:   Determine output related to the -oseq option, describing the existence and position
#             of the origin sequence in the genome.
#
# Args:       $accn:                accession we're working on
#             $origin_seq:          the origin sequence, or undef if -oseq not used
#             $origin_offset:       offset of origin in reference genome, or undef if -oseq not used
#             $origin_coords_HAR:   ref to origin_coords_HA hash of arrays
#
# Returns: 5 values:
#          $oseq_ct:       number of occurences of origin sequence found
#          $oseq_start:    start position of origin seq if $oseq_ct == 1, else '-'
#          $oseq_stop:     stop  position of origin seq if $oseq_ct == 1, else '-'
#          $oseq_offset:   offset position of origin seq if $oseq_ct == 1, else '-'
#          $oseq_passfail: 'P' if $oseq_ct is 1, else 'F'
# 
sub getOseqOutput {
  my $sub_name = "getOseqOutput";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($accn, $origin_seq, $origin_offset, $origin_coords_HAR) = @_;

  if(! defined $origin_seq)    { die "ERROR in $sub_name but origin sequence is undefined" }
  if(! defined $origin_offset) { die "ERROR in $sub_name but origin offset is undefined" }

  my $oseq_ct = (exists $origin_coords_HAR->{$accn}) ? scalar(@{$origin_coords_HAR->{$accn}}) : 0;;
  my $oseq_start    = "-"; # changed below if $oseq_ct == 1
  my $oseq_stop     = "-"; # changed below if $oseq_ct == 1
  my $oseq_offset   = "-"; # changed below if $oseq_ct == 1
  my $oseq_passfail = "F"; # changed below if $oseq_ct == 1

  if($oseq_ct == 1) { 
    ($oseq_start, $oseq_stop) = split(":", $origin_coords_HAR->{$accn}[0]);
    # print("accn: $accn oseq_start: $oseq_start oseq_stop: $oseq_stop\n");
    $oseq_offset = ($oseq_start < 0) ? ($oseq_start + $origin_offset) : ($oseq_start + $origin_offset - 1);
    # $oseq_offset is now number of nts to shift origin in counterclockwise direction
    if($oseq_offset > ($totlen_H{$accn} / 2)) { # simpler (shorter distance) to move origin clockwise
      $oseq_offset = $totlen_H{$accn} - $oseq_offset; # note, we don't add 1 here
    }
    else { # simpler to shift origin in counterclockwise direction, we denote this as a negative offset
      $oseq_offset *= -1;
    }
    $oseq_passfail = "P";
  }

  return ($oseq_ct, $oseq_start, $oseq_stop, $oseq_offset, $oseq_passfail);
}
# Subroutine: outputColumnHeaderExplanations()
# Args:       $out_header_exp_AR: ref to array of output explanation lines
#
# Synopsis:   Prints out output lines in @{$out_header_exp_AR} and exits.
#
# Returns:    void

sub outputColumnHeaderExplanations {
  my $sub_name = "outputColumnHeaderExplanations";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_header_exp_AR) = @_;

  foreach my $line (@{$out_header_exp_AR}) { 
    print $line;
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
# Synopsis:   Output gap information for all predicted features in a single sequence 
#             to a file handle.
#
# Args:       $perseq_FH:          output file handle to print per-sequence gap info to, undef to not print perseq info
#             $pergap_FH:          output file handle to print per-gap info to
#             $do_matpept:         '1' if we're in mat_peptide mode, else 0
#             $do_perseq_tbl:      '1' to output per sequence gaps as a table, '0' to print a list
#             $do_gap_all:         '1' to output per gap info for all gaps
#             $do_gap_not3:        '1' to output per gap info for gaps that are not a multiple of 3, not for all gaps
#             $do_gap_special:     '1' to output per gap info for special gaps that are possibly causative of a frameshift
#                                  Only 1 of $pergap_all, $pergap_not3, and $pergap_special can be '1'.
#             $mdl_AR:             reference to array of all model names
#             $seq_AR:             reference to array of all sequence names
#             $mdllen_HR:          ref to hash of model lengths
#             $mdl2mft_map_AR:     ref to array that maps each model to a feature
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
  my $nargs_exp = 13;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($perseq_FH, $pergap_FH, $do_matpept, $do_perseq_tbl, $do_gap_all, $do_gap_not3, $do_gap_special, $mdl_AR, $seq_AR, $mdllen_HR, $mdl2mft_map_AR, $refdel_HHAR, $refins_HHAR) = @_;

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
 
  my @gapstr_AA = ();           # [0..$nseq-1][0..$nmft-1]: string describing all gaps for this sequence and this feature
  my @w_gapstr_A = ();          # [0..$i..$nmft-1] max width of $gapstr_AA[$0..nseq-1][$i] for feature $i over all sequences

  my @tot_gap_length_AA = ();   # [0..$nseq-1][0..$nmft-1]: total number of gap positions for this sequence and this feature
  my @w_tot_gap_length_A = ();  # [0..$i..$nmft-1] max width of $tot_gap_length_AA[$0..nseq-1][$i] for feature $i over all sequences

  my @net_gap_length_AA = ();   # [0..$nseq-1][0..$nmft-1]: net number of gap positions for this sequence and this feature
  my @w_net_gap_length_A = ();  # [0..$i..$nmft-1] max width of $net_gap_length_AA[$0..nseq-1][$i] for feature $i over all sequences

  my @mft_gapstr_AH      = ();  # [0..$i..$nmft-1]: hash w/key: gapstring for a single position, value: number of times that gapstring occurs in any of @gapstr_AA for feature $i
  
  my $ch_gapstr         = "string";
  my $ch_tot_gap_length = "tot";
  my $ch_net_gap_length = "net";

  $w_gapstr_A[0]         = length($ch_gapstr);
  $w_tot_gap_length_A[0] = length($ch_tot_gap_length);
  $w_net_gap_length_A[0] = length($ch_net_gap_length);
  %{$mft_gapstr_AH[0]}   = ();

  my $width_seq = length("#accession");
  my $mft_idx = 0;

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
    my $nmdl = scalar(@{$mdl2mft_map_AR});
    my $gapstr = "";
    my $substr;
    my $tot_gap_length = 0;
    my $net_gap_length = 0;
    $mft_idx = 0;

    for(my $h = 0; $h < $nmdl; $h++) { 
      my $mdl = $mdl_AR->[$h];
      my $ndel = (exists $refdel_HHAR->{$mdl}{$seq}) ? scalar(@{$refdel_HHAR->{$mdl}{$seq}}) : 0;
      my $nins = (exists $refins_HHAR->{$mdl}{$seq}) ? scalar(@{$refins_HHAR->{$mdl}{$seq}}) : 0;
      my $mdl_is_final = (($h == ($nmdl-1)) || ($mdl2mft_map_AR->[($h+1)] != $mdl2mft_map_AR->[$h])) ? 1 : 0; 
      my $mdl_is_first = (($h == 0)         || ($mdl2mft_map_AR->[($h-1)] != $mdl2mft_map_AR->[$h])) ? 1 : 0; 
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
                $mft_gapstr_AH[$mft_idx]{$substr}++;
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
                $mft_gapstr_AH[$mft_idx]{$substr}++;
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
              $mft_gapstr_AH[$mft_idx]{$substr}++;
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
              $mft_gapstr_AH[$mft_idx]{$substr}++;
            }
          }
          $ins_idx++;
          $next_ins_str = ($ins_idx < $nins) ? $refins_HHAR->{$mdl}{$seq}[$ins_idx] : undef;
        }        
      } # end of 'while(defined $next_ins_str || defined $next_del_str) { 
      
      $offset += $mdllen_HR->{$mdl};
      
      if($mdl_is_final) {
        if($gapstr eq "") { $gapstr = "-"; }

        # important to update $gapstr here, before we store it if $do_gap_special is true
        if($do_gap_special) { 
          # printf("calling findSpecialGap with $gapstr\n");
          $gapstr = findSpecialGap($gapstr);
          $net_gap_length = $tot_gap_length;
          if($gapstr ne "?" && $gapstr ne "-") { 
            my @el_A = split(",", $gapstr); 
            foreach my $el (@el_A) { 
              $mft_gapstr_AH[$mft_idx]{$el}++;
            }
          }
        }

        push(@{$gapstr_AA[$i]}, $gapstr);
        if(length($gapstr) > $w_gapstr_A[$mft_idx]) { $w_gapstr_A[$mft_idx] = length($gapstr); }

        push(@{$tot_gap_length_AA[$i]}, $tot_gap_length);
        if(length($tot_gap_length) > $w_tot_gap_length_A[$mft_idx]) { $w_tot_gap_length_A[$mft_idx] = length($tot_gap_length); }
        $tot_gap_length = 0;

        push(@{$net_gap_length_AA[$i]}, $net_gap_length);
        if(length($net_gap_length) > $w_net_gap_length_A[$mft_idx]) { $w_net_gap_length_A[$mft_idx] = length($net_gap_length); }
        $net_gap_length = 0;

        $gapstr = "";
        $offset = 0;
        $mft_idx++;

        if(scalar(@w_gapstr_A)         <= $mft_idx) { $w_gapstr_A[$mft_idx]         = length($ch_gapstr); }
        if(scalar(@w_tot_gap_length_A) <= $mft_idx) { $w_tot_gap_length_A[$mft_idx] = length($ch_tot_gap_length); }
        if(scalar(@w_net_gap_length_A) <= $mft_idx) { $w_net_gap_length_A[$mft_idx] = length($ch_net_gap_length); }
        if(scalar(@mft_gapstr_AH)      <= $mft_idx) { %{$mft_gapstr_AH[$mft_idx]}   = (); }
      }
    }
  }
  my $nmft = $mft_idx;

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
    printf $perseq_FH ("# List of all gaps that may solely explain a %s's length not being a multiple of 3, for each sequence:\n#\n", ($do_matpept) ? "mat_peptide" : "CDS");
  }

  if(! $do_gap_special) { 
    printf $perseq_FH  ("%-*s  ", $width_seq, "#");
    for(my $c = 0; $c < $nmft; $c++) { 
      my $w_cur = $w_tot_gap_length_A[$c] + 2 + $w_net_gap_length_A[$c] + 2 + $w_gapstr_A[$c];
      if($c > 0) { print $perseq_FH "  "; }
      printf $perseq_FH ("%-*s", $w_cur, ($do_matpept) ? ("mat_peptide#" . ($c+1)) : ("CDS#" . ($c+1)));
    }
    print $perseq_FH "\n";
    
    # output line 2 (dashes under line 1 of column headers)
    printf $perseq_FH ("%-*s  ", $width_seq, "#");
    for(my $c = 0; $c < $nmft; $c++) { 
      my $w_cur = $w_tot_gap_length_A[$c] + 2 + $w_net_gap_length_A[$c] + 2 + $w_gapstr_A[$c];
      if($c > 0) { print $perseq_FH "  "; }
      printf $perseq_FH ("%-*s", $w_cur, monocharacterString($w_cur, "="));
    }    
    print $perseq_FH "\n";
  }

  # output line 3 of the column headers:
  printf $perseq_FH ("%-*s  ", $width_seq, "#accession");
  for(my $c = 0; $c < $nmft; $c++) { 
    if($c > 0) { print $perseq_FH "  "; }
    if(! $do_gap_special) { 
      printf $perseq_FH ("%-*s  ", $w_tot_gap_length_A[$c], $ch_tot_gap_length);
      printf $perseq_FH ("%-*s  ", $w_net_gap_length_A[$c], $ch_net_gap_length);
      printf $perseq_FH ("%-*s", $w_gapstr_A[$c], $ch_gapstr);
    }
    else { 
      printf $perseq_FH ("%-*s", $w_gapstr_A[$c], ($do_matpept) ? ("mat_peptide#" . ($c+1)) : ("CDS#" . ($c+1)));
    }
  }
  print $perseq_FH "\n";

  # output line 4 (dashes under line 3 of column headers)
  printf $perseq_FH ("%-*s  ", $width_seq, "#" . monocharacterString($width_seq-1, "-"));
  for(my $c = 0; $c < $nmft; $c++) { 
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
    for(my $c = 0; $c < $nmft; $c++) { 
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
    print  $perseq_FH ("# The table includes information on all gaps that exist between all pairwise alignments of\n");
    printf $perseq_FH ("# the reference %s and the predicted homologous %s for each sequence.\n", ($do_matpept) ? "mat_peptide" : "CDS", ($do_matpept) ? "mat_peptide" : "CDS");
  }
  elsif($do_gap_not3) { 
    print  $perseq_FH ("# The table includes information on all gaps of lengths that are not multiples of 3 that exist\n");
    printf $perseq_FH ("# between all pairwise alignments of the reference %s and the predicted homologous %s for each sequence.\n", ($do_matpept) ? "mat_peptide" : "CDS", ($do_matpept) ? "mat_peptide" : "CDS");
  }
  else { 
    print $perseq_FH ("# The table includes information on some gaps that can solely explain a feature (CDS or mat_peptide) not being a multiple of length 3.\n");
    print $perseq_FH ("# This is (probably) not an exhaustive list of all such gaps.\n");
    print $perseq_FH ("# Specifically it is only gaps X in a feature Y, such that the following criteria are met:\n");
    print $perseq_FH ("#   - length of feature Y is not a multiple of 3\n");
    print $perseq_FH ("#   - if you remove only X from list of all gaps, total feature length of Y becomes a multiple of 3\n");
    print $perseq_FH ("#       with length difference of D with reference feature\n");
    print $perseq_FH ("#    - there are no other gaps Z such that if you remove only Z then length of Y becomes a multiple\n");
    print $perseq_FH ("#       of 3 with length difference D2 from reference where D2 < D.\n");
  }
  print $perseq_FH ("#\n");
  if($do_gap_all || $do_gap_not3) { 
    printf $perseq_FH ("# There are 3 columns under each header \"%s#<n> (%s)\" named \"tot\", \"net\",\n", ($do_matpept) ? "mat_peptide" : "CDS", ($do_gap_all) ? "all gaps" : "gaps %3 != 0");
    print $perseq_FH ("# and \"string\".\n");
    print $perseq_FH ("# The \"tot\" columns list the total number of gap positions in either sequence in the pairwise alignment.\n");
    print $perseq_FH ("# The \"net\" columns list the net number of the listed gaps in the pairwise alignment; this is the number\n");
    print $perseq_FH ("#   of gaps in the reference sequence minus the number of gaps in the current sequence (inserts minus deletes)\n");
    print $perseq_FH ("# The \"string\" columns include a list of <n> tokens, each of which describes a gap of length >= 1 nucleotide.\n");
  }
  print  $perseq_FH ("#\n");
  print  $perseq_FH ("# Tokens are in the form: <char><position><length>\n");
  printf $perseq_FH ("#   <char>     is 'I' for an insertion relative to the reference %s (gap in reference sequence)\n", $do_matpept ? "mat_peptide" : "CDS");
  printf $perseq_FH ("#              or 'D' for a  deletion  relative to the reference %s (gap in current sequence)\n", $do_matpept ? "mat_peptide" : "CDS");
  print  $perseq_FH ("#   <position> is the nucleotide position of the gap in reference coordinates.\n");
  print  $perseq_FH ("#              For insertions this is the reference position after which the insertion occurs.\n");
  print  $perseq_FH ("#              For deletions  this is the first reference position for this deletion.\n");
  print  $perseq_FH ("#   <length>   length of the gap in nucleotides.\n");
  print  $perseq_FH ("#              For insertions this is the number of nucleotides inserted relative to the reference\n");
  print  $perseq_FH ("#              For deletions  this is the number of reference positions deleted.\n");
  print  $perseq_FH ("#\n");
  if($do_gap_special) { 
    print  $perseq_FH ("#\n");
    printf $perseq_FH ("# \"-\" tokens indicate the %s is a multiple of length 3\n", $do_matpept ? "mat_peptide" : "CDS");
    printf $perseq_FH ("# \"?\" tokens indicate the %s is not a multiple of length 3, but that no gaps that satisfy our criteria exist.\n", $do_matpept ? "mat_peptide" : "CDS");
    print  $perseq_FH ("#\n");
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
    print $perseq_FH ("# Specifically, these are counts of gaps X in a feature Y (CDS or mat_peptide), such that the following criteria are met:\n");
    print $perseq_FH ("#   - length of feature Y is not a multiple of 3\n");
    print $perseq_FH ("#   - if you remove only X from list of all gaps, total feature length of Y becomes a multiple of 3\n");
    print $perseq_FH ("#       with length difference of D with reference feature\n");
    print $perseq_FH ("#    - there are no other gaps Z such that if you remove only Z then length of Y becomes a multiple\n");
    print $perseq_FH ("#       of 3 with length difference D2 from reference where D2 < D.\n");
    print $perseq_FH ("#\n");
  }
  my $nprinted = 0;
  for(my $c = 0; $c < $nmft; $c++) { 
    if((scalar(keys %{$mft_gapstr_AH[$c]})) > 0) { 
      foreach my $key (sort keys %{$mft_gapstr_AH[$c]}) { 
        printf $pergap_FH ("%s#" . ($c+1) . " " . $key . " " . $mft_gapstr_AH[$c]{$key} . "\n", ($do_matpept) ? "mat_peptide" : "CDS");
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

# Subroutine: outputSeqAsColumnsPage()
#
# Synopsis:   Output a 'page' of annotation information with sequences as columns
#
# Args:       $header_AR:     reference to array of row headers
#             $out_AAR:       reference to the 2D array of output tokens for
#                             current sequences for the page
#             $ref_out_AR:    reference to array of output tokens for reference
#             $page_idx:      page number
#
# Returns:    void
#
sub outputSeqAsColumnsPage {
  my $sub_name = "outputSeqAsColumnsPage";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($header_AR, $out_AAR, $ref_out_AR, $page_idx) = @_;

  my $nseq = scalar(@{$out_AAR});
  my $AR; # reference to current array we are printing

  my @cwidth_A = (); # max width of each column (max width of all tokens for each sequence) 
  my ($ntok, $el);
  my $nrow = scalar(@{$header_AR});

  # first, make sure all sequences have same number of output tokens as we have row headers,
  # and determine maximum width of all tokens for each sequence
  for(my $i = 0; $i <= $nseq; $i++) { 
    # first column is reference, then come the other seqs
    $AR = ($i == 0) ? $ref_out_AR : \@{$out_AAR->[$i-1]}; 
    $cwidth_A[$i] = 0;
    $ntok = scalar(@{$AR});
    if($ntok != $nrow) { die sprintf("ERROR in $sub_name, we have $nrow headers, but sequence %s has $ntok tokens", $i+1, $ntok); }
    foreach $el (@{$AR}) { 
      $el =~ s/^\s+//; # remove leading whitespace
      $el =~ s/\s+$//; # remove trailing whitespace
      if(length($el) > $cwidth_A[$i]) { 
        $cwidth_A[$i] = length($el);
      }
    }
  }
  for(my $i = 0; $i <= $nseq; $i++) { 
    $cwidth_A[$i] += 2; # add 2 spaces for in-between columns
  }

  # determine max width of all row headers:
  my $hwidth = 0;
  for(my $r = 0; $r < $nrow; $r++) { 
    if(length($header_AR->[$r]) > $hwidth) { 
      $hwidth = length($header_AR->[$r]);
    }
  }
  $hwidth += 2;

  for(my $r = 0; $r < $nrow; $r++) { 
    printf("%-*s", $hwidth, $header_AR->[$r]);
    for(my $i = 0; $i <= $nseq; $i++) { 
      # first column is reference, then come the other seqs
      $AR = ($i == 0) ? $ref_out_AR : \@{$out_AAR->[$i-1]}; 
      $el = $AR->[$r];
      $el =~ s/\s+//g;
      printf("%*s", $cwidth_A[$i], $AR->[$r]);
    }
    printf("\n");
  }
  print "#\n";
  printf("# end of page %d\n", $page_idx);
  print "#\n";
  return;
}

###########################################################

