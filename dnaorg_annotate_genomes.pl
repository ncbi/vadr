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
my $esl_reformat   = "esl-reformat";
my $esl_epn_translate  = "/home/nawrocke/notebook/15_1118_dnaorg_annotate_genomes_translation/git-esl-epn-translate/esl-epn-translate.pl";
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
$usage .= "  -nocorrect   : do not correct annotations based on internal start/stop codons in predicted exons/CDS\n";
$usage .= "  -matpept <f> : read mat_peptide info in addition to CDS info, file <f> explains CDS:mat_peptide relationships\n";
$usage .= "  -strict      : require matching annotations to match CDS/exon index\n";
$usage .= "  -nodup       : do not duplicate genome seqs to identify features that span stop..start (for circular genomes)\n";
$usage .= "  -notexon     : do not use exon-specific models\n";
$usage .= "  -onlybuild   : exit after building reference models\n";
$usage .= "  -skipbuild   : skip the build and calibration step because you already did an -onlybuild run\n";
$usage .= "  -model <s>   : use model file <s>, instead of building one\n";
$usage .= "\n OPTIONS THAT ADD EXTRA ERROR TESTS:\n";
$usage .= "  -oseq <s>     : identify origin seq <s> in genomes, put \"|\" at site of origin, e.g. \"TAATATT\\|AC\"\n";
$usage .= "                  <s> must be a string consisting of only A,C,G,T and | characters. No regular expressions allowed.\n"; 
$usage .= "                  Note that \"|\" must be escaped, i.e. \"\\|\"; the \"|\" can be any char, incl. first or last\n";
$usage .= "                  This option is relevant only for circular genomes, e.g. Maize streak virus\n";
$usage .= "  -lfile <f>    : TO BE IMPLEMENTED\n";
$usage .= "  -idthresh <x> : TO BE IMPLEMENTED\n";
#$usage .= "  -lfile <f>    : read minimum and maximum length allowed per feature from file <f>\n";
#$usage .= "  -idthresh <x> : minimum allowed fractional identity per feature is <x>\n";
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

foreach my $x ($hmmbuild, $hmmpress, $nhmmscan, $cmbuild, $cmcalibrate, $cmfetch, $cmpress, $cmscan, $esl_epn_translate, $esl_ssplit) { 
  if(! -x $x) { die "ERROR executable file $x does not exist (or is not executable)"; }
}

# general options:
my $do_nocorrect   = 0; # set to '1' if -nocorrect  enabled, do not look for internal start/stop codons and update annotation if found
my $do_matpept     = 0; # set to '1' if -matpept    enabled, genome has a single polyprotein, use mat_peptide info, not CDS
my $matpept_infile = undef; # defined if -matpept   enabled, the input file that describes relationship between CDS and mat_peptides
my $do_strict      = 0; # set to '1' if -strict     enabled, matching annotations must be same index CDS+exon, else any will do
my $do_nodup       = 0; # set to '1' if -nodup      enabled, do not duplicate each genome, else do 
my $do_notexon     = 0; # set to '1' if -noexon     enabled, do not use exon-specific models, else do
my $do_onlybuild   = 0; # set to '1' if -onlybuild  enabled, exit after building the model
my $do_skipbuild   = 0; # set to '1' if -skipbuild  enabled, skip the build step
my $in_model_db    = undef; # defined if -model <s> enabled, use <s> as the model file instead of building one
# options for adding additional errors/tests
my $origin_seq     = undef; # defined if -oseq         enabled, the origin sequence string
my $length_infile  = undef; # defined if -lfile <s>    enabled, the input file defining allowed lengths
my $idthresh       = undef; # defined if -idthresh <x> enabled, the minimum allowed fractional identity
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

&GetOptions("nocorrect"  => \$do_nocorrect,
            "matpept=s"  => \$matpept_infile,
            "strict"     => \$do_strict,
            "nodup"      => \$do_nodup,
            "notexon"    => \$do_notexon,
            "onlybuild"  => \$do_onlybuild,
            "skipbuild"  => \$do_skipbuild,
            "model=s"    => \$in_model_db,
            "oseq=s"     => \$origin_seq,
            "lfile=s"    => \$length_infile,
            "idthresh=s" => \$idthresh,
            "c"          => \$do_concise,
            "seqrow"     => \$do_seqrow,
            "seqcol"     => \$do_seqcol,
            "nseqcol=s"  => \$nseqcol,
            "nomdlb"     => \$do_nomdlb,
            "noexist"    => \$do_noexist,
            "nobrack"    => \$do_nobrack,
            "nostop"     => \$do_nostop,
            "nofid"      => \$do_nofid,
            "noss3"      => \$do_noss3,
            "noolap"     => \$do_noolap,
            "noexp"      => \$do_noexp,
            "fullolap"   => \$do_fullolap,
            "fulladj"    => \$do_fulladj,
            "hmmer"      => \$do_hmmer,
            "hmmenv"     => \$do_hmmenv,
            "iglocal"    => \$do_iglocal,
            "cslow"      => \$do_cslow, 
            "cfarm"      => \$do_cfarm,
            "sfarm=s"    => \$sfarm_njobs,
            "swait=s"    => \$sfarm_wait,
            "skipfetch"  => \$do_skipfetch,
            "skipscan"   => \$do_skipscan,
            "skipaln"    => \$do_skipaln) ||
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
if(defined $matpept_infile) { 
  $do_matpept = 1;
  $opts_used_short .= "-matpept $matpept_infile";
  $opts_used_long  .= "# option:  using mat_peptide info, CDS:mat_peptide relationships explained in $matpept_infile [-matpept]\n";
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
if(defined $origin_seq) { 
  $opts_used_short .= "-oseq $origin_seq ";
  $opts_used_long  .= "# option:  searching for origin sequence of $origin_seq [-oseq]\n";
}
if(defined $length_infile) { 
  $opts_used_short .= "-lfile $length_infile ";
  $opts_used_long  .= "# option:  enforcing length minimums and maximums read from file [-lfile]\n";
}
if(defined $idthresh) { 
  $opts_used_short .= "-idthresh $idthresh ";
  $opts_used_long  .= "# option:  enforcing minimum fractional identity of $idthresh [-idthresh]\n";
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
  $opts_used_long  .= "# option:  waiting $sfarm_wait minutes for farm cmscan jobs to finish [-swait]\n";
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
# check for UNEXPECTED option value/combinations, we haven't coded for these yet, but we may be able to
if($do_matpept && (! $do_nodup)) { 
  die "ERROR -matpept requires -nodup (currently, this can be relaxed but you may want to look in code for 'matpept/nodup combo' comments)"
}

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

############################################
# parse the -matpept input file if necessary
###########################################
my @cds2pmatpept_AA = (); # 1st dim: cds index, 2nd dim: value array of primary matpept indices that comprise this CDS
my @cds2amatpept_AA = (); # 1st dim: cds index, 2nd dim: value array of all     matpept indices that comprise this CDS
if(defined $matpept_infile) { 
  matpeptParseInfile($matpept_infile, \@cds2pmatpept_AA, \@cds2amatpept_AA);
}

#####################
# parse the list file
# NOTE: the table and length files get created external to this program and have prescribed names
#####################
my @accn_A        = (); # array of accessions
my %accn_exists_H = (); # hash of accessions, key is accession, value is always '1', used only to check for duplicates
my $ref_accn      = undef; # changed to <s> with -ref <s>
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
$ref_accn = $accn_A[0];

##################################
# parse the table and length files
##################################
my %cds_tbl_HHA = ();   # CDS data from .cds.tbl file
                        # hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %mp_tbl_HHA = ();   # mat_peptide data from .matpept.tbl file
                        # hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %totlen_H = (); # key: accession, value length read from length file
my $mft_tbl_HHAR = undef; # 'main feature' (mft) table data, a reference to either cds_tbl_HHA or mft_tbl_HHA

parseLength($length_file, \%totlen_H);

if($do_matpept) { 
  parseTable($matpept_tbl_file, \%mp_tbl_HHA);
}
parseTable($cds_tbl_file, \%cds_tbl_HHA);

if($do_matpept) { 
  # validate the CDS:mat_peptide relationships that we read from $matpept_infile
  matpeptValidateCdsRelationships(\@cds2pmatpept_AA, \%{$cds_tbl_HHA{$ref_accn}}, \%{$mp_tbl_HHA{$ref_accn}});
}

if($do_matpept) { 
  $mft_tbl_HHAR = \%mp_tbl_HHA;
}
else { 
  $mft_tbl_HHAR = \%cds_tbl_HHA;
}

# if we're in matpept mode, determine which peptides should be adjacent to each other at 
# the nucleotide level

#######################
# variable declarations
#######################
my $strand_str;                # +/- string for all main features (CDS or mat_peptide) for an accession: e.g. '+-+': 1st and 3rd CDS are + strand, 2nd is -

# reference information on reference accession, first accession read in ntlist file
my $ref_label_str     = undef; # label string for reference accn
my $ref_nmft          = 0;     # number of main features (mft) (CDS or mat_peptide) in reference
my $ref_strand_str    = "";    # strand string for reference 
my @ref_mft_len_A     = ();    # [0..$i..$ref_nmft-1]: length of each reference main feature (CDS or mat_peptide)
#my @ref_mft_len_tol_A = ();   # [0..$i..$ref_nmft-1]: length tolerance, any gene that is within this fraction of the length of the ref gene is a match
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

# error codes variable declarations
my $nerrcodes_pf = 12;
my $nerrcodes_ps =  1;
my @err_pf_idx2code_A = (); # per-feature  map of error code array idx to code
my %err_pf_idx2msg_H  = (); # per-feature  map of error code array idx to error message
my %err_pf_code2idx_H = (); # per-feature  map of error code code to array idx
my @err_ps_idx2code_A = (); # per-sequence map of error code array idx to code
my %err_ps_idx2msg_H  = (); # per-sequence map of error code array idx to error message
my %err_ps_code2idx_H = (); # per-sequence map of error code code to array idx

$err_pf_idx2code_A[0]  = "nop";
$err_pf_idx2code_A[1]  = "nm3";
$err_pf_idx2code_A[2]  = "bd5";
$err_pf_idx2code_A[3]  = "bd3";
$err_pf_idx2code_A[4]  = "olp";
$err_pf_idx2code_A[5]  = "str";
$err_pf_idx2code_A[6]  = "stp";
$err_pf_idx2code_A[7]  = "ajb";
$err_pf_idx2code_A[8]  = "aja";
$err_pf_idx2code_A[9] = "trc";
$err_pf_idx2code_A[10] = "ext";
$err_pf_idx2code_A[11] = "ntr";

$err_pf_idx2msg_H{"nop"} = "unable to identify homologous feature";
$err_pf_idx2msg_H{"nm3"} = "length of nucleotide feature is not a multiple of 3";
$err_pf_idx2msg_H{"bd5"} = "alignment to reference does not extend to 5' boundary of reference";
$err_pf_idx2msg_H{"bd3"} = "alignment to reference does not extend to 5' boundary of reference";
$err_pf_idx2msg_H{"olp"} = "feature does not overlap with same set of features as in reference";
$err_pf_idx2msg_H{"str"} = "predicted CDS start position is not beginning of ATG start codon";
$err_pf_idx2msg_H{"stp"} = "predicted CDS stop  position is not end of valid stop codon (TAG|TAA|TGA)";
$err_pf_idx2msg_H{"ajb"} = "mature peptide is not adjacent to same set of mature peptides before it as in reference";
$err_pf_idx2msg_H{"aja"} = "mature peptide is not adjacent to same set of mature peptides after it as in reference";
$err_pf_idx2msg_H{"trc"} = "in-frame stop codon exists 5' of stop position predicted by homology to reference";
$err_pf_idx2msg_H{"ext"} = "first in-frame stop codon exists 3' of stop position predicted by homology to reference";
$err_pf_idx2msg_H{"ntr"} = "mature peptide is not translated because its CDS has an in-frame stop 5' of the mature peptide's predicted start";

my $e;
for($e = 0; $e < $nerrcodes_pf; $e++) { 
  $err_pf_code2idx_H{$err_pf_idx2code_A[$e]} = $e;
}

$err_ps_idx2code_A[0]    = "ori";
$err_ps_idx2msg_H{"ori"} = "there is not exactly 1 occurrence of origin sequence";
for($e = 0; $e < $nerrcodes_ps; $e++) { 
  $err_ps_code2idx_H{$err_ps_idx2code_A[$e]} = $e;
}

#####################################################
# Fetch all genome sequences, including the reference
#####################################################
my $naccn = scalar(@accn_A);
my $gnm_fetch_file = $out_root . ".fg.idfetch.in";
my $gnm_fasta_file = $out_root . ".fg.fa";
my @seq_accn_A = (); # [0..$naccn-1] name of genome fasta sequence for each accn
my $seq_accn;     # temp fasta sequence name
my %accn2seq_accn_H = (); # key: accession (from @accn_A), value, corresponding accession from @seq_accn_A
my $fetch_string = undef;
my $ref_seq_accn; # name of fasta sequence for reference
my $ref_totlen;   # total length of reference
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
  $accn2seq_accn_H{$accn} = $seq_accn;
  if($a == 0) { 
    $ref_seq_accn = $seq_accn; 
    $ref_totlen   = $totlen_H{$accn};
  }
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
if(! exists ($mft_tbl_HHAR->{$ref_accn})) { die sprintf("ERROR no %s information stored for reference accession", ($do_matpept) ? "mat_peptide" : "CDS"); }
(undef, undef, undef, undef, undef, $ref_strand_str) = getStrandStats($mft_tbl_HHAR, $ref_accn);
getLengthStatsAndCoordStrings(\%{$mft_tbl_HHAR->{$ref_accn}}, \@ref_mft_len_A, \@ref_mft_coords_A);
getQualifierValues($mft_tbl_HHAR, $ref_accn, "product", \@ref_mft_product_A);
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
my @mft2nexon_map_A = ();   # [0..h..$nmdl-1]: $nexons: model ($m+1) models an exon/mat_peptide that is one of $nexons exons/mat_peptides for this CDS
my @mdl_is_first_A = ();    # [0..h..$nmdl-1]: '1' if model ($m+1) is the first one for cds $mdl2mft_map_A[$h], else 0
my @mdl_is_final_A = ();    # [0..h..$nmdl-1]: '1' if model ($m+1) is the final one for cds $mdl2mft_map_A[$h], else 0
my @mdl_is_primary_A = ();  # [0..h..$nmdl-1]: '1' if model ($m+1) is for a primary matpept, else 0, only created/used if $do_matpept
my @mft2first_mdl_A = ();   # [0..$c..ncds-1]: $h, first exon of feature $f+1 is modeled by model $m+1
my @mft2final_mdl_A = ();   # [0..$c..ncds-1]: $h, final exon of feature $f+1 is modeled by model $m+1
my @mdl_A = ();             # [0..$nmdl-1]: array of model names, also name of stockholm alignments used to build those models
my @mft_out_tiny_A    = (); # [0..$ref_nmft-1]: array of very abbreviated model main feature names to print
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
  startStopsFromCoords($ref_mft_coords_A[$i], $ref_totlen, $do_nodup, \@starts_A, \@stops_A, undef, \$nexons);
  my $strand = substr($ref_strand_str, $i, 1);

  push(@ref_nexons_A, $nexons);
  $ref_tot_nexons += $nexons;

  # if we're on the negative strand, reverse the arrays, they'll be in the incorrect order
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
      my $tiny  = $short;
      $tiny =~ s/\s+//g; # remove whitespace
      if($act_nexons > 1) { $short .= sprintf(" [$act_nexons %s; $strand]", ($do_matpept) ? "segments" : "exons"); }
      else                { $short .= sprintf(" [single %s; $strand]",      ($do_matpept) ? "segment"  : "exon"); }
      push(@mft_out_tiny_A,    $tiny);
      push(@mft_out_short_A,   $short);
      push(@mft_out_product_A, $ref_mft_product_A[$i]);
    }
    push(@mdl_A, $cur_name_root);

    $mdllen_H{$cur_name_root} = $mdllen;

    # now append the named alignment to the growing stockholm alignment database $all_stk_file
    $cmd = "cat $cur_named_stkfile";
    if($nmdl == 0) { $cmd .= " >  $all_stk_file"; }
    else           { $cmd .= " >> $all_stk_file"; }
    runCommand($cmd, 0);
    push(@mdl2mft_map_A,  $i);
    push(@mft2exon_map_A, $e);
    push(@mft2nexon_map_A, $nexons);
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

# fill @mft2mdl_map_AA
my @mft2mdl_map_AA = ();  # [0..$c..$ref_nmft-1][0..$e..$nexons_for_this_feature-1]: value is index of mdl that encodes exon $e for feature $c
my $mft_idx = -1;
my $mdl_idx = 0;
for(my $h = 0; $h < $nmdl; $h++) { 
  if($mdl_is_first_A[$h]) {
    $mft_idx++;
    $mdl_idx = 0;
  }
  $mft2mdl_map_AA[$mft_idx][$mdl_idx] = $h;
  $mdl_idx++;
}

### TEMP
##for (my $z = 0; $z < $ref_nmft; $z++) { 
##  for (my $zz = 0; $zz < scalar(@{$mft2mdl_map_AA[$z]}); $zz++) { 
##    printf("mft2mdl_map_AA[$z][$zz] = $mft2mdl_map_AA[$z][$zz]\n");
##  }
##  printf("\n");
##}
##exit 0;

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
    #runCmscan($cmscan, $do_iglocal, ($do_skipscan || $do_skipaln), 0, $model_db, $gnm_fasta_file, $tblout, $stdout);
    runCmscan($cmscan, $do_iglocal, ($do_skipscan || $do_skipaln), 0, $model_db, $gnm_fasta_file, $tblout, undef); #undef stdout --> don't save stdout to save space
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
wrapperCombineExonsIntoCDS($dir, "predicted", \@mdl_A, \@accn_A, \@mdl2mft_map_A, \@mft2mdl_map_AA, \@pred_mft_fafile_A);

#####################################################
# COMBINE MULTI-MAT_PEPTIDE SEQUENCES INTO SINGLE CDS
#####################################################
my @matpept_cds_pass_fail_AH = ();
my @matpept_cds_stop_codon_AH = ();
my @pred_cds_mft_fafile_A = (); # array of predicted CDS sequence files, filled below
if($do_matpept) { 
  wrapperCombineExonsIntoCDS($dir, "predicted", \@mdl_A, \@accn_A, \@mdl2mft_map_A, \@cds2pmatpept_AA, \@pred_cds_mft_fafile_A);
  push(@pred_mft_fafile_A, @pred_cds_mft_fafile_A);

  # keep track of which CDS sequences are valid (pass all adjacency tests), we will only do the translation tests
  # (i.e. look for in-frame stops) for valid CDS
  my $ncds = scalar(@cds2pmatpept_AA);
  for(my $c = 0; $c < $ncds; $c++) { 
    %{$matpept_cds_pass_fail_AH[$c]} = ();
    for(my $a = 0; $a < $naccn; $a++) { 
      my $accn     = $accn_A[$a];
      my $seq_accn = $seq_accn_A[$a];
      my ($cds_start, $cds_stop, undef, undef, $cds_len, $cds_start_codon, $cds_stop_codon, $cds_pass_fail) = matpeptCheckCdsRelationships($sqfile, $seq_accn, $totlen_H{$accn}, \@mdl_A, \@{$cds2pmatpept_AA[$c]}, \%p_start_HH, \%p_stop_HH, undef, \%p_strand_HH, \@mft2first_mdl_A, \@mft2final_mdl_A);
      $matpept_cds_pass_fail_AH[$c]{$seq_accn}  = $cds_pass_fail;
      $matpept_cds_stop_codon_AH[$c]{$seq_accn} = $cds_stop_codon;
      # printf("0 matpept_cds_stop_codon_AH[$c]{$seq_accn}: $matpept_cds_stop_codon_AH[$c]{$seq_accn}\n");
    }
  }
}

###########################################################################
# Before we translate the sequences, we take some extra steps to correct
# stop predictions:
# - for all features with a valid start codon: look for in-frame stops
# - for all features with a valid start codon but no in-frame stops up
#   until predicted stop position (thus, stop is an invalid stop) look
#   downstream for first in-frame stop after predicted stop position
#
# CORRECT PREDICTIONS 
# - for all features: look for in-frame stops
# - for CDS/mature_peptides with 
# - combine new exons into new CDS
###########################################################################

($seconds, $microseconds) = gettimeofday();
$start_time = ($seconds + ($microseconds / 1000000.));
printf("%-65s ... ", sprintf("# Translating predicted %s to identify internal starts/stops", ($do_matpept) ? "mat_peptides" : "CDS"));

my $source_accn;   # name of CDS/mat_peptide sequence that was translated
my @corr_mft_stop_AH = ();   # [0..$i..nmft-1], each element is a hash with keys $key as sequence accessions and values 
                             # are number of nucleotides that the prediction of the stop coordinate should be corrected
                             # based on an esl-translate translation of the predicted CDS/mat_peptide sequence, values can be negative
                             # or positive
my %c_stop_HH  = ();        # corrected stop positions of hits,  start with a copy of p_stop_HH
my @did_corr_exon_stop_AH = ();   # [0..$i..ref_nexons-1], each element is a hash with keys $key as sequence accessions and values 
                                  # are '1' if this exon's stop position was corrected
my %corr_fafile_H = ();     # hash of names of fasta files for corrected exon sequences, keys: model name from @mdl_A, value: name of file
my @corr_mft_fafile_A = (); # array of corrected CDS/mat_peptide sequence files, filled below

my $ref_nmp_and_cds = $ref_nmft;
if($do_matpept) { 
  # add in CDS in matpept mode
  $ref_nmp_and_cds += scalar(@cds2pmatpept_AA); 
}
my @source_accn_ext_check_AA     = (); # [0..$c..$ref_nmp_and_cds-1][0..$j..$n] for each feature $c, array of source 
                                       # accessions to check for 'ext' error, b/c no in-frame stop within 
                                       # predicted start..stop
my @source_accn_trc_error_AH     = (); # [0..$c..$ref_nmp_and_cds-1][0..$j..$n] for each feature $c, array of source 
                                       # accessions that have 'trc' error, b/c in-frame stop exists between predicted start..stop
my @source_accn_ext_error_AH     = (); # [0..$c..$ref_nmp_and_cds-1][0..$j..$n] for each feature $c, array of source 
                                       # accessions that have 'ext' error, b/c no in-frame stop exists between predicted start..stop
my @source_accn_nst_error_AH     = (); # 1st dim [0..$c..$ref_nmp_and_cds-1], key: source accn, value: 1 if this 
                                       # feature has no in-frame stop
                                       # 2nd of 2 reasons we shouldn't translate it
                                       # these should be rare, these are a subset of the accessions in
                                       # @source_accn_ext_check_AA, the ones that had no valid stop in-frame
                                       # for the rest of the sequence (possibly duplicated)
my @source_accn_str_error_AH     = (); # 1st dim [0..$c..$ref_nmp_and_cds-1], key: source accn, value: 1 if this 
                                       # feature has an invalid start 
                                       # 1st of 2 reasons we shouldn't translate it
my @source_accn_stp_error_AH     = (); # 1st dim [0..$c..$ref_nmp_and_cds-1], key: source accn, value: 1 if this 
                                       # feature's predicted stop is not a valid stop codon
                                       # 1st of 2 reasons we shouldn't translate it

# Translate predicted CDS/mat_peptide sequences using esl-epn-translate to identify 
# in-frame stop codons.

for(my $c = 0; $c < $ref_nmp_and_cds; $c++) { 
  my $is_matpept = ($do_matpept && ($c < $ref_nmft)) ? 1 : 0; # if $is_matpept is 0, then we are a CDS
  #printf("HEYAD do_matpept: $do_matpept is_matpept: $is_matpept\n");

  # initialize data structures
  @{$source_accn_ext_check_AA[$c]} = ();
  %{$source_accn_trc_error_AH[$c]} = ();
  %{$source_accn_ext_error_AH[$c]} = ();
  %{$source_accn_nst_error_AH[$c]} = ();
  %{$source_accn_str_error_AH[$c]} = ();
  %{$source_accn_stp_error_AH[$c]} = ();
  %{$corr_mft_stop_AH[$c]}         = ();
  
  my $cur_fafile = $pred_mft_fafile_A[$c];
  
  # translate into AA sequences
  my $tmp_esl_epn_translate_output = $cur_fafile;
  $tmp_esl_epn_translate_output =~ s/\.mp/\.esl-epn-translate/; 
  $tmp_esl_epn_translate_output =~ s/\.cds/\.esl-epn-translate/;
  $tmp_esl_epn_translate_output =~ s/\.fa$//;
  
  my $other_options = "";
  $cmd = $esl_epn_translate . " $other_options -firststop $cur_fafile > $tmp_esl_epn_translate_output";
  runCommand($cmd, 0);
  
  # parse esl-epn-translate output
  open(IN, $tmp_esl_epn_translate_output) || die "ERROR unable to open $tmp_esl_epn_translate_output for reading";

  while(my $line = <IN>) { 
    # example line
    #HQ693465/1-306 1 1 304
    if($line =~ /^(\S+)\/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/) { 
      my ($source_accn, $coords, $start_is_valid, $stop_is_valid, $first_stop_first_posn) = ($1, $2, $3, $4, $5);
      
      # A wrinkle here related to CDS comprised of mat_peptide
      # sequences, the final codon of the final mat_peptide sequence
      # for the CDS *SHOULD NOT* encode a stop codon. This is NCBI
      # annotation protocol. The three nucleotides 5' of the final
      # mat_peptide stop are in fact the stop codon. This presents a
      # problem for us here, because we 'create' the full CDS by
      # simply concatenating all the mat_peptide coding sequences
      # together. However, when we do that the resulting CDS will
      # necessarily be missing the stop codon. NCBI CDS annotation
      # however, *does* include the stop codon, so this is an
      # inconsistency.
      #
      # What this means is that normally for a mat_peptide comprised
      # CDS, the $stop_is_valid value will be 0 because
      # esl-epn-translate has not checked the correct 3 nucleotides to
      # see if it is a valid stop codon. To deal with this, we've
      # previously stored whether the final mat_peptide's next 3
      # nucleotides encode a valid stop codon or not. We have that in
      # $matpept_cds_stop_codon_AH() which was filled following a call
      # to matpeptCheckCdsRelationships() which manually checks the
      # next 3 nucleotides 5' of the stop position. So we use that
      # value to override $stop_is_valid here.
      if($do_matpept && (! $is_matpept)) { # a mat_peptide comprised CDS
        my $seq_accn = $accn2seq_accn_H{$source_accn};
        my $cur_stop_codon = $matpept_cds_stop_codon_AH[($c - $ref_nmft)]{$seq_accn};
        #printf("REEXAMINING stop_is_valid, was: $stop_is_valid, now: ");
        $stop_is_valid = ((defined $cur_stop_codon) && 
                          ($cur_stop_codon eq "TAG" || $cur_stop_codon eq "TAA" || $cur_stop_codon eq "TGA"))
            ? 1 : 0;
        #printf("$stop_is_valid\n");
      }

      # determine what first_stop_first_posn would be if the final codon is the first in-frame stop
      my @start_stop_A = split(",", $coords);
      my $cds_len = 0;
      foreach my $start_stop (@start_stop_A) { 
        if($start_stop =~ /(\-?\d+)\-(\-?\d+)/) { 
          $cds_len += abs($1 - $2) + 1;
          my ($tmp_start, $tmp_stop) = ($1, $2);
          if(($tmp_start < 0 && $tmp_stop > 0) || ($tmp_start > 0 && $tmp_stop < 0)) { $cds_len--; } # correct for off-by-one error due to fact that posn 0 is not a nt
        }
        else { 
          die "ERROR unable to parse coordinates ($coords) read from esl-epn-translate.pl output";
        }
      }
      my $final_codon_first_posn = $cds_len - 2;
      my $early_inframe_stop     = (($first_stop_first_posn != 0) && ($first_stop_first_posn != $final_codon_first_posn)) ? 1 : 0;

      # We now have all of the relevant data on the current
      # CDS/mat_peptide sequence and we need to use that to determine
      # what errors each sequence should throw (at least for those
      # that can tell should be thrown at this point in the
      # processing) as well as any corrections to stop predictions
      # that we should make prior to translation (trc errors are
      # thrown when this happens).
      # 
      # There are 4 relevant variables that dictate which errors
      # should be thrown/checked for later and whether a stop
      # correction should be made. The table gives all possible
      # combinations of values of those variables and lists the
      # outcome of each possibility.
      #
      # Variables we know from earlier processing:
      # $is_matpept:     '1' if current feature is a mature peptide, '0' if it is a CDS
      #
      # Variables derived from esl-epn-translate output we're currently parsing
      # $start_is_valid:     '1' if current feature's first 3 nt encode a valid start codon
      # $stop_is_valid:      '1' if current feature's final 3 nt encode a valid stop codon and total feature length is multiple of 3
      # $early_inframe_stop: '1' if an inframe stop exists prior to predicted stop
      #
      # 7 possibilities, each with different outcome (P1-P7):
      #                                                                                       | make correction to |
      # idx | is_matpept | start_is_valid | stop_is_valid | early_inframe_stop ||   errors    | stop coordinate?   |
      # ----|------------|----------------|---------------|--------------------||-------------|--------------------|
      #  P1 |      false |          false |          any  |                any ||         str |                 no |
      #  P2 |      false |           true |        false  |              false ||     stp ext?|     maybe (if ext) |
      #  P3 |      false |           true |        false  |               true ||     stp trc |                yes |
      #  P4 |      false |           true |         true  |              false ||        none |                 no |
      #  P5 |      false |           true |         true  |               true ||         trc |                yes |
      # -----------------------------------------------------------------------||-----------------------------------
      #  P6 |       true |            any |          any  |              false ||        ntr? |                 no |
      #  P7 |       true |            any |          any  |               true ||    trc ntr? |                yes |
      # ------------------------------------------------------------------------------------------------------------
      # 
      # in table above:
      # '*' around value means that it is guaranteed to be that value (e.g. *true* means 'false' is impossible)
      # '?' after error code means that error is possible, we have to check for it later           
      # 'any' means that any value is possible, outcome is unaffected by value
      #
      if(! $is_matpept) { 
        if(! $start_is_valid) { # possibility 1 (P1)
          #printf("HEYA $source_accn $c possibility 1 (str)\n");
          $source_accn_str_error_AH[$c]{$source_accn} = 1;
        }
        else { 
          # $start_is_valid is 1
          if(! $stop_is_valid) { 
            if(! $early_inframe_stop) { 
              # possibility 2 (P2): stp error, need to check for ext error later
              #printf("HEYA $source_accn $c possibility 2 (stp, maybe ext)\n");
              $source_accn_stp_error_AH[$c]{$source_accn} = 1;
              push(@{$source_accn_ext_check_AA[$c]}, $source_accn);
            }
            else { # $early_inframe_stop is 1
              # possibility 3 (P3): stp and trc error
              #printf("HEYA $source_accn $c possibility 3 (trc and stp)\n");
              $source_accn_stp_error_AH[$c]{$source_accn} = 1;
              $source_accn_trc_error_AH[$c]{$source_accn} = 1;
              $corr_mft_stop_AH[$c]{$source_accn} = -1 * ($final_codon_first_posn - $first_stop_first_posn);
              if($final_codon_first_posn == $first_stop_first_posn) { die "ERROR trying to correct a stop by 0 nt"; }
              #printf("HEYAA set corr_mft_stop_AH[$c]{$source_accn} to $corr_mft_stop_AH[$c]{$source_accn}\n");
            }
          } # end of 'if(! $stop_is_valid)'
          else { # $stop_is_valid is 1
            if(! $early_inframe_stop) { 
              #printf("HEYA $source_accn $c possibility 4 (no errors)\n");
              ;# possibility 4 (P4): no errors, do nothing
            }
            else { 
              # possibility 5 (P5): trc error
              #printf("HEYA $source_accn $c possibility 5 (trc)\n");
              $source_accn_trc_error_AH[$c]{$source_accn} = 1;
              $corr_mft_stop_AH[$c]{$source_accn} = -1 * ($final_codon_first_posn - $first_stop_first_posn);
              if($final_codon_first_posn == $first_stop_first_posn) { die "ERROR trying to correct a stop by 0 nt"; }
              #printf("HEYAA set corr_mft_stop_AH[$c]{$source_accn} to $corr_mft_stop_AH[$c]{$source_accn}\n");
            }
          }              
        }
      } # end of 'if(! $is_matpept)'
      else { # $is_matpept is 1 
        if(! $early_inframe_stop) { 
          # possibility 6 (P6): maybe ntr error later, but can't check for it now, do nothing;
          #printf("HEYA $source_accn $c possibility 6 (no error)\n");
        }
        else { # $early_inframe_stop is '1'
          # possibility 7 (P7): trc error, maybe ntr error later, but can't check for it now
          #printf("HEYA $source_accn $c possibility 7 (trc)\n");
          $source_accn_trc_error_AH[$c]{$source_accn} = 1;
          $corr_mft_stop_AH[$c]{$source_accn} = -1 * ($final_codon_first_posn - $first_stop_first_posn);
          if($final_codon_first_posn == $first_stop_first_posn) { die "ERROR trying to correct a stop by 0 nt"; }
          #printf("HEYAA set corr_mft_stop_AH[$c]{$source_accn} to $corr_mft_stop_AH[$c]{$source_accn}\n");
        }
      }
    }
    else { 
      die "ERROR, unable to parse esl-epn-translate.pl output line: $line\n";
    }
  }
  close(IN);
}

# next step is to look at all accessions we just stored in
# @source_accn_ext_check_A, these are sequences that start with a valid start,
# but have no in-frame stop in them, look for an inframe stop codon past the 
# predicted stop. Two possibilities:
# (A) we find one
# (B) we do not find one
#
for(my $c = 0; $c < $ref_nmft; $c++) { 
  foreach $source_accn (@{$source_accn_ext_check_AA[$c]}) { 
    printf("HEYA NEED TO check $source_accn for extension...\n");
  }
}

# now we have all stop corrections, go back and create c_stop data structures based 
# on p_stop data structures and these corrections, this is tricky for multi-exon
# CDS because we have to make sure we correct the proper exon
# first copy p_stop_HH into c_stop_HH, any models which we *DO NOT* do corrections for
# will have c_stop_HH remain equal to p_stop_HH. 
for(my $a = 0; $a < $naccn; $a++) { 
  my $accn     = $accn_A[$a];
  my $seq_accn = $seq_accn_A[$a];
  my $mdl; # name of a model
  for(my $c = 0; $c < $ref_nmft; $c++) { # we don't do CDS in mat_peptide mode, since they will necessarily
                                         # be composed of mat_peptides that we'll correct individually
    my $corr_stop  = (exists $corr_mft_stop_AH[$c]{$accn}) ? $corr_mft_stop_AH[$c]{$accn} : 0;
    #printf("HEYAAA $accn mdl: $c, corr_stop: $corr_stop\n");

    # sanity check
    if($corr_stop  > 0) { die sprintf("ERROR corrected stop ($corr_stop) greater than 0, can't deal with that yet %s feature #%d", $accn,  $c+1);  }
    
    my @cur_mdl_A = ();
    @cur_mdl_A = @{$mft2mdl_map_AA[$c]};
    my $cur_nmdl = scalar(@cur_mdl_A);
    
    # if this is our first time setting c_stop_HH, copy p_stop_HH values to c_stop_HH, 
    # if no correction is made c_stop_HH value will remain as p_stop_HH
    for(my $h = 0; $h < $cur_nmdl; $h++) { 
      my $mdl = $mdl_A[$cur_mdl_A[$h]];
      if(exists $p_stop_HH{$mdl}{$seq_accn} && 
         ((! exists $c_stop_HH{$mdl}{$seq_accn}) || 
          (! defined $c_stop_HH{$mdl}{$seq_accn}))) { 
        $c_stop_HH{$mdl}{$seq_accn} = $p_stop_HH{$mdl}{$seq_accn};
        # printf("HEYAAAA initialized c_stop_HH{$mdl}{$seq_accn} to $c_stop_HH{$mdl}{$seq_accn}\n");
      }
    }
    
    # Two cases: (1) single exon CDS or mat_peptide (2) multiple exon CDS
    if($cur_nmdl == 1) { 
      # easy case: single exon CDS/mat_peptide:
      my $mdl_idx = $cur_mdl_A[0];
      my $mdl     = $mdl_A[$cur_mdl_A[0]];
      if(exists $p_strand_HH{$mdl}{$seq_accn}) { 
        if($corr_stop != 0) { 
          $did_corr_exon_stop_AH[$mdl_idx]{$accn} = 1;
          if($p_strand_HH{$mdl}{$seq_accn} eq "+") { 
            $c_stop_HH{$mdl}{$seq_accn}  = $p_stop_HH{$mdl}{$seq_accn}  + $corr_stop;  # note that $corr_stop  may be 0
            #printf("HEYAAAA c_stop_HH{$mdl}{$seq_accn}: $c_stop_HH{$mdl}{$seq_accn} p_stop_HH: $p_stop_HH{$mdl}{$seq_accn}, corr_stop: $corr_stop\n");
          }
          else { 
            $c_stop_HH{$mdl}{$seq_accn}  = $p_stop_HH{$mdl}{$seq_accn}  - $corr_stop;  # note that $corr_stop  may be 0
          }
        }
      }
    }
    else { # multi exon CDS
      # TODO: skip this if any of the exons are not predicted
      # get temporary array of exon lengths and determine full CDS length
      
      if($corr_stop < 0) { # determine which exon we need to correct the stop for, and correct it:
        my $found_exon = 0;
        my $len_so_far = 0;
        for(my $h = ($cur_nmdl-1); $h >= 0; $h--) { # we go backwards through all models
          my $mdl_idx = $cur_mdl_A[$h];
          my $mdl     = $mdl_A[$mdl_idx];
          my $strand  = $p_strand_HH{$mdl_A[$mdl_idx]}{$seq_accn};
          my $exon_len = ($strand eq "+") ? 
              ($p_stop_HH{$mdl_A[$mdl_idx]}{$seq_accn}  - $p_start_HH{$mdl_A[$mdl_idx]}{$seq_accn} + 1) : 
              ($p_start_HH{$mdl_A[$mdl_idx]}{$seq_accn} - $p_stop_HH{$mdl_A[$mdl_idx]}{$seq_accn}  + 1);
          if((-1 * $corr_stop) < ($len_so_far + $exon_len)) { # correction is in this exon
            if($strand eq "+") { 
              $c_stop_HH{$mdl}{$seq_accn} = $p_stop_HH{$mdl}{$seq_accn} - ((-1 * $corr_stop) - $len_so_far);
              printf("set positive strand c_stop_HH{$mdl}{$seq_accn} to $c_stop_HH{$mdl}{$seq_accn} (p_stop_HH: %d - (-1 * $corr_stop) - len_so_far:$len_so_far h: $h, mdl_idx: $mdl_idx, mdl: $mdl)\n", $p_stop_HH{$mdl}{$seq_accn});
            }
            else { 
              printf("INITIALLY: p_stop_HH{$mdl}{$seq_accn} is $p_stop_HH{$mdl}{$seq_accn}\n");
              $c_stop_HH{$mdl}{$seq_accn}  = $p_stop_HH{$mdl}{$seq_accn} + ((-1 * $corr_stop) - $len_so_far);
              printf("set negative strand c_stop_HH{$mdl}{$seq_accn} to $c_stop_HH{$mdl}{$seq_accn} (p_stop_HH{$mdl}{$seq_accn} is $p_stop_HH{$mdl}{$seq_accn}\n");
            }
            $h = -1; # breaks loop
            $did_corr_exon_stop_AH[$mdl_idx]{$accn} = 1;
            printf("set did_corr_exon_stop_AH[$mdl_idx]{$accn} to $did_corr_exon_stop_AH[$mdl_idx]{$accn}\n");
            $found_exon = 1;
          }
          $len_so_far += $exon_len;
        }
        if(! $found_exon) { 
          die "ERROR unable to find exon for corrected stop in multi-exon CDS/mature_peptide"; 
        }
      } # end of 'if($corr_stop < 0)'
    } # end of 'else' entered if we're a multi-exon CDS
  } # end of loop over CDS
} # end of loop over sequences
($seconds, $microseconds) = gettimeofday();
$stop_time = ($seconds + ($microseconds / 1000000.));
printf("done. [%.1f seconds]\n", ($stop_time - $start_time));

# fetch corrected hits into new files
fetchHits($sqfile, $do_skipaln, "corrected", \@mdl_A, \@seq_accn_A, \%totlen_H, \%p_start_HH, \%c_stop_HH, \%p_strand_HH, $out_root, \%corr_fafile_H);

# combine multi-exon sequences into CDS:
wrapperCombineExonsIntoCDS($dir, "corrected", \@mdl_A, \@accn_A, \@mdl2mft_map_A, \@mft2mdl_map_AA, \@corr_mft_fafile_A);

#########################################
# TRANSLATE PREDICTIONS INTO PROTEIN SEQS
#########################################
# TEMPORARILY SKIP THIS STEP IF WE'RE IN MATPEPT MODE
my @nfullprot_A   = ();  # [0..$i..nmft-1], number of accessions we have a full protein for, for CDS $i
my @fullprot_AH   = ();  # [0..$i..nmft-1], each element is a hash with keys $key as sequence accessions and values 
                         # of number of full length protein sequences we have for CDS/mat_peptide $i for key $key. Values should
                         # always be '1', more than '1' is an error, and we never create a value of '0'.
my @ntruncprot_A   = (); # [0..$i..nmft-1], number of accessions we do not have a full protein for, for CDS $i
my @truncprot_AH   = (); # [0..$i..nmft-1], each element is a hash with keys $key as sequence accessions and values 
                         # of number of non-full length protein sequences we have for CDS/mat_peptide $i for key $key. Values should
                         # always be '1', more than '1' is an error, and we never create a value of '0'.
my @aa_full_files_A = (); # array of protein sequence files we are about to create

($seconds, $microseconds) = gettimeofday();
$start_time = ($seconds + ($microseconds / 1000000.));
printf("%-65s ... ", "# Translating coding sequences into proteins/peptides ");
for(my $c = 0; $c < $ref_nmft; $c++) { 
  my $cur_fafile = ($do_nocorrect) ? $pred_mft_fafile_A[$c] : $corr_mft_fafile_A[$c];
  my $aa_full_fafile  = $cur_fafile;
  if($aa_full_fafile =~ m/\.mp/) { 
    $aa_full_fafile  =~ s/\.mp/\.aa.full/;
  }
  elsif($aa_full_fafile =~ m/\.cds/) { 
    $aa_full_fafile  =~ s/\.cds/\.aa.full/;
  }
  else { 
    die "ERROR, did not find \.mp or \.cds in fasta file name $aa_full_fafile when trying to translate coding sequences into proteins/peptides";
  }

  # translate into AA sequences
  my $opts = "";
  if(! $do_matpept) { 
    $opts = " -reqstart -reqstop ";
  }
  $cmd = $esl_epn_translate . " -endatstop -nostop $opts $cur_fafile > $aa_full_fafile";
  
  runCommand($cmd, 0);
  
#    printf("cmd: $cmd\n");
  
  push(@aa_full_files_A, $aa_full_fafile);
}
($seconds, $microseconds) = gettimeofday();
$stop_time = ($seconds + ($microseconds / 1000000.));
printf("done. [%.1f seconds]\n", ($stop_time - $start_time));

#################################
# CREATE MULTIPLE DNA ALIGNMENTS
#################################
if($do_hmmer) { 
  alignHits($hmmalign, $hmmfetch, $model_db, $do_skipaln, \@mdl_A, \@seq_accn_A, 
            (($do_nocorrect) ? \%pred_fafile_H : \%corr_fafile_H), 
            \%p_start_HH, \%p_fid2ref_HH, \%p_refdel_HHA, \%p_refins_HHA);
}
else { 
  alignHits($cmalign, $cmfetch, $model_db, $do_skipaln, \@mdl_A, \@seq_accn_A, 
            (($do_nocorrect) ? \%pred_fafile_H : \%corr_fafile_H), 
            \%p_start_HH, \%p_fid2ref_HH, \%p_refdel_HHA, \%p_refins_HHA);
}

#####################################
# CREATE MULTIPLE PROTEIN ALIGNMENTS
#####################################
# Make sure that the predictions in the reference of each feature are the same as the GenBank annotation
# If they're not we can't use the current method which fetches the CDS/mat_peptide sequence based on the
# predicion
# TEMPORARILIY SKIP THIS STEP IN MATPEPT MODE
#if(! $do_skipaln) {
if(! $do_skipaln) { 
  ($seconds, $microseconds) = gettimeofday();
  $start_time = ($seconds + ($microseconds / 1000000.));
  printf("%-65s ... ", sprintf("# Creating multiple alignments of protein sequences"));
  my @tmp_ref_act_exon_starts_AA = (); # [0..$nmft-1][0..$nexons-1] start positions of actual annotations of exons for ref accn, $nexons is main-feature (CDS or mat_peptide) specific
  my @tmp_ref_act_exon_stops_AA  = (); # [0..$nmft-1][0..$nexons-1] stop  positions of actual annotations of exons for ref accn, $nexons is main-feature (CDS or mat_peptide) specific
  my $tmp_ref_tot_nexons = 0;
  my $tmp_ref_nmft       = 0;
  
  getActualAnnotations($accn_A[0], $ref_totlen, $do_nodup, $mft_tbl_HHAR, \@tmp_ref_act_exon_starts_AA, \@tmp_ref_act_exon_stops_AA, \$tmp_ref_tot_nexons, \$tmp_ref_nmft);
  for(my $h = 0; $h < $nmdl; $h++) { 
    my $model   = $mdl_A[$h];
    my $mft_i   = $mdl2mft_map_A[$h];
    my $exon_i  = $mft2exon_map_A[$h];
    my $nexons  = $mft2nexon_map_A[$h];
    if(! exists $p_start_HH{$model}{$ref_seq_accn}) { die "ERROR no prediction in reference for feature $mft_i exon $exon_i.\n"; }
    my $start  = $p_start_HH{$model}{$ref_seq_accn};
    my $stop   = ($do_nocorrect) ? $p_stop_HH{$model}{$ref_seq_accn} : $c_stop_HH{$model}{$ref_seq_accn};
    
    my $ref_match = checkStrictBoundaryMatch(\@tmp_ref_act_exon_starts_AA, \@tmp_ref_act_exon_stops_AA, $mft_i, $exon_i, $start, $stop, ($do_nodup ? undef : $ref_totlen));
    if(! $ref_match) { die "ERROR, predicted reference feature $mft_i exon $exon_i does not match GenBank annotation."; }
    
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

#########################################
# OUTPUT ANNOTATION TABLE AND ERROR CODES
#########################################
# open files for writing
my $tblout_file = $out_root . ".tbl";
my $tblout_FH = undef;
open($tblout_FH, ">", $tblout_file) || die "ERROR, unable to open $tblout_file for writing";

my $fail_tblout_file = $out_root . ".fail.tbl";
my $fail_tblout_FH = undef;
open($fail_tblout_FH, ">", $fail_tblout_file) || die "ERROR, unable to open $fail_tblout_file for writing";
my $nfail_accn = 0;

my $err_tblout_file = $out_root . ".error.tbl";
my $err_tblout_FH = undef;
open($err_tblout_FH, ">", $err_tblout_file) || die "ERROR, unable to open $err_tblout_file for writing";
my $nerr_accn = 0;

my $errors_per_accn_file = $out_root . ".peraccn.errors";
my $errpaout_FH = undef;
open($errpaout_FH, ">", $errors_per_accn_file) || die "ERROR, unable to open $errors_per_accn_file for writing";
my $nlines_per_accn_file = 0;

my $errors_per_err_file = $out_root . ".all.errors";
my $errpeout_FH = undef;
open($errpeout_FH, ">", $errors_per_err_file) || die "ERROR, unable to open $errors_per_err_file for writing";
my $nlines_per_err_file = 0;

# output headers to $errpaout_FH file
printf $errpaout_FH ("# Each accession for which at least one error was found is printed below.\n");
printf $errpaout_FH ("# One line per accession. Each line is formatted as follows:\n");
printf $errpaout_FH ("#   <accession> <idxA>:<errorcodeA1>(,<errorcodeAM>) <idxN>:<errorcodeN1>(,<errorcodeNM>)\n");
printf $errpaout_FH ("# For indices (<idx>) A to N, each with M error codes.\n");
printf $errpaout_FH ("# Each index refers to a 'feature' in the reference accession as defined below.\n");
printf $errpaout_FH ("# If no index exists, then the error pertains to the entire sequence.\n");

# output header to $errpeout_FH file
printf $errpeout_FH ("# Each error encountered is printed below, one error per line.\n");
printf $errpeout_FH ("# Each line has four columns with the following labels:\n");
printf $errpeout_FH ("#   \"accn\"         : sequence accession\n");
printf $errpeout_FH ("#   \"idx\"          : feature index, full feature names are listed below for each index\n");
printf $errpeout_FH ("#   \"code\"         : 3 digit error code\n");
printf $errpeout_FH ("#   \"error-message\": error message, possibly with additional information at end enclosed in \"[]\"\n");
printf $errpeout_FH ("#\n");
printf $errpeout_FH ("# List of features:\n");

# list feature names to both $errpaout_FH and $errpeout_FH
my $mft_i = 0;
for($mft_i = 0; $mft_i < $ref_nmft; $mft_i++) { 
  printf $errpaout_FH ("# Feature \#%d: $mft_out_short_A[$mft_i] $mft_out_product_A[$mft_i]\n", ($mft_i+1));
  printf $errpeout_FH ("# Feature \#%d: $mft_out_short_A[$mft_i] $mft_out_product_A[$mft_i]\n", ($mft_i+1));
}
if($do_matpept) { 
  foreach (my $cds_idx = 0; $cds_idx < scalar(@cds2pmatpept_AA); $cds_idx++) { 
    printf $errpaout_FH ("# Feature \#%d: CDS \#%d\n", ($mft_i+1), $cds_idx+1);
    printf $errpeout_FH ("# Feature \#%d: CDS \#%d\n", ($mft_i+1), $cds_idx+1);
    $mft_i++;
  }
}
printf $errpeout_FH ("# \"N/A\" in feature index column (fidx) indicates error pertains to the entire sequence\n");
printf $errpeout_FH ("#       as opposed to a specific feature.\n");
printf $errpeout_FH ("#\n");
my $div_line = "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
if($do_matpept) { 
  outputCDSMaturePeptideRelationships($errpaout_FH, \@cds2pmatpept_AA, \@cds2amatpept_AA, $div_line);
  outputCDSMaturePeptideRelationships($errpeout_FH, \@cds2pmatpept_AA, \@cds2amatpept_AA, $div_line);
}
printf $errpeout_FH ("%-10s  %3s  %-5s  %4s  error-message\n", "#accn", "idx", "desc", "code");
printf $errpaout_FH ("#\n");

#######################

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
            \@mft_out_short_A, \@mft_out_product_A, \@cds2pmatpept_AA, ($do_seqrow) ? (\@out_col_header_AA) : undef,
            ($do_seqcol) ? (\@out_row_header_A)  : undef,
            (\@out_header_exp_A));


if($do_seqrow) { # output sequences as rows 
  for(my $i = 0; $i < 5; $i++) { 
    print $tblout_FH "#";
    print $fail_tblout_FH "#";
    print $err_tblout_FH "#";
    my $ncols = scalar(@{$out_col_header_AA[$i]});
    for(my $j = 0; $j < $ncols; $j++) { 
      print $tblout_FH $out_col_header_AA[$i][$j];
      print $fail_tblout_FH $out_col_header_AA[$i][$j];
      print $err_tblout_FH $out_col_header_AA[$i][$j];
    }
    print $tblout_FH "\n";
    print $fail_tblout_FH "\n";
    print $err_tblout_FH "\n";
  }
}

########################################################################
# Pass through all accessions, and gather and output annotation for each
########################################################################
my @ref_ol_AA     = (); # 2D array that describes the overlaps in the reference, $ref_ol_AA[$i][$j] is '1' if the exons modeled by model $i and $j overlap
my @ref_adj_AA    = (); # 2D array that describes the adjacencies in the reference, $ref_ol_AA[$i][$j] is '1' if the exons modeled by model $i and $j are adjacent,
my @ref_ol_str_A  = (); # array of strings that describes the overlaps in the reference, per feature
my @ref_ajb_str_A = (); # array of strings that describes the 'before' adjacencies in the reference, per feature
my @ref_aja_str_A = (); # array of strings that describes the 'after'  adjacencies in the reference, per feature

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
# analagous data structures for 'failure-only' version of the file
my @fail_page_accn_A = (); # [0..$cur_fail_pagesize]2D array, each element is an array of output tokens for one accession
my @fail_page_out_AA = (); # 2D array, [0..$a..$cur_fail_pagesize-1][0..(ntoks-1)] first dimension is of size $cur_pagesize, 
                           # each element is an array of output tokens for accession $fail_page_accn_A[$a]
my $cur_fail_pagesize = 0; # current number of accessions we have info for in fail_page_out_AA (size of first dimension in fail_page_out_AA)
                           # when this hits $nseqcol, we dump the output
my $nfail_pages = 0;       # number of fail pages output, only used if $do_seqcol
# analagous data structures for 'error-only' version of the file
my @err_page_accn_A = (); # [0..$cur_err_pagesize]2D array, each element is an array of output tokens for one accession
my @err_page_out_AA = (); # 2D array, [0..$a..$cur_err_pagesize-1][0..(ntoks-1)] first dimension is of size $cur_pagesize, 
                          # each element is an array of output tokens for accession $err_page_accn_A[$a]
my $cur_err_pagesize = 0; # current number of accessions we have info for in err_page_out_AA (size of first dimension in err_page_out_AA)
                          # when this hits $nseqcol, we dump the output
my $nerr_pages = 0;       # number of error pages output, only used if $do_seqcol


# data structures for reporting error codes
my @cur_err_ps_A        = (); # per-sequence errors [0..$nerrcodes_ps-1], reset for each accession
my @cur_err_pf_AA       = (); # per-feature  errors [0..$nerrcodes_pf-1][0..$nmft-1], reset for each accession
my @cur_err_extra_ps_A  = (); # per-sequence error extra information to print [0..$nerrcodes_pf-1], reset for each accession
my @cur_err_extra_pf_AA = (); # per-feature  error extra information to print [0..$nerrcodes_pf-1][0..$nmft-1], reset for each accession
for(my $a = 0; $a < $naccn; $a++) { 
  my $accn      = $accn_A[$a];
  my $seq_accn  = $seq_accn_A[$a];
  my @cur_out_A = (); # array of current tokens to print
  my $cur_nerr = 0;

  # reset error information to empty for this accession
  for(my $e = 0; $e < $nerrcodes_ps; $e++) { 
    $cur_err_ps_A[$e]       = 0;
    $cur_err_extra_ps_A[$e] = "";
  }
  my $nmft_for_errors = $ref_nmft;
  if($do_matpept) { 
    $nmft_for_errors += scalar(@cds2pmatpept_AA);
  }
  my $mft_i;
  for($mft_i = 0; $mft_i < $nmft_for_errors; $mft_i++) { 
    for(my $e = 0; $e < $nerrcodes_pf; $e++) { 
      $cur_err_pf_AA[$mft_i][$e]       = 0;
      $cur_err_extra_pf_AA[$mft_i][$e] = "";
    }
  }

  # sanity checks
  if(! exists $totlen_H{$accn}) { die "ERROR accession $accn does not exist in the length file $length_file"; }
  
  # Create the initial portion of the output line, the accession and length
  push(@cur_out_A, sprintf("%-5d  ", ($a+1)));
  push(@cur_out_A, sprintf("%-19s  ", $accn)); 
  push(@cur_out_A, sprintf("%6d ", $totlen_H{$accn}));

  #########################################################
  # Get information on the actual annotation of this genome
  #########################################################
  my @act_exon_starts_AA = (); # [0..$nmft-1][0..$nexons-1] start positions of actual annotations of exons for this accn, $nexons is main-feature (CDS or mat_peptide) specific
  my @act_exon_stops_AA  = (); # [0..$nmft-1][0..$nexons-1] stop  positions of actual annotations of exons for this accn, $nexons is main-feature (CDS or mat_peptide) specific
  my $tot_nexons = 0;

  getActualAnnotations($accn, $totlen_H{$accn}, $do_nodup, $mft_tbl_HHAR, \@act_exon_starts_AA, \@act_exon_stops_AA, \$tot_nexons, \$nmft);

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

    if($oseq_passfail eq "F") { 
      setErrorCode(\@cur_err_ps_A, \@cur_err_extra_ps_A, "ori", $oseq_ct . " occurrences", \%err_ps_code2idx_H, \$cur_nerr);
    }
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
      $start  = $p_start_HH{$model}{$seq_accn};
      $stop   = ($do_nocorrect) ? $p_stop_HH{$model}{$seq_accn}  : $c_stop_HH{$model}{$seq_accn};
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
    else { # 1st matpept does not start at nt 1 (normal case)
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
    my $nexons  = $mft2nexon_map_A[$h];
    
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
      my $stop_corrected_exon  = (exists $did_corr_exon_stop_AH[$h]{$accn})  ? 1 : 0;
      my $start    = $p_start_HH{$model}{$seq_accn}; # we never correct a start
      my $stop     = ($do_nocorrect) ? $p_stop_HH{$model}{$seq_accn} : $c_stop_HH{$model}{$seq_accn};
      my $hangover = $p_hangover_HH{$model}{$seq_accn};
      #printf("HEYAC h: $h accn: $accn stop_corrected_exon: $stop_corrected_exon $start..$stop do_nocorrect: $do_nocorrect\n");

      my ($hang5, $hang3) = split(":", $hangover);
      if($hang5 == 0) { 
        $hang5 = "."; 
      }
      else { 
        $at_least_one_fail = 1;
        setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "bd5", $hang5 . " nt from 5' end", \%err_pf_code2idx_H, \$cur_nerr);
        if($hang5 >  9) { $hang5 = "+"; } 
      }

      if($hang3 == 0) { 
        $hang3 = "."; 
      }
      else { 
        $at_least_one_fail = 1;
        setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "bd3", $hang3 . "nt from 3' end", \%err_pf_code2idx_H, \$cur_nerr);
        if($hang3 >  9) { $hang3 = "+"; } 
      }

      my $trc_extra_errmsg = ($nexons == 1) ? "" : sprintf("exon %d of %d ", $cur_nexons, $nexons);
      my $trc_errmsg = sprintf("homology search predicted $p_start_HH{$model}{$seq_accn}..$p_stop_HH{$model}{$seq_accn} %srevised to $start..$stop (stop shifted %d nt)", $trc_extra_errmsg, abs($stop-$p_stop_HH{$model}{$seq_accn}));
      if($stop_corrected_exon) { 
        setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "trc", $trc_errmsg, \%err_pf_code2idx_H, \$cur_nerr);
      }
      if($stop_corrected_exon) { 
        $hang3 = "c"; 
      }
      #printf("HEYAAA earlier mft_i: $mft_i e: %d ", $err_pf_code2idx_H{"trc"});
      #printf("%d\n", ((defined $cur_err_pf_AA[$mft_i][$err_pf_code2idx_H{"trc"}]) && ($cur_err_pf_AA[$mft_i][$err_pf_code2idx_H{"trc"}])) ? 1 : 0);

      my $boundary_match;
      $boundary_match = ($do_strict) ? 
          checkStrictBoundaryMatch   (\@act_exon_starts_AA, \@act_exon_stops_AA, $mft_i, $exon_i, $start, $stop, ($do_nodup ? undef : $totlen_H{$accn})) :
          checkNonStrictBoundaryMatch(\@act_exon_starts_AA, \@act_exon_stops_AA,                  $start, $stop, ($do_nodup ? undef : $totlen_H{$accn}));
      if($boundary_match) { $nmatch_boundaries++; }
 
      if($do_nobrack) { # set to '1' so brackets are never printed
        $boundary_match = 1;
      }
      
      $hit_length += abs($stop-$start) + 1;
      if(($stop < 0 && $start > 0) || 
         ($stop > 0 && $start < 0)) { 
        # correct for off-by-one induced by the way we use negative indices distance from -1..1 is 1 nt, not 2
        $hit_length -= 1;
      }

      # TODO: MODIFY ANNOTATION FOR EXONS WITHOUT CORRECTED STARTS OR STOPS IN FEATURES WITH
      #       OTHER EXONS THAT HAVE CORRECTED STARTS OR STOPS
      push(@cur_out_A, sprintf("  %8s ", ($boundary_match ? " " . $start . " " : "[" . $start . "]")));
      push(@cur_out_A, sprintf("%8s",    ($boundary_match ? " " . $stop .  " " : "[" . $stop . "]")));
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
        ($ol_pf_char, $ol_str, undef) = compareOverlapsOrAdjacencies(\@cur_name_A, "/", $h, 1, 1, \@ref_ol_AA, \@cur_ol_AA); # 1, 1 says check features before and after $h
        push(@cur_out_A, sprintf(" %10s", $ol_str));
        my $ol_str_errmsg = $ol_str;
        if($ol_pf_char eq "F") { 
          $ol_str_errmsg =~ s/^F://;
          if(! defined $ref_ol_str_A[$mft_i]) { # determine it
            $ref_ol_str_A[$mft_i] = getOverlapsOrAdjacenciesString(\@cur_name_A, $h, 1, 1, \@ref_ol_AA);
          }
          $ol_str_errmsg .= " != " . $ref_ol_str_A[$mft_i] . " (ref)";
        }
        
        if($ol_pf_char eq "F") { 
          $at_least_one_fail = 1;
          setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "olp", $ol_str_errmsg, \%err_pf_code2idx_H, \$cur_nerr);
        }
      }

      # adjacency string for this feature
      if($do_matpept) { 
        my $adj_pf_char  = undef;
        my $adj_str      = undef;
        my $fail_indices = undef;
        ($adj_pf_char, $adj_str) = compareOverlapsOrAdjacencies(\@cur_name_A, "|", $h, 1, 1, \@ref_adj_AA, \@cur_adj_AA);
        push(@cur_out_A, $adj_str);
        if($adj_pf_char eq "F") { 
          # determine whether we have a ajb or aja (adjacent before or adjacent after) error or both
          # first check for ajb error
          my ($ajb_pf_char, $ajb_str_errmsg) = compareOverlapsOrAdjacencies(\@cur_name_A, "|", $h, 1, 0, \@ref_adj_AA, \@cur_adj_AA); #1, 0: only consider features before $h
          if($ajb_pf_char eq 'F') { 
            $ajb_str_errmsg =~ s/^F://;
            if(! defined $ref_ajb_str_A[$mft_i]) { # determine it
              $ref_ajb_str_A[$mft_i] = getOverlapsOrAdjacenciesString(\@cur_name_A, $h, 1, 0, \@ref_adj_AA); # 1, 0: only consider features before $h (< $h)
            }
            $ajb_str_errmsg .= " != " . $ref_ajb_str_A[$mft_i] . " (ref)";
            setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "ajb", $ajb_str_errmsg, \%err_pf_code2idx_H, \$cur_nerr);
          }
          # now check for aja error
          my ($aja_pf_char, $aja_str_errmsg) = compareOverlapsOrAdjacencies(\@cur_name_A, "|", $h, 0, 1, \@ref_adj_AA, \@cur_adj_AA); #0, 1: only consider features after $h
          if($aja_pf_char eq 'F') { 
            $aja_str_errmsg =~ s/^F://;
            if(! defined $ref_aja_str_A[$mft_i]) { # determine it
              $ref_aja_str_A[$mft_i] = getOverlapsOrAdjacenciesString(\@cur_name_A, $h, 0, 1, \@ref_adj_AA); # 0, 1: only consider features after $h (< $h)
            }
            $aja_str_errmsg .= " != " . $ref_aja_str_A[$mft_i] . " (ref)";
            setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "aja", $aja_str_errmsg, \%err_pf_code2idx_H, \$cur_nerr);
          }
          $at_least_one_fail = 1;
        }
      }

      if($mdl_is_first_A[$h] && (! $do_matpept)) { # determine $start_codon_char
        $start_codon = fetchStartCodon($sqfile, $seq_accn, $start, $totlen_H{$accn}, $p_strand_HH{$model}{$seq_accn});
        $start_codon_char = ($start_codon eq "ATG") ? $ss3_yes_char : $ss3_no_char;
        if($start_codon_char eq $ss3_no_char) { 
          setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "str", $start_codon, \%err_pf_code2idx_H, \$cur_nerr);
          $at_least_one_fail = 1;
        }
      }

      if($mdl_is_final_A[$h]) { 
        if(! $do_matpept) { 
          $stop_codon = fetchStopCodon($sqfile, $seq_accn, $stop, $totlen_H{$accn}, $p_strand_HH{$model}{$seq_accn});
          $stop_codon_char = ($stop_codon eq "TAG" || $stop_codon eq "TAA" || $stop_codon eq "TGA") ? $ss3_yes_char : $ss3_no_char;
          if($stop_codon_char eq $ss3_no_char) { 
            $at_least_one_fail = 1;
            setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "stp", $stop_codon, \%err_pf_code2idx_H, \$cur_nerr);
          }
        }

        $multiple_of_3_char = (($hit_length % 3) == 0) ? $ss3_yes_char : $ss3_no_char;
        if($multiple_of_3_char eq $ss3_no_char) { 
          setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "nm3", $hit_length, \%err_pf_code2idx_H, \$cur_nerr);
          $at_least_one_fail = 1;
        }
        push(@cur_out_A, sprintf(" %6d", $hit_length));
        # TODO: put in length error check here!

        # add the ss3 (start/stop/multiple of 3 info) if we're not in mat_pept mode
        if((! $do_matpept) && (! $do_noss3)) { 
          push(@cur_out_A,  sprintf(" %s%s%s", $start_codon_char, $stop_codon_char, $multiple_of_3_char));
        }
        # add the stop codon if we're not in mat_pept mode
        if((! $do_matpept) && (! $do_nostop)) { 
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
      setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "nop", undef, \%err_pf_code2idx_H, \$cur_nerr);
      push(@cur_out_A, sprintf("  %8s ", "NP")); # start position
      push(@cur_out_A, sprintf("%8s",  "NP"));   # stop position
      if(! $do_nofid) { 
        push(@cur_out_A, sprintf(" %5s", "NP")); # fid
      }        
      if(! $do_nomdlb) { 
        push(@cur_out_A, "  NP"); # model boundaries
      }
      if(! $do_noolap) { 
        push(@cur_out_A, "  NP"); # overlaps
      }
      if($do_matpept) { 
        push(@cur_out_A, "  NP"); # adjacencies
      }
      if($mdl_is_final_A[$h]) { 
        push(@cur_out_A, sprintf(" %6s", "NP")); # length
        if((! $do_matpept) && (! $do_noss3)) { 
          push(@cur_out_A, "  NP"); # ss3
        }
        if((! $do_matpept) && (! $do_nostop)) { 
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

  # do CDS for matpept mode
  if($do_matpept) { 
    my $start_codon_char;   # one of the 3 characters in the ss3 field, indicating if the start codon is valid or not
    my $stop_codon_char;    # one of the 3 characters in the ss3 field, indicating if the stop  codon is valid or not
    my $multiple_of_3_char; # one of the 3 characters in the ss3 field, indicating if the length is a multiple of 3 or not
    my $mft_i = $ref_nmft-1;
    foreach (my $cds_idx = 0; $cds_idx < scalar(@cds2pmatpept_AA); $cds_idx++) { 
      $mft_i++;
      my ($cds_start, $cds_stop, $cds_corr_stop, $cds_corr_stop_mdl_idx, $cds_len, $cds_start_codon, $cds_stop_codon, $cds_pass_fail) = matpeptCheckCdsRelationships($sqfile, $seq_accn, $totlen_H{$accn}, \@mdl_A, \@{$cds2pmatpept_AA[$cds_idx]}, \%p_start_HH, \%p_stop_HH, \%c_stop_HH, \%p_strand_HH, \@mft2first_mdl_A, \@mft2final_mdl_A);
      my $orig_cds_stop = $cds_stop;

      # check if we have a truncated stop
      if(defined $cds_corr_stop) { 
        $cds_stop = $cds_corr_stop;
        $cds_len  = abs($cds_stop - $cds_start) + 1;
        # set 'trc' error
        my $trc_errmsg = sprintf("homology search predicted $cds_start..$orig_cds_stop revised to $cds_start..$cds_stop (stop shifted %d nt)", abs($cds_stop-$orig_cds_stop));
        setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "trc", $trc_errmsg, \%err_pf_code2idx_H, \$cur_nerr);

        # and set 'ntr' (not translated) errors for all mat peptides that are not translated due to this truncation:
        # first determine mft index (mat_peptide index) that $cds_corr_stop_mdl_idx corresponds to:
        my $cds_corr_stop_mft_idx = $mdl2mft_map_A[$cds_corr_stop_mdl_idx];
        my $ntr_errmsg = sprintf("early stop in mature peptide %d ending at position $cds_stop", $cds_corr_stop_mft_idx+1);
        for(my $i = 0; $i < scalar(@{$cds2amatpept_AA[$cds_idx]}); $i++) { # step through ALL peptides encoded by this CDS, not just primary ones
          my $mft_j = $cds2amatpept_AA[$cds_idx][$i];
          if($mft_j > $cds_corr_stop_mft_idx) { 
            setErrorCode(\@{$cur_err_pf_AA[$mft_j]}, \@{$cur_err_extra_pf_AA[$mft_j]}, "ntr", $ntr_errmsg, \%err_pf_code2idx_H, \$cur_nerr);
          }
        }
      }

      if(defined $cds_start) { 
        push(@cur_out_A, sprintf("  %6d", $cds_start)); 
      } 
      else { 
        push(@cur_out_A, sprintf("  %6s", "?")); 
      }
      if(defined $cds_stop) {
        push(@cur_out_A, sprintf("  %6d", $cds_stop)); 
      } 
      else { 
        push(@cur_out_A, sprintf("  %6s", "?")); 
      }

      if(defined $cds_len)   { 
        push(@cur_out_A, sprintf("  %6d", $cds_len)); 
        $multiple_of_3_char = (($cds_len % 3) == 0) ? $ss3_yes_char : $ss3_no_char;
      }
      else {
        push(@cur_out_A, sprintf("  %6s", "?"));
        $multiple_of_3_char = $ss3_no_char;
      }
      if($multiple_of_3_char ne $ss3_yes_char) { 
        setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "nm3", (defined $cds_len) ? $cds_len : "?", \%err_pf_code2idx_H, \$cur_nerr);
      }

      if(defined $cds_start_codon) { 
        push(@cur_out_A, sprintf("  %6s", $cds_start_codon)); 
        $start_codon_char = ($cds_start_codon eq "ATG") ? $ss3_yes_char : $ss3_no_char;
      } 
      else { 
        push(@cur_out_A, sprintf("  %6s", "?"));
        $start_codon_char = $ss3_unsure_char;
      }
      if($start_codon_char ne $ss3_yes_char) { 
        setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "str", ((defined $cds_start_codon) ? $cds_start_codon : "?"), \%err_pf_code2idx_H, \$cur_nerr);
      }

      if(defined $cds_stop_codon) { 
        if(! $do_nostop) { 
          push(@cur_out_A, sprintf("  %6s", $cds_stop_codon)); 
        }
        $stop_codon_char = ($cds_stop_codon eq "TAG" || $cds_stop_codon eq "TAA" || $cds_stop_codon eq "TGA") ? $ss3_yes_char : $ss3_no_char;
      } 
      else { 
        if(! $do_nostop) { 
          push(@cur_out_A, sprintf("  %6s", "?"));
        }
        $stop_codon_char = $ss3_unsure_char;
      }
      if($stop_codon_char ne $ss3_yes_char) { 
        setErrorCode(\@{$cur_err_pf_AA[$mft_i]}, \@{$cur_err_extra_pf_AA[$mft_i]}, "stp", ((defined $cds_stop_codon) ? $cds_stop_codon : "?"), \%err_pf_code2idx_H, \$cur_nerr);
      }

      if(! $do_noss3) { 
        push(@cur_out_A,  sprintf(" %s%s%s", $start_codon_char, $stop_codon_char, $multiple_of_3_char));
      }
      if($start_codon_char ne $ss3_yes_char || $stop_codon_char ne $ss3_yes_char || $multiple_of_3_char ne $ss3_yes_char) { 
        $cds_pass_fail = "F";
      }
      push(@cur_out_A, sprintf("  %3s", $cds_pass_fail)); 
      $pass_fail_str .= $cds_pass_fail;
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

  my $accn_failed = ($pass_fail_str =~ m/F/) ? 1 : 0;
  if($accn_failed) { $nfail_accn++; }
  my $result_str = ($accn_failed) ? "FAIL" : "PASS";
  $result_str .= " " . $pass_fail_str;
  push(@cur_out_A, sprintf("  %s", $result_str));

  # ACTUALLY OUTPUT THE INFORMATION TO RELEVANT FILE HANDLES
  if($do_seqrow) { 
    foreach my $el (@cur_out_A) { 
      print $tblout_FH $el;
    }
    print $tblout_FH "\n";
    if($accn_failed || $a == 0) { 
      foreach my $el (@cur_out_A) { 
        print $fail_tblout_FH $el;
      }
      print $fail_tblout_FH "\n";
    }
    if($cur_nerr > 0) { 
      foreach my $el (@cur_out_A) { 
        print $err_tblout_FH $el;
      }
      print $err_tblout_FH "\n";
    } 
  }
  elsif($do_seqcol) { 
    if($a == 0) { 
      # copy reference info if this is the reference
      @ref_out_A = @cur_out_A; 
    } 
    else { 
      push(@page_out_AA, [@cur_out_A]);
      $cur_pagesize++;
      if($accn_failed) { 
        push(@fail_page_out_AA, [@cur_out_A]);
        $cur_fail_pagesize++;
      }        
      if($cur_nerr > 0) { 
        push(@err_page_out_AA, [@cur_out_A]);
        $cur_err_pagesize++;
      }        
    }
    if(($cur_pagesize+1) == $nseqcol) { 
      $npages++;
      outputSeqAsColumnsPage($tblout_FH, \@out_row_header_A, \@page_out_AA, \@ref_out_A, $npages);
      @page_out_AA = ();
      $cur_pagesize = 0;
    }
    if($accn_failed && 
       (($cur_fail_pagesize+1) == $nseqcol)) { 
      $nfail_pages++;
      outputSeqAsColumnsPage($fail_tblout_FH, \@out_row_header_A, \@fail_page_out_AA, \@ref_out_A, $nfail_pages);
      @fail_page_out_AA = ();
      $cur_fail_pagesize = 0;
    }
    if(($cur_nerr > 0) && 
       (($cur_err_pagesize+1) == $nseqcol)) { 
      $nerr_pages++;
      outputSeqAsColumnsPage($err_tblout_FH, \@out_row_header_A, \@err_page_out_AA, \@ref_out_A, $nerr_pages);
      @err_page_out_AA = ();
      $cur_err_pagesize = 0;
    }
  }

  # OUTPUT to error file
  if($cur_nerr > 0) { 
    print $errpaout_FH $accn;
    # first print per-sequence errors, if any
    for($e = 0; $e < $nerrcodes_ps; $e++) { 
      if((defined $cur_err_ps_A[$e]) && ($cur_err_ps_A[$e])) { 
        print $errpaout_FH " " . $err_ps_idx2code_A[$e];
        printf $errpeout_FH ("%-10s  %3s  %4s  %s%s\n", $accn, "N/A", $err_ps_idx2code_A[$e], $err_ps_idx2msg_H{$err_ps_idx2code_A[$e]}, 
                        (defined $cur_err_extra_ps_A[$e]) ? " [" . $cur_err_extra_ps_A[$e]. "]" : "");
        $nlines_per_err_file++;
      }
    }        
    for($mft_i = 0; $mft_i < $nmft_for_errors; $mft_i++) { 
      my $nerr_printed = 0;
      for($e = 0; $e < $nerrcodes_pf; $e++) { 
        #printf("HEYAAA mft_i: $mft_i e: $e ");
        #printf("%d\n", ((defined $cur_err_pf_AA[$mft_i][$e]) && ($cur_err_pf_AA[$mft_i][$e])) ? 1 : 0);
        if((defined $cur_err_pf_AA[$mft_i][$e]) && ($cur_err_pf_AA[$mft_i][$e])) { 
          if($nerr_printed == 0) { 
            #print $errpaout_FH (" feature\#" . ($mft_i+1) . ":" . $err_pf_idx2code_A[$e]); 
            print $errpaout_FH (" " . ($mft_i+1) . ":" . $err_pf_idx2code_A[$e]); 
          }
          else { 
            print $errpaout_FH "," . $err_pf_idx2code_A[$e];
          }
          $nerr_printed++;
          my $desc = ($mft_i < $ref_nmft) ? $mft_out_tiny_A[$mft_i] : sprintf("CDS#%d", ($mft_i - $ref_nmft + 1));
          printf $errpeout_FH ("%-10s  %3s  %-5s  %4s  %s%s\n", $accn, ($mft_i+1), $desc, $err_pf_idx2code_A[$e], $err_pf_idx2msg_H{$err_pf_idx2code_A[$e]}, 
                               (defined $cur_err_extra_pf_AA[$mft_i][$e]) ? " [" . $cur_err_extra_pf_AA[$mft_i][$e]. "]" : "");
          $nlines_per_err_file++;
        }
      }
    }
    print $errpaout_FH "\n";
    $nlines_per_accn_file++;
  }
}
# print final page (if non-empty)
if(($do_seqcol) && ($cur_pagesize > 0)) { 
  $npages++;
  outputSeqAsColumnsPage($tblout_FH, \@out_row_header_A, \@page_out_AA, \@ref_out_A, $npages);
}
if(($do_seqcol) && ($cur_fail_pagesize > 0)) { 
  $nfail_pages++;
  outputSeqAsColumnsPage($fail_tblout_FH, \@out_row_header_A, \@fail_page_out_AA, \@ref_out_A, $nfail_pages);
}
if(($do_seqcol) && ($cur_err_pagesize > 0)) { 
  $nerr_pages++;
  outputSeqAsColumnsPage($err_tblout_FH, \@out_row_header_A, \@err_page_out_AA, \@ref_out_A, $nerr_pages);
}

print $tblout_FH $div_line;
print $fail_tblout_FH $div_line;
print $err_tblout_FH $div_line;

##############################################
# OUTPUT CDS:MAT_PEPTIDE RELATIONSHIPS, IF NEC 
##############################################
if($do_matpept) { 
  outputCDSMaturePeptideRelationships($tblout_FH,      \@cds2pmatpept_AA, \@cds2amatpept_AA, $div_line);
  outputCDSMaturePeptideRelationships($fail_tblout_FH, \@cds2pmatpept_AA, \@cds2amatpept_AA, $div_line);
  outputCDSMaturePeptideRelationships($err_tblout_FH,  \@cds2pmatpept_AA, \@cds2amatpept_AA, $div_line);
}

##########################
# OUTPUT EXPLANATORY TEXT 
##########################
if(! $do_noexp) { 
    outputColumnHeaderExplanations($tblout_FH, \@out_header_exp_A);
    outputColumnHeaderExplanations($fail_tblout_FH, \@out_header_exp_A);
    outputColumnHeaderExplanations($err_tblout_FH, \@out_header_exp_A);
    print $tblout_FH "#\n";
    print $tblout_FH $div_line;
    print $tblout_FH "#\n";
    print $fail_tblout_FH "#\n";
    print $fail_tblout_FH $div_line;
    print $fail_tblout_FH "#\n";
    print $err_tblout_FH "#\n";
    print $err_tblout_FH $div_line;
    print $err_tblout_FH "#\n";
}

# CLOSE TABULAR OUTPUT FILES
my $outfile_w = 85;
close($tblout_FH);
close($fail_tblout_FH);
close($err_tblout_FH);
printf ("%-*s %s\n", $outfile_w, "# Saved tabular annotation information on all $naccn accessions to file.", "[$tblout_file]");
printf ("%-*s %s\n", $outfile_w, "# Saved tabular annotation information on $nlines_per_accn_file accessions with >= 1 error to file.", "[$err_tblout_file]");
printf ("%-*s %s\n", $outfile_w, "# Saved tabular annotation information on $nfail_accn accessions that FAILed to file.", "[$fail_tblout_file]");

# CLOSE ERROR FILES
close($errpaout_FH);
close($errpeout_FH);
printf ("%-*s %s\n", $outfile_w, "# Saved information on $nlines_per_accn_file accessions with at least one error to file.", "[$errors_per_accn_file]");
printf ("%-*s %s\n", $outfile_w, "# Saved information on $nlines_per_err_file errors to file.", "[$errors_per_err_file]");

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

($seconds, $microseconds) = gettimeofday();
my $end_secs = ($seconds + ($microseconds / 1000000.));
printf("#\n");
printf("# Total time: %.1f seconds\n", ($end_secs - $start_secs));

print ("#[ok]\n");

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
# Args:       $tbl_HAR:   ref to hash of arrays for accession we're interested in
#             $len_AR:    ref to array to fill with lengths of features in %{$tbl_HAR}
#             $coords_AR: ref to array to fill with coordinates for each gene
# Returns:    void; fills @{$len_AR} and @{$coords_AR}
#
sub getLengthStatsAndCoordStrings { 
  my $sub_name = "getLengthStatsAndCoordStrings()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($tbl_HAR, $len_AR, $coords_AR) = @_;

  my $ngenes = scalar(@{$tbl_HAR->{"coords"}});

  if ($ngenes > 0) { 
    for(my $i = 0; $i < $ngenes; $i++) { 
      push(@{$len_AR},    lengthFromCoords($tbl_HAR->{"coords"}[$i]));
      push(@{$coords_AR}, $tbl_HAR->{"coords"}[$i]);
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
# Args:       $coords:     the coords string
#             $totlen:     total length of sequence, can be undef if we don't want to call checkForSpanningExon()
#             $do_nodup:   '1' if -nodup option enabled, can be undef if we don't want to call checkForSpanningExon()
#             $starts_AR:  ref to array to fill with start positions
#             $stops_AR:   ref to array to fill with stop positions
#             $strands_AR: ref to array to fill with strands of each exon, can be undef
#             $nexons_R:   ref to scalar that fill with the number of exons
#
# Returns:    void; but fills @{$starts_AR}, @{$stops_AR}, and $$nexons_R.
#
sub startStopsFromCoords { 
  my $sub_name = "startStopsFromCoords()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($coords, $totlen, $do_nodup, $starts_AR, $stops_AR, $strands_AR, $nexons_R) = @_;

  @{$starts_AR} = ();
  @{$stops_AR}  = ();
  $$nexons_R    = 0;
  
  my $orig_coords = $coords;
  # Examples:
  # complement(2173412..2176090)
  # complement(join(226623..226774, 226854..229725))

  # remove 'complement('  ')'
  my $strand = "+";
  if($coords =~ m/^complement\(/) { 
    $coords =~ s/^complement\(//;
    $strand = "-";
  }
  $coords =~ s/\)$//;

  # remove 'join('  ')'
  $coords =~ s/^join\(//;
  $coords =~ s/\)$//;

  my @el_A = split(/\s*\,\s*/, $coords);

  my $length = 0;
  my $cur_strand = $strand;
  foreach my $el (@el_A) { 
    # rare case: remove 'complement(' ')' that still exists:
    $cur_strand = $strand;
    if($el =~ m/^complement\(/) { 
      die "ERROR in $sub_name: found internal complement in coords string $coords, we assume all exons are on same strand..."; 
      $el =~ s/^complement\(//;
      if($cur_strand eq "-") { die "ERROR in $sub_name, found nested 'complement' annotations in coord string: $coords"; }
      $cur_strand = "-";
    }
    $el =~ s/\)$//;
    $el =~ s/\<//; # remove '<'
    $el =~ s/\>//; # remove '>'
    if($el =~ m/^(\d+)\.\.(\d+)$/) { 
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $2);
      if(defined $strands_AR) { push(@{$strands_AR}, $cur_strand); }
      $$nexons_R++;
    }
    else { 
      die "ERROR unable to parse $orig_coords in $sub_name"; 
    }
  }

  # if we're in a circular genome, and we've been passed in values for $do_nodup
  # and $totlen, then we want to check for a special case, where 
  # what looks like a 2-exon CDS is really a single exon that spans the stop..start boundary.
  if(defined $do_nodup && defined $totlen) { 
    if((! $do_nodup) && ($$nexons_R == 2)) { 
      checkForSpanningExon($starts_AR, $stops_AR, $nexons_R, $strand, $totlen);
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

  # note we pass in undef for $do_nodup and for $totlen because we don't need to check for 
  # exons that span stop..start in this case
  startStopsFromCoords($coords, undef, undef, \@starts_A, \@stops_A, undef, \$nexons);

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
#    if(! -s $stdout_file) { die "ERROR, with -skipscan or -skipaln nhmmscan output is expected to already exist, but $stdout_file does not or is empty"; }
    printf("done. [-skipscan or -skipaln]\n");
    return;
  }

  my $opts = " --noali --tblout $tblout_file ";

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
#    if(! -s $stdout_file) { die "ERROR, with -skipscan or -skipaln cmscan output is expected to already exist, but $stdout_file does not or is empty"; }
    printf("done. [-skipscan or -skipaln]\n");
    return;
  }

  my $opts = "";
  if($do_iglocal) { $opts .= "-g "; }
  # $opts .= " --cpu 0 --rfam --tblout $tblout_file --verbose --nohmmonly ";
  # opts: nhmmer F1, F2, F2b, F3 and F3b; Infernal --rfam F4, F4b, F5, and F6
  $opts .= " --noali --cpu 0 --F1 0.02 --F2 0.001 --F2b 0.001 --F3 0.00001 --F3b 0.00001 --F4 0.0002 --F4b 0.0002 --F5 0.0002 --F6 0.0001 --tblout $tblout_file --verbose --nohmmonly ";
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
# Synopsis:   Identify all exact occurrences of a sequence in a file
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
    # if $qseq_posn == -1, then no occurrences were found. In this case we don't store 
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
#             $totlen:        total length of (non-duplicated) sequence, undef if we don't want to 
#                             check if spanning match exists. (Spanning match is one predicted coords
#                             are -x..y, and actual coords are -x+L+1,y+L where l is total length
#                             of the sequence. This is the same coordinates but one is offset by L,
#                             which can occur in duplicated genomes.)
# 
# Returns:    '1' if $pstart == $act_start_AAR->[$mft_i][$exon_i] and 
#                    $pstop  == $act_stop_AAR->[$mft_i][$exon_i], else '0'
#
# 
sub checkStrictBoundaryMatch {
  my $sub_name = "checkStrictBoundaryMatch";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($act_start_AAR, $act_stop_AAR, $mft_i, $exon_i, $pstart, $pstop, $totlen) = @_;

  my $start1 = ((exists $act_start_AAR->[$mft_i]) && 
                (exists $act_start_AAR->[$mft_i][$exon_i])) ? 
                $act_start_AAR->[$mft_i][$exon_i] : undef;
  my $stop1  = ((exists $act_stop_AAR->[$mft_i]) && 
                (exists $act_stop_AAR->[$mft_i][$exon_i])) ? 
                $act_stop_AAR->[$mft_i][$exon_i] : undef;
  if(defined $start1 && (! defined $stop1)) { die "ERROR in $sub_name, start1 is defined but stop1 is not."; }
  my $start2 = undef;
  my $stop2  = undef;
  if(defined $totlen && defined $start1 && defined $stop1) { 
    $start2 = $start1 - $totlen;
    $stop2  = $stop1  - $totlen;
    if($start2 < 0) { $start2--; } # correct for off-by-one with negative coords (b/c '0' is not a valid posn)
    if($stop2  < 0) { $stop2--;  } # correct for off-by-one with negative coords (b/c '0' is not a valid posn)
  }

  #printf("in $sub_name\n");
  #printf("predicted start:   $pstart\n");
  #printf("predicted stop:    $pstop\n");
  #if(defined $start1) { printf("actual    start1:  $start1\n"); }
  #if(defined $stop1)  { printf("actual    stop1:   $stop1\n"); }
  #if(defined $start2) { printf("actual    start2:  $start2\n"); }
  #if(defined $stop2)  { printf("actual    stop2:   $stop2\n"); }

  if(defined $start1 && $start1 == $pstart && defined $stop1 && $stop1 == $pstop) { 
    return 1;
  }
  if(defined $start2 && $start2 == $pstart && defined $stop2 && $stop2 == $pstop) { 
    return 1;
  }
  return 0;
}

# Subroutine: checkNonStrictBoundaryMatch
#
# Synopsis:   Check if a given start..stop boundary pair matches any
#             annotation in the 2D array referred to by $act_AAR.
#
# Args:       $act_start_AAR: ref to 2D array [0..i..$nmft-1][0..e..$nexon-1], start for feature $i+1 exon $e+1
#     :       $act_stop_AAR:  ref to 2D array [0..i..$nmft-1][0..e..$nexon-1], stop for feature $i+1 exon $e+1
#             $pstart:        predicted start position
#             $pstop:         predicted stop position
#             $totlen:        total length of (non-duplicated) sequence, undef if we don't want to 
#                             check if spanning match exists. (Spanning match is one predicted coords
#                             are -x..y, and actual coords are -x+L+1,y+L where l is total length
#                             of the sequence. This is the same coordinates but one is offset by L,
#                             which can occur in duplicated genomes.)
#
# Returns:    Two values:
#             '1' if $pstart == $act_start_AAR->[$i][$e] and 
#                    $pstop  == $act_stop_AAR->[$i][$e] for any pair $i and $e
#             else '0'
#             That is both start and stop boundaries have to match the same start and stop.
#
sub checkNonStrictBoundaryMatch {
  my $sub_name = "checkNonStrictBoundaryMatch";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($act_start_AAR, $act_stop_AAR, $pstart, $pstop, $totlen) = @_;

  my $nmft = scalar(@{$act_start_AAR});
  for(my $i = 0; $i < $nmft; $i++) { 
    my $nexons = scalar(@{$act_start_AAR->[$i]});
    if(! exists $act_stop_AAR->[$i]) { die "ERROR in checkNonStrictBoundaryMatch() $i exists in first dimension of start coords, but not stop coords"; }
    for(my $e = 0; $e < $nexons; $e++) { 
      if(! exists $act_stop_AAR->[$i][$e]) { die "ERROR in checkNonStrictBoundaryMatch() $i $e exists in start coords, but not stop coords"; }
      if(checkStrictBoundaryMatch($act_start_AAR, $act_stop_AAR, $i, $e, $pstart, $pstop, $totlen)) {
        return 1;
      }
    }
  }

  # if we get here, we didn't find a match
  return 0;
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

  # printf("\nin $sub_name, $seqname $start $strand returning $seq\n");
  
  return $seq;
}

# Subroutine: fetchStartCodon()
#
# Synopsis:   Fetch a start codon given it's first position,
#             the strand and the total length of the sequence.
#             We need the total length if it's on the reverse 
#             strand.
#
# Args:       $sqfile:  Bio::Easel::SqFile object, open sequence
#                       file containing $seqname;
#             $seqname: name of sequence to fetch part of
#             $start:   start position of the start codon
#             $L:       total length of sequence
#             $strand:  strand we want ("+" or "-")
#
# Returns:    The start codon as a string
#
sub fetchStartCodon {
  my $sub_name = "fetchStartCodon";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seqname, $start, $L, $strand) = @_;

  my $start_codon_posn;
  if($strand eq "-") { 
    $start_codon_posn = (($start-2) < 0) ? $start + $L + 1 : $start;
  }
  else { 
    $start_codon_posn = ($start < 0) ? $start + $L + 1 : $start;
  }
  return fetchCodon($sqfile, $seqname, $start_codon_posn, $strand);
}

# Subroutine: fetchStopCodon()
#
# Synopsis:   Fetch a stop codon given it's final position,
#             the strand and the total length of the sequence.
#             We need the total length if it's on the reverse 
#             strand.
#
# Args:       $sqfile:  Bio::Easel::SqFile object, open sequence
#                       file containing $seqname;
#             $seqname: name of sequence to fetch part of
#             $stop:    final position of the stop codon
#             $L:       total length of sequence
#             $strand:  strand we want ("+" or "-")
#
# Returns:    The start codon as a string
#
sub fetchStopCodon {
  my $sub_name = "fetchStopCodon";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seqname, $stop, $L, $strand) = @_;

  my $stop_codon_posn;
  if($strand eq "-") { 
    $stop_codon_posn    = ($stop < 0) ? ($stop + $L) + 1 + 2 : $stop + 2;
  }
  else { 
    $stop_codon_posn = (($stop-2) <= 0) ? ($stop + $L) + 1 - 2 : $stop - 2;
  }
  # printf("in $sub_name, seqname $seqname, stop $stop\n");

  return fetchCodon($sqfile, $seqname, $stop_codon_posn, $strand);
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
#                              those two features are observed to overlap or be adjacecent
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
#             $do_before:    only relevant if $index is defined, if '1' check features that
#                            come before $index (< $index), if '0' do not check features < $index
#             $do_after:     only relevant if $index is defined, if '1' check features that
#                            come after $index (> $index), if '0' do not check features > $index
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
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($name_AR, $div_char, $index, $do_before, $do_after, $expected_AAR, $observed_AAR) = @_;

  if(defined $index && (! defined $do_before)) { die "ERROR in $sub_name() index is defined but do_before is not"; }
  if(defined $index && (! defined $do_after))  { die "ERROR in $sub_name() index is defined but do_before is not"; }
  if(defined $index && (! $do_before) && (! $do_after)) { die "ERROR in $sub_name, index is defined and both do_before and do_after are 0"; }

  my $size = scalar(@{$name_AR});
  my $ret_str = "";
  my $nfound = 0;
  my $first_j = undef; # first j index to check for an i, differs depending on values of $index, $do_before and $do_after
  my $final_j = undef; # final j index to check for an i, differs depending on values of $index, $do_before and $do_after
  # check @observed_AA against @{$expected_AAR}

  my $pass_fail_char = "P";
  for(my $i = 0; $i < $size; $i++) { 
    if((! defined $index) || ($index == $i)) { 
      if(! defined $index) { 
        $first_j = $i+1;
        $final_j = $size-1;
      }
      elsif($do_before && $do_after) { 
        $first_j = 0;
        $final_j = $size-1;
      }
      elsif($do_before) { 
        $first_j = 0;
        $final_j = $i-1;
      }
      elsif($do_after) { 
        $first_j = $i+1;
        $final_j = $size-1;
      }        
      else { 
        # checked above that at least one of $do_before or $do_after is true
        die "ERROR in $sub_name() unforeseen case of index, do_before and do_after";
      }
      for(my $j = $first_j; $j <= $final_j; $j++) { 
        if($observed_AAR->[$i][$j] ne $expected_AAR->[$i][$j]) { 
          $pass_fail_char = "F";
        }
        if($observed_AAR->[$i][$j]) { 
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

# Subroutine: getOverlapsOrAdjacenciesString()
#
# Synopsis:   Determine an overlap or adjacency string given 
#             a 2D array that describes it.
#
# Args:       $name_AR:      ref to array of short names for each annotation
#             $index:        the feature $index we check overlaps/adjacencies for
#             $do_before:    '1' to include feature indices in return string that are < $index
#                            '0' to not
#             $do_after:     '1' to include feature indices in return string that are > $index
#                            '0' to not
#             $ol_AAR:       ref to 2D array of overlaps $ol_AAR->[$i][$j] is '1' if 
#                            those two exons are expected to overlap or be adjacent, PRE-FILLED
#             
# Returns:    $ol_str:       string of overlap features, separated by commas
#
sub getOverlapsOrAdjacenciesString{
  my $sub_name = "getOverlapsOrAdjacencyString";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($name_AR, $index, $do_before, $do_after, $ol_AAR) = @_;

  if((! $do_before) && (! $do_after)) { die "ERROR in $sub_name, both do_before and do_after are 0"; }

  my $size = scalar(@{$name_AR});
  my $ret_str = "";
  my $nfound = 0;
  my $first_j = undef; # first j index to check for an i, differs depending on values of $index, $do_before and $do_after
  my $final_j = undef; # final j index to check for an i, differs depending on values of $index, $do_before and $do_after

  # check @observed_AA against @{$expected_AAR}
  my $pass_fail_char = "P";
  for(my $i = 0; $i < $size; $i++) { 
    if($index == $i) { 
      if($do_before && $do_after) { 
        $first_j = 0;
        $final_j = $size-1;
      }
      elsif($do_before) { 
        $first_j = 0;
        $final_j = $i-1;
      }
      elsif($do_after) { 
        $first_j = $i+1;
        $final_j = $size-1;
      }        
      # checked above that at least one of $do_before or $do_after is true

      for(my $j = $first_j; $j <= $final_j; $j++) { 
        if($ol_AAR->[$i][$j]) { 
          if($ret_str ne "") { 
            $ret_str .= ",";
          }
          $ret_str .= sprintf("%s", $name_AR->[$j]); 
        }
      }
    }
  }
  if($ret_str eq "") { $ret_str = "NONE"; }

  return $ret_str;
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

    if(! defined $len) { printf("gapstr: $gapstr, el: $el\n"); }

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
# Args:       $dir:             directory for output files
#             $key:             string for naming output files (e.g.: "predicted" or "corrected")
#             $mdl_AR:          ref to array of model names 
#             $acc_order_AR:    ref to array with desired order of accessions
#             $mdl2mft_map_AR:  ref to array mapping models to CDS
#             $mft2mdl_map_AAR: ref to array mapping features to models 
#                               [0..$c..$ref_nmft-1][0..$e..$nexons_for_this_feature-1]: value is index of mdl that encodes exon $e for feature $c
#             $outfile_AR:      ref to array of output files, filled here
#
# Returns:    void
# Dies:       if something unexpected happens when reading the exon fasta files
#
sub wrapperCombineExonsIntoCDS {
  my $sub_name = "combineExonsIntoCDS";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($dir, $key, $mdl_AR, $acc_order_AR,$mdl2mft_map_AR, $mft2mdl_map_AAR, $outfile_AR) = @_;

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));
  printf("%-65s ... ", "# Combining multiple exon $key CDS ");

  my @tmp_exon_fafile_A = (); # temporary array of exon fafiles for all exons in current CDS
  my $nmft = scalar(@{$mft2mdl_map_AAR}); 

  for(my $f = 0; $f < $nmft; $f++) { 
    my $nmdl = scalar(@{$mft2mdl_map_AAR->[$f]});
    my $out_fafile = $dir . "/" . $mdl_AR->[$mft2mdl_map_AAR->[$f][0]] . ".fa";
    $out_fafile  =~ s/ref/$key/;
    if($nmdl == 1) { 
      # a single exon gene, we should already have the sequence from alignHits
      if(! -s $out_fafile) { die sprintf("ERROR, expected output fasta file for CDS %s does not exist: $out_fafile", $mft_i+1); }
    }
    else { 
      # more than one exon/mat_peptide need to be combined to make this CDS/mat_peptide
      # if we're combining exons to make a CDS next line will remove exon substr
      $out_fafile =~ s/\.exon\.\d+//; 
      # if we're combining mature peptides to make a CDS next line will remove mp substr and replace with cds substr
      my $f2print = $f+1;
      $out_fafile =~ s/\.mp\.\d+/\.cds$f2print/; 
      my @tmp_exon_fafile_A = ();
      for(my $h = 0; $h < $nmdl; $h++) { 
        my $mdl_idx = $mft2mdl_map_AAR->[$f][$h];
        my $cur_fafile = $dir . "/" . $mdl_AR->[$mdl_idx] . ".fa";
        $cur_fafile =~ s/ref/$key/;
        push(@tmp_exon_fafile_A, $cur_fafile);
      }
      combineExonsIntoCDS(\@tmp_exon_fafile_A, $acc_order_AR, $out_fafile);
    }
    push(@{$outfile_AR}, $out_fafile);
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
#             Do this using BioEasel's SqFile module.
#
# Args:       $exon_fafile_AR: ref to array of all exon files for this CDS
#             $acc_order_AR:   ref to array of accessions in desired order
#             $cds_fafile:     name of CDS fasta file to create
#
# Returns:    void
# Dies:       if something unexpected happens when reading the exon fasta files
#
sub combineExonsIntoCDS {
  my $sub_name = "combineExonsIntoCDS";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($exon_fafile_AR, $acc_order_AR, $cds_fafile) = @_;
  
  #printf("HEYA in combineExonsIntoCDS(), combining files to create $cds_fafile\n");
  #for(my $f = 0; $f < scalar(@{$exon_fafile_AR}); $f++) { 
  #  printf("adding $exon_fafile_AR->[$f]\n");
  #}
  
  my @sqfile_A  = ();
  my @sqname_AA = ();
  my $nfiles = scalar(@{$exon_fafile_AR});
  my %acc_map_HA   = (); # key: accession, value: array [0..$f..$nfiles-1], value: ssi index of accession in file $f
  my %acc_ct_H     = (); # key: accession, value: number of files accession exists in
  my %acc_coords_H = (); # key: accession, value: concatenated set of coordinates for this accession
  for(my $f = 0; $f < $nfiles; $f++) { 
    if(! -s $exon_fafile_AR->[$f]) { 
      die "ERROR file $exon_fafile_AR->[$f] does not exist, we need to use it to combine exons into a CDS file.";
    }
    $sqfile_A[$f] = Bio::Easel::SqFile->new({ fileLocation => $exon_fafile_AR->[$f] });
    # get names all of sequences in each file
    for(my $i = 0; $i < $sqfile_A[$f]->nseq_ssi; $i++) { 
      $sqname_AA[$f][$i] = $sqfile_A[$f]->fetch_seq_name_given_ssi_number($i);
      my ($acc, $coords) = split("/", $sqname_AA[$f][$i]);
      if(! defined $coords || $coords !~ m/\-/) { 
        die "ERROR unable to parse sequence name $sqname_AA[$f][$i] into accession and coordinates";
      }
      $acc_map_HA{$acc}[$f] = $i;
      if(! exists $acc_ct_H{$acc}) { 
        $acc_ct_H{$acc} = 1;
        $acc_coords_H{$acc} = $coords;
      }
      else { 
        $acc_ct_H{$acc}++;
        $acc_coords_H{$acc} .= "," . $coords;
      }
    }
  }
  
  # now for each accession that exists in all files, fetch all exons for that accession
  # into a new sequence
  open(OUT, ">", $cds_fafile) || die "ERROR unable to open $cds_fafile for writing";
  foreach my $acc (@{$acc_order_AR}) { 
    if(exists $acc_ct_H{$acc} && $acc_ct_H{$acc} == $nfiles) { 
      print OUT ">" . $acc . "/" . $acc_coords_H{$acc} . "\n";
      for(my $f = 0; $f < $nfiles; $f++) { 
        my $sqname = $sqname_AA[$f][($acc_map_HA{$acc}[$f])];
        my $sqonly = $sqfile_A[$f]->fetch_seq_to_fasta_string($sqname);
        $sqonly =~ s/^\>.+\n//;
        print OUT $sqonly;
      }
    }
  }
  close(OUT);
  
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

# Subroutine: matpeptParseInfile()
#
# Synopsis:   Parse the input file that defines the relationship between each CDS and the mature peptides.
#
# Args:       $infile:           file to parse
#             $cds2pmatpept_AAR: ref to array of arrays to fill here, CDS: 'primary' mat_peptide relationships
#             $cds2amatpept_AAR: ref to array of arrays to fill here, CDS: 'all'     mat_peptide relationships
#
# Returns:    void
# Dies:       if problem reading $infile (in unexpected format)
#
sub matpeptParseInfile {
  my $sub_name = "matpeptParseInfile";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($infile, $cds2pmatpept_AAR, $cds2amatpept_AAR) = @_;
  my $ncds_read_primary = 0;
  my $ncds_read_all     = 0;
  my $cds_idx2store     = 0;
  my $max_cds_idx2store = 0;

  open(IN, $infile) || die "ERROR unable to open $infile for reading in $sub_name";

  while(my $line = <IN>) { 
    if($line !~ m/^\#/) { 
      # example input file:
      ## This file explains how CDS and mat_peptide annotation for NC_001477
      ## are related.
      ##
      ############beg of file
      ## Format of lines in this file:
      ## <CDS-idx> <'primary' OR 'all'> <mat_peptide-1-idx>:<mat_peptide-2-idx>:<mat_peptide-n-idx>
      ## 'primary' lines: these define the 'primary' peptides in order, the
      ##                  CDS <CDS-idx> is comprised of the peptides listed
      ##                  in final token, which are contiguous, start of 
      ##                  first mat_peptide to stop of final mat_peptide is
      ##                  one contiguous subsequence.
      ##
      ## 'all' lines:     these define all the peptides that are ultimately
      ##                  derived from CDS <CDS-idx>. It will be a superset
      ##                  of the primary line for this index but will
      ##                  additionally include mat_peptides that are
      ##                  secondarily cleaved from the primary mat_peptides.
      ##
      #1 primary 1:3:6:7:8:9:10:11:12:13:14
      #1 all     1:2:3:4:5:6:7:8:9:10:11:12:13:14
      ## 
      ############end of file
      # NOTE: in the input file CDS and matpept indices are in coordinate space 1..N, but we store them in 0..N-1
      # 
      # we need to have one 'primary' and one 'all' line for each  
      chomp $line;
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) != 3) { 
        die "ERROR in $sub_name, unable to parse matpept input file line: $line"; 
      }
      my ($cds_idx, $primary_or_all, $matpept_str) = ($el_A[0], $el_A[1], $el_A[2]);
      $primary_or_all =~ tr/A-Z/a-z/;
      if($primary_or_all ne "all" && $primary_or_all ne "primary") { 
        die "ERROR parsing $infile, second token of each non-comment line should be 'all' or 'primary'";
      }
      my $cds_idx2store = $cds_idx - 1;
      if($cds_idx2store < 0) { 
        die "ERROR in $sub_name, read CDS idx that is 0 or less ($cds_idx) in matpept input file"; 
      }
      if($cds_idx2store > $max_cds_idx2store) { 
        $max_cds_idx2store = $cds_idx2store; 
      }
      if($primary_or_all eq "primary") { 
        if(defined $cds2pmatpept_AAR->[$cds_idx2store] || exists $cds2pmatpept_AAR->[$cds_idx2store]) {
          die "ERROR in $sub_name, two primary lines for same CDS idx ($cds_idx) in matpept input file";
        }
        my @matpept_A = split(":", $matpept_str);
        @{$cds2pmatpept_AAR->[$cds_idx2store]} = ();
        foreach my $mp (@matpept_A) { 
          push(@{$cds2pmatpept_AAR->[$cds_idx2store]}, ($mp-1));
        }
        $ncds_read_primary++;
      }
      elsif($primary_or_all eq "all") { 
        if(defined $cds2amatpept_AAR->[$cds_idx2store] || exists $cds2amatpept_AAR->[$cds_idx2store]) {
          die "ERROR in $sub_name, two all lines for same CDS idx ($cds_idx) in matpept input file";
        }
        my @matpept_A = split(":", $matpept_str);
        @{$cds2amatpept_AAR->[$cds_idx2store]} = ();
        foreach my $mp (@matpept_A) { 
          push(@{$cds2amatpept_AAR->[$cds_idx2store]}, ($mp-1));
        }
        $ncds_read_all++;
      }
    }
  }
  close(IN);

  # three sanity checks:
  # 1: we should have stored all and primary info for CDS numbers 1..$max_cds_idx2stroe
  for(my $i = 0; $i <= $max_cds_idx2store; $i++) { 
    if((! defined $cds2pmatpept_AAR->[$cds_idx2store]) && (! exists $cds2pmatpept_AAR->[$cds_idx2store])) { 
      die "ERROR in $sub_name, did not read exactly $max_cds_idx2store primary and all lines in $infile\n";
    }
    if((! defined $cds2amatpept_AAR->[$cds_idx2store]) && (! exists $cds2amatpept_AAR->[$cds_idx2store])) { 
      die "ERROR in $sub_name, did not read exactly $max_cds_idx2store primary and all lines in $infile\n";
    }
  }
  # 2: we should have at least read at least one of each primary and all
  if($ncds_read_primary == 0) { 
    die "ERROR in $sub_name, no primary CDS:mat_peptide relationships read in matpept input file $infile"; 
  }
  if($ncds_read_all == 0) { 
    die "ERROR in $sub_name, no all CDS:mat_peptide relationships read in matpept input file $infile"; 
  }
  # 3: all info should be a superset of primary info
  for(my $i = 0; $i <= $max_cds_idx2store; $i++) { 
    foreach $mp (@{$cds2amatpept_AAR->[$i]}) { 
      # make sure it exists
      my $found_it = 0;
      foreach $mp2 (@{$cds2pmatpept_AAR->[$i]}) { 
        if($mp == $mp2) { $found_it = 1; }
      }
      if(! $found_it) { 
        die sprintf("ERROR in $sub_name, all information is not a superset of primary information: %d is in primary but not all", $mp2+1);
      }
    }
  }

  return;
}

# Subroutine: matpeptValidateCdsRelationships()
#
# Synopsis:   Validate that the CDS:mat_peptide relationships read from input matpept file.
#
# Args:       $cds2pmatpept_AAR: ref to array of arrays, 1st dim value is cds_idx, value is array
#                               of mat_peptide indices that comprise that CDS, PRE-FILLED
#             $ref_cds_tbl_HAR: ref to CDS table for reference
#             $ref_mp_tbl_HAR:  ref to mat_peptide table for reference
#
# Returns:    void 
# Dies:       if relationship is not valid
#
sub matpeptValidateCdsRelationships {
  my $sub_name = "matpeptValidateCdsRelationships()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cds2pmatpept_AAR, $ref_cds_tbl_HAR, $ref_mp_tbl_HAR) = @_;

  # get CDS and mat_peptide length and coordinate strings
  my @ref_cds_len_A    = ();
  my @ref_cds_coords_A = ();
  getLengthStatsAndCoordStrings($ref_cds_tbl_HAR, \@ref_cds_len_A, \@ref_cds_coords_A);

  my @ref_mp_len_A    = ();
  my @ref_mp_coords_A = ();
  getLengthStatsAndCoordStrings($ref_mp_tbl_HAR, \@ref_mp_len_A, \@ref_mp_coords_A);
  
  # validate
  my $ref_ncds = scalar(@ref_cds_len_A);
  my $ref_nmp  = scalar(@ref_mp_len_A);
  my $ncds2check = scalar(@{$cds2pmatpept_AAR});
  my $prv_stop = undef; # previous mat_peptide's stop position, 
  for(my $cds_idx = 0; $cds_idx < scalar(@{$cds2pmatpept_AAR}); $cds_idx++) { 
    if(($cds_idx < 0) || ($cds_idx >= $ref_ncds)) { 
      die "ERROR in $sub_name, cds_idx $cds_idx is out of range"; 
    }
    my @cds_starts_A  = ();
    my @cds_stops_A   = ();
    my @cds_strands_A = ();
    my $cds_nexons    = 0;
    startStopsFromCoords($ref_cds_coords_A[$cds_idx], undef, undef, \@cds_starts_A, \@cds_stops_A, \@cds_strands_A, \$cds_nexons);
    my $cds_start = $cds_starts_A[0];
    my $cds_stop  = $cds_stops_A[$cds_nexons-1];
    if($cds_nexons != 1) { 
      if($cds_nexons != 2) { 
        die "ERROR in $sub_name, triple or more exon CDS broken up to make mat_peptides, code for this does not yet exist.";
      }
      if($cds_strands_A[0] ne $cds_strands_A[1]) { 
        die "ERROR in $sub_name, double exon CDS with each exon on different strands to make mat_peptides, code for this does not yet exist.";
      }
      # two exon CDS, if any introns exist (any nt is not included between $cds_start..$cds_stop) then we can't handle it
      # example of a multi-'exon' CDS that we CAN handle is West Nile Virus CDS #2: NC_001563.2 join(97..3540,3540..3671)
      if($cds_strands_A[0] eq "+") { 
        if(($cds_starts_A[1] - $cds_stops_A[0] - 1) > 0) { 
          die "ERROR in $sub_name, multiple exon CDS with an intron broken up to make mat_peptides, code for this does not yet exist.";
        }
      }
      else { 
        if(($cds_stops_A[0] - $cds_starts_A[1] - 1) > 0) { 
          die "ERROR in $sub_name, multiple exon CDS with an intron broken up to make mat_peptides, code for this does not yet exist.";
        }
      }
    }

    # look at all mat_peptides that are supposed to comprise this CDS
    # and make sure that they do
    my $nmp2check = scalar(@{$cds2pmatpept_AAR->[$cds_idx]});
    for(my $x = 0; $x < $nmp2check; $x++) { 
      my $mp_idx = $cds2pmatpept_AAR->[$cds_idx][$x];
      if($mp_idx < 0 || $mp_idx >= $ref_nmp) { 
        die "ERROR in $sub_name, mp_idx $mp_idx for cds_idx $cds_idx is out of range"; 
      }
      my @mp_starts_A = ();
      my @mp_stops_A  = ();
      my $mp_nexons   = 0;
      startStopsFromCoords($ref_mp_coords_A[$mp_idx], undef, undef, \@mp_starts_A, \@mp_stops_A, undef, \$mp_nexons);
      if($x == 0) { # verify start matches with CDS start
        if($mp_starts_A[0] != $cds_start) { 
          die "ERROR in $sub_name, for cds_idx $cds_idx start of first mat_peptide doesn't match CDS start ($mp_starts_A[0] != $cds_start)"; 
        }
      }
      if($x > 0) { # check that this mat_peptide is adjacent to previous one
        if($mp_starts_A[0] != ($prv_stop+1)) { 
          die sprintf("ERROR in $sub_name, for mat_peptides %d and %d are not adjacent (%d != %d+1)", $x, $x-1, $mp_starts_A[0], $prv_stop);
        }
      }
      if($x == ($nmp2check-1)) { # verify stop matches with CDS stop-3
        if(($mp_stops_A[($mp_nexons-1)]+3) != $cds_stop) { 
          die sprintf("ERROR in $sub_name, for cds_idx $cds_idx stop of final mat_peptide doesn't match CDS stop (%d != %d)", $mp_stops_A[($mp_nexons-1)], $cds_stop);
        }
      }
      $prv_stop = $mp_stops_A[($mp_nexons-1)];
      # printf("checked mp $mp_idx %d..%d\n", $mp_starts_A[0], $mp_stops_A[($mp_nexons-1)]);
    }
  }
  return;
}


# Subroutine: matpeptCheckCdsRelationships()
#
# Synopsis:   Check that the CDS:mat_peptide relationships hold for accession $seq_accn.
#
# Args:       $sqfile:          Bio::Easel::SqFile object, the sequence file to fetch start/stop codons from
#             $seq_accn:        sequence accession, 2nd dim key for %{$start_HHR} and %{$stop_HHR}.
#             $totlen:          total length of $seq_accn
#             $mdl_AR:          ref to array of model names
#             $cds2matpept_AR:  ref to array of mat_peptides this CDS is comprised of
#             $start_HHR:       ref to 2D hash of predicted start positions for mat_peptides
#             $pstop_HHR:       ref to 2D hash of predicted stop  positions for mat_peptides
#             $cstop_HHR:       ref to 2D hash of corrected stop  positions for mat_peptides, can be undef
#             $strand_HHR:      ref to 2D hash of predicted strand for mat_peptides
#             $mp2first_mdl_AR: ref to array, [0..$i..$nmp-1] = $m, $m is first model for mat_peptide $i
#             $mp2final_mdl_AR: ref to array, [0..$i..$nmp-1] = $m, $m is final model for mat_peptide $i
#
# Returns:    7 values:
#             $start:         1st position of CDS, inferred from mat_peptide predictions, 
#                             undef if no first mat_peptide prediction
#             $stop:          final position of CDS, inferred from mat_peptide predictions
#                             this is 3 nt past stop position of final mat_peptide,
#                             undef if no final mat_peptide prediction
#             $corr_stop:     corrected stop, final position of CDS, earliest (5'-most) stop of
#                             any mat_peptide that comprises the CDS for which the stop was
#                             corrected, undef if no mat_peptide stop was corrected OR if
#                             input $cstop_HHR variable is undef
#             $corr_stop_idx: model index the corrected stop exists in, undef if $corr_stop is undef
#             $len:           length of predicted CDS, inferred from mat_peptide predictions, 
#                             undef if any mat_peptide predictions are missing
#                             or any are non-adjacent
#             $start_codon:   1st codon of CDS, inferred from mat_peptide predictions, 
#                             undef if no first mat_peptide prediction
#             $stop_codon:    final codon of CDS, inferred from mat_peptide predictions, 
#                             undef if no final mat_peptide prediction, if $corr_stop is
#                             defined, this will be the early stop codon.
#             $pass_fail:     'P' if CDS passes (all mat_peptides exist and are contiguous)
#                             'F' if CDS fails (not all mat_peptides exist, and/or they are not
#                             contiguous
#
# Dies:       If not all predicted mat_peptides are on the positive strand.
#
sub matpeptCheckCdsRelationships {
  my $sub_name = "matpeptCheckCdsRelationships()";
  my $nargs_exp = 11;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seq_accn, $totlen, $mdl_AR, $cds2matpept_AR, $start_HHR, $pstop_HHR, $cstop_HHR, $strand_HHR, $mft2first_mdl_AR, $mft2final_mdl_AR) = @_;

  # for each mat_peptide that comprises the cds, make sure it's start and stops are contiguous
  # with adjacent mat_peptide. If any mat_peptides are not predicted, we fail
  my %fail_H   = (); # we will populate this with the mat_peptide index of all failures
  my %nopred_H = (); # we will populate this with the mat_peptide index of any not predicted mat_peptides

  my $start       = undef; # start position of first mat_peptide, remains undef if no prediction for first mat_peptide
  my $stop        = undef; # stop  position of final mat_peptide, remains undef if no prediction for final mat_peptide
  my $length      = undef; # length of full CDS, only defined if we have adjacent predictions for all mat_peptides
  my $cur_length  = 0;     # current length of all mat_peptides seen thus far, $length is set to this all mat_peptides are adjacent
  my $start_codon = undef; # we'll fetch this when we get to it
  my $stop_codon  = undef; # we'll fetch this when we get to it
  my $nmp2check   = scalar(@{$cds2matpept_AR});
  my $prv_stop    = undef; # stop position of previous mat_peptide
  my $corrected_stop = undef; # earliest (5'-most) corrected stop of all mat_peptides that comprise this CDS, or undef if none are corrected
  my $corrected_stop_mdl_idx = undef; # model index that corrected_stop occurs in

  for(my $x = 0; $x < $nmp2check; $x++) { 
    my $mp_idx = $cds2matpept_AR->[$x];
    my $mdl1   = $mdl_AR->[$mft2first_mdl_AR->[$mp_idx]];
    my $mdl2   = $mdl_AR->[$mft2final_mdl_AR->[$mp_idx]];

    # check to see if the stop position for this mat_peptide was corrected due to an early stop
    for(my $tmp_mdl_idx = $mft2first_mdl_AR->[$mp_idx]; $tmp_mdl_idx <= $mft2final_mdl_AR->[$mp_idx]; $tmp_mdl_idx++) { 
      my $tmp_mdl = $mdl_AR->[$tmp_mdl_idx];
      if(defined $strand_HHR->{$tmp_mdl}{$seq_accn} && $strand_HHR->{$tmp_mdl}{$seq_accn} ne "+") { 
        die "ERROR in $sub_name, strand for mat_peptide $mp_idx is not positive! (seq_accn: $seq_accn)"; 
      }
      if((defined $cstop_HHR) && 
         (defined $cstop_HHR->{$tmp_mdl}{$seq_accn}) &&
         ($cstop_HHR->{$tmp_mdl}{$seq_accn} != $pstop_HHR->{$tmp_mdl}{$seq_accn})) { 
        if((! defined $corrected_stop) || 
           ($cstop_HHR->{$tmp_mdl}{$seq_accn} < $corrected_stop)) { 
          $corrected_stop = $cstop_HHR->{$tmp_mdl}{$seq_accn};
          $corrected_stop_mdl_idx = $tmp_mdl_idx;
        }
      }
    }     

    if((! exists $start_HHR->{$mdl1}{$seq_accn}) || (! exists $start_HHR->{$mdl2}{$seq_accn})) {
      $nopred_H{$mp_idx} = 1;
    }
    else { 
      my $cur_start  = $start_HHR->{$mdl1}{$seq_accn};
      my $cur_stop   = (defined $cstop_HHR) ? $cstop_HHR->{$mdl2}{$seq_accn} : $pstop_HHR->{$mdl2}{$seq_accn};
      # special case, if this is the final peptide, always use predicted stop, this way we'll return predicted 
      # stop in $stop return variable, and possibly corrected stop in $corrected_stop return variable, which 
      # is what caller should expect.
      if($x == $nmp2check-1) { 
        $cur_stop = $pstop_HHR->{$mdl2}{$seq_accn};
      }

      my $cur_length += abs($cur_stop - $cur_start) + 1;
      if(($strand_HHR->{$mdl1}{$seq_accn} ne "+") || ($strand_HHR->{$mdl2}{$seq_accn} ne "+")) { 
        die "ERROR in $sub_name, strand for mat_peptide $mp_idx is not positive! (seq_accn: $seq_accn)"; 
      }
      my $cur_strand = $strand_HHR->{$mdl1}{$seq_accn};

      if(defined $prv_stop) { 
        if($cur_start != ($prv_stop + 1)) { 
          $fail_H{$mp_idx} = 1;
        }
      }
      if($x == 0) { # get start codon
        $start_codon = fetchStartCodon($sqfile, $seq_accn, $cur_start, $totlen, $strand_HHR->{$mdl2}{$seq_accn});
        $start = $cur_start;
      }
      if($x == $nmp2check-1) { # get stop codon
        $stop = ($strand_HHR->{$mdl2}{$seq_accn} eq "+") ? $cur_stop + 3 : $cur_stop - 3;
        # since stop codon occurs 3 nt 3' of final peptide's stop, we have to make sure it's actually in the sequence first
        # (matpept/nodup option combo is currently required (that is -nodup is req'd with -matpept), but if its ever relaxed, 
        #  we may want to rethink this, do we want to allow the stop to wrap the stop/start boundary if -nodup is not used?)
        if(($strand_HHR->{$mdl2}{$seq_accn} eq "+" && ($stop > $totlen)) || # stop is off the end of the sequence on + strand
           ($strand_HHR->{$mdl2}{$seq_accn} eq "-" && ($stop < 1))) {       # stop is off the end of the sequence on - strand
          $stop       = undef;
          $stop_codon = undef;
          $fail_H{$mp_idx} = 1;
        }
        else { 
          $stop_codon = fetchStopCodon($sqfile, $seq_accn, $stop, $totlen, $strand_HHR->{$mdl2}{$seq_accn});
        }
      }
      $prv_stop = $cur_stop;
    }
  }
  # printf("FULL $start..$stop\n");
  my $pass_fail = "F";

  if((scalar(keys %fail_H) == 0) && (scalar(keys %nopred_H) == 0)) { 
    $length = $cur_length;
    $pass_fail = "P";
  }
  
  if(defined $corrected_stop) { 
    $stop_codon = fetchStopCodon($sqfile, $seq_accn, $corrected_stop, $totlen, "+");
  }

  return($start, $stop, $corrected_stop, $corrected_stop_mdl_idx, $length, $start_codon, $stop_codon, $pass_fail);
}

# Subroutine: checkForSpanningExon()
#
# Synopsis:   Check that two exons are really not just one that spans the stop/start boundary
#             in a circular genome.
#
# Args:       $starts_AR: ref of array of start positions to potentially overwrite (if we find this is really only one exon)
# Args:       $stops_AR:  ref of array of stop positions to potentially overwrite (if we find this is really only one exon)
#             $nexons_R:  ref to scalar of number of exons to overwrite (if we find this is really only one exon)
#             $strand:    strand the exons are on
#             $totlen:    total length of the sequence
#
# Returns:    Nothing, but potentially updates @{$starts_AR}, @{$stops_AR} and $$nexons_R.
#
sub checkForSpanningExon {
  my $sub_name = "checkForSpanningExon()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($starts_AR, $stops_AR, $nexons_R, $strand, $totlen) = @_;

  if($$nexons_R == 2) { 
    # if we're in a circular genome, we need to check for a special case, where 
    # what looks like a 2-exon CDS is really a single exon that spans the stop..start boundary.
    # [Note that if the stop..start boundary is spanned by an intron (i.e. exon i is before stop,
    # and i+1 is after start) then we don't need to modify anything, we'll still fetch the proper
    # sequence even in a duplicated genome].
    #
    # Example 1: single exon that spans stop..start boundary on positive strand
    # join(2309..3182,1..1625) in a seq of length 3182, this should really be a single exon
    # $starts_A[0] = 2309
    # $stops_A[0]  = 3182
    # $starts_A[1] = 1
    # $stops_A[1]  = 1625
    # $nexons = 2;
    # 
    # should become:
    # $starts_A[0] = 2309
    # $stops_A[0]  = 4807
    # $nexons = 1;
    # 
    # Example 2: single exon that spans stop..start boundary on negative strand
    # complement(join(2309..3182,1..1625))   in a seq of length 3182, this should really be a single exon
    # $starts_A[0] = 3182
    # $stops_A[0]  = 2309
    # $starts_A[1] = 1625
    # $stops_A[1]  = 1
    # $nexons = 2;
    # 
    # should become:
    # $starts_A[0] = 4807
    # $stops_A[0]  = 2309
    # $nexons = 1;
    #
    # we can easily check and fix these cases:
    my $tmp_start = undef;
    my $tmp_stop  = undef;
    if($strand eq "+" && $stops_AR->[0] == $totlen && $starts_AR->[1] == 1) { 
      $tmp_start = $starts_AR->[0];
      $tmp_stop  = $stops_AR->[1] + $ref_totlen;
    }
    elsif($strand eq "-" && $starts_AR->[0] == $totlen && $stops_AR->[1] == 1) { 
      my $tmp_start = $starts_AR->[1] + $totlen;
      my $tmp_stop  = $stops_AR->[0];
    }
    if(defined $tmp_start && defined $tmp_stop) { 
      @{$starts_AR} = ();
      @{$stops_AR} = ();
      $starts_AR->[0] = $tmp_start;
      $stops_AR->[0]  = $tmp_stop;
      $$nexons_R = 1;
    }    
  }
  return;
}

# Subroutine: setErrorCode()
#
# Synopsis:   Set an error code value to either '1' (has the error) or '0' (does not have the error)
#             given references to the relevant data structures.
#
# Args:       $err_AR:       ref to array to update one element of
#             $err_extra_AR: ref to array to add 'extra' information to
#             $errcode:      error code of string we want to update error info for
#             $errextra:     if defined, extra information to add to $err_extra_AR
#             $code2idx_HR:  ref to hash that maps $errcode to array index in @{$err_AR}
#             $errctr_R:     if defined, reference to error counter, to increment by 1
#
# Returns:    void
# Dies:       if $code2idx_HR->{$errcode} is not defined

sub setErrorCode {
  my $sub_name = "setErrorCode()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_AR, $err_extra_AR, $errcode, $errextra, $code2idx_HR, $errctr_R) = @_;

  #printf("HEYA in $sub_name\n");
  #printf("HEYA \terrcode:  $errcode\n");
  #printf("HEYA \terrctr:   $$errctr_R\n");
  #printf("HEYA \terrextra: $errextra\n");

  if(! defined $code2idx_HR->{$errcode}) { die "ERROR in $sub_name, unrecognized error code $errcode"; }
  $err_AR->[$code2idx_HR->{$errcode}] = 1;
  if(defined $err_extra_AR && defined $errextra) { 
    $err_extra_AR->[$code2idx_HR->{$errcode}] = $errextra;
  }
  if(defined $errctr_R) { 
    $$errctr_R++; 
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
#             $cds2pmatpept_AAR:    ref to @cds2pmatpept_AA, already filled
#             $out_col_header_AAR: ref to 2D array of column headers, filled here
#                                  undef unless $do_seqrow is '1'
#             $out_row_header_AR:  ref to 1D array of row headers, filled here
#                                  undef unless $do_seqcol is '1'
#             $out_header_exp_AR:  ref to 1D array of header explanations, each
#                                  element is a line to be printed in explanatory
#                                  section of the output; filled here
sub getHeadings {
  my $sub_name = "getHeadings";
  my $nargs_exp = 22;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($do_seqrow, $do_seqcol, $do_matpept, $do_nofid, $do_mdlb, $do_noss3, $do_nostop, $do_fullolap, $do_fulladj, $origin_seq, $ref_tot_nexons, $nmdl, $mdl2mft_map_AR, $mft2exon_map_AR, $mdl_is_first_AR, $mdl_is_final_AR, $mft_out_short_AR, $mft_out_product_AR, $cds2pmatpept_AAR, $out_col_header_AAR, $out_row_header_AR, $out_header_exp_AR) = @_;

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

  # determine how many features we have, crudely from the $mdl2mft_map_AR
  my $nmft = $mdl2mft_map_AR->[$nmdl-1] + 1; # the + 1 is to correct for the fact that mdl2mft_map_AR is in 0..$nmft-1 coordinate space not 1..$nmft

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

  my @pf_text_A = (); # array of lines to print about pass/fail strings to explanation at the end
  my $pf_idx = 1;

  my %need_to_define_H = (); # hash of definitions

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
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "number of occurrences of origin sequence (input with -oseq) in genome");

    # column/row #5: 'origin sequence:start'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %5s", "start");
    $tok5 = sprintf(" %5s", "-----");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "start position of lone occurrence of origin sequence (if only 1 exists)");

    # column/row #6: 'origin sequence:stop'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %5s", "stop");
    $tok5 = sprintf(" %5s", "-----");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "stop  position of lone occurrence of origin sequence (if only 1 exists)");

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
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "'P' (for PASS) if there is exactly 1 occurrence of the offset, else 'F' for FAIL");
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef); # adds a blank line
    push(@pf_text_A, "P/F character $pf_idx pertains to the origin sequence test");
    $pf_idx++;
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
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "start position of 5' UTR (inferred from other predictions, \"?\" if first mat_peptide is not predicted)");

    $tok4 = sprintf(" %6s", "stop");
    $tok5 = sprintf(" ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "stop  position of 5' UTR (inferred from other predictions, \"?\" if first mat_peptide is not predicted)");

    $tok4 = sprintf(" %6s", "length");
    $tok5 = sprintf(" ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "length of 5' UTR (inferred from other predictions, \"?\" if first mat_peptide is not predicted)");
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
      getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "start position of $exp_substr (\"NP\" if no prediction)");
      getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "enclosed in brackets \"\[e\]\" if start/stop different from all exon start/stops in existing GenBank annotation");
    }
    
    # stop, fid, and md rows take place for all exons
    # only token 4 changes
    $tok4 = sprintf(" %8s", sprintf("%s%s", "stop", $mft2exon_map_AR->[$h]+1));
    $exp_tok4 = "stop<j>";
    $tok5 = sprintf(" %8s", "--------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
    if($do_explanation) { 
      getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "stop  position of $exp_substr (\"NP\" if no prediction)");
      getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "enclosed in brackets \"\[e\]\" if start/stop different from all exon start/stops in existing GenBank annotation");
    }
    
    if(! $do_nofid) { 
      $tok4 = sprintf(" %5s", sprintf("%s%s", "fid", $mft2exon_map_AR->[$h]+1));
      $exp_tok4 = "fid<j>";
      $tok5 = sprintf(" %5s", "-----");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "fractional identity between $exp_substr and reference genome (\"NP\" if no prediction)");
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
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "annotation indicating if alignment to reference extends to 5' and 3' end of reference $exp_substr.");
        getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "first character pertains to 5' end and second character pertains to 3' end.");
        getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "possible values for each of the two characters:");
        getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "  \".\":   alignment extends to boundary of reference");
        getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "  \"<d>\": alignment truncates <d> nucleotides short of boundary of reference (1 <= <d> <= 9)");
        getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "  \"<d>\": alignment truncates <d> nucleotides short of boundary of reference (1 <= <d> <= 9)");
        getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "  \"+\":   alignment truncates >= 10 nucleotides short of boundary of reference");
        if(! $do_matpept) { 
          getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "  \"c\":   position has been corrected based on predicted protein sequence");
          getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "           3' position: stop coordinate has been adjusted to first in-frame stop");
          getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "  \"-\":   exon is not predicted due to stop codon in earlier exon");
        }
        getHeadingsExplanationHelper($out_header_exp_AR, undef,     undef,     undef, "  \"NP\":  (spanning both characters) no prediction");
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
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "\"NP\" if no prediction");
        $need_to_define_H{"overlap"} = 1;
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
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "\"NP\" if no prediction");      
      }
      $need_to_define_H{"adjacent"} = 1;
    }

    $exp_substr = $do_matpept ? "mat_peptide coding sequence" : "CDS";
    if($mdl_is_final_AR->[$h]) { 
      $tok4 = sprintf(" %6s", "length");
      $tok5 = sprintf(" %6s", "------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $tok4, undef, sprintf("length of $exp_substr #<i> (all %s summed)", $do_matpept ? "segments" : "exons"));
      }      

      if((! $do_matpept) && (! $do_noss3)) { # skip this in matpept mode, we don't check start/stop of mat_peptides, only CDS, later
        $tok4 = sprintf(" %3s", "ss3");
        $tok5 = sprintf(" %3s", "---");
        if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
        elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
        if($do_explanation) { 
          getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if predicted CDS has a valid start codon, stop codon and is a multiple of 3");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "first  character: '.' if predicted CDS has a valid start codon, else '!'");          
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "second character: '.' if predicted CDS has a valid stop  codon, else '!'");      
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "third  character: '.' if predicted CDS has a length which is a multiple of three, else '!'");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "\"NP\" if no prediction");
        }
      }
      if((! $do_matpept) && (! $do_nostop)) { # skip this in matpept mode, we only check stop of final mat_peptide, later
        $tok4 = sprintf(" %3s", "stp");
        $tok5 = sprintf(" %3s", "---");
        if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
        elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); }
        if($do_explanation) { 
          getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "the predicted stop codon for this CDS (\"NP\" if no prediction)");
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
          push(@pf_text_A, sprintf("P/F characters %d to %d pertain to each of the %d mature peptides, in order.", $pf_idx, $pf_idx + $nmft-1, $nmft));
          $pf_idx += $nmft;
          $need_to_define_H{"adjacent"} = 1;
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
          push(@pf_text_A, sprintf("P/F characters %d to %d pertain to each of the %d CDS, in order.", $pf_idx, $pf_idx + $nmft-1, $nmft));
          $pf_idx += $nmft;
          $need_to_define_H{"overlap"} = 1;
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
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "start position of 3' UTR (inferred from other predictions, \"?\" if final mat_peptide is not predicted)");

    $tok4 = sprintf(" %6s", "stop");
    $tok5 = sprintf(" ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "stop  position of 3' UTR (inferred from other predictions, \"?\" if final mat_peptide is not predicted)");

    $tok4 = sprintf(" %6s", "length");
    $tok5 = sprintf(" ------");
    if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
    elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
    getHeadingsExplanationHelper($out_header_exp_AR, $tok2, $tok4, undef, "length of 3' UTR (inferred from other predictions, \"?\" if final mat_peptide is not predicted)");

    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
  }

  # create columns for CDS annotation in mat_peptide mode, if nec
  if($do_matpept) { 
    $do_explanation = 1;
    my $ncds_for_matpept = scalar(@{$cds2pmatpept_AAR});
    foreach (my $cds_idx = 0; $cds_idx < scalar(@{$cds2pmatpept_AAR}); $cds_idx++) { 
      my $cds_idx2print = $cds_idx+1;
      $width = 6 + 1 + 6 + 1 + 6; #20
      $tok1     = sprintf("  %*s", $width, "");
      $tok2     = sprintf("  %*s", $width, "CDS #$cds_idx2print" . monocharacterString(($width-length("CDS #$cds_idx2print")), " "));
      my $exp_tok2 = "CDS #<i>";
      $tok3 = sprintf("  %s", monocharacterString($width, "-"));
      $tok4 = sprintf("  %6s", "start");
      $tok5 = sprintf("  ------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                $tok1, $tok2, $tok3, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok2, $tok4, undef, "start position of CDS #<i> (inferred from mat_peptides that comprise it,");
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok2, $tok4, undef, "\"?\" if first mat_peptide for this CDS is not predicted)");
      }

      $tok4 = sprintf(" %6s", "stop");
      $tok5 = sprintf(" %6s", "------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok2, $tok4, undef, "stop  position of CDS #<i> (inferred from mat_peptides that comprise it,");
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok2, $tok4, undef, "\"?\" if final mat_peptide for this CDS is not predicted)");
      }

      $tok4 = sprintf(" %6s", "length");
      $tok5 = sprintf(" %6s", "------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok2, $tok4, undef, "length of CDS #<i> (\"?\" if any of the mat_peptides that comprise this CDS are not predicted)");
      }

      $tok4 = sprintf(" %6s", "startc");
      $tok5 = sprintf(" %6s", "------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok2, $tok4, undef, "start codon of CDS #<i> (\"?\" if first mat_peptide for this CDS is not predicted)");
      }

      $tok4 = sprintf(" %6s", "stopc");
      $tok5 = sprintf(" %6s", "------");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok2, $tok4, undef, "stop codon of CDS #<i> (\"?\" if final mat_peptide for this CDS is not predicted)");
      }

      if(! $do_noss3) { 
        $tok4 = sprintf(" %3s", "ss3");
        $tok5 = sprintf(" %3s", "---");
        if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
        elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
        if($do_explanation) { 
          getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if predicted CDS has a valid start codon, stop codon and is a multiple of 3");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "first  character: '.' if predicted CDS has a valid start codon, '!' if not,");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "                  and '?' if first mat_peptide for this CDS is not predicted");          
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "second character: '.' if predicted CDS has a valid stop  codon, '!' if not,");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "                      and '?' if final mat_peptide for this CDS is not predicted");      
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "third  character: '.' if predicted CDS has a length which is a multiple of three, '!' if it is not a");
          getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "                  multiple of 3, and '?' if any of the mat_peptides that comprise it are not predicted.");
        }
      }

      $tok4 = sprintf(" %6s", "PF");
      $tok5 = sprintf(" %6s", "---");
      if   ($do_seqrow) { getHeadingsSeqRowHelper($out_col_header_AAR,                undef, undef, undef, $tok4, $tok5); }
      elsif($do_seqcol) { getHeadingsSeqColHelper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); }
      if($do_explanation) { 
        getHeadingsExplanationHelper($out_header_exp_AR, $exp_tok2, $tok4, undef, "annotation indicating if this CDS PASSED ('P') or FAILED ('F')");
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  a CDS sequence PASSES ('P') if and only if all of the following");
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  conditions are met (else it FAILS):");
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (1) it has a valid start codon at beginning of its first mat_peptide");
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (2) it has a valid stop  codon immediately after the end of its final mat_peptide");
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (3) its length is a multiple of 3");
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "  (4) all of the mat_peptides that comprise it are adjacent");
        getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
        if($ncds_for_matpept == 1) { 
          push(@pf_text_A, sprintf("P/F character %d pertains to the lone CDS.", $pf_idx));
        }
        else {
          push(@pf_text_A, sprintf("P/F characters %d to %d pertain to each of the %d CDS, in order.", $pf_idx, $pf_idx + $ncds_for_matpept-1, $ncds_for_matpept));
        }
        $pf_idx += $ncds_for_matpept;
      }
      $do_explanation = 0;
    }
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
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "first character:   'P' for PASS if predicted annotation for this accession has same overlaps as the reference");
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "                   'F' for FAIL if it does not");
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("second character:  total number of overlaps between any two %s", $do_matpept ? "mat_peptides" : "exons"));
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("remainder of line: text explaining which exons overlap", $do_matpept ? "mat_peptides" : "exons"));
    if($do_matpept) { 
      getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("  e.g.: \"3.1|4.1\" indicates segment #1 of mat_peptide #3 overlaps with exon #1 of mat_peptide #4 on either strand"));
    }
    else { 
      getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("  e.g.: \"3.2/4.1\" indicates exon #2 of CDS #3 overlaps with exon #1 of CDS #4 on either strand"));
    }
    push(@pf_text_A, sprintf("P/F character $pf_idx pertains to the overlap test of all %s",  $do_matpept ? "mat_peptides" : "exons"));
    $pf_idx++;
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
    $need_to_define_H{"overlap"} = 1;
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
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("two %s i and j are adjacent if i < j and final nt of i is 1 less than first nt of j", $do_matpept ? "mat_peptides" : "exons"));
      getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "first character:   'P' for PASS if predicted annotation for this accession has same adjacencies as the reference");
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "                   'F' for FAIL if it does not");
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("second character:  total number of adjacencies between any two %s", $do_matpept ? "mat_peptides" : "exons"));
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("remainder of line: text explaining which %s are adjacent", $do_matpept ? "mat_peptides" : "exons"));
    if($do_matpept) { 
      getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("  e.g.: \"3.1|4.1\" indicates segment #1 of mat_peptide #3 is adjacent to segment #1 of mat_peptide #4 on either strand"));
    }
    else { 
      getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("  e.g.: \"3.2/4.1\" indicates exon #2 of CDS #3 is adjacent with exon #1 of CDS #4 on either strand"));
    }
    push(@pf_text_A, sprintf("P/F character $pf_idx pertains to the full adjacency test of all %s", $do_matpept ? "mat_peptides" : "exons"));
    $pf_idx++;
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
    $need_to_define_H{"adjacent"} = 1;
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
  getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("as indicated in the \"PF\" %s.", ($do_seqcol) ? "rows" : "columns"));
  getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, sprintf("After that is a string of %d characters, these are the individual P/F results in order.", $pf_idx-1));
  foreach my $pf_text_str (@pf_text_A) { 
    getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, $pf_text_str);
  }
  getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, "See explanations of the individual P/F values above.");

  if((defined $need_to_define_H{"overlap"}) || (defined $need_to_define_H{"adjacent"})) {
     getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
     push(@{$out_header_exp_AR}, "# Definitions of non-obvious terms above:\n");
     getHeadingsExplanationHelper($out_header_exp_AR, undef, undef, undef, undef);
     if(defined $need_to_define_H{"overlap"}) { 
       push(@{$out_header_exp_AR}, "# overlap:  two features i and j overlap if they are on both on the same strand and overlap by >= 1 nt.\n");
     }
     if(defined $need_to_define_H{"adjacent"}) { 
       push(@{$out_header_exp_AR}, "# adjacent: two features i and j are adjacent if they are on the same strand and\n");
       push(@{$out_header_exp_AR}, "#           start_i < start_j and stop_i+1 == start_j.\n");
       push(@{$out_header_exp_AR}, "#           (Note that on the positive strand start_i <= stop_i for all i,\n");
       push(@{$out_header_exp_AR}, "#           (and that  on the negative strand start_i >= stop_i for all i)\n");
     }
  }
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
#             $totlen:              total length of accession
#             $do_nodup:            '1' if we're not duplicating this accession
#             $mft_tbl_HHAR:        ref to mft_tbl_HHA
#             $act_exon_starts_AAR: FILLED HERE [0..$nmft-1][0..$nexons-1] start positions of actual annotations of exons for this accn, $nexons is CDS specific
#             $act_exon_stops_AAR:  FILLED HERE [0..$nmft-1][0..$nexons-1] stop  positions of actual annotations of exons for this accn, $nexons is CDS specific
#             $tot_nexons_R:        FILLED HERE, total number of annotated exons for this accession
#             $nmft_R:              FILLED HERE, number of CDS for this accession
#
sub getActualAnnotations { 
  my $sub_name = "getActualAnnotations";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($accn, $totlen, $do_nodup, $mft_tbl_HHAR, $act_exon_starts_AAR, $act_exon_stops_AAR, $tot_nexons_R, $nmft_R) = @_;

  my ($nmft, $npos, $nneg, $nunc, $nbth, $strand_str, $tot_nexons); 

  if(exists ($mft_tbl_HHAR->{$accn})) { 
    ($nmft, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats($mft_tbl_HHAR, $accn);
    my @mft_len_A = ();
    my @mft_coords_A = ();
    my @mft_product_A = ();
    getLengthStatsAndCoordStrings($mft_tbl_HHAR->{$accn}, \@mft_len_A, \@mft_coords_A);
    getQualifierValues($mft_tbl_HHAR, $accn, "product", \@mft_product_A);
    for(my $i = 0; $i < $nmft; $i++) { 
      # determine start and stop positions of all exons
      my @starts_A = ();
      my @stops_A  = ();
      my $nexons   = 0;
      @{$act_exon_starts_AAR->[$i]} = ();
      @{$act_exon_stops_AAR->[$i]}  = ();
      startStopsFromCoords($mft_coords_A[$i], $totlen, $do_nodup, \@starts_A, \@stops_A, undef, \$nexons);
      
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
#          $oseq_ct:       number of occurrences of origin sequence found
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
# Args:       $FH:                file handle to print to
#             $out_header_exp_AR: ref to array of output explanation lines
#
# Synopsis:   Prints out output lines in @{$out_header_exp_AR} and exits.
#
# Returns:    void

sub outputColumnHeaderExplanations {
  my $sub_name = "outputColumnHeaderExplanations";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($FH, $out_header_exp_AR) = @_;

  foreach my $line (@{$out_header_exp_AR}) { 
    print $FH $line;
  }

  return;
}

# Subroutine: outputCDSMaturePeptideRelationships()
# Args:       $FH:                file handle to print to
#             $cds2pmatpept_AAR:   ref to array of arrays describing relatinoships
#             $div_line:          divider line
# Synopsis:   Prints out explanation of CDS:mature peptides relationships.
#
# Returns:    void

sub outputCDSMaturePeptideRelationships { 
  my $sub_name = "outputCDSMaturePeptideRelationships";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($FH, $cds2pmatpept_AAR, $cds2amatpept_AAR, $div_line) = @_;

  printf $FH ("#\n");
  printf $FH ("# CDS:MAT_PEPTIDE relationships:\n");
  printf $FH ("#\n");
  for(my $cds_idx = 0; $cds_idx < scalar(@{$cds2pmatpept_AAR}); $cds_idx++) { 
    printf $FH ("# CDS #%d is comprised of the following primary mat_peptides in order: ", ($cds_idx+1));
    for(my $i = 0; $i < scalar(@{$cds2pmatpept_AAR->[$cds_idx]}); $i++) { 
      if($i > 0) { print $FH (", "); }
      print $FH ($cds2pmatpept_AAR->[$cds_idx][$i] + 1);
    }
    print $FH "\n";
  }
  for(my $cds_idx = 0; $cds_idx < scalar(@{$cds2pmatpept_AAR}); $cds_idx++) { 
    printf $FH ("# CDS #%d encodes all of the following mat_peptides in order: ", ($cds_idx+1));
    for(my $i = 0; $i < scalar(@{$cds2amatpept_AAR->[$cds_idx]}); $i++) { 
      if($i > 0) { print $FH (", "); }
      print $FH ($cds2amatpept_AAR->[$cds_idx][$i] + 1);
    }
    print $FH "\n";
  }
  print $FH "#\n";
  print $FH $div_line;

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
        
        if(defined $next_ins_str && defined $next_del_str) { 
          if($next_del_rfpos <= $next_ins_rfpos) { # delete comes first, print it
            $substr = "D" . $next_del_rfpos . ":" . $next_del_count;
            if($do_gap_all || (($next_del_count % 3) != 0)) { 
              if($gapstr ne "") { $gapstr .= ","; }
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
              if($gapstr ne "") { $gapstr .= ","; }
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
            if($gapstr ne "") { $gapstr .= ","; }
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
            if($gapstr ne "") { $gapstr .= ","; }
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
# Args:       $FH:            file handle to print to
#             $header_AR:     reference to array of row headers
#             $out_AAR:       reference to the 2D array of output tokens for
#                             current sequences for the page
#             $ref_out_AR:    reference to array of output tokens for reference
#             $page_idx:      page number
#
# Returns:    void
#
sub outputSeqAsColumnsPage {
  my $sub_name = "outputSeqAsColumnsPage";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($FH, $header_AR, $out_AAR, $ref_out_AR, $page_idx) = @_;

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
    if($ntok != $nrow) { 
      die sprintf("ERROR in $sub_name, we have $nrow headers, but sequence %s has $ntok tokens", $i+1, $ntok); 
    }
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
    printf $FH ("%-*s", $hwidth, $header_AR->[$r]);
    for(my $i = 0; $i <= $nseq; $i++) { 
      # first column is reference, then come the other seqs
      $AR = ($i == 0) ? $ref_out_AR : \@{$out_AAR->[$i-1]}; 
      $el = $AR->[$r];
      $el =~ s/\s+//g;
      printf $FH ("%*s", $cwidth_A[$i], $AR->[$r]);
    }
    printf $FH ("\n");
  }
  print $FH "#\n";
  printf $FH ("# end of page %d\n", $page_idx);
  print $FH "#\n";
  return;
}

###########################################################

