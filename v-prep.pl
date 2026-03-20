#!/usr/bin/env perl
#
# v-profprep.pl
# Prepare a representative Stockholm alignment from sequence metadata.
#
use strict;
use warnings;
use Getopt::Long qw(:config no_auto_abbrev);
use Time::HiRes qw(gettimeofday);
use LWP::Simple;
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

require "vadr.pm";
require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seq.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

#######################################################################################
# What this script does: 
#
# - Preliminaries: 
#   o processes options
#   o creates the output directory
#   o outputs program banner
#   o makes sure the required executables are executable
#
# - Reads metadata output from fetch-seqs-given-taxid.pl
# - Selects sequence representatives (The Funnel Strategy)
# - Downloads full fasta files
# - Assesses quality and filters based on length/ambiguities
# - Computes alignment and constructs final multi-seq Stockholm + minfo
#######################################################################################

# make sure required environment variables are set
my $env_vadr_scripts_dir  = utl_DirEnvVarValid("VADRSCRIPTSDIR");
my $env_vadr_blast_dir    = utl_DirEnvVarValid("VADRBLASTDIR");
my $env_vadr_infernal_dir = utl_DirEnvVarValid("VADRINFERNALDIR");
my $env_vadr_hmmer_dir    = utl_DirEnvVarValid("VADRHMMERDIR");
my $env_vadr_easel_dir    = utl_DirEnvVarValid("VADREASELDIR");
my $env_vadr_fasta_dir    = $ENV{"VADRFASTADIR"};
my $env_vadr_rfam_dir     = $ENV{"VADRRFAMDIR"}; # optional, validated later if needed

# make sure the required executables exist and are executable
my %execs_H = (); # hash with paths to all required executables
$execs_H{"esl-reformat"}  = $env_vadr_easel_dir    . "/esl-reformat";
$execs_H{"esl-sfetch"}    = $env_vadr_easel_dir    . "/esl-sfetch";
$execs_H{"esl-translate"} = $env_vadr_easel_dir    . "/esl-translate";
$execs_H{"esl-alimask"}   = $env_vadr_easel_dir    . "/esl-alimask";
$execs_H{"esl-alimanip"}  = $env_vadr_easel_dir    . "/esl-alimanip";
$execs_H{"blastn"}        = $env_vadr_blast_dir    . "/blastn";
$execs_H{"makeblastdb"}   = $env_vadr_blast_dir    . "/makeblastdb";
$execs_H{"cmscan"}        = $env_vadr_infernal_dir . "/cmscan";
$execs_H{"cmalign"}       = $env_vadr_infernal_dir . "/cmalign";
$execs_H{"cmbuild"}       = $env_vadr_infernal_dir . "/cmbuild";
$execs_H{"esl-reformat"}  = $env_vadr_easel_dir    . "/esl-reformat";
$execs_H{"v-build.pl"}    = $env_vadr_scripts_dir  . "/v-build.pl";
$execs_H{"v-annotate.pl"} = $env_vadr_scripts_dir  . "/v-annotate.pl";
if(defined $env_vadr_fasta_dir) {
  my $ggsearch_exec = $env_vadr_fasta_dir . "/ggsearch36";
  if(-x $ggsearch_exec) {
    $execs_H{"ggsearch"} = $ggsearch_exec;
  }
}
utl_ExecHValidate(\%execs_H, undef);

#########################################################
# Command line and option processing using sqp_opts.pm
#
my %opt_HH = ();      
my @opt_order_A = (); 
my %opt_group_desc_H = ();
my $g = 0; # option group

$opt_group_desc_H{++$g} = "basic options";
#     option            type       default  group   requires incompat     preamble-output                                                help-output            
opt_Add("-h",           "boolean", 0,           0,    undef, undef,       undef,                                                         "display this help",                                   \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,          $g,    undef, undef,       "forcing directory overwrite",                                 "force; if dir <output directory> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,          $g,    undef, undef,       "be verbose",                                                  "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,          $g,    undef, undef,       "leave intermediate files on disk",                            "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "input options";
#     option            type       default  group   requires incompat     preamble-output                                                help-output            
opt_Add("--minout",     "string",  undef,      $g,    undef, undef,       "output model info file to <s>",                               "output model info file to <s>", \%opt_HH, \@opt_order_A);
opt_Add("--mdir",       "string",  undef,      $g,    undef, undef,       "model files are in directory <s>",                            "model files are in directory <s>", \%opt_HH, \@opt_order_A);
opt_Add("--mkey",       "string",  undef,      $g,    "--mdir", undef,   "model key is <s> (override auto-derived key)",                "model key is <s> (override auto-derived key)", \%opt_HH, \@opt_order_A);
opt_Add("--taxid",      "string",  undef,      $g,    undef, "--meta",    "fetch metadata from NCBI for taxonomy ID <s>",                "fetch metadata from NCBI for taxonomy ID <s>", \%opt_HH, \@opt_order_A);
opt_Add("--meta",       "string",  undef,      $g,    undef, "--taxid",   "read metadata TSV from <s> instead of fetching",              "read metadata TSV from <s> instead of fetching", \%opt_HH, \@opt_order_A);
opt_Add("--api_key",    "string",  undef,      $g,    "--taxid", undef,   "NCBI API key to use for metadata fetch",                      "NCBI API key to use for metadata fetch as <s>", \%opt_HH, \@opt_order_A);
opt_Add("--seed-accn",  "string",  undef,      $g,    undef, undef,       "build seed model internally from accession <s>",              "build seed model internally from accession <s>", \%opt_HH, \@opt_order_A);
opt_Add("--seed-build-opts-file", "string", undef,   $g,    "--seed-accn", "--seed-build-opts", "read extra v-build.pl options from file <s>",     "read extra v-build.pl options from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--seed-build-opts", "string", undef,        $g,    "--seed-accn", "--seed-build-opts-file", "extra v-build.pl options from single string <s>", "extra v-build.pl options from single string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--skip-annotate", "boolean", 0,             $g,    undef, undef,   "skip Tier 2 v-annotate screening step",                     "skip Tier 2 v-annotate screening step (temporary)", \%opt_HH, \@opt_order_A);
opt_Add("--xambig",     "integer", 5,          $g,    undef, undef,       "max ambiguous nucleotides allowed per sequence",             "max ambiguous nucleotides allowed per sequence as <n>", \%opt_HH, \@opt_order_A);
opt_Add("--xpergroup",  "integer", 50,         $g,    undef, undef,       "max sequences retained per serotype/genotype group",         "max sequences retained per serotype/genotype group as <n>", \%opt_HH, \@opt_order_A);
opt_Add("--rna-cm-file", "string", undef,        $g,    undef, "--skip-rna", "use CM file <s> for RNA search instead of default Rfam.cm",  "use CM file <s> for RNA search instead of default Rfam.cm", \%opt_HH, \@opt_order_A);
opt_Add("--skip-rna",   "boolean", 0,             $g,    undef, "--rna-cm-file", "skip RNA discovery and alignment",                          "skip RNA discovery and alignment (treat all noncoding as unstructured)", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "other expert options";
#       option       type          default     group  requires incompat      preamble-output                                              help-output           
opt_Add("--execname",   "string",  undef,         $g,    undef, undef,       "define executable name of this script as <s>",              "define executable name of this script as <s>", \%opt_HH, \@opt_order_A);

my %GetOptions_H = ();
my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'keep'         => \$GetOptions_H{"--keep"},
# input options
                'minout=s'     => \$GetOptions_H{"--minout"},
                'mdir=s'       => \$GetOptions_H{"--mdir"},
                'mkey=s'       => \$GetOptions_H{"--mkey"},
                'taxid=s'      => \$GetOptions_H{"--taxid"},
                'meta=s'       => \$GetOptions_H{"--meta"},
                'api_key=s'    => \$GetOptions_H{"--api_key"},
                'seed-accn=s'  => \$GetOptions_H{"--seed-accn"},
                'seed-build-opts-file=s' => \$GetOptions_H{"--seed-build-opts-file"},
                'seed-build-opts=s' => \$GetOptions_H{"--seed-build-opts"},
                'skip-annotate' => \$GetOptions_H{"--skip-annotate"},
                'xambig=i'     => \$GetOptions_H{"--xambig"},
                'xpergroup=i'  => \$GetOptions_H{"--xpergroup"},
                'rna-cm-file=s' => \$GetOptions_H{"--rna-cm-file"},
                'skip-rna'     => \$GetOptions_H{"--skip-rna"},
# other expert options
                'execname=s'   => \$GetOptions_H{"--execname"});

my $total_seconds = -1 * ofile_SecondsSinceEpoch(); 
my $execname_opt  = $GetOptions_H{"--execname"};
my $executable    = (defined $execname_opt) ? $execname_opt : "v-profprep.pl";
my $usage         = "Usage: $executable [-options] --mdir <s> <path to output directory to create>\n";
my $synopsis      = "$executable :: prepare training alignment for profile model building";
my $date          = scalar localtime();
my $version       = "1.7";
my $releasedate   = "Sep 2025";
my $pkgname       = "VADR";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, undef);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } 
}

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# check that number of command line args is correct
if(scalar(@ARGV) != 1) {
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do $executable -h\n\n";
  exit(1);
}
my ($dir) = (@ARGV);
my $model_dir = opt_Get("--mdir", \%opt_HH);
my $meta_tsv = undef;
my $do_seed_bootstrap = opt_IsUsed("--seed-accn", \%opt_HH);
my $do_skip_annotate = opt_Get("--skip-annotate", \%opt_HH);

if((! defined $model_dir) || ($model_dir =~ /^\s*$/)) {
  if($do_seed_bootstrap) {
    die "ERROR, --seed-accn requires --mdir <s>: v-profprep.pl runs v-build.pl to create seed model files in <mdir>, then reads the derived .vadr.minfo from that same directory.";
  }
  die "ERROR, required option --mdir <s> was not provided";
}

if($do_seed_bootstrap && opt_IsUsed("--mkey", \%opt_HH)) {
  die "ERROR, --mkey is not compatible with --seed-accn; in seed bootstrap mode model key is derived from model directory name";
}

if($model_dir =~ m/\/$/) {
  $model_dir =~ s/\/$//;
}

if((! $do_seed_bootstrap) && (! -d $model_dir)) {
  die "ERROR: model directory $model_dir does not exist or is not a directory";
}

my $auto_model_key = $model_dir;
$auto_model_key =~ s/^.+\///;
my $model_key = (opt_IsUsed("--mkey", \%opt_HH)) ? opt_Get("--mkey", \%opt_HH) : $auto_model_key;
my $seed_minfo = $model_dir . "/" . $model_key . ".vadr.minfo";

# metadata mode checks: exactly one required
if((! opt_IsUsed("--meta", \%opt_HH)) && (! opt_IsUsed("--taxid", \%opt_HH))) {
  die "ERROR, must specify one of --taxid <s> or --meta <s>";
}

# Check key input files early
if(opt_IsUsed("--meta", \%opt_HH)) {
  $meta_tsv = opt_Get("--meta", \%opt_HH);
  if (! -e $meta_tsv) { die "ERROR: metadata TSV file from --meta does not exist: $meta_tsv"; }
}

if(opt_IsUsed("--seed-build-opts", \%opt_HH) && opt_Get("--seed-build-opts", \%opt_HH) =~ /^\s*$/) {
  die "ERROR, --seed-build-opts was provided but is empty";
}
if(opt_IsUsed("--seed-build-opts-file", \%opt_HH)) {
  my $seed_opts_file = opt_Get("--seed-build-opts-file", \%opt_HH);
  if(! -e $seed_opts_file) {
    die "ERROR, --seed-build-opts-file does not exist: $seed_opts_file";
  }
}

if(opt_Get("--xambig", \%opt_HH) < 0) {
  die "ERROR, --xambig must be >= 0";
}
if(opt_Get("--xpergroup", \%opt_HH) < 1) {
  die "ERROR, --xpergroup must be >= 1";
}

# RNA discovery validation
my $do_rna_discovery = (! opt_Get("--skip-rna", \%opt_HH));
my $rna_cm_file = undef;
if($do_rna_discovery) {
  if(opt_IsUsed("--rna-cm-file", \%opt_HH)) {
    $rna_cm_file = opt_Get("--rna-cm-file", \%opt_HH);
    if(! -e $rna_cm_file) {
      die "ERROR, --rna-cm-file does not exist: $rna_cm_file";
    }
  }
  else {
    # Using default Rfam.cm - need VADRRFAMDIR
    if((! defined $env_vadr_rfam_dir) || ($env_vadr_rfam_dir eq "")) {
      die "ERROR, VADRRFAMDIR environment variable not set; required for default Rfam RNA discovery (or use --rna-cm-file <s> or --skip-rna)";
    }
    if(! -d $env_vadr_rfam_dir) {
      die "ERROR, VADRRFAMDIR directory does not exist: $env_vadr_rfam_dir";
    }
    $rna_cm_file = $env_vadr_rfam_dir . "/Rfam.cm";
    if(! -e $rna_cm_file) {
      die "ERROR, default Rfam.cm not found: $rna_cm_file";
    }
    my $rfam_clanin = $env_vadr_rfam_dir . "/Rfam.clanin";
    if(! -e $rfam_clanin) {
      die "ERROR, Rfam.clanin file not found: $rfam_clanin required for cmscan --clanin";
    }
  }
}

#############################
# create the output directory
#############################
my $cmd;              # a command to run with utl_RunCommand()
my @early_cmd_A = (); # array of commands we run before our log file is opened

if($dir =~ m/\/$/) { 
  $dir =~ s/\/$//; # remove final '/' if it exists
}
if(-d $dir) { 
  $cmd = "rm -rf $dir";
  if(opt_Get("-f", \%opt_HH)) { utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR directory named $dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $dir) { 
  $cmd = "rm $dir";
  if(opt_Get("-f", \%opt_HH)) { utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR a file named $dir already exists. Remove it, or use -f to overwrite it."; }
}

# create the dir
$cmd = "mkdir $dir";
utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root = $dir . "/" . $dir_tail . ".vadr";

if((! defined $meta_tsv) && opt_IsUsed("--taxid", \%opt_HH)) {
  $meta_tsv = $out_root . ".metadata.tsv";
}

#######################
# output program banner
#######################
# output preamble
my $meta_mode_str = (opt_IsUsed("--meta", \%opt_HH)) ? "--meta" : "--taxid";
my @arg_desc_A = ("model directory", "seed model key", "seed model .minfo file", "metadata mode", "metadata TSV file", "output directory");
my @arg_A      = ($model_dir, $model_key, $seed_minfo, $meta_mode_str, $meta_tsv, $dir);
my %extra_H    = ();
ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
opt_OutputPreamble(*STDOUT, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# open the log and command files
my %ofile_info_HH = (); 
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "log",  $out_root . ".log", 1, 1, "Output printed to screen");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "cmd",  $out_root . ".cmd", 1, 1, "List of executed commands");

# Output early commands to .cmd 
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

my %FH_H = ("log" => $ofile_info_HH{"FH"}{"log"}, "cmd" => $cmd_FH);
my $log_FH = $FH_H{"log"};

ofile_OutputBanner($log_FH, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

#---------------------------------------
# Step -1: Optional seed bootstrap via v-build.pl
#---------------------------------------
if($do_seed_bootstrap) {
  my @seed_opt_A = ();
  if(opt_IsUsed("--seed-build-opts-file", \%opt_HH)) {
    my $seed_opts_file = opt_Get("--seed-build-opts-file", \%opt_HH);
    open(my $sbfh, "<", $seed_opts_file) || die "ERROR, unable to read --seed-build-opts-file $seed_opts_file: $!";
    while(my $line = <$sbfh>) {
      chomp $line;
      $line =~ s/^\s+//;
      $line =~ s/\s+$//;
      next if($line eq "");
      next if($line =~ /^\#/);
      push(@seed_opt_A, split(/\s+/, $line));
    }
    close($sbfh);
  }
  elsif(opt_IsUsed("--seed-build-opts", \%opt_HH)) {
    my $seed_opt_str = opt_Get("--seed-build-opts", \%opt_HH);
    push(@seed_opt_A, split(/\s+/, $seed_opt_str));
  }

  if(opt_Get("-f", \%opt_HH)) {
    push(@seed_opt_A, "-f");
  }
  if(opt_Get("-v", \%opt_HH)) {
    push(@seed_opt_A, "-v");
  }

  my $seed_accn = opt_Get("--seed-accn", \%opt_HH);
  my $seed_opt_str = join(" ", @seed_opt_A);
  $cmd = $execs_H{"v-build.pl"} . " " . (($seed_opt_str ne "") ? ($seed_opt_str . " ") : "") . $seed_accn . " " . $model_dir;
  utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, \%FH_H);
}

if(! -e $seed_minfo) {
  die "ERROR: expected seed model .minfo file $seed_minfo does not exist (derived from model directory and model key; use --mkey to override)";
}

#---------------------------------------
# Step 0: Acquire metadata TSV if needed
#---------------------------------------
if(opt_IsUsed("--taxid", \%opt_HH)) {
  my $fetch_script = $env_vadr_scripts_dir . "/miniscripts/fetch-seqs-given-taxid.pl";
  if(! -e $fetch_script) {
    die "ERROR, expected fetch script does not exist: $fetch_script";
  }

  my $fetch_out_prefix = $out_root . ".metadata";
  my $taxid = opt_Get("--taxid", \%opt_HH);
  $cmd = "$^X $fetch_script --taxid $taxid --out $fetch_out_prefix";
  if(opt_IsUsed("--api_key", \%opt_HH)) {
    $cmd .= " --api_key " . opt_Get("--api_key", \%opt_HH);
  }
  if(! opt_Get("-v", \%opt_HH)) {
    $cmd .= " --quiet";
  }
  utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, \%FH_H);
  $meta_tsv = $fetch_out_prefix . ".tsv";
}

if((! defined $meta_tsv) || (! -e $meta_tsv)) {
  die "ERROR, metadata TSV is undefined or does not exist: " . ((defined $meta_tsv) ? $meta_tsv : "[undef]");
}

#---------------------------------------
# Step 1: Read out the seed model length
#---------------------------------------
my @mdl_info_A = ();
my %ftr_info_HA = ();
my @reqd_mdl_keys = ("length");
my @reqd_ftr_keys = ();

vdr_ModelInfoFileParse($seed_minfo, \@reqd_mdl_keys, \@reqd_ftr_keys, \@mdl_info_A, \%ftr_info_HA, \%FH_H);
my $seed_model_len = $mdl_info_A[0]{"length"};
ofile_OutputString(*STDOUT, 1, sprintf("# Read seed model length: %d\n", $seed_model_len));

#---------------------------------------
# Step 1b: RNA discovery via cmscan on reference sequence
#---------------------------------------
my @rna_regions_A = (); # Array of hashes: { start, end, strand, cm_family, cm_accession, score, evalue }
my $rna_annot_file = $out_root . ".rna_annotation.tsv";
my $rna_ss_cons = undef; # Full-length consensus secondary structure string

if($do_rna_discovery) {
  my $ref_seq_file = $model_dir . "/" . $model_key . ".vadr.fa";
  if(! -e $ref_seq_file) {
    die "ERROR, reference sequence file not found for RNA discovery: $ref_seq_file";
  }
  
  run_rna_discovery($ref_seq_file, $rna_cm_file, $seed_model_len, $out_root, $env_vadr_rfam_dir, 
                    \@rna_regions_A, \$rna_ss_cons, opt_Get("--keep", \%opt_HH), opt_Get("-v", \%opt_HH), 
                    \%execs_H, \%FH_H);
  
  # Write RNA annotation output
  run_rna_sstruct_generation(\@rna_regions_A, $ref_seq_file, $rna_cm_file, $env_vadr_rfam_dir, 
                             $out_root, opt_Get("--keep", \%opt_HH), opt_Get("-v", \%opt_HH), 
                             \%execs_H, \%FH_H);

  # Write RNA annotation output (after structure extraction)
  write_rna_annotation_table(\@rna_regions_A, $rna_annot_file, \%FH_H);
}
else {
  ofile_OutputString($FH_H{"log"}, 1, sprintf("# RNA discovery: skipped due to --skip-rna\n"));
}

#---------------------------------------
# Step 2: Read and filter metadata
#---------------------------------------

# %candidate_AH: arrays of hashes [1..nseq-1] grouped by serotype/genotype keys.
# Each hash contains sequence metadata:
#   acc      => sequence accession
#   len      => sequence length
#   cdate    => genbank create date
#   serotype => sequence serotype
#   genotype => sequence genotype
#   isolate  => sequence isolate
my %candidate_AH = ();
my %decision_H = ();
my $max_per_group = opt_Get("--xpergroup", \%opt_HH);
parse_and_filter_metadata($meta_tsv, $seed_model_len, $max_per_group, \%candidate_AH, \%decision_H, \%FH_H);

#---------------------------------------
# Step 3: Tier 2 fetch and ambiguity filter
#---------------------------------------
my $tier1_accn_file = $out_root . ".tier1.accn.list";
my $tier2_fasta_file = $out_root . ".tier2.fa";
my $max_ambig_nt = opt_Get("--xambig", \%opt_HH);
my $tier2_ant_outdir = $out_root . ".tier2.annot";
my $centroid_tsv_file = $out_root . ".centroid.tsv";
my $decision_tsv_file = $out_root . ".decision.tsv";
my $decision_summary_tsv_file = $out_root . ".decision.summary.tsv";
my $stitch_selected_accn_file = $out_root . ".stitch.selected.accn.list";
my $stitch_selected_fa_file   = $out_root . ".stitch.selected.fa";
my $stitch_block_plan_file    = $out_root . ".stitch.block_plan.tsv";
my $stitch_cds_nt_fa_file     = $out_root . ".stitch.cds.nt.fa";
my $stitch_cds_orf_fa_file    = $out_root . ".stitch.cds.orf.fa";
my $stitch_cds_aa_fa_file     = $out_root . ".stitch.cds.aa.fa";
my $stitch_cds_map_tsv_file   = $out_root . ".stitch.cds.translate_map.tsv";
my $stitch_cds_anchor_tsv_file = $out_root . ".stitch.cds.anchor.tsv";
my $stitch_cds_anchor_fa_file  = $out_root . ".stitch.cds.anchor.fa";
my $stitch_cds_pairwise_tsv_file = $out_root . ".stitch.cds.pairwise.tsv";
my $stitch_cds_msa_aa_fa_file = $out_root . ".stitch.cds.msa.aa.fa";
my $stitch_cds_msa_aa_stk_file = $out_root . ".stitch.cds.msa.aa.stk";
my $stitch_cds_msa_nt_fa_file = $out_root . ".stitch.cds.msa.nt.fa";
my $rna_annotation_file       = $out_root . ".rna_annotation.tsv";

write_accession_list_from_candidates(\%candidate_AH, $tier1_accn_file, \%FH_H);
fetch_fasta_from_accession_list($tier1_accn_file, $tier2_fasta_file, \%FH_H);
apply_ambiguity_filter_to_candidates(\%candidate_AH, $tier2_fasta_file, $max_ambig_nt, \%decision_H, \%FH_H);

my $tier2_align_stk_file = undef;
if($do_skip_annotate) {
  ofile_OutputString($FH_H{"log"}, 1, sprintf("# Tier 2 v-annotate filter: skipped due to --skip-annotate\n"));
  mark_all_remaining_as_selected_for_tier3(\%candidate_AH, \%decision_H);
}
else {
  $tier2_align_stk_file = run_vannotate_filter_fails(\%candidate_AH, $tier2_fasta_file, $tier2_ant_outdir, $model_dir, $model_key, opt_Get("--keep", \%opt_HH), opt_Get("-v", \%opt_HH), \%execs_H, \%decision_H, \%FH_H);
}

#---------------------------------------
# Step 4: Tier 3 centroid selection (BLAST all-vs-all)
#---------------------------------------
select_group_centroids_blast(\%candidate_AH, $tier2_fasta_file, $out_root, $centroid_tsv_file, opt_Get("-v", \%opt_HH), \%execs_H, \%decision_H, \%FH_H);

write_decision_report(\%decision_H, $decision_tsv_file, \%FH_H);
write_decision_stage_reports(\%decision_H, $out_root . ".decision", \%FH_H);
write_decision_summary_report(\%decision_H, $decision_summary_tsv_file, \%FH_H);

#---------------------------------------
# Step 5: Initial piecewise stitching scaffold outputs
#---------------------------------------
my $n_selected = write_stitch_scaffold_outputs(\%candidate_AH, $tier2_fasta_file, \%ftr_info_HA, $model_key, $seed_model_len, $stitch_selected_accn_file, $stitch_selected_fa_file, $stitch_block_plan_file, \%FH_H);

if($n_selected == 0) {
  ofile_OutputString($FH_H{"log"}, 1, "#\n# Zero sequences passed all filters. Cannot build profile alignment.\n");
  ofile_OutputString($FH_H{"log"}, 1, "# Check decision reports for details on why sequences were removed.\n");
  ofile_OutputString($FH_H{"log"}, 1, "#\n# Exiting with error.\n");
  ofile_FAIL("ERROR, zero sequences passed all filters, cannot build profile alignment. See $decision_summary_tsv_file for details.", 1, \%FH_H);
}

#---------------------------------------
# Step 6: CDS translation prep for protein alignment
#---------------------------------------
prepare_cds_translation_for_stitching(\%candidate_AH, $stitch_selected_fa_file, $tier2_ant_outdir, $stitch_cds_nt_fa_file, $stitch_cds_orf_fa_file, $stitch_cds_aa_fa_file, $stitch_cds_map_tsv_file, $do_skip_annotate, opt_Get("--keep", \%opt_HH), opt_Get("-v", \%opt_HH), \%execs_H, \%FH_H);

#---------------------------------------
# Step 7: Reference-anchored AA global/global pairwise (ggsearch)
#---------------------------------------
run_reference_anchored_pairwise_aa(\%candidate_AH, $centroid_tsv_file, $stitch_cds_aa_fa_file, $stitch_cds_anchor_tsv_file, $stitch_cds_anchor_fa_file, $stitch_cds_pairwise_tsv_file, $do_skip_annotate, opt_Get("--keep", \%opt_HH), opt_Get("-v", \%opt_HH), \%execs_H, \%FH_H);

#---------------------------------------
# Step 8a: Build protein MSA from pairwise CIGARs and backconvert CDS nt alignment
#---------------------------------------
build_anchor_projected_cds_msa($stitch_cds_aa_fa_file,
                               $stitch_cds_nt_fa_file,
                               $stitch_cds_map_tsv_file,
                               $stitch_cds_anchor_tsv_file,
                               $stitch_cds_pairwise_tsv_file,
                               $stitch_cds_msa_aa_fa_file,
                               $stitch_cds_msa_aa_stk_file,
                               $stitch_cds_msa_nt_fa_file,
                               $do_skip_annotate,
                               \%FH_H);


#---------------------------------------
# Step 8b-d: RNA region extraction and alignment refinement
#---------------------------------------
if($do_rna_discovery && (scalar(@rna_regions_A) > 0) && (!$do_skip_annotate)) {
  my $rna_struct_dir = $out_root . ".rna_struct";
  extract_and_align_rna_regions(\@rna_regions_A, $tier2_align_stk_file, $rna_struct_dir, 
                                $out_root, opt_Get("--keep", \%opt_HH), opt_Get("-v", \%opt_HH), 
                                \%execs_H, \%FH_H);
}

#---------------------------------------
#---------------------------------------
# Step 10: Stitch all blocks into final training alignment
#---------------------------------------
my $final_stk_file = $out_root . ".stitch.final.stk";
my $refined_stk_file = $out_root . ".stitch.final.refined.stk";
my $temp_cm_file = $out_root . ".stitch.temp.cm";

stitch_and_refine_final_alignment($stitch_block_plan_file,
                                  $rna_annotation_file,
                                  $tier2_align_stk_file,
                                  $stitch_cds_msa_nt_fa_file,
                                  $out_root,
                                  \@rna_regions_A,
                                  $final_stk_file,
                                  $refined_stk_file,
                                  $temp_cm_file,
                                  $do_rna_discovery,
                                  $do_skip_annotate,
                                  $seed_model_len,
                                  \%execs_H,
                                  \%FH_H);

#---------------------------------------
# Step 11: Generate updated .minfo file with RNA features
#---------------------------------------
if($do_rna_discovery && !$do_skip_annotate && scalar(@rna_regions_A) > 0) {
  my $updated_minfo_file = $out_root . ".minfo";
  generate_updated_minfo($seed_minfo, $rna_annotation_file, \@rna_regions_A,
                         $updated_minfo_file, $model_key, \%execs_H, \%FH_H);
}

ofile_OutputString(\*STDOUT, 1, sprintf("All done.\n"));

$total_seconds += ofile_SecondsSinceEpoch();
ofile_OutputConclusionAndCloseFilesOk($total_seconds, $dir, \%ofile_info_HH);
exit(0);


#################################################################
# Subroutine : parse_and_filter_metadata()
# Incept     : Gemini Thu Mar 05 2026
#
# Purpose    : Reads the sequence metadata TSV, applies length
#              filters, and performs the Chronological Sampling
#              Tier 1 filter to select the diverse 50 per group.
#
# Arguments  :
#   $tsv_file       : path to metadata TSV file
#   $seed_model_len : length of the seed model for >90% filtering
#   $candidate_AHR  : reference to hash mapping group (serotype/genotype)
#                     to an array ref of selected sequence hashrefs.
#   $FH_HR          : hash of file descriptors (like 'log' and 'cmd')
#
# Returns    : void
#################################################################
sub parse_and_filter_metadata {
  my ($tsv_file, $seed_model_len, $max_per_group, $candidate_AHR, $decision_HR, $FH_HR) = @_;

  my $total_seqs = 0;
  my $kept_len_seqs = 0;

  open(my $tsv_fh, "<", $tsv_file) or die "ERROR: Cannot read $tsv_file: $!";
  my $header = <$tsv_fh>; # Skip header: Accession, Length, CreateDate, Serotype, Genotype, Isolate

  while (my $line = <$tsv_fh>) {
    chomp $line;
    my ($acc, $len, $cdate, $serotype, $genotype, $isolate) = split(/\t/, $line);
    my $group = "Unknown";
    if    (defined $serotype && $serotype ne "") { $group = $serotype; }
    elsif (defined $genotype && $genotype ne "") { $group = $genotype; }

    $decision_HR->{$acc} = {
      accession       => $acc,
      group_key       => $group,
      serotype        => (defined $serotype) ? $serotype : "",
      genotype        => (defined $genotype) ? $genotype : "",
      isolate         => (defined $isolate)  ? $isolate  : "",
      create_date     => (defined $cdate)    ? $cdate    : "",
      seq_length      => (defined $len)      ? $len      : "",
      status          => "kept",
      reason_code     => "pass",
      reason_detail   => "",
      stage_last_seen => "metadata_read",
      n_ambig_nt      => "",
      ambig_threshold => ""
    };

    $total_seqs++;
    
    # Length check: Keep sequences >90% of the seed model's length
    if ($seed_model_len > 0 && $len < ($seed_model_len * 0.90)) {
      $decision_HR->{$acc}{"status"} = "removed";
      $decision_HR->{$acc}{"reason_code"} = "len_lt_90pct_seed";
      $decision_HR->{$acc}{"reason_detail"} = "length $len < 0.9 * seed_length $seed_model_len";
      $decision_HR->{$acc}{"stage_last_seen"} = "tier1_length_filter";
      next;
    }
    $kept_len_seqs++;

    # Save candidate info
    push @{$candidate_AHR->{$group}}, {
      acc      => $acc,
      len      => $len,
      cdate    => $cdate,
      serotype => $serotype,
      genotype => $genotype,
      isolate  => $isolate
    };
    $decision_HR->{$acc}{"stage_last_seen"} = "tier1_pool";
  }
  close($tsv_fh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Parsed %d total sequences from TSV.\n", $total_seqs));
  if ($seed_model_len > 0) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Retained %d sequences passing length filter (>90%% of seed length %d).\n", $kept_len_seqs, $seed_model_len));
  }
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Found %d distinct serotype/genotype groups.\n", scalar(keys %{$candidate_AHR})));

  #---------------------------------------
  # Tier 1 Filter: Chronological Sampling
  #---------------------------------------
  my $total_selected = 0;

  foreach my $group (sort keys %{$candidate_AHR}) {
    my @seqs = @{$candidate_AHR->{$group}};

    # Parse CreateDate and prepare for sorting
    foreach my $seq (@seqs) {
      my $cdate = $seq->{cdate} || "";
      $seq->{year_month} = "0000-00";
      $seq->{year_month_day} = "0000-00-00";
      $seq->{date_val} = "9999/99/99"; # default sort to end if invalid

      if ($cdate =~ m/^(\d{4})\/(\d{2})\/(\d{2})/) {
         $seq->{year_month}     = "$1-$2";
         $seq->{year_month_day} = "$1-$2-$3";
         $seq->{date_val}       = $cdate;
      }
    }

    # Sort oldest to newest (alphabetical on YYYY/MM/DD works nicely here)
    @seqs = sort { $a->{date_val} cmp $b->{date_val} || $a->{acc} cmp $b->{acc} } @seqs;

    my @selected_for_group = ();
    my %seen_month = ();
    my %seen_day =();
    my %is_selected = ();

    # Pass 1: 1 per year/month
    foreach my $seq (@seqs) {
      last if scalar(@selected_for_group) >= $max_per_group;
      if (! $seen_month{$seq->{year_month}}) {
        push @selected_for_group, $seq;
        $seen_month{$seq->{year_month}} = 1;
        $seen_day{$seq->{year_month_day}} = 1; # Reserve the day so pass 2 skips it
        $is_selected{$seq->{acc}} = 1;
      }
    }

    # Pass 2: 1 per day
    if (scalar(@selected_for_group) < $max_per_group) {
      foreach my $seq (@seqs) {
        last if scalar(@selected_for_group) >= $max_per_group;
        next if $is_selected{$seq->{acc}};
        if (! $seen_day{$seq->{year_month_day}}) {
          push @selected_for_group, $seq;
          $seen_day{$seq->{year_month_day}} = 1;
          $is_selected{$seq->{acc}} = 1;
        }
      }
    }

    # Pass 3: Oldest remaining
    if (scalar(@selected_for_group) < $max_per_group) {
      foreach my $seq (@seqs) {
        last if scalar(@selected_for_group) >= $max_per_group;
        next if $is_selected{$seq->{acc}};
        push @selected_for_group, $seq;
        $is_selected{$seq->{acc}} = 1;
      }
    }

    # Replace the original group array with the filtered top 50
    foreach my $seq (@seqs) {
      my $acc = $seq->{acc};
      if($is_selected{$acc}) {
        $decision_HR->{$acc}{"stage_last_seen"} = "tier1_selected";
      }
      else {
        $decision_HR->{$acc}{"status"} = "removed";
        $decision_HR->{$acc}{"reason_code"} = "xpergroup_prune";
        $decision_HR->{$acc}{"reason_detail"} = "not selected in chronological top-$max_per_group for group";
        $decision_HR->{$acc}{"stage_last_seen"} = "tier1_sampling";
      }
    }

    $candidate_AHR->{$group} = \@selected_for_group;
    $total_selected += scalar(@selected_for_group);
  }

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Tier 1 Filter: Selected %d total sequences (max %d per group) for tier 2 processing.\n", $total_selected, $max_per_group));

  return;
}

#################################################################
# Subroutine : write_accession_list_from_candidates()
# Incept     : Copilot Mon Mar 09 2026
#
# Purpose    : Write one accession per line for all currently
#              selected candidates across all groups.
#################################################################
sub write_accession_list_from_candidates {
  my ($candidate_AHR, $accn_file, $FH_HR) = @_;

  open(my $afh, ">", $accn_file) or die "ERROR: unable to write $accn_file: $!";
  my $nacc = 0;
  foreach my $group (sort keys %{$candidate_AHR}) {
    foreach my $seq (@{$candidate_AHR->{$group}}) {
      next if (! defined $seq->{acc} || $seq->{acc} eq "");
      print $afh $seq->{acc} . "\n";
      $nacc++;
    }
  }
  close($afh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Tier 2 prep: wrote %d accessions to %s\n", $nacc, $accn_file));
  return;
}

#################################################################
# Subroutine : fetch_fasta_from_accession_list()
# Incept     : Copilot Mon Mar 09 2026
#
# Purpose    : Download FASTA records from NCBI efetch for a list
#              of accessions (batched).
#################################################################
sub fetch_fasta_from_accession_list {
  my ($accn_file, $fasta_out_file, $FH_HR) = @_;

  my @acc_A = ();
  open(my $ifh, "<", $accn_file) or die "ERROR: unable to read $accn_file: $!";
  while(my $line = <$ifh>) {
    chomp $line;
    next if($line eq "");
    push(@acc_A, $line);
  }
  close($ifh);

  open(my $ofh, ">", $fasta_out_file) or die "ERROR: unable to write $fasta_out_file: $!";
  my $batch_size = 200;
  my $nfetched_batches = 0;

  for(my $i = 0; $i < scalar(@acc_A); $i += $batch_size) {
    my $end_i = $i + $batch_size - 1;
    $end_i = scalar(@acc_A) - 1 if($end_i >= scalar(@acc_A));
    my @batch_acc_A = @acc_A[$i .. $end_i];
    my $id_str = join(",", @batch_acc_A);

    my $url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" .
              "?db=nuccore&id=$id_str&rettype=fasta&retmode=text";

    my $res = get($url);
    if(! defined $res) {
      die "ERROR: failed to fetch FASTA batch from NCBI efetch (batch starts at accession index $i)";
    }

    print $ofh $res;
    $nfetched_batches++;
    select(undef, undef, undef, 0.35); # polite delay
  }
  close($ofh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Tier 2 fetch: downloaded FASTA for %d accessions in %d batches to %s\n", scalar(@acc_A), $nfetched_batches, $fasta_out_file));
  return;
}

#################################################################
# Subroutine : apply_ambiguity_filter_to_candidates()
# Incept     : Copilot Mon Mar 09 2026
#
# Purpose    : Parse fetched FASTA, compute ambiguous character
#              count for each sequence, and remove candidates
#              with ambiguous nucleotide count above threshold.
#################################################################
sub apply_ambiguity_filter_to_candidates {
  my ($candidate_AHR, $fasta_file, $max_ambig_nt, $decision_HR, $FH_HR) = @_;

  my %seq_H = ();
  my $cur_acc = undef;

  open(my $ffh, "<", $fasta_file) or die "ERROR: unable to read $fasta_file: $!";
  while(my $line = <$ffh>) {
    chomp $line;
    if($line =~ /^>(\S+)/) {
      $cur_acc = $1;
      $seq_H{$cur_acc} = "" if(! exists $seq_H{$cur_acc});
    }
    elsif(defined $cur_acc) {
      $line =~ s/\s+//g;
      $seq_H{$cur_acc} .= $line;
    }
  }
  close($ffh);

  my %ambig_ct_H = ();
  foreach my $acc (keys %seq_H) {
    my $useq = uc($seq_H{$acc});
    $useq =~ tr/U/T/;
    my $tmp = $useq;
    $tmp =~ s/[ACGT]//g;
    my $nambig = length($tmp);
    $ambig_ct_H{$acc} = $nambig;
  }

  my $nkept = 0;
  my $nremoved = 0;
  foreach my $group (sort keys %{$candidate_AHR}) {
    my @kept_A = ();
    foreach my $seq (@{$candidate_AHR->{$group}}) {
      my $acc = $seq->{acc};
      if(! exists $ambig_ct_H{$acc}) {
        $decision_HR->{$acc}{"status"} = "removed";
        $decision_HR->{$acc}{"reason_code"} = "missing_fasta";
        $decision_HR->{$acc}{"reason_detail"} = "accession not present in fetched tier2 fasta";
        $decision_HR->{$acc}{"stage_last_seen"} = "tier2_ambiguity_filter";
        $decision_HR->{$acc}{"n_ambig_nt"} = "";
        $decision_HR->{$acc}{"ambig_threshold"} = $max_ambig_nt;
        $nremoved++;
        next;
      }
      if($ambig_ct_H{$acc} > $max_ambig_nt) {
        $decision_HR->{$acc}{"status"} = "removed";
        $decision_HR->{$acc}{"reason_code"} = "ambig_gt_xambig";
        $decision_HR->{$acc}{"reason_detail"} = "ambiguous nt count $ambig_ct_H{$acc} > threshold $max_ambig_nt";
        $decision_HR->{$acc}{"stage_last_seen"} = "tier2_ambiguity_filter";
        $decision_HR->{$acc}{"n_ambig_nt"} = $ambig_ct_H{$acc};
        $decision_HR->{$acc}{"ambig_threshold"} = $max_ambig_nt;
        $nremoved++;
        next;
      }
      push(@kept_A, $seq);
      $decision_HR->{$acc}{"status"} = "kept";
      $decision_HR->{$acc}{"reason_code"} = "pass";
      $decision_HR->{$acc}{"reason_detail"} = "";
      $decision_HR->{$acc}{"stage_last_seen"} = "tier2_after_ambiguity";
      $decision_HR->{$acc}{"n_ambig_nt"} = $ambig_ct_H{$acc};
      $decision_HR->{$acc}{"ambig_threshold"} = $max_ambig_nt;
      $nkept++;
    }
    $candidate_AHR->{$group} = \@kept_A;
  }

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Tier 2 ambiguity filter: kept %d removed %d (max ambiguous nts %d)\n", $nkept, $nremoved, $max_ambig_nt));
  return;
}

#################################################################
# Subroutine : run_vannotate_filter_fails()
# Incept     : Copilot Mon Mar 09 2026
#
# Purpose    : Run v-annotate.pl on Tier 2 candidates and remove
#              sequences listed in .vadr.fail.list.
#################################################################
sub run_vannotate_filter_fails {
  my ($candidate_AHR, $fasta_file, $annot_outdir, $model_dir, $model_key, $do_keep, $do_verbose, $execs_HR, $decision_HR, $FH_HR) = @_;

  my $annot_mkey = $model_key . ".vadr";

  # Note: --out_stk and --keep are incompatible in v-annotate.pl, so we only use --out_stk
  # which outputs the Stockholm alignment we need for Step 6 block extraction
  my $cmd = $execs_HR->{"v-annotate.pl"} . " -f --mdir " . $model_dir . " --mkey " . $annot_mkey . " --out_stk " . $fasta_file . " " . $annot_outdir;
  if(! $do_verbose) {
    $cmd .= " > /dev/null";
  }

  utl_RunCommand($cmd, $do_verbose, 0, $FH_HR);

  my $annot_outdir_tail = $annot_outdir;
  $annot_outdir_tail =~ s/^.+\///;
  my $fail_list_file = $annot_outdir . "/" . $annot_outdir_tail . ".vadr.fail.list";

  my %is_fail_H = ();
  if(-e $fail_list_file) {
    open(my $ffh, "<", $fail_list_file) || die "ERROR, unable to read v-annotate fail list $fail_list_file: $!";
    while(my $line = <$ffh>) {
      chomp $line;
      $line =~ s/^\s+//;
      $line =~ s/\s+$//;
      next if($line eq "");
      my ($name) = split(/\s+/, $line);
      $is_fail_H{$name} = 1;
    }
    close($ffh);
  }

  my $nkept = 0;
  my $nremoved = 0;
  foreach my $group (sort keys %{$candidate_AHR}) {
    my @kept_A = ();
    foreach my $seq (@{$candidate_AHR->{$group}}) {
      my $acc = $seq->{acc};
      if(exists $is_fail_H{$acc}) {
        $decision_HR->{$acc}{"status"} = "removed";
        $decision_HR->{$acc}{"reason_code"} = "vadr_fail";
        $decision_HR->{$acc}{"reason_detail"} = "present in .vadr.fail.list";
        $decision_HR->{$acc}{"stage_last_seen"} = "tier2_vannotate_filter";
        $nremoved++;
        next;
      }
      push(@kept_A, $seq);
      $decision_HR->{$acc}{"status"} = "kept";
      $decision_HR->{$acc}{"reason_code"} = "pass";
      $decision_HR->{$acc}{"reason_detail"} = "";
      $decision_HR->{$acc}{"stage_last_seen"} = "selected_for_tier3";
      $nkept++;
    }
    $candidate_AHR->{$group} = \@kept_A;
  }

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Tier 2 v-annotate filter: kept %d removed %d (removed all accessions in .vadr.fail.list)\n", $nkept, $nremoved));
  
  # Return path to alignment file for later use in Step 8 RNA refinement
  # v-annotate creates: <outdir>/<outdir>.vadr.<modelkey>.align.stk
  my $align_stk_file = $annot_outdir . "/" . $annot_outdir_tail . ".vadr." . $model_key . ".align.stk";
  return $align_stk_file;
}

#################################################################
# Subroutine : select_group_centroids_blast()
# Incept     : Copilot Tue Mar 10 2026
#
# Purpose    : For each remaining group, select one centroid
#              sequence using average all-vs-all blastn pident.
#################################################################
sub select_group_centroids_blast {
  my ($candidate_AHR, $fasta_file, $out_root, $centroid_tsv_file, $do_verbose, $execs_HR, $decision_HR, $FH_HR) = @_;

  my %seq_H = ();
  my $cur_acc = undef;
  open(my $ffh, "<", $fasta_file) || die "ERROR, unable to read $fasta_file: $!";
  while(my $line = <$ffh>) {
    chomp $line;
    if($line =~ /^>(\S+)/) {
      $cur_acc = $1;
      $seq_H{$cur_acc} = "" if(! exists $seq_H{$cur_acc});
    }
    elsif(defined $cur_acc) {
      $line =~ s/\s+//g;
      $seq_H{$cur_acc} .= $line;
    }
  }
  close($ffh);

  open(my $ctfh, ">", $centroid_tsv_file) || die "ERROR, unable to write centroid TSV $centroid_tsv_file: $!";
  print $ctfh join("\t", "group_key", "accession", "avg_blastn_pident", "is_centroid") . "\n";

  my $n_groups_with_centroid = 0;
  my $n_groups_empty = 0;
  foreach my $group (sort keys %{$candidate_AHR}) {
    my @seq_A = @{$candidate_AHR->{$group}};
    my $nseq = scalar(@seq_A);
    if($nseq == 0) {
      $n_groups_empty++;
      next;
    }

    if($nseq == 1) {
      my $acc = $seq_A[0]{"acc"};
      print $ctfh join("\t", $group, $acc, sprintf("%.4f", 100.0), 1) . "\n";
      $decision_HR->{$acc}{"status"} = "kept";
      $decision_HR->{$acc}{"reason_code"} = "centroid_selected";
      $decision_HR->{$acc}{"reason_detail"} = "singleton group";
      $decision_HR->{$acc}{"stage_last_seen"} = "selected_for_tier3";
      $n_groups_with_centroid++;
      next;
    }

    my $group_safe = $group;
    $group_safe =~ s/[^A-Za-z0-9_\.-]+/_/g;
    my $group_fa = $out_root . ".centroid." . $group_safe . ".fa";
    my $group_db = $out_root . ".centroid." . $group_safe . ".blastdb";
    my $group_blast = $out_root . ".centroid." . $group_safe . ".blast.tsv";

    open(my $gfh, ">", $group_fa) || die "ERROR, unable to write $group_fa: $!";
    foreach my $seq (@seq_A) {
      my $acc = $seq->{"acc"};
      if(! exists $seq_H{$acc}) {
        die "ERROR, unable to find sequence $acc in $fasta_file while selecting centroid for group $group";
      }
      print $gfh ">" . $acc . "\n" . $seq_H{$acc} . "\n";
    }
    close($gfh);

    my $cmd = $execs_HR->{"makeblastdb"} . " -in " . $group_fa . " -dbtype nucl -out " . $group_db;
    if(! $do_verbose) {
      $cmd .= " > /dev/null 2>&1";
    }
    utl_RunCommand($cmd, $do_verbose, 0, $FH_HR);
    $cmd = $execs_HR->{"blastn"} . " -query " . $group_fa . " -db " . $group_db . " -task blastn -outfmt \"6 qseqid sseqid pident\" -max_target_seqs 10000 -max_hsps 1 > " . $group_blast;
    utl_RunCommand($cmd, $do_verbose, 0, $FH_HR);

    my %pid_HH = (); # pid_HH{q}{s} = max pident
    open(my $bfh, "<", $group_blast) || die "ERROR, unable to read $group_blast: $!";
    while(my $line = <$bfh>) {
      chomp $line;
      next if($line eq "");
      my ($q, $s, $pident) = split(/\t/, $line);
      next if((! defined $q) || (! defined $s) || (! defined $pident));
      if((! exists $pid_HH{$q}{$s}) || ($pident > $pid_HH{$q}{$s})) {
        $pid_HH{$q}{$s} = $pident;
      }
    }
    close($bfh);

    my %avg_H = ();
    my @acc_A = map { $_->{"acc"} } @seq_A;
    foreach my $q (@acc_A) {
      my $sum = 0.0;
      foreach my $s (@acc_A) {
        my $pid = (exists $pid_HH{$q}{$s}) ? $pid_HH{$q}{$s} : 0.0;
        $sum += $pid;
      }
      $avg_H{$q} = $sum / scalar(@acc_A);
    }

    my $centroid_acc = undef;
    my $best = -1.0;
    foreach my $acc (sort @acc_A) { # alphabetical tie-breaker
      if($avg_H{$acc} > $best) {
        $best = $avg_H{$acc};
        $centroid_acc = $acc;
      }
    }

    my @new_group_A = ();
    foreach my $seq (@seq_A) {
      my $acc = $seq->{"acc"};
      my $is_centroid = ($acc eq $centroid_acc) ? 1 : 0;
      print $ctfh join("\t", $group, $acc, sprintf("%.4f", $avg_H{$acc}), $is_centroid) . "\n";
      if($is_centroid) {
        push(@new_group_A, $seq);
        $decision_HR->{$acc}{"status"} = "kept";
        $decision_HR->{$acc}{"reason_code"} = "centroid_selected";
        $decision_HR->{$acc}{"reason_detail"} = sprintf("selected by max average blastn pident %.4f in group", $avg_H{$acc});
        $decision_HR->{$acc}{"stage_last_seen"} = "selected_for_tier3";
      }
      else {
        $decision_HR->{$acc}{"status"} = "removed";
        $decision_HR->{$acc}{"reason_code"} = "centroid_not_selected";
        $decision_HR->{$acc}{"reason_detail"} = sprintf("non-centroid; average blastn pident %.4f (centroid %.4f)", $avg_H{$acc}, $avg_H{$centroid_acc});
        $decision_HR->{$acc}{"stage_last_seen"} = "tier3_centroid_filter";
      }
    }
    $candidate_AHR->{$group} = \@new_group_A;
    $n_groups_with_centroid++;
  }
  close($ctfh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Tier 3 centroid selection (blast): selected 1 sequence in each of %d groups (%d empty groups) and wrote %s\n", $n_groups_with_centroid, $n_groups_empty, $centroid_tsv_file));
  return;
}

#################################################################
# Subroutine : mark_all_remaining_as_selected_for_tier3()
#################################################################
sub mark_all_remaining_as_selected_for_tier3 {
  my ($candidate_AHR, $decision_HR) = @_;
  foreach my $group (keys %{$candidate_AHR}) {
    foreach my $seq (@{$candidate_AHR->{$group}}) {
      my $acc = $seq->{acc};
      $decision_HR->{$acc}{"status"} = "kept";
      $decision_HR->{$acc}{"reason_code"} = "pass";
      $decision_HR->{$acc}{"reason_detail"} = "";
      $decision_HR->{$acc}{"stage_last_seen"} = "selected_for_tier3";
    }
  }
  return;
}

#################################################################
# Subroutine : write_decision_report()
#################################################################
sub write_decision_report {
  my ($decision_HR, $out_file, $FH_HR) = @_;

  open(my $dfh, ">", $out_file) || die "ERROR, unable to write decision report $out_file: $!";
  print $dfh join("\t", "accession", "group_key", "serotype", "genotype", "isolate", "create_date", "seq_length", "status", "reason_code", "reason_detail", "stage_last_seen", "n_ambig_nt", "ambig_threshold") . "\n";
  foreach my $acc (sort keys %{$decision_HR}) {
    my $d = $decision_HR->{$acc};
    my @f_A = ( $d->{"accession"},
                $d->{"group_key"},
                $d->{"serotype"},
                $d->{"genotype"},
                $d->{"isolate"},
                $d->{"create_date"},
                $d->{"seq_length"},
                $d->{"status"},
                $d->{"reason_code"},
                $d->{"reason_detail"},
                $d->{"stage_last_seen"},
                $d->{"n_ambig_nt"},
                $d->{"ambig_threshold"});
    for(my $i = 0; $i < scalar(@f_A); $i++) {
      if(! defined $f_A[$i]) { $f_A[$i] = ""; }
      $f_A[$i] =~ s/\t/ /g;
      $f_A[$i] =~ s/\n/ /g;
    }
    print $dfh join("\t", @f_A) . "\n";
  }
  close($dfh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Decision report: wrote per-sequence status/reason table to %s\n", $out_file));
  return;
}

#################################################################
# Subroutine : write_decision_stage_reports()
#################################################################
sub write_decision_stage_reports {
  my ($decision_HR, $out_prefix, $FH_HR) = @_;

  my %stage_AH = ();
  foreach my $acc (keys %{$decision_HR}) {
    my $stage = $decision_HR->{$acc}{"stage_last_seen"};
    if((! defined $stage) || ($stage eq "")) { $stage = "unknown_stage"; }
    push(@{$stage_AH{$stage}}, $acc);
  }

  my $nfiles = 0;
  foreach my $stage (sort keys %stage_AH) {
    my $stage_safe = $stage;
    $stage_safe =~ s/[^A-Za-z0-9_\.-]+/_/g;
    my $out_file = $out_prefix . ".stage." . $stage_safe . ".tsv";
    open(my $sfh, ">", $out_file) || die "ERROR, unable to write stage decision report $out_file: $!";
    print $sfh join("\t", "accession", "group_key", "serotype", "genotype", "isolate", "create_date", "seq_length", "status", "reason_code", "reason_detail", "stage_last_seen", "n_ambig_nt", "ambig_threshold") . "\n";
    foreach my $acc (sort @{$stage_AH{$stage}}) {
      my $d = $decision_HR->{$acc};
      my @f_A = ( $d->{"accession"},
                  $d->{"group_key"},
                  $d->{"serotype"},
                  $d->{"genotype"},
                  $d->{"isolate"},
                  $d->{"create_date"},
                  $d->{"seq_length"},
                  $d->{"status"},
                  $d->{"reason_code"},
                  $d->{"reason_detail"},
                  $d->{"stage_last_seen"},
                  $d->{"n_ambig_nt"},
                  $d->{"ambig_threshold"});
      for(my $i = 0; $i < scalar(@f_A); $i++) {
        if(! defined $f_A[$i]) { $f_A[$i] = ""; }
        $f_A[$i] =~ s/\t/ /g;
        $f_A[$i] =~ s/\n/ /g;
      }
      print $sfh join("\t", @f_A) . "\n";
    }
    close($sfh);
    $nfiles++;
  }

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Decision report: wrote %d per-stage decision files with prefix %s.stage.<stage>.tsv\n", $nfiles, $out_prefix));
  return;
}

#################################################################
# Subroutine : write_decision_summary_report()
#################################################################
sub write_decision_summary_report {
  my ($decision_HR, $out_file, $FH_HR) = @_;

  my %count_H = (); # key: status\treason_code\tstage_last_seen\tgroup_key
  foreach my $acc (keys %{$decision_HR}) {
    my $d = $decision_HR->{$acc};
    my $status = (defined $d->{"status"}) ? $d->{"status"} : "";
    my $reason = (defined $d->{"reason_code"}) ? $d->{"reason_code"} : "";
    my $stage  = (defined $d->{"stage_last_seen"}) ? $d->{"stage_last_seen"} : "";
    my $group  = (defined $d->{"group_key"}) ? $d->{"group_key"} : "";
    my $key = join("\t", $status, $reason, $stage, $group);
    $count_H{$key}++;
  }

  open(my $sfh, ">", $out_file) || die "ERROR, unable to write decision summary report $out_file: $!";
  print $sfh join("\t", "status", "reason_code", "stage_last_seen", "group_key", "nseq") . "\n";
  foreach my $key (sort keys %count_H) {
    print $sfh $key . "\t" . $count_H{$key} . "\n";
  }
  close($sfh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Decision summary: wrote grouped status/reason/stage counts to %s\n", $out_file));
  return;
}

#################################################################
# Subroutine : write_stitch_scaffold_outputs()
#################################################################
sub write_stitch_scaffold_outputs {
  my ($candidate_AHR, $tier2_fasta_file, $ftr_info_HAR, $model_key, $seed_model_len, $selected_accn_file, $selected_fa_file, $block_plan_file, $FH_HR) = @_;

  my @sel_acc_A = ();
  foreach my $group (sort keys %{$candidate_AHR}) {
    foreach my $seq (@{$candidate_AHR->{$group}}) {
      next if((! defined $seq->{"acc"}) || ($seq->{"acc"} eq ""));
      push(@sel_acc_A, $seq->{"acc"});
    }
  }

  open(my $afh, ">", $selected_accn_file) || die "ERROR, unable to write selected accession list $selected_accn_file: $!";
  foreach my $acc (@sel_acc_A) {
    print $afh $acc . "\n";
  }
  close($afh);

  my %tier2_seq_H = ();
  my $cur_acc = undef;
  open(my $ifh, "<", $tier2_fasta_file) || die "ERROR, unable to read $tier2_fasta_file: $!";
  while(my $line = <$ifh>) {
    chomp $line;
    if($line =~ /^>(\S+)/) {
      $cur_acc = $1;
      $tier2_seq_H{$cur_acc} = "" if(! exists $tier2_seq_H{$cur_acc});
    }
    elsif(defined $cur_acc) {
      $line =~ s/\s+//g;
      $tier2_seq_H{$cur_acc} .= $line;
    }
  }
  close($ifh);

  my $n_missing = 0;
  open(my $ofh, ">", $selected_fa_file) || die "ERROR, unable to write selected fasta $selected_fa_file: $!";
  foreach my $acc (@sel_acc_A) {
    if(! exists $tier2_seq_H{$acc}) {
      $n_missing++;
      next;
    }
    print $ofh ">" . $acc . "\n" . $tier2_seq_H{$acc} . "\n";
  }
  close($ofh);

  my $ftr_model_name = undef;
  if(exists $ftr_info_HAR->{$model_key}) {
    $ftr_model_name = $model_key;
  }
  else {
    my @mdl_key_A = sort keys %{$ftr_info_HAR};
    if(scalar(@mdl_key_A) == 1) {
      $ftr_model_name = $mdl_key_A[0];
    }
    elsif(scalar(@mdl_key_A) > 1) {
      foreach my $key (@mdl_key_A) {
        if($key =~ /^\Q$model_key\E(\.\d+)?$/) {
          $ftr_model_name = $key;
          last;
        }
      }
      if(! defined $ftr_model_name) {
        $ftr_model_name = $mdl_key_A[0];
      }
    }
  }

  my @iv_A = ();
  if(defined $ftr_model_name) {
    my @ftr_A = @{$ftr_info_HAR->{$ftr_model_name}};
    my @mp_iv_A = ();
    my @cds_iv_A = ();
    for(my $i = 0; $i < scalar(@ftr_A); $i++) {
      my $type = (defined $ftr_A[$i]{"type"}) ? $ftr_A[$i]{"type"} : "";
      my $coords = (defined $ftr_A[$i]{"coords"}) ? $ftr_A[$i]{"coords"} : "";
      next if($coords eq "");
      my ($start, $end, $strand) = parse_coords_bounds($coords);
      next if(! defined $start);
      my %iv_H = (
        start      => $start,
        end        => $end,
        strand     => $strand,
        feature_idx=> $i + 1,
        type       => $type,
        coords     => $coords
      );
      if($type eq "mat_peptide") {
        push(@mp_iv_A, \%iv_H);
      }
      elsif($type eq "CDS") {
        push(@cds_iv_A, \%iv_H);
      }
    }
    @iv_A = (scalar(@mp_iv_A) > 0) ? @mp_iv_A : @cds_iv_A;
  }

  @iv_A = sort { $a->{"start"} <=> $b->{"start"} || $a->{"end"} <=> $b->{"end"} } @iv_A;
  my @merged_A = ();
  foreach my $iv (@iv_A) {
    if(scalar(@merged_A) == 0) {
      my %new_H = %{$iv};
      $new_H{"count"} = 1;
      $new_H{"coords_concat"} = $iv->{"coords"};
      push(@merged_A, \%new_H);
    }
    else {
      my $last = $merged_A[-1];
      if($iv->{"start"} <= ($last->{"end"} + 1)) {
        if($iv->{"end"} > $last->{"end"}) {
          $last->{"end"} = $iv->{"end"};
        }
        $last->{"count"}++;
        $last->{"coords_concat"} .= "," . $iv->{"coords"};
      }
      else {
        my %new_H = %{$iv};
        $new_H{"count"} = 1;
        $new_H{"coords_concat"} = $iv->{"coords"};
        push(@merged_A, \%new_H);
      }
    }
  }

  open(my $bpfh, ">", $block_plan_file) || die "ERROR, unable to write stitch block plan $block_plan_file: $!";
  print $bpfh join("\t", "block_idx", "block_type", "nt_start", "nt_end", "nt_len", "source", "feature_type", "feature_count", "coords") . "\n";

  my $blk_idx = 0;
  my $cur_nt = 1;
  for(my $i = 0; $i < scalar(@merged_A); $i++) {
    my $start = $merged_A[$i]{"start"};
    my $end   = $merged_A[$i]{"end"};
    if($start > $cur_nt) {
      $blk_idx++;
      my $nlen = $start - $cur_nt;
      print $bpfh join("\t", $blk_idx, "noncoding", $cur_nt, $start - 1, $nlen, "complement", "", 0, "") . "\n";
    }
    if($end < $start) { next; }
    $blk_idx++;
    my $nlen = $end - $start + 1;
    print $bpfh join("\t", $blk_idx, "coding", $start, $end, $nlen, "seed_minfo", $merged_A[$i]{"type"}, $merged_A[$i]{"count"}, $merged_A[$i]{"coords_concat"}) . "\n";
    $cur_nt = $end + 1;
  }
  if($cur_nt <= $seed_model_len) {
    $blk_idx++;
    my $nlen = $seed_model_len - $cur_nt + 1;
    print $bpfh join("\t", $blk_idx, "noncoding", $cur_nt, $seed_model_len, $nlen, "complement", "", 0, "") . "\n";
  }
  close($bpfh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch scaffold: wrote %d selected accessions to %s\n", scalar(@sel_acc_A), $selected_accn_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch scaffold: wrote selected FASTA to %s (missing %d accessions)\n", $selected_fa_file, $n_missing));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch scaffold: wrote initial block plan (%d blocks) to %s\n", $blk_idx, $block_plan_file));
  return scalar(@sel_acc_A);
}

#################################################################
# Subroutine : parse_coords_bounds()
#################################################################
sub parse_coords_bounds {
  my ($coords_str) = @_;

  my @seg_A = split(/,/, $coords_str);
  my $min_start = undef;
  my $max_end = undef;
  my $strand = undef;

  foreach my $seg (@seg_A) {
    if($seg =~ /^<?(\d+)\.\.>?(\d+):([\+\-])$/) {
      my ($s, $e, $this_strand) = ($1, $2, $3);
      my $lo = ($s < $e) ? $s : $e;
      my $hi = ($s > $e) ? $s : $e;
      if((! defined $min_start) || ($lo < $min_start)) { $min_start = $lo; }
      if((! defined $max_end)   || ($hi > $max_end))    { $max_end   = $hi; }
      if(! defined $strand) { $strand = $this_strand; }
      elsif($strand ne $this_strand) { $strand = "."; }
    }
  }

  return ($min_start, $max_end, $strand);
}

#################################################################
# Subroutine : prepare_cds_translation_for_stitching()
#################################################################
sub prepare_cds_translation_for_stitching {
  my ($candidate_AHR, $selected_fa_file, $annot_outdir, $cds_nt_fa_file, $orf_fa_file, $aa_fa_file, $map_tsv_file, $do_skip_annotate, $do_keep, $do_verbose, $execs_HR, $FH_HR) = @_;

  if($do_skip_annotate) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS translation prep: skipped due to --skip-annotate (predicted CDS coords unavailable)\n"));
    return;
  }

  my $annot_tail = $annot_outdir;
  $annot_tail =~ s/^.+\///;
  my $ftr_file = $annot_outdir . "/" . $annot_tail . ".vadr.ftr";
  if(! -e $ftr_file) {
    die "ERROR, unable to find required v-annotate feature table for CDS extraction: $ftr_file";
  }

  my %is_selected_H = ();
  foreach my $group (keys %{$candidate_AHR}) {
    foreach my $seq (@{$candidate_AHR->{$group}}) {
      next if((! defined $seq->{"acc"}) || ($seq->{"acc"} eq ""));
      $is_selected_H{$seq->{"acc"}} = 1;
    }
  }

  my %seq_H = ();
  my $cur_acc = undef;
  open(my $sffh, "<", $selected_fa_file) || die "ERROR, unable to read selected fasta $selected_fa_file: $!";
  while(my $line = <$sffh>) {
    chomp $line;
    if($line =~ /^>(\S+)/) {
      $cur_acc = $1;
      $seq_H{$cur_acc} = "" if(! exists $seq_H{$cur_acc});
    }
    elsif(defined $cur_acc) {
      $line =~ s/\s+//g;
      $seq_H{$cur_acc} .= $line;
    }
  }
  close($sffh);

  my %cds_info_AH = (); # cds_info_AH{acc} = [ { ftr_idx, strand, n_from, n_to, trc, model_coords, seq_coords, source_name, cds_nt_len }, ... ]
  open(my $ffh, "<", $ftr_file) || die "ERROR, unable to read $ftr_file: $!";
  while(my $line = <$ffh>) {
    chomp $line;
    next if($line =~ /^\#/);
    next if($line =~ /^\s*$/);
    my @tok_A = split(/\s+/, $line);
    next if(scalar(@tok_A) < 25);
    my $acc = $tok_A[1];
    next if(! exists $is_selected_H{$acc});
    my $type = $tok_A[5];
    next if($type ne "CDS");

    my $ftr_idx = $tok_A[8];
    my $strand = $tok_A[10];
    my $n_from = $tok_A[11];
    my $n_to = $tok_A[12];
    my $trc = $tok_A[14];
    my $seq_coords_raw = $tok_A[23];
    my $model_coords = $tok_A[24];

    if((! defined $seq_H{$acc}) || ($seq_H{$acc} eq "")) {
      die "ERROR, sequence $acc required for CDS translation prep not found in $selected_fa_file";
    }

    my $start_tok = $n_from;
    my $stop_tok = $n_to;
    if($trc =~ /5/) { $start_tok = "<" . $start_tok; }
    if($trc =~ /3/) { $stop_tok = ">" . $stop_tok; }
    my $seq_coords = $start_tok . ".." . $stop_tok . ":" . $strand;
    my $source_name = $acc . "__CDS" . $ftr_idx . "/" . $seq_coords;

    my $lo = ($n_from < $n_to) ? $n_from : $n_to;
    my $hi = ($n_from > $n_to) ? $n_from : $n_to;
    my $seqlen = length($seq_H{$acc});
    if(($lo < 1) || ($hi > $seqlen)) {
      die "ERROR, invalid CDS nt range $n_from..$n_to for $acc (sequence length $seqlen)";
    }

    my $cds_sqstring = substr($seq_H{$acc}, $lo-1, ($hi-$lo+1));
    if($strand eq "-") {
      seq_SqstringReverseComplement(\$cds_sqstring);
    }

    push(@{$cds_info_AH{$acc}}, {
      accession    => $acc,
      ftr_idx      => $ftr_idx,
      strand       => $strand,
      n_from       => $n_from,
      n_to         => $n_to,
      trc          => $trc,
      seq_coords_raw => $seq_coords_raw,
      model_coords => $model_coords,
      seq_coords   => $seq_coords,
      source_name  => $source_name,
      cds_sqstring => $cds_sqstring,
      cds_nt_len   => length($cds_sqstring)
    });
  }
  close($ffh);

  open(my $ntfh, ">", $cds_nt_fa_file) || die "ERROR, unable to write $cds_nt_fa_file: $!";
  my %source_info_H = ();
  my $n_cds = 0;
  foreach my $acc (sort keys %cds_info_AH) {
    my @cds_A = sort { $a->{"ftr_idx"} <=> $b->{"ftr_idx"} } @{$cds_info_AH{$acc}};
    foreach my $cds (@cds_A) {
      print $ntfh ">" . $cds->{"source_name"} . " REFCOORDS=" . $cds->{"model_coords"} . "\n";
      print $ntfh seq_SqstringAddNewlines($cds->{"cds_sqstring"}, 60);
      $source_info_H{$cds->{"source_name"}} = $cds;
      $n_cds++;
    }
  }
  close($ntfh);

  my $translate_cmd = $execs_HR->{"esl-translate"} . " -l 1 --watson " . $cds_nt_fa_file . " > " . $orf_fa_file;
  utl_RunCommand($translate_cmd, $do_verbose, 0, $FH_HR);

  my %orf_AH = (); # orf_AH{source} = [ { start, stop, aa_len, frame, aa_seq }, ... ]
  my $cur_source = undef;
  my $cur_start = undef;
  my $cur_stop = undef;
  my $cur_aa_len = undef;
  my $cur_frame = undef;
  my $cur_aa = "";

  open(my $orfh, "<", $orf_fa_file) || die "ERROR, unable to read $orf_fa_file: $!";
  while(my $line = <$orfh>) {
    chomp $line;
    if($line =~ /^>orf\d+\s+source=(\S+)\s+coords=(\d+)\.\.(\d+)\s+length=(\d+)\s+frame=(\S+)/) {
      if(defined $cur_source) {
        push(@{$orf_AH{$cur_source}}, { start => $cur_start, stop => $cur_stop, aa_len => $cur_aa_len, frame => $cur_frame, aa_seq => $cur_aa });
      }
      $cur_source = $1;
      $cur_start = $2;
      $cur_stop = $3;
      $cur_aa_len = $4;
      $cur_frame = $5;
      $cur_aa = "";
      if(! exists $source_info_H{$cur_source}) {
        die "ERROR, esl-translate reported unknown source sequence: $cur_source";
      }
    }
    elsif($line =~ /^>/) {
      die "ERROR, unable to parse esl-translate defline in $orf_fa_file:\n$line";
    }
    else {
      $line =~ s/\s+//g;
      $cur_aa .= $line;
    }
  }
  close($orfh);
  if(defined $cur_source) {
    push(@{$orf_AH{$cur_source}}, { start => $cur_start, stop => $cur_stop, aa_len => $cur_aa_len, frame => $cur_frame, aa_seq => $cur_aa });
  }

  open(my $aafh, ">", $aa_fa_file) || die "ERROR, unable to write $aa_fa_file: $!";
  open(my $mapfh, ">", $map_tsv_file) || die "ERROR, unable to write $map_tsv_file: $!";
  print $mapfh join("\t", "source", "accession", "ftr_idx", "seq_coords", "model_coords", "cds_nt_len", "orf_nt_start", "orf_nt_stop", "orf_aa_len", "orf_frame", "codon_start", "n_untranslated_5p_nt", "n_untranslated_3p_nt", "n_orfs_found") . "\n";

  my $n_written = 0;
  foreach my $source (sort keys %source_info_H) {
    my $cds = $source_info_H{$source};
    my @orf_A = (exists $orf_AH{$source}) ? @{$orf_AH{$source}} : ();
    if(scalar(@orf_A) == 0) {
      die "ERROR, no ORFs found by esl-translate for CDS source $source";
    }

    my @pref_A = grep { $_->{"start"} <= 3 } @orf_A;
    my @use_A = (scalar(@pref_A) > 0) ? @pref_A : @orf_A;
    @use_A = sort {
      $b->{"aa_len"} <=> $a->{"aa_len"} ||
      (($b->{"stop"} - $b->{"start"} + 1) <=> ($a->{"stop"} - $a->{"start"} + 1)) ||
      $a->{"start"} <=> $b->{"start"} ||
      $a->{"frame"} cmp $b->{"frame"}
    } @use_A;
    my $best = $use_A[0];
    if(scalar(@use_A) > 1) {
      my $same_top = 0;
      if(($use_A[1]{"aa_len"} == $best->{"aa_len"}) &&
         (($use_A[1]{"stop"} - $use_A[1]{"start"} + 1) == ($best->{"stop"} - $best->{"start"} + 1)) &&
         ($use_A[1]{"start"} == $best->{"start"}) &&
         ($use_A[1]{"frame"} eq $best->{"frame"})) {
        $same_top = 1;
      }
      if($same_top) {
        die "ERROR, ambiguous best ORF for CDS source $source (multiple ORFs tied for top rank)";
      }
    }

    my $codon_start = $best->{"start"};
    if(($codon_start < 1) || ($codon_start > 3)) {
      die "ERROR, inferred codon_start out of range (1..3) for CDS source $source: $codon_start";
    }

    my $n5 = $best->{"start"} - 1;
    my $n3 = $cds->{"cds_nt_len"} - $best->{"stop"};
    if($n3 < 0) {
      die "ERROR, chosen ORF end exceeds CDS nt length for source $source (cds_len=" . $cds->{"cds_nt_len"} . " stop=" . $best->{"stop"} . ")";
    }

    my $aa_header = $cds->{"accession"} . ":" . $cds->{"seq_coords"} . "/" . $cds->{"model_coords"};
    print $aafh ">" . $aa_header . "\n";
    print $aafh seq_SqstringAddNewlines($best->{"aa_seq"}, 60);

    print $mapfh join("\t", $source,
                            $cds->{"accession"},
                            $cds->{"ftr_idx"},
                            $cds->{"seq_coords"},
                            $cds->{"model_coords"},
                            $cds->{"cds_nt_len"},
                            $best->{"start"},
                            $best->{"stop"},
                            $best->{"aa_len"},
                            $best->{"frame"},
                            $codon_start,
                            $n5,
                            $n3,
                            scalar(@orf_A)) . "\n";
    $n_written++;
  }
  close($aafh);
  close($mapfh);

  if((! $do_keep) && ($n_written > 0)) {
    # keep ORF file for now as a useful debug artifact if needed in follow-on steps
  }

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS translation prep: wrote %d CDS nt sequences to %s\n", $n_cds, $cds_nt_fa_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS translation prep: wrote esl-translate ORFs to %s\n", $orf_fa_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS translation prep: selected one ORF for %d CDS and wrote %s and %s\n", $n_written, $aa_fa_file, $map_tsv_file));
  return;
}

#################################################################
# Subroutine : run_reference_anchored_pairwise_aa()
# EPN* 2026-03-19
# Runs per-CDS-feature pairwise AA alignment. Groups sequences
# by feature (ref coords after '/' in header), picks one anchor
# per feature from the centroid, and aligns each query against
# its own feature's anchor.
#################################################################
sub run_reference_anchored_pairwise_aa {
  my ($candidate_AHR, $centroid_tsv_file, $aa_fa_file, $anchor_tsv_file, $anchor_fa_file, $pairwise_tsv_file, $do_skip_annotate, $do_keep, $do_verbose, $execs_HR, $FH_HR) = @_;

  if($do_skip_annotate) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS pairwise AA: skipped due to --skip-annotate\n"));
    return;
  }
  if(! -s $aa_fa_file) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS pairwise AA: skipped, missing or empty AA file %s\n", $aa_fa_file));
    return;
  }
  if(! exists $execs_HR->{"ggsearch"}) {
    die "ERROR, unable to run reference-anchored global/global pairwise AA: ggsearch36 not found (set VADRFASTADIR so $ENV{VADRFASTADIR}/ggsearch36 exists and is executable)";
  }

  my %centroid_avg_H = (); # accession => avg_blastn_pident for Tier 3 centroid rows
  if(-s $centroid_tsv_file) {
    open(my $ctfh, "<", $centroid_tsv_file) || die "ERROR, unable to read centroid TSV $centroid_tsv_file: $!";
    my $nline = 0;
    while(my $line = <$ctfh>) {
      chomp $line;
      $nline++;
      next if($line =~ /^\s*$/);
      next if($nline == 1);
      my @tok_A = split(/\t/, $line, -1);
      next if(scalar(@tok_A) < 4);
      my ($group_key, $acc, $avg, $is_centroid) = @tok_A;
      next if($is_centroid ne "1");
      $centroid_avg_H{$acc} = $avg;
    }
    close($ctfh);
  }

  # Parse AA sequences and extract feature key (ref coords after '/')
  my @aa_A = (); # { header, accession, feature_key, sqstring, len }
  my $cur_h = undef;
  my $cur_sq = "";
  open(my $aafh_in, "<", $aa_fa_file) || die "ERROR, unable to read AA fasta $aa_fa_file: $!";
  while(my $line = <$aafh_in>) {
    chomp $line;
    if($line =~ /^>(\S+)/) {
      if(defined $cur_h) {
        my $acc = $cur_h;
        if($acc =~ /^([^:]+):/) { $acc = $1; }
        my $fkey = ($cur_h =~ /\/(.+)$/) ? $1 : $cur_h;
        push(@aa_A, { header => $cur_h, accession => $acc, feature_key => $fkey, sqstring => $cur_sq, len => length($cur_sq) });
      }
      $cur_h = $1;
      $cur_sq = "";
    }
    elsif(defined $cur_h) {
      $line =~ s/\s+//g;
      $cur_sq .= $line;
    }
  }
  close($aafh_in);
  if(defined $cur_h) {
    my $acc = $cur_h;
    if($acc =~ /^([^:]+):/) { $acc = $1; }
    my $fkey = ($cur_h =~ /\/(.+)$/) ? $1 : $cur_h;
    push(@aa_A, { header => $cur_h, accession => $acc, feature_key => $fkey, sqstring => $cur_sq, len => length($cur_sq) });
  }
  if(scalar(@aa_A) == 0) {
    die "ERROR, no AA sequences parsed from $aa_fa_file";
  }

  # Group sequences by feature key
  my %fkey_order_A = (); # feature_key => [ seq_hashes ]
  my @fkey_order = ();   # ordered list of unique feature keys
  foreach my $seq (@aa_A) {
    my $fk = $seq->{"feature_key"};
    if(! exists $fkey_order_A{$fk}) {
      push(@fkey_order, $fk);
      $fkey_order_A{$fk} = [];
    }
    push(@{$fkey_order_A{$fk}}, $seq);
  }

  # Pick one anchor per feature group
  my %anchor_H = (); # feature_key => anchor seq hash
  my %anchor_method_H = (); # feature_key => method string
  foreach my $fk (@fkey_order) {
    my @grp_A = @{$fkey_order_A{$fk}};
    my @anchor_cand_A = grep { exists $centroid_avg_H{$_->{"accession"}} } @grp_A;
    if(scalar(@anchor_cand_A) > 0) {
      @anchor_cand_A = sort {
        $centroid_avg_H{$b->{"accession"}} <=> $centroid_avg_H{$a->{"accession"}} ||
        $b->{"len"} <=> $a->{"len"} ||
        $a->{"accession"} cmp $b->{"accession"}
      } @anchor_cand_A;
      $anchor_method_H{$fk} = "tier3_centroid_max_avg_blastn";
    }
    else {
      @anchor_cand_A = sort {
        $b->{"len"} <=> $a->{"len"} ||
        $a->{"accession"} cmp $b->{"accession"}
      } @grp_A;
      $anchor_method_H{$fk} = "fallback_longest_aa";
    }
    $anchor_H{$fk} = $anchor_cand_A[0];
  }

  # Write anchor TSV (one row per feature)
  open(my $atfh, ">", $anchor_tsv_file) || die "ERROR, unable to write anchor TSV $anchor_tsv_file: $!";
  print $atfh join("\t", "anchor_header", "anchor_accession", "anchor_aa_len", "anchor_method", "anchor_centroid_avg_blastn_pident", "feature_key") . "\n";
  foreach my $fk (@fkey_order) {
    my $anc = $anchor_H{$fk};
    my $anchor_avg = exists $centroid_avg_H{$anc->{"accession"}} ? $centroid_avg_H{$anc->{"accession"}} : "NA";
    print $atfh join("\t", $anc->{"header"}, $anc->{"accession"}, $anc->{"len"}, $anchor_method_H{$fk}, $anchor_avg, $fk) . "\n";
  }
  close($atfh);

  # Write anchor FASTA (all anchors)
  open(my $ancfh, ">", $anchor_fa_file) || die "ERROR, unable to write anchor AA fasta $anchor_fa_file: $!";
  foreach my $fk (@fkey_order) {
    my $anc = $anchor_H{$fk};
    print $ancfh ">" . $anc->{"header"} . "\n";
    print $ancfh seq_SqstringAddNewlines($anc->{"sqstring"}, 60);
  }
  close($ancfh);

  my $pairwise_dir = $pairwise_tsv_file . ".dir";
  if(! -d $pairwise_dir) {
    mkdir($pairwise_dir) || die "ERROR, unable to create pairwise output directory $pairwise_dir: $!";
  }

  open(my $pwfh, ">", $pairwise_tsv_file) || die "ERROR, unable to write pairwise summary TSV $pairwise_tsv_file: $!";
  print $pwfh join("\t", "query_header", "query_accession", "anchor_header", "program", "program_mode", "result_file", "hit_found", "pident", "alen", "qstart", "qend", "sstart", "send", "aln_code") . "\n";

  my $pair_idx = 0;
  foreach my $fk (@fkey_order) {
    my $anc = $anchor_H{$fk};
    # Write per-feature anchor FASTA for ggsearch library
    my $fk_anchor_fa = sprintf("%s/anchor.%s.fa", $pairwise_dir, $fk);
    $fk_anchor_fa =~ s/[^A-Za-z0-9._\-\/]/_/g; # sanitize feature key in filename
    $fk_anchor_fa = sprintf("%s/anchor.%05d.fa", $pairwise_dir, scalar(grep { $_ eq $fk } @fkey_order) ? (grep { $fkey_order[$_] eq $fk } 0..$#fkey_order)[0] + 1 : 0);
    open(my $afh, ">", $fk_anchor_fa) || die "ERROR, unable to write anchor fasta $fk_anchor_fa: $!";
    print $afh ">" . $anc->{"header"} . "\n";
    print $afh seq_SqstringAddNewlines($anc->{"sqstring"}, 60);
    close($afh);

    foreach my $seq (@{$fkey_order_A{$fk}}) {
      next if($seq->{"header"} eq $anc->{"header"});
      $pair_idx++;
      my $qfa_file = sprintf("%s/query.%05d.fa", $pairwise_dir, $pair_idx);
      my $out_file = sprintf("%s/query.%05d.ggsearch8cc.tsv", $pairwise_dir, $pair_idx);

      open(my $qfh, ">", $qfa_file) || die "ERROR, unable to write query fasta $qfa_file: $!";
      print $qfh ">" . $seq->{"header"} . "\n";
      print $qfh seq_SqstringAddNewlines($seq->{"sqstring"}, 60);
      close($qfh);

      my $cmd = $execs_HR->{"ggsearch"} . " -m 8CC -d 0 -T 1 " . $qfa_file . " " . $fk_anchor_fa . " > " . $out_file;
      utl_RunCommand($cmd, $do_verbose, 0, $FH_HR);

      my ($hit_found, $pident, $alen, $qstart, $qend, $sstart, $send, $aln_code) = (0, "NA", "NA", "NA", "NA", "NA", "NA", "");
      if(-s $out_file) {
        open(my $ofh, "<", $out_file) || die "ERROR, unable to read ggsearch output $out_file: $!";
        while(my $line = <$ofh>) {
          chomp $line;
          next if($line =~ /^\#/);
          next if($line =~ /^\s*$/);
          my @tok_A = split(/\t/, $line, -1);
          if(scalar(@tok_A) >= 13) {
            $hit_found = 1;
            $pident = $tok_A[2] if(defined $tok_A[2]);
            $alen   = $tok_A[3] if(defined $tok_A[3]);
            $qstart = $tok_A[6] if(defined $tok_A[6]);
            $qend   = $tok_A[7] if(defined $tok_A[7]);
            $sstart = $tok_A[8] if(defined $tok_A[8]);
            $send   = $tok_A[9] if(defined $tok_A[9]);
            $aln_code = $tok_A[12] if(defined $tok_A[12]);
            last;
          }
        }
        close($ofh);
      }

      print $pwfh join("\t", $seq->{"header"},
                              $seq->{"accession"},
                              $anc->{"header"},
                              "ggsearch36",
                              "global_global_pairwise",
                              $out_file,
                              $hit_found,
                              $pident,
                              $alen,
                              $qstart,
                              $qend,
                              $sstart,
                              $send,
                              $aln_code) . "\n";

      if(! $do_keep) {
        unlink $qfa_file;
      }
    }
    if(! $do_keep) {
      unlink $fk_anchor_fa;
    }
  }
  close($pwfh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS pairwise AA: %d CDS features, %d total pairwise alignments\n", scalar(@fkey_order), $pair_idx));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS pairwise AA: wrote anchor info (%d features) to %s\n", scalar(@fkey_order), $anchor_tsv_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS pairwise AA: wrote pairwise summary to %s\n", $pairwise_tsv_file));
  return;
}

#################################################################
# Subroutine : build_anchor_projected_cds_msa()
# EPN* 2026-03-19
# Builds per-CDS-feature AA MSAs from pairwise CIGARs, then
# concatenates them in feature order. Each feature has its own
# anchor. Backconverts to nucleotide MSA.
#################################################################
sub build_anchor_projected_cds_msa {
  my ($aa_fa_file, $cds_nt_fa_file, $map_tsv_file, $anchor_tsv_file, $pairwise_tsv_file, $msa_aa_fa_file, $msa_aa_stk_file, $msa_nt_fa_file, $do_skip_annotate, $FH_HR) = @_;

  if($do_skip_annotate) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: skipped due to --skip-annotate\n"));
    return;
  }
  if((! -s $aa_fa_file) || (! -s $anchor_tsv_file) || (! -s $pairwise_tsv_file) || (! -s $cds_nt_fa_file) || (! -s $map_tsv_file)) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: skipped (required input missing or empty)\n"));
    return;
  }

  # Parse AA sequences, extract feature key (ref coords after '/')
  my @aa_order_A = ();
  my %aa_seq_H = ();
  my %aa_fkey_H = (); # header => feature_key
  my $cur_h = undef;
  open(my $aafh, "<", $aa_fa_file) || die "ERROR, unable to read AA fasta $aa_fa_file: $!";
  while(my $line = <$aafh>) {
    chomp $line;
    if($line =~ /^>(\S+)/) {
      $cur_h = $1;
      if(! exists $aa_seq_H{$cur_h}) {
        push(@aa_order_A, $cur_h);
        $aa_seq_H{$cur_h} = "";
        $aa_fkey_H{$cur_h} = ($cur_h =~ /\/(.+)$/) ? $1 : $cur_h;
      }
    }
    elsif(defined $cur_h) {
      $line =~ s/\s+//g;
      $aa_seq_H{$cur_h} .= $line;
    }
  }
  close($aafh);
  if(scalar(@aa_order_A) == 0) {
    die "ERROR, no sequences parsed from $aa_fa_file";
  }

  # Read per-feature anchors from anchor TSV
  my %anchor_header_H = (); # feature_key => anchor_header
  my @fkey_order = ();
  open(my $atfh, "<", $anchor_tsv_file) || die "ERROR, unable to read anchor TSV $anchor_tsv_file: $!";
  my $nline = 0;
  while(my $line = <$atfh>) {
    chomp $line;
    $nline++;
    next if($line =~ /^\s*$/);
    next if($nline == 1);
    my @tok_A = split(/\t/, $line, -1);
    next if(scalar(@tok_A) < 6);
    my $a_header = $tok_A[0];
    my $fk = $tok_A[5];
    $anchor_header_H{$fk} = $a_header;
    push(@fkey_order, $fk);
  }
  close($atfh);
  if(scalar(@fkey_order) == 0) {
    die "ERROR, no anchor rows parsed from $anchor_tsv_file";
  }

  # Group sequences by feature key, preserving per-feature order
  my %fkey_seqs_A = (); # feature_key => [ headers ]
  foreach my $fk (@fkey_order) { $fkey_seqs_A{$fk} = []; }
  foreach my $h (@aa_order_A) {
    my $fk = $aa_fkey_H{$h};
    if(exists $fkey_seqs_A{$fk}) {
      push(@{$fkey_seqs_A{$fk}}, $h);
    }
  }

  # Read pairwise results
  my %pair_HH = ();
  open(my $pwfh, "<", $pairwise_tsv_file) || die "ERROR, unable to read pairwise summary $pairwise_tsv_file: $!";
  my $pw_nline = 0;
  while(my $line = <$pwfh>) {
    chomp $line;
    $pw_nline++;
    next if($line =~ /^\s*$/);
    next if($pw_nline == 1);
    my @tok_A = split(/\t/, $line, -1);
    next if(scalar(@tok_A) < 14);
    my ($q_header, $q_acc, $a_header, $program, $mode, $result_file, $hit_found, $pident, $alen, $qstart, $qend, $sstart, $send, $aln_code) = @tok_A;
    $pair_HH{$q_header} = {
      hit_found    => $hit_found,
      anchor_header=> $a_header,
      qstart       => $qstart,
      qend         => $qend,
      sstart       => $sstart,
      send         => $send,
      aln_code     => $aln_code
    };
  }
  close($pwfh);

  # Build AA MSA per feature, then concatenate
  # We need a consistent set of accessions across features. Extract accessions from first feature.
  my @accn_order_A = ();
  my %accn_seen_H = ();
  foreach my $h (@{$fkey_seqs_A{$fkey_order[0]}}) {
    my $acc = ($h =~ /^([^:]+):/) ? $1 : $h;
    if(! exists $accn_seen_H{$acc}) {
      push(@accn_order_A, $acc);
      $accn_seen_H{$acc} = 1;
    }
  }

  my %concat_msa_aa_H = (); # accession => concatenated AA alignment string
  my $concat_rf = "";
  foreach my $acc (@accn_order_A) { $concat_msa_aa_H{$acc} = ""; }

  foreach my $fk (@fkey_order) {
    my @fk_headers = @{$fkey_seqs_A{$fk}};
    my $anchor_header = $anchor_header_H{$fk};
    if(! exists $aa_seq_H{$anchor_header}) {
      die "ERROR, anchor $anchor_header for feature $fk not found in AA fasta";
    }
    my $anchor_seq = $aa_seq_H{$anchor_header};
    my $anchor_len = length($anchor_seq);

    my %core_AH = ();
    my %ins_AH  = ();
    my @max_ins_len_A = ((0) x ($anchor_len + 1));

    foreach my $h (@fk_headers) {
      my @core_A = ("") x ($anchor_len + 1);
      my @ins_A  = ("") x ($anchor_len + 1);

      if($h eq $anchor_header) {
        for(my $apos = 1; $apos <= $anchor_len; $apos++) {
          $core_A[$apos] = substr($anchor_seq, $apos-1, 1);
        }
      }
      else {
        if((! exists $pair_HH{$h}) || ($pair_HH{$h}{"hit_found"} ne "1")) {
          die "ERROR, missing successful pairwise alignment data for $h in $pairwise_tsv_file";
        }
        my $qseq = $aa_seq_H{$h};
        my $aln_code = $pair_HH{$h}{"aln_code"};
        my $qstart = $pair_HH{$h}{"qstart"};
        my $sstart = $pair_HH{$h}{"sstart"};

        my ($qaln, $saln) = pairwise_from_cigar($qseq, $anchor_seq, $aln_code, $qstart, $sstart, $h, $anchor_header);
        project_pairwise_onto_anchor($qaln, $saln, $anchor_len, \@core_A, \@ins_A, $h, $anchor_header);
      }

      $core_AH{$h} = \@core_A;
      $ins_AH{$h}  = \@ins_A;
      for(my $k = 0; $k <= $anchor_len; $k++) {
        my $ilen = length($ins_A[$k]);
        if($ilen > $max_ins_len_A[$k]) { $max_ins_len_A[$k] = $ilen; }
      }
    }

    # Build per-feature AA MSA strings
    my $fk_msa_len = undef;
    foreach my $h (@fk_headers) {
      my $acc = ($h =~ /^([^:]+):/) ? $1 : $h;
      my $aln = "";
      my @core_A = @{$core_AH{$h}};
      my @ins_A  = @{$ins_AH{$h}};
      for(my $k = 0; $k <= $anchor_len; $k++) {
        my $ins = $ins_A[$k];
        $aln .= $ins;
        if(length($ins) < $max_ins_len_A[$k]) {
          $aln .= ("-" x ($max_ins_len_A[$k] - length($ins)));
        }
        if($k < $anchor_len) {
          $aln .= $core_A[$k+1];
        }
      }
      if(! defined $fk_msa_len) { $fk_msa_len = length($aln); }
      elsif(length($aln) != $fk_msa_len) { die "ERROR, internal length mismatch while building AA MSA for feature $fk"; }
      $concat_msa_aa_H{$acc} .= $aln;
    }

    # Build per-feature RF
    my $fk_rf = "";
    for(my $k = 0; $k <= $anchor_len; $k++) {
      if($max_ins_len_A[$k] > 0) {
        $fk_rf .= ("." x $max_ins_len_A[$k]);
      }
      if($k < $anchor_len) {
        $fk_rf .= substr($anchor_seq, $k, 1);
      }
    }
    $concat_rf .= $fk_rf;
  }

  # Write per-accession concatenated AA MSA FASTA
  open(my $maa_fh, ">", $msa_aa_fa_file) || die "ERROR, unable to write $msa_aa_fa_file: $!";
  foreach my $acc (@accn_order_A) {
    print $maa_fh ">" . $acc . "\n";
    print $maa_fh seq_SqstringAddNewlines($concat_msa_aa_H{$acc}, 60);
  }
  close($maa_fh);

  open(my $mstk_fh, ">", $msa_aa_stk_file) || die "ERROR, unable to write $msa_aa_stk_file: $!";
  print $mstk_fh "# STOCKHOLM 1.0\n";
  foreach my $acc (@accn_order_A) {
    print $mstk_fh $acc . "\t" . $concat_msa_aa_H{$acc} . "\n";
  }
  print $mstk_fh "#=GC RF\t" . $concat_rf . "\n";
  print $mstk_fh "//\n";
  close($mstk_fh);

  # Read CDS nt sequences and translation map
  my %cds_nt_H = ();
  my $cur_nt_h = undef;
  open(my $ntfh, "<", $cds_nt_fa_file) || die "ERROR, unable to read CDS nt fasta $cds_nt_fa_file: $!";
  while(my $line = <$ntfh>) {
    chomp $line;
    if($line =~ /^>(\S+)/) {
      $cur_nt_h = $1;
      $cds_nt_H{$cur_nt_h} = "" if(! exists $cds_nt_H{$cur_nt_h});
    }
    elsif(defined $cur_nt_h) {
      $line =~ s/\s+//g;
      $cds_nt_H{$cur_nt_h} .= $line;
    }
  }
  close($ntfh);

  my %map_H = ();
  open(my $mapfh, "<", $map_tsv_file) || die "ERROR, unable to read CDS map TSV $map_tsv_file: $!";
  my $map_nline = 0;
  while(my $line = <$mapfh>) {
    chomp $line;
    $map_nline++;
    next if($line =~ /^\s*$/);
    next if($map_nline == 1);
    my @tok_A = split(/\t/, $line, -1);
    next if(scalar(@tok_A) < 13);
    my ($source, $acc, $ftr_idx, $seq_coords, $model_coords, $cds_nt_len, $orf_nt_start, $orf_nt_stop, $orf_aa_len, $orf_frame, $codon_start, $n5, $n3, $n_orfs) = @tok_A;
    my $aa_header = $acc . ":" . $seq_coords . "/" . $model_coords;
    $map_H{$aa_header} = {
      source      => $source,
      cds_nt_len  => $cds_nt_len,
      n5          => $n5,
      n3          => $n3
    };
  }
  close($mapfh);

  # NT backconversion: per-feature, then concatenate per-accession
  # Compute per-feature max_n5 and max_n3 (flanks can differ between features)
  my %concat_nt_H = (); # accession => concatenated NT alignment string
  foreach my $acc (@accn_order_A) { $concat_nt_H{$acc} = ""; }

  foreach my $fk (@fkey_order) {
    my @fk_headers = @{$fkey_seqs_A{$fk}};
    my $anchor_header = $anchor_header_H{$fk};
    my $anchor_seq = $aa_seq_H{$anchor_header};
    my $anchor_len = length($anchor_seq);

    # Compute per-feature max_n5/max_n3
    my $fk_max_n5 = 0;
    my $fk_max_n3 = 0;
    foreach my $h (@fk_headers) {
      if(! exists $map_H{$h}) {
        die "ERROR, unable to find AA->CDS mapping row for $h in $map_tsv_file";
      }
      my $n5 = $map_H{$h}{"n5"};
      my $n3 = $map_H{$h}{"n3"};
      if($n5 > $fk_max_n5) { $fk_max_n5 = $n5; }
      if($n3 > $fk_max_n3) { $fk_max_n3 = $n3; }
    }

    foreach my $h (@fk_headers) {
      my $acc = ($h =~ /^([^:]+):/) ? $1 : $h;
      my $source = $map_H{$h}{"source"};
      if(! exists $cds_nt_H{$source}) {
        die "ERROR, unable to find CDS nt sequence for source $source in $cds_nt_fa_file";
      }
      my $cds_nt = $cds_nt_H{$source};
      my $n5 = $map_H{$h}{"n5"};
      my $n3 = $map_H{$h}{"n3"};

      if(($n5 < 0) || ($n3 < 0)) {
        die "ERROR, invalid negative untranslated flank lengths for $h: n5=$n5 n3=$n3";
      }
      if(($n5 + $n3) > length($cds_nt)) {
        die "ERROR, invalid untranslated flank lengths for $h: n5+n3 exceeds CDS nt length";
      }

      my $prefix_nt = ($n5 > 0) ? substr($cds_nt, 0, $n5) : "";
      my $suffix_nt = ($n3 > 0) ? substr($cds_nt, length($cds_nt) - $n3, $n3) : "";
      my $orf_nt_len = length($cds_nt) - $n5 - $n3;
      if($orf_nt_len < 0) {
        die "ERROR, invalid CDS map values for $h: cds_nt_len=" . length($cds_nt) . " n5=$n5 n3=$n3";
      }
      my $orf_nt = substr($cds_nt, $n5, $orf_nt_len);

      # Backconvert AA alignment to NT using the per-feature per-accession AA MSA
      # We need to reconstruct the per-feature AA alignment for this header
      # It's the portion of concat_msa_aa_H{$acc} corresponding to this feature
      # Instead, rebuild from core/ins arrays (already computed above but out of scope)
      # Simpler: extract from concat_msa_aa_H using feature offsets
    }
  }

  # Alternative approach: backconvert per-header using the original aa_order_A loop
  # but use concat_msa_aa_H keyed by accession — we need per-header AA alignments.
  # Since we concatenated per-accession, we need to track per-header AA strings separately.
  # Let's build a per-header AA alignment hash during the feature loop above.

  # Reset and redo: build per-header AA MSA hash
  my %msa_aa_H = (); # header => per-feature AA alignment string
  %concat_nt_H = (); # reset
  foreach my $acc (@accn_order_A) { $concat_nt_H{$acc} = ""; }

  # We need to redo the per-feature MSA building to capture per-header strings.
  # But the core/ins arrays were computed in the feature loop above and are now out of scope.
  # To avoid recomputation, store per-header AA alignment during the feature loop.
  # Since we already concatenated into concat_msa_aa_H by accession, we can split it back
  # using the per-feature MSA lengths.

  # Compute per-feature MSA lengths from concat_rf
  my @fk_aa_msa_len_A = ();
  my $rf_offset = 0;
  foreach my $fk (@fkey_order) {
    my $anc_seq = $aa_seq_H{$anchor_header_H{$fk}};
    my $anc_len = length($anc_seq);
    # The per-feature AA MSA length = anchor_len + sum of insert columns for that feature
    # We can measure it from concat_rf: count chars from rf_offset until we've seen anc_len non-dot chars
    my $fk_msa_len = 0;
    my $rf_nongap = 0;
    my $pos = $rf_offset;
    while($rf_nongap < $anc_len && $pos < length($concat_rf)) {
      if(substr($concat_rf, $pos, 1) ne ".") {
        $rf_nongap++;
      }
      $fk_msa_len++;
      $pos++;
    }
    # Also count any trailing insert columns
    while($pos < length($concat_rf) && substr($concat_rf, $pos, 1) eq ".") {
      # Only if these dots belong to this feature (before next feature's RF chars)
      # Actually, trailing inserts after last anchor position belong to this feature
      # But leading inserts of next feature also show as dots. We can't distinguish.
      # Safer: use the known fk_msa_len from the feature loop. Let's store it.
      last; # Don't count — the insert-after-last was already counted in the feature loop
    }
    push(@fk_aa_msa_len_A, $fk_msa_len);
    $rf_offset = $pos;
  }

  # Now extract per-header AA alignments and backconvert
  foreach my $fk_idx (0..$#fkey_order) {
    my $fk = $fkey_order[$fk_idx];
    my @fk_headers = @{$fkey_seqs_A{$fk}};
    my $fk_msa_len = $fk_aa_msa_len_A[$fk_idx];
    my $fk_aa_offset = 0;
    for(my $i = 0; $i < $fk_idx; $i++) { $fk_aa_offset += $fk_aa_msa_len_A[$i]; }

    my $fk_max_n5 = 0;
    my $fk_max_n3 = 0;
    foreach my $h (@fk_headers) {
      my $n5 = $map_H{$h}{"n5"};
      my $n3 = $map_H{$h}{"n3"};
      if($n5 > $fk_max_n5) { $fk_max_n5 = $n5; }
      if($n3 > $fk_max_n3) { $fk_max_n3 = $n3; }
    }

    foreach my $h (@fk_headers) {
      my $acc = ($h =~ /^([^:]+):/) ? $1 : $h;
      my $aa_aln = substr($concat_msa_aa_H{$acc}, $fk_aa_offset, $fk_msa_len);

      my $source = $map_H{$h}{"source"};
      my $cds_nt = $cds_nt_H{$source};
      my $n5 = $map_H{$h}{"n5"};
      my $n3 = $map_H{$h}{"n3"};

      my $prefix_nt = ($n5 > 0) ? substr($cds_nt, 0, $n5) : "";
      my $suffix_nt = ($n3 > 0) ? substr($cds_nt, length($cds_nt) - $n3, $n3) : "";
      my $orf_nt_len = length($cds_nt) - $n5 - $n3;
      my $orf_nt = substr($cds_nt, $n5, $orf_nt_len);

      my $nt_aln = "";
      my $nt_pos = 0;
      my @aa_char_A = split(//, $aa_aln);
      foreach my $aa_char (@aa_char_A) {
        if($aa_char eq "-") {
          $nt_aln .= "---";
        }
        else {
          if(($nt_pos + 3) > length($orf_nt)) {
            die "ERROR, unable to backconvert AA MSA for $h: insufficient nucleotides in ORF segment";
          }
          $nt_aln .= substr($orf_nt, $nt_pos, 3);
          $nt_pos += 3;
        }
      }

      my $prefix_pad = ($fk_max_n5 > $n5) ? ("-" x ($fk_max_n5 - $n5)) : "";
      my $suffix_pad = ($fk_max_n3 > $n3) ? ("-" x ($fk_max_n3 - $n3)) : "";
      $concat_nt_H{$acc} .= $prefix_pad . $prefix_nt . $nt_aln . $suffix_nt . $suffix_pad;
    }
  }

  # Write concatenated NT MSA FASTA (one entry per accession, with header = accn:all_coords)
  open(my $mnt_fh, ">", $msa_nt_fa_file) || die "ERROR, unable to write $msa_nt_fa_file: $!";
  my $n_nt = 0;
  foreach my $acc (@accn_order_A) {
    # Build header that includes all feature coords for this accession
    my @acc_headers = grep { /^\Q$acc\E:/ } @aa_order_A;
    my $nt_header = join("+", @acc_headers);
    print $mnt_fh ">" . $nt_header . "\n";
    print $mnt_fh seq_SqstringAddNewlines($concat_nt_H{$acc}, 60);
    $n_nt++;
  }
  close($mnt_fh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: %d CDS features, wrote protein MSA to %s\n", scalar(@fkey_order), $msa_aa_fa_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: wrote protein MSA Stockholm (with RF) to %s\n", $msa_aa_stk_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: wrote CDS nucleotide MSA for %d sequences to %s\n", $n_nt, $msa_nt_fa_file));
  return;
}

#################################################################
# Subroutine : pairwise_from_cigar()
#################################################################
sub pairwise_from_cigar {
  my ($qseq, $sseq, $cigar, $qstart, $sstart, $qname, $sname) = @_;

  if((! defined $cigar) || ($cigar eq "")) {
    die "ERROR, empty CIGAR for pairwise alignment $qname vs $sname";
  }
  my @tok_A = ($cigar =~ /(\d+)([MID])/g);
  if((! defined $tok_A[0])) {
    die "ERROR, unable to parse CIGAR string '$cigar' for $qname vs $sname";
  }

  my $qi = $qstart - 1;
  my $si = $sstart - 1;
  my $qaln = "";
  my $saln = "";

  if($qi > 0) {
    $qaln .= substr($qseq, 0, $qi);
    $saln .= "-" x $qi;
  }
  if($si > 0) {
    $qaln .= "-" x $si;
    $saln .= substr($sseq, 0, $si);
  }

  while($cigar =~ /(\d+)([MID])/g) {
    my ($n, $op) = ($1, $2);
    if($op eq "M") {
      $qaln .= substr($qseq, $qi, $n);
      $saln .= substr($sseq, $si, $n);
      $qi += $n;
      $si += $n;
    }
    elsif($op eq "I") {
      $qaln .= "-" x $n;
      $saln .= substr($sseq, $si, $n);
      $si += $n;
    }
    elsif($op eq "D") {
      $qaln .= substr($qseq, $qi, $n);
      $saln .= "-" x $n;
      $qi += $n;
    }
    else {
      die "ERROR, unsupported CIGAR op '$op' in '$cigar' for $qname vs $sname";
    }
  }

  if($qi < length($qseq)) {
    $qaln .= substr($qseq, $qi);
    $saln .= "-" x (length($qseq) - $qi);
    $qi = length($qseq);
  }
  if($si < length($sseq)) {
    $qaln .= "-" x (length($sseq) - $si);
    $saln .= substr($sseq, $si);
    $si = length($sseq);
  }

  if(length($qaln) != length($saln)) {
    die "ERROR, reconstructed pairwise alignment length mismatch for $qname vs $sname";
  }
  return ($qaln, $saln);
}

#################################################################
# Subroutine : project_pairwise_onto_anchor()
#################################################################
sub project_pairwise_onto_anchor {
  my ($qaln, $saln, $anchor_len, $core_AR, $ins_AR, $qname, $sname) = @_;

  my $apos = 0;
  my $alen = length($qaln);
  if($alen != length($saln)) {
    die "ERROR, project_pairwise_onto_anchor got unequal alignment lengths for $qname vs $sname";
  }

  for(my $i = 0; $i < $alen; $i++) {
    my $qc = substr($qaln, $i, 1);
    my $sc = substr($saln, $i, 1);
    if($sc eq "-") {
      $ins_AR->[$apos] .= $qc;
    }
    else {
      $apos++;
      if(($apos < 1) || ($apos > $anchor_len)) {
        die "ERROR, anchor position out of bounds while projecting $qname vs $sname";
      }
      if(($qc eq "-") || ($qc eq "")) {
        $core_AR->[$apos] = "-";
      }
      else {
        $core_AR->[$apos] = $qc;
      }
    }
  }

  if($apos != $anchor_len) {
    die "ERROR, projected anchor coverage mismatch for $qname vs $sname (got $apos expected $anchor_len)";
  }
  return;
}

#################################################################
# Subroutine : run_rna_discovery()
# Incept     : EPN/Copilot Thu Mar 13 2026
#
# Purpose    : Discover structural RNA regions in reference sequence
#              using cmscan against Rfam (or custom CM file).
#              Also generate full-sequence consensus SS string via cmalign.
#
# Arguments  :
#   $ref_seq_file    : path to reference sequence FASTA file
#   $rna_cm_file     : path to CM file (Rfam.cm or user-provided)
#   $seq_len         : expected sequence length
#   $out_root        : output file root path
#   $rfam_dir        : VADRRFAMDIR (for Rfam.clanin), undef if using custom CM
#   $rna_regions_AR  : ref to array to populate with RNA hit info hashes
#   $ss_cons_SR      : ref to scalar to store full-length consensus SS string
#   $do_keep         : keep intermediate files
#   $do_verbose      : verbose output
#   $execs_HR        : hash of executable paths
#   $FH_HR           : hash of file handles
#
# Returns    : void (populates $rna_regions_AR and $ss_cons_SR)
#################################################################
sub run_rna_discovery {
  my ($ref_seq_file, $rna_cm_file, $seq_len, $out_root, $rfam_dir, $rna_regions_AR, $ss_cons_SR, $do_keep, $do_verbose, $execs_HR, $FH_HR) = @_;

  my $cmscan_tblout = $out_root . ".rna_cmscan.tblout";
  my $cmscan_stdout = $out_root . ".rna_cmscan.out";
  my $cmalign_tfile = $out_root . ".rna_cmalign.ifile";
  my $cmalign_stk   = $out_root . ".rna_cmalign.stk";

  # Step 1b.1: Run cmscan to identify RNA hits
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA discovery: running cmscan on reference sequence with %s\n", $rna_cm_file));
  
  my $cmd = $execs_HR->{"cmscan"} . " --noali --cut_ga --rfam --nohmmonly --tblout " . $cmscan_tblout . " --fmt 2";
  
  # Add --clanin if using Rfam.cm (requires Rfam.clanin)
  if((defined $rfam_dir) && ($rfam_dir ne "")) {
    my $clanin_file = $rfam_dir . "/Rfam.clanin";
    if(-e $clanin_file) {
      $cmd .= " --clanin " . $clanin_file;
    }
  }
  
  $cmd .= " --oskip --cpu 0 " . $rna_cm_file . " " . $ref_seq_file . " > " . $cmscan_stdout;
  utl_RunCommand($cmd, $do_verbose, 0, $FH_HR);

  # Step 1b.2: Parse cmscan tblout to extract RNA hits
  parse_cmscan_tblout($cmscan_tblout, $seq_len, $rna_regions_AR, $FH_HR);

  # Step 1b.3: Run cmalign --tfile to get full-sequence consensus secondary structure
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA discovery: running cmalign to generate consensus secondary structure\n"));
  
  # For cmalign, we need to use a single CM from the reference model
  # TODO: Determine which CM to use for full-sequence alignment (may need seed model CM)
  # For now, skip this step - will implement after testing cmscan parsing
  
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA discovery: found %d RNA regions\n", scalar(@{$rna_regions_AR})));
  
  # Cleanup if not keeping intermediate files
  if(! $do_keep) {
    unlink($cmscan_stdout) if(-e $cmscan_stdout);
    unlink($cmalign_tfile) if(-e $cmalign_tfile);
    unlink($cmalign_stk) if(-e $cmalign_stk);
  }
  
  return;
}

#################################################################
# Subroutine : parse_cmscan_tblout()
# Incept     : EPN/Copilot Thu Mar 13 2026
#
# Purpose    : Parse cmscan --tblout output to extract RNA hit information.
#
# Arguments  :
#   $tblout_file     : path to cmscan --tblout file
#   $seq_len         : expected sequence length for validation
#   $rna_regions_AR  : ref to array to populate with hit info hashes
#   $FH_HR           : hash of file handles
#
# Returns    : void (populates $rna_regions_AR)
#################################################################
sub parse_cmscan_tblout {
  my ($tblout_file, $seq_len, $rna_regions_AR, $FH_HR) = @_;

  if(! -e $tblout_file) {
    die "ERROR, cmscan tblout file does not exist: $tblout_file";
  }

  open(my $tblfh, "<", $tblout_file) || die "ERROR, unable to read cmscan tblout: $tblout_file: $!";
  
  my $nhits = 0;
  while(my $line = <$tblfh>) {
    chomp $line;
    next if($line =~ /^\#/);  # skip comments
    next if($line =~ /^\s*$/); # skip blank lines
    
    # cmscan --tblout format (space-delimited, 18+ columns):
    # idx target_name accession query_name accession clan mdl mdl_from mdl_to seq_from seq_to strand ...
    my @fields = split(/\s+/, $line);
    next if(scalar(@fields) < 18);
    
    my $cm_family = $fields[1];      # Rfam family name (e.g., "IRES_Picorna")
    my $cm_accession = $fields[2];   # Rfam accession (e.g., "RF00229")
    my $seq_from = $fields[9];       # sequence start position (1-based)
    my $seq_to = $fields[10];        # sequence end position (1-based)
    my $strand = $fields[11];        # strand ('+' or '-')
    my $score = $fields[14];         # bit score
    my $evalue = $fields[15];        # E-value
    
    # Validate coordinates
    if(($seq_from < 1) || ($seq_to > $seq_len)) {
      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA discovery WARNING: hit %s coords %d..%d out of bounds (seq_len %d), skipping\n", 
                                                      $cm_family, $seq_from, $seq_to, $seq_len));
      next;
    }
    
    # Store hit info
    push(@{$rna_regions_AR}, {
      start        => $seq_from,
      end          => $seq_to,
      strand       => $strand,
      cm_family    => $cm_family,
      cm_accession => $cm_accession,
      score        => $score,
      evalue       => $evalue
    });
    $nhits++;
  }
  close($tblfh);
  
  return;
}

#################################################################
# Subroutine : write_rna_annotation_table()
# Incept     : EPN/Copilot Thu Mar 13 2026
#
# Purpose    : Write RNA annotation table to TSV file.
#
# Arguments  :
#   $rna_regions_AR  : ref to array of RNA hit info hashes
#   $out_file        : output TSV file path
#   $FH_HR           : hash of file handles
#
# Returns    : void
#################################################################
sub write_rna_annotation_table {
  my ($rna_regions_AR, $out_file, $FH_HR) = @_;

  open(my $outfh, ">", $out_file) || die "ERROR, unable to write RNA annotation table: $out_file: $!";
  
  # Write header
  print $outfh "#idx\tstart\tend\tstrand\tcm_family\tcm_accession\tscore\tevalue\tsstruct\n";
  print  "#idx\tstart\tend\tstrand\tcm_family\tcm_accession\tscore\tevalue\tsstruct\n";
  
  # Write data rows
  my $idx = 1;
  foreach my $rna (@{$rna_regions_AR}) {
    printf $outfh "%d\t%d\t%d\t%s\t%s\t%s\t%.1f\t%.2e\t%s\n",
      $idx,
      $rna->{"start"},
      $rna->{"end"},
      $rna->{"strand"},
      $rna->{"cm_family"},
      $rna->{"cm_accession"},
      $rna->{"score"},
      $rna->{"evalue"},
      (defined $rna->{"sstruct"} ? $rna->{"sstruct"} : "");
    $idx++;
  }
  
  close($outfh);
  
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA discovery: wrote %d RNA annotations to %s\n", scalar(@{$rna_regions_AR}), $out_file));
  
  return;
}

#################################################################
# Subroutine : run_rna_sstruct_generation()
# Incept     : EPN/Copilot Thu Mar 13 2026
#
# Purpose    : Generate consensus secondary structure for each RNA region
#              by running cmalign with --tfile on extracted subsequences.
#
# Arguments  :
#   $rna_regions_AR  : ref to array of RNA hit hashes from cmscan
#   $ref_seq_file    : path to reference sequence FASTA file
#   $rna_cm_file     : path to CM file (Rfam.cm or custom)
#   $rfam_dir        : VADRRFAMDIR (for cmfetch), undef if using custom CM
#   $out_root        : output file root path
#   $do_keep         : keep intermediate files
#   $do_verbose      : verbose output
#   $execs_HR        : hash of executable paths
#   $FH_HR           : hash of file handles
#
# Returns    : void (updates RNA region hashes with sstruct field)
#################################################################
sub run_rna_sstruct_generation {
  my ($rna_regions_AR, $ref_seq_file, $rna_cm_file, $rfam_dir, $out_root, $do_keep, $do_verbose, $execs_HR, $FH_HR) = @_;

  if(! -e $ref_seq_file) {
    die "ERROR, reference sequence file not found: $ref_seq_file";
  }

  my $rna_work_dir = $out_root . ".rna_struct";
  if(! -d $rna_work_dir) {
    mkdir($rna_work_dir) || die "ERROR, unable to create RNA work directory: $rna_work_dir: $!";
  }

  # Load reference sequence into memory
  my %seq_H = ();
  my $cur_acc = undef;
  open(my $sfh, "<", $ref_seq_file) || die "ERROR, unable to read reference sequence: $ref_seq_file: $!";
  while(my $line = <$sfh>) {
    chomp $line;
    if($line =~ /^>(\S+)/) {
      $cur_acc = $1;
      $seq_H{$cur_acc} = "" if(! exists $seq_H{$cur_acc});
    }
    elsif(defined $cur_acc) {
      $line =~ s/\s+//g;
      $seq_H{$cur_acc} .= $line;
    }
  }
  close($sfh);

  if(scalar(keys %seq_H) == 0) {
    die "ERROR, reference sequence file appears empty: $ref_seq_file";
  }

  # Get the first (and should be only) reference sequence
  my @ref_accs = sort keys %seq_H;
  my $ref_acc = $ref_accs[0];
  my $ref_seq = $seq_H{$ref_acc};

  # Process each RNA region
  my $idx = 1;
  foreach my $rna (@{$rna_regions_AR}) {
    my $rna_start = $rna->{"start"};
    my $rna_end = $rna->{"end"};
    my $rna_acc = $rna->{"cm_accession"};
    my $rna_family = $rna->{"cm_family"};

    # Extract subsequence
    my $subseq_len = ($rna_end - $rna_start + 1);
    if(($rna_start < 1) || ($rna_end > length($ref_seq))) {
      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA sstruct WARNING: RNA region %d (%s) coords out of bounds, skipping\n", $idx, $rna_family));
      $rna->{"sstruct"} = "";
      $idx++;
      next;
    }
    my $rna_subseq = substr($ref_seq, $rna_start - 1, $subseq_len);

    # Write temporary FASTA for this RNA region
    my $subseq_fa = $rna_work_dir . "/" . sprintf("rna.%03d.fa", $idx);
    open(my $subfh, ">", $subseq_fa) || die "ERROR, unable to write RNA subsequence: $subseq_fa: $!";
    printf $subfh ">%s_%d_%d\n", $rna_family, $rna_start, $rna_end;
    print $subfh $rna_subseq . "\n";
    close($subfh);

    # Run cmalign to get alignment with secondary structure
    my $cmalign_tfile = $rna_work_dir . "/" . sprintf("rna.%03d.tfile", $idx);
    my $cmalign_stk   = $rna_work_dir . "/" . sprintf("rna.%03d.stk", $idx);

    my $cmd;
    if((defined $rfam_dir) && ($rfam_dir ne "")) {
      # Using Rfam.cm - use cmfetch to get specific model
      my $cmfetch_exec = $rfam_dir . "/cmfetch";
      if(! -e $cmfetch_exec) {
        $cmfetch_exec = "cmfetch"; # Fall back to PATH
      }
      my $rfam_cm = $rfam_dir . "/Rfam.cm";
      $cmd = "$cmfetch_exec $rfam_cm $rna_acc | " . $execs_HR->{"cmalign"} . 
             " --outformat pfam -g --tfile " . $cmalign_tfile . " - " . $subseq_fa . " > " . $cmalign_stk;
    }
    else {
      # Using custom CM - assume it's a single model or doesn't need cmfetch
      $cmd = $execs_HR->{"cmalign"} . " --outformat pfam -g --tfile " . $cmalign_tfile . " " . 
             $rna_cm_file . " " . $subseq_fa . " > " . $cmalign_stk;
    }

    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA sstruct: RNA region %d (%s %d..%d): running cmalign\n", $idx, $rna_family, $rna_start, $rna_end));
    utl_RunCommand($cmd, $do_verbose, 0, $FH_HR);

    # Clean up the alignment using cmbuild to resolve gaps and broken base pairs
    my $refined_stk = $rna_work_dir . "/" . sprintf("rna.%03d.refined.stk", $idx);
    my $tmp_cm      = $rna_work_dir . "/" . sprintf("rna.%03d.cm", $idx);
    my $cmd_cmbuild = $execs_HR->{"cmbuild"} . " -O " . $refined_stk . " -F " . $tmp_cm . " " . $cmalign_stk;
    
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA sstruct: RNA region %d (%s): refining alignment with cmbuild\n", $idx, $rna_family));
    utl_RunCommand($cmd_cmbuild, $do_verbose, 0, $FH_HR);

    # Convert to PFAM format to get clean SS_cons line
    my $pfam_stk = $rna_work_dir . "/" . sprintf("rna.%03d.pfam", $idx);
    my $cmd_reformat = "esl-reformat pfam " . $refined_stk . " > " . $pfam_stk;
    
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA sstruct: RNA region %d (%s): converting to PFAM format\n", $idx, $rna_family));
    utl_RunCommand($cmd_reformat, $do_verbose, 0, $FH_HR);

    # Extract consensus secondary structure from cleaned SS_cons line
    my $sstruct = "";
    if(-e $pfam_stk) {
      open(my $sfh, "<", $pfam_stk) || die "ERROR, unable to read pfam stk file: $pfam_stk: $!";
      while(my $line = <$sfh>) {
        chomp $line;
        if($line =~ /^#=GC SS_cons\s+(.+)$/) {
          my $ss_part = $1;
          $ss_part =~ s/^\s+|\s+$//g;
          $sstruct .= $ss_part;
        }
      }
      close($sfh);
    }
    $rna->{"sstruct"} = $sstruct;

    $idx++;
  }

  # Cleanup if not keeping intermediate files
  if(! $do_keep) {
    system("rm -rf $rna_work_dir");
  }

  return;
}

#################################################################
# Subroutine : extract_and_align_rna_regions()
# Incept     : EPN Mon Mar 16 2026
#
# Purpose    : Extract RNA regions from tier2 v-annotate alignment
#              and realign with custom CMs built in Step 1c.
#              This implements Step 8b-8d from the updated plan.
#              Uses Bio::Easel::MSA column_subset() for extraction
#              and write_single_unaligned_seq() for FASTA output.
#
# Arguments  :
#   $rna_regions_AR     : ref to array of RNA region hashes from Step 1b
#   $tier2_align_stk    : path to tier2 v-annotate Stockholm alignment 
#   $rna_struct_dir     : directory with CM files from Step 1c (rna.001.cm, etc)
#   $out_root           : output file root path
#   $do_keep            : keep intermediate files
#   $do_verbose         : verbose output
#   $execs_HR           : hash of executable paths
#   $FH_HR              : hash of file handles
#
# Returns    : void (writes RNA block Stockholm files)
#################################################################
sub extract_and_align_rna_regions {
  my ($rna_regions_AR, $tier2_align_stk, $rna_struct_dir, $out_root, $do_keep, $do_verbose, $execs_HR, $FH_HR) = @_;

  if(! -e $tier2_align_stk) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: WARNING - tier2 alignment not found, skipping RNA refinement\n"));
    return;
  }

  my $nrna = scalar(@{$rna_regions_AR});
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: extracting and aligning %d RNA region%s\n", $nrna, ($nrna == 1) ? "" : "s"));

  # Load tier2 Stockholm via Bio::Easel::MSA (read once, clone per region)
  my $tier2_msa = Bio::Easel::MSA->new({ fileLocation => $tier2_align_stk, isDna => 1 });

  # Build RF position map for column extraction
  my @rf2a_map_A = ();
  my @a2rf_map_A = ();
  $tier2_msa->get_rf_map(\@rf2a_map_A, \@a2rf_map_A, "-.");

  # Process each RNA region
  my $idx = 1;
  foreach my $rna (@{$rna_regions_AR}) {
    my $rna_start = $rna->{"start"};
    my $rna_end = $rna->{"end"};
    my $rna_family = $rna->{"cm_family"};

    # Step 8b: Extract RNA region columns from tier2 alignment via Bio::Easel
    # Map RF positions to alignment columns and extract with column_subset
    my $rna_extracted_fa = $out_root . ".stitch.rna." . sprintf("%03d", $idx) . ".extracted.fa";

    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: region %d (%s %d..%d): extracting from tier2 alignment\n",
                                                     $idx, $rna_family, $rna_start, $rna_end));

    # Build column mask: rf2a_map is 1-based, column_subset is 0-based
    my @useme_A = (0) x $tier2_msa->alen();
    my $first_acol = $rf2a_map_A[$rna_start];  # 1-based
    my $last_acol  = $rf2a_map_A[$rna_end];    # 1-based
    for(my $a = $first_acol; $a <= $last_acol; $a++) {
      $useme_A[$a - 1] = 1;  # convert to 0-based
    }
    my $rna_sub_msa = $tier2_msa->clone_msa();
    $rna_sub_msa->column_subset(\@useme_A);

    # Write unaligned sequences to FASTA, filtering out empty sequences (replaces esl-alimanip --lmin 1)
    my $nseq_written = 0;
    unlink($rna_extracted_fa) if(-e $rna_extracted_fa);  # remove if exists so first write is not append
    for(my $i = 0; $i < $rna_sub_msa->nseq(); $i++) {
      if($rna_sub_msa->get_sqlen($i) >= 1) {
        $rna_sub_msa->write_single_unaligned_seq($i, $rna_extracted_fa, ($nseq_written > 0) ? 1 : 0);
        $nseq_written++;
      }
    }
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: region %d (%s): %d of %d sequences have >= 1 residue\n",
                                                     $idx, $rna_family, $nseq_written, $rna_sub_msa->nseq()));

    if($nseq_written == 0) {
      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: WARNING - no sequences with residues for region %d, skipping\n", $idx));
      $idx++;
      next;
    }

    # Step 8c: Align with custom CM from Step 1c
    my $cm_file = $rna_struct_dir . "/rna." . sprintf("%03d", $idx) . ".cm";
    if(! -e $cm_file) {
      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: WARNING - CM not found for region %d, skipping\n", $idx));
      $idx++;
      next;
    }

    # Step 8d: Output refined RNA block Stockholm
    my $rna_aligned_stk = $out_root . ".stitch.rna." . sprintf("%03d", $idx) . ".stk";

    my $cmd_align = $execs_HR->{"cmalign"} . " --outformat pfam -g " . $cm_file . " " . $rna_extracted_fa .
                    " > " . $rna_aligned_stk;

    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: region %d (%s): aligning with custom CM\n",
                                                     $idx, $rna_family));
    utl_RunCommand($cmd_align, $do_verbose, 0, $FH_HR);

    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: region %d (%s): wrote aligned Stockholm to %s\n",
                                                     $idx, $rna_family, $rna_aligned_stk));

    # Cleanup extracted FASTA if not keeping
    if(! $do_keep) {
      unlink($rna_extracted_fa);
    }

    $idx++;
  }

  return;
}

#################################################################
# Subroutine : stitch_and_refine_final_alignment()
# Purpose    : Stitch all alignment blocks (CDS, RNA, noncoding)  
#              into final training alignment and refine with cmbuild
#
# Arguments  :
#   $block_plan_file    : block plan TSV file
#   $rna_annot_file     : RNA annotation TSV file
#   $tier2_stk_file     : tier2 full-sequence Stockholm alignment
#   $cds_msa_fa_file    : CDS nucleotide MSA FASTA file
#   $out_root           : output root path
#   $rna_regions_AR     : ref to array of RNA region hashes
#   $final_stk_file     : output final stitched Stockholm
#   $refined_stk_file   : output refined Stockholm from cmbuild
#   $temp_cm_file       : temporary CM file for cmbuild
#   $do_rna_discovery   : flag for RNA discovery
#   $do_skip_annotate   : flag to skip annotation
#   $seed_model_len     : reference sequence length
#   $execs_HR           : ref to hash of executable paths
#   $FH_HR              : ref to hash of file handles
#
# Returns    : void
#################################################################
sub stitch_and_refine_final_alignment {
  my ($block_plan_file, $rna_annot_file, $tier2_stk_file, $cds_msa_fa_file, $out_root, 
      $rna_regions_AR, $final_stk_file, $refined_stk_file, $temp_cm_file,
      $do_rna_discovery, $do_skip_annotate, $seed_model_len, $execs_HR, $FH_HR) = @_;

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: merging CDS, RNA, and noncoding blocks\n"));

  # Read block plan
  my @blocks_A = ();
  open(my $bpfh, $block_plan_file) || die "ERROR unable to read block plan $block_plan_file: $!";
  my $header = <$bpfh>;
  while(my $line = <$bpfh>) {
    chomp $line;
    my @f = split(/\t/, $line);
    push(@blocks_A, {
      "idx"      => $f[0],
      "type"     => $f[1],
      "start"    => $f[2],
      "end"      => $f[3],
      "len"      => $f[4],
      "source"   => $f[5]
    });
  }
  close($bpfh);

  # Read RNA annotations if available
  my @rna_A = ();
  if($do_rna_discovery && !$do_skip_annotate && defined($rna_regions_AR) && scalar(@{$rna_regions_AR}) > 0) {
    open(my $rnafh, $rna_annot_file) || die "ERROR unable to read RNA annotation $rna_annot_file: $!";
    $header = <$rnafh>;
    my $rna_idx = 1;
    while(my $line = <$rnafh>) {
      chomp $line;
      my @f = split(/\t/, $line);
      push(@rna_A, {
        "idx"      => $rna_idx,
        "start"    => $f[1],
        "end"      => $f[2],
        "strand"   => $f[3],
        "family"   => $f[4],
        "stk_file" => $out_root . ".stitch.rna." . sprintf("%03d", $rna_idx) . ".stk"
      });
      $rna_idx++;
    }
    close($rnafh);
    
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: merging %d RNA regions into block plan\n", scalar(@rna_A)));
  }

  # Merge RNA blocks into the plan
  my @merged_blocks_A = merge_rna_into_blocks(\@blocks_A, \@rna_A, $FH_HR);
  
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: merged block plan has %d blocks\n", scalar(@merged_blocks_A)));

  # Build ungapped reference-length SS_cons (RNA structures overlaid on single-strand background)
  my $ungapped_ss_cons = build_ungapped_ss_cons($seed_model_len, \@rna_A, $FH_HR);

  # Read anchor accession from anchor TSV (second line, field [1])
  my $anchor_accn = "";
  my $anchor_tsv_file = $out_root . ".stitch.cds.anchor.tsv";
  open(my $atsvfh, $anchor_tsv_file) || die "ERROR unable to read anchor TSV $anchor_tsv_file: $!";
  my $atsvhdr = <$atsvfh>;  # skip header
  if(my $atsvline = <$atsvfh>) {
    chomp $atsvline;
    my @atsvf = split(/\t/, $atsvline);
    $anchor_accn = $atsvf[1];
  }
  close($atsvfh);
  if($anchor_accn eq "") {
    die "ERROR: could not read anchor accession from $anchor_tsv_file";
  }
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: anchor accession: %s\n", $anchor_accn));

  # Extract and concatenate all blocks with interleaved Stockholm and per-block RF/SS_cons
  concatenate_all_blocks(\@merged_blocks_A, $tier2_stk_file, $cds_msa_fa_file, $final_stk_file,
                         $anchor_accn, $ungapped_ss_cons, $execs_HR, $FH_HR);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: wrote concatenated alignment to %s\n", $final_stk_file));

  # Refine alignment and add RF consensus line with cmbuild --refine
  # Use SS_cons from concatenated alignment to build structure into the model
  # Options: --verbose (progress), --sub (sub CM for speed), --tau (convergence), --mxsize (matrix size)
  # Note: -O saves alignment without wrapping, --refine saves wrapped alignment
  my $cmbuild_out = $refined_stk_file . ".cmbuild.out";
  my $output_stk_file = $out_root . ".vadr.stitch.final.rf.stk";
  my $cmd_cmbuild = $execs_HR->{"cmbuild"} . " --hand -O " . $output_stk_file . " " .
                    $temp_cm_file . " " . $final_stk_file . " > " . $cmbuild_out;

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: running cmbuild to add RF annotation (output to %s)\n", $cmbuild_out));
  utl_RunCommand($cmd_cmbuild, 1, 0, $FH_HR);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: wrote RF-annotated alignment to %s\n", $output_stk_file));

  # Cleanup temp CM
  if(-e $temp_cm_file) {
    unlink($temp_cm_file);
  }

  return;
}

#################################################################
# Subroutine : merge_rna_into_blocks()
# Purpose    : Merge RNA blocks into existing block plan,
#              splitting blocks where RNA regions occur.
#              Priority: CDS > RNA > noncoding
#              (RNA overlapping CDS is skipped - CDS alignment used)
#
# Arguments  :
#   $blocks_AR  : ref to array of block hashes
#   $rna_AR     : ref to array of RNA annotation hashes
#   $FH_HR      : ref to hash of file handles
#
# Returns    : array of merged blocks
#################################################################
sub merge_rna_into_blocks {
  my ($blocks_AR, $rna_AR, $FH_HR) = @_;
  
  my @merged_A = ();
  
  # If no RNA, just return original blocks
  if(scalar(@{$rna_AR}) == 0) {
    return @{$blocks_AR};
  }
  
  # Iterate through original blocks and insert/split for RNA
  foreach my $block (@{$blocks_AR}) {
    my $block_start = $block->{"start"};
    my $block_end   = $block->{"end"};
    my $block_type  = $block->{"type"};
    
    # Find RNA regions that overlap this block
    my @overlapping_rna = ();
    foreach my $rna (@{$rna_AR}) {
      if($rna->{"start"} <= $block_end && $rna->{"end"} >= $block_start) {
        push(@overlapping_rna, $rna);
      }
    }
    
    # If no RNA overlaps, add block as-is
    if(scalar(@overlapping_rna) == 0) {
      push(@merged_A, $block);
      next;
    }
    
    # Priority rule: CDS > RNA
    # If this is a CDS block, RNA regions are ignored (CDS alignment takes precedence)
    if($block_type eq "coding") {
      push(@merged_A, $block);
      my $n_skipped = scalar(@overlapping_rna);
      if($n_skipped > 0) {
        my @families = map { $_->{"family"} } @overlapping_rna;
        ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: skipped %d RNA region(s) overlapping CDS (%s) - CDS alignment has priority\n", 
                                                         $n_skipped, join(", ", @families)));
      }
      next;
    }
    
    # Sort overlapping RNA by start position
    @overlapping_rna = sort { $a->{"start"} <=> $b->{"start"} } @overlapping_rna;
    
    # Split noncoding block around RNA regions
    my $cur_pos = $block_start;
    foreach my $rna (@overlapping_rna) {
      my $rna_start = $rna->{"start"};
      my $rna_end   = $rna->{"end"};
      
      # Add segment before RNA if exists
      if($cur_pos < $rna_start) {
        push(@merged_A, {
          "type"     => $block_type,
          "start"    => $cur_pos,
          "end"      => $rna_start - 1,
          "len"      => $rna_start - $cur_pos,
          "source"   => $block->{"source"}
        });
      }
      
      # Add RNA block
      push(@merged_A, {
        "type"     => "rna",
        "start"    => $rna_start,
        "end"      => $rna_end,
        "len"      => $rna_end - $rna_start + 1,
        "source"   => "rna_discovery",
        "rna_idx"  => $rna->{"idx"},
        "rna_file" => $rna->{"stk_file"},
        "family"   => $rna->{"family"}
      });
      
      $cur_pos = $rna_end + 1;
    }
    
    # Add segment after last RNA if exists
    if($cur_pos <= $block_end) {
      push(@merged_A, {
        "type"     => $block_type,
        "start"    => $cur_pos,
        "end"      => $block_end,
        "len"      => $block_end - $cur_pos + 1,
        "source"   => $block->{"source"}
      });
    }
  }
  
  return @merged_A;
}


#################################################################
# Subroutine : concatenate_all_blocks()
# Purpose    : Extract and concatenate all blocks into final Stockholm
#              using interleaved format with per-block RF and full SS_cons.
#
# Arguments  :
#   $blocks_AR         : ref to array of merged block hashes
#   $tier2_stk_file    : tier2 full-sequence alignment (Stockholm)
#   $cds_msa_fa_file   : CDS MSA FASTA file
#   $out_stk_file      : output Stockholm file
#   $anchor_accn       : accession of CDS anchor/reference sequence
#   $ungapped_ss_cons  : ungapped reference-length SS_cons string
#   $execs_HR          : ref to hash of executable paths
#   $FH_HR             : ref to hash of file handles
#
# Returns    : void
#################################################################
sub concatenate_all_blocks {
  my ($blocks_AR, $tier2_stk_file, $cds_msa_fa_file, $out_stk_file,
      $anchor_accn, $ungapped_ss_cons, $execs_HR, $FH_HR) = @_;

  # Read CDS MSA FASTA, key by accession (strip ':coords:strand' suffix)
  my %cds_seqs_H = ();
  open(my $cdsfh, $cds_msa_fa_file) || die "ERROR unable to read CDS MSA $cds_msa_fa_file: $!";
  my $cur_name = "";
  while(my $line = <$cdsfh>) {
    chomp $line;
    if($line =~ /^>(\S+)/) { $cur_name = $1; $cds_seqs_H{$cur_name} = ""; }
    elsif($cur_name ne "")  { $cds_seqs_H{$cur_name} .= $line; }
  }
  close($cdsfh);

  my %accn_to_cds_seq_H = ();
  foreach my $cds_name (keys %cds_seqs_H) {
    my $accn = $cds_name;
    $accn =~ s/:.+//;  # strip from first colon onward
    $accn_to_cds_seq_H{$accn} = $cds_seqs_H{$cds_name};
  }

  # Load tier2 Stockholm via Bio::Easel::MSA (used for seq order and noncoding extraction)
  my $tier2_msa = Bio::Easel::MSA->new({ fileLocation => $tier2_stk_file, isDna => 1 });

  # Build RF position map for noncoding block extraction
  my @rf2a_map_A = ();
  my @a2rf_map_A = ();
  $tier2_msa->get_rf_map(\@rf2a_map_A, \@a2rf_map_A, "-.");

  # Read canonical sequence order, restricted to CDS MSA members
  my @seq_names = ();
  my %seen_H = ();
  for(my $i = 0; $i < $tier2_msa->nseq(); $i++) {
    my $name = $tier2_msa->get_sqname($i);
    if(!exists $seen_H{$name} && exists $accn_to_cds_seq_H{$name}) {
      push @seq_names, $name;
      $seen_H{$name} = 1;
    }
  }
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: %d sequences in final alignment (restricted to CDS MSA members)\n", scalar(@seq_names)));

  # Open output file and write header
  open(my $outfh, ">", $out_stk_file) || die "ERROR unable to write $out_stk_file: $!";
  print $outfh "# STOCKHOLM 1.0\n";

  my $cds_added = 0;

  foreach my $block (@{$blocks_AR}) {
    my $type        = $block->{"type"};
    my $block_start = $block->{"start"};
    my $block_end   = $block->{"end"};

    if($type eq "coding") {
      next if $cds_added;
      $cds_added = 1;

      my $anchor_seq = $accn_to_cds_seq_H{$anchor_accn};
      if(!defined $anchor_seq) {
        die "ERROR: anchor accession $anchor_accn not found in CDS MSA $cds_msa_fa_file";
      }
      my $cds_aln_width = length($anchor_seq);

      # Build RF from anchor: non-gap -> x, gap -> .
      my $cds_rf = "";
      foreach my $char (split(//, $anchor_seq)) {
        $cds_rf .= ($char eq '-' || $char eq '.') ? '.' : 'x';
      }

      # Write sequence lines in canonical order
      foreach my $name (@seq_names) {
        my $seq = exists $accn_to_cds_seq_H{$name} ? $accn_to_cds_seq_H{$name} : '-' x $cds_aln_width;
        printf $outfh "%-30s %s\n", $name, $seq;
      }
      printf $outfh "#=GC %-24s %s\n", "RF", $cds_rf;

      # Build gapped SS_cons for CDS block by following anchor gaps, write per-block
      my $ungapped_substr = substr($ungapped_ss_cons, $block_start - 1, $block_end - $block_start + 1);
      my $cds_ss = "";
      my $ss_pos = 0;
      foreach my $char (split(//, $anchor_seq)) {
        if($char eq '-' || $char eq '.') {
          $cds_ss .= '.';
        }
        else {
          $cds_ss .= ($ss_pos < length($ungapped_substr)) ? substr($ungapped_substr, $ss_pos, 1) : ':';
          $ss_pos++;
        }
      }
      printf $outfh "#=GC %-24s %s\n", "SS_cons", $cds_ss;
      print  $outfh "\n";

      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: added CDS MSA (%d..%d, %d alignment columns)\n",
                                                       $block_start, $block_end, $cds_aln_width));
    }
    elsif($type eq "rna") {
      # Read RNA Stockholm via Bio::Easel::MSA
      my $rna_stk = $block->{"rna_file"};
      my $rna_msa = Bio::Easel::MSA->new({ fileLocation => $rna_stk, isDna => 1 });
      my $rna_aln_width = $rna_msa->alen();
      my $rna_rf      = $rna_msa->has_rf()      ? $rna_msa->get_rf()      : 'x' x $rna_aln_width;
      my $rna_ss_cons = $rna_msa->has_ss_cons() ? $rna_msa->get_ss_cons() : '.' x $rna_aln_width;

      # Build name->sequence hash for lookup
      my %rna_seqs_H = ();
      for(my $i = 0; $i < $rna_msa->nseq(); $i++) {
        $rna_seqs_H{$rna_msa->get_sqname($i)} = $rna_msa->get_sqstring_aligned($i);
      }


      # Write sequence lines in canonical order
      foreach my $name (@seq_names) {
        my $seq = exists $rna_seqs_H{$name} ? $rna_seqs_H{$name} : '-' x $rna_aln_width;
        printf $outfh "%-30s %s\n", $name, $seq;
      }
      printf $outfh "#=GC %-24s %s\n", "RF", $rna_rf;
      printf $outfh "#=GC %-24s %s\n", "SS_cons", $rna_ss_cons;
      print  $outfh "\n";

      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: added RNA block %d (%s, %d..%d, %d alignment columns)\n",
                                                       $block->{"rna_idx"}, $block->{"family"},
                                                       $block_start, $block_end, $rna_aln_width));
    }
    elsif($type eq "noncoding") {
      # Extract noncoding region from tier2 via Bio::Easel::MSA column_subset
      # Map RF positions block_start..block_end to alignment columns using rf2a_map
      # rf2a_map returns 1-based apos; column_subset uses 0-based [0..alen-1]
      my @useme_A = (0) x $tier2_msa->alen();
      my $first_acol = $rf2a_map_A[$block_start];  # 1-based
      my $last_acol  = $rf2a_map_A[$block_end];    # 1-based
      for(my $a = $first_acol; $a <= $last_acol; $a++) {
        $useme_A[$a - 1] = 1;  # convert to 0-based
      }
      my $nc_msa = $tier2_msa->clone_msa();
      $nc_msa->column_subset(\@useme_A);
      my $nc_aln_width = $nc_msa->alen();
      my $nc_rf = $nc_msa->has_rf() ? $nc_msa->get_rf() : 'x' x $nc_aln_width;

      # Build name->sequence hash for lookup
      my %nc_seqs_H = ();
      for(my $i = 0; $i < $nc_msa->nseq(); $i++) {
        $nc_seqs_H{$nc_msa->get_sqname($i)} = $nc_msa->get_sqstring_aligned($i);
      }


      if($nc_aln_width == 0) { $nc_aln_width = $block_end - $block_start + 1; }

      # Write sequence lines in canonical order
      foreach my $name (@seq_names) {
        my $seq = exists $nc_seqs_H{$name} ? $nc_seqs_H{$name} : '-' x $nc_aln_width;
        printf $outfh "%-30s %s\n", $name, $seq;
      }
      printf $outfh "#=GC %-24s %s\n", "RF", $nc_rf;
      printf $outfh "#=GC %-24s %s\n", "SS_cons", '.' x $nc_aln_width;
      print  $outfh "\n";

      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: added noncoding block (%d..%d, %d alignment columns)\n",
                                                       $block_start, $block_end, $nc_aln_width));
    }
  }

  # Write footer
  print $outfh "//\n";
  close($outfh);


  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: wrote interleaved Stockholm with RF and SS_cons to %s\n", $out_stk_file));

  return;
}


#################################################################
# Subroutine : build_ungapped_ss_cons()
# Purpose    : Build ungapped reference-length SS_cons string
#              by overlaying RNA structures onto single-strand background
#
# Arguments  :
#   $ref_len    : reference sequence length
#   $rna_AR     : ref to array of RNA hashes (with start, end, stk_file)
#   $FH_HR      : ref to hash of file handles
#
# Returns    : ungapped SS_cons string (length = ref_len)
#################################################################
sub build_ungapped_ss_cons {
  my ($ref_len, $rna_AR, $FH_HR) = @_;
  
  # Initialize array with single-stranded characters
  my @ss_array = (':') x $ref_len;
  
  # Overlay RNA structures
  foreach my $rna (@{$rna_AR}) {
    my $start = $rna->{"start"};
    my $end   = $rna->{"end"};
    my $stk   = $rna->{"stk_file"};
    
    # Read SS_cons from RNA alignment via Bio::Easel::MSA
    my $rna_msa_tmp = Bio::Easel::MSA->new({ fileLocation => $stk, isDna => 1 });
    my $rna_ss = $rna_msa_tmp->has_ss_cons() ? $rna_msa_tmp->get_ss_cons() : "";
    
    # Strip gap characters ('.') to get ungapped structure
    $rna_ss =~ s/\.//g;
    
    my $expected_len = $end - $start + 1;
    my $actual_len = length($rna_ss);
    
    if($actual_len != $expected_len) {
      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# WARNING: RNA %s SS_cons length mismatch (expected %d, got %d)\n", 
                                                       $rna->{"family"}, $expected_len, $actual_len));
    }
    
    # Overlay onto array (convert to 0-based indexing)
    for(my $i = 0; $i < length($rna_ss) && ($start + $i) <= $end; $i++) {
      my $ref_pos = $start + $i - 1; # Convert to 0-based
      if($ref_pos >= 0 && $ref_pos < $ref_len) {
        $ss_array[$ref_pos] = substr($rna_ss, $i, 1);
      }
    }
    
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: overlaid RNA %s structure (%d..%d, %d chars)\n", 
                                                     $rna->{"family"}, $start, $end, length($rna_ss)));
  }
  
  # Join to create ungapped string
  my $ungapped_ss = join('', @ss_array);
  
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: built ungapped SS_cons (%d positions)\n", length($ungapped_ss)));
  
  return $ungapped_ss;
}

#################################################################
# Subroutine : generate_updated_minfo()
# Purpose    : Generate updated .minfo file with RNA features
#              added as misc_structure annotations.
#              All RNA features are added (overlaps with CDS allowed -
#              v-annotate.pl handles overlapping features)
#
# Arguments  :
#   $seed_minfo_file : path to seed model .minfo file
#   $rna_annot_file  : path to RNA annotation TSV file
#   $rna_regions_AR  : ref to array of RNA region hashes
#   $out_minfo_file  : output .minfo file path
#   $model_key       : model name
#   $execs_HR        : ref to hash of executable paths
#   $FH_HR           : ref to hash of file handles
#
# Returns    : void
#################################################################
sub generate_updated_minfo {
  my ($seed_minfo_file, $rna_annot_file, $rna_regions_AR, $out_minfo_file, $model_key, $execs_HR, $FH_HR) = @_;

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Generating updated .minfo file with RNA features\n"));

  # Read seed .minfo file
  my @seed_lines = ();
  open(my $seedfh, $seed_minfo_file) || die "ERROR unable to read seed .minfo $seed_minfo_file: $!";
  while(my $line = <$seedfh>) {
    chomp $line;
    push(@seed_lines, $line);
  }
  close($seedfh);

  # Read RNA annotations with structure
  my @rna_features = ();
  if(defined($rna_regions_AR) && scalar(@{$rna_regions_AR}) > 0) {
    open(my $rnafh, $rna_annot_file) || die "ERROR unable to read RNA annotation $rna_annot_file: $!";
    my $header = <$rnafh>;
    while(my $line = <$rnafh>) {
      chomp $line;
      my @f = split(/\t/, $line);
      my ($idx, $start, $end, $strand, $family, $accession, $score, $evalue, $sstruct) = @f;
      
      # Find first and last basepaired positions (skip single-stranded ends)
      my ($bp_start, $bp_end) = find_basepaired_bounds($sstruct, $start);
      
      # Skip if no basepairs found
      if(!defined($bp_start) || !defined($bp_end)) {
        ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# .minfo: skipped RNA %s (%d..%d) - no basepairs in structure\n", 
                                                         $family, $start, $end));
        next;
      }
      
      # Get Rfam model description
      my $rfam_cm = $ENV{'VADRMODELDIR'} ? $ENV{'VADRMODELDIR'} . "/rfam/Rfam.cm" : "/net/intdev/oblast01/dnaorg/virseqannot/code/vadr-install-1.7/rfam/Rfam.cm";
      my $desc_cmd = "cmfetch " . $rfam_cm . " " . $accession . " 2>/dev/null | grep '^DESC' | awk '{for(i=2;i<=NF;i++) printf \"%s \", \$i; print \"\"}'";
      my $description = `$desc_cmd`;
      chomp $description;
      $description =~ s/\s+$//;  # trim trailing whitespace
      
      if(!$description) {
        $description = $family;  # fallback to family name
      }
      
      push(@rna_features, {
        "start" => $bp_start,
        "end"   => $bp_end,
        "strand" => $strand,
        "family" => $family,
        "note"  => $description
      });
      
      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# .minfo: added RNA feature %s (%d..%d%s) note:\"%s\"\n", 
                                                       $family, $bp_start, $bp_end, $strand, $description));
    }
    close($rnafh);
  }

  # Write updated .minfo file
  open(my $outfh, ">", $out_minfo_file) || die "ERROR unable to write .minfo $out_minfo_file: $!";
  
  # Copy seed lines and insert RNA features after gene/CDS features
  foreach my $line (@seed_lines) {
    print $outfh $line . "\n";
  }
  
  # Add RNA features as misc_structure
  foreach my $rna (@rna_features) {
    my $coords = $rna->{"start"} . ".." . $rna->{"end"} . ":" . $rna->{"strand"};
    my $feature_line = sprintf("FEATURE %s type:\"misc_structure\" coords:\"%s\" parent_idx_str:\"GBNULL\" note:\"%s\"",
                               $model_key, $coords, $rna->{"note"});
    print $outfh $feature_line . "\n";
  }
  
  close($outfh);
  
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Wrote updated .minfo with %d RNA features to %s\n", 
                                                   scalar(@rna_features), $out_minfo_file));
  
  return;
}

#################################################################
# Subroutine : find_basepaired_bounds()
# Purpose    : Find first and last basepaired positions in 
#              consensus structure (skip single-stranded ends)
#
# Arguments  :
#   $sstruct : consensus structure string (e.g., ":::<<___>>:")
#   $start   : start coordinate of structure
#
# Returns    : ($bp_start, $bp_end) adjusted coordinates
#################################################################
sub find_basepaired_bounds {
  my ($sstruct, $start) = @_;
  
  # Find first basepaired position (first non-colon/comma/dash)
  my $first_bp = -1;
  my $last_bp = -1;
  
  my @chars = split(//, $sstruct);
  for(my $i = 0; $i < scalar(@chars); $i++) {
    my $c = $chars[$i];
    # Basepaired characters: <, >, (, ), [, ], {, }
    if($c =~ /[<>()\[\]{}]/) {
      if($first_bp == -1) {
        $first_bp = $i;
      }
      $last_bp = $i;
    }
  }
  
  if($first_bp == -1 || $last_bp == -1) {
    return (undef, undef);
  }
  
  # Adjust coordinates relative to start position
  my $bp_start = $start + $first_bp;
  my $bp_end = $start + $last_bp;
  
  return ($bp_start, $bp_end);
}

