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
my $env_vadr_muscle_dir   = $ENV{"VADRMUSCLEDIR"};  # muscle 3.8.31 (public domain)
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
if(defined $env_vadr_muscle_dir) {
  my $muscle_exec = $env_vadr_muscle_dir . "/muscle";
  if(-x $muscle_exec) {
    $execs_H{"muscle"} = $muscle_exec;
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
opt_Add("--refaccn",    "string",  undef,      $g,    undef, undef,       "reference accession is <s> (default: model key)",             "reference accession is <s> (default: model key)", \%opt_HH, \@opt_order_A);
opt_Add("--taxid",      "string",  undef,      $g,    undef, "--meta",    "fetch metadata from NCBI for taxonomy ID <s>",                "fetch metadata from NCBI for taxonomy ID <s>", \%opt_HH, \@opt_order_A);
opt_Add("--meta",       "string",  undef,      $g,    undef, "--taxid",   "read metadata TSV from <s> instead of fetching",              "read metadata TSV from <s> instead of fetching", \%opt_HH, \@opt_order_A);
opt_Add("--api_key",    "string",  undef,      $g,    "--taxid", undef,   "NCBI API key to use for metadata fetch",                      "NCBI API key to use for metadata fetch as <s>", \%opt_HH, \@opt_order_A);
opt_Add("--seed-accn",  "string",  undef,      $g,    undef, undef,       "build seed model internally from accession <s>",              "build seed model internally from accession <s>", \%opt_HH, \@opt_order_A);
opt_Add("--seed-build-opts-file", "string", undef,   $g,    "--seed-accn", "--seed-build-opts", "read extra v-build.pl options from file <s>",     "read extra v-build.pl options from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--seed-build-opts", "string", undef,        $g,    "--seed-accn", "--seed-build-opts-file", "extra v-build.pl options from single string <s>", "extra v-build.pl options from single string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--skip-annotate", "boolean", 0,             $g,    undef, undef,   "skip Tier 2 v-annotate screening step",                     "skip Tier 2 v-annotate screening step (temporary)", \%opt_HH, \@opt_order_A);
opt_Add("--xambig",     "integer", 5,          $g,    undef, undef,       "max ambiguous nucleotides allowed per sequence",             "max ambiguous nucleotides allowed per sequence as <n>", \%opt_HH, \@opt_order_A);
opt_Add("--xpergroup",  "integer", 50,         $g,    undef, undef,       "max sequences retained per serotype/genotype group",         "max sequences retained per serotype/genotype group as <n>", \%opt_HH, \@opt_order_A);
opt_Add("--group",      "string",  undef,      $g,    undef, undef,       "virus group name for #=GS GP annotations in output .stk",   "virus group name for #=GS GP annotations in output .stk as <s>", \%opt_HH, \@opt_order_A);
opt_Add("--rna-cm-file", "string", undef,        $g,    undef, "--skip-rna", "use CM file <s> for RNA search instead of default Rfam.cm",  "use CM file <s> for RNA search instead of default Rfam.cm", \%opt_HH, \@opt_order_A);
opt_Add("--skip-rna",   "boolean", 0,             $g,    undef, "--rna-cm-file", "skip RNA discovery and alignment",                          "skip RNA discovery and alignment (treat all noncoding as unstructured)", \%opt_HH, \@opt_order_A);
opt_Add("--no-auto-alt", "boolean", 0,           $g,    undef, "--alt-file",   "skip auto-detection of alternative CDS features and exceptions", "skip auto-detection of alternative CDS features and exceptions", \%opt_HH, \@opt_order_A);
opt_Add("--alt-file",   "string",  undef,         $g,    undef, "--no-auto-alt", "read alternative feature definitions from file <s>",           "read alternative feature definitions from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--alt-min-ind", "integer", 2,            $g,    undef, "--no-auto-alt", "min independent observations to add an alternative or exception", "min independent observations to add an alternative or exception as <n>", \%opt_HH, \@opt_order_A);
opt_Add("--vannot-opts-file", "string", undef,    $g,    undef, undef,          "read extra v-annotate.pl options from file <s>",              "read extra v-annotate.pl options from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--alt-max-fract", "real",    0.2,       $g,    undef, "--no-auto-alt", "max fractional length deviation for alternative CDS",         "max fractional length deviation for alternative CDS as <x>", \%opt_HH, \@opt_order_A);
opt_Add("--nper1grp",     "integer", 5,         $g,    undef, undef,          "number of seqs per group when 1 group",                       "number of seqs per group when 1 group as <n>", \%opt_HH, \@opt_order_A);

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
                'refaccn=s'    => \$GetOptions_H{"--refaccn"},
                'taxid=s'      => \$GetOptions_H{"--taxid"},
                'meta=s'       => \$GetOptions_H{"--meta"},
                'api_key=s'    => \$GetOptions_H{"--api_key"},
                'seed-accn=s'  => \$GetOptions_H{"--seed-accn"},
                'seed-build-opts-file=s' => \$GetOptions_H{"--seed-build-opts-file"},
                'seed-build-opts=s' => \$GetOptions_H{"--seed-build-opts"},
                'skip-annotate' => \$GetOptions_H{"--skip-annotate"},
                'xambig=i'     => \$GetOptions_H{"--xambig"},
                'xpergroup=i'  => \$GetOptions_H{"--xpergroup"},
                'group=s'      => \$GetOptions_H{"--group"},
                'rna-cm-file=s' => \$GetOptions_H{"--rna-cm-file"},
                'skip-rna'     => \$GetOptions_H{"--skip-rna"},
                'no-auto-alt'  => \$GetOptions_H{"--no-auto-alt"},
                'alt-file=s'   => \$GetOptions_H{"--alt-file"},
                'alt-min-ind=i' => \$GetOptions_H{"--alt-min-ind"},
                'vannot-opts-file=s' => \$GetOptions_H{"--vannot-opts-file"},
                'alt-max-fract=f' => \$GetOptions_H{"--alt-max-fract"},
                'nper1grp=i'   => \$GetOptions_H{"--nper1grp"},
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
my $group_name = opt_Get("--group", \%opt_HH);
if(! defined $group_name) { $group_name = $dir_tail; } # auto-guess from output dir basename

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

# open the log, command, and list files
my %ofile_info_HH = ();
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "log",  $out_root . ".log",      1, 1, "Output printed to screen");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "cmd",  $out_root . ".cmd",      1, 1, "List of executed commands");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".filelist", 1, 1, "List and description of all output files");
my $FH_HR  = $ofile_info_HH{"FH"};
my $log_FH = $FH_HR->{"log"};
my $cmd_FH = $FH_HR->{"cmd"};

# Output early commands to .cmd
foreach $cmd (@early_cmd_A) {
  print $cmd_FH $cmd . "\n";
}

ofile_OutputBanner($log_FH, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

#---------------------------------------
# Step 1: Optional seed bootstrap via v-build.pl
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
  utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
}

if(! -e $seed_minfo) {
  die "ERROR: expected seed model .minfo file $seed_minfo does not exist (derived from model directory and model key; use --mkey to override)";
}

#---------------------------------------
# Step 2: Acquire metadata TSV if needed
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
  utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
  $meta_tsv = $fetch_out_prefix . ".tsv";
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "metadata.tsv", $meta_tsv, 1, 1, "sequence metadata fetched from NCBI for taxid");
}

if((! defined $meta_tsv) || (! -e $meta_tsv)) {
  die "ERROR, metadata TSV is undefined or does not exist: " . ((defined $meta_tsv) ? $meta_tsv : "[undef]");
}

#---------------------------------------
# Read seed model length
#---------------------------------------
my @mdl_info_A = ();
my %ftr_info_HA = ();
my @reqd_mdl_keys = ("length");
my @reqd_ftr_keys = ();

vdr_ModelInfoFileParse($seed_minfo, \@reqd_mdl_keys, \@reqd_ftr_keys, \@mdl_info_A, \%ftr_info_HA, $FH_HR);
my $seed_model_len = $mdl_info_A[0]{"length"};
my $ref_accn = opt_IsUsed("--refaccn", \%opt_HH) ? opt_Get("--refaccn", \%opt_HH) : $model_key; # reference accession (default: model key, e.g. "NC_006232")
ofile_OutputString(*STDOUT, 1, sprintf("# Read seed model length: %d\n", $seed_model_len));
ofile_OutputString(*STDOUT, 1, sprintf("# Reference accession: %s\n", $ref_accn));

# Read extra v-annotate.pl options from file, if provided
my $vannot_extra_opts = "";
if(opt_IsUsed("--vannot-opts-file", \%opt_HH)) {
  my $vannot_opts_file = opt_Get("--vannot-opts-file", \%opt_HH);
  open(my $vofh, "<", $vannot_opts_file) || ofile_FAIL("ERROR, unable to open --vannot-opts-file $vannot_opts_file", 1, $FH_HR);
  while(my $voline = <$vofh>) {
    chomp $voline;
    $voline =~ s/^\s+//;
    $voline =~ s/\s+$//;
    next if($voline eq "" || $voline =~ /^\#/);
    $vannot_extra_opts .= " " . $voline;
  }
  close($vofh);
  ofile_OutputString(*STDOUT, 1, sprintf("# Extra v-annotate.pl options:%s\n", $vannot_extra_opts));
}

#---------------------------------------
# Step 3: RNA discovery via cmscan on reference sequence
#---------------------------------------
my @rna_regions_A = (); # Array of hashes: { start, end, strand, cm_family, cm_accession, score, evalue }
my $rna_annot_file = $out_root . ".rna_annotation.tsv";
my $rna_ss_cons = undef; # Full-length consensus secondary structure string
my $do_keep = opt_Get("--keep", \%opt_HH);
my @to_remove_A = (); # files to remove at end unless --keep

if($do_rna_discovery) {
  my $ref_seq_file = $model_dir . "/" . $model_key . ".vadr.fa";
  if(! -e $ref_seq_file) {
    die "ERROR, reference sequence file not found for RNA discovery: $ref_seq_file";
  }
  
  run_rna_discovery($ref_seq_file, $rna_cm_file, $seed_model_len, $out_root, $env_vadr_rfam_dir,
                    \@rna_regions_A, \$rna_ss_cons, $do_keep, opt_Get("-v", \%opt_HH),
                    \%execs_H, \%ofile_info_HH, \@to_remove_A, $FH_HR);

  # Step 3b: RNA structure generation via cmalign
  run_rna_sstruct_generation(\@rna_regions_A, $ref_seq_file, $rna_cm_file, $env_vadr_rfam_dir,
                             $out_root, $do_keep, opt_Get("-v", \%opt_HH),
                             \%execs_H, \%ofile_info_HH, $FH_HR);

  # Write RNA annotation output (after structure extraction)
  write_rna_annotation_table(\@rna_regions_A, $rna_annot_file, \%ofile_info_HH, $FH_HR);
}
else {
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA discovery: skipped due to --skip-rna\n"));
}

#---------------------------------------
# Step 4: Tier 1 metadata filtering
#---------------------------------------

# %candidate_AH: arrays of hashes [1..nseq-1] grouped by serotype/genotype keys.
# Each hash contains sequence metadata:
#   acc      => sequence accession
#   len      => sequence length
#   cdate    => genbank create date
#   serotype => sequence serotype
#   genotype => sequence genotype
#   isolate  => sequence isolate
# Resolve reference accession versioning BEFORE filtering.
# The model key may be unversioned (e.g. "NC_006232") while the metadata
# uses versioned accessions (e.g. "NC_006232.1"). Scan the metadata header
# to find the versioned match so all downstream filters use the correct accession.
{
  open(my $pre_fh, "<", $meta_tsv) || die "ERROR, unable to read metadata TSV $meta_tsv: $!";
  my $pre_hdr = <$pre_fh>; # skip header
  while(my $line = <$pre_fh>) {
    chomp $line;
    my ($acc) = split(/\t/, $line);
    if(defined $acc && $acc =~ /^\Q$ref_accn\E\.\d+$/) {
      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Reference accession %s matched versioned accession %s in metadata\n", $ref_accn, $acc));
      $ref_accn = $acc;
      last;
    }
    elsif(defined $acc && $acc eq $ref_accn) {
      last; # exact match, no versioning needed
    }
  }
  close($pre_fh);
}

my %candidate_AH = ();
my %decision_H = ();
my $max_per_group = opt_Get("--xpergroup", \%opt_HH);
parse_and_filter_metadata($meta_tsv, $seed_model_len, $max_per_group, \%candidate_AH, \%decision_H, $ref_accn, $FH_HR);

# Verify reference accession was found in metadata
if(! exists $decision_H{$ref_accn}) {
  ofile_FAIL("ERROR, reference accession $ref_accn not found in metadata TSV $meta_tsv.\n" .
             "The reference sequence must be present in the metadata for v-prep.pl to proceed.\n" .
             "If this accession is not returned by the NCBI taxonomy fetch, you can add it\n" .
             "manually to the metadata TSV file and re-run with --meta.", 1, $FH_HR);
}

#---------------------------------------
# Step 5: Tier 2 sequence fetch and filtering
#---------------------------------------
my $tier1_accn_file = $out_root . ".tier1.accn.list";
my $tier2_fasta_file = $out_root . ".tier2.fa";
my $max_ambig_nt = opt_Get("--xambig", \%opt_HH);
my $tier2_ant_outdir = $out_root . ".tier2.annot";
my $centroid_tsv_file = $out_root . ".centroid.tsv";
my $decision_tsv_file = $out_root . ".filter.seq.tsv";
my $decision_summary_tsv_file = $out_root . ".filter.sum.tsv";
my $stitch_selected_accn_file = $out_root . ".cds.selected.accn.list";
my $stitch_selected_fa_file   = $out_root . ".cds.selected.fa";
my $stitch_block_plan_file    = $out_root . ".block_plan.tsv";
my $stitch_cds_nt_fa_file     = $out_root . ".cds.nt.fa";
my $stitch_cds_orf_fa_file    = $out_root . ".cds.orf.fa";
my $stitch_cds_aa_fa_file     = $out_root . ".cds.aa.fa";
my $stitch_cds_map_tsv_file   = $out_root . ".cds.translate_map.tsv";
my $stitch_cds_anchor_tsv_file = $out_root . ".cds.anchor.tsv";
my $stitch_cds_anchor_fa_file  = $out_root . ".cds.anchor.fa";
my $stitch_cds_pairwise_tsv_file = $out_root . ".cds.pairwise.tsv";
my $stitch_cds_muscle_dir     = $out_root . ".cds.muscle.dir";
my $stitch_cds_msa_aa_fa_file = $out_root . ".cds.msa.aa.afa";
my $stitch_cds_msa_aa_stk_file = $out_root . ".cds.msa.aa.stk";
my $stitch_cds_msa_nt_fa_file = $out_root . ".cds.msa.nt.afa";
my $rna_annotation_file       = $out_root . ".rna_annotation.tsv";

write_accession_list_from_candidates(\%candidate_AH, $tier1_accn_file, $do_keep, \%ofile_info_HH, \@to_remove_A, $FH_HR);
fetch_fasta_from_accession_list($tier1_accn_file, $tier2_fasta_file, $do_keep, \%ofile_info_HH, \@to_remove_A, \%opt_HH, $FH_HR);
apply_ambiguity_filter_to_candidates(\%candidate_AH, $tier2_fasta_file, $max_ambig_nt, \%decision_H, $ref_accn, $FH_HR);

my $tier2_align_stk_file = undef;
if($do_skip_annotate) {
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Tier 2 v-annotate filter: skipped due to --skip-annotate\n"));
  mark_all_remaining_as_selected_for_tier3(\%candidate_AH, \%decision_H);
}
else {
  $tier2_align_stk_file = run_vannotate_filter_fails(\%candidate_AH, $tier2_fasta_file, $tier2_ant_outdir, $model_dir, $model_key, $do_keep, opt_Get("-v", \%opt_HH), \%execs_H, \%decision_H, $vannot_extra_opts, $FH_HR);
}

#---------------------------------------
# Step 5b: Auto-detect alternative CDS features and exceptions
#---------------------------------------
my $do_auto_alt = (! $do_skip_annotate) && (! opt_Get("--no-auto-alt", \%opt_HH));
if($do_auto_alt) {
  my $tier2_ant_tail = $tier2_ant_outdir;
  $tier2_ant_tail =~ s/^.+\///;
  my $tier2_alt_file = $tier2_ant_outdir . "/" . $tier2_ant_tail . ".vadr.alt";
  my $tier2_ftr_file = $tier2_ant_outdir . "/" . $tier2_ant_tail . ".vadr.ftr";
  my $min_independent = opt_Get("--alt-min-ind", \%opt_HH);

  if(-e $tier2_alt_file && -e $tier2_ftr_file) {
    # Detect alternative CDS features
    my $alt_groups_HHR = parse_alt_for_cds_boundary_alerts($tier2_alt_file, $FH_HR);
    my $max_fract_diff = opt_Get("--alt-max-fract", \%opt_HH);
    my $alt_features_AR = detect_alternative_features($alt_groups_HHR, $min_independent, $max_fract_diff,
                                                       \@{$ftr_info_HA{$model_key}}, $model_key, $FH_HR);

    # Detect exceptions
    my $exc_groups_HHR = parse_alt_for_exceptions($tier2_alt_file, $FH_HR);
    my $exceptions_AR = detect_exceptions($exc_groups_HHR, $min_independent, $FH_HR);

    my $n_alt = scalar(@{$alt_features_AR});
    my $n_exc = scalar(@{$exceptions_AR});
    ofile_OutputString($FH_HR->{"log"}, 1,
      sprintf("# Auto-alt detect: found %d alternative features and %d exceptions\n", $n_alt, $n_exc));

    if($n_alt > 0 || $n_exc > 0) {
      # Update the minfo file
      my $updated_minfo = $out_root . ".alt.minfo";
      add_alternatives_and_exceptions_to_minfo($seed_minfo, $updated_minfo,
                                                $alt_features_AR, $exceptions_AR,
                                                $model_key, $FH_HR);
      ofile_OutputString($FH_HR->{"log"}, 1,
        sprintf("# Auto-alt detect: wrote updated minfo to %s\n", $updated_minfo));

      # Translate alternative CDS proteins and update BLAST db
      if($n_alt > 0) {
        # Copy seed protein.fa to output dir so we can append to it
        my $updated_protein_fa = $out_root . ".alt.protein.fa";
        my $seed_protein_fa = $model_dir . "/" . $model_key . ".vadr.protein.fa";
        if(-e $seed_protein_fa) {
          utl_RunCommand("cp $seed_protein_fa $updated_protein_fa", opt_Get("-v", \%opt_HH), 0, $FH_HR);
          translate_alternative_cds_proteins($alt_features_AR, $tier2_ftr_file, $tier2_alt_file,
                                              $tier2_fasta_file, $updated_protein_fa,
                                              $model_key, \%execs_H, \%opt_HH, $FH_HR);
        }
      }

      # Selective v-annotate re-run: identify sequences that failed
      # ONLY due to alerts addressed by the new alternatives/exceptions,
      # and re-run v-annotate.pl on just those sequences
      my @rerun_accessions = identify_rerun_candidates($tier2_alt_file,
                                                        \%decision_H,
                                                        $alt_features_AR,
                                                        $exceptions_AR,
                                                        $FH_HR);

      if(scalar(@rerun_accessions) > 0) {
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# Auto-alt detect: re-running v-annotate on %d candidate sequences\n", scalar(@rerun_accessions)));

        # Create a FASTA subset for the re-run candidates
        my $rerun_fa = $out_root . ".alt.rerun.fa";
        my $rerun_accn_file = $out_root . ".alt.rerun.accn";
        open(my $accn_fh, ">", $rerun_accn_file) || ofile_FAIL("ERROR, unable to write $rerun_accn_file", 1, $FH_HR);
        foreach my $acc (@rerun_accessions) { print $accn_fh "$acc\n"; }
        close($accn_fh);
        # Index the tier2 FASTA for esl-sfetch if not already indexed
        if(! -e $tier2_fasta_file . ".ssi") {
          utl_RunCommand($execs_H{"esl-sfetch"} . " --index $tier2_fasta_file",
                         opt_Get("-v", \%opt_HH), 0, $FH_HR);
        }
        utl_RunCommand($execs_H{"esl-sfetch"} . " -f $tier2_fasta_file $rerun_accn_file > $rerun_fa",
                       opt_Get("-v", \%opt_HH), 0, $FH_HR);

        # Set up temp model dir with updated minfo + model files
        my $tmp_mdir = $out_root . ".alt.mdir";
        if(! -d $tmp_mdir) { utl_RunCommand("mkdir $tmp_mdir", 0, 0, $FH_HR); }
        utl_RunCommand("cp $updated_minfo $tmp_mdir/$model_key.vadr.minfo", 0, 0, $FH_HR);
        # Symlink CM and other model files
        my @model_exts = (".vadr.cm", ".vadr.cm.i1f", ".vadr.cm.i1i", ".vadr.cm.i1m", ".vadr.cm.i1p",
                          ".vadr.fa", ".vadr.fa.ssi",
                          ".vadr.fa.ndb", ".vadr.fa.nhr", ".vadr.fa.nin", ".vadr.fa.njs",
                          ".vadr.fa.not", ".vadr.fa.nsq", ".vadr.fa.ntf", ".vadr.fa.nto");
        foreach my $ext (@model_exts) {
          my $src = $model_dir . "/" . $model_key . $ext;
          my $dst = $tmp_mdir . "/" . $model_key . $ext;
          if(-e $src && ! -e $dst) { utl_RunCommand("ln -s $src $dst", 0, 0, $FH_HR); }
        }
        # Copy protein.fa into temp mdir and rebuild BLAST db
        # (must rebuild rather than copy db files because BLAST db
        # internally references the source filename)
        my $updated_protein_fa = $out_root . ".alt.protein.fa";
        my $src_protein_fa = (-e $updated_protein_fa) ? $updated_protein_fa :
                             $model_dir . "/" . $model_key . ".vadr.protein.fa";
        if(-e $src_protein_fa) {
          my $dst_protein_fa = $tmp_mdir . "/" . $model_key . ".vadr.protein.fa";
          utl_RunCommand("cp $src_protein_fa $dst_protein_fa", 0, 0, $FH_HR);
          sqf_BlastDbCreate($execs_H{"makeblastdb"}, "prot", $dst_protein_fa, \%opt_HH, $FH_HR);
        }

        # Run v-annotate.pl on the subset
        my $rerun_outdir = $out_root . ".vadr.tier2.annot.pass2";
        my $rerun_mkey = $model_key . ".vadr";
        my $cmd = $execs_H{"v-annotate.pl"} . " -f --mdir " . $tmp_mdir . " --mkey " . $rerun_mkey .
                  $vannot_extra_opts . " --out_stk " . $rerun_fa . " " . $rerun_outdir;
        if(! opt_Get("-v", \%opt_HH)) { $cmd .= " > /dev/null"; }
        utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);

        # Read pass2 fail list and update decisions
        my $rerun_outdir_tail = $rerun_outdir;
        $rerun_outdir_tail =~ s/^.+\///;
        my $pass2_fail_file = $rerun_outdir . "/" . $rerun_outdir_tail . ".vadr.fail.list";
        my %pass2_fail_H = ();
        if(-e $pass2_fail_file) {
          open(my $ffh, "<", $pass2_fail_file) || ofile_FAIL("ERROR, unable to read $pass2_fail_file", 1, $FH_HR);
          while(my $l = <$ffh>) { chomp $l; $l =~ s/^\s+//; $l =~ s/\s+$//; next if($l eq ""); my ($n) = split(/\s+/, $l); $pass2_fail_H{$n} = 1; }
          close($ffh);
        }

        my $n_rescued = 0;
        foreach my $acc (@rerun_accessions) {
          if(! exists $pass2_fail_H{$acc}) {
            # This sequence now passes — restore it to the candidate pool
            $decision_H{$acc}{"status"} = "kept";
            $decision_H{$acc}{"reason_code"} = "pass_alt";
            $decision_H{$acc}{"reason_detail"} = "rescued by auto-alt second pass";
            $decision_H{$acc}{"stage_last_seen"} = "selected_for_tier3";
            # Re-add to candidate_AH (need to find its group)
            # The sequence was removed from candidate_AH by run_vannotate_filter_fails
            # We need its group info from decision_H or the original candidate data
            $n_rescued++;
          }
        }
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# Auto-alt detect: pass2 rescued %d of %d re-run candidates\n",
                  $n_rescued, scalar(@rerun_accessions)));
      }
      else {
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# Auto-alt detect: no sequences eligible for re-run\n"));
      }

      # Update seed_minfo to the updated version for downstream steps
      $seed_minfo = $updated_minfo;

      # Re-parse feature info from updated minfo
      vdr_ModelInfoFileParse($seed_minfo, \@reqd_mdl_keys, \@reqd_ftr_keys,
                             \@mdl_info_A, \%ftr_info_HA, $FH_HR);
    }
  }
}

#---------------------------------------
# Step 6: Tier 3 centroid selection (BLAST all-vs-all)
#---------------------------------------
# Identify sequences with partial CDS coordinates so they are not chosen as
# centroids when non-partial alternatives exist in the same group.
my %partial_cds_H = ();
if(! $do_skip_annotate) {
  my $tier2_ant_tail = $tier2_ant_outdir;
  $tier2_ant_tail =~ s/^.+\///;
  my $tier2_ftr_file = $tier2_ant_outdir . "/" . $tier2_ant_tail . ".vadr.ftr";
  if(-e $tier2_ftr_file) {
    %partial_cds_H = get_partial_cds_accns_from_ftr($tier2_ftr_file);
  }
}
# Determine how many sequences to select per group based on number of non-empty groups
my $n_nonempty_groups = 0;
foreach my $group (keys %candidate_AH) {
  if(scalar(@{$candidate_AH{$group}}) > 0) { $n_nonempty_groups++; }
}
my $nper1grp = opt_Get("--nper1grp", \%opt_HH);
my $n_per_group = determine_seqs_per_group($n_nonempty_groups, $nper1grp);
ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Centroid selection: %d non-empty groups, selecting %d seqs per group\n", $n_nonempty_groups, $n_per_group));
select_group_centroids_blast(\%candidate_AH, $tier2_fasta_file, $out_root, $centroid_tsv_file, $n_per_group, $do_keep, opt_Get("-v", \%opt_HH), \%execs_H, \%decision_H, \%partial_cds_H, $ref_accn, \%ofile_info_HH, \@to_remove_A, $FH_HR);

write_decision_report(\%decision_H, $decision_tsv_file, \%ofile_info_HH, $FH_HR);
if($do_keep) {
  write_decision_stage_reports(\%decision_H, $out_root . ".filter", $do_keep, \%ofile_info_HH, \@to_remove_A, $FH_HR);
}
write_decision_summary_report(\%decision_H, $decision_summary_tsv_file, \%ofile_info_HH, $FH_HR);

#---------------------------------------
# Step 7: Initial piecewise stitching scaffold outputs
#---------------------------------------
my $n_selected = write_stitch_scaffold_outputs(\%candidate_AH, $tier2_fasta_file, \%ftr_info_HA, $model_key, $seed_model_len, $stitch_selected_accn_file, $stitch_selected_fa_file, $stitch_block_plan_file, $do_keep, \%ofile_info_HH, \@to_remove_A, $FH_HR);

if($n_selected == 0) {
  ofile_OutputString($FH_HR->{"log"}, 1, "#\n# Zero sequences passed all filters. Cannot build profile alignment.\n");
  ofile_OutputString($FH_HR->{"log"}, 1, "# Check decision reports for details on why sequences were removed.\n");
  ofile_OutputString($FH_HR->{"log"}, 1, "#\n# Exiting with error.\n");
  ofile_FAIL("ERROR, zero sequences passed all filters, cannot build profile alignment. See $decision_summary_tsv_file for details.", 1, $FH_HR);
}

#---------------------------------------
# Step 8: CDS translation prep for protein alignment
#---------------------------------------
prepare_cds_translation_for_stitching(\%candidate_AH, $stitch_selected_fa_file, $tier2_ant_outdir, $stitch_cds_nt_fa_file, $stitch_cds_orf_fa_file, $stitch_cds_aa_fa_file, $stitch_cds_map_tsv_file, $do_skip_annotate, $do_keep, opt_Get("-v", \%opt_HH), \%execs_H, \%ofile_info_HH, \@to_remove_A, $FH_HR);

#---------------------------------------
# Step 9: Multiple AA alignment with muscle (per CDS feature)
#---------------------------------------
run_muscle_aa_alignment($stitch_cds_aa_fa_file, $stitch_cds_map_tsv_file, $ref_accn,
                        $stitch_cds_anchor_tsv_file, $stitch_cds_muscle_dir,
                        $do_skip_annotate, $do_keep, opt_Get("-v", \%opt_HH),
                        \%execs_H, \%ofile_info_HH, \@to_remove_A, $FH_HR);

#---------------------------------------
# Step 10: Build CDS MSA from muscle alignment, backconvert to NT
#---------------------------------------
build_muscle_cds_msa($stitch_cds_aa_fa_file,
                     $stitch_cds_nt_fa_file,
                     $stitch_cds_map_tsv_file,
                     $stitch_cds_anchor_tsv_file,
                     $stitch_cds_muscle_dir,
                     $stitch_cds_msa_aa_fa_file,
                     $stitch_cds_msa_aa_stk_file,
                     $stitch_cds_msa_nt_fa_file,
                     $do_skip_annotate,
                     $do_keep, \%ofile_info_HH, \@to_remove_A, $FH_HR);


#---------------------------------------
# Step 11: RNA region extraction and alignment refinement
#---------------------------------------
if($do_rna_discovery && (scalar(@rna_regions_A) > 0) && (!$do_skip_annotate)) {
  my $rna_struct_dir = $out_root . ".rna_struct";
  extract_and_align_rna_regions(\@rna_regions_A, $tier2_align_stk_file, $rna_struct_dir,
                                $out_root, $do_keep, opt_Get("-v", \%opt_HH),
                                \%execs_H, \%ofile_info_HH, \@to_remove_A, $FH_HR);
}

#---------------------------------------
#---------------------------------------
# Step 12: Stitch all blocks into final training alignment
#---------------------------------------
my $final_stk_file = $out_root . ".final.stk";
my $temp_cm_file = $out_root . ".temp.cm";

stitch_and_refine_final_alignment($stitch_block_plan_file,
                                  $rna_annotation_file,
                                  $tier2_align_stk_file,
                                  $stitch_cds_msa_nt_fa_file,
                                  $out_root,
                                  \@rna_regions_A,
                                  $final_stk_file,
                                  $temp_cm_file,
                                  $do_rna_discovery,
                                  $do_skip_annotate,
                                  $seed_model_len,
                                  \%execs_H,
                                  $do_keep, \%ofile_info_HH, \@to_remove_A, $FH_HR);

#---------------------------------------
# Step 12b: Add #=GS GP/SG group/subgroup annotations to final alignment
#---------------------------------------
my $output_stk_file = $out_root . ".stk";
annotate_stk_group_subgroup($output_stk_file, $centroid_tsv_file, $group_name, $FH_HR);

#---------------------------------------
# Step 13: Generate updated .minfo file with RNA features
#---------------------------------------
if($do_rna_discovery && !$do_skip_annotate && scalar(@rna_regions_A) > 0) {
  my $updated_minfo_file = $out_root . ".minfo";
  generate_updated_minfo($seed_minfo, $rna_annotation_file, \@rna_regions_A,
                         $updated_minfo_file, $model_key, \%execs_H, \%ofile_info_HH, $FH_HR);
}

if(! $do_keep) {
  utl_FileRemoveList(\@to_remove_A, "v-prep.pl", \%opt_HH, $FH_HR);
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
  my ($tsv_file, $seed_model_len, $max_per_group, $candidate_AHR, $decision_HR, $ref_accn, $FH_HR) = @_;

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
    # Reference accession is always kept regardless of length
    if ($seed_model_len > 0 && $len < ($seed_model_len * 0.90) && $acc ne $ref_accn) {
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

    # Force-include reference accession even if not chronologically selected
    if(defined $ref_accn) {
      my $ref_in_group = 0;
      foreach my $seq (@seqs) {
        if($seq->{acc} eq $ref_accn) { $ref_in_group = 1; last; }
      }
      if($ref_in_group && !$is_selected{$ref_accn}) {
        foreach my $seq (@seqs) {
          if($seq->{acc} eq $ref_accn) {
            push @selected_for_group, $seq;
            $is_selected{$ref_accn} = 1;
            last;
          }
        }
      }
    }

    # Replace the original group array with the filtered top N
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
  my ($candidate_AHR, $accn_file, $do_keep, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

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

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "tier1.accn", $accn_file, $do_keep, $do_keep, "list of accessions selected for tier-2 download");
  if(! $do_keep) { push(@{$to_remove_AR}, $accn_file); }
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
  my ($accn_file, $fasta_out_file, $do_keep, $ofile_info_HHR, $to_remove_AR, $opt_HHR, $FH_HR) = @_;

  my @acc_A = ();
  open(my $ifh, "<", $accn_file) or die "ERROR: unable to read $accn_file: $!";
  while(my $line = <$ifh>) {
    chomp $line;
    next if($line eq "");
    push(@acc_A, $line);
  }
  close($ifh);

  my $api_key = (opt_IsUsed("--api_key", $opt_HHR)) ? opt_Get("--api_key", $opt_HHR) : undef;

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
    if(defined $api_key) { $url .= "&api_key=$api_key"; }

    my $res = get($url);
    if(! defined $res) {
      die "ERROR: failed to fetch FASTA batch from NCBI efetch (batch starts at accession index $i)";
    }

    print $ofh $res;
    $nfetched_batches++;
    select(undef, undef, undef, 0.35); # polite delay
  }
  close($ofh);

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "tier2.fa", $fasta_out_file, $do_keep, $do_keep, "FASTA sequences downloaded for tier-2 candidates");
  if(! $do_keep) { push(@{$to_remove_AR}, $fasta_out_file); }
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
  my ($candidate_AHR, $fasta_file, $max_ambig_nt, $decision_HR, $ref_accn, $FH_HR) = @_;

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
      # Reference accession is always kept regardless of ambiguity count
      if($ambig_ct_H{$acc} > $max_ambig_nt && $acc ne $ref_accn) {
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
  my ($candidate_AHR, $fasta_file, $annot_outdir, $model_dir, $model_key, $do_keep, $do_verbose, $execs_HR, $decision_HR, $vannot_extra_opts, $FH_HR) = @_;

  my $annot_mkey = $model_key . ".vadr";

  # Note: --out_stk and --keep are incompatible in v-annotate.pl, so we only use --out_stk
  # which outputs the Stockholm alignment we need for Step 8 block extraction
  my $cmd = $execs_HR->{"v-annotate.pl"} . " -f --mdir " . $model_dir . " --mkey " . $annot_mkey . $vannot_extra_opts . " --out_stk " . $fasta_file . " " . $annot_outdir;
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
      # Reference accession is always kept even if v-annotate reports a failure
      if(exists $is_fail_H{$acc} && $acc ne $model_key) {
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
  
  # Return path to alignment file for later use in Step 10 RNA refinement
  # v-annotate creates: <outdir>/<outdir>.vadr.<modelkey>.align.stk
  my $align_stk_file = $annot_outdir . "/" . $annot_outdir_tail . ".vadr." . $model_key . ".align.stk";
  return $align_stk_file;
}

#################################################################
# Subroutine : get_partial_cds_accns_from_ftr()
# Incept     : EPN* Tue Mar 24 2026
#
# Purpose    : Read a v-annotate .vadr.ftr file and return a hash
#              of accessions that have at least one CDS feature with
#              5' or 3' truncation (trc field != 'no').  These
#              sequences have partial CDS coordinates (<N or >N) and
#              should be deprioritised for centroid selection.
#
# Arguments  :
#   $ftr_file : path to .vadr.ftr file
#
# Returns    : hash (accession => 1) of partial-CDS accessions
#################################################################
sub get_partial_cds_accns_from_ftr {
  my ($ftr_file) = @_;
  my %partial_H = ();
  open(my $fh, "<", $ftr_file) || die "ERROR, unable to read ftr file $ftr_file: $!";
  while(my $line = <$fh>) {
    chomp $line;
    next if($line =~ /^\#/);
    next if($line =~ /^\s*$/);
    my @tok_A = split(/\s+/, $line);
    next if(scalar(@tok_A) < 15);
    my $acc  = $tok_A[1];
    my $type = $tok_A[5];
    my $trc  = $tok_A[14];
    next if($type ne "CDS");
    if($trc ne "no") {
      $partial_H{$acc} = 1;
    }
  }
  close($fh);
  return %partial_H;
}

#################################################################
# Subroutine : select_group_centroids_blast()
# Incept     : Copilot Tue Mar 10 2026
#
# Purpose    : For each remaining group, select one centroid
#              sequence using average all-vs-all blastn pident.
#              Sequences with partial CDS coordinates (< or > in
#              coords) are excluded from centroid candidacy unless
#              they are the only sequence in the group.
#################################################################
sub select_group_centroids_blast {
  my ($candidate_AHR, $fasta_file, $out_root, $centroid_tsv_file, $n_per_group, $do_keep, $do_verbose, $execs_HR, $decision_HR, $partial_cds_HR, $ref_accn, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

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

    # Filter out sequences with partial CDS coordinates from centroid candidacy.
    # Only apply when partial info is available and group has non-partial members.
    if(defined $partial_cds_HR) {
      my @non_partial_A = grep { ! $partial_cds_HR->{$_->{"acc"}} } @seq_A;
      if(scalar(@non_partial_A) > 0) {
        @seq_A = @non_partial_A;
      }
      # else: all sequences in group are partial — allow them all as candidates
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

    # Sort accessions by average pident (descending), alphabetical tie-breaker
    my @sorted_accs = sort { $avg_H{$b} <=> $avg_H{$a} || $a cmp $b } @acc_A;

    # Select top N sequences (or all if fewer than N available)
    my $n_to_select = ($n_per_group < scalar(@sorted_accs)) ? $n_per_group : scalar(@sorted_accs);
    my %selected_H = ();
    for(my $i = 0; $i < $n_to_select; $i++) {
      $selected_H{$sorted_accs[$i]} = 1;
    }
    # Also always keep the reference accession
    if(defined $ref_accn) {
      foreach my $seq (@seq_A) {
        if($seq->{"acc"} eq $ref_accn) { $selected_H{$ref_accn} = 1; }
      }
    }

    my $best_avg = $avg_H{$sorted_accs[0]};
    my @new_group_A = ();
    foreach my $seq (@seq_A) {
      my $acc = $seq->{"acc"};
      my $is_selected = exists $selected_H{$acc} ? 1 : 0;
      print $ctfh join("\t", $group, $acc, sprintf("%.4f", $avg_H{$acc}), $is_selected) . "\n";
      if($is_selected) {
        push(@new_group_A, $seq);
        $decision_HR->{$acc}{"status"} = "kept";
        if(defined $ref_accn && $acc eq $ref_accn && ! exists $selected_H{$acc}) {
          $decision_HR->{$acc}{"reason_code"} = "reference_forced";
          $decision_HR->{$acc}{"reason_detail"} = sprintf("reference accession (avg blastn pident %.4f)", $avg_H{$acc});
        }
        else {
          $decision_HR->{$acc}{"reason_code"} = "centroid_selected";
          $decision_HR->{$acc}{"reason_detail"} = sprintf("selected by avg blastn pident %.4f (rank %d of %d in group)",
                                                           $avg_H{$acc}, 1 + (grep { $avg_H{$_} > $avg_H{$acc} } @sorted_accs), scalar(@sorted_accs));
        }
        $decision_HR->{$acc}{"stage_last_seen"} = "selected_for_tier3";
      }
      else {
        $decision_HR->{$acc}{"status"} = "removed";
        $decision_HR->{$acc}{"reason_code"} = "centroid_not_selected";
        $decision_HR->{$acc}{"reason_detail"} = sprintf("not in top %d; avg blastn pident %.4f (best %.4f)", $n_per_group, $avg_H{$acc}, $best_avg);
        $decision_HR->{$acc}{"stage_last_seen"} = "tier3_centroid_filter";
      }
    }
    $candidate_AHR->{$group} = \@new_group_A;
    if(! $do_keep) {
      push(@{$to_remove_AR}, $group_fa);
      push(@{$to_remove_AR}, $group_blast);
      foreach my $ext ("nhr", "nin", "nsq", "nsi", "nsd", "not", "ntf", "nto", "ndb", "nos") {
        my $db_file = $group_db . "." . $ext;
        if(-e $db_file) { push(@{$to_remove_AR}, $db_file); }
      }
    }
    $n_groups_with_centroid++;
  }
  close($ctfh);

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "centroid.tsv", $centroid_tsv_file, 1, 1, "per-group centroid selection table (blastn average pident)");
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Tier 3 centroid selection (blast): selected up to %d seqs in each of %d groups (%d empty groups) and wrote %s\n", $n_per_group, $n_groups_with_centroid, $n_groups_empty, $centroid_tsv_file));
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
  my ($decision_HR, $out_file, $ofile_info_HHR, $FH_HR) = @_;

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

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "filter.seq.tsv", $out_file, 1, 1, "per-sequence filter status/reason table");
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Decision report: wrote per-sequence status/reason table to %s\n", $out_file));
  return;
}

#################################################################
# Subroutine : write_decision_stage_reports()
#################################################################
sub write_decision_stage_reports {
  my ($decision_HR, $out_prefix, $do_keep, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

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
    ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "filter.stage.$stage_safe", $out_file, $do_keep, $do_keep, "per-stage filter status/reason table for stage $stage");
    if(! $do_keep) { push(@{$to_remove_AR}, $out_file); }
    $nfiles++;
  }

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Decision report: wrote %d per-stage decision files with prefix %s.stage.<stage>.tsv\n", $nfiles, $out_prefix));
  return;
}

#################################################################
# Subroutine : write_decision_summary_report()
#################################################################
sub write_decision_summary_report {
  my ($decision_HR, $out_file, $ofile_info_HHR, $FH_HR) = @_;

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

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "filter.sum.tsv", $out_file, 1, 1, "filter status/reason summary counts table");
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Decision summary: wrote grouped status/reason/stage counts to %s\n", $out_file));
  return;
}

#################################################################
# Subroutine : write_stitch_scaffold_outputs()
#################################################################
sub write_stitch_scaffold_outputs {
  my ($candidate_AHR, $tier2_fasta_file, $ftr_info_HAR, $model_key, $seed_model_len, $selected_accn_file, $selected_fa_file, $block_plan_file, $do_keep, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

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

  open(my $bpfh, ">", $block_plan_file) || die "ERROR, unable to write block plan $block_plan_file: $!";
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

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.selected.accn", $selected_accn_file, $do_keep, $do_keep, "list of accessions selected for CDS stitching");
  if(! $do_keep) { push(@{$to_remove_AR}, $selected_accn_file); }
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.selected.fa", $selected_fa_file, $do_keep, $do_keep, "FASTA sequences selected for CDS stitching");
  if(! $do_keep) { push(@{$to_remove_AR}, $selected_fa_file); }
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "block.plan", $block_plan_file, $do_keep, $do_keep, "genome block plan (coding/noncoding partition)");
  if(! $do_keep) { push(@{$to_remove_AR}, $block_plan_file); }
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
  my ($candidate_AHR, $selected_fa_file, $annot_outdir, $cds_nt_fa_file, $orf_fa_file, $aa_fa_file, $map_tsv_file, $do_skip_annotate, $do_keep, $do_verbose, $execs_HR, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

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

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.nt.fa",          $cds_nt_fa_file, $do_keep, $do_keep, "CDS nucleotide sequences for stitching");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.orf.fa",         $orf_fa_file,    $do_keep, $do_keep, "esl-translate ORFs for CDS sequences");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.aa.fa",          $aa_fa_file,     $do_keep, $do_keep, "selected CDS amino acid sequences for stitching");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.translate_map.tsv", $map_tsv_file, $do_keep, $do_keep, "CDS coordinate/ORF mapping table for stitching");
  if(! $do_keep) {
    push(@{$to_remove_AR}, $cds_nt_fa_file);
    push(@{$to_remove_AR}, $orf_fa_file);
    push(@{$to_remove_AR}, $aa_fa_file);
    push(@{$to_remove_AR}, $map_tsv_file);
  }
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS translation prep: wrote %d CDS nt sequences to %s\n", $n_cds, $cds_nt_fa_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS translation prep: wrote esl-translate ORFs to %s\n", $orf_fa_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS translation prep: selected one ORF for %d CDS and wrote %s and %s\n", $n_written, $aa_fa_file, $map_tsv_file));
  return;
}

#################################################################
# Subroutine : run_muscle_aa_alignment()
# Incept:     EPN* Tue Mar 25 2026
#
# Purpose:    Run muscle 3.8.31 multiple alignment on CDS protein
#             sequences, one alignment per CDS feature group. The
#             reference sequence (identified by $ref_accn) is
#             included in each alignment and used to define RF
#             columns in the downstream MSA building step.
#
# Arguments:
#   $aa_fa_file:      CDS amino acid FASTA (from Step 8)
#   $map_tsv_file:    CDS translate_map.tsv (for ftr_idx grouping)
#   $ref_accn:        reference accession (for anchor.tsv)
#   $anchor_tsv_file: output anchor TSV (reference info per feature)
#   $muscle_dir:      output directory for per-feature muscle MSA files
#   $do_skip_annotate: skip flag
#   $do_keep:         keep intermediate files
#   $do_verbose:      verbose flag
#   $execs_HR:        executables hash (must contain "muscle")
#   $ofile_info_HHR:  output file info hash
#   $to_remove_AR:    files to remove
#   $FH_HR:           file handle hash
#
# Returns:    void
#################################################################
sub run_muscle_aa_alignment {
  my ($aa_fa_file, $map_tsv_file, $ref_accn, $anchor_tsv_file, $muscle_dir,
      $do_skip_annotate, $do_keep, $do_verbose, $execs_HR, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

  if($do_skip_annotate) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS muscle AA: skipped due to --skip-annotate\n"));
    return;
  }
  if(! -s $aa_fa_file) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS muscle AA: skipped, missing or empty AA file %s\n", $aa_fa_file));
    return;
  }
  if(! exists $execs_HR->{"muscle"}) {
    die "ERROR, unable to run muscle AA alignment: muscle not found (set VADRMUSCLEDIR so \$VADRMUSCLEDIR/muscle exists and is executable)";
  }

  # Read translate map to get ftr_idx for each AA header
  my %header_to_ftr_idx_H = ();
  if(-s $map_tsv_file) {
    open(my $mapfh, "<", $map_tsv_file) || die "ERROR, unable to read translate map $map_tsv_file: $!";
    my $maphdr = <$mapfh>;  # skip header
    while(my $line = <$mapfh>) {
      chomp $line;
      my @tok_A = split(/\t/, $line);
      my $acc = $tok_A[1];
      my $ftr_idx = $tok_A[2];
      my $seq_coords = $tok_A[3];
      my $model_coords = $tok_A[4];
      my $aa_header = $acc . ":" . $seq_coords . "/" . $model_coords;
      $header_to_ftr_idx_H{$aa_header} = $ftr_idx;
    }
    close($mapfh);
  }

  # Parse AA sequences and group by feature key (ftr_idx)
  my @aa_A = ();
  my $cur_h = undef;
  my $cur_sq = "";
  open(my $aafh_in, "<", $aa_fa_file) || die "ERROR, unable to read AA fasta $aa_fa_file: $!";
  while(my $line = <$aafh_in>) {
    chomp $line;
    if($line =~ /^>(\S+)/) {
      if(defined $cur_h) {
        my $acc = $cur_h;
        if($acc =~ /^([^:]+):/) { $acc = $1; }
        my $fkey = (exists $header_to_ftr_idx_H{$cur_h}) ? $header_to_ftr_idx_H{$cur_h} : (($cur_h =~ /\/(.+)$/) ? $1 : $cur_h);
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
    my $fkey = (exists $header_to_ftr_idx_H{$cur_h}) ? $header_to_ftr_idx_H{$cur_h} : (($cur_h =~ /\/(.+)$/) ? $1 : $cur_h);
    push(@aa_A, { header => $cur_h, accession => $acc, feature_key => $fkey, sqstring => $cur_sq, len => length($cur_sq) });
  }
  if(scalar(@aa_A) == 0) {
    die "ERROR, no AA sequences parsed from $aa_fa_file";
  }

  # Group sequences by feature key
  my %fkey_seqs_A = ();
  my @fkey_order = ();
  foreach my $seq (@aa_A) {
    my $fk = $seq->{"feature_key"};
    if(! exists $fkey_seqs_A{$fk}) {
      push(@fkey_order, $fk);
      $fkey_seqs_A{$fk} = [];
    }
    push(@{$fkey_seqs_A{$fk}}, $seq);
  }

  # Create output directory for per-feature muscle MSA files
  if(! -d $muscle_dir) {
    mkdir($muscle_dir) || die "ERROR, unable to create muscle output directory $muscle_dir: $!";
  }

  # Write anchor TSV with reference as the anchor for each feature
  open(my $atfh, ">", $anchor_tsv_file) || die "ERROR, unable to write anchor TSV $anchor_tsv_file: $!";
  print $atfh join("\t", "anchor_header", "anchor_accession", "anchor_aa_len", "anchor_method", "anchor_centroid_avg_blastn_pident", "feature_key") . "\n";

  my $n_alignments = 0;
  foreach my $fk (@fkey_order) {
    my @grp_A = @{$fkey_seqs_A{$fk}};

    # Find reference sequence in this feature group
    my $ref_seq = undef;
    foreach my $seq (@grp_A) {
      if($seq->{"accession"} eq $ref_accn) {
        $ref_seq = $seq;
        last;
      }
    }
    if(! defined $ref_seq) {
      die "ERROR, reference accession $ref_accn not found in CDS feature group $fk. " .
          "The reference must be in the training set for all CDS features.";
    }

    # Write anchor TSV row (reference = anchor for RF definition)
    print $atfh join("\t", $ref_seq->{"header"}, $ref_seq->{"accession"}, $ref_seq->{"len"}, "reference", "NA", $fk) . "\n";

    # Write per-feature input FASTA for muscle
    my $fk_input_fa  = sprintf("%s/ftr.%s.input.fa", $muscle_dir, $fk);
    my $fk_output_fa = sprintf("%s/ftr.%s.muscle.afa", $muscle_dir, $fk);

    open(my $fkfh, ">", $fk_input_fa) || die "ERROR, unable to write $fk_input_fa: $!";
    foreach my $seq (@grp_A) {
      print $fkfh ">" . $seq->{"header"} . "\n";
      print $fkfh seq_SqstringAddNewlines($seq->{"sqstring"}, 60);
    }
    close($fkfh);

    # Run muscle
    if(scalar(@grp_A) == 1) {
      # Single sequence: no alignment needed, just copy
      open(my $out1, ">", $fk_output_fa) || die "ERROR, unable to write $fk_output_fa: $!";
      print $out1 ">" . $grp_A[0]{"header"} . "\n";
      print $out1 seq_SqstringAddNewlines($grp_A[0]{"sqstring"}, 60);
      close($out1);
    }
    else {
      my $cmd = $execs_HR->{"muscle"} . " -in " . $fk_input_fa . " -out " . $fk_output_fa;
      utl_RunCommand($cmd, $do_verbose, 0, $FH_HR);
    }

    $n_alignments++;
    if(! $do_keep) { push(@{$to_remove_AR}, $fk_input_fa); }
  }
  close($atfh);

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.anchor.tsv", $anchor_tsv_file, 1, 1, "CDS anchor (reference) selection table");
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS muscle AA: %d CDS features, %d muscle alignments\n", scalar(@fkey_order), $n_alignments));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS muscle AA: per-feature MSA files in %s\n", $muscle_dir));
  return;
}

#################################################################
# Subroutine : build_muscle_cds_msa()
# Incept:     EPN* Tue Mar 25 2026
#
# Purpose:    Build concatenated CDS AA and NT MSAs from per-feature
#             muscle alignment output. For each feature, RF is
#             defined from the reference sequence's non-gap positions.
#             Back-translates AA alignment to NT using ORF nucleotides
#             and translation map.
#
# Arguments:
#   $aa_fa_file:       CDS amino acid FASTA (unaligned, for ordering)
#   $cds_nt_fa_file:   CDS nucleotide FASTA (for back-translation)
#   $map_tsv_file:     CDS translate_map.tsv
#   $anchor_tsv_file:  anchor TSV (reference info per feature)
#   $muscle_dir:       directory with per-feature muscle MSA files
#   $msa_aa_fa_file:   output: concatenated AA MSA FASTA
#   $msa_aa_stk_file:  output: concatenated AA MSA Stockholm (with RF)
#   $msa_nt_fa_file:   output: concatenated NT MSA FASTA
#   $do_skip_annotate: skip flag
#   $do_keep:          keep intermediate files
#   $ofile_info_HHR:   output file info hash
#   $to_remove_AR:     files to remove
#   $FH_HR:            file handle hash
#
# Returns:    void
#################################################################
sub build_muscle_cds_msa {
  my ($aa_fa_file, $cds_nt_fa_file, $map_tsv_file, $anchor_tsv_file, $muscle_dir,
      $msa_aa_fa_file, $msa_aa_stk_file, $msa_nt_fa_file,
      $do_skip_annotate, $do_keep, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

  if($do_skip_annotate) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: skipped due to --skip-annotate\n"));
    return;
  }
  if((! -s $aa_fa_file) || (! -s $anchor_tsv_file) || (! -s $cds_nt_fa_file) || (! -s $map_tsv_file)) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: skipped (required input missing or empty)\n"));
    return;
  }

  # Read translate map to get ftr_idx for each AA header
  my %header_to_ftr_idx_H = ();
  if(-s $map_tsv_file) {
    open(my $mapfh, "<", $map_tsv_file) || die "ERROR, unable to read translate map $map_tsv_file: $!";
    my $maphdr = <$mapfh>;  # skip header
    while(my $line = <$mapfh>) {
      chomp $line;
      my @tok_A = split(/\t/, $line);
      my $acc = $tok_A[1];
      my $ftr_idx = $tok_A[2];
      my $seq_coords = $tok_A[3];
      my $model_coords = $tok_A[4];
      my $aa_header = $acc . ":" . $seq_coords . "/" . $model_coords;
      $header_to_ftr_idx_H{$aa_header} = $ftr_idx;
    }
    close($mapfh);
  }

  # Parse original AA sequences (unaligned, for accession ordering)
  my @aa_order_A = ();
  my %aa_fkey_H = ();
  {
    my $cur_h = undef;
    open(my $aafh, "<", $aa_fa_file) || die "ERROR, unable to read AA fasta $aa_fa_file: $!";
    while(my $line = <$aafh>) {
      chomp $line;
      if($line =~ /^>(\S+)/) {
        $cur_h = $1;
        if(! exists $aa_fkey_H{$cur_h}) {
          push(@aa_order_A, $cur_h);
          $aa_fkey_H{$cur_h} = (exists $header_to_ftr_idx_H{$cur_h}) ? $header_to_ftr_idx_H{$cur_h} : (($cur_h =~ /\/(.+)$/) ? $1 : $cur_h);
        }
      }
    }
    close($aafh);
  }

  # Read per-feature anchor (reference) info from anchor TSV
  my %anchor_header_H = ();
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

  # Build accession order from first feature group
  my @accn_order_A = ();
  my %accn_seen_H = ();
  foreach my $h (@aa_order_A) {
    my $fk = $aa_fkey_H{$h};
    if($fk eq $fkey_order[0]) {
      my $acc = ($h =~ /^([^:]+):/) ? $1 : $h;
      if(! exists $accn_seen_H{$acc}) {
        push(@accn_order_A, $acc);
        $accn_seen_H{$acc} = 1;
      }
    }
  }

  # Read per-feature muscle MSA files, build concatenated AA MSA and RF
  my %concat_msa_aa_H = (); # accession => concatenated AA alignment string
  my $concat_rf = "";
  foreach my $acc (@accn_order_A) { $concat_msa_aa_H{$acc} = ""; }

  my @fk_aa_msa_len_A = (); # per-feature MSA width (for NT back-translation offset)
  foreach my $fk (@fkey_order) {
    my $fk_msa_file = sprintf("%s/ftr.%s.muscle.afa", $muscle_dir, $fk);
    if(! -s $fk_msa_file) {
      die "ERROR, per-feature muscle MSA file not found: $fk_msa_file";
    }

    # Read muscle MSA
    my %msa_seq_H = (); # header => aligned sequence
    my @msa_order_A = ();
    my $cur_h = undef;
    open(my $msafh, "<", $fk_msa_file) || die "ERROR, unable to read $fk_msa_file: $!";
    while(my $line = <$msafh>) {
      chomp $line;
      if($line =~ /^>(\S+)/) {
        $cur_h = $1;
        if(! exists $msa_seq_H{$cur_h}) {
          push(@msa_order_A, $cur_h);
          $msa_seq_H{$cur_h} = "";
        }
      }
      elsif(defined $cur_h) {
        $line =~ s/\s+//g;
        $msa_seq_H{$cur_h} .= $line;
      }
    }
    close($msafh);

    my $fk_msa_len = length($msa_seq_H{$msa_order_A[0]});
    push(@fk_aa_msa_len_A, $fk_msa_len);

    # Find reference sequence and build RF from its non-gap positions
    my $ref_header = $anchor_header_H{$fk};
    if(! exists $msa_seq_H{$ref_header}) {
      die "ERROR, reference header $ref_header not found in muscle MSA $fk_msa_file";
    }
    my $ref_aln = $msa_seq_H{$ref_header};
    my $fk_rf = "";
    for(my $i = 0; $i < length($ref_aln); $i++) {
      my $c = substr($ref_aln, $i, 1);
      $fk_rf .= ($c eq '-' || $c eq '.') ? '.' : 'x';
    }
    $concat_rf .= $fk_rf;

    # Build per-accession aligned AA strings from muscle MSA
    # Map header -> accession
    my %header_to_acc_H = ();
    foreach my $h (@msa_order_A) {
      my $acc = ($h =~ /^([^:]+):/) ? $1 : $h;
      $header_to_acc_H{$h} = $acc;
    }

    foreach my $acc (@accn_order_A) {
      # Find this accession's header in the MSA
      my $found_h = undef;
      foreach my $h (@msa_order_A) {
        if($header_to_acc_H{$h} eq $acc) {
          $found_h = $h;
          last;
        }
      }
      if(defined $found_h) {
        $concat_msa_aa_H{$acc} .= $msa_seq_H{$found_h};
      }
      else {
        # Accession not in this feature (e.g. partial CDS) — fill with gaps
        $concat_msa_aa_H{$acc} .= ('-' x $fk_msa_len);
      }
    }
  }

  # Write per-accession concatenated AA MSA FASTA
  open(my $maa_fh, ">", $msa_aa_fa_file) || die "ERROR, unable to write $msa_aa_fa_file: $!";
  foreach my $acc (@accn_order_A) {
    print $maa_fh ">" . $acc . "\n";
    print $maa_fh seq_SqstringAddNewlines($concat_msa_aa_H{$acc}, 60);
  }
  close($maa_fh);

  # Write AA MSA Stockholm with RF
  open(my $mstk_fh, ">", $msa_aa_stk_file) || die "ERROR, unable to write $msa_aa_stk_file: $!";
  print $mstk_fh "# STOCKHOLM 1.0\n";
  foreach my $acc (@accn_order_A) {
    print $mstk_fh $acc . "\t" . $concat_msa_aa_H{$acc} . "\n";
  }
  print $mstk_fh "#=GC RF\t" . $concat_rf . "\n";
  print $mstk_fh "//\n";
  close($mstk_fh);

  # ---------------------------------------------------------------
  # NT back-translation (same logic as original, reads from muscle MSA)
  # ---------------------------------------------------------------

  # Read CDS nt sequences
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

  # Read translation map (also extract strand per feature)
  my %map_H = ();
  my %fkey_strand_H = (); # feature_key => strand ("+" or "-")
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
    # Extract strand from seq_coords (e.g. "1734..1033:-" or "51..995:+")
    if($seq_coords =~ /:([\+\-])$/) {
      $fkey_strand_H{$ftr_idx} = $1 if(! exists $fkey_strand_H{$ftr_idx});
    }
  }
  close($mapfh);

  # Group headers by feature key
  my %fkey_headers_A = ();
  foreach my $h (@aa_order_A) {
    my $fk = $aa_fkey_H{$h};
    if(! exists $fkey_headers_A{$fk}) { $fkey_headers_A{$fk} = []; }
    push(@{$fkey_headers_A{$fk}}, $h);
  }

  # NT backconversion: per-feature, then concatenate per-accession
  my %concat_nt_H = ();
  foreach my $acc (@accn_order_A) { $concat_nt_H{$acc} = ""; }

  foreach my $fk_idx (0..$#fkey_order) {
    my $fk = $fkey_order[$fk_idx];
    my @fk_headers = exists $fkey_headers_A{$fk} ? @{$fkey_headers_A{$fk}} : ();
    my $fk_msa_len = $fk_aa_msa_len_A[$fk_idx];
    my $fk_aa_offset = 0;
    for(my $i = 0; $i < $fk_idx; $i++) { $fk_aa_offset += $fk_aa_msa_len_A[$i]; }

    # Compute per-feature max_n5 and max_n3
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
      my $fk_nt_aln = $prefix_pad . $prefix_nt . $nt_aln . $suffix_nt . $suffix_pad;

      # For minus-strand CDS features, reverse-complement the NT alignment
      # to convert from mRNA orientation back to forward genomic orientation.
      # Gaps ('-') are preserved; only residues are complemented; column order
      # is reversed so the alignment stays consistent across all sequences.
      my $fk_strand = (exists $fkey_strand_H{$fk}) ? $fkey_strand_H{$fk} : "+";
      if($fk_strand eq "-") {
        seq_SqstringReverseComplement(\$fk_nt_aln);
      }

      $concat_nt_H{$acc} .= $fk_nt_aln;
    }
  }

  # Write concatenated NT MSA FASTA
  open(my $mnt_fh, ">", $msa_nt_fa_file) || die "ERROR, unable to write $msa_nt_fa_file: $!";
  my $n_nt = 0;
  foreach my $acc (@accn_order_A) {
    my @acc_headers = grep { /^\Q$acc\E:/ } @aa_order_A;
    my $nt_header = join("+", @acc_headers);
    print $mnt_fh ">" . $nt_header . "\n";
    print $mnt_fh seq_SqstringAddNewlines($concat_nt_H{$acc}, 60);
    $n_nt++;
  }
  close($mnt_fh);

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.msa.aa.afa", $msa_aa_fa_file,  $do_keep, $do_keep, "concatenated CDS protein MSA (aligned FASTA)");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.msa.aa.stk", $msa_aa_stk_file, 1,        1,        "concatenated CDS protein MSA (Stockholm with RF)");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "cds.msa.nt.afa", $msa_nt_fa_file,  1,        1,        "concatenated CDS nucleotide MSA (aligned FASTA)");
  if(! $do_keep) { push(@{$to_remove_AR}, $msa_aa_fa_file); }
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: %d CDS features, wrote protein MSA to %s\n", scalar(@fkey_order), $msa_aa_fa_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: wrote protein MSA Stockholm (with RF) to %s\n", $msa_aa_stk_file));
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Stitch CDS MSA/backconvert: wrote CDS nucleotide MSA for %d sequences to %s\n", $n_nt, $msa_nt_fa_file));
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
  my ($ref_seq_file, $rna_cm_file, $seq_len, $out_root, $rfam_dir, $rna_regions_AR, $ss_cons_SR, $do_keep, $do_verbose, $execs_HR, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

  my $cmscan_tblout = $out_root . ".rna_cmscan.tblout";
  my $cmscan_stdout = $out_root . ".rna_cmscan.out";
  my $cmalign_tfile = $out_root . ".rna_cmalign.ifile";
  my $cmalign_stk   = $out_root . ".rna_cmalign.stk";

  # Step 3.1: Run cmscan to identify RNA hits
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

  # Step 3.2: Parse cmscan tblout to extract RNA hits
  parse_cmscan_tblout($cmscan_tblout, $seq_len, $rna_regions_AR, $FH_HR);

  # Step 3.3: Run cmalign --tfile to get full-sequence consensus secondary structure
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA discovery: running cmalign to generate consensus secondary structure\n"));
  
  # For cmalign, we need to use a single CM from the reference model
  # TODO: Determine which CM to use for full-sequence alignment (may need seed model CM)
  # For now, skip this step - will implement after testing cmscan parsing
  
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "rna.cmscan.tblout", $cmscan_tblout, 1,        1,        "cmscan tblout for RNA discovery");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "rna.cmscan.out",   $cmscan_stdout, $do_keep, $do_keep, "cmscan stdout for RNA discovery");
  if(! $do_keep) {
    push(@{$to_remove_AR}, $cmscan_stdout);
    push(@{$to_remove_AR}, $cmalign_tfile) if(-e $cmalign_tfile);
    push(@{$to_remove_AR}, $cmalign_stk)   if(-e $cmalign_stk);
  }
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA discovery: found %d RNA regions\n", scalar(@{$rna_regions_AR})));

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
  my ($rna_regions_AR, $out_file, $ofile_info_HHR, $FH_HR) = @_;

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

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "rna.annot.tsv", $out_file, 1, 1, "RNA annotation table");
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
  my ($rna_regions_AR, $ref_seq_file, $rna_cm_file, $rfam_dir, $out_root, $do_keep, $do_verbose, $execs_HR, $ofile_info_HHR, $FH_HR) = @_;

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
#              and realign with custom CMs built in Step 3b.
#              This implements Step 11 from the updated plan.
#              Uses Bio::Easel::MSA column_subset() for extraction
#              and write_single_unaligned_seq() for FASTA output.
#
# Arguments  :
#   $rna_regions_AR     : ref to array of RNA region hashes from Step 3
#   $tier2_align_stk    : path to tier2 v-annotate Stockholm alignment 
#   $rna_struct_dir     : directory with CM files from Step 3b (rna.001.cm, etc)
#   $out_root           : output file root path
#   $do_keep            : keep intermediate files
#   $do_verbose         : verbose output
#   $execs_HR           : hash of executable paths
#   $FH_HR              : hash of file handles
#
# Returns    : void (writes RNA block Stockholm files)
#################################################################
sub extract_and_align_rna_regions {
  my ($rna_regions_AR, $tier2_align_stk, $rna_struct_dir, $out_root, $do_keep, $do_verbose, $execs_HR, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

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

    # Step 11: Extract RNA region columns from tier2 alignment via Bio::Easel
    # Map RF positions to alignment columns and extract with column_subset
    my $rna_extracted_fa = $out_root . ".rna." . sprintf("%03d", $idx) . ".extracted.fa";

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

    # Step 11: Align with custom CM from Step 3b
    my $cm_file = $rna_struct_dir . "/rna." . sprintf("%03d", $idx) . ".cm";
    if(! -e $cm_file) {
      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: WARNING - CM not found for region %d, skipping\n", $idx));
      $idx++;
      next;
    }

    # Step 11: Output refined RNA block Stockholm
    my $rna_aligned_stk = $out_root . ".rna." . sprintf("%03d", $idx) . ".stk";

    my $cmd_align = $execs_HR->{"cmalign"} . " --outformat pfam -g " . $cm_file . " " . $rna_extracted_fa .
                    " > " . $rna_aligned_stk;

    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: region %d (%s): aligning with custom CM\n",
                                                     $idx, $rna_family));
    utl_RunCommand($cmd_align, $do_verbose, 0, $FH_HR);

    ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "rna." . sprintf("%03d", $idx) . ".stk",          $rna_aligned_stk,  1,        1,        sprintf("RNA region %d (%s) aligned Stockholm", $idx, $rna_family));
    ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "rna." . sprintf("%03d", $idx) . ".extracted.fa", $rna_extracted_fa, $do_keep, $do_keep, sprintf("RNA region %d (%s) extracted sequences", $idx, $rna_family));
    if(! $do_keep) { push(@{$to_remove_AR}, $rna_extracted_fa); }
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# RNA alignment: region %d (%s): wrote aligned Stockholm to %s\n",
                                                     $idx, $rna_family, $rna_aligned_stk));

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
      $rna_regions_AR, $final_stk_file, $temp_cm_file,
      $do_rna_discovery, $do_skip_annotate, $seed_model_len, $execs_HR,
      $do_keep, $ofile_info_HHR, $to_remove_AR, $FH_HR) = @_;

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
        "stk_file" => $out_root . ".rna." . sprintf("%03d", $rna_idx) . ".stk"
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
  my $anchor_tsv_file = $out_root . ".cds.anchor.tsv";
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
  my $cmbuild_out = $final_stk_file . ".cmbuild.out";
  my $output_stk_file = $out_root . ".stk";
  my $cmd_cmbuild = $execs_HR->{"cmbuild"} . " --hand -O " . $output_stk_file . " " .
                    $temp_cm_file . " " . $final_stk_file . " > " . $cmbuild_out;

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: running cmbuild to add RF annotation (output to %s)\n", $cmbuild_out));
  utl_RunCommand($cmd_cmbuild, 1, 0, $FH_HR);

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "final.stk",         $final_stk_file,   $do_keep, $do_keep, "pre-refinement stitched alignment (input to cmbuild)");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "final.cmbuild.out", $cmbuild_out,      $do_keep, $do_keep, "cmbuild output for final alignment refinement");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "out.stk",            $output_stk_file,  1,        1,        "final RF-annotated training alignment");
  if(! $do_keep) {
    push(@{$to_remove_AR}, $final_stk_file);
    push(@{$to_remove_AR}, $cmbuild_out);
    push(@{$to_remove_AR}, $temp_cm_file) if(-e $temp_cm_file);
  }
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: wrote RF-annotated alignment to %s\n", $output_stk_file));

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

  # Collect CDS coordinate ranges for clipping RNA blocks that straddle a CDS boundary.
  # When an RNA region only partially overlaps CDS, it gets its own block for the non-CDS
  # portion.  The CDS block handles the overlapping portion via the ungapped_ss_cons overlay
  # (build_ungapped_ss_cons overlays the full RNA structure, so the CDS block already carries
  # the RNA SS_cons characters for positions within the CDS).  We clip the RNA block here so
  # it only covers the non-CDS portion; base pairs that span the boundary will have one half
  # in the CDS block and the other in the RNA block, keeping them matched in the final SS_cons.
  my @cds_ranges_A = ();
  foreach my $block (@{$blocks_AR}) {
    if($block->{"type"} eq "coding") {
      push(@cds_ranges_A, { start => $block->{"start"}, end => $block->{"end"} });
    }
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

      # Clip RNA block to non-CDS portion if it straddles a CDS boundary.
      # The CDS-overlapping portion is already handled by build_ungapped_ss_cons overlaying
      # the full RNA structure into the CDS block's SS_cons.
      my $eff_rna_start = $rna_start;
      my $eff_rna_end   = $rna_end;
      foreach my $cds (@cds_ranges_A) {
        # RNA starts in noncoding and extends into CDS on the right
        if($eff_rna_start < $cds->{"start"} && $eff_rna_end >= $cds->{"start"}) {
          $eff_rna_end = $cds->{"start"} - 1;
        }
        # RNA starts inside CDS and extends into noncoding on the right
        if($eff_rna_start <= $cds->{"end"} && $eff_rna_end > $cds->{"end"}) {
          $eff_rna_start = $cds->{"end"} + 1;
        }
      }
      # 1-based RF position range within the RNA stk to extract for the clipped block
      my $rna_clip_rf_start = $eff_rna_start - $rna_start + 1;
      my $rna_clip_rf_end   = $eff_rna_end   - $rna_start + 1;
      my $is_clipped = ($eff_rna_start != $rna_start || $eff_rna_end != $rna_end);

      # Add segment before RNA if exists
      if($cur_pos < $eff_rna_start) {
        push(@merged_A, {
          "type"     => $block_type,
          "start"    => $cur_pos,
          "end"      => $eff_rna_start - 1,
          "len"      => $eff_rna_start - $cur_pos,
          "source"   => $block->{"source"}
        });
      }

      # Add RNA block (clipped to non-CDS portion)
      my %rna_block = (
        "type"     => "rna",
        "start"    => $eff_rna_start,
        "end"      => $eff_rna_end,
        "len"      => $eff_rna_end - $eff_rna_start + 1,
        "source"   => "rna_discovery",
        "rna_idx"  => $rna->{"idx"},
        "rna_file" => $rna->{"stk_file"},
        "family"   => $rna->{"family"}
      );
      if($is_clipped) {
        $rna_block{"rna_clip_rf_start"} = $rna_clip_rf_start;
        $rna_block{"rna_clip_rf_end"}   = $rna_clip_rf_end;
        ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: clipped RNA block %s to non-CDS portion (%d..%d, RF pos %d..%d of %d)\n",
                                                        $rna->{"family"}, $eff_rna_start, $eff_rna_end,
                                                        $rna_clip_rf_start, $rna_clip_rf_end,
                                                        $rna_end - $rna_start + 1));
      }
      push(@merged_A, \%rna_block);

      $cur_pos = $eff_rna_end + 1;
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

  # Pre-compute CDS MSA column slice [start, end] (0-based) for each coding block,
  # using the anchor sequence's non-gap positions as a ruler.
  # This ensures each coding block writes only its own columns in genomic order,
  # so that RF position n in the stitched alignment corresponds to genome position n.
  my $anchor_cds_seq = $accn_to_cds_seq_H{$anchor_accn};
  if(!defined $anchor_cds_seq) {
    die "ERROR in concatenate_all_blocks: anchor accession $anchor_accn not found in CDS MSA $cds_msa_fa_file";
  }
  my @cds_col_ranges_A = ();  # array of [start_col, end_col] (0-based) per coding block
  {
    my $col_cursor = 0;
    my $total_anchor_len = length($anchor_cds_seq);
    foreach my $blk (@{$blocks_AR}) {
      next unless $blk->{"type"} eq "coding";
      my $need = $blk->{"end"} - $blk->{"start"} + 1;  # genomic nt for this block
      my $start_col = $col_cursor;
      my $consumed = 0;
      while($col_cursor < $total_anchor_len && $consumed < $need) {
        my $c = substr($anchor_cds_seq, $col_cursor, 1);
        if($c ne '-' && $c ne '.') { $consumed++; }
        $col_cursor++;
      }
      push @cds_col_ranges_A, [$start_col, $col_cursor - 1];
    }
  }
  my $cds_blk_idx = 0;  # index into @cds_col_ranges_A, incremented per coding block

  # Open output file and write header
  open(my $outfh, ">", $out_stk_file) || die "ERROR unable to write $out_stk_file: $!";
  print $outfh "# STOCKHOLM 1.0\n";

  foreach my $block (@{$blocks_AR}) {
    my $type        = $block->{"type"};
    my $block_start = $block->{"start"};
    my $block_end   = $block->{"end"};

    if($type eq "coding") {
      # Extract only this block's slice of the CDS MSA
      my $slice_start = $cds_col_ranges_A[$cds_blk_idx][0];
      my $slice_end   = $cds_col_ranges_A[$cds_blk_idx][1];
      $cds_blk_idx++;
      my $slice_len = $slice_end - $slice_start + 1;

      my $anchor_slice = substr($anchor_cds_seq, $slice_start, $slice_len);

      # Build RF from anchor slice: non-gap -> x, gap -> .
      my $cds_rf = "";
      foreach my $char (split(//, $anchor_slice)) {
        $cds_rf .= ($char eq '-' || $char eq '.') ? '.' : 'x';
      }

      # Write sequence lines in canonical order (slice of each sequence)
      my $full_anchor_len = length($anchor_cds_seq);
      foreach my $name (@seq_names) {
        my $full_seq = exists $accn_to_cds_seq_H{$name} ? $accn_to_cds_seq_H{$name} : '-' x $full_anchor_len;
        my $seq = substr($full_seq, $slice_start, $slice_len);
        printf $outfh "%-30s %s\n", $name, $seq;
      }
      printf $outfh "#=GC %-24s %s\n", "RF", $cds_rf;

      # Build gapped SS_cons for this coding block
      my $ungapped_substr = substr($ungapped_ss_cons, $block_start - 1, $block_end - $block_start + 1);
      my $cds_ss = "";
      my $ss_pos = 0;
      foreach my $char (split(//, $anchor_slice)) {
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

      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Final stitching: added CDS MSA slice (%d..%d, %d alignment columns)\n",
                                                       $block_start, $block_end, $slice_len));
    }
    elsif($type eq "rna") {
      # Read RNA Stockholm via Bio::Easel::MSA
      my $rna_stk = $block->{"rna_file"};
      my $rna_msa = Bio::Easel::MSA->new({ fileLocation => $rna_stk, isDna => 1 });

      # If this RNA block was clipped to the non-CDS portion, extract only those columns.
      # The CDS-overlapping portion is already in the CDS block's SS_cons via the
      # ungapped_ss_cons overlay, so we only need the non-CDS columns here.
      # IMPORTANT: save SS_cons BEFORE column_subset — Easel removes unmatched bracket
      # symbols during column_subset, but we need to preserve them because their mates
      # reside in the CDS block (via ungapped_ss_cons).  Using the raw pre-clip string
      # ensures the final concatenated SS_cons has all base pairs matched.
      my $clipped_ss_cons = undef;
      if(defined $block->{"rna_clip_rf_start"}) {
        my $clip_rf_start = $block->{"rna_clip_rf_start"};  # 1-based first RF pos to keep
        my $clip_rf_end   = $block->{"rna_clip_rf_end"};    # 1-based last RF pos to keep
        my $full_rf = $rna_msa->has_rf() ? $rna_msa->get_rf() : 'x' x $rna_msa->alen();
        # Capture original SS_cons before column_subset can modify it
        my $orig_ss = $rna_msa->has_ss_cons() ? $rna_msa->get_ss_cons() : '.' x $rna_msa->alen();
        my @useme_A = (0) x $rna_msa->alen();
        my $rf_pos = 0;
        my $first_col = -1;
        my $last_col  = -1;
        for(my $a = 0; $a < $rna_msa->alen(); $a++) {
          my $is_consensus = (substr($full_rf, $a, 1) ne '.' && substr($full_rf, $a, 1) ne '-');
          if($is_consensus) { $rf_pos++; }
          if($is_consensus && $rf_pos == $clip_rf_start) { $first_col = $a; }
          if($is_consensus && $rf_pos == $clip_rf_end)   { $last_col  = $a; }
        }
        if($first_col >= 0 && $last_col >= $first_col) {
          # Save the SS_cons characters for these columns before Easel touches them
          $clipped_ss_cons = substr($orig_ss, $first_col, $last_col - $first_col + 1);
          for(my $a = $first_col; $a <= $last_col; $a++) { $useme_A[$a] = 1; }
          $rna_msa->column_subset(\@useme_A);
        }
      }

      my $rna_aln_width = $rna_msa->alen();
      my $rna_rf      = $rna_msa->has_rf() ? $rna_msa->get_rf() : 'x' x $rna_aln_width;
      # Use saved pre-clip SS_cons if available; otherwise use MSA's SS_cons
      my $rna_ss_cons = defined($clipped_ss_cons) ? $clipped_ss_cons :
                        ($rna_msa->has_ss_cons()  ? $rna_msa->get_ss_cons() : '.' x $rna_aln_width);

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

      # If column_subset produced a zero-width alignment (e.g., trailing
      # noncoding region beyond the tier2 alignment), use the expected
      # block width and fill with gaps
      if($nc_aln_width == 0) { $nc_aln_width = $block_end - $block_start + 1; }

      # Skip blocks with zero width entirely
      if($nc_aln_width == 0) { next; }

      # Write sequence lines in canonical order
      foreach my $name (@seq_names) {
        my $seq = (exists $nc_seqs_H{$name} && length($nc_seqs_H{$name}) > 0)
                  ? $nc_seqs_H{$name} : '-' x $nc_aln_width;
        printf $outfh "%-30s %s\n", $name, $seq;
      }
      # Ensure RF and SS_cons match the actual sequence width
      if(length($nc_rf) == 0 || length($nc_rf) != $nc_aln_width) {
        $nc_rf = 'x' x $nc_aln_width;
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
# Subroutine : annotate_stk_group_subgroup()
# Incept     : EPN* Wed Mar 25 2026
# Purpose    : Add #=GS GP and #=GS SG annotation lines to a
#              Stockholm alignment file.  GP is the user-supplied
#              or auto-guessed virus group name (same for all seqs).
#              SG is the per-sequence subgroup (group_key) from the
#              centroid TSV file (only is_centroid==1 rows are used).
#              Lines are inserted immediately after the last existing
#              #=GS line (typically the WT weight lines from cmbuild).
#
# Arguments  :
#   $stk_file       : path to Stockholm file to annotate (overwritten in place)
#   $centroid_tsv   : path to centroid.tsv (cols: group_key, accession, avg_pident, is_centroid)
#   $group_name     : GP value (same for all sequences)
#   $FH_HR          : ref to hash of file handles
#
# Returns    : void
#################################################################
sub annotate_stk_group_subgroup {
  my ($stk_file, $centroid_tsv, $group_name, $FH_HR) = @_;

  # Build accession -> subgroup map from centroid TSV (is_centroid==1 rows only)
  my %sg_H = ();
  if(-s $centroid_tsv) {
    open(my $cfh, "<", $centroid_tsv) || die "ERROR unable to read centroid TSV $centroid_tsv: $!";
    my $hdr = <$cfh>;  # skip header
    while(my $line = <$cfh>) {
      chomp $line;
      my @f = split(/\t/, $line);
      my ($group_key, $accn, $avg, $is_centroid) = @f;
      if(defined $is_centroid && $is_centroid eq "1") {
        $sg_H{$accn} = $group_key;
      }
    }
    close($cfh);
  }

  # Read Stockholm file
  open(my $infh, "<", $stk_file) || die "ERROR unable to read STK file $stk_file: $!";
  my @lines = <$infh>;
  close($infh);

  # Collect sequence names in alignment-order (first occurrence only)
  my @seqnames = ();
  my %seen_H = ();
  for my $line (@lines) {
    next if $line =~ /^#/;
    next if $line =~ /^\/\//;
    next if $line =~ /^\s*$/;
    if($line =~ /^(\S+)\s/) {
      my $seqname = $1;
      if(! $seen_H{$seqname}) {
        push(@seqnames, $seqname);
        $seen_H{$seqname} = 1;
      }
    }
  }

  # Find index of last #=GS line for insertion point
  my $last_gs_idx = -1;
  for(my $i = 0; $i < scalar(@lines); $i++) {
    if($lines[$i] =~ /^#=GS/) { $last_gs_idx = $i; }
  }
  # Fallback: insert after last #=GF line if no #=GS lines exist
  if($last_gs_idx == -1) {
    for(my $i = 0; $i < scalar(@lines); $i++) {
      if($lines[$i] =~ /^#=GF/) { $last_gs_idx = $i; }
    }
  }

  # Rewrite STK file with GP/SG lines inserted after insertion point
  open(my $outfh, ">", $stk_file) || die "ERROR unable to write STK file $stk_file: $!";
  for(my $i = 0; $i < scalar(@lines); $i++) {
    print $outfh $lines[$i];
    if($i == $last_gs_idx) {
      for my $seqname (@seqnames) {
        my $sg = (exists $sg_H{$seqname}) ? $sg_H{$seqname} : "unknown";
        print $outfh "#=GS $seqname GP $group_name\n";
        print $outfh "#=GS $seqname SG $sg\n";
      }
    }
  }
  close($outfh);

  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# Added GP/SG annotations for %d sequences to %s (GP: %s)\n",
                                                   scalar(@seqnames), $stk_file, $group_name));
  return;
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
  my ($seed_minfo_file, $rna_annot_file, $rna_regions_AR, $out_minfo_file, $model_key, $execs_HR, $ofile_info_HHR, $FH_HR) = @_;

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

  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "updated.minfo", $out_minfo_file, 1, 1, "updated model info file with RNA features");
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

#################################################################
# Subroutine: parse_alt_for_cds_boundary_alerts()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Parse a .vadr.alt file from v-annotate.pl output
#             and extract CDS boundary alerts (mutendex, cdsstopn,
#             cdsstopp, mutendcd, mutstart) grouped by feature and
#             alternative model coordinate.
#
# Arguments:
#   $alt_file: path to .vadr.alt file
#   $FH_HR:    REF to hash of file handles
#
# Returns: REF to hash of hash of arrays:
#   key1: "<ftr_name>:<ftr_idx>"
#   key2: "<mdl_coords>"
#   value: array ref of hashrefs with keys:
#          acc, alert_code, seq_coords, mdl_coords, fail, detail
#################################################################
sub parse_alt_for_cds_boundary_alerts {
  my $sub_name = "parse_alt_for_cds_boundary_alerts";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($alt_file, $FH_HR) = @_;

  my %target_alerts = map { $_ => 1 } qw(mutendex cdsstopn cdsstopp mutendcd mutstart);
  my %alt_groups = ();

  open(my $fh, "<", $alt_file) || ofile_FAIL("ERROR in $sub_name, unable to open $alt_file for reading", 1, $FH_HR);
  while(my $line = <$fh>) {
    chomp $line;
    next if($line =~ /^\#/);
    next if($line =~ /^\s*$/);

    my @tok = split(/\s+/, $line);
    next if(scalar(@tok) < 13);

    my $seq_name    = $tok[1];
    my $ftr_type    = $tok[3];
    my $ftr_name    = $tok[4];
    my $ftr_idx     = $tok[5] - 1;  # convert from 1-based (.vadr.alt) to 0-based (internal)
    my $alert_code  = $tok[6];
    my $fail        = $tok[7];
    my $seq_coords  = $tok[9];
    my $mdl_coords  = $tok[11];
    my $detail      = join(" ", @tok[13..$#tok]);

    next if(! exists $target_alerts{$alert_code});
    next if($ftr_type ne "CDS");
    next if($mdl_coords eq "-");

    my $ftr_key = $ftr_name . ":" . $ftr_idx;

    if(! exists $alt_groups{$ftr_key}{$mdl_coords}) {
      $alt_groups{$ftr_key}{$mdl_coords} = [];
    }
    push(@{$alt_groups{$ftr_key}{$mdl_coords}},
         { acc        => $seq_name,
           alert_code => $alert_code,
           seq_coords => $seq_coords,
           mdl_coords => $mdl_coords,
           fail       => $fail,
           detail     => $detail });
  }
  close($fh);

  return \%alt_groups;
}

#################################################################
# Subroutine: count_independent_observations()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Given an array of accession strings, count the number
#             of independent submissions by grouping accessions by
#             their alphabetic prefix (stripping version and trailing
#             digits). Accessions from the same submission tend to
#             have sequential numbers with the same prefix.
#
# Arguments:
#   $accessions_AR: REF to array of accession strings
#
# Returns: integer count of unique submission prefixes
#################################################################
sub count_independent_observations {
  my $sub_name = "count_independent_observations";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($accessions_AR) = @_;

  my %prefix_H = ();
  foreach my $acc (@{$accessions_AR}) {
    # Strip version suffix (e.g., ".1")
    my $bare = $acc;
    $bare =~ s/\.\d+$//;
    # Extract alphabetic prefix (e.g., "DQ" from "DQ649478", "K" from "K02990")
    my ($prefix) = ($bare =~ /^([A-Za-z]+)/);
    $prefix = "unknown" if(! defined $prefix);
    $prefix_H{$prefix} = 1;
  }

  return scalar(keys %prefix_H);
}

#################################################################
# Subroutine: detect_alternative_features()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Given grouped CDS boundary alerts from
#             parse_alt_for_cds_boundary_alerts(), apply the
#             independence threshold and determine which alternative
#             CDS features to add.
#
# Arguments:
#   $alt_groups_HHR:  REF to hash of hash of arrays from
#                     parse_alt_for_cds_boundary_alerts()
#   $min_independent: minimum number of independent observations
#   $max_fract_diff:  max fractional length deviation allowed
#                     (e.g., 0.2 means alt CDS must be within
#                     80%-120% of original CDS length)
#   $ftr_info_AHR:    REF to array of feature info hashes
#   $model_key:       model name/key
#   $FH_HR:           REF to hash of file handles
#
# Returns: REF to array of hashrefs, each with keys:
#          ftr_name, ftr_idx, alert_code, mdl_alt_coords,
#          original_coords, new_coords, n_seqs, n_independent,
#          accessions_AR
#################################################################
sub detect_alternative_features {
  my $sub_name = "detect_alternative_features";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($alt_groups_HHR, $min_independent, $max_fract_diff, $ftr_info_AHR, $model_key, $FH_HR) = @_;

  my @alt_features = ();

  foreach my $ftr_key (sort keys %{$alt_groups_HHR}) {
    my ($ftr_name, $ftr_idx) = split(/:/, $ftr_key);
    $ftr_idx = int($ftr_idx);

    # Get original CDS coords from feature info
    my $original_coords = $ftr_info_AHR->[$ftr_idx]{"coords"};
    if(! defined $original_coords) {
      ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# WARNING: $sub_name, unable to find coords for feature $ftr_key, skipping\n"));
      next;
    }

    foreach my $mdl_coords (sort keys %{$alt_groups_HHR->{$ftr_key}}) {
      my @entries = @{$alt_groups_HHR->{$ftr_key}{$mdl_coords}};
      my @accessions = map { $_->{"acc"} } @entries;
      my $n_ind = count_independent_observations(\@accessions);

      if($n_ind < $min_independent) {
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# Alt detect: %s %s n_seqs=%d n_ind=%d < min=%d, skipping\n",
                  $ftr_key, $mdl_coords, scalar(@entries), $n_ind, $min_independent));
        next;
      }

      # Determine alert type: is this an early stop, late stop, or alt start?
      # Use the first entry's alert_code as representative
      my $alert_code = $entries[0]{"alert_code"};
      my $alt_type;
      if($alert_code eq "cdsstopn" || $alert_code eq "cdsstopp") {
        $alt_type = "alt_stop_early";
      }
      elsif($alert_code eq "mutendex") {
        $alt_type = "alt_stop_late";
      }
      elsif($alert_code eq "mutstart") {
        $alt_type = "alt_start";
      }
      elsif($alert_code eq "mutendcd") {
        # mutendcd alone: the reference stop is invalid but no clear
        # alternative identified in this group — skip
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# Alt detect: %s %s mutendcd alone, skipping\n", $ftr_key, $mdl_coords));
        next;
      }
      else {
        next;
      }

      # Compute the new CDS coords with the alternative boundary
      my $new_coords = compute_alternative_cds_coords($original_coords, $mdl_coords, $alt_type, $FH_HR);
      if(! defined $new_coords) {
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# WARNING: $sub_name, unable to compute alt coords for %s %s %s, skipping\n",
                  $ftr_key, $mdl_coords, $alt_type));
        next;
      }

      # Don't add if the new coords are identical to the original
      if($new_coords eq $original_coords) {
        next;
      }

      # Don't add if the alternative has a different number of segments
      # than the original (e.g., original is multi-segment ribosomal
      # slippage but alternative is single-segment)
      my $orig_nsegs = scalar(split(/,/, $original_coords));
      my $new_nsegs  = scalar(split(/,/, $new_coords));
      if($orig_nsegs != $new_nsegs) {
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# Alt detect: %s %s nseg_orig=%d nseg_new=%d, segment count mismatch, skipping\n",
                  $ftr_key, $mdl_coords, $orig_nsegs, $new_nsegs));
        next;
      }

      # Don't add if the alternative CDS length deviates too much from the original
      my $orig_len = coords_total_length($original_coords);
      my $new_len  = coords_total_length($new_coords);
      if($orig_len > 0 && $new_len > 0) {
        my $fract_diff = abs($new_len - $orig_len) / $orig_len;
        if($fract_diff > $max_fract_diff) {
          ofile_OutputString($FH_HR->{"log"}, 1,
            sprintf("# Alt detect: %s %s fract_diff=%.3f > max=%.3f (orig_len=%d new_len=%d), skipping\n",
                    $ftr_key, $mdl_coords, $fract_diff, $max_fract_diff, $orig_len, $new_len));
          next;
        }
      }

      ofile_OutputString($FH_HR->{"log"}, 1,
        sprintf("# Alt detect: %s %s n_seqs=%d n_ind=%d alt_type=%s new_coords=%s ADDED\n",
                $ftr_key, $mdl_coords, scalar(@entries), $n_ind, $alt_type, $new_coords));

      push(@alt_features, {
        ftr_name        => $ftr_name,
        ftr_idx         => $ftr_idx,
        alert_code      => $alert_code,
        alt_type        => $alt_type,
        mdl_alt_coords  => $mdl_coords,
        original_coords => $original_coords,
        new_coords      => $new_coords,
        n_seqs          => scalar(@entries),
        n_independent   => $n_ind,
        accessions_AR   => \@accessions,
      });
    }
  }

  return \@alt_features;
}

#################################################################
# Subroutine: compute_alternative_cds_coords()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Given the original CDS coords and an alternative
#             stop/start position from a v-annotate alert, compute
#             the new CDS coords string.
#
#             For alt_stop_early/alt_stop_late: replace the 3' end
#             of the last segment with the alternative stop position.
#             For alt_start: replace the 5' start of the first
#             segment with the alternative start position.
#
#             Handles multi-segment CDS (ribosomal slippage) by
#             only modifying the relevant terminal segment.
#
# Arguments:
#   $original_coords: original CDS coords string (e.g., "1979..2442:+,2439..2490:+")
#   $mdl_alt_coords:  model coords of the alternative boundary (e.g., "2470..2472:+")
#   $alt_type:        "alt_stop_early", "alt_stop_late", or "alt_start"
#   $FH_HR:           REF to hash of file handles
#
# Returns: new coords string, or undef on failure
#################################################################
sub compute_alternative_cds_coords {
  my $sub_name = "compute_alternative_cds_coords";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($original_coords, $mdl_alt_coords, $alt_type, $FH_HR) = @_;

  # Parse mdl_alt_coords to get the alternative boundary position
  # Format: "2470..2472:+" — for a stop codon, the 3' position (2472) is the new CDS end
  # For a start codon, the 5' position (2470) is the new CDS start
  my ($alt_start, $alt_stop, $alt_strand);
  if($mdl_alt_coords =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
    ($alt_start, $alt_stop, $alt_strand) = ($1, $2, $3);
  }
  else {
    return undef;
  }

  # Parse original coords into segments
  my @segments = split(/,/, $original_coords);

  if($alt_type eq "alt_stop_early" || $alt_type eq "alt_stop_late") {
    # Modify the last segment's 3' end (for + strand) or 5' end (for - strand)
    my $last_seg = $segments[$#segments];
    if($last_seg =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
      my ($seg_start, $seg_stop, $seg_strand) = ($1, $2, $3);
      if($seg_strand eq "+") {
        # For + strand: the 3' end is the higher coordinate
        $segments[$#segments] = $seg_start . ".." . $alt_stop . ":" . $seg_strand;
      }
      else {
        # For - strand: the 3' end is the lower coordinate
        $segments[$#segments] = $alt_start . ".." . $seg_stop . ":" . $seg_strand;
      }
    }
    else {
      return undef;
    }
  }
  elsif($alt_type eq "alt_start") {
    # Modify the first segment's 5' start
    my $first_seg = $segments[0];
    if($first_seg =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
      my ($seg_start, $seg_stop, $seg_strand) = ($1, $2, $3);
      if($seg_strand eq "+") {
        # For + strand: the 5' start is the lower coordinate
        $segments[0] = $alt_start . ".." . $seg_stop . ":" . $seg_strand;
      }
      else {
        # For - strand: the 5' start is the higher coordinate
        $segments[0] = $seg_start . ".." . $alt_stop . ":" . $seg_strand;
      }
    }
    else {
      return undef;
    }
  }
  else {
    return undef;
  }

  # Validate: for each segment, start <= stop for + strand, start >= stop for - strand
  foreach my $seg (@segments) {
    if($seg =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
      my ($s, $e, $st) = ($1, $2, $3);
      if($st eq "+" && $s > $e) { return undef; }  # inverted + strand coords
      if($st eq "-" && $s < $e) { return undef; }  # inverted - strand coords
    }
  }

  return join(",", @segments);
}

#################################################################
# Subroutine: parse_alt_for_exceptions()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Parse a .vadr.alt file and extract alerts that could
#             warrant adding _exc (exception) keys to the minfo.
#             Groups alerts by feature and exception type, merging
#             overlapping model coordinate regions.
#
#             Alert-to-exception mapping:
#               fsthicft,fsthicfi,fstukcft,fstukcfi,
#                 fstlocft,fstlocfi              -> fst_exc (FEATURE)
#               insertnp,insertnn                -> insertn_exc (FEATURE)
#               deletinp,deletinn                -> deletin_exc (FEATURE)
#               lowsimic,lowsim5c,lowsim3c,
#                 lowsim5s,lowsim3s              -> lowsim_exc (MODEL)
#
# Arguments:
#   $alt_file: path to .vadr.alt file
#   $FH_HR:    REF to hash of file handles
#
# Returns: REF to hash of hash of arrays:
#   key1: "<exc_type>:<ftr_name>:<ftr_idx>" (or "<exc_type>:MODEL:-1" for lowsim)
#   key2: "<mdl_coords>"
#   value: array ref of accession strings
#################################################################
sub parse_alt_for_exceptions {
  my $sub_name = "parse_alt_for_exceptions";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($alt_file, $FH_HR) = @_;

  # Map alert codes to exception types
  my %alert_to_exc = (
    "fsthicft" => "fst_exc", "fsthicfi" => "fst_exc",
    "fstukcft" => "fst_exc", "fstukcfi" => "fst_exc",
    "fstlocft" => "fst_exc", "fstlocfi" => "fst_exc",
    "insertnp" => "insertn_exc", "insertnn" => "insertn_exc",
    "deletinp" => "deletin_exc", "deletinn" => "deletin_exc",
    "lowsimic" => "lowsim_exc", "lowsim5c" => "lowsim_exc",
    "lowsim3c" => "lowsim_exc", "lowsim5s" => "lowsim_exc",
    "lowsim3s" => "lowsim_exc",
  );

  my %exc_groups = ();

  open(my $fh, "<", $alt_file) || ofile_FAIL("ERROR in $sub_name, unable to open $alt_file for reading", 1, $FH_HR);
  while(my $line = <$fh>) {
    chomp $line;
    next if($line =~ /^\#/ || $line =~ /^\s*$/);

    my @tok = split(/\s+/, $line);
    next if(scalar(@tok) < 13);

    my $seq_name   = $tok[1];
    my $ftr_type   = $tok[3];
    my $ftr_name   = $tok[4];
    my $ftr_idx    = $tok[5] - 1;  # convert from 1-based (.vadr.alt) to 0-based (internal)
    my $alert_code = $tok[6];
    my $mdl_coords = $tok[11];

    next if(! exists $alert_to_exc{$alert_code});
    next if($mdl_coords eq "-");

    my $exc_type = $alert_to_exc{$alert_code};

    # lowsim_exc goes on MODEL line, others go on FEATURE line
    my $exc_key;
    if($exc_type eq "lowsim_exc") {
      $exc_key = "lowsim_exc:MODEL:-1";
    }
    else {
      $exc_key = $exc_type . ":" . $ftr_name . ":" . $ftr_idx;
    }

    # For insertn_exc and deletin_exc, extract the actual size from
    # the alert detail (format: [N>M] where N is actual, M is max allowed)
    my $exc_value = 0;
    if($exc_type eq "insertn_exc" || $exc_type eq "deletin_exc") {
      my $detail = join(" ", @tok[13..$#tok]);
      if($detail =~ /\[(\d+)>/) {
        $exc_value = $1;
      }
    }

    if(! exists $exc_groups{$exc_key}{$mdl_coords}) {
      $exc_groups{$exc_key}{$mdl_coords} = [];
    }
    push(@{$exc_groups{$exc_key}{$mdl_coords}}, { acc => $seq_name, value => $exc_value });
  }
  close($fh);

  return \%exc_groups;
}

#################################################################
# Subroutine: detect_exceptions()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Given grouped exception-worthy alerts from
#             parse_alt_for_exceptions(), apply the independence
#             threshold and determine which _exc keys to add.
#             Merges overlapping model coordinate regions for the
#             same feature and exception type.
#
# Arguments:
#   $exc_groups_HHR: REF to hash from parse_alt_for_exceptions()
#   $min_independent: minimum independent observations
#   $FH_HR:          REF to hash of file handles
#
# Returns: REF to array of hashrefs, each with keys:
#          exc_type, ftr_name, ftr_idx, exc_coords, n_seqs,
#          n_independent
#################################################################
sub detect_exceptions {
  my $sub_name = "detect_exceptions";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($exc_groups_HHR, $min_independent, $FH_HR) = @_;

  my @exceptions = ();

  foreach my $exc_key (sort keys %{$exc_groups_HHR}) {
    my ($exc_type, $ftr_name, $ftr_idx) = split(/:/, $exc_key);

    # Collect all regions and their accessions/values across all mdl_coords
    # for this exception type + feature, then merge overlapping regions
    my @all_regions = (); # array of [start, stop, strand, \@entries, max_value]

    foreach my $mdl_coords (sort keys %{$exc_groups_HHR->{$exc_key}}) {
      my @entries = @{$exc_groups_HHR->{$exc_key}{$mdl_coords}};

      if($mdl_coords =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
        # Track the max value across all entries for this coord
        my $max_val = 0;
        foreach my $e (@entries) { if($e->{"value"} > $max_val) { $max_val = $e->{"value"}; } }
        push(@all_regions, [$1, $2, $3, \@entries, $max_val]);
      }
    }

    # Merge overlapping regions (same strand)
    @all_regions = sort { $a->[0] <=> $b->[0] } @all_regions;
    my @merged = ();
    foreach my $region (@all_regions) {
      my ($rstart, $rstop, $rstrand, $rentries, $rmax_val) = @{$region};
      if(scalar(@merged) > 0 &&
         $merged[$#merged][2] eq $rstrand &&
         $rstart <= $merged[$#merged][1] + 1) {
        # Overlapping or adjacent: extend and merge
        if($rstop > $merged[$#merged][1]) {
          $merged[$#merged][1] = $rstop;
        }
        push(@{$merged[$#merged][3]}, @{$rentries});
        if($rmax_val > $merged[$#merged][4]) { $merged[$#merged][4] = $rmax_val; }
      }
      else {
        push(@merged, [$rstart, $rstop, $rstrand, [@{$rentries}], $rmax_val]);
      }
    }

    # Apply independence threshold to each merged region
    foreach my $region (@merged) {
      my ($rstart, $rstop, $rstrand, $rentries, $rmax_val) = @{$region};

      # Deduplicate accessions (same seq may trigger multiple alerts in region)
      my %unique_accs = map { $_->{"acc"} => 1 } @{$rentries};
      my @unique_acc_list = keys %unique_accs;
      my $n_ind = count_independent_observations(\@unique_acc_list);

      if($n_ind < $min_independent) {
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# Exc detect: %s %s %d..%d:%s n_seqs=%d n_ind=%d < min=%d, skipping\n",
                  $exc_type, $ftr_name, $rstart, $rstop, $rstrand,
                  scalar(@unique_acc_list), $n_ind, $min_independent));
        next;
      }

      # For insertn_exc and deletin_exc: use coords-value format (coords:maxlen)
      # For fst_exc and lowsim_exc: use coords-only format
      my $exc_coords;
      if($exc_type eq "insertn_exc" || $exc_type eq "deletin_exc") {
        $exc_coords = $rstart . ".." . $rstop . ":" . $rstrand . ":" . $rmax_val;
      }
      else {
        $exc_coords = $rstart . ".." . $rstop . ":" . $rstrand;
      }

      ofile_OutputString($FH_HR->{"log"}, 1,
        sprintf("# Exc detect: %s %s %s n_seqs=%d n_ind=%d ADDED\n",
                $exc_type, $ftr_name, $exc_coords, scalar(@unique_acc_list), $n_ind));

      push(@exceptions, {
        exc_type      => $exc_type,
        ftr_name      => $ftr_name,
        ftr_idx       => int($ftr_idx),
        exc_coords    => $exc_coords,
        n_seqs        => scalar(@unique_acc_list),
        n_independent => $n_ind,
      });
    }
  }

  return \@exceptions;
}

#################################################################
# Subroutine: add_alternatives_and_exceptions_to_minfo()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Read an existing minfo file and write a new one with
#             alternative CDS features and _exc exception keys added.
#
#             For each alternative CDS:
#             - Adds alternative_ftr_set:"<name>(cds)" to the original
#               CDS FEATURE line
#             - Inserts a new CDS FEATURE line with the alternative
#               coords and the same alternative_ftr_set value
#             - If a matching gene feature exists, adds alternative
#               gene features with alternative_ftr_set_subn links
#
#             For each exception:
#             - Appends the _exc key:"value" to the relevant
#               FEATURE or MODEL line
#
# Arguments:
#   $in_minfo_file:   path to input .minfo file
#   $out_minfo_file:  path to output .minfo file
#   $alt_features_AR: REF to array of alternative feature hashrefs
#                     from detect_alternative_features()
#   $exceptions_AR:   REF to array of exception hashrefs
#                     from detect_exceptions()
#   $model_key:       model name
#   $FH_HR:           REF to hash of file handles
#
# Returns: void
#################################################################
sub add_alternatives_and_exceptions_to_minfo {
  my $sub_name = "add_alternatives_and_exceptions_to_minfo";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($in_minfo_file, $out_minfo_file, $alt_features_AR, $exceptions_AR, $model_key, $FH_HR) = @_;

  # Read all lines from input minfo
  my @lines = ();
  open(my $infh, "<", $in_minfo_file) || ofile_FAIL("ERROR in $sub_name, unable to open $in_minfo_file for reading", 1, $FH_HR);
  while(my $line = <$infh>) {
    chomp $line;
    push(@lines, $line);
  }
  close($infh);

  # Index FEATURE lines by their 0-based feature index (order they appear)
  # and extract key info (type, gene name, coords)
  my @ftr_line_idx = ();  # ftr_line_idx[$ftr_idx] = index into @lines
  my @ftr_type = ();
  my @ftr_gene = ();
  my @ftr_coords = ();
  my $ftr_count = 0;
  for(my $i = 0; $i < scalar(@lines); $i++) {
    if($lines[$i] =~ /^FEATURE\s/) {
      $ftr_line_idx[$ftr_count] = $i;
      # Parse type
      if($lines[$i] =~ /type:"([^"]+)"/) { $ftr_type[$ftr_count] = $1; }
      else { $ftr_type[$ftr_count] = ""; }
      # Parse gene
      if($lines[$i] =~ /gene:"([^"]+)"/) { $ftr_gene[$ftr_count] = $1; }
      else { $ftr_gene[$ftr_count] = ""; }
      # Parse coords
      if($lines[$i] =~ /coords:"([^"]+)"/) { $ftr_coords[$ftr_count] = $1; }
      else { $ftr_coords[$ftr_count] = ""; }
      $ftr_count++;
    }
  }

  # --- Step 1: Apply exceptions to existing lines ---
  # Index for MODEL line
  my $model_line_idx = -1;
  for(my $i = 0; $i < scalar(@lines); $i++) {
    if($lines[$i] =~ /^MODEL\s/) { $model_line_idx = $i; last; }
  }

  foreach my $exc (@{$exceptions_AR}) {
    my $exc_type   = $exc->{"exc_type"};
    my $exc_coords = $exc->{"exc_coords"};
    my $ftr_idx    = $exc->{"ftr_idx"};

    if($exc_type eq "lowsim_exc") {
      # Append to MODEL line
      if($model_line_idx >= 0) {
        # Check if lowsim_exc already exists on the line
        if($lines[$model_line_idx] =~ /lowsim_exc:"([^"]*)"/) {
          # Append to existing value
          my $existing = $1;
          $lines[$model_line_idx] =~ s/lowsim_exc:"[^"]*"/lowsim_exc:"$existing,$exc_coords"/;
        }
        else {
          $lines[$model_line_idx] .= " lowsim_exc:\"$exc_coords\"";
        }
      }
    }
    else {
      # Append to FEATURE line at ftr_idx
      if($ftr_idx >= 0 && $ftr_idx < $ftr_count) {
        my $li = $ftr_line_idx[$ftr_idx];
        # Check if this exc_type already exists on the line
        if($lines[$li] =~ /$exc_type:"([^"]*)"/) {
          my $existing = $1;
          $lines[$li] =~ s/$exc_type:"[^"]*"/$exc_type:"$existing,$exc_coords"/;
        }
        else {
          $lines[$li] .= " $exc_type:\"$exc_coords\"";
        }
      }
    }
  }

  # --- Step 2: Build alternative feature insertions ---
  # Group alternatives by ftr_idx so we can handle multiple alternatives
  # for the same CDS
  my %alts_by_ftr_idx = ();
  foreach my $alt (@{$alt_features_AR}) {
    my $ftr_idx = $alt->{"ftr_idx"};
    if(! exists $alts_by_ftr_idx{$ftr_idx}) {
      $alts_by_ftr_idx{$ftr_idx} = [];
    }
    push(@{$alts_by_ftr_idx{$ftr_idx}}, $alt);
  }

  # For each CDS with alternatives, prepare:
  # 1. Modify the original CDS line to add alternative_ftr_set
  # 2. Create new CDS line(s) for each alternative
  # 3. Find and modify the parent gene (if any)
  # 4. Create new gene line(s)
  #
  # We build a map of: line_index -> [lines to insert AFTER this line]
  # and lines to modify in-place
  my %insert_after = ();  # line_idx -> [line1, line2, ...]
  my %used_set_names = (); # track used alternative_ftr_set names for uniqueness

  foreach my $ftr_idx (sort { $a <=> $b } keys %alts_by_ftr_idx) {
    my @alts = @{$alts_by_ftr_idx{$ftr_idx}};
    my $cds_line_idx = $ftr_line_idx[$ftr_idx];

    # Generate a unique alternative_ftr_set name from the gene/feature name
    my $gene_name = $ftr_gene[$ftr_idx];
    my $set_base = generate_alt_set_name($gene_name, $alts[0]{"ftr_name"});
    my $cds_set_name = $set_base . "(cds)";
    my $gene_set_name = $set_base . "(gene)";
    # Ensure uniqueness: if this name is already taken, append .N
    if(exists $used_set_names{$cds_set_name}) {
      my $n = 2;
      while(exists $used_set_names{$set_base . "." . $n . "(cds)"}) { $n++; }
      $set_base = $set_base . "." . $n;
      $cds_set_name = $set_base . "(cds)";
      $gene_set_name = $set_base . "(gene)";
    }
    $used_set_names{$cds_set_name} = 1;
    $used_set_names{$gene_set_name} = 1;

    # Add alternative_ftr_set to original CDS line (skip if already present)
    if($lines[$cds_line_idx] !~ /alternative_ftr_set:/) {
      $lines[$cds_line_idx] .= " alternative_ftr_set:\"$cds_set_name\"";
    }

    # Create new CDS lines for each alternative
    # Skip if the alternative coords match the original (no real change)
    # or match any other existing CDS coords in the minfo
    my %existing_coords = ();
    for(my $fi = 0; $fi < $ftr_count; $fi++) {
      if($ftr_type[$fi] eq "CDS") { $existing_coords{$ftr_coords[$fi]} = 1; }
    }
    my @new_cds_lines = ();
    foreach my $alt (@alts) {
      if(exists $existing_coords{$alt->{"new_coords"}}) {
        next;  # skip: coords already exist in minfo
      }
      my $new_line = $lines[$cds_line_idx];
      # Replace coords with alternative coords
      $new_line =~ s/coords:"[^"]*"/coords:"$alt->{"new_coords"}"/;
      push(@new_cds_lines, $new_line);
      $existing_coords{$alt->{"new_coords"}} = 1;  # track newly added too
    }

    # Insert new CDS lines after the original
    $insert_after{$cds_line_idx} = \@new_cds_lines;

    # Find the parent gene feature (same gene name, type "gene")
    my $gene_ftr_idx = -1;
    for(my $fi = 0; $fi < $ftr_count; $fi++) {
      if($ftr_type[$fi] eq "gene" && $ftr_gene[$fi] eq $gene_name && $gene_name ne "") {
        $gene_ftr_idx = $fi;
        last;
      }
    }

    if($gene_ftr_idx >= 0) {
      my $gene_line_idx = $ftr_line_idx[$gene_ftr_idx];

      # Check if the existing gene already encompasses ALL alternative CDS variants.
      # If so, no gene alternatives are needed (the gene spans all variants).
      my $gene_already_spans_all = 1;
      if($ftr_coords[$gene_ftr_idx] =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
        my ($gstart, $gstop, $gstrand) = ($1, $2, $3);
        foreach my $alt (@alts) {
          my $new_gene_coords = compute_gene_coords_for_alt($ftr_coords[$gene_ftr_idx],
                                                             $alt->{"original_coords"},
                                                             $alt->{"new_coords"},
                                                             $alt->{"alt_type"});
          if($new_gene_coords ne $ftr_coords[$gene_ftr_idx]) {
            $gene_already_spans_all = 0;
            last;
          }
        }
      }

      if(! $gene_already_spans_all) {
        # Gene coords need to change for some alternatives — add gene alternatives
        # Add alternative_ftr_set and subn to original gene line
        # The original gene corresponds to CDS alternative .1 (the original)
        if($lines[$gene_line_idx] !~ /alternative_ftr_set:/) {
          $lines[$gene_line_idx] .= " alternative_ftr_set:\"$gene_set_name\" alternative_ftr_set_subn:\"$cds_set_name.1\"";
        }

        # Create gene alternatives only for CDS alternatives that were
        # actually added (not skipped by duplicate-coords check).
        # The number of added CDS = scalar(@new_cds_lines), and
        # total CDS in set = 1 (original) + scalar(@new_cds_lines).
        # Gene subn indices must match: .1 = original, .2..N = alternatives.
        my @new_gene_lines = ();
        my $alt_num = 2;  # .1 is the original, alternatives start at .2
        foreach my $alt (@alts) {
          # Only create gene alternative if this CDS alt was actually added
          next if(exists $existing_coords{$alt->{"new_coords"}} &&
                  ! grep { /coords:"\Q$alt->{"new_coords"}\E"/ } @new_cds_lines);
          my $new_gene_line = $lines[$gene_line_idx];
          my $new_gene_coords = compute_gene_coords_for_alt($ftr_coords[$gene_ftr_idx],
                                                             $alt->{"original_coords"},
                                                             $alt->{"new_coords"},
                                                             $alt->{"alt_type"});
          $new_gene_line =~ s/coords:"[^"]*"/coords:"$new_gene_coords"/;
          $new_gene_line =~ s/alternative_ftr_set_subn:"[^"]*"/alternative_ftr_set_subn:"$cds_set_name.$alt_num"/;
          push(@new_gene_lines, $new_gene_line);
          $alt_num++;
        }

        # Insert new gene lines after the original gene line
        if(exists $insert_after{$gene_line_idx}) {
          push(@{$insert_after{$gene_line_idx}}, @new_gene_lines);
        }
        else {
          $insert_after{$gene_line_idx} = \@new_gene_lines;
        }
      } # end of if(! $gene_already_spans_all)
    } # end of if($gene_ftr_idx >= 0)
  }

  # --- Step 3: Write output minfo ---
  open(my $outfh, ">", $out_minfo_file) || ofile_FAIL("ERROR in $sub_name, unable to open $out_minfo_file for writing", 1, $FH_HR);
  for(my $i = 0; $i < scalar(@lines); $i++) {
    print $outfh $lines[$i] . "\n";
    if(exists $insert_after{$i}) {
      foreach my $new_line (@{$insert_after{$i}}) {
        print $outfh $new_line . "\n";
      }
    }
  }
  close($outfh);

  return;
}

#################################################################
# Subroutine: generate_alt_set_name()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Generate a short alternative_ftr_set base name from
#             the gene name and/or product name. Convention from
#             existing VADR models: lowercase abbreviation.
#
# Arguments:
#   $gene_name: gene qualifier value (e.g., "N", "SH", "P")
#   $ftr_name:  feature/product name (e.g., "nucleocapsid_protein")
#
# Returns: base name string (e.g., "n", "sh", "p")
#################################################################
sub generate_alt_set_name {
  my $sub_name = "generate_alt_set_name";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($gene_name, $ftr_name) = @_;

  # Prefer gene name if available
  if(defined $gene_name && $gene_name ne "") {
    return lc($gene_name);
  }

  # Fall back to first word of product name
  if(defined $ftr_name && $ftr_name ne "") {
    my $base = lc($ftr_name);
    $base =~ s/_.*//;  # take first word
    return $base;
  }

  return "alt";
}

#################################################################
# Subroutine: compute_gene_coords_for_alt()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Compute the gene coords for an alternative CDS by
#             adjusting the original gene coords to encompass the
#             alternative CDS coords.
#
# Arguments:
#   $gene_coords:    original gene coords (e.g., "146..1795:+")
#   $orig_cds_coords: original CDS coords
#   $alt_cds_coords:  alternative CDS coords
#   $alt_type:       "alt_stop_early", "alt_stop_late", or "alt_start"
#
# Returns: new gene coords string
#################################################################
sub compute_gene_coords_for_alt {
  my $sub_name = "compute_gene_coords_for_alt";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($gene_coords, $orig_cds_coords, $alt_cds_coords, $alt_type) = @_;

  # For simple cases: gene coords are single-segment
  # Adjust the gene end/start to match the CDS change
  if($gene_coords =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
    my ($gstart, $gstop, $gstrand) = ($1, $2, $3);

    # Parse the last segment of original and alternative CDS
    my @orig_segs = split(/,/, $orig_cds_coords);
    my @alt_segs = split(/,/, $alt_cds_coords);

    if($alt_type eq "alt_stop_early" || $alt_type eq "alt_stop_late") {
      # Get the 3' end of the last CDS segment
      if($alt_segs[$#alt_segs] =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
        my $alt_end = ($3 eq "+") ? $2 : $1;
        if($gstrand eq "+") {
          $gstop = $alt_end;
        }
        else {
          $gstart = $alt_end;
        }
      }
    }
    elsif($alt_type eq "alt_start") {
      # Get the 5' start of the first CDS segment
      if($alt_segs[0] =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
        my $alt_start = ($3 eq "+") ? $1 : $2;
        if($gstrand eq "+") {
          $gstart = $alt_start;
        }
        else {
          $gstop = $alt_start;
        }
      }
    }

    return $gstart . ".." . $gstop . ":" . $gstrand;
  }

  # Multi-segment gene coords: just return the original (rare case)
  return $gene_coords;
}

#################################################################
# Subroutine: translate_alternative_cds_proteins()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    For each detected alternative CDS feature, translate
#             the protein from a representative sequence and append
#             it to the seed model's protein.fa file. Then rebuild
#             the BLAST protein database.
#
#             For cdsstopn (early stop): uses the n_instp field from
#             .vadr.ftr which gives the exact sequence position of
#             the actual in-frame stop codon.
#
#             For mutendex (late stop): uses the seq_coords from
#             .vadr.alt which gives the stop codon position.
#
# Arguments:
#   $alt_features_AR:  REF to array of alternative feature hashrefs
#   $ftr_file:         path to .vadr.ftr from first v-annotate pass
#   $alt_file:         path to .vadr.alt from first v-annotate pass
#   $fasta_file:       path to tier 2 FASTA (indexed with .ssi)
#   $protein_fa_file:  path to seed model protein.fa (will be appended)
#   $model_key:        model name
#   $execs_HR:         REF to hash of executables
#   $opt_HHR:          REF to hash of option hashes
#   $FH_HR:            REF to hash of file handles
#
# Returns: void (appends to protein_fa_file and rebuilds BLAST db)
#################################################################
sub translate_alternative_cds_proteins {
  my $sub_name = "translate_alternative_cds_proteins";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($alt_features_AR, $ftr_file, $alt_file, $fasta_file, $protein_fa_file,
      $model_key, $execs_HR, $opt_HHR, $FH_HR) = @_;

  if(scalar(@{$alt_features_AR}) == 0) { return; }

  # Parse .vadr.ftr: for each (accession, ftr_idx), store:
  #   seq_coords (the multi-segment coords string)
  #   n_from, n_to (CDS start/end in sequence)
  #   n_instp (in-frame stop position, or "-" if none)
  my %ftr_info_HH = ();
  open(my $ftr_fh, "<", $ftr_file) || ofile_FAIL("ERROR in $sub_name, unable to open $ftr_file", 1, $FH_HR);
  while(my $line = <$ftr_fh>) {
    chomp $line;
    next if($line =~ /^\#/ || $line =~ /^\s*$/);
    my @tok = split(/\s+/, $line);
    next if(scalar(@tok) < 25);
    my $acc        = $tok[1];
    my $ftr_type   = $tok[5];
    my $ftr_idx    = $tok[8] - 1;  # convert from 1-based (.vadr.ftr) to 0-based (internal)
    my $n_from     = $tok[11];
    my $n_to       = $tok[12];
    my $n_instp    = $tok[13];
    my $seq_coords = $tok[23];
    next if($ftr_type ne "CDS");
    $ftr_info_HH{$acc}{$ftr_idx} = {
      seq_coords => $seq_coords,
      n_from     => $n_from,
      n_to       => $n_to,
      n_instp    => $n_instp,
    };
  }
  close($ftr_fh);

  # Parse .vadr.alt: for each (accession, ftr_idx, alert_code), store seq_coords
  # This gives the exact position of the alternative boundary in sequence coords
  my %alt_seq_coords_HH = ();
  open(my $alt_fh, "<", $alt_file) || ofile_FAIL("ERROR in $sub_name, unable to open $alt_file", 1, $FH_HR);
  while(my $line = <$alt_fh>) {
    chomp $line;
    next if($line =~ /^\#/ || $line =~ /^\s*$/);
    my @tok = split(/\s+/, $line);
    next if(scalar(@tok) < 13);
    my $acc        = $tok[1];
    my $ftr_idx    = $tok[5] - 1;  # convert from 1-based (.vadr.alt) to 0-based (internal)
    my $alert_code = $tok[6];
    my $seq_coords = $tok[9];
    # Store: key is "acc:ftr_idx:alert_code"
    $alt_seq_coords_HH{"$acc:$ftr_idx:$alert_code"} = $seq_coords;
  }
  close($alt_fh);

  # Index the FASTA file for esl-sfetch if not already indexed
  my $ssi_file = $fasta_file . ".ssi";
  if(! -e $ssi_file) {
    utl_RunCommand($execs_HR->{"esl-sfetch"} . " --index " . $fasta_file,
                   opt_Get("-v", $opt_HHR), 0, $FH_HR);
  }

  # Open protein.fa for appending
  open(my $prot_fh, ">>", $protein_fa_file) || ofile_FAIL("ERROR in $sub_name, unable to open $protein_fa_file for appending", 1, $FH_HR);

  foreach my $alt (@{$alt_features_AR}) {
    my $ftr_idx    = $alt->{"ftr_idx"};
    my $new_coords = $alt->{"new_coords"};
    my $alt_type   = $alt->{"alt_type"};
    my $alert_code = $alt->{"alert_code"};
    my @accessions = @{$alt->{"accessions_AR"}};

    # Find a representative sequence
    my $rep_acc = undef;
    foreach my $acc (@accessions) {
      if(exists $ftr_info_HH{$acc}{$ftr_idx}) {
        $rep_acc = $acc;
        last;
      }
    }
    if(! defined $rep_acc) {
      ofile_OutputString($FH_HR->{"log"}, 1,
        sprintf("# WARNING: $sub_name, no representative for alt CDS ftr_idx=%d, skipping\n", $ftr_idx));
      next;
    }

    my $ftr_info = $ftr_info_HH{$rep_acc}{$ftr_idx};
    my $orig_seq_coords = $ftr_info->{"seq_coords"};

    # Build the alternative sequence coords
    # Strategy: take the original seq_coords and replace the stop/start
    # boundary using the EXACT position from .vadr.ftr or .vadr.alt
    my $alt_seq_coords;
    if($alt_type eq "alt_stop_early") {
      # cdsstopn: n_instp from .vadr.ftr gives the exact early stop position
      my $instp = $ftr_info->{"n_instp"};
      if(!defined $instp || $instp eq "-") {
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# WARNING: $sub_name, no n_instp for %s ftr_idx=%d, skipping\n", $rep_acc, $ftr_idx));
        next;
      }
      $alt_seq_coords = replace_last_segment_end($orig_seq_coords, int($instp));
    }
    elsif($alt_type eq "alt_stop_late") {
      # mutendex: seq_coords from .vadr.alt gives the late stop codon position
      my $key = "$rep_acc:$ftr_idx:mutendex";
      if(! exists $alt_seq_coords_HH{$key}) {
        ofile_OutputString($FH_HR->{"log"}, 1,
          sprintf("# WARNING: $sub_name, no mutendex seq_coords for %s ftr_idx=%d, skipping\n", $rep_acc, $ftr_idx));
        next;
      }
      # Parse the stop codon end position from the alt seq_coords (e.g., "2536..2538:+")
      my $alt_stop_coords = $alt_seq_coords_HH{$key};
      if($alt_stop_coords =~ /^\d+\.\.(\d+):[\+\-]$/) {
        $alt_seq_coords = replace_last_segment_end($orig_seq_coords, $1);
      }
      else { next; }
    }
    elsif($alt_type eq "alt_start") {
      my $key = "$rep_acc:$ftr_idx:mutstart";
      if(! exists $alt_seq_coords_HH{$key}) { next; }
      my $alt_start_coords = $alt_seq_coords_HH{$key};
      if($alt_start_coords =~ /^(\d+)\.\.\d+:[\+\-]$/) {
        $alt_seq_coords = replace_first_segment_start($orig_seq_coords, $1);
      }
      else { next; }
    }
    else { next; }

    if(! defined $alt_seq_coords) {
      ofile_OutputString($FH_HR->{"log"}, 1,
        sprintf("# WARNING: $sub_name, failed to build alt seq coords for %s ftr_idx=%d, skipping\n", $rep_acc, $ftr_idx));
      next;
    }

    # Extract each segment of the CDS from the FASTA and concatenate
    my $cds_seq = "";
    my @segments = split(/,/, $alt_seq_coords);
    foreach my $seg (@segments) {
      if($seg =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
        my ($sstart, $sstop, $sstrand) = ($1, $2, $3);
        my $sfetch_coords = ($sstrand eq "+") ? "$sstart..$sstop" : "$sstop..$sstart";
        my $seg_seq = `$execs_HR->{"esl-sfetch"} -c $sfetch_coords $fasta_file $rep_acc 2>/dev/null`;
        $seg_seq =~ s/^>.*\n//;
        $seg_seq =~ s/\s//g;
        $cds_seq .= $seg_seq;
      }
    }

    if(length($cds_seq) == 0) {
      ofile_OutputString($FH_HR->{"log"}, 1,
        sprintf("# WARNING: $sub_name, empty CDS for %s ftr_idx=%d, skipping\n", $rep_acc, $ftr_idx));
      next;
    }

    # Translate: write to temp, run esl-translate, find full-length ORF
    my $tmp_cds = $protein_fa_file . ".tmp.alt.cds.fa";
    open(my $tmp_fh, ">", $tmp_cds) || ofile_FAIL("ERROR in $sub_name, unable to write $tmp_cds", 1, $FH_HR);
    print $tmp_fh ">$rep_acc/$alt_seq_coords\n$cds_seq\n";
    close($tmp_fh);

    my $tmp_prot = $protein_fa_file . ".tmp.alt.prot.fa";
    utl_RunCommand($execs_HR->{"esl-translate"} . " $tmp_cds > $tmp_prot",
                   opt_Get("-v", $opt_HHR), 0, $FH_HR);

    # Find ORF starting at position 1 (the full-length translation)
    my ($protein_seq, $protein_header);
    open(my $prot_in, "<", $tmp_prot) || ofile_FAIL("ERROR in $sub_name, unable to read $tmp_prot", 1, $FH_HR);
    my ($ch, $cs) = ("", "");
    while(my $pl = <$prot_in>) {
      chomp $pl;
      if($pl =~ /^>/) {
        if($ch =~ /coords=1\.\./) { $protein_seq = $cs; $protein_header = $ch; }
        $ch = $pl; $cs = "";
      } else { $cs .= $pl; }
    }
    if($ch =~ /coords=1\.\./ && ! defined $protein_seq) { $protein_seq = $cs; $protein_header = $ch; }
    close($prot_in);

    if(! defined $protein_seq) {
      ofile_OutputString($FH_HR->{"log"}, 1,
        sprintf("# WARNING: $sub_name, no full-length ORF for alt CDS %s ftr_idx=%d (cds_len=%d), skipping\n",
                $rep_acc, $ftr_idx, length($cds_seq)));
      next;
    }

    # Append to protein.fa
    my $prot_name = "$model_key/$new_coords";
    print $prot_fh ">$prot_name\n$protein_seq\n";

    ofile_OutputString($FH_HR->{"log"}, 1,
      sprintf("# Alt protein: %s from %s (%d aa)\n", $prot_name, $rep_acc, length($protein_seq)));

    unlink($tmp_cds);
    unlink($tmp_prot);
  }

  close($prot_fh);

  # Rebuild BLAST protein database
  sqf_BlastDbCreate($execs_HR->{"makeblastdb"}, "prot", $protein_fa_file, $opt_HHR, $FH_HR);

  return;
}

#################################################################
# Subroutine: replace_last_segment_end()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Replace the 3' end position of the last segment in a
#             VADR coords string with a new position.
#
# Arguments:
#   $coords:  VADR coords string (e.g., "1979..2442:+,2439..2490:+")
#   $new_end: new 3' end position for + strand (or 5' for - strand)
#
# Returns: updated coords string, or undef on failure
#################################################################
sub replace_last_segment_end {
  my ($coords, $new_end) = @_;
  my @segs = split(/,/, $coords);
  if($segs[$#segs] =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
    my ($ss, $se, $st) = ($1, $2, $3);
    if($st eq "+") { $segs[$#segs] = "$ss..$new_end:$st"; }
    else           { $segs[$#segs] = "$new_end..$se:$st"; }
    return join(",", @segs);
  }
  return undef;
}

#################################################################
# Subroutine: replace_first_segment_start()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Replace the 5' start position of the first segment
#             in a VADR coords string with a new position.
#
# Arguments:
#   $coords:    VADR coords string
#   $new_start: new 5' start position for + strand
#
# Returns: updated coords string, or undef on failure
#################################################################
sub replace_first_segment_start {
  my ($coords, $new_start) = @_;
  my @segs = split(/,/, $coords);
  if($segs[0] =~ /^(\d+)\.\.(\d+):([\+\-])$/) {
    my ($ss, $se, $st) = ($1, $2, $3);
    if($st eq "+") { $segs[0] = "$new_start..$se:$st"; }
    else           { $segs[0] = "$ss..$new_start:$st"; }
    return join(",", @segs);
  }
  return undef;
}

#################################################################
# Subroutine: identify_rerun_candidates()
# Incept:     EPN, Thu Mar 27 2026
#
# Purpose:    Identify sequences that failed the first v-annotate
#             pass ONLY due to alerts that are now addressed by
#             the detected alternatives and exceptions. These
#             sequences are candidates for a second v-annotate pass.
#
#             A sequence is a re-run candidate if ALL of its fatal
#             alerts are of types that would be resolved by the
#             new alternatives or exceptions.
#
# Arguments:
#   $alt_file:         path to .vadr.alt from first pass
#   $decision_HR:      REF to decision hash (has status for each accession)
#   $alt_features_AR:  REF to array of detected alternative features
#   $exceptions_AR:    REF to array of detected exceptions
#   $FH_HR:            REF to hash of file handles
#
# Returns: array of accession strings eligible for re-run
#################################################################
sub identify_rerun_candidates {
  my $sub_name = "identify_rerun_candidates";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($alt_file, $decision_HR, $alt_features_AR, $exceptions_AR, $FH_HR) = @_;

  # Build a set of (ftr_idx, alert_code) pairs that are addressed
  # by the new alternatives
  my %addressed_alt_H = ();
  foreach my $alt (@{$alt_features_AR}) {
    my $ftr_idx = $alt->{"ftr_idx"};
    # These alert codes on this feature are addressed by the new alternative
    foreach my $code (qw(mutendex cdsstopn cdsstopp mutendcd mutstart)) {
      $addressed_alt_H{"$ftr_idx:$code"} = 1;
    }
  }

  # Build a set of (ftr_idx, alert_code) pairs addressed by exceptions
  # Map exception types back to the alert codes they suppress
  my %exc_to_alerts = (
    "fst_exc"     => [qw(fsthicft fsthicfi fstukcft fstukcfi fstlocft fstlocfi)],
    "insertn_exc" => [qw(insertnp insertnn)],
    "deletin_exc" => [qw(deletinp deletinn)],
    "lowsim_exc"  => [qw(lowsimic lowsim5c lowsim3c lowsim5s lowsim3s lowsim5n lowsim3n lowsimin lowsimil lowsim5l lowsim3l)],
  );
  my %addressed_exc_H = ();
  foreach my $exc (@{$exceptions_AR}) {
    my $ftr_idx = $exc->{"ftr_idx"};
    my $exc_type = $exc->{"exc_type"};
    if(exists $exc_to_alerts{$exc_type}) {
      foreach my $code (@{$exc_to_alerts{$exc_type}}) {
        if($ftr_idx == -1) {
          # MODEL-level exception (lowsim): addresses this alert on ANY feature
          $addressed_exc_H{"*:$code"} = 1;
        }
        else {
          $addressed_exc_H{"$ftr_idx:$code"} = 1;
        }
      }
    }
  }

  # Parse .vadr.alt: for each failed sequence, check if ALL fatal alerts
  # are addressed
  my %seq_fatal_alerts = ();  # {accession} => [{ftr_idx, code, fail}, ...]
  open(my $fh, "<", $alt_file) || ofile_FAIL("ERROR in $sub_name, unable to open $alt_file", 1, $FH_HR);
  while(my $line = <$fh>) {
    chomp $line;
    next if($line =~ /^\#/ || $line =~ /^\s*$/);
    my @tok = split(/\s+/, $line);
    next if(scalar(@tok) < 13);
    my $acc  = $tok[1];
    my $fidx = $tok[5] - 1;  # convert from 1-based (.vadr.alt) to 0-based (internal)
    my $code = $tok[6];
    my $fail = $tok[7];
    next if($fail ne "yes");  # only care about fatal alerts
    if(! exists $seq_fatal_alerts{$acc}) {
      $seq_fatal_alerts{$acc} = [];
    }
    push(@{$seq_fatal_alerts{$acc}}, { ftr_idx => $fidx, code => $code });
  }
  close($fh);

  # Check each failed sequence
  my @rerun_candidates = ();
  foreach my $acc (sort keys %seq_fatal_alerts) {
    # Only consider sequences that were removed by the first pass
    next if(! exists $decision_HR->{$acc});
    next if($decision_HR->{$acc}{"status"} ne "removed");
    next if($decision_HR->{$acc}{"reason_code"} ne "vadr_fail");

    my $all_addressed = 1;
    foreach my $alert (@{$seq_fatal_alerts{$acc}}) {
      my $fidx = $alert->{"ftr_idx"};
      my $code = $alert->{"code"};
      my $key_specific = "$fidx:$code";
      my $key_wildcard = "*:$code";
      if(! exists $addressed_alt_H{$key_specific} &&
         ! exists $addressed_exc_H{$key_specific} &&
         ! exists $addressed_exc_H{$key_wildcard}) {
        $all_addressed = 0;
        last;
      }
    }

    if($all_addressed) {
      push(@rerun_candidates, $acc);
    }
  }

  return @rerun_candidates;
}

#################################################################
# Subroutine: coords_total_length()
# Incept:     EPN, Mon Mar 30 2026
#
# Purpose:    Calculate the total nucleotide length of a VADR coords
#             string by summing the lengths of all segments.
#
# Arguments:
#   $coords: VADR coords string (e.g., "1979..2442:+,2439..2490:+")
#
# Returns: total length (integer), or 0 on failure
#################################################################
sub coords_total_length {
  my ($coords) = @_;
  my $total = 0;
  foreach my $seg (split(/,/, $coords)) {
    if($seg =~ /^(\d+)\.\.(\d+):[\+\-]$/) {
      $total += abs($2 - $1) + 1;
    }
  }
  return $total;
}

#################################################################
# Subroutine: determine_seqs_per_group()
# Incept:     EPN, Mon Mar 31 2026
#
# Purpose:    Determine how many sequences to select per group
#             based on the total number of non-empty groups.
#             Fewer groups -> more sequences per group to ensure
#             adequate diversity in the final alignment.
#
#             Default schedule:
#               1 group:       $nper1grp (default 5)
#               2-5 groups:    3
#               6-10 groups:   2
#               11+ groups:    1
#
# Arguments:
#   $n_groups: number of non-empty groups
#   $nper1grp: number of seqs per group when 1 group (--nper1grp)
#
# Returns: integer number of sequences to select per group
#################################################################
sub determine_seqs_per_group {
  my ($n_groups, $nper1grp) = @_;

  if($n_groups <= 1)  { return $nper1grp; }
  if($n_groups <= 5)  { return 3; }
  if($n_groups <= 10) { return 2; }
  return 1;
}
