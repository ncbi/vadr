#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Mon Feb  1 15:07:43 2016 [dnaorg_build.pl split off from dnaorg_annotate_genomes.pl]
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

require "dnaorg.pm"; 
require "epn-options.pm";

#######################################################################################
# What this script does: 
#
# Preliminaries: 
#   - process options
#   - create the output directory
#   - output program banner and open output files
#   - parse the optional input files, if necessary
#   - make sure the required executables are executable
#
# Step 1. Gather and process information on reference genome using Edirect
#
# Step 2. Fetch and process the reference genome sequence
#
# Step 3. Build and calibrate models
#######################################################################################

# first, determine the paths to all modules, scripts and executables that we'll need
# we currently use hard-coded-paths for Infernal, HMMER and easel executables:
my $inf_exec_dir      = "/usr/local/infernal/1.1.1/bin/";
my $hmmer_exec_dir    = "/usr/local/hmmer/3.1/bin/";
my $esl_exec_dir      = "/usr/local/infernal/1.1.1/bin/";

# make sure the DNAORGDIR environment variable is set
my $dnaorgdir = $ENV{'DNAORGDIR'};
if(! exists($ENV{'DNAORGDIR'})) { 
    printf STDERR ("\nERROR, the environment variable DNAORGDIR is not set, please set it to the directory where you installed the dnaorg scripts and their dependencies.\n"); 
    exit(1); 
}
if(! (-d $dnaorgdir)) { 
    printf STDERR ("\nERROR, the dnaorg directory specified by your environment variable DNAORGDIR does not exist.\n"); 
    exit(1); 
}    
 
# determine other required paths to executables relative to $dnaorgdir
my $esl_fetch_cds     = $dnaorgdir . "/esl-fetch-cds/esl-fetch-cds.pl";

#########################################################
# Command line and option processing using epn-options.pm
#
# opt_HH: 2D hash:
#         1D key: option name (e.g. "-h")
#         2D key: string denoting type of information 
#                 (one of "type", "default", "group", "requires", "incompatible", "preamble", "help")
#         value:  string explaining 2D key:
#                 "type":          "boolean", "string", "int" or "real"
#                 "default":       default value for option
#                 "group":         integer denoting group number this option belongs to
#                 "requires":      string of 0 or more other options this option requires to work, each separated by a ','
#                 "incompatiable": string of 0 or more other options this option is incompatible with, each separated by a ','
#                 "preamble":      string describing option for preamble section (beginning of output from script)
#                 "help":          string describing option for help section (printed if -h used)
#                 "setby":         '1' if option set by user, else 'undef'
#                 "value":         value for option, can be undef if default is undef
#
# opt_order_A: array of options in the order they should be processed
# 
# opt_group_desc_H: key: group number (integer), value: description of group for help output
my %opt_HH = ();      
my @opt_order_A = (); 
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
$opt_group_desc_H{"1"} = "basic options";
#     option            type       default               group   requires incompat    preamble-output                          help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                   "display this help",                                  \%opt_HH, \@opt_order_A);
opt_Add("-c",           "boolean", 0,                        1,    undef, undef,      "genome is circular",                    "genome is circular",                                 \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",           "force; if dir <reference accession> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                            "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--dirout",     "string",  undef,                    1,    undef, undef,      "output directory specified as <s>",     "specify output directory as <s>, not <ref accession>", \%opt_HH, \@opt_order_A);
opt_Add("--matpept",    "string",  undef,                    1,    undef, undef,      "using pre-specified mat_peptide info",  "read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships", \%opt_HH, \@opt_order_A);
opt_Add("--nomatpept",  "boolean", 0,                        1,    undef,"--matpept", "ignore mat_peptide annotation",         "ignore mat_peptide information in reference annotation", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",    "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"2"} = "options affecting calibration of models";
#       option       type        default                group  requires incompat          preamble-output                                              help-output    
opt_Add("--slow",      "boolean", 0,                     2,    undef, undef,              "running cmcalibrate in slow mode",                           "use default cmcalibrate parameters, not parameters optimized for speed", \%opt_HH, \@opt_order_A);
opt_Add("--local",     "boolean", 0,                     2,    undef, undef,              "running cmcalibrate on local machine",                       "run cmcalibrate locally, do not submit calibration jobs for each CM to the compute farm", \%opt_HH, \@opt_order_A);
opt_Add("--wait",      "integer", 1800,                  2,    undef,"--local,--nosubmit","allow <n> minutes for cmcalibrate jobs on farm",             "allow <n> wall-clock minutes for cmscan jobs on farm to finish, including queueing time", \%opt_HH, \@opt_order_A);
opt_Add("--nosubmit",  "boolean", 0,                     2,    undef,"--local",           "do not submit cmcalibrate jobs to farm, run later",          "do not submit cmcalibrate jobs to farm, run later with qsub script", \%opt_HH, \@opt_order_A);
opt_Add("--errcheck",  "boolean", 0,                     2,    undef,"--local",           "consider any stderr output as indicating a job failure",     "consider any stderr output as indicating a job failure", \%opt_HH, \@opt_order_A);
opt_Add("--rammult",   "boolean", 0,                     2,    undef, undef,              "for all models, multiply RAM Gb by ncpu for mem_free",       "for all models, multiply RAM Gb by ncpu for mem_free", \%opt_HH, \@opt_order_A);
opt_Add("--bigthresh", "integer", "3000",                2,    undef, undef,              "set minimum length for a big model to <n>",                  "set minimum length for a big model to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--bigram",    "integer", "8",                   2,    undef, undef,              "for big models, set Gb RAM per core for calibration to <n>", "for big models, set Gb RAM per core for calibration to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--biglen",    "real",    "0.16",                2,    undef, undef,              "for big models, set length to search in Mb as <x>",          "for big models, set cmcalibrate length to search in Mb as <x>", \%opt_HH, \@opt_order_A);
opt_Add("--bigncpu",   "integer", "8",                   2,    undef, undef,              "for big models, set number of CPUs for calibration to <n>",  "for big models, set number of CPUs for calibration to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--bigtailp",  "real",    "0.30",                2,    undef, undef,              "for big models, set --tailp cmcalibrate parameter as <x>",   "for big models, set --tailp cmcalibrate parameter as <x>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"3"} = "optional output files";
#       option       type       default                group  requires incompat  preamble-output                          help-output    
opt_Add("--mdlinfo",    "boolean", 0,                        3,    undef, undef, "output internal model information",     "create file with internal model information",   \%opt_HH, \@opt_order_A);
opt_Add("--ftrinfo",    "boolean", 0,                        3,    undef, undef, "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"4"} = "options for skipping stages and using files from an earlier, identical run, primarily useful for debugging";
#     option               type       default               group   requires    incompat                  preamble-output                                            help-output    
opt_Add("--skipedirect",   "boolean", 0,                       4,   undef,      undef,                    "skip the edirect steps, use existing results",           "skip the edirect steps, use data from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipfetch",     "boolean", 0,                       4,   undef,      undef,                    "skip the sequence fetching steps, use existing results", "skip the sequence fetching steps, use files from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipbuild",     "boolean", 0,                       4,   undef,      undef,                    "skip the build/calibrate steps",                         "skip the model building/calibrating, requires --mdlinfo and/or --ftrinfo", \%opt_HH, \@opt_order_A);


# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_build.pl [-options] <reference accession>\n";
my $synopsis = "dnaorg_build.pl :: build homology models for features of a reference sequence";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'c'            => \$GetOptions_H{"-c"},
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'dirout=s'     => \$GetOptions_H{"--dirout"},
                'matpept=s'    => \$GetOptions_H{"--matpept"},
                'nomatpept'    => \$GetOptions_H{"--nomatpept"},
                'keep'         => \$GetOptions_H{"--keep"},
# calibration related options
                'slow'         => \$GetOptions_H{"--slow"},
                'local'        => \$GetOptions_H{"--local"},
                'wait=s'       => \$GetOptions_H{"--wait"},
                'nosubmit'     => \$GetOptions_H{"--nosubmit"},
                'errcheck'     => \$GetOptions_H{"--errcheck"},  
                'rammult'      => \$GetOptions_H{"--rammult"},
                'bigthresh=s'  => \$GetOptions_H{"--bigthresh"},
                'bigram=s'     => \$GetOptions_H{"--bigram"},
                'biglen=s'     => \$GetOptions_H{"--biglen"},
                'bigncpu=s'    => \$GetOptions_H{"--bigncpu"},
                'bigtailp=s'   => \$GetOptions_H{"--bigtailp"},
# options for skipping stages, using earlier results
                'skipedirect'  => \$GetOptions_H{"--skipedirect"},
                'skipfetch'    => \$GetOptions_H{"--skipfetch"},
                'skipbuild'    => \$GetOptions_H{"--skipbuild"},
# optional output files
                'mdlinfo'      => \$GetOptions_H{"--mdlinfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"});

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.14";
my $releasedate   = "Sep 2016";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 1) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do dnaorg_build.pl -h\n\n";
  exit(1);
}
my ($ref_accn) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# do checks that are too sophisticated for epn-options.pm
if((opt_Get("--skipbuild", \%opt_HH)) && 
   (! (opt_Get("--mdlinfo", \%opt_HH) || opt_Get("--ftrinfo", \%opt_HH)))) { 
  die "ERROR, --skipbuild requires one or both of --mdlinfo or --ftrinfo"; 
}

my $dir        = opt_Get("--dirout", \%opt_HH);          # this will be undefined unless -d set on cmdline
my $do_matpept = opt_IsOn("--matpept", \%opt_HH);

#############################
# create the output directory
#############################
my $cmd;              # a command to run with runCommand()
my @early_cmd_A = (); # array of commands we run before our log file is opened
# check if the $dir exists, and that it contains the files we need
# check if our output dir $symbol exists
if(! defined $dir) { 
  $dir = $ref_accn;
}
else { 
  if($dir !~ m/\/$/) { $dir =~ s/\/$//; } # remove final '/' if it exists
}
if(-d $dir) { 
  $cmd = "rm -rf $dir";
  if(opt_Get("-f", \%opt_HH)) { runCommand($cmd, opt_Get("-v", \%opt_HH), undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR directory named $dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $dir) { 
  $cmd = "rm $dir";
  if(opt_Get("-f", \%opt_HH)) { runCommand($cmd, opt_Get("-v", \%opt_HH), undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR a file named $dir already exists. Remove it, or use -f to overwrite it."; }
}

# create the dir
$cmd = "mkdir $dir";
runCommand($cmd, opt_Get("-v", \%opt_HH), undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root = $dir . "/" . $dir_tail . ".dnaorg_build";

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("reference accession");
my @arg_A      = ($ref_accn);
outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble(*STDOUT, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# open the log and command files:
# set output file names and file handles, and open those file handles
my %ofile_info_HH = ();  # hash of information on output files we created,
                         # 1D keys: 
                         #  "fullpath":  full path to the file
                         #  "nodirpath": file name, full path minus all directories
                         #  "desc":      short description of the file
                         #  "FH":        file handle to output to for this file, maybe undef
                         # 2D keys:
                         #  "log": log file of what's output to stdout
                         #  "cmd": command file with list of all commands executed

# open the log and command files 
openAndAddFileToOutputInfo(\%ofile_info_HH, "log", $out_root . ".log", 1, "Output printed to screen");
openAndAddFileToOutputInfo(\%ofile_info_HH, "cmd", $out_root . ".cmd", 1, "List of executed commands");
openAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".list", 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# open optional output files
if(opt_Get("--mdlinfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "mdlinfo", $out_root . ".mdlinfo", 1, "Model information (created due to --mdlinfo)");
}
if(opt_Get("--ftrinfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "ftrinfo", $out_root . ".ftrinfo", 1, "Feature information (created due to --ftrinfo)");
}

# now we have the log file open, output the banner there too
outputBanner($log_FH, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

########################################
# parse the optional input files, if nec
########################################
# -matpept <f>
my @cds2pmatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of primary matpept indices that comprise this CDS
my @cds2amatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of all     matpept indices that comprise this CDS
if($do_matpept) { 
  parseMatPeptSpecFile(opt_Get("--matpept", \%opt_HH), \@cds2pmatpept_AA, \@cds2amatpept_AA, $ofile_info_HH{"FH"});
}

###################################################
# make sure the required executables are executable
###################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmbuild"}       = $inf_exec_dir . "cmbuild";
$execs_H{"cmcalibrate"}   = $inf_exec_dir . "cmcalibrate";
$execs_H{"cmfetch"}       = $inf_exec_dir . "cmfetch";
$execs_H{"cmpress"}       = $inf_exec_dir . "cmpress";
$execs_H{"esl-afetch"}    = $esl_exec_dir . "esl-afetch";
$execs_H{"esl-reformat"}  = $esl_exec_dir . "esl-reformat";
$execs_H{"esl_fetch_cds"} = $esl_fetch_cds;
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################################################################
# Step 0. Output the consopts file that dnaorg_annotate.pl will use to 
#         make sure options used are consistent between dnaorg_build.pl and 
#         dnaorg_annotate.pl 
###########################################################################
my $progress_w = 80; # the width of the left hand column in our progress output, hard-coded
my $start_secs = outputProgressPrior("Outputting information on options used for future use with dnaorg_annotate.pl", $progress_w, $log_FH, *STDOUT);
output_consopts_file($out_root . ".consopts", \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###########################################################################
# Step 1. Gather and process information on reference genome using Edirect.
###########################################################################
$start_secs = outputProgressPrior("Gathering information on reference using edirect", $progress_w, $log_FH, *STDOUT);

my %cds_tbl_HHA = ();   # CDS data from .cds.tbl file, hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %mp_tbl_HHA = ();    # mat_peptide data from .matpept.tbl file, hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %seq_info_HA = ();   # hash of arrays, avlues are arrays [0..$nseq-1];
                        # 1st dim keys are "seq_name", "accn_name", "seq_len", "accn_len".
                        # $seq_info_HA{"accn_name"}[0] is our reference accession
@{$seq_info_HA{"accn_name"}} = ($ref_accn);

# Call the wrapper function that does the following:
#  1) creates the edirect .mat_peptide file, if necessary
#  2) creates the edirect .ftable file
#  3) creates the length file
#  4) parses the edirect .mat_peptide file, if necessary
#  5) parses the edirect .ftable file
#  6) parses the length file
wrapperGetInfoUsingEdirect(undef, $ref_accn, $out_root, \%cds_tbl_HHA, \%mp_tbl_HHA, \%seq_info_HA, \%ofile_info_HH, 
                           \%opt_HH, $ofile_info_HH{"FH"}); # 1st argument is undef because we are only getting info for $ref_accn

if($do_matpept) {  
  # validate the CDS:mat_peptide relationships that we read from the $matpept input file
  matpeptValidateCdsRelationships(\@cds2pmatpept_AA, \%{$cds_tbl_HHA{$ref_accn}}, \%{$mp_tbl_HHA{$ref_accn}}, opt_Get("-c", \%opt_HH), $seq_info_HA{"accn_len"}[0], $ofile_info_HH{"FH"});
}
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################
# Step 2. Fetch and process the reference genome sequence
##########################################################
$start_secs = outputProgressPrior("Fetching and processing the reference genome", $progress_w, $log_FH, *STDOUT);
my %mdl_info_HA = ();          # hash of arrays, values are arrays [0..$nmdl-1];
                               # see dnaorg.pm::validateModelInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAllSequencesAndProcessReferenceSequence()
my %ftr_info_HA = ();          # hash of arrays, values are arrays [0..$nftr-1], 
                               # see dnaorg.pm::validateFeatureInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAllSequencesAndProcessReferenceSequence()
my $sqfile = undef;            # pointer to the Bio::Easel::SqFile object we'll open in wrapperFetchAllSequencesAndProcessReferenceSequence()


# Call the wrapper function that does the following:
#   1) fetches the sequences listed in @{$seq_info_HAR->{"accn_name"} into a fasta file 
#      and indexes that fasta file, the reference sequence is $seq_info_HAR->{"accn_name"}[0].
#   2) determines information for each feature (strand, length, coordinates, product) in the reference sequence
#   3) determines type of each reference sequence feature ('cds-mp', 'cds-notmp', or 'mp')
#   4) fetches the reference sequence feature and populates information on the models and features
wrapperFetchAllSequencesAndProcessReferenceSequence(\%execs_H, \$sqfile, $out_root, $out_root, # yes, $out_root is passed in twice, on purpose
                                                    undef, undef, undef, undef,  # 4 variables used only if --infasta enabled in dnaorg_annotate.pl (irrelevant here)
                                                    \%cds_tbl_HHA,
                                                    ($do_matpept) ? \%mp_tbl_HHA      : undef, 
                                                    ($do_matpept) ? \@cds2pmatpept_AA : undef, 
                                                    ($do_matpept) ? \@cds2amatpept_AA : undef, 
                                                    \%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA,
                                                    \%opt_HH, \%ofile_info_HH);

# verify our model and feature info hashes are complete, 
# if validateFeatureInfoHashIsComplete() fails then the program will exit with an error message
my $nftr = validateFeatureInfoHashIsComplete(\%ftr_info_HA, undef, $ofile_info_HH{"FH"}); # nftr: number of features
my $nmdl = validateModelInfoHashIsComplete  (\%mdl_info_HA, undef, $ofile_info_HH{"FH"}); # nmdl: number of homology models

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#dumpInfoHashOfArrays("ftr_info", 0, \%ftr_info_HA, *STDOUT);

####################################
# Step 3. Build and calibrate models 
####################################
my $do_local      = opt_Get("--local", \%opt_HH); # are we running calibration locally
my $do_farm_now   = 0; # set to '1' if we are submitting to farm now and waiting for jobs to finish
my $do_farm_later = 0; # set to '1' if we are creating script to submit to farm later
if(! $do_local) { 
  $do_farm_now   = opt_Get("--nosubmit", \%opt_HH) ? 0 : 1; 
  $do_farm_later = opt_Get("--nosubmit", \%opt_HH) ? 1 : 0; 
}

if(! opt_Get("--skipbuild", \%opt_HH)) { 
  my $build_str = opt_Get("--local", \%opt_HH) ? "Building models locally" : "Submitting jobs to build models to compute farm and waiting for them to finish";
  $start_secs = outputProgressPrior($build_str, $progress_w, $log_FH, *STDOUT);
  run_cmbuild(\%execs_H, $out_root, $ofile_info_HH{"fullpath"}{"refstk"}, \@{$mdl_info_HA{"cmname"}}, \%opt_HH, \%ofile_info_HH);
  if(! opt_Get("--local", \%opt_HH)) { 
    outputString($log_FH, 1, "# ");
  }
  outputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);

  # calibrate models
  my $calibrate_str = "";
  if(opt_Get("--local", \%opt_HH)) { 
    $calibrate_str = "Calibrating models locally";
  }
  elsif(opt_Get("--nosubmit", \%opt_HH)) { 
    $calibrate_str  = "Creating script for calibrating models later (due to --nosubmit)";
  }
  else { 
    $calibrate_str  = "Submitting jobs to calibrate models to compute farm and waiting for them to finish";
  }
  $start_secs = outputProgressPrior($calibrate_str, $progress_w, $log_FH, *STDOUT);
  run_cmcalibrate_and_cmpress(\%execs_H, $out_root, $nmdl, \%opt_HH, \%ofile_info_HH);
  if((! opt_Get("--local", \%opt_HH)) && (! opt_Get("--nosubmit", \%opt_HH))) { 
    outputString($log_FH, 1, "# ");
  }
  for(my $i = 0; $i < $nmdl; $i++) { 
    addClosedFileToOutputInfo(\%ofile_info_HH, "cm$i", "$out_root.$i.cm", 1, 
                              sprintf("CM file #%d, %s%s", $i+1, $mdl_info_HA{"out_tiny"}[$i], 
                                      (opt_Get("--nosubmit", \%opt_HH)) ? " (needs to be calibrated later by running \"sh $out_root.cm.qsub\")" : ""));
  }
  if(opt_Get("--nosubmit", \%opt_HH)) { 
    addClosedFileToOutputInfo(\%ofile_info_HH, "qsub", "$out_root.cm.qsub", 1, "Shell script to submit cmcalibrate commands with (not executed yet, due to --nosubmit, execute it later)");
    addClosedFileToOutputInfo(\%ofile_info_HH, "postscript", "$out_root.cm.run_when_jobs_are_finished.sh", 1, "Shell script to run after cmcalibrate jobs submitted with $out_root.cm.qsub are **finished running**. This script prepares the CM files for use with dnaorg_annotate.pl.");
  }
  outputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);
}

##########
# Conclude
##########
# output optional output files
if(exists $ofile_info_HH{"FH"}{"mdlinfo"}) { 
  dumpInfoHashOfArrays("Model information (%mdl_info_HA)", 0, \%mdl_info_HA, $ofile_info_HH{"FH"}{"mdlinfo"});
}
if(exists $ofile_info_HH{"FH"}{"ftrinfo"}) { 
  dumpInfoHashOfArrays("Feature information (%ftr_info_HA)", 0, \%ftr_info_HA, $ofile_info_HH{"FH"}{"ftrinfo"});
}

# a quick note to the user about what to do next
if(! opt_Get("--skipbuild", \%opt_HH)) { 
  outputString($log_FH, 1, sprintf("#\n"));
  if(! opt_Get("--nosubmit", \%opt_HH)) { 
    outputString($log_FH, 1, "# You can now use dnaorg_annotate.pl to annotate genomes with the models that\n");
    outputString($log_FH, 1, "# you've created here.\n");
  }
  else { 
    outputString($log_FH, 1, "# The models you've created have not yet been calibrated so you cannot use\n");
    outputString($log_FH, 1, "# dnaorg_annotate.pl with them yet. First, you need to execute $out_root.cm.qsub\n");
    outputString($log_FH, 1, "# in the directory $dir. That will submit jobs to the compute farm. Wait until those\n");
    outputString($log_FH, 1, "# jobs have all finished running (monitor with qstat), and then execute\n");
    outputString($log_FH, 1, "# $out_root.cm.run_when_jobs_are_finished.sh in the directory $dir. That script\n");
    outputString($log_FH, 1, "# will prepare the calibrated CM files for use with dnaorg_annotate.pl. When that\n");
    outputString($log_FH, 1, "# script finishes you can then use dnaorg_annotate.pl with the models.\n");
  }    
  outputString($log_FH, 1, sprintf("#\n"));
}

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir, \%ofile_info_HH);
exit 0;


#################################################################
# Subroutine: output_consopts_file()
# Incept:     EPN, Fri May 27 14:20:28 2016
#
# Purpose:   Output a simple .consopts file that includes 
#            information on options that must be kept consistent
#            between dnaorg_build.pl and dnaorg_annotate.pl.
#            Currently this is only "--matpept", "--nomatpept",
#            and "-c".
#
# Arguments:
#  $consopts_file:     name of the file to create
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:    REF to output file hash
# 
# Returns:  void
# 
# Dies: If $consopts_file doesn't exist, or we can't parse it.
#
#       If an option enabled in dnaorg_build.pl that needs to 
#       be consistently used in dnaorg_annotate.pl is not, or
#       vice versa.
#
#       If an option enabled in dnaorg_build.pl that takes a file
#       as input has a different checksum for that file than 
#       
#################################################################
sub output_consopts_file { 
  my $sub_name = "output_consopts_files";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($consopts_file, $opt_HHR, $ofile_info_HHR) = @_;

  open(OUT, ">", $consopts_file) || fileOpenFailure($consopts_file, $sub_name, $!, "writing", $ofile_info_HHR->{"FH"});

  my $printed_any_options = 0;
  if(opt_Get("-c", $opt_HHR)) { 
    print OUT ("-c\n");
    $printed_any_options = 1;
  }
  if(opt_Get("--nomatpept", $opt_HHR)) { 
    print OUT ("--nomatpept\n");
    $printed_any_options = 1;
  }
  if(opt_IsUsed("--matpept", $opt_HHR)) { 
    my $matpept_file = opt_Get("--matpept", $opt_HHR);
    printf OUT ("--matpept %s %s\n", 
                $matpept_file, 
                md5ChecksumOfFile($matpept_file, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"}));
    $printed_any_options = 1;
  }
  
  if(! $printed_any_options) { 
    print OUT "none\n";
  }
  close(OUT);
  
  addClosedFileToOutputInfo($ofile_info_HHR, "consopts", "$consopts_file", 1, "File with list of options that must be kept consistent between dnaorg_build.pl and dnaorg_annotate.pl runs");
  return;
}

#################################################################
# Subroutine: run_cmbuild()
# Incept:     EPN, Wed Aug 31 13:30:02 2016
#
# Purpose:   Run cmbuild either locally or submit jobs to farm
#            and wait for them to finish.
#
# Arguments:
#   $execs_HR:       reference to hash with infernal executables, 
#                    e.g. $execs_HR->{"cmbuild"} is path to cmbuild, PRE-FILLED
#   $out_root:       string for naming output files
#   $stk_file_AR:    reference to array of stockholm files, one per model
#   $indi_name_AR:   ref to array of individual model names, we only use this if 
#                    $do_calib_local is 0 or undef, PRE-FILLED
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $ofile_info_HHR: REF to the 2D hash of output file information
# 
# Returns:  void
# 
# Dies:  if any of the files in @{$stk_file_AR} do not exist or are empty
#        if we can't determine consensus length from the model file
#
#################################################################
sub run_cmbuild { 
  my $sub_name = "run_cmbuild";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $out_root, $stk_file, $indi_name_AR, $opt_HHR, $ofile_info_HHR) = @_;

  # for convenience
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  if(! -s $stk_file)  { DNAORG_FAIL("ERROR in $sub_name, $stk_file file does not exist or is empty", 1, $FH_HR); }

  my $cmbuild = $execs_HR->{"cmbuild"};
  my $afetch  = $execs_HR->{"esl-afetch"};

  my $abs_out_root = getcwd() . "/" . $out_root;
  my $do_local = opt_Get("--local", $opt_HHR); # should we run locally (instead of on farm)?
  my ($cmbuild_opts, $cmbuild_cmd); # options and command for cmbuild

  my $nmodel = scalar(@{$indi_name_AR});

  my $out_tail    = $out_root;
  $out_tail       =~ s/^.+\///;

  my @out_A = ();
  my @err_A = ();
  for(my $i = 0; $i < $nmodel; $i++) { 
    $cmbuild_opts = "-F --informat stockholm";
    $cmbuild_cmd  = "$afetch $stk_file $indi_name_AR->[$i] | $cmbuild $cmbuild_opts $out_root.$i.cm - > $out_root.$i.cmbuild";
    if($do_local) { 
      runCommand($cmbuild_cmd, opt_Get("-v", $opt_HHR), $FH_HR);
    }
    else { # submit to farm
      my $jobname = "b." . $out_tail . $i;
      my $errfile = $abs_out_root . ".b." . $i . ".err";
      my $outfile = "$out_root.$i.cmbuild";
      my $farm_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n -l h_rt=288000,h_vmem=8G,mem_free=8G,reserve_mem=8G " . "\"" . $cmbuild_cmd . "\" > /dev/null\n";
      push(@out_A, $outfile);
      push(@err_A, $errfile);
      runCommand($farm_cmd, opt_Get("-v", $opt_HHR), $FH_HR);
    }
  }
  if(! $do_local) { # wait for farm jobs to finish
    my $njobs_finished = waitForFarmJobsToFinish(\@out_A, \@err_A, "# CPU time", opt_Get("--wait", $opt_HHR), opt_Get("--errcheck", $opt_HHR), $FH_HR);
    if($njobs_finished != $nmodel) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name only $njobs_finished of the $nmodel are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", $opt_HHR)), 1, $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine: run_cmcalibrate_and_cmpress()
# Incept:     EPN, Wed Aug 31 14:27:48 2016
# 
# Purpose:    Calibrate a set of CM files.
#
# Arguments:
#   $execs_HR:       reference to hash with infernal executables, 
#                    e.g. $execs_HR->{"cmcalibrate"} is path to cmcalibrate, PRE-FILLED
#   $out_root:       string for naming output files
#   $nmodel:         number of models
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $ofile_info_HHR: REF to the 2D hash of output file information
#                    
# Returns:    void
#
# Dies:       if $stk_file does not exist or is empty
#             if we can't determine consensus length from the model file
#################################################################
sub run_cmcalibrate_and_cmpress { 
  my $sub_name = "run_cmcalibrate_and_cmpress";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $out_root, $nmodel, $opt_HHR, $ofile_info_HHR) = @_;

  # for convenience
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $abs_out_root = getcwd() . "/" . $out_root;

  my $do_calib_slow  = opt_Get("--slow", $opt_HHR);  # should we run in 'slow' mode (instead of fast mode)?
  my $do_calib_local = opt_Get("--local", $opt_HHR); # should we run locally (instead of on farm)?
  my $do_calib_later = opt_Get("--nosubmit", $opt_HHR); # are we running calibrations later (creating shell script to do it)

  my ($df_cmcalibrate_opts,  $df_cmcalibrate_cmd);    # options and command for default cmcalibrate
  my ($big_cmcalibrate_opts, $big_cmcalibrate_cmd);   # options and command for big-model cmcalibrate
  my ($cmpress_opts,         $cmpress_cmd);           # options and command for cmpress
  my $df_qsub_opts;  # default   options for qsub for cmcalibrate job
  my $big_qsub_opts; # big-model options for qsub for cmcalibrate job

  my $cmcalibrate = $execs_HR->{"cmcalibrate"};
  my $cmpress     = $execs_HR->{"cmpress"};

  # set up cmcalibrate options for two scenarios: default and 'big' model
  my $df_ncpu         = 4;
  my $df_gb_per_core  = 8;
  my $df_gb_tot       = opt_Get("--rammult", $opt_HHR) ? $df_gb_per_core * $df_ncpu : $df_gb_per_core;
  my $df_time_and_mem_req = sprintf("-l h_rt=576000,h_vmem=%sG,mem_free=%sG,reserve_mem=%sG", $df_gb_tot, $df_gb_tot, $df_gb_per_core);
  my $big_ncpu        = opt_Get("--bigncpu", $opt_HHR);
  my $big_gb_per_core = opt_Get("--bigram", $opt_HHR);
  my $big_gb_tot      = opt_Get("--rammult", $opt_HHR) ? $big_gb_per_core * $big_ncpu : $big_gb_per_core;

  my $big_thresh      = opt_Get("--bigthresh", $opt_HHR);
  my $big_time_and_mem_req = sprintf("-l h_rt=576000,h_vmem=%sG,mem_free=%sG,reserve_mem=%sG", $big_gb_tot, $big_gb_tot, $big_gb_per_core);

  $df_cmcalibrate_opts  = " --cpu $df_ncpu ";
  $big_cmcalibrate_opts = " --cpu $big_ncpu ";
  if(! $do_calib_slow) { 
    $df_cmcalibrate_opts  .= " -L 0.04 ";
  }
  $big_cmcalibrate_opts .= " -L " . opt_Get("--biglen", $opt_HHR) . " --tailp " . opt_Get("--bigtailp", $opt_HHR) . " ";

  my $df_cmcalibrate_cmd_root   = "$cmcalibrate $df_cmcalibrate_opts";
  my $big_cmcalibrate_cmd_root  = "$cmcalibrate $big_cmcalibrate_opts";
  
  # run cmcalibrate either on farm or locally, one job for each CM file
  # split up model file into individual CM files, then submit a job to calibrate each one, and exit. 
  # the qsub commands will be submitted by executing a shell script with all of them
  # unless --nosubmit option is enabled, in which case we'll just create the file
  my $farm_cmd_file     = $out_root . ".cm.qsub";
  my $postfarm_cmd_file = $out_root . ".cm.run_when_jobs_are_finished.sh";
  my @tmp_out_file_A = (); # the array of cmcalibrate output files
  my @tmp_err_file_A = (); # the array of error files 
  my $nfarmjobs = 0; # number of jobs submitted to farm
  if(! $do_calib_local) { 
    open(FARMOUT, ">", $farm_cmd_file)     || fileOpenFailure($farm_cmd_file, $sub_name, $!, "writing", $FH_HR);
    open(POSTOUT, ">", $postfarm_cmd_file) || fileOpenFailure($postfarm_cmd_file, $sub_name, $!, "writing", $FH_HR);
  }
  for(my $i = 0; $i < $nmodel; $i++) { 
    my $abs_model_file = "$abs_out_root.$i.cm";
    my $model_file     = "$out_root.$i.cm";
    if(! -s $abs_model_file) { 
      DNAORG_FAIL("ERROR in $sub_name, $abs_model_file does not exist."); 
    }
    my $out_tail    = $out_root;
    $out_tail       =~ s/^.+\///;
    my $jobname     = "c." . $out_tail . $i;
    my $errfile     = $abs_out_root . "." . $i . ".err";
    
    # determine length of the model, if >= opt_HHR->bigthresh, use --tailp, and require double memory (16Gb for 4 threads instead of 8Gb)
    my $is_big = 0;
    my $clen = `grep ^CLEN $model_file`;
    chomp $clen;
    if($clen =~ /^CLEN\s+(\d+)/) { 
      $clen = $1;
      if($clen >= $big_thresh) { $is_big = 1; }
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, couldn't determine consensus length in CM file $model_file, got $clen", 1, $FH_HR);
    }
    
    my $actual_command = $df_cmcalibrate_cmd_root . " $abs_model_file > $abs_out_root.$i.cmcalibrate";
    if($is_big) { 
      $actual_command = $big_cmcalibrate_cmd_root . "$abs_model_file > $abs_out_root.$i.cmcalibrate";
    }
    if($do_calib_local) { 
      runCommand($actual_command, opt_Get("-v", $opt_HHR), $FH_HR);
    }
    else { # ! local, run job on farm
      my $farm_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n $df_time_and_mem_req -pe multicore $df_ncpu -R y " . "\"" . $actual_command . "\"" . " > /dev/null\n";
      if($is_big) { # rewrite it
        $farm_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n $big_time_and_mem_req -pe multicore $big_ncpu -R y " . "\"" . $actual_command . "\"" . " > /dev/null\n";
      }
      push(@tmp_out_file_A, "$abs_out_root.$i.cmcalibrate");
      push(@tmp_err_file_A, $errfile);
      $nfarmjobs++;
      print FARMOUT $farm_cmd;
    }
  } # end of 'for(my $i = 0; $i < $nmodel; $i++)' 

  if(! $do_calib_local) { 
    # if we didn't run jobs locally...
    close(FARMOUT);
  
    my $cmfile; # name of a CM file
    my $cmd; # a command
    if(! $do_calib_later) { 
      # run the cmcalibrate commands
      runCommand("sh $farm_cmd_file", opt_Get("-v", $opt_HHR), $FH_HR);
      # and wait for all jobs to finish
      my $njobs_finished = waitForFarmJobsToFinish(\@tmp_out_file_A, \@tmp_err_file_A, "[ok]", opt_Get("--wait", $opt_HHR), opt_Get("--errcheck", $opt_HHR), $FH_HR);
      if($njobs_finished != $nfarmjobs) { 
        DNAORG_FAIL(sprintf("ERROR in main() only $njobs_finished of the $nfarmjobs are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", $opt_HHR)), 1, $FH_HR);
      }
      # now all jobs are finished, now we need to press each file 
    }      

    for(my $i = 0; $i < $nmodel; $i++) { 
      $cmfile = "$out_root.$i.cm";
      $cmd = pressCmDb($cmfile, $cmpress, 1, (!$do_calib_later), $opt_HHR, $ofile_info_HHR);
      print POSTOUT $cmd . "\n";
    }
  }
    
  return;
}
