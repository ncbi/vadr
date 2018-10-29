#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Mon Feb  1 15:07:43 2016 [dnaorg_build.pl split off from dnaorg_annotate_genomes.pl]
# LES, Mon Jul 25 09:25    2016 [dnaorg_refseq_assign.pl split off from dnaorg_build.pl]
# EPN, Tue Nov 28 10:44:54 2017 [dnaorg_classify.pl created, renamed from dnaorg_refseq_assign.pl]
#
# Some code in here is also adopted from dnaorg.pm

use strict;
use warnings;
#use diagnostics;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

#use Data::Dumper;
#use File::Slurp;

require "dnaorg.pm"; 
require "epn-options.pm";

#######################################################################################
# What this script does: 
#
# Preliminaries: 
#   - process options
#   - create the output directory
#   - output program banner and open output files
#   - parse the input files
#
# If --onlybuild: Create an HMM library of RefSeqs, and then exit.
# 
# Else (--inlist of --infasta): 
#
# Step 1. Run nhmmscan for all accessions using HMM library from previous --onlybuild run.
#
# Step 2. Parse the nhmmscan output to determine the proper RefSeq for each accession
#
# Step 3. Generate ntlists and possibly fasta files for each RefSeq
#
#######################################################################################

# first, determine the paths to all modules, scripts and executables that we'll need
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
my $hmmer_exec_dir    = $dnaorgdir . "/hmmer-3.1b2/src/";
my $esl_exec_dir      = $dnaorgdir . "/infernal-1.1.2/easel/miniapps/";
my $bioeasel_exec_dir = $dnaorgdir . "/Bio-Easel/scripts/";
my $dnaorg_exec_dir   = $dnaorgdir . "/dnaorg_scripts/";

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
#     option            type       default               group   requires incompat    preamble-output                          help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                   "display this help",                                  \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"1"} = "REQUIRED options";
#     option            type       default               group   requires incompat    preamble-output                          help-output    
opt_Add("--dirout",     "string", undef,                     1,    undef, undef,       "REQUIRED name for output directory to create is <s>",  "REQUIRED: name for output directory to create is <s>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"2"} = "options for selecting mode to run in: build-mode (--onlybuild) or classify-mode (--inlist OR --infasta)";
#     option            type       default               group   requires       incompat                 preamble-output                                     help-output    
opt_Add("--onlybuild",  "string", undef,                    2,    undef,        "--infasta,--inlist",    "build an HMM library for seqs in <s>, then exit",  "build an HMM library for sequences listed in <s>, then exit", \%opt_HH, \@opt_order_A);
opt_Add("--inlist",     "string", undef,                    2,    "--dirbuild", "--onlybuild,--infasta", "list of sequence accessions to classify is <s>",   "list of sequence accessions to classify is <s>, requires --dirbuild", \%opt_HH, \@opt_order_A);
opt_Add("--infasta",    "string", undef,                    2,    "--dirbuild", "--onlybuild,--inlist",  "fasta file with sequences to classify is <s>",     "fasta file with sequences to classify is <s>, requires --dirbuild", \%opt_HH, \@opt_order_A);
opt_Add("--dirbuild",   "string", undef,                    2,    undef,        "--onlybuild",           "specify directory with HMM library is <s>",        "specify directory with HMM library is <s>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"3"} = "basic options";
#     option            type       default               group   requires incompat    preamble-output                                                help-output    
opt_Add("-f",           "boolean", 0,                        3,    undef, undef,      "forcing directory overwrite",                                 "force; if dir <reference accession> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        3,    undef, undef,      "be verbose",                                                  "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                        3,    undef, undef,      "leaving intermediate files on disk",                          "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);
opt_Add("--nkb",        "integer", 5,                        3,    undef,"--local",   "number of KB of sequence for each nhmmscan farm job",         "set target number of KB of sequences for each nhmscan farm job to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--maxnjobs",   "integer", 500,                      3,    undef,"--local",   "maximum allowed number of jobs for compute farm",             "set max number of jobs to submit to compute farm to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--wait",       "integer", 500,                      3,    undef,"--local",   "allow <n> minutes for nhmmscan jobs on farm",                 "allow <n> wall-clock minutes for nhmmscan jobs on farm to finish, including queueing time", \%opt_HH, \@opt_order_A);
opt_Add("--local",      "boolean", 0,                        3,    undef, undef,      "run nhmmscan locally instead of on farm",                     "run nhmmscan locally instead of on farm", \%opt_HH, \@opt_order_A);
opt_Add("--errcheck",   "boolean", 0,                        3,    undef,"--local",   "consider any farm stderr output as indicating a job failure", "consider any farm stderr output as indicating a job failure", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"4"} = "options for controlling what unexpected features cause sequences to PASS/FAIL";
#     option                type       default            group   requires incompat         preamble-output                                              help-output    
opt_Add("--allpass"  ,     "boolean",  0,                    4,   undef,   undef,           "all assigned seqs pass, only non-assigned seqs fail",        "all assigned sequences pass, only non-assigned seqs fail", \%opt_HH, \@opt_order_A);
opt_Add("--lowscpass",     "boolean",  0,                    4,   undef,"--allpass",        "seqs with LowScore     unexpected feature can PASS",         "sequences with LowScore     unexpected feature can PASS",  \%opt_HH, \@opt_order_A);
opt_Add("--vlowscpass",    "boolean",  0,                    4,"--lowscpass","--allpass",   "seqs with VeryLowScore unexpected feature can PASS",         "sequences with VeryLowScore unexpected feature can PASS",  \%opt_HH, \@opt_order_A);
opt_Add("--lowdiffpass",   "boolean",  0,                    4,   undef,"--allpass",        "seqs with LowDiff      unexpected feature can PASS",         "sequences with LowDiff      unexpected feature can PASS",  \%opt_HH, \@opt_order_A);
opt_Add("--vlowdiffpass",  "boolean",  0,                    4,"--lowdiffpass","--allpass", "seqs with VeryLowDiff  unexpected feature can PASS",         "sequences with VeryLowDiff  unexpected feature can PASS",  \%opt_HH, \@opt_order_A);
opt_Add("--biaspass",      "boolean",  0,                    4,   undef,"--allpass",        "seqs with HighBias     unexpected feature can PASS",         "sequences with HighBias     unexpected feature can PASS",  \%opt_HH, \@opt_order_A);
opt_Add("--minusfail",     "boolean",  0,                    4,   undef,"--allpass",        "seqs with MinusStrand  unexpected feature FAIL",             "sequences with MinusStrand  unexpected feature FAIL",      \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"5"} = "options for controlling reporting of unexpected features";
#     option                type         default            group   requires incompat     preamble-output                                                    help-output    
opt_Add("--lowscthresh",    "real",    0.3,                   5,   undef,   undef,        "bits per nucleotide threshold for LowScore is <x>",               "bits per nucleotide threshold for LowScore unexpected feature is <x>",         \%opt_HH, \@opt_order_A);
opt_Add("--vlowscthresh",   "real",    0.2,                   5,   undef,   undef,        "bits per nucleotide threshold for VeryLowScore is <x>",           "bits per nucleotide threshold for VeryLowScore unexpected feature is <x>",     \%opt_HH, \@opt_order_A);
opt_Add("--lowdiffthresh",  "real",    0.06,                  5,   undef,   undef,        "bits per nucleotide diff threshold for LowDiff is <x>",           "bits per nucleotide diff threshold for LowDiff unexpected feature is <x>",     \%opt_HH, \@opt_order_A);
opt_Add("--vlowdiffthresh", "real",    0.006,                 5,   undef,   undef,        "bits per nucleotide diff threshold for VeryLowDiff is <x>",       "bits per nucleotide diff threshold for VeryLowDiff unexpected feature is <x>", \%opt_HH, \@opt_order_A);
opt_Add("--biasfract",      "real",    0.25,                  5,   undef,   undef,        "fractional threshold for HighBias is <x>",                        "fractional threshold for HighBias unexpected feature is <x>",                  \%opt_HH, \@opt_order_A);
opt_Add("--lowscminlen",    "integer",  501,                   5,   undef,"--lowscpass",  "minimum seq len for which LowScore causes FAILure is <n>",        "set minimum length for which LowScore causes FAILure to <n>",                  \%opt_HH, \@opt_order_A);
opt_Add("--lowdiffminlen",  "integer", 1001,                  5,   undef,"--lowdiffpass", "minimum seq len for which LowDiff causes FAILure is <n>",         "set minimum length for which LowDiff causes FAILure to <n>",                   \%opt_HH, \@opt_order_A);
opt_Add("--nolowscminlen",  "boolean", 0,                     5,   undef,"--lowscpass",   "no minimum length for which LowScore causes a seq to FAIL",       "no minimum length for which LowScore causes a seq to FAIL",                    \%opt_HH, \@opt_order_A);
opt_Add("--nolowdiffminlen","boolean", 0,                     5,   undef,"--lowdiffpass", "no minimum length for which LowDiff causes a seq to FAIL",        "no minimum length for which LowDiff causes a seq to FAIL",                     \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"6"} = "options for automatically running dnaorg_annotate.pl for classified sequences";
#     option            type       default               group   requires       incompat          preamble-output                help-output    
opt_Add("-A",           "string", undef,                    6,    undef,        "--onlybuild",    "annotate after classifying using build dirs in dir <s>",  "annotate using dnaorg_build.pl build directories in <s> after classifying", \%opt_HH, \@opt_order_A);
opt_Add("--optsA",      "string", undef,                    6,    "-A",         "--onlybuild",    "read dnaorg_annotate.pl options from file <s>",           "read additional dnaorg_annotate.pl options from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--reflistA",   "string", undef,                    6,    "-A",         "--onlybuild",    "only annotate seqs that match to RefSeqs listed in <s>",  "only annotate seqs that match to RefSeqs listed in <s>", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: This script must be run in 1 of 3 modes:\n";
$usage      .= "\nBuild mode (--onlybuild): Build HMM library and exit (no classification).\nExample usage:\n\t";
$usage      .= "dnaorg_classify.pl [-options] --onlybuild <RefSeq list> --dirout <output directory to create with HMM library>\n";
$usage      .= "\nClassify mode given list (--inlist): Use a previously created HMM library to annotate sequence accessions listed in an input file.\nExample usage:\n\t";
$usage      .= "dnaorg_classify.pl [-options] --dirbuild <directory with HMM library to use> --dirout <output directory to create> --inlist <list of accessions to classify>\n";
$usage      .= "\nClassify mode given fasta (--infasta): Use a previously created HMM library to annotate sequences in an input fasta file.\nExample usage:\n\t";
$usage      .= "dnaorg_classify.pl [-options] --dirbuild <directory with HMM library to use> --dirout <output directory to create> --infasta <fasta file with sequences to classify>\n";
my $script_name = "dnaorg_classify.pl";
my $synopsis = "$script_name :: classify sequences using an HMM library of RefSeqs";

my $options_okay = 
    &GetOptions('h'                => \$GetOptions_H{"-h"}, 
# basic options
                'f'                => \$GetOptions_H{"-f"},
                'v'                => \$GetOptions_H{"-v"},
                'keep'             => \$GetOptions_H{"--keep"},
                'dirout=s'         => \$GetOptions_H{"--dirout"},
                'onlybuild=s'      => \$GetOptions_H{"--onlybuild"},
                'dirbuild=s'       => \$GetOptions_H{"--dirbuild"},
                'inlist=s'         => \$GetOptions_H{"--inlist"},
                'infasta=s'        => \$GetOptions_H{"--infasta"},
                'nkb=s'            => \$GetOptions_H{"--nkb"}, 
                'maxnjobs=s'       => \$GetOptions_H{"--maxnjobs"}, 
                'wait=s'           => \$GetOptions_H{"--wait"},
                'local'            => \$GetOptions_H{"--local"}, 
                'errcheck'         => \$GetOptions_H{"--errcheck"},       
                "allpass"          => \$GetOptions_H{"--allpass"},
                "lowscpass"        => \$GetOptions_H{"--lowscpass"},
                "vlowscpass"       => \$GetOptions_H{"--vlowscpass"},
                "lowdiffpass"      => \$GetOptions_H{"--lowdiffpass"},
                "vlowdiffpass"     => \$GetOptions_H{"--vlowdiffpass"},
                "biaspass"         => \$GetOptions_H{"--biaspass"},
                "minusfail"        => \$GetOptions_H{"--minusfail"},
                "lowscthresh=s"    => \$GetOptions_H{"--lowscthresh"},
                "vlowscthresh=s"   => \$GetOptions_H{"--vlowscthresh"},
                "lowdiffthresh=s"  => \$GetOptions_H{"--lowdiffthresh"},
                "vlowdiffthresh=s" => \$GetOptions_H{"--vlowdiffthresh"},
                'biasfract=s'      => \$GetOptions_H{"--biasfract"},  
                'lowscminlen=s'    => \$GetOptions_H{"--lowscminlen"},  
                'lowdiffminlen=s'  => \$GetOptions_H{"--lowdiffminlen"},  
                'nolowscminlen'    => \$GetOptions_H{"--nolowscminlen"},  
                'nolowdiffminlen'  => \$GetOptions_H{"--nolowdiffminlen"},  
                'A=s'              => \$GetOptions_H{"-A"},
                'optsA=s'          => \$GetOptions_H{"--optsA"},
                'reflistA=s'       => \$GetOptions_H{"--reflistA"},
                );

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.36";
my $releasedate   = "Sept 2018";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date, $dnaorgdir);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 0) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do $script_name -h\n\n";
  exit(1);
}

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

my $dir = opt_Get("--dirout", \%opt_HH);          # this will be undefined unless -d set on cmdline
if(! defined $dir) { 
  die "ERROR the --dirout <s> option is required to specify the name of the output directory";
}

# determine which 'mode' we are running in
my $onlybuild_mode = opt_IsUsed("--onlybuild", \%opt_HH);
my $inlist_mode    = opt_IsUsed("--inlist", \%opt_HH);
my $infasta_mode   = opt_IsUsed("--infasta", \%opt_HH);
if(($onlybuild_mode + $inlist_mode + $infasta_mode) != 1) { 
  die "ERROR exactly one of the options --onlybuild, --inlist, --infasta must be used.";
}

# determine if we are going to annotate after classifying
my $do_annotate = opt_IsUsed("-A", \%opt_HH);
my $annotate_dir = ($do_annotate) ? opt_Get("-A", \%opt_HH) : undef;
my $annotate_non_cons_opts = undef; # defined below, after output file hashes set up, if --optsA used
my $annotate_reflist_file  = opt_IsUsed("--reflistA", \%opt_HH) ? opt_Get("--reflistA", \%opt_HH) : undef; 
if((defined $annotate_dir) && (! -d $annotate_dir)) { 
  die "ERROR with -A <s>, directory <s> must exist";
}
if(defined $annotate_reflist_file) { 
  if(! -e $annotate_reflist_file) { 
    die "ERROR $annotate_reflist_file, used with --reflistA, does not exist"; 
  }
  if(! -s $annotate_reflist_file) { 
    die "ERROR $annotate_reflist_file, used with --reflistA, exists but is empty";
  }
}

# enforce that lowscthresh >= vlowscthresh
if(opt_IsUsed("--lowscthresh",\%opt_HH) || opt_IsUsed("--vlowscthresh",\%opt_HH)) { 
  if(opt_Get("--lowscthresh",\%opt_HH) < opt_Get("--vlowscthresh",\%opt_HH)) { 
    die sprintf("ERROR, with --lowscthresh <x> and --vlowscthresh <y>, <x> must be less than <y> (got <x>: %f, <y>: %f)\n", 
                opt_Get("--lowscthresh",\%opt_HH), opt_Get("--vlowscthresh",\%opt_HH)); 
  }
}
# enforce that lowdiffthresh >= vlowdiffthresh
if(opt_IsUsed("--lowdiffthresh",\%opt_HH) || opt_IsUsed("--vlowdiffthresh",\%opt_HH)) { 
  if(opt_Get("--lowdiffthresh",\%opt_HH) < opt_Get("--vlowdiffthresh",\%opt_HH)) { 
    die sprintf("ERROR, with --lowdiffthresh <x> and --vlowdiffthresh <y>, <x> must be less than <y> (got <x>: %f, <y>: %f)\n", 
                opt_Get("--lowdiffthresh",\%opt_HH), opt_Get("--vlowdiffthresh",\%opt_HH)); 
  }
}

#############################
# create the output directory
#############################
my $cmd;              # a command to run with runCommand()
my @early_cmd_A = (); # array of commands we run before our log file is opened
# check if the $dirout exists, and that it contains the files we need
# check if our output dir $symbol exists
if($dir !~ m/\/$/) { $dir =~ s/\/$//; } # remove final '/' if it exists
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
my $out_root = $dir . "/" . $dir_tail . ".dnaorg_classify";

my $dir_build = undef;
if(opt_IsUsed("--dirbuild", \%opt_HH)) { 
  $dir_build = opt_Get("--dirbuild", \%opt_HH);
}
else { 
  $dir_build = $dir;
}

$dir_build =~ s/\/$//; # remove final '/' if there is one
my $dir_build_tail = $dir_build;
$dir_build_tail =~ s/^.+\///; # remove all but last dir
my $build_root = $dir_build . "/" . $dir_build_tail . ".dnaorg_classify";

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ();
my @arg_A      = ();
outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date, $dnaorgdir);
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
                         #  "list": list and description of output files
                         #  Per accession '$seq':
                         #                       $seq.tbl
                         #                       $seq.results

# open the log and command files 
openAndAddFileToOutputInfo(\%ofile_info_HH, "log", $out_root . ".log", 1, "Output printed to screen");
openAndAddFileToOutputInfo(\%ofile_info_HH, "cmd", $out_root . ".cmd", 1, "List of executed commands");
openAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".list", 1, "List and description of all output files");

my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
# output files are all open, if we exit after this point, we'll need
# to close these first.


# now we have the log file open, output the banner there too
outputBanner($log_FH, $version, $releasedate, $synopsis, $date, $dnaorgdir);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $cmd_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

my $be_verbose = opt_Get("-v", \%opt_HH);     # is -v used?
my $do_keep    = opt_Get("--keep", \%opt_HH); # should we leave intermediates files on disk, instead of removing them?
my @files2rm_A = ();                          # will be filled with files to remove, --keep was not enabled

###################################################
# make sure the required executables are executable
###################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"nhmmscan"}     = $hmmer_exec_dir . "nhmmscan";
$execs_H{"hmmbuild"}     = $hmmer_exec_dir . "hmmbuild";
$execs_H{"hmmpress"}     = $hmmer_exec_dir . "hmmpress";
$execs_H{"esl-reformat"} = $esl_exec_dir   . "esl-reformat";
$execs_H{"esl-seqstat"}  = $esl_exec_dir   . "esl-seqstat";
$execs_H{"esl-sfetch"}   = $esl_exec_dir   . "esl-sfetch";
$execs_H{"esl-ssplit"}   = $bioeasel_exec_dir . "esl-ssplit.pl";
if($do_annotate) { 
  $execs_H{"dnaorg_annotate"} = $dnaorg_exec_dir . "dnaorg_annotate.pl";
}
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

my @annotate_summary_output_A = (); # the output that summarizes number of sequences that pass each dnaorg_annotate.pl call

#################################################################################
#
# LES Jul 26 2016
#
# ASSUMPTIONS: $ref_list    is the file name of a list of RefSeq accns
#              $cls_list    is a file which contains the accession numbers of all the
#                           sequences which are to be assigned to ntlists (one accn #
#                           per line)
#
#################################################################################

my $progress_w = 70; # the width of the left hand column in our progress output, hard-coded                                                     
my $start_secs = outputProgressPrior("Parsing RefSeq list", $progress_w, $log_FH, *STDOUT);

# open and parse list of RefSeqs
# we do this in all modes, because we use ref_list_seqname_A and ref_fasta_seqname_A in all modes
my @ref_list_seqname_A  = (); # array of accessions listed in RefSeq input file
my @ref_fasta_seqname_A = (); # array of names of RefSeqs fetched in fasta file (may include version)
my $onlybuild_file    = opt_Get("--onlybuild", \%opt_HH);
my $ref_library       = $build_root . ".hmm";
my $ref_fa            = $build_root . ".ref.fa";
my $ref_stk           = $build_root . ".ref.stk";
my $ref_list          = $build_root . ".ref.list";
my $ref_fasta_list    = $build_root . ".ref.fa.list"; # fasta names (may include version)
my $ref_fasta_seqname; # name of a sequence in the fasta file to classify
my $ref_list_seqname;  # name of a sequence in the list file to classify (this is what we'll output)
my %ref_list2fasta_seqname_H = (); # hash mapping a list sequence name (key) to a fasta sequence name (value)
my $cur_nseq = 0;
my $nseq_above_zero = 0; # number of refseq accessions assigned > 0 sequences

# copy the ref list to the build directory if --onlybuild
if($onlybuild_mode) { 
  validateFileExistsAndIsNonEmpty($onlybuild_file, "main", $ofile_info_HH{"FH"});
  runCommand("cp $onlybuild_file $ref_list", opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
  addClosedFileToOutputInfo(\%ofile_info_HH, "RefList", $ref_list, 1, "List of reference sequences used to build HMMs");
}
fileLinesToArray($ref_list, 1, \@ref_list_seqname_A, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
my $n_ref = scalar(@ref_list_seqname_A);

# initialize the hashes that store which RefSeqs we annotate for
my %annotate_ref_list_seqname_H  = (); # key: value from ref_list_seqname_A,  value '1' if we are going to annotate for this ref_list_seqname, '0' if not
# initialize to 0 for all, first
foreach $ref_list_seqname (@ref_list_seqname_A) { 
  $annotate_ref_list_seqname_H{$ref_list_seqname} = 0; 
}

# if $do_annotate, we check here to make sure that all dnaorg_build.pl output 
# directories we may need to annotate exist
my %buildopts_used_HH = ();
my $build_dir = undef;
if($do_annotate) { 
  $start_secs = outputProgressPrior("Verifying build directories exist (-A)", $progress_w, $log_FH, *STDOUT);
  # deal with --reflistA if it was used
  if(defined $annotate_reflist_file) { # --reflistA $annotate_reflist_file used
    # read $annotate_reflist_file and update $annotate_ref_list_seqname_H{$ref_list_seqname} to 1 
    # for all listed accessions, we will only annotate sequences that match these.
    my @tmp_ref_list_A = (); # temporary array of lines in $annotate_reflist_file
    fileLinesToArray($annotate_reflist_file, 1, \@tmp_ref_list_A, $ofile_info_HH{"FH"});
    foreach $ref_list_seqname (@tmp_ref_list_A) { 
      if(! exists $annotate_ref_list_seqname_H{$ref_list_seqname}) { 
        DNAORG_FAIL("ERROR in dnaorg_classify.pl::main(), sequence name $ref_list_seqname read from $annotate_reflist_file does not exist in reference list file $ref_list", 1, $ofile_info_HH{"FH"});
      }
      $annotate_ref_list_seqname_H{$ref_list_seqname} = 1;
    }
  }
  else { # --reflistA not used, we will annotate all seqs that match all RefSeqs
    foreach $ref_list_seqname (@ref_list_seqname_A) { 
      $annotate_ref_list_seqname_H{$ref_list_seqname} = 1; 
    }
  }
  
  # verify dnaorg_build directories we need to annotate exist
  $annotate_dir =~ s/\/*$//; # remove trailing '/'
  foreach $ref_list_seqname (@ref_list_seqname_A) { 
    if($annotate_ref_list_seqname_H{$ref_list_seqname} == 1) { 
      $build_dir = $annotate_dir . "/" . $ref_list_seqname;
      %{$buildopts_used_HH{$ref_list_seqname}} = ();
      validate_build_dir($build_dir, \%{$buildopts_used_HH{$ref_list_seqname}}, \%opt_HH, \%ofile_info_HH);
    }
  }
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

if($onlybuild_mode) { 
  $start_secs = outputProgressPrior("Creating RefSeq HMM Library", $progress_w, $log_FH, *STDOUT);

  # create fasta file of all refseqs
  list_to_fasta($ref_list, $ref_fa, \%opt_HH, $ofile_info_HH{"FH"});
  addClosedFileToOutputInfo(\%ofile_info_HH, "RefFasta", $ref_fa, 1, "Fasta file of all RefSeq sequences");

  # get a list file with their actual names 
  fasta_to_list_and_lengths($ref_fa, $ref_fasta_list, $ref_fasta_list . ".tmp", $execs_H{"esl-seqstat"}, undef, undef, \%opt_HH, $ofile_info_HH{"FH"});
  addClosedFileToOutputInfo(\%ofile_info_HH, "RefFastaList", $ref_fasta_list, 1, "List of RefSeq names from fasta file (may include version)");
  fileLinesToArray($ref_fasta_list, 1, \@ref_fasta_seqname_A, $ofile_info_HH{"FH"});

  # create the hash mapping the list sequence names to the fasta sequence names
  foreach $ref_fasta_seqname (@ref_fasta_seqname_A) {
    $ref_list_seqname = fetchedNameToListName($ref_fasta_seqname);
    $ref_list2fasta_seqname_H{$ref_list_seqname} = $ref_fasta_seqname;
  }

  if(! -e $ref_fa . ".ssi") { 
    $cmd = $execs_H{"esl-sfetch"} . " --index $ref_fa > /dev/null";
    runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
  }

  for(my $i = 0; $i < $n_ref; $i++) { 
    my $hmm_file = $out_root . ".$i.hmm";
    $ref_list_seqname  = $ref_list_seqname_A[$i];
    $ref_fasta_seqname = $ref_list2fasta_seqname_H{$ref_list_seqname};
    if(! defined $ref_fasta_seqname) { 
      DNAORG_FAIL("ERROR in dnaorg_classify.pl::main() 0, could not find mapping fasta sequence name for list sequence name $ref_list_seqname", 1, $ofile_info_HH{"FH"});
    }
    # build up the command, we'll fetch the sequence, pass it into esl-reformat to make a stockholm 'alignment', and pass that into hmmbuild
    $cmd  = $execs_H{"esl-sfetch"} . " $ref_fa $ref_fasta_seqname | "; # fetch the sequence
#    $cmd .= $execs_H{"esl-reformat"} . " --informat afa stockholm - | "; # reformat to stockholm
    $cmd .= $execs_H{"hmmbuild"} . " --informat afa -n $ref_list_seqname $hmm_file - > /dev/null"; # build HMM file
    runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
    $cmd = "cat $hmm_file >> $ref_library";   # adds each individual hmm to the hmm library
    runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
    
    if(! $do_keep) {
      push(@files2rm_A, $hmm_file);
    }
  }    
  addClosedFileToOutputInfo(\%ofile_info_HH, "HMMLib", $ref_library, 1, "Library of HMMs of RefSeqs");

  $cmd = $execs_H{"hmmpress"} . " $ref_library > /dev/null";
  runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
  addClosedFileToOutputInfo(\%ofile_info_HH, "HMMLib.press.m", $ref_library . ".h3m", 0, "HMM press index file (.h3m)");
  addClosedFileToOutputInfo(\%ofile_info_HH, "HMMLib.press.f", $ref_library . ".h3f", 0, "HMM press index file (.h3f)");
  addClosedFileToOutputInfo(\%ofile_info_HH, "HMMLib.press.p", $ref_library . ".h3p", 0, "HMM press index file (.h3p)");
  addClosedFileToOutputInfo(\%ofile_info_HH, "HMMLib.press.i", $ref_library . ".h3i", 0, "HMM press index file (.h3i)");

  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}
else { 
  # not --onlybuild, classify sequences
  # first, define some file names, copy files and make sure all required files exist and are nonempty
  my $infasta_file        = opt_Get("--infasta", \%opt_HH);
  my $inlist_file         = opt_Get("--inlist", \%opt_HH);
  my $cls_list            = $out_root . ".seqlist";
  my $cls_fa              = $out_root . ".fa";
  my %cls_fasta_seqlen_H  = (); # key is sequence name from fasta file, value is length of sequence
  my @cls_fasta_seqname_A = (); # array of sequences in $cls_fa,   may be in "accession-version" format
  my @cls_list_seqname_A  = (); # array of sequences in $cls_list, may not include version
  my $n_cls_seqs          = 0;  # number of sequences to classify, size of @cls_fasta_seqname_A and @cls_list_seqname_A
  my $ref_list_seqname;  # name of a reference sequence
  my $cls_fasta_seqname; # name of a sequence in the fasta file to classify
  my $cls_fasta_seqlen;  # length of $cls_fasta_seqname
  my $cls_list_seqname;  # name of a sequence in the list file to classify (this is what we'll output)
  my %cls_list2fasta_seqname_H = (); # hash mapping a list sequence name (key) to a fasta sequence name (value)

  # if we're annotating, read the dnaorg_annotate.pl options file, if nec
  if(($do_annotate) && (opt_IsUsed("--optsA", \%opt_HH))) { 
    $start_secs = outputProgressPrior("Parsing additional dnaorg_annotate.pl options (--optsA)", $progress_w, $log_FH, *STDOUT);
    my $failure_str   = "that option will automatically be set,\nas required, to be consistent with the relevant dnaorg_build.pl command used previously.";
    my $auto_add_opts = "--dirout,--dirbuild,--infasta,--refaccn";
    if($do_keep)    { $auto_add_opts .= ",--keep"; }
    if($be_verbose) { $auto_add_opts .= ",-v"; }
    $annotate_non_cons_opts = parseNonConsOptsFile(opt_Get("--optsA", \%opt_HH), "--optsA", $auto_add_opts, $failure_str, $ofile_info_HH{"FH"});
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }

  if($inlist_mode) { 
    $start_secs = outputProgressPrior("Fetching sequences to a fasta file", $progress_w, $log_FH, *STDOUT);
    # copy the inlist file to our output dir
    validateFileExistsAndIsNonEmpty($inlist_file, "main", $ofile_info_HH{"FH"});
    $cmd = "cp $inlist_file $cls_list";
    runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
    validateFileExistsAndIsNonEmpty($cls_list, "main", $ofile_info_HH{"FH"});
    addClosedFileToOutputInfo(\%ofile_info_HH, "ClassList", $cls_list, 1, "List file with sequences to classify (copy of $inlist_file)");

    # create the fasta file using efetch
    list_to_fasta($cls_list, $cls_fa, \%opt_HH, $ofile_info_HH{"FH"});
    validateFileExistsAndIsNonEmpty($cls_fa, "main", $ofile_info_HH{"FH"});
    addClosedFileToOutputInfo(\%ofile_info_HH, "SeqFasta", $cls_fa, 1, "Fasta file with sequences to classify");

    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

    # get sequence lengths 
    $start_secs = outputProgressPrior("Processing fasta file to get sequence lengths", $progress_w, $log_FH, *STDOUT);
    fasta_to_list_and_lengths($cls_fa, undef, $cls_list . ".tmp", $execs_H{"esl-seqstat"}, \@cls_fasta_seqname_A, \%cls_fasta_seqlen_H, \%opt_HH, $ofile_info_HH{"FH"});

    # store sequence names in list file
    fileLinesToArray($cls_list, 1, \@cls_list_seqname_A, $ofile_info_HH{"FH"});
    $n_cls_seqs = scalar(@cls_list_seqname_A);
    if($n_cls_seqs != (scalar(@cls_fasta_seqname_A))) { 
      DNAORG_FAIL(sprintf("ERROR in dnaorg_classify.pl::main(), fetched %d != %d seqs, when %d were listed in %s\n", scalar(@cls_fasta_seqname_A), $n_cls_seqs, $n_cls_seqs, $cls_list), 1, $ofile_info_HH{"FH"});
    }

    # create the hash mapping the list sequence names to the fasta sequence names
    foreach $cls_fasta_seqname (@cls_fasta_seqname_A) {
      $cls_list_seqname = fetchedNameToListName($cls_fasta_seqname);
      $cls_list2fasta_seqname_H{$cls_list_seqname} = $cls_fasta_seqname;
    }
    # make sure we have a valid mapping for all sequences in the list file
    foreach $cls_list_seqname (@cls_list_seqname_A) { 
      if(! defined $cls_list2fasta_seqname_H{$cls_list_seqname}) { 
        DNAORG_FAIL("ERROR in dnaorg_classify.pl::main() 1, could not find mapping fasta sequence name for list sequence name $cls_list_seqname", 1, $ofile_info_HH{"FH"});
      }
    }
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
  elsif($infasta_mode) { 
    $start_secs = outputProgressPrior("Processing fasta file to get sequence lengths", $progress_w, $log_FH, *STDOUT);
    # copy the fasta file to our output dir
    validateFileExistsAndIsNonEmpty($infasta_file, "main", $ofile_info_HH{"FH"});
    $cmd = "cp $infasta_file $cls_fa";
    runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
    validateFileExistsAndIsNonEmpty($cls_fa, "main", $ofile_info_HH{"FH"});
    addClosedFileToOutputInfo(\%ofile_info_HH, "SeqFasta", $cls_fa, 1, "Fasta file with sequences to classify (copy of $infasta_file)");

    if(! -e $cls_fa . ".ssi") { 
      $cmd = $execs_H{"esl-sfetch"} . " --index $cls_fa > /dev/null";
      runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
    }

    # create the list file and get sequence lengths using esl-seqstat    
    fasta_to_list_and_lengths($cls_fa, $cls_list, $cls_list . ".tmp", $execs_H{"esl-seqstat"}, \@cls_fasta_seqname_A, \%cls_fasta_seqlen_H, \%opt_HH, $ofile_info_HH{"FH"});
    validateFileExistsAndIsNonEmpty($cls_list, "main", $ofile_info_HH{"FH"});
    addClosedFileToOutputInfo(\%ofile_info_HH, "ClassList", $cls_list, 0, "List file with sequences to classify");

    @cls_list_seqname_A = @cls_fasta_seqname_A; # list seqname array is same as fasta seqname array if --infasta used
    $n_cls_seqs = scalar(@cls_list_seqname_A);

    # create the hash mapping the list sequence names to the fasta sequence names
    foreach $cls_fasta_seqname (@cls_fasta_seqname_A) {
      $cls_list2fasta_seqname_H{$cls_fasta_seqname} = $cls_fasta_seqname; # names are identical in both arrays
    }
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }

  # now regardless of whether --inlist or --infasta was chosen we're at the same spot
  # $cls_list:            list file of all seqs to classify
  # @cls_list_seqname_A:  array of all sequence names from $cls_list
  # $cls_fa:              fasta file of all seqs to classify
  # @cls_fasta_seqname_A: array of all sequence names from $cls_fa
  # %cls_list2fasta_seqname_H: hash mapping list sequence names to fasta sequence
  

  ########### RUN nhmmscan and generate output files ####################################################################################################
  my $tblout_file = $out_root . ".tblout"; # concatenated tblout file, created by concatenating all of the individual 
                                           # tblout files in cmscanOrNhmmscanWrapper()
  my @mdl_file_A = ($ref_library); # cmscanOrNhmmscanWrapper() needs an array of model files

  my $cls_tot_len_nt = sumHashValues(\%cls_fasta_seqlen_H);
  cmscanOrNhmmscanWrapper(\%execs_H, 0, $out_root, $cls_fa, $cls_tot_len_nt, $tblout_file, $progress_w, 
                          \@mdl_file_A, undef, \%opt_HH, \%ofile_info_HH); 
  # in above cmscanOrNhmmscanWrapper call: '0' means run nhmmscan, not cmscan, 'undef' is for the model length array, irrelevant b/c we're using nhmmscan

  ########################################################################################################################################################
  #
  # PARSE nhmmscan files
  # loop through sequences
  # loop through lines
  #
  # Chooses a RefSeq for each sequence
  # Creates: .tbl file giving information on each match
  # Creates: ntlists for each RefSeq
  #
  ########################################################################################################################################################

  $start_secs = outputProgressPrior("Creating tabular output file", $progress_w, $log_FH, *STDOUT);

  # determine which unexpected features cause a sequence fo tfail
  my $lowsc_fails     = (opt_Get("--allpass", \%opt_HH)) || (  opt_Get("--lowscpass",    \%opt_HH)) ? 0 : 1;
  my $vlowsc_fails    = (opt_Get("--allpass", \%opt_HH)) || (  opt_Get("--vlowscpass",   \%opt_HH)) ? 0 : 1;
  my $lowdiff_fails   = (opt_Get("--allpass", \%opt_HH)) || (  opt_Get("--lowdiffpass",  \%opt_HH)) ? 0 : 1;
  my $vlowdiff_fails  = (opt_Get("--allpass", \%opt_HH)) || (  opt_Get("--vlowdiffpass", \%opt_HH)) ? 0 : 1;
  my $bias_fails      = (opt_Get("--allpass", \%opt_HH)) || (  opt_Get("--biaspass",     \%opt_HH)) ? 0 : 1;
  my $minus_fails     = (opt_Get("--allpass", \%opt_HH)) || (! opt_Get("--minusfail",    \%opt_HH)) ? 0 : 1;

  my $lowsc_minlen    = opt_Get("--lowscminlen", \%opt_HH); 
  my $lowdiff_minlen  = opt_Get("--lowdiffminlen", \%opt_HH); 
  if(opt_Get("--nolowscminlen",   \%opt_HH)) { $lowsc_minlen   = 0; }
  if(opt_Get("--nolowdiffminlen", \%opt_HH)) { $lowdiff_minlen = 0; }

  my $query_width = 20;
  foreach $cls_fasta_seqname (keys %cls_fasta_seqlen_H) { 
    if(length($cls_fasta_seqname) > $query_width) { 
      $query_width = length($cls_fasta_seqname); 
    }
  }

  # Generate match information file header                                                                               
  my $match_file = $out_root . ".infotbl";
  open(MATCH_INFO, ">", $match_file) || fileOpenFailure($match_file, $0, $!, "writing", $ofile_info_HH{"FH"});
  print  MATCH_INFO "########################################################################################################################################\n";
  print  MATCH_INFO  "#\n";
  print  MATCH_INFO  "# query:       Accession number of the sequence\n";
  print  MATCH_INFO  "# qlen:        length of this sequence\n";
  print  MATCH_INFO  "# RefSeq:      RefSeq that the sequence was assigned to\n";
  print  MATCH_INFO  "# score:       (summed) bit score(s) of all hits to 'RefSeq'\n";
  print  MATCH_INFO  "# sc/nt:       'score' divided by 'qlen'\n";
  print  MATCH_INFO  "# E-val:       E-value of hit to 'RefSeq'\n";
  print  MATCH_INFO  "# coverage:    the percentage of the query that the hit to 'RefSeq' covers (Hit length)/(Query length)\n";
  print  MATCH_INFO  "# bias:        correction in bits for biased composition sequences\n";
  print  MATCH_INFO  "# #hits:       number of individual hits to this RefSeq (the program combines stats such as bit score and covg. from separate hits)\n";
  print  MATCH_INFO  "# strand:      strand of top hit to 'RefSeq'\n";
  print  MATCH_INFO  "# H2:RefSeq:   RefSeq that had the second strongest hit\n";
  print  MATCH_INFO  "# scdiff:      difference in bit score b/t 'RefSeq' hit(s) and 'H2:RefSeq' hit(s)\n";
  print  MATCH_INFO  "# scdiff/nt:   'scdiff' divided by 'qlen'\n";
  print  MATCH_INFO  "# covdiff:     amount by which the coverage of the 'RefSeq' hit is greater than that of the 'H2: RefSeq' hit\n";
  printf MATCH_INFO ("# CC:          'confidence class', first letter based on sc/nt: A: if sc/nt >= %.3f, B: if %.3f > sc/nt >= %.3f, C: if %.3f > sc_nt\n", opt_Get("--lowscthresh", \%opt_HH), opt_Get("--lowscthresh", \%opt_HH), opt_Get("--vlowscthresh", \%opt_HH), opt_Get("--vlowscthresh", \%opt_HH));
  printf MATCH_INFO ("#              second letter based on diff/nt: A: if diff/nt >= %.3f, B: if %.3f > diff/nt >= %.3f, C: if %.3f > diff_nt\n", opt_Get("--lowdiffthresh", \%opt_HH), opt_Get("--lowdiffthresh", \%opt_HH), opt_Get("--vlowdiffthresh", \%opt_HH), opt_Get("--vlowdiffthresh", \%opt_HH));
  printf MATCH_INFO ("# p/f:         'PASS' if sequence passes, 'FAIL' if it fails\n");
  print  MATCH_INFO  "# unexpected   \n";
  print  MATCH_INFO  "# features:    unexpected features for this sequence\n";
  print  MATCH_INFO  "#              Possible values in unexpected features column:\n";
  printf MATCH_INFO ("#              LowScore:     'sc/nt'   < %.3f (threshold settable with --lowscthresh)\n",    opt_Get("--lowscthresh", \%opt_HH));
  printf MATCH_INFO ("#              VeryLowScore: 'sc/nt'   < %.3f (threshold settable with --vlowscthresh)\n",   opt_Get("--vlowscthresh", \%opt_HH));
  printf MATCH_INFO ("#              LowDiff:      'diff/nt' < %.3f (threshold settable with --lowdiffthresh)\n",  opt_Get("--lowdiffthresh", \%opt_HH));
  printf MATCH_INFO ("#              VeryLowDiff:  'diff/nt' < %.3f (threshold settable with --vlowdiffthresh)\n", opt_Get("--vlowdiffthresh", \%opt_HH));
  printf MATCH_INFO ("#              MinusStrand:  top hit is on minus strand\n");
  printf MATCH_INFO ("#              HighBias:     'bias' > (%.3f * ('bias' + 'score')) (threshold settable with --biasfract)\n", opt_Get("--biasfract", \%opt_HH));
  print  MATCH_INFO  "########################################################################################################################################\n";
  print  MATCH_INFO "#\n";
  printf MATCH_INFO ("%-*s  %6s  %9s  %7s  %5s  %8s  %8s  %7s  %5s  %5s  %9s  %7s  %7s  %7s  %2s  %4s  %s\n", 
                     $query_width, "#query", "qlen", "RefSeq", "score", "sc/nt", "E-val", "coverage", "bias", "#hits", "strand", "H2:RefSeq", "scdiff", "diff/nt", "covdiff", "CC", "p/f", "unexpected-features");
  printf MATCH_INFO ("%-*s  %6s  %9s  %7s  %5s  %8s  %8s  %7s  %5s  %6s  %9s  %7s  %7s  %7s  %2s  %4s  %s\n", 
                     $query_width,
                     "#===================", 
                     "======",
                     "=========",
                     "=======",
                     "=====",
                     "========",
                     "========",
                     "=======",
                     "=====",
                     "======",
                     "=========",
                     "=======",
                     "=======",
                     "=======",
                     "==",
                     "====",
                     "===================");

  # Generate data structures to build ntlists from
  my %ntlist_HA = (); # Hash of arrays containing each RefSeq's ntlist
                       # Key:    RefSeq accession #
                      # Value:  Array of accessions for which have been assigned to this RefSeq
  foreach (@ref_list_seqname_A) {
    @{$ntlist_HA{$_}} = (); # initialize each RefSeq's value to an empty array
  }
  
  # initialize the non-assigned list to an empty array
  @{$ntlist_HA{"non-assigned"}} = ();

  my %pass_fail_H = (); # Hash of pass/fail values:
                        # Key: Sequence accession
                        # Value: "PASS" or "FAIL"
  
  # Generate data structures to build a sequence profile for further evaluation                       
  my %hit_info_HAA = (); # Hash of hash of arrays containing nhmmscan info for each hit for each sequence
                         # 1D Key:    Sequence accession #                                                                             
                         # 1D Array:  Order of hit (1st nhmmscan hit = index 0, 2nd = index 1, etc)
                         # Value:     Array w/ relevant nhmmscan info for the hit (see comments below)
  
  my @hit_output_A = (); # array containing all the output for this hit
  my %counter_H    = (); # key is sequence name, value is number of hits for this sequence
  my $h;                 # counter over hits
  my $i;                 # counter over sequences

  # indices in @hit_output_A:
  my $HOIDX_QUERY      = 0;
  my $HOIDX_QLEN       = 1;
  my $HOIDX_REFSEQ     = 2;
  my $HOIDX_BITSC      = 3;
  my $HOIDX_BITSCPNT   = 4;
  my $HOIDX_EVALUE     = 5;
  my $HOIDX_COVERAGE   = 6;
  my $HOIDX_BIAS       = 7;
  my $HOIDX_NUMHITS    = 8;
  my $HOIDX_STRAND     = 9;
  my $HOIDX_REFSEQ2    = 10;
  my $HOIDX_BITDIFF    = 11;
  my $HOIDX_BITDIFFPNT = 12;
  my $HOIDX_COVDIFF    = 13;

  open(NHMMSCANTBL, $tblout_file) || fileOpenFailure($tblout_file, $0, $!, "reading", $ofile_info_HH{"FH"});
  while(my $line = <NHMMSCANTBL>) {
    if($line !~ m/^\#/) {               #if this line is not commented
      my @hit_specs_A = split(/\s+/, $line);
      $ref_list_seqname = $hit_specs_A[0];
      $cls_fasta_seqname = $hit_specs_A[2];
      if(! exists $cls_fasta_seqlen_H{$cls_fasta_seqname}) { 
        DNAORG_FAIL("ERROR in dnaorg_classify.pl::main() sequence length for sequence $cls_fasta_seqname does not exist in the length hash", 1, $ofile_info_HH{"FH"}, 1, $ofile_info_HH{"FH"});
      }
      my $cls_fasta_seqlen = $cls_fasta_seqlen_H{$cls_fasta_seqname};

      ##########################################################################
      # In @hit_specs_A - array containing all the output info for the hit on this line
      #
      # index         feature
      #-----------------------
      # 0             target name
      # 1             accession (not sure what this refers to, shows up as a - for every hit in my output)
      # 2             query name (name from > line of fasta file)
      # 3             accession (again, ?)
      # 4             hmmfrom
      # 5             hmm to
      # 6             alifrom
      # 7             ali to
      # 8             env from
      # 9             env to
      # 10            modlen (entire length of model seq)
      # 11            strand (+,-)
      # 12            E-value
      # 13            bit score
      # 14            bias
      # 15            description of target (showing up as - for every hit in my output)
      ##########################################################################
        
      # ########################################################################
      # In @hit_output_A - array containing all the output for this hit
      #
      # index variable     feature
      #----------------------------
      # 0 $HOIDX_QUERY     Query accn
      # 1 $HOIDX_QLEN      Query accn sequence length
      # 2 $HOIDX_REFSEQ    RefSeq accn
      # 3 $HOIDX_BITSC     Bit-score
      # 4 $HOIDX_BITSCPNT  bits per nucleotide (Bit-score / sequence length)
      # 5 $HOIDX_EVALUE    E-value
      # 6 $HOIDX_COVERAGE  Coverage (Hit length / Query length)
      # 7 $HOIDX_BIAS      Bias
      # 8 $HOIDX_NUMHITS   Number of hits to this RefSeq
      # 9 $HODIX_STRAND    Strand of top hit (and all hits that combine to make-up coverage -- we skip any hits on the other strand)
      # 10-13 added below
      ###########################################################################
      # check if this RefSeq has appeared in a previous hit for this sequence                                              

      # initialize if necessary
      if(! exists($hit_info_HAA{$cls_fasta_seqname}) ) {
	@{$hit_info_HAA{$cls_fasta_seqname}} = ();
        $counter_H{$cls_fasta_seqname} = 0;
      }

      my $first_index = undef;  # the first occuring index of a previous hit, if one exists                                                             
      # if there are previous hits, search through them to see if this RefSeq has been hit before
      if( @{$hit_info_HAA{$cls_fasta_seqname}} != 0 ) {
        for($h=0;  $h < scalar(@{$hit_info_HAA{$cls_fasta_seqname}});  $h++) {
          if( $hit_info_HAA{$cls_fasta_seqname}[$h][1] eq $ref_list_seqname ) {       # if index $i's RefSeq is the same as this one    
            $first_index = $h;
            $h = scalar(@{$hit_info_HAA{$cls_fasta_seqname}});
          }
        }
      }
    
      # Deciding hit length from alito alifrom - will be used to find Coverage
      my $hit_len = abs($hit_specs_A[7] - $hit_specs_A[6]) + 1;

      # if this is the first hit to this RefSeq
      if(! defined($first_index)) {
        #debug
        #print "\tfirst occurence of $cls_fasta_seqname \n";
        
        @hit_output_A = (); # Array to be outputed to match info file
        push(@hit_output_A, $cls_fasta_seqname);
        push(@hit_output_A, sprintf("%d", $cls_fasta_seqlen));
        push(@hit_output_A, $ref_list_seqname);
        
        my $bitsc = $hit_specs_A[13];
        push(@hit_output_A, sprintf("%7.1f", $bitsc)); # add bit score

        my $bitsc_per_nt = $bitsc / $cls_fasta_seqlen;
        push(@hit_output_A, sprintf("%5.3f", $bitsc_per_nt)); # add bit score per nucleotide
        
        my $e_val = $hit_specs_A[12];
        push(@hit_output_A, $e_val); # add E-val                                                                                                            
        
        # add coverage				
        my $coverage = $hit_len / $cls_fasta_seqlen;
        $coverage = sprintf("%5.3f", $coverage);
        push(@hit_output_A, $coverage);
        
        # add bias
        my $bias = $hit_specs_A[14];
        push(@hit_output_A, $bias); # initialize 'Number of hits' to 1

        # initialize number of hits
        push(@hit_output_A, 1); # initialize 'Number of hits' to 1

        # add strand
        my $strand = $hit_specs_A[11];
        push(@hit_output_A, $hit_specs_A[11]);  # add strand
                    
        @{$hit_info_HAA{$cls_fasta_seqname}}[$counter_H{$cls_fasta_seqname}] = ();
        @{$hit_info_HAA{$cls_fasta_seqname}[$counter_H{$cls_fasta_seqname}]} = @hit_output_A;
          
        $counter_H{$cls_fasta_seqname}++;
          
      }else {    # if this is a hit to a sequence that has already been hit
        #@hit_output_A = \@{$hit_info_HAA{$cls_fasta_seqname}[$first_index]};
          
        #debug
        #print "\tsecond, third, etc hit of $cls_fasta_seqname \n";
          
        # TODO - check with Eric!
        if(($hit_specs_A[13] > 0) && # only add hits with positive bit scores 
           ($hit_specs_A[11] eq $hit_info_HAA{$cls_fasta_seqname}[$first_index][7])) { # only add hits on the same strand as the top hit
          @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[$HOIDX_BITSC] += $hit_specs_A[13]; # Add bit score
          @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[$HOIDX_BIAS]  += $hit_specs_A[14]; # Add bias
          @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[$HOIDX_NUMHITS] ++;                # increase 'Number of hits' by 1
          # recalculate bits per nucleotide
          @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[$HOIDX_BITSCPNT] = @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[$HOIDX_BITSC] / $cls_fasta_seqlen;

          ##############################################
          # Calculate and add coverage
          my $prev_covg = @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[4];
          my $prev_hit_len = $prev_covg*$cls_fasta_seqlen;
          printf("$cls_fasta_seqname prev_covg: $prev_covg prev_hit_len: $prev_hit_len\n");
          my $coverage = ($prev_hit_len + $hit_len) / $cls_fasta_seqlen;
          $coverage = sprintf("%.7f", $coverage);
          @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[4] = $coverage;
          #############################################
        }
      }  
    }
  }
  close NHMMSCANTBL || die "couldn't close NHMMSCANTBL. oops! \n";

    
  # now we have everything stored in $hit_info_HAA
  # go through and for each sequence create the output
  for($i = 0; $i < $n_cls_seqs; $i++) { 
    $cls_list_seqname  = $cls_list_seqname_A[$i];
    $cls_fasta_seqname = $cls_list2fasta_seqname_H{$cls_list_seqname};
    if(! defined $cls_fasta_seqname) { 
      DNAORG_FAIL("ERROR in dnaorg_classify.pl::main() 2, could not find mapping fasta sequence name for list sequence name $cls_list_seqname", 1, $ofile_info_HH{"FH"});
    }
    if(exists($hit_info_HAA{$cls_fasta_seqname})) { 
      # there was at least one valid hit for this sequence
      # parse and output the best hit's info

      # sort order of hits by bit score
      @{$hit_info_HAA{$cls_fasta_seqname}} = sort { $a->[$HOIDX_BITSC] <=> $b->[$HOIDX_BITSC] } @{$hit_info_HAA{$cls_fasta_seqname}};
      
      # Assigns the sequence to the refseq hit with the the highest combined bit score of all hits to that refseq
      #my $max_index = 0; # contains the index of the hit with the highest bit score - initialized at index 0
      #my $second_index = 0;
      #for( my $hit=0; $hit< scalar(@{$hit_info_HAA{$cls_fasta_seqname}}) ; $hit++  ) {
      #	if( @{$hit_info_HAA{$cls_fasta_seqname}[$hit]}[2] > @{$hit_info_HAA{$cls_fasta_seqname}[$max_index]}[2] ) {
      #	    $max_index = $hit;
      #	}
      #}
      
      # assign hit_output_A to the assigned refseq hit's info, assign $ref_list_seqname to the assigned refseq
      @hit_output_A = @{@{$hit_info_HAA{$cls_fasta_seqname}}[-1]}; 
      $ref_list_seqname = $hit_output_A[$HOIDX_REFSEQ];
      
      ##############################################################################
      # Add comparison indicies to @hit_output_A
      #
      # index variable        feature
      #--------------------------------
      # 10 $HOIDX_REFSEQ2     RefSeq accn # for hit 2
      # 11 $HOIDX_BITDIFF     difference between bit scores of hit 1 and hit 2
      # 12 $HOIDX_BITDIFFPNT  bit difference per nucleotide (bitdiff / sequence length)
      # 13 $HOIDX_COVDIFF     difference between coverage of hit 1 and hit 2
      ##############################################################################
      
      # TODO - hardcode-y 
      my $bit_diff_per_nt = undef;
      if( scalar(@{$hit_info_HAA{$cls_fasta_seqname}}) >= 2) {     # if there is more than one hit  
        push(@hit_output_A, @{@{$hit_info_HAA{$cls_fasta_seqname}}[-2]}[$HOIDX_REFSEQ]);
	
        my $bit_diff = $hit_info_HAA{$cls_fasta_seqname}[-1][$HOIDX_BITSC] - $hit_info_HAA{$cls_fasta_seqname}[-2][$HOIDX_BITSC];
        $bit_diff = sprintf("%7.1f", $bit_diff);
        $bit_diff_per_nt = sprintf("%5.3f", ($bit_diff / $cls_fasta_seqlen_H{$cls_fasta_seqname}));
        my $cov_diff = $hit_info_HAA{$cls_fasta_seqname}[-1][$HOIDX_COVERAGE] - $hit_info_HAA{$cls_fasta_seqname}[-2][$HOIDX_COVERAGE];
        $cov_diff = sprintf("%5.3f", $cov_diff);
        push(@hit_output_A, $bit_diff);
        push(@hit_output_A, $bit_diff_per_nt);
        push(@hit_output_A, $cov_diff);
 	
      } else { # if there was only one hit, there's no comparison
        push(@hit_output_A, "-----");
        push(@hit_output_A, "-----");
        push(@hit_output_A, "-----");
        push(@hit_output_A, "-----");
      }

      # determine ufeature string and whether the sequence passes or fails
      my $ufeature_str = "";
      my $top_bitsc      = $hit_info_HAA{$cls_fasta_seqname}[-1][$HOIDX_BITSC];
      my $top_bias       = $hit_info_HAA{$cls_fasta_seqname}[-1][$HOIDX_BIAS];
      my $top_strand     = $hit_info_HAA{$cls_fasta_seqname}[-1][$HOIDX_STRAND];
      my $top_bit_per_nt = $hit_info_HAA{$cls_fasta_seqname}[-1][$HOIDX_BITSCPNT];
      my $score_class = "A"; # set to 'B' or 'C' below if below threshold
      my $diff_class  = "A"; # set to 'B' or 'C' below if below threshold
      my $seqlen = $cls_fasta_seqlen_H{$cls_fasta_seqname};

      $pass_fail_H{$cls_fasta_seqname} = "PASS"; # will change to FAIL below if necessary
      
      if($top_bit_per_nt < (opt_Get("--vlowscthresh", \%opt_HH))) { 
        $ufeature_str .= "VeryLowScore(" . $top_bit_per_nt . "<" . sprintf("%.3f", opt_Get("--vlowscthresh", \%opt_HH)) . ");"; 
        $score_class = "C";
        if($vlowsc_fails) { $pass_fail_H{$cls_fasta_seqname} = "FAIL"; }
      }
      elsif($top_bit_per_nt < (opt_Get("--lowscthresh", \%opt_HH))) { 
        $ufeature_str .= "LowScore(" . $top_bit_per_nt . "<" . sprintf("%.3f", opt_Get("--lowscthresh", \%opt_HH)) . ");"; 
        $score_class = "B";
        if(($lowsc_fails) && ($seqlen > $lowsc_minlen)) { $pass_fail_H{$cls_fasta_seqname} = "FAIL"; }
      }
      if(defined $bit_diff_per_nt) { 
        if($top_bit_per_nt < (opt_Get("--vlowdiffthresh", \%opt_HH))) { 
          $ufeature_str .= "VeryLowDiff(" . $bit_diff_per_nt . "<" . sprintf("%.3f", opt_Get("--vlowdiffthresh", \%opt_HH)) . ");"; 
          $diff_class = "C";
          if($vlowdiff_fails) { $pass_fail_H{$cls_fasta_seqname} = "FAIL"; }
        }
        elsif($bit_diff_per_nt < (opt_Get("--lowdiffthresh", \%opt_HH))) { 
          $ufeature_str .= "LowDiff(" . $bit_diff_per_nt . "<" . sprintf("%.3f", opt_Get("--lowdiffthresh", \%opt_HH)) . ");"; 
          $diff_class = "B";
          if($lowdiff_fails && ($seqlen > $lowdiff_minlen)) { $pass_fail_H{$cls_fasta_seqname} = "FAIL"; }
        }
      }
      if($top_strand eq "-") { 
        $ufeature_str .= "MinusStrand;";
        if($minus_fails) { $pass_fail_H{$cls_fasta_seqname} = "FAIL"; }
      }
      if($top_bias > ((opt_Get("--biasfract", \%opt_HH)) * ($top_bitsc + $top_bias))) { 
        # $top_bitsc has already had bias subtracted from it so we need to add it back in before we compare with biasfract
        $ufeature_str .= "HighBias;";
        if($bias_fails) { $pass_fail_H{$cls_fasta_seqname} = "FAIL"; }
      }
      if($ufeature_str eq "") { 
        $ufeature_str = "-";
      }

      # add hit info to .infotbl file
      #print MATCH_INFO join("\t", (@hit_output_A, $ufeature_str));
      printf MATCH_INFO ("%-*s  %6s  %9s  %7s  %5s  %8s  %8s  %7s  %5s  %6s  %9s  %7s  %7s  %7s  %2s  %4s  %s\n", 
                         $query_width,
                         $hit_output_A[$HOIDX_QUERY], 
                         $hit_output_A[$HOIDX_QLEN], 
                         $hit_output_A[$HOIDX_REFSEQ], 
                         $hit_output_A[$HOIDX_BITSC], 
                         $hit_output_A[$HOIDX_BITSCPNT], 
                         $hit_output_A[$HOIDX_EVALUE], 
                         $hit_output_A[$HOIDX_COVERAGE],
                         $hit_output_A[$HOIDX_BIAS],
                         $hit_output_A[$HOIDX_NUMHITS],
                         ($hit_output_A[$HOIDX_STRAND] eq "+") ? "plus" : "minus",
                         $hit_output_A[$HOIDX_REFSEQ2],
                         $hit_output_A[$HOIDX_BITDIFF],
                         $hit_output_A[$HOIDX_BITDIFFPNT],
                         $hit_output_A[$HOIDX_COVDIFF], 
                         $score_class . $diff_class,
                         $pass_fail_H{$cls_fasta_seqname},
                         $ufeature_str);
      

      # Add this seq to the hash key that corresponds to its RefSeq
      push(@{$ntlist_HA{$ref_list_seqname}}, $cls_list_seqname);
      # IMPORTANT: push $cls_list_seqname, and not $cls_fasta_seqname (cls_list_seqname is from the input list 
      # (unless --infasta) and may not include the version
    } else { # if there were no eligible hits
      printf MATCH_INFO ("%-*s  %6s  %9s  %7s  %5s  %8s  %8s  %7s  %5s  %6s  %9s  %7s  %7s  %7s  %2s  %4s  %s\n", 
                         $query_width,
                         $cls_fasta_seqname,
                         $cls_fasta_seqlen_H{$cls_fasta_seqname},
                         "-",
                         "-",
                         "-",
                         "-",
                         "-",
                         "-",
                         "-", 
                         "-", 
                         "-", 
                         "-",
                         "-",
                         "-",
                         "--",
                         "FAIL",
                         "NoHits;");

      $pass_fail_H{$cls_fasta_seqname} = "FAIL";

      # add this sequence to the file that lists non-assigned sequences
      push(@{$ntlist_HA{"non-assigned"}}, $cls_list_seqname); 
      # IMPORTANT: push $cls_list_seqname, and not $cls_fasta_seqname (cls_list_seqname is from the input list 
      # (unless --infasta) and may not include the version
    }
  } # end of 'foreach $cls_fasta_seqname' loop

  addClosedFileToOutputInfo(\%ofile_info_HH, "match-info", $match_file, 1, "Table with statistics for each match");
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  ########################################################################
  $start_secs = outputProgressPrior("Creating ntlists and other output files", $progress_w, $log_FH, *STDOUT);

  my @tmp_output_A = (); # array of lines to output after we create the files we want to create
  
  # Generate ntlist files, one per RefSeq
  # If --infasta used: generate fasta files, one per RefSeq with at least 
  if($infasta_mode) { 
    if(! -s $cls_fa . ".ssi") { 
      DNAORG_FAIL("ERROR in dnaorg_classify.pl::main(), ssi file for $cls_fa does not exist when trying to fetch", 1, $ofile_info_HH{"FH"});
    }
  }

  push(@tmp_output_A, sprintf("#\n"));
  push(@tmp_output_A, sprintf("# Number of input sequences assigned to each RefSeq:\n"));
  push(@tmp_output_A, sprintf("#\n"));
  push(@tmp_output_A, sprintf("%-12s  %10s  %10s    %10s  %10s%s\n", "# RefSeq-acc", "num-PASS", "fract-PASS", "num-FAIL", "fract-FAIL", ($do_annotate) ? "    annotating?" : ""));
  push(@tmp_output_A, sprintf("%-12s  %10s  %10s    %10s  %10s%s\n", "#-----------", "----------", "----------", "----------", "----------", ($do_annotate) ? "    -----------" : ""));
  my $cur_nseq_pass = 0; # number of seqs for current refseq that pass
  my $cur_nseq_fail = 0; # number of seqs for current refseq that fail

  my %ntlist_file_H    = (); # two keys: PASS and FAIL, name of ntlist file for PASS/FAIL seqs
  my %dv_ntlist_file_H = (); # two keys: PASS and FAIL, name of deversioned ntlist file for PASS/FAIL seqs
  my %sub_fasta_file_H = (); # two keys: PASS and FAIL, name of sub fasta file for PASS/FAIL seqs 
  my %cur_nseq_H       = (); # two keys: PASS and FAIL, number of seqs that PASS/FAIL

  my @reordered_ref_list_seqname_A = (); # put RefSeqs we annotate for first
  foreach $ref_list_seqname (@ref_list_seqname_A) {
    if($annotate_ref_list_seqname_H{$ref_list_seqname}) { push(@reordered_ref_list_seqname_A, $ref_list_seqname); }
  }
  foreach $ref_list_seqname (@ref_list_seqname_A) {
    if(! $annotate_ref_list_seqname_H{$ref_list_seqname}) { push(@reordered_ref_list_seqname_A, $ref_list_seqname); }
  }

  my $total_assigned_pass   = 0; # total number of assigned seqs that PASS
  my $total_assigned_fail   = 0; # total number of assigned seqs that FAIL
  my $annot_assigned_pass   = 0; # total number of assigned seqs that PASS that will be annotated
  my $annot_assigned_fail   = 0; # total number of assigned seqs that FAIL that will be annotated
  my $noannot_assigned_pass = 0; # total number of assigned seqs that PASS that will not be annotated
  my $noannot_assigned_fail = 0; # total number of assigned seqs that FAIL that will not be annotated

  foreach my $ref_list_seqname (@reordered_ref_list_seqname_A) {
    $ntlist_file_H{"PASS"}    = $out_root . ".$ref_list_seqname.cp.ntlist";
    $dv_ntlist_file_H{"PASS"} = $out_root . ".$ref_list_seqname.cp.deversioned.ntlist";
    $sub_fasta_file_H{"PASS"} = $out_root . ".$ref_list_seqname.cp.fa";
    $ntlist_file_H{"FAIL"}    = $out_root . ".$ref_list_seqname.cf.ntlist";
    $dv_ntlist_file_H{"FAIL"} = $out_root . ".$ref_list_seqname.cf.deversion.ntlist";
    $sub_fasta_file_H{"FAIL"} = $out_root . ".$ref_list_seqname.cf.fa";
    $cur_nseq_H{"PASS"} = 0;
    $cur_nseq_H{"FAIL"} = 0;
    my $cur_list_seq;
    my $cur_fasta_seq;

    for(my $z = 0; $z < scalar(@{$ntlist_HA{$ref_list_seqname}}); $z++) { 
      $cur_list_seq   = $ntlist_HA{$ref_list_seqname}[$z];
      $cur_fasta_seq  = $cls_list2fasta_seqname_H{$cur_list_seq};
      if(! exists $pass_fail_H{$cur_fasta_seq}) { 
        DNAORG_FAIL("ERROR, pass_fail_H{$cur_fasta_seq} does not exist", 1, $ofile_info_HH{"FH"});
      }
      if(($pass_fail_H{$cur_fasta_seq} ne "PASS") && ($pass_fail_H{$cur_fasta_seq} ne "FAIL")) { 
        DNAORG_FAIL("ERROR, pass_fail_H{$cur_fasta_seq} does not equal PASS or FAIL but $pass_fail_H{$cur_fasta_seq}", 1, $ofile_info_HH{"FH"});
      }
      $cur_nseq_H{$pass_fail_H{$cur_fasta_seq}}++; 
    }

    my $annotate_field = "";
    if($do_annotate) { 
      $annotate_field = sprintf("    %11s", ($annotate_ref_list_seqname_H{$ref_list_seqname} ? "yes" : "no"));
    }
    push(@tmp_output_A, sprintf("%-12s  %10d  %10.4f    %10d  %10.4f%s\n", $ref_list_seqname, $cur_nseq_H{"PASS"}, $cur_nseq_H{"PASS"} / $n_cls_seqs, $cur_nseq_H{"FAIL"}, $cur_nseq_H{"FAIL"} / $n_cls_seqs, $annotate_field));

    $total_assigned_pass += $cur_nseq_H{"PASS"};
    $total_assigned_fail += $cur_nseq_H{"FAIL"};
    if($annotate_ref_list_seqname_H{$ref_list_seqname}) { 
      $annot_assigned_pass += $cur_nseq_H{"PASS"};
      $annot_assigned_fail += $cur_nseq_H{"FAIL"};
    }
    else { 
      $noannot_assigned_pass += $cur_nseq_H{"PASS"};
      $noannot_assigned_fail += $cur_nseq_H{"FAIL"};
    }

    # create files for sequences that PASS and FAIL
    foreach my $class ("PASS", "FAIL") { 
      if($cur_nseq_H{$class} > 0) { 
        if($annotate_ref_list_seqname_H{$ref_list_seqname} == 1) { $nseq_above_zero++; }
        open(NTLIST, "> $ntlist_file_H{$class}")  || fileOpenFailure($ntlist_file_H{$class}, $0, $!, "writing", $ofile_info_HH{"FH"});
        if(! $infasta_mode) { # only print ref accession if ! --infasta
          print NTLIST "$ref_list_seqname\n";
        }
        foreach $cur_list_seq (@{$ntlist_HA{$ref_list_seqname}}) {
          $cur_fasta_seq  = $cls_list2fasta_seqname_H{$cur_list_seq};
          if($pass_fail_H{$cur_fasta_seq} eq $class) { 
            print NTLIST "$cur_list_seq\n";
          }
        }
        close(NTLIST);
        addClosedFileToOutputInfo(\%ofile_info_HH, "$ref_list_seqname.$class.ntlist", $ntlist_file_H{$class}, 1, sprintf("ntlist with %sing seqs for $ref_list_seqname", $class));

        # if --infasta mode make another ntlist file without versions, and also use the original ntfiles file to fetch the sequences
        if($infasta_mode) { 
          open(DVNTLIST, "> $dv_ntlist_file_H{$class}")  || fileOpenFailure($dv_ntlist_file_H{$class}, $0, $!, "writing", $ofile_info_HH{"FH"});
          print DVNTLIST "$ref_list_seqname\n";
          foreach $cur_list_seq (@{$ntlist_HA{$ref_list_seqname}}) {
            $cur_fasta_seq  = $cls_list2fasta_seqname_H{$cur_list_seq};
            if($pass_fail_H{$cur_fasta_seq} eq $class) { 
              print DVNTLIST fetchedNameToListName($cur_list_seq) . "\n";
            }
          }
          close(DVNTLIST);
          addClosedFileToOutputInfo(\%ofile_info_HH, "$ref_list_seqname.$class.deversioned.ntlist", $dv_ntlist_file_H{$class}, 1, sprintf("De-versioned ntlist of %sing sequences for $ref_list_seqname", $class));
                                    
          # fetch the sequences into a new fasta file
          sleep(0.1); # make sure that NTLIST is closed
          $cmd  = "cat $ntlist_file_H{$class} |" . $execs_H{"esl-sfetch"} . " -f $cls_fa - > $sub_fasta_file_H{$class}"; # fetch the sequences
          runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
          addClosedFileToOutputInfo(\%ofile_info_HH, "$ref_list_seqname.$class.fa", $sub_fasta_file_H{$class}, 1, sprintf("Fasta file with %sing sequences assigned to $ref_list_seqname", $class));
        }
      }
    }
  }
  
  # Generate file that lists non-assigned sequences
  my $non_assigned_file = $out_root . ".non-assigned";
  
  my $cur_nseq = scalar(@{$ntlist_HA{"non-assigned"}});
  push(@tmp_output_A, sprintf("%-12s  %10d  %10.4f    %10d  %10.4f%s\n", "NON-ASSIGNED", $cur_nseq, $cur_nseq / $n_cls_seqs, $cur_nseq, $cur_nseq / $n_cls_seqs, ($do_annotate) ? "             no" : ""));
  if($cur_nseq > 0) { 
    open(NALIST, "> $non_assigned_file")  || fileOpenFailure($non_assigned_file, $0, $!, "writing", $ofile_info_HH{"FH"});
    foreach (@{$ntlist_HA{"non-assigned"}}) {
      print NALIST "$_\n";
    }
    addClosedFileToOutputInfo(\%ofile_info_HH, "non-assigned", $non_assigned_file, 1, "List of sequences not assigned to a RefSeq");
  }
  push(@tmp_output_A, sprintf("%-12s  %10s  %10s    %10s  %10s%s\n", "#-----------", "----------", "----------", "----------", "----------", ($do_annotate) ? "    -----------" : ""));
  if($do_annotate) { 
    push(@tmp_output_A, sprintf("%-12s  %10d  %10.4f    %10d  %10.4f%s\n", "SUM-ASSIGNED", $annot_assigned_pass, $annot_assigned_pass / $n_cls_seqs, $annot_assigned_fail, $annot_assigned_fail / $n_cls_seqs, ($do_annotate) ? "            yes" : ""));
    push(@tmp_output_A, sprintf("%-12s  %10d  %10.4f    %10d  %10.4f%s\n", "SUM-ASSIGNED", $noannot_assigned_pass, $noannot_assigned_pass / $n_cls_seqs, $noannot_assigned_fail, $noannot_assigned_fail / $n_cls_seqs, ($do_annotate) ? "             no" : ""));
  }
  push(@tmp_output_A, sprintf("%-12s  %10d  %10.4f    %10d  %10.4f%s\n", "SUM-ASSIGNED", $total_assigned_pass, $total_assigned_pass / $n_cls_seqs, $total_assigned_fail, $total_assigned_fail / $n_cls_seqs, ($do_annotate) ? "           both" : ""));
  

  # add $ref_list and $cls_list to output direcotry
  my $out_ref_list = $out_root . ".all.refseqs";
  my $out_cls_list = $out_root . ".all.seqs";
  $cmd = "cat $ref_list | grep . > $out_ref_list";
  runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
  $cmd = "cat $cls_list | grep . > $out_cls_list";
  runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
  
  addClosedFileToOutputInfo(\%ofile_info_HH, "RefSeqs", $out_ref_list, 1, "List of RefSeqs in the HMM library");
  addClosedFileToOutputInfo(\%ofile_info_HH, "ClsSeqs", $out_cls_list, 1, "List of sequences that were sorted into ntlists");
  
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  
  # now output statistics summary
  foreach my $tmp_line (@tmp_output_A) { 
    outputString($log_FH, 1, $tmp_line);
  }
##########################################################################################################################################################

#################################################################################
# If -A, run dnaorg_annotate.pl for each model that had >0 seqs classified to it,
# if not, 
#################################################################################
  if($do_annotate) { 
    outputString($log_FH, 1, "#\n#\n");
    $start_secs = outputProgressPrior("Running dnaorg_annotate.pl $nseq_above_zero time(s) to annotate sequences", $progress_w, $log_FH, *STDOUT);
    outputString($log_FH, 1, "\n");
    my $annotate_cmd = "";
    my $annotate_cons_opts = "";
    my $cur_out_dir  = "";
    my $cur_out_root = "";
    my $ctr = 1;
    my %filekey_H = ();
    $filekey_H{"PASS"} = "cp";
    $filekey_H{"FAIL"} = "cf";
    # for each family with >0 
    foreach my $ref_list_seqname (@ref_list_seqname_A) {
      if($annotate_ref_list_seqname_H{$ref_list_seqname} == 1) { 
        foreach my $class ("PASS", "FAIL") { 
          my $classkey = $filekey_H{$class};
          my $ntlist_file    = $out_root . ".$ref_list_seqname.$classkey.ntlist";
          my $dv_ntlist_file = $out_root . ".$ref_list_seqname.$classkey.deversioned.ntlist"; # same as $ntlist_file but accessions have been de-versioned
          my $sub_fasta_file = $out_root . ".$ref_list_seqname.$classkey.fa";
          $build_dir = $annotate_dir . "/" . $ref_list_seqname; # do not include $classkey here 
          $cur_out_root = $dir_tail . "-" . $ref_list_seqname  . ".$classkey";
          $cur_out_dir  = $dir . "/" . $cur_out_root;
          $cur_nseq = 0;
          for(my $z = 0; $z < scalar(@{$ntlist_HA{$ref_list_seqname}}); $z++) { 
            if($pass_fail_H{$ntlist_HA{$ref_list_seqname}[$z]} eq $class) { 
              $cur_nseq++; 
            }
          }
          if($cur_nseq > 0) { 
            $annotate_cons_opts = build_opts_hash_to_opts_string(\%{$buildopts_used_HH{$ref_list_seqname}});
            $annotate_cmd = $execs_H{"dnaorg_annotate"} . " " . $annotate_cons_opts . " --dirbuild $build_dir --dirout $cur_out_dir";
            if(defined $annotate_non_cons_opts) { $annotate_cmd .= " " . $annotate_non_cons_opts; }
            if($do_keep)    { $annotate_cmd .= " --keep"; }
            if($be_verbose) { $annotate_cmd .= " -v"; }
            if($infasta_mode) { 
              $annotate_cmd .= " --infasta $sub_fasta_file --refaccn $ref_list_seqname";
            }
            else { # not fasta mode
              $annotate_cmd .= " $ntlist_file";
            }
            runCommand($annotate_cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
            # now copy the sequin feature tables to this top level directory:
            my $src_pass_sqtable  = $cur_out_dir . "/" . $cur_out_root . ".dnaorg_annotate.ap.sqtable";
            my $src_fail_sqtable  = $cur_out_dir . "/" . $cur_out_root . ".dnaorg_annotate.af.sqtable";
            my $src_long_sqtable  = $cur_out_dir . "/" . $cur_out_root . ".dnaorg_annotate.long.af.sqtable";
            my $src_pass_list     = $cur_out_dir . "/" . $cur_out_root . ".dnaorg_annotate.ap.list";
            my $src_fail_list     = $cur_out_dir . "/" . $cur_out_root . ".dnaorg_annotate.af.list";
            my $dest_pass_sqtable = $dir . "/" . $cur_out_root . ".dnaorg_annotate.ap.sqtable";
            my $dest_fail_sqtable = $dir . "/" . $cur_out_root . ".dnaorg_annotate.af.sqtable";
            my $dest_long_sqtable = $dir . "/" . $cur_out_root . ".dnaorg_annotate.long.af.sqtable";
            my $dest_pass_list    = $dir . "/" . $cur_out_root . ".dnaorg_annotate.ap.list";
            my $dest_fail_list    = $dir . "/" . $cur_out_root . ".dnaorg_annotate.af.list";
            runCommand("cp $src_pass_sqtable $dest_pass_sqtable", opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
            runCommand("cp $src_fail_sqtable $dest_fail_sqtable", opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
            runCommand("cp $src_long_sqtable $dest_long_sqtable", opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
            runCommand("cp $src_pass_list    $dest_pass_list",    opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
            runCommand("cp $src_fail_list    $dest_fail_list",    opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
            addClosedFileToOutputInfo(\%ofile_info_HH, "pass-sqtbl" . $ctr++, $dest_pass_sqtable, 1, "annotation results for $ref_list_seqname sequences that pass dnaorg_annotate.pl");
            addClosedFileToOutputInfo(\%ofile_info_HH, "fail-sqtbl" . $ctr++, $dest_fail_sqtable, 1, "annotation results for $ref_list_seqname sequences that fail dnaorg_annotate.pl (minimal)");
            addClosedFileToOutputInfo(\%ofile_info_HH, "longsqtbl"  . $ctr++, $dest_long_sqtable, 1, "annotation results for $ref_list_seqname sequences that fail dnaorg_annotate.pl (verbose)");
            addClosedFileToOutputInfo(\%ofile_info_HH, "pass-list"  . $ctr++, $dest_pass_sqtable, 1, "list of $ref_list_seqname sequences that pass dnaorg_annotate.pl");
            addClosedFileToOutputInfo(\%ofile_info_HH, "fail-list"  . $ctr++, $dest_fail_sqtable, 1, "list of $ref_list_seqname sequences that fail dnaorg_annotate.pl (minimal)");

            # and finally, determine the number of passing sequences, by counting the 
            # number of lines in the pass_list file
            my $npass = countLinesInFile($dest_pass_list, $ofile_info_HH{"FH"});
            push(@annotate_summary_output_A, sprintf("%-12s  %6s  %10d    %10d  %10.4f    %10d  %10.4f\n", $ref_list_seqname, "C-" . $class,  $cur_nseq, $npass, $npass / $cur_nseq, $cur_nseq - $npass, ($cur_nseq - $npass) / $cur_nseq));
          }
        }
      }
    }
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
    
    my @output_A = ();
    push(@output_A, sprintf("#\n"));
    push(@output_A, sprintf("# Number of annotated sequences that PASSed/FAILed dnaorg_annotate.pl:\n"));
    push(@output_A, sprintf("#\n"));
    push(@output_A, sprintf("%-12s  %6s  %10s    %10s  %10s    %10s  %10s\n", "# RefSeq-acc", "C-p/f", "num-annot",  "num-A-PASS",   "frc-A-PASS", "num-A-FAIL",  "frc-A-FAIL"));
    push(@output_A, sprintf("%-12s  %6s  %10s    %10s  %10s    %10s  %10s\n", "#-----------", "------", "----------", "----------", "----------", "----------", "----------"));
    # output annotation summary
    foreach my $out_line (@output_A, @annotate_summary_output_A) { 
      outputString($log_FH, 1, $out_line);
    }
  }
  else { 
    if($inlist_mode) { 
      outputString($log_FH, 1, "#\n");
      outputString($log_FH, 1, "# *** The $out_root.dnaorg_classify.<s>.ntlist files\n");
      outputString($log_FH, 1, "# *** can be used as input to dnaorg_annotate.pl once you've run 'dnaorg_build.pl <s>'\n");
      outputString($log_FH, 1, "# *** to create models for RefSeq <s>.\n#\n");
    }
    elsif($infasta_mode) { 
      outputString($log_FH, 1, "#\n");
      outputString($log_FH, 1, "# *** The $out_root.dnaorg_classify.<s>.deversioned.ntlist files or\n");
      outputString($log_FH, 1, "# *** the $out_root.dnaorg_classify.<s>.fa files (with the --infasta option) can be\n");
      outputString($log_FH, 1, "# *** used as input to dnaorg_annotate.pl once you've run 'dnaorg_build.pl <s>'\n");
      outputString($log_FH, 1, "# *** to create models for RefSeq <s>.\n#\n");
    }
  }
} # end of 'else' entered if $onlybuild_mode is FALSE

##########
# Conclude
##########

# remove temp files, unless --keep
if(! opt_Get("--keep", \%opt_HH)) { 
  foreach my $file2rm (@files2rm_A) {
    if(-e $file2rm) { 
      removeFileUsingSystemRm($file2rm, "dnaorg_classify.pl::main", \%opt_HH, $ofile_info_HH{"FH"});
    }
  }
}

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir, \%ofile_info_HH);
exit 0;

#################################################################################
# SUBROUTINES:
#
# list_to_fasta():             given a list file, fetch sequences into a fasta file 
#                              using edirect
#
# fasta_to_list_and_lengths(): given a fasta file, create a list of sequences in it
#                              and get the sequence names and lengths to data structures
#
#
# Old subroutines no longer used but kept for reference:
#
# createFastas():              given a list of sequences create a single sequence fasta file
#                              for each sequence. Replaced with list_to_fasta() which creates
#                              one large fasta file.
#
# nhmmscanSeqs():              Calls nhmmscan for each sequence in input list. Replaced with
#                              dnaorg.pm::cmscanOrNhmmscanWrapper() a function also called 
#                              by dnaorg_annotate.pl which splits up the sequence file into
#                              an appropriate number of files (with >= 1 sequence), runs
#                              each job on the compute farm, and monitors the output.
#
#################################################################################
#
# Sub name:  list_to_fasta()
#
# Authors:    Lara Shonkwiler and Eric Nawrocki
# Date:       2016 Aug 01 [Lara] 2017 Nov 21 [Eric]
#
# Purpose:   Given a list of accessions, create a fasta file with all of them.
#            using efetch.
#
# Arguments: $list_file    list of all accessions for fastas to be made for
#            $fa_file      directory path for fasta file to create
#            $opt_HHR      reference to hash of options
#            $FH_HR        reference to hash of file handles
#
# Returns:   void
#
#################################################################################
sub list_to_fasta {
  my $sub_name  = "list_to_fasta()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }
  
  my ($list_file, $fa_file, $opt_HHR, $FH_HR) = (@_);
  
  $cmd = "cat $list_file | epost -db nuccore -format acc | efetch -format fasta > $fa_file";
  runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);

  return;
}

#################################################################################
#
# Sub name:  fasta_to_list_and_lengths()
#
# Author:     Eric Nawrocki
# Date:       2017 Nov 21
#
# Purpose:   Given a fasta file, create a list file with names of the sequences 
#            in that fasta file, and store the lengths of each sequence in
#            %{$seqlen_HR}.
#
# Arguments: $fa_file      directory path for fasta file to create
#            $list_file    list of all accessions for fastas to be made for, can be undef if undesired
#            $tmp_file     temporary file to print sequence names and lengths too
#            $seqstat      path to esl-seqstat executable
#            $seqname_AR:  ref to array of sequence names to fill, can be undef if undesired
#            $seqlen_HR:   ref to hash of sequence lengths to fill, can be undef if undesired
#            $opt_HHR      reference to hash of options
#            $FH_HR        reference to hash of file handles
#
# Returns:   void
#
#################################################################################
sub fasta_to_list_and_lengths { 
  my $sub_name  = "fasta_to_list()";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }
  
  my ($fa_file, $list_file, $tmp_file, $seqstat, $seqname_AR, $seqlen_HR, $opt_HHR, $FH_HR) = (@_);

  $cmd = $seqstat . " --dna -a $fa_file | grep ^\= | awk '{ printf(\"\%s \%s\\n\", \$2, \$3); }' > $tmp_file";
  runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);

  open(IN, $tmp_file) || fileOpenFailure($tmp_file, $0, $!, "reading", $FH_HR);
  if(defined $list_file) { 
    open(OUT, ">", $list_file) || fileOpenFailure($tmp_file, $0, $!, "reading", $FH_HR);
  }  
  while(my $line = <IN>) { 
    chomp $line;
    #print ("in $sub_name, input line $line\n");
    my ($seqname, $seqlength) = split(/\s+/, $line);
    if(defined $list_file) { 
      print OUT $seqname . "\n";
    }
    if(defined $seqname_AR) { 
      push(@{$seqname_AR}, $seqname); 
    }
    if(defined $seqlen_HR) { 
      #printf("\tsetting seqlen_HR->{$seqname} to $seqlength\n");
      $seqlen_HR->{$seqname} = $seqlength;
    }
  }
  close(IN);
  if(defined $list_file) { 
    close(OUT);
  }

  if((-e $tmp_file) && (! opt_Get("--keep", $opt_HHR))) { 
    removeFileUsingSystemRm($tmp_file, $sub_name, $opt_HHR, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: validate_build_dir()
# Incept:     EPN, Mon Feb 12 12:30:44 2018
#
# Purpose:   Validate that a build directory exists and has the files
#            that we need to annotate with. This actually can't check
#            that everything a downstream dnaorg_annotate.pl run needs
#            in a build dir (e.g. number and names of models) is kosher, 
#            but it can check for a few things that are required:
#            directory exists and contains a .consopts file.
#
# Arguments:
#  $build_dir:         the build directory we are validating
#  $buildopts_used_HR: REF to build options used when creating this build directory 
#                      with dnaorg_build.pl, key is option (e.g. --matpept), value is
#                      option argument, "" for none.
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:    REF to output file hash
# 
# Returns:  void
# 
# Dies: If $build_dir doesn't exist or .consopts file does not
#       exist or is unreadable or is in a bad format.
#       
#################################################################
sub validate_build_dir { 
  my $sub_name  = "validate_build_dir()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }
  
  my ($build_dir, $buildopts_used_HR, $opt_HHR, $ofile_info_HHR) = (@_);

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience
  
  if(! -d $build_dir) { 
    DNAORG_FAIL(sprintf("ERROR, the directory $build_dir should exist based on your -A option argument %s, but it does not.\nDid you run dnaorg_build.pl?\n", opt_Get("-A", $opt_HHR)), 1, $FH_HR);
  }

  my $build_root         = removeDirPath($build_dir);
  my $build_dir_and_root = $build_dir . "/" . $build_root;
  my $consopts_file      = $build_dir_and_root . ".dnaorg_build.consopts";
  if(! -e $consopts_file) { 
    DNAORG_FAIL(sprintf("ERROR, required file $consopts_file should exist based on your -A option argument %s, but it does not.\nDid you run dnaorg_build.pl?\n", opt_Get("-A", $opt_HHR)), 1, $FH_HR);
  }

  # parse the consopts file
  my %consopts_notused_H = (); # key option in consopts file, value argument used in dnaorg_build.pl
  my %consmd5_H  = ();         # key option used in dnaorg_build.pl in consopts file, md5 checksum value of the file name argument used in dnaorg_build.pl
  # we pass in $buildopts_used_HR, because caller needs that, but not %consopts_notused_H nor %consmd5_H
  parseConsOptsFile($consopts_file, $buildopts_used_HR, \%consopts_notused_H, \%consmd5_H, $FH_HR);

  # for any option that takes a file name as its argument, 
  # make sure that that file exists and its md5 matches what
  # is expected in the build dir
  my $opt;
  my $optfile;
  my $optfile_md5;
  my $optroot;
  foreach my $opt (sort keys (%consmd5_H)) { 
    if($consmd5_H{$opt} ne "") { 
      $optroot = $opt;
      $optroot =~ s/^\-*//; # remove trailing '-' values
      $optfile = $build_dir_and_root . ".dnaorg_build." . $optroot;
      if(! -s $optfile) { 
        DNAORG_FAIL("ERROR, the file $optfile required to run dnaorg_annotate.pl later for build directory $build_root does not exist.", 1, $FH_HR);
      }
      $optfile_md5 = md5ChecksumOfFile($optfile, $sub_name, $opt_HHR, $FH_HR);
      if($consmd5_H{$opt} ne $optfile_md5) { 
        DNAORG_FAIL("ERROR, the file $optfile required to run dnaorg_annotate.pl later for build directory $build_root\ndoes not appear to be identical to the file used with dnaorg_build.pl.\nThe md5 checksums of the two files differ: dnaorg_build.pl: " . $consmd5_H{$opt} . " $optfile: " . $optfile_md5, 1, $FH_HR);
      }
      # update the $buildopts_used_HR to use the local file
      $buildopts_used_HR->{$opt} = $optfile;
    }
  }
          
  return;
}

#################################################################
# Subroutine: build_opts_hash_to_opts_string()
# Incept:     EPN, Mon Feb 12 15:10:19 2018
#
# Purpose:   Given a hash of options, Validate that a build directory exists and has the files
#            that we need to annotate with. This actually can't check
#            that everything a downstream dnaorg_annotate.pl run needs
#            in a build dir (e.g. number and names of models) is kosher, 
#            but it can check for a few things that are required:
#            directory exists and contains a .consopts file.
#
# Arguments:
#  $buildopts_used_HR: REF to build options used when creating this build directory 
#                      with dnaorg_build.pl, key is option (e.g. --matpept), value is
#                      option argument, "" for none.
# 
# Returns:  string of options
# 
# Dies: Never
#       
#################################################################
sub build_opts_hash_to_opts_string { 
  my $sub_name  = "build_opts_hash_to_opts_string";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }
  
  my ($buildopts_used_HR) = (@_);

  my $ret_str = "";
  foreach my $opt (sort keys (%{$buildopts_used_HR})) { 
    if($ret_str ne "") { $ret_str .= " "; }
    $ret_str .= $opt;
    if($buildopts_used_HR->{$opt} ne "") { 
      $ret_str .= " " . $buildopts_used_HR->{$opt};
    }
  }

  return $ret_str;
}
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# OLD SUBROUTINES NO LONGER USED, KEPT HERE FOR REFERENCE
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# 
# Sub name:  createFastas()
# EPN, Wed Nov 22 15:09:41 2017: NO LONGER USED, list_to_fasta creates one
# large fasta file of all seqs, then it's split as necessary by other functions.
# 
# Author:    Lara Shonkwiler
# Date:      2016 Aug 01
#
# Purpose:   Given a list of accessions, creates a fasta file for each one
# 
# Arguments: $accn_list    list of all accessions for fastas to be made for
#            $out_root     directory path
#            $opt_HHR      reference to hash of options
#            $FH_HR        reference to hash of file handles
#
#
# Returns:   void
#
#################################################################################
sub createFastas {
    my $sub_name  = "createFastas()";
    my $nargs_expected = 4;
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

    my ($accn_list, $out_root, $opt_HHR, $FH_HR) = (@_);

    # create a file containing the concatenated fasta files for all refseqs
    my $all_refseq_fastas = $out_root . ".all.fa";
    $cmd = "cat $accn_list | epost -db nuccore -format acc | efetch -format fasta > $all_refseq_fastas";
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);

    # separate that file into individual fasta files
    # This section contains code adopted from break_fasta.pl (created by Alejandro Schaffer)
    # This section also contains code adopted from entrez_names.pl (created by Lara Shonkwiler)
    my $infile; #input FASTA file with many sequences                                                                                        
    my $outfile; #output FASTA file with one sequence                                                                                                                                                                 
    my $nofile = 1;
    my $nextline;
    my $same_sequence = 2;
    my $state = $nofile;
    my $id;
    my $new_file_name;

    open(SEQS, $all_refseq_fastas) or die "Cannot open $all_refseq_fastas\n";
    while(defined($nextline = <SEQS>)) {
	chomp($nextline);
	if ($nofile == $state) {
	    #$id = $nextline;
	    ($id) = ($nextline =~ m/^>(\S+)/);               # get rid of > and words after space                                                                                                                              
	    $id =~ s/^gi\|?\d+\|\w+\|//; # get rid of everything before the accession number                                                                                                                  
            # special case to deal with accessions that being with pdb\|, e.g. 'pdb|5TSN|T' which is the 
            # name of the sequence that gets fetched for the accession 5TSN|T
            if($id =~ /^pdb\|(\S+)\|(\S+)$/) { 
              $id = $1 . "_" . $2;
            }
	    else { 
              $id =~ s/\|$//;               # get rid of end | or                                                                                                                                  
              $id =~ s/\|\w+//;             # get rid of end | and everything after |
            }        
	    # gets rid of version number
            $id =~ s/\.\d+$//;

	    $new_file_name = $out_root . "." . $id . "." . "fasta";
	    open(ONEFILE, ">$new_file_name") or die "Cannot open $new_file_name\n";
	    $state = $same_sequence;
	    print ONEFILE "$nextline";
	    print ONEFILE "\n";
	}
	else {
	    if ($same_sequence == $state) {
		if ($nextline =~ m/>/ ) {
		    close(ONEFILE);
		    #$id = $nextline;
		    ($id) = ($nextline =~ m/^>(\S+)/);               # get rid of > and words after space                                                                                                                              
		    $id =~ s/^gi\|?\d+\|\w+\|//; # get rid of everything before the accession number                                                                                                     
                    # special case to deal with accessions that being with pdb\|, e.g. 'pdb|5TSN|T' which is the 
                    # name of the sequence that gets fetched for the accession 5TSN|T
                    if($id =~ /^pdb\|(\S+)\|(\S+)$/) { 
                      $id = $1 . "_" . $2;
                    }
                    else { 
                      $id =~ s/\|$//;               # get rid of end | or                                                                                                                                  
                      $id =~ s/\|\w+//;             # get rid of end | and everything after |
                    }        
		    
		    # gets rid of version number                                                                                                                                                                                                         
		    if($id =~ m/\.\d+$/) {
		     $id =~ s/\.\d+$//;
		    }
 
		    $new_file_name = $out_root . "." . $id . "." . "fasta";
		    open(ONEFILE, ">$new_file_name") or die "Cannot open $new_file_name\n";
		    $state = $same_sequence;
		}
		print ONEFILE "$nextline";
		print ONEFILE "\n";
	    }
	}
    }
    close(ONEFILE);
    close(SEQS);

    # get rid of concatenated fasta file
    $cmd = "rm $all_refseq_fastas";
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);
}


#######################################################################################
#
# Sub name:  nhmmscanSeqs()
#
# EPN, Wed Nov 22 15:08:08 2017: NO LONGER USED, dnaorg.pm:cmscanOrNhmmscanWrapper
#            handles this: splitting up sequence file as necessary, submitting
#            jobs (if necessary) and monitoring the farm until they're done.
#
# Author:    Lara Shonkwiler (adopted from code by Alejandro Schaffer)
# Date:      2016 Aug 01
#
# Purpose:   Calls nhmmscan for each seq against the given HMM library
# 
# Arguments: $nhmmscan            - path to hmmer-3.1b2 nhmmscan executable
#            $out_root            - path to the output directory
#            $cls_list_seqname_AR   - unique identifier that distinguishes this command from others
#            $ref_library         - HMM library of RefSeqs
#            $files2rm_AR         - REF to a list of files to remove at the end of the
#                                   program
#            $ofile_info_HHR      - REF to the output file info hash
#            
#            
#
# Returns:   void
#
#######################################################################################
sub nhmmscanSeqs {
  my $sub_name  = "nhmmscanSeqs()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }
  
  my ($nhmmscan, $out_root, $cls_list_seqname_AR, $ref_library, $files2rm_AR, $ofile_info_HHR) = (@_);
  
  my $qsub_script_name = $out_root . ".qsub.script";                                                                                                                                         
  open(QSUB, ">$qsub_script_name") || fileOpenFailure($qsub_script_name, $sub_name, $!, "writing", $ofile_info_HHR->{"FH"});
  foreach my $seq (@{$cls_list_seqname_AR}) {
    my $seq_tbl_file = $out_root . ".$seq.tblout";
    my $seq_results_file = $out_root . ".$seq.nhmmscan.out";
    my $fa_file = $out_root . ".$seq.fasta";
    
    # --noali , --cpu 0
    my $cmd = "$nhmmscan --noali --cpu 0 --tblout $seq_tbl_file $ref_library $fa_file > $seq_results_file";
    
    my $output_script_name = "$out_root" . "." . "$seq" . "\.qsub\.csh";                                                                                                                               
    open(SCRIPT, ">$output_script_name") or die "Cannot open 4 $output_script_name\n";
    print QSUB "chmod +x $output_script_name\n";
    print QSUB "qsub $output_script_name\n";
    print SCRIPT "source /etc/profile\n";
    print SCRIPT "\#!/bin/tcsh\n";
    print SCRIPT "\#\$ -P unified\n";
    print SCRIPT "\n";
    print SCRIPT "\# list resource request options\n";
    print SCRIPT "\#\$  -l h_rt=288000,h_vmem=32G,mem_free=32G,reserve_mem=32G,m_mem_free=32G  \n";
    # old version of above line: -l h_vmem=32G,reserve_mem=32G,mem_free=32G
    # new version of above line: -l h_rt=288000,h_vmem=32G,mem_free=32G,reserve_mem=32G,m_mem_free=32G
    print SCRIPT "\n";
    print SCRIPT "\# split stdout and stderr files (default is they are joined into one file)\n";
    print SCRIPT "\#\$ -j n\n";
    print SCRIPT "\n";
    print SCRIPT "\#define stderr file\n";
    my $error_file_name = "$out_root." . "$seq" . "\.qsub\.err";                                                                                                                               
    print SCRIPT "\#\$ -e $error_file_name\n";
    print SCRIPT "\# define stdout file\n";
    my $diagnostic_file_name = "$out_root" . "." . "$seq" . "\.qsub\.out";                                                                                                                               
    print SCRIPT "\#\$ -o $diagnostic_file_name\n";
    print SCRIPT "\n";
    print SCRIPT "\# job is re-runnable if SGE fails while it's running (e.g. the host reboots)\n";
    print SCRIPT "\#\$ -r y\n";
    print SCRIPT "\# stop email from being sent at the end of the job\n";
    print SCRIPT "\#\$ -m n\n";
    print SCRIPT "\n";
    print SCRIPT "\# trigger NCBI facilities so runtime enviroment is similar to login environment\n";
    print SCRIPT "\#\$ -v SGE_FACILITIES\n";
    print SCRIPT "\n";
    print SCRIPT "echo \"Running qsub\"\n";                                                                                                                               
    print SCRIPT "\n";
    print SCRIPT "$cmd";                                                                                                                               
    print SCRIPT "\n";
    close(SCRIPT);                                                                                                                               
    system("chmod +x $output_script_name");                                                                                                          
    
    addClosedFileToOutputInfo($ofile_info_HHR, "$seq.tbl", $seq_tbl_file, 1, "Table summarizing nhmmscan results for $seq against $ref_library");
    #addClosedFileToOutputInfo(\%ofile_info_HH, "$seq.results", $seq_results_file, 1, "nhmmscan results for searching $ref_library with $seq");                                         
    
    # clean up leftover fasta files and qsub files
    if(! $do_keep) {
      push(@{$files2rm_AR}, $fa_file);
      push(@{$files2rm_AR}, $seq_results_file);
      push(@{$files2rm_AR}, $output_script_name);
      push(@{$files2rm_AR}, $error_file_name);
      push(@{$files2rm_AR}, $diagnostic_file_name);
    }
    
  }
  
  close(QSUB);                                                                                                                               
  system("chmod +x $qsub_script_name");                                                                                                                               
  system("$qsub_script_name");
  
  # wait until all nhmmscan jobs are done before proceeding
  
#    # only checks last seq - causes null array ref error
#    my $last_seq = @{$cls_list_seqname_AR}[-1];
#    my $last_seq_file = $out_root . "." . $last_seq . ".tblout";
#    my $done = 0; # Boolean - set to false
#    my $last_line;
#
#   until( $done ) {
#	sleep(10);
#
#	if(-s $last_seq_file){
#	    $last_line = `tail -n1 $last_seq_file`;
#	    if( $last_line =~ m/\[ok\]/ ) {
#		$done = 1;
#	    }
#	}
#    }
  
  
#   slow, but accurate version
#    my $last_line;
#    
#  for my $seq (@{$cls_list_seqname_AR}) {
#	my $done = 0; # Boolean - set to false
#	
#	until( $done ) {
#	    sleep(10);
#	    #debug
#	    print "waiting on $seq";
#
#	    my $seq_file = $out_root . "." . $seq . ".tblout";
#	    if(-s $seq_file){
#		$last_line = `tail -n1 $seq_file`;
#		if( $last_line =~ m/\[ok\]/ ) {
#		    $done = 1;
#		}
#	    }
#	}
#    }
  
  
  
  # This one seems to work well
  my $last_line;
  
  my $counter = 0;
  my $end = scalar(@{$cls_list_seqname_AR});
  
  until($counter == $end){
    $counter = 0;
    
    for my $seq (@{$cls_list_seqname_AR}) {
      my $seq_file = $out_root . "." . $seq . ".tblout";
      my $full_seq_file = $out_root . "." . $seq . ".nhmmscan.out";
      if( (-s $seq_file) && (-s $full_seq_file)) {
        $last_line = `tail -n1 $seq_file`;
        if( $last_line =~ m/\[ok\]/ ) {
          $counter++;
        }
      }
    }
    
    sleep(10);
  }
}

