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
use LWP::Simple; 

require "dnaorg.pm"; 
require "epn-options.pm";

#######################################################################################
# What this script does: 
#
# - Preliminaries: 
#   o processes options
#   o creates the output directory
#   o outputs program banner
#   o makes sure the required executables are executable
#
# - Fetches GenBank file (if ! --gb)
# - Parses GenBank file
# - Prunes data read from GenBank file
# - Parses input stockholm file (if --stk)
# - Fills in feature and segment info
# - Translates CDS (if any) and creates BLAST db
# - Builds CM 
# - Writes model info file
# - Writes optional output files
# 
#######################################################################################
# make sure the DNAORGDIR environment variable is set
my $dnaorgdir      = verifyEnvVariableIsValidDir("DNAORGDIR");
my $dnaorgblastdir = verifyEnvVariableIsValidDir("DNAORGBLASTDIR");

my $inf_exec_dir      = $dnaorgdir . "/infernal-dev/src";
my $esl_exec_dir      = $dnaorgdir . "/infernal-dev/easel/miniapps";
my $blast_exec_dir    = $dnaorgblastdir;

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
my $g = 0; # option group

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
$opt_group_desc_H{++$g} = "basic options";
#     option            type       default  group   requires incompat     preamble-output                                                help-output    
opt_Add("-h",           "boolean", 0,           0,    undef, undef,       undef,                                                         "display this help",                                   \%opt_HH, \@opt_order_A);
opt_Add("-g",           "string", 0,           $g,    undef, undef,       "define model group for model info file as <s>",               "define model group for model info file as <s>", \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,          $g,    undef, undef,       "forcing directory overwrite",                                 "force; if dir <output directory> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,          $g,    undef, undef,       "be verbose",                                                  "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--stk",        "string",  undef,      $g,    undef,  undef,      "read single sequence stockholm 'alignment' from <s>",         "read single sequence stockholm 'alignment' from <s>", \%opt_HH, \@opt_order_A);
opt_Add("--gb",         "string",  undef,      $g,    undef,  undef,      "read genbank file from <s>, don't fetch it",                  "read genbank file from <s>, don't fetch it", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,          $g,    undef, undef,       "leave intermediate files on disk",                            "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling what feature types are stored in model info file";
#     option            type       default  group   requires incompat     preamble-output                                                 help-output    
opt_Add("--fall",       "boolean", 0,          $g,    undef,  undef,      "store info for all feature types",                             "store info for all feature types", \%opt_HH, \@opt_order_A);
opt_Add("--fadd",       "string",  undef,      $g,    undef,"--fall",     "also store features types in comma separated string <s>",      "also store feature types in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--fnocds",     "boolean", 0,          $g,    undef,"--fall",     "do not store info for CDS features",                           "do not store info for CDS features", \%opt_HH, \@opt_order_A);
opt_Add("--fnogene",    "boolean", 0,          $g,    undef,"--fall",     "do not store info for gene features",                          "do not store info for gene features", \%opt_HH, \@opt_order_A);
opt_Add("--fnomp",      "boolean", 0,          $g,    undef,"--fall",     "do not store info for mat_peptide features",                   "do not store info for mat_peptide features", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling what qualifiers are stored in model info file";
#     option            type       default  group   requires incompat     preamble-output                                                 help-output    
opt_Add("--qall",       "boolean",  0,        $g,    undef,  undef,       "store info for all qualifiers",                                "store info for all qualifiers", \%opt_HH, \@opt_order_A);
opt_Add("--qadd",       "string",   undef,    $g,    undef,"--qall",      "also store info for qualifiers in comma separated string <s>", "also store info for qualifiers in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--qnoproduct", "boolean",  0,        $g,    undef,"--qall",      "do not store info for product qualifier",                      "do not store info for product qualifier", \%opt_HH, \@opt_order_A);
opt_Add("--qnogene",    "boolean",  0,        $g,    undef,"--qall",      "do not store info for gene qualifier",                         "do not store info for gene qualifier", \%opt_HH, \@opt_order_A);
opt_Add("--qnoexc",     "boolean",  0,        $g,    undef,"--qall",      "do not store info for exception qualifier",                    "do not store info for exception features", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling CDS translation step";
#     option          type       default    group   requires    incompat   preamble-output                                             help-output    
opt_Add("--ttbl",     "integer", 1,            $g,  undef,    "--fnocds",  "use NCBI translation table <n> to translate CDS",          "use NCBI translation table <n> to translate CDS", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling cmbuild step";
#     option          type       default    group   requires    incompat   preamble-output                                             help-output    
opt_Add("--cmn",      "integer", 0,           $g,   undef, "--skipbuild",  "set number of seqs for glocal fwd HMM calibration to <n>", "set number of seqs for glocal fwd HMM calibration to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--cmp7ml",   "boolean", 0,           $g,   undef, "--skipbuild",  "set CM's filter p7 HMM as the ML p7 HMM",                  "set CM's filter p7 HMM as the ML p7 HMM",                  \%opt_HH, \@opt_order_A);
opt_Add("--cmere",    "real",    0,           $g,   undef,  "--skipbuild", "set CM relative entropy target to <x>",                    "set CM relative entropy target to <x>",                    \%opt_HH, \@opt_order_A);
opt_Add("--cmeset",   "real",    0,           $g,   undef,  "--skipbuild", "set CM eff seq # for CM to <x>",                           "set CM eff seq # for CM to <x>",                           \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for skipping stages";
#       option       type       default     group  requires     incompat  preamble-output                                    help-output    
opt_Add("--skipbuild",  "boolean", 0,         $g,    undef,     undef,    "skip the cmbuild step",                           "skip the cmbuild step", \%opt_HH, \@opt_order_A);
opt_Add("--onlyurl",    "boolean", 0,         $g,    undef,"--stk,--gb",  "output genbank file url for accession and exit",  "output genbank file url for accession and exit", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "optional output files";
#       option       type       default     group  requires     incompat  preamble-output                          help-output    
opt_Add("--ftrinfo",    "boolean", 0,         $g,    undef,     undef,    "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);
opt_Add("--sgminfo",    "boolean", 0,         $g,    undef,     undef,    "output internal segment information",   "create file with internal segment information", \%opt_HH, \@opt_order_A);


# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_build.pl [-options] <accession> <path to output directory to create>\n";
my $synopsis = "dnaorg_build.pl :: build homology model of a single sequence for feature annotation";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'g=s'          => \$GetOptions_H{"-g"},
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'stk=s'        => \$GetOptions_H{"--stk"},
                'gb=s'         => \$GetOptions_H{"--gb"},
                'keep'         => \$GetOptions_H{"--keep"},
# options for controlling what feature types are stored in model info file
                'fall'         => \$GetOptions_H{"--fall"},
                'fadd=s'       => \$GetOptions_H{"--fadd"},
                'fnocds'       => \$GetOptions_H{"--fnocds"},
                'fnogene'      => \$GetOptions_H{"--fnogene"},
                'fnomp'        => \$GetOptions_H{"--fnomp"},
# options for controlling what qualifiers are stored in model info file
                'qall'         => \$GetOptions_H{"--qall"},
                'qadd=s'       => \$GetOptions_H{"--qadd"},
                'qnoproduct'   => \$GetOptions_H{"--qnoproduct"},
                'qnogene'      => \$GetOptions_H{"--qnogene"},
                'qnoexc'       => \$GetOptions_H{"--qnoexc"},
# options for controlling CDS translation step
                'ttbl=s'       => \$GetOptions_H{"--ttbl"},
# options for controlling cmbuild step
                'cmn=s'        => \$GetOptions_H{"--cmn"},
                'cmp7ml'       => \$GetOptions_H{"--cmp7ml"},
                'cmere=s'      => \$GetOptions_H{"--cmere"},
                'cmeset=s'     => \$GetOptions_H{"--cmeset"},
# optional for skipping stages
                'skipbuild'    => \$GetOptions_H{"--skipbuild"},
                'onlyurl'      => \$GetOptions_H{"--onlyurl"},
# optional output files
                'sgminfo'      => \$GetOptions_H{"--sgminfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"});

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.45x";
my $releasedate   = "Feb 2019";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date, $dnaorgdir);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 2) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do dnaorg_build.pl -h\n\n";
  exit(1);
}
my ($mdl_name, $dir) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

############################################
# if --onlyurl used, output the url and exit
############################################
if(opt_Get("--onlyurl", \%opt_HH)) { 
  print eutilsFetchUrl($mdl_name, "gb") . "\n";
  exit 0;
}

#############################
# create the output directory
#############################
my $cmd;              # a command to run with runCommand()
my @early_cmd_A = (); # array of commands we run before our log file is opened

if($dir !~ m/\/$/) { $dir =~ s/\/$//; } # remove final '/' if it exists
if(-d $dir) { 
  $cmd = "rm -rf $dir";
  if(opt_Get("-f", \%opt_HH)) { runCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR directory named $dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $dir) { 
  $cmd = "rm $dir";
  if(opt_Get("-f", \%opt_HH)) { runCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR a file named $dir already exists. Remove it, or use -f to overwrite it."; }
}

# create the dir
$cmd = "mkdir $dir";
runCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root = $dir . "/" . $dir_tail . ".dnaorg_build";

#######################
# output program banner
#######################
# output preamble
my @arg_desc_A = ("accession/model name", "output directory");
my @arg_A      = ($mdl_name, $dir);
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

# open the log and command files 
openAndAddFileToOutputInfo(\%ofile_info_HH, "log", $out_root . ".log", 1, "Output printed to screen");
openAndAddFileToOutputInfo(\%ofile_info_HH, "cmd", $out_root . ".cmd", 1, "List of executed commands");
openAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".list", 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
my $FH_HR  = $ofile_info_HH{"FH"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# now we have the log file open, output the banner there too
outputBanner($log_FH, $version, $releasedate, $synopsis, $date, $dnaorgdir);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

#############################################################
# make sure the required executables exist and are executable
#############################################################
my %execs_H = (); # hash with paths to all required executables
if(! opt_Get("--skipbuild", \%opt_HH)) { 
  $execs_H{"cmbuild"}       = $inf_exec_dir . "/cmbuild";
}
$execs_H{"esl-reformat"}  = $esl_exec_dir . "/esl-reformat";
$execs_H{"esl-translate"} = $esl_exec_dir . "/esl-translate";
$execs_H{"makeblastdb"}   = $blast_exec_dir . "/makeblastdb";
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################################
# Fetch the genbank file (if --gb not used)
###########################################
my $progress_w = 60; # the width of the left hand column in our progress output, hard-coded
my $start_secs;
my $gb_file = undef;
if(opt_IsUsed("--gb", \%opt_HH)) { 
  $gb_file = opt_Get("--gb", \%opt_HH);
}
else { 
  # --gb not used, create gb file by fetching using eutils
  $start_secs = outputProgressPrior("Fetching GenBank file", $progress_w, $log_FH, *STDOUT);

  $gb_file = $out_root . ".gb";
  eutilsFetchToFile($gb_file, $mdl_name, "gb", 5, $ofile_info_HH{"FH"});  # number of attempts to fetch to make before dying
  addClosedFileToOutputInfo(\%ofile_info_HH, "gb", $gb_file, 1, "GenBank format file for $mdl_name");

  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

########################
# Parse the genbank file
########################
$start_secs = outputProgressPrior("Parsing GenBank file", $progress_w, $log_FH, *STDOUT);

my %ftr_info_HAH = (); # the feature info 
my %seq_info_HH  = (); # the sequence info 
genbankParse($gb_file, \%seq_info_HH, \%ftr_info_HAH, $FH_HR);
if((! exists $seq_info_HH{$mdl_name}) || (! defined $seq_info_HH{$mdl_name}{"seq"})) { 
  DNAORG_FAIL("ERROR parsing GenBank file $gb_file, did not read sequence for reference accession $mdl_name\n", 1, $FH_HR);
}
if(! exists $ftr_info_HAH{$mdl_name}) { 
  DNAORG_FAIL("ERROR parsing GenBank file $gb_file, did not read info for reference accession $mdl_name\n", 1, $FH_HR);
}

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#######################################################
# Prune data read from %ftr_info_HAH, only keeping what
# we want to output to the eventual model info file
#######################################################
$start_secs = outputProgressPrior("Pruning data read from GenBank file", $progress_w, $log_FH, *STDOUT);

my $ftr_idx;
my $key;
# determine what types of features we will store based on cmdline options
# --fall is incompatible with all other --f* options
my %ftype_H = ();
if(! opt_Get("--fnocds",  \%opt_HH)) { $ftype_H{"CDS"}         = 1; }
if(! opt_Get("--fnogene", \%opt_HH)) { $ftype_H{"gene"}        = 1; }
if(! opt_Get("--fnomp",   \%opt_HH)) { $ftype_H{"mat_peptide"} = 1; }
if(opt_IsUsed("--fadd", \%opt_HH)) { 
  my @fadd_A = split(",", opt_Get("--fadd", \%opt_HH));
  foreach my $f (@fadd_A) { $ftype_H{$f} = 1; }
}

# determine what qualifiers we will store based on cmdline options
# --qall is incompatible with all other --q* options
my %qual_H = ();
$qual_H{"type"}         = 1;
$qual_H{"location"}     = 1;
if(! opt_Get("--qnoproduct", \%opt_HH)) { $qual_H{"product"}   = 1; }
if(! opt_Get("--qnogene",    \%opt_HH)) { $qual_H{"gene"}      = 1; }
if(! opt_Get("--qnoexc",     \%opt_HH)) { $qual_H{"exception"} = 1; }
if(opt_IsUsed("--qadd", \%opt_HH)) { 
  my @qadd_A = split(",", opt_Get("--qadd", \%opt_HH));
  foreach my $q (@qadd_A) { $qual_H{$q} = 1; }
}

# remove all array elements with feature types not in %ftype_H, unless --fall used
if(! opt_Get("--fall", \%opt_HH)) { 
  my @ftr_idx_to_remove_A = ();
  for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$mdl_name}}); $ftr_idx++) { 
    my $ftype = $ftr_info_HAH{$mdl_name}[$ftr_idx]{"type"};
    if((! defined $ftype) || (! exists $ftype_H{$ftype})) { 
      splice(@{$ftr_info_HAH{$mdl_name}}, $ftr_idx, 1);
      $ftr_idx--; # this is about to be incremented
    }
  }
}

# remove any qualifier key/value pairs with keys not in %qual_H, unless --qall used
if(! opt_Get("--qall", \%opt_HH)) { 
  for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$mdl_name}}); $ftr_idx++) { 
    foreach $key (sort keys %{$ftr_info_HAH{$mdl_name}[$ftr_idx]}) { 
      if((! exists $qual_H{$key}) && (exists $ftr_info_HAH{$mdl_name}[$ftr_idx]{$key})) { 
        delete $ftr_info_HAH{$mdl_name}[$ftr_idx]{$key};
      }
    }
  }
}

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#####################################################################
# Parse the input stockholm file (if --stk) or create it (if ! --stk)
#####################################################################
my $stk_file = $out_root . ".stk";
my $fa_file  = $out_root . ".fa";
my $stk_has_ss = undef;
my $in_stk_file = opt_Get("--stk", \%opt_HH);
if(defined $in_stk_file) { 
  $start_secs = outputProgressPrior("Validating input Stockholm file", $progress_w, $log_FH, *STDOUT);

  $stk_has_ss = stockholm_validate_single_sequence_input($in_stk_file, $seq_info_HH{$mdl_name}{"seq"}, \%opt_HH, $FH_HR);

  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  
  $start_secs = outputProgressPrior("Reformatting Stockholm file to FASTA file", $progress_w, $log_FH, *STDOUT);

  runCommand("cp $in_stk_file $stk_file", opt_Get("-v", \%opt_HH), 0, $FH_HR);
  reformat_stockholm_file_to_unaligned_fasta_file($execs_H{"esl-reformat"}, $stk_file, $fa_file, \%opt_HH, $FH_HR);

  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}
else { 
  # --stk not used, we create it by first making a fasta file of the model
  # model sequence read from the gb file, then converting that fasta file 
  # to a stockholm file
  $start_secs = outputProgressPrior("Creating FASTA sequence file", $progress_w, $log_FH, *STDOUT);

  openAndAddFileToOutputInfo(\%ofile_info_HH, "fasta", $fa_file, 1, "fasta sequence file for $mdl_name");
  print_sequence_to_fasta_file($ofile_info_HH{"FH"}{"fasta"}, 
                               $seq_info_HH{$mdl_name}{"ver"}, 
                               $seq_info_HH{$mdl_name}{"def"}, 
                               $seq_info_HH{$mdl_name}{"seq"}, $FH_HR);
  close $ofile_info_HH{"FH"}{"fasta"};

  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  $start_secs = outputProgressPrior("Reformatting FASTA file to Stockholm file", $progress_w, $log_FH, *STDOUT);

  reformat_fasta_file_to_stockholm_file($execs_H{"esl-reformat"}, $fa_file, $stk_file, \%opt_HH, $FH_HR);
  addClosedFileToOutputInfo(\%ofile_info_HH, "stk", $stk_file, 1, "Stockholm alignment file for $mdl_name");

  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

######################################################################
# Finish populating @{$ftr_info_HAH{$mdl_name} and create @sgm_info_AH
######################################################################
$start_secs = outputProgressPrior("Finalizing feature information", $progress_w, $log_FH, *STDOUT);

featureInfoImputeCoords(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
featureInfoImputeLength(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
featureInfoImputeSourceIdx(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
featureInfoImputeParentIdx(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);

my @sgm_info_AH = (); # segment info, inferred from feature info
segmentInfoPopulate(\@sgm_info_AH, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###################################
# Translate the CDS, if we have any
###################################
my $ncds = featureInfoCountType(\@{$ftr_info_HAH{$mdl_name}}, "CDS");
my $cds_fa_file = undef;
my $protein_fa_file = undef;
if($ncds > 0) { 
  $start_secs = outputProgressPrior("Translating CDS and building BLAST DB", $progress_w, $log_FH, *STDOUT);

  $cds_fa_file  = $out_root . ".cds.fa";
  openAndAddFileToOutputInfo(\%ofile_info_HH, "cdsfasta", $cds_fa_file, 1, "fasta sequence file for CDS from $mdl_name");
  cdsFetchStockholmToFasta($ofile_info_HH{"FH"}{"cdsfasta"}, $stk_file, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  close $ofile_info_HH{"FH"}{"cdsfasta"};
  
  $protein_fa_file = $out_root . ".protein.fa";
  openAndAddFileToOutputInfo(\%ofile_info_HH, "proteinfasta", $protein_fa_file, 1, "fasta sequence file for translated CDS from $mdl_name");
  cdsTranslateToFastaFile($ofile_info_HH{"FH"}{"proteinfasta"}, $execs_H{"esl-translate"}, $cds_fa_file, 
                              $out_root, \@{$ftr_info_HAH{$mdl_name}}, \%opt_HH, $FH_HR);
  close $ofile_info_HH{"FH"}{"proteinfasta"};

  create_blast_protein_db($execs_H{"makeblastdb"}, $protein_fa_file, \%opt_HH, $FH_HR);

  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

##############
# Build the CM
##############
my $cm_file = undef;
if(! opt_Get("--skipbuild", \%opt_HH)) { 
  my $cmbuild_str = undef;
  my $clen_times_cmn = $seq_info_HH{$mdl_name}{"len"} * 200;
  if(opt_IsUsed("--cmn", \%opt_HH)) { 
    $clen_times_cmn *= (opt_Get("--cmn", \%opt_HH) / 200);
  }
  if   ($clen_times_cmn > 4000000) { $cmbuild_str = "(may take more than an hour)"; }
  elsif($clen_times_cmn > 3000000) { $cmbuild_str = "(should take roughly an hour)"; }
  elsif($clen_times_cmn > 2000000) { $cmbuild_str = "(should take roughly 20-40 minutes)"; }
  elsif($clen_times_cmn > 1000000) { $cmbuild_str = "(should take roughly 10-30 minutes)"; }
  elsif($clen_times_cmn >  500000) { $cmbuild_str = "(should take roughly 5-10 minutes)"; }
  else                             { $cmbuild_str = "(shouldn't take more than a few minutes)"; }

  $start_secs = outputProgressPrior("Building model $cmbuild_str", $progress_w, $log_FH, *STDOUT);

  my $cmbuild_opts = "-n $mdl_name --verbose ";
  if((! defined $stk_has_ss) || (! $stk_has_ss)) { $cmbuild_opts .= " --noss"; }
  if(opt_IsUsed("--cmn",    \%opt_HH)) { $cmbuild_opts .= " --EgfN " . opt_Get("--cmn", \%opt_HH); }
  if(opt_IsUsed("--cmp7ml", \%opt_HH)) { $cmbuild_opts .= " --p7ml"; }
  if(opt_IsUsed("--cmere",  \%opt_HH)) { $cmbuild_opts .= " --ere "  . opt_Get("--cmere", \%opt_HH); }
  if(opt_IsUsed("--cmeset", \%opt_HH)) { $cmbuild_opts .= " --eset " . opt_Get("--cmeset", \%opt_HH); }

  my $cmbuild_file = $out_root . ".cmbuild";
  $cm_file         = $out_root . ".cm";
  my $cmbuild_cmd  = $execs_H{"cmbuild"} . " " . $cmbuild_opts . " $cm_file $stk_file > $cmbuild_file";
  runCommand($cmbuild_cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});
  outputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);

  addClosedFileToOutputInfo(\%ofile_info_HH, "cm",      $cm_file, 1, "CM file");
  addClosedFileToOutputInfo(\%ofile_info_HH, "cmbuild", $cmbuild_file, 1, "cmbuild output file");
}

########################
# Output model info file
########################
$start_secs = outputProgressPrior("Creating model info file", $progress_w, $log_FH, *STDOUT);

# create @mdl_info_AH, and add info for our lone model
# modelInfoFileWrite() can output data for multiple models at once, 
# but we use it here only for a single model.
my @mdl_info_AH = (); 
%{$mdl_info_AH[0]} = ();
$mdl_info_AH[0]{"name"}   = $mdl_name;
$mdl_info_AH[0]{"length"} = $seq_info_HH{$mdl_name}{"len"};
if(defined $cm_file) { 
  $mdl_info_AH[0]{"cmfile"} = removeDirPath($cm_file);
}
if($ncds > 0) { 
  $mdl_info_AH[0]{"blastdb"} = removeDirPath($protein_fa_file);
  if((opt_IsUsed("--ttbl", \%opt_HH)) && (opt_Get("--ttbl", \%opt_HH) != 1))  { 
    $mdl_info_AH[0]{"transl_table"} = opt_Get("--ttbl", \%opt_HH);
  }
}
if(opt_IsUsed("-g", \%opt_HH)) { 
  $mdl_info_AH[0]{"group"} = opt_Get("-g", \%opt_HH); 
}
my $minfo_file  = $out_root . ".minfo";
modelInfoFileWrite($minfo_file, \@mdl_info_AH, \%ftr_info_HAH, $FH_HR);
addClosedFileToOutputInfo(\%ofile_info_HH, "minfo", $minfo_file, 1, "DNAORG 'model info' format file for $mdl_name");

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

##########
# Conclude
##########
# output optional output files
if(opt_Get("--sgminfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "sgminfo", $out_root . ".sgminfo", 1, "Model information (created due to --sgminfo)");
  dumpArrayOfHashes("Feature information (ftr_info_AH) for $mdl_name", \@{$ftr_info_HAH{$mdl_name}}, $ofile_info_HH{"FH"}{"ftrinfo"});
}
if(exists $ofile_info_HH{"FH"}{"sgminfo"}) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "ftrinfo", $out_root . ".ftrinfo", 1, "Feature information (created due to --ftrinfo)");
  dumpArrayOfHashes("Segment information (sgm_info_AH) for $mdl_name", \@sgm_info_AH, $ofile_info_HH{"FH"}{"sgminfo"});
}

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir, \%ofile_info_HH);
exit 0;

###############
# SUBROUTINES #
###############

#################################################################
# Subroutine: create_blast_protein_db
# Incept:     EPN, Wed Oct  3 16:31:38 2018
# 
# Purpose:    Create a protein blast database from a fasta file.
#
# Arguments:
#   $makeblastdb:    path to 'makeblastdb' executable
#   $fa_file:        FASTA file of protein sequences to make blast db from
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#                    
# Returns:    void
#
#################################################################
sub create_blast_protein_db {
  my $sub_name = "create_blast_protein_db";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($makeblastdb, $fa_file, $opt_HHR, $FH_HR) = @_;

  runCommand($makeblastdb . " -in $fa_file -dbtype prot > /dev/null", opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine: print_sequence_to_fasta_file()
# Incept:     EPN, Thu Mar 14 06:06:59 2019
#
# Synopsis: Print a sequence to a fasta file.
#
# Arguments:
#  $out_FH:    output file handle
#  $name:      sequence name
#  $def:       sequence definition, can be undef
#  $seq:       sequence string
#  $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $name or $seq is undef
#################################################################
sub print_sequence_to_fasta_file {
  my $sub_name = "print_sequence_to_fasta_file";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $name, $def, $seq, $FH_HR) = @_;

  if(! defined $name) { DNAORG_FAIL("ERROR in $sub_name, name is undefined", 1, $FH_HR); }
  if(! defined $seq)  { DNAORG_FAIL("ERROR in $sub_name, name is undefined", 1, $FH_HR); }

  # capitalize and DNAize $seq
  sqstringCapitalize(\$seq);
  sqstringDnaize(\$seq);
  printf $out_FH (">%s%s\n%s", 
                  $name, 
                  (defined $def) ? " " . $def : "",
                  sqstringAddNewlines($seq, 60));
  
  return;
}

#################################################################
# Subroutine: reformat_fasta_file_to_stockholm_file()
# Incept:     EPN, Sat Mar  9 10:24:14 2019
#
# Synopsis: Use esl-reformat to convert a fasta file to a stockholm file
#
# Arguments:
#  $esl_reformat: esl-reformat executable file
#  $fa_file:      fasta file
#  $stk_file:     stockholm file to create
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's a problem fetching the sequence file
#################################################################
sub reformat_fasta_file_to_stockholm_file { 
  my $sub_name = "reformat_fasta_file_to_stockholm_file";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_reformat, $fa_file, $stk_file, $opt_HHR, $FH_HR) = @_;

  my $cmd = $esl_reformat . " --informat afa stockholm $fa_file > $stk_file";
  runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine: reformat_stockholm_file_to_unaligned_fasta_file()
# Incept:     EPN, Fri Mar 15 12:56:04 2019
#
# Synopsis: Use esl-reformat to convert a stockholm file to unaligned fasta
#
# Arguments:
#  $esl_reformat: esl-reformat executable file
#  $stk_file:     stockholm file
#  $fa_file:      fasta file to create
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's a problem fetching the sequence file
#################################################################
sub reformat_stockholm_file_to_unaligned_fasta_file { 
  my $sub_name = "reformat_stockholm_file_to_unaligned_fasta_file";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_reformat, $stk_file, $fa_file, $opt_HHR, $FH_HR) = @_;

  my $cmd = $esl_reformat . " --informat stockholm fasta $stk_file > $fa_file";
  runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine: stockholm_validate_single_sequence_input()
# Incept:     EPN, Fri Mar 15 13:19:32 2019
#
# Synopsis: Validate an input Stockholm file is in the correct 
#           format, has exactly 1 sequence and no gaps.
#
# Arguments:
#  $in_stk_file:  input stockholm file to validate
#  $exp_sqstring: sequence we expect to be in the stockholm alignment
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    '1' if Stockholm file has SS_cons annotation, else '0'
#
# Dies:       if there's a problem parsing the file or 
#             a requirement is not met
#################################################################
sub stockholm_validate_single_sequence_input {
  my $sub_name = "stockholm_validate_single_sequence_input";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($in_stk_file, $exp_sqstring, $opt_HHR, $FH_HR) = @_;

  if(! -e $in_stk_file) { DNAORG_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file does not exist", 1, $FH_HR); }
  if(! -s $in_stk_file) { DNAORG_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file exists but is empty", 1, $FH_HR); }
  if(  -d $in_stk_file) { DNAORG_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file is actually a directory", 1, $FH_HR); }
  my $msa = Bio::Easel::MSA->new({ fileLocation => $in_stk_file, isDna => 1});
  my $nseq = $msa->nseq;
  if($nseq == 1) { 
    # single sequence, make sure there are no gaps
    if($msa->any_allgap_columns) { 
      DNAORG_FAIL("ERROR, read 1 sequence in --stk file $in_stk_file, but it has gaps, this is not allowed for single sequence 'alignments' (remove gaps with 'esl-reformat --mingap')", 1, $FH_HR);
    }
    # validate it matches $exp_sqstring
    my $fetched_sqstring = $msa->get_sqstring_unaligned(0);
    sqstringCapitalize(\$fetched_sqstring);
    sqstringCapitalize(\$exp_sqstring);
    sqstringDnaize(\$fetched_sqstring);
    sqstringDnaize(\$exp_sqstring);
    if($fetched_sqstring ne $exp_sqstring) { 
      my $summary_sqstring_diff_str = sqstringDiffSummary($fetched_sqstring, $exp_sqstring);
      DNAORG_FAIL("ERROR, read 1 sequence in --stk file $in_stk_file, but it does not match sequence read from GenBank file $gb_file:\n$summary_sqstring_diff_str", 1, $FH_HR); 
    }
  }
  else { # nseq != 1
    DNAORG_FAIL("ERROR, did not read exactly 1 sequence in --stk file $in_stk_file.\nTo use DNAORG with models built from alignments of multiple sequences,\nyou will have to build the CM with cmbuild and create the model info file manually.\n", 1, $FH_HR);
  }

  return $msa->has_ss;
}

