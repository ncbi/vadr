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
# Preliminaries: 
#   - process options
#   - create the output directory
#   - output program banner and open output files
#   - parse the optional input files, if necessary
#   - make sure the required executables are executable
#
# 
#######################################################################################

# first, determine the paths to all modules, scripts and executables that we'll need

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

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
$opt_group_desc_H{"1"} = "basic options";
#     option            type       default               group   requires incompat    preamble-output                                                help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                                         "display this help",                                   \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",                                 "force; if dir <output directory> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-m",           "boolean", 0,                        1,    undef, undef,      "first cmdline arg is an input .minfo file, not an accession", "first cmdline arg is an input .minfo file, not an accession", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                                  "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",                          "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"1"} = "options related to alternative modes:";
#       option          type       default                group  requires incompat        preamble-output                                            help-output    
opt_Add("--stk",        "string",  undef,                    2,    undef,  undef,         "read stockholm alignment with >= 1 sequences from <s>",   "read stockholm alignment with >= 1 sequences from <s>", \%opt_HH, \@opt_order_A);
opt_Add("--gb",         "string",  undef,                    2,    undef, "-m",           "read genbank file from <s>, don't fetch it",              "read genbank file from <s>, don't fetch it", \%opt_HH, \@opt_order_A);
opt_Add("--onlyurl",    "boolean", 0,                        2,    undef,"--stk,--gb,-m", "output genbank file url for accession and exit",          "output genbank file url for accession and exit", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"3"} = "optional output files";
#       option       type       default                group  requires incompat  preamble-output                          help-output    
opt_Add("--ftrinfo",    "boolean", 0,                        3,    undef, undef, "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);
opt_Add("--sgminfo",    "boolean", 0,                        3,    undef, undef, "output internal segment information",   "create file with internal segment information",   \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"4"} = "options for controlling cmbuild step";
#     option          type       default         group   requires    incompat   preamble-output                                             help-output    
opt_Add("--cmn",      "integer", 0,                  4,   undef,     undef,     "set number of seqs for glocal fwd HMM calibration to <n>", "set number of seqs for glocal fwd HMM calibration to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--cmp7ml",   "boolean", 0,                  4,   undef,     undef,     "set CM's filter p7 HMM as the ML p7 HMM",                  "set CM's filter p7 HMM as the ML p7 HMM",                  \%opt_HH, \@opt_order_A);
opt_Add("--cmere",    "real",    0,                  4,   undef,     undef,     "set CM relative entropy target to <x>",                    "set CM relative entropy target to <x>",                    \%opt_HH, \@opt_order_A);
opt_Add("--cmeset",   "real",    0,                  4,   undef,     "--cmere", "set CM eff seq # for CM to <x>",                           "set CM eff seq # for CM to <x>",                           \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_build.pl [-options] <accession (or input .minfo file if -m)> <path to output directory to create>\n";
my $synopsis = "dnaorg_build.pl :: build homology model for feature annotation";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'm'            => \$GetOptions_H{"-m"},
                'stk=s'        => \$GetOptions_H{"--stk"},
                'gb=s'         => \$GetOptions_H{"--gb"},
                'keep'         => \$GetOptions_H{"--keep"},
# options for alternative modes:
                'onlyurl'      => \$GetOptions_H{"--onlyurl"},
# optional output files
                'sgminfo'      => \$GetOptions_H{"--sgminfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"},
# options for controlling cmbuild step
                'cmn=s'        => \$GetOptions_H{"--cmn"},
                'cmp7ml'       => \$GetOptions_H{"--cmp7ml"},
                'cmere=s'      => \$GetOptions_H{"--cmere"},
                'cmeset=s'     => \$GetOptions_H{"--cmeset"});

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
my ($acc_or_minfo, $dir) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

############################################
# if --onlyurl used, output the url and exit
############################################
if(opt_Get("--onlyurl", \%opt_HH)) { 
  print eutilsFetchUrl($acc_or_minfo, "gb") . "\n";
  exit 0;
}

################################################
# validate .minfo file exists and parse it if -m
################################################
my $in_minfo_file = opt_Get("-m", \%opt_HH) ? $acc_or_minfo : undef;
my $model_name = undef;
my @mdl_info_AH  = (); # array of hashes with model info
my %ftr_info_HAH = (); # array of hashes with feature info 
if(defined $in_minfo_file) { 
  if(! -e $in_minfo_file) { 
    die "ERROR, -m enabled, model info file $in_minfo_file does not exist";
  }
  if(! -s $in_minfo_file) { 
    die "ERROR, -m enabled, model info file $in_minfo_file exists but is empty";
  }
  if(-d $in_minfo_file) { 
    die "ERROR, -m enabled, model info file $in_minfo_file is actually a directory";
  }
  # set required keys for models and features
  my @reqd_mdl_keys_A = ("name");
  my @reqd_ftr_keys_A = ("type", "coords");
  modelInfoFileParse($in_minfo_file, \@mdl_info_AH, \%ftr_info_HAH, \@reqd_mdl_keys_A, \@reqd_ftr_keys_A, undef);
  # ensure we read info only a single model from the input file
  if(scalar(@mdl_info_AH) == 0) { 
    die "ERROR did not read info for any models in $in_minfo_file";
  }
  if(scalar(@mdl_info_AH) > 1) { 
    die "ERROR read info for multiple models in $in_minfo_file, it must have info for only 1 model";
  }
  if(! defined $mdl_info_AH[0]{"name"}) { 
    die "ERROR did not read name of model in $in_minfo_file";
  }
  $model_name = $mdl_info_AH[0]{"name"};
}
else { 
  $model_name = $acc_or_minfo;
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

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("accession (or .minfo file if -m)", "output directory");
my @arg_A      = ($acc_or_minfo, $dir);
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
$execs_H{"cmbuild"}       = $inf_exec_dir . "/cmbuild";
$execs_H{"esl-reformat"}  = $esl_exec_dir . "/esl-reformat";
$execs_H{"esl-translate"} = $esl_exec_dir . "/esl-translate";
$execs_H{"makeblastdb"}   = $blast_exec_dir . "/makeblastdb";
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################################################
# Fetch the genbank file (if neither -m nor --gb used)
###########################################################
my $progress_w = 50; # the width of the left hand column in our progress output, hard-coded
my $start_secs;
my $gb_file = undef;
if(opt_IsUsed("--gb", \%opt_HH)) { 
  $gb_file = opt_Get("--gb", \%opt_HH);
}
elsif(opt_IsUsed("-m", \%opt_HH)) { 
  $gb_file = undef;
}
else { 
  # neither --gb nor -m used, create gb file by fetching using eutils
  $gb_file = $out_root . ".gb";
  $start_secs = outputProgressPrior("Fetching GenBank file", $progress_w, $log_FH, *STDOUT);

  eutilsFetchToFile($gb_file, $model_name, "gb", 5, $ofile_info_HH{"FH"});  # number of attempts to fetch to make before dying
  addClosedFileToOutputInfo(\%ofile_info_HH, "gb", $gb_file, 1, "GenBank format file for $model_name");

  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

#########################################
# Parse the genbank file (if -m not used)
#########################################
my %seq_info_HH  = ();
if(defined $gb_file) { 
  $start_secs = outputProgressPrior("Parsing GenBank file", $progress_w, $log_FH, *STDOUT);
  genbankParse($gb_file, \%seq_info_HH, \%ftr_info_HAH, $FH_HR);
  if(! exists $ftr_info_HAH{$model_name}) { 
    DNAORG_FAIL("ERROR parsing GenBank file $gb_file, did not read info for reference accession $model_name\n", 1, $FH_HR);
  }

  my $ftr_idx;
  my $key;
  my %ftype_H = ();
  $ftype_H{"CDS"}         = 1;
  $ftype_H{"gene"}        = 1;
  $ftype_H{"mat_peptide"} = 1;
  
  my %qual_H = ();
  $qual_H{"type"}         = 1;
  $qual_H{"location"}     = 1;
  $qual_H{"product"}      = 1;
  $qual_H{"gene"}         = 1;
  $qual_H{"exception"}    = 1;

  # first, remove all array elements with types not in %ftype_H
  my @ftr_idx_to_remove_A = ();
  for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$model_name}}); $ftr_idx++) { 
    my $ftype = $ftr_info_HAH{$model_name}[$ftr_idx]{"type"};
    if((! defined $ftype) || (! exists $ftype_H{$ftype})) { 
      splice(@{$ftr_info_HAH{$model_name}}, $ftr_idx, 1);
      $ftr_idx--; # this is about to be incremented
    }
  }

  # now go back through and remove any key/value pairs not in %qual_H
  for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$model_name}}); $ftr_idx++) { 
    foreach $key (sort keys %{$ftr_info_HAH{$model_name}[$ftr_idx]}) { 
      if((! exists $qual_H{$key}) && (exists $ftr_info_HAH{$model_name}[$ftr_idx]{$key})) { 
        delete $ftr_info_HAH{$model_name}[$ftr_idx]{$key};
      }
    }
  }
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

###################################################################################################
# Finish populating @{$ftr_info_HAH{$model_name} and create @sgm_info_AH and output model info file
###################################################################################################

$start_secs = outputProgressPrior("Creating model info file", $progress_w, $log_FH, *STDOUT);

# remove any features that are not of the type that we care about
# and from the features we do care about, remove keys we don't care about

# TODO, optionally prune to keep certain feature types, or all (?)
# TODO, optionally prune to keep certain qualifiers, or all (?)

# only impute coords if we didn't read them from an minfo file
if(! defined $in_minfo_file) { 
  featureInfoImputeCoords(\@{$ftr_info_HAH{$model_name}}, $FH_HR);
}
featureInfoImputeLength(\@{$ftr_info_HAH{$model_name}}, $FH_HR);
featureInfoImputeSourceIdx(\@{$ftr_info_HAH{$model_name}}, $FH_HR);
featureInfoImputeParentIdx(\@{$ftr_info_HAH{$model_name}}, $FH_HR);

my @sgm_info_AH = (); # segment info, inferred from feature info
segmentInfoPopulate(\@sgm_info_AH, \@{$ftr_info_HAH{$model_name}}, $FH_HR);

my $minfo_file  = $out_root . ".minfo";
my $cm_file     = $out_root . ".cm";
if(scalar(@mdl_info_AH) == 0) { # initialize
  %{$mdl_info_AH[0]} = ();
  $mdl_info_AH[0]{"name"} = $model_name;
}
$mdl_info_AH[0]{"cmfile"} = $cm_file;
modelInfoFileWrite($minfo_file, \@mdl_info_AH, \%ftr_info_HAH, $FH_HR);
addClosedFileToOutputInfo(\%ofile_info_HH, "minfo", $minfo_file, 1, "GenBank format file for $model_name");

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

exit 0;
######################################
# Output the fasta and stockholm files
######################################
# verify our model and feature info hashes are complete, 
# if validateFeatureInfoHashIsComplete() fails then the program will exit with an error message

$start_secs = outputProgressPrior("Creating FASTA and STOCHOLM sequence files", $progress_w, $log_FH, *STDOUT);

my $fa_file  = $out_root . ".fa";
openAndAddFileToOutputInfo(\%ofile_info_HH, "fasta", $fa_file, 1, "fasta sequence file for $model_name");
print_sequence_to_fasta_file($ofile_info_HH{"FH"}{"fasta"}, 
                             $seq_info_HH{$model_name}{"ver"}, 
                             $seq_info_HH{$model_name}{"def"}, 
                             $seq_info_HH{$model_name}{"seq"}, $FH_HR);
close $ofile_info_HH{"FH"}{"fasta"};

my $stk_file = $out_root . ".stk";
reformat_fasta_file_to_stockholm_file($execs_H{"esl-reformat"}, $fa_file, $stk_file, \%opt_HH, $FH_HR);
addClosedFileToOutputInfo(\%ofile_info_HH, "stk", $stk_file, 1, "Stockholm alignment file for $model_name");

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###################
# Translate the CDS
###################
my $ncds = featureInfoCountType(\@{$ftr_info_HAH{$model_name}}, "CDS");
if($ncds > 0) { 
  $start_secs = outputProgressPrior("Translating CDS and building BLAST DB", $progress_w, $log_FH, *STDOUT);

  my $cds_fa_file  = $out_root . ".cds.fa";
  openAndAddFileToOutputInfo(\%ofile_info_HH, "cdsfasta", $cds_fa_file, 1, "fasta sequence file for CDS from $model_name");
  fetch_cds_to_fasta_file($ofile_info_HH{"FH"}{"cdsfasta"}, $fa_file, $seq_info_HH{$model_name}{"ver"}, \@{$ftr_info_HAH{$model_name}}, $FH_HR);
  close $ofile_info_HH{"FH"}{"cdsfasta"};
  
  my $protein_fa_file  = $out_root . ".protein.fa";
  openAndAddFileToOutputInfo(\%ofile_info_HH, "proteinfasta", $protein_fa_file, 1, "fasta sequence file for translated CDS from $model_name");
  translate_cds_to_fasta_file($ofile_info_HH{"FH"}{"proteinfasta"}, $execs_H{"esl-translate"}, $cds_fa_file, 
                              $out_root, $seq_info_HH{$model_name}{"ver"}, \@{$ftr_info_HAH{$model_name}}, \%opt_HH, $FH_HR);
  close $ofile_info_HH{"FH"}{"proteinfasta"};

  create_blast_protein_db($execs_H{"makeblastdb"}, $protein_fa_file, \%opt_HH, $FH_HR);
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

##############
# Build the CM
##############
$start_secs = outputProgressPrior("Building model (this could take a while)", $progress_w, $log_FH, *STDOUT);

my $cmbuild_opts = "-n $model_name --noss --verbose -F";
if(opt_IsUsed("--cmn",    \%opt_HH)) { $cmbuild_opts .= " --EgfN " . opt_Get("--cmn", \%opt_HH); }
if(opt_IsUsed("--cmp7ml", \%opt_HH)) { $cmbuild_opts .= " --p7ml"; }
if(opt_IsUsed("--cmere",  \%opt_HH)) { $cmbuild_opts .= " --ere "  . opt_Get("--cmere", \%opt_HH); }
if(opt_IsUsed("--cmeset", \%opt_HH)) { $cmbuild_opts .= " --eset " . opt_Get("--cmeset", \%opt_HH); }

my $cmbuild_file = $out_root . ".cmbuild";
my $cmbuild_cmd  = $execs_H{"cmbuild"} . " " . $cmbuild_opts . " $cm_file $stk_file > $cmbuild_file";
#runCommand($cmbuild_cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);

addClosedFileToOutputInfo(\%ofile_info_HH, "cm",      $cm_file, 1, "CM file");
addClosedFileToOutputInfo(\%ofile_info_HH, "cmbuild", $cmbuild_file, 1, "cmbuild output file");

##########
# Conclude
##########
# output optional output files
if(opt_Get("--sgminfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "sgminfo", $out_root . ".sgminfo", 1, "Model information (created due to --sgminfo)");
  dumpArrayOfHashes("Feature information (ftr_info_AH) for $model_name", \@{$ftr_info_HAH{$model_name}}, $ofile_info_HH{"FH"}{"ftrinfo"});
}
if(exists $ofile_info_HH{"FH"}{"sgminfo"}) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "ftrinfo", $out_root . ".ftrinfo", 1, "Feature information (created due to --ftrinfo)");
  dumpArrayOfHashes("Segment information (sgm_info_AH) for $model_name", \@sgm_info_AH, $ofile_info_HH{"FH"}{"sgminfo"});
}

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir, \%ofile_info_HH);
exit 0;

#################################################################
# Subroutine: fetch_cds_to_fasta_file()
# Incept:     EPN, Thu Mar 14 12:30:33 2019
# 
# Purpose:    Fetch the protein translations of the CDS for the genome
#             and create multiple N+1 FASTA files, one with each
#             single sequence (N) and one with all sequences.
#             Fill @{$fa_file_AR} with the sequence file names.
#
# Arguments:
#   $out_FH:         output file handle
#   $fa_file:        fasta file with full length sequence of sequence $seqname
#   $seq_name:       sequence to fetch CDS from
#   $ftr_info_AHR:   REF to the feature info, pre-filled
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#                    
# Returns: void
#
# Dies:    if we have trouble fetching a sequence
#
#################################################################
sub fetch_cds_to_fasta_file { 
  my $sub_name = "fetch_cds_to_fasta_file";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $fa_file, $seq_name, $ftr_info_AHR, $FH_HR) = @_;

  my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $fa_file }); # the sequence file object
  my $nftr = scalar(@{$ftr_info_AHR});

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") { 
      my $cds_sqstring = "";
      my @start_A = ();
      my @stop_A = ();
      my @strand_A = ();
      featureStartStopStrandArrays($ftr_info_AHR->[$ftr_idx]{"coords"}, \@start_A, \@stop_A, \@strand_A, $FH_HR);
      foreach(my $sgm_idx = 0; $sgm_idx < scalar(@start_A); $sgm_idx++) { 
        $cds_sqstring .= $sqfile->fetch_subseq_to_sqstring($seq_name, $start_A[$sgm_idx], $stop_A[$sgm_idx], -1, ($strand_A[$sgm_idx] eq "-") ? 1 : 0);
      }
      print $out_FH(">" . $seq_name . "/" . $ftr_info_AHR->[$ftr_idx]{"coords"} . "\n" . sqstringAddNewlines($cds_sqstring, 60));
    }
  }
  return;
}

#################################################################
# Subroutine: translate_cds_to_fasta_file()
# Incept:     EPN, Thu Mar 14 12:30:28 2019
# 
# Purpose:    Use esl-translate to translate a fasta file with
#             CDS sequences pertaining to the CDS features in 
#             @{$ftr_info_AHR} into fasta protein files.
#
# Arguments:
#   $out_FH:         output file handle
#   $esl_translate:  path to esl-translate executable
#   $cds_fa_file:    fasta file with CDS sequences
#   $out_root:       string that is the 'root' for naming output files
#   $seq_name:       sequence to fetch CDS from
#   $ftr_info_AHR:   REF to the feature info, pre-filled
#   $opt_HHR:        command line options
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#                    
# Returns: void
#
# Dies:    if we have trouble fetching a sequence
#
#################################################################
sub translate_cds_to_fasta_file { 
  my $sub_name = "fetch_and_translate_proteins_into_fasta_files";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $esl_translate, $cds_fa_file, $out_root, $seq_name, $ftr_info_AHR, $opt_HHR, $FH_HR) = @_;

  my $tmp1_translate_fa_file = $out_root . ".cds.esl-translate.1.fa";
  my $tmp2_translate_fa_file = $out_root . ".cds.esl-translate.2.fa";
  my $translate_cmd = "$esl_translate -M -l 3 --watson $cds_fa_file > $tmp1_translate_fa_file";
  runCommand($translate_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # go through output fasta file and rewrite names, so we can fetch 
  open(IN,       $tmp1_translate_fa_file) || fileOpenFailure($tmp1_translate_fa_file, $sub_name, $!, "reading", $FH_HR);
  open(OUT, ">", $tmp2_translate_fa_file) || fileOpenFailure($tmp2_translate_fa_file, $sub_name, $!, "writing", $FH_HR);
  while(my $line = <IN>) { 
    if($line =~ m/^\>/) { 
      #>orf58 source=NC_039477.1/5..5104:+ coords=1..5097 length=1699 frame=1  
      chomp $line;
      if($line =~ /^\>orf\d+\s+(source\=\S+)\s+(coords\=\S+)\s+length\=\d+\s+frame\=\S+/) { 
        # rename as 'source=NC_039477.1/5..5104:+,coords=1..5097'
        print OUT (">" . $1 . "," . $2 . "\n");
      }
      else { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing esl-translate output file $tmp1_translate_fa_file, line:\n$line\n", 1, $FH_HR);
      }
    }
    else { 
      print OUT $line; 
    }
  }
  close(IN);
  close(OUT);

  # $tmp2_translate_fa_file now includes renamed translated sequences from esl-translate
  # fetch expected translated seqs and print to $out_FH
  my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $tmp2_translate_fa_file }); # the sequence file object
  
  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") { 
      my $fetch_name = "source=" . $seq_name . "/" . $ftr_info_AHR->[$ftr_idx]{"coords"} . ",coords=1.." . ($ftr_info_AHR->[$ftr_idx]{"length"} - 3); # subtract length of stop codon
      if(! $sqfile->check_seq_exists($fetch_name)) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, problem translating CDS feature (coords: %s)\n", $ftr_info_AHR->[$ftr_idx]{"coords"}), 1, $FH_HR);
      }
      print $out_FH $sqfile->fetch_seq_to_fasta_string($fetch_name, 60);
    }
  }

  # remove temporary files unless --keep
  if(! opt_Get("--keep", $opt_HHR)) { 
    removeFileUsingSystemRm($tmp1_translate_fa_file, $sub_name, $opt_HHR, $FH_HR);
    removeFileUsingSystemRm($tmp2_translate_fa_file, $sub_name, $opt_HHR, $FH_HR);
    removeFileUsingSystemRm($tmp2_translate_fa_file . ".ssi", $sub_name, $opt_HHR, $FH_HR);
  }

  return;
}

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
# Subroutine: modelInfoFileWrite()
# Incept:     EPN, Sat Mar  9 05:27:15 2019
#
# Synopsis: Output a model info file for model $mdlname based on 
#           feature information in @{$ftr_info_AHR}.
#
#           The following values must be set in @{$ftr_info_AHR}:
#             "type":   feature type, e.g. "mat_peptide", "CDS"
#             "coords": coordinates for this feature in the reference
#
# Arguments:
#  $out_file:      out file to create
#  $mdl_info_AHR:  REF to array of hashes with model info, pre-filled
#  $ftr_info_HAHR: REF to hash of array of hashes with information on the features, pre-filled
#  $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $ftr_info_HAR is not valid upon entering
#################################################################
sub modelInfoFileWrite { 
  my $sub_name = "modelInfoFileWrite";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_file, $mdl_info_AHR, $ftr_info_HAHR, $FH_HR) = @_;

  my $mdl_idx  = undef; # model index
  my $ftr_idx  = undef; # feature index
  my $nmdl     = undef; # number of models
  my $nftr     = undef; # number of features
  my $mdl_name = undef; # model name
  my $key      = undef; # a hash key
  my $value    = undef; # a hash value

  my @reqd_mdl_keys_A  = ("name", "cmfile");
  my @reqd_ftr_keys_A  = ("type", "coords");
  # validate all info first
  $nmdl = arrayOfHashesValidate($mdl_info_AHR, \@reqd_mdl_keys_A, "ERROR in $sub_name, mdl_info failed validation; missing required key(s)", $FH_HR);
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    arrayOfHashesValidate(\@{$ftr_info_HAHR->{$mdl_name}}, \@reqd_ftr_keys_A, "ERROR in $sub_name, ftr_info failed validation; missing required key(s)", $FH_HR);
  }

  # output 
  open(OUT, ">", $out_file) || fileOpenFailure($out_file, $sub_name, $!, "writing", $FH_HR);
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    print OUT ("MODEL $mdl_name");
    foreach $key (sort keys (%{$mdl_info_AHR->[$mdl_idx]})) { 
      $value = $mdl_info_AHR->[$mdl_idx]{$key};
      if($key =~ m/\:/) { 
        DNAORG_FAIL("ERROR in $sub_name, problem writing $out_file, illegal ':' character in model key $key for model $mdl_name", 1, $FH_HR);
      }
      if($value =~ m/\"/) { 
        DNAORG_FAIL("ERROR in $sub_name, problem writing $out_file, illegal '\"' character in model value $value for key $key for model $mdl_name", 1, $FH_HR);
      }
      if($key ne "name") { 
        print OUT (" $key:\"$value\"");
      }
    }
    print OUT ("\n");

    # define feature keys to ignore
    my %ftr_key_ignore_H = ();
    $ftr_key_ignore_H{"type"}        = 1; # this automatically gets added to @key_order_A, so it goes first
    $ftr_key_ignore_H{"coords"}      = 1; # this automatically gets added to @key_order_A, so it goes second
    $ftr_key_ignore_H{"length"}      = 1; # will be inferred from coords
    $ftr_key_ignore_H{"source_idx"}  = 1; # will be inferred from coords and type
    $ftr_key_ignore_H{"parent_idx"}  = 1; # will be inferred from coords and type
    $ftr_key_ignore_H{"3pa_ftr_idx"} = 1; # will be inferred from coords and type
    $ftr_key_ignore_H{"5p_sgm_idx"}  = 1; # will be inferred from coords, when sgm_info_HA is created
    $ftr_key_ignore_H{"3p_sgm_idx"}  = 1; # will be inferred from coords, when sgm_info_HA is created
    $ftr_key_ignore_H{"location"}    = 1; # *could* (but won't be) inferred from coords

    $nftr = scalar(@{$ftr_info_HAHR->{$mdl_name}});
    # determine order of keys for this feature
    my @ftr_key_order_A  = ("type", "coords");
    for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      foreach $key (sort keys %{$ftr_info_HAHR->{$mdl_name}[$ftr_idx]}) { 
        if(! exists $ftr_key_ignore_H{$key}) { 
          push(@ftr_key_order_A, $key);
          $ftr_key_ignore_H{$key} = 1; 
        }
      }
    }
    
    for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      print OUT ("FEATURE $mdl_name");
      foreach $key (@ftr_key_order_A) { 
        if(exists $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key}) { 
          $value = $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key};
          if($key =~ m/\:/) { 
            DNAORG_FAIL("ERROR in $sub_name, problem writing $out_file, illegal ':' character in feature key $key for model $mdl_name", 1, $FH_HR);
          }
          if($value =~ m/\"/) { 
            DNAORG_FAIL("ERROR in $sub_name, problem writing $out_file, illegal '\"' character in feature value $value for key $key for model $mdl_name", 1, $FH_HR);
          }
          print OUT (" $key:\"$value\"");
        }
      }
      print OUT ("\n");
    }
  }
  close(OUT);

  return;
}

#################################################################
# Subroutine: modelInfoFileParse()
# Incept:     EPN, Fri Mar 15 05:15:23 2019
#
# Synopsis: Parse a model info file for >= 1 models and collect 
#           feature information for each model $model in 
#           @{$ftr_info_HAHR->{$model}}.
#
#           The following keys must be defined for all features:
#             "type":   feature type, e.g. "mat_peptide", "CDS"
#             "coords": coordinates for this feature in the reference
#           We verify this at end of subroutine
# 
# Arguments:
#  $in_file:          input .minfo file to parse
#  $mdl_info_AHR:     REF to array of hashes of model information, filled here
#  $ftr_info_HAHR:    REF to hash of array of hashes with information 
#                     on the features per model, filled here
#  $reqd_mdl_keys_AR: REF to array of required keys read in MODEL lines
#  $reqd_ftr_keys_AR: REF to array of required keys read in FEATURE lines
#  $FH_HR:            REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if unable to parse $in_file
#             if a feature is defined without "type" or "coords" keys
#################################################################
sub modelInfoFileParse {
  my $sub_name = "modelInfoFileParse";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($in_file, $mdl_info_AHR, $ftr_info_HAHR, $reqd_mdl_keys_AR, $reqd_ftr_keys_AR, $FH_HR) = @_;
  
  my $format_str = "# DNAORG model info (.minfo) format specifications:\n";
  $format_str   .= "# Lines prefixed with '#' are ignored.\n";
  $format_str   .= "# All other lines must begin with either: 'MODEL' or 'FEATURE'\n";
  $format_str   .= "# followed by one or more whitespace characters and then the model\n";
  $format_str   .= "# name <modelname> which cannot include whitespace.\n";
  $format_str   .= "# On each line after <modelname>, both MODEL and FEATURE lines must\n";
  $format_str   .= "# contain 0 or more <key>:<value> pairs meeting the following criteria.\n";
  $format_str   .= "# <key> must not include any whitespace or ':' characters\n";
  $format_str   .= "# <value> must start *and* end with '\"' but include no other '\"'\n";
  $format_str   .= "# characters (but <value> may include whitespace characters).\n";
  $format_str   .= "# <key>:<value> pairs must be separated by one or more whitespace characters.\n";
  $format_str   .= "# <modelname> and the first <key>:<value> pair must be separated by one or\n";
  $format_str   .= "# more whitespace characters.\n";

  # example lines:
  #MODEL NC_039477 cmfile:"test/test.dnaorg_build.cm"
  #FEATURE NC_039477 type:"gene" coords:"5..5104:+" gene:"ORF1"
  #FEATURE NC_039477 type:"CDS" coords:"5..5104:+" gene:"ORF1" product:"nonstructural polyprotein"

  my $mdl_name   = undef; # name of current model
  my $ftr_idx    = undef; # index of current feature
  my $mdl_idx    = -1;    # index of current model
  my %mdl_read_H = ();    # keeps track of which model names we've seen MODEL lines for, to avoid duplicates
  open(IN, $in_file) || fileOpenFailure($in_file, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    if($line !~ /^#/) { 
      # not a comment line
      my $orig_line = $line;
      chomp $line;
      my $is_model_line = 0; # set to 1 if line we are parsing is a MODEL line, else it's a FEATURE line
      if($line =~ /^MODEL\s+(\S+)\s*/) { 
        $mdl_name = $1;
        if(exists $mdl_read_H{$mdl_name}) { 
          DNAORG_FAIL("ERROR in $sub_name, problem parsing $in_file: read multiple MODEL lines for $mdl_name, should only be 1; line:\n$orig_line\n", 1, $FH_HR);
        }
        $mdl_idx++;
        %{$mdl_info_AHR->[$mdl_idx]} = ();
        @{$ftr_info_HAHR->{$mdl_name}} = ();
        $mdl_info_AHR->[$mdl_idx]{"name"} = $mdl_name;
        $mdl_read_H{$mdl_name} = 1;

        $is_model_line = 1;
        $line =~ s/^MODEL\s+(\S+)\s*//; # remove MODEL and model value
      }
      elsif($line =~ /^FEATURE\s+(\S+)\s*/) { 
        $mdl_name = $1;
        if(! exists $mdl_read_H{$mdl_name}) { 
          DNAORG_FAIL("ERROR in $sub_name, problem parsing $in_file: read FEATURE line for model $mdl_name before a MODEL line for $mdl_name; line:\n$orig_line\n", 1, $FH_HR);
        }
        $ftr_idx = scalar(@{$ftr_info_HAHR->{$mdl_name}});
        # initialize ftr_info for this model/feature pair
        %{$ftr_info_HAHR->{$mdl_name}[$ftr_idx]} = (); 
        $line =~ s/^FEATURE\s+\S+\s*//; # remove FEATURE and model value
      }
      else { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $in_file, non-comment line does not start with 'MODEL <modelname>' or 'FEATURE <featurename>', line:\n$orig_line\n", 1, $FH_HR);
      }
      # if we get here we have either a MODEL or FEATURE line, parse the rest of it
      while($line ne "") { 
        if($line =~ /^([^\:\s]+)\:\"([^\"]+)\"\s*/) { 
          # key   must not include ':' or whitespace
          # value must begin and end with '"' but otherwise include on '"' characters
          my ($key, $value) = ($1, $2);
          if($is_model_line) { 
            if(exists $mdl_info_AHR->[$mdl_idx]{$key}) {
              DNAORG_FAIL("ERROR in $sub_name, problem parsing $in_file, read multiple values for key $key on MODEL line; line:\n$orig_line\n", 1, $FH_HR);
            }
            $mdl_info_AHR->[$mdl_idx]{$key} = $value;
          }
          else { # feature line
            if(exists $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key}) {
              DNAORG_FAIL("ERROR in $sub_name, problem parsing $in_file, read multiple values for key $key on MODEL line; line:\n$orig_line\n", 1, $FH_HR);
            }
            $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key} = $value;
            # printf("\tadded ftr_info_HAR->{$mdl_name}[$ftr_idx]{$key} as $value\n");
          }
          $line =~ s/^[^\:\s]+\:\"[^\"]+\"\s*//; # remove this key/value pair
        }
        else { 
          DNAORG_FAIL("ERROR in $sub_name, unable to parse $in_file, failed to parse key:value pairs in line:\n$orig_line\n$format_str\n", 1, $FH_HR);
        }
      } 
    }
  }
  close(IN);

  # verify we read what we need
  arrayOfHashesValidate($mdl_info_AHR, $reqd_mdl_keys_AR, "ERROR in $sub_name, problem parsing $in_file, required MODEL key missing", $FH_HR);
  my $nmdl = scalar(@{$mdl_info_AHR});
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    arrayOfHashesValidate($ftr_info_HAHR->{$mdl_name}, $reqd_ftr_keys_AR, "ERROR in $sub_name, problem parsing $in_file, required MODEL key missing for model " . $mdl_info_AHR->[$mdl_idx]{"name"}, $FH_HR);
  }

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

  # capitalize $seq
  $seq =~ tr/a-z/A-Z/;
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

