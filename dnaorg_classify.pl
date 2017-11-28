#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Mon Feb  1 15:07:43 2016 [dnaorg_build.pl split off from dnaorg_annotate_genomes.pl]
# LES, Mon Jul 25 09:25    2016 [dnaorg_refseq_assign.pl split off from dnaorg_build.pl]
# EPN, Tue Nov 28 10:44:54 2017 [dnaorg_classify.pl created, renamed from dnaorg_refseq_assign.pl]
#
# Some code in here is also adopted from dnaorg.pm

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

use Data::Dumper;
use File::Slurp;

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
opt_Add("--nseq",       "integer", 10,                       3,    undef,"--local",   "number of sequences for each nhmmscan farm job",              "set number of sequences for each nhmmscan farm job to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--maxnjobs",   "integer", 500,                      1,    undef,"--local",   "maximum allowed number of jobs for compute farm",             "set max number of jobs to submit to compute farm to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--wait",       "integer", 500,                      1,    undef,"--local",   "allow <n> minutes for nhmmscan jobs on farm",                 "allow <n> wall-clock minutes for nhmmscan jobs on farm to finish, including queueing time", \%opt_HH, \@opt_order_A);
opt_Add("--local",      "boolean", 0,                        3,    undef, undef,      "run nhmmscan locally instead of on farm",                     "run nhmmscan locally instead of on farm", \%opt_HH, \@opt_order_A);
opt_Add("--errcheck",   "boolean", 0,                        1,    undef,"--local",   "consider any farm stderr output as indicating a job failure", "consider any farm stderr output as indicating a job failure", \%opt_HH, \@opt_order_A);


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
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'keep'         => \$GetOptions_H{"--keep"},
                'dirout=s'     => \$GetOptions_H{"--dirout"},
                'onlybuild=s'  => \$GetOptions_H{"--onlybuild"},
                'dirbuild=s'   => \$GetOptions_H{"--dirbuild"},
                'inlist=s'     => \$GetOptions_H{"--inlist"},
                'infasta=s'    => \$GetOptions_H{"--infasta"},
                'nseq=s'       => \$GetOptions_H{"--nseq"}, 
                'maxnjobs=s'   => \$GetOptions_H{"--maxnjobs"}, 
                'wait=s'       => \$GetOptions_H{"--wait"},
                'local'        => \$GetOptions_H{"--local"}, 
                'errcheck'     => \$GetOptions_H{"--errcheck"},  
                );

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.19";
my $releasedate   = "Nov 2017";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date);
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
outputBanner($log_FH, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $cmd_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

my $do_keep = opt_Get("--keep", \%opt_HH); # should we leave intermediates files on disk, instead of removing them?
my @files2rm_A = ();                       # will be filled with files to remove, --keep was not enabled

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
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

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

# copy the ref list to the build directory if --onlybuild
if($onlybuild_mode) { 
  validateFileExistsAndIsNonEmpty($onlybuild_file, "main", $ofile_info_HH{"FH"});
  runCommand("cp $onlybuild_file $ref_list", opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
  addClosedFileToOutputInfo(\%ofile_info_HH, "RefList", $ref_list, 1, "List of reference sequences used to build HMMs");
}
fileLinesToArray($ref_list, 1, \@ref_list_seqname_A, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
my $n_ref = scalar(@ref_list_seqname_A);

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
      DNAORG_FAIL("ERROR in dnaorg_classify.pl::main(), could not find mapping fasta sequence name for list sequence name $ref_list_seqname", 1, $ofile_info_HH{"FH"});
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
        DNAORG_FAIL("ERROR in dnaorg_classify.pl::main(), could not find mapping fasta sequence name for list sequence name $cls_list_seqname", 1, $ofile_info_HH{"FH"});
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
  
  cmscanOrNhmmscanWrapper(\%execs_H, 0, $out_root, $cls_fa, $n_cls_seqs, $tblout_file, $progress_w, 
                          \@mdl_file_A, undef, \%opt_HH, \%ofile_info_HH); 
  # in above cmscanOrNhmmscanWrapper call: '0' means run nhmmscan, not cmscan, 'undef' is for the model length array, irrelevant b/c we're using nhmmscan

  ########################################################################################################################################################
  #
  # PARSE nhmmscan files
  # loop through sequences
  # loop through lines
  #
  # Chooses a RefSeq for each sequence
  # Creates: .matches.info file giving information on each match
  # Creates: ntlists for each RefSeq
  #
  ########################################################################################################################################################

  $start_secs = outputProgressPrior("Creating match info file", $progress_w, $log_FH, *STDOUT);

  # Generate match information file header                                                                               
  my $match_file = $out_root . ".matches.info";
  open(MATCH_INFO, ">", $match_file) || fileOpenFailure($match_file, $0, $!, "writing", $ofile_info_HH{"FH"});
  print MATCH_INFO "########################################################################################################################################\n";
  print MATCH_INFO "#\n";
  print MATCH_INFO "# Query:       Accession number of the sequence \n";
  print MATCH_INFO "# RefSeq:      The RefSeq that the sequence was assigned to\n";
  print MATCH_INFO "# Bit score:   Bit score of hit to 'RefSeq'\n";
  print MATCH_INFO "# E-val:       E-val of hit to 'RefSeq'\n";
  print MATCH_INFO "# Coverage:    The percentage of the query that the hit to 'RefSeq' covers (Hit length)/(Query length)\n";
  print MATCH_INFO "# Bias:        TODO\n";
  print MATCH_INFO "# # Hits:      The number of individual hits to this RefSeq (the program combines stats such as bit score and covg. from separate hits)\n";
  print MATCH_INFO "# H2: RefSeq:  The RefSeq that had the second strongest hit\n";
  print MATCH_INFO "# Bit Diff:    The amount by which the bit score of the 'RefSeq' hit is greater than that of the 'H2: RefSeq' hit\n";
  print MATCH_INFO "# Covg. Diff:  The amount by which the coverage of the 'RefSeq' hit is greater than that of the 'H2: RefSeq' hit\n";
#  print MATCH_INFO "# Num. Correct Hits: The amount of times 'Exp. RefSeq' produced a hit\n";
  print MATCH_INFO "#\n";
  print MATCH_INFO "########################################################################################################################################\n";
  print MATCH_INFO "\n";
  print MATCH_INFO "\n";
  print MATCH_INFO "#H ";
  my @header = ("Query   ","RefSeq   ","Bit score","E-val","Coverage","Bias","# Hits","H2: RefSeq","Bit Diff","Covg. Diff \n");
  print MATCH_INFO join("\t\t", @header);
  print MATCH_INFO "# ";
  print MATCH_INFO "--------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
  print MATCH_INFO "\n";
  
  # Generate data structures to build ntlists from
  my %ntlist_HA = (); # Hash of arrays containing each RefSeq's ntlist
                      # Key:    RefSeq accession #
                      # Value:  Array of accessions for which have been assigned to this RefSeq
  foreach (@ref_list_seqname_A) {
    @{$ntlist_HA{$_}} = (); # initialize each RefSeq's value to an empty array
  }
  
  # initialize the non-assigned list to an empty array
  @{$ntlist_HA{"non-assigned"}} = ();
  
  # Generate data structures to build a sequence profile for further evaluation                       
  my %hit_info_HAA = (); # Hash of hash of arrays containing nhmmscan info for each hit for each sequence
                         # 1D Key:    Sequence accession #                                                                             
                         # 1D Array:  Order of hit (1st nhmmscan hit = index 0, 2nd = index 1, etc)
                         # Value:     Array w/ relevant nhmmscan info for the hit (see comments below)
  
  my @hit_output_A = (); # array containing all the output for this hit
  my %counter_H    = (); # key is sequence name, value is number of hits for this sequence
  my $h;                 # counter over hits
  my $i;                 # counter over sequences

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
      # index         feature
      #--------------------------
      # 0             Query accn
      # 1             RefSeq accn
      # 2             Bit-score
      # 3             E-value
      # 4             Coverage (Hit length)/(Query length)
      # 5             Bias
      # 6             Number of hits to this RefSeq
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
      my $hit_len = $hit_specs_A[7] - $hit_specs_A[6];
      #printf("HEYA ref:$ref_list_seqname seq:$cls_fasta_seqname seqlen:$cls_fasta_seqlen hitlen:$hit_len\n");        

      # if this is the first hit to this RefSeq
      if(! defined($first_index)) {
        #debug
        #print "\tfirst occurence of $cls_fasta_seqname \n";
          
        @hit_output_A = (); # Array to be outputed to match info file                                                                                                                          
        push(@hit_output_A, $cls_fasta_seqname);
        push(@hit_output_A, $ref_list_seqname);
          
        push(@hit_output_A, $hit_specs_A[13]); # add bit score
        my $e_val = sprintf("%.7f", $hit_specs_A[12]);
        push(@hit_output_A, $e_val); # add E-val                                                                                                            
          
        # add coverage				
        my $coverage = $hit_len/$cls_fasta_seqlen;
        $coverage = sprintf("%.7f", $coverage);
        push(@hit_output_A, $coverage);
          
        push(@hit_output_A, $hit_specs_A[14]);  # add bias                                                                                                                            
        push(@hit_output_A, 1);                 # initialize 'Number of hits' to 1
          
          
        @{$hit_info_HAA{$cls_fasta_seqname}}[$counter_H{$cls_fasta_seqname}] = ();
        @{$hit_info_HAA{$cls_fasta_seqname}[$counter_H{$cls_fasta_seqname}]} = @hit_output_A;
          
        $counter_H{$cls_fasta_seqname}++;
          
      }else {    # if this is a hit to a sequence that has already been hit
        #@hit_output_A = \@{$hit_info_HAA{$cls_fasta_seqname}[$first_index]};
          
        #debug
        #print "\tsecond, third, etc hit of $cls_fasta_seqname \n";
          
        # TODO - check with Eric!
        if($hit_specs_A[13] > 0) { # only add hits with positive bit scores
          @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[2] += $hit_specs_A[13]; # Add bit score
###        @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[3] += $hit_specs_A[12]; # Add E-val
          @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[5] += $hit_specs_A[14]; # Add bias
          @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[6] ++;                  # increase 'Number of hits' by 1

          ##############################################
          # Calculate and add coverage
          my $prev_covg = @{$hit_info_HAA{$cls_fasta_seqname}[$first_index]}[4];
          my $prev_hit_len = $prev_covg*$cls_fasta_seqlen;
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
      DNAORG_FAIL("ERROR in dnaorg_classify.pl::main(), could not find mapping fasta sequence name for list sequence name $cls_list_seqname", 1, $ofile_info_HH{"FH"});
    }
    if(exists($hit_info_HAA{$cls_fasta_seqname})) { 
      # there was at least one valid hit for this sequence
      # parse and output the best hit's info

      # sort order of hits by bit score
      @{$hit_info_HAA{$cls_fasta_seqname}} = sort { $a->[2] <=> $b->[2] } @{$hit_info_HAA{$cls_fasta_seqname}};
      
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
      $ref_list_seqname = $hit_output_A[1];
      
      ##############################################################################
      # Add comparison indicies to @hit_output_A
      #
      # index        feature
      # ----------------------
      # 7            RefSeq accn # for hit 2
      # 8            difference between bit scores of hit 1 and hit 2
      # 9            difference between coverage of hit 1 and hit 2
      ##############################################################################
      
      # TODO - hardcode-y 
      if( scalar(@{$hit_info_HAA{$cls_fasta_seqname}}) >= 2) {     # if there is more than one hit  
        push(@hit_output_A, @{@{$hit_info_HAA{$cls_fasta_seqname}}[-2]}[1]);
	
        my $bit_diff = $hit_info_HAA{$cls_fasta_seqname}[-1][2] - $hit_info_HAA{$cls_fasta_seqname}[-2][2];
        $bit_diff = sprintf("%9.1f", $bit_diff);
        my $cov_diff = $hit_info_HAA{$cls_fasta_seqname}[-1][4] - $hit_info_HAA{$cls_fasta_seqname}[-2][4];
        $cov_diff = sprintf("%.7f", $cov_diff);
        push(@hit_output_A, $bit_diff);
        push(@hit_output_A, $cov_diff);
	
      } else { # if there was only one hit, there's no comparison
        push(@hit_output_A, "-----");
        push(@hit_output_A, "-----");
        push(@hit_output_A, "-----");
      }
      
      # add hit info to .matches.info file
      print MATCH_INFO join("\t\t", @hit_output_A);
      print MATCH_INFO "\n";
      
      # Add this seq to the hash key that corresponds to it's RefSeq
      push(@{$ntlist_HA{$ref_list_seqname}}, $cls_list_seqname);
      # IMPORTANT: push $cls_list_seqname, and not $cls_fasta_seqname (cls_list_seqname is from the input list 
      # (unless --infasta) and may not include the version
    } else { # if there were no eligible hits
      @hit_output_A = ($cls_fasta_seqname,"--------","--------","--------","--------","--------","--------","--------","--------","--------");
      print MATCH_INFO join("\t\t", @hit_output_A);
      print MATCH_INFO "\n";
      
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
  push(@tmp_output_A, sprintf("%-20s  %10s  %10s\n", "# RefSeq-accession", "num-seqs", "fract-seqs"));
  push(@tmp_output_A, sprintf("%-20s  %10s  %10s\n", "#-------------------", "----------", "----------"));
  foreach my $ref_list_seqname (@ref_list_seqname_A) {
    my $ntlist_file    = $out_root . ".$ref_list_seqname.ntlist";
    my $sub_fasta_file = $out_root . ".$ref_list_seqname.fa";
    my $cur_nseq = scalar(@{$ntlist_HA{$ref_list_seqname}});
    push(@tmp_output_A, sprintf("%-20s  %10d  %10.4f\n", $ref_list_seqname, $cur_nseq, $cur_nseq / $n_cls_seqs));
    if($cur_nseq > 0) { 
      open(NTLIST, "> $ntlist_file")  || fileOpenFailure($ntlist_file, $0, $!, "writing", $ofile_info_HH{"FH"});
      print NTLIST "$ref_list_seqname\n";
      foreach (@{$ntlist_HA{$ref_list_seqname}}) {
        print NTLIST "$_\n";
      }
      close(NTLIST);
      addClosedFileToOutputInfo(\%ofile_info_HH, "$ref_list_seqname.ntlist", $ntlist_file, 1, "ntlist for $ref_list_seqname");

      if($infasta_mode) { 
        # fetch the sequences into a new fasta file
        sleep(0.1); # make sure that NTLIST is closed
        $cmd  = "tail -n +2 $ntlist_file |" . $execs_H{"esl-sfetch"} . " -f $cls_fa - > $sub_fasta_file"; # fetch the sequences
        runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
        addClosedFileToOutputInfo(\%ofile_info_HH, "$ref_list_seqname.fa", $sub_fasta_file, 1, "fasta file with sequences assigned to $ref_list_seqname");
      }
    }
  }
  
  # Generate file that lists non-assigned sequences
  my $non_assigned_file = $out_root . ".non-assigned";
  
  my $cur_nseq = scalar(@{$ntlist_HA{"non-assigned"}});
  push(@tmp_output_A, sprintf("%-20s  %10d  %10.4f\n", "NON-ASSIGNED", $cur_nseq, $cur_nseq / $n_cls_seqs));
  if($cur_nseq > 0) { 
    open(NALIST, "> $non_assigned_file")  || fileOpenFailure($non_assigned_file, $0, $!, "writing", $ofile_info_HH{"FH"});
    foreach (@{$ntlist_HA{"non-assigned"}}) {
      print NALIST "$_\n";
    }
    addClosedFileToOutputInfo(\%ofile_info_HH, "non-assigned", $non_assigned_file, 1, "List of sequences not assigned to a RefSeq");
  }
  
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

} # end of 'else' entered if $onlybuild_mode is FALSE
##########################################################################################################################################################

##########
# Conclude
##########

# if not --keep, then remove unnessecary files (while checking to make sure that you're not removing the same file twice
my %seen;

foreach my $file2rm (@files2rm_A) {
  next if $seen{$file2rm}++;
  runCommand("rm $file2rm", 0, $ofile_info_HH{"FH"});
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
    print SCRIPT "\#\$  -l h_rt=288000,h_vmem=32G,mem_free=32G,reserve_mem=32G  \n";
    # old version of above line: -l h_vmem=32G,reserve_mem=32G,mem_free=32G
    # new version of above line: -l h_rt=288000,h_vmem=32G,mem_free=32G,reserve_mem=32G
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

