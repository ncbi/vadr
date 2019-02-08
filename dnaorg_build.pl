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
my $inf_exec_dir      = $dnaorgdir . "/infernal-dev/src/";
my $hmmer_exec_dir    = $dnaorgdir . "/hmmer-3.1b2/src/";
my $esl_exec_dir      = $dnaorgdir . "/infernal-dev/easel/miniapps/";
my $esl_fetch_cds     = $dnaorgdir . "/esl-fetch-cds/esl-fetch-cds.pl";
my $blast_exec_dir    = "/usr/bin/"; # HARD-CODED FOR NOW

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
#     option            type       default               group   requires incompat    preamble-output                                 help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                          "display this help",                                  \%opt_HH, \@opt_order_A);
opt_Add("-c",           "boolean", 0,                        1,    undef, undef,      "genome is circular",                           "genome is circular",                                 \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",                  "force; if dir <reference accession> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-n",           "integer", "4",                      1,    undef, undef,      "set number of CPUs for calibration to <n>",    "for non-big models, set number of CPUs for calibration to <n>", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                   "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--dirout",     "string",  undef,                    1,    undef, undef,      "output directory specified as <s>",            "specify output directory as <s>, not <ref accession>", \%opt_HH, \@opt_order_A);
opt_Add("--matpept",    "string",  undef,                    1,    undef, undef,      "using pre-specified mat_peptide info",         "read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships", \%opt_HH, \@opt_order_A);
opt_Add("--nomatpept",  "boolean", 0,                        1,    undef,"--matpept", "ignore mat_peptide annotation",                "ignore mat_peptide information in reference annotation", \%opt_HH, \@opt_order_A);
opt_Add("--xfeat",      "string",  undef,                    1,    undef, undef,      "build models of additional qualifiers",        "build models of additional qualifiers in string <s>", \%opt_HH, \@opt_order_A);  
opt_Add("--dfeat",      "string",  undef,                    1,    undef, undef,      "annotate additional qualifiers as duplicates", "annotate qualifiers in <s> from duplicates (e.g. gene from CDS)",  \%opt_HH, \@opt_order_A);  
opt_Add("--keep",       "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",           "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"2"} = "optional output files";
#       option       type       default                group  requires incompat  preamble-output                          help-output    
opt_Add("--mdlinfo",    "boolean", 0,                        2,    undef, undef, "output internal model information",     "create file with internal model information",   \%opt_HH, \@opt_order_A);
opt_Add("--ftrinfo",    "boolean", 0,                        2,    undef, undef, "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"3"} = "options for skipping stages and using files from an earlier, identical run, primarily useful for debugging";
#     option               type       default               group   requires    incompat                  preamble-output                                            help-output    
opt_Add("--skipedirect",   "boolean", 0,                       3,   undef,      undef,                    "skip the edirect steps, use existing results",           "skip the edirect steps, use data from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipfetch",     "boolean", 0,                       3,   undef,      undef,                    "skip the sequence fetching steps, use existing results", "skip the sequence fetching steps, use files from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipbuild",     "boolean", 0,                       3,   undef,      undef,                    "skip the build/calibrate steps",                         "skip the model building/calibrating, requires --mdlinfo and/or --ftrinfo", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"4"} = "options for building models for origin sequences";
#     option               type       default               group   requires               incompat                  preamble-output                                              help-output    
opt_Add("--orginput",      "string",  undef,                   4,   "-c,--orgstart,--orglen",  undef,                   "read training alignment for origin sequence from file <s>", "read training alignment for origin sequences from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--orgstart",      "integer", 0,                       4,   "-c,--orginput,--orglen",  undef,                   "origin sequence starts at position <n>",                    "origin sequence starts at position <n> in file <s> from --orginput <s>", \%opt_HH, \@opt_order_A);
opt_Add("--orglen",        "integer", 0,                       4,   "-c,--orginput,--orgstart",undef,                   "origin sequence is <n> nucleotides long",                   "origin sequence is <n> nucleotides long", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_build.pl [-options] <reference accession>\n";
my $synopsis = "dnaorg_build.pl :: build homology models for features of a reference sequence";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'c'            => \$GetOptions_H{"-c"},
                'f'            => \$GetOptions_H{"-f"},
                'n=s'          => \$GetOptions_H{"-n"},
                'v'            => \$GetOptions_H{"-v"},
                'dirout=s'     => \$GetOptions_H{"--dirout"},
                'matpept=s'    => \$GetOptions_H{"--matpept"},
                'nomatpept'    => \$GetOptions_H{"--nomatpept"},
                'xfeat=s'      => \$GetOptions_H{"--xfeat"},
                'dfeat=s'      => \$GetOptions_H{"--dfeat"},
                'keep'         => \$GetOptions_H{"--keep"},
# optional output files
                'mdlinfo'      => \$GetOptions_H{"--mdlinfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"},
# options for skipping stages, using earlier results
                'skipedirect'  => \$GetOptions_H{"--skipedirect"},
                'skipfetch'    => \$GetOptions_H{"--skipfetch"},
                'skipbuild'    => \$GetOptions_H{"--skipbuild"},
# options for building models for the origin sequence
                'orginput=s'   => \$GetOptions_H{"--orginput"},
                'orgstart=s'   => \$GetOptions_H{"--orgstart"},
                'orglen=s'     => \$GetOptions_H{"--orglen"});

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.45";
my $releasedate   = "Feb 2019";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date, $dnaorgdir);
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
my $do_origin  = opt_IsUsed("--matpept", \%opt_HH);

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
my @arg_desc_A = ("reference accession");
my @arg_A      = ($ref_accn);
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
outputBanner($log_FH, $version, $releasedate, $synopsis, $date, $dnaorgdir);
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
  my $matpept_optfile = opt_Get("--matpept", \%opt_HH);
  my $dest_matpept_optfile = $out_root . ".matpept";
  parseMatPeptSpecFile($matpept_optfile, \@cds2pmatpept_AA, \@cds2amatpept_AA, $ofile_info_HH{"FH"});
  # copy the matpept file to a special file name
  runCommand("cp $matpept_optfile $dest_matpept_optfile", opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});
}

###################################################
# make sure the required executables are executable
###################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmbuild"}       = $inf_exec_dir . "cmbuild";
$execs_H{"esl-reformat"}  = $esl_exec_dir . "esl-reformat";
$execs_H{"esl_fetch_cds"} = $esl_fetch_cds;
$execs_H{"makeblastdb"}   = $blast_exec_dir . "makeblastdb";
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################################################################
# Step 1. Output the consopts file that dnaorg_annotate.pl will use to 
#         make sure options used are consistent between dnaorg_build.pl and 
#         dnaorg_annotate.pl 
###########################################################################
my $progress_w = 80; # the width of the left hand column in our progress output, hard-coded
my $start_secs = outputProgressPrior("Outputting information on options used for future use with dnaorg_annotate.pl", $progress_w, $log_FH, *STDOUT);
output_consopts_file($out_root . ".consopts", \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###########################################################################
# Step 2. Validate the origin alignment file (if --orginput used).
###########################################################################
if(opt_IsUsed("--orginput", \%opt_HH)) { 
  $progress_w = 80; # the width of the left hand column in our progress output, hard-coded
  $start_secs = outputProgressPrior("Validating and processing origin input alignment", $progress_w, $log_FH, *STDOUT);
  process_origin_input_alignment($out_root, \%opt_HH, \%ofile_info_HH);
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

###########################################################################
# Step 3. Gather and process information on reference genome using Edirect.
###########################################################################
$start_secs = outputProgressPrior("Gathering information on reference using edirect", $progress_w, $log_FH, *STDOUT);

my %cds_tbl_HHA = ();    # CDS data from .cds.tbl file, hash of hashes of arrays, 
                         # 1D: key: accession
                         # 2D: key: column name in gene ftable file
                         # 3D: per-row values for each column
my %mp_tbl_HHA = ();     # mat_peptide data from .matpept.tbl file, hash of hashes of arrays, 
                         # 1D: key: accession
                         # 2D: key: column name in gene ftable file
                         # 3D: per-row values for each column
my %xfeat_tbl_HHHA = (); # xfeat (extra feature) data from feature table file, hash of hash of hashes of arrays
                         # 1D: qualifier name, e.g. 'gene'
                         # 2D: key: accession
                         # 3D: key: column name in gene ftable file
                         # 4D: per-row values for each column
my %dfeat_tbl_HHHA = (); # dfeat (duplicate feature) data from feature table file, hash of hash of hashes of arrays
                         # 1D: qualifier name, e.g. 'gene'
                         # 2D: key: accession
                         # 3D: key: column name in gene ftable file
                         # 4D: per-row values for each column
my %seq_info_HA = ();    # hash of arrays, avlues are arrays [0..$nseq-1];
                         # 1st dim keys are "seq_name", "accn_name", "seq_len", "accn_len".
                         # $seq_info_HA{"accn_name"}[0] is our reference accession
@{$seq_info_HA{"accn_name"}} = ($ref_accn);

# parse --xfeat option if necessary and initiate hash of hash of arrays for each comma separated value
my $do_xfeat = 0;
if(opt_IsUsed("--xfeat", \%opt_HH)) { 
  $do_xfeat = 1;
  my $xfeat_str = opt_Get("--xfeat", \%opt_HH);
  foreach my $xfeat (split(",", $xfeat_str)) { 
    %{$xfeat_tbl_HHHA{$xfeat}} = ();
  }
}
# parse --dfeat option if necessary
my $do_dfeat = 0;
if(opt_IsUsed("--dfeat", \%opt_HH)) { 
  $do_dfeat = 1;
  my $dfeat_str = opt_Get("--dfeat", \%opt_HH);
  foreach my $dfeat (split(",", $dfeat_str)) { 
    %{$dfeat_tbl_HHHA{$dfeat}} = ();
    if(exists $xfeat_tbl_HHHA{$dfeat}) {
      DNAORG_FAIL("ERROR, with --xfeat <s1> and --dfeat <s2>, no qualifier names can be in common between <s1> and <s2>, found $dfeat", 1, $ofile_info_HH{"FH"});
    }
  }
}

# Call the wrapper function that does the following:
#  1) creates the edirect .mat_peptide file, if necessary
#  2) creates the edirect .ftable file
#  3) creates the length file
#  4) parses the edirect .mat_peptide file, if necessary
#  5) parses the edirect .ftable file
#  6) parses the length file
wrapperGetInfoUsingEdirect(undef, $ref_accn, $out_root, \%cds_tbl_HHA, \%mp_tbl_HHA, \%xfeat_tbl_HHHA, \%dfeat_tbl_HHHA,
                           \%seq_info_HA, \%ofile_info_HH, 
                           \%opt_HH, $ofile_info_HH{"FH"}); # 1st argument is undef because we are only getting info for $ref_accn

if($do_matpept) {  
  # validate the CDS:mat_peptide relationships that we read from the $matpept input file
  matpeptValidateCdsRelationships(\@cds2pmatpept_AA, \%{$cds_tbl_HHA{$ref_accn}}, \%{$mp_tbl_HHA{$ref_accn}}, opt_Get("-c", \%opt_HH), $seq_info_HA{"accn_len"}[0], $ofile_info_HH{"FH"});
}
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################
# Step 4. Fetch and process the reference genome sequence
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
                                                    ($do_xfeat)   ? \%xfeat_tbl_HHHA  : undef,
                                                    ($do_dfeat)   ? \%dfeat_tbl_HHHA  : undef,
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

if(exists $ofile_info_HH{"FH"}{"mdlinfo"}) { 
  dumpInfoHashOfArrays("Model information (%mdl_info_HA)", 0, \%mdl_info_HA, $ofile_info_HH{"FH"}{"mdlinfo"});
}
if(exists $ofile_info_HH{"FH"}{"ftrinfo"}) { 
  dumpInfoHashOfArrays("Feature information (%ftr_info_HA)", 0, \%ftr_info_HA, $ofile_info_HH{"FH"}{"ftrinfo"});
}

#######################################################################
# Step 4B. Fetch the CDS protein translations and build BLAST database
#######################################################################
$start_secs = outputProgressPrior("Fetching protein translations of CDS and building BLAST DB", $progress_w, $log_FH, *STDOUT);
my @prot_fa_file_A = ();
fetch_proteins_into_fasta_files($out_root, $ref_accn, \%ftr_info_HA, \@prot_fa_file_A, \%opt_HH, \%ofile_info_HH);

foreach my $prot_fa_file (@prot_fa_file_A) { 
  create_blast_protein_db(\%execs_H, $prot_fa_file, \%opt_HH, \%ofile_info_HH);
}
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################
# Step 5. Build the model
#########################
if(! opt_Get("--skipbuild", \%opt_HH)) { 
  $start_secs = outputProgressPrior("Building model (this could take a while)", $progress_w, $log_FH, *STDOUT);

  my $cmbuild_cmd  = $execs_H{"cmbuild"} . " --verbose -F $out_root.cm " . $ofile_info_HH{"fullpath"}{"refstk"} . " > $out_root.cmbuild";
  runCommand($cmbuild_cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});

  outputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);

  addClosedFileToOutputInfo(\%ofile_info_HH, "cm", "$out_root.cm", 1, "CM file");
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
#            "--xfeat", "--dfeat", and "-c".
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
  if(opt_IsUsed("--xfeat", $opt_HHR)) { 
    printf OUT ("--xfeat %s\n", opt_Get("--xfeat", $opt_HHR));
    $printed_any_options = 1;
  }
  if(opt_IsUsed("--dfeat", $opt_HHR)) { 
    printf OUT ("--dfeat %s\n", opt_Get("--dfeat", $opt_HHR));
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
# Subroutine: process_origin_input_alignment()
# Incept:     EPN, Thu Jul  7 13:44:15 2016 [origin-hmms-01 branch]
#
# Purpose:   Validate an input alignment file that will be used
#            to identify origin sequences, and then 'process' it
#            by splitting it nearly in half, with an overlap between 
#            the two halves which is exactly the origin sequence.
#            The two halves will each serve as a training alignment for 
#            cmbuild in a later step.
#
# Arguments:
#  $out_root:          root for naming output files
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
sub process_origin_input_alignment {
  my $sub_name = "process_origin_input_alignment";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_root, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience

  my $input_file   = opt_Get("--orginput", $opt_HHR);
  my $origin_start = opt_Get("--orgstart", $opt_HHR);
  my $origin_stop  = $origin_start + opt_Get("--orglen", $opt_HHR) - 1; 
  my $origin_len   = $origin_stop - $origin_start + 1;

  my $msa1_outfile = $out_root . ".origin.5p.stk";
  my $msa2_outfile = $out_root . ".origin.3p.stk";

  # open and validate file
  my $msa1 = Bio::Easel::MSA->new({
    fileLocation => $input_file,
    isDna => 1
                                 });  

  # determine length of the alignment
  my $alen = $msa1->alen();

  if($alen < $origin_stop) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, with --orgstart, stop position of origin in alignment should be %d, but alignment is only %d columns long", $origin_stop, $alen),
                1, $FH_HR); 
  }

  # split the alignment in half, we need two versions of the alignment to do this:
  my $msa2 = Bio::Easel::MSA->new({
    fileLocation => $input_file,
    isDna => 1
                                 });  
  

  my @useme1_A = (); # 0..$i..$alen-1: '1' if column $i+1 should be kept in $msa1, '0' if not
  my @useme2_A = (); # 0..$i..$alen-1: '1' if column $i+1 should be kept in $msa2, '0' if not
  my $i;
  for($i = 0; $i < ($origin_start-1); $i++) { 
    $useme1_A[$i] = 1;
    $useme2_A[$i] = 0;
  }
  for($i = $origin_start-1; $i < $origin_stop; $i++) { 
    $useme1_A[$i] = 1;
    $useme2_A[$i] = 1;
  }
  for($i = $origin_stop; $i < $alen; $i++) { 
    $useme1_A[$i] = 0;
    $useme2_A[$i] = 1;
  }

  $msa1->column_subset_rename_nse(\@useme1_A, 0);
  $msa2->column_subset_rename_nse(\@useme2_A, 0);

  my $name1 = removeDirPath($out_root);
  my $name2 = removeDirPath($out_root);
  $name1 .= sprintf(".origin.5p.len%d", $origin_len);
  $name2 .= sprintf(".origin.3p.len%d", $origin_len);
  
  $msa1->set_name($name1);
  $msa2->set_name($name2);

  $msa1->write_msa($msa1_outfile);
  $msa2->write_msa($msa2_outfile);


  addClosedFileToOutputInfo($ofile_info_HHR, "origin5p", "$msa1_outfile", 1, sprintf("Subset of alignment in %s containing positions %d to %d. The origin and the flanking 5 prime sequence.", $input_file, 1, $origin_stop));
  addClosedFileToOutputInfo($ofile_info_HHR, "origin3p", "$msa2_outfile", 1, sprintf("Subset of alignment in %s containing positions %d to %d. The origin and the flanking 3 prime sequence.", $input_file, $origin_start, $alen));

  return;
}

#################################################################
# Subroutine: fetch_proteins_into_fasta_files()
# Incept:     EPN, Wed Oct  3 16:10:26 2018
# 
# Purpose:    Fetch the protein translations of CDS for the genome
#             and create multiple N+1 FASTA files, one with each
#             single sequence (N) and one with all sequences.
#             Fill @{$fa_file_AR} with the sequence file names.
#
# Arguments:
#   $out_root:       string for naming output files
#   $ref_accn:       reference accession
#   $ftr_info_HAR:   REF to the feature info, pre-filled
#   $fa_file_AR:     REF to array of fasta file names, filled here 
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $ofile_info_HHR: REF to the 2D hash of output file information
#                    
# Returns: void
# Dies:    if a fetched location for a feature does not match to any feature's "ref_coords" 
#
#################################################################
sub fetch_proteins_into_fasta_files { 
  my $sub_name = "fetch_proteins_into_fasta_files";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_root, $ref_accn, $ftr_info_HAR, $fa_file_AR, $opt_HHR, $ofile_info_HHR) = @_;
  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience

  my $efetch_out_file  = $out_root . ".prot.efetch";
  my $all_fa_out_file  = $out_root . ".prot.fa";
  runCommand("esearch -db nuccore -query $ref_accn | efetch -format gpc | xtract -insd CDS protein_id INSDFeature_location translation > $efetch_out_file", opt_Get("-v", $opt_HHR), 0, $FH_HR); 
  # NOTE: could get additional information to add to fasta defline, e.g. add 'product' after 'translation' above.

  # parse that file to create the fasta files
  open(IN,         $efetch_out_file) || fileOpenFailure($efetch_out_file,  $sub_name, $!, "reading", $FH_HR);
  open(ALLFA, ">", $all_fa_out_file) || fileOpenFailure($all_fa_out_file,      $sub_name, $!, "writing", $FH_HR);
  while(my $line = <IN>) { 
    chomp $line;
    my @el_A = split(/\t/, $line);
    if(scalar(@el_A) != 4) { 
      DNAORG_FAIL("ERROR in $sub_name, not exactly 4 tab delimited tokens in efetch output file line\n$line\n", 1, $ofile_info_HHR->{"FH"});
    }
    my ($read_ref_accn, $prot_accn, $location, $translation) = (@el_A);
    my $new_name = $prot_accn . "/" . $location;

    print ALLFA  (">$new_name\n$translation\n");

    # determine what feature this corresponds to, and create the individual fasta file for that
    my $ftr_idx = blastxDbSeqNameToFtrIdx($new_name, $ftr_info_HAR, $ofile_info_HHR->{"FH"}); # this will die if we can't find the feature with $location
    my $indi_fa_out_file = $out_root . ".f" . $ftr_idx . ".prot.fa";
    open(INDIFA, ">", $indi_fa_out_file) || fileOpenFailure($indi_fa_out_file, $sub_name, $!, "writing", $FH_HR);
    print INDIFA (">$new_name\n$translation\n");
    close(INDIFA);
    addClosedFileToOutputInfo($ofile_info_HHR, "prot-indi-f" . $ftr_idx . "-fa", $indi_fa_out_file, 0, "protein FASTA file with proteins for feature $ftr_idx");
    push(@{$fa_file_AR}, $indi_fa_out_file);
  }
  close(IN);
  close(ALLFA);

  addClosedFileToOutputInfo($ofile_info_HHR, "prot-all-fa", $all_fa_out_file, 0, "protein FASTA file with proteins for all features");
  push(@{$fa_file_AR}, $all_fa_out_file);

  addClosedFileToOutputInfo($ofile_info_HHR, "prot-fetch",   $efetch_out_file,  0, "efetch output with protein information");

  return;
}

#################################################################
# Subroutine: create_blast_protein_db
# Incept:     EPN, Wed Oct  3 16:31:38 2018
# 
# Purpose:    Create a protein blast database from a fasta file.
#
# Arguments:
#   $execs_HR:       reference to hash with infernal executables, 
#                    e.g. $execs_HR->{"cmcalibrate"} is path to cmcalibrate, PRE-FILLED
#   $prot_fa_file:   FASTA file of protein sequences to make blast db from
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $ofile_info_HHR: REF to the 2D hash of output file information
#                    
# Returns:    void
#
#################################################################
sub create_blast_protein_db {
  my $sub_name = "create_blast_protein_db";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $prot_fa_file, $opt_HHR, $ofile_info_HHR) = @_;

  runCommand($execs_HR->{"makeblastdb"} . " -in $prot_fa_file -dbtype prot > /dev/null", opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});

  return;
}
