#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Thu Feb 18 12:48:16 2016 [dnaorg_annotate.pl split off from dnaorg_annotate_genomes.pl]
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
#   - output program banner and open output files
#   - parse the optional input files, if necessary
#   - make sure the required executables are executable
#
# Step 1. Gather and process information on reference genome using Edirect
#
# Step 2. Fetch and process the reference genome sequence
#
#######################################################################################

# hard-coded-paths:
my $inf_exec_dir   = "/usr/local/infernal/1.1.1/bin/";
my $esl_exec_dir   = "/usr/local/infernal/1.1.1/bin/";
my $esl_fetch_cds  = "/panfs/pan1/dnaorg/programs/esl-fetch-cds.pl";

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
#     option            type       default               group   requires incompat    preamble-output                               help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                        "display this help",                                  \%opt_HH, \@opt_order_A);
opt_Add("-c",           "boolean", 0,                        1,    undef, undef,      "genome is circular",                         "genome is circular",                                 \%opt_HH, \@opt_order_A);
opt_Add("-d",           "string",  undef,                    1,    undef, undef,      "directory specified as",                     "specify output directory is <s1> (created with dnaorg_build.pl -d <s>), not <ref accession>", \%opt_HH, \@opt_order_A);
#opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",           "force; if dir <reference accession> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                 "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--matpept",    "string",  undef,                    1,    undef, undef,      "using pre-specified mat_peptide info",       "read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships", \%opt_HH, \@opt_order_A);
opt_Add("--nomatpept",  "boolean", 0,                        1,    undef,"--matpept", "ignore mat_peptide annotation",              "ignore mat_peptide information in reference annotation", \%opt_HH, \@opt_order_A);
opt_Add("--specstart",  "string",  undef,                    1,    undef, undef,      "using pre-specified alternate start codons", "read specified alternate start codons per CDS from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--dirty",      "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",         "do not remove intermediate files, leave them all on disk", \%opt_HH, \@opt_order_A);
opt_Add("--model",      "string",  undef,                    1,    undef, undef,      "use model in file",                          "use model file <s>", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_annotate.pl [-options] <file with list of accessions to annotate>\n";
my $synopsis = "dnaorg_annotate.pl :: annotate sequences based on a reference annotation";

my $options_okay = 
    &GetOptions('h'             => \$GetOptions_H{"-h"}, 
# basic options
                'c'             => \$GetOptions_H{"-c"},
                'd=s'           => \$GetOptions_H{"-d"},
                'f'             => \$GetOptions_H{"-f"},
                'v'             => \$GetOptions_H{"-v"},
                'matpept=s'     => \$GetOptions_H{"--matpept"},
                'nomatpept'     => \$GetOptions_H{"--nomatpept"},
                'specstart=s'   => \$GetOptions_H{"--specstart"},
                'dirty'         => \$GetOptions_H{"--dirty"},
                'model'         => \$GetOptions_H{"--model"});

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.1";
my $releasedate   = "Feb 2016";

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
  print "\nTo see more help on available options, do dnaorg_annotate.pl -h\n\n";
  exit(1);
}
my ($listfile) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

my $dir        = opt_Get("-d", \%opt_HH);          # this will be undefined unless -d set on cmdline
my $do_matpept = opt_IsOn("--matpept", \%opt_HH);  # this will be '0' unless --matpept set on cmdline 

###############
# Preliminaries
###############
# first, parse the list file, we need to do this first because we need
# to know what the reference accession <refaccn> is to check if the
# directory <refaccn> exists
my @accn_A = (); # array of accessions, $accn_A[0] is our reference
parseListFile($listfile, 1, \@accn_A, undef); # 1 
my $ref_accn = $accn_A[0];

my $dir_set_as_ref_accn = 0;
if(! defined $dir) { 
  $dir = $ref_accn;
  $dir_set_as_ref_accn = 1;
}

# make sure that $dir exists
if(! -d $dir) {
  DNAORG_FAIL(sprintf("ERROR, directory $dir %s does not exist", $dir_set_as_ref_accn ? "(first accession read from $listfile)" : "(specified with -d)"), 1, undef);
}

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root = $dir . "/" . $dir_tail . ".dnaorg_annotate";

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("file with list of accessions");
my @arg_A      = ($listfile);
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
openAndAddFileToOutputInfo(\%ofile_info_HH, "log", $out_root . ".log", "Output printed to screen");
openAndAddFileToOutputInfo(\%ofile_info_HH, "cmd", $out_root . ".cmd", "List of executed commands");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# now we have the log file open, output the banner there too
outputBanner($log_FH, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# now we have the log file open, output the banner there too
outputBanner($log_FH, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

########################################
# parse the optional input files, if nec
########################################
# --matpept <f>
my @cds2pmatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of primary matpept indices that comprise this CDS
my @cds2amatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of all     matpept indices that comprise this CDS
if($do_matpept) { 
  parseMatPeptSpecFile(opt_Get("--matpept", \%opt_HH), \@cds2pmatpept_AA, \@cds2amatpept_AA, $ofile_info_HH{"FH"});
}
# --specstart <f>
my @specstart_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of allowed start codons for this CDS
if(opt_IsOn("--specstart", \%opt_HH)) { 
  parseSpecStartFile(opt_Get("--specstart", \%opt_HH), \@specstart_AA, $ofile_info_HH{"FH"});
}

###################################################
# make sure the required executables are executable
###################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmsearch"}      = $inf_exec_dir . "cmsearch";
$execs_H{"cmalign"}       = $inf_exec_dir . "cmalign";
$execs_H{"cmpress"}       = $inf_exec_dir . "cmpress";
$execs_H{"esl-reformat"}  = $esl_exec_dir . "esl-reformat";
$execs_H{"esl_fetch_cds"} = $esl_fetch_cds;
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################################################################
# Step 1. Gather and process information on reference genome using Edirect.
###########################################################################
my $cmd;             # a command to run with runCommand()
my $progress_w = 60; # the width of the left hand column in our progress output, hard-coded
my $start_secs = outputProgressPrior(sprintf("Gathering information on %d sequences using edirect", scalar(@accn_A)), $progress_w, $log_FH, *STDOUT);

my %cds_tbl_HHA = ();   # CDS data from .cds.tbl file, hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %mp_tbl_HHA = ();    # mat_peptide data from .matpept.tbl file, hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %totlen_H   = ();    # key: accession, value: total length of the sequence for that accession


# Call the wrapper function that does the following:
#  1) creates the edirect .mat_peptide file, if necessary
#  2) creates the edirect .ftable file
#  3) creates the length file
#  4) parses the edirect .mat_peptide file, if necessary
#  5) parses the edirect .ftable file
#  6) parses the length file
wrapperGetInfoUsingEdirect($listfile, $ref_accn, $out_root, \%cds_tbl_HHA, \%mp_tbl_HHA, \%totlen_H, \%ofile_info_HH,
                           \%opt_HH, $ofile_info_HH{"FH"}); # 1st argument is undef because we are only getting info for $ref_accn

if($do_matpept) {  
    # validate the CDS:mat_peptide relationships that we read from the $matpept input file
    matpeptValidateCdsRelationships(\@cds2pmatpept_AA, \%{$cds_tbl_HHA{$ref_accn}}, \%{$mp_tbl_HHA{$ref_accn}}, opt_Get("-c", \%opt_HH), $totlen_H{$ref_accn}, $ofile_info_HH{"FH"});
}
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################
# Step 2. Fetch and process the reference genome sequence
##########################################################
$start_secs = outputProgressPrior("Fetching and processing the reference genome", $progress_w, $log_FH, *STDOUT);
my %mdl_info_HA = ();          # hash of arrays, values are arrays [0..$nmdl-1];
                               # see dnaorg.pm::validateModelInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAndProcessReferenceSequence()
my %ftr_info_HA = ();          # hash of arrays, values are arrays [0..$nftr-1], 
                               # see dnaorg.pm::validateFeatureInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAndProcessReferenceSequence()

# Call the wrapper function that does the following:
#   1) fetches the reference sequence into a fasta file and indexes that fasta file
#   2) determines reference information for each feature (strand, length, coordinates, product)
#   3) determines type of each reference feature ('cds-mp', 'cds-notmp', or 'mp')
#   4) fetches the reference feature sequences and populate information on the models and features
wrapperFetchAndProcessReferenceSequence(\@accn_A, $out_root, \%cds_tbl_HHA,
                                        ($do_matpept) ? \%mp_tbl_HHA      : undef, 
                                        ($do_matpept) ? \@cds2pmatpept_AA : undef, 
                                        \%totlen_H, \%ofile_info_HH,
                                        \%ftr_info_HA, \%mdl_info_HA, \%execs_H,
                                        \%opt_HH, $ofile_info_HH{"FH"});

# verify our model and feature info hashes are complete, 
# if validateFeatureInfoHashIsComplete() fails then the program will exit with an error message
my $nftr = validateFeatureInfoHashIsComplete(\%ftr_info_HA, undef, $ofile_info_HH{"FH"}); # nftr: number of features
my $nmdl = validateModelInfoHashIsComplete  (\%mdl_info_HA, undef, $ofile_info_HH{"FH"}); # nmdl: number of homology models

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

# verify that we have the model file we need
my $model_file = opt_Get("--model", \%opt_HH);     # this will be undefined unless --model set on cmdline
if(! defined $model_file) { 
  # -model not used, make sure we already have the CM DB file from
  # an earlier dnaorg_build.pl run. 
  $model_file = $out_root . ".dnaorg_build.ref.cm";
  if(! -s $model_file) { 
    # the model file does not (yet) exist. This is probably the first
    # time we've run dnaorg_annotate.pl and we have several individual
    # CM files because each was calibrated independently on the
    # farm. In this case, we concatenate the individual files to
    # create one CM database.
    $start_secs = outputProgressPrior("Creating CM database by concatenating individual CM files", $progress_w, $log_FH, *STDOUT);
    concatenate_individual_cm_files($model_file, $out_root, \%opt_HH, \%ofile_info_HH);
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
}

# validate that we have the CM binary index files from cmpress that we need, 
# if we don't, then create them
my $do_press = 0;
for my $suffix ("i1m", "i1i", "i1f", "i1p") { 
  my $file = $model_file . "." . $suffix;
  if(! -s $file) { $do_press = 1; }
}
if($do_press) { 
  # run cmpress to create the CM binary index files
  $start_secs = outputProgressPrior("Preparing the CM database for homology search using cmpress", $progress_w, $log_FH, *STDOUT);
  press_cm_database($model_file, $execs_H{"cmpress"}, \%opt_HH, \%ofile_info_HH);
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

# validate the CM we are about to use to annotate was actually created for the current reference
$start_secs = outputProgressPrior("Verifying CM database created for current reference $ref_accn", $progress_w, $log_FH, *STDOUT);
validate_cms_built_from_reference($model_file, \%mdl_info_HA, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

##########
# Conclude
##########
$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir, \%ofile_info_HH);

exit 0;


###############
# SUBROUTINES #
################################################################
# naming convention: words separated by underscores
# (e.g. 'get_value_from_array()') for subroutines in *this* file as
# opposed to a Perl module where they're named in camel caps
# (e.g. getValueFromArray()).
#################################################################

#################################################################
# Subroutine : concatenate_individual_cm_files()
# Incept:      EPN, Mon Feb 29 10:53:53 2016
#
# Purpose:     Concatenate individual CM files created and calibrated
#              by a previous call to dnaorg_build.pl into one 
#              CM file.
#
# Arguments: 
#   $model_file:     the full path to the concatenated model file to create
#   $out_root:       output root the individual CM files share
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $ofile_info_HHR: REF to the 2D hash of output file information
# 
# Returns:     void
# 
# Dies: If any of the expected individual CM files do not exist.
#
################################################################# 
sub concatenate_individual_cm_files {
  my $sub_name = "concatenate_individual_cm_files()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($model_file, $out_root, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # concatenate the files into a CM DB and run cmpress on it
  my $cat_cmd = "cat";
  my $out_root_no_dnaorg_annotate = $out_root;
  $out_root_no_dnaorg_annotate =~ s/\.dnaorg_annotate$//;
  for(my $i = 0; $i < $nmdl; $i++) { 
    my $indi_model = $out_root_no_dnaorg_annotate . ".dnaorg_build.ref.$i.cm";
    if(! -s ($indi_model)) { 
      DNAORG_FAIL("ERROR, model database file $model_file does not exist, nor does individual model file $indi_model.\nDid you run 'dnaorg_build.pl $ref_accn' -- it doesn't seem like you have.", 1, $FH_HR);
    }
    $cat_cmd .= " $indi_model ";
  }
  $cat_cmd .= " > $model_file";
  runCommand($cat_cmd, opt_Get("-v", $opt_HHR), $FH_HR);
  addClosedFileToOutputInfo($ofile_info_HHR, "cm", $model_file, "CM file (a concatenation of individual files created by dnaorg_build.pl)");

  # remove the binary index files if they exist, possibly from an earlier cmbuild/cmpress:
  for my $suffix ("i1m", "i1i", "i1f", "i1p") { 
    my $file = $model_file . "." . $suffix;
    if(-e $file) { 
      runCommand("rm $file", opt_Get("-v", $opt_HHR), $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine : press_cm_database()
# Incept:      EPN, Mon Feb 29 14:26:52 2016
#
# Purpose:     Run cmpress on a CM database file.
#
# Arguments: 
#   $model_file:     the full path to the concatenated model file to create
#   $cmpress:        path to cmpress executable
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $ofile_info_HHR: REF to the 2D hash of output file information
# 
# Returns:     void
# 
# Dies: If any of the expected individual CM files do not exist.
#
################################################################# 
sub press_cm_database {
  my $sub_name = "press_cm_database()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($model_file, $cmpress, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # remove the binary index files if they exist, possibly from an earlier cmbuild/cmpress:
  for my $suffix ("i1m", "i1i", "i1f", "i1p") { 
    my $file = $model_file . "." . $suffix;
    if(-e $file) { 
      runCommand("rm $file", opt_Get("-v", $opt_HHR), $FH_HR);
    }
  }

  my $cmpress_cmd = "$cmpress $model_file > /dev/null"; # output is irrelevant
  runCommand($cmpress_cmd, opt_Get("-v", $opt_HHR), $FH_HR);
  addClosedFileToOutputInfo($ofile_info_HHR, "cmi1m", $model_file.".i1m", "index file for the CM, created by cmpress");
  addClosedFileToOutputInfo($ofile_info_HHR, "cmi1i", $model_file.".i1i", "index file for the CM, created by cmpress");
  addClosedFileToOutputInfo($ofile_info_HHR, "cmi1f", $model_file.".i1f", "index file for the CM, created by cmpress");
  addClosedFileToOutputInfo($ofile_info_HHR, "cmi1p", $model_file.".i1p", "index file for the CM, created by cmpress");

  return;
}

#################################################################
# Subroutine : validate_cms_built_from_reference()
# Incept:      EPN, Mon Feb 29 11:21:11 2016
#
# Purpose:     Validate the CMs in a model file were built from
#              the current reference, with information in $mdl_info_HAR->{"cksum"}.
#
# Arguments: 
#  $model_file:      the full path to the concatenated model file to create
#  $mdl_info_HAR:    REF to hash of arrays with information on the models, PRE-FILLED
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:     void
# 
# Dies: If at least one CM was not built from the current reference.
#
################################################################# 
sub validate_cms_built_from_reference { 
  my $sub_name = "validate_cms_built_from_reference()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($model_file, $mdl_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # validate we have a complete model info hash
  validateModelInfoHashIsComplete($mdl_info_HAR, undef, $FH_HR);

  # get the checksum lines from the CM file into a file
  my $cksum_file  = $model_file . ".cksum";
  $cmd = "grep ^CKSUM $model_file | awk '{ print \$2 '} > $cksum_file";
  runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);

  # for each of the checksums in the CM file, make sure they match the checksum from
  # the alignment file used to build that CM, which is in $mdl_info_HAR->{"checksum"}[$i]
  open(CKSUM, $cksum_file) || fileOpenFailure($cksum_file, $!, "reading", $FH_HR);
  my $i = 0;
  while(my $cksum = <CKSUM>) { 
    chomp $cksum;
    if($cksum != $mdl_info_HAR->{"checksum"}[$i]) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, checksum mismatch for CM %d (CM: %d != alignment: %d)", $i+1, $cksum, $mdl_info_HAR->{"checksum"}[$i]), 1, $FH_HR);
    }
    $i++;
  }
  close(CKSUM);

  # clean up this file, unless -d used (in which case we don't remove intermediate files) 
  if(! opt_Get("-d", $opt_HHR)) { 
    unlink $cksum_file;
  }
  else { 
    addClosedFileToOutputInfo($ofile_info_HHR, "cmchecksum", $cksum_file, "Checksum lines from the CM file");
  }
  return;
}
