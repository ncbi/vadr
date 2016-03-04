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
my $esl_ssplit     = "/panfs/pan1/dnaorg/programs/Bio-Easel/scripts/esl-ssplit.pl";
my $esl_epn_translate  = "/home/nawrocke/notebook/15_1118_dnaorg_annotate_genomes_translation/git-esl-epn-translate/esl-epn-translate.pl";

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
opt_Add("-d",           "string",  undef,                    1,    undef, undef,      "directory specified as",                       "specify output directory is <s1> (created with dnaorg_build.pl -d <s>), not <ref accession>", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                   "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--matpept",    "string",  undef,                    1,    undef, undef,      "using pre-specified mat_peptide info",         "read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships", \%opt_HH, \@opt_order_A);
opt_Add("--nomatpept",  "boolean", 0,                        1,    undef,"--matpept", "ignore mat_peptide annotation",                "ignore mat_peptide information in reference annotation", \%opt_HH, \@opt_order_A);
opt_Add("--specstart",  "string",  undef,                    1,    undef, undef,      "using pre-specified alternate start codons",   "read specified alternate start codons per CDS from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",           "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);
opt_Add("--model",      "string",  undef,                    1,    undef, undef,      "use model in file",                            "use model file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--local",      "boolean", 0,                        1,    undef, undef,      "run cmscan locally instead of on farm",        "run cmscan locally instead of on farm", \%opt_HH, \@opt_order_A);
opt_Add("--nseq",       "integer", 5,                        1,    undef,"--local",   "number of sequences for each cmscan farm job", "set number of sequences for each cmscan farm job to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--wait",       "integer", 10,                       1,    undef,"--local",   "allow <n> minutes for cmscan jobs on farm",    "allow <n> minutes for cmscan jobs on farm to finish", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"2"} = "options for skipping stages, primarily useful for debugging";
#     option            type       default               group   requires incompat                 preamble-output            help-output    
opt_Add("--skipedirect", "boolean", 0,                       2,    undef,"--nseq,--local,--wait",   "skip the edirect steps",   "skip the edirect steps, use data from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipscan",   "boolean", 0,                       2,    undef,"--nseq,--local,--wait",   "skip the cmscan step",    "skip the cmscan step, use results from an earlier run of the script", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"3"} = "optional output files";
#       option       type       default                group  requires incompat  preamble-output                          help-output    
opt_Add("--mdlinfo",    "boolean", 0,                        1,    undef, undef, "output internal model information",     "create file with internal model information",   \%opt_HH, \@opt_order_A);
opt_Add("--ftrinfo",    "boolean", 0,                        1,    undef, undef, "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);
opt_Add("--errinfo",    "boolean", 0,                        1,    undef, undef, "output internal error information",     "create file with internal error information", \%opt_HH, \@opt_order_A);

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
                'keep'          => \$GetOptions_H{"--keep"},
                'model=s'       => \$GetOptions_H{"--model"}, 
                'local'         => \$GetOptions_H{"--local"}, 
                'nseq=s'        => \$GetOptions_H{"--nseq"}, 
                'wait=s'        => \$GetOptions_H{"--wait"},
# options for skipping stages
                'skipedirect'   => \$GetOptions_H{"--skipedirect"},
                'skipscan'      => \$GetOptions_H{"--skipscan"},
# optional output files
                'mdlinfo'      => \$GetOptions_H{"--mdlinfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"}, 
                'errinfo'      => \$GetOptions_H{"--errinfo"});

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.1";
my $releasedate   = "Feb 2016";

# make *STDOUT file handle 'hot' so it automatically flushes whenever
# it is printed to
select *STDOUT;
$| = 1;

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

# open optional output files
if(opt_Get("--mdlinfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "mdlinfo", $out_root . ".mdlinfo", "Model information (created due to --mdlinfo)");
}
if(opt_Get("--ftrinfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "ftrinfo", $out_root . ".ftrinfo", "Feature information (created due to --ftrinfo)");
}
if(opt_Get("--errinfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "errinfo", $out_root . ".errinfo", "Error information (created due to --errinfo)");
}

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
$execs_H{"cmscan"}        = $inf_exec_dir . "cmscan";
$execs_H{"cmalign"}       = $inf_exec_dir . "cmalign";
$execs_H{"cmpress"}       = $inf_exec_dir . "cmpress";
$execs_H{"esl-reformat"}  = $esl_exec_dir . "esl-reformat";
$execs_H{"esl_fetch_cds"} = $esl_fetch_cds;
$execs_H{"esl_ssplit"}    = $esl_ssplit;
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################################################################
# Step 1. Gather and process information on reference genome using Edirect.
###########################################################################
my $cmd;             # a command to run with runCommand()
my $progress_w = 70; # the width of the left hand column in our progress output, hard-coded
my $progress_str = (opt_Get("--skipedirect", \%opt_HH)) ? 
    sprintf("Processing information on %d sequences fetched earlier using edirect", scalar(@accn_A)) : 
    sprintf("Gathering information on %d sequences using edirect", scalar(@accn_A));
my $start_secs = outputProgressPrior($progress_str, $progress_w, $log_FH, *STDOUT);

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
my @seq_name_A   = ();          # actual name of all sequences in fasta file, after being fetched, same order as @accn_A
my %mdl_info_HA = ();          # hash of arrays, values are arrays [0..$nmdl-1];
                               # see dnaorg.pm::validateModelInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAndProcessReferenceSequence()
my %ftr_info_HA = ();          # hash of arrays, values are arrays [0..$nftr-1], 
                               # see dnaorg.pm::validateFeatureInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAndProcessReferenceSequence()
my $sqfile = undef;            # pointer to the Bio::Easel::SqFile object we'll open in wrapperFetchAllSequencesAndProcessReferenceSequence()

# Call the wrapper function that does the following:
#   1) fetches the sequences listed in @{$accn_AR} into a fasta file and indexes that fasta file,
#      the reference sequence is $accn_AR->[0].
#   2) determines information for each feature (strand, length, coordinates, product) in the reference sequence
#   3) determines type of each reference sequence feature ('cds-mp', 'cds-notmp', or 'mp')
#   4) fetches the reference sequence feature and populates information on the models and features
wrapperFetchAllSequencesAndProcessReferenceSequence(\@accn_A, \@seq_name_A, \$sqfile, $out_root, \%cds_tbl_HHA,
                                                    ($do_matpept) ? \%mp_tbl_HHA      : undef, 
                                                    ($do_matpept) ? \@cds2pmatpept_AA : undef, 
                                                    ($do_matpept) ? \@cds2amatpept_AA : undef, 
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
  $model_file = $out_root . ".ref.cm";
  $model_file =~ s/\dnaorg\_annotate/dnaorg\_build/; # we name this file dnaorg_build, not dnaorg_annotate
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

#########################################################
# Step 3. Perform homology searches
##########################################################
my $seq_file = $ofile_info_HH{"fullpath"}{"fasta"};
validateFileExistsAndIsNonEmpty($seq_file, undef, $ofile_info_HH{"FH"});

# cmscan output files
my $stdout_file = (opt_Get("--keep", \%opt_HH)) ? $out_root . ".stdout" : "/dev/null"; # only save the stdout if --keep used at cmdline
my $tblout_file = $out_root . ".tblout";

# determine how many jobs we need to run to satisfy <n> per job (where n is from --nseq <n>),
# if more than 1 and --local not used, we split up the sequence file and submit jobs to farm,
# if only 1 or --local used, we just run it locally
my $nfarmjobs = scalar(@accn_A) / opt_Get("--nseq", \%opt_HH); 

if(! opt_Get("--skipscan", \%opt_HH)) { 
  if(($nfarmjobs == 1) || (opt_Get("--local", \%opt_HH))) { 
    # run jobs locally
    $start_secs = outputProgressPrior("Running cmscan locally", $progress_w, $log_FH, *STDOUT);
    run_cmscan($execs_H{"cmscan"}, 1, $model_file, $seq_file, $stdout_file, $tblout_file, \%opt_HH, \%ofile_info_HH); # 1: run locally
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
  else { 
    # we need to split up the sequence file, and submit a separate cmscan job for each
    my @tmp_tblout_file_A = (); # the array of tblout files, we'll remove after we're done, unless --keep
    my @tmp_seq_file_A = ();    # the array of sequence files, we'll remove after we're done, unless --keep
    my @tmp_err_file_A = ();    # the array of error files, we'll remove after we're done, unless --keep
    my $nfasta_created = split_fasta_file($execs_H{"esl_ssplit"}, $seq_file, $nfarmjobs, \%opt_HH, \%ofile_info_HH);
    # we may redefined $nfarmjobs here, split_fasta_file will return the actual number of fasta files created, 
    # which can differ from the requested amount (which is $nfarmjobs) that we pass in
    
    # now submit a job for each
    $start_secs = outputProgressPrior("Submitting $nfasta_created cmscan jobs to the farm", $progress_w, $log_FH, *STDOUT);
    for(my $z = 1; $z <= $nfarmjobs; $z++) { 
      run_cmscan($execs_H{"cmscan"}, 0, $model_file,  # 0: do not run locally
                 $seq_file . "." . $z, 
                 ($stdout_file eq "/dev/null") ? "/dev/null" : $stdout_file . "." . $z,
                 $tblout_file . "." . $z, 
                 \%opt_HH, \%ofile_info_HH);
      push(@tmp_seq_file_A,    $seq_file    . "." . $z);
      push(@tmp_tblout_file_A, $tblout_file . "." . $z);
      push(@tmp_err_file_A,    $tblout_file . "." . $z . ".err"); # this will be the name of the error output file, set in run_cmscan
    }
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
    
    # wait for the jobs to finish
    $start_secs = outputProgressPrior(sprintf("Waiting a maximum of %d minutes for all farm jobs to finish", opt_Get("--wait", \%opt_HH)), 
                                      $progress_w, $log_FH, *STDOUT);
    my $njobs_finished = wait_for_farm_jobs_to_finish(\@tmp_tblout_file_A, "# [ok]", opt_Get("--wait", \%opt_HH));
    if($njobs_finished != $nfasta_created){ 
      DNAORG_FAIL(sprintf("ERROR in main() only $njobs_finished of the $nfasta_created are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", \%opt_HH)), 1, \%{$ofile_info_HH{"FH"}});
    }
    
    # concatenate all the results files into one 
    my $cat_cmd = "cat ";
    foreach my $tmp_file (@tmp_tblout_file_A) { 
      $cat_cmd .= $tmp_file . " ";
    }
    $cat_cmd .= "> $tblout_file";
    runCommand($cat_cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
    
    # remove temporary files, unless --keep
    if(! opt_Get("--keep", \%opt_HH)) { 
      my $rm_cmd = "rm ";
      foreach my $tmp_file (@tmp_seq_file_A, @tmp_tblout_file_A, @tmp_err_file_A) { 
        $rm_cmd .= $tmp_file . " ";
      }
      runCommand($rm_cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
    }
    
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  } # end of 'else' entered if($nfarmjobs > 1 && ! --local)
} # end of if(! opt_Get(--skipscan", \%opt_HH))

#########################################################
# Step 4. Parse homology search results into usable data structuers
##########################################################
#%results_AAH = (); # 1st dim: array, 0..$nmdl-1, one per model
#                   # 2nd dim: array, 0..$naccn-1, one per accessions
#                   # 3rd dim: hash, keys are "p_start", "p_stop", "p_strand", "p_5hangover", "p_3hangover", "p_evalue", "p_fid2ref"
$start_secs = outputProgressPrior("Parsing cmscan results", $progress_w, $log_FH, *STDOUT);
my $naccn = scalar(@accn_A);
# initialize the results AAH
my @results_AAH = ();
for(my $m = 0; $m < $nmdl; $m++) { 
  @{$results_AAH[$m]} = ();
  for(my $a = 0; $a < $naccn; $a++) { 
    %{$results_AAH[$m][$a]} = ();
  }
}

# parse the cmscan results
parse_cmscan_tblout($tblout_file, \%mdl_info_HA, \@seq_name_A, \@accn_A, \%totlen_H, \@results_AAH, $ofile_info_HH{"FH"});

# dump results
# dump_results(\@results_AAH, \%mdl_info_HA, \@seq_name_A, \@accn_A, \%totlen_H, *STDOUT);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################
# Step 5. Fetch predicted hits into fasta files.
##########################################################
my $out_root_and_key = $out_root . ".predicted";
my $results_prefix = "p_"; # we want to fetch the hits with information stored with 3rd dimension keys that begin with "p_" in @results_AAH
$start_secs = outputProgressPrior("Fetching cmscan predicted hits into fasta files", $progress_w, $log_FH, *STDOUT);
fetch_hits_given_results($sqfile, $out_root_and_key, \%mdl_info_HA, \@seq_name_A, \@results_AAH, $results_prefix, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#############################################################################
# Step 6. For features modelled by multiple models (e.g. multi-exon CDS)
#         combine all relevant hits into predicted feature sequences.
#############################################################################
$start_secs = outputProgressPrior("Combining predicted exons into CDS", $progress_w, $log_FH, *STDOUT);
combine_model_hits($out_root_and_key, $results_prefix, \@seq_name_A, \%mdl_info_HA, \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################
# Step 7. For features modelled by multiple other features (e.g. CDS comprised
#         of multiple mature peptides) combine the individual feature sequences into 
#         single sequences.
#########################################################

$start_secs = outputProgressPrior("Combining predicted mature peptides into CDS", $progress_w, $log_FH, *STDOUT);
combine_feature_hits($out_root_and_key, $results_prefix, \@seq_name_A, \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

# output optional output files
if(exists $ofile_info_HH{"FH"}{"mdlinfo"}) { 
  dumpInfoHashOfArrays("Model information (%mdl_info_HA)", 0, \%mdl_info_HA, $ofile_info_HH{"FH"}{"mdlinfo"});
}
if(exists $ofile_info_HH{"FH"}{"ftrinfo"}) { 
  dumpInfoHashOfArrays("Feature information (%ftr_info_HA)", 0, \%ftr_info_HA, $ofile_info_HH{"FH"}{"ftrinfo"});
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
$start_secs = outputProgressPrior("Identifying internal starts/stops in coding sequences", $progress_w, $log_FH, *STDOUT);

#my $accn;                    # name of source accession CDS/mat_peptide sequence came from
my @corr_mft_stop_trc_or_ext_AH = ();  # [0..$i..nmft-1], each element is a hash with keys $key as sequence accessions and values 
                                       # are number of nucleotides that the prediction of the stop coordinate should be corrected
                                       # based on an esl-translate translation of the predicted CDS/mat_peptide sequence, values can be negative
                                       # or positive. If negative: a trc error, if positive: a ext error
my %c_stop_HH  = ();                  # corrected stop positions of hits,  start with a copy of p_stop_HH
my @did_corr_exon_stop_trc_AH = ();   # [0..$i..ref_nexons-1], each element is a hash with keys $key as sequence accessions and values 
                                      # are '1' if this exon's stop position was corrected for a 'trc' error (early stop)
my @did_corr_exon_stop_ext_AH = ();   # [0..$i..ref_nexons-1], each element is a hash with keys $key as sequence accessions and values 
                                      # are '1' if this exon's stop position was corrected for a 'ext' error (late stop (past predicted stop))
 my %corr_fafile_H = ();     # hash of names of fasta files for corrected exon sequences, keys: model name from @mdl_A, value: name of file
 my @corr_mft_fafile_A = (); # array of corrected CDS/mat_peptide sequence files, filled below

my @accn_ext_check_AA     = (); # [0..$c..$ref_nmp_and_cds-1][0..$j..$n] for each feature $c, array of source 
#                                # accessions to check for 'ext' error, b/c no in-frame stop within 
#                                # predicted start..stop
#my @accn_trc_error_AH     = (); # [0..$c..$ref_nmp_and_cds-1][0..$j..$n] for each feature $c, array of source 
#                                # accessions that have 'trc' error, b/c in-frame stop exists between predicted start..stop
#my @accn_ext_error_AH     = (); # [0..$c..$ref_nmp_and_cds-1][0..$j..$n] for each feature $c, array of source 
#                                # accessions that have 'ext' error, b/c no in-frame stop exists between predicted start..stop
#my @accn_nst_error_AH     = (); # 1st dim [0..$c..$ref_nmp_and_cds-1], key: source accn, value: 1 if this 
#                                # feature has no in-frame stop
#                                # 2nd of 2 reasons we shouldn't translate it
#                                # these should be rare, these are a subset of the accessions in
#                                # @accn_ext_check_AA, the ones that had no valid stop in-frame
#                                # for the rest of the sequence (possibly duplicated)
#my @accn_str_error_AH     = (); # 1st dim [0..$c..$ref_nmp_and_cds-1], key: source accn, value: 1 if this 
#                                # feature has an invalid start 
#                                # 1st of 2 reasons we shouldn't translate it
#my @accn_stp_error_AH     = (); # 1st dim [0..$c..$ref_nmp_and_cds-1], key: source accn, value: 1 if this 
#                                # feature's predicted stop is not a valid stop codon
#                                # 1st of 2 reasons we shouldn't translate it

# initialize data structures
my %err_info_HA = (); 
initializeHardCodedErrorInfoHash(\%err_info_HA, $ofile_info_HH{"FH"});

if(exists $ofile_info_HH{"FH"}{"errinfo"}) { 
  dumpInfoHashOfArrays("Error information (%err_info_HA)", 0, \%err_info_HA, $ofile_info_HH{"FH"}{"errinfo"});
}

my @err_instances_AHH = ();
initialize_error_instances_AHH(\@err_instances_AHH, \%err_info_HA, \%ftr_info_HA, $ofile_info_HH{"FH"});

# Translate predicted CDS/mat_peptide sequences using esl-epn-translate to identify 
# in-frame stop codons.
for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
  my $is_matpept = ($ftr_info_HA{"type"}[$ftr_idx] eq "mp")     ? 1 : 0;
  my $is_cds_mp  = ($ftr_info_HA{"type"}[$ftr_idx] eq "cds-mp") ? 1 : 0;
  my $is_cds     = (! $is_matpept);
  my $type_idx   = $ftr_info_HA{"type_idx"}[$ftr_idx];

#  @{$accn_ext_check_AA[$ftr_idx]} = ();
#  %{$corr_mft_stop_trc_or_ext_AH[$ftr_idx]}  = ();
  
  my $ofile_info_key = $ftr_info_HA{"p_hits"}[$ftr_idx];
  my $cur_fafile     = $ofile_info_HH{"fullpath"}{$ofile_info_key};
  
  # use esl-epn-translate.pl to examine the start and stop codons in each feature sequence
  my $esl_epn_translate_outfile = $cur_fafile . ".esl-epn-translate";

  # deal with alternative starts here
  my $altstart_opt = "";
  if($is_cds && opt_IsUsed("--specstart", \%opt_HH)) { 
    $altstart_opt = "-altstart " . get_altstart_list(\@specstart_AA, $type_idx-1, $ofile_info_HH{"FH"}); 
  }
  $cmd = $esl_epn_translate . " $altstart_opt -startstop $cur_fafile > $esl_epn_translate_outfile";
  runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
  
  parse_esl_epn_translate_startstop_outfile($esl_epn_translate_outfile, $ftr_idx, \%ftr_info_HA, \@err_instances_AHH, $ofile_info_HH{"FH"});

  if(! opt_Get("--keep", \%opt_HH)) { 
    removeFileUsingSystemRm($esl_epn_translate_outfile, "dnaorg_annotate.pl:main", \%opt_HH, \%ofile_info_HH);
  }
  else { 
    addClosedFileToOutputInfo(\%ofile_info_HH, "p_esl_epn_translate_outfile.ftr." . $ftr_idx, $esl_epn_translate_outfile, sprintf("esl-epn-translate.pl output file for feature %s", $ftr_info_HA{"out_tiny"}[$ftr_idx]));
  }
}

dumpArrayOfHashesOfHashes("Error instances (%err_instances_AHH)", \@err_instances_AHH, *STDOUT);

exit 0;  

########################################################################
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
  my $out_root_no_dnaorg_annotate = removeScriptNameFromString($out_root); # this removes 'dnaorg_annotate' from $out_root
  for(my $i = 0; $i < $nmdl; $i++) { 
    my $indi_model = $out_root_no_dnaorg_annotate . ".dnaorg_build.ref.$i.cm"; # dnaorg_build created the CMs
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

#################################################################
# Subroutine : run_cmscan()
# Incept:      EPN, Mon Feb 29 15:09:22 2016
#
# Purpose:     Run Infernal's cmscan executable using $model_file
#              as the CM file on sequence file $seq_file.
#
# Arguments: 
#  $cmscan:          path to the cmscan executable file
#  $do_local:        '1' to run locally, '0' to submit job to farm
#  $model_file:      path to the CM file
#  $seq_file:        path to the sequence file
#  $stdout_file:     path to the stdout file to create, can be "/dev/null", or undef 
#  $tblout_file:     path to the cmscan --tblout file to create
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:     void
# 
# Dies: If at least one CM was not built from the current reference.
#
################################################################# 
sub run_cmscan { 
  my $sub_name = "run_cmscan()";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmscan, $do_local, $model_file, $seq_file, $stdout_file, $tblout_file, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  validateFileExistsAndIsNonEmpty($model_file, $sub_name, $FH_HR); 
  validateFileExistsAndIsNonEmpty($seq_file,   $sub_name, $FH_HR);

  # where I got these filter threshold options from: F1, F2, F2b, F3 and F3b from nhmmer, F4, F4b, F5, and F6 from --rfam
  my $opts .= " --noali --cpu 0 --F1 0.02 --F2 0.001 --F2b 0.001 --F3 0.00001 --F3b 0.00001 --F4 0.0002 --F4b 0.0002 --F5 0.0002 --F6 0.0001 --tblout $tblout_file --verbose --nohmmonly ";

  my $cmd = "$cmscan $opts $model_file $seq_file > $stdout_file";

  # run cmscan, either locally or by submitting jobs to the farm
  if($do_local) { 
    # run locally
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);
  }
  else { 
    # submit job to farm and return
    my $jobname = removeDirPath($seq_file);
    my $errfile = $tblout_file . ".err";
    if(-e $errfile) { 
      runCommand("rm $errfile", opt_Get("-v", $opt_HHR), $FH_HR);
    }
    my $farm_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n -l h_rt=288000,h_vmem=8G,mem_free=8G " . "\"" . $cmd . "\" > /dev/null\n";
    runCommand($farm_cmd, opt_Get("-v", $opt_HHR), $FH_HR);
  }

  return;
}

#################################################################
# Subroutine : wait_for_farm_jobs_to_finish()
# Incept:      EPN, Mon Feb 29 16:20:54 2016
#
# Purpose: Wait for jobs on the farm to finish by checking the final
#          line of their output files (in @{$outfile_AR}) to see
#          if the final line is exactly the string
#          $finished_string. We'll wait a maximum of $nmin
#          minutes, then return the number of jobs that have
#          finished. If all jobs finish before $nmin minutes we
#          return at that point.
#
# Arguments: 
#  $outfile_AR:      path to the cmscan executable file
#  $finished_str:    string that indicates a job is finished e.g. "[ok]"
#  $nmin:            number of minutes to wait
# 
# Returns:     Number of jobs (<= scalar(@{$outfile_AR})) that have
#              finished.
# 
# Dies: never.
#
################################################################# 
sub wait_for_farm_jobs_to_finish { 
  my $sub_name = "wait_for_farm_jobs_to_finish()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($outfile_AR, $finished_str, $nmin) = @_;

  my $njobs = scalar(@{$outfile_AR});
  my $nfinished      = 0;   # number of jobs finished
  my $cur_sleep_secs = 15;  # number of seconds to wait between checks, we'll double this until we reach $max_sleep, every $doubling_secs seconds
  my $doubling_secs  = 120; # number of seconds to wait before doublign $cur_sleep
  my $max_sleep_secs = 120; # maximum number of seconds we'll wait between checks
  my $secs_waited    = 0;   # number of total seconds we've waited thus far
  while($secs_waited < (($nmin * 60) + $cur_sleep_secs)) { # we add $cur_sleep so we check one final time before exiting after time limit is reached
    # check to see if jobs are finished, every $cur_sleep seconds
    sleep($cur_sleep_secs);
    $secs_waited += $cur_sleep_secs;
    if($secs_waited >= $doubling_secs) { 
      $cur_sleep_secs *= 2;
      if($cur_sleep_secs > $max_sleep_secs) { # reset to max if we've exceeded it
        $cur_sleep_secs = $max_sleep_secs;
      }
    }

    $nfinished = 0; # important to reset this
    for(my $i = 0; $i < $njobs; $i++) { 
      if(-s $outfile_AR->[$i]) { 
        my $final_line = `tail -n 1 $outfile_AR->[$i]`;
        chomp $final_line;
        if($final_line eq $finished_str) { 
          $nfinished++;
        }
      }
    }
    if($nfinished == $njobs) { 
      # we're done break out of it
      return $nfinished;
    }
  }
  
  return $nfinished;
}

#################################################################
# Subroutine : split_fasta_file()
# Incept:      EPN, Tue Mar  1 09:30:10 2016
#
# Purpose: Split up a fasta file into <n> smaller files by calling
#          the esl-ssplit perl script.
#
# Arguments: 
#  $esl_ssplit:      path to the esl-ssplit.pl script to use
#  $fasta_file:      fasta file to split up
#  $nfiles:          desired number of files to split $fasta_file into
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:    Number of files actually created (can differ from requested
#             amount (which is $nfiles)).
#
# Dies:       if esl-ssplit command fails
#
################################################################# 
sub split_fasta_file { 
  my $sub_name = "split_fasta_file()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_ssplit, $fasta_file, $nfiles, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $outfile = $fasta_file . ".esl-ssplit";
  my $cmd = "$esl_ssplit -v -n $fasta_file $nfiles > $outfile";
  runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);

  # parse output to determine exactly how many files were created:
  my $nfiles_created = countLinesInFile($outfile, $FH_HR);

  if(! opt_Get("--keep", $opt_HHR)) { 
    runCommand("rm $outfile", opt_Get("-v", $opt_HHR), $FH_HR);
  }

  return $nfiles_created;
}

#################################################################
# Subroutine : parse_cmscan_tblout()
# Incept:      EPN, Tue Mar  1 13:56:46 2016
#
# Purpose:    Parse Infernal 1.1 cmscan --tblout output and store
#             results in $results_AAH.
#
# Arguments: 
#  $tblout_file:   tblout file to parse
#  $mdl_info_HAR:  REF to hash of arrays with information on the models, PRE-FILLED
#  $seq_name_AR:    REF to array of sequence names, PRE-FILLED
#  $accn_AR:       REF to array of accessions, same order as @{$seq_name_AR}, PRE-FILLED
#  $totlen_HR:     REF to hash of total lengths of all accessions (keys are values from @{$accn_AR}), PRE-FILLED
#  $results_AAHR:  REF to results AAH, FILLED HERE
#  $FH_HR:         REF to hash of file handles
#
# Returns:    void
#
# Dies:       if we find a hit to a model or sequence that we don't
#             have stored in $mdl_info_HAR or $seq_name_AR
#
################################################################# 
sub parse_cmscan_tblout { 
  my $sub_name = "parse_cmscan_tblout()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($tblout_file, $mdl_info_HAR, $seq_name_AR, $accn_AR, $totlen_HR, $results_AAHR, $FH_HR) = @_;
  
  # make an 'order hash' for the model names and sequence names,
  my %mdlname_index_H = (); # mdlname_index_H{$model_name} = <n>, means that $model_name is the <n>th model in the @{$mdl_info_HAR{*}} arrays
  my %seqname_index_H = (); # seqname_order_H{$seq_name} = <n>, means that $seq_name is the <n>th sequence name in the @{$seq_name_AR}} array
  getIndexHashForArray($mdl_info_HAR->{"cmname"}, \%mdlname_index_H, $FH_HR);
  getIndexHashForArray($seq_name_AR, \%seqname_index_H, $FH_HR);

  open(IN, $tblout_file) || fileOpenFailure($tblout_file, $sub_name, $!, "reading", $FH_HR);

  my $did_field_check = 0; # set to '1' below after we check the fields of the file
  my $line_ctr = 0;
  while(my $line = <IN>) { 
    $line_ctr++;
    if(($line =~ m/^\#/) && (! $did_field_check)) { 
      # sanity check, make sure the fields are what we expect
      if($line !~ m/#target name\s+accession\s+query name\s+accession\s+mdl\s+mdl\s+from\s+mdl to\s+seq from\s+seq to\s+strand\s+trunc\s+pass\s+gc\s+bias\s+score\s+E-value inc description of target/) { 
        DNAORG_FAIL("ERROR in $sub_name, unexpected field names in $tblout_file\n$line\n", 1, $FH_HR);
      }
      $did_field_check = 1;
    }
    elsif($line !~ m/^\#/) { 
      chomp $line;
      # example line:
      #Maize-streak_r23.NC_001346.ref.mft.4        -         NC_001346:genome-duplicated:NC_001346:1:2689:+:NC_001346:1:2689:+: -          cm        1      819     2527     1709      -    no    1 0.44   0.2  892.0         0 !   -
      my @elA = split(/\s+/, $line);
      my ($mdlname, $seqname, $mod, $mdlfrom, $mdlto, $from, $to, $strand, $score, $evalue) = 
          ($elA[0], $elA[2], $elA[4], $elA[5], $elA[6], $elA[7], $elA[8], $elA[9], $elA[14], $elA[15]);

      if(! exists $mdlname_index_H{$mdlname}) { 
        DNAORG_FAIL("ERROR in $sub_name, do not have information for model $mdlname read in $tblout_file on line $line_ctr", 1, $FH_HR);
      }
      if(! exists $seqname_index_H{$seqname}) { 
        DNAORG_FAIL("ERROR in $sub_name, do not have information for sequence $seqname read in $tblout_file on line $line_ctr", 1, $FH_HR);
      }

      my $mdlidx = $mdlname_index_H{$mdlname}; # model    index for the hit in results_AAH (1st dim of results_AAH)
      my $seqidx = $seqname_index_H{$seqname}; # sequence index for the hit in results_AAH (2nd dim of results_AAH)
      my $mdllen = $mdl_info_HAR->{"length"}[$mdlidx]; # model length, used to determine how far hit is from boundary of the model
      if(! exists $totlen_HR->{$accn_AR->[$seqidx]}) { 
        DNAORG_FAIL("ERROR in $sub_name, do not have length information for sequence $seqname, accession $accn_AR->[$seqidx]", 1, $FH_HR);
      }
      my $seqlen = $totlen_HR->{$accn_AR->[$seqidx]}; # sequence length, used to exclude storing of hits that start and stop after $seqlen, 
                                                      # which can occur in circular genomes, where we've duplicated the sequence

      store_hit($results_AAHR, $mdlidx, $seqidx, $mdllen, $seqlen, $mdlfrom, $mdlto, $from, $to, $strand, $evalue, $FH_HR);
    }
  }
  close(IN);
  
  return;
}

#################################################################
# Subroutine : store_hit()
# Incept:      EPN, Tue Mar  1 14:33:38 2016
#
# Purpose:    Helper function for parse_cmscan_tblout().
#             Given info on a hit and a ref to the results AAH, 
#             store info on it. 
#
# Arguments: 
#  $results_AAHR:  REF to results AAH, FILLED HERE
#  $mdlidx:        model index, 1st dim index in results_AAH to store in
#  $seqidx:        sequence index, 2nd dim index in results_AAH to store in
#  $mdllen:        model length
#  $seqlen:        sequence length (before duplicating, if relevant)
#  $mdlfrom:       start position of hit
#  $mdlto:         stop position of hit
#  $seqfrom:       start position of hit
#  $seqto:         stop position of hit
#  $strand:        strand of hit
#  $evalue:        E-value of hit
#  $FH_HR:         REF to hash of file handles
#
# Returns:    void
#
# Dies:       never
#
################################################################# 
sub store_hit { 
  my $sub_name = "store_hit()";
  my $nargs_exp = 12;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($results_AAHR, $mdlidx, $seqidx, $mdllen, $seqlen, $mdlfrom, $mdlto, $seqfrom, $seqto, $strand, $evalue, $FH_HR) = @_;

  # only consider hits where either the start or end are less than the total length
  # of the genome. Since we sometimes duplicate all genomes, this gives a simple 
  # rule for deciding which of duplicate hits we'll store 
  if(($seqfrom <= $seqlen) || ($seqto <= $seqlen)) { 

# Code block below is how we used to modify start/stop if hit spanned
# stop..start boundary, now we just store it as it is, and then modify
# the start and/or stop when we output its coordinates. This removes a
# nasty off-by-one with negative indices that kept recurring in
# various situations.
# 
# This is kept here for reference
#
#    # deal with case where one but not both of from to is > L:
#    if($seqfrom > $L || $seqto > $L) { 
#      $seqfrom -= $L; 
#      $seqto   -= $L; 
#      if($seqfrom < 0)  { $seqfrom--; }
#      if($seqto   < 0)  { $seqto--; }
#    }
    
    # we only store the first hit we see, this is safe because we know 
    # that this will be the lowest E-value
    if(%{$results_AAHR->[$mdlidx][$seqidx]}) { 
      # a hit for this model:seq pair already exists, make sure it has a lower E-value than the current one
      if($results_AAHR->[$mdlidx][$seqidx]{"p_evalue"} > $evalue) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, already have hit stored for model index $mdlidx seq index $seqidx with higher evalue (%g > %g), this implies hits were not sorted by E-value...", $results_AAHR->[$mdlidx][$seqidx]{"evalue"}, $evalue), 1, $FH_HR); 
      }
    }
    else { 
      # no hit yet exists, make one
      %{$results_AAHR->[$mdlidx][$seqidx]} = ();
      $results_AAHR->[$mdlidx][$seqidx]{"p_start"}     = $seqfrom;
      $results_AAHR->[$mdlidx][$seqidx]{"p_stop"}      = $seqto;
      $results_AAHR->[$mdlidx][$seqidx]{"p_strand"}    = $strand;
      $results_AAHR->[$mdlidx][$seqidx]{"p_5hangover"} = ($mdlfrom - 1);
      $results_AAHR->[$mdlidx][$seqidx]{"p_3hangover"} = ($mdllen - $mdlto);
      $results_AAHR->[$mdlidx][$seqidx]{"p_evalue"}    = $evalue;
    }
  }

  return;
}

#################################################################
# Subroutine : dump_results()
# Incept:      EPN, Tue Mar  1 14:54:27 2016
#
# Purpose:    Dump results data structure to $FH. Probably only
#             useful for debugging.
#
# Arguments: 
#  $results_AAHR:  REF to results AAH, PRE-FILLED
#  $mdl_info_HAR:  REF to hash of arrays with information on the models, PRE-FILLED
#  $seq_name_AR:    REF to array of sequence names, PRE-FILLED
#  $accn_AR:       REF to array of accessions, same order as @{$seq_name_AR}, PRE-FILLED
#  $totlen_HR:     REF to hash of total lengths of all accessions (keys are values from @{$accn_AR}), PRE-FILLED
#  $FH:            file handle to output to
#
# Returns:    void
#
# Dies:       never
#
################################################################# 
sub dump_results {
  my $sub_name = "dump_results()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($results_AAHR, $mdl_info_HAR, $seq_name_AR, $accn_AR, $totlen_HR, $FH) = @_;

  my $nmdl = scalar(@{$results_AAHR});
  my $nseq = scalar(@{$seq_name_AR}); 
  
  for(my $m = 0; $m < $nmdl; $m++) { 
    for(my $s = 0; $s < $nseq; $s++) { 
      printf $FH ("model($m): %20s  seq($s): %20s  accn: %10s  len: %10d  ", 
                  $mdl_info_HAR->{"cmname"}[$m],
                  $seq_name_AR->[$s],
                  $accn_AR->[$s],
                  $totlen_HR->{$accn_AR->[$s]});
      foreach my $key (sort keys %{$results_AAHR->[$m][$s]}) { 
        printf $FH ("%s: %s ", $key, $results_AAHR->[$m][$s]{$key});
      }
      printf $FH ("\n")
    }
  }
  return;
}

#################################################################
# Subroutine:  fetch_hits_given_results()
# Incept:      EPN, Wed Mar  2 15:25:55 2016
#
# Purpose:    Given the results data structure, fetch all
#             hits to fasta files.
#             As an extra step, fetch any 'appended' sequence
#             for models in which $mdl_info_HAR->{"append_num"}[$mdl_idx]
#             is non-zero. This will go to a separate file,
#             and we'll append it to sequence hits from other 
#             models in combine_feature_hits.
# Arguments: 
#  $sqfile:            REF to Bio::Easel::SqFile object, open sequence file containing
#                      usually $out_root . ".predicted", or $out_root . ".corrected"
#  $out_root_and_key:  output root for the files we'll create here, 
#                      usually $out_root . ".predicted", or $out_root . ".corrected"
#  $mdl_info_HAR:      REF to hash of arrays with information on the models, PRE-FILLED
#  $seq_name_AR:       REF to array of sequence names, PRE-FILLED
#  $results_AAHR:      REF to results AAH, PRE-FILLED
#  $results_prefix:    prefix to results 3rd dim keys to use, e.g. "p_" to use "p_start" and "p_stop"
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#################################################################
sub fetch_hits_given_results { 
  my $sub_name = "fetch_hits_given_results";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $out_root_and_key, $mdl_info_HAR, $seq_name_AR, $results_AAHR, $results_prefix, $opt_HHR, $ofile_info_HHR) = @_;
    
  my $nmdl = scalar(@{$results_AAHR});
  my $nseq = scalar(@{$seq_name_AR}); 

  my $start_key  = $results_prefix . "start";  # key for 'start'  in 3rd dim of @{$results_AAHR}
  my $stop_key   = $results_prefix . "stop";   # key for 'stop'   in 3rd dim of @{$results_AAHR}
  my $strand_key = $results_prefix . "strand"; # key for 'strand' in 3rd dim of @{$results_AAHR}
  
  my $mdl_info_key        = $results_prefix . "hits";
  my $mdl_info_append_key = $results_prefix . "hits.append";

  for(my $m = 0; $m < $nmdl; $m++) { 
    my @fetch_AA        = ();
    my @fetch_append_AA = ();
    my $mdl_name = $mdl_info_HAR->{"cmname"}[$m];
    my $nseq2fetch        = 0; # number of sequences we'll fetch
    my $nseq2fetch_append = 0; # number of sequences we'll fetch for append file
    my $append_num = ($mdl_info_HAR->{"append_num"}[$m]);

    for(my $s = 0; $s < $nseq; $s++) { 
      my $seq_name = $seq_name_AR->[$s];
      if($s == 0 && (! %{$results_AAHR->[$m][$s]})) { 
        DNAORG_FAIL("ERROR in $sub_name(), no hit from model $mdl_name to the reference sequence $seq_name", 1, $ofile_info_HHR->{"FH"}); 
      }
      if(%{$results_AAHR->[$m][$s]}) { # hit exists
        if(! exists $results_AAHR->[$m][$s]{$start_key}) { 
          DNAORG_FAIL("ERROR in $sub_name(), no start value for hit of model $mdl_name to the sequence $seq_name", 1, $ofile_info_HHR->{"FH"}); 
        }
        if(! exists $results_AAHR->[$m][$s]{$stop_key}) { 
          DNAORG_FAIL("ERROR in $sub_name(), no stop value for hit of model $mdl_name to the sequence $seq_name", 1, $ofile_info_HHR->{"FH"}); 
        }
        my $start  = $results_AAHR->[$m][$s]{$start_key};
        my $stop   = $results_AAHR->[$m][$s]{$stop_key};
        my $strand = $results_AAHR->[$m][$s]{$strand_key};

        my $new_name = $seq_name . "/" . $start . "-" . $stop;
        push(@fetch_AA, [$new_name, $start, $stop, $seq_name]);

        $nseq2fetch++;

        if($append_num > 0) { 
          my $append_start;
          my $append_stop;
          if($strand eq "+") { 
            $append_start = $stop + 1;
            $append_stop  = $stop + $append_num;
          }
          else { 
            $append_start = $stop - 1;
            $append_stop  = $stop + $append_num;
          }
          my $append_new_name = $seq_name . "/" . $append_start . "-" . $append_stop;
          printf("added $append_new_name $append_start $append_stop $seq_name to fetch_append_AA\n");
          push(@fetch_append_AA, [$append_new_name, $append_start, $append_stop, $seq_name]);
        
          $nseq2fetch_append++;
        }
      }
    }

    if($nseq2fetch > 0) { 
      my $out_fafile = $out_root_and_key . ".". $mdl_info_HAR->{"filename_root"}[$m] . ".fa";
      $sqfile->fetch_subseqs(\@fetch_AA, undef, $out_fafile);
      # save information on this to the output file info hash
      my $ofile_info_key = $results_prefix . "hits" . ".model" . $m;
      addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_key, $out_fafile, "fasta file with predicted hits for model " . $mdl_info_HAR->{"out_tiny"}[$m]);
      # now save this file in the mdl_info_HAR
      $mdl_info_HAR->{$mdl_info_key}[$m] = $ofile_info_key;

      if($nseq2fetch_append > 0) { 
        my $out_fafile = $out_root_and_key . ".". $mdl_info_HAR->{"filename_root"}[$m] . ".append.fa";
        $sqfile->fetch_subseqs(\@fetch_append_AA, undef, $out_fafile);
        # save information on this to the output file info hash
        my $ofile_info_key = $results_prefix . "hits" . ".model" . $m . ".append";
        addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_key, $out_fafile, "fasta file with predicted appended hits for model " . $mdl_info_HAR->{"out_tiny"}[$m]);
        # now save this file in the mdl_info_HAR
        $mdl_info_HAR->{$mdl_info_append_key}[$m] = $ofile_info_key;
      }
      else { 
        $mdl_info_HAR->{$mdl_info_append_key}[$m] = ""; # no file to append for this model either
      }
    }
    else { 
      $mdl_info_HAR->{$mdl_info_key}[$m]        = ""; # no file for this model
      $mdl_info_HAR->{$mdl_info_append_key}[$m] = ""; # no file to append for this model either
    }
  }

  return;
}

#################################################################
# Subroutine:  combine_model_hits()
# Incept:      EPN, Thu Mar  3 11:39:17 2016
#
# Purpose:    For all features annotated by models 
#             ($ftr_info_HAR->{"annot_type"}[*] = "model")
#             assign (for single model features) or create
#             (for multiple model features) the hit fasta file
#             for each. Calls 'combine_sequences()' which does
#             much of the work.
#
# Arguments: 
#  $out_root_and_key:  output root for the files we'll create here, 
#                      usually $out_root . ".predicted", or $out_root . ".corrected"
#  $results_prefix:    prefix used for naming mdl fasta files in fetch_hits(), e.g. "p_" for 'predicted'
#  $seq_name_AR:       REF to array of sequence names, PRE-FILLED
#  $mdl_info_HAR:      REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:      REF to hash of arrays with information on the features, ADDED TO HERE
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If we have a problem reading the fasta files
#
################################################################# 
sub combine_model_hits { 
  my $sub_name = "combine_model_hits";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_root_and_key, $results_prefix, $seq_name_AR, $mdl_info_HAR, $ftr_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nftr: number of features
  my $nmdl = validateModelInfoHashIsComplete  ($mdl_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nmdl: number of homology models

  my $ftr_info_key        = $results_prefix . "hits";
  my $ftr_info_append_key = $results_prefix . "hits.append";
  my $mdl_info_key        = $results_prefix . "hits";
  my $mdl_info_append_key = $results_prefix . "hits.append";

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { # we only do this for features annotated by models
      my $mdl_idx        = $ftr_info_HAR->{"first_mdl"}[$ftr_idx];
      my $ofile_info_key = $mdl_info_HAR->{$mdl_info_key}[$mdl_idx];
      my $mdl_hit_fafile = $ofile_info_HH{"fullpath"}{$ofile_info_key};
      validateFileExistsAndIsNonEmpty($mdl_hit_fafile, $sub_name, $ofile_info_HHR->{"FH"});

      #######################################
      # single model (e.g. exon) features
      #######################################
      if($ftr_info_HAR->{"nmodels"}[$ftr_idx] == 1) { 
        # a single model (e.g. exon) gene, we should already have the sequence from fetch_hits
        $ftr_info_HAR->{$ftr_info_key}[$ftr_idx] = $ofile_info_key;
      }
      #######################################
      # multi model (e.g. exon) features
      #######################################
      else { 
        # more than one model's hit files need to be combined to make this feature 
        my $ftr_hit_fafile = $out_root_and_key . ".". $ftr_info_HAR->{"filename_root"}[$ftr_idx] . ".fa";
        my @tmp_hit_fafile_A = ($mdl_hit_fafile);
        for($mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx] + 1; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
          $mdl_info_key   = $results_prefix . "hits";
          $ofile_info_key = $mdl_info_HAR->{$mdl_info_key}[$mdl_idx];
          $mdl_hit_fafile = $ofile_info_HH{"fullpath"}{$ofile_info_key};
          validateFileExistsAndIsNonEmpty($mdl_hit_fafile, $sub_name, $ofile_info_HHR->{"FH"});
          push(@tmp_hit_fafile_A, $mdl_hit_fafile);
        }
        # combine the sequences into 1 file
        combine_sequences(\@tmp_hit_fafile_A, $seq_name_AR, $ftr_hit_fafile, $ofile_info_HHR->{"FH"});
        my $ofile_info_key = $results_prefix . "hits" . ".ftr" . $ftr_idx;
        addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_key, $ftr_hit_fafile, "fasta file with predicted hits for feature " . $ftr_info_HAR->{"out_tiny"}[$ftr_idx] . " from " . $ftr_info_HAR->{"nmodels"}[$ftr_idx] . " combined model predictions");
        # now save this file in the ftr_info_HAR
        $ftr_info_HAR->{$ftr_info_key}[$ftr_idx] = $ofile_info_key;
      }

      # check if there's a file to append
      my $final_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$ftr_idx];
      if($mdl_info_HAR->{$mdl_info_append_key}[$final_mdl_idx] ne "") { 
        # yes, there is
        my $ofile_info_append_key = $mdl_info_HAR->{$mdl_info_append_key}[$final_mdl_idx];
        my $mdl_hit_append_fafile = $ofile_info_HH{"fullpath"}{$ofile_info_append_key};
        validateFileExistsAndIsNonEmpty($mdl_hit_append_fafile, $sub_name, $ofile_info_HHR->{"FH"});
        $ftr_info_HAR->{$ftr_info_append_key}[$ftr_idx] = $ofile_info_append_key;
      }
      else { # no file to append, set to ""
        $ftr_info_HAR->{$ftr_info_append_key}[$ftr_idx] = "";
      }
    } # end of 'if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model")'
    else { # $ftr_info_HAR->{"annot_type"}[$ftr_idx] ne "model"
      # these never have anything appended
      $ftr_info_HAR->{$ftr_info_append_key}[$ftr_idx] = "";
    }
  }
  
  return;
}

#################################################################
# Subroutine:  combine_feature_hits()
# Incept:      EPN, Thu Mar  3 12:15:55 2016
#
# Purpose:    For all features that are annotated by combining
#             multiple other features ($ftr_info_HAR->{"annot_type"}[*] 
#             = "multifeature", e.g. CDS made up of mature peptides)
#             create the feature fasta file for each. Calls 
#             'combine_sequences()' which does much of the work.
#
# Arguments: 
#  $out_root_and_key:  output root for the files we'll create here, 
#                      usually $out_root . ".predicted", or $out_root . ".corrected"
#  $results_prefix:    prefix used for naming mdl fasta files in fetch_hits(), e.g. "p_" for 'predicted'
#  $seq_name_AR:       REF to array of sequence names, PRE-FILLED
#  $ftr_info_HAR:      REF to hash of arrays with information on the features, ADDED TO HERE
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If we have a problem reading the fasta files
#
################################################################# 
sub combine_feature_hits { 
  my $sub_name = "combine_feature_hits";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_root_and_key, $results_prefix, $seq_name_AR, $ftr_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nftr: number of features
  my $ftr_info_key        = $results_prefix . "hits";
  my $ftr_info_append_key = $results_prefix . "hits.append";

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") { # we only do this for features annotated by models
      # get the array of primary children feature indices for this feature
      my @primary_children_idx_A = (); # feature indices of the primary children of this feature
      getPrimaryChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, \@primary_children_idx_A, $ofile_info_HHR->{"FH"});
      my @tmp_hit_fafile_A = ();
      my $combined_ftr_hit_fafile = $out_root_and_key . ".". $ftr_info_HAR->{"filename_root"}[$ftr_idx] . ".fa";
      foreach $ftr_idx (@primary_children_idx_A) { 
        my $ofile_info_key = $ftr_info_HAR->{$ftr_info_key}[$ftr_idx];
        my $ftr_hit_fafile = $ofile_info_HH{"fullpath"}{$ofile_info_key};
        validateFileExistsAndIsNonEmpty($ftr_hit_fafile, $sub_name, $ofile_info_HHR->{"FH"});
        push(@tmp_hit_fafile_A, $ftr_hit_fafile);
        # check if this feature has a mandatory file to append
        if($ftr_info_HAR->{$ftr_info_append_key}[$ftr_idx] ne "") { 
          # it does, append it
          my $ofile_info_append_key = $ftr_info_HAR->{$ftr_info_append_key}[$ftr_idx];
          my $ftr_hit_append_fafile = $ofile_info_HH{"fullpath"}{$ofile_info_append_key};
          validateFileExistsAndIsNonEmpty($ftr_hit_append_fafile, $sub_name, $ofile_info_HHR->{"FH"});
          push(@tmp_hit_fafile_A, $ftr_hit_append_fafile);
        }
      }
      # combine the sequences into 1 file
      combine_sequences(\@tmp_hit_fafile_A, $seq_name_AR, $combined_ftr_hit_fafile, $ofile_info_HHR->{"FH"});
      my $ofile_info_key = $results_prefix . "hits" . ".ftr" . $ftr_idx;
      addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_key, $combined_ftr_hit_fafile, "fasta file with predicted hits for feature " . $ftr_info_HAR->{"out_tiny"}[$ftr_idx] . " from " . scalar(@primary_children_idx_A) . " combined feature predictions");

      # now save this file in the ftr_info_HAR
      $ftr_info_HAR->{$ftr_info_key}[$ftr_idx] = $ofile_info_key;
    }
  }
  return;
}

#################################################################
# Subroutine:  combine_sequences()
# Incept:      EPN, Wed Mar  2 16:11:40 2016
#
# Purpose:    Helper function for combine_multiple_hits(). Given an array of
#             fasta files, each with a different subsequence from the
#             same parent sequences, create a single new fasta file
#             that has the subsequences concatenated together.  An
#             example is stitching together exons into a CDS.  Uses
#             BioEasel's sqfile module.
#
# Arguments: 
#  $indi_file_AR: REF to array of fasta files to combine
#  $seq_name_AR:  REF to array with order of sequence names
#  $multi_file:   name of multi file to create
#  $FH_HR:        REF to hash of file handles
#
# Returns:    void
#
# Dies:       If we have a problem reading the fasta files
#
#################################################################
sub combine_sequences {
  my $sub_name = "combine_sequences";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($indi_file_AR, $seq_name_AR, $multi_file, $FH_HR) = @_;
  
  my @sqfile_A  = (); # array of open Bio::Easel::SqFile objects, one per indi_file_AR element
  my @sqname_AA = (); # 2D array of sequence names, read from the indi files
  my $nfiles = scalar(@{$indi_file_AR});
  my %rawname_map_HA   = (); # key: sequence name minus coords, value: array [0..$f..$nfiles-1], value: ssi index of accession in file $f
  my %rawname_ct_H     = (); # key: sequence name minus coords, value: number of files accession exists in
  my %rawname_coords_H = (); # key: sequence name minus coords, value: concatenated set of coordinates for this accession
  for(my $f = 0; $f < $nfiles; $f++) { 
    if(! -s $indi_file_AR->[$f]) { 
      DNAORG_FAIL("ERROR in $sub_name, file $indi_file_AR->[$f] does not exist, we need to use it to combine exons into a CDS file", 1, $FH_HR);
    }
    $sqfile_A[$f] = Bio::Easel::SqFile->new({ fileLocation => $indi_file_AR->[$f] });
    # get names all of sequences in each file
    for(my $i = 0; $i < $sqfile_A[$f]->nseq_ssi; $i++) { 
      $sqname_AA[$f][$i] = $sqfile_A[$f]->fetch_seq_name_given_ssi_number($i);
      my ($acc, $coords) = split("/", $sqname_AA[$f][$i]);
      if(! defined $coords || $coords !~ m/\-/) { 
        DNAORG_FAIL("ERROR in $sub_name, unable to parse sequence name $sqname_AA[$f][$i] into accession and coordinates", 1, $FH_HR);
      }
      $rawname_map_HA{$acc}[$f] = $i;
      if(! exists $rawname_ct_H{$acc}) { 
        $rawname_ct_H{$acc} = 1;
        $rawname_coords_H{$acc} = $coords;
      }
      else { 
        $rawname_ct_H{$acc}++;
        $rawname_coords_H{$acc} .= "," . $coords;
      }
    }
  }
  
  # now for each accession that exists in all files, fetch all exons for that accession
  # into a new sequence
  open(OUT, ">", $multi_file) || die "ERROR unable to open $multi_file for writing";
  foreach my $seq_name (@{$seq_name_AR}) { 
    if(exists $rawname_ct_H{$seq_name} && $rawname_ct_H{$seq_name} == $nfiles) { 
      print OUT ">" . $seq_name . "/" . $rawname_coords_H{$seq_name} . "\n";
      for(my $f = 0; $f < $nfiles; $f++) { 
        my $sqname = $sqname_AA[$f][($rawname_map_HA{$seq_name}[$f])];
        my $sqonly = $sqfile_A[$f]->fetch_seq_to_fasta_string($sqname);
        $sqonly =~ s/^\>.+\n//;
        print OUT $sqonly;
      }
    }
  }
  close(OUT);

  # clean up: remove all 'ssi' files we just created
  for(my $f = 0; $f < $nfiles; $f++) { 
    if(-e $indi_file_AR->[$f] . ".ssi") { 
      unlink $indi_file_AR->[$f] . ".ssi";
    }
  }

  return;
}

#################################################################
# Subroutine:  initialize_error_instances_AHH()
# Incept:      EPN, Fri Mar  4 12:26:42 2016
#
# Purpose:    Initialize the error instances array of arrays of 
#             2 dimensional hashes. The array is [0..$f..$nftr-1] 
#             where $f is a feature index in %{ftr_info_HA}. 
#             The key in the 1st hash dimension is an error code, 
#             the key in the 2nd hash dimension is a sequence 
#             name (from $accn_AR).
#
# Arguments: 
#  $err_instances_AHHR: REF to the array of 2D hashes, initialized here
#  $err_info_HAR:       REF to the error info hash of arrays, PRE-FILLED
#  $ftr_info_HAR:       REF to the feature info hash of arrays, PRE-FILLED
#  $FH_HR:              REF to hash of file handles
#
# Returns:    void
#
# Dies:       If err_info_HAR is not complete
#
#################################################################
sub initialize_error_instances_AHH { 
  my $sub_name = "initialize_error_instances_AHH";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_instances_AHHR, $err_info_HAR, $ftr_info_HAR, $FH_HR) = @_;
  
  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $FH_HR); 
  my $nerr = validateErrorInfoHashIsComplete($err_info_HAR, undef, $FH_HR); 

  @{$err_instances_AHHR} = ();
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    %{$err_instances_AHHR->[$ftr_idx]} = (); 
    for(my $err_idx = 0; $err_idx < $nerr; $err_idx++) { 
      %{$err_instances_AHHR->[$ftr_idx]{$err_info_HAR->{"code"}[$err_idx]}} = ();
    }
  }

  return;
}

#################################################################
# Subroutine:  parse_esl_epn_translate_startstop_outfile()
# Incept:      EPN, Fri Mar  4 13:56:56 2016
#
# Purpose:    Parse an output file from esl-epn-translate.pl run
#             with the --startstop option and store the relevant 
#             information we derive from it in @{$err_instances_AHHR}. 
#
# Arguments: 
#  $translate_outfile:  path to the file to parse
#  $ftr_idx:            feature index the translate output is for
#  $ftr_info_HAR:       REF to the feature info hash of arrays, PRE-FILLED
#  $err_instances_AHHR: REF to the error instances AAH, PRE-FILLED
#  $FH_HR:              REF to hash of file handles
#
# Returns:    void
#
# Dies:       If we have trouble parsing $translate_outfile
#
#################################################################
sub parse_esl_epn_translate_startstop_outfile { 
  my $sub_name = "parse_esl_epn_translate_startstop_outfile";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($translate_outfile, $ftr_idx, $ftr_info_HAR, $err_instances_AHHR, $FH_HR) = @_;
  
  # is this a mature peptide?
  my $is_matpept = ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp") ? 1 : 0;

  open(IN, $translate_outfile) || fileOpenFailure($translate_outfile, "main", $!, "reading", $FH_HR);

  while(my $line = <IN>) { 
    # example line
    #HQ693465/1-306 1 1 304
    if($line =~ /^(\S+)\/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/) { 
      my ($seq_name, $coords, $start_is_valid, $stop_is_valid, $first_stop_posn1) = ($1, $2, $3, $4, $5);
      
      # determine if we have an early stop
      my $cds_len            = dashCoordsStringCommaDelimitedToLength($coords, $sub_name, $FH_HR);
      my $final_codon_posn1  = $cds_len - 2; # final codon position 1 
      my $early_inframe_stop;
      if($first_stop_posn1 == $final_codon_posn1) { # first stop is the final codon
        $early_inframe_stop = 0;
      }
      elsif($first_stop_posn1 == 0) { # esl-epn-translate didn't find any in-frame stop
        $early_inframe_stop = 0;
      }
      else { # there is an early stop
        $early_inframe_stop = 1; 
      }

      # We now have all of the relevant data on the current
      # CDS/mat_peptide sequence and we need to use to determine
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
      # '?' after error code means that error is possible, we have to check for it later           
      # 'any' means that any value is possible, outcome is unaffected by value
      #
      if(! $is_matpept) { 
        if(! $start_is_valid) { # possibility 1 (P1)
          $err_instances_AHHR->[$ftr_idx]{"str"}{$seq_name} = "yes";
          # printf("in $sub_name, feature index $ftr_idx, seq $seq_name $c possibility 1 (str)\n");
        }
        else { 
          # $start_is_valid is 1
          if(! $stop_is_valid) { 
            if(! $early_inframe_stop) { 
              # possibility 2 (P2): stp error, need to check for ext error later
              $err_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name} = "yes";
              $err_instances_AHHR->[$ftr_idx]{"ext"}{$seq_name} = "maybe";
              #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 2 (stp, maybe ext)\n");
            }
            else { # $early_inframe_stop is 1
              # possibility 3 (P3): stp and trc error
              $err_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name} = "yes";
              $err_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name} = "yes" . "." . -1 * ($final_codon_posn1 - $first_stop_posn1);
              if($final_codon_posn1 == $first_stop_posn1) { 
                DNAORG_FAIL(sprintf("ERROR in $sub_name, trying to correct a stop codon by 0 nt for feature %s sequence $seq_name", $ftr_info_HAR->{"out_tiny"}[$ftr_idx]), 1, $FH_HR); 
              }
              #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 3 (trc and stp)\n");
            }
          } # end of 'if(! $stop_is_valid)'
          else { # $stop_is_valid is 1
            if(! $early_inframe_stop) { 
              ; 
              # possibility 4 (P4): no errors, do nothing
              #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 4 (no errors)\n");
            }
            else { 
              # possibility 5 (P5): trc error
              $err_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name} = "yes" . "." . -1 * ($final_codon_posn1 - $first_stop_posn1);
              if($final_codon_posn1 == $first_stop_posn1) { 
                DNAORG_FAIL(sprintf("ERROR in $sub_name, trying to correct a stop codon by 0 nt for feature %s sequence $seq_name", $ftr_info_HAR->{"out_tiny"}[$ftr_idx]), 1, $FH_HR); 
              }
              #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 5 (trc)\n");
            }
          }              
        }
      } # end of 'if(! $is_matpept)'
      else { # $is_matpept is 1 
        if(! $early_inframe_stop) { 
          ; 
          # possibility 6 (P6): maybe ntr error later, but can't check for it now, do nothing;
          #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 6 (no error)\n");
        }
        else { # $early_inframe_stop is '1'
          # possibility 7 (P7): trc error, maybe ntr error later, but can't check for it now
          #printf("HEYA $accn $ftr_idx possibility 7 (trc)\n");
          $err_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name} = "yes" . "." . -1 * ($final_codon_posn1 - $first_stop_posn1);
          if($final_codon_posn1 == $first_stop_posn1) { 
            DNAORG_FAIL(sprintf("ERROR in $sub_name, trying to correct a stop codon by 0 nt for feature %s sequence $seq_name", $ftr_info_HAR->{"out_tiny"}[$ftr_idx]), 1, $FH_HR); 
          }
          #printf("HEYAA set corr_mft_stop_trc_or_ext_AH[$ftr_idx]{$accn} to $corr_mft_stop_trc_or_ext_AH[$ftr_idx]{$accn}\n");
          #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 7 (trc error, maybe ntr error)\n");
        }
      }
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse esl-epn-translate.pl output line:\n$line\n", 1, $FH_HR);
    }
  }
  close(IN);
  
  return;
}
