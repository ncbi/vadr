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
opt_Add("--keep",      "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",           "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);
opt_Add("--model",      "string",  undef,                    1,    undef, undef,      "use model in file",                            "use model file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--local",      "boolean", 0,                        1,    undef, undef,      "run cmscan locally instead of on farm",        "run cmscan locally instead of on farm", \%opt_HH, \@opt_order_A);
opt_Add("--nseq",       "integer", 5,                        1,    undef,"--local",   "number of sequences for each cmscan farm job", "set number of sequences for each cmscan farm job to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--wait",       "integer", 10,                       1,    undef,"--local",   "allow <n> minutes for cmscan jobs on farm",    "allow <n> minutes for cmscan jobs on farm to finish", \%opt_HH, \@opt_order_A);

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
                'wait=s'        => \$GetOptions_H{"--wait"}); 

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
my $progress_w = 65; # the width of the left hand column in our progress output, hard-coded
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
my @seqname_A   = ();          # actual name of all sequences in fasta file, after being fetched, same order as @accn_A
my %mdl_info_HA = ();          # hash of arrays, values are arrays [0..$nmdl-1];
                               # see dnaorg.pm::validateModelInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAndProcessReferenceSequence()
my %ftr_info_HA = ();          # hash of arrays, values are arrays [0..$nftr-1], 
                               # see dnaorg.pm::validateFeatureInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAndProcessReferenceSequence()

# Call the wrapper function that does the following:
#   1) fetches the sequences listed in @{$accn_AR} into a fasta file and indexes that fasta file,
#      the reference sequence is $accn_AR->[0].
#   2) determines information for each feature (strand, length, coordinates, product) in the reference sequence
#   3) determines type of each reference sequence feature ('cds-mp', 'cds-notmp', or 'mp')
#   4) fetches the reference sequence feature and populates information on the models and features
wrapperFetchAllSequencesAndProcessReferenceSequence(\@accn_A, \@seqname_A, $out_root, \%cds_tbl_HHA,
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
}

#########################################################
# Step 4. Parse homology search results into usable data structuers
##########################################################
#%results_AAH = (); # 1st dim: array, 0..$nmdl-1, one per model
#                   # 2nd dim: array, 0..$naccn-1, one per accessions
#                   # 3rd dim: hash, keys are "p_start", "p_stop", "p_strand", "p_5hangover", "p_3hangover", "p_evalue", "p_fid2ref"
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
parse_cmscan_tblout($tblout_file, \%mdl_info_HA, \@seqname_A, \@accn_A, \%totlen_H, \@results_AAH, $ofile_info_HH{"FH"});

# dump results
dump_results(\@results_AAH, \%mdl_info_HA, \@seqname_A, \@accn_A, \%totlen_H, *STDOUT);

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
#  $seqname_AR:    REF to array of sequence names, PRE-FILLED
#  $accn_AR:       REF to array of accessions, same order as @{$seqname_AR}, PRE-FILLED
#  $totlen_HR:     REF to hash of total lengths of all accessions (keys are values from @{$accn_AR}), PRE-FILLED
#  $results_AAHR:  REF to results AAH, FILLED HERE
#  $FH_HR:         REF to hash of file handles
#
# Returns:    void
#
# Dies:       if we find a hit to a model or sequence that we don't
#             have stored in $mdl_info_HAR or $seqname_AR
#
################################################################# 
sub parse_cmscan_tblout { 
  my $sub_name = "parse_cmscan_tblout()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($tblout_file, $mdl_info_HAR, $seqname_AR, $accn_AR, $totlen_HR, $results_AAHR, $FH_HR) = @_;
  
  # make an 'order hash' for the model names and sequence names,
  my %mdlname_index_H = (); # mdlname_index_H{$model_name} = <n>, means that $model_name is the <n>th model in the @{$mdl_info_HAR{*}} arrays
  my %seqname_index_H = (); # seqname_order_H{$seq_name} = <n>, means that $seq_name is the <n>th sequence name in the @{$seqname_AR}} array
  getIndexHashForArray($mdl_info_HAR->{"cmname"}, \%mdlname_index_H, $FH_HR);
  getIndexHashForArray($seqname_AR, \%seqname_index_H, $FH_HR);

  open(IN, $tblout_file) || fileOpenFailure($tblout_file, $sub_name, $!, "reading", $FH_HR);

  my $did_field_check = 0; # set to '1' below after we check the fields of the file
  my $line_ctr = 0;
  while(my $line = <IN>) { 
    $line_ctr++;
    if($line =~ m/^\# \S/ && (! $did_field_check)) { 
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

      store_hit($results_AAHR, $mdlname, $seqname, $mdllen, $seqlen, $mdlfrom, $mdlto, $from, $to, $strand, $evalue, $FH_HR);
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
#  $seqname_AR:    REF to array of sequence names, PRE-FILLED
#  $accn_AR:       REF to array of accessions, same order as @{$seqname_AR}, PRE-FILLED
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

  my ($results_AAHR, $mdl_info_HAR, $seqname_AR, $accn_AR, $totlen_HR, $FH) = @_;

  my $nmdl = scalar(@{$results_AAHR});
  my $nseq = scalar(@{$seqname_AR}); 
  
  for(my $m = 0; $m < $nmdl; $m++) { 
    for(my $i = 0; $i < $nseq; $i++) { 
      printf $FH ("model($m): %20s  seq($i): %20s  accn: %10s  len: %10d  ", 
                  $mdl_info_HAR->{"cmname"}[$m],
                  $seqname_AR->[$i],
                  $accn_AR->[$i],
                  $totlen_HR->{$accn_AR->[$i]});
      foreach my $key (sort keys %{$results_AAHR->[$m][$i]}) { 
        printf $FH ("%s: %s\n", $key, $results_AAHR->[$m][$i]{$key});
      }
      printf $FH ("\n")
    }
  }
  return;
}
