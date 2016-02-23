#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Thu Feb 18 12:48:16 2016 [dnaorg_annotate.pl split off from dnaorg_annotate_genomes.pl]

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

require "dnaorg.pm"; 
require "epn-options.pm";

# hard-coded-paths:
my $inf_exec_dir       = "/usr/local/infernal/1.1.1/bin/";
my $esl_fetch_cds      = "/panfs/pan1/dnaorg/programs/esl-fetch-cds.pl";
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
#     option            type       default               group   requires incompat preamble-output                          help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,   undef,                                   "display this help",                                  \%opt_HH, \@opt_order_A);
opt_Add("-d",           "string",  0,                        1,    undef, undef,   "directory specified as",                "specify output directory is <s> (created with dnaorg_build.pl -d <s>), not <reference accession>\n", \%opt_HH, \@opt_order_AA);
opt_Add("-c",           "boolean", 0,                        1,    undef, undef,   "genome is circular",                    "genome is circular",                                 \%opt_HH, \@opt_order_A);
opt_Add("-matpept",     "string",  undef,                    1,    undef, undef,   "using pre-specified mat_peptide info",  "read mat_peptide info in addition to CDS info, file <f> explains CDS:mat_peptide relationships", \%opt_HH, \@opt_order_A);

#$opt_group_desc_H{"2"} = "options affecting window/hit definition";
##       option       type       default                group  requires incompat  preamble-output                          help-output    
#opt_Add("--b1",      "real",    undef,                     2,    undef, undef,   "bit score threshold for window defn",   "set bit score threshold for window definition to <x>", \%opt_HH, \@opt_order_A);
#opt_Add("--b2",      "real",    undef,             ,       2,    undef, "--e2",  "bit score threshold for hit defn",      "set bit score threshold for hit definition to <x>",    \%opt_HH, \@opt_order_A);
#opt_Add("--e2",      "real",    $RNAVORE_H{"DF_E2"},       2,    undef, "--b2",  "E-value threshold for hit defn",        "set E-value   threshold for hit definition to <x>",    \%opt_HH, \@opt_order_A);
#opt_Add("--toponly", "boolean", 0,                         2,    undef, undef,   "only allow hits on top-strand",         "only search top-strand of target sequence file",       \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_annotate.pl [-options] <file with list of accessions to annotate>\n";
my $synopsis = "dnaorg_annotate.pl :: annotate sequences based on a reference annotation";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'd=s'          => \$GetOptions_H{"-d"},
                'c'            => \$GetOptions_H{"-c"},
                'matpept=s'    => \$GetOptions_H{"-matpept"});


# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $synopsis);
  outputBanner(*STDOUT, $synopsis, $command, $date, $opts_used_long);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 1) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do rvr-align -h\n\n";
  exit(1);
}
my ($listfile) = (@ARGV);

###############
# Preliminaries
###############
# first, parse the list file, we need to do this first because we need
# to know what the reference accession <refaccn> is to check if the
# directory <refaccn> exists
my @accn_A        = (); # array of accessions
parseListFile($listfile, 1, \@accn_A, undef); # 1 
$ref_accn = $accn_A[0];

my $dir_set_as_ref_accn = 0;
if(! defined $dir) { 
  $dir = $ref_accn;
  $dir_set_as_ref_accn = 1;
}

# make sure that $dir exists
if(! -d $dir) {
  DNAORG_FAIL(sprintf("ERROR, directory $dir %s does not exist", $dir_set_as_ref_accn ? "(first accession read from $listfile)" : "(specified with -d)"), 1, undef);
}

#############################################
# output program banner and open output files
#############################################
my $synopsis = "dnaorg_annotate.pl :: annotate sequences using reference homology models";
my $command  = "$executable $opts_used_short $ref_accn";
my $date     = scalar localtime();
outputBanner(*STDOUT, $synopsis, $command, $date, $opts_used_long);

# open the log and command files:
# set output file names and file handles, and open those file handles
my %ofile_name_H          = (); # full name for (full path to) output files
my %ofile_sname_H         = (); # short name for (no dir path) output files
my %ofile_FH_H            = (); # file handle for output files, keys are in @ofile_keys_A
my @ofile_keys_A          = ("log", "cmd"); 
my %ofile_desc_H          = (); # description of each output file
$ofile_desc_H{"log"}      = "Output printed to screen";
$ofile_desc_H{"cmd"}      = "List of executed commands";

foreach my $key (@ofile_keys_A) { 
  $ofile_name_H{$key}  = $out_root . "." . $key;
  $ofile_sname_H{$key} = $dir_tail . ".dnaorg_annotate." . $key; # short name (lacks directory)
  if(! open($ofile_FH_H{$key}, ">", $ofile_name_H{$key})) { 
    printf STDERR ("ERROR, unable to open $ofile_name_H{$key} for writing.\n"); 
    exit(1);
  }
}
my $log_FH = $ofile_FH_H{"log"};
my $cmd_FH = $ofile_FH_H{"cmd"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# now we have the log file open, output the banner there too
outputBanner($log_FH, $synopsis, $command, $date, $opts_used_long);

############################################
# parse the optional input files, if nec
###########################################
# -matpept <f>
my @cds2pmatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of primary matpept indices that comprise this CDS
my @cds2amatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of all     matpept indices that comprise this CDS
if(defined $matpept_infile) { 
  parseMatPeptSpecFile($matpept_infile, \@cds2pmatpept_AA, \@cds2amatpept_AA, \%ofile_FH_H);
}
# -specstart <f>
my @specstart_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of allowed start codons for this CDS
if(defined $specstart_infile) { 
  parseSpecStartFile($specstart_infile, \@specstart_AA, \%ofile_FH_H);
}

###########################################################################
# Step 1. Gather and process information on all genomes using Edirect.
###########################################################################
my $progress_w = 60; # the width of the left hand column in our progress output, hard-coded
my $start_secs = outputProgressPrior("Gathering information on all sequences using edirect", $progress_w, $log_FH, *STDOUT);
my %cds_tbl_HHA = ();   # CDS data from .cds.tbl file, hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %mp_tbl_HHA = ();    # mat_peptide data from .matpept.tbl file, hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column

# 1) create the edirect .mat_peptide file, if necessary
# 2) create the edirect .ftable file
# 3) create the length file
# 
# We create the .mat_peptide file first because we will die with an
# error if mature peptide info exists and neither -matpept nor
# -nomatpept was used (and we want to die as early as possible in the
# script to save the user's time)
#
# 1) create the edirect .mat_peptide file, if necessary
my $mp_file = $out_root . ".mat_peptide";

#      if -nomatpept was   enabled we don't attempt to create a matpept file
# else if -matpept was     enabled we validate that the resulting matpept file is not empty
# else if -matpept was not enabled we validate that the resulting matpept file is     empty
if(! $do_nomatpept) { 
#  $cmd = "esearch -db nuccore -query $ref_accn | efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product > $mp_file";
  $cmd = "cat $listfile | epost -db nuccore -format acc | efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product > $mp_file";
  runCommand($cmd, $be_verbose, \%ofile_FH_H);
  
  if($do_matpept) { 
    if(! -s  $mp_file) { 
      DNAORG_FAIL("ERROR, -matpept enabled but no mature peptide information exists.", 1, \%ofile_FH_H); 
    }
    addClosedOutputFile(\@ofile_keys_A, \%ofile_name_H, \%ofile_sname_H, \%ofile_desc_H, "mp", $mp_file, "Mature peptide information obtained via edirect", \%ofile_FH_H);
  }
  else { # ! $do_matpept
    if(-s $mp_file) { 
      DNAORG_FAIL("ERROR, -matpept not enabled but mature peptide information exists, use -nomatpept to ignore it.", 1, \%ofile_FH_H); 
    }
    else { 
      # remove the empty file we just created
      runCommand("rm $mp_file", $be_verbose, \%ofile_FH_H);
    }
  }
}
                
# 2) create the edirect .ftable file
# create the edirect ftable file
my $ft_file  = $out_root . ".ftable";
#$cmd = "esearch -db nuccore -query $ref_accn | efetch -format ft > $ft_file";
$cmd = "cat $listfile | epost -db nuccore -format acc | efetch -format ft > $ft_file";
                runCommand($cmd, $be_verbose, \%ofile_FH_H);
addClosedOutputFile(\@ofile_keys_A, \%ofile_name_H, \%ofile_sname_H, \%ofile_desc_H, "ft", $ft_file, "Feature table obtained via edirect", \%ofile_FH_H);

# 3) create the length file
# create a file with total lengths of each accession
my $len_file  = $out_root . ".length";
my $len_file_created = $len_file . ".created";
                my $len_file_lost    = $len_file . ".lost";
#$cmd = "esearch -db nuccore -query $ref_accn | efetch -format gpc | xtract -insd INSDSeq_length | grep . | sort > $len_file";
                $cmd = "cat $listfile | epost -db nuccore -format acc | efetch -format gpc | xtract -insd INSDSeq_length | grep . | sort > $len_file";
                runCommand($cmd, $be_verbose, \%ofile_FH_H);
addClosedOutputFile(\@ofile_keys_A, \%ofile_name_H, \%ofile_sname_H, \%ofile_desc_H, "len", $len_file, "Sequence length file", \%ofile_FH_H);
if(! -s $len_file) { 
  DNAORG_FAIL("ERROR, no length information obtained using edirect.", 1, \%ofile_FH_H); 
}  

exit 0;

