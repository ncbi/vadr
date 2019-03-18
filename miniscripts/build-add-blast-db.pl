#!/usr/bin/env perl
# EPN, Mon Mar 18 09:51:47 2019
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
# - Reads in a model info file
# - Reads in a stockholm alignment file
# - Creates a fasta file of all CDS sequences in the sequences from the alignment
# - Translates the CDS sequences into protein sequences
# - Creates a blast db from the protein sequences
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
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                                  "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--ttbl",       "integer", 1,                        1,    undef, undef,      "use NCBI translation table <n> to translate CDS",             "use NCBI translation table <n> to translate CDS", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",                          "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $script_name = "build-add-blast-db.pl";
my $usage       = "Usage: $script_name [-options] <input model info file> <input stockholm alignment> <output root for naming output files>\n";
my $synopsis    = "$script_name :: create BLAST db for dnaorg model";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'v'            => \$GetOptions_H{"-v"},
                'ttbl=s'       => \$GetOptions_H{"--ttbl"},
                'keep'         => \$GetOptions_H{"--keep"});

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
if(scalar(@ARGV) != 3) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do dnaorg_build.pl -h\n\n";
  exit(1);
}
my ($in_minfo_file, $in_stk_file, $out_root) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("input model info file", "input stockholm alignment file", "output root for naming output files");
my @arg_A      = ($in_minfo_file, $in_stk_file, $out_root);
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

#############################################################
# make sure the required executables exist and are executable
#############################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"esl-reformat"}  = $esl_exec_dir . "/esl-reformat";
$execs_H{"esl-translate"} = $esl_exec_dir . "/esl-translate";
$execs_H{"makeblastdb"}   = $blast_exec_dir . "/makeblastdb";
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################
# Parse the model info file
###########################
my @mdl_info_AH  = (); # array of hashes with model info
my %ftr_info_HAH = (); # array of hashes with feature info 

my $progress_w = 50;
my $start_secs = outputProgressPrior("Parsing input model info file", $progress_w, $log_FH, *STDOUT);

if(! -e $in_minfo_file) { DNAORG_FAIL("ERROR, model info file $in_minfo_file does not exist", 1, $FH_HR); }
if(! -s $in_minfo_file) { DNAORG_FAIL("ERROR, model info file $in_minfo_file exists but is empty", 1, $FH_HR); }
if(  -d $in_minfo_file) { DNAORG_FAIL("ERROR, model info file $in_minfo_file is actually a directory", 1, $FH_HR); }

modelInfoFileParse($in_minfo_file, \@mdl_info_AH, \%ftr_info_HAH, $FH_HR);
# ensure we read info only a single model from the input file
if(scalar(@mdl_info_AH) == 0) { DNAORG_FAIL("ERROR did not read info for any models in $in_minfo_file", 1, $FH_HR); }
if(scalar(@mdl_info_AH) > 1)  { DNAORG_FAIL("ERROR read info for multiple models in $in_minfo_file, it must have info for only 1 model", 1, $FH_HR); }
if(! defined $mdl_info_AH[0]{"name"}) { 
  DNAORG_FAIL("ERROR did not read name of model in $in_minfo_file", 1, $FH_HR);
}
my $mdl_name = $mdl_info_AH[0]{"name"};
my $mdl_len  = $mdl_info_AH[0]{"length"};

my $ncds = featureInfoCountType(\@{$ftr_info_HAH{$mdl_name}}, "CDS");
if($ncds == 0) { 
  DNAORG_FAIL("ERROR did not read any CDS features in model info file $in_minfo_file", 1, $FH_HR);
}

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

################################
# Parse the input stockholm file
################################
my $fa_file  = $out_root . ".fa";

$start_secs  = outputProgressPrior("Validating input Stockholm file", $progress_w, $log_FH, *STDOUT);

my $stk_mdl_len = stockholm_validate_input($in_stk_file, \%opt_HH, $FH_HR);
if($stk_mdl_len != $mdl_len) { 
  DNAORG_FAIL(sprintf("ERROR, model length inferred from stockholm alignment ($stk_mdl_len) does not match length read from %s ($mdl_len)", 
                      (opt_IsUsed("-m", \%opt_HH) ? "model info file" : "GenBank file")), 1, $FH_HR);
}
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  
$start_secs = outputProgressPrior("Reformatting Stockholm file to FASTA file", $progress_w, $log_FH, *STDOUT);
fastaFileWriteFromStockholmFile($execs_H{"esl-reformat"}, $fa_file, $in_stk_file, \%opt_HH, $FH_HR);

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

######################################################################
# Finish populating @{$ftr_info_HAH{$mdl_name} and create @sgm_info_AH
######################################################################
my @sgm_info_AH = (); # segment info, inferred from feature info

$start_secs = outputProgressPrior("Finalizing feature information", $progress_w, $log_FH, *STDOUT);

featureInfoImputeLength(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
featureInfoImputeSourceIdx(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
featureInfoImputeParentIdx(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);

segmentInfoPopulate(\@sgm_info_AH, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###################
# Translate the CDS
###################
$start_secs = outputProgressPrior("Translating CDS and building BLAST DB", $progress_w, $log_FH, *STDOUT);

my $cds_fa_file  = $out_root . ".cds.fa";
openAndAddFileToOutputInfo(\%ofile_info_HH, "cdsfasta", $cds_fa_file, 1, "fasta sequence file for CDS from $mdl_name");
cdsFetchStockholmToFasta($ofile_info_HH{"FH"}{"cdsfasta"}, $in_stk_file, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
close $ofile_info_HH{"FH"}{"cdsfasta"};

my $protein_fa_file  = $out_root . ".protein.fa";
openAndAddFileToOutputInfo(\%ofile_info_HH, "proteinfasta", $protein_fa_file, 1, "fasta sequence file for translated CDS from $mdl_name");
cdsTranslateToFastaFile($ofile_info_HH{"FH"}{"proteinfasta"}, $execs_H{"esl-translate"}, $cds_fa_file, 
                        $out_root, \@{$ftr_info_HAH{$mdl_name}}, \%opt_HH, $FH_HR);
close $ofile_info_HH{"FH"}{"proteinfasta"};

blastDbProteinCreate($execs_H{"makeblastdb"}, $protein_fa_file, \%opt_HH, $FH_HR);

# add to mdl_info_AH
$mdl_info_AH[0]{"blastdb"} = $protein_fa_file;
if((opt_IsUsed("--ttbl", \%opt_HH)) && (opt_Get("--ttbl", \%opt_HH) != 1))  { 
  $mdl_info_AH[0]{"transl_table"} = opt_Get("--ttbl", \%opt_HH);
}

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

########################
# Output model info file
########################
$start_secs = outputProgressPrior("Creating model info file", $progress_w, $log_FH, *STDOUT);

my $out_minfo_file   = $out_root . ".minfo";
if($in_minfo_file eq $out_minfo_file) { 
  my $old_minfo_file = $in_minfo_file . ".old";
  runCommand("cp $in_minfo_file " . $in_minfo_file . ".old", opt_Get("-v", \%opt_HH), 0, $FH_HR);
  addClosedFileToOutputInfo(\%ofile_info_HH, "oldminfo", $old_minfo_file, 1, "Copy of input model info file");
}  
modelInfoFileWrite($out_minfo_file, \@mdl_info_AH, \%ftr_info_HAH, $FH_HR);
addClosedFileToOutputInfo(\%ofile_info_HH, "outminfo", $out_minfo_file, 1, "Output model info file with blastdb added");

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

##########
# Conclude
##########

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, "", \%ofile_info_HH);
exit 0;

#################################################################
# Subroutine: stockholm_validate_input()
# Incept:     EPN, Fri Mar 15 13:19:32 2019
#
# Synopsis: Validate an input Stockholm file is in the correct 
#           format and has all the information we need.
#
# Arguments:
#  $in_stk_file:  input stockholm file to validate
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    Length of the model that will be built from this alignment.
#             If RF annotation exists, return nongap RF length
#             If RF annotation does not exist, we verify that there is only 1 sequence and no gaps
#             so we return alignment length
#
# Dies:       if there's a problem parsing the file or 
#             a requirement is not met
#################################################################
sub stockholm_validate_input {
  my $sub_name = "stockholm_validate_input";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($in_stk_file, $opt_HHR, $FH_HR) = @_;

  if(! -e $in_stk_file) { DNAORG_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file does not exist", 1, $FH_HR); }
  if(! -s $in_stk_file) { DNAORG_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file exists but is empty", 1, $FH_HR); }
  if(  -d $in_stk_file) { DNAORG_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file is actually a directory", 1, $FH_HR); }
  my $msa = Bio::Easel::MSA->new({ fileLocation => $in_stk_file, isDna => 1});
  my $nseq = $msa->nseq;
  if($nseq > 1) { 
    # multiple sequences in the alignment
    # this only works if the alignment has RF annotation
    if(! $msa->has_rf) { 
      DNAORG_FAIL("ERROR, read more than 1 ($nseq) sequences in --stk file $in_stk_file, this is only allowed ifalignment has RF annotation (cmbuild -O will give you this)", 1, $FH_HR);
    }
  }
  elsif($nseq == 1) { 
    # a single sequence in the alignment, it must have zero gaps
    if($msa->any_allgap_columns) { 
      DNAORG_FAIL("ERROR, read 1 sequence in --stk file $in_stk_file, but it has gaps, this is not allowed for single sequence 'alignments' (remove gaps with 'esl-reformat --mingap')", 1, $FH_HR);
    }
  }
  else { # nseq == 0
    DNAORG_FAIL("ERROR, did not read any sequence data in --stk file $in_stk_file", 1, $FH_HR);
  }

  # determine return value
  my $ret_mdl_len = undef;
  if($msa->has_rf) { 
    $ret_mdl_len = $msa->get_rflen;
  }
  else { # no RF annotation, we checked above there was a single sequence with zero gaps
    $ret_mdl_len = $msa->alen;
  }
  return $ret_mdl_len;
}

