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
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",                  "force; if dir <reference accession> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                   "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--dirout",     "string",  undef,                    1,    undef, undef,      "output directory specified as <s>",            "specify output directory as <s>, not <ref accession>", \%opt_HH, \@opt_order_A);
opt_Add("--matpept",    "string",  undef,                    1,    undef, undef,      "using pre-specified mat_peptide info",         "read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships", \%opt_HH, \@opt_order_A);
opt_Add("--nomatpept",  "boolean", 0,                        1,    undef,"--matpept", "ignore mat_peptide annotation",                "ignore mat_peptide information in reference annotation", \%opt_HH, \@opt_order_A);
opt_Add("--xfeat",      "string",  undef,                    1,    undef, undef,      "build models of additional qualifiers",        "build models of additional qualifiers in string <s>", \%opt_HH, \@opt_order_A);  
opt_Add("--dfeat",      "string",  undef,                    1,    undef, undef,      "annotate additional qualifiers as duplicates", "annotate qualifiers in <s> from duplicates (e.g. gene from CDS)",  \%opt_HH, \@opt_order_A);  
opt_Add("--keep",       "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",           "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"2"} = "optional output files";
#       option       type       default                group  requires incompat  preamble-output                          help-output    
opt_Add("--seginfo",    "boolean", 0,                        2,    undef, undef, "output internal model information",     "create file with internal model information",   \%opt_HH, \@opt_order_A);
opt_Add("--ftrinfo",    "boolean", 0,                        2,    undef, undef, "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"3"} = "options for skipping stages and using files from an earlier, identical run, primarily useful for debugging";
#     option               type       default               group   requires    incompat                  preamble-output                                            help-output    
opt_Add("--skipedirect",   "boolean", 0,                       3,   undef,      undef,                    "skip the edirect steps, use existing results",           "skip the edirect steps, use data from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipfetch",     "boolean", 0,                       3,   undef,      undef,                    "skip the sequence fetching steps, use existing results", "skip the sequence fetching steps, use files from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipbuild",     "boolean", 0,                       3,   undef,      undef,                    "skip the build/calibrate steps",                         "skip the model building/calibrating, requires --seginfo and/or --ftrinfo", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_build.pl [-options] <reference accession>\n";
my $synopsis = "dnaorg_build.pl :: build homology models for features of a reference sequence";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'dirout=s'     => \$GetOptions_H{"--dirout"},
                'matpept=s'    => \$GetOptions_H{"--matpept"},
                'nomatpept'    => \$GetOptions_H{"--nomatpept"},
                'xfeat=s'      => \$GetOptions_H{"--xfeat"},
                'dfeat=s'      => \$GetOptions_H{"--dfeat"},
                'keep'         => \$GetOptions_H{"--keep"},
# optional output files
                'seginfo'      => \$GetOptions_H{"--seginfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"},
# options for skipping stages, using earlier results
                'skipedirect'  => \$GetOptions_H{"--skipedirect"},
                'skipfetch'    => \$GetOptions_H{"--skipfetch"},
                'skipbuild'    => \$GetOptions_H{"--skipbuild"});

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
   (! (opt_Get("--seginfo", \%opt_HH) || opt_Get("--ftrinfo", \%opt_HH)))) { 
  die "ERROR, --skipbuild requires one or both of --seginfo or --ftrinfo"; 
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
my $FH_HR  = $ofile_info_HH{"FH"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# open optional output files
if(opt_Get("--seginfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "seginfo", $out_root . ".seginfo", 1, "Model information (created due to --seginfo)");
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

###################################################
# make sure the required executables are executable
###################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmbuild"}       = $inf_exec_dir . "cmbuild";
$execs_H{"esl-reformat"}  = $esl_exec_dir . "esl-reformat";
#$execs_H{"esl_fetch_cds"} = $esl_fetch_cds;
$execs_H{"makeblastdb"}   = $blast_exec_dir . "makeblastdb";
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

########################
# Fetch the genbank file
########################
my $progress_w = 50; # the width of the left hand column in our progress output, hard-coded
my $start_secs = outputProgressPrior("Fetching GenBank file", $progress_w, $log_FH, *STDOUT);

my $gb_file  = $out_root . ".gb";
edirectFetchToFile($gb_file, $ref_accn, "gb", 5, $ofile_info_HH{"FH"});  # number of attempts to fetch to make before dying
addClosedFileToOutputInfo(\%ofile_info_HH, "gb", $gb_file, 1, "GenBank format file for $ref_accn");

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

########################
# Parse the genbank file
########################
my %seq_info_HH  = ();
my %full_ftr_info_HAH = ();
my @ftr_info_AH = ();

$start_secs = outputProgressPrior("Parsing GenBank file", $progress_w, $log_FH, *STDOUT);

genbankParse($gb_file, \%seq_info_HH, \%full_ftr_info_HAH, $FH_HR);
if(! exists $full_ftr_info_HAH{$ref_accn}) { 
  DNAORG_FAIL("ERROR parsing GenBank file $gb_file, did not read info for reference accession $ref_accn\n", 1, $FH_HR);
}

# remove any features that are not of the type that we care about
# and from the features we do care about, remove keys we don't care about

# TODO, optionally prune to keep certain feature types, or all (?)
# TODO, optionally prune to keep certain qualifiers, or all (?)

my $full_ftr_idx;
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

# copy subset of information from @{%full_ftr_info_HAH{$ref_accn}} to create %ftr_info_AH
$ftr_idx = -1;
for($full_ftr_idx = 0; $full_ftr_idx < scalar(@{$full_ftr_info_HAH{$ref_accn}}); $full_ftr_idx++) { 
  my $full_ftype = $full_ftr_info_HAH{$ref_accn}[$full_ftr_idx]{"type"};
  if((defined $full_ftype) && (exists $ftype_H{$full_ftype})) { 
    $ftr_idx++;
    %{$ftr_info_AH[$ftr_idx]} = ();
    foreach $key (sort keys %{$full_ftr_info_HAH{$ref_accn}[$full_ftr_idx]}) { 
      if((exists $qual_H{$key}) && (defined $full_ftr_info_HAH{$ref_accn}[$full_ftr_idx]{$key})) { 
        $ftr_info_AH[$ftr_idx]{$key} = $full_ftr_info_HAH{$ref_accn}[$full_ftr_idx]{$key};
      }
    }
  }
}
%full_ftr_info_HAH = (); # we don't need this any more

#############################################################
# Finish populating feature_info_HA and segment_info_HA hashes
# and output model info file.
#############################################################

my @seg_info_AH = ();          # hash of arrays, values are arrays [0..$nftr-1], 
                               # see dnaorg.pm::validateSegmentInfoHashIsComplete() for list of all keys

featureInfoImputeCoords(\@ftr_info_AH, \%opt_HH, $FH_HR);
featureInfoImputeSourceIdx(\@ftr_info_AH, \%opt_HH, $FH_HR);
featureInfoImputeParentIdx(\@ftr_info_AH, \%opt_HH, $FH_HR);

segmentInfoPopulate(\@seg_info_AH, \@ftr_info_AH, \%opt_HH, $FH_HR);

my $minfo_file  = $out_root . ".minfo";
output_model_info_file($minfo_file, $ref_accn, \@ftr_info_AH, $FH_HR);
addClosedFileToOutputInfo(\%ofile_info_HH, "minfo", $minfo_file, 1, "GenBank format file for $ref_accn");

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################
# Step 4. Fetch and process the reference genome sequence
##########################################################
# verify our model and feature info hashes are complete, 
# if validateFeatureInfoHashIsComplete() fails then the program will exit with an error message

$start_secs = outputProgressPrior("Creating FASTA and STOCHOLM sequence files", $progress_w, $log_FH, *STDOUT);

my $fa_file  = $out_root . ".fa";
openAndAddFileToOutputInfo(\%ofile_info_HH, "fasta", $fa_file, 1, "fasta sequence file for $ref_accn");
print_sequence_to_fasta_file($ofile_info_HH{"FH"}{"fasta"}, 
                             $seq_info_HH{$ref_accn}{"ver"}, 
                             $seq_info_HH{$ref_accn}{"def"}, 
                             $seq_info_HH{$ref_accn}{"seq"}, $FH_HR);
close $ofile_info_HH{"FH"}{"fasta"};

my $stk_file = $out_root . ".stk";
reformat_fasta_file_to_stockholm_file($execs_H{"esl-reformat"}, $fa_file, $stk_file, \%opt_HH, $FH_HR);
addClosedFileToOutputInfo(\%ofile_info_HH, "stk", $stk_file, 1, "Stockholm alignment file for $ref_accn");

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
exit 0;
#######################################################################
# Step 5. Fetch the CDS protein translations and build BLAST database
#######################################################################
$start_secs = outputProgressPrior("Fetching protein translations of CDS and building BLAST DB", $progress_w, $log_FH, *STDOUT);
my @prot_fa_file_A = ();
#fetch_proteins_into_fasta_files($out_root, $ref_accn, \@ftr_info_AH, \@prot_fa_file_A, \%opt_HH, \%ofile_info_HH);

foreach my $prot_fa_file (@prot_fa_file_A) { 
#  create_blast_protein_db(\%execs_H, $prot_fa_file, \%opt_HH, \%ofile_info_HH);
}
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################
# Step 6. Build the model
#########################
if(! opt_Get("--skipbuild", \%opt_HH)) { 
  $start_secs = outputProgressPrior("Building model (this could take a while)", $progress_w, $log_FH, *STDOUT);

  my $cmbuild_cmd  = $execs_H{"cmbuild"} . " --noss --verbose -F $out_root.cm " . $ofile_info_HH{"fullpath"}{"stk"} . " > $out_root.cmbuild";
  runCommand($cmbuild_cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});

  outputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);

  addClosedFileToOutputInfo(\%ofile_info_HH, "cm", "$out_root.cm", 1, "CM file");
}

##########
# Conclude
##########
# output optional output files
if(exists $ofile_info_HH{"FH"}{"ftrinfo"}) { 
  dumpArrayOfHashes("Feature information (@ftr_info_AH) for $ref_accn", $ofile_info_HH{"FH"}{"ftrinfo"});
}
if(exists $ofile_info_HH{"FH"}{"seginfo"}) { 
  dumpArrayOfHashes("Segment information (@seg_info_AH) for $ref_accn", $ofile_info_HH{"FH"}{"seginfo"});
}

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir, \%ofile_info_HH);
exit 0;

#################################################################
# Subroutine: translate_proteins_into_fasta_files()
# Incept:     EPN, Mon Mar 11 06:30:55 2019
# 
# Purpose:    Use esl-translate to create protein translations 
#             of each of the N CDS for an accession and create 
#             N+1 FASTA files, one with each single sequence
#             and one with all sequences.
#             Fill @{$fa_file_AR} with the sequence file names.
#
# Arguments:
#   $execs_HR:       hash with paths to esl-sfetch and esl-translate executables
#   $out_root:       string for naming output files
#   $ref_accn:       reference accession
#   $ftr_info_AHR:   REF to the feature info, pre-filled
#   $seg_info_HAR:   REF to the segment info, pre-filled
#   $fa_file_AR:     REF to array of fasta file names, filled here 
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $ofile_info_HHR: REF to the 2D hash of output file information
#                    
# Returns: void
# Dies:    if a fetched location for a feature does not match to any feature's "ref_coords" 
#
#################################################################
sub translate_into_protein_fasta_files { 
  my $sub_name = "translate_into_protein_fasta_files";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $out_root, $ref_accn, $ftr_info_HAR, $seg_info_HAR, $fa_file_AR, $opt_HHR, $ofile_info_HHR) = @_;
  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience

  my $nftr = getSizeOfInfoHashOfArrays(\@ftr_info_AH, "type", $FH_HR);

  my $sfetch_out_file  = $out_root . ".prot.sfetch";
  my $all_fa_out_file  = $out_root . ".prot.fa";

#  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
#    if(featureTypeIsCds($ftr_info_HAR, $ftr_idx)) { 
#      push(@fetch_name_A,  $ref_accn);
#      push(@fetch_coords_A,  $ftr_info_HAR->{"coords"}[$ftr_idx]);
#      # HERE HERE HERE 
#    }
#  }

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
  open(ALLFA, ">", $all_fa_out_file) || fileOpenFailure($all_fa_out_file,  $sub_name, $!, "writing", $FH_HR);
  while(my $line = <IN>) { 
    chomp $line;
    my @el_A = split(/\t/, $line);
    if(scalar(@el_A) != 4) { 
      DNAORG_FAIL("ERROR in $sub_name, not exactly 4 tab delimited tokens in efetch output file line\n$line\n", 1, $FH_HR);
    }
    my ($read_ref_accn, $prot_accn, $location, $translation) = (@el_A);
    $location =~ s/\.\./\-/g;
    my $new_name = $prot_accn . "/" . $location;

    print ALLFA  (">$new_name\n$translation\n");

    # determine what feature this corresponds to, and create the individual fasta file for that
    my $ftr_idx = blastxDbSeqNameToFtrIdx($new_name, $ftr_info_HAR, $FH_HR); # this will die if we can't find the feature with $location
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

#################################################################
# Subroutine:  fetch_reference_info_edirect()
# Incept:      EPN, Thu Mar  7 14:17:58 2019 [updated]
#              EPN, Tue Feb 23 13:00:23 2016
#
# Purpose:     Gets sequence information using Edirect and parse that
#              Edirect output into usable data structures.
#
#              Creates the following output files and stores
#              information on them in %{$ofile_info_HHR}
#              by calling the addClosedFileToOutputInfo() function:
#              - $out_root . ".mat_peptide": mature peptide info obtained via edirect
#              - $out_root . ".ftable":      feature table obtained via edirect
#                      
# Arguments: 
#   $ref_accn:              reference accession
#   $feature_str:           comma separated string of features to parse, e.g. "CDS,gene,RNA"
#   $out_root:              string that is the 'root' for naming output files
#   $feat_tbl_HHAR:         REF to hash of hash of arrays for other (non-CDS and non-MP) feature info, FILLED HERE
#   $ofile_info_HHR:        REF to 2D hash with output info, ADDED TO HERE
#   $opt_HHR:               REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#   $FH_HR:                 REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#
# Returns:     void
#
# Dies:        never
#     
################################################################# 
sub fetch_reference_info_edirect {
  my $sub_name = "fetch_reference_info_edirect";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name, entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($ref_accn, $feature_str, $out_root, $feat_tbl_HHAR, $ofile_info_HHR, $opt_HHR, $FH_HR) = @_;

  my $cmd; 
  my @feature_A = ();
  if((defined $feature_str) && ($feature_str ne "")) { 
    @feature_A = split(",", $feature_str);
  }  

  # mature peptides
  my $mp_file = $out_root . ".mat_peptide";
  $cmd = "esearch -db nuccore -query $ref_accn | efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product > $mp_file";
  runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  if(-s $mp_file) { 
    addClosedFileToOutputInfo($ofile_info_HHR, "mp", $mp_file, 0, "Mature peptide information obtained via edirect");
    edirect_ftable_or_matpept_to_single_feature_table_info($ref_accn, $mp_file, 1, "mat_peptide", \%{$feat_tbl_HHAR->{"mat_peptide"}}, $FH_HR); # 1: it is a mat_peptide file
  }
  else { 
    # remove the empty file we just created
    runCommand("rm $mp_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);
  }

  # feature table
  my $ft_file  = $out_root . ".fetched.ftable";
  $cmd = "esearch -db nuccore -query $ref_accn | efetch -format ft > $ft_file";
  runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  if(! -s $ft_file) { 
    DNAORG_FAIL("ERROR in $sub_name, unable to fetch feature table", 1, $FH_HR); 
  }
  addClosedFileToOutputInfo($ofile_info_HHR, "ft", $ft_file, 0, "Feature table obtained via edirect");

  # features from feature table
  foreach my $feat (@feature_A) { 
    edirect_ftable_or_matpept_to_single_feature_table_info($ref_accn, $ft_file, 0, $feat, \%{$feat_tbl_HHAR->{$feat}}, $FH_HR); # 0: it's not a mat_peptide file
  }    

  return 0;
}

#################################################################
# Subroutine: edirect_ftable_or_matpept_to_single_feature_table_info()
# Incept:     EPN, Fri Feb 12 09:48:00 2016
# 
# Purpose:    Given an edirect output feature table file or mature
#             peptide table file, parse it and store it's relevant
#             information in a 'single feature table' (%{$tbl_HAR}).
#
#             This is a wrapper subroutine that calls 
#             parse_edirect_mat_peptide_file() or 
#             parse_edirect_ftable_file() and then 
#             get_single_feature_table_info().
#
# Arguments:
#   $ref_accn:      reference accession
#   $edirect_file:  name of edirect output file to parse
#   $do_matpept:    '1' if edirect file is a mature peptide file, '0' or undef if it is a ftable file
#   $feature:       the feature we want to store info on in $tbl_HAR (e.g. "CDS" or "mat_peptide")
#   $tbl_HAR:       REF to hash of hash of arrays we'll fill with info on $qual_name:
#                   1D: key: qualifier (e.g. 'coords')
#                   2D: values for each qualifier, size will be number of features for this accession
#                       e.g. size of 5 means this accession has 5 CDS if $feature is CDS
#                       INITIALIZED AND FILLED HERE
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if we're not able to parse the edirect file
#################################################################
sub edirect_ftable_or_matpept_to_single_feature_table_info { 
  my $sub_name = "edirect_ftable_or_matpept_to_single_feature_table_info";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($ref_accn, $edirect_file, $do_matpept, $feature, $tbl_HAR, $FH_HR) = @_;

  # variables relevant to parsing the edirect ftable and edirect mat_peptide files
  my %quals_HHA       = ();      # values that will go into a table
                                 # key 1: feature name, e.g. "CDS"
                                 # key 2: 'fac' string: <full_accession><fac_seq><coordinates>
                                 # value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>
  my @faccn_A         = ();      # array of all full accessions read
  my %fac_HHA         = ();      # used to easily determine list of 2D keys ('fac's) in quals_HHA 
                                 # for a given feature and faccn.
                                 # key 1: feature name: e.g. "CDS"
                                 # key 2: 'faccn', full accession
                                 # value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
  my %faccn2accn_H    = ();      # key: full accession, value: short accession, used for convenience
  my %column_HA       = ();      # key: feature, e.g. "CDS"
                                 # value: array of qualifiers that exist for that feature in at least 1 faccn
  my $dummy_column    = "DuMmY"; # name of 'dummy' column that we don't actually print, for use with --matpept only
  my %sep_H           = ();      # hash of 'separation' values we use to concatenate multiple strings into 
                                 # one, to ease with their storage
  $sep_H{"qnqv"}      = "!!";    # string that separates 'qualifier name' and 'qualifier_values' 
  $sep_H{"fac"}       = "::";    # strings separating full accessions, start and stop coords
  $sep_H{"qval"}      = ";;";    # strings separating qualifier values  

  %{$tbl_HAR} = ();

  if((defined $do_matpept) && ($do_matpept)) { 
    parse_edirect_mat_peptide_file($edirect_file, $dummy_column, \%sep_H, \%quals_HHA, \@faccn_A, \%fac_HHA, \%faccn2accn_H, \%column_HA, $FH_HR);
  }
  else { # default, a ftable file
    parse_edirect_ftable_file($edirect_file, undef, \%sep_H, \%quals_HHA, \@faccn_A, \%fac_HHA, \%faccn2accn_H, \%column_HA, $FH_HR);
  }
  if(exists $quals_HHA{$feature}) { 
    get_single_feature_table_info($ref_accn, $dummy_column, \%sep_H, \%{$quals_HHA{$feature}}, \@faccn_A, \%{$fac_HHA{$feature}}, \%faccn2accn_H, \@{$column_HA{$feature}}, $tbl_HAR, $FH_HR);
  }
  
  # if $quals_HHA doesn't exist, return 
  return;
}

#################################################################
# Subroutine: parse_edirect_ftable_file()
# Incept:     EPN, Thu Feb 11 13:47:11 2016
# 
# Purpose:   Given a feature table file output from an 'efetch -format
#            ft' command, parse that file into usable data
#            structures. 
# 
#            Can be called by the wrapper subroutine: 
#             edirectFtableOrMatPept2SingleFeatureTableInfo().
#
#            Caller will commonly call get_single_feature_table_info() after calling this
#            subroutine.
#
#            Similar to (and shares some code with) parse_edirect_mat_peptide_file().
#
# Args:
#   $infile:          name of ftable file
#   $only_feature:    name of only feature we will parse, 'undef' to parse all features
#   $sep_HR:          ref to hash of 'separation' values, keys: "qnqv", "fac", and "qval", 
#                     used to split values that have been concatenated with these 'separation'
#                     values in between, PRE-FILLED
#   $quals_HHAR:      ref to 2D hash of arrays with values that could eventually go 
#                     into a single-feature table, FILLED HERE
#                     key 1: feature name, e.g. "CDS"
#                     key 2: 'fac' string: <full_accession><fac_seq><coordinates>
#                     value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>
#   $faccn_AR:        REF to array of full accessions read, FILLED HERE
#   $fac_HHAR:        REF to 2D hash of arrays that is used to easily determine
#                     list of 2D keys ('fac's) in quals_HHA for a given feature and faccn
#                     FILLED HERE
#                     key 1: feature name: e.g. "CDS"
#                     key 2: 'faccn', full accession
#                     value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
#   $faccn2accn_HR:   REF to hash, key: full accession, value: short accession, 
#                     used for convenience in output function, FILLED HERE
#   $column_HAR:      REF to hash of arrays, key: feature, e.g. "CDS"
#                     value: array of qualifiers that exist for that feature in at least 1 faccn
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
#
# Returns:  void
#
# Dies:     if <$infile> file is in unexpected format
#################################################################
sub parse_edirect_ftable_file {
  my $sub_name = "parse_edirect_ftable_file()";

  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $only_feature, $sep_HR, $quals_HHAR, $faccn_AR, $fac_HHAR, $faccn2accn_HR, $column_HAR, $FH_HR) = @_;

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  my $do_only_feature  = 0;
  my $cap_only_feature = undef; # capitalized version of $only_feature
  if(defined $only_feature) { 
    $do_only_feature = 1;
    $cap_only_feature = $only_feature;
    $cap_only_feature =~ tr/a-z/A-Z/;
  }

  # variables used to make sure input feature table is in 
  # expected format w.r.t to order of line types

  # define other variables we'll use
  my $faccn           = undef;   # full accession, e.g. ref|NC_000883.2|
  my $fac             = undef;   # full accession and coords, accession and each segment coords separated by $fac_sep
  my $faccn2          = undef;   # another full accession
  my $accn            = undef;   # full accession with extraneous characters removed: e.g. NC_000883.2
  my $start           = undef;   # start position
  my $stop            = undef;   # stop position
  my $qname           = undef;   # a qualifier name,  e.g. 'organism'
  my $qval            = undef;   # a qualifier value, e.g. 'Paramecia bursaria Chlorella virus 1'
  my $feature         = undef;   # a feature name, e.g. "CDS", or "gene"
  my $cap_feature     = undef;   # a feature name capitalized, necessary so -f option can be case insensitive
  my $coords          = undef;   # coordinates
  my $fac_sep  = $sep_HR->{"fac"};
  my $qnqv_sep = $sep_HR->{"qnqv"};
  my $qval_sep = $sep_HR->{"qval"};
  my $line_ctr        = 0;       # count of number of lines read in ftable
  my $accn_ctr        = 0;       # count of number of accessions read
  my %column_HH       = ();      # existence 2D hash, used to aid construction of %column_HA
                                 # key 1: feature, e.g. "CDS"
                                 # key 2: qualifier
                                 # value: '1' (always)
  my $prv_was_accn           = 0; # set to '1' if previous line was an accession line
  my $prv_was_coords_feature = 0; # set to '1' if previous line was a coordinates line with a feature name
  my $prv_was_coords_only    = 0; # set to '1' if previous line was a coordinates line without a feature name
  my $prv_was_quals          = 0; # set to '1' if previous line was a qualifier_name qualifier value line

  while(my $line = <IN>) { 
    $line_ctr++;
    chomp $line;
    if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
    if($line =~ m/\w/) { 
      # parse each of the 4 line types differently
      # -------------------------------------------------------
      if($line =~ /Feature\s+(\S+)$/) { 
        # ACCESSION LINE
        # example:
        #>Feature ref|NC_001359.1|    
        $faccn = $1;
        $accn_ctr++;
        # does line order make sense?
        if($accn_ctr == 1   ||          # first accession of the file
           $prv_was_quals ||          # previous line was a quals line (common)
           $prv_was_coords_feature || # previous line was a coords feature line (rare)
           $prv_was_coords_only    || # previous line was a coords line without a feature (rare)
           $prv_was_accn) {           # previous line was an accession line (rare)
          # line order makes sense, keep going...

          # determine short version of the accession, e.g. NC_001359 in above example
          if($faccn =~ /[^\|]*\|([^\|]*)\|/) { 
            $accn = $1;
            $faccn2accn_HR->{$faccn} = $accn;
            push(@{$faccn_AR}, $faccn);
          }
          else { 
            DNAORG_FAIL("ERROR in $sub_name, unable to parse Feature line $line", 1, $FH_HR); 
          }
          $feature     = undef; 
          $cap_feature = undef;
          $fac         = undef;

          # update '$prv_*' values that we use to make sure line order makes sense
          $prv_was_accn           = 1;
          $prv_was_coords_feature = 0;
          $prv_was_coords_only    = 0;
          $prv_was_quals          = 0;
          #printf("set prv_was_accn\n");
        }
        else { # line order is unexpected
          DNAORG_FAIL("ERROR in $sub_name, unexpected line order (accession line) at line $line_ctr: $line", 1, $FH_HR);
        }
      }
      # -------------------------------------------------------
      elsif($line =~ /^(\<?\d+\^?)\s+(\>?\d+)\s+(\S+)$/) { 
        # COORDINATES LINE WITH A FEATURE NAME
        # example:
        #230	985	gene
        ($start, $stop, $feature) = ($1, $2, $3);

        # does line order make sense?
        if($prv_was_accn  ||           # previous line was accession line (common)
           $prv_was_quals ||           # previous line was quals line (common)
           $prv_was_coords_feature ||  # previous line was coords line with a feature (rare)
           $prv_was_coords_only)    {  # previous line was coords line without a feature (rarer)
          # line order makes sense, keep going...
          $cap_feature = $feature;
          $cap_feature =~ tr/a-z/A-Z/;
          $fac = $faccn . $fac_sep . $start . $fac_sep . $stop;

          # update '$prv_*' values that we use to make sure line order makes sense
          $prv_was_accn           = 0;
          $prv_was_coords_feature = 1;
          $prv_was_coords_only    = 0;
          $prv_was_quals          = 0;
          #printf("set prv_was_coords_feature\n");
        }
        else { # line order is unexpected
          DNAORG_FAIL("ERROR in $sub_name, unexpected line order (coords_feature) at line $line_ctr: $line", 1, $FH_HR);
        }
      }
      # -------------------------------------------------------
      elsif($line =~ /^(\<?\d+)\s+(\>?\d+)$/) { 
        # COORDINATES LINE WITHOUT A FEATURE NAME
        # example:
        #1	54
        ($start, $stop) = ($1, $2);

        # does line order make sense?
        if($prv_was_coords_feature || # previous line was a coords line with a feature (common)
           $prv_was_coords_only) {    # previous line was a coords line without a feature (common)
          # line order makes sense, keep going...

          $fac .= $fac_sep . $start . $fac_sep . $stop;
          
          # update '$prv_*' values that we use to make sure line order makes sense
          $prv_was_accn           = 0;
          $prv_was_coords_feature = 0;
          $prv_was_coords_only    = 1;
          $prv_was_quals          = 0;
          #printf("set prv_was_coords_only\n");
        }
        else { # line order is unexpected
          DNAORG_FAIL("ERROR in $sub_name, unexpected line order (coords_only line) at line $line_ctr: $line", 1, $FH_HR);
        }
      }
      # -------------------------------------------------------
      elsif(($line =~ /^\s+\S+\s+.+$/) || 
            ($line =~ /^\s+\S+$/)) { 
        # QUALIFIER LINE
        # examples:
        #			gene	AR1
        #			locus_tag	PhyvvsAgp1

        # before parsing it, do two sanity checks
        if(! defined $fac)     { DNAORG_FAIL("ERROR in $sub_name, coordinates undefined at a qualifier line", 1, $FH_HR); }
        if(! defined $feature) { DNAORG_FAIL("ERROR didn't read feature line before line: $line", 1, $FH_HR); }
        # and determine if we even care about this feature
        if((! $do_only_feature) || ($cap_feature eq $cap_only_feature)) { 
          # does line order make sense?
          if($prv_was_coords_feature || 
             $prv_was_coords_only    ||
             $prv_was_quals) { 
            # line order makes sense, keep going...
            if(! $prv_was_quals) { 
              # first quals line for this feature
              # at this point, we know that we have the full coordinates for the feature
              # so initialize the information
              if(! exists $fac_HHAR->{$feature}) {
                %{$fac_HHAR->{$feature}} = (); 
              }
              if(! exists $fac_HHAR->{$feature}{$faccn}) { 
                @{$fac_HHAR->{$feature}{$faccn}} = ();
              }
              push(@{$fac_HHAR->{$feature}{$faccn}}, $fac);
              if(! exists $quals_HHAR->{$feature}) { 
                %{$quals_HHAR->{$feature}} = ();
                @{$quals_HHAR->{$feature}{$fac}} = ();
              }
            } # end of 'if(! $prv_was_quals)'

            # now parse the line;
            # examples:
            #			gene	AR1
            #			locus_tag	PhyvvsAgp1
            if($line =~ /^\s+(\S+)\s+(.+)$/) { 
              ($qname, $qval) = ($1, $2);
            }
            elsif($line =~ /^\s+(\S+)$/) { 
              $qname = $1;
              $qval = "<no_value>";
            }
            else { 
              DNAORG_FAIL("ERROR in $sub_name, didn't parse quals line on second pass: $line", 1, $FH_HR); 
            }
            if($qname =~ m/\Q$qnqv_sep/)   { DNAORG_FAIL("ERROR in $sub_name, qualifier_name $qname has the string $qnqv_sep in it", 1, $FH_HR); }
            if($qval  =~ m/\Q$qnqv_sep/)   { DNAORG_FAIL("ERROR in $sub_name, qualifier_value $qval has the string $qnqv_sep in it", 1, $FH_HR); }
            my $qnqv = $qname . $qnqv_sep . $qval; # this is how we store the qualifier name and value, as a concatenated string in quals_HHAR
            push(@{$quals_HHAR->{$feature}{$fac}}, $qnqv);

            # and update the column data structures which just keep info on names and order of columns
            if(! exists $column_HH{$feature}) { 
              %{$column_HH{$feature}} = ();
              @{$column_HAR->{$feature}} = (); 
            }
            if(! exists $column_HH{$feature}{$qname}) { 
              push(@{$column_HAR->{$feature}}, $qname);
              $column_HH{$feature}{$qname} = 1;
            }
          }
          else { # unexpected line order
            DNAORG_FAIL("ERROR in $sub_name, unexpected line order (quals line) at line $line_ctr: $line", 1, $FH_HR);
          }          
        } # end of 'if((! $do_only_feature) || ($cap_feature ne $cap_only_feature))'
        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 0;
        $prv_was_coords_feature = 0;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 1;
        #printf("set prv_was_quals\n");
      }
      # -------------------------------------------------------
      else { 
        DNAORG_FAIL("ERROR in $sub_name, unable to parse line $line_ctr: $line", 1, $FH_HR); 
      }
      # -------------------------------------------------------
    }
  }
}

#################################################################
# Subroutine: parse_edirect_mat_peptide_file()
# Incept:     EPN, Thu Feb 11 14:01:30 2016
#
# Purpose:   Given a mature peptide file output from a
#            'efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product'
#            command, parse that file into usable data structures.
#
#            Can be called by the wrapper subroutine: 
#            edirectFtableOrMatPept2SingleFeatureTableInfo().
#
#            Caller will commonly call get_single_feature_table_info() after calling this
#            subroutine.
#
#            Similar to (and shares some code with) parse_edirect_ftable_file().
#
# Args:
#   $infile:          name of ftable file
#   $dummy_column:    name for dummy columns that we won't output, can be undef
#   $sep_HR:          ref to hash of 'separation' values, keys: "qnqv", "fac", and "qval", 
#                     used to split values that have been concatenated with these 'separation'
#                     values in between, PRE-FILLED
#   $quals_HHAR:      ref to 2D hash of arrays with values that could eventually go 
#                     into a single-feature table, FILLED HERE
#                     key 1: feature name, e.g. "CDS"
#                     key 2: 'fac' string: <full_accession><fac_seq><coordinates>
#                     value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>    
#   $faccn_AR:        REF to array of full accessions read, FILLED HERE
#   $fac_HHAR:        REF to 2D hash of arrays that is used to easily determine
#                     list of 2D keys ('fac's) in quals_HHA for a given feature and faccn
#                     FILLED HERE
#                     key 1: feature name: e.g. "CDS"
#                     key 2: 'faccn', full accession
#                     value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
#   $faccn2accn_HR:   REF to hash, key: full accession, value: short accession, 
#                     used for convenience in output function, FILLED HERE
#   $column_HAR:      REF to hash of arrays, key: feature, e.g. "CDS"
#                     value: array of qualifiers that exist for that feature in at least 1 faccn
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
#
# Returns:  void
#
# Dies:     if <$infile> file is in unexpected format
#################################################################
sub parse_edirect_mat_peptide_file {
  my $sub_name = "parse_edirect_mat_peptide_file()";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $dummy_column, $sep_HR, $quals_HHAR, $faccn_AR, $fac_HHAR, $faccn2accn_HR, $column_HAR, $FH_HR) = @_;

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  # example lines of a .mat_peptide file
  # NC_009942.1	97..465	anchored capsid protein C	
  # NC_009942.1	97..411	capsid protein C	
  # NC_009942.1	466..966	membrane glycoprotein precursor M	
  # NC_009942.1	join(2470..3552,3552..3680)	nonstructural protein NS1'

  my $faccn     = undef;   # full accession, e.g. ref|NC_000883.2|
  my $fac       = undef;   # full accession and coords, accession and each segment coords separated by $fac_sep
  my $coords    = undef;   # coordinates
  my $faccn2    = undef;   # another full accession
  my %column_HH = ();      # existence 2D hash, used to aid construction of %column_HA
                           # key 1: feature, e.g. "CDS"
                           # key 2: qualifier
                           # value: '1' (always)

  my $fac_sep  = $sep_HR->{"fac"};
  my $qnqv_sep = $sep_HR->{"qnqv"};
  my $qval_sep = $sep_HR->{"qval"};

  my $feature = "mat_peptide";
  while(my $line = <IN>) { 
    if($line =~ m/\w/) { 
# example line:
# NC_001477.1	95..436	anchored capsid protein C
      if($line =~ /(\S+)\s+(\S+)\s*(.*)$/) { 
        my ($acc, $coords, $product) = ($1, $2, $3); 
        # with .mat_peptide files, the accession is not in the 'full accession' ($faccn, e.g. ref|NC_001477.1|) 
        # that appears in the .ftable files, but we still want to use similar downstream code to 
        # handle both the parsed .mat_peptide and .ftable files, so we still fill faccn2accn_HR
        # for .mat_peptide files, even though the keys and values are identical

        $product =~ s/\s+$//; # remove trailing whitespace
        $faccn = $acc; # in this case, full accession and accession are the same
        # we use both $faccn and $acc so that this code matches the code in this subroutine matches
        # with the code in parse_edirect_ftable_file() for consistency.

        if(! exists $faccn2accn_HR->{$acc}) { 
          push(@{$faccn_AR}, $acc);
          $faccn2accn_HR->{$acc} = $acc;
        }

        # now we have to break down coords just so we can store each segment separately in fac_HHA, 
        # because the feature table has each segment separated, and we want to process and OUTPUT
        # both feature tables and mat_peptide information the same way after we parse it here. At
        # that point we'll put the coordinates together again. 
        my @starts_A = ();
        my @stops_A  = ();
        my $nsegments = 0;
        startsStopsFromCoords($coords, \@starts_A, \@stops_A, \$nsegments, $FH_HR);
        
        # printf("coords:    $coords\n");
        # printf("nsegments: $nsegments\n");
        for(my $i = 0; $i < $nsegments; $i++) { 
          my $start = $starts_A[$i];
          my $stop  = $stops_A[$i];
          # printf("starts_A[$i]: $starts_A[$i]\n");
          # printf("stops_A[$i]:  $stops_A[$i]\n");
          if($i == 0) { 
            $fac = $faccn . $fac_sep . $start . $fac_sep . $stop;
          }
          else { 
            $fac .= $fac_sep . $start . $fac_sep . $stop;
          }
        }

        if(! exists $fac_HHAR->{"mat_peptide"}) {
          %{$fac_HHAR->{"mat_peptide"}} = (); 
        }
        push(@{$fac_HHAR->{"mat_peptide"}{$acc}}, $fac);
        
        # first add the 'dummy' qual
        # we need to do this so that we can store the coordinate information (in $fac)
        # for mature peptides that have no additional information besides coordinates
        # (.mat_peptide files often include only mature peptide coordinates and no
        # product information ($product below) 
        my $qname = $dummy_column;
        my $qval  = "<no_value>";
        if($qname =~ m/\Q$qnqv_sep/)   { DNAORG_FAIL("ERROR in $sub_name, qualifier_name $qname has the string $qnqv_sep in it", 1, $FH_HR); }
        if($qval  =~ m/\Q$qnqv_sep/)   { DNAORG_FAIL("ERROR qualifier_value $qval has the string $qnqv_sep in it", 1, $FH_HR); }
        my $qnqv = $qname . $qnqv_sep . $qval; # this is how we store the qualifier name and value, as a concatenated string in quals_HHAR
        push(@{$quals_HHAR->{"mat_peptide"}{$fac}}, $qnqv);
        
        # and update the column data structures which just keep info on names and order of columns
        if(! exists $column_HH{"mat_peptide"}) { 
          %{$column_HH{"mat_peptide"}} = ();
          @{$column_HAR->{"mat_peptide"}} = (); 
        }
        if(! exists $column_HH{"mat_peptide"}{$qname}) { 
          push(@{$column_HAR->{"mat_peptide"}}, $qname);
          $column_HH{"mat_peptide"}{$qname} = 1;
        }
        
        if(defined $product && $product ne "") { 
          # now if the product qualifier has a value, add that too
          $qname = "product";
          $qval  = $product;
          $qnqv = $qname . $qnqv_sep . $qval;
          push(@{$quals_HHAR->{"mat_peptide"}{$fac}}, $qnqv);
          
          if(! exists $column_HH{"mat_peptide"}{$qname}) { 
            push(@{$column_HAR->{"mat_peptide"}}, $qname);
            $column_HH{"mat_peptide"}{$qname} = 1;
          }
        }
      }
    } # end of 'if($line =~ m/\w/)
  } # end of '$line = <IN>' 
}

#################################################################
# Subroutine: get_single_feature_table_info()
# Incept:     EPN, Fri Feb 12 09:51:48 2016
#
# Purpose:  Given data structures collected from
#           parseEdirecFtableFile() or parseEdirectMatPeptFile(), fill
#           a hash of hash of arrays (%{$tbl_HAR}) with the
#           information for a specific feature ($feature_name,
#           e.g. 'CDS' or mat_peptide').
#
# Arguments:
#   $ref_accn:       reference accession (without version)
#   $feature_name:   name of feature we want information for (e.g. CDS or mat_peptide)
#   $dummy_column:   name for dummy columns that we won't output, can be undef
#   $sep_HR:         ref to hash of 'separation' values, keys: "qnqv", "fac", and "qval"
#                    used to split values that have been concatenated with these 'separation'
#                    values in between, PRE-FILLED
#   $quals_HAR:      ref to hash of arrays with values we are printing, PRE-FILLED
#                    key 1: 'fac' string: <full_accession><fac_seq><coordinates>
#                    value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>    
#   $faccn_AR:       REF to array of full accessions read, PRE-FILLED
#   $fac_HAR:        REF to hash of arrays that is used to easily determine
#                    list of keys ('fac's) in quals_HA for a given feature and faccn
#                    key 1: 'faccn', full accession
#                    value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
#                    PRE-FILLED
#   $faccn2accn_HR:  REF to hash, key: full accession, value: short accession, 
#                    used for convenience in output function, PRE-FILLED
#   $column_AR:      REF to array of qualifiers for feature we're printing table for, PRE-FILLED
#   $tbl_HAR:        REF to hash of hash of arrays we'll fill with info on $qual_name:
#                    1D: key: qualifier (e.g. 'coords')
#                    2D: values for each qualifier, size will be number of features for this accession
#                    e.g. size of 5 means this accession has 5 CDS if $feature is CDS
#                    FILLED HERE
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if an passed in data structure is empty or contains invalid information
#################################################################
sub get_single_feature_table_info { 
  my $sub_name = "get_single_feature_table_info()";
  my $nargs_expected = 10;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($ref_accn, $dummy_column, $sep_HR, $quals_HAR, $faccn_AR, $fac_HAR, $faccn2accn_HR, $column_AR, $tbl_HAR, $FH_HR) = @_;

  my $start           = undef;   # start position
  my $stop            = undef;   # stop position
  my $fac             = undef;   # full accession and coords, accession and each segment coords separated by $fac_sep
  my $faccn           = undef;   # full accession, e.g. ref|NC_000883.2|
  my $faccn2          = undef;   # another full accession
  my $coords          = undef;   # coordinates
  my $qname           = undef;   # a qualifier name,  e.g. 'organism'
  my $qval            = undef;   # a qualifier value, e.g. 'Paramecia bursaria Chlorella virus 1'

  my $column;
  my $strand; 
  my $sort_coord;

  my $fac_sep  = $sep_HR->{"fac"};
  my $qnqv_sep = $sep_HR->{"qnqv"};
  my $qval_sep = $sep_HR->{"qval"};

  if(! %{$quals_HAR}) { DNAORG_FAIL("ERROR in $sub_name, quals hash of arrays is not defined", 1, $FH_HR); }
  if(! %{$fac_HAR})   { DNAORG_FAIL("ERROR in $sub_name, fac hash of arrays is not defined", 1, $FH_HR); }
  if(! @{$column_AR}) { DNAORG_FAIL("ERROR in $sub_name, column hash of arrays is not defined", 1, $FH_HR); }

  # go through all full accessions in order they were read from feature table
  foreach $faccn (@{$faccn_AR}) { 
    if(exists $fac_HAR->{$faccn}) { # if this accession has >= 1 qualifiers for this feature
      my $accn = $faccn2accn_HR->{$faccn};
      stripVersion(\$accn);
      if($accn eq $ref_accn) { 
        foreach $fac (@{$fac_HAR->{$faccn}}) { # foreach 'fac', accession + set of coords
          
          ($faccn2, $coords, $strand) = helper_breakdown_fac($fac, $fac_sep, $FH_HR);
          if($faccn ne $faccn2) { DNAORG_FAIL("ERROR in $sub_name, inconsistent fac value: $faccn ne $faccn2", 1, $FH_HR); }
          
          if(exists $quals_HAR->{$fac}) { # if there's any qualifiers for this fac
            # printf("quals_HA feature: fac: $fac exists!\n"); 
            
            push(@{$tbl_HAR->{"coords"}}, $coords);
            
            # for all columns in the table
            foreach $column (@{$column_AR}) {
              if((! defined $dummy_column) || ($column ne $dummy_column)) { 
                my $save_str = ""; 
                
                # for all qualifier names and values 
                foreach my $qnqv (@{$quals_HAR->{$fac}}) { 
                  ($qname, $qval) = split($qnqv_sep, $qnqv);
                  ### printf("faccn: $faccn qnqv: $qnqv split into $qname $qval\n");
                  
                  # if this qname matches this column, then it's the appropriate value to save here
                  if($qname eq $column) { 
                    if($save_str eq "") { # first value in this cell
                      $save_str = $qval;  
                    }
                    else { 
                      if($qval =~ m/\Q$qval_sep/) { DNAORG_FAIL("ERROR in $sub_name, qualifier_name $qval has the string $qval_sep in it", 1, $FH_HR); }
                      $save_str .= $qval_sep . $qval; # not first value, concatenate onto previous values
                    }
                  }
                }
                # do not save values that are '-', as far as I can tell these are just empty placeholders,
                # as an example, see the first CDS in the feature table returned by this query (as of 02/13/18):
                # esearch -db nuccore -query NC_031324 | efetch -format ft 
                # It has a qualifier value of '-' for the 'gene' qualifier, which only occurs because the
                # GenBank flat file does not list a gene name (gene qualifier) for the first gene (all 
                # other genes have gene qualifiers and so their corresponding CDS features do not have 
                # 'gene' qualifiers.
                if($save_str eq "-") { $save_str = ""; } 
                push(@{$tbl_HAR->{$column}}, $save_str);
              }
            } 
          }
        }
      }
    }
  }
  
  return;
}
#################################################################
# Subroutine: helper_breakdown_fac()
# Incept      EPN, Thu Feb 11 14:14:23 2016
#
# Purpose:    Breakdown a 'fac' string into it's parts and 
#             create a string in NCBI coordinate format from it.
#             A 'fac' string has accessions and coordinates in it.
#           
#             A 'helper' function called by 'get_single_feature_table_info()'.
#
# Arguments:
#   $fac:     full accession and coordinates, concatenated together
#             can be multiple coordinates e.g. "NC_000001:100:300:350:400"
#             means two segments: 100..300 and 350..400.
#   $fac_sep: character in between each token in $fac (":" in above example)
#   $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:
#   Four values:
#   $faccn:       full accession
#   $coords:      coordinates in format 'start1-stop1,start2-stop2,...startn-stopn'
#   $strand:      '+' if all segments are on fwd strand
#                 '-' if all segments are on neg strand
#                 '?' if all segments are 1 nucleotide (strand is consequently uncertain)
#                 '!' if >= 1 segment on two or more of following: fwd strand, rev strand, uncertain
#
# Dies:       if there's not an odd number of tokens
#################################################################
sub helper_breakdown_fac {
  my $sub_name = "helper_breakdown_fac()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($fac, $fac_sep, $FH_HR) = @_;

  my @elA = split($fac_sep, $fac);
  my $nel = scalar(@elA);
  if($nel % 2 == 0) { DNAORG_FAIL("ERROR in $sub_name, unexpectedly there's an even number of tokens: fac: $fac", 1, $FH_HR); }

  my $faccn = $elA[0];
  my ($istart, $istop) = (1, 2);
  my $start  = undef;   # start position
  my $stop   = undef;   # stop position
  my $npos   = 0; # number of segments that are on forward strand
  my $nneg   = 0; # number of segments that are on reverse strand
  my $nunc   = 0; # number of segments for which strand is uncertain (start == stop)
  my $nseg   = 0;
  my $coords = "";
  while($istop <= $nel) { 
    my ($orig_start, $orig_stop) = ($elA[$istart], $elA[$istop]);
    $start = $orig_start;
    $stop  = $orig_stop;
    $start =~ s/^<//;
    $stop  =~ s/^>//;
    # a more complicated version of this function existed as dnaorg.pm::helperBreakdownFac in dnaorg_scripts v0.45 and earlier 
    # that could handle more complicated start and stop coordinates including start coords that end with '^'
    if($start !~ /^\d+$/) { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse start coordinate $orig_start", 1, $FH_HR); 
    }
    if($stop !~ /^\d+$/) { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse stop coordinate $orig_stop", 1, $FH_HR); 
    }
    if($start == $stop) { # special case, we can't tell if we're in reverse complement or not
      $nunc++;
    }
    elsif($start > $stop) {
      $nneg++;
    }
    else { # $start < $stop, not reverse complement
      $npos++;
    }
    if($coords ne "") { $coords .= ","; }
    $coords .= $start . "-" . $stop;
    $nseg++;
    $istart += 2; 
    $istop  += 2;
  }

  # determine strand
  my $ret_strand = undef;
  if   ($npos == $nseg)          { $ret_strand = "+"; }
  elsif($nneg == $nseg)          { $ret_strand = "-"; }
  elsif($nunc  > 0)              { $ret_strand = "?"; }
  elsif($npos  > 0 && $nneg > 0) { $ret_strand = "!"; }
  else                           { DNAORG_FAIL("ERROR in $sub_name, unable to determine strand for fac: $fac", 1, $FH_HR); }

  return($faccn, $coords, $ret_strand);
}
#################################################################

#################################################################
# Subroutine: initialize_feature_info_hash_from_ftable_data()
# Incept:     EPN, Tue Feb 16 14:05:51 2016
# 
# Purpose:    Fill "type", "coords" and other keys in $qual_str
#             arrays in the feature info hash of arrays (%{$ftr_info_HAR})
#             using $feat_tbl_HHAR:
#
#             "type":    feature type, e.g. "mat_peptide", "CDS"
#             "coords":  coordinates for this feature in the reference
#             "strands": strands for these features in the reference
#              And whatever else is in the comma separated string $qual_str
# 
# Arguments:
#   $feat_tbl_HHAR:  ref to feature information, PRE-FILLED
#   $ftr_info_HAR:   ref to hash of arrays with feature information, FILLED HERE
#   $qual_str:       comma separated string of additional qualifiers to add
#                    e.g. "product,gene,exception"
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       Never
#
#################################################################
sub initialize_feature_info_hash_from_ftable_data { 
  my $sub_name = "initialize_feature_info_hash_from_ftable_data";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($feat_tbl_HHAR, $ftr_info_HAR, $qual_str, $FH_HR) = @_;
  
  my $nftr  = 0; # number of features added so far
  my $ntype = 0; # number of features of a type 
  my $type;      # a feature type
  my $qual;      # a qualifier
  my $i;         # counter
  
  # create array of qualifiers we will add
  my @qual_A = ();
  if((defined $qual_str) && ($qual_str ne "")) { 
    @qual_A = split(",", $qual_str);
  }
  
  # create order of feature names
  my @type_order_A = ("mat_peptide", "CDS", "gene");
  foreach $type (sort keys %{$feat_tbl_HHAR}) { 
    if(($type ne "mat_peptide") && ($type ne "CDS") && ($type ne "gene")) { 
      push(@type_order_A, $type);
    }
  }

  # set "coords", "strands" and any other values from input %
  foreach $type (@type_order_A) { 
    if(exists $feat_tbl_HHAR->{$type}) { 
      foreach $qual ("coords", @qual_A) { 
        if($qual eq "coords") { 
          $ntype = scalar(@{$feat_tbl_HHAR->{$type}{$qual}});
          if($ntype == 0) {
            DNAORG_FAIL("ERROR in $sub_name, for $type no coords values exist", 1, $FH_HR); 
          }
          # fill strands
          foreach my $coords_str (@{$feat_tbl_HHAR->{$type}{$qual}}) { 
            my @seg_strand_A = ();
            strand_array_from_fetched_coords_str($coords_str, \@seg_strand_A, $FH_HR);
            
          }
        }
        if(exists $feat_tbl_HHAR->{$type}{$qual}) { 
          push(@{$ftr_info_HAR->{$qual}}, @{$feat_tbl_HHAR->{$type}{$qual}});
        }
        else { # doesn't exist, fill with "" values
          for($i = 0; $i < $ntype; $i++) { 
            push(@{$ftr_info_HAR->{$qual}}, "");
          }
        }
      }
      # finally, the type (e.g. "CDS")
      for($i = 0; $i < $ntype; $i++) { 
        push(@{$ftr_info_HAR->{"type"}}, $type);
      }
      $nftr += $ntype;
    }
  }

  return;
}

#################################################################
# Subroutine: output_model_info_file()
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
#  $out_file:     out file to create
#  $name:         model name
#  $ftr_info_AHR: REF to array of hashes with information on the features, pre-filed
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $ftr_info_HAR is not valid upon entering
#################################################################
sub output_model_info_file { 
  my $sub_name = "output_model_info_file";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_file, $name, $ftr_info_AHR, $FH_HR) = @_;

  # ftr_info_AHR should already have array data for keys "coords" and "type"
  my @reqd_keys_A  = ("type", "coords");
  my $nftr = arrayOfHashesValidate($ftr_info_AHR, \@reqd_keys_A, $FH_HR);


  # create order of keys
  my ($key, $value2print);
  my @key_order_A  = ("type", "coords");
  my %key_ignore_H = ();
  $key_ignore_H{"type"}        = 1; # already added this to @key_order_A, so it goes first
  $key_ignore_H{"coords"}      = 1; # already added this to @key_order_A, so it goes second
  $key_ignore_H{"source_idx"}  = 1; # will be inferred from coords and type
  $key_ignore_H{"parent_idx"}  = 1; # will be inferred from coords and type
  $key_ignore_H{"3pa_ftr_idx"} = 1; # will be inferred from coords and type
  $key_ignore_H{"5p_seg_idx"}  = 1; # will be inferred from coords, when seg_info_HA is created
  $key_ignore_H{"3p_seg_idx"}  = 1; # will be inferred from coords, when seg_info_HA is created
  $key_ignore_H{"location"}    = 1; # *could* (but won't be) inferred from coords

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    foreach $key (sort keys %{$ftr_info_AHR->[$ftr_idx]}) { 
      if(! exists $key_ignore_H{$key}) { 
        push(@key_order_A, $key);
        $key_ignore_H{$key} = 1; 
      }
    }
  }

  # open output file
  open(OUT, ">", $out_file) || fileOpenFailure($out_file, $sub_name, $!, "writing", $FH_HR);
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    printf OUT ("%s", $name);
    foreach $key (@key_order_A) { 
      $value2print = (defined $ftr_info_AHR->[$ftr_idx]{$key}) ? $ftr_info_AHR->[$ftr_idx]{$key} : "";
      if($value2print eq "") { $value2print = "-"; }
      printf OUT ("\t%s:%s", $key, $value2print);
    }
    print OUT "\n";
  }
  close(OUT);

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


  
#################################################################
# Subroutine: edirectFetchToFile()
# Incept:     EPN, Tue Mar 12 12:18:37 2019
#
# Synopsis: Fetch information for an accession using edirect.
#
# Arguments:
#  $out_file:  output file to create
#  $accn:      accession to fetch
#  $format:    format to fetch (e.g. "gpc", "ft", "fasta")
#  $nattempts: number of times to retry 
#  $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's a problem fetching the data
#################################################################
sub edirectFetchToFile { 
  my $sub_name = "edirectFetchToFile";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_file, $accn, $format, $nattempts, $FH_HR) = @_;
  if((! defined $nattempts) || ($nattempts < 1)) { $nattempts = 1; }

  my $url = edirectFetchUrl($accn, $format);

  my $n = 0;
  my $fetched_str = undef;
  while(($n < $nattempts) && (! defined $fetched_str)) { 
    $fetched_str = get($url);
    $n++;
    sleep(1);
  }
  if(! defined $fetched_str) { 
    DNAORG_FAIL("ERROR in $sub_name, problem fetching $accn (undefined)", 1, $FH_HR); 
  }

  open(OUT, ">", $out_file) || fileOpenFailure($out_file, $sub_name, $!, "writing", $FH_HR);
  print OUT $fetched_str;
  close(OUT);

  return;
}

#################################################################
# Subroutine: edirectFetchUrl()
# Incept:     EPN, Tue Mar 12 12:18:37 2019
#
# Synopsis: Return a url for an efetch command
#
# Arguments:
#  $accn:      accession to fetch
#  $format:    format to fetch (e.g. "gpc", "ft", "fasta")
#
# Returns:    void
#
# Dies:       never
#################################################################
sub edirectFetchUrl { 
  my $sub_name = "edirectFetchUrl";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($accn, $format) = @_;

  return sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=%s&retmode=text", $accn, $format);
}

#################################################################
# Subroutine: genbankParse()
# Incept:     EPN, Tue Mar 12 14:04:14 2019
#
# Synopsis: Parse a GenBank format file.
#
# Arguments:
#  $infile:   GenBank file to parse
#  $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if we have trouble parsing the file
#
# Reference: https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
#            http://www.insdc.org/documents/feature-table
#################################################################
sub genbankParse { 
  my $sub_name = "genbankParse";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $seq_info_HHR, $ftr_info_HAHR, $FH_HR) = @_;

  my $line_idx  = 0;     # line index of input file
  my $acc       = undef; # accession, read from LOCUS line
  my $tmp_acc   = undef; # accession, read from ACCESSION or VERSION line
  my $len       = undef; # length, read from LOCUS line
  my $def       = undef; # seq definition, read from DEFINITION line
  my $ver       = undef; # sequence version, read from VERSION line
  my $feature   = undef; # a feature   read from a feature/location line in the FEATURES section
  my $location  = undef; # a location  read from a feature/location line in the FEATURES section
  my $qualifier = undef; # a qualifier read from a qualifier/value  line in the FEATURES section
  my $value     = undef; # a value     read from a qualifier/value  line in the FEATURES section
  my $seq       = undef; # sequence, read from the ORIGIN section
  my $seqline   = undef; # single line of sequence
  my $seq_idx   = 0;     # number of sequences read
  my $ftr_idx   = -1;    # number of features read for current sequence
  my $line      = undef; # a line

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  $line = <IN>; 
  while(defined $line) { 
    chomp $line; $line_idx++;
    if($line =~ /^LOCUS\s+(\S+)\s+(\d+)/) { 
      #LOCUS       NC_039477               7567 bp    RNA     linear   VRL 22-FEB-2019
      if((defined $acc) || (defined $len)) { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple LOCUS lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      ($acc, $len) = ($1, $2);
      # initialize the array of hashes for this accession's features
      if(defined $ftr_info_HAHR->{$acc}) { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, trying to add feature info for accession $acc, but it already exists, line:\n$line\n", 1, $FH_HR);
      }
      @{$ftr_info_HAHR->{$acc}} = ();
      $line = <IN>; 
    }
    elsif($line =~ /^DEFINITION\s+(.*)$/) { 
      #DEFINITION  Norovirus GII isolate strain Hu/GBR/2016/GII.P16-GII.4_Sydney/226,
      #            complete genome.
      if(defined $def) { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple DEFINITION lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      $def = $1;
      # read remainder of the definition (>= 0 lines)
      $line = <IN>; 
      while((defined $line) && ($line =~ /^\s+(.+)$/)) {
        chomp $line; $line_idx++;
        $def .= $1;
        $line = <IN>; 
      }
    }
    elsif($line =~ /^ACCESSION\s+(\S+)$/) { 
      # ACCESSION   NC_039477
      # verify this matches what we read in the LOCUS line
      $tmp_acc = $1;
      if((! defined $acc) || ($tmp_acc ne $acc)) { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, accession mismatch for $tmp_acc, line:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>;
    }
    elsif($line =~ /^VERSION\s+(\S+)$/) { 
      #VERSION     NC_039477.1
      if(defined $ver) { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple VERSION lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      # verify this matches what we read in the LOCUS line
      $ver = $1;
      $tmp_acc = $ver;
      stripVersion(\$tmp_acc);
      if((! defined $acc) || ($tmp_acc ne $acc)) { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, version/accession mismatch for $tmp_acc, line:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>;
    }
    elsif($line =~ /^FEATURES\s+Location\/Qualifiers$/) { 
      # parse the features and then the sequence
      # FEATURES section
      # two types of line:
      # feature/location line
      #        example:      gene            5..5104
      #        example:      misc_feature    join(2682..2689,1..2)
      #        example:      misc_feature    join(161990..162784,complement(88222..88806),complement(86666..87448))
      # qualifier/value line type A, first line of a new qualifier
      #        example: /codon_start=1
      #        example: /gene="ORF1"
      # qualifier/value line type B, not the first line of a new qualifier, line 2 to N of a qualifier value
      #        example: QNVIDPWIRNNFVQAPGGEFTVSPRNAPGEILWSAPLGPDLNPYLSHLARMYNGYAGG
      #        example: IPPNGYFRFDSWVNQFYTLAPMGNGTGRRRVV"
      if($ftr_idx != -1) { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple FEATURES lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>;
      while((defined $line) && ($line !~ /^ORIGIN/)) { 
        chomp $line; $line_idx++;
        if($line =~ /^\s+\/(\S+)\=(.+)$/) { # first token must start with '/'
          # qualifier/value line type A, examples:
          #  /codon_start=1
          #  /gene="ORF1"
          #  /translation="MKMASNDATVAVACNNNNDKEKSSGEGLFTNMSSTLKKALGARP
          my ($save_qualifier, $save_value) = ($1, $2);
          if(defined $value) { # we are finished with previous value
            genbank_store_qualifier_value(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
          }
          ($qualifier, $value) = ($save_qualifier, $save_value);
        }
        elsif($line =~ /^\s+(\S+)\s+(\S+)$/) { 
          # NOTE: this will pass for a non-first line of a qualifier value that has whitespace in it:
          # e.g.                      KQP ASRDESQKPPRPPTPELVKRIPPPPPNGEEEEEPVIRYEVKSGISGLPELTTVPQ
          # But I think those are illegal, if they're not, then we'll set "KQP" as feature below, which is bad
          if(defined $value) { # we are finished with previous value
            genbank_store_qualifier_value(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
            ($qualifier, $value) = (undef, undef);
          }
          # feature/location line, examples:
          #   gene            5..5104
          ($feature, $location) = ($1, $2);
          $ftr_idx++;
          genbank_store_qualifier_value(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",     $feature,  $FH_HR);
          genbank_store_qualifier_value(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "location", $location, $FH_HR);
        }
        else { 
          # qualifier/value line type B
          #        example: QNVIDPWIRNNFVQAPGGEFTVSPRNAPGEILWSAPLGPDLNPYLSHLARMYNGYAGG
          #        example: IPPNGYFRFDSWVNQFYTLAPMGNGTGRRRVV"
          $line =~ s/^\s+//; # remove leading whitespace
          if(! defined $value) { 
            DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, in FEATURES section read qualifier value line without qualifier first, line:\n$line\n", 1, $FH_HR);
          }
          $value .= $line; 
        }
        $line = <IN>; chomp $line; $line_idx++;
      }
      if(! defined $line) { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, expected to read ORIGIN line after FEATURES but did not\n", 1, $FH_HR);
      }
      # if we get here we just read the ORIGIN line
      # first store final qualifier/value
      if(defined $value) { 
        genbank_store_qualifier_value($ftr_info_HAHR->{$acc}, $ftr_idx, $qualifier, $value, $FH_HR);
      }
      # parse the ORIGIN sequence
      $line = <IN>;
      # sanity check
      if(defined $seq) { 
        DNAORG_FAIL("ERROR in $sub_name, read multiple ORIGIN lines for single record ($acc), line:\n$line\n", 1, $FH_HR);
      }
      $seq = "";
      while((defined $line) && ($line !~ /^\/\/$/)) { 
        chomp $line; $line_idx++;
        # sequence lines
        # examples:
        # 7501 gtcacgggcg taatgtgaaa agacaaaact gattatcttt ctttttcttt agtgtctttt
        # 7561 aaaaaaa
        if($line =~ /^\s+\d+\s+(.+)$/) { 
          $seqline = $1;
          $seqline =~ s/\s+//g; # remove spaces
          $seq .= $seqline;
        }
        $line = <IN>;
      }
      if(! defined $line) { 
        DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, expected to find a // line after ORIGIN but did not, line $line_idx\n", 1, $FH_HR);
      }
      # if we get here we just read the // line
      # we are finished with this sequence, store the information
      if(! defined $acc) { DNAORG_FAIL(        "ERROR in $sub_name, failed to read accession, line: $line_idx\n", 1, $FH_HR); }
      if(! defined $len) { DNAORG_FAIL(sprintf("ERROR in $sub_name, failed to read length (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), 1, $FH_HR); }
      if(! defined $ver) { DNAORG_FAIL(sprintf("ERROR in $sub_name, failed to read version (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), 1, $FH_HR); }
      if(! defined $def) { DNAORG_FAIL(sprintf("ERROR in $sub_name, failed to read definition (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), 1, $FH_HR); }
      if(! defined $seq) { DNAORG_FAIL(sprintf("ERROR in $sub_name, failed to read sequence (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), 1, $FH_HR); }

      # store sequence info
      %{$seq_info_HH{$acc}} = ();
      $seq_info_HH{$acc}{"len"} = $len;
      $seq_info_HH{$acc}{"ver"} = $ver;
      $seq_info_HH{$acc}{"def"} = $def;
      $seq_info_HH{$acc}{"seq"} = $seq;

      # reset variables
      $seq = undef;
      $len = undef;
      $acc = undef;
      $ver = undef;
      $def = undef;
      $feature   = undef;
      $location  = undef;
      $qualifier = undef;
      $value     = undef;

      $seq_idx++;
      $ftr_idx = -1;
      $line = <IN>;
    } # end of 'elsif($line =~ /^FEATURES\s+Location\/Qualifiers$/) {' 
    else { 
      # not a line we will parse, read the next line
      $line = <IN>;
    }
  }

  if($seq_idx == 0) { 
    DNAORG_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, failed to read any sequence data\n", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: genbank_store_qualifier_value()
# Incept:     EPN, Wed Mar 13 09:42:22 2019
#
# Synopsis: Store a genbank qualifier and value.
#
# Arguments:
#  $ftr_info_AHR: REF to the array of hashes to store data in
#  $ftr_idx:      feature index
#  $qualifier:    qualifier
#  $value:        qualifier value
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#
# Returns:    '1' if $ftr_info_AHR->[$ftr_idx]{$qualifier} created
#             '0' if $ftr_info_AHR->[$ftr_idx]{$qualifier} exists upon entering function
#
# Dies:       If $value includes the string ":GPSEP:, which we use 
#             to separate multiple qualifier values for the same qualifier.
#             
#################################################################
sub genbank_store_qualifier_value { 
  my $sub_name = "genbank_store_qualifier_value";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($ftr_info_AHR, $ftr_idx, $qualifier, $value, $FH_HR) = @_;

  if($value =~ /\:GPSEP\:/) { 
    DNAORG_FAIL("ERROR in $sub_name, qualifier value $value includes the special string :GPSEP:, this is not allowed", 1, $FH_HR);
  }

  # remove leading and trailing " in the value, if they exist
  # GenBank format uses "" as a substitute for " in these strings
  $value =~ s/^\"//;
  $value =~ s/\"$//;

  # printf("in $sub_name q: $qualifier v: $value\n");
  if(! defined ($ftr_info_AHR->[$ftr_idx])) { 
    %{$ftr_info_AHR->[$ftr_idx]} = (); 
  }
  if(! defined $ftr_info_AHR->[$ftr_idx]{$qualifier}) { 
    $ftr_info_AHR->[$ftr_idx]{$qualifier} = $value;
  }
  else { 
    $ftr_info_AHR->[$ftr_idx]{$qualifier} .= ":GBSEP:" . $value;
  }

  return;
}


