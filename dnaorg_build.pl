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
opt_Add("--mdlinfo",    "boolean", 0,                        2,    undef, undef, "output internal model information",     "create file with internal model information",   \%opt_HH, \@opt_order_A);
opt_Add("--ftrinfo",    "boolean", 0,                        2,    undef, undef, "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"3"} = "options for skipping stages and using files from an earlier, identical run, primarily useful for debugging";
#     option               type       default               group   requires    incompat                  preamble-output                                            help-output    
opt_Add("--skipedirect",   "boolean", 0,                       3,   undef,      undef,                    "skip the edirect steps, use existing results",           "skip the edirect steps, use data from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipfetch",     "boolean", 0,                       3,   undef,      undef,                    "skip the sequence fetching steps, use existing results", "skip the sequence fetching steps, use files from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipbuild",     "boolean", 0,                       3,   undef,      undef,                    "skip the build/calibrate steps",                         "skip the model building/calibrating, requires --mdlinfo and/or --ftrinfo", \%opt_HH, \@opt_order_A);

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
                'mdlinfo'      => \$GetOptions_H{"--mdlinfo"},
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
#$execs_H{"esl_fetch_cds"} = $esl_fetch_cds;
$execs_H{"makeblastdb"}   = $blast_exec_dir . "makeblastdb";
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################################################################
# Step 1. Output the consopts file that dnaorg_annotate.pl will use to 
#         make sure options used are consistent between dnaorg_build.pl and 
#         dnaorg_annotate.pl 
###########################################################################
my $progress_w = 80; # the width of the left hand column in our progress output, hard-coded
my $start_secs;
#my $start_secs = outputProgressPrior("Outputting information on options used for future use with dnaorg_annotate.pl", $progress_w, $log_FH, *STDOUT);
#output_consopts_file($out_root . ".consopts", \%opt_HH, \%ofile_info_HH);
#outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###########################################################################
# Step 2. Gather and process information on reference genome using Edirect.
###########################################################################
$start_secs = outputProgressPrior("Gathering information on reference using edirect", $progress_w, $log_FH, *STDOUT);

my %feat_tbl_HHHA = ();  # feature data from feature table file, hash of hash of hashes of arrays
                         # 1D: qualifier name, e.g. 'CDS'
                         # 2D: key: accession
                         # 3D: key: column name in gene ftable file
                         # 4D: per-row values for each column

my $other_str = "gene";
fetch_reference_info_edirect($ref_accn, $other_str, $out_root, \%feat_tbl_HHHA,
                             \%ofile_info_HH, \%opt_HH, $ofile_info_HH{"FH"}); 

#
#if($do_matpept) {  
#  # validate the CDS:mat_peptide relationships that we read from the $matpept input file
#  matpeptValidateCdsRelationships(\@cds2pmatpept_AA, \%{$cds_tbl_HHA{$ref_accn}}, \%{$mp_tbl_HHA{$ref_accn}}, 0, $seq_info_HA{"len"}[0], $ofile_info_HH{"FH"});
#}
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

exit 0;

#########################################################
# Step 3. Fetch and process the reference genome sequence
##########################################################
$start_secs = outputProgressPrior("Fetching and processing the reference genome", $progress_w, $log_FH, *STDOUT);
my %mdl_info_HA = ();          # hash of arrays, values are arrays [0..$nmdl-1];
                               # see dnaorg.pm::validateModelInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAllSequencesAndProcessReferenceSequence()
my %ftr_info_HA = ();          # hash of arrays, values are arrays [0..$nftr-1], 
                               # see dnaorg.pm::validateFeatureInfoHashIsComplete() for list of all keys
                               # filled in wrapperFetchAllSequencesAndProcessReferenceSequence()

# convert the input stockholm alignment into a fasta file
#my $fa_file = $out_root . ".fa";
#stockholm_to_fasta($in_stk_file, $out_root . ".fa", $ofile_info_HH{"FH"});

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
# Step 3B. Fetch the CDS protein translations and build BLAST database
#######################################################################
$start_secs = outputProgressPrior("Fetching protein translations of CDS and building BLAST DB", $progress_w, $log_FH, *STDOUT);
my @prot_fa_file_A = ();
fetch_proteins_into_fasta_files($out_root, $ref_accn, \%ftr_info_HA, \@prot_fa_file_A, \%opt_HH, \%ofile_info_HH);

foreach my $prot_fa_file (@prot_fa_file_A) { 
  create_blast_protein_db(\%execs_H, $prot_fa_file, \%opt_HH, \%ofile_info_HH);
}
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################
# Step 4. Build the model
#########################
if(! opt_Get("--skipbuild", \%opt_HH)) { 
  $start_secs = outputProgressPrior("Building model (this could take a while)", $progress_w, $log_FH, *STDOUT);

  my $cmbuild_cmd  = $execs_H{"cmbuild"} . " --verbose -F $out_root.cm " . $ofile_info_HH{"fullpath"}{"refstk"} . " > $out_root.cmbuild";
  #runCommand($cmbuild_cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});

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
#   $ref_accn:              reference accession, first accession in $listfile (although this is 
#                           not enforced here, caller enforced it)
#   $other_str:             string of 'other' features to parse, e.g. "gene,RNA"
#   $out_root:              string that is the 'root' for naming output files
#   $feat_tbl_HHHAR:        REF to hash of hash of hash of arrays for other (non-CDS and non-MP) feature info, FILLED HERE
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

  my ($ref_accn, $other_str, $out_root, $feat_tbl_HHHAR, $ofile_info_HHR, $opt_HHR, $FH_HR) = @_;

  my $cmd; 
  my @other_A = ();
  if((defined $other_str) && ($other_str ne "")) { 
    @other_A = split(",", $other_str);
  }  

  # mature peptides
  my $mp_file = $out_root . ".mat_peptide";
  $cmd = "esearch -db nuccore -query $ref_accn | efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product > $mp_file";
  runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  if(-s $mp_file) { 
    addClosedFileToOutputInfo($ofile_info_HHR, "mp", $mp_file, 0, "Mature peptide information obtained via edirect");
    edirect_ftable_or_matpept_to_single_feature_table_info($mp_file, 1, "mat_peptide", $feat_tbl_HHHAR->{"mat_peptide"}, $FH_HR); # 1: it is a mat_peptide file
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

  # CDS from feature table
  edirect_ftable_or_matpept_to_single_feature_table_info($ft_file, 0, "CDS", $feat_tbl_HHHAR->{"CDS"}, $FH_HR); # 0: it's not a mat_peptide file

  # other features from feature table
  foreach my $ofeat (@other_A) { 
    edirect_ftable_or_matpept_to_single_feature_table_info($ft_file, 0, $ofeat, \%{$feat_tbl_HHHAR->{$ofeat}}, $FH_HR); # 0: it's not a mat_peptide file
  }    

  return 0;
}

#################################################################
# Subroutine: edirect_ftable_or_matpept_to_single_feature_table_info()
# Incept:     EPN, Fri Feb 12 09:48:00 2016
# 
# Purpose:    Given an edirect output feature table file or mature
#             peptide table file, parse it and store it's relevant
#             information in a 'single feature table' (%{$tbl_HHAR}).
#
#             This is a wrapper subroutine that calls 
#             parse_edirect_mat_peptide_file() or 
#             parse_edirect_ftable_file() and then 
#             get_single_feature_table_info().
#
# Arguments:
#   $edirect_file:  name of edirect output file to parse
#   $do_matpept:    '1' if edirect file is a mature peptide file, '0' or undef if it is a ftable file
#   $feature:       the feature we want to store info on in $tbl_HHAR (e.g. "CDS" or "mat_peptide")
#   $tbl_HHAR:      REF to hash of hash of arrays we'll fill with info on $qual_name:
#                   1D: key: accession
#                   2D: key: qualifier (e.g. 'coords')
#                   3D: values for each qualifier, size will be number of features for this accession
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
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($edirect_file, $do_matpept, $feature, $tbl_HHAR, $FH_HR) = @_;

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

  %{$tbl_HHAR} = ();

  if(defined $do_matpept && $do_matpept) { 
    parse_edirect_mat_peptide_file($edirect_file, $dummy_column, \%sep_H, \%quals_HHA, \@faccn_A, \%fac_HHA, \%faccn2accn_H, \%column_HA, $FH_HR);
  }
  else { # default, a ftable file
    parse_edirect_ftable_file($edirect_file, undef, \%sep_H, \%quals_HHA, \@faccn_A, \%fac_HHA, \%faccn2accn_H, \%column_HA, $FH_HR);
  }
  printf("CHECKING $feature\n");
  if(exists $quals_HHA{$feature}) { 
    get_single_feature_table_info($dummy_column, \%sep_H, \%{$quals_HHA{$feature}}, \@faccn_A, \%{$fac_HHA{$feature}}, \%faccn2accn_H, \@{$column_HA{$feature}}, $tbl_HHAR, $FH_HR);
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
#           a hash of hash of arrays (%{$tbl_HHAR}) with the
#           information for a specific feature ($feature_name,
#           e.g. 'CDS' or mat_peptide').
#
# Arguments:
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
#   $tbl_HHAR:       REF to hash of hash of arrays we'll fill with info on $qual_name:
#                    1D: key: accession
#                    2D: key: qualifier (e.g. 'coords')
#                    3D: values for each qualifier, size will be number of features for this accession
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
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($dummy_column, $sep_HR, $quals_HAR, $faccn_AR, $fac_HAR, $faccn2accn_HR, $column_AR, $tbl_HHAR, $FH_HR) = @_;

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
      if(%{$tbl_HHAR} && exists $tbl_HHAR->{$accn}) { 
        DNAORG_FAIL("ERROR in $sub_name, read data twice in $sub_name for accession $accn", 1, $FH_HR);
      }
      foreach $fac (@{$fac_HAR->{$faccn}}) { # foreach 'fac', accession + set of coords

        ($faccn2, $coords, $sort_coord, $strand) = helper_breakdown_fac($fac, $fac_sep, $FH_HR);
        if($faccn ne $faccn2) { DNAORG_FAIL("ERROR in $sub_name, inconsistent fac value: $faccn ne $faccn2", 1, $FH_HR); }
        
        if(exists $quals_HAR->{$fac}) { # if there's any qualifiers for this fac
          # printf("quals_HA feature: fac: $fac exists!\n"); 

          push(@{$tbl_HHAR->{$accn}{"coords"}},         $coords);
          push(@{$tbl_HHAR->{$accn}{"strand"}},         $strand);
          push(@{$tbl_HHAR->{$accn}{"min-coord"}},      $sort_coord);
          
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
              # old behavior:
              ## if there's no value for this qualifier, put '<empty>'
              ##if($save_str eq "") { $save_str = "<empty>"; }

              # do not save values that are '-', as far as I can tell these are just empty placeholders,
              # as an example, see the first CDS in the feature table returned by this query (as of 02/13/18):
              # esearch -db nuccore -query NC_031324 | efetch -format ft 
              # It has a qualifier value of '-' for the 'gene' qualifier, which only occurs because the
              # GenBank flat file does not list a gene name (gene qualifier) for the first gene (all 
              # other genes have gene qualifiers and so their corresponding CDS features do not have 
              # 'gene' qualifiers.
              if($save_str eq "-") { $save_str = ""; } 
              push(@{$tbl_HHAR->{$accn}{$column}}, $save_str);
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
#             A 'helper' function called by 'getSingleFeatureTableInfo()'.
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
#   $ncbi_coords: NCBI format for coordinates
#   $sort_coord:  minimum coordinate of all segments, possibly useful for ordering multiple features
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
  my $is_rev = 0; # true if current segment is on reverse strand
  my $nfwd   = 0; # number of segments that are on forward strand
  my $nrev   = 0; # number of segments that are on reverse strand
  my $nunc   = 0; # number of segments for which strand is uncertain (start == stop)
  my $some_revcomp  = undef;
  my @ncbi_coords_A = (); # [0..$nsegments-1] strings to concatenate to make ncbi coords
  my @is_rev_A      = (); # [0..$nsegments-1] '1' if current segment is on reverse strand
  my $nsegments     = 0;
  my $have_carrot   = 0; # set to '1' if $start has a '^' at the end, e.g. 804^ 805
  my $sort_coord = undef; # start position of first segment
  my $min_coord  = undef;
  while($istop <= $nel) { 
    my ($orig_start, $orig_stop) = ($elA[$istart], $elA[$istop]);
    $start = $orig_start;
    $stop  = $orig_stop;
    if($start =~ s/\^$//) { 
      $have_carrot = 1;
    }
    $start =~ s/^<//;
    $start =~ s/^>//;
    $stop  =~ s/^<//;
    $stop  =~ s/^>//;
    $min_coord = ($start < $stop) ? $start : $stop;
    if((! defined $sort_coord) || ($min_coord < $sort_coord)) { $sort_coord = $min_coord; }

    if($have_carrot) { 
      if(abs($start-$stop) > 1) { DNAORG_FAIL("ERROR found carrot but range is not exactly 1 nt: fac: $fac", 1, $FH_HR); }
    }

    if($start == $stop) { # special case, we can't tell if we're in reverse complement or not
      $nunc++;
      push(@is_rev_A, 0);
    }
    elsif($start > $stop) {
      $nrev++;
      $is_rev = 1;
      push(@is_rev_A, 1);
    }
    else { # $start < $stop, not reverse complement
      $nfwd++;
      push(@is_rev_A, 0);
    }
    if($have_carrot) { 
      push(@ncbi_coords_A, $orig_start . $orig_stop);
    }
    else { 
      if($is_rev) { push(@ncbi_coords_A, $orig_stop .  ".." . $orig_start); }
      else        { push(@ncbi_coords_A, $orig_start . ".." . $orig_stop); }
    }
    $nsegments++;
    $istart += 2; 
    $istop  += 2;
  }
  if($have_carrot) { 
    if($nsegments > 1) { DNAORG_FAIL("ERROR found carrot but more than one segment: fac: $fac", 1, $FH_HR); }
  }

  my $ncbi_coords = "";
  if($nfwd > 0 && $nrev > 0) { 
    # special case, we need to put 'complement(' ')' around each rev segment separately
    # e.g: join(161990..162784,complement(88222..88806),complement(86666..87448))
    for(my $i = 0; $i < $nsegments; $i++) { 
      if($i > 0) { $ncbi_coords .= ","; }
      if($is_rev_A[$i]) { 
        # first, deal with a non-obvious situation, where in the feature table
        # '>' and '<' characters indicating incompleteness are inverted relative
        # to how they are in the actual annotation. 
        # NC_007030.2 complement(4370..>4576)
        # is in the feature table as: <4576	4370	CDS
        $ncbi_coords_A[$i] =~ s/\</\!/;
        $ncbi_coords_A[$i] =~ s/\>/\</;
        $ncbi_coords_A[$i] =~ s/\!/\>/;
        $ncbi_coords .= "complement(" . $ncbi_coords_A[$i] . ")" 
      }
      else { 
        $ncbi_coords .= $ncbi_coords_A[$i]; 
      }
    }
  }
  else { # normal case, all exons/segments are on the same strand
    # if we're on the reverse strand, we need to reverse the order of the 
    # segments, because the order of reverse strand segments in a feature table is 
    # opposite what it is in Entrez, and our other scripts use Entrez
    # format, so we enforce that convention here.
    if($nrev > 0) { 
      for(my $i = $nsegments-1; $i >= 0; $i--) { 
        if($i < ($nsegments-1)) { $ncbi_coords .= ","; }
        $ncbi_coords .= $ncbi_coords_A[$i];
      }
    }
    else { # positive strand
      for(my $i = 0; $i < $nsegments; $i++) { 
        if($i > 0) { $ncbi_coords .= ","; }
        $ncbi_coords .= $ncbi_coords_A[$i];
      }
    }
  }
  if($nsegments > 1) { # more than one segment
    $ncbi_coords = "join(" . $ncbi_coords . ")";
  }
  # now add complement for cases where are segments/segments are on reverse strand
  # impt to do this after the join, so we get complement(join()) instead of
  # join(complement())
  if($nfwd == 0 && $nrev > 0) { # all segments are on reverse strand
    # first, deal with a non-obvious situation, where in the feature table
    # '>' and '<' characters indicating incompleteness are inverted relative
    # to how they are in the actual annotation. 
    # NC_007030.2 complement(4370..>4576)
    # is in the feature table as: <4576	4370	CDS
    $ncbi_coords =~ s/\</\!/g;
    $ncbi_coords =~ s/\>/\</g;
    $ncbi_coords =~ s/\!/\>/g;
    # now add the 'complement()'
    $ncbi_coords = "complement(" . $ncbi_coords . ")";
  }

  # printf("in $sub_name input: $fac, returning $faccn $coords ncbi_coords:$ncbi_coords\n");

  # determine strand
  my $ret_strand = undef;
  if   ($nfwd == $nsegments)     { $ret_strand = "+"; }
  elsif($nrev == $nsegments)     { $ret_strand = "-"; }
  elsif($nunc  > 0)              { $ret_strand = "?"; }
  elsif($nfwd  > 0 && $nrev > 0) { $ret_strand = "!"; }
  else                           { DNAORG_FAIL("ERROR in $sub_name, unable to determine strand for fac: $fac", 1, $FH_HR); }

  return($faccn, $ncbi_coords, $sort_coord, $ret_strand);
}
#################################################################
