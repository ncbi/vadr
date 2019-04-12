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
require "epn-ofile.pm";
require "epn-seq.pm";
require "epn-seqfile.pm";
require "epn-utils.pm";

#######################################################################################
# What this script does: 
#
# - Preliminaries: 
#   o processes options
#   o creates the output directory
#   o outputs program banner
#   o makes sure the required executables are executable
#
# - Fetches GenBank file (if ! --gb)
# - Parses GenBank file
# - Prunes data read from GenBank file
# - Parses input stockholm file (if --stk)
# - Fills in feature and segment info
# - Translates CDS (if any) and creates BLAST db
# - Builds CM 
# - Presses CM 
# - Writes model info file
# - Writes optional output files
# 
#######################################################################################
# make sure required environment variables are set
my $env_dnaorg_scripts_dir  = utl_DirEnvVarValid("DNAORGSCRIPTSDIR");
my $env_dnaorg_blast_dir    = utl_DirEnvVarValid("DNAORGBLASTDIR");
my $env_dnaorg_infernal_dir = utl_DirEnvVarValid("DNAORGINFERNALDIR");
my $env_dnaorg_easel_dir    = utl_DirEnvVarValid("DNAORGEASELDIR");

my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmbuild"}       = $env_dnaorg_infernal_dir . "/cmbuild";
$execs_H{"cmpress"}       = $env_dnaorg_infernal_dir . "/cmpress";
$execs_H{"esl-reformat"}  = $env_dnaorg_easel_dir    . "/esl-reformat";
$execs_H{"esl-translate"} = $env_dnaorg_easel_dir    . "/esl-translate";
$execs_H{"makeblastdb"}   = $env_dnaorg_blast_dir    . "/makeblastdb";
utl_ExecHValidate(\%execs_H, undef);

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
my $g = 0; # option group

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
$opt_group_desc_H{++$g} = "basic options";
#     option            type       default  group   requires incompat     preamble-output                                                help-output    
opt_Add("-h",           "boolean", 0,           0,    undef, undef,       undef,                                                         "display this help",                                   \%opt_HH, \@opt_order_A);
opt_Add("-g",           "string", 0,           $g,    undef, undef,       "define model group for model info file as <s>",               "define model group for model info file as <s>", \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,          $g,    undef, undef,       "forcing directory overwrite",                                 "force; if dir <output directory> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,          $g,    undef, undef,       "be verbose",                                                  "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--stk",        "string",  undef,      $g,    undef,  undef,      "read single sequence stockholm 'alignment' from <s>",         "read single sequence stockholm 'alignment' from <s>", \%opt_HH, \@opt_order_A);
opt_Add("--gb",         "string",  undef,      $g,    undef,  undef,      "read genbank file from <s>, don't fetch it",                  "read genbank file from <s>, don't fetch it", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,          $g,    undef, undef,       "leave intermediate files on disk",                            "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling what feature types are stored in model info file\n[default set is: CDS,gene,mat_peptide]";
#     option            type       default  group   requires incompat     preamble-output                                                      help-output    
opt_Add("--fall",       "boolean", 0,          $g,    undef,  undef,      "store info for all feature types (except those in --fskip)",        "store info for all feature types (except those in --fskip)", \%opt_HH, \@opt_order_A);
opt_Add("--fadd",       "string",  undef,      $g,    undef,"--fall",     "also store features types in comma separated string <s>",           "also store feature types in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--fskip",      "string",  undef,      $g,    undef,  undef,      "do not store info for feature types in comma separated string <s>",  "do not store info for feature types in comma separated string <s>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling what qualifiers are stored in model info file\n[default set is:product,gene,exception]";
#     option            type       default  group   requires incompat     preamble-output                                                   help-output    
opt_Add("--qall",       "boolean",  0,        $g,    undef,  undef,       "store info for all qualifiers (except those in --qskip)",        "store info for all qualifiers (except those in --qskip)", \%opt_HH, \@opt_order_A);
opt_Add("--qadd",       "string",   undef,    $g,    undef,"--qall",      "also store info for qualifiers in comma separated string <s>",   "also store info for qualifiers in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--qskip",      "string",   undef,    $g,    undef,  undef,       "do not store info for qualifiers in comma separated string <s>", "do not store info for qualifiers in comma separated string <s>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling CDS translation step";
#     option          type       default    group   requires    incompat   preamble-output                                             help-output    
opt_Add("--ttbl",     "integer", 1,            $g,  undef,         undef,  "use NCBI translation table <n> to translate CDS",          "use NCBI translation table <n> to translate CDS", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling cmbuild step";
#     option          type       default    group   requires    incompat   preamble-output                                             help-output    
opt_Add("--cmn",      "integer", 0,           $g,   undef, "--skipbuild",  "set number of seqs for glocal fwd HMM calibration to <n>", "set number of seqs for glocal fwd HMM calibration to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--cmp7ml",   "boolean", 0,           $g,   undef, "--skipbuild",  "set CM's filter p7 HMM as the ML p7 HMM",                  "set CM's filter p7 HMM as the ML p7 HMM",                  \%opt_HH, \@opt_order_A);
opt_Add("--cmere",    "real",    0,           $g,   undef,  "--skipbuild", "set CM relative entropy target to <x>",                    "set CM relative entropy target to <x>",                    \%opt_HH, \@opt_order_A);
opt_Add("--cmeset",   "real",    0,           $g,   undef,  "--skipbuild", "set CM eff seq # for CM to <x>",                           "set CM eff seq # for CM to <x>",                           \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for skipping stages";
#       option             type       default     group requires   incompat  preamble-output                                    help-output    
opt_Add("--skipbuild",     "boolean", 0,         $g,    undef,     undef,    "skip the cmbuild step",                           "skip the cmbuild step", \%opt_HH, \@opt_order_A);
opt_Add("--onlyurl",       "boolean", 0,         $g,    undef,"--stk,--gb",  "output genbank file url for accession and exit",  "output genbank file url for accession and exit", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "optional output files";
#       option       type       default     group  requires     incompat  preamble-output                          help-output    
opt_Add("--ftrinfo",    "boolean", 0,         $g,    undef,     undef,    "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);
opt_Add("--sgminfo",    "boolean", 0,         $g,    undef,     undef,    "output internal segment information",   "create file with internal segment information", \%opt_HH, \@opt_order_A);


# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_build.pl [-options] <accession> <path to output directory to create>\n";
my $synopsis = "dnaorg_build.pl :: build homology model of a single sequence for feature annotation";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'g=s'          => \$GetOptions_H{"-g"},
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'stk=s'        => \$GetOptions_H{"--stk"},
                'gb=s'         => \$GetOptions_H{"--gb"},
                'keep'         => \$GetOptions_H{"--keep"},
# options for controlling what feature types are stored in model info file
                'fall'         => \$GetOptions_H{"--fall"},
                'fadd=s'       => \$GetOptions_H{"--fadd"},
                'fskip=s'      => \$GetOptions_H{"--fskip"},
# options for controlling what qualifiers are stored in model info file
                'qall'         => \$GetOptions_H{"--qall"},
                'qadd=s'       => \$GetOptions_H{"--qadd"},
                'qskip=s'      => \$GetOptions_H{"--qskip"},
# options for controlling CDS translation step
                'ttbl=s'       => \$GetOptions_H{"--ttbl"},
# options for controlling cmbuild step
                'cmn=s'        => \$GetOptions_H{"--cmn"},
                'cmp7ml'       => \$GetOptions_H{"--cmp7ml"},
                'cmere=s'      => \$GetOptions_H{"--cmere"},
                'cmeset=s'     => \$GetOptions_H{"--cmeset"},
# optional for skipping stages
                'skipbuild'    => \$GetOptions_H{"--skipbuild"},
                'onlyurl'      => \$GetOptions_H{"--onlyurl"},
# optional output files
                'sgminfo'      => \$GetOptions_H{"--sgminfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"});

my $total_seconds = -1 * ofile_SecondsSinceEpoch(); # by multiplying by -1, we can just add another ofile_SecondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.9";
my $releasedate   = "Apr 2019";
my $pkgname       = "dnaorg";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, undef);
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
my ($mdl_name, $dir) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

############################################
# if --onlyurl used, output the url and exit
############################################
if(opt_Get("--onlyurl", \%opt_HH)) { 
  print dng_EutilsFetchUrl($mdl_name, "gb") . "\n";
  exit 0;
}

#############################
# create the output directory
#############################
my $cmd;              # a command to run with utl_RunCommand()
my @early_cmd_A = (); # array of commands we run before our log file is opened

if($dir !~ m/\/$/) { $dir =~ s/\/$//; } # remove final '/' if it exists
if(-d $dir) { 
  $cmd = "rm -rf $dir";
  if(opt_Get("-f", \%opt_HH)) { utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR directory named $dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $dir) { 
  $cmd = "rm $dir";
  if(opt_Get("-f", \%opt_HH)) { utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR a file named $dir already exists. Remove it, or use -f to overwrite it."; }
}

# create the dir
$cmd = "mkdir $dir";
utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root = $dir . "/" . $dir_tail . ".dnaorg_build";

#######################
# output program banner
#######################
# output preamble
my @arg_desc_A = ("accession/model name", "output directory");
my @arg_A      = ($mdl_name, $dir);
my %extra_H    = ();
$extra_H{"\$DNAORGSCRIPTSDIR"}  = $env_dnaorg_scripts_dir;
$extra_H{"\$DNAORGINFERNALDIR"} = $env_dnaorg_infernal_dir;
$extra_H{"\$DNAORGEASELDIR"}    = $env_dnaorg_easel_dir;
$extra_H{"\$DNAORGBLASTDIR"}    = $env_dnaorg_blast_dir;
ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
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
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, $pkgname, "log",  $out_root . ".log",  1, "Output printed to screen");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, $pkgname, "cmd",  $out_root . ".cmd",  1, "List of executed commands");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, $pkgname, "list", $out_root . ".list", 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
my $FH_HR  = $ofile_info_HH{"FH"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# open optional output files
if(opt_Get("--ftrinfo", \%opt_HH)) { 
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, $pkgname, "ftrinfo", $out_root . ".ftrinfo", 1, "Feature information (created due to --ftrinfo)");
}
if(opt_Get("--sgminfo", \%opt_HH)) { 
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, $pkgname, "sgminfo", $out_root . ".sgminfo", 1, "Segment information (created due to --sgminfo)");
}

# now we have the log file open, output the banner there too
ofile_OutputBanner($log_FH, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

#############################################################
# make sure the required executables exist and are executable
#############################################################

###########################################
# Fetch the genbank file (if --gb not used)
###########################################
my $progress_w = 60; # the width of the left hand column in our progress output, hard-coded
my $start_secs;
my $gb_file = undef;
if(opt_IsUsed("--gb", \%opt_HH)) { 
  $gb_file = opt_Get("--gb", \%opt_HH);
}
else { 
  # --gb not used, create gb file by fetching using eutils
  $start_secs = ofile_OutputProgressPrior("Fetching GenBank file", $progress_w, $log_FH, *STDOUT);

  $gb_file = $out_root . ".gb";
  dng_EutilsFetchToFile($gb_file, $mdl_name, "gb", 5, $ofile_info_HH{"FH"});  # number of attempts to fetch to make before dying
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "gb", $gb_file, 1, "GenBank format file for $mdl_name");

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

########################
# Parse the genbank file
########################
$start_secs = ofile_OutputProgressPrior("Parsing GenBank file", $progress_w, $log_FH, *STDOUT);

my %ftr_info_HAH = (); # the feature info 
my %seq_info_HH  = (); # the sequence info 
sqf_GenbankParse($gb_file, \%seq_info_HH, \%ftr_info_HAH, $FH_HR);
if((! exists $seq_info_HH{$mdl_name}) || (! defined $seq_info_HH{$mdl_name}{"seq"})) { 
  ofile_FAIL("ERROR parsing GenBank file $gb_file, did not read sequence for reference accession $mdl_name\n", "dnaorg", 1, $FH_HR);
}
if(! exists $ftr_info_HAH{$mdl_name}) { 
  ofile_FAIL("ERROR parsing GenBank file $gb_file, did not read info for reference accession $mdl_name\n", "dnaorg", 1, $FH_HR);
}

ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#######################################################
# Prune data read from %ftr_info_HAH, only keeping what
# we want to output to the eventual model info file
#######################################################
$start_secs = ofile_OutputProgressPrior("Pruning data read from GenBank file", $progress_w, $log_FH, *STDOUT);

# determine what types of features we will store based on cmdline options
# --fall is incompatible with --fadd
my %fdf_H   = (); # default feature types to keep
my %fadd_H  = (); # feature types to add
my %fskip_H = (); # feature types to skip
process_add_and_skip_options("CDS,gene,mat_peptide", "--fadd", "--fskip", \%fdf_H, \%fadd_H, \%fskip_H, \%opt_HH, $FH_HR);

# determine what qualifiers we will store based on cmdline options
# --qall is incompatible with --qadd
my %qdf_H   = (); # default qualifiers to keep
my %qadd_H  = (); # qualifiers to add
my %qskip_H = (); # qualifiers to skip
process_add_and_skip_options("type,location,product,gene,exception", "--qadd", "--qskip", \%qdf_H, \%qadd_H, \%qskip_H, \%opt_HH, $FH_HR);

# remove all features types we don't want
my $ftr_idx;
my @ftr_idx_to_remove_A = ();
for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$mdl_name}}); $ftr_idx++) { 
  my $ftype = $ftr_info_HAH{$mdl_name}[$ftr_idx]{"type"};
  # we skip this type and remove it from ftr_info_HAH
  # if all three of A1, A2, A3 OR B is satisfied
  # (A1) it's not a default feature type AND
  # (A2) it's not listed in --fadd AND
  # (A3) --fall not used
  # OR 
  # (B) it is listed in --fskip string 
  if(((! defined $fdf_H{$ftype})     && # (A1)
      (! defined $fadd_H{$ftype})    && # (A2)
      (! opt_Get("--fall", \%opt_HH)))  # (A3)
     || (defined $fskip_H{$ftype})) {   # (B)
    splice(@{$ftr_info_HAH{$mdl_name}}, $ftr_idx, 1);
    $ftr_idx--; # this is about to be incremented
  }
}

# remove any qualifier key/value pairs with keys not in %qual_H, unless --qall used
for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$mdl_name}}); $ftr_idx++) { 
  foreach my $qual (sort keys %{$ftr_info_HAH{$mdl_name}[$ftr_idx]}) { 
    # we skip this qualifier and remove it from ftr_info_HAH
    # if all three of A1, A2, A3 OR B is satisfied
    # (A1) it's not a default qualifier AND
    # (A2) it's not listed in --qadd AND
    # (A3) --qall not used
    # OR 
    # (B) it is listed in --qskip string 
    if(((! defined $qdf_H{$qual})        && # (A1)
        (! defined $qadd_H{$qual})       && # (A2)
        (! opt_Get("--qall", \%opt_HH)))    # (A3)
       || (defined $qskip_H{$qual})) {      # (B)
      delete $ftr_info_HAH{$mdl_name}[$ftr_idx]{$qual};
    }
  }
}

ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#####################################################################
# Parse the input stockholm file (if --stk) or create it (if ! --stk)
#####################################################################
my $stk_file = $out_root . ".stk";
my $fa_file  = $out_root . ".fa";
my $stk_has_ss = undef;
my $in_stk_file = opt_Get("--stk", \%opt_HH);
if(defined $in_stk_file) { 
  $start_secs = ofile_OutputProgressPrior("Validating input Stockholm file", $progress_w, $log_FH, *STDOUT);

  $stk_has_ss = stockholm_validate_single_sequence_input($in_stk_file, $seq_info_HH{$mdl_name}{"seq"}, \%opt_HH, $FH_HR);

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  
  $start_secs = ofile_OutputProgressPrior("Reformatting Stockholm file to FASTA file", $progress_w, $log_FH, *STDOUT);

  utl_RunCommand("cp $in_stk_file $stk_file", opt_Get("-v", \%opt_HH), 0, $FH_HR);
  sqf_EslReformatRun($execs_H{"esl-reformat"}, $stk_file, $fa_file, "stockholm", "fasta", \%opt_HH, $FH_HR);

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}
else { 
  # --stk not used, we create it by first making a fasta file of the model
  # model sequence read from the gb file, then converting that fasta file 
  # to a stockholm file
  $start_secs = ofile_OutputProgressPrior("Creating FASTA sequence file", $progress_w, $log_FH, *STDOUT);

  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, $pkgname, "fasta", $fa_file, 1, "fasta sequence file for $mdl_name");
  sqf_FastaWriteSequence($ofile_info_HH{"FH"}{"fasta"}, 
                         $seq_info_HH{$mdl_name}{"ver"}, 
                         $seq_info_HH{$mdl_name}{"def"}, 
                         $seq_info_HH{$mdl_name}{"seq"}, $FH_HR);
  close $ofile_info_HH{"FH"}{"fasta"};

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  $start_secs = ofile_OutputProgressPrior("Reformatting FASTA file to Stockholm file", $progress_w, $log_FH, *STDOUT);

  sqf_EslReformatRun($execs_H{"esl-reformat"}, $fa_file, $stk_file, "afa", "stockholm", \%opt_HH, $FH_HR);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "stk", $stk_file, 1, "Stockholm alignment file for $mdl_name");

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

######################################################################
# Finish populating @{$ftr_info_HAH{$mdl_name} and create @sgm_info_AH
######################################################################
$start_secs = ofile_OutputProgressPrior("Finalizing feature information", $progress_w, $log_FH, *STDOUT);

dng_FeatureInfoImputeCoords(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
dng_FeatureInfoImputeLength(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
dng_FeatureInfoImputeSourceIdx(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
dng_FeatureInfoImputeParentIndices(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
dng_FeatureInfoImputeOutname(\@{$ftr_info_HAH{$mdl_name}});

my @sgm_info_AH = (); # segment info, inferred from feature info
dng_SegmentInfoPopulate(\@sgm_info_AH, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);

ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###################################
# Translate the CDS, if we have any
###################################
my $ncds = dng_FeatureInfoCountType(\@{$ftr_info_HAH{$mdl_name}}, "CDS");
my $cds_fa_file = undef;
my $protein_fa_file = undef;
if($ncds > 0) { 
  $start_secs = ofile_OutputProgressPrior("Translating CDS and building BLAST DB", $progress_w, $log_FH, *STDOUT);

  $cds_fa_file  = $out_root . ".cds.fa";
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, $pkgname, "cdsfasta", $cds_fa_file, 1, "fasta sequence file for CDS from $mdl_name");
  dng_CdsFetchStockholmToFasta($ofile_info_HH{"FH"}{"cdsfasta"}, $stk_file, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  close $ofile_info_HH{"FH"}{"cdsfasta"};
  
  $protein_fa_file = $out_root . ".protein.fa";
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, $pkgname, "proteinfasta", $protein_fa_file, 1, "fasta sequence file for translated CDS from $mdl_name");
  sqf_EslTranslateCdsToFastaFile($ofile_info_HH{"FH"}{"proteinfasta"}, $execs_H{"esl-translate"}, $cds_fa_file, 
                                 $out_root, \@{$ftr_info_HAH{$mdl_name}}, \%opt_HH, $FH_HR);
  close $ofile_info_HH{"FH"}{"proteinfasta"};

  sqf_BlastDbProteinCreate($execs_H{"makeblastdb"}, $protein_fa_file, \%opt_HH, $FH_HR);

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

##############
# Build the CM
##############
my $cm_file = undef;
if(! opt_Get("--skipbuild", \%opt_HH)) { 
  my $cmbuild_str = undef;
  my $clen_times_cmn = $seq_info_HH{$mdl_name}{"len"} * 200;
  if(opt_IsUsed("--cmn", \%opt_HH)) { 
    $clen_times_cmn *= (opt_Get("--cmn", \%opt_HH) / 200);
  }
  if   ($clen_times_cmn > 4000000) { $cmbuild_str = "(may take more than an hour)"; }
  elsif($clen_times_cmn > 3000000) { $cmbuild_str = "(should take roughly an hour)"; }
  elsif($clen_times_cmn > 2000000) { $cmbuild_str = "(should take roughly 20-40 minutes)"; }
  elsif($clen_times_cmn > 1000000) { $cmbuild_str = "(should take roughly 10-30 minutes)"; }
  elsif($clen_times_cmn >  500000) { $cmbuild_str = "(should take roughly 5-10 minutes)"; }
  else                             { $cmbuild_str = "(shouldn't take more than a few minutes)"; }

  $start_secs = ofile_OutputProgressPrior("Building model $cmbuild_str", $progress_w, $log_FH, *STDOUT);

  my $cmbuild_opts = "-n $mdl_name --verbose ";
  if((! defined $stk_has_ss) || (! $stk_has_ss)) { $cmbuild_opts .= " --noss"; }
  if(opt_IsUsed("--cmn",    \%opt_HH)) { $cmbuild_opts .= " --EgfN " . opt_Get("--cmn", \%opt_HH); }
  if(opt_IsUsed("--cmp7ml", \%opt_HH)) { $cmbuild_opts .= " --p7ml"; }
  if(opt_IsUsed("--cmere",  \%opt_HH)) { $cmbuild_opts .= " --ere "  . opt_Get("--cmere", \%opt_HH); }
  if(opt_IsUsed("--cmeset", \%opt_HH)) { $cmbuild_opts .= " --eset " . opt_Get("--cmeset", \%opt_HH); }

  my $cmbuild_file = $out_root . ".cmbuild";
  $cm_file         = $out_root . ".cm";
  my $cmbuild_cmd  = $execs_H{"cmbuild"} . " " . $cmbuild_opts . " $cm_file $stk_file > $cmbuild_file";
  utl_RunCommand($cmbuild_cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});
  ofile_OutputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);

  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "cm",      $cm_file, 1, "CM file");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "cmbuild", $cmbuild_file, 1, "cmbuild output file");

  # press the file we just created 
  $start_secs = ofile_OutputProgressPrior("Pressing CM file", $progress_w, $log_FH, *STDOUT);
  my $cmpress_file = $out_root . ".cmpress";
  my $cmpress_cmd  = $execs_H{"cmpress"} . " $cm_file > $cmpress_file";
  utl_RunCommand($cmpress_cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});
  ofile_OutputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);

  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "i1m",     $cm_file . ".i1m", 1, "binary CM and p7 HMM filter file");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "i1i",     $cm_file . ".i1i", 1, "SSI index for binary CM file");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "i1f",     $cm_file . ".i1f", 1, "optimized p7 HMM filters (MSV part)");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "i1p",     $cm_file . ".i1p", 1, "optimized p7 HMM filters (remainder)");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "cmpress", $cmpress_file,     1, "cmpress output file");
}

########################
# Output model info file
########################
$start_secs = ofile_OutputProgressPrior("Creating model info file", $progress_w, $log_FH, *STDOUT);

# create @mdl_info_AH, and add info for our lone model
# modelInfoFileWrite() can output data for multiple models at once, 
# but we use it here only for a single model.
my @mdl_info_AH = (); 
%{$mdl_info_AH[0]} = ();
$mdl_info_AH[0]{"name"}   = $mdl_name;
$mdl_info_AH[0]{"length"} = $seq_info_HH{$mdl_name}{"len"};
if(defined $cm_file) { 
  $mdl_info_AH[0]{"cmfile"} = utl_RemoveDirPath($cm_file);
}
if($ncds > 0) { 
  $mdl_info_AH[0]{"blastdb"} = utl_RemoveDirPath($protein_fa_file);
  if((opt_IsUsed("--ttbl", \%opt_HH)) && (opt_Get("--ttbl", \%opt_HH) != 1))  { 
    $mdl_info_AH[0]{"transl_table"} = opt_Get("--ttbl", \%opt_HH);
  }
}
if(opt_IsUsed("-g", \%opt_HH)) { 
  $mdl_info_AH[0]{"group"} = opt_Get("-g", \%opt_HH); 
}
my $modelinfo_file  = $out_root . ".modelinfo";
dng_ModelInfoFileWrite($modelinfo_file, \@mdl_info_AH, \%ftr_info_HAH, $FH_HR);
ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $pkgname, "modelinfo", $modelinfo_file, 1, "DNAORG 'model info' format file for $mdl_name");

ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

##########
# Conclude
##########
# output optional output files
if(exists $ofile_info_HH{"FH"}{"ftrinfo"}) { 
  utl_AHDump("Feature information", \@{$ftr_info_HAH{$mdl_name}}, $ofile_info_HH{"FH"}{"ftrinfo"});
}
if(exists $ofile_info_HH{"FH"}{"sgminfo"}) { 
  utl_AHDump("Segment information", \@sgm_info_AH, $ofile_info_HH{"FH"}{"sgminfo"});
}

$total_seconds += ofile_SecondsSinceEpoch();
ofile_OutputConclusionAndCloseFiles($total_seconds, $pkgname, $dir, \%ofile_info_HH);
exit 0;

###############
# SUBROUTINES #
###############

#################################################################
# Subroutine: stockholm_validate_single_sequence_input()
# Incept:     EPN, Fri Mar 15 13:19:32 2019
#
# Synopsis: Validate an input Stockholm file is in the correct 
#           format, has exactly 1 sequence and no gaps.
#
# Arguments:
#  $in_stk_file:  input stockholm file to validate
#  $exp_sqstring: sequence we expect to be in the stockholm alignment
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    '1' if Stockholm file has SS_cons annotation, else '0'
#
# Dies:       if there's a problem parsing the file or 
#             a requirement is not met
#################################################################
sub stockholm_validate_single_sequence_input {
  my $sub_name = "stockholm_validate_single_sequence_input";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($in_stk_file, $exp_sqstring, $opt_HHR, $FH_HR) = @_;

  if(! -e $in_stk_file) { ofile_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file does not exist", "dnaorg", 1, $FH_HR); }
  if(! -s $in_stk_file) { ofile_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file exists but is empty", "dnaorg", 1, $FH_HR); }
  if(  -d $in_stk_file) { ofile_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file is actually a directory", "dnaorg", 1, $FH_HR); }
  my $msa = Bio::Easel::MSA->new({ fileLocation => $in_stk_file, isDna => 1});
  my $nseq = $msa->nseq;
  if($nseq == 1) { 
    # single sequence, make sure there are no gaps
    if($msa->any_allgap_columns) { 
      ofile_FAIL("ERROR, read 1 sequence in --stk file $in_stk_file, but it has gaps, this is not allowed for single sequence 'alignments' (remove gaps with 'esl-reformat --mingap')", "dnaorg", 1, $FH_HR);
    }
    # validate it matches $exp_sqstring
    my $fetched_sqstring = $msa->get_sqstring_unaligned(0);
    seq_SqstringCapitalize(\$fetched_sqstring);
    seq_SqstringCapitalize(\$exp_sqstring);
    seq_SqstringDnaize(\$fetched_sqstring);
    seq_SqstringDnaize(\$exp_sqstring);
    if($fetched_sqstring ne $exp_sqstring) { 
      my $summary_sqstring_diff_str = seq_SqstringDiffSummary($fetched_sqstring, $exp_sqstring);
      ofile_FAIL("ERROR, read 1 sequence in --stk file $in_stk_file, but it does not match sequence read from GenBank file $gb_file:\n$summary_sqstring_diff_str", "dnaorg", 1, $FH_HR); 
    }
  }
  else { # nseq != 1
    ofile_FAIL("ERROR, did not read exactly 1 sequence in --stk file $in_stk_file.\nTo use DNAORG with models built from alignments of multiple sequences,\nyou will have to build the CM with cmbuild and create the model info file manually.\n", "dnaorg", 1, $FH_HR);
  }

  return $msa->has_ss_cons;
}

#################################################################
# Subroutine: process_add_and_skip_options()
# Incept:     EPN, Mon Mar 18 06:29:21 2019
#
# Synopsis: Process cmdline --{f,q}add and --{f,q}skip options 
#           for features or qualifiers.
#
# Arguments:
#  $df_string:  comma separated string of default values (e.g. "CDS,gene,mat_peptide")
#  $add_opt:    name of add option (e.g. "--fadd")
#  $skip_opt:   name of skip option (e.g. "--fskip")
#  $df_HR:      ref to hash of default keys, filled here, values will all be '1'
#  $add_HR:     ref to hash of keys to add, filled here, values will all be '1'
#  $skip_HR:    ref to hash of keys to skip, filled here, values will all be '1'
#  $opt_HHR:    ref to hash of cmdline options
#  $FH_HR:      ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if a key is listed in both the $add_opt and $skip_opt option strings.
#################################################################
sub process_add_and_skip_options { 
  my $sub_name = "process_add_and_skip_options";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($df_string, $add_opt, $skip_opt, $df_HR, $add_HR, $skip_HR, $opt_HHR, $FH_HR) = @_;

  utl_ExistsHFromCommaSepString($df_HR, $df_string);
  if(opt_IsUsed($add_opt,  $opt_HHR)) { utl_ExistsHFromCommaSepString($add_HR,  opt_Get($add_opt,  $opt_HHR)); }
  if(opt_IsUsed($skip_opt, $opt_HHR)) { utl_ExistsHFromCommaSepString($skip_HR, opt_Get($skip_opt, $opt_HHR)); }
  # make sure $add_opt and $skip_opt have no values in common
  foreach my $key (sort keys (%{$add_HR})) { 
    if(defined $skip_HR->{$key}) { 
      ofile_FAIL("ERROR in $sub_name, processing $add_opt <s1> and $skip_opt <s2> options, $key exists in both <s1> and <s2>", "dnaorg", 1, $FH_HR);
    }
  }

  return;
}


