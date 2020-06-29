#!/usr/bin/env perl
# EPN, Wed May  1 10:48:51 2019 [renamed to vadr-build.pl]
# EPN, Mon Feb  1 15:07:43 2016 [dnaorg_build.pl split off from dnaorg_annotate_genomes.pl]
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
#
use strict;
use warnings;
use Getopt::Long qw(:config no_auto_abbrev);
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;
use LWP::Simple; 
use LWP::Protocol::https; 

require "vadr.pm";
require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seq.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

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
my $env_vadr_scripts_dir  = utl_DirEnvVarValid("VADRSCRIPTSDIR");
my $env_vadr_blast_dir    = utl_DirEnvVarValid("VADRBLASTDIR");
my $env_vadr_infernal_dir = utl_DirEnvVarValid("VADRINFERNALDIR");
my $env_vadr_hmmer_dir    = utl_DirEnvVarValid("VADRHMMERDIR");
my $env_vadr_easel_dir    = utl_DirEnvVarValid("VADREASELDIR");

# make sure the required executables exist and are executable
my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmbuild"}       = $env_vadr_infernal_dir . "/cmbuild";
$execs_H{"cmemit"}        = $env_vadr_infernal_dir . "/cmemit";
$execs_H{"cmfetch"}       = $env_vadr_infernal_dir . "/cmfetch";
$execs_H{"cmpress"}       = $env_vadr_infernal_dir . "/cmpress";
$execs_H{"hmmbuild"}      = $env_vadr_hmmer_dir    . "/hmmbuild";
$execs_H{"hmmpress"}      = $env_vadr_hmmer_dir    . "/hmmpress";
$execs_H{"esl-reformat"}  = $env_vadr_easel_dir    . "/esl-reformat";
$execs_H{"esl-sfetch"}    = $env_vadr_easel_dir    . "/esl-sfetch";
$execs_H{"esl-translate"} = $env_vadr_easel_dir    . "/esl-translate";
$execs_H{"makeblastdb"}   = $env_vadr_blast_dir    . "/makeblastdb";
utl_ExecHValidate(\%execs_H, undef);

#########################################################
# Command line and option processing using sqp_opts.pm
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
opt_Add("-f",           "boolean", 0,          $g,    undef, undef,       "forcing directory overwrite",                                 "force; if dir <output directory> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,          $g,    undef, undef,       "be verbose",                                                  "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--stk",        "string",  undef,      $g,    undef, undef,       "read single sequence stockholm 'alignment' from <s>",         "read single sequence stockholm 'alignment' from <s>", \%opt_HH, \@opt_order_A);
opt_Add("--infa",       "string",  undef,      $g,    undef, undef,       "read single sequence fasta file from <s>, don't fetch it",    "read single sequence fasta file from <s>, don't fetch it", \%opt_HH, \@opt_order_A);
opt_Add("--inft",       "string",  undef,      $g, "--inft", "--gb",      "read feature table file from <s>, don't fetch it",            "read feature table file from <s>, don't fetch it", \%opt_HH, \@opt_order_A);
opt_Add("--ftfetch1",   "boolean", 0,          $g,    undef, "--inft,--gb,--ftfetch2", "fetch feature table with efetch -format ft",      "fetch feature table with efetch -format ft", \%opt_HH, \@opt_order_A);
opt_Add("--ftfetch2",   "boolean", 0,          $g,    undef, "--inft,--gb,--ftfetch1", "fetch feature table with efetch -format gbc | xml2tbl", "fetch feature table with efetch -format gbc | xml2tbl", \%opt_HH, \@opt_order_A);
opt_Add("--gb",         "boolean", 0,          $g,    undef, undef,       "parse a genbank file, not a feature table file",              "parse a genbank file, not a feature table file", \%opt_HH, \@opt_order_A);
opt_Add("--ingb",       "string",  undef,      $g,   "--gb", undef,       "read genbank file from <s>, don't fetch it",                  "read genbank file from <s>, don't fetch it", \%opt_HH, \@opt_order_A);
opt_Add("--addminfo",   "string",  undef,      $g,    undef, undef,       "add feature info from model info file <s>",                   "add feature info from model info file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--forcelong",  "boolean", 0,          $g,    undef, undef,       "allow long models > 25Kb in length",                          "allow long models > 25Kb in length", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,          $g,    undef, undef,       "leave intermediate files on disk",                            "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling what feature types are stored in model info file\n[default set is: CDS,gene,mat_peptide]";
#     option            type       default  group   requires incompat     preamble-output                                                      help-output    
opt_Add("--fall",       "boolean", 0,          $g,    undef,  undef,      "store info for all feature types (except those in --fskip)",        "store info for all feature types (except those in --fskip)", \%opt_HH, \@opt_order_A);
opt_Add("--fadd",       "string",  undef,      $g,    undef,"--fall",     "also store features types in comma separated string <s>",           "also store feature types in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--fskip",      "string",  undef,      $g,    undef,  undef,      "do not store info for feature types in comma separated string <s>",  "do not store info for feature types in comma separated string <s>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling what qualifiers are stored in model info file\n[default set is:product,gene,exception]";
#     option            type       default  group   requires incompat     preamble-output                                                             help-output    
opt_Add("--qall",       "boolean",  0,        $g,    undef,  undef,       "store info for all qualifiers (except those in --qskip)",                  "store info for all qualifiers (except those in --qskip)", \%opt_HH, \@opt_order_A);
opt_Add("--qadd",       "string",   undef,    $g,    undef,"--qall",      "also store info for qualifiers in comma separated string <s>",             "also store info for qualifiers in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--qftradd",    "string",   undef,    $g,"--qadd",    undef,      "--qadd <s2> only applies for feature types in comma separated string <s>", "--qadd <s2> only applies for feature types in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--qskip",      "string",   undef,    $g,    undef,  undef,       "do not store info for qualifiers in comma separated string <s>",           "do not store info for qualifiers in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--noaddgene", "boolean",  0,        $g,    undef,  undef,       "do not add gene qualifiers from gene features to overlapping features",     "do not add gene qualifiers from gene features to overlapping features", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for including additional model attributes";
#     option           type       default    group   requires    incompat   preamble-output                       help-output    
opt_Add("--group",     "string",  undef,        $g,  undef,         undef,  "specify model group is <s>",         "specify model group is <s>",    \%opt_HH, \@opt_order_A);
opt_Add("--subgroup",  "string",  undef,        $g,  "--group",     undef,  "specify model subgroup is <s>",      "specify model subgroup is <s>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling CDS translation step";
#     option          type       default    group   requires    incompat   preamble-output                                             help-output    
opt_Add("--ttbl",     "integer", 1,            $g,  undef,         undef,  "use NCBI translation table <n> to translate CDS",          "use NCBI translation table <n> to translate CDS", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling cmbuild step";
#     option          type       default    group   requires    incompat              preamble-output                                             help-output    
opt_Add("--cmn",      "integer", undef,       $g,   undef, "--skipbuild,--cminfile",  "set number of seqs for glocal fwd HMM calibration to <n>", "set number of seqs for glocal fwd HMM calibration to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--cmp7ml",   "boolean", 0,           $g,   undef, "--skipbuild,--cminfile",  "set CM's filter p7 HMM as the ML p7 HMM",                  "set CM's filter p7 HMM as the ML p7 HMM",                  \%opt_HH, \@opt_order_A);
opt_Add("--cmere",    "real",    undef,       $g,   undef,  "--skipbuild,--cminfile", "set CM relative entropy target to <x>",                    "set CM relative entropy target to <x>",                    \%opt_HH, \@opt_order_A);
opt_Add("--cmeset",   "real",    undef,       $g,   undef,  "--skipbuild,--cminfile", "set CM eff seq # for CM to <x>",                           "set CM eff seq # for CM to <x>",                           \%opt_HH, \@opt_order_A);
opt_Add("--cmemaxseq","real",    undef,       $g,   undef,  "--skipbuild,--cminfile", "set CM maximum allowed eff seq # for CM to <x>",           "set CM maximum alowed eff seq # for CM to <x>",            \%opt_HH, \@opt_order_A);
opt_Add("--cminfile", "string",  undef,       $g,   undef,  "--skipbuild",            "read cmbuild options from file <s>",                       "read cmbuild options from file <s>",                       \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for skipping stages";
#       option             type       default     group requires   incompat  preamble-output                                    help-output    
opt_Add("--skipbuild",     "boolean", 0,         $g,    undef,     undef,    "skip the cmbuild step",                           "skip the cmbuild and blastn db creation steps", \%opt_HH, \@opt_order_A);
opt_Add("--onlyurl",       "boolean", 0,         $g,    undef,"--stk,--ingb,--inft",  "output genbank file url for accession and exit",  "output genbank file url for accession and exit", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "optional output files";
#       option       type       default     group  requires     incompat  preamble-output                          help-output    
opt_Add("--ftrinfo",    "boolean", 0,         $g,    undef,     undef,    "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);
opt_Add("--sgminfo",    "boolean", 0,         $g,    undef,     undef,    "output internal segment information",   "create file with internal segment information", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "other expert options";
#       option       type          default     group  requires incompat  preamble-output                                 help-output    
opt_Add("--execname",   "string",  undef,         $g,    undef, undef,   "define executable name of this script as <s>", "define executable name of this script as <s>", \%opt_HH, \@opt_order_A);        

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'stk=s'        => \$GetOptions_H{"--stk"},
                'infa=s'       => \$GetOptions_H{"--infa"},
                'inft=s'       => \$GetOptions_H{"--inft"},
                'ftfetch1'     => \$GetOptions_H{"--ftfetch1"},
                'ftfetch2'     => \$GetOptions_H{"--ftfetch2"},
                'gb'           => \$GetOptions_H{"--gb"},
                'ingb=s'       => \$GetOptions_H{"--ingb"},
                'addminfo=s'   => \$GetOptions_H{"--addminfo"},
                'forcelong'    => \$GetOptions_H{"--forcelong"},
                'keep'         => \$GetOptions_H{"--keep"},
# options for controlling what feature types are stored in model info file
                'fall'         => \$GetOptions_H{"--fall"},
                'fadd=s'       => \$GetOptions_H{"--fadd"},
                'fskip=s'      => \$GetOptions_H{"--fskip"},
# options for controlling what qualifiers are stored in model info file
                'qall'         => \$GetOptions_H{"--qall"},
                'qadd=s'       => \$GetOptions_H{"--qadd"},
                'qftradd=s'    => \$GetOptions_H{"--qftradd"},
                'qskip=s'      => \$GetOptions_H{"--qskip"},
                'noaddgene'    => \$GetOptions_H{"--noaddgene"},
# options for including additional model attributes
                'group=s'      => \$GetOptions_H{"--group"},
                'subgroup=s'   => \$GetOptions_H{"--subgroup"},
# options for controlling CDS translation step
                'ttbl=s'       => \$GetOptions_H{"--ttbl"},
# options for controlling cmbuild step
                'cmn=s'        => \$GetOptions_H{"--cmn"},
                'cmp7ml'       => \$GetOptions_H{"--cmp7ml"},
                'cmere=s'      => \$GetOptions_H{"--cmere"},
                'cmeset=s'     => \$GetOptions_H{"--cmeset"},
                'cmemaxseq=s'  => \$GetOptions_H{"--cmemaxseq"},
                'cminfile=s'   => \$GetOptions_H{"--cminfile"},
# options for skipping stages
                'skipbuild'    => \$GetOptions_H{"--skipbuild"},
                'onlyurl'      => \$GetOptions_H{"--onlyurl"},
# optional output files
                'sgminfo'      => \$GetOptions_H{"--sgminfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"},
# other expert options
                'execname=s'   => \$GetOptions_H{"--execname"});

my $total_seconds = -1 * ofile_SecondsSinceEpoch(); # by multiplying by -1, we can just add another ofile_SecondsSinceEpoch call at end to get total time
my $execname_opt  = $GetOptions_H{"--execname"};
my $executable    = (defined $execname_opt) ? $execname_opt : "v-build.pl";
my $usage         = "Usage: $executable [-options] <accession> <path to output directory to create>\n";
my $synopsis      = "$executable :: build homology model of a single sequence for feature annotation";
my $date          = scalar localtime();
my $version       = "1.1dev3";
my $releasedate   = "June 2020";
my $pkgname       = "VADR";

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
  print "\nTo see more help on available options, do $executable -h\n\n";
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
  if(opt_Get("--gb", \%opt_HH)) { 
    print vdr_EutilsFetchUrl($mdl_name, "nuccore", "gb") . "\n";
  }
  else { 
    print vdr_EutilsFetchUrl($mdl_name, "nuccore", "ft") . "\n";
  }
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
my $out_root = $dir . "/" . $dir_tail . ".vadr";

#######################
# output program banner
#######################
# output preamble
my @arg_desc_A = ("accession/model name", "output directory");
my @arg_A      = ($mdl_name, $dir);
my %extra_H    = ();
$extra_H{"\$VADRSCRIPTSDIR"}  = $env_vadr_scripts_dir;
$extra_H{"\$VADRINFERNALDIR"} = $env_vadr_infernal_dir;
$extra_H{"\$VADREASELDIR"}    = $env_vadr_easel_dir;
$extra_H{"\$VADRBLASTDIR"}    = $env_vadr_blast_dir;
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
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "log",  $out_root . ".log",      1, 1, "Output printed to screen");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "cmd",  $out_root . ".cmd",      1, 1, "List of executed commands");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".filelist", 1, 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
my $FH_HR  = $ofile_info_HH{"FH"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# open optional output files
if(opt_Get("--ftrinfo", \%opt_HH)) { 
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "ftrinfo", $out_root . ".ftrinfo", 1, 1, "Feature information (created due to --ftrinfo)");
}
if(opt_Get("--sgminfo", \%opt_HH)) { 
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "sgminfo", $out_root . ".sgminfo", 1, 1, "Segment information (created due to --sgminfo)");
}

# now we have the log file open, output the banner there too
ofile_OutputBanner($log_FH, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

#######################################################
# Parse the input minfo file if --addminfo file is used
#######################################################
my $progress_w = 60; # the width of the left hand column in our progress output, hard-coded
my $start_secs;
my $addminfo_file = opt_Get("--addminfo", \%opt_HH);
my @add_mdl_info_AH  = ();
my %add_ftr_info_HAH = ();

if(defined $addminfo_file) { 
  $start_secs = ofile_OutputProgressPrior("Processing --addminfo option", $progress_w, $log_FH, *STDOUT);

  my @reqd_mdl_keys_A = ("name");
  my @reqd_ftr_keys_A = ("type");
  utl_FileValidateExistsAndNonEmpty($addminfo_file, "--addminfo file", undef, 1, $FH_HR);
  vdr_ModelInfoFileParse($addminfo_file, \@reqd_mdl_keys_A, \@reqd_ftr_keys_A, \@add_mdl_info_AH, \%add_ftr_info_HAH, $FH_HR);

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

###################################################
# Fetch the fasta file (if necessary) and parse it
###################################################
my $fa_file = $out_root . ".fa";
my %seq_H = ();
if(opt_IsUsed("--infa", \%opt_HH)) { 
  utl_RunCommand("cp " . opt_Get("--infa", \%opt_HH) . " $fa_file", opt_Get("-v", \%opt_HH), 0, $FH_HR);
}
else { 
  $start_secs = ofile_OutputProgressPrior("Fetching FASTA file", $progress_w, $log_FH, *STDOUT);
  vdr_EutilsFetchToFile($fa_file, $mdl_name, "nuccore", "fasta", 5, $ofile_info_HH{"FH"});  # number of attempts to fetch to make before dying
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "fasta", $fa_file, 1, 1, "fasta file for $mdl_name");
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}
$start_secs = ofile_OutputProgressPrior("Parsing FASTA file", $progress_w, $log_FH, *STDOUT);
vdr_ParseSeqFileToSeqHash($fa_file, \%seq_H, $FH_HR);
my @fetched_seq_A = (sort keys %seq_H);
if(scalar(@fetched_seq_A) != 1) { 
  ofile_FAIL("ERROR did not fetch exactly 1 sequence from fasta file $fa_file\n", 1, $FH_HR);
}
my $mdl_name_ver = $fetched_seq_A[0];
# make sure it's the right sequence
if($mdl_name_ver =~ /(\S+)\.\d+/) { 
  if($1 ne $mdl_name) { 
    ofile_FAIL("ERROR did not fetch correct sequence from fasta file $fa_file (expected accession.version starting with $mdl_name, got $mdl_name_ver)\n", 1, $FH_HR);
  }
}
else {
  ofile_FAIL("ERROR did not fetch correct sequence from fasta file $fa_file (expected accession.version starting with $mdl_name, got $mdl_name_ver)\n", 1, $FH_HR);
}
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#######################################################################
# Verify our sequence is not longer than our hard-coded maximum of 25Kb 
#######################################################################
# allow any length if --forcelong
my $maxlen = 25000;
my $mdllen = length($seq_H{$mdl_name_ver});
if(! opt_Get("--forcelong", \%opt_HH)) { 
  if($mdllen > 25000) { 
    ofile_FAIL("ERROR, model length ($mdllen) exceeds maximum allowed length of $maxlen.\nYou can use --forcelong to bypass this at your own risk.\nUse of VADR on models > 25Kb is not recommended.\nModel building will be very slow and\ndownstream v-annotate.pl annotation will have large memory requirements.", 1, $FH_HR);
  }
}

###################################################################
# Fetch the feature table (ft) or GenBank (gb) file (if necessary)
# and parse it.
###################################################################
my $ft_file = undef;
my $gb_file = undef;
my %ftr_info_HAH = (); # the feature info 
if(! opt_IsUsed("--gb", \%opt_HH)) { 
  if(opt_IsUsed("--inft", \%opt_HH)) { 
    $ft_file = opt_Get("--inft", \%opt_HH);
  }
  else { 
    # --inft not used, create ft file by fetching using eutils
    $start_secs = ofile_OutputProgressPrior("Fetching feature table file", $progress_w, $log_FH, *STDOUT);
    $ft_file = $out_root . ".tbl";
    if(opt_Get("--ftfetch1", \%opt_HH)) { 
      utl_RunCommand("efetch -db nuccore -id $mdl_name -format ft > $ft_file", opt_Get("-v", \%opt_HH), 0, $FH_HR);
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "ft", $ft_file, 1, 1, "feature table format file for $mdl_name (--ftfetch1)");
    }
    elsif(opt_Get("--ftfetch2", \%opt_HH)) { 
      utl_RunCommand("efetch -db nuccore -id $mdl_name -format gbc | xml2tbl > $ft_file", opt_Get("-v", \%opt_HH), 0, $FH_HR);
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "ft", $ft_file, 1, 1, "feature table format file for $mdl_name (--ftfetch2)");
    }
    else { # default way of fetching a feature table
      vdr_EutilsFetchToFile($ft_file, $mdl_name, "nuccore", "ft", 5, $ofile_info_HH{"FH"});  # number of attempts to fetch to make before dying
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "ft", $ft_file, 1, 1, "feature table format file for $mdl_name");
    }
    ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
  # parse the feature table file
  $start_secs = ofile_OutputProgressPrior("Parsing feature table file", $progress_w, $log_FH, *STDOUT);
  sqf_FeatureTableParse($ft_file, \%ftr_info_HAH, $FH_HR);
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  # if we have any CDS features with protein_id qualifiers, fetch and parse those
  $start_secs = ofile_OutputProgressPrior("Fetching and parsing protein feature table file(s)", $progress_w, $log_FH, *STDOUT);
  fetch_and_parse_cds_protein_feature_tables(\@{$ftr_info_HAH{$mdl_name}}, $out_root, $FH_HR);
  if(! exists $ftr_info_HAH{$mdl_name}) { 
    ofile_FAIL("ERROR parsing GenBank file $gb_file, did not read info for reference accession $mdl_name\n", 1, $FH_HR);
  }
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
} # end of 'if(! opt_IsUsed("--gb", \%opt_HH))' { 
else { 
# If --gb (and not --ingb) used: fetch the genbank file 
  if(opt_IsUsed("--ingb", \%opt_HH)) { 
    $gb_file = opt_Get("--ingb", \%opt_HH);
  }
  else { 
    # --ingb not used, create gb file by fetching using eutils
    $start_secs = ofile_OutputProgressPrior("Fetching GenBank file", $progress_w, $log_FH, *STDOUT);
    
    $gb_file = $out_root . ".gb";
    vdr_EutilsFetchToFile($gb_file, $mdl_name, "nuccore", "gb", 5, $ofile_info_HH{"FH"});  # number of attempts to fetch to make before dying
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "gb", $gb_file, 1, 1, "GenBank format file for $mdl_name");
    
    ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
  # parse the genbank file
  $start_secs = ofile_OutputProgressPrior("Parsing GenBank file", $progress_w, $log_FH, *STDOUT);
  sqf_GenbankParse($gb_file, undef, \%ftr_info_HAH, $FH_HR);
  if(! exists $ftr_info_HAH{$mdl_name}) { 
    ofile_FAIL("ERROR parsing GenBank file $gb_file, did not read info for reference accession $mdl_name\n", 1, $FH_HR);
  }
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}
if(exists $ofile_info_HH{"FH"}{"ftrinfo"}) { 
  utl_AHDump("Feature information", \@{$ftr_info_HAH{$mdl_name}}, $ofile_info_HH{"FH"}{"ftrinfo"});
}

#######################################################
# Prune data read from %ftr_info_HAH, only keeping what
# we want to output to the eventual model info file
#######################################################
$start_secs = ofile_OutputProgressPrior("Pruning data read from GenBank", $progress_w, $log_FH, *STDOUT);

# determine what types of features we will store based on cmdline options
# --fall is incompatible with --fadd
my %fdf_H   = (); # default feature types to keep
my %fadd_H  = (); # feature types to add
my %fskip_H = (); # feature types to skip
process_add_and_skip_options("CDS,gene,mat_peptide", "--fadd", "--fskip", undef, \%fdf_H, \%fadd_H, \%fskip_H, undef, \%opt_HH, $FH_HR);

# determine what qualifiers we will store based on cmdline options
# --qall is incompatible with --qadd
my %qdf_H      = (); # default qualifiers to keep
my %qadd_H     = (); # qualifiers to add
my %qskip_H    = (); # qualifiers to skip
my %qftr_add_H = (); # if --qftradd, subset of features to add qualifiers in --qadd option for
process_add_and_skip_options("type,coords,location,product,gene,exception,parent_idx_str,5p_trunc,3p_trunc", "--qadd", "--qskip", "--qftradd", \%qdf_H, \%qadd_H, \%qskip_H, \%qftr_add_H, \%opt_HH, $FH_HR); 
# we only need ribosomal_slippage above so we can get the exception:ribosomal slippage 
# qualifier, if we switch to parsing feature tables instead of GenBank files, then
# "ribosomal_slippage" should be removed from the list.

# remove all features types we don't want
my $ftr_idx;
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

# deal with special case: remove any CDS features that have "trunc5"
# or "trunc3" keys set as 1 we can't deal with these because we
# don't know how to translate them in v-build.pl and (even if we did
# handle that based on codon_start) v-annotate.pl can't deal with
# these because a start/stop codon is not expected and all complete
# CDS are validated by looking for a start/stop
for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$mdl_name}}); $ftr_idx++) { 
  my $ftype = $ftr_info_HAH{$mdl_name}[$ftr_idx]{"type"};
  if(($ftype eq "CDS") && 
     (((defined $ftr_info_HAH{$mdl_name}[$ftr_idx]{"trunc5"}) && 
       ($ftr_info_HAH{$mdl_name}[$ftr_idx]{"trunc5"} == 1)) || 
      ((defined $ftr_info_HAH{$mdl_name}[$ftr_idx]{"trunc3"}) && 
       ($ftr_info_HAH{$mdl_name}[$ftr_idx]{"trunc3"} == 1)))) { 
    ofile_OutputString($log_FH, 1, "\n# WARNING: not modelling CDS feature with coords " . $ftr_info_HAH{$mdl_name}[$ftr_idx]{"coords"} . " because it is 5' and/or 3' truncated\n#          (e.g. incomplete, with a \"<\" or \">\" in its coordinates in the feature table.\n#\n# ");
    splice(@{$ftr_info_HAH{$mdl_name}}, $ftr_idx, 1);
    $ftr_idx--; # this is about to be incremented
  }
}

# remove any qualifier key/value pairs with keys not in %qual_H, unless --qall used
for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$mdl_name}}); $ftr_idx++) { 
  my $ftype = $ftr_info_HAH{$mdl_name}[$ftr_idx]{"type"};
  foreach my $qual (sort keys %{$ftr_info_HAH{$mdl_name}[$ftr_idx]}) { 
    # we skip this qualifier and remove it from ftr_info_HAH
    # if all three of A1, A2, A3 OR B is satisfied
    # (A1) it's not a default qualifier        AND
    # (A2) (it's not listed in --qadd OR 
    #       (--qftradd is used AND $ftype is not listed in --qftradd)) AND
    # (A3) --qall not used
    # OR 
    # (B) it is listed in --qskip string 
    if(((! defined $qdf_H{$qual})         && # (A1)
        ((! defined $qadd_H{$qual}) || 
         ((opt_IsUsed("--qftradd", \%opt_HH)) && 
          (! defined $qftr_add_H{$ftype})))   && # (A2)
        (! opt_Get("--qall", \%opt_HH)))     # (A3)
       || (defined $qskip_H{$qual})) {       # (B)
      delete $ftr_info_HAH{$mdl_name}[$ftr_idx]{$qual};
    }
  }
}

if(opt_Get("--gb", \%opt_HH)) { 
  # Deal with special case: we purposefully added 'ribosomal_slippage' qualifiers if they 
  # existed just so we could now add 'exception' qualifiers with 'ribosomal slippage' values
  # at this stage. This is ONLY to get around problem that GenBank format includes 'ribosomal_slippage'
  # qualifiers but not 'exception' qualifiers with 'ribosomal slippage' values, but 
  # only 'exception:ribosomal slippage' qualifier/values are desired in the 
  # output feature table. If we switch to parsing Entrez feature tables as input then
  # the need for this should go away because 'exception:ribosomal slippage' is in that
  # feature table file (along with the unwanted 'ribosomal_slippage' qualifier which
  # we can just ignore). 
  # 
  # If ribosomal_slippage qualifier exists: create a new "exception" 
  # qualifier with value of "ribosomal slippage"
  # 
  for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$mdl_name}}); $ftr_idx++) { 
    if((defined $ftr_info_HAH{$mdl_name}[$ftr_idx]) && 
       (defined $ftr_info_HAH{$mdl_name}[$ftr_idx]{"ribosomal_slippage"})) {
      if(defined $ftr_info_HAH{$mdl_name}[$ftr_idx]{"exception"}) { 
        $ftr_info_HAH{$mdl_name}[$ftr_idx]{"exception"} .= ":GBSEP:" . "ribosomal slippage";
      }
      else { 
        $ftr_info_HAH{$mdl_name}[$ftr_idx]{"exception"} = "ribosomal slippage";
      }
      # now remove the "ribosomal_slippage" qualifier UNLESS --qall used or $qadd_H{"ribosomal_slippage"} exists
      if((! opt_Get("--qall", \%opt_HH)) &&
         (! defined $qadd_H{"ribosomal_slippage"})) { 
        delete $ftr_info_HAH{$mdl_name}[$ftr_idx]{"ribosomal_slippage"};
      }
    }
  }
}
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###############################################
# Add in features read from --addminfo if used
###############################################
# if --addminfo was used, add the feature info read from that file
if(defined $addminfo_file) { 
  $start_secs = ofile_OutputProgressPrior("Adding feature info read from --addminfo file", $progress_w, $log_FH, *STDOUT);
  if(! defined $add_ftr_info_HAH{$mdl_name}) { 
    ofile_FAIL("ERROR with --addminfo <s>, <s> must include model $mdl_name, but it does not", 1, $FH_HR);
  }
  vdr_FeatureInfoMerge(\@{$add_ftr_info_HAH{$mdl_name}}, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

#####################################################################
# Parse the input stockholm file (if --stk) or create it (if ! --stk)
#####################################################################
my $stk_file = $out_root . ".stk";
my $stk_has_ss = undef;
my $in_stk_file = opt_Get("--stk", \%opt_HH);
if(defined $in_stk_file) { 
  $start_secs = ofile_OutputProgressPrior("Validating input Stockholm file", $progress_w, $log_FH, *STDOUT);

  $stk_has_ss = stockholm_validate_single_sequence_input($in_stk_file, $seq_H{$mdl_name_ver}, \%opt_HH, $FH_HR);

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  utl_RunCommand("cp $in_stk_file $stk_file", opt_Get("-v", \%opt_HH), 0, $FH_HR);
}
else { 
  # --stk not used, we create it from the fasta file we read
  $start_secs = ofile_OutputProgressPrior("Reformatting FASTA file to Stockholm file", $progress_w, $log_FH, *STDOUT);

  sqf_EslReformatRun($execs_H{"esl-reformat"}, undef, $fa_file, $stk_file, "afa", "stockholm", \%opt_HH, $FH_HR);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "stk", $stk_file, 1, 1, "Stockholm alignment file for $mdl_name");

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

######################################################################
# Finish populating @{$ftr_info_HAH{$mdl_name} and create @sgm_info_AH
######################################################################
$start_secs = ofile_OutputProgressPrior("Finalizing feature information", $progress_w, $log_FH, *STDOUT);

if(opt_Get("--gb", \%opt_HH)) { # we only need to derive 'coords' if we parsed the GenBank file
  vdr_FeatureInfoImputeCoords(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
}
vdr_FeatureInfoImputeLength(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
vdr_FeatureInfoInitializeParentIndexStrings(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);

# A special step only needed in v-build.pl (not needed in v-annotate.pl): 
# Convert parent_index_str values from the strings they were set to in
# fetch_and_parse_cds_protein_feature_tables to integers, now that all
# feature pruning is complete
integerize_parent_index_strings(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);

vdr_FeatureInfoImputeOutname(\@{$ftr_info_HAH{$mdl_name}});
# add 'gene' qualifiers to 'CDS' features
if((! opt_Get("--noaddgene", \%opt_HH)) && (! defined $qskip_H{"gene"})) { 
  vdr_FeatureInfoImputeByOverlap(\@{$ftr_info_HAH{$mdl_name}}, "gene", "gene", "CDS",        "gene", $FH_HR);
  vdr_FeatureInfoImputeByOverlap(\@{$ftr_info_HAH{$mdl_name}}, "gene", "gene", "mRNA",       "gene", $FH_HR);
  vdr_FeatureInfoImputeByOverlap(\@{$ftr_info_HAH{$mdl_name}}, "gene", "gene", "regulatory", "gene", $FH_HR);
}

my @sgm_info_AH = (); # segment info, inferred from feature info
vdr_SegmentInfoPopulate(\@sgm_info_AH, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);

ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###################################
# If any CDS: 
# - translate the CDS
# - build BLAST DB
# - build HMMER DB
###################################
my $ncds = vdr_FeatureInfoCountType(\@{$ftr_info_HAH{$mdl_name}}, "CDS");
my $cds_fa_file = undef;
my $protein_fa_file = undef;
my $hmm_file = undef;
my $hmmbuild_file = undef;
if($ncds > 0) { 
  # translate CDS
  $start_secs = ofile_OutputProgressPrior("Translating CDS ", $progress_w, $log_FH, *STDOUT);

  $cds_fa_file  = $out_root . ".cds.fa";
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "cdsfasta", $cds_fa_file, 1, 1, "fasta sequence file for CDS from $mdl_name");
  vdr_CdsFetchStockholmToFasta($ofile_info_HH{"FH"}{"cdsfasta"}, $stk_file, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  close $ofile_info_HH{"FH"}{"cdsfasta"};
  
  $protein_fa_file = $out_root . ".protein.fa";
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "proteinfasta", $protein_fa_file, 1, 1, "fasta sequence file for translated CDS from $mdl_name");
  sqf_EslTranslateCdsToFastaFile($ofile_info_HH{"FH"}{"proteinfasta"}, $execs_H{"esl-translate"}, $cds_fa_file, 
                                 $out_root, \@{$ftr_info_HAH{$mdl_name}}, \%opt_HH, $FH_HR);
  close $ofile_info_HH{"FH"}{"proteinfasta"};
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  # build blast db
  $start_secs = ofile_OutputProgressPrior("Building BLAST protein database ", $progress_w, $log_FH, *STDOUT);

  sqf_BlastDbCreate($execs_H{"makeblastdb"}, "prot", $protein_fa_file, \%opt_HH, $FH_HR);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-phr", $protein_fa_file . ".phr", 1, 1, "BLAST db .phr file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-pin", $protein_fa_file . ".pin", 1, 1, "BLAST db .pin file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-psq", $protein_fa_file . ".psq", 1, 1, "BLAST db .psq file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-pdb", $protein_fa_file . ".pdb", 1, 1, "BLAST db .pdb file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-pot", $protein_fa_file . ".pot", 1, 1, "BLAST db .pot file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-ptf", $protein_fa_file . ".ptf", 1, 1, "BLAST db .ptf file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-pto", $protein_fa_file . ".pto", 1, 1, "BLAST db .pto file for $mdl_name");

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  # build hmmer db, need to build one HMM per CDS and concatenate them
  $start_secs = ofile_OutputProgressPrior("Building HMMER protein database ", $progress_w, $log_FH, *STDOUT);

  # run esl-seqstat and parse it
  my $sfetch_cmd = $execs_H{"esl-sfetch"} . " --index $protein_fa_file";
  my $protein_sqfile = Bio::Easel::SqFile->new({ fileLocation => $protein_fa_file });
  my $nhmm = 0;
  my @hmm_file_A      = (); # array of names of individual HMM files to be concatenated to make the library
  my @hmmbuild_file_A = (); # array of names of individual hmmbuild output files to be concatenated to make the library
  for(my $i = 0; $i < $ncds; $i++) { 
    my $seq_name = $protein_sqfile->fetch_seq_name_given_ssi_number($i);
    my $hmm_name = undef;
    # remove version from $hmm_name
    if($seq_name =~ /^(.+)\.\d+(\/[^\/]+)$/) { 
      $hmm_name = $1 . $2;
    }
    else { 
      ofile_FAIL("ERROR, unable to parse protein sequence name $seq_name to make HMM model name", 1, $FH_HR);
    }
    my $tmp_hmm_file      = $out_root . "." . ($nhmm+1) . ".hmm";
    my $tmp_hmmbuild_file = $out_root . "." . ($nhmm+1) . ".hmmbuild";
    my $sfetch_to_hmmbuild_cmd = $execs_H{"esl-sfetch"} . " $protein_fa_file $seq_name | ";
    $sfetch_to_hmmbuild_cmd   .= $execs_H{"hmmbuild"} . " -n $hmm_name --informat afa $tmp_hmm_file - > $tmp_hmmbuild_file";
    utl_RunCommand($sfetch_to_hmmbuild_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
    push(@hmm_file_A,      $tmp_hmm_file);
    push(@hmmbuild_file_A, $tmp_hmmbuild_file);
    $nhmm++;
  }
  $hmm_file      = $out_root . ".protein.hmm";
  $hmmbuild_file = $out_root . ".protein.hmmbuild";
  utl_ConcatenateListOfFiles(\@hmm_file_A,      $hmm_file,      "v-build.pl main()", \%opt_HH, $FH_HR);
  utl_ConcatenateListOfFiles(\@hmmbuild_file_A, $hmmbuild_file, "v-build.pl main()", \%opt_HH, $FH_HR);

  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "proteinindex", $protein_fa_file . ".ssi", 0, 1, "esl-sfetch index file for $protein_sqfile");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "hmmdb",        $hmm_file,                 1, 1, "HMMER model db file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "hmmbuild",     $hmmbuild_file,            1, 1, "hmmbuild build output (concatenated)");

  # run hmmpress on it
  my $hmmpress_file = $out_root . ".hmmpress";
  my $hmmpress_cmd  = $execs_H{"hmmpress"} . " $hmm_file > $hmmpress_file";
  utl_RunCommand($hmmpress_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);

  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "h3m",      $hmm_file . ".h3m", 1, 1, "binary HMM and p7 HMM filter file");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "h3i",      $hmm_file . ".h3i", 1, 1, "SSI index for binary HMM file");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "h3f",      $hmm_file . ".h3f", 1, 1, "optimized p7 HMM filters (MSV part)");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "h3p",      $hmm_file . ".h3p", 1, 1, "optimized p7 HMM filters (remainder)");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "hmmpress", $hmmpress_file,     1, 1, "hmmpress output file");

  ofile_OutputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);
}

##############
# Build the CM
##############
my $cm_file = undef;
if(! opt_Get("--skipbuild", \%opt_HH)) { 
  my $cmbuild_str = undef;
  my $clen_times_cmn = length($seq_H{$mdl_name_ver}) * 200;
  if(opt_IsUsed("--cmn", \%opt_HH)) { 
    $clen_times_cmn *= (opt_Get("--cmn", \%opt_HH) / 200);
  }
  if   ($clen_times_cmn > 4000000) { $cmbuild_str = "(may take more than an hour)"; }
  elsif($clen_times_cmn > 3000000) { $cmbuild_str = "(should take roughly an hour)"; }
  elsif($clen_times_cmn > 2000000) { $cmbuild_str = "(should take roughly 20-40 minutes)"; }
  elsif($clen_times_cmn > 1000000) { $cmbuild_str = "(should take roughly 10-30 minutes)"; }
  elsif($clen_times_cmn >  500000) { $cmbuild_str = "(should take roughly 5-10 minutes)"; }
  else                             { $cmbuild_str = "(shouldn't take more than a few minutes)"; }

  $start_secs = ofile_OutputProgressPrior("Building CM $cmbuild_str", $progress_w, $log_FH, *STDOUT);

  my $cmbuild_occ_file = $out_root . ".cmbuild.occ";
  my $cmbuild_cp9occ_file = $out_root . ".cmbuild.cp9occ";
  my $cmbuild_fp7occ_file = $out_root . ".cmbuild.fp7occ";

  my $cmbuild_opts = "-n $mdl_name --verbose --occfile $cmbuild_occ_file --cp9occfile $cmbuild_cp9occ_file --fp7occfile $cmbuild_fp7occ_file ";
  if((! defined $stk_has_ss) || (! $stk_has_ss)) { $cmbuild_opts .= " --noss"; }
  if(opt_IsUsed("--cmn",       \%opt_HH)) { $cmbuild_opts .= " --EgfN "    . opt_Get("--cmn", \%opt_HH); }
  if(opt_IsUsed("--cmp7ml",    \%opt_HH)) { $cmbuild_opts .= " --p7ml"; }
  if(opt_IsUsed("--cmere",     \%opt_HH)) { $cmbuild_opts .= " --ere "     . opt_Get("--cmere", \%opt_HH); }
  if(opt_IsUsed("--cmeset",    \%opt_HH)) { $cmbuild_opts .= " --eset "    . opt_Get("--cmeset", \%opt_HH); }
  if(opt_IsUsed("--cmemaxseq", \%opt_HH)) { $cmbuild_opts .= " --emaxseq " . opt_Get("--cmemaxseq", \%opt_HH); }
  if(opt_IsUsed("--cminfile",  \%opt_HH)) { 
    my @cminfile_A = ();
    utl_FileLinesToArray(opt_Get("--cminfile", \%opt_HH), 1, \@cminfile_A, $FH_HR);
    foreach my $optline (@cminfile_A) { 
      chomp $optline;
      $cmbuild_opts .= " " . $optline . " ";
    }
  }
  # if model is big > 0.5 * 25Kb (maxlen), then use the --Egcmult option
  # this avoids problems and slowness with very large sequence lengths for glocal HMM calibration
  if($mdllen > (0.5 * $maxlen)) { 
    $cmbuild_opts .= " --Egcmult " . sprintf("%.5f", ($maxlen / $mdllen));
  }

  my $cmbuild_file = $out_root . ".cmbuild";
  $cm_file         = $out_root . ".cm";
  my $cmbuild_cmd  = $execs_H{"cmbuild"} . " " . $cmbuild_opts . " $cm_file $stk_file > $cmbuild_file";
  utl_RunCommand($cmbuild_cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});
  ofile_OutputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);

  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "cm",      $cm_file, 1, 1, "CM file");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "cmbuild", $cmbuild_file, 1, 1, "cmbuild output file");

  # press the file we just created 
  $start_secs = ofile_OutputProgressPrior("Pressing CM file", $progress_w, $log_FH, *STDOUT);
  my $cmpress_file = $out_root . ".cmpress";
  my $cmpress_cmd  = $execs_H{"cmpress"} . " $cm_file > $cmpress_file";
  utl_RunCommand($cmpress_cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});
  ofile_OutputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);

  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "i1m",     $cm_file . ".i1m", 1, 1, "binary CM and p7 HMM filter file");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "i1i",     $cm_file . ".i1i", 1, 1, "SSI index for binary CM file");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "i1f",     $cm_file . ".i1f", 1, 1, "optimized p7 HMM filters (MSV part)");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "i1p",     $cm_file . ".i1p", 1, 1, "optimized p7 HMM filters (remainder)");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "cmpress", $cmpress_file,     1, 1, "cmpress output file");
}

#####################################################################
# Build the blastn database, using the consensus sequence from the CM
#####################################################################
if(! opt_Get("--skipbuild", \%opt_HH)) { # we can only do this step if we built the CM
  
  $start_secs = ofile_OutputProgressPrior("Building BLAST nucleotide database of CM consensus ", $progress_w, $log_FH, *STDOUT);

  # emit the consensus sequence to a file
  my $cmemit_fa_file = $out_root . ".nt-cseq.fa";
  my $cmemit_cseq = vdr_CmemitConsensus(\%execs_H, $cm_file, $mdl_name, $cmemit_fa_file, \%opt_HH, \%ofile_info_HH);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "cmemit-fasta", $cmemit_fa_file, 0, opt_Get("--keep", \%opt_HH), "cmemit -c output");

  # make a new fasta file with same sequence but new name ($mdl_name)
  my $blastn_fa_file = $out_root . ".nt.fa";
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "blastn-fa",  $blastn_fa_file, 1, 1, "fasta sequence file with cmemit consensus sequence for $mdl_name");
  printf { $ofile_info_HH{"FH"}{"blastn-fa"} } ">" . $mdl_name . "\n" . seq_SqstringAddNewlines($cmemit_cseq, 60);
  close $ofile_info_HH{"FH"}{"blastn-fa"};

  sqf_BlastDbCreate($execs_H{"makeblastdb"}, "nucl", $blastn_fa_file, \%opt_HH, $FH_HR);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nhr", $blastn_fa_file . ".nhr", 1, 1, "BLAST db .nhr file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nin", $blastn_fa_file . ".nin", 1, 1, "BLAST db .nin file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nsq", $blastn_fa_file . ".nsq", 1, 1, "BLAST db .nsq file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-ndb", $blastn_fa_file . ".ndb", 1, 1, "BLAST db .ndb file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-not", $blastn_fa_file . ".not", 1, 1, "BLAST db .not file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-ntf", $blastn_fa_file . ".ntf", 1, 1, "BLAST db .ntf file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nto", $blastn_fa_file . ".nto", 1, 1, "BLAST db .nto file for $mdl_name");


  if(! opt_Get("--keep", \%opt_HH)) { 
    utl_FileRemoveUsingSystemRm($cmemit_fa_file, "v-build.pl main", \%opt_HH, $FH_HR);
  }

  ofile_OutputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);
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
$mdl_info_AH[0]{"length"} = length($seq_H{$mdl_name_ver});
if(defined $cm_file) { 
  $mdl_info_AH[0]{"cmfile"} = utl_RemoveDirPath($cm_file);
}
if($ncds > 0) { 
  $mdl_info_AH[0]{"blastdb"} = utl_RemoveDirPath($protein_fa_file);
  if((opt_IsUsed("--ttbl", \%opt_HH)) && (opt_Get("--ttbl", \%opt_HH) != 1))  { 
    $mdl_info_AH[0]{"transl_table"} = opt_Get("--ttbl", \%opt_HH);
  }
}
if(opt_IsUsed("--group", \%opt_HH)) { 
  $mdl_info_AH[0]{"group"} = opt_Get("--group", \%opt_HH); 
  if(opt_IsUsed("--subgroup", \%opt_HH)) { 
    $mdl_info_AH[0]{"subgroup"} = opt_Get("--subgroup", \%opt_HH); 
  }
}
my $modelinfo_file  = $out_root . ".minfo";
vdr_ModelInfoFileWrite($modelinfo_file, \@mdl_info_AH, \%ftr_info_HAH, $FH_HR);
ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "modelinfo", $modelinfo_file, 1, 1, "VADR 'model info' format file for $mdl_name");

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
ofile_OutputConclusionAndCloseFilesOk($total_seconds, $dir, \%ofile_info_HH);
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
#  $opt_HHR:      REF to 2D hash of option values, see top of sqp_opts.pm for description, PRE-FILLED
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

  if(! -e $in_stk_file) { ofile_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file does not exist", 1, $FH_HR); }
  if(! -s $in_stk_file) { ofile_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file exists but is empty", 1, $FH_HR); }
  if(  -d $in_stk_file) { ofile_FAIL("ERROR, --stk enabled, stockholm file $in_stk_file is actually a directory", 1, $FH_HR); }
  my $msa = Bio::Easel::MSA->new({ fileLocation => $in_stk_file, isDna => 1});
  my $nseq = $msa->nseq;
  if($nseq == 1) { 
    # single sequence, make sure there are no gaps
    if($msa->any_allgap_columns) { 
      ofile_FAIL("ERROR, read 1 sequence in --stk file $in_stk_file, but it has gaps, this is not allowed for single sequence 'alignments' (remove gaps with 'esl-reformat --mingap')", 1, $FH_HR);
    }
    # validate it matches $exp_sqstring
    my $fetched_sqstring = $msa->get_sqstring_unaligned(0);
    seq_SqstringCapitalize(\$fetched_sqstring);
    seq_SqstringCapitalize(\$exp_sqstring);
    seq_SqstringDnaize(\$fetched_sqstring);
    seq_SqstringDnaize(\$exp_sqstring);
    if($fetched_sqstring ne $exp_sqstring) { 
      my $summary_sqstring_diff_str = seq_SqstringDiffSummary($fetched_sqstring, $exp_sqstring);
      ofile_FAIL("ERROR, read 1 sequence in --stk file $in_stk_file, but it does not match sequence read from GenBank file $gb_file:\n$summary_sqstring_diff_str", 1, $FH_HR); 
    }
  }
  else { # nseq != 1
    ofile_FAIL("ERROR, did not read exactly 1 sequence in --stk file $in_stk_file.\nTo use VADR with models built from alignments of multiple sequences,\nyou will have to build the CM with cmbuild and create the model info file manually.\n", 1, $FH_HR);
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
#  $sub_opt:    name of option with subset of features add option applies to, or undef
#  $df_HR:      ref to hash of default keys, filled here, values will all be '1'
#  $add_HR:     ref to hash of keys to add, filled here, values will all be '1'
#  $skip_HR:    ref to hash of keys to skip, filled here, values will all be '1'
#  $add_sub_HR: ref to hash of feature keys $add_HR applies to, or undef if $sub_opt is undef
#  $opt_HHR:    ref to hash of cmdline options
#  $FH_HR:      ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if a key is listed in both the $add_opt and $skip_opt option strings.
#################################################################
sub process_add_and_skip_options { 
  my $sub_name = "process_add_and_skip_options";
  my $nargs_expected = 10;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($df_string, $add_opt, $skip_opt, $sub_opt, $df_HR, $add_HR, $skip_HR, $add_sub_HR, $opt_HHR, $FH_HR) = @_;

  utl_ExistsHFromCommaSepString($df_HR, $df_string);
  if(opt_IsUsed($add_opt,  $opt_HHR)) { utl_ExistsHFromCommaSepString($add_HR,  opt_Get($add_opt,  $opt_HHR)); }
  if(opt_IsUsed($skip_opt, $opt_HHR)) { utl_ExistsHFromCommaSepString($skip_HR, opt_Get($skip_opt, $opt_HHR)); }
  if((defined $sub_opt) && (opt_IsUsed($sub_opt, $opt_HHR))) { utl_ExistsHFromCommaSepString($add_sub_HR, opt_Get($sub_opt, $opt_HHR)); }
  # make sure $add_opt and $skip_opt have no values in common
  foreach my $key (sort keys (%{$add_HR})) { 
    if(defined $skip_HR->{$key}) { 
      ofile_FAIL("ERROR in $sub_name, processing $add_opt <s1> and $skip_opt <s2> options, $key exists in both <s1> and <s2>", 1, $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine: fetch_and_parse_cds_protein_feature_tables()
# Incept:     EPN, Tue May 21 20:49:40 2019
#
# Synopsis: Fetch and parse feature tables for proteins stored
#           as qualifier values for the "protein_id" qualifier
#           of "CDS" features in @{$ftr_info_AHR}. Features
#           and qualifiers read from these feature tables are
#           added to @{$ftr_info_AHR} after converting the
#           coordinates as necessary.
#
# Arguments:
#  $ftr_info_AHR:  ref to the feature info array of hashes
#  $out_root:      output root for the file names
#  $FH_HR:         ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's a problem fetching of parsing any
#             protein feature tables
#################################################################
sub fetch_and_parse_cds_protein_feature_tables { 
  my $sub_name = "fetch_and_parse_cds_protein_feature_tables()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($ftr_info_AHR, $out_root, $FH_HR) = @_;

  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my %prot_ftr_info_HAH = ();
    my $tmp_parent_idx_str = $ftr_info_AHR->[$ftr_idx]{"type"} . ":GBSEP:" . $ftr_info_AHR->[$ftr_idx]{"coords"}; # we will add this to any child's parent_idx_str value
    if(($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") && 
       (defined $ftr_info_AHR->[$ftr_idx]{"protein_id"})) { 
      my $protein_id = $ftr_info_AHR->[$ftr_idx]{"protein_id"};
      my $accver = undef;
      if($protein_id =~ /[^\|]*\|([^\|]+\.\d+)\|/) { 
        $accver = $1;
      }
      elsif($protein_id =~ /([^\|]+\.\d+)/) { 
        $accver = $1;
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, unable to parse protein_id $protein_id to get accession.version\n", 1, $FH_HR);
      }
      my $ft_file = $out_root . "." . $accver . ".tbl";
      vdr_EutilsFetchToFile($ft_file, $accver, "protein", "ft", 5, $ofile_info_HH{"FH"});  # number of attempts to fetch to make before dying
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "ft." . $accver, $ft_file, 1, 1, "feature table format file for $accver");

      # parse the file
      sqf_FeatureTableParse($ft_file, \%prot_ftr_info_HAH, $FH_HR);
    }
    # copy info from \%prot_ftr_info_HAH to %ftr_info_AHR,
    # but only if we don't already have that feature in %ftr_info_AHR
    foreach my $prot_accver (sort keys %prot_ftr_info_HAH) { 
      my $prot_nftr = utl_AHValidate(\@{$prot_ftr_info_HAH{$prot_accver}}, \@keys_A, "ERROR in $sub_name for accver $prot_accver", $FH_HR);
      for(my $prot_ftr_idx = 0; $prot_ftr_idx < $prot_nftr; $prot_ftr_idx++) { 
        # for non-CDS features, check to see if we already have a
        # feature with the same type and coords, if so, make sure all
        # data is consistent and skip it (if not all data is
        # consistent: die in error) if no other feature with same type
        # and coords exists, add it
        # (We skip all CDS because we should already have them from the nucleotide
        #  record, and because our check to see if an existing feature exists doesn't
        #  word because the coords will differ by 3 and the 3' end due to the stop
        #  codon coords being included in the nucleotide CDS record, but not the
        #  protein one.)
        # first, convert protein coords to nucleotide coords (before
        # checking if it already exists or not)
        if($prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"type"} ne "CDS") { 
          $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"coords"} = vdr_CoordsProteinRelativeToAbsolute($ftr_info_AHR->[$ftr_idx]{"coords"}, $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"coords"}, $FH_HR);
          my $found_ftr_idx = -1;
          for(my $chk_ftr_idx = 0; $chk_ftr_idx < scalar(@{$ftr_info_AHR}); $chk_ftr_idx++) { 
            if($found_ftr_idx == -1) { 
              if(($prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"type"}   eq $ftr_info_AHR->[$chk_ftr_idx]{"type"}) && 
                 ($prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"coords"} eq $ftr_info_AHR->[$chk_ftr_idx]{"coords"})) { 
                # add data from $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx] to ftr_info_AHR->[$chk_ftr_idx]
                foreach my $prot_key (sort keys (%{$prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]})) { 
                  if(! defined $ftr_info_AHR->[$chk_ftr_idx]{$prot_key}) { 
                    $ftr_info_AHR->[$chk_ftr_idx]{$prot_key} = $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{$prot_key};
                  }
                  else { 
                    if($ftr_info_AHR->[$chk_ftr_idx]{$prot_key} ne $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{$prot_key}) { 
                      # not equivalent, append 
                      $ftr_info_AHR->[$chk_ftr_idx]{$prot_key} .= ":GBSEP:" . $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{$prot_key};
                    }
                  }
                }
              }
            }
          }
          if($found_ftr_idx != -1) { # the feature already exists update its parent string
            if((! defined $ftr_info_AHR->[$found_ftr_idx]{"parent_idx_str"}) || 
               ($ftr_info_AHR->[$found_ftr_idx]{"parent_idx_str"} eq "GBNULL")) { 
              # set parent_idx_str to "parent's type" . ":GBSEP:" . "parent's coords", we need to do this because parent's ftr_idx may change when we prune unwanted features
              $ftr_info_AHR->[$found_ftr_idx]{"parent_idx_str"} = $tmp_parent_idx_str;
            }
            else { 
              # set parent_idx_str to "parent's type" . ":GBSEP:" . "parent's coords", we need to do this because parent's ftr_idx may change when we prune unwanted features
              $ftr_info_AHR->[$found_ftr_idx]{"parent_idx_str"} .= "!GBSEP!" . $tmp_parent_idx_str;
            }
          }
          else { # we didn't find this feature already in the feature info hash, add it
            # printf("adding feature " . $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"type"} . " with coords " . $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"coords"} . "\n");
            my $nxt_ftr_idx = scalar(@{$ftr_info_AHR});
            %{$ftr_info_AHR->[$nxt_ftr_idx]} = ();
            foreach my $prot_key (sort keys (%{$prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]})) { 
              $ftr_info_AHR->[$nxt_ftr_idx]{$prot_key} = $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{$prot_key};
            }
            # set parent_idx_str to "parent's type" . ":GBSEP:" . "parent's coords", we need to do this because parent's ftr_idx may change when we prune unwanted features
            $ftr_info_AHR->[$nxt_ftr_idx]{"parent_idx_str"} = $tmp_parent_idx_str;
          }
        } # end of 'if($prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"type"} ne "CDS") {'
      }
    }
  }
  return;
}

#################################################################
# Subroutine: integerize_parent_index_strings
# Incept:     EPN, Wed Aug 14 06:38:00 2019
# 
# Purpose:    Update "parent_idx_str" values that are set as 
#             N >= 1 "!GBSEP!" separated tokens of: 
#             <parent's type> . ":GBSEP:" . <parent's coords> to
#             a string of N integers separated by commas where
#             the integers are the parent's feature indices.
#             This is necessary because v-build.pl::fetch_and_parse_cds_protein_feature_tables
#             has to set them as "type:GBSEP:coords" instead of
#             just the feature indices because later steps in
#             v-build.pl may remove some features, thus making the
#             feature indices invalid, and because v-annotate.pl
#             expects feature index integers not "type:GBSEP:coords".
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub integerize_parent_index_strings { 
  my $sub_name = "integerize_parent_index_strings";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if((defined $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"}) && 
       ($ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} ne "GBNULL")) { 
      my @parent_type_coords_A = split("!GBSEP!", $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"});
      my $new_parent_idx_str = "";
      foreach my $parent_type_coords_str (@parent_type_coords_A) { 
        my @el_A = split(":GBSEP:", $parent_type_coords_str);
        if(scalar(@el_A) != 2) { 
          ofile_FAIL("ERROR in $sub_name, unable to parse temporary parent_idx_str $parent_type_coords_str\n", 1, $FH_HR);
        }
        my ($parent_type, $parent_coords) = ($el_A[0], $el_A[1]);
        my $parent_ftr_idx = undef;
        # find parent idx in ftr_info_AHR, if it exists
        for(my $ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
          if($ftr_idx2 ne $ftr_idx) { # a feature can't be the parent of itself
            if(($ftr_info_AHR->[$ftr_idx2]{"type"}   eq $parent_type) && 
               ($ftr_info_AHR->[$ftr_idx2]{"coords"} eq $parent_coords)) { 
              if(defined $parent_ftr_idx) { 
                ofile_FAIL("ERROR in $sub_name, found two features that qualify as parents of feature $ftr_idx with type $parent_type and coords $parent_coords: $parent_ftr_idx and $ftr_idx2", 1, $FH_HR);
              }
              $parent_ftr_idx = $ftr_idx2;
            }
          }
        }
        if(! defined $parent_ftr_idx) { 
          ofile_FAIL("ERROR in $sub_name, unable to find a feature parent of $ftr_idx, expected a parent feature with type $parent_type and coords $parent_coords", 1, $FH_HR);
        }
        if($new_parent_idx_str ne "") { $new_parent_idx_str .= ","; }
        $new_parent_idx_str .= $parent_ftr_idx;
      }
      $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} = $new_parent_idx_str;
    }
  }

  return;
}

