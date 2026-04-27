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
use Mozilla::CA;

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
$execs_H{"cmfetch"}       = $env_vadr_infernal_dir . "/cmfetch";
$execs_H{"cmemit"}        = $env_vadr_infernal_dir . "/cmemit";
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
opt_Add("--profile",    "boolean", 0,          $g,  "--stk,--minfoin", "--infa,--inft,--gb,--ingb,--addminfo,--onlyurl", "build a model from an input multi-sequence alignment (--stk) and model info (--minfoin)", "build a model from an input multi-sequence alignment (--stk) and model info (--minfoin)", \%opt_HH, \@opt_order_A);
opt_Add("--stk",        "string",  undef,      $g,    undef, undef,       "read stockholm alignment from <s> (single-seq unless --profile)", "read stockholm alignment from <s> (single-seq unless --profile)", \%opt_HH, \@opt_order_A);
opt_Add("--minfoin",    "string",  undef,      $g, "--profile", undef,   "read model info file from <s> (for --profile mode)",           "read model info file from <s> (for --profile mode)", \%opt_HH, \@opt_order_A);
opt_Add("--infa",       "string",  undef,      $g,    undef, undef,       "read single sequence fasta file from <s>, don't fetch it",    "read single sequence fasta file from <s>, don't fetch it", \%opt_HH, \@opt_order_A);
opt_Add("--inft",       "string",  undef,      $g, "--inft", "--gb",      "read feature table file from <s>, don't fetch it",            "read feature table file from <s>, don't fetch it", \%opt_HH, \@opt_order_A);
opt_Add("--ftfetch1",   "boolean", 0,          $g,    undef, "--inft,--gb,--ftfetch2", "fetch feature table with efetch -format ft",      "fetch feature table with efetch -format ft", \%opt_HH, \@opt_order_A);
opt_Add("--ftfetch2",   "boolean", 0,          $g,    undef, "--inft,--gb,--ftfetch1", "fetch feature table with efetch -format gbc | xml2tbl", "fetch feature table with efetch -format gbc | xml2tbl", \%opt_HH, \@opt_order_A);
opt_Add("--gb",         "boolean", 0,          $g,    undef, undef,       "parse a genbank file, not a feature table file",              "parse a genbank file, not a feature table file", \%opt_HH, \@opt_order_A);
opt_Add("--ingb",       "string",  undef,      $g,   "--gb", undef,       "read genbank file from <s>, don't fetch it",                  "read genbank file from <s>, don't fetch it", \%opt_HH, \@opt_order_A);
opt_Add("--addminfo",   "string",  undef,      $g,    undef, undef,       "add feature info from model info file <s>",                   "add feature info from model info file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--forcelong",  "boolean", 0,          $g,    undef, undef,       "allow long models > 25Kb in length",                          "allow long models > 25Kb in length", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,          $g,    undef, undef,       "leave intermediate files on disk",                            "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling what feature types are stored in model info file\n[default set is: CDS,gene,mat_peptide,sig_peptide]";
#     option            type       default  group   requires incompat     preamble-output                                                      help-output    
opt_Add("--fall",       "boolean", 0,          $g,    undef,  undef,      "store info for all feature types (except those in --fskip)",        "store info for all feature types (except those in --fskip)", \%opt_HH, \@opt_order_A);
opt_Add("--fadd",       "string",  undef,      $g,    undef,"--fall",     "also store features types in comma separated string <s>",           "also store feature types in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--fskip",      "string",  undef,      $g,    undef,  undef,      "do not store info for feature types in comma separated string <s>",  "do not store info for feature types in comma separated string <s>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling what qualifiers are stored in model info file\n[default set is:product,gene,exception]";
#     option            type       default  group   requires incompat     preamble-output                                                             help-output    
opt_Add("--qall",       "boolean",  0,        $g,    undef,  undef,       "store info for all qualifiers (except those in --qskip)",                  "store info for all qualifiers (except those in --qskip)", \%opt_HH, \@opt_order_A);
opt_Add("--qadd",       "string",   undef,    $g,    undef,"--qall",      "also store info for qualifiers in comma separated string <s>",             "also store info for qualifiers in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--qftradd",    "string",   undef,    $g,"--qadd",    undef,      "--qadd <s2> only applies for feature types in comma separated string <s>", "--qadd <s2> only applies for feature types in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--qskip",      "string",   undef,    $g,    undef,   undef,      "do not store info for qualifiers in comma separated string <s>",           "do not store info for qualifiers in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--noaddgene", "boolean",  0,         $g,    undef,   undef,      "do not add gene qualifiers from gene features to overlapping features",     "do not add gene qualifiers from gene features to overlapping features", \%opt_HH, \@opt_order_A);
opt_Add("--nosplice",  "boolean",  0,         $g,    undef,   undef,      "do not check and add valid splice sites qualifiers for CDS",                "do not check and add valid splice sites qualifiers for CDS", \%opt_HH, \@opt_order_A);
opt_Add("--ssplice",   "boolean",  0,         $g,    undef,"--nosplice",  "exit if any noncanonical splice sites exist in any CDS",                    "exit if any noncanonical splice sites exist in any CDS",  \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for including additional model attributes";
#     option           type       default    group   requires    incompat   preamble-output                       help-output    
opt_Add("--group",     "string",  undef,        $g,  undef,         undef,  "specify model group is <s>",         "specify model group is <s>",    \%opt_HH, \@opt_order_A);
opt_Add("--subgroup",  "string",  undef,        $g,  "--group",     undef,  "specify model subgroup is <s>",      "specify model subgroup is <s>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling CDS translation step";
#     option          type       default    group   requires    incompat   preamble-output                                             help-output    
opt_Add("--ttbl",     "integer", 1,            $g,  undef,         undef,  "use NCBI translation table <n> to translate CDS",          "use NCBI translation table <n> to translate CDS", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling optional HMMER profile HMM building";
#     option          type       default    group   requires    incompat     preamble-output                                         help-output
opt_Add("--addhmm",   "boolean", 0,            $g,  undef,     "--profile", "build HMMER profile HMM db for CDS (off by default)",  "build HMMER profile HMM db for CDS (off by default)", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling cmbuild step";
#     option          type       default    group   requires    incompat             preamble-output                                             help-output    
opt_Add("--cmn",      "integer", undef,       $g,   undef, "--skipbuild,--cminfile", "set number of seqs for glocal fwd HMM calibration to <n>", "set number of seqs for glocal fwd HMM calibration to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--cmp7ml",   "boolean", 0,           $g,   undef, "--skipbuild,--cminfile", "set CM's filter p7 HMM as the ML p7 HMM",                  "set CM's filter p7 HMM as the ML p7 HMM",                  \%opt_HH, \@opt_order_A);
opt_Add("--cmere",    "real",    undef,       $g,   undef, "--skipbuild,--cminfile", "set CM relative entropy target to <x>",                    "set CM relative entropy target to <x>",                    \%opt_HH, \@opt_order_A);
opt_Add("--cmeset",   "real",    undef,       $g,   undef, "--skipbuild,--cminfile", "set CM eff seq # for CM to <x>",                           "set CM eff seq # for CM to <x>",                           \%opt_HH, \@opt_order_A);
opt_Add("--cmemaxseq","real",    undef,       $g,   undef, "--skipbuild,--cminfile", "set CM maximum allowed eff seq # for CM to <x>",           "set CM maximum alowed eff seq # for CM to <x>",            \%opt_HH, \@opt_order_A);
opt_Add("--cmnoh3pri","boolean", 0,           $g,   undef, "--skipbuild,--cminfile", "do not use --noh3pri option with cmbuild",                 "do not use --noh3pri option with cmbuild",                 \%opt_HH, \@opt_order_A);
opt_Add("--cminfile", "string",  undef,       $g,   undef, "--skipbuild",            "read cmbuild options from file <s>, one per line",         "read cmbuild options from file <s>, one per line",         \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for skipping stages";
#       option             type       default     group requires   incompat  preamble-output                                    help-output    
opt_Add("--skipbuild",     "boolean", 0,         $g,    undef,     undef,    "skip the cmbuild step",                           "skip the cmbuild and blastn db creation steps", \%opt_HH, \@opt_order_A);
opt_Add("--onlyurl",       "boolean", 0,         $g,    undef,"--stk,--ingb,--inft",  "output genbank file url for accession and exit",  "output genbank file url for accession and exit", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "optional output files";
#       option       type       default     group  requires     incompat  preamble-output                          help-output    
opt_Add("--ftrinfo",    "boolean", 0,         $g,    undef,     undef,    "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);
opt_Add("--sgminfo",    "boolean", 0,         $g,    undef,     undef,    "output internal segment information",   "create file with internal segment information", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "other expert options";
#       option       type          default     group  requires incompat      preamble-output                                              help-output    
opt_Add("--execname",   "string",  undef,         $g,    undef, undef,       "define executable name of this script as <s>",              "define executable name of this script as <s>", \%opt_HH, \@opt_order_A);        
opt_Add("--nosig2mat",  "boolean", 0,             $g,    undef, undef,       "do not treat sig_peptide as mat_peptide",                   "do not treat sig_peptide as mat_peptide", \%opt_HH, \@opt_order_A);        
opt_Add("--intlen",     "integer", 40,            $g,    undef,"--nosplice", "set min length of intron to check for splice sites to <n>", "set min length of intron to check for splice sites to <n>", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
        'profile'      => \$GetOptions_H{"--profile"},
                'stk=s'        => \$GetOptions_H{"--stk"},
        'minfoin=s'    => \$GetOptions_H{"--minfoin"},
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
                'nosplice'     => \$GetOptions_H{"--nosplice"},
                'ssplice'      => \$GetOptions_H{"--ssplice"},
# options for including additional model attributes
                'group=s'      => \$GetOptions_H{"--group"},
                'subgroup=s'   => \$GetOptions_H{"--subgroup"},
# options for controlling CDS translation step
                'ttbl=s'       => \$GetOptions_H{"--ttbl"},
# options for controlling optional HMMER profile HMM building
                'addhmm'       => \$GetOptions_H{"--addhmm"},
# options for controlling cmbuild step
                'cmn=s'        => \$GetOptions_H{"--cmn"},
                'cmp7ml'       => \$GetOptions_H{"--cmp7ml"},
                'cmere=s'      => \$GetOptions_H{"--cmere"},
                'cmeset=s'     => \$GetOptions_H{"--cmeset"},
                'cmemaxseq=s'  => \$GetOptions_H{"--cmemaxseq"},
                'cmnoh3pr'     => \$GetOptions_H{"--cmnoh3pri"},
                'cminfile=s'   => \$GetOptions_H{"--cminfile"},
# options for skipping stages
                'skipbuild'    => \$GetOptions_H{"--skipbuild"},
                'onlyurl'      => \$GetOptions_H{"--onlyurl"},
# optional output files
                'sgminfo'      => \$GetOptions_H{"--sgminfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"},
# other expert options
                'execname=s'   => \$GetOptions_H{"--execname"},
                'nosig2mat'    => \$GetOptions_H{"--nosig2mat"},
                'intlen=s'     => \$GetOptions_H{"--intlen"});

my $total_seconds = -1 * ofile_SecondsSinceEpoch(); # by multiplying by -1, we can just add another ofile_SecondsSinceEpoch call at end to get total time
my $execname_opt  = $GetOptions_H{"--execname"};
my $executable    = (defined $execname_opt) ? $execname_opt : "v-build.pl";
my $usage         = "Usage: $executable [-options] <accession> <path to output directory to create>\n";
my $synopsis      = "$executable :: build homology model for feature annotation";
my $date          = scalar localtime();
my $version       = "1.7";
my $releasedate   = "Sep 2025";
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

if($mdl_name =~ /[\(\)]/) { 
  die "ERROR, accession cannot contain '(' or ')'";
}

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

if($dir =~ m/\/$/) { 
  $dir =~ s/\/$//; # remove final '/' if it exists
} 
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
my $out_root         = $dir . "/" . $dir_tail . ".vadr";

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

########################################
# Determine if we are in --profile mode
########################################
my $do_profile   = opt_Get("--profile", \%opt_HH);
my $minfoin_file = opt_Get("--minfoin", \%opt_HH);

##############################################
# Variables that exist for all modes
##############################################
my $fa_file = $out_root . ".fa";
my %seq_H = ();
my $mdl_name_ver = undef;
my $maxlen = 25000;
my $mdllen = undef;
my $ft_file = undef;
my $gb_file = undef;
my %ftr_info_HAH = (); # the feature info
my %minfoin_mdl_info_H = (); # model info from --minfoin for $mdl_name (only in --profile mode)

###################################################
# Fetch the fasta file (if necessary) and parse it
###################################################
if(! $do_profile) {
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
  $mdl_name_ver = $fetched_seq_A[0];
  # make sure it's the right sequence
  if($mdl_name_ver =~ /(\S+)\.\d+/) { 
    if($1 ne $mdl_name) { 
      ofile_FAIL("ERROR did not fetch correct sequence from fasta file $fa_file (expected accession.version starting with $mdl_name, got $mdl_name_ver)\n", 1, $FH_HR);
    }
  }
  else {
    ofile_FAIL("ERROR did not fetch correct sequence from fasta file $fa_file (expected accession.version starting with $mdl_name, got $mdl_name_ver)\n", 1, $FH_HR);
  }
  $mdllen = length($seq_H{$mdl_name_ver});
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}
else {
  # --profile mode: build from an input Stockholm alignment and an input .minfo file.
  $start_secs = ofile_OutputProgressPrior("Processing --profile inputs", $progress_w, $log_FH, *STDOUT);

  if(! defined opt_Get("--stk", \%opt_HH)) {
    ofile_FAIL("ERROR, --profile requires --stk <s>", 1, $FH_HR);
  }
  if(! defined $minfoin_file) {
    ofile_FAIL("ERROR, --profile requires --minfoin <s>", 1, $FH_HR);
  }
  utl_FileValidateExistsAndNonEmpty($minfoin_file, "--minfoin file", undef, 1, $FH_HR);
  utl_FileValidateExistsAndNonEmpty(opt_Get("--stk", \%opt_HH), "--stk file", undef, 1, $FH_HR);

  my @minfo_mdl_info_AH  = ();
  my %minfo_ftr_info_HAH = ();
  my @reqd_mdl_keys_A = ("name");
  my @reqd_ftr_keys_A = ("type", "coords");
  vdr_ModelInfoFileParse($minfoin_file, \@reqd_mdl_keys_A, \@reqd_ftr_keys_A, \@minfo_mdl_info_AH, \%minfo_ftr_info_HAH, $FH_HR);

  if(! defined $minfo_ftr_info_HAH{$mdl_name}) {
    ofile_FAIL("ERROR, --minfoin file $minfoin_file must include model $mdl_name, but it does not", 1, $FH_HR);
  }
  @{$ftr_info_HAH{$mdl_name}} = @{$minfo_ftr_info_HAH{$mdl_name}};
  # also capture any model-level info for this model, so we can preserve it
  for(my $i = 0; $i < scalar(@minfo_mdl_info_AH); $i++) {
    if((defined $minfo_mdl_info_AH[$i]{"name"}) && ($minfo_mdl_info_AH[$i]{"name"} eq $mdl_name)) {
      %minfoin_mdl_info_H = %{$minfo_mdl_info_AH[$i]};
      last;
    }
  }
  if(! defined $minfoin_mdl_info_H{"name"}) {
    ofile_FAIL("ERROR, --minfoin file $minfoin_file must include MODEL line for model $mdl_name, but it does not", 1, $FH_HR);
  }
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

#######################################################################
# Verify our sequence is not longer than our hard-coded maximum of 25Kb 
#######################################################################
# allow any length if --forcelong
if((defined $mdllen) && (! opt_Get("--forcelong", \%opt_HH))) { 
  if($mdllen > $maxlen) { 
    ofile_FAIL("ERROR, model length ($mdllen) exceeds maximum allowed length of $maxlen.\nYou can use --forcelong to bypass this at your own risk.\nUse of VADR on models > 25Kb is not recommended.\nModel building will be very slow and\ndownstream v-annotate.pl annotation will have large memory requirements.", 1, $FH_HR);
  }
}

###################################################################
# Fetch the feature table (ft) or GenBank (gb) file (if necessary)
# and parse it.
###################################################################
$ft_file = undef;
$gb_file = undef;
if($do_profile) {
  # in --profile mode, feature info already came from --minfoin
}
elsif(! opt_IsUsed("--gb", \%opt_HH)) { 
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
process_add_and_skip_options("CDS,gene,mat_peptide,sig_peptide", "--fadd", "--fskip", undef, \%fdf_H, \%fadd_H, \%fskip_H, undef, \%opt_HH, $FH_HR);

# determine what qualifiers we will store based on cmdline options
# --qall is incompatible with --qadd
my %qdf_H      = (); # default qualifiers to keep
my %qadd_H     = (); # qualifiers to add
my %qskip_H    = (); # qualifiers to skip
my %qftr_add_H = (); # if --qftradd, subset of features to add qualifiers in --qadd option for
process_add_and_skip_options("type,coords,location,product,gene,exception,parent_idx_str,5p_trunc,3p_trunc,alternative_ftr_set,alternative_ftr_set_subn", "--qadd", "--qskip", "--qftradd", \%qdf_H, \%qadd_H, \%qskip_H, \%qftr_add_H, \%opt_HH, $FH_HR);
# 'ribosomal_slippage' is intentionally absent from the default qualifier list:
# it is converted to 'exception:"ribosomal slippage"' below (before the qualifier
# filter), so callers see only the canonical 'exception' form in the minfo.

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

# Promote 'ribosomal_slippage' qualifiers to 'exception:"ribosomal slippage"'.
# Run this BEFORE the qualifier filter (below) so ribosomal_slippage is still
# present, and run it for both the --gb and feature-table parsing paths:
#   * GenBank format includes /ribosomal_slippage but not /exception, so the
#     conversion is required.
#   * NCBI's feature-table format normally includes both, in which case the
#     conversion is a no-op for the exception value, and the duplicate
#     ribosomal_slippage key is dropped by the filter below.
#   * Custom feature tables that include only /ribosomal_slippage (without
#     /exception) are also handled correctly here.
# Only 'exception:ribosomal slippage' qualifier/values are desired in the
# output minfo; the bare 'ribosomal_slippage' is dropped unless preserved
# explicitly via --qall or --qadd ribosomal_slippage.
for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$mdl_name}}); $ftr_idx++) {
  if((defined $ftr_info_HAH{$mdl_name}[$ftr_idx]) &&
     (defined $ftr_info_HAH{$mdl_name}[$ftr_idx]{"ribosomal_slippage"})) {
    my $exc_existing = $ftr_info_HAH{$mdl_name}[$ftr_idx]{"exception"};
    if((defined $exc_existing) && ($exc_existing =~ m/(?:^|:GBSEP:)ribosomal slippage(?::GBSEP:|$)/)) {
      ; # already present (e.g., NCBI FT carries both qualifiers); no-op
    }
    elsif(defined $exc_existing) {
      $ftr_info_HAH{$mdl_name}[$ftr_idx]{"exception"} = $exc_existing . ":GBSEP:" . "ribosomal slippage";
    }
    else {
      $ftr_info_HAH{$mdl_name}[$ftr_idx]{"exception"} = "ribosomal slippage";
    }
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
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

####################################################################################
# Change 'type' value for sig_peptide features to 'mat_peptide' (unless --nosig2mat)
####################################################################################
if(! opt_Get("--nosig2mat", \%opt_HH)) { 
  for($ftr_idx = 0; $ftr_idx < scalar(@{$ftr_info_HAH{$mdl_name}}); $ftr_idx++) { 
    if((defined $ftr_info_HAH{$mdl_name}[$ftr_idx]) && 
       (defined $ftr_info_HAH{$mdl_name}[$ftr_idx]{"type"}) && 
       ($ftr_info_HAH{$mdl_name}[$ftr_idx]{"type"} eq "sig_peptide")) { 
      $ftr_info_HAH{$mdl_name}[$ftr_idx]{"type"}     = "mat_peptide";
      $ftr_info_HAH{$mdl_name}[$ftr_idx]{"out_type"} = "sig_peptide";
    }
  }
}

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
my $stk_has_rf = 0;
my $in_stk_file = opt_Get("--stk", \%opt_HH);
if(defined $in_stk_file) { 
  $start_secs = ofile_OutputProgressPrior("Validating input Stockholm file", $progress_w, $log_FH, *STDOUT);

  if($do_profile) {
    my $msa = Bio::Easel::MSA->new({ fileLocation => $in_stk_file, isDna => 1});
    if(! $msa->has_rf) {
      ofile_FAIL("ERROR, --profile requires that the --stk Stockholm alignment includes RF annotation (#=GC RF)", 1, $FH_HR);
    }
    my $rf = $msa->get_rf;
    my $clen = 0;
    foreach my $rfchar (split(//, $rf)) {
      if($rfchar !~ /[\-\_\.\~]/) { $clen++; }
    }
    if($clen == 0) {
      ofile_FAIL("ERROR, --profile read RF annotation but computed 0 consensus columns; is RF all gaps?", 1, $FH_HR);
    }
    $mdllen = $clen;
    $stk_has_ss = $msa->has_ss_cons;
    $stk_has_rf = 1;
    undef $msa;
  }
  else {
    $stk_has_ss = stockholm_validate_single_sequence_input($in_stk_file, $seq_H{$mdl_name_ver}, \%opt_HH, $FH_HR);
  }

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

###########################
# Build the blastn database
###########################
# In --profile mode, the blastn db is built AFTER the CM (further below): the
# db sequence is generated from the final CM via 'cmemit -c' so the blastn
# subject sequence matches the consensus the CM emits at each match state.
# We also rewrite the RF annotation of the .stk file to match this consensus,
# so that v-annotate.pl's NN-based classification (which reads RF from the
# .stk) is consistent with the blastn db sequence.
my $tmp_blastn_fa_file = $out_root . ".fa.tmp";
my $blastn_fa_file     = $out_root . ".fa";

if(! $do_profile) {
  $start_secs = ofile_OutputProgressPrior("Building BLAST nucleotide database ", $progress_w, $log_FH, *STDOUT);

  sqf_EslReformatRun($execs_H{"esl-reformat"}, "-d -u", $fa_file, $tmp_blastn_fa_file, "fasta", "fasta", \%opt_HH, $FH_HR);
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "blastn-fa",  $blastn_fa_file, 1, 1, "nucleotide blastn db fasta sequence file for $mdl_name");
  my $nt_sqfile    = Bio::Easel::SqFile->new({ fileLocation => $tmp_blastn_fa_file });
  my $tmp_seq_name = $nt_sqfile->fetch_seq_name_given_ssi_number(0);
  my $tmp_sqstring = $nt_sqfile->fetch_seq_to_sqstring($tmp_seq_name);
  printf { $ofile_info_HH{"FH"}{"blastn-fa"} } ">" . $mdl_name . "\n" . seq_SqstringAddNewlines($tmp_sqstring, 60);
  close $ofile_info_HH{"FH"}{"blastn-fa"};

  # run makeblastdb
  sqf_BlastDbCreate($execs_H{"makeblastdb"}, "nucl", $blastn_fa_file, \%opt_HH, $FH_HR);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nhr", $blastn_fa_file . ".nhr", 1, 1, "BLAST db .nhr file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nin", $blastn_fa_file . ".nin", 1, 1, "BLAST db .nin file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nsq", $blastn_fa_file . ".nsq", 1, 1, "BLAST db .nsq file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-ndb", $blastn_fa_file . ".ndb", 1, 1, "BLAST db .ndb file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-not", $blastn_fa_file . ".not", 1, 1, "BLAST db .not file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-ntf", $blastn_fa_file . ".ntf", 1, 1, "BLAST db .ntf file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nto", $blastn_fa_file . ".nto", 1, 1, "BLAST db .nto file for $mdl_name");
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-njs", $blastn_fa_file . ".njs", 1, 1, "BLAST db .njs file for $mdl_name");

  if(! opt_Get("--keep", \%opt_HH)) {
    utl_FileRemoveUsingSystemRm($tmp_blastn_fa_file, "v-build.pl main", \%opt_HH, $FH_HR);
    if(-e $tmp_blastn_fa_file . ".ssi") {
      utl_FileRemoveUsingSystemRm($tmp_blastn_fa_file . ".ssi", "v-build.pl main", \%opt_HH, $FH_HR);
    }
  }
  # index the new file
  my $sfetch_cmd = $execs_H{"esl-sfetch"} . " --index $blastn_fa_file > /dev/null";
  utl_RunCommand($sfetch_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);

  ofile_OutputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);
} # end of 'if(! $do_profile)' block for non-profile blastn db creation

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
# - check for valid splice sites
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
  if($do_profile) {
    profile_CdsFetchStockholmToFasta($ofile_info_HH{"FH"}{"cdsfasta"}, $stk_file, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  }
  else {
    vdr_CdsFetchStockholmToFasta($ofile_info_HH{"FH"}{"cdsfasta"}, $stk_file, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  }
  close $ofile_info_HH{"FH"}{"cdsfasta"};
  
  $protein_fa_file = $out_root . ".protein.fa";
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "proteinfasta", $protein_fa_file, 1, 1, "fasta sequence file for translated CDS from $mdl_name");
  
  if($do_profile) {
    # For profile mode: translate to temp file, then rename headers
    my $tmp_protein_fa_file = $out_root . ".protein.tmp.fa";
    open(my $tmp_protein_FH, ">", $tmp_protein_fa_file) || ofile_FileOpenFailure($tmp_protein_fa_file, "v-build.pl::main()", $!, "writing", $FH_HR);
    sqf_EslTranslateCdsToFastaFile($tmp_protein_FH, $execs_H{"esl-translate"}, $cds_fa_file, 
                                   $out_root, \@{$ftr_info_HAH{$mdl_name}}, \%opt_HH, $FH_HR);
    close($tmp_protein_FH);
    
    # Rename protein headers from seqname/seqcoords to seqname:seqcoords/refcoords
    profile_RenameProteinHeaders($tmp_protein_fa_file, $ofile_info_HH{"FH"}{"proteinfasta"}, $cds_fa_file, $FH_HR);
    
    if(! opt_Get("--keep", \%opt_HH)) {
      utl_FileRemoveUsingSystemRm($tmp_protein_fa_file, "v-build.pl::main()", \%opt_HH, $FH_HR);
    }
  }
  else {
    sqf_EslTranslateCdsToFastaFile($ofile_info_HH{"FH"}{"proteinfasta"}, $execs_H{"esl-translate"}, $cds_fa_file, 
                                   $out_root, \@{$ftr_info_HAH{$mdl_name}}, \%opt_HH, $FH_HR);
  }
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

  # build hmmer db (optional), need to build one HMM per CDS and concatenate them
  if(opt_Get("--addhmm", \%opt_HH)) {
    $start_secs = ofile_OutputProgressPrior("Building HMMER protein database ", $progress_w, $log_FH, *STDOUT);

  # run esl-seqstat and parse it
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
      if($hmm_name =~ /[\(\)]/) { 
        ofile_FAIL("ERROR, illegal sequence name in $protein_fa_file, sequence names can't have ')' or '(' in them", 1, $FH_HR);
      }
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

  # check splice sites and add canon_splice_sites:"1" if all are valid
  if(! opt_Get("--nosplice", \%opt_HH)) { 
    $start_secs = ofile_OutputProgressPrior("Checking intron splice sites, if any", $progress_w, $log_FH, *STDOUT);
    check_and_add_cds_splice_sites($stk_file, \@{$ftr_info_HAH{$mdl_name}}, \@sgm_info_AH, \%opt_HH, $FH_HR);
    ofile_OutputProgressComplete($start_secs, undef,  $log_FH, *STDOUT);
  }
}

##############
# Build the CM
##############
my $cm_file = undef;
if(! opt_Get("--skipbuild", \%opt_HH)) { 
  my $cmbuild_str = undef;
  my $clen_times_cmn = $mdllen * 200;
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

#  my $cmbuild_occ_file = $out_root . ".cmbuild.occ";
#  my $cmbuild_cp9occ_file = $out_root . ".cmbuild.cp9occ";
#  my $cmbuild_fp7occ_file = $out_root . ".cmbuild.fp7occ";

#  my $cmbuild_opts = "-n $mdl_name --verbose --occfile $cmbuild_occ_file --cp9occfile $cmbuild_cp9occ_file --fp7occfile $cmbuild_fp7occ_file ";
  my $cmbuild_opts = "-n $mdl_name --verbose ";
  if(opt_IsUsed("--cminfile",  \%opt_HH))   { 
    my @cminfile_A = ();
    utl_FileLinesToArray(opt_Get("--cminfile", \%opt_HH), 1, \@cminfile_A, $FH_HR);
    foreach my $optline (@cminfile_A) { 
      chomp $optline;
      $cmbuild_opts .= " " . $optline . " ";
    }
  }
  else { 
    if((! defined $stk_has_ss) || (! $stk_has_ss)) { $cmbuild_opts .= " --noss"; }
    if(($do_profile) && ($stk_has_rf)) { $cmbuild_opts .= " --hand"; }
    if(opt_IsUsed("--cmn",       \%opt_HH))   { $cmbuild_opts .= " --EgfN "    . opt_Get("--cmn", \%opt_HH); }
    if(opt_IsUsed("--cmp7ml",    \%opt_HH))   { $cmbuild_opts .= " --p7ml"; }
    if(opt_IsUsed("--cmere",     \%opt_HH))   { $cmbuild_opts .= " --ere "     . opt_Get("--cmere", \%opt_HH); }
    if(opt_IsUsed("--cmeset",    \%opt_HH))   { $cmbuild_opts .= " --eset "    . opt_Get("--cmeset", \%opt_HH); }
    if(opt_IsUsed("--cmemaxseq", \%opt_HH))   { $cmbuild_opts .= " --emaxseq " . opt_Get("--cmemaxseq", \%opt_HH); }
    if(! opt_IsUsed("--cmnoh3pri", \%opt_HH)) { $cmbuild_opts .= " --noh3pri"; }
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

  ###############################################################
  # In --profile mode: build the blastn db from the final CM and
  # rewrite the RF annotation of the .stk to match the consensus
  ###############################################################
  if($do_profile) {
    $start_secs = ofile_OutputProgressPrior("Building BLAST nucleotide database from CM consensus", $progress_w, $log_FH, *STDOUT);

    # Run cmemit -c to get the consensus sequence emitted by the final CM.
    # The output sequence has length equal to the number of consensus (match)
    # columns of the model and characters equal to the most-likely emission
    # at each match state. This becomes the single subject sequence in the
    # blastn classification db (named after the model) AND is used to rewrite
    # the RF annotation of the .stk so that v-annotate.pl's NN-based
    # classification (which reads RF from the .stk) is consistent with it.
    my $cmemit_fa_file = $out_root . ".cmemit-c.fa";
    my $cmemit_cmd     = $execs_H{"cmemit"} . " -c $cm_file > $cmemit_fa_file";
    utl_RunCommand($cmemit_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);

    # read the single consensus sequence from the cmemit output
    my $cmemit_sqfile  = Bio::Easel::SqFile->new({ fileLocation => $cmemit_fa_file });
    my $cmemit_seqname = $cmemit_sqfile->fetch_seq_name_given_ssi_number(0);
    my $consensus      = $cmemit_sqfile->fetch_seq_to_sqstring($cmemit_seqname);
    undef $cmemit_sqfile;
    seq_SqstringCapitalize(\$consensus);
    seq_SqstringDnaize(\$consensus);
    # strip any whitespace/newlines just in case
    $consensus =~ s/\s+//g;
    if(length($consensus) != $mdllen) {
      ofile_FAIL(sprintf("ERROR, --profile: cmemit -c consensus length (%d) does not match model consensus length (%d)",
                         length($consensus), $mdllen), 1, $FH_HR);
    }

    # write the blastn db fasta: single sequence named after the model
    ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "blastn-fa", $blastn_fa_file, 1, 1, "nucleotide blastn db fasta sequence file for $mdl_name");
    printf { $ofile_info_HH{"FH"}{"blastn-fa"} } ">" . $mdl_name . "\n" . seq_SqstringAddNewlines($consensus, 60);
    close $ofile_info_HH{"FH"}{"blastn-fa"};

    # makeblastdb on the new file
    sqf_BlastDbCreate($execs_H{"makeblastdb"}, "nucl", $blastn_fa_file, \%opt_HH, $FH_HR);
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nhr", $blastn_fa_file . ".nhr", 1, 1, "BLAST db .nhr file for $mdl_name");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nin", $blastn_fa_file . ".nin", 1, 1, "BLAST db .nin file for $mdl_name");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nsq", $blastn_fa_file . ".nsq", 1, 1, "BLAST db .nsq file for $mdl_name");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-ndb", $blastn_fa_file . ".ndb", 1, 1, "BLAST db .ndb file for $mdl_name");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-not", $blastn_fa_file . ".not", 1, 1, "BLAST db .not file for $mdl_name");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-ntf", $blastn_fa_file . ".ntf", 1, 1, "BLAST db .ntf file for $mdl_name");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-nto", $blastn_fa_file . ".nto", 1, 1, "BLAST db .nto file for $mdl_name");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastdb-njs", $blastn_fa_file . ".njs", 1, 1, "BLAST db .njs file for $mdl_name");

    # esl-sfetch index
    my $sfetch_cmd = $execs_H{"esl-sfetch"} . " --index $blastn_fa_file > /dev/null";
    utl_RunCommand($sfetch_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);

    # rewrite the RF annotation of the .stk file so that the consensus
    # characters in RF match the cmemit -c consensus, position-for-position
    # (non-consensus columns -- gap chars in RF -- are left untouched)
    my $rf_msa = Bio::Easel::MSA->new({ fileLocation => $stk_file, isDna => 1 });
    my $old_rf = $rf_msa->get_rf;
    my @old_rf_A = split(//, $old_rf);
    my @cons_A   = split(//, $consensus);
    my $cidx = 0;
    my $new_rf = "";
    foreach my $rfchar (@old_rf_A) {
      if($rfchar =~ /[\-\_\.\~]/) {
        $new_rf .= $rfchar;
      }
      else {
        if($cidx >= scalar(@cons_A)) {
          ofile_FAIL("ERROR, --profile: ran out of consensus characters while rewriting RF annotation", 1, $FH_HR);
        }
        $new_rf .= $cons_A[$cidx];
        $cidx++;
      }
    }
    if($cidx != scalar(@cons_A)) {
      ofile_FAIL(sprintf("ERROR, --profile: RF rewrite consumed %d of %d consensus characters", $cidx, scalar(@cons_A)), 1, $FH_HR);
    }
    $rf_msa->set_rf($new_rf);
    # overwrite the .stk file in place with the updated RF annotation
    $rf_msa->write_msa($stk_file, "stockholm");
    undef $rf_msa;

    if(! opt_Get("--keep", \%opt_HH)) {
      utl_FileRemoveUsingSystemRm($cmemit_fa_file, "v-build.pl main", \%opt_HH, $FH_HR);
      if(-e $cmemit_fa_file . ".ssi") {
        utl_FileRemoveUsingSystemRm($cmemit_fa_file . ".ssi", "v-build.pl main", \%opt_HH, $FH_HR);
      }
    }
    else {
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "cmemit-c-fa", $cmemit_fa_file, 1, 1, "cmemit -c consensus fasta for $mdl_name");
    }

    ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
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
$mdl_info_AH[0]{"length"} = $mdllen;
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
elsif($do_profile) {
  # preserve group/subgroup from --minfoin unless overridden on cmdline
  if((defined $minfoin_mdl_info_H{"group"}) && (! opt_IsUsed("--group", \%opt_HH))) {
    $mdl_info_AH[0]{"group"} = $minfoin_mdl_info_H{"group"};
    if(defined $minfoin_mdl_info_H{"subgroup"}) {
      $mdl_info_AH[0]{"subgroup"} = $minfoin_mdl_info_H{"subgroup"};
    }
  }
  # preserve transl_table from --minfoin if present and not overridden
  if((defined $minfoin_mdl_info_H{"transl_table"}) && (! opt_IsUsed("--ttbl", \%opt_HH))) {
    $mdl_info_AH[0]{"transl_table"} = $minfoin_mdl_info_H{"transl_table"};
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
# Subroutine: profile_CdsFetchStockholmToFasta()
# Incept:     EPN?, Adapted for profile mode by GitHub Copilot
#
# Purpose:    Like vdr_CdsFetchStockholmToFasta(), but in --profile mode
#             we must name each CDS (and therefore the translated protein)
#             with BOTH:
#               (1) the dealigned (sequence) coordinates for that CDS in
#                   that specific sequence, derived from the alignment, and
#               (2) the model (reference) coordinates from the input .minfo.
#
#             Output FASTA header format per CDS per sequence:
#               <seqname>:<seqcoords>/<refcoords>
#             Example:
#               myseq:742..7323:+/745..7398:+
#
# Arguments:
#   $out_FH:         output file handle
#   $stk_file:       stockholm file with aligned full length sequences
#   $ftr_info_AHR:   REF to the feature info, pre-filled (model coords)
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns: void
#################################################################
sub profile_CdsFetchStockholmToFasta {
  my $sub_name = "profile_CdsFetchStockholmToFasta";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($out_FH, $stk_file, $ftr_info_AHR, $FH_HR) = (@_);

  my $msa = Bio::Easel::MSA->new({ fileLocation => $stk_file, isDna => 1});
  if(! $msa->has_rf) {
    ofile_FAIL("ERROR in $sub_name, --profile requires RF annotation in $stk_file", 1, $FH_HR);
  }

  # precompute start, stop, strand, for all features, so we don't have to redo this for each seq
  my @sgm_start_AA  = ();
  my @sgm_stop_AA   = ();
  my @sgm_strand_AA = ();
  vdr_FeatureInfoStartStopStrandArrays($ftr_info_AHR, \@sgm_start_AA, \@sgm_stop_AA, \@sgm_strand_AA, $FH_HR);

  my $nftr = scalar(@{$ftr_info_AHR});
  my $nseq = $msa->nseq;

  # Build map of alternative feature sets: for each set name, collect
  # all member ftr_idx values in order (first = primary)
  my %alt_set_members = ();  # set_name -> [ftr_idx1, ftr_idx2, ...]
  for(my $fi = 0; $fi < $nftr; $fi++) {
    if($ftr_info_AHR->[$fi]{"type"} eq "CDS" &&
       defined $ftr_info_AHR->[$fi]{"alternative_ftr_set"} &&
       $ftr_info_AHR->[$fi]{"alternative_ftr_set"} ne "") {
      my $afset = $ftr_info_AHR->[$fi]{"alternative_ftr_set"};
      if(! exists $alt_set_members{$afset}) {
        $alt_set_members{$afset} = [];
      }
      push(@{$alt_set_members{$afset}}, $fi);
    }
  }

  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) {
    my $sqname = $msa->get_sqname($seq_idx);
    my $aligned_sqstring = $msa->get_sqstring_aligned($seq_idx);
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) {
      if($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") {
        # For CDS with alternative_ftr_set: skip non-primary members
        # (the primary handles all alternatives for each sequence).
        # For the primary: try extracting with primary coords first,
        # then try each alternative if the primary fails validation.
        if(defined $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"} &&
           $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"} ne "") {
          my $afset = $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"};
          # Skip non-primary members
          if(exists $alt_set_members{$afset} && $alt_set_members{$afset}[0] != $ftr_idx) {
            next;
          }
          # This is the primary — try it and each alternative
          my @try_ftr_idxs = @{$alt_set_members{$afset}};
          my $success = 0;
          foreach my $try_fi (@try_ftr_idxs) {
            my ($cds_str, $header) = profile_ExtractOneCds($msa, $seq_idx, $sqname, $aligned_sqstring,
                                                           $try_fi, $ftr_info_AHR,
                                                           \@sgm_start_AA, \@sgm_stop_AA, \@sgm_strand_AA,
                                                           $FH_HR);
            if(defined $cds_str && profile_ValidateCdsIsComplete($cds_str)) {
              print $out_FH ">" . $header . "\n";
              print $out_FH seq_SqstringAddNewlines($cds_str, 60);
              $success = 1;
              last;
            }
          }
          if(! $success) {
            ofile_OutputString($FH_HR->{"log"}, 1,
              sprintf("# WARNING: profile_CdsFetchStockholmToFasta: no alternative CDS translated for %s in set %s, skipping\n",
                      $sqname, $afset));
          }
          next;  # skip the normal extraction below
        }
        # Non-alternative CDS: extract normally (fail on error)
        my $cds_sqstring = "";
        my @seq_sgm_coords_A = ();
        my $total_rf_offset_5p = 0; # cumulative 5' RF offset (nt to trim from 5' end)
        my $total_rf_offset_3p = 0; # cumulative 3' RF offset (nt to trim from 3' end)
        foreach(my $sgm_idx = 0; $sgm_idx < scalar(@{$sgm_start_AA[$ftr_idx]}); $sgm_idx++) {
          my $rfstart = $sgm_start_AA[$ftr_idx][$sgm_idx];
          my $rfstop  = $sgm_stop_AA[$ftr_idx][$sgm_idx];
          my $astart  = $msa->rfpos_to_aligned_pos($rfstart);
          my $astop   = $msa->rfpos_to_aligned_pos($rfstop);
          if($astart > $astop) { utl_Swap(\$astart, \$astop); }

          # derive dealigned (sequence) coords for this segment
          my ($ua_first, $ua_last) = profile_FirstAndLastUngappedPositionsInAlignedRange($aligned_sqstring, $astart, $astop);
          if(! defined $ua_first) {
            ofile_FAIL("ERROR in $sub_name, unable to determine dealigned coords for $sqname segment $astart..$astop (all gaps?)", 1, $FH_HR);
          }

          # fetch the unaligned segment sequence (in alignment coords)
          my $sgm_sqstring = $msa->get_sqstring_unaligned_and_truncated($seq_idx, $astart, $astop);
          if($sgm_strand_AA[$ftr_idx][$sgm_idx] eq "-") {
            seq_SqstringReverseComplement(\$sgm_sqstring);
          }

          # Compute RF-based offset: how many nucleotides at 5' end of this segment
          # correspond to gapped RF positions (i.e., before the first RF position that
          # is ungapped in this sequence).
          my ($rf_offset_5p, $rf_offset_3p) = profile_ComputeRfOffsets($aligned_sqstring, $astart, $astop, $rfstart, $rfstop, $msa, $sgm_strand_AA[$ftr_idx][$sgm_idx]);
          $total_rf_offset_5p += $rf_offset_5p;
          $total_rf_offset_3p += $rf_offset_3p;

          $cds_sqstring .= $sgm_sqstring;
        }

        # Don't trim the CDS sequence - pass it intact to esl-translate
        # esl-translate will find the correct ORF based on the coordinates we specify
        my $cds_len = length($cds_sqstring);
        if($total_rf_offset_5p + $total_rf_offset_3p >= $cds_len) {
          ofile_FAIL("ERROR in $sub_name, CDS for $sqname is entirely gapped in RF positions (5p_offset=$total_rf_offset_5p, 3p_offset=$total_rf_offset_3p, len=$cds_len)", 1, $FH_HR);
        }

        # The final CDS is the full extracted sequence
        my $final_cds = $cds_sqstring;
        my $final_len = length($final_cds);

        # Compute sequence coordinates based on the full (untrimmed) CDS
        my @final_seq_sgm_coords_A = ();

        foreach(my $sgm_idx = 0; $sgm_idx < scalar(@{$sgm_start_AA[$ftr_idx]}); $sgm_idx++) {
          my $rfstart = $sgm_start_AA[$ftr_idx][$sgm_idx];
          my $rfstop  = $sgm_stop_AA[$ftr_idx][$sgm_idx];
          my $astart  = $msa->rfpos_to_aligned_pos($rfstart);
          my $astop   = $msa->rfpos_to_aligned_pos($rfstop);
          if($astart > $astop) { utl_Swap(\$astart, \$astop); }

          my ($ua_first, $ua_last) = profile_FirstAndLastUngappedPositionsInAlignedRange($aligned_sqstring, $astart, $astop);

          # Use the full segment coordinates without trimming
          if($ua_first <= $ua_last) {
            my $seg_strand = $sgm_strand_AA[$ftr_idx][$sgm_idx];
            # Positive strand: low..high:+ (5'->3' = low->high)
            # Negative strand: high..low:- (5'->3' = high->low, VADR convention)
            if($seg_strand eq "-") {
              push(@final_seq_sgm_coords_A, $ua_last . ".." . $ua_first . ":" . $seg_strand);
            }
            else {
              push(@final_seq_sgm_coords_A, $ua_first . ".." . $ua_last . ":" . $seg_strand);
            }
          }
        }

        my $seq_coords_str = join(",", @final_seq_sgm_coords_A);
        my $ref_coords_str = $ftr_info_AHR->[$ftr_idx]{"coords"};
        
        # Determine truncation status:
        # 5': Check if (a) reference CDS was already 5'-truncated, OR
        #              (b) this sequence is missing RF positions at the 5' end (has RF offset > 0)
        # 3': Check if this sequence has ungapped residues at the last 3 RF positions
        #     (which should encode the stop codon). If any of those positions are gapped,
        #     then this sequence doesn't have a stop codon and should be 3'-truncated.
        my $ref_is_trunc5p = ($ref_coords_str =~ /</) ? 1 : 0;
        # Check if the biological 5' RF position of this CDS is gapped in this sequence.
        # This catches truncation even when total_rf_offset_5p==0 (i.e., the sequence starts
        # downstream of rfstart but happens to be in-frame with rfstart).
        my $is_trunc5p_for_this_seq = 0;
        if($total_rf_offset_5p > 0) {
          $is_trunc5p_for_this_seq = 1;
        }
        elsif(scalar(@{$sgm_start_AA[$ftr_idx]}) > 0) {
          my $sgm0_rfstart = $sgm_start_AA[$ftr_idx][0];
          my $sgm0_rfstop  = $sgm_stop_AA[$ftr_idx][0];
          my $sgm0_strand  = $sgm_strand_AA[$ftr_idx][0];
          # Biological 5' RF pos: max for minus strand, min for plus strand
          my $rf_5p = ($sgm0_strand eq "-") ?
            ($sgm0_rfstart > $sgm0_rfstop ? $sgm0_rfstart : $sgm0_rfstop) :
            ($sgm0_rfstart < $sgm0_rfstop ? $sgm0_rfstart : $sgm0_rfstop);
          my $apos_5p = $msa->rfpos_to_aligned_pos($rf_5p);
          if($apos_5p >= 1 && $apos_5p <= length($aligned_sqstring)) {
            my $char_5p = substr($aligned_sqstring, $apos_5p - 1, 1);
            if($char_5p =~ /[\-\_\.\~]/) { $is_trunc5p_for_this_seq = 1; }
          }
        }
        
        # For 3' truncation: check if this sequence has ungapped residues at the biological 3' end.
        # For positive strand: biological 3' end is at high RF positions (stop codon at end)
        # For negative strand: biological 3' end is at low RF positions (stop codon at beginning)
        my $is_trunc3p_for_this_seq = 0;
        my $last_sgm_idx = scalar(@{$sgm_start_AA[$ftr_idx]}) - 1;
        if($last_sgm_idx >= 0) {
          my $strand = $sgm_strand_AA[$ftr_idx][$last_sgm_idx];
          my $ref_start_rfpos = $sgm_start_AA[$ftr_idx][$last_sgm_idx];
          my $ref_stop_rfpos = $sgm_stop_AA[$ftr_idx][$last_sgm_idx];
          
          # For negative strand, start > stop (e.g., 81..1:-), so we need to swap for range checking
          my ($min_rfpos, $max_rfpos) = ($ref_start_rfpos < $ref_stop_rfpos) ? 
                                         ($ref_start_rfpos, $ref_stop_rfpos) : 
                                         ($ref_stop_rfpos, $ref_start_rfpos);
          
          # Determine which 3 RF positions should contain the stop codon
          my @stop_codon_rfpos_A = ();
          if($strand eq "+") {
            # Positive strand: biological 3' is at high RF positions (ref_stop_rfpos and two before it)
            push(@stop_codon_rfpos_A, $ref_stop_rfpos, $ref_stop_rfpos-1, $ref_stop_rfpos-2);
          } else {
            # Negative strand: biological 3' is at low RF positions (ref_stop_rfpos and two after it)
            # Note: for negative strand coords like 81..1:-, ref_stop_rfpos=1 (the biological 3' end)
            push(@stop_codon_rfpos_A, $ref_stop_rfpos, $ref_stop_rfpos+1, $ref_stop_rfpos+2);
          }
          
          # Check if this sequence has ungapped residues at all 3 stop codon RF positions
          my $stop_codon_rf_ungapped_count = 0;
          foreach my $check_rfpos (@stop_codon_rfpos_A) {
            if($check_rfpos >= $min_rfpos && $check_rfpos <= $max_rfpos) {
              my $apos = $msa->rfpos_to_aligned_pos($check_rfpos);
              my $c = substr($aligned_sqstring, $apos-1, 1);
              my $is_gap = ($c =~ /[\-\_\.\~]/) ? 1 : 0;
              if(! $is_gap) { $stop_codon_rf_ungapped_count++; }
            }
          }
          # If any of the stop codon RF positions are gapped, this sequence is 3'-truncated
          if($stop_codon_rf_ungapped_count < 3) {
            $is_trunc3p_for_this_seq = 1;
          }
        }
        
        # Compute codon_start: the position in the extracted sequence where
        # the in-frame translation starts (1-based)
        # With $total_rf_offset_5p out-of-frame nucleotides at the 5' end,
        # the translation frame starts at position ($total_rf_offset_5p + 1)
        my $codon_start = $total_rf_offset_5p + 1;
        
        # Format CDS header with truncation markers and codon_start
        # Rebuild coordinate string with < and > markers if truncated
        if($ref_is_trunc5p || $is_trunc5p_for_this_seq || $is_trunc3p_for_this_seq) {
          my @marked_coords_A = ();
          for(my $i = 0; $i < scalar(@final_seq_sgm_coords_A); $i++) {
            my $coord_str = $final_seq_sgm_coords_A[$i];
            if($coord_str =~ /^(\d+)\.\.(\d+):([+-])$/) {
              my ($seg_start, $seg_stop, $seg_strand) = ($1, $2, $3);
              
              # For positive strand: coords are low..high:+, so start=low=5', stop=high=3'
              # For negative strand: coords are high..low:-, so start=high=5', stop=low=3'
              
              if($seg_strand eq "+") {
                # Positive strand: 5' marker at start, 3' marker at stop
                if($i == 0 && ($ref_is_trunc5p || $is_trunc5p_for_this_seq)) {
                  $seg_start = "<" . $seg_start;
                }
                if($i == scalar(@final_seq_sgm_coords_A) - 1 && $is_trunc3p_for_this_seq) {
                  $seg_stop = ">" . $seg_stop;
                }
              } else {
                # Negative strand: coords are high..low:-, 5' marker at start (high), 3' marker at stop (low)
                if($i == 0 && ($ref_is_trunc5p || $is_trunc5p_for_this_seq)) {
                  $seg_start = "<" . $seg_start;
                }
                if($i == scalar(@final_seq_sgm_coords_A) - 1 && $is_trunc3p_for_this_seq) {
                  $seg_stop = ">" . $seg_stop;
                }
              }
              push(@marked_coords_A, $seg_start . ".." . $seg_stop . ":" . $seg_strand);
            }
            else {
              push(@marked_coords_A, $coord_str); # fallback if format doesn't match
            }
          }
          $seq_coords_str = join(",", @marked_coords_A);
        }
        
        my $cds_header = $sqname . "/" . $seq_coords_str;
        if($codon_start != 1) {
          $cds_header .= "/CS" . $codon_start;
        }
        # Add refcoords as a comment so we can reconstruct the final protein header later
        $cds_header .= " REFCOORDS=" . $ref_coords_str;

        # In profile mode, only include sequences whose CDS is complete
        # (ends with a stop codon and has no premature in-frame stops).
        # Sequences with truncated CDS (no stop codon — typically CDS-only
        # GenBank submissions) or frameshifts are skipped; only complete
        # proteins are needed in the protein BLAST db.
        if(! profile_ValidateCdsIsComplete($final_cds)) {
          ofile_OutputString($FH_HR->{"log"}, 1,
            sprintf("# WARNING: profile_CdsFetchStockholmToFasta: CDS for %s at %s is incomplete or invalid, skipping (only complete CDS are included in protein db)\n",
                    $sqname, $seq_coords_str));
          next;
        }

        print $out_FH ">" . $cds_header . "\n";
        print $out_FH seq_SqstringAddNewlines($final_cds, 60);
      }
    }
  }

  undef $msa;
  return;
}

#################################################################
# Subroutine: profile_ComputeRfOffsets()
#
# Purpose: For a given CDS segment, compute how many nucleotides at the
#          5' and 3' ends need to be trimmed to maintain the correct reading frame.
#          When a sequence has gapped RF positions at the 5' end of a CDS,
#          the extracted sequence starts out-of-frame relative to the reference.
#          We need to trim nucleotides to reach the first position that's in
#          the same frame as the reference CDS start.
#
# Arguments:
#   $aligned_sqstring: aligned sequence string for this sequence
#   $astart:           aligned start position of the segment (1-based)
#   $astop:            aligned stop position of the segment (1-based)
#   $rfstart:          RF start position of the segment (1-based)
#   $rfstop:           RF stop position of the segment (1-based)
#   $msa:              Bio::Easel::MSA object
#   $strand:           strand of the segment ("+" or "-")
#
# Returns: ($offset_5p, $offset_3p)
#          $offset_5p: number of nucleotides to trim from 5' end to get in-frame
#          $offset_3p: number of nucleotides to trim from 3' end
#################################################################

#################################################################
# Subroutine: profile_ExtractOneCds()
# Incept:     EPN, Fri Apr 04 2026
#
# Purpose:    Extract a CDS nucleotide sequence for one sequence
#             using the coords of a specified feature index.
#             This is the core extraction logic from
#             profile_CdsFetchStockholmToFasta(), refactored to
#             allow trying multiple alternative CDS coords per seq.
#
# Returns:    ($cds_string, $header_string) on success
#             (undef, undef) if the CDS region is all gaps
#################################################################
sub profile_ExtractOneCds {
  my $sub_name = "profile_ExtractOneCds";
  my ($msa, $seq_idx, $sqname, $aligned_sqstring, $ftr_idx, $ftr_info_AHR,
      $sgm_start_AAR, $sgm_stop_AAR, $sgm_strand_AAR, $FH_HR) = @_;

  my $cds_sqstring = "";
  my $total_rf_offset_5p = 0;
  my $total_rf_offset_3p = 0;

  foreach(my $sgm_idx = 0; $sgm_idx < scalar(@{$sgm_start_AAR->[$ftr_idx]}); $sgm_idx++) {
    my $rfstart = $sgm_start_AAR->[$ftr_idx][$sgm_idx];
    my $rfstop  = $sgm_stop_AAR->[$ftr_idx][$sgm_idx];
    my $astart  = $msa->rfpos_to_aligned_pos($rfstart);
    my $astop   = $msa->rfpos_to_aligned_pos($rfstop);
    if($astart > $astop) { utl_Swap(\$astart, \$astop); }

    my ($ua_first, $ua_last) = profile_FirstAndLastUngappedPositionsInAlignedRange($aligned_sqstring, $astart, $astop);
    if(! defined $ua_first) {
      return (undef, undef);  # all gaps in this segment
    }

    my $sgm_sqstring = $msa->get_sqstring_unaligned_and_truncated($seq_idx, $astart, $astop);
    if($sgm_strand_AAR->[$ftr_idx][$sgm_idx] eq "-") {
      seq_SqstringReverseComplement(\$sgm_sqstring);
    }

    my ($rf_offset_5p, $rf_offset_3p) = profile_ComputeRfOffsets($aligned_sqstring, $astart, $astop, $rfstart, $rfstop, $msa, $sgm_strand_AAR->[$ftr_idx][$sgm_idx]);
    $total_rf_offset_5p += $rf_offset_5p;
    $total_rf_offset_3p += $rf_offset_3p;

    $cds_sqstring .= $sgm_sqstring;
  }

  my $cds_len = length($cds_sqstring);
  if($cds_len == 0 || ($total_rf_offset_5p + $total_rf_offset_3p >= $cds_len)) {
    return (undef, undef);
  }

  # Compute sequence coordinates
  my @final_seq_sgm_coords_A = ();
  foreach(my $sgm_idx = 0; $sgm_idx < scalar(@{$sgm_start_AAR->[$ftr_idx]}); $sgm_idx++) {
    my $rfstart = $sgm_start_AAR->[$ftr_idx][$sgm_idx];
    my $rfstop  = $sgm_stop_AAR->[$ftr_idx][$sgm_idx];
    my $astart  = $msa->rfpos_to_aligned_pos($rfstart);
    my $astop   = $msa->rfpos_to_aligned_pos($rfstop);
    if($astart > $astop) { utl_Swap(\$astart, \$astop); }
    my ($ua_first, $ua_last) = profile_FirstAndLastUngappedPositionsInAlignedRange($aligned_sqstring, $astart, $astop);
    if(defined $ua_first && $ua_first <= $ua_last) {
      my $seg_strand = $sgm_strand_AAR->[$ftr_idx][$sgm_idx];
      if($seg_strand eq "-") {
        push(@final_seq_sgm_coords_A, $ua_last . ".." . $ua_first . ":" . $seg_strand);
      } else {
        push(@final_seq_sgm_coords_A, $ua_first . ".." . $ua_last . ":" . $seg_strand);
      }
    }
  }

  my $seq_coords_str = join(",", @final_seq_sgm_coords_A);
  my $ref_coords_str = $ftr_info_AHR->[$ftr_idx]{"coords"};

  # Build header with codon_start if needed
  my $codon_start = $total_rf_offset_5p + 1;
  my $cds_header = $sqname . "/" . $seq_coords_str;
  if($codon_start != 1) {
    $cds_header .= "/CS" . $codon_start;
  }
  $cds_header .= " REFCOORDS=" . $ref_coords_str;

  return ($cds_sqstring, $cds_header);
}

#################################################################
# Subroutine: profile_ValidateCdsIsComplete()
# Incept:     EPN* Sun Apr 13 2026
#
# Purpose:    Check if a CDS nucleotide sequence extracted from a
#             profile alignment is complete and suitable for
#             inclusion in the protein BLAST db:
#             - Must end with a stop codon (TAA/TAG/TGA)
#             - No in-frame stop codons before the terminal one
#             - Length minus 3 (stop) must be divisible by 3
#
#             Sequences that fail (truncated CDS, frameshifts,
#             premature stops) are excluded from the protein db
#             but remain in the nucleotide alignment.
#
# Arguments:
#   $cds_seq: CDS nucleotide sequence string
#
# Returns:    1 if complete and valid, 0 if not
#################################################################
sub profile_ValidateCdsIsComplete {
  my ($cds_seq) = @_;

  my $len = length($cds_seq);
  if($len < 6) { return 0; }

  my %stop_codons = ("TAA" => 1, "TAG" => 1, "TGA" => 1,
                     "taa" => 1, "tag" => 1, "tga" => 1);

  # Last codon must be a stop
  my $last_codon = substr($cds_seq, $len - 3, 3);
  if(! exists $stop_codons{$last_codon}) {
    return 0;  # no terminal stop — truncated/incomplete CDS
  }

  # Length minus stop must be divisible by 3
  if(($len - 3) % 3 != 0) {
    return 0;
  }

  # No premature in-frame stops before the terminal one
  for(my $i = 0; $i < $len - 3; $i += 3) {
    my $codon = substr($cds_seq, $i, 3);
    if(exists $stop_codons{$codon}) {
      return 0;  # premature stop
    }
  }

  return 1;
}

sub profile_ComputeRfOffsets {
  my $sub_name = "profile_ComputeRfOffsets";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($aligned_sqstring, $astart, $astop, $rfstart, $rfstop, $msa, $strand) = (@_);

  # For negative strand, rfstart > rfstop (e.g., 81..1:-), so we need to scan in reverse
  my ($min_rf, $max_rf) = ($rfstart <= $rfstop) ? ($rfstart, $rfstop) : ($rfstop, $rfstart);

  # Scan the RF span to find the first and last RF positions where this sequence is ungapped
  my $first_ungapped_rf = undef;
  my $last_ungapped_rf  = undef;

  for(my $rfpos = $min_rf; $rfpos <= $max_rf; $rfpos++) {
    my $apos = $msa->rfpos_to_aligned_pos($rfpos);
    my $c = substr($aligned_sqstring, $apos-1, 1);
    my $is_gap = ($c =~ /[\-\_\.\~]/) ? 1 : 0;
    if(! $is_gap) {
      if(! defined $first_ungapped_rf) { $first_ungapped_rf = $rfpos; }
      $last_ungapped_rf = $rfpos;
    }
  }

  if(! defined $first_ungapped_rf) {
    # Entire segment is gapped in this sequence
    return (0, 0);
  }

  # For 5' offset: we need to trim nucleotides until we reach a position
  # that's in the same reading frame as the reference CDS start (rfstart).
  # The reference frame is (rfstart % 3).
  # We want the first RF position where (rfpos % 3) == (rfstart % 3).
  # For negative strand, the biological 5' end is at the HIGH RF position (rfstart),
  # but in our scan we go low to high, so we need to find the first in-frame position
  # scanning FROM the biological 5' end.
  my $ref_frame = $rfstart % 3;
  my $offset_5p = 0;
  
  # Determine scan direction based on strand
  if($strand eq "+") {
    # Positive strand: scan forward from first ungapped
    for(my $rfpos = $first_ungapped_rf; $rfpos <= $max_rf; $rfpos++) {
      my $apos = $msa->rfpos_to_aligned_pos($rfpos);
      my $c = substr($aligned_sqstring, $apos-1, 1);
      my $is_gap = ($c =~ /[\-\_\.\~]/) ? 1 : 0;
      
      if(! $is_gap) {
        my $this_frame = $rfpos % 3;
        if($this_frame == $ref_frame) {
          # Found the first in-frame position, stop counting
          last;
        }
        $offset_5p++;
      }
    }
  } else {
    # Negative strand: scan backward from last ungapped (biological 5' end)
    for(my $rfpos = $last_ungapped_rf; $rfpos >= $min_rf; $rfpos--) {
      my $apos = $msa->rfpos_to_aligned_pos($rfpos);
      my $c = substr($aligned_sqstring, $apos-1, 1);
      my $is_gap = ($c =~ /[\-\_\.\~]/) ? 1 : 0;
      
      if(! $is_gap) {
        my $this_frame = $rfpos % 3;
        if($this_frame == $ref_frame) {
          # Found the first in-frame position (from biological 5' end), stop counting
          last;
        }
        $offset_5p++;
      }
    }
  }

  # For 3' offset: count ungapped nucleotides after the biological 3' end
  # For positive strand: biological 3' is after the last ungapped RF position
  # For negative strand: biological 3' is before the first ungapped RF position
  my $offset_3p = 0;
  if($strand eq "+") {
    for(my $rfpos = $last_ungapped_rf + 1; $rfpos <= $max_rf; $rfpos++) {
      my $apos = $msa->rfpos_to_aligned_pos($rfpos);
      my $c = substr($aligned_sqstring, $apos-1, 1);
      my $is_gap = ($c =~ /[\-\_\.\~]/) ? 1 : 0;
      if(! $is_gap) { $offset_3p++; }
    }
  } else {
    for(my $rfpos = $first_ungapped_rf - 1; $rfpos >= $min_rf; $rfpos--) {
      my $apos = $msa->rfpos_to_aligned_pos($rfpos);
      my $c = substr($aligned_sqstring, $apos-1, 1);
      my $is_gap = ($c =~ /[\-\_\.\~]/) ? 1 : 0;
      if(! $is_gap) { $offset_3p++; }
    }
  }

  return ($offset_5p, $offset_3p);
}

#################################################################
# Subroutine: profile_RenameProteinHeaders()
#
# Purpose: Rename protein FASTA headers from temporary format (seqname/seqcoords)
#          to final profile format (seqname:seqcoords/refcoords).
#          Reads REFCOORDS from CDS FASTA comments and removes truncation markers.
#
# Arguments:
#   $tmp_protein_file: temporary protein FASTA file (input)
#   $out_FH:           output file handle for final protein FASTA
#   $cds_fa_file:      CDS FASTA file (contains REFCOORDS comments)
#   $FH_HR:            REF to hash of file handles
#
# Returns: void
#################################################################
sub profile_RenameProteinHeaders {
  my $sub_name = "profile_RenameProteinHeaders";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($tmp_protein_file, $out_FH, $cds_fa_file, $FH_HR) = (@_);

  # Build mapping from CDS header (seqname/seqcoords) to refcoords
  my %refcoords_H = ();
  open(CDS, $cds_fa_file) || ofile_FileOpenFailure($cds_fa_file, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <CDS>) {
    if($line =~ /^>(\S+)\s+REFCOORDS=(\S+)/) {
      my ($cds_name, $refcoords) = ($1, $2);
      $refcoords_H{$cds_name} = $refcoords;
    }
  }
  close(CDS);

  # Build reverse lookup: protein_name (without markers or CS#) -> CDS key (with REFCOORDS)
  my %protein_to_cds_H = ();
  foreach my $cds_key (keys %refcoords_H) {
    my $protein_key = $cds_key;
    # Remove truncation markers < and >
    $protein_key =~ s/\<//g;
    $protein_key =~ s/\>//g;
    # Remove /CS# annotation (codon_start)
    $protein_key =~ s/\/CS\d$//;
    $protein_to_cds_H{$protein_key} = $cds_key;
  }

  # Read temp protein file and rename headers (preserving order)
  open(PROT, $tmp_protein_file) || ofile_FileOpenFailure($tmp_protein_file, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <PROT>) {
    if($line =~ /^>(\S+)/) {
      my $protein_name = $1;
      # protein_name is from sqf_EslTranslateCdsToFastaFile output: seqname/seqcoords (without <> markers)
      
      if(exists $protein_to_cds_H{$protein_name}) {
        my $cds_key = $protein_to_cds_H{$protein_name};
        my $refcoords = $refcoords_H{$cds_key};
        # Convert seqname/seqcoords to seqname:seqcoords/refcoords
        if($protein_name =~ /^([^\/]+)\/(.+)$/) {
          my ($seqname, $seqcoords) = ($1, $2);
          print $out_FH ">" . $seqname . ":" . $seqcoords . "/" . $refcoords . "\n";
        }
        else {
          ofile_FAIL("ERROR in $sub_name, unable to parse protein header: $protein_name\n", 1, $FH_HR);
        }
      }
      else {
        ofile_FAIL("ERROR in $sub_name, unable to find REFCOORDS for protein: $protein_name\n", 1, $FH_HR);
      }
    }
    else {
      print $out_FH $line;
    }
  }
  close(PROT);

  return;
}

#################################################################
# Subroutine: profile_EslTranslateCdsToFastaFile()
#
# Purpose: Translate CDS sequences to proteins for profile mode.
#          Adapted from sqf_EslTranslateCdsToFastaFile in sequip/sqp_seqfile.pm
#          to handle profile-mode header format: seqname:seqcoords/refcoords/CS#
#          (vs. single-seq format: seqname/coords)
#
# Arguments:
#   $out_FH:         output file handle for protein FASTA
#   $esl_translate:  path to esl-translate executable
#   $cds_fa_file:    input CDS FASTA file
#   $out_root:       output root for temp files
#   $ftr_info_AHR:   REF to feature info array
#   $opt_HHR:        REF to 2D hash of options
#   $FH_HR:          REF to hash of file handles
#
# Returns: void
#################################################################
sub profile_EslTranslateCdsToFastaFile {
  my $sub_name = "profile_EslTranslateCdsToFastaFile";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($out_FH, $esl_translate, $cds_fa_file, $out_root, $ftr_info_AHR, $opt_HHR, $FH_HR) = @_;

  my $tmp1_translate_fa_file  = $out_root . ".cds.esl-translate.1.fa";
  my $tmp2_translate_fa_file  = $out_root . ".cds.esl-translate.2.fa";
  my $tmp1_translate_ssi_file = $out_root . ".cds.esl-translate.1.fa.ssi";
  my $tmp2_translate_ssi_file = $out_root . ".cds.esl-translate.2.fa.ssi";
  if(-e $tmp1_translate_ssi_file) { unlink $tmp1_translate_ssi_file; }
  if(-e $tmp2_translate_ssi_file) { unlink $tmp2_translate_ssi_file; }

  my $c_opt = "";
  if((opt_IsUsed("--ttbl", $opt_HHR)) && (opt_Get("--ttbl", $opt_HHR) != 1)) {
    $c_opt = "-c " . opt_Get("--ttbl", $opt_HHR);
  }
  my $translate_cmd = "$esl_translate $c_opt -l 3 --watson $cds_fa_file > $tmp1_translate_fa_file";
  utl_RunCommand($translate_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # Rewrite names so we can fetch by source+coords
  open(IN,       $tmp1_translate_fa_file) || ofile_FileOpenFailure($tmp1_translate_fa_file, $sub_name, $!, "reading", $FH_HR);
  open(OUT, ">", $tmp2_translate_fa_file) || ofile_FileOpenFailure($tmp2_translate_fa_file, $sub_name, $!, "writing", $FH_HR);
  while(my $line = <IN>) {
    if($line =~ m/^\>/) {
      chomp $line;
      if($line =~ /^\>orf\d+\s+(source\=\S+)\s+(coords\=\S+)\s+length\=\d+\s+frame\=\S+/) {
        print OUT (">" . $1 . "," . $2 . "\n");
      }
      else {
        ofile_FAIL("ERROR in $sub_name, problem parsing esl-translate output file $tmp1_translate_fa_file, line:\n$line\n", 1, $FH_HR);
      }
    }
    else {
      print OUT $line;
    }
  }
  close(IN);
  close(OUT);

  # Fetch expected translated sequences
  my $cds_sqfile     = Bio::Easel::SqFile->new({ fileLocation => $cds_fa_file });
  my $protein_sqfile = Bio::Easel::SqFile->new({ fileLocation => $tmp2_translate_fa_file });

  for(my $seq_idx = 0; $seq_idx < $cds_sqfile->nseq_ssi; $seq_idx++) {
    my ($seq_name, $seq_length) = $cds_sqfile->fetch_seq_name_and_length_given_ssi_number($seq_idx);
    
    # Parse profile-mode header: seqname:seqcoords/refcoords/CS#
    # Split by : first to separate seqname from the rest
    my @colon_parts = split(":", $seq_name, 2);
    if(scalar(@colon_parts) < 2) {
      ofile_FAIL("ERROR in $sub_name, unable to parse profile CDS header (expected seqname:seqcoords/refcoords): $seq_name\n", 1, $FH_HR);
    }
    my $sqname = $colon_parts[0];
    my $coord_part = $colon_parts[1];
    
    # Split coord_part by / to get seqcoords, refcoords, and optional CS#
    my @slash_parts = split("/", $coord_part);
    if(scalar(@slash_parts) < 2) {
      ofile_FAIL("ERROR in $sub_name, unable to parse profile CDS coords (expected seqcoords/refcoords): $coord_part\n", 1, $FH_HR);
    }
    my $seq_coords_str = $slash_parts[0];
    my $ref_coords_str = $slash_parts[1];
    my $codon_start = 1;
    if(scalar(@slash_parts) >= 3 && $slash_parts[2] =~ /^CS(\d)$/) {
      $codon_start = $1;
    }
    
    # Determine truncation status from seq_coords_str
    my $is_trunc5 = ($seq_coords_str =~ /</) ? 1 : 0;
    my $is_trunc3 = ($seq_coords_str =~ />/) ? 1 : 0;
    
    # Compute expected translation coordinates
    my $expected_start = 1;
    my $expected_stop  = $seq_length;
    if(! $is_trunc3) {
      $expected_stop -= 3; # stop codon won't be translated
    }
    else {
      if($codon_start == 1)    { $expected_stop -= ($seq_length % 3); }
      elsif($codon_start == 2) { $expected_stop -= (($seq_length-1) % 3); }
      elsif($codon_start == 3) { $expected_stop -= (($seq_length-2) % 3); }
    }
    if($codon_start == 2) { $expected_start = 2; }
    if($codon_start == 3) { $expected_start = 3; }
    
    # Fetch the translation
    my $fetch_name = "source=" . $seq_name . ",coords=" . $expected_start . ".." . $expected_stop;
    if(! $protein_sqfile->check_seq_exists($fetch_name)) {
      ofile_FAIL("ERROR in $sub_name, problem translating CDS feature, unable to find expected translated sequence in $tmp2_translate_fa_file:\n\tseq: $seq_name\n\texpected sequence:$fetch_name\n", 1, $FH_HR);
    }
    
    # Build output protein header: seqname:seqcoords/refcoords
    my $protein_sqname = $sqname . ":" . $seq_coords_str . "/" . $ref_coords_str;
    # Remove truncation markers for protein header
    $protein_sqname =~ s/\<//g;
    $protein_sqname =~ s/\>//g;
    
    print $out_FH ">" . $protein_sqname . "\n";
    my $protein_sqstring = $protein_sqfile->fetch_seq_to_sqstring($fetch_name);
    if(! $is_trunc5) {
      if($protein_sqstring !~ m/^M/) {
        ofile_FAIL("ERROR in $sub_name, problem translating CDS feature, feature does not seem to be 5' truncated but translated protein does not start with an M:\n\tseq: $seq_name\n\texpected sequence:$fetch_name\n", 1, $FH_HR);
      }
    }
    print $out_FH seq_SqstringAddNewlines($protein_sqstring, 60);
  }

  # Cleanup temp files unless --keep
  if(! opt_Get("--keep", $opt_HHR)) {
    utl_FileRemoveUsingSystemRm($tmp1_translate_fa_file, $sub_name, $opt_HHR, $FH_HR);
    utl_FileRemoveUsingSystemRm($tmp2_translate_fa_file, $sub_name, $opt_HHR, $FH_HR);
    utl_FileRemoveUsingSystemRm($tmp2_translate_fa_file . ".ssi", $sub_name, $opt_HHR, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: profile_FirstAndLastUngappedPositionsInAlignedRange()
#
# Purpose: Given an aligned sequence string and an aligned coordinate
#          range (1..alen), return the first and last dealigned
#          (ungapped) positions spanned by any non-gap residues within
#          that aligned range.
#
# Returns: (undef, undef) if the range contains no residues.
#################################################################
sub profile_FirstAndLastUngappedPositionsInAlignedRange {
  my $sub_name = "profile_FirstAndLastUngappedPositionsInAlignedRange";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

  my ($aligned_sqstring, $astart, $astop) = (@_);

  my $alen = length($aligned_sqstring);
  if($astart < 1 || $astop < 1 || $astart > $alen || $astop > $alen) {
    return (undef, undef);
  }
  if($astart > $astop) { my $tmp = $astart; $astart = $astop; $astop = $tmp; }

  my $ua_pos = 0;
  my $ua_first = undef;
  my $ua_last  = undef;

  for(my $apos = 1; $apos <= $alen; $apos++) {
    my $c = substr($aligned_sqstring, $apos-1, 1);
    my $is_gap = ($c =~ /[\-\_\.\~]/) ? 1 : 0;
    if(! $is_gap) { $ua_pos++; }

    if(($apos >= $astart) && ($apos <= $astop) && (! $is_gap)) {
      if(! defined $ua_first) { $ua_first = $ua_pos; }
      $ua_last = $ua_pos;
    }
  }

  return ($ua_first, $ua_last);
}

#################################################################
# Subroutine: process_add_and_skip_options()
# Incept:     EPN, Mon Mar 18 06:29:21 2019
#
# Synopsis: Process cmdline --{f,q}add and --{f,q}skip options 
#           for features or qualifiers.
#
# Arguments:
#  $df_string:  comma separated string of default values (e.g. "CDS,gene,mat_peptide,sig_peptide")
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
            #printf("adding feature " . $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"type"} . " with coords " . $prot_ftr_info_HAH{$prot_accver}[$prot_ftr_idx]{"coords"} . "\n");
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
      my $parent_idx_str = $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"};
      $parent_idx_str =~ s/\s+//g;

      # If parent_idx_str already looks like a comma-separated list of integers,
      # nothing to do (this is the format expected by v-annotate.pl and is what
      # we read from .minfo files in --profile mode).
      if($parent_idx_str =~ /^\d+(?:,\d+)*$/) {
        $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} = $parent_idx_str;
        next;
      }

      my @parent_type_coords_A = split("!GBSEP!", $parent_idx_str);
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

#################################################################
# Subroutine: check_and_add_cds_splice_sites
# Incept:     EPN, Mon Jul 24 13:04:00 2023
# 
# Purpose:    Check any CDS features to see if they have any
#             introns, and if so if those introns have canonical
#             GT donor (5') and AG acceptor (3') sites. If all 
#             introns have canonical splice sites, then add 
#             canon_splice_sites="1" qualifier for that CDS.
# 
# Arguments:
#   $stk_file:      stockholm file with aligned full length sequences
#   $ftr_info_AHR:  REF to feature information, added to here
#   $opt_HHR:       REF to 2D hash of option values, see top of sqp_opts.pm for description, PRE-FILLED
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub check_and_add_cds_splice_sites { 
  my $sub_name = "check_and_add_cds_splice_sites";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($stk_file, $ftr_info_AHR, $sgm_info_AHR, $opt_HHR, $FH_HR) = @_;

  my $msa = Bio::Easel::MSA->new({ fileLocation => $stk_file, isDna => 1});
  my $msa_has_rf = $msa->has_rf;
  my $nseq = $msa->nseq;

  my $min_intron_length = opt_Get("--intlen", $opt_HHR);

  # precompute start, stop, strand, for all features, so we don't have to redo this for each seq
  my @sgm_start_AA  = ();
  my @sgm_stop_AA   = ();
  my @sgm_strand_AA = ();
  vdr_FeatureInfoStartStopStrandArrays($ftr_info_AHR, \@sgm_start_AA, \@sgm_stop_AA, \@sgm_strand_AA, $FH_HR);

  my $nftr = scalar(@{$ftr_info_AHR});

  my $canon_5p; # TRUE if all seqs have canonical 5' splice site (GT)
  my $canon_3p; # TRUE if all seqs have canonical 3' splice site (AG)
  my ($rfstart, $rfstop, $astart, $astop); # model/alignment positions 
  my ($nsgm, $next_sgm_idx, $strand);
  my $seq_idx;     # sequence index in the MSA (always 0 if 1 seq MSA)
  my $ss_sqstring; # the splice site string
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") { 
      $nsgm = scalar(@{$sgm_start_AA[$ftr_idx]});
      if($nsgm > 1) {
        # this function is authoritative for canon_splice_sites on multi-segment CDSs:
        # clear any stale qualifier inherited via --minfoin so we only re-set it below
        # when the splice sites actually validate
        delete $ftr_info_AHR->[$ftr_idx]{"canon_splice_sites"};
        $canon_5p = 1; # will set to 0 below if any 5' splice site for any intron is not GT
        $canon_3p = 1; # will set to 0 below if any 3' splice site for any intron is not AG
        my $nintron_checked = 0; # number of introns whose splice sites we actually inspected;
                                 # if 0 (e.g., joined-coord ribosomal-slippage CDSs whose
                                 # intersegment gap is below --intlen) do not set
                                 # canon_splice_sites:"1" vacuously

        # NOTE: as currently implemented, if there are more than one
        # introns, splice sites have to be canonical in all for
        # canon_splice_sites:"1" to get set, else we won't set it

        foreach(my $sgm_idx = 0; $sgm_idx < ($nsgm-1); $sgm_idx++) {
          $next_sgm_idx = $sgm_idx+1;
          # determine if the intron is >= 40 nt
          $strand = $sgm_strand_AA[$ftr_idx][$sgm_idx];
          if($strand ne $sgm_strand_AA[$ftr_idx][$next_sgm_idx]) { 
            ofile_FAIL(sprintf("ERROR in $sub_name, with --fss, all adjacent CDS segments must be on the same strand, but CDS with coordinates %s violates this", $ftr_info_AHR->[$ftr_idx]{"coords"}), 1, $FH_HR);
          }
          my $intron_length = vdr_FeatureLengthBetweenAdjacentSegments($ftr_info_AHR, $sgm_info_AHR, $ftr_idx, $sgm_idx, $FH_HR);

          if($intron_length >= $min_intron_length) {
            $nintron_checked++;
            # check 5' splice site
            $rfstart = ($strand eq "+") ? $sgm_stop_AA[$ftr_idx][$sgm_idx] + 1 : $sgm_stop_AA[$ftr_idx][$sgm_idx] - 1;
            $rfstop  = ($strand eq "+") ? $sgm_stop_AA[$ftr_idx][$sgm_idx] + 2 : $sgm_stop_AA[$ftr_idx][$sgm_idx] - 2;
            $astart  = ($msa_has_rf) ? $msa->rfpos_to_aligned_pos($rfstart) : $rfstart;
            $astop   = ($msa_has_rf) ? $msa->rfpos_to_aligned_pos($rfstop)  : $rfstop;
            if($astart > $astop) { utl_Swap(\$astart, \$astop); }
            for($seq_idx = 0; $seq_idx < $msa->nseq; $seq_idx++) { 
              $ss_sqstring = $msa->get_sqstring_unaligned_and_truncated($seq_idx, $astart, $astop);
              if($sgm_strand_AA[$ftr_idx][$sgm_idx] eq "-") { 
                seq_SqstringReverseComplement(\$ss_sqstring);
              }
              $ss_sqstring =~ tr/a-z/A-Z/; # convert to uppercase
              $ss_sqstring =~ tr/U/T/;     # convert to DNA
              if(length($ss_sqstring) != 2) {
                # exit if --strictss
                if(opt_Get("--strictss", $opt_HHR)) { 
                  ofile_FAIL(sprintf("ERROR in $sub_name, with --strictss, splice sites expected to be length 2, but got length of %d for 5' splice site for seq#%d (positions $astart..$astop:$sgm_strand_AA[$ftr_idx][$sgm_idx])", 
                                     length($ss_sqstring), $seq_idx), 1, $FH_HR);
                }
                $canon_5p = 0;
              }
              if($ss_sqstring ne "GT") { 
                if(opt_Get("--strictss", $opt_HHR)) { 
                  ofile_FAIL(sprintf("ERROR in $sub_name, with --strictss, 5' splice sites expected to be GT, but got %s for 5' splice site for seq#%d (positions $astart..$astop:$sgm_strand_AA[$ftr_idx][$sgm_idx])", 
                                     $ss_sqstring, $seq_idx), 1, $FH_HR);
                }
                $canon_5p = 0;
              }
            }
            
            # check 3' splice site
            $rfstart = ($strand eq "+") ? $sgm_start_AA[$ftr_idx][$next_sgm_idx] - 2 : $sgm_start_AA[$ftr_idx][$next_sgm_idx] + 2;
            $rfstop  = ($strand eq "+") ? $sgm_start_AA[$ftr_idx][$next_sgm_idx] - 1 : $sgm_start_AA[$ftr_idx][$next_sgm_idx] + 1;
            $astart  = ($msa_has_rf) ? $msa->rfpos_to_aligned_pos($rfstart) : $rfstart;
            $astop   = ($msa_has_rf) ? $msa->rfpos_to_aligned_pos($rfstop)  : $rfstop;
            if($astart > $astop) { utl_Swap(\$astart, \$astop); }
            for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
              $ss_sqstring = $msa->get_sqstring_unaligned_and_truncated($seq_idx, $astart, $astop);
              if($sgm_strand_AA[$ftr_idx][$sgm_idx] eq "-") { 
                seq_SqstringReverseComplement(\$ss_sqstring);
              }
              $ss_sqstring =~ tr/a-z/A-Z/; # convert to uppercase
              $ss_sqstring =~ tr/U/T/;     # convert to DNA
              if(length($ss_sqstring) != 2) {
                # exit if --strictss
                if(opt_Get("--strictss", $opt_HHR)) { 
                  ofile_FAIL(sprintf("ERROR in $sub_name, with --strictss, splice sites expected to be length 2, but got length of %d for 3' splice site for seq#%d (positions $astart..$astop:$sgm_strand_AA[$ftr_idx][$sgm_idx])", 
                                     length($ss_sqstring), $seq_idx), 1, $FH_HR);
                }
                $canon_3p = 0;
              }
              if($ss_sqstring ne "AG") { 
                if(opt_Get("--strictss", $opt_HHR)) { 
                  ofile_FAIL(sprintf("ERROR in $sub_name, with --strictss, 3' splice sites expected to be AG, but got %s for 5' splice site for seq#%d (positions $astart..$astop:$sgm_strand_AA[$ftr_idx][$sgm_idx])", 
                                     $ss_sqstring, $seq_idx), 1, $FH_HR);
                }
                $canon_3p = 0;
              }
            }
          } # end of if($intron_length >= $min_intron_length)
        } # end of for($sgm_idx...
        if($nintron_checked > 0 && $canon_5p && $canon_3p) {
          # set canon_splice_sites="1"
          $ftr_info_AHR->[$ftr_idx]{"canon_splice_sites"} = 1;
        }
      } # end of if($nsgm > 1)
    }
  }
  return;
}

