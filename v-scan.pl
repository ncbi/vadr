#!/usr/bin/env perl
# EPN, Wed May  1 10:18:55 2019 [renamed to v-annotate.pl]
# EPN, Thu Feb 18 12:48:16 2016 [dnaorg_annotate.pl split off from dnaorg_annotate_genomes.pl]
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
#
use strict;
use warnings;
use Getopt::Long qw(:config no_auto_abbrev);
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::Random;
use Bio::Easel::SqFile;

require "vadr.pm"; 
require "vadr_seed.pm"; 
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
#   o validates and parses input config file
#
# - Runs v-annotate.pl with --cls_only for each model library
#   listed in config file
# - Parses --cls_only output to determine which sequences to
#   annotate with each model library
# - For each model library with at least one sequence to annotate,
#   o fetches sequences to new fasta file
#   o runs v-annotate.pl again using options specified in config file
# - Summarizes results and exits
#
#######################################################################################

# make sure required environment variables are set
my $env_vadr_scripts_dir = utl_DirEnvVarValid("VADRSCRIPTSDIR");
my $env_vadr_easel_dir   = utl_DirEnvVarValid("VADREASELDIR");
my $env_vadr_config_file = (exists $ENV{"VADRCONFIGFILE"}) ? $ENV{"VADRCONFIGFILE"} : undef;

my %execs_H = (); # hash with paths to all required executables
$execs_H{"v-annotate.pl"} = $env_vadr_scripts_dir  . "/v-annotate.pl";
$execs_H{"esl-reformat"}  = $env_vadr_easel_dir    . "/esl-reformat";

utl_ExecHValidate(\%execs_H, undef);

#########################################################
# Command line and option processing using sqp_opts.pm
#
# opt_HH: 2D hash:
#         1D key: option name (e.g. "-h")
#         2D key: string denoting type of information 
#                 (one of "type", "default", "group", "requires", "incompatible", "preamble", "help")
#         value:  string explaining 2D key:
#                 "type":         "boolean", "string", "integer" or "real"
#                 "default":      default value for option
#                 "group":        integer denoting group number this option belongs to
#                 "requires":     string of 0 or more other options this option requires to work, each separated by a ','
#                 "incompatible": string of 0 or more other options this option is incompatible with, each separated by a ','
#                 "preamble":     string describing option for preamble section (beginning of output from script)
#                 "help":         string describing option for help section (printed if -h used)
#                 "setby":        '1' if option set by user, else 'undef'
#                 "value":        value for option, can be undef if default is undef
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
#     option            type       default group   requires incompat    preamble-output                                   help-output    
opt_Add("-h",           "boolean", 0,          0,    undef, undef,      undef,                                            "display this help",                                  \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "basic options";
#     option            type       default group   requires incompat    preamble-output                                                                            help-output    
opt_Add("-f",           "boolean", 0,         $g,    undef, undef,      "force directory overwrite",                                                               "force; if output dir exists, overwrite it",   \%opt_HH, \@opt_order_A);
opt_Add("-c",           "string",  0,         $g,    undef, undef,      "use config file <s> instead of default in \$VADRCONFIGFILE",                              "use config file <s> instead of default in \$VADRCONFIGFILE", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,         $g,    undef, undef,      "be verbose",                                                                              "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("-m",           "boolean", 0,         $g,    undef, undef,      "allow matches to multiple model libraries",                                               "allow matches to multiple model libraries", \%opt_HH, \@opt_order_A);
opt_Add("--lone",       "boolean", 0,         $g,    undef, undef,      "exit if at least one sequence matches to multiple libraries",                             "exit if at least one sequence matches to multiple libraries", \%opt_HH, \@opt_order_A);
opt_Add("--first",      "boolean", 0,         $g,    undef,"--lone",    "if a seq matches > 1 model library use first one [df: use best scoring]",                 "if a seq matches > 1 model library use first one [df: use best scoring]", \%opt_HH, \@opt_order_A);
opt_Add("--origfa",     "boolean", 0,         $g,    undef,   undef,    "do not copy fasta file prior to analysis, use original",                 "do not copy fasta file prior to analysis, use original", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,         $g,    undef, undef,      "leaving intermediate files on disk",                                                      "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);
#     option            type       default group   requires incompat    preamble-output                                                                            help-output    
$opt_group_desc_H{++$g} = "options for specifying which model libraries to use:";
opt_Add("--only",        "string", 0,         $g,   undef,"--skip",     "only use the model library(ies) in comma separated string <s>",                           "only use the model library(ies) in comma separated string <s>",   \%opt_HH, \@opt_order_A);
opt_Add("--skip",        "string", 0,         $g,   undef,"--only",     "do not use the model library(ies) in comma separated string <s>",                         "do nout use the model library(ies) in comma separated string <s>", \%opt_HH, \@opt_order_A);
#     option            type       default group   requires incompat    preamble-output                                                                            help-output    
$opt_group_desc_H{++$g} = "options related to the sampling of sequences for determining model library to use:";
opt_Add("--all",        "boolean", 0,         $g,   undef, undef,       "do not sample, pick model library(ies) based on all sequences (auto turned on if -m used)", "do not sample, pick model library(ies) based on all sequences (auto turned on if -m used)", \%opt_HH, \@opt_order_A);
opt_Add("--s_nseq",     "integer", 3,         $g,   undef, "--all",     "set the number of sequences to sample to <n>",                                            "set the number of sequences to sample to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--s_beg",      "boolean", 0,         $g,   undef, "--all",     "sample sequences from the beginning of the file, not randomly",                           "sample sequences from the beginning of the file, not randomly", \%opt_HH, \@opt_order_A);
opt_Add("--s_seed",     "integer", 181,       $g,undef,"--all,--s_beg", "set the random number generator seed to <n>",                                             "set the random number generator seed to <n>", \%opt_HH, \@opt_order_A);
#     option            type       default group   requires incompat    preamble-output                                                                            help-output    
$opt_group_desc_H{++$g} = "options for listing information on models and exiting:";
opt_Add("--l_all",        "boolean", 0,         $g,    undef, undef,    "list all info about all model libraries in the config file and exit",                     "list all info about all model libraries in the config file and exit", \%opt_HH, \@opt_order_A);
opt_Add("--l_lib",        "string", 0,          $g,    undef, undef,    "list all info about model library with key <s> in the config file and exit",              "list all info about model library with key <s> in the config file and exit", \%opt_HH, \@opt_order_A);
opt_Add("--l_dir",        "boolean", 0,         $g,    undef, undef,    "list all model library directories in the config file and exit",                          "list all model library directories in the config file and exit", \%opt_HH, \@opt_order_A);
opt_Add("--l_opt",        "boolean", 0,         $g,    undef, undef,    "list all model library v-annotate.pl options in the config file and exit",                "list all model library v-annotate.pl options in the config file and exit", \%opt_HH, \@opt_order_A);
opt_Add("--l_mdl",        "boolean", 0,         $g,    undef, undef,    "list all models in the model libraries in the config file and exit",                      "list all models in the model libraries in the config file and exit", \%opt_HH, \@opt_order_A);


# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay = 
    &GetOptions('h'        => \$GetOptions_H{"-h"}, 
# basic options
                'f'        => \$GetOptions_H{"-f"},
                'c=s'      => \$GetOptions_H{"-c"},
                'l'        => \$GetOptions_H{"-l"},
                'v'        => \$GetOptions_H{"-v"},
                'm'        => \$GetOptions_H{"-m"}, 
                'lone'     => \$GetOptions_H{"--lone"}, 
                'first'    => \$GetOptions_H{"--first"},
                'origfa'   => \$GetOptions_H{"--origfa"},
                'keep'     => \$GetOptions_H{"--keep"}, 
                'only=s'   => \$GetOptions_H{"--only"}, 
                'skip=s'   => \$GetOptions_H{"--skip"}, 
                'all'      => \$GetOptions_H{"--all"}, 
                's_nseq=s' => \$GetOptions_H{"--s_nseq"},
                's_beg'    => \$GetOptions_H{"--s_beg"},
                's_seed=s' => \$GetOptions_H{"--s_seed"},
                'l_all'    => \$GetOptions_H{"--l_all"},
                'l_lib=s'  => \$GetOptions_H{"--l_lib"},
                'l_dir'    => \$GetOptions_H{"--l_dir"},
                'l_opt'    => \$GetOptions_H{"--l_opt"},
                'l_mdl'    => \$GetOptions_H{"--l_mdl"});
                

my $total_seconds = -1 * ofile_SecondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $execname_opt  = $GetOptions_H{"--execname"};
my $executable    = (defined $execname_opt) ? $execname_opt : "v-scan.pl";
my $usage         = "Usage: $executable [-options] <fasta file to annotate> <output directory to create> <path to v-scan.pl config file>\n";
my $synopsis      = "$executable :: scan and annotate sequences against VADR model libraries ";
my $date          = scalar localtime();
my $version       = "1.7dev0";
my $releasedate   = "Feb 2025";
my $pkgname       = "VADR";

# make *STDOUT file handle 'hot' so it automatically flushes whenever we print to it
select *STDOUT;
$| = 1;

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, undef);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

my $do_verbose   = opt_Get("-v",       \%opt_HH);
my $do_multi     = opt_Get("-m",       \%opt_HH);
my $do_keep      = opt_Get("--keep",   \%opt_HH);

# related to sampling
my $sample_nseq   = opt_Get("--s_nseq", \%opt_HH);
my $do_sample_beg = opt_Get("--s_beg",  \%opt_HH);
my $rand_seed     = opt_Get("--s_seed", \%opt_HH);
my $do_all        = ($do_multi || opt_Get("--all", \%opt_HH)) ? 1 : 0;
my $do_sample     = ($do_all) ? 0 : 1;

# parse config file, we do this early so we can handle -l      
my $config_file = $env_vadr_config_file; # this may be undef
if(opt_IsUsed("-c", \%opt_HH)) {
  $config_file = opt_Get("-c", \%opt_HH);
  utl_FileValidateExistsAndNonEmpty($config_file, "config file specified with -c", undef, 1, undef); # '1' says: die if it doesn't exist or is empty
}
else {
  if(! defined $config_file) {
    die "\nERROR, the environment variable \$VADRCONFIGFILE is not set,\neither set it as the path to the v-scan.pl config file or use the -c option\n";
  }
  utl_FileValidateExistsAndNonEmpty($config_file, "config file defined by env variable \$VADRCONFIGFILE", undef, 1, undef); # '1' says: die if it doesn't exist or is empty
}

my @okey_A = ();      # array of model library keys, read from config file
my %okey_mdir_H   = (); # hash of model directories for each options key, read from config file, key is options key
my %okey_opts_H   = (); # hash of options for each option key, read from config file, key is model key
my %okey_mkey_H   = (); # hash of options for each option key, read from config file, key is model key
my %other_okey_HA = (); # hash of arrays, key is option key $okey, value is array of other okeys ($okey2) that
                        # use $okey as mkey, e.g. $okey = "flavi", @{$other_okey_HA{"flavi"} = ("dengue", "hcv")

parse_config_file($config_file, \@okey_A, \%okey_mdir_H, \%okey_opts_H, \%okey_mkey_H, \%opt_HH, undef);
validate_okey_mkey_values_and_fill_other_okey_HA(\@okey_A, \%okey_mdir_H, \%okey_mkey_H, \%other_okey_HA);

utl_HDump("okey_mkey_H",    \%okey_mkey_H, *STDOUT);
utl_HADump("other_okey_HA", \%other_okey_HA, *STDOUT);
exit 0;

# enforce that --only and --skip options are valid
if((opt_IsUsed("--only", \%opt_HH)) || (opt_IsUsed("--skip", \%opt_HH))) { 
  only_skip_options(\@okey_A, \%okey_mdir_H, \%okey_opts_H, \%opt_HH);
}

# handle --l (list) options, if any of these are selected we just output info and exit
# we do not run v-annotate.pl on any sequences
if(opt_IsUsed("--l_all", \%opt_HH) ||
   opt_IsUsed("--l_lib", \%opt_HH) ||
   opt_IsUsed("--l_dir", \%opt_HH) ||
   opt_IsUsed("--l_opt", \%opt_HH) ||
   opt_IsUsed("--l_mdl", \%opt_HH)) {
  list_options($config_file, \@okey_A, \%okey_mdir_H, \%okey_opts_H, $pkgname, $version, $releasedate, \%opt_HH);
  exit 0;
}

# check that number of command line args is correct
if(scalar(@ARGV) != 2) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do $executable -h\n\n";
  exit(1);
}

my ($orig_in_fa_file, $dir) = (@ARGV);

#############################
# create the output directory
#############################
my $cmd;               # a command to run with utl_RunCommand()
my @early_cmd_A = ();  # array of commands we run before our log file is opened

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
  else                       { die "ERROR a file named $dir already exists. Remove it, or use -f to overwrite it."; }
}

# create the dir
$cmd = "mkdir $dir";
utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root = $dir . "/" . $dir_tail . ".vadr";

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("sequence file", "output directory");
my @arg_A      = ($orig_in_fa_file, $dir);
my %extra_H    = ();
$extra_H{"\$VADRSCRIPTSDIR"}  = $env_vadr_scripts_dir;
$extra_H{"\$VADRCONFIGFILE"}  = (defined $env_vadr_config_file) ? $env_vadr_config_file : "undef";
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
                         # 2D keys (at least initially)
                         #  "log":  log file of what's output to stdout
                         #  "cmd":  command file with list of all commands executed
                         #  "list": file with list of all output files created

# open the log and command files 
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "log",  $out_root . ".log",      1, 1, "Output printed to screen");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "cmd",  $out_root . ".cmd",      1, 1, "List of executed commands");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".filelist", 1, 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
my $FH_HR  = $ofile_info_HH{"FH"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# now we have the log file open, output the banner there too
ofile_OutputBanner($log_FH, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

##############################################
# Validate that we have all the files we need
##############################################
my $progress_w = 60; # the width of the left hand column in our progress output, hard-coded
my $start_secs = ofile_OutputProgressPrior("Validating input", $progress_w, $log_FH, *STDOUT);

utl_FileValidateExistsAndNonEmpty($orig_in_fa_file, "input fasta sequence file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
my @to_remove_A   = (); # list of files to remove at end of subroutine, if --keep not used

ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###########################################
# Copy and validate the input sequence file
###########################################
my $in_fa_file = undef;
if(opt_Get("--origfa", \%opt_HH)) { 
  # --origfa: analyze original fasta file, do not copy it
  $in_fa_file = $orig_in_fa_file;
  if(-e $in_fa_file . ".ssi") { unlink $in_fa_file . ".ssi"}; # remove SSI file if it exists, it may be out of date
}
else { 
  # default: copy original fasta file and analyze that, but don't just copy it, 
  # use 'esl-reformat fasta', this was introduced to sidestep some mysterious 
  # SSI related issues
  $in_fa_file = $out_root . ".in.fa";
  utl_RunCommand($execs_H{"esl-reformat"} . " fasta $orig_in_fa_file > $in_fa_file", opt_Get("-v", \%opt_HH), 0, $FH_HR);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "cp.in.fasta", $in_fa_file, $do_keep, $do_keep, "copy of input fasta file");
  push(@to_remove_A, $in_fa_file);
  push(@to_remove_A, $in_fa_file . ".ssi");
}

my $in_sqfile = Bio::Easel::SqFile->new({ fileLocation => $in_fa_file }); # the sequence file object
my $in_nseq   = $in_sqfile->nseq_ssi;
my $sample_in_fa_file = undef;
my $rand = undef;

# if $do_sample, create the smaller file we'll use for classifying
if($do_sample) {
  if($sample_nseq >= $in_nseq) {
    # num to sample meets or exceeds number of seqs in file, look at all of them in original file
    $sample_in_fa_file = $in_fa_file;
    $sample_nseq = $in_nseq;
  }
  else { # we'll take a subset of all sequences for classification
    ofile_OutputProgressPrior("Sampling $sample_nseq sequences to use for classification", $progress_w, $log_FH, *STDOUT);
    $sample_in_fa_file = $out_root . ".sample.in.fa";
    if($do_sample_beg) {
      $in_sqfile->fetch_consecutive_seqs($sample_nseq, "", 60, $sample_in_fa_file);
    }
    else { # ! $do_sample_beg, sample randomly
      $rand = Bio::Easel::Random->new({ seed => $rand_seed }); # the RNG
      my %chosen_H = ();
      open(FA, ">", $sample_in_fa_file) || ofile_FileOpenFailure($sample_in_fa_file, "v-scan", $!, "writing", $FH_HR);
      my $nchosen = 0;
      my $nrolls  = 0;
      while($nchosen < $sample_nseq) {
        my $j = $rand->roll($in_nseq);
        if(! defined $chosen_H{$j}) {
          print FA $in_sqfile->fetch_seq_to_fasta_string_given_ssi_number($j, 60) . "\n";
          $chosen_H{$j} = 1; # so we don't pick same seq twice
          $nchosen++;
        }
        $nrolls++;
        if($nrolls > (100 * $sample_nseq)) {
          ofile_FAIL("ERROR, unexpectedly taking too many random rolls to pick $sample_nseq seqs, try a different strategy", 1, $FH_HR);
        }
      }
    }
    #push(@to_remove_A, $sample_in_fa_file);
    #push(@to_remove_A, $sample_in_fa_file . ".ssi");
    ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  } # end of else entered if ($sample_nseq < $in_nseq)
} # end of 'if($do_sample)'
else {
  $sample_nseq = $in_nseq;
}

##################################################
# For each model key, run v-annotate.pl --cls_only
##################################################
my %cls_outdir_H = (); # hash of output directories
my %sqc_H = ();     # hash of sqc files
my $okey;
my $n_okey = scalar(@okey_A);
my $clsonly_fa_file = ($do_sample) ? $sample_in_fa_file : $in_fa_file;
my @cls_outdir_A = ();
my $keep_opt   = ($do_keep) ? "--keep" : "";
my $okey_width = 0;
my $mkey_opt2use = "";
foreach $okey (@okey_A) {
  if(length($okey) > $okey_width) { $okey_width = length($okey); }
}
if($n_okey > 1) { 
  foreach $okey (@okey_A) {
    if($okey_mkey_H{$okey} eq $okey) { # we don't need to rerun clsonly mode when okey != mkey
      $cls_outdir_H{$okey} = $dir_tail . "/" . $dir_tail . ".clsonly." . $okey;
      push(@cls_outdir_A, $cls_outdir_H{$okey});
      $sqc_H{$okey} = $cls_outdir_H{$okey} . "/" . $dir_tail . ".clsonly." . $okey . ".vadr.sqc";
      # determine --okey option to use, this is --mkey $mkey unless specified in config file options string
      $mkey_opt2use = "--mkey " . mkey_from_opts($okey, $okey_opts_H{$okey});
      $cmd = $execs_H{"v-annotate.pl"} . " -f -s --origfa --cls_only $mkey_opt2use --mdir $okey_mdir_H{$okey} $keep_opt $clsonly_fa_file $cls_outdir_H{$okey}";
      if(! $do_verbose) { $cmd .= " > /dev/null"; }
      my $start_secs = ofile_OutputProgressPrior(sprintf("Scanning $sample_nseq sequences against %-*s library ", $okey_width, $okey), $progress_w, $log_FH, *STDOUT);
      utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
      ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
    }
  }
}
##################################################
# Parse sqc files to determine which seqs match best to each library, need
# to look at all sqc files before assigning sequences to a mkey
# because we may be determining best library based on score
##################################################
my %seq_H         = ();   # 'exists' hash, key is sequence name, value is always '1' 
my @seq_A         = ();   # array of sequence names
my %seq_okey_H    = ();   # key is seq name, value is best okey for this sequence
my %seq_mdl_H     = ();   # key is seq name, value is best model for this sequence
my %seq_sc_H      = ();   # key is seq name, value is score for best model for this sequence
if($n_okey > 1) { 
  foreach $okey (@okey_A) {
    if($okey_mkey_H{$okey} eq $okey) { 
      parse_sqc_clsonly_file($sqc_H{$okey}, $okey, \%seq_H, \@seq_A, \%seq_okey_H, \%seq_mdl_H, \%seq_sc_H, \%opt_HH, $FH_HR);
    }
  }
}
utl_HDump("seq_okey_H", \%seq_okey_H, *STDOUT);
exit 0;
    
# Fill per-mkey lists of sequences
my %seqlist_HA = (); # key is mkey, value is array of sequences that match to this mkey
my $n_okey_used = 0; # number of mkey (libraries) we have at least one sequence to rerun v-annotate.pl for
my %okey_ct_H  = (); # key is okey, value is number of seqs assigned to that okey, 'undef' if 0
if($n_okey > 1) {
  foreach my $seqname (@seq_A) {
    if(defined $seq_okey_H{$seqname}) {
      my $okey = $seq_okey_H{$seqname};
      my $mdl  = $seq_mdl_H{$seqname};
      if(! defined $seqlist_HA{$okey}) {
        @{$seqlist_HA{$okey}} = ();
        $okey_ct_H{$okey} = 0;
        $n_okey_used++;
      }
      push(@{$seqlist_HA{$okey}}, $seqname);
      $okey_ct_H{$okey}++;
    }
  } 
  
  if((! $do_multi) && ($n_okey_used > 1)) {
    my $okey_str = "";
    foreach $okey (sort keys %seqlist_HA) {
      if($okey_str ne "") { $okey_str .= ", "; }
      $okey_str .= $okey;
    }
    ofile_FAIL("ERROR, -m not used but found matches to multiple libraries: $okey_str", 1, $FH_HR);
  }
}
else {
  $n_okey_used = 1; # we didn't run in clsonly because we only have 1 library
}

###########################################################################
# Re-run v-annotate.pl for each model key that at least one seq matched to
###########################################################################
my @okey_used_A  = (); # array of the okeys with at least one sequence 
my @ant_outdir_A = (); # array of output directories we will create 
my @mdl_file_A   = (); # array of mdl files to output before exiting
my @alc_file_A   = (); # array of alc files to output before exiting
my @log_file_A   = (); # array of log files to process before exiting
my $okey_fa_file = undef;
my $progress_str = undef;
my $ant_okey_width = 0;   # max length of any okey we will annotate for
# first get max width of okey used
foreach $okey (@okey_A) {
  if(defined $seqlist_HA{$okey}) {
    if(length($okey) > $ant_okey_width) { $ant_okey_width = length($okey); }
  }
}
if($n_okey_used > 0) { 
  foreach $okey (@okey_A) {
    if((defined $seqlist_HA{$okey}) || ($n_okey == 1)) { # if $n_okey == 1, we didn't run --clsonly mode
      if($n_okey_used == 1) { 
        $okey_fa_file = $in_fa_file;
        $progress_str = "Annotating $in_nseq sequences with $okey model library ";
      }
      else {
        $okey_fa_file = $dir_tail . "/" . $dir_tail . "." . $okey . ".fa";
        $in_sqfile->fetch_seqs_given_names(\@{$seqlist_HA{$okey}}, 60, $okey_fa_file);
        $progress_str = sprintf("Annotating %*d %-*s sequences ", length($in_nseq), scalar(@{$seqlist_HA{$okey}}), $ant_okey_width, $okey);
      }
      my $ant_outdir = $dir_tail . "/" . $dir_tail . "." . $okey;
      push(@okey_used_A, $okey);
      push(@ant_outdir_A, $ant_outdir);
      push(@mdl_file_A, $ant_outdir . "/" . $dir_tail . "." . $okey . ".vadr.mdl");
      push(@alc_file_A, $ant_outdir . "/" . $dir_tail . "." . $okey . ".vadr.alc");
      push(@log_file_A, $ant_outdir . "/" . $dir_tail . "." . $okey . ".vadr.log");
      $mkey_opt2use = "--mkey " . mkey_from_opts($okey, $okey_opts_H{$okey});
      $cmd = $execs_H{"v-annotate.pl"} . " $mkey_opt2use --mdir $okey_mdir_H{$okey} $okey_opts_H{$okey} $keep_opt $okey_fa_file $ant_outdir";
      if(! $do_verbose) { $cmd .= " > /dev/null"; }
      my $start_secs = ofile_OutputProgressPrior($progress_str, $progress_w, $FH_HR->{"log"}, *STDOUT);
      utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
      ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
    }
  }
}

$start_secs = ofile_OutputProgressPrior("Generating tabular output", $progress_w, $log_FH, *STDOUT);

# create the @data_lib_AA
my $okey_idx = 1;
my $mdl_idx = 1;

# open files for writing
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "lib", $out_root . ".lib", 1, 1, "per-model library tabular summary file");
my @head_lib_AA = ();
my @data_lib_AA = ();
my @clj_lib_A   = ();
@{$head_lib_AA[0]} = ("",    "",        "num");
@{$head_lib_AA[1]} = ("idx", "library", "seqs");
@clj_lib_A         = (1,     1,         0);

foreach $okey (@okey_A) {
  if($okey_mkey_H{$okey} ne $okey) { 
    my $nseq2print = (defined $okey_ct_H{$okey}) ? $okey_ct_H{$okey} : 0;
    if(scalar(@okey_A) == 1) { # we didn't run clsonly mode, set nseq to '-'
      $nseq2print = "-";
    }
    push(@data_lib_AA, [$okey_idx, $okey, $nseq2print]);
    $okey_idx++;
  }
}

ofile_TableHumanOutput(\@data_lib_AA, \@head_lib_AA, \@clj_lib_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"lib"}, undef, $FH_HR);
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

############
# Conclude #
############
output_lib_mdl_and_alc_files_and_remove_temp_files($in_nseq, $sample_nseq, \@okey_used_A, \@mdl_file_A, \@alc_file_A, \@to_remove_A, \%opt_HH, \%ofile_info_HH);

my $z = 0;
if($do_keep) {
  # with --keep leave the files where they are
  for($z = 0; $z < scalar(@okey_A); $z++) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# All %-*s library --clsonly  output files can be found in directory $cls_outdir_A[$z]\n", $okey_width, $okey_A[$z]));
  }
  for($z = 0; $z < scalar(@okey_used_A); $z++) {
    ofile_OutputString($FH_HR->{"log"}, 1, sprintf("# All %-*s library annotation output files can be found in directory $ant_outdir_A[$z]\n", $ant_okey_width, $okey_used_A[$z]));
  }
}
else {
  # --keep not used, move all files in *annotation* subdirectories into $dir and remove subdirs
  # remove clsonly directories
  for($z = 0; $z < scalar(@okey_A); $z++) {
    utl_RunCommand("rm $cls_outdir_A[$z]/*; rmdir $cls_outdir_A[$z]", opt_Get("-v", \%opt_HH), 0, $FH_HR);
  }
  for($z = 0; $z < scalar(@okey_used_A); $z++) {
    # parse the .log file to determine which output files we want to list
    parse_log_file_for_out_files($log_file_A[$z], $okey_used_A[$z], $ant_okey_width, $FH_HR);
    utl_RunCommand("mv $ant_outdir_A[$z]/* ./$dir/; rmdir $ant_outdir_A[$z]", opt_Get("-v", \%opt_HH), 0, $FH_HR);
    if($z < (scalar(@okey_used_A) - 1)) {
      ofile_OutputString($FH_HR->{"log"}, 1, "#\n");
    }
  }
}
if($n_okey_used == 0) { # matches were found to zero libraries
  ofile_OutputString($FH_HR->{"log"}, 1, "# Zero sequences matched a model library so no annotations were performed.\n");
}

$total_seconds += ofile_SecondsSinceEpoch();
ofile_OutputConclusionAndCloseFilesOk($total_seconds, $dir, \%ofile_info_HH);


#################################################################
# Subroutine: parse_config_file
# Incept:     EPN, Thu Feb  6 18:13:40 2025
# 
# Purpose:    Parse the special v-scan.pl config file and store
#             the relevant info in 
# Arguments:
#  $config_file:   path to config file
#  $okey_AR:       REF to array of model library option keys
#  $okey_mdir_HR:  REF to hash, key is okey, value is model directory for this okey, filled here
#  $okey_opts_HR:  REF to hash, key is okey, value is options string to use when annotation for this okey, filled here
#  $okey_mkey_HR:  REF to hash, key is okey, value is --mkey used for annotation, filled here
#  $opt_HHR:       REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $FH_HR:         REF to hash of file handles
#
# Returns:  void
#           
# Dies:     if problem parsing config file
#
#################################################################
sub parse_config_file { 
  my $sub_name = "parse_config_file"; 
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($config_file, $okey_AR, $okey_mdir_HR, $okey_opts_HR, $okey_mkey_HR, $opt_HHR, $FH_HR) = (@_);

  open(CONFIG, $config_file) || ofile_FileOpenFailure($config_file, $sub_name, $!, "reading", $FH_HR);

  my $line;
  my ($okey, $mdir, $opts);
  while($line = <CONFIG>) {
    chomp $line;
    if(($line =~ m/\w/) && ($line !~ m/^\#/)) {
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) < 2) {
        ofile_FAIL("ERROR all non-comment lines should include at least two white space delimited fields: <outkey> <modeldir>\nread line:\n$line", 1, $FH_HR);
      }
      my ($okey, $mdir) = ($el_A[0], $el_A[1]);
      my $test_okey = $okey;
      $test_okey =~ s/[^a-z0-9]//g;
      if($test_okey ne $okey) { 
        ofile_FAIL("ERROR ready okey $okey, which includes some characters that are not lowercase or numeric, all okey values in field 1 must be all lowercase without any special non-alphanumeric characters", 1, $FH_HR);
      }
      my $opts = "";
      for(my $i = 2; $i < scalar(@el_A); $i++) {
        if($opts ne "") { $opts .= " "; }
        $opts .= $el_A[$i];
      }
      if(defined $okey_mdir_H{$okey}) {
        ofile_FAIL("ERROR read output model key $okey twice in config file", 1, $FH_HR);
      }
      push(@{$okey_AR}, $okey);
      $okey_mdir_HR->{$okey} = $mdir;
      $okey_opts_HR->{$okey} = $opts;
      $okey_mkey_HR->{$okey}  = mkey_from_opts($okey, $opts);
    }
  }
  close(CONFIG);

  # make sure all %{$okey_mkey_HR} values are valid mkeys
  foreach $okey (@{$okey_AR}) {
    if(! defined $okey_mdir_HR->{($okey_mkey_HR->{$okey})}) {
      ofile_FAIL("ERROR, in config file for okey $okey, options string includes --mkey $okey_mkey_HR->{$okey}\nbut $okey_mkey_HR->{$okey} does not have an entry in the config file.", 1, $FH_HR);
    }
  }
  return;
}

#################################################################
# Subroutine: parse_sqc_clsonly_file
# Incept:     EPN, Thu Feb  6 18:42:47 2025
# 
# Purpose:    Parse the .sqc file output from v-annotate.pl --cls_only
#
# Arguments:
#  $sqc_file:      name of sqc file to parse
#  $okey:          options key (e.g. flu) that this sqc file pertains to
#  $other_okey_AR: REF to array of other option keys that use this model library
#                  e.g. if $okey is 'flavi', @{$other_okey_AR} might be ('dengue', 'hcv')
#  $seq_HR:        REF to hash of sequence names, key is seq name, value is 1, to fill here
#  $seq_AR:        REF to array of sequence names, to fill here
#  $seq_okey_HR:   REF to hash, key is sequence name, value is winning mkey, to fill here
#  $seq_mdl_HR:    REF to hash, key is sequence name, value is winning model, to fill here
#  $seq_sc_HR:     REF to hash, key is sequence name, value is winning score, to fill here
#  $opt_HHR:       REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $FH_HR:         REF to hash of file handles
#
# Returns:  void
#           
# Dies:     if problem parsing sqc file
#
#################################################################
sub parse_sqc_clsonly_file { 
  my $sub_name = "parse_sqc_clsonly_file"; 
  my $nargs_exp = 10;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqc_file, $okey, $other_okey_AR, $seq_HR, $seq_AR, $seq_okey_HR, $seq_mdl_HR, $seq_sc_HR, $opt_HHR, $FH_HR) = (@_);

  my $do_lone  = opt_Get("--lone", $opt_HHR);
  my $do_first = opt_Get("--first", $opt_HHR);

  open(SQC, $sqc_file) || ofile_FileOpenFailure($sqc_file, $sub_name, $!, "reading", $FH_HR);

  #seq  seq                               seq                             sub                     seq  mdl         num                       sub     score  diff/  seq   
  #idx  name                              len  p/f   ant  model1    grp1  grp1    score  sc/nt    cov  cov  bias  hits  str  model2    grp2  grp2     diff     nt  alerts
  #---  -----------------------------  ------  ----  ---  --------  ----  ----  -------  -----  -----  ---  ----  ----  ---  --------  ----  ----  -------  -----  ------
  #1     gi|1205299962|gb|MF147925.1|     1765  FAIL  no   -         -     -           -      -      -    -     -     -    _  -         -     -           -      -  NO_ANNOTATION(noannotn)
  #10    KY654518                        15277  PASS  no   KY654518  RSV   A     28567.0  1.861  1.005    -     -     2    +  MZ516105  RSV   B     16488.0  0.791  -     
  while(my $line = <SQC>) { 
    if($line !~ m/^\#/) {
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) != 21) { 
        ofile_FAIL("ERROR problem parsing sqc file $sqc_file", 1, $FH_HR);
      }
      my ($seqname, $pf, $mdl, $grp, $subgrp, $score) = ($el_A[1], $el_A[3], $el_A[5], $el_A[6], $el_A[7], $el_A[8]);
      my $okey2use = $okey;
      my $mdl2use = $mdl;
      if(defined $other_okey_AR) {
        $grp     =~ tr/A-Z/a-z/;
        $subgrp  =~ tr/A-Z/a-z/;
        $mdl2use =~ tr/A-Z/a-z/;
        $grp     =~ s/[^a-z0-9]//g;
        $subgrp  =~ s/[^a-z0-9]//g;
        $mdl2use =~ s/[^a-z0-9]//g;
        foreach my $other_okey (@{$other_okey_AR}) {
          # other_okey will be lowercase without special characters, parse_config_file makes sure of this
          if($mdl2use eq $other_okey) { 
            $okey2use = $mdl2use;
          }
          elsif($grp eq $other_okey) { 
            $okey2use = $grp;
          }
          elsif($subgrp eq $other_okey) { 
            $okey2use = $subgrp
          }
        }
      }
      
      if(! defined $seq_HR->{$seqname}) {
        push(@seq_A, $seqname);
        $seq_HR->{$seqname} = 1;
      }
      if($el_A[3] eq "PASS") {
        my $keep_flag = 1;
        if(defined $seq_okey_HR->{$seqname}) {
          # this sequence already matched a model for a different $okey
          # we either:
          # 1) die with error message
          # 2) figure out best okey for this sequence
          #    either first okey seen, or okey that gave top scoring hit
          if($do_lone) {
            ofile_FAIL("ERROR sequence $seqname matched to two libraries: $seq_okey_HR->{$seqname} and $okey, omit --lone to allow this", 1, $FH_HR);
          }
          if($do_first) {
            $keep_flag = 0; # keep existing value in $seq_okey_HR->{$seqname}
          }
          else {
            # does it score better? it has to be better by at least 1
            # this means if you have identical models in multiple libraries,
            # list the library that you prefer to use earliest in the config file
            $keep_flag = (($score-1.) > $seq_sc_HR->{$seqname}) ? 1 : 0;
          }
        }
        if($keep_flag) {
          $seq_okey_HR->{$seqname} = $okey2use;
          $seq_mdl_HR->{$seqname}  = $mdl;
          $seq_sc_HR->{$seqname}   = $score;
        }
      }
    }
  }
  close(SQC);

  return;
}

#################################################################
# Subroutine: output_lib_mdl_and_alc_files_and_remove_temp_files()
# Incept:     EPN, Mon Feb 10 14:13:25 2025
#             based on v-annotate.pl:output_mdl_and_alc_files_and_remove_temp_files()
# Purpose:    Output the lib file and then for any library
#             with at least one sequence annotated, output the
#             mdl and alc files, then remove all files
#             in (@{$to_remove_A}), unless --keep. 
#
# Arguments:
#  $in_nseq;        number of sequences in input file
#  $sample_nseq:    number of sequences sampled
#  $okey_used_AR:   ref to array of option keys we want to output .mdl and .alc files for
#  $mdl_file_AR:    ref to array of .mdl files to output
#  $alc_file_AR:    ref to array of .alc files to output
#  $to_remove_AR:   ref to array of files to remove
#  $opt_HHR:        ref to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR: ref to 2D hash of output file information, added to here
#             
# Returns:  void
#
#################################################################
sub output_lib_mdl_and_alc_files_and_remove_temp_files { 
  my $sub_name = "output_lib_mdl_and_alc_files_and_remove_temp_files";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($in_nseq, $sample_nseq, $okey_used_AR, $mdl_file_AR, $alc_file_AR, $to_remove_AR, $opt_HHR, $ofile_info_HHR) = (@_);

  # close the two files we may output to stdout and the log
  close($ofile_info_HHR->{"FH"}{"lib"});
  
  my $FH_HR     = $ofile_info_HH{"FH"};
  my $do_multi  = opt_Get("-m", $opt_HHR);
  my $do_all    = ($do_multi || opt_Get("--all", $opt_HHR)) ? 1 : 0;
  my $do_sample = ($do_all) ? 0 : 1;
  my $sum_str   = "";

  # only output the lib file if -m was used 
  my @conclude_A = ();
  my @file_A = ();
  my ($okey, $mdl_file, $alc_file) = (undef, undef, undef);
  my $n_okey = scalar(@{$okey_used_AR});
  if($do_multi) { 
    if(($do_sample) && ($sample_nseq < $in_nseq)) {
      $sum_str = sprintf("# Summary of seqs matching each library (only %d of %d seqs scanned):", $sample_nseq, $in_nseq);
    }
    else {
      $sum_str = "# Summary of sequences matching each library:";
    }
    
    if($n_okey != scalar(@{$mdl_file_AR})) {
      ofile_FAIL("ERROR, in $sub_name, unexpected number of mdl files", 1, $FH_HR);
    }
    if(scalar(@{$mdl_file_AR}) != (scalar(@{$alc_file_AR}))) {
      ofile_FAIL("ERROR, in $sub_name, number of mdl and alc files differ", 1, $FH_HR);
    }
    
    push(@conclude_A, "#");
    push(@conclude_A, $sum_str);
    push(@conclude_A, "#");
    utl_FileLinesToArray($ofile_info_HHR->{"fullpath"}{"lib"}, 1, \@file_A, $FH_HR);
    push(@conclude_A, @file_A);
    push(@conclude_A, "#");
  }
  
  for(my $m = 0; $m < $n_okey; $m++) {
    $okey     = $okey_used_AR->[$m];
    $mdl_file = $mdl_file_AR->[$m];
    $alc_file = $alc_file_AR->[$m];

    push(@conclude_A, "#");
    push(@conclude_A, "# Summary of sequences matching $okey:");
    push(@conclude_A, "#");

    @file_A = ();
    utl_FileLinesToArray($mdl_file, 1, \@file_A, $FH_HR);
    push(@conclude_A, @file_A);
    push(@conclude_A, "#");

    @file_A = ();
    utl_FileLinesToArray($alc_file, 1, \@file_A, $FH_HR);
    if(scalar(@file_A == 3)) {
      push(@conclude_A, "# Zero alerts reported for seqs matching $okey.");
    }
    else {
      push(@conclude_A, "# Summary of reported alerts for seqs matching $okey:");
      push(@conclude_A, "#");
      push(@conclude_A, @file_A);
    }
    push(@conclude_A, "#");
  }

  foreach my $line (@conclude_A) { 
    ofile_OutputString($FH_HR->{"log"}, 1, $line . "\n");
  }
  
  # remove unwanted files, unless --keep
  if(! opt_Get("--keep", $opt_HHR)) { 
    my @to_actually_remove_A = (); # sanity check: make sure the files we're about to remove actually exist
    my %to_actually_remove_H = (); # sanity check: to make sure we don't try to delete 
    foreach my $to_remove_file (@to_remove_A) { 
      if((defined $to_remove_file) && (-e $to_remove_file) && (! defined $to_actually_remove_H{$to_remove_file})) { 
        push(@to_actually_remove_A, $to_remove_file); 
        $to_actually_remove_H{$to_remove_file} = 1; 
      }
    }
    utl_FileRemoveList(\@to_actually_remove_A, $sub_name, $opt_HHR, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine:  only_skip_options()
# Incept:      EPN, Wed Feb 12 10:48:27 2025
#
# Purpose:    Handle the --only and --skip options by 
#             parsing their strings, determining if they are valid
#             and updating the @{$okey_AR}, %{$okey_mdir_HR}
#             and %{$okey_opts_HR} data structures.
#
# Arguments: 
#  $okey_AR:      REF to array of all okeys read from config file, modified here
#  $okey_mdir_HR: REF to hash of directories for each option key, modified here
#  $okey_opts_HR: REF to hash of options for each option key, modified here
#  $okey_mkey_HR: REF to hash of mkeys to use for each option key, modified here
#  $opt_HHR:      REF to 2D hash of option values
#
# Returns:    void
#
# Dies:       if --only or --skip option strings are invalid
#
#################################################################
sub only_skip_options { 
  my $sub_name = "only_skip_options()"; 
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($okey_AR, $okey_mdir_HR, $okey_opts_HR, $okey_mkey_HR, $opt_HHR) = @_;
  
  my @only_A = ();  # array of values in --only arg
  my @skip_A = ();  # array of values in --skip arg
  my %only_H = ();  # 'exists' hash for values in --only arg
  my %skip_H = ();  # 'exists' hash for values in --skip arg
  my @new_okey_A = ();      # new array of mkeys we will replace @{$mkey_AR} with before returning
  my %new_okey_mdir_H = (); # new hash we will replace %{$mkey_mdir_HR} with before returning
  my %new_okey_opts_H = (); # new hash we will replace %{$mkey_opts_HR} with before returning
  my %new_okey_mkey_H = (); # new hash we will replace %{$mkey_mkey_HR} with before returning

  my $die_str = "";
  if(opt_IsUsed("--only", $opt_HHR)) { 
    @only_A = split(",", opt_Get("--only", $opt_HHR));
    foreach my $only_okey (@only_A) {
      $only_H{$only_okey} = 1;
      if(! defined $okey_mdir_HR->{$only_okey}) {
        $die_str .= "\t$only_okey specified in --only option but not listed in config file\n";
      }
    }
    foreach my $okey (@{$okey_AR}) {
      if(defined $only_H{$okey}) {
        push(@new_okey_A, $okey);
        $new_okey_mdir_H{$okey} = $okey_mdir_H{$okey};
        $new_okey_opts_H{$okey} = $okey_opts_H{$okey};
        $new_okey_mkey_H{$okey} = $okey_mkey_H{$okey};
      }
    }
  }
  if(opt_IsUsed("--skip", $opt_HHR)) { 
    if(scalar(@only_A) != 0) {
      # this should have been enforced by opt_ValidateSet() 
      ofile_FAIL("ERROR, in $sub_name, --only and --skip both used, pick one", 1, $FH_HR);
    }
    @skip_A = split(",", opt_Get("--skip", $opt_HHR));
    foreach my $skip_okey (@skip_A) {
      $skip_H{$skip_okey} = 1;
      if(! defined $okey_mdir_HR->{$skip_okey}) {
        $die_str .= "\t$skip_okey specified in --skip option but not listed in config file\n";
      }
    }
    foreach my $okey (@{$okey_AR}) {
      if(! defined $skip_H{$okey}) {
        push(@new_okey_A, $okey);
        $new_okey_mdir_H{$okey} = $okey_mdir_H{$okey};
        $new_okey_opts_H{$okey} = $okey_opts_H{$okey};
        $new_okey_mkey_H{$okey} = $okey_mkey_H{$okey};
      }
    }
  }

  if($die_str ne "") {
      ofile_FAIL("ERROR, in $sub_name:\n$die_str\n", 1, $FH_HR);
  }
  
  @{$okey_AR} = ();
  %{$okey_mdir_HR} = ();
  %{$okey_opts_HR} = ();
  %{$okey_mkey_HR} = ();

  # copy values
  @{$okey_AR} = @new_okey_A;
  foreach $okey (@{$okey_AR}) {
    $okey_mdir_HR->{$okey} = $new_okey_mdir_H{$okey};
    $okey_opts_HR->{$okey} = $new_okey_opts_H{$okey};
    $okey_opts_HR->{$okey} = $new_okey_mkey_H{$okey};
  }
  
  return;
}


#################################################################
# Subroutine:  list_options()
# Incept:      EPN, Wed Feb 19 16:50:35 2025
#
# Purpose:    Handle the --l_* options by printing out information
#             on the model libraries and exiting.
#
# Arguments: 
#  $config_file:  path to config file
#  $okey_AR:      REF to array of all mkeys read from config file, modified here
#  $okey_mdir_HR: REF to hash of directories for each output key, modified here
#  $okey_opts_HR: REF to hash of options for each output key, modified here
#  $okey_mkey_HR: REF to hash of mkeys for each output key, modified here
#  $pkgname:      package name
#  $version:      version
#  $releasedate:  release date for the package
#  $opt_HHR:      REF to 2D hash of option values
#
# Returns:    void
#
# Dies:       if --only or --skip option strings are invalid
#
#################################################################
sub list_options { 

  my $sub_name = "list_options()"; 
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($config_file, $okey_AR, $okey_mdir_HR, $okey_opts_HR, $okey_mkey_HR, $pkgname, $version, $releasedate, $opt_HHR) = @_;

  my $div_line = utl_StringMonoChar(60, "#", undef) . "\n";

  print $div_line;
  print "#\n";
  print "# $pkgname $version ($releasedate)\n";
  print "#\n";
  print "# config file: $config_file\n";
  print "#\n";

  my $do_lib = opt_IsUsed("--l_lib", $opt_HHR) ? 1 : 0;
  my $out_lib = undef;
  my $okey = undef;

  # if --l_lib doesn't exist, exit
  if($do_lib)  {
    $out_lib = opt_Get("--l_lib", $opt_HHR);
    if(! defined $okey_mdir_HR->{$out_lib}) {
      my $die_str = "ERROR, library key $out_lib specified with --l_lib does not exist in config file.\nExisting library keys are:\n";
      foreach $okey (@okey_A) {
        $die_str .= "\t$okey\n";
      }
      die $die_str;
    }
  }
    
  my @head_AA = ();
  my @data_AA = ();
  my @clj_A   = ();

  # model dir table:
  if(opt_IsUsed("--l_all", $opt_HHR) ||
     opt_IsUsed("--l_dir", $opt_HHR) || 
     opt_IsUsed("--l_lib", $opt_HHR)) { 
    if($do_lib) { 
      #print("# $out_lib library directory information:\n");
    }
    else { 
      print("# Model library directory information:\n#\n");
    }

    @{$head_AA[0]} = ("model key", "model dir");
    @clj_A         = (1,     1);
    foreach my $okey (@{$okey_AR}) {
      if((! $do_lib) || ($okey eq $out_lib)) { 
        push(@data_AA, [$okey, $okey_mdir_HR->{$okey}]);
      }
    }
    ofile_TableHumanOutput(\@data_AA, \@head_AA, \@clj_A, undef, undef, "  ", "-", "#", "#", "", 0, *STDOUT, undef, undef);
    print("#\n");
  }

  # model options table:
  if(opt_IsUsed("--l_all", $opt_HHR) ||
     opt_IsUsed("--l_opt", $opt_HHR) ||
     opt_IsUsed("--l_lib", $opt_HHR)) { 
    if($do_lib) { 
      #print("# $out_lib library options information:\n");
    }
    else {
      print("# Model library options information:\n#\n");
    }

    @data_AA = ();
    @{$head_AA[0]} = ("model key", "v-annotate.pl options");
    @clj_A         = (1,     1);
    foreach my $okey (@{$okey_AR}) {
      if((! $do_lib) || ($okey eq $out_lib)) { 
        push(@data_AA, [$okey, $okey_opts_HR->{$okey}]);
      }
    }
    ofile_TableHumanOutput(\@data_AA, \@head_AA, \@clj_A, undef, undef, "  ", "-", "#", "#", "", 0, *STDOUT, undef, undef);
    print("#\n");
  }

  # models table
  my @reqd_mdl_keys_A = ("name", "length");
  if(opt_IsUsed("--l_all", $opt_HHR) ||
     opt_IsUsed("--l_mdl", $opt_HHR) ||
     opt_IsUsed("--l_lib", $opt_HHR)) { 
    if($do_lib) { 
      #print("# List of models in $out_lib model library:\n");
    }
    else {
      print("# List of models in each model library:\n#\n");
    }

    @data_AA = ();
    @{$head_AA[0]} = ("idx", "model key", "model name", "length", "group", "subgroup");
    @clj_A         = (0,     1,           1,            0,        1,       1);
    my @reqd_ftr_keys_A = ("type", "coords");
    my $minfo_file = undef;
    my @mdl_info_AH = ();
    my %ftr_info_HAH = ();
    for(my $k = 0; $k < scalar(@{$okey_AR}); $k++) {
      my $okey = $okey_A[$k];
      $minfo_file = $okey_mdir_HR->{$okey} . "/" . $okey . ".minfo";
      @mdl_info_AH = ();
      %ftr_info_HAH = ();
      utl_FileValidateExistsAndNonEmpty($minfo_file, "$okey model info file", undef, 1, undef);
      vdr_ModelInfoFileParse($minfo_file, \@reqd_mdl_keys_A, \@reqd_ftr_keys_A, \@mdl_info_AH, \%ftr_info_HAH, undef);
      my $nmdl = scalar(@mdl_info_AH);
      if((! $do_lib) && ($k > 0)) { push(@data_AA, []); } # blank line
      for(my $m = 0; $m < $nmdl; $m++) {
        if((! $do_lib) || ($okey eq $out_lib)) { 
          push(@data_AA,
               [(sprintf("%d.%d", ($do_lib ? 1 : ($k+1)), ($m+1))), 
                $okey,
                $mdl_info_AH[$m]{"name"},
                $mdl_info_AH[$m]{"length"},
                ((defined $mdl_info_AH[$m]{"group"})    ? $mdl_info_AH[$m]{"group"} : "-"), 
                ((defined $mdl_info_AH[$m]{"subgroup"}) ? $mdl_info_AH[$m]{"subgroup"} : "-")]);
        }
      }
    }    
    ofile_TableHumanOutput(\@data_AA, \@head_AA, \@clj_A, undef, undef, "  ", "-", "#", "#", "", 0, *STDOUT, undef, undef);
  }

  return;
}

#################################################################
# Subroutine:  parse_log_file_for_out_files()
# Incept:      EPN, Thu Feb 20 13:30:54 2025
#
# Purpose:    Parse a log file and output the lines that list
#             the output files that were created.
#
# Arguments: 
#  $log_file:     path to config file
#  $FH_HR:        REF to hash of file handles
#
# Returns:    void
#
# Dies:       if there's a problem parsing the log file
#
#################################################################
sub parse_log_file_for_out_files { 
  my $sub_name = "parse_log_file_for_out_files";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($log_file, $okey, $okey_width, $FH_HR) = (@_);
  #
  # Zero alerts were reported.
  #
  # Output printed to screen saved in:                              vs-test.flu.vadr.log
  # List of executed commands saved in:                             vs-test.flu.vadr.cmd
  # List and description of all output files saved in:              vs-test.flu.vadr.filelist
  # esl-seqstat -a output for input fasta file saved in:            vs-test.flu.vadr.seqstat
  # 5 column feature table output for passing sequences saved in:   vs-test.flu.vadr.pass.tbl
  # 5 column feature table output for failing sequences saved in:   vs-test.flu.vadr.fail.tbl
  # list of passing sequences saved in:                             vs-test.flu.vadr.pass.list
  # list of failing sequences saved in:                             vs-test.flu.vadr.fail.list
  # list of alerts in the feature tables saved in:                  vs-test.flu.vadr.alt.list
  # alignment doctoring tabular summary file saved in:              vs-test.flu.vadr.dcr
  # replaced stretches of Ns summary file (-r) saved in:            vs-test.flu.vadr.rpn
  #
  # All output files created in directory ./vs-test/vs-test.flu/
  open(LOG, $log_file) || ofile_FileOpenFailure($log_file, $sub_name, $!, "reading", $FH_HR);

  my $print_flag = 0;
  my $did_print = 0;
  while(my $line = <LOG>) {
    chomp $line;
    if($line =~ m/^# Output printed to screen saved in/) {
      $print_flag = 1;
    }
    if(($print_flag) && ($line =~ m/^#\s*$/)) {
      $print_flag = 0;
    }
    if($print_flag) {
      $line =~ s/^# //;
      $line = sprintf("# %-*s library ", $okey_width, $okey) . lcfirst($line);
      ofile_OutputString($FH_HR->{"log"}, 1, $line . "\n");
      $did_print = 1;
    }
  }
  close(LOG);

  if(! $did_print) {
    ofile_FAIL("ERROR in $sub_name, unable to find any output files listed in $log_file\n", 1, $FH_HR);
  }
  return;
}

#################################################################
# Subroutine:  mkey_from_opts()
# Incept:      EPN, Fri Feb 21 14:10:57 2025
#
# Purpose:    Return the mkey set in an options string, if none
#             is set, return $mkey.
#
# Arguments: 
#  $mkey:     path to config file
#  $opts:     REF to hash of file handles
#
# Returns:    void
#
# Dies:       if there's a problem parsing the log file
#
#################################################################
sub mkey_from_opts {
  my $sub_name = "mkey_from_opts";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mkey, $opts) = (@_);

  if($opts =~ /\s+\-\-mkey\s+(\S+)/) {
    return $1;
  }
  return $mkey;
}

#################################################################
# Subroutine:  validate_okey_mkey_values_and_fill_other_okey_HA()
# Incept:      EPN, Fri Feb 21 14:25:32 2025
#
# Purpose:    For any $okey_mkey_HR->{$okey}=$mkey values that do not equal
#             $okey, read the minfo file for $mkey and make sure at
#             least one model exists with a name, group or subgroup
#             that, when lowercased and had specials removed, equals
#             $okey. If not, exit with an error.
#
# Arguments: 
#  $okey_AR:        ref to array of okeys
#  $okey_mdir_HR:   hash with mdir values for each okey
#  $okey_mkey_HR:   hash with mkey used for classifying/annotating
#                   for this okey
#  $other_okey_HAR: ref to hash of arrays, key is $okey, value is
#                   array of all $okey2 != $okey for which
#                   $okey_mkey_HR{$okey2} = $okey. For example
#                   if $okey is "flavi", $other_okey_HAR->{"flavi"}
#                   could be ("dengue", "hcv").
#
# Returns:    void
#
# Dies:       see 'Purpose'
#
#################################################################
sub validate_okey_mkey_values_and_fill_other_okey_HA {
  my $sub_name = "validate_okey_mkey_values_and_fill_other_okey_HA";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($okey_AR, $okey_mdir_HR, $okey_mkey_HR, $other_okey_HAR) = (@_);

  my $k; 
  my $die_str = "";
  for($k = 0; $k < scalar(@{$okey_AR}); $k++) {
    my $okey = $okey_A[$k];
    my $mkey = $okey_mkey_HR->{$okey};
    if($okey ne $mkey) { 
      if(! defined $other_okey_HAR->{$mkey}) {
        $other_okey_HAR->{$mkey} = ();
      }
      push(@{$other_okey_HAR->{$mkey}}, $okey);
    }
  }

  my @mdl_info_AH = ();
  my %ftr_info_HAH = ();
  my @reqd_mdl_keys_A = ("name", "length");
  my @reqd_ftr_keys_A = ("type", "coords");
  my $other_okey;
  my %found_match_H = (); # key is $okey2 from @{$other_okey_HAR->{$okey}}, value is 1 if found match to $okey2
  for($k = 0; $k < scalar(@{$okey_AR}); $k++) {
    if(defined $other_okey_HAR->{$okey}) {
      # we need to find at least one model with name/group/subgroup that matches $other_okey
      foreach $other_okey (@{$other_okey_HAR->{$okey}}) {
        $found_match_H{$other_okey} = 0;
      }

      my $minfo_file = $okey_mdir_HR->{$okey} . "/" . $okey . ".minfo";
      @mdl_info_AH = ();
      %ftr_info_HAH = ();
      utl_FileValidateExistsAndNonEmpty($minfo_file, "$okey model info file", undef, 1, undef);
      vdr_ModelInfoFileParse($minfo_file, \@reqd_mdl_keys_A, \@reqd_ftr_keys_A, \@mdl_info_AH, \%ftr_info_HAH, undef);
      my $nmdl = scalar(@mdl_info_AH);
      for(my $m = 0; $m < $nmdl; $m++) {
        # this model may only match to 0 or 1 of the $other_okey_HAR->{$okey} values
        my $nmatch_this_model = 0;
        my $okey_matches_this_model = "";
        my $mdl2use    = $mdl_info_AH[$m]{"name"};
        my $grp2use    = (defined $mdl_info_AH[$m]{"group"})    ? $mdl_info_AH[$m]{"group"}    : undef; 
        my $subgrp2use = (defined $mdl_info_AH[$m]{"subgroup"}) ? $mdl_info_AH[$m]{"subgroup"} : undef; 
        $mdl2use =~ tr/A-Z/a-z/;
        $mdl2use =~ s/[^a-z0-9]//g;
        if(defined $grp2use) {
          $grp2use =~ tr/A-Z/a-z/;
          $grp2use =~ s/[^a-z0-9]//g;
        }
        if(defined $subgrp2use) { 
          $subgrp2use =~ tr/A-Z/a-z/;
          $subgrp2use =~ s/[^a-z0-9]//g;
        }
        foreach $other_okey (@{$other_okey_HAR->{$okey}}) {
          if(($mdl2use    eq $okey) ||
             ((defined $grp2use)    && ($grp2use    eq $okey)) ||
             ((defined $subgrp2use) && ($subgrp2use eq $okey))) { 
            if($nmatch_this_model > 0) { 
              $nmatch_this_model++;
              if($okey_matches_this_model ne "") { $okey_matches_this_model .= ", "; }
              $okey_matches_this_model .= $other_okey;
            }
            $found_match_H{$other_okey} = 1;
          }
        }
        if($nmatch_this_model > 1) {
          $die_str .= sprintf("Model %s with group: %s and subgroup: %s in library $okey\nmatches more than one other mkey that uses $okey library: $okey_matches_this_model\n\n",
                              $mdl_info_AH[$m]{"name"},
                              (defined $mdl_info_AH[$m]{"group"}    ? $mdl_info_AH[$m]{"group"}    : "undef"),
                              (defined $mdl_info_AH[$m]{"subgroup"} ? $mdl_info_AH[$m]{"subgroup"} : "undef"));
        }
      } # end of (for(my $m = 0; $m < $nmdl; $m++)
    } # end of if(defined $other_okey_HAR->{$okey}
  } # end of for ($k = 0; $k < scalar(@{$okey_AR})    

  if($die_str ne "") {
    die "ERROR in $sub_name:\n$die_str"; 
  }
  
  return;
}

