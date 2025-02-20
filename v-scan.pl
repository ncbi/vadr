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
$opt_group_desc_H{++$g} = "options for choosing a model library based on only a subset of input sequences:";
opt_Add("-p",           "boolean", 0,         $g,   undef,  "-m",       "peek only at a few seqs for picking model library to use",                               "peek only at a few seqs for picking model library to use",   \%opt_HH, \@opt_order_A);
opt_Add("--p_nseq",     "integer", 3,         $g,    "-p", undef,       "with -p, set the number of sequences to peek at to <n>",                                  "with -p, set the number of sequences to peek at to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--p_beg",      "boolean", 0,         $g,    "-p", undef,       "with -p, do not select seqs randomly, use seqs from beginning of file",                   "with -p, do not select seqs randomly, use seqs from beginning of file", \%opt_HH, \@opt_order_A);
opt_Add("--p_seed",     "integer", 181,       $g,    "-p","--p_beg",    "with -p, set the random number generator seed to <n>",                                    "with -p, set the random number generator seed to <n>", \%opt_HH, \@opt_order_A);
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
                'p'        => \$GetOptions_H{"-p"}, 
                'p_nseq=s' => \$GetOptions_H{"--p_nseq"},
                'p_beg'    => \$GetOptions_H{"--p_beg"},
                'p_seed=s' => \$GetOptions_H{"--p_seed"},
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
my $do_peek      = opt_Get("-p",       \%opt_HH);
my $peek_nseq    = opt_Get("--p_nseq", \%opt_HH);
my $do_peek_beg  = opt_Get("--p_beg",  \%opt_HH);
my $rand_seed    = opt_Get("--p_seed", \%opt_HH);

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

my @mkey_A = ();      # array of model library keys, read from config file
my %mkey_mdir_H = (); # hash of model directories for each model library key, read from config file, key is model key
my %mkey_opts_H = (); # hash of options for each model library key, read from config file, key is model key
parse_config_file($config_file, \@mkey_A, \%mkey_mdir_H, \%mkey_opts_H, \%opt_HH, undef);

# enforce that --only and --skip options are valid
if((opt_IsUsed("--only", \%opt_HH)) || (opt_IsUsed("--skip", \%opt_HH))) { 
  only_skip_options(\@mkey_A, \%mkey_mdir_H, \%mkey_opts_H, \%opt_HH);
}

# handle --l (list) options, if any of these are selected we just output info and exit
# we do not run v-annotate.pl on any sequences
if(opt_IsUsed("--l_all", \%opt_HH) ||
   opt_IsUsed("--l_lib", \%opt_HH) ||
   opt_IsUsed("--l_dir", \%opt_HH) ||
   opt_IsUsed("--l_opt", \%opt_HH) ||
   opt_IsUsed("--l_mdl", \%opt_HH)) {
  list_options($config_file, \@mkey_A, \%mkey_mdir_H, \%mkey_opts_H, $pkgname, $version, $releasedate, \%opt_HH);
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
my $peek_in_fa_file = undef;
my $rand = undef;

# if $do_peek, create the smaller file we'll 
if($do_peek) {
  if($peek_nseq >= $in_nseq) {
    # num to peek meets or exceeds number of seqs in file, look at all of them in original file
    $peek_in_fa_file = $in_fa_file;
    $peek_nseq = $in_nseq;
  }
  else { # we'll take a subset of all files
    $peek_in_fa_file = $out_root . ".peek.in.fa";
    if($do_peek_beg) {
      $in_sqfile->fetch_consecutive_seqs($peek_nseq, "", 60, $peek_in_fa_file);
    }
    else { # ! $do_peek_beg
      $rand = Bio::Easel::Random->new({ seed => $rand_seed }); # the RNG
      my %chosen_H = ();
      open(FA, ">", $peek_in_fa_file) || ofile_FileOpenFailure($peek_in_fa_file, "v-scan", $!, "writing", $FH_HR);
      my $nchosen = 0;
      my $nrolls  = 0;
      while($nchosen < $peek_nseq) {
        my $j = $rand->roll($in_nseq);
        if(! defined $chosen_H{$j}) {
          print FA $in_sqfile->fetch_seq_to_fasta_string_given_ssi_number($j, 60) . "\n";
          $chosen_H{$j} = 1; # so we don't pick same seq twice
          $nchosen++;
        }
        $nrolls++;
        if($nrolls > (100 * $peek_nseq)) {
          ofile_FAIL("ERROR, unexpectedly taking too many random rolls to pick $peek_nseq seqs, try a different strategy", 1, $FH_HR);
        }
      }
    }
    #push(@to_remove_A, $peek_in_fa_file);
    #push(@to_remove_A, $peek_in_fa_file . ".ssi");
  } # end of else entered if ($peek_nseq < $in_nseq)
}

##################################################
# For each model key, run v-annotate.pl --cls_only
##################################################
my %out_dir_H = (); # hash of output directories
my %sqc_H = ();     # hash of sqc files
my $mkey;
my $nmkey = scalar(@mkey_A);
my $clsonly_fa_file = ($do_peek) ? $peek_in_fa_file : $in_fa_file;
if($nmkey > 1) { 
  foreach $mkey (@mkey_A) {
    $out_dir_H{$mkey} = $dir_tail . "/" . $mkey . ".0";
    $sqc_H{$mkey} = $out_dir_H{$mkey} . "/" . $mkey . ".0.vadr.sqc";
    $cmd = $execs_H{"v-annotate.pl"} . " -f -s --origfa --cls_only --mkey $mkey --mdir $mkey_mdir_H{$mkey} $clsonly_fa_file $out_dir_H{$mkey}";
    if(! $do_verbose) { $cmd .= " > /dev/null"; }
    my $start_secs = ofile_OutputProgressPrior(sprintf("Scanning sequences against %s library ... ", $mkey), $progress_w, $log_FH, *STDOUT);
    utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
    ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
}
##################################################
# Parse sqc files to determine which seqs match best to each library, need
# to look at all sqc files before assigning sequences to a mkey
# because we may be determining best library based on score
##################################################
my %seq_H         = ();   # 'exists' hash, key is sequence name, value is always '1' 
my @seq_A         = ();   # array of sequence names
my %seq_mkey_H    = ();   # key is seq name, value is best mkey for this sequence
my %seq_mdl_H     = ();   # key is seq name, value is best model for this sequence
my %seq_sc_H      = ();   # key is seq name, value is score for best model for this sequence
if($nmkey > 1) { 
  foreach $mkey (@mkey_A) {
    parse_sqc_clsonly_file($sqc_H{$mkey}, $mkey, \%seq_H, \@seq_A, \%seq_mkey_H, \%seq_mdl_H, \%seq_sc_H, \%opt_HH, $FH_HR);
  }
}

# Fill per-mkey lists of sequences
my %seqlist_HA = (); # key is mkey, value is array of sequences that match to this mkey
my $nmkey_used = 0;  # number of mkey (libraries) we have at least one sequence to rerun v-annotate.pl for
my %mkey_ct_H  = (); # key is mkey, value is number of seqs assigned to that mkey, 'undef' if 0
if($nmkey > 1) {
  foreach my $seqname (@seq_A) {
    if(defined $seq_mkey_H{$seqname}) {
      my $mkey = $seq_mkey_H{$seqname};
      my $mdl  = $seq_mdl_H{$seqname};
      if(! defined $seqlist_HA{$mkey}) {
        @{$seqlist_HA{$mkey}} = ();
        $mkey_ct_H{$mkey} = 0;
        $nmkey_used++;
      }
      push(@{$seqlist_HA{$mkey}}, $seqname);
      $mkey_ct_H{$mkey}++;
    }
  } 
  
  if((! $do_multi) && ($nmkey_used > 1)) {
    my $mkey_str = "";
    foreach $mkey (sort keys %seqlist_HA) {
      if($mkey_str ne "") { $mkey_str .= ", "; }
      $mkey_str .= $mkey;
    }
    ofile_FAIL("ERROR, -m not used but found matches to multiple libraries: $mkey_str", 1, $FH_HR);
  }
}
else {
  $nmkey_used = 1; # we didn't run in clsonly because we only have 1 library
}

###########################################################################
# Re-run v-annotate.pl for each model key that at least one seq matched to
###########################################################################
my @mkey_used_A = (); # array of the mkeys with at least one sequence 
my @mdl_file_A  = (); # array of mdl files to output before exiting
my @alc_file_A  = (); # array of alc files to output before exiting
my $mkey_fa_file = undef;
my $progress_str = undef;
if($nmkey_used > 0) { 
  foreach $mkey (@mkey_A) {
    if((defined $seqlist_HA{$mkey}) || ($nmkey == 1)) { # if $nmkey == 1, we didn't run --clsonly mode
      if($nmkey_used == 1) { 
        $mkey_fa_file = $in_fa_file;
        $progress_str = "Annotating all sequences with $mkey model library ... ";
      }
      else {
        $mkey_fa_file = $dir_tail . "/" . $mkey . ".fa";
        $in_sqfile->fetch_seqs_given_names(\@{$seqlist_HA{$mkey}}, 60, $mkey_fa_file);
        $progress_str = sprintf("Annotating %s sequences (%d) ... ", $mkey, scalar(@{$seqlist_HA{$mkey}}));
      }
      my $out_dir = $dir_tail . "/" . $mkey;
      push(@mkey_used_A, $mkey);
      push(@mdl_file_A, $dir_tail . "/" . $mkey . "/" . $mkey . ".vadr.mdl");
      push(@alc_file_A, $dir_tail . "/" . $mkey . "/" . $mkey . ".vadr.alc");

      $cmd = $execs_H{"v-annotate.pl"} . " --mkey $mkey --mdir $mkey_mdir_H{$mkey} $mkey_opts_H{$mkey} $mkey_fa_file $out_dir";
      if(! $do_verbose) { $cmd .= " > /dev/null"; }
      my $start_secs = ofile_OutputProgressPrior($progress_str, $progress_w, $FH_HR->{"log"}, *STDOUT);
      utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
      ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
    }
  }
}

$start_secs = ofile_OutputProgressPrior("Generating tabular output", $progress_w, $log_FH, *STDOUT);

# create the @data_lib_AA
my $mkey_idx = 1;
my $mdl_idx = 1;

# open files for writing
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "lib", $out_root . ".lib", 1, 1, "per-model library tabular summary file");
my @head_lib_AA = ();
my @data_lib_AA = ();
my @clj_lib_A   = ();
@{$head_lib_AA[0]} = ("",    "",        "num");
@{$head_lib_AA[1]} = ("idx", "library", "seqs");
@clj_lib_A         = (1,     1,         0);

foreach $mkey (@mkey_A) {
  my $nseq2print = (defined $mkey_ct_H{$mkey}) ? $mkey_ct_H{$mkey} : 0;
  if(scalar(@mkey_A) == 1) { # we didn't run clsonly mode, set nseq to '-'
    $nseq2print = "-";
  }
  push(@data_lib_AA, [$mkey_idx, $mkey, $nseq2print]);
  $mkey_idx++;
}

ofile_TableHumanOutput(\@data_lib_AA, \@head_lib_AA, \@clj_lib_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"lib"}, undef, $FH_HR);
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

############
# Conclude #
############
output_lib_mdl_and_alc_files_and_remove_temp_files($in_nseq, $peek_nseq, \@mkey_used_A, \@mdl_file_A, \@alc_file_A, \@to_remove_A, \%opt_HH, \%ofile_info_HH);

if($nmkey_used == 0) { # matches were found to zero libraries
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
#  $config_file:  path to config file
#  $mkey_AR:      REF to array of model library keys, one per model library hashes with
#                 information on the features, PRE-FILLED
#  $mkey_mdir_HR: REF to array of hashes with information on the features, PRE-FILLED
#  $mkey_opts_HR: REF to array of hashes with information on the features, PRE-FILLED
#  $opt_HHR:      REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $FH_HR:        REF to hash of file handles
#
# Returns:  void
#           
# Dies:     if problem parsing config file
#
#################################################################
sub parse_config_file { 
  my $sub_name = "parse_config_file"; 
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($config_file, $mkey_AR, $mkey_mdir_HR, $mkey_opts_HR, $opt_HHR, $FH_HR) = (@_);

  open(CONFIG, $config_file) || ofile_FileOpenFailure($config_file, $sub_name, $!, "reading", $FH_HR);

  my $line;
  my ($mkey, $mdir, $opts);
  while($line = <CONFIG>) {
    chomp $line;
    if(($line =~ m/\w/) && ($line !~ m/^\#/)) {
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) < 2) {
        ofile_FAIL("ERROR all non-comment lines should include at least two white space delimited fields: <modemkey> <modeldir>\nread line:\n$line", 1, $FH_HR);
      }
      my ($mkey, $mdir) = ($el_A[0], $el_A[1]);
      my $opts = "";
      for(my $i = 2; $i < scalar(@el_A); $i++) {
        if($opts ne "") { $opts .= " "; }
        $opts .= $el_A[$i];
      }
      if(defined $mkey_mdir_H{$mkey}) {
        ofile_FAIL("ERROR read model key $mkey twice in config file", 1, $FH_HR);
      }
      push(@{$mkey_AR}, $mkey);
      $mkey_mdir_HR->{$mkey} = $mdir;
      $mkey_opts_HR->{$mkey} = $opts;
    }
  }
  close(CONFIG);

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
#  $mkey:          REF model key (e.g. flu) that this sqc file pertains to
#  $seq_HR:        REF to hash of sequence names, key is seq name, value is 1, to fill here
#  $seq_AR:        REF to array of sequence names, to fill here
#  $seq_mkey_HR:   REF to hash, key is sequence name, value is winning mkey, to fill here
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
  my $nargs_exp = 9;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqc_file, $mkey, $seq_HR, $seq_AR, $seq_mkey_HR, $seq_mdl_HR, $seq_sc_HR, $opt_HHR, $FH_HR) = (@_);

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
        ofile_FAIL("ERROR problem parsing sqc file $sqc_H{$mkey}", 1, $FH_HR);
      }
      my ($seqname, $pf, $mdl, $score) = ($el_A[1], $el_A[3], $el_A[5], $el_A[8]);

      if(! defined $seq_HR->{$seqname}) {
        push(@seq_A, $seqname);
        $seq_HR->{$seqname} = 1;
      }
      if($el_A[3] eq "PASS") {
        my $keep_flag = 1;
        if(defined $seq_mkey_HR->{$seqname}) {
          # this sequence already matched a model for a different $mkey
          # we either:
          # 1) die with error message
          # 2) figure out best model for this sequence
          #    either first mkey seen, or top scoring mkey
          if($do_lone) {
            ofile_FAIL("ERROR sequence $seqname matched to two libraries: $seq_mkey_HR->{$seqname} and $mkey, omit --lone to allow this", 1, $FH_HR);
          }
          if($do_first) {
            $keep_flag = 0; # keep existing value in $seq_mkey_HR->{$seqname}
          }
          else {
            $keep_flag = ($score > $seq_sc_HR->{$seqname}) ? 1 : 0;
          }
        }
        if($keep_flag) {
          $seq_mkey_HR->{$seqname} = $mkey;
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
#  $peek_nseq:      number of sequences peeked at, only relevant if -p
#  $mkey_used_AR:   ref to array of model keys we want to output .mdl and .alc files for
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

  my ($in_nseq, $peek_nseq, $mkey_used_AR, $mdl_file_AR, $alc_file_AR, $to_remove_AR, $opt_HHR, $ofile_info_HHR) = (@_);

  # close the two files we may output to stdout and the log
  close($ofile_info_HHR->{"FH"}{"lib"});
  
  my $FH_HR   = $ofile_info_HH{"FH"};
  my $do_peek = opt_Get("-p", $opt_HHR);
  my $sum_str = "";

  if(($do_peek) && ($peek_nseq < $in_nseq)) {
    $sum_str = sprintf("# Summary of seqs matching each library (only %d of %d seqs scanned due to -p):", $peek_nseq, $in_nseq);
  }
  else {
    $sum_str = "# Summary of sequences matching each library:";
  }
  
  my $nmkey = scalar(@{$mkey_used_AR});
  if($nmkey != scalar(@{$mdl_file_AR})) {
    ofile_FAIL("ERROR, in $sub_name, unexpected number of mdl files", 1, $FH_HR);
  }
  if(scalar(@{$mdl_file_AR}) != (scalar(@{$alc_file_AR}))) {
    ofile_FAIL("ERROR, in $sub_name, number of mdl and alc files differ", 1, $FH_HR);
  }

  my @conclude_A = ();
  my @file_A = ();
  my ($mkey, $mdl_file, $alc_file) = (undef, undef, undef);
  push(@conclude_A, "#");
  push(@conclude_A, $sum_str);
  push(@conclude_A, "#");
  utl_FileLinesToArray($ofile_info_HHR->{"fullpath"}{"lib"}, 1, \@file_A, $FH_HR);
  push(@conclude_A, @file_A);
  push(@conclude_A, "#");

  for(my $m = 0; $m < $nmkey; $m++) {
    $mkey     = $mkey_used_AR->[$m];
    $mdl_file = $mdl_file_AR->[$m];
    $alc_file = $alc_file_AR->[$m];

    push(@conclude_A, "#");
    push(@conclude_A, "# Summary of sequences matching the $mkey model library:");
    push(@conclude_A, "#");

    @file_A = ();
    utl_FileLinesToArray($mdl_file, 1, \@file_A, $FH_HR);
    push(@conclude_A, @file_A);
    push(@conclude_A, "#");

    @file_A = ();
    utl_FileLinesToArray($alc_file, 1, \@file_A, $FH_HR);
    if(scalar(@file_A == 3)) {
      push(@conclude_A, "# Zero alerts reported for seqs matching to the $mkey library.");
    }
    else {
      push(@conclude_A, "# Summary of reported alerts for seqs matching the $mkey library:");
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
#             and updating the @{$mkey_AR}, %{$mkey_mdir_HR}
#             and %{$mkey_opts_HR} data structures.
#
# Arguments: 
#  $mkey_AR:      REF to array of all mkeys read from config file, modified here
#  $mkey_mdir_HR: REF to hash of directories for each model key, modified here
#  $mkey_opts_HR: REF to hash of options for each model key, modified here
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
  
  my ($mkey_AR, $mkey_mdir_HR, $mkey_opts_HR, $opt_HHR) = @_;
  
  my @only_A = ();  # array of values in --only arg
  my @skip_A = ();  # array of values in --skip arg
  my %only_H = ();  # 'exists' hash for values in --only arg
  my %skip_H = ();  # 'exists' hash for values in --skip arg
  my @new_mkey_A = ();      # new array of mkeys we will replace @{$mkey_AR} with before returning
  my %new_mkey_mdir_H = (); # new hash  we will replace %{$mkey_mdir_HR} with before returning
  my %new_mkey_opts_H = (); # new hash  we will replace %{$mkey_opts_HR} with before returning

  my $die_str = "";
  if(opt_IsUsed("--only", $opt_HHR)) { 
    @only_A = split(",", opt_Get("--only", $opt_HHR));
    foreach my $only_mkey (@only_A) {
      $only_H{$only_mkey} = 1;
      if(! defined $mkey_mdir_HR->{$only_mkey}) {
        $die_str .= "\t$only_mkey specified in --only option but not listed in config file\n";
      }
    }
    foreach my $mkey (@{$mkey_AR}) {
      if(defined $only_H{$mkey}) {
        push(@new_mkey_A, $mkey);
        $new_mkey_mdir_H{$mkey} = $mkey_mdir_H{$mkey};
        $new_mkey_opts_H{$mkey} = $mkey_opts_H{$mkey};
      }
    }
  }
  if(opt_IsUsed("--skip", $opt_HHR)) { 
    if(scalar(@only_A) != 0) {
      # this should have been enforced by opt_ValidateSet() 
      ofile_FAIL("ERROR, in $sub_name, --only and --skip both used, pick one", 1, $FH_HR);
    }
    @skip_A = split(",", opt_Get("--skip", $opt_HHR));
    foreach my $skip_mkey (@skip_A) {
      $skip_H{$skip_mkey} = 1;
      if(! defined $mkey_mdir_HR->{$skip_mkey}) {
        $die_str .= "\t$skip_mkey specified in --skip option but not listed in config file\n";
      }
    }
    foreach my $mkey (@{$mkey_AR}) {
      if(! defined $skip_H{$mkey}) {
        push(@new_mkey_A, $mkey);
        $new_mkey_mdir_H{$mkey} = $mkey_mdir_H{$mkey};
        $new_mkey_opts_H{$mkey} = $mkey_opts_H{$mkey};
      }
    }
  }

  if($die_str ne "") {
      ofile_FAIL("ERROR, in $sub_name:\n$die_str\n", 1, $FH_HR);
  }
  
  @{$mkey_AR} = ();
  %{$mkey_mdir_HR} = ();
  %{$mkey_opts_HR} = ();

  # copy values
  @{$mkey_AR} = @new_mkey_A;
  foreach $mkey (@{$mkey_AR}) {
    $mkey_mdir_HR->{$mkey} = $new_mkey_mdir_H{$mkey};
    $mkey_opts_HR->{$mkey} = $new_mkey_opts_H{$mkey};
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
#  $mkey_AR:      REF to array of all mkeys read from config file, modified here
#  $mkey_mdir_HR: REF to hash of directories for each model key, modified here
#  $mkey_opts_HR: REF to hash of options for each model key, modified here
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
  
  my ($config_file, $mkey_AR, $mkey_mdir_HR, $mkey_opts_HR, , $pkgname, $version, $releasedate, $opt_HHR) = @_;

  my $div_line = utl_StringMonoChar(60, "#", undef) . "\n";

  print $div_line;
  print "#\n";
  print "# $pkgname $version ($releasedate)\n";
  print "#\n";
  print "# config file: $config_file\n";
  print "#\n";

  my $do_lib = opt_IsUsed("--l_lib", $opt_HHR) ? 1 : 0;
  my $out_lib = undef;
  my $mkey = undef;

  # if --l_lib doesn't exist, exit
  if($do_lib)  {
    $out_lib = opt_Get("--l_lib", $opt_HHR);
    if(! defined $mkey_mdir_HR->{$out_lib}) {
      my $die_str = "ERROR, model library $out_lib specified with --l_lib does not exist in config file.\nExisting libraries are:\n";
      foreach $mkey (@mkey_A) {
        $die_str .= "\t$mkey\n";
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
    foreach my $mkey (@{$mkey_AR}) {
      if((! $do_lib) || ($mkey eq $out_lib)) { 
        push(@data_AA, [$mkey, $mkey_mdir_HR->{$mkey}]);
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
    foreach my $mkey (@{$mkey_AR}) {
      if((! $do_lib) || ($mkey eq $out_lib)) { 
        push(@data_AA, [$mkey, $mkey_opts_HR->{$mkey}]);
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
    for(my $k = 0; $k < scalar(@{$mkey_AR}); $k++) {
      my $mkey = $mkey_A[$k];
      $minfo_file = $mkey_mdir_HR->{$mkey} . "/" . $mkey . ".minfo";
      @mdl_info_AH = ();
      %ftr_info_HAH = ();
      utl_FileValidateExistsAndNonEmpty($minfo_file, "$mkey model info file", undef, 1, undef);
      vdr_ModelInfoFileParse($minfo_file, \@reqd_mdl_keys_A, \@reqd_ftr_keys_A, \@mdl_info_AH, \%ftr_info_HAH, undef);
      my $nmdl = scalar(@mdl_info_AH);
      if((! $do_lib) && ($k > 0)) { push(@data_AA, []); } # blank line
      for(my $m = 0; $m < $nmdl; $m++) {
        if((! $do_lib) || ($mkey eq $out_lib)) { 
          push(@data_AA,
               [(sprintf("%d.%d", ($do_lib ? 1 : ($k+1)), ($m+1))), 
                $mkey,
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
