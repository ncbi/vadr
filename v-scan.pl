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
my $env_vadr_scripts_dir  = utl_DirEnvVarValid("VADRSCRIPTSDIR");
my $env_vadr_easel_dir    = utl_DirEnvVarValid("VADREASELDIR");

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
opt_Add("-v",           "boolean", 0,         $g,    undef, undef,      "be verbose",                                                                              "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--lone",       "boolean", 0,         $g,    undef, undef,      "exit if at least one sequence matches to multiple libraries",                             "exit if at least one sequence matches to multiple libraries", \%opt_HH, \@opt_order_A);
opt_Add("--first",      "boolean", 0,         $g,    undef,"--lone",    "if a seq matches > 1 model library use first one [df: use best scoring]",                 "if a seq matches > 1 model library use first one [df: use best scoring]", \%opt_HH, \@opt_order_A);
opt_Add("--origfa",     "boolean", 0,         $g,    undef,   undef,    "do not copy fasta file prior to analysis, use original",                 "do not copy fasta file prior to analysis, use original", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,         $g,    undef, undef,      "leaving intermediate files on disk",                                                      "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay = 
    &GetOptions('h'        => \$GetOptions_H{"-h"}, 
# basic options
                'f'        => \$GetOptions_H{"-f"},
                'v'        => \$GetOptions_H{"-v"},
                'lone'     => \$GetOptions_H{"--lone"}, 
                'first'    => \$GetOptions_H{"--first"},
                'origfa'   => \$GetOptions_H{"--origfa"},
                'keep'     => \$GetOptions_H{"--keep"});

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

my $do_verbose = opt_Get("-v",     \%opt_HH);
my $do_keep    = opt_Get("--keep", \%opt_HH);

# check that number of command line args is correct
if(scalar(@ARGV) != 3) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do $executable -h\n\n";
  exit(1);
}

my ($orig_in_fa_file, $dir, $in_config_file) = (@ARGV);

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
my @arg_desc_A = ("sequence file", "output directory", "config file");
my @arg_A      = ($orig_in_fa_file, $dir, $in_config_file);
my %extra_H    = ();
$extra_H{"\$VADRSCRIPTSDIR"}  = $env_vadr_scripts_dir;
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

my $progress_w = 60; # the width of the left hand column in our progress output, hard-coded
my $start_secs = ofile_OutputProgressPrior("Validating input", $progress_w, $log_FH, *STDOUT);

my @to_remove_A   = (); # list of files to remove at end of subroutine, if --keep not used

###########################################
# Validate that we have all the files we need:
# fasta file
utl_FileValidateExistsAndNonEmpty($orig_in_fa_file, "input fasta sequence file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
utl_FileValidateExistsAndNonEmpty($in_config_file,  "input config file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty

my @mkey_A = ();      # array of model library keys, read from config file
my %mkey_mdir_H = (); # hash of model directories for each model library key, read from config file, key is model key
my %mkey_opts_H = (); # hash of options for each model library key, read from config file, key is model key
parse_config_file($in_config_file, \@mkey_A, \%mkey_mdir_H, \%mkey_opts_H, \%opt_HH, $FH_HR);

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

#################################################
# foreach model key, run v-annotate.pl --cls_only
my %out_dir_H = (); # hash of output directories
my %sqc_H = ();     # hash of sqc files
my $mkey;
foreach $mkey (@mkey_A) {
  $out_dir_H{$mkey} = $dir_tail . "/" . $mkey . ".0";
  $sqc_H{$mkey} = $out_dir_H{$mkey} . "/" . $mkey . ".0.vadr.sqc";
  $cmd = $execs_H{"v-annotate.pl"} . " -f -s --origfa --cls_only --mkey $mkey --mdir $mkey_mdir_H{$mkey} $in_fa_file $out_dir_H{$mkey}";
  if(! $do_verbose) { $cmd .= " > /dev/null"; }
  my $start_secs = ofile_OutputProgressPrior(sprintf("Scanning for %s sequences ... ", $mkey), $progress_w, $log_FH, *STDOUT);
  utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

# Parse sqc files to determine which seqs match best to each library, need
# to look at all sqc files before assigning sequences to a mkey
# because we may be determining best library based on score
my %seq_H         = ();   # 'exists' hash, key is sequence name, value is always '1' 
my @seq_A         = ();   # array of sequence names
my %seq_mkey_H    = ();   # key is seq name, value is best mkey for this sequence
my %seq_mdl_H     = ();   # key is seq name, value is best model for this sequence
my %seq_grp_H     = ();   # key is seq name, value is group of best model for this sequence
my %seq_subgrp_H  = ();   # key is seq name, value is subgroup of best model for this sequence
my %seq_sc_H      = ();   # key is seq name, value is score for best model for this sequence
foreach $mkey (@mkey_A) {
  parse_sqc_clsonly_file($sqc_H{$mkey}, $mkey, \%seq_H, \@seq_A, \%seq_mkey_H, \%seq_mdl_H, \%seq_grp_H, \%seq_subgrp_H, \%seq_sc_H, \%opt_HH, $FH_HR);
}

# Fill per-mkey lists of sequences
my %seqlist_HA = (); # key is mkey, value is array of sequences that match to this mkey
my $nmkey = 0;       # number of mkey (libraries) we have at least one sequence to rerun v-annotate.pl for

my %mkey_mdl_ct_HH = ();
my %mkey_mdl_order_HA = ();
my %grp_H    = (); # key is mdl name, value is group
my %subgrp_H = (); # key is mdl name, value is subgroup
foreach my $seqname (@seq_A) {
  if(defined $seq_mkey_H{$seqname}) {
    my $mkey   = $seq_mkey_H{$seqname};
    my $mdl    = $seq_mdl_H{$seqname};
    my $grp    = $seq_grp_H{$seqname};
    my $subgrp = $seq_subgrp_H{$seqname};
    if(! defined $seqlist_HA{$mkey}) {
      @{$seqlist_HA{$mkey}} = ();
      %{$mkey_mdl_ct_HH{$mkey}} = ();
      @{$mkey_mdl_order_HA{$mkey}} = ();
      $nmkey++;
    }
    if(! defined $mkey_mdl_ct_HH{$mkey}{$mdl}) {
      $mkey_mdl_ct_HH{$mkey}{$mdl} = 0;
      push(@{$mkey_mdl_order_HA{$mkey}}, $mdl);
      $grp_H{$mdl}    = $grp;
      $subgrp_H{$mdl} = $subgrp;
    }
    push(@{$seqlist_HA{$mkey}}, $seqname);
    $mkey_mdl_ct_HH{$mkey}{$mdl}++;
  }
} 

if($nmkey > 0) { 
  foreach $mkey (@mkey_A) {
    if(defined $seqlist_HA{$mkey}) {
      my $mkey_fasta_file = $dir_tail . "/" . $mkey . ".fa";
      my $out_dir = $dir_tail . "/" . $mkey;
      $in_sqfile->fetch_seqs_given_names(\@{$seqlist_HA{$mkey}}, 60, $mkey_fasta_file);
      $cmd = $execs_H{"v-annotate.pl"} . " --mkey $mkey --mdir $mkey_mdir_H{$mkey} $mkey_opts_H{$mkey} $mkey_fasta_file $out_dir";
      if(! $do_verbose) { $cmd .= " > /dev/null"; }
      my $start_secs = ofile_OutputProgressPrior(sprintf("Annotating %s sequences (%d) ... ", $mkey, scalar(@{$seqlist_HA{$mkey}})), $progress_w, $FH_HR->{"log"}, *STDOUT);
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
@{$head_lib_AA[0]} = ("",    "",        "",      "",      "",         "num");
@{$head_lib_AA[1]} = ("idx", "library", "model", "group", "subgroup", "seqs");
@clj_lib_A         = (1,     1,         1,       1,       1,          0);

foreach $mkey (@mkey_A) {
  my $mdl_idx = 1;
  my $mkey_mdl_idx = sprintf("%d.%d", $mkey_idx, $mdl_idx);
  if(! defined $mkey_mdl_ct_HH{$mkey}) {
    push(@data_lib_AA, [$mkey_mdl_idx, $mkey, "-", "-", "-", 0 ]);
  }
  else {
    foreach my $mdl (@{$mkey_mdl_order_HA{$mkey}}) {
      push(@data_lib_AA, [$mkey_mdl_idx, $mkey, $mdl, $grp_H{$mdl}, $subgrp_H{$mdl}, $mkey_mdl_ct_HH{$mkey}{$mdl} ]);
      $mdl_idx++;
    }
  }
  $mkey_idx++;
}

ofile_TableHumanOutput(\@data_lib_AA, \@head_lib_AA, \@clj_lib_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"lib"}, undef, $FH_HR);
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

############
# Conclude #
############
output_lib_file_and_remove_temp_files($zero_mdl, \@to_remove_A, \%opt_HH, \%ofile_info_HH);
                                      
# remove unwanted files, unless --keep
if(! opt_Get("--keep", \%opt_HH)) { 
  my @to_actually_remove_A = (); # sanity check: make sure the files we're about to remove actually exist
  my %to_actually_remove_H = (); # sanity check: to make sure we don't try to delete 
  foreach my $to_remove_file (@to_remove_A) { 
    if((defined $to_remove_file) && (-e $to_remove_file) && (! defined $to_actually_remove_H{$to_remove_file})) { 
      push(@to_actually_remove_A, $to_remove_file); 
      $to_actually_remove_H{$to_remove_file} = 1; 
    }
  }
  utl_FileRemoveList(\@to_actually_remove_A, "v-scan.pl", \%opt_HH, $FH_HR);
}

# output lib file to stdout
my @file_A = ();
my $line;
utl_FileLinesToArray($ofile_info_HH{"fullpath"}{"lib"}, 0, \@file_A, $FH_HR);
foreach $line (@file_A) {
  print $line . "\n";
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
#  $seq_grp_HR:    REF to hash, key is sequence name, value is group of winning model, "-" if undef, to fill here
#  $seq_subgrp_HR: REF to hash, key is sequence name, value is subgroup of winning model, "-" if undef, to fill here
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
  my $nargs_exp = 11;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqc_file, $mkey, $seq_HR, $seq_AR, $seq_mkey_HR, $seq_mdl_HR, $seq_grp_HR, $seq_subgrp_HR, $seq_sc_HR, $opt_HHR, $FH_HR) = (@_);

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
      my ($seqname, $pf, $mdl, $grp, $subgrp, $score) = ($el_A[1], $el_A[3], $el_A[5], $el_A[6], $el_A[7], $el_A[8]);

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
          $seq_mkey_HR->{$seqname}   = $mkey;
          $seq_mdl_HR->{$seqname}    = $mdl;
          $seq_grp_HR->{$seqname}    = $grp;
          $seq_subgrp_HR->{$seqname} = $subgrp;
          $seq_sc_HR->{$seqname}     = $score;
        }
      }
    }
  }
  close(SQC);

  return;
}

#################################################################
# Subroutine: output_lib_and_remove_temp_files()
# Incept:     EPN, Mon Feb 10 14:13:25 2025
#             based on v-annotate.pl:output_mdl_and_alc_files_and_remove_temp_files()
# Purpose:    Output the lib file and remove all files
#             if (@{$to_remove_A}) unless --keep. 
#
# Arguments:
#  $to_remove_AR:   ref to array of files to remove
#  $opt_HHR:        ref to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR: ref to 2D hash of output file information, added to here
#             
# Returns:  void
#
#################################################################
sub output_lib_file_and_remove_temp_files { 
  my $sub_name = "output_lib_file_and_remove_temp_files";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($to_remove_AR, $opt_HHR, $ofile_info_HHR) = (@_);

  # close the two files we may output to stdout and the log
  close($ofile_info_HHR->{"FH"}{"lib"});
  
  my $FH_HR  = $ofile_info_HH{"FH"};
  
  my @conclude_A = ();
  push(@conclude_A, "#");
  push(@conclude_A, "# Summary of classified sequences:");
  push(@conclude_A, "#");
  my @file_A = ();
  utl_FileLinesToArray($ofile_info_HHR->{"fullpath"}{"mdl"}, 1, \@file_A, $FH_HR);
  push(@conclude_A, @file_A);
  push(@conclude_A, "#");
  if($do_clsonly) {
    push(@conclude_A, "# Only classification-related alerts detected due to --cls_only.");
  }
  elsif($zero_alt) { 
    push(@conclude_A, "# Zero alerts were reported.");
  }
  else { 
    push(@conclude_A, "# Summary of reported alerts:");
    push(@conclude_A, "#");
    my @file_A = ();
    utl_FileLinesToArray($ofile_info_HHR->{"fullpath"}{"alc"}, 1, \@file_A, $FH_HR);
    push(@conclude_A, @file_A);
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
