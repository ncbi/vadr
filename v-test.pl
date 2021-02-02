#!/usr/bin/env perl
#
use strict;
use warnings;
use Getopt::Long qw(:config no_auto_abbrev);
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

require "vadr.pm";
require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_utils.pm";

#########################################################
# Command line and option processing using sqp_opts.pm
#
# opt_HH: 2D hash:
#         1D key: option name (e.g. "-h")
#         2D key: string denoting type of information 
#                 (one of "type", "default", "group", "requires", "incompatible", "preamble", "help")
#         value:  string explaining 2D key:
#                 "type":         "boolean", "string", "int" or "real"
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

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
$opt_group_desc_H{"1"} = "basic options";
#     option            type       default               group   requires incompat    preamble-output                       help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                "display this help",                                  \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",        "force; if dir <output directory> exists, overwrite it",   \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                         "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("-s",           "boolean", 0,                        1,    undef, "--rmout",  "skip commands, they were already run, just compare files",  "skip commands, they were already run, just compare files",   \%opt_HH, \@opt_order_A);
opt_Add("-m",           "boolean", 0,                        1,"--noteamcity", undef, "benchmark mode: compare certain fields, do not diff", "benchmark mode: compare certain fields, do not diff files",   \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"2"} = "options for defining variables in testing files";
#       option       type        default                group  requires incompat          preamble-output                                              help-output    
#opt_Add("--dirbuild",   "string",  undef,                    2,   undef, undef,       "build directory, replaces !dirbuild! in test file with <s>", "build directory, replaces !dirbuild! in test file with <s>", \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"3"} = "other options";
opt_Add("--rmout",      "boolean", 0,                        3,    undef, "-s",       "if output files listed in testin file already exist, remove them", "if output files listed in testin file already exist, remove them", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                        3,    undef, undef,      "leaving intermediate files on disk", "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);
opt_Add("--noteamcity", "boolean", 0,                        3,    undef, undef,      "do not output teamcity test info",   "do not output teamcity test info", \%opt_HH, \@opt_order_A);
opt_Add("--skipmsg",    "boolean", 0,                        3,    undef, undef,      "do not compare errors and warnings", "do not compare errors and warning lines", \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"4"} = "other expert options";
#       option       type          default     group  requires incompat  preamble-output                                 help-output    
opt_Add("--execname",   "string",  undef,         4,    undef, undef,   "define executable name of this script as <s>", "define executable name of this script as <s>", \%opt_HH, \@opt_order_A);        

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"},
                'v'            => \$GetOptions_H{"-v"},
                'f'            => \$GetOptions_H{"-f"},
                's'            => \$GetOptions_H{"-s"},
                'm'            => \$GetOptions_H{"-m"},
#                'dirbuild=s'   => \$GetOptions_H{"--dirbuild"},
                'rmout'        => \$GetOptions_H{"--rmout"},
                'keep'         => \$GetOptions_H{"--keep"},
                'noteamcity'   => \$GetOptions_H{"--noteamcity"},
                'skipmsg'      => \$GetOptions_H{"--skipmsg"},
                'execname=s'   => \$GetOptions_H{"--execname"});

my $total_seconds = -1 * ofile_SecondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $execname_opt  = $GetOptions_H{"--execname"};
my $executable    = (defined $execname_opt) ? $execname_opt : "v-test.pl";
my $usage         = "Usage: $executable [-options] <input test file e.g. testfiles/testin.1> <output directory to create>\n";
my $synopsis      = "$executable :: test VADR scripts [TEST SCRIPT]";
my $date          = scalar localtime();
my $version       = "1.1.3dev3";
my $releasedate   = "Feb 2021";
my $pkgname       = "VADR";

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

# check that number of command line args is correct
if(scalar(@ARGV) != 2) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do $executable -h\n\n";
  exit(1);
}

my ($test_file, $dir) = (@ARGV);

#############################
# create the output directory
#############################
my $cmd;               # a command to run with utl_RunCommand()
my @early_cmd_A = ();  # array of commands we run before our log file is opened

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

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("test file", "output directory");
my @arg_A      = ($test_file, $dir);
ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, undef);
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
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "log",  $out_root . ".log",  1, 1, "Output printed to screen");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "cmd",  $out_root . ".cmd",  1, 1, "List of executed commands");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".list", 1, 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
my $FH_HR  = $ofile_info_HH{"FH"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# now we have the log file open, output the banner there too
ofile_OutputBanner($log_FH, $pkgname, $version, $releasedate, $synopsis, $date, undef);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

# are we running in 'benchmark' mode?
my $do_mark = opt_Get("-m", \%opt_HH);

# read in the test file
my $progress_w = 50; # the width of the left hand column in our progress output, hard-coded
my $start_secs = ofile_OutputProgressPrior("Parsing test file", $progress_w, $log_FH, *STDOUT);
my @cmd_A      = (); # array of the commands to run
my @desc_A     = (); # array of the descriptions for the commands
my @outfile_AA = (); # array of arrays of output files to compare for each command
my @expfile_AA = (); # array of arrays of expected files to compare output to for each command
my @cmpstr_AA  = (); # array of arrays of 'compare' strings, to compare output to for each command, used only in benchmark mode (-m)
my @rmdir_AA   = (); # array of directories to remove after each command is completed
my $do_teamcity = opt_Get("--noteamcity", \%opt_HH) ? 0 : 1;
my $ncmd = parse_test_or_mark_file($test_file, \@cmd_A, \@desc_A, \@outfile_AA, \@expfile_AA, \@cmpstr_AA, \@rmdir_AA, \%opt_HH, $ofile_info_HH{"FH"});
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

my $npass = 0;
my $nfail = 0;
my $w_cmpcat = 20;
my @data_pertest_AA = (); # pertest file data, only used if -m
my @data_perline_AA = (); # perline file data, only used if -m
my @data_perval_AA  = (); # perval file data, only used if -m

for(my $i = 1; $i <= $ncmd; $i++) { 
  my $cmd  = $cmd_A[($i-1)];
  my $desc = $desc_A[($i-1)];
  if($do_teamcity) { 
    ofile_OutputString($log_FH, 1, sprintf("##teamcity[testStarted name='$desc' captureStandardOutput='true']\n"));
  }
  my $outfile_AR = \@{$outfile_AA[($i-1)]};
  my $expfile_AR = ($do_mark) ? undef : \@{$expfile_AA[($i-1)]};
  my $cmpstr_AR  = ($do_mark) ? \@{$cmpstr_AA[($i-1)]} : undef;
  my $rmdir_AR   = \@{$rmdir_AA[($i-1)]};
  my $npass_i    = 0; # number of tests that passed for current cmd
  my $nfail_i    = 0; # number of tests that failed for current cmd
  if((opt_IsUsed("-s", \%opt_HH)) && (opt_Get("-s", \%opt_HH))) { 
    # -s used, we aren't running commands, just comparing files
    $start_secs = ofile_OutputProgressPrior(sprintf("Skipping command %2d [%20s]", $i, $desc_A[($i-1)]), $progress_w, $log_FH, *STDOUT);
  }
  else { 
    # -s not used, run command
    $start_secs = ofile_OutputProgressPrior(sprintf("Running command %2d [%20s]", $i, $desc_A[($i-1)]), $progress_w, $log_FH, *STDOUT);
    utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 1, $ofile_info_HH{"FH"}); # 1: do not fail if command fails
  }
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  my $nout = scalar(@{$outfile_AR});
  for(my $j = 0; $j < $nout; $j++) { 
    if($do_mark) { 
      my @outfields_A = ();
      my @cmpfields_A = ();
      my ($cmpcat, $cmpfile) = parse_cmpstr($cmpstr_AR->[$j], \@outfields_A, \@cmpfields_A, $ofile_info_HH{"FH"});
      my $desccat = $desc . ":" . $cmpcat;
      my ($nlines, $nidentical, $ndifferent) = compare_two_files($desccat, $i.".".($j+1), $outfile_AR->[$j], $cmpfile, \@outfields_A, \@cmpfields_A, \@data_perline_AA, \@data_perval_AA, $ofile_info_HH{"FH"});
      my $idx = $i . "." . ($j+1);
      push(@data_pertest_AA, [$idx, $cmpcat, $nlines, $nidentical, $ndifferent, 
                              sprintf("%7.5f", ($nidentical / $nlines)),
                              sprintf("%7.5f", ($ndifferent / $nlines))]);
    }
    else { # default -m not enabled
      my $diff_file = $out_root . "." . $i . "." . ($j+1) . ".diff";
      my $pass = diff_two_files($outfile_AR->[$j], $expfile_AR->[$j], $diff_file, \%opt_HH, $ofile_info_HH{"FH"});
      if($pass) { $npass++; $npass_i++; }
      else      { $nfail++; $nfail_i++; }
    }
  }

  if(($nfail_i == 0) && (! opt_Get("--keep", \%opt_HH))) { # only remove dir if no tests failed
    my $nrmdir = (defined $rmdir_AR) ? scalar(@{$rmdir_AR}) : 0;
    for(my $k = 0; $k < $nrmdir; $k++) { 
      ofile_OutputString($log_FH, 1, sprintf("#\t%-60s ... ", "removing directory $rmdir_AR->[$k]"));
      utl_RunCommand("rm -rf $rmdir_AR->[$k]", opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"}); 
      ofile_OutputString($log_FH, 1, "done\n");
    }
  }

  if($do_teamcity) { 
    if($nfail_i > 0) { 
      ofile_OutputString($log_FH, 1, sprintf("##teamcity[testFailed name='$desc' message='v-test.pl failure']\n"));
    }
    ofile_OutputString($log_FH, 1, sprintf("##teamcity[testFinished name='$desc']\n"));
  }
}

##########
# Conclude
##########

if(! $do_mark) { 
  # summarize number of files checked and passed
  my $overall_pass = ($nfail == 0) ? 1 : 0;
  ofile_OutputString($log_FH, 1, "#\n#\n");
  if($overall_pass) { 
    ofile_OutputString($log_FH, 1, "# PASS: all $npass files were created correctly.\n");
    $total_seconds += ofile_SecondsSinceEpoch();
    ofile_OutputConclusionAndCloseFilesOk($total_seconds, $dir, \%ofile_info_HH);
  }
  else { 
    ofile_OutputString($log_FH, 1, sprintf("# FAIL: %d of %d files were not created correctly.\n", $nfail, $npass+$nfail));
    $total_seconds += ofile_SecondsSinceEpoch();
    ofile_OutputConclusionAndCloseFilesFail($total_seconds, $dir, \%ofile_info_HH);
    ofile_FAIL("ERROR, at least one test FAILed", 1, undef);
  }
}
else { 
  # -m enabled, output summary files
  my @head_pertest_AA = ();
  @{$head_pertest_AA[0]} = ("test",     "",         "number",      "number",    "number",    "fraction",  "fraction");
  @{$head_pertest_AA[1]} = ("index",    "category", "comparisons", "identical", "different", "identical", "different");
  my @clj_pertest_A      = (1,          1,          0,             0,           0,           0,           0);
  # data already added to @data_pertest_AA above

  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "pertest", $out_root . ".pertest", 1, 1, "per-test tabular summary file");
  ofile_TableHumanOutput(\@data_pertest_AA, \@head_pertest_AA, \@clj_pertest_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"pertest"}, undef, $FH_HR);
  close($ofile_info_HH{"FH"}{"pertest"}); 

  my @head_perline_AA = ();
  @{$head_perline_AA[0]} = ("line",     "description+",    "anchor",      "expected",  "observed",  "identical/");
  @{$head_perline_AA[1]} = ("index",    "category",        "data",        "data",      "data",      "different");
  my @clj_perline_A      = (1,          1,                 1,             1,           1,           1);
  # data already added to @data_perline_AA above

  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "perline", $out_root . ".perline", 1, 1, "per-line tabular summary file");
  ofile_TableHumanOutput(\@data_perline_AA, \@head_perline_AA, \@clj_perline_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"perline"}, undef, $FH_HR);
  close($ofile_info_HH{"FH"}{"perline"}); 

  my @head_perval_AA = ();
  @{$head_perval_AA[0]} = ("value",    "description+",    "",            "number",      "number",    "number",    "fraction",  "fraction");
  @{$head_perval_AA[1]} = ("index",    "category",        "value",       "comparisons", "identical", "different", "identical", "different");
  my @clj_perval_A      = (1,          1,                 1,             0,             0,           0,           0,           0);
  # data already added to @data_perval_AA above

  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "perval", $out_root . ".perval", 1, 1, "per-value tabular summary file");
  ofile_TableHumanOutput(\@data_perval_AA, \@head_perval_AA, \@clj_perval_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"perval"}, undef, $FH_HR);
  close($ofile_info_HH{"FH"}{"perval"}); 

  my @conclude_A = ();
  push(@conclude_A, "#");
  push(@conclude_A, "# Summary of comparison tests:");
  push(@conclude_A, "#");
  my @file_A = ();
  utl_FileLinesToArray($ofile_info_HH{"fullpath"}{"pertest"}, 1, \@file_A, $FH_HR);
  push(@conclude_A, @file_A);
  push(@conclude_A, "#");
  foreach my $line (@conclude_A) { 
    ofile_OutputString($log_FH, 1, $line . "\n");
  }
}

exit 0;
#################################################################
# Subroutine:  parse_test_or_mark_file()
# Incept:      EPN, Wed May 16 16:03:40 2018
#
# Purpose:     Parse an input test or mark file and store the relevant
#              information in passed in array references.
#
# Arguments:
#   $testfile:    name of file to parse
#   $cmd_AR:      ref to array of commands to fill here
#   $desc_AR:     ref to array of descriptions to fill here
#   $outfile_AAR: ref to 2D array of output files to fill here
#   $expfile_AAR: ref to 2D array of expected files to fill here
#   $cmpstr_AAR:  ref to 2D array of compare strings to fill here, only if -m
#   $rmdir_AAR:   ref to 2D array of directories to remove after calling each command
#   $opt_HHR:     ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $FH_HR:       ref to hash of file handles, including "log" and "cmd"
#
# Returns:    number of commands read in $testfile
#
# Dies:       - if any of the expected files do not exist
#             - if number of expected files is not equal to number of output files
#               for any command
#             - if there are 0 output files for a given command
#             - if output file already exists
#             - if any desc values are duplicates
#################################################################
sub parse_test_or_mark_file { 
  my $sub_name = "parse_test_or_mark_file";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($testfile, $cmd_AR, $desc_AR, $outfile_AAR, $expfile_AAR, $cmpstr_AAR, $rmdir_AAR, $opt_HHR, $FH_HR) = @_;

  my $do_mark = opt_Get("-m", $opt_HHR);

  open(IN, $testfile) || ofile_FileOpenFailure($testfile, $sub_name, $!, "reading", $FH_HR);
  my $ncmd = 0;
  my $ndesc = 0;
  my $outfile;
  my $expfile;
  my $rmdir;
  my %desc_H = (); # used to prevent two desc from being identical

  while(my $line = <IN>) { 
    if(($line !~ m/^\#/) && ($line =~ m/\w/)) { 
      # example input file:
      # # comment line (ignored)
      #command: perl $VADRSCRIPTSDIR/v-annotate.pl -f $VADRSCRIPTSDIR/testfiles/dengue.r5.fa va-dengue.r5 > va-dengue.r5.out
      #desc: annotate-dengue-5-local
      #out: va-dengue.r5/va-dengue.r5.vadr.pass.ft 
      #exp: @VADRSCRIPTSDIR@/testfiles/expected-files/va-dengue.r5/va-dengue.r5.vadr.pass.ft 
      chomp $line;
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists

      if($line =~ s/^command\:\s+//) { 
        my $cmd = $line;
        if($ncmd > 0) { 
          # make sure we have read >= 1 outfiles and expfiles for previous command
          if(! (@{$outfile_AAR->[($ncmd-1)]})) { ofile_FAIL("ERROR did not read any out: lines for command " . ($ncmd+1), 1, $FH_HR); }
          my $nout_prv = scalar(@{$outfile_AAR->[($ncmd-1)]});
          if($do_mark) { 
            if(! (@{$cmpstr_AAR->[($ncmd-1)]})) { ofile_FAIL("ERROR did not read any cmp: lines for command " . ($ncmd+1), 1, $FH_HR); }
            my $ncmp_prv = scalar(@{$cmpstr_AAR->[($ncmd-1)]});
            if($nout_prv != $ncmp_prv) { 
              ofile_FAIL("ERROR different number of output and compare lines for command " . ($ncmd+1), 1, $FH_HR);
            }
          }
          else { 
            if(! (@{$expfile_AAR->[($ncmd-1)]})) { ofile_FAIL("ERROR did not read any exp: lines for command " . ($ncmd+1), 1, $FH_HR); }
            my $nexp_prv = scalar(@{$expfile_AAR->[($ncmd-1)]});
            if($nout_prv != $nexp_prv) { 
              ofile_FAIL("ERROR different number of output and expected lines for command " . ($ncmd+1), 1, $FH_HR);
            }
          }
        }
        # replace !<s>! with value of --<s> from options, die if it wasn't set or doesn't exist
        while($cmd =~ /\!(\w+)\!/) { 
          my $var = $1;
          my $varopt = "--" . $var;
          if(! opt_Exists($varopt, $opt_HHR)) { 
            ofile_FAIL("ERROR trying to replace !$var! in test file but option --$var does not exist in command line options", 1, $FH_HR); 
          }
          if(! opt_IsUsed($varopt, $opt_HHR)) { 
            ofile_FAIL("ERROR trying to replace !$var! in test file but option --$var was not specified on the command line, please rerun with --$var", 1, $FH_HR); 
          }
          my $replacevalue = opt_Get($varopt, $opt_HHR);
          $cmd =~ s/\!$var\!/$replacevalue/g;
        }
        push(@{$cmd_AR}, $cmd); 
        $ncmd++;
      }
      elsif($line =~ s/^desc\:\s+//) { 
        my $desc = $line;
        if(defined $desc_H{$desc}) { 
          ofile_FAIL("ERROR read desc $desc twice in $testfile", 1, $FH_HR);
        }          
        push(@{$desc_AR}, $desc); 
        $desc_H{$desc} = 1;
        $ndesc++;
      }
      elsif($line =~ s/^out\:\s+//) { 
        $outfile = $line;
        $outfile =~ s/^\s+//;
        $outfile =~ s/\s+$//;
        if($outfile =~ m/\s/) { ofile_FAIL("ERROR output file has spaces: $outfile", 1, $FH_HR); }
        if(scalar(@{$outfile_AAR}) < $ncmd) { 
          @{$outfile_AAR->[($ncmd-1)]} = ();
        }
        push(@{$outfile_AAR->[($ncmd-1)]}, $outfile);
        if((opt_IsUsed("-s", $opt_HHR)) && (opt_Get("-s", $opt_HHR))) { 
          # -s used, we aren't running commands, just comparing files, output files must already exist
          if(! -e $outfile) { ofile_FAIL("ERROR, output file $outfile does not already exist (and -s used)", 1, $FH_HR); }
        }
        elsif((opt_IsUsed("--rmout", $opt_HHR)) && (opt_Get("--rmout", $opt_HHR))) { 
          # --rmout used, remove any file that already exists
          if(-e $outfile) { utl_FileRemoveUsingSystemRm($outfile, "v-test.pl", $opt_HHR, $FH_HR); }
        }
        else { # -s not used, --rmout not used
          if(-e $outfile) { ofile_FAIL("ERROR, output file $outfile already exists (and -s/--rmout not used)", 1, $FH_HR); }
        }
      }
      elsif($line =~ s/^exp\:\s+//) { 
        if($do_mark) { ofile_FAIL("ERROR, read exp line but -m enabled", 1, $FH_HR); }
        $expfile = $line;
        $expfile =~ s/^\s+//;
        $expfile =~ s/\s+$//;
        # replace @<s>@ with value of $ENV{'<s>'}
        while($expfile =~ /\@(\w+)\@/) { 
          my $envvar = $1;
          my $replacevalue = $ENV{"$envvar"};
          $expfile =~ s/\@$envvar\@/$replacevalue/g;
        }
        if($expfile =~ m/\s/) { ofile_FAIL("ERROR expected file has spaces: $expfile", 1, $FH_HR) }
        if(scalar(@{$expfile_AAR}) < $ncmd) { 
          @{$expfile_AAR->[($ncmd-1)]} = ();
        }
        push(@{$expfile_AAR->[($ncmd-1)]}, $expfile);
        if(! -e $expfile) { ofile_FAIL("ERROR, expected file $expfile does not exist", 1, $FH_HR); }
      }
      elsif($line =~ s/^cmp\:\s+//) { 
        if(! $do_mark) { ofile_FAIL("ERROR, read cmp line but -m not enabled", 1, $FH_HR); }
        # expect 4 tokens:
        # cmp: <category>:<comma-separated-out-fields>:<comma-separated-cmp-fields>:<cmp-file-path>
        my $cmpstr = $line;
        $cmpstr =~ s/^\s+//;
        $cmpstr =~ s/\s+$//;
        # replace @<s>@ with value of $ENV{'<s>'}
        while($cmpstr =~ /\@(\w+)\@/) { 
          my $envvar = $1;
          my $replacevalue = $ENV{"$envvar"};
          $cmpstr =~ s/\@$envvar\@/$replacevalue/g;
        }

        if($cmpstr =~ m/\s/) { ofile_FAIL("ERROR compare string has spaces: $expfile", 1, $FH_HR) }
        if(scalar(@{$cmpstr_AAR}) < $ncmd) { 
          @{$cmpstr_AAR->[($ncmd-1)]} = ();
        }
        push(@{$cmpstr_AAR->[($ncmd-1)]}, $cmpstr);

        # validate the line, we don't store this info, we reparse it when needed
        my ($cmpcat, $cmpfile) = parse_cmpstr($line, undef, undef, $FH_HR);
      }
      elsif($line =~ s/^rmdir\:\s+//) { 
        $rmdir = $line;
        $rmdir =~ s/^\s+//;
        $rmdir =~ s/\s+$//;
        if(scalar(@{$rmdir_AAR}) < $ncmd) { 
          @{$rmdir_AAR->[($ncmd-1)]} = ();
        }
        push(@{$rmdir_AAR->[($ncmd-1)]}, $rmdir);
      }
      else { 
        ofile_FAIL("ERROR unable to parse line $line in $testfile", 1, $FH_HR);
      }
    }
  }
  close(IN);

  if($ndesc != $ncmd) { ofile_FAIL("ERROR did not read same number of descriptions and commands", 1, $FH_HR); }

  # for final command, check that number of exp and out files is equal
  if(! (@{$outfile_AAR->[($ncmd-1)]})) { ofile_FAIL("ERROR did not read any out: lines for command " . ($ncmd+1), 1, $FH_HR); }
  my $nout_prv = scalar(@{$outfile_AAR->[($ncmd-1)]});
  if($do_mark) { 
    if(! (@{$cmpstr_AAR->[($ncmd-1)]})) { ofile_FAIL("ERROR did not read any cmp: lines for command " . ($ncmd+1), 1, $FH_HR); }
    my $ncmp_prv = scalar(@{$cmpstr_AAR->[($ncmd-1)]});
    if($nout_prv != $ncmp_prv) { 
      ofile_FAIL("ERROR different number of output and compare lines for command " . ($ncmd+1), 1, $FH_HR);
    }
  }
  else { 
    if(! (@{$expfile_AAR->[($ncmd-1)]})) { ofile_FAIL("ERROR did not read any exp: lines for command " . ($ncmd+1), 1, $FH_HR); }
    my $nexp_prv = scalar(@{$expfile_AAR->[($ncmd-1)]});
    if($nout_prv != $nexp_prv) { 
      ofile_FAIL("ERROR different number of output and expected lines for command " . ($ncmd+1), 1, $FH_HR);
    }
  }

  return $ncmd;
}


#################################################################
# Subroutine:  diff_two_files()
# Incept:      EPN, Thu May 17 14:24:06 2018
#
# Purpose:     Diff two files, and output whether they are identical or not.
#
# Arguments:
#   $out_file:    name of output file
#   $exp_file:    name of expected file
#   $diff_file:   output file for diff command
#   $opt_HHR:     REF to 2D hash of option values, see top of sqp_opts.pm for description
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:    '1' if $outfile is identical to $expfile as determined by diff
#
# Dies:       If an expected file does not exist or is empty.
#
#################################################################
sub diff_two_files { 
  my $sub_name = "diff_two_files";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_file, $exp_file, $diff_file, $opt_HHR, $FH_HR) = @_;

  my $out_file_exists   = (-e $out_file) ? 1 : 0;
  my $exp_file_exists   = (-e $exp_file) ? 1 : 0;
  my $out_file_nonempty = (-s $out_file) ? 1 : 0;
  my $exp_file_nonempty = (-s $exp_file) ? 1 : 0;

  my $conclusion = "";
  my $pass = 0;

  if(! $exp_file_exists) { 
    ofile_FAIL("ERROR in $sub_name, expected file $exp_file does not exist", 1, $FH_HR) ;
  }
    
  ofile_OutputString($FH_HR->{"log"}, 1, sprintf("#\tchecking %-100s ... ", $out_file));

  if($out_file_exists) { 
    utl_RunCommand("diff -U 0 $out_file $exp_file > $diff_file", opt_Get("-v", $opt_HHR), 1, $FH_HR); # 1: don't die if command fails
    if(-s $diff_file) { 
      # copy the two files here:
      my $copy_of_out_file = $diff_file . ".out";
      my $copy_of_exp_file = $diff_file . ".exp";
      utl_RunCommand("cp $out_file $copy_of_out_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);
      utl_RunCommand("cp $exp_file $copy_of_exp_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);
      # analyze the diff file and print out how many lines 
      if($out_file =~ m/\.sqtable/ && $exp_file =~ m/\.sqtable/) { 
        if($out_file_nonempty && $exp_file_nonempty) { 
          my $sqtable_diff_file = $diff_file . ".man";
          compare_two_sqtable_files($out_file, $exp_file, $sqtable_diff_file, \%opt_HH, $FH_HR);
          $conclusion = "FAIL [files differ, see $sqtable_diff_file]";
        }
        else {
          if(! $exp_file_nonempty) { 
            $conclusion = "FAIL [files differ, expected file is empty, but output file is not (see $diff_file)]";
          }
          elsif(! $out_file_nonempty) { 
            $conclusion = "FAIL [files differ, output file is empty, but expected file is not (see $diff_file)]";
          }
          else { 
            ofile_FAIL("ERROR in $sub_name, unexpected case exp_file_nonempty: $exp_file_nonempty out_file_nonempty: $out_file_nonempty", 1, $FH_HR) ;
          }
        }
      }
      else { 
        $conclusion = "FAIL [files differ, see $diff_file]";
      }
    }
    else { 
      $conclusion = "pass";
      $pass = 1;
    }
  }
  else { 
    $conclusion = "FAIL [output file does not exist]";
  }

  ofile_OutputString($FH_HR->{"log"}, 1, "$conclusion\n");

  return $pass;
}

#################################################################
# Subroutine:  compare_two_files()
# Incept:      EPN, Tue Jul  7 11:20:56 2020
#
# Purpose:     Compares ceratin fields of two files
#
# Arguments:
#   $desccat:          description and category concatenated together, to push
#                      to @{$data_perline_AAR}.
#   $prefix:           index prefix
#   $out_file:         name of output file
#   $cmp_file:         name of compare file
#   $outfields_AR:     ref to array of fields from out file to compare
#   $cmpfields_AR:     ref to array of fields from cmp file to compare
#   $data_perline_AAR: ref to 2D array of perline data
#   $data_perval_AAR:  ref to 2D array of perval data
#   $FH_HR:            ref to hash of file handles, including "log" and "cmd"
#
# Returns:    '1' if $outfile is identical to $expfile as determined by diff
#
# Dies:       If an expected file does not exist or is empty.
#
#################################################################
sub compare_two_files { 
  my $sub_name = "compare_two_files";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($desccat, $prefix, $out_file, $cmp_file, $outfields_AR, $cmpfields_AR, $data_perline_AAR, $data_perval_AAR, $FH_HR) = @_;
  
  if(! -s $cmp_file) { 
    ofile_FAIL("ERROR in $sub_name, compare file $cmp_file does not exist or is empty", 1, $FH_HR) ;
  }
  if(! -s $out_file) { 
    ofile_FAIL("ERROR in $sub_name, output file $out_file does not exist or is empty", 1, $FH_HR) ;
  }
    
  # ofile_OutputString($FH_HR->{"log"}, 1, sprintf("#\tchecking %-100s ... ", $out_file));

  # read in fields of out file
  my @outdata0_A = (); # value  from field 0 from outfile
  my @outdataN_A = (); # values from fields 1->N from outfile, concatenated and separated by ';'
  my @cmpdata0_A = (); # value  from field 0 from cmpfile
  my @cmpdataN_A = (); # values from fields 1->N from cmpfile, concatenated and separated by ';'
  get_field_data_from_file($out_file, $outfields_AR, \@outdata0_A, \@outdataN_A);
  get_field_data_from_file($cmp_file, $cmpfields_AR, \@cmpdata0_A, \@cmpdataN_A);

  my $noutlines = scalar(@outdata0_A);
  my $ncmplines = scalar(@cmpdata0_A);
  if($noutlines ne $ncmplines) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, when trying to compare $out_file and $cmp_file, number of lines differ: out: %d != cmp: %d\n", scalar(@outdata0_A), scalar(@cmpdata0_A)), 1, $FH_HR);
  }

  my $nidentical = 0; 
  my $ndifferent = 0; 
  my $id_or_df = "";
  my %nidentical_perval_H = (); # key is a outdataN_A[$i] value, value is number of identities with $cmpdataN_A[$i]
  my %ndifferent_perval_H = (); # key is a outdataN_A[$i] value, value is number of differents with $cmpdataN_A[$i]
  my %ntotal_perval_H     = (); # key is a outdataN_A[$i] value, value is number of comparisons with $cmpdataN_A[$i]
  for (my $i = 0; $i < $noutlines; $i++) { 
    # first field should always be identical
    if($outdata0_A[$i] ne $cmpdata0_A[$i]) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, when trying to compare $out_file and $cmp_file, line %d, first token differs: out: %s != cmp: %s\n", $i+1, $outdata0_A[$i], $cmpdata0_A[$i]), 1, $FH_HR);
    }

    # remaining fields (which were concatenated together (after separating by ;) by get_field_data_from_file, may be different
    my $perval_key = $cmpdataN_A[$i];
    if($outdataN_A[$i] eq $cmpdataN_A[$i]) { 
      $nidentical++;
      if(! defined $ntotal_perval_H{$perval_key}) { 
        $nidentical_perval_H{$perval_key} = 1;
        $ndifferent_perval_H{$perval_key} = 0;
        $ntotal_perval_H{$perval_key} = 1;
      }
      else { 
        $nidentical_perval_H{$perval_key}++;
        $ntotal_perval_H{$perval_key}++;
      }
      $id_or_df = "identical";
    }
    else { 
      $ndifferent++;
      if(! defined $ndifferent_perval_H{$perval_key}) { 
        $ndifferent_perval_H{$perval_key} = 1;
        $nidentical_perval_H{$perval_key} = 0;
        $ntotal_perval_H{$perval_key} = 1;
      }
      else { 
        $ndifferent_perval_H{$perval_key}++;
        $ntotal_perval_H{$perval_key}++;
      }
      $id_or_df = "different";
    }
    push(@{$data_perline_AAR}, [sprintf("%s.%d", $prefix, ($i+1)), $desccat, $cmpdata0_A[$i], $cmpdataN_A[$i], $outdataN_A[$i], $id_or_df]);
  }
  push(@{$data_perline_AAR}, []); # empty array -> blank line

  # add perval data
  my $i = 1;
  foreach my $val (sort { $ntotal_perval_H{$b} <=> $ntotal_perval_H{$a} } keys %ntotal_perval_H) { 
    push(@{$data_perval_AAR}, [sprintf("%s.%d", $prefix, $i), $desccat, $val, $ntotal_perval_H{$val}, 
                               $nidentical_perval_H{$val}, $ndifferent_perval_H{$val}, 
                               sprintf("%7.5f", $nidentical_perval_H{$val} / $ntotal_perval_H{$val}), 
                               sprintf("%7.5f", $ndifferent_perval_H{$val} / $ntotal_perval_H{$val})]);
    $i++;
  }
  push(@{$data_perval_AAR}, []); # empty array -> blank line

  return ($noutlines, $nidentical, $ndifferent);
}

#################################################################
# Subroutine:  get_field_data_from_file()
# Incept:      EPN, Tue Jul  7 11:52:33 2020
#
# Purpose:     Fill arrays with data from specific fields from a file
#
# Arguments:
#   $in_file:    name of output file
#   $fields_AR:  ref to array with field indices
#   $data0_AR:   REF to array to fill with per-line data from $fields_AR[0]
#   $dataN_AR:   REF to array to fill with per-line data from all fields $fields_AR[$i] with $i>0, 
#
# Returns:    void, fills @{$data0_AR} and @{$dataN_AR}
#
# Dies:       If a non-# comment line of the file does not have enough fields
#
#################################################################
sub get_field_data_from_file {
  my $sub_name = "get_field_data_from_file";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($in_file, $fields_AR, $data0_AR, $dataN_AR) = @_;

  if(! -s $in_file) { 
    ofile_FAIL("ERROR in $sub_name, file $in_file does not exist or is empty", 1, $FH_HR) ;
  }

  open(IN, $in_file) || ofile_FileOpenFailure($in_file, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    if($line !~ m/^\#/) { 
      chomp $line;
      my @el_A = split(/\s+/, $line);
      my $nfields = scalar(@el_A);
      my $data0 = "";
      my $dataN = "";
      for(my $i = 0; $i < scalar(@{$fields_AR}); $i++) { 
        my $act_field = $fields_AR->[$i]-1;
        if($act_field >= $nfields) { 
          ofile_FAIL("ERROR in $sub_name, trying to get token $fields_AR->[$i] but only read $nfields on line $line of $in_file", 1, $FH_HR); 
        }
        if($i == 0) { 
          $data0 = $el_A[$act_field]; 
        }
        else { 
          if($dataN ne "") { $dataN .= ";"; }
          $dataN .= $el_A[$act_field];
        }
      }
      push(@{$data0_AR}, $data0);
      push(@{$dataN_AR}, $dataN);
    }
  }
  close(IN);

  return;
}
  
#################################################################
# Subroutine:  parse_cmpstr()
# Incept:      EPN, Tue Jul  7 09:57:23 2020
#
# Purpose:     Validate and parse a 'compare string'. 
#              Format should be 4 ':' delimited tokens:
#              <category>:<comma-separated-out-fields>:<comma-separated-cmp-fields>:<cmp-file-path>
#              Checks that <comma-separated-out-fields> and 
#              <comma-separated-cmp-fields> have the same number 
#              of fields. Validates that the file <cmp-file-path> 
#              exists and is non empty. 
#              Fill @{$outfields_AR} with field indices from 
#              <comma-separated-out-fields> [0..nfields-1]
#              Fill @{$cmpfields_AR} with field indices from 
#              <comma-separated-cmp-fields> [0..nfields-1]
#              
#
# Arguments:
#   $cmpstr:       name of output sqtable file
#   $outfields_AR: ref to array of out file field indices
#   $cmpfields_AR: ref to array of cmp file field indices
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:    Two values: 
#             <category>:      first ':'-delimited token from $cmpstr
#             <cmp-file-path>: path to cmpfile 
#
# Dies:       If cmpstr is not parseable or invalid
#             If <cmp-file-path> file does not exist or is empty.
#
#################################################################
sub parse_cmpstr { 
  my $sub_name = "parse_cmpstr";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmpstr, $outfields_AR, $cmpfields_AR, $FH_HR) = @_;

  my @el_A = split(":", $cmpstr);
  if(scalar(@el_A) != 4) { 
    ofile_FAIL(sprintf("ERROR unable to parse cmp: string $cmpstr, expected 4 ':' delimited tokens got: %d", scalar(@el_A)), 1, $FH_HR);
  }
  my ($category, $outfields_str, $cmpfields_str, $cmpfile) = (@el_A);

  my @outfields_A = split(",", $outfields_str); 
  my @cmpfields_A = split(",", $cmpfields_str); 
  my $nout = scalar(@outfields_A);
  my $ncmp = scalar(@cmpfields_A);
  if($nout == 1) { 
    ofile_FAIL("ERROR unable to parse cmp: string $cmpstr, read only 1 field in outfields str $outfields_str, must be at least 2", 1, $FH_HR);
  }
  if($ncmp == 1) { 
    ofile_FAIL("ERROR unable to parse cmp: string $cmpstr, read only 1 field in cmpfields str $cmpfields_str, must be at least 2", 1, $FH_HR);
  }
  if($nout != $ncmp) { 
    ofile_FAIL("ERROR unable to parse cmp: string $cmpstr, different numbers of comma-delimited tokens in outfields str $outfields_str and cmpfields str $cmpfields_str", 1, $FH_HR);
  }
  for(my $i = 0; $i < $nout; $i++) { 
    if($outfields_A[$i] !~ m/^\d+$/) { 
      ofile_FAIL("ERROR unable to parse cmp: string $cmpstr, found non-positive-integer in outfields str $outfields_str: $outfields_A[$i]", 1, $FH_HR);
    }
    if($cmpfields_A[$i] !~ m/^\d+$/) { 
      ofile_FAIL("ERROR unable to parse cmp: string $cmpstr, found non-positive-integer in cmpfields str $cmpfields_str: $cmpfields_A[$i]", 1, $FH_HR);
    }
  }
  if(defined $outfields_AR) { @{$outfields_AR} = @outfields_A; }
  if(defined $cmpfields_AR) { @{$cmpfields_AR} = @cmpfields_A; }

  if(! -e $cmpfile) { ofile_FAIL("ERROR, problem with cmp: string $cmpstr, expected file $cmpfile does not exist", 1, $FH_HR); }

  return ($category, $cmpfile); 
}

#################################################################
# Subroutine:  compare_two_sqtable_files()
# Incept:      EPN, Mon Jun 11 09:42:11 2018
#
# Purpose:     Compare two sqtable files outputting the number of 
#              lost, added, and changed features.
#
# Arguments:
#   $out_file:    name of output sqtable file
#   $exp_file:    name of expected sqtable file
#   $diff_file:   name of file to create with differences
#   $opt_HHR:     ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       If sequences are not in the same order in the two files
#
#################################################################
sub compare_two_sqtable_files { 
  my $sub_name = "compare_two_sqtable_files";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_file, $exp_file, $diff_file, $opt_HHR, $FH_HR) = @_;

  my $out_file_exists   = (-e $out_file) ? 1 : 0;
  my $exp_file_exists   = (-e $exp_file) ? 1 : 0;
  my $out_file_nonempty = (-s $out_file) ? 1 : 0;
  my $exp_file_nonempty = (-s $exp_file) ? 1 : 0;

  my $conclusion = "";
  my $pass = 0;

  my $skip_msg_lines = opt_Get("--skipmsg", $opt_HHR) ? 1 : 0;

  if(! $exp_file_exists) { 
    ofile_FAIL("ERROR in $sub_name, expected file $exp_file does not exist", 1, $FH_HR) ;
  }
  if(! $exp_file_nonempty) { 
    ofile_FAIL("ERROR in $sub_name, expected file $exp_file exists but is empty", 1, $FH_HR);
  }
  if(! $out_file_exists) { 
    ofile_FAIL("ERROR in $sub_name, output file $out_file does not exist", 1, $FH_HR) ;
  }
  if(! $out_file_nonempty) { 
    ofile_FAIL("ERROR in $sub_name, output file $out_file exists but is empty", 1, $FH_HR);
  }
    
  my @out_line_A = ();  # array of all lines in out file
  my @exp_line_A = ();  # array of all lines in exp file
  my $out_line;         # single line from out file
  my $exp_line;         # single line from exp file

  my %exp_seq_lidx_H = (); # key: sequence name, value line index where sequence annotation starts for key in $exp_file
  my @exp_seq_A      = (); # array of all sequence name keys in order in %exp_seq_lidx_H

  my %out_seq_lidx_H = (); # key: sequence name, value line index where sequence annotation starts for key in $out_file
  my @out_seq_A      = (); # array of all sequence name keys in order in %out_seq_lidx_H

  my %seq_exists_H = ();   # key is a sequence name, value is '1' if sequence exists in either @exp_seq_A or @out_seq_A

  my $lidx = 0; # line index
  my $seq;      # a sequence name         

  my $nseq = 0; # total number of sequences
  my $tot_nseq_ftr_identical = 0; # number of sequences with all features identically annotated between output and expected
  my $tot_nseq_ftr_diff      = 0; # number of sequences with >=1 feature differently annotated between output and expected
  my $tot_nseq_ftr_out_only  = 0; # number of sequences with >=1 feature only in the output file (not in the expected file)
  my $tot_nseq_ftr_exp_only  = 0; # number of sequences with >=1 feature only in the expected file (not in the output file)

  my $tot_nseq_note_identical = 0; # number of sequences with all notes identically annotated between output and expected
  my $tot_nseq_note_diff      = 0; # number of sequences with >=1 note differently annotated between output and expected
  my $tot_nseq_note_out_only  = 0; # number of sequences with >=1 note only in the output file (not in the expected file)
  my $tot_nseq_note_exp_only  = 0; # number of sequences with >=1 note only in the expected file (not in the output file)

  if($out_file_nonempty && $exp_file_nonempty) { 
    open(DIFFOUT, ">", $diff_file) || ofile_FileOpenFailure($diff_file, $sub_name, $!, "writing", $FH_HR);
    my %seq_H     = (); # key is sequence name, 1 is if sequence exists in @seq_A or not
    my @seq_A     = (); # array of sequence names in order
    my %seq_ftr_HH = (); # key1: sequence, key2: feature, value always 1
    my %seq_ftr_HA = (); # key1: sequence, array of features in both out and exp for this sequence
    my @file_A    = ($out_file, $exp_file);
    my %file_seq_ftr_HHH = (); # key 1: "out" or "exp", key 2: sequence name, key 3: a single feature for that sequence in output file, value: always 1
    my @filekey_A = ("out", "exp");
    my ($line, $file, $seq, $ftr, $filekey); 
    for(my $f = 0; $f < 2; $f++) { # loop over out and exp files
      $file = $file_A[$f];
      $filekey = $filekey_A[$f];
      $seq = undef;
      $ftr = ""; # current feature string (multiple lines)
      my $prv_line_seq = 0;
      %{$file_seq_ftr_HHH{$filekey}} = ();
      open(IN, $file) || ofile_FileOpenFailure($file, $sub_name, $!, "reading", $FH_HR);
      while(my $line = <IN>) { 
        if($line !~ m/\S+/) { # skip blank lines
          ; # do nothing
        }
        elsif($skip_msg_lines && ($line =~ m/^Additional/)) { 
          ; # do nothing
        }
        elsif($skip_msg_lines && ($line =~ m/^ERROR/)) { 
          ; # do nothing
        }
        elsif($skip_msg_lines && ($line =~ m/^WARNING/)) { 
          ; # do nothing
        }
        elsif($line =~ m/^\>/) { # sequence line
          if(defined $seq) { # if this isn't our first sequence
            # store final feature string from previous sequence
            $file_seq_ftr_HHH{$filekey}{$seq}{$ftr} = 1;
            if(! exists $seq_ftr_HH{$seq}{$ftr}) { 
              push(@{$seq_ftr_HA{$seq}}, $ftr); 
              $seq_ftr_HH{$seq}{$ftr} = 1;
            }
          }
          $ftr = "";
          $seq = $line;
          chomp $seq;
          if(! exists $seq_H{$seq}) { 
            push(@seq_A, $seq);
            $seq_H{$seq} = 1;
          }
          $prv_line_seq = 1;
          %{$file_seq_ftr_HHH{$filekey}{$seq}} = ();
          if(! exists $seq_ftr_HA{$seq}) { 
            @{$seq_ftr_HA{$seq}} = ();
          }
        }
        elsif($line =~ m/^\<?\d+\s+\d+\s+\S+/) { # coordinate plus feature line
          if(! $prv_line_seq) { # previous line was not a sequence line so
                                # we know we just ended a feature, store it
            $file_seq_ftr_HHH{$filekey}{$seq}{$ftr} = 1;
            if(! exists $seq_ftr_HH{$seq}{$ftr}) { 
              push(@{$seq_ftr_HA{$seq}}, $ftr); 
              $seq_ftr_HH{$seq}{$ftr} = 1;
            }
            $ftr = $line; # reset this to new coordinate line
          }
          else { # previous line was a sequence line, this will be first part of $ftr
            $ftr .= $line
          }
        }
        elsif($line =~ m/^\<?\d+\s+\d+/) { # coordinate line with no feature
          # add to current feature
          $ftr .= $line; 
          $prv_line_seq = 0;
        }
        else { # not a sequence line and not a coordinate line
          # just add to current feature string
          $ftr .= $line;
          $prv_line_seq = 0;
        }
      }
      # finished reading all lines, store the final feature
      if($ftr ne "") { 
        $file_seq_ftr_HHH{$filekey}{$seq}{$ftr} = 1;
        if(! exists $seq_ftr_HH{$seq}{$ftr}) { 
          push(@{$seq_ftr_HA{$seq}}, $ftr); 
          $seq_ftr_HH{$seq}{$ftr} = 1;
        }
      }
    }

    # now just compare all feature lines, and determine how many are in common
    # and how many are unique to either OUT or EXP
    foreach $seq (@seq_A) { 
      my $seq_id_ct = 0;
      my $seq_out_uniq_ct = 0;
      my $seq_exp_uniq_ct = 0;
      my $prefix = "";
      my @out_A = ();
      if((exists $file_seq_ftr_HHH{"out"}{$seq}) && 
         (exists $file_seq_ftr_HHH{"exp"}{$seq})) { 
        foreach $ftr (@{$seq_ftr_HA{$seq}}) { 
          if((exists $file_seq_ftr_HHH{"out"}{$seq}{$ftr}) && 
             (exists $file_seq_ftr_HHH{"exp"}{$seq}{$ftr})) { 
            $seq_id_ct++; 
            $prefix = "FTR-IDENTICAL: "; 
          }
          elsif(exists $file_seq_ftr_HHH{"out"}{$seq}{$ftr}) { 
              $seq_out_uniq_ct++; 
              $prefix = "FTR-OUT-ONLY: ";
          }
          elsif(exists $file_seq_ftr_HHH{"exp"}{$seq}{$ftr}) { 
              $seq_exp_uniq_ct++; 
              $prefix = "FTR-EXP-ONLY: ";
          }
          else { 
            ofile_FAIL("ERROR in $sub_name, second pass feature not in either file:\n$ftr", 1, $FH_HR);
          }
          foreach $line (split("\n", $ftr)) { 
            push(@out_A, $prefix . $line . "\n");
          }
        }
        # mark completely identical sequences with " *" at the end of the line summarizing them
        my $fully_identical_str = "";
        if(scalar(keys %{$file_seq_ftr_HHH{"out"}{$seq}} == scalar(keys %{$file_seq_ftr_HHH{"exp"}{$seq}})) && 
           scalar(keys %{$file_seq_ftr_HHH{"out"}{$seq}} == $seq_id_ct)) { 
          $fully_identical_str = " *IDENTICAL*"; 
        }
        printf DIFFOUT ("$seq FTR-OUT: %2d FTR-EXP: %2d FTR-ID: %2d FTR-OUT-UNIQUE: %2d FTR-EXP-UNIQUE: %2d%s\n", 
                        scalar(keys %{$file_seq_ftr_HHH{"out"}{$seq}}), 
                        scalar(keys %{$file_seq_ftr_HHH{"exp"}{$seq}}), 
                        $seq_id_ct, $seq_out_uniq_ct, $seq_exp_uniq_ct, 
                        $fully_identical_str);
        if($fully_identical_str eq "") { 
          foreach $line (@out_A) { print DIFFOUT $line; }
        }
      }
      elsif(exists $file_seq_ftr_HHH{"out"}{$seq}) { 
        print DIFFOUT "$seq OUT-ONLY\n";
      }
      elsif(exists $file_seq_ftr_HHH{"exp"}{$seq}) { 
        print DIFFOUT "$seq EXP-ONLY\n";
      }
    }
  }
  close(DIFFOUT);
  return;
}
    
    

