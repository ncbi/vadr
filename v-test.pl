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
#     option            type       default               group   requires incompat    preamble-output                       help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                "display this help",                                  \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",        "force; if dir <output directory> exists, overwrite it",   \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                         "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("-s",           "boolean", 0,                        1,    undef, "--rmout",  "skip commands, they were already run, just compare files",  "skip commands, they were already run, just compare files",   \%opt_HH, \@opt_order_A);
$opt_group_desc_H{"2"} = "options for defining variables in testing files";
#       option       type        default                group  requires incompat          preamble-output                                              help-output    
opt_Add("--dirbuild",   "string",  undef,                    2,   undef, undef,       "build directory, replaces !dirbuild! in test file with <s>", "build directory, replaces !dirbuild! in test file with <s>", \%opt_HH, \@opt_order_A);
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
                'dirbuild=s'   => \$GetOptions_H{"--dirbuild"},
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
my $version       = "1.1dev4";
my $releasedate   = "July 2020";
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

# read in the test file
my $progress_w = 50; # the width of the left hand column in our progress output, hard-coded
my $start_secs = ofile_OutputProgressPrior("Parsing test file", $progress_w, $log_FH, *STDOUT);
my @cmd_A      = (); # array of the commands to run
my @desc_A     = (); # array of the descriptions for the commands
my @outfile_AA = (); # array of arrays of output files to compare for each command
my @expfile_AA = (); # array of arrays of expected files to compare output to for each command
my @rmdir_AA   = (); # array of directories to remove after each command is completed
my $do_teamcity = opt_Get("--noteamcity", \%opt_HH) ? 0 : 1;
my $ncmd = parse_test_file($test_file, \@cmd_A, \@desc_A, \@outfile_AA, \@expfile_AA, \@rmdir_AA, \%opt_HH, $ofile_info_HH{"FH"});
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

my $npass = 0;
my $nfail = 0;
for(my $i = 1; $i <= $ncmd; $i++) { 
  my $cmd  = $cmd_A[($i-1)];
  my $desc = $desc_A[($i-1)];
  if($do_teamcity) { 
    ofile_OutputString($log_FH, 1, sprintf("##teamcity[testStarted name='$desc' captureStandardOutput='true']\n"));
  }
  my $outfile_AR = \@{$outfile_AA[($i-1)]};
  my $expfile_AR = \@{$expfile_AA[($i-1)]};
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
    my $diff_file = $out_root . "." . $i . "." . ($j+1) . ".diff";
    my $pass = diff_two_files($outfile_AR->[$j], $expfile_AR->[$j], $diff_file, \%opt_HH, $ofile_info_HH{"FH"});
    if($pass) { $npass++; $npass_i++; }
    else      { $nfail++; $nfail_i++; }
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

exit 0;
#################################################################
# Subroutine:  parse_test_file()
# Incept:      EPN, Wed May 16 16:03:40 2018
#
# Purpose:     Parse an input test file and store the relevant information
#              in passed in array references.
#
# Arguments:
#   $testfile:    name of file to parse
#   $cmd_AR:      ref to array of commands to fill here
#   $desc_AR:     ref to array of descriptions to fill here
#   $outfile_AAR: ref to 2D array of output files to fill here
#   $expfile_AAR: ref to 2D array of expected files to fill here
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
#################################################################
sub parse_test_file { 
  my $sub_name = "parse_test_file";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($testfile, $cmd_AR, $desc_AR, $outfile_AAR, $expfile_AAR, $rmdir_AAR, $opt_HHR, $FH_HR) = @_;

  open(IN, $testfile) || ofile_FileOpenFailure($testfile, $sub_name, $!, "reading", $FH_HR);
  my $ncmd = 0;
  my $ndesc = 0;
  my $outfile;
  my $expfile;
  my $rmdir;

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
          if(! (@{$expfile_AAR->[($ncmd-1)]})) { ofile_FAIL("ERROR did not read any exp: lines for command " . ($ncmd+1), 1, $FH_HR); }
          my $nout_prv = scalar(@{$outfile_AAR->[($ncmd-1)]});
          my $nexp_prv = scalar(@{$expfile_AAR->[($ncmd-1)]});
          if($nout_prv != $nexp_prv) { 
            ofile_FAIL("ERROR different number of output and expected lines for command " . ($ncmd+1), 1, $FH_HR);
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
        push(@{$desc_AR}, $desc); 
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
          if(-e $outfile) { utl_FileRemoveUsingSystemRm($outfile, "v-test.pl", \%opt_HH, $FH_HR); }
        }
        else { # -s not used, --rmout not used
          if(-e $outfile) { ofile_FAIL("ERROR, output file $outfile already exists (and -s/--rmout not used)", 1, $FH_HR); }
        }
      }
      elsif($line =~ s/^exp\:\s+//) { 
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
  if(! (@{$expfile_AAR->[($ncmd-1)]})) { ofile_FAIL("ERROR did not read any exp: lines for command " . ($ncmd+1), 1, $FH_HR); }
  my $nout_prv = scalar(@{$outfile_AAR->[($ncmd-1)]});
  my $nexp_prv = scalar(@{$expfile_AAR->[($ncmd-1)]});
  if($nout_prv != $nexp_prv) { 
    ofile_FAIL("ERROR different number of output and expected lines for command " . ($ncmd+1), 1, $FH_HR);
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

#################################################################
# Subroutine:  OLD_compare_two_sqtable_files()
# Incept:      EPN, Mon Jun 11 09:42:11 2018
#
# Purpose:     Compare two sqtable files outputting the number of 
#              lost, added, and changed features.
#
# Arguments:
#   $out_file:    name of output sqtable file
#   $exp_file:    name of expected sqtable file
#   $diff_file:   name of file to create with differences
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       If sequences are not in the same order in the two files
#
#################################################################
sub OLD_compare_two_sqtable_files { 
  my $sub_name = "OLD_compare_two_sqtable_files";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_file, $exp_file, $diff_file, $FH_HR) = @_;

  my $out_file_exists   = (-e $out_file) ? 1 : 0;
  my $exp_file_exists   = (-e $exp_file) ? 1 : 0;
  my $out_file_nonempty = (-s $out_file) ? 1 : 0;
  my $exp_file_nonempty = (-s $exp_file) ? 1 : 0;

  my $conclusion = "";
  my $pass = 0;

  if(! $exp_file_exists) { 
    ofile_FAIL("ERROR in $sub_name, expected file $exp_file does not exist", 1, $FH_HR) ;
  }
  if(! $exp_file_nonempty) { 
    ofile_FAIL("ERROR in $sub_name, expected file $exp_file exists but is empty", 1, $FH_HR);
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

  if($out_file_nonempty) { 
    open(OUT, $out_file) || ofile_FileOpenFailure($out_file, $sub_name, $!, "reading", $FH_HR);
    $lidx = 0;
    while($out_line = <OUT>) { 
      if($out_line =~ m/^\>/) { 
        $seq = $out_line;
        chomp $seq;
        push(@out_seq_A, $seq);
        $out_seq_lidx_H{$seq} = $lidx;
        if(! exists $seq_exists_H{$seq}) { 
          $seq_exists_H{$seq} = 1;
        }
      }
      push(@out_line_A, $out_line);
      $lidx++;
    }
    close(OUT);

    open(EXP, $exp_file) || ofile_FileOpenFailure($exp_file, $sub_name, $!, "reading", $FH_HR);
    $lidx = 0;
    while($exp_line = <EXP>) { 
      if($exp_line =~ m/^\>/) { 
        $seq = $exp_line;
        chomp $seq;
        push(@exp_seq_A, $seq);
        $exp_seq_lidx_H{$seq} = $lidx;
        if(! exists $seq_exists_H{$seq}) { 
          $seq_exists_H{$seq} = 1;
        }
      }
      push(@exp_line_A, $exp_line);
      $lidx++;
    }
    close(EXP);

    # make sure all sequences existed in both files
    foreach $seq (sort keys %seq_exists_H) { 
      if(! exists $out_seq_lidx_H{$seq}) { 
        ofile_FAIL("ERROR in $sub_name, $seq not in out_seq_lidx_H", 1, $FH_HR);
      }
      if(! exists $exp_seq_lidx_H{$seq}) { 
        ofile_FAIL("ERROR in $sub_name, $seq not in exp_seq_lidx_H", 1, $FH_HR);
      }
    }
    # make sure all sequences are in the same order in both files
    $nseq = scalar(@out_seq_A);
    my $s;
    for($s = 0; $s < $nseq; $s++) { 
      if($out_seq_A[$s] ne $exp_seq_A[$s]) { 
        ofile_FAIL("ERROR in $sub_name, $seq not in same order in both files", 1, $FH_HR);
      }
      printf("HEYA $out_seq_A[$s] " . $out_seq_lidx_H{$out_seq_A[$s]} . " " . $exp_seq_lidx_H{$out_seq_A[$s]} . "\n");
    }
    
    open(DIFFOUT, ">", $diff_file) || ofile_FileOpenFailure($diff_file, $sub_name, $!, "writing", $FH_HR);

    # for each sequence, compare the annotations
    my $lidx;        # counter over lines
    my $first_lidx;  # first line index of annotation for a sequence
    my $final_lidx;  # final line index of annotation for a sequence
    my %cur_out_H = ();    # key: feature line, value: current annotation in out file for that feature
    my %cur_exp_H = ();    # key: feature line, value: current annotation in exp file for that feature
    my %cur_exists_H = (); # key: feature line, value '1'
    my @cur_A = ();        # array of feature lines in order they were read

    for($s = 0; $s < $nseq; $s++) { 
      $seq = $out_seq_A[$s];
      printf("\ns: $s\n");
      %cur_exists_H = ();
      @cur_A = ();

      # parse information from out_file
      %cur_out_H = ();
      $first_lidx = $out_seq_lidx_H{$seq} + 1; # first line after a sequence is the beginning of its annotation
      if($s < ($nseq-1)) { $final_lidx = $out_seq_lidx_H{$out_seq_A[($s+1)]} - 1; }
      else               { $final_lidx = scalar(@out_line_A) - 1; }
      my $fline = undef;
      my $cline = undef;
      my $cur_cline = "";
      printf("OUT $first_lidx..$final_lidx\n");
      for($lidx = $first_lidx; $lidx <= $final_lidx; $lidx++) { 
        printf("OUT line: $out_line_A[$lidx]\n");
        if($out_line_A[$lidx] =~ m/^\<?\d+/) { 
          # coordinate line
          $cline = $out_line_A[$lidx];
          $cur_cline .= $cline;
        }
        else { 
          # feature line
          $fline = $out_line_A[$lidx];
          $cur_out_H{$fline} = $cur_cline;
          printf("HEYA-OUT read fline: $fline added cline: $cur_cline\n");
          if(! exists $cur_exists_H{$fline}) { 
            push(@cur_A, $fline);
            $cur_exists_H{$fline} = 1;
          }
          $cur_cline = "";
        }
      }
      exit 0;

      # parse information from exp_file
      %cur_exp_H = ();
      $first_lidx = $exp_seq_lidx_H{$seq} + 1; # first line after a sequence is the beginning of its annotation
      if($s < ($nseq-1)) { $final_lidx = $exp_seq_lidx_H{$exp_seq_A[($s+1)]} - 1; }
      else               { $final_lidx = scalar(@exp_line_A) - 1; }
      $fline = undef;
      $cline = undef;
      $cur_cline = "";
      printf("EXP $first_lidx..$final_lidx\n");
      for($lidx = $first_lidx; $lidx <= $final_lidx; $lidx++) { 
        printf("EXP line: $exp_line_A[$lidx]\n");
        if($exp_line_A[$lidx] =~ m/^\<?\d+/) { 
          # coordinate line
          $cline = $exp_line_A[$lidx];
          $cur_cline .= $cline;
        }
        else { 
          # feature line
          $fline = $exp_line_A[$lidx];
          $cur_exp_H{$fline} = $cur_cline;
          if(! exists $cur_exists_H{$fline}) { 
            push(@cur_A, $fline);
            $cur_exists_H{$fline} = 1;
          }
          $cur_cline = "";
        }
      }
      
      # compare annotation from exp_file and out_file
      my $cur_ftr_nidentical  = 0; # number of features that are identical (all lines) between output and expected
      my $cur_ftr_ndiff       = 0; # number of features that are different (at least one line different) between output and expected
      my $cur_ftr_nout_only   = 0; # number of features that are only in output file
      my $cur_ftr_nexp_only   = 0; # number of features that are only in expected file
      my $cur_nftr            = 0; # number of features in the sequence
      my $cur_outstr          = ""; # output string

      my $cur_note_nidentical  = 0; # number of 'notes' that are identical (all lines) between output and expected
      my $cur_note_ndiff       = 0; # number of 'notes' that are different (at least one line different) between output and expected
      my $cur_note_nout_only   = 0; # number of 'notes' that are only in output file
      my $cur_note_nexp_only   = 0; # number of 'notes' that are only in expected file
      my $cur_nnote            = 0; # number of 'notes' in the sequence

      my $is_note = 0; # '1' if this line is a note, '0' if it is a feature

      foreach $fline (@cur_A) { 
        # determine if it is a 'note' or not, we treat them differently
        $is_note = ($fline =~ m/^\s+note/) ? 1 : 0;

        if($is_note) { 
          $cur_nnote++;
        }
        else { 
          $cur_nftr++;
        }

        if((! exists $cur_exp_H{$fline}) && (! exists $cur_out_H{$fline})) { 
          ofile_FAIL("ERROR in $sub_name, $fline does not exist in either exp or out", 1, $FH_HR);
        }
        elsif((exists $cur_exp_H{$fline}) && (! exists $cur_out_H{$fline})) { 
          if($is_note) { $cur_note_nexp_only++; }
          else         { $cur_ftr_nexp_only++; }
          $cur_outstr .= "EXP only: $fline$cur_exp_H{$fline}"; 
        }
        elsif((! exists $cur_exp_H{$fline}) && (exists $cur_out_H{$fline})) { 
          if($is_note) { $cur_note_nout_only++; }
          else         { $cur_ftr_nout_only++; }
          $cur_outstr .= "OUT only: $fline$cur_out_H{$fline}"; 
        }
        elsif($cur_exp_H{$fline} ne $cur_out_H{$fline}) { 
          if($is_note) { $cur_note_ndiff++; }
          else         { $cur_ftr_ndiff++; }
          if(! exists $cur_exp_H{$fline}) { 
            printf("exp does not exists $fline\n");
          }
          if(! exists $cur_out_H{$fline}) { 
            printf("out does not exists $fline\n");
          }
          $cur_outstr .= "DIFFERENT: " . $fline . "OUT:\n$cur_out_H{$fline}EXP:\n$cur_exp_H{$fline}"; 
        }
        else { 
          if($is_note) { $cur_note_nidentical++; }
          else         { $cur_ftr_nidentical++; }
        }
      }
      printf DIFFOUT ("%s %s %-50s %2d features [%2d %2d %2d %2d] %2d notes [%2d %2d %2d %2d]\n", 
                      ($cur_ftr_nidentical == $cur_nftr)   ? "FTR-IDENTICAL" : "FTR-DIFFERENT", 
                      ($cur_note_nidentical == $cur_nnote) ? "NOTE-IDENTICAL" : "NOTE-DIFFERENT", 
                      $seq,
                      $cur_nftr, $cur_ftr_nidentical, $cur_ftr_ndiff, $cur_ftr_nout_only, $cur_ftr_nexp_only, 
                      $cur_nnote, $cur_note_nidentical, $cur_note_ndiff, $cur_note_nout_only, $cur_note_nexp_only);
      if($cur_outstr ne "") { 
        printf DIFFOUT $cur_outstr;
      }
      if($cur_ftr_nidentical == $cur_nftr) { $tot_nseq_ftr_identical++; }
      if($cur_ftr_ndiff > 0)               { $tot_nseq_ftr_diff++; }
      if($cur_ftr_nout_only > 0)           { $tot_nseq_ftr_out_only++; }
      if($cur_ftr_nexp_only > 0)           { $tot_nseq_ftr_exp_only++; }

      if($cur_note_nidentical == $cur_nftr) { $tot_nseq_note_identical++; }
      if($cur_note_ndiff > 0)               { $tot_nseq_note_diff++; }
      if($cur_note_nout_only > 0)           { $tot_nseq_note_out_only++; }
      if($cur_note_nexp_only > 0)           { $tot_nseq_note_exp_only++; }
    }
  }

  printf DIFFOUT ("SUMMARY %5d sequences FEATURES: [identical:%5d;  different:%5d;  out-only:%5d;  exp-only:%5d;] NOTES: [identical:%5d;  different:%5d;  out-only:%5d;  exp-only:%5d;]\n", $nseq,
                  $tot_nseq_ftr_identical, $tot_nseq_ftr_diff, $tot_nseq_ftr_out_only, $tot_nseq_ftr_exp_only, 
                  $tot_nseq_note_identical, $tot_nseq_note_diff, $tot_nseq_note_out_only, $tot_nseq_note_exp_only);

  close(DIFFOUT);
  return;
}
    
    

