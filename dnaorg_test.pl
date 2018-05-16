#!/usr/bin/env perl
# 
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

require "dnaorg.pm"; 
require "epn-options.pm";

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
my $inf_exec_dir      = $dnaorgdir . "/infernal-1.1.2/src/";
my $hmmer_exec_dir    = $dnaorgdir . "/hmmer-3.1b2/src/";
my $esl_exec_dir      = $dnaorgdir . "/infernal-1.1.2/easel/miniapps/";
my $esl_fetch_cds     = $dnaorgdir . "/esl-fetch-cds/esl-fetch-cds.pl";
my $esl_epn_translate = $dnaorgdir . "/esl-epn-translate/esl-epn-translate.pl";
my $esl_ssplit        = $dnaorgdir . "/Bio-Easel/scripts/esl-ssplit.pl";

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
#     option            type       default               group   requires incompat    preamble-output                          help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                   "display this help",                                  \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                   "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",  "force; if dir <output directory> exists, overwrite it",   \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_test.pl [-options] <input test file e.g. testfiles/testin.1> <output directory to create>\n";
my $synopsis = "dnaorg_test.pl :: test dnaorg scripts [TEST SCRIPT]";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"},
                'v'            => \$GetOptions_H{"-v"},
                'f'            => \$GetOptions_H{"-f"}); 

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.31";
my $releasedate   = "May 2018";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date);
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
my ($test_file, $dir_out) = (@ARGV);

if(defined $dir_out) { 
  $dir_out =~ s/\/$//; # remove final '/' if there is one
}
my $dir_out_tail   = $dir_out;
$dir_out_tail   =~ s/^.+\///; # remove all but last dir
my $out_root   = $dir_out .   "/" . $dir_out_tail   . ".dnaorg_test";

my $cmd;
my @early_cmd_A = ();  # array of commands we run before our log file is opened

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

if(-d $dir_out) { 
  $cmd = "rm -rf $dir_out";
  if(opt_Get("-f", \%opt_HH)) { # -f used, always remove it
    runCommand($cmd, opt_Get("-v", \%opt_HH), undef); push(@early_cmd_A, $cmd); 
  }
  else { # dirout exists but -f not used
    die "ERROR directory named $dir_out (specified with --dirout) already exists. Remove it, or use -f to overwrite it."; 
  }
}

if(! -d $dir_out) {
  $cmd = "mkdir $dir_out";
  runCommand($cmd, 0, undef);
}

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("test file");
my @arg_A      = ($test_file);
outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date);
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

# now we have the log file open, output the banner there too
outputBanner($log_FH, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

# read in the test file
my @cmd_A      = (); # array of the commands to run
my @outfile_AA = (); # array of arrays of output files to compare for each command
my @expfile_AA = (); # array of arrays of expected files to compare output to for each command
my @rmdir_AA   = (); # array of directories to remove after each command is completed
my $ncmd = parse_test_file($test_file, \@cmd_A, \@outfile_AA, \@expfile_AA, \@rmdir_AA, $ofile_info_HH{"FH"});

for(my $i = 1; $i <= $ncmd; $i++) { 
  my $cmd = $cmd_A[($i-1)];
  my $outfile_AR = \@{$outfile_AA[($i-1)]};
  my $expfile_AR = \@{$expfile_AA[($i-1)]};
  my $progress_w = 50; # the width of the left hand column in our progress output, hard-coded
  my $start_secs = outputProgressPrior("Running cmd $i ...\n" . $cmd_A[($i-1)], $progress_w, $log_FH, *STDOUT);
  #runCommand($cmd, 0, $ofile_info_HHR->{"FH"});
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  my $nout = scalar(@{$outfile_AR});
  for(my $j = 0; $j < $nout; $j++) { 
    #diff_files($outfile_AR->[$j], $expfile_AR->[$j]);
    printf("diff_files($outfile_AR->[$j], $expfile_AR->[$j]);\n");
  }
}

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
#   $outfile_AAR: ref to 2D array of output files to fill here
#   $expfile_AAR: ref to 2D array of expected files to fill here
#   $rmdir_AAR:   ref to 2D array of directories to remove after calling each command
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
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
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($testfile, $cmd_AR, $outfile_AAR, $expfile_AAR, $rmdir_AAR, $FH_HR) = @_;

  open(IN, $testfile) || fileOpenFailure($testfile, $sub_name, $!, "reading", $FH_HR);
  my $ncmd = 0;
  my $outfile;
  my $expfile;
  my $rmdir;

  while(my $line = <IN>) { 
    if($line !~ m/^\#/) { 
      # example input file:
      # # comment line (ignored)
      # command: perl $DNAORGDIR/dnaorg_scripts/dnaorg_classify.pl -f -A /panfs/pan1/dnaorg/virseqannot/dnaorg-build-directories/norovirus-builds --infasta testfiles/noro.9.fa --dirbuild /panfs/pan1/dnaorg/virseqannot/dnaorg-build-directories/norovirus-builds --dirout test-noro.9
      # out: test-noro.9/test-noro.9-NC_001959.dnaorg_annotate.sqtable 
      # exp: testfiles/testout.1/test-noro.9/test-noro.9-NC_001959.dnaorg_annotate.sqtable 
      chomp $line;
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists

      if($line =~ s/^command\:\s+//) { 
        my $cmd = $line;
        if($ncmd > 0) { 
          # make sure we have 
          if(! (@{$outfile_AAR->[($ncmd-1)]})) { DNAORG_FAIL("ERROR did not read any out: lines for command " . ($ncmd+1), 1, $FH_HR); }
          if(! (@{$expfile_AAR->[($ncmd-1)]})) { DNAORG_FAIL("ERROR did not read any exp: lines for command " . ($ncmd+1), 1, $FH_HR); }
          my $nout_prv = scalar(@{$outfile_AAR->[($ncmd-1)]});
          my $nexp_prv = scalar(@{$expfile_AAR->[($ncmd-1)]});
          if($nout_prv != $nexp_prv) { 
            DNAORG_FAIL("ERROR different number of output and expected lines for command " . ($ncmd+1), 1, $FH_HR);
          }
        }
        push(@{$cmd_AR}, $cmd); 
        $ncmd++;
      }
      elsif($line =~ s/^out\:\s+//) { 
        $outfile = $line;
        $outfile =~ s/^\s+//;
        $outfile =~ s/\s+$//;
        if($outfile =~ m/\s/) { DNAORG_FAIL("ERROR output file has spaces: $outfile", 1, $FH_HR); }
        if(scalar(@{$outfile_AAR}) < $ncmd) { 
          @{$outfile_AAR->[($ncmd-1)]} = ();
        }
        push(@{$outfile_AAR->[($ncmd-1)]}, $outfile);
        if(-e $outfile) { DNAORG_FAIL("ERROR, output file $outfile already exists", 1, $FH_HR); }
      }
      elsif($line =~ s/^exp\:\s+//) { 
        $expfile = $line;
        $expfile =~ s/^\s+//;
        $expfile =~ s/\s+$//;
        if($expfile =~ m/\s/) { DNAORG_FAIL("ERROR expected file has spaces: $expfile", 1, $FH_HR) }
        if(scalar(@{$expfile_AAR}) < $ncmd) { 
          @{$expfile_AAR->[($ncmd-1)]} = ();
        }
        push(@{$expfile_AAR->[($ncmd-1)]}, $expfile);
        # TEMP if(! -e $expfile) { DNAORG_FAIL("ERROR, expected file $expfile does not exist", 1, $FH_HR); }
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
        DNAORG_FAIL("ERROR unable to parse line $line in $testfile", 1, $FH_HR);
      }
    }
  }
  close(IN);

  # for final command, check that number of exp and out files is equal
  if(! (@{$outfile_AAR->[($ncmd-1)]})) { DNAORG_FAIL("ERROR did not read any out: lines for command " . ($ncmd+1), 1, $FH_HR); }
  if(! (@{$expfile_AAR->[($ncmd-1)]})) { DNAORG_FAIL("ERROR did not read any exp: lines for command " . ($ncmd+1), 1, $FH_HR); }
  my $nout_prv = scalar(@{$outfile_AAR->[($ncmd-1)]});
  my $nexp_prv = scalar(@{$expfile_AAR->[($ncmd-1)]});
  if($nout_prv != $nexp_prv) { 
    DNAORG_FAIL("ERROR different number of output and expected lines for command " . ($ncmd+1), 1, $FH_HR);
  }

  return $ncmd;
}
