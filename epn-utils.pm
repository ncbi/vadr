#!/usr/bin/perl
#
# epn-utils.pm
# Eric Nawrocki
# EPN, Tue Mar 19 13:35:06 2019
# version: 0.00
#
use strict;
use warnings;
use Time::HiRes qw(gettimeofday);

require "epn-ofile.pm";

#####################################################################
# Data structures used in this module:
#
#####################################################################
# List of subroutines:
# 
#################################################################
# Subroutine: utl_RunCommand()
# Incept:     EPN, Thu Feb 11 13:32:34 2016 [dnaorg.pm]
#
# Purpose:     Runs a command using system() and exits in error 
#              if the command fails. If $be_verbose, outputs
#              the command to stdout. If $FH_HR->{"cmd"} is
#              defined, outputs command to that file handle.
#
# Arguments:
#   $cmd:         command to run, with a "system" command;
#   $be_verbose:  '1' to output command to stdout before we run it, '0' not to
#   $do_failok:   '1' to NOT exit if command fails, '0' to exit if command fails
#   $FH_HR:       REF to hash of file handles, including "cmd"
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails and $do_failok is '0'
#################################################################
sub utl_RunCommand {
  my $sub_name = "utl_RunCommand()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $be_verbose, $do_failok, $FH_HR) = @_;
  
  my $cmd_FH = undef;
  if(defined $FH_HR && defined $FH_HR->{"cmd"}) { 
    $cmd_FH = $FH_HR->{"cmd"};
  }

  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  if(defined $cmd_FH) { 
    print $cmd_FH ("$cmd\n");
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));

  system($cmd);

  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if(($? != 0) && (! $do_failok)) { 
    ofile_FAIL("ERROR in $sub_name, the following command failed:\n$cmd\n", $?, $FH_HR); 
  }

  return ($stop_time - $start_time);
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
