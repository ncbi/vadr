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

#################################################################
# Subroutine: utl_ArrayOfHashesToArray()
# Incept:     EPN, Wed Mar 20 09:07:06 2019
#
# Purpose:    Fill @{$AR} with all values in $AHR->[]{$key}.
# Arguments:
#   $AHR:      REF to array of hashes
#   $AR:       REF to array to add to
#   $key:      key of interest
# 
# Returns: number of elements added to @{$AR}
#
#################################################################
sub utl_ArrayOfHashesToArray {
  my $sub_name = "utl_ArrayOfHashesToArray()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($AHR, $AR, $key) = (@_);

  my $ret_n = 0;
  my $n = scalar(@{$AHR});
  for(my $i = 0; $i < $n; $i++) { 
    if(defined $AHR->[$i]{$key}) { 
      push(@{$AR}, $AHR->[$i]{$key}); 
      $ret_n++;
    }
  }
  
  return $ret_n;
}

#################################################################
# Subroutine:  utl_ConcatenateListOfFiles()
# Incept:      EPN, Sun Apr 24 08:08:15 2016
#
# Purpose:     Concatenate a list of files into one file.
#              If the list has more than 500 files, split
#              up job into concatenating 500 at a time.
# 
#              We remove all files that we concatenate unless
#              --keep option is on in %{$opt_HHR}.
#
# Arguments: 
#   $file_AR:          REF to array of all files to concatenate
#   $outfile:          name of output file to create by concatenating
#                      all files in @{$file_AR}.
#   $caller_sub_name:  name of calling subroutine (can be undef)
#   $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If one of the cat commands fails.
#              If $outfile is in @{$file_AR}
#              If @{$file_AR} contains more than 800*800 files
#              (640K) if so, we may need to call this function
#              recursively twice (that is, recursive call will
#              also call itself recursively) and we don't have 
#              a sophisticated enough temporary file naming
#              strategy to handle that robustly.
################################################################# 
sub utl_ConcatenateListOfFiles { 
  my $nargs_expected = 5;
  my $sub_name = "utl_ConcatenateListOfFiles()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($file_AR, $outfile, $caller_sub_name, $opt_HHR, $FH_HR) = @_;

  if(utl_ArrayFindNonNumericValue($file_AR, $outfile, $FH_HR) != -1) { 
    ofile_FAIL(sprintf("ERROR in $sub_name%s, output file name $outfile exists in list of files to concatenate", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), "dnaorg", 1, $FH_HR);
  }

  # first, convert @{$file_AR} array into a 2D array of file names, each of which has 
  # a max of 800 elements, we'll concatenate each of these lists separately
  my $max_nfiles = 800;
  my $nfiles = scalar(@{$file_AR});

  if($nfiles > ($max_nfiles * $max_nfiles)) { 
    ofile_FAIL(sprintf("ERROR in $sub_name%s, trying to concatenate %d files, our limit is %d", 
                       (defined $caller_sub_name) ? " called by $caller_sub_name" : "", $nfiles, $max_nfiles * $max_nfiles), 
               undef, 1, $FH_HR);
  }
    
  my ($idx1, $idx2); # indices in @{$file_AR}, and of secondary files
  my @file_AA = ();
  $idx2 = -1; # get's incremented to 0 in first loop iteration
  for($idx1 = 0; $idx1 < $nfiles; $idx1++) { 
    if(($idx1 % $max_nfiles) == 0) { 
      $idx2++; 
      @{$file_AA[$idx2]} = (); # initialize
    }
    push(@{$file_AA[$idx2]}, $file_AR->[$idx1]);
  }
  
  my $nconcat = scalar(@file_AA);
  my @tmp_outfile_A = (); # fill this with names of temporary files we create
  my $tmp_outfile; # name of an output file we'll create
  for($idx2 = 0; $idx2 < $nconcat; $idx2++) { 
    if($nconcat == 1) { # special case, we don't need to create any temporary files
      $tmp_outfile = $outfile;
    }
    else { 
      $tmp_outfile = $outfile . ".tmp" . ($idx2+1); 
      # make sure this file does not exist in @{$file_AA[$idx2]} to avoid klobbering
      # if it does, continue to append .tmp($idx2+1) until it doesn't
      while(utl_ArrayFindNonNumericValue($file_AA[$idx2], $tmp_outfile, $FH_HR) != -1) { 
        $tmp_outfile .= ".tmp" . ($idx2+1); 
      }
    }
    # create the concatenate command
    my $cat_cmd = "cat ";
    foreach my $tmp_file (@{$file_AA[$idx2]}) {
      $cat_cmd .= $tmp_file . " ";
    }
    $cat_cmd .= "> $tmp_outfile";

    # execute the command
    utl_RunCommand($cat_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

    # add it to the array of temporary files
    push(@tmp_outfile_A, $tmp_outfile); 
  }

  if(scalar(@tmp_outfile_A) > 1) { 
    # we created more than one temporary output file, concatenate them
    # by calling this function again
    utl_ConcatenateListOfFiles(\@tmp_outfile_A, $outfile, (defined $caller_sub_name) ? $caller_sub_name . ":" . $sub_name : $sub_name, $opt_HHR, $FH_HR);
  }

  if(! opt_Get("--keep", $opt_HHR)) { 
    # remove all of the original files, be careful to not remove @tmp_outfile_A
    # because the recursive call will handle that
    foreach my $file_to_remove (@{$file_AR}) { 
      utl_FileRemoveUsingSystemRm($file_to_remove, 
                                  (defined $caller_sub_name) ? $caller_sub_name . ":" . $sub_name : $sub_name, 
                                  $opt_HHR, $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine:  utl_ArrayFindNonNumericValue()
# Incept:      EPN, Tue Feb 16 10:40:57 2016 [dnaorg.pm]
#
# Purpose:     Returns (first) index in @{$AR} that has the 
#              nonnumeric value $value. Returns -1 
#              if it does not exist.
#
# Arguments: 
#   $AR:       REF to array 
#   $value:    the value we're checking exists in @{$AR}
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns:     index ($i) '1' if $value exists in @{$AR}, '-1' if not
#
# Dies:        if $value is numeric, or @{$AR} is not defined.
################################################################# 
sub utl_ArrayFindNonNumericValue { 
  my $nargs_expected = 3;
  my $sub_name = "utl_ArrayFindNonNumericValue()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR, $value, $FH_HR) = @_;

  if(verify_real($value)) { 
    ofile_FAIL("ERROR in $sub_name, value $value seems to be numeric, we can't compare it for equality", undef, 1, $FH_HR);
  }

  if(! defined $AR) { 
    ofile_FAIL("ERROR in $sub_name, array reference is not defined", undef, 1, $FH_HR);
  }

  for(my $i = 0; $i < scalar(@{$AR}); $i++) {
    if($AR->[$i] eq $value) { 
      return $i; 
    }
  }

  return -1; # did not find it
}

#################################################################
# Subroutine: utl_FileRemoveUsingSystemRm
# Incept:     EPN, Fri Mar  4 15:57:25 2016 [dnaorg.pm]
#
# Purpose:    Remove a file from the filesystem by using
#             the system rm command.
# Arguments:
#   $file:            file to remove
#   $caller_sub_name: name of caller, can be undef
#   $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    - if the file does not exist
#
#################################################################
sub utl_FileRemoveUsingSystemRm {
  my $sub_name = "utl_FileRemoveUsingSystemRm";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($file, $caller_sub_name, $opt_HHR, $FH_HR) = (@_);
  
  if(! -e $file) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, %s trying to remove file $file but it does not exist", 
                (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), "dnaorg", 1, $FH_HR); 
  }

  utl_RunCommand("rm $file", opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine:  utl_RemoveDirPath()
# Incept:      EPN, Mon Nov  9 14:30:59 2009 [ssu-align]
#
# Purpose:     Given a full path of a file remove the directory path.
#              For example: "foodir/foodir2/foo.stk" becomes "foo.stk".
#
# Arguments: 
#   $fullpath: name of original file
# 
# Returns:     The string $fullpath with dir path removed.
#
################################################################# 
sub utl_RemoveDirPath {
  my $sub_name = "utl_RemoveDirPath()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($fullpath) = (@_);

  $fullpath =~ s/^.+\///;

  return $fullpath;
}


#################################################################
# Subroutine:  utl_HHFromAH()
# Incept:      EPN, Thu Mar 21 06:18:38 2019
#
# Purpose:     Create a 2D hash %{$HHR} from an array of hashes
#              @{$AHR} using $AHR->[]{$key} as the 1D key in
#              %{$HHR}.
#              
# Arguments: 
#   $HHR:      ref to 2D hash to create
#   $AHR:      ref to array of hashes to copy from
#   $key:      key to use as 1st dim value
#   $FH_HR:    ref to hash of file handles, including "log" and "cmd"
# 
# Returns:     void
# 
# Dies:        If not all elements of @{$AHR} have @{$AHR->[]{$key}}
#              defined. If more than one elements of @{$AHR} have same
#              value for @{$AHR->[]{$key}}.
#
################################################################# 
sub utl_HHFromAH {
  my $sub_name = "utl_HHFromAH()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HHR, $AHR, $key, $caller_sub_name, $FH_HR);

  %{$HHR} = ();
  my $n = scalar(@{$AHR});
  for(my $i = 0; $i < $n; $i++) { 
    if(! defined $AHR->[$i]{$key}) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, %s trying to remove file $file but it does not exist", 
                         (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), "dnaorg", 1, $FH_HR); 
      

  return $fullpath;
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
