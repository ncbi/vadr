#!/usr/bin/perl
#
# epn-seq.pm
# Eric Nawrocki
# EPN, Tue Mar 19 13:29:58 2019
# version: 0.00
#
use strict;
use warnings;
use Time::HiRes qw(gettimeofday);

require "epn-utils.pm";
require "epn-ofile.pm";
require "epn-options.pm";

#####################################################################
# Data structures used in this module:
#
#####################################################################
#
# List of subroutines:
# 
#################################################################
# Subroutine : seq_ProcessSequenceFile()
# Incept:      EPN, Fri May 12 10:08:47 2017 [ribo.pm]
#
# Purpose:     Use esl-seqstat to get the lengths of all sequences in a
#              FASTA or Stockholm formatted sequence file and fill
#              @{$seq_order_AR} and %{$seq_len_HR}.
#
#              
# Arguments: 
#   $seqstat_exec:   path to esl-seqstat executable
#   $seq_file:       sequence file to process
#   $seqstat_file:   path to esl-seqstat output to create
#   $seq_name_AR:    ref to array of sequence names in order to fill here
#   $seq_len_HR:     ref to hash of sequence names
#   $opt_HHR:        reference to 2D hash of cmdline options
#   $ofile_info_HHR: ref to the ofile info 2D hash
# 
# Returns:     total number of nucleotides in all sequences read, 
#
# Dies:        If the sequence file has two sequences with identical names.
#              Error message will list all duplicates.
#              If no sequences were read.
#
################################################################# 
sub seq_ProcessSequenceFile { 
  my $nargs_expected = 7;
  my $sub_name = "seq_ProcessSequenceFile()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($seqstat_exec, $seq_file, $seqstat_file, $seq_name_AR, $seq_len_HR, $opt_HHR, $ofile_info_HHR) = (@_);

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience

  utl_RunCommand($seqstat_exec . " --dna -a $seq_file > $seqstat_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);
  if(defined $ofile_info_HHR) { 
    ofile_AddClosedFileToOutputInfo($ofile_info_HHR, undef, "seqstat", $seqstat_file, 0, "esl-seqstat -a output for $seq_file");
  }

  return seq_ParseSeqstatAFile($seqstat_file, $seq_name_AR, $seq_len_HR, $FH_HR);
}

#################################################################
# Subroutine : seq_ParseSeqstatAFile()
# Incept:      EPN, Wed Dec 14 16:16:22 2016 [ribo.pm]
#
# Purpose:     Parse an esl-seqstat -a output file.
#              
# Arguments: 
#   $seqstat_file:  file to parse
#   $seq_name_AR:   REF to array of sequences in order to fill here
#   $seq_len_HR     REF to hash of sequence names to fill here
#   $FH_HR:         REF to hash of file handles, including "cmd"
#
# Returns:     Total number of nucleotides read (summed length of all sequences). 
# 
# Dies:        If the sequence file has two sequences with identical names.
#              Error message will list all duplicates.
#              If no sequences were read.
#
################################################################# 
sub seq_ParseSeqstatAFile { 
  my $nargs_expected = 4;
  my $sub_name = "ribo_ParseSeqstatFile";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seqstat_file, $seq_name_AR, $seq_len_HR, $FH_HR) = @_;

  open(IN, $seqstat_file) || ofile_FileOpenFailure($seqstat_file, undef, $sub_name, $!, "reading", $FH_HR);

  my $nread = 0;            # number of sequences read
  my $tot_length = 0;       # summed length of all sequences
  my $seq_name;             # a sequence name
  my $length;               # length of a seq
  my %seq_exists_H = ();    # key is name of a read sequence, value is always '1'
  my %seq_dups_H = ();      # key is a sequence name that exists more than once in seq file, value is number of occurences
  my $at_least_one_dup = 0; # set to 1 if we find any duplicate sequence names
  @{$seq_name_AR} = ();
  %{$seq_len_HR} = ();

  # parse the seqstat -a output 
  # sequences must have non-empty names (else esl-seqstat call would have failed)
  # lengths must be >= 0 (lengths of 0 are okay)
  while(my $line = <IN>) { 
    # = lcl|dna_BP331_0.3k:467     1232 
    # = lcl|dna_BP331_0.3k:10     1397 
    # = lcl|dna_BP331_0.3k:1052     1414 
    chomp $line;
    #print $line . "\n";
    if($line =~ /^\=\s+(\S+)\s+(\d+)/) { 
      ($seq_name, $length) = ($1, $2);

      if(exists($seq_exists_H{$seq_name})) { 
        if(exists($seq_dups_H{$seq_name})) { 
          $seq_dups_H{$seq_name}++; 
        }
        else { 
          $seq_dups_H{$seq_name} = 2;
        }
        $at_least_one_dup = 1;
      }
      else { 
        $seq_exists_H{$seq_name} = 1; 
      }

      push(@{$seq_name_AR}, $seq_name);
      $seq_len_HR->{$seq_name} = $length;
      $tot_length += $length;
      $nread++;
    }
  }
  close(IN);
  if($nread == 0) { 
    ofile_FAIL("ERROR in $sub_name, did not read any sequence lengths in esl-seqstat file $seqstat_file, did you use -a option with esl-seqstat", "RIBO", 1, $FH_HR);
  }
  if($at_least_one_dup) { 
    my $i = 1;
    my $die_string = "\nERROR, not all sequences in input sequence file have a unique name. They must.\nList of sequences that occur more than once, with number of occurrences:\n";
    foreach $seq_name (sort keys %seq_dups_H) { 
      $die_string .= "\t($i) $seq_name $seq_dups_H{$seq_name}\n";
      $i++;
    }
    $die_string .= "\n";
    ofile_FAIL($die_string, undef, 1, $FH_HR);
  }

  return $tot_length;
}

#################################################################
# Subroutine : seq_FastaRemoveDescriptions()
# Incept:      EPN, Mon Mar 25 11:30:16 2019
#
# Purpose:     Given an FASTA file, create a new one that is
#              identical but with sequence descriptions removed.
#              DOES NOT VALIDATE INPUT FILE IS IN PROPER FASTA FORMAT.
#              
# Arguments: 
#   $in_file:        input FASTA file
#   $out_file:       output FASTA file, with sequence descriptions removed
#   $ofile_info_HHR: ref to the ofile info 2D hash, can be undef
# 
# Returns:     void
#
# Dies:        If unable to open $in_file for reading or $out_file
#              for writing.
#
################################################################# 
sub seq_FastaRemoveDescriptions {
  my $nargs_expected = 3;
  my $sub_name = "seq_FastaRemoveDescriptions()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($in_file, $out_file, $ofile_info_HHR) = (@_);

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef; # for convenience

  open(IN,       $in_file)  || ofile_FileOpenFailure($in_file,  undef, $sub_name, $!, "reading", $FH_HR);
  open(OUT, ">", $out_file) || ofile_FileOpenFailure($out_file, undef, $sub_name, $!, "writing", $FH_HR);

  while(my $line = <IN>) { 
    chomp $line;
    if($line =~ /^>(\S+)/) { 
      print OUT ">" . $1 . "\n";
    }
    else { 
      print OUT $line . "\n";
    }
  }
  close(IN);
  close(OUT);

  return;
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
