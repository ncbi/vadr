#!/usr/bin/perl
# 
# version: 1.0.4 [March 2020]
#
# vadr-fast.pm
# Eric Nawrocki
# EPN, Fri Mar 27 10:28:01 2020
# 
# Perl module used by v-annotate.pl script related to 
# the --fast option to use blastn for the classification,
# and coverage determination stages, and to accelerate
# the alignment stage. 
#
#########################
# Naming conventions used in this file:
#
# - Hash data structures use the abbreviation 'H', typically at the end of 
#   a variable name (e.g. %length_H)
#
# - Array data structures use the abbreviation 'A', typically at the end of 
#   a variable name (e.g. @accession_A)
#
# - Multidimensional arrays/hashes use concatenated 'H's and 'A's, e.g.:
#   @size_AA is a two dimensional array, and %info_HHA is a two dimensional 
#   hash of arrays
# 
# - File handles typically start or end with 'FH', e.g. $log_FH, or a hash
#   of file handles that is used commonly is %FH_HR.
# 
#########################
# Common data structures used in this file:
#
# - $ftr_info_AHR: reference to an array of hashes with feature information for a single model.
#                   
# - $sgm_info_AHR: reference to an array of hashes with model-segment information for a single model.
#                   
# - $alt_info_HHR: reference to a hash of hashes with alert information.
#
########################################################################################
#
# List of subroutines in this file, divided into categories. 
#
use strict;
use warnings;
use Cwd;
use LWP::Simple; 

require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seq.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

###########################################################################
# the next line is critical, a perl module must return a true value
return 1;
###########################################################################

#################################################################
# Subroutine:  blastn_run()
# Incept:      EPN, Fri Mar 27 11:11:24 2020
#
# Purpose:     Run blastn locally.
#
# Arguments: 
#  $execs_HR:        ref to executables with "esl-ssplit" and "cmsearch"
#                    defined as keys
#  $db_file:         name of blast db file to use
#  $seq_file:        name of sequence file with all sequences to run against
#  $out_root:        string for naming output files
#  $nseq:            number of sequences in $seq_file
#  $progress_w:      width for outputProgressPrior output
#  $opt_HHR:         REF to 2D hash of option values, see top of sqp-opts.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
#
# Returns:     void
# 
# Dies: If blastn executable doesn't exist or command fails
################################################################# 
sub blastn_run { 
  my $sub_name = "blastn_run";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $db_file, $seq_file, $out_root, 
      $nseq, $progress_w, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"};
  my $log_FH = $FH_HR->{"log"}; # for convenience

  my $do_keep = opt_Get("--keep", $opt_HHR);

  my $start_secs = ofile_OutputProgressPrior(sprintf("Classifying sequences with blastn ($nseq seq%s)", ($nseq > 1) ? "s" : ""), $progress_w, $log_FH, *STDOUT);

  my $blastn_out_file = $out_root . ".r1.blastn.out";
  my $opt_str = "-num_alignments 10 -query $seq_file -db $db_file -out $blastn_out_file";
  my $blastn_cmd = $execs_HR->{"blastn"} . " $opt_str";
  
  utl_RunCommand($blastn_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  return;
}
