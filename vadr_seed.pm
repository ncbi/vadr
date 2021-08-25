#!/usr/bin/env perl
# 
# version: 1.3 [Aug 2021]
#
# vadr_seed.pm
# Eric Nawrocki
# EPN, Fri Mar 27 10:28:01 2020
# 
# Perl module used by v-annotate.pl script related to 
# the -s option to use blastn for the classification,
# and coverage determination stages, and to accelerate
# the alignment stage. Alignment acceleration achieved
# by using the maximum length ungapped region in the top
# as a fixed alignment and only aligning the flanking terminal
# 5' and/or 3' regions with cmalign, if necessary.
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
## 
#
# See vadr.pm for additional notes
#
use strict;
use warnings;

require "vadr.pm"; 
require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seq.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

#########################################################################################
#
# List of subroutines in this file, divided into categories. 
#
# Subroutines related to running and parsing blastn:
# run_blastn_and_summarize_output()
# parse_blastn_results()
# blastn_pretblout_to_tblout()
# parse_blastn_indel_strings()
# parse_blastn_indel_token()
# parse_blastn_indel_file_to_get_subseq_info()
# 
# Subroutines related to joining alignments:
# join_alignments_and_add_unjoinbl_alerts()
# join_alignments_helper()
# update_overflow_info_for_joined_alignments
# 
#########################################################################################
# 
# List of subroutines in this file
# 
#################################################################
# Subroutine:  run_blastn_and_summarize_output()
# Incept:      EPN, Fri Mar 27 11:11:24 2020
#
# Purpose:     Run all input sequences as queries against a 
#              blastn db of the model sequences and summarize
#              the output with parse_blast.pl.
#
# Arguments: 
#  $execs_HR:        ref to executables with "esl-ssplit" and "cmsearch"
#                    defined as keys
#  $db_file:         name of blast db file to use
#  $seq_file:        name of sequence file with all sequences to run against
#  $out_root:        string for naming output files
#  $stg_key:         stage key, "rpn.cls" or "std.cls"
#  $nseq:            number of sequences in $seq_file
#  $progress_w:      width for outputProgressPrior output
#  $opt_HHR:         REF to 2D hash of option values, see top of sqp-opts.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
#
# Returns:     void
# 
# Dies: If blastn executable doesn't exist or command fails
################################################################# 
sub run_blastn_and_summarize_output { 
  my $sub_name = "run_blastn_and_summarize_output";
  my $nargs_expected = 9;

  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $db_file, $seq_file, $out_root, $stg_key,
      $nseq, $progress_w, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"};
  my $log_FH = $FH_HR->{"log"}; # for convenience
  my $ncpu = opt_Get("--cpu", $opt_HHR);
  if($ncpu == 0) { $ncpu = 1; }

  if(($stg_key ne "rpn.cls") && ($stg_key ne "std.cls")) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized stage key: $stg_key, should be rpn.cls or std.cls", 1, $FH_HR);
  }

  my $do_keep = opt_Get("--keep", $opt_HHR);
  my $stg_desc = ($stg_key eq "rpn.cls") ? 
      sprintf("Preprocessing for N replacement: blastn classification ($nseq seq%s)", (($nseq > 1) ? "s" : "")) :
      sprintf("Classifying sequences with blastn ($nseq seq%s)", (($nseq > 1) ? "s" : ""));
  my $start_secs = ofile_OutputProgressPrior($stg_desc, $progress_w, $log_FH, *STDOUT);
  my $blastn_out_file = $out_root . ".$stg_key.blastn.out";
  my $opt_str = "-num_threads $ncpu -query $seq_file -db $db_file -out $blastn_out_file -word_size " . opt_Get("--s_blastnws", $opt_HHR); 
  my $blastn_cmd = $execs_HR->{"blastn"} . " $opt_str";
  
  utl_RunCommand($blastn_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "$stg_key.blastn.out", $blastn_out_file, 0, $do_keep, "blastn output");

  # now summarize its output
  # use --splus to skip alignment parsing of hits to negative strand of subject
  # that's important because parse_blast.pl doesn't have the code to deal with
  # such alignments (it was written for blastx originally for which subject is
  # always +) but that's okay because those hits are to negative strand of
  # the sequence (actually negative strand of the subject/model but blastn revcomps
  # the subject instead of the query like cmsearch/cmscan would do), and we don't care
  # about negative strand hit indel info.
  my $blastn_summary_file = $out_root . ".$stg_key.blastn.summary.txt";
  my $parse_cmd = $execs_HR->{"parse_blast"} . " --program n --input $blastn_out_file --splus > $blastn_summary_file";
  utl_RunCommand($parse_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "$stg_key.blastn.summary", $blastn_summary_file, 0, $do_keep, "parsed (summarized) blastn output");

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  return;
}

#################################################################
# Subroutine:  parse_blastn_results()
# Incept:      EPN, Thu Oct  4 15:25:00 2018
#              [based on parse_blastx_results() which was 
#               modified from parse_blastx_results() subroutine
#               Alejandro Schaffer in compare_predictions.pl]
#
# Purpose:    Parse blastn summary file and output a summary file
#             and create one of two sets of output files:
#
#             Output mode 1: 1 file, produced if $seq2mdl_HR is undef
#
#             "blastn.{rpn.cls,std.cls}.pretblout" file: cmsearch --trmF3 output
#             format file with each hit on a separate line and
#             individual hit scores reported for each hit. This is
#             later processed by blastn_pretblout_to_tblout() to sum
#             bit scores for all hits with the same model/seq/strand
#             so we can classify sequences the same way we do in
#             default mode (cmsearch based classification).
#
#             Output mode 2: 2 files produced per model with >= 1 matching
#             sequence, produced if $seq2mdl_HR is defined.
#
#             "search.{rpn.cls,std.cls,rpn.cdt,std.cdt}.<mdlname>.tblout":
#             cmsearch --tblout format file with each hit for a
#             sequence on + strand that is classified to model
#             <mdlname>.
#
#             "blastn.{rpn.cls,std.cls,rpn.cdt,std.cdt}.<mdlname>.indel.txt":
#             one line per sequence with all inserts and deletes in
#             all blastn hit alignments for each sequence that is
#             classified to <mdlname> on strand +.
#
# Arguments: 
#  $blastn_summary_file: path to blastn summary file to parse
#  $seq_len_HR:          REF to hash of sequence lengths
#  $seq2mdl_HR:          REF to hash mapping each sequence to the model
#                        it is classified to, if undef serves as flag
#                        that we will create set 1 of output files
#  $mdl_name_AR:         REF to array of model names that are keys in
#                        %{$seq2mdl_HR}, can be undef if $seq2mdl_HR is undef
#  $out_root:            output root for the file names
#  $stg_key:             stage key, "rpn.cls" or "std.cls" for classification
#                        or "rpn.cdt" or "std.cdt" for coverage determination
#  $opt_HHR:             REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:      REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If blastn fails.
#
################################################################# 
sub parse_blastn_results { 
  my $sub_name = "parse_blastn_results";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastn_summary_file, $seq_len_HR, $seq2mdl_HR, $mdl_name_AR, 
      $out_root, $stg_key, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  if(($stg_key ne "rpn.cls") && ($stg_key ne "rpn.cdt") && 
     ($stg_key ne "std.cls") && ($stg_key ne "std.cdt")) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized stage key: $stg_key, should be rpn.cls, rpn.cdt, std.cls, or std.cdt", 1, $FH_HR);
  }
  if((($stg_key eq "rpn.cdt") || ($stg_key eq "std.cdt")) && 
     (! defined $seq2mdl_HR)) { 
    ofile_FAIL("ERROR in $sub_name, stage key is $stg_key but seq2mdl_HR is undef", 1, $FH_HR);
  }

  my $pretblout_FH = undef; # defined if output mode 1 (if ! defined $seq2mdl_HR)
  my %tblout_FH_H  = ();    # defined if output mode 2 (if   defined $seq2mdl_HR)
  my %indel_FH_H   = ();    # defined if output mode 2 (if   defined $seq2mdl_HR)
  my $outfile_key  = undef; # a key for an output file in %{$ofile_info_HHR}
  my $small_value  = 0.000001;
  my $min_bitsc    = opt_Get("--s_blastnsc", $opt_HHR) - $small_value;
  my $do_keep      = opt_Get("--keep", $opt_HHR) ? 1 : 0;
  my $mdl_name = undef;
  if(! defined $seq2mdl_HR) { 
    # output mode 1, open the pretblout output file 
    # 
    # We write to this file as we parse the $blastn_summary_file but we
    # have to post-process it in blastn_pretblout_to_tblout so that the
    # top hit per model/sequence/strand trio includes the *summed* score
    # for that model/strand/strand instead of just the hit score. This
    # way we will match the cmsearch --trmF3 output downstream steps
    # expect.
    ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, "$stg_key.blastn.pretblout", $out_root . ".$stg_key.blastn.pretblout",  0, $do_keep, "blastn output converted to cmsearch --trmF3 tblout format (hit scores)");
    $pretblout_FH = $ofile_info_HHR->{"FH"}{"$stg_key.blastn.pretblout"}; 
    printf $pretblout_FH ("%-30s  %-30s  %8s  %9s  %9s  %6s  %6s  %3s  %11s\n", 
                          "#modelname/subject", "sequence/query", "bitscore", "start", "end", "strand", "bounds", "ovp", "seqlen");
  }
  else { # $seq2mdl_HR is defined
    # output mode 2, for each model in @{$mdl_name_AR}, open the 
    # coverage determination tblout files in cmsearch --tblout 
    # format (not --trmF3 output format) and the indel files 
    foreach $mdl_name (@{$mdl_name_AR}) { 
      $outfile_key = "$stg_key.$mdl_name.tblout";
      ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, $outfile_key, $out_root . "." . $outfile_key,  0, $do_keep, "blastn output converted to cmsearch tblout format for model $mdl_name");
      $tblout_FH_H{$mdl_name} = $ofile_info_HHR->{"FH"}{$outfile_key};

      $outfile_key = "$stg_key.$mdl_name.indel";
      ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, $outfile_key, $out_root . "." . $outfile_key,  0, $do_keep, "blastn indel information for model $mdl_name");
      $indel_FH_H{$mdl_name} = $ofile_info_HHR->{"FH"}{$outfile_key};
    }
  }

  # open and parse input blastn summary file
  utl_FileValidateExistsAndNonEmpty($blastn_summary_file, "blastn summary file", $sub_name, 1, $FH_HR);
  open(IN, $blastn_summary_file) || ofile_FileOpenFailure($blastn_summary_file, $sub_name, $!, "reading", $FH_HR);
  my $line_idx   = 0;
  my $seq_name   = undef; # sequence name this hit corresponds to (query)
  my $seq_len    = undef; # length of query sequence 
  my $cur_FH     = undef; # current file handle 
  my %cur_H = ();     # values for current hit (HSP)
  my %scsum_HHH = (); # 3D hash with summed scores for model/seq/strand trios:
                      # key 1: model/subject
                      # key 2: sequence/query
                      # key 3: strand ("+" or "-")
                      # value: summed bit score for all hits for this model/sequence/strand trio
  
  # 
  # Order of lines in <IN>:
  # -----per-query/target-block---
  # QACC
  # QDEF   ignored
  # QLEN   
  # MATCH  ignored
  # HACC  
  # HDEF   ignored
  # SLEN   
  # ------per-HSP-block------
  # HSP   
  # BITSCORE 
  # RAWSCORE ignored
  # EVALUE 
  # HLEN   ignored
  # IDENT  ignored
  # GAPS   ignored
  # STRAND
  # STOP   
  # DEL    
  # MAXDE  
  # INS    
  # MAXINS 
  # QRANGE
  # SRANGE 
  # ------per-HSP-block------
  # END_MATCH

  while(my $line = <IN>) { 
    chomp $line;
    $line_idx++;
    if($line ne "END_MATCH") { 
      my @el_A = split(/\t/, $line);
      if(scalar(@el_A) != 2) { 
        ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, did not read exactly 2 tab-delimited tokens in line $line", 1, $FH_HR);
      }
      my ($key, $value) = (@el_A);
      if($key eq "QACC") { 
        $cur_H{$key} = $value;
        if(! defined $seq_len_HR->{$value}) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, unexpected sequence name $value in line $line", 1, $FH_HR);
        }
      }
      elsif($key eq "QLEN") { 
        if(! defined $cur_H{"QACC"}) { 
          ofile_FAIL("ERROR in $sub_name, ing $blastn_summary_file, read QLEN line before QACC line (line: $line_idx)", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
        if($value !~ /^\d+$/) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, unable to parse query length in line $line", 1, $FH_HR);
        }
        if($value != $seq_len_HR->{$cur_H{"QACC"}}) {
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read query length of $value from QLEN line for query " . $cur_H{"QACC"} . ", but expected " . $seq_len_HR->{$cur_H{"QACC"}} . " from seq_len_HR, on line: $line", 1, $FH_HR);
        }
      }
      elsif($key eq "HACC") { 
        if(! defined $cur_H{"QACC"}) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read HACC line $line before QACC line (line: $line_idx)\n", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
      }
      elsif($key eq "HLEN") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read HLEN line $line before one or both of QACC and HACC lines (line: $line_idx)\n", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
        if($value !~ /^\d+$/) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, unable to parse subject length in line $line", 1, $FH_HR);
        }
      }
      elsif($key eq "HSP") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read HSP line before one or both of QACC and HACC lines (line: $line_idx)\n", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
        if($value !~ /^(\d+)$/) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, unable to parse blastn summary HSP line $line", 1, $FH_HR);
        }
      }
      elsif($key eq "BITSCORE") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"}) || (! defined $cur_H{"HSP"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read BITSCORE line before one or more of QACC, HACC, or HSP lines (line: $line_idx)\n", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
        if(($value !~ /^\d+\.\d+$/) && ($value !~ /^\d+$/) && ($value !~ /^\d+\.\d+e[+-]\d+/)) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, unable to parse blastn summary BITSCORE line $line", 1, $FH_HR);
        }
      }
      elsif($key eq "EVALUE") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"}) || (! defined $cur_H{"HSP"}) || (! defined $cur_H{"BITSCORE"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read EVALUE line before one or more of QACC, HACC, HSP, or BITSCORE lines (line: $line_idx)\n", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
      }
      elsif($key =~ m/^[SQ]STRAND$/) { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"}) || (! defined $cur_H{"HSP"}) || (! defined $cur_H{"BITSCORE"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read QSTRAND OR SSTRAND line before one or more of QACC, HACC, HSP, or BITSCORE lines (line: $line_idx)\n", 1, $FH_HR);
        }
        if($value =~ /^([\+\-])$/) { 
          $cur_H{$key} = $1;
          # all query strands should be + (unless blastn is reverse complementing queries and running those too)
          if(($key eq "QSTRAND") && ($cur_H{$key} ne "+")) { 
            ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, query strand not + line $line", 1, $FH_HR);
          }            
        }
        else { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, unable to parse blastn summary $key line $line ($key $value)", 1, $FH_HR);
        }
      }
      elsif(($key eq "STOP") || ($key eq "DEL") || ($key eq "INS")) { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"}) || (! defined $cur_H{"HSP"}) || (! defined $cur_H{"BITSCORE"}) || (! defined $cur_H{"QSTRAND"}) || (! defined $cur_H{"SSTRAND"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read $key line before one or more of QACC, HACC, HSP, BITSCORE or FRAME lines (line: $line_idx)\n", 1, $FH_HR);
        }
        if(($value ne "") && ($value ne "BLASTNULL")) { 
          $cur_H{$key} = $value;
        } 
      }
      elsif($key eq "SRANGE") { 
        # we don't require all of QACC, HACC, HSP, BITSCORE and STRAND even though we should have them
        # sometimes we don't (may be a bug in parse-blastn.pl), we only require QACC
        if(! defined $cur_H{"QACC"}) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read $key line before QACC line (line: $line_idx)\n", 1, $FH_HR);
        }
        if($value eq "..") { # special case, no hits, silently move on
          ;
        }
        elsif(! defined $cur_H{"BITSCORE"}) { # special case, no BITSCORE lines yet seen (may be a bug in parse-blastn.pl?), silently move on
          ;
        }
        elsif($value =~ /^(\d+)..(\d+)$/) { 
          ($cur_H{"SRANGESTART"}, $cur_H{"SRANGESTOP"}) = ($1, $2);
        }
        else { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, unable to parse blastn summary QRANGE line $line", 1, $FH_HR);
        }
      }
      elsif($key eq "QRANGE") { 
        # we don't require all of QACC, HACC, HSP, BITSCORE and STRAND even though we should have them
        # sometimes we don't (may be a bug in parse-blastn.pl), we only require QACC
        if(! defined $cur_H{"QACC"}) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read $key line before QACC line (line: $line_idx)\n", 1, $FH_HR);
        }
        if($value eq "..") { # special case, no hits, silently move on
          ;
        }
        elsif(! defined $cur_H{"BITSCORE"}) { # special case, no BITSCORE lines yet seen (may be a bug in parse-blastn.pl?), silently move on
          ;
        }
        elsif(($value =~ /^(\d+)..(\d+)$/) && ($cur_H{"BITSCORE"} >= $min_bitsc)) { 
          ($cur_H{"QRANGESTART"}, $cur_H{"QRANGESTOP"}) = ($1, $2);

          # output data in cmscan --trmF3 format
          if((! defined $cur_H{"SRANGESTART"}) || (! defined $cur_H{"SRANGESTOP"})) { 
            ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read $key line before SRANGE line (line: $line_idx)\n", 1, $FH_HR);
          }
          my $cur_seq_name  = $cur_H{"QACC"};
          my $cur_seq_len   = $seq_len_HR->{$cur_seq_name}; # we checked this is defined when we read QACC line
          my $cur_mdl_name  = $cur_H{"HACC"};
          my $cur_bit_score = $cur_H{"BITSCORE"}; 
          my $cur_evalue    = $cur_H{"EVALUE"}; 
          my ($cur_seq_start, $cur_seq_stop, $cur_seq_strand, $cur_mdl_start, $cur_mdl_stop) = (undef, undef, undef, undef, undef);
          # if subject strand is negative (and subject is the model)
          # we need to revcomp both the query and subject coords, to match cmscan behavior
          # which is the inverse of blastn, whereas cmscan revcomps the sequence (blastn equivalent of query)
          # and searches it against the always positive model (blastn equivalent of the subject), blastn
          # revcomps the subject and searches it against the always positive query
          if($cur_H{"SSTRAND"} eq "+") { 
            $cur_seq_start  = $cur_H{"QRANGESTART"};
            $cur_seq_stop   = $cur_H{"QRANGESTOP"};
            $cur_seq_strand = "+";
            $cur_mdl_start  = $cur_H{"SRANGESTART"};
            $cur_mdl_stop   = $cur_H{"SRANGESTOP"};
          }
          else { # SSTRAND is -
            $cur_seq_start  = $cur_H{"QRANGESTOP"};
            $cur_seq_stop   = $cur_H{"QRANGESTART"};
            $cur_seq_strand = "-";
            $cur_mdl_start  = $cur_H{"SRANGESTOP"};
            $cur_mdl_stop   = $cur_H{"SRANGESTART"};
          }
          if(defined $pretblout_FH) { 
            printf $pretblout_FH ("%-30s  %-30s  %8.1f  %9d  %9d  %6s  %6s  %3s  %11s\n", 
                                  $cur_mdl_name, 
                                  $cur_seq_name,
                                  $cur_bit_score,
                                  $cur_seq_start, 
                                  $cur_seq_stop, 
                                  $cur_seq_strand, 
                                  sprintf("    %s%s", 
                                          ($cur_seq_start == 1)            ? "[" : ".",
                                          ($cur_seq_stop  == $cur_seq_len) ? "]" : "."),
                                  "?",
                                  $cur_seq_len);
            
            # update summed score in %scsum_HHH for this model/seq/strand trio
            if(! defined $scsum_HHH{$cur_mdl_name}) { 
              %{$scsum_HHH{$cur_mdl_name}} = ();
            }
            if(! defined $scsum_HHH{$cur_mdl_name}{$cur_seq_name}) { 
              %{$scsum_HHH{$cur_mdl_name}{$cur_seq_name}} = ();
            }
            if(! defined $scsum_HHH{$cur_mdl_name}{$cur_seq_name}{$cur_seq_strand}) { 
              $scsum_HHH{$cur_mdl_name}{$cur_seq_name}{$cur_seq_strand} = 0.;
            }
            $scsum_HHH{$cur_mdl_name}{$cur_seq_name}{$cur_seq_strand} += $cur_H{"BITSCORE"};
          }            
          else { 
            # in output mode 2, but we only output if this hit is 
            # a <seq>/<model> pair where <seq> is classified to <model>
            if((defined $seq2mdl_HR->{$cur_seq_name}) &&  
               ($seq2mdl_HR->{$cur_seq_name} eq $cur_mdl_name)) { 
              #target name         accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
              if(! defined $tblout_FH_H{$cur_mdl_name}) { 
                ofile_FAIL("ERROR in $sub_name, read unexpected model name: $cur_mdl_name", 1, $FH_HR);
              }
              $cur_FH = $tblout_FH_H{$cur_mdl_name};
              printf $cur_FH ("%-s  -  %-s  -  blastn  %d  %d  %d  %d  %s  -  -  -  0.0  %8.1f  %s  ?  -\n", 
                              $cur_seq_name,
                              $cur_mdl_name,
                              $cur_mdl_start,
                              $cur_mdl_stop,
                              $cur_seq_start, 
                              $cur_seq_stop, 
                              $cur_seq_strand, 
                              $cur_bit_score, 
                              $cur_evalue);
              # and output indel info to a separate file
              $cur_FH = $indel_FH_H{$cur_mdl_name};
              printf $cur_FH ("%s  %s  %s  %s  %s  %s  %s  %s\n",
                              $cur_mdl_name,
                              $cur_seq_name,
                              vdr_CoordsSegmentCreate($cur_mdl_start, $cur_mdl_stop, "+", $FH_HR),
                              $cur_H{"HLEN"},
                              vdr_CoordsSegmentCreate($cur_seq_start, $cur_seq_stop, $cur_seq_strand, $FH_HR),
                              $cur_seq_len,
                              (defined $cur_H{"INS"}) ? $cur_H{"INS"} : "BLASTNULL",
                              (defined $cur_H{"DEL"}) ? $cur_H{"DEL"} : "BLASTNULL");
            }
          }
        } # end of 'else' entered if SRANGE is NOT ".."
        # reset variables
        my $save_qacc = $cur_H{"QACC"}; 
        my $save_hacc = $cur_H{"HACC"};
        %cur_H = (); 
        $cur_H{"QACC"} = $save_qacc; # save QACC
        $cur_H{"HACC"} = $save_hacc; # save HACC
      } # end of 'elsif $key eq "SRANGE"
    } # end of 'if($line ne "END_MATCH")'
    else { # $line eq "END_MATCH"
      # reset variables
      my $save_qacc = $cur_H{"QACC"};
      %cur_H = (); 
      $cur_H{"QACC"} = $save_qacc; # save QACC
      # we don't update seq_name seq_len until we see a new QACC
    }
  } # end of 'while($my $line = <IN>)'
  close(IN);

  # close files depending on output mode
  if(defined $pretblout_FH) { 
    # output mode 1
    # close the pretblout file
    # then call convert it pretblout to tblout (cmscan --trmF3 format)
    # which will have scores summed for each seq/mdl/strand trio
    close $ofile_info_HHR->{"FH"}{"$stg_key.blastn.pretblout"};
    blastn_pretblout_to_tblout($ofile_info_HHR->{"fullpath"}{"$stg_key.blastn.pretblout"}, 
                               \%scsum_HHH, $out_root, $stg_key, $opt_HHR, $ofile_info_HHR);
  }
  else { 
    # output mode 2:
    # close the per model files:
    foreach $mdl_name (@{$mdl_name_AR}) { 
      $outfile_key = "$stg_key.$mdl_name.tblout";
      close $ofile_info_HHR->{"FH"}{$outfile_key};
      $outfile_key = "$stg_key.$mdl_name.indel";
      close $ofile_info_HHR->{"FH"}{$outfile_key};
    }
  }
  return;
}

#################################################################
# Subroutine:  blastn_pretblout_to_tblout()
# Incept:      EPN, Sat Mar 28 14:04:18 2020
#
# Purpose:     Given a 'pretblout' file created by parse_blastn_results
#              generate the 'tblout' file that is identical to the
#              'pretblout' except that scores are summed over all hits
#              for each sequence/model/strand trio for the top hit for
#              each sequence/model/strand trio and 0 for all other hits
#              for that trio. This unusual step is necessary so the 
#              highest scoring hit in the tblout file can be used to 
#              define the winning model (as it is for the cmscan --trmF3
#              tblout file).
#
# Arguments: 
#  $blastn_pretblout_file: path to blastn summary file to parse
#  $scsum_HHHR:            ref to 3D hash, 
#                          key 1: model/subject
#                          key 2: sequence/query
#                          key 3: strand ("+" or "-")
#                          value: summed bit score for all hits for this model/sequence/strand trio
#                          NOTE: all values in this 3D hash are set to 0. by this subroutine!
#  $out_root:              output root for the file names
#  $stg_key:               stage key, "rpn.cls" or "std.cls" for classification (cmscan) ,
#                          or "rpn.cdt" or "std.cdt" for coverage determination (cmsearch)
#  $opt_HHR:               REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:        REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       if unable to parse a pretblout file
#
################################################################# 
sub blastn_pretblout_to_tblout { 
  my $sub_name = "blastn_pretblout_to_tblout";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastn_pretblout_file, $scsum_HHHR, $out_root, $stg_key, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $do_keep = opt_Get("--keep", $opt_HHR);

  if(($stg_key ne "rpn.cls") && ($stg_key ne "rpn.cdt") && 
     ($stg_key ne "std.cls") && ($stg_key ne "std.cdt")) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized stage key: $stg_key, should be rpn.cls, rpn.cdt, std.cls, or std.cdt", 1, $FH_HR);
  }

  ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, "$stg_key.tblout", $out_root . ".$stg_key.tblout",  0, $do_keep, "blastn output converted to cmscan --trmF3 tblout format (summed hit scores)");
  my $tblout_FH = $FH_HR->{"$stg_key.tblout"}; 

  # open and parse input blastn summary file
  utl_FileValidateExistsAndNonEmpty($blastn_pretblout_file, "blastn pretblout file", $sub_name, 1, $FH_HR);
  open(IN, $blastn_pretblout_file) || ofile_FileOpenFailure($blastn_pretblout_file, $sub_name, $!, "reading", $FH_HR);

  while(my $line = <IN>) { 
    chomp $line;
    if($line !~ m/^#/) { 
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) != 9) { 
        ofile_FAIL("ERROR in $sub_name, unable to parse $blastn_pretblout_file line (wrong number of space-delimited tokens):\n$line\n", 1, $FH_HR);
      }
      my ($model, $seq, $bitsc, $start, $end, $strand, $bounds, $ovp, $seqlen) = (@el_A);
      if(! defined $scsum_HHHR->{$model}{$seq}{$strand}) { 
        ofile_FAIL("ERROR in $sub_name, read model/seq/strand trio not in input scsum hash: model:$model, seq:$seq, strand:$strand on line:\n$line\n", 1, $FH_HR);
      }
      printf $tblout_FH ("%-30s  %-30s  %8.1f  %9d  %9d  %6s  %6s  %3s  %11s\n", 
                         $seq, $model, $scsum_HHHR->{$model}{$seq}{$strand}, 
                         $start, $end, $strand, $bounds, $ovp, $seqlen);
      # set scsum to zero for all subsequent hits to this trio      
      $scsum_HHHR->{$model}{$seq}{$strand} = 0.; 
    }
    else { 
      print $tblout_FH ($line . "\n");
    }
  }
  close(IN);
  close $ofile_info_HHR->{"FH"}{"$stg_key.tblout"};
}

#################################################################
# Subroutine:  parse_blastn_indel_strings()
# Incept:      EPN, Mon Mar 30 06:53:35 2020
#
# Purpose:     Given information on where insertions and deletions
#              are in a blastn alignment, determine the
#              max ungapped string and optional return the alignment
#              of the sequence as a string.
#
# Arguments: 
#  $in_mdl_coords_sgm: vadr coords segment indicating model boundaries of alignment
#  $in_seq_coords_sgm: vadr coords segment indicating sequence boundaries of alignment
#  $blastn_ins_str:    insert string from parse_blast.pl run on blastn output indicating
#                      where insertions in query (sequence) with respect to model (subject)
#                         are, in format:
#                         <instok1>...<instokn> (note: final token does not have ';' after it)
#                         <instok> = Q<d1>:S<d2>+<d3>
#                         <d1> is query   (sequence) position after which insertion of length <d3> nt exists
#                         <d2> is subject (model)    position after which insertion of length <d3> nt exists
#  $blastn_del_str:     delete string from parse_blast.pl indicating where deletions are 
#                       where insertions are in format:
#                         <deltok1>;<deltok2>;...<deltokn> (note: final token does not have ';' after it)
#                         <deltok> = Q<d1>:S<d2>+<d3>
#                         <d1> is query   (sequence) position after which insertion of length <d3> nt exists
#                         <d2> is subject (model)    position after which insertion of length <d3> nt exists
#  $ofile_info_HHR:      REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    Two values:
#             $ugp_mdl_coords: vadr coords string with all segments of ungapped alignments between
#                              query and subject in model (subject) coords, between model coords
#                              in $in_mdl_coords_sgm.
#             $ugp_seq_coords: vadr coords string with all segments of ungapped alignments between
#                              query and subject in sequence (query) coords, between sequence coords
#                              in $in_seq_coords_sgm.
#
# Dies:       if unable to parse the indel strings
#
################################################################# 
sub parse_blastn_indel_strings { 
  my $sub_name = "parse_blastn_indel_strings";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($in_mdl_coords_sgm, $in_seq_coords_sgm, $blastn_ins_str, $blastn_del_str, $ofile_info_HHR) = @_;
  
  my $FH_HR = (defined $ofile_info_HHR) ? $ofile_info_HHR->{"FH"} : undef; # undef helpful for tests

  my @ins_tok_A = ();
  if($blastn_ins_str ne "BLASTNULL") { 
    @ins_tok_A = split(";", $blastn_ins_str);
  }
  my $nins = scalar(@ins_tok_A);
  
  my @del_tok_A = ();
  if($blastn_del_str ne "BLASTNULL") { 
    @del_tok_A = split(";", $blastn_del_str);
  }
  my $ndel = scalar(@del_tok_A);
  # printf("nins: $nins ndel: $ndel\n");
  
  my ($mdl_start, $mdl_stop, $mdl_strand) = vdr_CoordsSegmentParse($in_mdl_coords_sgm, $FH_HR);
  my ($seq_start, $seq_stop, $seq_strand) = vdr_CoordsSegmentParse($in_seq_coords_sgm, $FH_HR);

  if($mdl_strand ne "+") { ofile_FAIL("ERROR in $sub_name, model (subject) strand in model coords: $in_mdl_coords_sgm is not +", 1, $FH_HR); }
  if($seq_strand ne "+") { ofile_FAIL("ERROR in $sub_name, sequence (query) strand in model coords: $in_seq_coords_sgm is not +", 1, $FH_HR); }
  
  # construct a vadr coords string of ungapped aligned segments, segment by segment
  my $ugp_mdl_coords = "";
  my $ugp_seq_coords = "";
  my $cur_mdl_start = $mdl_start;
  my $cur_seq_start = $seq_start;

  my $ii = 0; # insert index in @ins_tok_A
  my $di = 0; # delete index in @del_tok_A
  my ($ins_seq_pos, $ins_mdl_pos, $ins_len) = (undef, undef, undef);
  my ($del_seq_pos, $del_mdl_pos, $del_len) = (undef, undef, undef);
  my $update_ins = 1; # flag set when we should update our insert token
  my $update_del = 1; # flag set when we should update our delete token
  my $ins_tok = undef;
  my $del_tok = undef;
  while(($ii < $nins) || ($di < $ndel)) { 
    my $ins_tok = ($ii < $nins) ? $ins_tok_A[$ii] : "NULL";
    my $del_tok = ($di < $ndel) ? $del_tok_A[$di] : "NULL";
    # printf("\tin loop ii: $ii ins_tok: $ins_tok, di: $di del_tok: $del_tok\n");
    if($update_ins) { 
      if($ii < $nins) { 
        $ins_tok = $ins_tok_A[$ii];
        ($ins_mdl_pos, $ins_seq_pos, $ins_len) = parse_blastn_indel_token($ins_tok, "insert", $FH_HR);
      }
      else {
        ($ins_tok, $ins_mdl_pos, $ins_seq_pos, $ins_len) = (undef, undef, undef, undef);
      }
    }
    if($update_del) {
      if($di < $ndel) { 
        $del_tok = $del_tok_A[$di];
        ($del_mdl_pos, $del_seq_pos, $del_len) = parse_blastn_indel_token($del_tok, "delete", $FH_HR);
      }
      else {
        ($del_tok, $del_mdl_pos, $del_seq_pos, $del_len) = (undef, undef, undef, undef);
      }
    }
    $update_ins = 0; # set to 1 below, if our next indel is an insert ($next_is_insert set to 1) and we need to read the next one
    $update_del = 0; # set to 1 below, if our next indel is a  delete ($next_is_delete set to 1) and we need to read the next one
    my $next_is_insert = 0; # set to 1 below if we determine next token is insert, not delete
    my $next_is_delete = 0; # set to 1 below if we determine next token is delete, not insert
    if((! defined $ins_seq_pos) && (! defined $del_seq_pos)) {
      # shouldn't happen
      ofile_FAIL("ERROR in $sub_name, ran out of insert and delete tokens, but insert token index is $ii out of $nins and delete token index is $di out of ndel\nblastn_ins_str: $blastn_ins_str\nblastn_del_str: $blastn_del_str\n", 1, $FH_HR);
    }      
    elsif(! defined $ins_seq_pos) {
      # $del_seq_pos must be defined
      $next_is_delete = 1;
    }
    elsif(! defined $del_seq_pos) {
      # $ins_seq_pos must be defined
      $next_is_insert = 1;
    }
    else {
      # both $del_seq_pos and $ins_seq_pos are defined
      # determine which comes first
      if(($ins_seq_pos < $del_seq_pos) && ($ins_mdl_pos <= $del_mdl_pos)) {
        # insert comes first
        $next_is_insert = 1;
      }
      elsif(($del_seq_pos < $ins_seq_pos) && ($del_mdl_pos <= $ins_mdl_pos)) {
        # delete comes first
        $next_is_delete = 1;
      }
      elsif(($ins_seq_pos == $del_seq_pos) && ($ins_mdl_pos == $del_mdl_pos)) {
        # shouldn't happen
        ofile_FAIL("ERROR in $sub_name, insert and delete tokens have identical query and subject positions, " . ((defined $ins_tok) ? $ins_tok : "UNDEFINED") . " and " . ((defined $del_tok) ? $del_tok : "UNDEFINED"), 1, $FH_HR);
      }
      else { 
        # shouldn't happen, seq and mdl are not in the same order
        ofile_FAIL("ERROR in $sub_name, insert and delete tokens imply query and subject coordinates out of order, " . ((defined $ins_tok) ? $ins_tok : "UNDEFINED") . " and " . ((defined $del_tok) ? $del_tok : "UNDEFINED"), 1, $FH_HR);
      }
    } # end of else entered if both $del_seq_pos and $ins_seq_pos are defined

    if($next_is_insert) {
      # sanity check that the segments we are about to add to mdl and seq
      # coords are the same length (they represent a chunk of ungapped alignment)
      if(($ins_mdl_pos - $cur_mdl_start + 1) != ($ins_seq_pos - $cur_seq_start + 1)) {
        # print("\t\tugp_mdl_coords: $ugp_mdl_coords\n");
        # print("\t\tugp_seq_coords: $ugp_seq_coords\n");
        ofile_FAIL("ERROR in $sub_name, trying to add ungapped segments before next insert, but lengths don't match up: mdl: $cur_mdl_start .. $ins_mdl_pos, seq: $cur_seq_start .. $ins_seq_pos", 1, $FH_HR);
      }
      $ugp_mdl_coords = vdr_CoordsAppendSegment($ugp_mdl_coords, vdr_CoordsSegmentCreate($cur_mdl_start, $ins_mdl_pos, "+", $FH_HR));
      $cur_mdl_start  = $ins_mdl_pos + 1;
      $ugp_seq_coords = vdr_CoordsAppendSegment($ugp_seq_coords, vdr_CoordsSegmentCreate($cur_seq_start, $ins_seq_pos, "+", $FH_HR));
      $cur_seq_start  = $ins_seq_pos + $ins_len + 1;
      $update_ins = 1;
      $ii++;
    }
    elsif($next_is_delete) {
      # sanity check that the segments we are about to add to mdl and seq
      # coords are the same length (they represent a chunk of ungapped alignment)
      if(($del_mdl_pos - $cur_mdl_start + 1) != ($del_seq_pos - $cur_seq_start + 1)) {
        ofile_FAIL("ERROR in $sub_name, trying to add ungapped segments before next delete, but lengths don't match up: mdl: $cur_mdl_start .. $del_mdl_pos, seq: $cur_seq_start .. $del_seq_pos", 1, $FH_HR);
      }
      $ugp_mdl_coords = vdr_CoordsAppendSegment($ugp_mdl_coords, vdr_CoordsSegmentCreate($cur_mdl_start, $del_mdl_pos, "+", $FH_HR));
      $cur_mdl_start  = $del_mdl_pos + $del_len + 1;
      $ugp_seq_coords = vdr_CoordsAppendSegment($ugp_seq_coords, vdr_CoordsSegmentCreate($cur_seq_start, $del_seq_pos, "+", $FH_HR));
      $cur_seq_start  = $del_seq_pos + 1;
      $update_del = 1;
      $di++;
    }
    else {
      # shouldn't happen, can't figure out what token is next
      ofile_FAIL("ERROR in $sub_name, unable to determine which insert and delete tokens imply query and subject coordinates out of order, $ins_tok_A[$ii] and $del_tok_A[$di]", 1, $FH_HR);
    }
  } # end of 'while(($ii < $nins) && ($di < $ndel))'
  # add the final segment
  if($cur_mdl_start > $mdl_stop) {
    ofile_FAIL("ERROR in $sub_name, trying to append final ungapped segment but out of model (subject) positions, trying to add $cur_mdl_start..$mdl_stop:+ to $ugp_mdl_coords", 1, $FH_HR);
  }
  if($cur_seq_start > $seq_stop) {
    ofile_FAIL("ERROR in $sub_name, trying to append final ungapped segment but out of sequence (query) positions, trying to add $cur_seq_start..$seq_stop:+ to $ugp_seq_coords", 1, $FH_HR);
  }
  # sanity check that the segments we are about to add to mdl and seq
  # coords are the same length (they represent a chunk of ungapped alignment)
  if(($mdl_stop - $cur_mdl_start + 1) != ($seq_stop - $cur_seq_start + 1)) {
    ofile_FAIL("ERROR in $sub_name, trying to add ungapped segments as final segment, but lengths don't match up: mdl: $cur_mdl_start .. $mdl_stop, seq: $cur_seq_start .. $seq_stop", 1, $FH_HR);
  }
  $ugp_mdl_coords = vdr_CoordsAppendSegment($ugp_mdl_coords, vdr_CoordsSegmentCreate($cur_mdl_start, $mdl_stop, "+", $FH_HR));
  $ugp_seq_coords = vdr_CoordsAppendSegment($ugp_seq_coords, vdr_CoordsSegmentCreate($cur_seq_start, $seq_stop, "+", $FH_HR));

  # print("in $sub_name, returning mdl: $ugp_mdl_coords, seq: $ugp_seq_coords\n");
  
  return ($ugp_mdl_coords, $ugp_seq_coords);
}
  
#################################################################
# Subroutine:  parse_blastn_indel_token()
# Incept:      EPN, Mon Mar 30 10:15:05 2020
#
# Purpose:     Parse a token from an INS or DEL output line created
#              by parse_blast.pl --program n, and return its constituent
#              parts.
#
#              if $type eq "insert": = Q<d1>:S<d2>+<d3>
#              if $type eq "delete": = Q<d1>:S<d2>-<d3>
#
# Arguments: 
#  $indel_tok:  token from parse_blast.pl INS or DEL line
#  $type:       "insert" or "delete" (just for checking if + or - is correct)
#  $FH_HR:      REF to hash of file handles, including "log" and "cmd"
#                         
# Returns:    Three values:
#             $mdl_pos: model (subject) position 
#             $seq_pos: sequence (query) position 
#             $len:     length of insert or delete
#
# Dies:       if unable to parse the token
#             if $type is "insert" and we read a - not a +
#             if $type is "delete" and we read a + not a -
#             if $type is not "insert" or "delete"
#
################################################################# 
sub parse_blastn_indel_token { 
  my $sub_name = "parse_blastn_indel_token";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($indel_tok, $type, $FH_HR) = @_;

  if(! defined $indel_tok) {
    ofile_FAIL("ERROR in $sub_name, indel_tok is undefined", 1, $FH_HR);
  }    
  if(($type ne "insert") && ($type ne "delete")) {
    ofile_FAIL("ERROR in $sub_name, type is $type but should be \"insert\" or \"delete\"", 1, $FH_HR);
  }  

  my $exp_plus_or_minus = ($type eq "insert") ? "+" : "-";

  # printf("in $sub_name indel_tok: $indel_tok\n");
  if($indel_tok =~ /^Q(\d+)\:S(\d+)([\+\-])(\d+)$/) { 
    my ($seq_pos, $mdl_pos, $plus_or_minus, $len) = ($1, $2, $3, $4);
    if($plus_or_minus ne $exp_plus_or_minus) {
      ofile_FAIL("ERROR in $sub_name, type is $type, expected $exp_plus_or_minus but read $plus_or_minus in input token $indel_tok", 1, $FH_HR);
    }
    return ($mdl_pos, $seq_pos, $len); # yes, mdl and seq are supposed to be reversed order relative to $indel_tok
  }

  # if we get here we couldn't parse the token
  ofile_FAIL("ERROR in $sub_name, unable to parse token $indel_tok", 1, $FH_HR);

  return; # NEVER REACHED
}

#################################################################
# Subroutine:  parse_blastn_indel_file_to_get_subseq_info()
# Incept:      EPN, Tue Mar 31 07:22:48 2020
#
# Purpose:     Parse a blastn indel file created by parse_blastn_results,
#              determine blastn aligned regions of each sequence that we will
#              trust, and subsequences around those regions that we will
#              align with cmalign. Return information on those subsequences
#              by filling a 2D array in @{$seqseq_AAR}, where each array
#              has 4 elements: <newname> <start> <end> <source>
#              <newname>: name to give subsequence after fetching
#              <start>:   start position of subsequence
#              <end>:     end position of subsequence
#              <source>:  name of source sequence to fetch subseq from 
#
# Arguments: 
#  $indel_file:      blastn indel file to parse, created by 
#                    parse_blastn_results() for a single model 
#  $seq_name_AR:     REF to array of sequences we want to parse indel info for
#  $seq_len_HR:      REF to hash of sequence lengths
#  $exp_mdl_name:    name of model we expect on all lines of $indel_file
#  $subseq_AAR:      REF to 2D array with subseq info, FILLED HERE
#  $ugp_mdl_HR:      REF to hash, key is <seq_name>, value is mdl coords
#                    segment of max ungapped blast aln, FILLED HERE
#  $ugp_seq_HR:      REF to hash, key is <seq_name>, value is seq coords
#                    segment of max ungapped blast aln, FILLED HERE
#  $seq2subseq_HAR:  REF to hash of arrays, key is <seq_name>,
#                    value is array of names of subsequences pertaining to
#                    <seq_name>, FILLED HERE
#  $subseq2seq_HR:   REF to hash, key is subseq name, value is <seq_name>
#                    it derives from, FILLED HERE
#  $subseq_len_HR:   REF to hash of lengths of subsequences, FILLED HERE
#  $opt_HHR:         REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information, ADDED TO HERE
#                         
# Returns:    void
#
# Dies:       if unable to parse $indel_file
#
################################################################# 
sub parse_blastn_indel_file_to_get_subseq_info { 
  my $sub_name = "parse_blastn_indel_file_to_get_subseq_info";
  my $nargs_exp = 12;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($indel_file, $seq_name_AR, $seq_len_HR, $exp_mdl_name, $subseq_AAR, $ugp_mdl_HR, $ugp_seq_HR,
      $seq2subseq_HAR, $subseq2seq_HR, $subseq_len_HR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"};
  my $nt_overhang = opt_Get("--s_overhang", $opt_HHR);

  my %processed_H = (); # key: sequence name we want indel info for, 
                        # value: 0 if we have not processed an HSP for this sequence
                        #        1 if we have
  foreach my $seq_name (@{$seq_name_AR}) { 
    $processed_H{$seq_name} = 0; 
  }
  my $nprocessed = 0;

  open(IN, $indel_file) || ofile_FileOpenFailure($indel_file, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    if($line !~ m/^#/) { 
      chomp $line;
      #model_name seq_name    mdl_coords mdl_len seq_coords seq_len  ins_str                 del_str
      #NC_039477  AB541310.1  5..7513:+  7535    1..7509:+  7509     Q41:S46+1;Q107:S111+1;  Q37:S41-1;Q113:S116-1;
      chomp $line;
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) != 8) {
        ofile_FAIL("ERROR in $sub_name, unable to parse indel line, unexpected number of tokens:\n$line\n", 1, $FH_HR);
      }
      my ($mdl_name, $seq_name, $mdl_coords, $mdl_len, $seq_coords, $seq_len, $ins_str, $del_str) = (@el_A);
      if(! defined $seq_len_HR->{$seq_name}) {
        ofile_FAIL("ERROR in $sub_name, unrecognized sequence $seq_name on line:\n$line\n", 1, $FH_HR);
      }          
      if($mdl_name ne $exp_mdl_name) { 
        ofile_FAIL("ERROR in $sub_name, unexpected model $mdl_name (expected $exp_mdl_name) on line:\n$line\n", 1, $FH_HR);
      }          
      if((defined $processed_H{$seq_name}) && ($processed_H{$seq_name} == 0)) {
        # top hit for this sequence
        my $seq_len = $seq_len_HR->{$seq_name};
        my ($ugp_mdl_coords, $ugp_seq_coords) = parse_blastn_indel_strings($mdl_coords, $seq_coords,
                                                                           $ins_str, $del_str, $FH_HR);
        my ($argmax_ugp_mdl_sgm, $max_ugp_mdl_sgm_len) = vdr_CoordsMaxLengthSegment($ugp_mdl_coords, $FH_HR);
        my ($argmax_ugp_seq_sgm, $max_ugp_seq_sgm_len) = vdr_CoordsMaxLengthSegment($ugp_seq_coords, $FH_HR);
        # sanity check: these should be the same length
        if($max_ugp_mdl_sgm_len != $max_ugp_seq_sgm_len) {
          ofile_FAIL("ERROR in $sub_name, max ungapped model segment length and ungapped sequence length differ ($max_ugp_mdl_sgm_len != $max_ugp_seq_sgm_len)", 1, $FH_HR);
        }

        # determine subseq info
        $ugp_mdl_HR->{$seq_name} = $argmax_ugp_mdl_sgm;
        $ugp_seq_HR->{$seq_name} = $argmax_ugp_seq_sgm;
        my ($ugp_seq_start, $ugp_seq_stop, $ugp_seq_strand) = vdr_CoordsSegmentParse($argmax_ugp_seq_sgm, $FH_HR);
        if(($ugp_seq_start == 1) && ($ugp_seq_stop == $seq_len)) {
          ; # do nothing, full sequence is covered by the max
          # length ungapped blast hit, no need to fetch or align with cmalign
        }
        else {
          # normal case: full seq not covered, fetch 5' end and/or 3' end
          my $start_5p = 1;
          my $stop_5p  = $ugp_seq_start + $nt_overhang - 1;
          my $start_3p = $ugp_seq_stop  - $nt_overhang + 1;
          my $stop_3p  = $seq_len;
          my $subseq_name = undef;
          if($stop_5p >= $start_3p) { # the two regions overlap, just fetch the full sequence
            $subseq_name = $seq_name . "/1-" . $seq_len;
            push(@{$subseq_AAR}, [ $subseq_name, 1, $seq_len, $seq_name ]);
            @{$seq2subseq_HAR->{$seq_name}} = ($subseq_name);
            $subseq2seq_HR->{$subseq_name} = $seq_name;
            $subseq_len_HR->{$subseq_name} = $seq_len;
          }
          else { # two regions do not overlap, fetch 5' and/or 3' ends
            @{$seq2subseq_HAR->{$seq_name}} = ();

            if($ugp_seq_start != 1) { 
              $subseq_name = $seq_name . "/" . $start_5p . "-" . $stop_5p; 
              push(@{$subseq_AAR}, [ $subseq_name, $start_5p, $stop_5p, $seq_name ]);
              push(@{$seq2subseq_HAR->{$seq_name}}, $subseq_name);
              $subseq2seq_HR->{$subseq_name} = $seq_name;
              $subseq_len_HR->{$subseq_name} = $stop_5p;
            }
            if($ugp_seq_stop != $seq_len) { 
              $subseq_name = $seq_name . "/" . $start_3p . "-" . $stop_3p; 
              push(@{$subseq_AAR}, [ $subseq_name, $start_3p, $stop_3p, $seq_name ]);
              push(@{$seq2subseq_HAR->{$seq_name}}, $subseq_name);
              $subseq2seq_HR->{$subseq_name} = $seq_name;
              $subseq_len_HR->{$subseq_name} = $stop_3p - $start_3p + 1;
            }
          }
        }
        $processed_H{$seq_name} = 1;
        $nprocessed++;
      }
    }
  }
  close(IN);

  # sanity check: we should have processed each sequence
  if($nprocessed != scalar(@{$seq_name_AR})) { 
    if($nprocessed < scalar(@{$seq_name_AR})) { 
      # shouldn't happen
      my $err_str = "ERROR in $sub_name, did not process indel strings for following sequences:\n";
      foreach my $seq_name (@{$seq_name_AR}) { 
        if($processed_H{$seq_name} != 1) { 
          $err_str .= "\t$seq_name\n";
        }
      }
      ofile_FAIL($err_str, 1, $FH_HR);
    }
    else { 
      # really shouldn't happen
      ofile_FAIL("ERROR in $sub_name, processed indel strings for more sequences than expected", 1, $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine:  join_alignments_and_add_unjoinbl_alerts()
# Incept:      EPN, Wed Apr  1 09:37:53 2020
#
# Purpose:     Join all alignments of subsequences with their ungapped
#              blastn alignments.
#
#              Report unjoinbl alerts for any sequence for which
#              we are unable to join the alignments (should be rare
#              if overhang (--s_overhang) is long enough (>50-100nt))
#
# Arguments: 
#  $sqfile_R:              REF to Bio::Easel::SqFile object, open sequence file containing the full input seqs
#  $blastn_db_sqfile_R:    REF to Bio::Easel::SqFile object, open blastn db file containing the full model seqs
#  $execs_HR:              REF to hash with paths to executables (for cmemit)
#  $do_glsearch:           '1' if we're running glsearch not cmalign
#  $cm_file:               name of CM file
#  $seq_name_AR:           REF to array of original (non subseq) sequence names
#  $seq_len_HR:            REF to hash of sequence lengths
#  $mdl_info_AHR:          REF to model info array of hashes, possibly added to here 
#  $mdl_idx:               index of model in @{mdl_info_AHR} these sequences were assigned to
#  $ugp_mdl_HR:            REF to hash, key is <seq_name>, value is mdl coords
#                          segment of max ungapped blast aln, already filled
#  $ugp_seq_HR:            REF to hash, key is <seq_name>, value is mdl coords
#                          segment of max ungapped blast aln, already filled
#  $seq2subseq_HAR:        REF to hash of arrays, key is <seq_name>,
#                          value is array of names of subsequences pertaining to
#                          <seq_name>, already filled
#  $subseq_len_HR:         REF to hash with lengths of subsequences, already filled
#  $in_stk_file_AR:        REF to array of existing stockholm files, already filled
#  $out_stk_file_AR:       REF to array of new stockholm files created here, FILLED HERE
#  $sda_output_HHR:        REF to 2D hash with information to output to .sda tabular file, ADDED TO HERE
#  $alt_seq_instances_HHR: REF to 2D hash with per-sequence alerts, PRE-FILLED
#  $alt_info_HHR:          REF to the alert info hash of arrays, PRE-FILLED
#  $unjoinbl_seq_name_AR:  REF to array of sequences with unjoinbl alerts, FILLED HERE
#  $out_root:              output root for the file names
#  $opt_HHR:               REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:        REF to 2D hash of output file information, ADDED TO HERE
#                         
# Returns:    void
#
# Dies:       if unable to parse $indel_file
#
################################################################# 
sub join_alignments_and_add_unjoinbl_alerts { 
  my $sub_name = "join_alignments_and_add_unjoinbl_alerts";
  my $nargs_exp = 22;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($sqfile_R, $blastn_db_sqfile_R, $execs_HR, $do_glsearch, $cm_file, 
      $seq_name_AR, $seq_len_HR, 
      $mdl_info_AHR, $mdl_idx, $ugp_mdl_HR, $ugp_seq_HR, 
      $seq2subseq_HAR, $subseq_len_HR, $in_stk_file_AR, 
      $out_stk_file_AR, $sda_output_HHR, 
      $alt_seq_instances_HHR, $alt_info_HHR, $unjoinbl_seq_name_AR,
      $out_root, $opt_HHR, $ofile_info_HHR) = @_;

  # printf("in $sub_name for mdl $mdl_name\n");
  
  my $FH_HR  = $ofile_info_HHR->{"FH"};
  my $do_keep = opt_Get("--keep", $opt_HHR);
  my $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
  my $mdl_len  = $mdl_info_AHR->[$mdl_idx]{"length"};
  my $mdl_consensus_sqstring;
  if(opt_Get("--keep", $opt_HHR) || opt_Get("--out_stk", $opt_HHR) || opt_Get("--out_afa", $opt_HHR) || opt_Get("--out_rpstk", $opt_HHR) || opt_Get("--out_rpafa", $opt_HHR)) { 
    # We may already have cseq from the blastn -r alignment, but even
    # if we do if we are outputting an alignment eventually then we
    # want to recreate it from cmemit because we are going to
    # eventually merge the alignments in which case we need all RF
    # annotation to be identical (equal to cmemit seq) Open all of the
    # input stk files and fetch the aligned sequence strings for all
    # sequences
    $mdl_consensus_sqstring = undef;
  }
  else { # if we're not going to output, then use it if we already have it
    $mdl_consensus_sqstring = (defined $mdl_info_AHR->[$mdl_idx]{"cseq"}) ? $mdl_info_AHR->[$mdl_idx]{"cseq"} : undef;
  }
  my $ninstk = scalar(@{$in_stk_file_AR});
  my %subseq2stk_idx_H = (); # key is subseq name, value is index of stockholm file name in @{$in_stk_file_AR}
  my %ali_subseq_H = ();     # key is subseq name, value is aligned sqstring for that subseq
  my %ali_subpp_H  = ();     # key is subseq name, value is aligned PP sqstring for that subseq, undef if $do_glsearch
  my @rf_subseq_A  = ();     # array: value $i is RF line from stockholm alignment $in_stk_file_AR->[$i]
  # insert 2D hash, filled from parsing stk file if $do_glsearch, or lower from parsing ifile if (! $do_glsearch)
  my %subseq_inserts_HH = (); # key 1: sequence name
                              # key 2: one of 'spos', 'epos', 'ins'
                              # $seq_inserts_HHR->{}{"spos"} is starting model position of alignment
                              # $seq_inserts_HHR->{}{"epos"} is ending model position of alignment
                              # $seq_inserts_HHR->{}{"ins"} is the insert string in the format:
                              # <mdlpos_1>:<uapos_1>:<inslen_1>;...<mdlpos_n>:<uapos_n>:<inslen_n>;
                              # for n inserts, where insert x is defined by:
                              # <mdlpos_x> is model position after which insert occurs 0..mdl_len (0=before first pos)
                              # <uapos_x> is unaligned sequence position of the first aligned nt
                              # <inslen_x> is length of the insert
  for(my $stk_idx = 0; $stk_idx < $ninstk; $stk_idx++) {
    my $stk_file = $in_stk_file_AR->[$stk_idx];
    my $msa = Bio::Easel::MSA->new({
      fileLocation => $stk_file,
      isDna => 1});
    my $cur_nseq = $msa->nseq();
    my $cur_rf = $msa->get_rf();
    my $cur_alen = $msa->alen();
    my @is_rf_A = split("", $cur_rf); # [0..$apos..$cur_alen-1] 1 if alignment position $apos is nongap RF position, 0 if gap RF position
    for(my $apos = 0; $apos < $cur_alen; $apos++) {
      $is_rf_A[$apos] = ($is_rf_A[$apos] =~ m/[\.\-\~]/) ? 0 : 1;
    }
    push(@rf_subseq_A, $cur_rf);
    for(my $i = 0; $i < $cur_nseq; $i++) {
      my $subseq_name = $msa->get_sqname($i);
      $subseq2stk_idx_H{$subseq_name} = $stk_idx;
      $ali_subseq_H{$subseq_name} = $msa->get_sqstring_aligned($i);
      $ali_subpp_H{$subseq_name}  = ($do_glsearch) ? undef : $msa->get_ppstring_aligned($i);
    }
    $msa = undef;
  }

  # parse the ifile for this model, we may not have one if $ninstk is 0,
  # this will happen if and only if all seqs for this model had ungapped
  # blastn hits that spanned the full sequence
  my $in_ifile = $out_root . "." . $mdl_name . ".align.ifile";
  if($ninstk > 0) { 
    vdr_CmalignParseInsertFile($in_ifile, \%subseq_inserts_HH, undef, undef, undef, undef, $FH_HR);
  }


  # define variables for the output insert file we will create
  my $out_ifile = $out_root . "." . $mdl_name . ".jalign.ifile";
  my $out_ifile_key = $mdl_name . ".jalign.ifile";
  my %seq_inserts_HH = (); # same format as subseq_inserts_HH but we fill this for joined sequences
  
  # For each sequence, determine which of the following three cases
  # (stored as $seq_case) it is:
  # Case 1: entire sequence was aligned with cmalign
  # Case 2: 5' and/or 3' ends of sequence were aligned with cmalign
  #         part of sequence covered with blastn alignment
  # Case 3: none of the sequence was aligned with cmalign,
  #         entire sequence covered by blastn alignment
  # 
  my $out_stk_idx = 0;
  foreach my $seq_name (@{$seq_name_AR}) {
    # sanity checks
    if(! defined $seq_len_HR->{$seq_name}) { ofile_FAIL("ERROR in $sub_name, no seq_len entry for sequence $seq_name", 1, $FH_HR); }
    if(! defined $ugp_mdl_HR->{$seq_name}) { ofile_FAIL("ERROR in $sub_name, no ugp_mdl entry for sequence $seq_name", 1, $FH_HR); }
    if(! defined $ugp_seq_HR->{$seq_name}) { ofile_FAIL("ERROR in $sub_name, no ugp_seq entry for sequence $seq_name", 1, $FH_HR); }

    my $seq_len = $seq_len_HR->{$seq_name};
    my $ali_seq_line = ""; # aligned full sequence string 
    my $ali_mdl_line = ""; # aligned full model string (#=GC RF annotation)
    my $ali_pp_line  = ""; # aligned full posterior probability string (#=GR PP annotation)
    my $seq_case = 3; # set to 1 or 2 if nec below
    my $alt_msg = ""; # filled if we can't join the alignments for this sequence
    
    # variables we may need to fill and send to join_alignments_helper() for cases 2 and 3
    my $ali_5p_mdl        = undef; # aligned RF (mdl) string of 5' end, if it exists
    my $ali_5p_seq        = undef; # aligned sequence string of 5' end, if it exists
    my $ali_5p_pp         = undef; # aligned PP string of 5' end, if it exists
    my $ali_5p_mdl_coords = undef; # mdl start/end points of 5' alignment, if it exists
    my $ali_5p_seq_coords = undef; # seq start/end points of 5' alignment, if it exists
    my $ali_3p_mdl        = undef; # aligned RF (mdl) string of 3' end, if it exists (from %subseq_inserts_HH from ifile)
    my $ali_3p_seq        = undef; # aligned sequence string of 3' end, if it exists
    my $ali_3p_pp         = undef; # aligned PP string of 3' end, if it exists
    my $ali_3p_mdl_coords = undef; # mdl start/end points of 5' alignment, if it exists (from %subseq_inserts_HH from ifile)
    my $ali_3p_seq_coords = undef; # seq start/end points of 5' alignment, if it exists
    my $ali_3p_seq_start  = undef; # updated if we aligned a 3' region with cmalign

    # variables we need for defining insert info for joined alignment
    %{$seq_inserts_HH{$seq_name}} = (); 
    $seq_inserts_HH{$seq_name}{"spos"} = undef; # initialize
    $seq_inserts_HH{$seq_name}{"epos"} = undef; # initialize
    $seq_inserts_HH{$seq_name}{"ins"}  = ""; # initialize

    # ungapped blastn region variables
    my ($ugp_seq_start, $ugp_seq_stop, $ugp_seq_strand) = (undef, undef, undef);
    my ($ugp_mdl_start, $ugp_mdl_stop, $ugp_mdl_strand) = (undef, undef, undef);
    my $ugp_seq = undef;
    if(defined $ugp_seq_HR->{$seq_name}) {
      ($ugp_seq_start, $ugp_seq_stop, $ugp_seq_strand) = vdr_CoordsSegmentParse($ugp_seq_HR->{$seq_name}, $FH_HR);
      ($ugp_mdl_start, $ugp_mdl_stop, $ugp_mdl_strand) = vdr_CoordsSegmentParse($ugp_mdl_HR->{$seq_name}, $FH_HR);
      $ugp_seq = $$sqfile_R->fetch_subseq_to_sqstring($seq_name, $ugp_seq_start, $ugp_seq_stop);
    }
    
    my $full_seq_idx = undef; # set to subseq idx if we have a alignment of the full sequence (case 1)
    my $ali_5p_idx   = undef; # set to subseq idx if we have a subsequence alignment on the 5' end (case 2)
    my $ali_3p_idx   = undef; # set to subseq idx if we have a subsequence alignment on the 3' end (case 2)
    if(defined $seq2subseq_HAR->{$seq_name}) { 
      # $seq2subseq_HAR->{$seq_name} is defined, this means that 
      # the longest blastn ungapped region did not span the full seq,
      # so either case 1 or 2
      $seq_case = 2; # set to 1 below if nec
      my $nsubseq = scalar(@{$seq2subseq_HAR->{$seq_name}});
      for(my $s = 0; $s < $nsubseq; $s++) { 
        my $subseq_name = $seq2subseq_HAR->{$seq_name}[$s];
        if($subseq_name =~ /^(\S+)\/(\d+)\-(\d+)$/) {
          my ($orig_seq_name, $subseq_start, $subseq_stop) = ($1, $2, $3);
          if($orig_seq_name ne $seq_name) {
            ofile_FAIL("ERROR in $sub_name, unexpected sequence name in subsequence name $subseq_name for seq $seq_name", 1, $FH_HR);
          }
          # check if we have case 1
          if(($subseq_start == 1) && ($subseq_stop == $seq_len)) { 
            $full_seq_idx = $s; # case 1
            $seq_case = 1;
            if($nsubseq != 1) {
              ofile_FAIL("ERROR in $sub_name, subsequence $subseq_name looks like full seq for seq $seq_name len $seq_len, but more than one subseqs for this seq exist", 1, $FH_HR);
            }
          }
          elsif($subseq_start == 1) {
            # check if this sequence is case 2, and this is the 5' end
            if(defined $ali_5p_idx) {
              ofile_FAIL("ERROR in $sub_name, read two aligned subseqs for 5' end of $orig_seq_name", 1, $FH_HR);
            }
            $ali_5p_idx = $s;
            $ali_5p_mdl_coords = vdr_CoordsSegmentCreate(1, $subseq_inserts_HH{$subseq_name}{"epos"}, "+", $FH_HR);
            $ali_5p_seq_coords = vdr_CoordsSegmentCreate($subseq_start, $subseq_stop, "+", $FH_HR);
          }
          elsif($subseq_stop == $seq_len) {
            # check if this sequence is case 2, and this is the 3' end
            if(defined $ali_3p_idx) {
              ofile_FAIL("ERROR in $sub_name, read two aligned subseqs for 3' end of $orig_seq_name", 1, $FH_HR);
            }
            $ali_3p_idx = $s;
            $ali_3p_mdl_coords = vdr_CoordsSegmentCreate($subseq_inserts_HH{$subseq_name}{"spos"}, $mdl_len, "+", $FH_HR);
            $ali_3p_seq_coords = vdr_CoordsSegmentCreate($subseq_start, $subseq_stop, "+", $FH_HR);
            $ali_3p_seq_start = $subseq_start;
          }
          else {
            # not 5' or 3' end or full seq, shouldn't happen
            ofile_FAIL("ERROR in $sub_name, have unexpected subseq that is none of full seq,  5' end, or 3' end, subseq name $subseq_name for seq $seq_name", 1, $FH_HR);
          }
        }
        else { # unable to parse $subseq_name
          ofile_FAIL("ERROR in $sub_name, unable to parse subsequence name $subseq_name for seq $seq_name", 1, $FH_HR);
        }
      } # end of 'for(my $s = 0; $s < $nsubseq; $s++) {' over subseqs

      # we know which case, and which aligned (sub)sequences pertain to this sequence (if case 2 or 3)
      my $stk_idx = undef;
      my $subseq_name = undef;
      if($seq_case == 1) { 
        # case 1: cmalign was used to align the full sequence
        #         we already have the fully aligned sequence and RF lines
        if(! defined $full_seq_idx) { 
          ofile_FAIL("ERROR in $sub_name, case 1, but full_seq_idx is undef on second pass for seq $seq_name", 1, $FH_HR);
        }
        $subseq_name = $seq2subseq_HAR->{$seq_name}[$full_seq_idx];
        $stk_idx = $subseq2stk_idx_H{$subseq_name};
        $ali_seq_line = $ali_subseq_H{$subseq_name};
        $ali_mdl_line = $rf_subseq_A[$stk_idx];
        $ali_pp_line  = ($do_glsearch) ? undef : $ali_subpp_H{$subseq_name};
        $seq_inserts_HH{$seq_name}{"spos"} = $subseq_inserts_HH{$subseq_name}{"spos"};
        $seq_inserts_HH{$seq_name}{"epos"} = $subseq_inserts_HH{$subseq_name}{"epos"};
        $seq_inserts_HH{$seq_name}{"ins"}  = $subseq_inserts_HH{$subseq_name}{"ins"};
        # printf("seq_case 1 $seq_name updated seq_inserts_HH to spos:" . $seq_inserts_HH{$seq_name}{"spos"} . " epos:" . $seq_inserts_HH{$seq_name}{"epos"} . " ins:" . $seq_inserts_HH{$seq_name}{"ins"} . "\n");
      }
      elsif($seq_case == 2) {
        # case 2: we don't have the full sequence aligned by cmalign,
        #         we have the 5' end, the 3' end or both, and we can infer the ungapped region
        #         get the 5' and 3' ends of the seq/mdl in preparation for a join_alignments_helper() call 
        if((! defined $ali_5p_idx) && (! defined $ali_3p_idx)) {
          ofile_FAIL("ERROR in $sub_name, unable to find subseq that is the full seq, the 5' end or the 3' end of $seq_name", 1, $FH_HR);
        }
        if(defined $ali_5p_idx) { 
          $subseq_name = $seq2subseq_HAR->{$seq_name}[$ali_5p_idx];
          $stk_idx     = $subseq2stk_idx_H{$subseq_name};
          # remove final $mdl_len - $subseq_inserts_HH{$seq_name}{"epos"}, 
          # these will certainly all be gaps, but we still may have gaps remaining after removing these
          # if there were more than one sequence in the stockholm alignment. So we do an additional step
          # of removing all terminal gaps from $ali_5p_seq (we could just do this instead of doign both
          # the substr and this step, but I think doing both is more efficient since the bulk of the gaps
          # will be removed by the substr() call which must be more efficient than a s//, right?)
          my $min_len_to_remove_at_3p_end = ($mdl_len - $subseq_inserts_HH{$subseq_name}{"epos"});
          $ali_5p_seq  = substr($ali_subseq_H{$subseq_name}, 0, (-1 * $min_len_to_remove_at_3p_end));
          $ali_5p_seq =~ s/[.\-\~]*$//; # remove trailing gaps leftover
          # now we know length of 5' segment and we can use it with substr for mdl and pp
          my $ali_5p_seq_len = length($ali_5p_seq);
          $ali_5p_pp   = ($do_glsearch) ? undef : substr($ali_subpp_H{$subseq_name},  0, $ali_5p_seq_len);
          $ali_5p_mdl  = substr($rf_subseq_A[$stk_idx],      0, $ali_5p_seq_len);
          # set insert info
          $seq_inserts_HH{$seq_name}{"spos"} = $subseq_inserts_HH{$subseq_name}{"spos"};
          if(defined $subseq_inserts_HH{$subseq_name}{"ins"}) { 
            $seq_inserts_HH{$seq_name}{"ins"} .= $subseq_inserts_HH{$subseq_name}{"ins"};
          }
          # printf("ali_5p_idx defined $seq_name updated seq_inserts_HH to spos:" . $seq_inserts_HH{$seq_name}{"spos"} . " ins:" . $seq_inserts_HH{$seq_name}{"ins"} . "\n");
        }
        else {
          $seq_inserts_HH{$seq_name}{"spos"} = $ugp_mdl_start;
          # $seq_inserts_HH{$seq_name}{"ins"} is not affected - no inserts in ungapped region
          # printf("ali_5p_idx undefined $seq_name updated seq_inserts_HH to spos:" . $seq_inserts_HH{$seq_name}{"spos"} . "\n");
        }
        if(defined $ali_3p_idx) {
          $subseq_name  = $seq2subseq_HAR->{$seq_name}[$ali_3p_idx];
          $stk_idx      = $subseq2stk_idx_H{$subseq_name};
          # remove first ($spos - 1) positions, these should be all gaps in sequence, nongap in model
          # these will certainly all be gaps, but we still may have gaps remaining after removing these
          # see comments for analogous code for 5' end above for an explanation of the logic here
          my $min_len_to_remove_at_5p_end = ($subseq_inserts_HH{$subseq_name}{"spos"} - 1);
          $ali_3p_seq  = substr($ali_subseq_H{$subseq_name}, $min_len_to_remove_at_5p_end);
          $ali_3p_seq =~ s/^[.\-\~]*//; # remove leading gaps leftover
          my $len_to_remove_at_5p_end = length($ali_subseq_H{$subseq_name}) - length($ali_3p_seq);
          $ali_3p_pp   = ($do_glsearch) ? undef : substr($ali_subpp_H{$subseq_name},  $len_to_remove_at_5p_end);
          $ali_3p_mdl  = substr($rf_subseq_A[$stk_idx],      $len_to_remove_at_5p_end);
          # set insert info
          $seq_inserts_HH{$seq_name}{"epos"} = $subseq_inserts_HH{$subseq_name}{"epos"};
          if(defined $subseq_inserts_HH{$subseq_name}{"ins"}) { 
            my $ualen_to_add = $ali_3p_seq_start -1;
            my @subseq_ins_tok_A = split(";", $subseq_inserts_HH{$subseq_name}{"ins"});
            foreach my $subseq_ins_tok (@subseq_ins_tok_A) {
              my ($subseq_mdl_pos, $subseq_uapos, $subseq_inslen) = split(":", $subseq_ins_tok);
              $seq_inserts_HH{$seq_name}{"ins"} .= $subseq_mdl_pos . ":" . ($subseq_uapos + $ali_3p_seq_start - 1) . ":" . $subseq_inslen . ";";
            }
          }
          # printf("ali_3p_idx defined $seq_name updated seq_inserts_HH to epos:" . $seq_inserts_HH{$seq_name}{"epos"} . " ins:" . $seq_inserts_HH{$seq_name}{"ins"} . "\n");
        }
        else {
          $seq_inserts_HH{$seq_name}{"epos"} = $ugp_mdl_stop;
          # $seq_inserts_HH{$seq_name}{"ins"} is not affected - no inserts in ungapped region
          # printf("ali_3p_idx undefined $seq_name updated seq_inserts_HH to epos:" . $seq_inserts_HH{$seq_name}{"epos"} . "\n");
        }
      }
      else {
        ofile_FAIL("ERROR in $sub_name, unable to determine case for sequence $seq_name", 1, $FH_HR);
      }
    }

    # finally if the ungapped region aligned by blastn was the full sequence, set spos and epos accordingly
    if($seq_case == 3) {
      $seq_inserts_HH{$seq_name}{"spos"} = $ugp_mdl_start;
      $seq_inserts_HH{$seq_name}{"epos"} = $ugp_mdl_stop;
    }
    
    # sanity check, only way we should have the aligned sequence is case 1
    if($ali_seq_line ne "") {
      if($seq_case != 1) { 
        ofile_FAIL("ERROR in $sub_name, for seq $seq_name, we have alignment prematurely", 1, $FH_HR);
      }
    }
    else {
      # case is 2 or 3
      # if case 2: ali_5p_{seq_coords,seq,mdl} and/or ali_3p_{seq_coords,seq,mdl} variables will be defined
      # if case 3: ali_5p_{seq_coords,seq,mdl} and ali_3p_{seq_coords,seq,mdl} variables will be undefined
      # 
      # for case 2 or 3, we call join_alignments_helper()
      # in case 2, this will join together the 5' and/or 3' cmalign alignments
      #            with the blastn alignment
      # in case 3, it will return the blastn alignment and construct the
      #            ungapped model/RF alignment
      # first, fetch the ungapped region of the sequence
      if(! defined $mdl_consensus_sqstring) { 
        if(! $do_glsearch) { 
          my $cseq_fa_file = $out_root . "." . $mdl_name . ".cseq.fa";
          $mdl_info_AHR->[$mdl_idx]{"cseq"} = vdr_CmemitConsensus($execs_HR, $cm_file, $mdl_name, $cseq_fa_file, $opt_HHR, $ofile_info_HHR);
          $mdl_consensus_sqstring = $mdl_info_AHR->[$mdl_idx]{"cseq"};
          ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".cseq.fa", $cseq_fa_file, 0, opt_Get("--keep", $opt_HHR), "fasta with consensus sequence for model $mdl_name");
        }
        else { 
          $mdl_consensus_sqstring = $$blastn_db_sqfile_R->fetch_seq_to_sqstring($mdl_name);
        }
      }
      ($ali_seq_line, $ali_mdl_line, $ali_pp_line) =
          join_alignments_helper($do_glsearch,
                                 $ali_5p_seq_coords, $ali_5p_mdl_coords, $ali_5p_seq, $ali_5p_mdl, $ali_5p_pp,
                                 $ali_3p_seq_coords, $ali_3p_mdl_coords, $ali_3p_seq, $ali_3p_mdl, $ali_3p_pp,
                                 $ugp_seq_HR->{$seq_name}, $ugp_mdl_HR->{$seq_name}, $ugp_seq, $mdl_consensus_sqstring, 
                                 $seq_len, $mdl_len, $ofile_info_HHR);
      if(! defined $ali_seq_line) {
        # this means something went wrong when we tried to join the alignments,
        # because the overhanging region of one or both of the cmalign alignments
        # did not match flush with the ungapped blastn alignment. This should be
        # very rare given long enough overhangs (50-100nt) but it can happen, and
        # if it does we report an unexdivg alert
        # unjoinbl: 
        $alt_msg = $ali_pp_line; # join_alignments_helper() put the alert msg into the third return value
        alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "unjoinbl", $seq_name, $alt_msg, $FH_HR);
        push(@{$unjoinbl_seq_name_AR}, $seq_name);
        # remove insert info for this seq
        delete($seq_inserts_HH{$seq_name});
      }
    }
    # now we have ali_seq_line and ali_mdl_line for this sequence, unless $alt_msg ne "" in which case
    # we added an alert, and won't align this sequence

    if($alt_msg eq "") { 
      # output the alignment to a file
      my $out_stk_file = $out_root . "." . $mdl_name . ".align.r3.s" . $out_stk_idx . ".stk";
      # printf("Writing $out_stk_file: $seq_name, ali length %d %d\n", length($ali_seq_line), length($ali_mdl_line));
      push(@{$out_stk_file_AR}, $out_stk_file);
      open(OUT, ">", $out_stk_file) || ofile_FileOpenFailure($out_stk_file, $sub_name, $!, "writing", $FH_HR);
      if($do_glsearch) { # do not include PPs
        print OUT ("# STOCKHOLM 1.0\n$seq_name $ali_seq_line\n#=GC RF $ali_mdl_line\n//\n");
      }
      else { # include PPs
        print OUT ("# STOCKHOLM 1.0\n$seq_name $ali_seq_line\n#=GR $seq_name PP $ali_pp_line\n#=GC RF $ali_mdl_line\n//\n");
      }
      close(OUT);
      $out_stk_idx++;
    }      

    # fill sda_output_HHR, to output later to .sda file in output_tabular
    # do this even if we didn't create the alignment (if $alt_msg ne "")
    $sda_output_HHR->{$seq_name}{"ugp_seq"} = $ugp_seq_HR->{$seq_name};
    $sda_output_HHR->{$seq_name}{"ugp_mdl"} = $ugp_mdl_HR->{$seq_name};
    $sda_output_HHR->{$seq_name}{"5p_seq"}  = (defined $ali_5p_seq_coords) ? $ali_5p_seq_coords : undef;
    $sda_output_HHR->{$seq_name}{"3p_seq"}  = (defined $ali_3p_seq_coords) ? $ali_3p_seq_coords : undef;
    $sda_output_HHR->{$seq_name}{"3p_mdl"}  = (defined $ali_3p_mdl_coords) ? $ali_3p_mdl_coords : undef;

    # mdl coords are not just $ali_{5p,3p}_mdl_coords, because that contains coordinates of mdl line
    # from subseq alignment but what we want is spos..epos
    $sda_output_HHR->{$seq_name}{"5p_mdl"}  = (defined $ali_5p_idx) ? 
        vdr_CoordsSegmentCreate($subseq_inserts_HH{$seq2subseq_HAR->{$seq_name}[$ali_5p_idx]}{"spos"},
                                $subseq_inserts_HH{$seq2subseq_HAR->{$seq_name}[$ali_5p_idx]}{"epos"},
                                "+", $FH_HR) : undef;
    $sda_output_HHR->{$seq_name}{"3p_mdl"}  = (defined $ali_3p_idx) ? 
        vdr_CoordsSegmentCreate($subseq_inserts_HH{$seq2subseq_HAR->{$seq_name}[$ali_3p_idx]}{"spos"},
                                $subseq_inserts_HH{$seq2subseq_HAR->{$seq_name}[$ali_3p_idx]}{"epos"},
                                "+", $FH_HR) : undef;
  } # end of 'foreach $seq_name (@{$seq_name_AR})'

  # write the new insert file, careful not to add insert info for unjoinbl seqs
  vdr_CmalignWriteInsertFile($out_ifile, 0, # do_append = 0
                             $mdl_name, $mdl_len, $seq_name_AR, $seq_len_HR,
                             \%seq_inserts_HH, $FH_HR);
    
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $out_ifile_key, $out_ifile, 0, $do_keep, sprintf("align ifile file%s", (defined $mdl_name) ? "for model $mdl_name" : ""));

  return;
}

#################################################################
# Subroutine:  join_alignments_helper()
# Incept:      EPN, Wed Apr  1 13:24:10 2020
#
# Purpose:     Figures out how to join an alignment based on
#              input and returns the joined strings.
#
# Arguments: 
#  $do_glsearch:            '1' if we used glsearch to create alignment in this
#                           case we do not have PP info
#  $ali_5p_seq_coords:      aligned 5' region sequence coords string
#                           undef if none, ugp_seq_start must be 1 in this case 
#  $ali_5p_mdl_coords:      aligned 5' region model coords string
#                           undef if none, ugp_seq_start must be 1 in this case 
#  $ali_5p_seq:             aligned 5' end of sequence from cmalign, 
#                           undef if none, ugp_seq_start must be 1 in this case 
#  $ali_5p_mdl:             aligned 5' end of model (RF) from cmalign, 
#                           undef if none, ugp_seq_start must be 1 in this case 
#  $ali_5p_pp:              aligned 5' end of posterior probability annotation from cmalign, 
#                           undef if none, ugp_seq_start must be 1 or $do_glsearch must be 1, in this case
#  $ali_3p_seq_coords:      aligned 3' region sequence coords string
#                           undef if none, ugp_seq_stop must be $seq_len in this case 
#  $ali_3p_mdl_coords:      aligned 3' region model coords string
#                           undef if none, ugp_seq_stop must be $seq_len in this case 
#  $ali_3p_pp:              aligned 3' end of posterior probability annotation from cmalign, 
#                           undef if none, ugp_seq_stop must be $seq_len or $do_glsearch must be 1, in this case
#  $ali_3p_seq:             aligned 3' end of sequence from cmalign, 
#                           undef if none, ugp_seq_stop must be $seq_len in this case 
#  $ali_3p_mdl:             aligned 3' end of model (RF) from cmalign, 
#                           undef if none, ugp_seq_stop must be $seq_len in this case 
#  $ugp_seq_coords:         ungapped region sequence coords string
#  $ugp_mdl_coords:         ungapped region model coords string
#  $ugp_seq:                ungapped region sequence string
#  $mdl_consensus_sqstring: the model consensus sequence, as a string
#  $seq_len:                total sequence length
#  $mdl_len:                total model length
#  $ofile_info_HHR:         REF to 2D hash of output file information, ADDED TO HERE
#                         
# Returns:    3 values:
#             $joined_seq: joined sequence as a string
#             $joined_mdl: joined RF as a string
#             $joined_pp:  joined PP as a string, undef if $do_glsearch
#
#             If we can't join the sequence because 3' endpoint of the 5' cmalign alignment
#             of the 5' endpoint of the cmalign alignment are inconsistent with the ungapped
#             seed region, then caller will add a unjoinbl alert, return:
#             undef,
#             undef,
#             $alert_msg: string describing the unjoinbl alert with coordinates
#    
# Dies:       if $ali_{5p,3p}_{coords,seq,mdl,pp} is undef but $ugp_seq_coords indicates it should be defined
#
################################################################# 
sub join_alignments_helper { 
  my $sub_name = "join_alignments_helper";
  my $nargs_exp = 18;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($do_glsearch, $ali_5p_seq_coords, $ali_5p_mdl_coords, $ali_5p_seq, $ali_5p_mdl, $ali_5p_pp,
      $ali_3p_seq_coords, $ali_3p_mdl_coords, $ali_3p_seq, $ali_3p_mdl, $ali_3p_pp, 
      $ugp_seq_coords, $ugp_mdl_coords, $ugp_seq, $mdl_consensus_sqstring, 
      $seq_len, $mdl_len, $ofile_info_HHR) = @_;

  my $FH_HR  = (defined $ofile_info_HHR) ? $ofile_info_HHR->{"FH"} : undef;

  my ($have_5p)  = ((defined $ali_5p_seq_coords) && (defined $ali_5p_mdl_coords) && (defined $ali_5p_seq) && (defined $ali_5p_mdl) && ((defined $ali_5p_pp) || $do_glsearch)) ? 1 : 0;
  my ($have_3p)  = ((defined $ali_3p_seq_coords) && (defined $ali_3p_mdl_coords) && (defined $ali_3p_seq) && (defined $ali_3p_mdl) && ((defined $ali_3p_pp) || $do_glsearch)) ? 1 : 0;
  my ($have_ugp) = ((defined $ugp_seq_coords) && (defined $ugp_seq)) ? 1 : 0;

  # parse coords
  my ($ugp_seq_start, $ugp_seq_stop, $ugp_seq_strand) = vdr_CoordsSegmentParse($ugp_seq_coords, $FH_HR);
  my ($ugp_mdl_start, $ugp_mdl_stop, $ugp_mdl_strand) = vdr_CoordsSegmentParse($ugp_mdl_coords, $FH_HR);
  my ($ali_5p_seq_start, $ali_5p_seq_stop, $ali_5p_seq_strand) = ($have_5p) ?
      vdr_CoordsSegmentParse($ali_5p_seq_coords, $FH_HR) :
      (undef, undef, undef);
  my ($ali_5p_mdl_start, $ali_5p_mdl_stop, $ali_5p_mdl_strand) = ($have_5p) ?
      vdr_CoordsSegmentParse($ali_5p_mdl_coords, $FH_HR) :
      (undef, undef, undef);
  my ($ali_3p_seq_start, $ali_3p_seq_stop, $ali_3p_seq_strand) = ($have_3p) ?
      vdr_CoordsSegmentParse($ali_3p_seq_coords, $FH_HR) :
      (undef, undef, undef);
  my ($ali_3p_mdl_start, $ali_3p_mdl_stop, $ali_3p_mdl_strand) = ($have_3p) ?
      vdr_CoordsSegmentParse($ali_3p_mdl_coords, $FH_HR) :
      (undef, undef, undef);
  
  # sanity checks
  if(($have_5p)   && ($ugp_seq_start == 1)) {
    ofile_FAIL("ERROR in $sub_name, ugp_seq_coords include 1st position but we have 5p segment info", 1, $FH_HR);
  }
  if((! $have_5p) && ($ugp_seq_start != 1)) {
    ofile_FAIL("ERROR in $sub_name, ugp_seq_coords do not include 1st position but we don't have 5p segment info", 1, $FH_HR);
  }
  if(($have_3p)   && ($ugp_seq_stop == $seq_len)) {
    ofile_FAIL("ERROR in $sub_name, ugp_seq_coords include final position but we have 3p segment info", 1, $FH_HR);
  }
  if((! $have_3p) && ($ugp_seq_stop != $seq_len)) {
    ofile_FAIL("ERROR in $sub_name, ugp_seq_coords do not include final position but we don't have 3p segment info", 1, $FH_HR);
  }

  # strand sanity checks
  if(($ugp_seq_strand ne "+") ||
     ($ugp_mdl_strand ne "+") || 
     ((defined $ali_5p_seq_strand) && ($ali_5p_seq_strand ne "+")) ||
     ((defined $ali_3p_seq_strand) && ($ali_3p_seq_strand ne "+"))) { 
    ofile_FAIL("ERROR in $sub_name, not all strands are positive", 1, $FH_HR);
  }

  # more sanity checks
  if(($have_5p) && ($ali_5p_mdl_start != 1)) {
    ofile_FAIL("ERROR in $sub_name, 5' model start is not 1 but $ali_5p_mdl_start", 1, $FH_HR);
  }
  if(($have_3p) && ($ali_3p_mdl_stop != $mdl_len)) {
    ofile_FAIL("ERROR in $sub_name, 3' model stop is not mdl_len ($mdl_len) but $ali_3p_mdl_stop", 1, $FH_HR);
  }
  # another sanity check: there must be some overhang between 5p and ugp and ugp and 3p
  if(($have_5p) && ($ali_5p_seq_stop < $ugp_seq_start)) {
    ofile_FAIL("ERROR in $sub_name, no overhang between ugp region ($ugp_seq_start .. $ugp_seq_stop) and 5p region ($ali_5p_seq_start .. $ali_5p_seq_stop)", 1, $FH_HR);
  }
  if(($have_3p) && ($ali_3p_seq_start > $ugp_seq_stop)) {
    ofile_FAIL("ERROR in $sub_name, no overhang between ugp region ($ugp_seq_start .. $ugp_seq_stop) and 3p region ($ali_3p_seq_start .. $ali_3p_seq_stop)", 1, $FH_HR);
  }

  # debugging print statements
  # if($have_5p) { print("ali_5p_seq: $ali_5p_seq_start .. $ali_5p_seq_stop\nali_5p_seq: $ali_5p_seq\n"); }
  # if($have_3p) { print("ali_3p_seq: $ali_3p_seq_start .. $ali_3p_seq_stop\nali_3p_seq: $ali_3p_seq\n"); }

  # if($have_5p) { print("ali_5p_mdl: $ali_5p_mdl_start .. $ali_5p_mdl_stop\n"); }
  # if($have_3p) { print("ali_3p_mdl: $ali_3p_mdl_start .. $ali_3p_mdl_stop\n"); }

  # print("ugp_seq: $ugp_seq_start .. $ugp_seq_stop\n");
  # print("ugp_mdl: $ugp_mdl_start .. $ugp_mdl_stop\n");

  ################################################
  # Check to make sure aligned overlapping (overhang) regions
  # were aligned as expected. More specifically ensure that:
  # 
  # ($ali_5p_seq_stop == ($ali_5p_mdl_stop - $ugp_seq_mdl_diff)):
  # If so, this means offset between model stop and sequence stop is
  # same as offset between model start and sequence stooo which means
  # the overlapping ('overhang' region between the 5' region and the
  # ungapped region was aligned as expected and ended with the final
  # overhanging position in the sequence aligned to the final overhanging
  # position in the model.
  #
  # Also make sure the inverse is true on the 3' end, if nec:
  # That is, that:
  # ($ali_3p_seq_start == ($ali_3p_mdl_start - $ugp_seq_mdl_diff))
  # 
  # variables starting with 'fetch' are in relative coordinate space for whatever they pertain to:
  # either ali_5p_{seq,mdl}, ali_3p_{seq,mdl} or ugp_seq
  my $ugp_seq_mdl_diff = $ugp_mdl_start - $ugp_seq_start; # offset between model start and sequence start
  # printf("ugp_seq_mdl_diff: $ugp_seq_mdl_diff\n");
  my $fetch_ali_5p_seq_start = 1;                    
  my $fetch_ali_5p_seq_stop  = length($ali_5p_seq);
  my $fetch_ali_3p_seq_start = 1;                   
  my $fetch_ali_3p_seq_stop  = length($ali_3p_seq);
  my $fetch_ugp_seq_start    = undef;
  my $fetch_ugp_seq_stop     = undef;
  my $fetch_ugp_mdl_start    = undef;
  my $fetch_ugp_mdl_stop     = undef;
  my $ugp_seq_len = ($ugp_seq_stop - $ugp_seq_start + 1);
  my $ugp_mdl_len = ($ugp_mdl_stop - $ugp_mdl_start + 1);

  my $alt_msg = ""; # added to if we can't join on 5' and/or 3' end
  if($have_5p) { 
    # usual case, final position of aligned 5' region is just 1 sequence position
    # *and* 1 model position prior to ungapped region
    if($ali_5p_seq_stop == ($ali_5p_mdl_stop - $ugp_seq_mdl_diff)) {
      $fetch_ugp_seq_start = ($ali_5p_seq_stop - $ugp_seq_start + 1) + 1; # one position past 5' overhang
      $fetch_ugp_mdl_start = $ali_5p_mdl_stop + 1; # one position past 5' overhang
    }
    else {
      $alt_msg .= "5' aligned region (mdl:$ali_5p_mdl_coords, seq:$ali_5p_seq_coords) unjoinable with seed (mdl:$ugp_mdl_coords; seq:$ugp_seq_coords);";
    }
  }
  else { # $have_5p == 0
    $fetch_ugp_seq_start = 1;
    $fetch_ugp_mdl_start = $ugp_mdl_start; # because mdl_consensus_sqstring is always length of model, and length of sequence is not
  }
  if($have_3p) {
    # usual case, first position of aligned 3' region is just 1 sequence position
    # *and* 1 model position after ungapped region
    if($ali_3p_seq_start == ($ali_3p_mdl_start - $ugp_seq_mdl_diff)) {
      $fetch_ugp_seq_stop = $ugp_seq_len - ($ugp_seq_stop - $ali_3p_seq_start + 1); # one position prior to 3' overhang
      $fetch_ugp_mdl_stop = $ali_3p_mdl_start - 1; # one position prior to 3' overhang
    }
    else {
      $alt_msg .= "3' aligned region (mdl:$ali_3p_mdl_coords, seq:$ali_3p_seq_coords) unjoinable with seed (mdl:$ugp_mdl_coords; seq:$ugp_seq_coords);";
    }
  }
  else { # $have_3p == 0
    $fetch_ugp_seq_stop = $ugp_seq_len;
    $fetch_ugp_mdl_stop = $ugp_mdl_stop; # because mdl_consensus string is always length of model, and length of sequence is not
  }

  if($alt_msg ne "") {
    # we couldn't join the alignment send the alt_msg back to caller who will add the unjoinbl alert
    return (undef, undef, $alt_msg);
  }

  ####################################################
  # Below is pseudocode for a more complicated algorithm we *could*
  # use if above checks are not satisfied
  # But in practice it's probably not necessary with 
  # long enough overhangs (50-100nt)
  # 
  # below, let ugp_seq_mdl_diff is ugp_mdl_start - ugp_seq_start
  # 
  # 5' end:
  # find X = max position <= ali_5p_stop where
  # c1: seqpos == mdlpos - ugp_seq_mdl_diff
  # where seqpos = ali_5p_seq_stop - number of nongap seq positions seen since ali_5p_seq_stop
  # where mdlpos = ali_5p_mdl_stop - number of nongap mdl positions seen since ali_5p_mdl_stop
  #
  # then
  # end ali_5p_seq at position X
  # end ali_5p_mdl at position X
  # start ungapped seq at position X+1 (chop off first X-ugp_start positions)
  #
  # 3' end:
  # find Y = min position >= ali_3p_start where
  # c2: seqpos == mdlpos - ugp_seq_mdl_diff
  # where seqpos = ali_3p_seq_start - number of nongap seq positions seen since ali_3p_seq_start
  # where mdlpos = ali_3p_mdl_start - number of nongap mdl positions seen since ali_3p_mdl_start
  # 
  # If we get to 1 and can't satisfy c1, we'd need to doctor the alignment
  # If we get to L and can't satisfy c2, we'd need to doctor the alignment
  #####################################################

  # fetch the sequence, model (RF) and PP from the 5' and 3' alignments
  # and create them for the ungapped region:
  my $apos;
  my $joined_seq = "";
  my $joined_mdl = "";
  my $joined_pp  = ($do_glsearch ? undef : "");
  if($have_5p) {
    # $fetch_ali_5p_seq_start == 1, but included below for consistency with 3p calls
    # printf("fetching 5p %d to %d from %s\n", $fetch_ali_5p_seq_start, $fetch_ali_5p_seq_stop, $ali_5p_seq_coords);
    $joined_seq .= substr($ali_5p_seq, $fetch_ali_5p_seq_start - 1, ($fetch_ali_5p_seq_stop - $fetch_ali_5p_seq_start + 1));
    $joined_mdl .= substr($ali_5p_mdl, $fetch_ali_5p_seq_start - 1, ($fetch_ali_5p_seq_stop - $fetch_ali_5p_seq_start + 1));
    if(! $do_glsearch) { 
      $joined_pp  .= substr($ali_5p_pp,  $fetch_ali_5p_seq_start - 1, ($fetch_ali_5p_seq_stop - $fetch_ali_5p_seq_start + 1));
    }
  }
  else {
    # we did not align the 5' end with cmalign, add all gap 5' chunk
    # of seq and model up to the start of ungapped blast alignment
    if($ugp_mdl_start != 1) {
      # seed-with-gaps need this to know about gaps in HSP - reparse blastn indel string?
      $joined_seq .= utl_StringMonoChar(($ugp_mdl_start - 1), "-", $FH_HR);
      #$joined_mdl .= utl_StringMonoChar(($ugp_mdl_start - 1), "x", $FH_HR);
      $joined_mdl .= substr($mdl_consensus_sqstring, 0, ($ugp_mdl_start - 1));
      if(! $do_glsearch) { 
        $joined_pp  .= utl_StringMonoChar(($ugp_mdl_start - 1), ".", $FH_HR);
      }
    }
  }
  
  # printf("fetching ugp seq %d to %d from %s\n", $fetch_ugp_seq_start, $fetch_ugp_seq_stop, $ugp_seq_coords);
  # printf("fetching ugp mdl %d to %d from %s\n", $fetch_ugp_mdl_start, $fetch_ugp_mdl_stop, $ugp_mdl_coords);
  my $fetch_ugp_seq_len = ($fetch_ugp_seq_stop - $fetch_ugp_seq_start + 1);
  my $fetch_ugp_mdl_len = ($fetch_ugp_mdl_stop - $fetch_ugp_mdl_start + 1);
  $joined_seq .= substr($ugp_seq, $fetch_ugp_seq_start - 1, $fetch_ugp_seq_len);
  #$joined_mdl .= utl_StringMonoChar($fetch_ugp_seq_len, "x", $FH_HR);
  $joined_mdl .= substr($mdl_consensus_sqstring, $fetch_ugp_mdl_start - 1, $fetch_ugp_mdl_len);
  if(! $do_glsearch) { 
    $joined_pp  .= utl_StringMonoChar($fetch_ugp_seq_len, "*", $FH_HR);
  }

  if($have_3p) {
    # printf("fetching 3p %d to %d from %s\n", $fetch_ali_3p_seq_start, $fetch_ali_3p_seq_stop, $ali_3p_seq_coords);
    $joined_seq .= substr($ali_3p_seq, $fetch_ali_3p_seq_start - 1, ($fetch_ali_3p_seq_stop - $fetch_ali_3p_seq_start + 1));
    $joined_mdl .= substr($ali_3p_mdl, $fetch_ali_3p_seq_start - 1, ($fetch_ali_3p_seq_stop - $fetch_ali_3p_seq_start + 1));
    if(! $do_glsearch) { 
      $joined_pp  .= substr($ali_3p_pp,  $fetch_ali_3p_seq_start - 1, ($fetch_ali_3p_seq_stop - $fetch_ali_3p_seq_start + 1));
    }
  }
  else { 
    # we did not align the 3' end with cmalign, add all gap 3' chunk
    # of seq and model after the end of the ungapped blast alignment
    if($ugp_mdl_stop != $mdl_len) {
      $joined_seq .= utl_StringMonoChar(($mdl_len - $ugp_mdl_stop), "-", $FH_HR);
      #$joined_mdl .= utl_StringMonoChar(($mdl_len - $ugp_mdl_stop), "x", $FH_HR);
      $joined_mdl .= substr($mdl_consensus_sqstring, $ugp_mdl_stop);
      if(! $do_glsearch) { 
        $joined_pp  .= utl_StringMonoChar(($mdl_len - $ugp_mdl_stop), ".", $FH_HR);
      }
    }
  }
  
  return ($joined_seq, $joined_mdl, $joined_pp);
}

#################################################################
# Subroutine: update_overflow_info_for_joined_alignments
# Incept:     EPN, Wed Apr  8 08:18:06 2020
# Purpose:    Given data in @{$overflow_{seq,mxsize}_AR} filled by cmalign_wrapper()
#             for subsequences of full seqs aligned due to -s, update
#             the values so they pertain to full sequences, given the
#             map from subsequences to full sequences in %{$subseq2seq_HR}.
#
# Arguments:
#  $sub_overflow_seq_AR:      REF to array of subseq names we had overflows for, ALREADY FILLED
#  $sub_overflow_mxsize_AR:   REF to array of mxsizes of overflows, ALREADY FILLED
#  $subseq2seq_HR:            REF to hash mapping subsequence names to full seq names, ALREADY FILLED
#  $full_overflow_seq_AR:     REF to array of full seq names we have overflows for, FILLED HERE
#  $full_overflow_mxsize_AR:  REF to array of mxsizes of overflows for full seqs, FILLED HERE
# 
# Returns:  void, fills @{$full_overflow_seq_AR} and @{$full_overflow_mxsize_AR}
#
# Dies:     never
#
#################################################################
sub update_overflow_info_for_joined_alignments { 
  my $sub_name = "update_overflow_info_for_joined_alignments";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sub_overflow_seq_AR, $sub_overflow_mxsize_AR, $subseq2seq_HR, $full_overflow_seq_AR, $full_overflow_mxsize_AR) = (@_);

  my %added_H = (); # used so we don't add overflow for same full seq twice
  for(my $i = 0; $i < scalar(@{$sub_overflow_seq_AR}); $i++) {
    my $subseq_name = $sub_overflow_seq_AR->[$i];
    if(defined $subseq2seq_HR->{$subseq_name}) { 
      my $seq_name = $subseq2seq_HR->{$subseq_name};
      if(! defined $added_H{$seq_name}) { 
        push(@{$full_overflow_seq_AR},    $seq_name);
        push(@{$full_overflow_mxsize_AR}, $sub_overflow_mxsize_AR->[$i]);
        $added_H{$seq_name} = 1;
      }
    }
  }

  return;
}

###########################################################################
# the next line is critical, a perl module must return a true value
return 1;
###########################################################################

    
