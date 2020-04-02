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

require "vadr.pm"; 

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
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $db_file, $seq_file, $out_root, 
      $nseq, $progress_w, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"};
  my $log_FH = $FH_HR->{"log"}; # for convenience

  my $do_keep = opt_Get("--keep", $opt_HHR);

  my $start_secs = ofile_OutputProgressPrior(sprintf("Classifying sequences with blastn ($nseq seq%s)", ($nseq > 1) ? "s" : ""), $progress_w, $log_FH, *STDOUT);

  my $blastn_out_file = $out_root . ".r1.blastn.out";
  my $opt_str = "-query $seq_file -db $db_file -out $blastn_out_file -word_size " . opt_Get("--blastnws", $opt_HHR); 
  my $blastn_cmd = $execs_HR->{"blastn"} . " $opt_str";
  
  utl_RunCommand($blastn_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "blastn-out", $blastn_out_file, 0, $do_keep, "blastn output");

  # now summarize its output
  # use --splus to skip alignment parsing of hits to negative strand of subject
  # that's important because parse_blast.pl doesn't have the code to deal with
  # such alignments (it was written for blastx originally for which subject is
  # always +) but that's okay because those hits are to negative strand of
  # the sequence (actually negative strand of the subject/model but blastn revcomps
  # the subject instead of the query like cmscan would do), and we don't care
  # about negative strand hit indel info.
  my $blastn_summary_file = $out_root . ".r1.blastn.summary.txt";
  my $parse_cmd = $execs_HR->{"parse_blast"} . " --program n --input $blastn_out_file --splus > $blastn_summary_file";
  utl_RunCommand($parse_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "blastn-summary", $blastn_summary_file, 0, $do_keep, "parsed (summarized) blastn output");

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
#             "blastn.r1.pretblout" file: cmscan --trmF3 output format
#             file with each hit on a separate line and individual hit
#             scores reported for each hit. This is later processed by
#             blastn_pretblout_to_tblout() to sum bit scores for all
#             hits with the same model/seq/strand so we can classify
#             sequences the same way we do in default mode (cmscan
#             based classification).
#
#             Output model 2: 2 files produced per model with >= 1 matching
#             sequence, produced if $seq2mdl_HR is defined.
#
#             "search.r2.<mdlname>.tblout": cmsearch --tblout format
#             file with each hit for a sequence on + strand that is 
#             classified to model <mdlname>. 
#
#             "blastn.r2.<mdlname>.indel.txt": one line per sequence
#             with all inserts and deletes in all blastn hit
#             alignments for each sequence that is classified to 
#             <mdlname> on strand +.
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
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastn_summary_file, $seq_len_HR, $seq2mdl_HR, $mdl_name_AR, 
      $out_root, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $pretblout_FH = undef; # defined if output mode 1 (if ! defined $seq2mdl_HR)
  my %tblout_FH_H  = ();    # defined if output mode 2 (if   defined $seq2mdl_HR)
  my %indel_FH_H   = ();    # defined if output mode 2 (if   defined $seq2mdl_HR)
  my $outfile_key  = undef; # a key for an output file in %{$ofile_info_HHR}
  my $small_value  = 0.000001;
  my $min_bitsc    = opt_Get("--blastnsc", $opt_HHR) - $small_value;

  my $mdl_name = undef;
  if(! defined $seq2mdl_HR) { 
    # output mode 1, open the pretblout output file 
    # 
    # We write to this file as we parse the $blastn_summary_file but we
    # have to post-process it in blastn_pretblout_to_tblout so that the
    # top hit per model/sequence/strand trio includes the *summed* score
    # for that model/strand/strand instead of just the hit score. This
    # way we will match the cmscan --trmF3 output downstream steps
    # expect.
    ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, "blastn.r1.pretblout", $out_root . ".blastn.r1.pretblout",  0, opt_Get("--keep", $opt_HHR), "blastn output converted to cmscan --trmF3 tblout format (hit scores)");
    $pretblout_FH = $ofile_info_HHR->{"FH"}{"blastn.r1.pretblout"}; 
    printf $pretblout_FH ("%-30s  %-30s  %8s  %9s  %9s  %6s  %6s  %3s  %11s\n", 
                          "#modelname/subject", "sequence/query", "bitscore", "start", "end", "strand", "bounds", "ovp", "seqlen");
  }
  else { # $seq2mdl_HR is defined
    # output mode 2, for each model in @{$mdl_name_AR}, open the 
    # coverage determination tblout files in cmsearch --tblout 
    # format (not --trmF3 output format) and the indel files 
    foreach $mdl_name (@{$mdl_name_AR}) { 
      $outfile_key = "search.r2.$mdl_name.tblout";
      ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, $outfile_key, $out_root . "." . $outfile_key,  0, opt_Get("--keep", $opt_HHR), "blastn output converted to cmsearch tblout format for model $mdl_name");
      $tblout_FH_H{$mdl_name} = $ofile_info_HHR->{"FH"}{$outfile_key};

      $outfile_key = "search.r2.$mdl_name.indel";
      ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, $outfile_key, $out_root . "." . $outfile_key,  0, opt_Get("--keep", $opt_HHR), "blastn indel information for model $mdl_name");
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
        if(($value !~ /^\d+\.\d+$/) && ($value !~ /^\d+$/)) { 
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
        elsif($value =~ /^(\d+)..(\d+)$/) { 
          ($cur_H{"QRANGESTART"}, $cur_H{"QRANGESTOP"}) = ($1, $2);
        }
        else { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, unable to parse blastn summary QRANGE line $line", 1, $FH_HR);
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
        elsif(($value =~ /^(\d+)..(\d+)$/) && ($cur_H{"BITSCORE"} >= $min_bitsc)) { 
          ($cur_H{"SRANGESTART"}, $cur_H{"SRANGESTOP"}) = ($1, $2);

          # output data in cmscan --trmF3 format
          if((! defined $cur_H{"QRANGESTART"}) || (! defined $cur_H{"QRANGESTOP"})) { 
            ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read $key line before QRANGE line (line: $line_idx)\n", 1, $FH_HR);
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
    # output mode 1:
    # close the pretblout file
    # then call convert it pretblout to tblout (cmscan --trmF3 format)
    # which will have scores summed for each seq/mdl/strand trio
    close $ofile_info_HHR->{"FH"}{"blastn.r1.pretblout"};
    blastn_pretblout_to_tblout($ofile_info_HHR->{"fullpath"}{"blastn.r1.pretblout"}, 
                               \%scsum_HHH, $out_root, $opt_HHR, $ofile_info_HHR);
  }
  else { 
    # output mode 2:
    # close the per model files:
    foreach $mdl_name (@{$mdl_name_AR}) { 
      $outfile_key = "search.r2.$mdl_name.tblout";
      close $ofile_info_HHR->{"FH"}{$outfile_key};
      $outfile_key = "search.r2.$mdl_name.indel";
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
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastn_pretblout_file, $scsum_HHHR, $out_root, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, "scan.r1.tblout", $out_root . ".blastn.r1.tblout",  0, opt_Get("--keep", $opt_HHR), "blastn output converted to cmscan --trmF3 tblout format (summed hit scores)");
  my $tblout_FH = $FH_HR->{"scan.r1.tblout"}; 

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
                         $model, $seq, $scsum_HHHR->{$model}{$seq}{$strand}, 
                         $start, $end, $strand, $bounds, $ovp, $seqlen);
      # set scsum to zero for all subsequent hits to this trio      
      $scsum_HHHR->{$model}{$seq}{$strand} = 0.; 
    }
    else { 
      print $tblout_FH ($line . "\n");
    }
  }
  close(IN);
  close $ofile_info_HHR->{"FH"}{"scan.r1.tblout"};
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
#  $seq_len_HR:      REF to hash of sequence lengths
#  $exp_mdl_name:    name of model we expect on all lines of $indel_file
#  $subseq_AAR:      REF to 2D array with subseq info, FILLED HERE
#  $ugp_mdl_HR:      REF to hash, key is <seq_name>, value is mdl coords
#                    segment of max ungapped blast aln, FILLED HERE
#  $ugp_seq_HR:      REF to hash, key is <seq_name>, value is mdl coords
#                    segment of max ungapped blast aln, FILLED HERE
#  $seq2subseq_HAR:  REF to hash of arrays, key is <seq_name>,
#                    value is array of names of subsequences pertaining to
#                    <seq_name>, FILLED HERE
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
  my $nargs_exp = 10;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($indel_file, $seq_len_HR, $exp_mdl_name, $subseq_AAR, $ugp_mdl_HR, $ugp_seq_HR,
      $seq2subseq_HAR, $subseq_len_HR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"};
  my $nt_overhang = opt_Get("--overhang", $opt_HHR);

  my %seen_H = (); # key: sequence name, value: 1 if we've already processed an HSP for this sequence
  
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
      if(! defined $seen_H{$seq_name}) {
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
          if($stop_5p > $start_3p) { # the two regions overlap, just fetch the full sequence
            $subseq_name = $seq_name . "/1-" . $seq_len;
            push(@{$subseq_AAR}, [ $subseq_name, 1, $seq_len, $seq_name ]);
            @{$seq2subseq_HAR->{$seq_name}} = ($subseq_name);
            $subseq_len_HR->{$seq_name} = $seq_len;
          }
          else { # two regions do not overlap, fetch 5' and/or 3' ends
            @{$seq2subseq_HAR->{$seq_name}} = ();

            if($ugp_seq_start != 1) { 
              $subseq_name = $seq_name . "/" . $start_5p . "-" . $stop_5p; 
              push(@{$subseq_AAR}, [ $subseq_name, $start_5p, $stop_5p, $seq_name ]);
              push(@{$seq2subseq_HAR->{$seq_name}}, $subseq_name);
              $subseq_len_HR->{$subseq_name} = $stop_5p;
            }
            if($ugp_seq_stop != $seq_len) { 
              $subseq_name = $seq_name . "/" . $start_3p . "-" . $stop_3p; 
              push(@{$subseq_AAR}, [ $subseq_name, $start_3p, $stop_3p, $seq_name ]);
              push(@{$seq2subseq_HAR->{$seq_name}}, $subseq_name);
              $subseq_len_HR->{$subseq_name} = $stop_3p - $start_3p + 1;
            }
          }
        }
        $seen_H{$seq_name} = 1;
      }
    }
  }
  close(IN);

  return;
}

#################################################################
# Subroutine:  join_alignments()
# Incept:      EPN, Wed Apr  1 09:37:53 2020
#
# Purpose:     Join all alignments of subsequences with their ungapped
#              blastn alignments.
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
#  $sqfile:          REF to Bio::Easel::SqFile object, open sequence file containing the full input seqs
#  $seq_name_AR:     REF to array of original (non subseq) sequence names
#  $seq_len_HR:      REF to hash of sequence lengths
#  $mdl_name:        name of model these sequences were assigned to
#  $mdl_len:         length of model these sequences were assigned to
#  $ugp_mdl_HR:      REF to hash, key is <seq_name>, value is mdl coords
#                    segment of max ungapped blast aln, already filled
#  $ugp_seq_HR:      REF to hash, key is <seq_name>, value is mdl coords
#                    segment of max ungapped blast aln, already filled
#  $seq2subseq_HAR:  REF to hash of arrays, key is <seq_name>,
#                    value is array of names of subsequences pertaining to
#                    <seq_name>, already filled
#  $subseq_len_HR:   REF to hash with lengths of subsequences, already filled
#  $stk_file_AR:     ref to array of stockholm files, already filled
#  $progress_w:      width for outputProgressPrior output
#  $opt_HHR:         REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information, ADDED TO HERE
#                         
# Returns:    void
#
# Dies:       if unable to parse $indel_file
#
################################################################# 
sub join_alignments { 
  my $sub_name = "join_alignments";
  my $nargs_exp = 13;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($sqfile, $seq_name_AR, $seq_len_HR, $mdl_name, $mdl_len, $ugp_mdl_HR, $ugp_seq_HR, $seq2subseq_HAR, $subseq_len_HR,
      $stk_file_AR, $progress_w, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"};

  my $nseq = scalar(@{$seq_name_AR});
  my $start_secs = ofile_OutputProgressPrior(sprintf("Joining alignments from cmalign and blastn for model $mdl_name ($nseq seq%s)", ($nseq > 1) ? "s" : ""), $progress_w, $FH_HR->{"log"}, *STDOUT);

  # get the aligned sequence strings for all sequences in the stockholm files
  my $nstk = scalar(@{$stk_file_AR});
  my %subseq2stk_idx_H = (); # key is subseq name, value is index of stockholm file name in @{$stk_file_AR}
  my %asubseq_H = ();        # key is subseq name, value is aligned sqstring for that subseq
  my @rf_seq_A = ();         # array: value $i is RF line from stockholm alignment $stk_file_AR->[$i]
  for(my $stk_idx = 0; $stk_idx < $nstk; $stk_idx++) {
    my $stk_file = $stk_file_AR->[$stk_idx];
    my $msa = Bio::Easel::MSA->new({
      fileLocation => $stk_file,
      isDna => 1});
    my $nseq = $msa->nseq();
    push(@rf_seq_A, $msa->get_rf());
    for(my $i = 0; $i < $nseq; $i++) {
      my $subseq_name = $msa->get_sqname($i);
      $subseq2stk_idx_H{$subseq_name} = $stk_idx;
      $asubseq_H{$subseq_name} = $msa->get_sqstring_aligned($i);
    }
    $msa = undef;
  }

  # join alignments for each sequence
  foreach my $seq_name (@{$seq_name_AR}) {
    if(! defined $seq_len_HR->{$seq_name})     { ofile_FAIL("ERROR in $sub_name, no seq_len entry for sequence $seq_name", 1, $FH_HR); }
    if(! defined $ugp_mdl_HR->{$seq_name})     { ofile_FAIL("ERROR in $sub_name, no ugp_mdl entry for sequence $seq_name", 1, $FH_HR); }
    if(! defined $ugp_seq_HR->{$seq_name})     { ofile_FAIL("ERROR in $sub_name, no ugp_seq entry for sequence $seq_name", 1, $FH_HR); }
    my $seq_len = $seq_len_HR->{$seq_name};
    my $seq_line = "";
    my $rf_line = "";

    if(defined $seq2subseq_HAR->{$seq_name}) {
      # longest ungapped blastn alignment did not cover the full sequence
      my $subseq_5p_idx = -1; # set to subseq idx if we have a subsequence alignment on the 5' end
      my $subseq_3p_idx = -1; # set to subseq idx if we have a subsequence alignment on the 3' end
      my $full_seq_idx  = -1; # set to subseq idx if we have a alignment of the full sequence
      my $nsubseq = scalar(@{$seq2subseq_HAR->{$seq_name}});
      for(my $s = 0; $s < $nsubseq; $s++) { 
        my $subseq_name = $seq2subseq_HAR->{$seq_name}[$s];
        if($subseq_name =~ /^(\S+)\/(\d+)\-(\d+)$/) {
          my ($orig_seq_name, $subseq_start, $subseq_stop) = ($1, $2, $3);
          if($orig_seq_name ne $seq_name) {
            ofile_FAIL("ERROR in $sub_name, unexpected sequence name in subsequence name $subseq_name for seq $seq_name", 1, $FH_HR);
          }
          if(($subseq_start == 1) && ($subseq_stop == $seq_len)) { 
            $full_seq_idx = $s;
            if($nsubseq != 1) { 
              ofile_FAIL("ERROR in $sub_name, subsequence $subseq_name looks like full seq for seq $seq_name len $seq_len, but more than one subseqs for this seq exist", 1, $FH_HR);
            }
          }
          elsif($subseq_start == 1) {
            $subseq_5p_idx = $s;
          }
          elsif($subseq_stop == $seq_len) {
            $subseq_3p_idx = $s;
          }
          else {
            ofile_FAIL("ERROR in $sub_name, have unexpected subseq that is none of full seq,  5' end, or 3' end, subseq name $subseq_name for seq $seq_name", 1, $FH_HR);
          }
        }
        else {
          ofile_FAIL("ERROR in $sub_name, unable to parse subsequence name $subseq_name for seq $seq_name", 1, $FH_HR);
        }
      }# end of 'for(my $s = 0; $s < $nsubseq; $s++) {' over subseqs

      # create the alignment for this sequence
      my $stk_idx = undef;
      my $subseq_name = undef;
      if($full_seq_idx != -1) {
        $subseq_name = $seq2subseq_HAR->{$seq_name}[$full_seq_idx];
        $seq_line .= $asubseq_H{$subseq_name};

        $stk_idx = $subseq2stk_idx_H{$subseq_name};
        $rf_line .= $rf_seq_A[$stk_idx];
      }
      else {
        if($subseq_5p_idx != -1) {
          $subseq_name = $seq2subseq_HAR->{$seq_name}[$subseq_5p_idx];
          $seq_line .= $asubseq_H{$subseq_name};

          $stk_idx = $subseq2stk_idx_H{$subseq_name};
          $rf_line .= $rf_seq_A[$stk_idx];
        }
        
        # add ungapped region detected by blastn
        my ($ugp_seq_start, $ugp_seq_stop, $ugp_seq_strand) = vdr_CoordsSegmentParse($ugp_seq_HR->{$seq_name}, $FH_HR);
        if($ugp_seq_strand ne "+") { ofile_FAIL("ERROR in $sub_name, ungapped sequence segment $ugp_seq_start .. $ugp_seq_stop is not + strand", 1, $FH_HR); }
        $seq_line .= $sqfile->fetch_subseq_to_sqstring($seq_name, $ugp_seq_start, $ugp_seq_stop);
        $rf_line  .= utl_StringMonoChar(abs($ugp_seq_stop - $ugp_seq_start) + 1, "x", $FH_HR);
        
        if($subseq_3p_idx != -1) {
          $subseq_name = $seq2subseq_HAR->{$seq_name}[$subseq_3p_idx];
          $seq_line .= $asubseq_H{$subseq_name};
          
          $stk_idx = $subseq2stk_idx_H{$subseq_name};
          $rf_line .= $rf_seq_A[$stk_idx];
        }
      }
    } # end of 'if(defined $seq2subseq_HAR->{$seq_name})'
    else {
      # longest ungapped blastn alignment did cover the full sequence
        my ($ugp_seq_start, $ugp_seq_stop, $ugp_seq_strand) = vdr_CoordsSegmentParse($ugp_seq_HR->{$seq_name}, $FH_HR);
        my ($ugp_mdl_start, $ugp_mdl_stop, $ugp_mdl_strand) = vdr_CoordsSegmentParse($ugp_mdl_HR->{$seq_name}, $FH_HR);
        if($ugp_seq_strand ne "+") { ofile_FAIL("ERROR in $sub_name, ungapped sequence segment $ugp_seq_start .. $ugp_seq_stop is not + strand", 1, $FH_HR); }
        if($ugp_mdl_strand ne "+") { ofile_FAIL("ERROR in $sub_name, ungapped sequence segment $ugp_mdl_start .. $ugp_mdl_stop is not + strand", 1, $FH_HR); }
        if($ugp_mdl_start > 1) {
          $rf_line  .= utl_StringMonoChar($ugp_mdl_start - 1, "x", $FH_HR);
          $seq_line .= utl_StringMonoChar($ugp_mdl_start - 1, ".", $FH_HR);
        }
        $rf_line  .= utl_StringMonoChar($ugp_mdl_stop - $ugp_mdl_start + 1, "x", $FH_HR);
        $seq_line .= $sqfile->fetch_subseq_to_sqstring($seq_name, $ugp_seq_start, $ugp_seq_stop);
        if($ugp_mdl_stop < $mdl_len) { 
          $rf_line  .= utl_StringMonoChar($mdl_len - $ugp_mdl_stop - 1, "x", $FH_HR);
          $seq_line .= utl_StringMonoChar($mdl_len - $ugp_mdl_stop - 1, ".", $FH_HR);
        }
    }
    printf("# STOCKHOLM 1.0\n");
    printf("$seq_name $seq_line\n");
    printf("#=GC RF $rf_line\n");
  } # end of 'foreach $seq_name (@{$seq_name_AR})'

  ofile_OutputProgressComplete($start_secs, undef, $FH_HR->{"log"}, *STDOUT);

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
#  $ali_5p_seq_coords:      aligned 5' region coords string
#                           undef if none, ugp_seq_start must be 1 in this case 
#  $ali_5p_seq:             aligned 5' end of sequence from cmalign, 
#                           undef if none, ugp_seq_start must be 1 in this case 
#  $ali_5p_mdl:             aligned 5' end of model (RF) from cmalign, 
#                           undef if none, ugp_seq_start must be 1 in this case 
#  $ali_3p_seq_coords:      aligned 3' region coords string
#                           undef if none, ugp_seq_stop must be $seq_len in this case 
#  $ali_3p_seq:             aligned 3' end of sequence from cmalign, 
#                           undef if none, ugp_seq_stop must be $seq_len in this case 
#  $ali_3p_mdl:             aligned 3' end of model (RF) from cmalign, 
#                           undef if none, ugp_seq_stop must be $seq_len in this case 
#  $ugp_seq_coords:         ungapped region sequence coords string
#  $ugp_mdl_coords:         ungapped region model coords string
#  $ugp_seq:                ungapped region sequence string
#  $seq_len:                total sequence length
#  $mdl_len:                total model length
#  $ofile_info_HHR:         REF to 2D hash of output file information, ADDED TO HERE
#                         
# Returns:    2 values:
#             $joined_seq: joined sequence as a string
#             $joined_mdl: joined RF as a string
#
# Dies:       if $ali_{5p,3p}_{coords,seq,mdl} is undef but $ugp_seq_coords indicates it should be defined
#
################################################################# 
sub join_alignments_helper { 
  my $sub_name = "join_alignments_helper";
  my $nargs_exp = 12;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ali_5p_seq_coords, $ali_5p_seq, $ali_5p_mdl,
      $ali_3p_seq_coords, $ali_3p_seq, $ali_3p_mdl,
      $ugp_seq_coords, $ugp_mdl_coords, $ugp_seq,
      $seq_len, $mdl_len, $ofile_info_HHR) = @_;

  my $FH_HR  = (defined $ofile_info_HHR) ? $ofile_info_HHR->{"FH"} : undef;

  my ($have_5p)  = ((defined $ali_5p_seq_coords) && (defined $ali_5p_seq) && (defined $ali_5p_mdl)) ? 1 : 0;
  my ($have_3p)  = ((defined $ali_3p_seq_coords) && (defined $ali_3p_seq) && (defined $ali_3p_mdl)) ? 1 : 0;
  my ($have_ugp) = ((defined $ugp_seq_coords) && (defined $ugp_seq)) ? 1 : 0;

  # parse coords
  my ($ugp_seq_start, $ugp_seq_stop, $ugp_seq_strand) = vdr_CoordsSegmentParse($ugp_seq_coords, $FH_HR);
  my ($ugp_mdl_start, $ugp_mdl_stop, $ugp_mdl_strand) = vdr_CoordsSegmentParse($ugp_mdl_coords, $FH_HR);
  my ($ali_5p_seq_start, $ali_5p_seq_stop, $ali_5p_seq_strand) = ($have_5p) ?
      vdr_CoordsSegmentParse($ali_5p_seq_coords, $FH_HR) :
      (undef, undef, undef);
  my ($ali_3p_seq_start, $ali_3p_seq_stop, $ali_3p_seq_strand) = ($have_3p) ?
      vdr_CoordsSegmentParse($ali_3p_seq_coords, $FH_HR) :
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

  # determine model length in 5' and 3' segments
  my $ali_5p_mdl_start = ($have_5p) ? 1        : undef;
  my $ali_5p_mdl_stop  = ($have_5p) ? ($ali_5p_mdl =~ tr/[.\-~]//c) : undef;
  my $ali_3p_mdl_start = ($have_5p) ? ($mdl_len - ($ali_3p_mdl =~ tr/[.\-~]//c) + 1) : undef;
  my $ali_3p_mdl_stop  = ($have_3p) ? $mdl_len : undef;

  # another sanity check: there must be some overhang between 5p and ugp and ugp and 3p
  if(($have_5p) && ($ali_5p_seq_stop < $ugp_seq_start)) {
    ofile_FAIL("ERROR in $sub_name, no overhang between ugp region ($ugp_seq_start .. $ugp_seq_stop) and 5p region ($ali_5p_seq_start .. $ali_5p_seq_stop)", 1, $FH_HR);
  }
  if(($have_3p) && ($ali_3p_seq_start > $ugp_seq_stop)) {
    ofile_FAIL("ERROR in $sub_name, no overhang between ugp region ($ugp_seq_start .. $ugp_seq_stop) and 3p region ($ali_3p_seq_start .. $ali_3p_seq_stop)", 1, $FH_HR);
  }
  
  print("ali_5p_mdl: $ali_5p_mdl_start .. $ali_5p_mdl_stop\n");
  print("ali_3p_mdl: $ali_3p_mdl_start .. $ali_3p_mdl_stop\n");

  my $ugp_seq_mdl_diff = $ugp_mdl_start - $ugp_seq_start;
  printf("ugp_seq_mdl_diff: $ugp_seq_mdl_diff\n");

  # where ugp_seq_mdl_diff is ugp_mdl_start - ugp_seq_start
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
  # WHAT IF WE GET TO 1 AND CAN'T SATISFY c1? NEED TO DOCTOR ALIGNMENT
  # WHAT IF WE GET TO L AND CAN'T SATISFY c2? NEED TO DOCTOR ALIGNMENT
  # 

  # find 3' end of 5' end of ungapped segment
  my $apos;

  # variables starting with 'fetch' are in relative coordinate space for whatever they pertain to:
  # either ali_5p_{seq,mdl}, ali_3p_{seq,mdl} or ugp_{seq}
  my $fetch_ali_5p_seq_start = 1;                    
  my $fetch_ali_5p_seq_stop  = length($ali_5p_seq);
  my $fetch_ali_3p_seq_start = 1;                   
  my $fetch_ali_3p_seq_stop  = length($ali_3p_seq);
  my $fetch_ugp_seq_start    = undef;
  my $fetch_ugp_seq_stop     = undef;
  my $ugp_seq_len = ($ugp_seq_stop - $ugp_seq_start + 1);

  if($have_5p) { 
    # usual case, final position of aligned 5' region is just 1 sequence position
    # *and* 1 model position prior to ungapped region
    if($ali_5p_seq_stop == ($ali_5p_mdl_stop - $ugp_seq_mdl_diff)) {
      $fetch_ugp_seq_start = ($ali_5p_seq_stop - $ugp_seq_start + 1) + 1; # one position past 5' overhang
    }
    else {
      die "NOT YET IMPLEMENTED, mismatch at 5' end";
    }
  }
  else { # $have_5p == 0
    $fetch_ugp_seq_start = $ugp_seq_start; # this is 1
  }
  if($have_3p) {
    # usual case, first position of aligned 3' region is just 1 sequence position
    # *and* 1 model position after ungapped region
    if($ali_3p_seq_start == ($ali_3p_mdl_start - $ugp_seq_mdl_diff)) {
      $fetch_ugp_seq_stop = $ugp_seq_len - ($ugp_seq_stop - $ali_3p_seq_start + 1); # one position prior to 3' overhang
    }
    else {
      die "NOT YET IMPLEMENTED, mismatch at 3' end";
    }
  }
  else { # $have_3p == 0
    $fetch_ugp_seq_stop = $ugp_seq_len - $fetch_ugp_seq_start + 1;
  }

  my $joined_seq = "";
  my $joined_mdl = "";
  if($have_5p) {
    # $fetch_ali_5p_seq_start == 1, but included below for consistency with 3p calls
    printf("fetching 5p %d to %d from %s\n", $fetch_ali_5p_seq_start, $fetch_ali_5p_seq_stop, $ali_5p_seq_coords);
    $joined_seq .= substr($ali_5p_seq, $fetch_ali_5p_seq_start - 1, ($fetch_ali_5p_seq_stop - $fetch_ali_5p_seq_start + 1));
    $joined_mdl .= substr($ali_5p_mdl, $fetch_ali_5p_seq_start - 1, ($fetch_ali_5p_seq_stop - $fetch_ali_5p_seq_start + 1));
  }

  printf("fetching ugp %d to %d from %s\n", $fetch_ugp_seq_start, $fetch_ugp_seq_stop, $ugp_seq_coords);
  my $fetch_ugp_seq_len = ($fetch_ugp_seq_stop - $fetch_ugp_seq_start + 1);
  $joined_seq .= substr($ugp_seq, $fetch_ugp_seq_start - 1, $fetch_ugp_seq_len);
  $joined_mdl .= utl_StringMonoChar($fetch_ugp_seq_len, "x", $FH_HR);

  if($have_3p) {
    printf("fetching 3p %d to %d from %s\n", $fetch_ali_3p_seq_start, $fetch_ali_3p_seq_stop, $ali_3p_seq_coords);
    $joined_seq .= substr($ali_3p_seq, $fetch_ali_3p_seq_start - 1, ($fetch_ali_3p_seq_stop - $fetch_ali_3p_seq_start + 1));
    $joined_mdl .= substr($ali_3p_mdl, $fetch_ali_3p_seq_start - 1, ($fetch_ali_3p_seq_stop - $fetch_ali_3p_seq_start + 1));
  }
  
  return ($joined_seq, $joined_mdl);
}

###########################################################################
# the next line is critical, a perl module must return a true value
return 1;
###########################################################################

    
