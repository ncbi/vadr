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
  my $opt_str = "-query $seq_file -db $db_file -out $blastn_out_file";
  if(opt_IsUsed("--blastnws", $opt_HHR)) { $opt_str .= " -word_size " . opt_Get("--blastnws", $opt_HHR); }
  my $blastn_cmd = $execs_HR->{"blastn"} . " $opt_str";
  
  utl_RunCommand($blastn_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "blastn-out", $blastn_out_file, 0, $do_keep, "blastn output");

  # now summarize its output
  my $blastn_summary_file = $out_root . ".r1.blastn.summary.txt";
  my $parse_cmd = $execs_HR->{"parse_blast"} . " --program n --input $blastn_out_file > $blastn_summary_file";
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
  # SLEN   ignored
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
      }
      elsif($key eq "HACC") { 
        if(! defined $cur_H{"QACC"}) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read HACC line $line before QACC line (line: $line_idx)\n", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
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
        elsif($value =~ /^(\d+)..(\d+)$/) { 
          ($cur_H{"SRANGESTART"}, $cur_H{"SRANGESTOP"}) = ($1, $2);

          # output data in cmscan --trmF3 format
          if((! defined $cur_H{"QRANGESTART"}) || (! defined $cur_H{"QRANGESTOP"})) { 
            ofile_FAIL("ERROR in $sub_name, reading $blastn_summary_file, read $key line before QRANGE line (line: $line_idx)\n", 1, $FH_HR);
          }
          my $cur_seq_name = $cur_H{"QACC"};
          my $cur_seq_len  = $seq_len_HR->{$cur_seq_name}; # we checked this is defined when we read QACC line
          my $cur_mdl_name = $cur_H{"HACC"};
          if(defined $pretblout_FH) { 
            printf $pretblout_FH ("%-30s  %-30s  %8.1f  %9d  %9d  %6s  %6s  %3s  %11s\n", 
                                  $cur_mdl_name, 
                                  $cur_seq_name,
                                  $cur_H{"BITSCORE"},
                                  $cur_H{"QRANGESTART"}, 
                                  $cur_H{"QRANGESTOP"}, 
                                  $cur_H{"QSTRAND"},
                                  sprintf("    %s%s", 
                                          ($cur_H{"QRANGESTART"} == 1)            ? "[" : ".",
                                          ($cur_H{"QRANGESTOP"}  == $cur_seq_len) ? "]" : "."),
                                  "?",
                                  $cur_seq_len);
            
            # update summed score in %scsum_HHH for this model/seq/strand trio
            if(! defined $scsum_HHH{$cur_H{"HACC"}}) { 
              %{$scsum_HHH{$cur_H{"HACC"}}} = ();
            }
            if(! defined $scsum_HHH{$cur_H{"HACC"}}{$cur_H{"QACC"}}) { 
              %{$scsum_HHH{$cur_H{"HACC"}}{$cur_H{"QACC"}}} = ();
            }
            if(! defined $scsum_HHH{$cur_H{"HACC"}}{$cur_H{"QACC"}}{$cur_H{"QSTRAND"}}) { 
              $scsum_HHH{$cur_H{"HACC"}}{$cur_H{"QACC"}}{$cur_H{"QSTRAND"}} = 0.;
            }
            $scsum_HHH{$cur_H{"HACC"}}{$cur_H{"QACC"}}{$cur_H{"QSTRAND"}} += $cur_H{"BITSCORE"};
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
                              $cur_H{"SRANGESTART"},
                              $cur_H{"SRANGESTOP"},
                              $cur_H{"QRANGESTART"},
                              $cur_H{"QRANGESTOP"},
                              $cur_H{"QSTRAND"},
                              $cur_H{"BITSCORE"},
                              $cur_H{"EVALUE"});
              $cur_FH = $indel_FH_H{$cur_mdl_name};
              printf $cur_FH ("%s\n", $cur_seq_name);
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
