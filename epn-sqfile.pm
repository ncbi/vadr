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
# Subroutine : sqfile_EslSeqstatOptAParse()
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
sub sqfile_EslSeqstatOptAParse { 
  my $nargs_expected = 4;
  my $sub_name = "sqfile_EslSeqstatOptAParse";
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
# Subroutine : sqfile_FastaFileRemoveDescriptions()
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
sub sqfile_FastaFileRemoveDescriptions {
  my $nargs_expected = 3;
  my $sub_name = "sqfile_FastaFileRemoveDescriptions()";
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

#################################################################
# Subroutine: sqfile_GenbankParse()
# Incept:     EPN, Tue Mar 12 14:04:14 2019
#
# Synopsis: Parse a GenBank format file.
#
# Arguments:
#  $infile:   GenBank file to parse
#  $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if we have trouble parsing the file
#
# Reference: https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
#            http://www.insdc.org/documents/feature-table
#################################################################
sub sqfile_GenbankParse { 
  my $sub_name = "sqfile_GenbankParse";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $seq_info_HHR, $ftr_info_HAHR, $FH_HR) = @_;

  my $line_idx  = 0;     # line index of input file
  my $acc       = undef; # accession, read from LOCUS line
  my $tmp_acc   = undef; # accession, read from ACCESSION or VERSION line
  my $len       = undef; # length, read from LOCUS line
  my $def       = undef; # seq definition, read from DEFINITION line
  my $ver       = undef; # sequence version, read from VERSION line
  my $feature   = undef; # a feature   read from a feature/location line in the FEATURES section
  my $location  = undef; # a location  read from a feature/location line in the FEATURES section
  my $qualifier = undef; # a qualifier read from a qualifier/value  line in the FEATURES section
  my $value     = undef; # a value     read from a qualifier/value  line in the FEATURES section
  my $seq       = undef; # sequence, read from the ORIGIN section
  my $seqline   = undef; # single line of sequence
  my $seq_idx   = 0;     # number of sequences read
  my $ftr_idx   = -1;    # number of features read for current sequence
  my $line      = undef; # a line

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  $line = <IN>; 
  while(defined $line) { 
    chomp $line; $line_idx++;
    if($line =~ /^LOCUS\s+(\S+)\s+(\d+)/) { 
      #LOCUS       NC_039477               7567 bp    RNA     linear   VRL 22-FEB-2019
      if((defined $acc) || (defined $len)) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple LOCUS lines for single record ($acc), line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
      ($acc, $len) = ($1, $2);
      # initialize the array of hashes for this accession's features
      if(defined $ftr_info_HAHR->{$acc}) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, trying to add feature info for accession $acc, but it already exists, line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
      @{$ftr_info_HAHR->{$acc}} = ();
      $line = <IN>; 
    }
    elsif($line =~ /^DEFINITION\s+(.*)$/) { 
      #DEFINITION  Norovirus GII isolate strain Hu/GBR/2016/GII.P16-GII.4_Sydney/226,
      #            complete genome.
      if(defined $def) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple DEFINITION lines for single record ($acc), line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
      $def = $1;
      # read remainder of the definition (>= 0 lines)
      $line = <IN>; 
      while((defined $line) && ($line =~ /^\s+(.+)$/)) {
        chomp $line; $line_idx++;
        $def .= $1;
        $line = <IN>; 
      }
    }
    elsif($line =~ /^ACCESSION\s+(\S+)$/) { 
      # ACCESSION   NC_039477
      # verify this matches what we read in the LOCUS line
      $tmp_acc = $1;
      if((! defined $acc) || ($tmp_acc ne $acc)) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, accession mismatch for $tmp_acc, line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
      $line = <IN>;
    }
    elsif($line =~ /^VERSION\s+(\S+)$/) { 
      #VERSION     NC_039477.1
      if(defined $ver) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple VERSION lines for single record ($acc), line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
      # verify this matches what we read in the LOCUS line
      $ver = $1;
      $tmp_acc = $ver;
      dng_StripVersion(\$tmp_acc);
      if((! defined $acc) || ($tmp_acc ne $acc)) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, version/accession mismatch for $tmp_acc, line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
      $line = <IN>;
    }
    elsif($line =~ /^FEATURES\s+Location\/Qualifiers$/) { 
      # parse the features and then the sequence
      # FEATURES section
      # two types of line:
      # feature/location line
      #        example:      gene            5..5104
      #        example:      misc_feature    join(2682..2689,1..2)
      #        example:      misc_feature    join(161990..162784,complement(88222..88806),complement(86666..87448))
      # qualifier/value line type A, first line of a new qualifier
      #        example: /codon_start=1
      #        example: /gene="ORF1"
      # qualifier/value line type B, not the first line of a new qualifier, line 2 to N of a qualifier value
      #        example: QNVIDPWIRNNFVQAPGGEFTVSPRNAPGEILWSAPLGPDLNPYLSHLARMYNGYAGG
      #        example: IPPNGYFRFDSWVNQFYTLAPMGNGTGRRRVV"
      if($ftr_idx != -1) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple FEATURES lines for single record ($acc), line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
      $line = <IN>;
      while((defined $line) && ($line !~ /^ORIGIN/)) { 
        chomp $line; $line_idx++;
        if($line =~ /^\s+\/(\S+)\=(.+)$/) { # first token must start with '/'
          # qualifier/value line type A, examples:
          #  /codon_start=1
          #  /gene="ORF1"
          #  /translation="MKMASNDATVAVACNNNNDKEKSSGEGLFTNMSSTLKKALGARP
          my ($save_qualifier, $save_value) = ($1, $2);
          if(defined $value) { # we are finished with previous value
            sqfile_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
          }
          ($qualifier, $value) = ($save_qualifier, $save_value);
        }
        elsif($line =~ /^\s+(\S+)\s+(\S+)$/) { 
          # NOTE: this will pass for a non-first line of a qualifier value that has whitespace in it:
          # e.g.                      KQP ASRDESQKPPRPPTPELVKRIPPPPPNGEEEEEPVIRYEVKSGISGLPELTTVPQ
          # But I think those are illegal, if they're not, then we'll set "KQP" as feature below, which is bad
          if(defined $value) { # we are finished with previous value
            sqfile_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
            ($qualifier, $value) = (undef, undef);
          }
          # feature/location line, examples:
          #   gene            5..5104
          ($feature, $location) = ($1, $2);
          $ftr_idx++;
          sqfile_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",     $feature,  $FH_HR);
          sqfile_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "location", $location, $FH_HR);
        }
        else { 
          # qualifier/value line type B
          #        example: QNVIDPWIRNNFVQAPGGEFTVSPRNAPGEILWSAPLGPDLNPYLSHLARMYNGYAGG
          #        example: IPPNGYFRFDSWVNQFYTLAPMGNGTGRRRVV"
          $line =~ s/^\s+//; # remove leading whitespace
          if(! defined $value) { 
            ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, in FEATURES section read qualifier value line without qualifier first, line:\n$line\n", "dnaorg", 1, $FH_HR);
          }
          $value .= $line; 
        }
        $line = <IN>; chomp $line; $line_idx++;
      }
      if(! defined $line) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, expected to read ORIGIN line after FEATURES but did not\n", "dnaorg", 1, $FH_HR);
      }
      # if we get here we just read the ORIGIN line
      # first store final qualifier/value
      if(defined $value) { 
        sqfile_GenbankStoreQualifierValue($ftr_info_HAHR->{$acc}, $ftr_idx, $qualifier, $value, $FH_HR);
      }
      # parse the ORIGIN sequence
      $line = <IN>;
      # sanity check
      if(defined $seq) { 
        ofile_FAIL("ERROR in $sub_name, read multiple ORIGIN lines for single record ($acc), line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
      $seq = "";
      while((defined $line) && ($line !~ /^\/\/$/)) { 
        chomp $line; $line_idx++;
        # sequence lines
        # examples:
        # 7501 gtcacgggcg taatgtgaaa agacaaaact gattatcttt ctttttcttt agtgtctttt
        # 7561 aaaaaaa
        if($line =~ /^\s+\d+\s+(.+)$/) { 
          $seqline = $1;
          $seqline =~ s/\s+//g; # remove spaces
          $seq .= $seqline;
        }
        $line = <IN>;
      }
      if(! defined $line) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, expected to find a // line after ORIGIN but did not, line $line_idx\n", "dnaorg", 1, $FH_HR);
      }
      # if we get here we just read the // line
      # we are finished with this sequence, store the information
      if(! defined $acc) { ofile_FAIL(        "ERROR in $sub_name, failed to read accession, line: $line_idx\n", "dnaorg", 1, $FH_HR); }
      if(! defined $len) { ofile_FAIL(sprintf("ERROR in $sub_name, failed to read length (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), "dnaorg", 1, $FH_HR); }
      if(! defined $ver) { ofile_FAIL(sprintf("ERROR in $sub_name, failed to read version (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), "dnaorg", 1, $FH_HR); }
      if(! defined $def) { ofile_FAIL(sprintf("ERROR in $sub_name, failed to read definition (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), "dnaorg", 1, $FH_HR); }
      if(! defined $seq) { ofile_FAIL(sprintf("ERROR in $sub_name, failed to read sequence (accn: %s), line: $line_idx\n", (defined $acc ? $acc : "undef")), "dnaorg", 1, $FH_HR); }

      # store sequence info
      %{$seq_info_HHR->{$acc}} = ();
      $seq_info_HHR->{$acc}{"len"} = $len;
      $seq_info_HHR->{$acc}{"ver"} = $ver;
      $seq_info_HHR->{$acc}{"def"} = $def;
      $seq_info_HHR->{$acc}{"seq"} = $seq;

      # reset variables
      $seq = undef;
      $len = undef;
      $acc = undef;
      $ver = undef;
      $def = undef;
      $feature   = undef;
      $location  = undef;
      $qualifier = undef;
      $value     = undef;

      $seq_idx++;
      $ftr_idx = -1;
      $line = <IN>;
    } # end of 'elsif($line =~ /^FEATURES\s+Location\/Qualifiers$/) {' 
    else { 
      # not a line we will parse, read the next line
      $line = <IN>;
    }
  }

  if($seq_idx == 0) { 
    ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, failed to read any sequence data\n", "dnaorg", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: sqfile_GenbankStoreQualifierValue()
# Incept:     EPN, Wed Mar 13 09:42:22 2019
#
# Synopsis: Store a genbank qualifier and value.
#
# Arguments:
#  $ftr_info_AHR: REF to the array of hashes to store data in
#  $ftr_idx:      feature index
#  $qualifier:    qualifier
#  $value:        qualifier value
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#
# Returns:    '1' if $ftr_info_AHR->[$ftr_idx]{$qualifier} created
#             '0' if $ftr_info_AHR->[$ftr_idx]{$qualifier} exists upon entering function
#
# Dies:       If $value includes the string ":GPSEP:, which we use 
#             to separate multiple qualifier values for the same qualifier.
#             
#################################################################
sub sqfile_GenbankStoreQualifierValue { 
  my $sub_name = "sqfile_GenbankStoreQualifierValue";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($ftr_info_AHR, $ftr_idx, $qualifier, $value, $FH_HR) = @_;

  if($value =~ /\:GPSEP\:/) { 
    ofile_FAIL("ERROR in $sub_name, qualifier value $value includes the special string :GPSEP:, this is not allowed", "dnaorg", 1, $FH_HR);
  }

  # remove leading and trailing " in the value, if they exist
  # GenBank format uses "" as a substitute for " in these strings
  $value =~ s/^\"//;
  $value =~ s/\"$//;

  # printf("in $sub_name q: $qualifier v: $value\n");
  if(! defined ($ftr_info_AHR->[$ftr_idx])) { 
    %{$ftr_info_AHR->[$ftr_idx]} = (); 
  }
  if(! defined $ftr_info_AHR->[$ftr_idx]{$qualifier}) { 
    $ftr_info_AHR->[$ftr_idx]{$qualifier} = $value;
  }
  else { 
    $ftr_info_AHR->[$ftr_idx]{$qualifier} .= ":GBSEP:" . $value;
  }

  return;
}

#################################################################
# Subroutine: dng_CdsFetchStockholmToFasta()
# Incept:     EPN, Thu Mar 14 12:30:33 2019
# 
# Purpose:    Given coordinates of all CDS features in %{$ftr_info_AHR}
#             fetch all the CDS for all sequences in the Stockholm alignment
#             and create a new output fasta file with just the CDS features.
#
#             We don't really need both the stockholm and fasta file 
#             if there are no gaps in the stockholm alignment (as is the
#             case in dnaorg_build.pl (which requires a single sequence 
#             'alignment' with no gaps), but this implemenation works for
#             alignments with gaps too.
#
# Arguments:
#   $out_FH:         output file handle
#   $stk_file:       stockholm file with aligned full length sequences
#   $ftr_info_AHR:   REF to the feature info, pre-filled
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#                    
# Returns: void
#
# Dies:    if we have trouble fetching a sequence
#
#################################################################
sub dng_CdsFetchStockholmToFasta { 
  my $sub_name = "dng_CdsFetchStockholmToFasta";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $stk_file, $ftr_info_AHR, $FH_HR) = @_;

  my $msa = Bio::Easel::MSA->new({ fileLocation => $stk_file, isDna => 1});
  my $msa_has_rf = $msa->has_rf;

  # precompute start, stop, strand, for all features, so we don't have to redo this for each seq
  my @sgm_start_AA  = ();
  my @sgm_stop_AA   = ();
  my @sgm_strand_AA = ();
  dng_FeatureInfoStartStopStrandArrays($ftr_info_AHR, \@sgm_start_AA, \@sgm_stop_AA, \@sgm_strand_AA, $FH_HR);

  my $nftr = scalar(@{$ftr_info_AHR});
  my $nseq = $msa->nseq;
  my $ftr_idx = undef; # feature index
  my $seq_idx = undef; # feature index
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") { 
        my $cds_sqstring = "";
        foreach(my $sgm_idx = 0; $sgm_idx < scalar(@{$sgm_start_AA[$ftr_idx]}); $sgm_idx++) { 
          my $rfstart = $sgm_start_AA[$ftr_idx][$sgm_idx];
          my $rfstop  = $sgm_stop_AA[$ftr_idx][$sgm_idx];
          my $astart  = ($msa_has_rf) ? $msa->rfpos_to_aligned_pos($rfstart) : $rfstart;
          my $astop   = ($msa_has_rf) ? $msa->rfpos_to_aligned_pos($rfstop)  : $rfstop;
          my $sgm_sqstring = $msa->get_sqstring_unaligned_and_truncated($seq_idx, $astart, $astop);
          if($sgm_strand_AA[$ftr_idx][$sgm_idx] eq "-") { 
            sqstringReverseComplement(\$sgm_sqstring);
          }
          $cds_sqstring .= $sgm_sqstring;
        }
        print $out_FH(">" . $msa->get_sqname($seq_idx) . "/" . $ftr_info_AHR->[$ftr_idx]{"coords"} . "\n" . seq_SqstringAddNewlines($cds_sqstring, 60));
      }
    }
  }
  return;
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
