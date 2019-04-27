#!/usr/bin/perl
#
# epn-seqfile.pm
# Eric Nawrocki
# EPN, Wed Apr  3 06:13:49 2019
# version: 0.00
#
use strict;
use warnings;
use Time::HiRes qw(gettimeofday);

require "epn-ofile.pm";
require "epn-options.pm";
require "epn-seq.pm";
require "epn-utils.pm";

#####################################################################
# Data structures used in this module:
#
#####################################################################
#
# List of subroutines:
# 
#################################################################
# Subroutine : sqf_EslSeqstatOptAParse()
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
sub sqf_EslSeqstatOptAParse { 
  my $nargs_expected = 4;
  my $sub_name = "sqf_EslSeqstatOptAParse";
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
# Subroutine : sqf_FastaFileRemoveDescriptions()
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
sub sqf_FastaFileRemoveDescriptions {
  my $nargs_expected = 3;
  my $sub_name = "sqf_FastaFileRemoveDescriptions()";
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
# Subroutine: sqf_GenbankParse()
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
sub sqf_GenbankParse { 
  my $sub_name = "sqf_GenbankParse";
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
      # 4 types of line:
      # feature/location line
      #        example1:      gene            5..5104
      #        example2:      misc_feature    join(2682..2689,1..2)
      #        example3:      misc_feature    join(161990..162784,complement(88222..88806),complement(86666..87448))
      # qualifier/value type A, single line of controlled vocab or enumerated value
      #        example: /codon_start=1
      # qualifier/value type B, single line of free text
      #        example1: /gene="ORF1"
      # qualifier/value type C, multiple lines of free text
      #       example1: /translation="MMMASKDVVPTAASSENANNNSSIKSRLLARLKGSGGATSPPNS
      #       example1: IKITNQDMALGLIGQVPAPKATSVDVPKQQRDRPPRTVAEVQQNLRWTERPQDQNVKT
      #       example1: QNVIDPWIRNNFVQAPGGEFTVSPRNAPGEILWSAPLGPDLNPYLSHLARMYNGYAGG
      #       example1: IPPNGYFRFDSWVNQFYTLAPMGNGTGRRRVV"
      #
      #       example2: /note="The mature peptides were added by the NCBI staff
      #       example2: following publications Liu et all (1996, 1999), that
      #       example2: (tentatively) determined processing map by site-directed
      #       example2: mutagenesis or by direct sequencing of the cleavage
      #       example2: products for closely related Southampton virus (a
      #       example2: Norwalk-like virus); ORF1; sequence homologies to 2C
      #       example2: helicase, 3C protease, and 3D RNA-dependent RNA polymerase
      #       example2: of picornavirus"
      #
      if($ftr_idx != -1) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read multiple FEATURES lines for single record ($acc), line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
      $line = <IN>;
      while((defined $line) && ($line !~ /^ORIGIN/)) { 
        chomp $line; $line_idx++;
        if($line =~ /^\s+\/(\S+)\=(.+)$/) { # first token must start with '/'
          my ($save_qualifier, $save_value) = ($1, $2);
          # determine type of qualifier value
          # A: single line of controlled vocabulary or enumerated value, does not start with '"'
          # B: single line of free text: starts and ends with single '"'
          # C: first line of multiple lines of free text: starts with single '"', does not end with single '"'
          if($save_value  =~ /^[^\"].*/) { 
            # A: single line of controlled vocabulary or enumerated value, does not start with '"'
            # example:
            # /codon_start=1
            ; # do nothing, we just need this if to catch a valid line (so we don't get to the else below that leads to a failure)
          }
          elsif(($save_value =~ /^\"[^\"]\"$/) || # single non-'"' character in between two '"' characters
                ($save_value =~ /^\"[^\"].*[^\"]\"$/)) { # '"' then a non-'"', then >= 0 characters, then a non-'"' character then ends with a '"' character
            # B: single line of free text: starts and ends with single '"'
            # example:
            # /product="nonstructural polyprotein"
            ; # do nothing, we just need this elsif to catch a valid line (so we don't get to the else below that leads to a failure)
          }
          elsif(($save_value =~ /^\"[^\"]$/) || # single non-'"' character after '"' character
                ($save_value =~ /^\"[^\"].*[^\"]$/)) { # '"' then a non-'"', then >= 0 characters, then ends with a non-'"' character
            # C: first line of multiple lines of free text: starts with single '"', does not end with single '"'
            # example:
            # /note="The mature peptides were added by the NCBI staff
            # 
            # continue to read lines until we read one that ends with a single '"'
            #
            # example of remaining lines in value:
            #  following publications Liu et all (1996, 1999), that
            #  (tentatively) determined processing map by site-directed
            #  mutagenesis or by direct sequencing of the cleavage
            #  products for closely related Southampton virus (a
            #  Norwalk-like virus); ORF1; sequence homologies to 2C
            #  helicase, 3C protease, and 3D RNA-dependent RNA polymerase
            #  of picornavirus"
            $line = <IN>; chomp $line; $line_idx++;
            while((defined $line) && ($line !~ /[^\"]\"$/)) { # ends with a non-'"' character followed by '"' character
              $line =~ s/^\s+//; # remove leading whitespace
              $save_value .= $line;
              $line = <IN>; chomp $line; $line_idx++;
            }
            # add final line we read
            $line =~ s/^\s+//; # remove leading whitespace
            $save_value .= $line;
          }
          else { 
            ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, in FEATURES section read unparseable qualifier value line:\n$line\n", "dnaorg", 1, $FH_HR);
          }
          if(defined $value) { # we are finished with previous value
            sqf_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
          }
          ($qualifier, $value) = ($save_qualifier, $save_value);
        }
        elsif($line =~ /^\s+(\S+)\s+(\S+)$/) { 
          if(defined $value) { # we are finished with previous value
            sqf_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
            ($qualifier, $value) = (undef, undef);
          }
          # feature/location line, examples:
          #   gene            5..5104
          ($feature, $location) = ($1, $2);
          $ftr_idx++;
          sqf_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",     $feature,  $FH_HR);
          sqf_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "location", $location, $FH_HR);
        }
        $line = <IN>;
      }
      if(! defined $line) { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, expected to read ORIGIN line after FEATURES but did not\n", "dnaorg", 1, $FH_HR);
      }
      # if we get here we just read the ORIGIN line
      # first store final qualifier/value
      if(defined $value) { 
        sqf_GenbankStoreQualifierValue($ftr_info_HAHR->{$acc}, $ftr_idx, $qualifier, $value, $FH_HR);
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
# Subroutine: sqf_GenbankStoreQualifierValue()
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
sub sqf_GenbankStoreQualifierValue { 
  my $sub_name = "sqf_GenbankStoreQualifierValue";
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
# Subroutine: sqf_BlastDbProteinCreate
# Incept:     EPN, Mon Mar 18 09:40:28 2019
# 
# Purpose:    Create a protein blast database from a fasta file.
#
# Arguments:
#   $makeblastdb:    path to 'makeblastdb' executable
#   $fa_file:        FASTA file of protein sequences to make blast db from
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#                    
# Returns:    void
#
#################################################################
sub sqf_BlastDbProteinCreate {
  my $sub_name = "sqf_BlastDbProteinCreate";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($makeblastdb, $fa_file, $opt_HHR, $FH_HR) = @_;

  utl_RunCommand($makeblastdb . " -in $fa_file -dbtype prot > /dev/null", opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine: sqf_EslTranslateCdsToFastaFile()
# Incept:     EPN, Thu Mar 14 12:30:28 2019
# 
# Purpose:    Use esl-translate to translate a fasta file with
#             CDS sequences pertaining to the CDS features in 
#             @{$ftr_info_AHR} into fasta protein files.
#
# Arguments:
#   $out_FH:         output file handle to print to 
#   $esl_translate:  path to esl-translate executable
#   $cds_fa_file:    fasta file with CDS sequences
#   $out_root:       string that is the 'root' for naming output files
#   $ftr_info_AHR:   REF to the feature info, pre-filled
#   $opt_HHR:        command line options
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#                    
# Returns: void
#
# Dies:    if we have trouble fetching a sequence
#
#################################################################
sub sqf_EslTranslateCdsToFastaFile { 
  my $sub_name = "sqf_EslTranslateCdsTranslateToFastaFile";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $esl_translate, $cds_fa_file, $out_root, $ftr_info_AHR, $opt_HHR, $FH_HR) = @_;

  my $tmp1_translate_fa_file  = $out_root . ".cds.esl-translate.1.fa";
  my $tmp2_translate_fa_file  = $out_root . ".cds.esl-translate.2.fa";
  my $tmp1_translate_ssi_file = $out_root . ".cds.esl-translate.1.fa.ssi";
  my $tmp2_translate_ssi_file = $out_root . ".cds.esl-translate.2.fa.ssi";
  if(-e $tmp1_translate_ssi_file) { unlink $tmp1_translate_ssi_file; }
  if(-e $tmp2_translate_ssi_file) { unlink $tmp2_translate_ssi_file; }

  my $c_opt = "";
  if((opt_IsUsed("--ttbl", $opt_HHR)) && (opt_Get("--ttbl", $opt_HHR) != 1)) { 
    $c_opt = "-c " . opt_Get("--ttbl", $opt_HHR);
  }
  my $translate_cmd = "$esl_translate $c_opt -M -l 3 --watson $cds_fa_file > $tmp1_translate_fa_file";
  utl_RunCommand($translate_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # go through output fasta file and rewrite names, so we can fetch 
  open(IN,       $tmp1_translate_fa_file) || fileOpenFailure($tmp1_translate_fa_file, $sub_name, $!, "reading", $FH_HR);
  open(OUT, ">", $tmp2_translate_fa_file) || fileOpenFailure($tmp2_translate_fa_file, $sub_name, $!, "writing", $FH_HR);
  while(my $line = <IN>) { 
    if($line =~ m/^\>/) { 
      #>orf58 source=NC_039477.1/5..5104:+ coords=1..5097 length=1699 frame=1  
      chomp $line;
      if($line =~ /^\>orf\d+\s+(source\=\S+)\s+(coords\=\S+)\s+length\=\d+\s+frame\=\S+/) { 
        # rename as 'source=NC_039477.1/5..5104:+,coords=1..5097'
        print OUT (">" . $1 . "," . $2 . "\n");
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, problem parsing esl-translate output file $tmp1_translate_fa_file, line:\n$line\n", "dnaorg", 1, $FH_HR);
      }
    }
    else { 
      print OUT $line; 
    }
  }
  close(IN);
  close(OUT);

  # $tmp2_translate_fa_file now includes renamed translated sequences from esl-translate
  # fetch expected translated seqs and print to $out_FH
  my $cds_sqfile     = Bio::Easel::SqFile->new({ fileLocation => $cds_fa_file });
  my $protein_sqfile = Bio::Easel::SqFile->new({ fileLocation => $tmp2_translate_fa_file });

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $seq_idx = 0; $seq_idx < $cds_sqfile->nseq_ssi; $seq_idx++) { 
    my ($seq_name, $seq_length) = $cds_sqfile->fetch_seq_name_and_length_given_ssi_number($seq_idx);
    my $fetch_name = "source=" . $seq_name . ",coords=1.." . ($seq_length - 3); # subtract length of stop codon
    if(! $protein_sqfile->check_seq_exists($fetch_name)) { 
      ofile_FAIL("ERROR in $sub_name, problem translating CDS feature, unable to find expected translated sequence in $tmp2_translate_fa_file:\n\tseq: $seq_name\n\texpected sequence:$fetch_name\n", "dnaorg", 1, $FH_HR);
    }
    print $out_FH ">" . $seq_name . "\n";
    print $out_FH seq_SqstringAddNewlines($protein_sqfile->fetch_seq_to_sqstring($fetch_name), 60);
  }
  # remove temporary files unless --keep
  if(! opt_Get("--keep", $opt_HHR)) { 
    utl_FileRemoveUsingSystemRm($tmp1_translate_fa_file, $sub_name, $opt_HHR, $FH_HR);
    utl_FileRemoveUsingSystemRm($tmp2_translate_fa_file, $sub_name, $opt_HHR, $FH_HR);
    utl_FileRemoveUsingSystemRm($tmp2_translate_fa_file . ".ssi", $sub_name, $opt_HHR, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: sqf_FastaWriteSequence()
# Incept:     EPN, Thu Mar 14 06:06:59 2019
#
# Synopsis: Print a sequence to a fasta file.
#
# Arguments:
#  $out_FH:    output file handle
#  $name:      sequence name
#  $def:       sequence definition, can be undef
#  $seq:       sequence string
#  $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $name or $seq is undef
#################################################################
sub sqf_FastaWriteSequence {
  my $sub_name = "sqf_FastaWriteSequence";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $name, $def, $seq, $FH_HR) = @_;

  if(! defined $name) { ofile_FAIL("ERROR in $sub_name, name is undefined", "dnaorg", 1, $FH_HR); }
  if(! defined $seq)  { ofile_FAIL("ERROR in $sub_name, name is undefined", "dnaorg", 1, $FH_HR); }

  # capitalize and DNAize $seq
  seq_SqstringCapitalize(\$seq);
  seq_SqstringDnaize(\$seq);
  printf $out_FH (">%s%s\n%s", 
                  $name, 
                  (defined $def) ? " " . $def : "",
                  seq_SqstringAddNewlines($seq, 60));
  
  return;
}

#################################################################
# Subroutine: sqf_EslReformatRun()
# Incept:     EPN, Fri Mar 15 12:56:04 2019
#
# Synopsis: Use esl-reformat to convert a file from one format to another
#
# Arguments:
#  $esl_reformat: esl-reformat executable file
#  $in_file:      input file for esl-reformat
#  $out_file:     output file from esl-reformat
#  $in_format:    input format
#  $out_format:   output format
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's a problem fetching the sequence file
#################################################################
sub sqf_EslReformatRun { 
  my $sub_name = "sqf_EslReformatRun";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_reformat, $in_file, $out_file, $in_format, $out_format, $opt_HHR, $FH_HR) = @_;

  my $cmd = $esl_reformat . " --informat $in_format $out_format $in_file > $out_file";
  utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # remove a .ssi file for newly created file if it exists
  my $ssi_file = $out_file . ".ssi";
  if(-e $ssi_file) { unlink $ssi_file; }

  return;
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
