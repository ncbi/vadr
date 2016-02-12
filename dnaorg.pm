#!/usr/bin/perl
#
# dnaorg.pm
# Eric Nawrocki
# EPN, Mon Feb  1 15:23:00 2016
# 
# Some subroutines are taken from rnavore.pm [EPN, Tue Sep 23 09:22:55 2014]
# and ssu.pm [EPN, Thu Nov  5 05:39:37 2009].
#
# Perl module used by dnaorg_build.pl and dnaorg_annotate.pl scripts
# which contains subroutines called by both of those scripts.
#
#########################
# Naming conventions used in this file:
# - Hash  data structures use the abbreviation 'H', typically at the end of 
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
# - TODO: FH_HR
########################################################################################
#
# List of subroutines in this file:
#
# runCommand():                            run and time a command and die if it fails
# parseMatPeptSpecFile():                  parse input file specifying relationship between CDS and mature peptides
# parseLengthFile():                       parse a file with accessions and lengths into a hash
# parseEdirectFtableFile():                parse a ftable file output from 'efetch' into usable data structures
# parseEdirectMatPeptideFile():            parse a file output from 'efetch' with mature peptide info into usable data structures
# helperBreakdownFac():                    helps breakdown a specially formatted '$fac' string (accession and coordinates)
# startsStopsFromCoords():                 extract start and stop positions from a coordinates string
# stripVersion():                          remove version from an accession.version string
# fetchReferenceFeatureSequences():        fetch sequences of features for a reference sequence
# mapFeaturesToModels():                   fill a 2D array that maps each feature to the models that model it 
# determineFeatureTypes():                 determine what type each feature is (e.g. "mp", "cds-notmp", "cds-mp")
# getNumFeaturesWithModels():              determine number of features that are modelled
# startsStopsStrandsFromCoordsLength():    determine start, stop and strands from a coordinate string and length
# dumpInfoHashOfArrays():                  print out all values in an 'info' hash of arrays (e.g. %ftr_info_HA or %mdl_info_HA)
# validateAndGetSizeOfInfoHashOfArrays():  validate that an 'info' hash of arrays is valid (all arrays are same size)
# getLengthsAndCoords():                   fill arrays with lengths and coordinate strings of features
# lengthFromCoords():                      determine length of a feature from its coordinate string
# fetchSequencesUsingEslFetchCds():        fetch sequences using the esl_fetch_cds.pl script
# addNameAndBlankSsToStockholmAlignment(): add name and blank SS_cons annotation to a Stockholm alignment
# getQualifierValues():                    extract qualifier values from a hash of hash of arrays
# createCmDb():                            create a CM database from a stockholm alignment file
# matpeptValidateCdsRelationships():       validate relationships between CDS and the mature peptides that comprise them
# checkForSpanningSequenceSegments():      check if two sequence segments span the stop..start boundary in a circular genome
# edirectFtableOrMatPept2SingleFeatureTableInfo(): extract 'feature table information' from a ftable or matpept edirect file
# getSingleFeatureTableInfo():             get single feature table information from a 'qualifiers' hash of arrays
# outputBanner():                          output the program banner, commands, date and options used to a file handle
# getStrandStats():                        return stats on strands for a given accession
# DNAORG_FAIL():                           exit with an error, after printing that error to the sum and log files and closing those files
# fileOpenFailure():                       exit with an error related to opening a file for reading or writing
use strict;
use warnings;

#################################################################
# Subroutine:  runCommand()
# Incept:      EPN, Thu Feb 11 13:32:34 2016
#
# Purpose:     Runs a command using system() and exits in error 
#              if the command fails. If $be_verbose, outputs
#              the command to stdout. If $FH_HR->{"cmd"} is
#              defined, outputs command to that file handle.
#
# Arguments:
#   $cmd:         command to run, with a "system" command;
#   $be_verbose:  '1' to output command to stdout before we run it, '0' not to
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails
#################################################################
sub runCommand {
  my $sub_name = "runCommand()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $be_verbose, $FH_HR) = @_;
  
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

  if($? != 0) { 
    DNAORG_FAIL("ERROR in $sub_name, the following command failed:\n$cmd\n", $?, $FH_HR); 
  }

  return ($stop_time - $start_time);
}

#################################################################
# Subroutine:  parseMatPeptSpecFile()
# Incept:      EPN, Thu Feb 11 13:21:02 2016
#
# Purpose:     Parse the input specifications file that defines the
#              relationship between each CDS and the mature peptides.
#
# Arguments:
#   $infile:           file to parse
#   $cds2pmatpept_AAR: ref to array of arrays to fill here, CDS: 'primary' mat_peptide relationships
#   $cds2amatpept_AAR: ref to array of arrays to fill here, CDS: 'all'     mat_peptide relationships
#   $FH_HR:            REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if problem reading $infile (in unexpected format)
#
#################################################################
sub parseMatPeptSpecFile {
  my $sub_name = "parseMatPeptSpecFile";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $cds2pmatpept_AAR, $cds2amatpept_AAR, $FH_HR) = @_;
  my $ncds_read_primary = 0;
  my $ncds_read_all     = 0;
  my $cds_idx2store     = 0;
  my @cds_idx_read_A    = (); # $cds_idx_read_A[$i] = 1 if we read info for CDS $i+1
  my $max_cds_idx2store = 0;

  open(IN, $infile) || fileOpenFailure($infile, $FH_HR, $!, "reading");

  while(my $line = <IN>) { 
    if($line !~ m/^\#/) { 
      # example input file:
      ## This file explains how CDS and mat_peptide annotation for NC_001477
      ## are related.
      ##
      ############beg of file
      ## Format of lines in this file:
      ## <CDS-idx> <'primary' OR 'all'> <mat_peptide-1-idx>:<mat_peptide-2-idx>:<mat_peptide-n-idx>
      ## 'primary' lines: these define the 'primary' peptides in order, the
      ##                  CDS <CDS-idx> is comprised of the peptides listed
      ##                  in final token, which are contiguous, start of 
      ##                  first mat_peptide to stop of final mat_peptide is
      ##                  one contiguous subsequence.
      ##
      ## 'all' lines:     these define all the peptides that are ultimately
      ##                  derived from CDS <CDS-idx>. It will be a superset
      ##                  of the primary line for this index but will
      ##                  additionally include mat_peptides that are
      ##                  secondarily cleaved from the primary mat_peptides.
      ##
      #1 primary 1:3:6:7:8:9:10:11:12:13:14
      #1 all     1:2:3:4:5:6:7:8:9:10:11:12:13:14
      ## 
      ############end of file
      # NOTE: in the input file CDS and matpept indices are in coordinate space 1..N, but we store them in 0..N-1
      # 
      # we need to have one 'primary' and one 'all' line for each  
      chomp $line;
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) != 3) { 
        DNAORG_FAIL("ERROR in $sub_name, unable to parse matpept input file line: $line", 1, $FH_HR); 
      }
      my ($cds_idx, $primary_or_all, $matpept_str) = ($el_A[0], $el_A[1], $el_A[2]);
      $primary_or_all =~ tr/A-Z/a-z/;
      if($primary_or_all ne "all" && $primary_or_all ne "primary") { 
        DNAORG_FAIL("ERROR in $sub_name, parsing $infile, second token of each non-comment line should be 'all' or 'primary'", 1, $FH_HR);
      }
      my $cds_idx2store = $cds_idx - 1;
      if($cds_idx2store < 0) { 
        DNAORG_FAIL("ERROR in $sub_name, read CDS idx that is 0 or less ($cds_idx) in matpept input file", 1, $FH_HR); 
      }
      $cds_idx_read_A[$cds_idx2store] = 1;
      if($cds_idx2store > $max_cds_idx2store) { 
        $max_cds_idx2store = $cds_idx2store; 
      }
      if($primary_or_all eq "primary") { 
        if(defined $cds2pmatpept_AAR->[$cds_idx2store] || exists $cds2pmatpept_AAR->[$cds_idx2store]) {
          DNAORG_FAIL("ERROR in $sub_name, two primary lines for same CDS idx ($cds_idx) in matpept input file", 1, $FH_HR);
        }
        my @matpept_A = split(":", $matpept_str);
        @{$cds2pmatpept_AAR->[$cds_idx2store]} = ();
        foreach my $mp (@matpept_A) { 
          push(@{$cds2pmatpept_AAR->[$cds_idx2store]}, ($mp-1));
        }
        $ncds_read_primary++;
      }
      elsif($primary_or_all eq "all") { 
        if(defined $cds2amatpept_AAR->[$cds_idx2store] || exists $cds2amatpept_AAR->[$cds_idx2store]) {
          DNAORG_FAIL("ERROR in $sub_name, two all lines for same CDS idx ($cds_idx) in matpept input file", 1, $FH_HR);
        }
        my @matpept_A = split(":", $matpept_str);
        @{$cds2amatpept_AAR->[$cds_idx2store]} = ();
        foreach my $mp (@matpept_A) { 
          push(@{$cds2amatpept_AAR->[$cds_idx2store]}, ($mp-1));
        }
        $ncds_read_all++;
      }
    }
  }
  close(IN);

  # three sanity checks:
  # 1: we should have stored all and primary info for any CDS $i for which $cds_idx_read_A[$i-1] is 1, and nothing else
  for(my $i = 0; $i <= $max_cds_idx2store; $i++) { 
    if($cds_idx_read_A[$i]) { 
     if((! defined $cds2pmatpept_AAR->[$i]) && (! exists $cds2pmatpept_AAR->[$i])) { 
       DNAORG_FAIL(sprintf("ERROR in $sub_name, did not properly read info for cds %d in $infile", $i+1), 1, $FH_HR);
     }
     if((! defined $cds2amatpept_AAR->[$i]) && (! exists $cds2amatpept_AAR->[$i])) { 
       DNAORG_FAIL(sprintf("ERROR in $sub_name, did not properly read info for cds %d in $infile\n", $i+1), 1, $FH_HR);
     }
    }
    else { # we didn't read this one
      if(defined $cds2pmatpept_AAR->[$i] || exists $cds2pmatpept_AAR->[$i]) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, improperly read non-existent info for cds %d in $infile\n", $i+1), 1, $FH_HR);
      }
      if(defined $cds2amatpept_AAR->[$i] || exists $cds2amatpept_AAR->[$i]) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, improperly read non-existent info for cds %d in $infile\n", $i+1), 1, $FH_HR);
      }
    }
  }
  # 2: we should have at least read at least one of each 'primary' and 'all'
  if($ncds_read_primary == 0) { 
    DNAORG_FAIL("ERROR in $sub_name, no primary CDS:mat_peptide relationships read in matpept input file $infile", 1, $FH_HR); 
  }
  if($ncds_read_all == 0) { 
    DNAORG_FAIL("ERROR in $sub_name, no all CDS:mat_peptide relationships read in matpept input file $infile", 1, $FH_HR); 
  }
  # 3: all info should be a superset of primary info
  for(my $i = 0; $i <= $max_cds_idx2store; $i++) { 
    if($cds_idx_read_A[$i]) { 
      foreach my $mp (@{$cds2pmatpept_AAR->[$i]}) { 
        # make sure this primary peptide exists in all
        my $found_it = 0;
        foreach my $mp2 (@{$cds2amatpept_AAR->[$i]}) { 
          if($mp == $mp2) { $found_it = 1; }
        }
        if(! $found_it) { 
          DNAORG_FAIL(sprintf("ERROR in $sub_name, all information is not a superset of primary information: %d is in primary but not all", $mp+1), 1, $FH_HR);
        }
      }
    }
  }

  return;
}

#################################################################
# Subroutine:  parseLengthFile()
# Incept:      EPN, Thu Feb 11 13:29:31 2016 
#
# Purpose:     Parses a length file and stores the lengths read
#              into %{$len_HR}. Each line has two tokens separated
#              by whitespace: <accession> <length>
#
# Arguments: 
#   $lenfile: full path to a length file
#   $len_HR:  ref to hash of lengths, key is accession
#   $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's problem parsing $lenfile
#################################################################
sub parseLengthFile {
  my $sub_name = "parseLengthFile()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($lenfile, $len_HR, $FH_HR) = @_;

  open(LEN, $lenfile) || fileOpenFailure($lenfile, $FH_HR, $!, "reading");

  while(my $line = <LEN>) { 
    # example line
    #HM448898.1	2751

    chomp $line;
    my ($accn, $length) = split(/\s+/, $line);
    if($length !~ m/^\d+$/) { 
      DNAORG_FAIL("ERROR in $sub_name, couldn't parse length file line: $line", 1, $FH_HR); 
    } 

    stripVersion(\$accn);
    $len_HR->{$accn} = $length;
  }
  close(LEN);

  return;
}

#################################################################
# Subroutine: parseEdirectFtableFile()
# Incept:     EPN, Thu Feb 11 13:47:11 2016
# 
# Purpose:   Given a feature table file output from an 'efetch -format
#            ft' command, parse that file into usable data
#            structures. 
# 
#            Can be called by the wrapper subroutine: 
#             edirectFtableOrMatPept2SingleFeatureTableInfo().
#
#            Caller will commonly call outputSingleFeatureTable() 
#            or getSingleFeatureTableInfo() after calling this
#            subroutine.
#
#            Similar to (and shares some code with) parseEdirectMatPeptideFile().
#
# Args:
#   $infile:          name of ftable file
#   $only_feature:    name of only feature we will parse, 'undef' to parse all features
#   $sep_HR:          ref to hash of 'separation' values, keys: "qnqv", "fac", and "qval", 
#                     used to split values that have been concatenated with these 'separation'
#                     values in between, PRE-FILLED
#   $quals_HHAR:      ref to 2D hash of arrays with values that could eventually go 
#                     into a single-feature table, FILLED HERE
#                     key 1: feature name, e.g. "CDS"
#                     key 2: 'fac' string: <full_accession><fac_seq><coordinates>
#                     value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>
#   $faccn_AR:        REF to array of full accessions read, FILLED HERE
#   $fac_HHAR:        REF to 2D hash of arrays that is used to easily determine
#                     list of 2D keys ('fac's) in quals_HHA for a given feature and faccn
#                     FILLED HERE
#                     key 1: feature name: e.g. "CDS"
#                     key 2: 'faccn', full accession
#                     value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
#   $faccn2accn_HR:   REF to hash, key: full accession, value: short accession, 
#                     used for convenience in output function, FILLED HERE
#   $column_HAR:      REF to hash of arrays, key: feature, e.g. "CDS"
#                     value: array of qualifiers that exist for that feature in at least 1 faccn
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
#
# Returns:  void
#
# Dies:     if <$infile> file is in unexpected format
#################################################################
sub parseEdirectFtableFile {
  my $sub_name = "parseEdirectFtableFile()";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $only_feature, $sep_HR, $quals_HHAR, $faccn_AR, $fac_HHAR, $faccn2accn_HR, $column_HAR, $FH_HR) = @_;

  open(IN, $infile) || fileOpenFailure($infile, $FH_HR, $!, "reading");

  my $do_only_feature  = 0;
  my $cap_only_feature = undef; # capitalized version of $only_feature
  if(defined $only_feature) { 
    $do_only_feature = 1;
    $cap_only_feature = $only_feature;
    $cap_only_feature =~ tr/a-z/A-Z/;
  }

  # variables used to make sure input feature table is in 
  # expected format w.r.t to order of line types

  # define other variables we'll use
  my $faccn           = undef;   # full accession, e.g. ref|NC_000883.2|
  my $fac             = undef;   # full accession and coords, accession and each segment coords separated by $fac_sep
  my $faccn2          = undef;   # another full accession
  my $accn            = undef;   # full accession with extraneous characters removed: e.g. NC_000883.2
  my $start           = undef;   # start position
  my $stop            = undef;   # stop position
  my $qname           = undef;   # a qualifier name,  e.g. 'organism'
  my $qval            = undef;   # a qualifier value, e.g. 'Paramecia bursaria Chlorella virus 1'
  my $feature         = undef;   # a feature name, e.g. "CDS", or "gene"
  my $cap_feature     = undef;   # a feature name capitalized, necessary so -f option can be case insensitive
  my $coords          = undef;   # coordinates
  my $fac_sep  = $sep_HR->{"fac"};
  my $qnqv_sep = $sep_HR->{"qnqv"};
  my $qval_sep = $sep_HR->{"qval"};
  my $line_ctr        = 0;       # count of number of lines read in ftable
  my $accn_ctr        = 0;       # count of number of accessions read
  my %column_HH       = ();      # existence 2D hash, used to aide construction of %column_HA
                                 # key 1: feature, e.g. "CDS"
                                 # key 2: qualifier
                                 # value: '1' (always)
  my $prv_was_accn           = 0; # set to '1' if previous line was an accession line
  my $prv_was_coords_feature = 0; # set to '1' if previous line was a coordinates line with a feature name
  my $prv_was_coords_only    = 0; # set to '1' if previous line was a coordinates line without a feature name
  my $prv_was_quals          = 0; # set to '1' if previous line was a qualifier_name qualifier value line

  while(my $line = <IN>) { 
    $line_ctr++;
    chomp $line;
    if($line =~ m/\w/) { 
      # parse each of the 4 line types differently
      # -------------------------------------------------------
      if($line =~ /Feature\s+(\S+)$/) { 
        # ACCESSION LINE
        # example:
        #>Feature ref|NC_001359.1|    
        $faccn = $1;
        $accn_ctr++;
        # does line order make sense?
        if($accn_ctr == 1   ||          # first accession of the file
           $prv_was_quals ||          # previous line was a quals line (common)
           $prv_was_coords_feature || # previous line was a coords feature line (rare)
           $prv_was_coords_only    || # previous line was a coords line without a feature (rare)
           $prv_was_accn) {           # previous line was an accession line (rare)
          # line order makes sense, keep going...

          # determine short version of the accession, e.g. NC_001359 in above example
          if($faccn =~ /[^\|]*\|([^\|]*)\|/) { 
            $accn = $1;
            $faccn2accn_HR->{$faccn} = $accn;
            push(@{$faccn_AR}, $faccn);
          }
          else { 
            DNAORG_FAIL("ERROR in $sub_name, unable to parse Feature line $line", 1, $FH_HR); 
          }
          $feature     = undef; 
          $cap_feature = undef;
          $fac         = undef;

          # update '$prv_*' values that we use to make sure line order makes sense
          $prv_was_accn           = 1;
          $prv_was_coords_feature = 0;
          $prv_was_coords_only    = 0;
          $prv_was_quals          = 0;
          #printf("set prv_was_accn\n");
        }
        else { # line order is unexpected
          DNAORG_FAIL("ERROR in $sub_name, unexpected line order (accession line) at line $line_ctr: $line", 1, $FH_HR);
        }
      }
      # -------------------------------------------------------
      elsif($line =~ /^(\<?\d+\^?)\s+(\>?\d+)\s+(\S+)$/) { 
        # COORDINATES LINE WITH A FEATURE NAME
        # example:
        #230	985	gene
        ($start, $stop, $feature) = ($1, $2, $3);

        # does line order make sense?
        if($prv_was_accn  ||           # previous line was accession line (common)
           $prv_was_quals ||           # previous line was quals line (common)
           $prv_was_coords_feature ||  # previous line was coords line with a feature (rare)
           $prv_was_coords_only)    {  # previous line was coords line without a feature (rarer)
          # line order makes sense, keep going...
          $cap_feature = $feature;
          $cap_feature =~ tr/a-z/A-Z/;
          $fac = $faccn . $fac_sep . $start . $fac_sep . $stop;

          # update '$prv_*' values that we use to make sure line order makes sense
          $prv_was_accn           = 0;
          $prv_was_coords_feature = 1;
          $prv_was_coords_only    = 0;
          $prv_was_quals          = 0;
          #printf("set prv_was_coords_feature\n");
        }
        else { # line order is unexpected
          DNAORG_FAIL("ERROR in $sub_name, unexpected line order (coords_feature) at line $line_ctr: $line", 1, $FH_HR);
        }
      }
      # -------------------------------------------------------
      elsif($line =~ /^(\<?\d+)\s+(\>?\d+)$/) { 
        # COORDINATES LINE WITHOUT A FEATURE NAME
        # example:
        #1	54
        ($start, $stop) = ($1, $2);

        # does line order make sense?
        if($prv_was_coords_feature || # previous line was a coords line with a feature (common)
           $prv_was_coords_only) {    # previous line was a coords line without a feature (common)
          # line order makes sense, keep going...

          $fac .= $fac_sep . $start . $fac_sep . $stop;
          
          # update '$prv_*' values that we use to make sure line order makes sense
          $prv_was_accn           = 0;
          $prv_was_coords_feature = 0;
          $prv_was_coords_only    = 1;
          $prv_was_quals          = 0;
          #printf("set prv_was_coords_only\n");
        }
        else { # line order is unexpected
          DNAORG_FAIL("ERROR in $sub_name, unexpected line order (coords_only line) at line $line_ctr: $line", 1, $FH_HR);
        }
      }
      # -------------------------------------------------------
      elsif(($line =~ /^\s+\S+\s+.+$/) || 
            ($line =~ /^\s+\S+$/)) { 
        # QUALIFIER LINE
        # examples:
        #			gene	AR1
        #			locus_tag	PhyvvsAgp1

        # before parsing it, do two sanity checks
        if(! defined $fac)     { DNAORG_FAIL("ERROR in $sub_name, coordinates undefined at a qualifier line", 1, $FH_HR); }
        if(! defined $feature) { DNAORG_FAIL("ERROR didn't read feature line before line: $line", 1, $FH_HR); }
        # and determine if we even care about this feature
        if((! $do_only_feature) || ($cap_feature eq $cap_only_feature)) { 
          # does line order make sense?
          if($prv_was_coords_feature || 
             $prv_was_coords_only    ||
             $prv_was_quals) { 
            # line order makes sense, keep going...
            if(! $prv_was_quals) { 
              # first quals line for this feature
              # at this point, we know that we have the full coordinates for the feature
              # so initialize the information
              if(! exists $fac_HHAR->{$feature}) {
                %{$fac_HHAR->{$feature}} = (); 
              }
              if(! exists $fac_HHAR->{$feature}{$faccn}) { 
                @{$fac_HHAR->{$feature}{$faccn}} = ();
              }
              push(@{$fac_HHAR->{$feature}{$faccn}}, $fac);
              if(! exists $quals_HHAR->{$feature}) { 
                %{$quals_HHAR->{$feature}} = ();
                @{$quals_HHAR->{$feature}{$fac}} = ();
              }
            } # end of 'if(! $prv_was_quals)'

            # now parse the line;
            # examples:
            #			gene	AR1
            #			locus_tag	PhyvvsAgp1
            if($line =~ /^\s+(\S+)\s+(.+)$/) { 
              ($qname, $qval) = ($1, $2);
            }
            elsif($line =~ /^\s+(\S+)$/) { 
              $qname = $1;
              $qval = "<no_value>";
            }
            else { 
              DNAORG_FAIL("ERROR in $sub_name, didn't parse quals line on second pass: $line", 1, $FH_HR); 
            }
            if($qname =~ m/\Q$qnqv_sep/)   { DNAORG_FAIL("ERROR in $sub_name, qualifier_name $qname has the string $qnqv_sep in it", 1, $FH_HR); }
            if($qval  =~ m/\Q$qnqv_sep/)   { DNAORG_FAIL("ERROR in $sub_name, qualifier_value $qval has the string $qnqv_sep in it", 1, $FH_HR); }
            my $qnqv = $qname . $qnqv_sep . $qval; # this is how we store the qualifier name and value, as a concatenated string in quals_HHAR
            push(@{$quals_HHAR->{$feature}{$fac}}, $qnqv);

            # and update the column data structures which just keep info on names and order of columns
            if(! exists $column_HH{$feature}) { 
              %{$column_HH{$feature}} = ();
              @{$column_HAR->{$feature}} = (); 
            }
            if(! exists $column_HH{$feature}{$qname}) { 
              push(@{$column_HAR->{$feature}}, $qname);
              $column_HH{$feature}{$qname} = 1;
            }
          }
          else { # unexpected line order
            DNAORG_FAIL("ERROR in $sub_name, unexpected line order (quals line) at line $line_ctr: $line", 1, $FH_HR);
          }          
        } # end of 'if((! $do_only_feature) || ($cap_feature ne $cap_only_feature))'
        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 0;
        $prv_was_coords_feature = 0;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 1;
        #printf("set prv_was_quals\n");
      }
      # -------------------------------------------------------
      else { 
        DNAORG_FAIL("ERROR in $sub_name, unable to parse line $line_ctr: $line", 1, $FH_HR); 
      }
      # -------------------------------------------------------
    }
  }
}

#################################################################
# Subroutine: parseEdirectMatPeptideFile()
# Incept:     EPN, Thu Feb 11 14:01:30 2016
#
# Purpose:   Given a mature peptide file output from a
#            'efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product'
#            command, parse that file into usable data structures.
#
#            Can be called by the wrapper subroutine: 
#             edirectFtableOrMatPept2SingleFeatureTableInfo().
#
#            Caller will commonly call outputSingleFeatureTable() 
#            or getSingleFeatureTableInfo() after calling this
#            subroutine.
#
#            Similar to (and shares some code with) parseEdirectFtableFile().
#
# Args:
#   $infile:          name of ftable file
#   $dummy_column:    name for dummy columns that we won't output, can be undef
#   $sep_HR:          ref to hash of 'separation' values, keys: "qnqv", "fac", and "qval", 
#                     used to split values that have been concatenated with these 'separation'
#                     values in between, PRE-FILLED
#   $quals_HHAR:      ref to 2D hash of arrays with values that could eventually go 
#                     into a single-feature table, FILLED HERE
#                     key 1: feature name, e.g. "CDS"
#                     key 2: 'fac' string: <full_accession><fac_seq><coordinates>
#                     value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>    
#   $faccn_AR:        REF to array of full accessions read, FILLED HERE
#   $fac_HHAR:        REF to 2D hash of arrays that is used to easily determine
#                     list of 2D keys ('fac's) in quals_HHA for a given feature and faccn
#                     FILLED HERE
#                     key 1: feature name: e.g. "CDS"
#                     key 2: 'faccn', full accession
#                     value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
#   $faccn2accn_HR:   REF to hash, key: full accession, value: short accession, 
#                     used for convenience in output function, FILLED HERE
#   $column_HAR:      REF to hash of arrays, key: feature, e.g. "CDS"
#                     value: array of qualifiers that exist for that feature in at least 1 faccn
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
#
# Returns:  void
#
# Dies:     if <$infile> file is in unexpected format
#################################################################
sub parseEdirectMatPeptideFile {
  my $sub_name = "parseEdirectMatPeptideFile()";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $dummy_column, $sep_HR, $quals_HHAR, $faccn_AR, $fac_HHAR, $faccn2accn_HR, $column_HAR, $FH_HR) = @_;

  open(IN, $infile) || fileOpenFailure($infile, $FH_HR, $!, "reading");

  # example lines of a .mat_peptide file
  # NC_009942.1	97..465	anchored capsid protein C	
  # NC_009942.1	97..411	capsid protein C	
  # NC_009942.1	466..966	membrane glycoprotein precursor M	
  # NC_009942.1	join(2470..3552,3552..3680)	nonstructural protein NS1'

  my $faccn     = undef;   # full accession, e.g. ref|NC_000883.2|
  my $fac       = undef;   # full accession and coords, accession and each segment coords separated by $fac_sep
  my $coords    = undef;   # coordinates
  my $faccn2    = undef;   # another full accession
  my %column_HH = ();      # existence 2D hash, used to aide construction of %column_HA
                           # key 1: feature, e.g. "CDS"
                           # key 2: qualifier
                           # value: '1' (always)

  my $fac_sep  = $sep_HR->{"fac"};
  my $qnqv_sep = $sep_HR->{"qnqv"};
  my $qval_sep = $sep_HR->{"qval"};

  my $feature = "mat_peptide";
  while(my $line = <IN>) { 
    if($line =~ m/\w/) { 
      if($line =~ /(\S+)\s+(\S+)\s*(.*)$/) { 
        my ($acc, $coords, $product) = ($1, $2, $3);
        $product =~ s/\s+$//; # remove trailing whitespace
        $faccn = $acc;
        if(! exists $faccn2accn_HR->{$acc}) { 
          push(@{$faccn_AR}, $acc);
          $faccn2accn_HR->{$acc} = $acc;
        }

        # now we have to break down coords just so we can store each segment separately in fac_HHA, 
        # because the feature table has each segment separated, and we want to process and OUTPUT
        # both feature tables and mat_peptide information the same way after we parse it here. At
        # that point we'll put the coordinates together again. 
        my @starts_A = ();
        my @stops_A  = ();
        my $nsegments = 0;
        startsStopsFromCoords($coords, \@starts_A, \@stops_A, \$nsegments, $FH_HR);
        
        # printf("coords:    $coords\n");
        # printf("nsegments: $nsegments\n");
        for(my $i = 0; $i < $nsegments; $i++) { 
          my $start = $starts_A[$i];
          my $stop  = $stops_A[$i];
          # printf("starts_A[$i]: $starts_A[$i]\n");
          # printf("stops_A[$i]:  $stops_A[$i]\n");
          if($i == 0) { 
            $fac = $faccn . $fac_sep . $start . $fac_sep . $stop;
          }
          else { 
            $fac .= $fac_sep . $start . $fac_sep . $stop;
          }
        }

        if(! exists $fac_HHAR->{"mat_peptide"}) {
          %{$fac_HHAR->{"mat_peptide"}} = (); 
        }
        push(@{$fac_HHAR->{"mat_peptide"}{$acc}}, $fac);
        
        # first add the 'dummy' qual
        my $qname = $dummy_column;
        my $qval  = "<no_value>";
        if($qname =~ m/\Q$qnqv_sep/)   { DNAORG_FAIL("ERROR in $sub_name, qualifier_name $qname has the string $qnqv_sep in it", 1, $FH_HR); }
        if($qval  =~ m/\Q$qnqv_sep/)   { DNAORG_FAIL("ERROR qualifier_value $qval has the string $qnqv_sep in it", 1, $FH_HR); }
        my $qnqv = $qname . $qnqv_sep . $qval; # this is how we store the qualifier name and value, as a concatenated string in quals_HHAR
        push(@{$quals_HHAR->{"mat_peptide"}{$fac}}, $qnqv);
        
        # and update the column data structures which just keep info on names and order of columns
        if(! exists $column_HH{"mat_peptide"}) { 
          %{$column_HH{"mat_peptide"}} = ();
          @{$column_HAR->{"mat_peptide"}} = (); 
        }
        if(! exists $column_HH{"mat_peptide"}{$qname}) { 
          push(@{$column_HAR->{"mat_peptide"}}, $qname);
          $column_HH{"mat_peptide"}{$qname} = 1;
        }
        
        if(defined $product && $product ne "") { 
          # now if the product qualifier has a value, add that too
          $qname = "product";
          $qval  = $product;
          $qnqv = $qname . $qnqv_sep . $qval;
          push(@{$quals_HHAR->{"mat_peptide"}{$fac}}, $qnqv);
          
          if(! exists $column_HH{"mat_peptide"}{$qname}) { 
            push(@{$column_HAR->{"mat_peptide"}}, $qname);
            $column_HH{"mat_peptide"}{$qname} = 1;
          }
        }
      }
    } # end of 'if($line =~ m/\w/)
  } # end of '$line = <IN>' 
}


#################################################################
# Subroutine: helperBreakdownFac()
# Incept      EPN, Thu Feb 11 14:14:23 2016
#
# Purpose:    Breakdown a 'fac' string into it's parts and 
#             create a string in NCBI coordinate format from it.
#             A 'fac' string has accessions and coordinates in it.
#           
#             A 'helper' function called by 'outputSingleFeatureTable()'
#             and 'getSingleFeatureTableInfo()'.
#
# Arguments:
#   $fac:     full accession and coordinates, concatenated together
#             can be multiple coordinates e.g. "NC_000001:100:300:350:400"
#             means two segments: 100..300 and 350..400.
#   $fac_sep: character in between each token in $fac (":" in above example)
#   $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:
#   Four values:
#   $faccn:       full accession
#   $ncbi_coords: NCBI format for coordinates
#   $sort_coord:  minimum coordinate of all segments, possibly useful for ordering multiple features
#   $strand:      '+' if all segments are on fwd strand
#                 '-' if all segments are on neg strand
#                 '?' if all segments are 1 nucleotide (strand is consequently uncertain)
#                 '!' if >= 1 segment on two or more of following: fwd strand, rev strand, uncertain
#
# Dies:       if there's not an odd number of tokens
#################################################################
sub helperBreakdownFac {
  my $sub_name = "helperBreakdownFac()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($fac, $fac_sep, $FH_HR) = @_;

  my @elA = split($fac_sep, $fac);
  my $nel = scalar(@elA);
  if($nel % 2 == 0) { DNAORG_FAIL("ERROR in $sub_name, unexpectedly there's an even number of tokens: fac: $fac", 1, $FH_HR); }

  my $faccn = $elA[0];
  my ($istart, $istop) = (1, 2);
  my $start  = undef;   # start position
  my $stop   = undef;   # stop position
  my $is_rev = 0; # true if current segment is on reverse strand
  my $nfwd   = 0; # number of segments that are on forward strand
  my $nrev   = 0; # number of segments that are on reverse strand
  my $nunc   = 0; # number of segments for which strand is uncertain (start == stop)
  my $some_revcomp  = undef;
  my @ncbi_coords_A = (); # [0..$nsegments-1] strings to concatenate to make ncbi coords
  my @is_rev_A      = (); # [0..$nsegments-1] '1' if current segment is on reverse strand
  my $nsegments     = 0;
  my $have_carrot   = 0; # set to '1' if $start has a '^' at the end, e.g. 804^ 805
  my $sort_coord = undef; # start position of first segment
  my $min_coord  = undef;
  while($istop <= $nel) { 
    my ($orig_start, $orig_stop) = ($elA[$istart], $elA[$istop]);
    $start = $orig_start;
    $stop  = $orig_stop;
    if($start =~ s/\^$//) { 
      $have_carrot = 1;
    }
    $start =~ s/^<//;
    $start =~ s/^>//;
    $stop  =~ s/^<//;
    $stop  =~ s/^>//;
    $min_coord = ($start < $stop) ? $start : $stop;
    if((! defined $sort_coord) || ($min_coord < $sort_coord)) { $sort_coord = $min_coord; }

    if($have_carrot) { 
      if(abs($start-$stop) > 1) { DNAORG_FAIL("ERROR found carrot but range is not exactly 1 nt: fac: $fac", 1, $FH_HR); }
    }

    if($start == $stop) { # special case, we can't tell if we're in reverse complement or not
      $nunc++;
      push(@is_rev_A, 0);
    }
    elsif($start > $stop) {
      $nrev++;
      $is_rev = 1;
      push(@is_rev_A, 1);
    }
    else { # $start < $stop, not reverse complement
      $nfwd++;
      push(@is_rev_A, 0);
    }
    if($have_carrot) { 
      push(@ncbi_coords_A, $orig_start . $orig_stop);
    }
    else { 
      if($is_rev) { push(@ncbi_coords_A, $orig_stop .  ".." . $orig_start); }
      else        { push(@ncbi_coords_A, $orig_start . ".." . $orig_stop); }
    }
    $nsegments++;
    $istart += 2; 
    $istop  += 2;
  }
  if($have_carrot) { 
    if($nsegments > 1) { DNAORG_FAIL("ERROR found carrot but more than one segment: fac: $fac", 1, $FH_HR); }
  }

  my $ncbi_coords = "";
  if($nfwd > 0 && $nrev > 0) { 
    # special case, we need to put 'complement(' ')' around each rev segment separately
    # e.g: join(161990..162784,complement(88222..88806),complement(86666..87448))
    for(my $i = 0; $i < $nsegments; $i++) { 
      if($i > 0) { $ncbi_coords .= ","; }
      if($is_rev_A[$i]) { 
        # first, deal with a non-obvious situation, where in the feature table
        # '>' and '<' characters indicating incompleteness are inverted relative
        # to how they are in the actual annotation. 
        # NC_007030.2 complement(4370..>4576)
        # is in the feature table as: <4576	4370	CDS
        $ncbi_coords_A[$i] =~ s/\</\!/;
        $ncbi_coords_A[$i] =~ s/\>/\</;
        $ncbi_coords_A[$i] =~ s/\!/\>/;
        $ncbi_coords .= "complement(" . $ncbi_coords_A[$i] . ")" 
      }
      else { 
        $ncbi_coords .= $ncbi_coords_A[$i]; 
      }
    }
  }
  else { # normal case, all exons/segments are on the same strand
    # if we're on the reverse strand, we need to reverse the order of the 
    # exons, because the order of reverse strand exons in a feature table is 
    # opposite what it is in Entrez, and our other scripts use Entrez
    # format, so we enforce that convention here.
    if($nrev > 0) { 
      for(my $i = $nsegments-1; $i >= 0; $i--) { 
        if($i < ($nsegments-1)) { $ncbi_coords .= ","; }
        $ncbi_coords .= $ncbi_coords_A[$i];
      }
    }
    else { # positive strand
      for(my $i = 0; $i < $nsegments; $i++) { 
        if($i > 0) { $ncbi_coords .= ","; }
        $ncbi_coords .= $ncbi_coords_A[$i];
      }
    }
  }
  if($nsegments > 1) { # more than one segment
    $ncbi_coords = "join(" . $ncbi_coords . ")";
  }
  # now add complement for cases where are exons/segments are on reverse strand
  # impt to do this after the join, so we get complement(join()) instead of
  # join(complement())
  if($nfwd == 0 && $nrev > 0) { # all segments are on reverse strand
    # first, deal with a non-obvious situation, where in the feature table
    # '>' and '<' characters indicating incompleteness are inverted relative
    # to how they are in the actual annotation. 
    # NC_007030.2 complement(4370..>4576)
    # is in the feature table as: <4576	4370	CDS
    $ncbi_coords =~ s/\</\!/g;
    $ncbi_coords =~ s/\>/\</g;
    $ncbi_coords =~ s/\!/\>/g;
    # now add the 'complement()'
    $ncbi_coords = "complement(" . $ncbi_coords . ")";
  }

  # printf("in $sub_name input: $fac, returning $faccn $coords ncbi_coords:$ncbi_coords\n");

  # determine strand
  my $ret_strand = undef;
  if   ($nfwd == $nsegments)     { $ret_strand = "+"; }
  elsif($nrev == $nsegments)     { $ret_strand = "-"; }
  elsif($nunc  > 0)              { $ret_strand = "?"; }
  elsif($nfwd  > 0 && $nrev > 0) { $ret_strand = "!"; }
  else                           { DNAORG_FAIL("ERROR in $sub_name, unable to determine strand for fac: $fac", 1, $FH_HR); }

  return($faccn, $ncbi_coords, $sort_coord, $ret_strand);
}

#################################################################
# Subroutine: startsStopsFromCoords()
# Incept:     EPN, Thu Feb 11 14:22:54 2016
#
# Purpose:    Extract the starts and stops from a coords string.
#
# Args:
#   $coords:    the coords string
#   $starts_AR: REF to array to fill with start positions, FILLED HERE
#   $stops_AR:  REF to array to fill with stop positions, FILLED HERE
#   $nexons_R:  REF to scalar that fill with the number of exons, FILLED HERE
#   $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:      void; but fills @{$starts_AR}, @{$stops_AR}, and $$nexons_R.
#
# Dies:         if we can't parse $coords
#################################################################
sub startsStopsFromCoords { 
  my $sub_name = "startsStopsFromCoords()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($coords, $starts_AR, $stops_AR, $nexons_R, $FH_HR) = @_;

  @{$starts_AR} = ();
  @{$stops_AR}  = ();
  $$nexons_R    = 0;
  
  my $orig_coords = $coords;
  # Examples:
  # complement(2173412..2176090)
  # complement(join(226623..226774, 226854..229725))

  # remove 'complement('  ')'
  $coords =~ s/^complement\(//;
  $coords =~ s/\)$//;

  # remove 'join('  ')'
  $coords =~ s/^join\(//;
  $coords =~ s/\)$//;

  my @el_A = split(/\s*\,\s*/, $coords);

  my $length = 0;
  foreach my $el (@el_A) { 
    # rare case: remove 'complement(' ')' that still exists:
    $el =~ s/^complement\(//;
    $el =~ s/\)$//;
    $el =~ s/\<//; # remove '<'
    $el =~ s/\>//; # remove '>'
    if($el =~ m/^(\d+)\.\.(\d+)$/) { 
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $2);
      $$nexons_R++;
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse $orig_coords", 1, $FH_HR); 
    }
  }

  return;
}

#################################################################
# Subroutine: stripVersion()
# Incept:     EPN, Thu Feb 11 14:25:52 2016
#
# Purpose:    Given a ref to an accession.version string, remove the version.
#
# Arguments: 
#   $accver_R: ref to accession version string to remove version from
#
# Returns:    Nothing, $$accver_R has version removed
#
# Dies:       never
#################################################################
sub stripVersion {
  my $sub_name  = "stripVersion()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($accver_R) = (@_);

  $$accver_R =~ s/\.[0-9]*$//; # strip version

  return;
}


#################################################################
# Subroutine: fetchReferenceFeatureSequences()
# Incept:     EPN, Thu Feb 11 14:35:31 2016
#
# Synopsis: Fetch the reference sequences that models will be built
#           from and populate data structures with information on
#           features and models.
#
# Arguments:
#  $sqfile:            Bio::Easel::SqFile object, the sequence file we'll fetch from
#  $ref_seq_accn:      sequence accession of reference
#  $ref_totlen:        length of reference 
#  $do_circular:       true if genome is circularized and we allow features to span stop..start boundary
#  $mdl_info_HAR:      ref to hash of arrays with information on the models, FILLED HERE
#  $ftr_info_HAR:      ref to hash of arrays with information on the features, ADDED TO HERE
#  $all_stk_file:      name of output file we will write the stockholm single sequence 
#                      alignments of all features to
#  $FH_HR:             REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void; fills @{$values_HHAR}
#
# Dies:       if $ftr_info_HAR is not valid upon entering
#################################################################
sub fetchReferenceFeatureSequences {
  my $sub_name = "fetchReferenceFeatureSequences()";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($sqfile, $ref_seq_accn, $ref_totlen, $out_root, $do_circular, $mdl_info_HAR, $ftr_info_HAR, $all_stk_file, $FH_HR) = @_;

  # contract check
  # ftr_info_HAR should have array data for keys "ref_coords", "ref_strand"
  my @reqd_ftr_info = ("ref_coords", "ref_strand", "has_models");
  my $nftr             = validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, \@reqd_ftr_info, $FH_HR);
  my $nftr_with_models = getNumFeaturesWithModels($ftr_info_HAR, $FH_HR);

  my $cur_out_root;
  my $cur_name_root;
  my $dir_tail = $out_root;
  $dir_tail =~ s/^.+\///;
  my $fetch_input;
  my $fetch_output;
  my $nmdl = 0;               # number of HMMs (and alignments used to build those HMMs)
  my @mdl_A = ();             # [0..$nmdl-1]: array of model names, also name of stockholm alignments used to build those models
  my %mdllen_H          = (); # key: model name from @mdl_A, value is model length
  my @ref_nexons_A      = (); # [0..$c..$ref_nmft-1]: number of exons in CDS or mat_peptide $c+1
  my $ref_tot_nexons    = 0;  # total number of exons in all CDS or mat_peptides
  my @indi_ref_name_A   = (); # [0..$nmdl-1]: name of individual stockholm alignments and models
  my @indi_cksum_stk_A  = (); # [0..$nmdl-1]: checksum's of each named individual stockholm alignment

  for(my $i = 0; $i < $nftr; $i++) { 
    my $do_model = $ftr_info_HAR->{"has_models"}[$i];

    if($do_model) { 
      my $cds_or_mp = ($ftr_info_HAR->{"type"}[$i] eq "mp") ? "mp" : "cds";
      my $i2print   = $ftr_info_HAR->{"type_idx"}[$i];

      # determine start and stop positions of all exons/segments
      my @starts_A = ();
      my @stops_A  = ();
      my $nexons   = 0;
      startsStopsStrandsFromCoordsLength($ftr_info_HAR->{"ref_coords"}[$i], $ref_totlen, $do_circular, \@starts_A, \@stops_A, undef, \$nexons, $FH_HR);
      $ftr_info_HAR->{"nexons"}[$i] = $nexons;
      my $strand = $ftr_info_HAR->{"ref_strand"}[$i];

      # if we're on the negative strand, reverse the arrays, they'll be in the incorrect order
      if($strand eq "-") { 
        @starts_A = reverse @starts_A;
        @stops_A  = reverse @stops_A;
      }

      for(my $e = 0; $e < $nexons; $e++) { 
        if($nexons > 1) { 
          $cur_out_root  = $out_root . ".ref." . $cds_or_mp . "." . $i2print . ".exon." . ($e+1);
          $cur_name_root = $dir_tail . ".ref." . $cds_or_mp . "." . $i2print . ".exon." . ($e+1);
        }
        else { 
          $cur_out_root  = $out_root . ".ref." . $cds_or_mp . "." . $i2print;
          $cur_name_root = $dir_tail . ".ref." . $cds_or_mp . "." . $i2print;
        }
        printf("out_root: $out_root\n");
        printf("dir_tail: $dir_tail\n");
        printf("cds_or_mp: $cds_or_mp\n");
        printf("i2print: $i2print\n");
    
        # determine start and stop of the region we are going to fetch
        my $start = $starts_A[$e];
        my $stop  = $stops_A[$e];
        if($strand eq "-") { # swap start and stop
          my $tmp = $start;
          $start = $stop;
          $stop  = $tmp;
        }
        my @fetch_AA = ();
    
        push(@fetch_AA, [$cur_name_root, $start, $stop, $ref_seq_accn]);
    
        # fetch the sequence
        my $cur_fafile = $cur_out_root . ".fa";
        $sqfile->fetch_subseqs(\@fetch_AA, undef, $cur_fafile);

        # reformat to stockholm
        my $cur_stkfile = $cur_out_root . ".stk";
        my $cmd = "esl-reformat --informat afa stockholm $cur_fafile > $cur_stkfile";
        runCommand($cmd, 0, $FH_HR);
    
        # annotate the stockholm file with a blank SS and with a name
        my $cur_named_stkfile = $cur_out_root . ".named.stk";
        my ($mdllen, $cksum) = addNameAndBlankSsToStockholmAlignment($cur_name_root, 1, $cur_stkfile, $cur_named_stkfile, $FH_HR); # 1: add blank SS_cons line
        $mdl_info_HAR->{"cmname"}[$nmdl]   = $cur_name_root;
        $mdl_info_HAR->{"checksum"}[$nmdl] = $cksum;
        $mdl_info_HAR->{"length"}[$nmdl]   = $mdllen;

        # store information on this model's name for output purposes
        $mdl_info_HAR->{"out_tiny"}[$nmdl] = sprintf("%s#%s", 
                                                     (($cds_or_mp eq "mp") ? "MP" : "CDS"), 
                                                     ($nexons == 1) ? sprintf("%d", $i2print) : sprintf("%d.%d", $i2print, ($e+1)));
        if($e == ($nexons-1)) { 
          my $short = ($cds_or_mp eq "mp") ? sprintf("MP #%d", $i2print) : sprintf("CDS #%d", $i2print);
          my $tiny  = $short;
          $tiny =~ s/\s+//g; # remove whitespace
          if($nexons > 1) { $short .= sprintf(" [$nexons %s; $strand]", ($cds_or_mp eq "mp") ? "segments" : "exons"); }
          else            { $short .= sprintf(" [single %s; $strand]",  ($cds_or_mp eq "mp") ? "segment"  : "exon"); }
          $ftr_info_HAR->{"out_tiny"}[$i]    = $tiny;
          $ftr_info_HAR->{"out_short"}[$i]   = $short;
        }
        
        # now append the named alignment to the growing stockholm alignment database $all_stk_file
        $cmd = "cat $cur_named_stkfile";
        if($nmdl == 0) { $cmd .= " >  $all_stk_file"; }
        else           { $cmd .= " >> $all_stk_file"; }
        runCommand($cmd, 0, $FH_HR);

        $mdl_info_HAR->{"map_mft"}[$nmdl]  = $i;
        $mdl_info_HAR->{"is_first"}[$nmdl] = ($e == 0)           ? 1 : 0;
        $mdl_info_HAR->{"is_final"}[$nmdl] = ($e == ($nexons-1)) ? 1 : 0;

        $ftr_info_HAR->{"map_exon"}[$i]  = $e;
        $ftr_info_HAR->{"map_nexon"}[$i] = $nexons;
        if($e == 0) { 
          $ftr_info_HAR->{"first_mdl"}[$i] = $nmdl;
        }
        if($e == ($nexons-1)) { 
          $ftr_info_HAR->{"final_mdl"}[$i] = $nmdl;
        }
        $nmdl++;
      }
    } # end of 'if($do_model)'
    else { 
      # this feature is NOT modelled ($ftr_info_HAR->{"has_models"}[$i] == 0)
      $ftr_info_HAR->{"nexons"}[$i]    = 0;
      $ftr_info_HAR->{"out_tiny"}[$i]  = "?";
      $ftr_info_HAR->{"out_short"}[$i] = "?";
      $ftr_info_HAR->{"map_exon"}[$i]  = -1;
      $ftr_info_HAR->{"map_nexon"}[$i] = -1;
      $ftr_info_HAR->{"first_mdl"}[$i] = -1;
      $ftr_info_HAR->{"final_mdl"}[$i] = -1;
    }
  }

  return;
}

#################################################################
# Subroutine: mapFeaturesToModels()
# Incept:     EPN, Thu Feb 11 14:48:20 2016
#
# Purpose:    Fill the two-dimensional array that maps features 
#             to model indices.
#
# Arguments:
#   $ftr2mdl_map_AAR: ref to 2D array that maps features to models
#                     1D:    feature index
#                     2D:    number of exons/segments for this feature ($ftr_info_HAR->{"nexon"})
#                     value: model index in $mdl_info_HAR
#   $ftr_info_HAR:    ref to hash of arrays with info on features, PRE-FILLED
#   $mdl_info_HAR:    ref to hash of arrays with info on model, PRE-FILLED
#   $do_print:        '1' to print out map after creating map
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $ftr_info_HAR or $mdl_info_HAR is invalid upon entering
#################################################################
sub mapFeaturesToModels { 
  my $sub_name  = "mapFeaturesToModels()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr2mdl_map_AAR, $ftr_info_HAR, $mdl_info_HAR, $do_print, $FH_HR) = (@_);

  # contract check, make sure that %{$mdl_info_HAR} and %{$ftr_info_HAR} are valid and have the info we need
  my @reqd_mdl_info = ("is_first", "is_final");
  my $nmdl = validateAndGetSizeOfInfoHashOfArrays($mdl_info_HAR, \@reqd_mdl_info);

  my @reqd_ftr_info = ("nexon");
  my $nftr = validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, \@reqd_ftr_info);

  my $ftr_idx = -1;
  my $mdl_idx = 0;
  for(my $h = 0; $h < $nmdl; $h++) { 
    if($mdl_info_HAR->{"is_first"}[$h]) { 
      $ftr_idx++;
      $mdl_idx = 0;
    }
    $ftr2mdl_map_AAR->[$ftr_idx][$mdl_idx] = $h;
    if($do_print) { 
      printf("ftr2mdl_map_AAR->[$ftr_idx][$mdl_idx]: $h\n");
    }
    $mdl_idx++;
    
    # and verify that this agrees with %{$ftr_info_HAR}
    if($mdl_info_HAR->{"is_final"}[$h]) { 
      if($ftr_info_HAR->{"nexon"}[$ftr_idx] != $mdl_idx) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, observed $mdl_idx total models for feature $ftr_idx, but expected %d", $ftr_info_HAR->{"nexon"}{$ftr_idx}), 1, $FH_HR);
      }
    }
  }
  return;
}

#################################################################
# Subroutine: determineFeatureTypes()
# Incept:     EPN, Thu Feb 11 14:50:53 2016
#
# Purpose:    Determine the type of each feature.
#
# Arguments:
#   $nmp:              number of mature peptides, may be 0
#   $cds2pmatpept_AAR: 1st dim: cds index (-1, off-by-one), 
#                      2nd dim: value array of primary matpept indices that comprise this CDS, 
#                      OR undefined if all features are CDS and there are no mature peptides; 
#                      PRE-FILLED
#   $ftr_info_HAR:     ref to hash of arrays to fill "type" and "type_idx" arrays for
#                      "type"     values are one of: "cds-notmp", "cds-mp", or "mp",
#                      "type_idx" values are indices of each feature for its type
#                      so a '2' means it's the 2nd of its type.
#   $FH_HR:            REF to hash of file handles, including "log" and "cmd"
#             
# Returns:    void
#
# Dies:       if $ftr_info_HAR is invalid upon entering
#################################################################
sub determineFeatureTypes {
  my $sub_name  = "determineFeatureTypes()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($nmp, $cds2pmatpept_AAR, $ftr_info_HAR, $FH_HR) = (@_);

  my $c; # counter

  my $nftr = validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, undef, $FH_HR);

  # initialize arrays
  @{$ftr_info_HAR->{"type"}}       = ();
  @{$ftr_info_HAR->{"type_idx"}}   = ();
  @{$ftr_info_HAR->{"has_models"}} = ();

  if(! defined $cds2pmatpept_AAR) { 
    # trivial, all features are CDS
    for($c = 0; $c < $nftr; $c++) { 
      $ftr_info_HAR->{"type"}[$c]       = "cds-notmp";
      $ftr_info_HAR->{"has_models"}[$c] = 1; # this feature does have homology models
      $ftr_info_HAR->{"type_idx"}[$c]   = $c+1;
    }
  }
  else { 
    # some features are mature peptides, some are CDS, some of the CDS are comprised of mature peptides, 
    # and (possibly) some of the CDS are not

    # first comes the $nmp mature peptides
    my $ncds = $nftr - $nmp;
    for($c = 0; $c < $nmp; $c++) { 
      $ftr_info_HAR->{"type"}[$c]       = "mp";
      $ftr_info_HAR->{"has_models"}[$c] = 1; # this feature does have homology models
      $ftr_info_HAR->{"type_idx"}[$c]   = $c+1;
    }
    # then come the CDS, but we need to consult $cds2pmatpept_AAR to see if they are comprised of mature peptides ("cds-mp") 
    # or not ("cds-notmp")
    my $nnotmp_idxed = 0; # we increment this as we index cds-notmp 
    my $nnotmp_extra = 0; # # of cds-notmp that occur after final cds-mp (idx >= scalar(@cds2pmatpept_AA))
    for(my $cds_idx = 0; $cds_idx < $ncds; $cds_idx++) { 
      $c = $nmp + $cds_idx;
      if((defined $cds2pmatpept_AAR->[$cds_idx]) && (exists $cds2pmatpept_AAR->[$cds_idx])) { 
        # a CDS that is comprised of mature peptides
        $ftr_info_HAR->{"type"}[$c]         = "cds-mp";
        $ftr_info_HAR->{"type_idx"}->[$c]   = $cds_idx+1;
        $ftr_info_HAR->{"has_models"}->[$c] = 0; # this feature does not have homology models
      }
      else { 
        # a CDS that is NOT comprised of mature peptides
        $ftr_info_HAR->{"type"}[$c] = "cds-notmp";
        $ftr_info_HAR->{"type_idx"}->[$c]   = $cds_idx+1;
        $ftr_info_HAR->{"has_models"}->[$c] = 1; # this feature does have homology models
      }
    }
  }
  return;
}

#################################################################
# Subroutine: getNumFeaturesWithModels()
# Incept:     EPN, Thu Feb 11 14:50:53 2016
#
# Purpose:    Return number of features for which "has_models" array values are 1.
#
# Arguments: 
#   $ftr_info_HAR: ref to hash of arrays to fill "type" and "type_idx" arrays for
#                  must have at least the "has_models" array, else we die 
#   $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    number of features for which the has_models array values are 1
#
# Dies:       if $ftr_info_HAR is not valid up on entering.
#################################################################
sub getNumFeaturesWithModels { 
  my $sub_name  = "getNumFeaturesWithModels()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_HAR, $FH_HR) = (@_);

  my @reqd_ftr_info = ("has_models");
  my $nftr = validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, \@reqd_ftr_info, $FH_HR);
  
  my $nftr_with_models = 0;
  for(my $i = 0; $i < $nftr; $i++) { 
    if($ftr_info_HAR->{"has_models"}[$i]) { $nftr_with_models++; }
  }
  return $nftr_with_models;
}

#################################################################
# Subroutine: startsStopsStrandsFromCoordsLength()
# Incept:     EPN, Thu Feb 11 14:57:22 2016
#
# Purpose:    Determine the starts, stops and strands from a coords string
#             and length.
# 
# Arguments:
#   $coords:      the coords string
#   $totlen:      total length of sequence
#   $do_circular: '1' if we're searching a circular (duplicated) genome and we'll allow stop..start boundary spans
#   $starts_AR:   ref to array to fill with start positions
#   $stops_AR:    ref to array to fill with stop positions
#   $strands_AR:  ref to array to fill with strands of each exon, can be undef
#   $nexons_R:    ref to scalar that fill with the number of exons
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void; but fills @{$starts_AR}, @{$stops_AR}, and $$nexons_R.
#
# Dies:       if we see a feature that spans stop..start but $do_circular is 0
#################################################################
sub startsStopsStrandsFromCoordsLength { 
  my $sub_name = "startsStopsStrandsFromCoordsLength()";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($coords, $totlen, $do_circular, $starts_AR, $stops_AR, $strands_AR, $nexons_R, $FH_HR) = @_;

  # zero/initialize what we will be determining in this subroutine
  @{$starts_AR} = ();
  @{$stops_AR}  = ();
  $$nexons_R    = 0;
  
  my $orig_coords = $coords;
  # Examples:
  # complement(2173412..2176090)
  # complement(join(226623..226774, 226854..229725))

  # remove 'complement('  ')'
  my $strand = "+";
  if($coords =~ m/^complement\(/) { 
    $coords =~ s/^complement\(//;
    $strand = "-";
  }
  $coords =~ s/\)$//;

  # remove 'join('  ')'
  $coords =~ s/^join\(//;
  $coords =~ s/\)$//;

  my @el_A = split(/\s*\,\s*/, $coords);

  my $length = 0;
  my $cur_strand = $strand;
  foreach my $el (@el_A) { 
    # rare case: remove 'complement(' ')' that still exists:
    $cur_strand = $strand;
    if($el =~ m/^complement\(/) { 
      DNAORG_FAIL("ERROR in $sub_name: found internal complement in coords string $coords, we assume all exons are on same strand...", 1, $FH_HR); 
      $el =~ s/^complement\(//;
      if($cur_strand eq "-") { DNAORG_FAIL("ERROR in $sub_name, found nested 'complement' annotations in coord string: $coords", 1, $FH_HR); }
      $cur_strand = "-";
    }
    $el =~ s/\)$//;
    $el =~ s/\<//; # remove '<'
    $el =~ s/\>//; # remove '>'
    if($el =~ m/^(\d+)\.\.(\d+)$/) { 
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $2);
      if(defined $strands_AR) { push(@{$strands_AR}, $cur_strand); }
      $$nexons_R++;
    }
    else { 
      DNAORG_FAIL("ERROR unable to parse $orig_coords in $sub_name", 1, $FH_HR); 
    }
  }

  # check if we have a spanning exon (that spans stop..start) and if we do
  # and (! $do_circular) then die, because that shouldn't happen.
  my $have_spanning_exon = checkForSpanningSequenceSegments($starts_AR, $stops_AR, $nexons_R, 0, $strand, $totlen); # 0 says: don't correct the spanning exon
  if($have_spanning_exon) { 
    if(! $do_circular) { 
      DNAORG_FAIL("ERROR in $sub_name, found exon that spanned stop..start boundary, but we're not allowing circular genomes...", 1, $FH_HR); 
    }
    else { 
      # fix it
      checkForSpanningSequenceSegments($starts_AR, $stops_AR, $nexons_R, 1, $strand, $totlen); # 0 says: don't correct the spanning exon
    }
  }

  return;
}

#################################################################
# Subroutine: dumpInfoHashOfArrays()
# Incept:     EPN, Thu Feb 11 15:06:31 2016
# Synopsis:   Print an 'info' hash of arrays, for debugging purposes.
#
# Args:       $name2print:  name of hash of arrays
#             $by_array:    determines order of output, 
#                           '1' to print all keys for each array element, 
#                           '0' to print all array values for each key
#             $info_HAR:    ref of the hash of arrays to print
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
#################################################################
sub dumpInfoHashOfArrays { 
  my $sub_name = "dumpInfoHashOfArrays()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $by_array, $info_HAR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  # determine max number of array elements over all keys
  my $max_nel = 0;
  my $key;
  foreach $key (keys %{$info_HAR}) { 
    if(scalar(@{$info_HAR->{$key}}) > $max_nel) { 
      $max_nel = scalar(@{$info_HAR->{$key}});
    }
  }

  if($by_array) { 
    for(my $a_ctr = 0; $a_ctr < $max_nel; $a_ctr++) { 
      printf $FH ("el %2d\n", ($a_ctr+1)); 
      foreach $key (sort keys %{$info_HAR}) { 
        printf $FH ("\t$key: %s\n", (defined $info_HAR->{$key}[$a_ctr]) ? $info_HAR->{$key}[$a_ctr] : "undef"); 
      }
      printf $FH ("\n");
    }
  }
  else { # printing all array values for each key
    my $key_ctr = 1;
    my $nkeys = scalar(keys %{$info_HAR});
    
    foreach $key (sort keys %{$info_HAR}) { 
      my $nel = (@{$info_HAR->{$key}}) ? scalar(@{$info_HAR->{$key}}) : 0;
      printf $FH ("key %2d of %2d: $key, $nel elements\n", $key_ctr++, $nkeys);
      for(my $a_ctr = 0; $a_ctr < $nel; $a_ctr++) { 
        printf $FH ("\tel %2d: %s\n", ($a_ctr+1), (defined $info_HAR->{$key}[$a_ctr]) ? $info_HAR->{$key}[$a_ctr] : "undef"); 
      }
      printf $FH ("\n");
    }
  }
  
  return;
}

#################################################################
# Subroutine: validateAndGetSizeOfInfoHashOfArrays()
# Incept:     EPN, Thu Feb 11 15:06:40 2016
#
# Purpose:    Validate that the arrays in a hash of arrays are all the
#             same size, and return the number of keys (first
#             dimension) and size of all arrays (second dimension).
#
# Arguments:
#   $HAR:          REF to hash of array to validate
#   $reqd_keys_AR: REF to array of keys that must be in $HAR, else we die
#   $FH_HR:        REF to hash of file handles, including "log" and "cmd"
# 
# Returns: number of elements (scalar(@{$$HAR->{$key}})) for all keys $key
#          because we validate that all arrays are the same size, and die
#          if that is not true
#
# Dies:    if not all of the arrays (for all keys) are not the same size
#          OR if $reqd_keys_AR is defined and $HAR->{$key} does not exist for
#          >= 1 of the keys in $reqd_keys_AR.
#################################################################
sub validateAndGetSizeOfInfoHashOfArrays {
  my $sub_name = "validateAndGetSizeOfInfoHashOfArrays()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HAR, $reqd_keys_AR, $FH_HR) = (@_);
  
  #  make sure that all required keys exist
  if(defined $reqd_keys_AR && @{$reqd_keys_AR}) { 
    foreach my $reqd_key (@{$reqd_keys_AR}) { 
      if(! exists $HAR->{$reqd_key}) { 
        DNAORG_FAIL("ERROR in $sub_name, required key $reqd_key does not exist in the hash of arrays", 1, $FH_HR); 
      }
    }
  }

  # check consistency: each array should be the same size
  my $nel = -1;
  my $nkey = scalar(keys %{$HAR});
  my $nel_key = undef;
  foreach my $key (sort keys %{$HAR}) { 
    if($nel == -1) { 
      $nel = scalar(@{$HAR->{$key}});
      $nel_key = $key;
    }
    else { 
      if($nel != scalar(@{$HAR->{$key}})) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, expected number of elements in array for key $key is $nel (from key $nel_key) but %d exist!", scalar(@{$HAR->{$key}})), 1, $FH_HR);
      }
    }
    # make sure all values are defined
    for(my $i = 0; $i < $nel; $i++) { 
      if(! defined $HAR->{$key}[$i]) { 
        DNAORG_FAIL("ERROR in $sub_name, undefined value: key: $key index: $i", 1, $FH_HR); 
      }        
    }
  }

  return $nel;
}

#################################################################
# Subroutine: getLengthsAndCoords()
# Incept:     EPN, Thu Feb 11 13:32:34 2016
#
# Purpose:    For a given accession, retreive lengths and coorindate
#             strings of all features.
#
# Arguments:
#   $tbl_HAR:   REF to hash of arrays for accession we're interested in, PRE-FILLED
#   $len_AR:    REF to array to fill with lengths of features in %{$tbl_HAR}, FILLED HERE
#   $coords_AR: REF to array to fill with coordinates for each gene, FILLED HERE
#   $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void; fills @{$len_AR} and @{$coords_AR}
#
#################################################################
sub getLengthsAndCoords { 
  my $sub_name = "getLengthsAndCoords()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($tbl_HAR, $len_AR, $coords_AR, $FH_HR) = @_;

  my $ngenes = scalar(@{$tbl_HAR->{"coords"}});

  if ($ngenes > 0) { 
    for(my $i = 0; $i < $ngenes; $i++) { 
      push(@{$len_AR},    lengthFromCoords($tbl_HAR->{"coords"}[$i], $FH_HR));
      push(@{$coords_AR}, $tbl_HAR->{"coords"}[$i]);
    }
  }

  return;
}

#################################################################
# Subroutine: lengthFromCoords()
# Incept:     EPN, Thu Feb 11 15:28:20 2016
#
# Purpose:    Determine the length of a region from its coords in NCBI format.
#
# Arguments:
#   $coords:  the coords string in NCBI format
#   $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:    length in nucleotides implied by $coords  
#################################################################
sub lengthFromCoords { 
  my $sub_name = "lengthFromCoords()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($coords, $FH_HR) = @_;

  my @starts_A = ();
  my @stops_A  = ();
  my $nexons   = 0;

  startsStopsFromCoords($coords, \@starts_A, \@stops_A, \$nexons, $FH_HR);

  my $length = 0;
  for(my $i = 0; $i < $nexons; $i++) { 
    $length += abs($starts_A[$i] - $stops_A[$i]) + 1;
  }

  return $length;
}

#################################################################
# Subroutine: fetchSequencesUsingEslFetchCds()
# Incept:     EPN, Thu Feb 11 15:30:07 2016
# 
# Purpose:    Determine the length of a region give its coords in NCBI format.
#
# Arguments:
#   $esl_fetch_cds:    path to esl-fetch-cds executable
#   $do_force:         '1' to overwrite $out_fetch_file and $out_fasta_file if they exist,
#                      '0' to die if they exist
#   $out_fetch_file:   name of file to use as input to $esl_fetch_cds, created here
#   $out_fasta_file:   name of fasta file to create, created here
#   $do_circular:      '1' to duplicate genome, '0' not to
#   $accn_AR:          ref to array of accessions to fetch, PRE-FILLED
#   $totlen_HR:        ref to hash with total lengths of each sequence, PREF-FILLED
#   $seq_accn_AR:      ref to array of sequence names in newly created sequence file, FILLED HERE
#   $accn2seq_accn_HR: ref to hash, keys are accessions from @{$accn_AR}, values are
#                      sequence accessions from @{$seq_accn_AR}, can be undef if unwanted
#                      FILLED HERE (if defined)
#   $seq_accn2accn_HR: ref to hash, keys are accessions from @{$seq_accn_AR}, values are
#                      sequence accessions from @{$accn_AR}, can be undef if unwanted
#                      FILLED HERE (if defined)
#   $FH_HR:            REF to hash of file handles, including "log" and "cmd"
# 
# Returns:    void
#          
# Dies: if $esl_fetch_cds executable does not exist
#       if $esl_fetch_cds command fails
#################################################################
sub fetchSequencesUsingEslFetchCds { 
  my $sub_name = "fetchSequencesUsingEslFetchCds";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($esl_fetch_cds, $do_force, $fetch_file, $fasta_file, $do_circular, $accn_AR, $totlen_HR, $seq_accn_AR, $accn2seq_accn_HR, $seq_accn2accn_HR, $FH_HR) = @_;

  # contract check
  if(! -e $esl_fetch_cds) { 
    DNAORG_FAIL("ERROR in $sub_name, required executable $esl_fetch_cds does not exist", 1, $FH_HR); 
  }

  my @out_files_A = ($fetch_file, $fasta_file, $fasta_file.".ssi");
  if($do_force) { 
    foreach my $out_file (@out_files_A) { 
      if(-e $out_file) { 
        unlink $out_file;
      }
    }
  }
  else { 
    foreach my $out_file (@out_files_A) { 
      if(-e $out_file) { 
        DNAORG_FAIL("ERROR in $sub_name, output file $out_file already exists, and do_force option is off", 1, $FH_HR);
      }
    }
  }

  my $seq_accn;     # sequence accession
  my $fetch_string; # string output to the input file for esl-fetch-cds.pl
  
  # create the esl-fetch-cds input file
  open(OUT, ">", $fetch_file) || fileOpenFailure($fetch_file, $FH_HR, $!, "writing");
  my $naccn = scalar(@{$accn_AR});
  for(my $a = 0; $a < $naccn; $a++) { 
#  print OUT $accn_A[$a] . "\n";
    my $accn = $accn_AR->[$a];
    if(! exists $totlen_HR->{$accn}) { DNAORG_FAIL("ERROR no total length read for accession $accn", 1, $FH_HR); } 

    if($do_circular) { 
      $fetch_string = "join(" . $accn . ":1.." . $totlen_HR->{$accn} . "," . $accn . ":1.." . $totlen_HR->{$accn} . ")\n";
      print OUT $accn . ":" . "genome-duplicated" . "\t" . $fetch_string;
      $seq_accn = $accn . ":genome-duplicated:" . $accn . ":1:" . $totlen_HR->{$accn} . ":+:" . $accn . ":1:" . $totlen_HR->{$accn} . ":+:";
    }
    else { # $do_circular is FALSE
      $fetch_string = $accn . ":1.." . $totlen_HR->{$accn} . "\n";
      print OUT $accn . ":" . "genome" . "\t" . $fetch_string;
      $seq_accn = $accn . ":genome:" . $accn . ":1:" . $totlen_HR->{$accn} . ":+:";
    }
    if(defined $seq_accn_AR)      { push(@{$seq_accn_AR},        $seq_accn); }
    if(defined $accn2seq_accn_HR) { $accn2seq_accn_HR->{$accn} = $seq_accn;  }
    if(defined $seq_accn2accn_HR) { $seq_accn2accn_HR->{$seq_accn} = $accn;  }
  }
  close(OUT);
  
  # execute $esl_fetch_cds to fetch it
  my $cmd = "perl $esl_fetch_cds -nocodon $fetch_file > $fasta_file";
  runCommand($cmd, 0, $FH_HR);

  return;
}

#################################################################
# Subroutine: addNameAndBlankSsToStockholmAlignment()
# Incept:     EPN, Thu Feb 11 15:38:00 2016
#
# Synopsis:   Read in a stockholm alignment ($in_file), and 
#             add a name ($name) to it, then optionally add
#             a blank SS (#=GC SS_cons) annotation to it
#             and output a new file ($out_file) that is
#             identical to it but with the name annotation
#             and possibly a blank SS.
#
# Arguments:
#   $name:         name to add to alignment
#   $do_blank_ss:  '1' to add a blank SS, else '0'
#   $in_file:      input stockholm alignment
#   $out_file:     output stockholm alignment to create
#   $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    Two values:
#             alignment length 
#             checksum of alignment
#################################################################
sub addNameAndBlankSsToStockholmAlignment {
  my $sub_name = "addNameAndBlankSsToStockholmAlignment";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($name, $do_blank_ss, $in_file, $out_file, $FH_HR) = @_;

  if(! -s $in_file) { 
    DNAORG_FAIL("ERROR in $sub_name, file $in_file does not exist."); 
  }

  # open and validate file
  my $msa = Bio::Easel::MSA->new({
    fileLocation => $in_file,
                                 });  
  $msa->set_name($name);
  if($do_blank_ss) { 
    $msa->set_blank_ss_cons;
  }
  $msa->write_msa($out_file);

  return ($msa->alen, $msa->checksum);
}

#################################################################
# Subroutine: getQualifierValues()
# Incept:     EPN, Fri Feb 12 05:33:19 2016
#
# Purpose:    Retreive values for the qualifier $qualifier 
#             for accession $accn in the given %tbl_HHAR 
#             and return the values 
#             in $values_AR.
#
# Arguments:
#   $tbl_HHAR:  REF to hash of hash of arrays to retrieve from
#   $accn:      accession we're interested in
#   $qualifier: qualifier we're interested in (e.g. 'Product')
#   $values_AR: REF to array to fill with values of $qualifier
#               ADDED TO HERE
#   $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void; fills @{$values_AR}
# 
# Dies: if $accn doesn't exist in %{$tbl_HHAR}
#################################################################
sub getQualifierValues {
  my $sub_name = "getQualifierValues()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($tbl_HHAR, $accn, $qualifier, $values_AR, $FH_HR) = @_;

  if(! exists $tbl_HHAR->{$accn}) { DNAORG_FAIL("ERROR in $sub_name, no data for accession: $accn", 1, $FH_HR); }

  if(! @{$tbl_HHAR->{$accn}{$qualifier}}) { return; } # no annotation for $qualifier, do not update array

  my $nvalues = scalar(@{$tbl_HHAR->{$accn}{$qualifier}});

  if ($nvalues > 0) { 
    for(my $i = 0; $i < $nvalues; $i++) { 
      push(@{$values_AR}, $tbl_HHAR->{$accn}{$qualifier}[$i]);
    }
  }

  return;
}

#################################################################
# Subroutine: createCmDb()
# Incept:     EPN, Fri Feb 12 09:06:52 2016
# 
# Purpose:    Create a CM database from a stockholm database file
#             for use with Infernal 1.1 using the $cmbuild 
#             executable. If $cmcalibrate is defined, also run 
#             cmcalibrate. 
#
# Arguments:
#   $execs_HR:       reference to hash with infernal executables, 
#                    e.g. $execs_HR->{"cmbuild"} is path to cmbuild, PRE-FILLED
#   $do_calib_slow:  '1' to calibrate using default parameters instead of
#                    options to make it go much faster
#   $do_calib_local: '1' to run cmcalibrate locally, '0' to do it on the farm
#   $stk_file:       stockholm DB file
#   $out_root:       string for naming output files
#   $indi_name_AR:   ref to array of individual model names, we only use this if 
#                    $do_calib_local is 0 or undef, PRE-FILLED
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#                    
# Returns:    void
#
# Dies:       if $stk_file does not exist or is empty
#################################################################
sub createCmDb { 
  my $sub_name = "createCmDb()";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $do_calib_slow, $do_calib_local, $stk_file, $out_root, $indi_name_AR, $FH_HR) = @_;

  if(! -s $stk_file)  { DNAORG_FAIL("ERROR in $sub_name, $stk_file file does not exist or is empty", 1, $FH_HR); }

  # remove the binary files for the CM from an earlier cmbuild/cmpress, if they exist:
  for my $suffix ("i1m", "i1i", "i1f", "i1p") { 
    my $file = $out_root . ".cm." . $suffix;
    if(-e $file) { unlink $file; }
  }

  my ($cmbuild_opts,     $cmbuild_cmd);        # options and command for cmbuild
  my ($cmcalibrate_opts, $cmcalibrate_cmd);    # options and command for cmcalibrate
  my ($cmpress_opts,     $cmpress_cmd);        # options and command for cmpress

  my $cmbuild     = $execs_HR->{"cmbuild"};
  my $cmcalibrate = $execs_HR->{"cmcalibrate"};
  my $cmpress     = $execs_HR->{"cmpress"};
  my $cmfetch     = $execs_HR->{"cmfetch"};

  my $nmodel = scalar(@{$indi_name_AR});
  if($nmodel == 0) { DNAORG_FAIL("ERROR no model names in $sub_name", 1, $FH_HR); }

  # build step:
  $cmbuild_opts = "-F";
  $cmbuild_cmd  = "$cmbuild $cmbuild_opts $out_root.cm $stk_file > $out_root.cmbuild";
#  printf("%-71s ... ", sprintf("# Running cmbuild to build %d CMs", $nmodel));
####  printf("calling $cmbuild_cmd\n");
  my $secs_elapsed = runCommand($cmbuild_cmd, 0, $FH_HR);
#  printf("done. [%.1f seconds]\n", $secs_elapsed);

  # calibration step:
  $cmcalibrate_opts = " --cpu 4 ";
  if(! $do_calib_slow) { $cmcalibrate_opts .= " -L 0.04 "; }
  $cmcalibrate_cmd  = "$cmcalibrate $cmcalibrate_opts $out_root.cm > $out_root.cmcalibrate";
  
  if($do_calib_local) { 
    # calibrate the model locally
#    printf("%-71s ... ", "# Running cmcalibrate");
####    printf("calling $cmcalibrate_cmd\n");
    $secs_elapsed = runCommand($cmcalibrate_cmd, 0, $FH_HR);
#    printf("done. [%.1f seconds]\n", $secs_elapsed);

    # press the model
    $cmpress_cmd = "$cmpress $out_root.cm > $out_root.cmpress";
#    printf("%-71s ... ", "# Running cmpress");
    $secs_elapsed = runCommand($cmpress_cmd, 0, $FH_HR);
####printf("calling $cmpress_cmd\n");
#    printf("done [%.1f seconds]\n", $secs_elapsed);
  } # end of 'else' entered if $do_calib_farm is false
  else { 
    # run cmcalibrate on farm, one job for each CM file
    # split up model file into individual CM files, then submit a job to calibrate each one, and exit. 
    for(my $i = 0; $i < $nmodel; $i++) { 
      my $cmfetch_cmd = "$cmfetch $out_root.cm $indi_name_AR->[$i] > $out_root.$i.cm";
      runCommand($cmfetch_cmd, 0, $FH_HR);
      my $out_tail    = $out_root;
      $out_tail       =~ s/^.+\///;
      my $jobname     = "c." . $out_tail . $i;
      my $errfile     = $out_root . "." . $i . ".err";
      $cmcalibrate_cmd  = "$cmcalibrate $cmcalibrate_opts $out_root.$i.cm > $out_root.$i.cmcalibrate";
      my $farm_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n -l h_rt=288000,h_vmem=8G,mem_free=8G -pe multicore 4 -R y " . "\"" . $cmcalibrate_cmd . "\"\n";
####printf("calling $farm_cmd\n");
      runCommand($farm_cmd, 0, $FH_HR);
    }
    # final step, remove the master CM file, so we can create a new one after we're done calibrating
    runCommand("rm $out_root.cm", 0, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: matpeptValidateCdsRelationships()
# Incept:     EPN, Fri Feb 12 09:15:58 2016
#
# Purpose:    Validate the CDS:mat_peptide relationships 
#             in @{$cds2pmatpept_AAR}, which were probably 
#             ready from a file in parsed in pasreMatPeptSpecFile().
#
# Arguments:
#   $cds2pmatpept_AAR: REF to array of arrays, 1st dim value is cds_idx, value is array
#                      of mat_peptide indices that comprise that CDS, PRE-FILLED
#   $ref_cds_tbl_HAR:  REF to CDS table for reference genome
#   $ref_mp_tbl_HAR:   REF to mat_peptide table for reference genome
#   $do_circular:      '1' if we're allowing circular genomes
#   $ref_totlen:       length of reference genome
#   $FH_HR:            REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void 
#
# Dies:       if relationship is not valid
#
#################################################################
sub matpeptValidateCdsRelationships {
  my $sub_name = "matpeptValidateCdsRelationships()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cds2pmatpept_AAR, $ref_cds_tbl_HAR, $ref_mp_tbl_HAR, $do_circular, $ref_totlen, $FH_HR) = @_;

  # get CDS and mat_peptide length and coordinate strings
  my @ref_cds_len_A    = ();
  my @ref_cds_coords_A = ();
  getLengthsAndCoords($ref_cds_tbl_HAR, \@ref_cds_len_A, \@ref_cds_coords_A, $FH_HR);

  my @ref_mp_len_A    = ();
  my @ref_mp_coords_A = ();
  getLengthsAndCoords($ref_mp_tbl_HAR, \@ref_mp_len_A, \@ref_mp_coords_A, $FH_HR);
  
  # validate
  my $ref_ncds = scalar(@ref_cds_len_A);
  my $ref_nmp  = scalar(@ref_mp_len_A);
  my $ncds2check = scalar(@{$cds2pmatpept_AAR});
  my $prv_stop = undef; # previous mat_peptide's stop position, 
  for(my $cds_idx = 0; $cds_idx < scalar(@{$cds2pmatpept_AAR}); $cds_idx++) { 
    if(($cds_idx < 0) || ($cds_idx >= $ref_ncds)) { 
      DNAORG_FAIL("ERROR in $sub_name, cds_idx $cds_idx is out of range", 1, $FH_HR); 
    }
    my @cds_starts_A  = ();
    my @cds_stops_A   = ();
    my @cds_strands_A = ();
    my $cds_nexons    = 0;
    startsStopsStrandsFromCoordsLength($ref_cds_coords_A[$cds_idx], $ref_totlen, $do_circular, \@cds_starts_A, \@cds_stops_A, \@cds_strands_A, \$cds_nexons, $FH_HR);
    my $cds_start = $cds_starts_A[0];
    my $cds_stop  = $cds_stops_A[$cds_nexons-1];
    if($cds_nexons != 1) { 
      if($cds_nexons != 2) { 
        DNAORG_FAIL("ERROR in $sub_name, triple or more exon CDS broken up to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
      }
      if($cds_strands_A[0] ne $cds_strands_A[1]) { 
        DNAORG_FAIL("ERROR in $sub_name, double exon CDS with each exon on different strands to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
      }
      # two exon CDS, if any introns exist (any nt is not included between $cds_start..$cds_stop) then we can't handle it
      # example of a multi-'exon' CDS that we CAN handle is West Nile Virus CDS #2: NC_001563.2 join(97..3540,3540..3671)
      if($cds_strands_A[0] eq "+") { 
        if(($cds_starts_A[1] - $cds_stops_A[0] - 1) > 0) { 
          DNAORG_FAIL("ERROR in $sub_name, multiple exon CDS with an intron broken up to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
        }
      }
      else { 
        if(($cds_stops_A[0] - $cds_starts_A[1] - 1) > 0) { 
          DNAORG_FAIL("ERROR in $sub_name, multiple exon CDS with an intron broken up to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
        }
      }
    }

    # look at all mat_peptides that are supposed to comprise this CDS
    # and make sure that they do
    my $nmp2check = scalar(@{$cds2pmatpept_AAR->[$cds_idx]});
    for(my $x = 0; $x < $nmp2check; $x++) { 
      my $mp_idx = $cds2pmatpept_AAR->[$cds_idx][$x];
      if($mp_idx < 0 || $mp_idx >= $ref_nmp) { 
        DNAORG_FAIL("ERROR in $sub_name, mp_idx $mp_idx for cds_idx $cds_idx is out of range", 1, $FH_HR); 
      }
      my @mp_starts_A = ();
      my @mp_stops_A  = ();
      my $mp_nexons   = 0;
      startsStopsStrandsFromCoordsLength($ref_mp_coords_A[$mp_idx], $ref_totlen, $do_circular, \@mp_starts_A, \@mp_stops_A, undef, \$mp_nexons, $FH_HR);
      if($x == 0) { # verify start matches with CDS start
        if($mp_starts_A[0] != $cds_start) { 
          DNAORG_FAIL("ERROR in $sub_name, for cds_idx $cds_idx start of first mat_peptide doesn't match CDS start ($mp_starts_A[0] != $cds_start)", 1, $FH_HR); 
        }
      }
      if($x > 0) { # check that this mat_peptide is adjacent to previous one
        if($mp_starts_A[0] != ($prv_stop+1)) { 
          DNAORG_FAIL(sprintf("ERROR in $sub_name, for mat_peptides %d and %d are not adjacent (%d != %d+1)", $x, $x-1, $mp_starts_A[0], $prv_stop), 1, $FH_HR);
        }
      }
      if($x == ($nmp2check-1)) { # verify stop matches with CDS stop-3
        if(($mp_stops_A[($mp_nexons-1)]+3) != $cds_stop) { 
          DNAORG_FAIL(sprintf("ERROR in $sub_name, for cds_idx $cds_idx stop of final mat_peptide doesn't match CDS stop (%d != %d)", $mp_stops_A[($mp_nexons-1)], $cds_stop), 1, $FH_HR);
        }
      }
      $prv_stop = $mp_stops_A[($mp_nexons-1)];
      # printf("checked mp $mp_idx %d..%d\n", $mp_starts_A[0], $mp_stops_A[($mp_nexons-1)]);
    }
  }
  return;
}

#################################################################
# Subroutine: checkForSpanningSequenceSegments()
#
# Synopsis:   Check if two sequence segments (e.g. exons) are really 
#             just one that spans the stop/start boundary in a circular 
#             genome. (More explanatory notes in comments within
#             code of the subroutine.)
#
# Arguments:
#   $starts_AR: ref of array of start positions to potentially overwrite (if we find this is really only one exon)
#   $stops_AR:  ref of array of stop positions to potentially overwrite (if we find this is really only one exon)
#   $nexons_R:  ref to scalar of number of exons to overwrite (if we find this is really only one exon)
#   $do_update: '1' to update $starts_AR, $stops_AR and $nexons_R if we find two exons that span stop..start
#   $strand:    strand the exons are on
#   $totlen:    total length of the sequence
#           
# Returns:    '1' if we found two exons that spanned stop..start 
#             (and if ($do_update) then we also updated @{$starts_AR}, @{$stops_AR} and $$nexons_R)
# 
# Dies:       Never.
# 
#################################################################
sub checkForSpanningSequenceSegments {
  my $sub_name = "checkForSpanningSequenceSegments()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($starts_AR, $stops_AR, $nexons_R, $do_update, $strand, $totlen) = @_;

  my $found_spanning_exons = 0;
  if($$nexons_R == 2) { 
    # if we're in a circular genome, we need to check for a special case, where 
    # what looks like a 2-exon CDS is really a single exon that spans the stop..start boundary.
    # [Note that if the stop..start boundary is spanned by an intron (i.e. exon i is before stop,
    # and i+1 is after start) then we don't need to modify anything, we'll still fetch the proper
    # sequence even in a duplicated genome].
    #
    # Example 1: single exon that spans stop..start boundary on positive strand
    # join(2309..3182,1..1625) in a seq of length 3182, this should really be a single exon
    # $starts_A[0] = 2309
    # $stops_A[0]  = 3182
    # $starts_A[1] = 1
    # $stops_A[1]  = 1625
    # $nexons = 2;
    # 
    # should become:
    # $starts_A[0] = 2309
    # $stops_A[0]  = 4807
    # $nexons = 1;
    # 
    # Example 2: single exon that spans stop..start boundary on negative strand
    # complement(join(2309..3182,1..1625))   in a seq of length 3182, this should really be a single exon
    # $starts_A[0] = 3182
    # $stops_A[0]  = 2309
    # $starts_A[1] = 1625
    # $stops_A[1]  = 1
    # $nexons = 2;
    # 
    # should become:
    # $starts_A[0] = 4807
    # $stops_A[0]  = 2309
    # $nexons = 1;
    #
    # we can easily check and fix these cases:
    my $tmp_start = undef;
    my $tmp_stop  = undef;
    if($strand eq "+" && $stops_AR->[0] == $totlen && $starts_AR->[1] == 1) { 
      $tmp_start = $starts_AR->[0];
      $tmp_stop  = $stops_AR->[1] + $totlen;
    }
    elsif($strand eq "-" && $starts_AR->[0] == $totlen && $stops_AR->[1] == 1) { 
      my $tmp_start = $starts_AR->[1] + $totlen;
      my $tmp_stop  = $stops_AR->[0];
    }
    if(defined $tmp_start && defined $tmp_stop) { 
      @{$starts_AR} = ();
      @{$stops_AR} = ();
      $starts_AR->[0] = $tmp_start;
      $stops_AR->[0]  = $tmp_stop;
      $$nexons_R = 1;
      $found_spanning_exons = 1;
    }    
  }
  return $found_spanning_exons;
}


#################################################################
# Subroutine: edirectFtableOrMatPept2SingleFeatureTableInfo()
# Incept:     EPN, Fri Feb 12 09:48:00 2016
# 
# Purpose:    Given an edirect output feature table file or mature
#             peptide table file, parse it and store it's relevant
#             information in a 'single feature table' (%{$tbl_HHAR}).
#
#             This is a wrapper subroutine that calls 
#             parseEdirectMatPeptideFile() or 
#             parseEdirectFtableFile() and then 
#             getSingleFeatureTableInfo().
#
# Arguments:
#   $edirect_file:  name out edirect output file to parse
#   $do_matpept:    '1' if edirect file is a mature peptide file, '0' or undef if it is a ftable file
#   $feature:       the feature we want to store info on in $tbl_HHAR (e.g. "CDS" or "mat_peptide")
#   $tbl_HHAR:      REF to hash of hash of arrays we'll fill with info on $qual_name:
#                   1D: key: accession
#                   2D: key: qualifier (e.g. 'coords')
#                   3D: values for each qualifier, size will be number of features for this accession
#                       e.g. size of 5 means this accession has 5 CDS if $feature is CDS
#                       FILLED HERE
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if we're not able to parse the edirect file
#################################################################
sub edirectFtableOrMatPept2SingleFeatureTableInfo { 
  my $sub_name = "edirectFtableOrMatPept2SingleFeatureTableInfo()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($edirect_file, $do_matpept, $feature, $tbl_HHAR, $FH_HR) = @_;

  # variables relevant to parsing the edirect ftable and edirect mat_peptide files
  my %quals_HHA       = ();      # values that will go into a table
                                 # key 1: feature name, e.g. "CDS"
                                 # key 2: 'fac' string: <full_accession><fac_seq><coordinates>
                                 # value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>
  my @faccn_A         = ();      # array of all full accessions read
  my %fac_HHA         = ();      # used to easily determine list of 2D keys ('fac's) in quals_HHA 
                                 # for a given feature and faccn.
                                 # key 1: feature name: e.g. "CDS"
                                 # key 2: 'faccn', full accession
                                 # value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
  my %faccn2accn_H    = ();      # key: full accession, value: short accession, used for convenience
  my %column_HA       = ();      # key: feature, e.g. "CDS"
                                 # value: array of qualifiers that exist for that feature in at least 1 faccn
  my $dummy_column    = "DuMmY"; # name of 'dummy' column that we don't actually print, for use with -matpept only
  my %sep_H           = ();      # hash of 'separation' values we use to concatenate multiple strings into 
                                 # one, to ease with their storage
  $sep_H{"qnqv"}      = "!!";    # string that separates 'qualifier name' and and 'qualifier_values' 
  $sep_H{"fac"}       = "::";    # strings separating full accessions, start and stop coords
  $sep_H{"qval"}      = ";;";    # strings separating qualifier values  

  if(defined $do_matpept && $do_matpept) { 
    parseEdirectMatPeptideFile($edirect_file, $dummy_column, \%sep_H, \%quals_HHA, \@faccn_A, \%fac_HHA, \%faccn2accn_H, \%column_HA, $FH_HR);
  }
  else { # default, a ftable file
    parseEdirectFtableFile($edirect_file, undef, \%sep_H, \%quals_HHA, \@faccn_A, \%fac_HHA, \%faccn2accn_H, \%column_HA, $FH_HR);
  }
  if(! exists $quals_HHA{$feature}) { 
    DNAORG_FAIL("ERROR, in $sub_name, trying to create info HHA for feature $feature, but it doesn't exist in the feature table file $edirect_file", 1, $FH_HR);
  }
  getSingleFeatureTableInfo($dummy_column, \%sep_H, \%{$quals_HHA{$feature}}, \@faccn_A, \%{$fac_HHA{$feature}}, \%faccn2accn_H, \@{$column_HA{$feature}}, $tbl_HHAR, $FH_HR);
  
  return;
}

#################################################################
# Subroutine: getSingleFeatureTableInfo()
# Incept:     EPN, Fri Feb 12 09:51:48 2016
#
# Purpose:  Given data structures collected from
#           parseEdirecFtableFile() or parseEdirectMatPeptFile(), fill
#           a hash of hash of arrays (%{$tbl_HHAR}) with the
#           information for a specific feature ($feature_name,
#           e.g. 'CDS' or mat_peptide').
#
# Arguments:
#   $feature_name:   name of feature we want information for (e.g. CDS or mat_peptide)
#   $dummy_column:   name for dummy columns that we won't output, can be undef
#   $sep_HR:         ref to hash of 'separation' values, keys: "qnqv", "fac", and "qval"
#                    used to split values that have been concatenated with these 'separation'
#                    values in between, PRE-FILLED
#   $quals_HAR:      ref to hash of arrays with values we are printing, PRE-FILLED
#                    key 1: 'fac' string: <full_accession><fac_seq><coordinates>
#                    value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>    
#   $faccn_AR:       REF to array of full accessions read, PRE-FILLED
#   $fac_HAR:        REF to hash of arrays that is used to easily determine
#                    list of keys ('fac's) in quals_HA for a given feature and faccn
#                    key 1: 'faccn', full accession
#                    value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
#                    PRE-FILLED
#   $faccn2accn_HR:  REF to hash, key: full accession, value: short accession, 
#                    used for convenience in output function, PRE-FILLED
#   $column_AR:      REF to array of qualifiers for feature we're printing table for, PRE-FILLED
#   $tbl_HHAR:       REF to hash of hash of arrays we'll fill with info on $qual_name:
#                    1D: key: accession
#                    2D: key: qualifier (e.g. 'coords')
#                    3D: values for each qualifier, size will be number of features for this accession
#                    e.g. size of 5 means this accession has 5 CDS if $feature is CDS
#                    FILLED HERE
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if an passed in data structure is empty or contains invalid information
#################################################################
sub getSingleFeatureTableInfo { 
  my $sub_name = "getSingleFeatureTableInfo()";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($dummy_column, $sep_HR, $quals_HAR, $faccn_AR, $fac_HAR, $faccn2accn_HR, $column_AR, $tbl_HHAR, $FH_HR) = @_;

  my $start           = undef;   # start position
  my $stop            = undef;   # stop position
  my $fac             = undef;   # full accession and coords, accession and each segment coords separated by $fac_sep
  my $faccn           = undef;   # full accession, e.g. ref|NC_000883.2|
  my $faccn2          = undef;   # another full accession
  my $coords          = undef;   # coordinates
  my $qname           = undef;   # a qualifier name,  e.g. 'organism'
  my $qval            = undef;   # a qualifier value, e.g. 'Paramecia bursaria Chlorella virus 1'

  my $column;
  my $strand; 
  my $sort_coord;

  my $fac_sep  = $sep_HR->{"fac"};
  my $qnqv_sep = $sep_HR->{"qnqv"};
  my $qval_sep = $sep_HR->{"qval"};

  if(! %{$quals_HAR}) { DNAORG_FAIL("ERROR in $sub_name, quals hash of arrays is not defined", 1, $FH_HR); }
  if(! %{$fac_HAR})   { DNAORG_FAIL("ERROR in $sub_name, fac hash of arrays is not defined", 1, $FH_HR); }
  if(! @{$column_AR}) { DNAORG_FAIL("ERROR in $sub_name, column hash of arrays is not defined", 1, $FH_HR); }

  # go through all full accessions in order they were read from feature table
  foreach $faccn (@{$faccn_AR}) { 
    if(exists $fac_HAR->{$faccn}) { # if this accession has >= 1 qualifiers for this feature
      my $accn = $faccn2accn_HR->{$faccn};
      printf("accn: $accn\n");
      stripVersion(\$accn);
      if(%{$tbl_HHAR} && %{$tbl_HHAR->{$accn}}) { 
        DNAORG_FAIL("ERROR read data twice in $sub_name for accession $accn", 1, $FH_HR);
      }
      foreach $fac (@{$fac_HAR->{$faccn}}) { # foreach 'fac', accession + set of coords

        ($faccn2, $coords, $sort_coord, $strand) = helperBreakdownFac($fac, $fac_sep, $FH_HR);
        if($faccn ne $faccn2) { DNAORG_FAIL("ERROR inconsistent fac value: $faccn ne $faccn2", 1, $FH_HR); }
        
        if(exists $quals_HAR->{$fac}) { # if there's any qualifiers for this fac
          # printf("quals_HA feature: fac: $fac exists!\n"); 

          push(@{$tbl_HHAR->{$accn}{"coords"}},         $coords);
          push(@{$tbl_HHAR->{$accn}{"strand"}},         $strand);
          push(@{$tbl_HHAR->{$accn}{"min-coord"}},      $sort_coord);
          
          # for all columns in the table
          foreach $column (@{$column_AR}) {
            if((! defined $dummy_column) || ($column ne $dummy_column)) { 
              my $save_str = ""; 
              
              # for all qualifier names and values 
              foreach my $qnqv (@{$quals_HAR->{$fac}}) { 
                ($qname, $qval) = split($qnqv_sep, $qnqv);
                ### printf("faccn: $faccn qnqv: $qnqv split into $qname $qval\n");
                  
                # if this qname matches this column, then it's the appropriate value to save here
                if($qname eq $column) { 
                  if($save_str eq "") { # first value in this cell
                    $save_str = $qval;  
                  }
                  else { 
                    if($qval =~ m/\Q$qval_sep/) { DNAORG_FAIL("ERROR qualifier_name $qval has the string $qval_sep in it", 1, $FH_HR); }
                    $save_str .= $qval_sep . $qval; # not first value, concatenate onto previous values
                  }
                }
              }
              # if there's no value for this qualifier, put '-'
              if($save_str eq "") { $save_str = "-"; }
              push(@{$tbl_HHAR->{$accn}{$column}}, $save_str);
            }
          } 
        }
      }
    }
  }
  
  return;
}

#####################################################################
# Subroutine: outputBanner()
# Incept:     EPN, Thu Oct 30 09:43:56 2014 (rnavore)
# 
# Purpose:    Output the dnaorg banner.
#
# Arguments: 
#    $FH:                file handle to print to
#    $synopsis:          string reporting the date
#    $command:           command used to execute the script
#    $date:              date information to print
#    $opts_used:         string describing options used 
#
# Returns:    Nothing, if it returns, everything is valid.
# 
# Dies: never
####################################################################
sub outputBanner {
  my $nargs_expected = 5;
  my $sub_name = "outputBanner()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($FH, $synopsis, $command, $date, $opts_used) = @_;

  print $FH ("\# $synopsis\n");
  print $FH ("\# dnaorg 0.1 (February 2016)\n");
#  print $FH ("\# Copyright (C) 2014 HHMI Janelia Research Campus\n");
#  print $FH ("\# Freely distributed under the GNU General Public License (GPLv3)\n");
  print $FH ("\# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  if(defined $command) { print $FH ("# command: $command\n"); }
  if(defined $date)    { print $FH ("# date:    $date\n"); }
  if(defined $opts_used && $opts_used ne "") { 
    print $FH $opts_used;
  }
  printf $FH ("#\n");

  return;
}

#################################################################
# Subroutine : DNAORG_FAIL()
# Incept:      EPN, Wed Nov 11 06:22:59 2009 (rnavore)
#
# Purpose:     Print an error message to STDERR and sum and 
#              log files in $FH_HR (ref to hash of file handles)
#              then close those file handles and exit.
#
# Arguments: 
#   $errmsg:       the error message to write
#   $status:       error status to exit with
#   $FH_HR:        ref to hash of file handles, including "log" and "cmd"
# 
# Returns:     Nothing, this function will exit the program.
#
################################################################# 
sub DNAORG_FAIL { 
  my $nargs_expected = 3;
  my $sub_name = "DNAORG_FAIL()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, DNAORG_FAIL() entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }
  my ($errmsg, $status, $FH_HR) = @_;
  
  if($errmsg !~ m/\n$/) { $errmsg .= "\n\n"; }
  else                  { $errmsg .= "\n"; }
  if($errmsg !~ m/^\n/) { $errmsg = "\n" . $errmsg; }
  
  if(defined $FH_HR) { 
    my $cmd_FH = $FH_HR->{"cmd"};
    my $log_FH = $FH_HR->{"log"};
    if(defined $cmd_FH) { 
      print $cmd_FH $errmsg;
      print $cmd_FH "# DNAORG-FAILURE\n";
    }
    if(defined $log_FH) {
      print $log_FH $errmsg;
      print $log_FH "# DNAORG-FAILURE\n";
    }
    # close each file handle
    foreach my $key (keys %{$FH_HR}) { 
      close $FH_HR->{$key};
    }
  }
  
  printf STDERR $errmsg; 
  exit($status);
}

#################################################################
# Subroutine: getStrandStats()
# Incept:     EPN, Thu Feb 11 15:14:09 2016
# 
# Purpose:    Retrieve strand stats from a tbl_HHA.
#
# Arguments:
#   $tbl_HHAR:  ref to hash of hash of arrays
#   $accn:      1D key to get strand info for
#   $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:    6 values:
#             $nfeatures:  number of features
#             $npos:       number of genes with all segments on positive strand
#             $nneg:       number of genes with all segments on negative strand
#             $nunc:       number of genes with all segments on unknown strand 
#             $nbth:       number of genes with that don't fit above 3 categories
#             $strand_str: strand string, summarizing strand of all genes, in order
#
# Dies: if 'strand' doesn't exist as a key in $tbl_HHAR
#       if we can't parse a strand value
#################################################################
sub getStrandStats {
  my $sub_name = "getStrandStats()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($tbl_HHAR, $accn, $FH_HR) = @_;

  my $nfeatures; # number of genes in this genome
  my $npos = 0;  # number of genes on positive strand 
  my $nneg = 0;  # number of genes on negative strand 
  my $nbth = 0;  # number of genes with >= 1 segment on both strands (usually 0)
  my $nunc = 0;  # number of genes with >= 1 segments that are uncertain (usually 0)
  my $strand_str = "";

  if(! exists $tbl_HHAR->{$accn}{"strand"}) { DNAORG_FAIL("ERROR in $sub_name, didn't read strand information for accn: $accn", 1, $FH_HR); }

  $nfeatures = scalar(@{$tbl_HHAR->{$accn}{"coords"}});
  if ($nfeatures > 0) { 
    for(my $i = 0; $i < $nfeatures; $i++) { 

      if   ($tbl_HHAR->{$accn}{"strand"}[$i] eq "+") { $npos++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "-") { $nneg++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "!") { $nbth++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "?") { $nunc++; }
      else { DNAORG_FAIL(sprintf("ERROR in $sub_name, unable to parse strand for feature %d for $accn\n", $i+1), 1, $FH_HR); }
      $strand_str .= $tbl_HHAR->{$accn}{"strand"}[$i];
    }
  }

  return ($nfeatures, $npos, $nneg, $nunc, $nbth, $strand_str);
}

#################################################################
# Subroutine : fileOpenFailure()
# Incept:      EPN, Wed Nov 11 05:39:56 2009 (rnavore)
#
# Purpose:     Called if a open() call fails on a file.
#              Print an informative error message
#              to $FH_HR->{"cmd"} and $FH_HR->{"log"}
#              and to STDERR, then exit with <$status>.
#
# Arguments: 
#   $filename: file that we couldn't open
#   $status:   error status
#   $action:   "reading", "writing", "appending"
#   $FH_HR:    ref to hash of open file handles to close
# 
# Returns:     Nothing, this function will exit the program.
#
################################################################# 
sub fileOpenFailure { 
  my $nargs_expected = 4;
  my $sub_name = "fileOpenFailure()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $status, $action, $FH_HR) = @_;

  my $errmsg;
  if(($action eq "reading") && (! (-e $filename))) { 
    $errmsg = "\nERROR, could not open $filename for reading. It does not exist.\n";
  }
  else { 
    $errmsg = "\nERROR, could not open $filename for $action.\n";
  }
  DNAORG_FAIL($errmsg, $status, $FH_HR);
}

###########################################################################
# the next line is critical, a perl module must return a true value
return 1;
###########################################################################
