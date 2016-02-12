#!/usr/bin/perl
#
# dnaorg_notused.pm
# Eric Nawrocki
# EPN, Thu Feb 11 15:17:16 2016
# 
# Subroutines written for the dnaorg project, that are currently
# not called by dnaorg_build.pl or dnaorg_annotate.pl.
#
# Kept here for reference, and for possible future use.
#
# List of subroutines in this file:
#
# getStrandStats():            return stats on strands for a given accession
# outputSingleFeatureTable():  output a file with information on a single feature obtained from an edirect file
# parseSingleFeatureTableFile():           parse a single feature table output from outputSingleFeatureTable()
# dumpHashOfHashOfArrays():                print the contents of a hash of hash of arrays, possibly for debugging purposes
# compareTwoHahsOfHashOfArrays():          compare two hash of hash of arrays

use strict;
use warnings;

#################################################################
# Subroutine: getStrandStats()
# Incept:     EPN, Thu Feb 11 15:14:09 2016
# 
# Purpose:    Retrieve strand stats from a tbl_HHA.
#
# Arguments:
#   $tbl_HHAR:  ref to hash of hash of arrays
#   $accn:      1D key to get strand info for
#   $FH_HR:     REF to hash of file handles, including "log" and "sum"
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
# Subroutine: outputSingleFeatureTable()
# Incept:     EPN, Thu Feb 11 14:01:30 2016
#
# Purpose: Given data structures collected from
#          parseEdirecFtableFile() or parseEdirectMatPeptideFile(),
#          output a 'single feature' table with all information for
#          a single 'feature', e.g. 'CDS' or 'mat_peptide'.
#
# Args:       $outfile:        name of output file to create
#             $dummy_column:   name for dummy columns that we won't output, can be undef
#             $sep_HR:         ref to hash of 'separation' values, keys: "qnqv", "fac", and "qval"
#                              used to split values that have been concatenated with these 'separation'
#                              values in between, PRE-FILLED
#             $quals_HAR:      ref to hash of arrays with values we are printing, PRE-FILLED
#                              key 1: 'fac' string: <full_accession><fac_seq><coordinates>
#                              value: array of 'qnqv' strings: <qualifier_name><qnqv_sep><qualifier_value>    
#             $faccn_AR:       REF to array of full accessions read, PRE-FILLED
#             $fac_HAR:        REF to hash of arrays that is used to easily determine
#                              list of keys ('fac's) in quals_HA for a given feature and faccn
#                              key 1: 'faccn', full accession
#                              value: array of 'fac' strings: <full_accession><fac_seq><coordinates>
#                              PRE-FILLED
#             $faccn2accn_HR:  REF to hash, key: full accession, value: short accession, 
#                              used for convenience in output function, PRE-FILLED
#             $column_AR:      REF to array of qualifiers for feature we're printing table for
#             $FH_HR:          REF to hash of file handles, including "log" and "sum"
#
# Returns:    void
#
# Dies:       if we can't open $outfile, or we see something inconsistent in the input
#################################################################
sub outputSingleFeatureTable {
  my $sub_name = "outputSingleFeatureTable()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($outfile, $dummy_column, $sep_HR, $quals_HAR, $faccn_AR, $fac_HAR, $faccn2accn_HR, $column_AR) = @_;

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

  if(! %{$quals_HAR}) { DNAORG_FAIL("ERROR in $sub_name, quals hash of arrays is not defined",  1, $FH_HR); }
  if(! %{$fac_HAR})   { DNAORG_FAIL("ERROR in $sub_name, fac hash of arrays is not defined",    1, $FH_HR); }
  if(! @{$column_AR}) { DNAORG_FAIL("ERROR in $sub_name, column hash of arrays is not defined", 1, $FH_HR); }

  open(OUT, ">", $outfile) || fileOpenFailure($outfile, $FH_HR, $!, "writing");

  # print header line with names of column
  print OUT "#full-accession\taccession\tcoords\tstrand\tmin-coord";
  foreach $column (@{$column_AR}) { 
    if((! defined $dummy_column) || ($column ne $dummy_column)) { 
      print OUT "\t$column"; 
    }
  }
  print OUT "\n";

  # go through all full accessions in order they were read from feature table
  foreach $faccn (@{$faccn_AR}) { 
    if(exists $fac_HAR->{$faccn}) { # if this accession has >= 1 qualifiers for this feature
      foreach $fac (@{$fac_HAR->{$faccn}}) { # foreach 'fac', accession + set of coords
        my @output_A = ();
        ($faccn2, $coords, $sort_coord, $strand) = parseHelperBreakdownFac($fac, $fac_sep, $FH_HR);
        if($faccn ne $faccn2) { DNAORG_FAIL("ERROR in $sub_name, inconsistent fac value: $faccn ne $faccn2", 1, $FH_HR); }
        
        if(exists $quals_HAR->{$fac}) { # if there's any qualifiers for this fac
          # printf("quals_HA feature: fac: $fac exists!\n"); 
          print OUT $faccn. "\t" . $faccn2accn_HR->{$faccn} . "\t" . $coords . "\t" . $strand . "\t" . $sort_coord;
          
          # for all columns in the table
          foreach $column (@{$column_AR}) {
            if((! defined $dummy_column) || ($column ne $dummy_column)) { 
              ### printf("\n\n");
              my $column_str = ""; 
              
              # for all qualifier names and values 
              foreach my $qnqv (@{$quals_HAR->{$fac}}) { 
                ($qname, $qval) = split($qnqv_sep, $qnqv);
                ### printf("faccn: $faccn qnqv: $qnqv split into $qname $qval\n");
                  
                # if this qname matches this column, then it's the appropriate value to output here
                if($qname eq $column) { 
                  if($column_str eq "") { # first value in this cell
                    $column_str = $qval;  
                  }
                  else { 
                    if($qval =~ m/\Q$qval_sep/) { DNAORG_FAIL("ERROR in $sub_name, qualifier_name $qval has the string $qval_sep in it", 1, $FH_HR); }
                    $column_str .= $qval_sep . $qval; # not first value, concatenate onto previous values
                  }
                }
              }
              # if there's no value for this qualifier, put '-'
              if($column_str eq "") { $column_str = "-"; }
              print OUT "\t$column_str";
            }
          }
          print OUT "\n";
        }
      }
    }
  }
  close(OUT);
  
  return;
}
#################################################################
# Subroutine: parseSingleFeatureTable()
# Incept:     EPN, Thu Feb 11 14:27:17 2016
#
# Synopsis:   Parses a single feature table file and stores the 
#             relevant info in it in $values_HHAR.
#
# Arguments:
#   $tblfile:      full path to a table file to parse
#   $values_HHAR:  REF to hash of hash of arrays to, FILLED HERE
#   $FH_HR:        REF to hash of file handles, including "log" and "sum"
#
# Returns:    void; fills @{$values_HHAR}
#
# Dies:       if there's problem parsing $tblfile
#################################################################
sub parseSingleFeatureTable {
  my $sub_name = "parseSingleFeatureTable()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 

  my ($tblfile, $values_HHAR) = @_;

  ##full-accession	accession	coords	strand	min-coord	gene
  #gb|HM448898.1|	HM448898.1	129..476	+	129	AV2

  open(TBL, $tblfile) || fileOpenFailure($tblfile, $FH_HR, $!, "reading");

  # get column header line:
  my $line_ctr = 0;
  my @colnames_A = ();
  my $line = <TBL>;
  my $ncols = undef;
  $line_ctr++;
  if(! defined $line) { DNAORG_FAIL("ERROR in $sub_name, did not read any lines from file $tblfile", 1, $FH_HR); }
  chomp $line;
  if($line =~ s/^\#//) { 
    @colnames_A = split(/\t/, $line);
    $ncols = scalar(@colnames_A);
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, first line of $tblfile did not start with \"#\"", 1, $FH_HR);
  }
  if($colnames_A[0] ne "full-accession") { DNAORG_FAIL("ERROR in $sub_name, first column name is not full-accession", 1, $FH_HR); }
  if($colnames_A[1] ne "accession")      { DNAORG_FAIL("ERROR in $sub_name, second column name is not accession",     1, $FH_HR); }
  if($colnames_A[2] ne "coords")         { DNAORG_FAIL("ERROR in $sub_name, third column name is not coords",         1, $FH_HR); }

  # read remaining lines
  while($line = <TBL>) { 
    chomp $line;
    $line_ctr++;
    if($line =~ m/^\#/) { DNAORG_FAIL("ERROR, line $line_ctr of $tblfile begins with \"#\"", 1, $FH_HR); }
    my @el_A = split(/\t/, $line);
    if(scalar(@el_A) != $ncols) { 
      DNAORG_FAIL("ERROR in $sub_name, read wrong number of columns in line $line_ctr of file $tblfile", 1, $FH_HR);
    }
    my $prv_min_coord = 0;
    # get accession
    my $accn = $el_A[1]; 
    stripVersion(\$accn);
    if(! exists $values_HHAR->{$accn}) { 
      %{$values_HHAR->{$accn}} = (); 
    }

    for(my $i = 0; $i < $ncols; $i++) { 
      my $colname = $colnames_A[$i];
      my $value   = $el_A[$i]; 
     if($colname eq "min-coord") { 
        if($value < $prv_min_coord) { 
          DNAORG_FAIL("ERROR, minimum coordinates out of order at line $line_ctr and previous line of file $tblfile", 1, $FH_HR); 
        }
        $prv_min_coord = $value; 
        # printf("prv_min_coord: $prv_min_coord\n");
      }
      
      if(! exists $values_HHAR->{$accn}{$colname}) { 
        @{$values_HHAR->{$accn}{$colname}} = ();
      }
      push(@{$values_HHAR->{$accn}{$colname}}, $el_A[$i]);
      #printf("pushed $accn $colname $el_A[$i]\n");
    }
  }
  close(TBL);
  return;
}
#################################################################
# Subroutine: dumpHashOfHashOfArrays()
# Incept:     EPN, Fri Feb 12 09:38:01 2016
#
# Purpose:    Print a hash of hash of arrays, possibly for 
#             debugging purposes.
# 
# Arguments:
#   $name2print: name to print for data structure, can be undef
#   $HHAR:       ref of the hash of hash of arrays to print
#   $FH:         file handle to print to 
#
# Returns:    void
#
# Dies:       Never.
#################################################################
sub dumpHashOfHashOfArrays { 
  my $sub_name = "dumpHashOfHashOfArrays()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
 
  my ($name2print, $HHAR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $key1;
  my $key2;
  my $i = 0;
  foreach $key1 (sort keys (%{$HHAR})) { 
    foreach $key2 (sort keys (%{$HHAR->{$key1}})) { 
      my $nel = scalar(@{$HHAR->{$key1}{$key2}});
      for($i = 0; $i < $nel; $i++) { 
        printf $FH ("HHA{$key1}{$key2}{$i] = $HHAR->{$key1}{$key2}[$i]\n");
      }
      printf $FH ("\n");
    }
    printf $FH ("\n");
  }
  return;
}

#################################################################
# Subroutine: compareTwoHashOfHashOfArrays()
# Incept:     EPN, Fri Feb 12 09:40:24 2016
#
# Purpose:    Compare two hash of hash of arrays and print any
#             differences (if $FH is defined).
#
# Arguments:
#   $name1:    name to print for first HHA
#   $name2:    name to print for second HHA
#   $HHAR1:    ref to first  hash of hash of arrays to print
#   $HHAR2:    ref to second hash of hash of arrays to print
#   $FH:       if defined, a file handle to print any 
#              differences to
#   $FH_HR:    REF to hash of file handles, including "log" and "sum"
#              differences to
#
# Returns:    '1' if the two hash of hash of arrays are the same, else '0'
#
# Dies:        if HHAR1 or HHAR2 is undefined
#################################################################
sub compareTwoHashOfHashOfArrays { 
  my $sub_name = "compareTwoHashOfHashOfArrays()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $narg_expected); exit(1); } 
 
  my ($name1, $name2, $HHAR1, $HHAR2, $FH, $FH_HR) = @_;

  if(! %{$HHAR1}) { DNAORG_FAIL("ERROR in $sub_name, first  hash of hash of arrays is undefined", 1, $FH_HR); }
  if(! %{$HHAR2}) { DNAORG_FAIL("ERROR in $sub_name, second hash of hash of arrays is undefined", 1, $FH_HR); }

  if(defined $FH && (defined $name1 || defined $name2)) { 
    printf $FH ("in $sub_name, comparing %s and %s\n", (defined $name1) ? $name1 : "undef", (defined $name2) ? $name2 : "undef");
  }

  my $key1;
  my $key2;
  my $i = 0;
  my $found_diff = 0;

  # look for any data in HHAR1 but not in HHAR2
  foreach $key1 (sort keys (%{$HHAR1})) { 
    if(! exists $HHAR2->{$key1}) {
      if(defined $FH) { printf $FH ("1st dim key $key1 is in data structure 1 but not in data structure 2\n"); }
      $found_diff = 1;
    }
    else { 
      foreach $key2 (sort keys (%{$HHAR1->{$key1}})) { 
        if(! exists $HHAR2->{$key1}{$key2}) {
          if(defined $FH) { printf $FH ("1st/2nd dim key combo {$key1}{$key2} is in data structure 1 but not in data structure 2\n"); }
          $found_diff = 1;
        }
        else { 
          my $nel1 = scalar(@{$HHAR1->{$key1}{$key2}});
          my $nel2 = scalar(@{$HHAR2->{$key1}{$key2}});
          for($i = 0; $i < $nel1; $i++) { 
            if($i > $nel2 || $HHAR1->{$key1}{$key2}[$i] ne $HHAR1->{$key1}{$key2}[$i]) { 
              if(defined $FH) { printf $FH ("element {$key1}{$key2}[$i]: $HHAR1->{$key1}{$key2}[$i] != %s\n", ($i < $nel2) ? $HHAR2->{$key1}{$key2}[$i] : "(undef)"); }
              $found_diff = 1;
            }
          }
        }
      }
    }
  }

  # now look for any data in HHAR2 but not in HHAR1
  foreach $key1 (sort keys (%{$HHAR2})) { 
    if(! exists $HHAR1->{$key1}) {
      if(defined $FH) { printf $FH ("2nd dim key $key1 is in data structure 2 but not in data structure 1\n"); }
      $found_diff = 1;
    }
    else { 
      foreach $key2 (sort keys (%{$HHAR2->{$key1}})) { 
        if(! exists $HHAR1->{$key1}{$key2}) {
          if(defined $FH) { printf $FH ("1st/2nd dim key combo {$key1}{$key2} is in data structure 2 but not in data structure 1\n"); }
          $found_diff = 1;
        }
        # if we get there then @{$HHAR1->{$key1}{$key2}} and @{$HHAR2->{$key1}{$key2}} both exist and 
        # we've already checked above whether they are identical or not, so we're done
      }
    }
  }

  return ($found_diff) ? 0 : 1;
}
