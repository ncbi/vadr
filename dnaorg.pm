#!/usr/bin/perl
# 
# version: 0.45 [Feb 2019]
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

#########################################################################################
#
# Subroutines related to features or segments:
# dng_FeatureInfoImputeCoords
# dng_FeatureInfoImputeLength
# dng_FeatureInfoImputeSourceIdx
# dng_FeatureInfoImputeParentIdx
# dng_FeatureInfoImputeOutname()
# dng_FeatureInfoImpute3paFtrIdx
# dng_FeatureInfoStartStopStrandArrays()
# dng_FeatureInfoCountType
# dng_FeatureInfoValidateCoords
# dng_FeatureInfoChildrenArrayOfArrays()
#
# dng_SegmentInfoPopulate()
# 
# dng_FeatureTypeAndTypeIndexString()
# dng_FeatureTypeIsCds()
# dng_FeatureTypeIsMatPeptide()
# dng_FeatureTypeIsGene()
# dng_FeatureTypeIsCdsOrMatPeptide()
# dng_FeatureTypeIsCdsOrMatPeptideOrGene()
# dng_FeatureChildrenArray()
# dng_FeatureNumSegments()
# dng_FeatureIsDuplicate()
# dng_Feature5pMostPosition()
# dng_Feature3pMostPosition()
# dng_FeatureSummarizeSegment()
# dng_FeatureStartStopStrandArrays()
# dng_FeatureSummaryStrand
# dng_FeaturePositionSpecificValueBreakdown()
# 
# Subroutines related to alerts:
# dng_AlertInfoInitialize()
# dng_AlertInfoAdd()
# dng_AlertInfoSetFTableInvalidatedBy
# 
# Subroutines related to parallelization on the compute farm:
# dng_SubmitJob()
# dng_WaitForFarmJobsToFinish()
#
# Subroutines related to sequence and model coordinates: 
# dng_CoordsTokenParse()
# dng_CoordsLength()
# dng_CoordsFromLocation()
# dng_CoordsComplement()
#
# Subroutines related to eutils:
# dng_EutilsFetchToFile()
# dng_EutilsFetchUrl()
# 
# Subroutines related to model info files:
# dng_ModelInfoFileWrite()
# dng_ModelInfoFileParse()
# 
# Miscellaneous subroutines:
# dng_SplitFastaFile()
# dng_SplitNumSeqFiles()
# dng_CdsFetchStockholmToFasta()
# dng_StripVersion()
#
#################################################################
# Subroutine: dng_FeatureInfoImputeCoords
# Incept:     EPN, Wed Mar 13 13:15:33 2019
# 
# Purpose:    Fill "coords" values in %{$ftr_info_AHR}
# 
# Arguments:
#   $ftr_info_AHR:   REF to feature information, added to here
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub dng_FeatureInfoImputeCoords { 
  my $sub_name = "dng_FeatureInfoImputeCoords";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "location"
  my @keys_A = ("type", "location");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"coords"} = dng_CoordsFromLocation($ftr_info_AHR->[$ftr_idx]{"location"}, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: dng_FeatureInfoImputeLength
# Incept:     EPN, Thu Mar 14 12:07:16 2019
# 
# Purpose:    Fill "length" values in @{$ftr_info_AHR}
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub dng_FeatureInfoImputeLength { 
  my $sub_name = "dng_FeatureInfoImputeSourceIdx";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  # go through all features and determine length by parsing the 
  # "coords" value
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $len = 0;
    my @sgm_start_A  = (); # array of starts, one per segment
    my @sgm_stop_A   = (); # array of stops, one per segment
    dng_FeatureStartStopStrandArrays($ftr_info_AHR->[$ftr_idx]{"coords"}, \@sgm_start_A, \@sgm_stop_A, undef, $FH_HR);
    for(my $s = 0; $s < scalar(@sgm_start_A); $s++) { 
      $len += abs($sgm_start_A[$s] - $sgm_stop_A[$s]) + 1;
    }
    $ftr_info_AHR->[$ftr_idx]{"length"} = $len;
  }

  return;
}

#################################################################
# Subroutine: dng_FeatureInfoImputeSourceIdx
# Incept:     EPN, Wed Mar 13 13:20:01 2019
# 
# Purpose:    Fill "source_idx" values in @{$ftr_info_AHR}
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub dng_FeatureInfoImputeSourceIdx { 
  my $sub_name = "dng_FeatureInfoImputeSourceIdx";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  # go through all features and determine duplicates (set 'source_idx')
  # 
  # $ftr_info_AHR->{"source_idx"}[$ftr_idx] set to $ftr_idx2 if
  # - $ftr_idx type is gene
  # - $ftr_idx2 type is CDS
  # - $ftr_idx and $ftr_idx2 have identical coords (for all segments)
  #
  # else "-1" if no $ftr_idx2 exists for $ftr_idx that satisfies above
  #
  # dies if more than one $ftr_idx2 satisfies above
  my ($ftr_idx, $ftr_idx2); # feature indices
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"source_idx"} = $ftr_idx; # initialize
    for($ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
      if($ftr_idx != $ftr_idx2) {
        if(($ftr_info_AHR->[$ftr_idx]{"type"}  eq "gene") && 
           ($ftr_info_AHR->[$ftr_idx2]{"type"} eq "CDS") && 
           ($ftr_info_AHR->[$ftr_idx]{"coords"} eq $ftr_info_AHR->[$ftr_idx2]{"coords"})) { 
          if($ftr_info_AHR->[$ftr_idx]{"source_idx"} != $ftr_idx) { 
            ofile_FAIL(sprintf("ERROR in $sub_name, unable to determine source (two choices) for duplicate feature of type %s and coords %s\n", 
                                $ftr_info_AHR->[$ftr_idx]{"type"}, $ftr_info_AHR->[$ftr_idx]{"coords"}), "dnaorg", 1, $FH_HR);
          }
          $ftr_info_AHR->[$ftr_idx]{"source_idx"} = $ftr_idx2;
        }
      }
    }
  }

  return;
}

#################################################################
# Subroutine: dng_FeatureInfoImputeParentIdx
# Incept:     EPN, Wed Mar 13 13:33:33 2019
# 
# Purpose:    Fill "parent_idx" values in @{$ftr_info_AHR}
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub dng_FeatureInfoImputeParentIdx {
  my $sub_name = "dng_FeatureInfoImputeParentIdx";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  # go through all features and determine parents (set 'parent_idx')
  # 
  # $ftr_info_AHR->{"parent_idx"}[$ftr_idx] set to $ftr_idx2 if:
  # - $ftr_idx type is mat_peptide
  # - $ftr_idx2 type is CDS
  # - start of $ftr_idx2 is 5' of or equal to start of $ftr_idx
  # - stop  of $ftr_idx2 is 3' of or equal to stop of $ftr_idx
  # - $ftr_idx and $ftr_idx2 are both "+" or both "-" strands
  # 
  # else "-1" if no $ftr_idx2 exists for $ftr_idx that satisfies above
  #
  # dies if more than one $ftr_idx2 satisfies above
  my ($ftr_idx, $ftr_idx2); # feature indices
  my $ftr_5p_pos;  # 3'-most position for feature $ftr_idx
  my $ftr_3p_pos;  # 3'-most position for feature $ftr_idx
  my $ftr_5p_pos2; # 5'-most position for feature $ftr_idx2
  my $ftr_3p_pos2; # 5'-most position for feature $ftr_idx2
  my $ftr_strand;  # strand for feature $ftr_idx
  my $ftr_strand2; # strand for feature $ftr_idx2
  my $found_parent; # flag for if we found a parent or not
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"parent_idx"} = -1; # initialize
    if($ftr_info_AHR->[$ftr_idx]{"type"} eq "mat_peptide") { 
      $ftr_5p_pos = dng_Feature5pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
      $ftr_3p_pos = dng_Feature3pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
      $ftr_strand = dng_FeatureSummaryStrand($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
      for($ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
        $ftr_5p_pos2 = dng_Feature5pMostPosition($ftr_info_AHR->[$ftr_idx2]{"coords"}, $FH_HR);
        $ftr_3p_pos2 = dng_Feature3pMostPosition($ftr_info_AHR->[$ftr_idx2]{"coords"}, $FH_HR);
        $ftr_strand2 = dng_FeatureSummaryStrand($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
        $found_parent = 0;
        if(($ftr_idx != $ftr_idx2) && 
           ($ftr_info_AHR->[$ftr_idx2]{"type"} eq "CDS") && 
           ($ftr_strand eq $ftr_strand2)) { 
          if(($ftr_strand eq "+") && 
             ($ftr_5p_pos2 <= $ftr_5p_pos) &&
             ($ftr_3p_pos2 >= $ftr_3p_pos)) { 
            $found_parent = 1;
          }
          if(($ftr_strand eq "-") && 
             ($ftr_5p_pos2 >= $ftr_5p_pos) &&
             ($ftr_3p_pos2 <= $ftr_3p_pos)) { 
            $found_parent = 1;
          }
          if($found_parent) { 
            if($ftr_info_AHR->[$ftr_idx]{"parent_idx"} != -1) { 
              printf("ftr_5p_pos:  $ftr_5p_pos\n");
              printf("ftr_3p_pos:  $ftr_3p_pos\n");
              printf("ftr_5p_pos2: $ftr_5p_pos2\n");
              printf("ftr_3p_pos2: $ftr_3p_pos2\n");
              ofile_FAIL(sprintf("ERROR in $sub_name, unable to determine parent of mature peptide with coords %s (multiple CDS cover it with coords %s and %s)\n", 
                                  $ftr_info_AHR->[$ftr_idx]{"coords"},
                                  $ftr_info_AHR->[($ftr_info_AHR->[$ftr_idx]{"parent_idx"})]{"coords"}, 
                                  $ftr_info_AHR->[$ftr_idx2]{"coords"}), "dnaorg", 1, $FH_HR);
            }
            $ftr_info_AHR->[$ftr_idx]{"parent_idx"} = $ftr_idx2;
          }
        }
      }
    }
  }   
  return 0;
}
      
#################################################################
# Subroutine: dng_FeatureInfoImputeOutname()
# Incept:     EPN, Mon Apr  1 06:49:20 2019
#
# Purpose:    Fill "outname" values in @{$ftr_info_AHR}
#             This is defined as:
#                  $ftr_info_AHR->[$ftr_idx]{"product"} if defined,
#             else $ftr_info_AHR->[$ftr_idx]{"gene"} if defined,
#             else string of type and type index (e.g. CDS.1)
#
# Arguments: 
#   $ftr_info_AHR:  REF to array of hashes of feature info
#
# Returns:    Feature name string
#
# Dies: Never, nothing is validated
# 
#################################################################
sub dng_FeatureInfoImputeOutname { 
  my $sub_name  = "dng_FeatureInfoImputeOutname";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_AHR, $ftr_idx) = (@_);

  my $nftr = scalar(@{$ftr_info_AHR}); 
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(defined $ftr_info_AHR->[$ftr_idx]{"product"}) { 
      $ftr_info_AHR->[$ftr_idx]{"outname"} = $ftr_info_AHR->[$ftr_idx]{"product"}; 
    }
    elsif(defined $ftr_info_AHR->[$ftr_idx]{"gene"}) { 
      $ftr_info_AHR->[$ftr_idx]{"outname"} = $ftr_info_AHR->[$ftr_idx]{"gene"}; 
    }
    else { 
      $ftr_info_AHR->[$ftr_idx]{"outname"} = dng_FeatureTypeAndTypeIndexString($ftr_info_AHR, $ftr_idx, ".");
    }
  }

  return;
}

#################################################################
# Subroutine: dng_FeatureInfoImpute3paFtrIdx
# Incept:     EPN, Wed Mar 13 13:39:34 2019
# 
# Purpose:    Fill "3pa_ftr_idx" values in @{$ftr_info_AHR}
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub dng_FeatureInfoImpute3paFtrIdx {
  my $sub_name = "dng_FeatureInfoImpute3paFtrIdx";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  # go through all features and determine adjacent mat_peptides (set '3pa_ftr_idx')
  # 
  # $ftr_info_AHR->{"3pa_ftr_idx"}[$ftr_idx] set to $ftr_idx2 if:
  # - $ftr_idx type is mat_peptide
  # - $ftr_idx2 type is mat_peptide
  # - $ftr_idx starts at 1 position 3' of stop position of $ftr_idx
  # - $ftr_idx and $ftr_idx2 are both "+" or both "-" strands
  # 
  # else "-1" if no $ftr_idx2 exists for $ftr_idx that satisfies above
  #
  # dies if more than one $ftr_idx2 satisfies above
  my ($ftr_idx, $ftr_idx2); # feature indices
  my $ftr_3p_pos;  # 3'-most position for feature $ftr_idx
  my $ftr_5p_pos2; # 5'-most position for feature $ftr_idx2
  my $ftr_strand;  # strand for feature $ftr_idx
  my $ftr_strand2; # strand for feature $ftr_idx2
  my $found_adj = 0;
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"3pa_ftr_idx"} = -1;
    if($ftr_info_AHR->[$ftr_idx]{"type"} eq "mat_peptide") { 
      $ftr_3p_pos = dng_Feature3pMostPosition($ftr_info_AHR, $ftr_idx, $FH_HR);
      $ftr_strand = dng_FeatureSummaryStrand($ftr_info_AHR, $ftr_idx, $FH_HR);
      for($ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
        $ftr_5p_pos2 = dng_Feature5pMostPosition($ftr_info_AHR, $ftr_idx2, $FH_HR);
        $ftr_strand2 = dng_FeatureSummaryStrand($ftr_info_AHR, $ftr_idx2, $FH_HR);
        $found_adj = 0;
        if(($ftr_idx != $ftr_idx2) && 
           ($ftr_info_AHR->[$ftr_idx2]{"type"} eq "mat_peptide") &&
           ($ftr_strand eq $ftr_strand2)) { 
          if(($ftr_strand eq "+") && (($ftr_3p_pos+1) == ($ftr_5p_pos2))) { 
            $found_adj = 1;
          }
          if(($ftr_strand eq "-") && (($ftr_3p_pos-1) == ($ftr_5p_pos2))) { 
            $found_adj = 1; 
          }
          if($found_adj) { 
            if($ftr_info_AHR->[$ftr_idx]{"3pa_ftr_idx"} != -1) { 
              ofile_FAIL(sprintf("ERROR in $sub_name, unable to determine 3' mature peptide of mature peptide with coords (multiple mature peptides satisfy criteria)\n", 
                                  $ftr_info_AHR->[$ftr_idx]{"coords"}), "dnaorg", 1, $FH_HR);
            }
            $ftr_info_AHR->[$ftr_idx]{"3pa_ftr_idx"} = $ftr_idx2; 
          }
        }
      }
    }
  }

  return;
}

#################################################################
# Subroutine: dng_FeatureInfoStartStopStrandArrays()
# Incept:     EPN, Fri Mar 15 15:39:35 2019
#
# Synopsis: For all features in a @{$ftr_info_AHR}, validate 
#           "coords" values and fill @{$start_AAR}, @{$stop_AAR} and
#           @{$strand_AAR} based on them.
# 
# Arguments:
#  $ftr_info_AHR:  REF to the feature info array of hashes
#  $start_AAR:     REF to array of start position array to fill here, FILLED here, can be undef
#  $stop_AAR:      REF to array of stop position array to fill here, FILLED here, can be undef
#  $strand_AAR:    REF to array of strand array to fill here with "+" or "-", FILLED here, can be undef
#  $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $coords_str
#
#################################################################
sub dng_FeatureInfoStartStopStrandArrays {
  my $sub_name = "dng_FeatureInfoStartStopStrandArrays";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($ftr_info_AHR, $start_AAR, $stop_AAR, $strand_AAR, $FH_HR) = @_;

  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  my @start_AA  = ();
  my @stop_AA   = ();
  my @strand_AA = ();
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    @{$start_AA[$ftr_idx]}  = ();
    @{$stop_AA[$ftr_idx]}   = ();
    @{$strand_AA[$ftr_idx]} = ();
    dng_FeatureStartStopStrandArrays($ftr_info_AHR->[$ftr_idx]{"coords"}, \@{$start_AA[$ftr_idx]}, \@{$stop_AA[$ftr_idx]}, \@{$strand_AA[$ftr_idx]}, $FH_HR);
  }
  if(defined $start_AAR)  { @{$start_AAR}   = @start_AA;   }
  if(defined $stop_AAR)   { @{$stop_AAR}    = @stop_AA;    }
  if(defined $strand_AAR) { @{$strand_AAR}  = @strand_AA;  }

  return;
}

#################################################################
# Subroutine: dng_FeatureInfoCountType
# Incept:     EPN, Thu Mar 14 12:16:26 2019
# 
# Purpose:    Count number of elements in @{$ftr_info_AHR} 
#             have type of $type.
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $type:          type we are interested in
#
# Returns:    void
# 
# Dies:       never, nothing is validated
#
#################################################################
sub dng_FeatureInfoCountType { 
  my $sub_name = "dng_FeatureInfoCountType";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $type) = @_;

  my $ntype = 0;
  my $nftr = scalar(@{$ftr_info_AHR});

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if((defined $ftr_info_AHR->[$ftr_idx]{"type"}) && 
       ($ftr_info_AHR->[$ftr_idx]{"type"} eq $type)) { 
      $ntype++;
    }
  }
  
  return $ntype;
}

#################################################################
# Subroutine: dng_FeatureInfoValidateCoords
# Incept:     EPN, Fri Mar 15 14:31:36 2019
# 
# Purpose:    Validate that "coords" values are in the proper
#             format and all less than or equal to $length.
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $length:        type we are interested in
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       If a "coords" value is in the incorrect format or
#             if a position in a "coords" value exceeds $length
#
#################################################################
sub dng_FeatureInfoValidateCoords { 
  my $sub_name = "dng_FeatureInfoValidateCoords";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $length, $FH_HR) = @_;

  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);
  my $fail_str = ""; # added to if any elements are out of range

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my @start_A  = (); # array of starts, one per segment
    my @stop_A   = (); # array of stops, one per segment
    # this sub dng_Will die if $ftr_info_AHR->[$ftr_idx]{"coords"} is in incorrect format
    dng_FeatureStartStopStrandArrays($ftr_info_AHR->[$ftr_idx]{"coords"}, \@start_A, \@stop_A, undef, $FH_HR); 
    foreach my $start (@start_A) { if($start > $length) { $fail_str .= "ftr_idx: $ftr_idx, start position $start > $length\n"; } }
    foreach my $stop  (@stop_A)  { if($stop  > $length) { $fail_str .= "ftr_idx: $ftr_idx, stop  position $stop  > $length\n"; } }
  }

  if($fail_str ne "") { 
    ofile_FAIL("ERROR in $sub_name, some coordinates exceed model length ($length):\n$fail_str\n", "dnaorg", 1, $FH_HR);
  }
  
  return;
}

#################################################################
# Subroutine:  dng_FeatureInfoChildrenArrayOfArrays()
# Incept:      EPN, Sun Mar 10 06:22:49 2019
#
# Purpose:     Fill @{$AAR} with arrays of children (feature indices)
#              for each feature in %{$ftr_info_AHR}.
# 
# Arguments: 
#   $ftr_info_AHR:   REF to hash of arrays with information on the features, PRE-FILLED
#   $AAR:            REF to array of arrays of children feature indices, FILLED HERE
#   $FH_HR:          REF to hash of file handles
# 
# Returns:     Nothing.
# 
#
################################################################# 
sub dng_FeatureInfoChildrenArrayOfArrays { 
  my $nargs_expected = 3;
  my $sub_name = "dng_FeatureInfoChildrenArrayOfArrays";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $AAR, $FH_HR) = @_;

  @{$AAR} = ();
  my $nftr = scalar(@{$ftr_info_AHR});

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    dng_FeatureChildrenArray($ftr_info_AHR, $ftr_idx, $nftr, \@{$AAR->[$ftr_idx]}, $FH_HR);
  }
  
  return;
}

#################################################################
# Subroutine: dng_SegmentInfoPopulate()
# Incept:     EPN, Wed Mar 13 13:55:56 2019
#
# Synopsis: Fill @{$sgm_info_AHR} and add to @{$ftr_info_AHR}
#           based on existing information in @{$ftr_info_AHR}.
#
#           The following values are added to %{$sgm_info_HAR}:
#                "start":    start position of segment in the reference genome
#                "stop":     stop position of segment in the reference genome
#                "strand":   strand of segment in the reference genome
#                "map_ftr":  the feature index (array index in ftr_info_AHR) 
#                            this segment is for
#                "is_5p":    '1' if this segment is the 5'-most segment for its feature
#                            (when the segments are joined to make the feature, not 
#                            necessarily in reference genome)
#                "is_3p":    '1' if this segment is the 3'-most model for its feature
#                            (when the segments are joined to make the feature, not 
#                            necessarily in reference genome)
#
#           The following values are added to %{$ftr_info_AHR}:
#                "5p_sgm_idx":   index (in arrays of %sgm_info_HA) of 5'-most segment for this feature
#                "3p_sgm_idx":   index (in arrays of %sgm_info_HA) of 3'-most segment for this feature
# Arguments:
#  $sgm_info_AHR:      ref to array of hashes with information on the segments, FILLED HERE
#  $ftr_info_AHR:      ref to array of hashes with information on the features, ADDED TO HERE
#  $FH_HR:             REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if @{$ftr_info_AHR} is not valid upon entry
#
#################################################################
sub dng_SegmentInfoPopulate {
  my $sub_name = "dng_SegmentInfoPopulate";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($sgm_info_AHR, $ftr_info_AHR, $FH_HR) = @_;

  # ftr_info_AHR should already have array data for keys "type", "coords", "source_idx"
  my @keys_A = ("type", "coords", "source_idx");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  # initialize new %{$ftr_info_AHR} values
  my ($ftr_idx, $ftr_idx2, $sgm_idx, $sgm_idx2); # feature and segment indices
  my ($sgm_start, $sgm_stop, $sgm_strand); # start, stop and strand for a segment
  my $nseg = 0; 
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} = -1; # remains -1 if $ftr_idxs_dup
    $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"} = -2; # remains -2 if $ftr_idxs_dup
    my $ftr_type   = $ftr_info_AHR->[$ftr_idx]{"type"};
    my $ftr_is_dup = ($ftr_info_AHR->[$ftr_idx]{"source_idx"} == $ftr_idx) ? 0 : 1;

    if(! $ftr_is_dup) { 
      # determine start and stop positions of all segments
      my @sgm_start_A  = (); # array of starts, one per segment
      my @sgm_stop_A   = (); # array of stops, one per segment
      my @sgm_strand_A = (); # array of strands ("+", "-"), one per segment
      dng_FeatureStartStopStrandArrays($ftr_info_AHR->[$ftr_idx]{"coords"}, \@sgm_start_A, \@sgm_stop_A, \@sgm_strand_A, $FH_HR);
      my $cur_nseg = scalar(@sgm_start_A);
      for(my $s = 0; $s < $cur_nseg; $s++) { 
        $sgm_info_AHR->[$nseg]{"start"}   = $sgm_start_A[$s];
        $sgm_info_AHR->[$nseg]{"stop"}    = $sgm_stop_A[$s];
        $sgm_info_AHR->[$nseg]{"strand"}  = $sgm_strand_A[$s];
        $sgm_info_AHR->[$nseg]{"map_ftr"} = $ftr_idx;
        if($s == 0) { 
          $sgm_info_AHR->[$nseg]{"is_5p"} = 1;
          $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} = $nseg; 
        }
        else { 
          $sgm_info_AHR->[$nseg]{"is_5p"} = 0;
        }
        if($s == ($cur_nseg-1)) { 
          $sgm_info_AHR->[$nseg]{"is_3p"} = 1;
          $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"} = $nseg;
        }
        else { 
          $sgm_info_AHR->[$nseg]{"is_3p"} = 0;
        }
        $nseg++;
      }
    }
  }

  return;
}

#################################################################
# Subroutine:  dng_FeatureTypeAndTypeIndexString()
# Incept:      EPN, Sun Mar 10 06:41:53 2019
#
# Purpose:     Return a string giving feature type and type index
#              (e.g. "CDS#2") for feature $ftr_idx.
# 
# Arguments: 
#   $ftr_info_AHR:   REF to hash of arrays with information on the features, PRE-FILLED
#   $ftr_idx:        index we are interested in
#   $sep_char:       character to separate type and index in return string
#                    e.g. "#" or "."
# 
# Returns:     String
#
################################################################# 
sub dng_FeatureTypeAndTypeIndexString { 
  my $nargs_expected = 3;
  my $sub_name = "dng_FeatureChildrenArray";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $ftr_idx, $sep_char) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  my $type = $ftr_info_AHR->[$ftr_idx]{"type"};
  my $type_idx = 1;
  for(my $ftr_idx2 = 0; $ftr_idx2 < $ftr_idx; $ftr_idx2++) { 
    if($ftr_info_AHR->[$ftr_idx2]{"type"} eq $type) { $type_idx++; }
  }
  
  return $type . $sep_char . $type_idx;
}

#################################################################
# Subroutine: dng_FeatureTypeIsCds()
# Incept:     EPN, Mon Mar 25 11:07:05 2019
#
# Purpose:    Is feature $ftr_idx a CDS?
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $ftr_idx:        feature index
#
# Returns:    1 or 0
#
# Dies:       never; does not validate anything.
#
################################################################# 
sub dng_FeatureTypeIsCds { 
  my $sub_name = "dng_FeatureTypeIsCds";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return ($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") ? 1 : 0;
}

#################################################################
# Subroutine: dng_FeatureTypeIsGene()
# Incept:     
#
# Purpose:    Is feature $ftr_idx a gene?
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $ftr_idx:        feature index
#
# Returns:    1 or 0
#
# Dies:       never; does not validate anything.
#
################################################################# 
sub dng_FeatureTypeIsGene { 
  my $sub_name = "dng_FeatureTypeIsGene";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return ($ftr_info_AHR->[$ftr_idx]{"type"} eq "gene") ? 1 : 0;
}


#################################################################
# Subroutine: dng_FeatureTypeIsCdsOrMatPeptide()
# Incept:     EPN, Mon Feb 25 14:30:34 2019
#
# Purpose:    Is feature $ftr_idx a CDS or mature peptide?
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $ftr_idx:        feature index
#
# Returns:    1 or 0
#
# Dies:       never; does not validate anything.
#
################################################################# 
sub dng_FeatureTypeIsCdsOrMatPeptide { 
  my $sub_name = "dng_FeatureTypeIsCdsOrMatPeptide";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return (($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") || 
          ($ftr_info_AHR->[$ftr_idx]{"type"} eq "mat_peptide")) ? 1 : 0;
}

#################################################################
# Subroutine: dng_FeatureTypeIsCdsOrMatPeptideOrGene()
# Incept:     EPN, Wed Apr  3 14:31:33 2019
#
# Purpose:    Is feature $ftr_idx a CDS or mature peptide or gene?
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $ftr_idx:        feature index
#
# Returns:    1 or 0
#
# Dies:       never; does not validate anything.
#
################################################################# 
sub dng_FeatureTypeIsCdsOrMatPeptideOrGene { 
  my $sub_name = "dng_FeatureTypeIsCdsOrMatPeptideOrGene";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return (($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") || 
          ($ftr_info_AHR->[$ftr_idx]{"type"} eq "mat_peptide") ||
          ($ftr_info_AHR->[$ftr_idx]{"type"} eq "gene")) ? 1 : 0;
}

#################################################################
# Subroutine:  dng_FeatureChildrenArray()
# Incept:      EPN, Sun Mar 10 06:26:03 2019
#
# Purpose:     Fill @{$AAR} with array of children (feature indices)
#              for feature $ftr_idx in %{$ftr_info_AHR}.
# 
# Arguments: 
#   $ftr_info_AHR:   REF to hash of arrays with information on the features, PRE-FILLED
#   $ftr_idx:        index we are interested in
#   $nftr:           number of features in %{$ftr_info_AHR}
#   $AR:             REF to array children feature indices, FILLED HERE
#   $FH_HR:          REF to hash of file handles
# 
# Returns:     Nothing.
# 
#
################################################################# 
sub dng_FeatureChildrenArray { 
  my $nargs_expected = 5;
  my $sub_name = "dng_FeatureChildrenArray";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $ftr_idx, $nftr, $AR, $FH_HR) = @_;

  @{$AR} = ();

  for(my $ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
    if(($ftr_idx2 != $ftr_idx) && 
       ($ftr_info_AHR->[$ftr_idx2]{"parent_idx"} == $ftr_idx)) { 
      push(@{$AR}, $ftr_idx2);
    }
  }
  
  return;
}

#################################################################
# Subroutine: dng_FeatureNumSegments()
# Incept:     EPN, Tue Mar  5 13:05:38 2019
#
# Purpose:    Return number of segments in feature $ftr_idx.
#
# Arguments: 
#   $ftr_info_AHR:  REF to array of hashes of feature info
#   $ftr_idx:       feature index we are interested in
#
# Returns:    Number of segments for $ftr_idx.
#
# Dies: Never, nothing is validated
# 
#################################################################
sub dng_FeatureNumSegments { 
  my $sub_name  = "featureNumSegments";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_AHR, $ftr_idx) = (@_);

  return ($ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"} - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} + 1);
}



#################################################################
# Subroutine: dng_FeatureIsDuplicate()
# Incept:      EPN, Sun Mar 10 07:04:24 2019
#
# Purpose:    Is feature $ftr_idx a duplicate of another feature?
#             This is true if $ftr_info_AHR->[$ftr_idx]{"source_idx"} != $ftr_idx
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $ftr_idx:        feature index
#
# Returns:    1 or 0 
# 
# Dies:       never; does not validate anything.
#
################################################################# 
sub dng_FeatureIsDuplicate { 
  my $sub_name = "dng_FeatureIsDuplicate";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return(($ftr_info_AHR->[$ftr_idx]{"source_idx"} != $ftr_idx) ? 1 : 0);
}

#################################################################
# Subroutine: dng_Feature5pMostPosition()
# Incept:      EPN, Fri Mar  8 12:57:21 2019
#
# Purpose:    Return 5'-most position in all segments for a feature.
#
# Arguments: 
#  $coords:  coords value from feature info
#  $FH_HR:   ref to hash of file handles
# 
# Returns:   5'-most position
#
# Dies:      if $coords is not parseable.
#
################################################################# 
sub dng_Feature5pMostPosition { 
  my $sub_name = "dng_Feature5pMostPosition";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords, $FH_HR) = @_;
  
  if($coords =~ /^(\d+)\.\.\d+/) { 
    return $1;
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse ftr_info_HA coords string " . $coords, "dnaorg", 1, $FH_HR); 
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine: dng_Feature3pMostPosition()
# Incept:      EPN, Fri Mar  8 13:00:31 2019
#
# Purpose:    Return 3'-most position in all segments for a feature.
#
# Arguments: 
#  $coords:  coords value from feature info
#  $FH_HR:   ref to hash of file handles
# 
# Returns:   3'-most position
#
# Dies:      if $coords is not parseable.
#
################################################################# 
sub dng_Feature3pMostPosition { 
  my $sub_name = "dng_Feature3pMostPosition";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords, $FH_HR) = @_;
  
  if($coords =~ /\d+\.\.(\d+)\:[\+\-]$/) { 
    return $1;
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse ftr_info_HA coords string " . $coords, "dnaorg", 1, $FH_HR); 
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine: dng_FeatureSummarizeSegment()
# Incept:      EPN, Fri Mar  1 12:36:36 2019
#
# Purpose:    Return a string indicating what model this is
#             for features that are covered by multiple model spans.
#
# Arguments: 
#  $ftr_info_AHR: ref to feature info array of hashes, PRE-FILLED
#  $sgm_info_AHR: ref to segment info array of hashes, PRE-FILLED
#  $sgm_idx:      model index
#
# Returns:    "" if this is the only model for this feature
#             string like ", model 1 of 2", if not
# 
# Dies:       never; does not validate anything.
#
################################################################# 
sub dng_FeatureSummarizeSegment { 
  my $sub_name = "dng_FeatureSummarizeSegment";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $sgm_info_AHR, $sgm_idx) = @_;

  my $ftr_idx = $sgm_info_AHR->[$sgm_idx]{"map_ftr"};
  my $nmdl = ($ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"} - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}) + 1;
  if($nmdl > 1) { 
    return sprintf(", segment %d of %d", ($sgm_idx - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}) + 1, $nmdl);
  }

  return ""; # return "" if $nmdl == 1;
}



#################################################################
# Subroutine: dng_FeatureStartStopStrandArrays()
# Incept:     EPN, Sat Mar  9 05:50:10 2019
#
# Synopsis: Given a comma separated coords string, parse it, 
#           validate it, and fill @{$start_AR}, @{$stop_AR} and
#           @{$strand_AR} based on it.
# 
# Arguments:
#  $coords:       coordinate string
#  $start_AR:     REF to start position array to fill here, FILLED here, can be undef
#  $stop_AR:      REF to stop position array to fill here, FILLED here, can be undef
#  $strand_AR:    REF to strand array to fill here with "+" or "-", FILLED here, can be undef
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $coords
#
#################################################################
sub dng_FeatureStartStopStrandArrays {
  my $sub_name = "dng_FeatureStartStopStrandArrays";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $start_AR, $stop_AR, $strand_AR, $FH_HR) = @_;
  if(! defined $coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", "dnaorg", 1, $FH_HR); 
  }

  my @start_A  = ();
  my @stop_A   = ();
  my @strand_A = ();
  my ($start, $stop, $strand, $sgm_idx);
  my @coords_A  = split(",", $coords);
  my $nsgm = scalar(@coords_A);
  for($sgm_idx = 0; $sgm_idx < $nsgm; $sgm_idx++) { 
    ($start, $stop, $strand) = dng_CoordsTokenParse($coords_A[$sgm_idx], $FH_HR);
    # dng_CoordsTokenParse() will fail if unable to parse $coords_A[$sgm_idx]
    push(@start_A,  $start);
    push(@stop_A,   $stop);
    push(@strand_A, $strand); 
  }

  if(defined $start_AR)  { @{$start_AR}   = @start_A;  }
  if(defined $stop_AR)   { @{$stop_AR}    = @stop_A;   }
  if(defined $strand_AR) { @{$strand_AR}  = @strand_A;  }

  return;
}

#################################################################
# Subroutine: dng_FeatureSummaryStrand
# Incept:     EPN, Wed Mar 13 15:38:06 2019
# 
# Purpose:    Summarize the strandedness of segments for a feature
#             by parsing the "coords" value.
# 
# Arguments:
#   $coords:   coords string to complement
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    "+" if all segments in $coords are "+"
#             "-" if all segments in $coords are "-"
#             "!" if >= 1 segment in $coords is "+"
#                 and >= 1 segment in $coords is "-"
#
# Dies:      if unable to parse $coords
#
#################################################################
sub dng_FeatureSummaryStrand { 
  my $sub_name = "dng_FeatureSummaryStrand";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($coords, $FH_HR) = @_;

  # Examples we can parse: 
  # $coords                  return value
  # -----------------------  -----------------
  # 1-200:+                  200-1:-
  # 1-200:+,300-400:+        400-300:-,200-1:-

  my @strand_A = ();
  dng_FeatureStartStopStrandArrays($coords, undef, undef, \@strand_A, $FH_HR);

  my $npos = 0;
  my $nneg = 0;
  foreach my $strand (@strand_A) { 
    if   ($strand eq "+") { $npos++; }
    elsif($strand eq "-") { $nneg++; }
    else { ofile_FAIL("ERROR in $sub_name, unable to determine strands in coords $coords", "dnaorg", 1, $FH_HR); }
  }

  if(($npos >  0) && ($nneg == 0)) { return "+"; }
  if(($npos == 0) && ($nneg >  0)) { return "-"; }
  if(($npos == 0) && ($nneg == 0)) { 
    ofile_FAIL("ERROR in $sub_name, unable to determine strands in coords $coords", "dnaorg", 1, $FH_HR); 
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine: dng_FeaturePositionSpecificValueBreakdown()
# Incept:     EPN, Tue Apr  2 10:22:16 2019
#
# Purpose:    Breakdown a list of position specific values
#             from a string in %{$ftr_info_AHR->[$ftr_idx]}
#             and fill %HR with key/value pairs.
# 
#             String must be in format of one or more tokens
#             of: "<d>:<s>" separated by ";" if more than one.
#
#             If $ftr_info_AHR->[$ftr_idx] does not exist just
#             return.
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $ftr_idx:        feature index
#  $key:            key in $ftr_info_AHR->[$ftr_idx]
#  $HR:             ref to hash to fill
#  $FH_HR:          ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $ftr_info_AHR->[$ftr_idx] exists but cannot
#             be parsed.
#
################################################################# 
sub dng_FeaturePositionSpecificValueBreakdown { 
  my $sub_name = "dng_FeaturePositionSpecificValueBreakdown";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx, $key, $HR, $FH_HR) = @_;
 
  if(defined $ftr_info_AHR->[$ftr_idx]{$key}) { 
    my @tok_A = split(";", $ftr_info_AHR->[$ftr_idx]{$key});
    foreach my $tok (@tok_A) { 
      if($tok =~ /^(\d+)\:(\S+)$/) { 
        $HR->{$1} = $2;
      }
      else { 
        ofile_FAIL("ERROR, in $sub_name, unable to parse token $tok parsed out of " . $ftr_info_AHR->[$ftr_idx]{$key}, 1, $FH_HR);
      }
    }
  }

  return;
}

#################################################################
# Subroutine: dng_AlertInfoInitialize()
# Incept:     EPN, Fri Mar  4 12:56:43 2016
#
# Purpose:    Set the initial values in an alert info hash or hashes,
#             using the hardcoded information in this
#             function.
#
# Arguments:
#   $alt_info_HHR:  REF to hash of arrays of alert information, FILLED HERE
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $alt_info_HHR already has keys upon entering this function
#
#################################################################
sub dng_AlertInfoInitialize { 
  my $sub_name = "dng_AlertInfoInitialize";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $FH_HR) = (@_);

  if(scalar (keys(%{$alt_info_HHR})) > 0) { 
    ofile_FAIL("ERROR in $sub_name, alert info hash of arrays already has at least one key", "dnaorg", 1, $FH_HR);
  }

  # add each alert code, this function will die if we try to add the same code twice, or if something is wrong 
  # with how we try to add it (args to dng_AlertInfoAdd don't pass the contract check)

  # classification errors
  dng_AlertInfoAdd($alt_info_HHR, "c_noa", "sequence",
                   "No Annotation", # short description
                   "no significant similarity detected", # long  description
                   1, 1, 1, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "c_mst", "sequence",
                   "Minus Strand", # short description
                   "sequence appears to be reverse complemented", # long description
                   1, 1, 1, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "c_usg", "sequence",
                   "Unexpected Subgroup Classification", # short description
                   "score difference too large between best overall model and best expected subgroup model", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "c_ugr", "sequence",
                   "Unexpected Group Classification", # short description
                   "score difference too large between best overall model and best expected group model", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "c_loc", "sequence",
                   "Low Coverage", # short description, 
                   "low sequence fraction with significant similarity to homology model", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "c_lod", "sequence",
                   "Low Score Difference", # short description
                   "low score difference between best overall model and second best model (not in best model's subgroup)", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "c_vld", "sequence",
                   "Very Low Score Difference", # description
                   "very low score difference between best overall model and second best model (not in best model's subgroup)", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "c_los", "sequence",
                   "Low Score", # short description
                   "score to homology model below low threshold", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "c_vls", "sequence",
                   "Very Low Score", # short description
                   "score to homology model below very low threshold", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "c_hbi", "sequence",
                   "Biased Sequence", # short description
                   "high fraction of score attributed to biased sequence composition", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "n_div", "sequence",
                   "Unexpected Divergence", # short description
                   "sequence is too divergent to confidently assign nucleotide-based annotation", # long description
                   1, 1, 1, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "b_zft", "sequence",
                   "No Features Annotated", # short description
                   "sequence similarity to homology model does not overlap with any features", # long description
                   1, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HHR, "n_str", "feature",
                   "Mutation at Start", # short description
                   "expected start codon could not be identified", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "n_stp", "feature",
                   "Mutation at End", # short description
                   "expected stop codon could not be identified, predicted CDS stop by homology is invalid", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "n_nst", "feature",  
                   "Mutation at End", # short description
                   "expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted valid start codon", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "n_ext", "feature",
                   "Mutation at End", # short description
                   "expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "n_nm3", "feature",  
                   "Unexpected Length", # short description
                   "length of complete coding (CDS or mat_peptide) feature is not a multiple of 3", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "n_trc", "feature",
                   "Unexpected Stop Codon", # short description
                   "in-frame stop codon exists 5' of stop position predicted by homology to reference", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "p_trc", "feature",
                   "Unexpected Stop Codon", # short description
                   "stop codon in protein-based alignment", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "b_per", "feature",
                   "Peptide Translation Problem", # short description
                   "mat_peptide may not be translated because its CDS has a problem", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "b_non", "feature",
                   "Indefinite Annotation", # short description
                   "protein-based search identifies CDS not identified in nucleotide-based search", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "b_nop", "feature",
                   "Indefinite Annotation", # short description
                   "nucleotide-based search identifies CDS not identified in protein-based search", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "n_gp5", "feature",
                   "Indefinite Annotation at Start", # short description
                   "alignment to homology model is a gap at 5' boundary", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "n_lp5", "feature",
                   "Indefinite Annotation at Start", # short description
                   "alignment to homology model has low confidence at 5' boundary", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "b_p5l", "feature",
                   "Indefinite Annotation at Start", # short description
                   "protein-based alignment extends past nucleotide-based alignment at 5' end", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "b_p5s", "feature",
                   "Indefinite Annotation at Start", , # short description
                   "protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "n_gp3", "feature",
                   "Indefinite Annotation at End", # short description
                   "alignment to homology model is a gap at 3' boundary", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "n_lp3", "feature",
                   "Indefinite Annotation at End", # short description
                   "alignment to homology model has low confidence at 3' boundary", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "b_p3l", "feature",
                   "Indefinite Annotation at End", # short description
                   "protein-based alignment extends past nucleotide-based alignment at 3' end", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "b_p3s", "feature",
                   "Indefinite Annotation at End", # short description
                   "protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "b_cst", "feature",
                   "Indefinite Strand", # short description
                   "strand mismatch between protein-based and nucleotide-based predictions", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "p_lin", "feature",
                   "Insertion of Nucleotides", # short description
                   "too large of an insertion in protein-based alignment", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HHR, "p_lde", "feature",
                   "Deletion of Nucleotides", # short description
                   "too large of a deletion in protein-based alignment", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  # define the ftbl_invalid_by values, these are one-sided, any error code listed in the 
  # 3rd argument invalidates the 2nd argument error code, but not vice versa

  # n_trc, n_ext and n_nst are preferred to n_stp
  dng_AlertInfoSetFTableInvalidatedBy($alt_info_HHR, "n_stp", "n_trc,n_ext,n_nst", $FH_HR); 

  # n_div is preferred to b_zft
  dng_AlertInfoSetFTableInvalidatedBy($alt_info_HHR, "b_zft", "n_div", $FH_HR);

  # validate the alert info hash
  #dng_AlertInfoValidate($alt_info_HHR, undef, $FH_HR); 

  return;
}

#################################################################
# Subroutine: dng_AlertInfoAdd()
# Incept:     EPN, Fri Mar  4 13:09:52 2016
#
# Purpose:    Add an element to the error info hash.
#             Die if the same error code already exists.
#
# Arguments:
#   $alt_info_HHR:    REF to array of hashes of error information, FILLED HERE
#   $code:            the code of the element we are adding
#   $pertype:         the 'per-type' of the element we are adding, "sequence" or "feature"
#   $sdesc:           short description of the alert we are adding
#   $ldesc:           long  description of the alert we are adding
#   $always_fails:    '1' if this alert *always* causes its sequence to FAIL, '0' if not
#   $causes_failure:  '1' if this alert causes its sequence to FAIL by default, '0' if not
#   $prevents_annot:  '1' if this alert prevents its sequence from being annotated, '0' if not
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $alt_info_HHR->{"$code"} already exists
#          if $type ne "feature and ne "sequence"
#          if $type ne "sequence" and $prevents_annot == 1 (not allowed)
#
######################p###########################################
sub dng_AlertInfoAdd { 
  my $sub_name = "dng_AlertInfoAdd";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $code, $pertype, $sdesc, $ldesc, $always_fails, $causes_failure, $prevents_annot, $FH_HR) = (@_);

  # make sure $pertype is valid
  if(($pertype ne "feature") && ($pertype ne "sequence")) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code with per-type $pertype that is not neither \"feature\" nor \"sequence\".", "dnaorg", 1, $FH_HR); 
  }
  
  # make sure $always_fails is valid
  if((! defined $causes_failure) || (($causes_failure != 0) && ($causes_failure != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but causes_failure is undefined or not 0 or 1", "dnaorg", 1, $FH_HR);
  }

  # make sure $causes_failure is valid, and makes sense with $always_fail
  if((! defined $causes_failure) || (($causes_failure != 0) && ($causes_failure != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but causes_failure is undefined or not 0 or 1", "dnaorg", 1, $FH_HR);
  }
  if($always_fails && (! $causes_failure)) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but always_fails is 1 and causes_failure 0", "dnaorg", 1, $FH_HR);
  }
  
  # make sure $prevents_annot is valid
  if((! defined $prevents_annot) || (($prevents_annot != 0) && ($prevents_annot != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but prevents_annot is undefined or not 0 or 1", "dnaorg", 1, $FH_HR);
  }

  # make sure $prevents_annot is only 1 if $pertype is "sequence"
  if(($prevents_annot == 1) && ($pertype ne "sequence")) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but prevents_annot is 1 and pertype is feature", "dnaorg", 1, $FH_HR);
  }
  
  # check if $code already exists
  if(defined $alt_info_HHR->{$code}) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code, but it already exists in the error info hash", "dnaorg", 1, $FH_HR);
  }

  # find highest "order" index so far in the 2D hash
  my $order = scalar(keys (%{$alt_info_HHR}));
  
  # initialize
  %{$alt_info_HHR->{$code}} = ();

  $alt_info_HHR->{$code}{"order"}           = $order; # first "order" value will be '0', second will be '1', etc.
  $alt_info_HHR->{$code}{"pertype"}         = $pertype;
  $alt_info_HHR->{$code}{"sdesc"}           = $sdesc;
  $alt_info_HHR->{$code}{"ldesc"}           = $ldesc;
  $alt_info_HHR->{$code}{"always_fails"}    = $always_fails;
  $alt_info_HHR->{$code}{"causes_failure"}  = $causes_failure;
  $alt_info_HHR->{$code}{"prevents_annot"}  = $prevents_annot;
  $alt_info_HHR->{$code}{"ftbl_invalid_by"} = ""; # initialized to no invalid_by's, possibly added to later with setFTableInvalidatedByErrorInfoHash()

  return;
}

#################################################################
# Subroutine: dng_AlertInfoSetFTableInvalidatedBy
# Incept:     EPN, Thu Nov  1 10:10:03 2018
#
# Purpose:    Add to the ftbl_invalid_by value for an error code $code1 given
#             a string of other error codes $code2. ftr_invalid_by values 
#             are uni-directional.
#
# Arguments:
#   $alt_info_HHR:  REF to hash of hashes of error information, FILLED HERE
#   $code1:         the code of the element we are adding ftbl_invalid_by values for
#   $code2str:      the codes $code1 is invalidated by, separated by commas
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if one of the error codes in $code1 or $code2str do not
#          exist in %{$alt_info_HHR}.
#
#################################################################
sub dng_AlertInfoSetFTableInvalidatedBy {
  my $sub_name = "dng_AlertInfoSetFTableInvalidatedBy";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $code1, $code2str, $FH_HR) = (@_);

  # verify the codes in $code2str
  my @code2_A = split(',', $code2str);
  foreach my $code2 (@code2_A) { 
    if(! defined $alt_info_HHR->{$code2}) { 
      ofile_FAIL("ERROR in $sub_name, trying to add invalidated by relationship between codes $code1 and $code2, but $code2 does not exist in the error info hash", "dnaorg", 1, $FH_HR);
    }
    if($code1 eq $code2) { 
      ofile_FAIL("ERROR in $sub_name, trying to add invalidated by relationship between a code and itself: $code1 and $code2", "dnaorg", 1, $FH_HR);
    }
  }

  # set the value
  $alt_info_HHR->{$code1}{"ftbl_invalid_by"} = $code2str;

  return;
}

#################################################################
# Subroutine: dng_AlertInfoSetCausesFailure
# Incept:     EPN, Wed Apr  3 12:58:58 2019
#
# Purpose:    Set the "causes_failure" value for %{$ftr_info_HHR->{$code}
#
# Arguments:
#   $alt_info_HHR: REF to hash of hashes of error information, FILLED HERE
#   $code:         the code of the element we are adding ftbl_invalid_by values for
#   $value:        value we are setting $ftbl_info_HHR->{$code}{"causes_failure"} to 
#   $FH_HR:        REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $code does not exist in %{$atl_info_HHR}
#          if $value is not '1' or '0'
#
#################################################################
sub dng_AlertInfoSetCausesFailure {
  my $sub_name = "dng_AlertInfoSetCausesFailure";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $code, $value, $FH_HR) = (@_);

  if(! defined $alt_info_HHR->{$code}) { 
    ofile_FAIL("ERROR in $sub_name, trying to set causes_failure for invalid code $code", "dnaorg", 1, $FH_HR);
  }
  if(($value ne "1") && ($value ne "0")) { 
    ofile_FAIL("ERROR in $sub_name, trying to set causes_failure to invalid value $value (must be 1 or 0)", "dnaorg", 1, $FH_HR);
  }

  $alt_info_HHR->{$code}{"causes_failure"} = $value;

  return;
}

#################################################################
# Subroutine:  dng_AlertInfoDump()
# Incept:      EPN, Wed Apr  3 11:11:26 2019
#
# Purpose:    Output all the information in the alert info 2D hash.
#
# Arguments: 
#  $alt_info_HHR:  REF to the alert info hash of arrays, PRE-FILLED
#  $FH:            file handle to print to
#
# Returns:    void
#
# Dies:       never
#
#################################################################
sub dng_AlertInfoDump { 
  my $sub_name = "dng_AlertInfoDump";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($alt_info_HHR, $FH) = @_;

  my $w_sdesc   = utl_HHMaxLengthValueGiven2DKey($alt_info_HHR, "sdesc");
  my $w_ldesc   = utl_HHMaxLengthValueGiven2DKey($alt_info_HHR, "ldesc");
  my $w_pertype = utl_HHMaxLengthValueGiven2DKey($alt_info_HHR, "pertype");
  my $w_afails  = length("always");
  my $w_fails   = length("fails");
  my $w_annot   = length("prevents");
  my $w_invalid = utl_HHMaxLengthValueGiven2DKey($alt_info_HHR, "ftbl_invalid_by");
  $w_invalid = utl_Max($w_invalid, length("invalidated"));

  # determine order of codes to print
  my @code_A = ();
  foreach my $code (sort keys (%{$alt_info_HHR})) { 
    $code_A[($alt_info_HHR->{$code}{"order"})] = $code;
  }

  printf $FH ("%-4s  %5s  %-*s  %*s  %*s  %*s  %*s  %-*s\n", 
              "#", "", $w_sdesc, "short", $w_afails, "fails", $w_fails, "", $w_annot, "prevents", $w_invalid, "invalidated", $w_ldesc, "long");
  printf $FH ("%-4s  %5s  %-*s  %*s  %*s  %*s  %*s  %-*s\n", 
              "#idx", "code", $w_sdesc, "desc", $w_afails, "always", $w_fails, "fails", $w_annot, "annot", $w_invalid, "by", $w_ldesc, "desc");
  printf $FH ("%-4s  %5s  %-*s  %*s  %*s  %*s  %*s  %-*s\n", 
              "#---", "-----", 
              $w_sdesc,   utl_StringMonoChar($w_sdesc,   "-", undef), 
              $w_afails,  utl_StringMonoChar($w_afails,  "-", undef), 
              $w_fails,   utl_StringMonoChar($w_fails,   "-", undef), 
              $w_annot,   utl_StringMonoChar($w_annot,   "-", undef), 
              $w_invalid, utl_StringMonoChar($w_invalid, "-", undef), 
              $w_ldesc,   utl_StringMonoChar($w_ldesc,   "-", undef));

  my $idx = 0;
  foreach my $code (@code_A) { 
    $idx++;
    printf $FH ("%-4s  %5s  %-*s  %*s  %*s  %*s  %*s  %-*s\n", 
                $idx, $code, 
                $w_sdesc,  helper_tabular_replace_spaces($alt_info_HHR->{$code}{"sdesc"}), 
                $w_afails, $alt_info_HHR->{$code}{"always_fails"}   ? "yes" : "no",
                $w_fails,  $alt_info_HHR->{$code}{"causes_failure"} ? "yes" : "no",
                $w_annot,  $alt_info_HHR->{$code}{"prevents_annot"} ? "yes" : "no",
                $w_invalid, $alt_info_HHR->{$code}{"ftbl_invalid_by"}, 
                $w_ldesc, $alt_info_HHR->{$code}{"ldesc"});
  }

  return;
}

#################################################################
# Subroutine: dng_SubmitJob()
# Incept:      EPN, Wed Feb  6 12:35:04 2019
#
# Purpose:     Submits a job to sge.
#
# Arguments:
#   $cmd:            command to run
#   $job_name:       name for job
#   $alt_file:       name of err file to create, can be "/dev/null"
#   $mem_gb:         number of Gb of memory required
#   $nsecs:          maximum number of seconds to allow jobs to take
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#   $ofile_info_HHR: REF to the 2D hash of output file information, ADDED TO HERE 
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if qsub dng_$cmd fails
#################################################################
sub dng_SubmitJob {
  my $sub_name = "dng_SubmitJob()";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $job_name, $alt_file, $mem_gb, $nsecs, $opt_HHR, $ofile_info_HHR) = @_;
  
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  if(($alt_file ne "/dev/null") && (-e $alt_file)) { 
    utl_FileRemoveUsingSystemRm($alt_file, $sub_name, $opt_HHR, $ofile_info_HHR); 
  }
  my $submit_cmd = sprintf("qsub -N $job_name -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $alt_file -m n -l h_rt=%d,h_vmem=%dG,mem_free=%dG,reserve_mem=%dG,m_mem_free=%dG " . "\"" . $cmd . "\" > /dev/null\n", $nsecs, $mem_gb, $mem_gb, $mem_gb, $mem_gb);

  utl_RunCommand($submit_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine:  dng_WaitForFarmJobsToFinish()
# Incept:      EPN, Mon Feb 29 16:20:54 2016
#              EPN, Wed Aug 31 09:07:05 2016 [moved from dnaorg_annotate.pl to dnaorg.pm]
#
# Purpose: Wait for jobs on the farm to finish by checking the final
#          line of their output files (in @{$outfile_AR}) to see
#          if the final line is exactly the string
#          $finished_string. We'll wait a maximum of $nmin
#          minutes, then return the number of jobs that have
#          finished. If all jobs finish before $nmin minutes we
#          return at that point.
#
#          If($do_cmalign) behavior changes to a more detailed check
#          of the cmalign output file, see cmalignCheckStdOutput().
#          
#          When $do_errcheck is 1, this function considers any output
#          written to stderr output files in @{$errfile_AR} to mean
#          that the corresponding job has 'failed' and should be
#          considered to have finished. When $do_errchecks is 0
#          we don't look at the err files.
# 
#
# Arguments: 
#  $do_cmalign:      '1' if we're running cmalign, which requires special care because we
#                    handle two cases: finish successfully or die with a specific error
#  $out_file_AHR:    ref to array of hashes of output files that will be created by jobs we are waiting for
#  $success_AR:      ref to array of success values, FILLED HERE, can be undef if ! $do_cmalign
#                    these will always all be '1' unless $do_cmalign
#                    if($do_cmalign) some may be '0'
#  $mxsize_AR:       ref to array of required matrix sizes, CAN BE UNDEF
#                    $mxsize_AR->[$j] set to value readh from cmalign output, if $success_AR->[$j] == 0
#                    else set to '0'
#  $finished_str:    string that indicates a job is finished e.g. "[ok]"
#  $opt_HHR:         REF to options hash
#  $FH_HR:           REF to hash of file handles
#
# Returns:     Number of jobs (<= scalar(@outfile_A)) that have
#              finished.
# 
# Dies: never.
#
################################################################# 
sub dng_WaitForFarmJobsToFinish { 
  my $sub_name = "dng_WaitForFarmJobsToFinish()";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($do_cmalign, $out_file_AHR, $success_AR, $mxsize_AR, $finished_str, $opt_HHR, $FH_HR) = @_;

  my $log_FH = $FH_HR->{"log"};
  my $nmin = opt_Get("--wait", $opt_HHR);
  my $do_errcheck = opt_Get("--errcheck", $opt_HHR);

  # contract check
  if(($do_cmalign) && (! exists $out_file_AHR->[0]{"stdout"})) { 
    ofile_FAIL("ERROR in $sub_name, cmalign mode, no stdout files in out_file_AHR", "dnaorg", 1, $FH_HR);
  }
  if((! $do_cmalign) && (! exists $out_file_AHR->[0]{"tblout"})) { 
    ofile_FAIL("ERROR in $sub_name, cmsearch mode, no stdout files in out_file_AHR", "dnaorg", 1, $FH_HR);
  }
  if(! exists $out_file_AHR->[0]{"err"}) { 
    ofile_FAIL("ERROR in $sub_name, no err files in out_file_AHR", "dnaorg", 1, $FH_HR);
  }

  my $outkey = ($do_cmalign) ? "stdout" : "tblout";
  my @outfile_A = ();
  my @errfile_A = ();
  utl_ArrayOfHashesToArray($out_file_AHR, \@outfile_A, $outkey);
  utl_ArrayOfHashesToArray($out_file_AHR, \@errfile_A, "err");

  my $njobs = scalar(@outfile_A);
  if($njobs != scalar(@errfile_A)) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, number of elements in outfile array ($njobs) differ from number of jobs in errfile array (%d)", scalar(@errfile_A)), "dnaorg", 1, $FH_HR);
  }
  my @is_finished_A  = ();  # $is_finished_A[$i] is 1 if job $i is finished (either successfully or having failed), else 0
  my @is_failed_A    = ();  # $is_failed_A[$i] is 1 if job $i has finished and failed (all failed jobs are considered 
                            # to be finished), else 0. We only use this array if the --errcheck option is enabled.
  my $nfinished      = 0;   # number of jobs finished
  my $nfail          = 0;   # number of jobs that have failed
  my $cur_sleep_secs = 15;  # number of seconds to wait between checks, we'll double this until we reach $max_sleep, every $doubling_secs seconds
  my $doubling_secs  = 120; # number of seconds to wait before doublign $cur_sleep
  my $max_sleep_secs = 120; # maximum number of seconds we'll wait between checks
  my $secs_waited    = 0;   # number of total seconds we've waited thus far

  # initialize @is_finished_A to all '0's
  for(my $i = 0; $i < $njobs; $i++) { 
    $is_finished_A[$i] = 0;
    $is_failed_A[$i] = 0;
  }

  my $keep_going = 1;  # set to 0 if all jobs are finished
  ofile_OutputString($log_FH, 1, "\n");
  while(($secs_waited < (($nmin * 60) + $cur_sleep_secs)) && # we add $cur_sleep so we check one final time before exiting after time limit is reached
        ($keep_going)) { 
    # check to see if jobs are finished, every $cur_sleep seconds
    sleep($cur_sleep_secs);
    $secs_waited += $cur_sleep_secs;
    if($secs_waited >= $doubling_secs) { 
      $cur_sleep_secs *= 2;
      if($cur_sleep_secs > $max_sleep_secs) { # reset to max if we've exceeded it
        $cur_sleep_secs = $max_sleep_secs;
      }
    }

    for(my $i = 0; $i < $njobs; $i++) { 
      if(! $is_finished_A[$i]) { 
        if(-s $outfile_A[$i]) { 
          if($do_cmalign) { 
            my $success = dng_CmalignCheckStdOutput($outfile_A[$i], 
                                                    (defined $mxsize_AR) ? \$mxsize_AR->[$i] : undef,
                                                    $FH_HR);
            if($success == 0 || $success == 1) { 
              if(defined $success_AR) { $success_AR->[$i] = $success; }
              $is_finished_A[$i] = 1;
              $nfinished++;
            }
          }
          else { 
            my $final_line = `tail -n 1 $outfile_A[$i]`;
            chomp $final_line;
            if($final_line =~ m/\r$/) { chop $final_line; } # remove ^M if it exists
            if($final_line =~ m/\Q$finished_str\E/) { 
              if(defined $success_AR) { $success_AR->[$i] = 1; } # if we're not running cmalign, if we see $finished_str, job was successful
              $is_finished_A[$i] = 1;
              $nfinished++;
            }
          }
        }
        if(($do_errcheck) && (-s $errfile_A[$i])) { # errfile exists and is non-empty, this is a failure, even if we saw $finished_str above
          if(! $is_finished_A[$i]) { 
            $nfinished++;
          }
          $is_finished_A[$i] = 1;
          $is_failed_A[$i] = 1;
          $nfail++;
        }
      }
    }

    # output update
    ofile_OutputString($log_FH, 1, sprintf("#\t%4d of %4d jobs finished (%.1f minutes spent waiting)\n", $nfinished, $njobs, $secs_waited / 60.));

    if($nfinished == $njobs) { 
      # we're done break out of it
      $keep_going = 0;
    }
  }

  if($nfail > 0) { 
    # construct error message
    my $errmsg = "ERROR in $sub_name, $nfail of $njobs finished in error (output to their respective error files).\n";
    $errmsg .= "Specifically the jobs that were supposed to create the following output and err files:\n";
    for(my $i = 0; $i < $njobs; $i++) { 
      if($is_failed_A[$i]) { 
        $errmsg .= "\t$outfile_A[$i]\t$errfile_A[$i]\n";
      }
    }
    ofile_FAIL($errmsg, "dnaorg", 1, $FH_HR);
  }

  # if we get here we have no failures
  return $nfinished;
}

#################################################################
# Subroutine: dng_CoordsTokenParse()
# Incept:     EPN, Tue Mar 26 06:15:09 2019
#
# Synopsis: Given a single coords token, validate it, 
#           and return its start, stop, strand values.
# 
# Arguments:
#  $coords_tok:   coordinate token
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    3 values:
#             $start:  start position
#             $stop:   stop position
#             $strand: strand
#
# Dies: if unable to parse $coords
#
#################################################################
sub dng_CoordsTokenParse {
  my $sub_name = "dng_CoordsTokenParse";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords_tok, $FH_HR) = @_;
  if(! defined $coords_tok) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", "dnaorg", 1, $FH_HR); 
  }
  if($coords_tok =~ /^\<?(\d+)\.\.\>?(\d+)\:([\+\-])$/) { 
    return ($1, $2, $3);
  }
  ofile_FAIL("ERROR in $sub_name, unable to parse coords token $coords_tok", "dnaorg", 1, $FH_HR); 

  return; # NEVER REACHED
}

#################################################################
# Subroutine: dng_CoordsLength()
# Incept:     EPN, Tue Mar 26 05:56:08 2019
#
# Synopsis: Given a comma separated coords string, parse it, 
#           validate it, and return its length.
# 
# Arguments:
#  $coords:       coordinate string
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $coords
#
#################################################################
sub dng_CoordsLength {
  my $sub_name = "dng_CoordsLength";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $FH_HR) = @_;
  if(! defined $coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", "dnaorg", 1, $FH_HR); 
  }

  # if there's no comma, we should have a single span
  if($coords !~ m/\,/) { 
    my ($start, $stop, undef) = dng_CoordsTokenParse($coords, $FH_HR);
    return abs($start - $stop) + 1;
  }
  # else, split it up and sum length of all
  my @coords_A  = split(",", $coords);
  my ($start, $stop);
  my $ret_len = 0;
  foreach my $coords_tok (@coords_A) { 
    ($start, $stop, undef) = dng_CoordsTokenParse($coords_tok, $FH_HR);
    $ret_len += abs($start - $stop) + 1;
  }

  return $ret_len;
}

#################################################################
# Subroutine: dng_CoordsFromLocation
# Incept:     EPN, Wed Mar 13 14:17:08 2019
# 
# Purpose:    Convert a GenBank file 'location' value to 
#             a coords string in the format:
#             <start1>-<stop2>:<strand1>,<start2>-<stop2>:<strand2>,...,<startN>-<stopN>:<strandN>
# 
#             This function has to call itself recursively in some
#             cases.
# 
# Arguments:
#   $location: GenBank file location string
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if unable to parse $location
#
# Ref: GenBank release notes (release 230.0) as of this writing
#      and
#      https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
#################################################################
sub dng_CoordsFromLocation { 
  my $sub_name = "dng_CoordsFromLocation";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($location, $FH_HR) = @_;

  # Examples we can parse: 
  # $location                          return value
  # ---------------------------------  -----------------
  # 1..200                             1..200:+
  # <1..200                            1..200:+
  # 100..200>                          100..200:+
  # <1..200>                           1..200:+
  # complement(1..200)                 200..1:-
  # join(1..200,300..400)              1..200:+,300..400:+
  # complement(join(1..200,300..400))  400..300:-,200..1:-
  # join(1..200,complement(300..400))  1..200:+,400..300:- ! NOT SURE IF THIS IS CORRECT !
  # join(complement(300..400),1..200)  400..300:-,1..200:+ ! NOT SURE IF THIS IS CORRECT !

  my $ret_val = "";
  if($location =~ /^join\((.+)\)$/) { 
    my $location_to_join = $1;
    $ret_val = dng_CoordsFromLocation($location_to_join, $FH_HR);
  }
  elsif($location =~ /^complement\((.+)\)$/) { 
    my $location_to_complement = $1;
    my $coords_to_complement = dng_CoordsFromLocation($location_to_complement, $FH_HR);
    $ret_val = dng_CoordsComplement($coords_to_complement, $FH_HR);
  }
  elsif($location =~ /\,/) { 
    # not wrapped in join() or complement(), but multiple segments
    foreach my $location_el (split(",", $location)) { 
      if($ret_val ne "") { $ret_val .= ","; }
      $ret_val .= dng_CoordsFromLocation($location_el, $FH_HR);
    }
  }
  elsif($location =~ /^(\<?\d+\.\.\>?\d+)$/) { 
    $ret_val = $1 . ":+"; # a recursive call due to the complement() may complement this
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse location token $location", "dnaorg", 1, $FH_HR);
  }

  return $ret_val;
}

#################################################################
# Subroutine: dng_CoordsComplement
# Incept:     EPN, Wed Mar 13 15:00:24 2019
# 
# Purpose:    Complement a coords string by complementing all
#             elements within it.
# 
# Arguments:
#   $coords:   coords string to complement
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    complemented $coords
# 
# Dies:       if unable to parse $coords, or any segment in $coords
#             is already on the negative strand.
#
#################################################################
sub dng_CoordsComplement { 
  my $sub_name = "dng_CoordsComplement";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($coords, $FH_HR) = @_;

  # Examples we can parse: 
  # $coords                  return value
  # -----------------------  -----------------
  # 1-200:+                  200-1:-
  # 1-200:+,300-400:+        400-300:-,200-1:-

  my $ret_val = "";
  my @el_A = split(",", $coords);
  for(my $i = scalar(@el_A)-1; $i >= 0; $i--) { 
    if($el_A[$i] =~ /^(\<?)(\d+)\.\.(\>?)(\d+)\:\+/) { 
      my ($start_carrot, $start, $stop_carrot, $stop) = ($1, $2, $3, $4);
      if($start_carrot eq "<") { $start_carrot = ">"; }
      if($stop_carrot  eq ">") { $stop_carrot  = "<"; }
      if($ret_val ne "") { $ret_val .= ","; }
      $ret_val .= $stop_carrot . $stop . ".." . $start_carrot . $start . ":-";
    }
    else { 
      ofile_FAIL("ERROR in $sub_name, unable to parse coords token $coords", "dnaorg", 1, $FH_HR);
    }
  }

  # printf("\tin $sub_name, coords: $coords ret_val: $ret_val\n");

  return $ret_val;
}


#################################################################
# Subroutine: dng_EutilsFetchToFile()
# Incept:     EPN, Tue Mar 12 12:18:37 2019
#
# Synopsis: Fetch information for an accession using edirect.
#
# Arguments:
#  $out_file:  output file to create
#  $accn:      accession to fetch
#  $format:    format to fetch (e.g. "gpc", "ft", "fasta")
#  $nattempts: number of times to retry 
#  $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's a problem fetching the data
#################################################################
sub dng_EutilsFetchToFile { 
  my $sub_name = "dng_EutilsFetchToFile";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_file, $accn, $format, $nattempts, $FH_HR) = @_;
  if((! defined $nattempts) || ($nattempts < 1)) { $nattempts = 1; }

  my $url = dng_EutilsFetchUrl($accn, $format);

  my $n = 0;
  my $fetched_str = undef;
  while(($n < $nattempts) && (! defined $fetched_str)) { 
    $fetched_str = get($url);
    $n++;
    sleep(1);
  }
  if(! defined $fetched_str) { 
    ofile_FAIL("ERROR in $sub_name, problem fetching $accn (undefined)", "dnaorg", 1, $FH_HR); 
  }

  open(OUT, ">", $out_file) || fileOpenFailure($out_file, $sub_name, $!, "writing", $FH_HR);
  print OUT $fetched_str;
  close(OUT);

  return;
}

#################################################################
# Subroutine: dng_EutilsFetchUrl()
# Incept:     EPN, Tue Mar 12 12:18:37 2019
#
# Synopsis: Return a url for an efetch command
#
# Arguments:
#  $accn:      accession to fetch
#  $format:    format to fetch (e.g. "gpc", "ft", "fasta")
#
# Returns:    void
#
# Dies:       never
#################################################################
sub dng_EutilsFetchUrl { 
  my $sub_name = "dng_EutilsFetchUrl";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($accn, $format) = @_;

  return sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=%s&retmode=text", $accn, $format);
}


#################################################################
# Subroutine: dng_ModelInfoFileWrite()
# Incept:     EPN, Sat Mar  9 05:27:15 2019
#
# Synopsis: Output a model info file for all models in @{$mdl_info_AHR}
#           and features in %{$ftr_info_HAHR}.
#
#           The following keys must be set in %{$mdl_info_AHR->[*]}:
#             "name":   name of the model
#             "length": length of the model
#
#           The following keys must be set in %{$ftr_info_HAHR->{*}[*]}:
#             "type":   feature type, e.g. "mat_peptide", "CDS"
#             "coords": coordinates for this feature in the reference
#
# Arguments:
#  $out_file:      out file to create
#  $mdl_info_AHR:  REF to array of hashes with model info, pre-filled
#  $ftr_info_HAHR: REF to hash of array of hashes with information on the features, pre-filled
#  $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $ftr_info_HAHR is not valid upon entering
#################################################################
sub dng_ModelInfoFileWrite { 
  my $sub_name = "dng_ModelInfoFileWrite";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_file, $mdl_info_AHR, $ftr_info_HAHR, $FH_HR) = @_;

  my $mdl_idx  = undef; # model index
  my $ftr_idx  = undef; # feature index
  my $nmdl     = undef; # number of models
  my $nftr     = undef; # number of features
  my $mdl_name = undef; # model name
  my $key      = undef; # a hash key
  my $value    = undef; # a hash value

  my @reqd_mdl_keys_A  = ("name", "length");
  my @reqd_ftr_keys_A  = ("type", "coords");
  # validate all info first
  $nmdl = utl_AHValidate($mdl_info_AHR, \@reqd_mdl_keys_A, "ERROR in $sub_name, mdl_info failed validation; missing required key(s)", $FH_HR);
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    utl_AHValidate(\@{$ftr_info_HAHR->{$mdl_name}}, \@reqd_ftr_keys_A, "ERROR in $sub_name, ftr_info failed validation; missing required key(s)", $FH_HR);
  }

  # verify feature coords make sense
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    dng_FeatureInfoValidateCoords($ftr_info_HAHR->{$mdl_name}, $mdl_info_AHR->[$mdl_idx]{"length"}, $FH_HR); 
  }

  # output 
  open(OUT, ">", $out_file) || fileOpenFailure($out_file, $sub_name, $!, "writing", $FH_HR);
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    print OUT ("MODEL $mdl_name");
    foreach $key (sort keys (%{$mdl_info_AHR->[$mdl_idx]})) { 
      $value = $mdl_info_AHR->[$mdl_idx]{$key};
      if($key =~ m/\:/) { 
        ofile_FAIL("ERROR in $sub_name, problem writing $out_file, illegal ':' character in model key $key for model $mdl_name", "dnaorg", 1, $FH_HR);
      }
      if($value =~ m/\"/) { 
        ofile_FAIL("ERROR in $sub_name, problem writing $out_file, illegal '\"' character in model value $value for key $key for model $mdl_name", "dnaorg", 1, $FH_HR);
      }
      if($key ne "name") { 
        print OUT (" $key:\"$value\"");
      }
    }
    print OUT ("\n");

    # define feature keys to ignore
    my %ftr_key_ignore_H = ();
    $ftr_key_ignore_H{"type"}        = 1; # this automatically gets added to @key_order_A, so it goes first
    $ftr_key_ignore_H{"coords"}      = 1; # this automatically gets added to @key_order_A, so it goes second
    $ftr_key_ignore_H{"length"}      = 1; # will be inferred from coords
    $ftr_key_ignore_H{"source_idx"}  = 1; # will be inferred from coords and type
    $ftr_key_ignore_H{"parent_idx"}  = 1; # will be inferred from coords and type
    $ftr_key_ignore_H{"3pa_ftr_idx"} = 1; # will be inferred from coords and type
    $ftr_key_ignore_H{"5p_sgm_idx"}  = 1; # will be inferred from coords, when sgm_info_HA is created
    $ftr_key_ignore_H{"3p_sgm_idx"}  = 1; # will be inferred from coords, when sgm_info_HA is created
    $ftr_key_ignore_H{"location"}    = 1; # *could* (but won't be) inferred from coords

    $nftr = scalar(@{$ftr_info_HAHR->{$mdl_name}});
    # determine order of keys for this feature
    my @ftr_key_order_A  = ("type", "coords");
    for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      foreach $key (sort keys %{$ftr_info_HAHR->{$mdl_name}[$ftr_idx]}) { 
        if(! exists $ftr_key_ignore_H{$key}) { 
          push(@ftr_key_order_A, $key);
          $ftr_key_ignore_H{$key} = 1; 
        }
      }
    }
    
    for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      print OUT ("FEATURE $mdl_name");
      foreach $key (@ftr_key_order_A) { 
        if(exists $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key}) { 
          $value = $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key};
          if($key =~ m/\:/) { 
            ofile_FAIL("ERROR in $sub_name, problem writing $out_file, illegal ':' character in feature key $key for model $mdl_name", "dnaorg", 1, $FH_HR);
          }
          if($value =~ m/\"/) { 
            ofile_FAIL("ERROR in $sub_name, problem writing $out_file, illegal '\"' character in feature value $value for key $key for model $mdl_name", "dnaorg", 1, $FH_HR);
          }
          print OUT (" $key:\"$value\"");
        }
      }
      print OUT ("\n");
    }
  }
  close(OUT);

  return;
}

#################################################################
# Subroutine: dng_ModelInfoFileParse()
# Incept:     EPN, Fri Mar 15 05:15:23 2019
#
# Synopsis: Parse a model info file for >= 1 models and collect 
#           feature information for each model $model in 
#           @{$ftr_info_HAHR->{$model}}.
#
#           The following keys must be defined for all features:
#             "type":   feature type, e.g. "mat_peptide", "CDS"
#             "coords": coordinates for this feature in the reference
#           We verify this at end of subroutine
# 
# Arguments:
#  $in_file:          input .minfo file to parse
#  $mdl_info_AHR:     REF to array of hashes of model information, filled here
#  $ftr_info_HAHR:    REF to hash of array of hashes with information 
#                     on the features per model, filled here
#  $FH_HR:            REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if unable to parse $in_file
#             if a feature is defined without "type" or "coords" keys
#################################################################
sub dng_ModelInfoFileParse {
  my $sub_name = "dng_ModelInfoFileParse";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($in_file, $mdl_info_AHR, $ftr_info_HAHR, $FH_HR) = @_;
  
  my $format_str = "# DNAORG model info (.minfo) format specifications:\n";
  $format_str   .= "# Lines prefixed with '#' are ignored.\n";
  $format_str   .= "# All other lines must begin with either: 'MODEL' or 'FEATURE'\n";
  $format_str   .= "# followed by one or more whitespace characters and then the model\n";
  $format_str   .= "# name <modelname> which cannot include whitespace.\n";
  $format_str   .= "# On each line after <modelname>, both MODEL and FEATURE lines must\n";
  $format_str   .= "# contain 0 or more <key>:<value> pairs meeting the following criteria.\n";
  $format_str   .= "# <key> must not include any whitespace or ':' characters\n";
  $format_str   .= "# <value> must start *and* end with '\"' but include no other '\"'\n";
  $format_str   .= "# characters (but <value> may include whitespace characters).\n";
  $format_str   .= "# <key>:<value> pairs must be separated by one or more whitespace characters.\n";
  $format_str   .= "# <modelname> and the first <key>:<value> pair must be separated by one or\n";
  $format_str   .= "# more whitespace characters.\n";

  # example lines:
  #MODEL NC_039477 cmfile:"test/test.dnaorg_build.cm"
  #FEATURE NC_039477 type:"gene" coords:"5..5104:+" gene:"ORF1"
  #FEATURE NC_039477 type:"CDS" coords:"5..5104:+" gene:"ORF1" product:"nonstructural polyprotein"

  my $mdl_name   = undef; # name of current model
  my $ftr_idx    = undef; # index of current feature
  my $mdl_idx    = -1;    # index of current model
  my %mdl_read_H = ();    # keeps track of which model names we've seen MODEL lines for, to avoid duplicates
  open(IN, $in_file) || fileOpenFailure($in_file, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    if($line !~ /^#/) { 
      # not a comment line
      my $orig_line = $line;
      chomp $line;
      my $is_model_line = 0; # set to 1 if line we are parsing is a MODEL line, else it's a FEATURE line
      if($line =~ /^MODEL\s+(\S+)\s*/) { 
        $mdl_name = $1;
        if(exists $mdl_read_H{$mdl_name}) { 
          ofile_FAIL("ERROR in $sub_name, problem parsing $in_file: read multiple MODEL lines for $mdl_name, should only be 1; line:\n$orig_line\n", "dnaorg", 1, $FH_HR);
        }
        $mdl_idx++;
        %{$mdl_info_AHR->[$mdl_idx]} = ();
        @{$ftr_info_HAHR->{$mdl_name}} = ();
        $mdl_info_AHR->[$mdl_idx]{"name"} = $mdl_name;
        $mdl_read_H{$mdl_name} = 1;

        $is_model_line = 1;
        $line =~ s/^MODEL\s+(\S+)\s*//; # remove MODEL and model value
      }
      elsif($line =~ /^FEATURE\s+(\S+)\s*/) { 
        $mdl_name = $1;
        if(! exists $mdl_read_H{$mdl_name}) { 
          ofile_FAIL("ERROR in $sub_name, problem parsing $in_file: read FEATURE line for model $mdl_name before a MODEL line for $mdl_name; line:\n$orig_line\n", "dnaorg", 1, $FH_HR);
        }
        $ftr_idx = scalar(@{$ftr_info_HAHR->{$mdl_name}});
        # initialize ftr_info for this model/feature pair
        %{$ftr_info_HAHR->{$mdl_name}[$ftr_idx]} = (); 
        $line =~ s/^FEATURE\s+\S+\s*//; # remove FEATURE and model value
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $in_file, non-comment line does not start with 'MODEL <modelname>' or 'FEATURE <featurename>', line:\n$orig_line\n", "dnaorg", 1, $FH_HR);
      }
      # if we get here we have either a MODEL or FEATURE line, parse the rest of it
      while($line ne "") { 
        if($line =~ /^([^\:\s]+)\:\"([^\"]+)\"\s*/) { 
          # key   must not include ':' or whitespace
          # value must begin and end with '"' but otherwise include on '"' characters
          my ($key, $value) = ($1, $2);
          if($is_model_line) { 
            if(exists $mdl_info_AHR->[$mdl_idx]{$key}) {
              ofile_FAIL("ERROR in $sub_name, problem parsing $in_file, read multiple values for key $key on MODEL line; line:\n$orig_line\n", "dnaorg", 1, $FH_HR);
            }
            $mdl_info_AHR->[$mdl_idx]{$key} = $value;
          }
          else { # feature line
            if(exists $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key}) {
              ofile_FAIL("ERROR in $sub_name, problem parsing $in_file, read multiple values for key $key on MODEL line; line:\n$orig_line\n", "dnaorg", 1, $FH_HR);
            }
            $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key} = $value;
            # printf("\tadded ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key} as $value\n");
          }
          $line =~ s/^[^\:\s]+\:\"[^\"]+\"\s*//; # remove this key/value pair
        }
        else { 
          ofile_FAIL("ERROR in $sub_name, unable to parse $in_file, failed to parse key:value pairs in line:\n$orig_line\n$format_str\n", "dnaorg", 1, $FH_HR);
        }
      } 
    }
  }
  close(IN);

  # verify we read what we need
  my @reqd_mdl_keys_A = ("name", "length");
  my @reqd_ftr_keys_A = ("type", "coords");
  utl_AHValidate($mdl_info_AHR, \@reqd_mdl_keys_A, "ERROR in $sub_name, problem parsing $in_file, required MODEL key missing", $FH_HR);
  my $nmdl = scalar(@{$mdl_info_AHR});
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    utl_AHValidate($ftr_info_HAHR->{$mdl_name}, \@reqd_ftr_keys_A, "ERROR in $sub_name, problem parsing $in_file, required MODEL key missing for model " . $mdl_info_AHR->[$mdl_idx]{"name"}, $FH_HR);
  }

  # verify feature coords make sense
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    dng_FeatureInfoValidateCoords($ftr_info_HAHR->{$mdl_name}, $mdl_info_AHR->[$mdl_idx]{"length"}, $FH_HR); 
  }
  return;
}

#################################################################
# Subroutine:  dng_SplitFastaFile()
# Incept:      EPN, Tue Mar  1 09:30:10 2016
#
# Purpose: Split up a fasta file into <n> smaller files by calling
#          the esl-ssplit perl script.
#
# Arguments: 
#  $esl_ssplit:      path to the esl-ssplit.pl script to use
#  $fasta_file:      fasta file to split up
#  $nfiles:          desired number of files to split $fasta_file into, -1 for one file for each sequence
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:    Number of files actually created (can differ from requested
#             amount (which is $nfiles)).
#
# Dies:       if esl-ssplit command fails
#
################################################################# 
sub dng_SplitFastaFile { 
  my $sub_name = "dng_SplitFastaFile()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_ssplit, $fasta_file, $nfiles, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to ofile_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $outfile = $fasta_file . ".esl-ssplit";
  my $cmd = undef;
  if($nfiles == -1) { # special case: put 1 file per sequence
    $cmd = "$esl_ssplit -v $fasta_file 1 > $outfile";
  }
  else { 
    $cmd = "$esl_ssplit -v -r -n $fasta_file $nfiles > $outfile";
  }
  utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # parse output to determine exactly how many files were created:
  # $esl_ssplit will have output exactly 1 line per fasta file it created
  my $nfiles_created = utl_FileCountLines($outfile, $FH_HR);

  if(! opt_Get("--keep", $opt_HHR)) { 
    utl_RunCommand("rm $outfile", opt_Get("-v", $opt_HHR), 0, $FH_HR);
  }

  return $nfiles_created;
}

#################################################################
# Subroutine: dng_SplitNumSeqFiles()
# Incept:     EPN, Mon Mar 18 15:01:44 2019
#
# Synopsis: Return number of sequence files we need to split a sequence
#           file of $tot_nt nucleotides into based on --nkb and 
#           --maxnjobs options in %{$opt_HHR}.
#
# Arguments:
#  $tot_nt:   total number of nucleotides in the file
#  $opt_HHR:  REF to 2D hash of option values
#
# Returns:    Number of sequence files (>= 1).
#
# Dies:       If --nkb or --maxnjobs options do not exist in %{$opt_HHR}.
#################################################################
sub dng_SplitNumSeqFiles { 
  my $sub_name = "dng_SplitNumSeqFiles";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($tot_nt, $opt_HHR) = @_;

  my $ret_nseqfiles = int($tot_nt / (opt_Get("--nkb", $opt_HHR) * 1000)); 
  # int() takes the floor, so there can be a nonzero remainder. We don't add 1 though, 
  # because splitFastaFile() will return the actual number of sequence files created
  # and we'll use that as the number of jobs subsequently. $targ_nseqfiles is currently only
  # the 'target' number of sequence files that we pass into splitFastaFile().
  # make sure we won't exceed our max number of jobs (from --maxnjobs)
  if(($ret_nseqfiles) > (opt_Get("--maxnjobs", $opt_HHR))) { 
    $ret_nseqfiles = int(opt_Get("--maxnjobs", $opt_HHR));
  }
  if($ret_nseqfiles == 0) { $ret_nseqfiles = 1; }

  return $ret_nseqfiles;
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
          if($astart > $astop) { utl_Swap(\$astart, \$astop); }
          my $sgm_sqstring = $msa->get_sqstring_unaligned_and_truncated($seq_idx, $astart, $astop);
          if($sgm_strand_AA[$ftr_idx][$sgm_idx] eq "-") { 
            seq_SqstringReverseComplement(\$sgm_sqstring);
          }
          $cds_sqstring .= $sgm_sqstring;
        }
        print $out_FH(">" . $msa->get_sqname($seq_idx) . "/" . $ftr_info_AHR->[$ftr_idx]{"coords"} . "\n" . seq_SqstringAddNewlines($cds_sqstring, 60));
      }
    }
  }
  return;
}

#################################################################
# Subroutine: dng_StripVersion()
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
sub dng_StripVersion {
  my $sub_name  = "dng_StripVersion()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($accver_R) = (@_);

  $$accver_R =~ s/\.[0-9]*$//; # strip version

  return;
}

###########################################################################
# the next line is critical, a perl module must return a true value
return 1;
###########################################################################
