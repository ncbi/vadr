#!/usr/bin/perl
# 
# version: 1.0.5 [March 2020]
#
# vadr.pm
# Eric Nawrocki
# EPN, Mon Feb  1 15:23:00 2016
# 
# Some subroutines are taken from rnavore.pm [EPN, Tue Sep 23 09:22:55 2014]
# and ssu.pm [EPN, Thu Nov  5 05:39:37 2009].
#
# Perl module used by vadr-build.pl and vadr-annotate.pl scripts
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
use LWP::Simple; 

require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seq.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

#########################################################################################
#
# Subroutines related to features or segments:
# vdr_FeatureInfoImputeCoords()
# vdr_FeatureInfoImputeLength()
# vdr_FeatureInfoInitializeParentIndexStrings()
# vdr_FeatureInfoImputeOutname()
# vdr_FeatureInfoImpute3paFtrIdx()
# vdr_FeatureInfoImputeByOverlap()
# vdr_FeatureInfoStartStopStrandArrays()
# vdr_FeatureInfoCountType()
# vdr_FeatureInfoValidateCoords()
# vdr_FeatureInfoValidateParentIndexStrings()
# vdr_FeatureInfoChildrenArrayOfArrays()
# vdr_FeatureInfoMapFtrTypeIndicesToFtrIndices()
#
# vdr_SegmentInfoPopulate()
# 
# vdr_FeatureTypeAndTypeIndexString()
# vdr_FeatureTypeIndex()
# vdr_FeatureTypeIsCds()
# vdr_FeatureTypeIsMatPeptide()
# vdr_FeatureTypeIsGene()
# vdr_FeatureTypeIsCdsOrMatPeptide()
# vdr_FeatureTypeIsCdsOrMatPeptideOrGene()
# vdr_FeatureNumSegments()
# vdr_FeatureRelativeSegmentIndex()
# vdr_Feature5pMostPosition()
# vdr_Feature3pMostPosition()
# vdr_FeatureSummarizeSegment()
# vdr_FeatureParentIndex()
# vdr_FeatureStartStopStrandArrays()
# vdr_FeatureSummaryStrand()
# vdr_FeaturePositionSpecificValueBreakdown()
# 
# Subroutines related to alerts:
# vdr_AlertInfoInitialize()
# vdr_AlertInfoAdd()
# vdr_AlertInfoSetFTableInvalidatedBy()
# vdr_AlertInfoSetCausesFailure()
# vdr_AlertInfoDump()
# 
# Subroutines related to parallelization on the compute farm:
# vdr_ParseQsubFile()
# vdr_SubmitJob()
# vdr_WaitForFarmJobsToFinish()
#
# Subroutines related to sequence and model coordinates: 
# vdr_CoordsSegmentParse()
# vdr_CoordsSegmentCreate()
# vdr_CoordsAppendSegment()
# vdr_CoordsLength()
# vdr_CoordsFromLocation()
# vdr_CoordsReverseComplement()
# vdr_CoordsSegmentReverseComplement()
# vdr_CoordsMin()
# vdr_CoordsMax()
# vdr_CoordsMissing()
# vdr_CoordsCheckIfSpans()
# vdr_CoordsSegmentOverlap()
# vdr_CoordsRelativeToAbsolute()
# vdr_CoordsRelativeSegmentToAbsolute()
# vdr_CoordsProteinRelativeToAbsolute()
# vdr_CoordsProteinToNucleotide()
# vdr_CoordsMergeAllAdjacentSegments()
# vdr_CoordsMergeTwoSegmentsIfAdjacent()
# vdr_CoordsMaxLengthSegment()
#
# Subroutines related to eutils:
# vdr_EutilsFetchToFile()
# vdr_EutilsFetchUrl()
# 
# Subroutines related to model info files:
# vdr_ModelInfoFileWrite()
# vdr_ModelInfoFileParse()
#
# Subroutines related to cmalign output files:
# vdr_CmalignParseIfile()
# vdr_CmalignWriteIfile()
# 
# Miscellaneous subroutines:
# vdr_SplitFastaFile()
# vdr_SplitNumSeqFiles()
# vdr_CdsFetchStockholmToFasta()
#
#################################################################
# Subroutine: vdr_FeatureInfoImputeCoords
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
sub vdr_FeatureInfoImputeCoords { 
  my $sub_name = "vdr_FeatureInfoImputeCoords";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }

  my $do_carrots = 0; # do not include carrots in coords strings, even if they exist in location strings

  my ($ftr_info_AHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "location"
  my @keys_A = ("type", "location");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"coords"} = vdr_CoordsFromLocation($ftr_info_AHR->[$ftr_idx]{"location"}, $do_carrots, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoImputeLength
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
sub vdr_FeatureInfoImputeLength { 
  my $sub_name = "vdr_FeatureInfoImputeLength";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  # go through all features and determine length by parsing the 
  # "coords" value
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"length"} = vdr_CoordsLength($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoInitializeParentIndexStrings
# Incept:     EPN, Wed Mar 13 13:33:33 2019
# 
# Purpose:    Set "parent_idx_str" value to "GBNULL" for any feature 
#             in which it is not already defined in @{$ftr_info_AHR}
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
sub vdr_FeatureInfoInitializeParentIndexStrings {
  my $sub_name = "vdr_FeatureInfoInitializeParentIndexStrings";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(! defined $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"}) { 
      $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} = "GBNULL";
    }
  }

  return;
}
      
#################################################################
# Subroutine: vdr_FeatureInfoImputeOutname()
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
sub vdr_FeatureInfoImputeOutname { 
  my $sub_name  = "vdr_FeatureInfoImputeOutname";
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
      $ftr_info_AHR->[$ftr_idx]{"outname"} = vdr_FeatureTypeAndTypeIndexString($ftr_info_AHR, $ftr_idx, ".");
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoImpute3paFtrIdx
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
sub vdr_FeatureInfoImpute3paFtrIdx {
  my $sub_name = "vdr_FeatureInfoImpute3paFtrIdx";
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
      $ftr_3p_pos = vdr_Feature3pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
      $ftr_strand = vdr_FeatureSummaryStrand($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
      for($ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
        $ftr_5p_pos2 = vdr_Feature5pMostPosition($ftr_info_AHR->[$ftr_idx2]{"coords"}, $FH_HR);
        $ftr_strand2 = vdr_FeatureSummaryStrand($ftr_info_AHR->[$ftr_idx2]{"coords"}, $FH_HR);
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
            if($ftr_info_AHR->[$ftr_idx]{"3pa_ftr_idx"} == -1) { 
              # if multiple features are 3'-adjacent, keep only the first one
              $ftr_info_AHR->[$ftr_idx]{"3pa_ftr_idx"} = $ftr_idx2; 
            }
          }
        }
      }
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoImputeByOverlap
# Incept:     EPN, Tue May 21 06:28:07 2019
# 
# Purpose:    Add a qualifier and value to instances of one feature
#             type based on qualifier and values of another feature
#             type and overlap of nucleotides
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $src_type:      source feature type
#   $src_key:       source feature key (value in %{$ftr_info_AH[$src_idx]})
#   $dst_type:      destination feature type
#   $dst_key:       destination feature key (value in %{$ftr_info_AH[$dst_idx]})
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub vdr_FeatureInfoImputeByOverlap {
  my $sub_name = "vdr_FeatureInfoImputeByOverlap";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $src_type, $src_key, $dst_type, $dst_key, $FH_HR) = @_;

  # printf("in $sub_name, src_type: $src_type, src_key: $src_key, dst_type: $dst_type, dst_key: $dst_key\n");
  # contract check
  if($src_type eq $dst_type) { 
    ofile_FAIL("ERROR in $sub_name, source type and destination type are equal ($src_type)\n", 1, $FH_HR);
  }

  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  # - go through all features and find F1 of type $dst_type
  #    - if F1 does not have a value for $dst_key
  #      find other features F2 of type $src_type that 'span' F1
  #      where spans means every position in F1 also exists in F2
  #      - set F1 $dst_key to value in F2 for $src_key
  # 
  # if more than one F2 exists for an F1 that meet criteria above:
  #     - set F1 $dst_key to value in F2 for $src_key for shortest F2
  #       
  # if more than one F2 exists with identical coords for an F1 that meet criteria above
  #     - die
  #
  my ($dst_ftr_idx, $src_ftr_idx); # feature indices
  for($dst_ftr_idx = 0; $dst_ftr_idx < $nftr; $dst_ftr_idx++) { 
    if(($ftr_info_AHR->[$dst_ftr_idx]{"type"} eq $dst_type) && 
       (! defined $ftr_info_AHR->[$dst_ftr_idx]{$dst_key})) { 
      my $found_src_ftr_idx = undef;
      for($src_ftr_idx = 0; $src_ftr_idx < $nftr; $src_ftr_idx++) { 
        if(($dst_ftr_idx != $src_ftr_idx) && 
           ($ftr_info_AHR->[$src_ftr_idx]{"type"} eq $src_type) && 
           (defined $ftr_info_AHR->[$src_ftr_idx]{$src_key})) { 
          # determine if $dst_ftr_idx is completely spanned by $src_ftr_idx overlap completely
          # printf("checking overlap between src " . $ftr_info_AHR->[$src_ftr_idx]{"coords"} . " and dst " . $ftr_info_AHR->[$dst_ftr_idx]{"coords"} . "\n");
          if(vdr_CoordsCheckIfSpans($ftr_info_AHR->[$src_ftr_idx]{"coords"}, $ftr_info_AHR->[$dst_ftr_idx]{"coords"}, $FH_HR)) { 
            # $src_ftr_idx completely spans $dst_ftr_idx
            my $do_update = 1; # set to '0' below if $found_src_ftr_idx is defined and is shorter than $src_ftr_idx
            if(defined $found_src_ftr_idx) { 
              # this is okay, pick the shorter of the two
              # (but if they have identical coords values --> die)
              my $existing_coords = $ftr_info_AHR->[$found_src_ftr_idx]{"coords"};
              my $current_coords  = $ftr_info_AHR->[$src_ftr_idx]{"coords"};
              if($existing_coords eq $current_coords) { 
                ofile_FAIL("ERROR in $sub_name, more than one feature of type $src_type (with key $src_key) of identical coords ($existing_coords) spans feature idx $dst_ftr_idx of type $dst_type with coords " . $ftr_info_AHR->[$dst_ftr_idx]{"coords"} . "\n1: idx: $found_src_ftr_idx, value: " . $ftr_info_AHR->[$found_src_ftr_idx]{$src_key} . "\n2: idx: $src_ftr_idx, value: " . $ftr_info_AHR->[$src_ftr_idx]{$src_key} . "\n", 1, $FH_HR);
              }
              # if we get here, the coords differ, only update if current is smaller than existing
              if((vdr_CoordsLength($current_coords, $FH_HR)) >= (vdr_CoordsLength($existing_coords, $FH_HR))) { 
                # printf("in $sub_name, skipping update of ftr_info_AHR->[$dst_ftr_idx]{$dst_key} of type $dst_type to feature idx $src_ftr_idx with value $ftr_info_AHR->[$src_ftr_idx]{$src_key} and coords " . $ftr_info_AHR->[$src_ftr_idx]{"coords"} . "\nbecause shorter overlap exists: feature $found_src_ftr_idx with value $ftr_info_AHR->[$found_src_ftr_idx]{$src_key} and coords . " . $ftr_info_AHR->[$found_src_ftr_idx]{"coords"} . "\n");
                $do_update = 0;
              }
            } # end of 'if(defined $found_src_ftr_idx)'
            if($do_update) { # will be '1' if (! defined $found_src_ftr_idx) or if $src_ftr_idx is shorter than $found_src_ftr_idx
              # printf("in $sub_name, setting ftr_info_AHR->[$dst_ftr_idx]{$dst_key} of type $dst_type to $ftr_info_AHR->[$src_ftr_idx]{$src_key}\n");
              $ftr_info_AHR->[$dst_ftr_idx]{$dst_key} = $ftr_info_AHR->[$src_ftr_idx]{$src_key};
              $found_src_ftr_idx = $src_ftr_idx;
            }
          }
        }
      }
    }
  }
  
  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoStartStopStrandArrays()
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
sub vdr_FeatureInfoStartStopStrandArrays {
  my $sub_name = "vdr_FeatureInfoStartStopStrandArrays";
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
    vdr_FeatureStartStopStrandArrays($ftr_info_AHR->[$ftr_idx]{"coords"}, \@{$start_AA[$ftr_idx]}, \@{$stop_AA[$ftr_idx]}, \@{$strand_AA[$ftr_idx]}, $FH_HR);
  }
  if(defined $start_AAR)  { @{$start_AAR}   = @start_AA;   }
  if(defined $stop_AAR)   { @{$stop_AAR}    = @stop_AA;    }
  if(defined $strand_AAR) { @{$strand_AAR}  = @strand_AA;  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoCountType
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
sub vdr_FeatureInfoCountType { 
  my $sub_name = "vdr_FeatureInfoCountType";
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
# Subroutine: vdr_FeatureInfoValidateCoords
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
sub vdr_FeatureInfoValidateCoords { 
  my $sub_name = "vdr_FeatureInfoValidateCoords";
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
    # this sub vdr_Will die if $ftr_info_AHR->[$ftr_idx]{"coords"} is in incorrect format
    vdr_FeatureStartStopStrandArrays($ftr_info_AHR->[$ftr_idx]{"coords"}, \@start_A, \@stop_A, undef, $FH_HR); 
    foreach my $start (@start_A) { if($start > $length) { $fail_str .= "ftr_idx: $ftr_idx, start position $start > $length\n"; } }
    foreach my $stop  (@stop_A)  { if($stop  > $length) { $fail_str .= "ftr_idx: $ftr_idx, stop  position $stop  > $length\n"; } }
  }

  if($fail_str ne "") { 
    ofile_FAIL("ERROR in $sub_name, some coordinates exceed model length ($length):\n$fail_str\n", 1, $FH_HR);
  }
  
  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoValidateParentIndexStrings
# Incept:     EPN, Wed Feb 19 11:44:28 2020
# 
# Purpose:    Validate "parent_idx_str" values are either "GBNULL"
#             or a valid feature index [0..$nftr-1]. Should probably
#             be called after vdr_FeatureInfoInitializeParentIndexStrings().
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
sub vdr_FeatureInfoValidateParentIndexStrings {
  my $sub_name = "vdr_FeatureInfoValidateParentIndexStrings";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
  
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  my $nftr = scalar(@{$ftr_info_AHR});
  my $fail_str = ""; # added to if any elements are out of range
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(! defined $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"}) { 
      $fail_str .= "ftr_idx: $ftr_idx, undefined\n"; 
    }
    elsif($ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} ne "GBNULL") { 
      if($ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} < 0) { 
        $fail_str .= "ftr_idx: $ftr_idx, " . $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} . " < 0\n"; 
      }
      elsif($ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} >= $nftr) { 
        $fail_str .= "ftr_idx: $ftr_idx, " . $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} . " >= $nftr (num features, should be 0.." . ($nftr-1) . ")\n";
      }
    }
  }
  
  if($fail_str ne "") { 
    ofile_FAIL("ERROR in $sub_name, some parent index strings are undefined or don't make sense:\n$fail_str\n", 1, $FH_HR);
  }

  return;
  
}

#################################################################
# Subroutine:  vdr_FeatureInfoChildrenArrayOfArrays()
# Incept:      EPN, Sun Mar 10 06:22:49 2019
#
# Purpose:     Fill @{$AAR} with arrays of children (feature indices)
#              for each feature in %{$ftr_info_AHR}.
# 
# Arguments: 
#   $ftr_info_AHR:   REF to hash of arrays with information on the features, PRE-FILLED
#   $type_or_undef:  feature type of children (e.g. mat_peptide) we want information on
#                    caller should set as 'undef' to get information on all types of children
#   $AAR:            REF to array of arrays of children feature indices, FILLED HERE
#   $FH_HR:          REF to hash of file handles
# 
# Returns:     Nothing.
# 
#
################################################################# 
sub vdr_FeatureInfoChildrenArrayOfArrays { 
  my $nargs_expected = 4;
  my $sub_name = "vdr_FeatureInfoChildrenArrayOfArrays";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $type_or_undef, $AAR, $FH_HR) = @_;

  @{$AAR} = ();
  my $nftr = scalar(@{$ftr_info_AHR});

  my ($parent_ftr_idx, $child_ftr_idx);

  # initialize
  for($parent_ftr_idx = 0; $parent_ftr_idx < $nftr; $parent_ftr_idx++) { 
    @{$AAR->[$parent_ftr_idx]} = ();
  }

  # fill
  for($child_ftr_idx = 0; $child_ftr_idx < $nftr; $child_ftr_idx++) { 
    if((! defined $type_or_undef) || ($ftr_info_AHR->[$child_ftr_idx]{"type"} eq $type_or_undef)) { 
      if($ftr_info_AHR->[$child_ftr_idx]{"parent_idx_str"} ne "GBNULL") { 
        my @parent_ftr_idx_A = split(",", $ftr_info_AHR->[$child_ftr_idx]{"parent_idx_str"});
        foreach $parent_ftr_idx (@parent_ftr_idx_A) { 
          push(@{$AAR->[$parent_ftr_idx]}, $child_ftr_idx);
        }
      }
    }
  }
  return;
}

#################################################################
# Subroutine:  vdr_FeatureInfoMapFtrTypeIndicesToFtrIndices()
# Incept:      EPN, Wed Mar 18 08:15:06 2020
#
# Purpose:    Map ftr_type_idx values (e.g. CDS.4) to feature indices
#             (ftr_idx) in \@{$ftr_info_AHR}.
#
# Arguments: 
#  $ftr_info_AHR:            REF to array of hashes with information on the features, PRE-FILLED
#  $ftr_type_idx2ftr_idx_HR: REF to hash to fill, key is ftr_type_idx (e.g. CDS.4) value is 
#  $FH_HR:                   REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#
# Returns:    void
#
################################################################# 
sub vdr_FeatureInfoMapFtrTypeIndicesToFtrIndices {
  my $sub_name = "vdr_FeatureInfoMapFtrTypeIndicesToFtrIndices";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $ftr_type_idx2ftr_idx_HR, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $ftr_type_idx = vdr_FeatureTypeAndTypeIndexString($ftr_info_AHR, $ftr_idx, ".");
    $ftr_type_idx2ftr_idx_HR->{$ftr_type_idx} = $ftr_idx;
  }
 
  return;
}

#################################################################
# Subroutine: vdr_SegmentInfoPopulate()
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
sub vdr_SegmentInfoPopulate {
  my $sub_name = "vdr_SegmentInfoPopulate";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($sgm_info_AHR, $ftr_info_AHR, $FH_HR) = @_;

  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR in $sub_name", $FH_HR);

  # initialize new %{$ftr_info_AHR} values
  my ($ftr_idx, $ftr_idx2, $sgm_idx, $sgm_idx2); # feature and segment indices
  my ($sgm_start, $sgm_stop, $sgm_strand); # start, stop and strand for a segment
  my $nseg = 0; 
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} = -1; 
    $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"} = -2; 
    my $ftr_type   = $ftr_info_AHR->[$ftr_idx]{"type"};

    # determine start and stop positions of all segments
    my @sgm_start_A  = (); # array of starts, one per segment
    my @sgm_stop_A   = (); # array of stops, one per segment
    my @sgm_strand_A = (); # array of strands ("+", "-"), one per segment
    vdr_FeatureStartStopStrandArrays($ftr_info_AHR->[$ftr_idx]{"coords"}, \@sgm_start_A, \@sgm_stop_A, \@sgm_strand_A, $FH_HR);
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

  return;
}

#################################################################
# Subroutine:  vdr_FeatureTypeAndTypeIndexString()
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
sub vdr_FeatureTypeAndTypeIndexString { 
  my $nargs_expected = 3;
  my $sub_name = "vdr_FeatureTypeAndTypeIndexString";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $ftr_idx, $sep_char) = @_;

  my $type = $ftr_info_AHR->[$ftr_idx]{"type"};

  return $type . $sep_char . vdr_FeatureTypeIndex($ftr_info_AHR, $ftr_idx);
}

#################################################################
# Subroutine:  vdr_FeatureTypeIndex()
# Incept:      EPN, Tue Feb 11 11:49:05 2020
#
# Purpose:     Return the type index for $ftr_idx.
# 
# Arguments: 
#   $ftr_info_AHR:   REF to hash of arrays with information on the features, PRE-FILLED
#   $ftr_idx:        index we are interested in
# 
# Returns:     Type index, string
#
################################################################# 
sub vdr_FeatureTypeIndex { 
  my $nargs_expected = 2;
  my $sub_name = "vdr_FeatureTypeIndex";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $ftr_idx) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  my $type = $ftr_info_AHR->[$ftr_idx]{"type"};
  my $type_idx = 1;
  for(my $ftr_idx2 = 0; $ftr_idx2 < $ftr_idx; $ftr_idx2++) { 
    if($ftr_info_AHR->[$ftr_idx2]{"type"} eq $type) { $type_idx++; }
  }
  
  return $type_idx;
}

#################################################################
# Subroutine: vdr_FeatureTypeIsCds()
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
sub vdr_FeatureTypeIsCds { 
  my $sub_name = "vdr_FeatureTypeIsCds";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return ($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") ? 1 : 0;
}

#################################################################
# Subroutine: vdr_FeatureTypeIsMatPeptide()
# Incept:     EPN, Fri Apr  5 14:32:28 2019
#
# Purpose:    Is feature $ftr_idx a mat_peptide?
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
sub vdr_FeatureTypeIsMatPeptide { 
  my $sub_name = "vdr_FeatureTypeIsMatPeptide";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return ($ftr_info_AHR->[$ftr_idx]{"type"} eq "mat_peptide") ? 1 : 0;
}

#################################################################
# Subroutine: vdr_FeatureTypeIsGene()
# Incept:     EPN, Fri Apr  5 14:32:54 2019
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
sub vdr_FeatureTypeIsGene { 
  my $sub_name = "vdr_FeatureTypeIsGene";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return ($ftr_info_AHR->[$ftr_idx]{"type"} eq "gene") ? 1 : 0;
}

#################################################################
# Subroutine: vdr_FeatureTypeIsCdsOrMatPeptide()
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
sub vdr_FeatureTypeIsCdsOrMatPeptide { 
  my $sub_name = "vdr_FeatureTypeIsCdsOrMatPeptide";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return (($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") || 
          ($ftr_info_AHR->[$ftr_idx]{"type"} eq "mat_peptide")) ? 1 : 0;
}

#################################################################
# Subroutine: vdr_FeatureTypeIsCdsOrMatPeptideOrGene()
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
sub vdr_FeatureTypeIsCdsOrMatPeptideOrGene { 
  my $sub_name = "vdr_FeatureTypeIsCdsOrMatPeptideOrGene";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return (($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") || 
          ($ftr_info_AHR->[$ftr_idx]{"type"} eq "mat_peptide") ||
          ($ftr_info_AHR->[$ftr_idx]{"type"} eq "gene")) ? 1 : 0;
}

#################################################################
# Subroutine: vdr_FeatureNumSegments()
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
sub vdr_FeatureNumSegments { 
  my $sub_name  = "vdr_FeatureNumSegments";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_AHR, $ftr_idx) = (@_);

  return ($ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"} - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} + 1);
}


#################################################################
# Subroutine: vdr_FeatureRelativeSegmentIndex()
# Incept:     EPN, Tue Feb 11 11:53:02 2020
#
# Purpose:    Return the relative index of segment $sgm_idx 
#             in $ftr_idx.
#
# Arguments: 
#   $ftr_info_AHR:  REF to array of hashes of feature info
#   $ftr_idx:       feature index we are interested in [0..nftr-1]
#   $sgm_idx:       segment index we are interested in [0..nsgm-1]
#
# Returns:    Index of segment $sgm_idx [1..vdr_FeatureNumSegements($ftr_idx)]
#             -1 if segment $sgm_idx is not a segment of $ftr_idx.
# 
# Dies: Never, nothing is validated.
# 
#################################################################
sub vdr_FeatureRelativeSegmentIndex { 
  my $sub_name  = "vdr_FeatureRelativeSegmentIndex";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_AHR, $ftr_idx, $sgm_idx) = (@_);

  if(($sgm_idx >= $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}) && 
     ($sgm_idx <= $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"})) { 
    return $sgm_idx - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} + 1;
  }

  return -1;
}


#################################################################
# Subroutine: vdr_Feature5pMostPosition()
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
sub vdr_Feature5pMostPosition { 
  my $sub_name = "vdr_Feature5pMostPosition";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords, $FH_HR) = @_;
  
  if($coords =~ /^(\d+)\.\.\d+/) { 
    return $1;
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse ftr_info_HA coords string " . $coords, 1, $FH_HR); 
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine: vdr_Feature3pMostPosition()
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
sub vdr_Feature3pMostPosition { 
  my $sub_name = "vdr_Feature3pMostPosition";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords, $FH_HR) = @_;
  
  if($coords =~ /\d+\.\.(\d+)\:[\+\-]$/) { 
    return $1;
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse ftr_info_HA coords string " . $coords, 1, $FH_HR); 
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine: vdr_FeatureParentIndex
# Incept:     EPN, Tue Mar 24 16:40:27 2020
# 
# Purpose:    Return "parent_idx_str" if it is not "GBNULL"
#             else return "1"
#             be called after vdr_FeatureInfoInitializeParentIndexStrings().
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $ftr_idx:       feature index
#
# Returns:    void
# 
# Dies:       Never; does not validate $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"}
#
#################################################################
sub vdr_FeatureParentIndex {
  my $sub_name = "vdr_FeatureParentIndex";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
  
  my ($ftr_info_AHR, $ftr_idx) = @_;
  
  if($ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} eq "GBNULL") { 
    return -1;
  }
  return $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"};
}

#################################################################
# Subroutine: vdr_FeatureSummarizeSegment()
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
sub vdr_FeatureSummarizeSegment { 
  my $sub_name = "vdr_FeatureSummarizeSegment";
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
# Subroutine: vdr_FeatureStartStopStrandArrays()
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
# Returns:    number of segments
#
# Dies: if unable to parse $coords
#
#################################################################
sub vdr_FeatureStartStopStrandArrays {
  my $sub_name = "vdr_FeatureStartStopStrandArrays";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $start_AR, $stop_AR, $strand_AR, $FH_HR) = @_;
  if(! defined $coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", 1, $FH_HR); 
  }

  my @start_A  = ();
  my @stop_A   = ();
  my @strand_A = ();
  my ($start, $stop, $strand, $sgm_idx);
  my @coords_A  = split(",", $coords);
  my $nsgm = scalar(@coords_A);
  for($sgm_idx = 0; $sgm_idx < $nsgm; $sgm_idx++) { 
    ($start, $stop, $strand) = vdr_CoordsSegmentParse($coords_A[$sgm_idx], $FH_HR);
    # vdr_CoordsSegmentParse() will fail if unable to parse $coords_A[$sgm_idx]
    push(@start_A,  $start);
    push(@stop_A,   $stop);
    push(@strand_A, $strand); 
  }

  if(defined $start_AR)  { @{$start_AR}   = @start_A;  }
  if(defined $stop_AR)   { @{$stop_AR}    = @stop_A;   }
  if(defined $strand_AR) { @{$strand_AR}  = @strand_A;  }

  return $nsgm;
}

#################################################################
# Subroutine: vdr_FeatureSummaryStrand
# Incept:     EPN, Wed Mar 13 15:38:06 2019
# 
# Purpose:    Summarize the strandedness of segments for a feature
#             by parsing the "coords" value.
# 
# Arguments:
#   $coords:   coords string
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
sub vdr_FeatureSummaryStrand { 
  my $sub_name = "vdr_FeatureSummaryStrand";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($coords, $FH_HR) = @_;

  # Examples we can parse: 
  # $coords                  return value
  # -----------------------  -----------------
  # 1-200:+                  200-1:-
  # 1-200:+,300-400:+        400-300:-,200-1:-

  my @strand_A = ();
  vdr_FeatureStartStopStrandArrays($coords, undef, undef, \@strand_A, $FH_HR);

  my $npos = 0;
  my $nneg = 0;
  foreach my $strand (@strand_A) { 
    if   ($strand eq "+") { $npos++; }
    elsif($strand eq "-") { $nneg++; }
    else { ofile_FAIL("ERROR in $sub_name, unable to determine strands in coords $coords", 1, $FH_HR); }
  }

  if(($npos >  0) && ($nneg >  0)) { return "!"; }
  if(($npos >  0) && ($nneg == 0)) { return "+"; }
  if(($npos == 0) && ($nneg >  0)) { return "-"; }
  if(($npos == 0) && ($nneg == 0)) { 
    ofile_FAIL("ERROR in $sub_name, unable to determine strands in coords $coords", 1, $FH_HR); 
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine: vdr_FeaturePositionSpecificValueBreakdown()
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
sub vdr_FeaturePositionSpecificValueBreakdown { 
  my $sub_name = "vdr_FeaturePositionSpecificValueBreakdown";
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
# Subroutine: vdr_AlertInfoInitialize()
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
sub vdr_AlertInfoInitialize { 
  my $sub_name = "vdr_AlertInfoInitialize";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $FH_HR) = (@_);

  if(scalar (keys(%{$alt_info_HHR})) > 0) { 
    ofile_FAIL("ERROR in $sub_name, alert info hash of arrays already has at least one key", 1, $FH_HR);
  }

  # add each alert code, this function will die if we try to add the same code twice, or if something is wrong 
  # with how we try to add it (args to vdr_AlertInfoAdd don't pass the contract check)

  vdr_AlertInfoAdd($alt_info_HHR, "noannotn", "sequence",
                   "NO_ANNOTATION", # short description
                   "no significant similarity detected", # long  description
                   1, 1, 1, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "revcompl", "sequence",
                   "REVCOMPLEM", # short description
                   "sequence appears to be reverse complemented", # long description
                   1, 1, 1, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "qstsbgrp", "sequence",
                   "QUESTIONABLE_SPECIFIED_SUBGROUP", # short description
                   "best overall model is not from specified subgroup", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "qstgroup", "sequence",
                   "QUESTIONABLE_SPECIFIED_GROUP", # short description
                   "best overall model is not from specified group", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "incsbgrp", "sequence",
                   "INCORRECT_SPECIFIED_SUBGROUP", # short description
                   "score difference too large between best overall model and best specified subgroup model", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "incgroup", "sequence",
                   "INCORRECT_SPECIFIED_GROUP", # short description
                   "score difference too large between best overall model and best specified group model", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "lowcovrg", "sequence",
                   "LOW_COVERAGE", # short description, 
                   "low sequence fraction with significant similarity to homology model", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "indfclas", "sequence",
                   "INDEFINITE_CLASSIFICATION", # short description
                   "low score difference between best overall model and second best model (not in best model's subgroup)", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "lowscore", "sequence",
                   "LOW_SCORE", # short description
                   "score to homology model below low threshold", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "biasdseq", "sequence",
                   "BIASED_SEQUENCE", # short description
                   "high fraction of score attributed to biased sequence composition", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "dupregin", "sequence",
                   "DUPLICATE_REGIONS", # short description
                   "similarity to a model region occurs more than once", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "discontn", "sequence",
                   "DISCONTINUOUS_SIMILARITY", # short description
                   "not all hits are in the same order in the sequence and the homology model", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indfstrn", "sequence",
                   "INDEFINITE_STRAND", # short description
                   "significant similarity detected on both strands", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim5s", "sequence",
                   "LOW_SIMILARITY_START", # short description
                   "significant similarity not detected at 5' end of the sequence", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim3s", "sequence",
                   "LOW_SIMILARITY_END", # short description
                   "significant similarity not detected at 3' end of the sequence", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsimis", "sequence",
                   "LOW_SIMILARITY", # short description
                   "internal region without significant similarity", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "unexdivg", "sequence",
                   "UNEXPECTED_DIVERGENCE", # short description
                   "sequence is too divergent to confidently assign nucleotide-based annotation", # long description
                   1, 1, 1, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "noftrann", "sequence",
                   "NO_FEATURES_ANNOTATED", # short description
                   "sequence similarity to homology model does not overlap with any features", # long description
                   1, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "mutstart", "feature",
                   "MUTATION_AT_START", # short description
                   "expected start codon could not be identified", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "mutendcd", "feature",
                   "MUTATION_AT_END", # short description
                   "expected stop codon could not be identified, predicted CDS stop by homology is invalid", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "mutendns", "feature",  
                   "MUTATION_AT_END", # short description
                   "expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted valid start codon", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "mutendex", "feature",
                   "MUTATION_AT_END", # short description
                   "expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "unexleng", "feature",  
                   "UNEXPECTED_LENGTH", # short description
                   "length of complete coding (CDS or mat_peptide) feature is not a multiple of 3", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "cdsstopn", "feature",
                   "CDS_HAS_STOP_CODON", # short description
                   "in-frame stop codon exists 5' of stop position predicted by homology to reference", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "fsthicnf", "feature",
                   "POSSIBLE_FRAMESHIFT_HIGH_CONF", # short description
                   "high confidence potential frameshift in CDS", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "fstlocnf", "feature",
                   "POSSIBLE_FRAMESHIFT_LOW_CONF", # short description
                   "low confidence potential frameshift in CDS", # long description
                   0, 0, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "cdsstopp", "feature",
                   "CDS_HAS_STOP_CODON", # short description
                   "stop codon in protein-based alignment", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "peptrans", "feature",
                   "PEPTIDE_TRANSLATION_PROBLEM", # short description
                   "mat_peptide may not be translated because its parent CDS has a problem", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "pepadjcy", "feature",
                   "PEPTIDE_ADJACENCY_PROBLEM", # short description
                   "predictions of two mat_peptides expected to be adjacent are not adjacent", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indfantp", "feature",
                   "INDEFINITE_ANNOTATION", # short description
                   "protein-based search identifies CDS not identified in nucleotide-based search", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indfantn", "feature",
                   "INDEFINITE_ANNOTATION", # short description
                   "nucleotide-based search identifies CDS not identified in protein-based search", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf5gap", "feature",
                   "INDEFINITE_ANNOTATION_START", # short description
                   "alignment to homology model is a gap at 5' boundary", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf5loc", "feature",
                   "INDEFINITE_ANNOTATION_START", # short description
                   "alignment to homology model has low confidence at 5' boundary", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf5plg", "feature",
                   "INDEFINITE_ANNOTATION_START", # short description
                   "protein-based alignment extends past nucleotide-based alignment at 5' end", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf5pst", "feature",
                   "INDEFINITE_ANNOTATION_START", # short description
                   "protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf3gap", "feature",
                   "INDEFINITE_ANNOTATION_END", # short description
                   "alignment to homology model is a gap at 3' boundary", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf3loc", "feature",
                   "INDEFINITE_ANNOTATION_END", # short description
                   "alignment to homology model has low confidence at 3' boundary", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf3plg", "feature",
                   "INDEFINITE_ANNOTATION_END", # short description
                   "protein-based alignment extends past nucleotide-based alignment at 3' end", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf3pst", "feature",
                   "INDEFINITE_ANNOTATION_END", # short description
                   "protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indfstrp", "feature",
                   "INDEFINITE_STRAND", # short description
                   "strand mismatch between protein-based and nucleotide-based predictions", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "insertnp", "feature",
                   "INSERTION_OF_NT", # short description
                   "too large of an insertion in protein-based alignment", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "deletinp", "feature",
                   "DELETION_OF_NT", # short description
                   "too large of a deletion in protein-based alignment", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim5f", "feature",
                   "LOW_FEATURE_SIMILARITY_START", # short description
                   "region within annotated feature at 5' end of sequence lacks significant similarity", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim3f", "feature",
                   "LOW_FEATURE_SIMILARITY_END", # short description
                   "region within annotated feature at 3' end of sequence lacks significant similarity", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsimif", "feature",
                   "LOW_FEATURE_SIMILARITY", # short description
                   "region within annotated feature lacks significant similarity", # long description
                   0, 1, 0, # always_fails, causes_failure, prevents_annot
                   $FH_HR);

  # define the ftbl_invalid_by values, these are one-sided, any error code listed in the 
  # 3rd argument invalidates the 2nd argument error code, but not vice versa

  # cdsstopn, mutendex, mutendns are preferred to mutendcd
  vdr_AlertInfoSetFTableInvalidatedBy($alt_info_HHR, "mutendcd", "cdsstopn,mutendex,mutendns", $FH_HR); 

  # unexdivg is preferred to noftrann
  vdr_AlertInfoSetFTableInvalidatedBy($alt_info_HHR, "noftrann", "unexdivg", $FH_HR);

  # validate the alert info hash
  #vdr_AlertInfoValidate($alt_info_HHR, undef, $FH_HR); 

  return;
}

#################################################################
# Subroutine: vdr_AlertInfoAdd()
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
sub vdr_AlertInfoAdd { 
  my $sub_name = "vdr_AlertInfoAdd";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $code, $pertype, $sdesc, $ldesc, $always_fails, $causes_failure, $prevents_annot, $FH_HR) = (@_);

  # make sure $pertype is valid
  if(($pertype ne "feature") && ($pertype ne "sequence")) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code with per-type $pertype that is not neither \"feature\" nor \"sequence\".", 1, $FH_HR); 
  }
  
  # make sure $always_fails is valid
  if((! defined $causes_failure) || (($causes_failure != 0) && ($causes_failure != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but causes_failure is undefined or not 0 or 1", 1, $FH_HR);
  }

  # make sure $causes_failure is valid, and makes sense with $always_fail
  if((! defined $causes_failure) || (($causes_failure != 0) && ($causes_failure != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but causes_failure is undefined or not 0 or 1", 1, $FH_HR);
  }
  if($always_fails && (! $causes_failure)) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but always_fails is 1 and causes_failure 0", 1, $FH_HR);
  }
  
  # make sure $prevents_annot is valid
  if((! defined $prevents_annot) || (($prevents_annot != 0) && ($prevents_annot != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but prevents_annot is undefined or not 0 or 1", 1, $FH_HR);
  }

  # make sure $prevents_annot is only 1 if $pertype is "sequence"
  if(($prevents_annot == 1) && ($pertype ne "sequence")) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but prevents_annot is 1 and pertype is feature", 1, $FH_HR);
  }
  
  # check if $code already exists
  if(defined $alt_info_HHR->{$code}) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code, but it already exists in the error info hash", 1, $FH_HR);
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
# Subroutine: vdr_AlertInfoSetFTableInvalidatedBy
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
sub vdr_AlertInfoSetFTableInvalidatedBy {
  my $sub_name = "vdr_AlertInfoSetFTableInvalidatedBy";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $code1, $code2str, $FH_HR) = (@_);

  # verify the codes in $code2str
  my @code2_A = split(',', $code2str);
  foreach my $code2 (@code2_A) { 
    if(! defined $alt_info_HHR->{$code2}) { 
      ofile_FAIL("ERROR in $sub_name, trying to add invalidated by relationship between codes $code1 and $code2, but $code2 does not exist in the error info hash", 1, $FH_HR);
    }
    if($code1 eq $code2) { 
      ofile_FAIL("ERROR in $sub_name, trying to add invalidated by relationship between a code and itself: $code1 and $code2", 1, $FH_HR);
    }
  }

  # set the value
  $alt_info_HHR->{$code1}{"ftbl_invalid_by"} = $code2str;

  return;
}

#################################################################
# Subroutine: vdr_AlertInfoSetCausesFailure
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
sub vdr_AlertInfoSetCausesFailure {
  my $sub_name = "vdr_AlertInfoSetCausesFailure";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $code, $value, $FH_HR) = (@_);

  if(! defined $alt_info_HHR->{$code}) { 
    ofile_FAIL("ERROR in $sub_name, trying to set causes_failure for invalid code $code", 1, $FH_HR);
  }
  if(($value ne "1") && ($value ne "0")) { 
    ofile_FAIL("ERROR in $sub_name, trying to set causes_failure to invalid value $value (must be 1 or 0)", 1, $FH_HR);
  }

  $alt_info_HHR->{$code}{"causes_failure"} = $value;

  return;
}

#################################################################
# Subroutine:  vdr_AlertInfoDump()
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
sub vdr_AlertInfoDump { 
  my $sub_name = "vdr_AlertInfoDump";
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
# Subroutine : vdr_ParseQsubFile()
# Incept:      EPN, Mon Jul  9 10:30:41 2018 [ribovore:ribo.pm]
#
# Purpose:     Parse a file that specifies the qsub command to use
#              when submitting jobs to the farm. The file should 
#              have exactly 2 non-'#' prefixed lines. Chomp each
#              and return them.
#              
# Arguments: 
#   $qsub_file:  file to parse
#   $FH_HR:      REF to hash of file handles
#
# Returns:     2 values:
#              $qsub_prefix: string that is the qsub command prior to the 
#                            actual cmsearch/cmalign command
#              $qsub_suffix: string that is the qsub command after the 
#                            actual cmsearch/cmalign command
# 
# Dies:        If we can't parse the qsub file because it is not
#              in the correct format.
#
################################################################# 
sub vdr_ParseQsubFile { 
  my $nargs_expected = 2;
  my $sub_name = "vdr_ParseQsubFile";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($qsub_file, $FH_HR) = @_;

  open(IN, $qsub_file) || ofile_FileOpenFailure($qsub_file, $sub_name, $!, "reading", $FH_HR);

  my $qsub_prefix_line = undef;
  my $qsub_suffix_line = undef;
  while(my $line = <IN>) { 
    if($line !~ m/^\#/) { 
      chomp $line;
      if   (! defined $qsub_prefix_line) { $qsub_prefix_line = $line; }
      elsif(! defined $qsub_suffix_line) { $qsub_suffix_line = $line; }
      else { # both $qsub_prefix_line and $qsub_suffix_line are defined, this shouldn't happen
        ofile_FAIL("ERROR in $sub_name, read more than 2 non-# prefixed lines in file $qsub_file:\n$line\n", 1, $FH_HR);
      }
    }
  }
  close(IN);
  
  if(! defined $qsub_prefix_line) { 
    ofile_FAIL("ERROR in $sub_name, read zero non-# prefixed lines in file $qsub_file, but expected 2", 1, $FH_HR);
  }
  if(! defined $qsub_suffix_line) { 
    ofile_FAIL("ERROR in $sub_name, read only one non-# prefixed lines in file $qsub_file, but expected 2", 1, $FH_HR);
  }

  return($qsub_prefix_line, $qsub_suffix_line);
}

#################################################################
# Subroutine: vdr_SubmitJob()
# Incept:      EPN, Wed Feb  6 12:35:04 2019
#
# Purpose:     Submits a job to sge.
#
# Arguments:
#   $cmd:            command to run
#   $qsub_prefix:    qsub command prefix to use when submitting to farm, undef if running locally
#   $qsub_suffix:    qsub command suffix to use when submitting to farm, undef if running locally
#   $job_name:       name for job
#   $err_file:       name of err file to create, can be "/dev/null"
#   $mem_gb:         number of Gb of memory required
#   $nsecs:          maximum number of seconds to allow jobs to take
#   $opt_HHR:        REF to 2D hash of option values, see top of sqp_opts.pm for description, PRE-FILLED
#   $ofile_info_HHR: REF to the 2D hash of output file information, ADDED TO HERE 
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if qsub vdr_$cmd fails
#################################################################
sub vdr_SubmitJob {
  my $sub_name = "vdr_SubmitJob()";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $qsub_prefix, $qsub_suffix, $job_name, $err_file, $mem_gb, $nsecs, $opt_HHR, $ofile_info_HHR) = @_;
  
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  if(($err_file ne "/dev/null") && (-e $err_file)) { 
    utl_FileRemoveUsingSystemRm($err_file, $sub_name, $opt_HHR, $ofile_info_HHR); 
  }
  my $submit_cmd = $qsub_prefix . $cmd . $qsub_suffix;
  # replace changeable parts of qsub suffix and prefix
  $submit_cmd =~ s/\!\[jobname\]\!/$job_name/g;
  $submit_cmd =~ s/\!\[errfile\]\!/$err_file/g;
  $submit_cmd =~ s/\!\[memgb\]\!/$mem_gb/g;
  $submit_cmd =~ s/\!\[nsecs\]\!/$nsecs/g;

  utl_RunCommand($submit_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine:  vdr_WaitForFarmJobsToFinish()
# Incept:      EPN, Mon Feb 29 16:20:54 2016
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
sub vdr_WaitForFarmJobsToFinish { 
  my $sub_name = "vdr_WaitForFarmJobsToFinish()";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($do_cmalign, $out_file_AHR, $success_AR, $mxsize_AR, $finished_str, $opt_HHR, $FH_HR) = @_;

  my $log_FH = $FH_HR->{"log"};
  my $nmin = opt_Get("--wait", $opt_HHR);
  my $do_errcheck = opt_Get("--errcheck", $opt_HHR);

  # contract check
  if(($do_cmalign) && (! exists $out_file_AHR->[0]{"stdout"})) { 
    ofile_FAIL("ERROR in $sub_name, cmalign mode, no stdout files in out_file_AHR", 1, $FH_HR);
  }
  if((! $do_cmalign) && (! exists $out_file_AHR->[0]{"tblout"})) { 
    ofile_FAIL("ERROR in $sub_name, cmsearch mode, no stdout files in out_file_AHR", 1, $FH_HR);
  }
  if(! exists $out_file_AHR->[0]{"err"}) { 
    ofile_FAIL("ERROR in $sub_name, no err files in out_file_AHR", 1, $FH_HR);
  }

  my $outkey = ($do_cmalign) ? "stdout" : "tblout";
  my @outfile_A = ();
  my @errfile_A = ();
  utl_ArrayOfHashesToArray($out_file_AHR, \@outfile_A, $outkey);
  utl_ArrayOfHashesToArray($out_file_AHR, \@errfile_A, "err");

  my $njobs = scalar(@outfile_A);
  if($njobs != scalar(@errfile_A)) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, number of elements in outfile array ($njobs) differ from number of jobs in errfile array (%d)", scalar(@errfile_A)), 1, $FH_HR);
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
            my $success = vdr_CmalignCheckStdOutput($outfile_A[$i], 
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
    ofile_FAIL($errmsg, 1, $FH_HR);
  }

  # if we get here we have no failures
  return $nfinished;
}

#################################################################
# Subroutine: vdr_CoordsSegmentParse()
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
sub vdr_CoordsSegmentParse {
  my $sub_name = "vdr_CoordsSegmentParse";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords_tok, $FH_HR) = @_;
  if(! defined $coords_tok) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", 1, $FH_HR); 
  }
  if($coords_tok =~ /^\<?(\d+)\.\.\>?(\d+)\:([\+\-])$/) { 
    return ($1, $2, $3);
  }
  ofile_FAIL("ERROR in $sub_name, unable to parse coords token $coords_tok", 1, $FH_HR); 

  return; # NEVER REACHED
}

#################################################################
# Subroutine: vdr_CoordsSegmentCreate()
# Incept:     EPN, Mon Apr 29 14:07:26 2019
#
# Synopsis: Create a coords token from a given start, stop, strand
# 
# Arguments:
#  $start:    start position
#  $stop:     stop position
#  $strand:   strand ("+" or "-")
#  $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    coordinate token <start>..<stop>:<strand>
#
# Dies:  if $start or $stop is invalid
#        if $strand is not "+" or "-"
#
#################################################################
sub vdr_CoordsSegmentCreate {
  my $sub_name = "vdr_CoordsSegmentCreate";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($start, $stop, $strand, $FH_HR) = @_;
  if($start !~ /^\<?(\d+)$/) { ofile_FAIL("ERROR in $sub_name, start is invalid ($start)", 1, $FH_HR); }
  if($stop  !~ /^\>?(\d+)$/) { ofile_FAIL("ERROR in $sub_name, stop is invalid ($stop)", 1, $FH_HR); }
  if(($strand ne "+") && ($strand ne "-")) { ofile_FAIL("ERROR in $sub_name, strand is invalid ($strand)", 1, $FH_HR); }

  return $start . ".." . $stop . ":" . $strand;
}

#################################################################
# Subroutine: vdr_CoordsAppendSegment()
#
# Incept:     EPN, Fri Mar 20 09:11:03 2020
#
# Synopsis: Append a segment $coords_sgm to $coords and return 
#           the result, after potentially adding a ',' before
#           $coords_sgm. If $coords is undef, return $coords_sgm.
#
# Arguments:
#  $coords:     existing coords string, can be undef or ""
#  $coords_sgm: segment to append, not checked for validity
#
# Returns:   Coords string with $coords_sgm appended.
#
# Dies: never
#
#################################################################
sub vdr_CoordsAppendSegment { 
  my $sub_name = "vdr_CoordsAppendSegment";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $coords_sgm) = @_;

  if((! defined $coords) || ($coords eq "")) { 
    return $coords_sgm;
  }
  return $coords . "," . $coords_sgm;
}

#################################################################
# Subroutine: vdr_CoordsLength()
# Incept:     EPN, Tue Mar 26 05:56:08 2019
#
# Synopsis: Given a comma separated coords string, parse it, 
#           validate it, and return its length.
# 
# Arguments:
#  $coords:  coordinate string
#  $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:   total length of segments in $coords
#
# Dies: if unable to parse $coords
#
#################################################################
sub vdr_CoordsLength {
  my $sub_name = "vdr_CoordsLength";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $FH_HR) = @_;
  if(! defined $coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", 1, $FH_HR); 
  }

  # if there's no comma, we should have a single span
  if($coords !~ m/\,/) { 
    my ($start, $stop, undef) = vdr_CoordsSegmentParse($coords, $FH_HR);
    return abs($start - $stop) + 1;
  }
  # else, split it up and sum length of all
  my @coords_A  = split(",", $coords);
  my ($start, $stop);
  my $ret_len = 0;
  foreach my $coords_tok (@coords_A) { 
    ($start, $stop, undef) = vdr_CoordsSegmentParse($coords_tok, $FH_HR);
    $ret_len += abs($start - $stop) + 1;
  }

  return $ret_len;
}

#################################################################
# Subroutine: vdr_CoordsFromLocation
# Incept:     EPN, Wed Mar 13 14:17:08 2019
# 
# Purpose:    Convert a GenBank file 'location' value to 
#             a coords string in the format:
#             <start1>-<stop2>:<strand1>,<start2>-<stop2>:<strand2>,...,<startN>-<stopN>:<strandN>
# 
#             <startN>: may begin with "<" carrot.
#             <stopN>: may begin with ">" carrot.
#
#             If $do_carrots is 1: include them in the
#             appropriate position in the returned
#             $coords string
#             Else: remove carrots
#
#             This function has to call itself recursively in some
#             cases.
# 
# Examples:
#                                      $do_carrots=0        $do_carrots=1
# $location                            return value         return value
# ---------------------------------    -----------------    ----------------
# 1..200                               1..200:+             1..200:+
# <1..200                              1..200:+             <1..200:+
# 100..>200                            100..200:+           100..>200:+
# <1..>200                             1..200:+             <1..>200:+
# complement(1..200)                   200..1:-             200..1:-           
# join(1..200,300..400)                1..200:+,300..400:+  1..200:+,300..400:+  
# complement(join(1..200,300..400))    400..300:-,200..1:-  400..300:-,200..1:-
# join(1..200,complement(<300..400))   1..200:+,400..300:-  1..200:+,400..>300:-
# join(complement(300..>400),<1..>200) 400..300:-,1..200:+  <400..300:-,<1..>200:+
#
# See t/01-coords.t for additional examples
# 
# Arguments:
#   $location:   GenBank file location string
#   $do_carrots: '1' to have return values include carrots, '0' not to
#   $FH_HR:      REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if unable to parse $location
#
# Ref: GenBank release notes (release 230.0) as of this writing
#      and
#      https://www.ncbi.nlm.nih.gov/genbank/samplerecord/
#################################################################
sub vdr_CoordsFromLocation { 
  my $sub_name = "vdr_CoordsFromLocation";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($location, $do_carrots, $FH_HR) = @_;

  my $ret_val = "";
  if($location =~ /^join\((.+)\)$/) { 
    my $location_to_join = $1;
    $ret_val = vdr_CoordsFromLocation($location_to_join, $do_carrots, $FH_HR);
  }
  elsif($location =~ /^complement\((.+)\)$/) { 
    my $location_to_complement = $1;
    my $coords_to_complement = vdr_CoordsFromLocation($location_to_complement, $do_carrots, $FH_HR);
    $ret_val = vdr_CoordsReverseComplement($coords_to_complement, $do_carrots, $FH_HR);
  }
  elsif($location =~ /\,/) { 
    # not wrapped in join() or complement(), but multiple segments
    foreach my $location_el (split(",", $location)) { 
      $ret_val = vdr_CoordsAppendSegment($ret_val, vdr_CoordsFromLocation($location_el, $do_carrots, $FH_HR));
    }
  }
  elsif($do_carrots) { 
    # only difference with this block and following block ($do_carrot == 0)
    # is that the carrots are possibly included in the $ret_val in this block
    if($location =~ /^(\<?\d+\.\.\>?\d+)$/) { 
      $ret_val = $1 . ":+"; # a recursive call due to the complement() may complement this
    }
    elsif($location =~ /^(\<?\d+)$/) { # single nucleotide
      $ret_val = $1 . ".." . $1 . ":+"; # a recursive call due to the complement() may complement this
    }
    elsif($location =~ /^(\>?\d+)$/) { # single nucleotide
      $ret_val = $1 . ".." . $1 . ":+"; # a recursive call due to the complement() may complement this
    }
    else { 
      ofile_FAIL("ERROR in $sub_name, unable to parse location token $location", 1, $FH_HR);
    }
  }
  else { # $do_carrot is 0
    if($location =~ /^\<?(\d+)\.\.\>?(\d+)$/) { 
      $ret_val = $1 . ".." . $2 . ":+"; # a recursive call due to the complement() may complement this
    }
    elsif($location =~ /^\<?(\d+)$/) { # single nucleotide
      $ret_val = $1 . ".." . $1 . ":+"; # a recursive call due to the complement() may complement this
    }
    elsif($location =~ /^\>?(\d+)$/) { # single nucleotide
      $ret_val = $1 . ".." . $1 . ":+"; # a recursive call due to the complement() may complement this
    }
    else { 
      ofile_FAIL("ERROR in $sub_name, unable to parse location token $location", 1, $FH_HR);
    }
  }

  return $ret_val;
}

#################################################################
# Subroutine: vdr_CoordsReverseComplement
# Incept:     EPN, Wed Mar 13 15:00:24 2019
# 
# Purpose:    Reverse complement a coords string by reverse complementing
#             all tokens within it (by calling 
#             vdr_CoordsSegmentReverseComplement()), then reversing 
#             the order of the resulting tokens, concatenating them
#             and returning that concatenated string.
#
#             If $do_carrots is 1: include "<" before start and
#             ">" before stop if they exist. 
#             Else: remove carrots
# 
# 
# Examples:
#                          $do_carrots=0         $do_carrots=1
# $coords                  return value          return value
# ---------------------    -----------------     ----------------
# 1..200:+                 200..1:-              200..1:-
# <1..200:+                200..1:-              200..>1:-
# 100..>200:+              200..100:-            <200..100:-
# <1..>200                 200..1:-              <200..>1:-
# 200..1:-                 1..200:+              1..200:+
# 1..200:+,300..400:+      400..300:-,200..1:-   400..300:-,200..1:-
# 400..300:-,200..1:-      1..200:+,300..400:+   1..200:+,300..400:+
# 1..200:+,400..>300:-     300..400:+,200..1:-   <300..400:+,200..1:-
# <400..300:-,1..200:+     200..1:-,300.. 400:+  200..1:-,300..>400:+
# 
# See t/01-coords.t for additional examples
#
# Arguments:
#   $coords:     coords string to complement
#   $do_carrots: '1' to have return values include carrots, '0' not to
#   $FH_HR:      REF to hash of file handles, including "log" and "cmd"
#
# Returns:    reverse complemented $coords
# 
# Dies:       if unable to parse $coords
#
#################################################################
sub vdr_CoordsReverseComplement { 
  my $sub_name = "vdr_CoordsReverseComplement";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($coords, $do_carrots, $FH_HR) = @_;

  my @tok_A = split(",", $coords);
  my $ntok = scalar(@tok_A);
  my @ret_coords_tok_A = (); 

  for(my $i = 0; $i < $ntok; $i++) { 
    push(@ret_coords_tok_A, vdr_CoordsSegmentReverseComplement($tok_A[$i], $do_carrots, $FH_HR));
  }
  
  # concatenate the tokens in reverse order
  my $ret_val = "";
  for(my $i = $ntok-1; $i >= 0; $i--) { 
    $ret_val = vdr_CoordsAppendSegment($ret_val, $ret_coords_tok_A[$i]);
  }
  # printf("\tin $sub_name, coords: $coords ret_val: $ret_val\n");
      
  return $ret_val;
}

#################################################################
# Subroutine: vdr_CoordsSegmentReverseComplement
# Incept:     EPN, Thu Mar 19 15:29:22 2020
# 
# Purpose:    Reverse complement a single coords tokenn by 
#             reverse complementing it. 
#
#             If $do_carrots is 1: include "<" before start and
#             ">" before stop if they exist. 
#             Else: remove carrots
# 
# 
# Examples:
#                  $do_carrots=0         $do_carrots=1
# $coords          return value          return value
# -------------    -----------------     ----------------
# 1..200:+         200..1:-              200..1:-
# <1..200:+        200..1:-              200..>1:-
# 100..>200:+      200..100:-            <200..100:-
# <1..>200         200..1:-              <200..>1:-
# 200..1:-         1..200:+              1..200:+
# 
# See t/01-coords.t for additional examples
#
# Arguments:
#   $coords_tok: coords string to complement
#   $do_carrots: '1' to have return values include carrots, '0' not to
#   $FH_HR:      REF to hash of file handles, including "log" and "cmd"
#
# Returns:    reverse complemented $coords_tok
# 
# Dies:       if we are unable to parse $coords_tok (including if it 
#             has a ',' in it implying it is multiple segments)
#
#################################################################
sub vdr_CoordsSegmentReverseComplement { 
  my $sub_name = "vdr_CoordsSegmentReverseComplement";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($coords_tok, $do_carrots, $FH_HR) = @_;

  if($coords_tok =~ /^(\<?)(\d+)\.\.(\>?)(\d+)\:([\+\-])$/) { 
    my ($start_carrot, $start, $stop_carrot, $stop, $strand) = ($1, $2, $3, $4, $5);
    if($do_carrots) { 
      if($start_carrot ne "") { $start_carrot = ">"; } # swap it from < to >
      if($stop_carrot  ne "") { $stop_carrot  = "<"; } # swap it from > to <
    }
    else { 
      $start_carrot = "";
      $stop_carrot  = "";
    }    
    my $ret_strand = "?";
    if   ($strand eq "+") { $ret_strand = "-"; }
    elsif($strand eq "-") { $ret_strand = "+"; }

    return $stop_carrot . $stop . ".." . $start_carrot . $start . ":" . $ret_strand;
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse input coords token: $coords_tok", 1, $FH_HR); 
  }
  return ""; # NEVER REACHED
}


#################################################################
# Subroutine: vdr_CoordsMin()
# Incept:     EPN, Mon Apr 29 13:49:55 2019
#
# Synopsis: Given a comma separated coords string, return the 
#           minimum position that it corresponds to.
# 
# Arguments:
#  $coords:  coordinate string
#  $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $coords
#
#################################################################
sub vdr_CoordsMin {
  my $sub_name = "vdr_CoordsMin";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $FH_HR) = @_;
  if(! defined $coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", 1, $FH_HR); 
  }

  # if there's no comma, we should have a single span
  if($coords !~ m/\,/) { 
    my ($start, $stop, undef) = vdr_CoordsSegmentParse($coords, $FH_HR);
    return utl_Min($start, $stop);
  }
  # else, split it up and find minimum
  my @coords_A  = split(",", $coords);
  my ($start, $stop);
  my $ret_min = undef;
  foreach my $coords_tok (@coords_A) { 
    ($start, $stop, undef) = vdr_CoordsSegmentParse($coords_tok, $FH_HR);
    if(! defined $ret_min) { 
      $ret_min = utl_Min($start, $stop);
    }
    else { 
      $ret_min = utl_Min($ret_min, utl_Min($start, $stop));
    }
  }

  return $ret_min;
}

#################################################################
# Subroutine: vdr_CoordsMax()
# Incept:     EPN, Mon Apr 29 13:57:57 2019
#
# Synopsis: Given a comma separated coords string, return the 
#           maximum position that it corresponds to.
# 
# Arguments:
#  $coords:  coordinate string
#  $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $coords
#
#################################################################
sub vdr_CoordsMax {
  my $sub_name = "vdr_CoordsMax";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $FH_HR) = @_;
  if(! defined $coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", 1, $FH_HR); 
  }

  # if there's no comma, we should have a single span
  if($coords !~ m/\,/) { 
    my ($start, $stop, undef) = vdr_CoordsSegmentParse($coords, $FH_HR);
    return utl_Max($start, $stop);
  }
  # else, split it up and find maximum
  my @coords_A  = split(",", $coords);
  my ($start, $stop);
  my $ret_max = undef;
  foreach my $coords_tok (@coords_A) { 
    ($start, $stop, undef) = vdr_CoordsSegmentParse($coords_tok, $FH_HR);
    if(! defined $ret_max) { 
      $ret_max = utl_Max($start, $stop);
    }
    else { 
      $ret_max = utl_Max($ret_max, utl_Max($start, $stop));
    }
  }

  return $ret_max;
}

#################################################################
# Subroutine: vdr_CoordsMissing()
# Incept:     EPN, Mon Apr 29 13:57:57 2019
#
# Synopsis: Given a comma separated coords string, a strand
#           ("+" or "-") and the total length, return a 
#           comma separated coords string with each interval
#           on strand $in_strand that is missing in $in_coords.
# 
# Arguments:
#  $in_coords: coordinate string
#  $in_strand: strand we are interested in
#  $in_length: length of sequence
#  $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $in_coords
#       if $in_coords has a position that is < 0 or exceeds $in_length
#################################################################
sub vdr_CoordsMissing {
  my $sub_name = "vdr_CoordsMissing";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($in_coords, $in_strand, $in_length, $FH_HR) = @_;

  if(! defined $in_coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", 1, $FH_HR); 
  }
  if(! defined $in_strand) { 
    ofile_FAIL("ERROR in $sub_name, strand is undefined", 1, $FH_HR); 
  }
  if(! defined $in_length) { 
    ofile_FAIL("ERROR in $sub_name, length is undefined", 1, $FH_HR); 
  }

  my @coords_A  = split(",", $in_coords); # the tokens in $in_coords
  my $ret_coords = ""; # return coordinates value
  my @covered_A = ();  # 1..$i..$in_length: '1' if a coordinate token in $in_coords covers position $i on strand $in_strand, else '0'
  my $i;               # sequence position 1..$in_length
  my ($start, $stop, $strand); # start, stop and strand of a coords token

  # initialize
  for($i = 0; $i <= $in_length; $i++) { 
    $covered_A[$i] = 0;
  }

  # fill @covered_A based on @coords_A
  foreach my $coords_tok (@coords_A) { 
    ($start, $stop, $strand) = vdr_CoordsSegmentParse($coords_tok, $FH_HR);
    if(($start < 0) || ($start > $in_length)) { 
      ofile_FAIL("ERROR in $sub_name, start is invalid ($start in_length: $in_length)", 1, $FH_HR); 
    }
    if(($stop < 0) || ($stop > $in_length)) { 
      ofile_FAIL("ERROR in $sub_name, stop is invalid ($stop, in_length: $in_length)", 1, $FH_HR); 
    }
    if($strand eq $in_strand) { 
      my $min = utl_Min($start, $stop);
      my $max = utl_Max($start, $stop);
      for($i = $min; $i <= $max; $i++) { $covered_A[$i] = 1; }
    }
  }
  
  # go back and create return coords
  for($i = 1; $i <= $in_length; $i++) { 
    if($covered_A[$i] == 0) { 
      $start = $i;
      while((($i+1) <= $in_length) && ($covered_A[($i+1)] == 0)) { $i++; }
      $stop = $i; 
      $ret_coords = vdr_CoordsAppendSegment($ret_coords, vdr_CoordsSegmentCreate($start, $stop, $in_strand, $FH_HR));
    }
  }

  return $ret_coords;
}

#################################################################
# Subroutine: vdr_CoordsCheckIfSpans()
# Incept:     EPN, Tue May 21 09:59:24 2019
#
# Synopsis: Check if coords value $coords1 completely spans
#           $coords2. This occurs if every segment in $coords2
#           is spanned completely by at least 1 segment in 
#           $coords2.
# 
# Arguments:
#  $coords1: coordinate string 1
#  $coords2: coordinate string 2
#  $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:   '1' if $coords1 completely spans $coords2
#            '0' if not
#
# Dies: if unable to parse $coords
#
#################################################################
sub vdr_CoordsCheckIfSpans {
  my $sub_name = "vdr_CoordsCheckIfSpans";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords1, $coords2, $FH_HR) = @_;
  if(! defined $coords1) { ofile_FAIL("ERROR in $sub_name, coords1 is undefined", 1, $FH_HR); }
  if(! defined $coords2) { ofile_FAIL("ERROR in $sub_name, coords2 is undefined", 1, $FH_HR); }

  my @coords1_A = split(",", $coords1);
  my @coords2_A = split(",", $coords2);

  foreach my $coords2_tok (@coords2_A) { 
    my $found_overlap = 0;
    my $coords2_tok_len = vdr_CoordsLength($coords2_tok, $FH_HR);
    foreach my $coords1_tok (@coords1_A) { 
      if(! $found_overlap) { 
        my ($noverlap, undef) = vdr_CoordsSegmentOverlap($coords2_tok, $coords1_tok, $FH_HR);
        if($noverlap == $coords2_tok_len) { 
          $found_overlap = 1;
        }
      }
    }
    if(! $found_overlap) { # no token in $coords2 completely contains $coords1_tok
      # printf("in $sub_name, coords1: $coords1 coords2: $coords2, returning 0\n");
      return 0;
    }
  }

  # if we get here, there exists a token in $coords2 that completely 
  # spans each token in $coords1 
  # printf("in $sub_name, coords1: $coords1 coords2: $coords2, returning 1\n");
  return 1;
}

#################################################################
# Subroutine: vdr_CoordsSegmentOverlap()
# Incept:     EPN, Tue May 21 10:06:26 2019
#
# Synopsis: Return number of positions of overlap between 
#           <$coords_tok1> and <$coords_tok2> on the same 
#           strand.
# 
# Arguments:
#  $coords_tok1: coordinate token 1
#  $coords_tok2: coordinate token 2
#  $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:   Number of positions of overap between <$coords1_tok>
#            and <$coords2_tok> on the same strand.
#            '0' if no overlap or the two tokens are on opposite strands.
#
# Dies: if unable to parse $coords_tok1 or $coords_tok2
#
#################################################################
sub vdr_CoordsSegmentOverlap { 
  my $sub_name = "vdr_CoordsSegmentOverlap";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords_tok1, $coords_tok2, $FH_HR) = @_;
  if(! defined $coords_tok1) { ofile_FAIL("ERROR in $sub_name, coords_tok1 is undefined", 1, $FH_HR); }
  if(! defined $coords_tok2) { ofile_FAIL("ERROR in $sub_name, coords_tok2 is undefined", 1, $FH_HR); }

  my ($start1, $stop1, $strand1) = vdr_CoordsSegmentParse($coords_tok1, $FH_HR);
  my ($start2, $stop2, $strand2) = vdr_CoordsSegmentParse($coords_tok2, $FH_HR);

  if($strand1 ne $strand2) { # strand mismatch
    return 0;
  }

  if($strand1 eq "-") { # $strand2 must be "-" too
    utl_Swap(\$start1, \$stop1); 
    utl_Swap(\$start2, \$stop2);
  }

  return seq_Overlap($start1, $stop1, $start2, $stop2, $FH_HR);
}

#################################################################
# Subroutine: vdr_CoordsRelativeToAbsolute()
#             formerly vdr_CoordsProtToNuc() (pre v1.0.5)
#
# Incept:     EPN, Wed May 22 09:49:30 2019
#
# Synopsis: Return absolute nucleotide coordinates that correspond to
#           the relative nucleotide coordinates in <$rel_coords>.
#           Work is done by calling vdr_CoordsRelativeSegmentToAbsolute()
#           for each segment in <$rel_coords>, concatenating all the 
#           returned coords, and then condensing them to combine any
#           adjacent segments.
#
# Arguments:
#  $abs_coords:  nucleotide coordinates in full sequence [1..seqlen]
#  $rel_coords:  relative nucleotide coordinates within $abs_coords
#  $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:   Absolute coordinates coords string corresponding to $rel_coords.
#
# Dies: if $rel_coords has a position that is longer than
#       total length of absolute coords 
#       (dies within vdr_CoordsRelativeSegmentToAbsolute())
#
#################################################################
sub vdr_CoordsRelativeToAbsolute { 
  my $sub_name = "vdr_CoordsRelativeToAbsolute";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($abs_coords, $rel_coords, $FH_HR) = @_;

  # check to make sure we have identical strands for all segments in both
  # $abs_coords and $rel_coords
  # the code may work even if we don't but I didn't add tests for these
  # mixed strand cases, so I'm disallowing them here.
  # If you want to allow mixed strand in the future, add tests in
  # 01-coords.t (see note in that file as well dated Fri Mar 20, 2020)
  my $abs_summary_strand = vdr_FeatureSummaryStrand($abs_coords, $FH_HR); # will be + if 
  my $rel_summary_strand = vdr_FeatureSummaryStrand($rel_coords, $FH_HR);
  # $*_summary_strand values will be + if all segments are +, - if all segments are -, else !
  if($abs_summary_strand eq "!") { 
    ofile_FAIL("ERROR in $sub_name, in abs_coords $abs_coords not all segments are the same strand", 1, $FH_HR);
  }
  if($rel_summary_strand eq "!") { 
    ofile_FAIL("ERROR in $sub_name, in rel_coords $rel_coords not all segments are the same strand", 1, $FH_HR);
  }

  my @rel_coords_sgm_A = split(",", $rel_coords);
  my $ret_coords = "";
  foreach my $rel_coords_sgm (@rel_coords_sgm_A) { 
    $ret_coords = vdr_CoordsAppendSegment($ret_coords, vdr_CoordsRelativeSegmentToAbsolute($abs_coords, $rel_coords_sgm, $FH_HR));
  }
  $ret_coords = vdr_CoordsMergeAllAdjacentSegments($ret_coords, $FH_HR);
  
  return $ret_coords;
}

#################################################################
# Subroutine: vdr_CoordsRelativeSegmentToAbsolute()
#             formerly vdr_CoordsProtToNuc() (pre v1.0.5)
#
# Incept:     EPN, Wed May 22 09:49:30 2019
#
# Synopsis: Return absolute nucleotide coordinates that correspond to
#           the relative nucleotide coordinates segment in <$rel_coords_tok>.
#           with nucleotide sequence with absolute coords <$abs_coords>.
#
#           Examples:
# 
#           abs_coords     rel_coords  returns
#           "11..100:+"    "6..38:+"   "16..48:+"     
#
# Arguments:
#  $abs_coords:     nucleotide coordinates in full sequence [1..seqlen]
#  $rel_coords_sgm: relative nucleotide coordinates (single segment) within $abs_coords
#  $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:   Absolute coordinates coords string corresponding to $rel_coords_sgm.
#
# Dies: if $rel_coords_sgm has a position that is longer than
#       total length of absolute coords
#
#################################################################
sub vdr_CoordsRelativeSegmentToAbsolute { 
  my $sub_name = "vdr_CoordsRelativeSegmentToAbsolute";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($abs_coords, $rel_coords_sgm, $FH_HR) = @_;

  # printf("in $sub_name, abs_coords: $abs_coords, rel_coords_sgm: $rel_coords_sgm\n");

  my $ret_coords = ""; # return value 

  # $rel_coords_sgm should not be multiple segments (i.e have any ',' characters)
  if($rel_coords_sgm =~ m/\,/) { 
    ofile_FAIL("ERROR in $sub_name, rel_coords_sgm: $rel_coords_sgm seems to be multiple segments (has at least one ',')", 1, $FH_HR); 
  }

  # breakdown the absolute and relative coords
  my $abs_len = vdr_CoordsLength($abs_coords, $FH_HR);
  my $rel_len = vdr_CoordsLength($rel_coords_sgm, $FH_HR);
  my ($orig_rel_start, $orig_rel_stop, $orig_rel_strand) = vdr_CoordsSegmentParse($rel_coords_sgm, $FH_HR);

  my @abs_start_A  = (); # array of starts  for $abs_coords segments
  my @abs_stop_A   = (); # array of stops   for $abs_coords segments
  my @abs_strand_A = (); # array of strands for $abs_coords segments
  my $nsgm_abs = vdr_FeatureStartStopStrandArrays($abs_coords, \@abs_start_A, \@abs_stop_A, \@abs_strand_A, $FH_HR);

  # if $orig_rel_strand is -, reverse complement the relative coords, we will reverse complement comp $ret_coords back at the end prior to return
  my ($rel_start, $rel_stop, $rel_strand); 
  if   ($orig_rel_strand eq "+") { ($rel_start, $rel_stop, $rel_strand) = ($orig_rel_start, $orig_rel_stop, "+"); }
  elsif($orig_rel_strand eq "-") { ($rel_start, $rel_stop, $rel_strand) = ($orig_rel_stop, $orig_rel_start, "-"); }
  else                           { ofile_FAIL("ERROR in $sub_name, relative coords segment $rel_coords_sgm is not + or - strand", 1, $FH_HR); }

  # make sure our relative start and stop are <= $abs_len
  if($rel_start > $abs_len) { ofile_FAIL("ERROR in $sub_name, relative coords segment has a position $rel_start that exceeds absolute coords ($abs_coords) length: $abs_len", 1, $FH_HR); }
  if($rel_stop  > $abs_len) { ofile_FAIL("ERROR in $sub_name, relative coords segment has a position $rel_stop  that exceeds absolute coords ($abs_coords) length: $abs_len", 1, $FH_HR); }
  
  # do the conversion, one absolute segment at a time
  my $cur_len_to_convert = abs($rel_stop - $rel_start) + 1; # number of positions left to convert from relative to absolute coords
  my $cur_abs_offset = $rel_start - 1; # number of nucleotides to skip from start of $abs_coords_tok when converting coords
  my $conv_start = undef; # a start position converted from relative to absolute coords
  my $conv_stop  = undef; # a stop  position converted from relative to absolute coords
  for(my $a = 0; $a < $nsgm_abs; $a++) { 
    if($cur_len_to_convert > 0) { # only need to look at this token if we still have sequence left to convert
      my ($abs_start, $abs_stop, $abs_strand) = ($abs_start_A[$a], $abs_stop_A[$a], $abs_strand_A[$a]);
      my $abs_sgm_len = abs($abs_start - $abs_stop) + 1;
      if($cur_abs_offset < $abs_sgm_len) { 
        # this abs token has >= 1 nt corresponding to $rel_coords_sgm
        if($abs_strand eq "+") { 
          $conv_start = $abs_start + $cur_abs_offset;
          $conv_stop  = $conv_start + $cur_len_to_convert - 1;
          # make sure we didn't go off the end of the abs token
          if($conv_stop > $abs_stop) { $conv_stop = $abs_stop; }
        }
        elsif($abs_strand eq "-") { 
          $conv_start = $abs_start - $cur_abs_offset;
          $conv_stop  = $conv_start - $cur_len_to_convert + 1;
          # make sure we didn't go off the end of the nt token
          if($conv_stop < $abs_stop) { $conv_stop = $abs_stop; }
        }          
        else { 
          ofile_FAIL("ERROR in $sub_name, abs_coords token $abs_start..$abs_stop:$abs_strand has strand $abs_strand that is not either + or -", 1, $FH_HR); 
        }
        # update length of $rel_coords_sgm we still have to convert
        $cur_len_to_convert -= abs($conv_stop - $conv_start) + 1;
        
        # append converted token to return coords string
        if($ret_coords ne "") { $ret_coords .= ","; }
        $ret_coords .= vdr_CoordsSegmentCreate($conv_start, $conv_stop, $abs_strand, $FH_HR);
      }
      # adjust nt offset so that for next nt coords token we start at the proper position
      $cur_abs_offset -= $abs_sgm_len;
      if($cur_abs_offset < 0) { $cur_abs_offset = 0; }
    } # end of 'if($cur_len_to_convert > 0)
  } # end of 'for(my $a = 0; $a < $nsgm_abs; $a++)'

  if($orig_rel_strand eq "-") { 
    $ret_coords = vdr_CoordsReverseComplement($ret_coords, 0, $FH_HR); # 0: do not include carrots
  }

  # printf("in $sub_name, returning $ret_coords\n");
  return $ret_coords;
}

#################################################################
# Subroutine: vdr_CoordsProteinRelativeToAbsoluate()
#
# Incept:     EPN, Fri Mar 20 07:12:19 2020
#
# Synopsis: Return absolute nucleotide coordinates that correspond to
#           the relative protein coordinates in the coords string
#           <$rel_coords> for a protein encoded by the nucleotide
#           sequence with absolute coordinates <$abs_coords>
#  
#           abs_nt_coords          rel_pt_coords     returns
#           "11..100:+"            "2..11:+"         "14..43:+"     
#           "100..11:-"            "2..11:+"         "97..68:-"
#           "11..40:+,42..101:+"   "2..11:+"         "14..40:+,42..44:+"
#           "11..100:+"            "2..3:+,5..11:+"  "14..19:+,23..43:+"
#           "100..11:-"            "2..3:+,5..11:+"  "97..92:-,88..68:-"
#           "11..40:+,42..101:+"   "2..3:+,5..11:+"  "14..19:+,23..40:+,42..44:+"
#
# See t/01-coords.t for additional examples
#
# Arguments:
#  $abs_nt_coords:  nucleotide coordinates in full sequence [1..seqlen]
#  $rel_pt_coords:  relative protein coordinates 
#  $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:   Coords string corresponding to $rel_pt_coords in nucleotide 
#            coordinates relative to $abs_nt_coords.
#
# Dies: if $rel_pt_coords coordinates imply positions outside 
#       $abs_nt_coords nucleotide coordinates (protein is too long)
#
#################################################################
sub vdr_CoordsProteinRelativeToAbsolute { 
  my $sub_name = "vdr_CoordsProteinRelativeToAbsolute";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($abs_nt_coords, $rel_pt_coords, $FH_HR) = @_;

  my $rel_nt_coords = vdr_CoordsProteinToNucleotide($rel_pt_coords, $FH_HR);
  return vdr_CoordsRelativeToAbsolute($abs_nt_coords, $rel_nt_coords, $FH_HR);
}

#################################################################
# Subroutine: vdr_CoordsProteinToNucleotide()
#
# Incept:     EPN, Fri Mar 20 07:40:19 2020
#
# Synopsis: Return nucleotide coordinates that correspond to the
#           protein coordinates in <$pt_coords>. 
#
#           Restricted to only work for plus strand segments, because
#           it doesn't make sense to have negative strand proteins but
#           code should be easily adaptable to work for minus strand
#           if desired.
#
#           A protein start position corresponds to the first
#           nucleotide position of the codon that encodes the first
#           amino acid of the protein. A protein stop position
#           corresponds to the third nucleotide position of the codon
#           that encodes the final amino acid of the protein.
#  
#           pt_coords            returns 
#           "1..10:+"            "1..30:+"
#           "1..10:+,15..30:+"   "1..30:+,43..90");
#
# See t/01-coords.t for additional examples
#
# Arguments:
#  $pt_coords:  protein coordinates 
#  $FH_HR:      REF to hash of file handles, including "log" and "cmd"
#
# Returns:   Coords string corresponding to $rel_pt_coords in nucleotide 
#            coordinates relative to $abs_nt_coords.
#
# Dies: if a segment in $pt_coords has a strand that is not "+"
#
#################################################################
sub vdr_CoordsProteinToNucleotide { 
  my $sub_name = "vdr_CoordsProteinToNucleotide";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($pt_coords, $FH_HR) = @_;

  my @start_A  = ();
  my @stop_A   = ();
  my @strand_A = ();
  vdr_FeatureStartStopStrandArrays($pt_coords, \@start_A, \@stop_A, \@strand_A, $FH_HR);
  my $nt_coords = "";
  my $start_coord = undef;
  my $stop_coord  = undef;
  my $nseg = scalar(@start_A);
  for(my $s = 0; $s < $nseg; $s++) { 
    if($nt_coords ne "") { $nt_coords .= ","; }
    if($strand_A[$s] ne "+") { 
      ofile_FAIL("ERROR in $sub_name, protein has segment that is not on + strand in coords string: $pt_coords", 1, $FH_HR);
    }
    my $start_coord = ($start_A[$s] * 3) - 2;
    my $stop_coord  = ($stop_A[$s] * 3);
    $nt_coords .= $start_coord . ".." . $stop_coord . ":+";
  }

  return $nt_coords;
}

#################################################################
# Subroutine: vdr_CoordsMergeAllAdjacentSegments()
#
# Incept:     EPN, Fri Mar 20 09:17:06 2020
#
# Synopsis: Merge any 'adjacent' segments in <$coords_str> and
#           return the possibly modified resulting coords string.
#           Actual work is done by vdr_CoordsMergeTwoSegmentsIfAdjacent()
#
# Arguments:
#  $coords: coords string to merge adjacent segments in
#  $FH_HR:  REF to hash of file handles, including "log" and "cmd"
#
# Returns:   Coords string with all adjacent segments merged together
#
# Dies: never
#
#################################################################
sub vdr_CoordsMergeAllAdjacentSegments { 
  my $sub_name = "vdr_CoordsMergeAllAdjacentSegments";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $FH_HR) = @_;

  my @sgm_A = split(",", $coords);
  my $nsgm = scalar(@sgm_A);

  # only 1 segment, can't merge anything 
  if($nsgm == 1) { return $coords; }

  my $ret_coords = "";
  my $cur_sgm = $sgm_A[0];
  my $merged_sgm = undef;
  for(my $i = 1; $i < $nsgm; $i++) { 
    $merged_sgm = vdr_CoordsMergeTwoSegmentsIfAdjacent($cur_sgm, $sgm_A[$i], $FH_HR);
    if($merged_sgm ne "") { 
      # we merged two segments, don't append yet
      $cur_sgm = $merged_sgm;
      # printf("in $sub_name, cur_sgm is merged_sum: $cur_sgm\n");
    }
    else { # did not merge, append $cur_sgm, then update $cur_sgm
      $ret_coords = vdr_CoordsAppendSegment($ret_coords, $cur_sgm);
      $cur_sgm = $sgm_A[$i];
    }
  }
  # add the final segment
  $ret_coords = vdr_CoordsAppendSegment($ret_coords, $cur_sgm);

  # printf("in $sub_name, input coords: $coords returning $ret_coords\n");
  return $ret_coords;
}

#################################################################
# Subroutine: vdr_CoordsMergeTwoSegmentsIfAdjacent()
#
# Incept:     EPN, Fri Mar 20 09:39:35 2020
#
# Synopsis: Merge two segments into one if they are 'adjacent'.
#           Two segments 1 and 2 are adjacent if they are either:
#           a) both are same strand and have fwd direction and
#              $coords_str and start_2 == (stop_1+1)
#           OR
#           b) both are same strand and have bck direction and
#              $coords_str and start_2 == (stop_1-1)
#
#           Where direction is defined as:
#           if strand is +
#             "fwd" if start_i <= stop_i
#             "bck" if start_i >  stop_i
#           if strand is -
#             "fwd" if start_i >= stop_i
#             "bck" if start_i <  stop_i
#       
#
# Arguments:
#  $sgm1:   segment 1
#  $sgm2:   segment 2
#  $FH_HR:  REF to hash of file handles, including "log" and "cmd"
#
# Returns:  A merged segment of segment 1 and 2 if they are adjacent
#           else "" if they are not adjacent
#
# Dies: If either $sgm1 or $sgm2 is not parseable
#
#################################################################
sub vdr_CoordsMergeTwoSegmentsIfAdjacent { 
  my $sub_name = "vdr_CoordsMergeTwoSegmentsIfAdjacent";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($sgm1, $sgm2, $FH_HR) = @_;

  my ($start1, $stop1, $strand1) = vdr_CoordsSegmentParse($sgm1, $FH_HR);
  my ($start2, $stop2, $strand2) = vdr_CoordsSegmentParse($sgm2, $FH_HR);
  my $dir1 = undef;
  my $dir2 = undef;
  if($strand1 eq "+") { 
    $dir1 = ($start1 <= $stop1) ? "fwd" : "bck";
  }
  else { 
    $dir1 = ($start1 >= $stop1) ? "bck" : "fwd";
  }
  if($strand2 eq "+") { 
    $dir2 = ($start2 <= $stop2) ? "fwd" : "bck";
  }
  else { 
    $dir2 = ($start2 >= $stop2) ? "bck" : "fwd";
  }

  # printf("in $sub_name, sgm1: $sgm1 sgm2: $sgm2 dir1: $dir1 dir2: $dir2\n");

  # return quick if we can
  if(($strand1 ne $strand2) || 
     ($dir1    ne $dir2)) { 
    return "";
  }

  # if we get here, $strand1 eq $strand2 and $dir1 eq $dir2
  if((($dir1 eq "fwd") && ($start2 == ($stop1+1))) ||
     (($dir1 eq "bck") && ($start2 == ($stop1-1)))) { 
    return vdr_CoordsSegmentCreate($start1, $stop2, $strand1, $FH_HR);
  }

  return "";  
}

#################################################################
# Subroutine: vdr_CoordsMaxLengthSegment()
#
# Incept:     EPN, Mon Mar 30 17:32:32 2020
#
# Synopsis: Return the maximum length segment and its length
#           from a coords string of one or more segments.
#           If multiple coords segments are tied for the maximum
#           length, return the first one. 
#
# Arguments:
#  $coords: segment 1
#  $FH_HR:  REF to hash of file handles, including "log" and "cmd"
#
# Returns:  Two values:
#           1. $argmax_sgm: he coords segment of maximum length 
#           2. $max_sgm_len: length of $argmax_coords_sgm
#
# Dies: If unable to parse $coords
#
#################################################################
sub vdr_CoordsMaxLengthSegment { 
  my $sub_name = "vdr_CoordsMaxLengthSegment";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $FH_HR) = @_;

  my @coords_tok_A  = split(",", $coords);
  my $nsgm = scalar(@coords_tok_A);
  my $argmax_sgm = $coords_tok_A[0];
  my $max_sgm_len = vdr_CoordsLength($coords_tok_A[0], $FH_HR);
  for(my $i = 1; $i < $nsgm; $i++) {
    my $sgm_len = vdr_CoordsLength($coords_tok_A[$i], $FH_HR);
    if($sgm_len > $max_sgm_len) {
      $max_sgm_len = $sgm_len;
      $argmax_sgm = $coords_tok_A[$i];
    }
  }
  return ($argmax_sgm, $max_sgm_len);
}

#################################################################
# Subroutine: vdr_EutilsFetchToFile()
# Incept:     EPN, Tue Mar 12 12:18:37 2019
#
# Synopsis: Fetch information for an accession using edirect.
#
# Arguments:
#  $out_file:  output file to create
#  $accn:      accession to fetch
#  $db:        database to fetch from (e.g. "nuccore", "protein")
#  $format:    format to fetch (e.g. "gpc", "ft", "fasta")
#  $nattempts: number of times to retry 
#  $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's a problem fetching the data
#################################################################
sub vdr_EutilsFetchToFile { 
  my $sub_name = "vdr_EutilsFetchToFile";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_file, $accn, $db, $format, $nattempts, $FH_HR) = @_;
  if((! defined $nattempts) || ($nattempts < 1)) { $nattempts = 1; }

  my $url = vdr_EutilsFetchUrl($accn, $db, $format);

  my $n = 0;
  my $fetched_str = undef;

  while(($n < $nattempts) && (! defined $fetched_str)) { 
    $fetched_str = get($url);
    $n++;
    sleep(3);
  }
  if(! defined $fetched_str) { 
    ofile_FAIL("ERROR in $sub_name, problem fetching $accn (undefined)", 1, $FH_HR); 
  }

  open(OUT, ">", $out_file) || ofile_FileOpenFailure($out_file, $sub_name, $!, "writing", $FH_HR);
  print OUT $fetched_str;
  close(OUT);

  return;
}

#################################################################
# Subroutine: vdr_EutilsFetchUrl()
# Incept:     EPN, Tue Mar 12 12:18:37 2019
#
# Synopsis: Return a url for an efetch command
#
# Arguments:
#  $accn:      accession to fetch
#  $db:        database to fetch from (e.g. "nuccore", "protein")
#  $format:    format to fetch (e.g. "gpc", "ft", "fasta")
#
# Returns:    void
#
# Dies:       never
#################################################################
sub vdr_EutilsFetchUrl { 
  my $sub_name = "vdr_EutilsFetchUrl";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($accn, $db, $format) = @_;

  return sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=%s&id=%s&rettype=%s&retmode=text", $db, $accn, $format);
}


#################################################################
# Subroutine: vdr_ModelInfoFileWrite()
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
sub vdr_ModelInfoFileWrite { 
  my $sub_name = "vdr_ModelInfoFileWrite";
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
    vdr_FeatureInfoValidateCoords($ftr_info_HAHR->{$mdl_name}, $mdl_info_AHR->[$mdl_idx]{"length"}, $FH_HR); 
  }

  # output 
  open(OUT, ">", $out_file) || ofile_FileOpenFailure($out_file, $sub_name, $!, "writing", $FH_HR);
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    print OUT ("MODEL $mdl_name");
    foreach $key (sort keys (%{$mdl_info_AHR->[$mdl_idx]})) { 
      $value = $mdl_info_AHR->[$mdl_idx]{$key};
      if($key =~ m/\:/) { 
        ofile_FAIL("ERROR in $sub_name, problem writing $out_file, illegal ':' character in model key $key for model $mdl_name", 1, $FH_HR);
      }
      if($value =~ m/\"/) { 
        ofile_FAIL("ERROR in $sub_name, problem writing $out_file, illegal '\"' character in model value $value for key $key for model $mdl_name", 1, $FH_HR);
      }
      if($key ne "name") { 
        print OUT (" $key:\"$value\"");
      }
    }
    print OUT ("\n");

    # define feature keys to ignore
    my %ftr_key_ignore_H = ();
    $ftr_key_ignore_H{"type"}           = 1; # this automatically gets added to @key_order_A, just so it goes first
    $ftr_key_ignore_H{"coords"}         = 1; # this automatically gets added to @key_order_A, just so it goes second
    $ftr_key_ignore_H{"parent_idx_str"} = 1; # this automatically gets added to @key_order_A, just so it goes third
    $ftr_key_ignore_H{"length"}         = 1; # will be inferred from coords
    $ftr_key_ignore_H{"3pa_ftr_idx"}    = 1; # will be inferred from coords and type
    $ftr_key_ignore_H{"outname"}        = 1; # will be inferred from product and gene (or lack of)
    $ftr_key_ignore_H{"5p_sgm_idx"}     = 1; # will be inferred from coords, when sgm_info_HA is created
    $ftr_key_ignore_H{"3p_sgm_idx"}     = 1; # will be inferred from coords, when sgm_info_HA is created
    $ftr_key_ignore_H{"location"}       = 1; # *could* (but won't be) inferred from coords

    $nftr = scalar(@{$ftr_info_HAHR->{$mdl_name}});
    # determine order of keys for this feature
    my @ftr_key_order_A  = ("type", "coords", "parent_idx_str");
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
            ofile_FAIL("ERROR in $sub_name, problem writing $out_file, illegal ':' character in feature key $key for model $mdl_name", 1, $FH_HR);
          }
          if($value =~ m/\"/) { 
            ofile_FAIL("ERROR in $sub_name, problem writing $out_file, illegal '\"' character in feature value $value for key $key for model $mdl_name", 1, $FH_HR);
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
# Subroutine: vdr_ModelInfoFileParse()
# Incept:     EPN, Fri Mar 15 05:15:23 2019
#
# Synopsis: Parse a model info file for >= 1 models and collect 
#           feature information for each model $model in 
#           @{$ftr_info_HAHR->{$model}}.
#
#           This subroutine validates that keys in @{$reqd_mdl_keys_AR}
#           are read and stored in $mdl_info_AHR, and that keys in 
#           @{$reqd_ftr_keys_AR} are read and stored in $ftr_info_HAHR.
# 
# Arguments:
#  $in_file:          input .minfo file to parse
#  $reqd_mdl_keys_AR: REF to array of required model   keys, e.g. ("name", "length")
#  $reqd_ftr_keys_AR: REF to array of required feature keys, e.g. ("type", "coords")
#  $mdl_info_AHR:     REF to array of hashes of model information, filled here
#  $ftr_info_HAHR:    REF to hash of array of hashes with information 
#                     on the features per model, filled here
#  $FH_HR:            REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if unable to parse $in_file
#             if a required mdl or ftr key does not exist
#################################################################
sub vdr_ModelInfoFileParse {
  my $sub_name = "vdr_ModelInfoFileParse";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($in_file, $reqd_mdl_keys_AR, $reqd_ftr_keys_AR, $mdl_info_AHR, $ftr_info_HAHR, $FH_HR) = @_;
  
  my $format_str = "# VADR model info (.minfo) format specifications:\n";
  $format_str   .= "# Lines prefixed with '#' are ignored.\n";
  $format_str   .= "# All other lines must begin with either: 'MODEL' or 'FEATURE'\n";
  $format_str   .= "# followed by one or more whitespace characters and then the model\n";
  $format_str   .= "# name <modelname> which cannot include whitespace.\n";
  $format_str   .= "# On each line after <modelname>, both MODEL and FEATURE lines must\n";
  $format_str   .= "# contain 0 or more <key>:<value> pairs meeting the following criteria.\n";
  $format_str   .= "# <key> must not include any whitespace or ':' characters\n";
  $format_str   .= "# <value> must start *and* end with '\"' but include no other '\"'\n";
  $format_str   .= "# characters (but <value> may include whitespace characters).\n";
  $format_str   .= "# To include multiple values for the same <key>, e.g. multiple 'note' qualifiers\n";
  $format_str   .= "# separate each qualifier value by the string ':GBSEP:' in the <value> field.\n";
  $format_str   .= "# <key>:<value> pairs must be separated by one or more whitespace characters.\n";
  $format_str   .= "# <modelname> and the first <key>:<value> pair must be separated by one or\n";
  $format_str   .= "# more whitespace characters.\n";

  # example lines:
  #MODEL NC_039477 cmfile:"test/test.vadr.cm"
  #FEATURE NC_039477 type:"gene" coords:"5..5104:+" gene:"ORF1"
  #FEATURE NC_039477 type:"CDS" coords:"5..5104:+" gene:"ORF1" product:"nonstructural polyprotein"

  my $mdl_name   = undef; # name of current model
  my $ftr_idx    = undef; # index of current feature
  my $mdl_idx    = -1;    # index of current model
  my %mdl_read_H = ();    # keeps track of which model names we've seen MODEL lines for, to avoid duplicates
  open(IN, $in_file) || ofile_FileOpenFailure($in_file, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    if($line !~ /^#/) { 
      # not a comment line
      my $orig_line = $line;
      chomp $line;
      my $is_model_line = 0; # set to 1 if line we are parsing is a MODEL line, else it's a FEATURE line
      if($line =~ /^MODEL\s+(\S+)\s*/) { 
        $mdl_name = $1;
        if(exists $mdl_read_H{$mdl_name}) { 
          ofile_FAIL("ERROR in $sub_name, problem parsing $in_file: read multiple MODEL lines for $mdl_name, should only be 1; line:\n$orig_line\n", 1, $FH_HR);
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
          ofile_FAIL("ERROR in $sub_name, problem parsing $in_file: read FEATURE line for model $mdl_name before a MODEL line for $mdl_name; line:\n$orig_line\n", 1, $FH_HR);
        }
        $ftr_idx = scalar(@{$ftr_info_HAHR->{$mdl_name}});
        # initialize ftr_info for this model/feature pair
        %{$ftr_info_HAHR->{$mdl_name}[$ftr_idx]} = (); 
        $line =~ s/^FEATURE\s+\S+\s*//; # remove FEATURE and model value
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $in_file, non-comment line does not start with 'MODEL <modelname>' or 'FEATURE <featurename>', line:\n$orig_line\n", 1, $FH_HR);
      }
      # if we get here we have either a MODEL or FEATURE line, parse the rest of it
      while($line ne "") { 
        if($line =~ /^([^\:\s]+)\:\"([^\"]+)\"\s*/) { 
          # key   must not include ':' or whitespace
          # value must begin and end with '"' but otherwise include no '"' characters
          my ($key, $value) = ($1, $2);
          if($is_model_line) { 
            if(exists $mdl_info_AHR->[$mdl_idx]{$key}) {
              ofile_FAIL("ERROR in $sub_name, problem parsing $in_file, read multiple values for key $key on MODEL line; line:\n$orig_line\n", 1, $FH_HR);
            }
            $mdl_info_AHR->[$mdl_idx]{$key} = $value;
          }
          else { # feature line
            if(exists $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key}) {
              ofile_FAIL("ERROR in $sub_name, problem parsing $in_file, read multiple values for key $key on FEATURE line; line:\n$orig_line\n", 1, $FH_HR);
            }
            $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key} = $value;
            # printf("\tadded ftr_info_HAHR->{$mdl_name}[$ftr_idx]{$key} as $value\n");
          }
          $line =~ s/^[^\:\s]+\:\"[^\"]+\"\s*//; # remove this key/value pair
        }
        else { 
          ofile_FAIL("ERROR in $sub_name, unable to parse $in_file, failed to parse key:value pairs in line:\n$orig_line\n$format_str\n", 1, $FH_HR);
        }
      } 
    }
  }
  close(IN);

  # verify we read what we need
  utl_AHValidate($mdl_info_AHR, $reqd_mdl_keys_AR, "ERROR in $sub_name, problem parsing $in_file, required MODEL key missing", $FH_HR);
  my $nmdl = scalar(@{$mdl_info_AHR});
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    utl_AHValidate($ftr_info_HAHR->{$mdl_name}, $reqd_ftr_keys_AR, "ERROR in $sub_name, problem parsing $in_file, required MODEL key missing for model " . $mdl_info_AHR->[$mdl_idx]{"name"}, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoMerge()
# Incept:     EPN, Sun May  5 10:32:40 2019
#
# Synopsis: Add info in one feature information arrays to another
#           after validating that they can be merged.
#
#           Features to be merged are identified as "consistent"
#           features between src_ftr_info_AHR and dst_ftr_info_AHR,
#           defined as those for which there is a subset of >=1
#           features that are in common (identical qualifier name and
#           value) between the two.
#
#           If any "inconsistent" features are identified between
#           src_ftr_info_AHR and dst_ftr_info_AHR the subroutine will
#           die. "Inconsistent" features are those for which a subset
#           of >=1 qualifiers have identical values but another subset
#           of >=1 qualifiers have different values.
# 
# Arguments:
#  $src_ftr_info_AHR:   REF to source feature info array of hashes to 
#                       add to $
#  $dst_ftr_info2_AHR:   REF to hash of array of hashes with information 
#                    on the features to add to  $ftr_info1_HAHR
#  $FH_HR:           REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if something is inconsistent between the two mdl_info_AHRs or
#             ftr_info_HAHRs that prevent merging
#################################################################
sub vdr_FeatureInfoMerge { 
  my $sub_name = "vdr_FeatureInfoMerge";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($src_ftr_info_AHR, $dst_ftr_info_AHR, $FH_HR) = @_;
  
  # for each feature in src_ftr_info_AHR, find the single consistent feature in dst_ftr_info_AHR
  my $src_nftr = scalar(@{$src_ftr_info_AHR});
  my $dst_nftr = scalar(@{$dst_ftr_info_AHR});
  my $src_ftr_key;
  for(my $src_ftr_idx = 0; $src_ftr_idx < $src_nftr; $src_ftr_idx++) { 
    my $found_consistent = 0;
    for(my $dst_ftr_idx = 0; $dst_ftr_idx < $dst_nftr; $dst_ftr_idx++) { 
      my $nconsistent   = 0;
      my $ninconsistent = 0;
      foreach $src_ftr_key (sort keys (%{$src_ftr_info_AHR->[$src_ftr_idx]})) { 
        if(defined $dst_ftr_info_AHR->[$dst_ftr_idx]{$src_ftr_key}) { 
          if($src_ftr_info_AHR->[$src_ftr_idx]{$src_ftr_key} eq
             $dst_ftr_info_AHR->[$dst_ftr_idx]{$src_ftr_key}) { 
            $nconsistent++;
          }
          else { 
            $ninconsistent++;
          }
        }
      }
      if(($nconsistent > 0) && ($ninconsistent == 0)) { 
        # we found a match, merge them
        if($found_consistent) { 
          ofile_FAIL("ERROR in $sub_name, trying to merge feature number " . ($src_ftr_idx + 1) . " but found more than one feature consistent with it", 1, $FH_HR);
        }
        foreach $src_ftr_key (sort keys (%{$src_ftr_info_AHR->[$src_ftr_idx]})) { 
          if(! defined $dst_ftr_info_AHR->[$dst_ftr_idx]{$src_ftr_key}) { 
            $dst_ftr_info_AHR->[$dst_ftr_idx]{$src_ftr_key} = 
                $src_ftr_info_AHR->[$src_ftr_idx]{$src_ftr_key};
          }
        }
        $found_consistent = 1;
      }
    }
    if(! $found_consistent) { 
      ofile_FAIL("ERROR in $sub_name, trying to merge feature number " . ($src_ftr_idx + 1) . " but did not find any features consistent with it (we can't add new features like this)", 1, $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine:  vdr_SplitFastaFile()
# Incept:      EPN, Tue Mar  1 09:30:10 2016
#
# Purpose: Split up a fasta file into <n> smaller files by calling
#          the esl-ssplit perl script.
#
# Arguments: 
#  $esl_ssplit:      path to the esl-ssplit.pl script to use
#  $fasta_file:      fasta file to split up
#  $nfiles:          desired number of files to split $fasta_file into, -1 for one file for each sequence
#  $opt_HHR:         REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:    Number of files actually created (can differ from requested
#             amount (which is $nfiles)).
#
# Dies:       if esl-ssplit command fails
#
################################################################# 
sub vdr_SplitFastaFile { 
  my $sub_name = "vdr_SplitFastaFile()";
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
# Subroutine: vdr_SplitNumSeqFiles()
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
sub vdr_SplitNumSeqFiles { 
  my $sub_name = "vdr_SplitNumSeqFiles";
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
# Subroutine: vdr_CdsFetchStockholmToFasta()
# Incept:     EPN, Thu Mar 14 12:30:33 2019
# 
# Purpose:    Given coordinates of all CDS features in %{$ftr_info_AHR}
#             fetch all the CDS for all sequences in the Stockholm alignment
#             and create a new output fasta file with just the CDS features.
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
sub vdr_CdsFetchStockholmToFasta { 
  my $sub_name = "vdr_CdsFetchStockholmToFasta";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $stk_file, $ftr_info_AHR, $FH_HR) = @_;

  my $msa = Bio::Easel::MSA->new({ fileLocation => $stk_file, isDna => 1});
  my $msa_has_rf = $msa->has_rf;

  # precompute start, stop, strand, for all features, so we don't have to redo this for each seq
  my @sgm_start_AA  = ();
  my @sgm_stop_AA   = ();
  my @sgm_strand_AA = ();
  vdr_FeatureInfoStartStopStrandArrays($ftr_info_AHR, \@sgm_start_AA, \@sgm_stop_AA, \@sgm_strand_AA, $FH_HR);

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
# Subroutine:  vdr_CmalignCheckStdOutput()
# Incept:      EPN, Wed Feb  6 14:18:59 2019
#
# Purpose:     Check cmalign output to see if it indicates that 
#              a cmalign run finished successfully, in error, or 
#              has not yet finished.
#              
# Arguments: 
#  $stdout_file:      path to the stdout file we will check
#  $ret_mxsize_R:     REF to required matrix size, only filled meaningfully if return value is '0'
#  $FH_HR:            REF to hash of file handles
# 
# Returns:     '1' if $stdout_file indicates cmalign job finished successfully
#              '0' if $stdout_file indicates cmalign job finished in error but in
#                  a way that is allowed, fills $$ret_mxsize_R
#             '-1' if $stdout_file indicates cmalign job is not yet finished
#                  or failed in some way we aren't looking for
#
# Dies: If $stdout_file does not exist or is empty
# 
################################################################# 
sub vdr_CmalignCheckStdOutput { 
  my $sub_name = "vdr_CmalignCheckStdOutput";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($stdout_file, $ret_mxsize_R, $FH_HR) = @_;
  if(defined $ret_mxsize_R) { 
    $$ret_mxsize_R = 0; # overwritten below if nec
  }

  if(! -e $stdout_file) { 
    ofile_FAIL("ERROR in $sub_name, cmalign stdout file $stdout_file does not exist", 1, $FH_HR);
  }
  if(! -s $stdout_file) { 
    ofile_FAIL("ERROR in $sub_name, cmalign stdout file $stdout_file exists but is empty", 1, $FH_HR);
  }

  # if we get here, the file exists and is non-empty
  my $final_line = `tail -n 1 $stdout_file`;
  chomp $final_line;
  if($final_line =~ m/\r$/) { chop $final_line; } # remove ^M if it exists
  if($final_line =~ m/\Q# CPU time/) { 
    return 1; 
  }
  else { 
    # job did NOT finish successfully, check for mx overflow error
    my $error_line = `grep ^Error $stdout_file | tail -n 1`;
    if($error_line =~ m/\r$/) { chop $error_line; } # remove ^M if it exists
    if($error_line =~ /Error: .+ alignment mxes need (\d+\.\d+)/) { 
      if(defined $ret_mxsize_R) { 
        $$ret_mxsize_R = $1;
      }
      return 0;
    }
    else { 
      return -1;
    }
  }
    
  return -1; # NEVER REACHED
}

#################################################################
# Subroutine: vdr_ParseSeqFileToSeqHash()
# Incept:     EPN, Mon May 20 12:19:29 2019
# 
# Purpose:    Parse an input sequence file using Bio::Easel::SqFile
#             and fill %{$seq_HR}. 
#
# Arguments:
#   $infile:   input file
#   $seq_HR:   sequence hash, key is sequence name, value is string of sequence
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#                    
# Returns: void
#
# Dies:    if we have trouble parsing the file
#
#################################################################
sub vdr_ParseSeqFileToSeqHash { 
  my $sub_name = "vdr_ParseSeqFileToSeqHash";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $seq_HR, $FH_HR) = @_;

  my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $infile }); # the sequence file object
  my $nseq = $sqfile->nseq_ssi;

  for(my $sidx = 0; $sidx < $nseq; $sidx++) {
    my ($sqname, undef) = $sqfile->fetch_seq_name_and_length_given_ssi_number($sidx);
    $seq_HR->{$sqname} = $sqfile->fetch_next_seq_to_sqstring();
  }
  
  return;
}

#################################################################
# Subroutine : vdr_CmalignParseIfile()
# Incept:      EPN, Thu Jan 31 13:06:54 2019
#
# Purpose:    Parse Infernal 1.1 cmalign --ifile output and store
#             results in %{$seq_inserts_HHR}.
#
#             %{$seq_inserts_HHR} is a 2D hash, with
#             key 1: sequence name
#             key 2: one of 'spos', 'epos', 'ins'
#             $seq_inserts_HHR->{}{"spos"} is starting model position of alignment
#             $seq_inserts_HHR->{}{"epos"} is ending model position of alignment
#             $seq_inserts_HHR->{}{"ins"} is the insert string in the format:
#             <mdlpos_1>:<uapos_1>:<inslen_1>;...<mdlpos_n>:<uapos_n>:<inslen_n>;
#             for n inserts, where insert x is defined by:
#             <mdlpos_x> is model position after which insert occurs 0..mdl_len (0=before first pos)
#             <uapos_x> is unaligned sequence position of the first aligned nt
#             <inslen_x> is length of the insert
#             
# Arguments: 
#  $ifile_file:       ifile file to parse
#  $mdl_name_R:       REF to model name read from ifile, defined here
#  $mdl_len_R:        REF to model length read from ifile, defined here
#  $seq_inserts_HHR:  REF to hash of hashes with insert information, added to here
#  $FH_HR:            REF to hash of file handles
#
# Returns:    void
#
# Dies:       if unable to parse the ifile
#
################################################################# 
sub vdr_CmalignParseIfile { 
  my $sub_name = "vdr_CmalignParseIfile()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ifile_file, $mdl_name_R, $mdl_len_R, $seq_inserts_HHR, $FH_HR) = @_;
  
  open(IN, $ifile_file) || ofile_FileOpenFailure($ifile_file, $sub_name, $!, "reading", $FH_HR);

  my $line_ctr = 0;  # counts lines in ifile_file
  while(my $line = <IN>) { 
    $line_ctr++;
    if($line !~ m/^\#/) { 
      chomp $line;
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
      # 2 types of lines, those with 2 tokens, and those with 4 or more tokens
      #norovirus.NC_039477 7567
      #gi|669176088|gb|KM198574.1| 7431 17 7447  2560 2539 3  2583 2565 3
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) == 2) {
        ($mdl_name, $mdl_len) = (@_);
        
      } 
      elsif(scalar(@el_A) >= 4) { 
        my $nel = scalar(@el_A); 
        if((($nel - 4) % 3) != 0) { # check number of elements makes sense
          ofile_FAIL("ERROR in $sub_name, unexpected number of elements ($nel) in ifile line in $ifile_file on line $line_ctr:\n$line\n", 1, $FH_HR);
        }          
        my ($seqname, $seqlen, $spos, $epos) = ($el_A[0], $el_A[1], $el_A[2], $el_A[3]);
        if(! defined $seq_inserts_HHR->{$seqname}) { 
          # initialize
          %{$seq_inserts_HHR->{$seqname}} = ();
        }
        # create the insert string
        my $insert_str = "";
        for(my $el_idx = 4; $el_idx < scalar(@el_A); $el_idx += 3) { 
          $insert_str .= $el_A[$el_idx] . ":" . $el_A[$el_idx+1] . ":" . $el_A[$el_idx+2] . ";"; 
        }
        $seq_inserts_HHR->{$seqname}{"spos"} = $spos;
        $seq_inserts_HHR->{$seqname}{"epos"} = $epos;
        $seq_inserts_HHR->{$seqname}{"ins"}  = $insert_str;
      }
    }
  }
  close(IN);
  
  return;
}

#################################################################
# Subroutine : vdr_CmalignWriteIfile()
# Incept:      EPN, Fri Apr  3 11:17:49 2020
#
# Purpose:    Write an Infernal 1.1 cmalign --ifile given
#             insert information in %{$seq_inserts_HHR}.
#
#             See vdr_CmalignParseIfile()'s 'Purpose' for
#             details on format of the cmalign ifile.
#
# Arguments: 
#  $ifile_file:       ifile file to parse
#  $mdl_name:  
#  $seq_inserts_HHR:  REF to hash of hashes with insert information, added to here
#  $FH_HR:            REF to hash of file handles
#
# Returns:    void
#
# Dies:       if unable to parse the ifile
#
################################################################# 
sub vdr_CmalignParseIfile { 
  my $sub_name = "vdr_CmalignParseIfile()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ifile_file, $seq_inserts_HHR, $FH_HR) = @_;
  
  open(OUT, ">", $ifile_file) || ofile_FileOpenFailure($ifile_file, $sub_name, $!, "writing", $FH_HR);

  my $line_ctr = 0;  # counts lines in ifile_file
  
  while(my $line = <IN>) { 
    $line_ctr++;
    if($line !~ m/^\#/) { 
      chomp $line;
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
      # 2 types of lines, those with 2 tokens, and those with 4 or more tokens
      #norovirus.NC_039477 7567
      #gi|669176088|gb|KM198574.1| 7431 17 7447  2560 2539 3  2583 2565 3
      my @el_A = split(/\s+/, $line);
      if   (scalar(@el_A) == 2) { ; } # ignore these lines
      elsif(scalar(@el_A) >= 4) { 
        my $nel = scalar(@el_A); 
        if((($nel - 4) % 3) != 0) { # check number of elements makes sense
          ofile_FAIL("ERROR in $sub_name, unexpected number of elements ($nel) in ifile line in $ifile_file on line $line_ctr:\n$line\n", 1, $FH_HR);
        }          
        my ($seqname, $seqlen, $spos, $epos) = ($el_A[0], $el_A[1], $el_A[2], $el_A[3]);
        if(! defined $seq_inserts_HHR->{$seqname}) { 
          # initialize
          %{$seq_inserts_HHR->{$seqname}} = ();
        }
        # create the insert string
        my $insert_str = "";
        for(my $el_idx = 4; $el_idx < scalar(@el_A); $el_idx += 3) { 
          $insert_str .= $el_A[$el_idx] . ":" . $el_A[$el_idx+1] . ":" . $el_A[$el_idx+2] . ";"; 
        }
        $seq_inserts_HHR->{$seqname}{"spos"} = $spos;
        $seq_inserts_HHR->{$seqname}{"epos"} = $epos;
        $seq_inserts_HHR->{$seqname}{"ins"}  = $insert_str;
      }
    }
  }
  close(IN);
  
  return;
}

###########################################################################
# the next line is critical, a perl module must return a true value
return 1;
###########################################################################
