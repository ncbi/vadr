#!/usr/bin/env perl
# 
# version: 1.6 [Nov 2023]
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
use strict;
use warnings;
use Cwd;
use LWP::Simple; 
# use LWP::Protocol::https; 
# above line is purposefully commented out, b/c LWP::Protocol::https is only
# needed for v-build.pl, and many users only use v-annotate.pl only, which can run without it
# when a better solution for checking for modules during installation is found, uncomment it

require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seq.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

#########################################################################################
#
# List of subroutines in this file, divided into categories. 
#
# Subroutines related to features or segments:
# vdr_FeatureInfoImputeCoords()
# vdr_FeatureInfoImputeLength()
# vdr_FeatureInfoInitializeParentIndexStrings()
# vdr_FeatureInfoImputeOutname()
# vdr_FeatureInfoImpute3paFtrIdx()
# vdr_FeatureInfoImputeByOverlap()
# vdr_FeatureInfoInitializeMiscNotFailure()
# vdr_FeatureInfoInitializeIsDeletable()
# vdr_FeatureInfoInitializeAlternativeFeatureSet()
# vdr_FeatureInfoInitializeAlternativeFeatureSetSubstitution()
# vdr_FeatureInfoValidateMiscNotFailure()
# vdr_FeatureInfoValidateIsDeletable()
# vdr_FeatureInfoValidateAlternativeFeatureSet()
# vdr_FeatureInfoValidateAndConvertAlternativeFeatureSetSubstitution()
# vdr_FeatureInfoValidateCanonSpliceSites()
# vdr_FeatureInfoValidateExceptionKeys()
# vdr_FeatureInfoStartStopStrandArrays()
# vdr_FeatureInfoCountType()
# vdr_FeatureInfoValidateCoords()
# vdr_FeatureInfoValidateParentIndexStrings()
# vdr_FeatureInfoChildrenArrayOfArrays()
# vdr_FeatureInfoMapFtrTypeIndicesToFtrIndices()
# vdr_FeatureInfoMerge()
# vdr_FeatureInfoCdsStartStopCodonCoords()
# vdr_FeatureInfoMaxNumCdsSegments()
#
# vdr_SegmentInfoPopulate()
# 
# vdr_FeatureTypeAndTypeIndexString()
# vdr_FeatureTypeIndex()
# vdr_FeatureTypeIsCds()
# vdr_FeatureTypeIsCdsOrIdStartAndStop()
# vdr_FeatureTypeIsMatPeptide()
# vdr_FeatureTypeIsGene()
# vdr_FeatureTypeIsCdsOrGene()
# vdr_FeatureTypeIsCdsOrMatPeptide()
# vdr_FeatureTypeIsCdsOrMatPeptideOrIdCoords()
# vdr_FeatureTypeIsCdsOrMatPeptideOrIdStartAndStop()
# vdr_FeatureTypeIsCdsOrMatPeptideOrGene()
# vdr_FeatureTypeCanBecomeMiscFeature()
# vdr_FeatureOutType()
# vdr_FeatureNumSegments()
# vdr_FeatureRelativeSegmentIndex()
# vdr_Feature5pMostPosition()
# vdr_Feature3pMostPosition()
# vdr_FeatureParentIndex()
# vdr_FeatureSummarizeSegment()
# vdr_FeatureStartStopStrandArrays()
# vdr_FeatureSummaryStrand()
# vdr_FeaturePositionSpecificValueBreakdown()
# vdr_FeatureCoordsListBreakdown()
#
# vdr_SegmentStartIdenticalToCds()
# vdr_SegmentStopIdenticalToCds()
#
# Subroutines related to alerts:
# vdr_AlertInfoInitialize()
# vdr_AlertInfoAdd()
# vdr_AlertInfoSetFTableInvalidatedBy()
# vdr_AlertInfoSetCausesFailure()
# vdr_AlertInfoSetMiscNotFailure()
# vdr_AlertInfoDump()
# 
# Subroutines related to feature-specific alerts (misc_not_failure support):
# vdr_FeatureAlertCausesFailure()
# vdr_FeatureAlertIsMiscNotFailure()
# 
# Subroutines related to alert exceptions:
# vdr_ExceptionCoordsAndValuesToSegmentsAndValues()
# vdr_ExceptionSegmentAndValueToPositionsAndValues()
# vdr_ExceptionSegmentAndValueParse()
# vdr_ExceptionCoordsAndValuesValidate()
# 
# Subroutines related to parallelization on the compute farm:
# vdr_ParseQsubFile()
# vdr_SubmitJob()
# vdr_WaitForFarmJobsToFinish()
#
# Subroutines related to sequence and model coordinates: 
# vdr_CoordsToSegments()
# vdr_CoordsSegmentParse()
# vdr_CoordsSegmentCreate()
# vdr_CoordsValidate()
# vdr_CoordsSegmentValidate()
# vdr_CoordsSinglePositionSegmentCreate()
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
# vdr_CoordsRelativeSingleCoordToAbsolute()
# vdr_CoordsProteinRelativeToAbsolute()
# vdr_CoordsProteinToNucleotide()
# vdr_CoordsMergeAllAdjacentSegments()
# vdr_CoordsMergeTwoSegmentsIfAdjacent()
# vdr_CoordsMaxLengthSegment()
# vdr_CoordsFromStartStopStrandArrays()
# vdr_CoordsSegmentActualToFractional()
# vdr_CoordsSegmentFractionalToActual()
#
# Subroutines related to eutils:
# vdr_EutilsFetchToFile()
# vdr_EutilsFetchUrl()
# 
# Subroutines related to model info files:
# vdr_ModelInfoFileWrite()
# vdr_ModelInfoFileParse()
# vdr_ModelInfoCoordListValueBreakdown()
# vdr_ModelInfoValidateExceptionKeys()
#
# Subroutines related to cmalign output:
# vdr_CmalignCheckStdOutput()
# vdr_CmalignParseInsertFile()
# vdr_CmalignWriteInsertFile()
#
# Subroutines related to the --split option:
# vdr_CmalignCheckStdOutput()
# vdr_CmalignParseInsertFile()
# vdr_CmalignWriteInsertFile()
# 
# Other subroutines related to running infernal programs
# vdr_CmemitConsensus()
# 
# Subroutines related to merging output (--split):
# vdr_MergeOutputConcatenateOnly()
# vdr_MergeOutputConcatenatePreserveSpacing()
# vdr_MergeOutputGetFileList()
# vdr_MergeOutputMdlTabularFile()
# vdr_MergeOutputAlcTabularFile()
# vdr_MergeFrameshiftStockholmFiles()
# vdr_MergeAlignments()
# 
# Subroutines supporting N-replacement (-r):
# vdr_ReplacePseudoCoordsStringCreate()
# vdr_ReplacePseudoCoordsStringParse()
#
# Subroutines related to backwards compatibility:
# vdr_BackwardsCompatibilityExceptions()
#
# Miscellaneous subroutines:
# vdr_SplitFastaFile()
# vdr_SplitNumSeqFiles()
# vdr_CdsFetchStockholmToFasta()
# vdr_ParseSeqFileToSeqHash()
# vdr_FrameAdjust()
# vdr_WriteCommandScript()
# vdr_GlsearchFormat3And9CToStockholmAndInsertFile()
# vdr_CigarToInsertsHash()
# vdr_CigarToPositionMap()
# vdr_UpdateInsertTokenInInsertString()
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
# Subroutine: vdr_FeatureInfoInitializeMiscNotFailure
# Incept:     EPN, Fri Feb  5 11:44:11 2021
# 
# Purpose:    Set "misc_not_failure" value to 0 for any feature 
#             in which it is not already defined in @{$ftr_info_AHR}.
#             If $force_zero, set all values to 0 even if they are
#             already defined.
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $force_zero:    '1' to set values to '0' for all features, even if already defined
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub vdr_FeatureInfoInitializeMiscNotFailure {
  my $sub_name = "vdr_FeatureInfoInitializeMiscNotFailure";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $force_zero, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(($force_zero) || (! defined $ftr_info_AHR->[$ftr_idx]{"misc_not_failure"})) { 
      $ftr_info_AHR->[$ftr_idx]{"misc_not_failure"} = 0;
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoInitializeIsDeletable
# Incept:     EPN, Tue Sep 28 21:09:22 2021
# 
# Purpose:    Set "is_deletable" value to 0 for any feature 
#             in which it is not already defined in @{$ftr_info_AHR}.
#             If $force_zero, set all values to 0 even if they are
#             already defined.
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $force_zero:    '1' to set values to '0' for all features, even if already defined
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub vdr_FeatureInfoInitializeIsDeletable {
  my $sub_name = "vdr_FeatureInfoInitializeIsDeletable";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $force_zero, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(($force_zero) || (! defined $ftr_info_AHR->[$ftr_idx]{"is_deletable"})) { 
      $ftr_info_AHR->[$ftr_idx]{"is_deletable"} = 0;
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoInitializeAlternativeFeatureSet
# Incept:     EPN, Tue Oct 12 19:43:46 2021
# 
# Purpose:    Set "alternative_ftr_set" value to "" for any feature 
#             in which it is not already defined in @{$ftr_info_AHR}.
#             If $force_empty, set all values to "" even if they are
#             already defined.
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $force_empty:   '1' to set values to "" for all features, even if already defined
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub vdr_FeatureInfoInitializeAlternativeFeatureSet {
  my $sub_name = "vdr_FeatureInfoInitializeAlternativeFeatureSet";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $force_empty, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(($force_empty) || (! defined $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"})) { 
      $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"} = "";
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoInitializeAlternativeFeatureSetSubstitution
# Incept:     EPN, Thu Oct 14 21:22:20 2021
# 
# Purpose:    Set "alternative_ftr_set_subn" value to "" for any feature 
#             in which it is not already defined in @{$ftr_info_AHR}.
#             If $force_empty, set all values to "" even if they are
#             already defined.
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $force_empty:   '1' to set values to "" for all features, even if already defined
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub vdr_FeatureInfoInitializeAlternativeFeatureSetSubstitution {
  my $sub_name = "vdr_FeatureInfoInitializeAlternativeFeatureSetSubstitution";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $force_empty, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(($force_empty) || (! defined $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set_subn"})) { 
      $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set_subn"} = "";
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoInitializeCanonSpliceSites
# Incept:     EPN, Mon Aug  7 13:34:57 2023
# 
# Purpose:    Set "canon_splice_sites" value to 0 for any feature 
#             in which it is not already defined in @{$ftr_info_AHR}.
#             If $force_empty, set all values to 0 even if they are
#             already defined.
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $force_one:     '1' to set values to 1 for all features, even if already defined as 0
#   $force_zero:    '1' to set values to 0 for all features, even if already defined as 1
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub vdr_FeatureInfoInitializeCanonSpliceSites {
  my $sub_name = "vdr_FeatureInfoInitializeCanonSpliceSites";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $force_one, $force_zero, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($force_one) { 
      $ftr_info_AHR->[$ftr_idx]{"canon_splice_sites"} = 1;
    }
    if(($force_zero) || (! defined $ftr_info_AHR->[$ftr_idx]{"canon_splice_sites"})) { 
      $ftr_info_AHR->[$ftr_idx]{"canon_splice_sites"} = 0;
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoValidateMiscNotFailure
# Incept:     EPN, Fri Feb  5 11:45:29 2021
# 
# Purpose:    Validate "misc_not_failure" values are either 0 or 1.
#             Should probably be called after vdr_FeatureInfoInitializeMiscNotFailure()
# 
#             I considered enforcing that "misc_not_failure" could only exist
#             for feature types that can become misc_features but decided 
#             against it so that 'gene' features can have alerts that would
#             be fatal except that a misc_not_failure=1 value causes them not
#             to be fatal. These won't be turned into misc_features because
#             output_feature_table() won't let them be, but that's ok.
#
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if misc_not_failure is invalid for any feature
#
#################################################################
sub vdr_FeatureInfoValidateMiscNotFailure {
  my $sub_name = "vdr_FeatureInfoValidateMiscNotFailure";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
  
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  my $nftr = scalar(@{$ftr_info_AHR});
  my $fail_str = ""; # added to if any elements are out of range
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(! defined $ftr_info_AHR->[$ftr_idx]{"misc_not_failure"}) { 
      $fail_str .= "ftr_idx: $ftr_idx, undefined\n"; 
    }
    elsif(($ftr_info_AHR->[$ftr_idx]{"misc_not_failure"} != 1) && ($ftr_info_AHR->[$ftr_idx]{"misc_not_failure"} != 0)) { 
      $fail_str .= "ftr_idx: $ftr_idx, " . $ftr_info_AHR->[$ftr_idx]{"misc_not_failure"} . " != 0 and != 1\n"; 
    }
# decided against this enforcement (see Purpose):
#    elsif(($ftr_info_AHR->[$ftr_idx]{"misc_not_failure"} == 1) && (! vdr_FeatureTypeCanBecomeMiscFeature($ftr_info_AHR, $ftr_idx))) { 
#      $fail_str .= "ftr_idx: $ftr_idx, " . $ftr_info_AHR->[$ftr_idx]{"misc_not_failure"} . " is 1 but type is " . $ftr_info_AHR->[$ftr_idx]{"type"} . " and that type cannot become a misc_feature (hard-coded)\n";
#    }
  }
  
  if($fail_str ne "") { 
    ofile_FAIL("ERROR in $sub_name, some misc_not_failure values are invalid or don't make sense:\n$fail_str\n", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoValidateIsDeletable
# Incept:     EPN, Tue Sep 28 21:09:10 2021
# 
# Purpose:    Validate "is_deletable" values are either 0 or 1.
#             Should probably be called after vdr_FeatureInfoInitializeIsDeletable()
#
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if is_deletable is invalid for any feature
#
#################################################################
sub vdr_FeatureInfoValidateIsDeletable {
  my $sub_name = "vdr_FeatureInfoValidateIsDeletable";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
  
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  my $nftr = scalar(@{$ftr_info_AHR});
  my $fail_str = ""; # added to if any elements are out of range
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(! defined $ftr_info_AHR->[$ftr_idx]{"is_deletable"}) { 
      $fail_str .= "ftr_idx: $ftr_idx, undefined\n"; 
    }
    elsif(($ftr_info_AHR->[$ftr_idx]{"is_deletable"} != 1) && ($ftr_info_AHR->[$ftr_idx]{"is_deletable"} != 0)) { 
      $fail_str .= "ftr_idx: $ftr_idx, " . $ftr_info_AHR->[$ftr_idx]{"is_deletable"} . " != 0 and != 1\n"; 
    }
  }
  
  if($fail_str ne "") { 
    ofile_FAIL("ERROR in $sub_name, some is_deletable values are invalid or don't make sense:\n$fail_str\n", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoValidateAlternativeFeatureSet
# Incept:     EPN, Tue Sep 28 21:09:10 2021
# 

# Purpose:    Validate "alternative_ftr_set" values are either "" or another
#             string. If another string, each other string must be the
#             value for "alternative_ftr_set" in more than one
#             feature. Also ensure that for any sets that have >= 1
#             children, all the features in that set are all the
#             children of the same parent.
#           
#             Should probably be called after
#             vdr_FeatureInfoInitializeAlternativeFeatureSet() 
#             and 
#             vdr_FeatureInfoValidateParentIndexStrings()
#
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    '1' if there are any 'alternative_ftr_set' values ne ""
#             '0' if all 'alternative_ftr_set' values are ""
# 
# Dies:       if any alternative_ftr_set values are undefined
#             if any alternative_ftr_set values exist only once 
#
#################################################################
sub vdr_FeatureInfoValidateAlternativeFeatureSet {
  my $sub_name = "vdr_FeatureInfoValidateAlternativeFeatureSet";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
  
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  my $nftr     = scalar(@{$ftr_info_AHR});
  my $ret_val  = 0; # set to '1' if we see any values ne ""
  my $fail_str = ""; # added to if any elements are out of range
  my %set_HA = (); # key is set value, array is feature indices in that set
  my $ftr_idx = undef;

  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(! defined $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"}) {
      $fail_str .= "ftr_idx: $ftr_idx, undefined\n"; 
    }
    else { 
      my $value = $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"};
      if($value ne "") { 
        $ret_val = 1;
        if(! defined $set_HA{$value}) { 
          @{$set_HA{$value}} = ();
        }
        push(@{$set_HA{$value}}, $ftr_idx);
      }
    }    
  }

  # make sure that each alternative_ftr_set has >= 2 members
  # and that for any set that has >= 1 children, all members are children with the same parent
  foreach my $key (sort keys (%set_HA)) { 
    my $nset = scalar(@{$set_HA{$key}});
    if($nset == 1) { 
      $fail_str .= "alternative_ftr_set value: $key exists only once, each value must exist at least twice\n"; 
    }
    else { 
      my $nchildren = 0;
      my $common_parent_idx = undef;
      foreach $ftr_idx (@{$set_HA{$key}}) { 
        my $parent_idx = $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"};
        if((defined $parent_idx) && ($parent_idx ne "GBNULL")) { 
          $nchildren++;
          if(! defined $common_parent_idx) { 
            $common_parent_idx = $parent_idx;
          }
          elsif($parent_idx != $common_parent_idx) { 
            $fail_str .= "ftr_idx: $ftr_idx is child of parent $parent_idx but >= 1 other feature(s) in same set ($key) have a different parent ($common_parent_idx), this is not allowed\n";
          }
        }
      }
      # make sure if any members are children, then all members are children
      if(($nchildren != 0) && ($nchildren != $nset)) { 
        $fail_str .= "for alternative_ftr_set with key $key, some but not all members are children of $common_parent_idx\n";
      }        
    }
  }

  if($fail_str ne "") { 
    ofile_FAIL("ERROR in $sub_name, some alternative_ftr_set values are invalid or don't make sense:\n$fail_str\n", 1, $FH_HR);
  }

  return $ret_val;
}

#################################################################
# Subroutine: vdr_FeatureInfoValidateAndConvertAlternativeFeatureSetSubstitution
# Incept:     EPN, Fri Oct 15 10:07:54 2021
# 


# Purpose:    Validate "alternative_ftr_set_subn" values are either "",
#             "<s>.<d1>" or "<d2>" where <s> is a valid
#             "alternative_ftr_set" for a set other than the set of
#             the current feature, and "<d1>" is an index 1..n where n
#             is the size of the set "<s>".  Alternatively, for
#             backwards compatibility, "<d2>" can be a feature index
#             <0..nftr-1> as long as it isn't self idx.
#
#             Should probably be called after
#             vdr_FeatureInfoInitializeAlternativeFeatureSetSubstitution() 
#
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if any alternative_ftr_set_subn values are undefined
#             if any alternative_ftr_set_subn values are invalid
#
#################################################################
sub vdr_FeatureInfoValidateAndConvertAlternativeFeatureSetSubstitution {
  my $sub_name = "vdr_FeatureInfoValidateAndConvertAlternativeFeatureSetSubstitution";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
  
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  my $nftr     = scalar(@{$ftr_info_AHR});
  my $fail_str = ""; # added to if any elements are out of range

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $subn_val = $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set_subn"};
    if(! defined $subn_val) { 
      $fail_str .= "ftr_idx: $ftr_idx, undefined\n"; 
    }
    elsif($subn_val ne "") { 
      if($subn_val =~ /^(\S+)\.(\d+)$/) { # e.g. "attachment(cds).2"
        my ($set, $set_idx) = ($1, $2);
        # make sure $set is not the set that this ftr belongs to
        if((defined $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"}) && 
           ($ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"} eq $set)) { 
          $fail_str .= "ftr_idx: $ftr_idx, " . $subn_val . " implies the substitution set is the same as the ftr set ($set) for this feature.\n";
        }
        else { 
          # find the ftr $set.$set_idx corresponds to:
          my $cur_set_idx = 0;
          my $found_flag = 0; # set to 1 if we find it
          for(my $ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
            if((defined $ftr_info_AHR->[$ftr_idx2]{"alternative_ftr_set"}) && 
               ($ftr_info_AHR->[$ftr_idx2]{"alternative_ftr_set"} eq $set)) { 
              $cur_set_idx++;
              if($cur_set_idx == $set_idx) { 
                # rewrite value
                $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set_subn"} = $ftr_idx2;
                $found_flag = 1;
              }
            }
          }
          if(! $found_flag) { 
            $fail_str .= "ftr_idx: $ftr_idx, " . $subn_val . " implies substitution ftr is from alternative_ftr_set $set index $set_idx, but found $cur_set_idx features with that alternative_ftr_set value.\n";
          }
        }
      }
      elsif($subn_val =~ /^\d+$/) { 
        if($subn_val < 0) { 
          $fail_str .= "ftr_idx: $ftr_idx, " . $subn_val . " < 0\n"; 
        }
        elsif($subn_val >= $nftr) { 
          $fail_str .= "ftr_idx: $ftr_idx, " . $subn_val . " >= $nftr (num features, should be 0.." . ($nftr-1) . ")\n";
        }
        elsif($subn_val == $ftr_idx) { 
          $fail_str .= "ftr_idx: $ftr_idx, is its own substitute, this is not allowed\n";
        }
      }
      # else valid feature index
    }
  }

  if($fail_str ne "") { 
    ofile_FAIL("ERROR in $sub_name, some alternative_ftr_set_subn index strings are undefined or don't make sense:\n$fail_str\n", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoValidateCanonSpliceSites
# Incept:     EPN, Mon Aug  7 13:37:15 2023
# 
# Purpose:    Validate "canon_splice_sites" values are either 0 or 1.
#             Should probably be called after vdr_FeatureInfoInitializeCanonSpliceSites()
#
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if canon_splice_sites is invalid for any feature
#
#################################################################
sub vdr_FeatureInfoValidateCanonSpliceSites {
  my $sub_name = "vdr_FeatureInfoValidateCanonSpliceSites";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
  
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  my $nftr = scalar(@{$ftr_info_AHR});
  my $fail_str = ""; # added to if any elements are out of range
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(! defined $ftr_info_AHR->[$ftr_idx]{"canon_splice_sites"}) { 
      $fail_str .= "ftr_idx: $ftr_idx, undefined\n"; 
    }
    elsif(($ftr_info_AHR->[$ftr_idx]{"canon_splice_sites"} != 1) && ($ftr_info_AHR->[$ftr_idx]{"canon_splice_sites"} != 0)) { 
      $fail_str .= "ftr_idx: $ftr_idx, " . $ftr_info_AHR->[$ftr_idx]{"canon_splice_sites"} . " != 0 and != 1\n"; 
    }
  }
  
  if($fail_str ne "") { 
    ofile_FAIL("ERROR in $sub_name, some canon_splice_sites values are invalid or don't make sense:\n$fail_str\n", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoValidateExceptionsKeys()
# Incept:     EPN, Tue Sep 12 13:51:19 2023
#
# Purpose:    Validate any feature info keys that pertain to
#             alert exceptions.
# 
#             Uses 'exc_key' value in alt_info_HH
#             to validate all ftr_info_AH keys that end in 
#             '_exc', and dies if any such keys are not valid.
#             
# Arguments: 
#  $ftr_info_AHR: ref to the feature info array of hashes
#  $alt_info_HHR: ref to the alert info hash of hashes
#  $FH_HR:        ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $ftr_info_AHR->[]{"*_exc"} = $exc_str is not valid 
#             (if no alt_info_HH{<code>}{"exc_key"} eq $exc_str)
#
################################################################# 
sub vdr_FeatureInfoValidateExceptionKeys {
  my $sub_name = "vdr_FeatureInfoValidateExceptionKeys";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ftr_info_AHR, $alt_info_HHR, $FH_HR) = @_;
  
  # small optimization, make temporary hash of all $alt_code_HH{<code>}{"exc_key"}
  my %tmp_key_H = ();
  foreach my $code (sort keys %{$alt_info_HHR}) { 
    if(defined $alt_info_HHR->{$code}{"exc_key"}) { 
      $tmp_key_H{$alt_info_HHR->{$code}{"exc_key"}} = $alt_info_HHR->{$code}{"exc_type"};
    }
  }
  
  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    foreach my $exc_key (sort keys %{$ftr_info_AHR->[$ftr_idx]}) { 
      if($exc_key =~ /^.+\_exc$/) { 
        if(! defined $tmp_key_H{$exc_key}) { 
          ofile_FAIL("ERROR, in $sub_name, invalid exception key $exc_key read in model info file", 1, $FH_HR);
        }
        # now validate format of value
        if($tmp_key_H{$exc_key} eq "coords-only") { 
          if(! vdr_CoordsValidate($ftr_info_AHR->[$ftr_idx]{$exc_key}, $FH_HR)) { 
            ofile_FAIL(sprintf("ERROR, in $sub_name, invalid format for exception key: $exc_key, type 'coords-only', read %s.\nBut expected coords string, with spans separated by commas, e.g. \"11..15:+\" or \"11..15:+,18..18:+\".", $ftr_info_AHR->[$ftr_idx]{$exc_key}), 1, $FH_HR);
          }
        }
        else { # coords-value
          # validate it by parsing it
          my %sgm_value_H = (); # we won't use this but need to pass it to vdr_ExceptionCoordsAndValuesToSegmentsAndValues()
          my $errmsg = sprintf("ERROR, in $sub_name, invalid format for exception key: $exc_key, type 'coords-value', read %s.\nBut expected coords-value string, with span separated by commas, e.g. \"11..15:+:23\" or \"11..15:+:23,18..18:+:36\".", $ftr_info_AHR->[$ftr_idx]{$exc_key});
          vdr_ExceptionCoordsAndValuesToSegmentsAndValues($ftr_info_AHR->[$ftr_idx]{$exc_key}, $errmsg, \%sgm_value_H, $FH_HR);
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
    # this sub vdr_FeatureStartStopStrandArrays will die if $ftr_info_AHR->[$ftr_idx]{"coords"} is in incorrect format
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
#             or a valid feature index [0..$nftr-1]. Also make sure that
#             all parents are not themselves children of another feature.
#             This subroutine should probably be called after 
#             vdr_FeatureInfoInitializeParentIndexStrings().
#             
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#             if any feature that is a parent is also a child of a different feature
#             if any feature is a child of itself
#################################################################
sub vdr_FeatureInfoValidateParentIndexStrings {
  my $sub_name = "vdr_FeatureInfoValidateParentIndexStrings";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
  
  my ($ftr_info_AHR, $FH_HR) = @_;
  
  my $nftr = scalar(@{$ftr_info_AHR});
  my $fail_str = ""; # added to if any elements are out of range
  my @parent_idx_A = (); # array of all features that are parents >= 1 children
  my %parent_idx_H = (); # hash used to avoid duplicates in @parent_A
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $parent_idx_str = $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"};
    if(! defined $parent_idx_str) { 
      $fail_str .= "ftr_idx: $ftr_idx, undefined\n"; 
    }
    elsif($parent_idx_str ne "GBNULL") { 
      if($parent_idx_str !~ /^\d+$/) { 
        $fail_str .= "ftr_idx: $ftr_idx, " . $parent_idx_str . " not an integer\n"; 
      }
      elsif($parent_idx_str < 0) { 
        $fail_str .= "ftr_idx: $ftr_idx, " . $parent_idx_str . " < 0\n"; 
      }
      elsif($parent_idx_str >= $nftr) { 
        $fail_str .= "ftr_idx: $ftr_idx, " . $parent_idx_str . " >= $nftr (num features, should be 0.." . ($nftr-1) . ")\n";
      }
      elsif($parent_idx_str == $ftr_idx) { 
        $fail_str .= "ftr_idx: $ftr_idx, is its own parent, this is not allowed\n";
      }
      else { # valid feature index
        # printf("in $sub_name, adding $parent_idx_str to parent_idx_A\n");
        if(! defined $parent_idx_H{$parent_idx_str}) { 
          push(@parent_idx_A, $parent_idx_str);
          $parent_idx_H{$parent_idx_str} = 1;
        }
      }
    }
  }

  # make sure all parents are not children themselves
  foreach my $parent_idx (@parent_idx_A) { 
    my $parent_idx_str = $ftr_info_AHR->[$parent_idx]{"parent_idx_str"};
    if((defined $parent_idx_str) && ($parent_idx_str ne "GBNULL")) { 
      $fail_str .= "ftr_idx: $parent_idx, is a parent but is also a child of parent " . $parent_idx_str . ", this is not allowed\n";
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
#   $i_am_child_AR:  REF to array of 1/0 values for each feature, '1' if child of type $type_or_undef,
#                    '0' if not, FILLED HERE, can be undef
#   $children_AAR:   REF to array of arrays of children feature indices, FILLED HERE, can be undef
#   $FH_HR:          REF to hash of file handles
# 
# Returns:     Number of features that are children of type $type_or_undef (or any type if $type_or_undef is not defined)
# 
################################################################# 
sub vdr_FeatureInfoChildrenArrayOfArrays { 
  my $nargs_expected = 5;
  my $sub_name = "vdr_FeatureInfoChildrenArrayOfArrays";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $type_or_undef, $i_am_child_AR, $children_AAR, $FH_HR) = @_;

  my @children_AA  = ();
  my @i_am_child_A = ();
  my $nftr = scalar(@{$ftr_info_AHR});

  my ($parent_ftr_idx, $child_ftr_idx);

  # initialize
  for($parent_ftr_idx = 0; $parent_ftr_idx < $nftr; $parent_ftr_idx++) { 
    @{$children_AA[$parent_ftr_idx]} = ();
  }

  # fill
  my $ret_val = 0;
  for($child_ftr_idx = 0; $child_ftr_idx < $nftr; $child_ftr_idx++) { 
    $i_am_child_A[$child_ftr_idx] = 0;
    if((! defined $type_or_undef) || ($ftr_info_AHR->[$child_ftr_idx]{"type"} eq $type_or_undef)) { 
      if($ftr_info_AHR->[$child_ftr_idx]{"parent_idx_str"} ne "GBNULL") { 
        $parent_ftr_idx = $ftr_info_AHR->[$child_ftr_idx]{"parent_idx_str"};
        push(@{$children_AA[$parent_ftr_idx]}, $child_ftr_idx);
        $i_am_child_A[$child_ftr_idx] = 1;
        $ret_val++;
      }
    }
  }
  
  if(defined $i_am_child_AR) { 
    @{$i_am_child_AR} = (@i_am_child_A);
  }
  if(defined $children_AAR) { 
    @{$children_AAR} = (@children_AA);
  }

  return $ret_val;
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
#                       add to $dst_ftr_info_AHR
#  $dst_ftr_info_AHR:   REF to hash of array of hashes we are adding to
#  $FH_HR:              REF to hash of file handles, including "log" and "cmd"
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
          ofile_FAIL("ERROR, in $sub_name, trying to merge feature number " . ($src_ftr_idx + 1) . " but found more than one feature consistent with it", 1, $FH_HR);
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
      ofile_FAIL("ERROR, in $sub_name, trying to merge feature number " . ($src_ftr_idx + 1) . " but did not find any features consistent with it (we can't add new features like this)", 1, $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureInfoCdsStartStopCodonCoords()
# Incept:     EPN, Mon Sep 13 17:55:41 2021
#
# Synopsis: Return coords strings of for start and stop codon
#           coordinates for all CDS. Start/stop codon coords will be
#           concatenated into long coords strings with multiple
#           segments. If a single start or stop codon is actually
#           N>1 segments, it will contribute N segments to the 
#           returned string.
# 
# Arguments:
#  $ftr_info_AHR: REF to array of hashes with information on the features, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    Two values:
#             $start_codon_coords: >= 1 segment coords string with coordinates of all start codons,
#                                  "" if no CDS exist in @{$ftr_info_AHR}
#             $stop_codon_coords:  >= 1 segment coords string with coordinates of all stop codons
#                                  "" if no CDS exist in @{$ftr_info_AHR}
#
# Dies:       If a CDS has length < 6
#################################################################
sub vdr_FeatureInfoCdsStartStopCodonCoords {
  my $sub_name = "vdr_FeatureInfoCdsStartStopCodonCoords";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_AHR, $FH_HR) = @_;

  my $ret_start_codon_coords = "";
  my $ret_stop_codon_coords = "";
  # for each feature in src_ftr_info_AHR, find the single consistent feature in dst_ftr_info_AHR
  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") { 
      my $coords = $ftr_info_AHR->[$ftr_idx]{"coords"};
      my $length = vdr_CoordsLength($coords, $FH_HR);
      if($length < 6) { 
        ofile_FAIL("ERROR, in $sub_name, CDS (ftr_idx: $ftr_idx) coords length is below 6", 1, $FH_HR); 
      }
      if($ret_start_codon_coords ne "") { 
        $ret_start_codon_coords .= ",";
        $ret_stop_codon_coords  .= ",";
      }
      $ret_start_codon_coords .= vdr_CoordsRelativeToAbsolute($coords, "1..3:+", $FH_HR);
      $ret_stop_codon_coords  .= vdr_CoordsRelativeToAbsolute($coords, ($length-2) . ".." . $length . ":+", $FH_HR);
    }
  }

  return($ret_start_codon_coords, $ret_stop_codon_coords);
}

#################################################################
# Subroutine: vdr_FeatureInfoMaxNumCdsSegments()
# Incept:     EPN, Wed Jun  7 12:18:03 2023
#
# Synopsis: Return the maximum number of segments in any
#           CDS feature. 
# 
# Arguments:
#  $ftr_info_AHR: REF to array of hashes with information on the features, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    Maximum number of segments in any CDS.
#################################################################
sub vdr_FeatureInfoMaxNumCdsSegments { 
  my $sub_name = "vdr_FeatureInfoMaxNumCdsSegments";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_AHR, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  my $ret_val = 0;
  my $nsgm = 0;
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") { 
      $nsgm = vdr_FeatureNumSegments($ftr_info_AHR, $ftr_idx);
      if($nsgm > $ret_val) { 
        $ret_val = $nsgm;
      }
    }
  }

  return $ret_val;
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
  my $nftr = utl_AHValidate($ftr_info_AHR, \@keys_A, "ERROR, in $sub_name", $FH_HR);

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
# Subroutine: vdr_FeatureTypeIsCdsOrIdStartAndStop()
# Incept:     EPN, Thu Dec  9 10:22:21 2021
#
# Purpose:    Is feature $ftr_idx either a CDS or does it have 
#             same start/stop coords as another feature
#             that is a CDS?
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
sub vdr_FeatureTypeIsCdsOrIdStartAndStop { 
  my $sub_name = "vdr_FeatureTypeIsCdsOrIdStartAndStop";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) { 
    return 1;
  }
  else { 
    my $nftr = scalar(@{$ftr_info_AHR});
    my $start1 = vdr_Feature5pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, undef);
    my $stop1  = vdr_Feature3pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, undef);
    for(my $ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
      if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx2)) { 
        my $start2 = vdr_Feature5pMostPosition($ftr_info_AHR->[$ftr_idx2]{"coords"}, undef);
        my $stop2  = vdr_Feature3pMostPosition($ftr_info_AHR->[$ftr_idx2]{"coords"}, undef);
        if(($start1 == $start2) && ($stop1 == $stop2)) { 
          return 1; 
        }
      }
    }
  }
  return 0;
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
# Subroutine: vdr_FeatureTypeIsCdsOrGene()
# Incept:     EPN, Tue Feb  9 15:16:21 2021
#
# Purpose:    Is feature $ftr_idx a CDS or gene?
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
sub vdr_FeatureTypeIsCdsOrGene { 
  my $sub_name = "vdr_FeatureTypeIsCdsOrGene";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return (($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") || 
          ($ftr_info_AHR->[$ftr_idx]{"type"} eq "gene")) ? 1 : 0;
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
# Subroutine: vdr_FeatureTypeIsCdsOrMatPeptideOrIdCoords()
# Incept:     EPN, Fri Feb  5 14:23:53 2021
#
# Purpose:    Is feature $ftr_idx either a CDS or mature peptide
#             or does it have equivalent coords to another feature
#             that is a CDS or mature peptide?
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
sub vdr_FeatureTypeIsCdsOrMatPeptideOrIdCoords { 
  my $sub_name = "vdr_FeatureTypeIsCdsOrMatPeptideOrIdCoords";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  if(vdr_FeatureTypeIsCdsOrMatPeptide($ftr_info_AHR, $ftr_idx)) { 
    return 1;
  }
  else { 
    my $nftr = scalar(@{$ftr_info_AHR});
    for(my $ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
      if((vdr_FeatureTypeIsCdsOrMatPeptide($ftr_info_AHR, $ftr_idx2)) && 
         ($ftr_info_AHR->[$ftr_idx]{"coords"} eq $ftr_info_AHR->[$ftr_idx2]{"coords"})) { 
        return 1; 
      }
    }
  }
  return 0;
}

#################################################################
# Subroutine: vdr_FeatureTypeIsCdsOrMatPeptideOrIdStartAndStop()
# Incept:     EPN, Wed Sep  1 15:42:19 2021
#
# Purpose:    Is feature $ftr_idx either a CDS or mature peptide
#             or does it have same start/stop coords as another feature
#             that is a CDS or mature peptide?
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
sub vdr_FeatureTypeIsCdsOrMatPeptideOrIdStartAndStop { 
  my $sub_name = "vdr_FeatureTypeIsCdsOrMatPeptideOrIdStartAndStop";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  if(vdr_FeatureTypeIsCdsOrMatPeptide($ftr_info_AHR, $ftr_idx)) { 
    return 1;
  }
  else { 
    my $nftr = scalar(@{$ftr_info_AHR});
    my $start1 = vdr_Feature5pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, undef);
    my $stop1  = vdr_Feature3pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, undef);
    for(my $ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
      if(vdr_FeatureTypeIsCdsOrMatPeptide($ftr_info_AHR, $ftr_idx2)) { 
        my $start2 = vdr_Feature5pMostPosition($ftr_info_AHR->[$ftr_idx2]{"coords"}, undef);
        my $stop2  = vdr_Feature3pMostPosition($ftr_info_AHR->[$ftr_idx2]{"coords"}, undef);
        if(($start1 == $start2) && ($stop1 == $stop2)) { 
          return 1; 
        }
      }
    }
  }
  return 0;
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
# Subroutine: vdr_FeatureTypeCanBecomeMiscFeature()
# Incept:     EPN, Tue Feb  9 13:15:31 2021
#
# Purpose:    Can feature $ftr_idx become a misc_feature?
#             Currently the definition of which feature types
#             cannot become misc_features is hard-coded in this
#             subroutine.
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
sub vdr_FeatureTypeCanBecomeMiscFeature { 
  my $sub_name = "vdr_FeatureTypeCanBecomeMiscFeature";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return (($ftr_info_AHR->[$ftr_idx]{"type"} ne "gene") && 
          ($ftr_info_AHR->[$ftr_idx]{"type"} ne "5'UTR") && 
          ($ftr_info_AHR->[$ftr_idx]{"type"} ne "3'UTR") && 
          ($ftr_info_AHR->[$ftr_idx]{"type"} ne "operon")) ? 1 : 0;
}

#################################################################
# Subroutine: vdr_FeatureOutType()
# Incept:     EPN, Thu Jun 29 09:38:05 2023
#
# Purpose:    Return the 'output type' for a feature. This is
#             the value for the 'out_type' key if it is defined
#             in %ftr_info_AHR, else it is the value for the 
#             'type' key.
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $ftr_idx:        feature index
#
# Returns:    $ftr_info_AHR->[$ftr_idx]{"out_type"} if defined
#             else $ftr_info_AHR->[$ftr_idx]{"type"}
#
# Dies:       never; does not validate anything.
#
################################################################# 
sub vdr_FeatureOutType { 
  my $sub_name = "vdr_FeatureOutType";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return (defined $ftr_info_AHR->[$ftr_idx]{"out_type"}) ? 
      $ftr_info_AHR->[$ftr_idx]{"out_type"} : $ftr_info_AHR->[$ftr_idx]{"type"};
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
    ofile_FAIL("ERROR, in $sub_name, unable to parse ftr_info_HA coords string " . $coords, 1, $FH_HR); 
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
    ofile_FAIL("ERROR, in $sub_name, unable to parse ftr_info_HA coords string " . $coords, 1, $FH_HR); 
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
# Purpose:    Return a string indicating what segment this is
#             for features with multiple segments and return 
#             the empty string for features with 1 segment.
#
# Arguments: 
#  $ftr_info_AHR: ref to feature info array of hashes, PRE-FILLED
#  $sgm_info_AHR: ref to segment info array of hashes, PRE-FILLED
#  $sgm_idx:      model index
#
# Returns:    "" if this is the only segment for this feature
#             string like ", segment 1 of 2", if not
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
  my $nsgm = ($ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"} - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}) + 1;
  if($nsgm > 1) { 
    return sprintf(", segment %d of %d", ($sgm_idx - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}) + 1, $nsgm);
  }

  return ""; # return "" if $nsgm == 1;
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
    ofile_FAIL("ERROR, in $sub_name, coords is undefined", 1, $FH_HR); 
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
    else { ofile_FAIL("ERROR, in $sub_name, unable to determine strands in coords $coords", 1, $FH_HR); }
  }

  if(($npos >  0) && ($nneg >  0)) { return "!"; }
  if(($npos >  0) && ($nneg == 0)) { return "+"; }
  if(($npos == 0) && ($nneg >  0)) { return "-"; }
  if(($npos == 0) && ($nneg == 0)) { 
    ofile_FAIL("ERROR, in $sub_name, unable to determine strands in coords $coords", 1, $FH_HR); 
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
#             of: "<d>:<s>" or "<d>..<d>:<s>" (for a range)
#             separated by ";" if more than one.
#
#             If $ftr_info_AHR->[$ftr_idx]{$key} does not exist just
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
      elsif($tok =~ /^(\d+)\.\.(\d+)\:(\S+)$/) { 
        my ($start, $stop, $value) = ($1, $2, $3);
        if($start > $stop) { 
          ofile_FAIL("ERROR, in $sub_name, unable to parse token $tok (stop > start) parsed out of " . $ftr_info_AHR->[$ftr_idx]{$key}, 1, $FH_HR);
        }
        for(my $i = $start; $i <= $stop; $i++) { 
          $HR->{$i} = $value;
        }
      }
      else { 
        ofile_FAIL("ERROR, in $sub_name, unable to parse token $tok parsed out of " . $ftr_info_AHR->[$ftr_idx]{$key}, 1, $FH_HR);
      }
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureCoordsListValueBreakdown()
# Incept:     EPN, Mon Oct 18 14:26:17 2021
#
# Purpose:    Breakdown a list of coords values 
#             from a string in %{$ftr_info_AHR->[$ftr_idx]}
#             and fill @{$AR} with key/value pairs.
# 
#             String must be in format of one or more single
#             segment coords tokens, "<d>..<d>:[+-]" separated 
#             by ";" if more than one.
#
#             If $ftr_info_AHR->[$ftr_idx]{$key} does not exist just
#             return.
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $ftr_idx:        feature index
#  $key:            key in $ftr_info_AHR->[$ftr_idx]
#  $AR:             ref to array to fill
#  $FH_HR:          ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $ftr_info_AHR->[$ftr_idx] exists but cannot
#             be parsed.
#
################################################################# 
sub vdr_FeatureCoordsListValueBreakdown { 
  my $sub_name = "vdr_FeatureCoordsListValueBreakdown";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx, $key, $AR, $FH_HR) = @_;
 
  if(defined $ftr_info_AHR->[$ftr_idx]{$key}) { 
    my @tok_A = split(";", $ftr_info_AHR->[$ftr_idx]{$key});
    foreach my $tok (@tok_A) { 
      if($tok =~ /^\d+\.\.\d+\:[\+\-]$/) { 
        push(@{$AR}, $tok);
      }
      else { 
        ofile_FAIL("ERROR, in $sub_name, unable to parse coords token $tok parsed out of " . $ftr_info_AHR->[$ftr_idx]{$key}, 1, $FH_HR);
      }
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureLengthBetweenAdjacentSegments()
# Incept:     EPN, Fri Jul 28 14:08:20 2023
#
# Synopsis: Return number of positions in between two adjacent 
#           segments for a feature.
# 
# Arguments:
#  $ftr_info_AHR: ref to feature info array of hashes, PRE-FILLED
#  $sgm_info_AHR: ref to segment info array of hashes, PRE-FILLED
#  $ftr_idx:      feature index we are interested in
#  $sgm_idx_5p:   segment index *within feature* (0 for feature's first segment)
#                 we will return distance between this segment and next segment
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:  number of positions between $sgm_idx_5p and ($sgm_idx_5p+1)
#
# Dies: if $sgm_idx_5p is the final segment for the feature
#       if $sgm_idx_5p is higher than number of segments for this feature
#       if the segment and the next segment are on different strands
#################################################################
sub vdr_FeatureLengthBetweenAdjacentSegments { 
  my $sub_name = "vdr_FeatureLengthBetweenAdjacentSegments";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($ftr_info_AHR, $sgm_info_AHR, $ftr_idx, $sgm_idx_5p, $FH_HR) = @_;
  
  my $nsgm = $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"} - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} + 1;
  if($sgm_idx_5p >= $nsgm) { 
    ofile_FAIL(sprintf("ERROR, in $sub_name, ftr_idx $ftr_idx has $nsgm segments, and requesting distance between segments %d and %d", $sgm_idx_5p, ($sgm_idx_5p+1)), 1, $FH_HR); 
  }

  my $sgm_idx     = $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} + $sgm_idx_5p;
  my $nxt_sgm_idx = $sgm_idx + 1;
  my $strand = $sgm_info_AHR->[$sgm_idx]{"strand"};
  if($sgm_info_AHR->[$nxt_sgm_idx]{"strand"} ne $strand) { 
    ofile_FAIL("ERROR, in $sub_name, segment and next segment have different strands", 1, $FH_HR); 
  }
    
  # determine start/stop and length of region between the two segments
  # we can't use strand agnostic 
  my $region_start  = ($strand eq "+") ? $sgm_info_AHR->[$sgm_idx]{"stop"}+1      : $sgm_info_AHR->[$sgm_idx]{"stop"}-1;
  my $region_stop   = ($strand eq "+") ? $sgm_info_AHR->[$nxt_sgm_idx]{"start"}-1 : $sgm_info_AHR->[$nxt_sgm_idx]{"start"}+1;
  my $region_length = ($strand eq "+") ? ($region_stop - $region_start) + 1       : ($region_start - $region_stop) + 1;

  # printf("in $sub_name, returning $region_length\n");

  return $region_length;
}

#################################################################
# Subroutine: vdr_SegmentStartIdenticalToCds()
# Incept:     EPN, Mon Feb  1 15:58:25 2021
#
# Purpose:    Determine if the start of a segment is identical
#             to the start position of any CDS feature. Return '1'
#             if yes, '0' if no.
# 
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $sgm_info_AHR:   ref to the segment info array of hashes 
#  $sgm_idx:        segment index
#  $FH_HR:          ref to hash of file handles, including "log" and "cmd"
#
# Returns:    1 or 0 (see 'Purpose')
#
# Dies:       If $sgm_info_AHR->[$sgm_idx] does not exist
#
################################################################# 
sub vdr_SegmentStartIdenticalToCds { 
  my $sub_name = "vdr_SegmentStartIdenticalToCds";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ftr_info_AHR, $sgm_info_AHR, $sgm_idx, $FH_HR) = @_;
  
  my $nftr = scalar(@{$ftr_info_AHR});
  my $nsgm = scalar(@{$sgm_info_AHR});

  if(($sgm_idx < 0) || ($sgm_idx >= $nsgm)) { 
    ofile_FAIL("ERROR, in $sub_name, invalid sgm idx: $sgm_idx", 1, $FH_HR);
  }
     
  my $sgm_start = $sgm_info_AHR->[$sgm_idx]{"start"};
  
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) {
    if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) { 
      my $sgm_5p_idx = $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"};
      my $ftr_start = $sgm_info_AHR->[$sgm_5p_idx]{"start"};
      if($ftr_start == $sgm_start) { return 1; }
    }
  }

  return 0;
}

#################################################################
# Subroutine: vdr_SegmentStopIdenticalToCds()
# Incept:     EPN, Mon Feb  1 15:58:25 2021
#
# Purpose:    Determine if the stop of a segment is identical
#             to the stop position of any CDS feature. Return '1'
#             if yes, '0' if no.
# 
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $sgm_info_AHR:   ref to the segment info array of hashes 
#  $sgm_idx:        segment index
#  $FH_HR:          ref to hash of file handles, including "log" and "cmd"
#
# Returns:    1 or 0 (see 'Purpose')
#
# Dies:       If $sgm_info_AHR->[$sgm_idx] does not exist
#
################################################################# 
sub vdr_SegmentStopIdenticalToCds { 
  my $sub_name = "vdr_SegmentStopIdenticalToCds";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ftr_info_AHR, $sgm_info_AHR, $sgm_idx, $FH_HR) = @_;
  
  my $nftr = scalar(@{$ftr_info_AHR});
  my $nsgm = scalar(@{$sgm_info_AHR});

  if(($sgm_idx < 0) || ($sgm_idx >= $nsgm)) { 
    ofile_FAIL("ERROR, in $sub_name, invalid sgm idx: $sgm_idx", 1, $FH_HR);
  }
     
  my $sgm_stop = $sgm_info_AHR->[$sgm_idx]{"stop"};
  
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) {
    if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) { 
      my $sgm_3p_idx = $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"};
      my $ftr_stop = $sgm_info_AHR->[$sgm_3p_idx]{"stop"};
      if($ftr_stop == $sgm_stop) { return 1; }
    }
  }

  return 0;
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
    ofile_FAIL("ERROR, in $sub_name, alert info hash of arrays already has at least one key", 1, $FH_HR);
  }

  # add each alert code, this function will die if we try to add the same code twice, or if something is wrong 
  # with how we try to add it (args to vdr_AlertInfoAdd don't pass the contract check)

  vdr_AlertInfoAdd($alt_info_HHR, "noannotn", "sequence",
                   "NO_ANNOTATION", # short description
                   "no significant similarity detected", # long  description
                   1, 1, 1, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "revcompl", "sequence",
                   "REVCOMPLEM", # short description
                   "sequence appears to be reverse complemented", # long description
                   1, 1, 1, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "qstsbgrp", "sequence",
                   "QUESTIONABLE_SPECIFIED_SUBGROUP", # short description
                   "best overall model is not from specified subgroup", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "qstgroup", "sequence",
                   "QUESTIONABLE_SPECIFIED_GROUP", # short description
                   "best overall model is not from specified group", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "incsbgrp", "sequence",
                   "INCORRECT_SPECIFIED_SUBGROUP", # short description
                   "score difference too large between best overall model and best specified subgroup model", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "incgroup", "sequence",
                   "INCORRECT_SPECIFIED_GROUP", # short description
                   "score difference too large between best overall model and best specified group model", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "ambgnt5s", "sequence",
                   "AMBIGUITY_AT_START", # short description
                   "first nucleotide of the sequence is an ambiguous nucleotide", # long  description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "ambgnt3s", "sequence",
                   "AMBIGUITY_AT_END", # short description
                   "final nucleotide of the sequence is an ambiguous nucleotide", # long  description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "lowcovrg", "sequence",
                   "LOW_COVERAGE", # short description, 
                   "low sequence fraction with significant similarity to homology model", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "indfclas", "sequence",
                   "INDEFINITE_CLASSIFICATION", # short description
                   "low score difference between best overall model and second best model (not in best model's subgroup)", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "lowscore", "sequence",
                   "LOW_SCORE", # short description
                   "score to homology model below low threshold", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "biasdseq", "sequence",
                   "BIASED_SEQUENCE", # short description
                   "high fraction of score attributed to biased sequence composition", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "dupregin", "sequence",
                   "DUPLICATE_REGIONS", # short description
                   "similarity to a model region occurs more than once", # long description
                   0, 1, 0, 0, "dupregin_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "discontn", "sequence",
                   "DISCONTINUOUS_SIMILARITY", # short description
                   "not all hits are in the same order in the sequence and the homology model", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indfstrn", "sequence",
                   "INDEFINITE_STRAND", # short description
                   "significant similarity detected on both strands", # long description
                   0, 1, 0, 0, "indfstr_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim5s", "sequence",
                   "LOW_SIMILARITY_START", # short description
                   "significant similarity not detected at 5' end of the sequence", # long description
                   0, 1, 0, 0, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim3s", "sequence",
                   "LOW_SIMILARITY_END", # short description
                   "significant similarity not detected at 3' end of the sequence", # long description
                   0, 1, 0, 0, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsimis", "sequence",
                   "LOW_SIMILARITY", # short description
                   "internal region without significant similarity", # long description
                   0, 1, 0, 0, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "extrant5", "sequence",
                   "EXTRA_SEQUENCE_START", # short description
                   "extra sequence detected 5' of expected sequence start", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "extrant3", "sequence",
                   "EXTRA_SEQUENCE_END", # short description
                   "extra sequence detected 3' of expected sequence end", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);
 
  vdr_AlertInfoAdd($alt_info_HHR, "unexdivg", "sequence",
                   "UNEXPECTED_DIVERGENCE", # short description
                   "sequence is too divergent to confidently assign nucleotide-based annotation", # long description
                   1, 1, 1, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "unjoinbl", "sequence",
                   "UNJOINABLE_SUBSEQ_ALIGNMENTS", # short description
                   "inconsistent alignment of overlapping region between seed and flanking region", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "noftrann", "sequence",
                   "NO_FEATURES_ANNOTATED", # short description
                   "sequence similarity to homology model does not overlap with any features", # long description
                   1, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "noftrant", "sequence",
                   "NO_FEATURES_ANNOTATED", # short description
                   "all annotated features are too short to output to feature table", # long description
                   1, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "nmiscftr", "sequence",
                   "TOO_MANY_MISC_FEATURES", # short description
                   "too many features are reported as misc_features", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "ftskipfl", "sequence",
                   "UNREPORTED_FEATURE_PROBLEM", # short description
                   "only fatal alerts are for feature(s) not output to feature table", # long description
                   1, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "deletins", "sequence",
                   "DELETION_OF_FEATURE", # short description
                   "internal deletion of a complete feature", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "deletina", "sequence",
                   "DELETION_OF_FEATURE", # short description
                   "allowed internal deletion of a complete feature", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "ambgntrp", "sequence",
                   "N_RICH_REGION_NOT_REPLACED", # short description
                   "N-rich region of unexpected length not replaced", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "mutstart", "feature",
                   "MUTATION_AT_START", # short description
                   "expected start codon could not be identified", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "mutendcd", "feature",
                   "MUTATION_AT_END", # short description
                   "expected stop codon could not be identified, predicted CDS stop by homology is invalid", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "mutendns", "feature",  
                   "MUTATION_AT_END", # short description
                   "expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted start codon", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "mutendex", "feature",
                   "MUTATION_AT_END", # short description
                   "expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "unexleng", "feature",  
                   "UNEXPECTED_LENGTH", # short description
                   "length of complete coding (CDS or mat_peptide) feature is not a multiple of 3", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "cdsstopn", "feature",
                   "CDS_HAS_STOP_CODON", # short description
                   "in-frame stop codon exists 5' of stop position predicted by homology to reference", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "cdsstopp", "feature",
                   "CDS_HAS_STOP_CODON", # short description
                   "stop codon in protein-based alignment", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "fsthicft", "feature",
                   "POSSIBLE_FRAMESHIFT_HIGH_CONF", # short description
                   "high confidence possible frameshift in CDS (frame not restored before end)", # long description
                   0, 1, 0, 1, "fst_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "fsthicfi", "feature",
                   "POSSIBLE_FRAMESHIFT_HIGH_CONF", # short description
                   "high confidence possible frameshift in CDS (frame restored before end)", # long description
                   0, 1, 0, 1, "fst_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "fstlocft", "feature",
                   "POSSIBLE_FRAMESHIFT_LOW_CONF", # short description
                   "low confidence possible frameshift in CDS (frame not restored before end)", # long description
                   0, 0, 0, 1, "fst_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "fstlocfi", "feature",
                   "POSSIBLE_FRAMESHIFT_LOW_CONF", # short description
                   "low confidence possible frameshift in CDS (frame restored before end)", # long description
                   0, 0, 0, 1, "fst_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "fstukcft", "feature",
                   "POSSIBLE_FRAMESHIFT", # short description
                   "possible frameshift in CDS (frame not restored before end)", # long description
                   0, 1, 0, 1, "fst_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "fstukcfi", "feature",
                   "POSSIBLE_FRAMESHIFT", # short description
                   "possible frameshift in CDS (frame restored before end)", # long description
                   0, 1, 0, 1, "fst_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "mutspst5", "feature",
                   "MUTATION_AT_SPLICE_SITE", # short description
                   "expected splice site at 5' end of intron (GT) could not be identified", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "mutspst3", "feature",
                   "MUTATION_AT_SPLICE_SITE", # short description
                   "expected splice site at 3' end of intron (AG) could not be identified", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "peptrans", "feature",
                   "PEPTIDE_TRANSLATION_PROBLEM", # short description
                   "mat_peptide may not be translated because its parent CDS has a problem", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "pepadjcy", "feature",
                   "PEPTIDE_ADJACENCY_PROBLEM", # short description
                   "predictions of two mat_peptides expected to be adjacent are not adjacent", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indfantp", "feature",
                   "INDEFINITE_ANNOTATION", # short description
                   "protein-based search identifies CDS not identified in nucleotide-based search", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indfantn", "feature",
                   "INDEFINITE_ANNOTATION", # short description
                   "nucleotide-based search identifies CDS not identified in protein-based search", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf5gap", "feature",
                   "INDEFINITE_ANNOTATION_START", # short description
                   "alignment to homology model is a gap at 5' boundary", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf5lcc", "feature",
                   "INDEFINITE_ANNOTATION_START", # short description
                   "alignment to homology model has low confidence at 5' boundary for feature that is or matches a CDS", # long description
                   0, 0, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf5lcn", "feature",
                   "INDEFINITE_ANNOTATION_START", # short description
                   "alignment to homology model has low confidence at 5' boundary for feature that does not match a CDS", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf5plg", "feature",
                   "INDEFINITE_ANNOTATION_START", # short description
                   "protein-based alignment extends past nucleotide-based alignment at 5' end", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf5pst", "feature",
                   "INDEFINITE_ANNOTATION_START", # short description
                   "protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf3gap", "feature",
                   "INDEFINITE_ANNOTATION_END", # short description
                   "alignment to homology model is a gap at 3' boundary", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf3lcc", "feature",
                   "INDEFINITE_ANNOTATION_END", # short description
                   "alignment to homology model has low confidence at 3' boundary for feature that is or matches a CDS", # long description
                   0, 0, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf3lcn", "feature",
                   "INDEFINITE_ANNOTATION_END", # short description
                   "alignment to homology model has low confidence at 3' boundary for feature that does not match a CDS", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf3plg", "feature",
                   "INDEFINITE_ANNOTATION_END", # short description
                   "protein-based alignment extends past nucleotide-based alignment at 3' end", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indf3pst", "feature",
                   "INDEFINITE_ANNOTATION_END", # short description
                   "protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint", # long description
                   0, 1, 0, 1, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "indfstrp", "feature",
                   "INDEFINITE_STRAND", # short description
                   "strand mismatch between protein-based and nucleotide-based predictions", # long description
                   0, 1, 0, 0, "indfstr_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "insertnp", "feature",
                   "INSERTION_OF_NT", # short description
                   "too large of an insertion in protein-based alignment", # long description
                   0, 1, 0, 0, "insertn_exc", "coords-value", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "insertnn", "feature",
                   "INSERTION_OF_NT", # short description
                   "too large of an insertion in nucleotide-based alignment of CDS feature", # long description
                   0, 0, 0, 0, "insertn_exc", "coords-value", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "deletinp", "feature",
                   "DELETION_OF_NT", # short description
                   "too large of a deletion in protein-based alignment", # long description
                   0, 1, 0, 1, "deletin_exc", "coords-value", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "deletinn", "feature",
                   "DELETION_OF_NT", # short description
                   "too large of a deletion in nucleotide-based alignment of CDS feature", # long description
                   0, 0, 0, 1, "deletin_exc", "coords-value", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "deletinf", "feature",
                   "DELETION_OF_FEATURE_SECTION", # short description
                   "internal deletion of complete section in multi-section feature with other section(s) annotated", # long description
                   0, 1, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim5c", "feature",
                   "LOW_FEATURE_SIMILARITY_START", # short description
                   "region overlapping annotated feature that is or matches a CDS at 5' end of sequence lacks significant similarity", # long description
                   0, 0, 0, 0, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim5n", "feature",
                   "LOW_FEATURE_SIMILARITY_START", # short description
                   "region overlapping annotated feature that does not match a CDS at 5' end of sequence lacks significant similarity", # long description
                   0, 1, 0, 1, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim5l", "feature",
                   "LOW_FEATURE_SIMILARITY_START", # short description
                   "long region overlapping annotated feature that does not match a CDS at 5' end of sequence lacks significant similarity", # long description
                   0, 1, 0, 0, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim3c", "feature",
                   "LOW_FEATURE_SIMILARITY_END", # short description
                   "region overlapping annotated feature that is or matches a CDS at 3' end of sequence lacks significant similarity", # long description
                   0, 0, 0, 0, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim3n", "feature",
                   "LOW_FEATURE_SIMILARITY_END", # short description
                   "region overlapping annotated feature that does not match a CDS at 3' end of sequence lacks significant similarity", # long description
                   0, 1, 0, 1, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsim3l", "feature",
                   "LOW_FEATURE_SIMILARITY_END", # short description
                   "long region overlapping annotated feature that does not match a CDS at 3' end of sequence lacks significant similarity", # long description
                   0, 1, 0, 0, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsimic", "feature",
                   "LOW_FEATURE_SIMILARITY", # short description
                   "region overlapping annotated feature that is or matches a CDS lacks significant similarity", # long description
                   0, 0, 0, 0, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsimin", "feature",
                   "LOW_FEATURE_SIMILARITY", # short description
                   "region overlapping annotated feature that does not match a CDS lacks significant similarity", # long description
                   0, 1, 0, 1, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "lowsimil", "feature",
                   "LOW_FEATURE_SIMILARITY", # short description
                   "long region overlapping annotated feature that does not match a CDS lacks significant similarity", # long description
                   0, 1, 0, 0, "lowsim_exc", "coords-only", # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR);

  vdr_AlertInfoAdd($alt_info_HHR, "ambgnt5f", "feature",
                   "AMBIGUITY_AT_FEATURE_START", # short description
                   "first nucleotide of non-CDS feature is an ambiguous nucleotide", # long  description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "ambgnt3f", "feature",
                   "AMBIGUITY_AT_FEATURE_END", # short description
                   "final nucleotide of non-CDS feature is an ambiguous nucleotide", # long  description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "ambgnt5c", "feature",
                   "AMBIGUITY_AT_CDS_START", # short description
                   "first nucleotide of CDS is an ambiguous nucleotide", # long  description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "ambgnt3c", "feature",
                   "AMBIGUITY_AT_CDS_END", # short description
                   "final nucleotide of CDS is an ambiguous nucleotide", # long  description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "ambgcd5c", "feature",
                   "AMBIGUITY_IN_START_CODON", # short description
                   "5' complete CDS starts with canonical nt but includes ambiguous nt in its start codon", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  vdr_AlertInfoAdd($alt_info_HHR, "ambgcd3c", "feature",
                   "AMBIGUITY_IN_STOP_CODON", # short description
                   "3' complete CDS ends with canonical nt but includes ambiguous nt in its stop codon", # long description
                   0, 0, 0, 0, undef, undef, # always_fails, causes_failure, prevents_annot, misc_not_failure, exc_key, exc_type
                   $FH_HR); 

  # define the ftbl_invalid_by values, these are one-sided, any alert code listed in the 
  # 3rd argument invalidates the 2nd argument alert code, but not vice versa
  # 
  # Note: only fatal alerts can be invalidated, and only then by other fatal alerts.
  # This is because alert-invalidation is only relevant to the feature table and the 
  # feature table only outputs fatal alerts. See v-annotate.pl:helper_ftable_process_sequence_alerts().

  # mutendcd is invalidated by cdsstopn, mutendex, mutendns
  vdr_AlertInfoSetFTableInvalidatedBy($alt_info_HHR, "mutendcd", "cdsstopn,mutendex,mutendns", $FH_HR); 

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
#   $alt_info_HHR:      REF to array of hashes of error information, FILLED HERE
#   $code:              the code of the element we are adding
#   $pertype:           the 'per-type' of the element we are adding, "sequence" or "feature"
#   $sdesc:             short description of the alert we are adding
#   $ldesc:             long  description of the alert we are adding
#   $always_fails:      '1' if this alert *always* causes its sequence to FAIL, '0' if not
#   $causes_failure:    '1' if this alert causes its sequence to FAIL by default, '0' if not
#   $prevents_annot:    '1' if this alert prevents its sequence from being annotated, '0' if not
#   $misc_not_failure:  '1' if this alert does not cause failure if ftr is "misc_not_failure" from .minfo file
#   $exc_key:           .minfo file key used for exceptions for this alert if it can handle exceptions, 
#                       undef if not
#   $exc_type:          type of exception string for this alert, must be undef if $exc_key is undef, 
#                       and defined if $exc_key is defined, if defined must be 'coords-only', or 'coords-value'
#   $FH_HR:             REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $alt_info_HHR->{"$code"} already exists
#          if $type ne "feature and ne "sequence"
#          if $type ne "sequence"  and $prevents_annot   == 1
#          if $type eq "sequence"  and $misc_not_failure == 1
#          if $always_fails == 1   and $causes_failure   != 1
#          if $always_fails == 1   and $misc_not_failure == 1
#          if (  defined $exc_key) and (! defined $exc_type) 
#          if (! defined $exc_key) and (  defined $exc_type) 
#          if (defined $exc_type)  and ($exc_type ne "coords-only) and ($exc_type ne "coords-value")
#          if an existing code with same exc_key has different exc_type
#
######################p###########################################
sub vdr_AlertInfoAdd { 
  my $sub_name = "vdr_AlertInfoAdd";
  my $nargs_expected = 12;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $code, $pertype, $sdesc, $ldesc, $always_fails, $causes_failure, $prevents_annot, $misc_not_failure, $exc_key, $exc_type, $FH_HR) = (@_);

  # make sure $pertype is valid
  if(($pertype ne "feature") && ($pertype ne "sequence")) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code with per-type $pertype that is not neither \"feature\" nor \"sequence\".", 1, $FH_HR); 
  }
  
  # make sure $always_fails is valid
  if((! defined $always_fails) || (($always_fails != 0) && ($always_fails != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but always_fails is undefined or not 0 or 1", 1, $FH_HR);
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

  # make sure $misc_not_failure is valid, makes sense with $always_fail, and is only true if type is feature
  if((! defined $misc_not_failure) || (($misc_not_failure != 0) && ($misc_not_failure != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but misc_not_failure is undefined or not 0 or 1", 1, $FH_HR);
  }
  if($always_fails && $misc_not_failure) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but always_fails is 1 and misc_not_failure is 1", 1, $FH_HR);
  }
  if(($misc_not_failure == 1) && ($pertype ne "feature")) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but misc_not_failure is 1 and pertype is sequence", 1, $FH_HR);
  }

  # make sure if exc_key is defined, exc_type is too, and vice versa
  if((  defined $exc_key) && (! defined $exc_type)) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but exc_key is defined and exc_key is not defined", 1, $FH_HR);
  }
  if((! defined $exc_key) && (  defined $exc_type)) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but exc_key is not defined and exc_key is defined", 1, $FH_HR);
  }
  # make sure exc_type is valid
  if((defined $exc_type) && (($exc_type ne "coords-only") && ($exc_type ne "coords-value"))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but exc_type ($exc_type) is invalid, it must be either 'coords-only' or 'coords-value'", 1, $FH_HR);
  }
  # check that all exisitng codes with same $exc_key (if any) are of type $exc_type 
  if(defined $exc_key) { 
    foreach my $existing_code (sort keys (%{$alt_info_HHR})) { 
      if((defined $alt_info_HHR->{$existing_code}{"exc_key"}) && 
         ($alt_info_HHR->{$existing_code}{"exc_key"}  eq $exc_key) && 
         ($alt_info_HHR->{$existing_code}{"exc_type"} ne $exc_type)) { 
        ofile_FAIL("ERROR in $sub_name, trying to add code $code with exc_key ($exc_key) and exc_type ($exc_type) but code $existing_code has same exc_key but different exc_type (" . $alt_info_HHR->{$existing_code}{"exc_type"} . ")", 1, $FH_HR);
      }
    }
  }
         
  # check if $code already exists
  if(defined $alt_info_HHR->{$code}) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code, but it already exists in the error info hash", 1, $FH_HR);
  }

  # find highest "order" index so far in the 2D hash
  my $order = scalar(keys (%{$alt_info_HHR}));
  
  # initialize
  %{$alt_info_HHR->{$code}} = ();

  $alt_info_HHR->{$code}{"order"}            = $order; # first "order" value will be '0', second will be '1', etc.
  $alt_info_HHR->{$code}{"pertype"}          = $pertype;
  $alt_info_HHR->{$code}{"sdesc"}            = $sdesc;
  $alt_info_HHR->{$code}{"ldesc"}            = $ldesc;
  $alt_info_HHR->{$code}{"always_fails"}     = $always_fails;
  $alt_info_HHR->{$code}{"causes_failure"}   = $causes_failure;
  $alt_info_HHR->{$code}{"prevents_annot"}   = $prevents_annot;
  $alt_info_HHR->{$code}{"misc_not_failure"} = $misc_not_failure;
  $alt_info_HHR->{$code}{"exc_key"}          = (defined $exc_key)  ? $exc_key  : undef;
  $alt_info_HHR->{$code}{"exc_type"}         = (defined $exc_type) ? $exc_type : undef;
  $alt_info_HHR->{$code}{"ftbl_invalid_by"}  = ""; # initialized to no invalid_by's, possibly added to later with setFTableInvalidatedByErrorInfoHash()

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
# Dies:    - if one of the error codes in $code1 or $code2str do not
#            exist in %{$alt_info_HHR}
#          - if $alt_info_HHR->{$code1}{"ftbl_invalid_by"} is already set
#
#################################################################
sub vdr_AlertInfoSetFTableInvalidatedBy {
  my $sub_name = "vdr_AlertInfoSetFTableInvalidatedBy";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $code1, $code2str, $FH_HR) = (@_);

  if((defined $alt_info_HHR->{$code1}{"ftbl_invalid_by"}) && 
     ($alt_info_HHR->{$code1}{"ftbl_invalid_by"} ne "")) { 
    ofile_FAIL("ERROR in $sub_name, trying to add invalidated by relationship for code $code1, but it is already set", 1, $FH_HR);
  }

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
#   $code:         the code of the element we are setting "causes_failure" value for
#   $value:        value we are setting $ftbl_info_HHR->{$code}{"causes_failure"} to 
#   $FH_HR:        REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $code does not exist in %{$alt_info_HHR}
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
# Subroutine: vdr_AlertInfoSetMiscNotFailure
# Incept:     EPN, Tue Feb  9 13:50:16 2021
#
# Purpose:    Set the "misc_not_failure" value for %{$ftr_info_HHR->{$code}
#
# Arguments:
#   $alt_info_HHR: REF to hash of hashes of error information, FILLED HERE
#   $code:         the code of the element we are setting "misc_not_failure" value for
#   $value:        value we are setting $ftbl_info_HHR->{$code}{"misc_not_failure"} to 
#   $FH_HR:        REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $code does not exist in %{$alt_info_HHR}
#          if $value is not '1' or '0'
#
#################################################################
sub vdr_AlertInfoSetMiscNotFailure {
  my $sub_name = "vdr_AlertInfoSetMiscNotFailure";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HHR, $code, $value, $FH_HR) = (@_);

  if(! defined $alt_info_HHR->{$code}) { 
    ofile_FAIL("ERROR in $sub_name, trying to set misc_not_failure for invalid code $code", 1, $FH_HR);
  }
  if(($value ne "1") && ($value ne "0")) { 
    ofile_FAIL("ERROR in $sub_name, trying to set misc_not_failure to invalid value $value (must be 1 or 0)", 1, $FH_HR);
  }

  $alt_info_HHR->{$code}{"misc_not_failure"} = $value;

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
  my $w_misc    = length("misc not");
  my $w_invalid = utl_HHMaxLengthValueGiven2DKey($alt_info_HHR, "ftbl_invalid_by");
  $w_invalid = utl_Max($w_invalid, length("invalidated"));

  # determine order of codes to print
  my @code_A = ();
  foreach my $code (sort keys (%{$alt_info_HHR})) { 
    $code_A[($alt_info_HHR->{$code}{"order"})] = $code;
  }

  printf $FH ("%-4s  %8s  %-*s  %*s  %*s  %*s  %*s  %-*s  %-*s\n", 
              "#", "", $w_sdesc, "short", $w_afails, "fails", $w_fails, "", $w_annot, "prevents", $w_misc, "misc not", $w_invalid, "invalidated", $w_ldesc, "long");
  printf $FH ("%-4s  %8s  %-*s  %*s  %*s  %*s  %*s  %-*s  %-*s\n", 
              "#idx", "code", $w_sdesc, "desc", $w_afails, "always", $w_fails, "fails", $w_annot, "annot", $w_misc, "failure", $w_invalid, "by", $w_ldesc, "desc");
  printf $FH ("%-4s  %8s  %-*s  %*s  %*s  %*s  %*s  %-*s  %-*s\n", 
              "#---", "--------", 
              $w_sdesc,   utl_StringMonoChar($w_sdesc,   "-", undef), 
              $w_afails,  utl_StringMonoChar($w_afails,  "-", undef), 
              $w_fails,   utl_StringMonoChar($w_fails,   "-", undef), 
              $w_annot,   utl_StringMonoChar($w_annot,   "-", undef), 
              $w_misc,    utl_StringMonoChar($w_misc,    "-", undef), 
              $w_invalid, utl_StringMonoChar($w_invalid, "-", undef), 
              $w_ldesc,   utl_StringMonoChar($w_ldesc,   "-", undef));

  my $idx = 0;
  foreach my $code (@code_A) { 
    $idx++;
    printf $FH ("%-4s  %5s  %-*s  %*s  %*s  %*s  %*s  %*s  %-*s\n", 
                $idx, $code, 
                $w_sdesc,   helper_tabular_replace_spaces($alt_info_HHR->{$code}{"sdesc"}), 
                $w_afails,  $alt_info_HHR->{$code}{"always_fails"}     ? "yes" : "no",
                $w_fails,   $alt_info_HHR->{$code}{"causes_failure"}   ? "yes" : "no",
                $w_annot,   $alt_info_HHR->{$code}{"prevents_annot"}   ? "yes" : "no",
                $w_misc,    $alt_info_HHR->{$code}{"misc_not_failure"} ? "yes" : "no",
                $w_invalid, $alt_info_HHR->{$code}{"ftbl_invalid_by"}, 
                $w_ldesc,   $alt_info_HHR->{$code}{"ldesc"});
  }

  return;
}

#################################################################
# Subroutine: vdr_FeatureAlertCausesFailure()
# Incept:     EPN, Fri Feb  5 11:35:58 2021
#
# Purpose:    Does $alt_code for feature $ftr_idx cause failure?
#             Considers both "causes_failure" and "misc_not_failure"
#             fields.
#
# Arguments: 
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $alt_info_HHR:  REF to the alert info hash of arrays, PRE-FILLED
#  $ftr_idx:       feature index
#  $alt_code:      alert code
#
# Returns:    1 or 0
#
# Dies:       never; does not validate anything.
#
################################################################# 
sub vdr_FeatureAlertCausesFailure { 
  my $sub_name = "vdr_FeatureAlertCausesFailure";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $alt_info_HHR, $ftr_idx, $alt_code) = @_;

  return (($alt_info_HHR->{$alt_code}{"causes_failure"}) && 
          (! vdr_FeatureAlertIsMiscNotFailure($ftr_info_AHR, $alt_info_HHR, $ftr_idx, $alt_code))) 
      ? 1 : 0;
}       

#################################################################
# Subroutine: vdr_FeatureAlertIsMiscNotFailure()
# Incept:     EPN, Fri Feb  5 11:41:11 2021
#
# Purpose:    Does $alt_code for feature $ftr_idx have the 
#             "misc_not_failure" attribute?
#
# Arguments: 
#  $ftr_info_AHR:  ref to the feature info array of hashes 
#  $alt_info_HHR:  REF to the alert info hash of arrays, PRE-FILLED
#  $ftr_idx:       feature index
#  $alt_code:      alert code
#
# Returns:    1 or 0
#
# Dies:       never; does not validate anything.
#
################################################################# 
sub vdr_FeatureAlertIsMiscNotFailure { 
  my $sub_name = "vdr_FeatureAlertIsMiscNotFailure";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $alt_info_HHR, $ftr_idx, $alt_code) = @_;

  return (($ftr_info_AHR->[$ftr_idx]{"misc_not_failure"}) && 
          ($alt_info_HHR->{$alt_code}{"misc_not_failure"})) 
      ? 1 : 0; 
}

#################################################################
# Subroutine: vdr_ExceptionCoordsAndValuesToSegmentsAndValues()
# Incept:     EPN, Thu Sep 14 14:55:00 2023
#
# Purpose:    Given a coords plus value string like we expect 
#             in exception strings of the type 'coords-value'.
#             Update a provided hash with keys equal to each segment
#             and values equal to the parsed value.
#
#             Example:
#             coords_value: "11..13:+:36,40..27:-:23";
#             $ret_sgm_value_HR->{"11..13:+"} = 36;
#             $ret_sgm_value_HR->{"40..27:-"} = 23;
#
# Arguments: 
#  $coords_value:      the coords segment plus value
#  $extra_errmsg:      upon failure, extra string to add to error message (can be undef)
#  $ret_sgm_value_HR:  ref to hash to update
#  $FH_HR:             ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $coords_value (invalid format)
#       if any of the segments overlap 
# 
################################################################# 
sub vdr_ExceptionCoordsAndValuesToSegmentsAndValues {
  my $sub_name = "vdr_ExceptionCoordsAndValuesToSegmentsAndValues";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords_value, $extra_errmsg, $ret_sgm_value_HR, $FH_HR) = @_;
  if(! defined $coords_value) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, coords_value is undefined.%s", (defined $extra_errmsg) ? ("\n" . $extra_errmsg) : ""), 1, $FH_HR);
  }
  my @sgm_value_A = split(",", $coords_value);
  foreach my $sgm_value (@sgm_value_A) { 
    my ($start, $stop, $strand, $value) = vdr_ExceptionSegmentAndValueParse($sgm_value, $extra_errmsg, $FH_HR); # this will die if format is invalid
    $ret_sgm_value_HR->{(vdr_CoordsSegmentCreate($start, $stop, $strand, $FH_HR))} = $value;
  }

  # check for any overlap segments, die if so
  my $noverlap;
  foreach my $key1 (sort keys (%{$ret_sgm_value_HR})) { 
    foreach my $key2 (sort keys (%{$ret_sgm_value_HR})) { 
      if($key1 ne $key2) { 
        ($noverlap, undef) = vdr_CoordsSegmentOverlap($key1, $key2, $FH_HR);
        if($noverlap > 0) { 
          ofile_FAIL(sprintf("ERROR in $sub_name, two coords segments overlap in exception coords: $key1 and $key2.%s", (defined $extra_errmsg) ? ("\n" . $extra_errmsg) : ""), 1, $FH_HR);
        }
      }
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_ExceptionSegmentAndValueParse()
# Incept:     EPN, Thu Sep 14 14:31:48 2023
#
# Purpose:    Given a single coords token plus a value in format
#             return its start, stop, strand and value.
#
# Arguments: 
#  $coords_value:  the coords segment plus value
#  $extra_errmsg:  upon failure, extra string to add to error message (can be undef)
#  $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    4 values:
#             $start:  start position
#             $stop:   stop position
#             $strand: strand
#             $value:  value
#
# Dies: if unable to parse $coords_value (invalid format)
#
################################################################# 
sub vdr_ExceptionSegmentAndValueParse {
  my $sub_name = "vdr_ExceptionSegmentAndValueParse";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords_value, $extra_errmsg, $FH_HR) = @_;
  if(! defined $coords_value) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, coords_value is undefined.%s", (defined $extra_errmsg) ? ("\n" . $extra_errmsg) : ""), 1, $FH_HR);
  }
  if($coords_value =~ /^\<?(\d+)\.\.\>?(\d+)\:([\+\-]):(\S+)$/) { 
    return ($1, $2, $3, $4);
  }
  ofile_FAIL(sprintf("ERROR in $sub_name, unable to parse coords-value token $coords_value%s", (defined $extra_errmsg) ? ("\n" . $extra_errmsg) : ""), 1, $FH_HR); 

  return; # NEVER REACHED
}

#################################################################
# Subroutine: vdr_ExceptionSegmentsAndValuesToPositionsAndValues()
# Incept:     EPN, Thu Sep 14 15:20:48 2023
#
# Purpose:    Given a hash with keys of coords segments, fill a 
#             different hash with per position keys, one for each
#             position of the coords range.
#
#             Example:
#             $coords_value_HR->{"11..13:+"} = "36";
#             $coords_value_HR->{"17..17:+"} = "42";
#
#             $ret_posn_value_HR->{"11"} = "36";
#             $ret_posn_value_HR->{"12"} = "36";
#             $ret_posn_value_HR->{"13"} = "36";
#             $ret_posn_value_HR->{"17"} = "42";
#
# Arguments: 
#  $coords_value_HR:   the coords segment plus value
#  $allow_updates:     1 if we can update value for an already set key in %{$ret_posn_value_HR}
#                      0 if we should die if we try to update a value for an already set key
#  $ret_posn_value_HR: ref to hash to update
#  $FH_HR:             ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $coords_value (invalid format)
#       if we try to set a value for a key in $ret_posn_value_HR->{key} that is already defined
#
################################################################# 
sub vdr_ExceptionSegmentsAndValuesToPositionsAndValues {
  my $sub_name = "vdr_ExceptionSegmentsAndValuesToPositionsAndValues";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords_value_HR, $allow_updates, $ret_posn_value_HR, $FH_HR) = @_;
  foreach my $coords_key (sort keys (%{$coords_value_HR})) { 
    vdr_ExceptionOneSegmentAndValueToPositionsAndValues($coords_key, $coords_value_HR->{$coords_key}, $allow_updates, $ret_posn_value_HR, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: vdr_ExceptionOneSegmentAndValueToPositionsAndValues()
# Incept:     EPN, Thu Sep 14 14:38:42 2023
#
# Purpose:    Given a single coords token plus a value 
#             update a provided hash with key equal to a position
#             and value.
#
#             Example:
#             coords: "11..13:+";
#             value:  "36"
#             $ret_posn_value_HR->{"11"} = 36;
#             $ret_posn_value_HR->{"12"} = 36;
#             $ret_posn_value_HR->{"13"} = 36;
#
# Arguments: 
#  $coords:            the coords segment
#  $value:             the value
#  $allow_updates:     1 if we can update value for an already set key in %{$ret_posn_value_HR}
#                      0 if we should die if we try to update a value for an already set key
#  $ret_posn_value_HR: ref to hash to update
#  $FH_HR:             ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies: if unable to parse $coords_value (invalid format)
#       if we try to set a value for a key in $ret_posn_value_HR->{key} that is already defined
#
################################################################# 
sub vdr_ExceptionOneSegmentAndValueToPositionsAndValues {
  my $sub_name = "vdr_ExceptionOneSegmentAndValueToPositionsAndValues";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $value, $allow_updates, $ret_posn_value_HR, $FH_HR) = @_;
  if(! defined $coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", 1, $FH_HR);
  }
  if(! vdr_CoordsSegmentValidate($coords, $FH_HR)) { 
    ofile_FAIL("ERROR in $sub_name, coords is in incorrect format", 1, $FH_HR);
  }
  my ($start, $stop, $strand) = vdr_CoordsSegmentParse($coords, $FH_HR);
  if($start >= $stop) { 
    my $tmp = $start;
    $start = $stop;
    $stop  = $tmp;
  }
  for(my $posn = $start; $posn <= $stop; $posn++) { 
    if((! $allow_updates) && (defined $ret_posn_value_HR->{$posn})) { 
      ofile_FAIL("ERROR in $sub_name, trying to update already set key ($posn) and allow_updates is false.\nThis could be because two coords segments in an exception string overlap.", 1, $FH_HR); 
    }
    $ret_posn_value_HR->{$posn} = $value;
  }

  return;
}


#################################################################
# Subroutine:  vdr_ParseQsubFile()
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
# Subroutine: vdr_SubmitJobAsScript()
# Incept:      EPN, Thu Jan 28 17:19:51 2021
#
# Purpose:     Creates a script that will execute a job and sumits a job to sge that will 
#              execute that script.
#
# Arguments:
#   $cmd:            command to run
#   $qsub_prefix:    qsub command prefix to use when submitting to farm, undef if running locally
#   $qsub_suffix:    qsub command suffix to use when submitting to farm, undef if running locally
#   $job_name:       name for job
#   $sh_file:        name of shell script file to create with command to run
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
sub vdr_SubmitJobAsScript {
  my $sub_name = "vdr_SubmitJobAsScript()";
  my $nargs_expected = 10;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $qsub_prefix, $qsub_suffix, $job_name, $sh_file, $err_file, $mem_gb, $nsecs, $opt_HHR, $ofile_info_HHR) = @_;
  
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  if(($err_file ne "/dev/null") && (-e $err_file)) { 
    utl_FileRemoveUsingSystemRm($err_file, $sub_name, $opt_HHR, $ofile_info_HHR); 
  }

  my $submit_cmd = $qsub_prefix . "sh $sh_file" . $qsub_suffix;
  # replace changeable parts of qsub suffix and prefix
  $submit_cmd =~ s/\!\[jobname\]\!/$job_name/g;
  $submit_cmd =~ s/\!\[errfile\]\!/$err_file/g;
  $submit_cmd =~ s/\!\[memgb\]\!/$mem_gb/g;
  $submit_cmd =~ s/\!\[nsecs\]\!/$nsecs/g;

  # create the shell script file with the cmsearch/cmalign/rRNA_sensor command $cmd
  vdr_WriteCommandScript($sh_file, $cmd, $FH_HR);
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
#          When --errcheck is used, this function considers any output
#          written to stderr output files in @{$errfile_AR} to mean
#          that the corresponding job has 'failed' and should be
#          considered to have finished. When --errcheck is not used we
#          don't look at the err files.
# 
#
# Arguments: 
#  $do_cmalign:       '1' if we're running cmalign, which requires special care because we
#                     handle two cases: finish successfully or die with a specific error
#  $do_errcheck:      '1' if we should fail if any error file is written to
#  $outkey:           key in second dimension of out_file_AHR we'll check to see if job is finished
#  $init_sleep_secs:  number of seconds to wait initially before checking 
#  $input_chunk_secs: number of seconds to increment wait by until we reach 120s
#  $out_file_AHR:     ref to array of hashes of output files that will be created by jobs we are waiting for
#  $success_AR:       ref to array of success values, FILLED HERE, can be undef if ! $do_cmalign
#                     these will always all be '1' unless $do_cmalign
#                     if($do_cmalign) some may be '0'
#  $mxsize_AR:        ref to array of required matrix sizes, CAN BE UNDEF
#                     $mxsize_AR->[$j] set to value readh from cmalign output, if $success_AR->[$j] == 0
#                     else set to '0'
#  $finished_str:     string that indicates a job is finished e.g. "[ok]"
#  $opt_HHR:          ref to options hash
#  $FH_HR:            ref to hash of file handles
#
# Returns:     Number of jobs (<= scalar(@outfile_A)) that have
#              finished.
# 
# Dies: never.
#
################################################################# 
sub vdr_WaitForFarmJobsToFinish { 
  my $sub_name = "vdr_WaitForFarmJobsToFinish()";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($do_cmalign, $do_errcheck, $outkey, $init_sleep_secs, $input_chunk_secs, $out_file_AHR, $success_AR, $mxsize_AR, $finished_str, $opt_HHR, $FH_HR) = @_;
 
  my $log_FH = $FH_HR->{"log"};
  my $nmin = opt_Get("--wait", $opt_HHR);

  # contract check
  if(! exists $out_file_AHR->[0]{$outkey}) { 
    ofile_FAIL("ERROR in $sub_name, no $outkey files in out_file_AHR", 1, $FH_HR);
  }
  if(! exists $out_file_AHR->[0]{"err"}) { 
    ofile_FAIL("ERROR in $sub_name, no err files in out_file_AHR", 1, $FH_HR);
  }

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

  # variables related to wait times between checks
  my $cur_sleep_secs = $init_sleep_secs;  # number of seconds to wait before first check
  my $chunk_secs     = $input_chunk_secs; # number of seconds to increment wait time by until we reach 
                                          # $doubling_secs at which point we'll double until we reach $max_sleep
  my $doubling_secs  = 120;               # number of seconds to wait before starting to double $cur_sleep_secs
  my $max_sleep_secs = 120;               # maximum number of seconds we'll wait between checks
  my $secs_waited    = 0;                 # number of total seconds we've waited thus far

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
    $secs_waited += $chunk_secs;
    if($secs_waited >= $doubling_secs) { 
      $cur_sleep_secs *= 2;
    }
    else { 
      $cur_sleep_secs += $chunk_secs;
    }
    if($cur_sleep_secs > $max_sleep_secs) { # reset to max if we've exceeded it
      $cur_sleep_secs = $max_sleep_secs;
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
# Subroutine: vdr_CoordsToSegments()
# Incept:     EPN, Thu Sep 14 13:53:42 2023
#
# Synopsis: Given a coords string, validate it is in the
#           correct format and fill an array of segments.
#
#           Example: 
#           $coords        = "30..100:+,102..105:+,200..333:+"
#           @{$ret_sgm_AR} = ("30..100:+","102..105:+","200..333:+");
# 
# Arguments:
#  $coords:     coordinate string
#  $ret_sgm_AR: ref to array of coords segment to fill
#  $FH_HR:      ref to hash of file handles, including "log" and "cmd"
#
# Returns:    number of segments, fills @{$ret_sgm_AR} with segments
#
# Dies: if unable to parse $coords (incorrect format)
#
#################################################################
sub vdr_CoordsToSegments {
  my $sub_name = "vdr_CoordsToSegments";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $ret_sgm_AR, $FH_HR) = @_;
  if(! defined $coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", 1, $FH_HR); 
  }
  @{$ret_sgm_AR} = split(",", $coords);
  foreach my $coords_tok (@{$ret_sgm_AR}) { 
    vdr_CoordsSegmentParse($coords_tok, $FH_HR); # this will fail if token is not in correct format
  }

  return scalar(@{$ret_sgm_AR});
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
# Subroutine: vdr_CoordsValidate()
# Incept:     EPN, Fri Sep 15 10:45:06 2023
#
# Synopsis: Return '1' if a given coords string is valid
#           else '0'. Useful if caller wants to check format
#           but not die if format is incorrect, possibly so 
#           it can report an informative error message.
# 
# Arguments:
#  $coords:   coordinate string
#  $FH_HR:    REF to has of file handles
#
# Returns:    '1' if coordinate string is in the correct format
#             '0' if not
#
# Dies: never
#
#################################################################
sub vdr_CoordsValidate {
  my $sub_name = "vdr_CoordsValidate";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $FH_HR) = @_;
  if(! defined $coords) { 
    ofile_FAIL("ERROR in $sub_name, coords is undefined", 1, $FH_HR); 
  }

  my @coords_A  = split(",", $coords);
  foreach my $coords_tok (@coords_A) {
    if(! vdr_CoordsSegmentValidate($coords_tok, $FH_HR)) { 
      return 0;
    }
  }

  return 1;
}

#################################################################
# Subroutine: vdr_CoordsSegmentValidate()
# Incept:     EPN, Fri Sep 15 10:49:19 2023
#
# Synopsis: Return '1' if a given coords segment string is valid
#           else '0'. Useful if caller wants to check format
#           but not die if format is incorrect, possibly so 
#           it can report an informative error message.
# 
# Arguments:
#  $coords_tok: coordinate string
#  $FH_HR:      REF to has of file handles
#
# Returns:    '1' if coordinate string is in the correct format
#             '0' if not
#
# Dies: never
#
#################################################################
sub vdr_CoordsSegmentValidate {
  my $sub_name = "vdr_CoordsSegmentValidate";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords_tok, $FH_HR) = @_;
  if(! defined $coords_tok) { 
    ofile_FAIL("ERROR in $sub_name, coords_tok is undefined", 1, $FH_HR); 
  }

  if($coords_tok =~ /^\<?(\d+)\.\.\>?(\d+)\:([\+\-])$/) { 
    return 1;
  }

  return 0;
}

#################################################################
# Subroutine: vdr_CoordsSinglePositionSegmentCreate()
# Incept:     EPN, Wed May 26 07:32:43 2021
#
# Synopsis: Create a length 1 coords token from a given position
#           and strand. Removes any carrots in $posn.
# 
# Arguments:
#  $posn:     start and stop position
#  $strand:   strand ("+" or "-")
#  $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    coordinate token <posn>..<posn>:<strand>
#
# Dies:  if $posn is invalid
#        if $strand is not "+" or "-"
#
#################################################################
sub vdr_CoordsSinglePositionSegmentCreate {
  my $sub_name = "vdr_CoordsSinglePositionSegmentCreate";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($posn, $strand, $FH_HR) = @_;
  if(($posn !~ /^\<?(\d+)$/) && ($posn !~ /^\>?(\d+)$/)) { ofile_FAIL("ERROR in $sub_name, posn is invalid ($posn)", 1, $FH_HR); }
  if(($strand ne "+")        && ($strand ne "-"))        { ofile_FAIL("ERROR in $sub_name, strand is invalid ($strand)", 1, $FH_HR); }
  $posn =~ s/^\>//;
  $posn =~ s/^\<//;

  return $posn . ".." . $posn . ":" . $strand;
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
# Returns:  Two values:
#           $noverlap:    Number of nucleotides of overlap between <$coords_tok1>
#                         and <$coords_tok2> on the same strand, 0 if none
#           $overlap_reg: region of overlap, "" if none
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
    return (0, "");
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
#           within nucleotide sequence with absolute coords <$abs_coords>.
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
# Subroutine: vdr_CoordsRelativeSingleCoordToAbsolute()
#
# Incept:     EPN, Wed Apr 29 06:27:45 2020
#
# Synopsis: Return absolute nucleotide coordinate that corresponds to
#           a single relative nucleotide coordinate <$rel_coord>.
#           within nucleotide sequence with absolute coords <$abs_coords>.
#           Single relative coordinate $rel_coord strand is irrelevant
#           because it is a single position.
#
#           Examples:
# 
#           abs_coords     rel_coord   returns
#           "11..100:+"    "6"         "16"
#           "11..100:+"    "38"        "48"     
#
# Arguments:
#  $abs_coords:     nucleotide coordinates in full sequence [1..seqlen]
#  $rel_coord:      relative nucleotide coordinate
#  $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:   Absolute coordinate corresponding to $rel_coord.
#
# Dies: if $rel_coord is longer than total length of absolute coords
#
#################################################################
sub vdr_CoordsRelativeSingleCoordToAbsolute { 
  my $sub_name = "vdr_CoordsRelativeSingleCoordToAbsolute";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($abs_coords, $rel_coord, $FH_HR) = @_;

  my $ret_coords = vdr_CoordsRelativeSegmentToAbsolute($abs_coords, vdr_CoordsSegmentCreate($rel_coord, $rel_coord, "+", $FH_HR), $FH_HR);
  if($ret_coords =~ /^(\d+)\.\.\d+/) { 
    return $1;
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, problem converting relative coordinate $rel_coord to absolute coords (abs_coords: $abs_coords)", 1, $FH_HR);
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine: vdr_CoordsProteinRelativeToAbsolute()
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
# Subroutine: vdr_CoordsFromStartStopStrandArrays()
#
# Incept:     EPN, Wed Sep  8 18:32:20 2021
#
# Synopsis: Given references to arrays of start, stop and strand
#           values for 1 or more segments, combine them to 
#           create a coords string with 1 or more segments.
#
# Arguments:
#  $start_AR:  reference to array of start positions
#  $stop_AR:   reference to array of stop positions
#  $strand_AR: reference to array of strand values
#  $FH_HR:     REF to hash of file handles, including "log" and "cmd"
#
# Returns:  coords string with all segments combined
#
# Dies: If unable to parse any segment 
#
#################################################################
sub vdr_CoordsFromStartStopStrandArrays {
  my $sub_name = "vdr_CoordsFromStartStopStrandArrays";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($start_AR, $stop_AR, $strand_AR, $FH_HR) = @_;

  my $coords = "";
  my $nsgm = scalar(@{$start_AR});
  if($nsgm != scalar(@{$stop_AR}))   { ofile_FAIL("ERROR in $sub_name, different number of start and stop elements",   1, $FH_HR); }
  if($nsgm != scalar(@{$strand_AR})) { ofile_FAIL("ERROR in $sub_name, different number of start and strand elements", 1, $FH_HR); }

  for(my $s = 0; $s < $nsgm; $s++) { 
    if($s > 0) { 
      $coords .= ","; 
    }
    $coords .= vdr_CoordsSegmentCreate($start_AR->[$s], $stop_AR->[$s], $strand_AR->[$s], $FH_HR);
  }

  return $coords;
}

#################################################################
# Subroutine: vdr_CoordsSegmentActualToFractional()
# Incept:     EPN, Tue Jan  4 09:00:48 2022
#
# Purpose:    Given two coords segments, one of a larger region and 
#             one of smaller region within that larger region, calculate
#             the fractional position of start and stop of the smaller
#             region within the larger one.
#
# Arguments:
#  $full_coords:     coordinates of full region
#  $subseq_coords:   coordinates of subsequence
#  $FH_HR:           ref to hash of file handles
#
# Returns:  2 values:
#           $fract_start: fractional start position [0.0 to 1.0], undef if $full_coords and subseq_coords are not same strand or do not overlap
#           $fract_stop:  fractional start position [0.0 to 1.0], undef if $full_coords and subseq_coords are not same strand or do not overlap
#           Always true that $fract_start <= $fract_stop
#
# Dies: If $full_coords or $subseq_coords is not a single parseable coords segment
#       If $strand from $full_coords or $subseq_coords is not "+" or "-"
#       If $strand from $full_coords or $subseq_coords is "+" and start > stop
#       If $strand from $full_coords or $subseq_coords is "-" and start < stop
#
#################################################################
sub vdr_CoordsSegmentActualToFractional { 
  my $sub_name = "vdr_CoordsActualToFractional";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($full_coords, $subseq_coords, $FH_HR) = (@_);

  # printf("in $sub_name, full_coords: $full_coords subseq_coords: $subseq_coords\n");

  my ($full_start,   $full_stop,   $full_strand)   = vdr_CoordsSegmentParse($full_coords, $FH_HR);
  my ($subseq_start, $subseq_stop, $subseq_strand) = vdr_CoordsSegmentParse($subseq_coords, $FH_HR);

  my $ret_fract_start = undef;
  my $ret_fract_stop  = undef;

  my $small_value = 0.0000000001;

  if(($full_strand ne "+") && ($full_strand ne "-")) { 
    ofile_FAIL("ERROR in $sub_name, full_coords $full_coords strand is not + or -", 1, $FH_HR);
  }
  if(($subseq_strand ne "+") && ($subseq_strand ne "-")) { 
    ofile_FAIL("ERROR in $sub_name, subseq_coords $subseq_coords strand is not + or -", 1, $FH_HR);
  }
  if(($full_strand eq "+") && ($full_start > $full_stop)) { 
    ofile_FAIL("ERROR in $sub_name, full_coords $full_coords + strand but start > stop", 1, $FH_HR);
  }
  if(($full_strand eq "-") && ($full_start < $full_stop)) { 
    ofile_FAIL("ERROR in $sub_name, full_coords $full_coords - strand but start < stop", 1, $FH_HR);
  }
  if(($subseq_strand eq "+") && ($subseq_start > $subseq_stop)) { 
    ofile_FAIL("ERROR in $sub_name, subseq_coords $subseq_coords + strand but start > stop", 1, $FH_HR);
  }
  if(($subseq_strand eq "-") && ($subseq_start < $subseq_stop)) { 
    ofile_FAIL("ERROR in $sub_name, subseq_coords $subseq_coords - strand but start < stop", 1, $FH_HR);
  }

  if($full_strand ne $subseq_strand) { 
    return (undef, undef); 
  }

  if($full_strand eq "+") { # positive strand
    if(($full_start > $full_stop) || ($subseq_start > $subseq_stop) || 
       ($subseq_start < $full_start) || ($subseq_start > $full_stop) || 
       ($subseq_stop  < $full_start) || ($subseq_stop  > $full_stop)) { 
      return (undef, undef); 
    }
    my $full_len = $full_stop - $full_start + 1;
    $ret_fract_start = ($subseq_start - $full_start + 1) / $full_len;
    $ret_fract_stop  = ($subseq_stop  - $full_start + 1) / $full_len;
  }
  else { # negative strand
    if(($full_start < $full_stop) || ($subseq_start < $subseq_stop) || 
       ($subseq_start > $full_start) || ($subseq_start < $full_stop) || 
       ($subseq_stop  > $full_start) || ($subseq_stop  < $full_stop)) { 
      return (undef, undef); 
    }
    my $full_len = $full_start - $full_stop + 1;
    $ret_fract_start = ($full_start - $subseq_start + 1) / $full_len;
    $ret_fract_stop  = ($full_start - $subseq_stop  + 1) / $full_len;
  }

  if(($ret_fract_start < (0 - $small_value)) || 
     ($ret_fract_stop  > (1 + $small_value))) { 
      return (undef, undef); 
  }

  # printf("in $sub_name returning %.5f .. %.5f\n", $ret_fract_start, $ret_fract_stop);
  return ($ret_fract_start, $ret_fract_stop);
}

#################################################################
# Subroutine: vdr_CoordsSegmentFractionalToActual()
# Incept:     EPN, Tue Jan  4 11:08:23 2022
#
# Purpose:    Given a coords segment and two fractional positions [0..1]
#             within it estimate the actual positions those fractional 
#             positions correspond to.
#
# Arguments:
#  $coords:   coordinates of full region (1 segment)
#  $fstart:   fractional start position
#  $fstop:    fractional stop  position
#  $FH_HR:    ref to hash of file handles
#
# Returns:  3 values:
#           $start:  start position (rounded)
#           $stop:   stop position (rounded)
#           $strand: strand 
#
# Dies: If $coords is not a single parseable coords segment
#       If $strand from $coords is not "+" or "-"
#       If $strand from $coords is "+" and start > stop
#       If $strand from $coords is "-" and start < stop
#################################################################
sub vdr_CoordsSegmentFractionalToActual {
  my $sub_name = "vdr_CoordsFractionalToActual";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords, $fstart, $fstop, $FH_HR) = (@_);

  # printf("in $sub_name, coords: $coords fstart: $fstart fstop: $fstop\n");

  my ($full_start, $full_stop, $full_strand) = vdr_CoordsSegmentParse($coords, $FH_HR);

  if(($full_strand ne "+") && ($full_strand ne "-")) { 
    ofile_FAIL("ERROR in $sub_name, coords $coords strand is not + or -", 1, $FH_HR);
  }
  if(($full_strand eq "+") && ($full_start > $full_stop)) { 
    ofile_FAIL("ERROR in $sub_name, coords $coords + strand but start > stop", 1, $FH_HR);
  }
  if(($full_strand eq "-") && ($full_start < $full_stop)) { 
    ofile_FAIL("ERROR in $sub_name, coords $coords - strand but start < stop", 1, $FH_HR);
  }

  my $ret_start = undef;
  my $ret_stop  = undef;

  my $full_len = abs($full_stop - $full_start) + 1;
  if($full_strand eq "+") { # positive strand 
    $ret_start = int($full_start + ($fstart * $full_len) - 1);
    if($ret_start < $full_start) { $ret_start = $full_start; }

    $ret_stop = int($full_start + ($fstop * $full_len) - 1);
    if($ret_stop > $full_stop) { $ret_stop = $full_stop; }

  }
  else { # negative strand
    $ret_start = int($full_start - ($fstart * $full_len) + 1);
    if($ret_start > $full_start) { $ret_start = $full_start; }

    $ret_stop = int($full_start - ($fstop * $full_len) + 1);
    if($ret_stop < $full_stop) { $ret_stop = $full_stop; }
  }

  # printf("in $sub_name returning $ret_start..$ret_stop:$full_strand\n");
  return($ret_start, $ret_stop, $full_strand);
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
        if($mdl_name =~ /[\)\(]/) { 
          ofile_FAIL("ERROR in $sub_name, model info file has model named $mdl_name which contains '(' and/or ')', which are not allowed in model names", 1, $FH_HR);
        }
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
# Subroutine: vdr_ModelInfoValidateExceptionsKeys()
# Incept:     EPN, Tue Sep 12 13:15:22 2023
#
# Purpose:    Validate any model info keys that pertain to
#             alert exceptions.
# 
#             Uses 'exc_key' value in alt_info_HH
#             to validate all mdl_info_H keys that end in 
#             '_exc', and dies if any such keys are not valid.
#             
# Arguments: 
#  $mdl_info_HR:    ref to the model info hash (for one model)
#  $alt_info_HHR: ref to the alert info hash of hashes
#  $FH_HR:          ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $mdl_info_HR->{"*_exc"} = $exc_str is not valid 
#             (if no alt_info_HH{<code>}{"exc_key"} eq $exc_str)
#
################################################################# 
sub vdr_ModelInfoValidateExceptionKeys {
  my $sub_name = "vdr_ModelInfoValidateExceptionKeys";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_info_HR, $alt_info_HHR, $FH_HR) = @_;

  # small optimization, make temporary hash of all $alt_info_HH{<code>}{"exc_key"}
  my %tmp_key_H = ();
  foreach my $code (sort keys %{$alt_info_HHR}) { 
    if(defined $alt_info_HHR->{$code}{"exc_key"}) { 
      $tmp_key_H{$alt_info_HHR->{$code}{"exc_key"}} = $alt_info_HHR->{$code}{"exc_type"};
    }
  }
 
  foreach my $exc_key (sort keys %{$mdl_info_HR}) { 
    if($exc_key =~ /^.+\_exc$/) { 
      if(! defined $tmp_key_H{$exc_key}) { 
        ofile_FAIL("ERROR, in $sub_name, no alert codes use the exception key $exc_key read in model info file", 1, $FH_HR);
      }
      # now validate format of value
      if($tmp_key_H{$exc_key} eq "coords-only") { 
        if(! vdr_CoordsValidate($mdl_info_HR->{$exc_key}, $FH_HR)) { 
          ofile_FAIL(sprintf("ERROR, in $sub_name, invalid format for exception key: $exc_key, type 'coords-only', read %s.\nBut expected coords string, with spans separated by commas, e.g. \"11..15:+\" or \"11..15:+,18..18:+\".", $mdl_info_HR->{$exc_key}), 1, $FH_HR);
        }
      }
      else { # coords-value
        # validate it by parsing it
        my %sgm_value_H = (); # we won't use this but need to pass it to vdr_ExceptionCoordsAndValuesToSegmentsAndValues()
        my $errmsg = sprintf("ERROR, in $sub_name, invalid format for exception key: $exc_key, type 'coords-value', read %s.\nBut expected coords-value string, with span separated by commas, e.g. \"11..15:+:23\" or \"11..15:+:23,18..18:+:36\".", $mdl_info_HR->{$exc_key});
        vdr_ExceptionCoordsAndValuesToSegmentsAndValues($mdl_info_HR->{$exc_key}, $errmsg, \%sgm_value_H, $FH_HR);
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
#              a cmalign run finished successfully,
#              in error, or has not yet finished.
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
    ofile_FAIL("ERROR in $sub_name, cmalign $stdout_file exists but is empty. v-annotate.pl may have run out of available memory, especially if you see a 'Killed' message.", 1, $FH_HR);
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
    if($error_line =~ /Error: .+ mxes need (\d+\.\d+)/) { 
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
# Subroutine:  vdr_CmalignParseInsertFile()
# Incept:      EPN, Thu Jan 31 13:06:54 2019
#
# Purpose:    Parse Infernal 1.1 cmalign --ifile output and store
#             results in %{$seq_inserts_HHR}.
#
#             
# Arguments: 
#  $ifile_file:       ifile file to parse
#  $seq_inserts_HHR:  REF to hash of hashes with insert information, added to here, must be defined
#                     key 1: sequence name
#                     key 2: one of 'spos', 'epos', 'ins'
#                     $seq_inserts_HHR->{}{"spos"} is starting model position of alignment
#                     $seq_inserts_HHR->{}{"epos"} is ending model position of alignment
#                     $seq_inserts_HHR->{}{"ins"} is the insert string in the format:
#                     <mdlpos_1>:<uapos_1>:<inslen_1>;...<mdlpos_n>:<uapos_n>:<inslen_n>;
#                     for n inserts, where insert x is defined by:
#                     <mdlpos_x> is model position after which insert occurs 0..mdl_len (0=before first pos)
#                     <uapos_x> is unaligned sequence position of the first aligned nt
#                     <inslen_x> is length of the insert
#  $mdl_name_HR:      REF to hash of model name hash to fill, added to here, can be undef
#                     key: sequence name, value: name of model this sequence was aligned to
#  $seq_name_AR:      REF to array of sequence names, in order read, added to here, can be undef
#  $seq_len_HR:       REF to hash of sequence lengths to fill, added to here, can be undef
#  $mdl_len_HR:       REF to hash of model name hash to fill, added to here, can be undef
#                     key: *model* name, value: length of model
#  $FH_HR:            REF to hash of file handles
#
# Returns:    void
#
# Dies:       if unable to parse the ifile
#             if $seq_inserts_HHR is undef upon entry
#             if the same model exists twice in the ifile with
#             different lengths and $mdl_len_HR is defined
#
################################################################# 
sub vdr_CmalignParseInsertFile { 
  my $sub_name = "vdr_CmalignParseInsertFile()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ifile_file, $seq_inserts_HHR, $mdl_name_HR, $seq_name_AR, $seq_len_HR, $mdl_len_HR, $FH_HR) = @_;
  
  open(IN, $ifile_file) || ofile_FileOpenFailure($ifile_file, $sub_name, $!, "reading", $FH_HR);

  if(! defined $seq_inserts_HHR) {
    ofile_FAIL("ERROR in $sub_name, seq_inserts_HHR is undef upon entry", 1, $FH_HR);
  }
    
  my $line_ctr = 0;  # counts lines in ifile_file
  my $mdl_name = undef;
  my $mdl_len  = undef;
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
        ($mdl_name, $mdl_len) = (@el_A);
        if(defined $mdl_len_HR) {
          if((defined $mdl_len_HR->{$mdl_name}) &&
             ($mdl_len_HR->{$mdl_name} != $mdl_len)) {
            ofile_FAIL("ERROR in $sub_name, model $mdl_name appears multiple times with different lengths $mdl_len and " . $mdl_len_HR->{$mdl_name}, 1, $FH_HR);
          }
          $mdl_len_HR->{$mdl_name} = $mdl_len;
        }
      }
      elsif(scalar(@el_A) >= 4) { 
        if((! defined $mdl_name) || (! defined $mdl_len)) { 
          ofile_FAIL("ERROR in $sub_name, read sequence line for sequence " . $el_A[0] . "before model line", 1, $FH_HR);
        }
        my $nel = scalar(@el_A); 
        if((($nel - 4) % 3) != 0) { # check number of elements makes sense
          ofile_FAIL("ERROR in $sub_name, unexpected number of elements ($nel) in ifile line in $ifile_file on line $line_ctr:\n$line\n", 1, $FH_HR);
        }          
        my ($seq_name, $seq_len, $spos, $epos) = ($el_A[0], $el_A[1], $el_A[2], $el_A[3]);
        if(! defined $seq_inserts_HHR->{$seq_name}) { 
          # initialize
          %{$seq_inserts_HHR->{$seq_name}} = ();
          if(defined $mdl_name_HR) { $mdl_name_HR->{$seq_name} = $mdl_name; }
          if(defined $seq_name_AR) { push(@{$seq_name_AR}, $seq_name); }
          if(defined $seq_len_HR)  { $seq_len_HR->{$seq_name} = $seq_len; }
        }
        # create the insert string
        my $insert_str = "";
        for(my $el_idx = 4; $el_idx < scalar(@el_A); $el_idx += 3) { 
          $insert_str .= $el_A[$el_idx] . ":" . $el_A[$el_idx+1] . ":" . $el_A[$el_idx+2] . ";"; 
        }
        $seq_inserts_HHR->{$seq_name}{"spos"} = $spos;
        $seq_inserts_HHR->{$seq_name}{"epos"} = $epos;
        $seq_inserts_HHR->{$seq_name}{"ins"}  = $insert_str;
      }
    }
  }
  close(IN);
  
  return;
}

#################################################################
# Subroutine:  vdr_CmalignWriteInsertFile()
# Incept:      EPN, Fri Apr  3 11:17:49 2020
#
# Purpose:    Write an Infernal 1.1 cmalign --ifile given
#             insert information in %{$seq_inserts_HHR}.
#
# Arguments: 
#  $ifile_file:       ifile file to write to
#  $do_append:        '1' to append to if file exists, '0' to create new file
#  $mdl_name:         name of model for model line, all seqs should have insert info relative to this model
#  $mdl_len:          length of model for model line,
#  $seq_name_AR:      REF to array of sequence names to output data for from seq_inserts_HHR
#  $seq_len_HR:       REF to hash of sequence lengths
#  $seq_inserts_HHR:  REF to hash of hashes with insert information, already filled
#                     key 1: sequence name
#                     key 2: one of 'spos', 'epos', 'ins'
#                     $seq_inserts_HHR->{}{"spos"} is starting model position of alignment
#                     $seq_inserts_HHR->{}{"epos"} is ending model position of alignment
#                     $seq_inserts_HHR->{}{"ins"} is the insert string in the format:
#                     <mdlpos_1>:<uapos_1>:<inslen_1>;...<mdlpos_n>:<uapos_n>:<inslen_n>;
#                     for n inserts, where insert x is defined by:
#                     <mdlpos_x> is model position after which insert occurs 0..mdl_len (0=before first pos)
#                     <uapos_x> is unaligned sequence position of the first aligned nt
#                     <inslen_x> is length of the insert
#  $FH_HR:            REF to hash of file handles
#
# Returns:    void
#
# Dies:       if unable to parse the ifile
#
################################################################# 
sub vdr_CmalignWriteInsertFile { 
  my $sub_name = "vdr_CmalignWriteInsertFile()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ifile_file, $do_append, $mdl_name, $mdl_len, $seq_name_AR, $seq_len_HR, $seq_inserts_HHR, $FH_HR) = @_;

  # print("in $sub_name\n");
  # utl_HHDump("seq_inserts_HH", $seq_inserts_HHR, *STDOUT);
  
  my $out_mode = (($do_append) && (-s $ifile_file)) ? ">>" : ">";
  open(OUT, $out_mode, $ifile_file) || ofile_FileOpenFailure($ifile_file, $sub_name, $!, (($out_mode eq ">") ? "writing" : "appending"), $FH_HR);

  my $line_ctr = 0;  # counts lines in ifile_file

  # model line: 2 tokens
  #<modelname> <modellen>
  #norovirus.NC_039477 7567
  print OUT $mdl_name . " " . $mdl_len . "\n";

  # per-sequence lines, each has at least 3 tokens:
  # <seqname> <spos> <epos>
  # and optionally 
  foreach my $seq_name (@{$seq_name_AR}) {
    if(defined $seq_inserts_HHR->{$seq_name}) { 
      print OUT ($seq_name . " " . $seq_len_HR->{$seq_name} . " " . $seq_inserts_HHR->{$seq_name}{"spos"} . " " . $seq_inserts_HHR->{$seq_name}{"epos"});
      if((defined $seq_inserts_HHR->{$seq_name}{"ins"}) &&
         $seq_inserts_HHR->{$seq_name}{"ins"} ne "") {
        my @ins_A = split(";", $seq_inserts_HHR->{$seq_name}{"ins"});
        foreach my $ins (@ins_A) {
          # example line:
          #gi|669176088|gb|KM198574.1| 7431 17 7447  2560 2539 3  2583 2565 3
          if($ins =~ /^(\d+)\:(\d+)\:(\d+)/) { 
            print OUT ("  " . $1 . " " . $2 . " " . $3);
          }
          else {
            ofile_FAIL("ERROR in $sub_name, unable to parse insert string " . $seq_inserts_HHR->{$seq_name}{"ins"} . " at token $ins", 1, $FH_HR);
          }
        }
      }
      print OUT "\n";
    }
  }
  print OUT "//\n";

  close OUT;
}

#################################################################
# Subroutine:  vdr_CmemitConsensusToFile()
# Incept:      EPN, Wed Apr 15 09:31:54 2020
#
# Purpose:    Run cmemit -c for model $mdl_name fetched from $cm_file
#             and return the consensus sequence.
#
# Arguments: 
#  $execs_HR:        REF to a hash with "blastx" and "parse_blastx.pl""
#  $cm_file:         CM file to fetch from
#  $mdl_name:        name of CM file to fetch
#  $cseq_fa_file:    name of output fasta file to create
#  $opt_HHR:         REF to options 2D hash
#  $ofile_info_HHR:  REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    the consensus sequence as a string
#
# Dies:       If cmemit fails or we can't read any sequence from 
#             the output fasta file
#
################################################################# 
sub vdr_CmemitConsensus {
  my $sub_name = "vdr_CmemitConsensusToFile";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $cm_file, $mdl_name, $cseq_fa_file, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $cmd = $execs_HR->{"cmfetch"} . " $cm_file $mdl_name | " . $execs_HR->{"cmemit"} . " -c - > $cseq_fa_file";
  utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), opt_Get("--keep", $opt_HHR), $FH_HR);

  # fetch the sequence
  my @file_lines_A = ();
  utl_FileLinesToArray($cseq_fa_file, 1, \@file_lines_A, $FH_HR);
  my $ret_cseq = "";
  my $nlines = scalar(@file_lines_A);
  if($nlines <= 1) { 
    ofile_FAIL("ERROR in $sub_name, read 0 seq data from $cseq_fa_file", 1, $FH_HR); 
  }
  for(my $i = 1; $i < $nlines; $i++) { # start at $i == 1, skip header line
    my $seq_line = $file_lines_A[$i];
    chomp $seq_line;
    $ret_cseq .= $seq_line;
  }
  
  return $ret_cseq;
}

#################################################################
# Subroutine:  vdr_MergeOutputConcatenateOnly()
# Incept:      EPN, Fri Mar 19 09:19:11 2021
#
# Purpose:    With --split, merge output files from multiple output
#             directories in @{$chunk_outdir_AR} into a single file.
#
# Arguments: 
#   $out_root_no_vadr:  root name for output file names, without '.vadr' suffix
#   $out_sfx:           output name suffix
#   $ofile_key:         key for %{$ofile_info_HHR}
#   $ofile_desc:        description for %{$ofile_info_HHR}, "" to not add the file to %{$ofile_info_HHR}
#   $do_check_exists:   '1' to check if all files to merge exist before concatenating and fail if not
#   $chunk_outdir_AR:   ref to array of output directories with files we are merging
#   $opt_HHR:           ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $ofile_info_HHR:    ref to the 2D hash of output file information, ADDED TO HERE 
#
# Returns:     name of the merged file created
# 
# Dies: if $check_exists is 1 and a file to merge does not exist
# 
################################################################# 
sub vdr_MergeOutputConcatenateOnly { 
  my $nargs_exp = 8;
  my $sub_name = "vdr_MergeOutputConcatenateOnly";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($out_root_no_vadr, $out_sfx, $ofile_key, $ofile_desc, $do_check_exists, $chunk_outdir_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my @filelist_A = (); # array of files to concatenate to make $merged_file
  vdr_MergeOutputGetFileList($out_root_no_vadr, $out_sfx, $do_check_exists, \@filelist_A, $chunk_outdir_AR, $FH_HR);

  my $merged_file = $out_root_no_vadr . ".vadr" . $out_sfx; # merged file to create by concatenating files in chunk dirs
  if(scalar(@filelist_A) > 0) { 
    utl_ConcatenateListOfFiles(\@filelist_A, $merged_file, $sub_name, $opt_HHR, $FH_HR);
    if($ofile_desc ne "") { 
      ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $ofile_key, $merged_file, 1, 1, $ofile_desc);
    }
  }
  elsif($do_check_exists) { 
    ofile_FAIL("ERROR in $sub_name, zero files from chunk dir to concatenate to make $merged_file", 1, $FH_HR);
  }
  return $merged_file;
}

#################################################################
# Subroutine:  vdr_MergeOutputConcatenatePreserveSpacing()
# Incept:      EPN, Sat May 22 09:16:03 2021
#
# Purpose:    With --split, merge output files from multiple output
#             directories in @{$chunk_outdir_AR} into a single file
#             and preserve spacing on <$np> at least consecutive
#             lines.
#
# Arguments: 
#   $out_root_no_vadr:  root name for output file names, without '.vadr' suffix
#   $out_sfx:           output name suffix
#   $ofile_key:         key for %{$ofile_info_HHR}
#   $ofile_desc:        description for %{$ofile_info_HHR}, "" to not add the file to %{$ofile_info_HHR}
#   $do_check_exists:   '1' to check if all files to merge exist before concatenating and fail if not
#   $np:                number of lines to preserve spacing for, -1 to preserve spacing for all lines
#   $csep:              column separator string, often "  "
#   $empty_flag:        '1' to output empty lines for empty arrays of data, '0'
#                       to output empty lines as header separation lines
#   $head_AAR:          header values
#   $cljust_AR:         ref to '1'/'0' array of indicating if a column is left justified or not
#   $chunk_outdir_AR:   ref to array of output directories with files we are merging
#   $opt_HHR:           ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $ofile_info_HHR:    ref to the 2D hash of output file information
#
# Returns:     name of the merged file created
# 
# Dies: if $check_exists is 1 and a file to merge does not exist
# 
################################################################# 
sub vdr_MergeOutputConcatenatePreserveSpacing { 
  my $nargs_exp = 13;
  my $sub_name = "vdr_MergeOutputConcatenatePreserveSpacing";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($out_root_no_vadr, $out_sfx, $ofile_key, $ofile_desc, $do_check_exists, $np, $csep, $empty_flag, $head_AAR, $cljust_AR, $chunk_outdir_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $ncol  = scalar(@{$cljust_AR});
  
  my @filelist_A = (); # array of files to concatenate to make $merged_file
  vdr_MergeOutputGetFileList($out_root_no_vadr, $out_sfx, $do_check_exists, \@filelist_A, $chunk_outdir_AR, $FH_HR);

  my $merged_file = $out_root_no_vadr . ".vadr" . $out_sfx; # merged file to create by concatenating files in chunk dirs
  my $out_FH = undef; # output file handle
  open($out_FH, ">", $merged_file) || ofile_FileOpenFailure($merged_file, $sub_name, $!, "writing", $FH_HR);

  my $line; 
  my $nline     = 0; # line number for current file
  my $nline_tot = 0; # line number over all files
  my $ncol2print = $ncol; # updated below by adding 1 if nec
  my @data_AA = (); # [0..$i..$nline-1][0..$j..$ncol2print-1], data read from input to output again with correct spacing
  my $seen_noncomment_line = 0; # used to skip header lines
  my $nout = 0; # number of times we output a chunk
  my $j;
  
  if(scalar(@filelist_A) > 0) { 
    foreach my $file (@filelist_A) {
      $seen_noncomment_line = 0; # used to skip header lines
      open(IN, $file) || ofile_FileOpenFailure($file, $sub_name, $!, "reading", $FH_HR);
      while($line = <IN>) { 
        chomp $line;
        my @el_A = split(/\s+/, $line);
        my $nel = scalar(@el_A);
        if($line !~ m/^\#/) { # a non-comment line
          $seen_noncomment_line = 1;
          @{$data_AA[$nline]} = ();
          for($j = 0; $j < ($ncol-1); $j++) {
            push(@{$data_AA[$nline]}, $el_A[$j]);
          }
          # combine all columns after $ncol into one, separated by single space
          if($nel >= $ncol) {
            my $toadd = "";
            $ncol2print = $ncol + 1; 
            for($j = ($ncol-1); $j < ($nel-1); $j++) {
              $toadd .= $el_A[$j] . " ";
            }
            $toadd .= $el_A[($nel-1)];
            push(@{$data_AA[$nline]}, $toadd);
          }
          $nline++;
          $nline_tot++;
        }
        else { # a comment-line
          if($seen_noncomment_line || $line eq "#") { 
            push(@data_AA, []);  # push empty array --> blank line 
            $nline++;
            $nline_tot++;
          }
        }
      }
      if($nline >= $np) {
        ofile_TableHumanOutput(\@data_AA, $head_AAR, $cljust_AR, undef, undef, $csep, undef, undef, undef, undef, $empty_flag, $out_FH, undef, $FH_HR);
        undef @data_AA;
        @data_AA = ();
        $nline = 0;
        $nout++;
      }
    }
    # output remaining lines
    if(($nline > 0) || ($nout == 0)) {
      ofile_TableHumanOutput(\@data_AA, $head_AAR, $cljust_AR, undef, undef, $csep, undef, undef, undef, undef, $empty_flag, $out_FH, undef, $FH_HR);
    }
    close($out_FH);
    if($ofile_desc ne "") { 
      ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $ofile_key, $merged_file, 1, 1, $ofile_desc);
    }
  }
  elsif($do_check_exists) { 
    ofile_FAIL("ERROR in $sub_name, zero files from chunk dir to concatenate to make $merged_file", 1, $FH_HR);
  }
  
  return $merged_file;
}

#################################################################
# Subroutine:  vdr_MergeOutputGetFileList()
# Incept:      EPN, Fri Mar 19 16:08:22 2021
#
# Purpose:    With --split, fill an array with a list of files
#             to merge.
#
# Arguments: 
#   $out_root_no_vadr:  root name for output file names, without '.vadr' suffix
#   $out_sfx:           output name suffix
#   $do_check_exists:   '1' to check if all files to merge exist before concatenating and fail if not
#   $filelist_AR:       list of files to merge   
#   $chunk_outdir_AR:   ref to array of output directories with files we are merging
#   $FH_HR:             ref to hash of file handles, including "cmd"
#
# Returns:     void
# 
# Dies: if $check_exists is 1 and a file to merge does not exist
# 
################################################################# 
sub vdr_MergeOutputGetFileList {
  my $nargs_exp = 6;
  my $sub_name = "vdr_MergeOutputGetFileList";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($out_root_no_vadr, $out_sfx, $do_check_exists, $filelist_AR, $chunk_outdir_AR, $FH_HR) = @_;

  @{$filelist_AR} = ();

  # make list of files to concatenate
  my $nchunk = scalar(@{$chunk_outdir_AR});
  my $out_dir_tail = utl_RemoveDirPath($out_root_no_vadr);
  for(my $i = 1; $i <= $nchunk; $i++) { 
    my $chunk_file = $chunk_outdir_AR->[($i-1)] . "/" . $out_dir_tail . "." . $i . ".vadr" . $out_sfx;
    if(-e $chunk_file) { 
      push(@{$filelist_AR}, $chunk_file);
    }
    elsif($do_check_exists) { # file does not exist
      ofile_FAIL("ERROR in $sub_name, expected file to concatenate $chunk_file does not exist", 1, $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine:  vdr_MergeOutputMdlTabularFile()
# Incept:      EPN, Fri Mar 19 13:27:00 2021
#
# Purpose:    With --split, merge .mdl tabular output files from 
#             multiple output directories in @{$chunk_outdir_AR} 
#             into a single file.
#
# Arguments: 
#   $out_root_no_vadr:  root name for output file names, without '.vadr' suffix
#   $ofile_desc:        description for %{$ofile_info_HHR}
#   $do_check_exists:   '1' to check if all files to merge exist before concatenating and fail if not
#   $chunk_outdir_AR:   ref to array of output directories with files we are merging
#   $opt_HHR:           ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $ofile_info_HHR:    ref to the 2D hash of output file information, ADDED TO HERE 
#
# Returns:     void
# 
# Dies: if $check_exists is 1 and a file to merge does not exist
# 
################################################################# 
sub vdr_MergeOutputMdlTabularFile { 
  my $nargs_exp = 6;
  my $sub_name = "vdr_MergeOutputMdlTabularFile";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($out_root_no_vadr, $ofile_desc, $do_check_exists, $chunk_outdir_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $out_sfx   = ".mdl";

  my $out_dir_tail = utl_RemoveDirPath($out_root_no_vadr);

  # make list of files to concatenate
  my @filelist_A = (); # array of files to concatenate to make $merged_file
  vdr_MergeOutputGetFileList($out_root_no_vadr, $out_sfx, $do_check_exists, \@filelist_A, $chunk_outdir_AR, $FH_HR);

  # th head_* definitions should be (manually) kept consistent with output_tabular()
  # alternatively we could parse the header lines in the files we want to merge,
  # but not doing that currently
  my @head_mdl_AA = ();
  my @data_mdl_AA = ();
  @{$head_mdl_AA[0]} = ("",    "",      "",      "",         "num",  "num",  "num");
  @{$head_mdl_AA[1]} = ("idx", "model", "group", "subgroup", "seqs", "pass", "fail");
  my @clj_mdl_A      = (1,     1,       1,       1,          0,      0,      0);

  # read each .mdl file and store info in it
  my ($idx, $model, $group, $subgroup, $num_seqs, $num_pass, $num_fail);
  my %group_H    = (); # key: model name, value: group
  my %subgroup_H = (); # key: model name, value: subgroup
  my %num_seqs_H = (); # key: model name, value: num seqs
  my %num_pass_H = (); # key: model name, value: num passing seqs
  my %num_fail_H = (); # key: model name, value: num failing seqs
  for(my $i = 0; $i < scalar(@filelist_A); $i++) { 
    open(IN, $filelist_A[$i]) || ofile_FileOpenFailure($filelist_A[$i], $sub_name, $!, "reading", $FH_HR);
    while(my $line = <IN>) { 
      ##                                                    num   num   num
      ##idx  model               group         subgroup    seqs  pass  fail
      ##---  ------------------  ------------  ----------  ----  ----  ----
      #1     NC_045512           Sarbecovirus  SARS-CoV-2     2     1     1
      #2     NC_045512-MW422255  Sarbecovirus  SARS-CoV-2     1     1     0
      ##---  ------------------  ------------  ----------  ----  ----  ----
      #-     *all*               -             -              3     2     1
      #-     *none*              -             -              0     0     0
      ##---  ------------------  ------------  ----------  ----  ----  ----
      if($line !~ m/^\#/) { 
        chomp $line;
        if($line =~ m/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)$/) { 
          ($idx, $model, $group, $subgroup, $num_seqs, $num_pass, $num_fail) = ($1, $2, $3, $4, $5, $6, $7);
        }
        else { 
          ofile_FAIL("ERROR in $sub_name unable to parse $filelist_A[$i] file line:\n$line\n", 1, $FH_HR);
        }
        if(! defined $group_H{$model}) { 
          $group_H{$model} = $group;
          # initialize other hashes for this model
          $num_seqs_H{$model} = 0;
          $num_pass_H{$model} = 0;
          $num_fail_H{$model} = 0;
        }
        elsif($group_H{$model} ne $group) { 
          ofile_FAIL("ERROR in $sub_name read more than one distinct group for model $model: $group_H{$model} and $group", 1, $FH_HR);
        }
        if(! defined $subgroup_H{$model}) { 
          $subgroup_H{$model} = $subgroup;
        }
        elsif($subgroup_H{$model} ne $subgroup) { 
          ofile_FAIL("ERROR in $sub_name read more than one distinct group for model $model: $group_H{$model} and $group", 1, $FH_HR);
        }
        $num_seqs_H{$model} += $num_seqs;
        $num_pass_H{$model} += $num_pass;
        $num_fail_H{$model} += $num_fail;
      }
    }
  }

  my @tmp_mdl_tbl_order_A = (sort { $num_seqs_H{$b} <=> $num_seqs_H{$a} or 
                                        $a cmp $b 
                             } keys (%num_seqs_H));

  # remove special "*all*" and "*none*" lines from @tmp_mdl_order_A to make @mdl_order_A
  my @mdl_tbl_order_A = ();
  foreach $model (@tmp_mdl_tbl_order_A) { 
    if(($model ne "*all*") && ($model ne "*none*")) { 
      push(@mdl_tbl_order_A, $model);
    }
  }

  # push data to @data_mdl_AA
  # the following block should be (manually) kept consistent with output_tabular()
  my $mdl_tbl_idx = 0;
  foreach $model (@mdl_tbl_order_A) { 
    if($num_seqs_H{$model} > 0) { 
      $mdl_tbl_idx++;
      push(@data_mdl_AA, [$mdl_tbl_idx, $model, $group_H{$model}, $subgroup_H{$model}, $num_seqs_H{$model}, $num_pass_H{$model}, $num_fail_H{$model}]);
    }
  }
  # add mdl summary line
  push(@data_mdl_AA, []); # separator line
  $model = "*all*";
  push(@data_mdl_AA, ["-", $model, $group_H{$model}, $subgroup_H{$model}, $num_seqs_H{$model}, $num_pass_H{$model}, $num_fail_H{$model}]);
  $model = "*none*";
  push(@data_mdl_AA, ["-", $model, $group_H{$model}, $subgroup_H{$model}, $num_seqs_H{$model}, $num_pass_H{$model}, $num_fail_H{$model}]);
  push(@data_mdl_AA, []); # separator line

  my $merged_file = $out_root_no_vadr . ".vadr" . $out_sfx; # merged file to create by concatenating files in chunk dirs
  ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, "mdl", $merged_file, 1, 1, "per-model tabular summary file");
  ofile_TableHumanOutput(\@data_mdl_AA, \@head_mdl_AA, \@clj_mdl_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"mdl"}, undef, $FH_HR);

  return;
}

#################################################################
# Subroutine:  vdr_MergeOutputAlcTabularFile()
# Incept:      EPN, Fri Mar 19 16:19:19 2021
#
# Purpose:    With --split, merge .alc tabular output files from 
#             multiple output directories in @{$chunk_outdir_AR} 
#             into a single file.
#
# Arguments: 
#   $out_root_no_vadr:  root name for output file names, without '.vadr' suffix
#   $alt_info_HHR:      ref to the alert info hash of arrays, PRE-FILLED
#   $ofile_desc:        description for %{$ofile_info_HHR}
#   $do_check_exists:   '1' to check if all files to merge exist before concatenating and fail if not
#   $chunk_outdir_AR:   ref to array of output directories with files we are merging
#   $opt_HHR:           ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $ofile_info_HHR:    ref to the 2D hash of output file information, ADDED TO HERE 
#
# Returns:     '1' if zero alerts were reported, else '0'
# 
# Dies: if $check_exists is 1 and a file to merge does not exist
# 
################################################################# 
sub vdr_MergeOutputAlcTabularFile { 
  my $nargs_exp = 7;
  my $sub_name = "vdr_MergeOutputAlcTabularFile";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($out_root_no_vadr, $alt_info_HHR, $ofile_desc, $do_check_exists, $chunk_outdir_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $out_sfx   = ".alc";

  my $out_dir_tail = utl_RemoveDirPath($out_root_no_vadr);

  # make list of files to concatenate
  my @filelist_A = (); # array of files to concatenate to make $merged_file
  vdr_MergeOutputGetFileList($out_root_no_vadr, $out_sfx, $do_check_exists, \@filelist_A, $chunk_outdir_AR, $FH_HR);

  # th head_* definitions should be (manually) kept consistent with output_tabular()
  # alternatively we could parse the header lines in the files we want to merge,
  # but not doing that currently
  my @head_alc_AA = ();
  my @data_alc_AA = ();
  @{$head_alc_AA[0]} = ("",    "alert",  "causes",  "short",       "per",  "num",   "num",  "long");
  @{$head_alc_AA[1]} = ("idx", "code",   "failure", "description", "type", "cases", "seqs", "description");
  my @clj_alc_A      = (1,     1,        1,          1,            0,      0,      0,        1);

  # read each .alc file and store info in it
  my ($idx, $alt_code, $causes_failure, $short_description, $per_type, $num_cases, $num_seqs, $long_description);

  my %data_HH = (); # key 1: alert code, key 2: column name, value: value from column
  my @invariant_keys_A = ("causes_failure", "short_description", "per_type", "long_description");
  my @sum_keys_A       = ("num_cases", "num_seqs"); 
  my $key;
  my %line_H = ();
  for(my $i = 0; $i < scalar(@filelist_A); $i++) { 
    open(IN, $filelist_A[$i]) || ofile_FileOpenFailure($filelist_A[$i], $sub_name, $!, "reading", $FH_HR);
    while(my $line = <IN>) { 
       ##---  --------  -------  ---------------------------  --------  -----  ----  -----------
       ##     alert     causes   short                       per    num   num  long       
       ##idx  code      failure  description                type  cases  seqs  description
       ##---  --------  -------  ----------------------  -------  -----  ----  -----------
       #1     ambgnt5f  no       N_AT_FEATURE_START      feature      3     1  first nucleotide of non-CDS feature is an N
       #2     ambgnt3f  no       N_AT_FEATURE_END        feature      3     1  final nucleotide of non-CDS feature is an N
       #3     ambgnt5c  no       N_AT_CDS_START          feature      2     1  first nucleotide of CDS is an N
       #4     ambgnt3c  no       N_AT_CDS_END            feature      2     1  final nucleotide of CDS is an N
       #---  --------  -------  ----------------------  -------  -----  ----  -----------
       #5     unexleng  yes*     UNEXPECTED_LENGTH       feature      1     1  length of complete coding (CDS or mat_peptide) feature is not a multiple of 3
       #6     cdsstopn  yes*     CDS_HAS_STOP_CODON      feature      1     1  in-frame stop codon exists 5' of stop position predicted by homology to reference
       #7     fstukcfi  yes*     POSSIBLE_FRAMESHIFT     feature      1     1  possible frameshift in CDS (internal)
       #8     indfantn  yes      INDEFINITE_ANNOTATION   feature      1     1  nucleotide-based search identifies CDS not identified in protein-based search
       #9     lowsimic  yes      LOW_FEATURE_SIMILARITY  feature      2     1  region within annotated feature that is or matches a CDS lacks significant similarity
       ##---  --------  -------  ----------------------  -------  -----  ----  -----------
      if($line !~ m/^\#/) { 
        chomp $line;
        if($line =~ m/^(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(.+)$/) { 
          ($idx, $alt_code, 
           $line_H{"causes_failure"}, 
           $line_H{"short_description"}, 
           $line_H{"per_type"}, 
           $line_H{"num_cases"}, 
           $line_H{"num_seqs"}, 
           $line_H{"long_description"}) = ($1, $2, $3, $4, $5, $6, $7, $8);
        }
        else { 
          ofile_FAIL("ERROR in $sub_name unable to parse $filelist_A[$i] file line:\n$line\n", 1, $FH_HR);
        }
        if(! defined $data_HH{$alt_code}) { 
          %{$data_HH{$alt_code}}= ();
          # init counters
          foreach $key (@sum_keys_A) { 
            $data_HH{$alt_code}{$key} = 0;
          }
        }
        foreach $key (@invariant_keys_A) { 
          if(! defined $data_HH{$alt_code}{$key}) { 
            $data_HH{$alt_code}{$key} = $line_H{$key};
          }
          elsif($data_HH{$alt_code}{$key} ne $line_H{$key}) { 
            ofile_FAIL("ERROR in $sub_name read more than one distinct value for $key for alert code $alt_code: $data_HH{$alt_code}{$key} and $line_H{$key}", 1, $FH_HR);
          }
        }
        foreach $key (@sum_keys_A) { 
          $data_HH{$alt_code}{$key} += $line_H{$key};
        }
      }
    }
  }

  # determine order of alert codes to print
  my @alt_code_A = (); # all alerts in output order
  alert_order_arrays($alt_info_HHR, \@alt_code_A, undef, undef);

  # push data to @data_mdl_AA
  # the following block should be (manually) kept consistent with output_tabular()
  my $alt_idx = 0;
  my $zero_alerts = 1; # set to '0' below if we have >= 1 alerts
  my $alc_sep_flag = 0;
  foreach my $alt_code (@alt_code_A) { 
    if(defined $data_HH{$alt_code}) { 
      if($data_HH{$alt_code}{"num_cases"} > 0) { 
        if(! $alt_info_HHR->{$alt_code}{"causes_failure"}) { 
          $alc_sep_flag = 1; 
        }
        if(($alt_info_HHR->{$alt_code}{"causes_failure"}) && $alc_sep_flag) { 
          # print separation line between alerts that cause and do not cause failure
          push(@data_alc_AA, []); # separator line
          $alc_sep_flag = 0; 
        }
        $alt_idx++;
        # don't need to check if alert is 'misc_not_feature' already have '*' added to 'causes_failure' column values
        # don't need to check if alert is 'causes_failure' already know that from tables we parsed
        push(@data_alc_AA, [$alt_idx, $alt_code, 
                            $data_HH{$alt_code}{"causes_failure"},
                            $data_HH{$alt_code}{"short_description"},
                            $data_HH{$alt_code}{"per_type"}, 
                            $data_HH{$alt_code}{"num_cases"}, 
                            $data_HH{$alt_code}{"num_seqs"}, 
                            $data_HH{$alt_code}{"long_description"}]);
        $zero_alerts = 0;
      }
    }
  }
  if(! $zero_alerts) { 
    push(@data_alc_AA, []); # separator line
  }

  my $merged_file = $out_root_no_vadr . ".vadr" . $out_sfx; # merged file to create by concatenating files in chunk dirs
  ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, "alc", $merged_file, 1, 1, "alert count tabular summary file");
  ofile_TableHumanOutput(\@data_alc_AA, \@head_alc_AA, \@clj_alc_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"alc"}, undef, $FH_HR);

  return $zero_alerts;
}

#################################################################
# Subroutine:  vdr_MergePerFeatureFastaFiles()
# Incept:      EPN, Mon Mar 22 06:28:21 2021
#
# Purpose:    With --out_allfasta or --keept merge per-feature fasta
#             files for each model.
#
# Arguments: 
#   $out_root_no_vadr:  root name for output file names, without '.vadr' suffix
#   $mdl_info_AHR:      ref to the model info hash of arrays, PRE-FILLED
#   $ftr_info_HAHR:     ref to hash of array of hashes with info on features per model, PRE-FILLED
#   $chunk_outdir_AR:   ref to array of output directories with files we are merging
#   $opt_HHR:           ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $ofile_info_HHR:    ref to the 2D hash of output file information, ADDED TO HERE 
#
# Returns:     void
# 
# Dies: if problem concatenating files
# 
################################################################# 
sub vdr_MergePerFeatureFastaFiles { 
  my $nargs_exp = 6;
  my $sub_name = "vdr_MergePerFeatureFastaFiles";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($out_root_no_vadr, $mdl_info_AHR, $ftr_info_HAHR, $chunk_outdir_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nmdl = scalar(@{$mdl_info_AHR});
  for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    my $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    my $nftr = scalar(@{$ftr_info_HAHR->{$mdl_name}});
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      my $ftr_out_sfx    = "." . $mdl_name . "." . vdr_FeatureTypeAndTypeIndexString($ftr_info_HAHR->{$mdl_name}, $ftr_idx, ".") . ".fa";
      my $ftr_ofile_key  = $mdl_name . ".pfa." . $ftr_idx;
      my $ftr_ofile_desc = "model " . $mdl_name . " feature " . vdr_FeatureTypeAndTypeIndexString($ftr_info_HAHR->{$mdl_name}, $ftr_idx, "#") . " predicted seqs";
      vdr_MergeOutputConcatenateOnly($out_root_no_vadr, $ftr_out_sfx, $ftr_ofile_key, $ftr_ofile_desc, 0, $chunk_outdir_AR, $opt_HHR, $ofile_info_HHR);
    }
  }
  return;
}

#################################################################
# Subroutine:  vdr_MergeFrameshiftStockholmFiles()
# Incept:      EPN, Mon Mar 22 10:56:39 2021
#
# Purpose:    With --out_fsstk or --keep merge per-segment Stockholm
#             files for each model.
#
# Arguments: 
#   $out_root_no_vadr:  root name for output file names, without '.vadr' suffix
#   $mdl_info_AHR:      ref to the model info hash of arrays, PRE-FILLED
#   $ftr_info_HAHR:     ref to hash of array of hashes with info on features per model, PRE-FILLED
#   $sgm_info_HAHR:     ref to hash of array of hashes with info on segments per model, PRE-FILLED
#   $chunk_outdir_AR:   ref to array of output directories with files we are merging
#   $opt_HHR:           ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $ofile_info_HHR:    ref to the 2D hash of output file information, ADDED TO HERE 
#
# Returns:     void
# 
# Dies: if problem concatenating files
# 
################################################################# 
sub vdr_MergeFrameshiftStockholmFiles { 
  my $nargs_exp = 7;
  my $sub_name = "vdr_MergeFrameshiftStockholmFiles";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($out_root_no_vadr, $mdl_info_AHR, $ftr_info_HAHR, $sgm_info_HAHR, $chunk_outdir_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nmdl = scalar(@{$mdl_info_AHR});
  for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    my $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    my $ftr_info_AHR = $ftr_info_HAHR->{$mdl_name}; # for convenience
    my $nftr = scalar(@{$ftr_info_AHR});
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) { 
        for(my $sgm_idx = $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}; $sgm_idx <= $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"}; $sgm_idx++) { 
          my $cds_and_sgm_idx = vdr_FeatureTypeAndTypeIndexString($ftr_info_AHR, $ftr_idx, ".") . "." . ($sgm_idx - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} + 1);
          my $sgm_out_sfx    = "." . $mdl_name . "." . $cds_and_sgm_idx . ".frameshift.stk";
          my $sgm_ofile_key  = $mdl_name . "." . $cds_and_sgm_idx . ".frameshift.stk";
          my $sgm_ofile_desc = "Stockholm file for >= 1 possible frameshifts for $cds_and_sgm_idx for model $mdl_name";
          vdr_MergeOutputConcatenateOnly($out_root_no_vadr, $sgm_out_sfx, $sgm_ofile_key, $sgm_ofile_desc, 0, $chunk_outdir_AR, $opt_HHR, $ofile_info_HHR);
        }
      }
    }
  }

  return;
}

#################################################################
# Subroutine:  vdr_MergeAlignments()
# Incept:      EPN, Mon Mar 22 07:26:57 2021
#
# Purpose:    With --out_*{stk,afa} or --keep merge alignment 
#             files for each model. If aligned fasta, we need to
#             convert to pfam first. 
#
# Arguments: 
#   $out_root_no_vadr:  root name for output file names, without '.vadr' suffix
#   $execs_HR:          ref to a hash with "esl-reformat" and "esl-alimerge"
#   $mdl_info_AHR:      ref to the model info hash of arrays, PRE-FILLED
#   $chunk_outdir_AR:   ref to array of output directories with files we are merging
#   $opt_HHR:           ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $ofile_info_HHR:    ref to the 2D hash of output file information, ADDED TO HERE 
#
# Returns:     void
# 
# Dies: if problem merging alignments
# 
################################################################# 
sub vdr_MergeAlignments { 
  my $nargs_exp = 6;
  my $sub_name = "vdr_MergeAlignments";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($out_root_no_vadr, $execs_HR, $mdl_info_AHR, $chunk_outdir_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # determine which type of alignment files we will merge:
  my $do_keep      = opt_Get("--keep", $opt_HHR);
  my $do_out_stk   = $do_keep || opt_Get("--out_stk",   $opt_HHR) ? 1 : 0;
  my $do_out_afa   = $do_keep || opt_Get("--out_afa",   $opt_HHR) ? 1 : 0;
  my $do_out_rpstk = $do_keep || opt_Get("--out_rpstk",   $opt_HHR) ? 1 : 0;
  my $do_out_rpafa = $do_keep || opt_Get("--out_rpafa",   $opt_HHR) ? 1 : 0;
  # NOTE: v-annotate.pl should have required that if --out_afa is used, --out_stk was also used
  # and that if --out_rpafa is used, --out_rpstk was also used. These are required because
  # we can't merge afa files (due to lack of RF annotation) so we need the stockholm equivalents
  # but we check here again to be safe.
  if(($do_out_afa) && (! $do_out_stk)) { 
    ofile_FAIL("ERROR in $sub_name, trying to merge afa files but don't have stk files", 1, $FH_HR);
  }
  if(($do_out_rpafa) && (! $do_out_rpstk)) { 
    ofile_FAIL("ERROR in $sub_name, trying to merge rpafa files but don't have rpstk files", 1, $FH_HR);
  }

  my $nmdl = scalar(@{$mdl_info_AHR});
  my @filelist_A = (); # array of alignment files to merge
  my $out_stk = undef;
  for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    @filelist_A = ();
    my $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
    my $afa_key = undef;
    my $do_stk  = 0;
    my $do_afa  = 0;
    my $desc    = "";
    foreach my $stk_key ("stk", "rpstk") { 
      if($stk_key eq "stk") { 
        $afa_key = "afa";
        $do_stk  = $do_out_stk;
        $do_afa  = $do_out_afa;
        $desc = "model " . $mdl_name . " full original sequence alignment";
      }
      else { # stk_key eq "rpstk"
        $afa_key = "rpafa";
        $do_stk  = $do_out_rpstk;
        $do_afa  = $do_out_rpafa;
        $desc = "model " . $mdl_name . " full replaced sequence alignment";
      }
      my $stk_sfx = "." . $mdl_name . ".align." . $stk_key;
      my $afa_sfx = "." . $mdl_name . ".align." . $afa_key;
      vdr_MergeOutputGetFileList($out_root_no_vadr, $stk_sfx, 0, \@filelist_A, $chunk_outdir_AR, $FH_HR);
      if(scalar(@filelist_A) > 0) { 
        my $stk_file = $out_root_no_vadr . ".vadr" . $stk_sfx;
        my $afa_file = $out_root_no_vadr . ".vadr" . $afa_sfx;
        my $stk_list_file = $stk_file . ".list";
        ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".align." . $stk_key . ".list", $stk_list_file, $do_keep, $do_keep, 
                                        "model " . $mdl_name . " alignment list ($stk_key)");
        utl_AToFile(\@filelist_A, $stk_list_file, 1, $FH_HR);
        sqf_EslAlimergeListRun($execs_HR->{"esl-alimerge"}, $stk_list_file, "", $stk_file, "stockholm", $opt_HHR, $FH_HR);

        # remember if we are outputting afa we are also outputting stk, see 'NOTE:' in comment above

        # need to use esl-alimerge and esl-alimanip to merge and add RF column numbering to the alignment
        # in the future, if esl-alimerge will keep RFCOL columns in the alignment, we can fall back to 
        # using only esl-alimerge
        my $merge_and_manip_cmd  = $execs_HR->{"esl-alimerge"} . " --list --outformat stockholm --informat stockholm --dna $stk_list_file | ";
        $merge_and_manip_cmd    .= $execs_HR->{"esl-alimanip"} . " --num-rf --outformat stockholm --informat stockholm --dna - > $stk_file";
        utl_RunCommand($merge_and_manip_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

        ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".align." . $stk_key, $stk_file, 1, 1, $desc . "(stk)");
        if($do_afa) { 
          sqf_EslReformatRun($execs_HR->{"esl-reformat"}, "", $stk_file, $afa_file, "stockholm", "afa", $opt_HHR, $FH_HR);
          ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".align." . $afa_key, $afa_file, 1, 1, $desc . "(afa)");
        }
        else {
          
        }
      }
    }
  }

  return;
}

#################################################################
# Subroutine: vdr_ReplacePseudoCoordsStringCreate()
# Incept:     EPN, Sat Dec 18 08:31:23 2021
#
# Purpose:    Format the "pseudo_coords" string for v-annotate.pl's
#             \%{$rpn_output_HHR} given the necessary information.
#
# Args:
#  $missing_seq_start:         position in model where missing model region starts
#  $missing_seq_stop:          position in model where missing model region stops
#  $missing_mdl_start:         position in model where missing model region starts
#  $missing_mdl_stop:          position in model where missing model region stops
#  $count_n:                   number of Ns in missing seq region
#  $nmatch:                    number of non-Ns in missing seq region that matched expected nt, can be undef
#  $nmismatch:                 number of non-Ns in missing seq region that did not match expected nt, can be undef
#  $flush_direction:           '5p', '3p' or '-'
#  $replaced_flag              '1' if region was replaced, '0' if not
#
# Returns: substring to append to $rpn_output_HHR->{$seq_name}{"coords"}
# Dies:    never
#
#################################################################
sub vdr_ReplacePseudoCoordsStringCreate {
  my $sub_name = "vdr_ReplacePseudoCoordsStringCreate";
  my $nargs_exp = 9;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($missing_seq_start, $missing_seq_stop, $missing_mdl_start, $missing_mdl_stop, 
      $count_n, $nmatch, $nmismatch, $flush_direction, $replaced_flag) = @_;
  
  my $missing_seq_len = abs($missing_seq_stop - $missing_seq_start) + 1;
  my $missing_mdl_len = abs($missing_mdl_stop - $missing_mdl_start) + 1;
  my $missing_diff    = $missing_seq_len - $missing_mdl_len;
  
  my $ret_str = "";
  $ret_str .= "[S:" . $missing_seq_start . ".." . $missing_seq_stop . ","; # S: seq coords of missing region
  $ret_str .= "M:" . $missing_mdl_start . ".." . $missing_mdl_stop . ","; # M: mdl coords of missing region
  # add '!' if missing_diff != 0;
  if($missing_diff != 0) { 
    $ret_str .= "D:" . $missing_diff . "!,";      # D: missing mdl len - missing seq len
  }
  else { 
    $ret_str .= "D:0,";
  }
  $ret_str .= "N:" . $count_n . "/" . $missing_seq_len . ",";             # N: num Ns in missing seq region / missing seq len
  if((defined $nmatch) && (defined $nmismatch)) { 
    my $missing_non_n   = $nmatch + $nmismatch;
    $ret_str .= "E:" . $nmatch  . "/" . $missing_non_n . ",";             # E: num non-Ns that match expected / num non Ns
  }
  else { 
    $ret_str .= "E:?/?,";
  }

  $ret_str .= "F:" . $flush_direction . ",";

  if($replaced_flag) { 
    $ret_str .= "R:Y];";                                                   # R: Y if replaced, N if not
  }
  else { 
    $ret_str .= "R:N];";
  }

  return $ret_str;
}

#################################################################
# Subroutine: vdr_ReplacePseudoCoordsStringParse()
# Incept:     EPN, Fri Jul  1 09:24:10 2022
#
# Purpose:    Given a .rpn "pseudo_coords" string from v-annotate.pl's
#             \%{$rpn_output_HHR} created by vdr_ReplacePseudoCoordsStringCreate(),
#             parse it and return its components. Caller should pass 'undef'
#             for unwanted components.
#
# Args:
#  $pseudo_coords:     the rpn_output_HHR->{}{"pseudo_coords"} string 
#  $scoords_sgm_AR:    REF to array to return VADR seq coords in, one segment per array element always + strand, FILLED here, can be undef
#  $mcoords_sgm_AR:    REF to array to return VADR mdl coords in, one segment per array element always + strand, FILLED here, can be undef
#  $diff_AR:           REF to array to return diff values in (e.g. "D:0", "D:-4!", D:2!"), FILLED here, can be undef
#  $ncount_AR:         REF to array to return N count values in (e.g. "20/20", "0/30"), FILLED here, can be undef
#  $ecount_AR:         REF to array to return 'expected' values in (e.g. "10/20", "0/0", "?/?"), FILLED here, can be undef
#  $flush_AR:          REF to array to return 'flush' values in (e.g. "-", "5'", "3'"), FILLED here, can be undef
#  $replaced_AR:       REF to array to return 'replaced' values in (e.g. "Y", "N"), FILLED here, can be undef
#  $FH_HR:             REF to hash of file handles
#
# Returns: number of elements in $pseudo_coords, this will be size of returned arrays
#          also fills any defined $*_AR passed in
#
#################################################################
sub vdr_ReplacePseudoCoordsStringParse { 
  my $sub_name = "vdr_ReplacePseudoCoordsStringParse";
  my $nargs_exp = 9;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($pseudo_coords, $scoords_sgm_AR, $mcoords_sgm_AR, $diff_AR, $ncount_AR, 
      $ecount_AR, $flush_AR, $replaced_AR, $FH_HR) = @_;

  # printf("in $sub_name, pseudo_coords: $pseudo_coords\n");
  
  my $ntok = 0;
  if($pseudo_coords eq "-") { 
    return 0;
  }

  my @ret_scoords_sgm_A = ();
  my @ret_mcoords_sgm_A = ();
  my @ret_diff_A        = ();
  my @ret_ncount_A      = ();
  my @ret_ecount_A      = ();
  my @ret_flush_A       = ();
  my @ret_replaced_A    = ();

  my @pseudo_tok_A = split(";", $pseudo_coords); 
  foreach my $pseudo_tok (@pseudo_tok_A) { 
    $ntok++;
    if($pseudo_tok =~ /^\[S\:(\d+\.\.\d+)\,M\:(\d+\.\.\d+),D\:(\-?\d+\!?),N\:(\d+\/\d+),E\:([^\,]+),F:([^\,]+),R:([YN])\]$/) { 
      my($scoords, $mcoords, $diff, $ncount, $ecount, $flush, $replaced) = ($1, $2, $3, $4, $5, $6, $7);
      $scoords .= ":+"; # always + strand
      $mcoords .= ":+"; # always + strand
      push(@ret_scoords_sgm_A, $scoords);
      push(@ret_mcoords_sgm_A, $mcoords);
      push(@ret_diff_A,        $diff);
      push(@ret_ncount_A,      $ncount);
      push(@ret_ecount_A,      $ecount);
      push(@ret_flush_A,       $flush);
      push(@ret_replaced_A,    $replaced);
    }
    else { 
      ofile_FAIL("ERROR in $sub_name, unable to parse pseudo_coords token $pseudo_tok", 1, $FH_HR);
    }
  }

  if(defined $scoords_sgm_AR)  { @{$scoords_sgm_AR} = @ret_scoords_sgm_A; }
  if(defined $mcoords_sgm_AR)  { @{$mcoords_sgm_AR} = @ret_mcoords_sgm_A; }
  if(defined $diff_AR)         { @{$diff_AR}        = @ret_diff_A; }
  if(defined $ncount_AR)       { @{$ncount_AR}      = @ret_ncount_A; }
  if(defined $ecount_AR)       { @{$ecount_AR}      = @ret_ecount_A; }
  if(defined $flush_AR)        { @{$flush_AR}       = @ret_flush_A; }
  if(defined $replaced_AR)     { @{$replaced_AR}    = @ret_replaced_A; }

  return $ntok;
}

#################################################################
# Subroutine: vdr_BackwardsCompatibilityExceptions()
# Incept:     EPN, Fri Sep 15 12:58:40 2023
#
# Purpose:    Update any 'exception' related keys and values in 
#             @{$ftr_info_AH} that were accepted in pre-v1.6 versions
#             of vadr. We will make two types of changes:
#
#             - change separator character for multiple coords strings 
#               any values of keys that end with "_exc" from ";" (old)
#               to "," (new)
#               example: "1..3:+;4..6:+" (old) to "1..3:+,4..6:+" (new)
#
#             - for "nmax{ins,del}_exc" exception keys, convert values from 
#               old format <posn>:<value> to new format: <coords_segment>:<value>
#               example: "333:40" (old) to "333..333:+:40" (new)
#               
#             - for "xmax{ins,del}_exc" exception keys, convert values from 
#               old format <posn>:<value> to new format: <coords_segment>:<value>
#               *and* convert from relative protein (amino acid) coords to 
#               global nucleotide coords example: 
# 
#               "1:40" (old) in protein with coords "100..129:+" to "102..102:+:40"
#               
#             Old {n,x}max{ins,del}_exc exception key values are single position, so
#             it could be possible to convert those to multiple position ranges, 
#             but we don't do that. We leave those as single position segments.
#       
#             This subroutine must be called prior to calling 
#               vdr_ModelInfoValidateExceptionKeys()
#               vdr_FeatureInfoValidateExceptionKeys()
# 
# modelinfo
# dupregin_exc - ; to ,
# indfstrn_exc - ; to ,
#
# featureinfo
# nmaxins_exc - convert to coords-value
# nmaxdel_exc - convert to coords-value
# xmaxins_exc - convert to coords-value
# xmaxdel_exc - convert to coords-value
# frameshift_exc - ; to ,
# 
# Args:
#  $mdl_info_HR:       ref to the model info hash (for one model)
#  $ftr_info_AHR:      ref to the feature info array of hashes
#  $alt_info_HHR:      ref to the alert info hash of hashes
#  $FH_HR:             ref to hash of file handles
#
# Returns: void
#
#################################################################
sub vdr_BackwardsCompatibilityExceptions {
  my $sub_name = "vdr_BackwardsCompatibilityExceptions";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_info_HR, $ftr_info_AHR, $alt_info_HHR, $FH_HR) = @_;
  
  if(! defined $mdl_info_HR->{"length"}) { 
    ofile_FAIL("ERROR in $sub_name, model length undefined", 1, $FH_HR);
  }
  my $mdl_len = $mdl_info_HR->{"length"};

  my $exc_key; 
  my $new_key = undef;
  foreach $exc_key (sort keys %{$mdl_info_HR}) { 
    # only 2 possibilities:
    # dupregin_exc
    # indfstrn_exc --> indfstr_exc
    # for both of these we need to swap ';' for ','
    if($exc_key =~ /^.+\_exc$/) { 
      # swap ; with ,
      $mdl_info_HR->{$exc_key} =~ s/\;/\,/g;
      $mdl_info_HR->{$exc_key} =~ s/\,$//; # remove final ',' if any
    }
    # and update indfstrn exc-key: indfstrn_exc --> indfstr_exc
    if($exc_key eq "indfstrn_exc") { 
      $new_key = $alt_info_HHR->{"indfstrn"}{"exc_key"};
      $mdl_info_HR->{$new_key} = $mdl_info_HR->{$exc_key};
      delete($mdl_info_HR->{$exc_key});
    }
  }

  # move onto ftr_info_AHR
  my $nftr = scalar(@{$ftr_info_AHR});
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    foreach my $exc_key (sort keys %{$ftr_info_AHR->[$ftr_idx]}) { 
      if($exc_key =~ /^.+\_exc$/) { 
        # swap ; with ,
        $ftr_info_AHR->[$ftr_idx]{$exc_key} =~ s/\;/\,/g;
        $ftr_info_AHR->[$ftr_idx]{$exc_key} =~ s/\,$//; # remove final ',' if any
      }
      if(($exc_key eq "nmaxins_exc") || ($exc_key eq "nmaxdel_exc") || ($exc_key eq "xmaxins_exc") || ($exc_key eq "xmaxdel_exc")) { 
        my $new_value = "";
        my @posn_value_A = split(",", $ftr_info_AHR->[$ftr_idx]{$exc_key});
        foreach my $posn_value (@posn_value_A) { 
          if($posn_value =~ /(\d+)\:(\d+)/) { 
            my ($posn, $value) = ($1, $2);
            if($new_value ne "") { $new_value .= ","; }
            my $nt_posn = $posn; # nt position
            if(($exc_key eq "xmaxins_exc") || ($exc_key eq "xmaxdel_exc")) { 
              # convert from protein to nucleotide coords
              $nt_posn = vdr_Feature3pMostPosition(vdr_CoordsProteinRelativeToAbsolute($ftr_info_AHR->[$ftr_idx]{"coords"}, 
                                                                                       vdr_CoordsSinglePositionSegmentCreate($posn, "+", $FH_HR),
                                                                                       $FH_HR), $FH_HR);
            }
            $new_value .= vdr_CoordsSinglePositionSegmentCreate($nt_posn, "+", $FH_HR) . ":" . $value;
          }
          else { 
            ofile_FAIL("ERROR, in $sub_name, trying to update old exception key $exc_key, but unable to parse value: $posn_value", 1, $FH_HR);
          }
        }
        if($new_value eq "") { 
          ofile_FAIL("ERROR, in $sub_name, trying to update old exception key $exc_key, but unable to parse value: " . $ftr_info_AHR->[$ftr_idx]{$exc_key}, 1, $FH_HR);
        }
        if($exc_key eq "nmaxins_exc") { 
          $new_key = $alt_info_HHR->{"insertnn"}{"exc_key"};
        }
        elsif($exc_key eq "nmaxdel_exc") { 
          $new_key = $alt_info_HHR->{"deletinn"}{"exc_key"};
        }
        if($exc_key eq "xmaxins_exc") { 
          $new_key = $alt_info_HHR->{"insertnp"}{"exc_key"};
        }
        if($exc_key eq "xmaxdel_exc") { 
          $new_key = $alt_info_HHR->{"deletinp"}{"exc_key"};
        }
        if(! defined $new_key) { 
          ofile_FAIL("ERROR, in $sub_name, trying to update old exception key $exc_key, but unable to determine new key", 1, $FH_HR);
        }
        $ftr_info_AHR->[$ftr_idx]{$new_key} = $new_value;
        delete($ftr_info_AHR->[$ftr_idx]{$exc_key});
      }
      if($exc_key eq "frameshift_exc") { 
        $new_key = $alt_info_HHR->{"fstukcfi"}{"exc_key"};
        $ftr_info_AHR->[$ftr_idx]{$new_key} = $ftr_info_AHR->[$ftr_idx]{$exc_key};
        delete($ftr_info_AHR->[$ftr_idx]{$exc_key});
      }
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
#  $esl_ssplit:       path to the esl-ssplit.pl script to use
#  $fasta_file:       fasta file to split up
#  $nfiles:           desired number of files to split $fasta_file into, -1 for one file for each sequence
#  $nseq_per_file_AR: [0..$nfiles-1]: number of seqs in each file, can be undef if unwanted
#  $opt_HHR:          REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:   REF to 2D hash of output file information
# 
# Returns:    Number of files actually created (can differ from requested
#             amount (which is $nfiles)).
#
# Dies:       if esl-ssplit command fails, or unable to open esl-ssplit output file
#
################################################################# 
sub vdr_SplitFastaFile { 
  my $sub_name = "vdr_SplitFastaFile()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_ssplit, $fasta_file, $nfiles, $nseq_per_file_AR, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to ofile_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $outfile = $fasta_file . ".esl-ssplit";
  my $cmd = undef;
  if($nfiles == -1) { # special case: put 1 sequence per file
    $cmd = "$esl_ssplit -v $fasta_file 1 > $outfile";
  }
  else { 
    $cmd = "$esl_ssplit -v -r -n $fasta_file $nfiles > $outfile";
  }
  utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # parse output to determine exactly how many files were created:
  # $esl_ssplit will have output exactly 1 line per fasta file it created
  my $nfiles_created = utl_FileCountLines($outfile, $FH_HR);

  if(defined $nseq_per_file_AR) { 
    open(IN, $outfile) || ofile_FileOpenFailure($outfile, $sub_name, $!, "reading", $FH_HR);
    while(my $line = <IN>) { 
      #va-r400/va-r400.vadr.in.fa.1 finished (11 seqs, 325770 residues)
      if($line =~ /^\S+\s+finished\s+\((\d+)\s+seqs/) { 
        push(@{$nseq_per_file_AR}, $1);
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, unable to parse esl-ssplit.pl output file $outfile line:\n$line\n", 1, $FH_HR);
      }
    }
    close(IN);
  }

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
# Subroutine: vdr_FrameAdjust
# Incept:     EPN, Fri May  1 17:46:54 2020
# Purpose:    Return the frame (1, 2, or 3) implied by a difference 
#             of <nt_diff> nt relative to <orig_frame>. I find this 
#             calculation unnervingly difficult to wrap my head around
#             for some reason.
#
# Arguments:
#  $orig_frame: 1, 2, or 3, you can also think of this as a codon_start
#               value:
#               '1': first position of sequence is first position to translate
#                  (first position of first codon)
#               '2': second position of sequence is first position to translate
#                  (first position of first codon, this implies that first 
#                  position of sequence is in frame *3* (third position of 
#                  'previous' codon)
#               '3': third position of sequence is first position to translate
#                  (first position of first codon, this implies that first 
#                  position of sequence is in frame *2* (second position of 
#                  'previous' codon)
#  $nt_diff:    number of nucleotides to adjust by; if positive, adjusting 
#               'forwards' (downstream, 3' direction); if negative, adjusting 
#               'backwards' (upstream, 5' direction)
#  $FH_HR:      ref to hash of file handles, including "cmd"
#
# Returns:  <ret_frame>: $orig_frame adjusted by $nt_diff nt.
#
#           <orig_frame> <nt_diff>  <ret_frame>
#           1           -2          3
#           1           -1          2
#           1            0          1
#           1            1          3
#           1            2          2
#
#           2           -2          1
#           2           -1          3
#           2            0          2
#           2            1          1
#           2            2          3
#
#           3           -1          1
#           3           -2          2
#           3            0          3
#           3            1          2
#           3            2          1
#
# Dies:     if $orig_frame is not 1, 2, or 3
#
#################################################################
sub vdr_FrameAdjust { 
  my $sub_name = "frame_adjust";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($orig_frame, $nt_diff, $FH_HR) = (@_);

  if(($orig_frame ne "1") && ($orig_frame ne "2") && ($orig_frame ne "3")) { 
    ofile_FAIL("ERROR in $sub_name, orig_frame must be 1, 2, or 3, got $orig_frame", 1, $FH_HR);
  }

  return (($orig_frame - $nt_diff - 1) % 3) + 1;
}

#################################################################
# Subroutine:  vdr_WriteCommandScript()
# Incept:      EPN, Fri Nov  9 14:26:07 2018 (ribo_WriteCommandScript)
#
# Purpose  : Create a new file to be executed as a job created by 
#            a qsub call.
# 
# Arguments: 
#   $file:  name of the file to create
#   $cmd:   the command to put in the file
#   $FH_HR: ref to hash of file handles, including "cmd"
#
# Returns:     void
# 
# Dies:        Never
#
################################################################# 
sub vdr_WriteCommandScript {
  my $nargs_exp = 3;
  my $sub_name = "vdr_WriteCommandScript";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($file, $cmd, $FH_HR) = @_;

  open(OUT, ">", $file) || ofile_FileOpenFailure($file, $sub_name, $!, "writing", $FH_HR);

  print OUT ("#!/bin/bash\n");
  print OUT ("#filename: $file\n");
  print OUT $cmd . "\n";

  close(OUT);

  return;
}

#################################################################
# Subroutine:  vdr_GlsearchFormat3And9CToStockholmAndInsertFile()
# Incept:      EPN, Wed Feb 17 09:21:08 2021
#
# Purpose  : Convert a FASTA package glsearch output file with 
#            format 3 (fasta) and 9C ("-m 3,9C") and >=1 query/target
#            alignments to a single Stockholm format multiple alignment
#            and create a corresponding insert file by parsing the cigar
#            strings. The glsearch output must have been run with:
#            "-m 3,9C": to specify format
#            "-z -1":   to turn off significance calculations
#            "-n":      to specify query is nucleotide
#            "-3":      to specify only top strand of query is searched
#            "-d 1":    to specify max number of alignments displayed is 1
# 
# Arguments: 
#   $gls_file:           name of output file from glsearch
#   $stk_file:           root name of stockholm files to create, we'll make one per seq
#   $insert_file:        name of insert file for all seqs to write
#   $blastn_db_sqfile_R: ref to open Bio:Easel:SqFile with model/target sequence
#   $exp_mdl_name:       expected single target sequence name
#   $opt_HHR:            ref to 2D hash of option values, see top of sqp_opts.pm for description
#   $ofile_info_HHR:     ref to 2D hash of output file information
#
# Returns:     Number of individual sequence $stk files created,
#              same as number of queries read.
# 
# Dies:        If there's a problem parsing the glsearch output.
#
################################################################# 
sub vdr_GlsearchFormat3And9CToStockholmAndInsertFile {
  my $nargs_exp = 7;
  my $sub_name = "vdr_GlsearchFormat3And9CToStockholmAndInsertFile";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($gls_file, $stk_file, $insert_file, $blastn_db_sqfile_R, $exp_mdl_name, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  #printf("in $sub_name\n\tgls_file: $gls_file\n\tstk_file: $stk_file\n\tinsert_file: $insert_file\nexp_mdl_name: $exp_mdl_name\n\n");
  open(IN,   $gls_file)  || ofile_FileOpenFailure($gls_file,  $sub_name, $!, "reading", $FH_HR);

  # fetch the model sequence, so we can use it to add to RF in alignments

  my $t_uaseq = $$blastn_db_sqfile_R->fetch_seq_to_sqstring($exp_mdl_name);

  my $q_name;         # name of query sequence
  my $q_len;          # length of query sequence
  my $nq = 0;         # number of queries read
  my $t_name;         # name of target sequence
  my ($an0, $ax0);    # start/stop position of alignment in query
  my ($an1, $ax1);    # start/stop position of alignment in library (target)
  my ($pn0, $px0);    # start/stop position of displayed query
  my ($pn1, $px1);    # start/stop position of displayed library (target)
  my $cigar;          # CIGAR string 
  my %q_len_H = ();   # key is query/sequence name, value is length
  my @q_name_A = ();  # array of query names
  # hash for storing insert info we will write to insert_file
  my %q_inserts_HH = (); # key 1: sequence name
                         # key 2: one of 'spos', 'epos', 'ins'
                         # $q_inserts_HHR->{}{"spos"} is starting model position of alignment
                         # $q_inserts_HHR->{}{"epos"} is ending model position of alignment
                         # $q_inserts_HHR->{}{"ins"} is the insert string in the format:
                         # <mdlpos_1>:<uapos_1>:<inslen_1>;...<mdlpos_n>:<uapos_n>:<inslen_n>;
                         # for n inserts, where insert x is defined by:
                         # <mdlpos_x> is model position after which insert occurs 0..mdl_len (0=before first pos)
                         # <uapos_x> is unaligned sequence position of the first aligned nt
                         # <inslen_x> is length of the insert

  # First 4 lines should look like this:
  ## /panfs/pan1/infernal/notebook/21_0213_vadr_hmmalign/fasta-experimenting-20210216/fasta-36.3.8h/bin/glsearch36 -z -1 -T 1 -3 -m 9C,3 -d 1 va-gls-cdc5/va-gls-cdc5.vadr.NC_045512.a.subseq.fa va-gls-cdc5/va-gls-cdc5.vadr.NC_045512.glsearch.fa
  #GLSEARCH performs a global-query/local-library search
  # version 36.3.8h May, 2020
  #Query: va-gls-cdc5/va-gls-cdc5.vadr.NC_045512.a.subseq.fa

  # in rare cases there's a line like the following before the Query line:
  # cannot read format 95 != lib_type 0

  # validate line 1
  ## /panfs/pan1/infernal/notebook/21_0213_vadr_hmmalign/fasta-experimenting-20210216/fasta-36.3.8h/bin/glsearch36 -z -1 -T 1 -3 -m 9C,3 -d 1 va-gls-cdc5/va-gls-cdc5.vadr.NC_045512.a.subseq.fa va-gls-cdc5/va-gls-cdc5.vadr.NC_045512.glsearch.fa
  my $line_ctr = 0;
  my $line = undef;
  $line = <IN>; $line_ctr++;
  chomp $line;
  if($line =~ m/^\#(.+)$/) { 
    my $first_line = $1;
    if(($first_line !~ m/\-m 3,9C/) && ($first_line !~ m/\-m 9C,3/)) { 
      ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, did not find \"-m 3,9C\" or \"-m 9C,3\" in first line:\n$line\n", 1, $FH_HR);
    }
  }
  else { 
    ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, first line is in unexpected format\n$line\n", 1, $FH_HR);
  }

  # validate line 2
  #GLSEARCH performs a global-query/local-library search
  $line = <IN>; $line_ctr++;
  chomp $line;
  if($line !~ m/^GLSEARCH/) { 
    ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, second line does not start with \"GLSEARCH\":\n$line\n", 1, $FH_HR);
  }

  # validate line 3
  # version 36.3.8h May, 2020
  $line = <IN>; $line_ctr++;
  chomp $line;
  if($line !~ m/version/) { 
    ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, third line did not include \"version\"\n", 1, $FH_HR);
  }

  # validate line 4
  #Query: va-gls-cdc5/va-gls-cdc5.vadr.NC_045512.a.subseq.fa
  $line = <IN>; $line_ctr++;
  chomp $line;
  while(($line !~ m/^Query/) && ($line =~ m/^\s+cannot read.*format/)) { 
    $line = <IN>; $line_ctr++;
    chomp $line;
  }
  if($line !~ m/^Query/) { 
    ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, fourth line did not start with \"Query\" (nor did it contain 'cannot read.*format')\n", 1, $FH_HR);
  }

  my $mdl_name;    # name of single target seq
  my $mdl_len;     # length of single target seq
  my $cur_mdl_len; # length of single target seq
  my $nseq;        # number of target seqs, should only be 1
  my $keep_going = 1;
  while((defined ($line = <IN>)) && ($keep_going)) { 
    $line_ctr++;
    chomp $line;
    #print("line: $line\n");
    if($line =~ /^\>\>\>\/\/\/$/) { 
      # end of all alignments
      $keep_going = 0;
    }
    elsif($line =~ /^\s*\d+\>\>\>(\S+).*(\d+)\s+nt/) { 
      # 1>>>lcl|SARS-CoV-2/human/USA/IN-CDC-LC00002770/2021/17579-27826 - 10248 nt (forward-only)
      #start of new query
      ($q_name, $q_len) = ($1, $2);
      push(@q_name_A, $q_name);
      $nq++;
      $q_len_H{$q_name} = $q_len;
      # parse next two lines
      $line = <IN>; $line_ctr++;
      if($line !~ m/^Library/) { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, first line after >>> line (line $line_ctr) does not start with Library\n$line\n", 1, $FH_HR);
      }
      $line = <IN>; $line_ctr++;
      if($line =~ /^\s*(\d+)\s+residues\s+in\s*(\d+)\s+sequences/) { 
        ($cur_mdl_len, $nseq) = ($1, $2);
        if(! defined $mdl_len) { 
          $mdl_len = $cur_mdl_len; 
        }
        elsif($cur_mdl_len != $mdl_len) { 
          ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, single target seq length $cur_mdl_len differs from previously read length $mdl_len", 1, $FH_HR);
        }
        if($nseq ne "1") { 
          ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, not exactly 1 sequence in target on line: $line_ctr\n$line\n", 1, $FH_HR);
        }
      }
      # validate next 5 lines:
      # <blank line>
      #Algorithm: Global/Local affine Needleman-Wunsch (SSE2, Michael Farrar 2010) (6.0 April 2007)
      #Parameters: +5/-4 matrix (5:-4), open/ext: -12/-4
      # <blank line>
      #The best scores are:                                                n-w	%_id  %_sim  gnw  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs  aln_code
      #NC_045512                                                (29903) [f] 21396	0.516 0.516 21396 10248    1 10248    1 10248 17633 27872    1 29903 4952   8   0	2949M1D5129M7D2162M
      $line = <IN>; $line_ctr++; # blank line
      $line = <IN>; $line_ctr++;
      if($line !~ /^Algorithm/) { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, does not begin with Algorithm", 1, $FH_HR);
      }
      $line = <IN>; $line_ctr++;
      if($line !~ /^Parameters/) { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, does not begin with Parameters", 1, $FH_HR);
      }
      $line = <IN>; $line_ctr++; # blank line
      $line = <IN>; $line_ctr++;
      if($line !~ /^The best scores/) { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, does not begin with The best scores", 1, $FH_HR);
      }
      $line = <IN>; $line_ctr++;
      # line with the info on the alignment we need to parse
      # Two examples below ("The best scores" line kept only to show field names
      # note the space in the parentheses before the model name, we have to take special care to deal with that
      #The best scores are:                                                n-w	%_id  %_sim  gnw  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs  aln_code
      #NC_045512                                                (29903) [f] 21396	0.516 0.516 21396 10248    1 10248    1 10248 17633 27872    1 29903 4952   8   0	2949M1D5129M7D2162M
      #entoy100a-dcr                                            ( 100) [f]  500	1.000 1.000  500  100    1  100    1  100    1  100    1  100   0   0   0	100M
      if($line =~ /^(\S+)\s*\(\s*\d+\)\s+\S+\s+/) { 
        # matches "entoy100a-dcr                                            ( 100) [f]"
        ($mdl_name) = $1;
        my $subline = $line;
        $subline =~ s/^\S+\s*\(\s*\d+\)\s+\S+\s+//;
        my @el_A = split(/\s+/, $subline);
        if(scalar(@el_A) != 17) { 
          ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, line $line_ctr, unable to parse best scores (1) line:\n$line\n", 1, $FH_HR);
        }
        ($an0, $ax0, $pn0, $px0, $an1, $ax1, $pn1, $px1, $cigar) = 
            ($el_A[5], $el_A[6], $el_A[7], $el_A[8], $el_A[9], $el_A[10], $el_A[11], $el_A[12], $el_A[16]);
        if($mdl_name ne $exp_mdl_name) { 
          ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, line $line_ctr, expected single target sequence name $mdl_name but read $mdl_name", 1, $FH_HR);
        }
        # parse cigar to get inserts in query to later write to insert_file
        vdr_CigarToInsertsHash(\%{$q_inserts_HH{$q_name}}, $cigar, $an0, $an1, $FH_HR);
      }
      else { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, unable to parse best scores (2) line:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>; $line_ctr++; # blank line
      $line = <IN>; $line_ctr++;
      if($line =~ /^\>\>\>(\S+)\,\s*/) { 
          if($1 ne $q_name) { 
          ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, expected >>>$q_name line preceding alignment but got:\n$line\n", 1, $FH_HR);
        }
      }
      else { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, expected >>>$q_name line preceding alignment but got:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>; $line_ctr++; # blank line
      $line = <IN>; $line_ctr++;
      if($line !~ /^global\/local score/) {   
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, expected line beginning with \"global/local\", but got:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>; $line_ctr++;
      my $q_afa = "";
      my $nspace_5p = 0;
      my $nspace_3p = 0;
      if($line !~ /^\>/) { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, expected line beginning with \">\" indicating beginning of q seq, but got:\n$line\n", 1, $FH_HR);
      }
      $line = <IN>; $line_ctr++;
      while($line !~ m/^\>/) { 
        if(! defined $line) { 
          ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, ran out of lines before target aligned seq", 1, $FH_HR);
        }
        chomp $line; 
        # do not remove leading/trailing whitespace, we deal with this after
        # we've read the full seq
        $q_afa .= $line;
        
        $line = <IN>; $line_ctr++;
      }
      # currently line is ">" indicating start of target alignment
      if($line =~ /^\>\>\>\<\<\</) { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, read end of alignment before aligned target\n", 1, $FH_HR);
      }
      $line = <IN>; $line_ctr++;
      my $t_afa = "";
      while($line !~ /^\>\>\>\<\<\</) { 
        if(! defined $line) { 
          ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, ran out of lines before end of target aligned seq", 1, $FH_HR);
        }
        chomp $line; 
        $t_afa .= $line;
        $line = <IN>; $line_ctr++;
      }
      # currently line is ">>><<<" indicating end of alignment for this query seq and target seq
      # make sure only spaces are at beginning and end of query seq,
      # count how many, and remove them and corresponding number of nt
      # from beginning/end of target seq too.
      if($q_afa =~ m/^(\s*)\S+(\s*)$/) { 
        $nspace_5p = length($1);
        $nspace_3p = length($2);
        $q_afa =~ s/^\s+//; # remove leading  whitespace
        $q_afa =~ s/\s+$//; # remove trailing whitespace
      }
      else { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, did not read query aligned sequence correctly:\n$q_afa\n", 1, $FH_HR);
      }
      if($t_afa !~ m/^\S+$/) { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, on line $line_ctr, did not read target aligned sequence correctly:\n$q_afa\n", 1, $FH_HR);
      }
      if($nspace_5p > 0) { 
        $t_afa = substr($t_afa, $nspace_5p);
      }
      if($nspace_3p > 0) { 
        $t_afa = substr($t_afa, 0, -1 * $nspace_3p);
      }
      my $q_len = length($q_afa);
      my $t_len = length($t_afa);
      
      if($q_len != $t_len) { 
        ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, at line $line_ctr; aligned query length $q_len differs from aligned target length $t_len", 1, $FH_HR);
      }

      # in target, which will become RF, replace - characters with '.' following hmmer/infernal convention
      $t_afa =~ s/\-/\./g; 

      # add 5' and 3' ends of target, if nec
      my $t_5p = "";
      my $q_5p = "";
      if($an1 > 1) { 
        $t_5p = substr($t_uaseq, 0, ($an1-1));
        $q_5p = utl_StringMonoChar(($an1-1), "-", undef); 
      }
      my $t_3p = "";
      my $q_3p = "";
      if($ax1 < $mdl_len) { 
        $t_3p = substr($t_uaseq, ($ax1 - $mdl_len));
        $q_3p = utl_StringMonoChar(($mdl_len - $ax1), "-", undef); 
      }
      $q_afa = $q_5p . $q_afa . $q_3p;
      $t_afa = $t_5p . $t_afa . $t_3p;

      my $q_name_len = length($q_name);

      my $cur_stk_file = $stk_file . "." . $nq;
      open(OUT, ">", $cur_stk_file) || ofile_FileOpenFailure($cur_stk_file, $sub_name, $!, "writing", $FH_HR);
      printf OUT ("# STOCKHOLM 1.0\n\n");
      printf OUT ("%-*s  %s\n", $q_name_len, $q_name, $q_afa);
      printf OUT ("%-*s  %s\n", $q_name_len, "#=GC RF",   $t_afa);
      print  OUT ("//\n");
      close(OUT);
    }
    else { 
      ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, at line $line_ctr, expected line beginning with \\d+>>> indicating next query or >>>/// line indicating end of alignments, got:\n$line\n", 1, $FH_HR);
    }
  }
  if($nq == 0) { 
    ofile_FAIL("ERROR, in $sub_name, parsing $gls_file, did not read any alignments\n", 1, $FH_HR);
  }

  # write insert file
  vdr_CmalignWriteInsertFile($insert_file, 0, $exp_mdl_name, $mdl_len, \@q_name_A, \%q_len_H, \%q_inserts_HH, $FH_HR);

  return $nq;
}

#################################################################
# Subroutine:  vdr_CigarToInsertString()
# Incept:      EPN, Wed Feb 17 18:50:56 2021
#
# Purpose:    Given a CIGAR string where one sequence in the alignment
#             is a model sequence, determine insert information and
#             add it to %{$inserts_HR}, where keys are:
#               "spos" is starting model position of aligned sequence
#               "epos" is ending model position of aligned sequence
#               "ins"  is the insert string in the format:
#                      <mdlpos_1>:<uapos_1>:<inslen_1>;...<mdlpos_n>:<uapos_n>:<inslen_n>;
#                      for n inserts, where insert x is defined by:
#                      <mdlpos_x> is model position after which insert occurs 0..mdl_len (0=before first pos)
#                      <uapos_x> is unaligned sequence position of the first aligned nt
#                      <inslen_x> is length of the insert
#             CIGAR is in formation (\d+[MID])+
#             where \d+ indicates length
#             M indicates matches (no inserts or deletes)
#             I indicates insertion in target/model, so deletion  in query
#             D indicates deletion  in target/model, so insertion in query (and so stored in %{$inserts_HR}{"ins"})
#
# Reference:  https://en.wikipedia.org/wiki/Sequence_alignment#Representations
#             https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/
# 
# Arguments: 
#   $inserts_HR:  ref to hash to fill, see 'Purpose' for keys
#   $cigar:       CIGAR string
#   $seqstart:    first sequence position of alignment (typically 1)
#   $mdlstart:    first model RF position of alignment (varies)
#   $FH_HR:       ref to hash of file handles, including "cmd"
#
# Returns:     void
# 
# Dies:        If unable to parse $cigar string
#
################################################################# 
sub vdr_CigarToInsertsHash { 
  my $nargs_exp = 5;
  my $sub_name = "vdr_CigarToInsertHash";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($inserts_HR, $cigar, $seqstart, $mdlstart, $FH_HR) = @_;
  
  # printf("in $sub_name, cigar: $cigar, seqstart: $seqstart mdlstart: $mdlstart\n");

  my $seqpos = $seqstart;
  my $mdlpos = $mdlstart;
  my $orig_cigar = $cigar;
  my $spos = undef;
  my $epos = undef;
  my $ins_str = "";
  while($cigar ne "") { 
    if($cigar =~ /^(\d+)([MID])/) {
      my ($len, $type) = ($1, $2);
      if($type eq "M") { 
        $seqpos += $len;
        if(! defined $spos) { $spos = $mdlpos; }
        $mdlpos += $len;
        $epos = $mdlpos - 1;
      }
      if($type eq "I") { 
        $mdlpos += $len;
      }
      if($type eq "D") { 
        $ins_str .= ($mdlpos-1) . ":" . $seqpos . ":" . $len . ";";
        $seqpos += $len;
      }
      $cigar =~ s/^\d+[MID]//;
    }
    else { 
      ofile_FAIL("ERROR, in $sub_name, unable to parse cigar string $orig_cigar", 1, $FH_HR);
    }
  }
  if(! defined $spos) { 
    ofile_FAIL("ERROR, in $sub_name, unable to determine spos for mdlstart: $mdlstart and cigar: $cigar", 1, $FH_HR);
  }
  if(! defined $epos) { 
    ofile_FAIL("ERROR, in $sub_name, unable to determine spos for mdlstart: $mdlstart and cigar: $cigar", 1, $FH_HR);
  }
  # printf("returning spos: $spos epos: $epos ins: $ins_str\n");
  $inserts_HR->{"spos"} = $spos; 
  $inserts_HR->{"epos"} = $epos; 
  $inserts_HR->{"ins"} = $ins_str;

  return;
}

#################################################################
# Subroutine:  vdr_CigarToPositionMap()
# Incept:      EPN, Mon Jul 26 12:36:38 2021
#
# Purpose:    Given a CIGAR string, model length <mlen> and sequence 
#             length <slen>, fill an array <map_AR> (1..<mlen>)
#             indicating which sequence position each model position
#             maps to.
#             map_AR->[mpos] = spos    indicates that sequence position 
#                                      spos is aligned to model position mpos
#             map_AR->[mpos] = -spos   indicates that model position mpos
#                                      is aligned to a gap in the sequence
#                                      and 5' most sequence position seen so
#                                      far (5'-most prior to gap) is spos
#             map_AR->[mpos] = 0       indicates that model position mpos
#                                      is aligned to a gap in the sequence
#                                      and no positions have been seen so far
#                                      in the sequence 
#             CIGAR string $cigar is in format (\d+[MID])+
#             where \d+ indicates length
#             M indicates matches (no inserts or deletes)
#             I indicates insertion in target/model, so deletion(gap) in query/sequence
#             D indicates deletion  in target/model, so insertion     in query/sequence (gap in target/model)
#
# Reference:  https://en.wikipedia.org/wiki/Sequence_alignment#Representations
#             https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/
# 
# Arguments: 
#   $map_AR:      ref to array to fill with map, see 'Purpose' for more info
#   $cigar:       CIGAR string
#   $mdllen:      length of model
#   $seqlen:      length of sequence
#   $FH_HR:       ref to hash of file handles, including "cmd", can be undef
#
# Returns:     void
# 
# Dies:        If unable to parse $cigar string
#              If cigar string implies model length different than $mdllen
#              If cigar string implies sequence length different than $seqlen
#
################################################################# 
sub vdr_CigarToPositionMap { 
  my $nargs_exp = 5;
  my $sub_name = "vdr_CigarToPositionMap";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($map_AR, $cigar, $mdllen, $seqlen, $FH_HR) = @_;
  
  #printf("in $sub_name, cigar: $cigar, mdllen: $mdllen, seqlen: $seqlen\n");

  my $seqpos = 0;
  my $mdlpos = 0;
  my $orig_cigar = $cigar;
  my $i;
  $map_AR->[0] = -2; # irrelevant
  while($cigar ne "") { 
    if($cigar =~ /^(\d+)([MID])/) {
      my ($len, $type) = ($1, $2);
      if($type eq "M") { 
        for($i = 1; $i <= $len; $i++) { 
          $map_AR->[($mdlpos+$i)] = $seqpos+$i;
        }
        $seqpos += $len;
        $mdlpos += $len;
      }
      if($type eq "I") { 
        # gap in sequence
        for($i = 1; $i <= $len; $i++) { 
          if($seqpos == 0) { 
            $map_AR->[($mdlpos+$i)] = 0;
          }
          else { 
            $map_AR->[($mdlpos+$i)] = -1 * $seqpos;
          }
        }
        $mdlpos += $len;
      }
      if($type eq "D") { 
        # don't need to update @{$map_AR}
        $seqpos += $len;
      }
      $cigar =~ s/^\d+[MID]//;
    }
    else { 
      ofile_FAIL("ERROR, in $sub_name, unable to parse cigar string $orig_cigar", 1, $FH_HR);
    }
  }
  if($mdlpos != $mdllen) { 
    ofile_FAIL("ERROR, in $sub_name, model length after mapping ($mdlpos) not equal to input model length ($mdllen)", 1, $FH_HR);
  }
  if($seqpos != $seqlen) { 
    ofile_FAIL("ERROR, in $sub_name, sequence length after mapping ($seqpos) not equal to input sequence length ($seqlen)", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine:  vdr_UpdateInsertTokenInInsertString()
# Incept:      EPN, Mon Mar 15 12:22:21 2021
#
# Purpose:    Given $ins_str an 'insert string' in the format:
#              <mdlpos_1>:<uapos_1>:<inslen_1>;...<mdlpos_n>:<uapos_n>:<inslen_n>;
#              for n inserts, where insert x is defined by:
#              <mdlpos_x> is model position after which insert occurs 0..mdl_len (0=before first pos)
#              <uapos_x> is unaligned sequence position of the first aligned nt
#              <inslen_x> is length of the insert
#             Find the 'insert token' <mdlpos_n>:<uapos_n>:<inslen_n> equal to $orig_ins_tok 
#             and update it with $new_ins_tok.
#
# Arguments: 
#   $ins_str:      insert string
#   $orig_ins_tok: insert token that should exist in $ins_str to replace
#   $new_ins_tok:  insert token to replace $orig_ins_tok with
#   $FH_HR:        ref to hash of file handles, including "cmd"
#
# Returns:     new insert string with $orig_ins_tok replaced with $new_ins_tok
# 
# Dies:        If we can't parse $ins_str
#              If $orig_ins_tok does not exist in $ins_str
#              If $orig_ins_tok exists more than once in $ins_str
#
################################################################# 
sub vdr_UpdateInsertTokenInInsertString { 
  my $nargs_exp = 4;
  my $sub_name = "vdr_UpdateInsertTokenInInsertString";
  if(scalar(@_) != $nargs_exp) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_exp); exit(1); } 

  my ($ins_str, $orig_ins_tok, $new_ins_tok, $FH_HR) = @_;

  # contract checks
  # $ins_str should be defined and not empty
  if(! defined $ins_str) { ofile_FAIL("ERROR in $sub_name, insert string is undefined", 1, $FH_HR); }
  if($ins_str eq "")     { ofile_FAIL("ERROR in $sub_name, insert string is empty", 1, $FH_HR); }

  # $orig_ins_tok should defined and valid
  if(! defined $orig_ins_tok) { ofile_FAIL("ERROR in $sub_name, orig_ins_tok string is undefined", 1, $FH_HR); }
  if($orig_ins_tok !~ /^(\d+)\:(\d+)\:(\d+)/) { ofile_FAIL("ERROR in $sub_name, unable to parse orig_ins_tok $orig_ins_tok", 1, $FH_HR); }
  
  # $new_ins_tok  should defined and valid
  if(! defined $new_ins_tok) { ofile_FAIL("ERROR in $sub_name, new_ins_tok string is undefined", 1, $FH_HR); }
  if($new_ins_tok !~ /^(\d+)\:(\d+)\:(\d+)/) { ofile_FAIL("ERROR in $sub_name, unable to parse new_ins_tok $new_ins_tok", 1, $FH_HR); }
  # end contract checks
  
  my @ins_A = split(";", $ins_str);
  my $found_orig_ins_tok = 0;
  my $ret_ins_str = "";
  foreach my $ins_tok (@ins_A) {
    if($ins_tok =~ /^(\d+)\:(\d+)\:(\d+)/) { 
      if($ins_tok eq $orig_ins_tok) { 
        if($found_orig_ins_tok) { 
          ofile_FAIL("ERROR in $sub_name, found original token $orig_ins_tok twice in insert string $ins_str", 1, $FH_HR);
        }
        $ret_ins_str .= $new_ins_tok . ";";
        $found_orig_ins_tok = 1;
      }
      else { 
        $ret_ins_str .= $ins_tok . ";";
      }
    }
    else {
      ofile_FAIL("ERROR in $sub_name, unable to parse insert string $ins_str at at token $ins_tok", 1, $FH_HR);
    }
  }
  if(! $found_orig_ins_tok) { 
    ofile_FAIL("ERROR in $sub_name, unable to find orig_ins_tok $orig_ins_tok in insert string $ins_str", 1, $FH_HR);
  }

  return $ret_ins_str;
}

###########################################################################
# the next line is critical, a perl module must return a true value
return 1;
###########################################################################
