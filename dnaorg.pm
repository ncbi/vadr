#!/usr/bin/perl
# 
# version: 0.37 [Nov 2018]
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
# - $FH_HR: a reference to a hash of file handles. Important keys are "log" and "cmd",
#           the log and command files we are outputting to. This data structure is passed
#           into nearly all functions because it is also passed into DNAORG_FAIL() which 
#           can be called from nearly all functions. DNAORG_FAIL() outputs an error message
#           both the summary and log files before exiting, and appends a # DNAORG-FAILURE
#           to those files before closing them and exiting. This is done so the user
#           has a record of the reason the execution of the program failed, and it is also
#           useful for debugging purposes.
#
# - $ofile_info_HHR: reference to an 'info hash of hashes' in which all 2nd dim hashes are the same
#                    size and have the same set of keys, and contain information on an output file
#                    created by the main script (e.g. dnaorg_build.pl). The 1st dim hash keys describe
#                    the type of information, e.g. "fullpath" is the full path of the file, and 2nd
#                    dim keys pertain to which file, e.g. "log" for the log file which contains 
#                    all the output printed to stdout during the course of the execution of the script.
#                    A special 1d key is 'order' which is used for keeping track of the order that
#                    the files are added in, mainly so we can output information on them in the same
#                    order. The HH values for "order" are 1..$num_ofiles, where $num_ofiles is the
#                    number of total output files (scalar(keys %{$ofile_info_HHR{"order"}})).
#                    See validateOutputFileInfoHashOfHashes() for a list and explanation of all
#                    of the keys.
# 
# - $ftr_info_HAR: reference to an 'info hash of arrays' in which all arrays are the same 
#                  size, and contain information on a 'feature' (e.g. CDS or mature peptide)
#                  The hash has specific keys (e.g. "ref_len"), each of which points to 
#                  an array with the relevant information. See validateFeatureInfoHashIsComplete()
#                  for a list and explanation of the keys. The contents of this data structure
#                  can be printed to file <f> by dnaorg_build.pl or dnaorg_annotate.pl with the
#                  --ftrinfo <f> command line option.
#                   
# - $mdl_info_HAR: similar to ${ftr,seq,err}_info_HAR, except contains information pertaining 
#                  to each model, >= 1 of which will model a single feature (1 model for single
#                  exon CDS, 2 models for dual exon CDS, etc.). See 
#                  validateModelInfoHashIsComplete() for a list and explanation of the keys.
#                  The contents of this data structure can be printed to file <f> by dnaorg_build.pl 
#                  or dnaorg_annotate.pl with the --mdlinfo <f> command line option.
#
# - $seq_info_HAR: similar to ${ftr,mdl,err}_info_HAR, except contains information pertaining 
#                  to each sequence. See validateSequenceInfoHashIsComplete()
#                  for a list and explanation of the keys. The contents of this data structure
#                  can be printed to file <f> by dnaorg_annotate.pl with the --seqinfo <f>
#                  command line option.
#                   
# - $err_info_HAR: similar to ${ftr,mdl}_info_HAR, except contains information pertaining to each 
#                  error code. See validateErrorInfoHashIsComplete() for a list and explanation 
#                  of the keys. The contents of this data structure can be printed to file <f>
#                  by dnaorg_annotate.pl with the --errinfo <f> command line option.
#                   
########################################################################################
#
# List of subroutines in this file, divided into categories. 
#
# The names of the subroutines are meant to be descriptive. The
# convention for subroutine names are to use camel-caps
# (e.g. determineFeatureTypes()), as opposed to local subroutine names
# in the scripts dnaorg_build.pl and dnaorg_annotate.pl which use
# underscores, e.g. determine_feature_types().
#
# Subroutines related to the feature and model info hash data
# structures:
#   getPrimaryOrAllChildrenFromFeatureInfo()
#   determineFeatureTypes()
#   getNumFeaturesAnnotatedByModels()
#   getReferenceFeatureInfo()
#   fetchReferenceFeatureSequences()
#   annotateAppendFeatures()
#   annotateOverlapsAndAdjacencies()
#
# Subroutines related to the output info hash:
#   openAndAddFileToOutputInfo()
#   addClosedFileToOutputInfo()
#   helperAddFileToOutputInfo()
#
# Subroutines related to the error info hash:
#   initializeHardCodedErrorInfoHash()
#   addToErrorInfoHash()
#   setIncompatibilityErrorInfoHash()
#   setRequiredErrorInfoHash()
#
# Subroutines related to the feature table error exception array of hashes:
#   initializeHardCodedFTableErrorExceptions()
#   addFTableErrorException()
#
# Subroutines related to overlaps and adjacencies:
#   overlapsAndAdjacenciesHelper()
#   getOverlap()
#   compareTwoOverlapOrAdjacencyIndexStrings()
#   checkForIndexInOverlapOrAdjacencyIndexString()
#
# Massive wrapper subroutines that call other subroutines:
#   wrapperGetInfoUsingEdirect()
#   wrapperFetchAllSequencesAndProcessReferenceSequence()
#
# Subroutines related to feature tables output from edirect:
#   edirectFtableOrMatPept2SingleFeatureTableInfo()
#   getSingleFeatureTableInfo()
#   helperBreakdownFac()
#
# Subroutines for parsing different file types:
#   parseMatPeptSpecFile()
#   parseLengthFile()
#   parseEdirectFtableFile()
#   parseEdirectMatPeptideFile()
#   parseListFile()
#   parseSpecStartFile()
#   parseConsOptsFile()
#   parseNonConsOptsFile()
#
# Subroutines related to parsing NCBI coordinate strings:
#   getStrandStats()
#   startsStopsStrandsFromCoordsLength()
#   startsStopsFromCoords()
#   getLengthsAndCoords()
#   lengthFromCoords()
#
# Subroutines related to parsing dash-coords (non-NCBI) coordinate strings:
#   dashCoordsStringCommaDelimitedToLength()
#   dashCoordsToLength()
#
# Subroutines related to output:
#   outputProgressPrior()
#   outputProgressComplete()
#   outputConclusionAndCloseFiles()
#   outputTiming()
#   outputString()
#   outputBanner()
#   outputDividingLine()
# 
# Subroutines for dumping data structures, usually for debugging:
#   dumpInfoHashOfArrays()
#   dumpArrayOfHashesOfHashes()
#   dumpArrayOfHashes()
#
# Subroutines for validating the special data structures:
#   validateExecutableHash()
#   validateFeatureInfoHashIsComplete()
#   validateModelInfoHashIsComplete()
#   validateSequenceInfoHashIsComplete()
#   validateErrorInfoHashIsComplete()
#   validateInfoHashOfArraysIsComplete()
#   validateOutputFileInfoHashOfHashes()
#   validateAndGetSizeOfInfoHashOfArrays()
#   getConsistentSizeOfInfoHashOfArrays()
#   validateFTableErrorExceptions()
#
# Subroutines related to codons:
#   fetchStopCodon()
#   fetchStartCodon()
#   fetchCodon()
#   validateStopCodon()
#
# Subroutines related to timings:
#   secondsSinceEpoch()
#   formatTimeString()
#
# Simple utility subroutines for hashes and arrays:
#   findNonNumericValueInArray()
#   numNonNumericValueInArray()
#   maxLengthScalarValueInHash()
#   maxLengthScalarValueInArray()
#   findValueInArray()
#
# Simple utility subroutines:
#   DNAORG_FAIL()
#   fileOpenFailure()
#   runCommand()
#   removeDirPath()
#   removeScriptNameFromString()
#   removeFileUsingSystemRm()
#   getMonocharacterString()
#   countLinesInFile()
#   fileLinesToArray()
#   arrayToNewlineDelimitedString()
#   hashKeysToNewlineDelimitedString()
#   hashValuesToNewlineDelimitedString()
#   validateFileExistsAndIsNonEmpty()
#   concatenateListOfFiles()
#   md5ChecksumOfFile()
#   nseBreakdown()
#
# Miscellaneous subroutines that don't fall into one of the above
# categories:
#   stripVersion()
#   fetchedNameToListName()
#   fetchSequencesUsingEslFetchCds()
#   addNameAndBlankSsToStockholmAlignment()
#   getQualifierValues()
#   createCmDb()
#   matpeptValidateCdsRelationships()
#   checkForSpanningSequenceSegments()
#   getIndexHashForArray()
#   waitForFarmJobsToFinish()
#   splitFastaFile()
#
use strict;
use warnings;
use Cwd;

#################################################################
#################################################################
#
# Subroutines related to the feature and model info hash data
# structures:
#   getPrimaryOrAllChildrenFromFeatureInfo()
#   determineFeatureTypes()
#   getNumFeaturesAnnotatedByModels()
#   getReferenceFeatureInfo()
#   fetchReferenceFeatureSequences()
#   annotateAppendFeatures()
#   annotateOverlapsAndAdjacencies()
#
#################################################################
# Subroutine : getPrimaryOrAllChildrenFromFeatureInfo()
# Incept:      EPN, Thu Mar  3 12:27:38 2016
#
# Purpose:     Fill @{$AR} with the space delimited tokens of 
#              @{$ftr_info_HAR->{"primary_children_ftr_str"}[$ftr_idx]}
#              (if $primary_or_all eq "primary") or of
#              @{$ftr_info_HAR->{"all_children_ftr_str"}[$ftr_idx]}
#              (if $primary_or_all eq "all") and return.
#              
# Arguments: 
#   $ftr_info_HAR:   REF to hash of arrays with information on the features, PRE-FILLED
#   $ftr_idx:        index of feature we're interested in
#   $primary_or_all: "primary" to fill array with primary children ("primary_children_ftr_str")
#                    "all" to fill array with all children ("all_children_ftr_str")
#   $AR:             REF to array to fill, FILLED HERE
#   $FH_HR:          REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: if feature info hash is not 'complete', or if
#       there are no primary/all children for feature index $ftr_idx,
#       or if it is not of type 'multifeature', or
#       if $primary_or_all is neither "primary" nor "all"
#
################################################################# 
sub getPrimaryOrAllChildrenFromFeatureInfo { 
  my $nargs_expected = 5;
  my $sub_name = "getPrimaryOrAllChildrenFromFeatureInfo";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_HAR, $ftr_idx, $primary_or_all, $AR, $FH_HR) = @_;

  if(($primary_or_all ne "primary") && ($primary_or_all ne "all")) { 
    DNAORG_FAIL("ERROR in $sub_name, special primary_or_all variable is neither \"primary\" nor \"all\"", 1, $FH_HR);
  }    

  if($ftr_info_HAR->{"annot_type"}[$ftr_idx] ne "multifeature") { 
    DNAORG_FAIL("ERROR in $sub_name, requesting info for feature index $ftr_idx, but it is not annotated type of multifeature.", 1, $FH_HR);
  }

  @{$AR} = ();
  if($primary_or_all eq "primary") { 
    @{$AR} = split(/\s+/, $ftr_info_HAR->{"primary_children_ftr_str"}[$ftr_idx]);
  }
  elsif($primary_or_all eq "all") { 
    @{$AR} = split(/\s+/, $ftr_info_HAR->{"all_children_ftr_str"}[$ftr_idx]);
  }

  if(scalar(@{$AR}) == 0) { 
    DNAORG_FAIL("ERROR in $sub_name, requesting info for feature index $ftr_idx, but it has no $primary_or_all children.", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: determineFeatureTypes()
# Incept:     EPN, Thu Feb 11 14:50:53 2016
#
# Purpose:    Determine the type of each feature, either
#             'cds-notmp': CDS not comprised of mature peptides
#             'cds-mp':    CDS comprised of mature peptides
#             'mp':        mature peptide
#             'xfeat'      extra feature annotated by its own model (e.g. gene)
#             'dfeat'      extra feature annotated as a duplicate of another feature (e.g. gene)
# 
#              And populate the %{$ftr_info_HAR} with the following
#              information (1D keys):
#              "type"        values are one of: "cds-notmp", "cds-mp", "mp", 'xfeat', or 'dfeat'
#              "type_idx"    values are indices of each feature for its type
#                            so a '2' means it's the 2nd of its type.
#              "annot_type": annotation type, either "model", "multifeature", or "duplicate"
#                            "model": a homology model's prediction annotates 
#                            these features
#                            "multifeature": derived from multiple "model" annot_type
#                            features, these are usually type: cds-mp, a CDS comprised
#                            of mature peptides.
#                            "duplicate": the annotation of this feature is copied from 
#                            another. Example is a 'gene' that is a copy of a 'CDS'.
#                            The feature that is the source to copy is stored in the
#                            "source_idx" key.
#               "primary_children_ftr_str": for features with annot_type eq "multifeature",
#                            string of feature indices that are the primary children of this
#                            feature (when the primary children are concatenated they make up
#                            the CDS that encodes all the primary children). Each index is
#                            separated by a space. Empty for features with annot_type eq 
#                            "model".
#               "primary_children_ftr_num": number of primary children (0 for features with
#                            annot_type eq "model").
#               "all_children_ftr_str": for features with annot_type eq "multifeature",
#                            string of all feature indices that are encoded by this CDS, 
#                            includes all primary children plus any secondarily processed
#                            mature peptides. Each index is separated by a space. Empty 
#                            for features with annot_type eq "model."
#               "all_children_ftr_num": number of all children (0 for features with
#                            annot_type eq "model").
#               "parent_ftr": index of parent feature for 'mp' types, this is the feature
#                            index of the annot_type eq "multifeature" of which the current
#                            feature is a child.
#
# Arguments:
#   $nmp:              number of mature peptides, may be 0
#   $ncds:             number of CDS features, may be 0
#   $nxfeat:           number of extra features, modeled independently
#   $ndfeat:           number of duplicate features, modeled by other features
#   $cds2pmatpept_AAR: 1st dim: cds index (-1, off-by-one), 
#                      2nd dim: value array of primary matpept indices that comprise this CDS, 
#                      OR undefined if all features are CDS and there are no mature peptides; 
#                      PRE-FILLED
#   $cds2amatpept_AAR: 1st dim: cds index (-1, off-by-one), 
#                      2nd dim: value array of all matpept indices that comprise this CDS, 
#                      OR undefined if all features are CDS and there are no mature peptides; 
#                      PRE-FILLED
#   $ftr_info_HAR:     ref to hash of arrays to add information to
#   $FH_HR:            REF to hash of file handles, including "log" and "cmd"
#             
# Returns:    void
#
# Dies:       if $ftr_info_HAR is invalid upon entering
#################################################################
sub determineFeatureTypes {
  my $sub_name  = "determineFeatureTypes()";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ncds, $nmp, $nxfeat, $ndfeat, $cds2pmatpept_AAR, $cds2amatpept_AAR, $ftr_info_HAR, $FH_HR) = (@_);

  my $c; # counter

  my $nftr = validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, undef, $FH_HR);

  # initialize arrays
  @{$ftr_info_HAR->{"type"}}                     = ();
  @{$ftr_info_HAR->{"type_idx"}}                 = ();
  @{$ftr_info_HAR->{"annot_type"}}               = ();
  @{$ftr_info_HAR->{"primary_children_ftr_str"}} = ();
  @{$ftr_info_HAR->{"primary_children_ftr_num"}} = ();
  @{$ftr_info_HAR->{"all_children_ftr_str"}}     = ();
  @{$ftr_info_HAR->{"all_children_ftr_num"}}     = ();
  @{$ftr_info_HAR->{"parent_ftr"}}               = ();

  if(! defined $cds2pmatpept_AAR) { 
    # first $ncds features are CDS, remainder (if any) are afeat features
    for($c = 0; $c < $ncds; $c++) { 
      $ftr_info_HAR->{"type"}[$c]                     = "cds-notmp";
      $ftr_info_HAR->{"annot_type"}[$c]               = "model"; # this feature is annotated by homology models
      $ftr_info_HAR->{"type_idx"}[$c]                 = $c+1;
      $ftr_info_HAR->{"primary_children_ftr_str"}[$c] = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"primary_children_ftr_num"}[$c] = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_str"}[$c]     = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_num"}[$c]     = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"parent_ftr"}[$c]               = -1;  # changed to a valid value for certain types (e.g. 'mp')
    }
    for(my $xfeat_idx = 0; $xfeat_idx < $nxfeat; $xfeat_idx++) { 
      $c = $ncds + $xfeat_idx;
      $ftr_info_HAR->{"type"}[$c]                     = "xfeat";
      $ftr_info_HAR->{"annot_type"}[$c]               = "model"; # this feature is annotated by homology models
      $ftr_info_HAR->{"type_idx"}[$c]                 = $xfeat_idx+1;
      $ftr_info_HAR->{"primary_children_ftr_str"}[$c] = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"primary_children_ftr_num"}[$c] = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_str"}[$c]     = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_num"}[$c]     = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"parent_ftr"}[$c]               = -1;  # changed to a valid value for certain types (e.g. 'mp')
    }
    for(my $dfeat_idx = 0; $dfeat_idx < $ndfeat; $dfeat_idx++) { 
      $c = $ncds + $nxfeat + $dfeat_idx;
      $ftr_info_HAR->{"type"}[$c]                     = "dfeat";
      $ftr_info_HAR->{"annot_type"}[$c]               = "duplicate"; # this feature is modeled as a duplicate of another feature
      $ftr_info_HAR->{"type_idx"}[$c]                 = $dfeat_idx+1;
      $ftr_info_HAR->{"primary_children_ftr_str"}[$c] = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"primary_children_ftr_num"}[$c] = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_str"}[$c]     = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_num"}[$c]     = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"parent_ftr"}[$c]               = -1;  # changed to a valid value for certain types (e.g. 'mp')
    }
  }
  else { 
    # some features are mature peptides, some are CDS, some of the CDS are comprised of mature peptides, 
    # and (possibly) some of the CDS are not, and possibly some are xfeat

    # first comes the $nmp mature peptides
    for($c = 0; $c < $nmp; $c++) { 
      $ftr_info_HAR->{"type"}[$c]                     = "mp";
      $ftr_info_HAR->{"annot_type"}[$c]               = "model"; # this feature is annotated by homology models
      $ftr_info_HAR->{"type_idx"}[$c]                 = $c+1;
      $ftr_info_HAR->{"primary_children_ftr_str"}[$c] = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"primary_children_ftr_num"}[$c] = 0;   # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_str"}[$c]     = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_num"}[$c]     = 0;   # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"parent_ftr"}[$c]               = -1;  # changed to a valid value for certain types (e.g. 'mp')
    }

    # then come the CDS, but we need to consult $cds2pmatpept_AAR to see if they are comprised of mature peptides ("cds-mp") 
    # or not ("cds-notmp")
    my $nnotmp_idxed = 0; # we increment this as we index cds-notmp 
    my $nnotmp_extra = 0; # # of cds-notmp that occur after final cds-mp (idx >= scalar(@cds2pmatpept_AA))
    for(my $cds_idx = 0; $cds_idx < $ncds; $cds_idx++) { 
      $c = $nmp + $cds_idx;
      if(defined $cds2pmatpept_AAR->[$cds_idx]) {
        if(! defined $cds2amatpept_AAR->[$cds_idx]) { 
          DNAORG_FAIL("ERROR in $sub_name, cds_idx has 'primary' mature peptide information but not 'all' mature peptide information.", 1, $FH_HR);
        }
        # a CDS that is comprised of mature peptides
        $ftr_info_HAR->{"type"}[$c]                     = "cds-mp";
        $ftr_info_HAR->{"annot_type"}[$c]               = "multifeature"; # this feature is annotated by multiple other features
        $ftr_info_HAR->{"type_idx"}[$c]                 = $cds_idx+1;
        $ftr_info_HAR->{"primary_children_ftr_str"}[$c] = ""; # added to below when we go through $cds2pmatpept_AAR
        $ftr_info_HAR->{"primary_children_ftr_num"}[$c] = 0;  # added to below when we go through $cds2pmatpept_AAR
        $ftr_info_HAR->{"all_children_ftr_str"}[$c]     = ""; # added to below when we go through $cds2amatpept_AAR
        $ftr_info_HAR->{"all_children_ftr_num"}[$c]     = 0;  # added to below when we go through $cds2amatpept_AAR
        $ftr_info_HAR->{"parent_ftr"}[$c]               = -1; # changed to a valid value below

        # step through @{$cds2pmatpept_AAR} and create the primary_children_ftr_str for this CDS
        for(my $z = 0; $z < scalar(@{$cds2pmatpept_AAR->[$cds_idx]}); $z++) { 
          if($ftr_info_HAR->{"primary_children_ftr_str"}[$c] ne "") { 
            $ftr_info_HAR->{"primary_children_ftr_str"}[$c] .= " ";
          }
          my $mp_idx = $cds2pmatpept_AAR->[$cds_idx][$z];
          $ftr_info_HAR->{"primary_children_ftr_str"}[$c] .= $mp_idx;
          $ftr_info_HAR->{"primary_children_ftr_num"}[$c]++;
        }
        # step through @{$cds2amatpept_AAR} and create the all_children_ftr_str for this CDS
        for(my $z = 0; $z < scalar(@{$cds2amatpept_AAR->[$cds_idx]}); $z++) { 
          if($ftr_info_HAR->{"all_children_ftr_str"}[$c] ne "") { 
            $ftr_info_HAR->{"all_children_ftr_str"}[$c] .= " ";
          }
          my $mp_idx = $cds2amatpept_AAR->[$cds_idx][$z];
          $ftr_info_HAR->{"all_children_ftr_str"}[$c] .= $mp_idx;
          $ftr_info_HAR->{"all_children_ftr_num"}[$c]++;
          $ftr_info_HAR->{"parent_ftr"}[$mp_idx] = $c;
        }          
      }
      else { 
        # a CDS that is NOT comprised of mature peptides
        $ftr_info_HAR->{"type"}[$c]                     = "cds-notmp";
        $ftr_info_HAR->{"type_idx"}[$c]                 = $cds_idx+1;
        $ftr_info_HAR->{"annot_type"}[$c]               = "model"; # this feature is annotated by homology models
        $ftr_info_HAR->{"primary_children_ftr_str"}[$c] = "";  # will remain "", only changed for annot_type eq "multifeature"
        $ftr_info_HAR->{"primary_children_ftr_num"}[$c] = 0;   # will remain 0, only changed for annot_type eq "multifeature"
        $ftr_info_HAR->{"all_children_ftr_str"}[$c]     = "";  # will remain "", only changed for annot_type eq "multifeature"
        $ftr_info_HAR->{"all_children_ftr_num"}[$c]     = 0;   # will remain 0, only changed for annot_type eq "multifeature"
        $ftr_info_HAR->{"parent_ftr"}[$c]               = -1;  # changed to a valid value for certain types (e.g. 'mp')
      }
    }
    # add 'xfeat' features, if any
    for(my $xfeat_idx = 0; $xfeat_idx < $nxfeat; $xfeat_idx++) { 
      $c = $nmp + $ncds + $xfeat_idx;
      $ftr_info_HAR->{"type"}[$c]                     = "xfeat";
      $ftr_info_HAR->{"annot_type"}[$c]               = "model"; # this feature is annotated by homology models
      $ftr_info_HAR->{"type_idx"}[$c]                 = $xfeat_idx+1;
      $ftr_info_HAR->{"primary_children_ftr_str"}[$c] = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"primary_children_ftr_num"}[$c] = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_str"}[$c]     = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_num"}[$c]     = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"parent_ftr"}[$c]               = -1;  # changed to a valid value for certain types (e.g. 'mp')
    }
    # add 'dfeat' features, if any
    for(my $dfeat_idx = 0; $dfeat_idx < $ndfeat; $dfeat_idx++) { 
      $c = $nmp + $ncds + $nxfeat + $dfeat_idx;
      $ftr_info_HAR->{"type"}[$c]                     = "dfeat";
      $ftr_info_HAR->{"annot_type"}[$c]               = "duplicate"; # this feature is modeled as a duplicate of another feature
      $ftr_info_HAR->{"type_idx"}[$c]                 = $dfeat_idx+1;
      $ftr_info_HAR->{"primary_children_ftr_str"}[$c] = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"primary_children_ftr_num"}[$c] = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_str"}[$c]     = "";  # will remain "", only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"all_children_ftr_num"}[$c]     = 0;   # will remain 0, only changed for annot_type eq "multifeature"
      $ftr_info_HAR->{"parent_ftr"}[$c]               = -1;  # changed to a valid value for certain types (e.g. 'mp')
    }
  }

  return;
}

#################################################################
# Subroutine: getNumFeaturesAnnotatedByModels()
# Incept:     EPN, Thu Feb 11 14:50:53 2016
#
# Purpose:    Return number of features for which "annot_type" array values are "model".
#
# Arguments: 
#   $ftr_info_HAR: ref to hash of arrays with feature info.
#                  must have at least the "annot_types" array, else we die 
#   $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    number of features for which the has_models array values are 1
#
# Dies:       if $ftr_info_HAR is not valid upon entering.
#################################################################
sub getNumFeaturesAnnotatedByModels { 
  my $sub_name  = "getNumFeaturesAnnotatedByModels()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_HAR, $FH_HR) = (@_);

  my @reqd_ftr_info_A = ("annot_type");
  my $nftr = validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, \@reqd_ftr_info_A, $FH_HR);
  
  my $nftr_with_models = 0;
  for(my $i = 0; $i < $nftr; $i++) { 
    if($ftr_info_HAR->{"annot_type"}[$i] eq "model") { $nftr_with_models++; }
  }
  return $nftr_with_models;
}

#################################################################
# Subroutine: getReferenceFeatureInfo()
# Incept:     EPN, Tue Feb 16 14:05:51 2016
# 
# Purpose:    Fill "ref_strand", "ref_len", "ref_coords" and "out_product"
#             arrays in the feature info hash of arrays (%{$ftr_info_HAR})
#             using $cds_tbl_HHAR and possibly $mp_tbl_HHAR:
#
#             "ref_strand":    strand of this feature in the reference
#             "ref_len":       length (in nt) of this feature in the reference
#             "ref_coords":    coordinates for this feature in the reference in NCBI format
#             "out_product":   'product' value or type_ftable value (if xfeat or dfeat) for this feature, in NCBI annotation
#             "out_gene":      'gene' value or type_ftable value (if xfeat or dfeat) for this feature, in NCBI annotation
#             "out_exception": 'gene' value or type_ftable value (if xfeat or dfeat) for this feature, in NCBI annotation
#             "type_fname":    string for naming output files related to this feature, e.g. "mp", "cds", "$xfeat" or "$dfeat"
#             "type_ftable":   output feature table feature name, e.g. "mat_peptide", "CDS", "$xfeat" or "$dfeat"
# 
# Arguments:
#   $cds_tbl_HHAR:   ref to CDS hash of hash of arrays, PRE-FILLED
#   $mp_tbl_HHAR:    ref to mature peptide hash of hash of arrays, can be undef, else PRE-FILLED
#   $xfeat_tbl_HHAR: ref to extra feature hash of hash of hash of arrays, can be undef, else PRE-FILLED
#   $dfeat_tbl_HHAR: ref to extra feature hash of hash of hash of arrays, can be undef, else PRE-FILLED
#   $ftr_info_HAR:   ref to hash of arrays with feature information, FILLED HERE
#   $ref_accn:       reference accession
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       Never
#
#################################################################
sub getReferenceFeatureInfo { 
  my $sub_name = "getReferenceFeatureInfo()";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($cds_tbl_HHAR, $mp_tbl_HHAR, $xfeat_tbl_HHHAR, $dfeat_tbl_HHHAR, $ftr_info_HAR, $ref_accn, $FH_HR) = @_;

  # initialize the arrays we are about to fill
  @{$ftr_info_HAR->{"ref_strand"}}    = ();
  @{$ftr_info_HAR->{"ref_len"}}       = ();
  @{$ftr_info_HAR->{"ref_coords"}}    = ();
  @{$ftr_info_HAR->{"out_product"}}   = ();
  @{$ftr_info_HAR->{"out_gene"}}      = ();
  @{$ftr_info_HAR->{"out_exception"}} = ();
  @{$ftr_info_HAR->{"type_fname"}}    = ();
  @{$ftr_info_HAR->{"type_ftable"}}   = ();

  my $do_matpept = (defined $mp_tbl_HHAR)     ? 1 : 0;
  my $do_xfeat   = (defined $xfeat_tbl_HHHAR) ? 1 : 0;
  my $do_dfeat   = (defined $dfeat_tbl_HHHAR) ? 1 : 0;

  my $ref_mp_strand_str         = ""; 
  my $ref_cds_strand_str        = "";
  my $ref_xfeat_strand_cur_str  = "";
  my $ref_xfeat_strand_full_str = "";
  my $ref_dfeat_strand_cur_str  = "";
  my $ref_dfeat_strand_full_str = "";
  my $i; # ctr;
  my $nadded_strand;    # number of elements just added to @{$ftr_info_HAR->{"strand"}}
  my $nadded_product;   # number of elements just added to @{$ftr_info_HAR->{"out_product"}}
  my $nadded_gene;      # number of elements just added to @{$ftr_info_HAR->{"out_gene"}}
  my $nadded_exception; # number of elements just added to @{$ftr_info_HAR->{"out_exception"}}

  if($do_matpept) { 
    (undef, undef, undef, undef, undef, $ref_mp_strand_str) = getStrandStats($mp_tbl_HHAR, $ref_accn, $FH_HR);
    getLengthsAndCoords(\%{$mp_tbl_HHAR->{$ref_accn}}, \@{$ftr_info_HAR->{"ref_len"}}, \@{$ftr_info_HAR->{"ref_coords"}}, $FH_HR);
    $nadded_product   = getQualifierValues($mp_tbl_HHAR, $ref_accn, "product",   \@{$ftr_info_HAR->{"out_product"}}, $FH_HR);
    $nadded_gene      = getQualifierValues($mp_tbl_HHAR, $ref_accn, "gene",      \@{$ftr_info_HAR->{"out_gene"}}, $FH_HR);
    $nadded_exception = getQualifierValues($mp_tbl_HHAR, $ref_accn, "exception", \@{$ftr_info_HAR->{"out_exception"}}, $FH_HR);
    $nadded_strand    = length($ref_mp_strand_str);
    for($i = 0; $i < $nadded_strand; $i++) { 
      push(@{$ftr_info_HAR->{"type_fname"}},  "mp"); 
      push(@{$ftr_info_HAR->{"type_ftable"}}, "mat_peptide"); 
    }
    if($nadded_product != $nadded_strand) { 
      DNAORG_FAIL("ERROR in $sub_name, for mature peptides, did not product info for all mature peptides", 1, $FH_HR);
    }
    if(($nadded_gene != 0) && ($nadded_gene != $nadded_strand)) {
      DNAORG_FAIL("ERROR in $sub_name, for mature peptides, added gene info for $nadded_gene mature peptides (should be 0 or $nadded_strand)", 1, $FH_HR);
    }
    if(($nadded_exception != 0) && ($nadded_exception != $nadded_strand)) {
      DNAORG_FAIL("ERROR in $sub_name, for mature peptides, added exception info for $nadded_exception mature peptides (should be 0 or $nadded_strand)", 1, $FH_HR);
    }
    if($nadded_gene == 0) { 
      for($i = 0; $i < $nadded_strand; $i++) { 
        push(@{$ftr_info_HAR->{"out_gene"}}, "");
      }
    }
    if($nadded_exception == 0) { 
      for($i = 0; $i < $nadded_strand; $i++) { 
        push(@{$ftr_info_HAR->{"out_exception"}}, "");
      }
    }
  }

  (undef, undef, undef, undef, undef, $ref_cds_strand_str) = getStrandStats($cds_tbl_HHAR, $ref_accn, $FH_HR);
  getLengthsAndCoords(\%{$cds_tbl_HHAR->{$ref_accn}}, \@{$ftr_info_HAR->{"ref_len"}}, \@{$ftr_info_HAR->{"ref_coords"}}, $FH_HR);
  $nadded_product   = getQualifierValues($cds_tbl_HHAR, $ref_accn, "product",   \@{$ftr_info_HAR->{"out_product"}}, $FH_HR);
  $nadded_gene      = getQualifierValues($cds_tbl_HHAR, $ref_accn, "gene",      \@{$ftr_info_HAR->{"out_gene"}}, $FH_HR);
  $nadded_exception = getQualifierValues($cds_tbl_HHAR, $ref_accn, "exception", \@{$ftr_info_HAR->{"out_exception"}}, $FH_HR);
  $nadded_strand    = length($ref_cds_strand_str);
  for($i = 0; $i < $nadded_strand; $i++) { 
    push(@{$ftr_info_HAR->{"type_fname"}},  "cds");
    push(@{$ftr_info_HAR->{"type_ftable"}}, "CDS");
  }
  if($nadded_product != $nadded_strand) { 
    DNAORG_FAIL("ERROR in $sub_name, for CDS, did not product info for all CDS", 1, $FH_HR);
  }
  if(($nadded_gene != 0) && ($nadded_gene != $nadded_strand)) {
    DNAORG_FAIL("ERROR in $sub_name, for CDS, added gene info for $nadded_gene CDS (should be 0 or $nadded_strand)", 1, $FH_HR);
  }
  if(($nadded_exception != 0) && ($nadded_exception != $nadded_strand)) {
    DNAORG_FAIL("ERROR in $sub_name, for CDS, added exception gene info for $nadded_exception CDS (should be 0 or $nadded_strand)", 1, $FH_HR);
  }
  if($nadded_gene == 0) { 
    for($i = 0; $i < $nadded_strand; $i++) { 
      push(@{$ftr_info_HAR->{"out_gene"}}, "");
    }
  }
  if($nadded_exception == 0) { 
    for($i = 0; $i < $nadded_strand; $i++) { 
      push(@{$ftr_info_HAR->{"out_exception"}}, "");
    }
  }

  if($do_xfeat) {
    foreach my $xfeat (sort keys (%{$xfeat_tbl_HHHAR})) { 
      (undef, undef, undef, undef, undef, $ref_xfeat_strand_cur_str) = getStrandStats(\%{$xfeat_tbl_HHHAR->{$xfeat}}, $ref_accn, $FH_HR);
      $ref_xfeat_strand_full_str .= $ref_xfeat_strand_cur_str;
      getLengthsAndCoords(\%{$xfeat_tbl_HHHAR->{$xfeat}{$ref_accn}}, \@{$ftr_info_HAR->{"ref_len"}}, \@{$ftr_info_HAR->{"ref_coords"}}, $FH_HR);
      $nadded_product   = getQualifierValues(\%{$xfeat_tbl_HHHAR->{$xfeat}}, $ref_accn, "product",   \@{$ftr_info_HAR->{"out_product"}}, $FH_HR);
      $nadded_gene      = getQualifierValues(\%{$xfeat_tbl_HHHAR->{$xfeat}}, $ref_accn, "gene",      \@{$ftr_info_HAR->{"out_gene"}}, $FH_HR);
      $nadded_exception = getQualifierValues(\%{$xfeat_tbl_HHHAR->{$xfeat}}, $ref_accn, "exception", \@{$ftr_info_HAR->{"out_exception"}}, $FH_HR);
      $nadded_strand    = length($ref_xfeat_strand_cur_str);
      for($i = 0; $i < $nadded_strand; $i++) { 
        push(@{$ftr_info_HAR->{"type_fname"}},  $xfeat);
        push(@{$ftr_info_HAR->{"type_ftable"}}, $xfeat);
      }
      if(($nadded_product != 0) && ($nadded_product != $nadded_strand)) {
        DNAORG_FAIL("ERROR in $sub_name, for $xfeat features, added gene info for $nadded_product features (should be 0 or $nadded_strand)", 1, $FH_HR);
      }
      if($nadded_product == 0) { 
        for($i = 0; $i < $nadded_strand; $i++) { 
          push(@{$ftr_info_HAR->{"out_product"}}, "");
        }
      }
      if(($nadded_gene != 0) && ($nadded_gene != $nadded_strand)) {
        DNAORG_FAIL("ERROR in $sub_name, for $xfeat features, added gene info for $nadded_gene features (should be 0 or $nadded_strand)", 1, $FH_HR);
      }
      if(($nadded_exception != 0) && ($nadded_exception != $nadded_strand)) {
        DNAORG_FAIL("ERROR in $sub_name, for $xfeat features, added exception info for $nadded_exception features (should be 0 or $nadded_strand)", 1, $FH_HR);
      }
      if($nadded_gene == 0) { 
        for($i = 0; $i < $nadded_strand; $i++) { 
          push(@{$ftr_info_HAR->{"out_gene"}}, "");
        }
      }
      if($nadded_exception == 0) { 
        for($i = 0; $i < $nadded_strand; $i++) { 
          push(@{$ftr_info_HAR->{"out_exception"}}, "");
        }
      }
    }
  }

  if($do_dfeat) {
    foreach my $dfeat (sort keys (%{$dfeat_tbl_HHHAR})) { 
      (undef, undef, undef, undef, undef, $ref_dfeat_strand_cur_str) = getStrandStats(\%{$dfeat_tbl_HHHAR->{$dfeat}}, $ref_accn, $FH_HR);
      $ref_dfeat_strand_full_str .= $ref_dfeat_strand_cur_str;
      getLengthsAndCoords(\%{$dfeat_tbl_HHHAR->{$dfeat}{$ref_accn}}, \@{$ftr_info_HAR->{"ref_len"}}, \@{$ftr_info_HAR->{"ref_coords"}}, $FH_HR);
      $nadded_product   = getQualifierValues(\%{$dfeat_tbl_HHHAR->{$dfeat}}, $ref_accn, "product",   \@{$ftr_info_HAR->{"out_product"}}, $FH_HR);
      $nadded_gene      = getQualifierValues(\%{$dfeat_tbl_HHHAR->{$dfeat}}, $ref_accn, "gene",      \@{$ftr_info_HAR->{"out_gene"}}, $FH_HR);
      $nadded_exception = getQualifierValues(\%{$dfeat_tbl_HHHAR->{$dfeat}}, $ref_accn, "exception", \@{$ftr_info_HAR->{"out_exception"}}, $FH_HR);
      $nadded_strand    = length($ref_dfeat_strand_cur_str);
      for($i = 0; $i < $nadded_strand; $i++) { 
        push(@{$ftr_info_HAR->{"type_fname"}},  $dfeat);
        push(@{$ftr_info_HAR->{"type_ftable"}}, $dfeat);
      }
      if(($nadded_product != 0) && ($nadded_product != $nadded_strand)) {
        DNAORG_FAIL("ERROR in $sub_name, for $dfeat features, added gene info for $nadded_product features (should be 0 or $nadded_strand)", 1, $FH_HR);
      }
      if($nadded_product == 0) { 
        for($i = 0; $i < $nadded_strand; $i++) { 
          push(@{$ftr_info_HAR->{"out_product"}}, "");
        }
      }
      if(($nadded_gene != 0) && ($nadded_gene != $nadded_strand)) {
        DNAORG_FAIL("ERROR in $sub_name, for $dfeat features, added gene info for $nadded_gene features (should be 0 or $nadded_strand)", 1, $FH_HR);
      }
      if(($nadded_exception != 0) && ($nadded_exception != $nadded_strand)) {
        DNAORG_FAIL("ERROR in $sub_name, for $dfeat features, added exception info for $nadded_exception features (should be 0 or $nadded_strand)", 1, $FH_HR);
      }
      if($nadded_gene == 0) { 
        for($i = 0; $i < $nadded_strand; $i++) { 
          push(@{$ftr_info_HAR->{"out_gene"}}, "");
        }
      }
      if($nadded_exception == 0) { 
        for($i = 0; $i < $nadded_strand; $i++) { 
          push(@{$ftr_info_HAR->{"out_exception"}}, "");
        }
      }
    }
  }


  @{$ftr_info_HAR->{"ref_strand"}} = split("", $ref_mp_strand_str . $ref_cds_strand_str . $ref_xfeat_strand_full_str . $ref_dfeat_strand_full_str);

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
#           As a special case, if --infasta enabled, we don't actually
#           fetch any of the reference sequences, but we do populate
#           the data structures with information on the features
#           and models.
#
#           The following values are filled in the %{$mdl_info_HAR}:
#                "checksum":      checksum of the 'alignment' (single sequence) file the model was built from
#                "cmname":        name of the model, used in infernal output 
#                "length":        length, in nucleotides, of the model
#                "ref_start":     start position of modelled region in the reference genome
#                "ref_stop":      stop position of modelled region in the reference genome
#                "ref_strand":    strand of modelled region in the reference genome
#                "filename_root": 'root' string for output file names related to this model: 
#                "out_tiny":      output value: very short name for this model (e.g. "CDS#4.2")
#                "map_ftr":       the feature index (array index in ftr_info_HAR) this model models
#                "is_final":      '1' if this model is the final model (e.g. final exon) for the feature it models ("map_ftr")
#                "is_first":      '1' if this model is the first model (e.g. final exon) for the feature it models ("map_ftr")
#                "map_exon":      the exon index this model models (1.."map_nexon" value)
#                "map_nexon":     the number of exons the feature this model models has 
#                "out_idx":       output value: feature index and segment index this model (e.g. "4.2")
#
#           The following values are filled in the %{$ftr_info_HAR}:
#                "final_mdl":   index (in arrays of %mdl_info_HA) of final model for this feature
#                "first_mdl":   index (in arrays of %mdl_info_HA) of first model for this feature
#                "nmodels":     number of models for this feature (e.g. number of exons or segments) for this feature, 
#                               0 if 'annotation' eq 'multifeature'
#                "out_short":   output value: short name for this feature (e.g. "CDS #4 [1 exon; +]")
#                "out_tiny":    output value: very short name for this feature (e.g. "CDS#4")
#
# Arguments:
#  $execs_HR:          reference to hash with executables, the key "esl-reformat"
#                      must be defined and the value must be a valid path to an 
#                      esl-reformat executable, PRE-FILLED
#  $sqfile:            Bio::Easel::SqFile object, the sequence file we'll fetch from, already opened by caller
#  $ref_seq_accn:      sequence accession of reference
#  $ref_totlen:        length of reference 
#  $out_root:          root of output file names
#  $build_root:        root of file names created by dnaorg_build.pl
#  $mdl_info_HAR:      ref to hash of arrays with information on the models, FILLED HERE
#  $ftr_info_HAR:      ref to hash of arrays with information on the features, ADDED TO HERE
#  $all_stk_file:      name of output file we will write the stockholm single sequence 
#                      alignments of all features to
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:             REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void; fills @{$values_HHAR}
#
# Dies:       if $ftr_info_HAR is not valid upon entering
#################################################################
sub fetchReferenceFeatureSequences {
  my $sub_name = "fetchReferenceFeatureSequences()";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $sqfile, $ref_seq_accn, $ref_totlen, $out_root, $build_root, $mdl_info_HAR, $ftr_info_HAR, $all_stk_file, $opt_HHR, $FH_HR) = @_;

  # contract check

  # ftr_info_HAR should have array data for keys "ref_coords", "ref_strand"
  my @reqd_ftr_info_A  = ("ref_coords", "ref_strand", "annot_type");
  my $nftr             = validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, \@reqd_ftr_info_A, $FH_HR);
  my $nftr_with_models = getNumFeaturesAnnotatedByModels($ftr_info_HAR, $FH_HR);

  # determine if we read the sequences from a provided fasta file
  # if this is true, we won't fetch the reference sequences
  my $do_infasta   = (opt_Exists("--infasta",   $opt_HHR) && opt_Get("--infasta",   $opt_HHR)) ? 1 : 0;

  my $do_keep     = opt_Get("--keep", $opt_HHR); # should we leave intermediates files on disk, instead of removing them?
  my $do_circular = opt_Get("-c", $opt_HHR); # are we allowing circular genomes?
  my $esl_reformat = $execs_HR->{"esl-reformat"};

  my $cur_out_root;        # output root plus the name of a current file 
  my $cur_out_name_root;   # same as $cur_out_root minus the directory path
  my $out_dir_tail = removeScriptNameFromString(removeDirPath($out_root)); # remove the directory path and any 'script-specific name', e.g. 'dnaorg_build'

  my $cur_build_root;      # build root plus the name of a current file 
  my $cur_build_name_root; # same as $cur_build_root minus the directory path
  my $build_dir_tail = removeScriptNameFromString(removeDirPath($build_root)); # remove the directory path and any 'script-specific name', e.g. 'dnaorg_build'
  my $fetch_input;
  my $fetch_output;
  my $nmdl = 0;               # number of HMMs (and alignments used to build those HMMs)
  my @mdl_A = ();             # [0..$nmdl-1]: array of model names, also name of stockholm alignments used to build those models
  my %mdllen_H          = (); # key: model name from @mdl_A, value is model length
  my @ref_nexons_A      = (); # [0..$c..$ref_nmft-1]: number of exons in CDS or mat_peptide $c+1
  my $ref_tot_nexons    = 0;  # total number of exons in all CDS or mat_peptides
  my @indi_ref_name_A   = (); # [0..$nmdl-1]: name of individual stockholm alignments and models
  my @indi_cksum_stk_A  = (); # [0..$nmdl-1]: checksum's of each named individual stockholm alignment

  my @files2rm_A = ();  # array of file names to remove at end of this function (remains empty if $do_keep)

  for(my $i = 0; $i < $nftr; $i++) { 
    my $type_fname    = $ftr_info_HAR->{"type_fname"}[$i];
    my $do_model      = ($ftr_info_HAR->{"annot_type"}[$i] eq "model") ? 1 : 0;
    my $ftr_type      = $ftr_info_HAR->{"type"}[$i];
    my $ftr_type_idx  = $ftr_info_HAR->{"type_idx"}[$i];
    $ftr_info_HAR->{"first_mdl"}[$i]  = -1; # remains -1 if $do_model is FALSE
    $ftr_info_HAR->{"final_mdl"}[$i]  = -1; # remains -1 if $do_model is FALSE
    $ftr_info_HAR->{"nmodels"}[$i]    = 0;  # remains 0 if $do_model is FALSE
    $ftr_info_HAR->{"append_num"}[$i] = 0;  # maybe changed later in determineFeatureTypes()
    $ftr_info_HAR->{"source_idx"}[$i] = -1; # maybe changed later in determineSourcesOfDuplicateFeatures()

    if($do_model) { 
      # determine start and stop positions of all exons/segments
      my @starts_A = ();
      my @stops_A  = ();
      my $nexons   = 0;
      startsStopsStrandsFromCoordsLength($ftr_info_HAR->{"ref_coords"}[$i], $ref_totlen, $do_circular, \@starts_A, \@stops_A, undef, \$nexons, $FH_HR);
      $ftr_info_HAR->{"nmodels"}[$i] = $nexons;
      my $strand = $ftr_info_HAR->{"ref_strand"}[$i];

      # if we're on the negative strand, reverse the arrays, they'll be in the incorrect order
      if($strand eq "-") { 
        @starts_A = reverse @starts_A;
        @stops_A  = reverse @stops_A;
      }

      for(my $e = 0; $e < $nexons; $e++) { 
        if($nexons > 1) { 
          $cur_out_root        = $out_root .       ".ref." . $type_fname . "." . $ftr_type_idx . ".exon." . ($e+1);
          $cur_out_name_root   = $out_dir_tail .   ".ref." . $type_fname . "." . $ftr_type_idx . ".exon." . ($e+1);
          $cur_build_root      = $build_root .     ".ref." . $type_fname . "." . $ftr_type_idx . ".exon." . ($e+1);
          $cur_build_name_root = $build_dir_tail . ".ref." . $type_fname . "." . $ftr_type_idx . ".exon." . ($e+1);
        }
        else { 
          $cur_out_root        = $out_root .       ".ref." . $type_fname . "." . $ftr_type_idx;
          $cur_out_name_root   = $out_dir_tail .   ".ref." . $type_fname . "." . $ftr_type_idx;
          $cur_build_root      = $build_root .     ".ref." . $type_fname . "." . $ftr_type_idx;
          $cur_build_name_root = $build_dir_tail . ".ref." . $type_fname . "." . $ftr_type_idx;
        }
    
        # determine start and stop of the region we are going to fetch
        my $start = $starts_A[$e];
        my $stop  = $stops_A[$e];
        if($strand eq "-") { # swap start and stop
          my $tmp = $start;
          $start = $stop;
          $stop  = $tmp;
        }

        my @fetch_AA = ();
        my $cksum  = undef;
        my $mdllen = undef;
        if($do_infasta) { 
          # do not actually fetch any sequence
          $cksum  = -1; 
          $mdllen = abs($stop-$start)+1;
        }
        else { 
          # $do_infasta is false, fetch the sequence
          push(@fetch_AA, [$cur_out_name_root, $start, $stop, $ref_seq_accn]);

          # fetch the sequence
          my $cur_fafile = $cur_out_root . ".fa";
          $sqfile->fetch_subseqs(\@fetch_AA, undef, $cur_fafile);
          if(! $do_keep) { push(@files2rm_A, $cur_fafile); }
          
          # reformat to stockholm
          my $cur_stkfile = $cur_out_root . ".stk";
          my $cmd = "$esl_reformat --informat afa stockholm $cur_fafile > $cur_stkfile";
          runCommand($cmd, 0, $FH_HR);
          if(! $do_keep) { push(@files2rm_A, $cur_stkfile); }
    
          # annotate the stockholm file with a blank SS and with a name
          my $cur_named_stkfile = $cur_out_root . ".named.stk";
          ($mdllen, $cksum) = addNameAndBlankSsToStockholmAlignment($cur_out_name_root, 1, $cur_stkfile, $cur_named_stkfile, $FH_HR); # 1: add blank SS_cons line
          if(! $do_keep) { push(@files2rm_A, $cur_named_stkfile); }
          
          # now append the named alignment to the growing stockholm alignment database $all_stk_file
          $cmd = "cat $cur_named_stkfile";
          if($nmdl == 0) { $cmd .= " >  $all_stk_file"; }
          else           { $cmd .= " >> $all_stk_file"; }
          runCommand($cmd, 0, $FH_HR);
        }

        $mdl_info_HAR->{"cmname"}[$nmdl]     = $cur_build_name_root;
        $mdl_info_HAR->{"checksum"}[$nmdl]   = $cksum;
        $mdl_info_HAR->{"length"}[$nmdl]     = $mdllen;
        $mdl_info_HAR->{"ref_start"}[$nmdl]  = $start;
        $mdl_info_HAR->{"ref_stop"}[$nmdl]   = $stop;
        $mdl_info_HAR->{"ref_strand"}[$nmdl] = $strand;

        # store information on this model's name for output purposes
        $mdl_info_HAR->{"filename_root"}[$nmdl] = sprintf("$ftr_type.%s", 
                                                          ($nexons == 1) ? sprintf("%d", $ftr_type_idx) : sprintf("%d.%d", $ftr_type_idx, ($e+1)));
        $mdl_info_HAR->{"out_tiny"}[$nmdl]  = sprintf("%s#%s", 
                                                      $type_fname,
                                                      ($nexons == 1) ? sprintf("%d", $ftr_type_idx) : sprintf("%d.%d", $ftr_type_idx, ($e+1)));
        

        $mdl_info_HAR->{"map_ftr"}[$nmdl]    = $i;
        $mdl_info_HAR->{"is_first"}[$nmdl]   = ($e == 0)           ? 1 : 0;
        $mdl_info_HAR->{"is_final"}[$nmdl]   = ($e == ($nexons-1)) ? 1 : 0;
        $mdl_info_HAR->{"map_exon"}[$nmdl]   = $e;
        $mdl_info_HAR->{"map_nexon"}[$nmdl]  = $nexons;
        $mdl_info_HAR->{"append_num"}[$nmdl] = 0; # maybe changed later in determineFeatureTypes()
        $mdl_info_HAR->{"out_idx"}[$nmdl] = sprintf("%d.%d", 
                                                    $mdl_info_HAR->{"map_ftr"}[$nmdl]+1, $mdl_info_HAR->{"map_exon"}[$nmdl]+1);

        if($e == 0) { 
          $ftr_info_HAR->{"first_mdl"}[$i] = $nmdl;
        }
        if($e == ($nexons-1)) { 
          $ftr_info_HAR->{"final_mdl"}[$i] = $nmdl;
        }
        $nmdl++;
      }
    } # end of 'if($do_model)'
                         
    # store information on this feature's name for output purposes
    my $short = "";
    if($ftr_type eq "mp") { 
      $short = sprintf("MP #%d", $ftr_type_idx);
    }
    elsif($ftr_type eq "cds-mp") { 
      $short = sprintf("CDS(MP) #%d", $ftr_type_idx);
    }
    elsif($ftr_type eq "cds") { 
      $short = sprintf("CDS #%d", $ftr_type_idx);
    }
    else {
      $short = sprintf("%s #%d", $ftr_info_HAR->{"type_ftable"}[$i], $ftr_type_idx);
    }
    my $tiny  = $short;
    $tiny =~ s/\s+//g; # remove whitespace
    if($ftr_info_HAR->{"annot_type"}[$i] eq "model") { 
      if($ftr_info_HAR->{"nmodels"}[$i] > 1) { $short .= sprintf(" [%d %s; %s]", $ftr_info_HAR->{"nmodels"}[$i], featureTypeIsCds($ftr_info_HAR->{"type"}[$i]) ? "exons" : "segments", $ftr_info_HAR->{"ref_strand"}[$i]); }
      else                                   { $short .= sprintf(" [single %s; %s]",  featureTypeIsCds($ftr_info_HAR->{"type"}[$i]) ? "exon" : "segment", $ftr_info_HAR->{"ref_strand"}[$i]); }
    }
    elsif($ftr_info_HAR->{"annot_type"}[$i] eq "multifeature") { 
      if($ftr_type eq "cds-mp") { 
        $short .= sprintf(" [%d mature_peptide(s)]", $ftr_info_HAR->{"primary_children_ftr_num"}[$i]);
      }
      else { 
        DNAORG_FAIL("ERROR in $sub_name, unexpected feature type for multifeature feature: $ftr_type.", 1, $FH_HR);
      }
    }
    $ftr_info_HAR->{"out_tiny"}[$i]      = $tiny;
    $ftr_info_HAR->{"out_short"}[$i]     = $short;
    $ftr_info_HAR->{"filename_root"}[$i] = $ftr_type . "." . $ftr_type_idx;
  }

  # clean up
  foreach my $file2rm (@files2rm_A) { 
    runCommand("rm $file2rm", 0, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine : annotateAppendFeatures()
# Incept:      EPN, Thu Mar  3 16:36:09 2016
#
# Purpose:     Look for and add information about special features
#              and their corresponding models where we want to 
#              annotate *past* the 3' end for the multifeature
#              that is comprised of other features. The specific
#              case we want to handle is to annotate the 3 nt 3'
#              of the end of the final mature peptide of a mature
#              peptide derived CDS, which is the stop codon.
#              
#              Fills for both %{$ftr_info_HAR} and %{$mdl_info_HAR}.
#                "append_num":    number of nucleotides to append for this model, usually 0, unless this
#                                 model is the final mature peptide of a CDS, in which case it is usually 3,
#                                 this informs other parts of the pipeline that we need to append the 3 
#                                 downstream nucleotides: the stop codon, which is not included in the 
#                                 annotation of the final mature peptide in NCBI, but is included in the
#                                 parent CDS annotation.
#
# Arguments: 
#   $ftr_info_HAR: REF to hash of arrays with information on the features, ADDED TO HERE
#   $mdl_info_HAR: REF to hash of arrays with information on the models, ADDED TO HERE
#   $FH_HR:        REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: 
################################################################# 
sub annotateAppendFeatures { 
  my $nargs_expected = 3;
  my $sub_name = "annotateAppendFeatures()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_HAR, $mdl_info_HAR, $FH_HR) = @_;

  my $nftr = getConsistentSizeOfInfoHashOfArrays($ftr_info_HAR, $FH_HR);
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"type"}[$ftr_idx]       eq "cds-mp" && 
       $ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") { 
      # find final feature 
      my @primary_children_idx_A = (); # feature indices of the primary children of this feature
      getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "primary", \@primary_children_idx_A, $FH_HR);
      my $final_ftr_idx = $primary_children_idx_A[(scalar(@primary_children_idx_A)-1)];
      my $final_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$final_ftr_idx];
      $mdl_info_HAR->{"append_num"}[$final_mdl_idx] = 3; # we want to append the 3 nt 3' of this model 
      $ftr_info_HAR->{"append_num"}[$ftr_idx] = 3; # we want to append the 3 nt 3' of this feature
    }
  }

  return;
}


#################################################################
# Subroutine : determineSourcesOfDuplicateFeatures()
# Incept:      EPN, Sun Jul 22 19:51:24 2018
#
# Purpose:     For any features of type 'duplicate', find the
#              other features that will be the source of their
#              annotation boundaries.
#              
#              Fills for %{$ftr_info_HAR}:
#                "source_idx": feature index of the feature that this feature will use
#                              as the source of its annotation boundaries. It will simply
#                              copy those boundaries.
#
# Arguments: 
#   $ftr_info_HAR: REF to hash of arrays with information on the features, ADDED TO HERE
#   $FH_HR:        REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: 
################################################################# 
sub determineSourcesOfDuplicateFeatures { 
  my $nargs_expected = 2;
  my $sub_name = "determineSourcesOfDuplicateFeatures";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_HAR, $FH_HR) = @_;

  my $nftr = getConsistentSizeOfInfoHashOfArrays($ftr_info_HAR, $FH_HR);
  my $match_idx = undef;
  
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"type"}[$ftr_idx]       eq "dfeat" && 
       $ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "duplicate") { 
      # find match to set as the source for this feature's boundary annotation
      $match_idx = undef;
      for(my $ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) {
        if(($ftr_info_HAR->{"type"}[$ftr_idx2] ne "dfeat") && 
           ($ftr_info_HAR->{"ref_coords"}[$ftr_idx2] eq $ftr_info_HAR->{"ref_coords"}[$ftr_idx])) {
          if(defined $match_idx) {
            DNAORG_FAIL("ERROR in $sub_name, more than one feature has identical coordinates (" . $ftr_info_HAR->{"ref_coords"}[$ftr_idx] . ") to a " . $ftr_info_HAR->{"type_ftable"}{$ftr_idx} . " feature", 1, $FH_HR);
          }
          $match_idx = $ftr_idx2;
        }
      }
      if(! defined $match_idx) {
        DNAORG_FAIL("ERROR in $sub_name, zero other features have the same coordinates as the " . $ftr_info_HAR->{"type_ftable"}{$ftr_idx} . " feature requiring a duplicate, with coordinates " . $ftr_info_HAR->{"ref_coords"}[$ftr_idx], 1, $FH_HR);
      }
      $ftr_info_HAR->{"source_idx"}[$ftr_idx] = $match_idx;
    }
    # else $ftr_info_HAR->{"source_idx"}[$ftr_idx] will stay as it was initialized (-1)
  }
  
  return;
}

#################################################################
# Subroutine : annotateOverlapsAndAdjacencies()
# Incept:      EPN, Fri Mar 11 20:24:46 2016
#
# Purpose:     Add information about which models overlap and
#              are adjacent to one another in the reference.
#              
#              Populates the %{$mdl_info_HAR} with the following
#              information (1D keys):
#              "idx_ajb_str":   string of model indices (each separated by a single ",") that are adjacent 'before'
#                               this model index
#              "idx_aja_str":   string of model indices (each separated by a single ",") that are adjacent 'after'
#                               this model index
#              "idx_olp_str":   string of model indices (each separated by a single ",") that overlap by >= 1 nt on
#                               same strand with this model index
#              "out_ajb_str":   string of short model descriptions that are adjacent 'before'
#                               this model index
#              "out_aja_str":   string of short model descriptions that are adjacent 'after'
#                               this model index
#              "out_olp_str":   string of short model descriptions that overlap by >= 1 nt on
#                               same strand with this model index
# Arguments: 
#   $ref_accnlen:  length of the reference accession 
#   $ref_seqlen:   length of the reference sequence (possibly 2*$ref_accnlen)
#   $mdl_info_HAR: REF to hash of arrays with information on the models, ADDED TO HERE
#   $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:        REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: 
################################################################# 
sub annotateOverlapsAndAdjacencies { 
  my $nargs_expected = 5;
  my $sub_name = "annotateOverlapsAndAdjacencies()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ref_accnlen, $ref_seqlen, $mdl_info_HAR, $opt_HHR, $FH_HR) = @_;

  # initialize
  @{$mdl_info_HAR->{"idx_ajb_str"}} = ();
  @{$mdl_info_HAR->{"idx_aja_str"}} = ();
  @{$mdl_info_HAR->{"idx_olp_str"}} = ();
  @{$mdl_info_HAR->{"out_ajb_str"}} = ();
  @{$mdl_info_HAR->{"out_aja_str"}} = ();
  @{$mdl_info_HAR->{"out_olp_str"}} = ();

  overlapsAndAdjacenciesHelper($mdl_info_HAR, 
                               $mdl_info_HAR->{"ref_start"}, $mdl_info_HAR->{"ref_stop"}, $mdl_info_HAR->{"ref_strand"}, 
                               $ref_seqlen, $ref_accnlen,
                               $mdl_info_HAR->{"idx_ajb_str"}, $mdl_info_HAR->{"idx_aja_str"}, $mdl_info_HAR->{"idx_olp_str"}, 
                               $mdl_info_HAR->{"out_ajb_str"}, $mdl_info_HAR->{"out_aja_str"}, $mdl_info_HAR->{"out_olp_str"}, 
                               $opt_HHR, $FH_HR);

  return;
}
#################################################################
#################################################################
#
# Subroutines related to the output info hash:
#   openAndAddFileToOutputInfo()
#   addClosedFileToOutputInfo()
#   helperAddFileToOutputInfo()
#
#################################################################
# Subroutine: openAndAddFileToOutputInfo()
# Incept:     EPN, Fri Feb 26 11:11:09 2016
# 
# Purpose:    Add information about an output file to the
#             %{$ofile_info_HHR} and open that output file. Eventually
#             we'll output information about this file with
#             outputConclusionAndCloseFiles().
#
#             Most of the work is done by helperAddFileToOutputInfo().
#
# Arguments:
#   $ofile_info_HHR:        REF to the 2D hash of output file information, ADDED TO HERE 
#   $key2d:                 2D key for the file we're adding and opening, e.g. "log"
#   $fullpath:              full path to the file we're adding and opening
#   $mainout:               '1' to always output description of this file to 'main' when script ends
#                           '0' to only output a description of this file to the "list" file
#   $desc:                  description of the file we're adding and opening
#
# Returns:    void
# 
# Dies:       If $ofile_info_HHR{*}{$key} already exists.
#
#################################################################
sub openAndAddFileToOutputInfo { 
  my $sub_name = "openAndAddFileToOutputInfo";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($ofile_info_HHR, $key2d, $fullpath, $mainout, $desc) = @_;

  # this helper function does everything but open the file handle
  helperAddFileToOutputInfo($ofile_info_HHR, $key2d, $fullpath, $mainout, $desc);

  # and open the file handle
  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  if(! open($ofile_info_HHR->{"FH"}{$key2d}, ">", $fullpath)) { 
    DNAORG_FAIL("ERROR in $sub_name, unable to open $fullpath for writing.", 1, $FH_HR); 
  }

  return;
}
#################################################################
# Subroutine: addClosedFileToOutputInfo()
# Incept:     EPN, Tue Feb 16 14:22:36 2016
# 
# Purpose:    Add information about a created output file (not open) to
#             the %{$ofile_info_HHR data structure, for eventual
#             output in outputConclusionAndCloseFiles().
#
#             Most of the work is done by helperAddFileToOutputInfo().
#
# Arguments:
#   $ofile_info_HHR:        REF to the 2D hash of output file information, ADDED TO HERE 
#                           for 1D key $key
#   $key2d:                 2D key for the file we're adding and opening, e.g. "fasta"
#   $fullpath:              full path to the closed file we're adding
#   $mainout:               '1' to always output description of this file to 'main' when script ends
#                           '0' to only output a description of this file to the "list" file
#   $desc:                  description of the closed file we're adding 
#
# Returns:    void
# 
# Dies:       If $ofile_desc_HR->{$key} already exists.
#
#################################################################
sub addClosedFileToOutputInfo { 
  my $sub_name = "addClosedFileToOutputInfo()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($ofile_info_HHR, $key2d, $fullpath, $mainout, $desc) = @_;

  # this helper function does everything but set the file handle ("FH") value
  helperAddFileToOutputInfo($ofile_info_HHR, $key2d, $fullpath, $mainout, $desc);

  # set FH value to undef
  $ofile_info_HHR->{"FH"}{$key2d} = undef;

  return;
}

#################################################################
# Subroutine: helperAddFileToOutputInfo()
# Incept:     EPN, Fri Feb 26 14:35:36 2016
# 
# Purpose:    Add information about an output file to the %{$ofile_info_HHR}
#             data structure. Helper function that's called by both 
#             openAndAddFileToOutputInfo() and addClosedFileToOutputInfo().
#             Also, if $ofile_info_HHR->{"FH"}{"list"} is defined, 
#             output the description of this file to the list file.
#
# Arguments:
#   $ofile_info_HHR:        REF to the 2D hash of output file information, ADDED TO HERE 
#                           for 1D key $key
#   $key2d:                 2D key for the file we're adding and opening, e.g. "log"
#   $fullpath:              full path to the file we're adding and opening
#   $mainout:               '1' to always output description of this file to 'main' when script ends
#                           '0' to only output a description of this file to the "list" file
#   $desc:                  description of the file we're adding and opening
#
# Returns:    void
# 
# Dies:       If $ofile_info_HHR{*}{$key} already exists.
#
#################################################################
sub helperAddFileToOutputInfo { 
  my $sub_name = "helperAddFileToOutputInfo";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($ofile_info_HHR, $key2d, $fullpath, $mainout, $desc) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # make sure $mainout value is 0 or 1
  if(($mainout ne "0") && ($mainout ne "1")) { 
    DNAORG_FAIL("ERROR in $sub_name, entered with invalid 'mainout' value of $mainout (should be 0 or 1)", 1, $FH_HR);
  }

  # make sure we don't already have any information for this 2nd dim key $key2d:
  foreach my $key1d (keys (%{$ofile_info_HHR})) { 
    if(exists $ofile_info_HHR->{$key1d}{$key2d}) { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add file $fullpath with key $key2d, but that key already exists for first dim key $key1d", 1, $FH_HR);
    }
  }

  # set the values of the 2D hash
  my $nodirpath = removeDirPath($fullpath);
  my $nidx = (defined $ofile_info_HHR->{"order"}) ? (scalar(keys %{$ofile_info_HHR->{"order"}})) : 0;
  $ofile_info_HHR->{"order"}{$key2d}     = $nidx+1; # first 2d key added will be '1', 2nd will be '2', etc.
  $ofile_info_HHR->{"fullpath"}{$key2d}  = $fullpath;
  $ofile_info_HHR->{"nodirpath"}{$key2d} = $nodirpath;
  $ofile_info_HHR->{"desc"}{$key2d}      = $desc;
  $ofile_info_HHR->{"mainout"}{$key2d}   = $mainout;

  # output the description of this file to the list file
  my $list_FH = ((defined $ofile_info_HHR) && (defined $ofile_info_HHR->{"FH"}) && (exists $ofile_info_HHR->{"FH"}{"list"})) ? 
      $ofile_info_HHR->{"FH"}{"list"} : undef;

  if(defined $list_FH) { 
    my $width_desc = length("# ") + maxLengthScalarValueInHash($ofile_info_HHR->{"desc"}) + length(" saved in:");
    if($width_desc < 80) { $width_desc = 80; }
    outputString($list_FH, 0, sprintf("# %-*s %s\n", $width_desc, $ofile_info_HHR->{"desc"}{$key2d} . " saved in:", $ofile_info_HHR->{"nodirpath"}{$key2d}));
  }

  # validate that we've correctly updated the output info 2D hash
  validateOutputFileInfoHashOfHashes($ofile_info_HHR);

  return;
}

#################################################################
#################################################################
#
# Subroutines related to the error info hash:
#   initializeHardCodedErrorInfoHash()
#   addToErrorInfoHash()
#   setIncompatibilityErrorInfoHash()
#   setRequiredErrorInfoHash()
#   setRequiredErrorInfoHash()
#   setFTableInvalidatedByErrorInfoHash()
#   processFeatureErrorsForFTable()
#   populateFTableNoteOrError()
#
#################################################################
# Subroutine: initializeHardCodedErrorInfoHash()
# Incept:     EPN, Fri Mar  4 12:56:43 2016
#
# Purpose:    Set the initial values in an error info hash,
#             using the hardcoded information in this
#             function.
#
# Arguments:
#   $err_info_HAR:  REF to hash of arrays of error information, FILLED HERE
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $err_info_HAR already has keys upon entering this function
#
#################################################################
sub initializeHardCodedErrorInfoHash { 
  my $sub_name = "initializeHardCodedErrorInfoHash";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($err_info_HAR, $FH_HR) = (@_);

  if(scalar (keys(%{$err_info_HAR})) > 0) { 
    DNAORG_FAIL("ERROR in $sub_name, error info hash of arrays already has at least one key", 1, $FH_HR);
  }

  # add each error code, this function will die if we try to add the same code twice, or if something is wrong 
  # with how we try to add it (args to addToErrorInfoHash don't pass the contract check)

  # errors that are not valid in the feature table: do not affect feature table output
  # nop, b5e, b3e, m5e, m3e, olp, ajb, aja, aji, inp (10)
  addToErrorInfoHash($err_info_HAR, "nop", "feature",  0,
                     "unable to identify homologous feature", # description
                     0, 0, "", "", # feature table info: valid, pred_stop, note, err,
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b5e", "feature",  0,
                     "alignment to reference does not extend to 5' boundary of reference but does extend to 5' boundary of target", # description
                     0, 0, "", "", # feature table info: valid, pred_stop, note, err,
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b3e", "feature",  0,
                     "alignment to reference does not extend to 3' boundary of reference but does extend to 3' boundary of target", # description
                     0, 0, "", "", # feature table info: valid, pred_stop, note, err,
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "m5e", "feature",  0,
                     "parent CDS has a b5e error", # description
                     0, 0, "", "", # feature table info: valid, pred_stop, note, err,
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "m3e", "feature",  0,
                     "parent CDS has a b3e error", # description
                     0, 0, "", "", # feature table info: valid, pred_stop, note, err,
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "olp", "feature",  0,
                     "feature does not overlap with same set of features as in reference", # description
                     0, 0, "", "", # feature table info: valid, pred_stop, note, err,
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "ajb", "feature",  0,
                     "feature (MP or CDS) is not adjacent to same set of features before it as in reference", # description
                     0, 0, "", "", # feature table info: valid, pred_stop, note, err,
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "aja", "feature",  0,
                     "feature (MP or CDS) is not adjacent to same set of features after it as in reference", # description
                     0, 0, "", "", # feature table info: valid, pred_stop, note, err,
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "inp", "feature",  0,
                     "CDS comprised of mat_peptides is incomplete: at least one primary mat_peptide is not identified (nop)", # description
                     # NOTE: currently ftbl_valid == 0, because any feature with this error has 
                     # no valid start/stop and so won't show up in the feature table
                     0, 0, "", "", # feature table info: valid, pred_stop, note, err,
                     $FH_HR);


  # errors that can be invalidated by other errors in feature table output, many of these have to do with premature stop codons
  # nm3, stp, trc, ext, ntr, ctr, int (7)
  addToErrorInfoHash($err_info_HAR, "nm3", "feature",  0,
                     "length of nucleotide feature is not a multiple of 3", # description
                     1, 0, "similar to !out_product,out_gene!; length is not a multiple of 3", # feature table info: valid, pred_stop, note
                     "Unexpected Length: (!out_product,out_gene!) length is not a multiple of 3; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "stp", "feature",  1,
                     "predicted CDS stop by homology is invalid; there may be a valid stop in a different location due to truncation (trc) or extension (ext) (TAG|TAA|TGA)", # description
                     1, 1, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified on !out_product,out_gene!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "trc", "feature",  0,
                     "in-frame stop codon exists 5' of stop position predicted by homology to reference", # description
                     # NOTE: invalidated by int and ctr because int or ctr will handle the feature table note/err
                     1, 1, "similar to !out_product,out_gene!; contains premature stop codon", # feature table info: valid, pred_stop, note
                     "!FEATURE_TYPE! Has Stop Codon: (!out_product,out_gene!) contains unexpected stop codon; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "ext", "feature",  1,
                     "first in-frame stop codon exists 3' of stop position predicted by homology to reference", # description
                     1, 1, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop note
                     "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "ntr", "feature",  0,
                     "mat_peptide may not be translated because its CDS has an in-frame stop 5' of the mat_peptide's predicted start", # description
                     1, 0, "similar to !out_product,out_gene!; polyprotein may not be translated", # feature table info: valid, pred_stop, note
                     "Peptide Translation Problem: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "ctr", "feature",  0,
                     "CDS comprised of mat_peptides includes a mat_peptide that may be truncated (has a trc error)", # description
                     # NOTE: currently valid is 1 and this doesn't invalidate trc, but maybe it should?
                     1, 1, "similar to !out_product,out_gene!; contains premature stop codon", # feature table info: valid, pred_stop, note
                     "!FEATURE_TYPE! Has Stop Codon: (!out_product,out_gene!) contains unexpected stop codon; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "int", "feature",  0,
                     "CDS comprised of mat_peptides is incomplete: at least one primary mat_peptide is not translated due to early stop (ntr)", # description
                     1, 1, "similar to !out_product,out_gene!; contains premature stop codon", # feature table info: valid, pred_stop, note
                     "!FEATURE_TYPE! Has Stop Codon: (!out_product,out_gene!) contains unexpected stop codon; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "nst", "feature",  1,
                     "no in-frame stop codon exists 3' of predicted valid start codon", # description
                     1, 1, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "str", "feature",  0,
                     "predicted CDS start position is not beginning of ATG start codon", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Mutation at Start: (!out_product,out_gene!) expected start codon could not be identified", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b5u", "feature",  0,
                     "alignment to reference does not extend to 5' boundary of reference or target", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b3u", "feature",  0,
                     "alignment to reference does not extend to 3' boundary of reference or target", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "lsc", "feature",  0,
                     "low homology score", # description
                     1, 0, "!DESC!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "xnn", "feature",  0,
                     "blastx identifies protein not identified in nucleotide-based search", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "xnh", "feature",  0,
                     "blastx protein validation failure, no blastx hits", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "xos", "feature",  0,
                     "blastx protein validation failure, strand mismatch between protein and nucleotide predictions", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "x5l", "feature",  0,
                     "blastx protein validation failure, protein alignment extends past nucleotide alignment at 5' end", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "x5s", "feature",  0,
                     "blastx protein validation failure, protein alignment does not extend close enough to nucleotide alignment 5' endpoint", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "x3l", "feature",  0,
                     "blastx protein validation failure, protein alignment extends past nucleotide alignment at 3' end", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation at Stop: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "x3s", "feature",  0,
                     "blastx protein validation failure, protein alignment does not extend close enough to nucleotide alignment 3' endpoint", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation at Stop: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "xin", "feature",  0,
                     "blastx protein validation failure, too large of an insert", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Insertion of Nucleotides: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "xtr", "feature",  0,
                     "blastx protein validation failure, stop codon in protein alignment", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "!FEATURE_TYPE! Has Stop Codon: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "mip", "feature",  0,
                     "mat_peptide may not be translated because its CDS has a blastx protein validation failure", # description
                     1, 0, "similar to !out_product,out_gene!; polyprotein may not be translated", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "dup", "feature",  0,
                     "more than one homologous region", # description
                     1, 0, "!DESC!", # feature table info: valid, pred_stop, note
                     "Duplicate Feature: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "aji", "feature",  0,
                     "CDS comprised of mat_peptides has at least one adjacency inconsistency between 2 primary mat_peptides", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Peptide Adjacency Problem: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "mtr", "feature",  0,
                     "mat_peptide may not be translated because its CDS has a problem", # description
                     1, 0, "similar to !out_product,out_gene!; polyprotein may not be translated", # feature table info: valid, pred_stop, note
                     "Peptide Translation Problem: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "ost", "feature",  0,
                     "predicted feature is on opposite strand from reference", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Reverse Complement: (!out_product,out_gene!) appears to be reverse complemented, sequence may be misassembled", # feature table error
                     $FH_HR);

  # errors that cannot be invalidated by other errors in feature table output
  # ori, str, b5u, b3u, ost, lsc, dup, xip, mpi, xnn, mtr (12)
  addToErrorInfoHash($err_info_HAR, "ori", "sequence", 0, # code, per-type, maybe-allowed
                     "there is not exactly 1 occurrence of origin sequence", # description
                     1, 0, "", "Duplicate Origin: (*sequence*) !DESC!", # feature table info: valid, pred_stop, note, err
                     $FH_HR); 

  # define the incompatibilities; these are two-sided, any error code listed in the 3rd arg is incompatible with the 2nd argument, and vice versa
  setIncompatibilityErrorInfoHash($err_info_HAR, "nop", "aji,int,inp", $FH_HR); 
  setIncompatibilityErrorInfoHash($err_info_HAR, "str", "stp,trc,ext,b5e", $FH_HR);
  setIncompatibilityErrorInfoHash($err_info_HAR, "trc", "ext,nst,aji,inp,b5e", $FH_HR);
  # nop incompatibilities were set as in following line up until 05/16/16, when I realized that
  # a nop in one exon or segment of a multi-model feature (e.g. one exon of a 2-exon feature) 
  # could have a nop and the other could have the other errors. May want to revisit this
  # at some point.
  #  setIncompatibilityErrorInfoHash($err_info_HAR, "nop", "nm3,b5e,b5u,b3e,b3u,str,stp,trc,ext,ntr,nst,aji,int,inp", $FH_HR); # only olp, aja and ajb are compatible with nop

  # define the required combinations, these are one-sided, error code arg 2 requires error code arg 3, but error code arg 3 does not require err code arg 2
  #
  # Previously these were set: 
  # setRequiredErrorInfoHash($err_info_HAR, "ext", "stp", $FH_HR); 
  # setRequiredErrorInfoHash($err_info_HAR, "nst", "stp", $FH_HR);
  #
  # But 'stp' error is reported only if predicted final 3 nts of a CDS are not a valid stop codon,
  # regardless of the frame they are in. 'ext' errors occur if no valid *in-frame* stop codon
  # exists between predicted start and stop, and 'nst' error occurs if no valid *in-frame* stop
  # codon exists 3' of predicted start. So there can be instances of 'ext' and not 'stp', and there
  # can be instances of 'nst' and not 'stp'.

  # define the ftbl_invalid_by values, these are one-sided, any error code listed in the 
  # 3rd argument invalidates the 2nd argument error code, but not vice versa
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "nm3", "b5e,b3e,m5e,m3e",             $FH_HR);
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "ext", "b5e,b3e,m5e,m3e",             $FH_HR);
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "ntr", "b5e,b3e,m5e,m3e",             $FH_HR);
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "ctr", "b5e,b3e,m5e,m3e",             $FH_HR);
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "nst", "b5e,b3e,m5e,m3e",             $FH_HR);

  # trc, ext and nst are preferred to stp
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "stp", "b5e,b3e,m5e,m3e,trc,ext,nst", $FH_HR); 

  # int and ctr are preffered to trc
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "trc", "b5e,b3e,m5e,m3e,int,ctr",     $FH_HR);

  # trc is preffered to ctr
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "ctr", "b5e,b3e,m5e,m3e,trc",         $FH_HR);

  # mip and ntr are preferred to mtr
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "mtr", "mip,ntr", $FH_HR);

  # validate the error info hash
  validateErrorInfoHashIsComplete($err_info_HAR, undef, $FH_HR); 

  return;
}

#################################################################
# Subroutine: addToErrorInfoHash
# Incept:     EPN, Fri Mar  4 13:09:52 2016
#
# Purpose:    Add an element to the error info hash.
#             Die if the same error code already exists.
#
# Arguments:
#   $err_info_HAR:    REF to hash of arrays of error information, FILLED HERE
#   $code:            the code of the element we are adding
#   $pertype:         the 'per-type' of the element we are adding, "sequence" or "feature"
#   $maybe_allowed:   '1' if we're allowing this error to be set as 'maybe', to facilitate
#                     reexamination later.
#   $desc:            the error description/message of the element we are adding
#   $ftbl_valid:      '1' if this error is valid/relevant to feature table output
#                     '0' if it is not valid/relevant
#   $ftbl_pred_stop:  '1' to use predicted stop (instead of corrected one) in 
#                     feature table when this error exists
#   $ftbl_note:       note message for feature table
#                     must eq "" if $ftbl_valid is 0 or per-type is "sequence"
#                     must ne "" if $ftbl_valid is 1 and per-type is "feature"
#   $ftbl_err:        ERROR message for feature table
#                     if $ftbl_valid is 0: must eq ""
#                     if $ftbl_valid is 1: must ne ""
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $err_info_HAR->{"$code"} already exists
#          if $type ne "feature and ne "sequence"
#
#################################################################
sub addToErrorInfoHash { 
  my $sub_name = "addToErrorInfoHash";
  my $nargs_expected = 10;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($err_info_HAR, $code, $pertype, $maybe_allowed, $desc, 
      $ftbl_valid, $ftbl_pred_stop, $ftbl_note, $ftbl_err,
      $FH_HR) = (@_);

  # make sure $pertype is valid
  if(($pertype ne "feature") && ($pertype ne "sequence")) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with per-type $pertype that is not neither \"feature\" nor \"sequence\".", 1, $FH_HR); 
  }
  
  # make sure $maybe_allowed is valid
  if(($maybe_allowed ne "1") && ($maybe_allowed ne "0")) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with invalid maybe_allowed value of $maybe_allowed", 1, $FH_HR);
  }
  
  # make sure $ftbl_valid is valid
  if(($ftbl_valid ne "1") && ($ftbl_valid ne "0")) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with invalid ftbl_valid value of $ftbl_valid", 1, $FH_HR);
  }
  
  # make sure $ftbl_pred_stop is valid
  if(($ftbl_pred_stop ne "1") && ($ftbl_pred_stop ne "0")) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with invalid ftbl_pred_stop value of $ftbl_pred_stop", 1, $FH_HR);
  }
  if(($ftbl_valid eq "0") && ($ftbl_pred_stop eq "1")) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with ftbl_valid of 0, but ftbl_pred_stop of 1", 1, $FH_HR);
  }
  
  # make sure $ftbl_note is valid
  if(! defined $ftbl_note) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add code $code but ftbl_note is undefined", 1, $FH_HR);
  }
  if($pertype eq "sequence") { 
    if($ftbl_note ne "") { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with pertype of sequence but ftbl_note is not empty", 1, $FH_HR);
    }
  }
  elsif($pertype eq "feature") { 
    if(($ftbl_valid eq "0") && ($ftbl_note ne "")) { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with ftbl_valid of 0 and pertype of feature but ftbl_note is not empty", 1, $FH_HR);
    }
    if(($ftbl_valid eq "1") && ($ftbl_note eq "")) { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with ftbl_valid of 1 and pertype of feature but ftbl_note is empty", 1, $FH_HR);
    }
  }

  # make sure $ftbl_err is valid
  # if $ftbl_valid is '0': $ftbl_err must eq ""
  # if $ftbl_valid is '1': $ftbl_err must ne ""
  if(! defined $ftbl_err) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add code $code but ftbl_err is undefined", 1, $FH_HR);
  }
  if($ftbl_valid eq "0") { 
    if($ftbl_err ne "") { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with ftbl_valid of 0 but ftbl_err is not empty", 1, $FH_HR);
    }
  }
  else { # ftbl_valid is not "0"
    if($ftbl_err eq "") { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with ftbl_valid of 1 but ftbl_err is empty", 1, $FH_HR);
    }
  }
  
  # check if $code already exists
  if(exists $err_info_HAR->{"code"}) { 
    my $nerr = scalar(@{$err_info_HAR->{"code"}});
    for(my $err_idx = 0; $err_idx < $nerr; $err_idx++) { 
      my $other_code = $err_info_HAR->{"code"}[$err_idx]; 
      if($code eq $other_code) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, trying to add code $code, but it already exists as element %d in the error info hash", $err_idx+1), 1, $FH_HR);
      }
    }
  }
  
  # initialize, if necessary
  if(! exists $err_info_HAR->{"code"})            { @{$err_info_HAR->{"code"}}            = (); }
  if(! exists $err_info_HAR->{"pertype"})         { @{$err_info_HAR->{"pertype"}}         = (); }
  if(! exists $err_info_HAR->{"maybe_allowed"})   { @{$err_info_HAR->{"maybe_allowed"}}   = (); }
  if(! exists $err_info_HAR->{"desc"})            { @{$err_info_HAR->{"desc"}}            = (); }
  if(! exists $err_info_HAR->{"incompat"})        { @{$err_info_HAR->{"incompat"}}        = (); }
  if(! exists $err_info_HAR->{"requires"})        { @{$err_info_HAR->{"requires"}}        = (); }
  if(! exists $err_info_HAR->{"ftbl_valid"})      { @{$err_info_HAR->{"ftbl_valid"}}      = (); }
  if(! exists $err_info_HAR->{"ftbl_invalid_by"}) { @{$err_info_HAR->{"ftbl_invalid_by"}} = (); }
  if(! exists $err_info_HAR->{"ftbl_pred_stop"})  { @{$err_info_HAR->{"ftbl_pred_stop"}}  = (); }
  if(! exists $err_info_HAR->{"ftbl_note"})       { @{$err_info_HAR->{"ftbl_note"}}       = (); }
  if(! exists $err_info_HAR->{"ftbl_err"})        { @{$err_info_HAR->{"ftbl_err"}}        = (); }
  
  push(@{$err_info_HAR->{"code"}},             $code); 
  push(@{$err_info_HAR->{"pertype"}},          $pertype); 
  push(@{$err_info_HAR->{"maybe_allowed"}},    $maybe_allowed); 
  push(@{$err_info_HAR->{"desc"}},             $desc); 
  push(@{$err_info_HAR->{"incompat"}},         ""); # initialized to no incompatabilities, possibly added to later with setIncompatibilityErrorInfoHash()
  push(@{$err_info_HAR->{"requires"}},         ""); # initialized to no requirements,      possibly added to later with setRequiredErrorInfoHash()
  push(@{$err_info_HAR->{"ftbl_valid"}},       $ftbl_valid);
  push(@{$err_info_HAR->{"ftbl_invalid_by"}},  ""); # initialized to no invalid_by's, possibly added to later with setFTableInvalidatedByErrorInfoHash()
  push(@{$err_info_HAR->{"ftbl_pred_stop"}},   $ftbl_pred_stop);
  push(@{$err_info_HAR->{"ftbl_note"}},        $ftbl_note);
  push(@{$err_info_HAR->{"ftbl_err"}},         $ftbl_err);

  return;
}

#################################################################
# Subroutine: setIncompatibilityErrorInfoHash
# Incept:     EPN, Tue Mar  8 11:18:14 2016
#
# Purpose:    Add to the incompatibility value for an error code $code1 given
#             a string of other error codes $code2. Incompatibilities
#             are bi-directional, so we add an incompatibility between 
#             $code1 and $code2 and between $code2 and $code1.
#
# Arguments:
#   $err_info_HAR:  REF to hash of arrays of error information, FILLED HERE
#   $code1:         the code of the element we are adding incompatibility for
#   $code2str:      the codes $code1 is incompatible with, separated by commas
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if one of the error codes in $code1 or $code2str do not
#          exist in %{$err_info_HAR}.
#
#################################################################
sub setIncompatibilityErrorInfoHash { 
  my $sub_name = "setIncompatibilityErrorInfoHash";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($err_info_HAR, $code1, $code2str, $FH_HR) = (@_);

  my $idx1 = findNonNumericValueInArray($err_info_HAR->{"code"}, $code1, $FH_HR);
  if($idx1 == -1) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add incompatibility for code $code1, but it does not exist in the error info hash", 1, $FH_HR);
  }

  my @code2_A = split(',', $code2str);
  foreach my $code2 (@code2_A) { 
    my $idx2 = findNonNumericValueInArray($err_info_HAR->{"code"}, $code2, $FH_HR);
    if($idx2 == -1) { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add incompatibility between codes $code1 and $code2, but $code2 does not exist in the error info hash", 1, $FH_HR);
    }
    if($idx1 == $idx2) { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add incompatibility between a code and itself: $code1 and $code2", 1, $FH_HR);
    }

    # add ',' if necessary
    if($err_info_HAR->{"incompat"}[$idx1] ne "") { $err_info_HAR->{"incompat"}[$idx1] .= ","; }
    if($err_info_HAR->{"incompat"}[$idx2] ne "") { $err_info_HAR->{"incompat"}[$idx2] .= ","; }

    # this is a bi-directional relationship
    $err_info_HAR->{"incompat"}[$idx1] .= $idx2;
    $err_info_HAR->{"incompat"}[$idx2] .= $idx1;
  }

  return;
}

#################################################################
# Subroutine: setRequiredErrorInfoHash
# Incept:     EPN, Tue Mar  8 11:38:26 2016
#
# Purpose:    Add to the required value for an error code $code1 given
#             a string of other error codes $code2. Required values
#             are uni-directional, so we add only a requirement that
#             $code1 needs $code2, but not that $code2 needs $code1.
#
# Arguments:
#   $err_info_HAR:  REF to hash of arrays of error information, FILLED HERE
#   $code1:         the code of the element we are adding requirement for
#   $code2str:      the codes $code1 requires, separated by a comma
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if one of the error codes in $code1 or $code2str do not
#          exist in %{$err_info_HAR}.
#
#################################################################
sub setRequiredErrorInfoHash { 
  my $sub_name = "setRequiredErrorInfoHash";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($err_info_HAR, $code1, $code2str, $FH_HR) = (@_);

  my $idx1 = findNonNumericValueInArray($err_info_HAR->{"code"}, $code1, $FH_HR);
  if($idx1 == -1) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add requirement for code $code1, but it does not exist in the error info hash", 1, $FH_HR);
  }

  my @code2_A = split(',', $code2str);
  foreach my $code2 (@code2_A) { 
    my $idx2 = findNonNumericValueInArray($err_info_HAR->{"code"}, $code2, $FH_HR);
    if($idx2 == -1) { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add requirement between codes $code1 and $code2, but $code2 does not exist in the error info hash", 1, $FH_HR);
    }
    if($idx1 == $idx2) { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add requirement between a code and itself: $code1 and $code2", 1, $FH_HR);
    }

    # add ',' if necessary
    if($err_info_HAR->{"requires"}[$idx1] ne "") { $err_info_HAR->{"requires"}[$idx1] .= ","; }

    # this is a uni-directional relationship
    $err_info_HAR->{"requires"}[$idx1] .= $idx2;
  }

  return;
}

#################################################################
# Subroutine: setFTableInvalidatedByErrorInfoHash
# Incept:     EPN, Thu Nov  1 10:10:03 2018
#
# Purpose:    Add to the ftbl_invalid_by value for an error code $code1 given
#             a string of other error codes $code2. ftr_invalid_by values 
#             are uni-directional.
#
# Arguments:
#   $err_info_HAR:  REF to hash of arrays of error information, FILLED HERE
#   $code1:         the code of the element we are adding ftbl_invalid_by values for
#   $code2str:      the codes $code1 is invalidated by, separated by commas
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if one of the error codes in $code1 or $code2str do not
#          exist in %{$err_info_HAR}.
#
#################################################################
sub setFTableInvalidatedByErrorInfoHash { 
  my $sub_name = "setFTableInvalidatedByErrorInfoHash";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($err_info_HAR, $code1, $code2str, $FH_HR) = (@_);

  my $idx1 = findNonNumericValueInArray($err_info_HAR->{"code"}, $code1, $FH_HR);
  if($idx1 == -1) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add ftbl_invalid_by for code $code1, but it does not exist in the error info hash", 1, $FH_HR);
  }
  if($err_info_HAR->{"ftbl_valid"}[$idx1] == 0) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add ftbl_invalid_by for code $code1, but its ftbl_valid value is 0", 1, $FH_HR);
  }

  # verify the codes in $code2str
  my @code2_A = split(',', $code2str);
  foreach my $code2 (@code2_A) { 
    my $idx2 = findNonNumericValueInArray($err_info_HAR->{"code"}, $code2, $FH_HR);
    if($idx2 == -1) { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add invalidated by relationship between codes $code1 and $code2, but $code2 does not exist in the error info hash", 1, $FH_HR);
    }
    if($idx1 == $idx2) { 
      DNAORG_FAIL("ERROR in $sub_name, trying to add invalidated by relationship between a code and itself: $code1 and $code2", 1, $FH_HR);
    }
  }

  # set the value
  $err_info_HAR->{"ftbl_invalid_by"}[$idx1] = $code2str;

  return;
}

#################################################################
# Subroutine: processFeatureErrorsForFTable()
# Incept:     EPN, Thu Nov  1 12:10:34 2018
#
# Purpose:    Given a string of errors that correspond to a specific
#             sequence and feature, use the %{$err_info_HAR} and
#             process that string to determine what (if any) notes,
#             and errors should be added to the feature table
#             for this seq/feature pair, also determine if the stop
#             coordinate should be the predicted stop instead of a
#             possibly corrected one.  error exceptions apply for this
#             sequence.  Die if more than one apply, that's supposed
#             to be impossible.
#
# Arguments:
#   $err_code_str:           string of errors, comma separated, can be ""
#   $seq_name:               name of sequence
#   $ftr_idx:                feature index
#   $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
#   $err_info_HAR:           REF to hash of arrays with information on the errors, PRE-FILLED
#   $err_ftr_instances_AHHR: REF to array of 2D hashes with per-feature errors, PRE-FILLED
#   $ret_note_AR:            REF to array of notes, possibly CREATED (not added to) here
#   $ret_error_AR:           REF to array of errors, possibly added to here (not created)
#   $FH_HR:                  REF to hash of file handles, including "log" and "cmd"
# 
# Returns: $do_pred_stop: '1' if predicted stop should be used for this feature in the feature table
#
# Dies: Never
#################################################################
sub processFeatureErrorsForFTable { 
  my $sub_name = "processFeatureErrorsForFTable";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($err_code_str, $seq_name, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, $ret_note_AR, $ret_error_AR, $FH_HR) = (@_);

  if($err_code_str eq "") { 
    return 0; 
  }

  # printf("HEYA in $sub_name $seq_name $ftr_idx, $err_code_str\n");

  # create a hash of all errors in the input $err_str, and also verify they are all valid errors
  my %input_err_code_H = (); # $input_err_code_H{$err_code} = 1 if $err_code is in $err_code_str
  my @err_idx_A = ();
  my $err_code; 
  my $err_idx; 
  foreach $err_code (split(",", $err_code_str)) { 
    $err_idx = findNonNumericValueInArray($err_info_HAR->{"code"}, $err_code, $FH_HR);
    if($err_idx == -1) { 
      DNAORG_FAIL("ERROR in $sub_name, input error of $err_code in string $err_code_str is invalid", 1, $FH_HR);
    }
    $input_err_code_H{$err_code} = 1; 
    push(@err_idx_A, $err_idx);
  }

  my @tmp_note_A = (); # holds all notes
  my $nerr  = scalar(@err_idx_A);
  my $ret_pred_stop = 0;
  my $valid = 0;
  for(my $e = 0; $e < $nerr; $e++) { 
    $err_idx = $err_idx_A[$e];
    # printf("\terr_idx: $err_idx " . $err_info_HAR->{"code"}[$err_idx] . " valid: " . $err_info_HAR->{"ftbl_valid"}[$err_idx] . " checking...\n");
    if($err_info_HAR->{"ftbl_valid"}[$err_idx]) { 
      $valid = 1; # may be set to '0' below
      if($err_info_HAR->{"ftbl_invalid_by"}[$err_idx] ne "") { 
        # printf("\t\tinvalid_by is " . $err_info_HAR->{"ftbl_invalid_by"}[$err_idx] . "\n");
        my @invalid_by_err_code_A = split(",", $err_info_HAR->{"ftbl_invalid_by"}[$err_idx]);
        foreach my $err_code2 (@invalid_by_err_code_A) {
          if(exists $input_err_code_H{$err_code2}) { 
            $valid = 0; # $err_idx is invalidated by $err_code2, which is also present in $err_str
            # printf("\t\t\tinvalidated by $err_code2\n");
          }
        }
      }
      if($valid) { 
        # printf("\t\tvalid\n");
        # valid error that will impact output of feature table
        if($err_info_HAR->{"ftbl_pred_stop"}[$err_idx]) { # if any of the errors are valid and have this set, we return it as 1
          $ret_pred_stop = 1;
        }
        # add notes and errors

        my $note_str = populateFTableNoteOrError("ftbl_note", $err_idx, $seq_name, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, $FH_HR);
        if($note_str ne "") { 
          push(@tmp_note_A, $note_str); # we will prune this array and populate @{$ret_note_AR} before returning
        }

        my $error_str = populateFTableNoteOrError("ftbl_err", $err_idx, $seq_name, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, $FH_HR);
        if($error_str ne "") { 
          # only add the error, if an identical error does not already exist in @{$ret_error_AR}
          my $idx = findNonNumericValueInArray($ret_error_AR, $error_str, $FH_HR);
          if($idx == -1) { 
            push(@{$ret_error_AR}, $error_str); 
          }
        }
      }
    }
  }

  # Create the @{$ret_note_AR}, we do not create the @{$ret_error_AR}
  # this is because @{$ret_note_AR} is per-sequence/feature pair, and caller will treat it as such
  # whereas @{$ret_error_AR} is per-sequence and caller will treat as such
  @{$ret_note_AR} = ();
  my $nnote = scalar(@tmp_note_A);
  my $ret_nnote = 0;
  for(my $n = 0; $n < $nnote; $n++) { 
    my $tmp_note = $tmp_note_A[$n];
    my $keep_flag = 1; # possibly set to 0 below
    # only keep this note if:
    # - we haven't already added an identical string to @{$ret_note_AR}
    # - there is not a superstring of it in @tmp_note_A, that we either 
    #   have already or will eventually add to @{$ret_note_AR}

    # check if we already have an identical string in @{$ret_note_AR}
    for(my $rn = 0; $rn < $ret_nnote; $rn++) { 
      my $ret_note = $ret_note_AR->[$rn];
      if($tmp_note eq $ret_note) { 
        $keep_flag = 0;
      }
    }

    # check if there is a superstring of $tmp_note in @tmp_note_A that we
    # either already have added or will add to @{$ret_note_AR}
    for(my $n2 = 0; $n2 < $nnote; $n2++) { 
      my $tmp_note2 = $tmp_note_A[$n2];
      if(($tmp_note ne $tmp_note2) && 
         ($tmp_note2 =~ m/^$tmp_note/)) { # $tmp_note2 is a superstring of $tmp_note
        $keep_flag = 0;
      }
    }
    
    if($keep_flag) { 
      push(@{$ret_note_AR}, $tmp_note);
      $ret_nnote++;
    }
  }

  return $ret_pred_stop;
}

#################################################################
# Subroutine: populateFTableNoteOrError
# Incept:     EPN, Thu Feb  8 14:31:16 2018
#
# Purpose:    Create notes and errors for the feature table for a specific
#             error, feature, and sequence trio.
#
# Arguments:
#   $ekey:                   either "ftbl_note" or "ftbl_err"
#   $err_idx:                index of current error in %{$err_info_HAR} arrays
#   $seq_name:               name of sequence
#   $ftr_idx:                feature index
#   $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
#   $err_info_HAR:           REF to hash of arrays with information on the errors, PRE-FILLED
#   $err_ftr_instances_AHHR: REF to array of 2D hashes with per-feature errors, PRE-FILLED
#   $FH_HR:                    REF to hash of file handles, including "log" 
#                              and "cmd"
# 
# Returns: string with the feature table note for the current sequence/feature combo
#
# Dies:    If ftbl_err_exceptions_HR doesn't have information we need
#          or has invalid information
#
#################################################################
sub populateFTableNoteOrError { 
  my $sub_name = "populateFTableNoteOrError";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ekey, $err_idx, $seq_name, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, $FH_HR) = (@_);
  
  if(! exists $err_info_HAR->{$ekey}) { 
    DNAORG_FAIL("ERROR in $sub_name, $ekey value is undefined in error info hash", 1, $FH_HR);
  }

  my $msg = $err_info_HAR->{$ekey}[$err_idx];

  if(! defined $msg) { 
    DNAORG_FAIL("ERROR in $sub_name, error $err_idx is invalid in error info hash", 1, $FH_HR);
  }

  if($msg eq "") { 
    return "";
  }

  my $err_code = $err_info_HAR->{"code"}[$err_idx];

  my $orig_msg = $msg;
  my $ret_msg  = $msg;
  my $idx;
  # replace !DESC! with description of the error
  if($ret_msg =~ /!DESC!/) { 
    if(exists $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}) { 
      my $desc_str = sprintf("%s%s", 
                             $err_info_HAR->{"desc"}[$err_idx], 
                             ($err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} eq "") ? "" : " [" . $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} . "]"); 
      $ret_msg =~ s/!DESC!/$desc_str/g;
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, trying to return $ekey message for $err_code and sequence $seq_name feature $ftr_idx, but no error instance exists", 1, $FH_HR);
    }
  }
  # replace !FEATURE_TYPE! with 
  if($ret_msg =~ /!FEATURE_TYPE!/) { 
    my $feature_type_str = $ftr_info_HAR->{"type_ftable"}[$ftr_idx];
    $ret_msg =~ s/!FEATURE_TYPE!/$feature_type_str/g;
  }
  # check if there is an internal !$key_str! string, where $key_str is either $key or $key_1,$key_2,...,$key_n for some number n,
  # which is replaced by the value: $ftr_info_HAR->{$key}[$ftr_idx]); for the first $key with a valid value
  if($ret_msg =~ /\!([^\!]*)\!/) {
    my $key_str = $1; 
    my @value_A = split(",", $key_str); 
    my $nvalue = scalar(@value_A);
    my $value = "?";
    for(my $v = 0; $v < $nvalue; $v++) { 
      my $key = $value_A[$v];
      if((exists $ftr_info_HAR->{$key}[$ftr_idx]) && ($ftr_info_HAR->{$key}[$ftr_idx] ne "")) { 
        $value = $ftr_info_HAR->{$key}[$ftr_idx];
        $v = $nvalue; # breaks loop
      }
    }
    $ret_msg =~ s/\![^\!]*\!/$value/;
  }

  return $ret_msg;
}


#################################################################
#
# Subroutines related to overlaps and adjacencies:
#   overlapsAndAdjacenciesHelper()
#   getOverlap()
#   compareTwoOverlapOrAdjacencyIndexStrings()
#   checkForIndexInOverlapOrAdjacencyIndexString()
#
#################################################################
# Subroutine : overlapsAndAdjacenciesHelper()
# Incept:      EPN, Sat Mar 12 10:27:59 2016
#
# Purpose:     Helper function for annotateOverlapsAndAdjacencies()
#              and dnaorg_annotate.pl:calculate_results_overlaps_and_adjacencies().
#
#              Given refs to three arrays, one with start positions,
#              one with stop positions and one with strands, 
#              determine which indices in the arrays are adjacent
#              to each other and overlap with each other
#              are adjacent to one another in the reference.
#              
# Arguments: 
#   $mdl_info_HAR:   REF to hash of arrays with information on the models, PRE-FILLED
#                    used only to get "out_idx" values to use for $out_*_str_AR values.
#   $start_AR:       REF to array with start positions 
#   $stop_AR:        REF to array with stop positions 
#   $strand_AR:      REF to array with strands
#   $seq_len:        total length of the sequence in our fetched sequence file, so we can 
#                    check if we're adjacent $len..1 if -c enabled
#   $accn_len:       length of the accession in GenBank (maybe half of $seq_len if -c used),
#                    we need this so we make a special check for overlaps if -c enabled
#   $idx_ajb_str_AR: REF to array to fill with strings of 'before' adjacency model indices
#   $idx_aja_str_AR: REF to array to fill with strings of 'after' adjacency model indices
#   $idx_olp_str_AR: REF to array to fill with strings of overlaps model indices
#   $out_ajb_str_AR: REF to array to fill with strings of 'before' adjacency model descriptions
#   $out_aja_str_AR: REF to array to fill with strings of 'after' adjacency model descriptions
#   $out_olp_str_AR: REF to array to fill with strings of overlaps
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:          REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: 
################################################################# 
sub overlapsAndAdjacenciesHelper { 
  my $nargs_expected = 14;
  my $sub_name = "overlapsAndAdjacenciesHelper()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($mdl_info_HAR, $start_AR, $stop_AR, $strand_AR, $seq_len, $accn_len, $idx_ajb_str_AR, $idx_aja_str_AR, $idx_olp_str_AR, 
      $out_ajb_str_AR, $out_aja_str_AR, $out_olp_str_AR, $opt_HHR, $FH_HR) = @_;

  my $nmdl = scalar(@{$start_AR});

  # the 2D arrays that hold info on which models overlap with each other and are adjacent to each other
  my @adj_AA = (); # [0..$i..$nmdl-1][0..$j..$nmdl-1], value is '1' if $i and $j are adjacent, else 0 (also, 0 if $i==$j)
  my @olp_AA = (); # [0..$i..$nmdl-1][0..$j..$nmdl-1], value is '1' if $i and $j overlap, else 0 (also, 0 if $i==$j)
  # diagonal values ($i == $j) for @adj_AA and @olp_AA are initialized to '0' and are 
  # skipped in main loop below which fills the matrices for all other values

  my $mdl_idx1; # counter over models
  my $mdl_idx2; # counter over models

  # initialize
  for($mdl_idx1 = 0; $mdl_idx1 < $nmdl; $mdl_idx1++) { 
    for($mdl_idx2 = 0; $mdl_idx2 < $nmdl; $mdl_idx2++) { 
      $adj_AA[$mdl_idx1][$mdl_idx2] = 0;
      $olp_AA[$mdl_idx1][$mdl_idx2] = 0;
    }
  }

  # fill the matrices
  for($mdl_idx1 = 0; $mdl_idx1 < $nmdl; $mdl_idx1++) { 
    my $start1  = $start_AR->[$mdl_idx1];
    my $stop1   = $stop_AR->[$mdl_idx1];
    my $strand1 = $strand_AR->[$mdl_idx1];
    if($start1 != -1) { # a flag that we don't have a start..stop, strand for this model
      for($mdl_idx2 = 0; $mdl_idx2 < $nmdl; $mdl_idx2++) { 
        my $start2  = $start_AR->[$mdl_idx2];
        my $stop2   = $stop_AR->[$mdl_idx2];
        my $strand2 = $strand_AR->[$mdl_idx2];
        #printf("mdl_idx1 $start1..$stop1 $strand1  mdl_idx2 $start2..$stop2 $strand2\n");
        if($start2 != -1 && ($mdl_idx1 != $mdl_idx2)) { 
          # we'll skip this if $mdl_idx1 == $mdl_idx2
          # $start2 == -1 is a flag that we don't have a start..stop, strand for this model
          if($strand1 eq $strand2) { 
            # strands match
            if((($strand1 eq "+") && (($stop1+1) == $start2)) || 
               (($strand1 eq "-") && (($stop1-1) == $start2)) || 
               (($strand1 eq "+") && (($stop2+1) == $start1)) || 
               (($strand1 eq "-") && (($stop2-1) == $start1))) { 
              # $mdl_idx1 and $mdl_idx2 are adjacent
              $adj_AA[$mdl_idx1][$mdl_idx2] = 1;
              #printf("HEYAA $mdl_idx1 and $mdl_idx2 are adjacent\n");
            }
            if((opt_Get("-c", $opt_HHR)) && 
               ((($strand1 eq "+") && ($stop1 == $seq_len) && ($start2 == 1)) || 
                (($strand1 eq "-") && ($stop1 == 1)        && ($start2 == $seq_len)) || 
                (($strand1 eq "+") && ($stop2 == $seq_len) && ($start1 == 1)) || 
                (($strand1 eq "-") && ($stop2 == 1)        && ($start1 == $seq_len)))) { 
              # $mdl_idx1 and $mdl_idx2 are adjacent, one ends the seq, the other begins it
              $adj_AA[$mdl_idx1][$mdl_idx2] = 1;
            }            
            # determine if we have an overlap
            # if -c enabled we have to do potentially 2 calls
            my $nres_overlap1 = 0;
            my $nres_overlap2 = 0;
            $nres_overlap1 = 
                getOverlap(($start1 < $stop1) ? $start1 : $stop1,
                           ($start1 < $stop1) ? $stop1  : $start1,
                           ($start2 < $stop2) ? $start2 : $stop2,
                           ($start2 < $stop2) ? $stop2  : $start2, $FH_HR);
            if((opt_Get("-c", $opt_HHR)) && $nres_overlap1 == 0) { 
              # do a second call iff only one of either start1..stop1
              # or start2..stop2 straddles $accn_len
              # example of why we need to do this is HBV ref NC_003977, CDS1 and CDS4:
              # search results are CDS1: 2309..4807, CDS4: 157..837
              # accn_len is 3182
              my $straddles1 = (($start1 < $accn_len && $stop1 > $accn_len) || 
                                ($start1 > $accn_len && $stop1 < $accn_len)) ? 1 : 0;
              my $straddles2 = (($start2 < $accn_len && $stop2 > $accn_len) || 
                                ($start2 > $accn_len && $stop2 < $accn_len)) ? 1 : 0;
              if($straddles1 && (! $straddles2)) { # start1..stop1 straddles $accn_len, but not start2..stop2
                my $adj_start2 = $start2;
                my $adj_stop2  = $stop2;
                if($start2 < $accn_len) { # this means $stop2 < $accn_len also
                  $adj_start2 += $accn_len;
                  $adj_stop2  += $accn_len;
                }
                else { # start2 and stop2 are > accn_len
                  $adj_start2 -= $accn_len;
                  $adj_stop2  -= $accn_len;
                }
                $nres_overlap2 = 
                    getOverlap(($start1 < $stop1) ? $start1 : $stop1,
                               ($start1 < $stop1) ? $stop1  : $start1,
                               ($adj_start2 < $adj_stop2) ? $adj_start2 : $adj_stop2,
                               ($adj_start2 < $adj_stop2) ? $adj_stop2  : $adj_start2, $FH_HR);
              }
              elsif((! $straddles1) && $straddles2) { # start1..stop1 straddles $accn_len, but not start2..stop2
                my $adj_start1 = $start1;
                my $adj_stop1  = $stop1;
                if($start1 < $accn_len) { # this means $stop1 < $accn_len also
                  $adj_start1 += $accn_len;
                  $adj_stop1  += $accn_len;
                }
                else { 
                  $adj_start1 -= $accn_len;
                  $adj_stop1  -= $accn_len;
                }
                $nres_overlap2 = 
                    getOverlap(($adj_start1 < $adj_stop1) ? $adj_start1 : $adj_stop1,
                               ($adj_start1 < $adj_stop1) ? $adj_stop1  : $adj_start1,
                               ($start2 < $stop2) ? $start2 : $stop2,
                               ($start2 < $stop2) ? $stop2  : $start2, $FH_HR);
              }
            } # end of if statement entered if -c enabled and nres_overlap1 == 0
            if(($nres_overlap1 > 0) || ($nres_overlap2 > 0)) { 
              $olp_AA[$mdl_idx1][$mdl_idx2] = 1;
              #printf("HEYAA $mdl_idx1 and $mdl_idx2 overlap\n");
            }
          }
        }
      }
    }
  }
  # make sure the matrices are symmetric
  for($mdl_idx1 = 0; $mdl_idx1 < $nmdl; $mdl_idx1++) { 
    for($mdl_idx2 = 0; $mdl_idx2 < $nmdl; $mdl_idx2++) { 
      if($adj_AA[$mdl_idx1][$mdl_idx2] != $adj_AA[$mdl_idx2][$mdl_idx1]) { 
        DNAORG_FAIL("ERROR in $sub_name, problem filling adjacency matrix, it is not symmetric [$mdl_idx1][$mdl_idx2] ($adj_AA[$mdl_idx1][$mdl_idx2]) != [$mdl_idx2][$mdl_idx1] ($adj_AA[$mdl_idx2][$mdl_idx1]", 1, $FH_HR);
      }
      if($olp_AA[$mdl_idx1][$mdl_idx2] != $olp_AA[$mdl_idx2][$mdl_idx1]) { 
        DNAORG_FAIL("ERROR in $sub_name, problem filling overlap matrix, it is not symmetric [$mdl_idx1][$mdl_idx2] ($olp_AA[$mdl_idx1][$mdl_idx2]) != [$mdl_idx2][$mdl_idx1] ($olp_AA[$mdl_idx2][$mdl_idx1]", 1, $FH_HR);
      }
    }
  }

  # convert 2D arrays into a series of strings, one per model
  for($mdl_idx1 = 0; $mdl_idx1 < $nmdl; $mdl_idx1++) { 
    $idx_ajb_str_AR->[$mdl_idx1] = "";
    $idx_aja_str_AR->[$mdl_idx1] = "";
    $idx_olp_str_AR->[$mdl_idx1] = "";
    $out_ajb_str_AR->[$mdl_idx1] = "";
    $out_aja_str_AR->[$mdl_idx1] = "";
    $out_olp_str_AR->[$mdl_idx1] = "";
    my $mdl_out1 = $mdl_info_HAR->{"out_idx"}[$mdl_idx1];
    for($mdl_idx2 = 0; $mdl_idx2 < $nmdl; $mdl_idx2++) { 
      my $mdl_out2 = $mdl_info_HAR->{"out_idx"}[$mdl_idx2];
      if($adj_AA[$mdl_idx1][$mdl_idx2]) { 
        if($mdl_idx1 < $mdl_idx2) { 
          $idx_aja_str_AR->[$mdl_idx1] .= ($idx_aja_str_AR->[$mdl_idx1] eq "") ? $mdl_idx2 : "," . $mdl_idx2;
          $out_aja_str_AR->[$mdl_idx1] .= ($out_aja_str_AR->[$mdl_idx1] eq "") ? $mdl_out2 : "," . $mdl_out2;
        }
        else { 
          $idx_ajb_str_AR->[$mdl_idx1] .= ($idx_ajb_str_AR->[$mdl_idx1] eq "") ? $mdl_idx2 : "," . $mdl_idx2;
          $out_ajb_str_AR->[$mdl_idx1] .= ($out_ajb_str_AR->[$mdl_idx1] eq "") ? $mdl_out2 : "," . $mdl_out2;
        }
      }
      if($olp_AA[$mdl_idx1][$mdl_idx2]) { 
        $idx_olp_str_AR->[$mdl_idx1] .= ($idx_olp_str_AR->[$mdl_idx1] eq "") ? $mdl_idx2 : "," . $mdl_idx2;
        $out_olp_str_AR->[$mdl_idx1] .= ($out_olp_str_AR->[$mdl_idx1] eq "") ? $mdl_out2 : "," . $mdl_out2;
      }
    }
  }

  return;
}

#################################################################
# Subroutine: getOverlap()
# Incept:     EPN, Mon Mar 14 13:47:57 2016
#
# Purpose:    Calculate number of nucleotides of overlap between
#             two regions.
#
# Args:
#  $start1: start position of hit 1 (must be <= $end1)
#  $end1:   end   position of hit 1 (must be >= $end1)
#  $start2: start position of hit 2 (must be <= $end2)
#  $end2:   end   position of hit 2 (must be >= $end2)
#  $FH_HR:  REF to hash of file handles, including "log" and "cmd"
#
# Returns:  Number of nucleotides of overlap between hit1 and hit2,
#           0 if none
#
# Dies:     if $end1 < $start1 or $end2 < $start2.

sub getOverlap {
  my $sub_name = "getOverlap";
  my $nargs_exp = 5;
  if(scalar(@_) != 5) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($start1, $end1, $start2, $end2, $FH_HR) = @_; 

  # printf("in $sub_name $start1..$end1 $start2..$end2\n");

  if($start1 > $end1) { DNAORG_FAIL("ERROR in $sub_name start1 > end1 ($start1 > $end1)", 1, $FH_HR); }
  if($start2 > $end2) { DNAORG_FAIL("ERROR in $sub_name start2 > end2 ($start2 > $end2)", 1, $FH_HR); }

  # Given: $start1 <= $end1 and $start2 <= $end2.
  
  # Swap if nec so that $start1 <= $start2.
  if($start1 > $start2) { 
    my $tmp;
    $tmp   = $start1; $start1 = $start2; $start2 = $tmp;
    $tmp   =   $end1;   $end1 =   $end2;   $end2 = $tmp;
  }
  
  # 3 possible cases:
  # Case 1. $start1 <=   $end1 <  $start2 <=   $end2  Overlap is 0
  # Case 2. $start1 <= $start2 <=   $end1 <    $end2  
  # Case 3. $start1 <= $start2 <=   $end2 <=   $end1
  if($end1 < $start2) { return 0; }                      # case 1
  if($end1 <   $end2) { return ($end1 - $start2 + 1); }  # case 2
  if($end2 <=  $end1) { return ($end2 - $start2 + 1); }  # case 3
  die "ERROR in $sub_name, unforeseen case in $start1..$end1 and $start2..$end2";

  return; # NOT REACHED
}

#################################################################
# Subroutine: compareTwoOverlapOrAdjacencyIndexStrings()
# Incept:     EPN, Sat Mar 12 12:40:54 2016
#
# Purpose:    Given two strings that represent adjacencies: 
#             indices separated by commas, compare them and 
#             return an array indicating differences.
#
# Arguments:
#   $str1:    first string, format "<idx1>,<idx2>,...,<idxN>"
#   $str2:    second string
#   $max:     max possible index that can be in $str1 and $str2
#   $diff_AR: ref to array to fill with differences
#             $diff_AR->[$idx] = 0  means $idx exists in both $str1 and $str2
#             $diff_AR->[$idx] = -1 means $idx exists only in $str1
#             $diff_AR->[$idx] = 1  means $idx exists only in $str2
#   $FH_HR:  REF to hash of file handles, including "log" and "cmd"
#
# Returns:  void, fills @{$diff_AR}
# 
# Dies:     If an index in $str1 or $str2 exceeds $max or is not
#           a valid integer >= 0.
#
#################################################################
sub compareTwoOverlapOrAdjacencyIndexStrings { 
  my $sub_name = "compareTwoOverlapOrAdjacencyIndexStrings";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($str1, $str2, $max, $diff_AR, $FH_HR) = @_;

  my $el; # an array element

  # turn them into arrays
  my @str1_A = split(",", $str1);
  my @str2_A = split(",", $str2);
  
  @{$diff_AR} = ();
  for(my $i = 0; $i <= $max; $i++) { 
    $diff_AR->[$i] = 0;
  }

  foreach $el (@str1_A) { 
    if(! verify_integer($el)) { 
      DNAORG_FAIL("ERROR in $sub_name, element $el in $str1 is not an integer", 1, $FH_HR);
    }
    if($el < 0) {  
      DNAORG_FAIL("ERROR in $sub_name, element $el in $str1 is negative, should be between 0 and $max", 1, $FH_HR);
    }
    if($el > $max) {  
      DNAORG_FAIL("ERROR in $sub_name, element $el in $str1 is above $max, should be between 0 and $max", 1, $FH_HR);
    }
    $diff_AR->[$el]--;
  }
  foreach $el (@str2_A) { 
    if(! verify_integer($el)) { 
      DNAORG_FAIL("ERROR in $sub_name, element $el in $str2 is not an integer", 1, $FH_HR);
    }
    if($el < 0) {  
      DNAORG_FAIL("ERROR in $sub_name, element $el in $str2 is negative, should be between 0 and $max", 1, $FH_HR);
    }
    if($el > $max) {  
      DNAORG_FAIL("ERROR in $sub_name, element $el in $str2 is above $max, should be between 0 and $max", 1, $FH_HR);
    }
    $diff_AR->[$el]++;
  }
  # now $diff_AR->[$el] == 0  if either $el doesn't exist in both of $str1 and $str2, or exists in both of them
  #     $diff_AR->[$el] == -1 if $el only exists $str1
  #     $diff_AR->[$el] == 1  if $el only exists $str2

  return;
}

#################################################################
# Subroutine: checkForIndexInOverlapOrAdjacencyIndexString
# Incept:     EPN, Tue Mar 15 10:20:22 2016
#
# Purpose:    Given a 'index' string that describes overlaps
#             or adjacencies, check if the index $idx is in it.
#             Return '1' if it is, else return '0'.
#
# Arguments:
#   $str:     string, format "<idx1>,<idx2>,...,<idxN>"
#   $idx:     index we're checking for
#   $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#
# Returns:  '1' if $idx is in $str, else '0'
# 
# Dies:     If $idx is not integer >= 0
#
#################################################################
sub checkForIndexInOverlapOrAdjacencyIndexString { 
  my $sub_name = "checkForIndexInOverlapOrAdjacencyIndexString()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($str, $idx, $FH_HR) = @_;

  if((! verify_integer($idx)) || ($idx < 0)) { 
    DNAORG_FAIL("ERROR in $sub_name, expected index to be integer greater than or equal to 0, got $idx", 1, $FH_HR);
  }
  if($str eq "") { 
    return 0;
  }
  if($str =~ m/\,$idx\,/) { 
    return 1;
  }
  if($str =~ m/^$idx$/) { 
    return 1;
  }
  if($str =~ m/^$idx\,/) { 
    return 1;
  }
  if($str =~ m/\,$idx$/) { 
    return 1;
  }

  return 0;
}

#################################################################
#################################################################
#
# Massive wrapper subroutines that call other subroutines:
#   wrapperGetInfoUsingEdirect()
#   wrapperFetchAllSequencesAndProcessReferenceSequence()
#
#################################################################
# Subroutine:  wrapperGetInfoUsingEdirect()
# Incept:      EPN, Tue Feb 23 13:00:23 2016
#
# Purpose:     A large block of code that is called once by each
#              dnaorg_build.pl and dnaorg_annotate.pl to gather
#              sequence information using Edirect and parse that
#              Edirect output into usable data structures, by
#              doing the following:
#
#              1) create the edirect .mat_peptide file, if necessary
#              2) create the edirect .ftable file
#              3) create the length file
#              4) parse the edirect .mat_peptide file, if necessary
#              5) parse the edirect .ftable file
#              6) parse the length file
#
#              Creates the following output files and stores
#              information on them in %{$ofile_info_HHR}
#              by calling the addClosedFileToOutputInfo() function:
#              - $out_root . ".mat_peptide": mature peptide info obtained via edirect
#              - $out_root . ".ftable":      feature table obtained via edirect
#              - $out_root . ".length":      sequence length file
#                      
#              As a special case, if $opt_HHR->{"--skipedirect"} is 'on',
#              then we skip the edirect steps and use pre-existing files.
#              
# Arguments: 
#   $listfile:              name of list file with all accessions, can be undef, in which case we 
#                           only gather information for the reference
#   $ref_accn:              reference accession, first accession in $listfile (although this is 
#                           not enforced here, caller enforced it)
#   $out_root:              string that is the 'root' for naming output files
#   $cds_tbl_HHAR:          REF to CDS hash of hash of arrays, FILLED HERE
#   $mp_tbl_HHAR:           REF to mature peptide hash of hash of arrays, can be undef, else FILLED HERE
#   $xfeat_tbl_HHHAR:       REF to hash of hash of hash of arrays for extra qualifier info, can be undef, else FILLED HERE
#   $dfeat_tbl_HHHAR:       REF to hash of hash of hash of arrays for duplicate qualifier info, can be undef, else FILLED HERE
#   $seq_info_HAR:          REF to 2D hash with sequence information, ADDED TO HERE
#   $ofile_info_HHR:        REF to 2D hash with output info, ADDED TO HERE
#   $opt_HHR:               REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#   $FH_HR:                 REF to hash of file handles, including "log" and "cmd", can be undef, PRE-FILLED
#
# Returns:     void
#
# Dies:        if $listfile is defined but it does not exist or is empty
#              if any of the edirect commands exit in error
#     
################################################################# 
sub wrapperGetInfoUsingEdirect {
  my $sub_name = "wrapperGetInfoUsingEdirect";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name, entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($listfile, $ref_accn, $out_root, $cds_tbl_HHAR, $mp_tbl_HHAR, $xfeat_tbl_HHHAR, $dfeat_tbl_HHHAR, $seq_info_HAR, $ofile_info_HHR, $opt_HHR, $FH_HR) = @_;

  my $cmd; # a command to be run by runCommand()
  my $do_matpept = opt_IsOn("--matpept", $opt_HHR);

  # should we skip the edirect calls?
  # Yes, if either --skipedirect or --infasta options have been used.
  my $do_skipedirect = (opt_Exists("--skipedirect", $opt_HHR) && opt_Get("--skipedirect", $opt_HHR)) ? 1 : 0;
  my $do_infasta     = (opt_Exists("--infasta", $opt_HHR)     && opt_Get("--infasta", $opt_HHR))     ? 1 : 0;
  my $do_skip = ($do_skipedirect || $do_infasta) ? 1 : 0;

  my $have_listfile = (defined $listfile) ? 1 : 0;
  if(defined $listfile) { 
    if(! -e $listfile) { DNAORG_FAIL("ERROR in $sub_name, $listfile does not exist", 1, $FH_HR); }
    if(! -s $listfile) { DNAORG_FAIL("ERROR in $sub_name, $listfile exists but is empty", 1, $FH_HR); }
  }

  # We create the .mat_peptide file first because we will die with an
  # error if mature peptide info exists and neither --matpept nor
  # --nomatpept was used (and we want to die as early as possible in the
  # script to save the user's time)
  #
  # 1) create the edirect .mat_peptide file, if necessary
  my $mp_file = $out_root . ".mat_peptide";

  #      if --nomatpept was   enabled we don't attempt to create a matpept file
  # else if --matpept was     enabled we validate that the resulting matpept file is not empty
  # else if --matpept was not enabled we validate that the resulting matpept file is     empty
  if((! $do_skip) && (! opt_Get("--nomatpept", $opt_HHR))) { 
    if($have_listfile) {
      $cmd = "cat $listfile | epost -db nuccore -format acc";
    }
    else { 
      $cmd = "esearch -db nuccore -query $ref_accn";
    }
    $cmd .= " | efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product > $mp_file";
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);
  
    if($do_matpept) { 
      if(! -s  $mp_file) { 
        DNAORG_FAIL("ERROR, in $sub_name, --matpept enabled but no mature peptide information exists.", 1, $FH_HR);
      }
      addClosedFileToOutputInfo($ofile_info_HHR, "mp", $mp_file, 0, "Mature peptide information obtained via edirect");
    }
    else { # ! $do_matpept
      if(-s $mp_file) { 
        # determine all the accessions in the mat_peptide file
        my $accn_list = `grep . $mp_file | awk \'{ print \$1 }\' | sort | uniq`; 
        # We allow this, the 'consopts' file enforces that build and annotate were run with the same options
        # related to mature peptide annotation, and that's good enough. We may see a situation where a 
        # non-reference accession has mat_peptide annotation, and so if we left in the following FAIL
        # statement we would fail because --nomatpept was not used, and the user would then have to 
        # rerun build with --nomatpept to get dnaorg_annotate.pl to run (or doctor the consopts file)
        # either way that's bad, so we let this go.
        # DNAORG_FAIL("ERROR, in $sub_name, --matpept not enabled but mature peptide information exists for the accessions printed below.\nUse --nomatpept to ignore it.\nAccessions with at least 1 mat_peptide annotation:\n$accn_list", 1, $FH_HR); 
        runCommand("rm $mp_file", opt_Get("-v", $opt_HHR), $FH_HR);
      }
      else { 
        # remove the empty file we just created
        runCommand("rm $mp_file", opt_Get("-v", $opt_HHR), $FH_HR);
      }
    }
  }

  # 2) create the edirect .fetched.ftable file
  # create the edirect ftable file
  my $ft_file  = $out_root . ".fetched.ftable";
  if(! $do_skip) { 
    if($have_listfile) { 
      $cmd = "cat $listfile | epost -db nuccore -format acc";
    }
    else { 
      $cmd = "esearch -db nuccore -query $ref_accn";
    }
    $cmd .= " | efetch -format ft > $ft_file";
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);
    addClosedFileToOutputInfo($ofile_info_HHR, "ft", $ft_file, 0, "Feature table obtained via edirect");
  }

  # 3) create the length file
  # create a file with total lengths of each accession
  my $len_file  = $out_root . ".length";
  if(! $do_skip) { 
    if($have_listfile) { 
      $cmd = "cat $listfile | epost -db nuccore -format acc";
    }
    else { 
      $cmd = "esearch -db nuccore -query $ref_accn";
    }
    $cmd .= " | efetch -format gpc | xtract -insd INSDSeq_length | grep . | sort > $len_file";
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);
    addClosedFileToOutputInfo($ofile_info_HHR, "len", $len_file, 0, "Sequence length file");
  }
  if(! -s $len_file) { 
    if($do_skip) { 
      DNAORG_FAIL("ERROR, no length information in file $len_file, are you sure you should be using --skipedirect?", 1, $FH_HR); 
    }
    else { 
      DNAORG_FAIL("ERROR, no length information obtained using edirect.", 1, $FH_HR); 
    }
  }

  # 4) parse the edirect .mat_peptide file, if necessary
  if($do_matpept) {  
    edirectFtableOrMatPept2SingleFeatureTableInfo($mp_file, 1, "mat_peptide", $mp_tbl_HHAR, $FH_HR); # 1: it is a mat_peptide file
    if (! exists ($mp_tbl_HHAR->{$ref_accn})) { 
      DNAORG_FAIL("ERROR in $sub_name, --matpept enabled, but no mature peptide information stored for reference accession", 1, $FH_HR); 
    }
  }

  # 5) parse the edirect .ftable file
  edirectFtableOrMatPept2SingleFeatureTableInfo($ft_file, 0, "CDS", $cds_tbl_HHAR, $FH_HR); # 0: it's not a mat_peptide file
  if(! exists ($cds_tbl_HHAR->{$ref_accn})) { 
    DNAORG_FAIL("ERROR in $sub_name, no CDS information stored for reference accession", 1, $FH_HR); 
  }

  foreach my $xfeat (keys (%{$xfeat_tbl_HHHAR})) { 
    edirectFtableOrMatPept2SingleFeatureTableInfo($ft_file, 0, $xfeat, \%{$xfeat_tbl_HHHAR->{$xfeat}}, $FH_HR); # 0: it's not a mat_peptide file
    if(! exists ($xfeat_tbl_HHHAR->{$xfeat}{$ref_accn})) { 
      DNAORG_FAIL("ERROR in $sub_name, no $xfeat information stored for reference accession", 1, $FH_HR); 
    }
  }    

  foreach my $dfeat (keys (%{$dfeat_tbl_HHHAR})) { 
    edirectFtableOrMatPept2SingleFeatureTableInfo($ft_file, 0, $dfeat, \%{$dfeat_tbl_HHHAR->{$dfeat}}, $FH_HR); # 0: it's not a mat_peptide file
    if(! exists ($dfeat_tbl_HHHAR->{$dfeat}{$ref_accn})) { 
      DNAORG_FAIL("ERROR in $sub_name, no $dfeat information stored for reference accession", 1, $FH_HR); 
    }
  }    

  # 6) parse the length file, and store accession lengths in $seq_info_HAR
  #    if --infasta was enabled, caller will likely have passed in a special $seq_info_HAR
  #    here (something like $ref_seq_info_HAR) which will only hold information on 
  #    the reference sequence, which will be read from the dnaorg_build length file.
  my %accn_len_H = ();
  parseLengthFile($len_file, \%accn_len_H, $FH_HR);
  my $nseq = 0;
  if(! $do_infasta) {  
    # we will have already stored all sequences in $seq_info_HAR->{"accn_name"}
    $nseq = scalar(@{$seq_info_HAR->{"accn_name"}});
    if($nseq == 0) { 
      DNAORG_FAIL("ERROR in $sub_name, no accessions in seq_info_HAR", 1, $FH_HR); 
    }
    for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
      if(! exists $accn_len_H{$accn_name}) { 
        DNAORG_FAIL("ERROR in $sub_name, problem fetching length of accession $accn_name", 1, $FH_HR); 
      }
      $seq_info_HAR->{"accn_len"}[$seq_idx] = $accn_len_H{$accn_name};
    }
  }
  else { # $do_infasta is true, $seq_info_HAR->{"accn_name"} is empty, but we should have
         # read exactly 1 length in parseLengthFile
    $nseq = scalar(keys %accn_len_H);
    if($nseq != 1) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, not exactly 1 accession read from length file $len_file (%d lengths read)", scalar(keys %accn_len_H)), 1, $FH_HR); 
    }
    my ($ref_accn) = (keys %accn_len_H);
    $seq_info_HAR->{"accn_name"}[0] = $ref_accn;
    $seq_info_HAR->{"accn_len"}[0]  = $accn_len_H{$ref_accn};
    $seq_info_HAR->{"seq_len"}[0]   = (opt_Exists("-c", $opt_HHR) && opt_Get("-c", $opt_HHR)) ? 2 * $accn_len_H{$ref_accn} : $accn_len_H{$ref_accn};
  }

  return 0;
}

#################################################################
# Subroutine:  wrapperFetchAllSequencesAndProcessReferenceSequence()
# Incept:      EPN, Fri Feb 26 09:27:39 2016
#
# Purpose:     A large block of code that is called once by each
#              dnaorg_build.pl and dnaorg_annotate.pl to fetch
#              either the reference sequence (if dnaorg_build.pl is caller) 
#              or all sequences (if dnaorg_annotate.pl is caller), 
#              and process information about the reference sequence.
#              This function does the following:
# 
#              1) fetches the sequences listed in @{$seq_info_HAR->{"accn_name"}} into a 
#                 fasta file and indexes that fasta file, the reference sequence is $seq_info_HAR->{"accn_name"}[0].
#              2) determines information for each feature (strand, length, coordinates, product) in the reference sequence
#              3) determines type of each reference sequence feature ('cds-mp', 'cds-notmp', or 'mp')
#              4) fetches the reference sequence feature and populates information on the models and features
#              5) look for special cases, where we want to append the 3 nt 3' of the final mature peptide in a cds-mp feature type
#              6) add information on the overlaps and adjacencies
#
#              Creates the following output files and stores
#              information on them in the ofile* data structures
#              by calling the addClosedFileToOutputInfo() function:
#              - $out_root . ".ft.idfetch.in": input file for esl-fetch-cds.pl
#              - $out_root . ".fg.fa":         sequence file with reference genome 
#              - $out_root . ".ref.all.stk":   Stockholm alignment file with reference features
#
#              The $out_root. ".ref.all.stk" will be used by the caller differently 
#              depending on if the caller is dnaorg_build.pl or dnaorg_annotate.pl:
#
#              dnaorg_build.pl:    .stk file will be used to create homology models (CMs)
#              
#              dnaorg_annotate.pl: .stk file will be used to verify the CM model we are about
#                                  to use (from a prior dnaorg_build.pl run) was built from the
#                                  current reference sequences. That is, that the homology models
#                                  are still current with respect to the GenBank reference annotation.
#                                  This check is done using a checksum value stored in the CM file
#                                  derived from the .stk file used to create the CMs.
#              
# Arguments: 
#   $execs_HR:              REF to hash with executables, the key "esl_fetch_cds"
#                           must be defined and the value must be a valid path to an 
#                           esl_fetch_cds Perl script, PRE-FILLED
#   $sqfile_R:              REF to a Bio::Easel::SqFile object, CREATED HERE
#   $out_root:              string that is the 'root' for naming output files
#   $build_root:            string that is the 'root' for files possibly already created with dnaorg_build.pl (may be same as $out_root)
#   $in_ref_accn:           reference accession, undef unless --infasta enabled
#   $in_ref_accnlen:        length of the reference accession, undef unless --infasta enabled
#   $in_ref_seqlen:         length of the reference sequence (possibly double $in_ref_accnlen), undef unless --infasta enabled
#   $infasta_file:          name of input fasta file, should be undef unless --infasta enabled
#   $cds_tbl_HHAR:          REF to CDS hash of hash of arrays, PRE-FILLED
#   $mp_tbl_HHAR:           REF to mature peptide hash of hash of arrays, can be undef, else PRE-FILLED
#   $xfeat_tbl_HHHAR:       REF to eXtra feature hash of hash of hash of arrays, can be undef, else PRE-FILLED
#   $dfeat_tbl_HHHAR:       REF to Duplicate feature hash of hash of hash of arrays, can be undef, else PRE-FILLED
#   $cds2pmatpept_AAR:      REF to 2D array, 1st dim: cds index (-1, off-by-one), 
#                           2nd dim: value array of primary matpept indices that comprise this CDS, 
#                           may be undef if $mp_tbl_HHAR is undef
#   $cds2amatpept_AAR:      1st dim: cds index (-1, off-by-one), 
#                           2nd dim: value array of all matpept indices that comprise this CDS, 
#                           OR undefined if all features are CDS and there are no mature peptides; 
#                           PRE-FILLED
#   $mdl_info_HAR:          REF to hash of arrays with information on the models, FILLED HERE
#   $ftr_info_HAR:          REF to hash of arrays with information on the features, FILLED HERE
#   $seq_info_HAR:          REF to 2D hash with sequence information, ADDED TO HERE 
#   $opt_HHR:               REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#   $ofile_info_HHR:        REF to 2D hash with output info, ADDED TO HERE
#
# Returns:     void
#
# Dies:        if $listfile is defined but it does not exist or is empty
#              if any of the edirect commands exit in error
#     
################################################################# 
sub wrapperFetchAllSequencesAndProcessReferenceSequence { 
  my $sub_name = "wrapperFetchAllSequencesAndProcessReferenceSequence()";
  my $nargs_expected = 19;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name, entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $sqfile_R, $out_root, $build_root, $in_ref_accn, $in_ref_accnlen, $in_ref_seqlen, $infasta_file, $cds_tbl_HHAR, $mp_tbl_HHAR, $xfeat_tbl_HHHAR, $dfeat_tbl_HHHAR, $cds2pmatpept_AAR, $cds2amatpept_AAR, $mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;
  
  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience

  # contract check
  if((defined $mp_tbl_HHAR) && (! defined $cds2pmatpept_AAR)) { 
    DNAORG_FAIL("ERROR in $sub_name, contract violation, mp_tbl_HHAR variable is defined but cds2pmatpept_AAR variable is not", 1, $FH_HR);
  }
  if((! defined $mp_tbl_HHAR) && (defined $cds2pmatpept_AAR)) { 
    DNAORG_FAIL("ERROR in $sub_name, contract violation, mp_tbl_HHAR variable is not defined but cds2pmatpept_AAR variable is", 1, $FH_HR);
  }

  # determine which one of the three possible modes we are in:
  # (1) default:     fetching sequences using esl-fetch-cds.pl (which calls idfetch)
  # (2) --skipfetch: not fetching, instead using sequences fetched in a previous run 
  # (3) --infasta:   not fetching, instead using sequences in $infasta_file

  my $do_skipfetch = (opt_Exists("--skipfetch", $opt_HHR) && opt_Get("--skipfetch", $opt_HHR)) ? 1 : 0;
  my $do_infasta   = (opt_Exists("--infasta",   $opt_HHR) && opt_Get("--infasta",   $opt_HHR)) ? 1 : 0;

  # make sure that at most one of --skipfetch and --infasta are enabled
  if($do_skipfetch && $do_infasta) { 
    DNAORG_FAIL("ERROR in $sub_name, contract violation, both --skipfetch and --infasta options are enabled, only one should be", 1, $FH_HR);
  }

  # if $do_infasta, make sure $infasta_file is defined, if (! $do_infasta), make sure that it is not defined
  if($do_infasta) { 
    if(! defined $infasta_file) { 
      DNAORG_FAIL("ERROR in $sub_name, contract violation, --infasta option enabled but <infasta_file> is not defined", 1, $FH_HR);
    }
  }
  elsif(defined $infasta_file) { 
    DNAORG_FAIL("ERROR in $sub_name, contract violation, --infasta option NOT enabled but <infasta_file> is defined", 1, $FH_HR);
  }

  # $in_ref_accn, $in_ref_accnlen and $in_ref_seqlen must be defined if $do_infasta, and must not be defined if ! $do_infasta
  if($do_infasta && ((! defined $in_ref_accn) || (! defined $in_ref_accnlen) || (! defined $in_ref_seqlen))) { 
    DNAORG_FAIL("ERROR in $sub_name, contract violation, --infasta option enabled but in_ref_accn and/or in_ref_accnlen and/or in_ref_seqlen variables are not defined", 1, $FH_HR);
  }
  elsif((! $do_infasta) && ((defined $in_ref_accn || defined $in_ref_accnlen || defined $in_ref_seqlen))) { 
    DNAORG_FAIL("ERROR in $sub_name, contract violation, --infasta option NOT enabled but in_ref_accn and/or in_ref_accnlen and/or in_ref_seqlen variables are defined", 1, $FH_HR);
  }
    
  # 1) fetch the sequences into a fasta file and index that fasta file
  my $nseq = scalar(@{$seq_info_HAR->{"accn_name"}});
  my $have_only_ref = ($nseq == 1) ? 1 : 0; # '1' if we only have the reference
  my $ref_accn      = ($do_infasta) ? $in_ref_accn    : $seq_info_HAR->{"accn_name"}[0]; # the reference accession
  my $ref_accnlen   = ($do_infasta) ? $in_ref_accnlen : $seq_info_HAR->{"accn_len"}[0]; # length of the accession
  my $fetch_file    = undef;
  my $fasta_file    = undef;

  if($do_infasta) { 
    $fasta_file = $infasta_file;
    # this may not actually be the input fasta file provided on the dnaorg_annotate.pl command line,
    # if -c was enabled then the process_input_fasta_file() function will have created a new file
    # that includes the duplicate of all the sequences in that user provided file
  }
  else { # $do_infasta is false 
    $fetch_file = sprintf($out_root . "%s.fg.idfetch.in", ($have_only_ref) ? ".ref" : "");  # the esl_fetch_cds.pl input file we are about to create
    $fasta_file = sprintf($out_root . "%s.fg.fa",         ($have_only_ref) ? ".ref" : "");  # the fasta file we are about to create
  }

  @{$seq_info_HAR->{"seq_name"}} = ();
  if((! $do_skipfetch) && (! $do_infasta)) { 
    # fetch the sequences into a fasta file
    fetchSequencesUsingEslFetchCds($execs_HR->{"esl_fetch_cds"}, $fetch_file, $fasta_file, opt_Get("-c", $opt_HHR), $seq_info_HAR, $FH_HR);
    addClosedFileToOutputInfo($ofile_info_HHR, "fetch", $fetch_file, 0, "Input file for esl-fetch-cds.pl");
    addClosedFileToOutputInfo($ofile_info_HHR, "fasta", $fasta_file, 0, "Sequence file (fetched with esl-fetch-cds.pl)");
  }
  else { # we should already have the fasta file that we need 
    validateFileExistsAndIsNonEmpty($fasta_file, $sub_name, $FH_HR);
    if($do_skipfetch) { 
      addClosedFileToOutputInfo($ofile_info_HHR, "fasta", $fasta_file, 0, "Sequence file (fetched on previous run (--skipfetch))");
    }
    else { # $do_infasta is true
      if(opt_Get("-c", $opt_HHR)) { 
        addClosedFileToOutputInfo($ofile_info_HHR, "fasta", $fasta_file, 0, "Sequence file (same sequences as supplied with --infasta, but all seqs are duplicated)");
      }
      else { 
        addClosedFileToOutputInfo($ofile_info_HHR, "fasta", $fasta_file, 0, "Sequence file (supplied with --infasta), NOT CREATED BY THE SCRIPT");
      }
    }
  }  

  # open the sequence file using Bio-Easel
  # remove the .ssi file first, if it exists
  my $ssi_file = $fasta_file . ".ssi";
  if(-e $ssi_file) { 
    runCommand("rm $ssi_file", opt_Get("-v", $opt_HHR), $FH_HR);
  }
  $$sqfile_R = Bio::Easel::SqFile->new({ fileLocation => $fasta_file }); # the sequence file object

  if($do_skipfetch || $do_infasta) { 
    # --skipfetch or --infasta enabled, so we never called 
    # fetchSequencesUsingEslFetchCds() above.
    # We need to fill $seq_info_HAR->{"seq_name"} and $seq_info_HAR->{"seq_len"}
    # get index hash for @{$seq_info_HAR->{"seq_accn"}} array
    # this simplifies determining sequence index in @{%seq_info_HAR->{}}
    # arrays for a given accession name.
    my %accn_name_idx_H = (); # key: $accn_name, value: idx of $accn_name in @{$seq_info_HAR->{"accn_name"}}
    getIndexHashForArray($seq_info_HAR->{"accn_name"}, \%accn_name_idx_H, $FH_HR);
    for(my $sqfile_seq_idx = 0; $sqfile_seq_idx < $nseq; $sqfile_seq_idx++) { 
      my ($seq_name, $seq_len) = $$sqfile_R->fetch_seq_name_and_length_given_ssi_number($sqfile_seq_idx);
      my $accn_name = undef;
      if($do_infasta && (! opt_Get("-c", $opt_HHR))) { 
        # if --infasta and ! -c, then we are fetching directly from the --infasta file, else 
        # we created the file we are fetching from, and so we need to determine the
        # accession from the sequence name using accn_name_from_seq_name().
        $accn_name = $seq_name;
      }
      else { 
        $accn_name = accn_name_from_seq_name($seq_name, $FH_HR);
      }
      if(! exists $accn_name_idx_H{$accn_name}) { 
        DNAORG_FAIL("ERROR in $sub_name, accession $accn_name derived from sqfile seq name: $seq_name does not exist in accn_name_idx_H", 1, $ofile_info_HHR->{"FH"});
      }
      $seq_info_HAR->{"seq_name"}[$accn_name_idx_H{$accn_name}] = $seq_name;
      $seq_info_HAR->{"seq_len"}[$accn_name_idx_H{$accn_name}]  = $seq_len;
    }
  }

  # 2) determine reference information for each feature (strand, length, coordinates, product)
  getReferenceFeatureInfo($cds_tbl_HHAR, $mp_tbl_HHAR, $xfeat_tbl_HHHAR, $dfeat_tbl_HHHAR, $ftr_info_HAR, $ref_accn, $FH_HR); # $mp_tbl_HHAR may be undefined and that's okay
  my @reqd_keys_A = ("ref_strand", "ref_len", "ref_coords", "out_product", "out_gene", "out_exception", "type_fname", "type_ftable");
  validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, \@reqd_keys_A, $FH_HR);

  # 3) determine type of each reference feature ('cds-mp', 'cds-notmp', 'mp', 'xfeat', or 'dfeat')
  my $ncds   = (defined $cds_tbl_HHAR) ? scalar(@{$cds_tbl_HHAR->{$ref_accn}{"coords"}}) : 0; # number of CDS features
  my $nmp    = (defined $mp_tbl_HHAR)  ? scalar(@{$mp_tbl_HHAR->{$ref_accn}{"coords"}})  : 0; # number of mature peptides
  my $nxfeat = (defined $xfeat_tbl_HHHAR) ? getNumExtraOrDuplicateFeatures($xfeat_tbl_HHHAR, $FH_HR) : 0;
  my $ndfeat = (defined $dfeat_tbl_HHHAR) ? getNumExtraOrDuplicateFeatures($dfeat_tbl_HHHAR, $FH_HR) : 0;
  determineFeatureTypes($ncds, $nmp, $nxfeat, $ndfeat, $cds2pmatpept_AAR, $cds2amatpept_AAR, $ftr_info_HAR, $FH_HR); # $cds2pmatpept_AAR may be undef and that's okay

  # 4) fetch the reference feature sequences and populate information on the models and features
  #    we won't actually fetch the reference sequence if --infasta is used (which is why $ref_seqname is undef if --infasta)
  my $ref_seqname  = ($do_infasta) ? undef          : $seq_info_HAR->{"seq_name"}[0]; # the reference sequence name the fetched sequence file $fasta_file
  my $all_stk_file = $out_root . ".ref.all.stk";     # name of output alignment file we are about to create, each single reference feature 
                                                     # sequence is a separate 'Stockholm (stk) alignment', and this single file contains all such 
                                                     # separate alignments, one per feature
  fetchReferenceFeatureSequences($execs_HR, $$sqfile_R, $ref_seqname, $ref_accnlen, $out_root, $build_root, $mdl_info_HAR, $ftr_info_HAR, $all_stk_file, $opt_HHR, $FH_HR); 
  if(! $do_infasta) { 
    addClosedFileToOutputInfo($ofile_info_HHR, "refstk", $all_stk_file, 0, "Stockholm alignment file with reference features");
  }

  # 5) look for special cases, where we want to append the 3 nt 3' of the final mature peptide in a cds-mp feature type
  annotateAppendFeatures($ftr_info_HAR, $mdl_info_HAR, $FH_HR);

  # 6) add information on the overlaps and adjacencies
  my $ref_seqlen = ($do_infasta) ? $in_ref_seqlen : $seq_info_HAR->{"seq_len"}[0]; # length of the sequence
  annotateOverlapsAndAdjacencies($ref_accnlen, $ref_seqlen, $mdl_info_HAR, $opt_HHR, $FH_HR);

  # 6) determine the source "dupsource" of any duplicate feature annotations, for
  #    all other non 'duplicate' features the "source_idx" will be -1
  determineSourcesOfDuplicateFeatures($ftr_info_HAR, $FH_HR);
  
  return 0;
}

#################################################################
#################################################################
#
# Subroutines related to feature tables output from edirect:
#   edirectFtableOrMatPept2SingleFeatureTableInfo()
#   getSingleFeatureTableInfo()
#   helperBreakdownFac()
#
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
#   $edirect_file:  name of edirect output file to parse
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
  my $dummy_column    = "DuMmY"; # name of 'dummy' column that we don't actually print, for use with --matpept only
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
      stripVersion(\$accn);
      if(%{$tbl_HHAR} && exists $tbl_HHAR->{$accn}) { 
        DNAORG_FAIL("ERROR in $sub_name, read data twice in $sub_name for accession $accn", 1, $FH_HR);
      }
      foreach $fac (@{$fac_HAR->{$faccn}}) { # foreach 'fac', accession + set of coords

        ($faccn2, $coords, $sort_coord, $strand) = helperBreakdownFac($fac, $fac_sep, $FH_HR);
        if($faccn ne $faccn2) { DNAORG_FAIL("ERROR in $sub_name, inconsistent fac value: $faccn ne $faccn2", 1, $FH_HR); }
        
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
                    if($qval =~ m/\Q$qval_sep/) { DNAORG_FAIL("ERROR in $sub_name, qualifier_name $qval has the string $qval_sep in it", 1, $FH_HR); }
                    $save_str .= $qval_sep . $qval; # not first value, concatenate onto previous values
                  }
                }
              }
              # old behavior:
              ## if there's no value for this qualifier, put '<empty>'
              ##if($save_str eq "") { $save_str = "<empty>"; }

              # do not save values that are '-', as far as I can tell these are just empty placeholders,
              # as an example, see the first CDS in the feature table returned by this query (as of 02/13/18):
              # esearch -db nuccore -query NC_031324 | efetch -format ft 
              # It has a qualifier value of '-' for the 'gene' qualifier, which only occurs because the
              # GenBank flat file does not list a gene name (gene qualifier) for the first gene (all 
              # other genes have gene qualifiers and so their corresponding CDS features do not have 
              # 'gene' qualifiers.
              if($save_str eq "-") { $save_str = ""; } 
              push(@{$tbl_HHAR->{$accn}{$column}}, $save_str);
            }
          } 
        }
      }
    }
  }
  
  return;
}

#################################################################
# Subroutine: helperBreakdownFac()
# Incept      EPN, Thu Feb 11 14:14:23 2016
#
# Purpose:    Breakdown a 'fac' string into it's parts and 
#             create a string in NCBI coordinate format from it.
#             A 'fac' string has accessions and coordinates in it.
#           
#             A 'helper' function called by 'getSingleFeatureTableInfo()'.
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
#################################################################
#
# Subroutines for parsing different file types:
#   parseMatPeptSpecFile()
#   parseLengthFile()
#   parseEdirectFtableFile()
#   parseEdirectMatPeptideFile()
#   parseListFile()
#   parseSpecStartFile()
#   parseConsOptsFile()
#   parseNonConsOptsFile()
#
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
  my $ncds_read_primary = 0;  # number of 'primary' CDS lines read
  my $ncds_read_all     = 0;  # number of 'all' CDS lines read
  my $cds_idx2store     = 0;  # index of CDS to store
  my @cds_idx_read_A    = (); # $cds_idx_read_A[$i] = 1 if we read info for CDS $i+1
  my $max_cds_idx2store = 0;  # maximum $cds_idx2store seen

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

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
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) != 3) { 
        DNAORG_FAIL("ERROR in $sub_name, unable to parse matpept input file line: $line\nline should have three tokens separated by whitespace", 1, $FH_HR); 
      }
      my ($cds_idx, $primary_or_all, $matpept_str) = ($el_A[0], $el_A[1], $el_A[2]);
      if($cds_idx !~ m/^\d+$/ || $cds_idx eq "0") { 
        DNAORG_FAIL("ERROR in $sub_name, parsing $infile, first token should be a positive integer, read $cds_idx on line $line", 1, $FH_HR);
      }
      $primary_or_all =~ tr/A-Z/a-z/;
      if($primary_or_all ne "all" && $primary_or_all ne "primary") { 
        DNAORG_FAIL("ERROR in $sub_name, parsing $infile, second token of each non-comment line should be 'all' or 'primary', found $primary_or_all in line $line", 1, $FH_HR);
      }
      my $cds_idx2store = $cds_idx - 1; # we verified this was a positive integer above
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

  open(LEN, $lenfile) || fileOpenFailure($lenfile, $sub_name, $!, "reading", $FH_HR);

  while(my $line = <LEN>) { 
    # example line
    #HM448898.1	2751

    chomp $line;
    if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
    my ($orig_accn, $length) = split(/\s+/, $line);
    if($length !~ m/^\d+$/) { 
      DNAORG_FAIL("ERROR in $sub_name, couldn't parse length file line: $line", 1, $FH_HR); 
    } 

    my $accn = $orig_accn;
    stripVersion(\$accn);
    if(exists $len_HR->{$accn}) { 
      DNAORG_FAIL("ERROR in $sub_name, accession $orig_accn exists more than once in $lenfile (possibly with different versions)", 1, $FH_HR); 
    } 
      
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
#            Caller will commonly call getSingleFeatureTableInfo() after calling this
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

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

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
  my %column_HH       = ();      # existence 2D hash, used to aid construction of %column_HA
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
    if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
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
#            edirectFtableOrMatPept2SingleFeatureTableInfo().
#
#            Caller will commonly call getSingleFeatureTableInfo() after calling this
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

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  # example lines of a .mat_peptide file
  # NC_009942.1	97..465	anchored capsid protein C	
  # NC_009942.1	97..411	capsid protein C	
  # NC_009942.1	466..966	membrane glycoprotein precursor M	
  # NC_009942.1	join(2470..3552,3552..3680)	nonstructural protein NS1'

  my $faccn     = undef;   # full accession, e.g. ref|NC_000883.2|
  my $fac       = undef;   # full accession and coords, accession and each segment coords separated by $fac_sep
  my $coords    = undef;   # coordinates
  my $faccn2    = undef;   # another full accession
  my %column_HH = ();      # existence 2D hash, used to aid construction of %column_HA
                           # key 1: feature, e.g. "CDS"
                           # key 2: qualifier
                           # value: '1' (always)

  my $fac_sep  = $sep_HR->{"fac"};
  my $qnqv_sep = $sep_HR->{"qnqv"};
  my $qval_sep = $sep_HR->{"qval"};

  my $feature = "mat_peptide";
  while(my $line = <IN>) { 
    if($line =~ m/\w/) { 
# example line:
# NC_001477.1	95..436	anchored capsid protein C
      if($line =~ /(\S+)\s+(\S+)\s*(.*)$/) { 
        my ($acc, $coords, $product) = ($1, $2, $3); 
        # with .mat_peptide files, the accession is not in the 'full accession' ($faccn, e.g. ref|NC_001477.1|) 
        # that appears in the .ftable files, but we still want to use similar downstream code to 
        # handle both the parsed .mat_peptide and .ftable files, so we still fill faccn2accn_HR
        # for .mat_peptide files, even though the keys and values are identical

        $product =~ s/\s+$//; # remove trailing whitespace
        $faccn = $acc; # in this case, full accession and accession are the same
        # we use both $faccn and $acc so that this code matches the code in this subroutine matches
        # with the code in parseEdirectFtableFile() for consistency.

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
        # we need to do this so that we can store the coordinate information (in $fac)
        # for mature peptides that have no additional information besides coordinates
        # (.mat_peptide files often include only mature peptide coordinates and no
        # product information ($product below) 
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
# Subroutine : parseListFile()
# Incept:      EPN, Thu Feb 18 13:05:30 2016
#
# Purpose:     Given a file name remove the directory path.
#              For example: "foodir/foodir2/foo.stk" becomes "foo.stk".
#
# Arguments: 
#   $listfile: name of list file to parse
#   $do_accn:  '1' if lines are accessions, else '0'.
#              If '1', this changes behavior 
#              of the function in two ways: 
#              1) we should strip version from each accession
#              2) if any accessions are duplicated in the list
#                 file, we should print all duplicates and die
#   $line_AR:  REF to array to fill, each element will be a line 
#              of $listfile with newline removed. FILLED HERE
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd", can be undef
#
# Returns:     void, fills @{$line_AR}
#
# Dies:        if $listfile does not exist or is not readable, or
#              if $listfile exists as a directory.
#              if $do_accn is 1 and >=1 accession occurs more 
#              than once in $listfile
################################################################# 
sub parseListFile {
  my $nargs_expected = 4;
  my $sub_name = "parseListFile()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name, entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($infile, $do_accn, $line_AR, $FH_HR) = @_;

  # a hash that keeps track of counts of accessions, only used if $do_accn = 1
  my %accn_ct_H    = (); # key: an accession, value: number of times accession occurs (possibly with different versions)
  my @accn_order_A = (); # the accessions read in the input file, in order

  if(-d $infile) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to read list file $infile, but a directory of the same name exists.", 1, $FH_HR);
  }

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  while(my $line = <IN>) { 
    if($line =~ m/\w/) {  # skip blank lines
      chomp $line;
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
      if($do_accn) { 
        my $accn = $line;
        stripVersion(\$accn); # remove version from $accn
        $accn_ct_H{$accn}++; 
        if($accn_ct_H{$accn} == 1) { 
          push(@accn_order_A, $accn);
        }
      }
      push(@{$line_AR}, $line);
    }
  }

  # check if we need to exit because we read >= 1 accessions 
  # more than once
  if($do_accn) { 
    my $errmsg = ""; # we may fill this below
    foreach my $accn (@accn_order_A) { 
      if($accn_ct_H{$accn} > 1) { 
        $errmsg .= "$accn\n";
      }
    }
    if($errmsg ne "") { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, the following accessions occur on multiple lines, possibly with different versions:\n%s", $errmsg), 1, $FH_HR);
    }
  }
    
  close(IN); 

  return;
}

#################################################################
# Subroutine:  parseSpecStartFile()
# Incept:      EPN, Thu Feb 18 15:45:49 2016
#
# Purpose:     Parse the input file that defines non-standard start 
#              codons for >= 1 CDS.
#
# Arguments: 
#   $infile:        file to parse
#   $specstart_AAR: ref to array of arrays to fill here, allowed start codons for each CDS
#                   if array doesn't exist for a CDS, ATG is only allowed start
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd", can be undef
#
# Returns:     void, fills @{$specstart_AAR}
#
# Dies:        if $infile does not exist or is not readable,
#              or we have some problem parsing it.
################################################################# 
sub parseSpecStartFile { 
  my $sub_name = "parseSpecStartFile";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name, entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $specstart_AAR, $FH_HR) = @_;

  my $ncds_read      = 0;
  my $cds_idx2store  = 0;
  my @cds_idx_read_A = (); # $cds_idx_read_A[$i] = 1 if we read info for CDS $i+1
  my $max_cds_idx2store = 0;

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  while(my $line = <IN>) { 
    if($line !~ m/^\#/) { 
      ## example input file:
      ## This file explains alternative start codons that are expected 
      ## West Nile lineage 1 CDS #3 (WARF4)
      ##
      ## Format of lines in this file:
      ## <CDS-idx> <alternate start codon 1>:<alternate start codon 2>:<alternate start codon n>
      #3 GGC
      #####################
      # NOTE: in the input file CDS and matpept indices are in the coordinate range 1..N, but we 
      # store them in the coordinate range 0..N-1
      chomp $line;
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) != 2) { 
        DNAORG_FAIL("ERROR in $sub_name, unable to parse specstart input file line: $line", 1, $FH_HR); 
      }
      my ($cds_idx, $codon_str) = ($el_A[0], $el_A[1]);
      my $cds_idx2store = $cds_idx - 1;
      if($cds_idx2store < 0) { 
        DNAORG_FAIL("ERROR in $sub_name, read CDS idx that is 0 or less ($cds_idx) in matpept input file in line $line", 1, $FH_HR); 
      }
      $cds_idx_read_A[$cds_idx2store] = 1;
      if($cds_idx2store > $max_cds_idx2store) { 
        $max_cds_idx2store = $cds_idx2store; 
      }
      my @codon_A = split(":", $codon_str);
      @{$specstart_AAR->[$cds_idx2store]} = ();
      foreach my $codon (@codon_A) { 
        $codon =~ tr/a-z/A-Z/;
        $codon =~ tr/U/T/;
        push(@{$specstart_AAR->[$cds_idx2store]}, $codon);
      }
      $ncds_read++;
    }
  }
  close(IN);

  # Two sanity checks:
  # 1: we should have stored all and primary info for any CDS $i for which $cds_idx_read_A[$i-1] is 1, and nothing else
  for(my $i = 0; $i <= $max_cds_idx2store; $i++) { 
    if($cds_idx_read_A[$i]) { 
     if((! defined $specstart_AAR->[$i]) && (! exists $specstart_AAR->[$i])) { 
       DNAORG_FAIL(sprintf("ERROR in $sub_name, did not properly read info for cds %d in $infile\n", $i+1), 1, $FH_HR);
     }
    }
    else { # we didn't read this one
      if(defined $specstart_AAR->[$i] || exists $specstart_AAR->[$i]) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, improperly read non-existent info for cds %d in $infile\n", $i+1), 1, $FH_HR);
      }
    }
  }
  # 2: we should have at least read at least one CDS info
  if($ncds_read == 0) { 
    DNAORG_FAIL("ERROR in $sub_name, no CDS start codon specifications read in matpept input file $infile", 1, $FH_HR); 
  }

  return;
}

#################################################################
# Subroutine: parseConsOptsFile()
# Incept:     EPN, Mon Feb 12 12:51:35 2018
#
# Purpose:   Parse the special .consopts file created by 
#            dnaorg_build.pl, which lists options that need
#            to be consistently used between dnaorg_build.pl and
#            dnaorg_annotate.pl.
#
# Arguments:
#  $consopts_file:       name of the dnaorg_build consopts file
#  $consopts_used_HR:    REF to hash of options in the consopts file
#                        key is option (e.g. --matpept), value is
#                        option argument, "" for none FILLED HERE 
#  $consopts_notused_HR: REF to options not used in the consopts file
#                        that could have been used, key is option,
#                        value is always 1.
#  $consmd5_HR:          REF to hash of MD5 values for files in the 
#                        consopts file, key is option (e.g. --matpept)
#                        that was used in dnaorg_build.pl, value is
#                        MD5 checksum value for the file associated 
#                        with the option, "" if option has no file associated
#                        with it. 
#                        FILLED HERE 
#  $FH_HR:               REF to hash of file handles
# 
# Returns:  void
# 
# Dies: If $consopts_file doesn't exist, or we can't parse it
#       because it includes an option that we don't expect
#       or is in an invalid format.
#
#################################################################
sub parseConsOptsFile { 
  my $sub_name = "parseConsOptsFile";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($consopts_file, $consopts_used_HR, $consopts_notused_HR, $consmd5_HR, $FH_HR) = @_;

  # read the consopts file
  if(! -e $consopts_file) { 
    DNAORG_FAIL("ERROR in $sub_name, consopts file $consopts_file does not exist.\nThis file should have been created by dnaorg_build.pl.\nYou probably need to rerun dnaorg_build.pl if it was run before May 31, 2016.", 1, $FH_HR);
  }
  open(IN, $consopts_file) || fileOpenFailure($consopts_file, $sub_name, $!, "reading", $FH_HR);
  my $line_ct = 0;

  # initialize the hashes
  %{$consopts_used_HR}    = (); 
  %{$consopts_notused_HR} = (); 
  %{$consmd5_HR}          = ();

  while(my $line = <IN>) { 
    chomp $line;
    $line_ct++;
    if(($line eq "none") && ($line_ct == 1)) { 
      ; # this is fine, none of the options that need to be consistent were set by dnaorg_build.pl
    }
    elsif($line =~ /^\-c$/) { 
      $consopts_used_HR->{"-c"} = "";
      $consmd5_HR->{"-c"}  = "";
    }
    elsif($line =~ /^\-\-nomatpept$/) { 
      $consopts_used_HR->{"--nomatpept"} = "";
      $consmd5_HR->{"--nomatpept"}  = "";
    }
    elsif($line =~ /^\-\-matpept\s+(\S+)\s+(\S+)$/) { # first string is file name, second is md5 checksum (obtained with 'md5sum' executable)
      ($consopts_used_HR->{"--matpept"}, $consmd5_HR->{"--matpept"}) = ($1, $2);
    }
    elsif($line =~ /^\-\-xfeat\s+(\S+)$/) { # first string is file name, second is argument
      $consopts_used_HR->{"--xfeat"} = $1;
      $consmd5_HR->{"--xfeat"}  = "";
    }
    elsif($line =~ /^\-\-dfeat\s+(\S+)$/) { # first string is file name, second is argument
      $consopts_used_HR->{"--dfeat"} = $1;
      $consmd5_HR->{"--dfeat"}  = "";
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse line from consopts file $consopts_file:\n$line\n", 1, $FH_HR);
    }
  }
  close(IN);

  # fill %{$consopts_notused_HR}
  foreach my $opt ("-c", "--nomatpept", "--matpept", "--xfeat", "--dfeat") { 
    if(! exists $consopts_used_HR->{$opt}) { 
      $consopts_notused_HR->{$opt} = 1;
    }
  }

  # if we get here, all options were recognized
  return;
}


#################################################################
# Subroutine: parseNonConsOptsFile()
# Incept:     EPN, Tue Feb 13 05:29:04 2018
#
# Purpose:   Read a file that contains command line options 
#            that are not required to be consistent between
#            different dnaorg scripts, and die if any 
#            options that are required to be consistent
#            are included. Also die if any options included
#            in the string $auto_opts are included. 
#            Return a string of verified, compatible options. 
# 
#            Format of file is a single line of options.
#
# Arguments:
#  $opts_file:      file with options
#  $cmdline_opt:    command line option that triggered this 
#                   subroutine call (e.g. "--Aopts")
#  $auto_opts_str:  string of options, separated by commas, which will be
#                   automatically added and so cannot listed in the file
#  $failure_str:    failure message (e.g. "that option will automatically be set, as required 
#                   to be consistent with relevant dnaorg_build.pl command used previously")
#  $FH_HR:          REF to hash of file handles
# 
# Returns:  void
# 
# Dies: If any of the options in $opts_file are 
#       those that are required to be consistent between
#       dnaorg_scripts.pl.
#       
#################################################################
sub parseNonConsOptsFile { 
  my $sub_name = "parseNonConsOptsFile";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($opts_file, $cmdline_opt, $auto_opts_str, $failure_str, $FH_HR) = @_;

  my $opts_str = "";

  # read the opts file
  if(! -e $opts_file) { 
    DNAORG_FAIL("ERROR in $sub_name, options file $opts_file does not exist.", 1, $FH_HR);
  }
  open(IN, $opts_file) || fileOpenFailure($opts_file, $sub_name, $!, "reading", $FH_HR);
  my $line_ct = 0;

  while(my $line = <IN>) { 
    if($line =~ m/\w/) { 
      chomp $line;
      $line_ct++;
      $opts_str .= $line;
    }
  }
  if($line_ct != 1) { 
    DNAORG_FAIL("ERROR in $sub_name, file $opts_file is required to have a single line (multiple options separated by spaces), read $line_ct.", 1, $FH_HR);
  }        
  close(IN);

  # check for the hard-coded illegal options that must be 
  # consistent between dnaorg_annotate.pl and dnaorg_build.pl
  if($opts_str =~ m/\s*\-c/)           { DNAORG_FAIL("ERROR with $cmdline_opt, command-line options cannot include -c, $failure_str", 1, $FH_HR); }
  if($opts_str =~ m/\s*\-\-nomatpept/) { DNAORG_FAIL("ERROR with $cmdline_opt, command-line options cannot include --nomatpept, $failure_str", 1, $FH_HR); }
  if($opts_str =~ m/\s*\-\-matpept/)   { DNAORG_FAIL("ERROR with $cmdline_opt, command-line options cannot include --matpept, $failure_str", 1, $FH_HR); }
  if($opts_str =~ m/\s*\-\-xfeat/)     { DNAORG_FAIL("ERROR with $cmdline_opt, command-line options cannot include --xfeat, $failure_str", 1, $FH_HR); }
  if($opts_str =~ m/\s*\-\-dfeat/)     { DNAORG_FAIL("ERROR with $cmdline_opt, command-line options cannot include --dfeat, $failure_str", 1, $FH_HR); }

  # check for the options that are automatically set or not set, 
  # these are illegal too
  foreach my $auto_opt (split(",", $auto_opts_str)) { 
    my $esc_auto_opt = quotemeta($auto_opt);
    if($opts_str =~ m/\s*$esc_auto_opt/) { 
      DNAORG_FAIL("ERROR with $cmdline_opt, command-line options cannot include $auto_opt, it will be automatically set by the program", 1, $FH_HR); 
    }
  }

  return $opts_str;
}
#################################################################
#################################################################
#
# Subroutines related to parsing NCBI coordinate strings:
#   getStrandStats()
#   startsStopsStrandsFromCoordsLength()
#   startsStopsFromCoords()
#   getLengthsAndCoords()
#   lengthFromCoords()
#
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
      else { DNAORG_FAIL(sprintf("ERROR in $sub_name, unable to parse strand (%s) for feature %d for $accn\n", $tbl_HHAR->{$accn}{"strand"}[$i], $i+1), 1, $FH_HR); }
      $strand_str .= $tbl_HHAR->{$accn}{"strand"}[$i];
    }
  }

  return ($nfeatures, $npos, $nneg, $nunc, $nbth, $strand_str);
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
    elsif($el =~ m/^(\d+)$/) { # a single nucleotide
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $1);
      if(defined $strands_AR) { push(@{$strands_AR}, $cur_strand); }
      $$nexons_R++;
    }
    else { 
      DNAORG_FAIL("ERROR unable to parse $orig_coords in $sub_name", 1, $FH_HR); 
    }
  }

  # check if we have a spanning exon (that spans stop..start) and if we do
  # and (! $do_circular) then die, because that shouldn't happen.
  my $have_spanning_exon = checkForSpanningSequenceSegments($starts_AR, $stops_AR, $nexons_R, 0, $strand, $totlen); # 1 says: do correct the spanning exon
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
    elsif($el =~ m/^(\d+)$/) { # a single nucleotide
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $1);
      $$nexons_R++;
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse coordinates $orig_coords", 1, $FH_HR); 
    }
  }

  return;
}

#################################################################
# Subroutine: getLengthsAndCoords()
# Incept:     EPN, Thu Feb 11 13:32:34 2016
#
# Purpose:    For a given accession, retreive lengths and coorindate
#             strings of all features.
#
# Arguments:
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
#################################################################
#
# Subroutines related to parsing dash-coords (non-NCBI) coordinate strings:
#   dashCoordsStringCommaDelimitedToLength()
#   dashCoordsToLength()
#
#################################################################
# Subroutine: dashCoordsStringCommaDelimitedToLength
# Incept:     EPN, Fri Mar  4 15:24:25 2016
#
# Purpose:    Given a string with >= 1 'dash coordinate' strings
#             separated by commas (e.g. "10-100,101-250") return the
#             total length of all coordinates summed together.  Die if
#             the string is not in the correct format.
#
# Arguments:
#   $dash_coords_str: the dashed coords string
#   $caller_sub_name: name of caller, can be undef
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
# 
# Returns: total length
#
# Dies:    - if we can't parse $dash_coords_str because it's
#            in an unexpected format
#          - if any coordinate is negative
#
#################################################################
sub dashCoordsStringCommaDelimitedToLength {
  my $sub_name = "dashCoordsStringCommaDelimitedToLength()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($dash_coords_str, $caller_sub_name, $FH_HR) = (@_);
  
  my $len = 0;
  my @start_stop_A = split(",", $dash_coords_str);
  foreach my $start_stop (@start_stop_A) { 
    $len += dashCoordsToLength($start_stop, $caller_sub_name, $FH_HR); 
    # dashCoordsToLength() will die if $start_stop has a negative index or is incorrectly formatted in some other way
  }
  
  return $len;
}

#################################################################
# Subroutine: dashCoordsToLength
# Incept:     EPN, Fri Mar  4 15:29:07 2016
#
# Purpose:    Given a string with 1 'dash coordinate' 
#             (e.g. "10-100), return the length implied.
#             Die if either start or stop is negative or
#             if there's something else is wrong with the
#             format of $dash_coords
# Arguments:
#   $start_stop: the dashed coords string
#   $caller_sub_name: name of caller, can be undef
#   $FH_HR:      REF to hash of file handles, including "log" and "cmd"
# 
# Returns: length
#
# Dies:    - if we can't parse $start_stop because it's
#            in an unexpected format
#          - if any coordinate is negative
#
#################################################################
sub dashCoordsToLength { 
  my $sub_name = "dashCoordsToLength";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($start_stop, $caller_sub_name, $FH_HR) = (@_);

  my $len = undef;

  # make sure that we $start_stop is not valid except for the fact that
  # $start or $stop is a negative position (we do not allow negative positions
  # in this function)
  if($start_stop =~ m/^\-+\d+\-\-\d+$/) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, %sstart and stop positions are negative in coords string $start_stop", 
                        (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), 1, $FH_HR); 
  }
  elsif($start_stop =~ m/^\-\d+\-\d+$/) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, %sstart position is negative in coords string $start_stop", 
                        (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), 1, $FH_HR); 
  }
  elsif($start_stop =~ m/^\d+\-\-\d+$/) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, %sstop position is negative in coords string $start_stop", 
                        (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), 1, $FH_HR); 
  }

  # if we get here, $start_stop is either valid, or invalid for a reason other
  # than having a negative position
  if($start_stop =~ m/^(\d+)\-(\d+)$/) { 
    # $start_stop is valid
    $len = (abs($1 - $2) + 1);
  }
  else { 
    # $start_stop is not valid, for some reason other than just having a negative position
    DNAORG_FAIL("ERROR in $sub_name, called by $caller_sub_name, unable to parse start-stop string: $start_stop", 1, $FH_HR); 
  }

  return $len;
}

#################################################################
#################################################################
#
# Subroutines related to output:
#   outputProgressPrior()
#   outputProgressComplete()
#   outputConclusionAndCloseFiles()
#   outputTiming()
#   outputString()
#   outputBanner()
# 
#################################################################
# Subroutine : outputProgressPrior()
# Incept:      EPN, Fri Feb 12 17:22:24 2016
#
# Purpose:      Output to $FH1 (and possibly $FH2) a message indicating
#               that we're about to do 'something' as explained in
#               $outstr.  
#
#               Caller should call *this* function, then do
#               the 'something', then call outputProgressComplete().
#
#               We return the number of seconds since the epoch, which
#               should be passed into the downstream
#               outputProgressComplete() call if caller wants to
#               output running time.
#
# Arguments: 
#   $outstr:     string to print to $FH
#   $progress_w: width of progress messages
#   $FH1:        file handle to print to
#   $FH2:        another file handle to print to, can be undef
# 
# Returns:     Number of seconds and microseconds since the epoch.
#
################################################################# 
sub outputProgressPrior { 
  my $nargs_expected = 4;
  my $sub_name = "outputProgressPrior()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($outstr, $progress_w, $FH1, $FH2) = @_;

  if(defined $FH1) { printf $FH1 ("# %-*s ... ", $progress_w, $outstr); }
  if(defined $FH2) { printf $FH2 ("# %-*s ... ", $progress_w, $outstr); }

  return secondsSinceEpoch();
}

#################################################################
# Subroutine : outputProgressComplete()
# Incept:      EPN, Fri Feb 12 17:28:19 2016
#
# Purpose:     Output to $FH1 (and possibly $FH2) a 
#              message indicating that we've completed 
#              'something'.
#
#              Caller should call *this* function,
#              after both a call to outputProgressPrior()
#              and doing the 'something'.
#
#              If $start_secs is defined, we determine the number
#              of seconds the step took, output it, and 
#              return it.
#
# Arguments: 
#   $start_secs:    number of seconds either the step took
#                   (if $secs_is_total) or since the epoch
#                   (if !$secs_is_total)
#   $extra_desc:    extra description text to put after timing
#   $FH1:           file handle to print to
#   $FH2:           another file handle to print to, can be undef
# 
# Returns:     Number of seconds the step took (if $secs is defined,
#              else 0)
#
################################################################# 
sub outputProgressComplete { 
  my $nargs_expected = 4;
  my $sub_name = "outputProgressComplete()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($start_secs, $extra_desc, $FH1, $FH2) = @_;

  my $total_secs = undef;
  if(defined $start_secs) { 
    $total_secs = secondsSinceEpoch() - $start_secs;
  }

  if(defined $FH1) { printf $FH1 ("done."); }
  if(defined $FH2) { printf $FH2 ("done."); }

  if(defined $total_secs || defined $extra_desc) { 
    if(defined $FH1) { printf $FH1 (" ["); }
    if(defined $FH2) { printf $FH2 (" ["); }
  }
  if(defined $total_secs) { 
    if(defined $FH1) { printf $FH1 (sprintf("%.1f seconds%s", $total_secs, (defined $extra_desc) ? ", " : "")); }
    if(defined $FH2) { printf $FH2 (sprintf("%.1f seconds%s", $total_secs, (defined $extra_desc) ? ", " : "")); }
  }
  if(defined $extra_desc) { 
    if(defined $FH1) { printf $FH1 $extra_desc };
    if(defined $FH2) { printf $FH2 $extra_desc };
  }
  if(defined $total_secs || defined $extra_desc) { 
    if(defined $FH1) { printf $FH1 ("]"); }
    if(defined $FH2) { printf $FH2 ("]"); }
  }

  if(defined $FH1) { printf $FH1 ("\n"); }
  if(defined $FH2) { printf $FH2 ("\n"); }
  
  return (defined $total_secs) ? $total_secs : 0.;
}

#######################################################################
# Subroutine: outputConclusionAndCloseFiles()
# Incept:     EPN, Thu Nov  5 18:25:31 2009 [ssu-align]
# 
# Purpose:    Output a list of the main output files created 
#             and the final few lines of output and optionally the 
#             run time timing to the summary file. Print date and
#             system information to the log file. 
#
#             Close all open file handles.
#
# Arguments: 
#  $total_secs:            total number of seconds, "" to not print timing
#  $odir:                  output directory, if "", files were put in cwd
#  $ofile_info_HHR:        REF to the 2D hash of output file information
#
# Returns:   Nothing.
# 
# Dies:      Never.
#
####################################################################
sub outputConclusionAndCloseFiles { 
  my $nargs_expected = 3;
  my $sub_name = "outputConclusionAndCloseFiles()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($total_secs, $odir, $ofile_info_HHR) = @_;

  validateOutputFileInfoHashOfHashes($ofile_info_HHR);

  my $log_FH  = $ofile_info_HHR->{"FH"}{"log"};
  my $cmd_FH  = $ofile_info_HHR->{"FH"}{"cmd"};
  my $list_FH = $ofile_info_HHR->{"FH"}{"list"};

  my $key2d; # a key in the 2nd dimension of $ofile_info_HHR

  # output the list of files that we created for which the 'mainout' variable 
  # is 1 to the $log file and stdout, we already printed the descriptions
  # to the list file in helperAddFileToOutputInfo().
  if(defined $log_FH) { 
    outputString($log_FH, 1, sprintf("#\n"));
    # create a temporary array with description of files with 'outmain' set to 1 (we'll only print these)
    # so we get pretty formatting
    my @tmp_A = ();
    foreach $key2d (keys (%{$ofile_info_HHR->{"desc"}})) { 
      if($ofile_info_HHR->{"mainout"}{$key2d}) { 
        push(@tmp_A, $ofile_info_HHR->{"desc"}{$key2d});
      }
    }
    my $width_desc = length("# ") + maxLengthScalarValueInArray(\@tmp_A) + length(" saved in:");
    my $cur_idx = 1;
    my $num_ofile = validateOutputFileInfoHashOfHashes($ofile_info_HHR); # this function validates we have exactly 1 of each "order" value of 1..$num_ofile
    for(my $i = 1; $i <= $num_ofile; $i++) { 
      foreach $key2d (keys (%{$ofile_info_HHR->{"order"}})) { 
        if(($ofile_info_HHR->{"order"}{$key2d} == $i) && 
           ($ofile_info_HHR->{"mainout"}{$key2d})) { # only print out files for which the special "mainout" value is '1'
          outputString($log_FH, 1, sprintf("# %-*s %s\n", $width_desc, $ofile_info_HHR->{"desc"}{$key2d} . " saved in:", $ofile_info_HHR->{"nodirpath"}{$key2d}));
        }
      }
    }
    outputString($log_FH, 1, sprintf("#\n"));
    outputString($log_FH, 1, sprintf("# All output files created in %s\n", ($odir eq "") ? "the current working directory" : "directory \.\/$odir\/"));
    outputString($log_FH, 1, sprintf("#\n"));
    if($total_secs ne "") { # don't print this if rvr-align is caller
      outputTiming("# CPU time: ", $total_secs, 1, $log_FH); 
      outputString($log_FH, 1, "#            hh:mm:ss\n");
      outputString($log_FH, 1, "# \n");
      outputString($log_FH, 1, "# DNAORG-SUCCESS\n");
    }
  }
  if(defined $cmd_FH) { 
    outputString($cmd_FH, 0, "# " . `date`);      # prints date,        e.g.: 'Mon Feb 22 16:37:09 EST 2016'
    outputString($cmd_FH, 0, "# " . `uname -a`);  # prints system info, e.g.: 'Linux cbbdev13 2.6.32-573.7.1.el6.x86_64 #1 SMP Tue Sep 22 22:00:00 UTC 2015 x86_64 x86_64 x86_64 GNU/Linux'
    if($total_secs ne "") { # don't print this if rvr-align is caller
      outputString($cmd_FH, 0, "# DNAORG-SUCCESS\n");
    }
  }

  # close any open file handles
  foreach $key2d (keys (%{$ofile_info_HHR->{"FH"}})) { 
    if(defined $ofile_info_HHR->{"FH"}{$key2d}) { 
      close $ofile_info_HHR->{"FH"}{$key2d};
    }
  }

  return;
}

#####################################################################
# Subroutine: outputTiming()
# Incept:     EPN, Tue Jun 16 08:52:08 2009 [ssu-align]
# 
# Purpose:    Output elapsed time in hhhh:mm:ss format.
# 
# Arguments:
#   $prefix:               string to print before the hhhh:mm:ss time info.
#   $inseconds:            number of seconds
#   $print_to_stdout:      '1' to print to stdout, '0' not to
#   $FH:                   file handle to print to
#
# Returns:    Nothing, if it returns, everything is valid.
# 
####################################################################
sub outputTiming { 
  my $nargs_expected = 4;
  my $sub_name = "outputTiming()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($prefix, $inseconds, $print_to_stdout, $FH) = @_;

  my $time_str = formatTimeString($inseconds);
  outputString($FH, $print_to_stdout, $prefix . " " . $time_str . "\n"); # this will always end with a newline
  
  return;
}

###########################################################
# Subroutine: outputString()
# Incept: EPN, Wed Oct 29 20:42:16 2014 [rnavore]
#
# Purpose: Given a string and an open file handle <$FH>, 
#          output the string to the file handle 
#          and potentially to stdout as well. 
#          If <$FH> is not defined then do not 
#          print to a file. 
#
# Arguments:
#   $FH:              file handle to output to, can be undef
#   $print_to_stdout: if '1' also output string to stdout
#   $string:          the string to output
#
# Returns: Nothing. 
#
# Dies:    Never.
#
###########################################################
sub outputString {
  my $nargs_expected = 3;
  my $sub_name = "outputString()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($FH, $print_to_stdout, $string) = @_;

  if(defined $FH)      { print $FH $string; }
  if($print_to_stdout) { print     $string; }

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
#    $version:           version of dnaorg
#    $releasedate:       month/year of version (e.g. "Feb 2016")
#    $synopsis:          string reporting the date
#    $date:              date information to print
#    $dnaorgdir:         dnaorg directory
#
# Returns:    Nothing, if it returns, everything is valid.
# 
# Dies: never
####################################################################
sub outputBanner {
  my $nargs_expected = 6;
  my $sub_name = "outputBanner()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($FH, $version, $releasedate, $synopsis, $date, $dnaorgdir) = @_;

  print $FH ("\# $synopsis\n");
  print $FH ("\# dnaorg $version ($releasedate)\n");
#  print $FH ("\# Copyright (C) 2014 HHMI Janelia Research Campus\n");
#  print $FH ("\# Freely distributed under the GNU General Public License (GPLv3)\n");
  print $FH ("\# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  if(defined $date)      { print $FH ("# date:       $date\n"); }
  if(defined $dnaorgdir) { print $FH ("# \$DNAORGDIR: $dnaorgdir\n"); }
  printf $FH ("#\n");

  return;
}

#################################################################
# Subroutine : outputDividingLine()
# Incept:      EPN, Tue Apr 12 14:40:13 2016
#
# Purpose:     Print a line of dashes followed by single spaces
#              with $ndash dashes to file handle $FH.
#              if $ndash is undefined, set it to 66.
#
# Arguments: 
#   $ndashes:  number of dashes in output dividing line
#   $FH:       file handle to print to
# 
# Returns:     Nothing.
# 
# Dies:        Never.
#
################################################################# 
sub outputDividingLine { 
  my $nargs_expected = 2;
  my $sub_name = "outputDividingLine()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ndashes, $FH) = @_;

  if(! defined $ndashes) { 
    $ndashes = 66;
  }

  my $div_line = "#";
  for(my $i = 0; $i < $ndashes; $i++) { 
    $div_line .= " -";
  }
  $div_line .= "\n";

  print $FH $div_line;
  
  return;
}

#################################################################
#################################################################
#
# Subroutines for dumping data structures, usually for debugging:
#   dumpInfoHashOfArrays()
#   dumpArrayOfHashesOfHashes()
#   dumpArrayOfHashes()
#
#################################################################
# Subroutine: dumpInfoHashOfArrays()
# Incept:     EPN, Thu Feb 11 15:06:31 2016
#
# Purpose:    Print an 'info' hash of arrays, probably for 
#             debugging purposes.
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
# Subroutine: dumpArrayOfHashesOfHashes()
# Incept:     EPN, Fri Mar  4 16:02:28 2016
#
# Purpose:    Dump the contents of an array of hashes of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $AHHR:        ref of the array of hashes of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub dumpArrayOfHashesOfHashes { 
  my $sub_name = "dumpArrayOfHashesOfHashes()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $AHHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $nel1 = scalar(@{$AHHR});
  for(my $i1 = 0; $i1 < $nel1; $i1++) { 
    printf $FH ("*A*HH el %2d\n", ($i1+1)); 
    my $nel2 = scalar(keys %{$AHHR->[$i1]}); 
    my $i2 = 0;
    foreach my $key2 (sort keys %{$AHHR->[$i1]}) { 
      printf("\tA*H*H el %2d key: $key2\n", ($i2+1)); 
      $i2++;
      my $nel3 = scalar(keys %{$AHHR->[$i1]{$key2}});
      my $i3 = 0;
      foreach my $key3 (sort keys %{$AHHR->[$i1]{$key2}}) { 
        printf("\tAH*H* el %2d key: $key3 value: %s\n", ($i3+1), $AHHR->[$i1]{$key2}{$key3}); 
        $i3++;
      }
      printf $FH ("\n");
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: dumpArrayOfHashes()
# Incept:     EPN, Thu Feb  8 11:01:29 2018
#
# Purpose:    Dump the contents of an array of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $AHR:         ref of the array of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub dumpArrayOfHashes { 
  my $sub_name = "dumpArrayOfHashes()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $AHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $nel1 = scalar(@{$AHR});
  for(my $i1 = 0; $i1 < $nel1; $i1++) { 
    printf $FH ("*A*H el %2d\n", ($i1+1)); 
    my $nel2 = scalar(keys %{$AHR->[$i1]}); 
    my $i2 = 0;
    foreach my $key2 (sort keys %{$AHR->[$i1]}) { 
      printf("\tA*H* el %2d key: $key2 value: %s\n", ($i2+1), $AHR->[$i1]{$key2}); 
      $i2++;
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
#################################################################
#
# Subroutines for validating the special data structures:
#   validateExecutableHash()
#   validateFeatureInfoHashIsComplete()
#   validateModelInfoHashIsComplete()
#   validateSequenceInfoHashIsComplete()
#   validateErrorInfoHashIsComplete()
#   validateInfoHashOfArraysIsComplete()
#   validateOutputFileInfoHashOfHashes()
#   validateAndGetSizeOfInfoHashOfArrays()
#   getConsistentSizeOfInfoHashOfArrays()
#   validateFTableErrorExceptions()
#
#################################################################
# Subroutine : validateExecutableHash()
# Incept:      EPN, Sat Feb 13 06:27:51 2016
#
# Purpose:     Given a reference to a hash in which the 
#              values are paths to executables, validate
#              those files are executable.
#
# Arguments: 
#   $execs_HR: REF to hash, keys are short names to executable
#              e.g. "cmbuild", values are full paths to that
#              executable, e.g. "/usr/local/infernal/1.1.1/bin/cmbuild"
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns:     void
#
# Dies:        if one or more executables does not exist#
#
################################################################# 
sub validateExecutableHash { 
  my $nargs_expected = 2;
  my $sub_name = "validateExecutableHash()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($execs_HR, $FH_HR) = @_;

  my $fail_str = undef;
  foreach my $key (sort keys %{$execs_HR}) { 
    if(! -e $execs_HR->{$key}) { 
      $fail_str .= "\t$execs_HR->{$key} does not exist.\n"; 
    }
    elsif(! -x $execs_HR->{$key}) { 
      $fail_str .= "\t$execs_HR->{$key} exists but is not an executable file.\n"; 
    }
  }
  
  if(defined $fail_str) { 
    DNAORG_FAIL("ERROR in $sub_name(),\n$fail_str", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: validateFeatureInfoHashIsComplete()
# Incept:     EPN, Tue Feb 16 10:15:45 2016
#
# Purpose:    Validate that a 'feature info' hash is valid and complete.
#             'Complete' means it has all the expected keys, each of which is an identically sized array.
#             The expected keys are:
#                "annot_type":    annotation type:
#                                 "model":        if this feature's annotation is derived from models (type: 'mp' or 'cds-notmp'
#                                 "multifeature": if this feature's annotation is derived from multiple other features (type: 'cds-mp')
#                                 "duplicate":    if this feature's annotation is derived by copying data from another (e.g., type: 'gene')
#                "append_num":    number of nucleotides to append for this model, usually 0, unless this
#                                 model is the final mature peptide of a CDS, in which case it is usually 3,
#                                 this informs other parts of the pipeline that we need to append the 3 
#                                 downstream nucleotides: the stop codon, which is not included in the 
#                                 annotation of the final mature peptide in NCBI, but is included in the
#                                 parent CDS annotation.
#                "filename_root": 'root' string for output file names related to this model: 
#                "final_mdl":     index (in arrays of %mdl_info_HA) of final model for this feature
#                "first_mdl":     index (in arrays of %mdl_info_HA) of first model for this feature
#                "nmodels":       number of models for this feature (e.g. number of exons or segments) for this feature, 
#                                 0 if 'annotation' eq 'multifeature'
#                "out_product":   output value: name of product for this feature (e.g. "replication-associated protein")
#                "out_gene":      output value: name of gene for this feature (e.g. "ORF2")
#                "out_exception": output value: name of exception for this feature (e.g. "RNA editing")
#                "out_short":     output value: short name for this feature (e.g. "CDS #4 [1 exon; +]")
#                "out_tiny":      output value: very short name for this feature (e.g. "CDS#4")
#                "ref_coords":    coordinates for this feature in the reference
#                "ref_len":       length of this feature in the reference
#                "ref_strand":    strand for this feature in the reference
#                "source_idx":    if "annot_type" is "duplicate", the index of the feature that
#                                 that is the source to copy for this feature's annotation.
#                "type":          type of feature: "mp", "cds-notmp", "cds-mp", "xfeat", or "dfeat"
#                "type_fname":    string for naming output files related to this feature, e.g. "mp", "cds", "$xfeat", or "$dfeat"
#                "type_ftable":   output feature table feature name, e.g. "mat_peptide", "CDS", "$xfeat", or "$dfeat"
#                "type_idx":      index for the type that this feature is (e.g. 4, if this is the 4th "cds-notmp" feature
#              
#             If @{exceptions_AR} is non-empty, then keys in 
#             in that array need not be in %{$ftr_info_HAR}.
#
# Arguments:
#   $ftr_info_HAR:  REF to hash of arrays of feature information
#   $exceptions_AR: REF to array of keys that may be excluded from the hash
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: Number of elements in each and every array of %{$ftr_info_HAR}
#
# Dies:    - if one of the expected keys (listed above and not in @{$exceptions_AR})
#            does not exist in $ftr_info_HAR
#          - if two arrays in $ftr_info_HAR are of different sizes
#          - if any key listed in @{$exceptions_AR} is not one of the expected keys
#
#################################################################
sub validateFeatureInfoHashIsComplete { 
  my $sub_name = "validateFeatureInfoHashIsComplete()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($ftr_info_HAR, $exceptions_AR, $FH_HR) = (@_);
  
  my @expected_keys_A = ("append_num", "final_mdl", "first_mdl", "annot_type", "nmodels", 
                         "out_product", "out_gene", "out_exception", 
                         "out_tiny", "out_short", "out_tiny", "ref_coords",
                         "ref_len", "ref_strand", "type", "type_idx", 
                         "parent_ftr", "primary_children_ftr_str", "primary_children_ftr_num", 
                         "all_children_ftr_str", "all_children_ftr_num", "source_idx");

  return validateInfoHashOfArraysIsComplete($ftr_info_HAR, \@expected_keys_A, $exceptions_AR, $FH_HR);
}

#################################################################
# Subroutine: validateModelInfoHashIsComplete()
# Incept:     EPN, Tue Feb 16 11:00:03 2016
#
# Purpose:    Validate that a 'model info' hash is valid and complete.
#             'Complete' means it has all the expected keys, each of which is an identically sized array.
#             The expected keys are:
#                "append_num":    number of nucleotides to append for this model, usually 0, unless this
#                                 model is the final mature peptide of a CDS, in which case it is usually 3,
#                                 this informs other parts of the pipeline that we need to append the 3 
#                                 downstream nucleotides: the stop codon, which is not included in the 
#                                 annotation of the final mature peptide in NCBI, but is included in the
#                                 parent CDS annotation.
#                "checksum":      checksum of the 'alignment' (single sequence) file the model was built from
#                "cmname":        name of the model, used in infernal output 
#                "filename_root": 'root' string for output file names related to this model: 
#                "idx_ajb_str":   string of model indices (each separated by a single ",") that are adjacent 'before'
#                                 this model index
#                "idx_aja_str":   string of model indices (each separated by a single ",") that are adjacent 'after'
#                                 this model index
#                "idx_olp_str":   string of model indices (each separated by a single ",") that overlap by >= 1 nt on
#                                 same strand with this model index
#                "is_final":      '1' if this model is the final model (e.g. final exon) for the feature it models ("map_ftr")
#                "is_first":      '1' if this model is the first model (e.g. final exon) for the feature it models ("map_ftr")
#                "length":        length, in nucleotides, of the model
#                "ref_start":     start position of modelled region in the reference genome
#                "ref_stop":      stop position of modelled region in the reference genome
#                "ref_strand":    strand of modelled region in the reference genome
#                "map_exon":      the exon index this model models (1.."map_nexon" value)
#                "map_ftr":       the feature index (array index in ftr_info_HAR) this model models
#                "map_nexon":     the number of exons the feature this model models has 
#                "out_ajb_str":   string of short model descriptions that are adjacent 'before'
#                                 this model index
#                "out_aja_str":   string of short model descriptions that are adjacent 'after'
#                                 this model index
#                "out_olp_str":   string of short model descriptions that overlap by >= 1 nt on
#                                 same strand with this model index
#                "out_tiny":      output value: very short name for this model (e.g. "CDS#4.2")
#                "out_idx":       output value: feature index and segment index this model (e.g. "4.2")
#
#             If @{exceptions_AR} is non-empty, then keys in 
#             in that array need not be in %{$mdl_info_HAR}.
#
# Arguments:
#   $mdl_info_HAR:  REF to hash of arrays of model information
#   $exceptions_AR: REF to array of keys that may be excluded from the hash
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: Number of elements in each and every array of %{$mdl_info_HAR}
#
# Dies:    - if one of the expected keys (listed above and not in @{$exceptions_AR})
#            does not exist in $mdl_info_HAR
#          - if two arrays in $mdl_info_HAR are of different sizes
#          - if any key listed in @{$exceptions_AR} is not one of the expected keys
#
#################################################################
sub validateModelInfoHashIsComplete { 
  my $sub_name = "validateModelInfoHashIsComplete()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($mdl_info_HAR, $exceptions_AR, $FH_HR) = (@_);
  
  my @expected_keys_A = ("append_num", "checksum", "cmname", "filename_root", "is_final", "is_first", "idx_ajb_str", 
                         "idx_aja_str", "idx_olp_str", "length", "out_aja_str", "out_ajb_str", "out_olp_str",
                         "ref_start", "ref_stop", "ref_strand", "map_exon", "map_ftr", "map_nexon", "out_tiny", "out_idx");

  return validateInfoHashOfArraysIsComplete($mdl_info_HAR, \@expected_keys_A, $exceptions_AR, $FH_HR);
}

#################################################################
# Subroutine: validateSequenceInfoHashIsComplete()
# Incept:     EPN, Tue Mar 15 05:35:22 2016
#
# Purpose:    Validate that a 'sequence info' hash is valid and complete.
#             'Complete' means it has all the expected keys, each of which is an identically sized array.
#             The expected keys are:
#                "seq_name":    name of the sequence in the sequence file we create and search in
#                "seq_len":     length of the sequence with name in "seq_name"in the sequence file we create and search in
#                "accn_name":   accession of the sequence in GenBank
#                "accn_len":    length of the sequence in GenBank, will be same as value in "seq_len" only if
#                               the -c option is not used to indicate a circular genome, if the -c option is
#                               used, then this will equal 2*seq_len.
#
#             This function also validates that one of the following is true:
#             1) -c is 'off' in %{$opt_HHR} and all "seq_len" and "accn_len" values are equal
#             2) -c is 'on'  in %{$opt_HHR} and all "seq_len" values are 2 * the "accn_len" values
#
#             If @{exceptions_AR} is non-empty, then keys in 
#             in that array need not be in %{$seq_info_HAR}.
#
# Arguments:
#   $seq_info_HAR:  REF to hash of arrays of sequence information
#   $exceptions_AR: REF to array of keys that may be excluded from the hash
#   $opts_HHR:      REF to the 2D hash of command line options
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: Number of elements in each and every array of %{$seq_info_HAR}
#
# Dies:    - if one of the expected keys (listed above and not in @{$exceptions_AR})
#            does not exist in $seq_info_HAR
#          - if two arrays in $seq_info_HAR are of different sizes
#          - if any key listed in @{$exceptions_AR} is not one of the expected keys
#          - if -c enabled     but not all seq_len values are 2X corresponding accn_len values
#          - if -c not enabled but not all seq_len values are equal to corresponding accn_len values
#
#################################################################
sub validateSequenceInfoHashIsComplete { 
  my $sub_name = "validateSequenceInfoHashIsComplete()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($seq_info_HAR, $exceptions_AR, $opt_HHR, $FH_HR) = (@_);
  
  my @expected_keys_A = ("seq_name", "seq_len", "accn_name", "accn_len");

  my $nseq = validateInfoHashOfArraysIsComplete($seq_info_HAR, \@expected_keys_A, $exceptions_AR, $FH_HR);
  # above call will die if we are invalid

  # make sure seq_len is either 2X accn_len (if -c) or equal to accn_len (if ! -c)
  if(opt_Get("-c", $opt_HHR)) { 
    # -c option on, seq_len should be 2X accn_len
    for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      if($seq_info_HAR->{"seq_len"}[$seq_idx] != (2 * $seq_info_HAR->{"accn_len"}[$seq_idx])) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, -c option on, but for sequence %s, seq_len value of %d not equal to 2X accn_len of %d", 
                            $seq_info_HAR->{"seq_name"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $seq_info_HAR->{"accn_len"}[$seq_idx]),
                    1, $FH_HR);
      }
    }
  }
  else { 
    # -c option off, seq_len should be accn_len
    for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      if($seq_info_HAR->{"seq_len"}[$seq_idx] != $seq_info_HAR->{"accn_len"}[$seq_idx]) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, -c option off, but for sequence %s, seq_len value of %d not equal to accn_len of %d", 
                            $seq_info_HAR->{"seq_name"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $seq_info_HAR->{"accn_len"}[$seq_idx]),
                    1, $FH_HR);
      }
    }
  }

  # make sure we do not have any duplicate accessions/names
  my %name_H = ();
  my %accn_H = ();
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
    my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
    if(exists $name_H{$seq_name}) { 
      DNAORG_FAIL("ERROR in $sub_name, sequence name $seq_name exists twice", 1, $FH_HR);
    }
    if(exists $accn_H{$accn_name}) { 
      DNAORG_FAIL("ERROR in $sub_name, accession $accn_name exists twice", 1, $FH_HR);
    }
    $name_H{$seq_name} = 1;
    $accn_H{$accn_name} = 1;
  }

  return $nseq;
}

#################################################################
# Subroutine: validateErrorInfoHashIsComplete()
# Incept:     EPN, Fri Mar  4 12:56:43 2016
#
# Purpose:    Validate that a 'error info' hash is valid and complete.
#             'Complete' means it has all the expected keys, each of which is an identically sized array.
#             The expected keys are:
#                "code":            the error code, e.g. "nop"
#                "pertype":         type of this error, either "feature" (error occurs for a single feature) or "sequence"
#                                   (error occurs for an entire sequence)
#                "maybe_allowed":   '1' if a value of 'maybe' is allowed for this error code, '0' if not
#                "desc":            the description/message for this error, e.g. ""unable to identify homologous feature"
#                "incompat":        string that lists all error code indices (separated by commas) the current error code is incompatible with
#                "requires":        string that lists all error code indices (separated by commas) the current error code requires
#             Information relevant  to feature table output:
#                "ftbl_valid":      '1' if this error is valid/relevant to feature table output
#                                   '0' if it is not valid/relevant and does not change feature table output
#                "ftbl_invalid_by": string that lists all error code indices (separated by commas) the current error code is invalidated by
#                                   for feature table output (if any of those errors are also present, this error becomes invalid and does
#                                   not impact feature table output)
#                "ftbl_pred_stop":  '1' to use predicted stop (instead of corrected one) in feature table when this error exists
#                "ftbl_note":       note message for feature table
#                "ftbl_err"         ERROR message for feature table
#                
#               
#             If @{exceptions_AR} is non-empty, then keys in 
#             in that array need not be in %{$err_info_HAR}.
#
# Arguments:
#   $err_info_HAR:  REF to hash of arrays of error information
#   $exceptions_AR: REF to array of keys that may be excluded from the hash
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: Number of elements in each and every array of %{$err_info_HAR}
#
# Dies:    - if one of the expected keys (listed above and not in @{$exceptions_AR})
#            does not exist in $err_info_HAR
#          - if two arrays in $err_info_HAR are of different sizes
#          - if any other key other than those listed above exist in ${%err_info_HAR}
#          - if any key listed in @{$exceptions_AR} is not one of the expected keys
#
#################################################################
sub validateErrorInfoHashIsComplete { 
  my $sub_name = "validateErrorInfoHashIsComplete()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($err_info_HAR, $exceptions_AR, $FH_HR) = (@_);
  
  my @expected_keys_A = ("code", "pertype", "maybe_allowed", "desc", "incompat", "requires", "ftbl_valid", "ftbl_invalid_by", "ftbl_pred_stop", "ftbl_note", "ftbl_err");

  # we do not do any more extensive checking because the function that adds elements to the error info hash
  # (addToErrorInfoHash()) does a lot of checking to validate the data before it is added

  return validateInfoHashOfArraysIsComplete($err_info_HAR, \@expected_keys_A, $exceptions_AR, $FH_HR);
}

#################################################################
# Subroutine: validateInfoHashOfArraysIsComplete()
# Incept:     EPN, Tue Feb 16 11:10:40 2016
#
# Purpose:    Validate that a 'info' hash is valid and complete.
#             'Complete' means it has all the expected keys, each of 
#             which is an identically sized array.
#             The expected keys are passed in in @{$expected_AR}.
#             If @{exceptions_AR} is non-empty, then keys in 
#             in that array do not need to be in %{$info_HAR},
#             but they can be.
#
# Arguments:
#   $info_HAR:         REF to hash of arrays of information
#   $expected_keys_AR: REF to array of keys that are expected to be in the hash
#   $exceptions_AR:    REF to array of keys in @{$expected_keys_AR} that may be excluded from the hash
#   $FH_HR:            REF to hash of file handles, including "log" and "cmd"
# 
# Returns: Number of elements in each and every array of %{$info_HAR}
#
# Dies:    - if one of the expected keys (listed above and not in @{$exceptions_AR})
#            does not exist in $info_HAR
#          - if two arrays in $info_HAR are of different sizes
#          - if any other key other than those listed above exist in ${%info_HAR}
#          - if any key listed in @{$exceptions_AR} is not one of the expected keys
#################################################################
sub validateInfoHashOfArraysIsComplete { 
  my $sub_name = "validateInfoHashOfArraysIsComplete()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($info_HAR, $expected_keys_AR, $exceptions_AR, $FH_HR) = (@_);
  
  # make sure our exceptions are actually in the expected array
  if(defined $exceptions_AR) { 
    foreach my $key (@{$exceptions_AR}) { 
      if(findNonNumericValueInArray($expected_keys_AR, $key, $FH_HR) == -1) { 
        DNAORG_FAIL("ERROR in $sub_name, excepted value $key is not an expected key in the feature info hash", 1, $FH_HR);
      }
    }
  }

  # make the list of keys we'll require, this is all expected keys minus any exceptions in @{$exceptions_AR}
  my @reqd_keys_A = ();
  foreach my $key (@{$expected_keys_AR}) { 
    if((! defined $exceptions_AR) || (findNonNumericValueInArray($key, $exceptions_AR, $FH_HR) == -1)) { 
      push(@reqd_keys_A, $key);
    }
  }
                         
  my $nftr = validateAndGetSizeOfInfoHashOfArrays($info_HAR, \@reqd_keys_A, $FH_HR);  

  return $nftr;
}

#################################################################
# Subroutine: validateOutputFileInfoHashOfHashes()
# Incept:     EPN, Mon Feb 29 09:21:51 2016
#
# Purpose:    Validate an 'output file' info hash of hashes.
#             A valid info hash of hashes has the same set of 2d
#             keys for each 1d key except for "FH". The set of 1d keys is 
#             "order": integer, the order in which this element was added
#             "fullpath":  full path to the file 
#             "nodirpath": file name, "fullpath" minus directories
#             "mainout":   '1' if this file should be listed in the main output,
#                          '0' if it should only be listed in the .list file
#             "desc":      short description of the file
#             "FH":        open file handle for this file, or undef             
#
#             For "FH", the value for each $ofile_info_HH{"FH"}{$key2d} can
#             be either defined or not defined.
#
# Arguments:
#   $ofile_info_HHR:  REF to hash of hashes of output file information
# 
# Returns: Number of elements in each and every 2d hash (except possibly %{$ofile_info_HHR->{"FH"}})
#
# Dies:    - if one of the expected keys (listed above) does not exist in $ofile_info_HHR
#          - if two 2d hashes in $ofile_info_HHR (besides %{$ofile_info_HHR->{"FH"}}) are of different sizes
#          - if two 2d hashes in $ofile_info_HHR (besides %{$ofile_info_HHR->{"FH"}}) have different set of keys
#################################################################
sub validateOutputFileInfoHashOfHashes { 
  my $sub_name = "validateOutputFileInfoHashOfHashes()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($ofile_info_HHR) = (@_);
  
  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my @same_keys1d_A = ("order", "fullpath", "nodirpath", "mainout", "desc"); # all of these 2nd dim hashes should have same set of keys
  my @all_keys1d_A   = (@same_keys1d_A, "FH");             # all 1d keys
  my $i;     # a counter
  my $key1d; # a 1st dim key
  my $key2d; # a 2nd dim key

  # make sure we don't have any extra 1d keys we don't expect
  foreach $key1d (keys %{$ofile_info_HHR}) { 
    my $found_it = 0;
    foreach my $tmp_key1d (@all_keys1d_A) { 
      if($key1d eq $tmp_key1d) { 
        $found_it = 1;
      }
    }
    if($found_it == 0) { 
      DNAORG_FAIL("ERROR in $sub_name, unexpected 1d key $key1d exists.", 1, $FH_HR);
    }     
  } 
  
  # make sure all 2nd dim keys for all 1st dim keys are the same as the 2nd dim keys for 1st dim key "order"
  if(! defined $ofile_info_HHR->{"order"}) { 
    DNAORG_FAIL("ERROR in $sub_name, expected 1d key order does not exist.", 1, $FH_HR);
  }
  foreach my $key1d (@same_keys1d_A) { 
    if($key1d ne "order") { # skip "order"
      if(! defined $ofile_info_HHR->{$key1d}) { 
        DNAORG_FAIL("ERROR in $sub_name, expected 1d key $key1d does not exist.", 1, $FH_HR);
      }
      # we make sure the set of 2d keys in $ofile_info_HHR->{"order"} and $ofile_info_HHR->{$key1d} are 
      # identical in 2 steps:
      # 1) make sure all 2d keys from $ofile_info_HHR->{"order"} are also in $ofile_info_HHR->{"order"}
      foreach $key2d (keys %{$ofile_info_HHR->{"order"}}) { 
        if(! defined $ofile_info_HHR->{$key1d}{$key2d}) { 
          DNAORG_FAIL("ERROR in $sub_name, 2nd dim key $key2d exists for ofile_info_HHR->{order} but not for ofile_info_HHR->{$key1d}", 1, $FH_HR); 
        }
      }
      # 2) make sure all the 2d keys in $ofile_info_HHR->{$key1d} are also in $ofile_info_HHR->{"order"}
      foreach $key2d (keys %{$ofile_info_HHR->{$key1d}}) { 
        if(! defined $ofile_info_HHR->{"order"}{$key2d}) { 
          DNAORG_FAIL("ERROR in $sub_name, 2nd dim key $key2d exists for ofile_info_HHR->{order} but not for ofile_info_HHR->{$key1d}", 1, $FH_HR); 
        }
      }
    }
  }

  # make sure that $ofile_info_HHR->{"order"} has all values 1..$nkey2d
  my $nkey2d = scalar(keys %{$ofile_info_HHR->{"order"}});
  my @check_A = (); 
  for ($i = 1; $i <= $nkey2d; $i++) { 
    $check_A[$i] = 0; # changed to 1 when we observe it below
  }
  foreach $key2d (keys %{$ofile_info_HHR->{"order"}}) { 
    $check_A[$ofile_info_HHR->{"order"}{$key2d}] = 1;
  }
  for ($i = 1; $i <= $nkey2d; $i++) { 
    if($check_A[$i] != 1) { 
      DNAORG_FAIL("ERROR in $sub_name, invalid values for ofile_info_HH{order}, $nkey2d 2nd dim keys, but value $i does not exist", 1, $FH_HR);
    }
  }

  return $nkey2d;
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
# Returns: number of elements (scalar(@{$HAR->{$key}})) for all keys $key
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

  # get size and check consistency
  my $nel = getConsistentSizeOfInfoHashOfArrays($HAR, $FH_HR); 
  # this will die if not all elements are the same size, or any values are undefined


  return $nel;
}

#################################################################
# Subroutine: getConsistentSizeOfInfoHashOfArrays()
# Incept:     EPN, Fri Mar  4 08:57:17 2016
#
# Purpose:    Return number of hash keys ($nkey) in a 'info hash',
#             which is actually a hash of arrays.
#             No keys are required, but at least one must exist.
#             Die if all existing keys do not have the exact
#             same number of ($nkey) elements in their arrays.
#             Die if any values are undefined.
#
# Arguments:
#   $info_HAR:  REF to hash of arrays of hash
#   $FH_HR:     REF to hash of file handles, including "log" and "cmd"
# 
# Returns: Number of elements in each and every array of %{$info_HAR}
#
# Dies:    - if two arrays in %{$HAR} are of different sizes
#          - if %{$HAR} has no keys
#          - if any value of $HAR{$key}[] is undef 
#################################################################
sub getConsistentSizeOfInfoHashOfArrays { 
  my $sub_name = "getConsistentSizeOfInfoHashOfArrays()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HAR, $FH_HR) = (@_);
  
  if((! %{$HAR}) || (scalar(keys %{$HAR}) == 0)) { 
    DNAORG_FAIL("ERROR in $sub_name, hash of arrays does not exist or has no keys", 1, $FH_HR);
  }

  # get size and check consistency: each array should be the same size
  my $nel = -1;
  my $nkey = scalar(keys %{$HAR});
  my $nel_key = undef;
  my $nel_key_values = ""; # the element values for $HAR->{$nel_key}, only used if we encounter an error
  foreach my $key (sort keys %{$HAR}) { 
    if($nel == -1) { 
      $nel = scalar(@{$HAR->{$key}});
      $nel_key = $key;
      foreach my $el (@{$HAR->{$key}}) { 
        $nel_key_values .= "$el ";
      }
    }
    else { 
      if($nel != scalar(@{$HAR->{$key}})) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, expected number of elements in array for key $key is $nel (from key $nel_key) but %d exist for key $key!\nElement values from key $nel_key are: %s", scalar(@{$HAR->{$key}}), $nel_key_values), 1, $FH_HR);
      }
    }
    # check that all values are defined
    for(my $i = 0; $i < $nel; $i++) { 
      if(! defined $HAR->{$key}[$i]) { 
        DNAORG_FAIL("ERROR in $sub_name, undefined value: key: $key index: $i", 1, $FH_HR); 
      }        
    }
  }

  return $nel;
}

#################################################################
# Subroutine: validateFTableErrorExceptions()
# Incept:     EPN, Thu Feb  8 12:36:37 2018
#
# Purpose:    Validate that every hash element of the array of hashes
#             @{$ftbl_err_exceptions_AHR} has all required keys, which 
#             includes all error codes from %{$err_info_HAR}.
#
# Arguments:
#   $ftbl_err_exceptions_AHR:  REF to array of hashes of feature table error
#                              exceptions
#   $err_info_HAR:             REF to hash of arrays of error information
#   $FH_HR:                    REF to hash of file handles, including "log" 
#                              and "cmd"
# 
# Returns: Number of elements in @{$ftbl_err_exceptions_HAR}
#
# Dies:    If any error combination matches more than one error exception.
#          Actually dies in checkErrorsAgainstFTableErrorExceptions().
#
#################################################################
sub validateFTableErrorExceptions { 
  my $sub_name = "validateFTableErrorExceptions";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($ftbl_err_exceptions_AHR, $err_info_HAR, $FH_HR) = (@_);

  my $nerr = validateErrorInfoHashIsComplete($err_info_HAR, undef, $FH_HR); 
  my @other_reqd_keys_A = ("misc_feature", "pred_stop", "note");

  my $nexc = scalar(@{$ftbl_err_exceptions_AHR});
  # make sure all error codes are hash keys with valid values ("R", "D", or "A")
  for(my $i = 0; $i < $nexc; $i++) { 
    for(my $err_idx = 0; $err_idx < $nerr; $err_idx++) { 
      my $err_code = $err_info_HAR->{"code"}[$err_idx];
      if(! exists $ftbl_err_exceptions_AHR->[$i]{$err_code}) { 
        DNAORG_FAIL("ERROR in $sub_name, array element $i does not have a hash key of $err_code", 1, $FH_HR);
      }
      if($ftbl_err_exceptions_AHR->[$i]{$err_code} ne "R" &&
         $ftbl_err_exceptions_AHR->[$i]{$err_code} ne "A" &&
         $ftbl_err_exceptions_AHR->[$i]{$err_code} ne "D") { 
        DNAORG_FAIL("ERROR in $sub_name, array element $i does not have valid value of \"R\" \"A\" or \"D\", but rather $ftbl_err_exceptions_AHR->[$i]{$err_code}", 1, $FH_HR);
      }
    }
  }

  # make sure all other required keys exist as well
  for(my $i = 0; $i < $nexc; $i++) { 
    for(my $j = 0; $j < scalar(@other_reqd_keys_A); $j++) { 
      if(! exists $ftbl_err_exceptions_AHR->[$i]{$other_reqd_keys_A[$j]}) { 
        DNAORG_FAIL("ERROR in $sub_name, array element $i does not have a hash key of $other_reqd_keys_A[$j]", 1, $FH_HR);
      }
    }
  }

  return $nexc;
}
#################################################################
#################################################################
#
# Subroutines related to codons:
#   fetchStopCodon()
#   fetchStartCodon()
#   fetchCodon()
#   validateStopCodon()
#
#################################################################
# Subroutine: fetchStopCodon()
# Incept:     EPN, Mon Mar 14 13:34:02 2016
#
# Synopsis:   Fetch a stop codon given its final position,
#             and strand.
#
# Args:       $sqfile:   Bio::Easel::SqFile object, open sequence
#                        file containing $seqname;
#             $seq_name: name of sequence to fetch part of
#             $stop:     final position of the stop codon
#             $strand:   strand we want ("+" or "-")
#             $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#             
# Returns:    The stop codon as a string
# 
# Dies:       If $stop (or $stop - 2) is negative.
#################################################################
sub fetchStopCodon {
  my $sub_name = "fetchStopCodon";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seq_name, $stop, $strand, $FH_HR) = @_;

  my $stop_codon_posn;
  if($strand eq "-") { 
    $stop_codon_posn = $stop + 2;
  }
  else { 
    $stop_codon_posn = $stop - 2;
  }
  if($stop_codon_posn < 0) { 
    DNAORG_FAIL("ERROR in $sub_name(), trying to fetch stop codon for $seq_name positions $stop_codon_posn..$stop, but we expect positive positions", 1, $FH_HR);
  }
  # printf("in $sub_name, seqname $seqname, stop $stop\n");

  return fetchCodon($sqfile, $seq_name, $stop_codon_posn, $strand, $FH_HR);
}

#################################################################
# Subroutine: fetchStartCodon()
# Incept:     EPN, Tue Mar 15 10:18:34 2016
#
# Synopsis:   Fetch a start codon given its first position,
#             and strand.
#
# Args:       $sqfile:   Bio::Easel::SqFile object, open sequence
#                        file containing $seqname;
#             $seq_name: name of sequence to fetch part of
#             $stop:     final position of the stop codon
#             $strand:   strand we want ("+" or "-")
#             $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#             
# Returns:    The stop codon as a string
# 
# Dies:       If $start is negative.
#################################################################
sub fetchStartCodon {
  my $sub_name = "fetchStartCodon";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seq_name, $start, $strand, $FH_HR) = @_;

  if($start < 0) { 
    DNAORG_FAIL("ERROR in $sub_name(), trying to fetch start codon for $seq_name at position $start on strand $strand, but we expect positive positions", 1, $FH_HR);
  }

  return fetchCodon($sqfile, $seq_name, $start, $strand, $FH_HR);
}

#################################################################
# Subroutine: fetchCodon()
# Incept:     EPN, Mon Mar 14 13:37:31 2016
# 
# Synopsis:   Fetch a codon given it's first position
#             and the strand and a Bio::Easel::SqFile object
#             that is the open sequence file with the desired
#             sequence.
#
# Args:       $sqfile:   Bio::Easel::SqFile object, open sequence
#                        file containing $seqname;
#             $seq_name: name of sequence to fetch part of
#             $start:    start position of the codon
#             $strand:   strand we want ("+" or "-")
#             $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    The codon as a string
#
#################################################################
sub fetchCodon {
  my $sub_name = "fetchCodon";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seqname, $start, $strand, $FH_HR) = @_;

  my $codon_start = $start;
  my $codon_stop  = ($strand eq "-") ? $start - 2 : $start + 2; 
  my $newname = $seqname . "/" . $codon_start . "-" . $codon_stop;

  my @fetch_AA = ();
  # if $seqname does not exist in $sqfile, the following line
  # will cause the script to fail, ungracefully
  my $seqlen = $sqfile->fetch_seq_length_given_name($seqname);
  if(($codon_start < 1)       || ($codon_stop < 1) || 
     ($codon_start > $seqlen) || ($codon_stop > $seqlen)) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, trying to fetch codon with coordinates %d..%d for sequence $seqname, valid coordinates are %d..%d\n", 
                        $codon_start, $codon_stop, 1, $seqlen), 1, $FH_HR);
  }

  push(@fetch_AA, [$newname, $codon_start, $codon_stop, $seqname]);

  my $faseq = $sqfile->fetch_subseqs(\@fetch_AA, -1);

  my ($header, $seq) = split("\n", $faseq);

  # printf("\nin $sub_name, $seqname $start $strand returning $seq\n");
  
  return $seq;
}

#################################################################
# Subroutine: validateStopCodon()
# Incept:     EPN, Mon Mar 14 13:47:57 2016
# 
# Purpose:    Given a codon, return '1' if it's a valid stop codon,
#             else return 0.
#
# Args:
#  $codon:  the codon
#
# Returns:    The codon as a string
#
#################################################################
sub validateStopCodon {
  my $sub_name = "validateStopCodon";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($codon) = @_;
  
  $codon =~ tr/a-z/A-Z/;
  $codon =~ tr/U/T/;

  if($codon eq "TAA" || 
     $codon eq "TGA" || 
     $codon eq "TAG" || 
     $codon eq "TAR") { 
    return 1;
  }

  return 0;
}

#################################################################
#
# Subroutines related to timings:
#   secondsSinceEpoch()
#   formatTimeString()
#
#################################################################
# Subroutine : secondsSinceEpoch()
# Incept:      EPN, Sat Feb 13 06:17:03 2016
#
# Purpose:     Return the seconds and microseconds since the 
#              Unix epoch (Jan 1, 1970) using 
#              Time::HiRes::gettimeofday().
#
# Arguments:   NONE
# 
# Returns:     Number of seconds and microseconds
#              since the epoch.
#
################################################################# 
sub secondsSinceEpoch { 
  my $nargs_expected = 0;
  my $sub_name = "secondsSinceEpoch()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seconds, $microseconds) = gettimeofday();
  return ($seconds + ($microseconds / 1000000.));
}

#####################################################################
# Subroutine: formatTimeString()
# Incept:     EPN, Fri Oct 24 13:18:23 2014 [rnavore]
# 
# Purpose:    Get a timing in hhhh:mm:ss format.
# 
# Arguments:
# $inseconds: number of seconds
#
# Returns:    string that describes time in hhhh:mm:ss format
# 
# Dies:       Never.
#
####################################################################
sub formatTimeString { 
  my $nargs_expected = 1;
  my $sub_name = "formatTimeString()";
   if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($inseconds) = @_;
  my ($i, $hours, $minutes, $seconds, $thours, $tminutes, $tseconds, $ndig_hours);

  $hours = int($inseconds / 3600);
  $inseconds -= ($hours * 3600);
  $minutes = int($inseconds / 60);
  $inseconds -= ($minutes * 60);
  $seconds = $inseconds;
  $thours   = sprintf("%02d", $hours);
  $tminutes = sprintf("%02d", $minutes);
  $ndig_hours = length($hours);
  if($ndig_hours < 2) { $ndig_hours = 2; }

  $tseconds = sprintf("%05.2f", $seconds);
  my $ret_str = sprintf("%*s:%2s:%5s", $ndig_hours, $thours, $tminutes, $tseconds);
  # %*s covers two of the arguments: $ndig_hours specifies width of string, $thours is the string
  
  return $ret_str;
}

#################################################################
#################################################################
#
# Simple utility subroutines for hashes and arrays:
#   findNonNumericValueInArray()
#   numNonNumericValueInArray()
#   maxLengthScalarKeyInHash()
#   maxLengthScalarValueInHash()
#   maxLengthScalarValueInArray()
#   findValueInArray()
#   sumArray()
#   sumHashValues()

#
#################################################################
# Subroutine : findNonNumericValueInArray()
# Incept:      EPN, Tue Feb 16 10:40:57 2016
#
# Purpose:     Returns (first) index in @{$AR} that has the 
#              nonnumeric value $value. Returns -1 
#              if it does not exist.
#
# Arguments: 
#   $AR:       REF to array 
#   $value:    the value we're checking exists in @{$AR}
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns:     index ($i) '1' if $value exists in @{$AR}, '-1' if not
#
# Dies:        if $value is numeric, or @{$AR} is not defined.
################################################################# 
sub findNonNumericValueInArray { 
  my $nargs_expected = 3;
  my $sub_name = "findNonNumericValueInArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR, $value, $FH_HR) = @_;

  if(verify_real($value)) { 
    DNAORG_FAIL("ERROR in $sub_name, value $value seems to be numeric, we can't compare it for equality", 1, $FH_HR);
  }

  if(! defined $AR) { 
    DNAORG_FAIL("ERROR in $sub_name, array reference is not defined", 1, $FH_HR);
  }

  for(my $i = 0; $i < scalar(@{$AR}); $i++) {
    if($AR->[$i] eq $value) { 
      return $i; 
    }
  }

  return -1; # did not find it
}

#################################################################
# Subroutine : numNonNumericValueInArray()
# Incept:      EPN, Fri Mar 11 06:34:51 2016
#
# Purpose:     Returns number of times nonnumeric value 
#              $value exists in @{$AR}. Returns 0 if
#              it doesn't exist.
#
# Arguments: 
#   $AR:       REF to array 
#   $value:    the value we're looking for in @{$AR}
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns:     Number of occurrences of $value in @{$AR}.
#
# Dies:        if $value is numeric, or @{$AR} is not defined.
################################################################# 
sub numNonNumericValueInArray { 
  my $nargs_expected = 3;
  my $sub_name = "numNonNumericValueInArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR, $value, $FH_HR) = @_;

  if(verify_real($value)) { 
    DNAORG_FAIL("ERROR in $sub_name, value $value seems to be numeric, we can't compare it for equality", 1, $FH_HR);
  }

  if(! defined $AR) { 
    DNAORG_FAIL("ERROR in $sub_name, array reference is not defined", 1, $FH_HR);
  }

  my $ct = 0;
  for(my $i = 0; $i < scalar(@{$AR}); $i++) {
    if($AR->[$i] eq $value) { 
      $ct++;
    }
  }

  return $ct;
}

#################################################################
# Subroutine : maxLengthScalarKeyInHash()
# Incept:      EPN, Thu Dec 13 15:52:09 2018
# 
# Purpose:     Return the maximum length of a scalar key
#              in a hash.
#
# Arguments: 
#   $HR: reference to the hash
# 
# Returns:     The length of the maximum length scalar key.
#
################################################################# 
sub maxLengthScalarKeyInHash { 
  my $nargs_expected = 1;
  my $sub_name = "maxLengthScalarKeyInHash()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = $_[0];

  my $max = 0;
  my $len = 0;
  foreach my $key (keys (%{$HR})) { 
    $len = length($key);
    if($len > $max) { $max = $len; }
  }
  return $max;
}

#################################################################
# Subroutine : maxLengthScalarValueInHash()
# Incept:      EPN, Mon Nov  3 09:09:59 2014 [rnavore]
# 
# Purpose:     Return the maximum length of a scalar value
#              in a hash.
#
# Arguments: 
#   $HR: reference to the hash
# 
# Returns:     The length of the maximum length scalar.
#
################################################################# 
sub maxLengthScalarValueInHash { 
  my $nargs_expected = 1;
  my $sub_name = "maxLengthScalarValueInHash()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = $_[0];

  my $max = 0;
  my $len = 0;
  foreach my $key (keys (%{$HR})) { 
    $len = length($HR->{$key});
    if($len > $max) { $max = $len; }
  }
  return $max;
}

#################################################################
# Subroutine : maxLengthScalarValueInArray()
# Incept:      EPN, Thu Mar 17 12:38:53 2016
# 
# Purpose:     Return the maximum length of a scalar value
#              in an array.
#
# Arguments: 
#   $AR: reference to the array
# 
# Returns:     The length of the maximum length scalar.
#
################################################################# 
sub maxLengthScalarValueInArray { 
  my $nargs_expected = 1;
  my $sub_name = "maxLengthScalarValueInArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR) = $_[0];

  my $max = 0;
  my $len = 0;
  foreach my $el (@{$AR}) { 
    $len = length($el);
    if($len > $max) { $max = $len; }
  }
  return $max;
}

#################################################################
# Subroutine: findValueInArray()
# Incept:     EPN, Tue Mar  8 11:26:03 2016
# Synopsis:   Look for a value in an array and return the index
#             of it, if found, else return -1. If it exists more than
#             once, return the minimum index.
#
# Arguments:
#  $value:   value to look for
#  $AR:      array to look in
#
# Returns:    void
# 
#################################################################
sub findValueInArray { 
  my $sub_name = "findValueInArray()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($value, $AR, $FH_HR) = @_;

  if(verify_real($value)) { # compare with ==
    for(my $i = 0; $i < scalar(@{$AR}); $i++) { 
      my $el = $AR->[$i];
      if(verify_real($el) && ($value == $el)) { 
        return $i;
      }
    }
  }
  else { # compare with 'eq'
    for(my $i = 0; $i < scalar(@{$AR}); $i++) { 
      my $el = $AR->[$i];
      if((! verify_real($el)) && ($value eq $el)) { 
        return $i;
      }
    }
  }
  return -1;
}  

#################################################################
# Subroutine : sumArray()
# Incept:      EPN, Wed Feb  7 13:53:07 2018
# 
# Purpose:     Sum the scalar values in an array
#
# Arguments: 
#   $AR: reference to the array
# 
# Returns:     Sum of all values in the array
#
################################################################# 
sub sumArray {
  my $nargs_expected = 1;
  my $sub_name = "sumArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR) = $_[0];

  my $sum = 0;
  foreach my $el (@{$AR}) { 
    $sum += $el; 
  }
  return $sum;
}

#################################################################
# Subroutine : sumHashValues()
# Incept:      EPN, Wed Feb  7 13:58:25 2018
# 
# Purpose:     Sum the values for all keys in a hash
#
# Arguments: 
#   $HR: reference to the hash
# 
# Returns:     Sum of all values in the hash
#
################################################################# 
sub sumHashValues {
  my $nargs_expected = 1;
  my $sub_name = "sumHashValues()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = $_[0];

  my $sum = 0;
  foreach my $key (keys (%{$HR})) { 
    $sum += $HR->{$key};
  }
  return $sum;
}

#################################################################
#################################################################
#
# Simple utility subroutines:
#   DNAORG_FAIL()
#   fileOpenFailure()
#   runCommand()
#   removeDirPath()
#   removeScriptNameFromString()
#   removeFileUsingSystemRm()
#   getMonocharacterString()
#   countLinesInFile()
#   validateFileExistsAndIsNonEmpty()
#   concatenateListOfFiles()
#   md5ChecksumOfFile()
#   nseBreakdown()
#
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
  if(scalar(@_) != $nargs_expected) { 
    if(scalar(@_) > 0) { 
      printf STDERR ("ERROR, DNAORG_FAIL() entered with %d != %d input arguments.\n(errmsg: $_[0])\n\n", scalar(@_), $nargs_expected); 
    }
    else { 
      printf STDERR ("ERROR, DNAORG_FAIL() entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); 
    }
    exit(1); 
  }
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
      if(defined $FH_HR->{$key}) { 
        close $FH_HR->{$key};
      }
    }
  }
  
  printf STDERR $errmsg; 
  exit($status);
}

#################################################################
# Subroutine : fileOpenFailure()
# Incept:      EPN, Wed Nov 11 05:39:56 2009 (rnavore)
#
# Purpose:     Called if an open() call fails on a file.
#              Print an informative error message
#              to $FH_HR->{"cmd"} and $FH_HR->{"log"}
#              and to STDERR, then exit with <$status>.
#
# Arguments: 
#   $filename:   file that we couldn't open
#   $c_sub_name: name of calling subroutine name
#   $status:     error status
#   $action:     "reading", "writing", "appending"
#   $FH_HR:      ref to hash of open file handles to close
# 
# Returns:     Nothing, this function will exit the program.
#
################################################################# 
sub fileOpenFailure { 
  my $nargs_expected = 5;
  my $sub_name = "fileOpenFailure()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $c_sub_name, $status, $action, $FH_HR) = @_;

  if(($action eq "reading") && (! (-e $filename))) { 
    DNAORG_FAIL(sprintf("ERROR, could not open %s%s for reading. It does not exist.", $filename, (defined $c_sub_name) ? " in subroutine $c_sub_name" : ""), $status, $FH_HR);
  }
  else { 
    DNAORG_FAIL(sprintf("ERROR, could not open %s%s for %s", $filename, (defined $c_sub_name) ? " in subroutine $c_sub_name" : "", $action), $status, $FH_HR);
  }

  return; # never reached
}

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
#   $FH_HR:       REF to hash of file handles, including "cmd"
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
# Subroutine : removeDirPath()
# Incept:      EPN, Mon Nov  9 14:30:59 2009 [ssu-align]
#
# Purpose:     Given a full path of a file remove the directory path.
#              For example: "foodir/foodir2/foo.stk" becomes "foo.stk".
#
# Arguments: 
#   $fullpath: name of original file
# 
# Returns:     The string $fullpath with dir path removed.
#
################################################################# 
sub removeDirPath {
  my $sub_name = "removeDirPath()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my $fullpath = $_[0];

  $fullpath =~ s/^.+\///;

  return $fullpath;
}

#################################################################
# Subroutine : removeScriptNameFromString()
# Incept:      EPN, Wed Mar  2 14:04:14 2016
#
# Purpose:     Remove a 'dnaorg_*' script name from a string and
#              return the modified version. 
#              For example: "NC_001346.dnaorg_annotate' becomes 'NC_001346'
#
# Arguments: 
#   $string: string to potentially modify
# 
# Returns:     The string $string with any script names removed.
#
################################################################# 
sub removeScriptNameFromString { 
  my $sub_name = "removeScriptNameFromString()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my $string = $_[0];

  my @script_A = ("dnaorg_build", "dnaorg_annotate");

  foreach my $script (@script_A) { 
    $string =~ s/\.$script\./\./g; # if bordered by two '.', leave one of them
    $string =~ s/\.$script$//g;    # if this ends the string and a '.' precedes it, remove ".$string"
    $string =~ s/$script//g;       # if it's a substr, but not caught by the two commands above, remove it
  }

  return $string;
}

#################################################################
# Subroutine: removeFileUsingSystemRm
# Incept:     EPN, Fri Mar  4 15:57:25 2016
#
# Purpose:    Remove a file from the filesystem by using
#             the system rm command.
# Arguments:
#   $file:            file to remove
#   $caller_sub_name: name of caller, can be undef
#   $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    - if the file does not exist
#
#################################################################
sub removeFileUsingSystemRm { 
  my $sub_name = "removeFileUsingSystemRm";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($file, $caller_sub_name, $opt_HHR, $FH_HR) = (@_);
  
  if(! -e $file) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, %s trying to remove file $file but it does not exist", 
                (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), 1, $FH_HR); 
  }

  runCommand("rm $file", opt_Get("-v", $opt_HHR), $FH_HR);

  return;
}

#################################################################
# Subroutine: getMonocharacterString()
# Incept:     EPN, Thu Mar 10 21:02:35 2016
#
# Purpose:    Return a string of length $len of repeated instances
#             of the character $char.
#
# Arguments:
#   $len:   desired length of the string to return
#   $char:  desired character
#   $FH_HR: REF to hash of file handles, including "log" and "cmd"
#
# Returns:  A string of $char repeated $len times.
# 
# Dies:     if $len is not a positive integer
#
#################################################################
sub getMonocharacterString {
  my $sub_name = "getMonocharacterString";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($len, $char, $FH_HR) = @_;

  if(! verify_integer($len)) { 
    DNAORG_FAIL("ERROR in $sub_name, passed in length ($len) is not a non-negative integer", 1, $FH_HR);
  }
  if($len < 0) { 
    DNAORG_FAIL("ERROR in $sub_name, passed in length ($len) is a negative integer", 1, $FH_HR);
  }
    
  my $ret_str = "";
  for(my $i = 0; $i < $len; $i++) { 
    $ret_str .= $char;
  }

  return $ret_str;
}

#################################################################
# Subroutine : countLinesInFile()
# Incept:      EPN, Tue Mar  1 09:36:56 2016
#
# Purpose:     Count the number of lines in a file
#              by opening it and reading in all lines.
#
# Arguments: 
#   $filename:         file that we are checking on
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If $filename does not exist or cannot be opened for reading.
#
################################################################# 
sub countLinesInFile { 
  my $nargs_expected = 2;
  my $sub_name = "countLinesInFile()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $FH_HR) = @_;

  my $nlines = 0;
  open(IN, $filename) || fileOpenFailure($filename, $sub_name, $!, "reading", $FH_HR);
  while(<IN>) { 
    $nlines++;
  }
  close(IN);

  return $nlines;
}

#################################################################
# Subroutine : fileLinesToArray()
# Incept:      EPN, Tue Nov 21 10:26:58 2017
#
# Purpose:     Store each non-blank line in a file as an element
#              in an array, after removing newline.
#
# Arguments: 
#   $filename:                   file that we are parsing
#   $remove_trailing_whitespace: '1' to remove trailing whitespace in each line, '0' not to
#   $AR:                         ref to array to add to
#   $FH_HR:                      ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If $filename does not exist or cannot be opened for reading.
#
################################################################# 
sub fileLinesToArray { 
  my $nargs_expected = 4;
  my $sub_name = "fileLinesToArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $remove_trailing_whitespace, $AR, $FH_HR) = @_;

  open(IN, $filename) || fileOpenFailure($filename, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    if($line =~ /\S/) { 
      chomp $line;
      if($remove_trailing_whitespace) { $line =~ s/\s*$//; }
      push(@{$AR}, $line);
    }
  }
  close(IN);

  return;
}

#################################################################
# Subroutine : arrayToNewlineDelimitedString()
# Incept:      EPN, Fri Dec 14 09:21:41 2018
#
# Purpose:     Return a newline delimited string with all values of an array.
#
# Arguments: 
#   $AR:     ref to array
# 
# Returns:     string
# 
# Dies:        Never.
#
################################################################# 
sub arrayToNewlineDelimitedString {
  my $nargs_expected = 1;
  my $sub_name = "arrayToNewlineDelimitedString";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR) = @_;

  my $retstr = "";
  foreach my $el (@{$AR}) { 
    $retstr .= $el . "\n";
  }
  if($retstr eq "") { 
    $retstr = "\n";
  }
  return $retstr;
}

#################################################################
# Subroutine : hashKeysToNewlineDelimitedString()
# Incept:      EPN, Fri Dec 14 09:25:24 2018
#
# Purpose:     Return a newline delimited string with all (sorted) keys in a hash.
#
# Arguments: 
#   $HR:     ref to hash
# 
# Returns:     string
# 
# Dies:        Never.
#
################################################################# 
sub hashKeysToNewlineDelimitedString {
  my $nargs_expected = 1;
  my $sub_name = "hashKeysToNewlineDelimitedString";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = @_;

  my $retstr = "";
  foreach my $el (sort keys %{$HR}) { 
    $retstr .= $el . "\n";
  }
  if($retstr eq "") { 
    $retstr = "\n";
  }
  return $retstr;
}

#################################################################
# Subroutine : hashValuesToNewlineDelimitedString()
# Incept:      EPN, Fri Dec 14 09:26:16 2018
#
# Purpose:     Return a newline delimited string with all values in a hash.
#
# Arguments: 
#   $HR:     ref to hash
# 
# Returns:     string
# 
# Dies:        Never.
#
################################################################# 
sub hashValuesToNewlineDelimitedString {
  my $nargs_expected = 1;
  my $sub_name = "hashValuesToNewlineDelimitedString";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($HR) = @_;

  my $retstr = "";
  foreach my $el (sort keys %{$HR}) { 
    $retstr .= $HR->{$el} . "\n";
  }
  if($retstr eq "") { 
    $retstr = "\n";
  }
  return $retstr;
}


#################################################################
# Subroutine : validateFileExistsAndIsNonEmpty()
# Incept:      EPN, Mon Feb 29 16:16:07 2016
#
# Purpose:     Check if a file exists and is non-empty.
#              If it does not exist or it is empty,
#              die via DNAORG_FAIL().
#
# Arguments: 
#   $filename:         file that we are checking on
#   $calling_sub_name: name of calling subroutine (can be undef)
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If $filename does not exist or is empty.
#
################################################################# 
sub validateFileExistsAndIsNonEmpty { 
  my $nargs_expected = 3;
  my $sub_name = "validateFileExistsAndIsNonEmpty()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $caller_sub_name, $FH_HR) = @_;

  if(! -e $filename) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, %sfile $filename does not exist.", (defined $caller_sub_name ? "called by $caller_sub_name," : "")), 1, $FH_HR);
  }
  elsif(! -s $filename) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, %sfile $filename exists but is empty.", (defined $caller_sub_name ? "called by $caller_sub_name," : "")), 1, $FH_HR);
  }
  
  return;
}

#################################################################
# Subroutine : concatenateListOfFiles()
# Incept:      EPN, Sun Apr 24 08:08:15 2016
#
# Purpose:     Concatenate a list of files into one file.
#              If the list has more than 500 files, split
#              up job into concatenating 500 at a time.
# 
#              We remove all files that we concatenate unless
#              --keep option is on in %{$opt_HHR}.
#
# Arguments: 
#   $file_AR:          REF to array of all files to concatenate
#   $outfile:          name of output file to create by concatenating
#                      all files in @{$file_AR}.
#   $caller_sub_name:  name of calling subroutine (can be undef)
#   $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If one of the cat commands fails.
#              If $outfile is in @{$file_AR}
#              If @{$file_AR} contains more than 800*800 files
#              (640K) if so, we may need to call this function
#              recursively twice (that is, recursive call will
#              also call itself recursively) and we don't have 
#              a sophisticated enough temporary file naming
#              strategy to handle that robustly.
################################################################# 
sub concatenateListOfFiles { 
  my $nargs_expected = 5;
  my $sub_name = "concatenateListOfFiles()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($file_AR, $outfile, $caller_sub_name, $opt_HHR, $FH_HR) = @_;

  if(findNonNumericValueInArray($file_AR, $outfile, $FH_HR) != -1) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name%s, output file name $outfile exists in list of files to concatenate", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), 1, $FH_HR);
  }

  # first, convert @{$file_AR} array into a 2D array of file names, each of which has 
  # a max of 800 elements, we'll concatenate each of these lists separately
  my $max_nfiles = 800;
  my $nfiles = scalar(@{$file_AR});

  if($nfiles > ($max_nfiles * $max_nfiles)) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name%s, trying to concatenate %d files, our limit is %d", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : "", $nfiles, $max_nfiles * $max_nfiles), 
                1, $FH_HR);
  }
    
  my ($idx1, $idx2); # indices in @{$file_AR}, and of secondary files
  my @file_AA = ();
  $idx2 = -1; # get's incremented to 0 in first loop iteration
  for($idx1 = 0; $idx1 < $nfiles; $idx1++) { 
    if(($idx1 % $max_nfiles) == 0) { 
      $idx2++; 
      @{$file_AA[$idx2]} = (); # initialize
    }
    push(@{$file_AA[$idx2]}, $file_AR->[$idx1]);
  }
  
  my $nconcat = scalar(@file_AA);
  my @tmp_outfile_A = (); # fill this with names of temporary files we create
  my $tmp_outfile; # name of an output file we'll create
  for($idx2 = 0; $idx2 < $nconcat; $idx2++) { 
    if($nconcat == 1) { # special case, we don't need to create any temporary files
      $tmp_outfile = $outfile;
    }
    else { 
      $tmp_outfile = $outfile . ".tmp" . ($idx2+1); 
      # make sure this file does not exist in @{$file_AA[$idx2]} to avoid klobbering
      # if it does, continue to append .tmp($idx2+1) until it doesn't
      while(findNonNumericValueInArray($file_AA[$idx2], $tmp_outfile, $FH_HR) != -1) { 
        $tmp_outfile .= ".tmp" . ($idx2+1); 
      }
    }
    # create the concatenate command
    my $cat_cmd = "cat ";
    foreach my $tmp_file (@{$file_AA[$idx2]}) {
      $cat_cmd .= $tmp_file . " ";
    }
    $cat_cmd .= "> $tmp_outfile";

    # execute the command
    runCommand($cat_cmd, opt_Get("-v", $opt_HHR), $FH_HR);

    # add it to the array of temporary files
    push(@tmp_outfile_A, $tmp_outfile); 
  }

  if(scalar(@tmp_outfile_A) > 1) { 
    # we created more than one temporary output file, concatenate them
    # by calling this function again
    concatenateListOfFiles(\@tmp_outfile_A, $outfile, (defined $caller_sub_name) ? $caller_sub_name . ":" . $sub_name : $sub_name, $opt_HHR, $FH_HR);
  }

  if(! opt_Get("--keep", $opt_HHR)) { 
    # remove all of the original files, be careful to not remove @tmp_outfile_A
    # because the recursive call will handle that
    foreach my $file_to_remove (@{$file_AR}) { 
      removeFileUsingSystemRm($file_to_remove, 
                              (defined $caller_sub_name) ? $caller_sub_name . ":" . $sub_name : $sub_name, 
                              $opt_HHR, $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine : md5ChecksumOfFile()
# Incept:      EPN, Fri May 27 14:02:30 2016
#
# Purpose:     Use md5sum to get a checksum of a file, return
#              the checksum. Not efficient. Creates a temporary
#              file and then deletes it.
# 
# Arguments: 
#   $file:             REF to array of all files to concatenate
#   $caller_sub_name:  name of calling subroutine (can be undef)
#   $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     md5sum of the file.
# 
# Dies:        If the file doesn't exist or the command fails.
#
################################################################# 
sub md5ChecksumOfFile { 
  my $nargs_expected = 4;
  my $sub_name = "md5ChecksumOfFile()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($file, $caller_sub_name, $opt_HHR, $FH_HR) = @_;

  if(! -s $file) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name%s, file to get md5 checksum of ($file) does no exist or is empty", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), 1, $FH_HR);
  }

  my $out_file = removeDirPath($file . ".md5sum");
  runCommand("md5sum $file > $out_file", opt_Get("-v", $opt_HHR), $FH_HR);

  open(MD5, $out_file) || fileOpenFailure($out_file, $sub_name, $!, "reading", $FH_HR);
  #194625f7c3e2a5129f9880c7e29f63de  wnv.lin2.matpept.in
  my $md5sum = <MD5>;
  chomp $md5sum;
  if($md5sum =~ /^(\S+)\s+(\S+)$/) { 
    $md5sum = $1;
  }
  else { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name%s, unable to parse md5sum output: $md5sum", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), 1, $FH_HR);
  }
  close(MD5);

  removeFileUsingSystemRm($out_file, $caller_sub_name, $opt_HHR, $FH_HR);

  return $md5sum;
}

#################################################################
# Subroutine : nseBreakdown()
# Incept:      EPN, Wed Jan 30 09:50:07 2013 [rfam-family-pipeline:Utils.pm]
#
# Purpose  : Checks if $nse is of format "name/start-end" and if so
#          : breaks it down into $n, $s, $e, $str (see 'Returns' section)
# 
# Arguments: 
#   $seqname:  sequence name, possibly in "name/start-end" format
# 
# Returns:     5 values:
#              '1' if seqname was of "name/start-end" format, else '0'
#              $n:   name ("" if seqname does not match "name/start-end")
#              $s:   start, maybe <= or > than $e (0 if seqname does not match "name/start-end")
#              $e:   end,   maybe <= or > than $s (0 if seqname does not match "name/start-end")
#              $str: strand, 1 if $s <= $e, else -1
# 
# Dies:        Never
#
################################################################# 
sub nseBreakdown {
  my $nargs_expected = 1;
  my $sub_name = "nseBreakdown()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($sqname) = $_[0];

  my $n;       # sqacc
  my $s;       # start, from seq name (can be > $end)
  my $e;       # end,   from seq name (can be < $start)
  my $str;     # strand, 1 if $start <= $end, else -1

  if($sqname =~ m/^(\S+)\/(\d+)\-(\d+)\s*/) {
    ($n, $s, $e) = ($1,$2,$3);
    $str = ($s <= $e) ? 1 : -1; 
    return (1, $n, $s, $e, $str);
  }
  return (0, "", 0, 0, 0); # if we get here, $sqname is not in name/start-end format
}

#################################################################
#
# Miscellaneous subroutines that don't fall into one of the above
# categories:
#   stripVersion()
#   fetchedNameToListName()
#   fetchSequencesUsingEslFetchCds()
#   addNameAndBlankSsToStockholmAlignment()
#   getQualifierValues()
#   createCmDb()
#   matpeptValidateCdsRelationships()
#   checkForSpanningSequenceSegments()
#   annotateAppendFeatures()
#   getIndexHashForArray()
#
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
# Subroutine: fetchedNameToListName()
# Incept:     EPN, Thu Feb 11 14:25:52 2016
#
# Purpose:    Convert a fetched sequence name via efetch to the 
#             name that was used as input to efetch, with a version.
#
# Arguments: 
#   $fetched_name: name from fasta file efetch returned
#
# Returns:    $list_name, processed fetched_name
#
# Dies:       never
#################################################################
sub fetchedNameToListName { 
  my $sub_name  = "fetchedNameToListName";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($fetched_name) = (@_);

  my $list_name = $fetched_name;
  
  $list_name =~ s/^gi\|?\d+\|\w+\|//; # get rid of everything before the accession number
  if($list_name =~ /^pdb\|(\S+)\|(\S+)$/) {  # special case of PDB accessions
    $list_name = $1 . "_" . $2;
  }
  else { 
    $list_name =~ s/\|$//;              # get rid of end | or
    $list_name =~ s/\|\w+//;            # get rid of end | and everything after |
  }
  $list_name =~ s/\.[0-9]*$//;        # strip version

  return $list_name;
}

#################################################################
# Subroutine: fetchSequencesUsingEslFetchCds()
# Incept:     EPN, Thu Feb 11 15:30:07 2016
# 
# Purpose:    Fetch sequences by calling the esl_fetch_cds program
#             with the path $esl_fetch_cds. 
#
# Arguments:
#   $esl_fetch_cds:    path to esl-fetch-cds executable
#   $out_fetch_file:   name of file to use as input to $esl_fetch_cds, created here
#   $out_fasta_file:   name of fasta file to create, created here
#   $do_circular:      '1' to duplicate genome, '0' not to
#   $seq_info_HAR:     REF to 2D hash with sequence information, ADDED TO HERE
#                      "seq_name", and "seq_len" arrays filled here
#   $FH_HR:            REF to hash of file handles, including "log" and "cmd"
# 
# Returns:    void
#          
# Dies: if $esl_fetch_cds executable does not exist
#       if $esl_fetch_cds command fails
#################################################################
sub fetchSequencesUsingEslFetchCds { 
  my $sub_name = "fetchSequencesUsingEslFetchCds";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($esl_fetch_cds, $fetch_file, $fasta_file, $do_circular, $seq_info_HAR, $FH_HR) = @_;

  # contract check
  if(! -e $esl_fetch_cds) { 
    DNAORG_FAIL("ERROR in $sub_name, required executable $esl_fetch_cds does not exist", 1, $FH_HR); 
  }

  my ($accn_name, $accn_len, $seq_name, $seq_len); # accession name/length (from GenBank) sequence name/length in the sequence file we create

  @{$seq_info_HAR->{"seq_name"}} = (); # initialize
  @{$seq_info_HAR->{"seq_len"}}  = (); # initialize
  my $seq_name_AR  = \@{$seq_info_HAR->{"seq_name"}};
  my $seq_len_AR = \@{$seq_info_HAR->{"seq_len"}};

  my $fetch_string; # string output to the input file for esl-fetch-cds.pl
  
  # create the esl-fetch-cds input file
  open(OUT, ">", $fetch_file) || fileOpenFailure($fetch_file, $sub_name, $!, "writing", $FH_HR);
  my $naccn = scalar(@{$seq_info_HAR->{"accn_name"}});
  for(my $a = 0; $a < $naccn; $a++) { 
    my $accn_name = $seq_info_HAR->{"accn_name"}[$a];
    my $accn_len  = $seq_info_HAR->{"accn_len"}[$a];

    if($do_circular) { 
      $fetch_string = "join(" . $accn_name . ":1.." . $accn_len . "," . $accn_name . ":1.." . $accn_len . ")\n";
      print OUT $accn_name . ":" . "dnaorg-duplicated" . "\t" . $fetch_string;
      $seq_name = $accn_name . ":dnaorg-duplicated:" . $accn_name . ":1:" . $accn_len . ":+:" . $accn_name . ":1:" . $accn_len . ":+:";
      $seq_len  = 2 * $accn_len;
    }
    else { # $do_circular is FALSE
      $fetch_string = $accn_name . ":1.." . $accn_len . "\n";
      print OUT $accn_name . ":" . "dnaorg" . "\t" . $fetch_string;
      $seq_name = $accn_name . ":dnaorg:" . $accn_name . ":1:" . $accn_len . ":+:";
      $seq_len  = $accn_len;
    }
    push(@{$seq_info_HAR->{"seq_name"}}, $seq_name);
    push(@{$seq_info_HAR->{"seq_len"}}, $seq_len);
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
    DNAORG_FAIL("ERROR in $sub_name, file $in_file does not exist.", 1, $FH_HR); 
  }

  # open and validate file
  my $msa = Bio::Easel::MSA->new({
    fileLocation => $in_file,
    isDna => 1
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
# Returns:    number of values added; fills @{$values_AR}
# 
# Dies: if $accn doesn't exist in %{$tbl_HHAR}
#################################################################
sub getQualifierValues {
  my $sub_name = "getQualifierValues()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($tbl_HHAR, $accn, $qualifier, $values_AR, $FH_HR) = @_;

  if(! exists $tbl_HHAR->{$accn}) { DNAORG_FAIL("ERROR in $sub_name, no data for accession: $accn", 1, $FH_HR); }

  if((! defined $tbl_HHAR->{$accn}{$qualifier}) ||  
     (scalar(@{$tbl_HHAR->{$accn}{$qualifier}}) == 0)) { 
    return 0; # no annotation for $qualifier, do not update array
  }

  my $nvalues = scalar(@{$tbl_HHAR->{$accn}{$qualifier}});

  if ($nvalues > 0) { 
    for(my $i = 0; $i < $nvalues; $i++) { 
      push(@{$values_AR}, $tbl_HHAR->{$accn}{$qualifier}[$i]);
    }
  }

  return $nvalues;
}


#################################################################
# Subroutine : pressCmDb()
# Incept:      EPN, Mon Feb 29 14:26:52 2016
#
# Purpose:     Run cmpress on a CM database file.
#
# Arguments: 
#   $model_file:            the full path to the concatenated model file to create
#   $cmpress:               path to cmpress executable
#   $do_add_to_output_info: '1' to add file info to output info, '0' not to
#   $do_run:                '1' to actually run the command, 0 to only return it
#   $opt_HHR:               REF to 2D hash of option values, see top of epn-options.pm for description
#   $ofile_info_HHR:        REF to the 2D hash of output file information
# 
# Returns:     cmpress command
# 
# Dies: If any of the expected individual CM files do not exist.
#
################################################################# 
sub pressCmDb {
  my $sub_name = "pressCmDb()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($model_file, $cmpress, $do_add_to_output_info, $do_run, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $cmpress_cmd = "$cmpress -F $model_file > /dev/null"; # output is irrelevant
  if($do_run) {
    runCommand($cmpress_cmd, opt_Get("-v", $opt_HHR), $FH_HR);

    if($do_add_to_output_info) { 
      addClosedFileToOutputInfo($ofile_info_HHR, "$model_file.i1m", $model_file.".i1m", 0, "index file for the CM, created by cmpress");
      addClosedFileToOutputInfo($ofile_info_HHR, "$model_file.i1i", $model_file.".i1i", 0, "index file for the CM, created by cmpress");
      addClosedFileToOutputInfo($ofile_info_HHR, "$model_file.i1f", $model_file.".i1f", 0, "index file for the CM, created by cmpress");
      addClosedFileToOutputInfo($ofile_info_HHR, "$model_file.i1p", $model_file.".i1p", 0, "index file for the CM, created by cmpress");
    }
  }
  return $cmpress_cmd;
}

#################################################################
# Subroutine : concatenateIndividualCmFiles()
# Incept:      EPN, Mon Feb 29 10:53:53 2016
#
# Purpose:     Concatenate individual CM files created and calibrated
#              by a previous call to dnaorg_build.pl into one 
#              CM file.
#
# Arguments: 
#   $combined_model_file:   the full path to the concatenated model file to create
#   $out_root:              output root the individual CM files share
#   $nmdl:                  number of models
#   $do_run:                '1' to actually run the command, 0 to only return it
#   $opt_HHR:               REF to 2D hash of option values, see top of epn-options.pm for description
#   $ofile_info_HHR:        REF to the 2D hash of output file information
# 
# Returns:     the cat command
# 
# Dies: If any of the expected individual CM files do not exist.
#
################################################################# 
sub concatenateIndividualCmFiles {
  my $sub_name = "concatenateIndividualCmFiles()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($combined_model_file, $out_root, $nmdl, $do_run, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # concatenate the files into a CM DB
  my $cat_cmd = "cat";
  for(my $i = 0; $i < $nmdl; $i++) { 
    my $indi_model = $out_root . ".$i.cm"; # dnaorg_build created the CMs
    if(($do_run) && (! -s ($indi_model))) { 
      DNAORG_FAIL("ERROR, individual model file $indi_model does not exist.", 1, $FH_HR);
    }
    $cat_cmd .= " $indi_model ";
  }
  $cat_cmd .= " > $combined_model_file";
  if($do_run) { 
    runCommand($cat_cmd, opt_Get("-v", $opt_HHR), $FH_HR);
    addClosedFileToOutputInfo($ofile_info_HHR, "cm", $combined_model_file, 1, "CM file (a concatenation of individual files created by dnaorg_build.pl)");
  }

  return $cat_cmd;
}

#################################################################
# Subroutine: matpeptValidateCdsRelationships()
# Incept:     EPN, Fri Feb 12 09:15:58 2016
#
# Purpose:    Validate the CDS:primary-mat_peptide relationships 
#             in @{$cds2pmatpept_AAR}, which were probably 
#             ready from a file in parsed in parseMatPeptSpecFile().
#             These 'primary' relationships are valid if
#             the mat_peptides that comprise the CDS
#             do in fact completely comprise it (when 
#             concatenated the mat_peptide coordinates completely
#             cover the CDS coordinates, in order, with no
#             overlaps). For example a CDS that spans 1..300
#             is comprised by three mature peptides with 
#             coordinates 1..150, 151..223, and 224..300.
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
  
  # validate cds2pmatpept_AA
  my $ref_ncds = scalar(@ref_cds_len_A);
  my $ref_nmp  = scalar(@ref_mp_len_A);
  my $ncds2check = scalar(@{$cds2pmatpept_AAR});
  my $prv_stop = undef; # previous mat_peptide's stop position, 
  for(my $cds_idx = 0; $cds_idx < scalar(@{$cds2pmatpept_AAR}); $cds_idx++) { 
    if(($cds_idx < 0) || ($cds_idx >= $ref_ncds)) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, cds_idx $cds_idx is out of range", $cds_idx+1), 1, $FH_HR); 
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
#          DNAORG_FAIL("ERROR in $sub_name, multiple exon CDS with an intron broken up to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
        }
      }
      else { # negative strand
        if(($cds_stops_A[0] - $cds_starts_A[1] - 1) > 0) { 
#          DNAORG_FAIL("ERROR in $sub_name, multiple exon CDS with an intron broken up to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
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
# Incept:     EPN, Mon Apr 11 14:36:45 2016
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

  # if $$nexons_R is not '2', we skip this block and return '0'.
  # 'spanning' exons can only occur if there are exactly 2 exons that
  # we are checking. 
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
    # remember if we get here, we know we only have 2 exons, i.e scalar(@{$starts_AR}) and scalar(@{$stops_AR}) is 2
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
# Subroutine : getIndexHashForArray()
# Incept:      EPN, Tue Mar  1 14:11:38 2016
#
# Purpose:     Create an 'index hash' %{$index_HR} for an array
#              @{$AR}, such that $index_HR{$value} = $n, if 
#              $AR->[$n] = $value. 
#              
# Arguments: 
#   $AR:       REF to the array, PRE-FILLED
#   $index_HR: REF to the hash to fill, FILLED HERE
#   $FH_HR:    REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If there are any duplicate values in array as
#              measured by eq, or if any elements in the array
#              are numeric (as determined by verify_real() from
#              epn-options.pm)
#
################################################################# 
sub getIndexHashForArray { 
  my $nargs_expected = 3;
  my $sub_name = "getIndexHashForArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR, $index_HR, $FH_HR) = @_;

  # initialize
  %{$index_HR} = (); 

  for(my $i = 0; $i < scalar(@{$AR}); $i++) { 
    my $el = $AR->[$i];
    # verify it's not a number, we don't want numbers because
    # then we can't be sure that $index_HR->{$num} will be 
    # testable using 'eq' due to precision issues with storing
    # numbers. 
    if(verify_real($el)) { 
      DNAORG_FAIL("ERROR in $sub_name(), value $el is numeric"); 
    }
    if(exists $index_HR->{$el}) { 
      DNAORG_FAIL("ERROR in $sub_name(), the value $el appears twice in the array"); 
    }
    $index_HR->{$el} = $i;
  }

  return;
}

#################################################################
# Subroutine : waitForFarmJobsToFinish()
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
#          A job is considered 'finished in error' if it outputs
#          anything to its err file in @{$errfile_AR}. (We do not kill
#          the job, although for the jobs we are monitoring with this
#          subroutine, it should already have died (though we don't
#          guarantee that in anyway).) If any jobs 'finish in error'
#          this subroutine will continue until all jobs have finished
#          or we've waited $nmin minutes and then it will cause the
#          program to exit in error and output an error message
#          listing the jobs that have 'finished in error'
#          
#          When $do_errcheck is 1, this function considers any output
#          written to stderr output files in @{$errfile_AR} to mean
#          that the corresponding job has 'failed' and should be
#          considered to have finished. When $do_errchecks is 0
#          we don't look at the err files.
# 
#
# Arguments: 
#  $outfile_AR:      ref to array of output files that will be created by jobs we are waiting for
#  $errfile_AR:      ref to array of err files that will be created by jobs we are waiting for if 
#                    any stderr output is created
#  $finished_str:    string that indicates a job is finished e.g. "[ok]"
#  $nmin:            number of minutes to wait
#  $do_errcheck:     '1' to consider output to an error file a 'failure' of a job, '0' not to.
#  $FH_HR:           REF to hash of file handles
#
# Returns:     Number of jobs (<= scalar(@{$outfile_AR})) that have
#              finished.
# 
# Dies: never.
#
################################################################# 
sub waitForFarmJobsToFinish { 
  my $sub_name = "waitForFarmJobsToFinish()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($outfile_AR, $errfile_AR, $finished_str, $nmin, $do_errcheck, $FH_HR) = @_;

  my $log_FH = $FH_HR->{"log"};

  my $njobs = scalar(@{$outfile_AR});
  if($njobs != scalar(@{$errfile_AR})) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, number of elements in outfile array ($njobs) differ from number of jobs in errfile array (%d)", scalar(@{$errfile_AR})), 1, $FH_HR);
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
  outputString($log_FH, 1, "\n");
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
        if(-s $outfile_AR->[$i]) { 
          my $final_line = `tail -n 1 $outfile_AR->[$i]`;
          chomp $final_line;
          if($final_line =~ m/\r$/) { chop $final_line; } # remove ^M if it exists
          if($final_line =~ m/\Q$finished_str\E/) { 
            $is_finished_A[$i] = 1;
            $nfinished++;
          }
        }
        if(($do_errcheck) && (-s $errfile_AR->[$i])) { # errfile exists and is non-empty, this is a failure, even if we saw $finished_str above
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
    outputString($log_FH, 1, sprintf("#\t%4d of %4d jobs finished (%.1f minutes spent waiting)\n", $nfinished, $njobs, $secs_waited / 60.));

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
        $errmsg .= "\t$outfile_AR->[$i]\t$errfile_AR->[$i]\n";
      }
    }
    DNAORG_FAIL($errmsg, 1, $FH_HR);
  }

  # if we get here we have no failures
  return $nfinished;
}

#################################################################
# Subroutine : splitFastaFile()
# Incept:      EPN, Tue Mar  1 09:30:10 2016
#
# Purpose: Split up a fasta file into <n> smaller files by calling
#          the esl-ssplit perl script.
#
# Arguments: 
#  $esl_ssplit:      path to the esl-ssplit.pl script to use
#  $fasta_file:      fasta file to split up
#  $nfiles:          desired number of files to split $fasta_file into
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:    Number of files actually created (can differ from requested
#             amount (which is $nfiles)).
#
# Dies:       if esl-ssplit command fails
#
################################################################# 
sub splitFastaFile { 
  my $sub_name = "splitFastaFile()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_ssplit, $fasta_file, $nfiles, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $outfile = $fasta_file . ".esl-ssplit";
  my $cmd = "$esl_ssplit -v -r -n $fasta_file $nfiles > $outfile";
  runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);

  # parse output to determine exactly how many files were created:
  # $esl_ssplit will have output exactly 1 line per fasta file it created
  my $nfiles_created = countLinesInFile($outfile, $FH_HR);

  if(! opt_Get("--keep", $opt_HHR)) { 
    runCommand("rm $outfile", opt_Get("-v", $opt_HHR), $FH_HR);
  }

  return $nfiles_created;
}

#################################################################
# Subroutine:  cmscanOrNhmmscanWrapper()
# Incept:      EPN, Tue Nov 21 15:54:20 2017
#
# Purpose:     Run one or more cmscan or nhmmscan jobs on the farm
#              or locally, after possibly splitting up the input
#              sequence file. 
#              The following must all be valid options in opt_HHR:
#              --nseq, --local, --wait, --keep, -v
#              See dnaorg_annotate.pl for examples of these options.
#
# Arguments: 
#  $execs_HR:        ref to executables with "esl-ssplit" and either "cmscan" or "nhmmscan" 
#                    defined as keys
#  $do_cmscan:       '1' if we are running cmscan, '0' if we are running nhmmscan
#  $out_root:        string for naming output files
#  $seq_file:        name of sequence file with all sequences to run against
#  $tot_len_nt:      total length of all nucleotides in $seq_file
#  $tblout_file:     string for name of concatnated tblout file to create
#  $progress_w:      width for outputProgressPrior output
#  $mdl_filename_AR: ref to array of model file names
#  $mdl_len_AR:      ref to array of model lengths, can be undef, but only if $do_cmscan == 0
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
#
# Returns:     void, updates $$nfa_created_R with number of
#              fasta files created.
# 
# Dies: If an executable doesn't exist, or cmscan or nhmmscan or esl-ssplit
#       command fails if we're running locally
################################################################# 
sub cmscanOrNhmmscanWrapper { 
  my $sub_name = "cmscanOrNhmmscanWrapper";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $do_cmscan, $out_root, $seq_file, $tot_len_nt, $tblout_file, $progress_w,
      $mdl_filename_AR, $mdl_len_AR, $opt_HHR, $ofile_info_HHR) = @_;

  # contract check
  if($do_cmscan) { 
    if(! defined $mdl_len_AR) { 
      DNAORG_FAIL("ERROR in $sub_name, do_cmscan is TRUE but mdl_len_AR is undefined.", 1, $ofile_info_HHR->{"FH"});
    }
  }

  my $program_choice = ($do_cmscan) ? "cmscan" : "nhmmscan";
  my @tmp_seq_file_A = ();    # array of sequence files created by esl-ssplit, we'll remove after we're done unless --keep
  my @tmp_tblout_file_A = (); # array of tblout files, we'll remove after we're done, unless --keep
  my @tmp_stdout_file_A = (); # array of stdout files created by cmscan/nhmmscan, only created if -v
  my @tmp_err_file_A = ();    # array of error files we'll remove after we're done, unless --keep (empty if we run locally)
  my $nfasta_created = 0;     # number of fasta files created by esl-ssplit
  my $log_FH = $ofile_info_HHR->{"FH"}{"log"}; # for convenience
  my $mdl_filename = undef; # file name of current model file
  my $mdl_len      = undef; # length of current model
  # filter threshold settings relevant only if $do_cmscan is TRUE
  my $do_max = 0; # if $do_cmscan, possibly set to true based on model length
  my $do_mid = 0; # if $do_cmscan, possibly set to true based on model length
  my $do_big = 0; # if $do_cmscan, possibly set to true based on model length
  my $start_secs; # timing start
  my $nmdl = scalar(@{$mdl_filename_AR});

  # determine how many sequence files we need to split $seq_file into to satisfy <n> sequences
  # per job (where n is from --nseq <n>). If more than 1 and --local not used, we split up 
  # the sequence file and submit jobs to farm, if only 1 or --local used, we just run it locally
  my $targ_nseqfiles = int($tot_len_nt / (opt_Get("--nkb", $opt_HHR) * 1000)); 
  # int() takes the floor, so there can be a nonzero remainder. We don't add 1 though, 
  # because splitFastaFile() will return the actual number of sequence files created
  # and we'll use that as the number of jobs subsequently. $nfarmjobs is currently only
  # the 'target' number of sequence files that we pass into splitFastaFile().
  # make sure we won't exceed our max number of jobs (from --maxnjobs)
  if(($targ_nseqfiles * $nmdl) > (opt_Get("--maxnjobs", $opt_HHR))) { 
    $targ_nseqfiles = int((opt_Get("--maxnjobs", $opt_HHR)) / $nmdl);
  }
  if($targ_nseqfiles == 0) { $targ_nseqfiles = 1; }

  # if($targ_nseqfiles == 1) || (opt_Get("--local", $opt_HHR))) {  # to run single jobs locally, used to be default 
  if(opt_Get("--local", $opt_HHR)) {  # to run single jobs locally, used to be default 
    # run jobs locally
    $start_secs = outputProgressPrior("Running $program_choice locally", $progress_w, $log_FH, *STDOUT);
    # run a different cmscan/nhmmscan run for each model file
    for(my $m = 0; $m < $nmdl; $m++) { 
      my $tmp_tblout_file = $out_root . ".m" . ($m+1) . ".tblout";
      my $tmp_stdout_file = opt_Get("-v", $opt_HHR) ? $out_root . ".m" . ($m+1) . ".cmscan" : "/dev/null";
      $mdl_filename = $mdl_filename_AR->[$m];
      $mdl_len      = (defined $mdl_len_AR) ? $mdl_len_AR->[$m] : 0;
      if($do_cmscan) { 
        ($do_max, $do_mid, $do_big) = determineCmscanFilterSettings($mdl_len, $opt_HHR);
      }
      runCmscanOrNhmmscan($execs_HR->{"$program_choice"}, $do_cmscan, 1, $do_max, $do_mid, $do_big, $mdl_filename,
                          $seq_file, $tmp_stdout_file, $tmp_tblout_file, $opt_HHR, $ofile_info_HHR); # 1: run locally
      push(@tmp_tblout_file_A, $tmp_tblout_file);
      push(@tmp_err_file_A,    $tmp_tblout_file . ".err"); # this will be the name of the error output file, set in run_cmscan
      if($tmp_stdout_file ne "/dev/null") { push(@tmp_stdout_file_A, $tmp_stdout_file); }
    }
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
  else { 
    # we need to split up the sequence file, and submit a separate set of nhmmscan/cmscan jobs (one per model file) for each
    my $nfasta_created = 0;
    if($targ_nseqfiles > 1) { 
      $nfasta_created = splitFastaFile($execs_HR->{"esl-ssplit"}, $seq_file, $targ_nseqfiles, $opt_HHR, $ofile_info_HHR);
      # splitFastaFile will return the actual number of fasta files created, 
      # which can differ from the requested amount (which is $targ_nseqfiles) that we pass in. 
    }
    else { 
      my $cp_command = "cp " . $seq_file . " " . $seq_file . ".1";
      runCommand($cp_command, opt_Get("-v", $opt_HHR), $ofile_info_HHR->{"FH"});
      $nfasta_created = 1;
    }
    my $nfarmjobs = $nfasta_created * $nmdl;
    # now submit a job for each
    $start_secs = outputProgressPrior("Submitting $nfarmjobs $program_choice jobs to the farm", $progress_w, $log_FH, *STDOUT);
    for(my $s = 1; $s <= $nfasta_created; $s++) { 
      my $tmp_seq_file = $seq_file . "." . $s;
      push(@tmp_seq_file_A, $tmp_seq_file);
      for(my $m = 0; $m < $nmdl; $m++) { 
        my $tmp_tblout_file = $out_root . ".m" . ($m+1) . ".s" . $s . ".tblout";
        my $tmp_stdout_file = opt_Get("-v", $opt_HHR) ? $out_root . ".m" . ($m+1) . ".s" . $s . "." . $program_choice : "/dev/null";
        $mdl_filename = $mdl_filename_AR->[$m];
        $mdl_len      = (defined $mdl_len_AR) ? $mdl_len_AR->[$m] : 0;
        if($do_cmscan) { 
          ($do_max, $do_mid, $do_big) = determineCmscanFilterSettings($mdl_len, $opt_HHR);
        }
        runCmscanOrNhmmscan($execs_HR->{"$program_choice"}, $do_cmscan, 0, $do_max, $do_mid, $do_big, $mdl_filename,
                            $tmp_seq_file, $tmp_stdout_file, $tmp_tblout_file, $opt_HHR, $ofile_info_HHR);   # 0: do not run locally
        push(@tmp_tblout_file_A, $tmp_tblout_file);
        push(@tmp_err_file_A,    $tmp_tblout_file . ".err"); # this will be the name of the error output file, set in run_cmscan
        if($tmp_stdout_file ne "/dev/null") { push(@tmp_stdout_file_A, $tmp_stdout_file); }
      }
    }
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
    
    # wait for the jobs to finish
    $start_secs = outputProgressPrior(sprintf("Waiting a maximum of %d minutes for all farm jobs to finish", opt_Get("--wait", $opt_HHR)), 
                                      $progress_w, $log_FH, *STDOUT);
    my $njobs_finished = waitForFarmJobsToFinish(\@tmp_tblout_file_A, \@tmp_err_file_A, "[ok]", opt_Get("--wait", $opt_HHR), opt_Get("--errcheck", $opt_HHR), $ofile_info_HHR->{"FH"});
    if($njobs_finished != $nfarmjobs) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name only $njobs_finished of the $nfarmjobs are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", $opt_HHR)), 1, $ofile_info_HHR->{"FH"});
    }
    if(! opt_Get("--local", $opt_HHR)) { 
      outputString($log_FH, 1, "# "); # necessary because waitForFarmJobsToFinish() creates lines that summarize wait time and so we need a '#' before 'done' printed by outputProgressComplete()
    }
  } # end of 'else' entered if($nfarmjobs > 1 && ! --local)
    
  # concatenate all the tblout files into one 
  concatenateListOfFiles(\@tmp_tblout_file_A, $tblout_file, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
  # remove temporary files if --keep not enabled
  if(! opt_Get("--keep", $opt_HHR)) { 
    foreach my $tmp_file (@tmp_seq_file_A, @tmp_tblout_file_A, @tmp_err_file_A) { 
      if(-e $tmp_file) { 
        removeFileUsingSystemRm($tmp_file, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
        # we don't remove stdout files, they'll only exist if -v, so we assume user wants them
      }
    }
  }

  return;
}

#################################################################
# Subroutine:  determineCmscanFilterSettings()
# Incept:      EPN, Wed Apr  5 10:52:38 2017
#
# Purpose:     Given a model length and command line options, 
#              determine the filtering mode to run cmscan in.
#
# Arguments: 
#  $model_len:       length of model
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#
# Returns:     3 values, exactly 0 or 1 of which will be '1'
#              others will be '0'
#              $do_max: '1' to run in max sensitivity mode
#              $do_mid: '1' to run in mid sensitivity mode
#              $do_big: '1' to run in fast (min sensitivity) mode
# 
# Dies: Never
################################################################# 
sub determineCmscanFilterSettings { 
  my $sub_name = "determineCmscanFilterSettings()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($model_len, $opt_HHR) = @_;

  my $do_max = 0;
  my $do_mid = 0;
  my $do_big = 0;

  if($model_len <= (opt_Get("--smallthresh", $opt_HHR))) { 
    $do_max = 1; # set filter thresholds to --max for small models
  } 
  if(! $do_max) { 
    if($model_len <= (opt_Get("--midthresh", $opt_HHR))) { 
      $do_mid = 1; # set filter thresholds to --mid for kind of small models
    } 
  }
  if((opt_IsUsed("--hmmonly", $opt_HHR)) && (! $do_max) && (! $do_mid)) {
    if($model_len >= (opt_Get("--hmmonly", $opt_HHR))) { 
      $do_big = 1; # use HMM mode for big models
    }
  }

  return ($do_max, $do_mid, $do_big);
}

#################################################################
# Subroutine : runCmscanOrNhmmscan()
# Incept:      EPN, Mon Feb 29 15:09:22 2016
#
# Purpose:     Run Infernal's cmscan executable or HMMER's nhmmscan 
#              executable using $model_file as the model file on sequence 
#              file $seq_file.
#
# Arguments: 
#  $executable:      path to the cmscan or nhmmscan executable file
#  $do_cmscan:       '1' if we are running cmscan, '0' if we are running nhmmscan
#  $do_local:        '1' to run locally, '0' to submit job to farm
#  $do_max:          '1' to run with --max option, '0' to run with default options
#                    '1' is usually used only when $model_file contains a single 
#                    very short model
#  $do_mid:          '1' to run with --mid option, '0' to not
#  $do_big:          '1' to run with --mxsize 6144 option, '0' to run with default options
#  $model_file:      path to the CM file
#  $seq_file:        path to the sequence file
#  $stdout_file:     path to the stdout file to create, can be "/dev/null", or undef 
#  $tblout_file:     path to the cmscan --tblout file to create
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:     void
# 
# Dies: If at least one CM was not built from the current reference.
#       If $do_max and $do_big are both true.
################################################################# 
sub runCmscanOrNhmmscan { 
  my $sub_name = "runCmscanOrNhmmscan()";
  my $nargs_expected = 12;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($executable, $do_cmscan, $do_local, $do_max, $do_mid, $do_big, $model_file, $seq_file, $stdout_file, $tblout_file, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # contract checks
  if($do_big && $do_max) { 
    DNAORG_FAIL("ERROR in $sub_name, do_big and do_max are both true, only one should be.", 1, $FH_HR);
  }
  if($do_big && $do_mid) { 
    DNAORG_FAIL("ERROR in $sub_name, do_big and do_mid are both true, only one should be.", 1, $FH_HR);
  }
  if($do_max && $do_mid) { 
    DNAORG_FAIL("ERROR in $sub_name, do_max and do_mid are both true, only one should be.", 1, $FH_HR);
  }
  validateFileExistsAndIsNonEmpty($model_file, $sub_name, $FH_HR); 
  validateFileExistsAndIsNonEmpty($seq_file,   $sub_name, $FH_HR);

  my $program_choice = ($do_cmscan) ? "cmscan" : "nhmmscan";

  my $opts = (opt_Get("-v", $opt_HHR)) ? " " : " --noali ";
  $opts .= " --cpu 0 --tblout $tblout_file";

  # cmscan specific options
  if($program_choice eq "cmscan") { 
    $opts .= " --verbose";
  
    if($do_max) { # no filtering
      $opts .= " --max -E 0.1 "; # with --max, a lot more FPs get through the filter, so we enforce an E-value cutoff,
                               # but we need to keep it high so very short models achieve it, 
                               # e.g. NC_002549 (Ebola) has a length 10 model that gets a 0.068 E-value in the reference (self-hit)
#      $opts .= " --max -E 0.01 "; # with --max, a lot more FPs get through the filter, so we enforce an E-value cutoff
    }
    elsif($do_mid) { 
#      $opts .= " --mid -E 0.1 --noF6 --olonepass --tau 0.001 --cyk --acyk -g --mxsize 1028."; # with --mid, more FPs get through the filter, so we enforce an E-value cutoff
      $opts .= " --mid -E 0.1 --noF6 --olonepass --cyk --acyk -g --mxsize 1028."; # with --mid, more FPs get through the filter, so we enforce an E-value cutoff
    }
    else { 
#      $opts .= " --F1 0.02 --F2 0.001 --F2b 0.001 --F3 0.00001 --F3b 0.00001 --F4 0.0002 --F4b 0.0002 --F5 0.0002 --noF6 --olonepass --tau 0.001 --cyk --acyk -g --mxsize 1028. ";
      $opts .= " --F1 0.02 --F2 0.001 --F2b 0.001 --F3 0.00001 --F3b 0.00001 --F4 0.0002 --F4b 0.0002 --F5 0.0002 --noF6 --olonepass --cyk --acyk -g --mxsize 1028. ";
    }
    # finally add --nohmmonly if we're not a big model
    if(! $do_big) { # do not use hmm unless model is big
      $opts .= " --nohmmonly ";
    }
  }

  my $cmd = "$executable $opts $model_file $seq_file > $stdout_file";

  # remove the tblout file if it exists, this is important because we'll use the existence and
  # final line of this file to determine when the jobs are finished, if it already exists, we'll
  # think the job is finished before it actual is.
  if(-e $tblout_file) { removeFileUsingSystemRm($tblout_file, $sub_name, $opt_HHR, $ofile_info_HHR); }

  # run cmscan, either locally or by submitting jobs to the farm
  if($do_local) { 
    # run locally
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);
  }
  else { 
    # submit job to farm and return
    my $jobname = "s" . removeDirPath($seq_file);
    my $errfile = $tblout_file . ".err";
    if(-e $errfile) { removeFileUsingSystemRm($errfile, $sub_name, $opt_HHR, $ofile_info_HHR); }
    my $farm_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n -l h_rt=288000,h_vmem=16G,mem_free=16G,reserve_mem=16G,m_mem_free=16G " . "\"" . $cmd . "\" > /dev/null\n";
    runCommand($farm_cmd, opt_Get("-v", $opt_HHR), $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: featureInfoKeyToFeatureTableQualifierName()
# Incept:     EPN, Tue Dec  5 14:22:25 2017
#
# Purpose:    Given a key from the ftr_info_HA{"type"} array, 
#             convert it into a string for a qualifier name
#             in a feature table.
#
#             Input           Return value
#             'out_product'   "product"
#             'out_gene'      "gene"
#             'out_exception' "exception"
#
# Arguments:
#   $in_key:  input key from %ftr_info_HA
#   $FH_HR:   REF to hash of file handles, including "log" and "cmd"
#             
# Returns:    qualifier name as a string for the feature table
#
# Dies:       if $in_key is unrecognized
#################################################################
sub featureInfoKeyToFeatureTableQualifierName { 
  my $sub_name  = "featureInfoKeyToFeatureTableQualifierName";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_key, $FH_HR) = (@_);

  if($in_key eq "out_product") { 
    return "product";
  }
  elsif($in_key eq "out_gene") { 
    return "gene";
  }
  elsif($in_key eq "out_exception") { 
    return "exception";
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unrecognized input key string: $in_key.", 1, $FH_HR);
  }

  return ""; # NEVERREACHED
}

#################################################################
# Subroutine: featureTypeIsCds()
# Incept:     EPN, Tue Feb  6 10:47:58 2018
#
# Purpose:    Is a feature type a CDS? 
#
#             Input        Return value
#             'cds-notmp'  1
#             'cds-mp'     1
#             'mp':        0
#             other:       0
#
# Arguments:
#   $feature_type:  input feature
#             
# Returns:    1 or 0
#
#################################################################
sub featureTypeIsCds() { 
  my $sub_name  = "featureTypeIsCds";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_feature) = (@_);

  if($in_feature eq "cds-notmp") { 
    return 1;
  }
  elsif($in_feature eq "cds-mp") { 
    return 1;
  }
  else { 
    return 0;
  }

  return ""; # NEVERREACHED
}
#################################################################
# Subroutine: featureTypeIsMaturePeptide()
# Incept:     EPN, Tue Feb  6 11:50:25 2018
#
# Purpose:    Is a feature type a mature peptide? 
#
#             Input        Return value
#             'mp':        1
#             'cds-notmp'  0
#             'cds-mp'     0
#             other:       0
#
# Arguments:
#   $feature_type:  input feature
#             
# Returns:    1 or 0
#
#################################################################
sub featureTypeIsMaturePeptide() { 
  my $sub_name  = "featureTypeIsMaturePeptide";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_feature) = (@_);

  if($in_feature eq "mp") { 
    return 1;
  }
  else { 
    return 0;
  }

  return ""; # NEVERREACHED
}
#################################################################
# Subroutine: featureTypeIsExtraFeature()
# Incept:     EPN, Tue Feb  6 12:17:29 2018
#
# Purpose:    Is a feature type an 'extra feature'?
#
#             Input        Return value
#             'xfeat':     1
#             'dfeat':     0
#             'mp':        0
#             'cds-notmp'  0
#             'cds-mp'     0
#             other:       0
#
# Arguments:
#   $feature_type:  input feature
#             
# Returns:    1 or 0
#
#################################################################
sub featureTypeIsExtraFeature() { 
  my $sub_name  = "featureTypeIsExtraFeature";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_feature) = (@_);

  if($in_feature eq "xfeat") { 
    return 1;
  }
  else { 
    return 0;
  }

  return ""; # NEVERREACHED
}

#################################################################
# Subroutine: featureTypeIsDuplicateFeature()
# Incept:     EPN, Sun Jul 22 18:38:57 2018
#
# Purpose:    Is a feature type a 'duplicate feature'?
#
#             Input        Return value
#             'xfeat':     0
#             'dfeat':     1
#             'mp':        0
#             'cds-notmp'  0
#             'cds-mp'     0
#             other:       0
#
# Arguments:
#   $feature_type:  input feature
#             
# Returns:    1 or 0
#
#################################################################
sub featureTypeIsDuplicateFeature() { 
  my $sub_name  = "featureTypeIsDuplicateFeature";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_feature) = (@_);

  if($in_feature eq "dfeat") { 
    return 1;
  }
  else { 
    return 0;
  }

  return ""; # NEVERREACHED
}

#################################################################
# Subroutine: getNumExtraOrDuplicateFeatures()
# Incept:     EPN, Sun Jul 22 13:10:33 2018
#
# Purpose:    Given a reference to a 3D hash of arrays, return
#             the number of elements in each of the arrays. 
#             (This will be the same number for all such arrays.)
#
# Arguments:
#   $tbl_HHHAR: 3D hash of arrays
#   $FH_HR:     ref to hash of file handles, including 'log'
#             
# Returns:    Number of elements in all arrays
#
# Dies: If not all arrays have the same number of elements.
#
#################################################################
sub getNumExtraOrDuplicateFeatures() { 
  my $sub_name  = "getNumExtraOrDuplicateFeatures()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($tbl_HHHAR, $FH_HR) = (@_);

  my $ret_size      = undef;
  my $ret_size_key1 = undef;
  my $ret_size_key2 = undef;
  my $ret_size_key3 = undef;
  my $size          = undef;
  foreach my $key1 (sort keys %{$tbl_HHHAR}) { 
    foreach my $key2 (sort keys %{$tbl_HHHAR->{$key1}}) { 
      foreach my $key3 (sort keys %{$tbl_HHHAR->{$key1}{$key2}}) {
        $size = scalar(@{$tbl_HHHAR->{$key1}{$key2}{$key3}});
        if(! defined $ret_size) {
          $ret_size = scalar(@{$tbl_HHHAR->{$key1}{$key2}{$key3}});
          $ret_size_key1 = $key1;
          $ret_size_key2 = $key2;
          $ret_size_key3 = $key3;
        }
        elsif($size != $ret_size) {
          DNAORG_FAIL("ERROR in $sub_name, not all arrays have the same number of elements.\ntbl_HHHAR->{$ret_size_key1}{$ret_size_key2}{$ret_size_key3} has $ret_size elements but tbl_HHHAR->{$key1}{$key2}{$key3} has $size elements", 1, $FH_HR);
        }
      }
    }
  }

  return $ret_size;
}

#################################################################
# Subroutine: formatTabDelimitedStringForErrorListFile()
# Incept:     EPN, Wed Dec 12 10:57:59 2018
#
# Purpose:    Given a sequence name and a string <error_str> that
#             describes an error, return a tab-delimited one that is
#             ready for output to an 'errors.list' file.
#
#             That return string will have 4 tokens:
#             <sequence-name>
#             <error-name>
#             <feature-name>
#             <error-description>
#
#             The input string must be in the following format:
#             <error-name>
# Arguments:
#   $seqname:   name of sequence
#   $errstr:    error string to convert
#   $FH_HR:     ref to hash of file handles, including 'log'
#             
# Returns:    $err_tab_str: tab delimited string in format described above.
#
# Dies: If $errstr is not in required format
#
#################################################################
sub formatTabDelimitedStringForErrorListFile() {
  my $sub_name  = "formatTabDelimitedStringForErrorListFile";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($seqname, $errstr, $FH_HR) = (@_);

  printf("HEYA 0 $seqname $errstr\n");

  # first replace any '_' characters with single spaces
  chomp $errstr;
  $errstr =~ s/\_/ /g;

  printf("HEYA 1 $seqname $errstr\n");

  my $error_name   = undef;
  my $feature_name = undef;
  my $error_desc   = undef;

  if($errstr =~ /^([^\:]+)\:\s*(.+)$/) {
    ($error_name, $error_desc) = ($1, $2);
    printf("HEYA 2 error_name: $error_name error_desc: $error_desc\n");
    if($error_desc =~ /^\([^\)]+\)\s*(.*)$/) { 
      $feature_name = $1;
      $error_desc   = $2;
    }
    if(! defined $feature_name) { 
      $feature_name = "*sequence*";
    }
  }
  elsif($errstr =~ /^([^\[]+)\[([^\]]+)\]\;$/) {
    # example from dnaorg_classify.pl:
    #Unexpected Classification[NC 001959,NC 029647 was specified, but NC 039476 is predicted];
    ($error_name, $error_desc) = ($1, $2);
    printf("HEYA 3 error_name: $error_name error_desc: $error_desc\n");
    $feature_name = "*sequence*";
  }
  elsif($errstr =~ /^([^\[\:]+)\;$/) {
    # example from dnaorg_classify.pl:
    #No Annotation;
    ($error_name) = ($1);
    printf("HEYA 4 error_name: $error_name error_desc: $error_desc\n");
    $feature_name = "*sequence*";
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unable to parse input error string: $errstr", 1, $FH_HR);
  }

  if($error_desc eq "") { 
    $error_desc = "-";
  }
  return $seqname . "\t" . $error_name . "\t" . $feature_name . "\t" . $error_desc;
}

###########################################################################
# the next line is critical, a perl module must return a true value
    return 1;
###########################################################################
