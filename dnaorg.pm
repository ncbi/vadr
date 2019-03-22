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
# - $FH_HR: a reference to a hash of file handles. Important keys are "log" and "cmd",
#           the log and command files we are outputting to. This data structure is passed
#           into nearly all functions because it is also passed into ofile_FAIL() which 
#           can be called from nearly all functions. ofile_FAIL() outputs an error message
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
# - $sgm_info_HAR: similar to ${ftr,seq,err}_info_HAR, except contains information pertaining 
#                  to each model, >= 1 of which will model a single feature (1 model for single
#                  exon CDS, 2 models for dual exon CDS, etc.). See 
#                  validateSegmentInfoHashIsComplete() for a list and explanation of the keys.
#                  The contents of this data structure can be printed to file <f> by dnaorg_build.pl 
#                  or dnaorg_annotate.pl with the --mdlinfo <f> command line option.
#
# - $seq_info_HAR: similar to ${ftr,mdl,err}_info_HAR, except contains information pertaining 
#                  to each sequence. See validateSequenceInfoHashIsComplete()
#                  for a list and explanation of the keys. The contents of this data structure
#                  can be printed to file <f> by dnaorg_annotate.pl with the --seqinfo <f>
#                  command line option.
#                   
# - $alt_info_HAR: similar to ${ftr,mdl}_info_HAR, except contains information pertaining to each 
#                  error code. See validateAlertInfoHashIsComplete() for a list and explanation 
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
#   getReferenceFeatureInfo()
#   fetchReferenceFeatureSequences()
#
# Subroutines related to the output info hash:
#   openAndAddFileToOutputInfo()
#   addClosedFileToOutputInfo()
#   helperAddFileToOutputInfo()
#
# Subroutines related to the error info hash:
#   dng_InitializeHardCodedAlertInfoHash()
#   addToAlertInfoHash()
#   setIncompatibilityAlertInfoHash()
#
# Subroutines related to the feature table error exception array of hashes:
#   initializeHardCodedFTableAlertExceptions()
#   addFTableAlertException()
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
#   ofile_OutputString()
#   outputBanner()
#   outputDividingLine()
# 
# Subroutines for dumping data structures, usually for debugging:
#   dumpInfoHashOfArrays()
#   dumpHashOfHashes()
#   dumpArrayOfHashesOfHashes()
#   dumpArrayOfHashes()
#
# Subroutines for validating the special data structures:
#   validateExecutableHash()
#   validateFeatureInfoHashIsComplete()
#   validateSegmentInfoHashIsComplete()
#   validateSequenceInfoHashIsComplete()
#   validateAlertInfoHashIsComplete()
#   validateInfoHashOfArraysIsComplete()
#   validateOutputFileInfoHashOfHashes()
#   validateAndGetSizeOfInfoHashOfArrays()
#   getConsistentSizeOfInfoHashOfArrays()
#   validateFTableAlertExceptions()
#
# Subroutines related to codons:
#   fetchStopCodon()
#   fetchStartCodon()
#   fetchCodon()
#   validateStopCodon()
#
# Subroutines related to timings:
#   ofile_SecondsSinceEpoch()
#   formatTimeString()
#
# Simple utility subroutines for hashes and arrays:
#   findNonNumericValueInArray()
#   numNonNumericValueInArray()
#   maxLengthScalarValueInHash()
#   maxLengthScalarValueInArray()
#   findValueInArray()
#   numberOfDigits()
# 
# Simple utility subroutines:
#   ofile_FAIL()
#   fileOpenFailure()
#   utl_RunCommand()
#   removeDirPath()
#   removeScriptNameFromString()
#   dng_RemoveFileUsingSystemRm()
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
#   dng_StripVersion()
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
#   determineFeatureTypes()
#   getNumFeaturesAnnotatedByModels()
#   getReferenceFeatureInfo()
#   fetchReferenceFeatureSequences()
#   featureHasChildren()
#   featureHasParent()
#


#################################################################
# Subroutine:  dng_FeaturesGetChildrenArrayOfArray()
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
sub dng_FeaturesGetChildrenArrayOfArrays { 
  my $nargs_expected = 3;
  my $sub_name = "dng_FeaturesGetChildrenArrayOfArrays";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $AAR, $FH_HR) = @_;

  @{$AAR} = ();
  my $nftr = scalar(@{$ftr_info_AHR});

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    dng_FeatureChildrenArray($ftr_info_AHR, $ftr_idx, $nftr, $AAR->[$ftr_idx], $FH_HR);
  }
  
  return;
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
# Subroutine:  dng_FeatureTypeAndTypeIndexString()
# Incept:      EPN, Sun Mar 10 06:41:53 2019
#
# Purpose:     Return a string giving feature type and type index
#              (e.g. "CDS#2") for feature $ftr_idx.
# 
# Arguments: 
#   $ftr_info_AHR:   REF to hash of arrays with information on the features, PRE-FILLED
#   $ftr_idx:        index we are interested in
#   $FH_HR:          REF to hash of file handles
# 
# Returns:     String
# 
#
################################################################# 
sub dng_FeatureTypeAndTypeIndexString { 
  my $nargs_expected = 3;
  my $sub_name = "dng_FeatureChildrenArray";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $ftr_idx, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  my $type = $ftr_info_AHR->[$ftr_idx]{"type"};
  my $type_idx = 1;
  for(my $ftr_idx2 = 0; $ftr_idx2 < $ftr_idx; $ftr_idx2++) { 
    if($ftr_info_AHR->[$ftr_idx2]{"type"} eq $type) { $type_idx++; }
  }
  
  return $type . "#" . $type_idx;
}


#################################################################
#################################################################
#
# Subroutines related to the error info hash:
#   dng_AlertInfoInitialize()
#   addToAlertInfoHash()
#   setIncompatibilityAlertInfoHash()
#   setFTableInvalidatedByAlertInfoHash()
#   processFeatureAlertsForFTable()
#   populateFTableNoteOrAlert()
#
#################################################################
# Subroutine: dng_AlertInfoInitialize()
# Incept:     EPN, Fri Mar  4 12:56:43 2016
#
# Purpose:    Set the initial values in an alert info hash,
#             using the hardcoded information in this
#             function.
#
# Arguments:
#   $alt_info_HAR:  REF to hash of arrays of alert information, FILLED HERE
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $alt_info_HAR already has keys upon entering this function
#
#################################################################
sub dng_AlertInfoInitialize { 
  my $sub_name = "dng_AlertInfoInitialize";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HAR, $FH_HR) = (@_);

  if(scalar (keys(%{$alt_info_HAR})) > 0) { 
    ofile_FAIL("ERROR in $sub_name, alert info hash of arrays already has at least one key", "dnaorg", 1, $FH_HR);
  }

  # add each alert code, this function will die if we try to add the same code twice, or if something is wrong 
  # with how we try to add it (args to dng_AlertInfoAdd don't pass the contract check)

  # classification errors
  dng_AlertInfoAdd($alt_info_HAR, "c_noa", "sequence",
                   "no significant similarity detected", # description
                   1, 1, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "No Annotation: (*sequence*) no significant similarity detected", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "c_usg", "sequence",
                   "unexpected classification", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Unexpected Classification: (*sequence*) !DESC!", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "c_ugr", "sequence",
                   "unexpected taxonomy", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Unexpected Taxonomy: (*sequence*) !DESC!", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "c_lcv", "sequence",
                   "low coverage", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Low Coverage: (*sequence*) !DESC!", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "c_lod", "sequence",
                   "low difference", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Low Difference: (*sequence*) !DESC!", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "c_vld", "sequence",
                   "very low difference", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Very Low Difference: (*sequence*) !DESC!", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "c_los", "sequence",
                   "low score", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Low Score: (*sequence*) !DESC!", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "c_vls", "sequence",
                   "very low score", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Very Low Score: (*sequence*) !DESC!", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "c_hbi", "sequence",
                   "high bias", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Biased Sequence: (*sequence*) !DESC!", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "c_mst", "sequence",
                   "minus strand", # description
                   1, 1, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Minus Strand: (*sequence*) !DESC!", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "n_nst", "feature",  
                   "no in-frame stop codon exists 3' of predicted valid start codon", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified; !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "n_nm3", "feature",  
                   "length of nucleotide feature is not a multiple of 3", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!; length is not a multiple of 3", # feature table note
                   "Unexpected Length: (!out_product,out_gene!) length is not a multiple of 3; !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "n_stp", "feature",
                   "predicted CDS stop by homology is invalid; there may be a valid stop in a different location due to truncation (trc) or extension (ext) (TAG|TAA|TGA)", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified on !out_product,out_gene!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "n_trc", "feature",
                   "in-frame stop codon exists 5' of stop position predicted by homology to reference", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!; contains premature stop codon", # feature table note
                   "CDS has Stop Codon: (!out_product,out_gene!) contains unexpected stop codon; !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "n_ext", "feature",
                   "first in-frame stop codon exists 3' of stop position predicted by homology to reference", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified; !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "b_per", "feature",
                   "mat_peptide may not be translated because its CDS has a problem", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note
                   "Peptide Translation Problem: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "n_str", "feature",
                   "predicted CDS start position is not beginning of start codon", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!; no start codon", # feature table note
                   "Mutation at Start: (!out_product,out_gene!) expected start codon could not be identified", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "n_gp5", "feature",
                   "alignment to reference is a gap at 5' boundary", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "n_gp3", "feature",
                   "alignment to reference is a gap at 3' boundary", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation at End: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "n_lp5", "feature",
                   "alignment to reference has low confidence at 5' boundary", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "n_lp3", "feature",
                   "alignment to reference has low confidence at 3' boundary", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation at End: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "b_non", "feature",
                   "protein-based search identifies CDS not identified in nucleotide-based search", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "b_nop", "feature",
                   "nucleotide-based search identifies CDS not identified in protein-based search", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "b_cst", "feature",
                   "strand mismatch between protein and nucleotide predictions", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "b_p5l", "feature",
                   "protein alignment extends past nucleotide alignment at 5' end", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "b_p5s", "feature",
                   "protein alignment does not extend close enough to nucleotide alignment 5' endpoint", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "b_p3l", "feature",
                   "protein alignment extends past nucleotide alignment at 3' end", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation at End: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "b_p3s", "feature",
                   "protein alignment does not extend close enough to nucleotide alignment 3' endpoint", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Indefinite Annotation at End: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "p_lin", "feature",
                   "too large of an insertion in protein alignment", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Insertion of Nucleotides: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "p_lde", "feature",
                   "too large of a deletion in protein alignment", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "Deletion of Nucleotides: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "p_trc", "feature",
                   "stop codon in protein alignment", # description
                   1, 0, # causes_failure, prevents_annot
                   "similar to !out_product,out_gene!", # feature table note
                   "CDS has Stop Codon: (!out_product,out_gene!) !DESC!", # feature table alert
                   $FH_HR);

  dng_AlertInfoAdd($alt_info_HAR, "b_zft", "sequence",
                   "zero features annotated", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "No Features Annotated: (*sequence*) zero annotated features", # feature table err
                   $FH_HR); 

  dng_AlertInfoAdd($alt_info_HAR, "n_div", "sequence",
                   "sequence too distant from reference to annotate", # description
                   1, 0, # causes_failure, prevents_annot
                   "", # feature table note, irrelevant for per-sequence errors
                   "Unexpected Divergence: (*sequence*) sequence is too divergent to confidently assign nucleotide-based annotation !DESC!", # feature table err
                   $FH_HR); 

  # define the ftbl_invalid_by values, these are one-sided, any error code listed in the 
  # 3rd argument invalidates the 2nd argument error code, but not vice versa

  # trc, ext and nst are preferred to stp
  dng_AlertInfoSetFTableInvalidatedBy($alt_info_HAR, "n_stp", "n_trc,n_ext,n_nst", $FH_HR); 

  # n_div is preferred to zft
  dng_AlertInfoSetFTableInvalidatedBy($alt_info_HAR, "b_zft", "n_div", $FH_HR);

  # validate the alert info hash
  dng_AlertInfoValidate($alt_info_HAR, undef, $FH_HR); 

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
#   $alt_info_HAR:    REF to array of hashes of error information, FILLED HERE
#   $code:            the code of the element we are adding
#   $pertype:         the 'per-type' of the element we are adding, "sequence" or "feature"
#   $desc:            the error description/message of the element we are adding
#   $causes_failure:  '1' if this alert causes its sequence to FAIL, '0' if not
#   $prevents_annot:  '1' if this alert prevents its sequence from being annotated, '0' if not
#   $ftbl_note:       note message for feature table
#                     must eq "" if per-type is "sequence"
#                     must ne "" if per-type is "feature"
#   $ftbl_alert:      alert message for feature table
#   $FH_HR:           REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if $alt_info_HAR->{"$code"} already exists
#          if $type ne "feature and ne "sequence"
#
######################p###########################################
sub dng_AlertInfoAdd { 
  my $sub_name = "dng_AlertInfoAdd";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HAR, $code, $pertype, $desc, $causes_failure, $prevents_annot, $fails, $ftbl_note, $ftbl_alert, $FH_HR) = (@_);

  # make sure $pertype is valid
  if(($pertype ne "feature") && ($pertype ne "sequence")) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code with per-type $pertype that is not neither \"feature\" nor \"sequence\".", "dnaorg", 1, $FH_HR); 
  }
  
  # make sure $ftbl_alert is valid
  if((! defined $ftbl_alert) || $ftbl_alert eq "") { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but ftbl_alert is undefined or empty", "dnaorg", 1, $FH_HR);
  }
  
  # make sure $causes_failure is valid
  if((! defined $causes_failure) || (($causes_failure != 0) && ($causes_failure != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but causes_failure is undefined or not 0 or 1", "dnaorg", 1, $FH_HR);
  }
  
  # make sure $prevents_annot is valid
  if((! defined $prevents_annot) || (($prevents_annot != 0) && ($prevents_annot != 1))) { 
    ofile_FAIL("ERROR in $sub_name, trying to add code $code but prevents_annot is undefined or not 0 or 1", "dnaorg", 1, $FH_HR);
  }
  
  # check if $code already exists
  if(exists $alt_info_HAR->{"code"}) { 
    my $nerr = scalar(@{$alt_info_HAR->{"code"}});
    for(my $alt_idx = 0; $alt_idx < $nerr; $alt_idx++) { 
      my $other_code = $alt_info_HAR->{"code"}[$alt_idx]; 
      if($code eq $other_code) { 
        ofile_FAIL(sprintf("ERROR in $sub_name, trying to add code $code, but it already exists as element %d in the error info hash", $alt_idx+1), "dnaorg", 1, $FH_HR);
      }
    }
  }
  
# initialize, if necessary
  if(! exists $alt_info_HAR->{"code"})            { @{$alt_info_HAR->{"code"}}            = (); }
  if(! exists $alt_info_HAR->{"pertype"})         { @{$alt_info_HAR->{"pertype"}}         = (); }
  if(! exists $alt_info_HAR->{"desc"})            { @{$alt_info_HAR->{"desc"}}            = (); }
  if(! exists $alt_info_HAR->{"causes_failure"})  { @{$alt_info_HAR->{"causes_failure"}}  = (); }
  if(! exists $alt_info_HAR->{"prevents_annot"})  { @{$alt_info_HAR->{"prevents_annot"}}  = (); }
  if(! exists $alt_info_HAR->{"ftbl_invalid_by"}) { @{$alt_info_HAR->{"ftbl_invalid_by"}} = (); }
  if(! exists $alt_info_HAR->{"ftbl_note"})       { @{$alt_info_HAR->{"ftbl_note"}}       = (); }
  if(! exists $alt_info_HAR->{"ftbl_alert"})      { @{$alt_info_HAR->{"ftbl_alert"}}      = (); }
  
  push(@{$alt_info_HAR->{"code"}},             $code); 
  push(@{$alt_info_HAR->{"pertype"}},          $pertype); 
  push(@{$alt_info_HAR->{"desc"}},             $desc); 
  push(@{$alt_info_HAR->{"causes_failure"}},   $causes_failure);
  push(@{$alt_info_HAR->{"prevents_annot"}},   $prevents_annot);
  push(@{$alt_info_HAR->{"ftbl_invalid_by"}},  ""); # initialized to no invalid_by's, possibly added to later with setFTableInvalidatedByErrorInfoHash()
  push(@{$alt_info_HAR->{"ftbl_note"}},        $ftbl_note);
  push(@{$alt_info_HAR->{"ftbl_alert"}},       $ftbl_alert);

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
#   $alt_info_HAR:  REF to hash of arrays of error information, FILLED HERE
#   $code1:         the code of the element we are adding ftbl_invalid_by values for
#   $code2str:      the codes $code1 is invalidated by, separated by commas
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies:    if one of the error codes in $code1 or $code2str do not
#          exist in %{$alt_info_HAR}.
#
#################################################################
sub dng_AlertInfoSetFTableInvalidatedBy {
  my $sub_name = "dng_AlertInfoSetFTableInvalidatedBy";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HAR, $code1, $code2str, $FH_HR) = (@_);

  my $idx1 = utl_ArrayFindNonNumericValue($alt_info_HAR->{"code"}, $code1, $FH_HR);
  if($idx1 == -1) { 
    ofile_FAIL("ERROR in $sub_name, trying to add ftbl_invalid_by for code $code1, but it does not exist in the error info hash", "dnaorg", 1, $FH_HR);
  }

  # verify the codes in $code2str
  my @code2_A = split(',', $code2str);
  foreach my $code2 (@code2_A) { 
    my $idx2 = utl_ArrayFindNonNumericValue($alt_info_HAR->{"code"}, $code2, $FH_HR);
    if($idx2 == -1) { 
      ofile_FAIL("ERROR in $sub_name, trying to add invalidated by relationship between codes $code1 and $code2, but $code2 does not exist in the error info hash", "dnaorg", 1, $FH_HR);
    }
    if($idx1 == $idx2) { 
      ofile_FAIL("ERROR in $sub_name, trying to add invalidated by relationship between a code and itself: $code1 and $code2", "dnaorg", 1, $FH_HR);
    }
  }

  # set the value
  $alt_info_HAR->{"ftbl_invalid_by"}[$idx1] = $code2str;

  return;
}

#################################################################
# Subroutine: dng_ProcessFeatureAlertsForFTable()
# Incept:     EPN, Thu Nov  1 12:10:34 2018
#
# Purpose:    Given a string of errors that correspond to a specific
#             sequence and feature, use the %{$alt_info_HAR} and
#             process that string to determine what (if any) notes,
#             and errors should be added to the feature table
#             for this seq/feature pair, also determine if the stop
#             coordinate should be the predicted stop instead of a
#             possibly corrected one.  
#
# Arguments:
#   $alt_code_str:           string of errors, comma separated, can be ""
#   $seq_name:               name of sequence
#   $ftr_idx:                feature index
#   $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
#   $alt_info_HAR:           REF to hash of arrays with information on the errors, PRE-FILLED
#   $alt_ftr_instances_AHHR: REF to array of 2D hashes with per-feature errors, PRE-FILLED
#   $ret_note_AR:            REF to array of notes, possibly CREATED (not added to) here
#   $ret_alert_AR:           REF to array of errors, possibly added to here (not created)
#   $FH_HR:                  REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies: Never
#################################################################
sub dng_ProcessFeatureAlertsForFTable { 
  my $sub_name = "dng_ProcessFeatureAlertsForFTable";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_code_str, $seq_name, $ftr_idx, $ftr_info_HAR, $alt_info_HAR, $alt_ftr_instances_AHHR, $ret_note_AR, $ret_alert_AR, $FH_HR) = (@_);

  if($alt_code_str eq "") { 
    return 0; 
  }

  # printf("HEYA in $sub_name $seq_name $ftr_idx, $alt_code_str\n");

  # create a hash of all errors in the input $alt_str, and also verify they are all valid errors
  my %input_err_code_H = (); # $input_err_code_H{$alt_code} = 1 if $alt_code is in $alt_code_str
  my @alt_idx_A = ();
  my $alt_code; 
  my $alt_idx; 
  foreach $alt_code (split(",", $alt_code_str)) { 
    $alt_idx = utl_ArrayFindNonNumericValue($alt_info_HAR->{"code"}, $alt_code, $FH_HR);
    if($alt_idx == -1) { 
      ofile_FAIL("ERROR in $sub_name, input error of $alt_code in string $alt_code_str is invalid", "dnaorg", 1, $FH_HR);
    }
    $input_err_code_H{$alt_code} = 1; 
    push(@alt_idx_A, $alt_idx);
  }

  my @tmp_note_A = (); # holds all notes
  my $nerr  = scalar(@alt_idx_A);
  my $valid = 0;
  for(my $e = 0; $e < $nerr; $e++) { 
    $alt_idx = $alt_idx_A[$e];
    $valid = 1; # may be set to '0' below
    if($alt_info_HAR->{"ftbl_invalid_by"}[$alt_idx] ne "") { 
      # printf("\t\tinvalid_by is " . $alt_info_HAR->{"ftbl_invalid_by"}[$alt_idx] . "\n");
      my @invalid_by_err_code_A = split(",", $alt_info_HAR->{"ftbl_invalid_by"}[$alt_idx]);
      foreach my $alt_code2 (@invalid_by_err_code_A) {
        if(exists $input_err_code_H{$alt_code2}) { 
          $valid = 0; # $alt_idx is invalidated by $alt_code2, which is also present in $alt_str
          # printf("\t\t\tinvalidated by $alt_code2\n");
        }
      }
    }
    if($valid) { 
      # valid error that will impact output of feature table
      # add notes and errors
      my $note_str = populateFTableNoteOrAlert("ftbl_note", $alt_idx, $seq_name, $ftr_idx, $ftr_info_HAR, $alt_info_HAR, $alt_ftr_instances_AHHR, undef, $FH_HR);
      if($note_str ne "") { 
        push(@tmp_note_A, $note_str); # we will prune this array and populate @{$ret_note_AR} before returning
      }
      
      my $error_str = populateFTableNoteOrAlert("ftbl_alert", $alt_idx, $seq_name, $ftr_idx, $ftr_info_HAR, $alt_info_HAR, $alt_ftr_instances_AHHR, undef, $FH_HR);
      if($error_str ne "") { 
        # only add the error, if an identical error does not already exist in @{$ret_alert_AR}
        my $idx = utl_ArrayFindNonNumericValue($ret_alert_AR, $error_str, $FH_HR);
        if($idx == -1) { 
          push(@{$ret_alert_AR}, $error_str); 
        }
      }
    }
  }

  # Create the @{$ret_note_AR}, we do not create the @{$ret_alert_AR}
  # this is because @{$ret_note_AR} is per-sequence/feature pair, and caller will treat it as such
  # whereas @{$ret_alert_AR} is per-sequence and caller will treat as such
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

  return;
}

#################################################################
# Subroutine: dng_ProcessSequenceAlertsForFTable()
# Incept:     EPN, Thu Jan 24 12:09:24 2019
#
# Purpose:    Given a string of per-sequence errors that correspond
#             to a specific sequence, use the %{$alt_info_HAR} and
#             process that string to determine what (if any) 
#             errors should be added to the feature table
#             for this sequence. Note that we do not add any 'notes'
#             as we possibly could in processFeatureAlertsForFTable() 
#             because we are dealing with the full sequence and not
#             a feature for a sequence.
#
# Arguments:
#   $alt_code_str:           string of errors, comma separated, can be ""
#   $seq_name:               name of sequence
#   $alt_info_HAR:           REF to hash of arrays with information on the errors, PRE-FILLED
#   $alt_seq_instances_HHR:  REF to 2D hashes with per-sequence errors, PRE-FILLED
#   $ret_alert_AR:           REF to array of errors, possibly added to here (not created)
#   $FH_HR:                  REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies: Never
#################################################################
sub dng_ProcessSequenceAlertsForFTable { 
  my $sub_name = "dng_ProcessSequenceAlertsForFTable";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_code_str, $seq_name, $alt_info_HAR, $alt_seq_instances_HHR, $ret_alert_AR, $FH_HR) = (@_);

  if($alt_code_str eq "") { 
    return 0; 
  }

  #printf("HEYA in $sub_name $seq_name, $alt_code_str\n");

  # NOTE: there's some code duplication in this sub with
  # processFeatureAlertsForFtable(), possibly a chance for additional
  # subroutines

  # create a hash of all errors in the input $alt_str, and also verify they are all valid errors
  my %input_err_code_H = (); # $input_err_code_H{$alt_code} = 1 if $alt_code is in $alt_code_str
  my @alt_idx_A = ();
  my $alt_code; 
  my $alt_idx; 
  foreach $alt_code (split(",", $alt_code_str)) { 
    $alt_idx = utl_ArrayFindNonNumericValue($alt_info_HAR->{"code"}, $alt_code, $FH_HR);
    if($alt_idx == -1) { 
      ofile_FAIL("ERROR in $sub_name, input error of $alt_code in string $alt_code_str is invalid", "dnaorg", 1, $FH_HR);
    }
    $input_err_code_H{$alt_code} = 1; 
    push(@alt_idx_A, $alt_idx);
  }

  my $nerr  = scalar(@alt_idx_A);
  my $valid = 0;
  for(my $e = 0; $e < $nerr; $e++) { 
    $alt_idx = $alt_idx_A[$e];
    $valid = 1; # may be set to '0' below
    if($alt_info_HAR->{"ftbl_invalid_by"}[$alt_idx] ne "") { 
      # printf("\t\tinvalid_by is " . $alt_info_HAR->{"ftbl_invalid_by"}[$alt_idx] . "\n");
      my @invalid_by_err_code_A = split(",", $alt_info_HAR->{"ftbl_invalid_by"}[$alt_idx]);
      foreach my $alt_code2 (@invalid_by_err_code_A) {
        if(exists $input_err_code_H{$alt_code2}) { 
          $valid = 0; # $alt_idx is invalidated by $alt_code2, which is also present in $alt_str
          # printf("\t\t\tinvalidated by $alt_code2\n");
        }
      }
    }
    if($valid) { 
      # add errors
      my $error_str = populateFTableNoteOrAlert("ftbl_alert", $alt_idx, $seq_name, -1, undef, $alt_info_HAR, undef, $alt_seq_instances_HHR, $FH_HR);
      if($error_str ne "") { 
        # only add the error, if an identical error does not already exist in @{$ret_alert_AR}
        my $idx = utl_ArrayFindNonNumericValue($ret_alert_AR, $error_str, $FH_HR);
        if($idx == -1) { 
          push(@{$ret_alert_AR}, $error_str); 
        }
      }
    }
  }

  return;
}

#################################################################
# Subroutine: dng_PopulateFTableNoteOrAlert
# Incept:     EPN, Thu Feb  8 14:31:16 2018
#
# Purpose:    Create notes and errors for the feature table for a specific
#             error, feature, and sequence trio.
#
# Arguments:
#   $ekey:                   either "ftbl_note" or "ftbl_alert"
#   $alt_idx:                index of current error in %{$alt_info_HAR} arrays
#   $seq_name:               name of sequence
#   $ftr_idx:                feature index, -1 if this is a per-sequence error
#   $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
#                            must be undefined if $ftr_idx == -1
#                            must be defined   if $ftr_idx != -1
#   $alt_info_HAR:           REF to hash of arrays with information on the errors, PRE-FILLED
#   $alt_ftr_instances_AHHR: REF to array of 2D hashes with per-feature errors, PRE-FILLED
#                            must be undefined if $ftr_idx == -1
#                            must be defined   if $ftr_idx != -1
#   $alt_seq_instances_HHR:  REF to array of 2D hashes with per-feature errors, PRE-FILLED
#                            must be undefined if $ftr_idx != -1
#                            must be defined   if $ftr_idx == -1
#   $FH_HR:                  REF to hash of file handles, including "log" 
#                            and "cmd"
# 
# Returns: string with the feature table note for the current sequence/feature combo
#
# Dies:    If alt_info_HAR doesn't have information we need
#          or has invalid information
#
#################################################################
sub dng_PopulateFTableNoteOrAlert { 
  my $sub_name = "dng_PopulateFTableNoteOrAlert";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ekey, $alt_idx, $seq_name, $ftr_idx, $ftr_info_HAR, $alt_info_HAR, $alt_ftr_instances_AHHR, $alt_seq_instances_HHR, $FH_HR) = (@_);

  if(! exists $alt_info_HAR->{$ekey}) { 
    ofile_FAIL("ERROR in $sub_name, $ekey value is undefined in error info hash", "dnaorg", 1, $FH_HR);
  }
  # check that combination of $ftr_idx and $alt_ftr_instances_AHHR and $alt_seq_instances_HHR is valid
  if($ftr_idx != -1 && (! defined $alt_ftr_instances_AHHR)) { 
    ofile_FAIL("ERROR in $sub_name, ftr_idx is not -1 but err_ftr_instances_AHHR is not defined", "dnaorg", 1, $FH_HR);
  }
  if($ftr_idx == -1 && (defined $alt_ftr_instances_AHHR)) { 
    ofile_FAIL("ERROR in $sub_name, ftr_idx is -1 but err_ftr_instances_AHHR is defined", "dnaorg", 1, $FH_HR);
  }
  if($ftr_idx != -1 && (! defined $ftr_info_HAR)) { 
    ofile_FAIL("ERROR in $sub_name, ftr_idx is not -1 but ftr_info_HAR is not defined", "dnaorg", 1, $FH_HR);
  }
  if($ftr_idx == -1 && (defined $ftr_info_HAR)) { 
    ofile_FAIL("ERROR in $sub_name, ftr_idx is -1 but ftr_info_HAR is defined", "dnaorg", 1, $FH_HR);
  }
  if($ftr_idx == -1 && (! defined $alt_seq_instances_HHR)) { 
    ofile_FAIL("ERROR in $sub_name, ftr_idx is -1 but err_seq_instances_AHHR is not defined", "dnaorg", 1, $FH_HR);
  }
  if($ftr_idx != -1 && (defined $alt_seq_instances_HHR)) { 
    ofile_FAIL("ERROR in $sub_name, ftr_idx is not -1 but err_ftr_instances_AHHR is defined", "dnaorg", 1, $FH_HR);
  }

  my $msg = $alt_info_HAR->{$ekey}[$alt_idx];

  if(! defined $msg) { 
    ofile_FAIL("ERROR in $sub_name, error $alt_idx is invalid in error info hash", "dnaorg", 1, $FH_HR);
  }

  if($msg eq "") { 
    return "";
  }

  my $alt_code = $alt_info_HAR->{"code"}[$alt_idx];

  my $orig_msg = $msg;
  my $ret_msg  = $msg;
  my $idx;
  # replace !DESC! with description of the error
  if($ret_msg =~ /!DESC!/) { 
    if(($ftr_idx != -1) && (exists $alt_ftr_instances_AHHR->[$ftr_idx]{$alt_code}{$seq_name})) { 
      my $desc_str = sprintf("%s%s", 
                             $alt_info_HAR->{"desc"}[$alt_idx], 
                             ($alt_ftr_instances_AHHR->[$ftr_idx]{$alt_code}{$seq_name} eq "") ? "" : " [" . $alt_ftr_instances_AHHR->[$ftr_idx]{$alt_code}{$seq_name} . "]"); 
      $ret_msg =~ s/!DESC!/$desc_str/g;
    }
    elsif(($ftr_idx == -1) && (exists $alt_seq_instances_HHR->{$alt_code}{$seq_name})) { 
      my $desc_str = sprintf("%s%s", 
                             $alt_info_HAR->{"desc"}[$alt_idx], 
                             ($alt_seq_instances_HHR->{$alt_code}{$seq_name} eq "") ? "" : " [" . $alt_seq_instances_HHR->{$alt_code}{$seq_name} . "]"); 
      $ret_msg =~ s/!DESC!/$desc_str/g;
    }
    else { 
      ofile_FAIL("ERROR in $sub_name, trying to return $ekey message for $alt_code and sequence $seq_name feature $ftr_idx, but no error instance exists", "dnaorg", 1, $FH_HR);
    }
  }
  # replace !FEATURE_TYPE! with 
  if(($ftr_idx != -1) && ($ret_msg =~ /!FEATURE_TYPE!/)) { 
    my $feature_type_str = $ftr_info_HAR->{"type_ftable"}[$ftr_idx];
    $ret_msg =~ s/!FEATURE_TYPE!/$feature_type_str/g;
  }
  # check if there is an internal !$key_str! string, where $key_str is either $key or $key_1,$key_2,...,$key_n for some number n,
  # which is replaced by the value: $ftr_info_HAR->{$key}[$ftr_idx]); for the first $key with a valid value
  if(($ftr_idx != -1) && ($ret_msg =~ /\!([^\!]*)\!/)) {
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
# Subroutine:  dng_ParseListFile()
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
sub dng_ParseListFile {
  my $nargs_expected = 4;
  my $sub_name = "dng_ParseListFile()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name, entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($infile, $do_accn, $line_AR, $FH_HR) = @_;

  # a hash that keeps track of counts of accessions, only used if $do_accn = 1
  my %accn_ct_H    = (); # key: an accession, value: number of times accession occurs (possibly with different versions)
  my @accn_order_A = (); # the accessions read in the input file, in order

  if(-d $infile) { 
    ofile_FAIL("ERROR in $sub_name, trying to read list file $infile, but a directory of the same name exists.", "dnaorg", 1, $FH_HR);
  }

  open(IN, $infile) || fileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  while(my $line = <IN>) { 
    if($line =~ m/\w/) {  # skip blank lines
      chomp $line;
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
      if($do_accn) { 
        my $accn = $line;
        dng_StripVersion(\$accn); # remove version from $accn
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
      ofile_FAIL(sprintf("ERROR in $sub_name, the following accessions occur on multiple lines, possibly with different versions:\n%s", $errmsg), "dnaorg", 1, $FH_HR);
    }
  }
    
  close(IN); 

  return;
}

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
# Subroutine: dng_GetStrandStats()
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
sub dng_GetStrandStats {
  my $sub_name = "dng_GetStrandStats()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($tbl_HHAR, $accn, $FH_HR) = @_;

  my $nfeatures; # number of genes in this genome
  my $npos = 0;  # number of genes on positive strand 
  my $nneg = 0;  # number of genes on negative strand 
  my $nbth = 0;  # number of genes with >= 1 segment on both strands (usually 0)
  my $nunc = 0;  # number of genes with >= 1 segments that are uncertain (usually 0)
  my $strand_str = "";

  if(! exists $tbl_HHAR->{$accn}{"strand"}) { ofile_FAIL("ERROR in $sub_name, didn't read strand information for accn: $accn", "dnaorg", 1, $FH_HR); }

  $nfeatures = scalar(@{$tbl_HHAR->{$accn}{"coords"}});
  if ($nfeatures > 0) { 
    for(my $i = 0; $i < $nfeatures; $i++) { 
      if   ($tbl_HHAR->{$accn}{"strand"}[$i] eq "+") { $npos++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "-") { $nneg++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "!") { $nbth++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "?") { $nunc++; }
      else { ofile_FAIL(sprintf("ERROR in $sub_name, unable to parse strand (%s) for feature %d for $accn\n", $tbl_HHAR->{$accn}{"strand"}[$i], $i+1), "dnaorg", 1, $FH_HR); }
      $strand_str .= $tbl_HHAR->{$accn}{"strand"}[$i];
    }
  }

  return ($nfeatures, $npos, $nneg, $nunc, $nbth, $strand_str);
}
#################################################################
# Subroutine: dng_StartsStopsStrandsFromCoordsLength()
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
#   $strands_AR:  ref to array to fill with strands of each segment, can be undef
#   $nsegments_R: ref to scalar that fill with the number of segments
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void; but fills @{$starts_AR}, @{$stops_AR}, and $$nsegments_R.
#
# Dies:       if we see a feature that spans stop..start but $do_circular is 0
#################################################################
sub dng_StartsStopsStrandsFromCoordsLength { 
  my $sub_name = "dng_StartsStopsStrandsFromCoordsLength()";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($coords, $totlen, $do_circular, $starts_AR, $stops_AR, $strands_AR, $nsegments_R, $FH_HR) = @_;

  # zero/initialize what we will be determining in this subroutine
  @{$starts_AR} = ();
  @{$stops_AR}  = ();
  $$nsegments_R    = 0;
  
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
      ofile_FAIL("ERROR in $sub_name: found internal complement in coords string $coords, we assume all segments are on same strand...", "dnaorg", 1, $FH_HR); 
      $el =~ s/^complement\(//;
      if($cur_strand eq "-") { ofile_FAIL("ERROR in $sub_name, found nested 'complement' annotations in coord string: $coords", "dnaorg", 1, $FH_HR); }
      $cur_strand = "-";
    }
    $el =~ s/\)$//;
    $el =~ s/\<//; # remove '<'
    $el =~ s/\>//; # remove '>'
    if($el =~ m/^(\d+)\.\.(\d+)$/) { 
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $2);
      if(defined $strands_AR) { push(@{$strands_AR}, $cur_strand); }
      $$nsegments_R++;
    }
    elsif($el =~ m/^(\d+)$/) { # a single nucleotide
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $1);
      if(defined $strands_AR) { push(@{$strands_AR}, $cur_strand); }
      $$nsegments_R++;
    }
    else { 
      ofile_FAIL("ERROR unable to parse $orig_coords in $sub_name", "dnaorg", 1, $FH_HR); 
    }
  }

  # check if we have a spanning segment (that spans stop..start) and if we do
  # and (! $do_circular) then die, because that shouldn't happen.
  my $have_spanning_segment = checkForSpanningSequenceSegments($starts_AR, $stops_AR, $nsegments_R, 0, $strand, $totlen); # 1 says: do correct the spanning segment
  if($have_spanning_segment) { 
    if(! $do_circular) { 
      ofile_FAIL("ERROR in $sub_name, found segment that spanned stop..start boundary, but we're not allowing circular genomes...", "dnaorg", 1, $FH_HR); 
    }
    else { 
      # fix it
      checkForSpanningSequenceSegments($starts_AR, $stops_AR, $nsegments_R, 1, $strand, $totlen); # 0 says: don't correct the spanning segment
    }
  }

  return;
}
#################################################################
# Subroutine: dng_StartsStopsFromCoords()
# Incept:     EPN, Thu Feb 11 14:22:54 2016
#
# Purpose:    Extract the starts and stops from a coords string.
#
# Args:
#   $coords:      the coords string
#   $starts_AR:   REF to array to fill with start positions, FILLED HERE
#   $stops_AR:    REF to array to fill with stop positions, FILLED HERE
#   $nsegments_R: REF to scalar that fill with the number of segments, FILLED HERE
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:      void; but fills @{$starts_AR}, @{$stops_AR}, and $$nsegments_R.
#
# Dies:         if we can't parse $coords
#################################################################
sub dng_StartsStopsFromCoords { 
  my $sub_name = "dng_StartsStopsFromCoords()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($coords, $starts_AR, $stops_AR, $nsegments_R, $FH_HR) = @_;

  @{$starts_AR} = ();
  @{$stops_AR}  = ();
  $$nsegments_R    = 0;
  
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
      $$nsegments_R++;
    }
    elsif($el =~ m/^(\d+)$/) { # a single nucleotide
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $1);
      $$nsegments_R++;
    }
    else { 
      ofile_FAIL("ERROR in $sub_name, unable to parse coordinates $orig_coords", "dnaorg", 1, $FH_HR); 
    }
  }

  return;
}

#################################################################
# Subroutine: dng_GetLengthsAndCoords()
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
sub dng_GetLengthsAndCoords { 
  my $sub_name = "dng_GetLengthsAndCoords()";
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
# Subroutine: dng_LengthFromCoords()
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
sub dng_LengthFromCoords { 
  my $sub_name = "dng_LengthFromCoords()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($coords, $FH_HR) = @_;

  my @starts_A = ();
  my @stops_A  = ();
  my $nsegments   = 0;

  startsStopsFromCoords($coords, \@starts_A, \@stops_A, \$nsegments, $FH_HR);

  my $length = 0;
  for(my $i = 0; $i < $nsegments; $i++) { 
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
# Subroutine: dng_DashCoordsStringCommaDelimitedToLength
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
sub dng_DashCoordsStringCommaDelimitedToLength {
  my $sub_name = "dng_DashCoordsStringCommaDelimitedToLength()";
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
# Subroutine: dng_DashCoordsToLength
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
sub dng_DashCoordsToLength { 
  my $sub_name = "dng_DashCoordsToLength";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($start_stop, $caller_sub_name, $FH_HR) = (@_);

  my $len = undef;

  # make sure that we $start_stop is not valid except for the fact that
  # $start or $stop is a negative position (we do not allow negative positions
  # in this function)
  if($start_stop =~ m/^\-+\d+\-\-\d+$/) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, %sstart and stop positions are negative in coords string $start_stop", 
                        (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), "dnaorg", 1, $FH_HR); 
  }
  elsif($start_stop =~ m/^\-\d+\-\d+$/) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, %sstart position is negative in coords string $start_stop", 
                        (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), "dnaorg", 1, $FH_HR); 
  }
  elsif($start_stop =~ m/^\d+\-\-\d+$/) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, %sstop position is negative in coords string $start_stop", 
                        (defined $caller_sub_name) ? "called by $caller_sub_name," : 0), "dnaorg", 1, $FH_HR); 
  }

  # if we get here, $start_stop is either valid, or invalid for a reason other
  # than having a negative position
  if($start_stop =~ m/^(\d+)\-(\d+)$/) { 
    # $start_stop is valid
    $len = (abs($1 - $2) + 1);
  }
  else { 
    # $start_stop is not valid, for some reason other than just having a negative position
    ofile_FAIL("ERROR in $sub_name, called by $caller_sub_name, unable to parse start-stop string: $start_stop", "dnaorg", 1, $FH_HR); 
  }

  return $len;
}
#################################################################
#
# Subroutines for dumping data structures, usually for debugging:
#   dumpInfoHashOfArrays()
#   dumpHashOfHashes()
#   dumpArrayOfHashesOfHashes()
#   dumpArrayOfHashes()
#
#################################################################
# Subroutine: dng_DumpInfoHashOfArrays()
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
sub dng_DumpInfoHashOfArrays { 
  my $sub_name = "dng_DumpInfoHashOfArrays()";
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
      printf $FH ("el %2d\n", ($a_ctr)); 
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
        printf $FH ("\tel %2d: %s\n", ($a_ctr), (defined $info_HAR->{$key}[$a_ctr]) ? $info_HAR->{$key}[$a_ctr] : "undef"); 
      }
      printf $FH ("\n");
    }
  }
  
  return;
}

#################################################################
# Subroutine: dng_DumpHashOfHashes()
# Incept:     EPN, Thu Dec 20 13:36:00 2018
#
# Purpose:    Dump the contents of  hashes of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $HHR:         ref of the hash of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub dng_DumpHashOfHashes { 
  my $sub_name = "dng_DumpHashOfHashes()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");
  
  foreach my $key1 (sort keys %{$HHR}) { 
    printf("*H*H key: $key1\n");
    my $nel = scalar(keys %{$HHR->{$key1}});
    foreach my $key2 (sort keys %{$HHR->{$key1}}) { 
      printf("\tH*H* key: $key2 value: %s\n", $HHR->{$key1}{$key2}); 
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: dng_DumpArrayOfHashesOfHashes()
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
sub dng_DumpArrayOfHashesOfHashes { 
  my $sub_name = "dng_DumpArrayOfHashesOfHashes()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $AHHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $nel1 = scalar(@{$AHHR});
  for(my $i1 = 0; $i1 < $nel1; $i1++) { 
    printf $FH ("*A*HH el %2d\n", ($i1)); 
    my $nel2 = scalar(keys %{$AHHR->[$i1]}); 
    my $i2 = 0;
    foreach my $key2 (sort keys %{$AHHR->[$i1]}) { 
      printf("\tA*H*H el %2d key: $key2\n", ($i2)); 
      $i2++;
      my $nel3 = scalar(keys %{$AHHR->[$i1]{$key2}});
      my $i3 = 0;
      foreach my $key3 (sort keys %{$AHHR->[$i1]{$key2}}) { 
        printf("\tAH*H* el %2d key: $key3 value: %s\n", ($i3), $AHHR->[$i1]{$key2}{$key3}); 
        $i3++;
      }
      printf $FH ("\n");
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: dng_DumpArrayOfHashes()
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
sub dng_DumpArrayOfHashes { 
  my $sub_name = "dng_DumpArrayOfHashes()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $AHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  my $nel1 = scalar(@{$AHR});
  for(my $i1 = 0; $i1 < $nel1; $i1++) { 
    printf $FH ("*A*H el %2d\n", ($i1)); 
    my $nel2 = scalar(keys %{$AHR->[$i1]}); 
    my $i2 = 0;
    foreach my $key2 (sort keys %{$AHR->[$i1]}) { 
      printf("\tA*H* el %2d key: $key2 value: %s\n", ($i2), $AHR->[$i1]{$key2}); 
      $i2++;
    }
    printf $FH ("\n");
  }

  return;
}


#################################################################
# Subroutine: dng_DumpHashOfArraysOfHashes()
# Incept:     EPN, Tue Mar 19 12:30:24 2019
#
# Purpose:    Dump the contents of a hash of arrays of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of array of hashes of hashes
#             $HAHR:        ref of the hash of array of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub dng_DumpHashOfArraysOfHashes { 
  my $sub_name = "dng_DumpHashOfArraysOfHashes()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HAHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  foreach my $key1 (sort keys %{$HAHR}) { 
    printf $FH ("*H*AH key $key1\n");
    my $nel2 = scalar(@{$HAHR->{$key1}});
    for (my $i2 = 0; $i2 < $nel2; $i2++) { 
      printf("\tH*A*H key: $key1 el: $i2:\n", $i2);
      foreach my $key3 (sort keys %{$HAHR->{$key1}[$i2]}) { 
        printf $FH ("\t\tHA*H* key: $key1 el: $i2 key: $key3 value: %s\n", $HAHR->{$key1}[$i2]{$key3});
      }
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: dng_DumpHashOfHashesOfHashes()
# Incept:     EPN, Wed Mar 20 14:44:54 2019
#
# Purpose:    Dump the contents of a hash of hashes of hashes,
#             probably for debugging purposes.
#
# Args:       $name2print:  name of hash of hashes of hashes
#             $HHHR:        ref of the hash of hashes of hashes
#             $FH:          file handle to print (often *STDOUT)
#
# Returns:    void
# 
#################################################################
sub dng_DumpHashOfHashesOfHashes { 
  my $sub_name = "dng_DumpHashOfHashesOfHashes()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($name2print, $HHHR, $FH) = @_;

  printf $FH ("in $sub_name, printing %s:\n", (defined $name2print) ? $name2print : "undefined");

  foreach my $key1 (sort keys %{$HHHR}) { 
    printf $FH ("*H*HH key1: $key1\n");
    foreach my $key2 (sort keys %{$HHHR->{$key1}}) { 
      printf $FH ("\tH*H*H key2: $key1 key2: $key2\n");
      foreach my $key3 (sort keys %{$HHHR->{$key1}{$key2}}) { 
        printf $FH ("\t\tHH*H* key: $key1 key2: $key2 key3: $key3 value: %s\n", $HHHR->{$key1}{$key2}{$key3});
      }
      printf $FH ("\n");
    }
    printf $FH ("\n");
  }

  return;
}

#################################################################
# Subroutine: dng_ArrayOfHashesKeepKeys()
# Incept:     EPN, Wed Mar 13 11:33:26 2019
#
# Purpose:    Delete all key/value pairs from all hashes
#             in an array of hashes, except those in @{$keys_AR}.
#
# Args:       $AHR:      ref of the array of hashes
#             $keys_HR:  ref to hash with key values to keep
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub dng_ArrayOfHashesKeepKeys { 
  my $sub_name = "dng_ArrayOfHashesKeepKeys";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($AHR, $keys_HR) = @_;

  my $nel1 = scalar(@{$AHR});
  for(my $i1 = 0; $i1 < $nel1; $i1++) { 
    foreach my $key (keys %{$AHR->[$i1]}) { 
      if(! exists $keys_HR->{$key}) { 
        delete $AHR->[$i1]{$key};
      }
    }
  }

  return;
}

#################################################################
#################################################################
#
# Subroutines for validating the special data structures:
#   validateExecutableHash()
#   validateFeatureInfoHashIsComplete()
#   validateSegmentInfoHashIsComplete()
#   validateSequenceInfoHashIsComplete()
#   validateAlertInfoHashIsComplete()
#   validateInfoHashOfArraysIsComplete()
#   validateOutputFileInfoHashOfHashes()
#   validateAndGetSizeOfInfoHashOfArrays()
#   getConsistentSizeOfInfoHashOfArrays()
#
#################################################################
# Subroutine:  dng_ValidateExecutableHash()
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
sub dng_ValidateExecutableHash { 
  my $nargs_expected = 2;
  my $sub_name = "dng_ValidateExecutableHash()";
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
    ofile_FAIL("ERROR in $sub_name(),\n$fail_str", "dnaorg", 1, $FH_HR);
  }

  return;
}


#################################################################
# Subroutine: dng_ValidateSequenceInfoHashIsComplete()
# Incept:     EPN, Tue Mar 15 05:35:22 2016
#
# Purpose:    Validate that a 'sequence info' hash is valid and complete.
#             'Complete' means it has all the expected keys, each of which is an identically sized array.
#             The expected keys are:
#                "seq_name":    name of the sequence in the sequence file we create and search in
#                "len":         length of the sequence with name in "seq_name"in the sequence file we create and search in
#                "accn_name":   accession of the sequence in GenBank
#
#             Optional keys are:
#                "ifile_spos":  starting model position for this aligned sequence, read from cmalign --ifile output file, 
#                               -1 if aligned sequence spans 0 model positions
#                "ifile_epos":  ending   model position for this aligned sequence, read from cmalign --ifile output file,
#                               -1 if aligned sequence spans 0 model positions
#                "ifile_ins":   string of all inserts in this aligned sequence, read from cmalign --ifile output file 
#                               "" if no inserts; else format is : 1 or more "<c_x>:<u_x>:<i_x>;" where
#                               <c_x> is a model position; if 0: inserts before 1st consensus posn)
#                               <u_x> is the *unaligned* sequence position of the first inserted residue after <c_x>.
#                               <i_x> is the number of inserted residues after position <c_x>
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
#
#################################################################
sub dng_ValidateSequenceInfoHashIsComplete { 
  my $sub_name = "dng_ValidateSequenceInfoHashIsComplete()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($seq_info_HAR, $exceptions_AR, $opt_HHR, $FH_HR) = (@_);
  
  my @expected_keys_A = ("seq_name", "len", "accn_name", "len");

  my $nseq = validateInfoHashOfArraysIsComplete($seq_info_HAR, \@expected_keys_A, $exceptions_AR, $FH_HR);
  # above call will die if we are invalid

  # make sure we do not have any duplicate accessions/names
  my %name_H = ();
  my %accn_H = ();
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
    my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
    if(exists $name_H{$seq_name}) { 
      ofile_FAIL("ERROR in $sub_name, sequence name $seq_name exists twice", "dnaorg", 1, $FH_HR);
    }
    if(exists $accn_H{$accn_name}) { 
      ofile_FAIL("ERROR in $sub_name, accession $accn_name exists twice", "dnaorg", 1, $FH_HR);
    }
    $name_H{$seq_name} = 1;
    $accn_H{$accn_name} = 1;
  }

  return $nseq;
}

#################################################################
# Subroutine: dng_AlertInfoValidate()
# Incept:     EPN, Fri Mar  4 12:56:43 2016
#
# Purpose:    Validate that a 'error info' hash is valid and complete.
#             'Complete' means it has all the expected keys, each of which is an identically sized array.
#             The expected keys are:
#                "code":            the error code, e.g. "b5e"
#                "pertype":         type of this error, either "feature" (error occurs for a single feature) or "sequence"
#                                   (error occurs for an entire sequence)
#                "desc":            the description/message for this error, e.g. ""unable to identify homologous feature"
#             Information relevant  to feature table output:
#                "ftbl_invalid_by": string that lists all error code indices (separated by commas) the current error code is invalidated by
#                                   for feature table output (if any of those errors are also present, this error becomes invalid and does
#                                   not impact feature table output)
#                "ftbl_note":       note message for feature table
#                "ftbl_alert"       alert message for feature table
#                
#               
#             If @{exceptions_AR} is non-empty, then keys in 
#             in that array need not be in %{$alt_info_HAR}.
#
# Arguments:
#   $alt_info_HAR:  REF to hash of arrays of error information
#   $exceptions_AR: REF to array of keys that may be excluded from the hash
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
# 
# Returns: Number of elements in each and every array of %{$alt_info_HAR}
#
# Dies:    - if one of the expected keys (listed above and not in @{$exceptions_AR})
#            does not exist in $alt_info_HAR
#          - if two arrays in $alt_info_HAR are of different sizes
#          - if any other key other than those listed above exist in ${%alt_info_HAR}
#          - if any key listed in @{$exceptions_AR} is not one of the expected keys
#
#################################################################
sub dng_AlertInfoValidate { 
  my $sub_name = "dng_AlertInfoValidate()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_info_HAR, $exceptions_AR, $FH_HR) = (@_);
  
  my @expected_keys_A = ("code", "pertype", "desc", "causes_failure", "prevents_annot", "ftbl_invalid_by", "ftbl_note", "ftbl_alert");

  # we do not do any more extensive checking because the function that adds elements to the error info hash
  # (dng_AlertInfoAdd()) does a lot of checking to validate the data before it is added

  return dng_ValidateInfoHashOfArraysIsComplete($alt_info_HAR, \@expected_keys_A, $exceptions_AR, $FH_HR);
}

#################################################################
# Subroutine: dng_ValidateInfoHashOfArraysIsComplete()
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
sub dng_ValidateInfoHashOfArraysIsComplete { 
  my $sub_name = "dng_ValidateInfoHashOfArraysIsComplete()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($info_HAR, $expected_keys_AR, $exceptions_AR, $FH_HR) = (@_);
  
  # make sure our exceptions are actually in the expected array
  if(defined $exceptions_AR) { 
    foreach my $key (@{$exceptions_AR}) { 
      if(utl_ArrayFindNonNumericValue($expected_keys_AR, $key, $FH_HR) == -1) { 
        ofile_FAIL("ERROR in $sub_name, excepted value $key is not an expected key in the feature info hash", "dnaorg", 1, $FH_HR);
      }
    }
  }

  # make the list of keys we'll require, this is all expected keys minus any exceptions in @{$exceptions_AR}
  my @reqd_keys_A = ();
  foreach my $key (@{$expected_keys_AR}) { 
    if((! defined $exceptions_AR) || (utl_ArrayFindNonNumericValue($key, $exceptions_AR, $FH_HR) == -1)) { 
      push(@reqd_keys_A, $key);
    }
  }
                         
  my $nftr = dng_ValidateAndGetSizeOfInfoHashOfArrays($info_HAR, \@reqd_keys_A, $FH_HR);  

  return $nftr;
}

#################################################################
# Subroutine: dng_ValidateOutputFileInfoHashOfHashes()
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
sub dng_ValidateOutputFileInfoHashOfHashes { 
  my $sub_name = "dng_ValidateOutputFileInfoHashOfHashes()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($ofile_info_HHR) = (@_);
  
  # we can only pass $FH_HR to ofile_FAIL if that hash already exists
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
      ofile_FAIL("ERROR in $sub_name, unexpected 1d key $key1d exists.", "dnaorg", 1, $FH_HR);
    }     
  } 
  
  # make sure all 2nd dim keys for all 1st dim keys are the same as the 2nd dim keys for 1st dim key "order"
  if(! defined $ofile_info_HHR->{"order"}) { 
    ofile_FAIL("ERROR in $sub_name, expected 1d key order does not exist.", "dnaorg", 1, $FH_HR);
  }
  foreach my $key1d (@same_keys1d_A) { 
    if($key1d ne "order") { # skip "order"
      if(! defined $ofile_info_HHR->{$key1d}) { 
        ofile_FAIL("ERROR in $sub_name, expected 1d key $key1d does not exist.", "dnaorg", 1, $FH_HR);
      }
      # we make sure the set of 2d keys in $ofile_info_HHR->{"order"} and $ofile_info_HHR->{$key1d} are 
      # identical in 2 steps:
      # 1) make sure all 2d keys from $ofile_info_HHR->{"order"} are also in $ofile_info_HHR->{"order"}
      foreach $key2d (keys %{$ofile_info_HHR->{"order"}}) { 
        if(! defined $ofile_info_HHR->{$key1d}{$key2d}) { 
          ofile_FAIL("ERROR in $sub_name, 2nd dim key $key2d exists for ofile_info_HHR->{order} but not for ofile_info_HHR->{$key1d}", "dnaorg", 1, $FH_HR); 
        }
      }
      # 2) make sure all the 2d keys in $ofile_info_HHR->{$key1d} are also in $ofile_info_HHR->{"order"}
      foreach $key2d (keys %{$ofile_info_HHR->{$key1d}}) { 
        if(! defined $ofile_info_HHR->{"order"}{$key2d}) { 
          ofile_FAIL("ERROR in $sub_name, 2nd dim key $key2d exists for ofile_info_HHR->{order} but not for ofile_info_HHR->{$key1d}", "dnaorg", 1, $FH_HR); 
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
      ofile_FAIL("ERROR in $sub_name, invalid values for ofile_info_HH{order}, $nkey2d 2nd dim keys, but value $i does not exist", "dnaorg", 1, $FH_HR);
    }
  }

  return $nkey2d;
}

#################################################################
# Subroutine: dng_ValidateAndGetSizeOfInfoHashOfArrays()
# Incept:     EPN, Thu Feb 11 15:06:40 2016
#
# Purpose:    Validate that the arrays in a hash of arrays are all the
#             same size, and return that size.
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
sub dng_ValidateAndGetSizeOfInfoHashOfArrays {
  my $sub_name = "dng_ValidateAndGetSizeOfInfoHashOfArrays()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HAR, $reqd_keys_AR, $FH_HR) = (@_);
  
  #  make sure that all required keys exist
  if(defined $reqd_keys_AR && @{$reqd_keys_AR}) { 
    foreach my $reqd_key (@{$reqd_keys_AR}) { 
      if(! exists $HAR->{$reqd_key}) { 
        ofile_FAIL("ERROR in $sub_name, required key $reqd_key does not exist in the hash of arrays", "dnaorg", 1, $FH_HR); 
      }
    }
  }

  # get size and check consistency
  my $nel = dng_GetConsistentSizeOfInfoHashOfArrays($HAR, $FH_HR); 
  # this will die if not all elements are the same size, or any values are undefined


  return $nel;
}

#################################################################
# Subroutine: dng_GetConsistentSizeOfInfoHashOfArrays()
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
sub dng_GetConsistentSizeOfInfoHashOfArrays { 
  my $sub_name = "dng_GetConsistentSizeOfInfoHashOfArrays()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($HAR, $FH_HR) = (@_);
  
  if((! %{$HAR}) || (scalar(keys %{$HAR}) == 0)) { 
    ofile_FAIL("ERROR in $sub_name, hash of arrays does not exist or has no keys", "dnaorg", 1, $FH_HR);
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
        ofile_FAIL(sprintf("ERROR in $sub_name, expected number of elements in array for key $key is $nel (from key $nel_key) but %d exist for key $key!\nElement values from key $nel_key are: %s", scalar(@{$HAR->{$key}}), $nel_key_values), "dnaorg", 1, $FH_HR);
      }
    }
    # check that all values are defined
    for(my $i = 0; $i < $nel; $i++) { 
      if(! defined $HAR->{$key}[$i]) { 
        ofile_FAIL("ERROR in $sub_name, undefined value: key: $key index: $i", "dnaorg", 1, $FH_HR); 
      }        
    }
  }

  return $nel;
}

#################################################################
# Subroutine: dng_GetSizeOfInfoHashOfArrays()
# Incept:     EPN, Tue Mar  5 08:00:34 2019
#
# Purpose:    Quickly determine the size of 'info hash' without
#             validating it, using key $key.
#
# Arguments:
#   $info_HAR:  REF to info hash of arrays
#   $key:       key to use to get size
#   $FH_HR:     REF to hash of file handles, including "log" and "cmd"
# 
# Returns: Number of elements in @{$info_HAR->{$key}}
#
# Dies:    - if $info_HAR is undefined or $info_HAR->{$key} does not exist
#
#################################################################
sub dng_GetSizeOfInfoHashOfArrays { 
  my $sub_name = "dng_GetSizeOfInfoHashOfArrays()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($info_HAR, $key, $FH_HR) = (@_);
  
  if(! defined $info_HAR) { 
    ofile_FAIL("ERROR in $sub_name, input info_HAR undefined", "dnaorg", 1, $FH_HR);
  }
  if(! exists $info_HAR->{$key}) { 
    ofile_FAIL("ERROR in $sub_name, key $key does not exist in info_HAR hash", "dnaorg", 1, $FH_HR);
  }

  return(scalar(@{$info_HAR->{$key}}));
}

#################################################################
# Subroutine: dng_ValidateCapitalizedDnaStopCodon()
# Incept:     EPN, Mon Mar 14 13:47:57 2016
# 
# Purpose:    Given an already capitalized DNA codon, return '1' 
#             if it's a valid stop codon, else return 0.
#
# Args:
#  $codon:  the codon
#
# Returns:    The codon as a string
#
#################################################################
sub dng_ValidateCapitalizedDnaStopCodon {
  my $sub_name = "dng_ValidateCapitaliedDnaStopCodon";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($codon) = @_;
  
  if($codon eq "TAA" || 
     $codon eq "TGA" || 
     $codon eq "TAG" || 
     $codon eq "TAR") { 
    return 1;
  }

  return 0;
}


#################################################################
# Subroutine: dng_ValidateCapitalizedDnaStartCodon()
# Incept:     EPN, Sat Feb 23 10:01:55 2019
# 
# Purpose:    Given an already capitalized DNA codon, return '1' 
#             if it's a valid start codon, else return 0.
#
# Args:
#  $codon:  the codon
#
# Returns:    The codon as a string
#
#################################################################
sub dng_ValidateCapitalizedDnaStartCodon {
  my $sub_name = "dng_ValidateCapitaliedDnaStartCodon";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($codon) = @_;
  
  if($codon eq "ATG") { 
    return 1;
  }

  return 0;
}

#################################################################
#################################################################
#
# Simple utility subroutines for hashes and arrays:
#   utl_ArrayFindNonNumericValue()
#   numNonNumericValueInArray()
#   maxLengthScalarKeyInHash()
#   maxLengthScalarValueInHash()
#   maxLengthScalarValueInArray()
#   numberOfDigits()
#   findValueInArray()
#   sumArray()
#   sumHashValues()

#

#################################################################
# Subroutine:  dng_NumNonNumericValueInArray()
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
sub dng_NumNonNumericValueInArray { 
  my $nargs_expected = 3;
  my $sub_name = "dng_NumNonNumericValueInArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR, $value, $FH_HR) = @_;

  if(verify_real($value)) { 
    ofile_FAIL("ERROR in $sub_name, value $value seems to be numeric, we can't compare it for equality", "dnaorg", 1, $FH_HR);
  }

  if(! defined $AR) { 
    ofile_FAIL("ERROR in $sub_name, array reference is not defined", "dnaorg", 1, $FH_HR);
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
# Subroutine:  dng_MaxLengthScalarKeyInHash()
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
sub dng_MaxLengthScalarKeyInHash { 
  my $nargs_expected = 1;
  my $sub_name = "dng_MaxLengthScalarKeyInHash()";
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
# Subroutine:  dng_MaxLengthScalarValueInHash()
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
sub dng_MaxLengthScalarValueInHash { 
  my $nargs_expected = 1;
  my $sub_name = "dng_MaxLengthScalarValueInHash()";
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
# Subroutine:  dng_MaxLengthScalarValueInArray()
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
sub dng_MaxLengthScalarValueInArray { 
  my $nargs_expected = 1;
  my $sub_name = "dng_MaxLengthScalarValueInArray()";
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
# Subroutine:  dng_NumberOfDigits()
# Incept:      EPN, Tue May  9 11:33:50 2017 [ribovore]
#              EPN, Fri Nov 13 06:17:25 2009 [ssu-align:ssu.pm:NumberOfDigits()]
# 
# Purpose:     Return the number of digits in a number before
#              the decimal point. (ex: 1234.56 would return 4).
# Arguments:
# $num:        the number
# 
# Returns:     the number of digits before the decimal point
#
################################################################# 
sub dng_NumberOfDigits { 
    my $nargs_expected = 1;
    my $sub_name = "dng_NumberOfDigits()";
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

    my ($num) = (@_);

    my $ndig = 1; 
    while($num >= 10) { $ndig++; $num /= 10.; }

    return $ndig;
}

#################################################################
# Subroutine: dng_FindValueInArray()
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
sub dng_FindValueInArray { 
  my $sub_name = "dng_FindValueInArray()";
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
# Subroutine:  dng_SumArray()
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
sub dng_SumArray {
  my $nargs_expected = 1;
  my $sub_name = "dng_SumArray()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR) = $_[0];

  my $sum = 0;
  foreach my $el (@{$AR}) { 
    $sum += $el; 
  }
  return $sum;
}

#################################################################
# Subroutine:  dng_SumHashValues()
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
sub dng_SumHashValues {
  my $nargs_expected = 1;
  my $sub_name = "dng_SumHashValues()";
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
#   utl_RunCommand()
#   removeDirPath()
#   removeScriptNameFromString()
#   utl_FileRemoveUsingSystemRm()
#   getMonocharacterString()
#   countLinesInFile()
#   validateFileExistsAndIsNonEmpty()
#   concatenateListOfFiles()
#   md5ChecksumOfFile()
#   nseBreakdown()
#
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
# Subroutine:  dng_RemoveScriptNameFromString()
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
sub dng_RemoveScriptNameFromString { 
  my $sub_name = "dng_RemoveScriptNameFromString()";
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
# Subroutine:  dng_RemoveListOfFiles()
# Incept:      EPN, Fri Oct 19 12:44:05 2018 [ribovore]
#
# Purpose:     Remove each file in an array of file
#              names. If there are more than 100 files, then
#              remove 100 at a time.
# 
# Arguments: 
#   $files2remove_AR:  REF to array with list of files to remove
#   $caller_sub_name:  name of calling subroutine (can be undef)
#   $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:            ref to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies:        If one of the rm -rf commands fails.
#
################################################################# 
sub dng_RemoveListOfFiles { 
  my $nargs_expected = 4;
  my $sub_name = "dng_RemoveListOfFiles()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($files2remove_AR, $caller_sub_name, $opt_HHR, $FH_HR) = @_;

  my $i = 0; 
  my $nfiles = scalar(@{$files2remove_AR});

  while($i < $nfiles) { 
    my $file_list = "";
    my $up = $i+100;
    if($up > $nfiles) { $up = $nfiles; }
    for(my $j = $i; $j < $up; $j++) { 
      $file_list .= " " . $files2remove_AR->[$j];
    }
    my $rm_cmd = "rm $file_list"; 
    utl_RunCommand($rm_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
    $i = $up;
  }
  
  return;
}

#################################################################
# Subroutine: dng_GetMonocharacterString()
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
sub dng_GetMonocharacterString {
  my $sub_name = "dng_GetMonocharacterString";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($len, $char, $FH_HR) = @_;

  if(! verify_integer($len)) { 
    ofile_FAIL("ERROR in $sub_name, passed in length ($len) is not a non-negative integer", "dnaorg", 1, $FH_HR);
  }
  if($len < 0) { 
    ofile_FAIL("ERROR in $sub_name, passed in length ($len) is a negative integer", "dnaorg", 1, $FH_HR);
  }
    
  my $ret_str = "";
  for(my $i = 0; $i < $len; $i++) { 
    $ret_str .= $char;
  }

  return $ret_str;
}

#################################################################
# Subroutine:  dng_CountLinesInFile()
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
sub dng_CountLinesInFile { 
  my $nargs_expected = 2;
  my $sub_name = "dng_CountLinesInFile()";
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
# Subroutine:  dng_FileLinesToArray()
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
sub dng_FileLinesToArray { 
  my $nargs_expected = 4;
  my $sub_name = "dng_FileLinesToArray()";
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
# Subroutine:  dng_ArrayToNewlineDelimitedString()
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
sub dng_ArrayToNewlineDelimitedString {
  my $nargs_expected = 1;
  my $sub_name = "dng_ArrayToNewlineDelimitedString";
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
# Subroutine:  dng_ArrayToString()
# Incept:      EPN, Tue Mar 12 13:07:03 2019
#
# Purpose:     Return a string with all values of an array delimited
#              by $delim_char;
#
# Arguments: 
#   $AR:         ref to array
#   $delim_char: character to delimit with
# 
# Returns:     string
# 
# Dies:        Never.
#
################################################################# 
sub dng_ArrayToString { 
  my $nargs_expected = 2;
  my $sub_name = "dng_ArrayToString";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($AR, $delim_char) = (@_);

  my $retstr = "";
  if(@{$AR}) { 
    $retstr .= $AR->[0];
    for(my $i = 1; $i < scalar(@{$AR}); $i++) { 
      $retstr .= $delim_char . $AR->[$i];
    }
  }
  return $retstr;
}

#################################################################
# Subroutine:  dng_hashKeysToNewlineDelimitedString()
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
sub dng_HashKeysToNewlineDelimitedString {
  my $nargs_expected = 1;
  my $sub_name = "dng_HashKeysToNewlineDelimitedString";
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
# Subroutine:  dng_hashValuesToNewlineDelimitedString()
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
sub dng_HashValuesToNewlineDelimitedString {
  my $nargs_expected = 1;
  my $sub_name = "dng_HashValuesToNewlineDelimitedString";
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
# Subroutine : dng_ValidateFileExistsAndIsNonEmpty()
# Incept:      EPN, Thu May  4 09:30:32 2017
#
# Purpose:     Check if a file exists and is non-empty. 
#
# Arguments: 
#   $filename:         file that we are checking on
#   $filedesc:         description of file
#   $calling_sub_name: name of calling subroutine (can be undef)
#   $do_die:           '1' if we should die if it does not exist.  
#   $FH_HR:            ref to hash of file handles, can be undef
# 
# Returns:     Return '1' if it does and is non empty
#              Return '0' if it does not exist (and ! $do_die)
#              Return '-1' if it exists but is empty (and ! $do_die)
#              Return '-2' if it exists as a directory (and ! $do_die)
#
# Dies:        If file does not exist or is empty and $do_die is 1.
# 
################################################################# 
sub dng_ValidateFileExistsAndIsNonEmpty { 
  my $nargs_expected = 5;
  my $sub_name = "dng_ValidateFileExistsAndIsNonEmpty()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($filename, $filedesc, $calling_sub_name, $do_die, $FH_HR) = @_;

  if(-d $filename) {
    if($do_die) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, %sfile $filename%s exists but is a directory.", 
                         (defined $calling_sub_name ? "called by $calling_sub_name," : ""),
                         (defined $filedesc         ? " ($filedesc)" : "")),
                 "RIBO", 1, $FH_HR); 
    }
    return -2;
  }
  elsif(! -e $filename) { 
    if($do_die) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, %sfile $filename%s does not exist.", 
                         (defined $calling_sub_name ? "called by $calling_sub_name," : ""),
                         (defined $filedesc         ? " ($filedesc)" : "")),
                 "RIBO", 1, $FH_HR); 
    }
    return 0;
  }
  elsif(! -s $filename) { 
    if($do_die) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, %sfile $filename%s exists but is empty.", 
                         (defined $calling_sub_name ? "called by $calling_sub_name," : ""),
                         (defined $filedesc         ? " ($filedesc)" : "")),
                 "RIBO", 1, $FH_HR); 
    }
    return -1;
  }
  
  return 1;
}


#################################################################
# Subroutine:  dng_md5ChecksumOfFile()
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
sub dng_Md5ChecksumOfFile { 
  my $nargs_expected = 4;
  my $sub_name = "dng_Md5ChecksumOfFile()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($file, $caller_sub_name, $opt_HHR, $FH_HR) = @_;

  if(! -s $file) { 
    ofile_FAIL(sprintf("ERROR in $sub_name%s, file to get md5 checksum of ($file) does no exist or is empty", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), "dnaorg", 1, $FH_HR);
  }

  my $out_file = removeDirPath($file . ".md5sum");
  utl_RunCommand("md5sum $file > $out_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);

  open(MD5, $out_file) || fileOpenFailure($out_file, $sub_name, $!, "reading", $FH_HR);
  #194625f7c3e2a5129f9880c7e29f63de  wnv.lin2.matpept.in
  my $md5sum = <MD5>;
  chomp $md5sum;
  if($md5sum =~ /^(\S+)\s+(\S+)$/) { 
    $md5sum = $1;
  }
  else { 
    ofile_FAIL(sprintf("ERROR in $sub_name%s, unable to parse md5sum output: $md5sum", 
                        (defined $caller_sub_name) ? " called by $caller_sub_name" : ""), "dnaorg", 1, $FH_HR);
  }
  close(MD5);

  utl_FileRemoveUsingSystemRm($out_file, $caller_sub_name, $opt_HHR, $FH_HR);

  return $md5sum;
}

#################################################################
# Subroutine:  dng_nseBreakdown()
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
sub dng_NseBreakdown {
  my $nargs_expected = 1;
  my $sub_name = "dng_NseBreakdown()";
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
#   dng_StripVersion()
#   fetchedNameToListName()
#   fetchSequencesUsingEslFetchCds()
#   addNameAndBlankSsToStockholmAlignment()
#   getQualifierValues()
#   createCmDb()
#   matpeptValidateCdsRelationships()
#   checkForSpanningSequenceSegments()
#   getIndexHashForArray()
#
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

#################################################################
# Subroutine: dng_FetchedNameToListName()
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
sub dng_FetchedNameToListName { 
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
# Subroutine:  dng_getIndexHashForArray()
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
sub dng_GetIndexHashForArray { 
  my $nargs_expected = 3;
  my $sub_name = "dng_GetIndexHashForArray()";
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
      ofile_FAIL("ERROR in $sub_name(), value $el is numeric"); 
    }
    if(exists $index_HR->{$el}) { 
      ofile_FAIL("ERROR in $sub_name(), the value $el appears twice in the array"); 
    }
    $index_HR->{$el} = $i;
  }

  return;
}

#################################################################
# Subroutine:  dng_waitForFarmJobsToFinish()
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
  my $nargs_expected = 9;
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
            my $success = dng_cmalignCheckStdOutput($outfile_A[$i], 
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
  my $nfiles_created = dng_CountLinesInFile($outfile, $FH_HR);

  if(! opt_Get("--keep", $opt_HHR)) { 
    utl_RunCommand("rm $outfile", opt_Get("-v", $opt_HHR), 0, $FH_HR);
  }

  return $nfiles_created;
}


#################################################################
# Subroutine:  dng_CmsearchWrapper()
# Incept:      EPN, Mon Mar 18 14:44:46 2019
#
# Purpose:     Run one or more cmsearch jobs on the farm
#              or locally, after possibly splitting up the input
#              sequence file with dng_SplitFastaFile and 
#              then calling dng_CmalignOrCmsearchWrapperHelper(). 
#
# Arguments: 
#  $execs_HR:              ref to executables with "esl-ssplit" and "cmsearch"
#                          defined as keys
#  $round:                 round index (1 or 2)
#  $out_root:              string for naming output files
#  $seq_file:              name of sequence file with all sequences to run against
#  $tot_len_nt:            total length of all nucleotides in $seq_file
#  $progress_w:            width for outputProgressPrior output
#  $mdl_file:              name of model file to use
#  $mdl_name:              name of model to fetch from $mdl_file (undef to not fetch)
#  $opt_HHR:               REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:        REF to 2D hash of output file information
#
# Returns:     void, updates $$nfa_created_R with number of
#              fasta files created.
# 
# Dies: If an executable doesn't exist, or cmalign or nhmmscan or esl-ssplit
#       command fails if we're running locally
################################################################# 
sub dng_CmsearchWrapper { 
  my $sub_name = "dng_CmsearchWrapper";
  my $nargs_expected = 10;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $round, $out_root, $seq_file, $tot_len_nt, $progress_w,
      $mdl_file, $mdl_name, $opt_HHR, $ofile_info_HHR) = @_;

  my $nfasta_created = 0; # number of fasta files created by esl-ssplit
  my $log_FH = $ofile_info_HHR->{"FH"}{"log"}; # for convenience
  my $start_secs; # timing start
  my $do_parallel = opt_Get("-p", $opt_HHR);

  # set up output file names
  my @concat_keys_A = (); # %out_file_HAR keys we are going to concatenate files for
  my %concat_HA = ();     # hash of arrays of all files to concatenate together
  my $out_key;            # key for an output file: e.g. "stdout", "ifile", "tfile", "tblout", "err"
  @concat_keys_A = ("tblout");
  if($do_parallel) { push(@concat_keys_A, "err"); }
  foreach my $out_key (@concat_keys_A) { 
    @{$concat_HA{$out_key}} = ();
  }    

  my @seq_file_A  = (); # [0..$nr-1]: name of sequence file for this run
  my %out_file_HA = (); # hash of arrays ([0..$nr-1]) of output files for cmsearch runs (filled by dng_CmalignOrCmsearchWrapperHelper())
  my $nr = 1;

  if($do_parallel) { 
    # -p used: we need to split up the sequence file, and submit a separate 
    # cmsearch job for each
    my $targ_nseqfiles = dng_SplitNumSeqFiles($tot_len_nt, $opt_HHR);
    if($targ_nseqfiles > 1) { # we are going to split up the fasta file 
      $nr = dng_SplitFastaFile($execs_HR->{"esl-ssplit"}, $seq_file, $targ_nseqfiles, $opt_HHR, $ofile_info_HHR);
      # dng_SplitFastaFile will return the actual number of fasta files created, 
      # which can differ from the requested amount (which is $targ_nseqfiles) that we pass in. 
      for(my $i = 0; $i < $nr; $i++) { 
        $seq_file_A[$i] = $seq_file . "." . ($i+1);
      }
    }
    else { # targ_nseqfiles is 1, no need to split
      $seq_file_A[0] = $seq_file;
    }
  }
  else { # -p not used
    $seq_file_A[0] = $seq_file;
  }
    
  dng_CmalignOrCmsearchWrapperHelper($execs_HR, $out_root . ".r" . $round, $progress_w, $mdl_file, $mdl_name, 
                                     \@seq_file_A, \%out_file_HA, undef, undef, $opt_HHR, $ofile_info_HHR); 
  
  # concatenate all the tblout and possibly .err files into one 
  foreach $out_key (@concat_keys_A) { 
    my $concat_file = $out_root . ".r" . $round . ".cmsearch." . $out_key;
    dng_ConcatenateListOfFiles($out_file_HA{$out_key}, $concat_file, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
    # dng_ConcatenateListOfFiles() removes individual files unless --keep enabled
    my $out_root_key = sprintf("r%s.concat.%s$out_key", $round, (defined $mdl_name) ? $mdl_name . "." : "");
    ofile_AddClosedFileToOutputInfo($ofile_info_HHR, "dnaorg", $out_root_key, $concat_file, 0, "round $round search $out_key file");
  }

  # remove sequence files if we created any
  if(($do_parallel) && (! opt_Get("--keep", $opt_HHR))) { 
    dng_RemoveListOfFiles(\@seq_file_A, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
  }

  return;
}


#################################################################
# Subroutine:  dng_CmalignOrCmsearchWrapperHelper()
# Incept:      EPN, Mon Mar 18 14:37:54 2019
#
# Purpose:     Run one or more cmalign or cmsearch jobs on the farm
#              or locally.
#
#              Helper subroutine for cmalignWrapper() and
#              cmsearchWrapper() subroutines.
#              see those sub's "Purpose" for more details.
#
# Arguments: 
#  $execs_HR:              ref to hash with paths to cmalign, cmsearch and cmfetch
#  $out_root:              string for naming output files
#  $progress_w:            width for outputProgressPrior output
#                          also serves as flag for what 'round' (1 or 2) we are in
#                          $progress_w will be > 0 for round 1 and < 0 for round 2
#                          in round 2 $progress_w is actually -1 * $progress_w
#  $mdl_file:              name of model file to use
#  $mdl_name:              name of model to fetch from $mdl_file (undef to not fetch)
#  $seq_file_AR:           ref to array of sequence file names for each cmalign/nhmmscan call, PRE-FILLED
#  $out_file_HAR:          ref to hash of arrays of output file names, FILLED HERE 
#  $success_AR:            ref to array of success values, can be undef if $executable is "cmsearch"
#                          $success_AR->[$j] set to '1' if job finishes successfully
#                                            set to '0' if job fails due to mx overflow (must be cmalign)
#  $mxsize_AR:             ref to array of required matrix sizes, can be undef if $executable is "cmsearch"
#                          $mxsize_AR->[$j] set to value readh from cmalign output, if $success_AR->[$j] == 0
#                                           else set to '0'
#  $opt_HHR:               REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:        REF to 2D hash of output file information
#
# Returns:     void
# 
# Dies: If an executable doesn't exist, or command fails (and its not a cmalign allowed failure)
#
################################################################# 
sub dng_CmalignOrCmsearchWrapperHelper { 
  my $sub_name = "dng_CmalignOrCmsearchWrapperHelper";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $out_root, $progress_w, $mdl_file, $mdl_name, $seq_file_AR, $out_file_HAR, $success_AR, $mxsize_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $do_cmalign     = (defined $mxsize_AR) ? 1 : 0;
  my $do_parallel    = opt_Get("-p", $opt_HHR) ? 1 : 0;
  my $program_choice = ($do_cmalign) ? "cmalign" : "cmsearch";
  my $nfasta_created = 0;     # number of fasta files created by esl-ssplit
  my $log_FH         = $ofile_info_HHR->{"FH"}{"log"}; # for convenience
  my $nseq_files     = scalar(@{$seq_file_AR});

  # determine description of the runs we are about to do, 
  # depends on $do_parallel, $round, and ($progress_w < 0), and 
  my $desc = "";
  if($do_parallel) { 
    if($progress_w < 0) { # flag for 'rerunning'
      $desc = "Resubmitting $nseq_files $program_choice jobs to the farm to find seqs too divergent to align";
      $progress_w *= -1;
    }
    else { 
      $desc = "Submitting $nseq_files $program_choice jobs to the farm"
    }
  }      
  else { 
    if($progress_w < 0) { # flag for 'rerunning'
      $desc = "Rerunning $program_choice locally to find sequences that are too divergent to align";
      $progress_w *= -1;
    }
    else { 
      $desc = "Running $program_choice locally";
    }
  }
  my $start_secs = ofile_OutputProgressPrior($desc, $progress_w, $log_FH, *STDOUT);

  # define output file names
  %{$out_file_HAR} = ();
  my @out_keys_A = ();
  my $s;
  my $key;
  if($do_cmalign) { 
    @out_keys_A = ("stdout", "err", "ifile", "tfile", "stk");
  }
  else { 
    @out_keys_A = ("stdout", "err", "tblout");
  }
  foreach $key (@out_keys_A) { 
    @{$out_file_HAR->{$key}} = ();
  }

  for($s = 0; $s < $nseq_files; $s++) { 
    foreach $key (@out_keys_A) { 
      $out_file_HAR->{$key}[$s] = $out_root . ".s" . $s . "." . $key;
    }
    if($do_cmalign) { 
      $success_AR->[$s] = dng_RunCmalign($execs_HR, $mdl_file, $mdl_name, $seq_file_AR->[$s], 
                                         $out_file_HAR->{"stdout"}[$s], $out_file_HAR->{"ifile"}[$s], $out_file_HAR->{"tfile"}[$s], 
                                         $out_file_HAR->{"stk"}[$s], $out_file_HAR->{"err"}[$s],
                                         (defined $mxsize_AR) ? \$mxsize_AR->[$s] : undef, 
                                         $opt_HHR, $ofile_info_HHR);   
    }
    else { 
      dng_RunCmsearch($execs_HR, $mdl_file, $mdl_name, $seq_file_AR->[$s], 
                      $out_file_HAR->{"stdout"}[$s], $out_file_HAR->{"tblout"}[$s], $out_file_HAR->{"err"}[$s],
                      $opt_HHR, $ofile_info_HHR);   
      # we don't bother with @{$success_AR->[$s]} for cmsearch, if dng_runCmsearch() returns the run was successful
    }
    # if we are running parallel, ignore the return values from the run{Cmalign,Cmsearch} subroutines
    # dng_WaitForFarmJobsToFinish() will fill these later
    if(($do_parallel) && (defined $success_AR)) { $success_AR->[$s] = 0; }
  }
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  if($do_parallel) { 
    if((opt_Exists("--skipalign", $opt_HHR)) && (opt_Get("--skipalign", $opt_HHR))) { 
      for($s = 0; $s < $nseq_files; $s++) { 
        $success_AR->[$s] = 1; 
      }
    }
    else { 
      # --skipalign not enabled
      # wait for the jobs to finish
      $start_secs = ofile_OutputProgressPrior(sprintf("Waiting a maximum of %d minutes for all farm jobs to finish", opt_Get("--wait", $opt_HHR)), 
                                              $progress_w, $log_FH, *STDOUT);
      my $njobs_finished = dng_WaitForFarmJobsToFinish($do_cmalign, 
                                                       ($do_cmalign) ? $out_file_HAR->{"stdout"} : $out_file_HAR->{"tblout"}, 
                                                       $out_file_HAR->{"err"}, 
                                                       $success_AR, # this may be undef
                                                       $mxsize_AR,  # this may be undef
                                                       ($do_cmalign) ? "" : "[ok]", # value is irrelevant for cmalign
                                                       opt_Get("--wait", $opt_HHR), opt_Get("--errcheck", $opt_HHR), $ofile_info_HHR->{"FH"});
      if($njobs_finished != $nseq_files) { 
        ofile_FAIL(sprintf("ERROR in $sub_name only $njobs_finished of the $nseq_files are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", $opt_HHR)), 1, $ofile_info_HHR->{"FH"});
      }
      ofile_OutputString($log_FH, 1, "# "); # necessary because waitForFarmJobsToFinish() creates lines that summarize wait time and so we need a '#' before 'done' printed by ofile_OutputProgressComplete()
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
# Subroutine: dng_FormatTabDelimitedStringForAlertListFile()
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
# Returns:    $alt_tab_str: tab delimited string in format described above.
#
# Dies: If $errstr is not in required format
#
#################################################################
sub dng_FormatTabDelimitedStringForAlertListFile() {
  my $sub_name  = "formatTabDelimitedStringForAlertListFile";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($seqname, $errstr, $FH_HR) = (@_);

  # first replace any '_' characters with single spaces
  chomp $errstr;
  $errstr =~ s/\_/ /g;

  my $error_name   = undef;
  my $feature_name = undef;
  my $error_desc   = undef;

  if($errstr =~ /^([^\:]+)\:\s*(.+)$/) {
    ($error_name, $error_desc) = ($1, $2);
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
    $feature_name = "*sequence*";
  }
  elsif($errstr =~ /^([^\[\:]+)\;$/) {
    # example from dnaorg_classify.pl:
    #No Annotation;
    ($error_name) = ($1);
    $feature_name = "*sequence*";
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse input error string: $errstr", "dnaorg", 1, $FH_HR);
  }

  if($error_desc eq "") { 
    $error_desc = "-";
  }
  return $seqname . "\t" . $error_name . "\t" . $feature_name . "\t" . $error_desc;
}

#################################################################
# Subroutine: dng_BlastxDbSeqnameToFtrIdx()
# Incept:      EPN, Tue Dec 18 13:27:50 2018
#
# Purpose:    Find the feature $ftr_idx that corresponds to the blastx
#             db sequence that was named with the convention:
#
#             <protein-accession>/<coords-str>
#
#             Where <coords-str> is identical to $ftr_info_HAR->{"ref_coords"}[$ftr_idx].
#
# Arguments: 
#  $blastx_seqname: sequence name
#  $ftr_info_HAR:   ref to the feature info array of hashes 
#  $FH_HR:          ref to hash of file handles
#
# Returns:    <$ftr_idx>
#
# Dies:       If we find zero features that match to this sequence
#             If we find more than 1 features that match to this sequence
#
################################################################# 
sub dng_BlastxDbSeqNameToFtrIdx { 
  my $sub_name = "dng_BlastxDbSeqNameToFtrIdx";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastx_seqname, $ftr_info_HAR, $FH_HR) = @_;

  my $nftr = getInfoHashSize($ftr_info_HAR, "type", $FH_HR); # nftr: number of features

  my $ret_ftr_idx = undef;
  if($blastx_seqname =~ /(\S+)\/(\S+)/) { 
    my ($accn, $coords) = ($1, $2);
    # find it in @{$ftr_info_HAR->{"coords"}}
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if(($ftr_info_HAR->{"type"}[$ftr_idx] eq "CDS")) { 
        if($ftr_info_HAR->{"coords"}[$ftr_idx] eq $coords) { 
          if(defined $ret_ftr_idx) { # found more than 1 features that match
            ofile_FAIL("ERROR in $sub_name, found blastx db sequence with coords that match two features, ftr_idx: $ftr_idx and $ret_ftr_idx", "dnaorg", 1, $FH_HR);
          }                  
          $ret_ftr_idx = $ftr_idx;
        }
      }
    }
    if(! defined $ret_ftr_idx) { # did not find match
      ofile_FAIL("ERROR in $sub_name, did not find matching feature for blastx db sequence $blastx_seqname", "dnaorg", 1, $FH_HR);
    }
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse blastx db sequence name $blastx_seqname", "dnaorg", 1, $FH_HR); 
  }

  return $ret_ftr_idx;
}

#################################################################
# Subroutine: dng_ValidateBlastDbExists()
# Incept:      EPN, Tue Dec 18 15:32:50 2018
#
# Purpose:    Validate that a blast database exists.
#
# Arguments: 
#  $blastdb_name:  name of the blast db
#  $FH_HR:         ref to hash of file handles
#
# Returns:    void
#
# Dies:       If any of the required files for a blast db do not exist.
#
################################################################# 
sub dng_ValidateBlastDbExists {
  my $sub_name = "dng_ValidateBlastDbExists";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastdb_name, $FH_HR) = @_;

  foreach my $sfx (".phr", ".pin", ".psq") { 
    if(! -s $blastdb_name . $sfx) { 
      ofile_FAIL("ERROR in $sub_name, required blast DB file " . $blastdb_name . $sfx . " does not exist or is empty", "dnaorg", 1, $FH_HR); 
    }
  }

  return;
}


#################################################################
# Subroutine: dng_FeatureTypeIsCdsOrMaturePeptide()
# Incept:      EPN, Mon Feb 25 14:30:34 2019
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
sub dng_FeatureTypeIsCdsOrMaturePeptide { 
  my $sub_name = "dng_FeatureTypeIsCdsOrMaturePeptide";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return((featureTypeIsCds($ftr_info_AHR, $ftr_idx)) || 
         (featureTypeIsMaturePeptide($ftr_info_AHR, $ftr_idx)));
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

  my ($ftr_info_HAR, $ftr_idx) = @_;

  return(($ftr_info_HAR->{"source_idx"}[$ftr_idx] != $ftr_idx) ? 1 : 0);
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
# Subroutine: dng_FeatureInfoSummarizeSegment()
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
sub dng_FeatureInfoSummarizeSegment { 
  my $sub_name = "dng_FeatureInfoSummarizeSegment";
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
# Subroutine: dng_FeatureInfoImputeCoords
# Incept:     EPN, Wed Mar 13 13:15:33 2019
# 
# Purpose:    Fill "coords" values in %{$ftr_info_HAR}
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
    $ftr_info_AHR->[$ftr_idx]{"coords"} = dng_FeatureCoordsFromLocation($ftr_info_AHR->[$ftr_idx]{"location"}, $FH_HR);
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
# Subroutine: dng_FeatureInfoStartStopStrandArrays()
# Incept:     EPN, Fri Mar 15 15:39:35 2019
#
# Synopsis: For all features in a @{$ftr_info_AHR}, validate 
#           "coords" values and fill @{$start_AAR}, @{$stop_AAR} and
#           @{$strand_AAR} based on them.
# 
# Arguments:
#  $ftr_AHR:       REF to the feature info array of hashes
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
  my $nseg = scalar(@coords_A);
  for($sgm_idx = 0; $sgm_idx < $nseg; $sgm_idx++) { 
    if($coords_A[$sgm_idx] =~ /^\<?(\d+)\.\.\>?(\d+)\:([\+\-])$/) { 
      ($start, $stop, $strand) = ($1, $2, $3);
      push(@start_A,  $start);
      push(@stop_A,   $stop);
      push(@strand_A, $strand); 
    }
    else { 
      ofile_FAIL("ERROR in $sub_name, unable to parse coords token $coords_A[$sgm_idx]", "dnaorg", 1, $FH_HR); 
    }
  }

  if(defined $start_AR)  { @{$start_AR}   = @start_A;  }
  if(defined $stop_AR)   { @{$stop_AR}    = @stop_A;   }
  if(defined $strand_AR) { @{$strand_AR}  = @strand_A;  }

  return;
}




#################################################################
# Subroutine: dng_FeatureCoordsFromLocation
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
sub dng_FeatureCoordsFromLocation { 
  my $sub_name = "dng_FeatureCoordsFromLocation";
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
    $ret_val = dng_FeatureCoordsFromLocation($location_to_join, $FH_HR);
  }
  elsif($location =~ /^complement\((.+)\)$/) { 
    my $location_to_complement = $1;
    my $coords_to_complement = dng_FeatureCoordsFromLocation($location_to_complement, $FH_HR);
    $ret_val = dng_FeatureCoordsComplement($coords_to_complement, $FH_HR);
  }
  elsif($location =~ /\,/) { 
    # not wrapped in join() or complement(), but multiple segments
    foreach my $location_el (split(",", $location)) { 
      if($ret_val ne "") { $ret_val .= ","; }
      $ret_val .= dng_FeatureCoordsFromLocation($location_el, $FH_HR);
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
# Subroutine: dng_FeatureCoordsComplement
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
sub dng_FeatureCoordsComplement { 
  my $sub_name = "dng_FeatureCoordsComplement";
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

  printf("\tin $sub_name, coords: $coords ret_val: $ret_val\n");

  return $ret_val;
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
# Subroutine:  dng_sqstringAddNewlines()
# Incept:      EPN, Thu Mar 14 06:12:11 2019
#
# Purpose:     Add newlines to $sqstring after every $linelen
#              characters and return result.
#
# Arguments: 
#   $sqstring: the sequence string
#   $linelen:  interval for newlines
# 
# Returns:     $sqstring with newlines inserted every $linelen 
#              characters and at end of string.
# 
# Dies:        Never.
#
################################################################# 
sub dng_SqstringAddNewlines { 
  my $nargs_expected = 2;
  my $sub_name = "dng_SqstringAddNewlines";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($sqstring, $linelen) = @_;

  my $retstr = "";
  my $sqpos = 0;
  my $sqlen = length($sqstring);

  while($sqpos < $sqlen) { 
    $retstr .= substr($sqstring, $sqpos, $linelen) . "\n";
    $sqpos += $linelen;
  }

  return $retstr;
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
# Subroutine: dng_GenbankParse()
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
sub dng_GenbankParse { 
  my $sub_name = "dng_GenbankParse";
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
            dng_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
          }
          ($qualifier, $value) = ($save_qualifier, $save_value);
        }
        elsif($line =~ /^\s+(\S+)\s+(\S+)$/) { 
          # NOTE: this will pass for a non-first line of a qualifier value that has whitespace in it:
          # e.g.                      KQP ASRDESQKPPRPPTPELVKRIPPPPPNGEEEEEPVIRYEVKSGISGLPELTTVPQ
          # But I think those are illegal, if they're not, then we'll set "KQP" as feature below, which is bad
          if(defined $value) { # we are finished with previous value
            dng_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qualifier, $value, $FH_HR);
            ($qualifier, $value) = (undef, undef);
          }
          # feature/location line, examples:
          #   gene            5..5104
          ($feature, $location) = ($1, $2);
          $ftr_idx++;
          dng_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",     $feature,  $FH_HR);
          dng_GenbankStoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "location", $location, $FH_HR);
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
        dng_GenbankStoreQualifierValue($ftr_info_HAHR->{$acc}, $ftr_idx, $qualifier, $value, $FH_HR);
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
# Subroutine: dng_GenbankStoreQualifierValue()
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
sub dng_GenbankStoreQualifierValue { 
  my $sub_name = "dng_GenbankStoreQualifierValue";
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
# Subroutine:  dng_verifyEnvVariableIsValidDir()
# Incept:      EPN, Wed Oct 25 10:09:28 2017 [ribo.pm]
#
# Purpose:     Verify that the environment variable $envvar exists 
#              and that it is a valid directory. Return directory path.
#              
# Arguments: 
#   $envvar:  environment variable
#
# Returns:    directory path $ENV{'$envvar'}
#
################################################################# 
sub dng_VerifyEnvVariableIsValidDir { 
  my $nargs_expected = 1;
  my $sub_name = "dng_VerifyEnvVariableIsValidDir()";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($envvar) = $_[0];

  if(! exists($ENV{"$envvar"})) { 
    die "ERROR, the environment variable $envvar is not set";
    # it's okay this isn't ofile_FAIL because this is called before ofile_info_HH is set-up
  }
  my $envdir = $ENV{"$envvar"};
  if(! (-d $envdir)) { 
    die "ERROR, the directory specified by your environment variable $envvar does not exist.\n"; 
    # it's okay this isn't ofile_FAIL because this is called before ofile_info_HH is set-up
  }    

  return $envdir;
}


#################################################################
# Subroutine: dng_SqstringCapitalize
# Incept:     EPN, Fri Mar 15 13:32:36 2019
# 
# Purpose:    Capitalize a string in place.
# 
# Arguments:
#   $sqstring_R: REF to sequence string to capitalize
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub dng_SqstringCapitalize {
  my $sub_name = "dng_SqstringCapitalize";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($sqstring_R) = @_;
  
  $$sqstring_R =~ tr/a-z/A-Z/;
  return;
}

#################################################################
# Subroutine: dng_SqstringDnaize
# Incept:     EPN, Fri Mar 15 13:33:39 2019
# 
# Purpose:    Convert a RNA/DNA sqstring to DNA in place.
# 
# Arguments:
#   $sqstring_R: REF to sequence string to capitalize
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub dng_SqstringDnaize {
  my $sub_name = "dng_SqstringDnaize";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($sqstring_R) = @_;
  
  $$sqstring_R =~ tr/Uu/Tt/;
  return;
}


#################################################################
# Subroutine: dng_SqstringReverseComplement
# Incept:     EPN, Fri Mar 15 15:35:10 2019
# 
# Purpose:    Reverse complement a RNA/DNA sqstring in place.
# 
# Arguments:
#   $sqstring_R: REF to reverse complement
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub dng_SqstringReverseComplement {
  my $sub_name = "dng_SqstringReverseComplement";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($sqstring_R) = @_;

  # DNA-ize it
  sqstringDnaize($sqstring_R);
  # reverse it 
  $$sqstring_R = reverse $$sqstring_R;
  # complement it
  $$sqstring_R =~ tr/ACGTRYMKHBVDacgtrymkhbvd/TGCAYRKMDVBHtgcayrkmdvbh/;
  # see esl_alphabet.c::set_complementarity()
  # note that S, W, N are omitted as they are their own complements

  return;
}


#################################################################
# Subroutine: dng_SqstringDiffSummary
# Incept:     EPN, Fri Mar 15 13:35:28 2019
# 
# Purpose:    Return a string summarizes the differences between
#             two sqstrings.
# 
# Arguments:
#   $sqstring1: sqstring 1
#   $sqstring2: sqstring 2
#
# Returns:    String with N newlines for N differences.
#             "" if sqstrings are identical.
# 
# Dies:       never
#
#################################################################
sub dng_SqstringDiffSummary {
  my $sub_name = "dng_SqstringDiffSummary";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($sqstring1, $sqstring2) = @_;

  if(! defined $sqstring1) { 
    return "sequence 1 is undefined\n"; 
  }
  if(! defined $sqstring2) { 
    return "sequence 2 is undefined\n"; 
  }
  dng_SqstringCapitalize(\$sqstring1);
  dng_SqstringCapitalize(\$sqstring2);
  dng_SqstringDnaize(\$sqstring1);
  dng_SqstringDnaize(\$sqstring2);
  if($sqstring1 eq $sqstring2) { 
    return "";
  }
  my $n1 = length($sqstring1); 
  my $n2 = length($sqstring2); 
  if($n1 != $n2) {
    return "sequence lengths mismatch ($n1 != $n2)\n"; 
  }
  my $ret_str = "";
  my @sqstring1_A = split("", $sqstring1); 
  my @sqstring2_A = split("", $sqstring2); 
  my $n = ($n1 > $n2) ? $n1 : $n2;
  for(my $i = 0; $i < $n; $i++) { 
    if($i >= $n1) { 
      $ret_str .= " char " . ($i+1) . " seq1: off-end seq2: " . $sqstring2_A[$i] . "\n";
    }
    elsif($i >= $n2) { 
      $ret_str .= " char " . ($i+1) . " seq1: " . $sqstring1_A[$i] . " seq2: off-end\n";
    }
    elsif($sqstring1_A[$i] ne $sqstring2_A[$i]) { 
      $ret_str .= " char " . ($i+1) . " seq1: " . $sqstring1_A[$i] . " seq2: " . $sqstring2_A[$i] . "\n";
    }
  }

  return $ret_str;
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
        print $out_FH(">" . $msa->get_sqname($seq_idx) . "/" . $ftr_info_AHR->[$ftr_idx]{"coords"} . "\n" . dng_SqstringAddNewlines($cds_sqstring, 60));
      }
    }
  }
  return;
}

#################################################################
# Subroutine: dng_CdsTranslateToFastaFile()
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
sub dng_CdsTranslateToFastaFile { 
  my $sub_name = "dng_CdsTranslateToFastaFiles";
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
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") { 
        my $fetch_name = "source=" . $seq_name . ",coords=1.." . ($seq_length - 3); # subtract length of stop codon
        if(! $protein_sqfile->check_seq_exists($fetch_name)) { 
          ofile_FAIL(sprintf("ERROR in $sub_name, problem translating CDS feature, unable to find expected translated sequence in $tmp2_translate_fa_file:\n\tseq: $seq_name\n\tftr_idx: $ftr_idx\n\tcoords: %s\n\texpected sequence:$fetch_name\n", $ftr_info_AHR->[$ftr_idx]{"coords"}), "dnaorg", 1, $FH_HR);
        }
        print $out_FH $protein_sqfile->fetch_seq_to_fasta_string($fetch_name, 60);
      }
    }
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
# Dies:       if $ftr_info_HAR is not valid upon entering
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
            # printf("\tadded ftr_info_HAR->{$mdl_name}[$ftr_idx]{$key} as $value\n");
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
# Subroutine: dng_HashFromCommaSeparatedString
# Incept:     EPN, Mon Mar 18 06:52:21 2019
#
# Synopsis: Given a hash reference and a comma separated string
#           fill the hash with keys for each token in the string,
#           with all values set as 1.
#
# Arguments:
#  $HR:      hash reference
#  $string:  comma separated string
#
# Returns:    void
#
# Dies:       never
#################################################################
sub dng_HashFromCommaSeparatedString {
  my $sub_name = "dng_HashFromCommaSeparatedString";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($HR, $string) = @_;

  %{$HR} = ();
  my @key_A = split(",", $string);
  foreach my $key (@key_A) { 
    $HR->{$key} = 1; 
  }

  return;
}

#################################################################
# Subroutine: dng_BlastDbProteinCreate
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
sub dng_BlastDbProteinCreate {
  my $sub_name = "dng_BlastDbProteinCreate";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($makeblastdb, $fa_file, $opt_HHR, $FH_HR) = @_;

  utl_RunCommand($makeblastdb . " -in $fa_file -dbtype prot > /dev/null", opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}

#################################################################
# Subroutine: dng_FastaWriteSequence()
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
sub dng_FastaWriteSequence {
  my $sub_name = "dng_FastaWriteSequence";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($out_FH, $name, $def, $seq, $FH_HR) = @_;

  if(! defined $name) { ofile_FAIL("ERROR in $sub_name, name is undefined", "dnaorg", 1, $FH_HR); }
  if(! defined $seq)  { ofile_FAIL("ERROR in $sub_name, name is undefined", "dnaorg", 1, $FH_HR); }

  # capitalize and DNAize $seq
  dng_SqstringCapitalize(\$seq);
  dng_SqstringDnaize(\$seq);
  printf $out_FH (">%s%s\n%s", 
                  $name, 
                  (defined $def) ? " " . $def : "",
                  dng_SqstringAddNewlines($seq, 60));
  
  return;
}


#################################################################
# Subroutine: dng_FastaFileWriteFromStockholmFile()
# Incept:     EPN, Fri Mar 15 12:56:04 2019
#
# Synopsis: Use esl-reformat to convert a stockholm file to unaligned fasta
#
# Arguments:
#  $esl_reformat: esl-reformat executable file
#  $fa_file:      fasta file to create
#  $stk_file:     stockholm file
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's a problem fetching the sequence file
#################################################################
sub dng_FastaFileWriteFromStockholmFile { 
  my $sub_name = "dng_FastaFileWriteFromStockholmFile";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_reformat, $fa_file, $stk_file, $opt_HHR, $FH_HR) = @_;

  my $cmd = $esl_reformat . " --informat stockholm fasta $stk_file > $fa_file";
  utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # remove a .ssi file if it exists
  my $ssi_file = $fa_file . ".ssi";
  if(-e $ssi_file) { unlink $ssi_file; }

  return;
}

#################################################################
# Subroutine: dng_StockholmFileWriteFromFastaFile()
# Incept:     EPN, Sat Mar  9 10:24:14 2019
#
# Synopsis: Use esl-reformat to convert a fasta file to a stockholm file
#
# Arguments:
#  $esl_reformat: esl-reformat executable file
#  $fa_file:      fasta file
#  $stk_file:     stockholm file to create
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:        REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if there's a problem fetching the sequence file
#################################################################
sub dng_StockholmFileWriteFromFastaFile { 
  my $sub_name = "dng_StockholmFileWriteFromFastaFile";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_reformat, $fa_file, $stk_file, $opt_HHR, $FH_HR) = @_;

  my $cmd = $esl_reformat . " --informat afa stockholm $fa_file > $stk_file";
  utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
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
# Subroutine:  dng_CmalignCheckStdOutput()
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
sub dng_CmalignCheckStdOutput { 
  my $sub_name = "dng_CmalignCheckStdOutput";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($stdout_file, $ret_mxsize_R, $FH_HR) = @_;
  if(defined $ret_mxsize_R) { 
    $$ret_mxsize_R = 0; # overwritten below if nec
  }

  if(! -e $stdout_file) { 
    ofile_FAIL("ERROR in $sub_name, cmalign stdout file $stdout_file does not exist", "dnaorg", 1, $FH_HR);
  }
  if(! -s $stdout_file) { 
    ofile_FAIL("ERROR in $sub_name, cmalign stdout file $stdout_file exists but is empty", "dnaorg", 1, $FH_HR);
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
    if($error_line =~ m/\r$/) { chop $final_line; } # remove ^M if it exists
    if($error_line =~ /Error: HMM banded truncated alignment mxes need (\d+\.\d+)/) { 
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


###########################################################################
# the next line is critical, a perl module must return a true value
return 1;
###########################################################################
