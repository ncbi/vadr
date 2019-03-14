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
# - $seg_info_HAR: similar to ${ftr,seq,err}_info_HAR, except contains information pertaining 
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
#   getReferenceFeatureInfo()
#   fetchReferenceFeatureSequences()
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
#
# Subroutines related to the feature table error exception array of hashes:
#   initializeHardCodedFTableErrorExceptions()
#   addFTableErrorException()
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
#   dumpHashOfHashes()
#   dumpArrayOfHashesOfHashes()
#   dumpArrayOfHashes()
#
# Subroutines for validating the special data structures:
#   validateExecutableHash()
#   validateFeatureInfoHashIsComplete()
#   validateSegmentInfoHashIsComplete()
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
#   numberOfDigits()
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
#   determineFeatureTypes()
#   getNumFeaturesAnnotatedByModels()
#   getReferenceFeatureInfo()
#   fetchReferenceFeatureSequences()
#   featureHasChildren()
#   featureHasParent()
#


#################################################################
# Subroutine : featuresGetChildrenArrayOfArray()
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
sub featuresGetChildrenArrayOfArrays { 
  my $nargs_expected = 3;
  my $sub_name = "featuresGetChildrenArrayOfArrays";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_AHR, $AAR, $FH_HR) = @_;

  @{$AAR} = ();
  my $nftr = scalar(@{$ftr_info_AHR});

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    featureGetChildrenArray($ftr_info_AHR, $ftr_idx, $nftr, $AAR->[$ftr_idx], $FH_HR);
  }
  
  return;
}

#################################################################
# Subroutine : featureGetChildrenArray()
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
sub featureGetChildrenArray { 
  my $nargs_expected = 5;
  my $sub_name = "featureGetChildrenArray";
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
# Subroutine : featureGetTypeAndTypeIndexString()
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
sub featureGetTypeAndTypeIndexString { 
  my $nargs_expected = 3;
  my $sub_name = "featureGetChildrenArray";
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
  addToErrorInfoHash($err_info_HAR, "n_nst", "feature",  
                     "no in-frame stop codon exists 3' of predicted valid start codon", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_nm3", "feature",  
                     "length of nucleotide feature is not a multiple of 3", # description
                     "similar to !out_product,out_gene!; length is not a multiple of 3", # feature table note
                     "Unexpected Length: (!out_product,out_gene!) length is not a multiple of 3; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_stp", "feature",
                     "predicted CDS stop by homology is invalid; there may be a valid stop in a different location due to truncation (trc) or extension (ext) (TAG|TAA|TGA)", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified on !out_product,out_gene!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_trc", "feature",
                     "in-frame stop codon exists 5' of stop position predicted by homology to reference", # description
                     "similar to !out_product,out_gene!; contains premature stop codon", # feature table note
                     "CDS has Stop Codon: (!out_product,out_gene!) contains unexpected stop codon; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_ext", "feature",
                     "first in-frame stop codon exists 3' of stop position predicted by homology to reference", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_per", "feature",
                     "mat_peptide may not be translated because its CDS has a problem", # description
                     "", # feature table note
                     "Peptide Translation Problem: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_str", "feature",
                     "predicted CDS start position is not beginning of start codon", # description
                     "similar to !out_product,out_gene!; no start codon", # feature table note
                     "Mutation at Start: (!out_product,out_gene!) expected start codon could not be identified", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_gp5", "feature",
                     "alignment to reference is a gap at 5' boundary", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_gp3", "feature",
                     "alignment to reference is a gap at 3' boundary", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation at End: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_lp5", "feature",
                     "alignment to reference has low confidence at 5' boundary", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_lp3", "feature",
                     "alignment to reference has low confidence at 3' boundary", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation at End: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_non", "feature",
                     "protein-based search identifies CDS not identified in nucleotide-based search", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_nop", "feature",
                     "nucleotide-based search identifies CDS not identified in protein-based search", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_cst", "feature",
                     "strand mismatch between protein and nucleotide predictions", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_p5l", "feature",
                     "protein alignment extends past nucleotide alignment at 5' end", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_p5s", "feature",
                     "protein alignment does not extend close enough to nucleotide alignment 5' endpoint", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_p3l", "feature",
                     "protein alignment extends past nucleotide alignment at 3' end", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation at End: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_p3s", "feature",
                     "protein alignment does not extend close enough to nucleotide alignment 3' endpoint", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Indefinite Annotation at End: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "p_lin", "feature",
                     "too large of an insertion in protein alignment", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Insertion of Nucleotides: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "p_lde", "feature",
                     "too large of a deletion in protein alignment", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "Deletion of Nucleotides: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "p_trc", "feature",
                     "stop codon in protein alignment", # description
                     "similar to !out_product,out_gene!", # feature table note
                     "CDS has Stop Codon: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_zft", "sequence",
                     "zero features annotated", # description
                     "", # feature table note, irrelevant for per-sequence errors
                     "No Features Annotated: (*sequence*) zero annotated features", # feature table err
                     $FH_HR); 

  addToErrorInfoHash($err_info_HAR, "n_div", "sequence",
                     "sequence too distant from reference to annotate", # description
                     "", # feature table note, irrelevant for per-sequence errors
                     "Unexpected Divergence: (*sequence*) sequence is too divergent to confidently assign nucleotide-based annotation !DESC!", # feature table err
                     $FH_HR); 

  # define the ftbl_invalid_by values, these are one-sided, any error code listed in the 
  # 3rd argument invalidates the 2nd argument error code, but not vice versa

  # trc, ext and nst are preferred to stp
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "n_stp", "n_trc,n_ext,n_nst", $FH_HR); 

  # n_div is preferred to zft
  setFTableInvalidatedByErrorInfoHash($err_info_HAR, "b_zft", "n_div", $FH_HR);

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
#   $desc:            the error description/message of the element we are adding
#   $ftbl_note:       note message for feature table
#                     must eq "" if per-type is "sequence"
#                     must ne "" if per-type is "feature"
#   $ftbl_err:        ERROR message for feature table
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
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($err_info_HAR, $code, $pertype, $desc, 
      $ftbl_note, $ftbl_err,
      $FH_HR) = (@_);

  # make sure $pertype is valid
  if(($pertype ne "feature") && ($pertype ne "sequence")) { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add code $code with per-type $pertype that is not neither \"feature\" nor \"sequence\".", 1, $FH_HR); 
  }
  
  # make sure $ftbl_err is valid
  if((! defined $ftbl_err) || $ftbl_err eq "") { 
    DNAORG_FAIL("ERROR in $sub_name, trying to add code $code but ftbl_err is undefined or empty", 1, $FH_HR);
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
  if(! exists $err_info_HAR->{"desc"})            { @{$err_info_HAR->{"desc"}}            = (); }
  if(! exists $err_info_HAR->{"ftbl_invalid_by"}) { @{$err_info_HAR->{"ftbl_invalid_by"}} = (); }
  if(! exists $err_info_HAR->{"ftbl_note"})       { @{$err_info_HAR->{"ftbl_note"}}       = (); }
  if(! exists $err_info_HAR->{"ftbl_err"})        { @{$err_info_HAR->{"ftbl_err"}}        = (); }
  
  push(@{$err_info_HAR->{"code"}},             $code); 
  push(@{$err_info_HAR->{"pertype"}},          $pertype); 
  push(@{$err_info_HAR->{"desc"}},             $desc); 
  push(@{$err_info_HAR->{"ftbl_invalid_by"}},  ""); # initialized to no invalid_by's, possibly added to later with setFTableInvalidatedByErrorInfoHash()
  push(@{$err_info_HAR->{"ftbl_note"}},        $ftbl_note);
  push(@{$err_info_HAR->{"ftbl_err"}},         $ftbl_err);

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
#             possibly corrected one.  
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
# Returns: void
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
  my $valid = 0;
  for(my $e = 0; $e < $nerr; $e++) { 
    $err_idx = $err_idx_A[$e];
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
      # valid error that will impact output of feature table
      # add notes and errors
      my $note_str = populateFTableNoteOrError("ftbl_note", $err_idx, $seq_name, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, undef, $FH_HR);
      if($note_str ne "") { 
        push(@tmp_note_A, $note_str); # we will prune this array and populate @{$ret_note_AR} before returning
      }
      
      my $error_str = populateFTableNoteOrError("ftbl_err", $err_idx, $seq_name, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, undef, $FH_HR);
      if($error_str ne "") { 
        # only add the error, if an identical error does not already exist in @{$ret_error_AR}
        my $idx = findNonNumericValueInArray($ret_error_AR, $error_str, $FH_HR);
        if($idx == -1) { 
          push(@{$ret_error_AR}, $error_str); 
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

  return;
}

#################################################################
# Subroutine: processSequenceErrorsForFTable()
# Incept:     EPN, Thu Jan 24 12:09:24 2019
#
# Purpose:    Given a string of per-sequence errors that correspond
#             to a specific sequence, use the %{$err_info_HAR} and
#             process that string to determine what (if any) 
#             errors should be added to the feature table
#             for this sequence. Note that we do not add any 'notes'
#             as we possibly could in processFeatureErrorsForFTable() 
#             because we are dealing with the full sequence and not
#             a feature for a sequence.
#
# Arguments:
#   $err_code_str:           string of errors, comma separated, can be ""
#   $seq_name:               name of sequence
#   $err_info_HAR:           REF to hash of arrays with information on the errors, PRE-FILLED
#   $err_seq_instances_HHR:  REF to 2D hashes with per-sequence errors, PRE-FILLED
#   $ret_error_AR:           REF to array of errors, possibly added to here (not created)
#   $FH_HR:                  REF to hash of file handles, including "log" and "cmd"
# 
# Returns: void
#
# Dies: Never
#################################################################
sub processSequenceErrorsForFTable { 
  my $sub_name = "processSequenceErrorsForFTable";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($err_code_str, $seq_name, $err_info_HAR, $err_seq_instances_HHR, $ret_error_AR, $FH_HR) = (@_);

  if($err_code_str eq "") { 
    return 0; 
  }

  #printf("HEYA in $sub_name $seq_name, $err_code_str\n");

  # NOTE: there's some code duplication in this sub with
  # processFeatureErrorsForFtable(), possibly a chance for additional
  # subroutines

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

  my $nerr  = scalar(@err_idx_A);
  my $valid = 0;
  for(my $e = 0; $e < $nerr; $e++) { 
    $err_idx = $err_idx_A[$e];
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
      # add errors
      my $error_str = populateFTableNoteOrError("ftbl_err", $err_idx, $seq_name, -1, undef, $err_info_HAR, undef, $err_seq_instances_HHR, $FH_HR);
      if($error_str ne "") { 
        # only add the error, if an identical error does not already exist in @{$ret_error_AR}
        my $idx = findNonNumericValueInArray($ret_error_AR, $error_str, $FH_HR);
        if($idx == -1) { 
          push(@{$ret_error_AR}, $error_str); 
        }
      }
    }
  }

  return;
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
#   $ftr_idx:                feature index, -1 if this is a per-sequence error
#   $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
#                            must be undefined if $ftr_idx == -1
#                            must be defined   if $ftr_idx != -1
#   $err_info_HAR:           REF to hash of arrays with information on the errors, PRE-FILLED
#   $err_ftr_instances_AHHR: REF to array of 2D hashes with per-feature errors, PRE-FILLED
#                            must be undefined if $ftr_idx == -1
#                            must be defined   if $ftr_idx != -1
#   $err_seq_instances_HHR:  REF to array of 2D hashes with per-feature errors, PRE-FILLED
#                            must be undefined if $ftr_idx != -1
#                            must be defined   if $ftr_idx == -1
#   $FH_HR:                  REF to hash of file handles, including "log" 
#                            and "cmd"
# 
# Returns: string with the feature table note for the current sequence/feature combo
#
# Dies:    If ftbl_err_exceptions_HR doesn't have information we need
#          or has invalid information
#
#################################################################
sub populateFTableNoteOrError { 
  my $sub_name = "populateFTableNoteOrError";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ekey, $err_idx, $seq_name, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, $err_seq_instances_HHR, $FH_HR) = (@_);

  if(! exists $err_info_HAR->{$ekey}) { 
    DNAORG_FAIL("ERROR in $sub_name, $ekey value is undefined in error info hash", 1, $FH_HR);
  }
  # check that combination of $ftr_idx and $err_ftr_instances_AHHR and $err_seq_instances_HHR is valid
  if($ftr_idx != -1 && (! defined $err_ftr_instances_AHHR)) { 
    DNAORG_FAIL("ERROR in $sub_name, ftr_idx is not -1 but err_ftr_instances_AHHR is not defined", 1, $FH_HR);
  }
  if($ftr_idx == -1 && (defined $err_ftr_instances_AHHR)) { 
    DNAORG_FAIL("ERROR in $sub_name, ftr_idx is -1 but err_ftr_instances_AHHR is defined", 1, $FH_HR);
  }
  if($ftr_idx != -1 && (! defined $ftr_info_HAR)) { 
    DNAORG_FAIL("ERROR in $sub_name, ftr_idx is not -1 but ftr_info_HAR is not defined", 1, $FH_HR);
  }
  if($ftr_idx == -1 && (defined $ftr_info_HAR)) { 
    DNAORG_FAIL("ERROR in $sub_name, ftr_idx is -1 but ftr_info_HAR is defined", 1, $FH_HR);
  }
  if($ftr_idx == -1 && (! defined $err_seq_instances_HHR)) { 
    DNAORG_FAIL("ERROR in $sub_name, ftr_idx is -1 but err_seq_instances_AHHR is not defined", 1, $FH_HR);
  }
  if($ftr_idx != -1 && (defined $err_seq_instances_HHR)) { 
    DNAORG_FAIL("ERROR in $sub_name, ftr_idx is not -1 but err_ftr_instances_AHHR is defined", 1, $FH_HR);
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
    if(($ftr_idx != -1) && (exists $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name})) { 
      my $desc_str = sprintf("%s%s", 
                             $err_info_HAR->{"desc"}[$err_idx], 
                             ($err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} eq "") ? "" : " [" . $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} . "]"); 
      $ret_msg =~ s/!DESC!/$desc_str/g;
    }
    elsif(($ftr_idx == -1) && (exists $err_seq_instances_HHR->{$err_code}{$seq_name})) { 
      my $desc_str = sprintf("%s%s", 
                             $err_info_HAR->{"desc"}[$err_idx], 
                             ($err_seq_instances_HHR->{$err_code}{$seq_name} eq "") ? "" : " [" . $err_seq_instances_HHR->{$err_code}{$seq_name} . "]"); 
      $ret_msg =~ s/!DESC!/$desc_str/g;
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, trying to return $ekey message for $err_code and sequence $seq_name feature $ftr_idx, but no error instance exists", 1, $FH_HR);
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
#   $strands_AR:  ref to array to fill with strands of each segment, can be undef
#   $nsegments_R: ref to scalar that fill with the number of segments
#   $FH_HR:       REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void; but fills @{$starts_AR}, @{$stops_AR}, and $$nsegments_R.
#
# Dies:       if we see a feature that spans stop..start but $do_circular is 0
#################################################################
sub startsStopsStrandsFromCoordsLength { 
  my $sub_name = "startsStopsStrandsFromCoordsLength()";
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
      DNAORG_FAIL("ERROR in $sub_name: found internal complement in coords string $coords, we assume all segments are on same strand...", 1, $FH_HR); 
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
      $$nsegments_R++;
    }
    elsif($el =~ m/^(\d+)$/) { # a single nucleotide
      push(@{$starts_AR}, $1);
      push(@{$stops_AR},  $1);
      if(defined $strands_AR) { push(@{$strands_AR}, $cur_strand); }
      $$nsegments_R++;
    }
    else { 
      DNAORG_FAIL("ERROR unable to parse $orig_coords in $sub_name", 1, $FH_HR); 
    }
  }

  # check if we have a spanning segment (that spans stop..start) and if we do
  # and (! $do_circular) then die, because that shouldn't happen.
  my $have_spanning_segment = checkForSpanningSequenceSegments($starts_AR, $stops_AR, $nsegments_R, 0, $strand, $totlen); # 1 says: do correct the spanning segment
  if($have_spanning_segment) { 
    if(! $do_circular) { 
      DNAORG_FAIL("ERROR in $sub_name, found segment that spanned stop..start boundary, but we're not allowing circular genomes...", 1, $FH_HR); 
    }
    else { 
      # fix it
      checkForSpanningSequenceSegments($starts_AR, $stops_AR, $nsegments_R, 1, $strand, $totlen); # 0 says: don't correct the spanning segment
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
sub startsStopsFromCoords { 
  my $sub_name = "startsStopsFromCoords()";
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
#   dumpHashOfHashes()
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
# Subroutine: dumpHashOfHashes()
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
sub dumpHashOfHashes { 
  my $sub_name = "dumpHashOfHashes()";
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
# Subroutine: arrayOfHashesKeepKeys()
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
sub arrayOfHashesKeepKeys { 
  my $sub_name = "arrayOfHashesKeepKeys";
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
#   validateErrorInfoHashIsComplete()
#   validateInfoHashOfArraysIsComplete()
#   validateOutputFileInfoHashOfHashes()
#   validateAndGetSizeOfInfoHashOfArrays()
#   getConsistentSizeOfInfoHashOfArrays()
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
# Subroutine: validateSequenceInfoHashIsComplete()
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
sub validateSequenceInfoHashIsComplete { 
  my $sub_name = "validateSequenceInfoHashIsComplete()";
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
#                "code":            the error code, e.g. "b5e"
#                "pertype":         type of this error, either "feature" (error occurs for a single feature) or "sequence"
#                                   (error occurs for an entire sequence)
#                "desc":            the description/message for this error, e.g. ""unable to identify homologous feature"
#             Information relevant  to feature table output:
#                "ftbl_invalid_by": string that lists all error code indices (separated by commas) the current error code is invalidated by
#                                   for feature table output (if any of those errors are also present, this error becomes invalid and does
#                                   not impact feature table output)
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
  
  my @expected_keys_A = ("code", "pertype", "desc", "ftbl_invalid_by", "ftbl_note", "ftbl_err");

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
# Subroutine: getSizeOfInfoHashOfArrays()
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
sub getSizeOfInfoHashOfArrays { 
  my $sub_name = "getSizeOfInfoHashOfArrays()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($info_HAR, $key, $FH_HR) = (@_);
  
  if(! defined $info_HAR) { 
    DNAORG_FAIL("ERROR in $sub_name, input info_HAR undefined", 1, $FH_HR);
  }
  if(! exists $info_HAR->{$key}) { 
    DNAORG_FAIL("ERROR in $sub_name, key $key does not exist in info_HAR hash", 1, $FH_HR);
  }

  return(scalar(@{$info_HAR->{$key}}));
}

#################################################################
# Subroutine: validateCapitalizedDnaStopCodon()
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
sub validateCapitalizedDnaStopCodon {
  my $sub_name = "validateCapitaliedDnaStopCodon";
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
# Subroutine: validateCapitalizedDnaStartCodon()
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
sub validateCapitalizedDnaStartCodon {
  my $sub_name = "validateCapitaliedDnaStartCodon";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($codon) = @_;
  
  if($codon eq "ATG") { 
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
#   numberOfDigits()
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
# Subroutine : numberOfDigits()
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
sub numberOfDigits { 
    my $nargs_expected = 1;
    my $sub_name = "numberOfDigits()";
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

    my ($num) = (@_);

    my $ndig = 1; 
    while($num >= 10) { $ndig++; $num /= 10.; }

    return $ndig;
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
#   $do_failok:   '1' to NOT exit if command fails, '0' to exit if command fails
#   $FH_HR:       REF to hash of file handles, including "cmd"
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails and $do_failok is '0'
#################################################################
sub runCommand {
  my $sub_name = "runCommand()";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $be_verbose, $do_failok, $FH_HR) = @_;
  
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

  if(($? != 0) && (! $do_failok)) { 
    DNAORG_FAIL("ERROR in $sub_name, the following command failed:\n$cmd\n", $?, $FH_HR); 
  }

  return ($stop_time - $start_time);
}

#################################################################
# Subroutine:  submitJob()
# Incept:      EPN, Wed Feb  6 12:35:04 2019
#
# Purpose:     Submits a job to sge.
#
# Arguments:
#   $cmd:            command to run
#   $job_name:       name for job
#   $err_file:       name of err file to create, can be "/dev/null"
#   $mem_gb:         number of Gb of memory required
#   $nsecs:          maximum number of seconds to allow jobs to take
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#   $ofile_info_HHR: REF to the 2D hash of output file information, ADDED TO HERE 
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if qsub $cmd fails
#################################################################
sub submitJob {
  my $sub_name = "submitJob()";
  my $nargs_expected = 7;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $job_name, $err_file, $mem_gb, $nsecs, $opt_HHR, $ofile_info_HHR) = @_;
  
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  if(($err_file ne "/dev/null") && (-e $err_file)) { 
    removeFileUsingSystemRm($err_file, $sub_name, $opt_HHR, $ofile_info_HHR); 
  }
  my $submit_cmd = sprintf("qsub -N $job_name -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $err_file -m n -l h_rt=%d,h_vmem=%dG,mem_free=%dG,reserve_mem=%dG,m_mem_free=%dG " . "\"" . $cmd . "\" > /dev/null\n", $nsecs, $mem_gb, $mem_gb, $mem_gb, $mem_gb);

  runCommand($submit_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
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

  runCommand("rm $file", opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}


#################################################################
# Subroutine : removeListOfFiles()
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
sub removeListOfFiles { 
  my $nargs_expected = 4;
  my $sub_name = "removeListOfFiles()";
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
    runCommand($rm_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
    $i = $up;
  }
  
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
# Subroutine : arrayToString()
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
sub arrayToString { 
  my $nargs_expected = 2;
  my $sub_name = "arrayToString";
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
    runCommand($cat_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

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
  runCommand("md5sum $file > $out_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);

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
#  $outfile_AR:      ref to array of output files that will be created by jobs we are waiting for
#  $errfile_AR:      ref to array of err files that will be created by jobs we are waiting for if 
#                    any stderr output is created
#  $success_AR:      ref to array of success values, FILLED HERE
#                    these will always all be '1' unless $do_cmalign
#                    if($do_cmalign) some may be '0'
#  $mxsize_AR:       ref to array of required matrix sizessuccess values, CAN BE UNDEF
#                    $mxsize_AR->[$j] set to value readh from cmalign output, if $success_AR->[$j] == 0
#                    else set to '0'
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
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($do_cmalign, $outfile_AR, $errfile_AR, $success_AR, $mxsize_AR, $finished_str, $nmin, $do_errcheck, $FH_HR) = @_;

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
          if($do_cmalign) { 
            my $success = cmalignCheckStdOutput($outfile_AR->[$i], 
                                                (defined $mxsize_AR) ? \$mxsize_AR->[$i] : undef,
                                                $FH_HR);
            if($success == 0 || $success == 1) { 
              $success_AR->[$i] = $success;
              $is_finished_A[$i] = 1;
              $nfinished++;
            }
          }
          else { 
            my $final_line = `tail -n 1 $outfile_AR->[$i]`;
            chomp $final_line;
            if($final_line =~ m/\r$/) { chop $final_line; } # remove ^M if it exists
            if($final_line =~ m/\Q$finished_str\E/) { 
              $success_AR->[$i] = 1; # if we're not running cmalign, if we see $finished_str, job was successful
              $is_finished_A[$i] = 1;
              $nfinished++;
            }
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
sub splitFastaFile { 
  my $sub_name = "splitFastaFile()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_ssplit, $fasta_file, $nfiles, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $outfile = $fasta_file . ".esl-ssplit";
  my $cmd = undef;
  if($nfiles == -1) { # special case: put 1 file per sequence
    $cmd = "$esl_ssplit -v $fasta_file 1 > $outfile";
  }
  else { 
    $cmd = "$esl_ssplit -v -r -n $fasta_file $nfiles > $outfile";
  }
  runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  # parse output to determine exactly how many files were created:
  # $esl_ssplit will have output exactly 1 line per fasta file it created
  my $nfiles_created = countLinesInFile($outfile, $FH_HR);

  if(! opt_Get("--keep", $opt_HHR)) { 
    runCommand("rm $outfile", opt_Get("-v", $opt_HHR), 0, $FH_HR);
  }

  return $nfiles_created;
}

#################################################################
# Subroutine:  cmalignOrNhmmscanWrapper()
# Incept:      EPN, Tue Nov 21 15:54:20 2017
#              EPN, Tue Jan 29 14:30:05 2019 [switched from cmscan to cmalign]
#
# Purpose:     Run one or more cmalign or nhmmscan jobs on the farm
#              or locally, after possibly splitting up the input
#              sequence file by calling 
#              cmalignOrNhmmscanWrapperHelper(). 
#
#              If we're doing cmalign we may have to do two rounds
#              of sequence file splitting and job running/submission
#              because there is an error case in cmalign that we want
#              to be able to detect. That error case is when the sequence
#              requires too much memory to align. In order to catch 
#              those error cases we need to run each offending sequence
#              individually, so our strategy for cmalign is:
#
#              Split full fasta file up using default method and run
#              >= 1 cmalign jobs. If any of those runs R fail, then 
#              split up run R's sequence file into >= 1 files with
#              exactly 1 sequence in them. One or more of those should
#              fail and that reveals which specific sequences are
#              causing the memory overflow.
# 
#              If we are doing nhmmscan then we expect that all 
#              jobs/runs will pass and so the above complications
#              are not necessary.
#
#              The following must all be valid options in opt_HHR:
#              --nseq, --local, --wait, --keep, -v
#              See dnaorg_annotate.pl for examples of these options.
#
# Arguments: 
#  $execs_HR:              ref to executables with "esl-ssplit" and either "cmalign" or "nhmmscan" 
#                          defined as keys
#  $do_cmalign:            '1' if we are running cmalign, '0' if we are running nhmmscan
#  $out_root:              string for naming output files
#  $seq_file:              name of sequence file with all sequences to run against
#  $tot_len_nt:            total length of all nucleotides in $seq_file
#  $progress_w:            width for outputProgressPrior output
#  $mdl_filename:          name of model file to use
#  $stk_file_AR:           ref to array of stockholm files created here, FILLED HERE if $do_cmalign
#  $overflow_seq_AR:       ref to array of sequences that failed due to matrix overflows, FILLED HERE if $do_cmalign
#  $overflow_mxsize_AR:    ref to array of required matrix sizes for each sequence that failed due to matrix overflows, FILLED HERE if $do_cmalign
#  $opt_HHR:               REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:        REF to 2D hash of output file information
#
# Returns:     void, updates $$nfa_created_R with number of
#              fasta files created.
# 
# Dies: If an executable doesn't exist, or cmalign or nhmmscan or esl-ssplit
#       command fails if we're running locally
################################################################# 
sub cmalignOrNhmmscanWrapper { 
  my $sub_name = "cmalignOrNhmmscanWrapper";
  my $nargs_expected = 12;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $do_cmalign, $out_root, $seq_file, $tot_len_nt, $progress_w,
      $mdl_filename, $stk_file_AR, $overflow_seq_AR, $overflow_mxsize_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $program_choice = ($do_cmalign) ? "cmalign" : "nhmmscan";
  my $nfasta_created = 0; # number of fasta files created by esl-ssplit
  my $log_FH = $ofile_info_HHR->{"FH"}{"log"}; # for convenience
  my $start_secs; # timing start
  my $do_local = opt_Get("--local", $opt_HHR);
  @{$overflow_seq_AR} = (); # we will fill this with names of sequences that fail cmalign because
                            # the matrix required to align them is too big

  # set up output file names
  my @concat_keys_A = (); # %r{1,2}_out_file_HAR keys we are going to concatenate files for
  my %concat_HA = ();     # hash of arrays of all files to concatenate together
  my $out_key;            # key for an output file: e.g. "stdout", "ifile", "tfile", "tblout", "err"
  if($do_cmalign) { 
    @concat_keys_A = ("stdout", "ifile", "tfile"); 
  }
  else { 
    @concat_keys_A = ("tblout");
  }
  if(! $do_local) { push(@concat_keys_A, "err"); }
  foreach my $out_key (@concat_keys_A) { 
    @{$concat_HA{$out_key}} = ();
  }    

  # determine how many sequence files we need to split $seq_file into to satisfy <n> sequences
  # per job (where n is from --nseq <n>). If more than 1 and --local not used, we split up 
  # the sequence file and submit jobs to farm, if only 1 or --local used, we just run it locally
  my $targ_nseqfiles = int($tot_len_nt / (opt_Get("--nkb", $opt_HHR) * 1000)); 
  # int() takes the floor, so there can be a nonzero remainder. We don't add 1 though, 
  # because splitFastaFile() will return the actual number of sequence files created
  # and we'll use that as the number of jobs subsequently. $targ_nseqfiles is currently only
  # the 'target' number of sequence files that we pass into splitFastaFile().
  # make sure we won't exceed our max number of jobs (from --maxnjobs)
  if(($targ_nseqfiles) > (opt_Get("--maxnjobs", $opt_HHR))) { 
    $targ_nseqfiles = int(opt_Get("--maxnjobs", $opt_HHR));
  }
  if($targ_nseqfiles == 0) { $targ_nseqfiles = 1; }

  # we need to split up the sequence file, and submit a separate set of nhmmscan/cmalign jobs (one per model file) for each
  my $nr1 = 0; # number of runs in round 1 (one per sequence file we create)
  my %r1_out_file_HA = (); # hash of arrays ([0..$nr1-1]) of output files for cmalign/nhmmscan round 1 runs
  my @r1_success_A   = (); # [0..$nr1-1]: '1' if this run finishes successfully, '0' if not
  my @r1_mxsize_A    = (); # [0..$nr1-1]: if $r1_success_A[$r1_i] is '0', required size for dp mx, else '0'
  my @r1_seq_file_A  = (); # [0..$nr1-1]: name of sequence file for this run
  my $r1_i;                 # counter over round 1 runs
  if($targ_nseqfiles == 1) { 
    my $cp_command = "cp " . $seq_file . " " . $seq_file . ".1";
    runCommand($cp_command, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
    $nr1 = 1;
  }
  else { # we are going to split up the fasta file 
    $nr1 = splitFastaFile($execs_HR->{"esl-ssplit"}, $seq_file, $targ_nseqfiles, $opt_HHR, $ofile_info_HHR);
    # splitFastaFile will return the actual number of fasta files created, 
    # which can differ from the requested amount (which is $targ_nseqfiles) that we pass in. 
  }
  for($r1_i = 0; $r1_i < $nr1; $r1_i++) { # update sequence file names
    $r1_seq_file_A[$r1_i] = $seq_file . "." . ($r1_i+1);
  }
  
  cmalignOrNhmmscanWrapperHelper($execs_HR, $do_cmalign, $out_root, $progress_w, 
                                 \@r1_seq_file_A, \%r1_out_file_HA, \@r1_success_A, \@r1_mxsize_A,
                                 $mdl_filename, $opt_HHR, $ofile_info_HHR); 
  
  my $nr2            = 0;  # number of round 2 runs (sequence files)
  my %r2_out_file_HA = (); # hash of arrays ([0..$nr2-1]) of output files for cmalign/nhmmscan round 2 runs
  my @r2_success_A   = (); # [0..$nr2-1]: '1' if this run finishes successfully, '0' if not
  my @r2_mxsize_A    = (); # [0..$nr2-1]: if $r2_success_A[$r2_i] is '0', required size for dp mx, else '0'
  my @r2_seq_file_A  = (); # [0..$nr2-1]: name of sequence file for this run
  my $r2_i;                # counter over round 2 runs

  # go through each run:
  # if it finished successfully record its output files to concatenate later
  # if it did not finish successfully rerun all of its sequences (if $do_cmalign)
  for($r1_i = 0; $r1_i < $nr1; $r1_i++) { 
    if($r1_success_A[$r1_i]) { 
      # run finished successfully
      foreach $out_key (@concat_keys_A) { 
        push(@{$concat_HA{$out_key}}, $r1_out_file_HA{$out_key}[$r1_i]);
      }
      if($do_cmalign) { 
        push(@{$stk_file_AR}, $r1_out_file_HA{"stk"}[$r1_i]);
      }
    }
    else { 
      # run did not finish successfully
      if(! $do_cmalign) { 
        DNAORG_FAIL("ERROR in $sub_name a nhmmscan job failed.", 1, $ofile_info_HHR->{"FH"});
      }
      # if we get here, we know that $do_cmalign is 1
      # split this sequence file up into multiple files with only 1 sequence each, 
      # remember which 
      my $cur_nr2 = splitFastaFile($execs_HR->{"esl-ssplit"}, $r1_seq_file_A[$r1_i], -1, $opt_HHR, $ofile_info_HHR);
      if($cur_nr2 == 1) { 
        # special case, r1 sequence file had only 1 sequence, so we know the culprit
        # and don't need to rerun cmalign
        cmalignStoreOverflow($r1_seq_file_A[$r1_i], $r1_mxsize_A[$r1_i], $overflow_seq_AR, $overflow_mxsize_AR, $ofile_info_HHR->{"FH"}); 
        if(! opt_Get("--keep", $opt_HHR)) { 
          runCommand("rm " . $r1_seq_file_A[$r1_i] . ".1", opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
        }
      }
      else { 
        # r1 sequence file had > 1 sequence, we need to run each sequence independently through cmalign
        for($r2_i = 0; $r2_i < $cur_nr2; $r2_i++) { 
          push(@r2_seq_file_A, $r1_seq_file_A[$r1_i] . "." . ($r2_i+1));
        }
        $nr2 += $cur_nr2;
      }
    }
  }

  # do all round 2 runs
  if($nr2 > 0) { 
    cmalignOrNhmmscanWrapperHelper($execs_HR, $do_cmalign, 
                                   $out_root . ".r2",
                                   -1 * $progress_w,          # passing in $progress_w < 0 is our flag that we are 'rerunning'
                                   \@r2_seq_file_A, \%r2_out_file_HA, \@r2_success_A, \@r2_mxsize_A, 
                                   $mdl_filename, $opt_HHR, $ofile_info_HHR);
    # go through all round 2 runs: 
    # if it finished successfully record its output files to concatenate later
    # if it did not finish successfully, record the name of the sequence and mxsize required
    for($r2_i = 0; $r2_i < $nr2; $r2_i++) { 
      if($r2_success_A[$r2_i]) { 
        # run finished successfully
        foreach my $out_key (@concat_keys_A) { 
          push(@{$concat_HA{$out_key}}, $r2_out_file_HA{$out_key}[$r2_i]);
        }
        push(@{$stk_file_AR}, $r2_out_file_HA{"stk"}[$r2_i]);
      }
      else { 
        # run did not finish successfully
        cmalignStoreOverflow($r2_seq_file_A[$r2_i], $r2_mxsize_A[$r2_i], $overflow_seq_AR, $overflow_mxsize_AR, $ofile_info_HHR->{"FH"}); 
      }
    }
    # remove sequence files if --keep not used
    if(! opt_Get("--keep", $opt_HHR)) { 
      removeListOfFiles(\@r2_seq_file_A, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
    }
  }
    
  # concatenate all the stdout/ifile/tblout files into one 
  foreach $out_key (@concat_keys_A) { 
    my $concat_file = $out_root . "." . $program_choice . "." . $out_key;
    concatenateListOfFiles($concat_HA{$out_key}, $concat_file, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
    # concatenateListOfFiles() removes individual files unless --keep enabled
    addClosedFileToOutputInfo($ofile_info_HHR, "concat." . $out_key, $concat_file, 0, "concatenated $out_key file");
  }
  # remove sequence files 
  if(! opt_Get("--keep", $opt_HHR)) { 
    removeListOfFiles(\@r1_seq_file_A, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
  }

  return;
}

#################################################################
# Subroutine:  cmalignOrNhmmscanWrapperHelper()
# Incept:      EPN, Tue Feb  5 09:15:49 2019
#
# Purpose:     Run one or more cmalign or nhmmscan jobs on the farm
#              or locally.
#
#              Helper subroutine for cmalignOrNhmmscanWrapper()
#              see that sub's "Purpose" for more details.
#
#              The following must all be valid options in opt_HHR:
#              --nseq, --local, --wait, --keep, -v
#              See dnaorg_annotate.pl for examples of these options.
#
# Arguments: 
#  $execs_HR:              ref to executables with "esl-ssplit" and either "cmalign" or "nhmmscan" 
#                          defined as keys
#  $do_cmalign:            '1' if we are running cmalign, '0' if we are running nhmmscan
#  $out_root:              string for naming output files
#  $progress_w:            width for outputProgressPrior output, -1 * $progress_w if we are calling this to rerun cmalign again
#  $seq_file_AR:           ref to array of sequence file names for each cmalign/nhmmscan call, PRE-FILLED
#  $out_file_HAR:          ref to hash of arrays of output file names, FILLED HERE 
#  $success_AR:            ref to array of success values, 
#                          $success_AR->[$j] set to '1' if job finishes successfully
#                                            set to '0' if job fails due to mx overflow (must be cmalign)
#  $mxsize_AR:             ref to array of required matrix sizessuccess values, CAN BE UNDEF
#                          $mxsize_AR->[$j] set to value readh from cmalign output, if $success_AR->[$j] == 0
#                                           else set to '0'
#  $mdl_filename:          name of model file to use
#  $opt_HHR:               REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:        REF to 2D hash of output file information
#
# Returns:     void
# 
# Dies: If an executable doesn't exist, or command fails (and its not a cmalign allowed failure)
#
################################################################# 
sub cmalignOrNhmmscanWrapperHelper { 
  my $sub_name = "cmalignOrNhmmscanWrapperHelper";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $do_cmalign, $out_root, $progress_w, $seq_file_AR, $out_file_HAR, $success_AR, $mxsize_AR, $mdl_filename, $opt_HHR, $ofile_info_HHR) = @_;

  my $program_choice = ($do_cmalign) ? "cmalign" : "nhmmscan";
  my $do_local       = opt_Get("--local", $opt_HHR) ? 1 : 0;
  my $nfasta_created = 0;     # number of fasta files created by esl-ssplit
  my $log_FH         = $ofile_info_HHR->{"FH"}{"log"}; # for convenience
  my $nseq_files     = scalar(@{$seq_file_AR});

  my $desc = "";
  if($do_local) { 
    if($progress_w < 0) { # flag for 'rerunning'
      $desc = "Rerunning $program_choice locally to find sequences that are too divergent to align";
      $progress_w *= -1;
    }
    else { 
      $desc = "Running $program_choice locally";
    }
  }      
  else { 
    if($progress_w < 0) { # flag for 'rerunning'
      $desc = "Resubmitting $nseq_files $program_choice jobs to the farm to find seqs too divergent to align";
      $progress_w *= -1;
    }
    else { 
      $desc = "Submitting $nseq_files $program_choice jobs to the farm"
    }
  }
  my $start_secs = outputProgressPrior($desc, $progress_w, $log_FH, *STDOUT);

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
      $success_AR->[$s] = runCmalign($execs_HR->{"$program_choice"}, $mdl_filename, $seq_file_AR->[$s], 
                                     $out_file_HAR->{"stdout"}[$s], $out_file_HAR->{"ifile"}[$s], $out_file_HAR->{"tfile"}[$s], 
                                     $out_file_HAR->{"stk"}[$s], $out_file_HAR->{"err"}[$s],
                                     (defined $mxsize_AR) ? \$mxsize_AR->[$s] : undef, 
                                     $opt_HHR, $ofile_info_HHR);   
    }
    else { 
      runNhmmscan($execs_HR->{"$program_choice"}, $mdl_filename, $seq_file_AR->[$s], 
                  $out_file_HAR->{"stdout"}[$s], $out_file_HAR->{"tblout"}[$s], $out_file_HAR->{"err"}[$s],
                  $opt_HHR, $ofile_info_HHR);   
      $success_AR->[$s] = 1; # if runNhmmscan() returns the run was successful
    }
    # if we are not running local, ignore the return values from the run{Cmalign,Nhmmscan} subroutines
    if(! $do_local) { $success_AR->[$s] = 0; }
  }
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  if(! $do_local) { 
    if((opt_Exists("--skipalign", $opt_HHR)) && (opt_Get("--skipalign", $opt_HHR))) { 
      for($s = 0; $s < $nseq_files; $s++) { 
        $success_AR->[$s] = 1; 
      }
    }
    else { 
      # --skipalign not enabled
      # wait for the jobs to finish
      $start_secs = outputProgressPrior(sprintf("Waiting a maximum of %d minutes for all farm jobs to finish", opt_Get("--wait", $opt_HHR)), 
                                        $progress_w, $log_FH, *STDOUT);
      my $njobs_finished = waitForFarmJobsToFinish($do_cmalign, 
                                                   ($do_cmalign) ? $out_file_HAR->{"stdout"} : $out_file_HAR->{"tblout"}, 
                                                   $out_file_HAR->{"err"}, 
                                                   $success_AR,  
                                                   $mxsize_AR, # this may be undef
                                                   ($do_cmalign) ? "" : "[ok]", # value is irrelevant for cmalign
                                                   opt_Get("--wait", $opt_HHR), opt_Get("--errcheck", $opt_HHR), $ofile_info_HHR->{"FH"});
      if($njobs_finished != $nseq_files) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name only $njobs_finished of the $nseq_files are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", $opt_HHR)), 1, $ofile_info_HHR->{"FH"});
      }
      outputString($log_FH, 1, "# "); # necessary because waitForFarmJobsToFinish() creates lines that summarize wait time and so we need a '#' before 'done' printed by outputProgressComplete()
    }
  }
  
  return;
}

#################################################################
# Subroutine : runCmalign()
# Incept:      EPN, Wed Feb  6 12:30:08 2019
#
# Purpose:     Run Infernal's cmalign executable using $model_file
#              as the model file on sequence file $seq_file, either
#              locally or on the farm.
#              
#              If job does not finish successfully, we need to 
#              parse the stderr output (which we redirect to stdout)
#              and see if it failed because of a specific type of
#              error, because the required DP matrix exceeded the
#              size limit. That error looks like this:
#        
#              Error: HMM banded truncated alignment mxes need 60.75 Mb > 2.00 Mb limit.
#              Use --mxsize, --maxtau or --tau.
#              
# Arguments: 
#  $executable:       path to the cmscan or nhmmscan executable file
#  $model_file:       path to the CM file
#  $seq_file:         path to the sequence file
#  $stdout_file:      path to the stdout file to create
#  $ifile_file:       path to the cmalign --ifile file to create
#  $tfile_file:       path to the cmalign --tfile file to create
#  $stk_file:         path to the cmalign stockholm alignment output file to create
#  $err_file:         path to the error output file (only used if --local not used)
#  $ret_mxsize_R:     REF to required matrix size, only filled if return value is '0'
#  $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:   REF to 2D hash of output file information
# 
# Returns:     '1' upon success, which occurs if
#                  job run on farm and submission went ok
#                  job run locally and finished without error
#              '0' upon allowed failure, which occurs if
#                  job run locally and fails because of too big a required matrix
#
# Dies: upon unallowed failure, which occurs if
#                  job run on farm and submission failed
#                  job run locally and finished with unallowed failure
# 
################################################################# 
sub runCmalign { 
  my $sub_name = "runCmalign()";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($executable, $model_file, $seq_file, $stdout_file, $ifile_file, $tfile_file, $stk_file, $err_file, $ret_mxsize_R, $opt_HHR, $ofile_info_HHR) = @_;

  if(defined $ret_mxsize_R) { 
    $$ret_mxsize_R = 0; # overwritten below if nec
  }

  my $FH_HR    = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $do_local = opt_Get("--local", $opt_HHR) ? 1 : 0;

  if($stdout_file eq "/dev/null") { 
    DNAORG_FAIL("ERROR in $sub_name, stdout_file is /dev/null", 1, $FH_HR);
  }

  validateFileExistsAndIsNonEmpty($model_file, $sub_name, $FH_HR); 
  validateFileExistsAndIsNonEmpty($seq_file,   $sub_name, $FH_HR);

  # determine cmalign options based on command line options
  my $opts = sprintf(" --verbose --cpu 0 --ifile $ifile_file --tfile $tfile_file -o $stk_file --tau %s --mxsize %s", opt_Get("--tau", $opt_HHR), opt_Get("--mxsize", $opt_HHR));
  # add --sub and --notrunc unless --nosub used
  if(! opt_Get("--nosub", $opt_HHR)) { 
    $opts .= " --sub --notrunc"; 
  }
  # add -g unless --noglocal used
  if(! opt_Get("--noglocal", $opt_HHR)) { 
    $opts .= " -g"; 
  }
  if(! opt_Get("--nofixedtau", $opt_HHR)) { 
    $opts .= " --fixedtau"; 
  }
  
  # remove the tblout file or ifile file if it exists, this is important because we'll use the existence and
  # final line of this file to determine when the jobs are finished, if it already exists, we'll
  # think the job is finished before it actual is.
  if((! opt_Exists("--skipalign", $opt_HHR)) || (! opt_Get("--skipalign", $opt_HHR))) { 
    if(-e $stdout_file) { removeFileUsingSystemRm($stdout_file, $sub_name, $opt_HHR, $ofile_info_HHR); }
    if(-e $ifile_file)  { removeFileUsingSystemRm($ifile_file,  $sub_name, $opt_HHR, $ofile_info_HHR); }
    if(-e $tfile_file)  { removeFileUsingSystemRm($tfile_file,  $sub_name, $opt_HHR, $ofile_info_HHR); }
    if(-e $stk_file)    { removeFileUsingSystemRm($stk_file,    $sub_name, $opt_HHR, $ofile_info_HHR); }
    if(-e $err_file)    { removeFileUsingSystemRm($err_file,    $sub_name, $opt_HHR, $ofile_info_HHR); }
  }

  my $cmd = "$executable $opts $model_file $seq_file > $stdout_file 2>&1";

  if($do_local) { 
    if((! opt_Exists("--skipalign", $opt_HHR)) || (! opt_Get("--skipalign", $opt_HHR))) { 
      runCommand($cmd, opt_Get("-v", $opt_HHR), 1, $FH_HR); # 1 says: it's okay if job fails
    }
  }
  else { 
    my $job_name = "J" . removeDirPath($seq_file);
    my $nsecs  = opt_Get("--wait", $opt_HHR) * 60.;
    my $mem_gb = (opt_Get("--mxsize", $opt_HHR) / 1000.) * 2; # multiply --mxsize Gb by 2 to be safe
    if($mem_gb < 8.) { $mem_gb = 8.; } # set minimum of 8 Gb
    if((! opt_Exists("--skipalign", $opt_HHR)) || (! opt_Get("--skipalign", $opt_HHR))) { 
      submitJob($cmd, $job_name, $err_file, $mem_gb, $nsecs, $opt_HHR, $ofile_info_HHR);
    }
  }
  
  my $success = 1;
  if($do_local) { 
    # command has completed, check for the error in the stdout, or a final line of 'CPU' indicating that it worked.
    $success = cmalignCheckStdOutput($stdout_file, $ret_mxsize_R, $FH_HR);
    if($success == -1) { # indicates job did not finish properly, this shouldn't happen because runCommand() didn't die
      DNAORG_FAIL("ERROR in $sub_name, cmalign failed in a bad way, see $stdout_file for error output", 1, $ofile_info_HHR->{"FH"});
    }
  }
  
  return $success; 
}

#################################################################
# Subroutine : runNhmmscan()
# Incept:      EPN, Wed Feb  6 12:38:11 2019
#
# Purpose:     Run HMMER's nhmmscan executable using $model_file
#              as the model file on sequence file $seq_file, either
#              locally or on the farm.
#
# Arguments: 
#  $executable:       path to the cmscan or nhmmscan executable file
#  $model_file:       path to the CM file
#  $seq_file:         path to the sequence file
#  $stdout_file:      path to the stdout file to create, can be "/dev/null"
#  $tblout_file:      path to the nhmmscan --tblout file to create
#  $err_file:         path to the error output file (only used if --local not used)
#  $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:   REF to 2D hash of output file information
# 
# Returns:     void
# 
################################################################# 
sub runNhmmscan {
  my $sub_name = "runNhmmscan()";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($executable, $model_file, $seq_file, $stdout_file, $tblout_file, $err_file, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR    = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $do_local = opt_Get("--local", $opt_HHR) ? 1 : 0;

  validateFileExistsAndIsNonEmpty($model_file, $sub_name, $FH_HR); 
  validateFileExistsAndIsNonEmpty($seq_file,   $sub_name, $FH_HR);

  my $opts = " --cpu 0 --tblout $tblout_file";
  if(opt_Get("-v", $opt_HHR)) { $opts .= " --noali"; }

  # remove the tblout file or stdout files if they exist
  if(-e $stdout_file) { removeFileUsingSystemRm($stdout_file, $sub_name, $opt_HHR, $ofile_info_HHR); }
  if(-e $tblout_file) { removeFileUsingSystemRm($tblout_file, $sub_name, $opt_HHR, $ofile_info_HHR); }

  my $cmd = "$executable $opts $model_file $seq_file > $stdout_file";

  if($do_local) { 
    runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  }
  else { 
    my $job_name = "J" . removeDirPath($seq_file);
    my $nsecs  = opt_Get("--wait", $opt_HHR) * 60.;
    my $mem_gb = 8.; # hardcoded for nhmmscan
    submitJob($cmd, $job_name, $err_file, $mem_gb, $nsecs, $opt_HHR, $ofile_info_HHR);
  }
  return; 
}

#################################################################
# Subroutine : cmalignCheckStdOutput()
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
sub cmalignCheckStdOutput { 
  my $sub_name = "cmalignCheckStdOutput";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($stdout_file, $ret_mxsize_R, $FH_HR) = @_;
  if(defined $ret_mxsize_R) { 
    $$ret_mxsize_R = 0; # overwritten below if nec
  }

  if(! -e $stdout_file) { 
    DNAORG_FAIL("ERROR in $sub_name, cmalign stdout file $stdout_file does not exist", 1, $FH_HR);
  }
  if(! -s $stdout_file) { 
    DNAORG_FAIL("ERROR in $sub_name, cmalign stdout file $stdout_file exists but is empty", 1, $FH_HR);
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

#################################################################
# Subroutine : cmalignStoreOverflow()
# Incept:      EPN, Wed Feb 13 16:04:53 2019
#
# Purpose:     Store information on a sequence that has caused
#              a DP matrix memory overflow. 
#              
# Arguments: 
#  $seq_file:           the sequence file with the single sequence that failed in it
#  $mxsize:             matrix size to add to @{$overflow_mxsize_AR}
#  $overflow_seq_AR:    ref to array of sequences that failed due to matrix overflows, to add to
#  $overflow_mxsize_AR: ref to array of required matrix sizes for each sequence that failed due to matrix overflows, to add to
#  $FH_HR:              ref to file handle hash
# 
# Returns:     void
#
# Dies: if there's some problem opening the sequence file
#
################################################################# 
sub cmalignStoreOverflow { 
  my $sub_name = "cmalignStoreOverflow";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seq_file, $mxsize, $overflow_seq_AR, $overflow_mxsize_AR, $FH_HR) = @_;

  my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $seq_file }); # the sequence file object

  my $r2_seqstring = $sqfile->fetch_consecutive_seqs(1, "", -1); 
  if($r2_seqstring =~ /^\>(\S+)/) { 
    my $overflow_seq = $1;
    push(@{$overflow_seq_AR},    $overflow_seq);
    push(@{$overflow_mxsize_AR}, $mxsize);
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name failed to parse sequence name from matrix overflow sequence:\n$r2_seqstring", 1, $FH_HR);
  }
  $sqfile = undef;

  return;
}

#################################################################
# Subroutine: featureTypeIsCds()
# Incept:     EPN, Tue Feb  6 10:47:58 2018
#
# Purpose:    Is feature $ftr_idx a CDS?
#
# Arguments:
#   $ftr_info_AHR:  REF to array of hashes with info on the features, PRE-FILLED
#   $ftr_idx:       feature index we are interested in
#             
# Returns:    1 or 0
#
# Dies:       never; does not validate anything.
#
#################################################################
sub featureTypeIsCds { 
  my $sub_name  = "featureTypeIsCds";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_AHR, $ftr_idx) = (@_);

  return($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS") ? 1 : 0;
}

#################################################################
# Subroutine: featureTypeIsMaturePeptide()
# Incept:     EPN, Sun Mar 10 06:57:01 2019
#
# Purpose:    Is feature $ftr_idx a mat_peptide?
#
# Arguments:
#   $ftr_info_AHR:  REF to array of hashes with info on the features, PRE-FILLED
#   $ftr_idx:       feature index we are interested in
#             
# Returns:    1 or 0
#
# Dies:       never; does not validate anything.
#
#################################################################
sub featureTypeIsMaturePeptide { 
  my $sub_name  = "featureTypeIsMaturePeptide";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_AHR, $ftr_idx) = (@_);

  return($ftr_info_AHR->[$ftr_idx]{"type"} eq "mat_peptide") ? 1 : 0;
}


#################################################################
# Subroutine: featureNumSegments()
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
sub featureNumSegments { 
  my $sub_name  = "featureNumSegments";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_AHR, $ftr_idx) = (@_);

  return ($ftr_info_AHR->[$ftr_idx]{"3p_seg_idx"} - $ftr_info_AHR->[$ftr_idx]{"5p_seg_idx"} + 1);
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
    DNAORG_FAIL("ERROR in $sub_name, unable to parse input error string: $errstr", 1, $FH_HR);
  }

  if($error_desc eq "") { 
    $error_desc = "-";
  }
  return $seqname . "\t" . $error_name . "\t" . $feature_name . "\t" . $error_desc;
}

#################################################################
# Subroutine:  blastxDbSeqnameToFtrIdx()
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
sub blastxDbSeqNameToFtrIdx { 
  my $sub_name = "blastxDbSeqNameToFtrIdx";
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
            DNAORG_FAIL("ERROR in $sub_name, found blastx db sequence with coords that match two features, ftr_idx: $ftr_idx and $ret_ftr_idx", 1, $FH_HR);
          }                  
          $ret_ftr_idx = $ftr_idx;
        }
      }
    }
    if(! defined $ret_ftr_idx) { # did not find match
      DNAORG_FAIL("ERROR in $sub_name, did not find matching feature for blastx db sequence $blastx_seqname", 1, $FH_HR);
    }
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unable to parse blastx db sequence name $blastx_seqname", 1, $FH_HR); 
  }

  return $ret_ftr_idx;
}

#################################################################
# Subroutine:  validateBlastDbExists()
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
sub validateBlastDbExists {
  my $sub_name = "validateBlastDbExists";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastdb_name, $FH_HR) = @_;

  foreach my $sfx (".phr", ".pin", ".psq") { 
    if(! -s $blastdb_name . $sfx) { 
      DNAORG_FAIL("ERROR in $sub_name, required blast DB file " . $blastdb_name . $sfx . " does not exist or is empty", 1, $FH_HR); 
    }
  }

  return;
}


#################################################################
# Subroutine:  featureTypeIsCdsOrMaturePeptide()
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
sub featureTypeIsCdsOrMaturePeptide { 
  my $sub_name = "featureTypeIsCdsOrMaturePeptide";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_idx) = @_;

  return((featureTypeIsCds($ftr_info_AHR, $ftr_idx)) || 
         (featureTypeIsMaturePeptide($ftr_info_AHR, $ftr_idx)));
}


#################################################################
# Subroutine:  featureIsDuplicate()
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
sub featureIsDuplicate { 
  my $sub_name = "featureIsDuplicate";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $ftr_idx) = @_;

  return(($ftr_info_HAR->{"source_idx"}[$ftr_idx] != $ftr_idx) ? 1 : 0);
}

#################################################################
# Subroutine:  featureGet5pMostPosition()
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
sub featureGet5pMostPosition { 
  my $sub_name = "featureGet5pMostPosition";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords, $FH_HR) = @_;
  
  if($coords =~ /^(\d+)\.\.\d+/) { 
    return $1;
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unable to parse ftr_info_HA coords string " . $coords, 1, $FH_HR); 
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine:  featureGet3pMostPosition()
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
sub featureGet3pMostPosition { 
  my $sub_name = "featureGet3pMostPosition";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords, $FH_HR) = @_;
  
  if($coords =~ /\d+\.\.(\d+)\:[\+\-]$/) { 
    return $1;
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unable to parse ftr_info_HA coords string " . $coords, 1, $FH_HR); 
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine:  featureInfoSummarizeSegment()
# Incept:      EPN, Fri Mar  1 12:36:36 2019
#
# Purpose:    Return a string indicating what model this is
#             for features that are covered by multiple model spans.
#
# Arguments: 
#  $ftr_info_AHR: ref to feature info array of hashes, PRE-FILLED
#  $seg_info_AHR: ref to segment info array of hashes, PRE-FILLED
#  $seg_idx:      model index
#
# Returns:    "" if this is the only model for this feature
#             string like ", model 1 of 2", if not
# 
# Dies:       never; does not validate anything.
#
################################################################# 
sub featureInfoSummarizeSegment { 
  my $sub_name = "featureInfoSummarizeSegment";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $seg_info_AHR, $seg_idx) = @_;

  my $ftr_idx = $seg_info_AHR->[$seg_idx]{"map_ftr"};
  my $nmdl = ($ftr_info_AHR->[$ftr_idx]{"3p_seg_idx"} - $ftr_info_AHR->[$ftr_idx]{"5p_seg_idx"}) + 1;
  if($nmdl > 1) { 
    return sprintf(", segment %d of %d", ($seg_idx - $ftr_info_AHR->[$ftr_idx]{"5p_seg_idx"}) + 1, $nmdl);
  }

  return ""; # return "" if $nmdl == 1;
}

#################################################################
# Subroutine: featureInfoImputeCoords
# Incept:     EPN, Wed Mar 13 13:15:33 2019
# 
# Purpose:    Fill "coords" values in %{$ftr_info_HAR}
# 
# Arguments:
#   $ftr_info_AHR:   REF to feature information, added to here
#   $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub featureInfoImputeCoords { 
  my $sub_name = "featureInfoImputeCoords";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $opt_HHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "location"
  my @keys_A = ("type", "location");
  my $nftr = arrayOfHashesValidate($ftr_info_AHR, \@keys_A, $FH_HR);

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"coords"} = featureInfoCoordsFromLocation($ftr_info_AHR->[$ftr_idx]{"location"}, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: featureInfoImputeSourceIdx
# Incept:     EPN, Wed Mar 13 13:20:01 2019
# 
# Purpose:    Fill "source_idx" values in @{$ftr_info_AHR}
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $opt_HHR:       REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub featureInfoImputeSourceIdx { 
  my $sub_name = "featureInfoImputeSourceIdx";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $opt_HHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = arrayOfHashesValidate($ftr_info_AHR, \@keys_A, $FH_HR);

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
            DNAORG_FAIL(sprintf("ERROR in $sub_name, unable to determine source (two choices) for duplicate feature of type %s and coords %s\n", 
                                $ftr_info_AHR->[$ftr_idx]{"type"}, $ftr_info_AHR->[$ftr_idx]{"coords"}), 1, $FH_HR);
          }
          $ftr_info_AHR->[$ftr_idx]{"source_idx"} = $ftr_idx2;
        }
      }
    }
  }

  return;
}

#################################################################
# Subroutine: featureInfoImputeParentIdx
# Incept:     EPN, Wed Mar 13 13:33:33 2019
# 
# Purpose:    Fill "parent_idx" values in @{$ftr_info_AHR}
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $opt_HHR:       REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub featureInfoImputeParentIdx {
  my $sub_name = "featureInfoImputeParentIdx";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $opt_HHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = arrayOfHashesValidate($ftr_info_AHR, \@keys_A, $FH_HR);

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
      $ftr_5p_pos = featureGet5pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
      $ftr_3p_pos = featureGet3pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
      $ftr_strand = featureSummaryStrand($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
      for($ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
        $ftr_5p_pos2 = featureGet5pMostPosition($ftr_info_AHR->[$ftr_idx2]{"coords"}, $FH_HR);
        $ftr_3p_pos2 = featureGet3pMostPosition($ftr_info_AHR->[$ftr_idx2]{"coords"}, $FH_HR);
        $ftr_strand2 = featureSummaryStrand($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
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
              DNAORG_FAIL(sprintf("ERROR in $sub_name, unable to determine parent of mature peptide with coords %s (multiple CDS cover it with coords %s and %s)\n", 
                                  $ftr_info_AHR->[$ftr_idx]{"coords"},
                                  $ftr_info_AHR->[($ftr_info_AHR->[$ftr_idx]{"parent_idx"})]{"coords"}, 
                                  $ftr_info_AHR->[$ftr_idx2]{"coords"}), 1, $FH_HR);
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
# Subroutine: featureInfoImpute3paFtrIdx
# Incept:     EPN, Wed Mar 13 13:39:34 2019
# 
# Purpose:    Fill "3pa_ftr_idx" values in @{$ftr_info_AHR}
# 
# Arguments:
#   $ftr_info_AHR:  REF to feature information, added to here
#   $opt_HHR:       REF to 2D hash of option values, see top of epn-options.pm for description
#   $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
# 
# Dies:       if $ftr_info_AHR is invalid upon entry
#
#################################################################
sub featureInfoImpute3paFtrIdx {
  my $sub_name = "featureInfoImpute3paFtrIdx";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($ftr_info_AHR, $opt_HHR, $FH_HR) = @_;
  
  # ftr_info_AHR should already have array data for keys "type", "coords"
  my @keys_A = ("type", "coords");
  my $nftr = arrayOfHashesValidate($ftr_info_AHR, \@keys_A, $FH_HR);

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
      $ftr_3p_pos = featureGet3pMostPosition($ftr_info_AHR, $ftr_idx, $FH_HR);
      $ftr_strand = featureSummaryStrand($ftr_info_AHR, $ftr_idx, $FH_HR);
      for($ftr_idx2 = 0; $ftr_idx2 < $nftr; $ftr_idx2++) { 
        $ftr_5p_pos2 = featureGet5pMostPosition($ftr_info_AHR, $ftr_idx2, $FH_HR);
        $ftr_strand2 = getFeatureStrand($ftr_info_AHR, $ftr_idx2, $FH_HR);
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
              DNAORG_FAIL(sprintf("ERROR in $sub_name, unable to determine 3' mature peptide of mature peptide with coords (multiple mature peptides satisfy criteria)\n", 
                                  $ftr_info_AHR->[$ftr_idx]{"coords"}), 1, $FH_HR);
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
# Subroutine: segmentInfoPopulate()
# Incept:     EPN, Wed Mar 13 13:55:56 2019
#
# Synopsis: Fill @{$seg_info_AHR} and add to @{$ftr_info_AHR}
#           based on existing information in @{$ftr_info_AHR}.
#
#           The following values are added to %{$seg_info_HAR}:
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
#                "5p_seg_idx":   index (in arrays of %seg_info_HA) of 5'-most segment for this feature
#                "3p_seg_idx":   index (in arrays of %seg_info_HA) of 3'-most segment for this feature
# Arguments:
#  $seg_info_AHR:      ref to array of hashes with information on the segments, FILLED HERE
#  $ftr_info_AHR:      ref to array of hashes with information on the features, ADDED TO HERE
#  $opt_HHR:           REF to 2D hash of option values
#  $FH_HR:             REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if @{$ftr_info_AHR} is not valid upon entry
#
#################################################################
sub segmentInfoPopulate {
  my $sub_name = "segmentInfoPopulate";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seg_info_AHR, $ftr_info_AHR, $opt_HHR, $FH_HR) = @_;

  # ftr_info_AHR should already have array data for keys "type", "coords", "source_idx"
  my @keys_A = ("type", "coords", "source_idx");
  my $nftr = arrayOfHashesValidate($ftr_info_AHR, \@keys_A, $FH_HR);

  # initialize new %{$ftr_info_AHR} values
  my ($ftr_idx, $ftr_idx2, $seg_idx, $seg_idx2); # feature and segment indices
  my ($seg_start, $seg_stop, $seg_strand); # start, stop and strand for a segment
  my $nseg = 0; 
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $ftr_info_AHR->[$ftr_idx]{"5p_seg_idx"} = -1; # remains -1 if $ftr_idxs_dup
    $ftr_info_AHR->[$ftr_idx]{"3p_seg_idx"} = -2; # remains -2 if $ftr_idxs_dup
    my $ftr_type   = $ftr_info_AHR->[$ftr_idx]{"type"};
    my $ftr_is_dup = ($ftr_info_AHR->[$ftr_idx]{"source_idx"} == $ftr_idx) ? 0 : 1;

    if(! $ftr_is_dup) { 
      # determine start and stop positions of all segments
      my @seg_start_A  = (); # array of starts, one per segment
      my @seg_stop_A   = (); # array of stops, one per segment
      my @seg_strand_A = (); # array of strands ("+", "-"), one per segment
      featureStartStopStrandArrays($ftr_info_AHR->[$ftr_idx]{"coords"}, \@seg_start_A, \@seg_stop_A, \@seg_strand_A, $FH_HR);
      my $cur_nseg = scalar(@seg_start_A);
      for(my $s = 0; $s < $cur_nseg; $s++) { 
        $seg_info_AHR->[$nseg]{"start"}   = $seg_start_A[$s];
        $seg_info_AHR->[$nseg]{"stop"}    = $seg_stop_A[$s];
        $seg_info_AHR->[$nseg]{"strand"}  = $seg_strand_A[$s];
        $seg_info_AHR->[$nseg]{"map_ftr"} = $ftr_idx;
        if($s == 0) { 
          $seg_info_AHR->[$nseg]{"is_5p"} = 1;
          $ftr_info_AHR->[$ftr_idx]{"5p_seg_idx"} = $nseg; 
        }
        else { 
          $seg_info_AHR->[$nseg]{"is_5p"} = 0;
        }
        if($s == ($cur_nseg-1)) { 
          $seg_info_AHR->[$nseg]{"is_3p"} = 1;
          $ftr_info_AHR->[$ftr_idx]{"3p_seg_idx"} = $nseg;
        }
        else { 
          $seg_info_AHR->[$nseg]{"is_3p"} = 0;
        }
        $nseg++;
      }
    }
  }

  return;
}

#################################################################
# Subroutine: featureStartStopStrandArrays()
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
# Dies: if unable to parse $coords_str
#
#################################################################
sub featureStartStopStrandArrays {
  my $sub_name = "featureStartStopStrandArrays";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($coords, $start_AR, $stop_AR, $strand_AR, $FH_HR) = @_;

  my @start_A  = ();
  my @stop_A   = ();
  my @strand_A = ();
  my ($start, $stop, $strand, $seg_idx);
  my @coords_A  = split(",", $coords);
  my $nseg = scalar(@coords_A);
  for($seg_idx = 0; $seg_idx < $nseg; $seg_idx++) { 
    if($coords_A[$seg_idx] =~ /^(\d+)\.\.(\d+)\:([\+\-])$/) { 
      ($start, $stop, $strand) = ($1, $2, $3);
      push(@start_A,  $start);
      push(@stop_A,   $stop);
      push(@strand_A, $strand); 
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse coords token $coords_A[$seg_idx]", 1, $FH_HR); 
    }
  }

  if(defined $start_AR)  { @{$start_AR}   = @start_A;  }
  if(defined $stop_AR)   { @{$stop_AR}    = @stop_A;   }
  if(defined $strand_AR) { @{$strand_AR}  = @strand_A;  }

  return;
}

#################################################################
# Subroutine: arrayOfHashesValidate()
# Incept:     EPN, Wed Mar 13 13:24:38 2019
#
# Purpose:    Validate an array of hashes, by making sure it
#             includes a key/value for all keys in @{$keys_AR}.
# Arguments:
#   $AHR:      REF to array of hashes to validate
#   $keys_AR:  REF to array of keys that may be excluded from the hash
#   $FH_HR:    REF to hash of file handles, including "log" and "cmd"
# 
# Returns: scalar(@{$AHR});
#
# Dies:    - if one of the keys in @{$keys_AR} does not exist in all hashes
#            of the array
#
#################################################################
sub arrayOfHashesValidate {
  my $sub_name = "arrayOfHashesValidate()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($AHR, $keys_AR, $FH_HR) = (@_);

  my $n = scalar(@{$AHR});

  for(my $i = 0; $i < $n; $i++) { 
    foreach my $key (@{$keys_AR}) { 
      if(! exists $AHR->[$i]{$key}) { 
        DNAORG_FAIL("ERROR in $sub_name, array element $i does not have key $key\n"); 
      }
    }
  }

  return $n;
}


#################################################################
# Subroutine: featureInfoCoordsFromLocation
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
sub featureInfoCoordsFromLocation { 
  my $sub_name = "featureInfoCoordsFromLocation";
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
    $ret_val = featureInfoCoordsFromLocation($location_to_join, $FH_HR);
  }
  elsif($location =~ /^complement\((.+)\)$/) { 
    my $location_to_complement = $1;
    my $coords_to_complement = featureInfoCoordsFromLocation($location_to_complement, $FH_HR);
    $ret_val = featureInfoCoordsComplement($coords_to_complement, $FH_HR);
  }
  elsif($location =~ /\,/) { 
    # not wrapped in join() or complement(), but multiple segments
    foreach my $location_el (split(",", $location)) { 
      if($ret_val ne "") { $ret_val .= ","; }
      $ret_val .= featureInfoCoordsFromLocation($location_el, $FH_HR);
    }
  }
  elsif($location =~ /\<?(\d+)\.\.(\d+)\>?/) { 
    $ret_val = $1 . ".." . $2 . ":" . "+"; # a recursive call due to the complement() may complement this
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unable to parse location token $location", 1, $FH_HR);
  }

  return $ret_val;
}

#################################################################
# Subroutine: featureInfoCoordsComplement
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
sub featureInfoCoordsComplement { 
  my $sub_name = "featureInfoCoordsComplement";
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
    if($el_A[$i] =~ /(\d+)\.\.(\d+)\:\+/) { 
      if($ret_val ne "") { $ret_val .= ","; }
      $ret_val .= $2 . ".." . $1. ":-";
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse coords token $coords", 1, $FH_HR);
    }
  }

  printf("\tin $sub_name, coords: $coords ret_val: $ret_val\n");

  return $ret_val;
}

#################################################################
# Subroutine: featureSummaryStrand
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
sub featureSummaryStrand { 
  my $sub_name = "featureSummaryStrand";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($coords, $FH_HR) = @_;

  # Examples we can parse: 
  # $coords                  return value
  # -----------------------  -----------------
  # 1-200:+                  200-1:-
  # 1-200:+,300-400:+        400-300:-,200-1:-

  my @strand_A = ();
  featureStartStopStrandArrays($coords, undef, undef, \@strand_A, $FH_HR);

  my $npos = 0;
  my $nneg = 0;
  foreach my $strand (@strand_A) { 
    if   ($strand eq "+") { $npos++; }
    elsif($strand eq "-") { $nneg++; }
    else { DNAORG_FAIL("ERROR in $sub_name, unable to determine strands in coords $coords", 1, $FH_HR); }
  }

  if(($npos >  0) && ($nneg == 0)) { return "+"; }
  if(($npos == 0) && ($nneg >  0)) { return "-"; }
  if(($npos == 0) && ($nneg == 0)) { 
    DNAORG_FAIL("ERROR in $sub_name, unable to determine strands in coords $coords", 1, $FH_HR); 
  }

  return; # NEVER REACHED
}
  
#################################################################
# Subroutine : sqstringAddNewlines()
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
sub sqstringAddNewlines { 
  my $nargs_expected = 2;
  my $sub_name = "sqstringAddNewlines";
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
    

###########################################################################
# the next line is critical, a perl module must return a true value
return 1;
###########################################################################
