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
#   determineFeatureTypes()
#   getNumFeaturesAnnotatedByModels()
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
#   setRequiredErrorInfoHash()
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
#   determineFeatureTypes()
#   getNumFeaturesAnnotatedByModels()
#   getReferenceFeatureInfo()
#   fetchReferenceFeatureSequences()
#   featureHasChildren()
#   featureHasParent()
#

#################################################################
# Subroutine: determineFeatureTypes()
# Incept:     EPN, Thu Feb 11 14:50:53 2016
#
# Purpose:    Determine the type of each feature, either
#             'cds':       CDS, possibly comprised of mature peptides
#             'mp':        mature peptide
#             'xfeat'      extra feature annotated by its own model (e.g. gene)
#             'dfeat'      extra feature annotated as a duplicate of another feature (e.g. gene)
# 
#              And populate the %{$ftr_info_HAR} with the following
#              information (1D keys):
#              "type"        values are one of: "cds", "mp", "xfeat" or "dfeat"
#              "type_idx"    values are indices of each feature for its type
#                            so a '2' means it's the 2nd of its type.
#              "annot_type": annotation type, either "model" or "duplicate"
#                            "model": the homology model's prediction annotates 
#                            these features
#                            "duplicate": the annotation of this feature is copied from 
#                            another. Example is a 'gene' that is a copy of a 'CDS'.
#                            The feature that is the source to copy is stored in the
#                            "source_idx" key.
#               "primary_children_ftr_str": for features with children (e.g. cds comprised of
#                            mature peptides) a string of feature indices that are the primary 
#                            children of this feature (when the primary children are concatenated 
#                            they make up the CDS that encodes all the primary children). Each 
#                            index is separated by a space. Empty for features with annot_type eq 
#                            "model".
#               "primary_children_ftr_num": number of primary children (0 for features with
#                            annot_type eq "model").
#               "all_children_ftr_str": for features with children (e.g. cds comprised of
#                            mature peptides) a string of feature indices that are encoded by 
#                            this CDS includes all primary children plus any secondarily processed
#                            mature peptides. Each index is separated by a space. Empty 
#                            for features with annot_type eq "model."
#               "all_children_ftr_num": number of all children (0 for features with
#                            annot_type eq "model").
#               "parent_ftr": index of parent feature for 'mp' types, this is the feature
#                            index of the feature of which the current feature is a child.
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

  # order: 
  # Mature peptides
  # CDS
  # xfeat
  # dfeat
  for($c = 0; $c < $nmp; $c++) { 
    determineFeatureTypesHelper($ftr_info_HAR, "mp", undef, undef, $FH_HR);
  }
  for($c = 0; $c < $ncds; $c++) { 
    # determine type of CDS
    if((defined $cds2pmatpept_AAR) && (defined $cds2pmatpept_AAR->[$c])) {
      # CDS comprised of mature peptides
      determineFeatureTypesHelper($ftr_info_HAR, "cds", $cds2pmatpept_AAR->[$c], $cds2amatpept_AAR->[$c], $FH_HR);
    }
    else { 
      # CDS NOT comprised of mature peptides
      determineFeatureTypesHelper($ftr_info_HAR, "cds", undef, undef, $FH_HR);
    }
  }
  for(my $xfeat_idx = 0; $xfeat_idx < $nxfeat; $xfeat_idx++) { 
    determineFeatureTypesHelper($ftr_info_HAR, "xfeat", undef, undef, $FH_HR);
  }
  for(my $dfeat_idx = 0; $dfeat_idx < $ndfeat; $dfeat_idx++) { 
    determineFeatureTypesHelper($ftr_info_HAR, "dfeat", undef, undef, $FH_HR);
  }
  return;
}

#################################################################
# Subroutine: determineFeatureTypesHelper()
# Incept:     EPN, Fri Feb 15 11:23:58 2019
#
# Purpose:    Initialize element $c of the feature info hash of arrays
#             of type $type.
#
# Arguments: 
#   $ftr_info_HAR:     ref to hash of arrays with feature info.
#                      must have at least the "annot_types" array, else we die 
#   $type:             type of the feature: "cds-notmp", "xfeat", "dfeat", "mp", or "cds-mp"
#   $pchildren_AR:     array of child indices for this feature OR undefined
#                      if no primary children (PRE-FILLED)
#   $achildren_AR:     array of child indices for this feature OR undefined
#                      if no children (PRE-FILLED)
#   $FH_HR:            REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if $ftr_info_HAR is not valid upon entering.
#################################################################
sub determineFeatureTypesHelper { 
  my $sub_name  = "determineFeatureTypesHelper()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_HAR, $type, $pchildren_AR, $achildren_AR, $FH_HR) = (@_);

  my @reqd_ftr_info_A = ();
  # we can't use validateAndGetSizeOfInfoHashOfArrays() because we probably will not have the 
  # same number of elements for all values, e.g. "type" may be > 0, but "type" may be 0
  my $nftr_type = scalar(@{$ftr_info_HAR->{"type"}});

  # determine type index
  my $type_idx = 1;
  my $ftr_idx;
  for($ftr_idx = 0; $ftr_idx < $nftr_type; $ftr_idx++) { 
    if($ftr_info_HAR->{"type"}[$ftr_idx] eq $type) { 
      $type_idx++;
    }
  }
  $ftr_idx = $nftr_type; 

  # initialize
  $ftr_info_HAR->{"type"}[$ftr_idx]                     = $type;
  $ftr_info_HAR->{"type_idx"}[$ftr_idx]                 = $type_idx;
  $ftr_info_HAR->{"primary_children_ftr_str"}[$ftr_idx] = "";  # possibly created later 
  $ftr_info_HAR->{"primary_children_ftr_num"}[$ftr_idx] = 0;   # possibly created later
  $ftr_info_HAR->{"all_children_ftr_str"}[$ftr_idx]     = "";  # possibly created later
  $ftr_info_HAR->{"all_children_ftr_num"}[$ftr_idx]     = 0;   # possibly created later
  $ftr_info_HAR->{"parent_ftr"}[$ftr_idx]               = -1;  # changed to a valid value for certain types (e.g. 'mp')

  # set annot_type, children strings and validate $type at the same time:
  if(($type eq "cds")   ||
     ($type eq "xfeat") ||
     ($type eq "mp")) { 
    $ftr_info_HAR->{"annot_type"}[$ftr_idx]  = "model"; # this feature is annotated by homology models
  }
  elsif($type eq "dfeat") { 
    $ftr_info_HAR->{"annot_type"}[$ftr_idx]  = "duplicate"; # this feature is modeled as a duplicate of another feature
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, invalid type: $type.", 1, $FH_HR);
  }

  my $child_idx;
  my $z;
  if(defined $achildren_AR) { 
    # step through @{$achildren_AR} and create the all_children_ftr_str for this CDS
    my $na = scalar(@{$achildren_AR});
    for($z = 0; $z < $na; $z++) { 
      if($ftr_info_HAR->{"all_children_ftr_str"}[$ftr_idx] ne "") { 
        $ftr_info_HAR->{"all_children_ftr_str"}[$ftr_idx] .= " ";
      }
      $child_idx = $achildren_AR->[$z];
      $ftr_info_HAR->{"all_children_ftr_str"}[$ftr_idx] .= $child_idx;
      $ftr_info_HAR->{"all_children_ftr_num"}[$ftr_idx]++;
      # set the parent ftr index for this child
      $ftr_info_HAR->{"parent_ftr"}[$child_idx] = $ftr_idx;
    }          
  }
  if(defined $pchildren_AR) { 
    # step through @{$pchildren_AR} and create the primary_children_ftr_str for this CDS
    my $np = scalar(@{$pchildren_AR});
    for($z = 0; $z < $np; $z++) { 
      if($ftr_info_HAR->{"primary_children_ftr_str"}[$ftr_idx] ne "") { 
        $ftr_info_HAR->{"primary_children_ftr_str"}[$ftr_idx] .= " ";
      }
      $child_idx = $pchildren_AR->[$z];
      $ftr_info_HAR->{"primary_children_ftr_str"}[$ftr_idx] .= $child_idx;
      $ftr_info_HAR->{"primary_children_ftr_num"}[$ftr_idx]++;
      # make sure this child was in the all_children array by checking "parent_ftr"
      if($ftr_info_HAR->{"parent_ftr"}[$child_idx] == -1) { 
        DNAORG_FAIL("ERROR in $sub_name, child $child_idx is primary children list but not in all children list", 1, $FH_HR);
      }
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
    if($ref_mp_strand_str !~ /^[+-]+$/) { 
      DNAORG_FAIL("ERROR in $sub_name, at least one MP is not + or - strand (must be single nt MP or mixed strand MP), ref_mp_strand_str: $ref_mp_strand_str", 1, $FH_HR);
    }
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
  # CDS strand must be all + or -, we can't deal with unknown strands or mixed strands because we need to 
  # identify starts and stops in CDS 
  if($ref_cds_strand_str !~ /^[+-]+$/) { 
    DNAORG_FAIL("ERROR in $sub_name, at least one CDS is not + or - strand (must be single nt CDS or mixed strand CDS), ref_cds_strand_str: $ref_cds_strand_str", 1, $FH_HR);
  }
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
#                "is_final":      '1' if this model is the final model (e.g. final segment/exon) for the feature it models ("map_ftr")
#                "is_first":      '1' if this model is the first model (e.g. final segment/exon) for the feature it models ("map_ftr")
#                "map_segment":   the segment index this model models (1.."map_nsegment" value)
#                "map_nsegment":  the number of segments the feature this model models has 
#                "out_idx":       output value: feature index and segment index this model (e.g. "4.2")
#
#           The following values are filled in the %{$ftr_info_HAR}:
#                "final_mdl":   index (in arrays of %mdl_info_HA) of final model for this feature
#                "first_mdl":   index (in arrays of %mdl_info_HA) of first model for this feature
#                "nmodels":     number of models for this feature (e.g. number of segments) for this feature, 
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
#  $stk_file:          name of output file we will write the stockholm single sequence 
#                      'alignment' to
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

  my ($execs_HR, $sqfile, $ref_seq_accn, $ref_totlen, $out_root, $build_root, $mdl_info_HAR, $ftr_info_HAR, $stk_file, $opt_HHR, $FH_HR) = @_;

  # contract check

  # ftr_info_HAR should have array data for keys "ref_coords", "ref_strand"
  my @reqd_ftr_info_A  = ("ref_coords", "ref_strand", "annot_type");
  my $nftr             = validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, \@reqd_ftr_info_A, $FH_HR);
  my $nftr_with_models = getNumFeaturesAnnotatedByModels($ftr_info_HAR, $FH_HR);

  # determine if we read the sequences from a provided fasta file
  # if this is true, we won't fetch the reference sequences
  my $do_infasta   = (opt_Exists("--infasta",   $opt_HHR) && opt_Get("--infasta",   $opt_HHR)) ? 1 : 0;

  my $do_keep     = opt_Get("--keep", $opt_HHR); # should we leave intermediates files on disk, instead of removing them?
  my $do_circular = 0; # placeholder for possibly circular genomes in the future, always FALSE for now
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
  my @ref_nsegments_A   = (); # [0..$c..$ref_nmft-1]: number of segments in CDS or mat_peptide $c+1
  my $ref_tot_nsegments = 0;  # total number of segments in all CDS or mat_peptides
  my @indi_ref_name_A   = (); # [0..$nmdl-1]: name of individual stockholm alignments and models
  my @indi_cksum_stk_A  = (); # [0..$nmdl-1]: checksum's of each named individual stockholm alignment

  my @files2rm_A = ();  # array of file names to remove at end of this function (remains empty if $do_keep)

  my $cksum  = -1;
  my $mdllen = undef;
  my $model_name = removeScriptNameFromString(removeDirPath($build_root)); 

  if(! $do_infasta) { 
    # fetch the sequence and convert it to stockholm
    my $fa_file = $out_root . ".ref.fa";
    if($do_circular) { 
      # placeholder for circulare genomes, in that case, we build a model of *duplicated* full sequence
      DNAORG_FAIL("ERROR in $sub_name, do_circular is TRUE", 1, $FH_HR);
      my $fa_string = $sqfile->fetch_seq_to_fasta_string($ref_seq_accn, -1);
      my @lines_A = split("\n", $fa_string);
      if(scalar(@lines_A) != 2) { 
        DNAORG_FAIL("ERROR in $sub_name, error fetching $ref_seq_accn into a fasta string", 1, $FH_HR);
      }
      open(FA, ">", $fa_file) || fileOpenFailure($fa_file, $sub_name, $!, "writing", $FH_HR);
      print FA $lines_A[0] . "\n"; # fasta header line
      print FA $lines_A[1] . "\n"; # fasta sequence line
      print FA $lines_A[1] . "\n"; # fasta sequence line again (to duplicate)
      close(FA);
    }
    else { 
      # ! $do_circular, normal case
      my @fetch_A = ($ref_seq_accn);
      $sqfile->fetch_seqs_given_names(\@fetch_A, -1, $fa_file);
    }
    my $tmp_stk_file = $stk_file . ".tmp";
    my $cmd = "$esl_reformat --informat afa stockholm $fa_file > $tmp_stk_file";
    runCommand($cmd, 0, 0, $FH_HR);
    if(! $do_keep) { push(@files2rm_A, $tmp_stk_file); }
    
    # annotate the stockholm file with a blank SS and with a name
    my (undef, $cksum) = addNameAndBlankSsToStockholmAlignment($model_name, 1, $tmp_stk_file, $stk_file, $FH_HR); # 1: add blank SS_cons line
  }

  for(my $i = 0; $i < $nftr; $i++) { 
    my $type_fname    = $ftr_info_HAR->{"type_fname"}[$i];
    my $do_model      = ($ftr_info_HAR->{"annot_type"}[$i] eq "model") ? 1 : 0;
    my $ftr_type      = $ftr_info_HAR->{"type"}[$i];
    my $ftr_type_idx  = $ftr_info_HAR->{"type_idx"}[$i];
    $ftr_info_HAR->{"first_mdl"}[$i]  = -1; # remains -1 if $do_model is FALSE
    $ftr_info_HAR->{"final_mdl"}[$i]  = -1; # remains -1 if $do_model is FALSE
    $ftr_info_HAR->{"nmodels"}[$i]    = 0;  # remains 0 if $do_model is FALSE
    $ftr_info_HAR->{"source_idx"}[$i] = -1; # maybe changed later in determineSourcesOfDuplicateFeatures()

    if($do_model) { 
      # determine start and stop positions of all exons/segments
      my @starts_A = ();
      my @stops_A  = ();
      my $nsegments   = 0;
      startsStopsStrandsFromCoordsLength($ftr_info_HAR->{"ref_coords"}[$i], $ref_totlen, $do_circular, \@starts_A, \@stops_A, undef, \$nsegments, $FH_HR);
      $ftr_info_HAR->{"nmodels"}[$i] = $nsegments;
      my $strand = $ftr_info_HAR->{"ref_strand"}[$i];

      # if we're on the negative strand, reverse the arrays, they'll be in the incorrect order
      if($strand eq "-") { 
        @starts_A = reverse @starts_A;
        @stops_A  = reverse @stops_A;
      }

      for(my $e = 0; $e < $nsegments; $e++) { 
        if($nsegments > 1) { 
          $cur_out_root        = $out_root .       ".ref." . $type_fname . "." . $ftr_type_idx . ".segment." . ($e+1);
          $cur_out_name_root   = $out_dir_tail .   ".ref." . $type_fname . "." . $ftr_type_idx . ".segment." . ($e+1);
          $cur_build_root      = $build_root .     ".ref." . $type_fname . "." . $ftr_type_idx . ".segment." . ($e+1);
          $cur_build_name_root = $build_dir_tail . ".ref." . $type_fname . "." . $ftr_type_idx . ".segment." . ($e+1);
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

        $mdllen = abs($stop-$start)+1;

        $mdl_info_HAR->{"cmname"}[$nmdl]     = $model_name; # will be the same for all models
        $mdl_info_HAR->{"checksum"}[$nmdl]   = $cksum;      # will be the same for all models
        $mdl_info_HAR->{"length"}[$nmdl]     = $mdllen;
        $mdl_info_HAR->{"ref_start"}[$nmdl]  = $start;
        $mdl_info_HAR->{"ref_stop"}[$nmdl]   = $stop;
        $mdl_info_HAR->{"ref_strand"}[$nmdl] = $strand;

        # store information on this model's name for output purposes
        $mdl_info_HAR->{"filename_root"}[$nmdl] = sprintf("$ftr_type.%s", 
                                                          ($nsegments == 1) ? sprintf("%d", $ftr_type_idx) : sprintf("%d.%d", $ftr_type_idx, ($e+1)));
        $mdl_info_HAR->{"out_tiny"}[$nmdl]  = sprintf("%s#%s", 
                                                      $type_fname,
                                                      ($nsegments == 1) ? sprintf("%d", $ftr_type_idx) : sprintf("%d.%d", $ftr_type_idx, ($e+1)));
        
        $mdl_info_HAR->{"map_ftr"}[$nmdl]       = $i;
        $mdl_info_HAR->{"is_first"}[$nmdl]      = ($e == 0)           ? 1 : 0;
        $mdl_info_HAR->{"is_final"}[$nmdl]      = ($e == ($nsegments-1)) ? 1 : 0;
        $mdl_info_HAR->{"map_segment"}[$nmdl]   = $e;
        $mdl_info_HAR->{"map_nsegment"}[$nmdl]  = $nsegments;
        $mdl_info_HAR->{"out_idx"}[$nmdl]       = sprintf("%d.%d", 
                                                          $mdl_info_HAR->{"map_ftr"}[$nmdl]+1, $mdl_info_HAR->{"map_segment"}[$nmdl]+1);

        if($e == 0) { 
          $ftr_info_HAR->{"first_mdl"}[$i] = $nmdl;
        }
        if($e == ($nsegments-1)) { 
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
    $ftr_info_HAR->{"out_tiny"}[$i]      = $tiny;
    $ftr_info_HAR->{"out_short"}[$i]     = $short;
    $ftr_info_HAR->{"filename_root"}[$i] = $ftr_type . "." . $ftr_type_idx;
  }

  # clean up
  foreach my $file2rm (@files2rm_A) { 
    runCommand("rm $file2rm", 0, 0, $FH_HR);
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
# Subroutine : featureHasChildren()
# Incept:      EPN, Sat Feb 16 06:43:26 2019
#
# Purpose:     Return '1' if a feature has >= 1 children features.
#              Return '0' if a feature has    0 children features.
#              
# Arguments: 
#   $ftr_info_HAR: REF to hash of arrays with information on the features, ADDED TO HERE
#   $ftr_idx:      index of feature we are interested in
#   $FH_HR:        REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: If any of the following are undefined or have 
#       fewer than $ftr_idx elements.
#       @{$ftr_info_HAR->{"primary_children_ftr_str"}}
#       @{$ftr_info_HAR->{"primary_children_ftr_num"}}
#       @{$ftr_info_HAR->{"all_children_ftr_str"}}
#       @{$ftr_info_HAR->{"all_children_ftr_num"}}
# 
################################################################# 
sub featureHasChildren { 
  my $nargs_expected = 3;
  my $sub_name = "featureHasChildren";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_HAR, $ftr_idx, $FH_HR) = @_;

  foreach my $key ("primary_children_ftr_str", 
                "primary_children_ftr_num",
                "all_children_ftr_str",
                "all_children_ftr_num") { 
    if(! defined $ftr_info_HAR->{$key}) { 
      DNAORG_FAIL("ERROR in $sub_name, ftr_info_HAR->{$key} is not defined", 1, $FH_HR);
    }
    if(scalar(@{$ftr_info_HAR->{$key}}) < $ftr_idx) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, ftr_idx: $ftr_idx, ftr_info_HAR->{$key} only has %d elements", scalar(@{$ftr_info_HAR->{$key}})), 1, $FH_HR);
    }
  }        

  if(featureHasPrimaryChildren($ftr_info_HAR, $ftr_idx, $FH_HR)) { 
    return 1; 
  }
  if(featureHasAllChildren($ftr_info_HAR, $ftr_idx, $FH_HR)) { 
    return 1; 
  }

  return 0;
}

#################################################################
# Subroutine : featureHasPrimaryChildren()
# Incept:      EPN, Sat Feb 16 06:51:52 2019
#
# Purpose:     Return '1' if a feature has >= 1 primary children features.
#              Return '0' if a feature has    0 primary children features.
#              
# Arguments: 
#   $ftr_info_HAR: REF to hash of arrays with information on the features, ADDED TO HERE
#   $ftr_idx:      index of feature we are interested in
#   $FH_HR:        REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: If any of the following are undefined or have 
#       fewer than $ftr_idx elements.
#       @{$ftr_info_HAR->{"primary_children_ftr_str"}}
#       @{$ftr_info_HAR->{"primary_children_ftr_num"}}
#
################################################################# 
sub featureHasPrimaryChildren { 
  my $nargs_expected = 3;
  my $sub_name = "featureHasPrimaryChildren";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_HAR, $ftr_idx, $FH_HR) = @_;

  foreach my $key ("primary_children_ftr_str", 
                   "primary_children_ftr_num") { 
    if(! defined $ftr_info_HAR->{$key}) { 
      DNAORG_FAIL("ERROR in $sub_name, ftr_info_HAR->{$key} is not defined", 1, $FH_HR);
    }
    if(scalar(@{$ftr_info_HAR->{$key}}) < $ftr_idx) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, ftr_idx: $ftr_idx, ftr_info_HAR->{$key} only has %d elements", scalar(@{$ftr_info_HAR->{$key}})), 1, $FH_HR);
    }
  }        

  if($ftr_info_HAR->{"primary_children_ftr_num"}[$ftr_idx] > 0) { 
    return 1;
  }

  return 0;
}

#################################################################
# Subroutine : featureHasAllChildren()
# Incept:      EPN, Sat Feb 16 06:51:52 2019
#
# Purpose:     Return '1' if a feature has >= 1 all children features.
#              Return '0' if a feature has    0 all children features.
#              
# Arguments: 
#   $ftr_info_HAR: REF to hash of arrays with information on the features, ADDED TO HERE
#   $ftr_idx:      index of feature we are interested in
#   $FH_HR:        REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: If any of the following are undefined or have 
#       fewer than $ftr_idx elements.
#       @{$ftr_info_HAR->{"all_children_ftr_str"}}
#       @{$ftr_info_HAR->{"all_children_ftr_num"}}
#
################################################################# 
sub featureHasAllChildren { 
  my $nargs_expected = 3;
  my $sub_name = "featureHasAllChildren";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_HAR, $ftr_idx, $FH_HR) = @_;

  foreach my $key ("all_children_ftr_str", 
                   "all_children_ftr_num") { 
    if(! defined $ftr_info_HAR->{$key}) { 
      DNAORG_FAIL("ERROR in $sub_name, ftr_info_HAR->{$key} is not defined", 1, $FH_HR);
    }
    if(scalar(@{$ftr_info_HAR->{$key}}) < $ftr_idx) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, ftr_idx: $ftr_idx, ftr_info_HAR->{$key} only has %d elements", scalar(@{$ftr_info_HAR->{$key}})), 1, $FH_HR);
    }
  }        

  if($ftr_info_HAR->{"all_children_ftr_num"}[$ftr_idx] > 0) { 
    return 1;
  }

  return 0;
}

#################################################################
# Subroutine : featureGetPrimaryChildren()
# Incept:      EPN, Sat Feb 16 07:00:13 2019
#
# Purpose:     Fill @{$AR} with the space delimited tokens of 
#              @{$ftr_info_HAR->{"primary_children_ftr_str"}[$ftr_idx]}
#              and return;
# 
# Arguments: 
#   $ftr_info_HAR:   REF to hash of arrays with information on the features, PRE-FILLED
#   $ftr_idx:        index of feature we're interested in
#   $AR:             REF to array to fill, FILLED HERE
#   $FH_HR:          REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: If there are no primary children for feature index $ftr_idx,
#       If number of primary children in $ftr_info_HAR->{"primary_children_ftr_str"}
#       disagrees with $ftr_info_HAR->{"primary_children_ftr_num"}.
#
################################################################# 
sub featureGetPrimaryChildren { 
  my $nargs_expected = 4;
  my $sub_name = "featureGetPrimaryChildren";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_HAR, $ftr_idx, $AR, $FH_HR) = @_;

  if(! featureHasPrimaryChildren($ftr_info_HAR, $ftr_idx, $FH_HR)) { 
    DNAORG_FAIL("ERROR in $sub_name, feature $ftr_idx has no primary children", 1, $FH_HR);
  }

  @{$AR} = ();
  @{$AR} = split(/\s+/, $ftr_info_HAR->{"primary_children_ftr_str"}[$ftr_idx]);

  if(scalar(@{$AR}) == 0) { 
    DNAORG_FAIL("ERROR in $sub_name, feature $ftr_idx has no primary children (on second pass -- shouldn't happen)", 1, $FH_HR);
  }
  if(scalar(@{$AR}) != $ftr_info_HAR->{"primary_children_ftr_num"}[$ftr_idx]) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, feature $ftr_idx has %d children in ftr_info_HAR->{primary_children_ftr_str} but ftr_info_HAR->{primary_children_ftr_num} is %d", 
                        scalar(@{$AR}), $ftr_info_HAR->{"primary_children_ftr_num"}[$ftr_idx]), 1, $FH_HR); 
  }

  return;
}

#################################################################
# Subroutine : featureGetAllChildren()
# Incept:      EPN, Sat Feb 16 07:00:08 2019
#
# Purpose:     Fill @{$AR} with the space delimited tokens of 
#              @{$ftr_info_HAR->{"all_children_ftr_str"}[$ftr_idx]}
#              and return;
# 
# Arguments: 
#   $ftr_info_HAR:   REF to hash of arrays with information on the features, PRE-FILLED
#   $ftr_idx:        index of feature we're interested in
#   $AR:             REF to array to fill, FILLED HERE
#   $FH_HR:          REF to hash of file handles
# 
# Returns:     Nothing.
# 
# Dies: If there are no all children for feature index $ftr_idx,
#       If number of all children in $ftr_info_HAR->{"all_children_ftr_str"}
#       disagrees with $ftr_info_HAR->{"all_children_ftr_num"}.
#
################################################################# 
sub featureGetAllChildren { 
  my $nargs_expected = 4;
  my $sub_name = "featureGetAllChildren";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($ftr_info_HAR, $ftr_idx, $AR, $FH_HR) = @_;

  if(! featureHasAllChildren($ftr_info_HAR, $ftr_idx, $FH_HR)) { 
    DNAORG_FAIL("ERROR in $sub_name, feature $ftr_idx has no all children", 1, $FH_HR);
  }

  @{$AR} = ();
  @{$AR} = split(/\s+/, $ftr_info_HAR->{"all_children_ftr_str"}[$ftr_idx]);

  if(scalar(@{$AR}) == 0) { 
    DNAORG_FAIL("ERROR in $sub_name, feature $ftr_idx has no all children (on second pass -- shouldn't happen)", 1, $FH_HR);
  }
  if(scalar(@{$AR}) != $ftr_info_HAR->{"all_children_ftr_num"}[$ftr_idx]) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, feature $ftr_idx has %d children in ftr_info_HAR->{all_children_ftr_str} but ftr_info_HAR->{all_children_ftr_num} is %d", 
                        scalar(@{$AR}), $ftr_info_HAR->{"all_children_ftr_num"}[$ftr_idx]), 1, $FH_HR); 
  }

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
  addToErrorInfoHash($err_info_HAR, "n_nst", "feature",  1,
                     "no in-frame stop codon exists 3' of predicted valid start codon", # description
                     1, 1, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_nm3", "feature",  0,
                     "length of nucleotide feature is not a multiple of 3", # description
                     1, 0, "similar to !out_product,out_gene!; length is not a multiple of 3", # feature table info: valid, pred_stop, note
                     "Unexpected Length: (!out_product,out_gene!) length is not a multiple of 3; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_stp", "feature",  1,
                     "predicted CDS stop by homology is invalid; there may be a valid stop in a different location due to truncation (trc) or extension (ext) (TAG|TAA|TGA)", # description
                     1, 1, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified on !out_product,out_gene!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_trc", "feature",  0,
                     "in-frame stop codon exists 5' of stop position predicted by homology to reference", # description
                     # NOTE: invalidated by int and ctr because int or ctr will handle the feature table note/err
                     1, 1, "similar to !out_product,out_gene!; contains premature stop codon", # feature table info: valid, pred_stop, note
                     "!FEATURE_TYPE! Has Stop Codon: (!out_product,out_gene!) contains unexpected stop codon; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_ext", "feature",  1,
                     "first in-frame stop codon exists 3' of stop position predicted by homology to reference", # description
                     1, 1, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop note
                     "Mutation at End: (!out_product,out_gene!) expected stop codon could not be identified; !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_pe", "feature",  0,
                     "mat_peptide may not be translated because its CDS has an in-frame stop 5' of the mat_peptide's predicted start", # description
                     1, 0, "similar to !out_product,out_gene!; polyprotein may not be translated", # feature table info: valid, pred_stop, note
                     "Peptide Translation Problem: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_str", "feature",  0,
                     "predicted CDS start position is not beginning of start codon", # description
                     1, 0, "similar to !out_product,out_gene!; no start codon", # feature table info: valid, pred_stop, note
                     "Mutation at Start: (!out_product,out_gene!) expected start codon could not be identified", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_gp5", "feature",  0,
                     "alignment to reference does not extend to 5' boundary of reference or target", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "n_gp3", "feature",  0,
                     "alignment to reference does not extend to 3' boundary of reference or target", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_non", "feature",  0,
                     "blastx identifies protein not identified in nucleotide-based search", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_xnh", "feature",  0,
                     "blastx protein validation failure, no blastx hits", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_cst", "feature",  0,
                     "blastx protein validation failure, strand mismatch between protein and nucleotide predictions", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_p5l", "feature",  0,
                     "blastx protein validation failure, protein alignment extends past nucleotide alignment at 5' end", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_p5s", "feature",  0,
                     "blastx protein validation failure, protein alignment does not extend close enough to nucleotide alignment 5' endpoint", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation at Start: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_p3l", "feature",  0,
                     "blastx protein validation failure, protein alignment extends past nucleotide alignment at 3' end", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation at Stop: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_p3s", "feature",  0,
                     "blastx protein validation failure, protein alignment does not extend close enough to nucleotide alignment 3' endpoint", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation at Stop: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "p_lin", "feature",  0,
                     "blastx protein validation failure, too large of an insert", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Insertion of Nucleotides: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "p_lde", "feature",  0,
                     "blastx protein validation failure, too large of a deletion", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "Deletion of Nucleotides: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "p_trc", "feature",  0,
                     "blastx protein validation failure, stop codon in protein alignment", # description
                     1, 0, "similar to !out_product,out_gene!", # feature table info: valid, pred_stop, note
                     "!FEATURE_TYPE! Has Stop Codon: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "p_per", "feature",  0,
                     "mat_peptide may not be translated because its CDS has a blastx protein validation failure", # description
                     1, 0, "similar to !out_product,out_gene!; polyprotein may not be translated", # feature table info: valid, pred_stop, note
                     "Indefinite Annotation: (!out_product,out_gene!) !DESC!", # feature table error
                     $FH_HR);

  addToErrorInfoHash($err_info_HAR, "b_zft", "sequence", 0, # code, per-type, maybe-allowed
                     "zero features annotated", # description
                     1, 0, "", "No Features Annotated: (*sequence*) zero annotated features", # feature table info: valid, pred_stop, note, err
                     $FH_HR); 

  addToErrorInfoHash($err_info_HAR, "n_div", "sequence", 0, # code, per-type, maybe-allowed
                     "sequence too distant from reference to annotate", # description
                     1, 0, "", "Unexpected Divergence: (*sequence*) sequence is too divergent to confidently assign nucleotide-based annotation !DESC!", # feature table info: valid, pred_stop, note, err
                     $FH_HR); 

  # define the incompatibilities; these are two-sided, any error code listed in the 3rd arg is incompatible with the 2nd argument, and vice versa
  #setIncompatibilityErrorInfoHash($err_info_HAR, "str", "stp,trc,ext,b5e", $FH_HR);
  #setIncompatibilityErrorInfoHash($err_info_HAR, "trc", "ext,nst,b5e",     $FH_HR);
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
    runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  
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
        runCommand("rm $mp_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);
      }
      else { 
        # remove the empty file we just created
        runCommand("rm $mp_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);
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
    runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
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
    runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
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
  my %len_H = ();
  parseLengthFile($len_file, \%len_H, $FH_HR);
  my $nseq = 0;
  if(! $do_infasta) {  
    # we will have already stored all sequences in $seq_info_HAR->{"accn_name"}
    $nseq = scalar(@{$seq_info_HAR->{"accn_name"}});
    if($nseq == 0) { 
      DNAORG_FAIL("ERROR in $sub_name, no accessions in seq_info_HAR", 1, $FH_HR); 
    }
    for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
      if(! exists $len_H{$accn_name}) { 
        DNAORG_FAIL("ERROR in $sub_name, problem fetching length of accession $accn_name", 1, $FH_HR); 
      }
      $seq_info_HAR->{"len"}[$seq_idx] = $len_H{$accn_name};
    }
  }
  else { # $do_infasta is true, $seq_info_HAR->{"accn_name"} is empty, but we should have
         # read exactly 1 length in parseLengthFile
    $nseq = scalar(keys %len_H);
    if($nseq != 1) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, not exactly 1 accession read from length file $len_file (%d lengths read)", scalar(keys %len_H)), 1, $FH_HR); 
    }
    my ($ref_accn) = (keys %len_H);
    $seq_info_HAR->{"accn_name"}[0] = $ref_accn;
    $seq_info_HAR->{"len"}[0]       = $len_H{$ref_accn};
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
#              3) determines type of each reference sequence feature
#              4) fetches the reference sequence feature and populates information on the models and features
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
#   $in_ref_len:            length of the reference sequence 
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
  my $nargs_expected = 18;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name, entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $sqfile_R, $out_root, $build_root, $in_ref_accn, $in_ref_len, $infasta_file, $cds_tbl_HHAR, $mp_tbl_HHAR, $xfeat_tbl_HHHAR, $dfeat_tbl_HHHAR, $cds2pmatpept_AAR, $cds2amatpept_AAR, $mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;
  
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

  # $in_ref_accn, $in_ref_len must be defined if $do_infasta, and must not be defined if ! $do_infasta
  if($do_infasta && ((! defined $in_ref_accn) || (! defined $in_ref_len))) { 
    DNAORG_FAIL("ERROR in $sub_name, contract violation, --infasta option enabled but in_ref_accn and/or in_ref_len variables are not defined", 1, $FH_HR);
  }
  elsif((! $do_infasta) && ((defined $in_ref_accn || defined $in_ref_len))) { 
    DNAORG_FAIL("ERROR in $sub_name, contract violation, --infasta option NOT enabled but in_ref_accn and/or in_ref_len variables are defined", 1, $FH_HR);
  }
    
  # 1) fetch the sequences into a fasta file and index that fasta file
  my $nseq = scalar(@{$seq_info_HAR->{"accn_name"}});
  my $have_only_ref = ($nseq == 1) ? 1 : 0; # '1' if we only have the reference
  my $ref_accn      = ($do_infasta) ? $in_ref_accn : $seq_info_HAR->{"accn_name"}[0]; # the reference accession
  my $ref_len       = ($do_infasta) ? $in_ref_len  : $seq_info_HAR->{"len"}[0]; # length of the accession
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
    fetchSequencesUsingEslFetchCds($execs_HR->{"esl_fetch_cds"}, $fetch_file, $fasta_file, $seq_info_HAR, $FH_HR);
    addClosedFileToOutputInfo($ofile_info_HHR, "fetch", $fetch_file, 0, "Input file for esl-fetch-cds.pl");
    addClosedFileToOutputInfo($ofile_info_HHR, "fasta", $fasta_file, 0, "Sequence file (fetched with esl-fetch-cds.pl)");
  }
  else { # we should already have the fasta file that we need 
    validateFileExistsAndIsNonEmpty($fasta_file, $sub_name, $FH_HR);
    if($do_skipfetch) { 
      addClosedFileToOutputInfo($ofile_info_HHR, "fasta", $fasta_file, 0, "Sequence file (fetched on previous run (--skipfetch))");
    }
    else { # $do_infasta is true
      addClosedFileToOutputInfo($ofile_info_HHR, "fasta", $fasta_file, 0, "Sequence file (supplied with --infasta), NOT CREATED BY THE SCRIPT");
    }
  }  
  # make a duplicate of the fasta file with no descriptions, because blastx adds the 
  # description on lines immediately after the query, and sometimes breaks up sequence
  # name across multiple lines, so to be sure that we have the query name and only the
  # query name we need to read it from >=1 lines when blastx run on an input file with
  # no descriptions for each sequence. Example from blastx output:
  #Query= gi|337255678|gb|JF830576.1|
  #Norovirus/Hu/GII.4/757/10/17-Mar-2010/Slovenia strain 757/10
  #RNA-dependent RNA polymerase (RdRp) gene, partial cds
  my $fasta_nodesc_file = $out_root . ".nodesc.fa"; 
  open(FAOUT, ">", $fasta_nodesc_file)                || fileOpenFailure($fasta_nodesc_file, $sub_name, $!, "writing", $FH_HR);
  open(FAIN,  $ofile_info_HHR->{"fullpath"}{"fasta"}) || fileOpenFailure($ofile_info_HHR->{"fullpath"}{"fasta"}, $sub_name, $!, "reading", $FH_HR);
  while(my $fasta_line = <FAIN>) { 
    if($fasta_line =~ m/^\>(\S+).*/) { 
      print FAOUT (">" . $1 . "\n"); # only print sequence name, no desc, if there was one
    }
    else { 
      print FAOUT $fasta_line;
    }
  }
  close(FAIN);
  close(FAOUT);
  addClosedFileToOutputInfo($ofile_info_HHR, "fastanodesc", $fasta_nodesc_file, 0, "Sequence file with sequence descriptions removed for blastx");

  # open the sequence file using Bio-Easel
  # remove the .ssi file first, if it exists
  my $ssi_file = $fasta_file . ".ssi";
  if(-e $ssi_file) { 
    runCommand("rm $ssi_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);
  }
  $$sqfile_R = Bio::Easel::SqFile->new({ fileLocation => $fasta_file }); # the sequence file object

  if($do_skipfetch || $do_infasta) { 
    # --skipfetch or --infasta enabled, so we never called 
    # fetchSequencesUsingEslFetchCds() above.
    # We need to fill $seq_info_HAR->{"seq_name"} and $seq_info_HAR->{"len"}
    # get index hash for @{$seq_info_HAR->{"seq_accn"}} array
    # this simplifies determining sequence index in @{%seq_info_HAR->{}}
    # arrays for a given accession name.
    my %accn_name_idx_H = (); # key: $accn_name, value: idx of $accn_name in @{$seq_info_HAR->{"accn_name"}}
    getIndexHashForArray($seq_info_HAR->{"accn_name"}, \%accn_name_idx_H, $FH_HR);
    for(my $sqfile_seq_idx = 0; $sqfile_seq_idx < $nseq; $sqfile_seq_idx++) { 
      my ($seq_name, $seq_len) = $$sqfile_R->fetch_seq_name_and_length_given_ssi_number($sqfile_seq_idx);
      my $accn_name = undef;
      if($do_infasta) { 
        # if --infasta, then we are fetching directly from the --infasta file, else 
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
      $seq_info_HAR->{"len"}[$accn_name_idx_H{$accn_name}]      = $seq_len;
    }
  }

  # 2) determine reference information for each feature (strand, length, coordinates, product)
  getReferenceFeatureInfo($cds_tbl_HHAR, $mp_tbl_HHAR, $xfeat_tbl_HHHAR, $dfeat_tbl_HHHAR, $ftr_info_HAR, $ref_accn, $FH_HR); # $mp_tbl_HHAR may be undefined and that's okay
  my @reqd_keys_A = ("ref_strand", "ref_len", "ref_coords", "out_product", "out_gene", "out_exception", "type_fname", "type_ftable");
  validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, \@reqd_keys_A, $FH_HR);

  # 3) determine type of each reference feature
  my $ncds   = (defined $cds_tbl_HHAR) ? scalar(@{$cds_tbl_HHAR->{$ref_accn}{"coords"}}) : 0; # number of CDS features
  my $nmp    = (defined $mp_tbl_HHAR)  ? scalar(@{$mp_tbl_HHAR->{$ref_accn}{"coords"}})  : 0; # number of mature peptides
  my $nxfeat = (defined $xfeat_tbl_HHHAR) ? getNumExtraOrDuplicateFeatures($xfeat_tbl_HHHAR, $FH_HR) : 0;
  my $ndfeat = (defined $dfeat_tbl_HHHAR) ? getNumExtraOrDuplicateFeatures($dfeat_tbl_HHHAR, $FH_HR) : 0;
  determineFeatureTypes($ncds, $nmp, $nxfeat, $ndfeat, $cds2pmatpept_AAR, $cds2amatpept_AAR, $ftr_info_HAR, $FH_HR); # $cds2pmatpept_AAR may be undef and that's okay

  # 4) fetch the reference feature sequences and populate information on the models and features
  #    we won't actually fetch the reference sequence if --infasta is used (which is why $ref_seqname is undef if --infasta)
  my $ref_seqname = ($do_infasta) ? undef  : $seq_info_HAR->{"seq_name"}[0]; # the reference sequence name the fetched sequence file $fasta_file
  my $stk_file    = $out_root . ".ref.stk";     # name of output alignment file we are about to create, with the single full sequence as the 'alignment'
  fetchReferenceFeatureSequences($execs_HR, $$sqfile_R, $ref_seqname, $ref_len, $out_root, $build_root, $mdl_info_HAR, $ftr_info_HAR, $stk_file, $opt_HHR, $FH_HR); 
  if(! $do_infasta) { 
    addClosedFileToOutputInfo($ofile_info_HHR, "refstk", $stk_file, 0, "Stockholm alignment file with reference sequence");
  }

  # 5) determine the source "dupsource" of any duplicate feature annotations, for
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
    # segments, because the order of reverse strand segments in a feature table is 
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
  # now add complement for cases where are segments/segments are on reverse strand
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
  foreach my $opt ("--nomatpept", "--matpept", "--xfeat", "--dfeat") { 
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
#                                 "duplicate":    if this feature's annotation is derived by copying data from another (e.g., type: 'gene')
#                "filename_root": 'root' string for output file names related to this model: 
#                "final_mdl":     index (in arrays of %mdl_info_HA) of final model for this feature
#                "first_mdl":     index (in arrays of %mdl_info_HA) of first model for this feature
#                "nmodels":       number of models for this feature (e.g. number of segments) for this feature, 
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
#                "type":          type of feature: e.g. "cds"
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
  
  my @expected_keys_A = ("final_mdl", "first_mdl", "annot_type", "nmodels", 
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
#                "checksum":      checksum of the 'alignment' (single sequence) file the model was built from
#                "cmname":        name of the model, used in infernal output 
#                "filename_root": 'root' string for output file names related to this model: 
#                "is_final":      '1' if this model is the final model (e.g. final segment) for the feature it models ("map_ftr")
#                "is_first":      '1' if this model is the first model (e.g. final segment) for the feature it models ("map_ftr")
#                "length":        length, in nucleotides, of the model
#                "ref_start":     start position of modelled region in the reference genome
#                "ref_stop":      stop position of modelled region in the reference genome
#                "ref_strand":    strand of modelled region in the reference genome
#                "map_segment":   the segment index this model models (1.."map_nsegment" value)
#                "map_ftr":       the feature index (array index in ftr_info_HAR) this model models
#                "map_nsegment":  the number of segments the feature this model models has 
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
  
  my @expected_keys_A = ("checksum", "cmname", "filename_root", "is_final", "is_first", 
                         "length", "ref_start", "ref_stop", "ref_strand", "map_segment", 
                         "map_ftr", "map_nsegment", "out_tiny", "out_idx");

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
#             $short_ok: '1' if it's okay to return a codon less than 3 nucleotides
#                        if the desired codon runs off the end of the sequence
#             $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#             
# Returns:    The stop codon as a string
# 
# Dies:       If $stop (or $stop - 2) is negative.
#################################################################
sub fetchStopCodon {
  my $sub_name = "fetchStopCodon";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seq_name, $stop, $strand, $short_ok, $FH_HR) = @_;

  my $seqlen = $sqfile->fetch_seq_length_given_name($seq_name);

  # printf("\nin $sub_name, $seq_name stop:$stop strand:$strand short_ok: $short_ok\n");

  if($seqlen == -1) { 
    DNAORG_FAIL("ERROR in $sub_name, unable to fetch sequence $seq_name", 1, $FH_HR);
  }

  if(($stop < 1) || ($stop > $seqlen)) { 
    DNAORG_FAIL("ERROR in $sub_name, input stop coordinate ($stop) invalid (length of sequence ($seq_name) is $seqlen", 1, $FH_HR); 
  }    

  my $stop_codon_posn;
  if($strand eq "-") { 
    $stop_codon_posn = $stop + 2;
    if($stop_codon_posn > $seqlen) { $stop_codon_posn = $seqlen; }
  }
  else { 
    $stop_codon_posn = $stop - 2;
    if($stop_codon_posn < 1) { $stop_codon_posn = 1; }
  }
  # printf("in $sub_name, seq_name $seq_name, stop $stop\n");

  return fetchCodon($sqfile, $seq_name, $stop_codon_posn, $strand, $short_ok, $FH_HR);
}

#################################################################
# Subroutine: fetchStartCodon()
# Incept:     EPN, Tue Mar 15 10:18:34 2016
#
# Synopsis:   Fetch a start codon given its first position,
#             and strand.
#
# Args:       $sqfile:   Bio::Easel::SqFile object, open sequence
#                        file containing $seq_name;
#             $seq_name: name of sequence to fetch part of
#             $stop:     final position of the stop codon
#             $strand:   strand we want ("+" or "-")
#             $short_ok: '1' if it's okay to return a codon less than 3 nucleotides
#                        if the desired codon runs off the end of the sequence
#             $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#             
# Returns:    The stop codon as a string
# 
# Dies:       If $start is negative.
#################################################################
sub fetchStartCodon {
  my $sub_name = "fetchStartCodon";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args (expected $nargs_exp)"; }

  my ($sqfile, $seq_name, $start, $strand, $short_ok, $FH_HR) = @_;

  my $seqlen = $sqfile->fetch_seq_length_given_name($seq_name);

  # printf("\nin $sub_name, $seq_name start:$start strand:$strand short_ok: $short_ok\n");

  if(($start < 1) || ($start > $seqlen)) { 
    DNAORG_FAIL("ERROR in $sub_name, input start coordinate ($start) invalid (length of sequence ($seq_name) is $seqlen", 1, $FH_HR); 
  }    

  return fetchCodon($sqfile, $seq_name, $start, $strand, $short_ok, $FH_HR);
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
#                        file containing $seq_name;
#             $seq_name: name of sequence to fetch part of
#             $start:    start position of the codon
#             $strand:   strand we want ("+" or "-")
#             $short_ok: '1' if it's okay to return a codon less than 3 nucleotides
#                        if the desired codon runs off the end of the sequence
#             $FH_HR:    REF to hash of file handles, including "log" and "cmd"
#
# Returns:    The codon as a string
#
#################################################################
sub fetchCodon {
  my $sub_name = "fetchCodon";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seq_name, $start, $strand, $short_ok, $FH_HR) = @_;

  my $seqlen = $sqfile->fetch_seq_length_given_name($seq_name);

  # printf("\nin $sub_name, $seq_name start:$start strand:$strand short_ok: $short_ok\n");

  if($seqlen == -1) { 
    DNAORG_FAIL("ERROR in $sub_name, unable to fetch sequence $seq_name", 1, $FH_HR);
  }

  if(($start < 1) || ($start > $seqlen)) { 
    DNAORG_FAIL("ERROR in $sub_name, input start coordinate ($start) invalid (length of sequence ($seq_name) is $seqlen", 1, $FH_HR); 
  }    

  my $codon_start = $start;
  my $codon_stop  = undef;
  my $valid_coords = 1; # changed to 0 below if coords are invalid

  if($strand eq "+") { 
    $codon_stop  = $start + 2; 
    if($codon_stop > $seqlen) { 
      if($short_ok) { $codon_stop = $seqlen; }
      else          { $valid_coords = 0; }
    }
  }
  else { # strand is '-'
    $codon_stop = $start - 2;
    if($codon_stop < 1) { 
      if($short_ok) { $codon_stop = 1; } 
      else          { $valid_coords = 0; }
    }
  }
  if(! $valid_coords) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, trying to fetch codon with coordinates %d..%d for sequence $seq_name, valid coordinates are %d..%d\n", 
                        $codon_start, $codon_stop, 1, $seqlen), 1, $FH_HR);
  }

  my $newname = $seq_name . "/" . $codon_start . "-" . $codon_stop;

  my @fetch_AA = ();
  # if $seq_name does not exist in $sqfile, the following line
  # will cause the script to fail, ungracefully

  push(@fetch_AA, [$newname, $codon_start, $codon_stop, $seq_name]);

  my $faseq = $sqfile->fetch_subseqs(\@fetch_AA, -1);

  my ($header, $seq) = split("\n", $faseq);

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
#   $seq_info_HAR:     REF to 2D hash with sequence information, ADDED TO HERE
#                      "seq_name", and "len" arrays filled here
#   $FH_HR:            REF to hash of file handles, including "log" and "cmd"
# 
# Returns:    void
#          
# Dies: if $esl_fetch_cds executable does not exist
#       if $esl_fetch_cds command fails
#################################################################
sub fetchSequencesUsingEslFetchCds { 
  my $sub_name = "fetchSequencesUsingEslFetchCds";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($esl_fetch_cds, $fetch_file, $fasta_file, $seq_info_HAR, $FH_HR) = @_;

  # contract check
  if(! -e $esl_fetch_cds) { 
    DNAORG_FAIL("ERROR in $sub_name, required executable $esl_fetch_cds does not exist", 1, $FH_HR); 
  }

  my ($accn_name, $seq_name, $seq_len); # accession name/length (from GenBank) sequence name/length in the sequence file we create

  @{$seq_info_HAR->{"seq_name"}} = (); # initialize
  my $seq_name_AR = \@{$seq_info_HAR->{"seq_name"}};

  my $fetch_string; # string output to the input file for esl-fetch-cds.pl
  
  # create the esl-fetch-cds input file
  open(OUT, ">", $fetch_file) || fileOpenFailure($fetch_file, $sub_name, $!, "writing", $FH_HR);
  my $naccn = scalar(@{$seq_info_HAR->{"accn_name"}});
  for(my $a = 0; $a < $naccn; $a++) { 
    my $accn_name = $seq_info_HAR->{"accn_name"}[$a];
    my $seq_len   = $seq_info_HAR->{"len"}[$a];

    $fetch_string = $accn_name . ":1.." . $seq_len . "\n";
    print OUT $accn_name . ":" . "dnaorg" . "\t" . $fetch_string;
    $seq_name = $accn_name . ":dnaorg:" . $accn_name . ":1:" . $seq_len . ":+:";

    push(@{$seq_info_HAR->{"seq_name"}}, $seq_name);
  }
  close(OUT);
  
  # execute $esl_fetch_cds to fetch it
  my $cmd = "perl $esl_fetch_cds -nocodon $fetch_file > $fasta_file";
  runCommand($cmd, 0, 0, $FH_HR);

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
    runCommand($cmpress_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

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
    runCommand($cat_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
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
    my $cds_nsegments = 0;
    startsStopsStrandsFromCoordsLength($ref_cds_coords_A[$cds_idx], $ref_totlen, $do_circular, \@cds_starts_A, \@cds_stops_A, \@cds_strands_A, \$cds_nsegments, $FH_HR);
    my $cds_start = $cds_starts_A[0];
    my $cds_stop  = $cds_stops_A[$cds_nsegments-1];
    if($cds_nsegments != 1) { 
      if($cds_nsegments != 2) { 
        DNAORG_FAIL("ERROR in $sub_name, triple or more segment CDS broken up to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
      }
      if($cds_strands_A[0] ne $cds_strands_A[1]) { 
        DNAORG_FAIL("ERROR in $sub_name, double segment CDS with each segment on different strands to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
      }
      # two segment CDS, if any introns exist (any nt is not included between $cds_start..$cds_stop) then we can't handle it
      # example of a multi-'segment' CDS that we CAN handle is West Nile Virus CDS #2: NC_001563.2 join(97..3540,3540..3671)
      if($cds_strands_A[0] eq "+") { 
        if(($cds_starts_A[1] - $cds_stops_A[0] - 1) > 0) { 
#          DNAORG_FAIL("ERROR in $sub_name, multiple segment CDS with an intron broken up to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
        }
      }
      else { # negative strand
        if(($cds_stops_A[0] - $cds_starts_A[1] - 1) > 0) { 
#          DNAORG_FAIL("ERROR in $sub_name, multiple segment CDS with an intron broken up to make mat_peptides, code for this does not yet exist.", 1, $FH_HR);
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
      my $mp_nsegments   = 0;
      startsStopsStrandsFromCoordsLength($ref_mp_coords_A[$mp_idx], $ref_totlen, $do_circular, \@mp_starts_A, \@mp_stops_A, undef, \$mp_nsegments, $FH_HR);
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
        if(($mp_stops_A[($mp_nsegments-1)]+3) != $cds_stop) { 
          DNAORG_FAIL(sprintf("ERROR in $sub_name, for cds_idx $cds_idx stop of final mat_peptide doesn't match CDS stop (%d != %d)", $mp_stops_A[($mp_nsegments-1)], $cds_stop), 1, $FH_HR);
        }
      }
      $prv_stop = $mp_stops_A[($mp_nsegments-1)];
      # printf("checked mp $mp_idx %d..%d\n", $mp_starts_A[0], $mp_stops_A[($mp_nsegments-1)]);
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
#   $starts_AR:   ref of array of start positions to potentially overwrite (if we find this is really only one segment)
#   $stops_AR:    ref of array of stop positions to potentially overwrite (if we find this is really only one segment)
#   $nsegments_R: ref to scalar of number of segments to overwrite (if we find this is really only one segment)
#   $do_update:   '1' to update $starts_AR, $stops_AR and $nsegments_R if we find two segments that span stop..start
#   $strand:      strand the segments are on
#   $totlen:      total length of the sequence
#           
# Returns:    '1' if we found two segments that spanned stop..start 
#             (and if ($do_update) then we also updated @{$starts_AR}, @{$stops_AR} and $$nsegments_R)
# 
# Dies:       Never.
# 
#################################################################
sub checkForSpanningSequenceSegments {
  my $sub_name = "checkForSpanningSequenceSegments()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($starts_AR, $stops_AR, $nsegments_R, $do_update, $strand, $totlen) = @_;

  my $found_spanning_segments = 0;

  # if $$nsegments_R is not '2', we skip this block and return '0'.
  # 'spanning' segments can only occur if there are exactly 2 segments that
  # we are checking. 
  if($$nsegments_R == 2) { 
    # if we're in a circular genome, we need to check for a special case, where 
    # what looks like a 2-segment CDS is really a single segment that spans the stop..start boundary.
    # [Note that if the stop..start boundary is spanned by an intron (i.e. segment i is before stop,
    # and i+1 is after start) then we don't need to modify anything, we'll still fetch the proper
    # sequence even in a duplicated genome].
    #
    # Example 1: single segment that spans stop..start boundary on positive strand
    # join(2309..3182,1..1625) in a seq of length 3182, this should really be a single segment
    # $starts_A[0] = 2309
    # $stops_A[0]  = 3182
    # $starts_A[1] = 1
    # $stops_A[1]  = 1625
    # $nsegments = 2;
    # 
    # should become:
    # $starts_A[0] = 2309
    # $stops_A[0]  = 4807
    # $nsegments = 1;
    # 
    # Example 2: single segment that spans stop..start boundary on negative strand
    # complement(join(2309..3182,1..1625))   in a seq of length 3182, this should really be a single segment
    # $starts_A[0] = 3182
    # $stops_A[0]  = 2309
    # $starts_A[1] = 1625
    # $stops_A[1]  = 1
    # $nsegments = 2;
    # 
    # should become:
    # $starts_A[0] = 4807
    # $stops_A[0]  = 2309
    # $nsegments = 1;
    #
    # we can easily check and fix these cases:
    my $tmp_start = undef;
    my $tmp_stop  = undef;
    # remember if we get here, we know we only have 2 segments, i.e scalar(@{$starts_AR}) and scalar(@{$stops_AR}) is 2
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
      $$nsegments_R   = 1;
      $found_spanning_segments = 1;
    }    
  }
  return $found_spanning_segments;
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
#             'cds'        1
#             other:       0
#
# Arguments:
#   $feature_type:  input feature
#             
# Returns:    1 or 0
#
#################################################################
sub featureTypeIsCds { 
  my $sub_name  = "featureTypeIsCds";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_feature) = (@_);

  if($in_feature eq "cds") { 
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
#             'mp":        1
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
#             'dfeat':     1
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
#  $ftr_info_HAR:   ref to the feature info hash of arrays 
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

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $FH_HR); # nftr: number of features

  my $ret_ftr_idx = undef;
  if($blastx_seqname =~ /(\S+)\/(\S+)/) { 
    my ($accn, $coords) = ($1, $2);
    # find it in @{$ftr_info_HAR->{"ref_coords"}}
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if(($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds")) { 
        if($ftr_info_HAR->{"ref_coords"}[$ftr_idx] eq $coords) { 
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
# Subroutine:  countFeatureType()
# Incept:      EPN, Sat Feb 23 07:29:49 2019
#
# Purpose:    Return number of features of type $type.
#             Does not check that $type is a valid type.
#
# Arguments: 
#  $ftr_info_HAR:   ref to the feature info hash of arrays 
#  $type:           feature 'type'
#
# Returns:    Number of features of type $type.
#
# Dies:       Never (does not validate $ftr_info_HAR or enforce that $type 
#             is a valid type)
#
################################################################# 
sub countFeatureType { 
  my $sub_name = "countFeatureType";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $type) = @_;

  my $ret_n = 0;
  foreach my $ftype (@{$ftr_info_HAR->{"type"}}) { 
    if($ftype eq $type) { $ret_n++; }
  }
  return $ret_n;
}

#################################################################
# Subroutine:  countFeatureTypeAndStrand()
# Incept:      EPN, Sat Feb 23 07:32:39 2019
#
# Purpose:    Return number of features of type $type
#             *and* ref_strand $strand.
#             Does not check that $type is a valid type
#             or $strand is a valid strand value.
#
# Arguments: 
#  $ftr_info_HAR:   ref to the feature info hash of arrays 
#  $type:           feature 'type'
#  $strand:         feature 'ref_strand'
#  $FH_HR:          ref to hash of file handles
#
# Returns:    Number of features of type $type.
#
# Dies:       If number of "type" values differs from number of 
#             "ref_strand" values.
#
################################################################# 
sub countFeatureTypeAndStrand { 
  my $sub_name = "countFeatureTypeAndStrand";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $type, $strand, $FH_HR) = @_;

  my @type_A   = @{$ftr_info_HAR->{"type"}};
  my @strand_A = @{$ftr_info_HAR->{"ref_strand"}};

  my $ntype   = scalar(@type_A);
  my $nstrand = scalar(@strand_A);
  if($ntype != $nstrand) { 
    DNAORG_FAIL("ERROR in $sub_name, number of types ($ntype) != number of strands ($nstrand)", 1, $FH_HR); 
  }
  my $ret_n = 0;
  for(my $i = 0; $i < $ntype; $i++) { 
    if(($type_A[$i]   eq $type) &&
       ($strand_A[$i] eq $strand)) { 
      $ret_n++;
    }
  }
  return $ret_n;
}

#################################################################
# Subroutine:  checkIfFeatureIsCdsOrMp()
# Incept:      EPN, Mon Feb 25 14:30:34 2019
#
# Purpose:    Return '1' if feature is type is 'cds' or 'mp', else return '0'.
#
# Arguments: 
#  $ftr_info_HAR:   ref to the feature info hash of arrays 
#  $ftr_idx:        feature index
#
# Returns:    '1' if $ftr_info_HAR->{"type"}[$ftr_idx] is "cds" or "mp"
#             else '0'
#
# Dies:       never; does not validate anything.
#
################################################################# 
sub checkIfFeatureIsCdsOrMp { 
  my $sub_name = "checkIfFeatureIsCdsOrMp";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $ftr_idx) = @_;

  if(($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds") ||
     ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp")) { 
    return 1; 
  }
  return 0;
}

#################################################################
# Subroutine:  checkIfFeatureIsCds()
# Incept:      EPN, Wed Feb 27 10:38:34 2019
#
# Purpose:    Return '1' if feature type is 'cds', else return '0'.
#
# Arguments: 
#  $ftr_info_HAR:   ref to the feature info hash of arrays 
#  $ftr_idx:        feature index
#
# Returns:    '1' if $ftr_info_HAR->{"type"}[$ftr_idx] is "cds"
#             else '0'
#
# Dies:       never; does not validate anything.
#
################################################################# 
sub checkIfFeatureIsCds { 
  my $sub_name = "checkIfFeatureIsCds";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $ftr_idx) = @_;

  if($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds") {
    return 1; 
  }
  return 0;
}

#################################################################
# Subroutine:  checkIfFeatureIsDuplicate()
# Incept:      EPN, Mon Feb 25 14:30:34 2019
#
# Purpose:    Return '1' if feature type is 'dfeat', else return '0'.
#
# Arguments: 
#  $ftr_info_HAR:   ref to the feature info hash of arrays 
#  $ftr_idx:        feature index
#
# Returns:    '1' if $ftr_info_HAR->{"type"}[$ftr_idx] is "dfeat"
#             else '0'
# 
# Dies:       never; does not validate anything.
#
################################################################# 
sub checkIfFeatureIsDuplicate { 
  my $sub_name = "checkIfFeatureIsDuplicate"; 
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $ftr_idx) = @_;

  if($ftr_info_HAR->{"type"}[$ftr_idx] eq "dfeat") { 
    return 1; 
  }
  return 0;
}

###########################################################################
# the next line is critical, a perl module must return a true value
    return 1;
###########################################################################
