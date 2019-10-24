#!/bin/bash

# make sure correct number of cmdline arguments were used, exit if not
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <path to v-annotate.pl output directory>"
  exit 1
fi

# make sure the specified v-annotate.pl output directory exists
if [ ! -d $1 ]; then
   echo "ERROR: specified v-annotate.pl output directory $1 does not exist"
   exit 1
fi

# to get list of errors: v-annotate.pl --alt_list | grep -v ^\# | awk '{ print $4 }' | sort | uniq | awk '{ printf("%s \\\n", $1); }'
for e in \
BIASED_SEQUENCE \
CDS_HAS_STOP_CODON \
DELETION_OF_NT \
DISCONTINUOUS_SIMILARITY \
DUPLICATE_REGIONS \
INCORRECT_SPECIFIED_GROUP \
INCORRECT_SPECIFIED_SUBGROUP \
INDEFINITE_ANNOTATION \
INDEFINITE_ANNOTATION_END \
INDEFINITE_ANNOTATION_START \
INDEFINITE_CLASSIFICATION \
INDEFINITE_STRAND \
INSERTION_OF_NT \
LOW_COVERAGE \
LOW_FEATURE_SIMILARITY \
LOW_FEATURE_SIMILARITY_END \
LOW_FEATURE_SIMILARITY_START \
LOW_SCORE \
LOW_SIMILARITY \
LOW_SIMILARITY_END \
LOW_SIMILARITY_START \
MUTATION_AT_END \
MUTATION_AT_START \
NO_ANNOTATION \
NO_FEATURES_ANNOTATED \
PEPTIDE_ADJACENCY_PROBLEM \
PEPTIDE_TRANSLATION_PROBLEM \
QUESTIONABLE_SPECIFIED_GROUP \
QUESTIONABLE_SPECIFIED_SUBGROUP \
REVCOMPLEM \
UNEXPECTED_DIVERGENCE \
UNEXPECTED_LENGTH \
; do 
    grep -v ^\# $1/*.alt | grep $e > $1/$e.alt.txt
    if [ -s $1/$e.alt.txt ]
    then
        cat $1/$e.alt.txt | awk '{ print $2 }' | sort | uniq > $1/$e.alt.acclist.txt
    else
        rm $1/$e.alt.txt
    fi
done


