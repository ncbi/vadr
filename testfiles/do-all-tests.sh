#!/bin/bash

RETVAL=0;

# Run tests in t/ subdir
for t in \
    do-prove-all-tests.sh \
    ; do
    sh $VADRSCRIPTSDIR/t/$t teamcity
    if [ "$?" -ne 0 ]; then
        RETVAL=1;
    fi   
done

# Run sequip tests in sequip/t/ subdir
for t in \
    do-prove-all-tests-from-vadr.sh \
    ; do
    sh $VADRSEQUIPDIR/t/$t teamcity
    if [ "$?" -ne 0 ]; then
        RETVAL=1;
    fi   
done

# If you want to test -p option for parallelization, add
# do-install-tests-parallel.sh to the following for loop.
# Note: this test requires qsub is in your path and qsub options are
# configured similarly to ncbi cluster, email eric.nawrocki@nih.gov
# for information on how to configure for different clusters
for t in \
    do-install-tests-local.sh \
    do-fs-tests.sh \
    do-replace-tests.sh \
    do-seed-tests.sh \
    do-hmmer-tests.sh \
    do-nindel-tests.sh \
    do-outaln-tests.sh \
    do-mxsize-tests.sh \
    do-uj-tests.sh \
    do-noftr-tests.sh \
    do-nends-tests.sh \
    do-ftskipfl-tests.sh \
    do-sub-tests.sh \
    do-lowsim-tests.sh \
    do-mnf-tests.sh \
    do-dcr-tests.sh \
    do-glsearch-tests.sh \
    do-split-tests.sh \
    do-af-tests.sh \
    do-indfstrn-tests.sh \
    do-dupregin-tests.sh \
    do-minimap2-tests.sh \
    do-lowsimexc-tests.sh \
    do-ss-tests.sh \
    do-extrant-tests.sh \
    github-issues/do-issue-tests.sh \
    ; do
    sh $VADRSCRIPTSDIR/testfiles/$t
    if [ "$?" -ne 0 ]; then
        RETVAL=1;
    fi   
done

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-all-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-all-tests.sh]"
   exit 1
fi
