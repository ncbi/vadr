#!/bin/bash

RETVAL=0;
CURRETVAL=0;

# we require 0 or 1 args, if 1 arg, it must be 'teamcity'
do_teamcity=0
if [ $# != 0 ]; then
    if [ $# -gt 1 ] || [ $1 != "teamcity" ]; then 
        echo "Usage:"
        echo "$0"
        echo "OR"
        echo "$0 teamcity"
        exit 1;
    fi
    # if we get here, there's 1 arg and it's 'teamcity'
    do_teamcity=1;
fi

if [ -z "${VADRSCRIPTSDIR}" ] ; then
    echo "VADRSCRIPTSDIR environment variable is not set, set it to top-level vadr/ dir and rerun"
    exit 1
fi
if [ -z "${VADRINSTALLDIR}" ] ; then
    echo "VADRINSTALLDIR environment variable is not set, set it to top-level vadr/ dir and rerun"
    exit 1
fi
if [ -z "${VADRMODELDIR}" ] ; then
    echo "VADRMODELDIR environment variable is not set, set it to top-level vadr/ dir and rerun"
    exit 1
fi

for test in \
    01-coords.t \
    02-seed.t \
    03-parseblast.t \
    04-ifile.t \
    05-annot-mdlopts.t \
    06-iss12-longseqnames.t \
    07-iss22-modelname-parantheses.t \
; do
    if [ $do_teamcity == 1 ]; then
        echo "##teamcity[testStarted name=\"$test\" captureStandardOutput='true']"
    fi

    prove -v $VADRSCRIPTSDIR/t/$test;
    CURRETVAL=$?

    if [ $do_teamcity == 1 ]; then 
        if [ $CURRETVAL != 0 ]; then
            echo "##teamcity[testFailed name=\"$test\" message=\"v-test.pl failure\"]"
        fi
        echo "##teamcity[testFinished name=\"$test\"]"
    fi

    if [ $CURRETVAL != 0 ]; then
        RETVAL=1
    fi
done

if [ $RETVAL == 0 ]; then
   echo "Success: all tests passed"
   exit 0
else 
   echo "FAIL: at least one test failed"
   exit 1
fi
