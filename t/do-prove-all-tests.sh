#!/bin/bash

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
; do
    if [ $do_teamcity == 1 ]; then
        echo "##teamcity[testStarted name=\"$test\" captureStandardOutput='true']"
    fi

    prove -v $VADRSCRIPTSDIR/t/$test;
    ret_val=$?

    if [ $do_teamcity == 1 ]; then 
        if [ $ret_val != 0 ]; then
            echo "##teamcity[testFailed name=\"$test\" message=\"v-test.pl failure\"]"
        fi
        echo "##teamcity[testFinished name=\"$test\"]"
    fi
done
