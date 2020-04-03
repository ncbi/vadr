#!/bin/bash
set -e 

if [ ! "$BASH_VERSION" ] ; then
    echo "Please do not use sh to run this script ($0), just execute it directly" 1>&2
    exit 1
fi

if [ -z "${VADRSCRIPTSDIR}" ] ; then
    echo "VADRSCRIPTSDIR environment variable is not set, set it as instructed during VADR installation and rerun"
    exit 1
fi

$VADRSCRIPTSDIR/testfiles/github-issues/iss3/do-iss3-tests.sh
