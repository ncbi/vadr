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

for test in \
    noro.r10.local.testin \
    dengue.r5.local.testin \
    noro.fail.testin \
; do
    if [ $do_teamcity == 1 ]; then
        echo "##teamcity[testStarted name='$test' captureStandardOutput='true']"
    fi

    # actual v-test call
    $VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/$test vt-$test
    ret_val=$?

    if [ $do_teamcity == 1 ]; then 
        if [ $ret_val != 0 ]; then
            echo "##teamcity[testFailed name='$test' message='v-test.pl failure']"
        fi
        echo "##teamcity[testFinished name='$test']"
    fi
done


