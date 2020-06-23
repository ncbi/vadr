#!/bin/bash

RETVAL=0;

for i in \
    iss3 \
    iss4 \
    iss5 \
    iss6 \
    iss13 \
    ; do
    $VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/github-issues/$i/$i.testin vt-$i
    if [ $? != 0 ]; then
        RETVAL=1;
    fi   
done

if [ $RETVAL == 0 ]; then
   echo "Success: all tests passed"
   exit 0
else 
   echo "FAIL: at least one test failed"
   exit 1
fi
