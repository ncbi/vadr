#!/bin/bash

RETVAL=0;

# noro r10
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.split.testin vt-n10-split
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-split-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-split-tests.sh]"
   exit 1
fi
