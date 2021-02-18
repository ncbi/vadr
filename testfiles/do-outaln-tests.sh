#!/bin/bash

RETVAL=0;

# noro r1 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r3.outaln.testin vt-n3-outaln
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-outaln-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-outaln-tests.sh]"
   exit 1
fi

