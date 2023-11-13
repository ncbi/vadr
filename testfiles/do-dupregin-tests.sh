#!/bin/bash

RETVAL=0;

# noro.dupregin
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.dupregin.testin noro.dupregin
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-dupregin-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-dupregin-tests.sh]"
   exit 1
fi
