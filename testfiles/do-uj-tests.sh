#!/bin/bash

RETVAL=0;

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.uj.testin vt-noro-uj
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-uj-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-uj-tests.sh]"
   exit 1
fi
