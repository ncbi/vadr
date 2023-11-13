#!/bin/bash

RETVAL=0;

# noro.lowsimexc
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.lowsimexc.testin noro.lowsimexc
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-lowsimexc-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-lowsimexc-tests.sh]"
   exit 1
fi
