#!/bin/bash

RETVAL=0;

# noro.indfstrn
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.indfstrn.testin noro.indfstrn
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-indfstrn-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-indfstrn-tests.sh]"
   exit 1
fi
