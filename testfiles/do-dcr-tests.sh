#!/bin/bash

RETVAL=0;

# entoy100a-dcr
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-dcr.testin vt-entoy100a-dcr
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-dcr-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-dcr-tests.sh]"
   exit 1
fi
