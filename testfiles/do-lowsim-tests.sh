#!/bin/bash

RETVAL=0;

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-lowsim.testin vt-entoy100a-lowsim
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-lowsim-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-lowsim-tests.sh]"
   exit 1
fi
