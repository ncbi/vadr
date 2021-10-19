#!/bin/bash

RETVAL=0;

# entoy100a-af1
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-af.testin entoy100a-af
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-af-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-af-tests.sh]"
   exit 1
fi
