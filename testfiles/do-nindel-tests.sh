#!/bin/bash

RETVAL=0;

# entoy100a-nindel
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-nindel.testin entoy100a-nindel
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-nindel-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed do-nindel-tests.sh]"
   exit 1
fi
