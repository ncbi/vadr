#!/bin/bash

RETVAL=0;

# entoy100a-nindel
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-minfo.testin entoy100a-minfo
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-minfo-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed do-minfo-tests.sh]"
   exit 1
fi
