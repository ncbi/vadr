#!/bin/bash

RETVAL=0;

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.ftskipfl.testin vt-noro.ftskipfl
if [ "$?" -eq 0 ]; then
   echo "Success: all tests passed [do-ftskipfl-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-ftskipfl-tests.sh]"
   exit 1
fi
