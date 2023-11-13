#!/bin/bash

RETVAL=0;

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.noftr.testin vt-noro-noftr
if [ "$?" -eq 0 ]; then
   echo "Success: all tests passed [do-noftr-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-noftr-tests.sh]"
   exit 1
fi
