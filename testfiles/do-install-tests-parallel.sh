#!/bin/bash

RETVAL=0;

# noro r10
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.parallel.testin n10-parallel
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# dengue r5
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.parallel.testin d5-parallel
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-install-tests-parallel.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-install-tests-parallel.sh]"
   exit 1
fi
