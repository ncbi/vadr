#!/bin/bash

RETVAL=0;

# noro r10 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.seed.testin vt-n10-seed
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# dengue r5 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.seed.testin vt-d5-seed
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-seed-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-seed-tests.sh]"
   exit 1
fi
