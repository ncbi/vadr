#!/bin/bash

RETVAL=0;

# noro r10 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.minimap2.testin vt-n10-minimap2
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# dengue r5 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.minimap2.testin vt-d5-minimap2
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-minimap2-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-minimap2-tests.sh]"
   exit 1
fi
