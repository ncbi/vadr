#!/bin/bash

RETVAL=0;

# noro r10 replace
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.rpn.testin vt-n10-rpn
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# noro r2 replace
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r2.rpn.testin vt-n2-rpn
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   
   
# dengue r5 replace
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.rpn.testin vt-d5-rpn
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-replace-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed do-replace-tests.sh]"
   exit 1
fi
