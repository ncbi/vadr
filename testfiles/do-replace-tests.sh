#!/bin/bash

RETVAL=0;

# noro r10 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.rpn.testin vt-n10-rpn
if [ $? != 0 ]; then
   RETVAL=1;
fi   
   
# dengue r5 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.rpn.testin vt-d5-rpn
if [ $? != 0 ]; then
   RETVAL=1;
fi   

if [ $RETVAL == 0 ]; then
   echo "Success: all tests passed"
   exit 0
else 
   echo "FAIL: at least one test failed"
   exit 1
fi
