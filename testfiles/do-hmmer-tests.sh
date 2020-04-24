#!/bin/bash

RETVAL=0;

# noro r10 hmmer
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.hmmer.testin vt-n10-hmmer
if [ $? != 0 ]; then
   RETVAL=1;
fi   

# dengue r5 hmmer
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.hmmer.testin vt-d5-hmmer
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
