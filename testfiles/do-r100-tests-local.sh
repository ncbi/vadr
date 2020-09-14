#!/bin/bash

RETVAL=0;

# noro r100
$VADRSCRIPTSDIR/v-test.pl -f $VADRSCRIPTSDIR/testfiles/noro.r100.local.testin n100-local
if [ $? != 0 ]; then
   RETVAL=1;
fi   

# dengue r100
$VADRSCRIPTSDIR/v-test.pl -f $VADRSCRIPTSDIR/testfiles/dengue.r100.local.testin d100-local
if [ $? != 0 ]; then
   RETVAL=1;
fi   

if [ $RETVAL = 0 ]; then
   echo "Success: all tests passed"
   exit 0
else 
   echo "FAIL: at least one test failed"
   exit 1
fi


