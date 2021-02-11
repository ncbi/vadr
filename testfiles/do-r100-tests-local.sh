#!/bin/bash

RETVAL=0;

# noro r100
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r100.local.testin n100-local
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# dengue r100
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r100.local.testin d100-local
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-r100-tests-local.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed do-r100-tests-local.sh]"
   exit 1
fi


