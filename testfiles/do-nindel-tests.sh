#!/bin/bash

RETVAL=0;

# entoy100a-nindel
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-nindel.testin entoy100a-nindel
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
