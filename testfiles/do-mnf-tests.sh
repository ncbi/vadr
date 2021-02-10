#!/bin/bash

RETVAL=0;

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-mnf-fs.testin vt-entoy100a-mnf-fs
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-mnf-lowsim.testin vt-entoy100a-mnf-lowsim
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [$0]"
   exit 0
else 
   echo "FAIL: at least one test failed [$0]"
   exit 1
fi
