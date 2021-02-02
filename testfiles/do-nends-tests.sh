#!/bin/bash

RETVAL=0;

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a.nends.testin vt-entoy100a-nends
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-rev.nends.testin vt-entoy100a-rev-nends
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.nends.testin vt-noro-nends
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-nends-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-nends-tests.sh]"
   exit 1
fi
