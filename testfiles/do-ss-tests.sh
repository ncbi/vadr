#!/bin/bash

RETVAL=0;

# entoy200b-ss
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy200b-ss.testin entoy200b-ss
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# entoy200b-ss-partials
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy200b-ss-partials.testin entoy200b-ss-partials
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# entoy200b-rev-ss
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy200b-rev-ss.testin entoy200b-rev-ss
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# entoy200b-rev-ss-partials
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy200b-rev-ss-partials.testin entoy200b-rev-ss-partials
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-ss-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-ss-tests.sh]"
   exit 1
fi
