#!/bin/bash

RETVAL=0;

# entoy100a-fs{1,2,3}
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-fs.testin entoy100a-fs
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# entoy100a-rev-fs{1,2,3}
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-rev-fs.testin entoy100a-rev-fs
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# noro.fs
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.fs.testin noro.fs
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# noro.fs.multisgm
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.fs.multisgm.testin noro.fs.multisgm
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-fs-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-fs-tests.sh]"
   exit 1
fi
