#!/bin/bash

RETVAL=0;

# noro r10 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.glsearch.testin vt-n10-glsearch
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# dengue r5 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.glsearch.testin vt-d5-glsearch
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

if [ "$RETVAL" -eq 0 ]; then
   echo "Success: all tests passed [do-glsearch-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-glsearch-tests.sh]"
   exit 1
fi
