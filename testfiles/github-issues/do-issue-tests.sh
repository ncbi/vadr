#!/bin/bash

RETVAL=0;

for i in \
    iss3 \
    iss4 \
    iss5 \
    iss6 \
    iss13 \
    iss21-delcds \
    iss30-toomanyNs \
    iss33-longdesc \
    iss37-roverlap \
    iss47-productparantheses \
    iss58-cdsstopp \
    iss61-blaststrand \
    iss68-blastmaxde \
    iss69-shortfs \
    iss70-cdsstopn3p \
    ; do
    $VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/github-issues/$i/$i.testin vt-$i
    if [ $? != 0 ]; then
        RETVAL=1;
    fi   
done

if [ $RETVAL = 0 ]; then
   echo "Success: all tests passed [do-issue-tests.sh]"
   exit 0
else 
   echo "FAIL: at least one test failed [do-issue-tests.sh]"
   exit 1
fi
