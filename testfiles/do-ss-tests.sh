#!/bin/bash

RETVAL=0;

$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy200b-ss.testin entoy200b-ss
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   
