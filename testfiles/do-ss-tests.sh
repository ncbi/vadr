#!/bin/bash

RETVAL=0;

# entoy200b-ss
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy200b-ss.testin entoy200b-ss
if [ "$?" -ne 0 ]; then
   RETVAL=1;
fi   

# entoy200b-rev-ss
#$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy200b-rev-ss.testin entoy200b-rev-ss
#if [ "$?" -ne 0 ]; then
#   RETVAL=1;
#fi   
