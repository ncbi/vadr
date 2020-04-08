#!/bin/bash

# noro.r10
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-fs.testin entoy100a-fs

# noro.r10.1 --overhang 1
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/entoy100a-rev-fs.testin entoy100a-rev-fs

