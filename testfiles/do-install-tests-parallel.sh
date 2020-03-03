#!/bin/bash

# noro r10
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.parallel.testin n10-parallel

# dengue r5
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.parallel.testin d5-parallel
