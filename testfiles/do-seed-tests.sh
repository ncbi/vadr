#!/bin/bash

# noro r10 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.seed.testin vt-n10-seed

# dengue r5 seed
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.seed.testin vt-d5-seed
