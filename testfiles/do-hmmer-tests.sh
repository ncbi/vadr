#!/bin/bash

# noro r10 hmmer
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/noro.r10.hmmer.testin vt-n10-hmmer

# dengue r5 hmmer
$VADRSCRIPTSDIR/v-test.pl -f --rmout $VADRSCRIPTSDIR/testfiles/dengue.r5.hmmer.testin vt-d5-hmmer
