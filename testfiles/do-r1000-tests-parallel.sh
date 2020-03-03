#!/bin/bash

# noro r1000
$VADRSCRIPTSDIR/v-test.pl -f $VADRSCRIPTSDIR/testfiles/noro.r1000.parallel.testin n1000-parallel

# dengue r1000
$VADRSCRIPTSDIR/v-test.pl -f $VADRSCRIPTSDIR/testfiles/dengue.r1000.parallel.testin d1000-parallel


