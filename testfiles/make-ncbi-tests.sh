#!/bin/bash

# add the 'no option tests' first
for a in \
dengue.r5.local \
noro.r10.local \
; do 
    perl $VADR_BB_TEST_DIR/regression_tests/testin2ncbi.pl --noopts $VADRSCRIPTSDIR/testfiles/$a.testin
done

# add the testsuite tests, keeping options 
for a in \
dengue.r5.hmmer \
dengue.r5.rpn \
dengue.r5.seed \
entoy100a-fs \
entoy100a-nindel \
entoy100a-rev-fs \
noro.fs.multisgm \
noro.fs \
noro.r10.hmmer \
noro.r10.rpn \
noro.r10.seed \
noro.r3.outaln \
; do 
    perl $VADR_BB_TEST_DIR/regression_tests/testin2ncbi.pl $VADRSCRIPTSDIR/testfiles/$a.testin
done


