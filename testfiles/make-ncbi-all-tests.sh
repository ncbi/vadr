#!/bin/bash

if [ -e ncbi-all-tests.sh ]; then 
    rm ncbi-all-tests.sh
fi

for a in \
dengue.r5.hmmer \
dengue.r5.local \
dengue.r5.parallel \
dengue.r5.rpn \
dengue.r5.seed \
entoy100a-fs \
entoy100a-nindel \
entoy100a-rev-fs \
noro.fs.multisgm \
noro.fs \
noro.r10.hmmer \
noro.r10.local \
noro.r10.parallel \
noro.r10.rpn \
noro.r10.seed \
noro.r3.outaln \
; do 
    perl ./testin2ncbi.pl $a.testin >> ncbi-all-tests.sh
done


