#!/bin/bash
set -e

sh $VADRSCRIPTSDIR/testfiles/do-install-tests-local.sh
sh $VADRSCRIPTSDIR/testfiles/do-fs-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-replace-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-seed-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-outaln-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-nindel-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-hmmer-tests.sh

# to test if -p will work (requires qsub and qsub flags similar to NCBI internal set-up)
# sh $VADRSCRIPTSDIR/testfiles/do-install-tests-parallel.sh
