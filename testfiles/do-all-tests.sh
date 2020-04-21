#!/bin/bash
set -e

sh $VADRSCRIPTSDIR/testfiles/do-install-tests-local.sh
sh $VADRSCRIPTSDIR/testfiles/do-fs-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-replace-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-seed-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-outaln-tests.sh

sh $VADRSCRIPTSDIR/t/do-prove-all-tests.sh teamcity

# If you want to test -p option for parallelization, uncomment the 
# next 'sh' line.
# Note: this test requires qsub is in your path and qsub options are
# configured similarly to ncbi cluster, email eric.nawrocki@nih.gov
# for information on how to configure for different clusters

# sh $VADRSCRIPTSDIR/testfiles/do-install-tests-parallel.sh
