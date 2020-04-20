#!/bin/bash
set -e

sh $VADRSCRIPTSDIR/testfiles/do-install-tests-local.sh
sh $VADRSCRIPTSDIR/testfiles/do-fs-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-replace-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-seed-tests.sh
sh $VADRSCRIPTSDIR/testfiles/do-outaln-tests.sh
# sh $VADRSCRIPTSDIR/testfiles/do-install-tests-parallel.sh
# sh $VADRSCRIPTSDIR/testfiles/do-r100-tests-local.sh
# sh $VADRSCRIPTSDIR/testfiles/do-r1000-tests-parallel.sh
