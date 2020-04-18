#!/bin/bash
sh do-install-tests-local.sh
sh do-fs-tests.sh
sh do-replace-tests.sh
sh do-seed-tests.sh 
# sh do-install-tests-parallel.sh
# sh do-r100-tests-local.sh
# sh do-r1000-tests-parallel.sh
