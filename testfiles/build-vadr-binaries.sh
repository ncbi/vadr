#!/bin/bash
# The previous line forces this script to be run with bash, regardless of user's 
# shell.
#
# EPN, Tue Oct 15 14:44:22 2019
#
# build-vadr-binaries.sh
# A shell script for building VADR dependencies. 
# To be run after 'checkout-source.sh' in the same directory that
# 'checkout-source.sh' was run in.
# 
# usage: 
# build-vadr-binaries.sh
# 
# for example:
# build-vadr-binaries.sh
#
# The following line will make the script fail if any commands fail
set -e
#
# make sure correct number of cmdline arguments were used, exit if not
if [ "$#" -ne 0 ]; then
  echo "Usage: $0"
  exit 1
fi

# Build Bio-Easel:
echo "Building Bio-Easel ... "
cd Bio-Easel
perl Makefile.PL
make
make test
cd ..
echo "Finished building Bio-Easel."
echo "------------------------------------------------"

# Build Infernal:
echo "Building Infernal ... "
cd infernal-dev
sh ./configure 
make
cd ..
echo "Finished building Infernal "
echo "------------------------------------------------"

