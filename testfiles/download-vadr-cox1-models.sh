#!/bin/bash
# The previous line forces this script to be run with bash, regardless of user's 
# shell.
#
# EPN, Tue Oct 15 14:38:57 2019
#
# checkout-source.sh 
# A shell script for downloading VADR COX1 models.
# 
# usage: 
# download-vadr-cox1-models.sh <vadr cox1 model version (e.g. 0.971.1)>
# 
# for example:
# download-vadr-cox1-models.sh <0.971.1>
#
# The following line will make the script fail if any commands fail
set -e

# make sure correct number of cmdline arguments were used, exit if not
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <vadr cox1 model version (e.g. 0.971.1)>"
  exit 1
fi

MVERSION=$1

# download vadr-cox1-models 
echo "Downloading VADR COX1 models ... "
curl -k -L -o vadr-cox1-models.tar.gz http://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/cox1/CURRENT/vadr-cox1-models-$MVERSION-dir.tar.gz
gunzip vadr-cox1-models.tar.gz
tar xf vadr-cox1-models.tar
rm vadr-cox1-models.tar
echo "------------------------------------------------"


