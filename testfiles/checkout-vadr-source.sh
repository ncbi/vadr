#!/bin/bash
# The previous line forces this script to be run with bash, regardless of user's 
# shell.
#
# EPN, Tue Oct 15 14:38:57 2019
#
# checkout-source.sh 
# A shell script for downloading VADR and its dependencies.
# 
# usage: 
# checkout-source.sh <vadr version (e.g. 0.991)> <vadr model version (e.g. 0.991.1)> <blast version (e.g. 2.9.0)>"
# 
# for example:
# checkout-source.sh 0.991 0.991.1 2.9.0+
#
# The following line will make the script fail if any commands fail
set -e

# make sure correct number of cmdline arguments were used, exit if not
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <vadr version (e.g. 0.991)> <vadr model version (e.g. 0.991.1)> <blast version (e.g. 2.9.0)>"
  exit 1
fi

VERSION=$1
VVERSION="vadr-$VERSION"
MVERSION=$2
BVERSION=$3
# TEMPORARY 
BEVERSION=Bio-Easel-0.08

# vadr
echo "Downloading vadr ... "
curl -k -L -o vadr-$VERSION.zip https://github.com/nawrockie/vadr/archive/$VERSION.zip; unzip vadr-$VERSION.zip; mv vadr-$VERSION vadr; rm vadr-$VERSION.zip

# sequip and Bio-Easel
for m in sequip Bio-Easel; do 
    echo "Downloading $m ... "
    curl -k -L -o $m-$VVERSION.zip https://github.com/nawrockie/$m/archive/$VVERSION.zip; unzip $m-$VVERSION.zip; mv $m-$VVERSION $m; rm $m-$VVERSION.zip
done
cd Bio-Easel
mkdir src
(cd src; curl -k -L -o easel-$BEVERSION.zip https://github.com/EddyRivasLab/easel/archive/$BEVERSION.zip; unzip easel-$BEVERSION.zip; mv easel-$BEVERSION easel; rm easel-$BEVERSION.zip; cd easel; autoconf)
cd ..
echo "------------------------------------------------"

# download Infernal (TEMPORARY: develop branch, specific commit, will be v1.1.3 once that is released)
# UPDATE THIS TO USE curl TO DOWNLOAD 1.1.3 WHEN IT IS AVAILABLE
echo "Downloading Infernal (develop branch) ... "
git clone https://github.com/EddyRivasLab/infernal.git infernal-dev
cd infernal-dev
git checkout 7d93882
rm -rf .git
git clone https://github.com/EddyRivasLab/hmmer
(cd hmmer; git checkout 498ec7c; rm -rf .git)
git clone https://github.com/EddyRivasLab/easel
(cd easel; git checkout 5288a95; rm -rf git)
autoconf
cd ..
echo "------------------------------------------------"

# download blast binaries
# For Linux: 
echo "Downloading BLAST version $BVERSION for Linux"
curl -k -L -o blast.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BVERSION/ncbi-blast-$BVERSION+-x64-linux.tar.gz
gunzip blast.tar.gz
tar xf blast.tar
mv ncbi-blast-2.9.0+ ncbi-blast
rm blast.tar
echo "------------------------------------------------"

# download vadr-models 
echo "Downloading latest VADR models ... "
curl -k -L -o vadr-models.tar.gz http://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/CURRENT/vadr-models-$MVERSION.tar.gz
gunzip vadr-models.tar.gz
tar xf vadr-models.tar
rm vadr-models.tar
echo "------------------------------------------------"


