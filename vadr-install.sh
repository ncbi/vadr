#!/bin/bash
# EPN, Wed Nov 20 06:06:22 2019
#
# vadr-install.sh
# A shell script for downloading and installing VADR and its dependencies.
# 
# usage: 
#
# To download and install:
# vadr-install.sh <"linux" or "macosx">
#
# To only download: 
# vadr-install.sh <"linux" or "macosx"> download
#
# for example:
# vadr-install.sh linux
# 
# The following line will make the script fail if any commands fail
set -e

VADRINSTALLDIR=$PWD

# versions
VERSION="1.1"
# bio-easel
BEVERSION="Bio-Easel-0.12"
# blast+
BVERSION="2.10.0"
# infernal
IVERSION="1.1.3"
# hmmer
HVERSION="3.3"
# dependency git tag
VVERSION="vadr-$VERSION"
# vadr models
MVERSION="1.1-1"

# set defaults
INPUTSYSTEM="?"

########################
# Validate correct usage
########################
# make sure correct number of cmdline arguments were used, exit if not
if [ "$#" -ne 1 ]; then
   echo "Usage: $0 <\"linux\" or \"macosx\">"
   exit 1
fi

# make sure 1st argument is either "linux" or "macosx"
if [ $1 == "linux" ]; then
    INPUTSYSTEM="linux";
fi
if [ $1 == "macosx" ]; then
    INPUTSYSTEM="macosx";
fi
if [ $INPUTSYSTEM == "?" ]; then 
   echo "Usage: $0 <\"linux\" or \"macosx\">"
   exit 1
fi
########################################################
echo "------------------------------------------------"
echo "DOWNLOADING AND BUILDING VADR $VERSION"
echo "------------------------------------------------"
echo ""
echo "************************************************************"
echo "IMPORTANT: BEFORE YOU WILL BE ABLE TO RUN VADR SCRIPTS,"
echo "YOU NEED TO FOLLOW THE INSTRUCTIONS OUTPUT AT THE END"
echo "OF THIS SCRIPT TO UPDATE YOUR ENVIRONMENT VARIABLES."
echo "************************************************************"

echo ""
echo "Determining current directory ... "
echo "Set VADRINSTALLDIR as current directory ($VADRINSTALLDIR)."

###########################################
# Download section
###########################################
echo "------------------------------------------------"
# vadr
echo "Downloading vadr ... "
curl -k -L -o $VVERSION.zip https://github.com/nawrockie/vadr/archive/$VVERSION.zip; unzip $VVERSION.zip; mv vadr-$VVERSION vadr; rm $VVERSION.zip
# for a test build of a release, comment out above curl and uncomment block below
# ----------------------------------------------------------------------------
#git clone https://github.com/nawrockie/vadr.git vadr
#cd vadr
#git checkout release-$VERSION
#rm -rf .git
#cd ..
# ----------------------------------------------------------------------------
 
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

# download infernal binary distribution
# - to download source distribution and build, see 
#   'infernal block 2' below.

# ----- infernal block 1 start  -----
if [ $INPUTSYSTEM == "linux" ]; then
    echo "Downloading Infernal version $IVERSION for Linux"
    curl -k -L -o infernal.tar.gz http://eddylab.org/infernal/infernal-$IVERSION-linux-intel-gcc.tar.gz
else
    echo "Downloading Infernal version $IVERSION for Mac/OSX"
    curl -k -L -o infernal.tar.gz http://eddylab.org/infernal/infernal-$IVERSION-macosx-intel.tar.gz
fi
tar xfz infernal.tar.gz
rm infernal.tar.gz
if [ $INPUTSYSTEM == "linux" ]; then
    mv infernal-$IVERSION-linux-intel-gcc infernal
else
    mv infernal-$IVERSION-macosx-intel infernal
fi
# ----- infernal block 1 end -----

# if you'd rather download the source distro and build it yourself
# (maybe because the binaries aren't working for you for some reason)
# comment out 'infernal block 1' above and 
# uncomment 'infernal block 2' below
# ----- infernal block 2 start  -----
#echo "Downloading Infernal version $IVERSION src distribution"
#curl -k -L -o infernal.tar.gz http://eddylab.org/infernal/infernal-$IVERSION.tar.gz
#tar xfz infernal.tar.gz
#rm infernal.tar.gz
#echo "Building Infernal ... "
#mv infernal-$IVERSION infernal
#cd infernal
#mkdir binaries
#sh ./configure --bindir=$PWD/binaries --prefix=$PWD
#make
#make install
#(cd easel/miniapps; make install)
#cd ..
#echo "Finished building Infernal "
# ----- infernal block 2 end -----
echo "------------------------------------------------"

# download hmmer source distribution
# (precompiled binaries are no longer provided as of v3.3)
echo "Downloading HMMER version $HVERSION src distribution"
curl -k -L -o hmmer.tar.gz http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz
tar xfz hmmer.tar.gz
rm hmmer.tar.gz
echo "Building HMMER ... "
mv hmmer-$HVERSION hmmer
cd hmmer
mkdir binaries
sh ./configure --bindir=$PWD/binaries --prefix=$PWD
make
make install
cd ..
echo "Finished building HMMER "
echo "------------------------------------------------"

# download blast binaries
if [ $INPUTSYSTEM == "linux" ]; then
echo "Downloading BLAST version $BVERSION for Linux"
curl -k -L -o blast.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BVERSION/ncbi-blast-$BVERSION+-x64-linux.tar.gz
else 
echo "Downloading BLAST version $BVERSION for Mac/OSX"
curl -k -L -o blast.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BVERSION/ncbi-blast-$BVERSION+-x64-macosx.tar.gz
fi
tar xfz blast.tar.gz
rm blast.tar.gz
mv ncbi-blast-$BVERSION+ ncbi-blast
echo "------------------------------------------------"

# download vadr-models 
echo "Downloading latest VADR models ... "
curl -k -L -o vadr-models.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/CURRENT/vadr-models-$MVERSION.tar.gz
# for a test build of a release, or of the develop branch, you may want different models,
# such as those in the develop/ dir, in that case comment out above curl and uncomment
# and possibly modify the one below
# ----------------------------------------------------------------------------
#curl -k -L -o vadr-models.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/develop/vadr-models-$MVERSION.tar.gz
# ----------------------------------------------------------------------------
tar xfz vadr-models.tar.gz
rm vadr-models.tar.gz
echo "------------------------------------------------"

###########################################
# Build section
###########################################
if [ ! -d Bio-Easel ]; then
   echo "ERROR: Bio-Easel dir does not exist"
   exit 1
fi
# Build Bio-Easel:
echo "------------------------------------------------"
echo "Building Bio-Easel ... "
cd Bio-Easel
perl Makefile.PL
make
make test
cd ..
echo "Finished building Bio-Easel."
echo "------------------------------------------------"

###############################################
# Message about setting environment variables
###############################################
echo ""
echo ""
echo "********************************************************"
echo "The final step is to update your environment variables."
echo "(See https://github.com/nawrockie/vadr/blob/$VERSION/documentation/install.md for more information.)"
echo ""
echo "If you are using the bash or zsh shell (zsh is default in MacOS/X as"
echo "of v10.15 (Catalina)), add the following lines to the end of your"
echo "'.bashrc' or '.zshrc' file in your home directory:"
echo ""
echo "export VADRINSTALLDIR=\"$VADRINSTALLDIR\""
echo "export VADRSCRIPTSDIR=\"\$VADRINSTALLDIR/vadr\""
echo "export VADRMODELDIR=\"\$VADRINSTALLDIR/vadr-models\""
echo "export VADRINFERNALDIR=\"\$VADRINSTALLDIR/infernal/binaries\""
echo "export VADREASELDIR=\"\$VADRINSTALLDIR/infernal/binaries\""
echo "export VADRHMMERDIR=\"\$VADRINSTALLDIR/hmmer/binaries\""
echo "export VADRBIOEASELDIR=\"\$VADRINSTALLDIR/Bio-Easel\""
echo "export VADRSEQUIPDIR=\"\$VADRINSTALLDIR/sequip\""
echo "export VADRBLASTDIR=\"\$VADRINSTALLDIR/ncbi-blast/bin\""
echo "export PERL5LIB=\"\$VADRSCRIPTSDIR\":\"\$VADRSEQUIPDIR\":\"\$VADRBIOEASELDIR/blib/lib\":\"\$VADRBIOEASELDIR/blib/arch\":\"\$PERL5LIB\""
echo "export PATH=\"\$VADRSCRIPTSDIR\":\"\$PATH\""
echo ""
echo "After adding the export lines to your .bashrc or .zshrc file, source that file"
echo "to update your current environment with the command:"
echo ""
echo "source ~/.bashrc"
echo ""
echo "OR"
echo ""
echo "source ~/.zshrc"
echo ""
echo "---"
echo "If you are using the C shell, add the following"
echo "lines to the end of your '.cshrc' file in your home"
echo "directory:"
echo ""
echo "setenv VADRINSTALLDIR \"$VADRINSTALLDIR\""
echo "setenv VADRSCRIPTSDIR \"\$VADRINSTALLDIR/vadr\""
echo "setenv VADRMODELDIR \"\$VADRINSTALLDIR/vadr-models\""
echo "setenv VADRINFERNALDIR \"\$VADRINSTALLDIR/infernal/binaries\""
echo "setenv VADRHMMERDIR \"\$VADRHMMERDIR/hmmer/binaries\""
echo "setenv VADREASELDIR \"\$VADRINSTALLDIR/infernal/binaries\""
echo "setenv VADRBIOEASELDIR \"\$VADRINSTALLDIR/Bio-Easel\""
echo "setenv VADRSEQUIPDIR \"\$VADRINSTALLDIR/sequip\""
echo "setenv VADRBLASTDIR \"\$VADRINSTALLDIR/ncbi-blast/bin\""
echo "setenv PERL5LIB \"\$VADRSCRIPTSDIR\":\"\$VADRSEQUIPDIR\":\"\$VADRBIOEASELDIR/blib/lib\":\"\$VADRBIOEASELDIR/blib/arch\":\"\$PERL5LIB\""
echo "setenv PATH \"\$VADRSCRIPTSDIR\":\"\$PATH\""
echo ""
echo "After adding the setenv lines to your .cshrc file, source that file"
echo "to update your current environment with the command:"
echo ""
echo "source ~/.cshrc"
echo ""
echo "(To determine which shell you use, type: 'echo \$SHELL')"
echo ""
echo ""
echo "********************************************************"
echo ""
