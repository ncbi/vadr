#!/bin/bash
# The previous line forces this script to be run with bash, regardless of user's 
# shell.
#
# EPN, Thu Jan 24 13:22:39 2019
#
# A shell script for installing dnaorg_scripts its dependencies
# for viral sequence classification and annotation.
#
DNAORGINSTALLDIR=$PWD
VERSION="0.44"
DVERSION="dnaorg_scripts-$VERSION"

# The following line will make the script fail if any commands fail
set -e
#
echo "------------------------------------------------"
echo "INSTALLING DNAORG_SCRIPTS $VERSION"
echo "------------------------------------------------"
echo ""
echo "************************************************************"
echo "IMPORTANT: BEFORE YOU WILL BE ABLE TO RUN DNAORG"
echo "SCRIPTS, YOU NEED TO FOLLOW THE INSTRUCTIONS OUTPUT AT"
echo "THE END OF THIS SCRIPT TO UPDATE YOUR ENVIRONMENT VARIABLES."
echo "************************************************************"
echo ""
echo "Determining current directory ... "
echo "Set DNAORGINSTALLDIR as current directory ($DNAORGINSTALLDIR)."

echo "------------------------------------------------"
# Clone what we need from GitHub (these are all public)

# dnaorg_scripts
echo "Installing dnaorg_scripts ... "
curl -k -L -o dnaorg_scripts-$VERSION.zip https://github.com/nawrockie/dnaorg_scripts/archive/$VERSION.zip; unzip dnaorg_scripts-$VERSION.zip; mv dnaorg_scripts-$VERSION dnaorg_scripts; rm dnaorg_scripts-$VERSION.zip

# esl-epn-translate, esl-fetch-cds, epn-options, Bio-Easel
for m in esl-epn-translate esl-fetch-cds epn-options Bio-Easel; do 
    echo "Installing $m ... "
    curl -k -L -o $m-$DVERSION.zip https://github.com/nawrockie/$m/archive/$DVERSION.zip; unzip $m-$DVERSION.zip; mv $m-$DVERSION $m; rm $m-$DVERSION.zip
    # we can't leave these directories with $DVERSION because dnaorg scripts expect non-versioned names
done
echo "------------------------------------------------"

# dnaorg-build-directories, use git here so user can update it with 'git pull'
echo "Cloning updatable github repo with virus models ... "
git clone https://github.com/nawrockie/dnaorg-build-directories.git
echo "------------------------------------------------"

########## BUILD BIO-EASEL ###############
# Build Bio-Easel:
echo "Building Bio-Easel ... "
cd $DNAORGINSTALLDIR/Bio-Easel
# clone Easel subdir
mkdir src
(cd src; git clone https://github.com/EddyRivasLab/easel.git easel)
(cd src/easel; git checkout tags/Bio-Easel-0.06; rm -rf .git)
perl Makefile.PL
make
make test
echo "Finished building Bio-Easel."
echo "------------------------------------------------"
cd $DNAORGINSTALLDIR
########## BUILD BIO-EASEL ###############

##########BEGINNING OF LINES TO COMMENT OUT TO SKIP INFERNAL INSTALLATION##########################
# Install Infernal (develop branch, specific commit)
# UPDATE THIS TO USE curl TO DOWNLOAD 1.1.3 WHEN IT IS AVAILABLE
cd $DNAORGINSTALLDIR
echo "Installing Infernal (develop branch) ... "
git clone https://github.com/EddyRivasLab/infernal.git infernal-dev
cd infernal-dev
git checkout b748e2c
rm -rf .git
git clone https://github.com/EddyRivasLab/hmmer
(cd hmmer; git checkout 498ec7c; rm -rf .git)
git clone https://github.com/EddyRivasLab/easel
(cd easel; git checkout 5288a95; rm -rf git)
autoconf
sh ./configure 
make
# uncomment to install system wide
#make install
cd $DNAORGINSTALLDIR
echo "Finished installing Infernal "
echo "------------------------------------------------"
##########END OF LINES TO COMMENT OUT TO SKIP INFERNAL INSTALLATION##########################
# 
##########BEGINNING OF LINES TO COMMENT OUT TO SKIP HMMER INSTALLATION##########################
cd $DNAORGINSTALLDIR
echo "Installing HMMER 3.1b2 ... "
curl -k -L -o hmmer-3.1b2.tar.gz http://eddylab.org/software/hmmer/hmmer-3.1b2.tar.gz
tar xfz hmmer-3.1b2.tar.gz
cd hmmer-3.1b2
sh ./configure 
make
# uncomment to install system wide
#make install
cd $DNAORGINSTALLDIR
echo "Finished installing HMMER 3.1b2."
echo "------------------------------------------------"
##########END OF LINES TO COMMENT OUT TO SKIP HMMER INSTALLATION##########################

################
# Output the final message:
echo ""
echo ""
echo "********************************************************"
echo "The final step is to update your environment variables."
echo "(See dnaorg_scripts/README.txt for more information.)"
echo ""
echo "If you are using the bash shell, add the following"
echo "lines to the '.bashrc' file in your home directory:"
echo ""
echo "export DNAORGDIR=$DNAORGINSTALLDIR"
echo "export DNAORGBUILDDIR=\"\$DNAORGDIR\"/dnaorg-build-directories"
echo "export PERL5LIB=\"\$DNAORGDIR\"/dnaorg_scripts:\"\$DNAORGDIR\"/epn-options:\"\$DNAORGDIR\"/Bio-Easel/blib/lib:\"\$DNAORGDIR\"/Bio-Easel/blib/arch:\"\$PERL5LIB\""
echo "export PATH=\"\$DNAORGDIR\"/dnaorg_scripts:\"\$PATH\""
echo ""
echo "After adding the export lines to your .bashrc file, source that file"
echo "to update your current environment with the command:"
echo ""
echo "source ~/.bashrc"
echo ""
echo "---"
echo "If you are using the C shell, add the following"
echo "lines to the '.cshrc' file in your home directory:"
echo ""
echo "setenv DNAORGDIR $DNAORGINSTALLDIR"
echo "setenv DNAORGBUILDDIR \"\$DNAORGDIR\"/dnaorg-build-directories"
echo "setenv PERL5LIB \"\$DNAORGDIR\"/dnaorg_scripts:\"\$DNAORGDIR\"/epn-options:\"\$DNAORGDIR\"/Bio-Easel/blib/lib:\"\$DNAORGDIR\"/Bio-Easel/blib/arch:\"\$PERL5LIB\""
echo "setenv PATH \"\$DNAORGDIR\"/dnaorg_scripts:\"\$PATH\""
echo ""
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
