#!/bin/bash
# The previous line forces this script to be run with bash, regardless of user's 
# shell.
#
# EPN, Mon May  6 09:21:35 2019
#
# A shell script for installing VADR its dependencies
# for viral sequence classification and annotation.
#
VADRINSTALLDIR=$PWD
VERSION="0.97"
VVERSION="vadr-$VERSION"

# The following line will make the script fail if any commands fail
set -e
#
echo "------------------------------------------------"
echo "INSTALLING VADR $VERSION"
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

echo "------------------------------------------------"
# Clone what we need from GitHub (these are all public)

# vadr
# ONLY DIFFERENCE BETWEEN THIS SCRIPT AND install.sh IS THAT WE USE GIT CLONE TO GET vadr INSTEAD OF INSTALLING A SET VERSION
#echo "Installing vadr ... "
#curl -k -L -o vadr-$VERSION.zip https://github.com/nawrockie/vadr/archive/$VERSION.zip; unzip vadr-$VERSION.zip; mv vadr-$VERSION vadr; rm vadr-$VERSION.zip
echo "Installing 'git pull'-updatable copy of vadr ... "
git clone https://github.com/nawrockie/vadr.git

# epn-options, epn-ofile, and Bio-Easel
for m in epn-options epn-ofile Bio-Easel; do 
    echo "Installing $m ... "
    curl -k -L -o $m-$VVERSION.zip https://github.com/nawrockie/$m/archive/$VVERSION.zip; unzip $m-$VVERSION.zip; mv $m-$VVERSION $m; rm $m-$VVERSION.zip
    # we can't leave these directories with $VVERSION because vadr scripts expect non-versioned names
done
echo "------------------------------------------------"

# download vadr-models 
echo "Downloading latest VADR models ... "
curl -k -L -o vadr-models.tar.gz http://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/CURRENT/vadr-models-0.94.2.tar.gz
gunzip vadr-models.tar.gz
tar xf vadr-models.tar
rm vadr-models.tar
echo "------------------------------------------------"

########## BUILD BIO-EASEL ###############
# Build Bio-Easel:
echo "Building Bio-Easel ... "
cd $VADRINSTALLDIR/Bio-Easel
# clone Easel subdir
mkdir src
(cd src; git clone https://github.com/EddyRivasLab/easel.git easel)
(cd src/easel; git checkout tags/Bio-Easel-0.08; rm -rf .git)
perl Makefile.PL
make
make test
echo "Finished building Bio-Easel."
echo "------------------------------------------------"
cd $VADRINSTALLDIR
########## BUILD BIO-EASEL ###############

##########BEGINNING OF LINES TO COMMENT OUT TO SKIP INFERNAL INSTALLATION##########################
# Install Infernal (develop branch, specific commit)
# UPDATE THIS TO USE curl TO DOWNLOAD 1.1.3 WHEN IT IS AVAILABLE
cd $VADRINSTALLDIR
echo "Installing Infernal (develop branch) ... "
git clone https://github.com/EddyRivasLab/infernal.git infernal-dev
cd infernal-dev
git checkout 7d93882
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
cd $VADRINSTALLDIR
echo "Finished installing Infernal "
echo "------------------------------------------------"
##########END OF LINES TO COMMENT OUT TO SKIP INFERNAL INSTALLATION##########################
# 
################
# Output the final message:
echo ""
echo ""
echo "********************************************************"
echo "The final step is to update your environment variables."
echo "(See vadr/README.txt for more information.)"
echo ""
echo "If you are using the bash shell, add the following"
echo "lines to the end of your '.bashrc' file in your home"
echo "directory:"
echo ""
echo "export VADRINSTALLDIR=\"$VADRINSTALLDIR\""
echo "export VADRSCRIPTSDIR=\"\$VADRINSTALLDIR/vadr\""
echo "export VADRMODELDIR=\"$VADRINSTALLDIR/vadr-models\""
echo "export VADRINFERNALDIR=\"\$VADRINSTALLDIR/infernal-dev/src\""
echo "export VADREASELDIR=\"\$VADRINSTALLDIR/infernal-dev/easel/miniapps\""
echo "export VADRBIOEASELDIR=\"\$VADRINSTALLDIR/Bio-Easel\""
echo "export VADRBLASTDIR=\"/usr/bin\""
echo "export PERL5LIB=\"\$VADRSCRIPTSDIR\":\"\$VADRINSTALLDIR/epn-options\":\"\$VADRINSTALLDIR/epn-ofile\":\"\$VADRBIOEASELDIR/blib/lib\":\"\$VADRBIOEASELDIR/blib/arch\":\"\$PERL5LIB\""
echo "export PATH=\"\$VADRSCRIPTSDIR\":\"\$PATH\""
echo ""
echo "After adding the export lines to your .bashrc file, source that file"
echo "to update your current environment with the command:"
echo ""
echo "source ~/.bashrc"
echo ""
echo "---"
echo "If you are using the C shell, add the following"
echo "lines to the end of your '.cshrc' file in your home"
echo "directory:"
echo ""
echo "setenv VADRINSTALLDIR \"$VADRINSTALLDIR\""
echo "setenv VADRSCRIPTSDIR \"\$VADRINSTALLDIR/vadr\""
echo "setenv VADRMODELDIR \"$VADRINSTALLDIR/vadr-models\""
echo "setenv VADRINFERNALDIR \"\$VADRINSTALLDIR/infernal-dev/src\""
echo "setenv VADREASELDIR \"\$VADRINSTALLDIR/infernal-dev/easel/miniapps\""
echo "setenv VADRBIOEASELDIR \"\$VADRINSTALLDIR/Bio-Easel\""
echo "setenv VADRBLASTDIR \"/usr/bin\""
echo "setenv PERL5LIB \"\$VADRSCRIPTSDIR\":\"\$VADRINSTALLDIR/epn-options\":\"\$VADRINSTALLDIR/epn-ofile\":\"\$VADRBIOEASELDIR/blib/lib\":\"\$VADRBIOEASELDIR/blib/arch\":\"\$PERL5LIB\""
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
