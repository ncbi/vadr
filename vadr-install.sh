#!/bin/bash
# EPN, Thu Oct 24 14:48:19 2019
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
# To only build (only works if already run in download mode in same dir): 
# vadr-install.sh <"linux" or "macosx"> build
#
# for example:
# vadr-install.sh linux
# 
# OR
# 
# vadr-install.sh macosx download
# followed by 
# vadr-install.sh macosx build
#
# The following line will make the script fail if any commands fail
set -e

VADRINSTALLDIR=$PWD

# versions
VERSION="0.991"
BEVERSION="Bio-Easel-0.08"
BVERSION="2.9.0"
VVERSION="vadr-$VERSION"
MVERSION="0.991.1"

# set defaults
INPUTSYSTEM="?"
MODE="both"

########################################################
# Validate correct usage:
########################################################
# make sure correct number of cmdline arguments were used, exit if not
# and if 2 are used, that second one is "download" or "build"
if [ "$#" -ne 1 ]; then
    if [ "$#" -ne 2 ]; then
        echo "Usage:"
        echo ""
        echo "To download source and build:"
        echo ""
        echo "$0 <\"linux\" or \"macosx\">"
        echo ""
        echo "To only download source *or* build:"
        echo ""
        echo "$0 <\"linux\" or \"macosx\"> <\"download\" or \"build\">"
        exit 1
    else 
        if [ $2 == "download" ]; then
            MODE="download";
        fi
        if [ $2 == "build" ]; then
            MODE="build";
        fi
        if [ $MODE == "both" ]; then 
            echo "Usage:"
            echo "To download source and build:"
            echo "$0 <\"linux\" or \"macosx\">"
            echo "OR"
            echo "To only download source or build:"
            echo "$0 <\"linux\" or \"macosx\"> <\"download\" or \"build\">"
            exit 1
        fi
    fi
fi

# make sure 1st argument is either "linux" or "macosx"
if [ $1 == "linux" ]; then
    INPUTSYSTEM="linux";
fi
if [ $1 == "macosx" ]; then
    INPUTSYSTEM="macosx";
fi
if [ $INPUTSYSTEM == "?" ]; then 
    echo "Usage:"
    echo "To download source and build:"
    echo "$0 <\"linux\" or \"macosx\">"
    echo "OR"
    echo "To only download source or build:"
    echo "$0 <\"linux\" or \"macosx\"> <\"download\" or \"build\">"
    exit 1
fi
########################################################

#
if [ $MODE == "both" ]; then
    echo "------------------------------------------------"
    echo "DOWNLOADING AND BUILDING VADR $VERSION"
    echo "------------------------------------------------"
fi
if [ $MODE == "download" ]; then
    echo "------------------------------------------------"
    echo "DOWNLOADING VADR $VERSION"
    echo "------------------------------------------------"
fi
if [ $MODE == "build" ]; then
    echo "------------------------------------------------"
    echo "BUILDING VADR $VERSION"
    echo "------------------------------------------------"
fi

echo ""
if [ $MODE != "download" ]; then
    echo "************************************************************"
    echo "IMPORTANT: BEFORE YOU WILL BE ABLE TO RUN VADR SCRIPTS,"
    echo "YOU NEED TO FOLLOW THE INSTRUCTIONS OUTPUT AT THE END"
    echo "OF THIS SCRIPT TO UPDATE YOUR ENVIRONMENT VARIABLES."
    echo "************************************************************"
fi
if [ $MODE == "download" ]; then
    echo ""
    echo "************************************************************"
    echo "IMPORTANT: BEFORE YOU WILL BE ABLE TO RUN VADR SCRIPTS,"
    echo "YOU NEED TO BUILD THE VADR EXECUTABLE FILES BY RERUNNING"
    echo "THIS SCRIPT FROM THIS SAME DIRECTORY LIKE THIS"
    echo "./vadr-install.sh <\"linux\" or \"macosx\"> build"
    echo "************************************************************"
fi

echo ""
echo "Determining current directory ... "
echo "Set VADRINSTALLDIR as current directory ($VADRINSTALLDIR)."

###########################################
# Download section
###########################################
if [ $MODE != "build" ]; then
    echo "------------------------------------------------"
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
    git clone https://github.com/EddyRivasLab/infernal.git infernal
    cd infernal
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
    if [ $INPUTSYSTEM == "linux" ]; then
        echo "Downloading BLAST version $BVERSION for Linux"
        curl -k -L -o blast.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BVERSION/ncbi-blast-$BVERSION+-x64-linux.tar.gz
    else 
        echo "Downloading BLAST version $BVERSION for Mac/OSX"
        curl -k -L -o blast.tar.gz ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BVERSION/ncbi-blast-$BVERSION+-x64-macosx.tar.gz
    fi
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
fi

###########################################
# Build section
###########################################
if [ $MODE != "download" ]; then
    if [ ! -d Bio-Easel ]; then
        echo "ERROR: Bio-Easel dir does not exist"
        if [ $MODE == "build" ]; then 
            echo "Did you run 'vadr-install.sh <\"linux\" or \"macosx\"> download' first?"
        fi
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

    if [ ! -d infernal ]; then
        echo "ERROR: infernal dir does not exist"
        if [ $MODE == "build" ]; then 
            echo "Did you run 'vadr-install.sh <\"linux\" or \"macosx\"> download' first?"
        fi
        exit 1
    fi
    # Build Infernal:
    echo "Building Infernal ... "
    cd infernal
    sh ./configure 
    make
    cd ..
    echo "Finished building Infernal "
    echo "------------------------------------------------"
fi

###############################################
# Message about setting environment variables
###############################################
echo ""
echo ""
echo "********************************************************"
echo "The final step is to update your environment variables."
echo "(See https://github.com/nawrockie/vadr/blob/0.991/documentation/install.md for more information.)"
echo ""
echo "If you are using the bash shell, add the following"
echo "lines to the end of your '.bashrc' file in your home"
echo "directory:"
echo ""
echo "export VADRINSTALLDIR=\"$VADRINSTALLDIR\""
echo "export VADRSCRIPTSDIR=\"\$VADRINSTALLDIR/vadr\""
echo "export VADRMODELDIR=\"\$VADRINSTALLDIR/vadr-models\""
echo "export VADRINFERNALDIR=\"\$VADRINSTALLDIR/infernal/src\""
echo "export VADREASELDIR=\"\$VADRINSTALLDIR/infernal/easel/miniapps\""
echo "export VADRBIOEASELDIR=\"\$VADRINSTALLDIR/Bio-Easel\""
echo "export VADRSEQUIPDIR=\"\$VADRINSTALLDIR/sequip\""
echo "export VADRBLASTDIR=\"\$VADRINSTALLDIR/ncbi-blast/bin\""
echo "export PERL5LIB=\"\$VADRSCRIPTSDIR\":\"\$VADRSEQUIPDIR\":\"\$VADRBIOEASELDIR/blib/lib\":\"\$VADRBIOEASELDIR/blib/arch\":\"\$PERL5LIB\""
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
echo "setenv VADRMODELDIR \"\$VADRINSTALLDIR/vadr-models\""
echo "setenv VADRINFERNALDIR \"\$VADRINSTALLDIR/infernal/src\""
echo "setenv VADREASELDIR \"\$VADRINSTALLDIR/infernal/easel/miniapps\""
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
