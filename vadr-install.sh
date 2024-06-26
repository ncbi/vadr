#!/bin/bash
# EPN, Wed Nov 20 06:06:22 2019
#
# vadr-install.sh
# A shell script for downloading and installing VADR and its dependencies.
# 
# usage: 
# to download and build files:
# vadr-install.sh <"linux" or "macosx-silicon" or "macosx-intel">
#
# or to only download files:
# vadr-install.sh <"linux" or "macosx-silicon" or "macosx-intel"> download
#
# or to only build files (after running in 'download' mode):
# vadr-install.sh <"linux" or "macosx-silicon" or "macosx-intel"> download
# 
# for example:
# vadr-install.sh linux
# 
# or
# vadr-install.sh macosx-silicon download
# vadr-install.sh macosx-silicon build
# 
# The following line will make the script fail if any commands fail
set -e

VADRINSTALLDIR=$PWD

# versions
VERSION="1.6.4"
# bio-easel (need this version info here only so we can check out correct easel branch in Bio-Easel/src)
BEVERSION="Bio-Easel-0.16"
# blast+
BVERSION="2.15.0"
# infernal
IVERSION="1.1.5"
# fasta
FVERSION="36.3.8h"
FVERSIONGIT="v36.3.8h_04-May-2020"
FVERSIONGITNOV="36.3.8h_04-May-2020"
# minimap2
MM2VERSIONGIT="v2.26"
MM2VERSIONGITNOV="2.26"
# dependency git tag
VVERSION="vadr-$VERSION"
# vadr models
MVERSION="1.2-1"
# hmmer (not needed in this release, we can use infernal's hmmer executables)
#HVERSION="3.4"

# set defaults
INPUTSYSTEM="?"
DOWNLOADORBUILD="both"

########################
# Validate correct usage
########################
# make sure correct number of cmdline arguments were used, exit if not
if [ "$#" -ne 1 ]; then
    if [ "$#" -ne 2 ]; then
        echo "Usage:"
        echo "To download and build:"
        echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\">"
        echo ""
        echo "or to only download files:"
        echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> download"
        echo ""
        echo "or to only build the software (after running in download mode):"
        echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> build"
        echo ""
        exit 1
    fi
fi

# make sure 1st argument is either "linux" or "macosx-silicon" or "macosx-intel"
if [ "$1" = "linux" ]; then
    INPUTSYSTEM="linux";
fi
if [ "$1" = "macosx-silicon" ]; then
    INPUTSYSTEM="macosx-silicon";
fi
if [ "$1" = "macosx-intel" ]; then
    INPUTSYSTEM="macosx-intel";
fi
if [ "$INPUTSYSTEM" = "?" ]; then 
    echo "Usage:"
    echo "To download and build:"
    echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\">"
    echo ""
    echo "or to only download files:"
    echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> download"
    echo ""
    echo "or to only build the software (after running in download mode):"
    echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> build"
    echo ""
    exit 1
fi

# make sure 2nd argument (if we have one) is "download" or "build"
if [ "$#" -eq 2 ]; then
    if [ "$2" = "download" ]; then
        DOWNLOADORBUILD="download";
    fi
    if [ "$2" = "build" ]; then
        DOWNLOADORBUILD="build";
    fi
    if [ "$DOWNLOADORBUILD" = "both" ]; then 
        echo "Usage:"
        echo "To download and build:"
        echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\">"
        echo ""
        echo "or to only download files:"
        echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> download"
        echo ""
        echo "or to only build the software (after running in download mode):"
        echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> build"
        echo ""
        exit 1
    fi
fi

########################################################
if [ "$DOWNLOADORBUILD" = "both" ]; then 
    echo "------------------------------------------------------------"
    echo "DOWNLOADING AND BUILDING VADR $VERSION"
    echo "------------------------------------------------------------"
fi
if [ "$DOWNLOADORBUILD" = "download" ]; then 
    echo "------------------------------------------------------------"
    echo "DOWNLOADING VADR $VERSION"
    echo "------------------------------------------------------------"
fi
if [ "$DOWNLOADORBUILD" = "build" ]; then 
    echo "------------------------------------------------------------"
    echo "BUILDING VADR $VERSION"
    echo "------------------------------------------------------------"
fi
if [ "$DOWNLOADORBUILD" != "download" ]; then 
    echo ""
    echo "************************************************************"
    echo "IMPORTANT: BEFORE YOU WILL BE ABLE TO RUN VADR SCRIPTS,"
    echo "YOU NEED TO FOLLOW THE INSTRUCTIONS OUTPUT AT THE END"
    echo "OF THIS SCRIPT TO UPDATE YOUR ENVIRONMENT VARIABLES."
    echo "************************************************************"
    echo ""
    echo "Determining current directory ... "
    echo "Set VADRINSTALLDIR as current directory ($VADRINSTALLDIR)."
fi


###########################################
# Download section
###########################################
if [ "$DOWNLOADORBUILD" != "build" ]; then
    echo "------------------------------------------------------------"
    # vadr
    echo "Downloading vadr ... "
    curl -k -L -o $VVERSION.zip https://github.com/ncbi/vadr/archive/$VVERSION.zip; unzip $VVERSION.zip; mv vadr-$VVERSION vadr; rm $VVERSION.zip
    # for a test build of a release, comment out above curl and uncomment block below
    # ------------------------------------------------------------
    #git clone https://github.com/ncbi/vadr.git vadr
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
    echo "------------------------------------------------------------"

    echo "Downloading Infernal version $IVERSION src distribution"
    curl -k -L -o infernal.tar.gz http://eddylab.org/infernal/infernal-$IVERSION.tar.gz
    tar xfz infernal.tar.gz
    rm infernal.tar.gz
    echo "------------------------------------------------------------"

    # download fasta source distribution from github
    echo "Downloading FASTA version $FVERSIONGIT src distribution"
    curl -k -L -o $FVERSIONGIT.zip https://github.com/wrpearson/fasta36/archive/$FVERSIONGIT.zip; unzip $FVERSIONGIT.zip; mv fasta36-$FVERSIONGITNOV fasta; rm $FVERSIONGIT.zip
    # patch Makefile with vadr specific changes and copy to expected name so 'build' step is linux/osx agnostic
    if [ "$INPUTSYSTEM" = "linux" ]; then
        patch fasta/make/Makefile.linux vadr/fasta-mods/vadr-fasta-Makefile.linux.patch
        cp fasta/make/Makefile.linux fasta/make/Makefile.vadr_install
    else 
        patch fasta/make/Makefile.os_x86_64 vadr/fasta-mods/vadr-fasta-Makefile.os_x86_64.patch
        cp fasta/make/Makefile.os_x86_64 fasta/make/Makefile.vadr_install
    fi
    # patch defs.h with vadr specific changes
    patch fasta/src/defs.h vadr/fasta-mods/vadr-fasta-defs.patch
    echo "------------------------------------------------------------"

    # download minimap2 source distribution from github
    echo "Downloading minimap2 version $MM2VERSIONGIT src distribution"
    curl -k -L -o $MM2VERSIONGIT.zip https://github.com/lh3/minimap2/archive/$MM2VERSIONGIT.zip; unzip $MM2VERSIONGIT.zip; mv minimap2-$MM2VERSIONGITNOV minimap2; rm $MM2VERSIONGIT.zip
    echo "------------------------------------------------------------"
    
    # download blast binaries
    if [ "$INPUTSYSTEM" = "linux" ]; then
        echo "Downloading BLAST version $BVERSION for Linux"
        curl -k -L -o blast.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BVERSION/ncbi-blast-$BVERSION+-x64-linux.tar.gz
    else 
        echo "Downloading BLAST version $BVERSION for Mac/OSX"
        curl -k -L -o blast.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BVERSION/ncbi-blast-$BVERSION+-x64-macosx.tar.gz
    fi
    tar xfz blast.tar.gz
    rm blast.tar.gz
    mv ncbi-blast-$BVERSION+ ncbi-blast
    echo "------------------------------------------------------------"

    # download vadr models, calici and flavi model sets only
    for v in calici flavi; do 
        echo "Downloading VADR ${v}viridae models ($MVERSION) ... "
        curl -k -L -o vadr-models-$v.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/${v}viridae/$MVERSION/vadr-models-$v-$MVERSION.tar.gz
        tar xfz vadr-models-$v.tar.gz
        mv vadr-models-$v-$MVERSION vadr-models-$v
        rm vadr-models-$v.tar.gz
    done
    echo "------------------------------------------------------------"
fi

if [ "$DOWNLOADORBUILD" = "download" ]; then
    echo ""
    echo ""
    echo "********************************************************"
    echo "Downloads finished successfully."
    echo "You will need to build the software before you can use it."
    echo "To do that, run this script in 'build' mode with the command:"
    echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> build"
    echo "********************************************************"
    echo ""
fi

###########################################
# Build section
###########################################
if [ "$DOWNLOADORBUILD" != "download" ]; then

    # Build Bio-Easel:
    if [ ! -d Bio-Easel ]; then
        echo ""
        echo "ERROR: Bio-Easel dir does not exist"
        if [ "$DOWNLOADORBUILD" = "build" ]; then
            echo ""
            echo "This may be because you did not yet run this script in download mode from this directory,"
            echo "which is required prior to running in build mode. To do that, execute:"
            echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> download"
            echo ""
            exit 1
        fi
    fi
    echo "------------------------------------------------------------"
    echo "Building Bio-Easel ... "
    cd Bio-Easel
    perl Makefile.PL
    make
    make test
    cd ..
    echo "Finished building Bio-Easel."
    echo "------------------------------------------------------------"
    
    # Build infernal:
    if [ ! -d infernal-$IVERSION ]; then
        echo ""
        echo "ERROR: infernal-$IVERSION dir does not exist"
        if [ "$DOWNLOADORBUILD" = "build" ]; then
            echo ""
            echo "This may be because you did not yet run this script in download mode from this directory,"
            echo "which is required prior to running in build mode. To do that, execute:"
            echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> download"
            echo ""
            exit 1
        fi
        exit 1
    fi
    echo "------------------------------------------------------------"
    echo "Building Infernal ... "
    mv infernal-$IVERSION infernal
    cd infernal
    mkdir binaries
    sh ./configure --bindir=$PWD/binaries --prefix=$PWD
    make
    make install
    (cd easel/miniapps; make install)
    (cd hmmer; make install)
    cd ..
    echo "Finished building Infernal."
    echo "------------------------------------------------------------"

    # Build FASTA:
    if [ ! -d fasta ]; then
        echo ""
        echo "ERROR: fasta dir does not exist"
        if [ "$DOWNLOADORBUILD" = "build" ]; then
            echo ""
            echo "This may be because you did not yet run this script in download mode from this directory,"
            echo "which is required prior to running in build mode. To do that, execute:"
            echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> download"
            echo ""
            exit 1
        fi
        exit 1
    fi
    echo "------------------------------------------------------------"
    echo "Building FASTA ... "
    cd fasta/src
    # note, download step copied either ../make/Makefile.linux_sse2 or ../make/Makefile.os_x86_64 to ../make/Makefile.vadr_install
    make -f ../make/Makefile.vadr_install all
    cd ../../
    echo "Finished building FASTA."
    echo "------------------------------------------------------------"

    # Build minimap2:
    if [ ! -d minimap2 ]; then
        echo ""
        echo "ERROR: minimap2 dir does not exist"
        if [ "$DOWNLOADORBUILD" = "build" ]; then
            echo ""
            echo "This may be because you did not yet run this script in download mode from this directory,"
            echo "which is required prior to running in build mode. To do that, execute:"
            echo "  $0 <\"linux\" or \"macosx-silicon\" or \"macosx-intel\"> download"
            echo ""
            exit 1
        fi
        exit 1
    fi
    echo "------------------------------------------------------------"
    echo "Building minimap2 ... "
    cd minimap2
    if [ "$INPUTSYSTEM" = "macosx-silicon" ]; then
        make arm_neon=1 aarch64=1
    fi
    if [ "$INPUTSYSTEM" != "macosx-silicon" ]; then
        make
    fi
    cd ../../
    echo "Finished building minimap2."
    echo "------------------------------------------------------------"
    
    ###############################################
    # Message about setting environment variables
    ###############################################
    echo ""
    echo ""
    echo "********************************************************"
    echo "The final step is to update your environment variables."
    echo "(See https://github.com/ncbi/vadr/blob/$VERSION/documentation/install.md for more information.)"
    echo ""
    echo "If you are using the bash or zsh shell (zsh is default in MacOS/X as"
    echo "of v10.15 (Catalina)), add the following lines to the end of your"
    echo "'.bashrc' or '.zshrc' file in your home directory:"
    echo ""
    echo "export VADRINSTALLDIR=\"$VADRINSTALLDIR\""
    echo "export VADRSCRIPTSDIR=\"\$VADRINSTALLDIR/vadr\""
    echo "export VADRMODELDIR=\"\$VADRINSTALLDIR/vadr-models-calici\""
    echo "export VADRINFERNALDIR=\"\$VADRINSTALLDIR/infernal/binaries\""
    echo "export VADREASELDIR=\"\$VADRINSTALLDIR/infernal/binaries\""
    echo "export VADRHMMERDIR=\"\$VADRINSTALLDIR/infernal/binaries\""
    echo "export VADRBIOEASELDIR=\"\$VADRINSTALLDIR/Bio-Easel\""
    echo "export VADRSEQUIPDIR=\"\$VADRINSTALLDIR/sequip\""
    echo "export VADRBLASTDIR=\"\$VADRINSTALLDIR/ncbi-blast/bin\""
    echo "export VADRFASTADIR=\"\$VADRINSTALLDIR/fasta/bin\""
    echo "export VADRMINIMAP2DIR=\"\$VADRINSTALLDIR/minimap2\""
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
    echo "setenv VADRMODELDIR \"\$VADRINSTALLDIR/vadr-models-calici\""
    echo "setenv VADRINFERNALDIR \"\$VADRINSTALLDIR/infernal/binaries\""
    echo "setenv VADRHMMERDIR \"\$VADRINSTALLDIR/infernal/binaries\""
    echo "setenv VADREASELDIR \"\$VADRINSTALLDIR/infernal/binaries\""
    echo "setenv VADRBIOEASELDIR \"\$VADRINSTALLDIR/Bio-Easel\""
    echo "setenv VADRSEQUIPDIR \"\$VADRINSTALLDIR/sequip\""
    echo "setenv VADRBLASTDIR \"\$VADRINSTALLDIR/ncbi-blast/bin\""
    echo "setenv VADRFASTADIR \"\$VADRINSTALLDIR/fasta/bin\""
    echo "setenv VADRMINIMAP2DIR=\"\$VADRINSTALLDIR/minimap2\""
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
fi
