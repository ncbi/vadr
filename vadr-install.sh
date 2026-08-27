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
# or to only download files, and minimize number of models downloaded:
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
# Requirements: a C and C++ compiler, make, perl, curl, unzip and git.
# Installing R2DT, which is used only by 'v-annotate.pl --draw_r2dt',
# additionally requires python3 (see R2DTMINPYTHON below) and network access to
# PyPI. R2DT is the only optional dependency: if it cannot be installed, this
# script warns and installs the rest of VADR.
#
# The following line will make the script fail if any commands fail
set -e

VADRINSTALLDIR=$PWD

# versions
VERSION="1.7"
# bio-easel (need this version info here only so we can check out correct easel branch in Bio-Easel/src)
BEVERSION="Bio-Easel-0.17"
# blast+
BVERSION="2.17.0"
# infernal
IVERSION="1.1.5"
# fasta
FVERSION="36.3.8h"
FVERSIONGIT="v36.3.8h_04-May-2020"
FVERSIONGITNOV="36.3.8h_04-May-2020"
# minimap2
MM2VERSIONGIT="v2.30"
MM2VERSIONGITNOV="2.30"
# dependency git tag
VVERSION="vadr-$VERSION"
# vadr models
CALICIVERSION="1.2-1"
FLAVIVERSION="1.7-1"
CORONAVERSION="1.3-3"
SARSCOV2VERSION="1.3-2"
FLUVERSION="1.6.3-2"
RSVVERSION="1.5-2"
MPXVVERSION="1.4.2-1"
# hmmer (not needed in this release, we can use infernal's hmmer executables)
#HVERSION="3.4"
# R2DT (used only by 'v-annotate.pl --draw_r2dt')
#
# TODO: BEFORE VADR 1.7.1 IS RELEASED, REPLACE THIS COMMIT WITH AN R2DT RELEASE
# TODO: TAG. R2DT has no release that contains the drawing code VADR requires,
# TODO: so this pins a reviewed commit on the R2DT 'develop' branch instead.
# TODO: A release tag must replace it before 1.7.1 ships.
R2DTVERSION="a3ad3588e5015ff5cbbc7c56cf6e8f57dce2aa93"
# traveler, the program R2DT uses to render diagrams. This commit, and the patch
# applied on top of it, are both taken from R2DT's own base_image/Dockerfile and
# base_image/patches/, so they are properties of R2DTVERSION and must be updated
# together with it.
TVERSION="78c1cd63dcf40be83d97d07b8eb9eb6a9dc63f7f"
TPATCH="traveler-pr20-pk-per-bp.patch"
# jiffy-infernal-hmmer-scripts, also required by R2DT's drawing path, also
# pinned by R2DT's base_image/Dockerfile
JVERSION="31d6c3b826d432a30620507830749cee58e15e68"
# minimum python3 version R2DT's python requirements can be installed on
R2DTMINPYTHON="3.9"

# set defaults
INPUTSYSTEM="?"
DOWNLOADORBUILD="both"

# R2DT is the one dependency this script treats as optional: it is used only by
# 'v-annotate.pl --draw_r2dt', it is the only dependency that needs python3 and
# network access to PyPI, and everything else VADR does works without it. If any
# part of the R2DT install fails we report it loudly and finish installing the
# rest of VADR, rather than leaving the user with no VADR at all. R2DTFAILED is
# set to 1 if that happens, so that a warning can be repeated at the very end.
R2DTFAILED=0
R2DTDIR="$VADRINSTALLDIR/R2DT"
R2DTVENVDIR="$VADRINSTALLDIR/r2dt-venv"
R2DTTRAVELERDIR="$VADRINSTALLDIR/traveler"
R2DTJIFFYDIR="$VADRINSTALLDIR/jiffy-infernal-hmmer-scripts"
R2DTLOG="$VADRINSTALLDIR/r2dt-install.log"

# r2dt_failure(): report that the R2DT install did not complete, explain what
# still works and what does not, and make sure nothing is left behind that could
# be mistaken for a working R2DT install. In particular the site config file
# $R2DTDIR/r2dt-vadr-env.sh is only ever written once every step has succeeded,
# so a partial install cannot silently produce empty diagrams.
r2dt_failure () {
    R2DTFAILED=1
    echo ""
    echo "************************************************************"
    echo "WARNING: the R2DT install did not complete."
    echo ""
    echo "Reason: $1"
    echo ""
    echo "The rest of VADR is unaffected and this script will continue."
    echo "Everything except 'v-annotate.pl --draw_r2dt' will work normally."
    echo "R2DT is only used to draw secondary structure diagrams."
    echo ""
    echo "To retry after fixing the problem above, remove these if they exist"
    echo "and rerun this script:"
    echo "  $R2DTDIR"
    echo "  $R2DTVENVDIR"
    echo "  $R2DTTRAVELERDIR"
    echo "  $R2DTJIFFYDIR"
    echo "************************************************************"
    echo ""
    # remove the site config file if a previous run left one: without it,
    # v-annotate.pl --draw_r2dt fails with a clear error instead of running
    # r2dt.py in an environment that cannot work.
    rm -f "$R2DTDIR/r2dt-vadr-env.sh"
}

# r2dt_find_python(): echo the path of a python3 that R2DT's requirements can be
# installed on, or echo nothing if there is none. The user's own python3 is
# preferred when it is new enough, otherwise the newest suitable python3.X on
# PATH is used. Set VADRPYTHON to override this entirely.
r2dt_python_is_new_enough () {
    $1 -c "import sys; sys.exit(0 if sys.version_info[:2] >= tuple(int(x) for x in \"$R2DTMINPYTHON\".split(\".\")) else 1)" > /dev/null 2>&1
}
r2dt_find_python () {
    if [ "$VADRPYTHON" != "" ]; then
        # an explicitly requested interpreter is still checked, so that a stale
        # or wrong VADRPYTHON is reported as such instead of failing later with
        # a less obvious error
        if r2dt_python_is_new_enough "$VADRPYTHON"; then
            echo "$VADRPYTHON"
        fi
        return
    fi
    for p in python3 python3.13 python3.12 python3.11 python3.10 python3.9; do
        if command -v $p > /dev/null 2>&1; then
            if r2dt_python_is_new_enough $p; then
                command -v $p
                return
            fi
        fi
    done
}

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

    # download vadr models
    for v in calici; do 
        echo "Downloading VADR $v models ($CALICIVERSION) ... "
        curl -k -L -o vadr-models-$v.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/${v}viridae/$CALICIVERSION/vadr-models-$v-$CALICIVERSION.tar.gz
        tar xfz vadr-models-$v.tar.gz
        mv vadr-models-$v-$CALICIVERSION vadr-models-$v
        rm vadr-models-$v.tar.gz
    done
    for v in flavi; do 
        echo "Downloading VADR $v models ($FLAVIVERSION) ... "
        curl -k -L -o vadr-models-$v.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/${v}viridae/$FLAVIVERSION/vadr-models-$v-$FLAVIVERSION.tar.gz
        tar xfz vadr-models-$v.tar.gz
        mv vadr-models-$v-$FLAVIVERSION vadr-models-$v
        rm vadr-models-$v.tar.gz
    done
    for v in corona; do 
        echo "Downloading VADR $v models ($CORONAVERSION) ... "
        curl -k -L -o vadr-models-$v.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/${v}viridae/$CORONAVERSION/vadr-models-$v-$CORONAVERSION.tar.gz
        tar xfz vadr-models-$v.tar.gz
        mv vadr-models-$v-$CORONAVERSION vadr-models-$v
        rm vadr-models-$v.tar.gz
    done
    for v in sarscov2; do 
        echo "Downloading VADR $v models ($SARSCOV2VERSION) ... "
        curl -k -L -o vadr-models-$v.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/$v/$SARSCOV2VERSION/vadr-models-$v-$SARSCOV2VERSION.tar.gz
        tar xfz vadr-models-$v.tar.gz
        mv vadr-models-$v-$SARSCOV2VERSION vadr-models-$v
        rm vadr-models-$v.tar.gz
    done
    for v in flu; do 
        echo "Downloading VADR $v models ($FLUVERSION) ... "
        curl -k -L -o vadr-models-$v.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/$v/$FLUVERSION/vadr-models-$v-$FLUVERSION.tar.gz
        tar xfz vadr-models-$v.tar.gz
        mv vadr-models-$v-$FLUVERSION vadr-models-$v
        rm vadr-models-$v.tar.gz
    done
    for v in rsv; do 
        echo "Downloading VADR $v models ($RSVVERSION) ... "
        curl -k -L -o vadr-models-$v.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/$v/$RSVVERSION/vadr-models-$v-$RSVVERSION.tar.gz
        tar xfz vadr-models-$v.tar.gz
        mv vadr-models-$v-$RSVVERSION vadr-models-$v
        rm vadr-models-$v.tar.gz
    done
    for v in mpxv; do 
        echo "Downloading VADR $v models ($MPXVVERSION) ... "
        curl -k -L -o vadr-models-$v.tar.gz https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/$v/$MPXVVERSION/vadr-models-$v-$MPXVVERSION.tar.gz
        tar xfz vadr-models-$v.tar.gz
        mv vadr-models-$v-$MPXVVERSION vadr-models-$v
        rm vadr-models-$v.tar.gz
    done
    echo "------------------------------------------------------------"

    ###########################################
    # R2DT download (optional dependency)
    ###########################################
    # Everything that needs the network happens here, including installing
    # R2DT's python requirements into a virtual environment, so that a later
    # run in 'build' mode needs a compiler but not a network connection, the
    # same as for every other dependency.
    echo "Downloading R2DT (needed only by 'v-annotate.pl --draw_r2dt') ... "
    echo "  (see $R2DTLOG for the full log of this step)"
    R2DTPYTHON=`r2dt_find_python`
    if [ "$R2DTPYTHON" = "" ]; then
        if [ "$VADRPYTHON" != "" ]; then
            r2dt_failure "VADRPYTHON is set to '$VADRPYTHON', which is not a working python3 of version $R2DTMINPYTHON or later. R2DT's python requirements cannot be installed with it. Set VADRPYTHON to the full path of a newer python3, or unset it to let this script look for one, then rerun this script."
        else
            r2dt_failure "no python3 of version $R2DTMINPYTHON or later was found on your PATH. R2DT's python requirements cannot be installed without one. Install a newer python3, or set the VADRPYTHON environment variable to the full path of one, then rerun this script."
        fi
    else
        echo "  Using python: $R2DTPYTHON (`$R2DTPYTHON --version 2>&1`)"
        # NOTE: the subshell below is deliberately NOT run as the condition of
        # an 'if'. Inside an 'if' condition bash ignores errexit for the whole
        # command, including a 'set -e' the subshell runs itself, so a step
        # failing halfway through would not stop the steps after it. Disabling
        # errexit around a plain subshell and testing $? afterwards does behave
        # as intended.
        set +e
        ( set -e
             exec > "$R2DTLOG" 2>&1
             cd "$VADRINSTALLDIR"
             rm -rf "$R2DTDIR" "$R2DTVENVDIR" "$R2DTTRAVELERDIR" "$R2DTJIFFYDIR"

             echo "== cloning R2DT at $R2DTVERSION =="
             git clone https://github.com/r2dt-bio/R2DT.git R2DT
             (cd R2DT && git checkout $R2DTVERSION)

             echo "== cloning jiffy-infernal-hmmer-scripts at $JVERSION =="
             git clone https://github.com/nawrockie/jiffy-infernal-hmmer-scripts.git jiffy-infernal-hmmer-scripts
             (cd jiffy-infernal-hmmer-scripts && git checkout $JVERSION)
             # R2DT runs these as commands found on PATH, but not all of them
             # are executable in the repository, so make them all executable,
             # exactly as R2DT's own base_image/Dockerfile does. Without this,
             # r2dt.py fails with "Permission denied" on
             # ali-pfam-lowercase-rf-gap-columns.pl and draws nothing.
             chmod +x jiffy-infernal-hmmer-scripts/*.pl

             echo "== cloning traveler at $TVERSION =="
             git clone https://github.com/cusbg/traveler.git traveler
             (cd traveler && git checkout $TVERSION)
             echo "== applying R2DT's traveler patch $TPATCH =="
             (cd traveler && git apply --verbose "$R2DTDIR/base_image/patches/$TPATCH")
             echo "== installing the <bits/stdc++.h> portability shim =="
             mkdir -p traveler/src/include/bits
             cp vadr/traveler-mods/bits-stdcxx-shim.h traveler/src/include/bits/stdc++.h

             echo "== creating python virtual environment =="
             "$R2DTPYTHON" -m venv "$R2DTVENVDIR"
             echo "== installing R2DT python requirements =="
             "$R2DTVENVDIR/bin/python" -m pip install --upgrade pip
             "$R2DTVENVDIR/bin/python" -m pip install -r "$R2DTDIR/requirements-minimal.txt"
             "$R2DTVENVDIR/bin/python" -m pip list
        )
        R2DTSTATUS=$?
        set -e
        if [ $R2DTSTATUS -eq 0 ]; then
            echo "Finished downloading R2DT."
        else
            r2dt_failure "one of the R2DT download steps failed. The last lines of $R2DTLOG were:
`tail -20 "$R2DTLOG" 2>/dev/null`"
        fi
    fi
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

    ###########################################
    # R2DT build (optional dependency)
    ###########################################
    # Nothing here needs the network: the download step already cloned the
    # sources and populated the python virtual environment.
    cd "$VADRINSTALLDIR"
    if [ ! -d "$R2DTDIR" ] || [ ! -d "$R2DTTRAVELERDIR" ] || [ ! -d "$R2DTJIFFYDIR" ] || [ ! -d "$R2DTVENVDIR" ]; then
        r2dt_failure "the R2DT source directories are not all present, so the download step for R2DT either did not run or did not complete."
    else
        echo "------------------------------------------------------------"
        echo "Building traveler (the renderer R2DT uses) ... "
        echo "  (appending to $R2DTLOG)"
        set +e
        ( set -e
             exec >> "$R2DTLOG" 2>&1
             echo "== building traveler =="
             cd "$R2DTTRAVELERDIR/src"
             make build
             test -x "$R2DTTRAVELERDIR/bin/traveler"
             # R2DT locates traveler's utils/ directory (infernal2mapping.py,
             # enrich_json.py, json2svg.py) relative to the traveler executable
             # it finds on PATH, so bin/ and utils/ must stay siblings. They are,
             # because traveler is built in place inside its own checkout.
             test -f "$R2DTTRAVELERDIR/utils/infernal2mapping.py"
        )
        R2DTSTATUS=$?
        set -e
        if [ $R2DTSTATUS -ne 0 ]; then
            r2dt_failure "the traveler build failed. The last lines of $R2DTLOG were:
`tail -20 "$R2DTLOG" 2>/dev/null`"
        else
            echo "Finished building traveler."
            echo "------------------------------------------------------------"
            echo "Writing R2DT site configuration ... "
            # v-annotate.pl --draw_r2dt sources $R2DT_DIR/r2dt-vadr-env.sh, if it
            # exists, immediately before running r2dt.py. Writing it here is what
            # makes --draw_r2dt work with no further configuration by the user.
            # It is written last, and only after every other R2DT step has
            # succeeded, so that a partial install never looks like a working one.
            # It is sourced with POSIX '.', so keep it POSIX compatible.
            cat > "$R2DTDIR/r2dt-vadr-env.sh" << EOF
# Site configuration for 'v-annotate.pl --draw_r2dt', written by vadr-install.sh.
#
# v-annotate.pl --draw_r2dt sources this file with POSIX '.' immediately before
# invoking r2dt.py, so that the environment r2dt.py needs is supplied here
# rather than hardcoded in VADR. If this file is removed, --draw_r2dt will only
# work if 'python \$R2DT_DIR/r2dt.py' already works in the calling environment.
#
# The paths below were determined when VADR was installed. If this installation
# is moved, they must be updated.

# The python virtual environment holding R2DT's python requirements is put
# first, so that 'python' resolves to it. The virtual environment is put on PATH
# rather than activated, which keeps this file POSIX compatible: a virtual
# environment python works by virtue of where it sits on PATH.
#
# The remaining directories hold the external programs R2DT calls:
#   Infernal (cmalign, cmbuild, and the esl-* Easel miniapps), the Bio-Easel
#   scripts, the jiffy Infernal/HMMER scripts, and traveler.
export PATH="$R2DTVENVDIR/bin:$VADRINSTALLDIR/infernal/binaries:$VADRINSTALLDIR/Bio-Easel/scripts:$R2DTJIFFYDIR:$R2DTTRAVELERDIR/bin:\$PATH"

# The Bio-Easel perl modules, needed by the Bio-Easel and jiffy scripts above.
export PERL5LIB="$VADRINSTALLDIR/Bio-Easel/blib/lib:$VADRINSTALLDIR/Bio-Easel/blib/arch:\$PERL5LIB"
EOF
            echo "Wrote $R2DTDIR/r2dt-vadr-env.sh"
            echo "Finished installing R2DT."
            echo "------------------------------------------------------------"
        fi
    fi

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
    echo "export VADRCONFIGFILE=\"\$VADRSCRIPTSDIR/vadr.config\""
    echo "export VADRMODELDIR=\"\$VADRINSTALLDIR/vadr-models-calici\""
    echo "export VADRINFERNALDIR=\"\$VADRINSTALLDIR/infernal/binaries\""
    echo "export VADREASELDIR=\"\$VADRINSTALLDIR/infernal/binaries\""
    echo "export VADRHMMERDIR=\"\$VADRINSTALLDIR/infernal/binaries\""
    echo "export VADRBIOEASELDIR=\"\$VADRINSTALLDIR/Bio-Easel\""
    echo "export VADRSEQUIPDIR=\"\$VADRINSTALLDIR/sequip\""
    echo "export VADRBLASTDIR=\"\$VADRINSTALLDIR/ncbi-blast/bin\""
    echo "export VADRFASTADIR=\"\$VADRINSTALLDIR/fasta/bin\""
    echo "export VADRMINIMAP2DIR=\"\$VADRINSTALLDIR/minimap2\""
    if [ "$R2DTFAILED" = "0" ]; then
        echo "export R2DT_DIR=\"\$VADRINSTALLDIR/R2DT\""
    fi
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
    echo "setenv VADRCONFIGFILE \"\$VADRSCRIPTSDIR/vadr.config\""
    echo "setenv VADRMODELDIR \"\$VADRINSTALLDIR/vadr-models-calici\""
    echo "setenv VADRINFERNALDIR \"\$VADRINSTALLDIR/infernal/binaries\""
    echo "setenv VADRHMMERDIR \"\$VADRINSTALLDIR/infernal/binaries\""
    echo "setenv VADREASELDIR \"\$VADRINSTALLDIR/infernal/binaries\""
    echo "setenv VADRBIOEASELDIR \"\$VADRINSTALLDIR/Bio-Easel\""
    echo "setenv VADRSEQUIPDIR \"\$VADRINSTALLDIR/sequip\""
    echo "setenv VADRBLASTDIR \"\$VADRINSTALLDIR/ncbi-blast/bin\""
    echo "setenv VADRFASTADIR \"\$VADRINSTALLDIR/fasta/bin\""
    echo "setenv VADRMINIMAP2DIR=\"\$VADRINSTALLDIR/minimap2\""
    if [ "$R2DTFAILED" = "0" ]; then
        echo "setenv R2DT_DIR \"\$VADRINSTALLDIR/R2DT\""
    fi
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
    if [ "$R2DTFAILED" = "1" ]; then
        echo "************************************************************"
        echo "WARNING: VADR was installed, but R2DT was NOT."
        echo ""
        echo "'v-annotate.pl --draw_r2dt' will not work. Everything else"
        echo "will. See the warning earlier in this script's output, and"
        echo "$R2DTLOG, for the reason."
        echo "************************************************************"
        echo ""
    fi
fi
