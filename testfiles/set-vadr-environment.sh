#!/bin/bash
# The previous line forces this script to be run with bash, regardless of user's 
# shell.
#
# EPN, Wed Oct 16 10:42:43 2019
#
# set-vadr-environment.sh <directory where checkout-source.sh and build-binaries.sh were run from>
# 
# for example:
# set-vadr-environment.sh /home/nawrocki/vadr-install
#
# The following line will make the script fail if any commands fail
set -e

# make sure correct number of cmdline arguments were used, exit if not
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <directory where checkout-source.sh and build-binaries.sh were run from>"
  exit 1
fi

if [ ! "$BASH_VERSION" ] ; then
    echo "Please do not use sh to run this script ($0), just execute it directly" 1>&2
    exit 1
fi

VADRINSTALLDIR=$1
VADRSCRIPTDIR="$VADRINSTALLDIR/vadr"
VADRMODELDIR="$VADRINSTALLDIR/vadr-models"
VADRINFERNALDIR="$VADRINSTALLDIR/infernal-dev/src"
VADREASELDIR="$VADRINSTALLDIR/infernal-dev/easel/miniapps"
VADRBIOEASELDIR="$VADRINSTALLDIR/Bio-Easel"
VADRSEQUIPDIR="$VADRINSTALLDIR/sequip"
VADRBLASTDIR="$VADRINSTALLDIR/ncbi-blast/bin"

# make sure directories exist
if [ ! -d $VADRINSTALLDIR ]; then 
    echo "ERROR: VADRINSTALLDIR directory ($VADRINSTALLDIR) does not exist";
    exit 1
fi
if [ ! -d $VADRMODELDIR ]; then 
    echo "ERROR: VADRMODELDIR directory ($VADRMODELDIR) does not exist";
    exit 1
fi
if [ ! -d $VADRINFERNALDIR ]; then 
    echo "ERROR: VADRINFERNALDIR directory ($VADRINFERNALDIR) does not exist";
    exit 1
fi
if [ ! -d $VADREASELDIR ]; then 
    echo "ERROR: VADREASELDIR directory ($VADREASELDIR) does not exist";
    exit 1
fi
if [ ! -d $VADRBIOEASELDIR ]; then 
    echo "ERROR: VADRBIOEASELDIR directory ($VADRBIOEASELDIR) does not exist";
    exit 1
fi
if [ ! -d $VADRSEQUIPDIR ]; then 
    echo "ERROR: VADRSEQUIPDIR directory ($VADRSEQUIPDIR) does not exist";
    exit 1
fi
if [ ! -d $VADRBLASTDIR ]; then 
    echo "ERROR: VADRBLASTDIR directory ($VADRBLASTDIR) does not exist";
    exit 1
fi

export VADRINSTALLDIR=\"$VADRINSTALLDIR\"
export VADRSCRIPTSDIR=\"\$VADRSCRIPTSDIR\"
export VADRMODELDIR=\"$VADRMODELDIR\"
export VADRINFERNALDIR=\"\$VADRINFERNALDIR\"
export VADREASELDIR=\"\$VADREASELDIR\"
export VADRBIOEASELDIR=\"\$VADRBIOEASELDIR\"
export VADRSEQUIPDIR=\"\$VADRSEQUIPDIR\"
export VADRBLASTDIR=\"\$VADRBLASTDIR\"
export PERL5LIB=\"\$VADRSCRIPTSDIR\":\"\$VADRSEQUIPDIR\":\"\$VADRBIOEASELDIR/blib/lib\":\"\$VADRBIOEASELDIR/blib/arch\":\"\$PERL5LIB\"
export PATH=\"\$VADRSCRIPTSDIR\":\"\$PATH\"
