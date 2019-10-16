#!/bin/bash
# The previous line forces this script to be run with bash, regardless of user's 
# shell.
#
# EPN, Wed Oct 16 10:42:43 2019
#
# generate-vadr-environment-commands.sh <directory where checkout-source.sh and build-binaries.sh were run from> <your shell: 'bash' or 'csh'>
# 
# for example:
# generate-vadr-environment-commands.sh /home/nawrocki/vadr-install bash
#
# The following line will make the script fail if any commands fail
set -e

# make sure correct number of cmdline arguments were used, exit if not
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <directory where checkout-vadr-source.sh and build-vadr-binaries.sh were run from> <your shell: 'bash' or 'csh'"
  exit 1
fi

if [ ! "$BASH_VERSION" ] ; then
    echo "Please do not use sh to run this script ($0), just execute it directly" 1>&2
    exit 1
fi

# check 2nd cmdline arg is 'bash' or 'csh'
INPUTSHELL="?"
if [ $2 == "bash" ]; then
    INPUTSHELL="bash";
fi
if [ $2 == "csh" ]; then
    INPUTSHELL="csh";
fi
if [ $INPUTSHELL == "?" ]; then 
  echo "Usage: $0 <directory where checkout-vadr-source.sh and build-vadr-binaries.sh were run from> <your shell: 'bash' or 'csh'"
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

if [ $INPUTSHELL == "bash" ]; then
    echo "export VADRINSTALLDIR=\"$VADRINSTALLDIR\""
    echo "export VADRSCRIPTSDIR=\"\$VADRINSTALLDIR/vadr\""
    echo "export VADRMODELDIR=\"\$VADRINSTALLDIR/vadr-models\""
    echo "export VADRINFERNALDIR=\"\$VADRINSTALLDIR/infernal-dev/src\""
    echo "export VADREASELDIR=\"\$VADRINSTALLDIR/infernal-dev/easel/miniapps\""
    echo "export VADRBIOEASELDIR=\"\$VADRINSTALLDIR/Bio-Easel\""
    echo "export VADRSEQUIPDIR=\"\$VADRINSTALLDIR/sequip\""
    echo "export VADRBLASTDIR=\"\$VADRINSTALLDIR/ncbi-blast/bin\""
    echo "export PERL5LIB=\"\$VADRSCRIPTSDIR\":\"\$VADRSEQUIPDIR\":\"\$VADRBIOEASELDIR/blib/lib\":\"\$VADRBIOEASELDIR/blib/arch\":\"\$PERL5LIB\""
    echo "export PATH=\"\$VADRSCRIPTSDIR\":\"\$PATH\""
fi
if [ $INPUTSHELL == "csh" ]; then
    echo "setenv VADRINSTALLDIR \"$VADRINSTALLDIR\""
    echo "setenv VADRSCRIPTSDIR \"\$VADRINSTALLDIR/vadr\""
    echo "setenv VADRMODELDIR \"\$VADRINSTALLDIR/vadr-models\""
    echo "setenv VADRINFERNALDIR \"\$VADRINSTALLDIR/infernal-dev/src\""
    echo "setenv VADREASELDIR \"\$VADRINSTALLDIR/infernal-dev/easel/miniapps\""
    echo "setenv VADRBIOEASELDIR \"\$VADRINSTALLDIR/Bio-Easel\""
    echo "setenv VADRSEQUIPDIR \"\$VADRINSTALLDIR/sequip\""
    echo "setenv VADRBLASTDIR \"\$VADRINSTALLDIR/ncbi-blast/bin\""
    echo "setenv PERL5LIB \"\$VADRSCRIPTSDIR\":\"\$VADRSEQUIPDIR\":\"\$VADRBIOEASELDIR/blib/lib\":\"\$VADRBIOEASELDIR/blib/arch\":\"\$PERL5LIB\""
    echo "setenv PATH \"\$VADRSCRIPTSDIR\":\"\$PATH\""
fi
