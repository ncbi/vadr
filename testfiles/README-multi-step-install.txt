EPN, Wed Oct 16 10:45:13 2019

Instructions for installing VADR in multiple steps.

============
INSTALLATION
============

Choose a directory to install in. All source code, executables, and
VADR models (several directories) will be placed there. 

Move to that directory and copy the following files there:
checkout-vadr-source.sh
build-vadr-binaries.sh
set-vadr-environment.sh
generate-vadr-environment-commands.sh

Step 1. Checkout source with 'checkout-vadr-source.sh' supplying 3
command line arguments: <vadr version> <vadr model version> <blast+ version>

$ checkout-vadr-source.sh 0.991 0.991.1 2.9.0

Step 2. Download VADR cox1 models (if desired) with
'download-vadr-cox1-models.sh' with 1 command line argument: <vadr
cox1 model version (e.g. 0.971.1)> 

$ download-vadr-cox1-models.sh 0.971.1

Step 3. Build binaries with 'build-vadr-binaries.sh'

$ build-vadr-binaries.sh 

=============================
SETTING ENVIRONMENT VARIABLES 
=============================

Before you can use VADR you need to set specific environment
variables. 

To set the environment automatically whenever you login to your home
directory, you can add lines to the .bashrc or .cshrc files in your
home directory that get run whenever you login.

The script 'generate-environment-commands.sh' will output these
commands for you, and then you will need to copy them and paste them
in your .bashrc or .cshrc file. 'generate-environment-commands.sh'
takes two command line arguments. <path-to-current-directory> is the
same argument used above in 'Setting environment for current shell
session'. The second command line argument must be either 'bash' or
'csh' depending on whether you use the bash or C shell. To determine
what shell you use, execute:

$ echo $SHELL

Example command if you use bash shell:
$ generate-environment-commands.sh <path-to-current-directory> bash

Example command if you use C shell or tcsh:

$ generate-environment-commands.sh <path-to-current-directory> csh

The output of this command will be a bunch of lines like this:

export VADRMODELDIR="$VADRINSTALLDIR/vadr-models"

OR

setenv VADRMODELDIR "$VADRINSTALLDIR/vadr-models"

Take all of these lines and add them to your .bashrc or .cshrc
file. Then do 'source ~/.bashrc' or 'source ~/.cshrc'. After that, you
should be able to run the vadr scripts.

To test this, do:

$ which v-build.pl

The output should be:

<path-to-current-directory>/vadr/v-build.pl



