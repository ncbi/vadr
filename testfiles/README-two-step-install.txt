EPN, Wed Oct 16 10:45:13 2019

Instructions for installing VADR in two steps.

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

$ checkout-vadr-source.sh 0.991 0.991-1 2.9.0

Step 2. Build binaries with 'build-vadr-binaries.sh'

$ build-vadr-binaries.sh 

=============================
SETTING ENVIRONMENT VARIABLES 
=============================

Before you can use VADR you need to set specific environment
variables. 

----------------------------------------------
Setting environment for current shell session:
----------------------------------------------
This can be done for the current shell with the script
'set-vadr-environment.sh'. You need to supply the path to the current
directory (where you performed steps 1 and 2).

First, determine current directory:

$ pwd
<path-to-current-directory>

Then run 'set-vadr-environment.sh' with that path, like:

$ set-vadr-environment.sh <path-to-current-directory>.

-------------------------------------------
Setting environment for all shell sessions: 
-------------------------------------------

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


