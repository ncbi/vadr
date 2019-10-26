EPN, Wed Oct 16 10:45:13 2019

Instructions for installing VADR in multiple steps
--------------------------------------------------

Choose a directory to install in. All source code, executables, and
VADR models (several directories) will be placed there. 

Move to that directory and copy the following files there:
vadr/vadr-install.sh
vadr/testfiles/download-vadr-cox1-models.sh

To install on a Linux system (to install on Mac/OSX replace 'linux'
below with 'macosx'):

Step 1. Download packages with command:

$ ./vadr-install.sh linux download

Step 2. Build binaries with the command below. 

!!! Make sure you save or copy the end of the resulting output,
specifically the section beginning with the text 'The final step is to
update your environment variables.' to the end of the output. You will
use part of this text to set your environment variables later.

$ ./vadr-install.sh linux build

Step 3. Download VADR cox1 models (if desired) with the script
'download-vadr-cox1-models.sh' with 1 command line argument: <vadr
cox1 model version (e.g. 0.971.1)>:

$ ./download-vadr-cox1-models.sh 0.971.1

Follow the instructions here:
https://github.com/nawrockie/vadr/blob/master/documentation/install.md#environment

To set your environment variables.

Then, follow the instructions here:
https://github.com/nawrockie/vadr/blob/master/documentation/install.md#tests

To verify that installation was successful. For this step it is only
necessary to run `do-install-tests-local.sh`, you can skip the
parallel tests.

--
For help, contact: eric.nawrocki@nih.gov
--

