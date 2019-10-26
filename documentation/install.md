# VADR installation instructions

* [Installation using `vadr-install.sh`](#install.sh)
* [Setting environment variables](#environment)
* [Verifying successful installation](#tests)
* [Further information](#further)

---
## VADR installation using the `vadr-install.sh` script

The file `vadr-install.sh` is an executable file for installing VADR
and its dependencies. That file is located online at github.
To install the latest version of VADR download this file:

https://raw.githubusercontent.com/nawrockie/vadr/master/vadr-install.sh

To download any specific release/version, for example version 1.0 download
the corresponding `vadr-install.sh` file for that release/version
(prior to version 1.0, the name of the installation script was
`install.sh`, not `vadr-install.sh`):

https://raw.githubusercontent.com/nawrockie/vadr/1.0/vadr-install.sh

Copy the `vadr-install.sh` file into the directory in which you want
to install VADR. A good name for that directory is
`vadr-install-dir`. Then move into that directory and run one of the
following two commands depending on whether you are installing on a
Linux or Mac/OSX system:

```
./vadr-install.sh linux
```

OR

```
./vadr-install.sh macosx
```
The `linux` or `macosx` argument controls (only) the type of blast
executable files that will be installed.

The `vadr-install.sh` command will create several directories in the current directory.
It will download and install the software packages
VADR, and [Infernal](http://eddylab.org/infernal/), the required perl
module libraries [sequip](https://github.com/nawrockie/sequip),
[Bio-Easel](https://github.com/nawrockie/Bio-Easel), as well as the
binary executables of the NCBI BLAST package (for either Linux or
Mac/OSX).

When `vadr-install.sh` is finished running it will print important
instructions to the screen that explain how to modify your environment
variables so that you can run the VADR scripts, as discussed next.

If you have trouble with installation, email eric.nawrocki@nih.gov.

---
## Setting VADR environment variables <a name="environment"></a>

As mentioned above, when you run `vadr-install.sh`, instructions will be
output about how to change your environment variables so that you can
run the VADR scripts. Those instructions are also included below for
reference, but without the actual path to where you ran `vadr-install.sh`
(below it is replaced with `<full path to directory in which you ran
vadr-install.sh>`)

---
### **Note to internal NCBI users**
Contact Eric Nawrocki (eric.nawrock@nih.gov) to find out the
path to the centrally installed copy of VADR at NCBI (`<ncbi-vadr-install-dir>`)

Then, to set up your environment variables follow the instructions
below but replace:
`<full path to directory in which you ran vadr-install.sh>`
with
`<ncbi-vadr-install-dir>`
---

### Instructions for setting environment variables output by `vadr-install.sh`
```
The final step is to update your environment variables.
(See vadr/README.txt for more information.)

If you are using the bash shell, add the following
lines to the '.bashrc' file in your home directory:

export VADRINSTALLDIR=<full path to directory in which you ran install.sh>
export VADRSCRIPTSDIR="$VADRINSTALLDIR/vadr"
export VADRMODELDIR="$VADRINSTALLDIR/vadr-models"
export VADRINFERNALDIR="$VADRINSTALLDIR/infernal-dev/src"
export VADREASELDIR="$VADRINSTALLDIR/infernal-dev/easel/miniapps"
export VADRBIOEASELDIR="$VADRINSTALLDIR/Bio-Easel"
export VADRSEQUIPDIR="$VADRINSTALLDIR/sequip"
export VADRBLASTDIR="$VADRINSTALLDIR/ncbi-blast/bin"
export PERL5LIB="$VADRSCRIPTSDIR":"$VADRSEQUIPDIR":"$VADRBIOEASELDIR/blib/lib":"$VADRBIOEASELDIR/blib/arch":"$PERL5LIB"
export PATH="$VADRSCRIPTSDIR":"$PATH"

After adding the export lines to your .bashrc file, source that file
to update your current environment with the command:

source ~/.bashrc

---
If you are using the C shell, add the following
lines to the end of your '.cshrc' file in your home
directory:

setenv VADRINSTALLDIR <full path to directory in which you ran install.sh>
setenv VADRSCRIPTSDIR "$VADRINSTALLDIR/vadr"
setenv VADRMODELDIR "$VADRINSTALLDIR/vadr-models"
setenv VADRINFERNALDIR "$VADRINSTALLDIR/infernal-dev/src"
setenv VADREASELDIR "$VADRINSTALLDIR/infernal-dev/easel/miniapps"
setenv VADRBIOEASELDIR "$VADRINSTALLDIR/Bio-Easel"
setenv VADRSEQUIPDIR "$VADRINSTALLDIR/sequip"
setenv VADRBLASTDIR "$VADRINSTALLDIR/ncbi-blast/bin"
setenv PERL5LIB "$VADRSCRIPTSDIR":"$VADRSEQUIPDIR":"$VADRBIOEASELDIR/blib/lib":"$VADRBIOEASELDIR/blib/arch":"$PERL5LIB"
setenv PATH "$VADRSCRIPTSDIR":"$PATH"

After adding the setenv lines to your .cshrc file, source that file
to update your current environment with the command:

source ~/.cshrc

(To determine which shell you use, type: 'echo $SHELL')

```
---

### If you get an error about `PERL5LIB` being undefined, change the PERL5LIB
line to add to:

```
export PERL5LIB="$VADRSCRIPTSDIR":"$VADRSEQUIPDIR":"$VADRBIOEASELDIR/blib/lib":"$VADRBIOEASELDIR/blib/arch"
````

for `.bashrc`, *or*

```
setenv PERL5LIB "$VADRSCRIPTSDIR":"$VADRSEQUIPDIR":"$VADRBIOEASELDIR/blib/lib":"$VADRBIOEASELDIR/blib/arch"
```

for `.cshrc`. And then execute `source ~/.bashrc` or `source ~/.cshrc` again.

---
## Verifying successful installation with test runs

The VADR package includes some tests you can run to make sure that
your installation was successful and that your environment variables
are set-up correctly.

There are 2 shell scripts for running tests; with respect to the
installation directory they are:

1. `vadr/testfiles/do-install-tests-local.sh`
2. `vadr/testfiles/do-install-tests-parallel.sh`

As explained [here](annotate.md#exampleparallel), the VADR scripts can
be run locally on your computer or in parallel on a compute
farm. These two test files test each of those modes.  If you plan to
run the scripts locally at least some of the time, then run
`do-install-tests-local.sh`. If you plan to run the scripts on a
compute farm at least some of the time, then run
`do-install-tests-parallel.sh`.

You can run the scripts like this:

```
$VADRSCRIPTSDIR/testfiles/do-install-tests-local.sh
```
and 
```
$VADRSCRIPTSDIR/testfiles/do-install-tests-parallel.sh
```

These scripts can take up to several minutes to run. Please be patient.
If something goes wrong, the `local` script will exit quickly. If the 
compute farm is busy, the `parallel` script make take longer as the
relevant jobs wait to run.

Below is example output for do-install-tests-local.sh:

```
# v-test.pl :: test VADR scripts [TEST SCRIPT]
# VADR 0.991 (Aug 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Sat Oct 26 16:28:18 2019
#
# test file:                    /usr/local/vadr-install-dir/vadr/testfiles/noro.r10.local.testin
# output directory:             n10-local
# forcing directory overwrite:  yes [-f]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parsing test file                                  ... done. [    0.0 seconds]
# Running command  1 [annotate-noro-10-local]        ... done. [   45.2 seconds]
#	checking va-noro.r10/va-noro.r10.vadr.pass.tbl                                                                ... pass
#	checking va-noro.r10/va-noro.r10.vadr.fail.tbl                                                                ... pass
#	checking va-noro.r10/va-noro.r10.vadr.sqa                                                                     ... pass
#	checking va-noro.r10/va-noro.r10.vadr.sqc                                                                     ... pass
#	checking va-noro.r10/va-noro.r10.vadr.ftr                                                                     ... pass
#	checking va-noro.r10/va-noro.r10.vadr.sgm                                                                     ... pass
#	checking va-noro.r10/va-noro.r10.vadr.mdl                                                                     ... pass
#	checking va-noro.r10/va-noro.r10.vadr.alt                                                                     ... pass
#	checking va-noro.r10/va-noro.r10.vadr.alc                                                                     ... pass
#	removing directory va-noro.r10                               ... done
#
#
# PASS: all 9 files were created correctly.
#
#
# Output printed to screen saved in:                   n10-local.vadr.log
# List of executed commands saved in:                  n10-local.vadr.cmd
# List and description of all output files saved in:   n10-local.vadr.list
#
# All output files created in directory ./n10-local/
#
# Elapsed time:  00:00:45.73
#                hh:mm:ss
# 
[ok]
# v-test.pl :: test VADR scripts [TEST SCRIPT]
# VADR 0.991 (Aug 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Sat Oct 26 16:29:05 2019
#
# test file:                    /usr/local/vadr-install-dir/vadr/testfiles/dengue.r5.local.testin
# output directory:             d5-local
# forcing directory overwrite:  yes [-f]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parsing test file                                  ... done. [    0.0 seconds]
# Running command  1 [annotate-dengue-5-local]       ... done. [   44.2 seconds]
#	checking va-dengue.r5/va-dengue.r5.vadr.pass.tbl                                                              ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.fail.tbl                                                              ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.sqa                                                                   ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.sqc                                                                   ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.ftr                                                                   ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.sgm                                                                   ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.mdl                                                                   ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.alt                                                                   ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.alc                                                                   ... pass
#	removing directory va-dengue.r5                              ... done
#
#
# PASS: all 9 files were created correctly.
#
#
# Output printed to screen saved in:                   d5-local.vadr.log
# List of executed commands saved in:                  d5-local.vadr.cmd
# List and description of all output files saved in:   d5-local.vadr.list
#
# All output files created in directory ./d5-local/
#
# Elapsed time:  00:00:44.97
#                hh:mm:ss
# 
[ok]
```

The two most important lines are the lines that begins with `# PASS`

```
PASS: all 9 files were created correctly.
PASS: all 9 files were created correctly.
```

This means that the test has passed. You should see similar 
lines when you run the other tests. If you do not and need help
figuring out why, email me at eric.nawrocki@nih.gov.

---
## Further information

* [`v-annotate.pl` example usage and command-line options](annotate.md)
* [`v-build.pl` example usage and command-line options](build.md)
* [VADR output formats](formats.md)

---
#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.


