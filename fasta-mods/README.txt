EPN, Fri Jun 24 15:57:08 2022

How these files were created:

These are from fasta v36.3.8h_04-May-2020

vadr-fasta-defs.h:
- is manually modified version of fasta's src/defs.h

vadr-fasta-defs.patch: 
- created with command 'diff -u defs.h vadr-fasta-defs.h > vadr-fasta-defs.patch'

vadr-fasta-Makefile.linux: 
- is manually modified version of fasta's make/Makefile.linux

vadr-fasta-Makefile.patch: 
- created with command 'diff -u Makefile.linux vadr-fasta-Makefile.linux > vadr-fasta-Makefile.linux.patch'

vadr-fasta-Makefile.os_x86_64: 
- is manually modified version of fasta's make/Makefile.os_x86_64

vadr-fasta-Makefile.patch: 
- created with command 'diff -u Makefile.os_x86_64 vadr-fasta-Makefile.os_x86_64 > vadr-fasta-Makefile.os_x86_64.patch'

----------

vadr-fasta.patch:
- this file is used internally by NCBI in TeamCity to patch both
  fasta/src/defs.h and fasta/make/Makefile.linux at the same time
  using `patch -p1 < vadr-fasta.patch` from the 'fasta/' directory
  after downloading v36.3.8h_04-May-2020.

  The following steps show how this file was created (mostly logged
  here in case it needs to be made again in the future):
  1. Make a new directory called vadr-fasta-patch, and initialize it
     as a git repo:
     > mkdir vadr-fasta-patch
     > cd vadr-fasta-patch
     > git init .

  2. download fasta with command:
     > curl -k -L -o v36.3.8h_04-May-2020.zip https://github.com/wrpearson/fasta36/archive/v36.3.8h_04-May-2020.zip; unzip v36.3.8h_04-May-2020.zip; mv fasta36-36.3.8h_04-May-2020 fasta; rm v36.3.8h_04-May-2020.zip

  3. Make new directories called src and make, and copy fasta files
     into them:
     > mkdir src 
     > mkdir make
     > cp fasta/src/defs.h src/
     > cp fasta/make/Makefile.linux make/
     > git add src
     > git add make     
     > git commit -m "Adds original files"

  4. Copy modified versions of src/defs.h and make/Makefile.linux
     from vadr git repo over the original verions:
     > cp PATH-TO-VADR-REPO/fasta-mods/vadr-fasta-defs.h src/defs.h
     > cp PATH-TO-VADR-REPO/fasta-mods/vadr-fasta-Makefile.linux make/Makefile.linux 

  5. Use 'git diff' to make a patch file:
     > git diff > vadr-fasta.patch

  Then as an example of how to apply this patch, move into fasta and
  use the patch command: 
  > cd fasta
  > patch -p1 < ../vadr-fasta.patch 

  And to test if it worked, make sure the newly patched files are
  identical to their counterparts in the vadr git repo:
  > diff make/Makefile.linux PATH-TO-VADR-REPO/fasta-mods/vadr-fasta-Makefile.linux
  > diff src/defs.h PATH-TO-VADR-REPO/fasta-mods/vadr-fasta-defs.h 
  > 

You can now delete the 'vadr-fasta-patch' directory you created in
step 1.



