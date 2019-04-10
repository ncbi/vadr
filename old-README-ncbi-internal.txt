EPN, Fri Feb  8 10:42:31 2019

dnaorg_scripts 0.45 README-ncbi-internal.txt 

Organization of this file:

INTRODUCTION
SETTING UP ENVIRONMENT VARIABLES
VERIFYING YOU CAN RUN DNAORG SCRIPTS
USAGE AND OPTIONS OF DNAORG SCRIPTS
EXAMPLE RUNS FOR INTERNAL NCBI USERS
OPTIONS THAT MUST BE USED CONSISTENTLY IN BOTH dnaorg_build.pl AND dnaorg_annotate.pl

Questions:
email Eric Nawrocki: eric.nawrocki@nih.gov

See also README.txt which has additional information including a list
of possible errors that can be output in feature tables.

##############################################################################
INTRODUCTION

dnaorg_scripts version 0.45 is installed in system-wide directories
for internal use at NCBI. The top-level directory is:

/panfs/pan1/dnaorg/virseqannot/code/dnaorg_scripts

And other packages/modules that dnaorg_scripts requires to run are
installed here:

/panfs/pan1/dnaorg/virseqannot/code/Bio-Easel
/panfs/pan1/dnaorg/virseqannot/code/epn-options
/panfs/pan1/dnaorg/virseqannot/code/esl-epn-translate
/panfs/pan1/dnaorg/virseqannot/code/esl-fetch-cds
/panfs/pan1/dnaorg/virseqannot/code/hmmer-3.1b2
/panfs/pan1/dnaorg/virseqannot/code/infernal-dev
/panfs/pan1/dnaorg/virseqannot/code/dnaorg-build-directories

You will need to update your environment variables as described below
in order to run the dnaorg scripts. Additionally if you want to run
dnaorg script jobs in parallel on the Sun Grid Engine (SGE) compute farm 
you will need access to the farm. See the SGE quick start confluence
page for help on this:
https://confluence.ncbi.nlm.nih.gov/pages/viewpage.action?spaceKey=UGE&title=Grid+Engine+Quick+Start

The README.txt for dnaorg_scripts includes an example run of using
dnaorg_classify.pl to annotate Norovirus sequences. If that is what
you intend to use dnaorg scripts for, that example is more relevant 
than the examples used here. That README.txt is here:

/panfs/pan1/dnaorg/virseqannot/code/dnaorg_scripts/README.txt

Git repository for dnaorg_scripts:
https://github.com/nawrockie/dnaorg_scripts.git

(To be automatically updated about new releases, you can become a
watcher of the repo on github.)

authors: Eric Nawrocki and Alejandro Schaffer
contact: eric.nawrocki@nih.gov

##############################################################################
SETTING UP ENVIRONMENT VARIABLES

Before you can run any dnaorg scripts, you will need to update some of your
environment variables. To do this, add the following lines to
your .bashrc file (if you use bash shell) or .cshrc file (if you use C
shell or tcsh). The .bashrc or .cshrc file is in your home
directory. To determine what shell you use, type
> echo $SHELL
If this command returns '/bin/bash', then update your .bashrc file.
If this command returns'/bin/csh' or '/bin/tcsh' then update your .cshrc file.

The lines to add to your .bashrc file:
-----------
export DNAORGDIR=/panfs/pan1/dnaorg/virseqannot/code
export DNAORGBUILDDIR="$DNAORGDIR"/dnaorg-build-directories
export PERL5LIB="$DNAORGDIR"/dnaorg_scripts:"$DNAORGDIR"/epn-options:"$DNAORGDIR"/Bio-Easel/blib/lib:"$DNAORGDIR"/Bio-Easel/blib/arch:"$PERL5LIB"
export PATH="$DNAORGDIR"/dnaorg_scripts:"$PATH"
-------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The analogous lines to add to your .cshrc file:
-----------
setenv DNAORGDIR /panfs/pan1/dnaorg/virseqannot/code
setenv DNAORGBUILDDIR "$DNAORGDIR"/dnaorg-build-directories
setenv PERL5LIB "$DNAORGDIR"/dnaorg_scripts:"$DNAORGDIR"/epn-options:"$DNAORGDIR"/Bio-Easel/blib/lib:"$DNAORGDIR"/Bio-Easel/blib/arch:"$PERL5LIB"
setenv PATH "$DNAORGDIR"/dnaorg_scripts:"$PATH"
-----------

After adding the lines specified above, execute the command:
> source ~/.bashrc
or
> source ~/.cshrc

If you get an error about PERL5LIB being undefined, change the PERL5LIB
line to add to:
export PERL5LIB="$DNAORGDIR"/dnaorg_scripts:"$DNAORGDIR"/epn-options:"$DNAORGDIR"/Bio-Easel/blib/lib:"$DNAORGDIR"/Bio-Easel/blib/arch
for .bashrc, OR
setenv PERL5LIB "$DNAORGDIR"/dnaorg_scripts:"$DNAORGDIR"/epn-options:"$DNAORGDIR"/Bio-Easel/blib/lib:"$DNAORGDIR"/Bio-Easel/blib/arch
for .cshrc. And then do
> source ~/.bashrc
or
> source ~/.cshrc
again.

###########################################################################
VERIFYING YOU CAN RUN DNAORG SCRIPTS

The dnaorg package include some tests you can run to make sure that
your environment variables are set-up correctly.

There are 2 shell scripts for running tests:

     1	do-install-tests-local.sh
     2	do-install-tests-parallel.sh

The dnaorg scripts can be run locally on your computer or in parallel
on a compute farm. These two test files test each of those modes.  To
make sure your environment is set-up correctly and all the dnaorg
related programs are working as intended, you should run both of these
and ensure they finish successfully as explained below.

These scripts can take up to several minutes to run. Please be patient.
If something goes wrong the 'local' script will exit quickly. If the 
compute farm is busy, the 'parallel' script make take longer as the
relevant jobs wait to run.

Below is example output for do-install-tests-local.sh:
> sh $DNAORGDIR/testfiles/do-install-tests-local.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_test.pl :: test dnaorg scripts [TEST SCRIPT]
# dnaorg 0.45 (Feb 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Feb  8 10:38:51 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0205_virus_dengue_linda_premature_stop
#
# test file:                                                   /panfs/pan1/infernal/notebook/19_0205_virus_dengue_linda_premature_stop/dnaorg_scripts/testfiles/noro.r10.local.testin
# forcing directory overwrite:                                 yes [-f]
# build directory, replaces !dirbuild! in test file with <s>:  /panfs/pan1/infernal/notebook/19_0205_virus_dengue_linda_premature_stop/dnaorg-build-directories/norovirus-builds [--dirbuild]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Running command  1 [classify-noro-10-local]        ... done. [74.7 seconds]
#	checking dt-noro.r10.l/dt-noro.r10.l-NC_029646/dt-noro.r10.l-NC_029646.dnaorg_annotate.ap.sqtable             ... pass
#	checking dt-noro.r10.l/dt-noro.r10.l-NC_029646/dt-noro.r10.l-NC_029646.dnaorg_annotate.af.sqtable             ... pass
#	checking dt-noro.r10.l/dt-noro.r10.l-NC_039476/dt-noro.r10.l-NC_039476.dnaorg_annotate.ap.sqtable             ... pass
#	checking dt-noro.r10.l/dt-noro.r10.l-NC_039476/dt-noro.r10.l-NC_039476.dnaorg_annotate.af.sqtable             ... pass
#	checking dt-noro.r10.l/dt-noro.r10.l-NC_039477/dt-noro.r10.l-NC_039477.dnaorg_annotate.ap.sqtable             ... pass
#	checking dt-noro.r10.l/dt-noro.r10.l-NC_039477/dt-noro.r10.l-NC_039477.dnaorg_annotate.af.sqtable             ... pass
#	removing directory dt-noro.r10.l                             ... done
#
#
# PASS: all 6 files were created correctly.
#
#
# Output printed to screen saved in:                   n10-local.dnaorg_test.log
# List of executed commands saved in:                  n10-local.dnaorg_test.cmd
# List and description of all output files saved in:   n10-local.dnaorg_test.list
#
# All output files created in directory ./n10-local/
#
# CPU time:  00:01:16.79
#            hh:mm:ss
# 
# DNAORG-SUCCESS
# dnaorg_test.pl :: test dnaorg scripts [TEST SCRIPT]
# dnaorg 0.45 (Feb 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Feb  8 10:40:09 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0205_virus_dengue_linda_premature_stop
#
# test file:                                                   /panfs/pan1/infernal/notebook/19_0205_virus_dengue_linda_premature_stop/dnaorg_scripts/testfiles/dengue.r5.local.testin
# forcing directory overwrite:                                 yes [-f]
# build directory, replaces !dirbuild! in test file with <s>:  /panfs/pan1/infernal/notebook/19_0205_virus_dengue_linda_premature_stop/dnaorg-build-directories/dengue-builds [--dirbuild]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Running command  1 [classify-dengue-5-local]       ... done. [51.1 seconds]
#	checking dt-dengue.r5.l/dt-dengue.r5.l-NC_001474/dt-dengue.r5.l-NC_001474.dnaorg_annotate.ap.sqtable          ... pass
#	checking dt-dengue.r5.l/dt-dengue.r5.l-NC_001474/dt-dengue.r5.l-NC_001474.dnaorg_annotate.af.sqtable          ... pass
#	checking dt-dengue.r5.l/dt-dengue.r5.l-NC_001477/dt-dengue.r5.l-NC_001477.dnaorg_annotate.ap.sqtable          ... pass
#	checking dt-dengue.r5.l/dt-dengue.r5.l-NC_001477/dt-dengue.r5.l-NC_001477.dnaorg_annotate.af.sqtable          ... pass
#	removing directory dt-dengue.r5.l                            ... done
#
#
# PASS: all 4 files were created correctly.
#
#
# Output printed to screen saved in:                   d5-local.dnaorg_test.log
# List of executed commands saved in:                  d5-local.dnaorg_test.cmd
# List and description of all output files saved in:   d5-local.dnaorg_test.list
#
# All output files created in directory ./d5-local/
#
# CPU time:  00:00:51.60
#            hh:mm:ss
# 
# DNAORG-SUCCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The two most important lines are the lines that begins with "# PASS"

# PASS: all 6 files were created correctly.
# PASS: all 4 files were created correctly.

This means that the test has passed. You should see similar 
lines when you run the other tests. If you do not and need help
figuring out why, email me at eric.nawrocki@nih.gov.

##############################################################################
USAGE AND OPTIONS OF DNAORG SCRIPTS

The output of each program run with the '-h' option
is informative about how to use that script:

> dnaorg_classify.pl -h
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_classify.pl :: classify sequences using an HMM library of RefSeqs
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Thu Jan 24 21:22:17 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
Usage: This script must be run in 1 of 2 modes:

Build mode (--onlybuild): Build HMM library and exit (no classification).
Example usage:
	dnaorg_classify.pl [-options] --onlybuild <RefSeq list> --dirout <output directory to create with HMM library>

Classify mode (--infasta): Use a previously created HMM library to annotate sequences in an input fasta file.
Example usage:
	dnaorg_classify.pl [-options] --dirbuild <directory with HMM library to use> --dirout <output directory to create> --infasta <fasta file with sequences to classify>

REQUIRED options:
  --dirout <s> : REQUIRED: name for output directory to create is <s>

options for selecting mode to run in: build-mode (--onlybuild) or classify-mode:
  --onlybuild <s> : build an HMM library for sequences listed in <s>, then exit
  --infasta <s>   : fasta file with sequences to classify is <s>, requires --dirbuild
  --dirbuild <s>  : specify directory with HMM library is <s>

basic options:
  -f             : force; if dir <reference accession> exists, overwrite it
  -v             : be verbose; output commands to stdout as they're run
  --keep         : do not remove intermediate files, keep them all on disk
  --nkb <n>      : set target number of KB of sequences for each nhmscan farm job to <n> [5]
  --maxnjobs <n> : set max number of jobs to submit to compute farm to <n> [500]
  --wait <n>     : allow <n> wall-clock minutes for nhmmscan jobs on farm to finish, including queueing time [500]
  --local        : run nhmmscan locally instead of on farm
  --errcheck     : consider any farm stderr output as indicating a job failure

options for controlling what unexpected features cause sequences to PASS/FAIL:
  --lowcovpass   : sequences with low coverage can PASS
  --unexppass    : sequences with unexpected classification can PASS
  --allfail      : sequences with >=1 unexpected feature(s) FAIL
  --lowscfail    : sequences with LowScore     unexpected feature FAIL
  --vlowscfail   : sequences with VeryLowScore unexpected feature FAIL
  --lowdifffail  : sequences with LowDiff      unexpected feature FAIL
  --vlowdifffail : sequences with VeryLowDiff  unexpected feature FAIL
  --biasfail     : sequences with HighBias     unexpected feature FAIL
  --minusfail    : sequences with MinusStrand  unexpected feature FAIL

options for controlling reporting of unexpected features:
  --lowcovthresh <x>   : fractional coverage threshold for Low Coverage is <x> [0.9]
  --lowscthresh <x>    : bits per nucleotide threshold for LowScore unexpected feature is <x> [0.3]
  --vlowscthresh <x>   : bits per nucleotide threshold for VeryLowScore unexpected feature is <x> [0.2]
  --lowdiffthresh <x>  : bits per nucleotide diff threshold for LowDiff unexpected feature is <x> [0.06]
  --vlowdiffthresh <x> : bits per nucleotide diff threshold for VeryLowDiff unexpected feature is <x> [0.006]
  --biasfract <x>      : fractional threshold for HighBias unexpected feature is <x> [0.25]
  --lowscminlen <n>    : set minimum length for which LowScore causes FAILure to <n> [501]
  --lowdiffminlen <n>  : set minimum length for which LowDiff causes FAILure to <n> [1001]
  --nolowscminlen      : no minimum length for which LowScore causes a seq to FAIL
  --nolowdiffminlen    : no minimum length for which LowDiff causes a seq to FAIL

options for automatically running dnaorg_annotate.pl for classified sequences:
  -A <s>         : annotate using dnaorg_build.pl build directories in <s> after classifying
  --optsA <s>    : read additional dnaorg_annotate.pl options from file <s>
  --reflistA <s> : only annotate seqs that match to RefSeqs listed in <s>

in combination with -A, options for tuning protein validation with blastx (don't list these in --optsA <f> file):
  --xalntol <n>    : max allowed difference in nucleotides b/t nucleotide and blastx start/end postions is <n> [5]
  --xindeltol <n>  : max allowed nucleotide insertion and deletion length in blastx validation is <n> [27]
  --xlonescore <n> : minimum score for a lone blastx (not supported by a CM hit) to cause an error is <n> [80]

options for defining expected classifications:
  --ecmap <s>    : read map of model names to taxonomic names from file <s>
  --ecall <s>    : set expected classification of all seqs to model/tax <s>
  --eceach <s>   : read expected classification for each sequence from file <s>
  --ecthresh <x> : expected classification must be within <x> bits/nt of top match [0.3]
  --ectoponly    : top match must be expected classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As reported by the dnaorg_classify.pl -h output, dnaorg_classify.pl
must be run in 1 of 2 modes. The first time you run it, it must be run
in the --onlybuild mode, which builds an HMM library of the RefSeqs
that you are going to use to classify input sequences in subsequent
runs of the program in --infasta modes.  An example is given above in
"Example 1".

############

After you run dnaorg_classify.pl in the --infasta mode, you will have
a list of sequences that have been classified as most similar to each
of the RefSeqs listed in the input list of RefSeqs you specified when
you ran dnaorg_classify.pl with the --onlybuild option. You can then
use the dnaorg_build.pl and dnaorg_annotate.pl scripts for each of
those RefSeqs to annotate the sequences that have been classified to
each RefSeq.

Another way you can use dnaorg_classify.pl is to have it automatically
call dnaorg_annotate.pl for all of the sequences it classifies using
the appropriate RefSeq models. In order to do this, you must have
already ran dnaorg_build.pl for all of the possible RefSeqs that
dnaorg_classify.pl may classify a sequence to. An example of this is
example 3 below.

Here are the options and usage for dnaorg_build.pl:

> dnaorg_build.pl -h
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_build.pl :: build homology models for features of a reference sequence
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 06:17:10 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
Usage: dnaorg_build.pl [-options] <reference accession>

basic options:
  -c            : genome is circular
  -f            : force; if dir <reference accession> exists, overwrite it
  -n <n>        : for non-big models, set number of CPUs for calibration to <n> [4]
  -v            : be verbose; output commands to stdout as they're run
  --dirout <s>  : specify output directory as <s>, not <ref accession>
  --matpept <s> : read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships
  --nomatpept   : ignore mat_peptide information in reference annotation
  --xfeat <s>   : build models of additional qualifiers in string <s>
  --dfeat <s>   : annotate qualifiers in <s> from duplicates (e.g. gene from CDS)
  --keep        : do not remove intermediate files, keep them all on disk

options affecting calibration of models:
  --slow          : use default cmcalibrate parameters, not parameters optimized for speed
  --local         : run cmcalibrate locally, do not submit calibration jobs for each CM to the compute farm
  --wait <n>      : allow <n> wall-clock minutes for cmcalibrate jobs on farm to finish, including queueing time [1800]
  --nosubmit      : do not submit cmcalibrate jobs to farm, run later with qsub script
  --errcheck      : consider any stderr output as indicating a job failure
  --rammult       : for all models, multiply RAM Gb by ncpu for mem_free
  --bigthresh <n> : set minimum length for a big model to <n> [2500]
  --bigram <n>    : for big models, set Gb RAM per core for calibration to <n> [8]
  --biglen <x>    : for big models, set cmcalibrate length to search in Mb as <x> [0.16]
  --bigncpu <n>   : for big models, set number of CPUs for calibration to <n> [4]
  --bigtailp <x>  : for big models, set --tailp cmcalibrate parameter as <x> [0.30]

optional output files:
  --mdlinfo : create file with internal model information
  --ftrinfo : create file with internal feature information

options for skipping stages and using files from an earlier, identical run, primarily useful for debugging:
  --skipedirect : skip the edirect steps, use data from an earlier run of the script
  --skipfetch   : skip the sequence fetching steps, use files from an earlier run of the script
  --skipbuild   : skip the model building/calibrating, requires --mdlinfo and/or --ftrinfo

options for building models for origin sequences:
  --orginput <s> : read training alignment for origin sequences from file <s>
  --orgstart <n> : origin sequence starts at position <n> in file <s> from --orginput <s> [0]
  --orglen <n>   : origin sequence is <n> nucleotides long [0]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dnaorg_build.pl takes a single argument: the reference 
accession. It then builds models based for sequence features in the
reference sequence using the compute farm. 

dnaorg_annotate.pl should not be run until after dnaorg_build.pl is
finished.

dnaorg_annotate.pl is typically run with a single command line
argument, a .ntlist file that contains a list of accessions, one per
line. 

The accessions in the .ntlist file should be the same species
as the reference accession used in dnaorg_build.pl. The reference
accession used passed into dnaorg_build.pl should be the first
accession listed in the .ntlist file.

See examples 3 and 4 below for more detailed instructions.

Alternatively, with the --infasta option, the user can provide a
fasta file of sequences to annotate instead of a list of accessions
as the .ntlist file. This may be desirable because in this mode,
none of the NCBI databases are accessed (for sequence fetching,
etc.). In particular, the --infasta option makes it feasible to
annotate newly arriving sequences that are not in GenBank. See
example 5 below for an example of running dnaorg_annotate.pl in this
mode. 

When the dnaorg_classify.pl script is run with the --infasta
option, it will create fasta files for each of the RefSeqs that had
at least one sequence assigned to it. These fasta files can be 
used as input to dnaorg_annotate.pl with the --infasta option.

Here is the usage and option information for dnaorg_annotate.pl:

> dnaorg_annotate.pl -h
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_annotate.pl :: annotate sequences based on a reference annotation
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 06:19:13 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
Usage: dnaorg_annotate.pl [-options] <file with list of accessions to annotate>

       OR

       dnaorg_annotate.pl [-options] --refaccn <reference accession> --infasta <fasta sequence file with sequences to annotate>

basic options:
  -c              : genome is closed (a.k.a circular)
  -f              : force; if dir from --dirout exists, overwrite it
  -v              : be verbose; output commands to stdout as they're run
  --dirout <s>    : specify output directory as <s>, not <ref accession>
  --dirbuild <s>  : specify output directory used for dnaorg_build.pl as <s> (created with dnaorg_build.pl --dirout <s>), not <ref accession>
  --origin <s>    : identify origin seq <s> in genomes, put "|" at site of origin ("|" must be escaped, i.e. "\|"
  --matpept <s>   : read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships
  --nomatpept     : ignore mat_peptide information in reference annotation
  --xfeat <s>     : use models of additional qualifiers in string <s>
  --dfeat <s>     : annotate qualifiers in <s> from duplicates (e.g. gene from CDS)
  --specstart <s> : read specified alternate start codons per CDS from file <s>
  --keep          : do not remove intermediate files, keep them all on disk

options for alternative modes:
  --infasta     : single cmdline argument is a fasta file of sequences, not a list of accessions
  --refaccn <s> : specify reference accession is <s> (must be used in combination with --infasta)

options for tuning protein validation with blastx:
  --xalntol <n>    : max allowed difference in nucleotides b/t nucleotide and blastx start/end postions is <n> [5]
  --xindeltol <n>  : max allowed nucleotide insertion and deletion length in blastx validation is <n> [27]
  --xlonescore <n> : minimum score for a lone blastx (not supported by a CM hit) to cause an error is <n> [80]

options for modifying which errors are reported:
  --classerrors <s> : read per-sequence classification errors from <s> [0]
  --allolp          : report all olp errors, even when other feature is not predicted (nop error)
  --alladj          : report all aja/ajb errors, even when other feature is not predicted (nop error)

options for changing search sensitivity modes:
  --midthresh <n>   : set max model length for using mid sensitivity mode to <n> [75]
  --smallthresh <n> : set max model length for using max sensitivity mode to <n> [30]
  --hmmonly <n>     : run in HMM-only mode for models with >= <n> positions [0]

options related to parallelization on compute farm:
  --local        : run cmscan locally instead of on farm
  --errcheck     : consider any farm stderr output as indicating a job failure
  --nkb <n>      : set target number of KB of sequences for each cmscan farm job to <n> [50]
  --maxnjobs <n> : set max number of jobs to submit to compute farm to <n> [2500]
  --wait <n>     : allow <n> wall-clock minutes for cmscan jobs on farm to finish, including queueing time [500]

options for skipping/adding optional stages:
  --doalign    : create nucleotide and protein alignments
  --mxsize <n> : with --doalign, set --mxsize <n> for cmalign to <n> [2048]

options that modify the tabular output file:
  --tblfirst  : include annotation for first accession on each page of .tbl output file
  --tblnocomp : do not include information comparing predicted annotations to existing GenBank annotations

optional output files:
  --mdlinfo : create file with internal model information
  --ftrinfo : create file with internal feature information
  --seqinfo : create file with internal sequence information
  --errinfo : create file with internal error information

options for skipping stages and using files from earlier, identical runs, primarily useful for debugging:
  --skipedirect   : skip the edirect steps, use data from an earlier run of the script
  --skipfetch     : skip the sequence fetching steps, use files from an earlier run of the script
  --skipscan      : skip the cmscan step, use results from an earlier run of the script
  --skiptranslate : skip the translation steps, use results from an earlier run of the script

TEMPORARY options for the alternative method of identifying origin sequences:
  --aorgmodel <s>    : use alternative origin method with origin model in <s>
  --aorgstart <n>    : origin begins at position <n> in --aorgmodel model [0]
  --aorgoffset <n>   : first position of genome sequence is position <n> in origin sequence [0]
  --aorglen <n>      : length of origin sequence is <n> [0]
  --aorgethresh <x>  : E-value threshold for origin detection is <x> [1]
  --aorgppthresh <x> : average PP threshold for origin detection is <x> [0.6]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################################################
EXAMPLE RUNS FOR INTERNAL NCBI USERS

IMPORTANT: Currently, most of the functionality of the dnaorg scripts
is limited to their use internally at NCBI. None of the following
examples will work external to NCBI. See the file
'README.txt' for examples of how users external to NCBI can
use dnaorg scripts to classify and annotate Norovirus sequences.

There are currently three dnaorg scripts and they are typically
meant to be used in the following order:

1) dnaorg_build.pl: build homology models for a reference accession.
   Must be run once per reference accession you want to use as
   a reference for annotation.

2) dnaorg_classify.pl: classify sequences by comparing to a set of
   reference sequences and picking the most similar one. (The first
   time you run this for a set of possible RefSeqs, you must use
   the --onlybuild option which constructs models for the set of
   reference sequences. You only need to do the --onlybuild run
   once per set of reference sequences.) 

3) dnaorg_annotate.pl: use models built from dnaorg_build.pl to
   annotate other sequences of the same species as the reference.
   Can be run many times for each species provided that
   dnaorg_build.pl has been run for that species (and reference).
   You can also have dnaorg_classify.pl call dnaorg_annotate.pl
   automatically when it finishes classifying using the "-A"
   option as show in example 3 below.

The discussion below includes three examples, which together cover
each of the three scripts.

Example 1: classifying sequences as Maize Streak Virus or Dengue 
           using dnaorg_classify.pl. 

Example 2: annotating MSV sequences using dnaorg_build.pl and
           dnaorg_annotate.pl. 

Example 3. Classifying and annotating Norovirus sequences from
           different Norovirus genotypes using dnaorg_classify.pl.                                  
           (This is nearly identical to the example in README.txt)

The Maize streak virus, which has a circular ssDNA genome with 4 CDS
that code for proteins, one of which has 2 exons.  The Norovirus virus
has a linear ssRNA genome with 3 CDS that code for proteins, one of
which is cleaved into 6 smaller peptides, called "mature
peptides". Together, using the scripts on these two different species
demonstrate most of the features and versatility of the
dnaorg_classify.pl, dnaorg_build.pl and dnaorg_annotate.pl scripts.

###########################################################
Example 1: dnaorg_classify.pl to classify viral sequences #
###########################################################

===================================================================
Example 1, Step 1 of 2. Run dnaorg_classify.pl in --onlybuild mode
===================================================================

Here is an example of running dnaorg_classify.pl to classify 
sequences as either MSV, Dengue or neither.
 
First, dnaorg_classify.pl must be run in 'build mode' using the
--onlybuild option. From above, in the "USAGE AND OPTIONS" section
above, the script reported how to run itself in build mode when it
was called with only the -h option:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Build mode (--onlybuild): Build HMM library and exit (no classification).
Example usage:
	dnaorg_classify.pl [-options] --onlybuild <RefSeq list> --dirout <output directory to create with HMM library>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For this example we'll use the RefSeq list in the file
testfiles/test.refseq.ntlist. Take a look at that file
cat testfiles/test.refseq.ntlist
NC_001346
NC_001477

NC_001346 is the accession of the Maize streak virus RefSeq and
NC_001477 is the accession of the Dengue virus type 1 RefSeq

When we run dnaorg_classify.pl with this input file it will create
an HMM library with HMM models for each of these two RefSeqs. 
We need to specify an output directory name as well, let's use
'test-build': 

> dnaorg_classify.pl --onlybuild $DNAORGDIR/dnaorg_scripts/testfiles/test.refseq.ntlist --dirout test-build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_classify.pl :: classify sequences using an HMM library of RefSeqs
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 06:23:27 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
# REQUIRED name for output directory to create is <s>:  test-build [--dirout]
# build an HMM library for seqs in <s>, then exit:      /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg_scripts/testfiles/test.refseq.ntlist [--onlybuild]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parsing RefSeq list                                                    ... done. [0.0 seconds]
# Creating RefSeq HMM Library                                            ... done. [5.5 seconds]
#
# Output printed to screen saved in:                                     test-build.dnaorg_classify.log
# List of executed commands saved in:                                    test-build.dnaorg_classify.cmd
# List and description of all output files saved in:                     test-build.dnaorg_classify.list
# List of reference sequences used to build HMMs saved in:               test-build.dnaorg_classify.ref.list
# Fasta file of all RefSeq sequences saved in:                           test-build.dnaorg_classify.ref.fa
# List of RefSeq names from fasta file (may include version) saved in:   test-build.dnaorg_classify.ref.fa.list
# Library of HMMs of RefSeqs saved in:                                   test-build.dnaorg_classify.hmm
#
# All output files created in directory ./test-build/
#
# CPU time:  00:00:05.64
#            hh:mm:ss
# 
# DNAORG-SUCCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The script explains that it first parses the RefSeq list and 
then created the RefSeq HMM library. 

The script then reports that it has created several files. These are
all in the directory test-build/. The most important of these is the
HMM file test-build.dnaorg_classify.hmm. This file will be used by
subsequent runs to classify sequences.

#================================================================
Example 1, Step 2 of 2. Run dnaorg_classify.pl in --infasta mode
#================================================================
#
The next step is to run dnaorg_classify.pl again but in the --infasta
mode with a fasta file of sequences to classify.

The usage for the --infasta mode was shown above in the "USAGE AND
OPTIONS" section when the script was run with -h:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Classify mode (--infasta): Use a previously created HMM library to annotate sequences in an input fasta file.
Example usage:
	dnaorg_classify.pl [-options] --dirbuild <directory with HMM library to use> --dirout <output directory to create> --infasta <fasta file with sequences to classify>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We first need to specify the directory with the HMM library to use
with the --dirbuild option, which is 'test-build' from above.  We
also need to specify a new output directory with the --dirout
option, let's use "test-classify". 
Finally, we need a fasta file to use with --infasta. The
file testfiles/test.toclassify.fa is a fasta file with 8
sequences: 2 Dengue type 1 sequences, 4 Maize streak sequences and 2
Norovirus sequences.

To execute dnaorg_classify.pl:
> dnaorg_classify.pl --dirbuild test-build --dirout test-classify --infasta $DNAORGDIR/dnaorg_scripts/testfiles/test.toclassify.fa 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_classify.pl :: classify sequences using an HMM library of RefSeqs
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 06:35:09 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
# REQUIRED name for output directory to create is <s>:  test-classify [--dirout]
# fasta file with sequences to classify is <s>:         /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg_scripts/testfiles/test.toclassify.fa [--infasta]
# specify directory with HMM library is <s>:            test-build [--dirbuild]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parsing RefSeq list                                                    ... done. [0.0 seconds]
# Processing fasta file to get sequence lengths                          ... done. [0.1 seconds]
# Submitting 3 nhmmscan jobs to the farm                                 ... done. [0.3 seconds]
# Waiting a maximum of 500 minutes for all farm jobs to finish           ... 
#	   1 of    3 jobs finished (0.2 minutes spent waiting)
#	   3 of    3 jobs finished (0.5 minutes spent waiting)
# Creating tabular output file                                           ... done. [0.0 seconds]
# Creating seqlists and other output files                               ... done. [0.2 seconds]
#
# Number of input sequences assigned to each RefSeq:
#
#model           num-PASS    num-FAIL
#------------  ----------  ----------
NC_001346               4           0
NC_001477               2           0
NON-ASSIGNED            0           2
#------------  ----------  ----------
SUM-ASSIGNED            6           0
#
# *** The test-classify/test-classify.dnaorg_classify.dnaorg_classify.<s>.fa files (with the --infasta option) can be
# *** used as input to dnaorg_annotate.pl once you've run 'dnaorg_build.pl <s>'
# *** to create models for RefSeq <s>.
#
#
# Output printed to screen saved in:                                                                                                                                      test-classify.dnaorg_classify.log
# List of executed commands saved in:                                                                                                                                     test-classify.dnaorg_classify.cmd
# List and description of all output files saved in:                                                                                                                      test-classify.dnaorg_classify.list
# Fasta file with sequences to classify (copy of /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg_scripts/testfiles/test.toclassify.fa) saved in:   test-classify.dnaorg_classify.fa
# Per-sequence hit and classification information, human readable saved in:                                                                                               test-classify.dnaorg_classify.rdb.infotbl
# Per-sequence hit and classification information, tab-delimited saved in:                                                                                                test-classify.dnaorg_classify.tab.infotbl
# List of errors (unexpected features that cause failure) for all sequences saved in:                                                                                     test-classify.dnaorg_classify.all.errlist
# List of PASSing seqs for NC_001346 saved in:                                                                                                                            test-classify.dnaorg_classify.NC_001346.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_001346 saved in:                                                                                                            test-classify.dnaorg_classify.NC_001346.cp.fa
# List of all seqs for NC_001346 saved in:                                                                                                                                test-classify.dnaorg_classify.NC_001346.seqlist
# Fasta file with all seqs assigned to NC_001346 saved in:                                                                                                                test-classify.dnaorg_classify.NC_001346.fa
# List of PASSing seqs for NC_001477 saved in:                                                                                                                            test-classify.dnaorg_classify.NC_001477.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_001477 saved in:                                                                                                            test-classify.dnaorg_classify.NC_001477.cp.fa
# List of all seqs for NC_001477 saved in:                                                                                                                                test-classify.dnaorg_classify.NC_001477.seqlist
# Fasta file with all seqs assigned to NC_001477 saved in:                                                                                                                test-classify.dnaorg_classify.NC_001477.fa
# List of sequences not assigned to a RefSeq saved in:                                                                                                                    test-classify.dnaorg_classify.noassign.seqlist
# List of RefSeqs in the HMM library saved in:                                                                                                                            test-classify.dnaorg_classify.all.refseqs
# List of sequences that were sorted into seqlists saved in:                                                                                                              test-classify.dnaorg_classify.all.seqs
#
# All output files created in directory ./test-classify/
#
# CPU time:  00:00:31.81
#            hh:mm:ss
# 
# DNAORG-SUCCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This time, the script performs different steps than it did with the
--onlybuild option. From the output:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parsing RefSeq list                                                    ... done. [0.0 seconds]
# Processing fasta file to get sequence lengths                          ... done. [0.1 seconds]
# Submitting 3 nhmmscan jobs to the farm                                 ... done. [0.3 seconds]
# Waiting a maximum of 500 minutes for all farm jobs to finish           ... 
#	   1 of    3 jobs finished (0.2 minutes spent waiting)
#	   3 of    3 jobs finished (0.5 minutes spent waiting)
# Creating tabular output file                                           ... done. [0.0 seconds]
# Creating seqlists and other output files                               ... done. [0.2 seconds]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, the RefSeq list is parsed and sequences to be classified are
processed. The program 'nhmmscan' is used to compare each HMM in the
HMM library to each sequence. The nhmmscan runs are submitted to the
cluster and the script waits for them to finish. In this case they all
finish within 0.5 minutes. The output of nhmmscan is then processed to
get 'infotbl' files. The format of the 'infotbl' files is explained
below in the "FORMAT OF dnaorg_classify.pl's INFOTBL FILE" section.
Finally, .seqlist files were created, which are explained below.

--
IMPORTANT: If you are outside NCBI and your run finished in error, it
is likely because your system does not have a compute farm or does but
does not use Sun Grid Engine the same way that NCBI does. If this
happens, you will want to use the --local option for
dnaorg_classify.pl and dnaorg_annotate.pl. You can repeat the above
dnaorg_classify.pl --infasta command by adding '--local' to your
command like this:

> dnaorg_classify.pl --local --dirbuild test-build --dirout test-classify --infasta $DNAORGDIR/dnaorg_scripts/testfiles/test.toclassify.fa 

And you should get identical output. The only difference is that the
nhmmscan jobs are run locally instead of submitted to the compute
farm. 
--

The script reports on how many sequences were assigned to each
NC_001346 (4), NC_001477 (2) and how many were not assigned (2). 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Number of input sequences assigned to each RefSeq:
#
#model           num-PASS    num-FAIL
#------------  ----------  ----------
NC_001346               4           0
NC_001477               2           0
NON-ASSIGNED            0           2
#------------  ----------  ----------
SUM-ASSIGNED            6           0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As before, the script outputs a list of the output files that it
creates. Of importance here are:

test-classify.dnaorg_classify.NC_001346.cp.fa
and
test-classify.dnaorg_classify.NC_001477.cp.fa

These are sequence files with the sequences from 
testfiles/test.toclassify.fa that were assigned to NC_001346 and
NC_001477 respectively. Exact copies of these files will be used as
input to dnaorg_annotate.pl in examples 2 and 3 below.

#################################################################
Example 2. Annotating MSV sequences using dnaorg_build.pl and 
           dnaorg_annotate.pl                                 
#################################################################

====================================================
Example 2, Step 1 of 2. Run dnaorg_build.pl for MSV 
====================================================
The RefSeq accession for MSV is NC_001346, so to run the
dnaorg_build.pl script one should run:

> dnaorg_build.pl -c NC_001346
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
<[(testfiles)]> dnaorg_build.pl -c NC_001346
# dnaorg_build.pl :: build homology models for features of a reference sequence
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 06:59:31 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
# reference accession:  NC_001346
# genome is circular:   yes [-c]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Outputting information on options used for future use with dnaorg_annotate.pl    ... done. [0.0 seconds]
# Gathering information on reference using edirect                                 ... done. [3.7 seconds]
# Fetching and processing the reference genome                                     ... done. [1.9 seconds]
# Fetching protein translations of CDS and building BLAST DB                       ... done. [1.6 seconds]
# Submitting jobs to build models to compute farm and waiting for them to finish   ... 
#	   4 of    5 jobs finished (0.2 minutes spent waiting)
#	   5 of    5 jobs finished (0.5 minutes spent waiting)
# done. [30.4 seconds]
# Submitting jobs to calibrate models to compute farm and waiting for them to finish ... 
#	   0 of    5 jobs finished (0.2 minutes spent waiting)
#	   0 of    5 jobs finished (0.5 minutes spent waiting)
#	   0 of    5 jobs finished (0.8 minutes spent waiting)
#	   0 of    5 jobs finished (1.0 minutes spent waiting)
#	   0 of    5 jobs finished (1.2 minutes spent waiting)
#	   0 of    5 jobs finished (1.5 minutes spent waiting)
#	   0 of    5 jobs finished (1.8 minutes spent waiting)
#	   0 of    5 jobs finished (2.0 minutes spent waiting)
#	   1 of    5 jobs finished (2.5 minutes spent waiting)
#	   2 of    5 jobs finished (3.5 minutes spent waiting)
#	   3 of    5 jobs finished (5.5 minutes spent waiting)
#	   5 of    5 jobs finished (7.5 minutes spent waiting)
# done. [454.1 seconds]
#
# You can now use dnaorg_annotate.pl to annotate genomes with the models that
# you've created here.
#
#
# Output printed to screen saved in:                                                                                     NC_001346.dnaorg_build.log
# List of executed commands saved in:                                                                                    NC_001346.dnaorg_build.cmd
# List and description of all output files saved in:                                                                     NC_001346.dnaorg_build.list
# File with list of options that must be kept consistent between dnaorg_build.pl and dnaorg_annotate.pl runs saved in:   NC_001346.dnaorg_build.consopts
# CM file #1, cds#1 saved in:                                                                                            NC_001346.dnaorg_build.0.cm
# CM file #2, cds#2 saved in:                                                                                            NC_001346.dnaorg_build.1.cm
# CM file #3, cds#3.1 saved in:                                                                                          NC_001346.dnaorg_build.2.cm
# CM file #4, cds#3.2 saved in:                                                                                          NC_001346.dnaorg_build.3.cm
# CM file #5, cds#4 saved in:                                                                                            NC_001346.dnaorg_build.4.cm
#
# All output files created in directory ./NC_001346/
#
# CPU time:  00:08:11.72
#            hh:mm:ss
# 
# DNAORG-SUCCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We use the -c option because MSV is a circular genome; more
information on the annotation of circular genomes is given below.

The output printed to the screen is also saved to the file
NC_001346/NC_001346.dnaorg_build.log. It explains the three main steps
the script performs as they are being performed, and then outputs a
list of some of the output files created by the script.  For a
complete list see NC_001346.dnaorg_build.list.  All output files have
been created in the subdirectory 'NC_001346/'. A particularly
important file is the NC_001346.dnaorg_build.cmd file that includes
all of the commands executed by the script during its execution.

The three main steps performed by the script are:
1. 'Gathering information on reference using edirect'. In this step, 
   the edirect tools (esearch, efetch, xtract, etc.) are used to
   extract information from the NCBI databases for the reference
   accession NC_001346. The output files ending in .ftable and
   .length include this information.

2. 'Fetching and processing the reference genome'. In this step, the
   actual sequence for NC_001346 is fetched using a program called
   esl-fetch.cds.pl, and the CDS features are extracted from this
   sequence and saved separately to be used to create homology
   models in step 3. The output files with the .fa and .stk suffixes
   are created in this step.

3. 'Submitting jobs to build models to compute farm and waiting for
   them to finish' In this step, commands to build models for each
   of the reference CDS are submitted to the compute farm, and the
   program waits for them to finish. These commands can take a while
   if the farm is busy and they have to wait to run for a long time.

4. 'Submitting jobs to calibrate models to compute farm and waiting
   for them to finish' In this step, commands to 'calibrate' the models 
   is submitted to the farm and the script waits for them to
   finish. A 'calibration' is a simulation used to determine E-value
   parameters for the model. These commands can take a while, even
   once they start running.

During steps 3 and 4, the user can open a separate shell and check
the status of the jobs if they'd like, using the 'qstat' program. 
For example:

> qstat
#job-ID  prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID 
#------------------------------------------------------------------------------------------------------------------------------------------------
187690 1.50472 c.NC_00134 nawrocke     r     05/20/2016 11:59:34 unified@sge581.be-md.ncbi.nlm.                                    4        
187692 1.50472 c.NC_00134 nawrocke     r     05/20/2016 11:59:34 unified@sge151.be-md.ncbi.nlm.                                    4        
187694 1.50472 c.NC_00134 nawrocke     r     05/20/2016 11:59:34 unified@sge109.be-md.ncbi.nlm.                                    4        

The 'r' value in the 'state' column indicates that the jobs are
currently running on a compute farm node. The '4' in the 'slots'
column indicates that they are running on 4 processors. If you
notice that the jobs are in a 'qw' state, this means that they are
currently waiting in the queue for a computer to run on. If they are
stuck in this 'qw' state for longer than you would like, it is most
likely due to the fact that they are waiting for a computer with 4
free processors to run on. In this case, you can cancel the current
command by hitting 'Ctrl-C' on the shell and rerun the command
with the '-n 1' command line option (like this:
'dnaorg_build.pl -c -f -n 1 NC_001346')
to specify that only 1 processor is necessary, although the
calibrations will take roughly 4 times as long. The '-f' option in
the above command is necessary to allow the script to overwrite the
previous files created by the run you aborted.

=======================================================
Example 2, Step 2 of 2. Run dnaorg_annotate.pl for MSV 
=======================================================

Once the calibrations are finished, we can run dnaorg_annotate.pl
to annotate other MSV genomes. For this example, we will annotate
5 genomes, one of which is the reference genome. These are in 
in testfiles/msv.example.fa.

To run dnaorg_annotate.pl do:

> dnaorg_annotate.pl --infasta --refaccn NC_001346 --dirout da-msv.example --dirbuild NC_001346 -c $DNAORGDIR/dnaorg_scripts/testfiles/msv.example.fa
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_annotate.pl :: annotate sequences based on a reference annotation
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 09:59:43 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
# fasta file with sequences to annotate (--infasta):                               /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg_scripts/testfiles/msv.example.fa
# genome is closed (a.k.a. circular):                                              yes [-c]
# output directory specified as:                                                   da-msv.example [--dirout]
# output directory used for dnaorg_build.pl:                                       NC_001346 [--dirbuild]
# single cmdline argument is a fasta file of sequences, not a list of accessions:  yes [--infasta]
# specify reference accession is <s>:                                              NC_001346 [--refaccn]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Verifying options are consistent with options used for dnaorg_build.pl                ... done. [0.0 seconds]
# Processing input fasta file                                                           ... done. [0.0 seconds]
# Fetching all sequences and processing the reference genome                            ... done. [0.0 seconds]
# Skipping verification that CMs created for current reference NC_001346 (--infasta)    ... done. [0.0 seconds]
# Submitting 5 cmscan jobs to the farm                                                  ... done. [0.4 seconds]
# Waiting a maximum of 500 minutes for all farm jobs to finish                          ... 
#	   4 of    5 jobs finished (0.2 minutes spent waiting)
#	   5 of    5 jobs finished (0.5 minutes spent waiting)
# Parsing cmscan results                                                                ... done. [0.0 seconds]
# Calculating predicted feature lengths                                                 ... done. [0.0 seconds]
# Fetching cmscan predicted hits into fasta files                                       ... done. [0.0 seconds]
# Combining predicted exons into CDS                                                    ... done. [0.1 seconds]
# Combining predicted mature peptides into CDS                                          ... done. [0.0 seconds]
# Identifying errors associated with incomplete alignment to the model                  ... done. [0.0 seconds]
# Identifying internal starts/stops in coding sequences                                 ... done. [1.5 seconds]
# Correcting homology search stop codon predictions to account for observed stop codons ... done. [0.0 seconds]
# Identifying overlap and adjacency errors                                              ... done. [0.0 seconds]
# Finalizing annotations and validating error combinations                              ... done. [0.0 seconds]
# Fetching corrected matches into fasta files                                           ... done. [0.0 seconds]
# Combining corrected exons into CDS                                                    ... done. [0.1 seconds]
# Combining corrected mature peptides into CDS                                          ... done. [0.0 seconds]
# Running and parsing BLASTX                                                            ... done. [0.7 seconds]
# Translating corrected nucleotide features into protein sequences                      ... done. [1.7 seconds]
# Generating error code output                                                          ... done. [0.0 seconds]
# Generating tabular annotation output                                                  ... done. [0.0 seconds]
# Generating feature table output                                                       ... done. [0.0 seconds]
#
# Annotated 4 accessions:
#      2 PASS (0.500) listed in da-msv.example/da-msv.example.dnaorg_annotate.ap.seqlist
#      2 FAIL (0.500) listed in da-msv.example/da-msv.example.dnaorg_annotate.af.seqlist
#
# Output printed to screen saved in:                                              da-msv.example.dnaorg_annotate.log
# List of executed commands saved in:                                             da-msv.example.dnaorg_annotate.cmd
# List and description of all output files saved in:                              da-msv.example.dnaorg_annotate.list
# information on blast and CM hits for CDS features in tabular format saved in:   da-msv.example.dnaorg_annotate.blastx.tbl
# All annotations in tabular format saved in:                                     da-msv.example.dnaorg_annotate.tbl
# Summary of all annotations saved in:                                            da-msv.example.dnaorg_annotate.tbl.summary
# Annotations for all sequences with >= 1 failure in tabular format saved in:     da-msv.example.dnaorg_annotate.fail.tbl
# Annotations for all sequences with >= 1 error in tabular format saved in:       da-msv.example.dnaorg_annotate.error.tbl
# List of errors, one line per sequence saved in:                                 da-msv.example.dnaorg_annotate.peraccn.errors
# List of errors, one line per error saved in:                                    da-msv.example.dnaorg_annotate.all.errors
# Summary of all errors saved in:                                                 da-msv.example.dnaorg_annotate.errors.summary
# Sequin feature table output for passing sequences saved in:                     da-msv.example.dnaorg_annotate.ap.sqtable
# Sequin feature table output for failing sequences (minimal) saved in:           da-msv.example.dnaorg_annotate.af.sqtable
# Sequin feature table output for failing sequences (verbose) saved in:           da-msv.example.dnaorg_annotate.long.sqtable
# list of passing sequences saved in:                                             da-msv.example.dnaorg_annotate.ap.seqlist
# list of failing sequences saved in:                                             da-msv.example.dnaorg_annotate.af.seqlist
# list of errors in the sequence tables saved in:                                 da-msv.example.dnaorg_annotate.errlist
#
# All output files created in directory ./da-msv.example/
#
# CPU time:  00:00:34.99
#            hh:mm:ss
# 
# DNAORG-SUCCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is an explanation of all of the options used in the 
above command:

'--infasta' specifies that the one command line argument
(NC_001346.fa) is a fasta file of sequences to annotate instead of
a list of accessions.

'--refaccn' is a required option when --infasta is used. It must
be followed by the reference accession that was passed into
dnaorg_build.pl, 'NC_001346' in this case.

'--local' specifies that the cmscan homology search steps will be
run locally (not submitted to the compute farm).

'--dirout NC_001346-out1' specifies that a new directory should be
created for the output named 'NC_001346-out1'. 

'--dirbuild NC_001346' specifies that a directory named 'NC_001346'
already exists, and was created by an earlier dnaorg_build.pl run.
That dnaorg_build.pl run must be completed (any cmcalibrate steps
submitted to the compute farm by dnaorg_build.pl must have
successfully completed.)

'-c' specifies that the sequences are circular (this option was also
used in example run #1).

The dnaorg_annotate.pl output above begins with a description of the
many steps it performs, with timings of each step.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Verifying options are consistent with options used for dnaorg_build.pl                ... done. [0.0 seconds]
# Processing input fasta file                                                           ... done. [0.0 seconds]
# Fetching all sequences and processing the reference genome                            ... done. [0.0 seconds]
# Skipping verification that CMs created for current reference NC_001346 (--infasta)    ... done. [0.0 seconds]
# Submitting 5 cmscan jobs to the farm                                                  ... done. [0.4 seconds]
# Waiting a maximum of 500 minutes for all farm jobs to finish                          ... 
#	   4 of    5 jobs finished (0.2 minutes spent waiting)
#	   5 of    5 jobs finished (0.5 minutes spent waiting)
# Parsing cmscan results                                                                ... done. [0.0 seconds]
# Calculating predicted feature lengths                                                 ... done. [0.0 seconds]
# Fetching cmscan predicted hits into fasta files                                       ... done. [0.0 seconds]
# Combining predicted exons into CDS                                                    ... done. [0.1 seconds]
# Combining predicted mature peptides into CDS                                          ... done. [0.0 seconds]
# Identifying errors associated with incomplete alignment to the model                  ... done. [0.0 seconds]
# Identifying internal starts/stops in coding sequences                                 ... done. [1.5 seconds]
# Correcting homology search stop codon predictions to account for observed stop codons ... done. [0.0 seconds]
# Identifying overlap and adjacency errors                                              ... done. [0.0 seconds]
# Finalizing annotations and validating error combinations                              ... done. [0.0 seconds]
# Fetching corrected matches into fasta files                                           ... done. [0.0 seconds]
# Combining corrected exons into CDS                                                    ... done. [0.1 seconds]
# Combining corrected mature peptides into CDS                                          ... done. [0.0 seconds]
# Running and parsing BLASTX                                                            ... done. [0.7 seconds]
# Translating corrected nucleotide features into protein sequences                      ... done. [1.7 seconds]
# Generating error code output                                                          ... done. [0.0 seconds]
# Generating tabular annotation output                                                  ... done. [0.0 seconds]
# Generating feature table output                                                       ... done. [0.0 seconds]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step verifies that some special command line options used
with dnaorg_build.pl were also used with dnaorg_annotate.pl, e.g. 
'-c'. For details on this, see 'OPTIONS THAT MUST BE USED
CONSISTENTLY IN BOTH dnaorg_build.pl AND dnaorg_annotate.pl' below

The next several steps gather information about the reference
accession, the existing annotations of all accessions, and validate
the CM file in the NC_001346/ directory that is about to be used for
annotation was indeed made for the current reference sequence.

The script then runs the Infernal 'cmscan' program to annotate all
the 'features' of each of the accessions. It does this by submitting 5
cmscan jobs to the compute farm and waits for them to finish before
proceeding. To disable this parallelization, use the --local option.

The cmscan results are then parsed and the predicted hits are
fetched into fasta files. These sequence files are then examined for
the validity of their start and stop codons and also to identify any
early in-frame stop codons. The results of these examinations are
used to 'correct' predictions where early in-frame stop codons were
found. The 'corrected' sequences are then fetched, aligned,
translated, and then those translations are aligned.

Finally the script outputs feature tables with annotation
information. 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Annotated 4 accessions:
#      2 PASS (0.500) listed in da-msv.example/da-msv.example.dnaorg_annotate.ap.seqlist
#      2 FAIL (0.500) listed in da-msv.example/da-msv.example.dnaorg_annotate.af.seqlist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case, 2 sequences passed and 2 failed. If >= 1 sequences
passed a feature table file with all passing sequences will be created
witht he suffix .dnaorg_annotate.ap.sqtable. If >= 1 sequences
passed a feature table file with all failing sequences will be
created. The passing and failing sequences will be listed, one per
line, the files ending with .dnaorg_annotate.ap.seqlist and
.dnaorg_annotate.af.seqlist. 

For this example those files are:
 da-msv.example/da-msv.example.dnaorg_annotate.ap.seqlist
 da-msv.example/da-msv.example.dnaorg_annotate.af.seqlist
 da-msv.example/da-msv.example.dnaorg_annotate.ap.sqtable
 da-msv.example/da-msv.example.dnaorg_annotate.af.sqtable

#######################################################################
Example 3. Classifying and annotating Norovirus sequences from
           different Norovirus genotypes.                                  
#######################################################################

For this example run, we will use the scripts in an alternative way
by first building a set of multiple models for multiple RefSeqs, one
for each of the 6 Norovirus RefSeqs, and then using
dnaorg_classify.pl in a special mode to both classify sequences and
then annotate them (by internally calling dnaorg_annotate.pl).

#=======================================================
Example 3, Steps 1 to 3 of 4. DO NOT EXECUTE, JUST READ
                              THROUGH THIS SECTION
#=======================================================

Unlike the other examples, you won't repeat the first steps required 
to do this example, because they have already been done for your (and
they would take too long to repeat here). The files resulting from
these steps were installed by 'install.sh' here:

$DNAORGDIR/dnaorg-build-directories

And you will use those files to complete step 4 below.

(If somehow you did not install the 'dnaorg-build-directories'
directory, you can do that using this command:
git clone https://github.com/nawrockie/dnaorg-build-directories.git)

In the dnaorg-build-directories directory there is a subdirectory
named 'norovirus-builds'. That directory was created by following
the steps in the file noro-00RECREATE.sh that is itself within that 
'norovirus-builds' directory. Here is that file:

> cat $DNAORGDIR/dnaorg-build-directories/norovirus-builds/noro-00RECREATE.sh 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EPN, Wed Mar 14 15:13:52 2018 [dnaorg_scripts version 0.28 -->] current
# EPN, Tue Feb 13 12:58:45 2018 [dnaorg_scripts version 0.26] previous
#
# If you're updating the models. Do these steps first:
# > mv norovirus-builds old-norovirus-builds
# > cp old-norovirus-builds/noro-00RECREATE.sh .
# > cp old-norovirus-builds/norovirus_refseq_list.txt .
# 
# Step 1: create the norovirus-builds directory using dnaorg_classify.pl --onlybuild:
dnaorg_classify.pl --dirout norovirus-builds --onlybuild norovirus_refseq_list.txt 

# Step 2: create the matpept files we need for dnaorg_build.pl and copy
# them into the norovirus-builds directory:
for a in \
    NC_001959 \
    NC_008311 \
    NC_029645 \
    NC_029646 \
    NC_029647 \
    NC_031324 \
    NC_039475 \
    NC_039476 \
    NC_039477 \
    ; do 
    dnaorg_get_matpepts.pl -f $a
    cp $a.matpept.in norovirus-builds
done

# Step 3: create the individual dnaorg_build.pl output directories by calling dnaorg_build.pl a bunch of times,
# use the -n 1 :
cd norovirus-builds;
for a in \
    NC_001959 \
    NC_008311 \
    NC_029645 \
    NC_029646 \
    NC_029647 \
    NC_031324 \
    NC_039475 \
    NC_039476 \
    NC_039477 \
    ; do 
    qsub -N build.$a -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e build.$a.err -l h_rt=28800,mem_free=8G,h_vmem=16G -m n "dnaorg_build.pl -f -n 1 --matpept $a.matpept.in --dfeat gene $a > $a.dnaorg_build.out"
done
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 3 steps listed in the above file are the first 3 steps of this
example. You *could* run all three of them by executing the
noro-00RECREATE.sh script, but you do not and should not do that (that
file just exists for reference). 

In Step 1, dnaorg_classify.pl was run in --onlybuild mode to create
the 'norovirus-builds' directory, similarly to step 1 of example 1.

In Step 2, special 'matpept' files are created using the
dnaorg_get_matpepts.pl script for each of the six Norovirus RefSeq
accessions (NC_001959 NC_008311 NC_029645 NC_029646 NC_029647
NC_031324) These matpept files are required input to
dnaorg_build.pl. An example of a matpept file was used in Example 4,
where it was provided. The dnaorg_get_matpepts.pl script can generate
a matpete file for any RefSeq that has mature peptide annotation. An
example usage is: "dnaorg_get_matpepts.pl NC_001959', which will
create the file NC_001959.matpept.in.

In Step 3, dnaorg_build.pl is called for each of the six Norovirus
RefSeq accessions, using the 'qsub' program which submits each
dnaorg_build.pl run to the compute farm. This step will probably only
work at NCBI.  

After all three steps are done, the directory norovirus-builds is now
ready to be used with dnaorg_classify.pl in a special mode, using the
-A option, to classify and annotate Norovirus sequences.

This next step *is* recommended for the reader to try.

======================================================= 
Example 3, Step 4. Run dnaorg_classify.pl -A for Norovirus sequences of different
genotypes.  
=======================================================

For this step, we will use the -A option for dnaorg_classify.pl,
which is described if you use the -h option to
dnaorg_classify.pl. From above: 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
options for automatically running dnaorg_annotate.pl for classified sequences:
  -A <s>         : annotate using dnaorg_build.pl build directories in <s> after classifying
  --optsA <s>    : read additional dnaorg_annotate.pl options from file <s>
  --reflistA <s> : only annotate seqs that match to RefSeqs listed in <s>
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case, we want to use the dnaorg_build.pl directories that
exist in the norovirus-builds directory referred to above:

$DNAORGDIR/dnaorg-build-directories/norovirus-builds

We will use dnaorg_classify.pl to classify and annotate the 9
Norovirus sequences in the file testfiles/noro.9.fa with the
following command. Note that we need to supply the path to
norovirus-builds twice, once each with the -A and --dirbuild
options.

> dnaorg_classify.pl -A $DNAORGDIR/dnaorg-build-directories/norovirus-builds --infasta $DNAORGDIR/dnaorg_scripts/testfiles/noro.9.fa --dirbuild $DNAORGDIR/dnaorg-build-directories/norovirus-builds --dirout noro.9
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_classify.pl :: classify sequences using an HMM library of RefSeqs
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 10:15:47 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
# REQUIRED name for output directory to create is <s>:     noro.9 [--dirout]
# fasta file with sequences to classify is <s>:            /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg_scripts/testfiles/noro.9.fa [--infasta]
# specify directory with HMM library is <s>:               /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg-build-directories/norovirus-builds [--dirbuild]
# annotate after classifying using build dirs in dir <s>:  /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg-build-directories/norovirus-builds [-A]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parsing RefSeq list                                                    ... done. [0.0 seconds]
# Verifying build directories exist (-A)                                 ... done. [0.5 seconds]
# Processing fasta file to get sequence lengths                          ... done. [0.1 seconds]
# Submitting 2 nhmmscan jobs to the farm                                 ... done. [0.2 seconds]
# Waiting a maximum of 500 minutes for all farm jobs to finish           ... 
#	   0 of    2 jobs finished (0.2 minutes spent waiting)
#	   0 of    2 jobs finished (0.5 minutes spent waiting)
#	   1 of    2 jobs finished (0.8 minutes spent waiting)
#	   2 of    2 jobs finished (1.0 minutes spent waiting)
# Creating tabular output file                                           ... done. [0.0 seconds]
# Creating seqlists and other output files                               ... done. [0.4 seconds]
#
# Number of input sequences assigned to each RefSeq:
#
#model           num-PASS    num-FAIL    annotating?
#------------  ----------  ----------    -----------
NC_001959               2           0            yes
NC_008311               2           0            yes
NC_029645               2           0            yes
NC_029646               0           0            yes
NC_029647               0           0            yes
NC_031324               1           0            yes
NC_039475               0           0            yes
NC_039476               0           0            yes
NC_039477               2           0            yes
NON-ASSIGNED            0           0             no
#------------  ----------  ----------    -----------
SUM-ASSIGNED            9           0            yes
SUM-ASSIGNED            0           0             no
SUM-ASSIGNED            9           0           both
#
#
# Running dnaorg_annotate.pl 5 time(s) to annotate sequences             ... 
# dnaorg_annotate.pl :: annotate sequences based on a reference annotation
# dnaorg 0.44 (Jan 2019)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TRUNCATED OUTPUT

This output is similar to the dnaorg_classify.pl output above. Except
when the script is done classifying the sequences, it moves directly
into dnaorg_annotate.pl. At this point it outputs that it is now:
"Running dnaorg_annotate.pl 5 time(s) to annotate sequences". The
output of the five runs of dnaorg_annotate.pl is not shown. After they
are all finished, the output is:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DNAORG-SUCCESS
done. [114.6 seconds]
#
# Number of annotated sequences that PASSed/FAILed dnaorg_annotate.pl:
#
#model          num-annot  num-C-PASS  num-C-FAIL  num-A-PASS  num-A-FAIL  num-A-FAIL-CO
#------------  ----------  ----------   ---------  ----------  ----------  -------------
NC_001959               2           2           0           2           0              0
NC_008311               2           2           0           1           1              0
NC_029645               2           2           0           2           0              0
NC_031324               1           1           0           1           0              0
NC_039477               2           2           0           2           0              0
#------------  ----------  ----------   ---------  ----------  ----------  -------------
*ALL*                   9           9           0           8           1              0
#
# Output printed to screen saved in:                                                                                                                             noro.9.dnaorg_classify.log
# List of executed commands saved in:                                                                                                                            noro.9.dnaorg_classify.cmd
# List and description of all output files saved in:                                                                                                             noro.9.dnaorg_classify.list
# Fasta file with sequences to classify (copy of /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg_scripts/testfiles/noro.9.fa) saved in:   noro.9.dnaorg_classify.fa
# Per-sequence hit and classification information, human readable saved in:                                                                                      noro.9.dnaorg_classify.rdb.infotbl
# Per-sequence hit and classification information, tab-delimited saved in:                                                                                       noro.9.dnaorg_classify.tab.infotbl
# List of errors (unexpected features that cause failure) for all sequences saved in:                                                                            noro.9.dnaorg_classify.all.errlist
# List of errors (unexpected features that cause failure) for sequences that will not be annotated saved in:                                                     noro.9.dnaorg_classify.noannot.errlist
# List of sequences that will not be annotated saved in:                                                                                                         noro.9.dnaorg_classify.noannot.seqlist
# List of PASSing seqs for NC_001959 saved in:                                                                                                                   noro.9.dnaorg_classify.NC_001959.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_001959 saved in:                                                                                                   noro.9.dnaorg_classify.NC_001959.cp.fa
# List of all seqs for NC_001959 saved in:                                                                                                                       noro.9.dnaorg_classify.NC_001959.seqlist
# Fasta file with all seqs assigned to NC_001959 saved in:                                                                                                       noro.9.dnaorg_classify.NC_001959.fa
# List of PASSing seqs for NC_008311 saved in:                                                                                                                   noro.9.dnaorg_classify.NC_008311.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_008311 saved in:                                                                                                   noro.9.dnaorg_classify.NC_008311.cp.fa
# List of all seqs for NC_008311 saved in:                                                                                                                       noro.9.dnaorg_classify.NC_008311.seqlist
# Fasta file with all seqs assigned to NC_008311 saved in:                                                                                                       noro.9.dnaorg_classify.NC_008311.fa
# List of PASSing seqs for NC_029645 saved in:                                                                                                                   noro.9.dnaorg_classify.NC_029645.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_029645 saved in:                                                                                                   noro.9.dnaorg_classify.NC_029645.cp.fa
# List of all seqs for NC_029645 saved in:                                                                                                                       noro.9.dnaorg_classify.NC_029645.seqlist
# Fasta file with all seqs assigned to NC_029645 saved in:                                                                                                       noro.9.dnaorg_classify.NC_029645.fa
# List of PASSing seqs for NC_031324 saved in:                                                                                                                   noro.9.dnaorg_classify.NC_031324.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_031324 saved in:                                                                                                   noro.9.dnaorg_classify.NC_031324.cp.fa
# List of all seqs for NC_031324 saved in:                                                                                                                       noro.9.dnaorg_classify.NC_031324.seqlist
# Fasta file with all seqs assigned to NC_031324 saved in:                                                                                                       noro.9.dnaorg_classify.NC_031324.fa
# List of PASSing seqs for NC_039477 saved in:                                                                                                                   noro.9.dnaorg_classify.NC_039477.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_039477 saved in:                                                                                                   noro.9.dnaorg_classify.NC_039477.cp.fa
# List of all seqs for NC_039477 saved in:                                                                                                                       noro.9.dnaorg_classify.NC_039477.seqlist
# Fasta file with all seqs assigned to NC_039477 saved in:                                                                                                       noro.9.dnaorg_classify.NC_039477.fa
# List of RefSeqs in the HMM library saved in:                                                                                                                   noro.9.dnaorg_classify.all.refseqs
# List of sequences that were sorted into seqlists saved in:                                                                                                     noro.9.dnaorg_classify.all.seqs
# annotation results for NC_001959 sequences that pass dnaorg_annotate.pl saved in:                                                                              noro.9-NC_001959.dnaorg_annotate.ap.sqtable
# annotation results for NC_001959 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    noro.9-NC_001959.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      noro.9-NC_001959.dnaorg_annotate.long.sqtable
# list of NC_001959 sequences that pass dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_001959.dnaorg_annotate.ap.seqlist
# list of NC_001959 sequences that fail dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_001959.dnaorg_annotate.af.seqlist
# list of NC_001959 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            noro.9-NC_001959.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 noro.9-NC_001959.dnaorg_annotate.errlist
# annotation results for NC_008311 sequences that pass dnaorg_annotate.pl saved in:                                                                              noro.9-NC_008311.dnaorg_annotate.ap.sqtable
# annotation results for NC_008311 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    noro.9-NC_008311.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      noro.9-NC_008311.dnaorg_annotate.long.sqtable
# list of NC_008311 sequences that pass dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_008311.dnaorg_annotate.ap.seqlist
# list of NC_008311 sequences that fail dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_008311.dnaorg_annotate.af.seqlist
# list of NC_008311 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            noro.9-NC_008311.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 noro.9-NC_008311.dnaorg_annotate.errlist
# annotation results for NC_029645 sequences that pass dnaorg_annotate.pl saved in:                                                                              noro.9-NC_029645.dnaorg_annotate.ap.sqtable
# annotation results for NC_029645 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    noro.9-NC_029645.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      noro.9-NC_029645.dnaorg_annotate.long.sqtable
# list of NC_029645 sequences that pass dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_029645.dnaorg_annotate.ap.seqlist
# list of NC_029645 sequences that fail dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_029645.dnaorg_annotate.af.seqlist
# list of NC_029645 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            noro.9-NC_029645.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 noro.9-NC_029645.dnaorg_annotate.errlist
# annotation results for NC_031324 sequences that pass dnaorg_annotate.pl saved in:                                                                              noro.9-NC_031324.dnaorg_annotate.ap.sqtable
# annotation results for NC_031324 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    noro.9-NC_031324.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      noro.9-NC_031324.dnaorg_annotate.long.sqtable
# list of NC_031324 sequences that pass dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_031324.dnaorg_annotate.ap.seqlist
# list of NC_031324 sequences that fail dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_031324.dnaorg_annotate.af.seqlist
# list of NC_031324 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            noro.9-NC_031324.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 noro.9-NC_031324.dnaorg_annotate.errlist
# annotation results for NC_039477 sequences that pass dnaorg_annotate.pl saved in:                                                                              noro.9-NC_039477.dnaorg_annotate.ap.sqtable
# annotation results for NC_039477 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    noro.9-NC_039477.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      noro.9-NC_039477.dnaorg_annotate.long.sqtable
# list of NC_039477 sequences that pass dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_039477.dnaorg_annotate.ap.seqlist
# list of NC_039477 sequences that fail dnaorg_annotate.pl saved in:                                                                                             noro.9-NC_039477.dnaorg_annotate.af.seqlist
# list of NC_039477 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            noro.9-NC_039477.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 noro.9-NC_039477.dnaorg_annotate.errlist
#
# All output files created in directory ./noro.9/
#
# CPU time:  00:02:57.18
#            hh:mm:ss
# 
# DNAORG-SUCCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After finishing, the Sequin format feature table files
(e.g. noro-9-NC_001959.dnaorg_annotate.ap.sqtable) are in
the noro.9 output directory.

The output directories created by dnaorg_annotate.pl with all of the
output files from that script are within the noro.9 output
directory, e.g. noro.9/noro.9-NC_001959.pass/.

You can view the output of the six dnaorg_annotate.pl runs in the
files that end in .log in those directories,
e.g. noro.9/noro.9-NC_001959.pass/noro.9-NC_001959.pass.dnaorg_annotate.log

The dnaorg_classify.pl command in step 4 above could have been
repeated with a list of accessions instead of an input fasta file. 
The relevant command would be:
#perl /panfs/pan1/dnaorg/virseqannot/code/dnaorg_scripts/dnaorg_classify.pl -A /panfs/pan1/dnaorg/virseqannot/dnaorg-build-directories/norovirus-builds --inlist testfiles/noro.9.ntlist --dirbuild /panfs/pan1/dnaorg/virseqannot/dnaorg-build-directories/norovirus-builds --dirout noro.9.list

The --optsA <s> option can be used to supply additional options 
to dnaorg_annotate.pl by creating a file that contains a single
line with the additional options and arguments as necessary that
should be passed to all runs of dnaorg_annotate.pl. For example,
if that file were called annotate.opts, add "--optsA annotate.opts" 
to the dnaorg_classify.pl command.

##########################################################################
OPTIONS THAT MUST BE USED CONSISTENTLY IN BOTH dnaorg_build.pl AND 
dnaorg_annotate.pl                                                 

Some options exist for both dnaorg_build.pl and dnaorg_annotate.pl
and must be used consistently with both programs. These options are:

 -c             : genome is closed (a.k.a circular)
 --matpept <s>  : read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships
 --nomatpept    : ignore mat_peptide information in reference annotation
 --xfeat <s>    : build models of additional qualifiers in string <s>
 --dfeat <s>   : annotate qualifiers in <s> from duplicates (e.g. gene from CDS)

If you use either -c, --matpept <s>, --nomatpept, --xfeat <s>, or
--dfeat <s> with dnaorg_build.pl for reference sequence x, then
you must also use it
when you run dnaorg_annotate.pl to annotate sequences according
to reference sequence x. Further with --matpept <s>, the file
supplied as <s> must be the same exact file for both dnaorg_build.pl
and dnaorg_annotate.pl.

The dnaorg_annotate.pl script enforces that these options are used
consistently. dnaorg_annotate.pl will fail with an informative error
message if they are not used consistently. For example, if you run
dnaorg_build.pl like this (the command from example run 1, above):

'perl dnaorg_build.pl -c NC_001346'

and subsequently run dnaorg_annotate.pl without -c, like this:

'perl dnaorg_annotate.pl NC_001346.ntlist'

The program will quickly exit with the following output explaining
the problem:

ERROR, the -c option was used when dnaorg_build.pl was run (according to file NC_001346/NC_001346.dnaorg_build.consopts).
You must also use it with dnaorg_annotate.pl.

Similarly, if you use dnaorg_classify.pl with the -A option,
the dnaorg_classify.pl script will call dnaorg_annotate.pl with
the appropriate values for these special options based on 
the options used to create the dnaorg_build.pl directories.

