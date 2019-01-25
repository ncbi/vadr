EPN, Thu Jan 24 13:11:35 2019

dnaorg_scripts 0.44 README.txt 

Organization of this file:

INTRODUCTION
INSTALLATION
SETTING UP ENVIRONMENT VARIABLES
VERIFYING SUCCESSFUL INSTALLATION WITH TEST RUNS
USAGE AND OPTIONS OF dnaorg_classify.pl
EXAMPLE RUN OF dnaorg_classify.pl
POSSIBLE ERRORS
NAMING CONVENTION FOR SEQUENCES IN FASTA FILES 
FORMAT OF dnaorg_classify.pl's .infotbl OUTPUT FILES
GETTING MORE INFORMATION

Questions:
email Eric Nawrocki: eric.nawrocki@nih.gov

##############################################################################
INTRODUCTION

This is documentation for dnaorg_scripts, a suite of tools for
classifying and analyzing viral sequences based on viral reference
sequences. (It is not able to detect or analyze novel viruses.)

Author: Eric Nawrocki and Alejandro Schaffer

Git repository for dnaorg_scripts:
https://github.com/nawrockie/dnaorg_scripts.git

This file gives examples of using the dnaorg scripts
dnaorg_classify.pl, dnaorg_build.pl and
dnaorg_annotate.pl. dnaorg_classify.pl can be used to classify a given
sequence or accession as most similar to one of a set of "reference"
sequences. dnaorg_build.pl and dnaorg_annotate.pl can be used annotate
viral genomes of species S given the most closely related single
reference GenBank entry for S. In this context, the meaning of
"reference" sequence is that it has all the annotations expected and
for circular genomes that position 1 is consistent with the accepted
convention about where the origin should be. If there are multiple
reference sequences for S, one runs dnaorg_build.pl and
dnaorg_annotate.pl separately for each reference sequence. In some
places, the technical term "RefSeq" is used instead of the less formal
"reference sequence".

The overarching objectives of *dnaorg*.pl are to compute all information
needed to produce annotations and to output that information via
Sequin formatted feature tables and FASTA files;

The dnaorg scripts do not produce an updated GenBank flat file.

IMPORTANT: Currently, most of the functionality of the dnaorg scripts
is limited to their use internally at NCBI. We hope to change this in
future versions. The examples below in the EXAMPLE RUNS section can
NOT be run by users external to NCBI. However, most of this
README.txt, including the INSTALLATION section is still relevant to
all users. See the file 'README-norovirus.txt' for examples of how
users external to NCBI can use dnaorg scripts to classify and annotate
Norovirus sequences.

##############################################################################
INSTALLATION

The file 'install.sh' is an executable file for installing dnaorg_scripts
and its dependencies. Copy that file into the directory in which you 
want to install dnaorg_scripts. You may want to call that
'dnaorg-0.44'. Then move into that directory and run the
command './install.sh'. That will create several directories in the
current directory.

'install.sh' will download and install the software packages
dnaorg_scripts, Infernal, HMMER 3.1b2, esl-fetch-cds,
and esl-epn-translate, as well as the required perl modules
epn-options and Bio-Easel.

When install.sh is finished running it will print important
instructions to the screen that explain how to modify your environment
variables so that you can run the dnaorg_scripts scripts.
These are discussed further below in the SETTING UP ENVIRONMENT
VARIABLES section.

##############################################################################
SETTING UP ENVIRONMENT VARIABLES

As mentioned above, when you run 'install.sh' or any of the optional install
scripts, instructions will be output about how to change your environment
variables so that you can run the dnaorg scripts. Those instructions are
also included below for reference, but without the actual path to
where you ran install.sh (below it is replaced with 
"<full path to directory in which you ran install.sh>")

--
********************************************************
The final step is to update your environment variables.
(See dnaorg_scripts/README.txt for more information.)

If you are using the bash shell, add the following
lines to the '.bashrc' file in your home directory:

export DNAORGDIR=<full path to directory in which you ran install.sh>
export DNAORGBUILDDIR="$DNAORGDIR"/dnaorg-build-directories
export PERL5LIB="$DNAORGDIR"/dnaorg_scripts:"$DNAORGDIR"/epn-options:"$DNAORGDIR"/Bio-Easel/blib/lib:"$DNAORGDIR"/Bio-Easel/blib/arch:"$PERL5LIB"
export PATH="$DNAORGDIR"/dnaorg_scripts:"$PATH"

After adding the export lines to your .bashrc file, source that file
to update your current environment with the command:

source ~/.bashrc

---
If you are using the C shell, add the following
lines to the '.cshrc' file in your home directory:

setenv DNAORGDIR <full path to directory in which you ran install.sh>
setenv DNAORGBUILDDIR "$DNAORGDIR"/dnaorg-build-directories
setenv PERL5LIB "$DNAORGDIR"/dnaorg_scripts:"$DNAORGDIR"/epn-options:"$DNAORGDIR"/Bio-Easel/blib/lib:"$DNAORGDIR"/Bio-Easel/blib/arch:"$PERL5LIB"
setenv PATH "$DNAORGDIR"/dnaorg_scripts:"$PATH"


After adding the setenv lines to your .cshrc file, source that file
to update your current environment with the command:

source ~/.cshrc

(To determine which shell you use, type: 'echo $SHELL')


********************************************************
--

If you get an error about PERL5LIB being undefined, change the PERL5LIB
line to add to:
export PERL5LIB="$DNAORGDIR"/dnaorg_scripts-0.44:"$DNAORGDIR"/epn-options:"$DNAORGDIR"/Bio-Easel/blib/lib:"$DNAORGDIR"/Bio-Easel/blib/arch
for .bashrc, OR
setenv PERL5LIB "$DNAORGDIR"/dnaorg_scripts-0.44:"$DNAORGDIR"/epn-options:"$DNAORGDIR"/Bio-Easel/blib/lib:"$DNAORGDIR"/Bio-Easel/blib/arch
for .cshrc. And then do
> source ~/.bashrc
or
> source ~/.cshrc
again.

###########################################################################
VERIFYING SUCCESSFUL INSTALLATION WITH TEST RUNS

The dnaorg package include some tests you can run to make sure that
your installation was successful and that your environment variables 
are set-up correctly.

There are 2 shell scripts for running tests:

     1	do-install-tests-local.sh
     2	do-install-tests-parallel.sh

The dnaorg scripts can be run locally on your computer or in parallel
on a compute farm. These two test files test each of those modes. 
If you plan to run the scripts at locally at least some of the time,
then run 'do-install-tests-local.sh'. If you plan to run the scripts
on a compute farm at least some of the time, then run
'do-install-tests-parallel.sh' 

Below is example output for do-install-tests-local.sh:
> sh $DNAORGDIR/testfiles/do-install-tests-local.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_test.pl :: test dnaorg scripts [TEST SCRIPT]
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 06:05:52 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
# test file:                                                   /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg_scripts/testfiles/noro.r10.local.testin
# forcing directory overwrite:                                 yes [-f]
# build directory, replaces !dirbuild! in test file with <s>:  /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg-build-directories/norovirus-builds [--dirbuild]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Running command  1 [classify-noro-10-local]        ... done. [52.9 seconds]
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
# CPU time:  00:00:53.25
#            hh:mm:ss
# 
# DNAORG-SUCCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The most important line is the line that begins with "# PASS"

# PASS: all 6 files were created correctly.

This means that the test has passed. You should see similar 
lines when you run the other tests. If you do not and need help
figuring out why, email me at eric.nawrocki@nih.gov.

###############################################
USING dnaorg_classify.pl TO CLASSIFY AND ANNOTATE NOROVIRUS SEQUENCES

The dnaorg_classify.pl program can be used to classify and annotate
norovirus sequences. It actually can be used generally for any type of
virus given pre-computed models for that virus. However, it has
currently only been tested extensively for norovirus and the examples
here are limited to norovirus. 

dnaorg_classify.pl takes as input a fasta file and for each sequence
in that file, it classifies the sequence as norovirus by determining
the most similar norovirus RefSeq and then annotates that sequence
based on that RefSeq's annotation. 

The following command provides an example, using the input file 
sequence file noro.9.fa with 9 sequences in the 
$DNAORGDIR/dnaorg_scripts/testfiles directory.

> dnaorg_classify.pl --local -A $DNAORGDIR/dnaorg-build-directories/norovirus-builds --infasta $DNAORGDIR/dnaorg_scripts/testfiles/noro.9.fa --dirbuild   $DNAORGDIR/dnaorg-build-directories/caliciviridae-classify-build --reflistA $DNAORGDIR/dnaorg-build-directories/norovirus-builds/norovirus-builds.dnaorg_classify.ref.list --dirout dc-noro.9 

In this command we are using several command line options:

--local: this specifies that all programs are to be run locally on
         this computer. Remove this option to use the SGE compute farm
         to parallelize the classification and annotation steps, but
         that will likely only work internally at NCBI.

-A <s>:  this tells the script to annotate sequences after classifying
         them using the model files in directory <s>. This directory
         ($DNAORGDIR/dnaorg-build-directories/norovirus-builds) was
         installed by 'install.sh' when you installed the package.
         It includes pre-built models for Norovirus RefSeqs.

--infasta <s>: this specifies that the input sequence file is
               <f>. This file must be in FASTA format.

--dirbuild <s>: this specifies that the directory <s> includes the
                models used for the classification step. This directory 
                ($DNAORGDIR/dnaorg-build-directories/caliciviridae-classify-build) 
                contains models for all Calicivirus RefSeqs, the
                family that includes Norovirus. It was installed by
                'install.sh'. 

--reflistA <s>: this specifies a file that is a list of RefSeqs, only
                sequences that are classified to one of these RefSeqs
                will be annotated. This file
                ($DNAORGDIR/dnaorg-build-directories/norovirus-builds/norovirus-builds.dnaorg_classify.ref.list)
                includes the list of all Norovirus RefSeqs.

--dirout <s>: this specifies that the output directory is named
              <s>. This directory will be created by the script and
              all output files will be placed within it.

You should always use these options when annotation norovirus
sequences. Only two of the arguments ('<s>' above)  will change, those
for --infasta, which will change to the name of the sequence you are
annotating and for --dirout to the name of your desired output directory.

--
POSSIBLE ADDITIONAL OPTIONS FOR SPECIFYING GENOGROUP:

If you know what genogroup your sequences are, and you want
dnaorg_classify.pl to enforce that genogroup during its classification
stage (and throw an error for any sequence that seems to not be that
genogroup), you should use two additional options:

--ecall <s>: where <s> is one of "GI", "GII", "GIII", "GIV", or "GV"

--ecmap <s>: where <s> is $DNAORGDIR/dnaorg-build-directories/norovirus.tax2ref.map 
--

Here is the output from the above command (without the --ecall and
--ecmap options):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_classify.pl :: classify sequences using an HMM library of RefSeqs
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 11:39:18 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
# REQUIRED name for output directory to create is <s>:     dc-noro.9 [--dirout]
# fasta file with sequences to classify is <s>:            /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg_scripts/testfiles/noro.9.fa [--infasta]
# specify directory with HMM library is <s>:               /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg-build-directories/caliciviridae-classify-build [--dirbuild]
# run nhmmscan locally instead of on farm:                 yes [--local]
# annotate after classifying using build dirs in dir <s>:  /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg-build-directories/norovirus-builds [-A]
# only annotate seqs that match to RefSeqs listed in <s>:  /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg-build-directories/norovirus-builds/norovirus-builds.dnaorg_classify.ref.list [--reflistA]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parsing RefSeq list                                                    ... done. [0.0 seconds]
# Verifying build directories exist (-A)                                 ... done. [0.4 seconds]
# Processing fasta file to get sequence lengths                          ... done. [0.1 seconds]
# Running nhmmscan locally                                               ... done. [120.2 seconds]
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
NC_002551               0           0             no
NC_000940               0           0             no
NC_004542               0           0             no
NC_004064               0           0             no
NC_004541               0           0             no
NC_001481               0           0             no
NC_010624               0           0             no
NC_006269               0           0             no
NC_006554               0           0             no
NC_002615               0           0             no
NC_006875               0           0             no
NC_001543               0           0             no
NC_007916               0           0             no
NC_008580               0           0             no
NC_011050               0           0             no
NC_011704               0           0             no
NC_012699               0           0             no
NC_017936               0           0             no
NC_019712               0           0             no
NC_024031               0           0             no
NC_024078               0           0             no
NC_025676               0           0             no
NC_027026               0           0             no
NC_027122               0           0             no
NC_030793               0           0             no
NC_033081               0           0             no
NC_033776               0           0             no
NC_034444               0           0             no
NC_035675               0           0             no
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

The output lists each step the program performs, and then information
on what sequences it has classified to which of the Caliciviridae
RefSeqs. Note that the first 9 lines in the table have 'yes' for
annotation and the remainder have 'no', for example:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#model           num-PASS    num-FAIL    annotating?
#------------  ----------  ----------    -----------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NC_039477               2           0            yes
NC_002551               0           0             no
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This means that sequences assigned to NC_039477 (a norovirus RefSeq)
will be annotated, and those assigned to NC_002551 (a non-norovirus
Calicivirus) will not (although in this case zero sequences were
assigned to NC_002551.

After the classification stage, the program outputs:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Running dnaorg_annotate.pl 5 time(s) to annotate sequences             ... 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Annotation is performed by a different script dnaorg_annotate.pl that
is called internally by dnaorg_classify.pl, which calls
dnaorg_annotate.pl once per Norovirus RefSeq to which at least 1
sequence was assigned.

The output of dnaorg_annotate.pl lists the steps taken by that
program:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dnaorg_annotate.pl :: annotate sequences based on a reference annotation
# dnaorg 0.44 (Jan 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:       Fri Jan 25 11:41:21 2019
# $DNAORGDIR: /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44
#
# fasta file with sequences to annotate (--infasta):                               dc-noro.9/dc-noro.9.dnaorg_classify.NC_001959.fa
# output directory specified as:                                                   dc-noro.9/dc-noro.9-NC_001959 [--dirout]
# output directory used for dnaorg_build.pl:                                       /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg-build-directories/norovirus-builds/NC_001959 [--dirbuild]
# using pre-specified mat_peptide info:                                            /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg-build-directories/norovirus-builds/NC_001959/NC_001959.dnaorg_build.matpept [--matpept]
# annotate additional qualifiers as duplicates:                                    gene [--dfeat]
# single cmdline argument is a fasta file of sequences, not a list of accessions:  yes [--infasta]
# specify reference accession is <s>:                                              NC_001959 [--refaccn]
# read per-sequence classification errors from <s>:                                dc-noro.9/dc-noro.9.dnaorg_classify.all.errlist [--classerrors]
# run cmscan locally instead of on farm:                                           yes [--local]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Verifying options are consistent with options used for dnaorg_build.pl                ... done. [0.0 seconds]
# Processing input fasta file                                                           ... done. [0.0 seconds]
# Fetching all sequences and processing the reference genome                            ... done. [0.0 seconds]
# Skipping verification that CMs created for current reference NC_001959 (--infasta)    ... done. [0.0 seconds]
# Running cmscan locally                                                                ... done. [1.2 seconds]
# Parsing cmscan results                                                                ... done. [0.0 seconds]
# Calculating predicted feature lengths                                                 ... done. [0.0 seconds]
# Fetching cmscan predicted hits into fasta files                                       ... done. [0.0 seconds]
# Combining predicted exons into CDS                                                    ... done. [0.0 seconds]
# Combining predicted mature peptides into CDS                                          ... done. [0.0 seconds]
# Identifying errors associated with incomplete alignment to the model                  ... done. [0.0 seconds]
# Identifying internal starts/stops in coding sequences                                 ... done. [1.1 seconds]
# Correcting homology search stop codon predictions to account for observed stop codons ... done. [0.0 seconds]
# Identifying overlap and adjacency errors                                              ... done. [0.0 seconds]
# Finalizing annotations and validating error combinations                              ... done. [0.0 seconds]
# Fetching corrected matches into fasta files                                           ... done. [0.0 seconds]
# Combining corrected exons into CDS                                                    ... done. [0.0 seconds]
# Combining corrected mature peptides into CDS                                          ... done. [0.0 seconds]
# Running and parsing BLASTX                                                            ... done. [0.5 seconds]
# Translating corrected nucleotide features into protein sequences                      ... done. [1.1 seconds]
# Generating error code output                                                          ... done. [0.0 seconds]
# Generating tabular annotation output                                                  ... done. [0.0 seconds]
# Generating feature table output                                                       ... done. [0.0 seconds]
#
# Annotated 2 accessions:
#      2 PASS (1.000) listed in dc-noro.9/dc-noro.9-NC_001959/dc-noro.9-NC_001959.dnaorg_annotate.ap.seqlist
#      0 FAIL (0.000) listed in dc-noro.9/dc-noro.9-NC_001959/dc-noro.9-NC_001959.dnaorg_annotate.af.seqlist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

And then lists how many sequences 'PASSed' or 'FAILed'
annotation. (More on this below.)

After dnaorg_annotate.pl is run once per RefSeq (in this case 5
times), dnaorg_classify.pl outputs some summary statistics:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This indicates per-RefSeq, how many sequences were annotated
('num-annot'), passed classification ('num-C-PASS'), failed
classification ('num-C-FAIL'), passed annotation ('num-A-PASS'),
failed annotation ('num-A-FAIL') and failed annotation only because
they failed classification ('num-A-FAIL-CO').

And then a list of output files:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Output printed to screen saved in:                                                                                                                             dc-noro.9.dnaorg_classify.log
# List of executed commands saved in:                                                                                                                            dc-noro.9.dnaorg_classify.cmd
# List and description of all output files saved in:                                                                                                             dc-noro.9.dnaorg_classify.list
# Fasta file with sequences to classify (copy of /panfs/pan1/infernal/notebook/19_0124_virus_dnaorg_release_0p44/dnaorg_scripts/testfiles/noro.9.fa) saved in:   dc-noro.9.dnaorg_classify.fa
# Per-sequence hit and classification information, human readable saved in:                                                                                      dc-noro.9.dnaorg_classify.rdb.infotbl
# Per-sequence hit and classification information, tab-delimited saved in:                                                                                       dc-noro.9.dnaorg_classify.tab.infotbl
# List of errors (unexpected features that cause failure) for all sequences saved in:                                                                            dc-noro.9.dnaorg_classify.all.errlist
# List of errors (unexpected features that cause failure) for sequences that will not be annotated saved in:                                                     dc-noro.9.dnaorg_classify.noannot.errlist
# List of sequences that will not be annotated saved in:                                                                                                         dc-noro.9.dnaorg_classify.noannot.seqlist
# List of PASSing seqs for NC_001959 saved in:                                                                                                                   dc-noro.9.dnaorg_classify.NC_001959.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_001959 saved in:                                                                                                   dc-noro.9.dnaorg_classify.NC_001959.cp.fa
# List of all seqs for NC_001959 saved in:                                                                                                                       dc-noro.9.dnaorg_classify.NC_001959.seqlist
# Fasta file with all seqs assigned to NC_001959 saved in:                                                                                                       dc-noro.9.dnaorg_classify.NC_001959.fa
# List of PASSing seqs for NC_008311 saved in:                                                                                                                   dc-noro.9.dnaorg_classify.NC_008311.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_008311 saved in:                                                                                                   dc-noro.9.dnaorg_classify.NC_008311.cp.fa
# List of all seqs for NC_008311 saved in:                                                                                                                       dc-noro.9.dnaorg_classify.NC_008311.seqlist
# Fasta file with all seqs assigned to NC_008311 saved in:                                                                                                       dc-noro.9.dnaorg_classify.NC_008311.fa
# List of PASSing seqs for NC_029645 saved in:                                                                                                                   dc-noro.9.dnaorg_classify.NC_029645.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_029645 saved in:                                                                                                   dc-noro.9.dnaorg_classify.NC_029645.cp.fa
# List of all seqs for NC_029645 saved in:                                                                                                                       dc-noro.9.dnaorg_classify.NC_029645.seqlist
# Fasta file with all seqs assigned to NC_029645 saved in:                                                                                                       dc-noro.9.dnaorg_classify.NC_029645.fa
# List of PASSing seqs for NC_031324 saved in:                                                                                                                   dc-noro.9.dnaorg_classify.NC_031324.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_031324 saved in:                                                                                                   dc-noro.9.dnaorg_classify.NC_031324.cp.fa
# List of all seqs for NC_031324 saved in:                                                                                                                       dc-noro.9.dnaorg_classify.NC_031324.seqlist
# Fasta file with all seqs assigned to NC_031324 saved in:                                                                                                       dc-noro.9.dnaorg_classify.NC_031324.fa
# List of PASSing seqs for NC_039477 saved in:                                                                                                                   dc-noro.9.dnaorg_classify.NC_039477.cp.seqlist
# Fasta file with PASSing seqs assigned to NC_039477 saved in:                                                                                                   dc-noro.9.dnaorg_classify.NC_039477.cp.fa
# List of all seqs for NC_039477 saved in:                                                                                                                       dc-noro.9.dnaorg_classify.NC_039477.seqlist
# Fasta file with all seqs assigned to NC_039477 saved in:                                                                                                       dc-noro.9.dnaorg_classify.NC_039477.fa
# List of RefSeqs in the HMM library saved in:                                                                                                                   dc-noro.9.dnaorg_classify.all.refseqs
# List of sequences that were sorted into seqlists saved in:                                                                                                     dc-noro.9.dnaorg_classify.all.seqs
# annotation results for NC_001959 sequences that pass dnaorg_annotate.pl saved in:                                                                              dc-noro.9-NC_001959.dnaorg_annotate.ap.sqtable
# annotation results for NC_001959 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    dc-noro.9-NC_001959.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      dc-noro.9-NC_001959.dnaorg_annotate.long.sqtable
# list of NC_001959 sequences that pass dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_001959.dnaorg_annotate.ap.seqlist
# list of NC_001959 sequences that fail dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_001959.dnaorg_annotate.af.seqlist
# list of NC_001959 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            dc-noro.9-NC_001959.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 dc-noro.9-NC_001959.dnaorg_annotate.errlist
# annotation results for NC_008311 sequences that pass dnaorg_annotate.pl saved in:                                                                              dc-noro.9-NC_008311.dnaorg_annotate.ap.sqtable
# annotation results for NC_008311 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    dc-noro.9-NC_008311.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      dc-noro.9-NC_008311.dnaorg_annotate.long.sqtable
# list of NC_008311 sequences that pass dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_008311.dnaorg_annotate.ap.seqlist
# list of NC_008311 sequences that fail dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_008311.dnaorg_annotate.af.seqlist
# list of NC_008311 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            dc-noro.9-NC_008311.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 dc-noro.9-NC_008311.dnaorg_annotate.errlist
# annotation results for NC_029645 sequences that pass dnaorg_annotate.pl saved in:                                                                              dc-noro.9-NC_029645.dnaorg_annotate.ap.sqtable
# annotation results for NC_029645 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    dc-noro.9-NC_029645.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      dc-noro.9-NC_029645.dnaorg_annotate.long.sqtable
# list of NC_029645 sequences that pass dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_029645.dnaorg_annotate.ap.seqlist
# list of NC_029645 sequences that fail dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_029645.dnaorg_annotate.af.seqlist
# list of NC_029645 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            dc-noro.9-NC_029645.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 dc-noro.9-NC_029645.dnaorg_annotate.errlist
# annotation results for NC_031324 sequences that pass dnaorg_annotate.pl saved in:                                                                              dc-noro.9-NC_031324.dnaorg_annotate.ap.sqtable
# annotation results for NC_031324 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    dc-noro.9-NC_031324.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      dc-noro.9-NC_031324.dnaorg_annotate.long.sqtable
# list of NC_031324 sequences that pass dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_031324.dnaorg_annotate.ap.seqlist
# list of NC_031324 sequences that fail dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_031324.dnaorg_annotate.af.seqlist
# list of NC_031324 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            dc-noro.9-NC_031324.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 dc-noro.9-NC_031324.dnaorg_annotate.errlist
# annotation results for NC_039477 sequences that pass dnaorg_annotate.pl saved in:                                                                              dc-noro.9-NC_039477.dnaorg_annotate.ap.sqtable
# annotation results for NC_039477 sequences that fail dnaorg_annotate.pl (minimal) saved in:                                                                    dc-noro.9-NC_039477.dnaorg_annotate.af.sqtable
# annotation results for sequences that pass or fail dnaorg_annotate.pl (verbose) saved in:                                                                      dc-noro.9-NC_039477.dnaorg_annotate.long.sqtable
# list of NC_039477 sequences that pass dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_039477.dnaorg_annotate.ap.seqlist
# list of NC_039477 sequences that fail dnaorg_annotate.pl saved in:                                                                                             dc-noro.9-NC_039477.dnaorg_annotate.af.seqlist
# list of NC_039477 sequences that fail dnaorg_annotate.pl ONLY do to classification errors saved in:                                                            dc-noro.9-NC_039477.dnaorg_annotate.af-co.seqlist
# list of all sequence table errors in (tab-delimited) saved in:                                                                                                 dc-noro.9-NC_039477.dnaorg_annotate.errlist
#
# All output files created in directory ./dc-noro.9/
#
# CPU time:  00:02:44.13
#            hh:mm:ss
# 
# DNAORG-SUCCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All output files will be in the output directory 'dc-noro-9'. The most
important files are:

dc-noro.9-<REFSEQ-ACCN>.dnaorg_annotate.ap.seqlist: list of sequences
that passed for this RefSeq.

dc-noro.9-<REFSEQ-ACCN>.dnaorg_annotate.af.seqlist: list of sequences
that failed for this RefSeq.

dc-noro.9-<REFSEQ-ACCN>.dnaorg_annotate.ap.sqtable: feature table of
annotation for sequences that passed for this RefSeq.

dc-noro.9-<REFSEQ-ACCN>.dnaorg_annotate.af.sqtable: feature table of
annotation for sequences that failed for this RefSeq.

dc-noro.9-<REFSEQ-ACCN>.dnaorg_annotate.ap.sqtable: feature table of
annotation for sequefnces that passed for this RefSeq.

dc-noro.9-<REFSEQ-ACCN>.dnaorg_annotate.af.sqtable: feature table of
annotation for sequences that failed for this RefSeq.

dc-noro.9-<REFSEQ-ACCN>.dnaorg_annotate.ap.sqtable: feature table of
annotation for sequefnces that passed for this RefSeq.

dc-noro.9-<REFSEQ-ACCN>.dnaorg_annotate.errlist: list of all errors
for all failing sequences for this RefSeq.

dc-noro.9.dnaorg_classify.errlist: list of all classification errors

dc-noro.9.dnaorg_classify.noannot.errlist: list of classification
errors for sequences that were not annotated (did not match any
RefSeq or matched best to a non-Norovirus calicivirus).

dc-noro.9.dnaorg_classify.noannot.seqlist: list of sequences that were
not annotated.

Any sequence that includes 1 or more error will FAIL. Sequences with 0
errors PASS.

###############################################
POSSIBLE ERRORS

Below is a list of all possible errors and brief explanations of
each. These will appear in the .errlist and .sqtable files.

The 'per-feature annotation errors' are errors that occur for a
specific feature (e.g. CDS) during the annotation stage
(dnaorg_annotate.pl). 

The 'per-sequence annotation errors' are errors that occur for a
complete sequence during the annotation stage (dnaorg_annotate.pl). 

The 'per-sequence classification errors' are errors that occur for a
complete sequence during the classification stage (beginning stage of
dnaorg_classify.pl). 

(per-feature annotation errors)
Mutation at Start               Sequence contains a mutation at point where start codon for coding region should be located.
Mutation at End                 Sequence contains a mutation at point where stop codon for coding region should be located.
{CDS, Peptide} Has Stop Codon   The predicted coding region contains an internal stop codon.
Reverse Complement              Feature seems to be reverse complemented.
Insertion of Nucleotides        The predicted coding region contains an insertion longer than the maximum allowed (27)
Deletion of Nucleotides         The predicted coding region contains a deletion longer than the maximum allowed (27)
Indefinite Annotation           Sequence similarity is too low to add confident annotation.     
Unexpected Length               The length of the predicted coding region is not a multiple of 3.
Peptide Adjacency Problem       The predicted mature peptides are not perfectly adjacent.
Indefinite Annotation at Start  Sequence similarity is too low near the start to add confident annotation.
Indefinite Annotation at Stop   Sequence similarity is too low near the stop to add confident annotation.
Peptide Translation Problem     The predicted peptide may not be translated due to a problem in its CDS.
Duplicate Feature               There are unexpectedly multiple predictions of the coding region.

(per-sequence annotation errors)
No Feature Annotation           Sequence does not include any detectable features to annotate despite being classified as Norovirus.

(per-sequence classification errors)
No Annotation                   Sequence similarity is too low to add annotation. 
Unexpected Taxonomy             Sequence is not even classified as Norovirus, but rather as another calicivirus.        
Unexpected Classification       Sequence is classified as Norovirus, but to a different genogroup than specified by the user.
Low Coverage                    Sequence alignment to its best RefSeq match covers less than 90% of the sequence.

The "No Annotation" and "Unexpected Taxonomy" errors are special in
that they prevent a sequence from being annotated. These will only be
included in the 'dnaorg_classify.noannot.seqlist' and
'dnaorg_classify.noannot.errlist' files, and will never occur in any
dnaorg_annotate*seqlist or dnaorg_annotate*errlist files. 

Sequences that have "Unexpected Classification" or "Low Coverage"
errors will still be annotated, but will always FAIL annotation
because they have at least 1 error.

"Unexpected Classification" is only a possible error if the --ecall
and --ecmap options were used to dnaorg_classify.pl, as discussed
above. 

######################################################
NAMING CONVENTION FOR SEQUENCES IN FASTA FILES 
Single exon nucleotide sequences are named <s>/<d1>-<d2>: 
Where <s> is the sequence accession, <d1> is the coordinate of the
start, and <d2> is the coordinate of the stop. 

For example:
NC_001346/150-455

Multi-exon nucleotide sequences of <m> exons are similarly named,
but contain <m> sets of <d1>-<d2>, one per exon separated by a
single ','.

For example:
NC_001346/2527-1886,1793-1353

Protein sequences are similarly named with the coordinates
pertaining to nucleotide coordinates the protein was translated
from. Additionally, the string "-translated" is appended to the
name. 

For the two examples above:
NC_001346/150-455-translated
NC_001346/2527-1886,1793-1353-translated

############################################################
FORMAT OF dnaorg_classify.pl's .infotbl OUTPUT FILE 

dnaorg_classify.pl outputs two files that end in .infotbl, one that
ends in .tab.infotbl (file with tab delimited columns, one line per
input sequence) and .rdb.infotbl (same as .tab.infotbl but with spaces
instead of tabs to make the file more human readable). 

Both files begin with comment ('#'-prefixed) lines explaining each
column. 

Here is an example of a rdb.infotbl file:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Explanations of each column:
#  1. sequence:    Accession number of the sequence
#  2. seqlen:      length of this sequence
#  3. topmodel:    model/RefSeq that the sequence was assigned to (max score)
#  4. toptax:      taxonomic group for 'topmodel' from --ecmap ('-' if none)
#  5. score:       (summed) bit score(s) of all hits to 'topmodel'
#  6. sc/nt:       'score' divided by 'qlen'
#  7. E-val:       E-value of best hit to 'topmodel'
#  8. coverage:    the percentage of the sequence that all hits to 'topmodel' cover
#  9. bias:        correction in bits for biased composition sequences, summed for all hits to 'topmodel'
# 10. #hits:       number of hits to 'topmodel'
# 11. strand:      strand of best hit and all considered to 'topmodel'
# 12. scdmodel:    second best model/Refseq (2nd highest score)
# 13. scdtax:      taxonomic group for 'scdmodel' from --ecmap ('-' if none)
# 14. scdiff:      difference in summed bit score b/t 'topmodel' hit(s) and 'scdmodel' hit(s)
# 15. diff/nt:     'scdiff' divided by 'seqlen'
# 16. covdiff:     amount by which the coverage of the 'topmodel' hit(s) is greater than that of the 'scdmodel' hit(s)
# 17. CC:          'confidence class', first letter based on sc/nt: A: if sc/nt >= 0.300, B: if 0.300 > sc/nt >= 0.200, C: if 0.200 > sc_nt
#                  second letter based on diff/nt: A: if diff/nt >= 0.060, B: if 0.060 > diff/nt >= 0.006, C: if 0.006 > diff_nt
# 18. p/f:         'PASS' if sequence passes, 'FAIL' if it fails
# 19. unexpected
#     features:    unexpected features for this sequence
#                  Possible values in unexpected features column:
#                  Low Score:      'sc/nt'   < 0.300 (threshold settable with --lowscthresh)
#                  Very Low Score: 'sc/nt'   < 0.200 (threshold settable with --vlowscthresh)
#                  Low Diff:       'diff/nt' < 0.060 (threshold settable with --lowdiffthresh)
#                  Very Low Diff:  'diff/nt' < 0.006 (threshold settable with --vlowdiffthresh)
#                  Minus Strand:   top hit is on minus strand
#                  High Bias:     'bias' > (0.250 * ('bias' + 'score')) (threshold settable with --biasfract)
########################################################################################################################################
#
#sequence   seqlen  topmodel       toptax    score  sc/nt     E-val  coverage     bias  #hits  strand  scdmodel       scdtax   scdiff  diff/nt  covdiff  CC   p/f  unexpected-features
#=========  ======  =============  ======  =======  =====  ========  ========  =======  =====  ======  =============  ======  =======  =======  =======  ==  ====  ===================
KY887602.1    7547  NC_039477      -        8084.7  1.071         0     1.000     57.4      1    plus  NC_039476      -        1883.9    0.250    0.004  AA  PASS  -
KT818729.1     243  NC_001959      -         119.0  0.490   6.7e-39     0.988      3.7      1    plus  NC_031324      -          17.6    0.072    0.000  AA  PASS  -
EU437710.1     291  NC_001959      -         237.8  0.817   1.3e-74     0.997      3.8      1    plus  NC_031324      -          96.0    0.330    0.007  AA  PASS  -
DQ288307.1    1094  NC_029645      -         941.9  0.861  4.1e-286     0.998     26.0      1    plus  NC_001959      -         666.7    0.609    0.074  AA  PASS  -
AY237464.1     255  NC_039477      -         212.3  0.833   5.6e-67     0.996      0.8      1    plus  NC_039475      -         111.7    0.438    0.012  AA  PASS  -
KF475958.1     275  NC_031324      -         154.3  0.561   1.6e-49     0.993      0.9      1    plus  NC_001959      -           1.1    0.004    0.033  AC  PASS  Very Low Diff[0.004<0.006];
AB713840.1     347  NC_008311      -         311.3  0.897   1.2e-96     0.997     11.8      1    plus  NC_029645      -         245.7    0.708    0.182  AA  PASS  -
JN585032.1     286  NC_029645      -         226.6  0.792   3.4e-71     0.993      8.0      1    plus  NC_031324      -         157.2    0.550    0.003  AA  PASS  -
JN975492.1    7286  NC_008311      -        4625.9  0.635         0     1.000     39.6      1    plus  NC_039475      -        3472.1    0.477    0.079  AA  PASS  -
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

############################################################
GETTING MORE INFORMATION

The file README-ncbi-internal.txt has additional information and
examples. However these additional examples highlight functionality
that is only currently available from within NCBI due to its
dependence on internal databases. 

github has revision information:
https://github.com/nawrockie/dnaorg_scripts

Questions? contact eric.nawrocki@nih.gov
