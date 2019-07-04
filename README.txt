EPN, Wed Jul  3 21:23:36 2019

VADR 0.981 README.txt 

Organization of this file:

INTRODUCTION
INSTALLATION
SETTING UP ENVIRONMENT VARIABLES
VERIFYING SUCCESSFUL INSTALLATION WITH TEST RUNS
USING v-annotate.pl TO CLASSIFY AND ANNOTATE NOROVIRUS SEQUENCES
LIST OF POSSIBLE ALERTS
LIST OF v-annotate.pl OPTIONS     
VADR MODELS
NAMING CONVENTION FOR SEQUENCES IN OUTPUT FASTA FILES 
GETTING MORE INFORMATION

Questions:
email Eric Nawrocki: eric.nawrocki@nih.gov

Note: VADR was previously named dnaorg_scripts (version 0.01 to 0.91)

##############################################################################
INTRODUCTION

This is documentation for VADR: Viral Annotation DefineR, a suite
of tools for classifying and analyzing viral sequences based on viral
reference sequences. (It is not able to detect or analyze novel
viruses.)

authors: Eric Nawrocki and Alejandro Schaffer
contact: eric.nawrocki@nih.gov

Git repository for vadr:
https://github.com/nawrockie/vadr.git

(To be automatically updated about new releases, you can become a
watcher of the repo on github.)

This file gives examples of using the VADR script v-annotate.pl
v-annotate.pl is used to classify a sequence, by determining which in
a set of 'reference' models it is most similar to, and then annotate
that sequence based on that most similar model.  

Another VADR script, v-build.pl, is used to create the models from
NCBI RefSeq sequences or from input multiple sequence alignments,
potentially with secondary structure annotation. v-build.pl stores the
RefSeq feature annotation in the model, and v-annotate.pl maps that
annotation (e.g. CDS coordinates) onto the sequences it annotates.
VADR includes about 100 prebuilt models. See the section "VADR
MODELS/RefSeqs" below for more information on these.

v-annotate.pl identifies unexpected or divergent attributes of the
sequences it annotates (e.g. invalid or early stop codons in CDS
features) and reports them to the user in the form of 'alerts'.
Additional features users want to be annotated with v-annotate.pl can
be added manually to the models after the v-build.pl step. If
secondary structure annotation is used to build the models with
v-build.pl, it will be used to align the sequences to the models
during the annotation phase of v-annotate.pl.

v-annotate.pl determines if each sequence PASSes or FAILs and outputs
tabular information on the annotation as well as Sequin formatted
feature table files. A subset of the alerts reported by v-annotate.pl
cause a sequence to FAIL. A sequence PASSes if 0 alerts from this
subset are reported. 

VADR is used by GenBank staff to evaluate incoming sequence
submissions of some viruses: Norovirus and Dengue at the time of this
writing. Submitted sequences that PASS v-annotate.pl are accepted into
GenBank.

##############################################################################
INSTALLATION

The file 'install.sh' is an executable file for installing VADR
and its dependencies. That file is located online at github here:
https://raw.githubusercontent.com/nawrockie/vadr/master/install.sh

Copy that file into the directory in which you want to install
VADR. You may want to call that 'vadr-install'. Then move into that
directory and run the command './install.sh'. That will create several
directories in the current directory.

'install.sh' will download and install the software packages VADR, and
Infernal, and the required perl module libraries sequip and
Bio-Easel.

When install.sh is finished running it will print important
instructions to the screen that explain how to modify your environment
variables so that you can run the VADR scripts.  These are discussed
further below in the SETTING UP ENVIRONMENT VARIABLES section.

##############################################################################
SETTING UP ENVIRONMENT VARIABLES

As mentioned above, when you run 'install.sh' or any of the optional
install scripts, instructions will be output about how to change your
environment variables so that you can run the VADR scripts. Those
instructions are also included below for reference, but without the
actual path to where you ran install.sh (below it is replaced with
"<full path to directory in which you ran install.sh>")

**********************NOTE TO INTERNAL NCBI USERS***********************
If you are internal to NCBI, VADR is installed in a central location
here: 
/panfs/pan1/dnaorg/virseqannot/code/vadr-install
To set up your environment variables follow the instructions below but
replace:
<full path to directory in which you ran install.sh>
with
/panfs/pan1/dnaorg/virseqannot/code/vadr-install
***********************************************************************

********************************************************
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
export VADRBLASTDIR="/usr/bin"
export PERL5LIB="$VADRSCRIPTSDIR":"$VADRSEQUIPDIR":"$VADRBIOEASELDIR/blib/lib":"$VADRBIOEASELDIR/blib/arch":"$PERL5LIB"
export PATH="$VADRSCRIPTSDIR":"$PATH"

After adding the export lines to your .bashrc file, source that file
to update your current environment with the command:

source ~/.bashrc

---
If you are using the C shell, add the following
lines to the '.cshrc' file in your home directory:

setenv VADRINSTALLDIR <full path to directory in which you ran install.sh>
setenv VADRSCRIPTSDIR "$VADRINSTALLDIR/vadr"
setenv VADRMODELDIR "$VADRINSTALLDIR/vadr-models"
setenv VADRINFERNALDIR "$VADRINSTALLDIR/infernal-dev/src"
setenv VADREASELDIR "$VADRINSTALLDIR/infernal-dev/easel/miniapps"
setenv VADRBIOEASELDIR "$VADRINSTALLDIR/Bio-Easel"
setenv VADRSEQUIPDIR "$VADRINSTALLDIR/sequip"
setenv VADRBLASTDIR "/usr/bin"
setenv PERL5LIB "$VADRSCRIPTSDIR":"$VADRSEQUIPDIR":"$VADRBIOEASELDIR/blib/lib":"$VADRBIOEASELDIR/blib/arch":"$PERL5LIB"
setenv PATH "$VADRSCRIPTSDIR":"$PATH"


After adding the setenv lines to your .cshrc file, source that file
to update your current environment with the command:

source ~/.cshrc

(To determine which shell you use, type: 'echo $SHELL')


********************************************************
--

If you get an error about PERL5LIB being undefined, change the PERL5LIB
line to add to:
export PERL5LIB="$VADRSCRIPTSDIR":"$VADRSEQUIPDIR":"$VADRBIOEASELDIR/blib/lib":"$VADRBIOEASELDIR/blib/arch"
for .bashrc, OR
setenv PERL5LIB "$VADRSCRIPTSDIR":"$VADRSEQUIPDIR":"$VADRBIOEASELDIR/blib/lib":"$VADRBIOEASELDIR/blib/arch"
for .cshrc. And then do
> source ~/.bashrc
or
> source ~/.cshrc
again.

###########################################################################
VERIFYING SUCCESSFUL INSTALLATION WITH TEST RUNS

The VADR package includes some tests you can run to make sure that
your installation was successful and that your environment variables
are set-up correctly.

There are 2 shell scripts for running tests:

     1	do-install-tests-local.sh
     2	do-install-tests-parallel.sh

The VADR scripts can be run locally on your computer or in parallel
on a compute farm. These two test files test each of those modes. 
If you plan to run the scripts at locally at least some of the time,
then run 'do-install-tests-local.sh'. If you plan to run the scripts
on a compute farm at least some of the time, then run
'do-install-tests-parallel.sh' 

These scripts can take up to several minutes to run. Please be patient.
If something goes wrong the 'local' script will exit quickly. If the 
compute farm is busy, the 'parallel' script make take longer as the
relevant jobs wait to run.

Below is example output for do-install-tests-local.sh:
> sh $VADRSCRIPTSDIR/testfiles/do-install-tests-local.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# v-test.pl :: test VADR scripts [TEST SCRIPT]
# VADR 0.92 (May 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Mon May  6 12:12:22 2019
#
# test file:                    /panfs/pan1/dnaorg/virseqannot/code/vadr-install/vadr/testfiles/noro.r10.local.testin
# output directory:             n10-local
# forcing directory overwrite:  yes [-f]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parsing test file                                  ... done. [    0.0 seconds]
# Running command  1 [annotate-noro-10-local]        ... done. [   22.5 seconds]
#	checking va-noro.r10/va-noro.r10.vadr.pass.sqtable                                                            ... pass
#	checking va-noro.r10/va-noro.r10.vadr.fail.sqtable                                                            ... pass
#	checking va-noro.r10/va-noro.r10.vadr.sqa.tbl                                                                 ... pass
#	checking va-noro.r10/va-noro.r10.vadr.sqc.tbl                                                                 ... pass
#	checking va-noro.r10/va-noro.r10.vadr.ftr.tbl                                                                 ... pass
#	checking va-noro.r10/va-noro.r10.vadr.sgm.tbl                                                                 ... pass
#	checking va-noro.r10/va-noro.r10.vadr.mdl.tbl                                                                 ... pass
#	checking va-noro.r10/va-noro.r10.vadr.alt.tbl                                                                 ... pass
#	checking va-noro.r10/va-noro.r10.vadr.alc.tbl                                                                 ... pass
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
# Elapsed time:  00:00:22.92
#                hh:mm:ss
# 
[ok]
# v-test.pl :: test VADR scripts [TEST SCRIPT]
# VADR 0.92 (May 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Mon May  6 12:12:46 2019
#
# test file:                    /panfs/pan1/dnaorg/virseqannot/code/vadr-install/vadr/testfiles/dengue.r5.local.testin
# output directory:             d5-local
# forcing directory overwrite:  yes [-f]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parsing test file                                  ... done. [    0.0 seconds]
# Running command  1 [annotate-noro-10-local]        ... done. [   29.2 seconds]
#	checking va-dengue.r5/va-dengue.r5.vadr.pass.sqtable                                                          ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.fail.sqtable                                                          ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.sqa.tbl                                                               ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.sqc.tbl                                                               ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.ftr.tbl                                                               ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.sgm.tbl                                                               ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.mdl.tbl                                                               ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.alt.tbl                                                               ... pass
#	checking va-dengue.r5/va-dengue.r5.vadr.alc.tbl                                                               ... pass
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
# Elapsed time:  00:00:29.60
#                hh:mm:ss
# 
[ok]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The two most important lines are the lines that begins with "# PASS"

# PASS: all 9 files were created correctly.
# PASS: all 9 files were created correctly.

This means that the test has passed. You should see similar 
lines when you run the other tests. If you do not and need help
figuring out why, email me at eric.nawrocki@nih.gov.

###############################################
USING v-annotate.pl TO CLASSIFY AND ANNOTATE NOROVIRUS SEQUENCES

The v-annotate.pl program can be used to classify and annotate
norovirus sequences. It actually can be used generally for any type of
virus given pre-computed models for that virus. However, it has
currently only been tested extensively for norovirus and the examples
here are limited to norovirus. 

v-annotate.pl takes as input a fasta file and for each sequence
in that file, it classifies the sequence as norovirus by determining
the most similar norovirus RefSeq and then annotates that sequence
based on that RefSeq's annotation. 

The following command provides an example, using the input file 
sequence file noro.9.fa with 9 sequences in the 
$VADRSCRIPTSDIR/testfiles directory.

v-annotate.pl takes two command line arguments, the first is the fasta
file with the sequences to classify and annotate and the second is the
name of the output directory to create and populate with output files:

> v-annotate.pl --group Norovirus $VADRSCRIPTSDIR/testfiles/noro.9.fa va-noro.9

Additionally, if you know what genogroup your sequences are, and you want
v-annotate.pl to enforce that genogroup during its classification
stage (and throw an error for any sequence that seems to not be that
genogroup), you can use the --subgroup <s> option as well, where <s>
is one of the possible Norovirus genogroups: "GI", "GII", "GIII",
"GIV", or "GV", like this:

> v-annotate.pl --group Norovirus --subgroup GII $VADRSCRIPTSDIR/testfiles/noro.9.fa va-noro.9

Here is the output from the above command (without the --subgroup option):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
> v-annotate.pl --group Norovirus $VADRSCRIPTSDIR/testfiles/noro.9.fa va-noro.9

# v-annotate.pl :: classify and annotate sequences using a CM library
# VADR 0.92 (May 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:              Mon May  6 12:27:22 2019
# $VADRBIOEASELDIR:  /panfs/pan1/dnaorg/virseqannot/code/vadr-install/Bio-Easel
# $VADRBLASTDIR:     /usr/bin
# $VADREASELDIR:     /panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal-dev/easel/miniapps
# $VADRINFERNALDIR:  /panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal-dev/src
# $VADRMODELDIR:     /panfs/pan1/dnaorg/virseqannot/code/vadr-install/vadr-models
# $VADRSCRIPTSDIR:   /panfs/pan1/dnaorg/virseqannot/code/vadr-install/vadr
#
# sequence file:                                         /panfs/pan1/dnaorg/virseqannot/code/vadr-install/vadr/testfiles/noro.9.fa
# output directory:                                      va-noro.9
# set expected classification of all seqs to group <s>:  Norovirus [--group]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validating input                                             ... done. [    0.6 seconds]
# Classifying sequences (9 seqs)                               ... done. [   20.9 seconds]
# Determining sequence coverage (NC_001959: 2 seqs)            ... done. [    0.7 seconds]
# Determining sequence coverage (NC_008311: 2 seqs)            ... done. [    3.1 seconds]
# Determining sequence coverage (NC_029645: 2 seqs)            ... done. [    1.2 seconds]
# Determining sequence coverage (NC_031324: 1 seq)             ... done. [    0.8 seconds]
# Determining sequence coverage (NC_039477: 2 seqs)            ... done. [    3.1 seconds]
# Aligning sequences (NC_001959: 2 seqs)                       ... done. [    0.9 seconds]
# Aligning sequences (NC_008311: 2 seqs)                       ... done. [   13.1 seconds]
# Aligning sequences (NC_029645: 2 seqs)                       ... done. [    1.8 seconds]
# Aligning sequences (NC_031324: 1 seq)                        ... done. [    0.7 seconds]
# Aligning sequences (NC_039477: 2 seqs)                       ... done. [   13.8 seconds]
# Determining annotation                                       ... done. [    0.4 seconds]
# Running and parsing BLASTX                                   ... done. [    2.3 seconds]
# Generating tabular output                                    ... done. [    0.0 seconds]
# Generating feature table output                              ... done. [    0.0 seconds]
#
# Summary of classified sequences:
#
#                                      num   num   num       num
#idx  model      group      subgroup  seqs  pass  fail  notannot
#---  ---------  ---------  --------  ----  ----  ----  --------
1     NC_008311  Norovirus  GV           2     1     1         0
2     NC_001959  Norovirus  GI           2     2     0         0
3     NC_039477  Norovirus  GII          2     2     0         0
4     NC_029645  Norovirus  GIII         2     2     0         0
5     NC_031324  Norovirus  GI           1     1     0         0
#---  ---------  ---------  --------  ----  ----  ----  --------
-     *all*      -          -            9     8     1         0
#---  ---------  ---------  --------  ----  ----  ----  --------
#
# Summary of reported alerts:
#
#     alert  causes   short                               per    num   num  long
#idx  code   failure  description                        type  cases  seqs  description
#---  -----  -------  ------------------------------  -------  -----  ----  -----------
1     n_stp  yes      Mutation_at_End                 feature      1     1  expected stop codon could not be identified, predicted CDS stop by homology is invalid
2     n_trc  yes      Unexpected_Stop_Codon           feature      1     1  in-frame stop codon exists 5' of stop position predicted by homology to reference
3     p_trc  yes      Unexpected_Stop_Codon           feature      1     1  stop codon in protein-based alignment
4     b_p5s  yes      Indefinite_Annotation_at_Start  feature      1     1  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint
5     b_p3s  yes      Indefinite_Annotation_at_End    feature      1     1  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
#---  -----  -------  ------------------------------  -------  -----  ----  -----------
#
# Output printed to screen saved in:                                      va-noro.9.vadr.log
# List of executed commands saved in:                                     va-noro.9.vadr.cmd
# List and description of all output files saved in:                      va-noro.9.vadr.list
# predicted hits to model NC_001959 for feature CDS#1 saved in:           va-noro.9.vadr.NC_001959.CDS.1.fa
# predicted hits to model NC_001959 for feature mat_peptide#6 saved in:   va-noro.9.vadr.NC_001959.mat_peptide.6.fa
# predicted hits to model NC_001959 for feature CDS#2 saved in:           va-noro.9.vadr.NC_001959.CDS.2.fa
# predicted hits to model NC_008311 for feature CDS#1 saved in:           va-noro.9.vadr.NC_008311.CDS.1.fa
# predicted hits to model NC_008311 for feature mat_peptide#6 saved in:   va-noro.9.vadr.NC_008311.mat_peptide.6.fa
# predicted hits to model NC_008311 for feature CDS#2 saved in:           va-noro.9.vadr.NC_008311.CDS.2.fa
# predicted hits to model NC_008311 for feature CDS#3 saved in:           va-noro.9.vadr.NC_008311.CDS.3.fa
# predicted hits to model NC_008311 for feature mat_peptide#1 saved in:   va-noro.9.vadr.NC_008311.mat_peptide.1.fa
# predicted hits to model NC_008311 for feature mat_peptide#2 saved in:   va-noro.9.vadr.NC_008311.mat_peptide.2.fa
# predicted hits to model NC_008311 for feature mat_peptide#3 saved in:   va-noro.9.vadr.NC_008311.mat_peptide.3.fa
# predicted hits to model NC_008311 for feature mat_peptide#4 saved in:   va-noro.9.vadr.NC_008311.mat_peptide.4.fa
# predicted hits to model NC_008311 for feature mat_peptide#5 saved in:   va-noro.9.vadr.NC_008311.mat_peptide.5.fa
# predicted hits to model NC_008311 for feature CDS#4 saved in:           va-noro.9.vadr.NC_008311.CDS.4.fa
# predicted hits to model NC_029645 for feature CDS#1 saved in:           va-noro.9.vadr.NC_029645.CDS.1.fa
# predicted hits to model NC_029645 for feature mat_peptide#6 saved in:   va-noro.9.vadr.NC_029645.mat_peptide.6.fa
# predicted hits to model NC_029645 for feature CDS#2 saved in:           va-noro.9.vadr.NC_029645.CDS.2.fa
# predicted hits to model NC_031324 for feature CDS#2 saved in:           va-noro.9.vadr.NC_031324.CDS.2.fa
# predicted hits to model NC_039477 for feature CDS#1 saved in:           va-noro.9.vadr.NC_039477.CDS.1.fa
# predicted hits to model NC_039477 for feature mat_peptide#1 saved in:   va-noro.9.vadr.NC_039477.mat_peptide.1.fa
# predicted hits to model NC_039477 for feature mat_peptide#2 saved in:   va-noro.9.vadr.NC_039477.mat_peptide.2.fa
# predicted hits to model NC_039477 for feature mat_peptide#3 saved in:   va-noro.9.vadr.NC_039477.mat_peptide.3.fa
# predicted hits to model NC_039477 for feature mat_peptide#4 saved in:   va-noro.9.vadr.NC_039477.mat_peptide.4.fa
# predicted hits to model NC_039477 for feature mat_peptide#5 saved in:   va-noro.9.vadr.NC_039477.mat_peptide.5.fa
# predicted hits to model NC_039477 for feature mat_peptide#6 saved in:   va-noro.9.vadr.NC_039477.mat_peptide.6.fa
# predicted hits to model NC_039477 for feature CDS#2 saved in:           va-noro.9.vadr.NC_039477.CDS.2.fa
# predicted hits to model NC_039477 for feature CDS#3 saved in:           va-noro.9.vadr.NC_039477.CDS.3.fa
# per-sequence tabular annotation summary file saved in:                  va-noro.9.vadr.sqa.tbl
# per-sequence tabular classification summary file saved in:              va-noro.9.vadr.sqc.tbl
# per-feature tabular summary file saved in:                              va-noro.9.vadr.ftr.tbl
# per-model-segment tabular summary file saved in:                        va-noro.9.vadr.sgm.tbl
# per-model tabular summary file saved in:                                va-noro.9.vadr.mdl.tbl
# per-alert tabular summary file saved in:                                va-noro.9.vadr.alt.tbl
# alert count tabular summary file saved in:                              va-noro.9.vadr.alc.tbl
# Sequin feature table output for passing sequences saved in:             va-noro.9.vadr.pass.sqtable
# Sequin feature table output for failing sequences saved in:             va-noro.9.vadr.fail.sqtable
# list of passing sequences saved in:                                     va-noro.9.vadr.pass.seqlist
# list of failing sequences saved in:                                     va-noro.9.vadr.fail.seqlist
# list of alerts in the feature tables saved in:                          va-noro.9.vadr.alt.list
#
# All output files created in directory ./va-noro.9/
#
# Elapsed time:  00:01:04.91
#                hh:mm:ss
# 
[ok]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The output lists each step the program performs. The first step is to
validate the input sequences and then to classify them. Next, the
coverage of each sequence is determined and each sequence is aligned
to its best-matching model/RefSeq to determine annotation. These two
steps are done for each model/RefSeq to which one or more sequences
were classified. After the annotation is determined, BLASTX is run to
validate the CDS predictions and output is generated.

The output summarizes how many sequences were classified to each
model/RefSeq and how many PASSed/FAILed, and the number of and types
of alerts printed.

Finally, the progam outputs a list of output files created in the
output directory (va-noro.9), along with brief descriptions.

###############################################
LIST OF POSSIBLE ALERTS

The possible alerts can be output using the --alt_list option to
v-annotate.pl:

> v-annotate.pl --alt_list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################################################
#
# VADR 0.92 (May 2019)
#
# Alert codes that ALWAYS cause a sequence to FAIL, and cannot be
# listed in --alt_pass or --alt_fail option strings:
#
#     alert  short                  long       
#idx  code   description            description
#---  -----  ---------------------  -----------
1     c_noa  No_Annotation          no significant similarity detected
2     c_mst  Minus_Strand           sequence appears to be reverse complemented
3     n_div  Unexpected_Divergence  sequence is too divergent to confidently assign nucleotide-based annotation
4     b_zft  No_Features_Annotated  sequence similarity to homology model does not overlap with any features
#
############################################################
#
# Alert codes that cause a sequence to FAIL by default, but can be set
# to not FAIL a sequence by listing the code as part of a comma separated
# string of codes in <s> with the --alt_pass <s> option:
#
#     alert  short                            long       
#idx  code   description                      description
#---  -----  -------------------------------  -----------
1     c_isg  Incorrect_Specified_Subgroup     score difference too large between best overall model and best specified subgroup model
2     c_igr  Incorrect_Specified_Group        score difference too large between best overall model and best specified group model
3     c_loc  Low_Coverage                     low sequence fraction with significant similarity to homology model
4     c_hbi  Biased_Sequence                  high fraction of score attributed to biased sequence composition
5     c_dpr  Duplicate_Regions                similarity to a model region occurs more than once
6     c_dcs  Discontiguous_Similarity         not all hits are in the same order in the sequence and the homology model
7     c_bst  Indefinite_Strand                significant similarity detected on both strands
8     c_lss  Low_Similarity_at_Start          significant similarity not detected at 5' end of the sequence
9     c_lse  Low_Similarity_at_End            significant similarity not detected at 3' end of the sequence
10    c_lsi  Low_Similarity                   internal region without significant similarity
11    n_str  Mutation_at_Start                expected start codon could not be identified
12    n_stp  Mutation_at_End                  expected stop codon could not be identified, predicted CDS stop by homology is invalid
13    n_nst  Mutation_at_End                  expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted valid start codon
14    n_ext  Mutation_at_End                  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position
15    n_nm3  Unexpected_Length                length of complete coding (CDS or mat_peptide) feature is not a multiple of 3
16    n_trc  Unexpected_Stop_Codon            in-frame stop codon exists 5' of stop position predicted by homology to reference
17    p_trc  Unexpected_Stop_Codon            stop codon in protein-based alignment
18    b_per  Peptide_Translation_Problem      mat_peptide may not be translated because a CDS that spans it has a problem
19    n_adj  Peptide_Adjacency_Problem        predictions of two mat_peptides expected to be adjacent are not adjacent
20    b_non  Indefinite_Annotation            protein-based search identifies CDS not identified in nucleotide-based search
21    b_nop  Indefinite_Annotation            nucleotide-based search identifies CDS not identified in protein-based search
22    n_gp5  Indefinite_Annotation_at_Start   alignment to homology model is a gap at 5' boundary
23    n_lp5  Indefinite_Annotation_at_Start   alignment to homology model has low confidence at 5' boundary
24    b_p5l  Indefinite_Annotation_at_Start   protein-based alignment extends past nucleotide-based alignment at 5' end
25    b_p5s  Indefinite_Annotation_at_Start   protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint
26    n_gp3  Indefinite_Annotation_at_End     alignment to homology model is a gap at 3' boundary
27    n_lp3  Indefinite_Annotation_at_End     alignment to homology model has low confidence at 3' boundary
28    b_p3l  Indefinite_Annotation_at_End     protein-based alignment extends past nucleotide-based alignment at 3' end
29    b_p3s  Indefinite_Annotation_at_End     protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
30    b_cst  Indefinite_Strand                strand mismatch between protein-based and nucleotide-based predictions
31    p_lin  Insertion_of_Nucleotides         too large of an insertion in protein-based alignment
32    p_lde  Deletion_of_Nucleotides          too large of a deletion in protein-based alignment
33    x_fss  Low_Feature_Similarity_at_Start  region within annotated feature at 5' end of sequence lacks significant similarity
34    x_fse  Low_Feature_Similarity_at_End    region within annotated feature at 3' end of sequence lacks significant similarity
35    x_fsi  Low_Feature_Similarity           region within annotated feature lacks significant similarity
#
############################################################
#
# Alert codes that do not cause a sequence to FAIL by default, but can be set
# to FAIL a sequence by listing the code as part of a comma separated
# string of codes in <s> with the --alt_fail <s> option:
#
#     alert  short                            long       
#idx  code   description                      description
#---  -----  -------------------------------  -----------
1     c_qsg  Questionable_Specified_Subgroup  best overall model is not from specified subgroup
2     c_qgr  Questionable_Specified_Group     best overall model is not from specified group
3     c_idc  Indefinite_Classification        low score difference between best overall model and second best model (not in best model's subgroup)
4     c_los  Low_Score                        score to homology model below low threshold

######################################################
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################
LIST OF v-annotate.pl OPTIONS

To get a full list of possible options to VADR scripts, use the -h
option without any other command line arguments. The list of
v-annotate.pl options are below:

> v-annotate.pl -h
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# v-annotate.pl :: classify and annotate sequences using a CM library
# VADR 0.97 (Jun 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Mon Jun 10 11:50:55 2019
#
Usage: v-annotate.pl [-options] <fasta file to annotate> <output directory to create>


basic options:
  -f        : force; if output dir exists, overwrite it
  -v        : be verbose; output commands to stdout as they're run
  -s <n>    : seed for random number generator is <n> [181]
  -m <s>    : use CM file <s> instead of default
  -i <s>    : use model info file <s> instead of default
  -b <s>    : specify BLAST dbs are in dir <s>, instead of default
  -n <n>    : use <n> CPUs [0]
  --atgonly : only consider ATG a valid start codon
  --keep    : do not remove intermediate files, keep them all on disk

options for specifying classification:
  --group <s>    : set expected classification of all seqs to group <s>
  --subgroup <s> : set expected classification of all seqs to subgroup <s>

options for controlling which alerts cause a sequence to FAIL:
  --alt_list     : output summary of all alerts and exit
  --alt_pass <s> : specify that alert codes in comma-separated <s> DO NOT cause FAILure
  --alt_fail <s> : specify that alert codes in comma-separated <s> DO     cause FAILure

options for tuning classification alerts:
  --lowcov <x>      : 'Low Coverage' fractional coverage threshold is <x> [0.9]
  --lowsc <x>       : 'Low Score' bits per nucleotide threshold is <x> [0.3]
  --lowsimterm <n>  : 'Low Similarity at Start/End' minimum length is <n> [10]
  --lowsimint <n>   : 'Low Similarity' (internal) minimum length is <n> [1]
  --indefclass <x>  : 'Indefinite Classification' bits per nucleotide diff threshold is <x> [0.03]
  --biasfract <x>   : 'Biased Sequence' fractional threshold is <x> [0.25]
  --dupreg <n>      : 'Duplicate Regions' minimum model overlap is <n> [20]
  --indefstr <x>    : 'Indefinite Strand' minimum weaker strand bit score is <x> [20]
  --indefann <x>    : 'Indefinite Annotation at Start/End' min allowed post probability is <x> for non-mat_peptide features [0.8]
  --indefann_mp <x> : 'Indefinite Annotation at Start/End' min allowed post probability is <x> for mat_peptide features [0.6]
  --incspec <x>     : 'Incorrect Specified {Sub}Group' bits/nt threshold is <x> [0.2]

options for tuning protein validation with blastx:
  --xminntlen <n>  : min CDS/mat_peptide/gene length for feature table output and blastx analysis is <n> [30]
  --xalntol <n>    : max allowed difference in nucleotides b/t nucleotide and blastx start/end postions is <n> [5]
  --xmaxins <n>    : max allowed nucleotide insertion length in blastx validation is <n> [27]
  --xmaxdel <n>    : max allowed nucleotide deletion length in blastx validation is <n> [27]
  --xlonescore <n> : minimum score for a lone blastx (not supported by a CM hit) to cause an error is <n> [80]
  --xmatrix <s>    : use the matrix <s> with blastx (e.g. BLOSUM45)
  --xdrop <n>      : set the xdrop value for blastx to <n> [25]

options for modifying cmalign runs:
  --mxsize <n> : set max allowed dp matrix size --mxsize value for cmalign calls to <n> Mb [8000]
  --tau <x>    : set the initial tau value for cmalign to <x> [0.001]
  --nofixedtau : do not fix the tau value when running cmalign, allow it to decrease if nec
  --nosub      : use alternative alignment strategy for truncated sequences
  --noglocal   : do not run cmalign in glocal mode (run in local mode)

options related to parallelization on compute farm:
  -p             : parallelize cmscan/cmsearch/cmalign on a compute farm
  -q <s>         : use qsub info file <s> instead of default
  --nkb <n>      : number of KB of sequence for each farm job is <n> [10]
  --wait <n>     : allow <n> wall-clock minutes for jobs on farm to finish, including queueing time [500]
  --errcheck     : consider any farm stderr output as indicating a job failure
  --maxnjobs <n> : set max number of jobs to submit to compute farm to <n> [2500]

optional output files:
  --ftrinfo : create file with internal feature information
  --sgminfo : create file with internal segment information
  --altinfo : create file with internal error information

options for skipping stages and using files from earlier, identical runs, primarily useful for debugging:
  --skipalign : skip the cmalign step, use results from an earlier run of the script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################################################
VADR MODELS

VADR includes a set of about 100 viral RefSeq "models". These models
are in the $VADRMODELDIR directory. The vadr.minfo file in that
directory has information on the models. To get information on those
models, grep the "MODEL" lines from that file, like this:

> grep MODEL $VADRMODELDIR/vadr.minfo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MODEL NC_001959 blastdb:"NC_001959.vadr.protein.fa" cmfile:"NC_001959.vadr.cm" group:"Norovirus" length:"7654" subgroup:"GI"
MODEL NC_008311 blastdb:"NC_008311.vadr.protein.fa" cmfile:"NC_008311.vadr.cm" group:"Norovirus" length:"7382" subgroup:"GV"
MODEL NC_029645 blastdb:"NC_029645.vadr.protein.fa" cmfile:"NC_029645.vadr.cm" group:"Norovirus" length:"7313" subgroup:"GIII"
MODEL NC_029646 blastdb:"NC_029646.vadr.protein.fa" cmfile:"NC_029646.vadr.cm" group:"Norovirus" length:"7518" subgroup:"GII"
MODEL NC_029647 blastdb:"NC_029647.vadr.protein.fa" cmfile:"NC_029647.vadr.cm" group:"Norovirus" length:"7527" subgroup:"GIV"
MODEL NC_031324 blastdb:"NC_031324.vadr.protein.fa" cmfile:"NC_031324.vadr.cm" group:"Norovirus" length:"7753" subgroup:"GI"
MODEL NC_039475 blastdb:"NC_039475.vadr.protein.fa" cmfile:"NC_039475.vadr.cm" group:"Norovirus" length:"7556" subgroup:"GII"
MODEL NC_039476 blastdb:"NC_039476.vadr.protein.fa" cmfile:"NC_039476.vadr.cm" group:"Norovirus" length:"7536" subgroup:"GII"
MODEL NC_039477 blastdb:"NC_039477.vadr.protein.fa" cmfile:"NC_039477.vadr.cm" group:"Norovirus" length:"7567" subgroup:"GII"
MODEL NC_001477 blastdb:"NC_001477.vadr.protein.fa" cmfile:"NC_001477.vadr.cm" group:"Dengue" length:"10735" subgroup:"1"
MODEL NC_001474 blastdb:"NC_001474.vadr.protein.fa" cmfile:"NC_001474.vadr.cm" group:"Dengue" length:"10723" subgroup:"2"
MODEL NC_001475 blastdb:"NC_001475.vadr.protein.fa" cmfile:"NC_001475.vadr.cm" group:"Dengue" length:"10707" subgroup:"3"
MODEL NC_002640 blastdb:"NC_002640.vadr.protein.fa" cmfile:"NC_002640.vadr.cm" group:"Dengue" length:"10649" subgroup:"4"
MODEL NC_004102 blastdb:"NC_004102.vadr.protein.fa" cmfile:"NC_004102.vadr.cm" group:"HCV" length:"9646" subgroup:"1"
MODEL NC_038882 blastdb:"NC_038882.vadr.protein.fa" cmfile:"NC_038882.vadr.cm" group:"HCV" length:"9599" subgroup:"1"
MODEL NC_009823 blastdb:"NC_009823.vadr.protein.fa" cmfile:"NC_009823.vadr.cm" group:"HCV" length:"9711" subgroup:"2"
MODEL NC_009824 blastdb:"NC_009824.vadr.protein.fa" cmfile:"NC_009824.vadr.cm" group:"HCV" length:"9456" subgroup:"3"
MODEL NC_009825 blastdb:"NC_009825.vadr.protein.fa" cmfile:"NC_009825.vadr.cm" group:"HCV" length:"9355" subgroup:"4"
MODEL NC_009826 blastdb:"NC_009826.vadr.protein.fa" cmfile:"NC_009826.vadr.cm" group:"HCV" length:"9343" subgroup:"5"
MODEL NC_009827 blastdb:"NC_009827.vadr.protein.fa" cmfile:"NC_009827.vadr.cm" group:"HCV" length:"9628" subgroup:"6"
MODEL NC_030791 blastdb:"NC_030791.vadr.protein.fa" cmfile:"NC_030791.vadr.cm" group:"HCV" length:"9443" subgroup:"7"
MODEL NC_001802 blastdb:"NC_001802.vadr.protein.fa" cmfile:"NC_001802.vadr.cm" group:"HIV" length:"9181" subgroup:"1"
MODEL NC_000940 blastdb:"NC_000940.vadr.protein.fa" cmfile:"NC_000940.vadr.cm" length:"7320"
MODEL NC_000943 blastdb:"NC_000943.vadr.protein.fa" cmfile:"NC_000943.vadr.cm" length:"11014"
MODEL NC_001461 blastdb:"NC_001461.vadr.protein.fa" cmfile:"NC_001461.vadr.cm" length:"12573"
MODEL NC_001481 blastdb:"NC_001481.vadr.protein.fa" cmfile:"NC_001481.vadr.cm" length:"7683"
MODEL NC_001543 blastdb:"NC_001543.vadr.protein.fa" cmfile:"NC_001543.vadr.cm" length:"7437"
MODEL NC_001563 blastdb:"NC_001563.vadr.protein.fa" cmfile:"NC_001563.vadr.cm" length:"10962"
MODEL NC_001564 blastdb:"NC_001564.vadr.protein.fa" cmfile:"NC_001564.vadr.cm" length:"10682"
MODEL NC_001655 blastdb:"NC_001655.vadr.protein.fa" cmfile:"NC_001655.vadr.cm" length:"9399"
MODEL NC_001672 blastdb:"NC_001672.vadr.protein.fa" cmfile:"NC_001672.vadr.cm" length:"11141"
MODEL NC_001710 blastdb:"NC_001710.vadr.protein.fa" cmfile:"NC_001710.vadr.cm" length:"9392"
MODEL NC_001809 blastdb:"NC_001809.vadr.protein.fa" cmfile:"NC_001809.vadr.cm" length:"10871"
MODEL NC_001837 blastdb:"NC_001837.vadr.protein.fa" cmfile:"NC_001837.vadr.cm" length:"9550"
MODEL NC_002031 blastdb:"NC_002031.vadr.protein.fa" cmfile:"NC_002031.vadr.cm" length:"10862"
MODEL NC_002551 blastdb:"NC_002551.vadr.protein.fa" cmfile:"NC_002551.vadr.cm" length:"8284"
MODEL NC_002615 blastdb:"NC_002615.vadr.protein.fa" cmfile:"NC_002615.vadr.cm" length:"7442"
MODEL NC_002657 blastdb:"NC_002657.vadr.protein.fa" cmfile:"NC_002657.vadr.cm" length:"12301"
MODEL NC_003635 blastdb:"NC_003635.vadr.protein.fa" cmfile:"NC_003635.vadr.cm" length:"10600"
MODEL NC_003678 blastdb:"NC_003678.vadr.protein.fa" cmfile:"NC_003678.vadr.cm" length:"12602"
MODEL NC_003679 blastdb:"NC_003679.vadr.protein.fa" cmfile:"NC_003679.vadr.cm" length:"12333"
MODEL NC_003687 blastdb:"NC_003687.vadr.protein.fa" cmfile:"NC_003687.vadr.cm" length:"10839"
MODEL NC_003690 blastdb:"NC_003690.vadr.protein.fa" cmfile:"NC_003690.vadr.cm" length:"10943"
MODEL NC_004064 blastdb:"NC_004064.vadr.protein.fa" cmfile:"NC_004064.vadr.cm" length:"7453"
MODEL NC_004119 blastdb:"NC_004119.vadr.protein.fa" cmfile:"NC_004119.vadr.cm" length:"10690"
MODEL NC_004541 blastdb:"NC_004541.vadr.protein.fa" cmfile:"NC_004541.vadr.cm" length:"8289"
MODEL NC_004542 blastdb:"NC_004542.vadr.protein.fa" cmfile:"NC_004542.vadr.cm" length:"8513"
MODEL NC_005039 blastdb:"NC_005039.vadr.protein.fa" cmfile:"NC_005039.vadr.cm" length:"10857"
MODEL NC_005062 blastdb:"NC_005062.vadr.protein.fa" cmfile:"NC_005062.vadr.cm" length:"10787"
MODEL NC_005064 blastdb:"NC_005064.vadr.protein.fa" cmfile:"NC_005064.vadr.cm" length:"11375"
MODEL NC_006269 blastdb:"NC_006269.vadr.protein.fa" cmfile:"NC_006269.vadr.cm" length:"7429"
MODEL NC_006554 blastdb:"NC_006554.vadr.protein.fa" cmfile:"NC_006554.vadr.cm" length:"7476"
MODEL NC_006875 blastdb:"NC_006875.vadr.protein.fa" cmfile:"NC_006875.vadr.cm" length:"7453"
MODEL NC_007580 blastdb:"NC_007580.vadr.protein.fa" cmfile:"NC_007580.vadr.cm" length:"10940"
MODEL NC_007916 blastdb:"NC_007916.vadr.protein.fa" cmfile:"NC_007916.vadr.cm" length:"7454"
MODEL NC_008580 blastdb:"NC_008580.vadr.protein.fa" cmfile:"NC_008580.vadr.cm" length:"8380"
MODEL NC_008604 blastdb:"NC_008604.vadr.protein.fa" cmfile:"NC_008604.vadr.cm" length:"10837"
MODEL NC_008718 blastdb:"NC_008718.vadr.protein.fa" cmfile:"NC_008718.vadr.cm" length:"10510"
MODEL NC_008719 blastdb:"NC_008719.vadr.protein.fa" cmfile:"NC_008719.vadr.cm" length:"10793"
MODEL NC_009026 blastdb:"NC_009026.vadr.protein.fa" cmfile:"NC_009026.vadr.cm" length:"10815"
MODEL NC_009028 blastdb:"NC_009028.vadr.protein.fa" cmfile:"NC_009028.vadr.cm" length:"10755"
MODEL NC_009029 blastdb:"NC_009029.vadr.protein.fa" cmfile:"NC_009029.vadr.cm" length:"10874"
MODEL NC_009942 blastdb:"NC_009942.vadr.protein.fa" cmfile:"NC_009942.vadr.cm" length:"11029"
MODEL NC_010624 blastdb:"NC_010624.vadr.protein.fa" cmfile:"NC_010624.vadr.cm" length:"7458"
MODEL NC_011050 blastdb:"NC_011050.vadr.protein.fa" cmfile:"NC_011050.vadr.cm" length:"8305"
MODEL NC_011704 blastdb:"NC_011704.vadr.protein.fa" cmfile:"NC_011704.vadr.cm" length:"7422"
MODEL NC_012532 blastdb:"NC_012532.vadr.protein.fa" cmfile:"NC_012532.vadr.cm" length:"10794"
MODEL NC_012533 blastdb:"NC_012533.vadr.protein.fa" cmfile:"NC_012533.vadr.cm" length:"10723"
MODEL NC_012534 blastdb:"NC_012534.vadr.protein.fa" cmfile:"NC_012534.vadr.cm" length:"10941"
MODEL NC_012671 blastdb:"NC_012671.vadr.protein.fa" cmfile:"NC_012671.vadr.cm" length:"10865"
MODEL NC_012699 blastdb:"NC_012699.vadr.protein.fa" cmfile:"NC_012699.vadr.cm" length:"6434"
MODEL NC_012735 blastdb:"NC_012735.vadr.protein.fa" cmfile:"NC_012735.vadr.cm" length:"10814"
MODEL NC_012812 blastdb:"NC_012812.vadr.protein.fa" cmfile:"NC_012812.vadr.cm" length:"12337"
MODEL NC_012932 blastdb:"NC_012932.vadr.protein.fa" cmfile:"NC_012932.vadr.cm" length:"11064"
MODEL NC_016997 blastdb:"NC_016997.vadr.protein.fa" cmfile:"NC_016997.vadr.cm" length:"10791"
MODEL NC_017086 blastdb:"NC_017086.vadr.protein.fa" cmfile:"NC_017086.vadr.cm" length:"10733"
MODEL NC_017936 blastdb:"NC_017936.vadr.protein.fa" cmfile:"NC_017936.vadr.cm" length:"7696"
MODEL NC_018713 blastdb:"NC_018713.vadr.protein.fa" cmfile:"NC_018713.vadr.cm" length:"12292"
MODEL NC_019712 blastdb:"NC_019712.vadr.protein.fa" cmfile:"NC_019712.vadr.cm" length:"8427"
MODEL NC_020902 blastdb:"NC_020902.vadr.protein.fa" cmfile:"NC_020902.vadr.cm" length:"11197"
MODEL NC_021069 blastdb:"NC_021069.vadr.protein.fa" cmfile:"NC_021069.vadr.cm" length:"10865"
MODEL NC_021153 blastdb:"NC_021153.vadr.protein.fa" cmfile:"NC_021153.vadr.cm" length:"8879"
MODEL NC_021154 blastdb:"NC_021154.vadr.protein.fa" cmfile:"NC_021154.vadr.cm" length:"11279"
MODEL NC_023176 blastdb:"NC_023176.vadr.protein.fa" cmfile:"NC_023176.vadr.cm" length:"12656"
MODEL NC_024017 blastdb:"NC_024017.vadr.protein.fa" cmfile:"NC_024017.vadr.cm" length:"10891"
MODEL NC_024018 blastdb:"NC_024018.vadr.protein.fa" cmfile:"NC_024018.vadr.cm" length:"12273"
MODEL NC_024031 blastdb:"NC_024031.vadr.protein.fa" cmfile:"NC_024031.vadr.cm" length:"7443"
MODEL NC_024078 blastdb:"NC_024078.vadr.protein.fa" cmfile:"NC_024078.vadr.cm" length:"8013"
MODEL NC_024299 blastdb:"NC_024299.vadr.protein.fa" cmfile:"NC_024299.vadr.cm" length:"10878"
MODEL NC_025676 blastdb:"NC_025676.vadr.protein.fa" cmfile:"NC_025676.vadr.cm" length:"8489"
MODEL NC_026620 blastdb:"NC_026620.vadr.protein.fa" cmfile:"NC_026620.vadr.cm" length:"10125"
MODEL NC_026624 blastdb:"NC_026624.vadr.protein.fa" cmfile:"NC_026624.vadr.cm" length:"10242"
MODEL NC_027026 blastdb:"NC_027026.vadr.protein.fa" cmfile:"NC_027026.vadr.cm" length:"7521"
MODEL NC_027122 blastdb:"NC_027122.vadr.protein.fa" cmfile:"NC_027122.vadr.cm" length:"8376"
MODEL NC_027709 blastdb:"NC_027709.vadr.protein.fa" cmfile:"NC_027709.vadr.cm" length:"10870"
MODEL NC_027817 blastdb:"NC_027817.vadr.protein.fa" cmfile:"NC_027817.vadr.cm" length:"10893"
MODEL NC_027819 blastdb:"NC_027819.vadr.protein.fa" cmfile:"NC_027819.vadr.cm" length:"10938"
MODEL NC_027999 blastdb:"NC_027999.vadr.protein.fa" cmfile:"NC_027999.vadr.cm" length:"10761"
MODEL NC_029055 blastdb:"NC_029055.vadr.protein.fa" cmfile:"NC_029055.vadr.cm" length:"10290"
MODEL NC_030290 blastdb:"NC_030290.vadr.protein.fa" cmfile:"NC_030290.vadr.cm" length:"10864"
MODEL NC_030793 blastdb:"NC_030793.vadr.protein.fa" cmfile:"NC_030793.vadr.cm" length:"7484"
MODEL NC_033081 blastdb:"NC_033081.vadr.protein.fa" cmfile:"NC_033081.vadr.cm" length:"8176"
MODEL NC_033776 blastdb:"NC_033776.vadr.protein.fa" cmfile:"NC_033776.vadr.cm" length:"7475"
MODEL NC_034007 blastdb:"NC_034007.vadr.protein.fa" cmfile:"NC_034007.vadr.cm" length:"10131"
MODEL NC_034017 blastdb:"NC_034017.vadr.protein.fa" cmfile:"NC_034017.vadr.cm" length:"10884"
MODEL NC_034151 blastdb:"NC_034151.vadr.protein.fa" cmfile:"NC_034151.vadr.cm" length:"10937"
MODEL NC_034204 blastdb:"NC_034204.vadr.protein.fa" cmfile:"NC_034204.vadr.cm" length:"10897"
MODEL NC_034444 blastdb:"NC_034444.vadr.protein.fa" cmfile:"NC_034444.vadr.cm" length:"8399"
MODEL NC_035118 blastdb:"NC_035118.vadr.protein.fa" cmfile:"NC_035118.vadr.cm" length:"10646"
MODEL NC_035432 blastdb:"NC_035432.vadr.protein.fa" cmfile:"NC_035432.vadr.cm" length:"12614"
MODEL NC_035675 blastdb:"NC_035675.vadr.protein.fa" cmfile:"NC_035675.vadr.cm" length:"6564"
MODEL NC_035889 blastdb:"NC_035889.vadr.protein.fa" cmfile:"NC_035889.vadr.cm" length:"10808"
MODEL NC_038425 blastdb:"NC_038425.vadr.protein.fa" cmfile:"NC_038425.vadr.cm" length:"9531"
MODEL NC_038433 blastdb:"NC_038433.vadr.protein.fa" cmfile:"NC_038433.vadr.cm" length:"10479"
MODEL NC_038912 blastdb:"NC_038912.vadr.protein.fa" cmfile:"NC_038912.vadr.cm" length:"12298"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################################################
NAMING CONVENTION FOR SEQUENCES IN OUTPUT FASTA FILES 

Predicted single-segment feature sequences are named
<s1>/<s2>.<d1>/<d2>..<d3>:<s3>

Where <s1> is the sequence name from the input fasta file, <s2> is the
feature type, <d1> is the index for that feature, <d2> and <d3> are
the start and stop nucleotide coordinates, and <s3> is "+" for
positive strand, or "-" for negatiave strand.

For example:
DQ288307.1/CDS.2/197..1094:+

Multi-segment feature sequences of <m> segments are similarly named,
but contain <m> sets of <d2>..<d3>:<s3>, one per segment separated by a
single ','.

For example:
AY940619.1/CDS.2/200..227:+,229..686:+

############################################################
GETTING MORE INFORMATION

The file README-ncbi-internal.txt has additional information and
examples. However these additional examples highlight functionality
that is only currently available from within NCBI due to its
dependence on internal databases. 

github has revision information:
https://github.com/nawrockie/vadr

Questions? contact eric.nawrocki@nih.gov
