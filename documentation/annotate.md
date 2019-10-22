# `v-annotate.pl` usage and command-line options

* [`v-annotate.pl` example usage](#exampleusage)
* [`v-annotate.pl` command-line options](#options)

## `v-annotate.pl` example usage<a name="exampleusage"></a>

`v-annotate.pl` uses previously created VADR models from `v-build.pl`
and uses them to analyze and annotate sequences in an input sequence
file. As part of the analysis of the sequences, more than 40 types of
unexpected characteristics, or *alerts* are detected and reported in
the output. Most of the alerts are *fatal* in that if a sequence has
one or more fatal alerts, they will be designated as *failing*
sequences. Sequences with zero fatal alerts are designated as
*passing* sequences. The types of alerts are described further below.

To determine the command-line usage of 
`v-annotate.pl` (or any VADR script), use the `-h` option, like this:

```
v-annotate.pl -h 
```

You'll see something like the following output:
```
# v-annotate.pl :: classify and annotate sequences using a CM library
# VADR 0.991 (Aug 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Mon Oct 21 15:14:21 2019
#
Usage: v-annotate.pl [-options] <fasta file to annotate> <output directory to create>
```

The first few lines are the banner which show the name of the VADR
script being run along with the version and release date. This is
followed by the time and date the command was executed.  The `Usage:`
line details the expected command line arguments.  `v-annotate.pl` takes
as input two command line arguments, a fasta file with sequences to 
analyze and annotate (`<fasta file to annotate>`) and and the name of the output directory you
want it to create (`<output directory to create>`) and
populate with output files.

After that comes a list of all available command-line options. These
are explained in more detail [below](#options).

Here is an example `v-annotate.pl` command using the fasta file
[vadr/documentation/annotate-files/noro.9.fa](annotate-files/noro.9.fa)
with 9 norovirus sequences, and creating an output directory with
`va-noro.9`:

```
v-annotate.pl $VADRSCRIPTSDIR/documentation/annotate-files/noro.9.fa va-noro.9
```
The standard output of `v-annotate.pl` that is printed to the screen
(which is also output to the `.log` output file) begins with the
banner and date again followed by a list of relevant environment
variables, the command line arguments used and any command line
options used:

```
# date:              Mon Oct 21 15:41:06 2019
# $VADRBIOEASELDIR:  /home/nawrocki/vadr-install/Bio-Easel
# $VADRBLASTDIR:     /usr/bin
# $VADREASELDIR:     /home/nawrocki/vadr-install/infernal-dev/easel/miniapps
# $VADRINFERNALDIR:  /home/nawrocki/vadr-install/infernal-dev/src
# $VADRMODELDIR:     /home/nawrocki/vadr-install/vadr-models
# $VADRSCRIPTSDIR:   /home/nawrocki/vadr-install/vadr
#
# sequence file:     ../vadr/documentation/annotate-files/noro.9.fa
# output directory:  va-noro.9
```

No command line options were used in our example output, but if they
were information on them would have appeared after the `output
directory` line.

Next, information is output about each step the script is proceeding
through. When each step is completed, the elapsed time for that step
is output. 

`v-annotate.pl` will use the default VADR model library (CM file
`$VADRMODELDIR/vadr.cm`, model info file `$VADRMODELDIR/vadr.minfo`
and BLAST DBs in $VADRMODELDIR) to analyze the sequences in
`noro.9.fa`, and will create an output directory named `va-noro.9` and
populate it with [many output files](formats.md#annotate). 

There are four main stages to `v-annotate.pl`

1. Classification: each sequence `S` in the input file is compared
against all models in the model library to determine the best-matching
(highest-scoring) model `M(S)` for each sequence.

2. Coverage determination: each sequence `S` is compared again to
`M(S)` to determine the *coverage* of `S`, which is basically the
fraction of nucleotides in `S` that appear to be homologous
(sufficiently similar) to `M(S)`.

3. Alignment and feature mapping: each sequence `S` is aligned to
`M(S)` and based on that alignment, features of `M(S)` are mapped onto
`S`.

4. Protein validation of CDS features: for each sequence `S` that has 1 or
more predicted CDS features, `blastx` is used to compare the predicted
CDS and the full sequence `S` against the VADR library BLAST DB. 

The output of `v-annotate.pl` lists one or more steps per stage. The
first two steps are:

```
# Validating input                                                                    ... done. [   16.5 seconds]
# Classifying sequences (9 seqs)                                                      ... done. [   37.8 seconds]
```

The first step validates that the VADR library `.minfo` file being
used corresponds to the VADR library `.cm` file. Next, the
classification stage is performed. After that, each model that is 
`M(S)` for at least one `S` is run separately through the coverage
determination stage for all of its sequences:

```
# Determining sequence coverage (NC_001959: 2 seqs)                                   ... done. [    1.6 seconds]
# Determining sequence coverage (NC_008311: 2 seqs)                                   ... done. [    4.8 seconds]
# Determining sequence coverage (NC_029645: 2 seqs)                                   ... done. [    1.7 seconds]
# Determining sequence coverage (NC_031324: 1 seq)                                    ... done. [    1.4 seconds]
# Determining sequence coverage (NC_039477: 2 seqs)                                   ... done. [    5.0 seconds]
```

Next, the alignments are performed for each model, and used to map
feature annotation:

```
# Aligning sequences (NC_001959: 2 seqs)                                              ... done. [    1.7 seconds]
# Aligning sequences (NC_008311: 2 seqs)                                              ... done. [   23.0 seconds]
# Aligning sequences (NC_029645: 2 seqs)                                              ... done. [    2.3 seconds]
# Aligning sequences (NC_031324: 1 seq)                                               ... done. [    1.2 seconds]
# Aligning sequences (NC_039477: 2 seqs)                                              ... done. [   25.3 seconds]
# Determining annotation                                                              ... done. [    0.8 seconds]
```

The classification and alignment stages are typically the
slowest. The protein validation stage is usually relatively fast:

```
# Running and parsing BLASTX                                                          ... done. [    4.1 seconds]
```

The only remaining steps are to create the output files:

```
# Generating tabular output                                                           ... done. [    0.0 seconds]
# Generating feature table output                                                     ... done. [    0.0 seconds]
```

After the output files are generated, a summary of the results is
output, starting with a list of all models that were 
the best matching model to one or more sequences: 

```
# Summary of classified sequences:
#
#                                      num   num   num
#idx  model      group      subgroup  seqs  pass  fail
#---  ---------  ---------  --------  ----  ----  ----
1     NC_008311  Norovirus  GV           2     1     1
2     NC_001959  Norovirus  GI           2     2     0
3     NC_039477  Norovirus  GII          2     2     0
4     NC_029645  Norovirus  GIII         2     2     0
5     NC_031324  Norovirus  GI           1     1     0
#---  ---------  ---------  --------  ----  ----  ----
-     *all*      -          -            9     8     1
-     *none*     -          -            0     0     0
#---  ---------  ---------  --------  ----  ----  ----
```

The above table is also saved to a file with the suffix `.mdl`,
and the format is described in more detail [here](formats.md#mdl).

The nine input sequences collectively had five different best-matching
models, all Norovirus models with the subgroups shown above. (The
groups and subgroups shown above were input as [command-line
options](build.md#options-attributes) to `v-build.pl` when the models
were built. An example is [here](build.md#1.0library-noro).  Eight of
the nine sequences *passed* and one *failed*.

Next, the alerts reported for the sequences are summarized:

```
# Summary of reported alerts:
#
#     alert     causes   short                            per    num   num  long
#idx  code      failure  description                     type  cases  seqs  description
#---  --------  -------  ---------------------------  -------  -----  ----  -----------
1     mutendcd  yes      MUTATION_AT_END              feature      1     1  expected stop codon could not be identified, predicted CDS stop by homology is invalid
2     cdsstopn  yes      CDS_HAS_STOP_CODON           feature      1     1  in-frame stop codon exists 5' of stop position predicted by homology to reference
3     cdsstopp  yes      CDS_HAS_STOP_CODON           feature      1     1  stop codon in protein-based alignment
4     indf5pst  yes      INDEFINITE_ANNOTATION_START  feature      1     1  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint
5     indf3pst  yes      INDEFINITE_ANNOTATION_END    feature      1     1  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
#---  --------  -------  ---------------------------  -------  -----  ----  -----------
```

For the nine sequences, there were five different alert codes
reported, each with one case for a single sequence. This tabular
summary does not indicate which sequences the alerts were reported
for. We can find that out from other output files as discussed further
below.

Next, the list of output files created by `v-annotate.pl` is
printed, along with elapsed time:

```
# Output printed to screen saved in:                               va-noro.9.vadr.log
# List of executed commands saved in:                              va-noro.9.vadr.cmd
# List and description of all output files saved in:               va-noro.9.vadr.filelist
# esl-seqstat -a output for input fasta file saved in:             va-noro.9.vadr.seqstat
# model NC_001959 feature gene#1 predicted seqs saved in:          va-noro.9.vadr.NC_001959.gene.1.fa
# model NC_001959 feature CDS#1 predicted seqs saved in:           va-noro.9.vadr.NC_001959.CDS.1.fa
# model NC_001959 feature mat_peptide#6 predicted seqs saved in:   va-noro.9.vadr.NC_001959.mat_peptide.6.fa
# model NC_001959 feature gene#2 predicted seqs saved in:          va-noro.9.vadr.NC_001959.gene.2.fa
# model NC_001959 feature CDS#2 predicted seqs saved in:           va-noro.9.vadr.NC_001959.CDS.2.fa
# model NC_008311 feature gene#1 predicted seqs saved in:          va-noro.9.vadr.NC_008311.gene.1.fa
# model NC_008311 feature CDS#1 predicted seqs saved in:           va-noro.9.vadr.NC_008311.CDS.1.fa
# model NC_008311 feature gene#2 predicted seqs saved in:          va-noro.9.vadr.NC_008311.gene.2.fa
# model NC_008311 feature CDS#2 predicted seqs saved in:           va-noro.9.vadr.NC_008311.CDS.2.fa
# model NC_008311 feature gene#3 predicted seqs saved in:          va-noro.9.vadr.NC_008311.gene.3.fa
# model NC_008311 feature CDS#3 predicted seqs saved in:           va-noro.9.vadr.NC_008311.CDS.3.fa
# model NC_008311 feature mat_peptide#6 predicted seqs saved in:   va-noro.9.vadr.NC_008311.mat_peptide.6.fa
# model NC_008311 feature gene#4 predicted seqs saved in:          va-noro.9.vadr.NC_008311.gene.4.fa
# model NC_008311 feature CDS#4 predicted seqs saved in:           va-noro.9.vadr.NC_008311.CDS.4.fa
# model NC_008311 feature mat_peptide#1 predicted seqs saved in:   va-noro.9.vadr.NC_008311.mat_peptide.1.fa
# model NC_008311 feature mat_peptide#2 predicted seqs saved in:   va-noro.9.vadr.NC_008311.mat_peptide.2.fa
# model NC_008311 feature mat_peptide#3 predicted seqs saved in:   va-noro.9.vadr.NC_008311.mat_peptide.3.fa
# model NC_008311 feature mat_peptide#4 predicted seqs saved in:   va-noro.9.vadr.NC_008311.mat_peptide.4.fa
# model NC_008311 feature mat_peptide#5 predicted seqs saved in:   va-noro.9.vadr.NC_008311.mat_peptide.5.fa
# model NC_029645 feature gene#1 predicted seqs saved in:          va-noro.9.vadr.NC_029645.gene.1.fa
# model NC_029645 feature CDS#1 predicted seqs saved in:           va-noro.9.vadr.NC_029645.CDS.1.fa
# model NC_029645 feature gene#2 predicted seqs saved in:          va-noro.9.vadr.NC_029645.gene.2.fa
# model NC_029645 feature CDS#2 predicted seqs saved in:           va-noro.9.vadr.NC_029645.CDS.2.fa
# model NC_029645 feature mat_peptide#6 predicted seqs saved in:   va-noro.9.vadr.NC_029645.mat_peptide.6.fa
# model NC_031324 feature gene#2 predicted seqs saved in:          va-noro.9.vadr.NC_031324.gene.2.fa
# model NC_031324 feature CDS#2 predicted seqs saved in:           va-noro.9.vadr.NC_031324.CDS.2.fa
# model NC_039477 feature gene#1 predicted seqs saved in:          va-noro.9.vadr.NC_039477.gene.1.fa
# model NC_039477 feature CDS#1 predicted seqs saved in:           va-noro.9.vadr.NC_039477.CDS.1.fa
# model NC_039477 feature gene#2 predicted seqs saved in:          va-noro.9.vadr.NC_039477.gene.2.fa
# model NC_039477 feature CDS#2 predicted seqs saved in:           va-noro.9.vadr.NC_039477.CDS.2.fa
# model NC_039477 feature gene#3 predicted seqs saved in:          va-noro.9.vadr.NC_039477.gene.3.fa
# model NC_039477 feature CDS#3 predicted seqs saved in:           va-noro.9.vadr.NC_039477.CDS.3.fa
# model NC_039477 feature mat_peptide#1 predicted seqs saved in:   va-noro.9.vadr.NC_039477.mat_peptide.1.fa
# model NC_039477 feature mat_peptide#2 predicted seqs saved in:   va-noro.9.vadr.NC_039477.mat_peptide.2.fa
# model NC_039477 feature mat_peptide#3 predicted seqs saved in:   va-noro.9.vadr.NC_039477.mat_peptide.3.fa
# model NC_039477 feature mat_peptide#4 predicted seqs saved in:   va-noro.9.vadr.NC_039477.mat_peptide.4.fa
# model NC_039477 feature mat_peptide#5 predicted seqs saved in:   va-noro.9.vadr.NC_039477.mat_peptide.5.fa
# model NC_039477 feature mat_peptide#6 predicted seqs saved in:   va-noro.9.vadr.NC_039477.mat_peptide.6.fa
# per-sequence tabular annotation summary file saved in:           va-noro.9.vadr.sqa
# per-sequence tabular classification summary file saved in:       va-noro.9.vadr.sqc
# per-feature tabular summary file saved in:                       va-noro.9.vadr.ftr
# per-model-segment tabular summary file saved in:                 va-noro.9.vadr.sgm
# per-model tabular summary file saved in:                         va-noro.9.vadr.mdl
# per-alert tabular summary file saved in:                         va-noro.9.vadr.alt
# alert count tabular summary file saved in:                       va-noro.9.vadr.alc
# 5 column feature table output for passing sequences saved in:    va-noro.9.vadr.pass.tbl
# 5 column feature table output for failing sequences saved in:    va-noro.9.vadr.fail.tbl
# list of passing sequences saved in:                              va-noro.9.vadr.pass.list
# list of failing sequences saved in:                              va-noro.9.vadr.fail.list
# list of alerts in the feature tables saved in:                   va-noro.9.vadr.alt.list
#
# All output files created in directory ./va-noro.9/
#
# Elapsed time:  00:02:09.59
#                hh:mm:ss
# 
[ok]
```

All of these files were created in the newly created directory
`va-noro-9`.  What follows is a brief discussion of some of these
output file types.  More information on these files and their formats
can be found [here](formats.md#annotate).

The first three files are the [`.log` file](formats.md#log), which is
the same as the standard output printed to the screen currently being
discussed, the [`.cmd` file](formats.md#cmd), and the [`.filelist`
file](formats.md#filelist) which lists the output files created by
`v-annotate.pl`. Next comes a [`.seqstat file](annotate.md#seqstat)
with lengths for each sequence in the input file.

Then comes about 40 FASTA-formatted sequence files with subsequences
of the predicted features for each sequence, grouped by each
model. There will be one file for each model/feature pair that is
annotated in at least one sequence.  Each file will contain all
subsequences of input sequences that match best to the corresponding
model, and that have the corresponding feature annotated. An example
of this type of file is `va-noro.9.vadr.NC_031324.CDS.2.fa`, which
contains the predicted CDS sequences for CDS #2 for the sequences for
which the best matching model was `NC_031324`. The convention used for
naming the sequences in these files is explained
[here](formats#seqnames).

After the FASTA files are seven tabular summary files that end with
three letter suffixes:

| suffix | description | reference | 
|--------|-------------|-----------|
| `.alc` | per-alert code information (counts)     | [description of format](formats.md#alc) |
| `.alt` | per-alert instance information          | [description of format](formats.md#alt) |
| `.ftr` | per-feature information                 | [description of format](formats.md#ftr) |
| `.mdl` | per-model information                   | [description of format](formats.md#mdl) |
| `.sgm` | per-segment information                 | [description of format](formats.md#sgm) |
| `.sqa` | per-sequence annotation information     | [description of format](formats.md#sqa) |
| `.sqc` | per-sequence classification information | [description of format](formats.md#sqc) |

The contents of the `.mdl` and `.alc` files were already output by
`v-annotate.pl` as covered above. To get more information on each
sequence, see the `.sqa` and `.sqc` files. The `.sqc` file
(`va-noro.vadr.sqc`) includes
information on the classification of each sequence:

```
#seq  seq          seq                                   sub                    seq    mdl         num                             sub    score  diff/  seq   
#idx  name         len  p/f   ant  model1     grp1       grp1   score  sc/nt    cov    cov  bias  hits  str  model2     grp2       grp2    diff     nt  alerts
#---  ----------  ----  ----  ---  ---------  ---------  ----  ------  -----  -----  -----  ----  ----  ---  ---------  ---------  ----  ------  -----  ------
1     KY887602.1  7547  PASS  yes  NC_039477  Norovirus  GII   8142.8  1.079  1.000  0.997  12.5     1    +  NC_029647  Norovirus  GIV   5829.8  0.772  -     
2     KT818729.1   243  PASS  yes  NC_001959  Norovirus  GI     130.7  0.538  0.996  0.032     0     1    +  NC_029647  Norovirus  GIV     56.5  0.233  -     
3     EU437710.1   291  PASS  yes  NC_001959  Norovirus  GI     249.4  0.857  1.000  0.038     0     1    +  NC_029645  Norovirus  GIII   166.1  0.571  -     
4     DQ288307.1  1094  PASS  yes  NC_029645  Norovirus  GIII   973.4  0.890  1.000  0.150   6.2     1    +  NC_001959  Norovirus  GI     671.1  0.613  -     
5     AY237464.1   255  PASS  yes  NC_039477  Norovirus  GII    221.1  0.867  1.000  0.034     0     1    +  NC_029647  Norovirus  GIV    133.9  0.525  -     
6     KF475958.1   275  PASS  yes  NC_031324  Norovirus  GI     163.1  0.593  0.993  0.035     0     1    +  NC_029647  Norovirus  GIV     85.0  0.309  -     
7     AB713840.1   347  PASS  yes  NC_008311  Norovirus  GV     330.6  0.953  1.000  0.047   0.2     1    +  NC_029645  Norovirus  GIII   244.6  0.705  -     
8     JN585032.1   286  PASS  yes  NC_029645  Norovirus  GIII   242.3  0.847  0.997  0.039   0.1     1    +  NC_031324  Norovirus  GI     157.4  0.550  -     
9     JN975492.1  7286  FAIL  yes  NC_008311  Norovirus  GV    4666.2  0.640  1.000  0.987  16.8     1    +  NC_039475  Norovirus  GII   3478.9  0.477  -   ```

This file includes per-sequence information on whether each sequence
passed or failed, its best and second-best matching model, and scores
and coverage. The difference in score between the best and second-best
model gives an indication of how confidently classified the sequence
is to the best-matching model. (The `qstgroup` and `qstsbgrp` alerts
are reported if these scores are not sufficiently far apart to alert
the user that a sequence is not confidently classified.) Note that the
sequence `JN975492.1` is the only sequence that failed. It has no `seq
alerts` listed in the final field, so the fatal alerts that caused it
to fail must have been per-feature alerts. These can be seen in the
`.ftr` and `.alt` files. An explanation of all fields in the `.sqc`
file type is [here](formats.md#sqc).

Next, take a look at the first few lines of the `.ftr` file (`va-noro9.vadr.ftr`):

```
#     seq          seq                   ftr          ftr                         ftr  ftr                                                                                     seq         model  ftr   
#idx  name         len  p/f   model      type         name                        len  idx  str  n_from  n_to  n_instp  trc    p_from  p_to  p_instp  p_sc  nsa  nsn        coords        coords  alerts
#---  ----------  ----  ----  ---------  -----------  -------------------------  ----  ---  ---  ------  ----  -------  -----  ------  ----  -------  ----  ---  ---  ------------  ------------  ------
1.1   KY887602.1  7547  PASS  NC_039477  gene         ORF1                       5083    1    +       1  5083        -  5'          -     -        -     -    1    0     1..5083:+    22..5104:+  -     
1.2   KY887602.1  7547  PASS  NC_039477  CDS          nonstructural_polyprotein  5083    2    +       1  5083        -  5'          2  5080        -  9094    1    0     1..5083:+    22..5104:+  -     
1.3   KY887602.1  7547  PASS  NC_039477  gene         ORF2                       1623    3    +    5064  6686        -  no          -     -        -     -    1    0  5064..6686:+  5085..6707:+  -     
1.4   KY887602.1  7547  PASS  NC_039477  CDS          VP1                        1623    4    +    5064  6686        -  no       5064  6683        -  2878    1    0  5064..6686:+  5085..6707:+  -     
1.5   KY887602.1  7547  PASS  NC_039477  gene         ORF3                        807    5    +    6686  7492        -  no          -     -        -     -    1    0  6686..7492:+  6707..7513:+  -     
1.6   KY887602.1  7547  PASS  NC_039477  CDS          VP2                         807    6    +    6686  7492        -  no       6686  7486        -  1393    1    0  6686..7492:+  6707..7513:+  -     
1.7   KY887602.1  7547  PASS  NC_039477  mat_peptide  p48                         979    7    +       1   979        -  5'          -     -        -     -    1    0      1..979:+    22..1000:+  -     
1.8   KY887602.1  7547  PASS  NC_039477  mat_peptide  NTPase                     1098    8    +     980  2077        -  no          -     -        -     -    1    0   980..2077:+  1001..2098:+  -     
1.9   KY887602.1  7547  PASS  NC_039477  mat_peptide  p22                         531    9    +    2078  2608        -  no          -     -        -     -    1    0  2078..2608:+  2099..2629:+  -     
1.10  KY887602.1  7547  PASS  NC_039477  mat_peptide  VPg                         399   10    +    2609  3007        -  no          -     -        -     -    1    0  2609..3007:+  2630..3028:+  -     
1.11  KY887602.1  7547  PASS  NC_039477  mat_peptide  Pro                         543   11    +    3008  3550        -  no          -     -        -     -    1    0  3008..3550:+  3029..3571:+  -     
1.12  KY887602.1  7547  PASS  NC_039477  mat_peptide  RdRp                       1530   12    +    3551  5080        -  no          -     -        -     -    1    0  3551..5080:+  3572..5101:+  -     
#
```
This file includes information on each annotated feature, organized by
sequence. The lines above show the annotated features for the first
sequence `KY887602.1` and include information on the feature length,
strand, positions in the sequence and in the reference model, whether
it is truncated or not, and where the blastx protein validation stage
predicted coordinates are. The `seq coords` and `mdl coords` lines
near the end show the sequence and model coordinates of each feature
in the VADR coordinate string format, described
[here](formats.md#coords). 
The final column lists any alerts
pertaining to each feature. For the features in this sequence there
are none, but for the final sequence, some alerts are listed for the
second CDS. An explanation of all fields in the `.ftr`
file type is [here](formats.md#ftr).

Another important file is the `.alt` output file (`va-noro9.vadr.alt`)
which includes one line per alert reported:

```
#      seq                    ftr   ftr   ftr  alert           alert                        alert 
#idx   name        model      type  name  idx  code      fail  desc                         detail
#----  ----------  ---------  ----  ----  ---  --------  ----  ---------------------------  ------
9.1.1  JN975492.1  NC_008311  CDS   VF1     6  mutendcd  yes   MUTATION_AT_END              expected stop codon could not be identified, predicted CDS stop by homology is invalid [TCA ending at position 5685 on + strand]
9.1.2  JN975492.1  NC_008311  CDS   VF1     6  cdsstopn  yes   CDS_HAS_STOP_CODON           in-frame stop codon exists 5' of stop position predicted by homology to reference [revised to 5044..5277 (stop shifted 408 nt)]
9.1.3  JN975492.1  NC_008311  CDS   VF1     6  cdsstopp  yes   CDS_HAS_STOP_CODON           stop codon in protein-based alignment [stop codon(s) end at position(s) 5277]
9.1.4  JN975492.1  NC_008311  CDS   VF1     6  indf3pst  yes   INDEFINITE_ANNOTATION_END    protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint [36 > 5 (strand:+ CM:5685 blastx:5649, no valid stop codon in CM prediction)]
9.2.1  JN975492.1  NC_008311  CDS   VP2     8  indf5pst  yes   INDEFINITE_ANNOTATION_START  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint [54 > 5 (strand:+ CM:6656 blastx:6710)]
```     

All alerts are for the `JN975492.1` sequence. Four are for the VF1 CDS
and 1 is for the VP2 CDS. The alert codes are listed in the seventh
column, along with a brief description in the eight column followed by
a more detailed description, sometimes with positional information, at
the end of the line. All possible alerts are listed in the [alert
table](#alertable).

In addition to the tabular output files, `v-annotate.pl` also creates 
5-column tab-delimited feature table files that end with the suffix
`.tbl`. There is a separate file for passing
(`va-noro9.vadr.pass.tbl`) and failing (`va-noro9.vadr.fail.tbl`)
sequences. 
The format of the `.tbl` files is described here:
https://www.ncbi.nlm.nih.gov/Sequin/table.html.
These files contain some of the same information as the `.ftr` files
discussed above, with some additional information on *qualifiers* for
features read from the [model info file](formats.md#minfo). 
From the `va-noro9.vadr.pass.tbl` file, the feature table for
`KY887602` is: 

```
>Feature KY887602.1
<1	5083	gene
			gene	ORF1
<1	5083	CDS
			product	nonstructural polyprotein
			codon_start	2
<1	979	mat_peptide
			product	p48
980	2077	mat_peptide
			product	NTPase
2078	2608	mat_peptide
			product	p22
2609	3007	mat_peptide
			product	VPg
3008	3550	mat_peptide
			product	Pro
3551	5080	mat_peptide
			product	RdRp
5064	6686	gene
			gene	ORF2
5064	6686	CDS
			product	VP1
6686	7492	gene
			gene	ORF3
6686	7492	CDS
			product	VP2
```

In the `.fail.tbl` file, each sequence's feature table contains a
special section at the end of the table that lists errors
due to reported fatal alerts. For example, at the end of
`va-noro9.vadr.fail.tbl`:

```
Additional note(s) to submitter:
ERROR: CDS_HAS_STOP_CODON: (VF1) in-frame stop codon exists 5' of stop position predicted by homology to reference [revised to 5044..5277 (stop shifted 408 nt)]
ERROR: CDS_HAS_STOP_CODON: (VF1) stop codon in protein-based alignment [stop codon(s) end at position(s) 5277]
ERROR: INDEFINITE_ANNOTATION_END: (VF1) protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint [36 > 5 (strand:+ CM:5685 blastx:5649, no valid stop codon in CM prediction)]
ERROR: INDEFINITE_ANNOTATION_START: (VP2) protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint [54 > 5 (strand:+ CM:6656 blastx:6710)]
```

Note that this file lists only four ERRORs while the `.alt` output
file above listed five alerts. **Not all** fatal alerts will be
printed to this `.fail.tbl` file, because when specific pairs of
alerts occur, only one is output to reduce the number of overlapping
or redundant problems reported to the submitter/user. In this case the
`mutendcd` (`MUTATION_AT_END`) alert is not reported in the
`.fail.tbl` file because it occurs in combination with a `cdsstopn`
(`CDS_HAS_STOP_CODON`) alert.  More information on pairs of alerts
that cause one to be unreported in `.fail.tbl` and `.fail.list` files
is [below](#invalidated).

The final three output files all end in `.list`. The `.pass.list` and
`.fail.list` files are simply lists of all passing and failing
sequences, respectively, with each sequence listed on a separate
line. The `.alt.list` file lists all alerts that cause errors in the
`.fail.tbl` file in a four field  tab-delimited format described
[here](formats.md#altlist). The `va-noro9.vadr.alt.list` file lists
the four alerts/errors in the `.fail.tbl` file shown above:

```
#sequence	error	feature	error-description
JN975492.1	CDS_HAS_STOP_CODON	VF1	in-frame stop codon exists 5' of stop position predicted by homology to reference [revised to 5044..5277 (stop shifted 408 nt)]
JN975492.1	CDS_HAS_STOP_CODON	VF1	stop codon in protein-based alignment [stop codon(s) end at position(s) 5277]
JN975492.1	INDEFINITE_ANNOTATION_END	VF1	protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint [36 > 5 (strand:+ CM:5685 blastx:5649, no valid stop codon in CM prediction)]
JN975492.1	INDEFINITE_ANNOTATION_START	VP2	protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint [54 > 5 (strand:+ CM:6656 blastx:6710)]
```

---
## `v-annotate.pl` command-line options<a name="options"></a>

---
OTHER SECTIONS TODO :
1. alert table
2. example showing how to make alert non-fatal
3. invalidated alerts not in .tbl and .alt.list files.
4. -p



