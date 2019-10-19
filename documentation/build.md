# `v-build.pl` example usage and command-line options

`v-build.pl` takes as input two command line arguments, the name
of a GenBank accession (typically a RefSeq accession) 
and the name of the output directory you
want it to create and populate with output files.

Here is an example `v-build.pl` command using the RefSeq
accession `NC_039897` a Norovirus GI complete genome sequence.

```
`v-build.pl NC_039897 NC_039897`
```

The standard output of `v-build.pl` that is printed to the screen
(which is also output to the `.log` output file) begins with 
the name of the VADR
script being run and the VADR version and release date:

```
# v-build.pl :: build homology model of a single sequence for feature annotation
# VADR 0.991 (Aug 2019)
```

Next comes information on the time and date of execution, any relevant
environment variables, the command line arguments used and any command
line options used:

```
# date:              Fri Oct  4 13:11:22 2019
# $VADRBLASTDIR:     /usr/bin
# $VADREASELDIR:     /home/nawrocki/vadr-install/infernal-dev/easel/miniapps
# $VADRINFERNALDIR:  /home/nawrocki/vadr-install/infernal-dev/src
# $VADRSCRIPTSDIR:   /home/nawrocki/vadr-install/vadr
#
# accession/model name:           NC_039897
# output directory:               NC_039897
```

No command line options were used in our example output, but 

Next is the stage list section <a name="log-stagelist"></a> which lists
each stage the script proceeds through, along with the time that
elapsed during that stage: 

```
# Fetching FASTA file                                          ... done. [    3.9 seconds]
# Parsing FASTA file                                           ... done. [    0.0 seconds]
# Fetching feature table file                                  ... done. [    3.1 seconds]
# Parsing feature table file                                   ... done. [    0.0 seconds]
# Fetching and parsing protein feature table file(s)           ... done. [   13.5 seconds]
# Pruning data read from GenBank                               ... done. [    0.0 seconds]
# Reformatting FASTA file to Stockholm file                    ... done. [    0.2 seconds]
# Finalizing feature information                               ... done. [    0.0 seconds]
# Translating CDS and building BLAST DB                        ... done. [    0.2 seconds]
# Building model (should take roughly 10-30 minutes)           ... done. [  818.7 seconds]
# Pressing CM file                                             ... done. [    0.6 seconds]
# Creating model info file                                     ... done. [    0.0 seconds]
```

The summary section <a name="log-summary"></a> is optional.
`v-annotate.pl` log files will have a summary section, but
`v-build.pl` output files do not. In this section, `v-annotate.pl` log
files include information on the number of sequences classified to
each model (this is identical to the information output to the
`.mdl`(#mdl-format) output file) and the number of each type of
reported alert (this is identical to the information output to the
`.alc`(#alc-format) output file). An example of this section from the
command `v-annotate.pl $VADRSCRIPTSDIR/testfiles/noro.9.fa va-noro9`
is:

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
#
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

The summary section is followed by the output file list <a
name="log-outputfilelist"></a> which lists many (but not necessarily
all of) the files that were created by the script. This list is meant
to contain the files that are most relevant to the typical
user. Switching back to the above `v-build.pl` command example output: 

```
# Output printed to screen saved in:                                NC_039897.vadr.log
# List of executed commands saved in:                               NC_039897.vadr.cmd
# List and description of all output files saved in:                NC_039897.vadr.filelist
# fasta file for NC_039897 saved in:                                NC_039897.vadr.fa
# feature table format file for NC_039897 saved in:                 NC_039897.vadr.tbl
# feature table format file for YP_009538340.1 saved in:            NC_039897.vadr.YP_009538340.1.tbl
# feature table format file for YP_009538341.1 saved in:            NC_039897.vadr.YP_009538341.1.tbl
# feature table format file for YP_009538342.1 saved in:            NC_039897.vadr.YP_009538342.1.tbl
# Stockholm alignment file for NC_039897 saved in:                  NC_039897.vadr.stk
# fasta sequence file for CDS from NC_039897 saved in:              NC_039897.vadr.cds.fa
# fasta sequence file for translated CDS from NC_039897 saved in:   NC_039897.vadr.protein.fa
# BLAST db .phr file for NC_039897 saved in:                        NC_039897.vadr.protein.fa.phr
# BLAST db .pin file for NC_039897 saved in:                        NC_039897.vadr.protein.fa.pin
# BLAST db .psq file for NC_039897 saved in:                        NC_039897.vadr.protein.fa.psq
# CM file saved in:                                                 NC_039897.vadr.cm
# cmbuild output file saved in:                                     NC_039897.vadr.cmbuild
# binary CM and p7 HMM filter file saved in:                        NC_039897.vadr.cm.i1m
# SSI index for binary CM file saved in:                            NC_039897.vadr.cm.i1i
# optimized p7 HMM filters (MSV part) saved in:                     NC_039897.vadr.cm.i1f
# optimized p7 HMM filters (remainder) saved in:                    NC_039897.vadr.cm.i1p
# cmpress output file saved in:                                     NC_039897.vadr.cmpress
# VADR 'model info' format file for NC_039897 saved in:             NC_039897.vadr.minfo
#
# All output files created in directory ./NC_039897/
```

The final section <a name="log-timing"></a> lists how much time elapsed
while the script was executing. The final line of the output is either
`[ok]` if the script finished successfully without any unexpected
runtime errors. This final line will be `[fail]` if the script did not finish successfully
due to a runtime error. In the latter case, an error message will
occur just prior to the `[fail]` line. It may also be helpful to look
at the `.cmd`<a name="cmd-format"></a>' output file to see what the final
command was prior to failure. For the `v-build.pl` command this
section is:

```
# Elapsed time:  00:14:00.24
#                hh:mm:ss
# 
[ok]
```



- Example runs
- Explanation of options
- Common option usages 
  * different VADR library

