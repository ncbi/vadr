# VADR output file formats

VADR creates many different types of output files. You can find an explanation
of these formats below divided into three categories:

| category | description |
|--------|-----------------------|
| [generic VADR output files](#generic-formats) | files created by all VADR scripts (`v-annotate.pl` and `v-build.pl`) |
| [`v-build.pl` output files](#build-formats) | files created only by `v-build.pl` |
| [`v-annotate.pl` output files](#annotate-formats) | files created only by `v-annotate.pl` |

---
## Format of generic VADR output files created by all VADR scripts<a name="generic-formats"></a>

All VADR scripts (e.g. `v-build.pl` and `v-annotate.pl`) create a
common set of three output files. These files are named
`<outdir>.vadr.<suffix>` where `<suffix>` is either `log`, `cmd` or
`filelist` and `<outdir>` is the command line argument
that specifies the name of the output directory to create.
These files are the three types of generic files that are supported by
the [Sequip]:https://github.com/nawrockie/sequip 
ofile Perl module.

| suffix | description |
|--------|-----------------------|
| [`.log`](#logformat) | log of steps taken by the VADR script (this is identical to the standard output of the program) |
| [`.cmd`](#cmdformat) | list of the commands run using Perl's `system` command internally by the VADR script |
| [`.filelist`](#filelistformat) | list of output files created by the VADR script | 

Each format is explained in more detail below.

---
### Explanation of `.log`-suffixed output files<a name="logformat"></a>

The `.log` files include the same text that is printed to standard output. 
Log files have five or six sections: [banner](#log-banner), [input
information](#log-inputinformation), [stage
list](#log-stagelist), [summary](#log-summary) (optional), [output
file list](#log-outputfilelist), and [timing](#log-timing). 

The banner section <a name="log-banner"></a> lists the name of the VADR
script being run and the version and date. An example of this section
from the command `v-build.pl -f --group Norovirus --subgroup GI NC_039897 NC_039897` with VADR 0.991 is:

```
# v-build.pl :: build homology model of a single sequence for feature annotation
# VADR 0.991 (Aug 2019)
```

The input information section <a name="log-inputinformation"></a> includes
information on the time and date of execution, any relevant
environment variables, the command line arguments used and any command
line options used: 

```
# date:              Fri Oct  4 13:11:22 2019
# $VADRBLASTDIR:     /usr/bin
# $VADREASELDIR:     /panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal-dev/easel/miniapps
# $VADRINFERNALDIR:  /panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal-dev/src
# $VADRSCRIPTSDIR:   /panfs/pan1/dnaorg/virseqannot/code/vadr-install/vadr
#
# accession/model name:           NC_039897
# output directory:               NC_039897
# forcing directory overwrite:    yes [-f]
# specify model group is <s>:     Norovirus [--group]
# specify model subgroup is <s>:  GI [--subgroup]
```

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
`.mdl`(#mdlformat) output file) and the number of each type of
reported alert (this is identical to the information output to the
`.alc`(#alcformat) output file). An example of this section from the
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
at the `.cmd`<a name="cmdformat"></a>' output file to see what the final
command was prior to failure. For the `v-build.pl` command this
section is:

```
# Elapsed time:  00:14:00.24
#                hh:mm:ss
# 
[ok]
```

---
### Explanation of `.cmd`-suffixed output files<a name="cmdformat"></a>

The `.cmd` files simply list all the commands run by the Perl `system`
function internally by the VADR script, each separated by a newline. 
The final three lines are special. The third-to-last line lists the
date and time just before the script completed execution. The
second-to-last line lists system info (output by the `uname -a` unix
command). The final line is either `[ok]` or `[fail]`, depending on if
the script ended successfully (zero return status) or did not
(non-zero return status), respectively. If this line is
`[fail]` it will be followed by an error message. An example `.cmd`
output file for the example command `v-build.pl -f --group Norovirus
--subgroup GI NC_039897 NC_039897` is:

```
mkdir NC_039897
/panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal-dev/easel/miniapps/esl-reformat --informat afa stockholm NC_039897/NC_039897.vadr.fa > NC_039897/NC_039897.vadr.stk
/panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal-dev/easel/miniapps/esl-translate  -M -l 3 --watson NC_039897/NC_039897.vadr.cds.fa > NC_039897/NC_039897.vadr.cds.esl-translate.1.fa
rm NC_039897/NC_039897.vadr.cds.esl-translate.1.fa
rm NC_039897/NC_039897.vadr.cds.esl-translate.2.fa
rm NC_039897/NC_039897.vadr.cds.esl-translate.2.fa.ssi
/usr/bin/makeblastdb -in NC_039897/NC_039897.vadr.protein.fa -dbtype prot > /dev/null
/panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal-dev/src/cmbuild -n NC_039897 --verbose  --noss NC_039897/NC_039897.vadr.cm NC_039897/NC_039897.vadr.stk > NC_039897/NC_039897.vadr.cmbuild
/panfs/pan1/dnaorg/virseqannot/code/vadr-install/infernal-dev/src/cmpress NC_039897/NC_039897.vadr.cm > NC_039897/NC_039897.vadr.cmpress
# Fri Oct  4 13:25:23 EDT 2019
# Darwin ericsmac 18.7.0 Darwin Kernel Version 18.7.0: Tue Aug 20 16:57:14 PDT 2019; root:xnu-3273.123.2~2/RELEASE_X86_64 x86_64
[ok]
```
---
### Explanation of `.filelist`-suffixed output files<a name="filelistformat"></a>

The `.filelist` files list the output files created by the VADR
script. This list will typically include at least those files printed
in the file output file list section[log-outputfilelist] of the `.log`
file, and sometimes more. Each line includes a brief description of
each file. An example `.filelist` output file for the example command
`v-build.pl -f --group Norovirus --subgroup GI NC_039897 NC_039897`
is:

```
# fasta file for NC_039897 saved in:                                               NC_039897.vadr.fa
# feature table format file for NC_039897 saved in:                                NC_039897.vadr.tbl
# feature table format file for YP_009538340.1 saved in:                           NC_039897.vadr.YP_009538340.1.tbl
# feature table format file for YP_009538341.1 saved in:                           NC_039897.vadr.YP_009538341.1.tbl
# feature table format file for YP_009538342.1 saved in:                           NC_039897.vadr.YP_009538342.1.tbl
# Stockholm alignment file for NC_039897 saved in:                                 NC_039897.vadr.stk
# fasta sequence file for CDS from NC_039897 saved in:                             NC_039897.vadr.cds.fa
# fasta sequence file for translated CDS from NC_039897 saved in:                  NC_039897.vadr.protein.fa
# BLAST db .phr file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.phr
# BLAST db .pin file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.pin
# BLAST db .psq file for NC_039897 saved in:                                       NC_039897.vadr.protein.fa.psq
# CM file saved in:                                                                NC_039897.vadr.cm
# cmbuild output file saved in:                                                    NC_039897.vadr.cmbuild
# binary CM and p7 HMM filter file saved in:                                       NC_039897.vadr.cm.i1m
# SSI index for binary CM file saved in:                                           NC_039897.vadr.cm.i1i
# optimized p7 HMM filters (MSV part) saved in:                                    NC_039897.vadr.cm.i1f
# optimized p7 HMM filters (remainder) saved in:                                   NC_039897.vadr.cm.i1p
# cmpress output file saved in:                                                    NC_039897.vadr.cmpress
# VADR 'model info' format file for NC_039897 saved in:                            NC_039897.vadr.minfo
```

## Format of `v-build.pl` output files<a name="build-formats"></a>

`v-build.pl` creates many output files. The following table lists many
of the output files with a brief description and in some cases further
references on the file type/format. The `.minfo` file format is documented
further below. 

| file suffix | description | reference |
|--------|-----------------------|
| `.tbl`  | 5 column tab-delimited `feature table` | https://www.ncbi.nlm.nih.gov/Sequin/table.html | 


---
## Format of `v-annotate.pl` output files<a name="annotate-formats"></a>

`v-annotate.pl` creates many output files. The formats of many of these file
types are discussed below.

---
## `v-annotate.pl` tabular output files 

There are seven types of `v-annotate.pl` tabular output files with
fields separated by one or more spaces, that are meant to be easily
parseable. These files are named `<outdir>.vadr.<suffix>` where
`<outdir>` is the second command line argument given to
`v-annotate.pl`. The seven suffixes are:

| suffix | description |
|--------|-----------------------|
| [`.alc`](#alcformat) | per-alert code information (counts) |
| [`.alt`](#altformat) | per-alert instance information |
| [`.ftr`](#ftrformat) | per-feature information |
| [`.mdl`](#mdlformat) | per-model information |
| [`.sgm`](#sgmformat) | per-segment information |
| [`.sqa`](#sqaformat) | per-sequence annotation information |
| [`.sqc`](#sqcformat) | per-sequence classification information |

All seven types of tabular output files share the following
characteristics: 

1. fields are separated by whitespace (with the possible exception of
   the final field) 
2. comment lines begin with `#`
3. data lines begin with a non-whitespace character other than `#`
4. all lines are either comment lines or data lines

Each format is explained in more detail below.

---
### Explanation of `.alc`-suffixed output files<a name="alcformat"></a>

`.alc` data lines have 8 or more fields, the names of which appear in the first two
comment lines in each file. There is one data line for each alert code
that occurs at least once in the input sequence file that
`v-annotate.pl` processed.


| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of alert code |
|   2 | `alert code`          | 8 character VADR alert code |
|   3 | `causes failure`      | `yes` if this code is fatal and causes the associated input sequence to FAIL, `no` if this code is non-fatal |
|   4 | `short description`   | short description of the alert that often maps to error message from NCBI's submission system, multiple alert codes can have the same short description |
|   5 | `per type`            | `feature` if this alert pertains to a specific feature in a sequence, `sequence` if it does not |
|   6 | `num cases`           | number of instances of this alert in the output (number of rows for this alert in `.alt` file), can be more than 1 per sequence |
|   7 | `num seqs`            | number of input sequences with at least one instance of this alert |
|   8 to end | `long description`    |longer description of the alert, specific to each alert type; **this field contains whitespace** |

---
### Explanation of `.alt`-suffixed output files<a name="altformat"></a>

`.alt` data lines have 10 or more fields, the names of which appear in the first two
comment lines in each file. There is one data line for each **alert instance**
that occurs for each input sequence file that `v-annotate.pl` processed.

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of alert instance in format `<d1>.<d2>.<d3>`, where `<d1>` is the index of the sequence this alert instance pertains to in the input sequence file, `<d2>` is the index of the feature this alert instance pertains to (range 1..`<n>`, where `<n>` is the number of features in this sequence with at least 1 alert instance) and `<d3>` is the index of the alert instance for this sequence/feature pair |
|   2 | `seq name`            | sequence name | 
|   3 | `model`               | name of the best-matching model for this sequence |
|   4 | `ftr type`            | type of the feature this alert instance pertains to (e.g. CDS) |
|   5 | `ftr name`            | name of the feature this alert instance pertains to |
|   6 | `ftr idx`             | index (in input model info file) this alert instance pertains to |
|   7 | `alert code`          | 8 character VADR alert code |
|   8 | `fail`                | `yes` if this alert code is fatal (automatically causes the sequence to fail), `no` if not |
|   9 | `alert desc`          | short description of the alert code that often maps to error message from NCBI's submission system, multiple alert codes can have the same short description |
| 10 to end | `alert detail`  | detailed description of the alert instance, possibly with sequence position information; **this field contains whitespace** |

---
### Explanation of `.ftr`-suffixed output files<a name="ftrformat"></a>

`.ftr` data lines have 23 fields, the names of which appear in the first two
comment lines in each file. There is one data line for each
**feature** that is annotated for each input sequence file that
`v-annotate.pl` processed. The set of possible features for each
input sequence depend on its best-matching model, and can be found in
the model info file.

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of feature in format `<d1>.<d2>`, where `<d1>` is the index of the sequence in which this feature is annotated in the input sequence file, `<d2>` is the index of the feature (range 1..`<n>`, where `<n>` is the number of features annotated for this sequence |
|   2 | `seq name`            | sequence name in which this feature is annotated |
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alert instances) |
|   5 | `model`               | name of the best-matching model for this sequence |
|   6 | `ftr type`            | type of the feature (e.g. CDS) |
|   7 | `ftr name`            | name of the feature |
|   8 | `ftr len`             | length of the annotated feature in nucleotides in input sequence |
|   9 | `ftr idx`             | index (in input model info file) of this feature |
|  10 | `str`                 | strand on which the feature is annotated: `+` for positive/forward/Watson strand, `-` for negative/reverse/Crick strand |
|  11 | `n_from`              | nucleotide start position for this feature in input sequence |
|  12 | `n_to`                | nucleotide end position for this feature in input sequence |
|  13 | `n_instp`             | nucleotide position of stop codon not at `n_to`, or `-` if none, will be 5' of `n_to` if early stop (`cdsstopn` alert), or 3' of `n_to` if first stop is 3' of `n_to` (`mutendex` alert), or `?` if no in-frame stop exists 3' of `n_from`; will always be `-` if `trunc` is not `no`; |
|  14 | `trc`                 | indicates whether the feature is truncated or not, where one or both ends of the feature are missing due to a premature end to the sequence; possible values are `no` for not truncated; `5'` for truncated on the 5' end; `3'` for truncated on the 3' end; and `5'&3'` for truncated on both the 5' and 3' ends; |
|  15 | `p_from`              | nucleotide start position for this feature based on the blastx protein-validation step | 
|  16 | `p_to`                | nucleotide stop position for this feature based on the blastx protein-validation step | 
|  17 | `p_instp`             | nucleotide position of stop codon 5' of `p_to` if an in-frame stop exists before `p_to` |
|  18 | `p_sc`                | raw score of best blastx alignment |
|  19 | `nsa`                 | number of segments annotated for this feature |
|  20 | `nsn`                 | number of segments not annotated for this feature |
|  21 | `seq coords`          | sequence coordinates of feature, see [`format of coordinate strings`(#coordformat)] |
|  22 | `mdl coords`          | model coordinates of feature, see [`format of coordinate strings`(#coordformat)] |
|  23 | `ftr alerts`          | alerts that pertain to this feature, listed in format `SHORT_DESCRIPTION(alertcode)`, separated by commas if more than one, `-` if none |

---
### Explanation of `.mdl`-suffixed output files<a name="mdlformat"></a>

`.mdl` data lines have 7 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**model** that is the best-matching model for at least one sequence in
the input file, plus 2 additional lines, a line with `*all*` in the
`model` field reports the summed counts over all models, and a line
with `*none*` in the `model` field reports the summed counts for all
sequences that did not match any models. This information is also
included in the `.log` output file.

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of model | 
|   2 | `model`               | name of model | 
|   3 | `group`               | group of model, defined in model info file, or `-` if none |
|   4 | `subgroup`            | subgroup of model, defined in model info file, or `-`' if none | 
|   5 | `num seqs`            | number of sequences for which this model was the best-matching model |
|   6 | `num pass`            | number of sequences from `num seqs` that passed with 0 fatal alerts | 
|   7 | `num fail`            | number of sequences from `num seqs` that failed with >= 1 fatal alerts | 

### Explanation of `.sgm`-suffixed output files<a name="sgmformat"></a>

`.sgm` data lines have 21 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**segment** of a feature that is annotated for each input sequence
file that `v-annotate.pl` processed. Each feature is composed of
1 or more segments, as defined by the `coords` field in the model info
file. 

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `idx`                 | index of segment in format `<d1>.<d2>.<d3> where `<d1>` is the index of the sequence in which this segment is annotated in the input sequence file, `<d2>` is the index of the feature (range 1..`<n1>`, where `<n1>` is the number of features annotated for this sequence) and `<d3>` is the index of the segment annotated within that feature (range 1..`<n2>` where `<n2>` is the number of segments annotated for this feature | 
|   2 | `seq name`            | sequence name in which this feature is annotated |
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alert instances) |
|   5 | `model`               | name of the best-matching model for this sequence |
|   6 | `ftr type`            | type of the feature (e.g. CDS) |
|   7 | `ftr name`            | name of the feature |
|   8 | `ftr idx`             | index (in input model info file) of this feature |
|   9 | `num sgm`             | number of segments annotated for this sequence/feature pair |
|  10 | `sgm idx`             | index (in feature) of this segment |
|  11 | `seq from`            | nucleotide start position for this segment in input sequence, will be <= `seq to` if strand (`str`) `-` |
|  12 | `seq to`              | nucleotide end position for this segment in input sequence, will be >= `seq from` if strand (`str`) `-` |
|  13 | `mdl from`            | model start position for this segment, will be <= `mdl to` if strand (`str`) `-` | |
|  14 | `mdl to`              | model end position for this segment, will be >= `mdl from` if strand (`str`) `-` | |
|  15 | `sgm len`             | length, in nucleotides, for this annotated segment in the input sequence |
|  16 | `str`                 | strand (`+` or `-`) for this segment in the input sequence |
|  17 | `trc`                 | indicates whether the segment is truncated or not, where one or both ends of the segment are missing due to a premature end to the sequence; possible values are `no` for not truncated; `5'` for truncated on the 5' end; `3'` for truncated on the 3' end; and `5'&3'` for truncated on both the 5' and 3' ends; |
|  18 | `5' pp`               | posterior probability of the aligned nucleotide at the 5' boundary of the segment, or `-` if 5' boundary aligns to a gap (possibly due to a 5' truncation) }
|  19 | `3' pp`               | posterior probability of the aligned nucleotide at the 3' boundary of the segment, or `-` if 3' boundary aligns to a gap (possibly due to a 3' truncation) }
|  20 | `5' gap`              | `yes` if the 5' boundary of the segment is a gap (possibly due to a 5' truncation), else `no` }
|  21 | `3' gap`              | `yes` if the 3' boundary of the segment is a gap (possibly due to a 3' truncation), else `no` }

---
### Explanation of `.sqa`-suffixed output files<a name="sqaformat"></a>

`.sqa` data lines have 14 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**sequence** in the input sequence file file that `v-annotate.pl`
processed. `.sqa` files include **annotation** information for each
sequence. `.sqc` files include **classification** information for each
sequence. 


| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `seq idx`             | index of sequence in the input file |
|   2 | `seq name`            | sequence name | 
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alert instances) |
|   5 | `ant`                 | `yes` if this sequence was annotated, `no` if not, due to a per-sequence alert that prevents annotation |
|   6 | `best model`          | name of the best-matching model for this sequence |
|   7 | `grp`                 | group of model `best model`, defined in model info file, or `-` if none |
|   8 | `subgp`               | subgroup of model `best model`, defined in model info file, or `-`' if none | 
|   9 | `nfa`                 | number of features annotated for this sequence |
|  10 | `nfn`                 | number of features in model `best model` that are not annotated for this sequence |
|  11 | `nf5`                 | number of annotated features that are 5' truncated |
|  12 | `nf3`                 | number of annotated features that are 3' truncated |
|  13 | `nfalt`               | number of per-feature alerts reported for this sequence (does not count per-sequence alerts) |
|  14 | `seq alerts`          | per-sequence alerts that pertain to this sequence, listed in format `SHORT_DESCRIPTION(alertcode)`, separated by commas if more than one, `-` if none |

---
### Explanation of `.sqc`-suffixed output files<a name="sqcformat"></a>

`.sqc` data lines have 21 fields, the names of which appear in the
first two comment lines in each file. There is one data line for each
**sequence** in the input sequence file file that `v-annotate.pl`
processed. `.sqc` files include **classification** information for
each sequence.  `.sqa` files include **annotation** information for
each sequence. For more information on bit scores and `bias` see the Infernal User's Guide
(http://eddylab.org/infernal/Userguide.pdf)

| idx | field                 | description |
|-----|-----------------------|-------------|
|   1 | `seq idx`             | index of sequence in the input file |
|   2 | `seq name`            | sequence name | 
|   3 | `seq len`             | length of the sequence with name `seq name` | 
|   4 | `p/f`                 | `PASS` if this sequence passes, `FAIL` if it fails (has >= 1 fatal alert instances) |
|   5 | `ant`                 | `yes` if this sequence was annotated, `no` if not, due to a per-sequence alert that prevents annotation |
|   6 | `model1`              | name of the best-matching model for this sequence, this is the model with the top-scoring hit for this sequence in the classification stage |
|   7 | `grp1`                | group of model `model1`, defined in model info file, or `-` if none |
|   8 | `subgrp1`             | subgroup of model `model1`, defined in model info file, or `-` if none |
|   9 | `score`               | summed bit score for all hits on strand `str` to model `model1` for this sequence in the classification stage |
|  10 | `sc/nt`               | bit score per nucleotide; `score` divided by total length (in sequence positions) of all hits to model `model1` on strand `str` in the classification stage |
|  11 | `seq cov`             | fraction of sequence positions (`seq len`) covered by any hit to model `model1` on strand `str` in the coverage determination stage | 
|  12 | `mdl cov`             | fraction of model positions (model length - the number of reference positions in `model1`) covered by any hit to model `model1` on strand `str` in the coverage determination stage | 
|  13 | `bias`                | summed bit score due to biased composition (deviation from expected nucleotide frequencies) of all hits on strand `str` to model `model1` for this sequence in the coverage determination stage |
|  14 | `num hits`            | number of hits on strand `str` to model `model1` for this sequence in the coverage determination stage |
|  15 | `str`                 | strand with the top-scoring hit to `model1` for this sequence in the classification stage |
|  16 | `model2`              | name of the second best-matching model for this sequence, this is the model with the top-scoring hit for this sequence across all hits that are not to `model1` in the classification stage |
|  17 | `grp2`                | group of model `model2`, defined in model info file, or `-` if none |
|  18 | `subgrp2`             | subgroup of model `model2`, defined in model info file, or `-`' if none | 
|  19 | `score diff`          | bit score difference between summed bit score for all hits to `model1` on strand `str` and summed bit score for all hits to `model2` on strand with top-scoring hit to `model2` in the classification stage |
|  20 | `diff/nt`             | bit score difference per nucleotide; `sc/nt` minus sc2/nt where sc2/nt is summed bit score for all hits to `model2` on strand with top-scoring hit to `model2` in the classification stage |
|  21 | `seq alerts`          | per-sequence alerts that pertain to this sequence, listed in format `SHORT_DESCRIPTION(alertcode)`, separated by commas if more than one, `-` if none |

---
TODO:
* `v-build.pl` output formats (including `modelinfo`)
* `coords` field in modelinfo explanation
* `posterior probability` explanation
