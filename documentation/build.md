# `v-build.pl` example usage and command-line options

* [`v-build.pl` example usage](#exampleusage)
* [`v-build.pl` command-line options](#options)
  * [basic options](#options-basic)
  * [options for controlling what feature types are stored in the output model info file](#options-featuretypes) 
  * [options for controlling what qualifiers are stored in the output model info file](#options-qualifiers)
  * [options for including additional model attributes (`group` and `subgroup`)](#options-attributes) 
  * [options for controlling the CDS translation step (translation tables)](#options-tt) 
  * [options for controlling the `cmbuild` step](#options-cmbuild)
  * [options for skipping stages](#options-skip)
  * [options for additional output files](#options-output)
* [Building a VADR model library](#library)
* [How the VADR 1.0 model library was constructed](#1.0library)
  * [Norovirus models](#1.0library-noro)
  * [Dengue virus models](#1.0library-dengue)
  * [Hepatitis C virus (HCV) models](#1.0library-hcv)
  * [Other 173 *Caliciviridae* and *Flaviviridae* models](#1.0library-173)
  * [Concatenating files for 197 models to make the library](#1.0library-concat)

## `v-build.pl` example usage<a name="exampleusage"></a>

`v-build.pl` creates the VADR model files for a specified reference
GenBank (typically RefSeq) sequence necessary for running
`v-annotate.pl` to validate and annotate sequences similar to that
reference sequence. To determine the command-line usage of 
`v-build.pl` (or any VADR script), use the `-h` option, like this:

```
v-build.pl -h 
```

You'll see something like the following output:
```
# v-build.pl :: build homology model of a single sequence for feature annotation
# VADR 0.991 (Aug 2019)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Mon Oct 21 06:18:24 2019
#
Usage: v-build.pl [-options] <accession> <path to output directory to create>
```

The first few lines are the banner which show the name of the VADR
script being run along with the version and release date. This is
followed by the time and date the command was executed.  The `Usage:`
line details the expected command line arguments.  `v-build.pl` takes
as input two command line arguments, the name of the reference
accession (`<accession>`) and the name of the output directory you
want it to create (`<path to output directory to create>`) and
populate with output files, including model files for the accession to
use with `v-annotate.pl`.  

After that comes a list of all available
command-line options. These are explained in more detail [below](#build-options).

Here is an example `v-build.pl` command using the RefSeq
accession `NC_039897` a Norovirus GI complete genome sequence, and creating
an output directory with the same name as the accession:

```
`v-build.pl NC_039897 NC_039897`
```

The standard output of `v-build.pl` that is printed to the screen
(which is also output to the `.log` output file) begins with the
banner and date again followed by a list of relevant environment
variables, the command line arguments used and any command line
options used:

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

No command line options were used in our example output, but if they
were information on them would have appeared after the `output
directory` line.

Next, information is output about each step the script is proceeding
through. When each step is completed, the elapsed time for that step
is output. The first few steps for `v-build.pl` are to fetch and parse 
a FASTA sequence file and feature table file for the input accession
`NC_039897`. 

```
# Fetching FASTA file                                          ... done. [    3.9 seconds]
# Parsing FASTA file                                           ... done. [    0.0 seconds]
# Fetching feature table file                                  ... done. [    3.1 seconds]
# Parsing feature table file                                   ... done. [    0.0 seconds]
```

Next, any relevant feature tables for proteins referenced in the 
`NC_039897` feature table are fetched and parsed. The fetched
data from GenBank is then pruned to remove information that will
not be stored in the VADR model files:

```
# Fetching and parsing protein feature table file(s)           ... done. [   13.5 seconds]
# Pruning data read from GenBank                               ... done. [    0.0 seconds]
```

The fetched FASTA file is then reformatted to Stockholm for input 
to `cmbuild` to create the model CM file, feature information read
from the feature tables is finalized internally, and CDS features
are translated and used to create the BLAST database with `makeblastdb`:

```
# Reformatting FASTA file to Stockholm file                    ... done. [    0.2 seconds]
# Finalizing feature information                               ... done. [    0.0 seconds]
# Translating CDS and building BLAST DB                        ... done. [    0.2 seconds]
```

Next, the covariance model is built using `cmbuild`. This is by far
the slowest step of `v-build.pl` and a rough estimate of how long it
will take is output. The longer the sequence being modeled the longer
this step will take.  When that completes the cm file is pressed with
`cmpress` to prepare it for use with `v-annotate.pl` and the final
step is to create the model info file:

```
# Building model (should take roughly 10-30 minutes)           ... done. [  818.7 seconds]
# Pressing CM file                                             ... done. [    0.6 seconds]
# Creating model info file                                     ... done. [    0.0 seconds]
```

When all steps are complete, `v-build.pl` ends by outputting a list of
relevant output files with brief descriptions:

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

These files include the FASTA and and their formats are described more
[here](formats.md).

Only some of these files will be used by `v-annotate.pl`. These are
the files with the following suffixes:

* `vadr.cm`: the covariance model (CM) file
* `.vadr.cm.i1{i,m,f,p}`: the CM index files
* `.vadr.protein.fa.p{hr,in,sq}`: the BLAST DB files 
* `.vadr.minfo`: the VADR model info file that specifies locations and attributes of features 

You can use only this model for `NC_039897` when annotating sequences
with `v-annotate.pl`, by using the `-m`, `-b` and `-i` options as
explained in the [`v-annotate.pl` documentation](annotate.md), or you
can combine these files together with analogous files from additional
`v-build.pl` runs for other accessions to create a VADR model library,
and then use `v-annotate.pl` `-m`, `-b` and `-i` options to 
specify that library be used. This is explained in more detail 
[below](#library). The VADR 1.0 library was created in
this manner, as explained in [another section](#1.0library).

---
## `v-build.pl` command-line options<a name="options"></a>

To get a list of command-line options, execute:

`v-build.pl -h`

This will output the usage and available command-line options. 
Each option has a short description, but additional information on many 
of these options can be found below.
For `v-build.pl` the available options are split into 8 different categories, 
each explained in their own subsection below.

### `v-build.pl` basic options<a name="options-basic"></a>

The first category of options are the *basic* options:

| .......option....... | explanation | 
|--------|-------------|
| `-f`   | if `<output directory>` already exists, then using this option will cause it to be overwritten, otherwise the progam exits in error |
| `-v`   | *verbose* mode: all commands will be output to standard output as they are run | 
| `--stk <s>` | single sequence stockholm alignment will be read from `<s>`, possibly with secondary structure annotation |
| `--infa <s>` | instead of fetching the FASTA for this sequence from GenBank, use the sequence in file <s> |
| `--inft <s>` | instead of fetching the feature table for this sequence from GenBank, use the feature table in file <s> |
| `--ftfetch1` | use `efetch` program (must be in your `PATH`) to fetch feature table with `efetch -format ft` instead of default method of fetching from an `eutils` URL | 
| `--ftfetch2` | use `efetch` program (must be in your `PATH`) to fetch feature table with `efetch -format gbc | xml2tbl` instead of default method of fetching from an `eutils` URL | 
| `--gb` | fetch and parse a GenBank-format file from GenBank instead of a feature table | 
| `--ingb <s>` | read the GenBank-format file in `<s>` instead of a feature table file (requires `--gb`)| 
| `--addminfo <s>` | add arbitrary feature info in file `<s>` to output `.minfo` file, see an example [here](#1.0library-noro") | 
| `--keep` | keep additional output files that are normally removed |

### `v-build.pl` options for controlling what feature types are stored in the output model info file<a name="options-featuretypes"></a>

By default, only `CDS`, `gene` and `mat_peptide` feature types read from the GenBank feature table file
will be stored in the output `.minfo` file. This default set can be changed using the following three
command line options. For an example of using the `--fadd` option, see the construction of the dengue virus
RefSeq models for the VADR 1.0 model library [here](#1.0library-dengue").

| ......option...... | explanation | 
|--------|-------------| 
| `--fall`   | specify that all feature types (except those in `<s>` from `--fskip <s>`)  be added to the `.minfo` output file |
| `--fadd <s>`  | add feature types listed in `<s>` to the default set, where `<s>` is a comma-separated string with each feature type separated by a comma with no whitespace |
| `--fskip <s>`  | do not store information for feature types listed in `<s>`, where `<s>` is a comma-separated string with each feature type separated by a comma with no whitespace; `<s>` may contain feature types from the default set, or from other features (if `--fall` also used) |


### `v-build.pl` options for controlling what qualifiers are stored in the output model info file<a name="options-qualifiers"></a>

By default, only `product`, `gene` and `exception` qualifiers read
from the GenBank feature table file will be stored in the output
`.minfo` file, and then only for the feature types that will be
stored. This default set can be changed using the following five
command line options. 
For an example of using the `--qadd` and `--qftradd` options, see
the construction of the dengue virus RefSeq models for the VADR 1.0
model library [here](#1.0library-dengue").

| ......option...... | explanation | 
|--------|-------------| 
| `--qall`   | specify that all qualifiers (except those in `<s> from `--qskip <s>`) be added to the `.minfo` output file |
| `--qadd <s>`  | add qualifiers listed in `<s>` to the default set, where `<s>` is a comma-separated string with each qualifier separated by a comma with no whitespace |
| `--qftradd <s>`  | specify that the qualifiers listed in `<s2>` from `qadd <s2>` only apply for feature types in the string `<s>`, where `<s>` is a comma-separated string with each qualifier separated by a comma with no whitespace |
| `--qskip <s>`  | do not store information for qualifiers listed in `<s>`, where `<s>` is a comma-separated string with each qualifier separated by a comma with no whitespace; `<s>` may contain qualifiers from the default set, or from other qualifiers (if `--qall` also used) |
| `--noaddgene`  | do not automatically add `gene` qualifiers from `gene` features to any overlapping non-gene features | 

### `v-build.pl` options for including additional model attributes<a name="options-attributes"></a>

Besides qualifiers read from GenBank and information included in the input
`.minfo` file with the `--addminfo` option, two additional attributes
can be added using the command line options `--group` and
`--subgroup`. For an example of using these options see
construction of the Norovirus VADR 1.0 library model files [here](#1.0library-noro).

| .....option..... | explanation | 
|--------|-------------| 
| `--group <s>`   | specify that the model `group` attribute is `<s>`, e.g. `Norovirus` |
| `--subgroup <s>`   | specify that the model `subgroup` attribute is `<s>`, e.g. `GI`, requires `--group` |

### `v-build.pl` options for controlling the CDS translation step<a name="options-tt"></a>

By default, the NCBI translation table 1 is used to translate CDS
sequences into proteins. This can be changed to use translation table
`<s>` with the `--ttbl` option. A non-1 translation table will be
stored in the output `.minfo` file so that `v-annotate.pl` is aware of
the translation table to use when analyzing CDS predictions.

| ......option...... | explanation | 
|--------|-------------| 
| `--tt <s>`   | specify that NCBI translation table `<s>` be used instead of `1` |

Reference on NCBI translation tables:
https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes

### `v-build.pl` options for controlling the `cmbuild` step<a name="options-cmbuild"></a>

Several options exist for controlling the command-line options that will be passed
to Infernal's `cmbuild` program. For more information on these options and how 
they control `cmbuild`, see the Infernal 
User's Guide (http://eddylab.org/infernal/Userguide.pdf).

| ......option...... | explanation | 
|--------|-------------| 
| `--cmn <n>` | set the number of seqs for glocal forward profile HMM filter calibration to `<n>` (sets the `cmbuild --EgfN` option), default is to use default `cmbuild` value | 
| `--cmp7ml` | set CM's filter profile HMM as the maximum likelihood profile HMM (sets the `cmbuild --p7ml` option) |
| `--cmere` | set CM relative entropy target bits to position to `<x>` (sets the `cmbuild --ere <x>` option), default is to use default `cmbuild` value |
| `--cmemaxseq` | set CM maximum allowed effective sequence # for CM to `<x>` (sets the `cmbuild --emaxseq <x>` option) | 
| `--cminfile` | read `cmbuild` options from an input file `<s>`, the contents of the file (after removing newlines) will be supplied directly to `cmbuild` as an options string (possibly with more than one option separated by whitespace) |

### `v-build.pl` options for skipping stages<a name="options-skip"></a>

| ......option...... | explanation | 
|--------|-------------| 
| `--skipbuild` | skip the `cmbuild` step; this is mostly useful for debugging purposes, but also possibly for creating a different `.minfo` file using options like `--fadd` and/or `--qadd` for an already created model without the need to wait for the slow `cmbuild` step | 
| `--onlyurl` | output the url for the GenBank feature table file (or GenBank-format file if `--gb` also used, and exit; possibly helpful if `v-build.pl` is having trouble fetching from GenBank, you can manually download the feature table file from a browswer with the output url and use the downloaded file with `--inft` | 

### `v-build.pl` options for optional output files<a name="options-output"></a>

| .......option....... | explanation | 
|--------|-------------| 
| `--sgminfo <s>` | output information on the internal data structure used to keep track of segments of features to `<s>`, mainly useful for debugging |
| `--ftrinfo <s>` | output information on the internal data structure used to keep track of features to `<s>`, mainly useful for debugging |

---
## Building a VADR model library<a name="library"></a>

Follow these steps to build a VADR model library:

1. Run `v-build.pl` multiple (`N>1`) times for different accessions.
2. Concatenate all resulting `N` `.vadr.minfo` files into a single file, call it `my.vadr.minfo`, for example.
3. Concatenate all resulting `N` `.vadr.cm` files into a single file, call it `my.vadr.cm`, for example.
4. Move all resulting BLAST DB files (`.vadr.protein.fa.phr`, `.vadr.protein.fa.pin`, and
`.vadr.protein.fa.psq`) from all `N` runs into the same directory, call it `my-vadr-blastdb-dir`
for example.
5. Run `cmpress` on the `my.vadr.cm` file created in step 3 like this:
```
$VADRINFERNALDIR/cmpress my.vadr.cm
```

You can then use this library with `v-annotate.pl` as follows:
```
v-annotate.pl -m my.vadr.cm -i my.vadr.minfo -b my-vadr-blastdb-dir
```

Or substitute the full paths to the CM and minfo files and to the
BLAST DB directory if they are not in the current directory you are
running `v-annotate.pl` from.

If you move `my.vadr.cm` into a new directory, make sure you also move
the four `cmpress` output index files (`my.vadr.cm.i1i`, `my.vadr.cm.i1m`, 
`my.vadr.cm.i1f` and `my.vadr.cm.i1p`) into
the same directory. 

---
## How the VADR 1.0 model library was constructed <a name="1.0library"></a>

The VADR 1.0 model library consists of 197 VADR models. Nine of these
are Norovirus RefSeq models, listed in
[vadr/documentation/build-files/norovirus.9.list](build-files/norovirus.9.list).
Four of these are Dengue virus RefSeq models, listed in
[vadr/documentation/build-files/dengue.4.list](build-files/dengue.4.list).
Eight of these are Hepatitis C virus Refseq models, listed in
[vadr/documentation/build-files/hcv.8.list](build-files/hcv.8.list).
The remaining 173 are additional *Caliciviridae* and *Flaviviridae*
RefSeq models, listed in
[vadr/documentation/build-files/non-noro-dengue-hcv.173.list](build-files/non-noro-dengue-hcv.173.list).

### Building the VADR 1.0 library Norovirus models <a name="1.0library-noro"></a>

To build models for each of the nine norovirus RefSeqs listed in 
[vadr/documentation/build-files/norovirus.9.list](build-files/norovirus.9.list),
run `v-build.pl` nine separate times as follows: 

```
v-build.pl --group Norovirus --subgroup GI NC_001959 NC_001959
v-build.pl --group Norovirus --subgroup GV NC_008311 NC_008311
v-build.pl --group Norovirus --subgroup GIII NC_029645 NC_029645
v-build.pl --group Norovirus --subgroup GII --addminfo NC_029646.addminfo NC_029646 NC_029646
v-build.pl --group Norovirus --subgroup GIV NC_029647 NC_029647
v-build.pl --group Norovirus --subgroup GI NC_031324 NC_031324
v-build.pl --group Norovirus --subgroup GII --addminfo NC_039475.addminfo NC_039475 NC_039475
v-build.pl --group Norovirus --subgroup GII NC_039476 NC_039476
v-build.pl --group Norovirus --subgroup GII NC_039477 NC_039477
```

(The shell script
[vadr/documentation/build-files/norovirus.9.build.sh](build-files/norovirus.9.build.sh)
contains these commands.)

The `--group` and `--subgroup` options specify the group and subgroup
values that will be added to the output `.minfo` file. These will
enable `v-annotate.pl` to label sequences classified to these models
with the corresponding group and subgroup, as well as to fail
sequences users expect are norovirus sequences that actually are
classified best to a non-norovirus model.

The `--addminfo` options specify an input file that contains information
on additional feature attributes not from GenBank that are desired in the output
`.minfo` file. Those input files are located here:
[vadr/documentation/build-files/NC_029646.addminfo](build-files/NC_029646.addminfo)
and 
[vadr/documentation/build-files/NC_039475.addminfo](build-files/NC_039475.addminfo).
The `NC_029646.addminfo` file is in the [`.minfo` format](formats.md#minfo-format) and looks like this:

```
MODEL NC_029646 
FEATURE NC_029646 type:"CDS" coords:"5085..6692:+" xmaxins_exc:"297:36"
```

This file specifies that the additional `<key>:<value>` pair of
`xmaxins_exc:"297:36"` be added to the CDS feature with coordinates
`5085..6692:+`. (The VADR coordinate string format is described
[here](formats.md#coords-format)). The `NC_039475.addminfo` file is
similar except with the value `"295:36"`.  These two additions allow
the corresponding CDS features to have an exception to the default
maximum allowed insert length by `v-annotate.pl` without causing an
*insertnp* alert, setting it as 36 after position 297 in `NC_029646`
and after position 295 in `NC_039475`.  This change to the default was
allowed after GenBank indexers observed a common biologically valid
insertion at these positions of length 36 nucleotides (nt), which
exceeds the default maximum of 27 nt.

### Building the VADR 1.0 library Dengue virus models <a name="1.0library-dengue"></a>

The four Dengue RefSeq models are built using the `--stk` option 
to specify the secondary structure of structured regions of the 
genome at the 5' and 3' ends, as well as some additional options.
To build models for each of these RefSeqs, run
run `v-build.pl` from VADR v1.0 four separate times as follows: 

```
v-build.pl --stk NC_001477.v1.stk --qftradd stem_loop,ncRNA --qadd note,ncRNA_class --fadd stem_loop,ncRNA --group Dengue --subgroup 1 NC_001477 NC_001477
v-build.pl --stk NC_001474.v1.stk --qftradd stem_loop,ncRNA --qadd note,ncRNA_class --fadd stem_loop,ncRNA --group Dengue --subgroup 2 NC_001474 NC_001474
v-build.pl --stk NC_001475.v1.stk --qftradd stem_loop,ncRNA --qadd note,ncRNA_class --fadd stem_loop,ncRNA --group Dengue --subgroup 3 NC_001475 NC_001475
v-build.pl --stk NC_002640.v1.stk --qftradd stem_loop,ncRNA --qadd note,ncRNA_class --fadd stem_loop,ncRNA --group Dengue --subgroup 4 NC_002640 NC_002640
```

(The shell script
[vadr/documentation/build-files/dengue.4.build.sh](build-files/dengue.4.build.sh)
contains these commands.)

The `--qftradd`, `--qadd`, and `--fadd` options all take comma-separated strings
as arguments. For example, `--qftradd` takes the argument `stem_loop,ncRNA`.
These options specify that `v-build.pl` should include GenBank feature information for 
`stem_loop` and `ncRNA` features in addition to its default set of `CDS`, `gene`, and 
mat_peptide` features. Further, the `note` and `ncRNA_class` GenBank qualifiers should
be included in addition to the default set of `product`, `gene` and `exception` qualifiers
for these `stem_loop` and `ncRNA` features.

The `--group` and `--subgroup` options are used in a similar way to how they were used
to build the norovirus models.

### Building the VADR 1.0 library Hepatitis C virus models <a name="1.0library-hcv"></a>

To build models for each of the eight Hepatitis C RefSeqs listed in 
[vadr/documentation/build-files/hcv.8.list](build-files/hcv.8.list),
run `v-build.pl` eight separate times as follows: 

```
v-build.pl --group HCV --subgroup 1 NC_004102 NC_004102
v-build.pl --group HCV --subgroup 1 NC_038882 NC_038882
v-build.pl --group HCV --subgroup 2 NC_009823 NC_009823
v-build.pl --group HCV --subgroup 3 NC_009824 NC_009824
v-build.pl --group HCV --subgroup 4 NC_009825 NC_009825
v-build.pl --group HCV --subgroup 5 NC_009826 NC_009826
v-build.pl --group HCV --subgroup 6 NC_009827 NC_009827
v-build.pl --group HCV --subgroup 7 NC_030791 NC_030791
```

(The shell script
[vadr/documentation/build-files/hcv.8.build.sh](build-files/hcv.8.build.sh)
contains these commands.)

The `--group` and `--subgroup` options are used in a similar way to how they were used
to build the norovirus models.

### Building the VADR 1.0 library models for the 173 other RefSeqs <a name="1.0library-173"></a>

To build models for the other 173 *Caliciviridae* and *Flaviviridae* models listed in
[vadr/documentation/build-files/non-noro-dengue-hcv.173.list](build-files/non-noro-dengue-hcv.173.list)

Simply run `v-build.pl` from VADR v1.0 using default parameters for each accession.
For example:

```
  v-build.pl NC_034444 NC_034444 
```

(The shell script
[vadr/documentation/build-files/non-noro-dengue.hcv.173.build.sh](build-files/non-noro-dengue.hcv.173.build.sh)
will execute these 173 commands.)

Each of these commands takes roughtly between 10 minutes and an hour. 

Note: if the RefSeq annotation for any of these 197 VADR models has
changed since October 2019, then it may not be able to identically
reproduce the VADR 1.0 model library using the steps outlined
above. This is because `v-build.pl` fetches the current RefSeq
annotation data from GenBank when it is run. If necessary, contact
Eric Nawrocki at nawrocke@ncbi.nlm.nih.gov for additional files needed
to reproduce the library exactly.

### Concatenating files for the 197 models to create the VADR library <a name="1.0library-concat"></a>

After completing the steps above to make the 197 models, you can make
the VADR 1.0 library by following the instructions for creating a VADR
library [here](#library).
