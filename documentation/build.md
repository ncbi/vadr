#  <a name="top"></a> `v-build.pl` example usage and command-line options

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
  * [additional expert options](#options-expert)
* [Building a VADR model library](#library)
* [Practical tutorial: building a VADR RSV model library](#rsv)
* [How the VADR 1.0 model library was constructed](#1.0library)
  * [Norovirus models](#1.0library-noro)
  * [Dengue virus models](#1.0library-dengue)
  * [Hepatitis C virus (HCV) models](#1.0library-hcv)
  * [Other 173 *Caliciviridae* and *Flaviviridae* models](#1.0library-173)
  * [Concatenating files for 197 models to make the library](#1.0library-concat)

---

## `v-build.pl` example usage<a name="exampleusage"></a>

`v-build.pl` creates the VADR model files for a specified reference
GenBank (typically RefSeq) sequence necessary for running
`v-annotate.pl` to validate and annotate sequences similar to that
reference sequence. It is recommended to run `v-build.pl` only on sequences 
of length 25Kb (25,000 nucleotides) or less due to the prohibitively
large memory requirements of `v-annotate.pl` for larger models. To
determine the command-line usage of `v-build.pl` (or any VADR script),
use the `-h` option, like this:

```
v-build.pl -h 
```

You'll see something like the following output:
```
# v-build.pl :: build homology model of a single sequence for feature annotation
# VADR 1.4 (Dec 2021)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Tue Dec 21 10:59:56 2021
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
command-line options. These are explained in more detail [below](#options).

Here is an example `v-build.pl` command using the RefSeq
accession `NC_039897` a Norovirus GI complete genome sequence, and creating
an output directory with the same name as the accession:

```
v-build.pl NC_039897 NC_039897
```

***If, when you the above command, it fails and you see a (potentially
   very long) error that starts with `Can't locate
   LWP/Protocol/https.pm in @INC` or something similar. You will need
   to install one or more perl modules before `v-build.pl` will run as
   explained more [here](install.md#inline).***

The standard output of `v-build.pl` that is printed to the screen
(which is also output to the `.log` output file) begins with the
banner and date again followed by a list of relevant environment
variables, the command line arguments used and any command line
options used:

```
# $VADRBLASTDIR:     /home/nawrocki/vadr-install-dir/ncbi-blast
# $VADREASELDIR:     /home/nawrocki/vadr-install-dir/infernal/binaries
# $VADRINFERNALDIR:  /home/nawrocki/vadr-install-dir/infernal/binaries
# $VADRSCRIPTSDIR:   /home/nawrocki/vadr-install-dir/vadr
#
# accession/model name:  NC_039897
# output directory:      NC_039897
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
# Fetching FASTA file                                          ... done. [    3.3 seconds]
# Parsing FASTA file                                           ... done. [    0.0 seconds]
# Fetching feature table file                                  ... done. [    3.2 seconds]
# Parsing feature table file                                   ... done. [    0.0 seconds]
```

Next, any relevant feature tables for proteins referenced in the 
`NC_039897` feature table are fetched and parsed. The fetched
data from GenBank is then pruned to remove information that will
not be stored in the VADR model files:

```
# Fetching and parsing protein feature table file(s)           ... done. [    9.4 seconds]
# Pruning data read from GenBank                               ... done. [    0.0 seconds]
```

The fetched FASTA file is then reformatted to Stockholm for input to
`cmbuild` to create the model CM file, the BLAST nucleotide database
is construced, feature information read from the feature tables is
finalized internally, and CDS features are translated and used to
create the BLAST database with `makeblastdb`, and the HMMER protein
HMM database with `hmmbuild`.

```
# Reformatting FASTA file to Stockholm file                    ... done. [    0.3 seconds]
# Building BLAST nucleotide database                           ... done. [    1.2 seconds]
# Finalizing feature information                               ... done. [    0.0 seconds]
# Translating CDS                                              ... done. [    0.1 seconds]
# Building BLAST protein database                              ... done. [    0.2 seconds]
# Building HMMER protein database                              ... done. [    1.2 seconds]
```

Next, the covariance model is built using `cmbuild`. This is by far
the slowest step of `v-build.pl` and a rough estimate of how long it
will take is output. The longer the sequence being modeled the longer
this step will take.  When that completes the cm file is pressed with
`cmpress` to prepare it for use with `v-annotate.pl`, and a BLAST
nucleotide database is created from the CM consensus sequence. The final
step is to create the model info file:

```
# Building CM (should take roughly 10-30 minutes)              ... done. [  480.6 seconds]
# Pressing CM file                                             ... done. [    0.3 seconds]
# Building BLAST nucleotide database of CM consensus           ... done. [    0.3 seconds]
# Creating model info file                                     ... done. [    0.0 seconds]
```
When all steps are complete, `v-build.pl` ends by outputting a list of
relevant output files with brief descriptions:

```
# Output printed to screen saved in:                                 NC_039897.vadr.log
# List of executed commands saved in:                                NC_039897.vadr.cmd
# List and description of all output files saved in:                 NC_039897.vadr.filelist
# fasta file for NC_039897 saved in:                                 NC_039897.vadr.fa
# feature table format file for NC_039897 saved in:                  NC_039897.vadr.tbl
# feature table format file for YP_009538340.1 saved in:             NC_039897.vadr.YP_009538340.1.tbl
# feature table format file for YP_009538341.1 saved in:             NC_039897.vadr.YP_009538341.1.tbl
# feature table format file for YP_009538342.1 saved in:             NC_039897.vadr.YP_009538342.1.tbl
# Stockholm alignment file for NC_039897 saved in:                   NC_039897.vadr.stk
# nucleotide blastn db fasta sequence file for NC_039897 saved in:   NC_039897.vadr.fa
# BLAST db .nhr file for NC_039897 saved in:                         NC_039897.vadr.fa.nhr
# BLAST db .nin file for NC_039897 saved in:                         NC_039897.vadr.fa.nin
# BLAST db .nsq file for NC_039897 saved in:                         NC_039897.vadr.fa.nsq
# BLAST db .ndb file for NC_039897 saved in:                         NC_039897.vadr.fa.ndb
# BLAST db .not file for NC_039897 saved in:                         NC_039897.vadr.fa.not
# BLAST db .ntf file for NC_039897 saved in:                         NC_039897.vadr.fa.ntf
# BLAST db .nto file for NC_039897 saved in:                         NC_039897.vadr.fa.nto
# fasta sequence file for CDS from NC_039897 saved in:               NC_039897.vadr.cds.fa
# fasta sequence file for translated CDS from NC_039897 saved in:    NC_039897.vadr.protein.fa
# BLAST db .phr file for NC_039897 saved in:                         NC_039897.vadr.protein.fa.phr
# BLAST db .pin file for NC_039897 saved in:                         NC_039897.vadr.protein.fa.pin
# BLAST db .psq file for NC_039897 saved in:                         NC_039897.vadr.protein.fa.psq
# BLAST db .pdb file for NC_039897 saved in:                         NC_039897.vadr.protein.fa.pdb
# BLAST db .pot file for NC_039897 saved in:                         NC_039897.vadr.protein.fa.pot
# BLAST db .ptf file for NC_039897 saved in:                         NC_039897.vadr.protein.fa.ptf
# BLAST db .pto file for NC_039897 saved in:                         NC_039897.vadr.protein.fa.pto
# HMMER model db file for NC_039897 saved in:                        NC_039897.vadr.protein.hmm
# hmmbuild build output (concatenated) saved in:                     NC_039897.vadr.protein.hmmbuild
# binary HMM and p7 HMM filter file saved in:                        NC_039897.vadr.protein.hmm.h3m
# SSI index for binary HMM file saved in:                            NC_039897.vadr.protein.hmm.h3i
# optimized p7 HMM filters (MSV part) saved in:                      NC_039897.vadr.protein.hmm.h3f
# optimized p7 HMM filters (remainder) saved in:                     NC_039897.vadr.protein.hmm.h3p
# hmmpress output file saved in:                                     NC_039897.vadr.hmmpress
# CM file saved in:                                                  NC_039897.vadr.cm
# cmbuild output file saved in:                                      NC_039897.vadr.cmbuild
# binary CM and p7 HMM filter file saved in:                         NC_039897.vadr.cm.i1m
# SSI index for binary CM file saved in:                             NC_039897.vadr.cm.i1i
# optimized p7 HMM filters (MSV part) saved in:                      NC_039897.vadr.cm.i1f
# optimized p7 HMM filters (remainder) saved in:                     NC_039897.vadr.cm.i1p
# cmpress output file saved in:                                      NC_039897.vadr.cmpress
# VADR 'model info' format file for NC_039897 saved in:              NC_039897.vadr.minfo
#
# All output files created in directory ./NC_039897/
```

These files include the FASTA and and their formats are described more
[here](formats.md#build).

Only some of these files will be used by `v-annotate.pl`. These are
the files with the following suffixes:

| file suffix | description | reference |
|--------|-----------------------|-------------|
| `.protein.fa.p{hr,in,sq,db,ot,tf,to}` | BLAST protein database index files, created by `makeblastdb` | binary files, not meant to be human-readable |
| `.fa.n{hr,in,sq,db,ot,tf,to}` | BLAST nucleotide database index files, created by `makeblastdb` | binary files, not meant to be human-readable |
| `.cm` | Infernal 1.1x covariance model file | http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output formats") |
| `.cm.i1{m,i,f,p}` | Infernal 1.1x covariance model index files, created by `cmpress` | binary files, not meant to be human-readable |
| `.hmm` | HMMER 3.x HMM file | http://eddylab.org/software/hmmer/Userguide.pdf ("HMMER profile HMM files" section) |
| `.hmm.h3{m,i,f,p}` | HMMER 3.x covariance model index files, created by `hmmpress` | binary files, not meant to be human-readable |
| `.minfo`  | VADR model info file | [description of format](formats.md#minfo) |

You can use only this model to annotate sequences
with `v-annotate.pl` that are similar to `NC_039897` using the
`--mdir` and `--mkey` options with a command like `v-annotate.pl
--mdir NC_039897 --mkey NC_039897.vadr <fasta-seq-file> <output-directory>`,
explained more [here](annotate.md#options-modelfiles).

Alternatively, you can combine these files together with analogous
files from additional `v-build.pl` runs for other accessions to create
a VADR model library and then use `v-annotate.pl` `--mdir` and  `--mkey` 
options to specify that the library be used. This is explained in more
detail [below](#library). The VADR 1.0 library was created in this
manner, as explained in [another section](#1.0library).

---
## `v-build.pl` command-line options<a name="options"></a>

To get a list of command-line options, execute:

`v-build.pl -h`

This will output the usage and available command-line options. 
Each option has a short description, but additional information on some
of these options can be found below.
For `v-build.pl` the available options are split into eight different categories, 
each explained in their own subsection below.

In the tables describing options below, `<s>` represents a string
and `<n>` represents an integer. 

### `v-build.pl` basic options<a name="options-basic"></a>

The first category of options are the *basic* options:

| ........option........ | explanation | 
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
| `--addminfo <s>` | add arbitrary feature info in file `<s>` to output `.minfo` file, see an example [here](#1.0library-noro) | 
| `--forcelong` | *use at your own risk*; allow long models > 25Kb in length; by default `v-build.pl` will fail for any model more than 25Kb (25,000 nucleotides) because `v-build.pl` will be very slow and the memory requirements of `v-annotate.pl` will be prohibitively large
| `--keep` | keep additional output files that are normally removed |

### `v-build.pl` options for controlling what feature types are stored in the output model info file<a name="options-featuretypes"></a>

By default, only `CDS`, `gene` and `mat_peptide` feature types read from the GenBank feature table file
will be stored in the output `.minfo` file. This default set can be changed using the following three
command line options. For an example of using the `--fadd` option, see the construction of the dengue virus
RefSeq models for the VADR 1.0 model library [here](#1.0library-dengue).

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
model library [here](#1.0library-dengue).

| .......option....... | explanation | 
|--------|-------------| 
| `--qall`   | specify that all qualifiers (except those in `<s>` from `--qskip <s>`) be added to the `.minfo` output file |
| `--qadd <s>`  | add qualifiers listed in `<s>` to the default set, where `<s>` is a comma-separated string with each qualifier separated by a comma with no whitespace |
| `--qftradd <s>`  | specify that the qualifiers listed in `<s2>` from `qadd <s2>` only apply for feature types in the string `<s>`, where `<s>` is a comma-separated string with each qualifier separated by a comma with no whitespace |
| `--qskip <s>`  | do not store information for qualifiers listed in `<s>`, where `<s>` is a comma-separated string with each qualifier separated by a comma with no whitespace; `<s>` may contain qualifiers from the default set, or from other qualifiers (if `--qall` also used) |
| `--noaddgene`  | do not automatically add `gene` qualifiers from `gene` features to any overlapping non-gene features | 
| `--nosplice`   | do not automatically check for GT/AG intron splice sites and add `canon_splice_sites` qualifiers for introns with valid splice sites, an intron is defined as a gap between CDS segments >= `<n>` nucleotides from `--intlen` option, by default `<n>` is `40` |
| `--ssplice`    | exit if any introns in CDS are found that do not have valid GT/AG splice sites, introns defined as explained in `--nosplice` description above |

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
| `--ttbl <n>`   | specify that NCBI translation table `<n>` be used instead of `1` |

Reference on NCBI translation tables:
https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes

### `v-build.pl` options for controlling the `cmbuild` step<a name="options-cmbuild"></a>

Several options exist for controlling the command-line options that will be passed
to Infernal's `cmbuild` program. For more information on these options and how 
they control `cmbuild`, see the Infernal 
User's Guide manual page for `cmbuild` (section 8 of http://eddylab.org/infernal/Userguide.pdf) .

| ......option...... | explanation | 
|--------|-------------| 
| `--cmn <n>` | set the number of seqs for glocal forward profile HMM filter calibration to `<n>` (sets the `cmbuild --EgfN` option), default is to use default `cmbuild` value | 
| `--cmp7ml` | set CM's filter profile HMM as the maximum likelihood profile HMM (sets the `cmbuild --p7ml` option) |
| `--cmere` | set CM relative entropy target bits to position to `<x>` (sets the `cmbuild --ere <x>` option), default is to use default `cmbuild` value |
| `--cmeset` | set CM effective sequence number to '<x>` (sets the `cmbuild --eset <x>` option), default is to use default `cmbuild` value | 
| `--cmemaxseq` | set CM maximum allowed effective sequence # for CM to `<x>` (sets the `cmbuild --emaxseq <x>` option) | 
| `--cmnoh3pri` | do not use `--noh3pri` option with `cmbuild`, to allow CM to use prior from HMMER3 for zero basepair models |
| `--cminfile` | read `cmbuild` options from an input file `<s>`, the contents of the file (after removing newlines) will be supplied directly to `cmbuild` as an options string (possibly with more than one option separated by whitespace) |

### `v-build.pl` options for skipping stages<a name="options-skip"></a>

| ......option...... | explanation | 
|--------|-------------| 
| `--skipbuild` | skip the `cmbuild` step; this is mostly useful for debugging purposes, but also possibly for creating a different `.minfo` file using options like `--fadd` and/or `--qadd` for an already created model without the need to wait for the slow `cmbuild` step | 
| `--onlyurl` | output the url for the GenBank feature table file (or GenBank-format file if `--gb` also used, and exit; possibly helpful if `v-build.pl` is having trouble fetching from GenBank, you can manually download the feature table file from a browswer with the output url and use the downloaded file with `--inft` | 

### `v-build.pl` options for optional output files<a name="options-output"></a>

| .......option....... | explanation | 
|--------|-------------| 
| `--ftrinfo <s>` | output information on the internal data structure used to keep track of features to `<s>`, mainly useful for debugging |
| `--sgminfo <s>` | output information on the internal data structure used to keep track of segments of features to `<s>`, mainly useful for debugging |

### Other `v-build.pl` expert options<a name="options-expert"></a>

| .......option....... | explanation | 
|--------|-------------| 
| `--execname <s>` | in banner and usage output, replace `v-annotate.pl` with `<s>` |
| `--nosig2mat`    | do not treat `sig_peptide` features as `mat_peptide` in model info file |
| `--intlen <n>`   | define intron as any gap >= `<n>` nucleotides between segments in a CDS, only relevant for checking for canonical splice sites, the default value for `<n>` is `40` |


---
## Building a VADR model library<a name="library"></a>

Follow these steps to build a VADR model library:

1. Run `v-build.pl` multiple (`N>1`) times for different accessions.
If your model is larger than 30Kb, or you know that you will use the
`-s` and `--glsearch` options with `v-annotate.pl` as is recommended
for coronavirus annotation
(https://github.com/ncbi/vadr/wiki/Coronavirus-annotation), you should
also use the `--skipbuild` option with `v-build.pl`. If you do, no CM
files will be created so you can skip steps 4 and 8 below and ignore
any other steps involving a `.cm` file.

2. Create a new directory call it `my-vadr-model-dir`, for example.

3. Concatenate all resulting `N` `.vadr.minfo` files into a single
file, call it `my.vadr.minfo`, and move it to the
`my-vadr-model-dir` directory.

4. Concatenate all resulting `N` `.vadr.cm` files into a single file,
call it `my.vadr.cm`, for example, and move it to the
`my-vadr-model-dir` directory.

5. Concatenate all resulting `N` `.vadr.hmm` files into a single file,
call it `my.vadr.hmm`, for example, and move it to the
`my-vadr-model-dir` directory.

6. Concatenate all resulting `N` `.vadr.fa` files into a single file,
call it `my.vadr.fa`, for example, and move it to the
`my-vadr-model-dir` directory, and then create a BLAST nucleotide 
database from it with the command:
```
$VADRBLASTDIR/makeblastdb -dbtype nucl -in my-vadr-model-dir/my.vadr.fa
```

7. Create an Easel index `.ssi` file for the `.vadr.fa` file you created and
moved in step 6 with the command: 
```
$VADREASELDIR/esl-sfetch --index my-vadr-model-dir/my.vadr.fa
```

8. Move all resulting BLAST protein DB files (`.vadr.protein.fa`,
`.vadr.protein.fa.p{hr,in,sq,db,ot,tf,to}`) from all `N` runs into the 
`my-vadr-model-dir` directory.

9. Run `cmpress` on the `my.vadr.cm` file created in step 4 like this:

```
$VADRINFERNALDIR/cmpress my-vadr-model-dir/my.vadr.cm
```

10. Run `hmmpress` on the `my.vadr.hmm` file created in step 5 like this:

```
$VADRHMMERDIR/hmmpress my-vadr-model-dir/my.vadr.hmm
```

You can then use this library with `v-annotate.pl` as follows:
```
v-annotate.pl --mdir my-vadr-model-dir --mkey my.vadr <fasta-file>
<output directory>
```

Or substitute the full paths to `my-vadr-model-dir` if it is not a
subdirectory of the current directory. Optionally, you can supply the
paths to each of the relevant files or directories, possibly if they
are in different directories or are not consistently named, using the `-m`,
`-a`, `-i`, `-n`, and `-x` options as listed
[here](annotate.md#options-modelfiles), with a command like:
```
v-annotate.pl -m my-vadr-model-dir/my.vadr.cm -a my-vadr-model-dir/my.vadr.hmm -i my-vadr-model-dir/my.vadr.minfo -n
my-vadr-model-dir/my.vadr.fa -x my-vadr-model-dir <fasta-file> <output directory>
```

***If you used `--skipbuild` with `v-build.pl`, you will also have
to use the `-s` and `--glsearch` options with `v-annotate.pl`.***

If you ever move `.cm`, `.hmm`, or BLAST `.fa` files into new
directories, make sure you also move the corresponding index files 
(.cm.*, `.hmm.*`, `.fa.*`) along with them.

---
## Practical tutorial: Building a VADR RSV model library<a name="rsv"></a>

The `v-build.pl` program will create a model from a single INSDC
accession and include CDS, gene and mature peptide features. However,
often using a model built from a single accession is not general
enough to allow most high quality sequences from a viral species to pass. For
example, some other sequences may include an extended CDS that has a
different stop codon position from the sequence the model was built
from. It is possible to make VADR models more general but it requires
some manual effort. Below I outline some steps I took when building a
general VADR respiratory syncitial virus (RSV) model library.

1. Determine good reference sequence(s)

There are two RSV subtypes, RSV-A and RSV-B, so the first step I took
was to determine a good reference sequence for each subtype. A good
strategy is often to start with a RefSeq sequence if any are
available. In this case there are two RefSeq sequences which can be
found on the [NCBI virus
resource](#https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), by
selecting the "Search by virus" button and entering "RSV". The top
suggestion will be "Human orthopneumovirus (taxid 11250)", which is a
another name for RSV. At the time of writing, the top two sequences listed in the resulting
list will be RefSeq sequences `NC_038235` (subgroup A) and `NC_001781`
(subgroup B). You can also filter to only RefSeq sequences using the
"Sequence type" filter. 

2. Build initial models from reference sequence(s)

Next, use `v-build.pl` to build the two models, specifying the
`--group` and `--subgroup` options as below:

```
$ v-build.pl --group RSV --subgroup A NC_038235 NC_038235
$ v-build.pl --group RSV --subgroup B NC_001781 NC_001781
```

These commands will take a long time, up to one hour each. 

When they are finished combine the two models into a model library by
following the steps below ([also listed here](#library)).

```
# create a new directory
$ mkdir rsv-models

# concatenate .minfo, .cm .fa and .hmm files:
$ cat NC_001781/*.vadr.minfo > rsv-models/rsv.minfo
$ cat NC_001781/*.vadr.cm > rsv-models/rsv.cm
$ cat NC_001781/*.vadr.fa > rsv-models/rsv.fa
$ cat NC_001781/*.vadr.protein.hmm > rsv-models/rsv.hmm
$ cat NC_038235/*.vadr.minfo >> rsv-models/rsv.minfo
$ cat NC_038235/*.vadr.cm >> rsv-models/rsv.cm
$ cat NC_038235/*.vadr.fa >> rsv-models/rsv.fa
$ cat NC_038235/*.vadr.protein.hmm >> rsv-models/rsv.hmm

# copy the blastdb files:
$ cp NC_001781/*.vadr.protein.fa rsv-models/
$ cp NC_038235/*.vadr.protein.fa* rsv-models/

# prepare the library files:
$ $VADRINFERNALDIR/esl-sfetch --index rsv-models/rsv.fa
$ $VADRINFERNALDIR/cmpress rsv-models/rsv.cm
$ $VADRHMMERDIR/hmmpress rsv-models/rsv.hmm
$ $VADRBLASTDIR/makeblastdb -dbtype nucl -in rsv-models/rsv.fa
```

3. Construct a test set for testing your initial models

INTRODUCE CONCEPT OF IMPROVING/OPTIMIZING MODEL? OR SOME OTHER WORD? INTRODUCE
CONCEPT OF COMMON FAILURE MODES (or something like this).

When evaluating a VADR model it is critical to look at how it performs
when used with `v-annotate.pl` on example sequences. Ideally, you
would know what the expected result is (pass or fail status, and
specific alerts you expect to be or not be reported). For example, if
you have a set of high quality sequences that have been expertly
validated and annotated, you could use this as a positive control -
`v-annotate.pl` should pass all of those sequences and give annotation
matching what is expected. Also, you could intentionally introduce
sequencing errors (insertions, deletions, rearrangments) into high
quality sequences, and checking to see if `v-annotate.pl` detects
those problems and reports them. 

Often times however, the most readily available set of sequences is
simply INSDC sequences of the viral species you are modelling. In this
case, I often take the strategy of selecting a random subset of nearly
full length sequences, evaluating them with `v-annotate.pl` and
manually analyzing the results with this set of training
sequences. While many of these sequences will be high quality and so
are expected to pass, some may have problems with them such as early
stop codons that are either bonafide (biologically real) or the result
of sequencing error or some other artifact. In either case,
`v-annotate.pl` ideally will detect the unusual features and report
them as alerts. For some of these situations, such as an early stop
codon (relative to the reference model) that exists in a significantly
large subset of sequences at a particular position, we may want
`v-annotate.pl` to allow the early stop and not report an alert. For
these cases we will want to modify our model as explained in the
subsequent steps.

The decision to use only full length sequences is debatable, as by
doing it I am assuming that the sequence diversity represented by all
sequences, including partial sequences, is well represented by only
the full length sequences. In other words, if we optimize the
models for only full length sequences, they may perform poorly on
existing partial length sequences that are sufficiently divergent from
all full length sequences. Ideally, you would test on a set of full
length and partial sequences. Alternatively, you could optimize on
full length sequences first, and then as a sanity check, test the
performance of the resulting models on a random subset of partial
sequences to see if any new unexpected failure modes exist.

There are several ways to select a random subset of sequences of a
viral species. One way is to use the NCBI virus resource to download
an accession list of all available sequences, and then select a random
subset from that list using command-line tools, such as the
`esl-selectn` program from Sean Eddy's Easel library that in installed
as part of Infernal with installation of VADR. 

I usually define nearly full length as 95% the length of the RefSeq
sequence, or greater. In this case `NC_038235` is 15222 nucleotides
and `NC_001781` is 15225 nucleotides, so a 95% length cutoff is about
14460 nucleotides.  At the time of writing there are about 5500 nearly
full length RSV INSDC sequences by this definition. For our random subset it is important
to select enough sequences that you should identify any common failure
modes, to decide if they should be addressed by changing the model,
but not too many that your manual analysis will take too long. For
RSV, I chose to use 500 randomly chosen nearly full length sequences
for model improvement.

To download 500 randomly chosen RSV sequences of length
14460 or greater, from the RSV list page you
reached in step 1, use the "Sequence Length" filter to set a minimum
of 14460 and then click the "Download" button, then select
"Nucleotide" under "Sequence data (FASTA format)", then "Download a
randomized subset of all records" then enter "500", click "Next" and 
"Use Default" for the FASTA definition line, then finally click "Download". This should download a file
called something like "sequences_20231006_7993457.fasta". 

4. Use `v-annotate.pl` to validate and annotate test set sequences and
analyze results

The next step is to use our new models HERE HERE HERE 

4. Potentially choose alternative representative sequences

5. Rerun `v-annotate.pl` on test set (if models were updated)

6. Iterative refine models and rerun test set

Methods for updating models to improve performance

1. Add new proteins to blastx database

2. Add alternative features to model info file

3. Add alert exceptions 

4. Build new reference alignments and rebuild CM file

5. Use command-line options to `v-annotate.pl` to make some alerts
fatal or not fatal.



Follow these steps to build a VADR model library:

1. Run `v-build.pl` multiple (`N>1`) times for different accessions.
If your model is larger than 30Kb, or you know that you will use the
`-s` and `--glsearch` options with `v-annotate.pl` as is recommended
for coronavirus annotation
(https://github.com/ncbi/vadr/wiki/Coronavirus-annotation), you should
also use the `--skipbuild` option with `v-build.pl`. If you do, no CM
files will be created so you can skip steps 4 and 8 below and ignore
any other steps involving a `.cm` file.

2. Create a new directory call it `my-vadr-model-dir`, for example.

3. Concatenate all resulting `N` `.vadr.minfo` files into a single
file, call it `my.vadr.minfo`, and move it to the
`my-vadr-model-dir` directory.

4. Concatenate all resulting `N` `.vadr.cm` files into a single file,
call it `my.vadr.cm`, for example, and move it to the
`my-vadr-model-dir` directory.

5. Concatenate all resulting `N` `.vadr.hmm` files into a single file,
call it `my.vadr.hmm`, for example, and move it to the
`my-vadr-model-dir` directory.

6. Concatenate all resulting `N` `.vadr.fa` files into a single file,
call it `my.vadr.fa`, for example, and move it to the
`my-vadr-model-dir` directory, and then create a BLAST nucleotide 
database from it with the command:
```
$VADRBLASTDIR/makeblastdb -dbtype nucl -in my-vadr-model-dir/my.vadr.fa
```

7. Create an Easel index `.ssi` file for the `.vadr.fa` file you created and
moved in step 6 with the command: 
```
$VADREASELDIR/esl-sfetch --index my-vadr-model-dir/my.vadr.fa
```

8. Move all resulting BLAST protein DB files (`.vadr.protein.fa`,
`.vadr.protein.fa.p{hr,in,sq,db,ot,tf,to}`) from all `N` runs into the 
`my-vadr-model-dir` directory.

9. Run `cmpress` on the `my.vadr.cm` file created in step 4 like this:

```
$VADRINFERNALDIR/cmpress my-vadr-model-dir/my.vadr.cm
```

10. Run `hmmpress` on the `my.vadr.hmm` file created in step 5 like this:

```
$VADRHMMERDIR/hmmpress my-vadr-model-dir/my.vadr.hmm
```

You can then use this library with `v-annotate.pl` as follows:
```
v-annotate.pl --mdir my-vadr-model-dir --mkey my.vadr <fasta-file>
<output directory>
```

Or substitute the full paths to `my-vadr-model-dir` if it is not a
subdirectory of the current directory. Optionally, you can supply the
paths to each of the relevant files or directories, possibly if they
are in different directories or are not consistently named, using the `-m`,
`-a`, `-i`, `-n`, and `-x` options as listed
[here](annotate.md#options-modelfiles), with a command like:
```
v-annotate.pl -m my-vadr-model-dir/my.vadr.cm -a my-vadr-model-dir/my.vadr.hmm -i my-vadr-model-dir/my.vadr.minfo -n
my-vadr-model-dir/my.vadr.fa -x my-vadr-model-dir <fasta-file> <output directory>
```

***If you used `--skipbuild` with `v-build.pl`, you will also have
to use the `-s` and `--glsearch` options with `v-annotate.pl`.***

If you ever move `.cm`, `.hmm`, or BLAST `.fa` files into new
directories, make sure you also move the corresponding index files 
(.cm.*, `.hmm.*`, `.fa.*`) along with them.

---
## How the VADR 1.0 model library was constructed <a name="1.0library"></a>

The VADR 1.0 library was built with version 1.0 of VADR. It is included here as 
an example of how to build a VADR library, but also so it can be reproduced, because it is the model library
used in the [paper on VADR 1.0](../README.md#reference) (https://doi.org/10.1186/s12859-020-3537-3).
If you want to reproduce it exactly, you'll want to install
version 1.0. The install script for v1.0 is here:

https://raw.githubusercontent.com/ncbi/vadr/1.0/vadr-install.sh

Additionally, if the RefSeq annotation for any of these 197 VADR
models has changed since October 2019, then it may not be able to
identically reproduce the VADR 1.0 model library using the steps
outlined above. This is because `v-build.pl` fetches the current
RefSeq annotation data from GenBank when it is run. If necessary,
email eric.nawrocki@nih.gov for additional files
needed to reproduce the library exactly.

Also note that the library has changed since version 1.0. For example,
the default set of models included with version 1.1 has 205 total
models, not 197. To see a list of changes, see the `RELEASE-NOTES.txt`
file in the directory pointed to by your `$VADRMODELDIR` environment
variable after installing VADR. To reproduce the construction of the
1.1 library, you would run similar steps to those below but also
adding the additional models listed in the `RELEASE-NOTES.txt` file.

Additionally, as of version 1.2, the *Caliciviridae* models and
*Flaviviridae* models have been split up into two different model
sets, but both are installed by the VADR install script
`vadr-install.sh`. The *Caliciviridae* models are used by default with
`v-annotate.pl`. To use the *Flaviviridae* models, use the options `--mkey
flavi --mdir $VADRMODELDIR/vadr-models-flavi`.

### Models in the VADR 1.0 library

The VADR 1.0 model library consists of 197 VADR models. Nine of these
are Norovirus RefSeq models, listed in
[vadr/documentation/build-files/1p0-models/norovirus.9.list](build-files/1p0-models/norovirus.9.list).
Four of these are Dengue virus RefSeq models, listed in
[vadr/documentation/build-files/1p0-models/dengue.4.list](build-files/1p0-models/dengue.4.list).
Eight of these are Hepatitis C virus Refseq models, listed in
[vadr/documentation/build-files/1p0-models/hcv.8.list](build-files/1p0-models/hcv.8.list).
The remaining 173 are additional *Caliciviridae* and *Flaviviridae*
RefSeq models, listed in
[vadr/documentation/build-files/1p0-models/non-noro-dengue-hcv.173.list](build-files/1p0-models/non-noro-dengue-hcv.173.list).

### Building the VADR 1.0 library Norovirus models <a name="1.0library-noro"></a>

To build models for each of the nine norovirus RefSeqs listed in 
[vadr/documentation/build-files/1p0-models/norovirus.9.list](build-files/1p0-models/norovirus.9.list),
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
[vadr/documentation/build-files/1p0-models/norovirus.9.build.sh](build-files/1p0-models/norovirus.9.build.sh)
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
[vadr/documentation/build-files/1p0-models/NC_029646.addminfo](build-files/1p0-models/NC_029646.addminfo)
and 
[vadr/documentation/build-files/1p0-models/NC_039475.addminfo](build-files/1p0-models/NC_039475.addminfo).
The `NC_029646.addminfo` file is in the [`.minfo` format](formats.md#minfo) and looks like this:

```
MODEL NC_029646 
FEATURE NC_029646 type:"CDS" coords:"5085..6692:+" xmaxins_exc:"297:36"
```

This file specifies that the additional `<key>:<value>` pair of
`xmaxins_exc:"297:36"` be added to the CDS feature with coordinates
`5085..6692:+`. (The VADR coordinate string format is described
[here](formats.md#coords)). The `NC_039475.addminfo` file is
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
[vadr/documentation/build-files/1p0-models/dengue.4.build.sh](build-files/1p0-models/dengue.4.build.sh)
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

#### Manually edit the `.minfo` files for the Dengue models

Not all the annotation in the Dengue RefSeqs as of October 2019 is in
the VADR 1.0 library. To make your `.minfo` files match those in the
library exactly do the following:

1. Remove the `ncRNA` features with product names `sfRNA2`, `sfRNA3`,
   and `sfRNA4` from the `.minfo` files of `NC_001477`, `NC_001474`,
   and `NC_001475`.

2. Remove the `ncRNA` features with product names `sfRNA3` and
   `sfRNA4` from the `.minfo` file of `NC_002640`.

3. Remove all occurrences of the substring `gene:"-"` (only) from
   `FEATURE` lines in the `.minfo` files of `NC_001477`, `NC_001474`,
   `NC_001475` and `NC_002640`.

### Building the VADR 1.0 library Hepatitis C virus models <a name="1.0library-hcv"></a>

To build models for each of the eight Hepatitis C RefSeqs listed in 
[vadr/documentation/build-files/1p0-models/hcv.8.list](build-files/1p0-models/hcv.8.list),
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
[vadr/documentation/build-files/1p0-models/hcv.8.build.sh](build-files/1p0-models/hcv.8.build.sh)
contains these commands.)

The `--group` and `--subgroup` options are used in a similar way to how they were used
to build the norovirus models.

### Building the VADR 1.0 library models for the other 173 *Caliciviridae* and *Flaviviridae* viral species <a name="1.0library-173"></a>

To build models for the other 173 *Caliciviridae* and *Flaviviridae* models listed in
[vadr/documentation/build-files/1p0-models/non-noro-dengue-hcv.173.list](build-files/1p0-models/non-noro-dengue-hcv.173.list)

Simply run `v-build.pl` from VADR v1.0 using default parameters for each accession.
For example:

```
  v-build.pl NC_034444 NC_034444 
```

(The shell script
[vadr/documentation/build-files/1p0-models/non-noro-dengue.hcv.173.build.sh](build-files/1p0-models/non-noro-dengue.hcv.173.build.sh)
will execute these 173 commands.)

Each of these commands takes roughtly between 10 minutes and an hour. 

### Concatenating files for the 197 models to create the VADR library <a name="1.0library-concat"></a>

After completing the steps above to make the 197 models, you can make
the VADR 1.0 library by following the instructions for creating a VADR
library starting at step 2 [here](#library).

---
#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.

