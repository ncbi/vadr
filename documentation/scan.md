# <a name="top"></a> `v-scan.pl` example usage and command-line options

* [Quickstart `v-scan.pl` examples](#quickstart)
* [The `v-scan.pl` config file](#config)
  * [adding libraries to the config file](#add2config)
  * [printing information from the config file](#printinfo)
* [Walk-throughs of `v-scan.pl` examples](#longwalk)
  * [running `v-scan.pl` when all sequences are expected to match one unknown model library (mode 1)](#mode1)
  * [running `v-scan.pl` when sequences may match multiple models libraries (mode 2)](#mode2)
  * [running `v-scan.pl` when all sequences are expected to match one known model library (mode 3)](#mode3)
* [`v-scan.pl` command-line options](#options)
  * [basic options](#options-basic)
  * [options for specifying which model libraries to use](#options-libraries)
  * [options controlling the random sampling of sequences](#options-sampling)
  * [options for listing information from the config model or about models and exiting](#options-listing)
* [Using the list options](#listexamples)

---

## Quickstart `v-scan.pl` examples <a name="quickstart"></a>

`v-scan.pl` can be run in 3 modes:

1. model library is unknown but all sequences are expected to match to
   the same library, run like: `v-scan.pl <fasta file to annotate>
   <output directory to create>`. For example:

   ```
   v-scan.pl $VADRSCRIPTSDIR/documentation/scan-files/n5.fa vs-n5
   ```

   After running the output files will be in the newly created directory
   `vs-n5`.

2. model libraries are unknown and sequences may match to different
libraries, same as above but add the `-m` option:

   ```
   v-scan.pl -m $VADRSCRIPTSDIR/documentation/scan-files/m5.fa vs-m5
   ```

3. model library is known and all sequences are expected to match to
it, same as 1 but use the `--only` option:

   ```
   v-scan.pl --only norovirus $VADRSCRIPTSDIR/documentation/scan-files/n5.fa vs-n5-only
   ```

---

## `v-scan.pl` config file<a name="config"></a>

<a name="config"></a> `v-scan.pl` can be used to annotate sequences
that match to one or more VADR model libraries. `v-scan.pl`
determines which model library to use for the input sequences and then
calls `v-annotate.pl` to annotate all sequences that match to that
model library. It will supply `v-annotate.pl` with the command-line
options specific for that model library read from the input config
file.

 `v-scan.pl` requires a 'config' file that lists information on the
model libraries it will use. The config file that will be used 
will be stored in the `$VADRCONFIGFILE` environment variable after the
[installation procedure](install.md#top). If you want you can modify
this file or modify a copy of it, and specify that a different file
`<s>` be used by using the `-c <s>` option to `v-scan.pl` or by
changing the value of `$VADRCONFIGFILE` to point to a different file. 

Here is the config file that is included with VADR and is used by
default in [vadr/vadr.config](../vadr.config) with
comment lines removed for brevity (all lines that begin with a `#` are
comment lines):

```
dengue    $VADRINSTALLDIR/vadr-models-flavi    --split --cpu 1 --group Dengue --nomisc --noprotid --mkey flavi -r
hcv       $VADRINSTALLDIR/vadr-models-flavi    --split --cpu 1 -r --mkey flavi --group HCV
flavi     $VADRINSTALLDIR/vadr-models-flavi    --split --cpu 1 -r --nomisc
norovirus $VADRINSTALLDIR/vadr-models-calici   --split --cpu 1 --group Norovirus --nomisc --noprotid --mkey calici -r
calici    $VADRINSTALLDIR/vadr-models-calici   --split --cpu 1 -r --nomisc 
sarscov2  $VADRINSTALLDIR/vadr-models-sarscov2 --split --cpu 4 -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --glsearch
corona    $VADRINSTALLDIR/vadr-models-corona   --split --cpu 1 -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --glsearch 
flu       $VADRINSTALLDIR/vadr-models-flu      --split --cpu 2 -r --atgonly --alt_fail extrant5,extrant3 --xnocomp --nomisc --forcegene
rsv       $VADRINSTALLDIR/vadr-models-rsv      --split --cpu 1 -r --xnocomp --nomisc 
mpxv      $VADRINSTALLDIR/vadr-models-mpxv     --split --cpu 1 --glsearch --minimap2 -s -r --nomisc --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --s_overhang 150
```
(Each line prefixed with `#` is a comment line and is ignored by
`v-scan.pl`.) All other lines have 3 or more fields:

| idx      | field               | description |
|----------|---------------------|-------------| 
|   1      | `<options key>`     | name for this library, a unique key that will be used for naming `v-scan.pl` output files, cannot contain whitespace |
|   2      | `<model directory>` | path to the model directory that includes all model files for this library, the same model directory can be used for multiple `<options key>` values, cannot contain whitespace |
| 3 to end | `<options string>`  | the `v-annotate.pl` options that should be used for sequences matching this `<options key>` library during the annotation stage of `v-scan.pl`; this must contain `--mkey <s>` if the model files in the `<model directory>` are named with a key other than the `<options key>`, an example is the use of `--mkey calici` for the `norovirus` `<options key>` in the example file above; may contain whitespace |

Note that in the above example config file, both `dengue` and `hcv`
`<options key>` values use the `flavi` model library. We can tell this
because `--mkey flavi` exists in the `<options string>` for the lines
beginning with `dengue` and `hcv`, *and* because the `<model directory>` is the
same for all three of `dengue`, `hcv` and `flavi`. Similarly
`norovirus` uses the `calici` library.

Multiple `<options string>` values can use the same model libraries
because sometimes we may want to use different `v-annotate.pl` options
for different models within those libraries. In the example of
`dengue`, `hcv`, and `flavi`, you can see that `dengue` uses the
`--noprotid` and `--group Dengue` options, where as `hcv` and `flavi`
use different options. All sequences that match to `dengue` will be
annotated by `v-annotate.pl` using the `dengue` options, all sequences
that match to `hcv` will be annotated using the `hcv` options, and all
sequences that match to `flavi` will be annotated using the `flavi`
options.

The way that sequences are *matched* to an `<options key>` is as
follows: in the `v-scan.pl` classification stage, each sequence is
scanned against each unique model library from the config file. In
the example config file, this is only two model libraries:

1. the `flavi` model library, with files named with the key `flavi`
(e.g. `flavi.minfo`) in the directory
`$VADRINSTALLDIR/vadr-models-flavi`.

2. the `calici` model library, with files named with the key `calici` 
(e.g. `calici.minfo`) in the directory
`$VADRINSTALLDIR/vadr-models-calici`.

The `dengue` and `hcv` `<options key>` values use the `flavi` library due to
the `--mkey flavi` in their `<options string>`, and `norovirus` 
uses the `calici` library due to `--mkey calici` in its `<options string>`.

Then, when parsing the output for the scan against the `flavi`
library, sequences are matched to `dengue` or `hcv` by
checking if the best matching model for each sequence matches to
`dengue` or `hcv` (`flavi` itself is the default `<options key>`, used for any
sequence that matches neither `dengue` nor `hcv`, as described below). A
model matches to an `<options key>` if
its `model name`, `group` or `subgroup` is identical to that `<options key>`
*after lowercasing and removing all
special characters* from the name, group or subgroup. For example, if a
sequence's best matching model is `NC_001477` which has `group` defined
as `Dengue` in the `flavi.minfo` file (relevant line below)
```
MODEL NC_001477 blastdb:"NC_001477.vadr.protein.fa" group:"Dengue" length:"10735" subgroup:"1"
```
then that sequence will match to `dengue` and
that sequence will then be annotated with `v-annotate.pl` using the
`dengue` `<options string>`. Or, if a sequence matched to a model
named `HCV!` then it would match to `hcv` because `HCV!` becomes `hcv`
after making it lowercase and removing all special (non-alphanumeric)
characters. Any sequence that matches best to a model in the `flavi`
library that does not match either `dengue` or `hcv` in this way will
match to `flavi` and `v-annotate.pl` will be used with the `flavi`
options string to annotate it.

Similarly, when parsing the output for the scan against the `calici`
library, sequences are matched to either `norovirus` or `calici` in
the same way. 

Model libraries do not have to be nested in this way. Each `<options
key>` in the config file can pertain to its own unique model library. 

#### Modifying the `--split` and `--cpu <n>` options 

Note that in the default config file shown above, each `<options
string>` includes the `--split` and `--cpu <n>` options with `<n>`
varying between `1`, `2` and `4`. These options are used to
parallelize processing in `v-annotate.pl` across `<n>` threads. If you
prefer to run non-threaded, you can remove these two options from all
the lines of the config file. You should also potentially change the
number of threads that are used for each `<options key>` by modifying
`<n>` after considering how much RAM you have on your computer The
suggested amount of RAM **per thread** for each virus is:

| options key | recommended RAM per thread |
|-------------|----------------------------|
| sarscov2    | 2Gb |
| flu         | 4Gb |
| dengue, hcv, flavi, norovirus, calici, corona, mpxv | 8Gb |
| rsv         |  16Gb |

#### Adding libraries to the config file<a name="add2config"></a>

You can download or build additional vadr model libraries and add them
to the config file or make your own config file. To use a different
config file `<s>` use the `-c <s>` option. The list of available VADR
models and an example config file that uses them all is 
[here](https://github.com/ncbi/vadr/wiki/Available-VADR-model-files).

#### Printing information from the config file<a name="printinfo">

Several command-line options exist for outputting information from the 
config file, and on the models that are in the libraries listed in the
config file: `--l_all`, `--l_lib <s>`, `--l_dir`, `--l_opt` and
`--l_mdl`. Examples of these can be found [below.](#listexamples)

---

## Walk-throughs of `v-scan.pl` examples <a name="longwalk"></a>

This section includes more detailed information on how to use
`v-scan.pl` and understand its output.

To determine the command-line usage of 
`v-scan.pl` (or any VADR script), use the `-h` option, like this:

```
v-scan.pl -h 
```

You'll see something like the following output:
```
# v-scan.pl :: scan and annotate sequences against VADR model libraries 
# VADR 1.7 (Sep 2025)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Tue Sep  2 10:39:19 2025
#
Usage: v-scan.pl [-options] <fasta file to annotate> <output directory to create>
```

The first few lines are the banner which show the name of the VADR
script being run along with the version and release date. This is
followed by the time and date the command was executed.  The `Usage:`
line details the expected command line arguments.  `v-scan.pl` takes
as input two command line arguments, a fasta file with sequences to
analyze and annotate (`<fasta file to annotate>`) and the name of the
output directory you want it to create (`<output directory to
create>`) and populate with output files.

After that comes a list of all available command-line options. These
are explained in more detail [below](#options).

---

### Running `v-scan.pl` (mode 1)<a name="mode1"></a>

<a name="examplebasic"></a>Below is an example `v-scan.pl` command
run in mode 1, where all sequences are expected to match to a single
model library (using our default config file this would mean either
'flavi' or 'calici') but that library is unknown. For this example we
will use the sequence file [vadr/documentation/scan-files/n5.fa](scan-files/n5.fa)
that includes 5 norovirus sequences, and we will create the directly 
`vs-n5`:

```
v-scan.pl $VADRSCRIPTSDIR/documentation/scan-files/n5.fa vs-n5
```

The standard output of `v-scan.pl` that is printed to the screen
(which is also output to the `.log` output file) begins with the
banner and date again followed by a list of relevant environment
variables, the command line arguments used and any command line
options used:

```
# date:             Tue Sep  2 14:16:40 2025
# $VADRCONFIGFILE:  /home/nawrocki/vadr-install-dir/vadr.config
# $VADRSCRIPTSDIR:  /home/nawrocki/vadr-install-dir
#
# sequence file:     /home/nawrocki/vadr-install-dir/documentation/scan-files/n5.fa
# output directory:  vs-n5
```

No command line options were used in our example output, but if they
were information on them would have appeared after the `output
directory` line.

Next, information is output about each step the script is proceeding
through. When each step is completed, the elapsed time for that step
is output. 

`v-scan.pl` will use the default VADR config file
(`$VADRMODELDIR/default.vadr.config`) to analyze the sequences in
`n5.fa`, and will create an output directory named `vs-5` and populate
it with [many output files](formats.md#annotate).

The output of `v-scan.pl` lists the steps it takes: 

```
# Validating input                                             ... done. [    0.0 seconds]
# Sampling 3 sequences to use for classification               ... done. [    0.0 seconds]
# Scanning 3 sequences against flavi  library                  ... done. [    3.3 seconds]
# Scanning 3 sequences against calici library                  ... done. [    3.8 seconds]
```

After validating the input, `v-scan.pl` samples 3 sequences to use for
the classification stage. Only 3 sequences are used to make this stage
more efficient based on the assumption is that all the sequences will
be for the same model library. If `v-scan.pl` determines that more
than one model library is matched by the sampled sequences, then it
will fail with an error message indicating which libraries were
matched.  The `--s_nseq <n>` option will specify that `<n>` sequences
be sampled instead of 3, and the `--all` option will turn off sampling
and all sequences be used for classification. An example of allowing
multiple libraries to match with the `-m` option is below.

After the classification stage, all of the sequences will be annotated
using the matching library, which in this case is the `norovirus`
library:

```
# Annotating 5 sequences with norovirus model library          ... done. [   31.7 seconds]
# Generating tabular output                                    ... done. [    0.0 seconds]
```

```
# Summary of sequences matching norovirus:
#
#                                      num   num   num
#idx  model      group      subgroup  seqs  pass  fail
#---  ---------  ---------  --------  ----  ----  ----
1     NC_039475  Norovirus  GII          2     2     0
2     NC_039477  Norovirus  GII          2     2     0
3     NC_039476  Norovirus  GII          1     1     0
#---  ---------  ---------  --------  ----  ----  ----
-     *all*      -          -            5     5     0
-     *none*     -          -            0     0     0
#---  ---------  ---------  --------  ----  ----  ----
#
# Zero alerts reported for seqs matching norovirus.
```

And finally the output files are listed:

```
# norovirus library output printed to screen saved in:                              vs-n5.norovirus.vadr.log
# norovirus library list of executed commands saved in:                             vs-n5.norovirus.vadr.cmd
# norovirus library list and description of all output files saved in:              vs-n5.norovirus.vadr.filelist
# norovirus library esl-seqstat -a output for input fasta file saved in:            vs-n5.norovirus.vadr.seqstat
# norovirus library 5 column feature table output for passing sequences saved in:   vs-n5.norovirus.vadr.pass.tbl
# norovirus library 5 column feature table output for failing sequences saved in:   vs-n5.norovirus.vadr.fail.tbl
# norovirus library list of passing sequences saved in:                             vs-n5.norovirus.vadr.pass.list
# norovirus library list of failing sequences saved in:                             vs-n5.norovirus.vadr.fail.list
# norovirus library list of alerts in the feature tables saved in:                  vs-n5.norovirus.vadr.alt.list
# norovirus library fasta file with passing sequences saved in:                     vs-n5.norovirus.vadr.pass.fa
# norovirus library fasta file with failing sequences saved in:                     vs-n5.norovirus.vadr.fail.fa
# norovirus library per-sequence tabular classification summary file saved in:      vs-n5.norovirus.vadr.sqc
# norovirus library per-model tabular summary file saved in:                        vs-n5.norovirus.vadr.mdl
# norovirus library per-sequence tabular annotation summary file saved in:          vs-n5.norovirus.vadr.sqa
# norovirus library per-feature tabular summary file saved in:                      vs-n5.norovirus.vadr.ftr
# norovirus library per-model-segment tabular summary file saved in:                vs-n5.norovirus.vadr.sgm
# norovirus library per-alert tabular summary file saved in:                        vs-n5.norovirus.vadr.alt
# norovirus library alert count tabular summary file saved in:                      vs-n5.norovirus.vadr.alc
# norovirus library alignment doctoring tabular summary file saved in:              vs-n5.norovirus.vadr.dcr
# norovirus library replaced stretches of Ns summary file (-r) saved in:            vs-n5.norovirus.vadr.rpn
#
# Output printed to screen saved in:                   vs-n5.vadr.log
# List of executed commands saved in:                  vs-n5.vadr.cmd
# List and description of all output files saved in:   vs-n5.vadr.filelist
# per-model library tabular summary file saved in:     vs-n5.vadr.lib
#
# All output files created in directory ./vs-n5/
#
# Elapsed time:  00:00:38.82
#                hh:mm:ss
# 
[ok]
```

All of these files were created in the newly created directory
`vs-n5`. The files that include `norovirus` in their name were created
by `v-annotate.pl` and more information on those files and their
formats can be found [here](formats.md#annotate). The [`v-annotate.pl`
documentation](formats.md#lib) also includes a description of several of these files.

The final four files are the [`.log` file](formats.md#log), which is
the same as the standard output printed to the screen currently being
discussed, the [`.cmd` file](formats.md#cmd), and the [`.filelist`
file](formats.md#filelist) which lists the output files created by
`v-scan.pl`, and finally the [`.lib` file](formats.md#lib)
which explains how many sequences matched to each library in the
classification stage. 

---

### Running `v-scan.pl` (mode 2)<a name="mode2"></a>

If there may be sequences that match to multiple libraries in the
input fasta file, use the `-m` option, like this:

```
v-scan.pl -m $VADRSCRIPTSDIR/documentation/scan-files/m5.fa vs-m5
```

This will generate:
```
# Validating input                                             ... done. [    0.0 seconds]
# Scanning 5 sequences against flavi  library                  ... done. [    4.8 seconds]
# Scanning 5 sequences against calici library                  ... done. [    4.1 seconds]
# Annotating 1 dengue    sequences                             ... done. [   40.4 seconds]
# Annotating 1 flavi     sequences                             ... done. [   45.3 seconds]
# Annotating 2 norovirus sequences                             ... done. [   27.7 seconds]
# Annotating 1 calici    sequences                             ... done. [   17.3 seconds]
# Generating tabular output                                    ... done. [    0.0 seconds]
#
# Summary of sequences matching each library:
#
#     options    model    num
#idx  key        key     seqs
#---  ---------  ------  ----
1     dengue     flavi      1
2     hcv        flavi      0
3     flavi      flavi      1
4     norovirus  calici     2
5     calici     calici     1
```

You may notice a few differences between the output with `-m`. First,
there is no sampling step, all sequences must be used for the
classification step because with `-m` multiple libraries may be matched and
we want to determine which library to use for annotating each
sequence. 

Also, following the summary of sequences matching each library, you
will see per-library summary statistics for each of the four libraries
matched, and a list of output sequences for each library.

---

### Running `v-scan.pl` (mode 3)<a name="mode3"></a>

If you know which library your sequences will match to, you can use
the `--only` option, like this:

```
v-scan.pl --only norovirus $VADRSCRIPTSDIR/documentation/scan-files/n5.fa vs-n5-only
```

```
# Validating input                                             ... done. [    0.0 seconds]
# Annotating 5 sequences with norovirus model library          ... done. [   31.6 seconds]
#
# Summary of sequences matching norovirus:
#
#                                      num   num   num
#idx  model      group      subgroup  seqs  pass  fail
#---  ---------  ---------  --------  ----  ----  ----
1     NC_039475  Norovirus  GII          2     2     0
2     NC_039477  Norovirus  GII          2     2     0
3     NC_039476  Norovirus  GII          1     1     0
#---  ---------  ---------  --------  ----  ----  ----
-     *all*      -          -            5     5     0
-     *none*     -          -            0     0     0
#---  ---------  ---------  --------  ----  ----  ----
```
Note that the sampling step and classification step performed in mode
1 are skipped. This makes it slightly more efficient to use `--only`
if you know the library.

If you want to restrict the libraries that can be matched, but not
specify a single one, you can use `--only` with a list of comma
separated library keys, like this: `v-scan.pl --only norovirus,flavi
n5.fa vs-n5-only2`. Alternatively, you can list libraries that you want
`v-scan.pl` to skip (ignore) using the `--skip` option like this:
`v-scan.pl --skip dengue,flavi n5.fa vs-n5-skip`.

---

## `v-scan.pl` command-line options<a name="options"></a>

To get a list of command-line options, execute:

`v-scan.pl -h`

This will output the usage and available command-line options. 
Each option has a short description, but additional information on some
of these options can be found below.
For `v-scan.pl` the available options are split into four different categories, 
each explained in their own subsection below.

In the tables describing options below, `<s>` represents a string,
`<x>` indicates a floating point number and `<n>` represents an
integer. 

### `v-scan.pl` basic options<a name="options-basic"></a>

| ......option.... | explanation | 
|------------------|-------------|
| `-f`             | if `<output directory>` already exists, then using this option will cause it to be overwritten, otherwise the progam exits in error |
| `-m`             | multiple-library mode, allow matches to multiple model libraries, without this option matches to only one library are allowed and the program will exit if matches to multiple libraries are found |
| `-c <s>`         | use the config file `<s>` instead of the config file in $VADRCONFIGFILE |
| `-v`             | *verbose* mode: all commands will be output to standard output as they are run | 
| `--first`        | specify that if a sequence matches to more than one library, use the first one; by default the higher scoring match is used |
| `--cpu <n>`      | parallelize classification stage across <n> CPU workers, requires `-m`, only impacts the initial classification stage, parallelization of annotation stage can be controlled by adding `--split --cpu <n>` to options strings in the config file |
| `--lone`         | exit if at least one sequence matches to multiple libraries |
| `--origfa`       | do not copy the input fasta file into output directory prior to analysis, use the original |
| `--keep`         | keep [additional `v-annotate.pl` output files](formats.md#annotate-keep) that are normally removed |

### `v-scan.pl` options for specifying which model libraries to use<a name="options-libraries"></a>
| ..........option.......... | explanation | 
|--------|-------------| 
| `--only <s>`      | only use the model library(ies) with option keys (e.g. `flavi`) listed in the comma separated string `<s>`, all option keys must exist in config file | 
| `--skip <s>`      | do not use the model library(ies) with option keys (e.g. `flavi`) listed in the comma separated string `<s>`,  all option keys must exist in config file | 

### <a name="options-sampling"></a> `v-scan.pl` options related to the random sampling of sequences for determining model library to use (sampling is turned off if `-m` is used or only one model library is being used)

| ............option............ | explanation | 
|----------------------------|-------------| 
| `--all`              | do not sample, pick model library(ies) based on all sequences (automatically turned on if `-m` used) | 
| `--s_nseq <n>`       | set the number of sequences to sample to `<n>`, default value is `3` |
| `--s_beg`            | sample sequences from the beginning of the sequence file, not randomly |
| `--s_seed <n>`       | set the random number generator seed to `<n>`, default value is `181` |

### `v-scan.pl` options for listing information from the config file or about models and exiting<a name="options-listing"></a>

| ............option............ | explanation | 
|----------------------------|-------------| 
| `--l_all`       | list information about all models, model directories, and options strings in the config file and exit |
| `--l_lib <s>`   | list all information about the model library for options key `<s>` (e.g. `flavi`) in the config file and exit |
| `--l_dir`       | list all model directories in the config file and exit |
| `--l_opt`       | list `v-annotate.pl` options for each option key in the config file and exit |
| `--l_mdl`       | list information about all the models in all libraries in the config file and exit |

---

## Using the list options <a name="listexamples"></a>

The `v-scan.pl` options listed above beginning with `-l` can be
useful for listing information about the config file and the model
libraries listed in the config file. Here are some examples of using
these options when the environment variable `$VADRCONFIGFILE` points
to the default config file in `$VADRINSTALLDIR/default.vadr.config`,
which it should be default after following the [installation
instructions](install.md#top):

---

List the model directories in the config file, for all libraries:

```
v-scan.pl --l_dir
```

```
############################################################
#
# VADR 1.7 (Sep 2025)
#
# config file: /home/nawrocki/vadr-install-dir/default.vadr.config
#
# Model library directory information:
#
#options key  model key  model dir
#-----------  ---------  ---------
dengue        flavi      /home/nawrocki/vadr-install-dir/vadr-models-flavi
hcv           flavi      /home/nawrocki/vadr-install-dir/vadr-models-flavi
flavi         "          /home/nawrocki/vadr-install-dir/vadr-models-flavi
norovirus     calici     /home/nawrocki/vadr-install-dir/vadr-models-calici
calici        "          /home/nawrocki/vadr-install-dir/vadr-models-calici
sarscov2      "          /home/nawrocki/vadr-install-dir/vadr-models-sarscov2
corona        "          /home/nawrocki/vadr-install-dir/vadr-models-corona
flu           "          /home/nawrocki/vadr-install-dir/vadr-models-flu
rsv           "          /home/nawrocki/vadr-install-dir/vadr-models-rsv
mpxv          "          /home/nawrocki/vadr-install-dir/vadr-models-mpxv
#

#
```

---

List the `v-annotate.pl` options in the config file, for all libraries:

```
v-scan.pl --l_opt
```

```
# Options information:
#
#options key  v-annotate.pl options
#-----------  ---------------------
dengue        --split --cpu 1 --group Dengue --nomisc --noprotid --mkey flavi -r
hcv           --split --cpu 1 -r --mkey flavi --group HCV
flavi         --split --cpu 1 -r --nomisc
norovirus     --split --cpu 1 --group Norovirus --nomisc --noprotid --mkey calici -r
calici        --split --cpu 1 -r --nomisc
sarscov2      --split --cpu 4 -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --glsearch
corona        --split --cpu 1 -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --glsearch
flu           --split --cpu 2 -r --atgonly --alt_fail extrant5,extrant3 --xnocomp --nomisc --forcegene
rsv           --split --cpu 1 -r --xnocomp --nomisc
mpxv          --split --cpu 1 --glsearch --minimap2 -s -r --nomisc --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --s_overhang 150
#
```

---

List all models in the model libraries in the config file:

```
v-scan.pl --l_mdl
```
```
# List of models in each library:
#
#idx   model key  options key  model name  length  group          subgroup
#----  ---------  -----------  ----------  ------  -------------  --------
1.1    flavi      dengue       NC_001477    10735  Dengue         1       
1.2    flavi      dengue       NC_001474    10723  Dengue         2       
1.3    flavi      dengue       NC_001475    10707  Dengue         3       
1.4    flavi      dengue       NC_002640    10649  Dengue         4       
1.5    flavi      hcv          NC_004102     9646  HCV            1       
1.6    flavi      hcv          NC_038882     9599  HCV            1       
1.7    flavi      hcv          NC_009823     9711  HCV            2       
1.8    flavi      hcv          NC_009824     9456  HCV            3       
1.9    flavi      hcv          NC_009825     9355  HCV            4       
1.10   flavi      hcv          NC_009826     9343  HCV            5       
1.11   flavi      hcv          NC_009827     9628  HCV            6       
1.12   flavi      hcv          NC_030791     9443  HCV            7       
1.13   flavi      "            NC_040815     8684  -              -       
1.14   flavi      "            NC_040788    10311  -              -       
1.15   flavi      "            NC_040776    10794  -              -       
..snip..
1.153  flavi      "            NC_001837     9550  -              -       
1.154  flavi      "            NC_001710     9392  -              -       
1.155  flavi      "            NC_001461    12573  -              -       
1.156  flavi      "            NC_031327    10588  -              -       
#----  ---------  -----------  ----------  ------  -------------  --------
2.1    calici     norovirus    NC_001959     7654  Norovirus      GI      
2.2    calici     norovirus    NC_008311     7382  Norovirus      GV      
2.3    calici     norovirus    NC_029645     7313  Norovirus      GIII    
2.4    calici     norovirus    NC_029646     7518  Norovirus      GII     
..snip..
5.67   flu        "            CY125947      1426  fluA-seg6      N11     
5.68   flu        "            CY125948      1027  fluA-seg7      -       
5.69   flu        "            CY125949       895  fluA-seg8      -       
5.70   flu        "            ON637239      1686  fluA-seg4      H19     
#----  ---------  -----------  ----------  ------  -------------  --------
6.1    rsv        "            KY654518     15277  RSV            A       
6.2    rsv        "            MZ516105     15276  RSV            B       
#----  ---------  -----------  ----------  ------  -------------  --------
7.1    mpxv       "            NC_063383   197209  Orthopoxvirus  Monkeypox_virus
```

---

List all information for a particular library/options key:
```
v-scan.pl --l_lib norovirus
```

This will print all the above information (model directory, options,
model information), but only for `norovirus`.

---

List all information for all library/options keys:
```
v-scan.pl --l_all
```

This will print all the above information (model directory, options,
model information), for all libraries.

---

#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.


