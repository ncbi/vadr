# <a name="top"></a> `v-scan.pl` example usage and command-line options

* [Quickstart `v-scan.pl` examples](#quickstart)
* [Walk-throughs of `v-scan.pl` examples](#longwalk)
  * [config file](#config)
  * [running `v-scan.pl` when all sequences are expected to match one unknown model library (mode 1)](#mode1)
  * [running `v-scan.pl` when sequences may match multiple models libraries (mode 2)](#mode2)
  * [running `v-scan.pl` when all sequences are expected to match one known model library (mode 3)](#mode3)
* [`v-scan.pl` command-line options](#options)
  * [basic options](#options-basic)
  * [options for specifying which model libraries to use](#options-libraries)
  * [options controlling the random sampling of sequences](#options-sampling)
  * [options for listing information from the config model or about models and exiting(#options-listing)

---

## Quickstart `v-scan.pl` examples <a name="quickexamples"></a>

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

Another key option is the `-c <s>` option to specify a different config
file `<s>` besides the default one. See [here](#config) for more on config files.

## Walk-throughs of `v-scan.pl` examples <a name="longwalk"></a>

This section includes more detailed information on how to use `v-scan.pl`.
`v-scan.pl` is a wrapper script for `v-annotate.pl`. It first
determines which model library to use for the input sequences and then
calls `v-annotate.pl` for that model library. It will supply
`v-annotate.pl` with the command-line options specific for that model
library read from the input config file.

To determine the command-line usage of 
`v-scan.pl` (or any VADR script), use the `-h` option, like this:

```
v-scan.pl -h 
```

You'll see something like the following output:
```
# v-scan.pl :: scan and annotate sequences against VADR model libraries 
# VADR 1.7 (Mar 2025)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Wed Mar  5 10:39:19 2025
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

#### `v-scan.pl` config file<a name="config"></a>

<a name="config"></a> `v-scan.pl` requires a 'config' file that lists information on the
model libraries it will use. Here is the config file that is included
with VADR in [vadr/default.vadr.config](../default.vadr.config) with
comment lines removed for brevity (all lines that begin with a `#` are
comment lines):

```
dengue    $VADRINSTALLDIR/vadr-models-flavi  --split --cpu 1 --group Dengue --nomisc --noprotid --mkey flavi -r
hcv       $VADRINSTALLDIR/vadr-models-flavi  --split --cpu 4 -r --mkey flavi --group HCV
flavi     $VADRINSTALLDIR/vadr-models-flavi  --split --cpu 1 -r --nomisc
norovirus $VADRINSTALLDIR/vadr-models-calici --split --cpu 1 --group Norovirus --nomisc --noprotid --mkey calici -r
calici    $VADRINSTALLDIR/vadr-models-calici --split --cpu 1 -r --nomisc 
```

Note that in the above example config file, both `dengue` and `hcv`
`<options key>` values use the `flavi` model library: `--mkey flavi`
exists in the `<options string>` *and* the `<model directory>` is the
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

Then, when parsing the output for the the scan against the `flavi`
library, sequences are matched to either `dengue`, `hcv` or `flavi` by
checking if the best matching model for each sequence matches to
`dengue` or `hcv` or `flavi`. A model matches to an `<options key>` if
its name equals that `<options key>`, or its `group` or `subgroup`
equals that `<options key>` *after lowercasing and removing all
special characters from the name, group or subgroup. For example, if a
sequence's best matching model is `NC_001477` which has `group` defined
as `Dengue` in the `flavi.minfo` file (relevant line below)
```
MODEL NC_001477 blastdb:"NC_001477.vadr.protein.fa" group:"Dengue" length:"10735" subgroup:"1"
```
then that sequence will match to the `dengue` `<options string>` and
that sequence will then be annotated with `v-annotate.pl` using the
`dengue` `<options string>`. Or, if a sequence matched to a model
named `HCV!` then it would match to `hcv` because `HCV!` becomes `hcv`
after making it lowercase and removing all special (non-alphanumeric)
characters.

Similarly, when parsing the output for the scan against the `calici`
library, sequences are matched to either `norovirus` or `calici` in
the same way. 

You may be wondering why a user wouldn't just separate out all the
`<options key>` models into their own libraries so that each one has
its own unique `<model directory>` and model key. That will certainly
work and it may be preferred by some users, but one reason not to do
that is simply convenience: using a larger library like `flavi` for
`dengue`, `hcv` and other flaviviruses can be more convenient because
it requires less files, and less partitioning of files into separate
model directories. 

#### Adding libraries to the config file<a name="add2config"></a>

You can download or build additional vadr model libraries and add them
to the config file or make your own config file. To use a different
config file `<s>` use the `-c <s>` option. The list of available VADR
models and an example config file that uses them all is 
[here](https://github.com/ncbi/vadr/wiki/Available-VADR-model-files).

#### Running `v-scan.pl` (mode 1)<a name="mode1"></a>

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
# date:             Tue Mar  4 14:16:40 2025
# $VADRCONFIGFILE:  /home/nawrocki/vadr-install-dir/default.vadr.config
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
documentation](annotate.md#lib) also includes a description of several of these files.

The final four files are the [`.log` file](formats.md#log), which is
the same as the standard output printed to the screen currently being
discussed, the [`.cmd` file](formats.md#cmd), and the [`.filelist`
file](formats.md#filelist) which lists the output files created by
`v-scan.pl`, and finally the [`.lib` file](annotate.md#lib)
which explains how many sequences matched to each library in the
classification stage. 

#### Running `v-scan.pl` (mode 2)<a name="mode2"></a>

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

#### Running `v-scan.pl` (mode 3)<a name="mode3"></a>

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
vs-n5-only2`. Alternatively, you can list libraries that you want
`v-scan.pl` to skip (ignore) using the `--skip` option like this:
`v-scan.pl --skip dengue,flavi`.

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

#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.


