# `v-annotate.pl` usage and command-line options

* [`v-annotate.pl` example usage](#exampleusage)
  * classification stage 
  * coverage determination stage
  * alignment and annotation stage
  * protein validation stage
* [`v-annotate.pl` command-line options](#options)

## `v-annotate.pl` example usage<a name="exampleusage"></a>

`v-annotate.pl` uses previously created VADR models from `v-build.pl`
and uses them to annotate sequences in an input sequence file. 
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

1. Classification: each sequence `S` in the input file is compared against
all models in the model library to determine the highest-scoring model
`M(S)` for each sequence. 

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

Go through summary section, 

Then go through some output files



OTHER SECTIONS:
1. -p
2. example showing how to make alert non-fatal




