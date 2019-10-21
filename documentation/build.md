# `v-build.pl` example usage and command-line options

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
the files that end with the suffixes `.vadr.cm`,
`.vadr.cm.i1{i,m,f,p}`, `.vadr.protein.fa.p{hr,in,sq}`, and
`.vadr.minfo`. You can either specify that `v-annotate.pl` use only
this model when annotating sequences, by using the `-m`, `-b` and `-i`
options to `v-annotate.pl` as explained in the [`v-annotate.pl`
documentation] (annotate.md), or you can combine these files together
with analogous files from additional `v-build.pl` runs for other
accessions to create a VADR model library, and then use the `-m`, `-b`
and `-i` options to `v-annotate.pl` to specify that library be
used. This is explained in more detail in the [next
section](#general-library). The VADR 1.0 library was created in this
manner, as explained in [another section](#1p0-library).

