# <a name="top"></a> Advanced tutorial: building an RSV model library

---

The `v-build.pl` program will create a model from a single INSDC
accession and include CDS, gene and mature peptide features. However,
often using a model built from a single accession is not general
enough to allow most high quality sequences from a viral species to pass. For
example, some other sequences may include an extended CDS that has a
different stop codon position from the sequence the model was built
from. It is possible to make VADR models more general but it requires
some manual effort. Below I outline some steps I took when building a
general VADR respiratory syncitial virus (RSV) model library.

---

### 1. Determine good reference sequence(s)

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

--- 

### 2. Build initial models from reference sequence(s)

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
$ cp NC_001781/*.vadr.protein.fa* rsv-models/
$ cp NC_038235/*.vadr.protein.fa* rsv-models/

# prepare the library files:
$ $VADRINFERNALDIR/esl-sfetch --index rsv-models/rsv.fa
$ $VADRINFERNALDIR/cmpress rsv-models/rsv.cm
$ $VADRHMMERDIR/hmmpress rsv-models/rsv.hmm
$ $VADRBLASTDIR/makeblastdb -dbtype nucl -in rsv-models/rsv.fa
```

--- 

### 3. Construct a test set for testing your initial models

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

To download 500 randomly chosen RSV sequences of length 14460 or
greater, from the RSV list page you reached in step 1, use the
"Sequence Length" filter to set a minimum of 14460 and then click the
"Download" button, then select "Nucleotide" under "Sequence data
(FASTA format)", then "Download a randomized subset of all records"
then enter "500", click "Next" and "Use Default" for the FASTA
definition line, then finally click "Download". This should download a
file called something like "sequences_20231010_2146812.fasta".
(You can find the accession list for the 500 randomly selected
sequences I used in
[vadr/documentation/build-files/rsv.r500.list](build-files/rsv.r500.list).

---

### 4. Use `v-annotate.pl` to validate and annotate test set sequences and analyze results

The next step is to use our new models to annotate our set of 500
training sequences. The RSV genome is about 15Kb, putting on the
larger end of what VADR can be used to annotate. The memory and speed
requirements for `v-annotate.pl` don't scale very well, and annotation
of RSV can require up to 64Gb of RAM. VADR is used for SARS-CoV-2
(30Kb genome) with significantly lower memory requirements, but it
utilizes heuristics that take advantage of the extremely high (at the time
of writing) sequence identity between all SARS-CoV-2 genomes. I don't
recommend using those heuristics with viral sequences with expected
pairwise sequence similarity of less than 95%. RSV similarity is about
x%.

To use v-annotate.pl on our set of 500 sequences, we would execute:

```
$ v-annotate.pl --mdir rsv-models --mkey rsv sequences_20231010_2146812.fasta va-r500
```

Importantly, this command will take about 1 minute per sequence, so roughly
8 hours to complete. You can parallelize it if you have a lot of RAM
and multiple CPUs using the `--split` and `--cpu` options
as described [here](annotate.md#options-split).

From my set of training sequences 

Next, we want to analyze the results. Ideally we would 

---
### 5. Potentially choose alternative representative sequences

---
### 6. Rerun `v-annotate.pl` on test set (if models were updated)

---
### 7. Iterative refine models and rerun test set

---
## Methods for updating models to improve performance

1. Add new proteins to blastx database

2. Add alternative features to model info file

3. Add alert exceptions 

4. Build new reference alignments and rebuild CM file

5. Use command-line options to `v-annotate.pl` to make some alerts
fatal or not fatal.

---
#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.

