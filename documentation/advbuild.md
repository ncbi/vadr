# <a name="top"></a> Advanced tutorial: building an RSV model library

The `v-build.pl` program will create a model from a single INSDC
accession and include CDS, gene and mature peptide features. However,
often using a model built from a single accession is not general
enough to allow most high quality sequences from a viral species to pass. For
example, some other sequences may include an extended CDS that has a
different stop codon position from the sequence the model was built
from.  It is possible to make VADR models more general but it requires
some manual effort. 

A good strategy is to take an iterative approach where you:

1. Build a model(s) from representative and well annotated sequence(s)
   as a starting point. These may be RefSeq sequences.

2. Construct a training set of randomly chosen existing sequences for
   the viral species (or use the same training set from the previous
   iteration) and use `v-annotate.pl` to validate and annotate
   the training sequences using your models from step 1.

3. Analyze the results by looking for common failure modes, and
   investigate the sequence characteristics that are responsible for
   them. This will be characteristics not present in the reference
   sequences. Classify these characteristics into major and minor types: 

   *Major*: these are common to a majority or large fraction of
   sequences but don't exist in the reference sequence. To address
   these, we probably want to pick a new representative sequence(s)
   that includes these characteristics. This means we will return to
   step 1 for another iteration of the three steps.

   *Minor*: these exist in a smaller fraction of sequences, and
   don't exist in the reference sequence. For these characteristics,
   we may be able to update our models without rebuilding them from
   new reference sequences. 

If the initial sequences chosen in the first iteration turn out to
representative enough that there are no major characteristics in step
3, then we can likely stop at the end of iteration 1 and update the
model(s) to tolerate the minor characteristics. 

But if there are one or more major characteristics identified in step
3, we will want to do an additional iteration of all three steps. Two
iterations of all three steps are usually enough to build a nearly
optimal model(s).

This tutorial will walk you through the approach that I took
when building RSV models, using two iterations of the approach
outlined above. The tutorial is long and detailed, and not all of it
may be relevant for what you are interested in. For example, if you
are only interested in fine-tuning existing models to tolerate minor
characteristics of some viral sequences, you may only want to look at
the X.

# RSV model building tutorial, divided into six steps:

<details>

<summary>

## Iteration 1, step 1: build model(s) from initial reference sequence(s)

</summary>

### Determine good reference sequence(s) to use

There are two RSV subtypes, RSV-A and RSV-B, so the first step I took
was to pick reference sequences for each subtype. A good
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

### Build initial models from reference sequence(s)

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
$ mkdir rsv-models1

# concatenate .minfo, .cm .fa and .hmm files:
$ cat NC_038235/*.vadr.minfo > rsv-models1/rsv.minfo
$ cat NC_038235/*.vadr.cm > rsv-models1/rsv.cm
$ cat NC_038235/*.vadr.fa > rsv-models1/rsv.fa
$ cat NC_038235/*.vadr.protein.hmm > rsv-models1/rsv.hmm
$ cat NC_001781/*.vadr.minfo >> rsv-models1/rsv.minfo
$ cat NC_001781/*.vadr.cm >> rsv-models1/rsv.cm
$ cat NC_001781/*.vadr.fa >> rsv-models1/rsv.fa
$ cat NC_001781/*.vadr.protein.hmm >> rsv-models1/rsv.hmm

# copy the blastdb files:
$ cp NC_038235/*.vadr.protein.fa* rsv-models1/
$ cp NC_001781/*.vadr.protein.fa* rsv-models1/

# prepare the library files:
$ $VADRINFERNALDIR/esl-sfetch --index rsv-models1/rsv.fa
$ $VADRINFERNALDIR/cmpress rsv-models1/rsv.cm
$ $VADRHMMERDIR/hmmpress rsv-models1/rsv.hmm
$ $VADRBLASTDIR/makeblastdb -dbtype nucl -in rsv-models1/rsv.fa
```

Before we go on, it's a good idea to test the models by using them to
annotate the reference sequences they were built from. This will
check that we combined the models correctly, and also the sequences
should pass, although they are not guaranteed to do so (some reference
sequences may have special characteristics that cause them to fail
even the sequences they were built from). For these RSV models,
both sequences should pass:

```
$ v-annotate.pl --mdir rsv-models1 --mkey rsv rsv-models1/rsv.fa va-rsv
```

Eventual output:

```
#                                  num   num   num
#idx  model      group  subgroup  seqs  pass  fail
#---  ---------  -----  --------  ----  ----  ----
1     NC_001781  RSV    B            1     1     0
2     NC_038235  RSV    A            1     1     0
#---  ---------  -----  --------  ----  ----  ----
-     *all*      -      -            2     2     0
-     *none*     -      -            0     0     0
#---  ---------  -----  --------  ----  ----  ----
#
# Zero alerts were reported.
#
```

</details>

<details>

<summary>

## Iteration 1, step 2: construct a training set and run `v-annotate.pl` on it

</summary>

### Construct a training set for testing your initial models

INTRODUCE CONCEPT OF IMPROVING/OPTIMIZING MODEL? OR SOME OTHER WORD? INTRODUCE
CONCEPT OF COMMON FAILURE MODES (or something like this).

MAKE A NOTE ABOUT HOW MAKING A LESS BIASED TRAINING SET (x num seqs
per year or something) WOULD BE BETTER.

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

### Run `v-annotate.pl` to validate and annotate sequences in training set

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
$ v-annotate.pl --mdir rsv-models1 --mkey rsv rsv.r500.fasequences_20231010_2146812.fasta va-r500
```

Importantly, this command will take about 1 minute per sequence, so roughly
8 hours to complete. You can parallelize it if you have a lot of RAM
and multiple CPUs using the `--split` and `--cpu` options
as described [here](annotate.md#options-split).

</details>

<details>

<summary>

## Iteration 1, step 3: analyze the results

</summary>

Next, we want to analyze the results. From the v-annotate.pl output
(also saved to the va-r500/va-r500.vadr.log file), we can see that 286
sequences were classified as RSV A (matched best to the NC_038235
model) and 214 as RSV B (matched best to the NC_001781 model). And
only 6 out of the 500 total training sequences "pass":

```
#                                  num   num   num
#idx  model      group  subgroup  seqs  pass  fail
#---  ---------  -----  --------  ----  ----  ----
1     NC_038235  RSV    A          286     2   284
2     NC_001781  RSV    B          214     4   210
#---  ---------  -----  --------  ----  ----  ----
-     *all*      -      -          500     6   494
-     *none*     -      -            0     0     0
#---  ---------  -----  --------  ----  ----  ----
```

The output also contains a summary of reported alerts shows which alerts were most common:

```
# Summary of reported alerts:
#
#     alert     causes   short                               per    num   num  long
#idx  code      failure  description                        type  cases  seqs  description
#---  --------  -------  -----------------------------  --------  -----  ----  -----------
1     fstlocft  no       POSSIBLE_FRAMESHIFT_LOW_CONF    feature      9     9  low confidence possible frameshift in CDS (frame not restored before end)
2     fstlocfi  no       POSSIBLE_FRAMESHIFT_LOW_CONF    feature      3     3  low confidence possible frameshift in CDS (frame restored before end)
3     indf5lcc  no       INDEFINITE_ANNOTATION_START     feature     15    10  alignment to homology model has low confidence at 5' boundary for feature that is or matches a CDS
4     indf3lcc  no       INDEFINITE_ANNOTATION_END       feature     19    12  alignment to homology model has low confidence at 3' boundary for feature that is or matches a CDS
5     insertnn  no       INSERTION_OF_NT                 feature    404   345  too large of an insertion in nucleotide-based alignment of CDS feature
6     lowsim5c  no       LOW_FEATURE_SIMILARITY_START    feature      2     2  region overlapping annotated feature that is or matches a CDS at 5' end of sequence lacks significant similarity
7     lowsim3c  no       LOW_FEATURE_SIMILARITY_END      feature      2     2  region overlapping annotated feature that is or matches a CDS at 3' end of sequence lacks significant similarity
8     lowsimic  no       LOW_FEATURE_SIMILARITY          feature     98    45  region overlapping annotated feature that is or matches a CDS lacks significant similarity
9     ambgnt5f  no       AMBIGUITY_AT_FEATURE_START      feature     17    11  first nucleotide of non-CDS feature is an ambiguous nucleotide
10    ambgnt3f  no       AMBIGUITY_AT_FEATURE_END        feature     30    20  final nucleotide of non-CDS feature is an ambiguous nucleotide
11    ambgnt5c  no       AMBIGUITY_AT_CDS_START          feature     22    14  first nucleotide of CDS is an ambiguous nucleotide
12    ambgnt3c  no       AMBIGUITY_AT_CDS_END            feature     19    12  final nucleotide of CDS is an ambiguous nucleotide
13    ambgcd5c  no       AMBIGUITY_IN_START_CODON        feature      1     1  5' complete CDS starts with canonical nt but includes ambiguous nt in its start codon
14    ambgcd3c  no       AMBIGUITY_IN_STOP_CODON         feature      1     1  3' complete CDS ends with canonical nt but includes ambiguous nt in its stop codon
#---  --------  -------  -----------------------------  --------  -----  ----  -----------
15    lowcovrg  yes      LOW_COVERAGE                   sequence     13    13  low sequence fraction with significant similarity to homology model
16    dupregin  yes      DUPLICATE_REGIONS              sequence    344   344  similarity to a model region occurs more than once
17    deletins  yes      DELETION_OF_FEATURE            sequence      2     1  internal deletion of a complete feature
18    mutstart  yes      MUTATION_AT_START               feature    261   260  expected start codon could not be identified
19    mutendcd  yes      MUTATION_AT_END                 feature      8     8  expected stop codon could not be identified, predicted CDS stop by homology is invalid
20    mutendns  yes      MUTATION_AT_END                 feature      2     2  expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted start codon
21    mutendex  yes      MUTATION_AT_END                 feature      5     5  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position
22    unexleng  yes      UNEXPECTED_LENGTH               feature     19    17  length of complete coding (CDS or mat_peptide) feature is not a multiple of 3
23    cdsstopn  yes      CDS_HAS_STOP_CODON              feature    419   402  in-frame stop codon exists 5' of stop position predicted by homology to reference
24    cdsstopp  yes      CDS_HAS_STOP_CODON              feature    164   164  stop codon in protein-based alignment
25    fsthicft  yes      POSSIBLE_FRAMESHIFT_HIGH_CONF   feature     14    14  high confidence possible frameshift in CDS (frame not restored before end)
26    fsthicfi  yes      POSSIBLE_FRAMESHIFT_HIGH_CONF   feature      1     1  high confidence possible frameshift in CDS (frame restored before end)
27    indfantn  yes      INDEFINITE_ANNOTATION           feature      5     3  nucleotide-based search identifies CDS not identified in protein-based search
28    indf5gap  yes      INDEFINITE_ANNOTATION_START     feature      1     1  alignment to homology model is a gap at 5' boundary
29    indf5lcn  yes      INDEFINITE_ANNOTATION_START     feature     11     8  alignment to homology model has low confidence at 5' boundary for feature that does not match a CDS
30    indf5pst  yes      INDEFINITE_ANNOTATION_START     feature     36    32  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint
31    indf3gap  yes      INDEFINITE_ANNOTATION_END       feature     37    36  alignment to homology model is a gap at 3' boundary
32    indf3lcn  yes      INDEFINITE_ANNOTATION_END       feature   1057   320  alignment to homology model has low confidence at 3' boundary for feature that does not match a CDS
33    indf3pst  yes      INDEFINITE_ANNOTATION_END       feature    218   203  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
34    insertnp  yes      INSERTION_OF_NT                 feature    205   205  too large of an insertion in protein-based alignment
35    lowsim5n  yes      LOW_FEATURE_SIMILARITY_START    feature      1     1  region overlapping annotated feature that does not match a CDS at 5' end of sequence lacks significant similarity
36    lowsim5l  yes      LOW_FEATURE_SIMILARITY_START    feature      2     2  long region overlapping annotated feature that does not match a CDS at 5' end of sequence lacks significant similarity
37    lowsim3n  yes      LOW_FEATURE_SIMILARITY_END      feature      1     1  region overlapping annotated feature that does not match a CDS at 3' end of sequence lacks significant similarity
38    lowsim3l  yes      LOW_FEATURE_SIMILARITY_END      feature      1     1  long region overlapping annotated feature that does not match a CDS at 3' end of sequence lacks significant similarity
39    lowsimin  yes      LOW_FEATURE_SIMILARITY          feature     14    14  region overlapping annotated feature that does not match a CDS lacks significant similarity
40    lowsimil  yes      LOW_FEATURE_SIMILARITY          feature     91    31  long region overlapping annotated feature that does not match a CDS lacks significant similarity
#---  --------  -------  -----------------------------  --------  -----  ----  -----------
```

At this point, we need to determine if these results suggest that our
models should be changed, or if they are giving us the desired
behavior and thus are okay as they are. The fact that so many
sequences fail seems to indicate that the models should be modified,
but it could be that RSV viral sequences are so highly variable that a
high failure rate is expected, and no single sequence based model
would allow significantly more sequences to pass. The only way to know
for sure is to drill down deeper into the results.

To investigate we need to look in further detail at the reasons
sequences are failing. A good way to do this is to go through the most
commonly reported fatal alerts.  Of the 26 fatal alerts (numbers 15 to
40) that occur at least once, the ones that occur in the highest
number of sequences are:

```
23    cdsstopn  yes      CDS_HAS_STOP_CODON              feature    419   402  in-frame stop codon exists 5' of stop position predicted by homology to reference
16    dupregin  yes      DUPLICATE_REGIONS              sequence    344   344  similarity to a model region occurs more than once
32    indf3lcn  yes      INDEFINITE_ANNOTATION_END       feature   1057   320  alignment to homology model has low confidence at 3' boundary for feature that does not match a CDS
18    mutstart  yes      MUTATION_AT_START               feature    261   260  expected start codon could not be identified
34    insertnp  yes      INSERTION_OF_NT                 feature    205   205  too large of an insertion in protein-based alignment
33    indf3pst  yes      INDEFINITE_ANNOTATION_END       feature    218   203  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
```

Information on the individual alert instances can be found in the
`.alt` and `.alt.list` files. We'll go through each of these top six
most common alerts in detail next.

### Investigating common `cdsstopn` alerts

To see all the `cdsstopn` (early stop
codon) alerts in the `.alt` file, we can use `grep`. Here are the
first 10 lines that contain `cdsstopn`:

ADD LINK TO alerts.md cdsstopn example

```
<[(tutorial-20231006)]> grep cdsstopn va-rsv.r500/va-rsv.r500.vadr.alt
1.5.3   OR143220.1  NC_038235  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5615..5617:+      3              5579..5581:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TGA, shifted S:3,M:3]
2.4.1   KX655635.1  NC_038235  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5500..5502:+      3              5579..5581:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TGA, shifted S:3,M:3]
4.5.3   OR287871.1  NC_038235  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5551..5553:+      3              5579..5581:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TGA, shifted S:3,M:3]
5.4.2   OM857265.1  NC_038235  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5544..5546:+      3              5576..5578:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TAA, shifted S:6,M:6]
7.4.1   KJ627366.1  NC_038235  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5498..5500:+      3              5579..5581:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TGA, shifted S:3,M:3]
8.5.2   OR143199.1  NC_038235  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5615..5617:+      3              5579..5581:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TGA, shifted S:3,M:3]
9.2.2   MG431251.1  NC_001781  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5602..5604:+      3              5566..5568:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TAA, shifted S:21,M:21]
10.2.2  KY249668.1  NC_001781  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5611..5613:+      3              5566..5568:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TAA, shifted S:21,M:21]
12.4.1  MK109787.1  NC_038235  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5542..5544:+      3              5579..5581:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TGA, shifted S:3,M:3]
13.2.2  OR326741.1  NC_001781  CDS   attachment_glycoprotein   14  cdsstopn  yes   CDS_HAS_STOP_CODON                      5569..5571:+      3              5566..5568:+      3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TAA, shifted S:21,M:21]
```

The format for the `.alt` file is described
[here](formats.md#alt). The 5th field is the product name, and as you
can see at least the first ten are for the same CDS, the
`attachment_glycoprotein`. The 12th field is the model (reference)
coordinates for the alert, indicating which reference positions in the
model sequence (`NC_038235` or `NC_001781` as indicated in field 3)
the early stop codon aligns/maps to. Note that there are some
positions that occur multiple times, namely `5579..5581:+` for
`NC_038235` and `5566..5568:+` for `NC_001781`. 

We can easily count how many times each model and reference position
pair occurs in the full list using the `grep`, `awk`, `sort` and
`uniq` command line utilities:

```
$ grep cdsstopn va-rsv.r500/va-rsv.r500.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
    213 NC_038235 attachment_glycoprotein 5579..5581:+
    154 NC_001781 attachment_glycoprotein 5566..5568:+
     15 NC_038235 attachment_glycoprotein 5576..5578:+
     15 NC_001781 attachment_glycoprotein 5575..5577:+
      3 NC_038235 phosphoprotein 2459..2461:+
      1 NC_038235 polymerase_protein 9108..9110:+
      1 NC_038235 polymerase_protein 14718..14720:+
      1 NC_038235 polymerase_protein 13347..13349:+
      1 NC_038235 polymerase_protein 13234..13236:+
      1 NC_038235 polymerase_protein 11461..11463:+
      1 NC_038235 polymerase_protein 10096..10098:+
      1 NC_038235 attachment_glycoprotein 5555..5557:+
      1 NC_038235 attachment_glycoprotein 5321..5323:+
      1 NC_038235 M2-2_protein 8175..8177:+
      1 NC_038235 M2-1_protein 8170..8172:+
      1 NC_001781 polymerase_protein 9153..9155:+
      1 NC_001781 polymerase_protein 14986..14988:+
      1 NC_001781 polymerase_protein 14984..14986:+
      1 NC_001781 polymerase_protein 12674..12676:+
      1 NC_001781 polymerase_protein 11015..11017:+
      1 NC_001781 nucleoprotein 2285..2287:+
      1 NC_001781 fusion_glycoprotein 6265..6267:+
      1 NC_001781 attachment_glycoprotein 5554..5556:+
      1 NC_001781 attachment_glycoprotein 5172..5174:+
```

Based on this we can see that 213 of the 286 sequences that match
best to `NC_038235` have an early stop at reference positions
`5579..5581:+`. The reference stop codon position for the `attachment
glycoprotein` CDS can be determined by the stop coordinate for that
CDS in the `rsv.minfo` file:

```
$ grep NC\_038235 rsv-models1/rsv.minfo | grep attachment
FEATURE NC_038235 type:"CDS" coords:"4688..5584:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein"
```
So the most common position for the stop is actually 3 positions
upstream of the `NC_038235` stop codon.  (This can also be inferred
from the detailed alert message in the `.alt` output file above:
`[TGA, shifted S:3,M:3]`.)

Of the 214 RSV B sequences that match best to `NC_001781`, 153 of them
have an early stop at reference positions `5566..5568:+`, which is 21
positions upstream of the attachment glycoprotein CDS stop codon in
`NC_001781`:

```
$ grep NC\_001781 rsv-models1/rsv.minfo | grep attachment
FEATURE NC_001781 type:"CDS" coords:"4690..5589:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein"
```

The fact that more than half of the sequences in our training set have
different stop positions than the reference is a first hint that maybe
our reference models should be updated to a different accession. But
it makes sense at this point to look at all of the common failure
modes before making that decision and potentially looking for new
references. 

### Investigating common `dupregin` alerts

The second most common fatal alert is the `dupregin` alert that occurs
when `similarity to a model region occurs more than once`. 

```
16    dupregin  yes      DUPLICATE_REGIONS              sequence    344   344  similarity to a model region occurs more than once
```

We can use `grep` to look at some examples of the `dupregin` alert:
```
$ grep dupregin va-rsv.r500/va-rsv.r500.vadr.alt 
1.1.1   OR143220.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5502..15225:+,2..5501:+  15224  5466..15199:+,31..5537:+  15241  similarity to a model region occurs more than once [5466..5537:+ (len 72>=20) hits 1 (8569.8 bits) and 2 (4779.2 bits)]
4.1.1   OR287871.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5438..14973:+,1..5437:+  14973  5466..15003:+,93..5537:+  14983  similarity to a model region occurs more than once [5466..5537:+ (len 72>=20) hits 1 (8198.7 bits) and 2 (4724.2 bits)]
5.1.1   OM857265.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5434..14961:+,1..5433:+  14961  5466..14995:+,99..5537:+  14969  similarity to a model region occurs more than once [5466..5537:+ (len 72>=20) hits 1 (8465.0 bits) and 2 (4737.7 bits)]
8.1.1   OR143199.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5502..15222:+,2..5501:+  15221  5466..15199:+,31..5537:+  15241  similarity to a model region occurs more than once [5466..5537:+ (len 72>=20) hits 1 (8590.1 bits) and 2 (4779.3 bits)]
```

The `dupregin` alert is explained more
[here](alerts.md#example-dupregin). The `alert detail` at the end of
each line are very similar for the first four instances of `dupregin`:
`similarity to a model region occurs more than once [5466..5537:+ (len
72>=20)`. This suggests these alerts are all referring to the same
situation. If a genome has repetitive regions you may see this alert
for all sequences, but in that case you will likely see it for the
model sequence as well. If that occurs see SECTION ON EXCEPTION
ALERTS. Even though the alert details are similar note that the model
positions in field 12 are not identical. For example, in the first
alert, the model coordinates are `5466..15199:+,31..5537:+`, and in
the second the coordinates are `5466..15003:+,93..5537:+`. These
coordinates are showing the reference coordinates of the two hits that
overlap, whereas the alert detail shows the actual region of
overlap. So for this alert, in order to group together similar
situations we use `awk` to select different space delimited fields
than in the above `cdsstopn` alert. Here we are more interested in the
24th field (e.g. `[5466..5537:+`), then then 12th. 

Let's again use `grep`, `awk`, `sort` and `uniq` to group the
instances of these alerts together:

```
$ grep dupregin va-rsv.r500/va-rsv.r500.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $23); }' | sort | uniq -c | sort -rnk 1
    167 NC_001781 - [5410..5469:+
    159 NC_038235 - [5466..5537:+
      4 NC_001781 - [5411..5469:+
      3 NC_001781 - [5409..5468:+
      2 NC_038235 - [5466..5536:+
      2 NC_038235 - [5465..5536:+
      2 NC_001781 - [5412..5469:+
      1 NC_038235 - [5478..5518:+
      1 NC_038235 - [5468..5538:+
      1 NC_038235 - [5468..5537:+
      1 NC_038235 - [5467..5537:+
      1 NC_001781 - [5414..5469:+
```

About 56% (159/286) or RSV A and 80% (167/214) of RSV B sequences have
one of two duplicated regions. To get a better idea of this
duplicatedd region, we can select an example of one RSV A sequence and
one RSV B sequence that contain this alert, and run it back through
`v-annotate.pl`. To select two sequences:

```
# prepare the sequence file for the esl-sfetch program 
$ $VADREASELDIR/esl-sfetch --index rsv.r500.fa

# pick a sequence with a dupregin alert that contains the strings we are interested it
$ grep dupregin va-rsv.r500/va-rsv.r500.vadr.alt | grep 5410..5469 | awk '{ print $2 }' | esl-selectn 1 - > ex1.list
$ grep dupregin va-rsv.r500/va-rsv.r500.vadr.alt | grep 5466..5537 | awk '{ print $2 }' | esl-selectn 1 - > ex2.list
$ cat ex1.list 
MZ516003.1
$ cat ex2.list
ON237248.1

# fetch the sequences:
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex1.list > ex1.fa
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex2.list > ex2.fa

# run v-annotate.pl on these sequences with 
# the --out_stk option to save the output alignments
$ v-annotate.pl --out_stk --mdir rsv-models1 --mkey rsv ex1.fa va-ex1
$ v-annotate.pl --out_stk --mdir rsv-models1 --mkey rsv ex2.fa va-ex2
```

Let's look at the RSV B sequence first. The
`va-ex1/va-ex1.vadr.NC_001781.align.stk` file contains the ex1
sequence aligned to the `NC_001781` model in [Stockholm alignment file
format](https://en.wikipedia.org/wiki/Stockholm_format). With
reference positions numbered using rows with the header `#=GC
RFCOL`. Below is an excerpt of that alignment file with the duplicated
region annotated in an extra line that I've added labelled
`dupregin`. The positions marked with `1` show the first instance of
the repeated region, and those marked with `2` show the second
instance. Note that these positions correspond to positions
`5410..5466` of the reference model:

```
MZ516003.1         AATAAACCAAAGAAAAAACCAACTACAAAACCCACAAACAAACCACCTACCAAAACCACAAACAAAAGAGACCCCAAAACACTAGCCAAAACACCGAAAAAAGAAACCACCATTAACCCAACAAAAAAACCAACCCCC
#=GR MZ516003.1 PP ******************************************************************************************************************************************
#=GC SS_cons       ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF            AACAAACCAAAGAAGAAACCAACCATCAAACCCACAAACAAACCAACCACCAAAACCACAAACAAAAGAGACCCAAAAACACCAGCCAAAACGACGAAAAAAGAAACTACCACCAACCCAACAAAAAAACCAACCCTC
#=GC RFCOLX....    000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#=GC RFCOL.X...    555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
#=GC RFCOL..X..    222222222222222222222222222222222222233333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333334
#=GC RFCOL...X.    666666677777777778888888888999999999900000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
#=GC RFCOL....X    345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890

dupregin                    1111111111111111111111111111111111111111111111111111111111112222222222222222222222222222222222222222222222222222222222222
MZ516003.1         AAGACTACAGAAAGAGACACCAGCACCCCACAATCCACTGTGCTCGACATAACCACATcaaaacacacagaaagggacaccagcacctcacaatccattgtgcttgacacaaccgcatCAAAACACACAACCCAACAG
#=GR MZ516003.1 PP ****************99999999999888888888888888888887777777666511111111111111111111100000000000000000000000000000000000000045566777888999******
#=GC SS_cons       ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::............................................................::::::::::::::::::::
#=GC RF            ACGACCACAGAAAGAGACACCAGCACCTCACAATCCACTGTGCTCGACACAACCACAT............................................................TAGAACACACAATCCAACAG
#=GC RFCOLX....    0000000000000000000000000000000000000000000000000000000000............................................................00000000000000000000
#=GC RFCOL.X...    5555555555555555555555555555555555555555555555555555555555............................................................55555555555555555555
#=GC RFCOL..X..    4444444444444444444444444444444444444444444444444444444444............................................................44444444444444444444
#=GC RFCOL...X.    0000000001111111111222222222233333333334444444444555555555............................................................56666666666777777777
#=GC RFCOL....X    1234567890123456789012345678901234567890123456789012345678............................................................90123456789012345678
```

(The lines that include `PP` in the header indicate the alignment
confidence at each position as explained more [here](alerts.md#pp).

This region falls within the `attachment glycoprotein` CDS. To learn
more about this duplication we might search
[PubMed](https://pubmed.ncbi.nlm.nih.gov/) with the query `RSV
attachment glycoprotein duplicated region`, which returns several
articles discussion this region, including a paper entitled
["Functional Analysis of the 60-Nucleotide Duplication in the
Respiratory Syncytial Virus Buenos Aires Strain Attachment
Glycoprotein"](#https://pubmed.ncbi.nlm.nih.gov/26018171/) by Hotard
et al. 

A similar situation exists for the `ex2` sequence `ON237248.1`. Here
is an excerpt of the alignment
`va-ex2/va-ex2.vadr.NC_038235.align.stk`: 

```

                              1111111111111111111111111111111111111111111111111111111111111111111111112222222222222222222222222222222222222222222222222222222222222
ON237248.1         AGAACACACAAGTCAAGAGAAAACCCTCCACTCAACCACCTCCGAAGGCtatctaagcccatcccaagtctatacaacatccggtcaagaggaaaccctccactcaaccacctccgaaggcTATCTAAGCTCATCACAAGTCTA
#=GR ON237248.1 PP ************999887777777776666666666655555555555500000000000000000000000000000000000000000000000000000000000000000000000055555666666667778888888
#=GC SS_cons       :::::::::::::::::::::::::::::::::::::::::::::::::........................................................................:::::::::::::::::::::::
#=GC RF            AGAACTCACAAGTCAAATGGAAACCTTCCACTCAACTTCCTCCGAAGGC........................................................................AATCCAAGCCCTTCTCAAGTCTC
#=GC RFCOLX....    0000000000000000000000000000000000000000000000000........................................................................00000000000000000000000
#=GC RFCOL.X...    5555555555555555555555555555555555555555555555555........................................................................55555555555555555555555
#=GC RFCOL..X..    4444444444444444444444444444444444444444444445555........................................................................55555555555555555555555
#=GC RFCOL...X.    5555566666666667777777777888888888899999999990000........................................................................00000011111111112222222
#=GC RFCOL....X    5678901234567890123456789012345678901234567890123........................................................................45678901234567890123456

                   22222222222
ON237248.1         TACAACATCCGAGTACTTATCACAATCTCTATCTTCATCTAACACAACAAAATGATAGTCATTAAAAAGCGTATTGTTGCAAAAAGCCATGACCAAATCAAGCAGAATCAAAATCAACTCTGGGGCAAATAACAATGGAGTTGC
#=GR ON237248.1 PP 99999999****************************************************************************************************************************************
#=GC SS_cons       ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF            TACAACATCCGAGTACCCATCACAACCTTCATCTCCACCCAACACACCACGCCAGTAGTTACTTAAAAACATATTATCACAAAAAGCCATGACCAACTTAAACAGAATCAAAATAAACTCTGGGGCAAATAACAATGGAGTTGC
#=GC RFCOLX....    000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#=GC RFCOL.X...    555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
#=GC RFCOL..X..    555555555555555555555555555555555555555555555555555555555555555555555555566666666666666666666666666666666666666666666666666666666666666666666666
#=GC RFCOL...X.    222333333333344444444445555555555666666666677777777778888888888999999999900000000001111111111222222222233333333334444444444555555555566666666667
#=GC RFCOL....X    789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890

```

Based on the counts of these two specific instances of the `dupregin` alert
above, we know that many sequences have these exact, or highly
similar, duplications, and that they are not present in the
`NC_001781` and `NC_038235` reference sequences. It is starting to
look like we should change our reference sequences to ones that
include the more common stop position of the attachment glycoprotein,
and that include this duplicated region. 

But before we should choose a new reference there may be more
characteristics that we want to include, so first we should
investigate the other common alerts from above. The third most common
alert was `indf3lcn`:

### Investigating common `indf3lcn` alerts

```
32    indf3lcn  yes      INDEFINITE_ANNOTATION_END       feature   1057   320  alignment to homology model has low confidence at 3' boundary for feature that does not match a CDS
```

This alert occurs because the alignment confidence at the end
coordinate/position  of a feature is too low, indicating that it may
be incorrect, based on the parameters of the model.

We can use `grep` and `awk` again to group together the `indf3lcn`
alerts and see which features they correspond to, this time outputting
the model name, feature type (e.g. CDS or gene) and product name
(fields 3, 4 and 5 in the `.alt` file):

```
$ grep indf3lcn va-rsv.r500/va-rsv.r500.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $4, $5); }' | sort | uniq -c | sort -rnk 1
    265 NC_038235 gene F
    263 NC_038235 gene P
    234 NC_038235 gene N
    226 NC_038235 gene M
     18 NC_001781 gene P
     11 NC_001781 gene L
      8 NC_038235 gene SH
      8 NC_038235 gene G
      7 NC_038235 gene L
      3 NC_038235 gene NS2
      3 NC_038235 gene M2
      3 NC_001781 gene SH
      2 NC_038235 gene NS1
      2 NC_001781 gene F
      1 NC_001781 gene NS2
      1 NC_001781 gene N
      1 NC_001781 gene M2
      1 NC_001781 gene M
```

The vast majority of these alerts correspond to the NC\_038235 model
and all of them are for `gene` features. We can inspect the `gene`
features for that model in the `.minfo` file:

```
$ grep NC_038235 rsv-models1/rsv.minfo
MODEL NC_038235 blastdb:"NC_038235.vadr.protein.fa" cmfile:"NC_038235.vadr.cm" group:"RSV" length:"15222" subgroup:"A"
FEATURE NC_038235 type:"gene" coords:"45..576:+" parent_idx_str:"GBNULL" gene:"NS1"
FEATURE NC_038235 type:"CDS" coords:"99..518:+" parent_idx_str:"GBNULL" gene:"NS1" product:"nonstructural protein 1"
FEATURE NC_038235 type:"gene" coords:"596..1098:+" parent_idx_str:"GBNULL" gene:"NS2"
FEATURE NC_038235 type:"CDS" coords:"628..1002:+" parent_idx_str:"GBNULL" gene:"NS2" product:"nonstructural protein 2"
FEATURE NC_038235 type:"gene" coords:"1125..2327:+" parent_idx_str:"GBNULL" gene:"N"
FEATURE NC_038235 type:"CDS" coords:"1140..2315:+" parent_idx_str:"GBNULL" gene:"N" product:"nucleoprotein"
FEATURE NC_038235 type:"gene" coords:"2329..3242:+" parent_idx_str:"GBNULL" gene:"P"
FEATURE NC_038235 type:"CDS" coords:"2346..3071:+" parent_idx_str:"GBNULL" gene:"P" product:"phosphoprotein"
FEATURE NC_038235 type:"gene" coords:"3252..4209:+" parent_idx_str:"GBNULL" gene:"M"
FEATURE NC_038235 type:"CDS" coords:"3261..4031:+" parent_idx_str:"GBNULL" gene:"M" product:"matrix protein"
FEATURE NC_038235 type:"gene" coords:"4219..4628:+" parent_idx_str:"GBNULL" gene:"SH"
FEATURE NC_038235 type:"CDS" coords:"4303..4497:+" parent_idx_str:"GBNULL" gene:"SH" product:"small hydrophobic protein"
FEATURE NC_038235 type:"gene" coords:"4673..5595:+" parent_idx_str:"GBNULL" gene:"G"
FEATURE NC_038235 type:"CDS" coords:"4688..5584:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein"
FEATURE NC_038235 type:"gene" coords:"5648..7550:+" parent_idx_str:"GBNULL" gene:"F"
FEATURE NC_038235 type:"CDS" coords:"5661..7385:+" parent_idx_str:"GBNULL" gene:"F" product:"fusion glycoprotein"
FEATURE NC_038235 type:"gene" coords:"7597..8557:+" parent_idx_str:"GBNULL" gene:"M2"
FEATURE NC_038235 type:"CDS" coords:"7606..8190:+" parent_idx_str:"GBNULL" gene:"M2" product:"M2-1 protein"
FEATURE NC_038235 type:"CDS" coords:"8159..8431:+" parent_idx_str:"GBNULL" gene:"M2" product:"M2-2 protein"
FEATURE NC_038235 type:"gene" coords:"8489..15067:+" parent_idx_str:"GBNULL" gene:"L"
FEATURE NC_038235 type:"CDS" coords:"8498..14995:+" parent_idx_str:"GBNULL" gene:"L" product:"polymerase protein"
```

Note that the coordinates of the `gene` and `CDS` are different for
the same gene. For example for `F` the `gene` coordinates are
`5648..7550:+`, while the `fusion glycoprotein` CDS has coordinates
`5661..7385:+`. In this case, if we look at the [GenBank
record](https://www.ncbi.nlm.nih.gov/nuccore/NC_038235) we can see
that the gene coordinates correspond to the `mRNA` coordinates for
this gene. This is potentially useful information that could be
annotated on subsequent sequences, but the high number of `indf3lcn`
alerts indicates that VADR is not confident about the annotations of
many of these `gene` boundaries. Because I don't want VADR to get
these annotations wrong, or fail these sequences only because of that
low confidence in the boundaries, I decided to forego the distinct
`gene` boundaries and instead use the `CDS` boundaries as the `gene`
boundaries. The norovirus, dengue, and SARS-CoV-2 VADR models all
enforce this identity between CDS and gene start and end points. If
for your own models you desire the distinct boundaries, then you can
keep them, and possibly set `indf3lcn` alerts as non-fatal using the
`--alt_pass indf3lcn` option to `v-annotate.pl`. 

When we choose our new reference sequences in the next section we will
revisit this issue of differing gene and CDS boundaries.

### Investigating common `mutstart` alerts

Moving on to the fourth most common fatal alert:

```
18    mutstart  yes      MUTATION_AT_START               feature    261   260  expected start codon could not be identified
```

Which model and features does this pertain to?

```
$ grep mutstart va-rsv.r500/va-rsv.r500.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
    259 NC_038235 M2-2_protein 8159..8161:+
      1 NC_038235 phosphoprotein 2365..2367:+
      1 NC_001781 matrix_protein 3263..3265:+
```

All but two of the 261 `mutstart` alert instances pertain to the
M2-2_protein in `NC_038235`. There were only 286 total RSV A
sequences, this means more than 90% of them do not have the start
codon at positions `8159..8161`. To investigate what is going on here,
let's take a random sample of 10 of these 259 sequences, rerun
`v-annotate.pl` on them and look at their alignment to the `NC_038235`
model. 

```
# pick 10 random sequences with the mutstart alert
$ grep mutstart va-rsv.r500/va-rsv.r500.vadr.alt | grep 8159..8161 | awk '{ print $2 }' | esl-selectn 10 - > ex3.list
$ cat ex3.list
MN536997.1
KJ627263.1
MZ515675.1
MF001052.1
MH182035.1
MZ515887.1
MZ516120.1
OR143178.1
OR143202.1
KJ627322.1

# fetch the sequences:
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex3.list > ex3.fa

# run v-annotate.pl on these sequences with 
# the --out_stk option to save the output alignments
$ v-annotate.pl --out_stk --mdir rsv-models1 --mkey rsv ex3.fa va-ex3
```

Below is an excerpt of the resulting alignment in
`va-ex3/va-ex3.vadr_NC_038235.align.stk

```
       
                                                   111   222 
MN536997.1         TAACCCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAAATTCCATA
KJ627263.1         TAACCCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAAATTCCATA
MZ515675.1         TAACCCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAAATTCCATA
MF001052.1         TAACCCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAAATTCCATA
MH182035.1         TAACCCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGATAAATATCCTTGTAGTATAAATTCCATA
MZ515887.1         TAACCCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAAATTCCATA
MZ516120.1         TAACCCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAAATTCCATA
OR143178.1         TAACCCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAAATTCCATA
OR143202.1         TAACNCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAAATTCCATA
KJ627322.1         TAACCCAAAAGAATCAACTGTTAGTGATACGAACGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAAATTCCATA
#=GC SS_cons       ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF            CAACCCAAAAGAATCAACTGTTAGTGATACAAATGACCATGCCAAAAATAATGATACTACCTGACAAATATCCTTGTAGTATAACTTCCATA
#=GC RFCOLX....    00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#=GC RFCOL.X...    88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#=GC RFCOL..X..    11111111111111111111111111111111111111111111111111111111111111111111111112222222222222222222
#=GC RFCOL...X.    22233333333334444444444555555555566666666667777777777888888888899999999990000000000111111111
#=GC RFCOL....X    78901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678
                                                   111   222 
 
```

The `111` at the top indicates the position of the `ATG` for the
`M2-2` protein in the `NC_038325` model (the `RF` line in the
alignment). Note that all 10 sequences have `ACG` aligned at these
positions. The `222` six nucleotides downstream at positions
`8165..8167` indicate another `ATG` which is in-frame with the reference
with positions `8159..8161`. Because the majority of RSV A sequences
in our training set have the second `ATG` but not the first, we
may want our model to use that as the start position. If we do that
however, then the annotated start for the `NC_038325` model will be
annotated differently. There is a way to deal with this and allow
`v-annotate.pl` to pick either start position. We'll come back to this
later. 

### Investigating common `insertnp` alerts

The next most common alert in our training set is:
```
34    insertnp  yes      INSERTION_OF_NT                 feature    205   205  too large of an insertion in protein-based alignment
```

We can determine which model and features these pertain to with:

```
$ grep insertnp va-rsv.r500/va-rsv.r500.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
    147 NC_001781 attachment_glycoprotein 5442..5442:+
     14 NC_038235 attachment_glycoprotein 5509..5509:+
      9 NC_001781 attachment_glycoprotein 5460..5460:+
      7 NC_001781 attachment_glycoprotein 5439..5439:+
      5 NC_001781 attachment_glycoprotein 5445..5445:+
      4 NC_038235 attachment_glycoprotein 5494..5494:+
      3 NC_038235 attachment_glycoprotein 5482..5482:+
      2 NC_038235 attachment_glycoprotein 5524..5524:+
      2 NC_038235 attachment_glycoprotein 5515..5515:+
      2 NC_038235 attachment_glycoprotein 5473..5473:+
      2 NC_001781 attachment_glycoprotein 5448..5448:+
      1 NC_038235 attachment_glycoprotein 5503..5503:+
      1 NC_038235 attachment_glycoprotein 5497..5497:+
      1 NC_001781 attachment_glycoprotein 5472..5472:+
      1 NC_001781 attachment_glycoprotein 5469..5469:+
      1 NC_001781 attachment_glycoprotein 5463..5463:+
      1 NC_001781 attachment_glycoprotein 5457..5457:+
      1 NC_001781 attachment_glycoprotein 5451..5451:+
      1 NC_001781 attachment_glycoprotein 5433..5433:+
```

About 75% of these alerts are for large insertions in the `blastx`
protein alignment in RSV B (`NC_001781`) sequences at position
5442. This position is within the `attachment glycoprotein` region for
which `dupregin` alerts were also reported: 5410..5469. It makes sense
that a large duplication could result in a large insertion. In fact,
it is somewhat surprising that there aren't similar alerts for the
`NC_038235` model, although we will find out why soon enough. Dealing
with the `NC_001781` duplication in the next round of model building
will likely eliminate these alerts.

It is possible to look at these insertions in the actual blastx output
alignments, but you'll need to rerun `v-annotate.pl` using the
`--keep` option, which makes it save all intermediate files that are
usually deleted. To do that we would sample one sequences and
rerun them (you can sample more than one, but often it is easier to
find the relevant part of the blastx output if you have restricted
your input to a single sequence):

```
# pick a sequence with the insertnp alert:
$ grep insertnp va-rsv.r500/va-rsv.r500.vadr.alt | grep 5442 | awk '{ print $2 }' | esl-selectn 1 - > ex4.list
$ cat ex4.list
MH760706.1

# fetch the sequence
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex4.list > ex4.fa

# run v-annotate.pl on these sequences with 
# the --keep option to save all output files
$ v-annotate.pl --keep --mdir rsv-models1 --mkey rsv ex4.fa va-ex4
```

The relevant `blastx` output will be in the file
`va-ex4/va-ex4.vadr.NC_001781.blastx.out`. The section of the file you
are looking for is the results for when the predicted `attachment glycoprotein`
CDS is used as a query sequence. To find this we need to know what the
coordinates are for that prediction. We can find these in the `.ftr`
output file (format described [here](formats.md#ftr) or in the `.tbl`
file `va-ex4/va-ex4.vadr.fail.tbl`: 

```
4590	5543	misc_feature
			note	similar to attachment glycoprotein
```

If we search for `4590..5543` in the `blastx` output file, we will
find the results and alignments:

```
Query= MH760706.1/CDS.7/4590..5543:+

Length=954
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

NC_001781.1/4690..5589:+                                              363     6e-130
NC_001781.1/3263..4033:+                                              21.2    0.51  
NC_001781.1/8509..15009:+                                             20.4    0.97  
NC_001781.1/5666..7390:+                                              20.0    1.1   
NC_001781.1/99..518:+                                                 16.9    8.9   


>NC_001781.1/4690..5589:+
Length=299

 Score = 363 bits (932),  Expect = 6e-130, Method: Compositional matrix adjust.
 Identities = 265/319 (83%), Positives = 273/319 (86%), Gaps = 22/319 (7%)
 Frame = +1

Query  1    MSKNKNQRTARTLEKTWDTLNHLIVISSCLYKLNLKSIAQIALSVLAMIISTSLIIAAII  180
            MSK+KNQRTARTLEKTWDTLNHLIVISSCLY+LNLKSIAQIALSVLAMIISTSLIIAAII
Sbjct  1    MSKHKNQRTARTLEKTWDTLNHLIVISSCLYRLNLKSIAQIALSVLAMIISTSLIIAAII  60

Query  181  FIISANHKVTLTTVTVQTIKNHTEKNMTTYLTQVSPERVSPSKQPTATPPIHTNSATISP  360
            FIISANHKVTLTTVTVQTIKNHTEKN+TTYLTQV PERVS SKQPT T PIHTNSAT SP
Sbjct  61   FIISANHKVTLTTVTVQTIKNHTEKNITTYLTQVPPERVSSSKQPTTTSPIHTNSATTSP  120

Query  361  NTKSETHHTTAQTKGTTSTPTQNNKPSTEPRPKKPPK--KDDYHFEVFNFVPCSICGNNQ  534
            NTKSETHHTTAQTKG T+T TQ NKPST+PR K PPK  KDDYHFEVFNFVPCSICGNNQ
Sbjct  121  NTKSETHHTTAQTKGRTTTSTQTNKPSTKPRLKNPPKKPKDDYHFEVFNFVPCSICGNNQ  180

Query  535  LCKSICKTIPSNKPKKKPTTKPTNKPPTKTTNKRDPKTLAKTPKKENTINPTKKPTPKTT  714
            LCKSICKTIPSNKPKKKPT KPTNKP TKTTNKRDPKT AKT KKE T NPTKKPT  TT
Sbjct  181  LCKSICKTIPSNKPKKKPTIKPTNKPTTKTTNKRDPKTPAKTTKKETTTNPTKKPTLTTT  240

Query  715  ERDTSTPQSTVLDITTSKHTERDTSTSQSIALDTTTSKHTTQQQSLYSTTPENTPNSTQT  894
            ERDTST QSTV                    LDTTT +HT QQQSL+STTPENTPNSTQT
Sbjct  241  ERDTSTSQSTV--------------------LDTTTLEHTIQQQSLHSTTPENTPNSTQT  280

Query  895  PTASEPSTSNST*RLQSYA  951
            PTASEPSTSNST   QS+A
Sbjct  281  PTASEPSTSNSTQNTQSHA  299

```

Note the large insertion in the query around position 748 of the CDS.

---

### Investigating common `indf3pst` alerts

The sixth and final common alert that we will look at is:

```
33    indf3pst  yes      INDEFINITE_ANNOTATION_END       feature    217   202  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
```

To learn more about the model and features for this alert:

```
$ grep indf3pst va-rsv.r500/va-rsv.r500.vadr.alt | awk '{ printf ("%s %s\n", $3, $5); }' | sort | uniq -c | sort -rnk 1
    182 NC_038235 attachment_glycoprotein
     12 NC_001781 attachment_glycoprotein
      9 NC_038235 polymerase_protein
      7 NC_001781 polymerase_protein
      2 NC_038235 fusion_glycoprotein
      2 NC_038235 M2-2_protein
      1 NC_038235 phosphoprotein
      1 NC_038235 M2-1_protein
      1 NC_001781 small_hydrophobic_protein
      1 NC_001781 nucleoprotein
```

The `indf3pst` alert occurs when the `blastx` alignment in the protein
validation stage does not extend closely enough to the 3' end. As this
is again in the attachment glycoprotein CDS it may be related to the
duplicated region. Note that the vast majority of instances are for
RSV A (`NC_038235`) sequences. A possible explanation is that these
are predominantly made up of sequences that also have a `dupregin`
alert and `blastx` does not create an alignment that spans the
duplicated region as it did for the example RSV B sequence MH760706.1
with the `insertnp` alert above. The duplicated region was typically the same
size in both RSV A and RSV B sequences, so the difference between VADR
reporting an `insertnp` or `indf3pst` alert is probably due to the 
similarity in the 3' ends of the protein between the training
sequences and the references: for RSV B, the 3' ends are similar
enough that the best scoring `blastx` alignment extends across the
duplication, whereas for RSV A, the 3' end is not similar enough and
the `blastx` alignment stops before the duplication. To verify this we
can again look at the `blastx` output after rerunning an example
sequence:

```
# pick a sequence with the indf3pst alert:
$ grep indf3pst va-rsv.r500/va-rsv.r500.vadr.alt | grep attachment | grep NC_038235 | awk '{ print $2 }' | esl-selectn 1 - > ex5.list
$ cat ex5.list
MH181932.1

# fetch the sequence
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex5.list > ex5.fa

# run v-annotate.pl on these sequences with 
# the --keep option to save all output files
$ v-annotate.pl --keep --mdir rsv-models1 --mkey rsv ex5.fa va-ex5
```

Again, consult the `.tbl` file to determine attachment glycoprotein
coordinates: 
```
4637	5605	misc_feature
			note	similar to attachment glycoprotein
```

And in the `va-ex5/va-ex5.vadr.NC_038235.blastx.out` file:

```
Query= MH181932.1/CDS.7/4637..5605:+

Length=969
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

NC_038235.1/4688..5584:+                                              362     2e-129
NC_038235.1/1140..2315:+                                              24.6    0.045 
NC_038235.1/7606..8190:+                                              20.4    0.72  
NC_038235.1/8498..14995:+                                             20.0    1.3   
NC_038235.1/2346..3071:+                                              19.6    1.4   


>NC_038235.1/4688..5584:+
Length=298

 Score = 362 bits (929),  Expect = 2e-129, Method: Compositional matrix adjust.
 Identities = 251/283 (89%), Positives = 257/283 (91%), Gaps = 0/283 (0%)
 Frame = +1

Query  1    MSKTKDQRTAKTLERTWDTLNHLLFISSCLYKLNLKSIAQITLSILAMIISTSLIIAAII  180
            MSK KDQRTAKTLERTWDTLNHLLFISSCLYKLNLKS+AQITLSILAMIISTSLIIAAII
Sbjct  1    MSKNKDQRTAKTLERTWDTLNHLLFISSCLYKLNLKSVAQITLSILAMIISTSLIIAAII  60

Query  181  FIASANHKVTLTTAIIQDATNQIKNTTPTYLTQNPQLGISFSNLSGTTSQSTTILASTTP  360
            FIASANHKVT TTAIIQDAT+QIKNTTPTYLTQNPQLGIS SN S  TSQ TTILASTTP
Sbjct  61   FIASANHKVTPTTAIIQDATSQIKNTTPTYLTQNPQLGISPSNPSEITSQITTILASTTP  120

Query  361  SAESTPQSTTVKIKNITTTQILPSKPTTKQRQNKPQNKPNNDFHFEVFNFVPCSICSNNP  540
              +ST QSTTVK KN TTTQ  PSKPTTKQRQNKP +KPNNDFHFEVFNFVPCSICSNNP
Sbjct  121  GVKSTLQSTTVKTKNTTTTQTQPSKPTTKQRQNKPPSKPNNDFHFEVFNFVPCSICSNNP  180

Query  541  TCWAICKRIPNKKPGKKTTTKPTKKPTLKTTKKDPKPQTTKPKEVLTTKPTGKPTINTTK  720
            TCWAICKRIPNKKPGKKTTTKPTKKPTLKTTKKDPKPQTTK KEV TTKPT +PTINTTK
Sbjct  181  TCWAICKRIPNKKPGKKTTTKPTKKPTLKTTKKDPKPQTTKSKEVPTTKPTEEPTINTTK  240

Query  721  TNIRTILLTSNTKGNPEHTSQEETLHSTTSEGYPSPSQVYTTS  849
            TNI T LLTSNT GNPE TSQ ET HST+SEG PSPSQV TTS
Sbjct  241  TNIITTLLTSNTTGNPELTSQMETFHSTSSEGNPSPSQVSTTS  283
```

Note that there is no insertion in this alignment. The alignment ends
just after the first occurence of the duplicated region in
`MH181932.1`. If it extended all the way to the stop codon of the
prediction region `4688..5584` the end `Query` position would be 969
because `5584-4688+1=969`. This is why the length of the relevant sequence region reported in
the `alert detail` in the `alt` file is `120`:

```
<[(va-ex5)]> grep indf3pst va-ex5.vadr.alt
1.5.3  MH181932.1  NC_038235  CDS   attachment_glycoprotein   14  indf3pst  yes   INDEFINITE_ANNOTATION_END             5486..5605:+    120              5584..5584:+      1  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint [120>8, valid stop codon in nucleotide-based prediction]
```
--- 

### <a name="majorchar"></a>Lessons from investigating common alerts in our training set

We've learned that both `NC_038235` and `NC_001781` are lacking
several major characteristics that are (by definition) present in the
majority of RSV A and RSV B sequences, respectively. These include:

 * attachment glycoprotein stop codon at reference coordinates
   `5579..5581` in RSV A and `5566..5568` in RSV B.

 * duplicated region in attachment glycoprotein near reference
   reference positions `5466..5537` in RSV A and `5410..5469` in RSV
   B. 

 * M2-2 protein start codon at reference positions `8165..8167`

And further we've observed that:

  * `NC_038235` and `NC_001781` both have `gene` positional
    boundaries that differ from the CDS. Some of the `gene`
    boundaries lead to low alignment confidence related alerts in
    VADR. Typically for viral GenBank submissions based on VADR, the
    gene and CDS boundaries are kept consistent.

At this point, because we've found at least one major characteristic
present in the majority of the sequences but lacking in the reference
model for both models, we will start a new iteration of our model
building strategy, by first choosing representative sequences that
contain these major characteristics and building new models from them.

</details>

---

<details>

<summary>

## Iteration 2, step 1: build model(s) from new reference sequence(s)

</summary>

## Choosing new representative sequences

We will choose a new representative from our random set of 500
training sequences. Of course, it is possible, and probably even
likely, that there is a 'better' representative that is not in our
training set, but using a sequence from our training set of 500 is
probably good enough. First we need to identify the subset of our 500
sequences that contain the three major characteristics we identified
in iteration 1, summarized [above](#majorchar). We'll do this
separately for RSV A and RSV B:

```
# fetch out the sequences with the early stop codon (cdsstopn) at 5579..5581:
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep cdsstopn | grep 5579..5581 | awk '{ print $2 }' | wc -l
213
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep cdsstopn | grep 5579..5581 | awk '{ print $2 }' | sort > 213.list 

# fetch out the sequences with the dupregin at 5466..5537:
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep dupregin | grep 5466..5537 | awk '{ print $2 }' | wc -l
159
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep dupregin | grep 5466..5537 | awk '{ print $2 }' | sort > 158.list

# use the unix 'comm' command to get the subset of seqs common to both lists
$ comm -1 -2 213.list 159.list | wc -l
143
$ comm -1 -2 213.list 159.list > 143.list

# fetch out the sequences with the M2-2 protein start codon at reference positions `8165..8167
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep M2-2 | grep mutstart | grep 8159..8161 | awk '{ print $2 }' | wc -l
259
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep M2-2 | grep mutstart | grep 8159..8161 | awk '{ print $2 }' | sort > 259.list

# use the unix 'comm' command to get the subset of seqs common to both lists
$ comm -1 -2 143.list 259.list | wc -l
140
$ comm -1 -2 143.list 259.list > 140.list

# fetch the 139 sequences into a new fasta file:
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa 139.list > rsvA.139.fa
```
And then repeat the same for RSV B, skipping the M2-2 start codon step
which doesn't apply to RSV B:
```
# fetch out the sequences with the early stop codon (cdsstopn) at 5566..5568:
$ grep NC_001781 va-rsv.r500/va-rsv.r500.vadr.alt | grep cdsstopn | grep 5566..5568 | awk '{ print $2 }' | wc -l
154
$ grep NC_001781 va-rsv.r500/va-rsv.r500.vadr.alt | grep cdsstopn | grep 5566..5568 | awk '{ print $2 }' | sort > 153.list

# fetch out the sequences with the dupregin at 5410..5469
$ grep NC_001781 va-rsv.r500/va-rsv.r500.vadr.alt | grep dupregin | grep 5410..5469 | awk '{ print $2 }' | wc -l
167
$ grep NC_001781 va-rsv.r500/va-rsv.r500.vadr.alt | grep dupregin | grep 5410..5469 | awk '{ print $2 }' | sort > 167.list

# use the unix 'comm' command to get the subset of seqs common to both lists
$ comm -1 -2 153.list 167.list | wc -l
129
$ comm -1 -2 153.list 167.list > 129.list

# fetch the 129 sequences into a new fasta file:
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa 129.list > rsvB.129.fa
```

It is important that our reference sequences do not have too many
ambiguous nucleotides as make the models less specific. (Indeed we may
want our models to be less specific and more general in some areas of
the genome that are less highly conserved, but including ambiguous
nucleotides from a single reference sequence is not the best way to do
this. We'll revisit this topic briefly at the end of the tutorial.) 

Because we have many candidates, we can afford to remove sequences
with 1 or more ambiguous nucleotides. We can list all such sequences
using the `count-ambigs.pl` *miniscript* included with VADR:

```
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl rsvA.140.fa
KJ672441.1     0 15173 0.0000
KJ672451.1     0 15173 0.0000
KJ672457.1     0 15172 0.0000
KU839631.1     0 15172 0.0000
KU950492.1     0 15233 0.0000
KU950506.1     0 15202 0.0000
KU950524.1     0 15233 0.0000
KU950596.1     0 15228 0.0000
KU950627.1     0 15228 0.0000
KU950639.1     0 15232 0.0000
KU950666.1     0 15061 0.0000
..snip..
```

The second field is the number of ambiguous nucleotides in each
sequence. We can fetch out all lines that have 1 or more in this field
with: 

```
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl rsvA.140.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep -v " 0" 
LC530050.1 1
MG813982.1 46
MN078121.1 7
MN535098.1 217
MN536995.1 1238
MN536996.1 65
..snip..
```

And to only save the sequences with 0 ambiguous nucleotides:
```
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl rsvA.140.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | wc -l 
90
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl rsvA.140.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" > rsvA.90.list
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa rsvA.90.list > rsvA.90.fa
```

Repeating for RSV B:
```
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl rsvB.129.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | wc -l
113
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl rsvB.129.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" > rsvB.113.list
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa rsvB.113.list > rsvB.113.fa
```

We will choose our reference sequences from these candidate sets using
two criteria: 
* similarity to other candidate sequences
* length

We want a sequence that is representative and one way to do that is to
pick the sequence with a high average percent identity to all other
candidates based on an alignment. We can use `v-annotate.pl` to
generate multiple alignments of all candidates:

```
$ v-annotate.pl --out_stk --mdir rsv-models1 --mkey rsv rsvA.89.fa va-rsvA.89
$ v-annotate.pl --out_stk --mdir rsv-models1 --mkey rsv rsvB.112.fa va-rsvB.112
```

When these commands are completed the alignment will be in 
`va-rsvA.89/va-rsv.89.vadr.NC_038235.align.stk` and
`va-rsvB.112/va-rsv.112.vadr.NC_001781.align.stk`.

We can use the `esl-alipid` program that is installed with VADR to
determine the alignment percent identity between all pairs of
sequences. The `esl-alipid-per-seq-stats.pl` script can then be
used to output the average pairwise identities for each sequence,
which we can use to select a new representative:

```
# compute the pairwise ids with esl-alipid:
$ $VADREASELDIR/esl-alipid va-rsvA.90/va-rsvA.90.vadr.NC_038235.align.stk > rsvA.90.alipid
$ head rsvA.90.alipid 
# seqname1 seqname2 %id nid denomid %match nmatch denommatch
KJ672441.1 KJ672451.1  99.29  15066  15173  99.03  15099  15247
KJ672441.1 KJ672457.1  99.12  15038  15172  99.04  15099  15246
KJ672441.1 KU839631.1  99.19  15049  15172  99.05  15100  15245
KJ672441.1 KU950492.1  99.27  15062  15173  98.65  15100  15306
KJ672441.1 KU950506.1  99.03  15026  15173  98.74  15091  15284
KJ672441.1 KU950524.1  99.28  15064  15173  98.64  15099  15307
KJ672441.1 KU950596.1  99.22  15055  15173  98.66  15098  15303
KJ672441.1 KU950627.1  99.24  15057  15173  98.67  15099  15302
KJ672441.1 KU950639.1  99.19  15050  15173  98.66  15100  15305

# compute the per-seq average percent id 
$ perl $VADRSCRIPTSDIR/miniscripts/esl-alipid-per-seq-stats.pl rsvA.90.alipid > rsvA.90.alipid.perseq
$ head rsvA.90.alipid.perseq 
#seq        avgpid  minpidseq   minpid  maxpidseq   maxpid
KJ672441.1  98.818  OQ941773.1  95.040  KJ672451.1  99.290
KJ672451.1  98.860  OQ941773.1  94.920  KU950680.1  99.840
KJ672457.1  98.794  OQ941773.1  95.380  OR466347.1  99.680
KU839631.1  98.824  OQ941773.1  95.420  KU950639.1  99.950
KU950492.1  98.853  OQ941773.1  94.950  KY982516.1  99.700
KU950506.1  98.771  OQ941773.1  95.370  KY654518.1  99.630
KU950524.1  98.890  OQ941773.1  94.910  KU950627.1  99.660
KU950596.1  98.823  OQ941773.1  94.860  KU950627.1  99.990
KU950627.1  98.827  OQ941773.1  94.860  KU950596.1  99.990

# sort by average percent id:
$ grep -v ^\# rsvA.90.alipid.perseq | sort -rnk 2 | head > rsvA.top10.alipid.perseq
$ cat rsvA.top10.alipid.perseq 
KY654518.1  98.950  OQ941773.1  95.480  KU950639.1  99.670
KY982516.1  98.922  OQ941773.1  95.000  KU950492.1  99.700
KU950524.1  98.890  OQ941773.1  94.910  KU950627.1  99.660
KX655644.1  98.871  OQ941773.1  94.930  MH181953.1  99.780
KU950680.1  98.869  OQ941773.1  94.910  KJ672451.1  99.840
KJ672451.1  98.860  OQ941773.1  94.920  KU950680.1  99.840
KU950639.1  98.856  OQ941773.1  95.430  KU839631.1  99.950
KU950492.1  98.853  OQ941773.1  94.950  KY982516.1  99.700
KU950627.1  98.827  OQ941773.1  94.860  KU950596.1  99.990
KU839631.1  98.824  OQ941773.1  95.420  KU950639.1  99.950
```

And for RSV B:

```
# compute the pairwise ids with esl-alipid:
$ $VADREASELDIR/esl-alipid va-rsvB.113/va-rsvB.113.vadr.NC_001781.align.stk > rsvB.113.alipid

# compute the per-seq average percent id 
$ perl $VADRSCRIPTSDIR/miniscripts/esl-alipid-per-seq-stats.pl rsvB.113.alipid > rsvB.113.alipid.perseq

# sort by average percent id:
$ grep -v ^\# rsvB.113.alipid.perseq | sort -rnk 2 | head > rsvB.top10.alipid.perseq
$ cat rsvB.top10.alipid.perseq
ON237084.1  98.927  ON237101.1  97.950  MH760654.1  99.800
ON237174.1  98.902  ON237101.1  97.870  ON237173.1  99.930
MH760699.1  98.887  ON237101.1  97.860  ON237174.1  99.680
MH760654.1  98.875  ON237101.1  97.900  ON237084.1  99.800
ON237177.1  98.873  ON237101.1  97.850  ON237174.1  99.930
ON237173.1  98.873  ON237101.1  97.840  ON237174.1  99.930
OR496332.1  98.857  ON237101.1  97.770  MZ516105.1  99.760
MZ516105.1  98.851  ON237101.1  97.760  OR496332.1  99.760
MW160818.1  98.849  ON237101.1  97.740  OR496332.1  99.660
MH760707.1  98.837  ON237101.1  97.860  MH760706.1  99.850
```

The sequences with the highest average percent identities are good
candidates. The final criterion is the length, we don't want the
sequence to be too short. In this case, we know that the
NC_038235 and NC_001781 models are considered "full length" as
indicated in their GenBank annotation:

[NC_038235](https://www.ncbi.nlm.nih.gov/nuccore/NC_038235)
```
COMMENT     REVIEWED REFSEQ: This record has been curated by NCBI staff. The
            reference sequence is identical to M74568.
            COMPLETENESS: full length.
```

[NC_001781](https://www.ncbi.nlm.nih.gov/nuccore/NC_001781) RefSeq
```
COMMENT     REVIEWED REFSEQ: This record has been curated by NCBI staff. The
            reference sequence is identical to AF013254.
            COMPLETENESS: full length.

```

So our new representatives should align end to end with their
respective model RefSeqs. We can determine those that do by inspecting
the alignments. To make viewing the alignments easier, we can pull out
only the top 10 candidates from each using the `esl-alimanip` program
that is installed with VADR:

```
$ cat rsvA.top10.alipid.perseq | awk '{ print $1 }' > rsvA.top10.list
$ $VADREASELDIR/esl-alimanip --seq-k rsvA.top10.list va-rsvA.90/va-rsvA.90.vadr.NC_038235.align.stk > rsvA.top10.stk

$ cat rsvB.top10.alipid.perseq | awk '{ print $1 }' > rsvB.top10.list
$ $VADREASELDIR/esl-alimanip --seq-k rsvB.top10.list va-rsvB.113/va-rsvB.113.vadr.NC_001781.align.stk > rsvB.top10.stk
```

Sequences that align to the full length of the reference model will
have zero terminal gaps at reference (nongap positions in the `#=GC
RF` lines). Taking a look at the 5' end of the `rsvA.top10.stk` alignment:

```
KJ672451.1         ----------------------------------------------------------------------CACTTAAATTTAACTCCT
#=GR KJ672451.1 PP ......................................................................******************
KU839631.1         ----------------------------------------------------------------------CACTTAAATTTAACTCCT
#=GR KU839631.1 PP ......................................................................******************
KU950492.1         ---------------------AAACTTGCGTAAACCAAAAAAATGGGGCAAATAAGAATTTGATAAGTACCACTTAAATCTAACTCCT
#=GR KU950492.1 PP .....................*******************************************************************
KU950524.1         ---------------------AAACTTGCGTAAACCAAAAAAATGGGGCAAATAAGAATTTGATAAGTACCACTTAAATTTAACTCCT
#=GR KU950524.1 PP .....................*******************************************************************
KU950627.1         --------------------------TGCGTAAACCAAAAAAATGGGGCAAATAAGAATTTGATAAGTACCACTTAAATTTAACTCCT
#=GR KU950627.1 PP ..........................**************************************************************
KU950639.1         ---------------------AAACTTGCGTAAACCAAAAAAATGGGGCAAATAAGAATTTGATAAGTACCACTTAAATTTAACTCCT
#=GR KU950639.1 PP .....................*******************************************************************
KU950680.1         -----------------------------GTAAACCAAAAAAATGGGGCAAATAAGAATTTGATAAGTACCACTTAAATTTAACTCCT
#=GR KU950680.1 PP .............................***********************************************************
KX655644.1         -----------------------------------CAAAAAAATGGGGCAAATAAGAATTTGATAAGTACCACTTAAATTTAACTCCT
#=GR KX655644.1 PP ...................................*****************************************************
KY654518.1         ACGCGAAAAAATGCGTACAACAAACTTGCGTAAACCAAAAAAATGGGGCAAATAAGAATTTGATAAGTACCACTTAAATTTAACTCCT
#=GR KY654518.1 PP ****************************************************************************************
KY982516.1         -----------------------------------CAAAAAAATGGGGCAAATAAGAATTTGATAAGTACCACTTAAATTTAACTCCT
#=GR KY982516.1 PP ...................................*****************************************************
#=GC SS_cons       ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF            ACGCGAAAAAATGCGTACAACAAACTTGCATAAACCAAAAAAATGGGGCAAATAAGAATTTGATAAGTACCACTTAAATTTAACTCCC
```

Only one sequence `KY654518.1` extends to the 5' end of the reference
model. Fortunately, it also extends to the 3' end as well:

```
KJ672451.1         TTATATGTATATTAACTAAATT-----------------------------------
#=GR KJ672451.1 PP **********************...................................
KU839631.1         TTATATGTATATTAACTAAATT-----------------------------------
#=GR KU839631.1 PP **********************...................................
KU950492.1         TTATATGTATATTAACTAAATTACGAGATATTA------------------------
#=GR KU950492.1 PP *********************************........................
KU950524.1         TTATATGTATATTAACTAAATTACGAGATATTA------------------------
#=GR KU950524.1 PP *********************************........................
KU950627.1         TTATATGTATATTAACTAAATTACGAGATATTA------------------------
#=GR KU950627.1 PP *********************************........................
KU950639.1         TTATATGTATATTAACTAAATTACGAGATATTA------------------------
#=GR KU950639.1 PP *********************************........................
KU950680.1         TTATATGTATATTAACTAAATTACGAGATATTA------------------------
#=GR KU950680.1 PP *********************************........................
KX655644.1         TTATATGTATATTAACTAAATTACGAGATATTA------------------------
#=GR KX655644.1 PP *********************************........................
KY654518.1         TTATATGTATATTAACTAAATTACGAGATATTAGTTTTTGACACTTTT-TTTCTCGT
#=GR KY654518.1 PP ************************************************.********
KY982516.1         TTATATGTATATTAACTAAATTACGAGATATTA------------------------
#=GR KY982516.1 PP *********************************........................
#=GC SS_cons       ::::::::::::::::::::::::::::::::::::::::::::::::.::::::::
#=GC RF            TTATATGTGTATTAACTAAATTACGAGATATTAGTTTTTGACACTTTT.TTTCTCGT
//
```

So `KY654518.1` will be our new RSV A reference sequence. Repeat the
drill for the `rsvB.top10.stk` alignment:

```
MH760654.1         ----------------------------------------------------------------------------------------
#=GR MH760654.1 PP ........................................................................................
MH760699.1         ----------------------------------------------------------------------------------------
#=GR MH760699.1 PP ........................................................................................
MH760707.1         ----------------------------------------------------------------------------------------
#=GR MH760707.1 PP ........................................................................................
MW160818.1         -------------------------TTGCATACTCGAAAA-AAATGGGGCAAATAAGAATTTGATAAGTGCTATTTAAGTCTAACCTT
#=GR MW160818.1 PP .........................***************.***********************************************
MZ516105.1         ACGCGAAAAAATGCGTACTACAAACTTGCACACTCGGAAA-AAATGGGGCAAATAAGAATTTGATGAGTGCTATTTAAGTCTAACCTT
#=GR MZ516105.1 PP ****************************************.***********************************************
ON237084.1         ---------------------------------------------GGGGCAAATAAGAATTTGATAAGTGCTATTTAAGTCTAACCTT
#=GR ON237084.1 PP .............................................*******************************************
ON237173.1         --------------------------------------------TGGGGCAAATAAGAATTTGATAAGTGCTATTTAAGTCTAACCTT
#=GR ON237173.1 PP ............................................********************************************
ON237174.1         --------------------------------------------TGGGGCAAATAAGAATTTGATAAGTGCTATTTAAGTCTAACCTT
#=GR ON237174.1 PP ............................................********************************************
ON237177.1         ------------------------------------------AATGGGGCAAATAAGAATTTGATAAGTGCTATTTAAGTCTAACCTT
#=GR ON237177.1 PP ..........................................**********************************************
OR496332.1         ----------------------------------------------GGGCAAATAAGAATTTGATGAGTGCTATTTAAGTCTAACCTT
#=GR OR496332.1 PP ..............................................******************************************
#=GC SS_cons       ::::::::::::::::::::::::::::::::::::::::.:::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF            ACGCGAAAAAATGCGTACTACAAACTTGCACATTCGGAAA.AAATGGGGCAAATAAGAATTTGATAAGTGCTATTTAAGTCTAACCTT
```

Again, only 1 candidate `MZ516105.1`, which fortunately also extends
to the 3' end position as well. (Note that there is one gap at the end
of `MZ516105.1` but this is to a gap in the `RF` annotation. The
sequence does extend to the final nongap position of the `RF`
annotation.)

```
MH760654.1         -----------------------------------------------------------------------------
#=GR MH760654.1 PP .............................................................................
MH760699.1         -----------------------------------------------------------------------------
#=GR MH760699.1 PP .............................................................................
MH760707.1         -----------------------------------------------------------------------------
#=GR MH760707.1 PP .............................................................................
MW160818.1         GTCTAAAACTAACAATCACACATGTGCATTTGCAACACA--------------------------------------
#=GR MW160818.1 PP ***************************************......................................
MZ516105.1         GTCTAAAACTAACAATCACACATGTGCATTTACAACACAACGAGACATTAGTTTTTGACACTTTT-TTTC--TCGT-
#=GR MZ516105.1 PP *****************************************************************.****..****.
ON237084.1         GTCTAAAACTAACAATCACACATGTGCATTTACAACACAACGAGACATTA---------------------------
#=GR ON237084.1 PP **************************************************...........................
ON237173.1         GTCTAAAACTAACAATCACACATGTGCATTTACAACACAACGAGACATTA---------------------------
#=GR ON237173.1 PP **************************************************...........................
ON237174.1         GTCTAAAACTAACAATCACACATGTGCATTTACAACACAACGAGACATTA---------------------------
#=GR ON237174.1 PP **************************************************...........................
ON237177.1         GTCTAAAACTAACAATCACACATGTGCATTTACAACACAACGAGACATTA---------------------------
#=GR ON237177.1 PP **************************************************...........................
OR496332.1         GTCTAAAACTAACAATCACACATGTGCATTTAC--------------------------------------------
#=GR OR496332.1 PP *********************************............................................
#=GC SS_cons       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.::::..::::.
#=GC RF            GTCTAAAACTAACAATGATACATGTGCATTTACAACACAACGAGACATTAGTTTTTGACACTTTT.TTTC..TCGT.
//
```

### Build new models from new reference sequence(s)

To build our new models, we run `v-build.pl`:

```
$ v-build.pl --group RSV --subgroup A KY654518 KY654510
$ v-build.pl --group RSV --subgroup B MZ516105 MZ516105 
```

These models will take up to an hour to build.
When they're finished, we can combine them as before:

```
# create a new directory
$ mkdir rsv-models2

# concatenate .minfo, .cm .fa and .hmm files:
$ cat KY654518/*.vadr.minfo > rsv-models2/rsv.minfo
$ cat KY654518/*.vadr.cm > rsv-models2/rsv.cm
$ cat KY654518/*.vadr.fa > rsv-models2/rsv.fa
$ cat KY654518/*.vadr.protein.hmm > rsv-models2/rsv.hmm
$ cat MZ516105/*.vadr.minfo >> rsv-models2/rsv.minfo
$ cat MZ516105/*.vadr.cm >> rsv-models2/rsv.cm
$ cat MZ516105/*.vadr.fa >> rsv-models2/rsv.fa
$ cat MZ516105/*.vadr.protein.hmm >> rsv-models2/rsv.hmm

# copy the blastdb files:
$ cp KY654518/*.vadr.protein.fa* rsv-models2/
$ cp MZ516105/*.vadr.protein.fa* rsv-models2/

# prepare the library files:
$ $VADRINFERNALDIR/esl-sfetch --index rsv-models2/rsv.fa
$ $VADRINFERNALDIR/cmpress rsv-models2/rsv.cm
$ $VADRHMMERDIR/hmmpress rsv-models2/rsv.hmm
$ $VADRBLASTDIR/makeblastdb -dbtype nucl -in rsv-models2/rsv.fa
```

As in iteration 1, it's a good idea to run `v-annotate.pl` with the new models against the two model
sequences as a sanity check. For these two RSV models, both sequences
should pass:

```
$ v-annotate.pl --mdir rsv-models2 --mkey rsv rsv-models2/rsv.fa va-rsv2
```

```
# Summary of classified sequences:
#
#                                 num   num   num
#idx  model     group  subgroup  seqs  pass  fail
#---  --------  -----  --------  ----  ----  ----
1     KY654518  RSV    A            1     1     0
2     MZ516105  RSV    B            1     1     0
#---  --------  -----  --------  ----  ----  ----
-     *all*     -      -            2     2     0
-     *none*    -      -            0     0     0
#---  --------  -----  --------  ----  ----  ----
#
# Zero alerts were reported.
#
```

</details>

<details>

<summary>

## Iteration 2, step 2: run `v-annotate.pl` on our existing training set

</summary>

We can reuse the training set from iteration 1 here in the second
iteration. We could create a different training set but if our
original training set was sufficiently large and truly random, then we
probably don't need to. One upside of keeping the same training set is
that it makes differences between results with our new models and
results with our original models in iteration 1 easier to understand
and interpret.

We can repeat the `v-annotate.pl` command from iteration 1, but this
time we will use the `--out_stk` option to save multiple alignments
for reasons that will become clear below:

```
$ v-annotate.pl --out_stk --mdir rsv-models2 --mkey rsv rsv.r500fa va2-r500
```
</details>

<details>

<summary>

## Iteration 2, step 3: analyze the results and update the models accordingly

</summary>

This time, from the `v-annotate.pl` output we can tell that many more
sequences passed than in iteration 1, when only 6 of 500 passed:

```
# Summary of classified sequences:
#
#                                 num   num   num
#idx  model     group  subgroup  seqs  pass  fail
#---  --------  -----  --------  ----  ----  ----
1     KY654518  RSV    A          286   135   151
2     MZ516105  RSV    B          214   132    82
#---  --------  -----  --------  ----  ----  ----
-     *all*     -      -          500   267   233
-     *none*    -      -            0     0     0
#---  --------  -----  --------  ----  ----  ----
```

But there are still 233 sequences that do not pass. We know from
iteration 1 that we don't expect any of the fatal alerts remaining in
these 233 sequences to exist in a majority of sequences, but there
still may be some alerts that are somewhat common and exist in a
significant number of sequences. Many of these alerts may correspond to real
biological variability in RSV sequences that we'd rather not cause a
sequence to fail, and they may be addressable by modifying the model
without changing the underlying reference sequences like we did at the
end of iteration 1. There are several ways we can address these
situations, including:

1. *Add a new protein to the blastx protein library for a model to
   account for protein sequence variability.* This
   can help remove unnecessary alerts related to the protein validation stage, most
   commonly: `indfpst5`, `indfpst3`, `insertnp` and `deletinp`.

2. *Add **alternative features** to the model info file, along with
   corresponding proteins to the blastx protein library.* Some CDS may
   have multiple positions, with respect to the reference model, for the
   start and/or stop codon. We can deal with this by adding all of the
   acceptable alternatives as features to the model info file, and
   `v-annotate.pl` will attempt to annotate each of them and report
   annotation for the single alternative that yields the fewest fatal
   alerts. 

3. *Add an alert *exception* to the model info file.* For some alerts,
   exceptions can be added that prevent the reporting of alerts in
   specific model regions. For example, an exception can be added to
   prevent the reporting of a `dupregin` alert due to an expected
   repetitive region. 

4. *Rebuild the covariance model from a new input alignment.* We can
   still use the same reference sequence and coordinates, but by
   buildng a new model from an alignment with multiple sequences,
   certain alerts that caused by sequence differences with the
   reference model can be prevented.

5. *Specify features as non-essential*. By adding a a
   `misc_not_failure` flag to a feature in a feature info file, we can
   make it so that many types of fatal alerts do not cause a sequence
   to fail, but instead cause the associated feature to be annotated
   as a `misc_feature` in the output. 

6. *Specify command-line options when running `v-annotate.pl` to make
   some alerts fatal or non-fatal*. For some viruses, certain fatal alerts
   may be so common that they are expected for many sequences, yet we
   don't want them to cause a sequence to fail. We can make those
   alerts non-fatal using the command-line option `alt_pass`.

I will walk you through examples of some of these strategies. First, we
need to identify the minority characteristics that we would like to
address through model modification. Once again we can start by looking
at the most common types of reported alerts:


Here are the counts for all of the fatal alerts, from the
`va2-rsv.r500/va2-rsv.r500.vadr.alc` file:

```
#     alert     causes   short                               per    num   num  long       
#idx  code      failure  description                        type  cases  seqs  description
#---  --------  -------  -----------------------------  --------  -----  ----  -----------
16    lowcovrg  yes      LOW_COVERAGE                   sequence     13    13  low sequence fraction with significant similarity to homology model
17    lowsimis  yes      LOW_SIMILARITY                 sequence      3     3  internal region without significant similarity
18    mutstart  yes      MUTATION_AT_START               feature      6     6  expected start codon could not be identified
19    mutendcd  yes      MUTATION_AT_END                 feature    114   110  expected stop codon could not be identified, predicted CDS stop by homology is invalid
20    mutendns  yes      MUTATION_AT_END                 feature      2     2  expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted start codon
21    mutendex  yes      MUTATION_AT_END                 feature    110   109  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position
22    unexleng  yes      UNEXPECTED_LENGTH               feature     16    15  length of complete coding (CDS or mat_peptide) feature is not a multiple of 3
23    cdsstopn  yes      CDS_HAS_STOP_CODON              feature     36    35  in-frame stop codon exists 5' of stop position predicted by homology to reference
24    cdsstopp  yes      CDS_HAS_STOP_CODON              feature      3     3  stop codon in protein-based alignment
25    fsthicft  yes      POSSIBLE_FRAMESHIFT_HIGH_CONF   feature     12    12  high confidence possible frameshift in CDS (frame not restored before end)
26    fsthicfi  yes      POSSIBLE_FRAMESHIFT_HIGH_CONF   feature      1     1  high confidence possible frameshift in CDS (frame restored before end)
27    indfantn  yes      INDEFINITE_ANNOTATION           feature      5     3  nucleotide-based search identifies CDS not identified in protein-based search
28    indf5gap  yes      INDEFINITE_ANNOTATION_START     feature      4     2  alignment to homology model is a gap at 5' boundary
29    indf5pst  yes      INDEFINITE_ANNOTATION_START     feature     25    22  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint
30    indf3gap  yes      INDEFINITE_ANNOTATION_END       feature      6     3  alignment to homology model is a gap at 3' boundary
31    indf3pst  yes      INDEFINITE_ANNOTATION_END       feature    132   124  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
32    deletinp  yes      DELETION_OF_NT                  feature    138   138  too large of a deletion in protein-based alignment
```

Sorting by number of sequences, there are 9 alerts that occur in more
than 10 sequence (more than 2\% of sequences):

```
32    deletinp  yes      DELETION_OF_NT                  feature    138   138  too large of a deletion in protein-based alignment
31    indf3pst  yes      INDEFINITE_ANNOTATION_END       feature    132   124  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
19    mutendcd  yes      MUTATION_AT_END                 feature    114   110  expected stop codon could not be identified, predicted CDS stop by homology is invalid
21    mutendex  yes      MUTATION_AT_END                 feature    110   109  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position
23    cdsstopn  yes      CDS_HAS_STOP_CODON              feature     36    35  in-frame stop codon exists 5' of stop position predicted by homology to reference
29    indf5pst  yes      INDEFINITE_ANNOTATION_START     feature     25    22  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint
22    unexleng  yes      UNEXPECTED_LENGTH               feature     16    15  length of complete coding (CDS or mat_peptide) feature is not a multiple of 3
16    lowcovrg  yes      LOW_COVERAGE                   sequence     13    13  low sequence fraction with significant similarity to homology model
25    fsthicft  yes      POSSIBLE_FRAMESHIFT_HIGH_CONF   feature     12    12  high confidence possible frameshift in CDS (frame not restored before end)
```

We will investigate each of these types of alerts, starting with
`deletinp`. As in iteration 1, we can sort all the occurences of this
alert in the `.alt` file and group them:

```
$ grep deletinp va2-rsv.r500/va2-rsv.r500.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
     69 KY654518 attachment_glycoprotein 5461..5532:+
     23 KY654518 attachment_glycoprotein 5509..5574:+
     20 MZ516105 attachment_glycoprotein 5435..5494:+
     15 MZ516105 attachment_glycoprotein 5399..5458:+
      7 KY654518 attachment_glycoprotein 5515..5580:+
      1 MZ516105 attachment_glycoprotein 5444..5503:+
      1 MZ516105 attachment_glycoprotein 5390..5440:+
      1 KY654518 attachment_glycoprotein 5494..5547:+
      1 KY654518 attachment_glycoprotein 5458..5529:+
```

The deletions are occuring in the attachment glycoprotein CDS in a
region that may be familiar - it is close to the duplicate region we
observed in iteration 1. We attempted to address the `dupregin` alerts
in iteration 1 by rebuilding our models with new sequences that
*included* the duplicated region, but now it seems that we are observing the sequences
that do not have the duplication failing with a `deletinp`
alert. These are a minority of the sequences but still a significant
number in each model. 

The `deletinp` alert occurs when there is a region in the best
'blastx' alignment that includes a deletion that is longer than 9
amino acids (27 nucleotides). Let's take a look at one example of the
most common deletion span: `5461..5532:+` to the `KY654518` (RSV A)
model, by randoming selecting one sequence and rerunning
`v-annotate.pl` on it with the `--keep` option as we did earlier in
iteration 1. 

```
# pick a sequence with the deletinp alert:
$ grep deletinp va2-rsv.r500/va2-rsv.r500.vadr.alt | grep 5461..5532 | awk '{ print $2 }' | esl-selectn 1 - > ex6.list
$ cat ex6.list
KX655676.1

# fetch the sequence
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex6.list > ex6.fa

# run v-annotate.pl on these sequences with 
# --keep option to save all output files
$ $VADRSCRIPTSDIR/v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex6.fa va-ex6
```

```
# Summary of reported alerts:
#
#     alert     causes   short               per    num   num  long
#idx  code      failure  description        type  cases  seqs  description
#---  --------  -------  --------------  -------  -----  ----  -----------
1     deletinn  no       DELETION_OF_NT  feature      1     1  too large of a deletion in nucleotide-based alignment of CDS feature
#---  --------  -------  --------------  -------  -----  ----  -----------
2     deletinp  yes      DELETION_OF_NT  feature      1     1  too large of a deletion in protein-based alignment
#---  --------  -------  --------------  -------  -----  ----  -----------
```

The relevant file is the `blastx` output file
`va-ex6/va-ex6.vadr.KY654518.blastx.out`, and we are interested in the
`blastx` alignment for the attachment glycoprotein, which is positions
`4681..5646:+` as found in the `.minfo` file:

```
$ grep attachment rsv-models2/rsv.minfo | grep KY654518
FEATURE KY654518 type:"CDS" coords:"4681..5646:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein"
```

Here is the relevant alignment from 
`va-ex6/va-ex6.vadr.KY654518.blastx.out`:

```
>KY654518.1/4681..5646:+
Length=321

 Score = 545 bits (1404),  Expect = 6e-180, Method: Compositional matrix adjust.
 Identities = 287/321 (89%), Positives = 291/321 (91%), Gaps = 24/321 (7%)
 Frame = +1

Query  4609  MSKTKDQRTAKTLERTWDTLNHLLFISSCLYKLNLKSIAQITLSILAMIISTSLIIAAII  4788
             MSKTKDQRTAKTLERTWDTLNHLLFISSCLYKLNLKSIAQITLSILAMIISTSLIIAAII
Sbjct  1     MSKTKDQRTAKTLERTWDTLNHLLFISSCLYKLNLKSIAQITLSILAMIISTSLIIAAII  60

Query  4789  FIASANHKVTLTTAIIQDATNQIKNTTPTYLTQNPQLGISFTNLSGTTSKSTTILASTTP  4968
             FIASANHKVTLTTAIIQDATNQIKNTTPTYLTQNPQLGISF+NLSGTTS+STTILASTTP
Sbjct  61    FIASANHKVTLTTAIIQDATNQIKNTTPTYLTQNPQLGISFSNLSGTTSQSTTILASTTP  120

Query  4969  SAESTPQSTTVKIKNTTTTQIQPSKPTTKQRQNKPQNKPNNDFHFEVFNFVPCSICSNNP  5148
             SAESTPQSTTVKIKNTTTTQI PSKPTTKQRQNKPQNKPNNDFHFEVFNFVPCSICSNNP
Sbjct  121   SAESTPQSTTVKIKNTTTTQILPSKPTTKQRQNKPQNKPNNDFHFEVFNFVPCSICSNNP  180

Query  5149  TCWAICKRIPNKKPGKKTTTKPTKKPTIKTTKKDPKPQTTKPKEVLTTKPTEKPTIDTTK  5328
             TCWAICKRIPNKKPGKKTTTKPTKKPT+KTTKKDPKPQTTKPKEVLTTKPT KPTI+TTK
Sbjct  181   TCWAICKRIPNKKPGKKTTTKPTKKPTLKTTKKDPKPQTTKPKEVLTTKPTGKPTINTTK  240

Query  5329  TNIRTTPLTSNTTGNPEHTS------------------------QEETLHSTTSEGNLSP  5436
             TNIRTT LTSNT GNPEHTS                        QEETLHSTTSEG LSP
Sbjct  241   TNIRTTLLTSNTKGNPEHTSQEETLHSTTSEGYLSPSQVYTTSGQEETLHSTTSEGYLSP  300

Query  5437  SQVYTTSEYLSQSPSSSNTTK  5499
             SQVYTTSEYLSQS SSSNTTK
Sbjct  301   SQVYTTSEYLSQSLSSSNTTK  321
```

Note the 24 amino acid deletion in the query with respect to the
subject. This indicates that sequence `KX655676.1` has a 72nt deletion
relative to the reference sequence `KY654318.1`. The `blastx` library
created by `v-build.pl` only has a single protein sequence for the
attachment glycoprotein, the translation of the `KY654318.1` CDS from
positions `4681..5646:+`. But we can add additional proteins that will
then be used as additional subjects in the `blastx` protein validation
stage. 

One option for such a sequence would be to use the attachment
glycoprotein from `KX655676.1`, which would certainly fix the issue
for at least itself. A more well-principled way to find a better
sequence would be to identify a more representative sequence out of
those that include this deletion in our training set. (A yet more
well-principled way would be to inspect all existing RSV A sequences
instead of just those in our training set, but this would require
significantly more effort for a presumably small improvement over the
alternative of restricting possibilities to only our training set.)

We can extract the candidate sequences from the alignment that was
output from `v-annotate.pl` (due to the use of the `--out_stk`
option). We are interested in a representative example of the CDS with
the `deletinp` alert for the `KY654318` model, so it makes sense to
use a sequence with the most common deletion of reference positions
`5461..5532:+`. To get a list of the sequences and extract the relevant alignment subset: 

```
# get a list of the 69 candidates
$ cat va2-rsv.r500/va2-rsv.r500.vadr.alt | grep deletinp | grep KY654518 | grep attachment_glycoprotein | grep 5461..5532 | awk '{ printf("%s\n", $2); }' > ex7.list

# extract the 69 aligned sequences from the alignment:
$ $VADREASELDIR/esl-alimanip --seq-k ex7.list va2-rsv.r500/va2-rsv.r500.vadr.KY654518.align.stk > ex7.stk 

# extract the attachment glycoprotein region:
$ esl-alimask -t --t-rf ex7.stk 4681..5646 > ex7.ag.stk
```

Next, as we did in iteration 1 to find new reference sequences, we can
take the follow steps to identify our reference:

1. Remove all sequences with any ambiguous nucleotides using the
`count-ambigs.pl` script.

2. use the `esl-alipid` and the `esl-alipid-per-seq-stats.pl` script
   to find the sequence out of the candidates that has the highest average
   percent identity to all other candidate sequences.

```
# convert to fasta 
# $ $VADREASELDIR/esl-reformat fasta ex7.ag.stk > ex7.ag.fa
# remove any sequences with ambiguous nucleotides
# $ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl ex7.ag.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | awk '{ printf("%s\n", $1); }' > ex7.60.list
# $ $VADREASELDIR/esl-alimanip --seq-k ex7.60.list ex7.ag.stk > ex7.ag.60.stk

# determine pairwise percent identity
$ $VADREASELDIR/esl-alipid ex7.ag.60.stk > ex7.ag.60.alipid

# calculate average percent identity
$ perl $VADRSCRIPTSDIR/miniscripts/esl-alipid-per-seq-stats.pl ex7.ag.60.alipid > ex7.ag.60.alipid.perseq

# list top candidates
$ grep -v ^\# ex7.ag.alipid.perseq | sort -rnk 2 | head
KX510193.1  96.631  KF826854.1  91.390  KJ627315.1  99.890
KJ627315.1  96.593  KF826854.1  91.500  KX510193.1  99.890
KX510189.1  96.541  KF826854.1  91.280  KX510193.1  99.890
KJ627320.1  96.521  KF826854.1  91.280  KX510193.1  99.890
OK649655.1  96.400  KF826854.1  91.280  KX510193.1  99.660
KX510264.1  96.357  KF826854.1  91.050  KX510250.1  100.000
KX510250.1  96.357  KF826854.1  91.050  KX510264.1  100.000
KX510230.1  96.357  KF826854.1  91.050  KX510250.1  100.000
KX510195.1  96.357  KF826854.1  91.050  KX510250.1  100.000
KX510148.1  96.357  KF826854.1  91.050  KX510250.1  100.000
```

We will use the `KX510193.1` attachment glycoprotein sequence, which
has the highest average nucleotide percent identity with all other
candidates. (We could have used a protein alignment for this step,
especially since we are trying to optimize blastx alignments, but the
nucleotide based approach we've taken here should be a good proxy for
the protein-based approach.)

To actually add the new sequence to our blastx library, we will use
the `build-add-to-blast-db.pl` script. To determine how to use that
script, execute it with the `-h` option and no other command-line arguments:

```
$ perl $VADRSCRIPTSDIR/miniscripts/build-add-to-blast-db.pl -h
# build-add-to-blast-db.pl :: add a single protein to a VADR blastx protein database
# VADR 1.6 (Nov 2023)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# date:    Mon Oct 16 15:23:37 2023
#
Usage: build-add-to-blast-db.pl [-options]
	<path to .minfo file>
	<path to blast db dir>
	<model name>
	<nt-accn-to-add>
	<nt-coords-to-add>
	<model-CDS-feature-coords>
	<name for output directory>

basic options:
  -f         : force; if dir <output directory> exists, overwrite it
  -v         : be verbose; output commands to stdout as they're run
  --ttbl <n> : use NCBI translation table <n> to translate CDS [1]
  --keep     : do not remove intermediate files, keep them all on disk
```

This script takes seven command-line arguments, we already know all of
them except for one: `<nt-coords-to-add>`. These are the nucleotide
coordinates of the attachment glycoprotein CDS in the `KY510193.1`
sequence. We can find these in the `.ftr` or `.tbl` output files from
`v-annotate.pl`:

```
$ head -n 2 va2-rsv.r500/va2-rsv.r500.vadr.ftr 
#      seq           seq                  ftr   ftr                         ftr  ftr  par                                                                                                 seq          model  ftr   
#idx   name          len  p/f   model     type  name                        len  idx  idx  str  n_from   n_to  n_instp  trc  5'N  3'N  p_from   p_to  p_instp   p_sc  nsa  nsn         coords         coords  alerts
$ grep KX510193 va2-rsv.r500/va2-rsv.r500.vadr.ftr | grep attachment
286.14  KX510193.1  14701  FAIL  KY654518  CDS   attachment_glycoprotein     894   14   -1    +    4409   5302        -  no     0    0    4409   5299        -   1420    1    0   4409..5302:+   4681..5646:+  DELETION_OF_NT(deletinn),DELETION_OF_NT(deletinp)
```

The relevant field is the `seq coords` field, so the relevant value is
`4409..5302:+`.

So all of the relevant values for all arguments are:

| command-line argument    | value |
|--------------------------|-------------------------|
| `<path to .minfo file>`  | `rsv-models2/rsv.minfo` | 
| `<path to blast db dir>` | `rsv-models2/`          |
| `<model name>`           | `KY654518`              | 
| `<nt-accn-to-add>`       | `KX510193`              |
| `nt-coords-to-add>`      | `4409..5302:+`          |
| `<model-CDS-feature-coords> | `4681..5646:+`       | 
| <name for output directory> | `vb-ex7`             | 

To add the protein:
```
$ perl $VADRSCRIPTSDIR/miniscripts/build-add-to-blast-db.pl \
rsv-models2/rsv.minfo \
rsv-models2 \
KY654518 \
KX510193 \
4409..5302:+ \
4681..5646:+ \
vb-ex7
```

The script will output the steps it takes:
```
# input model info file:                      rsv-models2/rsv.minfo
# input blast db path:                        rsv-models2
# input model name:                           KY654518
# nucleotide accession to add:                KX510193
# nt coords of CDS to add:                    4409..5302:+
# CDS feature coords this CDS should map to:  4681..5646:+
# output directory:                           vb-ex7
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parsing input model info file                      ... done. [    0.0 seconds]
# Fetching the CDS source sequence                   ... done. [    4.2 seconds]
# Translating CDS                                    ... done. [    0.0 seconds]
# Adding to BLAST DB                                 ... done. [    0.2 seconds]
#
# Output printed to screen saved in:                           vb-ex7.vadr.log
# List of executed commands saved in:                          vb-ex7.vadr.cmd
# List and description of all output files saved in:           vb-ex7.vadr.filelist
# fasta file with source sequence (KX510193) saved in:         vb-ex7.vadr.source.fa
# fasta file with CDS from KX510193 saved in:                  vb-ex7.vadr.cds.fa
# fasta file with translated protein from KX510193 saved in:   vb-ex7.vadr.prot.fa
#
# All output files created in the current working directory
#
# Elapsed time:  00:00:04.53
#                hh:mm:ss
# 
[ok]
```

We can verify the new protein was added by examining the protein
database fasta file using the `esl-seqstat` program with the `-a`
option which lists each sequence and its length:

```
$ esl-seqstat -a rsv-models2/KY654518.vadr.protein.fa
= KY654518.1/1140..2315:+        391 
= KY654518.1/2347..3072:+        241 
= KY654518.1/3255..4025:+        256 
= KY654518.1/4295..4489:+         64 
= KY654518.1/4681..5646:+        321 
= KY654518.1/5726..7450:+        574 
= KY654518.1/628..1002:+         124 
= KY654518.1/7669..8253:+        194 
= KY654518.1/8228..8494:+         88 
= KY654518.1/8561..15058:+      2165 
= KY654518.1/99..518:+           139 
= KX510193.1:4409..5302:+/4681..5646:+      297 
Format:              FASTA
Alphabet type:       amino
Number of sequences: 12
Total # residues:    4854
Smallest:            64
Largest:             2165
Average length:      404.5
```

This file includes all the protein sequences for the
model. `v-annotate.pl` uses the coordinates at the end of each name to
determine which sequence pertains to which protein, by matching those
coordinates up with the coordinates read for each CDS feature in the
model info file. Note that now there are two sequences that end with
`4681..5646:+`: the sequence `KY654518.1/4681..5646:+`, which is the
original translated CDS from `KY654518.1` and
`KX510193.1:4409..5302:+/4681..5646:+` which is the sequence we just
added. 

We can then perform a sanity check to make sure that this added
sequence has the intended effect. Let's run `v-annotate.pl` on our
`ex6.fa` sequence. It should now pass. 

```
$ $VADRSCRIPTSDIR/v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex6.fa va-ex6
```

```
# Summary of classified sequences:
#
#                                 num   num   num
#idx  model     group  subgroup  seqs  pass  fail
#---  --------  -----  --------  ----  ----  ----
1     KY654518  RSV    A            1     1     0
#---  --------  -----  --------  ----  ----  ----
-     *all*     -      -            1     1     0
-     *none*    -      -            0     0     0
#---  --------  -----  --------  ----  ----  ----
#
# Summary of reported alerts:
#
#     alert     causes   short               per    num   num  long
#idx  code      failure  description        type  cases  seqs  description
#---  --------  -------  --------------  -------  -----  ----  -----------
1     deletinn  no       DELETION_OF_NT  feature      1     1  too large of a deletion in nucleotide-based alignment of CDS feature
#---  --------  -------  --------------  -------  -----  ----  -----------
```

The sequence now passes. We still have a `deletinn` alert letting us
know that there is still a long deletion in the nucleotide-based
alignment, but this is a non-fatal alert. The `deletinp` alert is now
gone.

At this point, we could continue to address the `deletinp` instances,
probably first for the 35 `MZ516105` alerts. To do that, we would
repeat the above procedure to find a suitable representative sequence
to add to the protein blast library. After that, the next step would
be to rerun all of the sequences that failed due to `deletinp` alerts
with the updated models. If a significant number of sequences still
fail due to `deletinp` alerts at that stage, then we could repeat the
process again. 

For the purposes of this tutorial, we will move on to the next most
common alert `indf3pst` to provide a slightly different example of
updating a model.


Let's take a look at one example of this alert:
```
$ cat va2-rsv.r500/*alt | head -n 3
#        seq                   ftr   ftr                      ftr  alert           alert                                     seq   seq             mdl   mdl  alert 
#idx     name        model     type  name                     idx  code      fail  description                            coords   len          coords   len  detail
#------  ----------  --------  ----  -----------------------  ---  --------  ----  -----------------------------  --------------  ----  --------------  ----  ------
$ cat va2-rsv.r500/*alt | grep indf3pst | head -n 1
3.1.4    KY674983.1  MZ516105  CDS   attachment_glycoprotein   14  indf3pst  yes   INDEFINITE_ANNOTATION_END        5513..5518:+     6    5620..5620:+     1  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint [6>5, no valid stop codon in nucleotide-based prediction]

```

This alert occurs if the "protein-based alignment does not extend
close enough to nucleotide-based alignment 3' endpoint". A relevant
field is field 26, which explains how far the protein-based endpoint
and nucleotide based endpoint is. In this case it is 6 nucleotides
which exceeds the maximum allowed without an alert: `6>5`.

When grouping instances of this alert we should output this value as
well. Any instances for the same feature and the same distance may be
able to be addressed with the same model modification:

```
$ grep indf3pst va2-rsv.r500/va2-rsv.r500.vadr.alt | awk '{ printf ("%s %s %s %s\n", $3, $5, $12, $26); }' | sort | uniq -c | sort -rnk 1
     49 MZ516105 attachment_glycoprotein 5620..5620:+ [6>5,
     25 KY654518 attachment_glycoprotein 5646..5646:+ [6>5,
     15 KY654518 attachment_glycoprotein 5646..5646:+ [9>8,
     10 KY654518 attachment_glycoprotein 5646..5646:+ [45>5,
      7 KY654518 attachment_glycoprotein 5646..5646:+ [45>8,
      1 MZ516105 small_hydrophobic_protein 4498..4498:+ [132>120,
      1 MZ516105 polymerase 15060..15060:+ [25>5]
      1 MZ516105 polymerase 15060..15060:+ [24>5,
      1 MZ516105 polymerase 15060..15060:+ [2316>5,
      1 MZ516105 polymerase 15060..15060:+ [22>5]
      1 MZ516105 polymerase 15060..15060:+ [2078>5,
..snip..
```

It turns out the example we looked at about in `attachment
glycoprotein` for model `MZ516105` with a difference of 6 nucleotides
is the most common one. Let's investigate one of those sequences further:

```
$ grep indf3pst va2-rsv.r500/va2-rsv.r500.vadr.alt | grep MZ516105 | grep 5620 | awk '{ print $2 }' | esl-selectn 1 - > ex8.list
$ cat ex8.list
OR326763.1

# re-run v-annotate.pl on:
$ $VADRSCRIPTSDIR/v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex8.fa va-ex8
```
# Summary of reported alerts:
#
#     alert     causes   short                          per    num   num  long
#idx  code      failure  description                   type  cases  seqs  description
#---  --------  -------  -------------------------  -------  -----  ----  -----------
1     mutendcd  yes      MUTATION_AT_END            feature      1     1  expected stop codon could not be identified, predicted CDS stop by homology is invalid
2     mutendex  yes      MUTATION_AT_END            feature      1     1  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position
3     indf3pst  yes      INDEFINITE_ANNOTATION_END  feature      1     1  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
#---  --------  -------  -------------------------  -------  -----  ----  -----------
```

In this case, we find that this sequence not only has the common
`indf3pst` alert, but also `mutendcd` and `mutendex` alerts that occur
when the expected stop codon is missing, and the closest stop codon
occurs downstream of the expected position. This suggests that perhaps
many of those 49 sequences with this same alert also have these
alerts. To test that we can run `v-annotate.pl` on all 49 of them, or
if we want to save time, on a subset of 10:

```
$ grep indf3pst va2-rsv.r500/va2-rsv.r500.vadr.alt | grep MZ516105 | grep 5620 | awk '{ print $2 }' | esl-selectn 10 - > ex8.10.list
$ $VADRSCRIPTSDIR/v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex8.10.fa va-ex8.10
```

After this finishes, we can see that these three alerts do tend to
co-occur by looking at the `.alc` file:
```
$ cat *alc
#     alert     causes   short                          per    num   num  long       
#idx  code      failure  description                   type  cases  seqs  description
#---  --------  -------  -------------------------  -------  -----  ----  -----------
1     deletinn  no       DELETION_OF_NT             feature      3     3  too large of a deletion in nucleotide-based alignment of CDS feature
#---  --------  -------  -------------------------  -------  -----  ----  -----------
2     mutendcd  yes      MUTATION_AT_END            feature     10    10  expected stop codon could not be identified, predicted CDS stop by homology is invalid
3     mutendex  yes      MUTATION_AT_END            feature     10    10  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position
4     indf3pst  yes      INDEFINITE_ANNOTATION_END  feature     10    10  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
5     deletinp  yes      DELETION_OF_NT             feature      3     3  too large of a deletion in protein-based alignment
#---  --------  -------  -------------------------  -------  -----  ----  -----------
```

Let's go back to our single example `OR326763.1`, and look at details
on the alerts in the `.alt` file:
```
<[(tutorial-20231006)]> cat va-ex8/va-ex8.vadr.alt
#      seq                   ftr   ftr                      ftr  alert           alert                               seq  seq           mdl  mdl  alert 
#idx   name        model     type  name                     idx  code      fail  description                      coords  len        coords  len  detail
#----  ----------  --------  ----  -----------------------  ---  --------  ----  -------------------------  ------------  ---  ------------  ---  ------
1.1.1  OR326763.1  MZ516105  CDS   attachment_glycoprotein   14  mutendcd  yes   MUTATION_AT_END            5569..5571:+    3  5618..5620:+    3  expected stop codon could not be identified, predicted CDS stop by homology is invalid [CAA]
1.1.2  OR326763.1  MZ516105  CDS   attachment_glycoprotein   14  mutendex  yes   MUTATION_AT_END            5590..5592:+    3  5639..5641:+    3  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position [TAG]
1.1.3  OR326763.1  MZ516105  CDS   attachment_glycoprotein   14  indf3pst  yes   INDEFINITE_ANNOTATION_END  5566..5571:+    6  5620..5620:+    1  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint [6>5, no valid stop codon in nucleotide-based prediction]
```

For `mutendcd` and `mutendex` the alert detail field explains that in
this sequence there is a `CAG` at the expected stop codon position of `5618..5620`,
and the first in-frame stop codon `TAG` occurs 21 nucleotides
downstream at reference positions `5639..5641`.
We can see this in the `va-ex8/va-ex8.vadr.MZ516105.align.stk` alignment:


```
                                    vvv                  vvv    
OR326763.1         CCATATCAAATTCCACCCAAATACTCCAGTCATATGCTTAGTTATTTAAAAACTACATC
#=GR OR326763.1 PP ***********************************************************
#=GC SS_cons       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF            CCACATCAAATTCTATCTAAAGACTCCAGTCATATGCTTAGTTATTTAAAAACTACATC
#=GC RFCOLX....    00000000000000000000000000000000000000000000000000000000000
#=GC RFCOL.X...    55555555555555555555555555555555555555555555555555555555555
#=GC RFCOL..X..    66666666666666666666666666666666666666666666666666666666666
#=GC RFCOL...X.    00000000011111111112222222222333333333344444444445555555555
#=GC RFCOL....X    12345678901234567890123456789012345678901234567890123456789
                   *** ********* * * ***** ***********************************
```

I've added `vvv` characters indicating the expected stop codon
position and also the existing stop codon 9 nucleotides
downstream. I've also added `*` characters at the bottom of the
alignment at positions where the sequence and reference model are
identical. Note that the end of the CDS has the highest number of
mismatches, which explains the `indf3pst` alert for this CDS.

Because this is such a common characteristic of RSV B sequences, we'd
like our model to allow for it and not report any fatal alerts when it
occurs. To do this, we can modify our model by adding an **alternative
feature** for the attachment glycoprotein CDS. 

Adding an alternative feature for a CDS requires two steps:

Step 1. Manually edit the model info file
`rsv-models2/rsv.minfo` in a text editor. 

Step 2. Add one (or more) protein sequences to the blastx protein
     library that correspond to the alternative CDS feature, like we
     did to address the common `deletinp` alerts above.

In step 1, we want to make an alternative stop position for the
attachment glycoprotein CDS. To do this we will add an additional CDS
feature to the model info file. Currently the two lines pertaining to
the CDS and associated gene feature are:


```
FEATURE MZ516105 type:"gene" coords:"4688..5620:+" parent_idx_str:"GBNULL" gene:"G"
FEATURE MZ516105 type:"CDS" coords:"4688..5620:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein"
```

To create an alternative CDS feature that starts at the same position
`4688` but ends at position `5641` we will add the line:

```
FEATURE MZ516105 type:"CDS" coords:"4688..5641:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)"
```

Note the different stop codon *and* the new key/value pair:
`alternative_ftr_set="attachment(cds)"`. We also need to update the
first line above to include this field. This will inform
`v-annotate.pl` that these two CDS are members of the same alternative
feature set, and only 1 of them should be annotated in the output
`.tbl` file. `v-annotate.pl` will select the CDS that has fewer fatal
alerts and annotate that one. If they have the same number of fatal
alerts, the one that occurs first in the model info file will be
annotated. 

We also want to add a new `gene` feature line that will be coupled
with the new `CDS` feature line (same coordinates) so that the correct
gene coordinates will be annotated based on which CDS is annotated. 

```
FEATURE MZ516105 type:"gene" coords:"4688..5641:+" parent_idx_str:"GBNULL" gene:"G" alternative_ftr_set="attachment(gene)" alternative_ftr_set_subn:"attachment(cds).1"
```

This new `gene` feature line includes 
`alternative_ftr_set="attachment(gene)"` (it is important that this is
different from the CDS feature set name), and an additional key/value
pair: `alternative_ftr_set_subn="attachment(CDS).2"`. This means that
this `gene` should only be annotated if the second feature in the
`alternative_ftr_set="attachment(CDS)"` group is annotated. 

Additionally, we need to update the two original lines for the
features at positions `4688..5620:+` by adding `alternative_ftr_set`
key/values to them as well. The four new lines should be: 

```
FEATURE MZ516105 type:"gene" coords:"4688..5620:+" parent_idx_str:"GBNULL" gene:"G" alternative_ftr_set:"attachment(gene)" alternative_ftr_set_subn:"attachment(cds).1"
FEATURE MZ516105 type:"gene" coords:"4688..5641:+" parent_idx_str:"GBNULL" gene:"G" alternative_ftr_set:"attachment(gene)" alternative_ftr_set_subn:"attachment(cds).2"
FEATURE MZ516105 type:"CDS" coords:"4688..5620:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)"
FEATURE MZ516105 type:"CDS" coords:"4688..5641:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)"
```

Note that the order is important because the index in the
`alternative_ftr_set` values 
`attachment(cds).1` and `attachment(cds).2` for the two `gene`
features correspond to specifically to the first and second CDS features
that have the `alternative_ftr_set` value `attachment(cds)`.

Now we can move onto step 2, adding a protein to the blastx protein
library. As we did above when addressing the common `deletinp` alert,
we want to find a representative sequence out of the CDS that stop at
position `5641`. I repeated the steps detailed above using the set of 
41 sequences that end at this position as candidates and ended up
choosing `OK654726.1` as the representative. I then added it to the
blastx library using the `build-add-to-blast-db.pl` script. The steps
are shown below (and discussed in more detail for the `deletinp`
alert example above). One difference here is that we need to make sure
that we include the stop position in `OK654726.1` that aligns to the
new stop position at reference position `5641` which differs from what
is reported in the `.ftr` file. We may need to consult the alignment
of `OK654726.1` to determine this (after maybe rerunning
`v-annotate.pl`). 

```
# determine number of sequences that have attachment glycoprotein ended at 5641 
$ grep 5641 va2-rsv.r500/*alt | grep MZ516105 | grep mutendex | wc -l

# fetch out 41 sequences and extract only the attachment glycoprotein alignment
$ $VADREASELDIR/esl-alimanip --seq-k ex8.41.list va2-rsv.r500/va2-rsv.r500.vadr.MZ516105.align.stk > ex8.41.stk 
$ esl-alimask -t --t-rf ex8.41.stk 4688..5641 > ex8.41.ag.stk

# convert to fasta and remove any seqs with ambiguous nts
$ $VADREASELDIR/esl-reformat fasta ex8.41.ag.stk > ex8.41.ag.fa
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl ex8.41.ag.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | awk '{ printf("%s\n", $1); }' > ex8.40.list
$ $VADREASELDIR/esl-alimanip --seq-k ex8.40.list ex8.41.ag.stk > ex8.40.stk

# determine average percent id and choose representative
$ $VADREASELDIR/esl-alipid ex8.40.stk > ex8.40.alipid
$ perl $VADRSCRIPTSDIR/miniscripts/esl-alipid-per-seq-stats.pl ex8.40.alipid > ex8.40.alipid.perseq
$ grep -v ^\# ex8.40.alipid.perseq | sort -rnk 2 | head
OK649726.1  97.727  MG642027.1  92.780  OK649687.1  99.480

# add representative to blastx db
$ perl $VADRSCRIPTSDIR/miniscripts/build-add-to-blast-db.pl  \
rsv-models2/rsv.minfo \
rsv-models2 \
MZ516105 \
OK649726 \
4680..5633:+ \
4688..5641:+ \
vb-ex8
```

As a sanity check, we can rerun `v-annotate.pl` on our `ex8` sequence
`OR326763.1`, it should now pass:

```
$ $VADRSCRIPTSDIR/v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex8.fa va-ex8.2
```

```
# Summary of classified sequences:
#
#                                 num   num   num
#idx  model     group  subgroup  seqs  pass  fail
#---  --------  -----  --------  ----  ----  ----
1     MZ516105  RSV    B            1     1     0
#---  --------  -----  --------  ----  ----  ----
-     *all*     -      -            1     1     0
-     *none*    -      -            0     0     0
#---  --------  -----  --------  ----  ----  ----
#
# Zero alerts were reported.

```

---
TOADD: 
- add limitations/criticisms/alternatives to this approach:
  * the way we select training seqs, they may be very redundant, could
    weight them somehow would be better
  * experts may have a favorite reference sequence, they should use
    that
- what about multiple models? If we have many minor characteristics
   and only some groups have some chracteristics, it may be better to
   make multiple models


#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.

