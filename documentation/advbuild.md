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

1. Build a model(s) from representative and well annotated sequences
   as a starting point. These may be RefSeq sequences.

2. Construct a training set of randomly chosen existing sequences for
   the viral species (or use the same training set from the previous
   iteration) and use `v-annotate.pl` to validate and annotate
   the training sequences using your models from step 1.

3. Analyze the results by looking for common failure modes, and
   investigate the sequence characteristics that are responsible for
   them. This will be characteristics not present in the reference
   sequences. Classify these characteristics into major and minor types: 

   Major: these are common to a majority or large fraction of
   sequences but don't exist in the reference sequence. To address
   these, we probably want to pick a new representative sequence(s)
   that includes these characteristics. This means we will return to
   step 1 for another iteration of the three steps.

   Minor: these exist in a smaller fraction of sequences, and
   don't exist in the reference sequence. For these characteristics,
   we may be able to update our models without rebuilding them from
   new reference sequences. 

If the initial sequences chosen in the first iteration turn out to
representative enough that there are no major characteristics in step
3, then we can likely stop at the end of iteration 1 and update the
model(s) to tolerate the minor characteristics. 

But if there are one or more major characteristics identified in step
3, we will want to do an additional iteration of all three steps.  Two
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
</details>

<details>

<summary>

## Iteration 1, step 2: construct a training set and run `v-annotate.pl` on it

</summary>

ADD STEP OF RUNNING MODELS AGAINST MODEL SEQUENCES AS A
SANITY CHECK. I WILL THEN REFER TO THIS IN THE dupregin ALERT SECTION 

--- 

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
$ v-annotate.pl --mdir rsv-models --mkey rsv sequences_20231010_2146812.fasta va-r500
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
(also saved to the `va-r500/va-r500.vadr.log` file), we can see that
286 sequences were classified as RSV A (matched best to the
`NC_038235` model) and 214 as RSV B (matched best to the `NC_001781`
model). And only 6 out of the 500 total training sequences "pass":

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
5     insertnn  no       INSERTION_OF_NT                 feature    403   344  too large of an insertion in nucleotide-based alignment of CDS feature
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
16    dupregin  yes      DUPLICATE_REGIONS              sequence    343   343  similarity to a model region occurs more than once
17    deletins  yes      DELETION_OF_FEATURE            sequence      2     1  internal deletion of a complete feature
18    mutstart  yes      MUTATION_AT_START               feature    261   260  expected start codon could not be identified
19    mutendcd  yes      MUTATION_AT_END                 feature      8     8  expected stop codon could not be identified, predicted CDS stop by homology is invalid
20    mutendns  yes      MUTATION_AT_END                 feature      2     2  expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted start codon
21    mutendex  yes      MUTATION_AT_END                 feature      5     5  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position
22    unexleng  yes      UNEXPECTED_LENGTH               feature     19    17  length of complete coding (CDS or mat_peptide) feature is not a multiple of 3
23    cdsstopn  yes      CDS_HAS_STOP_CODON              feature    418   401  in-frame stop codon exists 5' of stop position predicted by homology to reference
24    cdsstopp  yes      CDS_HAS_STOP_CODON              feature    163   163  stop codon in protein-based alignment
25    fsthicft  yes      POSSIBLE_FRAMESHIFT_HIGH_CONF   feature     14    14  high confidence possible frameshift in CDS (frame not restored before end)
26    fsthicfi  yes      POSSIBLE_FRAMESHIFT_HIGH_CONF   feature      1     1  high confidence possible frameshift in CDS (frame restored before end)
27    indfantn  yes      INDEFINITE_ANNOTATION           feature      5     3  nucleotide-based search identifies CDS not identified in protein-based search
28    indf5gap  yes      INDEFINITE_ANNOTATION_START     feature      1     1  alignment to homology model is a gap at 5' boundary
29    indf5lcn  yes      INDEFINITE_ANNOTATION_START     feature     11     8  alignment to homology model has low confidence at 5' boundary for feature that does not match a CDS
30    indf5pst  yes      INDEFINITE_ANNOTATION_START     feature     36    32  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint
31    indf3gap  yes      INDEFINITE_ANNOTATION_END       feature     37    36  alignment to homology model is a gap at 3' boundary
32    indf3lcn  yes      INDEFINITE_ANNOTATION_END       feature   1057   320  alignment to homology model has low confidence at 3' boundary for feature that does not match a CDS
33    indf3pst  yes      INDEFINITE_ANNOTATION_END       feature    217   202  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
34    insertnp  yes      INSERTION_OF_NT                 feature    205   205  too large of an insertion in protein-based alignment
35    lowsim5n  yes      LOW_FEATURE_SIMILARITY_START    feature      1     1  region overlapping annotated feature that does not match a CDS at 5' end of sequence lacks significant similarity
36    lowsim5l  yes      LOW_FEATURE_SIMILARITY_START    feature      2     2  long region overlapping annotated feature that does not match a CDS at 5' end of sequence lacks significant similarity
37    lowsim3n  yes      LOW_FEATURE_SIMILARITY_END      feature      1     1  region overlapping annotated feature that does not match a CDS at 3' end of sequence lacks significant similarity
38    lowsim3l  yes      LOW_FEATURE_SIMILARITY_END      feature      1     1  long region overlapping annotated feature that does not match a CDS at 3' end of sequence lacks significant similarity
39    lowsimin  yes      LOW_FEATURE_SIMILARITY          feature     14    14  region overlapping annotated feature that does not match a CDS lacks significant similarity
40    lowsimil  yes      LOW_FEATURE_SIMILARITY          feature     91    31  long region overlapping annotated feature that does not match a CDS lacks significant similarity
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
23    cdsstopn  yes      CDS_HAS_STOP_CODON              feature    418   401  in-frame stop codon exists 5' of stop position predicted by homology to reference
16    dupregin  yes      DUPLICATE_REGIONS              sequence    343   343  similarity to a model region occurs more than once
32    indf3lcn  yes      INDEFINITE_ANNOTATION_END       feature   1057   320  alignment to homology model has low confidence at 3' boundary for feature that does not match a CDS
18    mutstart  yes      MUTATION_AT_START               feature    261   260  expected start codon could not be identified
34    insertnp  yes      INSERTION_OF_NT                 feature    205   205  too large of an insertion in protein-based alignment
33    indf3pst  yes      INDEFINITE_ANNOTATION_END       feature    217   202  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
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
$ grep cdsstopn va-sequences_20231010_2146812.vadr.alt
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
$ grep cdsstopn va-r500/va-r500.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
    213 NC_038235 attachment_glycoprotein 5579..5581:+
    153 NC_001781 attachment_glycoprotein 5566..5568:+
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
$ grep NC\_038235 rsv-models/rsv.minfo | grep attachment
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
$ grep NC\_001781 rsv-models/rsv.minfo | grep attachment
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
16    dupregin  yes      DUPLICATE_REGIONS              sequence    343   343  similarity to a model region occurs more than once
```

We can use `grep` to look at some examples of the `dupregin` alert:
```
$ grep dupregin va-r500/va-r500.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
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
<[(tutorial-20231006)]> grep dupregin va-sequences_20231010_2146812/va-sequences_20231010_2146812.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $23); }' | sort | uniq -c | sort -rnk 1
    167 NC_001781 - [5410..5469:+
    158 NC_038235 - [5466..5537:+
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

About 55% (158/286) or RSV A and 80% (167/214) of RSV B sequences have
one of two duplicated regions. To get a better idea of this
duplicatedd region, we can select an example of one RSV A sequence and
one RSV B sequence that contain this alert, and run it back through
`v-annotate.pl`. To select two sequences:

```
# prepare the sequence file for the esl-sfetch program 
$ esl-sfetch --index rsv.r500.fa

# pick a sequence with a dupregin alert that contains the strings we are interested it
$ grep dupregin va-sequences_20231010_2146812/va-sequences_20231010_2146812.vadr.alt | grep 5410..5469 | awk '{ print $2 }' | esl-selectn 1 - > ex1.list
$ grep dupregin va-sequences_20231010_2146812/va-sequences_20231010_2146812.vadr.alt | grep 5466..5537 | awk '{ print $2 }' | esl-selectn 1 - > ex2.list
$ cat ex1.list 
MZ516003.1
$ cat ex2.list
ON237248.1

# fetch the sequences:
$ esl-sfetch -f rsv.r500.fa ex1.list > ex1.fa
$ esl-sfetch -f rsv.r500.fa ex2.list > ex2.fa

# run v-annotate.pl on these sequences with 
# the --out_stk option to save the output alignments
$ v-annotate.pl --out_stk --mdir rsv-models --mkey rsv ex1.fa va-ex1
$ v-annotate.pl --out_stk --mdir rsv-models --mkey rsv ex2.fa va-ex2
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
$ grep indf3lcn va-sequences_20231010_2146812/va-sequences_20231010_2146812.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $4, $5); }' | sort | uniq -c | sort -rnk 1
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
$ grep NC_038235 rsv-models/rsv.minfo
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
$ grep mutstart va-sequences_20231010_2146812/va-sequences_20231010_2146812.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
    259 NC_038235 M2-2_protein 8159..8161:+
      1 NC_038235 phosphoprotein 2365..2367:+
      1 NC_001781 matrix_protein 3263..3265:+
```

All but 2 of the 261 `mutstart` alert instances pertain to the
M2-2_protein in `NC_038235`. There were only 286 total RSV A
sequences, this means more than 90% of them do not have the start
codon at positions `8159..8161`. To investigate what is going on here,
let's take a random sample of 10 of these 259 sequences, rerun
`v-annotate.pl` on them and look at their alignment to the `NC_038235`
model. 

```
# pick 10 random sequences with the mutstart alert
$ grep mutstart va-sequences_20231010_2146812/va-sequences_20231010_2146812.vadr.alt | grep 8159..8161 | awk '{ print $2 }' | esl-selectn 10 - > ex3.list
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
$ esl-sfetch -f rsv.r500.fa ex3.list > ex3.fa

# run v-annotate.pl on these sequences with 
# the --out_stk option to save the output alignments
$ v-annotate.pl --out_stk --mdir rsv-models --mkey rsv ex3.fa va-ex3
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
$ grep insertnp va-sequences_20231010_2146812/va-sequences_20231010_2146812.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
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
$ grep insertnp va-sequences_20231010_2146812/va-sequences_20231010_2146812.vadr.alt | grep 5442 | awk '{ print $2 }' | esl-selectn 1 - > ex4.list
$ cat ex4.list
MH760706.1

# fetch the sequence
$ esl-sfetch -f rsv.r500.fa ex4.list > ex4.fa

# run v-annotate.pl on these sequences with 
# the --keep option to save all output files
$ v-annotate.pl --keep --mdir rsv-models --mkey rsv ex4.fa va-ex4
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
$ grep indf3pst va-sequences_20231010_2146812/va-sequences_20231010_2146812.vadr.alt | awk '{ printf ("%s %s\n", $3, $5); }' | sort | uniq -c | sort -rnk 1
    181 NC_038235 attachment_glycoprotein
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
$ esl-sfetch -f rsv.r500.fa ex5.list > ex5.fa

# run v-annotate.pl on these sequences with 
# the --keep option to save all output files
$ v-annotate.pl --keep --mdir rsv-models --mkey rsv ex5.fa va-ex5
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

### Lessons from investigating common alerts in our training set

To summarize what we've learned from manually investigating the top
six most common types of fatal alerts: 

   * The attachment glyoprotein CDS stop codon position in both
     `NC_038235` and `NC_001781` is downstream of the most common stop codon
     position in our RSV A training sequences, by 3 nt (`NC_038235`)
     or 21 nt (`NC_001781`).

   * The `NC_001781` attachment glyoprotein CDS stop codon at
     positions `5587..5589` is 21 nucleotides downstream of the
     most common stop codon position in our RSV B training sequences.

   * There is a duplicated region in many RSV A and B sequences in the
     attachment glycoprotein that is not present in the `NC_038235`
     and `NC_001781` reference sequences. This duplicated region leads
     to several types of VADR alerts including `dupregin`, `insertnp`
     and `indf3pst`.

   * `NC_038235` and `NC_001781` both have `gene` positional
     boundaries that differ from the CDS. Some of the `gene`
     boundaries lead to low alignment confidence related alerts in
     VADR. Typically for viral GenBank submissions based on VADR, the
     gene and CDS boundaries are kept consistent.

   * The `M2-2 protein` CDS in `NC_038235` has a start codon that
     differs from many other RSV A sequences by six positions.

### 5. Potentially choose alternative representative sequences

</details>

---

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

