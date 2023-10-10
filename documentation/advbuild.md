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
`.alt` and `.alt.list` files. To see all the `cdsstopn` (early stop
codon) alerts in the `.alt` file, we can use `grep`. Here are the
first 10 lines that contain `cdsstopn`:

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


We can investigate the 

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

