# <a name="top"></a> Advanced tutorial: building an RSV model library

The `v-build.pl` program will create a model from a single INSDC
accession and include CDS, gene and mature peptide features. However,
a model built from a single accession is often not general enough to
allow most high quality sequences for a viral species to *pass*. For
example, some other sequences may include an extended CDS that has a
different stop codon position from the sequence the model was built
from, and these sequences will *fail* due to fatal alerts related to
the different stop codon. If your goal is too have `v-annotate.pl`
pass the vast majority of sequences that are error-free (lacking
misassemblies, sequencing errors and other artifacts), then you may
want to spend some manual effort One good strategy for building and
refining a model library is:

[Step 1.](#step1) Build one or more models from representative and well annotated
       sequences as a starting point. These may be RefSeq sequences.

[Step 2.](#step2) Construct a training set of randomly chosen existing sequences
        for the viral species.

[Step 3.](#step3) Use `v-annotate.pl` to validate and annotate the training
        sequences using your models from step 1.

[Step 4.](#step4) Analyze the results by looking for common failure modes and
        investigate the sequence characteristics that are responsible
        for them. Based on this analysis, determine if the models from
        Step 1 are sufficient.

[Step 5.](#step5) (Potentially) build new models. If in step 4 there are one or
        more characteristics found that occur in the majority of
        training sequences, it makes sense to pick a new reference
        sequence that includes those characteristics and rebuild the
        models based on those new references. Then you'll want to rerun
        `v-annotate.pl` on the training set using the new models.

[Step 6.](#step6) Analyze results and update models to accomodate existing
        biological sequence and feature diversity.

In this tutorial, we will follow these six steps in building an RSV
library using software installed with VADR, unix command-line
utilities like `grep` and `awk`, as well as the [NCBI Virus
website](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/). 
This tutorial is long and detailed. It can be followed step
by step by rerunning the provided commands locally, or just by reading
through it, or just for reference to specific examples of how to
analyze `v-annotate.pl` results and update VADR models based on those
analyses.

[A discussion of some limitations and alternatives to this
approach](#limit) is included at the end of the tutorial.

---

## <a name="outline"></a>Tutorial outline:

* [Step 1: build model(s) from initial reference sequence(s)](#step1)
  * [Determine good reference sequence(s)](#step1-refseq)
  * [Build initial models](#step1-build)
* [Step 2: construct a training set](#step2)
* [Step 3: run `v-annotate.pl` on training set](#step3)
* [Step 4: analyze results to determine if model is sufficient](#step4)
  * [investigate common `dupregin` alerts](#step4-dupregin)
  * [investigate common `indf3lcn` alerts](#step4-indf3lcn)
  * [investigate common `mutstart` alerts](#step4-mutstart)
  * [investigate common `insertnp` alerts](#step4-insertnp)
  * [investigate common `indf3pst` alerts](#step4-indf3pst)
  * [lessons from investigating common alerts](#step4-lessons)
* [Step 5: (potentially) build new models](#step5)
  * [choose new representative sequences](#step5-chooserep)
  * [build new models from new representative sequences](#step5-build)
  * [rerun `v-annotate.pl` on training set using new models](#step5-rerun)
* [Step 6: analyze results and update models](#step6)
  * [strategies for updating models](#step6-strategies)
  * [adding a protein to a model blastx library](#step6-addblastx)
  * [adding *alternative features*](#step6-alternative)
  * [allowing an *alert exception*](#step6-exception)
  * [rebuilding CM with additional information](#step6-cm)
  * [treating a feature as non-essential](#step6-miscfeat)
  * [making an alert non-fatal](#step6-altpass)
  * [summary of model modifications](#step6-summary)
* [Limitations of and alternatives to this approach](#limit)
  * [reference sequence selection](#limit-ref)
  * [training sequence selection](#limit-train)
  * [replacing one model with multiple models](#limit-multiple)
  * [alignment-based models](#limit-align)
  * [incorporating secondary structure](#limit-secondary)

---

## <a name="step1"></a> Step 1: build model(s) from initial reference sequence(s)

### <a name="step1-refseq"></a> Determine good reference sequence(s) to use

For viruses with multiple subtypes, genotypes, or serotypes, it makes
sense to build at least one model for each. There are two RSV
subtypes, RSV-A and RSV-B, so the first step is to pick
reference sequences for each subtype. 

A good strategy is often to start with a RefSeq sequence if any are
available. In this case there are two RefSeq sequences which can be
found on the [NCBI virus
resource](#https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/), by
selecting the "Search by virus" button and entering "RSV". The top
suggestion will be "Human orthopneumovirus (taxid 11250)", which is a
another name for RSV. At the time of writing, the top two sequences
listed in the resulting list will be RefSeq sequences `NC_038235`
(subgroup A) and `NC_001781` (subgroup B). You can also filter to only
RefSeq sequences using the "Sequence type" filter.

### <a name="step1-build"></a> Build initial models from reference sequence(s)

Next, use `v-build.pl` to build the two models, specifying the
`--group` and `--subgroup` options as below:

```
$ v-build.pl --group RSV --subgroup A NC_038235 NC_038235
$ v-build.pl --group RSV --subgroup B NC_001781 NC_001781
```

These commands will take a long time, up to one hour each. 

When they are finished combine the two models into a model library by
following the steps below ([also listed here](build.md#library)).

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
$ $VADREASELDIR/esl-sfetch --index rsv-models1/rsv.fa
$ $VADRINFERNALDIR/cmpress rsv-models1/rsv.cm
$ $VADRHMMERDIR/hmmpress rsv-models1/rsv.hmm
$ $VADRBLASTDIR/makeblastdb -dbtype nucl -in rsv-models1/rsv.fa
```

At this point, it's a good idea to test the models by using them to
annotate the reference sequences they were built from. This will
check that we combined the models correctly. The sequences
should pass, although they are not guaranteed to do so (some reference
sequences may have special characteristics that cause them to fail
even using models built from themselves). For these RSV models,
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
---

## <a name="step2"></a> Step 2: construct a training set

When evaluating a VADR model it is critical to examine at how it performs
when used with `v-annotate.pl` on example sequences. Ideally, you
would know what the expected result is (pass or fail status, and
specific alerts you expect to be or not be reported) for each of the
example sequences.

For instance, if you have a set of high quality sequences that have
been expertly validated and annotated, you could use this as a
positive control `v-annotate.pl` should pass all of those sequences
and give annotation matching what is expected.

Alternatively, you could intentionally introduce sequencing errors
(insertions, deletions, rearrangements) into high quality sequences,
and check to see if `v-annotate.pl` detects those problems and
reports them.

Often times, however, the most readily available set of sequences is
simply INSDC sequences of the viral species you are modelling. Many of
these will be of high quality without any sequencing errors,
misassemblies or other artifacts, but some may have those types of errors.

For this tutorial, we will take the strategy of selecting a random subset
of nearly full length sequences, evaluating them with `v-annotate.pl`
and manually analyzing the results as a way of evaluating our models.

<a name="step2-length"></a>(The decision to use only full length sequences is debatable, as by
doing it we are assuming that the sequence diversity represented by all
sequences, including partial sequences, is well represented by only
the full length sequences. In other words, if we optimize the
models for only full length sequences, they may perform poorly on
existing partial length sequences that are sufficiently divergent from
all full length sequences. Two alternatives you might consider are:
allowing partial sequences in your training set (although this
somewhat complicates the analysis of the results), and optimizing
models on full length sequences first, and then testing on partial
sequences to see if any new failure modes occur.)

One way to select a random set of training sequences is to use the
*NCBI Virus* resource. I usually define nearly full length as 95% the
length of the RefSeq sequence, or greater. In this case `NC_038235` is
15222 nucleotides and `NC_001781` is 15225 nucleotides, so a 95%
length cutoff is about 14460 nucleotides.  At the time of writing
there are about 5500 nearly full length RSV INSDC sequences by this
definition. For our random subset it is important to select enough
sequences that you should identify any common failure modes, but not
too many that your manual analysis will take a prohibitively long
amount of time (VADR is slow). For RSV, I chose to use 500 randomly
chosen nearly full length sequences for model improvement.

To download 500 randomly chosen RSV sequences of length 14460 or
greater: 

1. go to the [NCBI Virus RSV list page](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Human%20orthopneumovirus,%20taxid:11250) you reached in step 1

2. use the "Sequence Length" filter to set a minimum of 14460 and then
click the "Download" button near the top left of the page

3. select "Nucleotide" under "Sequence data (FASTA format)"

4. select  "Download a randomized subset of all records",
enter "500", click "Next" and "Use Default" for the FASTA
definition line, then finally click "Download". 

This should download a file called something like `sequences_20231010_2146812.fasta`.
For convenience, rename the downloaded fasta file as `rsv.r500.fa`.
(You can find the accession list for the 500 randomly selected
sequences used in this tutorial in [`vadr/documentation/build-files/rsv.r500.list`](build-files/rsv.r500.list).)

---

## <a name="step3"></a> Step 3. Run `v-annotate.pl` to validate and annotate sequences in training set

Next, we'll use our new models to annotate our set of 500 training
sequences. The RSV genome is about 15Kb long, which is towards the
longer end of genome size that VADR is capable of annotating. The
memory and speed requirements for `v-annotate.pl` don't scale well,
and annotation of RSV can require up to 64Gb of RAM and take about 1
minute per full length sequence.  (VADR is used for GenBank screening
of SARS-CoV-2 (30Kb genome) sequences with significantly lower memory
requirements, but it utilizes heuristics that work particularly well
for SARS-CoV-2 and I don't recommend using those heuristics with RSV
sequences.)

To use v-annotate.pl on our set of 500 sequences, we would execute:

```
$ v-annotate.pl --mdir rsv-models1 --mkey rsv rsv.r500.fa va-rsv.r500
```

This will take about 1 minute per sequence, so roughly
8 hours to complete. You can parallelize it if you have a lot of RAM
and multiple CPUs using the `--split` and `--cpu` options
as described [here](annotate.md#options-split).

---

---

## <a name="step4"></a> Step 4: analyze the results and update models

Next, we want to analyze the results. From the `v-annotate.pl` output (also saved to
the `va-r500/va-r500.vadr.log` file), we can see that 286 sequences were
classified as RSV A (matched best to the `NC_038235` model) and 214 as
RSV B (matched best to the `NC_001781` model). And only 6 out of the 500
total training sequences pass:

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
behavior and so are okay as they are. The fact that so many
sequences fail seems to indicate that the models should be modified,
but it could be that RSV viral sequences are so highly variable that a
high failure rate is expected, and no single sequence based model
would allow significantly more sequences to pass. The only way to know
for sure is to drill down deeper into the results.

To investigate we need to look in further detail at the reasons
sequences are failing. A good way to do this is to go through the most
commonly reported fatal alerts.  Of the 26 fatal alerts (`yes` in
`causes failure` column, numbers 15 to
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
most common alerts in detail next. Documentation on the alerts can be
found [here](annotate.md#alerts) with additional examples [here](alerts.md#top).

### Investigate common `cdsstopn` alerts

To see all the `cdsstopn` (early stop
codon) alerts in the `.alt` file, we can use `grep` and `head`. Here are the
first 10 lines that contain `cdsstopn`. (The first `head` command is
used only to display the column headings.):


```
$ head -n 3 va-rsv.r500/va-rsv.r500.vadr.alt
#       seq                    ftr   ftr                      ftr  alert           alert                                            seq    seq                       mdl    mdl  alert 
#idx    name        model      type  name                     idx  code      fail  description                                   coords    len                    coords    len  detail
#-----  ----------  ---------  ----  -----------------------  ---  --------  ----  ---------------------------  -----------------------  -----  ------------------------  -----  ------

$ grep cdsstopn va-rsv.r500/va-rsv.r500.vadr.alt | head
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
`uniq` unix command line utilities:

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
`5579..5581:+`. Because this is causing a `cdsstopn` alert, it must
differ from the the model reference stop codon position for the `attachment
glycoprotein` CDS. We can be figure out what that is by looking at the 
`rsv.minfo` file:

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
different stop positions than the reference is strong indication that maybe
our reference models should be rebuilt based on different sequences
that include the most common attachment glycoprotein stop position. But
it makes sense at this point to look at all of the common failure
modes before making that decision and looking for new
references. We may find additional attributes that we want our new
references to have. 

---

### <a name="step4-dupregin"></a> Investigate common `dupregin` alerts

The second most common fatal alert is the `dupregin` alert that occurs
when "similarity to a model region occurs more than once". 

```
16    dupregin  yes      DUPLICATE_REGIONS              sequence    344   344  similarity to a model region occurs more than once
```

We can use `grep` to look at some examples of the `dupregin` alert:
```
$ head -n 3 va-rsv.r500/va-rsv.r500.vadr.alt
#       seq                    ftr   ftr                      ftr  alert           alert                                            seq    seq                       mdl    mdl  alert 
#idx    name        model      type  name                     idx  code      fail  description                                   coords    len                    coords    len  detail
#-----  ----------  ---------  ----  -----------------------  ---  --------  ----  ---------------------------  -----------------------  -----  ------------------------  -----  ------

$ grep dupregin va-rsv.r500/va-rsv.r500.vadr.alt | head
1.1.1   OR143220.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5502..15225:+,2..5501:+  15224  5466..15199:+,31..5537:+  15241  similarity to a model region occurs more than once [5466..5537:+ (len 72>=20) hits 1 (8569.8 bits) and 2 (4779.2 bits)]
4.1.1   OR287871.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5438..14973:+,1..5437:+  14973  5466..15003:+,93..5537:+  14983  similarity to a model region occurs more than once [5466..5537:+ (len 72>=20) hits 1 (8198.7 bits) and 2 (4724.2 bits)]
5.1.1   OM857265.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5434..14961:+,1..5433:+  14961  5466..14995:+,99..5537:+  14969  similarity to a model region occurs more than once [5466..5537:+ (len 72>=20) hits 1 (8465.0 bits) and 2 (4737.7 bits)]
8.1.1   OR143199.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5502..15222:+,2..5501:+  15221  5466..15199:+,31..5537:+  15241  similarity to a model region occurs more than once [5466..5537:+ (len 72>=20) hits 1 (8590.1 bits) and 2 (4779.3 bits)]
9.1.1   MG431251.1  NC_001781  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5446..15245:+,1..5445:+  15245  5410..15210:+,17..5469:+  15254  similarity to a model region occurs more than once [5410..5469:+ (len 60>=20) hits 1 (9138.1 bits) and 2 (4959.7 bits)]
10.1.1  KY249668.1  NC_001781  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5455..15267:+,1..5454:+  15267   5410..15223:+,1..5469:+  15283  similarity to a model region occurs more than once [5410..5469:+ (len 60>=20) hits 1 (9118.2 bits) and 2 (4945.9 bits)]
11.1.1  OR287899.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5478..15168:+,1..5477:+  15168  5465..15165:+,54..5536:+  15184  similarity to a model region occurs more than once [5465..5536:+ (len 72>=20) hits 1 (8586.5 bits) and 2 (4753.8 bits)]
13.1.1  OR326741.1  NC_001781  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5413..15191:+,2..5412:+  15190  5410..15189:+,50..5469:+  15200  similarity to a model region occurs more than once [5410..5469:+ (len 60>=20) hits 1 (9027.0 bits) and 2 (4873.7 bits)]
15.1.1  OR143187.1  NC_038235  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5502..15224:+,2..5501:+  15223  5466..15199:+,31..5537:+  15241  similarity to a model region occurs more than once [5466..5537:+ (len 72>=20) hits 1 (8573.2 bits) and 2 (4797.4 bits)]
17.1.1  OR326763.1  NC_001781  -     -                          -  dupregin  yes   DUPLICATE_REGIONS            5412..15191:+,2..5411:+  15190  5409..15189:+,50..5468:+  15200  similarity to a model region occurs more than once [5409..5468:+ (len 60>=20) hits 1 (8998.9 bits) and 2 (4888.7 bits)]
```

The `alert detail` at the end of each line are very similar for the
six of the first ten instances of `dupregin`: `similarity to a model
region occurs more than once [5466..5537:+ (len 72>=20)`. This
suggests these alerts are all referring to the same situation. If a
genome has repetitive regions you may see this alert for all
sequences, but in that case you will likely see it for the model
sequence as well. (If that occurs, refer to the section on [alert
exceptions below](#step6-exception). Even though the alert details are
similar, note that the model positions in field 12 are not
identical. For example, in the first alert, the model coordinates are
`5466..15199:+,31..5537:+`, and in the second the coordinates are
`5466..15003:+,93..5537:+`. These coordinates are showing the
reference coordinates of the two hits that overlap, whereas the alert
detail shows the actual region of overlap. So for this alert, in order
to group together similar situations we use `awk` to select different
fields than in the above `cdsstopn` alert. Here we are more interested
in the 23rd field (e.g. `[5466..5537:+`), then the 12th.

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
one of two duplicated regions. To gain a better understanding of this
duplicated region, we can select an example of one RSV A sequence and
one RSV B sequence that contain this alert and run those through
`v-annotate.pl` again. To select two sequences we'll use the
`esl-sfetch` program that is installed with VADR in the directory
pointed to by your `$VADREASELDIR` environment variable:

```
# prepare the sequence file for the esl-sfetch program 
$ $VADREASELDIR/esl-sfetch --index rsv.r500.fa

# pick a sequence with a dupregin alert that contains the strings we are interested it
$ grep dupregin va-rsv.r500/va-rsv.r500.vadr.alt | grep 5410..5469 | awk '{ print $2 }' | $VADREASELDIR/esl-selectn 1 - > ex1.list
$ grep dupregin va-rsv.r500/va-rsv.r500.vadr.alt | grep 5466..5537 | awk '{ print $2 }' | $VADREASELDIR/esl-selectn 1 - > ex2.list
$ cat ex1.list 
MZ516003.1
$ cat ex2.list
ON237248.1

# fetch the sequences:
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex1.list > ex1.fa
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex2.list > ex2.fa

# run v-annotate.pl on these sequences with the --out_stk option to save the output alignments
$ v-annotate.pl --out_stk --mdir rsv-models1 --mkey rsv ex1.fa va-ex1
$ v-annotate.pl --out_stk --mdir rsv-models1 --mkey rsv ex2.fa va-ex2
```

Let's look at the RSV B sequence first. The
`va-ex1/va-ex1.vadr.NC_001781.align.stk` file contains the `MZ516003.1`
sequence aligned to the `NC_001781` model in [Stockholm alignment file
format](https://en.wikipedia.org/wiki/Stockholm_format), with
reference positions numbered using rows with the header `#=GC
RFCOL`. Below is an excerpt of that alignment file with the duplicated
region annotated in an extra line that I've added labelled
`dupregin`. The positions marked with `1` show the first instance of
the repeated region, and those marked with `2` show the second
instance. Note that these positions correspond to positions
`5410..5469` of the reference model:

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

dupregin                    111111111111111111111111111111111111111111111111111111111111222222222222222222222222222222222222222222222222222222222222
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
confidence at each position as explained more [here](alerts.md#pp)).

This region falls within the `attachment glycoprotein` CDS. To learn
more about this duplication we might search
[PubMed](https://pubmed.ncbi.nlm.nih.gov/) with the query "RSV
attachment glycoprotein duplicated region", which returns several
articles related to this region, including a paper entitled
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
`NC_001781` and `NC_038235` reference sequences. This is another
indication that we should change our reference sequences to sequences that
include the more common stop position of the attachment glycoprotein,
and that also include this duplicated region. 

There may still be more characteristics that we want to include, so we should
continue to investigate the other common alerts from above. 

---

### <a name="step4-indf3lcn"></a> Investigate common `indf3lcn` alerts

The third most common alert was `indf3lcn`:

```
32    indf3lcn  yes      INDEFINITE_ANNOTATION_END       feature   1057   320  alignment to homology model has low confidence at 3' boundary for feature that does not match a CDS
```

This alert occurs because the alignment confidence at the end
coordinate/position  of a feature is relatively low, indicating that it may
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
boundaries. This is consistent with the norovirus, dengue, and
SARS-CoV-2 VADR models used by GenBank, which all have identical start
and end points for corresponding CDS and gene features. If for your
own models you desire the distinct boundaries then you can, of course,
keep them. In that case you may want to possibly define `indf3lcn` alerts
as non-fatal using the `--alt_pass indf3lcn` option to
`v-annotate.pl`.

When we choose our new reference sequences in the next section we will
revisit this issue of differing gene and CDS boundaries.

---

### <a name="step4-mutstart"></a> Investigate common `mutstart` alerts

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
M2-2_protein in `NC_038235`. There are only 286 total RSV A
sequences, so this means that more than 90% of them do not have the start
codon at positions `8159..8161`. To investigate this further, 
let's take a random sample of 10 of these 259 sequences, rerun
`v-annotate.pl` on them and look at their alignment to the `NC_038235`
model. 

```
# pick 10 random sequences with the mutstart alert
$ grep mutstart va-rsv.r500/va-rsv.r500.vadr.alt | grep 8159..8161 | awk '{ print $2 }' | $VADREASELDIR/esl-selectn 10 - > ex3.list
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

<a name="step4-mutstart-aln"></a> Below is an excerpt of the resulting alignment in
`va-ex3/va-ex3.vadr_NC_038235.align.stk`

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
positions. The three nucleotides labelled with `222` that occur six
nucleotides downstream at positions `8165..8167` are 
`ATG`, and are in-frame with the reference with positions
`8159..8161`. Because the majority of RSV A sequences in our training
set have the second `ATG` but not the first, we may want our model to
use that as the start position. If we do that however, then the
annotated start for the `NC_038325` model will be annotated
differently. There is a way to deal with this and allow
`v-annotate.pl` to pick either start position. We'll revisit this
[below](#step6-m2start). 

---

### <a name="step4-insertnp"></a> Investigate common `insertnp` alerts

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

About 75% of these alerts are for large insertions in the blastx
protein alignment in RSV B (`NC_001781`) sequences at position
5442. This position is within the attachment glycoprotein region for
which dupregin alerts were also reported: `5410..5469`. It makes sense
that a large duplication could result in a large insertion. In fact,
it is somewhat surprising that there aren't similar alerts for the
`NC_038235` model, although we will find out why soon enough. Dealing
with the `NC_00178`1 duplication in the next round of model building
will likely eliminate these alerts.

It is possible to look at these insertions in the actual blastx output
alignments, but you'll need to rerun `v-annotate.pl` using the
`--keep` option, which makes it save all intermediate files that are
usually deleted. To do that we would sample one sequences and
rerun it (you can sample more than one, but often it is easier to
find the relevant part of the blastx output if you have restricted
your input to a single sequence):

```
# pick a sequence with the insertnp alert:
$ grep insertnp va-rsv.r500/va-rsv.r500.vadr.alt | grep 5442 | awk '{ print $2 }' | $VADREASELDIR/esl-selectn 1 - > ex4.list
$ cat ex4.list
MH760706.1

# fetch the sequence
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex4.list > ex4.fa

# run v-annotate.pl on these sequences with the --keep option to save all output files
$ v-annotate.pl --keep --mdir rsv-models1 --mkey rsv ex4.fa va-ex4
```

The relevant blastx output will be in the file
`va-ex4/va-ex4.vadr.NC_001781.blastx.out`. The section of the file you
are looking for is the results for when the predicted attachment glycoprotein
CDS is used as a query sequence. To find this we need to know what the
coordinates are for that prediction. We can find these in the `.ftr`
output file (format described [here](formats.md#ftr) or in the `.tbl`
file `va-ex4/va-ex4.vadr.fail.tbl`: 

```
4590	5543	misc_feature
			note	similar to attachment glycoprotein
```

If we search for `4590..5543` in the blastx output file, we will
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

### <a name="step4-indf3pst"></a>Investigate common `indf3pst` alerts

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
validation stage does not extend close enough to the 3' end. As this
is again in the attachment glycoprotein CDS it may be related to the
duplicated region. Note that the vast majority of instances are for
RSV A (`NC_038235`) sequences. A possible explanation is that these
are predominantly made up of sequences that also have a `dupregin`
alert and blastx does not create an alignment that spans the
duplicated region as it did for the example RSV B sequence MH760706.1
with the `insertnp` alert above. The duplicated region was typically
the same size in both RSV A and RSV B sequences, so the difference
between VADR reporting an `insertnp` or `indf3pst` alert is probably
due to the similarity in the 3' ends of the protein between the
training sequences and the references: for RSV B, the 3' ends are
similar enough that the best scoring blastx alignment extends across
the duplication, whereas for RSV A, the 3' end is not similar enough
and the blastx alignment stops before the duplication. To check if
this is indeed what's happening we can again look at the blastx output
after rerunning an example sequence:

```
# pick a sequence with the indf3pst alert:
$ grep indf3pst va-rsv.r500/va-rsv.r500.vadr.alt | grep attachment | grep NC_038235 | awk '{ print $2 }' | $VADREASELDIR/esl-selectn 1 - > ex5.list
$ cat ex5.list
MH181932.1

# fetch the sequence
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex5.list > ex5.fa

# run v-annotate.pl on these sequences with the --keep option to save all output files
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
because `5584-4688+1=969`. This is why the length of the relevant
sequence region reported in the `alert detail` in the `alt` file is
`120` (`969-849=120`):

```
$ grep indf3pst va-ex5.vadr.alt
1.5.3  MH181932.1  NC_038235  CDS   attachment_glycoprotein   14  indf3pst  yes   INDEFINITE_ANNOTATION_END             5486..5605:+    120              5584..5584:+      1  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint [120>8, valid stop codon in nucleotide-based prediction]
```
--- 

### <a name="step4-lessons"></a>Lessons from investigating common alerts

We've learned that both `NC_038235` and `NC_001781` are lacking
several important characteristics that are present in the
majority of RSV A and RSV B sequences. These include:

 * attachment glycoprotein stop codon at reference coordinates
   `5579..5581` in RSV A and `5566..5568` in RSV B.

 * duplicated region in attachment glycoprotein near reference
   positions `5466..5537` in RSV A and `5410..5469` in RSV B.

 * M2-2 protein start codon at reference positions `8165..8167` in RSV
   A. 

And further we've observed that:

  * `NC_038235` and `NC_001781` both have gene positional boundaries
    that differ from the corresponding CDS boundaries. `v-annotate.pl`
    commonly reports low alignment confidence related alerts
    (e.g. `indf3lcn`) for the gene boundaries. Typically for viral
    GenBank submissions based on VADR, the gene and CDS boundaries are
    kept consistent.

At this point, because we've found, for both models, at least one
characteristic that causes fatal alert(s) and is present in the
majority of the sequences but lacking in the reference sequence, we
will identify new, more representative sequences from which to build a
new set of models. 

---

### <a name="step5"></a> Step 5. (Potentially) choose new representative sequences and build new models

It is not always necessary to build new models at this point. If all
of the failure modes above had been in a minority of sequences, then
it may have made more sense to stick with our RefSeq-based models and
simply modified them (as described below) to address the failure modes
we wanted to eliminate. But in this case, as explained above, we want
to choose new representatives and build new models.

#### <a name="step5-chooserep"></a> Identifying new representative sequences

We will choose a new representative from our random set of 500
training sequences. Of course, it is possible, and probably even
likely, that there exists a 'better' representative sequence that is
not in our training set, but using the most representative sequence
from our training set of 500 is probably good enough. First we need to
identify the subset of our 500 sequences that contain the three major
characteristics we've already identified, summarized
[above](#step4-lessons). We'll do this separately for RSV A and RSV B:

```
# fetch out the sequences with the early stop codon (cdsstopn) at 5579..5581:
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep cdsstopn | grep 5579..5581 | awk '{ print $2 }' | wc -l
213
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep cdsstopn | grep 5579..5581 | awk '{ print $2 }' | sort > 213.list 

# fetch out the sequences with the dupregin at 5466..5537:
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep dupregin | grep 5466..5537 | awk '{ print $2 }' | wc -l
159
$ grep NC_038235 va-rsv.r500/va-rsv.r500.vadr.alt | grep dupregin | grep 5466..5537 | awk '{ print $2 }' | sort > 159.list

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

# fetch the 140 sequences into a new fasta file:
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa 140.list > rsvA.140.fa
```
And then repeat the same for RSV B, skipping the M2-2 start codon step
which doesn't apply to RSV B:
```
# fetch out the sequences with the early stop codon (cdsstopn) at 5566..5568:
$ grep NC_001781 va-rsv.r500/va-rsv.r500.vadr.alt | grep cdsstopn | grep 5566..5568 | awk '{ print $2 }' | wc -l
154
$ grep NC_001781 va-rsv.r500/va-rsv.r500.vadr.alt | grep cdsstopn | grep 5566..5568 | awk '{ print $2 }' | sort > 154.list

# fetch out the sequences with the dupregin at 5410..5469
$ grep NC_001781 va-rsv.r500/va-rsv.r500.vadr.alt | grep dupregin | grep 5410..5469 | awk '{ print $2 }' | wc -l
167
$ grep NC_001781 va-rsv.r500/va-rsv.r500.vadr.alt | grep dupregin | grep 5410..5469 | awk '{ print $2 }' | sort > 167.list

# use the unix 'comm' command to get the subset of seqs common to both lists
$ comm -1 -2 154.list 167.list | wc -l
129
$ comm -1 -2 154.list 167.list > 129.list

# fetch the 129 sequences into a new fasta file:
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa 129.list > rsvB.129.fa
```

It is important that the reference sequences we choose do not have too
many ambiguous nucleotides because they make the models less
specific. (Indeed we may want our models to be less specific and more
general in some areas of the genome that are less highly conserved,
but including ambiguous nucleotides from a single reference sequence
is not the best way to do this. We'll revisit this topic briefly at
[the end of the tutorial](#limit-align).)

Because we have many candidates, we may be able to afford to removing
all sequences with 1 or more ambiguous nucleotides. We can determine
how many ambiguous characters are in each sequence using the
`count-ambigs.pl` *miniscript* included with VADR:

```
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl rsvA.140.fa | head
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
```

The second field is the number of ambiguous nucleotides in each
sequence. We can fetch out all lines that have 1 or more in this field
with: 

```
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl rsvA.140.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep -v " 0" | head 
LC530050.1 1
MG813982.1 46
MN078121.1 7
MN535098.1 217
MN536995.1 1238
MN536996.1 65
MN536997.1 103
MN536999.1 150
MZ515551.1 98
MZ515654.1 2780
```

And so to only save the sequences with 0 ambiguous nucleotides:
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

One way to pick a representative is to pick the sequence with a high
average percent identity to all other candidates based on an
alignment. We can use `v-annotate.pl` to generate multiple alignments
of all candidates:

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
candidates. The final criterion is the length - we don't want the
representative sequence to be too short. In this case, we know that the
`NC_038235` and `NC_001781` models are considered "full length" as
indicated in their GenBank annotation:

[`NC_038235`](https://www.ncbi.nlm.nih.gov/nuccore/NC_038235)
```
COMMENT     REVIEWED REFSEQ: This record has been curated by NCBI staff. The
            reference sequence is identical to M74568.
            COMPLETENESS: full length.
```

[`NC_001781`](https://www.ncbi.nlm.nih.gov/nuccore/NC_001781)
```
COMMENT     REVIEWED REFSEQ: This record has been curated by NCBI staff. The
            reference sequence is identical to AF013254.
            COMPLETENESS: full length.

```

So our new representatives should align end to end with their
respective model RefSeqs. We can determine those that do by inspecting
the alignments. To make viewing the alignments easier, we can pull out
only the top 10 candidates for each subtype using the `esl-alimanip` program
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

So `KY654518.1` will be our new RSV A reference sequence. We can
repeat the drill for the `rsvB.top10.stk` alignment:

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

The final step is to see if `KY654518` and `MZ516105` satisfy the
criteria that the corresponding `gene` and `CDS` coordinates are
identical, unlike in the RefSeqs. Inspecting the GenBank record pages for
[`KY654518`](https://www.ncbi.nlm.nih.gov/nuccore/KY654518) and 
[`MZ516105`]((https://www.ncbi.nlm.nih.gov/nuccore/MZ516105) allows us
to confirm that they do. (If they did not, we could manually modify the `gene`
feature boundaries in the model info file after running `v-build.pl`
so that they did.)

#### <a name="step5-build"></a> Building new models from our new representative sequences

To build our new models, we run `v-build.pl`:

```
$ v-build.pl --group RSV --subgroup A KY654518 KY654518
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

As in step 1, it's a good idea to run `v-annotate.pl` with the new models against the two model
sequences as a sanity check. For these two RSV models, both sequences
should pass:

```
$ v-annotate.pl --out_stk --mdir rsv-models2 --mkey rsv rsv-models2/rsv.fa va-rsv2
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

#### <a name="step5-rerun"></a> Rerun `v-annotate.pl` on our existing training set using new models

Next we want to evaluate the performance of our new models. 
We can repeat the `v-annotate.pl` command from step 2 using our
new models, but this time we will use the `--out_stk` option to save
multiple alignments for reasons that will become clear below:

```
$ v-annotate.pl --out_stk --mdir rsv-models2 --mkey rsv rsv.r500fa va2-r500
```

---

### <a name="step6"></a> Step 6: Analyze the results and update the models accordingly

This time, from the `v-annotate.pl` output we can tell that many more
sequences passed than with the original models, when only 6 out of 500 passed:

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

But there are still 233 sequences that do not pass. As we did
previously, we can walk through the most common types of alerts to
investigate the reasons for the failures. But this time, we will
modify our models so that they do not report alerts due to valid
biological diversity that we want `v-annotate.pl` to allow. It makes
sense to start modifying the model at this stage (as opposed to
building yet another model) because we know we are using good
representative sequences for our models, based on the work we did
analyzing the results of the initial RefSeq models.

### <a name="step6-strategies"></a> Strategies for modifying models

There are several ways we can modify or update our models to avoid reporting
alerts for expected biological characteristics, including:

1. **Add a new protein to the blastx protein library for a model to
   account for protein sequence variability.** This
   can help remove unnecessary alerts related to the protein validation stage, most
   commonly: `indfpst5`, `indfpst3`, `insertnp` and
   `deletinp`. [There are two examples of this strategy detailed below:
   [1](#step6-addblastx) and [2](#step6-alternative-blastx).

2. **Add **alternative features** to the model info file, along with
   corresponding proteins to the blastx protein library.** Some CDS
   may have multiple start and stop positions, with respect to the
   reference model.  We can deal with this by adding all of the
   acceptable alternatives as separate features to the model info
   file. `v-annotate.pl` will attempt to annotate each of the
   alternatives and report annotation for the single alternative that
   yields the fewest fatal alerts. There is an example of this
   strategy detailed [below](#step6-alternative).

3. **Add an alert *exception* to the model info file.** For some alert
   types, exceptions can be added that prevent the reporting of alerts
   in specific model regions. For example, an exception can be added
   to prevent the reporting of a `dupregin` alert due to an expected
   repetitive region.

4. **Rebuild the covariance model from a new input alignment.** We can
   still use the same reference sequence and positions, but by
   buildng a new model from an alignment with multiple sequences,
   certain alert instances caused by sequence differences with the
   reference sequence can be prevented.

5. **Specify features as non-essential**. By adding a a
   `misc_not_failure` flag to a feature in the model info file, we can
   specify that many types of fatal alerts *not* cause a sequence
   to fail, but instead cause the associated feature to be annotated
   as a `misc_feature` in the output. 

6. **Specify command-line options when running `v-annotate.pl` to make
   some alerts fatal or non-fatal**. For some viruses, certain fatal alerts
   may be so common that they are expected for many sequences, yet we
   don't want them to cause a sequence to fail. We can make those
   alerts non-fatal using the command-line option
   `alt_pass`. (Technically, this isn't a modification to the models,
   but rather a change in how `v-annotate.pl` is run.)

We will now walk through examples for the first four strategies, and
elaborate on five and six below.

First, we need to identify the characteristics that we would
like to address through model modification. Once again we can start by
looking at the most common types of reported alerts:

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
than 10 sequences (more than 2\% of the sequences):

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

### <a name="step6-addblastx"></a> Adding a protein to model blastx library

Let's examine the `deletinp` alerts. As above, we can sort
all the occurences of this alert in the `.alt` file and group them:

```
$ grep deletinp va3-rsv.r500/va3-rsv.r500.vadr.alt | awk '{ printf ("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
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
region that is close to the duplicate region [we observed in our
RefSeq model results](#step4-dupregin). We attempted to address those
`dupregin` alerts by rebuilding our models with new sequences that
*included* the duplicated region, but now it seems that we are
observing that sequences that do not have the duplication are failing with
a `deletinp` alert. These are a minority of the sequences but still a
significant number in each model.

The `deletinp` alert occurs when there is a region in the best
blastx alignment that includes a deletion that is longer than 9
amino acids (27 nucleotides). Let's take a look at one example of the
most common deletion span: `5461..5532:+` to the `KY654518` (RSV A)
model, by randoming selecting one sequence and rerunning
`v-annotate.pl` on it with the `--keep` option as we did earlier:

```
# pick a sequence with the deletinp alert:
$ grep deletinp va2-rsv.r500/va2-rsv.r500.vadr.alt | grep 5461..5532 | awk '{ print $2 }' | $VADREASELDIR/esl-selectn 1 - > ex6.list
$ cat ex6.list
KX655676.1

# fetch the sequence
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex6.list > ex6.fa

# run v-annotate.pl on these sequences with --keep option to save all output files
$ v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex6.fa va-ex6
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

The relevant file is the blastx output file
`va-ex6/va-ex6.vadr.KY654518.blastx.out`, and we are interested in the
blastx alignment for the attachment glycoprotein, which is positions
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
then be used as additional subjects in the blastx-based protein validation
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
output from `v-annotate.pl` (which was only output because we used the
`--out_stk` option). We are interested in a representative example of
the CDS with the `deletinp` alert for the `KY654318` model, so it
makes sense to use a sequence with the most common deletion of
reference positions `5461..5532:+`. To get a list of qualifying
sequences and extract the relevant alignment subset:

```
# get a list of the 69 candidates
$ cat va2-rsv.r500/va2-rsv.r500.vadr.alt | grep deletinp | grep KY654518 | grep attachment_glycoprotein | grep 5461..5532 | awk '{ printf("%s\n", $2); }' > ex7.list

# extract the 69 aligned sequences from the alignment:
$ $VADREASELDIR/esl-alimanip --seq-k ex7.list va2-rsv.r500/va2-rsv.r500.vadr.KY654518.align.stk > ex7.stk 
```

Now, since we are only interested in finding a representative
attachment glycoprotein CDS sequence (as opposed to the full
sequence), we can restrict the sequences to only that CDS using the
`esl-alimask` program that is installed with VADR:

```
# extract the attachment glycoprotein region:
$ $VADREASELDIR/esl-alimask -t --t-rf ex7.stk 4681..5646 > ex7.ag.stk
```

Next, as we did when looking for new reference sequences above, we can
take the following steps to identify a good representative attachment
glycoprotein CDS sequence:

1. Remove all sequences with any ambiguous nucleotides using the
`count-ambigs.pl` script.

2. use `esl-alipid` and the `esl-alipid-per-seq-stats.pl` script
   to find the candidate sequence that has the highest average
   percent identity to all other candidate sequences.

```
# convert to fasta 
# $ $VADREASELDIR/esl-reformat fasta ex7.ag.stk > ex7.ag.fa

# remove any sequences with ambiguous nucleotides
# $ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl ex7.ag.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | awk '{ printf("%s\n", $1); }' | wc -l
60
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl ex7.ag.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | awk '{ printf("%s\n", $1); }' > ex7.60.list
# $VADREASELDIR/esl-alimanip --seq-k ex7.60.list ex7.ag.stk > ex7.ag.60.stk

# determine pairwise percent identity
$ $VADREASELDIR/esl-alipid ex7.ag.60.stk > ex7.ag.60.alipid

# calculate average percent identity
$ perl $VADRSCRIPTSDIR/miniscripts/esl-alipid-per-seq-stats.pl ex7.ag.60.alipid > ex7.ag.60.alipid.perseq

# list top candidates
$ grep -v ^\# ex7.ag.60.alipid.perseq | sort -rnk 2 | head
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

We will use the `KX510193.1` attachment glycoprotein CDS sequence,
which has the highest average nucleotide percent identity with all
other candidates. (We could have used a protein alignment as the basis
for selecting a representative in this step, especially since we are
trying to optimize blastx alignments, but the nucleotide-based
approach we've taken here should be a good enough proxy for the
protein-based approach.)

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

So all of the relevant values for the `build-add-to-blast-db.pl`
script are:

| command-line argument    | value |
|--------------------------|-------------------------|
| `<path to .minfo file>`  | `rsv-models2/rsv.minfo` | 
| `<path to blast db dir>` | `rsv-models2/`          |
| `<model name>`           | `KY654518`              | 
| `<nt-accn-to-add>`       | `KX510193`              |
| `nt-coords-to-add>`      | `4409..5302:+`          |
| `<model-CDS-feature-coords>` | `4681..5646:+`       | 
| `<name for output directory>` | `vb-ex7`             | 

To add the protein to the blastx library:
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
$ $VADREASELDIR/esl-seqstat -a rsv-models2/KY654518.vadr.protein.fa
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
added. Both of these sequences are now possible blastx subject
sequences for the attachment glycoprotein CDS.

We can then perform a sanity check to make sure that this added
sequence has the intended effect. Let's run `v-annotate.pl` on our
`ex6.fa` sequence. It should now pass. 

```
$ v-annotate.pl -f --keep --mdir rsv-models2 --mkey rsv ex6.fa va-ex6
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
probably next for the 35 `MZ516105` alerts. To do that, we would
repeat the above procedure to find and add a suitable representative sequence
to add to the protein blast library. After that, we might want to 
to rerun all of the sequences that failed due to `deletinp` alerts
with the updated models. If a significant number of sequences still
fail due to `deletinp` alerts at that stage, then we could repeat the
process again. For the purposes of this tutorial, we will move on to the next most
common alert `indf3pst` to provide a different example of
updating a model.

### <a name="step6-alternative"></a>Adding an alternative CDS feature with different stop coordinate

Let's take a look at one example of the `indf3pst` alert:

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
field is space-delimited field 26, which is part of the `alert detail`
column. Field 26 for `indf3pst` alerts lists how far the protein-based endpoint
and nucleotide based endpoint is. In this case it is 6 nucleotides
which exceeds the maximum allowed without an alert: `6>5`.

When grouping instances of this alert we should output this field 26 value as
well because any instances of this alert for the same feature *and* the same distance may be
able to be addressed with the same model modification:

```
$ grep indf3pst va2-rsv.r500/va2-rsv.r500.vadr.alt | awk '{ printf ("%s %s %s %s\n", $3, $5, $12, $26); }' | sort | uniq -c | sort -rnk 1 | head
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
```

It turns out that the example we looked at above in attachment
glycoprotein for model `MZ516105` with a difference of 6 nucleotides
is the most common one. Let's investigate one of the sequences with
that particular alert further:

```
$ grep indf3pst va2-rsv.r500/va2-rsv.r500.vadr.alt | grep MZ516105 | grep 5620 | awk '{ print $2 }' | $VADREASELDIR/esl-selectn 1 - > ex8.list
$ cat ex8.list
OR326763.1
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex8.list > ex8.fa

# re-run v-annotate.pl on:
$ v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex8.fa va-ex8
```

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
to save some time, on a random subset of 10:

```
$ grep indf3pst va2-rsv.r500/va2-rsv.r500.vadr.alt | grep MZ516105 | grep 5620 | awk '{ print $2 }' | $VADREASELDIR/esl-selectn 10 - > ex8.10.list
$ v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex8.10.fa va-ex8.10
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
$ cat va-ex8/va-ex8.vadr.alt
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
position and also the existing stop codon 21 nucleotides
downstream. I've also added `*` characters at the bottom of the
alignment at positions where the sequence and reference model are
identical. Note that the final few codons of the CDS before the
expected stop have the highest number of mismatches, which explains
the `indf3pst` alert for this CDS - the blastx alignment stopped prior
to the final few codons.

Because this is such a common characteristic of RSV B sequences, we'd
like our model to allow for it and not report any fatal alerts when it
occurs. To do this, we can modify our model by adding an **alternative
feature** for the attachment glycoprotein CDS. 

Adding an alternative feature for a CDS requires two steps:

1. Manually edit the model info file
`rsv-models2/rsv.minfo` in a text editor. 

2. Add one (or more) protein sequences to the blastx protein
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

Note the different stop position *and* the new key/value pair:
`alternative_ftr_set="attachment(cds)"`. We also need to update the
first line above to include this new key/value pair. This will inform
`v-annotate.pl` that these two CDS are members of the same alternative
feature set, and that only one of them should be annotated in the output
`.tbl` file. `v-annotate.pl` will select the CDS feature that has fewer fatal
alerts and annotate only that one. If they have the same number of fatal
alerts, the one that occurs first in the model info file will be
annotated. 

We also want to add a new `gene` feature line that will be coupled
with the new `CDS` feature line (same coordinates) so that the correct
gene coordinates will be annotated based on which CDS is annotated. 
(Importantly: we only need to do this because we want the annotated `gene`
feature to have the same coordinates as the annotated `CDS` feature. If
this wasn't the case *and* the `gene` boundaries completely spanned the
`CDS` boundaries (potentially with extra nucleotides for the `gene`), then
we would **not** need to add a new `gene` feature line.)

```
FEATURE MZ516105 type:"gene" coords:"4688..5641:+" parent_idx_str:"GBNULL" gene:"G" alternative_ftr_set="attachment(gene)" alternative_ftr_set_subn:"attachment(cds).1"
```

This new `gene` feature line includes
`alternative_ftr_set="attachment(gene)"` (it is important that the
value in quotes ("attachment(gene)") is
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
features correspond specifically to the first and second CDS features
that have the `alternative_ftr_set` value `attachment(cds)`.

<a name="step6-alternative-blastx"></a>Now we can move onto step 2,
adding a protein to the blastx protein library. As we did above when
addressing the common `deletinp` alert, we want to find a
representative sequence out of all the CDS sequences that stop at
position `5641`. I repeated the steps detailed above using the set of
41 sequences that end at this position as candidates and ended up
choosing `OK654726.1` as the representative. I then added it to the
blastx library using the `build-add-to-blast-db.pl` script. The steps
are shown below (and discussed in more detail for the `deletinp` alert
example above). One difference here is that we need to make sure that
we include the stop position in `OK654726.1` that aligns to the new
stop position at reference position `5641` which differs from what is
reported in the `.ftr` file. We may need to consult the alignment of
`OK654726.1` to determine this (after maybe rerunning
`v-annotate.pl`).

```
# determine number of sequences that have attachment glycoprotein ended at 5641 
$ grep 5641 va2-rsv.r500/*alt | grep MZ516105 | grep mutendex | wc -l
41
$ grep 5641 va2-rsv.r500/*alt | grep MZ516105 | grep mutendex | awk '{ print $2 }' > ex8.41.list

# fetch out 41 sequences and extract only the attachment glycoprotein alignment
$ $VADREASELDIR/esl-alimanip --seq-k ex8.41.list va2-rsv.r500/va2-rsv.r500.vadr.MZ516105.align.stk > ex8.41.stk 
$ $VADREASELDIR/esl-alimask -t --t-rf ex8.41.stk 4688..5641 > ex8.41.ag.stk

# convert to fasta and remove any seqs with ambiguous nts
$ $VADREASELDIR/esl-reformat fasta ex8.41.ag.stk > ex8.41.ag.fa
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl ex8.41.ag.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | awk '{ printf("%s\n", $1); }' | wc -l
40
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl ex8.41.ag.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | awk '{ printf("%s\n", $1); }' > ex8.40.list
$ $VADREASELDIR/esl-alimanip --seq-k ex8.40.list ex8.41.ag.stk > ex8.40.stk

# determine average percent id and choose representative
$ $VADREASELDIR/esl-alipid ex8.40.stk > ex8.40.alipid
$ perl $VADRSCRIPTSDIR/miniscripts/esl-alipid-per-seq-stats.pl ex8.40.alipid > ex8.40.alipid.perseq
$ grep -v ^\# ex8.40.alipid.perseq | sort -rnk 2 | head -n 1
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
$ v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex8.fa va-ex8.2
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

And we can verify that the attachment glycoprotein CDS was annotated
using the alternative stop ending at `5641` in the `.ftr` file, by looking in the `model coords` column:

```
$ head -n 3 va-ex8.2/va-ex8.2.vadr.ftr
#     seq           seq                  ftr   ftr                         ftr  ftr  par                                                                                                seq          model  ftr   
#idx  name          len  p/f   model     type  name                        len  idx  idx  str  n_from   n_to  n_instp  trc  5'N  3'N  p_from  p_to  p_instp   p_sc  nsa  nsn         coords         coords  alerts
#---  ----------  -----  ----  --------  ----  -------------------------  ----  ---  ---  ---  ------  -----  -------  ---  ---  ---  ------  ----  -------  -----  ---  ---  -------------  -------------  ------
$ grep 5641 va-ex8.2/va-ex8.2.vadr.ftr
1.13  OR326763.1  15191  PASS  MZ516105  gene  G                           954   14   -1    +    4639   5592        -  no     0    0       -     -        -      -    1    0   4639..5592:+   4688..5641:+  -     
1.14  OR326763.1  15191  PASS  MZ516105  CDS   attachment_glycoprotein     954   16   -1    +    4639   5592        -  no     0    0    4639  5589        -   1537    1    0   4639..5592:+   4688..5641:+  -     
```

The 5641 stop codon was the most common alternative to the `5620` stop
codon in `MZ516105`, but it wasn't the only alternative. Looking again
at the list of remaining examples of other `mutendex` failures (by
removing any that include `5641` with `grep -v`:

```
$ cat va2-rsv.r500/va2-rsv.r500.vadr.alt | grep mutendex | grep -v 5641 | awk '{ printf("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
     49 KY654518 attachment_glycoprotein 5647..5649:+
     15 MZ516105 attachment_glycoprotein 5627..5629:+
      3 KY654518 large_polymerase 15059..15061:+
      1 KY654518 nonstructural_protein_2 1048..1050:+
      1 KY654518 M2-1_protein 8492..8494:+
```

The most common alternative now is for the `KY654518` model ending at
position `5649`, which is 3 nucleotide downstream of the expected stop
at `5646`:

```
$ grep attachment rsv-models2/rsv.minfo | grep KY654518
FEATURE KY654518 type:"CDS" coords:"4681..5646:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein"
```

To address this we can repeat the drill we just did for the `5641` endpoint
for the `MZ516105` model, and add an alternative feature by adding to
and modifying existing lines in the model info file, and then adding a
representative protein for the coordinates `4681..5649:+` to the
`KY654518` model blastx library.

<details>

<summary>
Click to expand commands used to address `KY654518` attachment glycoprotein ending at `5649`
</summary>

```
# manually edit the rsv.minfo file to include the alternative feature
# by replacing the two lines corresponding to the `KY654518` attachment 
# glycoprotein with these four lines:
FEATURE KY654518 type:"gene" coords:"4681..5646:+" parent_idx_str:"GBNULL" gene:"G" alternative_ftr_set:"attachment(gene)" alternative_ftr_set_subn:"attachment(cds).1"
FEATURE KY654518 type:"gene" coords:"4681..5649:+" parent_idx_str:"GBNULL" gene:"G" alternative_ftr_set:"attachment(gene)" alternative_ftr_set_subn:"attachment(cds).2"
FEATURE KY654518 type:"CDS" coords:"4681..5646:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)"
FEATURE KY654518 type:"CDS" coords:"4681..5649:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)"
        
# find representative protein sequence for the blastx protein library
# fetch out all sequence names with this stop codon
$ grep 5647..5649 va2-rsv.r500/*alt | grep KY654518 | grep mutendex | awk '{ print $2 }' > ex9.49.list
# extract the attachment glycoprotein alignment
$ $VADREASELDIR/esl-alimanip --seq-k ex9.49.list va2-rsv.r500/va2-rsv.r500.vadr.KY654518.align.stk > ex9.49.stk 
$ $VADREASELDIR/esl-alimask -t --t-rf ex9.49.stk 4681..5649 > ex9.49.ag.stk

# convert to fasta and remove any seqs with ambiguous nts
$ $VADREASELDIR/esl-reformat fasta ex9.49.ag.stk > ex9.49.ag.fa
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl ex9.49.ag.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | awk '{ printf("%s\n", $1); }' > ex9.47.list
$ $VADREASELDIR/esl-alimanip --seq-k ex9.47.list ex9.49.ag.stk > ex9.47.stk

# determine average percent id and choose representative
$ $VADREASELDIR/esl-alipid ex9.47.stk > ex9.47.alipid
$ perl $VADRSCRIPTSDIR/miniscripts/esl-alipid-per-seq-stats.pl ex9.47.alipid > ex9.47.alipid.perseq
$ grep -v ^\# ex9.47.alipid.perseq | sort -rnk 2 | head -n 1
KJ643564.1  95.374  KU316171.1  89.410  KJ643503.1  100.000

# add representative to blastx db
$ perl $VADRSCRIPTSDIR/miniscripts/build-add-to-blast-db.pl  \
rsv-models2/rsv.minfo \
rsv-models2 \
KY654518 \
KJ643564 \
4672..5568:+ \
4681..5649:+ \
vb-ex9
```

</details>

After that, we want to address the `MZ516105` stop codon ending at
`5629` that occurs in 15 of our training sequences, following the
same steps.

<details>

<summary>
Click to expand commands used to address `MZ516105` attachment glycoprotein ending at `5629`
</summary>

```
# manually edit the rsv.minfo file to include the alternative feature
# by adding two lines corresponding to the `MZ516105` attachment 
# glycoprotein to the existing four, resulting in these six:
FEATURE MZ516105 type:"gene" coords:"4688..5620:+" parent_idx_str:"GBNULL" gene:"G" alternative_ftr_set:"attachment(gene)" alternative_ftr_set_subn:"attachment(cds).1"
FEATURE MZ516105 type:"gene" coords:"4688..5641:+" parent_idx_str:"GBNULL" gene:"G" alternative_ftr_set:"attachment(gene)" alternative_ftr_set_subn:"attachment(cds).2"
FEATURE MZ516105 type:"gene" coords:"4688..5629:+" parent_idx_str:"GBNULL" gene:"G" alternative_ftr_set:"attachment(gene)" alternative_ftr_set_subn:"attachment(cds).3"
FEATURE MZ516105 type:"CDS" coords:"4688..5620:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)"
FEATURE MZ516105 type:"CDS" coords:"4688..5641:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)"
FEATURE MZ516105 type:"CDS" coords:"4688..5629:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)"

# find representative protein sequence for the blastx protein library
# fetch out all sequence names with this stop codon
$ grep MZ516105 va2-rsv.r500/*alt | grep mutendex | grep attachment | awk '{ printf("%s %s\n", $2, $12); }' | grep 5627..5629 | awk '{ print $1 }' > ex10.15.list
# extract the attachment glycoprotein alignment
$ $VADREASELDIR/esl-alimanip --seq-k ex10.15.list va2-rsv.r500/va2-rsv.r500.vadr.MZ516105.align.stk > ex10.15.stk 
$ $VADREASELDIR/esl-alimask -t --t-rf ex10.15.stk 4688..5629 > ex10.15.ag.stk

# convert to fasta and remove any seqs with ambiguous nts
$ $VADREASELDIR/esl-reformat fasta ex10.15.ag.stk > ex10.15.ag.fa
$ perl $VADRSCRIPTSDIR/miniscripts/count-ambigs.pl ex10.15.ag.fa | awk '{ printf("%s %s\n", $1, $2); }' | grep " 0" | awk '{ printf("%s\n", $1); }' > ex10.14.list
$ $VADREASELDIR/esl-alimanip --seq-k ex10.14.list ex10.15.ag.stk > ex10.14.stk

# determine average percent id and choose representative
$ $VADREASELDIR/esl-alipid ex10.14.stk > ex10.14.alipid
$ perl $VADRSCRIPTSDIR/miniscripts/esl-alipid-per-seq-stats.pl ex10.14.alipid > ex10.14.alipid.perseq
$ grep -v ^\# ex10.14.alipid.perseq | sort -rnk 2 | head -n 1
JX576758.1  96.315  KU316128.1  94.480  KP258713.1  98.200

# add representative to blastx db
$ perl $VADRSCRIPTSDIR/miniscripts/build-add-to-blast-db.pl \
rsv-models2/rsv.minfo \
rsv-models2 \
MZ516105 \
JX576758 \
4690..5637:+ \
4688..5629:+ \
vb-ex10
```
</details>

#### <a name="step6-m2start"></a> Modeling `NC_038235`'s M2-2 CDS alternative start position with an alternative feature. 
Above, when investigating `mutstart` alerts returned using the RefSeq
models, we determined that `NC_038235` had a start position for M2-2
that was six nucleotides upstream from the majority of our RSV A
trainig sequences. When we switched to using a model based on
`KY654518` we began annotating the more common start position at
position `8228`. But what if we wanted to annotate the earlier start
position for those sequences that had it? We could do that by adding
an alternative feature for the M2-2 protein CDS that started six
nucleotides upstream at position `8222`. The two existing lines
beginning with `FEATURE KY654518` in the `rsv-models2/rsv.minfo` file
that include `M2-2` would be replaced with these four lines:

```
FEATURE KY654518 type:"gene" coords:"8222..8494:+" parent_idx_str:"GBNULL" gene:"M2-2" alternative_ftr-set:"M2-2(gene)" alternative_ftr_set_subn:"M2-2(cds).1"
FEATURE KY654518 type:"gene" coords:"8228..8494:+" parent_idx_str:"GBNULL" gene:"M2-2" alternative_ftr-set:"M2-2(gene)" alternative_ftr_set_subn:"M2-2(cds).2"
FEATURE KY654518 type:"CDS" coords:"8222..8494:+" parent_idx_str:"GBNULL" gene:"M2-2" product:"M2-2 protein" alternative_ftr-set:"M2-2(cds)"
FEATURE KY654518 type:"CDS" coords:"8228..8494:+" parent_idx_str:"GBNULL" gene:"M2-2" product:"M2-2 protein" alternative_ftr-set:"M2-2(cds)"
```

Importantly, the new alternative starting at `8222` needs to go first,
because `v-annotate.pl` will choose to annotate the alternative which
has the fewest fatal alerts, and in the case of ties, it will choose
the alternative that comes first in the model info file. Because
`NC_038235` includes a valid start codon beginning at reference
positions `8222` and `8228` (using `KY654518` as a reference
coordinate system) as shown in the `RF` line in the [alignment
above](#step4-mutstart-aln), it's important the `8222` feature comes
before the `8228` feature in the model info file, so that the `8222` start
is chosen for `NC_038235`. 

We'd also need to add a new protein to the blastx database that
include the two additional amino acids at the beginning of the longer
protein. We could add `NC_038235`'s M2-2 protein. After doing that, if
we ran `v-annotate.pl` on `NC_038235` it should be annotated using the
earlier start at reference position `8222`.

### <a name="step6-exception"></a> Allowing an alert exception for specific reference positions

After addressing the attachment glycoprotein `mutendex` alerts above,
the next most common `mutendex` alert is for the `large polymerase` at
positions `15059..15061:+` for 3 sequences. At this stage I would stop
and move on to other types of alerts, because with only 3 sequences it
is not as obvious that this is a legitimate type of variability that
we want our models to allow. There is probably lower hanging fruit
that we can address. It makes sense to rerun all 500 of the training
sequences using the updated models at this point, to make sure that we
base the decision on which alert instances to investigate next based
on performance using the current (updated) models.

To rerun the training set:
```
v-annotate.pl --out_stk --mdir rsv-models2 --mkey rsv rsv.r500fa va3-r500
```

Our new results:
```
# Summary of classified sequences:
#
#                                 num   num   num
#idx  model     group  subgroup  seqs  pass  fail
#---  --------  -----  --------  ----  ----  ----
1     KY654518  RSV    A          286   223    63
2     MZ516105  RSV    B          214   167    47
#---  --------  -----  --------  ----  ----  ----
-     *all*     -      -          500   390   110
-     *none*    -      -            0     0     0
#---  --------  -----  --------  ----  ----  ----
```

Previously 233 sequences failed, so we've cut the number of failing
sequences in half with our recent modifications.

What are the most common alerts remaining? Ranking the lines with
fatal alerts in the 
`va3-rsv.r500/va2-rsv.r500.vadr.alc` file by number of sequences:

```
#     alert     causes   short                               per    num   num  long       
#idx  code      failure  description                        type  cases  seqs  description
#---  --------  -------  -----------------------------  --------  -----  ----  -----------
32    indf3pst  yes      INDEFINITE_ANNOTATION_END       feature     58    56  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint
33    deletinp  yes      DELETION_OF_NT                  feature     46    46  too large of a deletion in protein-based alignment
24    cdsstopn  yes      CDS_HAS_STOP_CODON              feature     36    35  in-frame stop codon exists 5' of stop position predicted by homology to reference
30    indf5pst  yes      INDEFINITE_ANNOTATION_START     feature     25    22  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint
23    unexleng  yes      UNEXPECTED_LENGTH               feature     16    15  length of complete coding (CDS or mat_peptide) feature is not a multiple of 3
16    lowcovrg  yes      LOW_COVERAGE                   sequence     13    13  low sequence fraction with significant similarity to homology model
26    fsthicft  yes      POSSIBLE_FRAMESHIFT_HIGH_CONF   feature     12    12  high confidence possible frameshift in CDS (frame not restored before end)
20    mutendcd  yes      MUTATION_AT_END                 feature      8     8  expected stop codon could not be identified, predicted CDS stop by homology is invalid
19    mutstart  yes      MUTATION_AT_START               feature      6     6  expected start codon could not be identified
22    mutendex  yes      MUTATION_AT_END                 feature      5     5  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position
25    cdsstopp  yes      CDS_HAS_STOP_CODON              feature      4     4  stop codon in protein-based alignment
17    lowsimis  yes      LOW_SIMILARITY                 sequence      3     3  internal region without significant similarity
28    indfantn  yes      INDEFINITE_ANNOTATION           feature      5     3  nucleotide-based search identifies CDS not identified in protein-based search
31    indf3gap  yes      INDEFINITE_ANNOTATION_END       feature      6     3  alignment to homology model is a gap at 3' boundary
18    ftskipfl  yes      UNREPORTED_FEATURE_PROBLEM     sequence      2     2  only fatal alerts are for feature(s) not output to feature table
21    mutendns  yes      MUTATION_AT_END                 feature      2     2  expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted start codon
29    indf5gap  yes      INDEFINITE_ANNOTATION_START     feature      4     2  alignment to homology model is a gap at 5' boundary
27    fsthicfi  yes      POSSIBLE_FRAMESHIFT_HIGH_CONF   feature      1     1  high confidence possible frameshift in CDS (frame restored before end)
```

The most common alert remaining is `indf3pst` with instances in 56
sequences followed by `deletinp` in 46
sequences. We've seen examples above on how to address each of these
alerts: by adding sequences to the blastx protein library. That
strategy is the only method for addressing `indf3pst` alerts, but
there is an additional way to address `deletinp` alerts, by adding an
alert exception to the model info file for specific reference model
regions. Instead of including another example of adding a protein to
the blastx library here, let's try specifying an alert exception.

Looking further at the 46 remaining `deletinp` alerts, we can use
the `.ftr` file to determine which regions are commonly deleted.

```
$ cat va3-rsv.r500/va3-rsv.r500.vadr.alt | grep attachment | grep deletinp | awk '{ printf("%s %s %s\n", $3, $5, $12); }' | sort | uniq -c | sort -rnk 1
     12 MZ516105 attachment_glycoprotein 5435..5494:+
      7 KY654518 attachment_glycoprotein 5461..5532:+
      6 MZ516105 attachment_glycoprotein 5399..5458:+
      5 MZ516105 attachment_glycoprotein 5447..5506:+
      4 MZ516105 attachment_glycoprotein 5471..5530:+
      4 MZ516105 attachment_glycoprotein 5405..5464:+
      3 MZ516105 attachment_glycoprotein 5372..5407:+
      2 MZ516105 attachment_glycoprotein 5378..5422:+
      2 KY654518 attachment_glycoprotein 5515..5580:+
      1 MZ516105 attachment_glycoprotein 5444..5503:+
```

Most of the alerts (37/46) are for RSV B sequences (`MZ516105`) and they
relate to deletions at the reference positions `5372..5530`. We can specify a
`deletinp` alert exception for a specific model coordinate range for
the start position of the deletion and the maximum length we want to
allow an exception for. To do that we need to know the
range of start positions for the deletion exception, and the length we
want to allow for the deletion (any deletions above this length
starting within the range of start positions will trigger a `deletinp`
alert). For `MZ516105` an exception
ranging from positions 5372..5471 with maximum length of 60
would span all of the 37 instances
in our training dataset. 

To add the exception we need to manually modify the
`rsv-models2/rsv.minfo` file in a text editor by adding the key:value
pair `deletin_exc:5435..5471:+:60` to the `FEATURE` lines for the
attachment glycoprotein CDS for the `MZ516105` model. Remember there
are currently three `FEATURE` lines because we have alternative
features for the attachment glycoprotein and we should add the
exception to all three. The `deletin_exc` alert pertains to both
`deletinp` and `deletinn` alerts, so `deletin` non-fatal alerts will
also be excepted in this region.

Alert exceptions are only possible for a subset of alerts, for the
list of possible exceptions and more information on how to use them,
see [here](annotate.md#exceptions).

The 3 relevant lines in the model info file should now be:

```
FEATURE MZ516105 type:"CDS" coords:"4688..5620:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)" deletin_exc:"5372..5471:+:60"
FEATURE MZ516105 type:"CDS" coords:"4688..5641:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)" deletin_exc:"5372..5471:+:60"
FEATURE MZ516105 type:"CDS" coords:"4688..5629:+" parent_idx_str:"GBNULL" gene:"G" product:"attachment glycoprotein" alternative_ftr_set:"attachment(cds)" deletin_exc:"5372..5471:+:60"
```

If we now rerun these 37 sequences the `deletinp` and `deletinn`
alerts should now be absent. Let's randomly sample one to check:

```
$ grep deletinp va3-rsv.r500/va3-rsv.r500.vadr.alt | grep MZ516105 | grep attachment | awk '{ print $2 }' | $VADREASELDIR/esl-selectn 1 - > ex11.list
$ cat ex11.list
OQ101882.1
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex11.list > ex11.fa
$ v-annotate.pl --keep --mdir rsv-models2 --mkey rsv ex11.fa va-ex11
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

This sequence now passes, because the `deletinp` alert was its only
fatal one. Other sequences may still have other alerts, including
`indf3pst` alerts, which could be addressed using the other strategies
above. 

### <a name="step6-cm"></a>Rebuilding the CM with additional information

Another way to update a model is to rebuild the underlying CM using a
different input alignment. The `v-build.pl` script will by default
build a CM from a single-sequence 'alignment', but can also take a
multiple alignment as input and build a profile model from it. That
profile will have position specific parameters, meaning that each
position will have a different probability distribution for each of
the four nucleotides, and a different probability of insertion and
deletion. Those parameters are learned from the input alignment. Input
of a single sequence alignment is a special case where all positions
have the same probability distributions. Additionally, the
`v-build.pl` input alignment can have secondary structure annotation,
and the resulting CM will model the expected secondary structure and
use it when aligning input sequences. In summary, some reasons to
provide an alignment to `v-build.pl` are:

1. to increase the sequence diversity the model can handle
2. to allow the model to put insertions and deletions in 'expected'
   places (example below)
3. to add secondary structure to the model so sequences will be
   aligned based on sequence and structure (example
   [here](#https://github.com/ncbi/vadr/wiki/Rfam-based-structural-annotation-of-a-viral-genome-sequence))

As an example, we can build a new CM for one of our RSV models that
does a more consistent job of modelling the deletion in the attachment
glycoprotein CDS. (I'm including this example not because it will
address any common fatal alert instances, but just to provide an
example of rebuilding the CM.)  Below is a doctored version of the
alignment in `va3-r500/va3-r500.align.KY654518.align.stk` created
in one of the steps above, with some sequences removed and truncated
to the reference positions `5450..5600` (with the command
`$VADREASELDIR/esl-alimask -t --t-rf
va3-rsv.r500/va3-rsv.r500.align.KY654518.align.stk 5450..5600`).  This
alignment demonstrates the variability in the placement of the
deletion that occurs in some sequences that do *not* include the
duplicated region that caused the *dupregin* alerts when we were
testing the original RefSeq-based models:


```
OR143220.1         AACACACAAGTCAAGAGGTAACCCTCCACTCAACCACCTCCGAAGGCTATCCAAACCCATCACAAGTCTATACAACATCCGGTCAAGAGGAAACCCTCCACTCAACTACTTCCGAAGACTATCCAAGCCCATCACAAGTCCATACAACATC
#=GR OR143220.1 PP *******************************************************************************************************************************************************
KX655635.1         AACACACAAGTCAAGAGGAAACCCTTCACTCAACCACCTCCGAAGGC------------------------------------------------------------------------AATCCAAGCCCATCACAATTCTATACAACATC
#=GR KX655635.1 PP ************99999999988888887777777766666666555........................................................................555566666777777778888889999999**
KJ627366.1         AAAACACAAGTCAAAAGGAAACCCTCCACTCAACTACTCCCGAAG------------------------------------------------------------------------GCAATCCAAGCCCTTCACAAGTCTATGCAACATC
#=GR KJ627366.1 PP ************999999998888888887777776666655555........................................................................5555666666777777788888888899999999
MK109787.1         AACACACAAGTCAAGAGGAAACCCTCCACTCAACCACCTCCGAAGGCA------------------------------------------------------------------------ACCCAAGCCCATCACAAGTCTATACAACATC
#=GR MK109787.1 PP *************99999999988888888877777777666666655........................................................................555566777777788888888899999999*
KJ627665.1         AACACACAAGTCAAGAGGAAACCCTCCATTCAACCTCCTCCGAAGGCAA------------------------------------------------------------------------TACAAGCCCTTCACAAATCTATACAACATC
#=GR KJ627665.1 PP *************999999999988888888777777766666666555........................................................................555666777777888888899999999***
MK810782.1         AACTCACAAGTCAAATGGAAACCTTCCACTCAACCTCCTCCGAAGGCAAT------------------------------------------------------------------------CTAAGCCCTTCTCAAGTCTCCACAACATC
#=GR MK810782.1 PP ************99999988888888877777777666666666665555........................................................................5555566666667777777788999999*
KJ627688.1         AACACACAAGTCAAGAGGAAACCCTCCATTCAACCTCCTCCGAAGGCAA------------------------------------------------------------------------TACAAGCCCTTCACAAATCTATACAACATC
#=GR KJ627688.1 PP *************999999999988888888777777766666666555........................................................................555666777777888888899999999***
ON237257.1         AACACACAAGTCAAGAGGAAACCCTCCACTCAACCACCTCCGAAGGC------------------------------------------------------------------------TATCTAAGCCCATCACAAGTCTATACAACATC
#=GR ON237257.1 PP ************99999999888888887777777666666665555........................................................................55555666666677777777888888899999
MG813989.1         AACTCACAAGTCAAATGGAAACCTTCCACTCAACTTCCTCCGAAGGTAATC------------------------------------------------------------------------CAAGCCCTTCTCAAGTCTCCATAACATC
#=GR MG813989.1 PP ************999999998888888888877777766666666655555........................................................................55556666666777777778889999**
#=GC SS_cons       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF            AACACACAAGTCAAGAGGAAACCCTCCACTCAACCACCTCCGAAGGCTATCTAAGCCCATCACAAGTCTATACAACATCCGGTCAAGAGGAAACCCTCCACTCAACCACCTCCGAAGGCTATCTAAGCCCATCACAAGTCTATACAACATC
#=GC RFCOLX....    0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#=GC RFCOL.X...    5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
#=GC RFCOL..X..    4444444444444444444444444444444444444444444444444455555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555556
#=GC RFCOL...X.    5555555555666666666677777777778888888888999999999900000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
#=GC RFCOL....X    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
```

The `cmalign` program's alignment algorithm optimizes the expected
accuracy of the alignment. The expected accuracy of each aligned
nucleotide is shown in the `PP` lines of the alignment with `*`
indicating the highest level of expected accuracy, with `9` the second
highest level, and `8` the third and so on (as explained more
[here](alerts.md#pp)). Note that near the deletion is the lowest
expected accuracy. This makes sense, because those nucleotides could
reasonably be aligned at the opposite end of the deletion because this
deletion is actually just a lack of a short duplicated region.  The
issue responsible for this alignment inconsistency is that the CM does
not have any information about where the deletion should occur,
because it is only based on the one `KY654518` sequence which does not
have the deletion at all. We can add information about where the
deletion should occur (and only that information) to the CM by
rebuilding it from an alignment of two sequences: the original
`KY654518` and a synthetic sequence that is a copy of `KY654518` but
with the common deletion in the specific positions we want it to be
placed in output `v-annotate.pl` alignments.

Because our goal is to make the alignment of this region more
consistent, it makes sense to find the *average* position span for
this deletion. We may find that the region from reference positions
5496 to 5567 (72 positions) is the average. We can manually
create a two sequence Stockholm alignment file with `KY654518`
duplicated by starting with the file
`va-rsv2/va-rsv2.vadr.KY654518.align.stk` that was created with the
`v-annotate.pl` command:

```
$ v-annotate.pl --out_stk --mdir rsv-models2 --mkey rsv rsv-models2/rsv.fa va-rsv2
```

We can reformat this to a special type of Stockholm format referred to
as Pfam in the Easel and Infernal codebases/documentation that has
only one line per sequence (as opposed to the standard interleaved
Stockholm format). This will make it easier to duplicate the sequence.

```
$VADREASELDIR/esl-reformat pfam va-rsv2/va-rsv2.vadr.KY654518.align.stk > KY654518.2.pfam
```

Next we need to open this file in a text editor, duplicate the
sequence line that begins with `KY654518` and rename the second
sequence something like `KY654518-5496del72`. Then remove the
sequence in reference positions `5496` to `5567` in this second sequence,
replacing those nucleotides with `-` characters. After that you can 
and remove the line that starts with `#=GR KY654518 PP` as that
is irrelevant for the `cmbuild` step. When you are finished, save the
file, and then reformat it back to interleaved Stockholm format with:

```
$VADREASELDIR/esl-reformat stockholm KY654518.2.pfam > KY654518.2.stk
```

A copy of the `KY654518.2.stk` alignment file can be found in
[`vadr/documentation/build-files/KY654518.2.stk`](build-files/KY654518.2.stk)

The next step is to use this alignment to build a new CM file. We can
do this using the `cmbuild` program which was called by `v-build.pl`
when we built the initial model. We'll
want to use similar command-line options to what `v-build.pl` used,
which we can find in the `.cmd` output file from `v-build.pl`:

```
$ grep cmbuild KY654518/KY654518.vadr.cmd
/usr/local/vadr-install/infernal/binaries/cmbuild -n KY654518 --verbose --noss --noh3pri --Egcmult 1.63645 KY654518/KY654518.vadr.cm KY654518/KY654518.vadr.stk > KY654518/KY654518.vadr.cmbuild
```

So we'll use the `-n KY654518 --verbose --noss --noh3pri --Egcmult
1.63645`  options and we'll add one more important one: `--hand` which
informs `cmbuild` to maintain the existing reference positions in the
alignment, which correspond to the `KY654518` sequence, instead of
inferring new ones. This is extremely important whenever rebuilding a
CM for a VADR model because the coordinates of all of the features in
the existing model info file are with respect to the `KY654518`
sequence. If the reference positions change, then the model info file
coordinates would need to be updated accordingly.

So the `cmbuild` command is:
```
$VADRINFERNALDIR/cmbuild -n KY654518 --verbose --noss --noh3pri --Egcmult 1.63645 --hand KY654518.2.cm KY654518.2.stk
```

This may take up to an hour to complete. 

The final step is to remake the `rsv.cm` model file by combining our
new model with the existing `MZ516105` model. We can use the `cmfetch`
program to help with this:

```
# create a temporary CM file `new.rsv.cm`
$ $VADRINFERNALDIR/cmfetch rsv-models2/rsv.cm MZ516105 > new.rsv.cm
$ cat KY654518.2.cm >> new.rsv.cm

# copy it over our previous model (after saving a copy just in case):
$ cp rsv-models2/rsv.cm ./old.rsv.cm
$ cp new.rsv.cm rsv-models2/rsv.cm

# and remember to re-press this new file:
$ rm rsv-models2/rsv.cm.*
$ $VADRINFERNALDIR/cmpress rsv-models2/rsv.cm
```

We can test out our new model on the set of sequences that had the
jagged alignment above:

```
$ cat ex13.list
OR143220.1
KX655635.1
KJ627366.1
MK109787.1
KJ627665.1
MK810782.1
KJ627688.1
ON237257.1
MG813989.1
$ $VADREASELDIR/esl-sfetch -f rsv.r500.fa ex13.list > ex13.fa
$ v-annotate.pl --out_stk --mdir rsv-models2 --mkey rsv ex13.fa va-ex13
```

Then if we look at the relevant region of the alignment: 

```
$ $VADREASELDIR/esl-alimask -t --t-rf va-ex13/va-ex13.vadr.KY654518.align.stk 5450..5600
```

```
OR143220.1         AACACACAAGTCAAGAGGTAACCCTCCACTCAACCACCTCCGAAGGCTATCCAAACCCATCACAAGTCTATACAACATCCGGTCAAGAGGAAACCCTCCACTCAACTACTTCCGAAGACTATCCAAGCCCATCACAAGTCCATACAACATC
#=GR OR143220.1 PP *******************************************************************************************************************************************************
KX655635.1         AACACACAAGTCAAGAGGAAACCCTTCACTCAACCACCTCCGAAGG------------------------------------------------------------------------CAATCCAAGCCCATCACAATTCTATACAACATC
#=GR KX655635.1 PP **********************************************........................................................................899******************************
KJ627366.1         AAAACACAAGTCAAAAGGAAACCCTCCACTCAACTACTCCCGAAGG------------------------------------------------------------------------CAATCCAAGCCCTTCACAAGTCTATGCAACATC
#=GR KJ627366.1 PP *********************************************9........................................................................899******************************
MK109787.1         AACACACAAGTCAAGAGGAAACCCTCCACTCAACCACCTCCGAAGG------------------------------------------------------------------------CAACCCAAGCCCATCACAAGTCTATACAACATC
#=GR MK109787.1 PP **********************************************........................................................................889999***************************
KJ627665.1         AACACACAAGTCAAGAGGAAACCCTCCATTCAACCTCCTCCGAAGG------------------------------------------------------------------------CAATACAAGCCCTTCACAAATCTATACAACATC
#=GR KJ627665.1 PP **********************************************........................................................................889999***************************
MK810782.1         AACTCACAAGTCAAATGGAAACCTTCCACTCAACCTCCTCCGAAGG------------------------------------------------------------------------CAATCTAAGCCCTTCTCAAGTCTCCACAACATC
#=GR MK810782.1 PP **********************************************........................................................................899******************************
KJ627688.1         AACACACAAGTCAAGAGGAAACCCTCCATTCAACCTCCTCCGAAGG------------------------------------------------------------------------CAATACAAGCCCTTCACAAATCTATACAACATC
#=GR KJ627688.1 PP **********************************************........................................................................889999***************************
ON237257.1         AACACACAAGTCAAGAGGAAACCCTCCACTCAACCACCTCCGAAGG------------------------------------------------------------------------CTATCTAAGCCCATCACAAGTCTATACAACATC
#=GR ON237257.1 PP **********************************************........................................................................9********************************
MG813989.1         AACTCACAAGTCAAATGGAAACCTTCCACTCAACTTCCTCCGAAGG------------------------------------------------------------------------TAATCCAAGCCCTTCTCAAGTCTCCATAACATC
#=GR MG813989.1 PP **********************************************........................................................................77999****************************
#=GC SS_cons       :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF            AACACACAAGTCAAGAGGAAACCCTCCACTCAACCACCTCCGAAGGCTATCTAAGCCCATCACAAGTCTATACAACATCCGGTCAAGAGGAAACCCTCCACTCAACCACCTCCGAAGGCTATCTAAGCCCATCACAAGTCTATACAACATC
#=GC RFCOLX....    0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#=GC RFCOL.X...    5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
#=GC RFCOL..X..    4444444444444444444444444444444444444444444444444455555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555556
#=GC RFCOL...X.    5555555555666666666677777777778888888888999999999900000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
#=GC RFCOL....X    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
```

The updated model now aligns all the deletions in the same place, and
with higher expected accuracy values. 

This is just an example of how updating the training alignment and
rebuilding the model can effect the output alignment. In this case, it
actually does not change any annotations or pass/fail outcomes for
these 9 sequences, but there are other situations for other viruses
where updating the model could have a more significant impact.

### <a name="step6-miscfeat"></a>Treating a feature as non-essential by allowing it to be a `misc_feature`

For some viruses, some features may be non-essential and so can
tolerate mutations that disrupt or modify the function, such as early
stop codons, or frameshifts. For these features, we may want to allow
alerts to be reported but not be fatal. VADR allows you to specify
features as non-essential by adding a `misc_not_failure:"1"` key/value
pair to the relevant `FEATURE` lines in the model info file as
explained more [here](annotate.md#mnf). Such features will be
annotated as `misc_feature` if they include a normally fatal alert,
but the sequence will not fail because of such alerts.

When building RSV models, we could have defined the attachment
glycoprotein CDS, which is responsible for many of the fatal alerts,
as non-essential. However, this would have meant that it would often
be annotated as a `misc_feature`. By modifying the model through
adding proteins to the blastx library as well as the other strategies
above, we have specified the range of possible variability we want to
allow in the attachment glycoprotein CDS without reporting an alert for it,
while still validating and annotating it at as a CDS. For other
viruses, treating some features as non-essential can be a useful
strategy. Since August 2021, at least up until the time of writing (October
2023), VADR-based annotation of SARS-CoV-2 sequences submitted
to GenBank treats ORF3a, ORF6, ORF7a, ORF7b, ORF8, and ORF10 CDS as
well as the Coronavirus 3' stem-loop II-like motif (s2m) as
non-essential using this `misc_not_failure` strategy.

### <a name="step6-altpass"></a>Making an alert non-fatal using the `--alt_pass` option

For some viruses, specific fatal alerts are so common that we may want
to make them non-fatal. For example, the Mpox genome has several
repetitive regions that cause `dupregin` and `discontn` alerts for
nearly all mpox sequences. We could try to define alert exceptions to
allow the `dupregin` alerts, but there are no exceptions supported for
`discontn`. An alternative strategy is to specify that these two alerts be
considered non-fatal using the `v-annotate.pl` option: `--alt_pass
dupregin,discontn`. An example of using this option can be found
[here](annotate.md#examplealtpass).

### <a name="step6-summary"></a> Final summary of RSV model modifications

I followed the six step procedure above when creating the [VADR RSV
models](https://bitbucket.org/nawrockie/vadr-models-rsv/src/master/),
which can be used with `v-annotate.pl` to validate and annotate RSV
sequences as described
[here](https://github.com/ncbi/vadr/wiki/RSV-annotation). Those RSV
models are based on the `KY654518` and `MZ516105` reference sequences,
but the exact steps I took when creating them are not identical to
those listed above (I created this tutorial several months after I had
finished building the RSV models).  For example, I added proteins to
the blastx databases and alternative features to the model info files
that are not listed above. The table below summarizes most of the
additions I made to the original `KY654518` and `MZ516105` models that
are built by `v-build.pl`:

| model    | type of modification | feature | detail | 
----------|---------------------- |---------|--------|
`KY654518`  | added proteins to blastx library (13) | attachment glycoprotein (CDS) | proteins added (name format: `source accession:source coordinates/model coordinates`): `OM857255.1:4629..5591:+/4681..5643:+`, `KU316164.1:4611..5504:+/4681..5646:+`,  `AF065254.1:16..909:+/4681..5646:+`, `OK649616.1:4670..5563:+/4681..5646:+`, `hybrid:KY654518.1:4681..4695:+:AF065410.1:1..879:+/4681..5646:+`, `KF826850.1:4675..5568:+/4681..5646:+`, `KU316092.1:4620..5516:+/4681..5646:+`, `NC_038235.1:4688..5584:+/4681..5649:+`, `MZ515659.1:4681..5649:+/4681..5649:+`, `HQ699266.1:1..897:+/4681..5649:+`, `KJ641590.1:4630..5526:+/4681..5649:+`, `OK649616.1:4670..5566:+/4681..5649:+`, `M17212.1:16..912:+/4681..5649:+` |
`KY654518` | added protein to blastx library (1) | M2-1(CDS) | protein added: `OM857351.1:7614..8180:+/7669..8235:+`
`KY654518` | added alternative features (2)  | attachment glycoprotein (CDS + gene) | alternative feature coordinates: `4681..5643:+`, `4681..5649:+` | 
`KY654518` | added alert exception (1)     | attachment glycoprotein (CDS) | key/value pair added to model info file: `deletin_exc:5457..5508:+:72` | 
`KY654518` | rebuilt CM           | full model                    | added duplicate `KY654518` sequence with 72nt deletion after position `5496` |
| | | | 
`MZ516105` | added proteins to blastx library (21) | attachment glycoprotein (CDS) | proteins added: `MG642047.1:4666..5565:+/4688..5620:+`, `MG431253.1:4674..5567:+/4688..5620:+`,  `LC474547.1:4663..5595:+/4688..5620:+`, `MZ962122.1:1..933:+/4688..5620:+`, `KC297442.1:1..933:+/4688..5620:+`, `LC311384.1:1..933:+/4688..5620:+`, `MZ515748.1:4689..5621:+/4688..5620:+`, `MT040088.1:4679..5572:+/4688..5620:+`, `KJ627364.1:4618..5550:+/4688..5620:+`, `MH760718.1:4597..5529:+/4688..5620:+`, `KP856962.1:4618..5505:+/4688..5629:+`, `KU950619.1:4663..5604:+/4688..5629:+`, `KP258745.1:4620..5507:+/4688..5629:+`, `KU316181.1:4618..5505:+/4688..5629:+`, `MF185751.1:4640..5527:+/4688..5629:+`, `KC297470.1:1..882:+/4688..5629:+`, `KJ627249.1:4618..5565:+/4688..5635:+`, `KU316144.1:4618..5517:+/4688..5641:+`, `MN365572.1:4676..5629:+/4688..5641:+`, `OK649740.1:4675..5574:+/4688..5641:+`, `NC_001781.1:4690..5589:+/4688..5641:+` |
`MZ516105` | added protein to blastx library (1) | RNA-dependent RNA polymerase (CDS)| proteins added: `LC474543.1:8538..15017:+/8560..15039:+` |
`MZ516105` | added alternative features (2)  | attachment glycoprotein (CDS + gene) | alternative feature coordinates: `4688..5629:+`, `4688..5635:+` | 
`MZ516105` | added alternative feature (1)  | RNA-dependent RNA polymerase (CDS + gene)| alternative feature coordinates: `8560..15039:+` |
`MZ516105` | added alert exceptions (2)      | attachment glycoprotein (CDS) | key/value pair added to model info file: `deletin_exc:5441..5441:+:60`, `insertn_exc:5392..5467:+:60` | 
`MZ516105` | rebuilt CM           | full model                         | added duplicate `MZ516105` sequence with 60nt deletion after position `5441` |

---

## <a name="limit"></a> Limitations of and alternatives to this approach

The above procedure is one possible strategy for building VADR models
for RSV. While the strategy is somewhat general, different strategies
may work better for other viruses. Below I discuss some of the
limitations of this strategy and ideas for other strategies.

### <a name="limit-ref"></a> Reference sequence selection

Above we started with the RefSeq sequences as the basis for the
original models, then determined that they were not very
representative of RSV sequences in the database. Alternatively, if we
had expert knowledge of good representative sequences, we could have
started with those. Or we could have tried to find 'centroid'
sequences that were maximally similar to all RSV sequences to begin
with. Or, if we had a favorite sequence that was extremely well
studied and annotated, we might want to start with that. The approach
above is one that is reasonable if very little is known beforehand
about the virus being modelled and its sequence diversity, but it
makes sense to take advantage of any expert knowledge you have when
picking the initial representative sequences.

### <a name="limit-training"></a> Training sequence selection

The steps above explain how to select a random subset of 500 from all
existing INSDC full length RSV sequences to use as a training
set. Alternatively, we could have removed redundancy from the set of
candidate sequences first, so that our set of 500 was not biased
towards those sequences that are overrepresented in the database. For
example, if there was a major RSV sequencing project in 2019 then
there may tend to be more sequences in the database from the
particular virus population that circulated in 2019 than of sequences
from other years. We could filter our candidate sequences by sequence
identity, or by year, or by something else, to try and deal with this
redundancy, and then choose a training set by taking a randomly subset
of that filtered set of candidate sequences.

We could also not restrict our model training to only full length
sequences, and instead use all sequences. Or we could have two
training sets: one full length and the other partial length
sequences. We could follow the procedure above based on full length
sequences, and then check how the models worked on partial length
sequences too. This is also [briefly discussed above](#step2-length).

### <a name="limit-multiple"></a> Replacing one model with multiple models

In the steps outlined above, there are examples of modifying
single-sequence based models to be more general. An alternative
strategy would be to add one or more new models built from new
sequences that would allow sequences with divergent features to
pass. This can work especially well if you are using one model for a
set of sequences that can be easily separated into two distinct
clusters based on sequence identity (or on an inferred evolutionary
tree). In that situation, one model per cluster may be the best
approach. Be careful though, when you expand the number of models you
may increase the number of common alerts due to acceptable sequence
diversity that are returned by `v-annotate.pl`, each of which you will
have to deal with through some kind of model modification. In other
words, adding models can lead to more work manually tweaking those
models.

### <a name="limit-align"></a> Alignment-based models

If you're familiar with [Infernal](http://eddylab.org/infernal/), the
software package that VADR uses to build a statistical model of a
virus, you may be confused as to why this advanced tutorial on model
building focuses on building single-sequence models, because Infernal
is typically used to build profile models from multiple sequence
alignments. Profiles have advantages over single-sequence-based
methods because they include position-specific information about the
expected nucleotide distribution and probability of insertions and
deletions at each position, whereas single sequence based models treat
all positions identically. I provided one [example above](#step6-cm)
of rebuilding a CM from multiple sequences, but even that example is
only to deal with a single deletion, and doesn't introduce any other
sequence variability into the alignment.

If you do build a multiple sequence alignment, you'll still want to
select one of the sequences to be the "reference sequence" to use for
the reference coordinate system. The coordinates/positions in the
model info file will correspond to positions in that reference
sequence. This means a reasonable approach is to first use
`v-build.pl` using the accession you've selected as the reference
sequence. Then construct your alignment, possibly by running
`v-annotate.pl` using your initial model and the `--out_stk` option. 
Then rebuild the CM using `cmbuild` with the `--hand` option
and your multiple alignment as input as explained in the [example
above](#step6-cm) and finally overwrite the `v-build.pl` created
single sequence CM with your newly created one. 

The performance difference between using multiple sequence
alignment-based models versus versus single sequence-based models will
largely depend on on how similar all of the sequences being annotated
are to the single-sequence model. If they are highly similar, it won't
make a huge difference; if there is a lot of sequence variability it
will make more of a difference.  In my experience, building CMs from
multiple alignments for well conserved viruses like RSV, does not
significantly improve the performance of `v-annotate.pl`. This is
likely because all the sequences are so similar that the position
specific parameters are not necessary to get the correct
alignment. Because single sequence models perform acceptably well for
norovirus, dengue virus and SARS-CoV-2, the VADR models used to screen
incoming GenBank submissions of those virus sequences are all based on
single sequences.

Multiple sequence alignment-based VADR models are used for one type of
sequences at GenBank - the Cytochrome C Oxidase 1 (COX1) mitochondrial
protein coding gene, which exhibits significantly more sequence
variability than any of the viruses VADR is used for. [The VADR
library of COX1
models](https://bitbucket.org/nawrockie/vadr-models-cox1) includes 78
profile models built from multiple alignments, 20 of which include
more than 100 aligned sequences.

### <a name="limit-secondary"></a> Incorporating secondary structure

If you're familiar with [Infernal](http://eddylab.org/infernal/), you
may be confused again as to why there is no step at the beginning of
the procedure above to add predicted (or known) secondary structure to
the CM. The main reason CM methods exist is to take advantage of
conserved RNA secondary structure for RNA sequence analysis. The
conserved structure information could be taken advantage
of during sequence validation and annotation. To do this, we could add a
step to create an alignment file annotated with secondary structure of
any hits found in our reference sequence using Infernal and the Rfam
database of RNA families, which includes hundreds of viral RNA
families. There is a separate documentation page on how to do that (on
the VADR GitHub wiki)
[here](https://github.com/ncbi/vadr/wiki/Rfam-based-structural-annotation-of-a-viral-genome-sequence).

For RSV, it turns out there are actually zero Rfam hits (as of release
14.9) in the RefSeq sequences. For other viruses though, if you have
the time and interest, please do try to add secondary structure. There
are at least two reasons you may want to try this: adding secondary
structure could improve the alignment accuracy, and it will allow you
to enrich your annotations by adding some structural RNA features that
are commonly absent from GenBank annotation.  The [dengue virus] and
[SARS-CoV-2 VADR
models](https://bitbucket.org/nawrockie/vadr-models-sarscov2) both
include some secondary structure in their CMs and structural RNA
features (e.g. `stem_loop`, `ncRNA`).

Once you've created a structure annotated stockholm alignment file,
you can either use it as input to `v-build.pl` using the `--stk`
option like in [this example for dengue
virus](build.md#1.0library-dengue), or you can rebuild the CM later
using `cmbuild` with the structure annotated alignment as input,
similar to the [example above](#step6-cm).

---
 
#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.

