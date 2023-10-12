# <a name="top"></a> Advanced tutorial: building an RSV model library

---

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

---

<details>
<summary>### [Iteration 1, step 1: build model(s) from initial reference sequence(s)]</summary>

<br>
### Determine good reference sequence(s) to use

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

### [Iteration 1, step 2: construct a training set and run `v-annotate.pl` on it](#advbuild-files/advbuild-i1.s2.md)

### [Iteration 1, step 3: analyze the results, and identify major and minor sequence characteristics responsible for common failure modes](#advbuild-files/advbuild-i1.s2.md)


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

