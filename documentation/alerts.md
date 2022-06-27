# <a name="top"></a> Explanations and examples of `v-annotate.pl` detailed alert and error messages

* [Output files with detailed alert and error messages](#files)
* [Explanation of sequence and model coordinate fields in `.alt` files](#coords)
* [`toy50` toy model used in examples of alert messages](#toy)
* [Examples of different alert types and corresponding `.alt` output](#examples)
  * [Frameshift](#example-frameshift)
  * [Large insertion](#example-insert)
  * [Large deletion](#example-delete)
  * [Mutated start codon](#example-start)
  * [Stop codon alerts](#example-stop)
  * [Indefinite annotation at the start of a sequence or a feature](#example-indf5)
  * [Indefinite annotation at the end of a sequence or a feature](#example-indf3)
  * [Ambiguous nucleotides at the start of a sequence or a feature](#example-ambg5)
  * [Ambiguous nucleotides at the end of a sequence or a feature](#example-ambg3)
  * [Ambiguous nucleotides in start/stop codon that starts/ends with canonical nt](#example-ambgcd)
  * [Mature peptide-specific alerts](#example-pep)
  * [Low similarity at the start of a sequence of a feature](#example-lowsim5)
  * [Low similarity at the end of a sequence of a feature](#example-lowsim3)
  * [Low similarity internal to a sequence of a feature](#example-lowsimi)
  * [Deleted feature](#example-delftr)
  * [Duplicate regions](#example-dupregin)
  * [Discontinuous similarity](#example-discontn)
  * [Indefinite strand](#example-indfstrn)
  * [Low coverage](#example-lowcovrg)
* [Posterior probability annotation in VADR output Stockholm alignments](#pp)

---

# <a name="files"></a> Output files with detailed alert and error messages

`v-annotate.pl` outputs two types of files with detailed alert/error
messages: 

  * `.alt` files
  * `.alt.list` files
  * GenBank submission portal detailed error report `.tsv` files

Documentation on the format of `.alt` files can be found
[here](formats.md#alt) and for `.alt.list` files can be found
[here](formats.md#altlist). The GenBank submission portal detailed
error report `.tsv` files are the same format as the `.alt.list`
files. 

This page includes examples of many of the different alerts and 
corresponding `.alt` file output [below](#examples).

# <a name="coords"></a> Explanation of sequence and model coordinate fields in `.alt` files

| alert code(s) | alert desc(s) | sequence coords description | model coords explanation | link to example | 
|---------------|---------------|-----------------------------|--------------------------|-----------------|
| *fsthicft*, *fsthicfi*, *fstlocft*, *fstlocfi*, *fstukcft*, *fstukcfi* | *POSSIBLE_FRAMESHIFT_HIGH_CONF*,  *POSSIBLE_FRAMESHIFT_LOW_CONF*, *POSSIBLE_FRAMESHIFT* | sequence positions of the frameshifted region | model (reference) positions of the frameshifted region, some nucleotides may be inserted **before or after** these positions | [frameshift example](#example-frameshift) | 
| *insertnn*, *insertnp* | *INSERTION_OF_NT* | sequence positions of inserted nucleotides with respect to the model | model (reference) position after which insertion occurs (always length 1) | [large insertion example](#example-insert) | 
| *deletinn*, *deletinp* | *DELETION_OF_NT*  | sequence position just prior to (5' of) deletion with respect to the model (always length 1) | model (reference) positions that are deleted in sequence | [large deletion example](#example-delete) | 
| *mutstart* | *MUTATION_AT_START*  | sequence positions of predicted start codon (length <= 3) | model (reference) positions that align to the predicted start codon | [mutated start codon example](#example-start) | 
| *mutendcd* | *MUTATION_AT_END*  | sequence positions of predicted stop codon (length <= 3) | model (reference) positions that align to the predicted stop codon | [stop codon alert examples](#example-stop) | 
| *mutendex* | *MUTATION_AT_END*  | sequence positions of 5'-most in-frame stop codon in the CDS, this stop codon will be 3' of expected stop codon position (always length 3) | model (reference) positions that align to stop codon in `sequence coords` | [stop codon alert examples](#example-stop) | 
| *mutendns* | *MUTATION_AT_END*  | will be blank (`-`) | will be blank (`-`) | [stop codon alert examples](#example-stop) | 
| *unexleng* | *UNEXPECTED_LENGTH* | sequence positions of the predicted CDS, the length of which is not a multiple of 3 | model (reference) positions that the predicted CDS align to, some nucleotides may be inserted *before or after* these positions | [stop codon alert examples](#example-stop) |
| *cdsstopn* | *CDS_HAS_STOP_CODON* | sequence positions of the 5'-most in-frame stop codon in the CDS, this stop will be 5' of expected stop codong position (always length 3) | model (reference) positions that align to stop codon in `sequence coords` | [stop codon alert examples](#example-stop) | 
| *indf5gap* | *INDEFINITE_ANNOTATION_START* | sequence position of first nucleotide aligned 3' of gap that aligns to the feature boundary (always length 1) | model (reference) position of the 5' feature boundary (always length 1) | [stop codon alert examples](#example-indf5) | 
| *indf5lcc* | *INDEFINITE_ANNOTATION_START* | sequence position of nucleotide aligned at the 5' feature boundary (always length 1) | model (reference) position of the 5' feature boundary (always length 1) | [examples of indefinite annotation at start](#example-indf5) | 
| *indf5pst* | *INDEFINITE_ANNOTATION_START* | sequence positions of the nucleotide alignment of the 5' end of the CDS **not** covered by the protein-based alignment | model (reference) positions the sequences positions in `sequence coords` are aligned to in the nucleotide alignment | [examples of indefinite annotation at start](#example-indf5) | 
| *indf5plg* | *INDEFINITE_ANNOTATION_START* | sequence positions of the protein-based alignment **not** covered by the nucleotide alignment at the 5' end of the CDS | model (reference) position of the 5' boundary of the CDS (always length 1) | [examples of indefinite annotation at start](#example-indf5) | 
| *indf3gap* | *INDEFINITE_ANNOTATION_END*   | sequence position of final nucleotide aligned 5' of gap that aligns to the feature boundary (always length 1) | model (reference) position of the 3' feature boundary (always length 1) | [examples of indefinite annotation at end](#example-indf3) | 
| *indf3lcc* | *INDEFINITE_ANNOTATION_START* | sequence position of nucleotide aligned at the 3' feature boundary (always length 1) | model (reference) position of the 3' feature boundary (always length 1) | [examples of indefinite annotation at end](#example-indf3) | 
| *indf3pst* | *INDEFINITE_ANNOTATION_START* | sequence positions of the nucleotide alignment of the 3' end of the CDS **not** covered by the protein-based alignment | model (reference) positions the sequences positions in `sequence coords` are aligned to in the nucleotide alignment | [examples of indefinite annotation at end](#example-indf3) |   
| *indf3plg* | *INDEFINITE_ANNOTATION_START* | sequence positions of the protein-based alignment **not** covered by the nucleotide alignment at the 3' end of the CDS | model (reference) position of the 3' boundary of the CDS (always length 1) |  [examples of indefinite annotation at end](#example-indf3) | 
| *ambgnt5f*, ambgnt5c* | *AMBIGUITY_AT_FEATURE_START*, *AMBIGUITY_AT_CDS_START* | sequence position(s) of stretch of 1 or more consecutive ambiguous (non-ACGTU) nts ending at the 3' end of a feature | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | [examples of ambiguous nucleotides at start](#example-ambg5) | 
| *ambgnt5s* | *AMBIGUITY_AT_END* | sequence position(s) of 1 or more consecutive ambiguous (non-ACGTU) nts starting at position 1 of the sequence | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | [examples of ambiguous nucleotides at start](#example-ambg5) | 
| *ambgnt3f*, ambgnt3c* | *AMBIGUITY_AT_FEATURE_END*, *AMBIGUITY_AT_CDS_END* | sequence position(s) of 1 or more consecutive ambiguous (non-ACGTU) nts starting at the predicted 5' end of a feature | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | [examples of ambiguous nucleotides at end](#example-ambg3) | 
| *ambgnt3s* | *AMBIGUITY_AT_END* | sequence position(s) of 1 or more consecutive ambiguous (non-ACGTU) nts ending at the final position of the sequence | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | [examples of ambiguous nucleotides at end](#example-ambg3) | 
| *ambgcd5c* | *AMBIGUITY_IN_START_CODON* | sequence position(s) of start codon | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | [examples of ambiguous nucleotides in start/stop codons](#example-ambgcd) | 
| *ambgcd3c* | *AMBIGUITY_IN_STOP_CODON* | sequence position(s) of stop codon | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | [examples of ambiguous nucleotides in start/stop codons](#example-ambgcd) | 
| *pepadjcy* | *PEPTIDE_ADJACENCY_PROBLEM* | sequence position(s) of nucleotides inserted between two mature peptide predictions that are expected to be adjacent | model (reference) position(s) corresponding to the end of the 5' mature peptide and the start of the 3' mature peptide (always length 2) | [mature peptide-specific alert examples](#example-pep) | 
| *peptrans* | *PEPTIDE_TRANSLATION_PROBLEM* | will be blank  (`-`) | will be blank (`-`) | [mature peptide-specific alert examples](#example-pep) | 
| *lowsim5c*, *lowsim5n*, *lowsim5l* | *LOW_FEATURE_SIMILARITY_START* | sequence position(s) at 5' end of sequence that have low similarity to the reference model (low similarity region overlaps with a predicted feature) | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | [examples of low similarity at start](#example-lowsim5) | 
| *lowsim5s* | *LOW_SIMILARITY_START* | sequence position(s) at 5' end of sequence (not overlapping with a feature) that have low similarity to the reference model | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to, or '`-`' if sequence is not aligned | [examples of low similarity at start](#example-lowsim5) | 
| *lowsim3c*, *lowsim3n*, *lowsim3l* | *LOW_FEATURE_SIMILARITY_END* | sequence position(s) at 3' end of sequence that have low similarity to the reference model (low similarity region overlaps with a predicted feature) | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | [examples of low similarity at end](#example-lowsim3) | 
| *lowsim3s* | *LOW_SIMILARITY_END* | sequence position(s) at 3' end of sequence (not overlapping with a feature) that have low similarity to the reference model | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to, or '`-`' if sequence is not aligned | [examples of low similarity at end](#example-lowsim3) | 
| *lowsimic*, *lowsimin*, *lowsimil* | *LOW_FEATURE_SIMILARITY* | sequence position(s) internal to a predicted feature (not including first or final position of the feature) that have low similarity to the reference model | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | [examples of internal low similarity](#example-lowsimi) | 
| *lowsimis* | *LOW_SIMILARITY* | sequence position(s) internal to a sequence (not including first or final position of the sequence) and not overlapping with a feature that have low similarity to the reference model | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to, or '`-`' if sequence is not aligned | [examples of internal low similarity](#example-lowsimi) | 
| *deletins* | *DELETION_OF_FEATURE* | will be blank (`-`) | model (reference) positions that correspond to the feature that is deleted in the sequence | [deleted feature examples](#example-deletin) | 
| *deletinf* | *DELETION_OF_FEATURE_SECTION* | will be blank (`-`) | model (reference) positions that correspond to the segment of the feature that is deleted in the sequence | [deleted features examples](#example-delftr) | 
| *dupregin* | *DUPLICATE_REGIONS* | *N* sets of sequence coordinates, in pairs, each pair is two hits that overlap in model coordinates, *N* will be a factor of 2 | *N* model (reference) coordinates, one for each of the hits in the coverage determination stage that correspond to each set of sequence coordinates 1 to *N* | [duplicate regions example](#example-dupregin) | 
| *discontn* | *DISCONTINUOUS_SIMILARITY* | *N* sets of sequence coordinates, one for each hit in the coverage determination stage | *N* model (reference) coordinates, one for each of the hits in the coverage determination stage that correspond to each set of sequence coordinates 1 to *N* | [discontinuous similarity example](#example-discontn) | 
| *indfstrn* | *INDEFINITE_STRAND* | sequence coordinates of the best hit on the **opposite** strand from the overall best hit for this sequence in the coverage determination stage | model (reference) coordinates for the hit pertaining to the sequence coordinates in `sequence coords` | [indefinite strand example](#example-indfstrn) | 
| *lowcovrg* | *LOW_COVERAGE* | one or more set of sequence coordinates that are **not** covered by any hit to the model on the top-scoring strand in the coverage determination stage | will be blank (`-`) | [low coverage example](#example-lowcovrg) | 

--- 

# <a name="toy"></a>`toy50` toy model used in the examples of alert messages below

The toy50 model is a toy example used to illustrate many of the
problems with sequences that VADR can detect using simple examples on
this page. The toy50 model is 50 nucleotides long and includes 1 CDS
feature from positions 11 to 31 with the name (product) of `protein
one` . That CDS is composed of two adjacent mature peptides: `protein
one mp1` from positions 11 to 22 and `protein one mp2` from positions
23 to 28. The final 3 nucleotides of the CDS, 29 to 31, are the stop
codon. The `model info` file for the toy50 model is shown below.

```
MODEL toy50 cmfile:"toy50.cm" group:"toy" length:"50" subgroup:"A" blastdb:"toy50.protein.fa"
FEATURE toy50 type:"CDS" coords:"11..31:+" parent_idx_str:"GBNULL" gene:"one" product:"protein one"
FEATURE toy50 type:"mat_peptide" coords:"11..22:+" parent_idx_str:"0" product:"protein one mp1"
FEATURE toy50 type:"mat_peptide" coords:"23..28:+" parent_idx_str:"0" product:"protein one mp2"
```

The reference sequence for the toy50 model is shown below, as a
Stockholm format *alignment* file (even though it has one sequence)
with special markup in the form of `#=GC` columns to show where the
CDS and mature peptide features are, as well as the sequence position
information:

```
# STOCKHOLM 1.0

toy50              GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC CDS1.11..31:+ ..........123123123123123123stp...................
#=GC MP1.11..22:+  ..........123123123123............................
#=GC MP2.23..28:+  ......................123123......................
#=GC COLX.         00000000011111111112222222222333333333344444444445
#=GC COL.X         12345678901234567890123456789012345678901234567890
//
```

Stockholm format is described in more detail at
https://en.wikipedia.org/wiki/Stockholm_format and
http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output
formats") 

--- 

# <a name="examples"></a>Examples of alerts and their corresponding `.alt` output lines

## <a name="example-frameshift"></a>Frameshift example

| relevant alert codes | corresponding alert descriptions | notes
|----------------------|----------------------------------|---------
| *fsthicft*, *fsthicfi* | *POSSIBLE_FRAMESHIFT_HIGH_CONF*| only possible if `--glsearch` *is not* used |
| *fstlocft*, *fstlocfi* | *POSSIBLE_FRAMESHIFT_LOW_CONF* | only possible if `--glsearch` *is not* used |
| *fstukcft*, *fstukcfi* | *POSSIBLE_FRAMESHIFT*          | only possible if `--glsearch` *is* used |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-frameshift.sh
  ```



  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-frameshift/va-example-frameshift.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr          ftr              ftr  alert           alert                               seq  seq       mdl  mdl  alert 
#idx   name       model  type         name             idx  code      fail  description                      coords  len    coords  len  detail
#----  ---------  -----  -----------  ---------------  ---  --------  ----  -----------------------------  --------  ---  --------  ---  ------
1.1.2  TOY50-FS1  toy50  CDS          protein_one        1  fsthicfi  yes   POSSIBLE_FRAMESHIFT_HIGH_CONF  13..25:+   13  14..23:+   10  high confidence possible frameshift in CDS (frame restored before end) [cause:insert,S:13..17(5),M:13; restore:delete,S:25,M:22..23(2); frame:1(3)1; length:3:(13):8; shifted_avgpp:0.825; exp_avgpp:0.892;]
```

  **Alignment of `TOY50-FS1` sequence to the toy50 model:** The output file 
  `va-example-frameshift/va-example-frameshift.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be
  helpful in understanding the frameshift, as explained more below.
  The `#=GC RF` line shows the toy50 reference model sequence.
  The `#=GR PP` line indicates confidence estimates for each aligned
  nucleotide as explained more [here](#pp).
  This alignment is only output when the `--keep` or
  `--out_stk` options are used with `v-annotate.pl`. 

```
TOY50-FS1         -AAATCACCGATGcccccGTGATCGC--TACCATAAATGAGCATTCTACGTGCAT
#=GR TOY50-FS1 PP .**********987777789***998..59*************************
#=GC SS_cons      :::::::::::::.....:::::::::::::::::::::::::::::::::::::
#=GC RF           GAAATCACCGATG.....GTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      0000000001111.....1111112222222222333333333344444444445
#=GC RFCOL.X      1234567890123.....4567890123456789012345678901234567890
```

  **How to interpret this alert based on the above output**: 
  As reported in the `.alt` file and `.alt.list` files as shown above, a
  possible frameshift exists in the CDS named `protein one` in the
  sequence named `TOY50-FS1` which matches best to the model named
  `toy50`. The `.alt` file contains details on the frameshift. The
  frameshifted region is sequence positions 13 to 25 (`seq coords:
  13..25:+` in the `.alt` file) which is aligned to the reference model
  positions 14 to 23 (`mdl coords: 14..23:+`).

  The `alert detail` field provides further information:
  the indels that "cause" the frameshifted region are an insertion of length 5 of nucleotides
  13 to 17 after model position 13 (`cause:insert,S:13..17(5),M:13;`).
  The mutation that restores the frame is a deletion of length 2 *after* nucleotide 25 corresponding to model
  positions 22 and 23 (`restore:delete:S:25,M:22..23(2);`). 
  The CDS starts in frame 1 for 3 nt, shifts into frame 3 for the frameshifted region
  for 13 nt and restores to frame 1 for 8 nt before the end of the CDS (`frame:1(3)1` and
  `length:3:(13):8`). The [table below](#framelengths) show how to interpret different `frame` and `length`
  strings. The frame that the CDS starts in is defined as the
  'expected frame' (this will always be 1 unless the CDS is truncated
  at the 5' end). 

  This frameshift is a high confidence
  frameshift in that the average posterior probability of the aligned
  nucleotides in the frameshifted region is `0.825`
  (`shifted_frame:0.825`) and of the expected region prior to (5' of) the
  frameshift is 0.892 (`exp_frame:0.892`) both of which exceed the
  threshold for high confidence (`0.8` by default). (The PP of both
  the shifted and expected regions must exceed the high confidence threshold for a frameshift
  to be defined as high confidence.) Other possible
  frameshifts with lower posterior probability values will be reported
  with the `POSSIBLE_FRAMESHIFT_LOW_CONF` error. If the `--glsearch`
  option is used with `v-annotate.pl`, as is recommended with
  SARS-CoV-2 analysis, posterior probability values are
  not calculated and so all frameshifts are reported with the
  `POSSIBLE_FRAMESHIFT` error. 

  A separate alignment file showing the CDS features that include
  possible frameshifts can be optionally output from `v-annotate.pl`
  using the `--out_fsstk` option. An example excerpt from such an
  alignment file for this possible frameshift is below (see
  `va-example-frameshift/va-exammle-frameshift.vadr.toy50.CDS.1.1.frameshift.stk`). The
  `#=GR PP` shows an estimate of the posterior probability of each
  aligned nucleotide as explained more [here](#pp).  The `#=GR CS`
  line shows the implied frame of each aligned nucleotide and have `i`
  for inserted nucleotides and `d` for deleted reference
  positions. The `#=GC RF` line shows the reference model sequence.

```
TOY50-FS1         ATGCCCCCGTGATCGC--TACCATAA
#=GR TOY50-FS1 PP *987777789***998..59******
#=GR TOY50-FS1 CS 111iiiii33333333dd11111111
#=GC SS_cons      ::::::::::::::::::::::::::
#=GC RF           ATG.....GTGATCGCTTTACCATAA
#=GC RFCOLX.      111.....111111222222222233
#=GC RFCOL.X      123.....456789012345678901
```

### <a name="framelengths"></a> Explanation of `frame` and `length` strings in POSSIBLE_FRAMESHIFT_* alert detail 

| frame string | ......length_string...... | 5' truncated? | 3' truncated? | frame of shifted region | explanation    | 
|--------------|---------------|---------------|---------------|-------------------------|----------------|
| `1(2)`         | `10:(50)`       | no            | no            | 2 | frame 1 for 10 nt, frame 2 for 50 nt |
| `<1(2)`        | `10:(50)`       | yes           | no            | 2 | frame 1 for 10 nt, frame 2 for 50 nt |
| `<1(2)>`       | `10:(50)`       | yes           | yes           | 2 | frame 1 for 10 nt, frame 2 for 50 nt |
| `1(23)131>`    | `10:(50:30)20:9:31` | no        | yes           | 2 then 3 | frame 1 for 10 nt, frame 2 for 50 nt, frame 3 for 30 nt, frame 1 for 20 nt, frame 3 for 9 nt, frame 1 for 31 nt |

---

## <a name="example-insert"></a>Large insertion example

| relevant alert codes | corresponding alert descriptions |
|----------------------|----------------------------------|
| *insertnn*, *insertnp* | *INSERTION_OF_NT* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-insert.sh
  ```
  
  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-insert/va-example-insert.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq              ftr   ftr          ftr  alert           alert                 seq  seq       mdl  mdl  alert 
#idx   name      model  type  name         idx  code      fail  description        coords  len    coords  len  detail
#----  --------  -----  ----  -----------  ---  --------  ----  ---------------  --------  ---  --------  ---  ------
1.1.1  TOY50-I1  toy50  CDS   protein_one    1  insertnn  no    INSERTION_OF_NT  23..28:+    6  23..23:+    1  too large of an insertion in nucleotide-based alignment of CDS feature [6>2]
```

  **Alignment of `TOY50-I1` sequence to the toy50 model:** The output file 
  `va-example-insert/va-example-insert.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be
  helpful in understanding the insertion, as explained more below.
  The `#=GC RF` line shows the toy50 reference model sequence.
  The `#=GR PP` line indicates confidence estimates for each aligned
  nucleotide as explained more [here](#pp).
  This alignment is only output when the `--keep` or
  `--out_stk` options are used with `v-annotate.pl`. 

```
TOY50-I1         -AAATCACCGATGGTGATCGCTTggggggTACCATAAATGAGCATTCTACGTGCAT
#=GR TOY50-I1 PP .********************9866666689*************************
#=GC SS_cons     :::::::::::::::::::::::......:::::::::::::::::::::::::::
#=GC RF          GAAATCACCGATGGTGATCGCTT......TACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.     00000000011111111112222......222222333333333344444444445
#=GC RFCOL.X     12345678901234567890123......456789012345678901234567890
//
```

  **How to interpret this alert based on the above output**: As
  reported in the `.alt` file shown above, in the CDS feature with
  name `protein one`, the nucleotides 23 to 28, length 6, (`seq
  coords:23..28:+`, `seq len:6`) on the + strand insert after
  reference model position 23 (`mdl coords:23..23:+`).  This length
  exceeds the minimum allowed length of 2 (set with the
  `v-annotate.pl` option `--nmaxins 2` option for purposes of this
  example, `alert detail:[6>2]`). You can see the `gggggg` insertion after model position 23
  in the above alignment.

  The `insertnn` alert is detected in the nucleotide alignment stage.
  A similar `insertnp` alert can be reported for insertions detected in the
  protein validation stage, and often you will see both alerts for the
  same insertion. 

---
## <a name="example-delete"></a>Large deletion example

| relevant alert codes | corresponding alert descriptions |
|----------------------|----------------------------------|
| *deletinn*, *deletinp* | *DELETION_OF_NT* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-delete.sh
  ```
  
  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-delete/va-example-delete.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq              ftr   ftr          ftr  alert           alert                seq  seq       mdl  mdl  alert 
#idx   name      model  type  name         idx  code      fail  description       coords  len    coords  len  detail
#----  --------  -----  ----  -----------  ---  --------  ----  --------------  --------  ---  --------  ---  ------
1.1.1  TOY50-D1  toy50  CDS   protein_one    1  deletinn  no    DELETION_OF_NT  15..15:+    1  17..19:+    3  too large of a deletion in nucleotide-based alignment of CDS feature [3>2]
```

  **Alignment of `TOY50-D1` sequence to the toy50 model:** The output file 
  `va-example-delete/va-example-delete.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be
  helpful in understanding the deletion, as explained more below.
  The `#=GC RF` line shows the toy50 reference model sequence.
  The `#=GR PP` line indicates confidence estimates for each aligned
  nucleotide as explained more [here](#pp).
  This alignment is only output when the `--keep` or
  `--out_stk` options are used with `v-annotate.pl`. 

```
TOY50-D1         -AAATCACCGATGGTG---GCTTTACCATAAATGAGCATTCTACGTGCAT
#=GR TOY50-D1 PP .***********9987...89*****************************
#=GC SS_cons     ::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF          GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.     00000000011111111112222222222333333333344444444445
#=GC RFCOL.X     12345678901234567890123456789012345678901234567890
//
```

  **How to interpret this alert based on the above output**: 
  As reported in the `.alt` file shown above,
  in the CDS feature with name `protein one`, the reference model 
  positions 17 to 19, length 3 (`mdl coords:17..19:+`, `mdl len:3`) are deleted in the sequence `TOY50-D1`. The
  deletion occurs after position 15 (`seq coords:15..15:+`) in the sequence.
  The deletion length of 3 exceeds the minimum allowed length of 2 (set with the
  `v-annotate.pl` option `--nmaxdel 2` option for purposes of this
  example, `alert detail:[3>2]`). You can see the three deleted positions (`---`) in model
  RF positions 17 to 19 in the above alignment.

  The `deletinn` alert is detected in the nucleotide alignment stage.
  A similar `deletnp` alert can be reported for deletions detected in the
  protein validation stage, and often you will see both alerts for the
  same deletion. 

---
## <a name="example-start"></a>Mutated start codon example

| relevant alert code  | corresponding alert description  |
|----------------------|----------------------------------|
| *mutstart* | *MUTATION_AT_START* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-start.sh
  ```
  
  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-start/va-example-start.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq              ftr          ftr          ftr  alert           alert                  seq  seq       mdl  mdl  alert 
#idx   name      model  type         name         idx  code      fail  description         coords  len    coords  len  detail
#----  --------  -----  -----------  -----------  ---  --------  ----  ----------------  --------  ---  --------  ---  ------
1.1.1  TOY50-S1  toy50  CDS          protein_one   1  mutstart  yes   MUTATION_AT_START  10..12:+    3  11..13:+    3  expected start codon could not be identified [ATT]
```

  **Alignment of `TOY50-S1` sequence to the toy50 model:** The output
  file `va-example-start/va-example-start.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, allows
  you to see the alignment of the predicted start codon and
  surrounding sequence.  The `#=GC RF` line
  shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                           vvv
TOY50-S1         -AAATCACCGATTGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GR TOY50-S1 PP .*************************************************
#=GC SS_cons     ::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF          GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.     00000000011111111112222222222333333333344444444445
#=GC RFCOL.X     12345678901234567890123456789012345678901234567890
```


  **How to interpret this alert based on the above output**: The first
  three nucleotides of any CDS feature that is not truncated on the 5'
  end are checked to see if they are a valid start codon, and if not,
  the *mutstart* alert is reported. For this specific example, as
  reported in the `.alt` file shown above, the CDS start codon is from
  model (reference) positions 11 to 13 (`mdl coords:11..13:+`), and the first 3 nucleotides of
  the predicted CDS are positions 10 to 12 (`seq coords:10..12:+`). The predicted invalid
  `ATT` start codon (`alert detail:[ATT]`) can be seen in the above alignment, marked by the
  `vvv` characters at the top of the alignment. (These `vvv`
  characters have been added here, and do not exist in the actual
  output alignment file.)

---
## <a name="example-stop"></a>Examples of stop codon problems

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *mutendcd*, *mutendns*, *mutendex* | *MUTATION_AT_END*    |
| *cdsstopn*, *cdsstopp*             | *CDS_HAS_STOP_CODON* |
| *unexleng*                         | *UNEXPECTED_LENGTH*  |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-stop.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-stop/va-example-stop.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr          ftr  alert           alert                     seq  seq       mdl  mdl  alert 
#idx   name       model  type  name         idx  code      fail  description            coords  len    coords  len  detail
#----  ---------  -----  ----  -----------  ---  --------  ----  -------------------  --------  ---  --------  ---  ------
1.1.1  TOY50-SP1  toy50  CDS   protein_one    1  mutendcd  yes   MUTATION_AT_END      28..30:+    3  29..31:+    3  expected stop codon could not be identified, predicted CDS stop by homology is invalid [TAC]
1.1.2  TOY50-SP1  toy50  CDS   protein_one    1  mutendns  yes   MUTATION_AT_END             -    -         -    -  expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted start codon [-]
#
2.1.1  TOY50-SP2  toy50  CDS   protein_one    1  mutendcd  yes   MUTATION_AT_END      27..29:+    3  28..30:+    3  expected stop codon could not be identified, predicted CDS stop by homology is invalid [ATA]
2.1.2  TOY50-SP2  toy50  CDS   protein_one    1  mutendex  yes   MUTATION_AT_END      31..33:+    3  33..35:+    3  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position [TGA]
2.1.3  TOY50-SP2  toy50  CDS   protein_one    1  unexleng  yes   UNEXPECTED_LENGTH    10..29:+   20  11..30:+   20  length of complete coding (CDS or mat_peptide) feature is not a multiple of 3 [20]
#
3.1.1  TOY50-SP3  toy50  CDS   protein_one    1  mutendcd  yes   MUTATION_AT_END      28..30:+    3  29..31:+    3  expected stop codon could not be identified, predicted CDS stop by homology is invalid [TAC]
3.1.2  TOY50-SP3  toy50  CDS   protein_one    1  cdsstopn  yes   CDS_HAS_STOP_CODON   22..24:+    3  23..25:+    3  in-frame stop codon exists 5' of stop position predicted by homology to reference [TAA, shifted S:6,M:6]
```

  **Alignment of `TOY50-SP1`, `TOY50-SP2` and `TOY50-SP3` sequences to the toy50 model:** The output
  file `va-example-stop/va-example-stop.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the stop-codon related alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                              vvv
TOY50-SP1         -AAATCACCGATGGTGATCGCTTTACCATACATGAGCATTCTACGTGCAT
#=GR TOY50-SP1 PP .*************************************************
TOY50-SP2         -AAATCACCGATGGTGATCGCTTTACCATA-CTGAGCATTCTACGTGCAT
#=GR TOY50-SP2 PP .***************************96.69*****************
TOY50-SP3         -AAATCACCGATGGTGATCGCTTAACCATACATGAGCATTCTACGTGCAT
#=GR TOY50-SP3 PP .*************************************************
#=GC SS_cons      ::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      00000000011111111112222222222333333333344444444445
#=GC RFCOL.X      12345678901234567890123456789012345678901234567890
                                        ^^^
```

  **How to interpret this alert based on the above output**: 
  Several checks are made on CDS features related to stop codons and these three sequences demonstrate
  failures of some of these checks. 

  * `TOY50-SP1` sequence: 

  The *mutendcd* alert with description *MUTATION_AT_END* is reported
  because it has an invalid stop codon `TAC` (`alert detail:[TAC]`) at positions 28 to 30
  (`seq coords:28..30:+`) which are aligned to model positions 29 to 31 (`mdl coords:29..31:+`, marked by `vvv`
  characters in the alignment above), which are the final three
  positions of the predicted `protein one` CDS feature. The `mutendns`
  alert with description *MUTATION_AT_END* is also reported because no
  valid in-frame stop codon (in the same frame as the predicted start
  codon) exists in the remainder of the sequence.

  * `TOY50-SP2` sequence: 

  The *mutendcd* alert with description *MUTATION_AT_END* is
  reported because it has an invalid stop codon `ATA` (`alert detail:[ATA]`) at positions 27 to 29 (`seq coords:27..29:+`) which are aligned to model positions 28
  to 30 (`mdl coords:28..30:+`), which are the final 3 positions of the predicted `protein one` 
  CDS feature. Unlike `TOY50-SP1`, there is an in-frame valid `TGA` stop codon 
  3' of the expected stop position (`alert detail:[TGA]`), at sequence positions 31 to 33 (`seq coords:31..33:+`) which 
  align to model positions 33 to 35 (`mdl coords:33..35:+`), which cause a `mutendex` alert with
  description *MUTATION_AT_END* to be reported. Finally, because the predicted CDS feature
  is length 20 (`alert detail:[20]`), from sequence positions 10 to 29 (`seq coords:10..29:+`) aligned to model positions 
  11 to 30 (`mdl coords:11..30:+`), and 20 is not a multiple of 3, the *unexleng* alert with description
  *UNEXPECTED_LENGTH* is reported.

  * `TOY50-SP3` sequence: 

  Like `TOY50-SP1`, the *mutendcd* alert with description *MUTATION_AT_END* is
  reported because it has an invalid stop
  codon `TAC` (`alert detail:[TAC]`) at positions 28 to 30 (`seq coords:28..30:+`) which are aligned to model positions 29
  to 31 (`mdl coords:29..31:+`), which are the final three positions of the predicted `protein one` 
  CDS feature. Unlike the other two sequences, there is an in-frame valid `TAA` 
  stop codon 5' of the expected stop position (`alert detail:[TAA`), at sequence positions 22 to 24 (`seq coords:22..24:+`)
  which align to model positions 23 to 25 (`mdl coords:23..25:+`, marked by `^^^` characters at the 
  bottom of the alignment above), which causes a `cdsstopn` alert with 
  description *CDS_HAS_STOP_CODON*. The early stop codon shifts the stop codon 6 positions in the
  sequence and the reference relative to the expected stop position (`alert detail:shifted S:6,M:6`).

  The `cdsstopn` alert is detected based on the nucleotide alignment.
  A similar `cdsstopp` alert with description *CDS_HAS_STOP_CODON* is
  reported when early in-frame stop codons are detected in the protein validation
  stage by blastx. Often you will see both alerts for the same early stop codon,
  but sometimes you will only see one or the other. 

---
## <a name="example-indf5"></a>Examples of indefinite annotation at the start of a sequence or a feature

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *indf5gap*, *indf5lcc*, indf5pst*, *indf5lcn* (not shown), *indf5plg* (not shown) | *INDEFINITE_ANNOTATION_START* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-indefstart.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-indefstart/va-example-indefstart.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr          ftr  alert           alert                             seq  seq       mdl  mdl  alert 
#idx   name       model  type  name         idx  code      fail  description                    coords  len    coords  len  detail
#----  ---------  -----  ----  -----------  ---  --------  ----  ---------------------------  --------  ---  --------  ---  ------
1.1.4  TOY50-IS1  toy50  CDS   protein_one    1  indf5gap  yes   INDEFINITE_ANNOTATION_START  10..10:+    1  11..11:+    1  alignment to homology model is a gap at 5' boundary [-]
#
2.1.1  TOY50-IS2  toy50  CDS   protein_one    1  indf5lcc  no    INDEFINITE_ANNOTATION_START    8..8:+    1  11..11:+    1  alignment to homology model has low confidence at 5' boundary for feature that is or matches a CDS [0.80<0.90]
```

  **Alignment of `TOY50-IS1` and `TOY50-IS2` sequences to the toy50 model:** The output
  file `va-example-indefstart/va-example-indefstart.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                            v
TOY50-IS1         -AAATCACCG-TGGTGATCGCTTTACCATAAATGAGCAT-----------
#=GR TOY50-IS1 PP .*******98.89**************************...........
TOY50-IS2         -AAATCAC--ATGGTGATCGCTTTACCATAAATGAGCAT-----------
#=GR TOY50-IS2 PP .*****98..89***************************...........
#=GC SS_cons      ::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      00000000011111111112222222222333333333344444444445
#=GC RFCOL.X      12345678901234567890123456789012345678901234567890
```

  **How to interpret these alerts based on the above output**: 
  Several checks are made on the first position of all features, and these
  two sequences demonstrate failures due to some of these checks. 

  * `TOY50-IS1` sequence: 

  The *indf5gap* alert with description *INDEFINITE_ANNOTATION_START*
  is reported because there is a gap at the start position (5'
  boundary) of the CDS named `protein one`. The first *non-gap*
  nucleotide in the predicted CDS is position 10 of the sequence (`seq
  coords:10..10:+`). The 5' boundary gap is position 11 in the model
  (`mdl coords:11..11:+`).

  * `TOY50-IS2` sequence: 

  The *indf5lcc* alert with description *INDEFINITE_ANNOTATION_START*
  is reported because the posterior probability of the aligned nucleotide
  at the start position (5' boundary) of the CDS named `protein one` is too low. 
  In this example the value is 0.8, and the minimum value to not report an 
  alert is 0.9 (`alert detail:[0.80<0.90]`). The minimum is actually
  usually 0.8, but it was changed to 0.9 for demonstration purposes in this example
  with the option `--indefann 0.9`. The nucleotide aligned at the 5' boundary is 
  sequence position 10 and it aligns to position 11 in the model.
  Posterior probabilities are explained more [here](#pp). 

  A similar `indf5lcn` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding
  features cause a sequence to fail or not.

  * Related `indf5pst` and `indf5plg` alerts (not shown in example sequences above):

  If the protein validation stage for a sequence results in too short of blastx
  alignment of the predicted CDS to the reference protein that does not 
  extend close enough to the 5' boundary of the predicted CDS, the *indf5pst*
  alert will be reported. This alert is not possible to generate with the *toy50*
  toy model because the  CDS in this example is too short for constructive use with `blastx`,
  but a fabricated example `.alt` file output line is shown below, 
  in which the blastx alignment  
  extended only to sequence position 19, leaving positions 10 to 18 of the predicted CDS (`seq coords:10..18:+`), which align to 
  model positions 11 to 19 (`mdl coords:11..19:+`), uncovered by the blastx alignment. This length difference of 9
  exceeds the maximum allowed difference of 5 (`alert detail:[9>5]`), so the alert is reported.
  A similar `indf5plg` alert exists for when the blastx alignment extends *longer* than the 
  predicted nucleotide alignment on the 5' end, but no example is shown here. 

```
#      seq              ftr   ftr          ftr  alert           alert                             seq  seq       mdl  mdl  alert 
#idx   name      model  type  name         idx  code      fail  description                    coords  len    coords  len  detail
#----  --------- -----  ----  -----------  ---  --------  ----  ---------------------------  --------  ---  --------  ---  ------
4.2.1  TOY50.F1  toy50  CDS   protein_one    1  indf5pst  yes   INDEFINITE_ANNOTATION_START  10..18:+    9  11..19:+    9  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint [9>5]
```

---
## <a name="example-indf3"></a>Examples of indefinite annotation at end of a feature

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *indf3gap*, *indf3lcc*, indf3pst*, *indf3lcn* (not shown), *indf3plg* (not shown) | *INDEFINITE_ANNOTATION_END* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-indefend.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-indefend/va-example-indefend.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr          ftr  alert           alert                             seq  seq       mdl  mdl  alert 
#idx   name       model  type  name         idx  code      fail  description                    coords  len    coords  len  detail
#----  ---------  -----  ----  -----------  ---  --------  ----  ---------------------------  --------  ---  --------  ---  ------
1.1.4  TOY50-IE1  toy50  CDS   protein_one    1  indf3gap  yes   INDEFINITE_ANNOTATION_END    28..28:+    1  31..31:+    1  alignment to homology model is a gap at 3' boundary [-]
#
2.1.1  TOY50-IE2  toy50  CDS   protein_one    1  indf3lcc  no    INDEFINITE_ANNOTATION_END    30..30:+    1  31..31:+    1  alignment to homology model has low confidence at 3' boundary for feature that is or matches a CDS [0.70<0.80]
```

  **Alignment of `TOY50-IE1` and `TOY50-IE2` sequences to the toy50 model:** The output
  file `va-example-indefend/va-example-indefend.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                                v
TOY50-IE1         -AAATCACCGATGGTGATCGCTTTACCAT--ATGAGCAT-----------
#=GR TOY50-IE1 PP .**************************98..59******...........
TOY50-IE2         -AAATCACCGATGGTGATCCCTCTAGCATAGA-GAGCAT-----------
#=GR TOY50-IE2 PP .***************************9876.9*****...........
#=GC SS_cons      ::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      00000000011111111112222222222333333333344444444445
#=GC RFCOL.X      12345678901234567890123456789012345678901234567890
```

  **How to interpret these alerts based on the above output**: 
  Several checks are made on the first position of all features, and these
  two sequences demonstrate failures to some of these checks. 

  * `TOY50-IE1` sequence: 

  The *indf3gap* alert with description *INDEFINITE_ANNOTATION_END*
  is reported because there is a gap at the end position (3'
  boundary) of the CDS named `protein one`. The final
  *non-gap* nucleotide in the predicted CDS is position 28 (`seq coords:28..28:+`) of the
  sequence. The 3' boundary gap is position 31 in the model (`mdl coords:31..31:+`).

  * `TOY50-IE2` sequence: 

  The *indf3lcc* alert with description *INDEFINITE_ANNOTATION_END*
  is reported because the posterior probability of the aligned nucleotide
  at the end position (3' boundary) of the CDS named `protein one` is too low. 
  In this example the value is 0.7, and the minimum value to not report an 
  alert is 0.8 (`alert detail:[0.70<0.80]`). The nucleotide aligned at the 3' boundary is 
  sequence position 30 (`seq coords:30..30:+`) and it aligns to position 31 (`mdl coords:31..31:+`) in the model.
  Posterior probabilities are explained more [here](#pp). 

  A similar `indf3lcn` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding
  features cause a sequence to fail or not.

  * Related `indf3pst` and `indf3plg` alerts (not shown in example sequences above):

  If the protein validation stage for a sequence results in too short of blastx
  alignment of the predicted CDS to the reference protein that does not 
  extend close enough to the 3' boundary of the predicted CDS, the *indf3pst*
  alert will be reported. This alert is not possible to generate with the *toy50*
  toy model because the CDS in this example is too short for constructive use with `blastx`. 
  But a fabricated example `.alt` file output line is shown below, 
  in which the blastx alignment  
  extended only to sequence position 21, leaving positions 22 to 30 of the predicted CDS (`seq coords:22..30:+`), which align to 
  model positions 23 to 31 (`mdl coords:23..31:+`), uncovered by the blastx alignment. This length difference of 9
  exceeds the maximum allowed difference of 8 (`alert detail:[9>8]`), so the alert is reported.
  A similar `indf3plg` alert exists for when the blastx alignment extends *longer* than the 
  predicted nucleotide alignment on the 3' end, but no example is shown here. 

```
#      seq              ftr   ftr          ftr  alert           alert                             seq  seq       mdl  mdl  alert 
#idx   name      model  type  name         idx  code      fail  description                    coords  len    coords  len  detail
#----  --------- -----  ----  -----------  ---  --------  ----  ---------------------------  --------  ---  --------  ---  ------
7.2.1  TOY50-F2  toy50  CDS   protein_one    1  indf3pst  yes   INDEFINITE_ANNOTATION_END    22..30:+    9  23..31:+    9  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint [9>8]
```

---
## <a name="example-ambg5"></a>Examples of ambiguous nucleotides at the start of a sequence or a feature

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *ambgnt5s* | *AMBIGUITY_AT_START* | 
| *ambgnt5c* | *AMBIGUITY_AT_CDS_START* | 
| *ambgnt5f* | *AMBIGUITY_AT_FEATURE_START* | 

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-ambigstart.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-ambigstart/va-example-ambigstart.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr          ftr              ftr  alert           alert                            seq  seq       mdl  mdl  alert 
#idx   name       model  type         name             idx  code      fail  description                   coords  len    coords  len  detail
#----  ---------  -----  -----------  ---------------  ---  --------  ----  --------------------------  --------  ---  --------  ---  ------
1.1.1  TOY50-AS1  toy50  CDS          protein_one        1  ambgnt5c  no    AMBIGUITY_AT_CDS_START      10..13:+    4  11..14:+    4  first nucleotide of CDS is an ambiguous nucleotide [-]
1.2.1  TOY50-AS1  toy50  mat_peptide  protein_one_mp1    2  ambgnt5f  no    AMBIGUITY_AT_FEATURE_START  10..13:+    4  11..14:+    4  first nucleotide of non-CDS feature is an ambiguous nucleotide [-]
#
2.1.1  TOY50-AS2  toy50  -            -                  -  ambgnt5s  no    AMBIGUITY_AT_START           1..13:+   13   2..14:+   13  first nucleotide of the sequence is an ambiguous nucleotide [-]
2.2.1  TOY50-AS2  toy50  CDS          protein_one        1  ambgnt5c  no    AMBIGUITY_AT_CDS_START      10..13:+    4  11..14:+    4  first nucleotide of CDS is an ambiguous nucleotide [-]
2.3.1  TOY50-AS2  toy50  mat_peptide  protein_one_mp1    2  ambgnt5f  no    AMBIGUITY_AT_FEATURE_START  10..13:+    4  11..14:+    4  first nucleotide of non-CDS feature is an ambiguous nucleotide [-]
```

  **Alignment of `TOY50-AS1` and `TOY50-AS2` sequences to the toy50 model:** The output
  file `va-example-ambigstart/va-example-ambigstart.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                            vvvv
TOY50-AS1         -AAATCACCGNNNNTGATCGCTTTACCATAAATGAGCAT-----------
#=GR TOY50-AS1 PP .**************************************...........
TOY50-AS2         -NNNNNNNNNNNNNTGATCGCTTTACCATAAATGAGCAT-----------
#=GR TOY50-AS2 PP .**************************************...........
#=GC SS_cons      ::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      00000000011111111112222222222333333333344444444445
#=GC RFCOL.X      12345678901234567890123456789012345678901234567890
```
  **How to interpret these alerts based on the above output**: 
  If the first position of each sequence or predicted feature 
  is an N, an alert is reported.

  * `TOY50-AS1` sequence: 

  The *ambgnt5c* alert with description *AMBIGUITY_AT_CDS_START*
  is reported because the predicted CDS named `protein one` begins
  with 4 consecutive Ns from sequence positions 10
  to 13 (`seq coords:10..13:+`) that align to the model (reference) positions 11 to 14 (`mdl coords:11..14:+`). These positions
  are marked with `vvvv` in the alignment above. This alert is specific
  to CDS features.

  Similarly, the *ambgnt5f* alert with description *AMBIGUITY_AT_FEATURE_START*
  is reported because the predicted mat_peptide named `protein one mp1` begins
  with 4 consecutive Ns from sequence positions 10
  to 13 (`seq coords:10..13:+`) that align to the model (reference) positions 11 to 14 (`mdl coords:11..14:+`). These positions
  are marked with `vvvv` in the alignment above. This alert is only reported
  for non-CDS features.

  * `TOY50-AS2` sequence: 

  This sequence is similar to `TOY50-AS1` except the stretch of Ns begins at 
  position 1 and extends to position 13 (`seq coords:1..13:+`). Because the first 4 nucleotides of
  the CDS and mat_peptide are Ns like in `TOY50-AS1`, the same *ambgnt5c* and
  *ambgnt5f* alerts are reported, but now with an additional *ambgnt5s* alert
  with description *AMBIGUITY_AT_START* because the beginning of the sequence begins
  with a consecutive string of 13 Ns.

  This sequence is similar to `TOY50-AS1` except the stretch of Ns
  begins at position 1 and extends to position 13 
  position of the sequence causing the *ambgnt5s*
  alert with description *AMBIGUITY_AT_START* to be reported because the sequence ends 
  with a stretch of 13 consecutive Ns. Because the first 4 nucleotides of
  the CDS and mat_peptide are Ns like in `TOY50-AS1`, the same *ambgnt5c* and
  *ambgnt5f* alerts are reported as well.

---

## <a name="example-ambg3"></a>Examples of ambiguous nucleotides at end of a sequence or a feature

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *ambgnt3s* | *AMBIGUITY_AT_END* | 
| *ambgnt3c* | *AMBIGUITY_AT_CDS_END* | 
| *ambgnt3f* | *AMBIGUITY_AT_FEATURE_END* | 

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-ambigend.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-ambigend/va-example-ambigend.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr          ftr  alert           alert                      seq  seq       mdl  mdl  alert 
#idx   name       model  type  name         idx  code      fail  description             coords  len    coords  len  detail
#----  ---------  -----  ----  -----------  ---  --------  ----  --------------------  --------  ---  --------  ---  ------
1.1.1  TOY50-AE1  toy50  CDS   protein_one    1  ambgnt3c  no    AMBIGUITY_AT_CDS_END  29..30:+    2  30..31:+    2  final nucleotide of CDS is an ambiguous nucleotide [-]
#
2.1.1  TOY50-AE2  toy50  -     -              -  ambgnt3s  no    AMBIGUITY_AT_END      29..38:+   10  30..39:+   10  final nucleotide of the sequence is an ambiguous nucleotide [-]
2.2.1  TOY50-AE2  toy50  CDS   protein_one    1  ambgnt3c  no    AMBIGUITY_AT_CDS_END  29..30:+    2  30..31:+    2  final nucleotide of CDS is an ambiguous nucleotide [-]
```

  **Alignment of `TOY50-AE1` and `TOY50-AE2` sequences to the toy50 model:** The output
  file `va-example-ambigend/va-example-ambigend.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                               vv
TOY50-AE1         -AAATCACCGATGGTGATCGCTTTACCATNNATGAGCAT-----------
#=GR TOY50-AE1 PP .**************************************...........
TOY50-AE2         -AAATCACCGATGGTGATCGCTTTACCATNNNNNNNNNN-----------
#=GR TOY50-AE2 PP .**************************************...........
#=GC SS_cons      ::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      00000000011111111112222222222333333333344444444445
#=GC RFCOL.X      12345678901234567890123456789012345678901234567890
```
  **How to interpret these alerts based on the above output**: 
  If the first position of each sequence or predicted feature 
  is an N, an alert is reported.

  * `TOY50-AE1` sequence: 

  The *ambgnt3c* alert with description *AMBIGUITY_AT_CDS_END*
  is reported because the predicted CDS named `protein one` ends
  with 2 consecutive Ns from sequence positions 29
  to 30 (`seq coords:29..30:+`) that align to the model (reference) positions 30 to 31 (`mdl coords:30..31:+`). These positions
  are marked with `vv` in the alignment above. This alert is specific
  to CDS features.

  A similar *ambgnt3f* alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here. Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus
  non-coding features cause a sequence to fail or not.

  * `TOY50-AE2` sequence: 

  This sequence is similar to `TOY50-AE1` except the stretch of Ns
  begins at position 29 and continues to the position 38 (`seq coords:29..38:+`) which is the final
  position of the sequence causing the *ambgnt3s*
  alert with description *AMBIGUITY_AT_END* to be reported because the sequence ends 
  with a stretch of 10 consecutive Ns. Because the
  final 2 nucleotides of the CDS are Ns like in `TOY50-AS1`, the same
  *ambgnt3c* alert is reported. 

---
## <a name="example-ambgcd"></a>Examples of ambiguous nucleotides in start/stop codons that start/stop with canonical nt

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *ambgcd5c* | *AMBIGUITY_IN_START_CODON* |
| *ambgcd3c* | *AMBIGUITY_IN_STOP_CODON* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-ambigcodon.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-ambigend/va-example-ambigend.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** *****You may have to scroll to the right to see the `alert detail` field.*****

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr          ftr  alert           alert                      seq  seq       mdl  mdl  alert 
#idx   name       model  type  name         idx  code      fail  description             coords  len    coords  len  detail
#----  ---------  -----  ----  -----------  ---  --------  ----  --------------------  --------  ---  --------  ---  ------
1.1.1  TOY50-AC1  toy50  CDS   protein_one    1  ambgcd5c  no    AMBIGUITY_IN_START_CODON  10..12:+    3  11..13:+    3  5' complete CDS starts with canonical nt but includes ambiguous nt in its start codon [ANG]
#
2.1.1  TOY50-AC2  toy50  CDS   protein_one    1  ambgcd3c  no    AMBIGUITY_IN_STOP_CODON   28..30:+    3  29..31:+    3  3' complete CDS ends with canonical nt but includes ambiguous nt in its stop codon [TNA]
```

  **Alignment of `TOY50-AC1` and `TOY50-AC2` sequences to the toy50 model:** The output
  file `va-example-ambigcodon/va-example-ambigcodon.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                            vvv               
TOY50-AC1         -AAATCACCGANGNTGATCGCTTTACCATAAATGAGCAT-----------
#=GR TOY50-AC1 PP .**************************************...........
TOY50-AC2         -AAATCACCGATGGTGATCGCTTTACCATNAATGAGCAT-----------
#=GR TOY50-AC2 PP .**************************************...........
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      00000000011111111112222222222333333333344444444445
#=GC RFCOL.X      12345678901234567890123456789012345678901234567890
                                              ^^^
```
  **How to interpret these alerts based on the above output**: 
  If the first position of each sequence or predicted feature 
  is an N, an alert is reported.

  * `TOY50-AC1` sequence: 

  The *ambgcd5c* alert with description *AMBIGUITY_IN_START_CODON*
  is reported because the predicted CDS named `protein one` begins
  with a start codon `ANG` that includes an ambiguous nucleotide *but*
  starts with a canonical nt at sequence positions 10 to 12 (`seq coords:10..12:+`)
  that align to the model (reference) positions 11 to 13 (`mdl coords:11..13:+`). These positions
  are marked with `vv` in the alignment above. This alert is specific
  to CDS features that are not truncated at the 5' end.

  * `TOY50-AC2` sequence: 

  The *ambgcd3c* alert with description *AMBIGUITY_IN_STOP_CODON*
  is reported because the predicted CDS named `protein one` ends
  with a stop codon `TNA` that includes an ambiguous nucleotide *but*
  ends with a canonical nt at sequence positions 28 to 30 (`seq coords:28..30:+`)
  that align to the model (reference) positions 29 to 31 (`mdl coords:29..31:+`). These positions
  are marked with `^^` in the alignment above. This alert is specific
  to CDS features that are not truncated at the 3' end.

---
## <a name="example-pep"></a>Examples of mature peptide-specific alerts

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *pepadjcy* | *PEPTIDE_ADJACENCY_PROBLEM* | 
| *peptrans* | *PEPTIDE_TRANSLATION_PROBLEM* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-matpep.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-matpep/va-example-matpep.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr          ftr              ftr  alert           alert                             seq  seq       mdl  mdl  alert 
#idx   name       model  type         name             idx  code      fail  description                    coords  len    coords  len  detail
#----  ---------  -----  -----------  ---------------  ---  --------  ----  ---------------------------  --------  ---  --------  ---  ------
1.1.1  TOY50-MP1  toy50  mat_peptide  protein_one_mp1    2  pepadjcy  yes   PEPTIDE_ADJACENCY_PROBLEM    22..24:+    3  22..23:+    2  predictions of two mat_peptides expected to be adjacent are not adjacent [-]
#
2.1.1  TOY50-MP2  toy50  CDS          protein_one        1  mutstart  yes   MUTATION_AT_START            10..12:+    3  11..13:+    3  expected start codon could not be identified [GTG]
2.2.1  TOY50-MP2  toy50  mat_peptide  protein_one_mp1    2  peptrans  yes   PEPTIDE_TRANSLATION_PROBLEM         -    -         -    -  mat_peptide may not be translated because its parent CDS has a problem [-]
2.3.1  TOY50-MP2  toy50  mat_peptide  protein_one_mp2    3  peptrans  yes   PEPTIDE_TRANSLATION_PROBLEM         -    -         -    -  mat_peptide may not be translated because its parent CDS has a problem [-]
```

  **Alignment of `TOY50-MP1` and `TOY50-MP2` sequences to the toy50 model:** The output
  file `va-example-matpep/va-example-matpep.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                       v   v
TOY50-MP1         -AAATCACCGATGGTGATCGCTgggTTACCATAAATGAGCAT-----------
#=GR TOY50-MP1 PP .******************99766589***************...........
TOY50-MP2         -AAATCACCGGTGGTGATCGCT...TTACCATAAATGAGCAT-----------
#=GR TOY50-MP2 PP .*********************...*****************...........
#=GC SS_cons      ::::::::::::::::::::::...::::::::::::::::::::::::::::
#=GC RF           GAAATCACCGATGGTGATCGCT...TTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      0000000001111111111222...2222222333333333344444444445
#=GC RFCOL.X      1234567890123456789012...3456789012345678901234567890
                            ^^^
```

  **How to interpret these alerts based on the above output**: 

  * `TOY50-MP1` sequence: 

  When the protein product of a CDS is cleaved into multiple mature
  peptides the mature peptide features typically are adjacent to each
  other like they are for the toy example `toy50` (details
  [here](#toy)) in which `protein one mp1` which spans nucleotide
  positions 11 to 22, is adjacent to `protein one mp2` which spans
  positions 23 to 38. For these situations, `v-annotate.pl` checks that the predicted
  mature peptide features in each sequence are adjacent, and if
  not the *pepadjcy* alert is reported.

  The *pepadjcy* alert with description *PEPTIDE_ADJACENCY_PROBLEM* is
  reported for the `TOY50-MP1` sequence because the predicted CDS
  named `protein one mp1` ends at position 21 and `protein one mp2`
  begins at position 25, so the positions 22 to 24 lie between the two
  mature peptide predictions (`seq coords:22..24:+`). 
  The model position that ends
  `protein one mp1` is 22, and that begins `protein one mp2` is 23 (`mdl coords:22..23:+`).

  * `TOY50-MP2` sequence: 

  When a parent CDS of a mature peptide has a fatal alert, a
  *peptrans* alert with description *PEPTIDE_TRANSLATION_PROBLEM* is
  reported for all of the predicted child mature peptides. The 
  `TOY50-MP2` sequence shows an example, which has a `mutstart` alert
  for the CDS `protein one` which causes *peptrans* alerts for 
  the two mature peptides.

---
## <a name="example-lowsim5"></a>Examples of low similarity to the model at the start of a sequence or a feature

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *lowsim5c*, *lowsim5n* (not shown), *lowsim5l* (not shown) | *LOW_FEATURE_SIMILARITY_START* |
| *lowsim5s*                         | *LOW_SIMILARITY_START* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-lowsimstart.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-lowsimstart/va-example-lowsimstart.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq                ftr   ftr          ftr  alert           alert                             seq  seq       mdl  mdl  alert 
#idx   name        model  type  name         idx  code      fail  description                    coords  len    coords  len  detail
#----  ----------  -----  ----  -----------  ---  --------  ----  ----------------------------  -------  ---  --------  ---  ------
1.2.3  TOY50-LSS1  toy50  CDS   protein_one    1  lowsim5c  no    LOW_FEATURE_SIMILARITY_START  1..13:+   13   5..13:+    9  region overlapping annotated feature that is or matches a CDS at 5' end of sequence lacks significant similarity [7 nt overlap b/t low similarity region of length 13 (1..13) and annotated feature (7..31)]
#
2.1.2  TOY50-LSS2  toy50  -     -              -  lowsim5s  yes   LOW_SIMILARITY_START          1..10:+   10   3..10:+    8  significant similarity not detected at 5' end of the sequence [low similarity region of length 10]
```

  **Alignment of `TOY50-LSS1` and `TOY50-LSS2` sequences to the toy50 model:** The output
  file `va-example-lowsimstart/va-example-lowsimstart.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                       vvvvvv  vvvvvvv                 
TOY50-LSS1         ----TTGTAG..GTTcgacGTGATCGCTTTACCATAAATGAGCAT-----------
#=GR TOY50-LSS1 PP ....998766..542223389************************...........
TOY50-LSS2         --GTTTAGTGgcATG....GTGATCGCTTTACCATAAATGAGCAT-----------
#=GR TOY50-LSS2 PP ..**9998775589*....**************************...........
#=GC SS_cons       ::::::::::..:::....:::::::::::::::::::::::::::::::::::::
#=GC RF            GAAATCACCG..ATG....GTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.       0000000001..111....1111112222222222333333333344444444445
#=GC RFCOL.X       1234567890..123....4567890123456789012345678901234567890
```

  **Explanation of `lowsim5c` alert**: The 13 nucleotides from
  positions 1 to 13 in the sequence `TOY50-LSS1` (`seq coords:1..13:+`) are at the 5' end
  of the sequence and overlap by 7 nucleotides with the predicted CDS `protein one` but are not similar to the reference
  model. These 13 nucleotides align to reference model positions 5 to 13 (`mdl coords:5..13:+`) and 
  include an insertion of 4 nucleotides after reference position 13. The 13 nucleotides
  are marked by `v` characters in the alignment above.          

  Two similar alerts, `lowsim5f` and `lowsim5l`, exist for non-coding (non-CDS and
  non-mature peptide) features, but no examples are shown here. `lowsim5f` is for relatively short regions 
  of low similarity, while `lowsim5l` is for longer regions (30nt or more). Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding features.

  **Explanation of `lowsim5s` alert**: The first ten nucleotides 
  in the sequence `TOY50-LSS2` are not similar to the reference model (`seq coords:1..10:+`)
  and do not overlap with any predicted features. These ten nucleotides align to reference model positions 3 to 10,
  (`mdl coords:3..10:+`) and include two inserted nucleotides after reference position 10.

  Regions of low similarity are detected in the *coverage
  determination* stage, as regions that are not covered by local
  alignment *hits* between the sequence and the model, not based on
  the global alignment determined in the alignment stage. In the
  alignment above, note that the nucleotides in the low similarity
  regions at the beginning of each sequence do not match well to the
  nucleotides in the reference model (`#=GC RF` line).

---
## <a name="example-lowsim3"></a>Examples of low similarity to the model at the end of a sequence or a feature

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *lowsim3c*, *lowsim3n* (not shown), *lowsim3l* (not shown) | *LOW_FEATURE_SIMILARITY_END* |
| *lowsim3s*                         | *LOW_SIMILARITY_END* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-lowsimend.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-lowsimend/va-example-lowsimend.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq                ftr          ftr              ftr  alert           alert                             seq  seq       mdl  mdl  alert 
#idx   name        model  type         name             idx  code      fail  description                    coords  len    coords  len  detail
#----  ----------  -----  -----------  ---------------  ---  --------  ----  ----------------------------  -------  ---  --------  ---  ------
1.2.2  TOY50-LSE1  toy50  CDS          protein_one        1  lowsim3c  no    LOW_FEATURE_SIMILARITY_END    26..37:+   12  26..31:+    6  region overlapping annotated feature that is or matches a CDS at 3' end of sequence lacks significant similarity [12 nt overlap b/t low similarity region of length 12 (26..37) and annotated feature (10..37)]
#
2.1.2  TOY50-LSE2  toy50  -            -                  -  lowsim3s  yes   LOW_SIMILARITY_END            40..53:+   14  35..46:+   12  significant similarity not detected at 3' end of the sequence [low similarity region of length 14]
```

  **Alignment of `TOY50-LSE1` and `TOY50-LSE2` sequences to the toy50 model:** The output
  file `va-example-lowsimend/va-example-lowsimend.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                            vvvvvvvvvvvvv                 
TOY50-LSE1         -AAATCACCGATGGTGATCGCTTTACgtcccgtCTTAA----........---------------
#=GR TOY50-LSE1 PP .***********************975544444689**...........................
TOY50-LSE2         -AAATCACCGATGGTGATCGCTTTAC.......CATAAATGAcgatacacGAACTGCACGA----
#=GR TOY50-LSE2 PP .*************************.......********922222222334466799**....
#=GC SS_cons       ::::::::::::::::::::::::::.......:::::::::........:::::::::::::::
#=GC RF            GAAATCACCGATGGTGATCGCTTTAC.......CATAAATGA........GCATTCTACGTGCAT
#=GC RFCOLX.       00000000011111111112222222.......222333333........333344444444445
#=GC RFCOL.X       12345678901234567890123456.......789012345........678901234567890
                                                            ^^^^^^^^^^^^^^^^^^^^
```

  **Explanation of `lowsim3c` alert**: The 12 nucleotides from
  positions 26 to 37 (`seq coords:26..37:+`) in the sequence `TOY50-LSE1` are the 3' end
  of the sequence and overlap with the predicted CDS `protein one` but are not similar to the reference
  model. These 12 nucleotides align to reference model positions 26 to 31
  (`mdl coords:26..31:+`) including an insertion of length 7 after model position 26.
  The 12 nucleotides
  are marked by `v` characters in the alignment above.
  
  Two similar alerts, `lowsim3f` and `lowsim3l`, exist for non-coding (non-CDS and
  non-mature peptide) features, but no examples are shown here. `lowsim3f` is for relatively short regions 
  of low similarity, while `lowsim3l` is for longer regions (30nt or more). Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding features.

  **Explanation of `lowsim3s` alert**: The final 14 nucleotides 
  from positions 40 to 53 (`seq coords:40..53:+`) in the sequence `TOY50-LSE2` are not similar to the reference model 
  and do not overlap with any predicted features. These 14 nucleotides align to reference model positions 35 to 46 (`mdl coords:35..46:+`),
  and include 8 inserted nucleotides after reference position 35.
  The alignment above shows the region of low similarity marked by `^` characters.

  Regions of low similarity are detected in the *coverage determination* 
  stage, as regions that are not covered by local alignment *hits* between the 
  sequence and the model, not based on the global alignment determined in the 
  alignment stage.

---
## <a name="example-lowsimi"></a>Examples of low similarity to the model in a region internal to a sequence or feature

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *lowsimic*, *lowsimin* (not shown), *lowsimil* (not shown) | *LOW_FEATURE_SIMILARITY* |
| *lowsimis*                         | *LOW_SIMILARITY* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-lowsimint.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-lowsimint/va-example-lowsimint.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq                ftr   ftr          ftr  alert           alert                         seq  seq       mdl  mdl  alert 
#idx   name        model  type  name         idx  code      fail  description                coords  len    coords  len  detail
#----  ----------  -----  ----  -----------  ---  --------  ----  ----------------------  ---------  ---  --------  ---  ------
1.2.3  TOY50-LSI1  toy50  CDS   protein_one  1    lowsimic  no    LOW_FEATURE_SIMILARITY  32..107:+   76  27..28:+    2  region overlapping annotated feature that is or matches a CDS lacks significant similarity [76 nt overlap b/t low similarity region of length 76 (32..107) and annotated feature (10..114)]
```

  **Alignment of `TOY50-LSI1` sequence to the toy50 model:** The output
  file `va-example-lowsimint/va-example-lowsimint.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                                  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
TOY50-LSI1         -AAATCACCGATGGTGATCGCTTTACCaaagcagtacaggcacatgacaaagcagtacaggca-catgacaaagcagtacaggcacatgacaaagcagtacaggcacatgacaTAAATGAGCATTCTACGTGCAT
#=GR TOY50-LSI1 PP .**********************9875222222222222222222222222222222222222.222222222222222222222222222222222222222222222222267899*****************
#=GC SS_cons       :::::::::::::::::::::::::::....................................:.................................................::::::::::::::::::::::
#=GC RF            GAAATCACCGATGGTGATCGCTTTACC....................................A.................................................TAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.       000000000111111111122222222....................................2.................................................2333333333344444444445
#=GC RFCOL.X       123456789012345678901234567....................................8.................................................9012345678901234567890
//
```

  **Explanation of `lowsimic` alert**: The 76 nucleotides from
  positions 32 to 107 (`seq coords:32..107:+`) in the sequence `TOY50-LSI1` occur in the middle
  of the predicted CDS `protein one` but are not similar to the reference
  model. These 76 nucleotides align to reference model positions 27 to 28 (`mdl coords:27..28:+`),
  but are all actually inserted after those two reference positions as marked with 
  `v` characters in the alignment shown above.

  A similar `lowsimin` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding

  Two similar alerts, `lowsimif` and `lowsimil`, exist for non-coding (non-CDS and
  non-mature peptide) features, but no examples are shown here. `lowsimif` is for relatively short regions 
  of low similarity, while `lowsimil` is for longer regions (30nt or more). Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding features.

  Regions of low similarity are detected in the *coverage determination* 
  stage, as regions that are not covered by local alignment *hits* between the 
  sequence and the model, not based on the global alignment determined in the 
  alignment stage. That is why the region marked with `v` characters above
  which corresponds to the region of low similarity does not include the 
  full insertions after reference positions 27 and 28. Enough similarity exists in 
  the beginning of the insertion after position 27 and at the end of the insertion 
  after position 28 for a local alignment hit to include those regions. 

---
## <a name="example-delftr"></a>Example of a deleted feature

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *deletins*, *deletina* (not shown) | *DELETION_OF_FEATURE* |
| *deletinf* (not shown)             | *DELETION_OF_FEATURE_SECTION* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-delftr.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-delftr/va-example-delftr.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr   ftr  alert           alert                   seq  seq       mdl  mdl  alert 
#idx   name       model  type  name  idx  code      fail  description          coords  len    coords  len  detail
#----  ---------  -----  ----  ----  ---  --------  ----  -------------------  ------  ---  --------  ---  ------
1.1.1  TOY50-DF1  toy50  -     -       -  deletins  yes   DELETION_OF_FEATURE       -    -  23..28:+    6  internal deletion of a complete feature [mat_peptide feature number 2: protein one mp2]
```

  **Alignment of `TOY50-DF1` sequence to the toy50 model:** The output
  file `va-example-delftr/va-example-delftr.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                        vvvvvv
TOY50-DF1         -AAATCACCGATGGTGATCGC--------AAATGAGCATTCTACGTGCAT
#=GR TOY50-DF1 PP .*****************997........79*******************
#=GC SS_cons      ::::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      00000000011111111112222222222333333333344444444445
#=GC RFCOL.X      12345678901234567890123456789012345678901234567890
```

  **Explanation of `deletins` alert**: The mature peptide feature named `protein one mp2` 
  is deleted in the sequence `TOY50-DF1`, as indicated in the `alert detail` field of the
  `.alt` file above. That mature peptide spans model positions 23 to 28 (`model coords:23..28:+`) which can be
  seen as gaps in the above alignment, and are marked with positions `v`.

  A similar `deletina` alert will be reported if the feature that is deleted
  has the `is_deletable` field set to `1` in the input `.minfo` file.
  No example of this alert is shown here.

A similar `deletinf` alert exists for *multi-segment* features for which
  one or more but not all segments are deleted. The associated value for this
  alert in the `alert desc` field is `DELETION_OF_FEATURE_SECTION`. 
  No example of this alert is shown here.
  
---
## <a name="example-dupregin"></a>Example of duplicate regions

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *dupregin* | *DUPLICATE_REGION* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-dupregin.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-dupregin/va-example-dupregin.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr   ftr  alert           alert                           seq  seq              mdl  mdl  alert 
#idx   name       model  type  name  idx  code      fail  description                  coords  len           coords  len  detail
#----  ---------  -----  ----  ----  ---  --------  ----  -----------------  ----------------  ---  ---------------  ---  ------
1.1.1  TOY50-DR1  toy50  -     -       -  dupregin  yes   DUPLICATE_REGIONS  1..49:+,50..86:+   86  2..50:+,7..43:+   86  similarity to a model region occurs more than once [7..43 (len 37>=20) hits 1 (39.5 bits) and 2 (27.0 bits)]
```

  **Alignment of `TOY50-DR1` sequence to the toy50 model:** The output
  file `va-example-dupregin/va-example-dupregin.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                                                    |positions 50 to 86 in the sequence-|
                                                                    vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
TOY50-DR1         -AAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCataccgatggtgatcgctttaccataaatgagcattctA-
#=GR TOY50-DR1 PP .********************************************987444444444444444444444444444444444444446.
#=GC SS_cons      ::::::::::::::::::::::::::::::::::::::::::::::::......................................::
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGC......................................AT
#=GC RFCOLX.      000000000111111111122222222223333333333444444444......................................45
#=GC RFCOL.X      123456789012345678901234567890123456789012345678......................................90
                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                        |--positions 7 to 43 in the model---|
```

  **Explanation of `dupregin` alert**: There are two hits in the
  coverage determination stage, which overlap in model coordinates by
  more than 20 positions which cause the `dupregin` alert to be reported.
  Specifically, hit 1 is from
  positions 1 to 49 in the sequence (`seq coords:1..49:+`) and 2 to 50 in the model (`mdl coords:2..50:+`) and hit 2
  is from positions 50 to 86 in the sequence (`seq coords:50..86:+`) and positions 7 to 43 in
  the model (`mdl coords;7..43)`).  The hits overlap by
  37 model positions from 7 to 43, which is above the minimum length
  threshold for reporting this alert of 20 (`alert detail`: `7..43 (len
  37 >= 20)`). The first hit has a bit score of 39.5 bits and the
  second hit of 27.0 bits (`alert detail`: `hits 1 (39.5 bits) and 2
  (27.0 bits)`).

  In the alignment of the sequence above sequence positions 50 to 86
  are marked marked by `v` characters at the top of the alignment,
  which match identically to positions 7 to 43 in the reference model,
  marked by `^` at the bottom of the alignment.

---
## <a name="example-discontn"></a>Example of discontinuous similarity

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *discontn* | *DISCONTINUOUS_SIMILARITY |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-discontn.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-discontn/va-example-discontn.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr   ftr  alert           alert                                  seq  seq               mdl  mdl  alert 
#idx   name       model  type  name  idx  code      fail  description                         coords  len            coords  len  detail
#----  ---------  -----  ----  ----  ---  --------  ----  ------------------------  ----------------  ---  ----------------  ---  ------
1.1.1  TOY50-DS1  toy50  -     -       -  discontn  yes   DISCONTINUOUS_SIMILARITY  1..26:+,27..49:+   49  25..50:+,2..24:+   49  not all hits are in the same order in the sequence and the homology model [seq order: 1,2, mdl order: 2,1]
```

  **Alignment of `TOY50-DS1` sequence to the toy50 model:** The output
  file `va-example-discontn/va-example-discontn.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                                                    |-27 to 49 in the seq-|
                                                                    vvvvvvvvvvvvvvvvvvvvvvv
TOY50-DS1         ------------------------ACCATAAATGAGCATTCTACGTGCAtaaatcaccgatggtgatcgcttT
#=GR TOY50-DS1 PP ........................**********************98666666667777777777777777*
#=GC SS_cons      :::::::::::::::::::::::::::::::::::::::::::::::::.......................:
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCA.......................T
#=GC RFCOLX.      0000000001111111111222222222233333333334444444444.......................5
#=GC RFCOL.X      1234567890123456789012345678901234567890123456789.......................0
                   ^^^^^^^^^^^^^^^^^^^^^^^
                   |--2-24 in the model--|
```

  **Explanation of `discontn` alert**: For the `TOY50-DS1` sequence,
  there are two hits in the coverage determination stage and they are
  not in the same order in sequence and model coordinates. The
  positional information for these hits are reported in the `seq
  coords`, `mdl coords` and `alert detail` fields of the `.alt` file
  above. Specifically, hit 1 is from positions 1 to 26 in the
  sequence (`seq coords:1..26:+`) and 25 to 50 in the model (`mdl coords:25..50:+`) and hit 2 is from positions 27 to
  49 in the sequence (`seq coords:27..49:+`) and positions 2 to 24 in the model (`mdl coords:2..24:+`).
  So hit 1 comes before hit 2 in the
  sequence, but hit 2 comes before hit 1 in the model. The order of
  the hits in the sequence and model is reported in the `alert detail`
  field: `seq order: 1,2, mdl order: 2,1`.

  In the alignment of the sequence above, sequence positions 27 to 49
  are marked marked by `v` characters at the top of the alignment.
  These sequence positions match identically to positions 2 to 24 in
  the reference model, marked by `^` at the bottom of the alignment,
  showing why the second hit is to the end of the sequence and the
  beginning of the model. The first hit pertains to the the beginning
  of the sequence which matches well to the end of the model (the
  region with `*` values in the `PP` line).

---
## <a name="example-indfstrn"></a>Example of indefinite strand

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *indfstrn*, *indfstrp* (not shown) | *INDEFINITE_STRAND* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-indfstrn.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-indfstrn/va-example-indfstrn.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr   ftr  alert           alert                   seq  seq       mdl  mdl  alert 
#idx   name       model  type  name  idx  code      fail  description          coords  len    coords  len  detail
#----  ---------  -----  ----  ----  ---  --------  ----  -----------------  --------  ---  --------  ---  ------
1.1.2  TOY50-IS1  toy50  -     -       -  indfstrn  yes   INDEFINITE_STRAND  49..25:-   25  26..50:+   25  significant similarity detected on both strands [score:12.7>12.0]

```

  **Alignment of `TOY50-IS1` sequence to the toy50 model:** The output
  file `va-example-indfstrn/va-example-indfstrn.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                                           |  25 to 49 in the seq  |
                                           vvvvvvvvvvvvvvvvvvvvvvvvv
TOY50-IS1         -AAATCACCGATGGTGATCGCTTTAATGCAcgtagAATGCTCATTTATGG-----
#=GR TOY50-IS1 PP .**********************98544332222389999**********.....
#=GC SS_cons      ::::::::::::::::::::::::::::::.....::::::::::::::::::::
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATA.....AATGAGCATTCTACGTGCAT
#=GC RFCOLX.      000000000111111111122222222223.....33333333344444444445
#=GC RFCOL.X      123456789012345678901234567890.....12345678901234567890
                                           ^^^^^     ^^^^^^^^^^^^^^^^^^^^
                                           |   26 to 50 in the model    |
```

  **Explanation of `indfstrn` alert**: For the `TOY50-IS1` sequence
  there are two hits in the coverage determination stage and they are
  not to the same strand. The
  positional information for these hits are reported in the `seq
  coords`, `mdl coords` and `alert detail` fields of the `.alt` file
  above. Specifically, the top-scoring hit is to the `+` strand,
  but the second hit with a score of 12.7 bits is on the `-` strand.
  The second hit occurs from positions to 49 to 25 on the negative strand (`seq coords:49..25:-`) to
  positions 26 to 50 in the reference model (`mdl coords:26..50:+`)
  The second hit therefore
  conflicts with the top-scoring hit and so it is the one that is
  described in the `seq coords` and `mdl coords` fields. Normally, a
  *indfstrn* alert is only reported if a hit on the opposite strand
  from the top hit exceeds 25 bits, but to construct this example, the
  minimum of 25 bits was lowered to 12 bits with the `--indefstr 12`
  flag to `v-annotate.pl`.

  In the alignment of the sequence above, sequence positions 25 to
  49 are marked marked by `v` characters at the top of the alignment,
  *after reverse complementing* these match identically to positions 26
  to 50 in the reference model, marked by `^` at the bottom of the
  alignment.

  The `indftrsn` alert is detected in the coverage determination stage
  which is based on the nucleotide sequence. A similar `indfstrp`
  alert can be reported when blastx hits occur on both strands as
  detected in the protein validation stage, but is not
  shonw. Sometimes you will see both alerts for the same sequence.

---
## <a name="example-lowcovrg"></a>Example of low coverage

| relevant alert codes  | corresponding alert descriptions  |
|----------------------|------------------------------------|
| *lowcovrg* | *LOW_COVERAGE* |

  **Instructions to reproduce this example and create the files discussed below:**
  ```
  > sh $VADRSCRIPTSDIR/documentation/alert-files/example-lowcovrg.sh
  ```

  **The above command will run `v-annotate.pl` and create many output files in a newly created directory. The output shown in the box below is from the `.alt` output file (`va-example-lowcovrg/va-example-lowcovrg.vadr.alt`) . Similar information with the `seq` and `mdl coords` fields can be found in the `.alt.list` output file and in the detailed error report `.tsv` file generated by the GenBank submission portal.** ***You may have to scroll to the right to see the `alert detail` field.***

#### Example lines from `.alt` file:

```
#      seq               ftr   ftr   ftr  alert           alert                                  seq  seq               mdl  mdl  alert 
#idx   name       model  type  name  idx  code      fail  description                         coords  len            coords  len  detail
#----  ---------  -----  ----  ----  ---  --------  ----  ------------------------  ----------------  ---  ----------------  ---  ------
1.1.1  TOY50-LC1  toy50  -     -       -  lowcovrg  yes   LOW_COVERAGE          1..25:+   25       -    -  low sequence fraction with significant similarity to homology model [0.545<0.900]
```

  **Alignment of `TOY50-LC1` sequence to the toy50 model:** The output
  file `va-example-lowcovrg/va-example-lowcovrg.vadr.toy50.align.stk`.
  includes the alignment shown below. Looking at this alignment, or an
  alignment of the sequence generated by a different program, can be 
  helpful in understanding the alerts.
  The `#=GC RF` line shows the toy50 reference model sequence.  The `#=GR PP` line
  indicates confidence estimates for each aligned nucleotide as
  explained more [here](#pp).  This alignment is only output when the
  `--keep` or `--out_stk` options are used with `v-annotate.pl`.

```
                  |1 to 25 in the sequence|
                  vvvvvvvvvvvvvvvvvvvvvvvvv
TOY50-LC1         GAtacatcgatcatcaatccgggacaAATCACCGATGGTGATCGCTTTACCATAA-------------------
#=GR TOY50-LC1 PP *64444444444444333333333336899*************************...................
#=GC SS_cons      ::........................::::::::::::::::::::::::::::::::::::::::::::::::
#=GC RF           GA........................AATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.      00........................000000011111111112222222222333333333344444444445
#=GC RFCOL.X      12........................345678901234567890123456789012345678901234567890
```

  **Explanation of `lowcovrg` alert**: For the `TOY50-LC1` sequence
  the nucleotides from 1 to 25 are not part of any hit to the model in
  the coverage determation stage (`seq coords:1..25:+`). These 25 nucleotides make up 45.5\%
  the total length of the sequence, meaning that only 54.5\% of the
  sequence *is covered* by hits to the model (`alert detail:[0.545<0.900]`). The `lowcovrg` alert is reported 
  for any sequence if less than 90\% of the sequence is covered by hits to the model
  (changeable to `<x>` with the `--lowcov <x>` option).

  In the alignment of the sequence above, note that positions 1 to 25, marked
  by `v` characters at the top of the alignment, do not match well to the model.

---
### <a name="toy50"></a>Toy example of a model: TOY50

The TOY50 model is a toy example used to illustrate many of the
problems with sequences that VADR can detect using simple examples on
this page. The TOY50 model is 50 nucleotides long and includes 1 CDS
feature from positions 11 to 31 with the name (product) of `protein
one` . That CDS is composed of two adjacent mature peptides: `protein
one mp1` from positions 11 to 22 and `protein one mp2` from positions
23 to 28. The final 3 nucleotides of the CDS, 29 to 31, are the stop
codon. The `model info` file for the TOY50 model is shown below.

```
MODEL TOY50 cmfile:"toy50.cm" group:"toy" length:"50" subgroup:"A" blastdb:"toy50.protein.fa"
FEATURE TOY50 type:"CDS" coords:"11..31:+" parent_idx_str:"GBNULL" gene:"one" product:"protein one"
FEATURE TOY50 type:"mat_peptide" coords:"11..22:+" parent_idx_str:"0" product:"protein one mp1"
FEATURE TOY50 type:"mat_peptide" coords:"23..28:+" parent_idx_str:"0" product:"protein one mp2"
```

The reference sequence for the TOY50 model is shown below, as a
Stockholm format *alignment* file (even though it has one sequence)
with special markup in the form of `#=GC` columns to show where the
CDS and mature peptide features are, as well as the sequence position
information:

```
# STOCKHOLM 1.0

TOY50              GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC CDS1.11..31:+ ..........123123123123123123stp...................
#=GC MP1.11..22:+  ..........123123123123............................
#=GC MP2.23..28:+  ......................123123......................
#=GC RFCOLX.       00000000011111111112222222222333333333344444444445
#=GC RFCOL.X       12345678901234567890123456789012345678901234567890
//
```

Stockholm format is described in more detail at
https://en.wikipedia.org/wiki/Stockholm_format and
http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output
formats") 

---
### <a name="pp"></a>Posterior probability annotation in VADR output Stockholm alignments

The `v-annotate.pl` script uses a probabilistic model called a
covariance model (CM) to calculate glocal alignments of input
sequences in the alignment stage (unless the `--glsearch` option is
used, in which case a non-probabilistic alignment algorithm is used).
As part of the CM alignment calculation, the posterior probability
(PP) that each aligned nucleotide appears at its assigned alignment
position is calculated. These posterior probabilities are confidence
estimates that each nucleotide is correctly aligned given the
parameters of the CM.  The Stockholm alignments output from
`v-annotate.pl` include annotation on the alignment that indicates
these PP values. An example is below:

```
# STOCKHOLM 1.0
#=GF AU Infernal 1.1.4

toy50-1         GGTATCACAGATGGGATCCGCTGACTCAT-AATGTGTGTTCAaAAGTGCAT
#=GR toy50-1 PP 999*********99888888877766766.567888887775257999999
#=GC RF         GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCT.ACGTGCAT
#=GC RFCOLX.    000000000111111111122222222223333333333444.44444445
#=GC RFCOL.X    123456789012345678901234567890123456789012.34567890

```

The top line beginning with `toy50-1` shows the alignment of the sequence named
`toy50-1`. `-` characters here indicate deletions with respect to the model (gaps).

The `#=GC RF` line shows the reference model sequence. `.` characters in this line indicate
positions that are insertions relative to the model in one or more of the aligned sequences.

The `#=GC RFCOLX.` line shows the tens value of the reference position of each column.
The `#=GC RFCOL.X` line shows the ones value of the reference position of each column.
(E.g. the final column is reference position 50.)

The posterior probability confidence estimates are shown in the `#=GR
toy50-1 PP` line. Characters in PP rows have 12 possible values:

| PP character | meaning | 
|--------------|---------|
|`.`| this position is a gap in the sequence |
|`0`| posterior probability between 0.00 and 0.05 | 
|`1`| posterior probability between 0.05 and 0.15 | 
|`2`| posterior probability between 0.15 and 0.25 | 
|`3`| posterior probability between 0.25 and 0.35 | 
|`4`| posterior probability between 0.35 and 0.45 | 
|`5`| posterior probability between 0.45 and 0.55 | 
|`6`| posterior probability between 0.55 and 0.65 | 
|`7`| posterior probability between 0.65 and 0.75 | 
|`8`| posterior probability between 0.75 and 0.85 | 
|`9`| posterior probability between 0.85 and 0.95 | 
|`*`| posterior probability between 0.95 and 1.00 | 

In the above alignment, the 9 positions 4 to 12 have `*` PP values, indicating that
those positions are very confidently aligned at the correct positions given the parameters of the model. 
As you might expect, these positions nearly exactly match the model nucleotides they are aligned to.

Positions surrounding insertions and deletions tend to have lower PP
values because there is more uncertainty as to where each specific
nucleotide should be placed in the alignment.  For example, reference
position 29 and 30 which are separated by an insertion in the sequence
have PP values of `6` and `5`.

PP values are used by `v-annotate.pl` in two contexts: 

1. To distinguish between low-confidence and high-confidence frameshift regions (*fsthicf** vs *fstlocf** alerts, see this [example](#example-frameshift)).
2. To report alerts when feature boundaries have low confidence (*indf5lc** and *indf3lc** alerts, see example at [5' end](#example-indf5) and [3' end](#example-indf3).

---
#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.


