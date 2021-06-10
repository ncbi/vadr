# <a name="top"></a> Explanations and examples of `v-annotate.pl` detailed alert and error messages

* [Output fields with detailed alert and error messages](#files)
* [Explanation of sequence and model coordinate fields in `.alt` files](#coords)
* [TOY50 toy model used in examples of alert messages](#toy)
* [Example `.alt` output for different alert types](#examples)
* [TOY50 toy model used in examples of alert messages](#toy)

---

## <a name="files"></a> Output files with detailed alert and error messages

`v-annotate.pl` outputs two types of files with detailed alert/error
messages: 

  * `.alt` files
  * `.alt.list` files

Documentation on the format of `.alt` files can be found
[here](formats.md#alt) and for `.alt.list` files can be found
[here](formats.md#altlist).

This page includes examples of many of the different alerts and 
corresponding `.alt` file output [below](#examples).

## <a name="coords"></a> Explanation of sequence and model coordinate fields in `.alt` files

| alert code(s) | alert desc(s) | sequence coords description | sequence coords length constraints | model coords explanation | model coords length constraints | link to example | 
|---------|---------------------|-----------------------------|--------------------|--------------------------|--------------------|---------|
| *fsthicf5*, *fsthicf3*, *fsthicfi*, *fstlocf5*, *fstlocf3*, *fstlocfi*, *fstukcf5*, *fstukcf3*, *fstukcfi* | *POSSIBLE_FRAMESHIFT_HIGH_CONF*,  *POSSIBLE_FRAMESHIFT_LOW_CONF*, *POSSIBLE_FRAMESHIFT* | sequence positions of the frameshifted region | none | model (reference) positions of the frameshifted region, some nucleotides may be inserted **before or after** these positions | none | [frameshift alert example](#example-frameshift) | 
| *insertnn*, *insertnp* | *INSERTION_OF_NT* | sequence positions of inserted nucleotides with respect to the model |  none | model (reference) position after which insertion occurs | always length 1 | [insert alert example](#example-insert) | 
| *deletinn*, *deletinp* | *DELETION_OF_NT*  | sequence position just prior to (5' of) deletion with respect to the model | always length 1 | model (reference) positions that are deleted in sequence | none | [delete alert example](#example-delete) | 
| *mutstart* | *MUTATION_AT_START*  | sequence positions of predicted start codon | length <= 3 | model (reference) positions that align to the predicted start codon | none | [mutstart alert example](#example-mutstart) | 
| *mutendcd* | *MUTATION_AT_END*  | sequence positions of predicted stop codon | length <= 3 | model (reference) positions that align to the predicted stop codon | none | [stop codon alert examples](#example-stop) | 
| *mutendex* | *MUTATION_AT_END*  | sequence positions of 5'-most in-frame stop codon in the CDS, this stop codon will be 3' of expected stop codon position | always length 3 | model (reference) positions that align to stop codon in `sequence coords` | none | [mutend* alert examples](#example-stop) | 
| *mutendns* | *MUTATION_AT_END*  | will be blank (`-`) | N/A | will be blank (`-`) | N/A | [mutend* alert examples](#example-stop) | 
| *unexleng* | *UNEXPECTED_LENGTH* | sequence positions of the predicted CDS, the length of which is not a multiple of 3 | none | model (reference) positions that the predicted CDS align to, some nucleotides may be inserted *before or after* these positions | none | [mutend* alert examples](#example-stop) |
| *cdsstopn* | *CDS_HAS_STOP_CODON* | sequence positions of the 5'-most in-frame stop codon in the CDS, this stop will be 5' of expected stop codong position | always length 3 | model (reference) positions that align to stop codon in `sequence coords` | none | [mutend* alert examples](#example-stop) | 
| *indf5gap* | *INDEFINITE_ANNOTATION_START* | sequence position of first nucleotide aligned 3' of gap that aligns to the feature boundary | always length 1 | model (reference) position of the 5' feature boundary | always length 1 | [indf* alert examples](#example-indf) | 
| *indf5lcc* | *INDEFINITE_ANNOTATION_START* | sequence position of nucleotide aligned at the 5' feature boundary | always length 1 | model (reference) position of the 5' feature boundary | always length 1 | [indf* alert examples](#example-indf) | 
| *indf5pst* | *INDEFINITE_ANNOTATION_START* | sequence positions of the nucleotide alignment of the 5' end of the CDS **not** covered by the protein-based alignment | none | model (reference) positions the sequences positions in `sequence coords` are aligned to in the nucleotide alignment | none | 
| *indf5plg* | *INDEFINITE_ANNOTATION_START* | sequence positions of the protein-based alignment **not** covered by the nucleotide alignment at the 5' end of the CDS | none | model (reference) position of the 5' boundary of the CDS | always length 1 | 
| *indf3gap* | *INDEFINITE_ANNOTATION_END*   | sequence position of final nucleotide aligned 5' of gap that aligns to the feature boundary | always length 1 | model (reference) position of the 3' feature boundary | always length 1 | [indf* alert examples](#example-indf) | 
| *indf3lcc* | *INDEFINITE_ANNOTATION_START* | sequence position of nucleotide aligned at the 3' feature boundary | always length 1 | model (reference) position of the 3' feature boundary | always length 1 | [indf* alert examples](#example-indf) | 
| *indf3pst* | *INDEFINITE_ANNOTATION_START* | sequence positions of the nucleotide alignment of the 3' end of the CDS **not** covered by the protein-based alignment | none | model (reference) positions the sequences positions in `sequence coords` are aligned to in the nucleotide alignment | none | 
| *indf3plg* | *INDEFINITE_ANNOTATION_START* | sequence positions of the protein-based alignment **not** covered by the nucleotide alignment at the 3' end of the CDS | none | model (reference) position of the 3' boundary of the CDS | always length 1 | 
| *ambgnt5f*, ambgnt5c* | *N_AT_FEATURE_START*, *N_AT_CDS_START* | sequence position(s) of stretch of 1 or more consecutive Ns ending at the 3' end of a feature | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none| [ambg3* alert examples](#example-ambg5) | 
| *ambgnt5s* | *N_AT_END* | sequence position(s) of 1 or more consecutive Ns starting at position 1 of the sequence | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none| [ambg5* alert examples](#example-ambg5) | 
| *ambgnt3f*, ambgnt3c* | *N_AT_FEATURE_END*, *N_AT_CDS_END* | sequence position(s) of 1 or more consecutive Ns starting at the predicted 5' end of a feature | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none| [ambg3* alert examples](#example-ambg3) | 
| *ambgnt3s* | *N_AT_END* | sequence position(s) of 1 or more consecutive Ns ending at the final position of the sequence | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none| [ambg3* alert examples](#example-ambg3) | 
| *pepadjcy* | *PEPTIDE_ADJACENCY_PROBLEM* | sequence position(s) of nucleotides inserted between two mature peptide predictions that are expected to be adjacent | none | model (reference) position(s) corresponding to the end of the 5' mature peptide and the start of the 3' mature peptide | always length 2 | [pep* alert examples](#example-pep) | 
| *peptrans* | *PEPTIDE_TRANSLATION_PROBLEM* | will be blank  (`-`) | N/A | will be blank (`-`) | N/A | [pep* alert examples](#example-pep) | 
| *lowsim5c*, *lowsim5n* | *LOW_FEATURE_SIMILARITY_START* | sequence position(s) at 5' end of predicted feature that have low similarity to the reference model | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none | [lowsim5* alert examples](#example-lowsim5) | 
| *lowsim5s* | *LOW_SIMILARITY_START* | sequence position(s) at 5' end of sequence (not overlapping with a feature) that have low similarity to the reference model | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none | [lowsim5* alert examples](#example-lowsim5) | 
| *lowsim3c*, *lowsim3n* | *LOW_FEATURE_SIMILARITY_START* | sequence position(s) at 3' end of predicted feature that have low similarity to the reference model | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none | [lowsim3* alert examples](#example-lowsim5) | 
| *lowsim3s* | *LOW_SIMILARITY_START* | sequence position(s) at 3' end of sequence (not overlapping with a feature) that have low similarity to the reference model | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none | [lowsim3* alert examples](#example-lowsim5) | 
| *lowsimic*, *lowsimin* | *LOW_FEATURE_SIMILARITY* | sequence position(s) internal to a predicted feature (not including first or final position of the feature) that have low similarity to the reference model | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none | [lowsimi* alert examples](#example-lowsimi) | 
| *lowsimis* | *LOW_SIMILARITY* | sequence position(s) internal to a sequence (not including first or final position of the sequence) and not overlapping with a feature that have low similarity to the reference model | none | model (reference) position(s) the sequence position(s) in `sequence coords` are aligned to | none | [lowsimi* alert examples](#example-lowsimi) | 
| *deletins* | *DELETION_OF_FEATURE* | will be blank (`-`) | N/A | model (reference) positions that correspond to the feature that is deleted in the sequence | none | [deleted feature examples*#example-deletin) | 
| *deletinf* | *DELETION_OF_FEATURE_SECTION* | will be blank (`-`) | N/A | model (reference) positions that correspond to the segment of the feature that is deleted in the sequence | none | [deleted feature examples](#example-deletin) | 
| *dupregin* | *DUPLICATE_REGIONS* | *N* sets of sequence coordinates, in pairs, each pair is two hits that overlap in model coordinates, *N* will be a factor of 2 | none | *N* model (reference) coordinates, one for each of the hits in the coverage determination stage that correspond to each set of sequence coordinates 1 to *N* | none | [duplicate regions examples](#example-dupregin) | 
| *discontn* | *DISCONTINUOUS_SIMILARITY* | *N* sets of sequence coordinates, one for each hit in the coverage determination stage | none | *N* model (reference) coordinates, one for each of the hits in the coverage determination stage that correspond to each set of sequence coordinates 1 to *N* | none | [discontinuous similarity example](#example-discontn) | 
| *indfstrn* | *INDEFINITE_STRAND* | sequence coordinates of the best hit on the **opposite** strand from the overall best hit for this sequence in the coverage determination stage | none | model (reference) coordinates for the hit pertaining to the sequence coordinates in `sequence coords` | none | [indefinite strand example](#example-indfstrn) | 
| *lowcovrg* | *LOW_COVERAGE* | one or more set of sequence coordinates that are **not** covered by any hit to the model on the top-scoring strand in the coverage determination stage | none | will be blank (`-`) | N/A | [low coverage example](#example-lowcovrg) | 

### <a name="toy"></a>TOY50 toy model used in examples of alert messages

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
#=GC COLX.         00000000011111111112222222222333333333344444444445
#=GC COL.X         12345678901234567890123456789012345678901234567890
//
```

Stockholm format is described in more detail at
https://en.wikipedia.org/wiki/Stockholm_format and
http://eddylab.org/infernal/Userguide.pdf (section 9: "File and output
formats") 

--- 

## <a name="examples"></a>Examples of different types of alerts and corresponding `.alt` output lines

### <a name="example-frameshift"></a>Example frameshift alert

#### alert codes: *fsthicf5*, *fsthicf3*, *fsthicfi*, *fstlocf5*, *fstlocf3*, *fstlocfi*, *fstukcf5*, *fstukcf3*, *fstukcfi* 

#### corresponding error messages: *POSSIBLE_FRAMESHIFT_HIGH_CONF*,  *POSSIBLE_FRAMESHIFT_LOW_CONF*, *POSSIBLE_FRAMESHIFT*

#### Example line from `.alt` file:

```
#       seq                        ftr          ftr              ftr  alert           alert                               seq     seq       mdl     mdl  alert 
#idx    name            model      type         name             idx  code      fail  desc                             coords  length    coords  length  detail
#-----  --------------  ---------  -----------  ---------------  ---  --------  ----  -----------------------------  --------  ------  --------  ------  ------
7.1.2   ENTOY100A-fs6   ENTOY100A  CDS          protein_one        2  fsthicfi  yes   POSSIBLE_FRAMESHIFT_HIGH_CONF  14..25:+      12  14..22:+       9  high confidence possible frameshift in CDS (internal) [nucleotide alignment of internal sequence positions 14..25 (12 nt, avgpp: 0.890) to model positions 14..22 (9 nt) on + strand are frame 3 (dominant frame is 1); inserts:S:14..18(5),M:13; deletes:S:25,M:21..22(2);]
```

  **Explanation**: a possible frameshift exists in the CDS named
  `protein one` in the sequence named `ENTOY100A-fs6` which matches
  best to the model named `ENTOY100A`. The frameshifted region is
  sequence positions 14 to 25 (`seq coords: 14..25:+`) which is
  aligned to the reference model positions 14 to 22 (`mdl coords:
  14..22:+`) and are in frame 3, while the dominant frame for the CDS
  (frame in which the most nucleotides are in) is frame 1. The indels
  that cause the frameshifted region are an insertion of length 5 of nucleotides
  14 to 18 after model position 13 (`inserts:S:14..18(5),M:13;`) and a
  deletion of length 2 *after* nucleotide 25 corresponding to model
  positions 21 and 22 (`deletes:S:25,M:21..22(2);`).  This frameshift is a high confidence
  frameshift in that the average posterior probability of the aligned
  nucleotides in the frameshifted region is `0.890` which exceeds the
  threshold for high confidence (`0.8` by default). Other possible
  frameshifts with lower posterior probability values will be reported
  with the `POSSIBLE_FRAMESHIFT_LOW_CONF` error. If the `--glsearch`
  option is used with `v-annotate.pl`, as is recommended with
  SARS-CoV-2 analysis, posterior probability values are
  not calculated and so all frameshifts are reported with the
  `POSSIBLE_FRAMESHIFT` error. 

  An alignment file showing the CDS features that include possible
  frameshifts can be optionally output from `v-annotate.pl` using the
  `--out_fsstk` option. An example excerpt from such an alignment file
  for this possible frameshift is below. The `#=GR PP` shows an
  estimate of the posterior probability of each aligned nucleotide,
  (`*` is the highest value, next highest is `9`, then `8`, then `7`,
  etc.). The `#=GR CS` line shows the implied frame of each aligned
  nucleotide and have `i` for inserted nucleotides and `d` for deleted
  reference positions. The `#=GC RF` line shows the reference model
  sequence.

```
ENTOY100A-fs6         ATGCCCCCGTGATCG--TTACCATAA
#=GR ENTOY100A-fs6 PP **9888889*****9..68*******
#=GR ENTOY100A-fs6 CS 111iiiii3333333dd111111111
#=GC RF               atG.....GTGatCGCTTTACCATAA
```

  Alignment files with complete aligned sequences can be output using
  the `--out_stk` or `--out_afa` options. An example excerpt from such
  an alignment file is below. The `#=GC RFCOL*` lines indicate the
  positions of each reference model position (`01` to `50`) in this example.

```
ENTOY100A-fs6          GAAATCACCGATGCCCCCGTGATCG--TTACCATAAATGAGCATTCTACGTGCAT
#=GR ENTOY100A-fs6  PP ************9888889*****9..68**************************
#=GC RF                GAAATCACCGatG.....GTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.           0000000001111.....1111112222222222333333333344444444445
#=GC RFCOL.X           1234567890123.....4567890123456789012345678901234567890
```
---

### <a name="example-insert"></a>Example insert alerts

#### alert codes: *insertnn*, *insertnp*

#### corresponding error messages: *INSERTION_OF_NT*

#### Example line from `.alt` file:

```
#       seq                          ftr          ftr                ftr  alert           alert                seq     seq       mdl     mdl  alert 
#idx    name              model      type         name               idx  code      fail  desc              coords  length    coords  length  detail
#-----  ----------------  ---------  -----------  -----------------  ---  --------  ----  ---------------  -------  ------  --------  ------  ------
19.1.2  ENTOY100A-2-fs18  ENTOY100A  CDS          protein_three        6  insertnn  no    INSERTION_OF_NT  88..93:+      5  86..86:+       1  too large of an insertion in nucleotide-based alignment of CDS feature [nucleotide alignment insert of length 6>2 after reference nucleotide posn 86 on strand +]
19.1.3  ENTOY100A-2-fs18  ENTOY100A  CDS          protein_three        6  insertnp  yes   INSERTION_OF_NT  88..93:+      5  86..86:+       1  too large of an insertion in protein-based alignment [blastx predicted insert of length nucleotide alignment insert of length 6>2 starting after reference amino acid position 5]
```

  **Explanation**: in sequence named `ENTOY100A-2-fs18`, the CDS
  feature with name `protein three`, the nucleotides 88 to 93 (length
  6) on the + strand insert after reference model position 86.  This
  length exceeds the minimum allowed length of 2 (set with the
  `v-annotate.pl` option `--nmaxins 2` option for purposes of this
  example). Both lines of the `.alt` file pertain to the same
  insertion, which is common. The `insertnn` alert is detected during
  the nucleotide alignment stage of the entire sequence. The
  `insertnp` alert is detected in the protein validation stage with
  `blastx`. The `alert detail` field for the `insertnp` alert reports
  the additional information that the insertion occurs after the 5th
  amino acid of `protein three`.

  The alignment of the sequence to the model (`#=GC RF` line) below
  shows the insertion of `TTTTTT` after position 86:

```
ENTOY100A-2-fs18         GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCATCTTGCGGTGCCATACAATGGTAGAAAATTGCCATTCATTTTTTCGTACGTAGCATCA
#=GR ENTOY100A-2-fs18 PP ************************************************************************9752799****999855555579************
#=GC RF                  GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCATCTTGCGGTGCCATACAATGGTAGAA.ATTGCCATTCA......CGTACGTAGCATCA
#=GC RFCOLX..            000000000000000000000000000000000000000000000000000000000000000000000000000.00000000000......00000000000001
#=GC RFCOL.X.            000000000111111111122222222223333333333444444444455555555556666666666777777.77778888888......88899999999990
#=GC RFCOL..X            123456789012345678901234567890123456789012345678901234567890123456789012345.67890123456......78901234567890
```

---
### <a name="example-delete"></a>Example delete alerts

#### alert codes: *deletinn*, *deletinp*

#### corresponding error messages: *DELETION_OF_NT*

#### Example line from `.alt` file:

```
#       seq                          ftr          ftr                ftr  alert           alert                                                 seq     seq                         mdl     mdl  alert 
#idx    name              model      type         name               idx  code      fail  desc                                               coords  length                      coords  length  detail
#-----  ----------------  ---------  -----------  -----------------  ---  --------  ----  -----------------------------  --------------------------  ------  --------------------------  ------  ------
22.4.1  ENTOY100A-2-fs21  ENTOY100A  CDS          protein_four         9  deletinn  no    DELETION_OF_NT                                   86..86:-       1                    89..87:-       2  too large of a deletion in nucleotide-based alignment of CDS feature [nucleotide alignment delete of length 3>2 starting at reference nucleotide posn 88 on strand -]
```

  **Explanation**: in sequence named `ENTOY100A-2-fs21`, for the CDS
  feature with name `protein_four`, the reference model positions 89
  to 87 (length 3) on the negative (-) strand are deleted. This
  deletion occurs 'after' (3' of) sequence position 86 (that is,
  sequence position 85 is the first nucleotide in the CDS that is 3'
  of the deletion, so the the nucleotides 85 and 86 bracket the
  deletion). This length 3 deletion exceeds the minimum allowed of 2
  (set with the `v-annotate.pl` option `--nmaxdel 2` option for purposes of this
  example). 
  Both lines of the `.alt` file pertain to the same
  deletion, which is common. The `deletinn` alert is detected during
  the nucleotide alignment stage of the entire sequence. The
  `deletinp` alert is detected in the protein validation stage with
  `blastx`. The `alert detail` field for the `deletinp` alert reports
  the additional information that the deletioninsertion occurs after the 3rd
  amino acid of `protein four`.

  The alignment of the sequence to the model (`#=GC RF`
  line) below shows the deletion at reference positions 87, 88 and 89.

```
ENTOY100A-2-fs21         GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCATCTTGCGGTGCCATACAATGGTAGAAATT-CCATTCA---ACGTAGCATCA
#=GR ENTOY100A-2-fs21 PP ***************************************************************************998.789998A...8**********
#=GC RF                  GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCATCTTGCGGTGCCATACAATGGTAGAAATTGCCATTCACGTACGTAGCATCA
#=GC RFCOLX..            0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
#=GC RFCOL.X.            0000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
#=GC RFCOL..X            1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
```
---
### <a name="example-mutstart"></a>Example *mutstart* alert

#### alert code: *mutstart*

#### corresponding error messages: *MUTATION_AT_START*

#### Example line from `.alt` file:

```
#      seq                   ftr          ftr              ftr  alert           alert                             seq     seq       mdl     mdl  alert 
#idx   name       model      type         name             idx  code      fail  desc                           coords  length    coords  length  detail
#----  ---------  ---------  -----------  ---------------  ---  --------  ----  ---------------------------  --------  ------  --------  ------  ------
1.1.1  ENTOY100A  ENTOY100A  CDS          protein_one        2  mutstart  yes   MUTATION_AT_START            10..12:+       3  11..13:+       3  expected start codon could not be identified [ATT starting at sequence position 10 (model position 11) on + strand is not a valid start]
```

  **Explanation**: The first three nucleotides of any CDS feature are checked to see if they 
  are a valid start codon, and if not, the *mutstart* alert is reported. For this specific example, the
  CDS starts at model (reference) position 11, and the first 3 nucleotides of the predicted CDS are positions 10 to 12.
  The alignment of the sequence `ENTOY100A` to the model (`#=GC RF`
  line) below shows the invalid `ATT` start codon aligned to reference positions 11 to 13.

```
                            vvv
ENTOY100A         -AAATCACCGATTGTGATCGCTTTACCATAAATGAGCATTCTACGTGCATCTTGCGGTGCCATACAATGGTAGAAATTGCCATTCACGTACGTAGCATCA
#=GR ENTOY100A PP ****************************************************************************************************
#=GC RF           GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCATCTTGCGGTGCCATACAATGGTAGAAATTGCCATTCACGTACGTAGCATCA
#=GC RFCOLX..     0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
#=GC RFCOL.X.     0000000001111111111222222222233333333334444444444555555555566666666667777777777888888888899999999990
#=GC RFCOL..X     1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
```

---
### <a name="example-stop"></a>Example *mutend**, *cdsstop**, and *unexleng* alerts

#### alert code: *mutendcd*, *mutendns*, *mutendex*, *cdsstopn*, *cdsstopp*, *unexleng*

#### corresponding error messages: *MUTATION_AT_END*, *UNEXPECTED_LENGTH*

#### Example lines from `.alt` file:

```
#      seq                         ftr   ftr          ftr  alert           alert                           seq     seq       mdl     mdl  alert 
#idx   name              model     type  name         idx  code      fail  desc                         coords  length    coords  length  detail
#----  ----------------  --------  ----  -----------  ---  --------  ----  -------------------------  --------  ------  --------  ------  ------
2.1.1  ENTOY50A.mutend1  ENTOY50A  CDS   protein_one    1  mutendcd  yes   MUTATION_AT_END            28..30:+       3  29..31:+       3  expected stop codon could not be identified, predicted CDS stop by homology is invalid [TAC ending at sequence position 30 (model position 31) on + strand is not a valid stop]
2.1.2  ENTOY50A.mutend1  ENTOY50A  CDS   protein_one    1  mutendns  yes   MUTATION_AT_END                   -       -         -       -  expected stop codon could not be identified, no in-frame stop codon exists 3' of predicted start codon
#
3.1.1  ENTOY50A.mutend2  ENTOY50A  CDS   protein_one    1  mutendcd  yes   MUTATION_AT_END            27..29:+       3  28..30:+       3  expected stop codon could not be identified, predicted CDS stop by homology is invalid [ATA ending at sequence position 29 (model position 30) on + strand is not a valid stop]
3.1.2  ENTOY50A.mutend2  ENTOY50A  CDS   protein_one    1  mutendex  yes   MUTATION_AT_END            31..33:+       3  33..35:+       3  expected stop codon could not be identified, first in-frame stop codon exists 3' of predicted stop position [sequence positions 31 to 33 (model positions 33 to 35) on + strand]
3.2.3  ENTOY50A.mutend2  ENTOY50A  CDS   protein_one    1  unexleng  yes   UNEXPECTED_LENGTH          10..29:+      20  11..30:+      20  length of complete coding (CDS or mat_peptide) feature is not a multiple of 3 [20]
#
4.2.1  ENTOY50A.mutend3  ENTOY50A  CDS   protein_one    1  cdsstopn  yes   CDS_HAS_STOP_CODON         22..24:+       3  23..25:+       3  in-frame stop codon exists 5' of stop position predicted by homology to reference [revised to 10..24 (stop shifted 6 nt)]
4.2.2  ENTOY50A.mutend3  ENTOY50A  CDS   protein_one    1  cdsstopn  yes   CDS_HAS_STOP_CODON         22..24:+       3  23..25:+       3  stop codon in protein-based alignment
```

  **Explanation**: The sequence `ENTOY50A.mutend1` has an invalid
  `TAC` stop codon for the CDS with name `protein one` at the
    expected position in the reference model, causing the *mutendcd*
  alert. There are zero valid in-frame stop codons (in the same frame as the predicted start codon) in the 
  remainder of the sequence, so a *mutendns* alert is also reported.

    The sequence `ENTOY50.mutend2` also has an invalid stop codon, the
    final 3 nt in the predicted CDS are `ATA`, which are sequence
    positions 27 to 29 and align to model reference positions 28 to 30
    (see alignment below). There is an in-frame stop codon 3' of this
    though, at sequence positions 31 to 33 which align to model
    reference positions 33 to 35, which result in the *mutendex*
    alert. Finally, the predicted CDS is length 20 (positions 10 to
    29) and aligns to positions 11 to 30, this CDS length is not a
    multiple of 3 so the *unexleng* alert is reported.

    The sequence `ENTOY50.earlystop` has a valid stop codon at the
    expected position (reference positions 29 to 31) but also has an
    earlier in-frame stop codon at reference positions 23 to 25, which
    correspond to sequence positions 22 to 24. This results in
    *cdsstopn* and *cdsstopp* alerts. The *cdsstopn* alert is reported
    when an early stop is identified during the nucleotide-based
    alignment stage, and the *cdsstopp* alert is reported when an
    early stop is detected in the protein-based alignment of the
    protein validation stage. They are often reported together but not
    always, because some cases of predicted earlier stops are
    detectable only in one of the two types of alignment.

    The alignment of the three sequences is below, with the three positions of the 
    stop codon in the reference model marked with `vvv`. The reference positions 
    that align to the early stop codon in the `ENTOY50.mutend3` sequence are 
    marked with `^^^`.

```
                                                       vvv
ENTOY50A.mutend1           -AAATCACCGATGGTGATCGCTTTACCATACATGAGCAT-----------
#=GR ENTOY50A.mutend1 PP   .**************************************...........
ENTOY50A.mutend2           -AAATCACCGATGGTGATCGCTTTACCATA-CTGAGCAT-----------
#=GR ENTOY50A.mutend2 PP   .***************************96.69******...........
ENTOY50A.earlystop         -AAATCACCGATGGTGATCGCTTAACCATAAATGAGCAT-----------
#=GR ENTOY50A.earlystop PP .**************************************...........
#=GC RF                    GAAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.               00000000011111111112222222222333333333344444444445
#=GC RFCOL.X               12345678901234567890123456789012345678901234567890
                                                 ^^^
```                          

---
### <a name="example-indf5"></a>Examples of INDEFINITE_ANNOTATION_START

#### alert codes: *indf5gap*, *indf5lcc*, *indf5pst*, *indf5lcn* (not shown), *indf5plg* (not shown)

#### corresponding error message: *INDEFINITE_ANNOTATION_START*

#### Example lines from `.alt` file:

```
#      seq                    ftr   ftr          ftr  alert           alert                             seq     seq       mdl     mdl  alert 
#idx   name         model     type  name         idx  code      fail  desc                           coords  length    coords  length  detail
#----  -----------  --------  ----  -----------  ---  --------  ----  ---------------------------  --------  ------  --------  ------  ------
2.2.4  ENTOY50A.5A  ENTOY50A  CDS   protein_one    1  indf5gap  yes   INDEFINITE_ANNOTATION_START  10..10:+       1  11..11:+       1  alignment to homology model is a gap at 5' boundary [RF position 11]
#
3.2.1  ENTOY50A.5B  ENTOY50A  CDS   protein_one    1  indf5lcc  no    INDEFINITE_ANNOTATION_START    8..8:+       1  11..11:+       1  alignment to homology model has low confidence at 5' boundary for feature that is or matches a CDS [0.70 < 0.80, RF position 11]
#
4.2.1  ENTOY50A.5C  ENTOY50A  CDS   protein_one    1  indf5pst  yes   INDEFINITE_ANNOTATION_START  10..18:+       9  11..19:+       9  protein-based alignment does not extend close enough to nucleotide-based alignment 5' endpoint [9>5]
```

  **Explanation of `indf5gap` alert**: The sequence `ENTOY50A.5A` has a gap at the 5'
  boundary of the CDS named `protein one` causing the `indf5gap`
  alert. The first *non-gap* nucleotide in the predicted CDS is
  position 10 of the sequence. The 5' boundary gap is position 11 in
  the model.

  **Explanation of `indf5lcc` alert**:
  The alignment of the sequence `ENTOY50A.5B` has low confidence in
  the alignment of the nucleotide at position 8 of the sequence to the
  5' boundary of the `protein one` CDS at model position 11.
   The `#=GR `ENTOY50A.5B PP` line in
  the alignment shows the per position posterior probabilities for
  each aligned nucleotide and position 30 has a value of `7`
  indicating that the posterior probability is between 0.65 and 0.75.

  A similar `indf5lcn` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding
  features cause a sequence to fail or not.

  **Explanation of `indf5pst` alert**:
  The protein validation stage for the `ENTOY50A.5C` sequence resulted in a blastx
  alignment of the predicted nucleotide CDS to the reference protein that did not 
  extend close enough to the 5' boundary of the predicted CDS. The blastx alignment
  extended to sequence position 19, leaving positions 10 to 18 of the predicted CDS, which align to 
  model positions 11 to 19, uncovered by the blastx alignment. In the alignment
  below three nucleotide changes relative to the model in the first three codons
  are marked with `^` at the bottom of the alignment.

  This alert is actually not reported by `v-annotate.pl` for this toy example because the 
  CDS in this example is too short for constructive use with `blastx`. This alert is 
  fabricated here for purposes of illustration. 

  A similar `indf5plg` alert exists for when the blastx alignment extends *longer* than the 
  predicted nucleotide alignment, but no example is shown here. 

The alignment of the three sequences with `indf5*` alerts is
below. The first position of the start codon, model reference position
11 is marked with a `v` at the top. The three positions with
nucleotide changes at the first positions of the first 3 
codons in the CDS are marked at the bottom with `^` characters. These
three mutations cause amino acid changes relative to the model
and cause the blastx alignment to stop prior to the 3' end of the
CDS leading to the *indf5pst* alert.

```
                              v
ENTOY50A.5A         -AAATCACCG-TGGTGATCGCTTTACCATAAATGAGCAT-----------
#=GR ENTOY50A.5A PP .********9.****************************...........
ENTOY50A.5B         -AAATCA--CATGGTGATCGCTTTACCATAAATGAGCAT-----------
#=GR ENTOY50A.5B PP .*****9..67****************************...........
ENTOY50A.5C         -AAATCACCGCTGCTGTTCGCTTTACCATAAATGAGCAT-----------
#=GR ENTOY50A.5C PP .**************************************...........
#=GC RF             GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.        00000000011111111112222222222333333333344444444445
#=GC RFCOL.X        12345678901234567890123456789012345678901234567890
                              ^  ^  ^                         
```                                             

### <a name="example-indf3"></a>Examples of INDEFINITE_ANNOTATION_END

#### alert codes: *indf3gap*, *indf3lcc*, *indf3pst*, *indf3lcn* (not shown), *indf3plg* (not shown)

#### corresponding error message: *INDEFINITE_ANNOTATION_END*

#### Example lines from `.alt` file:

```
#      seq                    ftr   ftr          ftr  alert           alert                             seq     seq       mdl     mdl  alert 
#idx   name         model     type  name         idx  code      fail  desc                           coords  length    coords  length  detail
#----  -----------  --------  ----  -----------  ---  --------  ----  ---------------------------  --------  ------  --------  ------  ------
5.2.4  ENTOY50A.3A  ENTOY50A  CDS   protein_one    1  indf3gap  yes   INDEFINITE_ANNOTATION_END    28..28:+       1  31..31:+       1  alignment to homology model is a gap at 3' boundary [RF position 31]
#
6.2.1  ENTOY50A.3B  ENTOY50A  CDS   protein_one    1  indf3lcc  no    INDEFINITE_ANNOTATION_END    30..30:+       1  31..31:+       1  alignment to homology model has low confidence at 3' boundary for feature that is or matches a CDS [0.70 < 0.80, RF position 31]
#
7.2.1  ENTOY50A.3C  ENTOY50A  CDS   protein_one    1  indf3pst  yes   INDEFINITE_ANNOTATION_END    22..30:+       9  23..31:+       9  protein-based alignment does not extend close enough to nucleotide-based alignment 3' endpoint [9>8]
```

  **Explanation of `indf3gap` alert**: The sequence `ENTOY50A.3A` has a gap at the 3'
  boundary of the CDS named `protein one` causing the `indf3gap`
  alert. The final *non-gap* nucleotide in the predicted CDS is
  position 28 of the sequence. The 3' boundary gap is position 31 in
  the reference model.

  **Explanation of `indf3lcc` alert**: The alignment of the sequence
  `ENTOY50A.3B` has low confidence in the alignment of the nucleotide
  at position 30 of the sequence to the 3' boundary of the `protein
  one` CDS at model position 31. The `#=GR `ENTOY50A.3B PP` line in
  the alignment shows the per position posterior probabilities for
  each aligned nucleotide and position 30 has a value of `7`
  indicating that the posterior probability is between 0.65 and 0.75.

  A similar `indf3lcn` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding
  features cause a sequence to fail or not.

  **Explanation of `indf3pst` alert**:
  The protein validation stage for the `ENTOY50A.3C` sequence resulted in a blastx
  alignment of the predicted nucleotide CDS to the reference protein that did not 
  extend close enough to the 3' boundary of the predicted CDS. The blastx alignment
  extended to sequence position 22, leaving positions 22 to 30 of the predicted CDS, which align to 
  model positions 23 to 31, uncovered by the blastx alignment. In the alignment
  below three nucleotide changes relative to the model in the final three codons
  are marked with `^` at the bottom of the alignment.

  This alert is actually not reported by `v-annotate.pl` for this toy example because the 
  CDS in this example is too short for constructive use with `blastx`. This alert is 
  fabricated here for purposes of illustration. 

  A similar `indf3plg` alert exists for when the blastx alignment extends *longer* than the 
  predicted nucleotide alignment, but no example is shown here. 

The alignment of the three sequences with `indf3*` alerts is
below. The final position of the stop codon, model reference position
31 is marked with a `v` at the top. The three positions with
nucleotide changes at the first positions of the final 3 non-stop
codons in the CDS are marked at the bottom with `^` characters. These
three mutations would cause amino acid changes relative to the model
and cause the blastx alignment to stop prior to the 3' end of the
CDS.

```                                             
                                                  v
ENTOY50A.3A         -AAATCACCGATGGTGATCGCTTTACCAT--ATGAGCAT-----------
#=GR ENTOY50A.3A PP .**************************98..59******...........
ENTOY50A.3B         -AAATCACCGATGGTGATCCCTCTAGCATAGA-GAGCAT-----------
#=GR ENTOY50A.3B PP .****************************976.9*****...........
ENTOY50A.3C         -AAATCACCGATGGTGATCGCTATAGCAGAAATGAGCAT-----------
#=GR ENTOY50A.3C PP .**************************************...........
#=GC SS_cons        :::::::::<<<____>>>:::::::::::::::::::::::::::::::
#=GC RF             GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.        00000000011111111112222222222333333333344444444445
#=GC RFCOL.X        12345678901234567890123456789012345678901234567890
                                       ^  ^  ^
```                                             

### <a name="example-ambg5"></a>Examples of *N_AT_START* 

#### alert codes: *ambgnt5s*, *ambgnt5c*, *ambgnt5f* 

#### corresponding error message: *N_AT_START*, *N_AT_CDS_START*, *N_AT_FEATURE_START*

#### Example lines from `.alt` file:

```
#      seq                     ftr   ftr          ftr  alert           alert                       seq     seq       mdl     mdl  alert 
#idx   name          model     type  name         idx  code      fail  desc                     coords  length    coords  length  detail
#----  ------------  --------  ----  -----------  ---  --------  ----  ---------------------  --------  ------  --------  ------  ------
2.1.1  ENTOY50A.N5L  ENTOY50A  -     -              -  ambgnt5s  no    N_AT_START              1..13:+      13   2..14:+      13  first nucleotide of the sequence is an N [first 13 positions are Ns, first non-N is position 14]
2.2.1  ENTOY50A.N5L  ENTOY50A  CDS   protein_one    1  ambgnt5c  no    N_AT_CDS_START         10..13:+       4  11..14:+       4  first nucleotide of CDS is an N [first 4 positions are Ns, first non-N is position 14]
#
3.2.1  ENTOY50A.N5S  ENTOY50A  CDS   protein_one    1  ambgnt5c  no    N_AT_CDS_START         10..13:+       4  11..14:+       4  first nucleotide of CDS is an N [first 4 positions are Ns, first non-N is position 14]
#
```
  **Explanation of `ambgnt5s` alert**: The sequence `ENTOY50A.N5L` begins with 13 consecutive Ns
  that align to model (reference) positions 2 to 14. 
  Any sequence that begins with 1 or consecutive Ns will cause this alert.
  
  **Explanation of `ambgnt5c` alert**: Both sequences `ENTOY50A.N5L`
  and `ENTOY50A.N5S` have 4 consecutive Ns from sequence positions 10
  to 13 that align to the first four positions of the CDS `protein
  one`, which is model (reference) positions 11 to 14. You can 
  see this in the alignment below.

  A similar `ambgnt5f` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding
  features cause a sequence to fail or not.

The alignment of the two sequences is below, the first four positions of the 
CDS are marked by `v` chacacters:

```                                             
                               vvvv
ENTOY50A.N5L         -NNNNNNNNNNNNNTGATCGCTTTACCATAAATGAGCAT-----------
#=GR ENTOY50A.N5L PP .**************************************...........
ENTOY50A.N5S         -AAATCACCGNNNNTGATCGCTTTACCATAAATGAGCAT-----------
#=GR ENTOY50A.N5S PP .**************************************...........
#=GC RF              GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.         00000000011111111112222222222333333333344444444445
#=GC RFCOL.X         12345678901234567890123456789012345678901234567890
```                                             

---
### <a name="example-ambg3"></a>Examples of *N_AT_END* 

#### alert codes: *ambgnt3s*, *ambgnt3c*, *ambgnt3f* 

#### corresponding error message: *N_AT_END*, *N_AT_CDS_END*, *N_AT_FEATURE_END*

#### Example lines from `.alt` file:

```
#      seq                     ftr   ftr          ftr  alert           alert                       seq     seq       mdl     mdl  alert 
#idx   name          model     type  name         idx  code      fail  desc                     coords  length    coords  length  detail
#----  ------------  --------  ----  -----------  ---  --------  ----  ---------------------  --------  ------  --------  ------  ------
4.1.1  ENTOY50A.N3L  ENTOY50A  -     -              -  ambgnt3s  no    N_AT_END               28..38:+      11  29..39:+      11  final nucleotide of the sequence is an N [final 11 positions are Ns, final non-N is position 27]
4.2.1  ENTOY50A.N3L  ENTOY50A  CDS   protein_one    1  ambgnt3c  no    N_AT_CDS_END           29..30:+       2  30..31:+       2  final nucleotide of CDS is an N [final 2 positions are Ns, final non-N is position 28]
#
5.2.1  ENTOY50A.N3S  ENTOY50A  CDS   protein_one    1  ambgnt3c  no    N_AT_CDS_END           29..30:+       2  30..31:+       2  final nucleotide of CDS is an N [final 2 positions are Ns, final non-N is position 28]
#
```
  **Explanation of `ambgnt3s` alert**: The sequence `ENTOY50A.N3L` ends with 11 consecutive Ns
  from positions 28 to 38, which align to model (reference) positions 29 to 39.
  Any sequence that ends with 1 or consecutive Ns will cause this alert.
  
  **Explanation of `ambgnt3c` alert**: Both sequences `ENTOY50A.N3L`
  and `ENTOY50A.N3S` have 2 consecutive Ns from sequence positions 29
  to 30 that align to the final two positions of the CDS `protein
  one`, which is model (reference) positions 30 to 31. You can 
  see this in the alignment below.

  A similar `ambgnt3f` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding
  features cause a sequence to fail or not.

The alignment of the two sequences is below, the final two positions of the 
CDS are marked by `v` chacacters:


```                                             
                                                  vv
ENTOY50A.N3S         -AAATCACCGATGGTGATCGCTTTACCATNNATGAGCAT-----------
#=GR ENTOY50A.N3S PP .**************************************...........
ENTOY50A.N3L         -AAATCACCGATGGTGATCGCTTTACCATNNNNNNNNNN-----------
#=GR ENTOY50A.N3L PP .**************************************...........
#=GC RF              GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.         00000000011111111112222222222333333333344444444445
#=GC RFCOL.X         12345678901234567890123456789012345678901234567890
```                                             

---
### <a name="example-pep"></a>Examples of *PEPTIDE_ADJACENCY_PROBLEM* and *PEPTIDE_TRANSLATION_PROBLEM*

#### alert codes: *pepadjcy*, *peptrans*

#### corresponding error message: *PEPTIDE_ADJACENCY_PROBLEM*, *PEPTIDE_TRANSLATION_PROBLEM*

#### Example lines from `.alt` file:

```
#      seq                     ftr          ftr              ftr  alert           alert                             seq     seq       mdl     mdl  alert 
#idx   name          model     type         name             idx  code      fail  desc                           coords  length    coords  length  detail
#----  ------------  --------  -----------  ---------------  ---  --------  ----  ---------------------------  --------  ------  --------  ------  ------
2.1.1  ENTOY50A.MP1  ENTOY50A  mat_peptide  protein_one_mp1    2  pepadjcy  yes   PEPTIDE_ADJACENCY_PROBLEM    22..24:+       3  22..23:+       2  predictions of two mat_peptides expected to be adjacent are not adjacent [abs(21 - 25) != 1 (strand:+)]
#
3.1.1  ENTOY50A.MP2  ENTOY50A  CDS          protein_one        1  mutstart  yes   MUTATION_AT_START            10..12:+       3  11..13:+       3  expected start codon could not be identified [GTG starting at sequence position 10 (model position 11) on + strand is not a valid start]
3.2.1  ENTOY50A.MP2  ENTOY50A  mat_peptide  protein_one_mp1    2  peptrans  yes   PEPTIDE_TRANSLATION_PROBLEM         -       -         -       -  mat_peptide may not be translated because its parent CDS has a problem
3.3.1  ENTOY50A.MP2  ENTOY50A  mat_peptide  protein_one_mp2    3  peptrans  yes   PEPTIDE_TRANSLATION_PROBLEM         -       -         -       -  mat_peptide may not be translated because its parent CDS has a problem
```

  **Explanation of `pepadjcy` alert**: The three nucleotides from
  positions 22 to 24 in the sequence `ENTOY50A.MP1` are inserted
  between the end of the `mat_peptide` called `protein one mp1` and
  the second `mat_peptide` called `protein one mp2` (not shown), which
  means that they are not annotated as adjacent as they are expected
  to be. The `protein one mp1` mature peptide ends at model
  (reference) position 22 and the `protein one mp2` mature peptide
  begins at model (reference) position 23. This alert is specific to
  only `mat_peptide` features. 

  **Explanation of `peptrans` alert**: The CDS `protein_one` which is
  the *parent* CDS of the `protein one mp1` and `protein one mp2`
  mature peptides has a *MUTATION_AT_START* error, which means that
  the children mature peptides may not be translated, so a `peptrans`
  alert is reported. This alert is specific to only `mat_peptide`
  features and occurs for all children mature peptides when the parent
  CDS has any fatal alert. The parent/child relationship is defined
  in the *modelinfo* file, an example of which is [here](formats.md#minfo).

The alignment of the two sequences is below, the three inserted `G`
nucleotides in `ENTOY50A.MP1` are marked by `v` chacacters at the top
of the alignment. The mutated and invalid `GTG` start codon in
`ENTOY50A.MP2` is marked by `^` chacacters at the bottom of the
alignment.

```                                             
                                           vvv
ENTOY50A.MP1         -AAATCACCGATGGTGATCGCTGGGTTACCATAAATGAGCAT-----------
#=GR ENTOY50A.MP1 PP .********************97779****************...........
ENTOY50A.MP2         -AAATCACCGGTGGTGATCGCT---TTACCATAAATGAGCAT-----------
#=GR ENTOY50A.MP2 PP .*********************...*****************...........
#=GC RF              GAAATCACCGatGGTGatCGCT...TTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.         0000000001111111111222...2222222333333333344444444445
#=GC RFCOL.X         1234567890123456789012...3456789012345678901234567890
                               ^^^   
```                                             

---
### <a name="example-lowsim5"></a>Examples of 5' low similarity problems

#### alert codes: *lowsim5c*, *lowsim5n* (not shown), *lowsim5s*

#### corresponding error message: *LOW_FEATURE_SIMILARITY_START*, *LOW_SIMILARITY_START*

#### Example lines from `.alt` file:

```
#      seq                     ftr          ftr              ftr  alert           alert                               seq     seq       mdl     mdl  alert 
#idx   name          model     type         name             idx  code      fail  desc                             coords  length    coords  length  detail
#----  ------------  --------  -----------  ---------------  ---  --------  ----  -----------------------------  --------  ------  --------  ------  ------
2.2.1  ENTOY50A.LS1  ENTOY50A  CDS          protein_one        1  lowsim5c  no    LOW_FEATURE_SIMILARITY_START   12..13:+       2  11..12:+       2  region within annotated feature that is or matches a CDS at 5' end of sequence lacks significant similarity [2 nt overlap b/t low similarity region of length 13 (1..13) and annotated feature (12..31), strand: +]
#
3.1.2  ENTOY50A.LS2  ENTOY50A  -            -                  -  lowsim5s  yes   LOW_SIMILARITY_START            1..10:+      10   3..10:+       8  significant similarity not detected at 5' end of the sequence [seq:1..10:+;mdl:3..10:+;low similarity region of length 10 (1..10)]
```

  **Explanation of `lowsim5c` alert**: The two nucleotides from
  positions 12 to 13 in the sequence `ENTOY50A.LS1` are the 5' end
  of the predicted CDS `protein one` but are not similar to the reference
  model. These two nucleotides align to reference model positions 11 and 12.
  The alignment below shows the region of low similarity actually extends from 
  sequence position 1 to position 13, as reported in the `alert detail` field.

  A similar `lowsim5f` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding

  **Explanation of `lowsim5s` alert**: The first ten nucleotides 
  in the sequence `ENTOY50A.LS2` are not similar to the reference model 
  and do not overlap with any predicted features.  These ten nucleotides align to reference model positions 3 to 10,
  and include insertions after reference position 10.
  The alignment below shows the region of low similarity.

  Regions of low similarity are detected in the *coverage determination* 
  stage, as regions that are not covered by local alignment *hits* between the 
  sequence and the model, not based on the global alignment determined in the 
  alignment stage.

The alignment of the two sequences is below, the two nucleotides that
overlap with the beginning of the CDS in `ENTOY50A.LS1` are marked with 
`v` characters at the top of the alignment. The two inserted nucleotides
in sequence `ENTOY50A.LS2` after model position 10 are marked by `^` characters at the bottom 
of the alignment. Note that the nucleotides in the low similarity regions at
the beginning of each sequence do not match the nucleotides in the reference
model (`#=GC RF` line).


```                                             
                                  vv
ENTOY50A.LS1         TTGTAGGTTCG--AC-GTGATCGCTTTACCATAAATGAGCAT-----------
#=GR ENTOY50A.LS1 PP 64156677889..97.79************************...........
ENTOY50A.LS2         ---GTTTAGTGGCATGGTGATCGCTTTACCATAAATGAGCAT-----------
#=GR ENTOY50A.LS2 PP ...******99889****************************...........
#=GC SS_cons         ::.:::::::<..<<____>>>:::::::::::::::::::::::::::::::
#=GC RF              GA.AATCACCG..atGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.         00.00000001..1111111112222222222333333333344444444445
#=GC RFCOL.X         12.34567890..1234567890123456789012345678901234567890
                                ^^
```                                             

---
### <a name="example-lowsim3"></a>Examples of 3' low similarity problems

#### alert codes: *lowsim3c*, *lowsim3n* (not shown), *lowsim3s*

#### corresponding error message: *LOW_FEATURE_SIMILARITY_END*, *LOW_SIMILARITY_END*

#### Example lines from `.alt` file:

```
#      seq                     ftr          ftr              ftr  alert           alert                               seq     seq       mdl     mdl  alert 
#idx   name          model     type         name             idx  code      fail  desc                             coords  length    coords  length  detail
#----  ------------  --------  -----------  ---------------  ---  --------  ----  -----------------------------  --------  ------  --------  ------  ------
4.2.1  ENTOY50A.LS3  ENTOY50A  CDS          protein_one        1  lowsim3c  no    LOW_FEATURE_SIMILARITY_END     26..30:+       5  27..31:+       5  region within annotated feature that is or matches a CDS at 3' end of sequence lacks significant similarity [5 nt overlap b/t low similarity region of length 12 (26..37) and annotated feature (10..30), strand: +]
#
5.1.2  ENTOY50A.LS4  ENTOY50A  -            -                  -  lowsim3s  yes   LOW_SIMILARITY_END             40..53:+      14  41..50:+      10  significant similarity not detected at 3' end of the sequence [seq:40..53:+;mdl:41..50:+;low similarity region of length 14 (40..53)]
```

  **Explanation of `lowsim3c` alert**: The five nucleotides from
  positions 26 to 30 in the sequence `ENTOY50A.LS3` are the 3' end
  of the predicted CDS `protein one` but are not similar to the reference
  model. These two nucleotides align to reference model positions 27 to 31.
  The alignment below shows the region of low similarity actually extends from 
  sequence position 26 to the end of the sequence, as reported in the `alert detail` field.

  A similar `lowsim3f` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding

  **Explanation of `lowsim3s` alert**: The final 14 nucleotides 
  in the sequence `ENTOY50A.LS4` are not similar to the reference model 
  and do not overlap with any predicted features. These 14 nucleotides align to reference model positions 41 to 50,
  and include two inserted nucleotides after reference position 50.
  The alignment below shows the region of low similarity.

  Regions of low similarity are detected in the *coverage determination* 
  stage, as regions that are not covered by local alignment *hits* between the 
  sequence and the model, not based on the global alignment determined in the 
  alignment stage.

The alignment of the two sequences is below, the five nucleotides that
overlap with the end of the CDS in `ENTOY50A.LS3` are marked with 
`v` characters at the top of the alignment. The two inserted nucleotides
in sequence `ENTOY50A.LS4` after model position 50 are marked by `^` characters at the bottom 
of the alignment. Note that the nucleotides in the low similarity regions at
the end of each sequence do not match the nucleotides in the reference
model (`#=GC RF` line).


```                                             
                                               vvvvv              
ENTOY50A.LS3         -AAATCACCGATGGTGATCGCTTTACGTCCCGTCTTAA----------------
#=GR ENTOY50A.LS3 PP .*************************************................
ENTOY50A.LS4         -AAATCACCGATGGTGATCGCTTTACCATAAATGACGATACACGAACTGCACGA
#=GR ENTOY50A.LS4 PP .**********************************8887651156667777566
#=GC SS_cons         :::::::::<<<____>>>::::::::::::::::::::::..:::::::::..
#=GC RF              GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTC..TACGTGCAT..
#=GC RFCOLX.         00000000011111111112222222222333333333344..444444445..
#=GC RFCOL.X         12345678901234567890123456789012345678901..234567890..
                                                                         ^^
```                                             
---
### <a name="example-lowsimi"></a>Examples of internal low similarity problems

#### alert codes: *lowsimic*, *lowsimin* (not shown), *lowsimis (not shown)*

#### corresponding error message: *LOW_FEATURE_SIMILARITY*, *LOW_SIMILARITY*

#### Example lines from `.alt` file:

```
#      seq                     ftr          ftr              ftr  alert           alert                               seq     seq       mdl     mdl  alert 
#idx   name          model     type         name             idx  code      fail  desc                             coords  length    coords  length  detail
#----  ------------  --------  -----------  ---------------  ---  --------  ----  -----------------------------  --------  ------  --------  ------  ------
6.2.3  ENTOY50A.LS5  ENTOY50A  CDS          protein_one        1  lowsimic  no    LOW_FEATURE_SIMILARITY          32..106:+      75  16..22:+       7  region within annotated feature that is or matches a CDS lacks significant similarity [75 nt overlap b/t low similarity region of length 75 (32..106) and annotated feature (10..114), strand: +]
```

  **Explanation of `lowsimic` alert**: The 75 nucleotides from
  positions 32 to 106 in the sequence `ENTOY50A.LS5` occur in the middle
  of the predicted CDS `protein one` but are not similar to the reference
  model. These 75 nucleotides align to reference model positions 16 to 22.
  Most of them are actually inserted after reference position 16.
  The alignment below shows the region of low similarity marked with `v` characters.

  A similar `lowsimin` alert exists for non-coding (non-CDS and
  non-mature peptide) features, but no example is shown here.  Having
  separate alerts for coding and non-coding features gives the user
  control over whether these types of alerts in coding versus non-coding

  A similar `lowsimic` alert exists when the region of low similarity is 
  not within an annotated feature, and also does not include the first or 
  final nucleotide in the sequence, but no example is shown here.

  Regions of low similarity are detected in the *coverage determination* 
  stage, as regions that are not covered by local alignment *hits* between the 
  sequence and the model, not based on the global alignment determined in the 
  alignment stage. That is why the region marked with `v` characters below
  which corresponds to the region of low similartity does not include the 
  full insertion after reference position 16. Enough similarity exists in 
  the beginning of that insertion for a local alignment hit to include that region. 

The alignment of the sequence is below, the 75 nucleotides in the 
low similarity region are marked with 
`v` characters at the top of the alignment. 


```
                                                        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv           
ENTOY50A.LS5         -A-AATCACCG--ATGGTGATCGCTTTACCAAAGCAGTACAGGCACATGACAAAGCAGTACAGGCACATGACAAAGCAGTACAGGCACATGACAAAGCAGTACAGGCACAT-GACATAAATGAGCATTC--TACGTGCAT--
#=GR ENTOY50A.LS5 PP .*.********..**987555566666666666666666666666666677777777777777777777788888888888888888888888886665555544455544.4589*************..*********..
#=GC SS_cons         ::.:::::::<..<<____.....................................................................................>>>::::::::::::::::::::::..:::::::::..
#=GC RF              GA.AATCACCG..atGGTG.....................................................................................atCGCTTTACCATAAATGAGCATTC..TACGTGCAT..
#=GC RFCOLX..        00.00000000..000000.....................................................................................0000000000000000000000000..000000000..
#=GC RFCOL.X.        00.00000001..111111.....................................................................................1112222222222333333333344..444444445..
#=GC RFCOL..X        12.34567890..123456.....................................................................................7890123456789012345678901..234567890..
                       
```                                             
---
### <a name="example-deletin"></a>Examples of deleted feature problems

#### alert codes: *deletins*, *deletinf* (not shown)

#### corresponding error message: *DELETION_OF_FEATURE*, *DELETION_OF_FEATURE_SECTION*

#### Example lines from `.alt` file:

```
#      seq                     ftr   ftr   ftr  alert           alert                     seq     seq       mdl     mdl  alert 
#idx   name          model     type  name  idx  code      fail  desc                   coords  length    coords  length  detail
#----  ------------  --------  ----  ----  ---  --------  ----  -------------------  --------  ------  --------  ------  ------
2.1.1  ENTOY50A.DF1  ENTOY50A  -     -       -  deletins  yes   DELETION_OF_FEATURE         -       -  23..28:+       6  internal deletion of a complete feature [mat_peptide feature number 2: protein one mp2]
```

  **Explanation of `deletins` alert**: The mature peptide feature named `protein one mp2` 
  is deleted in the sequence `ENTOY50A.DF1`. This is inferred based on the alignment of the sequence
  to the model in region from model positions 23 to 28 which correspond to that mature peptide.
  
  A similar `deletinf` alert exists for *multi-segment* features for which
  one or more but not all segments are deleted. The associated value for this
  alert in the `alert desc` field is `DELETION_OF_FEATURE_SECTION`. 
  No example of this alert is shown here.
  
The alignment of the sequence is below, note that the model positions 23 to 28, 
marked by `v` characters at the top of the alignment, are all gaps in the sequence

```
                                           vvvvvv
ENTOY50A.DF1         -AAATCACCGATGGTGATCGC--------AAATGAGCATTCTACGTGCAT
#=GR ENTOY50A.DF1 PP .*******************8........89*******************
#=GC RF              GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.         00000000011111111112222222222333333333344444444445
#=GC RFCOL.X         12345678901234567890123456789012345678901234567890
```                                             

---
### <a name="example-dupregin"></a>Example of duplicate region

#### alert codes: *dupregin*

#### corresponding error message: *DUPLICATE_REGION*

#### Example lines from `.alt` file:

```
#      seq                     ftr   ftr   ftr  alert           alert                           seq     seq              mdl     mdl  alert 
#idx   name          model     type  name  idx  code      fail  desc                         coords  length           coords  length  detail
#----  ------------  --------  ----  ----  ---  --------  ----  -----------------  ----------------  ------  ---------------  ------  ------
2.1.1  ENTOY50A.DR1  ENTOY50A  -     -       -  dupregin  yes   DUPLICATE_REGIONS  1..49:+,50..86:+      86  2..50:+,7..43:+      86  similarity to a model region occurs more than once [7-43 (len 37 >= 20) hits 1 (37.7 bits) and 2 (25.5 bits) (seq:1..49,50..86 mdl:2..50,7..43)]
```

  **Explanation of `dupregin` alert**: There are two hits in the coverage determination stage, which 
  overlap in model coordinates by more than 20 positions. Specifically, hit 1 is from positions 1 to 49 in the 
  sequence and 2 to 50 in the model and hit 2 is from positions 50 to 86 in the sequence and positions 7 
  to 43 in the model. So the two hits overlap by 37 model positions from 7 to 43. The first hit has a bit
  score of 37.7 bits and the second hit of 25.5 bits. 

The alignment of the sequence is below, and sequence positions 50 to 86 are marked 
marked by `v` characters at the top of the alignment, which match identically
to positions 7 to 43 in the reference model, marked by '^' at the bottom of
the alignment. 

```
                                                                       |positions 50 to 86 in the sequence-|
                                                                       vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
ENTOY50A.DR1         -AAATCACCGATGGTGATCGCTTTACCATAAATGAGCATTCTACGTGCATACCGATGGTGATCGCTTTACCATAAATGAGCATTCTA
#=GR ENTOY50A.DR1 PP .**************************************************************************************
#=GC RF              GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT.....................................
#=GC RFCOLX.         00000000011111111112222222222333333333344444444445.....................................
#=GC RFCOL.X         12345678901234567890123456789012345678901234567890.....................................
                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                           |--positions 7 to 43 in the model---|
```                                             
---
### <a name="example-discontn"></a>Example of discontinuous similarity

#### alert codes: *discontn*

#### corresponding error message: *DISCONTINUOUS SIMILARITY*

#### Example lines from `.alt` file:

```
#      seq                     ftr   ftr   ftr  alert           alert                                  seq     seq               mdl     mdl  alert 
#idx   name          model     type  name  idx  code      fail  desc                                coords  length            coords  length  detail
#----  ------------  --------  ----  ----  ---  --------  ----  ------------------------  ----------------  ------  ----------------  ------  ------
2.1.1  ENTOY50A.DC1  ENTOY50A  -     -       -  discontn  yes   DISCONTINUOUS_SIMILARITY  1..26:+,27..49:+      49  25..50:+,2..24:+      49  not all hits are in the same order in the sequence and the homology model [seq order: 1,2, mdl order: 2,1]
```

  **Explanation of `discontn` alert**: For the `ENTOY5A.DC1` sequence,
  there are two hits in the coverage determination stage and they are
  not in the same order in sequence and model
  coordinates. Specifically, hit 1 is from positions 1 to 26 in the
  sequence and 25 to 50 in the model and hit 2 is from positions 27 to
  49 in the sequence and positions 2 to 24 in the model. So hit 1
  comes before hit 2 in the sequence, but hit 2 comes before hit 1 in
  the model.

The alignment of the sequence is below, and sequence positions 27 to 49 are marked 
marked by `v` characters at the top of the alignment, which match identically
to positions 2 to 24 in the reference model, marked by '^' at the bottom of
the alignment. 

```
                                                                       |-27 to 49 in the seq-|
                                                                       vvvvvvvvvvvvvvvvvvvvvvv
ENTOY50A.DC1         ------------------------ACCATAAATGAGCATTCTACGTGCATAAATCACCGATGGTGATCGCTTT
#=GR ENTOY50A.DC1 PP ........................*************************************************
#=GC RF              GAAATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT.......................
#=GC RFCOLX.         00000000011111111112222222222333333333344444444445.......................
#=GC RFCOL.X         12345678901234567890123456789012345678901234567890.......................
                      ^^^^^^^^^^^^^^^^^^^^^^^
                      |--2-24 in the model--|
```                                             
---
### <a name="example-indfstrn"></a>Example of indefinite strand

#### alert codes: *indfstrn*

#### corresponding error message: *INDEFINITE_STRAND*

#### Example lines from `.alt` file:
```
#      seq                     ftr   ftr   ftr  alert           alert                              seq     seq       mdl     mdl  alert 
#idx   name          model     type  name  idx  code      fail  desc                            coords  length    coords  length  detail
#----  ------------  --------  ----  ----  ---  --------  ----  ----------------------------  --------  ------  --------  ------  ------
2.1.2  ENTOY50A.ID1  ENTOY50A  -     -       -  indfstrn  yes   INDEFINITE_STRAND             49..25:-      25  26..50:+      25  significant similarity detected on both strands [best hit is on + strand, but hit on - strand from 49 to 25 has score 12.8 > 12.0]
```

  **Explanation of `indfstrn` alert**: For the `ENTOY50A.ID1` sequence there are two hits in the coverage determination stage and they 
  are not to the same strand. The top-scoring hit is to the `+` strand, but the second hit with a score of 12.8 bits is on the `-` strand,
  that hit occurs from positions to 49 to 25 on the negative strand to positions 26 to 50 in the reference model. The second hit therefore
  conflicts with the top-scoring hit and so it is the one that is described in the `seq coords` and `mdl coords` fields. Normally, a 
  *indfstrn* alert is only reported if a hit on the opposite strand from the top hit exceeds 25 bits, but to construct this example, the
  minimum of 25 bits was lowered to 12 bits with the `--indefstr 12` flag to `v-annotate.pl`.

The alignment of the sequence is below, and sequence positions 25 to
49 are marked marked by `v` characters at the top of the alignment,
*after reverse complementing* these match identically to positions 2
to 27 in the reference model, marked by '^' at the bottom of the
alignment.

```
                                              |  25 to 49 in the seq  |
                                              vvvvvvvvvvvvvvvvvvvvvvvvv
ENTOY50A.ID1         -AAATCACCGATGGTGATCGCTTTAATGCACGTAGAATGCTCATTTATGG-----
#=GR ENTOY50A.ID1 PP .***********************8654433333499*************.....
#=GC RF              GAAATCACCGatGGTGatCGCTTTACCATA.....AATGAGCATTCTACGTGCAT
#=GC RFCOLX.         000000000111111111122222222223.....33333333344444444445
#=GC RFCOL.X         123456789012345678901234567890.....12345678901234567890
                                              ^^^^^     ^^^^^^^^^^^^^^^^^^^^
                                              |   27 to 50 in the model    |
```                                             
---
### <a name="example-lowcovrg"></a>Example of low coverage

#### alert codes: *lowcovrg*

#### corresponding error message: *LOW_COVERAGE*

#### Example lines from `.alt` file:
```
#      seq                     ftr   ftr   ftr  alert           alert                     seq     seq     mdl     mdl  alert 
#idx   name          model     type  name  idx  code      fail  desc                   coords  length  coords  length  detail
#----  ------------  --------  ----  ----  ---  --------  ----  --------------------  -------  ------  ------  ------  ------
2.1.1  ENTOY50A.LC1  ENTOY50A  -     -       -  lowcovrg  yes   LOW_COVERAGE          1..25:+      25       -       -  low sequence fraction with significant similarity to homology model [0.545 < 0.900]
```

  **Explanation of `lowcovrg` alert**: For the `ENTOY50A.LC1` sequence the nucleotides from 1 to 25 are not part of any hit to the model
  in the coverage determation stage. These 25 nucleotides make up 45.5\% the total length of the sequence, meaning that only 54.5\% of
  the sequence *is covered* by hits to the model. 

The alignment of the sequence is below. Note that positions 1..25 do
not match well to the model. These are marked by `v` characters at the top of the alignment.

```
                     |1 to 25 in the sequence|
                     vvvvvvvvvvvvvvvvvvvvvvvvv
ENTOY50A.LC1         GATACATAACCATAAAATCCGGGACAAATCACCGATGGTGATCGCTTTACCATAA-------------------
#=GR ENTOY50A.LC1 PP *6444444444444444444444444789**************************...................
#=GC RF              GA........................AATCACCGatGGTGatCGCTTTACCATAAATGAGCATTCTACGTGCAT
#=GC RFCOLX.         00........................000000011111111112222222222333333333344444444445
#=GC RFCOL.X         12........................345678901234567890123456789012345678901234567890
```                                             

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
#### Questions, comments or feature requests? Send a
mail to eric.nawrocki@nih.gov.


