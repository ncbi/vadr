# <a name="top"></a> Explanations and examples of `v-annotate.pl` detailed alert and error messages

* [Output fields with detailed alert and error messages](#files)
* [Meaning of sequence and model coordinates in `.alt` files](#coords)
* [Example `.alt` output for different alert types](#examples)

---

##<aname="files"></a> Output files with detailed alert and error messages

`v-annotate.pl` outputs two types of files with detailed alert/error
messages: 

  * `.alt` files
  * `.alt.list` files

Documentation on the format of `.alt` files can be found
[here](formats.md#alt) and for `.alt.list` files can be found
[here](formats.md#altlist).

This page includes example lines for all the different alerts from
`.alt` files. Much of the same information, in particular the `seq
coords`, `mdl coords`, and `alert detail` is also present in the
`.alt.list` files. 

## Explanation of sequence and model coordinate fields in `.alt` files

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
| - | - | - | - | - | - | - | 
| *peptrans* | **various** | will be blank (`-`) | N/A | will be blank (`-`) | N/A | - | 

## Example `.alt` output for different alert types:

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
3.2.1  ENTOY50A.5B  ENTOY50A  CDS   protein_one    1  indf5lcc  no    INDEFINITE_ANNOTATION_START    8..8:+       1  11..11:+       1  alignment to homology model has low confidence at 5' boundary for feature that is or matches a CDS [0.70 < 0.80, RF position 11]
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
6.2.1  ENTOY50A.3B  ENTOY50A  CDS   protein_one    1  indf3lcc  no    INDEFINITE_ANNOTATION_END    30..30:+       1  31..31:+       1  alignment to homology model has low confidence at 3' boundary for feature that is or matches a CDS [0.70 < 0.80, RF position 31]
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

---


#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.


