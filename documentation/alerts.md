# <a name="top"></a> Explanations and examples of `v-annotate.pl` detailed alert and error messages

* [Output fiels with detailed alert and error messages](#files)
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
| *fsthicf5*, *fsthicf3*, *fsthicfi*, *fstlocf5*, *fstlocf3*, *fstlocfi*, *fstukcf5*, *fstukcf3*, *fstukcfi* | *POSSIBLE_FRAMESHIFT_HIGH_CONF*,  *POSSIBLE_FRAMESHIFT_LOW_CONF*, *POSSIBLE_FRAMESHIFT* | sequence positions of the frameshifted region | none | model (reference) positions of the frameshifted region, may some nucleotides may be inserted **before or after** these positions | none | [frameshift alert example](#example-frameshift) | 
| *insertnn*, *insertnp* | *INSERTION_OF_NT* | sequence positions of a 'large' insertion with respect to the model |  none | model (reference) position after which insertion occurs | always length 1 | [insert alert example](#example-insert) | 
| *deletinn*, *deletinp* | *DELETION_OF_NT*  | sequence position just prior to (5' of) deletion with respect to the model | always length 1 | model (reference) positions that are deleted in sequence | none | [delete alert example](#example-delete) | 

## Example `.alt` output for different alert types:

###<a name="example-frameshift"></a>Example frameshift alert

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

###<a name="example-insert"></a>Example insert alerts

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
###<a name="example-delete"></a>Example delete alerts

#### alert codes: *insertnp*, *deletinp*

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


#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.


