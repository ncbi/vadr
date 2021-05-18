# <a name="top"></a> Examples and explanations of `v-annotate.pl` detailed alert and error messages

* [`v-annotate.pl` example usage](#exampleusage)
  * [example annotation of norovirus sequences](#examplebasic)
  * [example of using `--alt_pass` to change alerts from fatal to non-fatal](#examplealtpass)
  * [example of using `-p` to run in parallel mode](#exampleparallel)
* [`v-annotate.pl` command-line options](#options)
  * [basic options](#options-basic)
  * [options for specifying expected sequence classification](#options-classification)
  * [options for controlling which alerts are fatal](#options-fatal)
  * [options related to model files](#options-modelfiles)
  * [options for controlling output feature table](#options-featuretable)
  * [options for controlling alert thresholds](#options-alerts)
  * [options for controlling the alignment stage](#options-align)
  * [options for controlling the blastx protein validation stage](#options-blastx)
  * [options for using hmmer instead of blastx for protein validation](#options-hmmer)
  * [options related to blastn-based seeded alignment acceleration strategy](#options-seed)
  * [options related to pre-processing to replace Ns with expected nucleotides](#options-replace)
  * [options related to splitting input fasta file and multithreading](#options-split)
  * [options related to parallelization on a compute farm/cluster](#options-parallel)
  * [options related to both splitting input and parallelization on a compute farm/cluster](#options-split-and-parallel)
  * [options for skipping stages](#options-skip)
  * [options for additional output files](#options-output)
  * [additional expert options](#options-expert)
* [Basic Information on `v-annotate.pl` alerts](#alerts)
* [Additional information on `v-annotate.pl` alerts](#alerts2)
* [Expendable features: allowing sequences to pass despite fatal alerts for specific features](#mnf)
* [Limiting memory usage and multi-threading](#memory)

---

## Output files with detailed alert and error messages

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

## Example `.alt` output for different alert types:

### Frameshift alerts

#### alert codes: *fsthicf5*, *fsthicf3*, *fsthicfi*, *fstlocf5*, *fstlocf3*, *fstlocfi*, *fstukcf5*, *fstukcf3*, *fstukcfi* 

#### corresponding error messages: *POSSIBLE_FRAMESHIFT_HIGH_CONF*,  *POSSIBLE_FRAMESHIFT_LOW_CONF*, *POSSIBLE_FRAMESHIFT*

#### Example line from `.alt` file:

```
#       seq                        ftr          ftr              ftr  alert           alert                               seq       mdl  alert 
#idx    name            model      type         name             idx  code      fail  desc                             coords    coords  detail
#-----  --------------  ---------  -----------  ---------------  ---  --------  ----  -----------------------------  --------  --------  ------
7.1.2   ENTOY100A-fs6   ENTOY100A  CDS          protein_one        2  fsthicfi  yes   POSSIBLE_FRAMESHIFT_HIGH_CONF  14..25:+  14..22:+  high confidence possible frameshift in CDS (internal) [nucleotide alignment of internal sequence positions 14..25 (12 nt, avgpp: 0.890) to model positions 14..22 (9 nt) on + strand are frame 3 (dominant frame is 1); inserts:S:14..18(5),M:13; deletes:S:25,M:21..22(2);]
```


  **Explanation**: a possible frameshift exists in the CDS named
  `protein_one` in the sequence named `ENTOY100A-fs6` which matches
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
  option is used with `v-annotate.pl` posterior probability values are
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

#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.


