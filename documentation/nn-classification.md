# `v-annotate.pl` nearest-neighbor based classification mode

* [Overview](#overview)
* [Requirements](#requirements)
* [How It Works](#how-it-works)
* [Model Info File Setup](#minfo-setup)
* [Stockholm Alignment File Format](#stk-format)
* [Classification Region Specification](#region-spec)
* [Output Files](#output)
* [Related Alerts](#alerts)
* [Tutorial Example](#tutorial)

---

## Overview<a name="overview"></a>

VADR's nearest-neighbor based classification mode allows sequences to be classified based on their similarity to reference sequences in a training alignment, rather than solely on the best-scoring model match, which is the default mode. This is particularly useful for viruses with well-defined subtypes or serotypes where classification should be based on genetic similarity within specific genomic regions.

In this mode, after a sequence is aligned to its best-matching model, VADR compares it against all sequences in a reference alignment and assigns the group and subgroup of the sequence with the highest percent identity (the "nearest neighbor"). This classification can optionally be restricted to a specific model region, such as VP1 for enteroviruses.

## Requirements<a name="requirements"></a>

To enable nearest-neighbor based classification, you need:

1. **A Stockholm format seed alignment file** (`.stk`) - This must be the alignment used to build the covariance model (CM). The alignment must have the same number of reference (RF) positions as the model. Use `cmbuild -O` to output the seed alignment when building your CM.
2. **Group and subgroup annotations** in the seed alignment file using `#=GS` markup
3. **A modified model info file** (`.minfo`) that references the seed alignment file using special `:FILE:` syntax
4. **(Optional) Region boundaries** specified in the alignment file for region-specific classification

## How It Works<a name="how-it-works"></a>

The classification process follows these steps:

1. **Initial classification**: The sequence is classified to its best-matching model using standard VADR procedures
2. **Alignment**: The sequence is aligned to the model
3. **Nearest-neighbor search**: The aligned sequence is compared to all sequences in the reference alignment
4. **Percent identity calculation**: For each reference sequence, the fractional identity is calculated as the number of identical nucleotides divided by the number of **nongap reference (RF) positions**. **Crucially, gap RF positions (inserts) are completely ignored** - only RF (model) positions are used in the calculation.
5. **Assignment**: The group and subgroup of the reference sequence with the highest percent identity becomes the classification for the input sequence
6. **Second nearest-neighbor**: The reference sequence with the highest identity from a different subgroup is identified to assess classification confidence

If a classification region is specified (e.g., VP1 positions 2467-3393), classification is based only on that region. If the input sequence doesn't sufficiently cover the specified region (default minimum: 40nt), the full sequence is used instead and an alert is reported.

## Model Info File Setup<a name="minfo-setup"></a>

To enable nearest-neighbor classification, modify your model info (`.minfo`) file to use the `:FILE:` syntax for the `group` and `subgroup` fields:

```
MODEL evB group:":FILE:evB.stk" subgroup:":FILE:evB.stk" cmfile:"evB.cm" length:"7562" blastdb:"evB.vadr.protein.fa"
```

**Key points:**
- Both `group` and `subgroup` must reference the same alignment file
- The syntax is `:FILE:<filename>` where `<filename>` is the Stockholm alignment file
- The alignment file should be in the same directory as the `.minfo` file or in `$VADRMODELDIR`
- You can disable this feature with `--ignore_nnclass`

**Example `.minfo` file:**
```
MODEL evB group:":FILE:evB.stk" subgroup:":FILE:evB.stk" cmfile:"evB.cm" length:"7562" blastdb:"evB.vadr.protein.fa"
FEATURE evB type:"gene" coords:"742..7398:+" gene:"POLY" product:"polyprotein"
FEATURE evB type:"CDS" coords:"742..7398:+" gene:"POLY" product:"polyprotein"
FEATURE evB type:"mat_peptide" coords:"742..1518:+" gene:"POLY" product:"VP4"
FEATURE evB type:"mat_peptide" coords:"1519..2316:+" gene:"POLY" product:"VP2"
FEATURE evB type:"mat_peptide" coords:"2317..3036:+" gene:"POLY" product:"VP3"
FEATURE evB type:"mat_peptide" coords:"3037..3945:+" gene:"POLY" product:"VP1"
```

## Stockholm Alignment File Format<a name="stk-format"></a>

The seed alignment must be in Stockholm format with specific annotations. This should be the same alignment file used to build the covariance model.

**Important:** When building your CM with `cmbuild`, use the `--hand` option to preserve the RF annotation:

```bash
cmbuild --noss --hand toy-nn.cm toy-nn.stk
```

The `--hand` option ensures that the RF annotation (which columns are reference vs. insert) from your input alignment is preserved exactly. Without `--hand`, `cmbuild` may modify the RF definition based on conservation patterns.

### Required per-sequence annotations:

Group and subgroup must be specified for each sequence using `#=GS` markup:

```
#=GS <seqname> GP <group>
#=GS <seqname> SG <subgroup>
```

**Example:**
```
# STOCKHOLM 1.0

#=GS AB426608.1 GP EVB
#=GS AB426608.1 SG E30
#=GS HM777023.1 GP EVB
#=GS HM777023.1 SG E18
#=GS KF311743.1 GP EVB
#=GS KF311743.1 SG CV-B5
#=GS MK012537.1 GP EVB
#=GS MK012537.1 SG CV-B3

AB426608.1    TTAAAACAGCCTGTGGGTTGT...
HM777023.1    TTAAAACAGCCTGTGGGTTCT...
KF311743.1    TTAAAACAGCCTGTGGGTTCT...
MK012537.1    TTAAAACAGCCTGTGGGTTGT...
#=GC RF       TTAAAACAGCCTGTGGGTTgT...
//
```

**Key points:**
- `GP` annotations define the group (e.g., "EVB" for Enterovirus B)
- `SG` annotations define the subgroup (e.g., "E30", "CV-B5" for specific serotypes)
- The `#=GC RF` line defines reference (model) positions used for percent identity calculations:
  - Any **non-gap character** (letter, uppercase or lowercase) indicates a nongap RF position used in calculations
  - Any **gap character** (`.` or `-`) indicates an insert column (gap RF position) ignored in calculations
  - Uppercase vs. lowercase letters indicate conservation level (uppercase = more conserved)
- All sequences must have both GP and SG annotations

## Classification Region Specification<a name="region-spec"></a>

To restrict classification to a specific model region (recommended for viruses with hypervariable regions):

### Optional alignment file annotations:

```
#=GF VADR-classification-rf-start-pos <start>
#=GF VADR-classification-rf-stop-pos  <stop>
```

**Example (VP1 region for enterovirus B):**
```
# STOCKHOLM 1.0

#=GF VADR-classification-rf-start-pos 2467
#=GF VADR-classification-rf-stop-pos  3393

#=GS AB426608.1 GP EVB
#=GS AB426608.1 SG E30
...
```

**Key points:**
- Positions are 1-based model reference (RF) positions
- If specified, only these positions are used for calculating percent identity
- Sequences with fewer than 40nt (default) in this region will use the full sequence instead
- You can change the minimum length with `--nnregionlen <n>`
- You can ignore region specifications with `--ignore_nnregion`

## Output Files<a name="output"></a>

### `.scn` output file

When nearest-neighbor classification is enabled, VADR creates a `.scn` (sequence classification neighbor) output file with detailed classification information for each sequence.

**Example `.scn` output:**
```
#seq  seq         seq                           sub     fract                     sub      fract                  fid  nnregion_seqspan    nnregion  nnregion  seq   
#idx  name        len   p/f   ant  model  grp1  grp1      id1  seq1         grp2  grp2     id2     seq2          diff        mdl_coords  mdl_coords     covrg  alerts
#---  ----------  ----  ----  ---  -----  ----  -----  ------  -----------  ----  -------  ------  ----------  ------  ----------------  ----------  --------  ------
1     KX171337.1  7421  PASS  yes  evB    EVB   EV-B106 1.0000  KX171337.1  EVB   EV-B77   0.7168  AJ493062.2  0.2832  2467..3393:+  2467..3393:+      1.0000  -
2     JF416934.1   976  PASS  yes  evB    EVB   EV-B110 1.0000  JF416934.1  EVB   EV-B112  0.7440  KJ418244.1  0.2560  2467..3393:+  2467..3393:+      1.0000  -
```

**Field descriptions:**

| Field | Description |
|-------|-------------|
| `seq idx` | Index of sequence in input file |
| `seq name` | Sequence name |
| `seq len` | Sequence length |
| `p/f` | PASS or FAIL |
| `ant` | yes if annotated, no if not |
| `model` | Best-matching model name |
| `grp1` | Group of nearest-neighbor sequence |
| `sub grp1` | Subgroup of nearest-neighbor sequence |
| `fract id1` | Fractional identity to nearest neighbor (0-1) |
| `seq1` | Name of nearest-neighbor sequence |
| `grp2` | Group of second nearest-neighbor (different subgroup) |
| `sub grp2` | Subgroup of second nearest-neighbor |
| `fract id2` | Fractional identity to second nearest neighbor |
| `seq2` | Name of second nearest-neighbor sequence |
| `fid diff` | Difference: `fract id1 - fract id2` |
| `nnregion_seqspan mdl_coords` | Model positions sequence actually spans |
| `nnregion mdl_coords` | Model region used for classification |
| `nnregion covrg` | Coverage: sequence span / classification region |
| `seq alerts` | Per-sequence alerts for this sequence |

For complete `.scn` format documentation, see [formats.md](formats.md#scn).

### Modified `.mdl` and `.sqc` files

The `.mdl` and `.sqc` files show group/subgroup assignments based on nearest-neighbor classification when enabled.

## Related Alerts<a name="alerts"></a>

VADR reports several alerts specific to nearest-neighbor classification:

| Alert code | Long name | Description | Default fatal? |
|------------|-----------|-------------|----------------|
| `nnindfcl` | INDEFINITE_CLASSIFICATION_NN | Low difference between fractional identity of sequence and its nearest neighbor vs. 2nd nearest neighbor | never |
| `nnloidcl` | LOW_ID_CLASSIFICATION_NN | Low fractional identity of sequence and its nearest neighbor | never |
| `nnalrgcl` | ALT_REGION_CLASSIFICATION_NN | Alternative alignment region used to find nearest neighbor because sequence doesn't include specified region | never |
| `nnptrgcl` | PARTIAL_REGION_CLASSIFICATION_NN | Only part of the specified alignment region used to find nearest neighbor | never |

### Adjusting alert thresholds:

You can modify alert thresholds with command-line options:

```bash
v-annotate.pl --nn_indefclass 0.10    # fractional difference threshold (default: 0.05)
v-annotate.pl --nn_lowidclass 0.70    # fractional identity threshold (default: 0.75)
v-annotate.pl --nn_partclass 0.40     # fractional region coverage threshold (default: 0.50)
```

**Example:** To require 80% identity to nearest neighbor and 10% difference to second nearest:
```bash
v-annotate.pl --nn_lowidclass 0.80 --nn_indefclass 0.10 myseqs.fa output-dir
```

## Tutorial Example<a name="tutorial"></a>

This tutorial demonstrates nearest-neighbor classification using a minimal toy example with 3 enterovirus-like sequences.

### Step 1: Create the reference alignment

Create a file `toy-nn.stk` with three reference sequences from different serotypes (this file is in `documentation/annotate-files/toy-nn.stk`):

```stockholm
# STOCKHOLM 1.0

#=GF VADR-classification-rf-start-pos 7
#=GF VADR-classification-rf-stop-pos  18

#=GS refseq_E30 GP EVB
#=GS refseq_E30 SG E30
#=GS refseq_E18 GP EVB
#=GS refseq_E18 SG E18
#=GS refseq_CVB5 GP EVB
#=GS refseq_CVB5 SG CV-B5

refseq_E30    AAAAACGGGGCCCCCTTTTTAAAA..GGGGGGGG
refseq_E18    AAAAACGGGGTTTTTTTTTTAAAA..GGGGGGGG
refseq_CVB5   AAAAACGGGGAAAAATTTTTAAAAaaGGGGGGGG
#=GC RF       AAAAACGGGGCCCCCTTTTAAAAA..GGGGGGGG
//
```

**Explanation:**
- Three reference sequences representing E30, E18, and CV-B5 serotypes
- Classification region: RF positions 7-18 (the middle variable region)
- `#=GC RF` line shows reference positions:
  - **Any character (uppercase or lowercase letter)** = nongap RF (model) position that **is used** for percent identity calculation
  - **Uppercase letters** = more highly conserved nongap RF positions
  - **Lowercase letters** = less well conserved nongap RF positions
  - **Gaps (periods or dashes)** = insert columns (gap RF positions) that **are ignored** for percent identity calculation
- Insert columns appear at alignment positions 25-26 (`..`, `..`, `aa` in sequences, `..` in RF)
- All three sequences share identical nongap RF nucleotides in flanking regions
- They differ in the classification region (nongap RF positions 7-18)

### Step 2: Build the covariance model

Build the covariance model (CM) from the seed alignment using `cmbuild` (this CM file is in documentation/annotate-files/toy-nn.cm):

```bash
cmbuild --noss --hand toy-nn.cm toy-nn.stk
```

**Explanation:**
- `--noss`: Do not use secondary structure information (our toy example has none)
- `--hand`: **Critical option** - Preserve the RF annotation from the input alignment. Without this option, `cmbuild` may modify which columns are reference vs. insert columns. This option ensures the RF annotation you defined is maintained exactly in the model.
- `toy-nn.cm`: Output CM file
- `toy-nn.stk`: Input seed alignment (this same file will be referenced in the model info file)

**Important**: The input alignment file `toy-nn.stk` already contains all required annotations (`#=GS` for group/subgroup, `#=GF` for classification region, and `#=GC RF` for reference positions). This same file will be used by VADR for nearest-neighbor classification.

### Step 3: Create the model info file

Create `toy-nn.minfo` (this file is in documentation/annotate-files/toy-nn.minfo):

```
MODEL toy-nn group:":FILE:toy-nn.stk" subgroup:":FILE:toy-nn.stk" cmfile:"toy-nn.cm" length:"32"
FEATURE toy-nn type:"gene" coords:"1..32:+" gene:"TEST"
```

**Explanation:**
- The `:FILE:toy-nn.stk` syntax tells VADR to use nearest-neighbor classification
- Both `group` and `subgroup` reference the same seed alignment file (the one with `#=GS` and `#=GF` annotations)
- The alignment file must have the same RF structure as the CM built from it

### Step 4: Create test sequences

Create `test-nn.fa` with sequences to classify (this file is in documentation/annotate-files/test-nn.fa):

```fasta
>seq1_should_be_E30
AAAAAACGGGGCCCCCTTTTTAAAAGGGGGGGG
>seq2_should_be_E18
AAAAAACGGGGTTTTTTTTTTAAAAGGGGGGGG
>seq3_closest_to_CVB5
AAAAAACGGGGAAAAGTTTTTAAAAGGGGGGGG
>seq4_different_inserts
AAAAAACGGGGCCCCCTTTTTAAAATTGGGGGGGG
>seq5_partial_region
AAAAAACGGGGNNNNNNNNNNAAAAGGGGGGGG
```

**Explanation:**
- `seq1`: Identical to refseq_E30 in all RF positions → should classify as E30
- `seq2`: Identical to refseq_E18 in all RF positions → should classify as E18
- `seq3`: One mismatch from refseq_CVB5 in classification region (G vs A at RF position 12) → should classify as CV-B5
- `seq4`: Identical to refseq_E30 in RF positions but has different inserts (TT at insert positions 25-26) - inserts are ignored so should classify as E30 with 100% identity
- `seq5`: Has Ns in classification region → may trigger `nnptrgcl` or `nnalrgcl` alert

### Step 5: Run `v-annotate.pl`

```bash
v-annotate.pl -f --out_stk --mdir $VADRSCRIPTSDIR/documentation/annotate-files --mkey toy-nn --nnregionlen 12 test-nn.fa va-nn
```

### Step 6: Examine the alignment and .scn output

Look at `va-nn/va-nn.vadr.scn`:
```
# STOCKHOLM 1.0
#=GF AU Infernal 1.1.5

seq1_should_be_E30             AAAAACGGG..GCCCCCTTTTTAAAA..GGGGGGGG
#=GR seq1_should_be_E30     PP *********..***************..********
seq2_should_be_E18             AAAAACGGG..GTTTTTTTTTTAAAA..GGGGGGGG
#=GR seq2_should_be_E18     PP *********..***************..********
seq3_closest_to_CVB5           AAAAACGGG..GAAAAGTTTTTAAAA..GGGGGGGG
#=GR seq3_closest_to_CVB5   PP *********..***************..********
seq4_different_inserts         AAAAACGGGggGCCCCCTTTTTAAAAttGGGGGGGG
#=GR seq4_different_inserts PP ******87511678899*******************
seq5_partial_region            AAAAACGGG..GNNNNNNNNNNAAAA..GGGGGGGG
#=GR seq5_partial_region    PP *********..***************..********
#=GC SS_cons                   :::::::::..:::::::::::::::..::::::::
#=GC RF                        AAAAACGGG..GtttttTTTTTAAAA..GGGGGGGG
#=GC RFCOLX.                   000000000..111111111122222..22222333
#=GC RFCOL.X                   123456789..012345678901234..56789012
//
```

And look at `va-nn/va-nn.vadr.scn`:

```
#seq  seq                     seq                           sub     fract                     sub    fract                 fid  nnregion_seqspan    nnregion  nnregion  seq   
#idx  name                    len  p/f   ant  model   grp1  grp1      id1  seq1         grp2  grp2     id2  seq2          diff        mdl_coords  mdl_coords     covrg  alerts
#---  ----------------------  ---  ----  ---  ------  ----  -----  ------  -----------  ----  ----  ------  ----------  ------  ----------------  ----------  --------  ------
1     seq1_should_be_E30       32  PASS  yes  toy-nn  EVB   E30    1.0000  refseq_E30   EVB   E18   0.5833  refseq_E18  0.4167           7..18:+     7..18:+    1.0000  -     
2     seq2_should_be_E18       32  PASS  yes  toy-nn  EVB   E18    1.0000  refseq_E18   EVB   E30   0.5833  refseq_E30  0.4167           7..18:+     7..18:+    1.0000  -     
3     seq3_closest_to_CVB5     32  PASS  yes  toy-nn  EVB   CV-B5  0.9167  refseq_CVB5  EVB   E30   0.5833  refseq_E30  0.3334           7..18:+     7..18:+    1.0000  -     
4     seq4_different_inserts   36  PASS  yes  toy-nn  EVB   E30    1.0000  refseq_E30   EVB   E18   0.5833  refseq_E18  0.4167           7..18:+     7..18:+    1.0000  -     
5     seq5_partial_region      32  PASS  yes  toy-nn  EVB   E30    0.3333  refseq_E30   EVB   E18   0.3333  refseq_E18  0.0000           7..18:+     7..18:+    1.0000  INDEFINITE_CLASSIFICATION_NN(nnindfcl),LOW_ID_CLASSIFICATION_NN(nnloidcl)
```

**Interpretation:**
- `seq1`: 100% identity (12/12 RF positions) to refseq_E30 in classification region → classified as E30
- `seq2`: 100% identity (12/12 RF positions) to refseq_E18 → classified as E18  
- `seq3`: 91.67% (11/12) identity to refseq_CVB5 (one mismatch at RF position 12) → classified as CV-B5
- `seq4`: **100% identity to refseq_E30 despite different inserts** (gg at insert positions after RF position 9) - demonstrates that insert columns are completely ignored in identity calculation
- `seq5`: Used partial sequence for classification due to Ns, classified as E30 with low difference to E18 (triggers `nnindfcl` and `nnloidcl` alerts)

### Step 7: Understanding the results

**Key observations:**

1. **Fractional identity (`fract id1`)**: Calculated only over non-gap reference positions in the classification region (would have been full sequence if region not used)

2. **Region coverage (`nnregion covrg`)**: Shows what fraction of the classification region the sequence actually spans
   - Value of 1.0 = sequence fully spans the classification region
   - Value < 1.0 = partial coverage (may trigger alerts)

3. **Second nearest neighbor (`seq2`)**: Always from a different subgroup, used to assess classification confidence
   - Large difference (`fid diff`) = confident classification
   - Small difference = ambiguous classification (may trigger `nnindfcl` alert)

4. **Classification coordinates**: 
   - `nnregion_seqspan mdl_coords` shows which model positions the sequence actually covers
   - `nnregion mdl_coords` shows the full classification region defined in the alignment

### Common scenarios and alerts:

**Scenario 1: High confidence classification**
```
fract id1: 0.95, fract id2: 0.70, fid diff: 0.25
```
→ Clear classification, no alerts

**Scenario 2: Low identity to nearest neighbor**
```
fract id1: 0.68, threshold: 0.75
```
→ `nnloidcl` alert: sequence may be a novel serotype or contamination

**Scenario 3: Ambiguous classification**
```
fract id1: 0.85, fract id2: 0.83, fid diff: 0.02, threshold: 0.05
```
→ `nnindfcl` alert: sequence may be recombinant or misclassified

**Scenario 4: Partial region coverage**
```
nnregion covrg: 0.45, threshold: 0.50
```
→ `nnptrgcl` alert: sequence doesn't fully span classification region

**Scenario 5: Classification region not used**
```
nnregion covrg: 0.02 (< 40nt minimum)
nnregion mdl_coords: 1..33:+ (full sequence used instead of 10..25)
```
→ `nnalrgcl` alert: used alternative region for classification

---

## Additional Options

### Disabling nearest-neighbor classification:

```bash
v-annotate.pl --ignore_nnclass test-seqs.fa output-dir
```
This disables nearest-neighbor classification even if `:FILE:` syntax is present in the `.minfo` file.

### Ignoring classification region specifications:

```bash
v-annotate.pl --ignore_nnregion test-seqs.fa output-dir
```
This uses the full sequence for classification even if region boundaries are specified.

### Changing minimum region length:

```bash
v-annotate.pl --nnregionlen 60 test-seqs.fa output-dir
```
Sets minimum subsequence length for region-based classification to 60nt (default: 40).

---

## Tips and Best Practices

1. **Choose classification regions carefully**: Select regions with high variability between subtypes but conservation within subtypes (e.g., VP1 for enteroviruses)

2. **Include representative sequences**: Your reference alignment should include at least one sequence from each subgroup you want to detect

3. **Monitor alert thresholds**: Adjust `--nn_lowidclass` and `--nn_indefclass` based on your data:
   - More divergent virus families: lower thresholds (0.70-0.75)
   - Closely related strains: higher thresholds (0.85-0.90)

4. **Check classification confidence**: Always review sequences with:
   - Low `fract id1` values (< 0.75)
   - Small `fid diff` values (< 0.05)
   - Low `nnregion covrg` values (< 0.50)

5. **Validate with full-length sequences**: Test your classification system with complete genomes before applying to partial sequences

6. **Update reference alignments**: Periodically add newly characterized sequences to improve classification accuracy

---

## References and Additional Information

- Complete output format descriptions: [formats.md](formats.md)
- Alert descriptions and thresholds: [annotate.md](annotate.md#alerts)
- Model building instructions: [build.md](build.md)
- VADR main documentation: [README.md](README.md)

---

[Back to top](#top)
