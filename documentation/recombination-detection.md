# `v-annotate.pl` recombination detection

* [Overview](#overview)
* [How It Works](#how-it-works)
* [Options](#options)
* [Output](#output)
* [Toy Example](#toy-example)
* [Real-World Example](#real-world)
* [Caveats and Current Status](#caveats)

---

## Overview<a name="overview"></a>

VADR's experimental recombination detection feature is an extension of
[nearest-neighbor (NN) classification](nn-classification.md). When enabled
with `--do_rc`, it scans the alignment for a **breakpoint** position where the
query sequence transitions from being most similar to one reference subgroup
(the *left parent*) to being most similar to a different reference subgroup
(the *right parent*). If the evidence for a switch is strong enough, a
non-fatal `recombin` / `POSSIBLE_RECOMBINATION` alert is reported.

This feature is **experimental and off by default**. See
[Caveats and Current Status](#caveats).

---

## How It Works<a name="how-it-works"></a>

Recombination detection runs after NN classification has been performed and
the alignment of the query sequence to its best-matching model is available.
It proceeds in two passes over the alignment columns of the NN classification
region.

### Log-likelihood ratio (LLR) scoring

At each aligned position the algorithm computes a **per-position LLR score**
comparing a homology null model vs. a random-sequence null model:

```
  match:    score_i = log2(rc_match / (2 * freq_i^2))
  mismatch: score_i = log2(rc_mismatch / (2 * freq_seq_i * freq_mdl_i))
```

where `rc_match` is the assumed match probability (default: 0.95),
`rc_mismatch = (1 - rc_match) / 3` (≈ 0.017), and `freq_i` is the observed
nucleotide frequency at that alignment column from the model's seed alignment.
Positions where either the query or reference has a gap are skipped.
Match scores are positive; large mismatches are negative.

The **cumulative forward score** through position `k` is the sum of
per-position scores for the query vs. a given reference sequence from the
left end of the region up to `k`. The **cumulative backward score** from
position `k` to the right end is computed analogously.

### Breakpoint scan

For every pair of reference sequences belonging to **different subgroups**
(parent-L and parent-R), the algorithm finds the breakpoint `k` that maximizes:

```
  recomb_score(k) = fwd_LLR(parent-L, left_of_k) + bck_LLR(parent-R, right_of_k)
```

A `recombin` alert is reported when the best recomb_score exceeds:

```
  rc_thresh * (npos_left + npos_right)
```

where `rc_thresh` is the minimum required per-position improvement (default:
0.05 bits/position) and `npos_left`, `npos_right` are the non-gap position
counts on each side. Both sides must have at least `rc_minlen` (default: 10)
non-gap positions.

---

## Options<a name="options"></a>

All options are grouped under "options for experimental recombination detection"
in the `v-annotate.pl` usage output.

| Option | Default | Description |
|--------|---------|-------------|
| `--do_rc` | off | Enable recombination detection |
| `--rc_thresh <x>` | 0.05 | Min per-position LLR improvement to report recombination; increase to be more conservative |
| `--rc_match <x>` | 0.95 | Assumed per-position match probability for the homology model. Overridden by `VADR-DEFAULT-RC_MATCH` in the model info STK if `--rc_match` is not set explicitly |
| `--rc_minlen <n>` | 10 | Min non-gap positions required on each side of the breakpoint |
| `--rc_igself` | off | Ignore the query sequence itself as a parent candidate (useful when the query is also in the model alignment) |
| `--rc_iglist <s>` | — | Comma-separated list of group.subgroup strings to skip as parent candidates |

---

## Output<a name="output"></a>

### Alert in `.vadr.alt`

A `recombin` alert appears as a **sequence-level** (non-feature) row in the
`.vadr.alt` file. It always has `fail: no` — it is non-fatal and is flagged
for review only. The columns `seq coords` and `mdl coords` each give the
single-position breakpoint location in sequence and model (RF) coordinates.

**Example row** (wrapped for readability):

```
#idx   seq name    model  type  name  idx  code      fail  description           seq coords    mdl coords  detail
----   ----------  -----  ----  ----  ---  --------  ----  --------------------  ----------    ----------  ------
1.1.1  MZ268661.1  hrvA   -     -     -    recombin  no    POSSIBLE_RECOMBINATION 5250..5250:+  5331..5331:+ possible recombination detected in sequence [...]
```

### Detail string

The bracketed detail string has this structure:

```
[<L_gsg>,<seq_bp>,<mdl_bp>,<R_gsg>;L:<L_name>;id:<L_fid>(+<L_margin>);llr:<L_ppb>(+<L_ppb_margin>),R:<R_name>;id:<R_fid>(+<R_margin>);llr:<R_ppb>(+<R_ppb_margin>),S:<score>]
```

| Field | Meaning |
|-------|---------|
| `<L_gsg>` | group.subgroup of the left parent |
| `<seq_bp>` | Breakpoint position in input sequence coordinates  |
| `<mdl_bp>` | Breakpoint position in model (RF) coordinates |
| `<R_gsg>` | group.subgroup of the right parent |
| `L:<name>` | Accession of the left parent reference sequence |
| `id:<f>(+<d>)` | Fractional identity of query to left parent on the 5' segment; `+<d>` is the margin over using the right parent on that same segment |
| `llr:<s>(+<d>)` | Per-position (bits/nt) LLR score of query vs. left parent on the 5' segment; `+<d>` is the improvement over using the right parent on that segment |
| `R:<name>` | Accession of the right parent reference sequence |
| (id, llr for R) | Same fields for the right parent on the 3' segment |
| `S:<s>` | Total recombination score (bits/nt) averaged over both segments |

### `.vadr.scn` alert column

The `.vadr.scn` file's rightmost `seq alerts` column will show
`POSSIBLE_RECOMBINATION(recombin)` for any sequence with a recombination
alert. Other NN alerts such as `INDEFINITE_CLASSIFICATION_NN(nnindfcl)` often
co-appear for recombinant sequences because the overall nearest-neighbor
classification is ambiguous.

---

## Toy Example<a name="toy-example"></a>

This example uses the toy model files from `testing-20260109/` in the
notebook directory. The model has three reference sequences:

| Seq | Subgroup | Approximate identity to others |
|-----|----------|---------------------------------|
| Seq1 | A.1 | 90% to Seq2, 80% to Seq3 |
| Seq2 | A.2 | 90% to Seq1, 85% to Seq3 |
| Seq3 | B.1 | 80% to Seq1, 85% to Seq2 |

Test sequences in `test-rc.fa` include pure (non-recombinant) sequences and
designed recombinants.

**Command:**
```bash
v-annotate.pl -f --do_rc --mdir . --mkey toy-rc test-rc.fa va-doc-toy-rc
```

**Non-recombinant control (`test1_nonrecomb`, 95% match to Seq1):**

In `va-doc-toy-rc.vadr.scn`:
```
1  test1_nonrecomb  100  PASS  yes  toy-rc  A  A.1  0.9500  Seq1  A  A.2  0.8600  Seq2  0.0900  1..100:+  1..100:+  1.0000  -
```
No alert. Classified as A.1, 95% identical to Seq1, 9 percentage points above
the next-best subgroup (A.2). No breakpoint found.

**Recombinant (`recomb_seq1_seq3_50`, Seq1[1-50] + Seq3[51-100]):**

In `va-doc-toy-rc.vadr.scn`:
```
5  recomb_seq1_seq3_50  100  PASS  yes  toy-rc  B  B.1  0.8700  Seq3  A  A.1  0.8600  Seq1  0.0100  1..100:+  1..100:+  1.0000  INDEFINITE_CLASSIFICATION_NN(nnindfcl),POSSIBLE_RECOMBINATION(recombin)
```

In `va-doc-toy-rc.vadr.alt`:
```
5.1.2  recomb_seq1_seq3_50  toy-rc  -  -  -  recombin  no  POSSIBLE_RECOMBINATION  58..58:+  1  58..58:+  1
  possible recombination detected in sequence
  [A.A.1,58,58,B.B.1;L:Seq1;id:0.931(+0.138);llr:0.445(+0.649),R:Seq3;id:0.976(+0.214);llr:0.807(+1.417),S:2.066]
```

**Interpretation:**

- The sequence was designed as a 50/50 Seq1+Seq3 chimera, but the breakpoint
  was found at position 58 (instead of 50) due to the ~5% noise added to both
  halves — the algorithm finds the position that maximizes the LLR split, not
  necessarily the designed breakpoint.
- Left parent: Seq1 (subgroup A.1), 93.1% identity on positions 1–58,
  13.8 percentage points above Seq3 on that segment.
- Right parent: Seq3 (subgroup B.1), 97.6% identity on positions 59–100,
  21.4 percentage points above Seq1 on that segment.
- The overall NN classification called B.1 (best over the full sequence) — the
  recombinant contains more B.1-like sequence, but is flagged because the 5'
  half clearly resembles A.1.
- `nnindfcl` co-appears because the whole-sequence classification is close
  (B.1 at 0.870 vs A.1 at 0.860, difference 0.010 < threshold 0.050).

**Recombinant at an earlier breakpoint (`recomb_seq1_seq3_25`, Seq1[1-25] + Seq3[26-100]):**

In `va-doc-toy-rc.vadr.alt`:
```
10.1.1  recomb_seq1_seq3_25  toy-rc  -  -  -  recombin  no  POSSIBLE_RECOMBINATION  23..23:+  1  23..23:+  1
  possible recombination detected in sequence
  [A.A.1,23,23,B.B.1;L:Seq1;id:0.913(+0.217);llr:0.592(+1.051),R:Seq3;id:0.961(+0.143);llr:0.568(+0.976),S:2.027]
```

Same parents detected (A.1 left, B.1 right), breakpoint found near position
23, close to the designed breakpoint at 25.

---

## Real-World Example<a name="real-world"></a>

Goya et al. (2024, *J Infect Dis* 231:e154–e164) identified a clade of HRV-A
sequences from Washington State as A105/A21 inter-serotype recombinants, with
a breakpoint near nucleotide 5250 (95% CI: 5216–5271) in the 3C protease
gene. MZ268661.1 (strain RvA105/USA/2021/AHU4DQ) is the primary example. The
parents are:

| Parent | Accession | Serotype | Genome region |
|--------|-----------|----------|---------------|
| Left (5' side) | MZ542285.3 | A105 | VP4 through most of 3C |
| Right (3' side) | JN837693.1 | A21 | Remainder of 3C, 3D polymerase |

**Command** (run from `testing-20260120/`):
```bash
v-annotate.pl -f --do_rc --rc_igself --mdir ../vadr-models-hrv/hrvA --mkey hrvA MZ268661.fa va-doc-MZ268661
```

(`--rc_igself` prevents MZ268661.1 from matching itself if it is present in
the model alignment.)

**In `va-doc-MZ268661.vadr.scn`:**
```
1  MZ268661.1  7051  PASS  yes  hrvA  hrvA  A105  0.9276  MZ542285.3  hrvA  A57  0.8650  FJ445141.1  0.0626  44..7134:+  1..7240:+  0.9794  POSSIBLE_RECOMBINATION(recombin)
```

**In `va-doc-MZ268661.vadr.alt`:**
```
1.1.1  MZ268661.1  hrvA  -  -  -  recombin  no  POSSIBLE_RECOMBINATION  5250..5250:+  1  5331..5331:+  1
  possible recombination detected in sequence
  [hrvA.A105,5250,5331,hrvA.A21;L:MZ542285.3;id:0.956(+0.134);llr:0.780(+0.763),R:JN837693.1;id:0.946(+0.104);llr:0.656(+0.559),S:1.322]
```

**Interpretation:**

- Breakpoint at sequence position 5250, model (RF) position 5331 — falls
  within the 5216–5271 confidence interval from Goya et al.
- Left parent: MZ542285.3 (A105), 95.6% identical on the 5' segment,
  13.4 percentage points above the right parent on that segment.
- Right parent: JN837693.1 (A21), 94.6% identical on the 3' segment,
  10.4 percentage points above the left parent on that segment.
- The alert is non-fatal (`fail: no`).
- MZ268661.1 is classified as nearest to MZ542285.3 overall, consistent with
  the fact that the A105-like region comprises the majority of the genome.

---

## Caveats and Current Status<a name="caveats"></a>

- **Experimental, off by default.** The `recombin` alert is always non-fatal.
  The feature needs broader validation across more viruses before it could be
  enabled in production workflows.
- **Coverage-dependent.** Detection quality depends on having representative
  sequences for both parental subgroups in the model alignment. Missing or
  unannotated serotypes may cause the wrong parents to be reported or
  breakpoints to be missed entirely.
- **Only inter-subgroup recombination is detected.** Parent candidates must
  belong to different annotated subgroups. Intra-subgroup recombination is
  not reported.
- **Classification region only.** Detection operates within the NN
  classification region (e.g., VP1 for enteroviruses). Recombination outside
  that region is invisible to this algorithm.
- **No performance penalty when disabled.** With `--do_rc` off (the default),
  the code path is identical to running without the feature — only forward LLR
  scores over the classification region are computed, needed anyway for NN
  classification.

---

## See Also

- [Nearest-neighbor classification](nn-classification.md)
- [Options reference](annotate.md#options-recomb)
- [Alert descriptions](annotate.md#alerts)
