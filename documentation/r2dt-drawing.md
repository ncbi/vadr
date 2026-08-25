# <a name="top"></a> Drawing R2DT secondary structure figures with `v-annotate.pl`

* [Overview](#overview)
* [Requirements](#requirements)
* [Determining whether a model has R2DT templates](#templates)
* [Running `--draw_r2dt`](#running)
* [Output files](#output)
* [How to read a diagram](#reading)
* [Coverage: why some sequences get no diagram](#coverage)
* [Other failure modes and troubleshooting](#troubleshooting)
* [Adding templates for your own model](#adding)
* [References](#references)

---

## Overview<a name="overview"></a>

[R2DT](https://r2dt.bio/) draws RNA secondary structure diagrams in a
consistent, reproducible layout by mapping a sequence onto a *template*, a
fixed reference layout for a particular structured RNA. The actual drawing is
done by [Traveler](https://github.com/cusbg/traveler), which R2DT calls
internally.

The `v-annotate.pl` option `--draw_r2dt` connects VADR's alignments to R2DT. For
each sequence that VADR classifies to a model that has one or more R2DT
templates declared, VADR extracts the residues the sequence has at the model
positions the template covers, hands them to R2DT, and collects one SVG per
(sequence, template) pair. Because every diagram is drawn on the same template
layout, the same structural element is laid out the same way, at the same
scale, in every sequence's diagram, which is what makes them comparable.

**This is not structure prediction.** The structure shown comes from the
template, not from folding the sequence. What varies between diagrams is which
residues are present, which differ from the template consensus, and how R2DT
colors them. If you want a structure predicted for a sequence, this is the wrong
tool.

Two things worth knowing before reading further:

* Drawing is done for **every sequence classified to a model that has
  templates, whether the sequence passes or fails**. A failing sequence is often
  exactly the one worth looking at. Sequences that are not classified to any
  model are not drawn.
* A diagram covers **only the model positions the template declares**, and a
  sequence that does not reach those positions gets no diagram at all, with a
  warning rather than an error. See [Coverage](#coverage). This catches people
  out, and it is the single most important part of this page.

`--draw_r2dt` requires an R2DT installation, which VADR does not install. It
also requires a model package whose `.minfo` file declares templates. Most
VADR model packages do not, and adding them is a model-building task covered on
a separate page, [Adding R2DT templates for a VADR model](r2dt-templates.md#top).

Throughout this page, examples use a Zika virus model with two templates named
`zika-linear` and `zika-circular`. Those names, and the model positions they
cover, are properties of that particular model. Substitute your own.

## Requirements<a name="requirements"></a>

### The `R2DT_DIR` environment variable

`--draw_r2dt` requires the `R2DT_DIR` environment variable to be set to the root
of an R2DT installation, the directory that contains `r2dt.py`:

```bash
export R2DT_DIR=/path/to/R2DT
```

`v-annotate.pl` checks at startup that `R2DT_DIR` is set and that
`$R2DT_DIR/r2dt.py` exists, and exits with an error naming the problem if
either check fails. Both checks happen before any analysis, so a misconfigured
install fails immediately rather than after the alignment stage.

### The R2DT version

<!-- TODO: name the R2DT version installed by vadr-install.sh here, once
     vadr-install.sh installs R2DT. -->

`vadr-install.sh` does not currently install R2DT; you install it yourself,
following [R2DT's own instructions](https://r2dt.readthedocs.io/).

**VADR does not state a minimum R2DT version**, because it has not been tested
against a range of them. Use a current R2DT release. If a future
`vadr-install.sh` installs R2DT, use the version it installs.

### What R2DT itself needs for this code path

VADR invokes R2DT as `r2dt.py draw --force_template`, which skips R2DT's own
classification stage entirely. A minimal R2DT installation is therefore
sufficient. R2DT ships a `requirements-minimal.txt` listing the Python packages
this path needs, all of them pure Python. Beyond those, the path calls Traveler,
Infernal (`cmalign` and several `esl-` miniapps), and a handful of Bio-Easel and
jiffy Infernal/HMMER scripts. The larger set of dependencies R2DT needs for its
full classification pipeline is not required.

R2DT's [installation documentation](https://r2dt.readthedocs.io/) is the
authoritative source for how to install it, and for what its current
dependencies are. Do not treat the summary above as a substitute.

### The optional site configuration file `r2dt-vadr-env.sh`

An R2DT installation frequently needs environment set up before `r2dt.py` will
run, typically `PATH` additions for Infernal, Bio-Easel, the jiffy scripts and
Traveler, and often a virtualenv `python`. VADR does not hardcode any of this.
Instead, if the file

```
$R2DT_DIR/r2dt-vadr-env.sh
```

exists, `v-annotate.pl` sources it with POSIX `.` immediately before invoking
`r2dt.py`. If the file does not exist, VADR assumes that `python
$R2DT_DIR/r2dt.py` already works in your environment.

Because the file is sourced with `.` from `/bin/sh`, keep it POSIX compatible.
An annotated example is at
[r2dt-files/example-r2dt-vadr-env.sh](r2dt-files/example-r2dt-vadr-env.sh);
copy it to `$R2DT_DIR/r2dt-vadr-env.sh` and edit the paths for your site.

The only environment VADR sets on its own is thread count pinning. It sets
`OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS` and `MKL_NUM_THREADS` to `1` on the
`r2dt.py` command line, so that a numerical library inside R2DT does not
oversubscribe the machine.

### What you can do without an R2DT installation

You can read this page and
[r2dt-templates.md](r2dt-templates.md#top), understand the mechanism, and write
correct `R2DT_TEMPLATE` lines with no R2DT installed. You cannot produce a
diagram without one.

## Determining whether a model has R2DT templates<a name="templates"></a>

Templates are declared in the model package's model info (`.minfo`) file, one
`R2DT_TEMPLATE` line per (model, template) pair. To find out whether a package
declares any:

```bash
grep R2DT_TEMPLATE $VADRMODELDIR/vadr.minfo
```

A declaration looks like this:

```
R2DT_TEMPLATE name=zika-linear model=MG807646 ranges=1..210,10380..10807
```

`name` is the R2DT template name, `model` is the VADR model name, and `ranges`
is a comma separated list of inclusive, 1-indexed model (RF) position ranges
that feed the template. The syntax and its validation rules are documented in
[Adding R2DT templates for a VADR model](r2dt-templates.md#minfo) and in
[formats.md](formats.md#minfo).

Three outcomes are possible:

* **The `.minfo` has no `R2DT_TEMPLATE` lines at all.** `v-annotate.pl
  --draw_r2dt` exits with a fatal error naming the model info file. This is
  the fastest way to learn that a model package does not support the option.
* **The `.minfo` has `R2DT_TEMPLATE` lines, but not for the model your
  sequence is classified to.** Those sequences are recorded in the summary
  file with status `skipped`, and no diagram is drawn for them. This is not an
  error.
* **The `.minfo` declares templates for the relevant model.** Those sequences
  are drawn, subject to [coverage](#coverage).

### <a name="compat"></a>Backward compatibility warning for model package maintainers

**A `.minfo` file containing `R2DT_TEMPLATE` lines cannot be parsed by VADR 1.7
or earlier.** Older versions reject the file as malformed, and the failure is
not limited to drawing. The *entire model package* becomes unusable with that
VADR version, for every script and every sequence.

`R2DT_TEMPLATE` is first accepted by VADR 1.7.1. If you maintain a model package
whose users may still be on an older VADR, do one of the following rather than
adding the lines to the package's main `.minfo`:

* keep the `R2DT_TEMPLATE` lines in a separate copy of the `.minfo` file, in
  its own directory, and have users who want drawing point at it with
  `--mdir`; or
* ship a version of the package with the lines stripped for older VADR users.

## Running `--draw_r2dt`<a name="running"></a>

The general form is an ordinary `v-annotate.pl` command with `--draw_r2dt`
added:

```bash
export R2DT_DIR=/path/to/R2DT
v-annotate.pl --draw_r2dt --mdir <model directory> --mkey <model key> <fasta file> <output directory>
```

`--draw_r2dt` does not require, and is not incompatible with, any other option.
It runs after annotation is complete and does not change any other output.

A worked example, using a Zika model directory:

```bash
export R2DT_DIR=/path/to/R2DT
v-annotate.pl --draw_r2dt --mkey zikv-canonical --mdir /path/to/zika-models \
    zika.fa va-zika
```

### Runtime cost

`--draw_r2dt` runs `r2dt.py` once per (sequence, template) pair that has
residues to draw. Each invocation performs an Infernal alignment of the
extracted residues to the template covariance model and then a Traveler
rendering, so the cost is seconds per pair rather than milliseconds. On an input
of a few thousand sequences with two templates per model, the drawing stage can
easily dominate the total run time. Know this before running it on a large
input; there is no option to draw only a subset, so use a smaller input file if
you want a sample.

### `--keep`

With `--keep`, VADR retains R2DT's full per-pair output tree and the captured
`r2dt.py` standard output. Without it, both are deleted after the SVG has been
copied out, since the tree is large. The extracted input FASTA files and the
copied SVGs are kept either way. Retaining the tree is the main way to
diagnose an `r2dt.py` failure; see [Troubleshooting](#troubleshooting).

## Output files<a name="output"></a>

With `--draw_r2dt`, `v-annotate.pl` creates the following, where `<out_root>` is
`<output directory>/<output directory>.vadr`:

```
<out_root>.r2dt.tsv                             summary table, one row per (sequence, template)
<out_root>.r2dt.warn                            warnings, created only if at least one pair failed
<out_root>.r2dt/<seq>/<seq>-<template>.svg      the diagrams
<out_root>.r2dt-input/<seq>-<template>.fa       the extracted residues handed to r2dt.py
```

For the example command above, that is `va-zika/va-zika.vadr.r2dt.tsv`,
`va-zika/va-zika.vadr.r2dt/`, and so on.

### <a name="tsv"></a>`.r2dt.tsv`

A tab delimited summary, one row per (sequence, template) pair, with a comment
line naming the columns. **This is the file to look at first.** Real output from
a three sequence run:

```
#seq_id	pass_fail	template_name	r2dt_status	overlaps	output_svg
NC_035889.1	PASS	zika-linear	ok	0	va-zika/va-zika.vadr.r2dt/NC_035889.1/NC_035889.1-zika-linear.svg
NC_035889.1	PASS	zika-circular	ok	0	va-zika/va-zika.vadr.r2dt/NC_035889.1/NC_035889.1-zika-circular.svg
KF383047.1	PASS	zika-linear	ok	0	va-zika/va-zika.vadr.r2dt/KF383047.1/KF383047.1-zika-linear.svg
KF383047.1	PASS	zika-circular	fail	-	-
AB908162.1	PASS	zika-linear	fail	-	-
AB908162.1	PASS	zika-circular	fail	-	-
```

| column | meaning |
|--------|---------|
| `seq_id` | sequence name |
| `pass_fail` | `PASS` or `FAIL`, the sequence's overall VADR pass/fail status, the same status that determines whether it appears in `.pass.list` or `.fail.list`. Both are drawn. |
| `template_name` | the R2DT template, or `-` for a `skipped` row |
| `r2dt_status` | `ok`, `fail`, or `skipped` (see below) |
| `overlaps` | Traveler's count of colliding drawn elements for this diagram, or `-` if no diagram was produced |
| `output_svg` | path to the SVG, as `v-annotate.pl` wrote it, or `-` if no diagram was produced |

The three `r2dt_status` values are:

* **`ok`**: a diagram was produced.
* **`fail`**: no diagram was produced for this pair. Some of these are
  [coverage](#coverage) cases and some are genuine `r2dt.py` failures; the
  corresponding line in `.r2dt.warn` says which.
* **`skipped`**: the sequence was classified to a model that has no
  `R2DT_TEMPLATE` lines. Nothing was attempted. There is one such row per
  sequence, not one per template.

**`overlaps` is the field to scan for diagram quality.** It is Traveler's count
of drawn elements that collide with one another. A well matched template and
target give `0`. A nonzero count means the template layout and this particular
sequence disagree enough that the drawing is crowded, and the diagram is worth
looking at before it is used for anything.

### <a name="warn"></a>`.r2dt.warn`

Created only if at least one pair failed. One line per failure, saying which
sequence, which template, and why. Continuing the example above:

```
WARNING: sequence KF383047.1 has zero residues in template zika-circular's RF column range(s); skipping r2dt.py
WARNING: sequence AB908162.1 has zero residues in template zika-linear's RF column range(s); skipping r2dt.py
WARNING: sequence AB908162.1 has zero residues in template zika-circular's RF column range(s); skipping r2dt.py
```

The presence of this file is the signal that some sequences did not get all of
the diagrams you might have expected. **`v-annotate.pl` still exits with status
0.**

### <a name="svg"></a>`.r2dt/<seq>/<seq>-<template>.svg`

One SVG per drawn (sequence, template) pair, in a per-sequence subdirectory. The
subdirectory is created for every drawn sequence, so a sequence for which every
template failed leaves an empty directory behind rather than none.

These are the R2DT *colored* SVGs, in which residues are colored according to
how they relate to the template. They are self contained and open in any
browser. Note that some command line SVG rasterizers do not apply the
id-scoped CSS that R2DT uses to set those colors, and will render the diagram in
black; a browser is the reliable way to look at one.

Two examples ship with this documentation, both drawn on the same template:

* [NC_035889.1-zika-linear.svg](r2dt-files/NC_035889.1-zika-linear.svg):
  a full length sequence, `r2dt_status` `ok`, `overlaps` `0`. Every position the
  template covers is present.
* [KF383047.1-zika-linear.svg](r2dt-files/KF383047.1-zika-linear.svg):
  a partial sequence, also `ok` and `overlaps` `0`, but covering only part of
  the template. Compare it against the full length one. The whole first block of
  the template is absent, because this sequence has no residues there, and the
  drawing simply starts where the sequence starts.

### <a name="input"></a>`.r2dt-input/<seq>-<template>.fa`

The two line FASTA file VADR extracted and handed to `r2dt.py`, one per pair
that got as far as being attempted. These are kept whether or not `--keep` is
used, because inspecting one is the quickest way to understand a failure. A
short or absent file here **is** the coverage rule firing.

## How to read a diagram<a name="reading"></a>

**The layout is the template's, not the sequence's.** Where a residue is drawn
on the page is decided by the template, so the same helix appears in the same
orientation, at the same scale, in every diagram drawn from that template. Small
local adjustments are made where a sequence differs from the template
consensus, but the overall geometry does not move. This comparability is the
entire point of templated drawing.

The SVG canvas, however, is cropped to the residues actually drawn. A partial
sequence's diagram is a smaller image whose origin sits at its own first drawn
residue, not a full sized image with a blank region where the missing residues
would be. Two diagrams from one template are directly comparable in shape and
scale, but they are not overlayable without shifting one of them.

**The drawn region is only the template's declared model position ranges.**
A diagram is not a picture of the whole sequence. Everything outside those
ranges is absent by construction, and its absence carries no meaning.

**WARNING: the `5'` and `3'` labels mark the ends of the drawn region, not the ends of
the sequence and not the ends of the genome.** For a partial sequence this is
the one real risk of misreading a diagram. In the partial example above, the
residue labeled `5'` is at model position 10380, the start of the template's
second range, because that is simply where this sequence's drawn residues
begin. It is not the 5' end of anything. If the template
carries numbering labels tied to genome coordinates, as the example templates
do, those labels resolve the ambiguity, and checking them is worth the few
seconds it takes. Whether a template has such labels is a choice made when the
template is authored; see [r2dt-templates.md](r2dt-templates.md#numbering).

For what the residue colors mean, see R2DT's own
[documentation](https://r2dt.readthedocs.io/). The color scheme is R2DT's, and
VADR neither sets nor modifies it.

Finally, use the [`overlaps`](#tsv) column as the machine readable quality
signal. A diagram with a nonzero overlap count deserves a look before you trust
it.

## Coverage: why some sequences get no diagram<a name="coverage"></a>

**A sequence gets no diagram from a template unless it actually has residues at
that template's declared model position ranges.**

VADR builds a template's input by walking the declared ranges in order, taking
the sequence's aligned residue at each model position, concatenating them, and
stripping gaps. If the sequence covers none of those positions, the result is
empty and there is nothing to draw. If it covers only some of them, the result
is short, and R2DT may or may not manage to draw it.

The failure is quiet:

* `r2dt.py` is not run at all when the extraction is empty, and may fail when
  it is very short;
* a line is written to `.r2dt.warn`;
* the pair is recorded as `fail` in `.r2dt.tsv`;
* **the run continues and `v-annotate.pl` exits with status 0.**

A user who does not read `.r2dt.warn` or `.r2dt.tsv` sees only missing files.

**What to check, in order:** the `fail` rows in `.r2dt.tsv`, then the matching
lines in `.r2dt.warn`, then, for rows whose warning is not about zero residues,
the corresponding file in `.r2dt-input/`.

### Which sequences hit it

Partial sequences, and sequences with large deletions relative to the model, are
the cases that hit the coverage rule. Whether any given partial sequence is
affected depends entirely on which model positions the template covers, so it is
a property of the template, not of VADR.

> **Example, and these numbers are one model's, not a general rate.** The
> example Zika model's two templates each cover a block at the 5' end of the
> genome and a block at the 3' end, because that is where the structured
> elements of a flavivirus genome are. The roughly 10 kb of coding sequence
> between them is covered by neither. A fragment consisting only of coding
> sequence therefore overlaps neither block and gets no diagram from either
> template, no matter how long or how good it is. In one run of that model over
> a set of 1026 partial sequences, 255 produced at least one diagram and 758
> produced none.
>
> The arithmetic follows from the range choice. A model whose template covers
> one contiguous, commonly sequenced region would see a very different split.

### If you want diagrams for partial sequences

Choose ranges that the partial sequences you care about actually cover, or
accept that only sequences spanning the ranges will be drawn. Running
`--draw_r2dt` on a set of partial sequences is supported and produces correct
output; you just need to read `.r2dt.tsv` to know what you got. Choosing a
template's ranges is choosing its coverage, and that decision is discussed
further in [r2dt-templates.md](r2dt-templates.md#decisions).

## Other failure modes and troubleshooting<a name="troubleshooting"></a>

### Fatal errors, all detected at startup

| condition | what you see |
|-----------|--------------|
| `R2DT_DIR` not set | error naming the variable and how to set it |
| `$R2DT_DIR/r2dt.py` does not exist | error naming the path that was checked |
| a declared template directory `$R2DT_DIR/data/local_data/<name>/` does not exist | error naming the directory and the `R2DT_TEMPLATE` line that referenced it |
| the model info file has no `R2DT_TEMPLATE` lines at all | error naming the model info file |
| an `R2DT_TEMPLATE` line has an unrecognized field, or is missing `name=`, `model=` or `ranges=` | parse error quoting the offending line |
| `model=` names a model not in the model info file | parse error quoting the line |
| a range is not of the form `<start>..<end>` | parse error quoting the range and the line |
| a range has `start < 1`, or `end < start`, or `end` greater than the model's length | parse error naming the violation |
| ranges within one line overlap or are not in ascending order | parse error naming the range and the previous range's end |

All of these happen before any sequences are processed.

### Non-fatal, warn and continue

| condition | what happens |
|-----------|--------------|
| the sequence has zero residues at the template's ranges | warning, `fail` row, `r2dt.py` not run |
| `r2dt.py` exits nonzero, or produces no SVG | warning naming the retained stdout file, `fail` row |
| the sequence has no row in the model's alignment | warning, `fail` row for each of that model's templates |

### Reproducing an `r2dt.py` failure by hand

For each pair, `v-annotate.pl` runs the equivalent of:

```sh
cd $R2DT_DIR && \
{ if [ -f $R2DT_DIR/r2dt-vadr-env.sh ]; then . $R2DT_DIR/r2dt-vadr-env.sh; fi; } && \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python $R2DT_DIR/r2dt.py draw --force_template <template_name> \
    <absolute path to .r2dt-input/<seq>-<template>.fa> <absolute path to a per-pair run directory> \
    > <run directory>.stdout 2>&1
```

Points worth noting:

* VADR changes directory into `$R2DT_DIR`, because `r2dt.py` expects to run from
  its installation root. Every path VADR passes is therefore made absolute
  first. If you run the command by hand, do the same.
* `--force_template` bypasses R2DT's classification stage. R2DT does not choose
  the template; VADR names it.
* The file VADR looks for afterwards is
  `<run directory>/results/svg/<seq>-<template>.colored.svg`. If that file is
  missing or empty, the pair is recorded as `fail` even if `r2dt.py` exited 0.
* With `--keep`, the run directory and the captured stdout survive the run, and
  the stdout is where R2DT's own error message will be.

If `r2dt.py` fails on an input you can draw by hand, the difference is usually
environment: check that `r2dt-vadr-env.sh` provides the same `PATH` your
interactive shell does.

### The run is slow

See [Runtime cost](#running). The drawing stage is one `r2dt.py` process per
drawn pair and it is not parallelized.

## Adding templates for your own model<a name="adding"></a>

Adding R2DT drawing support to a model means installing one or more R2DT
templates and adding `R2DT_TEMPLATE` lines to the model's `.minfo` file. It does
not require rebuilding the covariance model or changing any feature
definitions.

It does require a curated secondary structure for the region you want drawn, and
a hands-on layout step in an external editor. The VADR side of the procedure,
and the parts of it that are easy to get silently wrong, are documented in
[Adding R2DT templates for a VADR model](r2dt-templates.md#top).

## References<a name="references"></a>

* **R2DT**: McCann H, Meade CD, Williams LD, Petrov AS, Johnson PZ, Simon AE,
  Hoksza D, Nawrocki EP, Chan PP, Lowe TM, *et al.* R2DT: a comprehensive
  platform for visualizing RNA secondary structure. *Nucleic Acids Research*
  2025;53(4):gkaf032. <https://doi.org/10.1093/nar/gkaf032>
* R2DT documentation: <https://r2dt.readthedocs.io/>
* R2DT repository: <https://github.com/r2dt-bio/R2DT>
* R2DT web application: <https://r2dt.bio/>
* **Traveler**, which renders the diagrams: Elias R, Hoksza D. TRAVeLer: a tool
  for template-based RNA secondary structure visualization. *BMC
  Bioinformatics* 2017;18:487. <https://doi.org/10.1186/s12859-017-1885-4>.
  Repository: <https://github.com/cusbg/traveler>
* **RNAcanvas**, one of the layout editors R2DT supports: Johnson PZ, Simon AE.
  RNAcanvas: interactive drawing and exploration of nucleic acid structures.
  *Nucleic Acids Research* 2023;51(W1):W501-W508.
  <https://doi.org/10.1093/nar/gkad302>. Web application:
  <https://rnacanvas.app>

---

#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.
