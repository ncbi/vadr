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
Structure that the template does not declare is never drawn, no matter what the
sequence could form.

One thing is worth knowing before reading further. A diagram covers **only the
model positions the template declares**, and a sequence that does not reach
those positions gets no diagram at all. That is recorded rather than raised as
an error. See [Coverage](#coverage).

`--draw_r2dt` requires an R2DT installation, which `vadr-install.sh` installs
for you. It also requires a model package whose `.minfo` file declares
templates. Most VADR model packages do not, and adding them is a
model-building task covered on a separate page,
[Adding R2DT templates for a VADR model](r2dt-templates.md#top).

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

`vadr-install.sh` installs R2DT for you, along with Traveler and the helper
scripts R2DT's drawing path needs. R2DT is the only optional dependency in that
script: if it cannot be installed, the rest of the VADR installation still
succeeds, and everything except `--draw_r2dt` works. The R2DT step needs
`python3` and network access to PyPI, which the rest of the installation does
not. See [install.md](install.md#r2dt) for the details and for what to do if
that step fails.

**NOTE:** `vadr-install.sh` pins a specific R2DT release tag, named in
`vadr-install.sh` itself in the `R2DTVERSION` variable.

If you would rather install R2DT yourself, follow
[R2DT's own instructions](https://r2dt.readthedocs.io/) and point `R2DT_DIR` at
the result. **VADR does not state a minimum R2DT version**, because it has not
been tested against a range of them. The release `vadr-install.sh` pins is the
one VADR is known to work with.

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

### <a name="runtime"></a>Runtime cost

`--draw_r2dt` runs `r2dt.py` once per (sequence, template) pair that has
residues to draw. Each invocation performs an Infernal alignment of the
extracted residues to the template covariance model and then a Traveler
rendering, so the cost is seconds per pair rather than milliseconds. On an input
of a few thousand sequences with two templates per model, the drawing stage can
easily dominate the total run time. Know this before running it on a large
input; there is no option to draw only a subset, so use a smaller input file if
you want a sample.

### `--keep`

Without `--keep`, VADR keeps only the summary table and the diagrams. With
`--keep` it also keeps `<out_root>.r2dt-input/`, which holds the residues
extracted for each pair, R2DT's own output directory for that pair, and the
captured `r2dt.py` standard output. That directory is scratch data VADR can
regenerate, and it is large, so it is removed at the end of a run unless you ask
for it. It is also where you look when `r2dt.py` fails; see
[Troubleshooting](#troubleshooting).

## Output files<a name="output"></a>

With `--draw_r2dt`, `v-annotate.pl` creates the following, where `<out_root>` is
`<output directory>/<output directory>.vadr`:

```
<out_root>.rdt                                    summary table, one row per (sequence, template)
<out_root>.r2dt-svg/<seq>-<template>.svg          the diagrams
<out_root>.r2dt-input/<seq>-<template>.fa         the extracted residues, --keep only
<out_root>.r2dt-input/<seq>-<template>.r2dt-out/  r2dt.py's own output directory, --keep only
```

For the example command above, that is `va-zika/va-zika.vadr.rdt`,
`va-zika/va-zika.vadr.r2dt-svg/`, and so on.

### <a name="rdt"></a>`.rdt`

A space delimited, column aligned table with one row per (sequence, template)
pair, in the same style as VADR's other per-run tables, with comment lines
naming the columns. **This is the file to look at first.** Real output from a
three sequence run:

```
#seq_id      pass_fail  template_name  r2dt_status  overlaps  covered_ranges                                              covered_pct  output_svg
#----------  ---------  -------------  -----------  --------  ----------------------------------------------------------  -----------  ----------
AY632535.2   FAIL       zika-linear    pass                0  1..73:+,75..137:+,139..210:+,10380..10710:+,10712..10807:+         99.5  va-example.vadr.r2dt-svg/AY632535.2-zika-linear.svg
AY632535.2   FAIL       zika-circular  pass                0  1..73:+,75..137:+,139..190:+,10666..10710:+,10712..10807:+         99.1  va-example.vadr.r2dt-svg/AY632535.2-zika-circular.svg
NC_012532.1  FAIL       zika-linear    pass                0  1..73:+,75..137:+,139..210:+,10380..10710:+,10712..10807:+         99.5  va-example.vadr.r2dt-svg/NC_012532.1-zika-linear.svg
NC_012532.1  FAIL       zika-circular  pass                0  1..73:+,75..137:+,139..190:+,10666..10710:+,10712..10807:+         99.1  va-example.vadr.r2dt-svg/NC_012532.1-zika-circular.svg
KF383047.1   PASS       zika-linear    pass                0  10380..10637:+                                                     40.4  va-example.vadr.r2dt-svg/KF383047.1-zika-linear.svg
KF383047.1   PASS       zika-circular  fail-nocov          -  -                                                                   0.0  -
```

| column | meaning |
|--------|---------|
| `seq_id` | sequence name |
| `pass_fail` | `PASS` or `FAIL`, the sequence's overall VADR pass/fail status, the same status that determines whether it appears in `.pass.list` or `.fail.list`. Both are drawn. |
| `template_name` | the R2DT template, or `-` for a `skipped` row |
| `r2dt_status` | what happened for this pair (see below) |
| `overlaps` | Traveler's count of colliding drawn elements for this diagram, or `-` if no diagram was produced |
| `covered_ranges` | the model (RF) sub-ranges within the template's declared ranges at which this sequence actually has residues, in VADR coords format, or `-` if none |
| `covered_pct` | what percentage of the template's declared length those sub-ranges amount to, or `0.0` if none |
| `output_svg` | path to the SVG, relative to the output directory, or `-` if no diagram was produced |

The `r2dt_status` values say what happened, so you do not have to look anywhere
else to find out why a pair produced nothing:

* **`pass`**: a diagram was produced.
* **`skipped`**: the sequence was classified to a model that has no
  `R2DT_TEMPLATE` lines. Nothing was attempted. There is one such row per
  sequence, not one per template.
* **`fail-nocov`**: the sequence has no residues at the template's declared
  ranges, so there was nothing to draw and `r2dt.py` was not run. This is the
  [coverage](#coverage) rule firing, and it is the common one.
* **`fail-noaln`**: the sequence has no row in the model's alignment.
* **`fail-r2dt`**: `r2dt.py` ran and either exited nonzero or produced no SVG.

**WARNING: a run in which every pair failed still exits with status 0.** Missing
diagrams are reported in this file and nowhere else, so a caller that does not
read it sees only absent files.

**`overlaps` is the field to scan for diagram quality.** It is Traveler's count
of drawn elements that collide with one another. A well matched template and
target give `0`. A nonzero count means the template layout and this particular
sequence disagree enough that the drawing is crowded, and the diagram is worth
looking at before it is used for anything.

**`covered_pct` is the field to scan for coverage.** The example above shows why
both matter: `KF383047.1` is a 3' fragment, so it covers 40.4% of `zika-linear`
and none of `zika-circular`, whose ranges start further along. It gets one
diagram, not two.

### <a name="svg"></a>`.r2dt-svg/<seq>-<template>.svg`

One SVG per drawn (sequence, template) pair, all in one flat directory. Only
pairs that produced a diagram appear, so the file count matches the number of
`pass` rows in the `.rdt` file.

These are the R2DT *colored* SVGs, in which residues are colored according to
how they relate to the template. They are self contained and open in any
browser. Note that some command line SVG rasterizers do not apply the
id-scoped CSS that R2DT uses to set those colors, and will render the diagram in
black; a browser is the reliable way to look at one.

Two examples ship with this documentation, both drawn on the same template:

* [NC_035889.1-zika-linear.svg](r2dt-files/NC_035889.1-zika-linear.svg):
  a full length sequence, `r2dt_status` `pass`, `overlaps` `0`. Every position
  the template covers is present.
* [KF383047.1-zika-linear.svg](r2dt-files/KF383047.1-zika-linear.svg):
  a partial sequence, also `pass` and `overlaps` `0`, but covering only part of
  the template. Compare it against the full length one. The whole first block of
  the template is absent, because this sequence has no residues there, and the
  drawing simply starts where the sequence starts.

### <a name="input"></a>`.r2dt-input/`

**Only written with `--keep`.** For each pair that got as far as being
attempted, this directory holds three things:

```
<seq>-<template>.fa                  the two line FASTA VADR extracted and handed to r2dt.py
<seq>-<template>.r2dt-out/           r2dt.py's own output directory for that pair
<seq>-<template>.r2dt-out.stdout     r2dt.py's captured standard output and standard error
```

The FASTA is the quickest way to understand a coverage failure, since a short
file here **is** the coverage rule firing. The other two are what you need when
`r2dt.py` itself failed, and the `.stdout` file is where R2DT's own error
message will be.

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
would be, so a partial diagram will not overlay a complete one. Diagrams of
sequences that cover the whole template usually do come out the same size, but
not always, because the local adjustments described above can move an outermost
residue and change the crop. Read the diagrams as comparable in shape and scale,
not as overlayable.

**The drawn region is only the template's declared model position ranges.**
A diagram is not a picture of the whole sequence. Everything outside those
ranges is absent by construction, and its absence carries no meaning.

**WARNING: the `5'` and `3'` labels mark the ends of the drawn region, not the ends of
the sequence and not the ends of the genome.** For a partial sequence this is
the one real risk of misreading a diagram. In the partial example above, the
residue labeled `5'` is simply where this sequence's drawn residues begin. It
is not the 5' end of anything.

**The periodic numbering ticks are always submitted-sequence positions.** VADR
rewrites every tick Traveler draws along a diagram from its position within the
drawn, extracted residues to the corresponding position in the sequence exactly
as you submitted it, for every template, whether or not the template's own
layout declares numbering of its own. A tick therefore tells you where you are
in *that sequence record*, not in the genome. For a complete genome numbered
from its own position 1 the two coincide; for a fragment, or a record numbered
some other way relative to the genome, they do not, and the tick is only a
position within the record you gave VADR. This is what resolves the `5'`/`3'`
ambiguity above: checking the nearest tick tells you where the drawn region
starts in the sequence you submitted, and is worth the few seconds it takes.

**A `//` marker with a nucleotide count marks a break between declared
ranges.** When a template declares more than one model position range and a
sequence's drawn residues span the boundary between two of them, VADR draws a
`//` between the two flanking residues, labeled with the exact number of
submitted-sequence nucleotides that lie between them and are not drawn. The
count is a property of that sequence, not of the template: two full length
Zika genomes in one run got `10,169 nt` and `10,157 nt` at the same boundary,
a real 12 nt difference between the two genomes' lengths, not an off-by-one.
A sequence whose drawn residues fall entirely within one declared range gets
no marker, because there is no boundary to cross; the partial-sequence example
above is this case, since its drawn residues never leave the template's second
range. **The count is not recoverable from the [`covered_ranges`](#rdt)
column**, which is in model (RF) coordinates: subtracting those gives the same
gap, `10,169`, for both sequences above, because both reach the same model
positions on either side of the boundary. The two counts differ because the
sequences themselves differ in length across the omitted stretch, something
`covered_ranges` does not report since that stretch is outside every declared
range. Treat the label as VADR's own count, not one you can check by hand from
the `.rdt` file.

For what the residue colors mean, see R2DT's own
[documentation](https://r2dt.readthedocs.io/). The color scheme is R2DT's, and
VADR neither sets nor modifies it.

Finally, use the [`overlaps`](#rdt) column as the machine readable quality
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
* the pair is recorded as `fail-nocov` in the `.rdt` file, with `covered_pct`
  `0.0`;
* **the run continues and `v-annotate.pl` exits with status 0.**

A user who does not read the `.rdt` file sees only missing files.

**What to check, in order:** the `fail-nocov` rows in the `.rdt` file, then the
`covered_ranges` and `covered_pct` columns of the rows that did draw, to see how
much of the template each sequence actually reached. For a `fail-r2dt` row,
rerun with `--keep` and read the pair's `.stdout` file in `.r2dt-input/`.

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
> a set of 1026 partial sequences, 1013 of which passed VADR and 13 of which
> failed, 261 produced at least one diagram and 765 produced none.
>
> The arithmetic follows from the range choice. A model whose template covers
> one contiguous, commonly sequenced region would see a very different split.

### If you want diagrams for partial sequences

Choose ranges that the partial sequences you care about actually cover, or
accept that only sequences spanning the ranges will be drawn. Running
`--draw_r2dt` on a set of partial sequences is supported and produces correct
output; you just need to read the `.rdt` file to know what you got. Choosing a
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

### Non-fatal, recorded and continue

| condition | what happens |
|-----------|--------------|
| the sequence has zero residues at the template's ranges | `fail-nocov` row, `r2dt.py` not run |
| `r2dt.py` exits nonzero, or produces no SVG | `fail-r2dt` row. Rerun with `--keep` and read the pair's `.stdout` file for R2DT's own error message |
| the sequence has no row in the model's alignment | `fail-noaln` row for each of that model's templates |

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
  missing or empty, the pair is recorded as `fail-r2dt` even if `r2dt.py` exited 0.
* With `--keep`, the run directory and the captured stdout survive the run, in
  `.r2dt-input/`, and the stdout is where R2DT's own error message will be.

If `r2dt.py` fails on an input you can draw by hand, the difference is usually
environment: check that `r2dt-vadr-env.sh` provides the same `PATH` your
interactive shell does.

### The run is slow

See [Runtime cost](#runtime). The drawing stage is one `r2dt.py` process per
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
