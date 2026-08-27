# <a name="top"></a> Adding R2DT templates for a VADR model

* [Scope of this page](#scope)
* [What you need before starting](#prereqs)
* [The two decisions that come first](#decisions)
* [What an R2DT template is: four files](#files)
* [Building the template](#building)
* [Installing the template and wiring the `.minfo`](#minfo)
* [Validation checklist](#checklist)
* [What NOT to look for here](#boundary)
* [Worked example](#example)
* [Optional refinements](#numbering)
* [References](#references)

---

## Scope of this page<a name="scope"></a>

This page covers the VADR side of adding R2DT drawing support to a VADR model:
the four files that make up a template, how to install it, how to declare it in
the model's `.minfo` file, how to test it, and the ways it can silently come out
wrong.

**Creating the layout itself is documented by R2DT, and this page links to that
documentation rather than restating it.** R2DT's template documentation is at
[docs/templates.md](https://github.com/r2dt-bio/R2DT/blob/main/docs/templates.md)
in the R2DT repository, and it is the authoritative source for the editors, the
file formats, and how to submit a template to R2DT upstream. The boundary is
stated explicitly in [What NOT to look for here](#boundary).

If you are trying to *use* templates that already exist rather than make new
ones, you want [Drawing R2DT secondary structure figures with
`v-annotate.pl`](r2dt-drawing.md#top) instead.

**NOTE: this is a hands-on procedure with a manual layout editing step in a
browser.** It is not a one-command build, and the layout step is the part that
takes the longest and is easiest to get wrong.

## What you need before starting<a name="prereqs"></a>

* **An R2DT installation you can write to**, or at least whose
  `data/local_data/` directory you can add a subdirectory or a symlink to. See
  [r2dt-drawing.md](r2dt-drawing.md#requirements).
* **Infernal**, specifically `cmbuild` and `cmemit`. A VADR installation
  already provides these.
* **A curated secondary structure** for the region or regions you want drawn, as
  a Stockholm alignment with an `SS_cons` line. This is the input everything
  else is derived from, and it is the part nobody else can do for you.
* **A browser**, for the layout step.
* A VADR model whose `.minfo` file you can edit.

You do **not** need to rebuild the VADR model's covariance model, and you do not
need to change any feature definitions. A model gains drawing support by adding
lines to its `.minfo` file and installing template directories. Nothing else
about the model changes.

## The two decisions that come first<a name="decisions"></a>

### 1. What to draw

A template is a fixed layout. It should cover a region, or a set of regions,
that is structurally comparable across the sequences you expect to annotate. A
region that is structured in some sequences and absent in others will produce
diagrams that are not comparable, which defeats the purpose.

Nothing caps the number of templates a model may have; one is the common case.
A model may declare several, and each is drawn independently.

### 2. The model position ranges

Each template is fed by one or more inclusive, 1-indexed ranges of the VADR
model's consensus (RF) positions. This is the same coordinate frame VADR uses
for model coordinates everywhere else, so a position you read off a VADR
alignment or an `.alt` file is directly usable.

To find the consensus positions for a region of your model, the practical route
is the model's own alignment together with the `#=GC RFCOLX.` and `#=GC RFCOL.X`
numbering annotation that VADR writes into its output Stockholm files, which
gives the model position of every column as a pair of stacked digit rows. There
are worked examples of reading those lines in
[alerts.md](alerts.md#examples).

The rules the ranges must obey, all enforced when the `.minfo` file is parsed:

* each range is `<start>..<end>` with `start >= 1` and `end >= start`;
* `end` must not exceed the model's consensus length;
* ranges on one line must be in ascending order and must not overlap.

Two consequences that are easy to miss:

* **Multiple ranges are concatenated in the order given, and the template's
  own consensus must be that same concatenation, in that same order.** The
  template does not know what the ranges are; it just receives a string of
  residues. If the template's consensus is not built from the same regions in
  the same order, you get a template that builds, aligns and draws, and draws
  nonsense.
* **Choosing the ranges is choosing the template's coverage.** A sequence
  that has no residues at these positions gets no diagram. If you want partial
  sequences drawn, the ranges have to be positions the partial sequences
  actually cover. See [r2dt-drawing.md](r2dt-drawing.md#coverage) for what the
  user sees when they do not.

## What an R2DT template is: four files<a name="files"></a>

A template is a directory `$R2DT_DIR/data/local_data/<name>/` containing four
files, all named after the directory:

| file | what it is | how R2DT uses it at draw time |
|------|------------|-------------------------------|
| `<name>.cm` | Infernal covariance model of the template consensus | aligns the target sequence to the template with `cmalign` |
| `<name>.sto` | the Stockholm alignment the CM was built from, carrying `SS_cons` | passed as `cmalign --mapali <name>.sto --mapstr`, which threads the template's `SS_cons` onto the aligned target. **This is where the target's drawn structure comes from.** |
| `<name>.fasta` | three line structure FASTA: header, consensus sequence, dot-bracket structure | passed to Traveler as the template structure. **This is where the template's structure comes from.** |
| `<name>.xml` | Traveler format layout: per-residue x/y coordinates | passed to Traveler as the reference layout |

**All four are read on every draw.** In particular the `.sto` is not merely an
authoring artifact you can delete once the CM is built; R2DT reads it at run
time. It is also the authoring source of truth, since it is the only file from
which the `.cm` and the `.fasta` can be regenerated consistently.

The `.xml` carries coordinates and nothing else. It contains no pair
information and no pseudoknot information at all.

### Pseudoknots are not in the covariance model

A covariance model is a stochastic context-free grammar and can only represent
nested pairs. `cmbuild` silently discards the crossing pairs in `SS_cons`, and
the resulting `.cm` carries no record of them. In the example template described
[below](#example), `SS_cons` has 176 nested pairs and 16 crossing pairs, and
`cmbuild` reports exactly 176 basepairs in the model.

**Pseudoknots reach the drawing through the `.sto` and the `.fasta`, not the
`.cm`.** Both retain the WUSS letter annotations (`A`/`a`, `B`/`b`, `C`/`c`, and
so on) verbatim: the `.sto` supplies them to the target through `cmalign
--mapali --mapstr`, and the `.fasta` supplies them to the template through
Traveler.

## Building the template<a name="building"></a>

The steps below are marked by who owns them. Steps marked **[R2DT]** are
documented upstream and are only outlined here.

| # | step | owner |
|---|------|-------|
| 1 | Curate the Stockholm alignment `<name>.sto` for the region(s), with `SS_cons` carrying the nested pairs and, if you have them, WUSS letter annotations for pseudoknots. This file ships in the template directory. | you |
| 2 | `cmbuild --hand <name>.cm <name>.sto` | you |
| 3 | `cmemit -c` the consensus, convert the `SS_cons` WUSS string to dot-bracket **preserving the pseudoknot letter anchors**, and write the three line `<name>.fasta` | you |
| 4 | Obtain an initial layout for that sequence and structure | **[R2DT]** |
| 5 | Hand-adjust the layout in a 2D editor and export an RNA 2D JSON Schema file | **[R2DT]** |
| 6 | Convert the layout JSON to the Traveler XML R2DT reads, as `<name>.xml` | **[R2DT]** |
| 7 | Install the directory at `$R2DT_DIR/data/local_data/<name>/` | you |
| 8 | Smoke test with `r2dt.py draw --force_template <name> <name>.fasta <outdir>` | you |
| 9 | Add the `R2DT_TEMPLATE` line to the model's `.minfo` file | you, VADR side |
| 10 | End-to-end test with `v-annotate.pl --draw_r2dt` on a handful of sequences | you, VADR side |

`cmbuild --hand` is important at step 2: it takes the consensus column
definition from the alignment's `RF` annotation rather than letting `cmbuild`
choose columns by conservation. Without it the template's consensus length can
differ from the length of the region you meant to cover, and the `ranges` you
declare will no longer line up with the template.

### The layout route, in outline

Steps 4 through 6 are R2DT's territory and are documented in
[docs/templates.md](https://github.com/r2dt-bio/R2DT/blob/main/docs/templates.md).
The route used to build the example templates was: run the sequence through the
[R2DT web server](https://r2dt.bio/) to get a starting layout, edit it in
[RNAcanvas](https://rnacanvas.app), export an
[RNA 2D JSON Schema](https://github.com/LDWLab/RNA2D-data-schema/) file, and
convert that to Traveler XML. R2DT also supports the
[Exornata](https://exornata.chemistry.gatech.edu/) editor, and a route that
starts from FASTA or BPSEQ plus a Traveler XML file.

R2DT's `r2dt.py generate-template <layout.json>` will build the `.xml`,
`.fasta`, `.cm` and `.sto` from a layout JSON in one step, and **that is the
simplest route when your template has no pseudoknots.** It derives the
structure from the base pairs recorded in the JSON and writes a dot-bracket
using only `(` and `)`, so a template built this way cannot carry pseudoknots
even if the underlying structure has them. Building the `.sto`, `.cm` and
`.fasta` yourself, as steps 1 through 3 above describe, is what preserves them.

### WARNING: two different WUSS conversions

The same `SS_cons` string is converted **twice, differently**, and confusing the
two is a silent failure. The layout tool cannot lay out crossing pairs, so the
file you hand it must have the pseudoknots removed. The template's shipped
`.fasta` must keep them.

| target file | `(` `[` `{` `<` | `)` `]` `}` `>` | pseudoknot letters `Aa Bb Cc ...` | unpaired `,` `_` `-` `:` `~` `.` |
|-------------|-----------------|-----------------|-----------------------------------|----------------------------------|
| the file uploaded to the layout editor (step 4) | `(` | `)` | **`.`** | `.` |
| the shipped `<name>.fasta` (step 3) | `(` | `)` | **kept verbatim** | `.` |

Nested pairs survive both conversions, including long-range ones. Only the
layout file needs the pseudoknots dropped.

**Get the shipped one wrong and every pseudoknot vanishes from every diagram,
with no error anywhere.** The template still builds, still aligns, and still
draws. Nothing in `.r2dt.tsv` reports it, and it is invisible in a rendered
diagram unless you already know which pseudoknots should be there. This is the
single most expensive mistake on this page, which is why the
[checklist](#checklist) below has an item for it.

## Installing the template and wiring the `.minfo`<a name="minfo"></a>

### Installing

Put the template directory at `$R2DT_DIR/data/local_data/<name>/`. A symlink to
a directory you own works, and keeps the template under your own version
control rather than inside the R2DT installation.

`v-annotate.pl --draw_r2dt` checks that every referenced directory exists at
startup and exits with an error naming the missing one, so a broken symlink is
caught immediately.

### Declaring the template

Add one line per (model, template) pair to the model's `.minfo` file:

```
R2DT_TEMPLATE name=<template_name> model=<vadr_model_name> ranges=<start>..<end>[,<start>..<end>]*
```

The three fields may appear in any order and are all required. `name` must match
the `local_data` directory name, `model` must match a `MODEL` line's name in the
same `.minfo` file, and `ranges` is as described in
[The two decisions](#decisions). Every validation failure is fatal and quotes
the offending line.

An annotated example file is at
[r2dt-files/example-r2dt.minfo](r2dt-files/example-r2dt.minfo). The key is also
documented in [formats.md](formats.md#minfo).

**WARNING: adding `R2DT_TEMPLATE` lines to a `.minfo` file makes the whole model
package unreadable by VADR 1.7 and earlier.** If your users may be on an older
VADR, see the
[backward compatibility warning](r2dt-drawing.md#compat) before you edit a
published package's `.minfo` file.

### Testing

Two tests, in order:

```bash
# 1. R2DT alone, drawing the template's own consensus. Must exit 0 and write an SVG.
cd $R2DT_DIR
python $R2DT_DIR/r2dt.py draw --force_template <name> \
    $R2DT_DIR/data/local_data/<name>/<name>.fasta /tmp/<name>-smoketest

# 2. VADR end to end, on a few sequences you know well.
export R2DT_DIR=/path/to/R2DT
v-annotate.pl --draw_r2dt --mdir <model dir> --mkey <model key> few.fa va-few
```

Then read `va-few/va-few.vadr.r2dt.tsv` and check that the rows you expect are
`ok` and that `overlaps` is `0`. A nonzero overlap count on a sequence that
closely matches the template means the layout and the structure disagree, and
is worth resolving before the template is used for real work.

## Validation checklist<a name="checklist"></a>

Every item here is a real failure mode, and most of them produce a template that
works well enough to look fine.

* **Consensus length equals `SS_cons` length equals the residue count in the
  layout JSON.** Three files, one number. If they disagree, stop.
* **The ranges sum to the template's consensus length.** The `ranges` you
  declare in the `.minfo` file must add up to the number of consensus positions
  the template has. In the example, `1..210` plus `10380..10807` is 210 plus 428
  is 638, and the CM's `CLEN` is 638.
* **The layout's base pair set matches `SS_cons`'s nested pair set, in both
  directions.** A layout editor exports only the pairs it actually drew. A pair
  that was silently dropped during editing renders perfectly well and is very
  easy to miss. Check that every nested pair in `SS_cons` is in the layout and
  that every pair in the layout is in `SS_cons`.
* **A layout-only edit should change coordinates and nothing else.** Merge
  the edited x/y values into the existing template rather than round-tripping
  the whole file through the editor. A whole-file round trip re-derives the
  structure from whatever the editor thinks the pairs are, which is how pairs
  get silently deleted. This applies to any layout editor, not to one in
  particular.
* **End labels are not counted as residues.** Editors commonly draw `5'` and
  `3'` as text elements, and some exports include them in the residue list. One
  or two extra "residues" shifts every downstream index by one or two.
* **The shipped `.fasta` dot-bracket still contains the pseudoknot letter
  anchors, and they match `SS_cons`.** No other check on this list catches a
  template whose pseudoknots were converted away, and no rendered diagram will
  tell you either. See [the two conversions](#building).
* **The smoke test exits 0, and the end-to-end test reports `overlaps` `0`.**
* If you serve the layout JSON to a browser-based editor from a local web
  server, note that a plain static file server may be rejected by the editor
  for cross-origin reasons, and the failure can be silent in the browser rather
  than reported. If an editor will not load a file that looks correct, check
  that first.

## What NOT to look for here<a name="boundary"></a>

This page deliberately does not document the following, because R2DT and
Traveler own them and would go out of date here. This list is the boundary,
stated explicitly rather than implied.

| topic | where it lives |
|-------|----------------|
| creating a layout, and the editors that do it (RNAcanvas, Exornata) | [R2DT docs/templates.md](https://github.com/r2dt-bio/R2DT/blob/main/docs/templates.md) and [docs/editors.md](https://github.com/r2dt-bio/R2DT/blob/main/docs/editors.md) |
| `r2dt.py generate-template`, `generatecm`, `generatemodelinfo` | [R2DT docs/templates.md](https://github.com/r2dt-bio/R2DT/blob/main/docs/templates.md) |
| the `data/local_data` convention itself, and R2DT's template library | [R2DT documentation](https://r2dt.readthedocs.io/) |
| submitting a template upstream so others can use it | [R2DT docs/templates.md](https://github.com/r2dt-bio/R2DT/blob/main/docs/templates.md) |
| R2DT's colour scheme and what its own outputs mean | [R2DT documentation](https://r2dt.readthedocs.io/) |
| the Traveler intermediate XML format, and how overlaps are counted | [Traveler repository](https://github.com/cusbg/traveler) |
| installing R2DT | [R2DT documentation](https://r2dt.readthedocs.io/) |

What VADR owns, and what this page and
[r2dt-drawing.md](r2dt-drawing.md#top) therefore do document, is the
`--draw_r2dt` option and the `R2DT_DIR` contract, the `R2DT_TEMPLATE` key and its
validation, the model position ranges and how extraction works, the coverage
rule, VADR's own output files, and the optional `r2dt-vadr-env.sh` site
configuration.

## Worked example<a name="example"></a>

Everything in this section is an example. A model with one contiguous structured
region would declare a single range and a single template, and the rest of this
page would apply unchanged.

A Zika virus model, VADR model name `MG807646`, consensus length 10807, declares
two templates:

```
R2DT_TEMPLATE name=zika-linear   model=MG807646 ranges=1..210,10380..10807
R2DT_TEMPLATE name=zika-circular model=MG807646 ranges=1..190,10666..10807
```

Points this example illustrates:

* **Two ranges, not one, because this genome's structured elements are at its
  two ends.** A flavivirus genome has structured elements in the 5' region and
  in the 3' untranslated region, and roughly 10 kb of coding sequence between
  them with nothing to draw. The block-at-each-end shape is a property of this
  genome, not of the mechanism.
* **Two templates for one model.** `zika-linear` draws the genome laid out as a
  line; `zika-circular` covers narrower blocks and lays them out so that the
  long-range pairing between the 5' and 3' ends is visible, 33 of its 92 pairs
  being between the two blocks. They are independent templates
  with independent ranges, and a sequence can be drawn by one and not the other.
  Two is a choice; one is the common case.
* **The ranges sum to the consensus length.** `zika-linear` covers 210 plus 428
  is 638 positions, and its CM has `CLEN` 638.
* **Coverage follows from the ranges.** `zika-circular` reaches only back to
  model position 10666, so a 3' fragment that stops before that gets no circular
  diagram even though it gets a linear one. That is exactly the `KF383047.1` row
  in the example `.r2dt.tsv` in
  [r2dt-drawing.md](r2dt-drawing.md#tsv).

The two example SVGs shipped with this documentation,
[NC_035889.1-zika-linear.svg](r2dt-files/NC_035889.1-zika-linear.svg) and
[KF383047.1-zika-linear.svg](r2dt-files/KF383047.1-zika-linear.svg), are both
drawn on `zika-linear`.

The two templates were built from one curated alignment: the circular
template's `.sto` was derived by slicing the linear one down to the narrower
ranges. If you build several templates over overlapping regions, deriving them
from a single curated alignment rather than curating each separately is worth
the effort, because it keeps the structures consistent between them.

## Optional refinements<a name="numbering"></a>

Both of these are optional, both are template properties chosen by whoever
authors the template, and both are more example than mechanism.

### Genome-coordinate numbering labels

Traveler draws periodic position labels, and the numbers it draws come from
`numbering-label` attributes on the `<point>` elements of the layout XML. If you
set them to coordinates in your model's frame rather than to positions within
the template, every diagram carries the information a reader
needs to place it. This is what stops a partial sequence's diagram from reading
as though it covers the whole molecule, which is the main misreading risk
described in [r2dt-drawing.md](r2dt-drawing.md#reading). If your templates cover
non-contiguous regions, this is worth doing.

### Colouring specific residue groups

R2DT's colored output supports a fixed set of text colour classes. If the
template's layout already carries them, specific residue groups can be
highlighted, which is useful for showing which residues pair with which across
a pseudoknot. The constraint is that only classes already present in the
template can be used; you cannot introduce a new one from the VADR side.

## References<a name="references"></a>

* **R2DT**: McCann H, Meade CD, Williams LD, Petrov AS, Johnson PZ, Simon AE,
  Hoksza D, Nawrocki EP, Chan PP, Lowe TM, *et al.* R2DT: a comprehensive
  platform for visualizing RNA secondary structure. *Nucleic Acids Research*
  2025;53(4):gkaf032. <https://doi.org/10.1093/nar/gkaf032>
* R2DT documentation: <https://r2dt.readthedocs.io/>
* R2DT repository: <https://github.com/r2dt-bio/R2DT>
* R2DT template creation documentation:
  <https://github.com/r2dt-bio/R2DT/blob/main/docs/templates.md>
* R2DT web application: <https://r2dt.bio/>
* **Traveler**: Elias R, Hoksza D. TRAVeLer: a tool for template-based RNA
  secondary structure visualization. *BMC Bioinformatics* 2017;18:487.
  <https://doi.org/10.1186/s12859-017-1885-4>. Repository:
  <https://github.com/cusbg/traveler>
* **RNAcanvas**: Johnson PZ, Simon AE. RNAcanvas: interactive drawing and
  exploration of nucleic acid structures. *Nucleic Acids Research*
  2023;51(W1):W501-W508. <https://doi.org/10.1093/nar/gkad302>. Web
  application: <https://rnacanvas.app>
* **RNA 2D JSON Schema**, the layout interchange format:
  <https://github.com/LDWLab/RNA2D-data-schema/>
* **Infernal**, for `cmbuild` and `cmemit`: <http://eddylab.org/infernal/>

---

#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.
