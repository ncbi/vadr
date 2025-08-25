# VADR - Viral Annotation DefineR <a name="top"></a>
#### Version 1.6.4; June 2024
#### https://github.com/ncbi/vadr.git

VADR is a suite of tools for classifying and analyzing sequences
homologous to a set of reference models of viral genomes or gene
families. It includes models that can be used to validate and annotate
Norovirus, Dengue virus, SARS-CoV-2 virus as well as other
flaviviruses, caliciviruses, and coronaviruses, plus influenza virus,
mpox virus, and respiratory syncitial virus (RSV). Additional models
are [available to download](#models) or can be created using the
`v-build.pl` program.

---
## Quick-start: install VADR and classify and annotate viral sequences using `v-scan.pl`<a name="quickstart"></a>

#### Install VADR:

Download this file:

```
https://raw.githubusercontent.com/ncbi/vadr/master/vadr-install.sh
```

possibly with a command like:
```
curl -o vadr-install.sh https://raw.githubusercontent.com/ncbi/vadr/master/vadr-install.sh
```

And execute it, with one of the following commands depending on your
system type:

```
sh ./vadr-install.sh linux
```

OR

```
sh ./vadr-install.sh macosx-silicon
```

OR

```
sh ./vadr-install.sh macosx-intel
```

Then follow the instructions output at the end of the installation for
updating your `.bashrc` or `.cshrc` file and defining important
environment variables that VADR relies on.

#### Run `v-scan.pl` to annotate viral sequences

Given a fasta sequence file called `my.fa` with any combination of flavivirus,
calicivirus, coronavirus, influenza, RSV, or Mpox sequences, run:

`v-scan.pl -m in.fa out`

This will list each stage of the processing and
ultimately create an output directory called `out` and fill it with 
output files. Short descriptions of the output files will be printed to the
screen. More detailed explanation of output file types can be found
[here](formats.md#annotate). For a more detailed walk-through example
of `v-scan.pl` see [this page](scan.md#longwalk).

---
## VADR programs

The VADR `v-scan.pl` script classifies and annotates sequences that
match to any of your VADR model libraries. 
Once `v-scan.pl` determines the library to use for a given set of
sequences, it runs a different VADR program called `v-annotate.pl`
which identifies the appropriate model in the library to use for each
sequence and defines the annotation based on that most similar model.
`v-scan.pl` will automatically run `v-annotate.pl` using the
recommended settings (`v-annotate.pl` command-line options) for each
library but alternatively, users can run the `v-annotate.pl`
separately. Example usage of `v-annotate.pl` can be found
[here](documentation/annotate.md#top).

Another VADR script, `v-build.pl`, is used to create the models from
individual sequences from GenBank or from input multiple sequence
alignments, potentially with secondary structure
annotation. `v-build.pl` stores the GenBank feature annotation in the
model, and `v-annotate.pl` maps that annotation (e.g. CDS coordinates)
onto the sequences it annotates.  Example usage of `v-build.pl` can be
found [here](documentation/build.md#top). An advanced tutorial on
building VADR models using RSV as an example can be found
[here](documentation/advbuild.md#top).

`v-annotate.pl` identifies unexpected or divergent attributes of the
sequences it annotates (e.g. invalid or early stop codons in CDS
features) and reports them to the user in the form of *alerts*.  A
subset of alerts are *fatal* and cause a sequence to *fail*. A
sequence *passes* if zero fatal alerts are reported for it.  VADR is
used by GenBank staff to evaluate incoming sequence submissions of
some viruses (currently Norovirus, Dengue virus, and SARS-CoV-2).
Submitted Norovirus, Dengue virus and SARS-CoV-2 sequences that pass
`v-annotate.pl` are accepted into GenBank.

The homology search and alignment components of VADR scripts, the most
computationally expensive steps, are performed by the Infernal, HMMER,
FASTA, MINIMAP2 and BLAST software packages, which are downloaded and installed
with [VADR installation](documentation/install.md#top).

---

## VADR model libraries <a name="models"></a>

[VADR installation](install.md#top) includes the following model libraries:

| library      | model key (short name) | rigorously tested? | number of models | notes | 
|--------------|------------------------|--------------------|------------------|-------|
| *Caliciviridae*      | calici | norovirus models only      | 49  | norovirus models used by GenBank |
| *Flaviviridae*       | flavi  | dengue and HCV models only | 156 | dengue models used by GenBank |
| *Coronaviridae*      | corona | SARS-CoV-2 only            | 55  | SARS-CoV-2 models used by GenBank |
| influenza            | flu    | yes                        | 70  | [described in Database article](https://pubmed.ncbi.nlm.nih.gov/39297389/) |
| Mpox                 | mpxv   | yes                        | 1   | | 
| respiratory syncitial virus (RSV) | rsv    | yes                        | 2   | | 

Additional models are available. See [this
page](https://github.com/ncbi/vadr/wiki/Available-VADR-model-files)
for a list of all available models and additional information.

---
## VADR documentation <a name="documentation"></a>

* [VADR installation instructions](documentation/install.md#top)
  * [Installation using `vadr-install.sh`](documentation/install.md#install)
  * [Setting environment variables](documentation/install.md#environment)
  * [Verifying successful installation](documentation/install.md#tests)
  * [Further information](documentation/install.md#further)
* [`v-build.pl` example usage and command-line options](documentation/build.md#top)
  * [`v-build.pl` example usage](documentation/build.md#exampleusage)
  * [`v-build.pl` command-line options](documentation/build.md#options)
  * [Building a VADR model library](documentation/build.md#library)
  * [How the VADR 1.0 model library was constructed](documentation/build.md#1.0library)
* [`v-annotate.pl` example usage, command-line options and alert information](documentation/annotate.md#top)
  * [`v-annotate.pl` example usage](documentation/annotate.md#exampleusage)
  * [Running `v-annotate.pl` inside the `v-scan.pl` wrapper](documentation/annotate.md#scan)
  * [`v-annotate.pl` command-line options](documentation/annotate.md#options)
  * [Basic Information on `v-annotate.pl` alerts](documentation/annotate.md#alerts)
  * [Additional information on `v-annotate.pl` alerts](documentation/annotate.md#alerts2)
* [`v-scan.pl` example usage and command-line options](documentation/scan.md#top)
  * [Quickstart `v-scan.pl` examples](documentation/scan.md#quickstart)
  * [The `v-scan.pl` config file](documentation/scan.md#config)
  * [Walk-throughs of `v-scan.pl` examples](documentation/scan.md#longwalk)
  * [`v-scan.pl` command-line options](documentation/scan.md#options)
* [***Advanced tutorial: building an RSV model library***](documentation/advbuild.md#top)
* [Explanations and examples of `v-annotate.pl` detailed alert and error messages](documentation/alerts.md#top)
  * [Output fields with detailed alert and error messages](documentation/alerts.md#files)
  * [Explanation of sequence and model coordinate fields in `.alt` files](documentation/alerts.md#coords)
  * [`toy50` toy model used in examples of alert messages](documentation/alerts.md#toy)
  * [Examples of different alert types and corresponding `.alt` output](documentation/alerts.md#examples)
  * [Posterior probability annotation in VADR output Stockholm alignments](documentation/alerts.md#pp)
* [VADR output file formats](documentation/formats.md#top)
  * [VADR output files created by all VADR scripts](documentation/formats.md#generic)
  * [`v-build.pl` output files](documentation/formats.md#build)
  * [`v-annotate.pl` output files](documentation/formats.md#annotate)
  * [VADR `coords` coordinate string format](documentation/formats.md#coords)
  * [VADR sequence naming conventions](documentation/formats.md#seqnames)
* [Available VADR model files (github wiki)](https://github.com/ncbi/vadr/wiki/Available-VADR-model-files)
* [SARS-CoV-2 annotation (github wiki)](https://github.com/ncbi/vadr/wiki/Coronavirus-annotation)
* [Rfam-based structural annotation of a viral genome sequence for use with VADR (github wiki)](https://github.com/ncbi/vadr/wiki/Rfam-based-structural-annotation-of-a-viral-genome-sequence)
* [Development notes and instructions (github wiki)](https://github.com/ncbi/vadr/wiki/Development-notes-and-instructions)

---
## Contributors <a name="contributors"></a>
* VADR includes contributions and input from current and former
  colleagues at NCBI, including:

  Rodney Brister
  
  Vince Calhoun
  
  Sergiy Gotvyanskyy
  
  Eneida Hatcher
  
  Sophia Hu
  
  Ilene Karsch-Mizrachi
  
  Rich McVeigh
  
  Susan Schafer
  
  Alejandro Schäffer
  
  Lara Shonkwiler
  
  Beverly Underwood
  
  Yuri Wolf
  
  Linda Yankie

---
## Reference <a name="reference"></a>

* The recommended citation for influenza analysis using VADR is:
  *Vincent C Calhoun, Eneida L Hatcher, Linda Yankie, Eric P Nawrocki; 
  Influenza sequence validation and annotation using VADR. Database.
  baae091. (2024).* https://doi.org/10.1093/database/baae091

* The recommended citation for using VADR for SARS-CoV-2 analysis:
  *Eric P Nawrocki; Faster SARS-CoV-2 sequence validation and
  annotation for GenBank using VADR. NAR Genom Bioinform. 2023 Jan
  20;5(1)::lqad002. (2023).* https://doi.org/10.1093/nargab/lqad002

* The recommended citation for all other uses of VADR is:
  *Alejandro A Schäffer, Eneida L Hatcher, Linda Yankie, Lara Shonkwiler,
  J Rodney Brister, Ilene Karsch-Mizrachi, Eric P Nawrocki; VADR:
  validation and annotation of virus sequence submissions to
  GenBank. BMC Bioinformatics 21, 211
  (2020).* https://doi.org/10.1186/s12859-020-3537-3

---
#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.
