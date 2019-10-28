# VADR - Viral Annotation DefineR <a name="top"></a>
#### Version 0.991dev; October 2019
#### https://github.com/nawrockie/vadr.git

VADR is a suite of tools for classifying and analyzing sequences
homologous to a set of reference models of viral genomes or gene
families. It has been mainly tested for analysis of Norovirus and
Dengue virus sequences in preparation for submission to the GenBank
database. 

The VADR `v-annotate.pl` script is used to classify a sequence, by
determining which in a set of reference models it is most similar
to, and then annotate that sequence based on that most similar model.
Example usage of `v-annotate.pl` can be found [here](annotate.md#top).
Another VADR script, `v-build.pl`, is used to create the models from
NCBI RefSeq sequences or from input multiple sequence alignments,
potentially with secondary structure annotation. `v-build.pl` stores
the RefSeq feature annotation in the model, and `v-annotate.pl` maps
that annotation (e.g. CDS coordinates) onto the sequences it
annotates.  VADR includes 197 prebuilt models of *Flaviviridae* and
*Caliciviridae* viral RefSeq genomes, described
[here](build.md#1.0library).  Example usage of `v-build.pl` can be
found [here](build.md#top).

`v-annotate.pl` identifies unexpected or divergent attributes of the
sequences it annotates (e.g. invalid or early stop codons in CDS
features) and reports them to the user in the form of *alerts*.  A
subset of alerts are *fatal* and cause a sequence to *fail*. A
sequence *passes* if zero fatal alerts are reported for it.  VADR is
used by GenBank staff to evaluate incoming sequence submissions of
some viruses (currently Norovirus and Dengue virus).  Submitted
sequences that pass `v-annotate.pl` are accepted into GenBank.

The homology search and alignment components of VADR scripts, the most
computationally expensive steps, are performed by the Infernal and
BLAST software packages, which are downloaded and installed with [VADR
installation](install.md#top).

---
## VADR documentation <a name="documentation"></a>

* [VADR installation instructions](install.md#top)
  * [Installation using `vadr-install.sh`](install.md#install)
  * [Setting environment variables](install.md#environment)
  * [Verifying successful installation](install.md#tests)
  * [Further information](install.md#further)
* [`v-build.pl` example usage and command-line options](build.md#top)
  * [`v-build.pl` example usage](build.md#exampleusage)
  * [`v-build.pl` command-line options](build.md#options)
  * [Building a VADR model library](build.md#library)
  * [How the VADR 1.0 model library was constructed](build.md#1.0library)
* [`v-annotate.pl` example usage, command-line options and alert information](annotate.md#top)
  * [`v-annotate.pl` example usage](annotate.md#exampleusage)
  * [`v-annotate.pl` command-line options](annotate.md#options)
  * [Basic Information on `v-annotate.pl` alerts](annotate.md#alerts)
  * [Additional information on `v-annotate.pl` alerts](annotate.md#alerts2)
* [VADR output file formats](formats.md#top)
  * [VADR output files created by all VADR scripts](formats.md#generic)
  * [`v-build.pl` output files](formats.md#build)
  * [`v-annotate.pl` output files](formats.md#annotate)
  * [VADR `coords` coordinate string format](formats.md#coords)
  * [VADR sequence naming conventions](formats.md#seqnames)

---
#### Questions, comments or feature requests? Send a mail to eric.nawrocki@nih.gov.
