# VADR 1.x release notes 

### VADR 1.1.3 release (Feb 2021): Minor update

  v-annotate.pl changes:
  * adds support for allowing sequences to pass with some
    feature-specific fatal alerts for features with
    'misc_not_failure:"1"' listed in input modelinfo file.
    Originally implemented to allow SARS-CoV-2 sequences with certain
    ORF8 errors to still pass. More information in
    documentation/annotate.md.
  * indf5loc and indf3loc alerts no longer reported for CDS or
    features with identical coordinates to a CDS (e.g. gene) because
    these are subject to more stringent tests involving start/stop
    codons 
  * improved detailed error messages for CDS_HAS_STOP_CODON errors
    when predicted CDS have lengths that are not a multiple of 3
  * outputs separate fasta files of all passing seqs and all failing
    seqs, can be turned off with --out_nofasta
  * fixes bug with -r when number of leading 5' Ns exceeds expected
    number as determined by alignment to best-matching model 
    (github issue #30)
  * with -p, qsub commands now execute a shell script with relevant
    command instead of including command in qsub call, allowing
    arbitrarily long commands and removing need for --longdir option.

  Other changes:
  * updates versions of dependencies installed with vadr to:
    - infernal 1.1.4 (allowing v-build.pl to build large >25Kb models,
      e.g. coronaviridae)
    - hmmer 3.3.2
    - ncbi-blast 2.11.0+
    - Bio-Easel 0.13
    - sequip 0.08

### VADR 1.1.2 release (Nov 2020): Minor update

  * adds v-annotate.pl option --mlist to specify only a subset of
    models in the model info file be used
  * adds v-annotate.pl option --msub to substitute models after the
    classification stage
  * adds v-annotate.pl option --xsub to substitute blast dbs for the
    blastx protein validation stage
  * adds v-annotate.pl options to separately modify minimum length for
    lowsim5s/lowsim5f and lowsim3s/lowsim3f alerts 
  * adds 'benchmark mode' for v-test.pl with -m that produces output
    detailing differences between expected and observed output files
  * enforces that model names not include '(' or ')' (github issue
    #22)
  * fixes bug in install script that prevented building on some Ubuntu
    OS versions (github issue #24)

### VADR 1.1.1 release (July 2020): Minor update

  * in feature table output (only) features that start and/or end with
    one or more ambiguous 'N' nucleotides now have their start and
    stop positions 'trimmed' to the first or final non-N to be
    consisent with how GenBank annotates such features, and adds
    related --noftrtrim and --notrim options to turn this off for some
    or all types of features
  * in .ftr output file, adds two new columns that list number of
    consecutive Ns at beginning and end of each feature (often 0)
  * adds ambgnt5c, ambgnt3c alerts for CDS features that start/end
    with one or more Ns, non-fatal by default
  * adds ambgnt5f, ambgnt3f alerts for non-CDS features that start/end
    with one or more Ns, non-fatal by default
  * adds ftskipfl alert for rare case that >= 1 per-feature fatal
    alerts exist, but none of those features are included in the
    output feature table
  * adds deletinf alert for a sequence that has a feature with a
    complete segment deleted, fatal by default
  * adds deletins alert for a sequence that has a complete feature
    deleted, fatal by default
  * sequence names must now be 50 characters or less to be consistent
    with GenBank submission rules, adds --noseqnamemax option to relax
    this restriction (github issue #12)
  * fixes bug with deleted features causing v-annotate.pl to fail
    (github issue #21) 
  * fixes bug with -s and duplicate blastn hits in models with
    repetitive regions (github issue #13)

---

### VADR 1.1 release (May 2020): Major update
  * adds -s option for accelerating v-annotate.pl, using fixed
    alignment regions derived from blastn, mainly useful for
    SARS-CoV-2 annotation
  * adds -r option for v-annotate.pl for replacing Ns with expected
    nucleotides where possible, motivated by the high fraction of Ns
    in many SARS-CoV-2 sequences
  * adds --hmmer option for profile HMM based protein validation
  * makes fsthicnf alert fatal by default
  * adds ambgnt5s and ambgnt3s alerts, non-fatal by default
  * several bug fixes and other less significant new options
  * additional tests 

---

### VADR 1.0.6 release (April 2020): Bug fix release
* protein_id qualifiers now accession-only unless --forceid used.
* --execname option added to v-{annotate,build,test}.pl

### VADR 1.0.5 release (March 2020): Minor update
* adds protein_id qualifiers to CDS and mat_peptide
  features in output feature tables.

### VADR 1.0.4 release (March 2020): Bug fix release
* fixes installation test (do-install-tests-{local,parallel}.sh),
  which failed in 1.0.3.

### VADR 1.0.3 release (March 2020): Minor update
* Adds frameshift detection capability with associated fsthicnf and
  fstlowcnf alerts, as non-fatal alerts.
* Several bug fixes to previously untested code related to negative
  strand and multisegment features.
* More tests (entoy100a model).

### VADR 1.0.2 release (January 2020): Bug fix release
* Minor update: CDS features that are 5' truncated but not 3'
  truncated are now inspected for a valid stop codon (mutendcd alert
  reported if stop is invalid). Also, adds more informative error
  message for mutstart alerts.

### VADR 1.0.1 release (December 2019): Minor update
* Minor update: adds `--nomisc` option to `v-annotate.pl` for
  preventing conversion of some types of features to `misc_feature` if
  they have fatal alerts. 

### VADR 1.0 release (November 2019): First major release
* The version used in the [VADR manuscript currently in press at BMC
  Bioinformatics.](https://www.biorxiv.org/content/10.1101/852657v2)

---

For more information, see the [git log for the develop
branch](https://github.com/ncbi/vadr/commits/develop).

