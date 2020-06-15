# VADR 1.x release notes 

### VADR 1.1.1 release (June 2020): Minor update

  * in feature table output (only) features that start and/or end with
    one or more ambiguous 'N' nucleotides now have their start and
    stop positions 'trimmed' to the first or final non-N to be
    consisent with how GenBank annotates such features, and adds
    related --noftrtrim and --notrim options to turn this off for some
    or all types of features
  * adds ambgnt5c, ambgnt3c alerts for CDS features that start/end
    with one or more Ns, fatal by default
  * adds ambgnt5f, ambgnt3f alerts for non-CDS features that start/end
    with one or more Ns, non-fatal by default
  * adds ftskipfl alert for rare case that >= 1 per-feature fatal
    alerts exist, but none of those features are included in the
    output feature table
  * sequence names must now be 50 characters or less to be consistent
    with GenBank submission rules, adds --noseqnamemax option to relax
    this restriction

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
branch](https://github.com/nawrockie/vadr/commits/develop).

