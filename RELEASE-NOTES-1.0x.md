# VADR 1.0x release notes 

### VADR 1.0.5 release (March 2020)
* Minor update: adds protein_id qualifiers to CDS and mat_peptide
  features in output feature tables.

### VADR 1.0.4 release (March 2020)
* Bug fix release: fixes installation test
  (do-install-tests-{local,parallel}.sh), which failed in 1.0.3.

### VADR 1.0.3 release (March 2020)
* Minor update: Adds frameshift detection capability with associated
  fsthicnf and fstlowcnf alerts, as non-fatal alerts.
* Several bug fixes to previously untested code related to
  negative strand and multisegment features. 
* More tests (entoy100a model).

### VADR 1.0.2 release (January 2020)
* Minor update: CDS features that are 5' truncated but not 3'
  truncated are now inspected for a valid stop codon (mutendcd alert
  reported if stop is invalid). Also, adds more informative error
  message for mutstart alerts.

### VADR 1.0.1 release (December 2019)
* Minor update: adds `--nomisc` option to `v-annotate.pl` for
  preventing conversion of some types of features to `misc_feature` if
  they have fatal alerts. 

### VADR 1.0 release (November 2019)
* First major release. The version used in the [VADR manuscript submitted
  to BMC Bioinformatics.](https://www.biorxiv.org/content/10.1101/852657v1)

For more information, see the [git log for the develop branch](https://github.com/nawrockie/vadr/commits/develop).
