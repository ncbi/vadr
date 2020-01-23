# VADR 1.0x release notes 

### VADR 1.0.2 release (Jan 2020)
* Minor update: CDS features that are 5' truncated but not 3'
  truncated are now inspected for a valid stop codon (mutendcd alert
  reported if stop is invalid). Also, adds more informative error
  message for mutstart alerts.

### VADR 1.0.1 release (Dec 2019)
* Minor update: adds `--nomisc` option to `v-annotate.pl` for
  preventing conversion of some types of features to `misc_feature` if
  they have fatal alerts. 

### VADR 1.0 release (Nov 2019)
* First major release. The version used in the [VADR manuscript submitted
  to BMC Bioinformatics.](https://www.biorxiv.org/content/10.1101/852657v1)

For more information, see the [git log for the develop branch](https://github.com/nawrockie/vadr/commits/develop).
