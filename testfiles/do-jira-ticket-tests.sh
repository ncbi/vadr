#BUILDDIR=/panfs/pan1/infernal/notebook/18_0722_virus_duplicate_features/dnaorg-build-directories
#BUILDDIR=/panfs/pan1/infernal/notebook/18_1002_virus_dnaorg_protein_validation/dnaorg-build-directories
BUILDDIR=/panfs/pan1/dnaorg/virseqannot/dnaorg-build-directories

# https://jira.ncbi.nlm.nih.gov/browse/VIV-664
# https://jira.ncbi.nlm.nih.gov/browse/VIV-665
# https://jira.ncbi.nlm.nih.gov/browse/VIV-697
# https://jira.ncbi.nlm.nih.gov/browse/VIV-698
# https://jira.ncbi.nlm.nih.gov/browse/VIV-789
# https://jira.ncbi.nlm.nih.gov/browse/VIV-792
# https://jira.ncbi.nlm.nih.gov/browse/VIV-793
# https://jira.ncbi.nlm.nih.gov/browse/VIV-818
$DNAORGDIR/dnaorg_scripts/dnaorg_test.pl --dirbuild $BUILDDIR/norovirus-builds -f $DNAORGDIR/dnaorg_scripts/testfiles/jira-viv.testin dt-jira-viv

