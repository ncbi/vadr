#BUILDDIR=/panfs/pan1/infernal/notebook/18_0722_virus_duplicate_features/dnaorg-build-directories
#BUILDDIR=/panfs/pan1/infernal/notebook/18_1002_virus_dnaorg_protein_validation/dnaorg-build-directories
BUILDDIR=/panfs/pan1/dnaorg/virseqannot/dnaorg-build-directories

# same as do-v7-eneida-example-tests.sh 
# except compares output of v7 tests with 'current' expected output, from most recent stable release
# instead of with Eneida's expected output.

# v7 norovirus:
$DNAORGDIR/dnaorg_scripts/dnaorg_test.pl --dirbuild $BUILDDIR/norovirus-builds/NC_029646 -f $DNAORGDIR/dnaorg_scripts/testfiles/noro.current.v7.nocdsgene.testin dt-current-noro-v7-test

# v7 ebolavirus:
$DNAORGDIR/dnaorg_scripts/dnaorg_test.pl --dirbuild $BUILDDIR/ebolavirus-builds/NC_002549 -f $DNAORGDIR/dnaorg_scripts/testfiles/ebola.current.v7.noutr.nocdsgene.testin dt-current-ebola-v7-test
