#BUILDDIR=/panfs/pan1/infernal/notebook/18_0722_virus_duplicate_features/dnaorg-build-directories
#BUILDDIR=/panfs/pan1/infernal/notebook/18_1002_virus_dnaorg_protein_validation/dnaorg-build-directories
BUILDDIR=/panfs/pan1/dnaorg/virseqannot/dnaorg-build-directories

# noro r1000
$DNAORGDIR/dnaorg_scripts/dnaorg_test.pl --dirbuild $BUILDDIR/norovirus-builds -f $DNAORGDIR/dnaorg_scripts/testfiles/testin.noro.r1000 n1000

# ebola r1000
$DNAORGDIR/dnaorg_scripts/dnaorg_test.pl --dirbuild $BUILDDIR/ebolavirus-builds -f $DNAORGDIR/dnaorg_scripts/testfiles/testin.ebola.r1000 e1000


