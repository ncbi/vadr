#BUILDDIR=/panfs/pan1/infernal/notebook/18_0722_virus_duplicate_features/dnaorg-build-directories
#BUILDDIR=/panfs/pan1/infernal/notebook/18_1002_virus_dnaorg_protein_validation/dnaorg-build-directories
BUILDDIR=/panfs/pan1/dnaorg/virseqannot/dnaorg-build-directories

# noro r100
$DNAORGDIR/dnaorg_scripts/dnaorg_test.pl --dirbuild $BUILDDIR/norovirus-builds -f $DNAORGDIR/dnaorg_scripts/testfiles/testin.noro.r100 n100

# ebola r100
$DNAORGDIR/dnaorg_scripts/dnaorg_test.pl --dirbuild $BUILDDIR/ebolavirus-builds -f $DNAORGDIR/dnaorg_scripts/testfiles/testin.ebola.r100 e100


