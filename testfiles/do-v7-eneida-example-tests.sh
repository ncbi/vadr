#BUILDDIR=/panfs/pan1/infernal/notebook/18_0722_virus_duplicate_features/dnaorg-build-directories
BUILDDIR=/panfs/pan1/dnaorg/virseqannot/dnaorg-build-directories

# same as do-v7-current-example-tests.sh 
# except compares output of v7 tests with Eneida's expected output instead of 
# 'current' expected output from most recent stable release

# v7 norovirus:
# https://confluence.ncbi.nlm.nih.gov/pages/viewpage.action?spaceKey=VG&title=Annotation+Test+Sets+and+5-Column+Feature+Table
# original file with expected sequence table:
# On https://jira.ncbi.nlm.nih.gov/browse/VIV-528 
# Noro_II_feature_tables_v7.tbl [April 18, 2018]
# Then to get noro.v7.nocdsgene.nocodonstart.sqtable, I removed all codon_start qualifiers and all gene qualifiers in CDS features, which 
# dnaorg pipeline does not attempt to annotate.
# See /panfs/pan1/infernal/notebook/18_0507_dnaorg_virus_example_feature_table_update/00LOG.txt for details
$DNAORGDIR/dnaorg_scripts/dnaorg_test.pl --skipmsg --dirbuild $BUILDDIR/norovirus-builds/NC_029646 -f $DNAORGDIR/dnaorg_scripts/testfiles/noro.eneida.v7.nocdsgene.nocodonstart.testin dt-eneida-noro-v7-test

# v7 ebolavirus:
# https://confluence.ncbi.nlm.nih.gov/pages/viewpage.action?spaceKey=VG&title=Annotation+Test+Sets+and+5-Column+Feature+Table
# original file with expected sequence table:
# On https://jira.ncbi.nlm.nih.gov/browse/VIV-528 
# Zaire_ebola_feature_tables_v7.tbl [April 18, 2018]
# Then to get ebola.v7.noutr.nocdsgene.nocodonstart.sqtable, I removed all utr annotation, all codon_start qualifiers and 
# all gene qualifiers in CDS features, which 
# dnaorg pipeline does not attempt to annotate.
# See /panfs/pan1/infernal/notebook/18_0507_dnaorg_virus_example_feature_table_update/00LOG.txt for details
$DNAORGDIR/dnaorg_scripts/dnaorg_test.pl --skipmsg --dirbuild $BUILDDIR/ebolavirus-builds/NC_002549 -f $DNAORGDIR/dnaorg_scripts/testfiles/ebola.eneida.v7.noutr.nocdsgene.nocodonstart.testin dt-eneida-ebola-v7-test

