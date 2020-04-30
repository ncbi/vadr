rm dengue.2.cm*
rm dengue.2.minfo
for a in NC_001474 NC_001477; do
    cp $VADRMODELDIR/$a.vadr.protein.fa* .
    $VADRINFERNALDIR/cmfetch $VADRMODELDIR/vadr.cm $a >> dengue.2.cm
    grep $a $VADRMODELDIR/vadr.minfo >> dengue.2.minfo
done
$VADRINFERNALDIR/cmpress dengue.2.cm
$VADRINFERNALDIR/cmemit -c dengue.2.cm > tmp.fa
$VADREASELDIR/esl-reformat -d fasta tmp.fa > tmp2.fa
cat tmp2.fa | sed 's/-cmconsensus//' > dengue.2.fa
$VADRBLASTDIR/makeblastdb -in dengue.2.fa -dbtype nucl
rm tmp.fa
rm tmp2.fa
