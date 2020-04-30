rm noro.3.cm*
rm noro.3.minfo
for a in NC_029646 NC_039476 NC_039477; do 
    cp $VADRMODELDIR/$a.vadr.protein.fa* .
    $VADRINFERNALDIR/cmfetch $VADRMODELDIR/vadr.cm $a >> noro.3.cm
    grep $a $VADRMODELDIR/vadr.minfo >> noro.3.minfo
done
$VADRINFERNALDIR/cmpress noro.3.cm
$VADRINFERNALDIR/cmemit -c noro.3.cm > tmp.fa
$VADREASELDIR/esl-reformat -d fasta tmp.fa > tmp2.fa
cat tmp2.fa | sed 's/-hmmconsensus//' > noro.3.fa
$VADRBLASTDIR/makeblastdb -in noro.3.fa -dbtype nucl
rm tmp.fa
rm tmp2.fa
