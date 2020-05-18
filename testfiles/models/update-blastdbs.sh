for f in NC_001959.vadr.protein.fa NC_029646.vadr.protein.fa NC_039476.vadr.protein.fa NC_039477.vadr.protein.fa NC_001474.vadr.protein.fa NC_001477.vadr.protein.fa; do
    $VADRBLASTDIR/makeblastdb -dbtype prot -in $f
done
for f in NC_001959 NC_039477 noro.3.fa dengue.2.fa entoy100a.fa; do
    $VADRBLASTDIR/makeblastdb -dbtype nucl -in $f
done
cd blastx-NC_001959-multisgm/
$VADRBLASTDIR/makeblastdb -dbtype prot -in NC_001959.vadr.protein.fa
cd ..
