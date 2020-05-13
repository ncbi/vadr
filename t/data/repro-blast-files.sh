# blastn
$VADRBLASTDIR/makeblastdb -in noro.20.fa -dbtype nucl
$VADRBLASTDIR/blastn -word_size 7 -out blastn.1.out -db noro.20.fa -query noro.1.fa -num_threads 1
perl $VADRSCRIPTSDIR/parse_blast.pl --splus --program n --input blastn.1.out > blastn.1.summary

# blastx
$VADRBLASTDIR/makeblastdb -in NC_039477.vadr.protein.fa -dbtype prot
$VADRBLASTDIR/blastx -num_threads 1 -query noro.3.pv.blastx.fa -db NC_039477.vadr.protein.fa -seg no -out blastx.1.out
perl $VADRSCRIPTSDIR/parse_blast.pl --program x --input blastx.1.out > blastx.1.summary
