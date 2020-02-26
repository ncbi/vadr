for a in fs1 fs2 fs3; do 
    perl revcomp-stk.pl entoy100a-$a.stk 100 ENTOY100A ENTOY100A-REV | esl-alimanip --num-rf - > tmp.stk
    perl ~/src/jiffy-infernal-hmmer-scripts/ali-pfam-lowercase-rf-gap-columns.pl tmp.stk > entoy100a-rev-$a.stk
    esl-reformat --informat stockholm fasta entoy100a-rev-$a.stk > entoy100a-rev-$a.fa
done

