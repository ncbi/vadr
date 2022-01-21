rm -rf va-iss58
mkdir va-iss58
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/esl-reformat fasta /home/nawrocke/git/nawrockie/vadr/testfiles/github-issues/iss58-cdsstopp/iss58.fa > va-iss58/va-iss58.vadr.in.fa
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/esl-seqstat --dna -a va-iss58/va-iss58.vadr.in.fa > va-iss58/va-iss58.vadr.seqstat
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/ncbi-blast/bin/blastn -num_threads 1 -query va-iss58/va-iss58.vadr.blastn.fa -db /home/nawrocke/git/nawrockie/vadr/testfiles/models/noro.3.fa -out va-iss58/va-iss58.vadr.rpn.cls.blastn.out -word_size 7 -reward 1 -penalty -2 -xdrop_gap_final 110 -gapopen 2 -gapextend 1
/home/nawrocke/git/nawrockie/vadr/parse_blast.pl --program n --input va-iss58/va-iss58.vadr.rpn.cls.blastn.out --splus > va-iss58/va-iss58.vadr.rpn.cls.blastn.summary.txt
grep -v ^# va-iss58/va-iss58.vadr.rpn.cls.tblout | sed 's/  */ /g' | sort -k 1,1 -k 3,3rn > va-iss58/va-iss58.vadr.rpn.cls.tblout.sort
cat va-iss58/va-iss58.vadr.rpn.cdt.NC_029646.tblout | grep -v ^# | sed 's/  */ /g' | sort -k 1,1 -k 15,15rn -k 16,16g > va-iss58/va-iss58.vadr.rpn.cdt.tblout.sort
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/cmsearch  -T -10 --cpu 1 --trmF3 --noali --hmmonly --tblout va-iss58/va-iss58.vadr.std.cls.s0.tblout /home/nawrocke/git/nawrockie/vadr/testfiles/models/noro.3.cm va-iss58/va-iss58.vadr.rpn.fa > va-iss58/va-iss58.vadr.std.cls.s0.stdout
cat va-iss58/va-iss58.vadr.std.cls.s0.stdout > va-iss58/va-iss58.vadr.std.cls.stdout
cat va-iss58/va-iss58.vadr.std.cls.s0.tblout > va-iss58/va-iss58.vadr.std.cls.tblout
grep -v ^# va-iss58/va-iss58.vadr.std.cls.tblout | sed 's/  */ /g' | sort -k 1,1 -k 3,3rn > va-iss58/va-iss58.vadr.std.cls.tblout.sort
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/cmfetch /home/nawrocke/git/nawrockie/vadr/testfiles/models/noro.3.cm NC_029646 | /panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/cmsearch  -T -10 --cpu 1 --hmmonly  --noali  --tblout va-iss58/va-iss58.vadr.std.cdt.s0.tblout - va-iss58/va-iss58.vadr.NC_029646.fa > va-iss58/va-iss58.vadr.std.cdt.s0.stdout
cat va-iss58/va-iss58.vadr.std.cdt.s0.stdout > va-iss58/va-iss58.vadr.std.cdt.NC_029646.stdout
cat va-iss58/va-iss58.vadr.std.cdt.s0.tblout > va-iss58/va-iss58.vadr.std.cdt.NC_029646.tblout
cat va-iss58/va-iss58.vadr.std.cdt.NC_029646.tblout | grep -v ^# | sed 's/  */ /g' | sort -k 1,1 -k 15,15rn -k 16,16g > va-iss58/va-iss58.vadr.std.cdt.tblout.sort
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/cmfetch /home/nawrocke/git/nawrockie/vadr/testfiles/models/noro.3.cm NC_029646 | /panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/cmalign  --dnaout --verbose --cpu 1 --ifile va-iss58/va-iss58.vadr.NC_029646.align.r1.s0.ifile -o va-iss58/va-iss58.vadr.NC_029646.align.r1.s0.stk --tau 0.001 --mxsize 4000.00 --sub --notrunc -g --fixedtau - va-iss58/va-iss58.vadr.NC_029646.a.fa > va-iss58/va-iss58.vadr.NC_029646.align.r1.s0.stdout 2>&1
cat va-iss58/va-iss58.vadr.NC_029646.align.r1.s0.stdout > va-iss58/va-iss58.vadr.NC_029646.align.stdout
cat va-iss58/va-iss58.vadr.NC_029646.align.r1.s0.ifile > va-iss58/va-iss58.vadr.NC_029646.align.ifile
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/esl-alimerge --list --outformat stockholm --dna va-iss58/va-iss58.vadr.NC_029646.align.stk.list > va-iss58/va-iss58.vadr.NC_029646.rfrna.align.rpstk
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/esl-reformat --informat stockholm -d pfam va-iss58/va-iss58.vadr.NC_029646.align.rpstk > va-iss58/va-iss58.vadr.NC_029646.align.rpstk.stockholm.pfam
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/esl-alimanip --informat stockholm --dna --outformat stockholm --num-rf va-iss58/va-iss58.vadr.NC_029646.align.rpstk.stockholm.new.pfam > va-iss58/va-iss58.vadr.NC_029646.align.stk
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/esl-alimerge --list --outformat afa --dna va-iss58/va-iss58.vadr.NC_029646.align.stk.list > va-iss58/va-iss58.vadr.NC_029646.align.rpafa
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/esl-reformat --informat afa -d pfam va-iss58/va-iss58.vadr.NC_029646.align.rpafa > va-iss58/va-iss58.vadr.NC_029646.align.rpafa.afa.pfam
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/infernal/binaries/esl-alimanip --informat stockholm --dna --outformat afa va-iss58/va-iss58.vadr.NC_029646.align.rpafa.afa.new.pfam > va-iss58/va-iss58.vadr.NC_029646.align.afa
cat va-iss58/va-iss58.vadr.NC_029646.CDS.1.pv.fa >> va-iss58/va-iss58.vadr.NC_029646.pv.blastx.fa
cat va-iss58/va-iss58.vadr.NC_029646.CDS.2.pv.fa >> va-iss58/va-iss58.vadr.NC_029646.pv.blastx.fa
/panfs/pan1/infernal/notebook/22_0106_vadr_1p4p1_release/post-release-install-20220106/ncbi-blast/bin/blastx -num_threads 1 -num_alignments 20 -query va-iss58/va-iss58.vadr.NC_029646.pv.blastx.fa -db /home/nawrocke/git/nawrockie/vadr/testfiles/models/NC_029646.vadr.protein.fa -seg no -out va-iss58/va-iss58.vadr.NC_029646.blastx.out
/home/nawrocke/git/nawrockie/vadr/parse_blast.pl --program x --input va-iss58/va-iss58.vadr.NC_029646.blastx.out > va-iss58/va-iss58.vadr.NC_029646.blastx.summary.txt
# Thu Jan 20 21:25:24 EST 2022
# Linux cbbdev13 3.10.0-1160.49.1.el7.x86_64 #1 SMP Tue Nov 30 15:51:32 UTC 2021 x86_64 x86_64 x86_64 GNU/Linux
[ok]
