rm -rf va-nn
mkdir va-nn
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/esl-reformat fasta test-nn.fa > va-nn/va-nn.vadr.in.fa
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/esl-seqstat --dna -a va-nn/va-nn.vadr.in.fa > va-nn/va-nn.vadr.seqstat
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmsearch  -T -10 --cpu 1 --trmF3 --noali --hmmonly --tblout va-nn/va-nn.vadr.std.cls.s0.tblout /net/intdev/oblast01/infernal/git/nawrockie/vadr/documentation/annotate-files/toy-nn.cm va-nn/va-nn.vadr.in.fa > va-nn/va-nn.vadr.std.cls.s0.stdout
cat va-nn/va-nn.vadr.std.cls.s0.stdout > va-nn/va-nn.vadr.std.cls.stdout
rm  va-nn/va-nn.vadr.std.cls.s0.stdout
cat va-nn/va-nn.vadr.std.cls.s0.tblout > va-nn/va-nn.vadr.std.cls.tblout
rm  va-nn/va-nn.vadr.std.cls.s0.tblout
grep -v ^# va-nn/va-nn.vadr.std.cls.tblout | sed 's/  */ /g' | sort -k 1,1 -k 3,3rn > va-nn/va-nn.vadr.std.cls.tblout.sort
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmfetch /net/intdev/oblast01/infernal/git/nawrockie/vadr/documentation/annotate-files/toy-nn.cm toy-nn | /net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmsearch  -T -10 --cpu 1 --hmmonly  --noali  --tblout va-nn/va-nn.vadr.std.cdt.s0.tblout - va-nn/va-nn.vadr.toy-nn.fa > va-nn/va-nn.vadr.std.cdt.s0.stdout
cat va-nn/va-nn.vadr.std.cdt.s0.stdout > va-nn/va-nn.vadr.std.cdt.toy-nn.stdout
rm  va-nn/va-nn.vadr.std.cdt.s0.stdout
cat va-nn/va-nn.vadr.std.cdt.s0.tblout > va-nn/va-nn.vadr.std.cdt.toy-nn.tblout
rm  va-nn/va-nn.vadr.std.cdt.s0.tblout
cat va-nn/va-nn.vadr.std.cdt.toy-nn.tblout | grep -v ^# | sed 's/  */ /g' | sort -k 1,1 -k 15,15rn -k 16,16g > va-nn/va-nn.vadr.std.cdt.tblout.sort
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmfetch /net/intdev/oblast01/infernal/git/nawrockie/vadr/documentation/annotate-files/toy-nn.cm toy-nn | /net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmalign  --dnaout --verbose --cpu 1 --ifile va-nn/va-nn.vadr.toy-nn.align.r1.s0.ifile -o va-nn/va-nn.vadr.toy-nn.align.r1.s0.stk --tau 0.001 --mxsize 4000.00 --sub --notrunc -g --fixedtau --flanktoins 0.1 --flankselfins 0.8 - va-nn/va-nn.vadr.toy-nn.a.fa > va-nn/va-nn.vadr.toy-nn.align.r1.s0.stdout 2>&1
cat va-nn/va-nn.vadr.toy-nn.align.r1.s0.stdout > va-nn/va-nn.vadr.toy-nn.align.stdout
rm  va-nn/va-nn.vadr.toy-nn.align.r1.s0.stdout
cat va-nn/va-nn.vadr.toy-nn.align.r1.s0.ifile > va-nn/va-nn.vadr.toy-nn.align.ifile
rm  va-nn/va-nn.vadr.toy-nn.align.r1.s0.ifile
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/esl-alimerge --list --outformat stockholm  va-nn/va-nn.vadr.toy-nn.align.stk.list > va-nn/va-nn.vadr.toy-nn.rfrna.align.stk
rm  va-nn/va-nn.vadr.in.fa va-nn/va-nn.vadr.in.fa.ssi va-nn/va-nn.vadr.std.cls.tblout va-nn/va-nn.vadr.std.cls.stdout va-nn/va-nn.vadr.std.cls.tblout.sort va-nn/va-nn.vadr.toy-nn.fa va-nn/va-nn.vadr.std.cdt.toy-nn.tblout va-nn/va-nn.vadr.std.cdt.toy-nn.stdout va-nn/va-nn.vadr.std.cdt.tblout.sort va-nn/va-nn.vadr.toy-nn.a.fa va-nn/va-nn.vadr.toy-nn.align.stdout va-nn/va-nn.vadr.toy-nn.align.ifile va-nn/va-nn.vadr.toy-nn.align.r1.s0.stk va-nn/va-nn.vadr.toy-nn.rfrna.align.stk va-nn/va-nn.vadr.toy-nn.align.stk.list
# Fri Jan  9 14:23:34 EST 2026
# Linux cbbdev13 4.18.0-553.85.1.el8_10.x86_64 #1 SMP Mon Nov 24 09:05:24 EST 2025 x86_64 x86_64 x86_64 GNU/Linux
[ok]
