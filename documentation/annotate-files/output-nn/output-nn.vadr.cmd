rm -rf output-nn
mkdir output-nn
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/esl-reformat fasta test-seqs.fa > output-nn/output-nn.vadr.in.fa
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/esl-seqstat --dna -a output-nn/output-nn.vadr.in.fa > output-nn/output-nn.vadr.seqstat
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmsearch  -T -10 --cpu 1 --trmF3 --noali --hmmonly --tblout output-nn/output-nn.vadr.std.cls.s0.tblout /net/intdev/oblast01/infernal/git/nawrockie/vadr/documentation/annotate-files/toy-ev.cm output-nn/output-nn.vadr.in.fa > output-nn/output-nn.vadr.std.cls.s0.stdout
cat output-nn/output-nn.vadr.std.cls.s0.stdout > output-nn/output-nn.vadr.std.cls.stdout
rm  output-nn/output-nn.vadr.std.cls.s0.stdout
cat output-nn/output-nn.vadr.std.cls.s0.tblout > output-nn/output-nn.vadr.std.cls.tblout
rm  output-nn/output-nn.vadr.std.cls.s0.tblout
grep -v ^# output-nn/output-nn.vadr.std.cls.tblout | sed 's/  */ /g' | sort -k 1,1 -k 3,3rn > output-nn/output-nn.vadr.std.cls.tblout.sort
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmfetch /net/intdev/oblast01/infernal/git/nawrockie/vadr/documentation/annotate-files/toy-ev.cm toy-ev | /net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmsearch  -T -10 --cpu 1 --hmmonly  --noali  --tblout output-nn/output-nn.vadr.std.cdt.s0.tblout - output-nn/output-nn.vadr.toy-ev.fa > output-nn/output-nn.vadr.std.cdt.s0.stdout
cat output-nn/output-nn.vadr.std.cdt.s0.stdout > output-nn/output-nn.vadr.std.cdt.toy-ev.stdout
rm  output-nn/output-nn.vadr.std.cdt.s0.stdout
cat output-nn/output-nn.vadr.std.cdt.s0.tblout > output-nn/output-nn.vadr.std.cdt.toy-ev.tblout
rm  output-nn/output-nn.vadr.std.cdt.s0.tblout
cat output-nn/output-nn.vadr.std.cdt.toy-ev.tblout | grep -v ^# | sed 's/  */ /g' | sort -k 1,1 -k 15,15rn -k 16,16g > output-nn/output-nn.vadr.std.cdt.tblout.sort
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmfetch /net/intdev/oblast01/infernal/git/nawrockie/vadr/documentation/annotate-files/toy-ev.cm toy-ev | /net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/cmalign  --dnaout --verbose --cpu 1 --ifile output-nn/output-nn.vadr.toy-ev.align.r1.s0.ifile -o output-nn/output-nn.vadr.toy-ev.align.r1.s0.stk --tau 0.001 --mxsize 4000.00 --sub --notrunc -g --fixedtau --flanktoins 0.1 --flankselfins 0.8 - output-nn/output-nn.vadr.toy-ev.a.fa > output-nn/output-nn.vadr.toy-ev.align.r1.s0.stdout 2>&1
cat output-nn/output-nn.vadr.toy-ev.align.r1.s0.stdout > output-nn/output-nn.vadr.toy-ev.align.stdout
rm  output-nn/output-nn.vadr.toy-ev.align.r1.s0.stdout
cat output-nn/output-nn.vadr.toy-ev.align.r1.s0.ifile > output-nn/output-nn.vadr.toy-ev.align.ifile
rm  output-nn/output-nn.vadr.toy-ev.align.r1.s0.ifile
/net/intdev/oblast01/infernal/notebook/25_0902_vadr_1p7_release/test-install1/infernal/binaries/esl-alimerge --list --outformat stockholm  output-nn/output-nn.vadr.toy-ev.align.stk.list > output-nn/output-nn.vadr.toy-ev.rfrna.align.stk
rm  output-nn/output-nn.vadr.in.fa output-nn/output-nn.vadr.in.fa.ssi output-nn/output-nn.vadr.std.cls.tblout output-nn/output-nn.vadr.std.cls.stdout output-nn/output-nn.vadr.std.cls.tblout.sort output-nn/output-nn.vadr.toy-ev.fa output-nn/output-nn.vadr.std.cdt.toy-ev.tblout output-nn/output-nn.vadr.std.cdt.toy-ev.stdout output-nn/output-nn.vadr.std.cdt.tblout.sort output-nn/output-nn.vadr.toy-ev.a.fa output-nn/output-nn.vadr.toy-ev.align.stdout output-nn/output-nn.vadr.toy-ev.align.ifile output-nn/output-nn.vadr.toy-ev.align.r1.s0.stk output-nn/output-nn.vadr.toy-ev.rfrna.align.stk output-nn/output-nn.vadr.toy-ev.align.stk.list
# Fri Jan  9 12:59:11 EST 2026
# Linux cbbdev13 4.18.0-553.85.1.el8_10.x86_64 #1 SMP Mon Nov 24 09:05:24 EST 2025 x86_64 x86_64 x86_64 GNU/Linux
[ok]
