mkdir entoy100a-minfo
perl $VADRSCRIPTSDIR/v-annotate.pl -f --alicheck --minpvlen 3 --pv_skip -m $VADRSCRIPTSDIR/testfiles/models/entoy100a.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a.qual.minfo $VADRSCRIPTSDIR/testfiles/entoy100a.minfo.fa va-entoy100a-minfo1 > va-entoy100a-minfo1.out
diff -U 0 va-entoy100a-minfo1/va-entoy100a-minfo1.vadr.pass.tbl /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo1/va-entoy100a-minfo1.vadr.pass.tbl > entoy100a-minfo/entoy100a-minfo.vadr.1.1.diff
diff -U 0 va-entoy100a-minfo1/va-entoy100a-minfo1.vadr.fail.tbl /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo1/va-entoy100a-minfo1.vadr.fail.tbl > entoy100a-minfo/entoy100a-minfo.vadr.1.2.diff
cp va-entoy100a-minfo1/va-entoy100a-minfo1.vadr.fail.tbl entoy100a-minfo/entoy100a-minfo.vadr.1.2.diff.out
cp /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo1/va-entoy100a-minfo1.vadr.fail.tbl entoy100a-minfo/entoy100a-minfo.vadr.1.2.diff.exp
perl $VADRSCRIPTSDIR/v-annotate.pl --forcegene -f --alicheck --minpvlen 3 --pv_skip -m $VADRSCRIPTSDIR/testfiles/models/entoy100a.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a.qual.minfo $VADRSCRIPTSDIR/testfiles/entoy100a.minfo.fa va-entoy100a-minfo2 > va-entoy100a-minfo2.out
diff -U 0 va-entoy100a-minfo2/va-entoy100a-minfo2.vadr.pass.tbl /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo2/va-entoy100a-minfo2.vadr.pass.tbl > entoy100a-minfo/entoy100a-minfo.vadr.2.1.diff
diff -U 0 va-entoy100a-minfo2/va-entoy100a-minfo2.vadr.fail.tbl /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo2/va-entoy100a-minfo2.vadr.fail.tbl > entoy100a-minfo/entoy100a-minfo.vadr.2.2.diff
cp va-entoy100a-minfo2/va-entoy100a-minfo2.vadr.fail.tbl entoy100a-minfo/entoy100a-minfo.vadr.2.2.diff.out
cp /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo2/va-entoy100a-minfo2.vadr.fail.tbl entoy100a-minfo/entoy100a-minfo.vadr.2.2.diff.exp
perl $VADRSCRIPTSDIR/v-annotate.pl -f --forcequal qual1 --alicheck --minpvlen 3 --pv_skip -m $VADRSCRIPTSDIR/testfiles/models/entoy100a.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a.qual.minfo $VADRSCRIPTSDIR/testfiles/entoy100a.minfo.fa va-entoy100a-minfo3 > va-entoy100a-minfo3.out
diff -U 0 va-entoy100a-minfo3/va-entoy100a-minfo3.vadr.pass.tbl /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo3/va-entoy100a-minfo3.vadr.pass.tbl > entoy100a-minfo/entoy100a-minfo.vadr.3.1.diff
diff -U 0 va-entoy100a-minfo3/va-entoy100a-minfo3.vadr.fail.tbl /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo3/va-entoy100a-minfo3.vadr.fail.tbl > entoy100a-minfo/entoy100a-minfo.vadr.3.2.diff
cp va-entoy100a-minfo3/va-entoy100a-minfo3.vadr.fail.tbl entoy100a-minfo/entoy100a-minfo.vadr.3.2.diff.out
cp /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo3/va-entoy100a-minfo3.vadr.fail.tbl entoy100a-minfo/entoy100a-minfo.vadr.3.2.diff.exp
perl $VADRSCRIPTSDIR/v-annotate.pl -f --forcequal qual2,rfam_id,qual3 --alicheck --minpvlen 3 --pv_skip -m $VADRSCRIPTSDIR/testfiles/models/entoy100a.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a.qual.minfo $VADRSCRIPTSDIR/testfiles/entoy100a.minfo.fa va-entoy100a-minfo4 > va-entoy100a-minfo4.out
diff -U 0 va-entoy100a-minfo4/va-entoy100a-minfo4.vadr.pass.tbl /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo4/va-entoy100a-minfo4.vadr.pass.tbl > entoy100a-minfo/entoy100a-minfo.vadr.4.1.diff
diff -U 0 va-entoy100a-minfo4/va-entoy100a-minfo4.vadr.fail.tbl /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo4/va-entoy100a-minfo4.vadr.fail.tbl > entoy100a-minfo/entoy100a-minfo.vadr.4.2.diff
cp va-entoy100a-minfo4/va-entoy100a-minfo4.vadr.fail.tbl entoy100a-minfo/entoy100a-minfo.vadr.4.2.diff.out
cp /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/expected-files/va-entoy100a-minfo4/va-entoy100a-minfo4.vadr.fail.tbl entoy100a-minfo/entoy100a-minfo.vadr.4.2.diff.exp
# Mon Jan  5 14:47:22 EST 2026
# Linux cbbdev13 4.18.0-553.85.1.el8_10.x86_64 #1 SMP Mon Nov 24 09:05:24 EST 2025 x86_64 x86_64 x86_64 GNU/Linux
[FAIL]
