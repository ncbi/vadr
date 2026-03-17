rm -rf profbuild
mkdir profbuild
perl $VADRSCRIPTSDIR/v-build.pl -f --profile --stk $VADRSCRIPTSDIR/testfiles/profile-build/test1_singlepos.stk --minfoin $VADRSCRIPTSDIR/testfiles/profile-build/test1_singlepos.minfo test1 test1_singlepos_out > test1_singlepos.out
diff -U 0 test1_singlepos_out/test1_singlepos_out.vadr.cds.fa /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/profile-build/test1_singlepos.cds.fa > profbuild/profbuild.vadr.1.1.diff
diff -U 0 test1_singlepos_out/test1_singlepos_out.vadr.protein.fa /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/profile-build/test1_singlepos.protein.fa > profbuild/profbuild.vadr.1.2.diff
diff -U 0 test1_singlepos_out/test1_singlepos_out.vadr.minfo /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/profile-build/test1_singlepos.minfo > profbuild/profbuild.vadr.1.3.diff
rm -rf test1_singlepos_out
perl $VADRSCRIPTSDIR/v-build.pl -f --profile --stk $VADRSCRIPTSDIR/testfiles/profile-build/test2_multipos.stk --minfoin $VADRSCRIPTSDIR/testfiles/profile-build/test2_multipos.minfo test2 test2_multipos_out > test2_multipos.out
diff -U 0 test2_multipos_out/test2_multipos_out.vadr.cds.fa /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/profile-build/test2_multipos.cds.fa > profbuild/profbuild.vadr.2.1.diff
diff -U 0 test2_multipos_out/test2_multipos_out.vadr.protein.fa /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/profile-build/test2_multipos.protein.fa > profbuild/profbuild.vadr.2.2.diff
rm -rf test2_multipos_out
perl $VADRSCRIPTSDIR/v-build.pl -f --profile --stk $VADRSCRIPTSDIR/testfiles/profile-build/test3_singleneg.stk --minfoin $VADRSCRIPTSDIR/testfiles/profile-build/test3_singleneg.minfo test3 test3_singleneg_out > test3_singleneg.out
diff -U 0 test3_singleneg_out/test3_singleneg_out.vadr.cds.fa /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/profile-build/test3_singleneg.cds.fa > profbuild/profbuild.vadr.3.1.diff
diff -U 0 test3_singleneg_out/test3_singleneg_out.vadr.protein.fa /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/profile-build/test3_singleneg.protein.fa > profbuild/profbuild.vadr.3.2.diff
rm -rf test3_singleneg_out
perl $VADRSCRIPTSDIR/v-build.pl -f --profile --stk $VADRSCRIPTSDIR/testfiles/profile-build/test4_multineg.stk --minfoin $VADRSCRIPTSDIR/testfiles/profile-build/test4_multineg.minfo test4 test4_multineg_out > test4_multineg.out
diff -U 0 test4_multineg_out/test4_multineg_out.vadr.cds.fa /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/profile-build/test4_multineg.cds.fa > profbuild/profbuild.vadr.4.1.diff
diff -U 0 test4_multineg_out/test4_multineg_out.vadr.protein.fa /net/intdev/oblast01/infernal/git/nawrockie/vadr/testfiles/profile-build/test4_multineg.protein.fa > profbuild/profbuild.vadr.4.2.diff
rm -rf test4_multineg_out
# Mon Jan  5 16:07:35 EST 2026
# Linux cbbdev13 4.18.0-553.85.1.el8_10.x86_64 #1 SMP Mon Nov 24 09:05:24 EST 2025 x86_64 x86_64 x86_64 GNU/Linux
[ok]
