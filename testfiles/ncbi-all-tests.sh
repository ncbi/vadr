### tests from .testin file dengue.r5.hmmer.testin
#annotate-dengue-5-hmmer
vadr_noro -- -f --hmmer --alicheck --mdir $VADRSCRIPTSDIR/testfiles/models --mkey dengue.2 $VADRSCRIPTSDIR/testfiles/dengue.r5.fa va-dengue-hmmer.r5 > va-dengue-hmmer.r5.out
### tests from .testin file dengue.r5.local.testin
#annotate-dengue-5-local
vadr_dengue -- -f --alicheck $VADRSCRIPTSDIR/testfiles/dengue.r5.fa va-dengue.r5 > va-dengue.r5.out
### tests from .testin file dengue.r5.parallel.testin
#annotate-dengue-5-parallel
vadr_dengue -- -p -f --alicheck $VADRSCRIPTSDIR/testfiles/dengue.r5.fa va-dengue.r5 > va-dengue.r5.out
### tests from .testin file dengue.r5.rpn.testin
#annotate-dengue-5-rpn
vadr_noro -- -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey dengue.2 --alicheck -r $VADRSCRIPTSDIR/testfiles/dengue.r5.rpn.fa va-dengue-rpn.r5 > va-dengue-rpn.r5.out
#annotate-dengue-5-rpn-prof
vadr_noro -- -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey dengue.2 --alicheck -r --r_prof $VADRSCRIPTSDIR/testfiles/dengue.r5.rpn.fa va-dengue-rpn-prof.r5 > va-dengue-rpn-prof.r5.out
### tests from .testin file dengue.r5.seed.testin
#annotate-dengue-5-seed
vadr_noro -- -s -f --alicheck --mdir $VADRSCRIPTSDIR/testfiles/models --mkey dengue.2 $VADRSCRIPTSDIR/testfiles/dengue.r5.fa va-dengue-seed.r5 > va-dengue-seed.r5.out
### tests from .testin file entoy100a-fs.testin
#annotate-entoy100a-fs1
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 --fstminnt 1 --minpvlen 3 --keep -v --skip_pv -m $VADRSCRIPTSDIR/testfiles/models/entoy100a.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a.minfo -f $VADRSCRIPTSDIR/testfiles/entoy100a-fs1.fa va-entoy100a-fs1 > va-entoy100a-fs1.out
#annotate-entoy100a-fs2
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 --fstminnt 1 --minpvlen 3 --keep -v --skip_pv -m $VADRSCRIPTSDIR/testfiles/models/entoy100a.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a.minfo -f $VADRSCRIPTSDIR/testfiles/entoy100a-fs2.fa va-entoy100a-fs2 > va-entoy100a-fs2.out
#annotate-entoy100a-fs3
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 --fstminnt 1 --minpvlen 3 --keep -v --skip_pv -m $VADRSCRIPTSDIR/testfiles/models/entoy100a.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a.minfo -f $VADRSCRIPTSDIR/testfiles/entoy100a-fs3.fa va-entoy100a-fs3 > va-entoy100a-fs3.out
### tests from .testin file entoy100a-nindel.testin
#annotate-entoy100a-nindel
vadr_noro -- -f --nmaxdel 1 --nmaxins 2 --alicheck --minpvlen 3 --skip_pv --mdir $VADRSCRIPTSDIR/testfiles/models --mkey entoy100a $VADRSCRIPTSDIR/testfiles/entoy100a-fs2.fa va-entoy100a-fs2-nindel > va-entoy100a-fs2-nindel.out
#annotate-entoy100a-nindel-nmaxins_exc
vadr_noro -- -f --nmaxdel 1 --nmaxins 2 --alicheck --minpvlen 3 --skip_pv -m $VADRSCRIPTSDIR/testfiles/models/entoy100a.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a.nmaxins_exc.minfo $VADRSCRIPTSDIR/testfiles/entoy100a-fs2.fa va-entoy100a-fs2-nindel-nmaxins_exc > va-entoy100a-fs2-nindel-nmaxins_exc.out
### tests from .testin file entoy100a-rev-fs.testin
#annotate-entoy100a-rev-fs1
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 --fstminnt 1 --minpvlen 3 --keep -v --skip_pv -m $VADRSCRIPTSDIR/testfiles/models/entoy100a-rev.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a-rev.minfo -f $VADRSCRIPTSDIR/testfiles/entoy100a-rev-fs1.fa va-entoy100a-rev-fs1 > va-entoy100a-rev-fs1.out
#annotate-entoy100a-rev-fs2
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 --fstminnt 1 --minpvlen 3 --keep -v --skip_pv -m $VADRSCRIPTSDIR/testfiles/models/entoy100a-rev.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a-rev.minfo -f $VADRSCRIPTSDIR/testfiles/entoy100a-rev-fs2.fa va-entoy100a-rev-fs2 > va-entoy100a-rev-fs2.out
#annotate-entoy100a-rev-fs3
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 --fstminnt 1 --minpvlen 3 --keep -v --skip_pv -m $VADRSCRIPTSDIR/testfiles/models/entoy100a-rev.cm -i $VADRSCRIPTSDIR/testfiles/models/entoy100a-rev.minfo -f $VADRSCRIPTSDIR/testfiles/entoy100a-rev-fs3.fa va-entoy100a-rev-fs3 > va-entoy100a-rev-fs3.out
### tests from .testin file noro.fs.multisgm.testin
#annotate-noro-fs-multisgm.1
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 -f -m $VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i $VADRSCRIPTSDIR/testfiles/models/NC_001959.multisgm.minfo -x $VADRSCRIPTSDIR/testfiles/models/blastx-NC_001959-multisgm $VADRSCRIPTSDIR/testfiles/noro.fs.1.fa va-noro.fs.multisgm.1 > va-noro.fs.multisgm.1.out
#annotate-noro-fs-multisgm.2
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 -f -m $VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i $VADRSCRIPTSDIR/testfiles/models/NC_001959.multisgm.minfo -x $VADRSCRIPTSDIR/testfiles/models/blastx-NC_001959-multisgm $VADRSCRIPTSDIR/testfiles/noro.fs.2.fa va-noro.fs.multisgm.2 > va-noro.fs.multisgm.2.out
### tests from .testin file noro.fs.testin
#annotate-noro-fs.1
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 -f -m $VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i $VADRSCRIPTSDIR/testfiles/models/NC_001959.minfo -x $VADRSCRIPTSDIR/testfiles/models $VADRSCRIPTSDIR/testfiles/noro.fs.1.fa va-noro.fs.1 > va-noro.fs.1.out
#annotate-noro-fs.2
vadr_noro -- --alicheck --alt_fail fsthicnf --fsthighthr 0.300001 -f -m $VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i $VADRSCRIPTSDIR/testfiles/models/NC_001959.minfo -x $VADRSCRIPTSDIR/testfiles/models $VADRSCRIPTSDIR/testfiles/noro.fs.2.fa va-noro.fs.2 > va-noro.fs.2.out
### tests from .testin file noro.r10.hmmer.testin
#annotate-noro-10-hmmer
vadr_noro -- -f --hmmer --alicheck --mdir $VADRSCRIPTSDIR/testfiles/models --mkey noro.3 $VADRSCRIPTSDIR/testfiles/noro.r10.fa va-noro-hmmer.r10 > va-noro-hmmer.r10.out
### tests from .testin file noro.r10.local.testin
#annotate-noro-10-local
vadr_noro -- -f --alicheck $VADRSCRIPTSDIR/testfiles/noro.r10.fa va-noro.r10 > va-noro.r10.out
### tests from .testin file noro.r10.parallel.testin
#annotate-noro-10-local
vadr_noro -- --alicheck -p -f $VADRSCRIPTSDIR/testfiles/noro.r10.fa va-noro.r10 > va-noro.r10.out
### tests from .testin file noro.r10.rpn.testin
#annotate-noro-10-rpn
vadr_noro -- --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey noro.3 -r $VADRSCRIPTSDIR/testfiles/noro.r10.rpn.fa va-noro-rpn.r10 > va-noro-rpn.r10.out
#annotate-noro-10-rpn-prof
vadr_noro -- --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey noro.3 -r --r_prof $VADRSCRIPTSDIR/testfiles/noro.r10.rpn.fa va-noro-rpn-prof.r10 > va-noro-rpn-prof.r10.out
### tests from .testin file noro.r10.seed.testin
#annotate-noro-10-seed
vadr_noro -- -s --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey noro.3 $VADRSCRIPTSDIR/testfiles/noro.r10.fa va-noro-seed.r10 > va-noro-seed.r10.out
### tests from .testin file noro.r3.outaln.testin
#annotate-noro-1-outaln-os
vadr_noro -- --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_stk $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-os.r3 > va-noro-outaln-os.r3.out
#annotate-noro-1-outaln-oa
vadr_noro -- --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_afa $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-oa.r3 > va-noro-outaln-oa.r3.out
#annotate-noro-1-outaln-osa
vadr_noro -- --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_stk --out_afa $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-osa.r3 > va-noro-outaln-osa.r3.out
#annotate-noro-1-outaln-r-os
vadr_noro -- -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_stk $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-os.r3 > va-noro-outaln-r-os.r3.out
#annotate-noro-1-outaln-r-oa
vadr_noro -- -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_afa $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-oa.r3 > va-noro-outaln-r-oa.r3.out
#annotate-noro-1-outaln-r-osa
vadr_noro -- -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_stk --out_afa $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-osa.r3 > va-noro-outaln-r-osa.r3.out
#annotate-noro-1-outaln-r-osrs
vadr_noro -- -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_stk --out_rpstk $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-osrs.r3 > va-noro-outaln-r-osrs.r3.out
#annotate-noro-1-outaln-r-oara
vadr_noro -- -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_afa --out_rpafa $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-oara.r3 > va-noro-outaln-r-oara.r3.out
#annotate-noro-1-outaln-r-osrsoara
vadr_noro -- -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_stk --out_afa --out_rpstk --out_rpafa $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-osrsoara.r3 > va-noro-outaln-r-osrsoara.r3.out
#annotate-noro-1-outaln-r-rs
vadr_noro -- -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_rpstk $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-rs.r3 > va-noro-outaln-r-rs.r3.out
#annotate-noro-1-outaln-r-ra
vadr_noro -- -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_rpafa $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-ra.r3 > va-noro-outaln-r-ra.r3.out
#annotate-noro-1-outaln-r-rsra
vadr_noro -- -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 --out_rpstk --out_rpafa $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-rsra.r3 > va-noro-outaln-r-rsra.r3.out
#annotate-noro-1-outaln-keep
vadr_noro -- --keep -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-keep.r3 > va-noro-outaln-keep.r3.out
#annotate-noro-1-outaln-r-keep
vadr_noro -- --keep -r --alicheck -f --mdir $VADRSCRIPTSDIR/testfiles/models --mkey NC_039477 $VADRSCRIPTSDIR/testfiles/noro.r3.rpn.fa va-noro-outaln-r-keep.r3 > va-noro-outaln-r-keep.r3.out
