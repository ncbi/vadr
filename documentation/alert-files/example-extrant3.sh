#!/bin/bash
perl $VADRSCRIPTSDIR/v-annotate.pl --minpvlen 3 --lowcov 0.8 --pv_skip --keep -m $VADRSCRIPTSDIR/documentation/alert-files/toy50.vadrprior.cm -i $VADRSCRIPTSDIR/documentation/alert-files/toy50.minfo -f $VADRSCRIPTSDIR/documentation/alert-files/example-extrant3.fa va-example-extrant3
