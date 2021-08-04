#!/bin/bash
perl $VADRSCRIPTSDIR/v-annotate.pl --indefstr 12 --minpvlen 3 --pv_skip --keep --mkey toy50 -mdir $VADRSCRIPTSDIR/documentation/alert-files -f $VADRSCRIPTSDIR/documentation/alert-files/example-indfstrn.fa va-example-indfstrn
