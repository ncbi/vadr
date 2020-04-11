use strict;
use warnings FATAL => 'all';
use Test::More tests => 23;

# make sure the VADRINSTALLDIR, VADRSCRIPTSDIR and VADRMODELDIR env variables are set
my $env_ok = exists $ENV{"VADRINSTALLDIR"} ? 1 : 0;
is($env_ok, 1, "VADRINSTALLDIR env variable set");

$env_ok = exists $ENV{"VADRSCRIPTSDIR"} ? 1 : 0;
is($env_ok, 1, "VADRSCRIPTSDIR env variable set");

$env_ok = exists $ENV{"VADRMODELDIR"} ? 1 : 0;
is($env_ok, 1, "VADRMODELDIR env variable set");

# set VADRMODELDIR to default
# uncomment this if you've set VADRMODELDIR to something other than default, and you want to test that
$ENV{'VADRMODELDIR'} = '$VADRINSTALLDIR/vadr-models';

# could check for required files here, but there are many, so I just 
# rely on them being there, if they're not the tests below that 
# require them will fail.

my @cmd_A   = ();
my @desc_A  = ();
my @fail_A  = ();
my @rmdir_A = ();

# -h: pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -h > /dev/null");
push(@desc_A, "v-annotate.pl -h");
push(@fail_A, "0");

# default: pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa default models");
push(@fail_A, "0");

# --mdir without vadr.cm etc.: fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f --mdir \$VADRSCRIPTSDIR/testfiles/models \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa --mdir");
push(@fail_A, "1");

# change VADRMODELDIR env variable then use --mkey: pass
push(@cmd_A,  "export VADRMODELDIR=\$VADRSCRIPTSDIR/testfiles/models; \$VADRSCRIPTSDIR/v-annotate.pl -f --mkey NC_001959 \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa change VADRMODELDIR --mkey");
push(@fail_A, "0");

# change VADRMODELDIR env variable to bogus dir, then use --mkey: fail
push(@cmd_A,  "export VADRMODELDIR=\$VADRSCRIPTSDIR/bogus; \$VADRSCRIPTSDIR/v-annotate.pl -f --mkey NC_001959 \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa bogus VADRMODELDIR --mkey");
push(@fail_A, "1");

# --mdir --mkey: pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f --mdir \$VADRSCRIPTSDIR/testfiles/models --mkey NC_001959 \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa --mdir --mkey");
push(@fail_A, "0");

# nonexistent --mdir : fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f --mdir \$VADRSCRIPTSDIR/bogus --mkey NC_001959 \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa nonexistent --mdir");
push(@fail_A, "1");

# existent but unequipped --mdir : fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f --mdir \$VADRSCRIPTSDIR --mkey NC_001959 \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa existent but unequipped --mdir");
push(@fail_A, "1");

# --mdir bad --mkey: fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f --mdir \$VADRSCRIPTSDIR/testfiles/models --mkey bogus \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa --mdir bad --mkey");
push(@fail_A, "1");

# --mdir --mkey -s: pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f -s --mdir \$VADRSCRIPTSDIR/testfiles/models --mkey NC_001959 \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa -s --mdir --mkey");
push(@fail_A, "0");

# --mdir --mkey --addhmmer --skipblastx: pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f --mdir \$VADRSCRIPTSDIR/testfiles/models --mkey NC_001959 --addhmmer --skipblastx \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa -addhmmer --skipblast --mdir --mkey");
push(@fail_A, "0");

# --mdir --mkey -s --addhmmer --skipblastx: pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f -s --mdir \$VADRSCRIPTSDIR/testfiles/models --mkey NC_001959 --addhmmer --skipblastx \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa -s -addhmmer --skipblast --mdir --mkey");
push(@fail_A, "0");

# valid -m -i -x: pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f -m \$VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i \$VADRSCRIPTSDIR/testfiles/models/NC_001959.minfo -x \$VADRSCRIPTSDIR/testfiles/models \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa -m -i -x");
push(@fail_A, "0");

# valid -m -i bad -x: fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f -m \$VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i \$VADRSCRIPTSDIR/testfiles/models/NC_001959.minfo -x \$VADRSCRIPTSDIR/bogus \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa -m -i bad -x");
push(@fail_A, "1");

# valid -m -x bad -i: fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f -m \$VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i \$VADRSCRIPTSDIR/testfiles/models/bogus.minfo -x \$VADRSCRIPTSDIR/testfiles/models \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa -m -x bad -i");
push(@fail_A, "1");

# valid -i -x bad -m: fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f -m \$VADRSCRIPTSDIR/testfiles/models/bogus.cm -i \$VADRSCRIPTSDIR/testfiles/models/NC_001959.minfo -x \$VADRSCRIPTSDIR/testfiles/models \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa -i -x bad -m");
push(@fail_A, "1");

# -s and valid -m -i -n -x: pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -s -f -m \$VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i \$VADRSCRIPTSDIR/testfiles/models/NC_001959.minfo -x \$VADRSCRIPTSDIR/testfiles/models -n \$VADRSCRIPTSDIR/testfiles/models/NC_001959.fa \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa -s -m -i -n -x");
push(@fail_A, "0");

# -s and valid -m -i -x bad -n: fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -s -f -m \$VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i \$VADRSCRIPTSDIR/testfiles/models/NC_001959.minfo -x \$VADRSCRIPTSDIR/testfiles/models -n \$VADRSCRIPTSDIR/testfiles/models/bogus.fa \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa -s -m -i -x bad -n");
push(@fail_A, "1");

# --addhmmer --skipblastx and valid -m -i -a: pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl --addhmmer --skipblastx -f -m \$VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i \$VADRSCRIPTSDIR/testfiles/models/NC_001959.minfo -a \$VADRSCRIPTSDIR/testfiles/models/NC_001959.hmm \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa --addhmmer --skipblastx -m -i -a");
push(@fail_A, "0");

# --addhmmer --skipblastx and valid -m -i bad -a: fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl --addhmmer --skipblastx -f -m \$VADRSCRIPTSDIR/testfiles/models/NC_001959.cm -i \$VADRSCRIPTSDIR/testfiles/models/NC_001959.minfo -a \$VADRSCRIPTSDIR/testfiles/models/bogus.hmm \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa --addhmmer --skipblastx -m -i bad -a");
push(@fail_A, "1");

push (@rmdir_A, "va-test");

my $ncmd = scalar(@cmd_A);
my $retval = undef;
for(my $i = 0; $i < $ncmd; $i++) { 
  $retval = system($cmd_A[$i]);
  if($retval != 0) { $retval = 1; }
  is($retval, $fail_A[$i], sprintf("%s: expected %s", $desc_A[$i], ($fail_A[$i] ? "fail" : "return code 0 (pass)")));
}

my $dir;
foreach $dir (@rmdir_A) { 
  system("rm -rf $dir");
}
