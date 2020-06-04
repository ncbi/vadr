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

# default with noro.subseq.fa with short name, should pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa short name, --noseqnamemax not used");
push(@fail_A, "0");

# default with noro.subseq.fa with short name and --noseqnamemax, should pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f --noseqnamemax \$VADRSCRIPTSDIR/testfiles/noro.subseq.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa short name, --noseqnamemax used");
push(@fail_A, "0");

# default with noro.subseq.fa with long name, should fail
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f \$VADRSCRIPTSDIR/testfiles/noro.subseq-longname.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa long name, --noseqnamemax not used");
push(@fail_A, "1");

# default with noro.subseq.fa with long name and --noseqnamemax, should pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -f --noseqnamemax \$VADRSCRIPTSDIR/testfiles/noro.subseq-longname.fa va-test > /dev/null 2>&1");
push(@desc_A, "v-annotate.pl noro.subseq.fa long name, --noseqnamemax used");
push(@fail_A, "0");

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
