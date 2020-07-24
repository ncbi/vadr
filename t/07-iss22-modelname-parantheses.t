use strict;
use warnings FATAL => 'all';
use Test::More tests => 5;

# make sure the VADRINSTALLDIR, VADRSCRIPTSDIR and VADRMODELDIR env variables are set
my $env_ok = exists $ENV{"VADRINSTALLDIR"} ? 1 : 0;
is($env_ok, 1, "VADRINSTALLDIR env variable set");

$env_ok = exists $ENV{"VADRSCRIPTSDIR"} ? 1 : 0;
is($env_ok, 1, "VADRSCRIPTSDIR env variable set");

$env_ok = exists $ENV{"VADRMODELDIR"} ? 1 : 0;
is($env_ok, 1, "VADRMODELDIR env variable set");

# test entoy100a-parantheses.minfo should fail

my @cmd_A     = ();
my @desc_A    = ();
my @fail_A    = ();
my @errfile_A = ();
my @rmdir_A   = ();

# default with noro.subseq.fa with short name, should pass
push(@cmd_A,  "\$VADRSCRIPTSDIR/v-annotate.pl -m \$VADRSCRIPTSDIR/testfiles/models/entoy100a-parantheses.cm -i \$VADRSCRIPTSDIR/testfiles/models/entoy100a-parantheses.minfo -f --skip_pv \$VADRSCRIPTSDIR/testfiles/entoy100a-fs1.fa va-entoy100a-parantheses > /dev/null 2> va-entoy100a-parantheses.err");
push(@errfile_A, "va-entoy100a-parantheses.err");
push(@desc_A, "model name with parantheses");
push(@fail_A, "1");

push (@rmdir_A, "va-entoy100a-parantheses");

my $ncmd = scalar(@cmd_A);
my $retval = undef;
for(my $i = 0; $i < $ncmd; $i++) { 
  $retval = system($cmd_A[$i]);
  if($retval != 0) { $retval = 1; }
  is($retval, $fail_A[$i], sprintf("%s: expected %s", $desc_A[$i], ($fail_A[$i] ? "fail" : "return code 0 (pass)")));
  my $errmsg = `cat $errfile_A[$i]`;
  my $errmatch = 0;
  if($errmsg =~ "not allowed in model names") { $errmatch = 1; }
  is($errmatch, 1, "error message as expected");
}

my $dir;
foreach $dir (@rmdir_A) { 
  system("rm -rf $dir");
}
