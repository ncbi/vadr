use strict;
use warnings FATAL => 'all';
use Test::More tests => 36;

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

#######################
# tests of maximum sequence length in input sequence file

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

###############################################################
#### tests of maximum protein id length in output feature table
# seqs in the noro.subseq-namelen50.fa file have sequence names of length, 48, 49 and 50
@cmd_A = ();
@desc_A = ();
@fail_A = ();
my @longok_A = (); # 0 if it is NOT okay to have long protein_ids > 50 chars, 1 if it is 

# default, protein_ids should have lengths <= 50
push(@cmd_A,    "\$VADRSCRIPTSDIR/v-annotate.pl -f \$VADRSCRIPTSDIR/testfiles/noro.subseq-namelen50.fa va-test > /dev/null 2>&1");
push(@desc_A,   "feature table check: v-annotate.pl default noro.subseq.fa name length 48, 49, 50");
push(@fail_A,   0);
push(@longok_A, 0);

# --forceprotid, protein_ids should have lengths >= 50
push(@cmd_A,    "\$VADRSCRIPTSDIR/v-annotate.pl -f --forceprotid \$VADRSCRIPTSDIR/testfiles/noro.subseq-namelen50.fa va-test > /dev/null 2>&1");
push(@desc_A,   "feature table check: v-annotate.pl --forceprotid noro.subseq.fa name length 48, 49, 50");
push(@fail_A,   0);
push(@longok_A, 1);

# --noseqnamemax, protein_ids should have lengths >= 50
push(@cmd_A,    "\$VADRSCRIPTSDIR/v-annotate.pl -f --noseqnamemax \$VADRSCRIPTSDIR/testfiles/noro.subseq-namelen50.fa va-test > /dev/null 2>&1");
push(@desc_A,   "feature table check: v-annotate.pl --noseqnamemax noro.subseq.fa name length 48, 49, 50");
push(@fail_A,   0);
push(@longok_A, 1);

# --forceprotid and --noseqnamemax, protein_ids should have lengths >= 50
push(@cmd_A,    "\$VADRSCRIPTSDIR/v-annotate.pl -f --forceprotid --noseqnamemax \$VADRSCRIPTSDIR/testfiles/noro.subseq-namelen50.fa va-test > /dev/null 2>&1");
push(@desc_A,   "feature table check: v-annotate.pl --forceprotid --noseqnamemax noro.subseq.fa name length 48, 49, 50");
push(@fail_A,   0);
push(@longok_A, 1);

$ncmd = scalar(@cmd_A);
$retval = undef;
my $ftbl_file = "va-test/va-test.vadr.pass.tbl";
for(my $i = 0; $i < $ncmd; $i++) { 
  $retval = system($cmd_A[$i]);
  if($retval != 0) { $retval = 1; }
  is($retval, $fail_A[$i], sprintf("%s: expected %s", $desc_A[$i], ($fail_A[$i] ? "fail" : "return code 0 (pass)")));
  # parse the pass.tbl to get the protein_id values
  open(IN, $ftbl_file) || die "ERROR unable to open $ftbl_file";
  my $j = 1;
  while(my $line = <IN>) { 
    if($line =~ /^\s+protein_id\s+(\S+)/) {
      my $prot_id_value = $1;
      my $len = length($prot_id_value);
      my $more_than_50 = ($len > 50) ? 1 : 0;
      my $exactly_50 = ($len == 50) ? 1 : 0;
      if($longok_A[$i]) { 
        if($prot_id_value =~ /48/) { 
          # 48 char name should go to 50
          is($exactly_50, 1, "$desc_A[$i] protein_id $prot_id_value, length $len == 50, as expected");
        }
        else { 
          # 49 and 50 char names should go to more than 50
          is($more_than_50, 1, "$desc_A[$i] protein_id $prot_id_value, length $len > 50, as expected");
        }
      }
      else { # longid not okay, all ids should be 50 chars
        is($exactly_50, 1, "$desc_A[$i] protein_id $prot_id_value, length $len == 50, as expected");
      }
    }
  }
  close(IN);
}

# add test with >= 10 protein_ids, or just manually test

my $dir;
foreach $dir (@rmdir_A) { 
#  system("rm -rf $dir");
}
