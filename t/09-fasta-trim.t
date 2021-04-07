use strict;
use warnings FATAL => 'all';
use Test::More tests => 18;

BEGIN {
    use_ok( 'vadr' )      || print "Bail out!\n";
    use_ok( 'vadr_seed' ) || print "Bail out!\n";
}

# make sure the VADRSCRIPTSDIR env variable is set
my $env_ok = exists $ENV{"VADRSCRIPTSDIR"} ? 1 : 0;
is($env_ok, 1, "VADRSCRIPTSDIR env variable set");

# make sure we have the required files
my @reqd_file_A = ();
my $trim_script_path         = $ENV{"VADRSCRIPTSDIR"} . "/miniscripts/fasta-trim-terminal-ambigs.pl";
my $trim_fasta_path          = $ENV{"VADRSCRIPTSDIR"} . "/t/data/trim1.fa";
my $exp_df_fasta_path        = $ENV{"VADRSCRIPTSDIR"} . "/t/data/trim1_trimmed.fa";
my $exp_minlen_fasta_path    = $ENV{"VADRSCRIPTSDIR"} . "/t/data/trim1_trimmed.min10.fa";
my $exp_maxlen_fasta_path    = $ENV{"VADRSCRIPTSDIR"} . "/t/data/trim1_trimmed.max10.fa";
my $exp_minmaxlen_fasta_path = $ENV{"VADRSCRIPTSDIR"} . "/t/data/trim1_trimmed.min10.max17.fa";
my $exp_sfx_fasta_path       = $ENV{"VADRSCRIPTSDIR"} . "/t/data/trim1_trimmed.sfx.fa";

push(@reqd_file_A, $trim_script_path);
push(@reqd_file_A, $trim_fasta_path);
push(@reqd_file_A, $exp_df_fasta_path);
push(@reqd_file_A, $exp_minlen_fasta_path);
push(@reqd_file_A, $exp_maxlen_fasta_path);
push(@reqd_file_A, $exp_minmaxlen_fasta_path);
push(@reqd_file_A, $exp_sfx_fasta_path);

foreach my $reqd_file (@reqd_file_A) {
  my $file_nonempty = (-s $reqd_file) ? 1 : 0;
  is($file_nonempty, 1, "Required file $reqd_file exists and is nonempty");
}

###################################
# 
# run script
my @to_remove_A = ();
my $tmp_out  = "out.fa.tmp";
my $diff_out = "out.fa.diff";
push(@to_remove_A, $tmp_out);

my $cmd;
my $exp_val;

############
# default 
$cmd = "perl $trim_script_path $trim_fasta_path > $tmp_out";
system($cmd);
if($? != 0) { die "ERROR command $cmd failed"; }

$cmd = "diff $tmp_out $exp_df_fasta_path > $diff_out";
system($cmd);

$exp_val = `cat $diff_out`;
is($exp_val, "", "trimming script reproduced expected output (default)");
#############

############
# --minlen 10 
$cmd = "perl $trim_script_path --minlen 10 $trim_fasta_path > $tmp_out";
system($cmd);
if($? != 0) { die "ERROR command $cmd failed"; }

$cmd = "diff $tmp_out $exp_minlen_fasta_path > $diff_out";
system($cmd);

$exp_val = `cat $diff_out`;
is($exp_val, "", "trimming script reproduced expected output (minlen)");
#############

############
# --maxlen 10 
$cmd = "perl $trim_script_path --maxlen 10 $trim_fasta_path > $tmp_out";
system($cmd);
if($? != 0) { die "ERROR command $cmd failed"; }

$cmd = "diff $tmp_out $exp_maxlen_fasta_path > $diff_out";
system($cmd);

$exp_val = `cat $diff_out`;
is($exp_val, "", "trimming script reproduced expected output (maxlen)");
#############

############
# --minlen 10 --maxlen 17 
$cmd = "perl $trim_script_path --minlen 10 --maxlen 17 $trim_fasta_path > $tmp_out";
system($cmd);
if($? != 0) { die "ERROR command $cmd failed"; }

$cmd = "diff $tmp_out $exp_minmaxlen_fasta_path > $diff_out";
system($cmd);

$exp_val = `cat $diff_out`;
is($exp_val, "", "trimming script reproduced expected output (minmaxlen)");
#############

############
# --sfx .test
$cmd = "perl $trim_script_path --sfx .test $trim_fasta_path > $tmp_out";
system($cmd);
if($? != 0) { die "ERROR command $cmd failed"; }

$cmd = "diff $tmp_out $exp_sfx_fasta_path > $diff_out";
system($cmd);

$exp_val = `cat $diff_out`;
is($exp_val, "", "trimming script reproduced expected output (sfx)");
#############

# expected failures:
############
# --min -1 
$cmd = "perl $trim_script_path --min -1 $trim_fasta_path 2>/dev/null";
system($cmd);
my $did_fail = ($? == 0) ? 0 : 1;
is($did_fail, 1, "trimming script failed with --min -1 as expected");

#############
# --min 10 --max 9 
$cmd = "perl $trim_script_path --min 10 --max 9 $trim_fasta_path 2>/dev/null";
system($cmd);
$did_fail = ($? == 0) ? 0 : 1;
is($did_fail, 1, "trimming script failed with --min 10 --max 9 as expected");

#############

#############
# --min 10 --strict
$cmd = "perl $trim_script_path --min 10 strict $trim_fasta_path 2>/dev/null";
system($cmd);
$did_fail = ($? == 0) ? 0 : 1;
is($did_fail, 1, "trimming script failed with --min 10 --strict as expected");

#############

foreach my $tmp_file (@to_remove_A) {
  unlink $tmp_file;
}

