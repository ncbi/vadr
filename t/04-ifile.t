use strict;
use warnings FATAL => 'all';
use Test::More tests => 13;

BEGIN {
    use_ok( 'vadr' )      || print "Bail out!\n";
}

# make sure the VADRSCRIPTSDIR env variable is set
my $env_ok = exists $ENV{"VADRSCRIPTSDIR"} ? 1 : 0;
is($env_ok, 1, "VADRSCRIPTSDIR env variable set");

# make sure we have the required files
my @reqd_file_A = ();
my $ifile_path = $ENV{"VADRSCRIPTSDIR"} . "/t/data/noro.ifile";
push(@reqd_file_A, $ifile_path);

foreach my $reqd_file (@reqd_file_A) {
  my $file_nonempty = (-s $reqd_file) ? 1 : 0;
  is($file_nonempty, 1, "Required file $reqd_file exists and is nonempty");
}

my %seq_inserts_HH = ();
my %mdl_name_H     = ();
my @seq_name_A     = ();
my %seq_len_H      = ();
my %mdl_len_H      = ();
my $mdl_name;
my $seq_name;
my %FH_H = ();

# parse the ifile
vdr_CmalignParseInsertFile($ifile_path, \%seq_inserts_HH, \%mdl_name_H, \@seq_name_A, \%seq_len_H, \%mdl_len_H, \%FH_H);
# now write it out again
my @mdl_order_A = ();
my %mdl_order_H = ();
my %seq_by_model_HA = ();
foreach $seq_name (@seq_name_A) {
  my $mdl_name = $mdl_name_H{$seq_name};
  if(! defined $mdl_order_H{$mdl_name}) { 
    push(@mdl_order_A, $mdl_name);
    $mdl_order_H{$mdl_name} = 1;
  }
  push(@{$seq_by_model_HA{$mdl_name}}, $seq_name);
}
my $ifile_out = "noro.out.ifile";
if(-s $ifile_out) { unlink $ifile_out; }
foreach $mdl_name (@mdl_order_A) {
  vdr_CmalignWriteInsertFile($ifile_out,
                             1, # $do_append
                             $mdl_name, $mdl_len_H{$mdl_name},
                             \@{$seq_by_model_HA{$mdl_name}}, \%seq_len_H,
                             \%seq_inserts_HH, \%FH_H);
}

# check the two files for equivalence

# first grep the non model lines, these should be equivalent, check number of lines then strings
my ($cur_val, $exp_val);
$exp_val = `grep -v ^\# $ifile_path | grep -v ^\/\/ | grep -v ^NC | wc -l`;
$cur_val = `grep -v ^\# $ifile_out  | grep -v ^\/\/ | grep -v ^NC | wc -l`;
is($cur_val, $exp_val, "ifile: parse/write test grep sequence line count match expected");

$exp_val = `grep -v ^\# $ifile_path | grep -v ^\/\/ | grep -v ^NC`;
$cur_val = `grep -v ^\# $ifile_out  | grep -v ^\/\/ | grep -v ^NC`;
is($cur_val, $exp_val, "ifile: parse/write test grep sequence lines match expected");

######
# now grep the model lines, and remove dups, these should be equivalent too
$exp_val = `grep -v ^\# $ifile_path | grep -v ^\/\/ | grep ^NC | sort | uniq `;
$cur_val = `grep -v ^\# $ifile_out  | grep -v ^\/\/ | grep ^NC | sort | uniq `;
is($cur_val, $exp_val, "ifile: parse/write test grep uniq model lines match expected");

# now check actual values for expected values
my $nseq = scalar(@seq_name_A);
is($nseq, 167, "ifile: number of sequences correct");

#NC_039477 7567
#JQ911595.1 7511 3 7513  2560 2553 3  2583 2579 3
$mdl_name = "NC_039477";
$seq_name = "JQ911595.1";
is($mdl_name_H{$seq_name}, "NC_039477", "ifile: anecdotal mdl name is expected");
is($seq_len_H{$seq_name},  7511,      "ifile: anecdotal mdl name is expected");
is($mdl_len_H{$mdl_name},  7567,      "ifile: anecdotal mdl len is expected");
is($seq_inserts_HH{$seq_name}{"spos"}, 3, "ifile: anecdotal spos is expected");
is($seq_inserts_HH{$seq_name}{"epos"}, 7513, "ifile: anecdotal spos is expected");
is($seq_inserts_HH{$seq_name}{"ins"}, "2560:2553:3;2583:2579:3;", "ifile: anecdotal ins is expected");

foreach my $tmp_file ($ifile_out) { 
  unlink $tmp_file;
}

