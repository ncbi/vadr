use strict;
use warnings FATAL => 'all';
use Test::More tests => 4275;

BEGIN {
    use_ok( 'vadr' )      || print "Bail out!\n";
    use_ok( 'vadr_seed' ) || print "Bail out!\n";
}

# make sure the VADRSCRIPTSDIR env variable is set
my $env_ok = exists $ENV{"VADRSCRIPTSDIR"} ? 1 : 0;
is($env_ok, 1, "VADRSCRIPTSDIR env variable set");

# make sure we have the required files
my @reqd_file_A = ();
my $parse_blast_path = $ENV{"VADRSCRIPTSDIR"} . "/parse_blast.pl";
push(@reqd_file_A, $parse_blast_path);

foreach my $reqd_file (@reqd_file_A) {
  my $file_nonempty = (-s $reqd_file) ? 1 : 0;
  is($file_nonempty, 1, "Required file $reqd_file exists and is nonempty");
}

###################################
# blastx testing
@reqd_file_A = ();
my $blastx_in_path   = $ENV{"VADRSCRIPTSDIR"} . "/t/data/blastx.1.out";
my $blastx_sum_path  = $ENV{"VADRSCRIPTSDIR"} . "/t/data/blastx.1.summary";
push(@reqd_file_A, $blastx_in_path);
push(@reqd_file_A, $blastx_sum_path);
foreach my $reqd_file (@reqd_file_A) {
  my $file_nonempty = (-s $reqd_file) ? 1 : 0;
  is($file_nonempty, 1, "Required file $reqd_file exists and is nonempty");
}

# run parse_blastx.pl
my @to_remove_A = ();
my $tmp_blastx_out = "blastx1.tmp";
push(@to_remove_A, $tmp_blastx_out);
my $cmd = "perl $parse_blast_path --in $blastx_in_path --program x > $tmp_blastx_out";
system($cmd);
if($? != 0) { die "ERROR command $cmd failed"; }

# grep all line types from the output file and compare with expected output file
my @line_types_A = ("BITSCORE", "DEL", "END_MATCH", "EVALUE", "FRAME",
                    "GAPS", "HACC", "HDEF", "HLEN", "HSP", "IDENT",
                    "INS", "MATCH", "MAXDE", "MAXIN", "QACC",
                    "QDEF", "QLEN", "QRANGE", "RAWSCORE", "SLEN",
                    "SRANGE", "STOP"); 

my $exp_val;
my $cur_val;
my $cur_nlines;
my $exp_nlines;
foreach my $line_type (@line_types_A) {
  $exp_val = `grep ^" . $line_type . " $blastx_sum_path`;
  $cur_val = `grep ^" . $line_type . " $tmp_blastx_out`;
  $cur_nlines = $cur_val =~ tr/\n//;
  $exp_nlines = $exp_val =~ tr/\n//;
  is($cur_nlines, $exp_nlines, "blastn test 1 $line_type num lines match expected");
  my @cur_A = split(/\n/, $cur_val);
  my @exp_A = split(/\n/, $exp_val);
  for(my $i = 0; $i < scalar(@exp_A); $i++) { 
    is($cur_A[$i], $exp_A[$i], "blastn test 1 $line_type lines $i match expected");
  }
}
  
####################################
# blastn testing

@reqd_file_A = ();
my $blastn_in_path   = $ENV{"VADRSCRIPTSDIR"} . "/t/data/blastn.1.out";
my $blastn_sum_path  = $ENV{"VADRSCRIPTSDIR"} . "/t/data/blastn.1.summary";
push(@reqd_file_A, $blastn_in_path);
push(@reqd_file_A, $blastn_sum_path);
foreach my $reqd_file (@reqd_file_A) {
  my $file_nonempty = (-s $reqd_file) ? 1 : 0;
  is($file_nonempty, 1, "Required file $reqd_file exists and is nonempty");
}

# run parse_blastn.pl
@to_remove_A = ();
my $tmp_blastn_out = "blastn1.tmp";
push(@to_remove_A, $tmp_blastn_out);
$cmd = "perl $parse_blast_path --in $blastn_in_path --program n --splus > $tmp_blastn_out";
system($cmd);
if($? != 0) { die "ERROR command $cmd failed"; }

# grep all line types from the output file and compare with expected output file
@line_types_A = ("BITSCORE", "DEL", "END_MATCH", "EVALUE", "GAPS",
                    "HACC", "HDEF", "HLEN", "HSP", "IDENT",
                    "INS", "MATCH", "MAXDE", "MAXIN", "QACC",
                    "QDEF", "QLEN", "QRANGE", "QSTRAND", "RAWSCORE",
                    "SLEN", "SRANGE", "SSTRAND", "STOP");

foreach my $line_type (@line_types_A) {
  $exp_val = `grep ^$line_type $blastn_sum_path`;
  $cur_val = `grep ^$line_type $tmp_blastn_out`;
  #print("EXPECTED\n$exp_val\n");
  #print("CURRENT\n$cur_val\n");
  $cur_nlines = $cur_val =~ tr/\n//;
  $exp_nlines = $exp_val =~ tr/\n//;
  is($cur_nlines, $exp_nlines, "blastn test 1 $line_type num lines match expected");
  my @cur_A = split(/\n/, $cur_val);
  my @exp_A = split(/\n/, $exp_val);
  for(my $i = 0; $i < scalar(@exp_A); $i++) { 
    is($cur_A[$i], $exp_A[$i], "blastn test 1 $line_type lines $i match expected");
  }
}

foreach my $tmp_file (@to_remove_A) {
  unlink $tmp_file;
}

