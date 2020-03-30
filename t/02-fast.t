use strict;
use warnings FATAL => 'all';
use Test::More tests => 24;

BEGIN {
    use_ok( 'vadr' )      || print "Bail out!\n";
    use_ok( 'vadr_fast' ) || print "Bail out!\n";
}

###########################################
# parse_blastn_indel_strings() tests
###########################################

my @desc_A = ();
my @in_mdl_coords_A = ();
my @in_seq_coords_A = ();
my @ins_str_A = ();
my @del_str_A = ();
my @exp_mdl_coords_A = ();
my @exp_seq_coords_A = ();

push(@desc_A,           "no indels");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..102:+");
push(@ins_str_A,        "BLASTNULL");
push(@del_str_A,        "BLASTNULL");
push(@exp_mdl_coords_A, "1..100:+");
push(@exp_seq_coords_A, "3..102:+");

push(@desc_A,           "1 insert len 1");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..103:+");
push(@ins_str_A,        "Q12:S10:+1");
push(@del_str_A,        "BLASTNULL");
push(@exp_mdl_coords_A, "1..10:+,11..100:+");
push(@exp_seq_coords_A, "3..12:+,14..103:+");

push(@desc_A,           "1 insert len 3");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..105:+");
push(@ins_str_A,        "Q12:S10:+3");
push(@del_str_A,        "BLASTNULL");
push(@exp_mdl_coords_A, "1..10:+,11..100:+");
push(@exp_seq_coords_A, "3..12:+,16..105:+");

push(@desc_A,           "1 delete len 1");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..101:+");
push(@ins_str_A,        "BLASTNULL");
push(@del_str_A,        "Q12:S10:-1");
push(@exp_mdl_coords_A, "1..10:+,12..100:+");
push(@exp_seq_coords_A, "3..12:+,13..101:+");

push(@desc_A,           "1 delete len 5");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..97:+");
push(@ins_str_A,        "BLASTNULL");
push(@del_str_A,        "Q12:S10:-5");
push(@exp_mdl_coords_A, "1..10:+,16..100:+");
push(@exp_seq_coords_A, "3..12:+,13..97:+");

push(@desc_A,           "2 inserts");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..110:+");
push(@ins_str_A,        "Q12:S10:+3;Q20:S15:+5");
push(@del_str_A,        "BLASTNULL");
push(@exp_mdl_coords_A, "1..10:+,11..15:+,16..100:+");
push(@exp_seq_coords_A, "3..12:+,16..20:+,26..110:+");

push(@desc_A,           "2 deletes");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..94:+");
push(@ins_str_A,        "BLASTNULL");
push(@del_str_A,        "Q12:S10:-5;Q20:S23:-3");
push(@exp_mdl_coords_A, "1..10:+,16..23:+,27..100:+");
push(@exp_seq_coords_A, "3..12:+,13..20:+,21..94:+");

push(@desc_A,           "1 insert len 3, 1 delete len 3");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..102:+");
push(@ins_str_A,        "Q12:S10:+3");
push(@del_str_A,        "Q22:S17:-3");
push(@exp_mdl_coords_A, "1..10:+,11..17:+,21..100:+");
push(@exp_seq_coords_A, "3..12:+,16..22:+,23..102:+");

push(@desc_A,           "1 delete len 3, 1 insert len 3");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..102:+");
push(@ins_str_A,        "Q22:S23:+3");
push(@del_str_A,        "Q12:S10:-3");
push(@exp_mdl_coords_A, "1..10:+,14..23:+,24..100:+");
push(@exp_seq_coords_A, "3..12:+,13..22:+,26..102:+");

push(@desc_A,           "1 insert len 3, 1 delete len 3, 1 insert len 1");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..103:+");
push(@ins_str_A,        "Q12:S10:+3;Q30:S28:+1",);
push(@del_str_A,        "Q22:S17:-3;");
push(@exp_mdl_coords_A, "1..10:+,11..17:+,21..28:+,29..100:+");
push(@exp_seq_coords_A, "3..12:+,16..22:+,23..30:+,32..103:+");

push(@desc_A,           "1 delete len 3, 1 insert len 3, 1 delete len 1");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..101:+");
push(@ins_str_A,        "Q22:S23:+3");
push(@del_str_A,        "Q12:S10:-3;Q30:S28:-1");
push(@exp_mdl_coords_A, "1..10:+,14..23:+,24..28:+,30..100:+");
push(@exp_seq_coords_A, "3..12:+,13..22:+,26..30:+,31..101:+");

my $ntests = scalar(@desc_A);
my $cur_mdl_coords;
my $cur_seq_coords;
my $i;
for($i = 0; $i < $ntests; $i++) { 
  ($cur_mdl_coords, $cur_seq_coords) =
      parse_blastn_indel_strings($in_mdl_coords_A[$i],
                                 $in_seq_coords_A[$i],
                                 $ins_str_A[$i],
                                 $del_str_A[$i],
                                 undef);
  is($cur_mdl_coords, $exp_mdl_coords_A[$i], "parse_blastn_indel_strings() ungapped model coords: $desc_A[$i]");
  is($cur_seq_coords, $exp_seq_coords_A[$i], "parse_blastn_indel_strings() ungapped seq coords: $desc_A[$i]");
}

