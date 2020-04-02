use strict;
use warnings FATAL => 'all';
use Test::More tests => 61;

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
push(@ins_str_A,        "Q12:S10+1");
push(@del_str_A,        "BLASTNULL");
push(@exp_mdl_coords_A, "1..10:+,11..100:+");
push(@exp_seq_coords_A, "3..12:+,14..103:+");

push(@desc_A,           "1 insert len 3");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..105:+");
push(@ins_str_A,        "Q12:S10+3");
push(@del_str_A,        "BLASTNULL");
push(@exp_mdl_coords_A, "1..10:+,11..100:+");
push(@exp_seq_coords_A, "3..12:+,16..105:+");

push(@desc_A,           "1 delete len 1");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..101:+");
push(@ins_str_A,        "BLASTNULL");
push(@del_str_A,        "Q12:S10-1");
push(@exp_mdl_coords_A, "1..10:+,12..100:+");
push(@exp_seq_coords_A, "3..12:+,13..101:+");

push(@desc_A,           "1 delete len 5");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..97:+");
push(@ins_str_A,        "BLASTNULL");
push(@del_str_A,        "Q12:S10-5");
push(@exp_mdl_coords_A, "1..10:+,16..100:+");
push(@exp_seq_coords_A, "3..12:+,13..97:+");

push(@desc_A,           "2 inserts");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..110:+");
push(@ins_str_A,        "Q12:S10+3;Q20:S15+5");
push(@del_str_A,        "BLASTNULL");
push(@exp_mdl_coords_A, "1..10:+,11..15:+,16..100:+");
push(@exp_seq_coords_A, "3..12:+,16..20:+,26..110:+");

push(@desc_A,           "2 deletes");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..94:+");
push(@ins_str_A,        "BLASTNULL");
push(@del_str_A,        "Q12:S10-5;Q20:S23-3");
push(@exp_mdl_coords_A, "1..10:+,16..23:+,27..100:+");
push(@exp_seq_coords_A, "3..12:+,13..20:+,21..94:+");

push(@desc_A,           "1 insert len 3, 1 delete len 3");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..102:+");
push(@ins_str_A,        "Q12:S10+3");
push(@del_str_A,        "Q22:S17-3");
push(@exp_mdl_coords_A, "1..10:+,11..17:+,21..100:+");
push(@exp_seq_coords_A, "3..12:+,16..22:+,23..102:+");

push(@desc_A,           "1 delete len 3, 1 insert len 3");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..102:+");
push(@ins_str_A,        "Q22:S23+3");
push(@del_str_A,        "Q12:S10-3");
push(@exp_mdl_coords_A, "1..10:+,14..23:+,24..100:+");
push(@exp_seq_coords_A, "3..12:+,13..22:+,26..102:+");

push(@desc_A,           "1 insert len 3, 1 delete len 3, 1 insert len 1");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..103:+");
push(@ins_str_A,        "Q12:S10+3;Q30:S28+1",);
push(@del_str_A,        "Q22:S17-3;");
push(@exp_mdl_coords_A, "1..10:+,11..17:+,21..28:+,29..100:+");
push(@exp_seq_coords_A, "3..12:+,16..22:+,23..30:+,32..103:+");

push(@desc_A,           "1 delete len 3, 1 insert len 3, 1 delete len 1");
push(@in_mdl_coords_A,  "1..100:+");
push(@in_seq_coords_A,  "3..101:+");
push(@ins_str_A,        "Q22:S23+3");
push(@del_str_A,        "Q12:S10-3;Q30:S28-1");
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

###########################################
# vdr_CoordsMaxLengthSegment() tests
###########################################

@desc_A             = ();
my @coords_A        = ();
my @exp_coords_A    = ();
my @exp_rc_coords_A = ();
my @exp_len_A       = ();

push(@desc_A,           "1 sgm");
push(@coords_A,         "1..3:+");
push(@exp_coords_A,     "1..3:+");   
push(@exp_rc_coords_A,  ""); # exp_rc_coords is just reverse complement of exp_coords
push(@exp_len_A,        "3");   

push(@desc_A,           "2 sgm 1");
push(@coords_A,         "1..3:+,1..2:+");
push(@exp_coords_A,     "1..3:+");   
push(@exp_rc_coords_A,  "");   
push(@exp_len_A,        "3");   

push(@desc_A,           "2 sgm 2");
push(@coords_A,         "1..3:+,1..4:+");
push(@exp_coords_A,     "1..4:+");   
push(@exp_rc_coords_A,  "");   
push(@exp_len_A,        "4");   

push(@desc_A,           "2 sgm 3");
push(@coords_A,         "1..3:+,1..4:-");
push(@exp_coords_A,     "1..4:-");   
push(@exp_rc_coords_A,  "");   
push(@exp_len_A,        "4");   

push(@desc_A,           "3 sgm 1");
push(@coords_A,         "1..3:+,1..4:+,20..21:+");
push(@exp_coords_A,     "1..4:+");   
push(@exp_rc_coords_A,  "");   
push(@exp_len_A,        "4");   

push(@desc_A,           "3 sgm 2");
push(@coords_A,         "1..3:+,4..1:-,20..21:+");
push(@exp_coords_A,     "4..1:-");   
push(@exp_rc_coords_A,  "");   
push(@exp_len_A,        "4");   

push(@desc_A,           "3 sgm 3");
push(@coords_A,         "1..3:+,4..6:+,20..21:+");
push(@exp_coords_A,     "1..3:+");   
push(@exp_rc_coords_A,  "6..4:-");   
push(@exp_len_A,        "3");   

$ntests = scalar(@desc_A);
my ($cur_coords, $cur_len, $rc_coords, $rc_exp_coords);
for($i = 0; $i < $ntests; $i++) { 
  ($cur_coords, $cur_len) = vdr_CoordsMaxLengthSegment($coords_A[$i], undef);
  is($cur_coords, $exp_coords_A[$i], "vdr_CoordsMaxLengthSegment() coords: $desc_A[$i]");
  is($cur_len, $exp_len_A[$i], "vdr_CoordsMaxLengthSegment() length: $desc_A[$i]");

  # sanity check, reverse complement it and check again
  $rc_coords     = vdr_CoordsReverseComplement($coords_A[$i], 0, undef);
  $rc_exp_coords = ($exp_rc_coords_A[$i] eq "") ?
      vdr_CoordsReverseComplement($exp_coords_A[$i], 0, undef) :
      $exp_rc_coords_A[$i];
  ($cur_coords, $cur_len) = vdr_CoordsMaxLengthSegment($rc_coords, undef);
  is($cur_coords, $rc_exp_coords, "vdr_CoordsMaxLengthSegment() coords reverse complement: $desc_A[$i]");
  is($cur_len, $exp_len_A[$i], "vdr_CoordsMaxLengthSegment() len reverse complement: $desc_A[$i]");

  # another sanity check, reverse complement again and make sure we get back orig
  $cur_coords = vdr_CoordsReverseComplement($rc_coords, 0, undef);
  is($cur_coords, $coords_A[$i], "vdr_CoordsMaxLengthSegment() coords double reverse complement 1: $desc_A[$i]");
}

###########################################
# join_alignments_helper() tests
###########################################

@desc_A                  = ();
my @ali_5p_seq_coords_A  = ();
my @ali_5p_seq_A         = ();
my @ali_5p_mdl_A         = ();
my @ali_3p_seq_coords_A  = ();
my @ali_3p_seq_A         = ();
my @ali_3p_mdl_A         = ();
my @ugp_mdl_coords_A     = ();
my @ugp_seq_coords_A     = ();
my @ugp_seq_A            = ();       
my @seq_len_A            = ();       
my @mdl_len_A            = ();       
my @exp_joined_seq_A     = ();
my @exp_joined_mdl_A     = ();

push(@desc_A,               "ungapped");    
push(@ali_5p_seq_coords_A,  "1..20:+");  
push(@ali_5p_seq_A,         "AAAAAAAAAACCCCCCCCCC");
push(@ali_5p_mdl_A,         "xxxxxxxxxxxxxxxxxxxx");
push(@ali_3p_seq_coords_A,  "41..70:+");  
push(@ali_3p_seq_A,         "GGGGGGGGGGUUUUUUUUUUAAAAAAAAAA");
push(@ali_3p_mdl_A,         "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
push(@ugp_seq_coords_A,     "11..50:+");  
push(@ugp_mdl_coords_A,     "11..50:+");  
push(@ugp_seq_A,            "CCCCCCCCCCACGUACGUACGUACGUACGUGGGGGGGGGG");
push(@seq_len_A,            "70");
push(@mdl_len_A,            "70");
push(@exp_joined_seq_A,     "AAAAAAAAAACCCCCCCCCCACGUACGUACGUACGUACGUGGGGGGGGGGUUUUUUUUUUAAAAAAAAAA");
push(@exp_joined_mdl_A,     "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

$ntests = scalar(@desc_A);
my ($cur_joined_seq, $cur_joined_mdl);
for($i = 0; $i < $ntests; $i++) { 
  ($cur_joined_seq, $cur_joined_mdl) =
      join_alignments_helper($ali_5p_seq_coords_A[$i],
                             $ali_5p_seq_A[$i],
                             $ali_5p_mdl_A[$i],
                             $ali_3p_seq_coords_A[$i],
                             $ali_3p_seq_A[$i],
                             $ali_3p_mdl_A[$i],
                             $ugp_seq_coords_A[$i],
                             $ugp_mdl_coords_A[$i],
                             $ugp_seq_A[$i],
                             $seq_len_A[$i],
                             $mdl_len_A[$i],
                             undef);
  is($cur_joined_seq, $exp_joined_seq_A[$i], "join_alignments() correctly joined seq: $desc_A[$i]");
  is($cur_joined_mdl, $exp_joined_mdl_A[$i], "join_alignments() correctly joined seq: $desc_A[$i]");
}
