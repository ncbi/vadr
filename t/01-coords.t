use strict;
use warnings FATAL => 'all';
use Test::More tests => 165;

BEGIN {
    use_ok( 'vadr' ) || print "Bail out!\n";
}

my @desc_A     = ();  # array of descriptions for each test
my @exp_val_A  = ();  # array of return values for each test
my @exp_val2_A = ();  # array of return values for each test
my $ntests     = 0;      # number of tests in a section
my $cur_val    = undef; # current return value
###########################################
# vdr_CoordsReverseComplement() tests
# vdr_CoordsTokenReverseComplement() tests
#
# We don't have negative strand versions
# of some of these because we test that
# by reverse complementing the reverse
# complement and make sure we get the
# original coords back again.
#
# vdr_CoordsReverseComplement() calls
# vdr_CoordsTokenReverseComplement() 
# internally also.
##########################################
my @coords_A = (); # coords value to reverse complement
@desc_A = ();
@exp_val_A = ();
@exp_val2_A = ();

push(@desc_A,     "+");
push(@coords_A,   "1..200:+");
push(@exp_val_A,  "200..1:-");
push(@exp_val2_A, "200..1:-");

push(@desc_A,     "<+");
push(@coords_A,   "<1..200:+");
push(@exp_val_A,  "200..1:-");
push(@exp_val2_A, "200..>1:-");

push(@desc_A,     "+>");
push(@coords_A,   "100..>200:+");
push(@exp_val_A,  "200..100:-");
push(@exp_val2_A, "<200..100:-");

push(@desc_A,     "<+>");
push(@coords_A,   "<1..>200:+");
push(@exp_val_A,  "200..1:-");
push(@exp_val2_A, "<200..>1:-");

push(@desc_A,     "++");
push(@coords_A,   "1..200:+,300..400:+");
push(@exp_val_A,  "400..300:-,200..1:-");
push(@exp_val2_A, "400..300:-,200..1:-");

push(@desc_A,     "+><+");
push(@coords_A,   "1..>200:+,<300..400:+");
push(@exp_val_A,  "400..300:-,200..1:-");
push(@exp_val2_A, "400..>300:-,<200..1:-");

push(@desc_A,     "+-");
push(@coords_A,   "1..200:+,400..300:-");
push(@exp_val_A,  "300..400:+,200..1:-");
push(@exp_val2_A, "300..400:+,200..1:-");

push(@desc_A,     "+><-");
push(@coords_A,   "1..>200:+,<400..300:-");
push(@exp_val_A,  "300..400:+,200..1:-");
push(@exp_val2_A, "300..>400:+,<200..1:-");

push(@desc_A,     "+++");
push(@coords_A,   "1..200:+,300..400:+,402..700:+");
push(@exp_val_A,  "700..402:-,400..300:-,200..1:-");
push(@exp_val2_A, "700..402:-,400..300:-,200..1:-");

push(@desc_A,     "<+><+><+>");
push(@coords_A,   "<1..>200:+,<300..>400:+,<402..>700:+");
push(@exp_val_A,  "700..402:-,400..300:-,200..1:-");
push(@exp_val2_A, "<700..>402:-,<400..>300:-,<200..>1:-");

push(@desc_A,     "+-+");
push(@coords_A,   "1..200:+,400..300:-,402..700:+");
push(@exp_val_A,  "700..402:-,300..400:+,200..1:-");
push(@exp_val2_A, "700..402:-,300..400:+,200..1:-");

push(@desc_A,     "<+><-><+>");
push(@coords_A,   "<1..>200:+,<400..>300:-,<402..>700:+");
push(@exp_val_A,  "700..402:-,300..400:+,200..1:-");
push(@exp_val2_A, "<700..>402:-,<300..>400:+,<200..>1:-");

push(@desc_A,     "+-+-");
push(@coords_A,   "1..200:+,400..300:-,402..700:+,1000..650:-");
push(@exp_val_A,  "650..1000:+,700..402:-,300..400:+,200..1:-");
push(@exp_val2_A, "650..1000:+,700..402:-,300..400:+,200..1:-");

push(@desc_A,     "+->+<-");
push(@coords_A,   "1..200:+,400..>300:-,402..700:+,<1000..650:-");
push(@exp_val_A,  "650..1000:+,700..402:-,300..400:+,200..1:-");
push(@exp_val2_A, "650..>1000:+,700..402:-,<300..400:+,200..1:-");

$ntests = scalar(@desc_A);
my $do_carrots = undef;
my $cur_exp_val = undef;
# for each coords string, make sure that:
# vdr_CoordsReverseComplement works with $do_carrots as 0 and 1
# vdr_CoordsTokenReverseComplement works with $do_carrots as 0 and 1 if it is a single segment
# And that reverse complementing the reverse complement gives you back the original.
for(my $i = 0; $i < $ntests; $i++) { 
  $do_carrots = 0;
  $cur_val = vdr_CoordsReverseComplement($coords_A[$i], $do_carrots, undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsReverseComplement() no carrots: $desc_A[$i]");
  # now reverse complement the reverse complement and you should get back orig coords, minus carrots
  $cur_val = vdr_CoordsReverseComplement($cur_val, $do_carrots, undef);
  $cur_exp_val = $coords_A[$i]; 
  $cur_exp_val =~ s/\>//g;
  $cur_exp_val =~ s/\<//g;
  is($cur_val, $cur_exp_val, "vdr_CoordsReverseComplement() no carrots squared: $desc_A[$i]");

  if($coords_A[$i] !~ m/,/) { # single segment, test vdr_CoordsTokenReverseComplement too
    $cur_val = vdr_CoordsTokenReverseComplement($coords_A[$i], $do_carrots, undef);
    is($cur_val, $exp_val_A[$i], "vdr_CoordsTokenReverseComplement() no carrots: $desc_A[$i]");
    $cur_val = vdr_CoordsTokenReverseComplement($cur_val, $do_carrots, undef);
    is($cur_val, $cur_exp_val, "vdr_CoordsTokenReverseComplement() no carrots squared: $desc_A[$i]");
  }

  $do_carrots = 1;
  $cur_val = vdr_CoordsReverseComplement($coords_A[$i], $do_carrots, undef);
  is($cur_val, $exp_val2_A[$i], "vdr_CoordsReverseComplement() with carrots: $desc_A[$i]");
  # now reverse complement the reverse complement and you should get back orig coords, minus carrots
  $cur_val = vdr_CoordsReverseComplement($cur_val, $do_carrots, undef);
  $cur_exp_val = $coords_A[$i]; 
  is($cur_val, $cur_exp_val, "vdr_CoordsReverseComplement() with carrots squared: $desc_A[$i]");

  if($coords_A[$i] !~ m/,/) { # single segment, test vdr_CoordsTokenReverseComplement too
    $cur_val = vdr_CoordsTokenReverseComplement($coords_A[$i], $do_carrots, undef);
    is($cur_val, $exp_val2_A[$i], "vdr_CoordsTokenReverseComplement() with carrots: $desc_A[$i]");
    $cur_val = vdr_CoordsTokenReverseComplement($cur_val, $do_carrots, undef);
    is($cur_val, $cur_exp_val, "vdr_CoordsTokenReverseComplement() with carrots squared: $desc_A[$i]");
  }
}

################################
# vdr_CoordsFromLocation() tests
################################
my @location_A = (); # location value from a GenBank file 
@desc_A = ();
@exp_val_A = ();
@exp_val2_A = ();

push(@desc_A,      "+");
push(@location_A,  "1..200");
push(@exp_val_A,   "1..200:+");
push(@exp_val2_A,  "1..200:+");

push(@desc_A,      "<+");
push(@location_A,  "<1..200");
push(@exp_val_A,   "1..200:+");
push(@exp_val2_A,  "<1..200:+");

push(@desc_A,     "+>");
push(@location_A, "100..>200");
push(@exp_val_A,  "100..200:+");
push(@exp_val2_A, "100..>200:+");

push(@desc_A,     "<+>");
push(@location_A, "<100..>200");
push(@exp_val_A,  "100..200:+");
push(@exp_val2_A, "<100..>200:+");

push(@desc_A,     "complement(+)");
push(@location_A, "complement(1..200)");
push(@exp_val_A,  "200..1:-");
push(@exp_val2_A, "200..1:-");

push(@desc_A,     "complement(<+)");
push(@location_A, "complement(<1..200)");
push(@exp_val_A,  "200..1:-");
push(@exp_val2_A, "200..>1:-");

push(@desc_A,     "complement(+>)");
push(@location_A, "complement(100..>200)");
push(@exp_val_A,  "200..100:-");
push(@exp_val2_A, "<200..100:-");

push(@desc_A,     "complement(<+>)");
push(@location_A, "complement(<100..>200)");
push(@exp_val_A,  "200..100:-");
push(@exp_val2_A, "<200..>100:-");

push(@desc_A,     "join(++)");
push(@location_A, "join(1..200,300..400)");
push(@exp_val_A,  "1..200:+,300..400:+");
push(@exp_val2_A, "1..200:+,300..400:+");

push(@desc_A,     "join(<++)");
push(@location_A, "join(<1..200,300..400)");
push(@exp_val_A,  "1..200:+,300..400:+");
push(@exp_val2_A, "<1..200:+,300..400:+");

push(@desc_A,     "join(<++>)");
push(@location_A, "join(<1..200,300..>400)");
push(@exp_val_A,  "1..200:+,300..400:+");
push(@exp_val2_A, "<1..200:+,300..>400:+");

push(@desc_A,     "join(<+>+)");
push(@location_A, "join(<1..200,<300..400)");
push(@exp_val_A,  "1..200:+,300..400:+");
push(@exp_val2_A, "<1..200:+,<300..400:+");

push(@desc_A,     "join(<+><+>)");
push(@location_A, "join(<1..>200,<300..>400)");
push(@exp_val_A,  "1..200:+,300..400:+");
push(@exp_val2_A, "<1..>200:+,<300..>400:+");

push(@desc_A,     "complement(join(++)");
push(@location_A, "complement(join(1..200,300..400))");
push(@exp_val_A,  "400..300:-,200..1:-");
push(@exp_val2_A, "400..300:-,200..1:-");

push(@desc_A,     "complement(join(<++>)");
push(@location_A, "complement(join(<1..200,300..>400))");
push(@exp_val_A,  "400..300:-,200..1:-");
push(@exp_val2_A, "<400..300:-,200..>1:-");

push(@desc_A,     "complement(join(<+<+)");
push(@location_A, "complement(join(<1..200,<300..400))");
push(@exp_val_A,  "400..300:-,200..1:-");
push(@exp_val2_A, "400..>300:-,200..>1:-");

push(@desc_A,     "join(+,complement(+))");
push(@location_A, "join(1..200,complement(300..400))");
push(@exp_val_A,  "1..200:+,400..300:-");
push(@exp_val2_A, "1..200:+,400..300:-");

push(@desc_A,     "join(<+,complement(+>))");
push(@location_A, "join(<1..200,complement(300..>400))");
push(@exp_val_A,  "1..200:+,400..300:-");
push(@exp_val2_A, "<1..200:+,<400..300:-");

push(@desc_A,     "join(complement(+),+)");
push(@location_A, "join(complement(300..400),100..200)");
push(@exp_val_A,  "400..300:-,100..200:+");
push(@exp_val2_A, "400..300:-,100..200:+");

push(@desc_A,     "join(complement(+>),<+)");
push(@location_A, "join(complement(300..>400),<100..200)");
push(@exp_val_A,  "400..300:-,100..200:+");
push(@exp_val2_A, "<400..300:-,<100..200:+");

$ntests = scalar(@desc_A);
for(my $i = 0; $i < $ntests; $i++) { 
  $do_carrots = 0;
  $cur_val = vdr_CoordsFromLocation($location_A[$i], $do_carrots, undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsFromLocation(): $desc_A[$i]");

  $do_carrots = 1;
  $cur_val = vdr_CoordsFromLocation($location_A[$i], $do_carrots, undef);
  is($cur_val, $exp_val2_A[$i], "vdr_CoordsFromLocationWithCarrots(): $desc_A[$i]");
}

############################################
# vdr_CoordsRelativeToAbsolute() tests
############################################
# ABSOLUTE: POSITIVE STRAND, 1-3 segments
# RELATIVE: POSITIVE STRAND, 1 segment
# positive strand, 1 segment
@desc_A                = ();
my @full_abs_nt_coords_A  = ();
my @rel_nt_or_aa_coords_A = ();
my @rel_is_aa_A           = ();
@exp_val_A             = ();

push(@desc_A,                "abs (+), rel (+)");
push(@full_abs_nt_coords_A,  "11..100:+");
push(@rel_nt_or_aa_coords_A, "6..38:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..48:+");   

# positive strand, 2 segments
push(@desc_A,                "abs (+)+, rel (+) ");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "6..25:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..35:+");   

push(@desc_A,                "abs +(+), rel (+)");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "31..33:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "42..44:+");   

push(@desc_A,                "abs (++), rel (+)");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "6..38:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..40:+,42..49:+");   

# positive strand, 3 segments
push(@desc_A,                "abs (+)++, rel (+)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "6..25:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..35:+");   

push(@desc_A,                "abs +(+)+, rel (+)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "31..33:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "42..44:+");   

push(@desc_A,                "abs ++(+), rel (+)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "90..105:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "105..120:+");   

push(@desc_A,                "abs (++)+, rel (+)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "6..38:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..40:+,42..49:+");   

push(@desc_A,                "abs +(++), rel (+)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "39..99:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "50..100:+,105..114:+");   

push(@desc_A,                "abs (+++), rel (+)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "19..99:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "29..40:+,42..100:+,105..114:+");   

# ABSOLUTE: POSITIVE STRAND, 1-3 segments
# RELATIVE: POSITIVE STRAND, 2 segments
# note: all relative segments must be 'spanned' in return segment

push(@desc_A,                "abs (+), rel (++)");
push(@full_abs_nt_coords_A,  "11..100:+");
push(@rel_nt_or_aa_coords_A, "6..38:+,45..50:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..48:+,55..60:+");   

push(@desc_A,                "abs (+), rel (+++)");
push(@full_abs_nt_coords_A,  "11..100:+");
push(@rel_nt_or_aa_coords_A, "6..38:+,45..50:+,53..59:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..48:+,55..60:+,63..69:+");   

push(@desc_A,                "abs (++), rel (++)");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "6..38:+,45..50:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..40:+,42..49:+,56..61:+");   

push(@desc_A,                "abs (++), rel (+++)");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "6..38:+,45..50:+,53..59:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..40:+,42..49:+,56..61:+,64..70:+");   

push(@desc_A,                "abs (+++), rel (++)");
push(@full_abs_nt_coords_A,  "11..40:+,42..58:+,62..104:+");
push(@rel_nt_or_aa_coords_A, "6..38:+,45..50:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..40:+,42..49:+,56..58:+,62..64:+");   

push(@desc_A,                "abs (+++), rel (+++)");
push(@full_abs_nt_coords_A,  "11..40:+,42..58:+,62..104:+");
push(@rel_nt_or_aa_coords_A, "6..38:+,45..50:+,53..59:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "16..40:+,42..49:+,56..58:+,62..64:+,67..73:+");   

# ABSOLUTE: POSITIVE STRAND, 1-3 segments
# RELATIVE: NEGATIVE STRAND, 1 segment
push(@desc_A,                "abs (+), rel (-)");
push(@full_abs_nt_coords_A,  "11..100:+");
push(@rel_nt_or_aa_coords_A, "38..6:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "48..16:-");   

# positive strand, 2 segments
push(@desc_A,                "abs (+)+, rel (-) ");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "25..6:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "35..16:-");   

push(@desc_A,                "abs +(+), rel (-)");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "33..31:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "44..42:-");   

push(@desc_A,                "abs (++), rel (-)");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "38..6:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "49..42:-,40..16:-");

# positive strand, 3 segments
push(@desc_A,                "abs (+)++, rel (-)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "25..6:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "35..16:-");   

push(@desc_A,                "abs +(+)+, rel (-)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "33..31:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "44..42:-");   

push(@desc_A,                "abs ++(+), rel (-)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "105..90:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "120..105:-");   

push(@desc_A,                "abs (++)+, rel (-)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "38..6:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "49..42:-,40..16:-");

push(@desc_A,                "abs +(++), rel (-)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "99..39:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "114..105:-,100..50:-");

push(@desc_A,                "abs (+++), rel (-)");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "99..19:-");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@exp_val_A,             "114..105:-,100..42:-,40..29:-");



$ntests = scalar(@desc_A);
for(my $i = 0; $i < $ntests; $i++) { 
  $cur_val = vdr_CoordsRelativeToAbsolute($full_abs_nt_coords_A[$i], 
                                             $rel_nt_or_aa_coords_A[$i], 
                                             $rel_is_aa_A[$i], undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsRelativeToAbsolute(): $desc_A[$i]");
  # sanity check:
  # make sure the length of the returned coords string is the same as the
  # length of the $rel_nt_or_aa_coords string
  my $rel_nt_length   = vdr_CoordsLength($rel_nt_or_aa_coords_A[$i], undef);
  my $cur_val_length  = vdr_CoordsLength($cur_val, undef);
  is($cur_val_length, $rel_nt_length, "vdr_CoordsRelativeToAbsolute() and vdr_CoordsLength(): length sanity check for $desc_A[$i]");
}
