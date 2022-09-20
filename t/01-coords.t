use strict;
use warnings FATAL => 'all';
use Test::More tests => 554;

BEGIN {
    use_ok( 'vadr' )      || print "Bail out!\n";
    use_ok( 'vadr_seed' ) || print "Bail out!\n";
}

my @desc_A     = ();  # array of descriptions for each test
my @exp_val_A  = ();  # array of return values for each test
my @exp_val2_A = ();  # array of return values for each test
my $ntests     = 0;      # number of tests in a section
my $cur_val    = undef; # current return value
my $i;
my @sgm1_A = ();   # segment 1 values
my @sgm2_A = ();   # segment 2 values
my @coords_A = (); # coords values
my $do_carrots = undef;
my $cur_exp_val = undef;
my @location_A = (); # location values
my ($cur_sgm1, $cur_sgm2, $cur_coords);
my @abs_coords_A  = ();
my @rel_coords_A = ();
my @rel_pt_coords_A  = ();
my ($rel_length, $cur_val_length, $cur_rel_coords, $cur_rel_coord);
my @full_coords_A   = ();
my @subseq_coords_A = ();
my ($fract_start, $fract_stop, $ret_start, $ret_stop, $ret_strand);

###########################################
# vdr_CoordsReverseComplement() tests
# vdr_CoordsSegmentReverseComplement() tests
#
# We don't have negative strand versions
# of some of these because we test that
# by reverse complementing the reverse
# complement and make sure we get the
# original coords back again.
#
# vdr_CoordsReverseComplement() calls
# vdr_CoordsSegmentReverseComplement() 
# internally also.
##########################################
@coords_A = (); # coords value to reverse complement
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
$do_carrots = undef;
$cur_exp_val = undef;
# for each coords string, make sure that:
# vdr_CoordsReverseComplement works with $do_carrots as 0 and 1
# vdr_CoordsSegmentReverseComplement works with $do_carrots as 0 and 1 if it is a single segment
# And that reverse complementing the reverse complement gives you back the original.
for($i = 0; $i < $ntests; $i++) { 
  $do_carrots = 0;
  $cur_val = vdr_CoordsReverseComplement($coords_A[$i], $do_carrots, undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsReverseComplement() no carrots: $desc_A[$i]");
  # now reverse complement the reverse complement and you should get back orig coords, minus carrots
  $cur_val = vdr_CoordsReverseComplement($cur_val, $do_carrots, undef);
  $cur_exp_val = $coords_A[$i]; 
  $cur_exp_val =~ s/\>//g;
  $cur_exp_val =~ s/\<//g;
  is($cur_val, $cur_exp_val, "vdr_CoordsReverseComplement() no carrots squared: $desc_A[$i]");

  if($coords_A[$i] !~ m/,/) { # single segment, test vdr_CoordsSegmentReverseComplement too
    $cur_val = vdr_CoordsSegmentReverseComplement($coords_A[$i], $do_carrots, undef);
    is($cur_val, $exp_val_A[$i], "vdr_CoordsSegmentReverseComplement() no carrots: $desc_A[$i]");
    $cur_val = vdr_CoordsSegmentReverseComplement($cur_val, $do_carrots, undef);
    is($cur_val, $cur_exp_val, "vdr_CoordsSegmentReverseComplement() no carrots squared: $desc_A[$i]");
  }

  $do_carrots = 1;
  $cur_val = vdr_CoordsReverseComplement($coords_A[$i], $do_carrots, undef);
  is($cur_val, $exp_val2_A[$i], "vdr_CoordsReverseComplement() with carrots: $desc_A[$i]");
  # now reverse complement the reverse complement and you should get back orig coords, minus carrots
  $cur_val = vdr_CoordsReverseComplement($cur_val, $do_carrots, undef);
  $cur_exp_val = $coords_A[$i]; 
  is($cur_val, $cur_exp_val, "vdr_CoordsReverseComplement() with carrots squared: $desc_A[$i]");

  if($coords_A[$i] !~ m/,/) { # single segment, test vdr_CoordsSegmentReverseComplement too
    $cur_val = vdr_CoordsSegmentReverseComplement($coords_A[$i], $do_carrots, undef);
    is($cur_val, $exp_val2_A[$i], "vdr_CoordsSegmentReverseComplement() with carrots: $desc_A[$i]");
    $cur_val = vdr_CoordsSegmentReverseComplement($cur_val, $do_carrots, undef);
    is($cur_val, $cur_exp_val, "vdr_CoordsSegmentReverseComplement() with carrots squared: $desc_A[$i]");
  }
}

################################
# vdr_CoordsFromLocation() tests
################################
@location_A = (); # location value from a GenBank file 
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
for($i = 0; $i < $ntests; $i++) { 
  $do_carrots = 0;
  $cur_val = vdr_CoordsFromLocation($location_A[$i], $do_carrots, undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsFromLocation(): $desc_A[$i]");

  $do_carrots = 1;
  $cur_val = vdr_CoordsFromLocation($location_A[$i], $do_carrots, undef);
  is($cur_val, $exp_val2_A[$i], "vdr_CoordsFromLocationWithCarrots(): $desc_A[$i]");
}

###########################################
# vdr_CoordsMergeAllAdjacentSegments() tests
# vdr_CoordsMergeTwoSegmentsIfAdjacent() tests
#
# We don't have negative strand versions
# of some of these because we test that
# by reverse complementing the reverse
# complement and make sure we get the
# original value back again.
##########################################
@sgm1_A = ();   # segment 1
@sgm2_A = ();   # segment 2
@desc_A = ();
@coords_A = ();
@exp_val_A = ();
@exp_val2_A = ();

push(@desc_A,     "+ +, not adj");
push(@coords_A,   "1..200:+,202..300:+");
push(@sgm1_A,     "1..200:+");
push(@sgm2_A,     "202..300:+");
push(@exp_val_A,  "1..200:+,202..300:+");
push(@exp_val2_A, "");

push(@desc_A,     "++, adj");
push(@coords_A,   "1..200:+,201..300:+");
push(@sgm1_A,     "1..200:+");
push(@sgm2_A,     "201..300:+");
push(@exp_val_A,  "1..300:+");
push(@exp_val2_A, "1..300:+");

push(@desc_A,     "+ + +, not adj");
push(@coords_A,   "1..200:+,202..300:+,302..310:+");
push(@sgm1_A,     "202..300:+");
push(@sgm2_A,     "302..310:+");
push(@exp_val_A,  "1..200:+,202..300:+,302..310:+");
push(@exp_val2_A, "");

push(@desc_A,     "++ +, 2 adj");
push(@coords_A,   "1..200:+,201..300:+,302..310:+");
push(@sgm1_A,     "201..300:+");
push(@sgm2_A,     "302..310:+");
push(@exp_val_A,  "1..300:+,302..310:+");
push(@exp_val2_A, "");

push(@desc_A,     "+++, 3 adj");
push(@coords_A,   "1..200:+,201..300:+,301..310:+");
push(@sgm1_A,     "201..300:+");
push(@sgm2_A,     "301..310:+");
push(@exp_val_A,  "1..310:+");
push(@exp_val2_A, "201..310:+");

push(@desc_A,     "+ ++ +, 2 adj");
push(@coords_A,   "1..200:+,202..300:+,301..310:+,500..700:+");
push(@sgm1_A,     "301..310:+");
push(@sgm2_A,     "500..700:+");
push(@exp_val_A,  "1..200:+,202..310:+,500..700:+");
push(@exp_val2_A, "");

push(@desc_A,     "+++++, 5 adj");
push(@coords_A,   "1..1:+,2..2:+,3..3:+,4..4:+,5..5:+");
push(@sgm1_A,     "3..3:+");
push(@sgm2_A,     "4..4:+");
push(@exp_val_A,  "1..5:+");
push(@exp_val2_A, "3..4:+");

# mixed strands
push(@desc_A,     "+ -, not adj");
push(@coords_A,   "1..200:+,300..202:-");
push(@sgm1_A,     "1..200:+");
push(@sgm2_A,     "300..202:-");
push(@exp_val_A,  "1..200:+,300..202:-");
push(@exp_val2_A, "");

push(@desc_A,     "+-, adj");
push(@coords_A,   "1..200:+,300..201:-");
push(@sgm1_A,     "1..200:+");
push(@sgm2_A,     "300..201:-");
push(@exp_val_A,  "1..200:+,300..201:-");
push(@exp_val2_A, "");

push(@desc_A,     "+ - +, not adj");
push(@coords_A,   "1..200:+,300..202:-,302..310:+");
push(@sgm1_A,     "300..202:-");
push(@sgm2_A,     "302..310:+");
push(@exp_val_A,  "1..200:+,300..202:-,302..310:+");
push(@exp_val2_A, "");

push(@desc_A,     "++ -, 2 adj");
push(@coords_A,   "1..200:+,201..300:+,310..301:-");
push(@sgm1_A,     "201..300:+");
push(@sgm2_A,     "310..301:-");
push(@exp_val_A,  "1..300:+,310..301:-");
push(@exp_val2_A, "");

push(@desc_A,     "++-++, 2 adj");
push(@coords_A,   "1..200:+,201..300:+,310..301:-,311..400:+,401..500:+");
push(@sgm1_A,     "310..301:-");
push(@sgm2_A,     "311..400:+");
push(@exp_val_A,  "1..300:+,310..301:-,311..500:+");
push(@exp_val2_A, "");

# different directions, same strand 
push(@desc_A,     "++ fwd bck");
push(@coords_A,   "1..200:+,201..1:+");
push(@sgm1_A,     "1..200:+");
push(@sgm2_A,     "201..1:+");
push(@exp_val_A,  "1..200:+,201..1:+");
push(@exp_val2_A, "");

push(@desc_A,     "++++ fwd bck");
push(@coords_A,   "1..200:+,201..1:+,2..300:+,301..302:+");
push(@sgm1_A,     "201..1:+");
push(@sgm2_A,     "2..300:+");
push(@exp_val_A,  "1..200:+,201..1:+,2..302:+");
push(@exp_val2_A, "");

push(@desc_A,     "++++ fwd bck");
push(@coords_A,   "1..200:+,201..1:+,2..300:+,301..302:+");
push(@sgm1_A,     "201..1:+");
push(@sgm2_A,     "2..300:+");
push(@exp_val_A,  "1..200:+,201..1:+,2..302:+");
push(@exp_val2_A, "");

push(@desc_A,     "++, bck bck");
push(@coords_A,   "6..5:+,4..3:+,2..1:+");
push(@sgm1_A,     "6..5:+");
push(@sgm2_A,     "4..3:+");
push(@exp_val_A,  "6..1:+");
push(@exp_val2_A, "6..3:+");

push(@desc_A,     "++, bck bck");
push(@coords_A,   "6..5:+,4..3:+,3..1:+");
push(@sgm1_A,     "4..3:+");
push(@sgm2_A,     "3..1:+");
push(@exp_val_A,  "6..3:+,3..1:+");
push(@exp_val2_A, "");


$ntests = scalar(@desc_A);
for($i = 0; $i < $ntests; $i++) { 
  $cur_val = vdr_CoordsMergeAllAdjacentSegments($coords_A[$i], undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsMergeAllAdjacentSegments(): $desc_A[$i]");
  # now reverse complement $coords_A[$i], run the subroutine again, revcomp the result and 
  # make sure we get the expected return value again
  $cur_coords = vdr_CoordsReverseComplement($coords_A[$i], 0, undef);
  $cur_val = vdr_CoordsReverseComplement(vdr_CoordsMergeAllAdjacentSegments($cur_coords, undef), 0, undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsMergeAllAdjacentSegments() rev comp: $desc_A[$i]");

  $cur_val = vdr_CoordsMergeTwoSegmentsIfAdjacent($sgm1_A[$i], $sgm2_A[$i], undef);
  is($cur_val, $exp_val2_A[$i], "vdr_CoordsMergeTwoSegmentsIfAdjacent(): $desc_A[$i]");
  # now reverse complement $sgm1_A[$i] and $sgm2_A[$i], run the subroutine again, revcomp the result and 
  # make sure we get the expected return value again
  $cur_sgm1 = vdr_CoordsReverseComplement($sgm1_A[$i], 0, undef);
  $cur_sgm2 = vdr_CoordsReverseComplement($sgm2_A[$i], 0, undef);
  # note: pass in cur_sgm2 first, then cur_sgm1 b/c they're reverse complemented
  $cur_val = vdr_CoordsReverseComplement(vdr_CoordsMergeTwoSegmentsIfAdjacent($cur_sgm2, $cur_sgm1, undef), 0, undef);
  is($cur_val, $exp_val2_A[$i], "vdr_CoordsMergeTwoSegmentsIfAdjacent() rev comp: $desc_A[$i]");
}

############################################
# vdr_CoordsRelativeToAbsolute() tests
############################################
@desc_A       = ();
@abs_coords_A = ();
@rel_coords_A = ();
@exp_val_A    = ();

push(@desc_A,       "abs (+), rel (+)");
push(@abs_coords_A, "11..100:+");
push(@rel_coords_A, "6..38:+");    
push(@exp_val_A,    "16..48:+");   

push(@desc_A,       "abs (+)+, rel (+) ");
push(@abs_coords_A, "11..40:+,42..101:+");
push(@rel_coords_A, "6..25:+");    
push(@exp_val_A,    "16..35:+");   

push(@desc_A,       "abs +(+), rel (+)");
push(@abs_coords_A, "11..40:+,42..101:+");
push(@rel_coords_A, "31..33:+");    
push(@exp_val_A,    "42..44:+");   

push(@desc_A,        "abs (++), rel (+)");
push(@abs_coords_A,  "11..40:+,42..101:+");
push(@rel_coords_A,  "6..38:+");    
push(@exp_val_A,     "16..40:+,42..49:+");   

push(@desc_A,        "abs (+)++, rel (+)");
push(@abs_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_coords_A,  "6..25:+");    
push(@exp_val_A,     "16..35:+");   

push(@desc_A,        "abs +(+)+, rel (+)");
push(@abs_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_coords_A,  "31..33:+");    
push(@exp_val_A,     "42..44:+");   

push(@desc_A,        "abs ++(+), rel (+)");
push(@abs_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_coords_A,  "90..105:+");    
push(@exp_val_A,     "105..120:+");   

push(@desc_A,        "abs (++)+, rel (+)");
push(@abs_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_coords_A,  "6..38:+");    
push(@exp_val_A,     "16..40:+,42..49:+");   

push(@desc_A,        "abs +(++), rel (+)");
push(@abs_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_coords_A,  "39..99:+");    
push(@exp_val_A,     "50..100:+,105..114:+");   

push(@desc_A,        "abs (+++), rel (+)");
push(@abs_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_coords_A,  "19..99:+");    
push(@exp_val_A,     "29..40:+,42..100:+,105..114:+");   

push(@desc_A,        "abs (+), rel (++)");
push(@abs_coords_A,  "11..100:+");
push(@rel_coords_A,  "6..38:+,45..50:+");    
push(@exp_val_A,     "16..48:+,55..60:+");   

push(@desc_A,        "abs (+), rel (+++)");
push(@abs_coords_A,  "11..100:+");
push(@rel_coords_A,  "6..38:+,45..50:+,53..59:+");    
push(@exp_val_A,     "16..48:+,55..60:+,63..69:+");   

push(@desc_A,        "abs (++), rel (++)");
push(@abs_coords_A,  "11..40:+,42..101:+");
push(@rel_coords_A,  "6..38:+,45..50:+");    
push(@exp_val_A,     "16..40:+,42..49:+,56..61:+");   

push(@desc_A,        "abs (++), rel (+++)");
push(@abs_coords_A,  "11..40:+,42..101:+");
push(@rel_coords_A,  "6..38:+,45..50:+,53..59:+");    
push(@exp_val_A,      "16..40:+,42..49:+,56..61:+,64..70:+");   

push(@desc_A,        "abs (+++), rel (++)");
push(@abs_coords_A,  "11..40:+,42..58:+,62..104:+");
push(@rel_coords_A,  "6..38:+,45..50:+");    
push(@exp_val_A,     "16..40:+,42..49:+,56..58:+,62..64:+");   

push(@desc_A,        "abs (+++), rel (+++)");
push(@abs_coords_A,  "11..40:+,42..58:+,62..104:+");
push(@rel_coords_A,  "6..38:+,45..50:+,53..59:+");    
push(@exp_val_A,     "16..40:+,42..49:+,56..58:+,62..64:+,67..73:+");   

# negative strand abs coords
push(@desc_A,       "abs (-), rel (+)");
push(@abs_coords_A, "100..11:-");
push(@rel_coords_A, "6..38:+");    
push(@exp_val_A,    "95..63:-");   

push(@desc_A,       "abs (-)-, rel (+) ");
push(@abs_coords_A, "101..42:-,40..11:-");
push(@rel_coords_A, "6..25:+");    
push(@exp_val_A,    "96..77:-");   

push(@desc_A,       "abs -(-), rel (+)");
push(@abs_coords_A, "101..42:-,40..11:-");
push(@rel_coords_A, "71..73:+");    
push(@exp_val_A,    "30..28:-");   

push(@desc_A,        "abs (--), rel (+)");
push(@abs_coords_A,  "101..42:-,40..11:-");
push(@rel_coords_A,  "6..78:+");    
push(@exp_val_A,     "96..42:-,40..23:-");   

push(@desc_A,        "abs (-)--, rel (+)");
push(@abs_coords_A,  "121..105:-,100..42:-,40..11:-");
push(@rel_coords_A,  "6..15:+");    
push(@exp_val_A,     "116..107:-");

push(@desc_A,        "abs -(-)-, rel (+)");
push(@abs_coords_A,  "121..105:-,100..42:-,40..11:-");
push(@rel_coords_A,  "71..73:+");
push(@exp_val_A,     "47..45:-");   

push(@desc_A,        "abs --(-), rel (+)");
push(@abs_coords_A,  "121..105:-,100..42:-,40..11:-");
push(@rel_coords_A,  "80..85:+");    
push(@exp_val_A,     "37..32:-");   

push(@desc_A,        "abs (--)-, rel (+)");
push(@abs_coords_A,  "121..105:-,100..42:-,40..11:-");
push(@rel_coords_A,  "6..25:+");    
push(@exp_val_A,     "116..105:-,100..93:-");   

push(@desc_A,        "abs -(--), rel (+)");
push(@abs_coords_A,  "121..105:-,100..42:-,40..11:-");
push(@rel_coords_A,  "39..99:+");    
push(@exp_val_A,     "79..42:-,40..18:-");   

push(@desc_A,        "abs (---), rel (+)");
push(@abs_coords_A,  "121..105:-,100..42:-,40..11:-");
push(@rel_coords_A,  "9..99:+");    
push(@exp_val_A,     "113..105:-,100..42:-,40..18:-");   

push(@desc_A,        "abs (+), rel (--)");
push(@abs_coords_A,  "100..11:-");
push(@rel_coords_A,  "6..38:+,45..50:+");    
push(@exp_val_A,     "95..63:-,56..51:-");

push(@desc_A,        "abs (+), rel (---)");
push(@abs_coords_A,  "100..11:-");
push(@rel_coords_A,  "6..38:+,45..50:+,53..59:+");    
push(@exp_val_A,     "95..63:-,56..51:-,48..42:-");

push(@desc_A,        "abs (++), rel (--)");
push(@abs_coords_A,  "101..42:-,40..11:-");
push(@rel_coords_A,  "6..38:+,45..70:+");    
push(@exp_val_A,     "96..64:-,57..42:-,40..31:-");

push(@desc_A,        "abs (++), rel (---)");
push(@abs_coords_A,  "101..42:-,40..11:-");
push(@rel_coords_A,  "6..38:+,45..70:+,73..80:+");    
push(@exp_val_A,     "96..64:-,57..42:-,40..31:-,28..21:-");

push(@desc_A,        "abs (+++), rel (--)");
push(@abs_coords_A,  "104..94:-,101..42:-,40..11:-");
push(@rel_coords_A,  "6..38:+,45..50:+");    
push(@exp_val_A,     "99..94:-,101..75:-,68..63:-");

push(@desc_A,        "abs (+++), rel (---)");
push(@abs_coords_A,  "104..94:-,101..42:-,40..11:-");
push(@rel_coords_A,  "6..38:+,45..50:+,53..59:+");    
push(@exp_val_A,     "99..94:-,101..75:-,68..63:-,60..54:-");

# Did not include these mixed strands tests yet, not sure if they'll
# ever be needed and I had already spent too much time on this I put a
# check in the code for mixed strand abs or rel coords (it fails if it
# finds them) and a note about not allowing it because I didn't write
# tests yet. 
# EPN, Fri Mar 20 16:12:17 2020
#push(@desc_A,        "abs (-++), rel (+++)");
#push(@desc_A,        "abs (+-+), rel (+++)");
#push(@desc_A,        "abs (++-), rel (+++)");
#push(@desc_A,        "abs (-+-), rel (+++)");
#push(@desc_A,        "abs (--+), rel (+++)");

#push(@desc_A,        "abs (-++), rel (+-+)");
#push(@desc_A,        "abs (+-+), rel (+-+)");
#push(@desc_A,        "abs (++-), rel (+-+)");
#push(@desc_A,        "abs (-+-), rel (+-+)");
#push(@desc_A,        "abs (--+), rel (+-+)");

#push(@desc_A,        "abs (-++), rel (+--)");
#push(@desc_A,        "abs (+-+), rel (+--)");
#push(@desc_A,        "abs (++-), rel (+--)");
#push(@desc_A,        "abs (-+-), rel (+--)");
#push(@desc_A,        "abs (--+), rel (+--)");

#push(@desc_A,        "abs (-++), rel (--+)");
#push(@desc_A,        "abs (+-+), rel (--+)");
#push(@desc_A,        "abs (++-), rel (--+)");
#push(@desc_A,        "abs (-+-), rel (--+)");
#push(@desc_A,        "abs (--+), rel (--+)");

$ntests = scalar(@desc_A);
for($i = 0; $i < $ntests; $i++) { 
  $cur_val = vdr_CoordsRelativeToAbsolute($abs_coords_A[$i], 
                                          $rel_coords_A[$i], undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsRelativeToAbsolute(): $desc_A[$i]");
  # sanity check:
  # make sure the length of the returned coords string is the same as the
  # length of the $rel_nt_or_aa_coords string
  $rel_length     = vdr_CoordsLength($rel_coords_A[$i], undef);
  $cur_val_length = vdr_CoordsLength($cur_val, undef);
  is($cur_val_length, $rel_length, "vdr_CoordsRelativeToAbsolute() and vdr_CoordsLength(): length sanity check for $desc_A[$i]");

  # test vdr_CoordsRelativeSingleCoordToAbsolute()
  # 5' most position
  $cur_rel_coord = vdr_Feature5pMostPosition($rel_coords_A[$i], undef);
  $cur_exp_val   = vdr_Feature5pMostPosition($exp_val_A[$i], undef);
  $cur_val       = vdr_CoordsRelativeSingleCoordToAbsolute($abs_coords_A[$i], $cur_rel_coord, undef);
  is($cur_val, $cur_exp_val, "vdr_CoordsRelativeSingleCoordToAbsolute() 5'-most position: $desc_A[$i]");
  # 3' most position
  $cur_rel_coord = vdr_Feature3pMostPosition($rel_coords_A[$i], undef);
  $cur_exp_val   = vdr_Feature3pMostPosition($exp_val_A[$i], undef);
  $cur_val       = vdr_CoordsRelativeSingleCoordToAbsolute($abs_coords_A[$i], $cur_rel_coord, undef);
  is($cur_val, $cur_exp_val, "vdr_CoordsRelativeSingleCoordToAbsolute() 3'-most position: $desc_A[$i]");

  # reverse complement $rel_coords_A
  $cur_rel_coords = vdr_CoordsReverseComplement($rel_coords_A[$i], 0, undef);
  $cur_val = vdr_CoordsRelativeToAbsolute($abs_coords_A[$i], 
                                          $cur_rel_coords, undef);
  $cur_val = vdr_CoordsReverseComplement($cur_val, 0, undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsRelativeToAbsolute() revcomp: $desc_A[$i]");

  # test vdr_CoordsRelativeSingleCoordToAbsolute() on revcomp
  # 5' most position
  $cur_rel_coord = vdr_Feature5pMostPosition($cur_rel_coords, undef);
  $cur_exp_val   = vdr_Feature3pMostPosition($exp_val_A[$i], undef);
  $cur_val       = vdr_CoordsRelativeSingleCoordToAbsolute($abs_coords_A[$i], $cur_rel_coord, undef);
  is($cur_val, $cur_exp_val, "vdr_CoordsRelativeSingleCoordToAbsolute() 5'-most position revcomp: $desc_A[$i]");
  # 3' most position
  $cur_rel_coord = vdr_Feature3pMostPosition($cur_rel_coords, undef);
  $cur_exp_val   = vdr_Feature5pMostPosition($exp_val_A[$i], undef);
  $cur_val       = vdr_CoordsRelativeSingleCoordToAbsolute($abs_coords_A[$i], $cur_rel_coord, undef);
  is($cur_val, $cur_exp_val, "vdr_CoordsRelativeSingleCoordToAbsolute() 3'-most position revcomp: $desc_A[$i]");
}
#############################################
# vdr_CoordsProteinRelativeToAbsolute() tests
#
# Relative protein coords segments must all
# be "+" strand.
#############################################
@desc_A          = ();
@abs_coords_A    = ();
@rel_pt_coords_A = ();
@exp_val_A       = ();

push(@desc_A,          "abs (+), rel (+)");
push(@abs_coords_A,    "11..100:+");
push(@rel_pt_coords_A, "6..10:+");    
push(@exp_val_A,       "26..40:+");   

push(@desc_A,          "abs (+), rel (++)");
push(@abs_coords_A,    "11..100:+");
push(@rel_pt_coords_A, "6..10:+,12..13:+");    
push(@exp_val_A,       "26..40:+,44..49:+");   

push(@desc_A,          "abs (++), rel (+)");
push(@abs_coords_A,    "11..30:+,30..100:+");
push(@rel_pt_coords_A, "6..10:+");    
push(@exp_val_A,       "26..30:+,30..39:+");

push(@desc_A,          "abs (++), rel (++)");
push(@abs_coords_A,    "11..30:+,30..100:+");
push(@rel_pt_coords_A, "6..10:+,12..13:+");    
push(@exp_val_A,       "26..30:+,30..39:+,43..48:+");

# The following test is excluded because current code requires all
# absolute coords segments to be the same strand as a contract check.
# existing code should work for it though if you ever want to 
# relax that requirement.
#push(@desc_A,          "abs (+-), rel (++)");
#push(@abs_coords_A,    "11..30:+,100..30:-");
#push(@rel_pt_coords_A, "6..10:+,12..13:+");    
#push(@exp_val_A,       "26..30:+,100..91:-,87..82:-");

$ntests = scalar(@desc_A);
for($i = 0; $i < $ntests; $i++) { 
  $cur_val = vdr_CoordsProteinRelativeToAbsolute($abs_coords_A[$i], 
                                                 $rel_pt_coords_A[$i], undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsProteinRelativeToAbsolute(): $desc_A[$i]");
  # sanity check:
  # make sure the length of the returned coords string is the same as the
  # length of the $rel_nt_or_aa_coords string
  $rel_length     = vdr_CoordsLength($rel_pt_coords_A[$i], undef);
  $cur_val_length = (vdr_CoordsLength($cur_val, undef) / 3);
  is($cur_val_length, $rel_length, "vdr_CoordsProteinRelativeToAbsolute() and vdr_CoordsLength(): length sanity check for $desc_A[$i]");
}

#############################################
# vdr_FrameAdjust() tests
# hard-coded
#############################################
my ($ret_frame1, $ret_frame2, $ret_frame3) = undef;
#                             <orig_frame> <nt_diff>
$ret_frame1 = vdr_FrameAdjust(1,           0,        undef);
$ret_frame2 = vdr_FrameAdjust(1,           1,        undef);
$ret_frame3 = vdr_FrameAdjust(1,           2,        undef);
is($ret_frame1, 1, "vdr_FrameAdjust(): orig_frame:1, nt_diff:0");
is($ret_frame2, 3, "vdr_FrameAdjust(): orig_frame:1, nt_diff:1");
is($ret_frame3, 2, "vdr_FrameAdjust(): orig_frame:1, nt_diff:2");
$ret_frame1 = vdr_FrameAdjust(1,          -1,        undef);
$ret_frame2 = vdr_FrameAdjust(1,          -2,        undef);
is($ret_frame1, 2, "vdr_FrameAdjust(): orig_frame:1, nt_diff:-2");
is($ret_frame2, 3, "vdr_FrameAdjust(): orig_frame:1, nt_diff:-1");

$ret_frame1 = vdr_FrameAdjust(2,           0,        undef);
$ret_frame2 = vdr_FrameAdjust(2,           1,        undef);
$ret_frame3 = vdr_FrameAdjust(2,           2,        undef);
is($ret_frame1, 2, "vdr_FrameAdjust(): orig_frame:2, nt_diff:0");
is($ret_frame2, 1, "vdr_FrameAdjust(): orig_frame:2, nt_diff:1");
is($ret_frame3, 3, "vdr_FrameAdjust(): orig_frame:2, nt_diff:2");
$ret_frame1 = vdr_FrameAdjust(2,          -1,        undef);
$ret_frame2 = vdr_FrameAdjust(2,          -2,        undef);
is($ret_frame1, 3, "vdr_FrameAdjust(): orig_frame:2, nt_diff:-2");
is($ret_frame2, 1, "vdr_FrameAdjust(): orig_frame:2, nt_diff:-1");

$ret_frame1 = vdr_FrameAdjust(3,           0,        undef);
$ret_frame2 = vdr_FrameAdjust(3,           1,        undef);
$ret_frame3 = vdr_FrameAdjust(3,           2,        undef);
is($ret_frame1, 3, "vdr_FrameAdjust(): orig_frame:3, nt_diff:0");
is($ret_frame2, 2, "vdr_FrameAdjust(): orig_frame:3, nt_diff:1");
is($ret_frame3, 1, "vdr_FrameAdjust(): orig_frame:3, nt_diff:2");
$ret_frame1 = vdr_FrameAdjust(3,          -1,        undef);
$ret_frame2 = vdr_FrameAdjust(3,          -2,        undef);
is($ret_frame1, 1, "vdr_FrameAdjust(): orig_frame:3, nt_diff:-2");
is($ret_frame2, 2, "vdr_FrameAdjust(): orig_frame:3, nt_diff:-1");

# a few with abs(nt_diff) above 2  
$ret_frame1 = vdr_FrameAdjust(1,         108,        undef);
$ret_frame2 = vdr_FrameAdjust(1,         109,        undef);
$ret_frame3 = vdr_FrameAdjust(1,         110,        undef);
is($ret_frame1, 1, "vdr_FrameAdjust(): orig_frame:1, nt_diff%3:0");
is($ret_frame2, 3, "vdr_FrameAdjust(): orig_frame:1, nt_diff%3:1");
is($ret_frame3, 2, "vdr_FrameAdjust(): orig_frame:1, nt_diff%3:2");
$ret_frame1 = vdr_FrameAdjust(1,        -109,        undef);
$ret_frame2 = vdr_FrameAdjust(1,        -110,        undef);
is($ret_frame1, 2, "vdr_FrameAdjust(): orig_frame:1, nt_diff%3:-2");
is($ret_frame2, 3, "vdr_FrameAdjust(): orig_frame:1, nt_diff%3:-1");

$ret_frame1 = vdr_FrameAdjust(2,          72,        undef);
$ret_frame2 = vdr_FrameAdjust(2,          73,        undef);
$ret_frame3 = vdr_FrameAdjust(2,          74,        undef);
is($ret_frame1, 2, "vdr_FrameAdjust(): orig_frame:2, nt_diff%3:0");
is($ret_frame2, 1, "vdr_FrameAdjust(): orig_frame:2, nt_diff%3:1");
is($ret_frame3, 3, "vdr_FrameAdjust(): orig_frame:2, nt_diff%3:2");
$ret_frame1 = vdr_FrameAdjust(2,         -73,        undef);
$ret_frame2 = vdr_FrameAdjust(2,         -74,        undef);
is($ret_frame1, 3, "vdr_FrameAdjust(): orig_frame:2, nt_diff%3:-2");
is($ret_frame2, 1, "vdr_FrameAdjust(): orig_frame:2, nt_diff%3:-1");

$ret_frame1 = vdr_FrameAdjust(3,         900,        undef);
$ret_frame2 = vdr_FrameAdjust(3,         901,        undef);
$ret_frame3 = vdr_FrameAdjust(3,         902,        undef);
is($ret_frame1, 3, "vdr_FrameAdjust(): orig_frame:3, nt_diff%3:0");
is($ret_frame2, 2, "vdr_FrameAdjust(): orig_frame:3, nt_diff%3:1");
is($ret_frame3, 1, "vdr_FrameAdjust(): orig_frame:3, nt_diff%3:2");
$ret_frame1 = vdr_FrameAdjust(3,        -901,        undef);
$ret_frame2 = vdr_FrameAdjust(3,        -902,        undef);
is($ret_frame1, 1, "vdr_FrameAdjust(): orig_frame:3, nt_diff%3:-2");
is($ret_frame2, 2, "vdr_FrameAdjust(): orig_frame:3, nt_diff%3:-1");

#######################################################
# vadr_seed.pm:prune_seed_of_terminal_short_segments()
#######################################################
@desc_A                = ();
@exp_val_A             = ();
my @orig_seq_coords_A  = ();  # this will double as mdl coords
my @min_term_sgm_len_A = ();
my @seq_len_A          = ();
my @seq_start_A        = ();
my @seq_stop_A         = ();
my @seq_strand_A       = ();
my @mdl_start_A        = ();
my @mdl_stop_A         = ();
my @mdl_strand_A       = ();
my $seq_coords         = undef;
my $mdl_coords         = undef;

push(@desc_A,             "no change");
push(@orig_seq_coords_A,  "2..102:+,103..205:+,206..400:+");
push(@exp_val_A,          "2..102:+,103..205:+,206..400:+");
push(@min_term_sgm_len_A, "100");
push(@seq_len_A,          "401");

push(@desc_A,             "remove one 5'");
push(@orig_seq_coords_A,  "2..102:+,103..205:+,206..400:+");
push(@exp_val_A,          "103..205:+,206..400:+");
push(@min_term_sgm_len_A, "103");
push(@seq_len_A,          "401");

push(@desc_A,             "remove none'");
push(@orig_seq_coords_A,  "1..102:+,103..205:+,206..400:+");
push(@exp_val_A,          "1..102:+,103..205:+,206..400:+");
push(@min_term_sgm_len_A, "103");
push(@seq_len_A,          "401");

push(@desc_A,             "remove two 5'");
push(@orig_seq_coords_A,  "2..102:+,103..205:+,206..400:+");
push(@exp_val_A,          "206..400:+");
push(@min_term_sgm_len_A, "106");
push(@seq_len_A,          "401");

push(@desc_A,             "remove all");
push(@orig_seq_coords_A,  "2..102:+,103..205:+,206..400:+");
push(@exp_val_A,          "");
push(@min_term_sgm_len_A, "300");
push(@seq_len_A,          "401");

push(@desc_A,             "remove all but last");
push(@orig_seq_coords_A,  "2..102:+,103..205:+,206..400:+");
push(@exp_val_A,          "");
push(@min_term_sgm_len_A, "300");
push(@seq_len_A,          "400");

$ntests = scalar(@desc_A);
for($i = 0; $i < $ntests; $i++) { 
  vdr_FeatureStartStopStrandArrays($orig_seq_coords_A[$i], \@seq_start_A, \@seq_stop_A, \@seq_strand_A, undef);
  vdr_FeatureStartStopStrandArrays($orig_seq_coords_A[$i], \@mdl_start_A, \@mdl_stop_A, \@mdl_strand_A, undef);
  prune_seed_of_terminal_short_segments(\@seq_start_A, \@seq_stop_A, \@seq_strand_A,
                                        \@mdl_start_A, \@mdl_stop_A, \@mdl_strand_A,
                                        $min_term_sgm_len_A[$i],
                                        $seq_len_A[$i], $seq_len_A[$i]);
  $cur_val = vdr_CoordsFromStartStopStrandArrays(\@seq_start_A, \@seq_stop_A, \@seq_strand_A, undef);
  is($cur_val, $exp_val_A[$i], "prune_seed_of_terminal_short_segments(): $desc_A[$i]");
}

########################################################
# vadr_seed.pm:prune_seed_given_minimum_length_segment()
########################################################
@desc_A            = ();
@exp_val_A         = ();
@orig_seq_coords_A = (); # this will double as mdl coords
my @min_sgm_len_A  = ();

push(@desc_A,             "no change");
push(@orig_seq_coords_A,  "1..10:+,11..21:+,22..33:+");
push(@exp_val_A,          "1..10:+,11..21:+,22..33:+");
push(@min_sgm_len_A,      "10");

push(@desc_A,             "remove one 5'");
push(@orig_seq_coords_A,  "1..10:+,11..21:+,22..33:+");
push(@exp_val_A,          "11..21:+,22..33:+");
push(@min_sgm_len_A,      "11");

push(@desc_A,             "remove two 5'");
push(@orig_seq_coords_A,  "1..10:+,11..21:+,22..33:+");
push(@exp_val_A,          "22..33:+");
push(@min_sgm_len_A,      "12");

push(@desc_A,             "keep longest");
push(@orig_seq_coords_A,  "1..10:+,11..21:+,22..33:+");
push(@exp_val_A,          "22..33:+");
push(@min_sgm_len_A,      "50");

push(@desc_A,             "no change (long)");
push(@orig_seq_coords_A,  "1..10:+,11..110:+,111..121:+,122..222:+,223..233:+,233..334:+,335..435:+,434..446:+,444..543:+");
push(@exp_val_A,          "1..10:+,11..110:+,111..121:+,122..222:+,223..233:+,233..334:+,335..435:+,434..446:+,444..543:+");
push(@min_sgm_len_A,      "10");

push(@desc_A,             "remove one 5' (long)");
push(@orig_seq_coords_A,  "1..10:+,11..110:+,111..121:+,122..222:+,223..233:+,233..334:+,335..435:+,434..446:+,444..543:+");
push(@exp_val_A,          "11..110:+,111..121:+,122..222:+,223..233:+,233..334:+,335..435:+,434..446:+,444..543:+");
push(@min_sgm_len_A,      "11");

push(@desc_A,             "remove one 5' (long)");
push(@orig_seq_coords_A,  "1..10:+,11..110:+,111..121:+,122..222:+,223..233:+,233..334:+,335..435:+,434..446:+,444..543:+");
push(@exp_val_A,          "233..334:+,335..435:+,434..446:+,444..543:+");
push(@min_sgm_len_A,      "12");

push(@desc_A,             "remove all but two (long)");
push(@orig_seq_coords_A,  "1..10:+,11..110:+,111..121:+,122..222:+,223..233:+,233..334:+,335..435:+,434..446:+,444..543:+");
push(@exp_val_A,          "233..334:+,335..435:+");
push(@min_sgm_len_A,      "14");

push(@desc_A,             "remove all but one (long)");
push(@orig_seq_coords_A,  "1..10:+,11..110:+,111..121:+,122..222:+,223..233:+,233..334:+,335..435:+,434..446:+,444..543:+");
push(@exp_val_A,          "233..334:+");
push(@min_sgm_len_A,      "102");

push(@desc_A,             "remove all, keep longest (long)");
push(@orig_seq_coords_A,  "1..10:+,11..110:+,111..121:+,122..222:+,223..233:+,233..334:+,335..435:+,434..446:+,444..543:+");
push(@exp_val_A,          "233..334:+");
push(@min_sgm_len_A,      "200");

$ntests = scalar(@desc_A);
for($i = 0; $i < $ntests; $i++) { 
  vdr_FeatureStartStopStrandArrays($orig_seq_coords_A[$i], \@seq_start_A, \@seq_stop_A, \@seq_strand_A, undef);
  vdr_FeatureStartStopStrandArrays($orig_seq_coords_A[$i], \@mdl_start_A, \@mdl_stop_A, \@mdl_strand_A, undef);
  prune_seed_given_minimum_length_segment(\@seq_start_A, \@seq_stop_A, \@seq_strand_A,
                                          \@mdl_start_A, \@mdl_stop_A, \@mdl_strand_A,
                                          $min_sgm_len_A[$i]);
  $cur_val = vdr_CoordsFromStartStopStrandArrays(\@seq_start_A, \@seq_stop_A, \@seq_strand_A, undef);
  is($cur_val, $exp_val_A[$i], "prune_seed_given_minimum_length_segment(): $desc_A[$i]");
}

#############################################
# vdr_CoordsSegmentActualToFractional() and 
# vdr_CoordsSegmentFractionalToActual() tests
#############################################
@desc_A          = ();
@full_coords_A   = ();
@subseq_coords_A = ();
@exp_val_A       = ();

push(@desc_A,          "full (+), subseq (+) 1");
push(@full_coords_A,   "1..100:+");
push(@subseq_coords_A, "11..90:+");    
push(@exp_val_A,       "0.110/0.900");

push(@desc_A,          "full (+), subseq (+) 2");
push(@full_coords_A,   "11..90:+");
push(@subseq_coords_A, "20..60:+");    
push(@exp_val_A,       "0.125/0.625");

push(@desc_A,          "full (+), subseq (+) 3");
push(@full_coords_A,   "11..90:+");
push(@subseq_coords_A, "20..91:+");    
push(@exp_val_A,       "undef");

push(@desc_A,          "full (+), subseq (+) 4");
push(@full_coords_A,   "11..90:+");
push(@subseq_coords_A, "10..90:+");    
push(@exp_val_A,       "undef");

push(@desc_A,          "full (+), subseq (+) 5");
push(@full_coords_A,   "11..110:+");
push(@subseq_coords_A, "11..110:+");    
push(@exp_val_A,       "0.010/1.000");

push(@desc_A,          "full (-), subseq (-) 1");
push(@full_coords_A,   "100..1:-");
push(@subseq_coords_A, "90..11:-");    
push(@exp_val_A,       "0.110/0.900");

push(@desc_A,          "full (-), subseq (-) 2");
push(@full_coords_A,   "90..11:-");
push(@subseq_coords_A, "81..41:-");    
push(@exp_val_A,       "0.125/0.625");

push(@desc_A,          "full (-), subseq (-) 3");
push(@full_coords_A,   "90..11:-");
push(@subseq_coords_A, "91..20:-");    
push(@exp_val_A,       "undef");

push(@desc_A,          "full (-), subseq (-) 4");
push(@full_coords_A,   "90..11:-");
push(@subseq_coords_A, "90..10:-");    
push(@exp_val_A,       "undef");

push(@desc_A,          "full (-), subseq (-) 5");
push(@full_coords_A,   "110..11:-");
push(@subseq_coords_A, "110..11:-");    
push(@exp_val_A,       "0.010/1.000");

push(@desc_A,          "full (+), subseq (-) 1");
push(@full_coords_A,   "1..100:+");
push(@subseq_coords_A, "90..11:-");    
push(@exp_val_A,       "undef");

push(@desc_A,          "full (-), subseq (+) 1");
push(@full_coords_A,   "100..1:-");
push(@subseq_coords_A, "11..90:+");    
push(@exp_val_A,       "undef");

$ntests = scalar(@desc_A);
for($i = 0; $i < $ntests; $i++) { 
  # vdr_CoordsSegmentActualToFractional
  ($fract_start, $fract_stop) = vdr_CoordsSegmentActualToFractional($full_coords_A[$i], $subseq_coords_A[$i], undef);
  if((defined $fract_start) && (defined $fract_stop)) { 
    $cur_val = sprintf("%.3f/%.3f", $fract_start, $fract_stop);
  }
  else { 
    $cur_val = "undef";
  }
  is($cur_val, $exp_val_A[$i], "vdr_CoordsSegmentActualToFractional(): $desc_A[$i]");

  if($cur_val ne "undef") { 
    # vdr_CoordsSegmentFractionalToActual
    ($ret_start, $ret_stop, $ret_strand) = vdr_CoordsSegmentFractionalToActual($full_coords_A[$i], $fract_start, $fract_stop, undef);
    $cur_val  = vdr_CoordsSegmentCreate($ret_start, $ret_stop, $ret_strand, undef);
    is($cur_val, $subseq_coords_A[$i], "vdr_CoordsSegmentFractionalToActual(): $desc_A[$i]");
  }
}

###########################################
# vdr_ReplacePseudoCoordsStringCreate() tests
# vdr_ReplacePseudoCoordsStringParse() tests
##########################################
my @scoords_A   = ("1..15:+", "30..50:+", "70..1000:+", "11..32:+", "11..11:+");
my @mcoords_A   = ("1..15:+", "30..53:+", "72..1000:+", "10..32:+", "10..10:+");
my @diff_A      = ("0",     "-3!",    "2!",      "-1!",     "0");
my @ncount_A    = ("10/15", "21/21",  "431/931",  "0/22",   "1/1");
my @ecount_A    = ("3/5",   "0/0",    "?/?",      "?/?",    "0/0");
my @flush_A     = ("5'",    "3'",     "-",         "3'",    "-");
my @replaced_A  = ("Y",     "N",      "Y",         "N",     "Y");

my @countn_A    = ("10", "21", "431", "0",   "1");
my @nmatch_A    = ("3",  "0",  undef, undef,  "0");
my @nmismatch_A = ("2",  "0",  undef, undef,  "0");

my @exp_pcoords_A    = ();
push(@exp_pcoords_A, "[S:1..15,M:1..15,D:0,N:10/15,E:3/5,F:5',R:Y];");
push(@exp_pcoords_A, "[S:30..50,M:30..53,D:-3!,N:21/21,E:0/0,F:3',R:N];");
push(@exp_pcoords_A, "[S:70..1000,M:72..1000,D:2!,N:431/931,E:?/?,F:-,R:Y];");
push(@exp_pcoords_A, "[S:11..32,M:10..32,D:-1!,N:0/22,E:?/?,F:3',R:N];");
push(@exp_pcoords_A, "[S:11..11,M:10..10,D:0,N:1/1,E:0/0,F:-,R:Y];");

my @cur_scoords_A  = ();
my @cur_mcoords_A  = ();
my @cur_diff_A     = ();
my @cur_ncount_A   = ();
my @cur_ecount_A   = ();
my @cur_flush_A    = ();
my @cur_replaced_A = ();

$ntests = scalar(@scoords_A);
for($i = 0; $i < $ntests; $i++) { 
  my ($sstart, $sstop, $mstart, $mstop) = (undef, undef, undef, undef);
  if($scoords_A[$i] =~ /^(\d+)\.\.(\d+)/) { 
    ($sstart, $sstop) = ($1, $2);
  }
  if($mcoords_A[$i] =~ /^(\d+)\.\.(\d+)/) { 
    ($mstart, $mstop) = ($1, $2);
  }
  $cur_val = vdr_ReplacePseudoCoordsStringCreate($sstart, $sstop, $mstart, $mstop, 
                                                 $countn_A[$i], $nmatch_A[$i], $nmismatch_A[$i], $flush_A[$i], 
                                                 ($replaced_A[$i] eq "Y") ? 1 : 0);
       
  is($cur_val, $exp_pcoords_A[$i], "vdr_ReplacePseudoCoordsStringCreate(): " . ($i+1) . " of $ntests");

  vdr_ReplacePseudoCoordsStringParse($cur_val, \@cur_scoords_A, \@cur_mcoords_A, \@cur_diff_A, \@cur_ncount_A,
                                     \@cur_ecount_A, \@cur_flush_A, \@cur_replaced_A, undef);

  is($cur_scoords_A[0],  $scoords_A[$i],  "vdr_ReplacePseudoCoordsStringParse(): " . ($i+1) . " of $ntests (scoords)");
  is($cur_mcoords_A[0],  $mcoords_A[$i],  "vdr_ReplacePseudoCoordsStringParse(): " . ($i+1) . " of $ntests (mcoords)");
  is($cur_diff_A[0],     $diff_A[$i],     "vdr_ReplacePseudoCoordsStringParse(): " . ($i+1) . " of $ntests (diff)");
  is($cur_ncount_A[0],   $ncount_A[$i],   "vdr_ReplacePseudoCoordsStringParse(): " . ($i+1) . " of $ntests (ncount)");
  is($cur_ecount_A[0],   $ecount_A[$i],   "vdr_ReplacePseudoCoordsStringParse(): " . ($i+1) . " of $ntests (ecount)");
  is($cur_flush_A[0],    $flush_A[$i],    "vdr_ReplacePseudoCoordsStringParse(): " . ($i+1) . " of $ntests (flush)");
  is($cur_replaced_A[0], $replaced_A[$i], "vdr_ReplacePseudoCoordsStringParse(): " . ($i+1) . " of $ntests (replaced)");
}

# same as above, but concatenate all 5, then parse them
$cur_val = "";
for($i = 0; $i < $ntests; $i++) { 
  my ($sstart, $sstop, $mstart, $mstop) = (undef, undef, undef, undef);
  if($scoords_A[$i] =~ /^(\d+)\.\.(\d+)/) { 
    ($sstart, $sstop) = ($1, $2);
  }
  if($mcoords_A[$i] =~ /^(\d+)\.\.(\d+)/) { 
    ($mstart, $mstop) = ($1, $2);
  }
  $cur_val .= vdr_ReplacePseudoCoordsStringCreate($sstart, $sstop, $mstart, $mstop, 
                                                  $countn_A[$i], $nmatch_A[$i], $nmismatch_A[$i], $flush_A[$i], 
                                                  ($replaced_A[$i] eq "Y") ? 1 : 0);
       
}

my $ntok = vdr_ReplacePseudoCoordsStringParse($cur_val, \@cur_scoords_A, \@cur_mcoords_A, \@cur_diff_A, \@cur_ncount_A,
                                              \@cur_ecount_A, \@cur_flush_A, \@cur_replaced_A, undef);

is($ntok, $ntests, "vdr_ReplacePseudoCoordsStringParse() count is correct");
for($i = 0; $i < $ntok; $i++) { 
  is($cur_scoords_A[$i],  $scoords_A[$i],  "vdr_ReplacePseudoCoordsStringParse(): concat (scoords)");
  is($cur_mcoords_A[$i],  $mcoords_A[$i],  "vdr_ReplacePseudoCoordsStringParse(): concat (mcoords)");
  is($cur_diff_A[$i],     $diff_A[$i],     "vdr_ReplacePseudoCoordsStringParse(): concat (diff)");
  is($cur_ncount_A[$i],   $ncount_A[$i],   "vdr_ReplacePseudoCoordsStringParse(): concat (ncount)");
  is($cur_ecount_A[$i],   $ecount_A[$i],   "vdr_ReplacePseudoCoordsStringParse(): concat (ecount)");
  is($cur_flush_A[$i],    $flush_A[$i],    "vdr_ReplacePseudoCoordsStringParse(): concat (flush)");
  is($cur_replaced_A[$i], $replaced_A[$i], "vdr_ReplacePseudoCoordsStringParse(): concat (replaced)");
}
