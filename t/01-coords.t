use strict;
use warnings FATAL => 'all';
use Test::More tests => 287;

BEGIN {
    use_ok( 'vadr' ) || print "Bail out!\n";
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
my ($rel_length, $cur_val_length, $cur_rel_coords);

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

  # reverse complement $rel_coords_A
  $cur_rel_coords = vdr_CoordsReverseComplement($rel_coords_A[$i], 0, undef);
  $cur_val = vdr_CoordsRelativeToAbsolute($abs_coords_A[$i], 
                                          $cur_rel_coords, undef);
  $cur_val = vdr_CoordsReverseComplement($cur_val, 0, undef);
  is($cur_val, $exp_val_A[$i], "vdr_CoordsRelativeToAbsolute() revcomp: $desc_A[$i]");
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

push(@desc_A,          "abs (+-), rel (++)");
push(@abs_coords_A,    "11..30:+,100..30:-");
push(@rel_pt_coords_A, "6..10:+,12..13:+");    
push(@exp_val_A,       "26..30:+,100..91:-,87..82:-");

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

