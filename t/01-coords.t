use strict;
use warnings FATAL => 'all';
use Test::More tests => 11;

BEGIN {
    use_ok( 'vadr' ) || print "Bail out!\n";
}

my @desc_A                = ();
my @full_abs_nt_coords_A  = ();
my @rel_nt_or_aa_coords_A = ();
my @rel_is_aa_A           = ();
my @ret_val_A             = ();

# positive strand, 1 segment
push(@desc_A,                "+ strand, single segment, relative coords are nt coords");
push(@full_abs_nt_coords_A,  "11..100:+");
push(@rel_nt_or_aa_coords_A, "6..38:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "16..48:+");   

# positive strand, 2 segments
push(@desc_A,                "+ strand, two segments, only first spanned");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "6..25:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "16..35:+");   

push(@desc_A,                "+ strand, two segments, only second spanned");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "31..33:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "42..44:+");   

push(@desc_A,                "+ strand, two segments, both spanned");
push(@full_abs_nt_coords_A,  "11..40:+,42..101:+");
push(@rel_nt_or_aa_coords_A, "6..38:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "16..40:+,42..49:+");   

# positive strand, 3 segments
push(@desc_A,                "+ strand, three segments, only first spanned");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "6..25:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "16..35:+");   

push(@desc_A,                "+ strand, three segments, only first spanned");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "6..25:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "16..35:+");   

push(@desc_A,                "+ strand, three segments, only second spanned");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "31..33:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "42..44:+");   

push(@desc_A,                "+ strand, three segments, only third spanned");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "90..105:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "105..120:+");   

push(@desc_A,                "+ strand, three segments, first two spanned");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "6..38:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "16..40:+,42..49:+");   

push(@desc_A,                "+ strand, three segments, last two spanned");
push(@full_abs_nt_coords_A,  "11..40:+,42..100:+,105..121:+");
push(@rel_nt_or_aa_coords_A, "39..99:+");    
push(@rel_is_aa_A,           "0");          # nt coords
push(@ret_val_A,             "50..100:+,105..114:+");   

my $n = scalar(@full_abs_nt_coords_A);
for(my $i = 0; $i < $n; $i++) { 
  my $cur_val = vdr_CoordsRelativeToAbsolute($full_abs_nt_coords_A[$i], 
                                             $rel_nt_or_aa_coords_A[$i], 
                                             $rel_is_aa_A[$i], undef);
  is($cur_val, $ret_val_A[$i], "vdr_CoordsRelativeToAbsolute(): $desc_A[$i]");
}
