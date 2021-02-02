use strict;
use warnings FATAL => 'all';
use Test::More tests => 36;

BEGIN {
    use_ok( 'vadr' ) || print "Bail out!\n";
}

# make sure the VADRSCRIPTSDIR env variable is set
my $env_ok = exists $ENV{"VADRSCRIPTSDIR"} ? 1 : 0;
is($env_ok, 1, "VADRSCRIPTSDIR env variable set");

# make sure we have the required files
my @reqd_file_A = ();
my $minfo_file = $ENV{"VADRSCRIPTSDIR"} . "/t/data/entoy100a.minfo";
push(@reqd_file_A, $minfo_file);

foreach my $reqd_file (@reqd_file_A) {
  my $file_nonempty = (-s $reqd_file) ? 1 : 0;
  is($file_nonempty, 1, "Required file $reqd_file exists and is nonempty");
}

###########################
# Parse the model info file
###########################
my @mdl_info_AH  = (); # array of hashes with model info
my %ftr_info_HAH = (); # hash of array of hashes with feature info 
my %sgm_info_HAH = (); # hash of array of hashes with segment info 

my @reqd_mdl_keys_A = ("name", "length");
my @reqd_ftr_keys_A = ("type", "coords");
vdr_ModelInfoFileParse($minfo_file, \@reqd_mdl_keys_A, \@reqd_ftr_keys_A, \@mdl_info_AH, \%ftr_info_HAH, undef);

# validate %mdl_info_AH
my @mdl_reqd_keys_A = ("name", "length");
my $nmdl = utl_AHValidate(\@mdl_info_AH, \@mdl_reqd_keys_A, "ERROR reading model info from $minfo_file", undef);

is($nmdl, 1, "correct number of models read");
my $mdl_name = $mdl_info_AH[0]{"name"};

vdr_FeatureInfoImputeLength(\@{$ftr_info_HAH{$mdl_name}}, undef);
vdr_FeatureInfoInitializeParentIndexStrings(\@{$ftr_info_HAH{$mdl_name}}, undef);
vdr_FeatureInfoValidateParentIndexStrings(\@{$ftr_info_HAH{$mdl_name}}, undef);
vdr_FeatureInfoImpute3paFtrIdx(\@{$ftr_info_HAH{$mdl_name}}, undef);
vdr_FeatureInfoImputeOutname(\@{$ftr_info_HAH{$mdl_name}});
vdr_SegmentInfoPopulate(\@{$sgm_info_HAH{$mdl_name}}, \@{$ftr_info_HAH{$mdl_name}}, undef);

my $nftr = scalar(@{$ftr_info_HAH{$mdl_name}});
is($nftr, 11, "correct number of features read");

my $nsgm = scalar(@{$sgm_info_HAH{$mdl_name}});
is($nsgm, 15, "correct number of segments created");

my @segment_start_identical_to_cds_A = ();
my @segment_stop_identical_to_cds_A  = ();

for(my $sgm_idx = 0; $sgm_idx < $nsgm; $sgm_idx++) {
  $segment_start_identical_to_cds_A[$sgm_idx] = vdr_SegmentStartIdenticalToCds(\@{$ftr_info_HAH{$mdl_name}}, \@{$sgm_info_HAH{$mdl_name}}, $sgm_idx, undef);
  $segment_stop_identical_to_cds_A[$sgm_idx]  = vdr_SegmentStopIdenticalToCds (\@{$ftr_info_HAH{$mdl_name}}, \@{$sgm_info_HAH{$mdl_name}}, $sgm_idx, undef);
}

#MODEL ENTOY100A cmfile:"entoy100a.cm" group:"toy" length:"100" subgroup:"A"
#FEATURE ENTOY100A type:"gene" coords:"8..33:+" parent_idx_str:"GBNULL" gene:"one"
#FEATURE ENTOY100A type:"CDS" coords:"11..31:+" parent_idx_str:"GBNULL" gene:"one" product:"protein one"
#FEATURE ENTOY100A type:"mat_peptide" coords:"11..22:+" parent_idx_str:"1" product:"protein one mp1"
#FEATURE ENTOY100A type:"mat_peptide" coords:"23..28:+" parent_idx_str:"1" product:"protein one mp2"
#FEATURE ENTOY100A type:"CDS" coords:"32..47:+,47..64:+,66..73:+" parent_idx_str:"GBNULL" gene:"two" product:"protein two"
#FEATURE ENTOY100A type:"CDS" coords:"67..79:+,85..95:+" parent_idx_str:"GBNULL" gene:"three" product:"protein three"
#FEATURE ENTOY100A type:"mat_peptide" coords:"67..75:+" parent_idx_str:"5" product:"protein three mp1"
#FEATURE ENTOY100A type:"mat_peptide" coords:"76..79:+,85..92:+" parent_idx_str:"5" product:"protein three mp2"
#FEATURE ENTOY100A type:"CDS" coords:"98..84:-" parent_idx_str:"GBNULL" gene:"four" product:"protein four"
#FEATURE ENTOY100A type:"ncRNA" coords:"5..22:+" parent_idx_str:"GBNULL" note:"ncRNA 1" product:"enRNA" ncRNA_class:"lncRNA"
#FEATURE ENTOY100A type:"stem_loop" coords:"10..19:+" parent_idx_str:"GBNULL" note:"stem-loop 1"

is($segment_start_identical_to_cds_A[0],  0,  "segment 0 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[1],  1,  "segment 1 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[2],  1,  "segment 2 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[3],  0,  "segment 3 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[4],  1,  "segment 4 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[5],  0,  "segment 5 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[6],  0,  "segment 6 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[7],  1,  "segment 7 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[8],  0,  "segment 8 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[9],  1,  "segment 9 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[10], 0, "segment 10 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[11], 0, "segment 11 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[12], 1, "segment 12 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[13], 0, "segment 13 start identical to cds calc'ed correctly");
is($segment_start_identical_to_cds_A[14], 0, "segment 14 start identical to cds calc'ed correctly");

is($segment_stop_identical_to_cds_A[0],  0,  "segment 0 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[1],  1,  "segment 1 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[2],  0,  "segment 2 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[3],  0,  "segment 3 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[4],  0,  "segment 4 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[5],  0,  "segment 5 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[6],  1,  "segment 6 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[7],  0,  "segment 7 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[8],  1,  "segment 8 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[9],  0,  "segment 9 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[10], 0, "segment 10 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[11], 0, "segment 11 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[12], 1, "segment 12 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[13], 0, "segment 13 stop identical to cds calc'ed correctly");
is($segment_stop_identical_to_cds_A[14], 0, "segment 14 stop identical to cds calc'ed correctly");
