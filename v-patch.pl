#!/usr/bin/env perl
# EPN* Thu May  7 2026
#
# v-patch.pl: post-hoc patcher for v-build --profile output directories.
#
# Applies a TSV-described set of in-place modifications to a model's
# .minfo (and where applicable .protein.fa + BLAST db). Lets users add
# CDS alternative_ftr_set members and/or feature/model exception
# qualifiers without re-running v-prep/v-build.
#
# The CM file is never modified — added alts share emission probabilities
# with the primary CDS, matching v-build's autoalt convention.
#
use strict;
use warnings;
use Cwd qw(abs_path);
use Getopt::Long qw(:config no_auto_abbrev);
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;

require "vadr.pm";
require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seq.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

#########################################################
# environment / exec discovery
my $env_vadr_blast_dir = utl_DirEnvVarValid("VADRBLASTDIR");
my $env_vadr_easel_dir = utl_DirEnvVarValid("VADREASELDIR");

my %execs_H = ();
$execs_H{"esl-translate"} = $env_vadr_easel_dir . "/esl-translate";
$execs_H{"makeblastdb"}   = $env_vadr_blast_dir . "/makeblastdb";
utl_ExecHValidate(\%execs_H, undef);

#########################################################
# options
my %opt_HH = ();
my @opt_order_A = ();
my %opt_group_desc_H = ();
my $g = 0;

$opt_group_desc_H{++$g} = "basic options";
#       option         type        default group  requires incompat preamble                                    help
opt_Add("-h",          "boolean",  0,          $g,    undef, undef, undef,                                      "display this help",                                    \%opt_HH, \@opt_order_A);
opt_Add("-f",          "boolean",  0,          $g,    undef, undef, "force; overwrite existing .vpatch.log",    "force overwrite of existing patch artifacts",          \%opt_HH, \@opt_order_A);
opt_Add("-v",          "boolean",  0,          $g,    undef, undef, "verbose; echo commands to stdout",         "be verbose",                                           \%opt_HH, \@opt_order_A);
opt_Add("--dry",       "boolean",  0,          $g,    undef, undef, "dry-run; parse and validate only",         "do not modify any files; just report planned ops",     \%opt_HH, \@opt_order_A);
opt_Add("--keep",      "boolean",  0,          $g,    undef, undef, "keep tmp files",                           "do not remove tmp files",                              \%opt_HH, \@opt_order_A);
opt_Add("--stk",       "string",   undef,      $g,    undef, undef, "alignment STK for new-alt protein source", "use <s> as source alignment for translating added alts", \%opt_HH, \@opt_order_A);

my %GetOptions_H = ();
my $synopsis = "v-patch.pl :: post-hoc minfo / protein-db patcher for v-build --profile output\n";
my $usage    = "Usage: v-patch.pl [-options] <vb-dir> <patch.tsv>\n";
my $options_okay = &GetOptions('h'      => \$GetOptions_H{"-h"},
                               'f'      => \$GetOptions_H{"-f"},
                               'v'      => \$GetOptions_H{"-v"},
                               'dry'    => \$GetOptions_H{"--dry"},
                               'keep'   => \$GetOptions_H{"--keep"},
                               'stk=s'  => \$GetOptions_H{"--stk"});

my $total_seconds = -1 * ofile_SecondsSinceEpoch();
my $date          = scalar localtime();
my $version       = "1.x-vpatch";
my $releasedate   = "May 2026";
my $pkgname       = "VADR";

if((! $options_okay) || ($GetOptions_H{"-h"})) {
  ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, undef);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  exit 0;
}

if(scalar(@ARGV) != 2) {
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do v-patch.pl -h\n\n";
  exit(1);
}

my ($vb_dir, $patch_tsv) = (@ARGV);

opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);
opt_ValidateSet(\%opt_HH, \@opt_order_A);

$vb_dir    = abs_path($vb_dir);
$patch_tsv = (-f $patch_tsv) ? abs_path($patch_tsv) : $patch_tsv;

if(! -d $vb_dir)    { die "ERROR, <vb-dir> '$vb_dir' is not a directory"; }
if(! -s $patch_tsv) { die "ERROR, <patch.tsv> '$patch_tsv' does not exist or is empty"; }

#########################################################
# locate model artifacts in vb_dir
my ($root, $minfo_file, $protein_fa_file, $cm_file) = vpatch_locate_artifacts($vb_dir);
my $log_file  = "$vb_dir/$root.vpatch.log";
my $bak_minfo = "$minfo_file.bak";
my $bak_prot  = "$protein_fa_file.bak";

if((! opt_Get("--dry", \%opt_HH)) && (-s $log_file) && (! opt_Get("-f", \%opt_HH))) {
  die "ERROR, $log_file already exists; use -f to overwrite";
}

#########################################################
# parse patch TSV
my @ops_A = ();
vpatch_parse_tsv($patch_tsv, \@ops_A);

#########################################################
# parse current minfo
my @mdl_info_AH = ();
my %ftr_info_HAH = ();
my %FH_H = ("log" => *STDOUT, "cmd" => *STDOUT);
vdr_ModelInfoFileParse($minfo_file, ["name","length"], ["type","coords"], \@mdl_info_AH, \%ftr_info_HAH, \%FH_H);

printf("# v-patch.pl on vb-dir: %s\n", $vb_dir);
printf("# model root:           %s\n", $root);
printf("# minfo:                %s (%d models, %d features in %s)\n",
       $minfo_file, scalar(@mdl_info_AH),
       scalar(@{$ftr_info_HAH{$mdl_info_AH[0]{"name"}}}),
       $mdl_info_AH[0]{"name"});
printf("# patch tsv:            %s (%d ops)\n", $patch_tsv, scalar(@ops_A));

#########################################################
# dispatch
my $nadd_alt = 0;
my $nadd_exc = 0;
my @applied_ops_A = ();   # for logging
for my $op_HR (@ops_A) {
  if($op_HR->{"op"} eq "add_exc") {
    vpatch_apply_add_exc($op_HR, \@mdl_info_AH, \%ftr_info_HAH);
    $nadd_exc++;
    push(@applied_ops_A, $op_HR);
  }
  elsif($op_HR->{"op"} eq "add_alt") {
    vpatch_apply_add_alt($op_HR, \@mdl_info_AH, \%ftr_info_HAH);
    $nadd_alt++;
    push(@applied_ops_A, $op_HR);
  }
  else {
    die "ERROR, line " . $op_HR->{"line"} . ": unsupported op '" . $op_HR->{"op"} . "'";
  }
}

printf("# add_alt ops:          %d\n", $nadd_alt);
printf("# add_exc ops:          %d\n", $nadd_exc);

if(opt_Get("--dry", \%opt_HH)) {
  print "# DRY RUN — no files written.\n";
  exit 0;
}

#########################################################
# write outputs
vpatch_finalize($vb_dir, $root, $minfo_file, $protein_fa_file,
                $bak_minfo, $bak_prot, $log_file,
                \@mdl_info_AH, \%ftr_info_HAH, \@applied_ops_A);

$total_seconds += ofile_SecondsSinceEpoch();
printf("# done in %.1f seconds.\n", $total_seconds);
exit 0;

#################################################################
# Subroutine: vpatch_locate_artifacts()
# Incept:     EPN* Thu May  7 2026
#
# Given a <vb-dir>, identify <root>.minfo, <root>.protein.fa, <root>.cm.
# Returns ($root, $minfo, $protein_fa, $cm).
#################################################################
sub vpatch_locate_artifacts {
  my ($dir) = @_;
  opendir(my $dh, $dir) or die "ERROR, cannot read $dir: $!";
  my @minfo = grep { /\.vadr\.minfo$/ } readdir($dh);
  closedir($dh);
  if(scalar(@minfo) != 1) {
    die sprintf("ERROR, expected exactly one *.vadr.minfo in %s, found %d", $dir, scalar(@minfo));
  }
  (my $root = $minfo[0]) =~ s/\.minfo$//;
  my $minfo_f = "$dir/$minfo[0]";
  my $prot_f  = "$dir/$root.protein.fa";
  my $cm_f    = "$dir/$root.cm";
  (-s $prot_f) or die "ERROR, missing $prot_f";
  (-s $cm_f)   or die "ERROR, missing $cm_f";
  return ($root, $minfo_f, $prot_f, $cm_f);
}

#################################################################
# Subroutine: vpatch_parse_tsv()
# Incept:     EPN* Thu May  7 2026
#
# Parse a TSV with a header line ("op<TAB>target_model<TAB>target_ftr_coords
# <TAB>key<TAB>value") plus one row per op. Stores each row as a hashref
# in @{$ops_AR}. Whitespace-separated (any whitespace, not strictly tabs).
#################################################################
sub vpatch_parse_tsv {
  my ($tsv, $ops_AR) = @_;
  open(my $fh, "<", $tsv) or die "ERROR, cannot read $tsv: $!";
  my $hdr;
  while(defined($hdr = <$fh>)) {
    chomp $hdr;
    next if $hdr =~ /^\s*$/;
    next if $hdr =~ /^\s*#/;
    last;
  }
  if(! defined $hdr) { die "ERROR, $tsv contains no header"; }
  $hdr =~ s/^\s+//;
  my @hdr_A = split(/\s+/, $hdr);
  if(scalar(@hdr_A) < 5 || $hdr_A[0] ne "op") {
    die "ERROR, $tsv header must start with: op target_model target_ftr_coords key value (whitespace separated)";
  }
  my $line_n = 1;
  while(my $line = <$fh>) {
    $line_n++;
    chomp $line;
    next if $line =~ /^\s*$/;
    next if $line =~ /^\s*#/;
    my @F = split(/\s+/, $line);
    if(scalar(@F) < 5) {
      die "ERROR, $tsv line $line_n: expected >=5 whitespace-separated fields, got " . scalar(@F);
    }
    my ($op, $model, $coords, $key, @rest) = @F;
    my $value = join(" ", @rest);
    push(@{$ops_AR}, {
      "op"     => $op,
      "model"  => $model,
      "coords" => $coords,
      "key"    => $key,
      "value"  => $value,
      "line"   => $line_n,
    });
  }
  close($fh);
}

#################################################################
# Subroutine: vpatch_apply_add_exc()
# Incept:     EPN* Thu May  7 2026
#
# Apply one add_exc op: append a key:value qualifier to either the
# MODEL line (if key is "lowsim_exc") or a single matching CDS FEATURE
# line (matched on coords). If the key already exists, append the new
# value to the existing one separated by ",". The minfo writer in
# vadr.pm preserves the in-memory map exactly, so order of existing
# keys is unchanged.
#################################################################
sub vpatch_apply_add_exc {
  my ($op_HR, $mdl_info_AHR, $ftr_info_HAHR) = @_;
  my $model  = $op_HR->{"model"};
  my $coords = $op_HR->{"coords"};
  my $key    = $op_HR->{"key"};
  my $value  = $op_HR->{"value"};
  my $line_n = $op_HR->{"line"};

  # locate model
  my $mdl_idx = -1;
  for(my $i = 0; $i < scalar(@{$mdl_info_AHR}); $i++) {
    if($mdl_info_AHR->[$i]{"name"} eq $model) { $mdl_idx = $i; last; }
  }
  if($mdl_idx < 0) { die "ERROR, add_exc line $line_n: model '$model' not found in minfo"; }

  # MODEL-level exceptions
  if($key eq "lowsim_exc") {
    if(defined $mdl_info_AHR->[$mdl_idx]{$key} && $mdl_info_AHR->[$mdl_idx]{$key} ne "") {
      $mdl_info_AHR->[$mdl_idx]{$key} .= "," . $value;
    }
    else {
      $mdl_info_AHR->[$mdl_idx]{$key} = $value;
    }
    printf("# applied add_exc MODEL %s key=%s value=%s\n", $model, $key, $value);
    return;
  }

  # FEATURE-level: find a single CDS feature with matching coords.
  my $ftr_AR = $ftr_info_HAHR->{$model};
  if(! defined $ftr_AR) { die "ERROR, add_exc line $line_n: model $model has no features in ftr_info"; }

  my @match_idxs = ();
  for(my $fi = 0; $fi < scalar(@{$ftr_AR}); $fi++) {
    if($ftr_AR->[$fi]{"coords"} eq $coords && $ftr_AR->[$fi]{"type"} eq "CDS") {
      push(@match_idxs, $fi);
    }
  }
  if(scalar(@match_idxs) == 0) {
    die "ERROR, add_exc line $line_n: no CDS feature with coords '$coords' on model '$model'";
  }
  if(scalar(@match_idxs) > 1) {
    die sprintf("ERROR, add_exc line %d: %d CDS features match coords %s on model %s; cannot disambiguate",
                $line_n, scalar(@match_idxs), $coords, $model);
  }
  my $fi = $match_idxs[0];
  if(defined $ftr_AR->[$fi]{$key} && $ftr_AR->[$fi]{$key} ne "") {
    $ftr_AR->[$fi]{$key} .= "," . $value;
  }
  else {
    $ftr_AR->[$fi]{$key} = $value;
  }
  printf("# applied add_exc FEATURE CDS coords=%s key=%s value=%s\n", $coords, $key, $value);
}

#################################################################
# Subroutine: vpatch_apply_add_alt()  -- stub for commit 1
# Incept:     EPN* Thu May  7 2026
#################################################################
sub vpatch_apply_add_alt {
  my ($op_HR, $mdl_info_AHR, $ftr_info_HAHR) = @_;
  printf("# [planned] add_alt %s coords=%s set=%s\n",
         $op_HR->{"model"}, $op_HR->{"coords"}, $op_HR->{"value"});
}

#################################################################
# Subroutine: vpatch_finalize()  -- stub for commit 1
# Incept:     EPN* Thu May  7 2026
#################################################################
sub vpatch_finalize {
  die "ERROR, finalize not yet implemented (will be in commit 4) — use --dry for now";
}
