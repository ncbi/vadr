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
opt_Add("--altprotein", "string",  undef,      $g,    undef, undef, "pre-translated protein fasta <s> for added alts", "fasta <s> of pre-built alt proteins (headers must end with /refcoords)", \%opt_HH, \@opt_order_A);

my %GetOptions_H = ();
my $synopsis = "v-patch.pl :: post-hoc minfo / protein-db patcher for v-build --profile output\n";
my $usage    = "Usage: v-patch.pl [-options] <vb-dir> <patch.tsv>\n";
my $options_okay = &GetOptions('h'            => \$GetOptions_H{"-h"},
                               'f'            => \$GetOptions_H{"-f"},
                               'v'            => \$GetOptions_H{"-v"},
                               'dry'          => \$GetOptions_H{"--dry"},
                               'keep'         => \$GetOptions_H{"--keep"},
                               'stk=s'        => \$GetOptions_H{"--stk"},
                               'altprotein=s' => \$GetOptions_H{"--altprotein"});

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
      $mdl_info_AHR->[$mdl_idx]{$key} = vpatch_merge_exc_value($key, $mdl_info_AHR->[$mdl_idx]{$key}, $value, $line_n);
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
    $ftr_AR->[$fi]{$key} = vpatch_merge_exc_value($key, $ftr_AR->[$fi]{$key}, $value, $line_n);
  }
  else {
    $ftr_AR->[$fi]{$key} = $value;
  }
  printf("# applied add_exc FEATURE CDS coords=%s key=%s value=%s\n", $coords, $key, $value);
}

#################################################################
# Subroutine: vpatch_merge_exc_value()
# Incept:     EPN* Tue May 12 2026
#
# Merge a new exception-coords value into an existing one for the same
# key, detecting and resolving overlaps between segments. Each segment
# is "<start>..<stop>:<strand>:<value>" (coords-value) or
# "<start>..<stop>:<strand>" (coords-only). v-annotate's validator
# (vdr_ExceptionCoordsAndValuesToSegmentsAndValues) rejects overlapping
# segments, so blind concatenation can produce an invalid minfo when
# v-build/v-prep has already emitted a segment that the patch TSV
# overlaps.
#
# Resolution: same-strand overlap with identical value collapses to the
# larger-containing span; partial overlap (neither contains) or
# conflicting values dies with a clear error.
#################################################################
sub vpatch_merge_exc_value {
  my ($key, $existing, $new, $line_n) = @_;
  my @exist_segs = split(/,/, $existing);
  for my $nseg (split(/,/, $new)) {
    my ($ns, $ne, $nstr, $nv) = vpatch_parse_exc_segment($nseg, $line_n);
    my ($nlo, $nhi) = ($ns <= $ne) ? ($ns, $ne) : ($ne, $ns);
    my $skip_new = 0;
    my @replace_idxs = ();
    for(my $i = 0; $i < scalar(@exist_segs); $i++) {
      my ($es, $ee, $estr, $ev) = vpatch_parse_exc_segment($exist_segs[$i], $line_n);
      next if $estr ne $nstr;
      my ($elo, $ehi) = ($es <= $ee) ? ($es, $ee) : ($ee, $es);
      next if($nhi < $elo || $ehi < $nlo); # no overlap
      if(defined $nv && defined $ev && $nv ne $ev) {
        die "ERROR, add_exc line $line_n: $key overlap with conflicting values: existing '$exist_segs[$i]' vs new '$nseg'";
      }
      my $new_contains = ($nlo <= $elo && $nhi >= $ehi);
      my $exist_contains = ($elo <= $nlo && $ehi >= $nhi);
      if($exist_contains) {
        printf("# add_exc line %d: dropping redundant %s span %s (contained in existing %s)\n", $line_n, $key, $nseg, $exist_segs[$i]);
        $skip_new = 1;
        last;
      }
      elsif($new_contains) {
        push(@replace_idxs, $i);
      }
      else {
        die "ERROR, add_exc line $line_n: $key partial overlap (neither contains the other): existing '$exist_segs[$i]' vs new '$nseg'";
      }
    }
    next if $skip_new;
    for my $idx (sort {$b<=>$a} @replace_idxs) {
      printf("# add_exc line %d: replacing %s segment %s with broader %s\n", $line_n, $key, $exist_segs[$idx], $nseg);
      splice(@exist_segs, $idx, 1);
    }
    push(@exist_segs, $nseg);
  }
  return join(",", @exist_segs);
}

#################################################################
# Subroutine: vpatch_parse_exc_segment()
# Incept:     EPN* Tue May 12 2026
#
# Parse one exception-coords segment of the form
# "<start>..<stop>:<strand>:<value>" (coords-value) or
# "<start>..<stop>:<strand>" (coords-only). Tolerates leading "<" or
# ">" on the start/stop tokens. Returns (start, stop, strand, value-or-undef).
#################################################################
sub vpatch_parse_exc_segment {
  my ($seg, $line_n) = @_;
  if($seg =~ /^\<?(\d+)\.\.\>?(\d+):([+\-]):(\S+)$/) {
    return ($1, $2, $3, $4);
  }
  if($seg =~ /^\<?(\d+)\.\.\>?(\d+):([+\-])$/) {
    return ($1, $2, $3, undef);
  }
  die "ERROR, add_exc line $line_n: cannot parse exception segment '$seg'";
}

#################################################################
# Subroutine: vpatch_apply_add_alt()
# Incept:     EPN* Thu May  7 2026
#
# Apply one add_alt op: locate the alternative_ftr_set's primary CDS,
# determine the matching gene-set name and the next subn index, append
# a new CDS entry and a new gene entry to ftr_info. Protein extraction
# is deferred to finalize (see vpatch_extract_alt_protein()).
#
# Only key supported: "alternative_ftr_set" (value is the set name).
# Refuses to add a CDS whose coords already exist in the same set.
# Inherits gene/product/parent_idx_str from the primary; does NOT
# inherit exceptions (those are alt-specific and added via add_exc).
#################################################################
sub vpatch_apply_add_alt {
  my ($op_HR, $mdl_info_AHR, $ftr_info_HAHR) = @_;
  my $model      = $op_HR->{"model"};
  my $new_coords = $op_HR->{"coords"};
  my $key        = $op_HR->{"key"};
  my $set_name   = $op_HR->{"value"};
  my $line_n     = $op_HR->{"line"};

  if($key ne "alternative_ftr_set") {
    die "ERROR, add_alt line $line_n: only key 'alternative_ftr_set' is supported (got '$key')";
  }
  my $ftr_AR = $ftr_info_HAHR->{$model};
  if(! defined $ftr_AR) { die "ERROR, add_alt line $line_n: no features for model $model"; }

  # locate existing members of the CDS alt set
  my @set_cds_idxs = ();
  for(my $fi = 0; $fi < scalar(@{$ftr_AR}); $fi++) {
    if($ftr_AR->[$fi]{"type"} eq "CDS"
       && defined $ftr_AR->[$fi]{"alternative_ftr_set"}
       && $ftr_AR->[$fi]{"alternative_ftr_set"} eq $set_name) {
      push(@set_cds_idxs, $fi);
    }
  }
  if(scalar(@set_cds_idxs) == 0) {
    die "ERROR, add_alt line $line_n: alternative_ftr_set '$set_name' has no existing CDS members on model $model";
  }
  # refuse duplicate
  for my $fi (@set_cds_idxs) {
    if($ftr_AR->[$fi]{"coords"} eq $new_coords) {
      die "ERROR, add_alt line $line_n: CDS with coords '$new_coords' already exists in set '$set_name' (ftr_idx=$fi)";
    }
  }
  my $primary_cds_HR = $ftr_AR->[$set_cds_idxs[0]];

  # locate corresponding gene set (gene with coords matching any CDS member that has alternative_ftr_set)
  my %cds_coords_in_set = ();
  for my $fi (@set_cds_idxs) { $cds_coords_in_set{$ftr_AR->[$fi]{"coords"}} = 1; }
  my $gene_set_name = undef;
  my $primary_gene_HR = undef;
  for(my $fi = 0; $fi < scalar(@{$ftr_AR}); $fi++) {
    next unless $ftr_AR->[$fi]{"type"} eq "gene";
    next unless exists $cds_coords_in_set{$ftr_AR->[$fi]{"coords"}};
    if(defined $ftr_AR->[$fi]{"alternative_ftr_set"} && $ftr_AR->[$fi]{"alternative_ftr_set"} ne "") {
      $gene_set_name //= $ftr_AR->[$fi]{"alternative_ftr_set"};
      if($ftr_AR->[$fi]{"coords"} eq $primary_cds_HR->{"coords"}) {
        $primary_gene_HR = $ftr_AR->[$fi];
      }
    }
  }

  # the new CDS will be appended at index (current size); its position in
  # the alt set becomes (count of existing CDS members + 1)
  my $new_subn_idx = scalar(@set_cds_idxs) + 1;

  # build new CDS feature: shallow copy of primary, swap coords + alt_set,
  # drop primary-specific exception qualifiers
  my %new_cds = ();
  for my $k (keys %{$primary_cds_HR}) { $new_cds{$k} = $primary_cds_HR->{$k}; }
  $new_cds{"coords"} = $new_coords;
  $new_cds{"alternative_ftr_set"} = $set_name;
  delete $new_cds{"alternative_ftr_set_subn"};
  for my $exc_k (qw(deletin_exc insertn_exc lowsim_exc cdsstopn_exc cdsstopp_exc mutendcd_exc indf3pst_exc indf5pst_exc)) {
    delete $new_cds{$exc_k};
  }

  # build new gene feature (if a primary gene was found)
  my $new_gene_HR = undef;
  if(defined $primary_gene_HR && defined $gene_set_name) {
    my %new_gene = ();
    for my $k (keys %{$primary_gene_HR}) { $new_gene{$k} = $primary_gene_HR->{$k}; }
    $new_gene{"coords"} = $new_coords;
    $new_gene{"alternative_ftr_set"} = $gene_set_name;
    $new_gene{"alternative_ftr_set_subn"} = $set_name . "." . $new_subn_idx;
    $new_gene_HR = \%new_gene;
  }

  # append gene first (matches v-build's gene-before-CDS ordering for alts)
  if(defined $new_gene_HR) { push(@{$ftr_AR}, $new_gene_HR); }
  push(@{$ftr_AR}, \%new_cds);

  # remember what we appended; finalize will extract a protein for this CDS
  $op_HR->{"_new_cds_coords"} = $new_coords;
  $op_HR->{"_new_cds_gene"}   = $primary_cds_HR->{"gene"};
  $op_HR->{"_set_subn"}       = $set_name . "." . $new_subn_idx;

  printf("# applied add_alt CDS+gene coords=%s set=%s subn=%s\n",
         $new_coords, $set_name, $op_HR->{"_set_subn"});
}

#################################################################
# Subroutine: vpatch_load_altprotein_fasta()
# Incept:     EPN* Sun May 11 2026
#
# Parse a fasta file of pre-translated alt proteins. Each header is
# expected to end with "/<refcoords>" (matching the convention v-build
# uses, e.g. ">MZ515594.1:4644..5576:+/4688..5620:+"). Stored as:
#   $by_coords{$refcoords} = [ $seedname, $protein_aa_string ]
# where $seedname is the part of the header BEFORE the final "/".
#
# If multiple entries map to the same coords, the first wins.
#################################################################
sub vpatch_load_altprotein_fasta {
  my ($file, $by_coords_HR) = @_;
  open(my $fh, "<", $file) or die "ERROR, cannot read --altprotein $file: $!";
  my $cur_hdr = undef;
  my $cur_seq = "";
  while(my $line = <$fh>) {
    chomp $line;
    if($line =~ /^>(.+)$/) {
      _altprot_store($cur_hdr, $cur_seq, $by_coords_HR) if defined $cur_hdr;
      $cur_hdr = $1;
      $cur_seq = "";
    } elsif(defined $cur_hdr) {
      $line =~ s/\s+//g;
      $cur_seq .= $line;
    }
  }
  _altprot_store($cur_hdr, $cur_seq, $by_coords_HR) if defined $cur_hdr;
  close($fh);
}

sub _altprot_store {
  my ($hdr, $seq, $by_coords_HR) = @_;
  return unless $hdr =~ /^(\S+)\/(\d+\.\.\d+:[\+\-])$/;
  my ($seedfull, $coords) = ($1, $2);
  return if exists $by_coords_HR->{$coords};
  $by_coords_HR->{$coords} = [ $seedfull, $seq ];
}

#################################################################
# Subroutine: vpatch_extract_alt_protein()
# Incept:     EPN* Thu May  7 2026
#
# Given a Bio::Easel::MSA object (with RF annotation) and reference
# coords for a new alt CDS, find the first seed sequence whose
# ungapped residues at those RF positions form a complete in-frame
# ORF (matches profile_ValidateCdsIsComplete semantics), translate
# it, and return:
#   ($seedname, $protein_string)
# or (undef, undef) if no seed validates.
#
# Only single-segment, single-strand coords are supported here (which
# is sufficient for CDS alts created from autoalt-style stop-extension).
#################################################################
sub vpatch_extract_alt_protein {
  my ($msa, $coords, $FH_HR) = @_;
  # parse coords like "4688..5641:+"
  if($coords !~ /^(\d+)\.\.(\d+)\:([\+\-])$/) {
    die "ERROR, vpatch_extract_alt_protein: unsupported coords '$coords' (single-segment +/- only)";
  }
  my ($rfstart, $rfstop, $strand) = ($1, $2, $3);
  my $astart = $msa->rfpos_to_aligned_pos($rfstart);
  my $astop  = $msa->rfpos_to_aligned_pos($rfstop);
  if($astart > $astop) { ($astart, $astop) = ($astop, $astart); }

  my $nseq = $msa->nseq;
  for(my $i = 0; $i < $nseq; $i++) {
    my $sqname = $msa->get_sqname($i);
    my $sub_str = $msa->get_sqstring_unaligned_and_truncated($i, $astart, $astop);
    next if (! defined $sub_str) || $sub_str eq "";
    if($strand eq "-") { seq_SqstringReverseComplement(\$sub_str); }
    $sub_str =~ s/[\.\-\~]//g;  # final safety: strip any remaining gap chars
    if(vpatch_cds_is_complete($sub_str)) {
      my $protein = vpatch_translate_cds($sub_str);
      return ($sqname, $protein);
    }
  }
  return (undef, undef);
}

#################################################################
# Subroutine: vpatch_cds_is_complete()
# Mirrors profile_ValidateCdsIsComplete in v-build.pl: a valid CDS
# (a) has length divisible by 3 (b) ends in a stop codon
# (c) has no premature in-frame stops.
#################################################################
sub vpatch_cds_is_complete {
  my ($cds) = @_;
  my $len = length($cds);
  if($len < 6 || ($len % 3 != 0)) { return 0; }
  my %stop = ("TAA" => 1, "TAG" => 1, "TGA" => 1);
  my $tail = uc(substr($cds, $len - 3, 3));
  if(! exists $stop{$tail}) { return 0; }
  for(my $i = 0; $i < $len - 3; $i += 3) {
    my $c = uc(substr($cds, $i, 3));
    if(exists $stop{$c}) { return 0; }
  }
  return 1;
}

#################################################################
# Subroutine: vpatch_translate_cds()
# Translate a complete CDS (standard genetic code) to a protein
# string, stripping the terminal stop codon.
#################################################################
sub vpatch_translate_cds {
  my ($cds) = @_;
  my %codon = (
    TTT=>'F', TTC=>'F', TTA=>'L', TTG=>'L',
    CTT=>'L', CTC=>'L', CTA=>'L', CTG=>'L',
    ATT=>'I', ATC=>'I', ATA=>'I', ATG=>'M',
    GTT=>'V', GTC=>'V', GTA=>'V', GTG=>'V',
    TCT=>'S', TCC=>'S', TCA=>'S', TCG=>'S',
    CCT=>'P', CCC=>'P', CCA=>'P', CCG=>'P',
    ACT=>'T', ACC=>'T', ACA=>'T', ACG=>'T',
    GCT=>'A', GCC=>'A', GCA=>'A', GCG=>'A',
    TAT=>'Y', TAC=>'Y', TAA=>'*', TAG=>'*',
    CAT=>'H', CAC=>'H', CAA=>'Q', CAG=>'Q',
    AAT=>'N', AAC=>'N', AAA=>'K', AAG=>'K',
    GAT=>'D', GAC=>'D', GAA=>'E', GAG=>'E',
    TGT=>'C', TGC=>'C', TGA=>'*', TGG=>'W',
    CGT=>'R', CGC=>'R', CGA=>'R', CGG=>'R',
    AGT=>'S', AGC=>'S', AGA=>'R', AGG=>'R',
    GGT=>'G', GGC=>'G', GGA=>'G', GGG=>'G',
  );
  my $up = uc($cds);
  my $aa = "";
  for(my $i = 0; $i + 3 <= length($up) - 3; $i += 3) {  # stop one codon before the terminal stop
    my $c = substr($up, $i, 3);
    $aa .= exists $codon{$c} ? $codon{$c} : 'X';
  }
  return $aa;
}

#################################################################
# Subroutine: vpatch_finalize()
# Incept:     EPN* Thu May  7 2026
#
# Commit in-memory changes to disk:
#   1) backup current minfo + protein.fa to *.bak (unless already there)
#   2) for each applied add_alt op, extract a protein from the seed STK
#      and append (header: "<seedname>/<refcoords>") to protein.fa
#   3) write patched minfo via vdr_ModelInfoFileWrite()
#   4) rebuild BLAST db (makeblastdb) over the appended protein.fa
#   5) emit a .vpatch.log recording the applied ops
#################################################################
sub vpatch_finalize {
  my ($vb_dir, $root, $minfo_file, $protein_fa_file,
      $bak_minfo, $bak_prot, $log_file,
      $mdl_info_AHR, $ftr_info_HAHR, $ops_AR) = @_;

  # backups
  if(! -s $bak_minfo) {
    utl_RunCommand("cp $minfo_file $bak_minfo", opt_Get("-v", \%opt_HH), 0, \%FH_H);
  }
  if(! -s $bak_prot) {
    utl_RunCommand("cp $protein_fa_file $bak_prot", opt_Get("-v", \%opt_HH), 0, \%FH_H);
  }

  # protein extraction for any add_alt ops
  my $nalt = 0;
  for my $op (@{$ops_AR}) { $nalt++ if $op->{"op"} eq "add_alt"; }
  if($nalt > 0) {
    # Optional pre-translated protein library — headers must end with /<refcoords>
    my %altprot_by_coords_H = ();
    if(opt_IsUsed("--altprotein", \%opt_HH)) {
      my $apf = opt_Get("--altprotein", \%opt_HH);
      vpatch_load_altprotein_fasta($apf, \%altprot_by_coords_H);
      printf("# finalize: loaded %d pre-translated proteins from %s\n",
             scalar(keys %altprot_by_coords_H), $apf);
    }

    # STK source (only loaded if needed)
    my $stk_file = opt_IsUsed("--stk", \%opt_HH)
                     ? opt_Get("--stk", \%opt_HH)
                     : "$vb_dir/$root.stk";
    my $msa = undef;
    my $msa_load_failed = 0;

    open(my $prot_fh, ">>", $protein_fa_file)
      or die "ERROR, cannot append to $protein_fa_file: $!";
    for my $op (@{$ops_AR}) {
      next unless $op->{"op"} eq "add_alt";
      my $new_coords = $op->{"_new_cds_coords"};
      my ($seed, $protein) = (undef, undef);

      # 1) library lookup
      if(exists $altprot_by_coords_H{$new_coords}) {
        ($seed, $protein) = @{$altprot_by_coords_H{$new_coords}};
      }
      # 2) STK extraction fallback
      if(! defined $protein) {
        if(! defined $msa && ! $msa_load_failed) {
          if(! -s $stk_file) { $msa_load_failed = 1; }
          else {
            $msa = Bio::Easel::MSA->new({ fileLocation => $stk_file, isDna => 1 });
            if(! $msa->has_rf) {
              warn "WARNING, $stk_file lacks RF annotation; STK extraction disabled\n";
              $msa_load_failed = 1; $msa = undef;
            }
          }
        }
        if(defined $msa) {
          ($seed, $protein) = vpatch_extract_alt_protein($msa, $new_coords, \%FH_H);
        }
      }

      if(! defined $protein) {
        die "ERROR, finalize: no protein found for new alt coords $new_coords. "
          . "Tried --altprotein library (if any) and STK '$stk_file'. "
          . "Provide --altprotein <file> with a fasta entry whose header ends with /$new_coords, "
          . "or supply a --stk that contains a seed with the matching ORF.";
      }
      my $hdr = "$seed/$new_coords";
      print $prot_fh ">$hdr\n";
      for(my $i = 0; $i < length($protein); $i += 60) {
        print $prot_fh substr($protein, $i, 60) . "\n";
      }
      $op->{"_protein_seed"} = $seed;
      $op->{"_protein_len"}  = length($protein);
      printf("# finalize: appended protein '%s' (%d aa, from seed %s)\n",
             $hdr, length($protein), $seed);
    }
    close($prot_fh);
  }

  # sanity-check patched feature coords are within model length
  for(my $mi = 0; $mi < scalar(@{$mdl_info_AHR}); $mi++) {
    my $mname = $mdl_info_AHR->[$mi]{"name"};
    vdr_FeatureInfoValidateCoords($ftr_info_HAHR->{$mname}, $mdl_info_AHR->[$mi]{"length"}, \%FH_H);
  }

  # write patched minfo
  vdr_ModelInfoFileWrite($minfo_file, $mdl_info_AHR, $ftr_info_HAHR, \%FH_H);
  printf("# finalize: wrote patched minfo %s (backup at %s)\n", $minfo_file, $bak_minfo);

  # rebuild blastdb (only if any add_alt was applied; add_exc doesn't touch protein.fa)
  if($nalt > 0) {
    # remove stale blastdb index files so makeblastdb writes a clean set
    for my $ext (qw(phr pin psq pdb pjs pot ptf pto)) {
      my $f = "$protein_fa_file.$ext";
      if(-e $f) { unlink($f); }
    }
    my $cmd = $execs_H{"makeblastdb"} . " -in $protein_fa_file -dbtype prot > /dev/null";
    utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, \%FH_H);
    printf("# finalize: rebuilt BLAST db for %s\n", $protein_fa_file);
  }

  # write .vpatch.log
  open(my $lfh, ">", $log_file) or die "ERROR, cannot write $log_file: $!";
  printf $lfh "# v-patch.pl log -- %s\n", scalar localtime();
  printf $lfh "# vb-dir:           %s\n", $vb_dir;
  printf $lfh "# patched minfo:    %s\n", $minfo_file;
  printf $lfh "# minfo backup:     %s\n", $bak_minfo;
  printf $lfh "# protein.fa:       %s\n", $protein_fa_file;
  printf $lfh "# protein backup:   %s\n", $bak_prot;
  printf $lfh "# applied ops (%d):\n", scalar(@{$ops_AR});
  for my $op (@{$ops_AR}) {
    if($op->{"op"} eq "add_alt") {
      printf $lfh "add_alt\t%s\t%s\t%s\t%s\tsubn=%s\tseed=%s\tlen=%d\n",
        $op->{"model"}, $op->{"coords"}, $op->{"key"}, $op->{"value"},
        ($op->{"_set_subn"}  // ""), ($op->{"_protein_seed"} // ""),
        ($op->{"_protein_len"} // 0);
    } else {
      printf $lfh "add_exc\t%s\t%s\t%s\t%s\n",
        $op->{"model"}, $op->{"coords"}, $op->{"key"}, $op->{"value"};
    }
  }
  close($lfh);
  printf("# finalize: wrote %s\n", $log_file);
}
