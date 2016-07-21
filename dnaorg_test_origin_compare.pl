use strict;

my $usage = "perl dnaorg_test_origin_compare.pl <dnaorg_test_origin_*.pl output 1> <dnaorg_test_origin_*.pl output 2> <origin length> <output root>";

if(scalar(@ARGV) != 4) { 
  die $usage;
}

my ($in1, $in2, $clen, $out_root) = @ARGV;

my %info1_HH = (); # file 1: 1st dim key: sequence name, 2nd dim key: "coords", "len", "cseq", "nmismatch", "strand", "passfail", "line";
my %info2_HH = (); # file 2: 1st dim key: sequence name, 2nd dim key: "coords", "len", "cseq", "nmismatch", "strand", "passfail", "line";

read_input_file($in1, \%info1_HH);
read_input_file($in2, \%info2_HH);

my $n1 = scalar(keys %info1_HH);
my $n2 = scalar(keys %info2_HH);
if($n1 != $n2) { die "ERROR different number of sequences in file 1 and file 2 ($n1 != $n2)"; }

my $nseq = 0;
my $neq_pred = 0;
my $nneq1_pred = 0;
my $nneq2_pred = 0;

my $neq_nopred = 0;
my $nneq1_nopred = 0;
my $nneq2_nopred = 0;

my $neq_wronglen = 0;
my $nneq1_wronglen = 0;
my $nneq2_wronglen = 0;

my $neq_correctlen = 0;
my $nneq1_correctlen = 0;
my $nneq2_correctlen = 0;

my %neq_mismatch_H = ();
my %nneq1_mismatch_H = ();
my %nneq2_mismatch_H = ();

my @neq_pred_A = ();
my %neq_mismatch_HA = ();
my @neq_correctlen_A = ();
my @neq_wronglen_A = ();

my @neq_nopred_A = ();

my %nneq1_mismatch_HA = ();
my %nneq2_mismatch_HA = ();
my @nneq1_correctlen_A = ();
my @nneq2_correctlen_A = ();
my @nneq1_wronglen_A = ();
my @nneq2_wronglen_A = ();

my @nneq1_pred_A = ();
my @nneq2_pred_A = ();
my @nneq1_nopred_A = ();
my @nneq2_nopred_A = ();

for(my $z = 0; $z <= $clen; $z++) { 
  @{$neq_mismatch_HA{$z}} = ();
  @{$nneq1_mismatch_HA{$z}} = ();
  @{$nneq2_mismatch_HA{$z}} = ();
}

foreach my $seqname (sort keys %info1_HH) { 
  $nseq++;
  if(! exists $info2_HH{$seqname}) { 
    die "ERROR seqname: $seqname in $in1 but not $in2"; 
  }
  my $have_pred1 = ($info1_HH{$seqname}{"coords"} eq "?") ? 0 : 1;
  my $have_pred2 = ($info2_HH{$seqname}{"coords"} eq "?") ? 0 : 1;
  my $is_equal   = ($info1_HH{$seqname}{"line"} eq $info2_HH{$seqname}{"line"}) ? 1 : 0;
  my $len1       = $info1_HH{$seqname}{"len"};
  my $len2       = $info2_HH{$seqname}{"len"};
  my $nmismatch1 = $info1_HH{$seqname}{"nmismatch"};
  my $nmismatch2 = $info2_HH{$seqname}{"nmismatch"};

  if($is_equal) { 
    if($have_pred1 && $have_pred2) { 
      # equal predictions
      $neq_pred++; 
      push(@neq_pred_A, $seqname);
      if($len1 == $clen) { 
        $neq_mismatch_H{$nmismatch1}++; 
        push(@{$neq_mismatch_HA{$nmismatch1}}, $seqname); 
        $neq_correctlen++; 
        push(@neq_correctlen_A, $seqname); 
      }
      else { 
        $neq_wronglen++; 
        push(@neq_wronglen_A, $seqname);
      }
    }
    else { 
      # no predictions (both)
      $neq_nopred++;
      push(@neq_nopred_A, $seqname);
      if(! ((! $have_pred1) && (! $have_pred2))) { 
        die "ERROR logic error 1"; 
      }
    }
  }
  else { # not equal (lines are different)
    # printf("%s\n%s\n\n", $info1_HH{$seqname}{"line"}, $info2_HH{$seqname}{"line"});

    if($have_pred1 && $have_pred2) { 
      # both have predictions that are not equal
      $nneq1_pred++; 
      push(@nneq1_pred_A, $seqname);
      $nneq2_pred++; 
      push(@nneq2_pred_A, $seqname);

      if($len1 == $clen) { 
        $nneq1_mismatch_H{$nmismatch1}++; 
        push(@{$nneq1_mismatch_HA{$nmismatch1}}, $seqname);
        $nneq1_correctlen++; 
        push(@nneq1_correctlen_A, $seqname); 
      }
      else { 
        $nneq1_wronglen++; 
        push(@nneq1_wronglen_A, $seqname);
      }

      if($len2 == $clen) { 
        $nneq2_mismatch_H{$nmismatch2}++; 
        push(@{$nneq2_mismatch_HA{$nmismatch2}}, $seqname);
        $nneq2_correctlen++; 
        push(@nneq2_correctlen_A, $seqname);
      }
      else { 
        $nneq2_wronglen++; 
        push(@nneq2_wronglen_A, $seqname);
      }
    }
    elsif((! $have_pred1) && (! $have_pred2)) { 
      # neither have predictions, but they're different, shouldn't happen
      die "ERROR logic error 2";
    }
    elsif(! $have_pred1) { # have_pred2 is '1'
      # no prediction for 1, but prediction for 2
      $nneq2_pred++; 
      push(@nneq2_pred_A, $seqname);
      $nneq1_nopred++; 
      push(@nneq1_nopred_A, $seqname);
      if($len2 == $clen) { 
        $nneq2_mismatch_H{$nmismatch2}++; 
        push(@{$nneq2_mismatch_HA{$nmismatch2}}, $seqname);
        $nneq2_correctlen++; 
        push(@nneq2_correctlen_A, $seqname);
      }
      else { 
        $nneq2_wronglen++; 
        push(@nneq2_wronglen_A, $seqname);
      }
    }
    elsif(! $have_pred2) { # have_pred1 is '1'
      # no prediction for 2, but prediction for 1
      $nneq1_pred++; 
      push(@nneq1_pred_A, $seqname);
      $nneq2_nopred++; 
      push(@nneq2_nopred_A, $seqname);
      if($len1 == $clen) { 
        $nneq1_mismatch_H{$nmismatch1}++; 
        push(@{$nneq1_mismatch_HA{$nmismatch1}}, $seqname);
        $nneq1_correctlen++; 
        push(@nneq1_correctlen_A, $seqname);
      }
      else { 
        $nneq1_wronglen++; 
        push(@nneq1_wronglen_A, $seqname);
      }
    }
  }
}

# output summary
# print summary
printf("#\n# Summary:\n#\n");
printf("# Number of sequences:                       %4d\n", $nseq);
printf("# Number of no predictions:                  %4d  %4d  %4d\n", $neq_pred,       $nneq1_pred,       $nneq2_pred);
printf("# Number of predictions of unexpected len:   %4d  %4d  %4d\n", $neq_wronglen,   $nneq1_wronglen,   $nneq2_wronglen);
printf("# Number of predictions of expected len:     %4d  %4d  %4d\n", $neq_correctlen, $nneq1_correctlen, $nneq2_correctlen);
for(my $z = 0; $z <= $clen; $z++) { 
  printf("# Number of predictions with %2d mismatches:  %4d  %4d  %4d\n", $z, 
         (exists $neq_mismatch_H{$z})   ? $neq_mismatch_H{$z} : 0, 
         (exists $nneq1_mismatch_H{$z}) ? $nneq1_mismatch_H{$z} : 0, 
         (exists $nneq2_mismatch_H{$z}) ? $nneq2_mismatch_H{$z} : 0);
}

# output sequences to files
printf("#\n");
output_info_to_file($out_root . ".neq_pred.compare",            \@neq_pred_A,         \%info1_HH, \%info2_HH);
output_info_to_file($out_root . ".nneq1_pred.compare",          \@nneq1_pred_A,       \%info1_HH, \%info2_HH);
output_info_to_file($out_root . ".nneq2_pred.compare",          \@nneq2_pred_A,       \%info1_HH, \%info2_HH);
output_info_to_file($out_root . ".neq_unexpected.compare",      \@neq_wronglen_A,     \%info1_HH, \%info2_HH);
output_info_to_file($out_root . ".nneq1_unexpectedlen.compare", \@nneq1_wronglen_A,   \%info1_HH, \%info2_HH);
output_info_to_file($out_root . ".nneq2_unexpectedlen.compare", \@nneq2_wronglen_A,   \%info1_HH, \%info2_HH);
output_info_to_file($out_root . ".neq_expectedlen.compare",     \@neq_correctlen_A,   \%info1_HH, \%info2_HH);
output_info_to_file($out_root . ".nneq1_expectedlen.compare",   \@nneq1_correctlen_A, \%info1_HH, \%info2_HH);
output_info_to_file($out_root . ".nneq2_expectedlen.compare",   \@nneq2_correctlen_A, \%info1_HH, \%info2_HH);
for(my $z = 0; $z <= $clen; $z++) { 
  output_info_to_file($out_root . ".neq_mismatch.$z.compare",   \@{$neq_mismatch_HA{$z}},   \%info1_HH, \%info2_HH);
  output_info_to_file($out_root . ".nneq1_mismatch.$z.compare", \@{$nneq1_mismatch_HA{$z}}, \%info1_HH, \%info2_HH);
  output_info_to_file($out_root . ".nneq2_mismatch.$z.compare", \@{$nneq2_mismatch_HA{$z}}, \%info1_HH, \%info2_HH);
}


exit 0;                                                                     

########################
sub read_input_file {
  my ($infile, $info_HHR) = @_;
  open(IN, $infile) || die "ERROR unable to open $infile for reading";
  while(my $line = <IN>) { 
    chomp $line;
    if($line !~ m/^\#/) { 
      ## example line:
      #AF003952:dnaorg-duplicated:AF003952:1:2690:+:AF003952:1:2690:+:                   2684..2692   9   TAATATTAC   0  + PASS
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) != 7) { die "ERROR unable to parse line from $infile:\n$line\n"; }
      my ($seqname, $coords, $len, $cseq, $nmismatch, $strand, $passfail) = @el_A;
      if(exists $info_HHR->{$seqname}) { 
        die "ERROR duplicate line for $seqname in $infile";
      }
      %{$info_HHR->{$seqname}} = ();
      $info_HHR->{$seqname}{"coords"}    = $coords;
      $info_HHR->{$seqname}{"len"}       = $len;
      $info_HHR->{$seqname}{"cseq"}      = $cseq;
      $info_HHR->{$seqname}{"nmismatch"} = $nmismatch;
      $info_HHR->{$seqname}{"strand"}    = $strand;
      $info_HHR->{$seqname}{"passfail"}  = $passfail;
      $info_HHR->{$seqname}{"line"}      = $line;
    }
  }
  return;
}

########################
sub output_info_to_file {
  my ($outfile, $AR, $info1_HHR, $info2_HHR) = @_;

  if(scalar(@{$AR}) == 0) { return; }

  open(OUT, ">", $outfile) || die "ERROR unable to open $outfile for writing";

  foreach my $seqname (@{$AR}) { 
    printf OUT ("%s\n%s\n\n", $info1_HHR->{$seqname}{"line"}, $info2_HHR->{$seqname}{"line"});
  }
  close(OUT);
  printf("# Information on %4d sequences printed to output file $outfile\n", scalar(@{$AR}));

  return;
}
