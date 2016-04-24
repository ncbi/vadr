#!/usr/bin/env perl
# 
# copy-evalue-parameters.pl: copy E-value parameters from one CM file into another
# EPN, Sun Apr 24 09:15:37 2016
# 
# Each CM file must only have 1 CM in it.
# Destination CM file must not have E value parameters in it.
# 
my $usage;
$usage  = "copy-evalue-parameters.pl <CM file with parameters to copy FROM> <CM file with parameters to copy TO>\n";

if(scalar(@ARGV) != 2) { die $usage; }
($src_cmfile, $dst_cmfile) = @ARGV;

open(SRC, $src_cmfile) || die "ERROR unable to open $src_cmfile for reading";

my @src_evalue_lines_A = ();
while($line = <SRC>) { 
  if($line =~ m/^ECM/) { 
    push(@src_evalue_lines_A, $line);
  }
}
close(SRC);

my $nsrc_evalue_lines = scalar(@src_evalue_lines_A);
if($nsrc_evalue_lines != 4) { 
  die "ERROR, expected to read 4 lines starting with ECM in $src_cmfile, but read $nsrc_evalue_lines";
}

my @dst_A = (); # array of all lines for new CM file output
my $nefp7_lines = 0;
open(DST, $dst_cmfile) || die "ERROR unable to open $dst_cmfile for reading";
while($line = <DST>) { 
  if($line =~ m/^ECM/) { 
    die "ERROR, read E-value parameter line in $dst_cmfile, didn't expect any: $line"; 
  }
  push(@dst_A, $line);
  if($line =~ m/EFP7GF/) { # we need to insert the E-value parameters after this line
    push(@dst_A, @src_evalue_lines_A);
    $nefp7_lines++;
  }
}
close(DST);

if($nefp7_lines == 0) { 
  die "ERROR, didn't read any EFP7GF lines in $dst_cmfile, expected exactly one";
}
if($nefp7_lines != 1) { 
  die "ERROR, read $nefp7_lines lines in $dst_cmfile, it probably had $nefp7_lines CMs in it. This script only works on CM files with 1 model in them";
}

# create new CM file
$cmd = "cp " . $dst_cmfile . " " . $dst_cmfile . ".old";
system("$cmd");
open(OUT, ">", $dst_cmfile) || die "ERROR unable to open $dst_cmfile for writing";
print("# Copied $dst_cmfile to $dst_cmfile.old\n");
foreach my $line (@dst_A) { 
  print OUT $line;
}
close(OUT);
print("# Created new $dst_cmfile with E-value parameters from $src_cmfile\n");

