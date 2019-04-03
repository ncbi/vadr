#!/usr/bin/perl
#
# epn-seq.pm
# Eric Nawrocki
# EPN, Tue Mar 19 13:29:58 2019
# version: 0.00
#
use strict;
use warnings;
use Time::HiRes qw(gettimeofday);

#################################################################
# Subroutine:  seq_SqstringAddNewlines()
# Incept:      EPN, Thu Mar 14 06:12:11 2019
#
# Purpose:     Add newlines to $sqstring after every $linelen
#              characters and return result.
#
# Arguments: 
#   $sqstring: the sequence string
#   $linelen:  interval for newlines
# 
# Returns:     $sqstring with newlines inserted every $linelen 
#              characters and at end of string.
# 
# Dies:        Never.
#
################################################################# 
sub seq_SqstringAddNewlines { 
  my $nargs_expected = 2;
  my $sub_name = "seq_SqstringAddNewlines";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($sqstring, $linelen) = @_;

  my $retstr = "";
  my $sqpos = 0;
  my $sqlen = length($sqstring);

  while($sqpos < $sqlen) { 
    $retstr .= substr($sqstring, $sqpos, $linelen) . "\n";
    $sqpos += $linelen;
  }

  return $retstr;
}

#################################################################
# Subroutine: seq_SqstringCapitalize
# Incept:     EPN, Fri Mar 15 13:32:36 2019
# 
# Purpose:    Capitalize a string in place.
# 
# Arguments:
#   $sqstring_R: REF to sequence string to capitalize
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub seq_SqstringCapitalize {
  my $sub_name = "seq_SqstringCapitalize";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($sqstring_R) = @_;
  
  $$sqstring_R =~ tr/a-z/A-Z/;
  return;
}

#################################################################
# Subroutine: seq_SqstringDnaize
# Incept:     EPN, Fri Mar 15 13:33:39 2019
# 
# Purpose:    Convert a RNA/DNA sqstring to DNA in place.
# 
# Arguments:
#   $sqstring_R: REF to sequence string to capitalize
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub seq_SqstringDnaize {
  my $sub_name = "seq_SqstringDnaize";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($sqstring_R) = @_;
  
  $$sqstring_R =~ tr/Uu/Tt/;
  return;
}


#################################################################
# Subroutine: seq_SqstringReverseComplement
# Incept:     EPN, Fri Mar 15 15:35:10 2019
# 
# Purpose:    Reverse complement a RNA/DNA sqstring in place.
# 
# Arguments:
#   $sqstring_R: REF to reverse complement
#
# Returns:    void
# 
# Dies:       never
#
#################################################################
sub seq_SqstringReverseComplement {
  my $sub_name = "seq_SqstringReverseComplement";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($sqstring_R) = @_;

  # DNA-ize it
  seq_SqstringDnaize($sqstring_R);
  # reverse it 
  $$sqstring_R = reverse $$sqstring_R;
  # complement it
  $$sqstring_R =~ tr/ACGTRYMKHBVDacgtrymkhbvd/TGCAYRKMDVBHtgcayrkmdvbh/;
  # see esl_alphabet.c::set_complementarity()
  # note that S, W, N are omitted as they are their own complements

  return;
}


#################################################################
# Subroutine: seq_SqstringDiffSummary
# Incept:     EPN, Fri Mar 15 13:35:28 2019
# 
# Purpose:    Return a string summarizes the differences between
#             two sqstrings.
# 
# Arguments:
#   $sqstring1: sqstring 1
#   $sqstring2: sqstring 2
#
# Returns:    String with N newlines for N differences.
#             "" if sqstrings are identical.
# 
# Dies:       never
#
#################################################################
sub seq_SqstringDiffSummary {
  my $sub_name = "seq_SqstringDiffSummary";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { die "ERROR $sub_name entered with wrong number of input args" }
 
  my ($sqstring1, $sqstring2) = @_;

  if(! defined $sqstring1) { 
    return "sequence 1 is undefined\n"; 
  }
  if(! defined $sqstring2) { 
    return "sequence 2 is undefined\n"; 
  }
  seq_SqstringCapitalize(\$sqstring1);
  seq_SqstringCapitalize(\$sqstring2);
  seq_SqstringDnaize(\$sqstring1);
  seq_SqstringDnaize(\$sqstring2);
  if($sqstring1 eq $sqstring2) { 
    return "";
  }
  my $n1 = length($sqstring1); 
  my $n2 = length($sqstring2); 
  if($n1 != $n2) {
    return "sequence lengths mismatch ($n1 != $n2)\n"; 
  }
  my $ret_str = "";
  my @sqstring1_A = split("", $sqstring1); 
  my @sqstring2_A = split("", $sqstring2); 
  my $n = ($n1 > $n2) ? $n1 : $n2;
  for(my $i = 0; $i < $n; $i++) { 
    if($i >= $n1) { 
      $ret_str .= " char " . ($i+1) . " seq1: off-end seq2: " . $sqstring2_A[$i] . "\n";
    }
    elsif($i >= $n2) { 
      $ret_str .= " char " . ($i+1) . " seq1: " . $sqstring1_A[$i] . " seq2: off-end\n";
    }
    elsif($sqstring1_A[$i] ne $sqstring2_A[$i]) { 
      $ret_str .= " char " . ($i+1) . " seq1: " . $sqstring1_A[$i] . " seq2: " . $sqstring2_A[$i] . "\n";
    }
  }

  return $ret_str;
}

#################################################################
# Subroutine: seq_CodonValidateStartCapDna()
# Incept:     EPN, Sat Feb 23 10:01:55 2019
# 
# Purpose:    Given an already capitalized DNA codon, return '1' 
#             if it's a valid start codon, else return 0.
#
# Args:
#  $codon:  the codon
#
# Returns:    '1' if $codon is "ATG" else '0'
#
#################################################################
sub seq_CodonValidateStartCapDna {
  my $sub_name = "seq_CodonValidateStartCapDna";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($codon) = @_;
  
  if($codon eq "ATG") { 
    return 1;
  }

  return 0;
}

#################################################################
# Subroutine: seq_CodonValidateStopCapDna()
# Incept:     EPN, Mon Mar 14 13:47:57 2016
# 
# Purpose:    Given an already capitalized DNA codon, return '1' 
#             if it's a valid stop codon, else return 0.
#
# Args:
#  $codon:  the codon
#
# Returns:    '1' if codon is valid stop, else '0'
#
#################################################################
sub seq_CodonValidateStopCapDna {
  my $sub_name = "seq_CodonValidateStopCapDna";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($codon) = @_;
  
  if($codon eq "TAA" || 
     $codon eq "TGA" || 
     $codon eq "TAG" || 
     $codon eq "TAR") { 
    return 1;
  }

  return 0;
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
