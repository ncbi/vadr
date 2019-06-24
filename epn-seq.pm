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

require "epn-ofile.pm";

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
#  $codon:    the codon
#  $tt:       the translation table ('1' for standard)
#  $atg_only: return '1' only for 'ATG'
#
# Returns:    '1' if $codon is valid start codon for translation table $tt
#             else '0'
#
#################################################################
sub seq_CodonValidateStartCapDna {
  my $sub_name = "seq_CodonValidateStartCapDna";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($codon, $tt, $atg_only) = @_;

  if($atg_only) { 
    if($codon eq "ATG") {
      return 1; 
    }
    else { 
      return 0;
    }
  }
  else { 
    if($tt == 1) { # TTG CTG ATG
      if($codon =~ m/^[TCAYWMH]TG$/) { 
        return 1;
      }
    }
    elsif($tt == 2) { # ATT ATC ATA ATG GTG
      if(($codon =~ m/^AT[ACGTRYSWKMBDHVN]$/) || # ATN
         ($codon =~ m/^[AGR]TG$/)) { 
        return 1;
      }
    }
    elsif($tt == 3) { # ATA ATG GTG
      if(($codon =~ m/^AT[AGR]$/) || 
         ($codon =~ m/^[AGR]TG$/)) { 
        return 1;
      }
    }
    elsif($tt == 4) { # TTA TTG CTG ATT ATC ATA ATG GTG
      if(($codon =~ m/^AT[ACGTRYSWKMBDHVN]$/) || # ATN
         ($codon =~ m/^TT[AGR]$/) || 
         ($codon =~ m/^[CGS]TG$/)) { 
        return 1;
      }
    }
    elsif($tt == 5) { # TTG ATT ATC ATA ATG GTG
      if(($codon =~ m/^AT[ACGTRYSWKMBDHVN]$/) || # ATN
         ($codon =~ m/^[TGK]TG$/)) { 
        return 1;
      }
    }
    elsif(($tt == 6) || ($tt == 10) || ($tt == 14) || ($tt == 16) || ($tt == 22) || ($tt == 27) || ($tt == 28) || ($tt == 29) || ($tt == 30) || ($tt == 31)) { #ATG
      if(($codon eq "ATG")) {
        return 1;
      }
    }
    elsif(($tt == 9) || ($tt == 21)) { # ATG GTG
      if($codon =~ m/^[AGR]TG$/) { 
        return 1;
      }
    }
    elsif($tt == 11) { # TTG CTG ATT ATC ATA ATG GTG
      if(($codon =~ m/^AT[ACGTRYSWKMBDHVN]$/) || # ATN
         ($codon =~ m/^[TCGYKSB]TG$/)) { 
        return 1;
      }
    }
    elsif(($tt == 12) || ($tt == 26)) { # CTG ATG
      if($codon =~ m/^[CAM]TG$/) { 
        return 1;
      }
    }
    elsif($tt == 13) { # TTG ATA ATG GTG
      if(($codon =~ m/^[TAGWKRD]TG$/) || 
         ($codon =~ m/^AT[AGR]$/)) { 
        return 1;
      }
    }
    elsif($tt == 23) { # ATT ATG GTG
      if(($codon =~ m/^AT[TGK]$/) || 
         ($codon =~ m/^[AGR]TG$/)) { 
        return 1;
      }
    }
    elsif(($tt == 24) || ($tt == 33)) { # TTG CTG ATG GTG
      if($codon =~ m/^[ACGTRYSWKMBDHVN]TG$/) { # NTG
        return 1;
      }
    }
    elsif($tt == 25) { # TTG ATG GTG
      if($codon =~ m/^[TAGWKRD]TG$/) { 
        return 1;
      }
    }
  }

  return 0;
}

#################################################################
# Subroutine: seq_CodonValidateStopCapDna()
# Incept:     EPN, Mon Mar 14 13:47:57 2016
# 
# Purpose:    Given an already capitalized DNA codon, return '1' 
#             if it's a valid stop codon, else return 0.
#             Any codon with an ambiguous nt will return 0.
#
# Args:
#  $codon:  the codon
#  $tt:     the translation table ('1' for standard)
#
# Returns:    '1' if codon is valid stop, else '0'
#
# Ref: https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes
#################################################################
sub seq_CodonValidateStopCapDna {
  my $sub_name = "seq_CodonValidateStopCapDna";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($codon, $tt) = @_;
  
  if(($tt == 1) || ($tt == 11) || ($tt == 12) || ($tt == 26) || ($tt == 28)) { # TAA TAG TGA
    if(($codon =~ m/^TA[AGR]$/) || 
       ($codon =~ m/^T[AGR]A$/)) { 
      return 1;
    }
  }
  elsif($tt == 2) { # TAA TAG AGA AGG
    if(($codon =~ m/^TA[AGR]$/) || 
       ($codon =~ m/^AG[AGR]$/)) { 
      return 1;
    }
  }
  elsif(($tt == 3) || ($tt == 4) || ($tt == 5) || ($tt == 9) || ($tt == 10) || ($tt == 13) || ($tt == 21) || ($tt == 24) || ($tt == 25) || ($tt == 31)) { # TAA TAG
    if($codon =~ m/^TA[AGR]$/) {
      return 1;
    }
  }
  elsif(($tt == 6) || ($tt == 27) || ($tt == 29) || ($tt == 30)) { # TGA
    if($codon eq "TGA") { 
      return 1;
    }
  }
  elsif(($tt == 14) || ($tt == 33)) { # TAG
    if($codon eq "TAG") { 
      return 1;
    }
  }
  elsif($tt == 16) { # TAA TGA
    if($codon =~ m/^T[AGR]A$/) { 
      return 1;
    }
  }
  elsif($tt == 22) { # TCA TAA TGA
    if($codon =~ m/^T[ACGMSRV]A$/) { 
      return 1;
    }
  }
  elsif($tt == 23) { # TTA TAA TAG TGA
    if(($codon =~ m/^T[TAGWKRD]A$/) || 
       ($codon =~ m/^TA[AGR]$/)) { 
      return 1;
    }
  }
  
  return 0;
}

#################################################################
# Subroutine: seq_Overlap()
# Incept:     EPN, Mon Mar 14 13:47:57 2016 [dnaorg.pm]
#
# Purpose:    Calculate number of nucleotides of overlap between
#             two regions.
#
# Args:
#  $start1: start position of hit 1 (must be <= $end1)
#  $end1:   end   position of hit 1 (must be >= $end1)
#  $start2: start position of hit 2 (must be <= $end2)
#  $end2:   end   position of hit 2 (must be >= $end2)
#  $FH_HR:  REF to hash of file handles, including "log" and "cmd"
#
# Returns:  Two values:
#           $noverlap:    Number of nucleotides of overlap between hit1 and hit2, 
#                         0 if none
#           $overlap_reg: region of overlap, "" if none
#
# Dies:     if $end1 < $start1 or $end2 < $start2.
#
#################################################################
sub seq_Overlap {
  my $sub_name = "seq_Overlap";
  my $nargs_exp = 5;
  if(scalar(@_) != 5) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($start1, $end1, $start2, $end2, $FH_HR) = @_; 

  # printf("in $sub_name $start1..$end1 $start2..$end2\n");

  if($start1 > $end1) { ofile_FAIL("ERROR in $sub_name start1 > end1 ($start1 > $end1)", 1, $FH_HR); }
  if($start2 > $end2) { ofile_FAIL("ERROR in $sub_name start2 > end2 ($start2 > $end2)", 1, $FH_HR); }

  # Given: $start1 <= $end1 and $start2 <= $end2.
  
  # Swap if nec so that $start1 <= $start2.
  if($start1 > $start2) { 
    my $tmp;
    $tmp   = $start1; $start1 = $start2; $start2 = $tmp;
    $tmp   =   $end1;   $end1 =   $end2;   $end2 = $tmp;
  }
  
  # 3 possible cases:
  # Case 1. $start1 <=   $end1 <  $start2 <=   $end2  Overlap is 0
  # Case 2. $start1 <= $start2 <=   $end1 <    $end2  
  # Case 3. $start1 <= $start2 <=   $end2 <=   $end1
  if($end1 < $start2) { return (0, ""); }                                           # case 1
  if($end1 <   $end2) { return (($end1 - $start2 + 1), ($start2 . "-" . $end1)); }  # case 2
  if($end2 <=  $end1) { return (($end2 - $start2 + 1), ($start2 . "-" . $end2)); }  # case 3
  die "ERROR in $sub_name, unforeseen case in $sub_name, $start1..$end1 and $start2..$end2";

  return; # NOT REACHED
}

#################################################################
# Subroutine: seq_StripVersion()
# Incept:     EPN, Thu Feb 11 14:25:52 2016
#
# Purpose:    Given a ref to an accession.version string, remove the version.
#
# Arguments: 
#   $accver_R: ref to accession version string to remove version from
#
# Returns:    Nothing, $$accver_R has version removed
#
# Dies:       never
#################################################################
sub seq_StripVersion {
  my $sub_name  = "seq_StripVersion()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($accver_R) = (@_);

  $$accver_R =~ s/\.[0-9]*$//; # strip version

  return;
}

####################################################################
# the next line is critical, a perl module must return a true value
return 1;
####################################################################
