#!usr/bin/env perl
my $usage = "perl revcomp-minfo.pl <minfo file> <model length>";
if(scalar(@ARGV) != 2) { die $usage; }

my ($minfo_file, $mlen) = (@ARGV);

open(MINFO, $minfo_file) || die "ERROR unable to open $minfo_file";
while($line = <MINFO>) { 
  if($line =~ m/^FEATURE/) {
    chomp $line;
    my $new_line = "";
    my @el_A = split(/\s+/, $line);
    foreach $el (@el_A) {
      if($el =~ /^coords\:\"(.+)\"/) { 
        my $new_coords_str = "coords:\"";
        $coords_str = $1;
        my @coords_token_A = split(",", $coords_str);
        foreach $coords_token (@coords_token_A) {
          if($coords_token =~ /^(\d+)\.\.(\d+)\:([+-])$/) {
            my ($start, $stop, $strand) = ($1, $2, $3);
            if($new_coords_str ne "coords:\"") { $new_coords_str .= ","; }
            $new_coords_str .= sprintf("%d..%d:%s", $mlen - $start + 1, $mlen - $stop + 1, $strand eq "+" ? "-" : "+");
          }
          else {
            die "ERROR unable to parse coords token $coords_token\n";
          }
        }
        if($new_line ne "") { $new_line .= " "; }
        $new_line .= "$new_coords_str\"";
      }
      else {
        if($new_line ne "") { $new_line .= " "; }
        $new_line .= $el;
      }
    }
    print $new_line . "\n";
  }
  else {
    print $line;
  }
}

      
