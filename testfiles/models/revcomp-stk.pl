#!usr/bin/env perl
my $usage = "perl revcomp-stk.pl <stk file> <model length> <old name: e.g. ENTOY100A> <new name: e.g. ENTOY100A-REV";
if(scalar(@ARGV) != 4) { die $usage; }

my ($stk_file, $mlen, $old_name, $new_name) = (@ARGV);

open(STK, $stk_file) || die "ERROR unable to open $minfo_file";
while($line = <STK>) { 
  chomp $line;
  if($line =~ /^#=GS(\s+)(\S+)(\s+)DE(\s+)(.+)$/) { 
    # DE line
    my ($space1, $name, $space2, $space3, $remainder) = ($1, $2, $3, $4, $5);
    $name =~ s/$old_name/$new_name/;
    printf("#=GS%s%s%sDE%s$remainder\n", $space1, $name, $space2, $space3);
  }
  elsif($line =~ /^\#=GC(\s+)(\S+)(\s+)(.+)$/) { 
    # GC annotation, reverse string and update name
    my ($space1, $name, $space2, $remainder) = ($1, $2, $3, $4);
    if($name =~ m/^SS\_cons/) { 
      $remainder =~ tr/\<\>/\>\</;
      printf("#=GC" . $space1 . $name . $space2 . reverse($remainder) . "\n");
    }
    elsif($name =~ m/^COL/) { 
      printf("#=GC%s%s%s%s\n", $space1, $name, $space2, $remainder);
    }
    elsif($name =~ /(.+\.)(\d+)\-(\d+)/) { 
      my ($beg, $start, $end) = ($1, $2, $3);
      $name = $beg . ($mlen - $start + 1) . "-" . ($mlen - $end + 1);
      $remainder =~ tr/\<\>/\>\</;
      printf("#=GC" . $space1 . $name . $space2 . reverse($remainder) . "\n");
      #print("GC line\n");
      #print("name: $name\n");
      #print("remainder: $remainder\n");
      #print("rev remainder: " . reverse($remainder) . "\n");
      #exit(1);
    }
    else { 
      die "ERROR unable to parse GC line $line";
    }
  }
  elsif($line =~ /^\#/) { 
    print $line . "\n";
  }
  elsif($line =~ /(\S+)(\s+)(\S+)/) { 
    # sequence line
    my ($name, $space, $seq) = ($1, $2, $3);
    $name =~ s/$old_name/$new_name/;
    $seq = reverse($seq);
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    printf("%s%s%s\n", $name, $space, $seq);
  }
  else { 
    print $line . "\n";
  }
}
close(STK);

      
