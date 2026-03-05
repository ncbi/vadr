my $usage = "perl fasta-rpn-add-model.pl <fa file> <model name>\n";
if(scalar(@ARGV) != 2) { 
  die $usage;
}
my ($fa_file, $model) = (@ARGV);

open(IN, $fa_file) || die "ERROR unable to open $fa_file for reading"; 

while($line = <IN>) { 
  chomp $line;
  if($line =~ /^\>(\S+).*/) { 
    print(">$1:MODEL:$model\n");
  }
  else { 
    print $line . "\n"; 
  }
}
