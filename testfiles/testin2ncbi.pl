#!/usr/bin/perl -w
# Convert a vadr .testin input file to v-test.pl to 
# a command file with vadr_* deployed artifact commands.

use strict;
use warnings;
use Getopt::Long;

my $usage;
$usage  = "testin2ncbi.pl\n\n";
$usage .= "Usage:\n";
$usage .= "\ttestin2ncbi.pl <v-test.pl .testin input file>\n";
#$usage .= "\ttestin2ncbi.pl [OPTIONS]\n";
#$usage .= "\t<v-test.pl .testin input file>\n";

#my $do_tab    = 0;   # set to '1' if -t used
#my $add_empty = "-"; # set to <s> if --empty <s>used
#&GetOptions( "t"       => \$do_tab,
#             "empty=s" => \$add_empty);

if(scalar(@ARGV) != 1) { die $usage; }
my ($in_testin) = @ARGV;

my @cmd_A  = ();
my @desc_A = ();

#my $exec_path = "$NCBI/bin/_latest/MSS/vadr";

open(IN, $in_testin) || die "ERROR, unable to open $in_testin for reading";
while(my $line = <IN>) { 
  chomp $line;
  if($line =~ /^command\:\s+(.+)$/) { 
    push(@cmd_A, $1);
  }
  if($line =~ /^desc\:\s+(.+)$/) { 
    push(@desc_A, $1); 
  }
}
close(IN);

my $ncmd  = scalar(@cmd_A);
my $ndesc = scalar(@desc_A);
if($ncmd != $ndesc) { 
  die "ERROR different number of commands and descs"; 
}
print "### tests from .testin file $in_testin\n";

# for each command line, create the corresponding vadr_* command
for(my $i = 0; $i < $ncmd; $i++) { 
  my $cmd  = $cmd_A[$i];
  my $desc = $desc_A[$i];
  my $is_noro   = ($desc =~ /noro/)   ? 1 : 0;
  my $is_dengue = ($desc =~ /dengue/) ? 1 : 0;

  my $short_cmd; 
  if($cmd =~ /^perl \$VADRSCRIPTSDIR\/v-annotate.pl\s+(.+)$/) { 
    $short_cmd = $1;
  }
  else { 
    die "ERROR unable to parse $cmd\n";
  }

  # default to vadr_noro unless model directory is not supplied and dengue is in the desc
  my $exec_name = "vadr_noro";
  if(($short_cmd !~ /\-\-mdir/) && 
     (($short_cmd !~ /\-m /) || ($short_cmd !~ /\-i /))) { # model directory not specified
    if($is_dengue) { 
      $exec_name = "vadr_dengue";
    }      
    elsif(! $is_noro) { 
      die "ERROR in $in_testin, no --mdir option but desc does not have noro or dengue for desc: $desc and cmd:$cmd\n";
    }
  }

  print "#" . $desc . "\n";
  #print $exec_path . "/" . $exec_name . " -- " . $short_cmd . "\n";
  print $exec_name . " -- " . $short_cmd . "\n";
}
