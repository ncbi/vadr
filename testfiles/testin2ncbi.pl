#!/usr/bin/perl -w
# Convert a vadr .testin input file to v-test.pl to 
# a command file with vadr_* deployed artifact commands.

use strict;
use warnings;
use Getopt::Long;

my $usage;
$usage  = "testin2ncbi.pl\n\n";
$usage .= "Usage:\n";
$usage .= "\ttestin2ncbi.pl [OPTIONS]\n";
$usage .= "\t<v-test.pl .testin input file>\n";

my $do_noopts = 0;
my $opt_execname = undef;
&GetOptions( "noopts"     => \$do_noopts,
             "execname=s" => \$opt_execname);

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

  my $short_cmd = undef; # command without executable
  my ($opts, $arg1, $arg2, $outfile); 
  # commands need to be in a specific format, or else --noopts won't work
  # we assume the final two args are the command line args (not options, for v-annotate.pl this will be fasta file and output dir)
  # and that the output is piped to a file
  if($cmd =~ /^perl \$VADRSCRIPTSDIR\/v-annotate.pl\s+(.+)\s+(\S+)\s+(\S+)\s+\>\s+(\S+)/) { 
    ($opts, $arg1, $arg2, $outfile) = ($1, $2, $3, $4);
    $short_cmd = sprintf("%s%s %s > %s", (($do_noopts) ? "" : $opts . " "), $arg1, $arg2, $outfile);
  }
  else { 
    die "ERROR unable to parse $cmd\n";
  }

  # determine name of executable (e.g. vadr_noro)
  # if --execname was used, use that
  my $execname = undef;
  if(defined $opt_execname) { 
    $execname = $opt_execname;
  }
  else { 
    # default to vadr_noro unless model directory is not supplied and dengue is in the desc
    $execname = "vadr_noro";
    if(($short_cmd !~ /\-\-mdir/) && 
       (($short_cmd !~ /\-m /) || ($short_cmd !~ /\-i /))) { # model directory not specified
      if($is_dengue) { 
        $execname = "vadr_dengue";
      }      
      elsif(! $is_noro) { 
        die "ERROR in $in_testin, no --mdir option but desc does not have noro or dengue for desc: $desc and cmd:$cmd\n";
      }
    }
  }

  printf("#%s%s\n", $desc, ($do_noopts ? " (with all options removed)" : ""));
  #print $exec_path . "/" . $execname . " -- " . $short_cmd . "\n";

  if($do_noopts) { 
    print $execname . " " .  $short_cmd . "\n";
  }
  else { 
    print $execname . " -- " . $short_cmd . "\n";
  }
}
