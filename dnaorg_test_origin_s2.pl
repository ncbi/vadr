#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Mon Feb  1 15:07:43 2016 [dnaorg_build.pl split off from dnaorg_annotate_genomes.pl]
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

require "dnaorg.pm"; 
require "epn-options.pm";

#######################################################################################
# What this script does: 
#
# Preliminaries: 
#   - process options
#   - create the output directory
#   - output program banner and open output files
#   - parse the optional input files, if necessary
#   - make sure the required executables are executable
#
# Step 1. Gather and process information on reference genome using Edirect
#
# Step 2. Fetch and process the reference genome sequence
#
# Step 3. Build and calibrate models
#######################################################################################

# first, determine the paths to all modules, scripts and executables that we'll need
# make sure the DNAORGDIR environment variable is set
my $dnaorgdir = $ENV{'DNAORGDIR'};
if(! exists($ENV{'DNAORGDIR'})) { 
    printf STDERR ("\nERROR, the environment variable DNAORGDIR is not set, please set it to the directory where you installed the dnaorg scripts and their dependencies.\n"); 
    exit(1); 
}
if(! (-d $dnaorgdir)) { 
    printf STDERR ("\nERROR, the dnaorg directory specified by your environment variable DNAORGDIR does not exist.\n"); 
    exit(1); 
}    
 
# determine other required paths to executables relative to $dnaorgdir
my $inf_exec_dir      = $dnaorgdir . "/infernal-dev/src/";
my $hmmer_exec_dir    = $dnaorgdir . "/hmmer-3.1b2/src/";
my $esl_exec_dir      = $dnaorgdir . "/infernal-dev/easel/miniapps/";
my $esl_fetch_cds     = $dnaorgdir . "/esl-fetch-cds/esl-fetch-cds.pl";

my $nnop  = 0; # number of sequences for which an origin is not predicted
my $npred = 0; # number of sequences for which an origin is predicted
my $npred_len = 0; # number of sequences for which an origin is predicted of expected len
my %nmismatch_H = ();

#########################################################
# Command line and option processing using epn-options.pm
#
# opt_HH: 2D hash:
#         1D key: option name (e.g. "-h")
#         2D key: string denoting type of information 
#                 (one of "type", "default", "group", "requires", "incompatible", "preamble", "help")
#         value:  string explaining 2D key:
#                 "type":          "boolean", "string", "int" or "real"
#                 "default":       default value for option
#                 "group":         integer denoting group number this option belongs to
#                 "requires":      string of 0 or more other options this option requires to work, each separated by a ','
#                 "incompatiable": string of 0 or more other options this option is incompatible with, each separated by a ','
#                 "preamble":      string describing option for preamble section (beginning of output from script)
#                 "help":          string describing option for help section (printed if -h used)
#                 "setby":         '1' if option set by user, else 'undef'
#                 "value":         value for option, can be undef if default is undef
#
# opt_order_A: array of options in the order they should be processed
# 
# opt_group_desc_H: key: group number (integer), value: description of group for help output
my %opt_HH = ();      
my @opt_order_A = (); 
my %opt_group_desc_H = ();

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
$opt_group_desc_H{"1"} = "basic options";
#     option            type       default               group   requires incompat    preamble-output                          help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                   "display this help",                                  \%opt_HH, \@opt_order_A);
opt_Add("--hmmonly",    "boolean", 0,                        1,    undef, undef,      "search with HMMs not CMs",              "search with HMMs not CMs",                           \%opt_HH, \@opt_order_A);


# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_test_origin_s2.pl [-options] <CM file with origin model> <start position of origin in CM> <fasta file> <output directory> <consensus sequence>\n";
my $synopsis = "dnaorg_test_origin_s2.pl :: search for origin sequences [TEST SCRIPT]";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
                'hmmonly'      => \$GetOptions_H{"--hmmonly"});

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.40";
my $releasedate   = "Dec 2018";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date, $dnaorgdir);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 5) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do dnaorg_build.pl -h\n\n";
  exit(1);
}
my ($model_file, $ori_start_rfpos, $fasta_file, $dir_out, $cons_seq) = (@ARGV);

if(defined $dir_out) { 
  $dir_out =~ s/\/$//; # remove final '/' if there is one
}
my $dir_out_tail   = $dir_out;
$dir_out_tail   =~ s/^.+\///; # remove all but last dir
my $out_root   = $dir_out .   "/" . $dir_out_tail   . ".dnaorg_test_origin_s2";

my $cmd;
if(! -d $dir_out) {
  $cmd = "mkdir $dir_out";
  runCommand($cmd, 0, undef);
}

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("model file with origin model", "origin start position", "fasta file", "output file root", "consensus origin sequence");
my @arg_A      = ($model_file, $ori_start_rfpos, $fasta_file, $dir_out, $cons_seq);
outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date, $dnaorgdir);
opt_OutputPreamble(*STDOUT, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# open the log and command files:
# set output file names and file handles, and open those file handles
my %ofile_info_HH = ();  # hash of information on output files we created,
                         # 1D keys: 
                         #  "fullpath":  full path to the file
                         #  "nodirpath": file name, full path minus all directories
                         #  "desc":      short description of the file
                         #  "FH":        file handle to output to for this file, maybe undef
                         # 2D keys:
                         #  "log": log file of what's output to stdout
                         #  "cmd": command file with list of all commands executed

# open the log and command files 
openAndAddFileToOutputInfo(\%ofile_info_HH, "log", $out_root . ".log", 1, "Output printed to screen");
openAndAddFileToOutputInfo(\%ofile_info_HH, "cmd", $out_root . ".cmd", 1, "List of executed commands");
openAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".list", 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# now we have the log file open, output the banner there too
outputBanner($log_FH, $version, $releasedate, $synopsis, $date, $dnaorgdir);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

my $cons_len = length($cons_seq);
my @cons_seq_A = split("", $cons_seq);
###################################################
# make sure the required executables are executable
###################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmscan"}       = $inf_exec_dir . "cmscan";
$execs_H{"cmalign"}      = $inf_exec_dir . "cmalign";
$execs_H{"esl-alimask"}  = $esl_exec_dir . "esl-alimask";
$execs_H{"esl-reformat"} = $esl_exec_dir . "esl-reformat";
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

#####################################################
# Determine length of all sequences in the fasta file
#####################################################
if(-e $fasta_file . ".ssi") { 
  unlink $fasta_file . ".ssi";
}
if(-e $fasta_file.".ssi") { unlink $fasta_file.".ssi"; }
my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $fasta_file }); # the sequence file object
my $nseq = $sqfile->nseq_ssi;
my %seqlen_H = ();
my @seq_order_A = ();
for(my $i = 0; $i < $nseq; $i++) { 
  my ($seqname, $seqlen) = $sqfile->fetch_seq_name_and_length_given_ssi_number($i);
  $seqlen_H{$seqname} = $seqlen / 2;
  push(@seq_order_A, $seqname);
}

###############################
# Run cmscan on all sequences
###############################
my $tblout_file = $out_root . ".tbl";
my $stdout_file = $out_root . ".cmscan";
my $opts = " --cpu 0 --tblout $tblout_file --verbose ";
if(! opt_Get("--hmmonly", \%opt_HH)) { 
  $opts .= " --nohmmonly --F1 0.02 --F2 0.001 --F2b 0.001 --F3 0.00001 --F3b 0.00001 --F4 0.0002 --F4b 0.0002 --F5 0.0002 --F6 0.0001 ";
}

$cmd = $execs_H{"cmscan"} . " " . $opts . " $model_file $fasta_file > $stdout_file";

runCommand($cmd, 0, $ofile_info_HH{"FH"});

addClosedFileToOutputInfo(\%ofile_info_HH, "tblout", "$tblout_file", 1, "cmscan tabular output");
addClosedFileToOutputInfo(\%ofile_info_HH, "stdout", "$stdout_file", 1, "cmscan standard output");

################################################
# Parse cmscan tabular output and output results
################################################

my %hit_HH  = (); # 2D hash of top hits, 1st dim key is sequence name, 2nd is attribute, e.g. "start"    
parse_cmscan_tblout_s2($tblout_file, \%seqlen_H, \%hit_HH, $ofile_info_HH{"FH"});

#################################################################
# Fetch all hits into a fasta file that we can align with cmalign
#################################################################
my @fetch_AA = ();
foreach my $seqname (@seq_order_A) { 
  if(exists $hit_HH{$seqname}) { 
    my $newname = $seqname . "/" . $hit_HH{$seqname}{"start"} . "-" . $hit_HH{$seqname}{"stop"};
    push(@fetch_AA, [ $newname, $hit_HH{$seqname}{"start"}, $hit_HH{$seqname}{"stop"}, $seqname ]);
    #printf("pushed $newname\n");
  }
}
my $out_fasta_file = $out_root . ".cmscan.fa";
$sqfile->fetch_subseqs(\@fetch_AA, undef, $out_fasta_file);
addClosedFileToOutputInfo(\%ofile_info_HH, "outfasta", "$out_fasta_file", 1, "cmscan hits in fasta format");
    
#############################
# Align all hits with cmalign
#############################

my $out_stk_file     = $out_root . ".cmalign.stk";
my $out_cmalign_file = $out_root . ".cmalign";
$cmd = $execs_H{"cmalign"} . " -o $out_stk_file $model_file $out_fasta_file > $out_cmalign_file";

runCommand($cmd, 0, $ofile_info_HH{"FH"});
addClosedFileToOutputInfo(\%ofile_info_HH, "outstk",     "$out_stk_file",     1, "alignment of cmscan hits");
addClosedFileToOutputInfo(\%ofile_info_HH, "outcmalign", "$out_cmalign_file", 1, "cmalign output");

#########################################
# Pull out the predicted origin sequences
#########################################
my $out_origin_fa_file = $out_root . ".origin.fa";
my $ori_stop_rfpos = ($ori_start_rfpos + $cons_len - 1);
my $ori_coords = $ori_start_rfpos . ".." . $ori_stop_rfpos;
$cmd = $execs_H{"esl-alimask"} . " --t-rf -t $out_stk_file $ori_coords | " . $execs_H{"esl-reformat"} . " -d --informat stockholm fasta - > $out_origin_fa_file";

runCommand($cmd, 0, $ofile_info_HH{"FH"});
addClosedFileToOutputInfo(\%ofile_info_HH, "originfasta",     "$out_origin_fa_file",     1, "fasta file of all predicted origin sequences");

################
# Output results
################
if(-e $out_origin_fa_file.".ssi") { unlink $out_origin_fa_file.".ssi"; }
my $ori_sqfile = Bio::Easel::SqFile->new({ fileLocation => $out_origin_fa_file }); 
my $ori_msa    = Bio::Easel::MSA->new   ({ fileLocation => $out_stk_file });
my $ori_start_apos = $ori_msa->rfpos_to_aligned_pos($ori_start_rfpos);
my $ori_stop_apos  = $ori_msa->rfpos_to_aligned_pos($ori_stop_rfpos);

foreach my $seqname (@seq_order_A) { 
  my $seqlen = $seqlen_H{$seqname};
  my $has_hit = (defined $hit_HH{$seqname}) ? 1 : 0;

  if($has_hit) { 
    my $start = $hit_HH{$seqname}{"start"};
    my $stop  = $hit_HH{$seqname}{"stop"};
    my $msa_seqname = $seqname . "/" . $start . "-" . $stop;
    my $msa_sqidx = $ori_msa->get_sqidx($msa_seqname);
    if($msa_sqidx == -1) { 
      DNAORG_FAIL("ERROR, unable to find $msa_seqname in $out_stk_file", 1, $ofile_info_HH{"FH"});
    }

    # determine the unaligned positions the origin spans in the alignment file $out_stk_file
    my ($cur_ori_start_uapos, $cur_ori_start_apos) = $ori_msa->aligned_to_unaligned_pos($msa_sqidx, $ori_start_apos, 1); # '1': if ori_start_apos is a gap, return position of 1st non-gap nucleotide after it 
    my ($cur_ori_stop_uapos,  $cur_ori_stop_apos)  = $ori_msa->aligned_to_unaligned_pos($msa_sqidx, $ori_stop_apos,  0); # '0': if ori_start_apos is a gap, return position of 1st non-gap nucleotide before it 

    # do we have at least 1 nucleotide predicted in the origin positions?
    if(($cur_ori_start_apos < $ori_stop_apos) && ($cur_ori_stop_apos > $ori_start_apos))  { 
      # yes, we do:
      my $cur_ori_len = $cur_ori_stop_uapos - $cur_ori_start_uapos + 1;
      # determine strand 
      if($start <= $stop) { 
        # positive strand
        my $cur_ori_coords = ($cur_ori_start_uapos + $start - 1) . ".." . ($cur_ori_stop_uapos + $start - 1);
        # fetch origin sequence, and make sure it matches what we fetched in $out_origin_fa_file
        #printf("seqname $seqname start: %d (%d+%d-1) stop: %d (%d+%d-1)\n", $cur_ori_start_uapos + $start - 1, $cur_ori_start_uapos, $start, $cur_ori_stop_uapos + $start - 1, $cur_ori_stop_uapos, $start);
        my $ori_fasta_seq1 = $sqfile->fetch_subseq_to_fasta_string($seqname, $cur_ori_start_uapos + $start - 1, $cur_ori_stop_uapos + $start -1, -1, 0);
        my $ori_fasta_seq2 = $ori_sqfile->fetch_seq_to_fasta_string($msa_seqname, -1);
        (my $cur_ori_seq1 = $ori_fasta_seq1) =~ s/^\>.+\n//;
        (my $cur_ori_seq2 = $ori_fasta_seq2) =~ s/^\>.+\n//;
        $cur_ori_seq1 =~ tr/a-z/A-Z/;
        $cur_ori_seq2 =~ tr/a-z/A-Z/;
        if($cur_ori_seq1 ne $cur_ori_seq2) { 
          DNAORG_FAIL("ERROR, seqname: $seqname fetched origin check failed:\n$cur_ori_seq1\nne\n$cur_ori_seq2\n", 1, $ofile_info_HH{"FH"});
        }
        chomp $cur_ori_seq1;
        my $nmismatch = compare_to_consensus($cur_ori_seq1, \@cons_seq_A);

        outputString($ofile_info_HH{"FH"}{"log"}, 1, sprintf("%-80s  %10s  %2d  %10s  %2d  + %s\n", $seqname, $cur_ori_coords, $cur_ori_len, $cur_ori_seq1, $nmismatch, ($nmismatch == 0) ? "PASS" : "FAIL"));
        $npred++;
        if($cur_ori_len == $cons_len) { 
          $npred_len++;
          $nmismatch_H{$nmismatch}++;
        }
      } # end of 'if($start < $stop)'
      else { 
        # negative strand
        my $cur_ori_coords = ($start - $cur_ori_start_uapos + 1) . ".." . ($start - $cur_ori_stop_uapos + 1);
        # fetch origin sequence, and make sure it matches what we fetched in $out_origin_fa_file
        #printf("seqname $seqname start: %d (%d+%d-1) stop: %d (%d+%d-1)\n", $cur_ori_start_uapos + $start - 1, $cur_ori_start_uapos, $start, $cur_ori_stop_uapos + $start - 1, $cur_ori_stop_uapos, $start);
        my $ori_fasta_seq1 = $sqfile->fetch_subseq_to_fasta_string($seqname, ($start - $cur_ori_start_uapos + 1), ($start - $cur_ori_stop_uapos + 1), -1, 0);
        my $ori_fasta_seq2 = $ori_sqfile->fetch_seq_to_fasta_string($msa_seqname, -1);
        (my $cur_ori_seq1 = $ori_fasta_seq1) =~ s/^\>.+\n//;
        (my $cur_ori_seq2 = $ori_fasta_seq2) =~ s/^\>.+\n//;
        $cur_ori_seq1 =~ tr/a-z/A-Z/;
        $cur_ori_seq2 =~ tr/a-z/A-Z/;
        if($cur_ori_seq1 ne $cur_ori_seq2) { 
          DNAORG_FAIL("ERROR, seqname: $seqname fetched origin check failed:\n$cur_ori_seq1\nne\n$cur_ori_seq2\n", 1, $ofile_info_HH{"FH"});
        }
        chomp $cur_ori_seq1;
        my $nmismatch = compare_to_consensus($cur_ori_seq1, \@cons_seq_A);

        outputString($ofile_info_HH{"FH"}{"log"}, 1, sprintf("%-80s  %10s  %2d  %10s  %2d  - %s\n", $seqname, $cur_ori_coords, $cur_ori_len, $cur_ori_seq1, $nmismatch, ($nmismatch == 0) ? "PASS" : "FAIL"));
        $npred++;
        if($cur_ori_len == $cons_len) { 
          $npred_len++;
          $nmismatch_H{$nmismatch}++;
        }
      }
    }
    else { 
      outputString($ofile_info_HH{"FH"}{"log"}, 1, sprintf("%-80s  %10s  %2s  %10s  %2d  ? FAIL\n", $seqname, "?", "?", "?", $cons_len));
      $nnop++;
    }
  } 
  else { 
    outputString($ofile_info_HH{"FH"}{"log"}, 1, sprintf("%-80s  %10s  %2s  %10s  %2d  ? FAIL\n", $seqname, "?", "?", "?", $cons_len));
    $nnop++;
  }
}


##########
# Conclude
##########

# print summary
outputString($ofile_info_HH{"FH"}{"log"}, 1, "#\n# Summary:\n#\n");
outputString($ofile_info_HH{"FH"}{"log"}, 1, sprintf("# Number of sequences:                       %4d\n", $nseq));
outputString($ofile_info_HH{"FH"}{"log"}, 1, sprintf("# Number of no predictions:                  %4d (%.3f)\n", $nnop, $nnop / $nseq));
outputString($ofile_info_HH{"FH"}{"log"}, 1, sprintf("# Number of predictions of unexpected len:   %4d (%.3f)\n", ($npred-$npred_len), ($npred-$npred_len) / $nseq));
outputString($ofile_info_HH{"FH"}{"log"}, 1, sprintf("# Number of predictions of expected len:     %4d (%.3f)\n", $npred_len, $npred_len / $nseq));
for(my $z = 0; $z <= $cons_len; $z++) { 
  my $cur_nmismatch = (exists $nmismatch_H{$z}) ? $nmismatch_H{$z} : 0;
  outputString($ofile_info_HH{"FH"}{"log"}, 1, sprintf("# Number of predictions with %2d mismatches:  %4d (%.3f)\n", $z, $cur_nmismatch, $cur_nmismatch / $npred_len));
}

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir_out, \%ofile_info_HH);
exit 0;


#################################################################
# Subroutine : parse_cmscan_tblout_s2()
# Incept:      EPN, Tue Jul 12 08:54:07 2016
#
# Arguments: 
#  $tblout_file: tblout file to parse
#  $seqlen_HR:    REF to hash, key is sequence name, value is length
#  $hit_HHR:      REF to 2D hash of top hits, 1st dim key is sequence name, 2nd is attribute, e.g. "start"    
#  $FH_HR:        REF to hash of file handles
#
# Returns:    void
#
#################################################################
sub parse_cmscan_tblout_s2 { 
  my $sub_name = "parse_cmscan_tblout_s2()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($tblout_file, $seqlen_HR, $hit_HHR, $FH_HR) = @_;
  
  open(IN, $tblout_file) || fileOpenFailure($tblout_file, $sub_name, $!, "reading", $FH_HR);

  my $did_field_check = 0; # set to '1' below after we check the fields of the file
  my $line_ctr = 0;  # counts lines in tblout_file
  while(my $line = <IN>) { 
    $line_ctr++;
    if(($line =~ m/^\#/) && (! $did_field_check)) { 
      # sanity check, make sure the fields are what we expect
      if($line !~ m/#target name\s+accession\s+query name\s+accession\s+mdl\s+mdl\s+from\s+mdl to\s+seq from\s+seq to\s+strand\s+trunc\s+pass\s+gc\s+bias\s+score\s+E-value inc description of target/) { 
        DNAORG_FAIL("ERROR in $sub_name, unexpected field names in $tblout_file\n$line\n", 1, $FH_HR);
      }
      $did_field_check = 1;
    }
    elsif($line !~ m/^\#/) { 
      chomp $line;
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
      # example line:
      # NC_001346.dnaorg_build.origin.5p -         KJ699341             -         hmm        1       59     2484     2542      +     -    6 0.59   0.1   78.5     2e-24 !   -
      my @elA = split(/\s+/, $line);
      my ($mdlname, $seqname, $mod, $mdlfrom, $mdlto, $seqfrom, $seqto, $strand, $score, $evalue) = 
          ($elA[0], $elA[2], $elA[4], $elA[5], $elA[6], $elA[7], $elA[8], $elA[9], $elA[14], $elA[15]);

      my $seqlen = $seqlen_HR->{$seqname};

      # only consider hits where either the start or end are less than the total length
      # of the genome. Since we sometimes duplicate all genomes, this gives a simple 
      # rule for deciding which of duplicate hits we'll store 
      if(($seqfrom <= $seqlen) || ($seqto <= $seqlen)) { 
        if(! exists $hit_HHR->{$seqname}) { 
          %{$hit_HHR->{$seqname}} = ();
          $hit_HHR->{$seqname}{"start"}  = $seqfrom;
          $hit_HHR->{$seqname}{"stop"}   = $seqto;
          $hit_HHR->{$seqname}{"score"}  = $score;
          $hit_HHR->{$seqname}{"evalue"} = $evalue;
        } 
      }
    }
  }
}

#################################################################
# Subroutine : compare_to_consensus()
# Incept:      EPN, Tue Jul 12 14:33:02 2016
#
# Purpose:    Given a predicted origin sequence, compare the
#             it to the consensus sequence, and report number
#             of mismatches.
#
# Arguments: 
#  $pred_seq:    predicted consensus sequence
#  $cons_seq_AR: REF to an array that is the consensus sequence, each element is an array element
#
# Returns:    void
#
#################################################################
sub compare_to_consensus { 
  my $sub_name = "compare_to_consensus()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($pred_seq, $cons_seq_AR) = @_;

  my @pred_A = split("", $pred_seq);
  my $cons_len = scalar(@{$cons_seq_AR});
  my $pred_len = scalar(@pred_A);

  my $nmatch = 0;
  my $min_len = ($cons_len < $pred_len) ? $cons_len : $pred_len;

  for(my $i = 0; $i < $min_len; $i++) { 
    if(uc($cons_seq_AR->[$i]) eq uc($pred_A[$i])) { 
      $nmatch++; 
    }
  }
  
  return $cons_len - $nmatch;
}  

sub tmp_aligned_to_unaligned_pos
{
  my ($msa, $sqidx, $apos, $do_after) = @_;

  if(! defined $do_after) { $do_after = 0; }

  $msa->_check_msa();
  $msa->_check_sqidx($sqidx);
  $msa->_check_ax_apos($apos);

  my $sqstring = $msa->get_sqstring_aligned($sqidx);
  my $ret_apos; # return apos
  my $uapos;    # return value, the unaligned position corresponding to $ret_apos

  #printf("sqstring: $sqstring sqidx: $sqidx apos: $apos\n");
  
  # is alignment position $apos a gap in $sqidx?
  my $apos_char = substr($sqstring, $apos-1, 1);
  #printf("apos_char: $apos_char\n");
  my $is_gap = ($apos_char =~ m/[a-zA-Z]/) ? 0 : 1;
  #print("is_gap: $is_gap\n");

  if(! $is_gap) { 
    # easy case
    # remove all non-alphabetic characters, to get unaligned length
    my $sqstring_to_apos_no_gaps = substr($sqstring, 0, $apos);
    $sqstring_to_apos_no_gaps =~ s/[^a-zA-Z]//g;
    $uapos = length($sqstring_to_apos_no_gaps);
    return ($uapos, $apos);
  }
  else { 
    # $apos is a gap for $sqidx:
    # determine last  position before $apos that is not a gap, if any (if ! $do_after)
    #        or first position after  $apos that is not a gap, if any (if $do_after)
    if(! $do_after) { 
      # $do_after is '0': determine last  position before $apos that is not a gap, if any (if ! $do_after)
      # first check if there are any characters that are not gaps:
      # remove all characters after apos, we don't care about them
      my $sqstring_to_apos = substr($sqstring, 0, $apos);
      if ($sqstring_to_apos =~ /[a-zA-Z]/) {
        # we have at least 1 non-gap
        (my $sqstring_to_apos_no_trailing_gaps = $sqstring_to_apos) =~ s/[^a-zA-Z]*$//;
        #printf("before sqstring_to_apos_no_trailing_gaps1: $sqstring_to_apos_no_trailing_gaps\n");
        $ret_apos = length($sqstring_to_apos_no_trailing_gaps);
        #printf("ret_apos: $ret_apos\n");
        (my $sqstring_to_apos_no_gaps = $sqstring_to_apos_no_trailing_gaps) =~ s/[^a-zA-Z]//g; # remove all gaps from sqstring_no_gaps
        $uapos = length($sqstring_to_apos_no_gaps); # length of substr_no_gaps gives us uapos
        #printf("uapos: $uapos\n");
      }
      else { # no alphabetic characters in the string
        $ret_apos = -1;
        $uapos    = -1;
      }
    }
    else { 
      # $do_after is '1': determine first position after  $apos that is not a gap, if any 
      my $sqstring_apos_to_alen = substr($sqstring, $apos-1); # we want to examine from $apos to $alen (remember apos is 1..alen, not 0..alen-1)
      if ($sqstring_apos_to_alen  =~ /[a-zA-Z]/) {
        # we have at least 1 non-gap
        (my $sqstring_apos_to_alen_no_leading_gaps = $sqstring_apos_to_alen) =~ s/^[^a-zA-Z]*//;
        #printf("after sqstring_apos_to_alen_no_leading_gaps1: $sqstring_apos_to_alen_no_leading_gaps\n");
        $ret_apos  = $msa->alen - length($sqstring_apos_to_alen_no_leading_gaps) + 1; # the +1 is to account for the fact that we didn't remove the first nt
        #printf("ret_apos: $ret_apos\n");
        (my $sqstring_apos_to_alen_no_gaps = $sqstring_apos_to_alen_no_leading_gaps) =~ s/[^a-zA-Z]//g; # remove all gaps from sqstring_apos_to_alen_no_gaps
        (my $sqstring_no_gaps = $sqstring) =~ s/[^a-zA-Z]//g;
        $uapos = length($sqstring_no_gaps) - length($sqstring_apos_to_alen_no_gaps) + 1; # again, +1 b/c we didn't remove the first nt;
        #printf("uapos: $uapos\n");
      }
      else { # no alphabetic characters in the string
        $ret_apos = -1;
        $uapos    = -1;
      }
    }    
    return ($uapos, $ret_apos);
  }
}
