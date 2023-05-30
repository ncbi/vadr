#!/usr/bin/env perl
# EPN, Fri May 19 14:01:42 2023
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;
use LWP::Simple; 

require "vadr.pm"; 
require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

#######################################################################################
# What this script does: 
#
# - Reads in a model info file
# - Reads in a stockholm alignment file
# - Creates a fasta file of all CDS sequences in the sequences from the alignment
# - Translates the CDS sequences into protein sequences
# - Creates a blast db from the protein sequences
# 
#######################################################################################

# first, determine the paths to all modules, scripts and executables that we'll need

# make sure the required environment variables are set
my $env_vadr_scripts_dir  = utl_DirEnvVarValid("VADRSCRIPTSDIR");
my $env_vadr_blast_dir    = utl_DirEnvVarValid("VADRBLASTDIR");
my $env_vadr_easel_dir    = utl_DirEnvVarValid("VADREASELDIR");

my %execs_H = (); # hash with paths to all required executables
$execs_H{"esl-translate"} = $env_vadr_easel_dir . "/esl-translate";
$execs_H{"makeblastdb"}   = $env_vadr_blast_dir . "/makeblastdb";
utl_ExecHValidate(\%execs_H, undef);

#########################################################
# Command line and option processing using sqp_opts.pm
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
#     option            type       default               group   requires incompat    preamble-output                                                help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                                         "display this help",                                   \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",                                 "force; if dir <output directory> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                                  "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--ttbl",       "integer", 1,                        1,    undef, undef,      "use NCBI translation table <n> to translate CDS",             "use NCBI translation table <n> to translate CDS", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",                          "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'ttbl=s'       => \$GetOptions_H{"--ttbl"},
                'keep'         => \$GetOptions_H{"--keep"});

my $total_seconds = -1 * ofile_SecondsSinceEpoch(); # by multiplying by -1, we can just add another ofile_SecondsSinceEpoch call at end to get total time
my $execname_opt  = $GetOptions_H{"--execname"};
my $executable    = "build-add-to-blast-db.pl";
my $usage         = "Usage: $executable [-options]\n\t<path to .minfo file>\n\t<path to blast db dir>\n\t<model name>\n\t<nt-accn-to-add>\n\t\<nt-coords-to-add>\n\t<model-CDS-feature-coords>\n\t<name for output directory>\n";
my $synopsis      = "$executable :: add a single protein to a VADR blastx protein database";
my $date          = scalar localtime();
my $version       = "1.5.1";
my $releasedate   = "Feb 2023";
my $pkgname       = "VADR";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, undef);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 7) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do $executable -h\n\n";
  exit(1);
}
my ($in_minfo_file, $in_blastdb_dir, $in_mdl_name, $in_nt_accn, $in_nt_cds_coords, $in_mdl_cds_coords, $dir) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# make sure that $in_blastdb_dir exists
my $mdl_blastdb_file = $in_blastdb_dir . "/" . $in_mdl_name . ".vadr.protein.fa";
if(! -d $in_blastdb_dir) { 
  die "ERROR, blast db directory $in_blastdb_dir does not exist";
}
if(! -e $mdl_blastdb_file) {
  die "ERROR, blast db file $mdl_blastdb_file does not exist";
}

#############################
# create the output directory
#############################
my $cmd;              # a command to run with utl_RunCommand()
my @early_cmd_A = (); # array of commands we run before our log file is opened

if($dir =~ m/\/$/) { 
  $dir =~ s/\/$//; # remove final '/' if it exists
} 
if(-d $dir) { 
  $cmd = "rm -rf $dir";
  if(opt_Get("-f", \%opt_HH)) { utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR directory named $dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $dir) { 
  $cmd = "rm $dir";
  if(opt_Get("-f", \%opt_HH)) { utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR a file named $dir already exists. Remove it, or use -f to overwrite it."; }
}

# create the dir
$cmd = "mkdir $dir";
utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root         = $dir . "/" . $dir_tail . ".vadr";

#######################
# output program banner
#######################
# output preamble
my @arg_desc_A = ("input model info file", "input blast db path", "input model name", "nucleotide accession to add", "nt coords of CDS to add", "CDS feature coords this CDS should map to", "output directory");
my @arg_A      = ($in_minfo_file, $in_blastdb_dir, $in_mdl_name, $in_nt_accn, $in_nt_cds_coords, $in_mdl_cds_coords, $dir);
my %extra_H    = ();
$extra_H{"\$VADRSCRIPTSDIR"}  = $env_vadr_scripts_dir;
$extra_H{"\$VADRBLASTDIR"}    = $env_vadr_blast_dir;
$extra_H{"\$VADREASELDIR"}    = $env_vadr_easel_dir;
ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
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
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "log", $out_root . ".log", 1, 1, "Output printed to screen");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "cmd", $out_root . ".cmd", 1, 1, "List of executed commands");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".filelist", 1, 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
my $FH_HR  = $ofile_info_HH{"FH"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# now we have the log file open, output the banner there too
ofile_OutputBanner($log_FH, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

###########################
# Parse the model info file
###########################
my @mdl_info_AH  = (); # array of hashes with model info
my %ftr_info_HAH = (); # hash of array of hashes with feature info 

my $progress_w = 50;
my $start_secs = ofile_OutputProgressPrior("Parsing input model info file", $progress_w, undef, *STDOUT);

my @reqd_mdl_keys_A = ("name", "length");
my @reqd_ftr_keys_A = ("type", "coords");
utl_FileValidateExistsAndNonEmpty($in_minfo_file, "model info file", undef, 1, $FH_HR);
vdr_ModelInfoFileParse($in_minfo_file, \@reqd_mdl_keys_A, \@reqd_ftr_keys_A, \@mdl_info_AH, \%ftr_info_HAH, $FH_HR);

# make sure the specified model exists
if(! defined $ftr_info_HAH{$in_mdl_name}) { 
  ofile_FAIL("ERROR did not read info for any models in $in_minfo_file", 1, $FH_HR); 
}
  
# make sure we can find the appropriate feature
my $nftr = scalar(@{$ftr_info_HAH{$in_mdl_name}});
my $ftr_idx = undef;
my $found_ftr_idx = undef;
for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
  if(vdr_FeatureTypeIsCds(\@{$ftr_info_HAH{$in_mdl_name}}, $ftr_idx)) { 
    #printf("comparing %s and %s\n", $ftr_info_HAH{$in_mdl_name}[$ftr_idx]{"coords"}, $in_nt_cds_coords);
    if($ftr_info_HAH{$in_mdl_name}[$ftr_idx]{"coords"} eq $in_mdl_cds_coords) { 
      $found_ftr_idx = $ftr_idx;
      $ftr_idx = $nftr; # breaks loop
    }
  }
}
if(! defined $found_ftr_idx) { 
  ofile_FAIL("ERROR did not find any CDS with coords $in_mdl_cds_coords for model $in_mdl_name in input file $in_minfo_file", 1, $FH_HR); 
}
ofile_OutputProgressComplete($start_secs, undef, undef, *STDOUT);

###############################################
# Fetch the source sequence and CDS subsequence
###############################################
$start_secs = ofile_OutputProgressPrior("Fetching the CDS source sequence ", $progress_w, undef, *STDOUT);
my $source_fa_file = $out_root . ".source.fa";
vdr_EutilsFetchToFile($source_fa_file, $in_nt_accn, "nuccore", "fasta", 5, $ofile_info_HH{"FH"});  # number of attempts to fetch to make before dying
ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "sourcefasta", $source_fa_file, 1, 1, "fasta file with source sequence ($in_nt_accn)");

my $sqfile  = Bio::Easel::SqFile->new({ fileLocation => $source_fa_file }); # the sequence file object
my @start_A  = ();
my @stop_A   = ();
my @strand_A = ();
vdr_FeatureStartStopStrandArrays($in_nt_cds_coords, \@start_A, \@stop_A, \@strand_A, $FH_HR);
my $nsgm = scalar(@start_A);
my $cds_sqstring = "";
my $seq_name = $sqfile->fetch_seq_name_given_ssi_number(0); # we know the ssi index must be 0 because there's only 1 seq
for(my $sgm_idx = 0; $sgm_idx < $nsgm; $sgm_idx++) { 
  $cds_sqstring .= $sqfile->fetch_subseq_to_sqstring($seq_name, $start_A[$sgm_idx], $stop_A[$sgm_idx], ($strand_A[$sgm_idx] eq "-"));
}

# output the CDS to a fasta file
my $cds_fa_file = $out_root . ".cds.fa";
open(FA, ">", $cds_fa_file) || ofile_FileOpenFailure($cds_fa_file, "main", $!, "writing", $FH_HR);
my $cds_seqname = sprintf("%s/%s", $seq_name, $in_nt_cds_coords);
printf FA (">$cds_seqname\n%s", seq_SqstringAddNewlines($cds_sqstring, 60));
close(FA);
ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "cdsfasta", $cds_fa_file, 1, 1, "fasta file with CDS from $in_nt_accn");
ofile_OutputProgressComplete($start_secs, undef, undef, *STDOUT);

###################
# Translate the CDS
###################
$start_secs = ofile_OutputProgressPrior("Translating CDS ", $progress_w, undef, *STDOUT);
my $tmp_translate_fa_file = $out_root . ".prot.fa.tmp";
my $translate_fa_file = $out_root . ".prot.fa";
my $db_seqname = sprintf("%s:%s/%s", $seq_name, $in_nt_cds_coords, $in_mdl_cds_coords);
my $cds_nt_len = vdr_CoordsLength($in_nt_cds_coords, $FH_HR);
if(($cds_nt_len % 3) != 0) { 
  ofile_FAIL("ERROR CDS length calculated as $cds_nt_len from coords $in_nt_cds_coords, which is not a multiple of 3", 1, $FH_HR);
}
my $exp_prot_len = ($cds_nt_len / 3) - 1;
my $c_opt = "";
if(opt_IsUsed("--ttbl", \%opt_HH)) { 
  $c_opt = "-c " . opt_Get("--ttbl", \%opt_HH);
}
my $translate_cmd = $execs_H{"esl-translate"} . " $c_opt -M -l $exp_prot_len --watson $cds_fa_file > $tmp_translate_fa_file";
utl_RunCommand($translate_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);

# go through the output fasta file and rewrite the seq name
my $seqname_seen = 0;
open(IN,       $tmp_translate_fa_file) || ofile_FileOpenFailure($tmp_translate_fa_file, "main", $!, "reading", $FH_HR);
open(OUT, ">", $translate_fa_file)     || ofile_FileOpenFailure($translate_fa_file,     "main", $!, "writing", $FH_HR);
while(my $line = <IN>) { 
  if($line =~ m/^\>/) { 
    #>orf58 source=NC_039477.1/5..5104:+ coords=1..5097 length=1699 frame=1  
    chomp $line;
    if($line =~ /^\>orf\d+\s+source\=\S+\s+coords\=(\d+)\.\.(\d+)\s+length\=\d+\s+frame\=\S+/) { 
      # rename as 'source=NC_039477.1/5..5104:+,coords=1..5097'
      my ($start, $stop) = ($1, $2);
      if($start != 1)               { ofile_FAIL("ERROR CDS not translated as expected, start is $start but should be 1",           1, $FH_HR); }
      if($stop  != ($cds_nt_len-3)) { ofile_FAIL("ERROR CDS not translated as expected, stop  is $stop  but should be $cds_nt_len", 1, $FH_HR); }
      if($seqname_seen == 1) { 
        ofile_FAIL("ERROR more than one sequence exists in $tmp_translate_fa_file", 1, $FH_HR);
      }
      print OUT (">$db_seqname\n");
      $seqname_seen = 1;
    } 
    else { 
      ofile_FAIL("ERROR problem parsing esl-translate output file $tmp_translate_fa_file, line:\n$line\n", 1, $FH_HR);
    }
  }
  else { 
    print OUT $line; 
  }
}
close(OUT);
close(IN);
ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "protfasta", $translate_fa_file, 1, 1, "fasta file with translated protein from $in_nt_accn");
ofile_OutputProgressComplete($start_secs, undef, undef, *STDOUT);

##########################################
# Concontenate the protein to the BLAST DB 
##########################################
$start_secs = ofile_OutputProgressPrior("Adding to BLAST DB", $progress_w, undef, *STDOUT);
my $concat_cmd = "cat $translate_fa_file >> $mdl_blastdb_file";
utl_RunCommand($concat_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);

my $makeblastdb_cmd = $execs_H{"makeblastdb"} . " -in $mdl_blastdb_file -dbtype prot > /dev/null";
utl_RunCommand($makeblastdb_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
ofile_OutputProgressComplete($start_secs, undef, undef, *STDOUT);

##########
# Conclude
##########

$total_seconds += ofile_SecondsSinceEpoch();
ofile_OutputConclusionAndCloseFilesOk($total_seconds, "", \%ofile_info_HH);
exit 0;
