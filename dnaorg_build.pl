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

# hard-coded-paths:
my $inf_exec_dir   = "/usr/local/infernal/1.1.1/bin/";
my $esl_fetch_cds  = "/panfs/pan1/dnaorg/programs/esl-fetch-cds.pl";

# The definition of $usage explains the script and usage:
my $usage = "\ndnaorg_build.pl <reference accession>\n";
$usage .= "\n"; 
$usage .= " FIXME: This script annotates genomes from the same species based\n";
$usage .= " on reference annotation. The reference accession is the\n";
$usage .= " first accession listed in <list file with all accessions>\n";
$usage .= "\n";
$usage .= " BASIC/COMMON OPTIONS:\n";
$usage .= "  -c           : genome is circular\n";
$usage .= "  -d <s>       : define output directory as <s>, not <reference accession>\n";
$usage .= "  -f           : force; if dir <reference accession> exists, overwrite it\n";
$usage .= "  -v           : be verbose; output commands to stdout as they're run\n";
$usage .= "  -matpept <f> : read mat_peptide info in addition to CDS info, file <f> explains CDS:mat_peptide relationships\n";
$usage .= "  -nomatpept   : ignore mat_peptide information for <reference accession> (do not model it)\n";
$usage .= "\n OPTIONS SPECIFIC TO INFERNAL:\n";
$usage .= "  -cslow  : use default cmcalibrate parameters, not parameters optimized for speed\n";
$usage .= "  -clocal : run cmcalibrate locally, do not submit calibration jobs for each CM to the compute farm\n";
$usage .= "\n";

my ($seconds, $microseconds) = gettimeofday();
my $start_secs      = ($seconds + ($microseconds / 1000000.));
my $executable      = $0;


###################################################
# process options
###################################################
# basic/common options:
my $dir            = undef; # set to a value <s> with -d <s>
my $do_circular    = 0;     # set to '1' with -c, genome is circular and features can span stop..start boundary
my $do_force       = 0;     # set to '1' with -f, overwrite output dir if it exists
my $be_verbose     = 0;     # set to '1' with -v, output commands as they're run
my $do_matpept     = 0;     # set to '1' if -matpept    enabled, genome has a single polyprotein, use mat_peptide info, not CDS
my $do_nomatpept   = 0;     # set to '1' if -nomatpept  enabled, we will ignore mat_peptide information if it exists for the reference
my $matpept_infile = undef; # defined if -matpept   enabled, the input file that describes relationship between CDS and mat_peptides
my $do_cslow       = 0; # set to '1' if -cslow   enabled, use default, slow, cmcalibrate parameters instead of speed optimized ones
my $do_clocal      = 0; # set to '1' if -clocal  enabled, do not submit cmcalibrate job to farm

&GetOptions("d=s"         => \$dir,
            "c"           => \$do_circular,
            "f"           => \$do_force, 
            "v"           => \$be_verbose,
            "matpept=s"   => \$matpept_infile,
            "cslow"       => \$do_cslow, 
            "clocal"      => \$do_clocal) ||
    die "Unknown option";

if(scalar(@ARGV) != 1) { die $usage; }
my ($ref_accn) = (@ARGV);

# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
if(defined $dir) { 
  $opts_used_short .= "-d $dir ";
  $opts_used_long  .= "# option:  output directory specified as $dir [-d]\n"; 
}
if(defined $do_circular) { 
  $opts_used_short .= "-c ";
  $opts_used_long  .= "# option:  genome is circular, so features can span stop..start [-c]\n"; 
}
if($do_force) { 
  $opts_used_short .= "-f ";
  $opts_used_long  .= sprintf("# option:  forcing overwrite of %s directory [-f]\n", (defined $dir ? $dir : $ref_accn)); 
}
if(defined $matpept_infile) { 
  $do_matpept = 1;
  $opts_used_short .= "-matpept $matpept_infile";
  $opts_used_long  .= "# option:  using mat_peptide info, CDS:mat_peptide relationships explained in $matpept_infile [-matpept]\n";
}
if($do_cslow) { 
  $opts_used_short .= "-cslow ";
  $opts_used_long  .= "# option:  run cmcalibrate in default (slow) mode [-cslow]\n";
}
if($do_clocal) { 
  $opts_used_short .= "-clocal ";
  $opts_used_long  .= "# option:  running calibration jobs locally instead of on the farm [-clocal]\n";
}

###############
# Preliminaries
###############
my $cmd;              # a command to run with runCommand()
my @early_cmd_A = (); # array of commands we run before our log file is opened
# check if the $dir exists, and that it contains the files we need
# check if our output dir $symbol exists
if(! defined $dir) { 
  $dir = $ref_accn;
}
else { 
  if($dir !~ m/\/$/) { $dir =~ s/\/$//; } # remove final '/' if it exists
}
if(-d $dir) { 
  $cmd = "rm -rf $dir";
  if($do_force) { runCommand($cmd, $be_verbose, undef); push(@early_cmd_A, $cmd); }
  else          { die "ERROR directory named $dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $dir) { 
  $cmd = "rm $dir";
  if($do_force) { runCommand($cmd, $be_verbose, undef); push(@early_cmd_A, $cmd); }
  else          { die "ERROR a file named $dir already exists. Remove it, or use -f to overwrite it."; }
}

# create the dir
$cmd = "mkdir $dir";
runCommand($cmd, $be_verbose, undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root = $dir . "/" . $dir_tail . ".dnaorg_build";

# output banner
my $synopsis = "dnaorg_build.pl :: build homology models for features of a reference sequence";
my $command  = "$executable $opts_used_short $ref_accn";
my $date     = scalar localtime();
outputBanner(*STDOUT, $synopsis, $command, $date, $opts_used_long);

# open the log and command files:
# set output file names and file handles, and open those file handles
my %ofile_name_H          = (); # full name for (full path to) output files
my %ofile_sname_H         = (); # short name for (no dir path) output files
my %ofile_FH_H            = (); # file handle for output files, keys are in @ofile_keys_A
my @ofile_keys_A          = ("log", "cmd"); 
my %ofile_desc_H          = (); # description of each output file
$ofile_desc_H{"log"}      = "Output printed to screen";
$ofile_desc_H{"cmd"}      = "List of executed commands";

foreach my $key (@ofile_keys_A) { 
  $ofile_name_H{$key}  = $out_root . "." . $key;
  $ofile_sname_H{$key} = $dir_tail . "." . $key; # short name
  if(! open($ofile_FH_H{$key}, ">", $ofile_name_H{$key})) { 
    printf STDERR ("ERROR, unable to open $ofile_name_H{$key} for writing.\n"); 
    exit(1);
  }
}
my $log_FH = $ofile_FH_H{"log"};
my $cmd_FH = $ofile_FH_H{"cmd"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# now we have the log file open, output the banner there too
outputBanner($log_FH, $synopsis, $command, $date, $opts_used_long);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

############################################
# parse the optional input files, if nec
###########################################
# -matpept <f>
my @cds2pmatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of primary matpept indices that comprise this CDS
my @cds2amatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of all     matpept indices that comprise this CDS
if(defined $matpept_infile) { 
  parseMatPeptSpecFile($matpept_infile, \@cds2pmatpept_AA, \@cds2amatpept_AA, \%ofile_FH_H);
}

#####################################################
# make sure we have the executables that we require
#####################################################
my %execs_H = (); # hash with paths to all required executables
my @reqd_inf_execs_A = ("cmbuild", "cmcalibrate", "cmfetch", "cmpress");
my $reqd_exec;
foreach $reqd_exec (@reqd_inf_execs_A) { 
  $execs_H{$reqd_exec} = $inf_exec_dir . $reqd_exec;
}
my @reqd_other_execs_A = ("esl_fetch_cds");
$execs_H{"esl_fetch_cds"} = $esl_fetch_cds;
foreach $reqd_exec ((@reqd_inf_execs_A), (@reqd_other_execs_A)) { 
  if(! -x $execs_H{$reqd_exec}) { die "ERROR executable file $execs_H{$reqd_exec} does not exist (or is not executable)"; }
}

#########################################################
# generate the data files we need for building the models
#########################################################
# create the edirect .mat_peptide file, if nec, we do this first, because we will
# die with an error if mature peptide info exists and neither -matpept nor -nomatpept
# was used (and we want to die as early as possible in the script to save the user's time)


my $nsecs;
my $mp_file = $out_root . ".matpept";
my $mp_cmd = "esearch -db nuccore -query $ref_accn | efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product";
if($do_matpept) { 
  $cmd = $mp_cmd . "> $mp_file";
  $nsecs = runCommand($cmd, $be_verbose, \%ofile_FH_H);
}
else { 
  # user did not specify -matpept
  if(! $do_nomatpept) { 
    # user also did not specify -nomatpept, ensure there are no mature peptides in the annotation,
    # if there are, we die and tell the user to use -nomatpept or -matpept
    my $mp_output = `$mp_cmd`;
    if($mp_output =~ m/\w/) { 
      die "ERROR $ref_accn has mat_peptide information in its annotation.\nEither use -matpept <f> with input file <f> to model it, or -nomatpept to ignore it.";
    }
  }
}
  
# create a file with total lengths of each accession
my $len_file  = $out_root . ".length";
my $len_file_created = $len_file . ".created";
my $len_file_lost    = $len_file . ".lost";
$cmd = "esearch -db nuccore -query $ref_accn | efetch -format gpc | xtract -insd INSDSeq_length | grep . | sort > $len_file";
$nsecs = runCommand($cmd, $be_verbose, \%ofile_FH_H);

# create the edirect ftable file
my $ft_file  = $out_root . ".ftable";
$cmd = "esearch -db nuccore -query $ref_accn | efetch -format ft > $ft_file";
$nsecs = runCommand($cmd, $be_verbose, \%ofile_FH_H);

###########################################################################################
# parse those edirect output data files that we just created to get their data into usable data structures
###########################################################################################
# This is a bit convoluted because we actually do two iterations. First, we parse the
# edirect files into usable data structures. Then we write new output files in our
# standardized tabular format. Then we parse *those* standardized format files into
# our final data structures that we'll use in the remainder of the script.

# parse the length file
my %totlen_H = ();
parseLengthFile($len_file, \%totlen_H, \%ofile_FH_H);
if(! exists $totlen_H{$ref_accn}) { 
  die "ERROR, problem fetching length of reference accession $ref_accn, not in $len_file created by cmd $cmd"; 
}
my $ref_totlen = $totlen_H{$ref_accn};
printf("ref_totlen: $ref_totlen\n");
# parse the edirect ftable and (optionally) mat_peptide file
my %cds_tbl_HHA = ();   # CDS data from .cds.tbl file
                        # hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %mp_tbl_HHA = ();    # mat_peptide data from .matpept.tbl file
                        # hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column

edirectFtableOrMatPept2SingleFeatureTableInfo($ft_file, 0, "CDS", \%cds_tbl_HHA, \%ofile_FH_H); # 0: it's not a mat_peptide file
if($do_matpept) {  
  edirectFtableOrMatPept2SingleFeatureTableInfo($mp_file, 1, "mat_peptide", \%mp_tbl_HHA, \%ofile_FH_H); # 1: it is a mat_peptide file

  # now that we have the mat_peptide tabular data,
  # validate the CDS:mat_peptide relationships that we read from $matpept_infile
  matpeptValidateCdsRelationships(\@cds2pmatpept_AA, \%{$cds_tbl_HHA{$ref_accn}}, \%{$mp_tbl_HHA{$ref_accn}}, $do_circular, $totlen_H{$ref_accn}, \%ofile_FH_H);
}

########################################################
# Gather information and sequence data on the reference.
# Use each reference CDS exon or mat_peptide (if $do_matpept)
# as a homology search model against all the genomes.
#######################################################
# fetch the reference genome
my $fetch_file = $out_root . ".ref.fg.idfetch.in";
my $fasta_file = $out_root . ".ref.fg.fa";
my @accn_A = ($ref_accn);
my @seq_accn_A = (); # filled below
fetchSequencesUsingEslFetchCds($execs_H{"esl_fetch_cds"}, 0, $fetch_file, $fasta_file, $do_circular, \@accn_A, \%totlen_H, \@seq_accn_A, undef, undef, \%ofile_FH_H);

my %ftr_info_HA = (); # hash of arrays, hash keys: "type", "type_idx", "exon_idx", "nexon_map", "first_mdl", "final_mdl", values are arrays [0..$nftr-1];
                      # TODO: make this list complete
if(! exists ($cds_tbl_HHA{$ref_accn})) { 
  die "ERROR no CDS information stored for reference accession"; 
}
if($do_matpept && (! exists ($mp_tbl_HHA{$ref_accn}))) { 
  die "ERROR -matpept enabled, but no mature peptide information stored for reference accession"; 
}

my $ref_mp_strand_str = "";
my $ref_cds_strand_str = "";
if($do_matpept) { 
  (undef, undef, undef, undef, undef, $ref_mp_strand_str) = getStrandStats(\%mp_tbl_HHA, $ref_accn, \%ofile_FH_H);
}
if(%cds_tbl_HHA) { 
  (undef, undef, undef, undef, undef, $ref_cds_strand_str) = getStrandStats(\%cds_tbl_HHA, $ref_accn, \%ofile_FH_H);
}
my $ref_strand_str = $ref_mp_strand_str . $ref_cds_strand_str;

@{$ftr_info_HA{"ref_strand"}} = split("", $ref_strand_str);

my @ref_ftr_len_A     = (); # [0..$i..$nftr-1]: length of each reference main feature (CDS or mat_peptide)
my @ref_ftr_coords_A  = (); # [0..$i..$nftr-1]: main feature coords (CDS or mat_peptide) for reference
@{$ftr_info_HA{"ref_len"}}     = ();
@{$ftr_info_HA{"ref_coords"}}  = ();
@{$ftr_info_HA{"out_product"}} = ();
if($do_matpept) { 
  getLengthsAndCoords(\%{$mp_tbl_HHA{$ref_accn}}, \@{$ftr_info_HA{"ref_len"}}, \@{$ftr_info_HA{"ref_coords"}}, \%ofile_FH_H);
  getQualifierValues(\%mp_tbl_HHA, $ref_accn, "product", \@{$ftr_info_HA{"out_product"}}, \%ofile_FH_H);
}
if(%cds_tbl_HHA) { 
  getLengthsAndCoords(\%{$cds_tbl_HHA{$ref_accn}}, \@{$ftr_info_HA{"ref_len"}}, \@{$ftr_info_HA{"ref_coords"}}, \%ofile_FH_H);
  getQualifierValues(\%cds_tbl_HHA, $ref_accn, "product", \@{$ftr_info_HA{"out_product"}}, \%ofile_FH_H);
}

my $nftr; # number of features
my $nmp;  # number of mature peptide features
$nftr = validateAndGetSizeOfInfoHashOfArrays(\%ftr_info_HA, undef, \%ofile_FH_H);
$nmp  = ($do_matpept) ? scalar(@{$mp_tbl_HHA{$ref_accn}{"coords"}}) : 0;

# determine types of each feature 
determineFeatureTypes($nmp, ((@cds2pmatpept_AA) ? \@cds2pmatpept_AA : undef), \%ftr_info_HA, \%ofile_FH_H);

# fetch the reference feature sequences and populate information on the models and features
my $all_stk_file = $out_root . ".ref.all.stk";
my %mdl_info_HA = (); # hash of arrays, hash keys: "ftr_idx",  "is_first",  "is_final",  values are arrays [0..$nmdl-1];
my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $fasta_file });
fetchReferenceFeatureSequences($sqfile, $seq_accn_A[0], $ref_totlen, $out_root, $do_circular, \%mdl_info_HA, \%ftr_info_HA, $all_stk_file, \%ofile_FH_H); # 0 is 'do_circular' which is irrelevant in this context

dumpInfoHashOfArrays("ftr_info", 1, \%ftr_info_HA, *STDOUT);
dumpInfoHashOfArrays("mdl_info", 1, \%mdl_info_HA, *STDOUT);
validateAndGetSizeOfInfoHashOfArrays(\%ftr_info_HA, undef, \%ofile_FH_H);
validateAndGetSizeOfInfoHashOfArrays(\%mdl_info_HA, undef, \%ofile_FH_H);

##################
# Build the models
##################
#createCmDb(\%execs_H, $do_cslow, $do_clocal, $all_stk_file, $out_root . ".ref", \@{$mdl_info_HA{"cmname"}}, \%ofile_FH_H);
printf("#\n# Model calibration %s. Exiting.\n", ($do_clocal) ? "complete" : "job(s) submitted");

exit 0;


