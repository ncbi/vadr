#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Thu Feb 18 12:48:16 2016 [dnaorg_annotate.pl split off from dnaorg_annotate_genomes.pl]
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
# Given a set of homology models built for a reference accession
# created by a previous run of dnaorg_build.pl, this script uses those
# homology models to annotate the features (CDS and mature peptides)
# in other accessions, provided as a list as the lone command line
# argument. This script outputs tabular annotations and a list of
# errors found. 
#
# A list of subroutines can be found after the main script before
# the subroutines, on or about line 1000.
# 
# Immediately below are a list of the steps performed by this
# script. Each has a corresponding code block in the main script.
#
# Preliminaries: 
#   - process options
#   - output program banner and open output files
#   - parse the optional input files, if necessary
#   - make sure the required executables are executable
#
# Step 1.  Gather and process information on reference genome using Edirect.
# Step 2.  Fetch and process the reference genome sequence.
# Step 3.  Verify we have the model file that we need to do homology searches.
#
#    Steps 1-2 are very similar to those done in dnaorg_build.pl, but
#    dnaorg_build.pl and dnaorg_annotate,pl are run separately, perhaps
#    with a big time lag. Therefore, this script, in step 3, verifies 
#    that the models built by dnaorg_build.pl are still up-to-date with 
#    respect to the reference used in running this script.
#    The reference sequence is required to be the first accession 
#    listed in the accession list file.
#
# Step 4.  (OPTIONAL) Search for origin sequences, if --origin used.
#          The term "origin sequences" refers to the sequence around the position
#          numbered 1 in circular genomes. 
# Step 5.  Perform homology searches.
# Step 6.  Parse homology search results into usable data structures.
# Step 7.  Define names of all per-model and per-feature output files 
#          we're about to create.
# Step 8.  Fetch predicted matches into fasta files.
# Step 9.  For features modeled by multiple models (e.g. multi-exon CDS),
#          combine all relevant matches into predicted feature sequences.
# Step 10. For features modeled by multiple other features (e.g. CDS comprised
#          of multiple mature peptides), cocatenate the individual predicted feature
#          sequences into single sequences.
# Step 11. Examine each predicted feature to see if:
#          - start codon is valid
#          - where first in-frame stop exists (for all features with valid start codon)
#          Store the relevant information in @err_ftr_instances_AHH to either output later
#          or examine more later.
# Step 12. For all features that have a valid start but no in-frame
#          stop codon, look for an in-frame stop later in the sequence (downstream).
# Step 13. For all features with either a 'trc' or 'ext'
#          error, add a new stop position to the results_AAH
#          data structure, the 'corrected' stop. We utilize
#          the 'cumlen' values in the @mdl_results_AAH to 
#          simplify this for multi-exon features. 
# Step 14. Identify overlap and adjacency errors.
# Step 15. Finalize the mdl_results, and fill ftr_results and check
#          for incompatible error combinations.
# Step 16. Refetch corrected matches into new files.
# Step 17. For features modeled by multiple models (e.g. multi-exon CDS),
#          combine all relevant matches into corrected feature sequences.
# Step 18. For features modelled by multiple other features (e.g. CDS
#          comprised of multiple mature peptides) concatenate the
#          individual corrected feature sequences into single
#          sequences.
# Step 19. Translate corrected feature sequences into putative proteins.
# Step 20. (OPTIONAL) Create multiple alignments of DNA sequences, if --doalign.
# Step 21. (OPTIONAL) Create multiple alignments of protein sequences, if --doalign.
# Step 22. Output annotations and errors.
#
#######################################################################################
#
# Error identification:
#
# This script identifies and outputs a list of all errors in each of
# the accessions it annotates. Each type of error has an associated
# three letter error 'code'. The list of error codes is in the
# dnaorg.pm perl module file. It can also be found here:
#
#  /panfs/pan1/dnaorg/virseqannot/error_code_documentation/errorcodes.v5d.documentation.txt
#
# The most complicated code in this script is related to identifying
# and storing these errors. Ideally, there would be one or a small set
# of functions that identified the errors. However, the current
# implementation has the identification of errors spread across 7
# functions. This is partly due to a imperfect design of the code, but
# also due to the requirement of identifying different types of errors
# at different stages of the analysis.
#
# For example, it is necessary to identify in-frame stops in the
# predicted hits (potential trc errors) before we can finalize our
# annotations (since the homology search software (Infernal) is not
# looking for start and stop codons). This means we need to identify
# potential trc errors early, and the identification of these errors
# will affect later errors, for example 'nm3' errors rely on the final
# annotation lengths and so cannot be computed until after all the
# 'trc' errors have been found.
#
# Another example is that CDS features that are comprised of mature
# peptides (type: cds-mp) cannot be examined for errors until all of
# the mature peptide annotations that comprise them are finalized.
# 
# The following table shows which errors are identified in which
# functions. There are two annotation types (annot_type) of features,
# which is why there are two columns under 'annot_type', namely the
# 'model' features which are annotated based on homology model
# predictions (type: 'mp' and 'cds-notmp') and 'multifeature' features
# which are based on the annotations of multiple other features (type:
# 'cds-mp'). The identification of errors for these two feature
# annotation types is done differently, and is often done in different
# functions.
# 
# The sequence column is only relevant for errors that are
# per-sequence, instead of per-feature. Currently the only
# per-sequence error is the ori error.
#
# List of functions in which errors may be detected: 
# 1. parse_esl_epn_translate_startstop_outfile()
# 2. results_calculate_corrected_stops()
# 3. results_calculate_overlaps_and_adjacencies() 
# 4. mdl_results_add_stp_nop_ost_b3e_b3u_errors()
# 5. ftr_results_calculate()
# 6. find_origin_sequences()
# 7. MAIN (not a function but rather the main body of the script):
# 8. mdl_results_add_b5e_b5u_errors()
#
#              annot_type
#          -------------------
# errcode  model  multifeature sequence
# -------  -----  ------------ --------
# nop      4      N/A          N/A
# nm3      5      5            N/A
# b5e      8      N/A          N/A
# b5u      8      N/A          N/A
# b3e      4      N/A          N/A
# b3u      4      N/A          N/A
# olp      3      N/A          N/A
# str      1,4    1,5          N/A
# stp      1,2    1,5          N/A
# ajb      3      N/A          N/A
# aja      3      N/A          N/A                           
# trc      1,2    1,2          N/A
# ext      1,7,2  1,7,2        N/A
# ntr      5      N/A          N/A
# nst      1,7    1,7          N/A
# ost      4      N/A          N/A
# aji      N/A    5            N/A
# int      N/A    5            N/A
# inp      N/A    5            N/A
# ori      N/A    N/A          6
# 
#######################################################################################

# first, determine the paths to all modules, scripts and executables that we'll need
# we currently use hard-coded-paths for Infernal, HMMER and easel executables:
my $inf_exec_dir      = "/usr/local/infernal/1.1.1/bin/";
my $hmmer_exec_dir    = "/usr/local/hmmer/3.1/bin/";
my $esl_exec_dir      = "/usr/local/infernal/1.1.1/bin/";

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
my $esl_fetch_cds     = $dnaorgdir . "/esl-fetch-cds/esl-fetch-cds.pl";
my $esl_epn_translate = $dnaorgdir . "/esl-epn-translate/esl-epn-translate.pl";
my $esl_ssplit        = $dnaorgdir . "/Bio-Easel/scripts/esl-ssplit.pl";

#########################################################
# Command line and option processing using epn-options.pm
#
# opt_HH: 2D hash:
#         1D key: option name (e.g. "-h")
#         2D key: string denoting type of information 
#                 (one of "type", "default", "group", "requires", "incompatible", "preamble", "help")
#         value:  string explaining 2D key:
#                 "type":          "boolean", "string", "integer" or "real"
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
#     option            type       default               group   requires incompat    preamble-output                                 help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                          "display this help",                                  \%opt_HH, \@opt_order_A);
opt_Add("-c",           "boolean", 0,                        1,    undef, undef,      "genome is closed (a.k.a. circular)",           "genome is closed (a.k.a circular)",                  \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,                        1,"--dirout",undef,      "forcing directory overwrite (with --dirout)",  "force; if dir from --dirout exists, overwrite it",   \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                                   "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--dirout",     "string",  undef,                    1,    undef, undef,   "output directory specified as",                "specify output directory as <s>, not <ref accession>", \%opt_HH, \@opt_order_A);
opt_Add("--dirbuild",   "string",  undef,                    1,"--dirout",   undef,   "output directory used for dnaorg_build.pl",    "specify output directory used for dnaorg_build.pl as <s> (created with dnaorg_build.pl --dirout <s>), not <ref accession>", \%opt_HH, \@opt_order_A);
opt_Add("--origin",     "string",  undef,                    1,     "-c", undef,      "identify origin seq <s> in genomes",           "identify origin seq <s> in genomes, put \"|\" at site of origin (\"|\" must be escaped, i.e. \"\\|\"", \%opt_HH, \@opt_order_A);
opt_Add("--matpept",    "string",  undef,                    1,    undef, undef,      "using pre-specified mat_peptide info",         "read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships", \%opt_HH, \@opt_order_A);
opt_Add("--nomatpept",  "boolean", 0,                        1,    undef,"--matpept", "ignore mat_peptide annotation",                "ignore mat_peptide information in reference annotation", \%opt_HH, \@opt_order_A);
opt_Add("--specstart",  "string",  undef,                    1,    undef, undef,      "using pre-specified alternate start codons",   "read specified alternate start codons per CDS from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",           "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);
opt_Add("--local",      "boolean", 0,                        1,    undef, undef,      "run cmscan locally instead of on farm",        "run cmscan locally instead of on farm", \%opt_HH, \@opt_order_A);
opt_Add("--errcheck",   "boolean", 0,                        1,    undef,"--local",   "consider any farm stderr output as indicating a job failure", "consider any farm stderr output as indicating a job failure", \%opt_HH, \@opt_order_A);
opt_Add("--nseq",       "integer", 5,                        1,    undef,"--local",   "number of sequences for each cmscan farm job", "set number of sequences for each cmscan farm job to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--wait",       "integer", 500,                      1,    undef,"--local",   "allow <n> minutes for cmscan jobs on farm",    "allow <n> wall-clock minutes for cmscan jobs on farm to finish, including queueing time", \%opt_HH, \@opt_order_A);
opt_Add("--bigthresh",  "integer", 4000,                     1,    undef, undef,      "set minimum model length for using HMM mode to <n>", "set minimum model length for using HMM mode to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--midthresh",  "integer", 100,                      1,    undef, undef,      "set max model length for using mid sensitivity mode to <n>", "set max model length for using mid sensitivity mode to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--smallthresh","integer", 30,                       1,    undef, undef,      "set max model length for using max sensitivity mode to <n>", "set max model length for using max sensitivity mode to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--mxsize",     "integer", 2048,                     1,"--doalign",undef,     "with --doalign, set --mxsize <n> to <n>",      "with --doalign, set --mxsize <n> for cmalign to <n>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"2"} = "options for alternative modes";
 #       option               type   default                group  requires incompat                        preamble-output                                                      help-output    
opt_Add("--infasta",     "boolean", 0,                       2,"--refaccn", "--skipedirect,--skipfetch",   "single cmdline argument is a fasta file of sequences, not a list of accessions", "single cmdline argument is a fasta file of sequences, not a list of accessions", \%opt_HH, \@opt_order_A);
opt_Add("--refaccn",     "string",  undef,                   2,"--infasta", "--skipedirect,--skipfetch",   "specify reference accession is <s>",                                "specify reference accession is <s> (must be used in combination with --infasta)", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"3"} = "options that modify the tabular output file";
#       option               type   default                group  requires incompat preamble-output                                               help-output    
opt_Add("--tblfirst",    "boolean", 0,                      3,    undef,   undef,   "put first accession first on each .tbl page",               "include annotation for first accession on each page of .tbl output file", \%opt_HH, \@opt_order_A);
opt_Add("--tblnocomp",   "boolean", 0,                      3,    undef,   undef,   "do not compare annotations to existing GenBank annotation", "do not include information comparing predicted annotations to existing GenBank annotations", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"4"} = "options for skipping/adding optional stages";
#       option               type   default                group  requires incompat preamble-output                             help-output    
opt_Add("--doalign",    "boolean", 0,                       4,    undef,   undef,   "create nucleotide and protein alignments", "create nucleotide and protein alignments", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"5"} = "optional output files";
#       option       type       default                  group  requires incompat  preamble-output                          help-output    
opt_Add("--mdlinfo",    "boolean", 0,                        5,    undef, undef, "output internal model information",     "create file with internal model information",   \%opt_HH, \@opt_order_A);
opt_Add("--ftrinfo",    "boolean", 0,                        5,    undef, undef, "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);
opt_Add("--seqinfo",    "boolean", 0,                        5,    undef, undef, "output internal sequence information",  "create file with internal sequence information", \%opt_HH, \@opt_order_A);
opt_Add("--errinfo",    "boolean", 0,                        5,    undef, undef, "output internal error information",     "create file with internal error information", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{"6"} = "options for skipping stages and using files from earlier, identical runs, primarily useful for debugging";
#     option               type       default               group   requires    incompat                    preamble-output                                            help-output    
opt_Add("--skipedirect",   "boolean", 0,                       6,   undef,      "-f,--nseq,--local,--wait", "skip the edirect steps, use existing results",           "skip the edirect steps, use data from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipfetch",     "boolean", 0,                       6,   undef,      "-f,--nseq,--local,--wait", "skip the sequence fetching steps, use existing results", "skip the sequence fetching steps, use files from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipscan",      "boolean", 0,                       6,   undef,      "-f,--nseq,--local,--wait", "skip the cmscan step, use existing results",             "skip the cmscan step, use results from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skiptranslate", "boolean", 0,                       6,"--skipscan",  undef,                      "skip the translation steps, use existing resutls",       "skip the translation steps, use results from an earlier run of the script", \%opt_HH, \@opt_order_A);


$opt_group_desc_H{"7"} = "TEMPORARY options for the alternative method of identifying origin sequences";
#     option               type       default               group   requires                                   incompat     preamble-output                                                         help-output    
opt_Add("--aorgmodel",     "string",  undef,                   7,   "-c,--aorgstart,--aorgoffset,--aorglen",   "--origin",  "use alternative origin method with model <s>",                         "use alternative origin method with origin model in <s>", \%opt_HH, \@opt_order_A);
opt_Add("--aorgstart",     "integer", 0,                       7,   "-c,--aorgmodel,--aorgoffset,--aorglen",   "--origin",  "origin begins at position <n> in --aorgmodel model",                   "origin begins at position <n> in --aorgmodel model",     \%opt_HH, \@opt_order_A);
opt_Add("--aorgoffset",    "integer", 0,                       7,   "-c,--aorgmodel,--aorgstart,--aorglen",    "--origin",  "first position of genome sequence is position <n> in origin sequence", "first position of genome sequence is position <n> in origin sequence", \%opt_HH, \@opt_order_A);
opt_Add("--aorglen",       "integer", 0,                       7,   "-c,--aorgmodel,--aorgstart,--aorgoffset", "--origin",  "length of origin sequence is <n>",                                     "length of origin sequence is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--aorgethresh",   "real",    1.0,                     7,   "-c,--aorgmodel,--aorgstart,--aorgoffset", "--origin",  "E-value threshold for origin detection is <x>",                        "E-value threshold for origin detection is <x>", \%opt_HH, \@opt_order_A);
opt_Add("--aorgppthresh",  "real",    0.6,                     7,   "-c,--aorgmodel,--aorgstart,--aorgoffset", "--origin",  "average PP threshold for origin detection is <x>",                     "average PP threshold for origin detection is <x>", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_annotate.pl [-options] <file with list of accessions to annotate>\n";
$usage      .= "\n";
$usage      .= "       OR\n";
$usage      .= "\n";
$usage      .= "       dnaorg_annotate.pl [-options] --refaccn <reference accession> --infasta <fasta sequence file with sequences to annotate>\n";
my $synopsis = "dnaorg_annotate.pl :: annotate sequences based on a reference annotation";
my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'c'            => \$GetOptions_H{"-c"},
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'dirout=s'     => \$GetOptions_H{"--dirout"},
                'dirbuild=s'   => \$GetOptions_H{"--dirbuild"},
                'origin=s'     => \$GetOptions_H{"--origin"},
                'matpept=s'    => \$GetOptions_H{"--matpept"},
                'nomatpept'    => \$GetOptions_H{"--nomatpept"},
                'specstart=s'  => \$GetOptions_H{"--specstart"},
                'keep'         => \$GetOptions_H{"--keep"},
                'local'        => \$GetOptions_H{"--local"}, 
                'errcheck'     => \$GetOptions_H{"--errcheck"},  
                'nseq=s'       => \$GetOptions_H{"--nseq"}, 
                'wait=s'       => \$GetOptions_H{"--wait"},
                'bigthresh=s'  => \$GetOptions_H{"--bigthresh"},
                'midthresh=s'  => \$GetOptions_H{"--midthresh"},
                'smallthresh=s'=> \$GetOptions_H{"--smallthresh"},
                'mxsize=s'     => \$GetOptions_H{"--mxsize"},
# options for alternative modes
                'infasta'      => \$GetOptions_H{"--infasta"},
                'refaccn=s'    => \$GetOptions_H{"--refaccn"},
# options that affect tabular output file
                'tblfirst'     => \$GetOptions_H{"--tblfirst"},
                'tblnocomp'    => \$GetOptions_H{"--tblnocomp"},
# options for skipping/adding optional stages
                'doalign'      => \$GetOptions_H{"--doalign"},
# optional output files
                'mdlinfo'      => \$GetOptions_H{"--mdlinfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"}, 
                'seqinfo'      => \$GetOptions_H{"--seqinfo"}, 
                'errinfo'      => \$GetOptions_H{"--errinfo"},
# options for skipping stages, using earlier results
                'skipedirect'   => \$GetOptions_H{"--skipedirect"},
                'skipfetch'     => \$GetOptions_H{"--skipfetch"},
                'skipscan'      => \$GetOptions_H{"--skipscan"},
                'skiptranslate' => \$GetOptions_H{"--skiptranslate"}, 
# options for alternative origin detection method
                'aorgmodel=s'   => \$GetOptions_H{"--aorgmodel"},
                'aorgstart=s'   => \$GetOptions_H{"--aorgstart"},
                'aorglen=s'     => \$GetOptions_H{"--aorglen"},
                'aorgoffset=s'  => \$GetOptions_H{"--aorgoffset"}, 
                'aorgethresh=s' => \$GetOptions_H{"--aorgethresh"}, 
                'aorgppthresh=s'=> \$GetOptions_H{"--aorgppthresh"});

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.15";
my $releasedate   = "Oct 2016";

# make *STDOUT file handle 'hot' so it automatically flushes whenever we print to it
# it is printed to
select *STDOUT;
$| = 1;

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 1) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do dnaorg_annotate.pl -h\n\n";
  exit(1);
}
my ($listfile) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# validate origin sequence if necessary
my $origin_offset = undef;
my $origin_seq = undef;
if(opt_IsUsed("--origin", \%opt_HH)) { 
  $origin_seq = opt_Get("--origin", \%opt_HH);
  $origin_seq =~ tr/a-z/A-Z/; # capitalize origin seq
  $origin_offset = validate_origin_seq($origin_seq);
  $origin_seq =~ s/\|//; # remove the single "|"
}

# validate alternative origin detection options if necessary
if(opt_IsUsed("--aorgmodel", \%opt_HH)) { 
  my $aorg_model  = opt_Get("--aorgmodel", \%opt_HH);
  if(! -s $aorg_model) { 
    die "ERROR with --aorgmodel <s>, $aorg_model does not exist or is empty"; 
  }    
  my $aorg_offset = opt_Get("--aorgoffset", \%opt_HH);
  my $aorg_len    = opt_Get("--aorglen",  \%opt_HH);
  if(($aorg_offset <= 0) || ($aorg_offset > $aorg_len)) { 
    die "ERROR with --aorgoffset <n>, <n> must be greater than 0 and less than or equal to <n> from --aorglen";
  }
}

my $dir_build  = opt_Get("--dirbuild", \%opt_HH);  # this will be undefined unless --dirbuild set on cmdline
my $dir_out    = opt_Get("--dirout",   \%opt_HH);  # this will be undefined unless --dirout set on cmdline
my $do_matpept = opt_IsOn("--matpept", \%opt_HH);  # this will be '0' unless --matpept set on cmdline 

if(defined $dir_out) { 
  $dir_out =~ s/\/$//; # remove final '/' if there is one
}
if(defined $dir_build) { 
  $dir_build =~ s/\/$//; # remove final '/' if there is one
}

# if --infasta used, $listfile is actually a fasta file
my $infasta_file = undef;
my $do_infasta = 0;
if(opt_Get("--infasta", \%opt_HH)) { 
  $infasta_file = $listfile;
  $listfile = undef;
  $do_infasta = 1;
}

# if --smallthresh or --midthresh or --bigthresh used, validate that the thresholds make sense:
# small < mid < big
if(opt_IsUsed("--smallthresh", \%opt_HH)) { 
  if(opt_Get("--smallthresh", \%opt_HH) >= opt_Get("--midthresh", \%opt_HH)) { 
    die sprintf("ERROR, with --smallthresh <x> and --midthresh <y>, <x> must be < <y>, got <x> = %d and <y> = %d\n", 
                opt_Get("--smallthresh", \%opt_HH), opt_Get("--midthresh", \%opt_HH));
  }
  if(opt_Get("--smallthresh", \%opt_HH) >= opt_Get("--bigthresh", \%opt_HH)) { 
    die sprintf("ERROR, with --smallthresh <x> and --bigthresh <y>, <x> must be < <y>, got <x> = %d and <y> = %d\n", 
                opt_Get("--smallthresh", \%opt_HH), opt_Get("--bigthresh", \%opt_HH));
  }
}
if(opt_IsUsed("--midthresh", \%opt_HH)) { 
  if(opt_Get("--midthresh", \%opt_HH) < opt_Get("--smallthresh", \%opt_HH)) { 
    die sprintf("ERROR, with --smallthresh <x> and --midthresh <y>, <x> must be < <y>, got <x> = %d and <y> = %d\n", 
                opt_Get("--smallthresh", \%opt_HH), opt_Get("--midthresh", \%opt_HH));
  }
  if(opt_Get("--midthresh", \%opt_HH) >= opt_Get("--bigthresh", \%opt_HH)) { 
    die sprintf("ERROR, with --midthresh <x> and --bigthresh <y>, <x> must be < <y>, got <x> = %d and <y> = %d\n", 
                opt_Get("--midthresh", \%opt_HH), opt_Get("--bigthresh", \%opt_HH));
  }
}
if(opt_IsUsed("--bigthresh", \%opt_HH)) { 
  if(opt_Get("--bigthresh", \%opt_HH) < opt_Get("--smallthresh", \%opt_HH)) { 
    die sprintf("ERROR, with --smallthresh <x> and --bigthresh <y>, <x> must be < <y>, got <x> = %d and <y> = %d\n", 
                opt_Get("--smallthresh", \%opt_HH), opt_Get("--bigthresh", \%opt_HH));
  }
  if(opt_Get("--bigthresh", \%opt_HH) < opt_Get("--midthresh", \%opt_HH)) { 
    die sprintf("ERROR, with --midthresh <x> and --bigthresh <y>, <x> must be < <y>, got <x> = %d and <y> = %d\n", 
                opt_Get("--midthresh", \%opt_HH), opt_Get("--bigthresh", \%opt_HH));
  }
}

###############
# Preliminaries
###############
# first, parse the list file, we need to do this first because we need
# to know what the reference accession <refaccn> is to check if the
# directory <refaccn> exists
my $cmd;               # a command to run with runCommand()
my @early_cmd_A = ();  # array of commands we run before our log file is opened
my %seq_info_HA = ();  # hash of arrays, values are arrays with index range [0..$nseq-1];
                       # 1st dim keys are "seq_name", "accn_name", "seq_len", "accn_len".
                       # $seq_info_HA{"accn_name"}[0] is our reference accession
@{$seq_info_HA{"accn_name"}} = ();

my %infasta_ref_seq_info_HA = ();  # hash of arrays, for reference sequence information. 
                                   # only used if --infasta used. Actually only stores information
                                   # on 1 sequence, so could be just a hash, but it is a hash of 
                                   # single element arrays so that it is the same type of data
                                   # structure as %seq_info_HA so we can pass it into 
                                   # functions (namely wrapperGetInfoUsingEdirect) in place
                                   # of %seq_info_HA.
                                   # 1st dim keys are "seq_name", "accn_name", "seq_len", "accn_len".
                                   # $infasta_ref_seq_info_HA{"accn_name"}[0] is our reference accession

my $nseq = 0;
my $ref_accn = undef;
if(! defined $infasta_file) { # default mode
  parseListFile($listfile, 1, $seq_info_HA{"accn_name"}, undef); # 1 
  $nseq = scalar(@{$seq_info_HA{"accn_name"}});
  $ref_accn = $seq_info_HA{"accn_name"}[0];
}
else { # --infasta used
  if(! (opt_Get("--refaccn", \%opt_HH))) { 
    # we should never get here, because epn-options.pm should have enforced that
    # --refaccn and --infasta were both used if one was, but we do a second check
    # here just to make sure
    DNAORG_FAIL("ERROR, --infasta requires --refaccn", 1, undef);
  }
  else { 
    $ref_accn = opt_Get("--refaccn", \%opt_HH);
    stripVersion(\$ref_accn);
    # initialize the sequence info hash of arrays
  }
}

# determine the directory in which dnaorg_build files are in ($dir_build)
# it will already be defined if --dirbuild was used
if(! defined $dir_build) { 
  $dir_build = $ref_accn;
  if(! -d $dir_build) {
    DNAORG_FAIL(sprintf("ERROR, directory $dir_build (%s) does not exist.\nDid you run \"dnaorg_build.pl $dir_build\" yet? If not, you need to do that first.", opt_IsUsed("--refaccn", \%opt_HH) ? "specified with --refaccn" : "first accession from in $listfile"), 1, undef);
  }
}
else { # --dirbuild was used on the command line
  # currently --dirbuild requires --dirout (enforced by epn-options.pm) but we double
  # check that if --dirbuild is on, then --dirout must be too here (just to be safe)
  if(! defined $dir_out) { 
    DNAORG_FAIL("ERROR the --dirbuild option requires the --dirout option be used also", 1, undef);
  }
  if(! -d $dir_build) {
    DNAORG_FAIL("ERROR, directory $dir_build (specified with --dirbuild) does not exist.\nDid you run \"dnaorg_build.pl --dirout $dir_build\" yet? If not, you need to do that first.", 1, undef);
  }
}

# determine the directory we'll put output files in,
# it will already be defined if --dirout was used
if(! defined $dir_out) { 
  $dir_out = $ref_accn;
}
else { 
  # --dirout was used to specify $dir_out
  # if it already exists (and it's not $ref_accn and it's not $dir_build) 
  # check if one of the skip options was used (begin with --skip) 
  # if so, try to use it. Else tell user to either rerun with -f
  # or delete it.
  if(($dir_out ne $ref_accn) && ($dir_out ne $dir_build)) { 
    if(-d $dir_out) { 
      $cmd = "rm -rf $dir_out";
      if(opt_Get("-f", \%opt_HH)) { # -f used, always remove it
        runCommand($cmd, opt_Get("-v", \%opt_HH), undef); push(@early_cmd_A, $cmd); 
      }
      else { # dirout exists but -f not used
        if(! ((opt_IsUsed("--skipedirect",   \%opt_HH)) || 
              (opt_IsUsed("--skipfetch",     \%opt_HH)) || 
              (opt_IsUsed("--skipscan",      \%opt_HH)) || 
              (opt_IsUsed("--skiptranslate", \%opt_HH)))) { 
          die "ERROR directory named $dir_out (specified with --dirout) already exists. Remove it, or use -f to overwrite it."; 
        }
        # if a --skip option is used, we just press on
      }
    }
    elsif(-e $dir_out) { 
      $cmd = "rm $dir_out";
      if(opt_Get("-f", \%opt_HH)) { runCommand($cmd, opt_Get("-v", \%opt_HH), undef); push(@early_cmd_A, $cmd); }
      else                        { die "ERROR a file named $dir_out (specified with --dirout) already exists. Remove it, or use -f to overwrite it."; }
    }
  }
}
# if $dir_out does not exist, create it
if(! -d $dir_out) {
  $cmd = "mkdir $dir_out";
  runCommand($cmd, opt_Get("-v", \%opt_HH), undef); push(@early_cmd_A, $cmd);
}

my $dir_out_tail   = $dir_out;
my $dir_build_tail = $dir_build;
$dir_out_tail   =~ s/^.+\///; # remove all but last dir
$dir_build_tail =~ s/^.+\///; # remove all but last dir
my $out_root   = $dir_out .   "/" . $dir_out_tail   . ".dnaorg_annotate";
my $build_root = $dir_build . "/" . $dir_build_tail . ".dnaorg_build";

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ();
my @arg_A      = ();

if(defined $listfile) { 
  push(@arg_desc_A, "file with list of accessions");
  push(@arg_A, $listfile);
}
elsif(defined $infasta_file) { 
  push(@arg_desc_A, "fasta file with sequences to annotate (--infasta)");
  push(@arg_A, $infasta_file);
}
else { 
  DNAORG_FAIL("ERROR, both listfile and infasta_file are undefined...", 1, undef);
}

outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble(*STDOUT, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# open the log and command files:
# set output file names and file handles, and open those file handles
my %ofile_info_HH = ();  # hash of information on output files we created,
                         # 1D keys: 
                         #  "fullpath":  full path to the file
                         #  "nodirpath": file name, full path minus all directories
                         #  "desc":      short description of the file
                         #  "FH":        file handle to output to for this file, maybe undef
                         # 2D keys (at least initially)
                         #  "log":  log file of what's output to stdout
                         #  "cmd":  command file with list of all commands executed
                         #  "list": file with list of all output files created

# open the log and command files 
openAndAddFileToOutputInfo(\%ofile_info_HH, "log",  $out_root . ".log",  1, "Output printed to screen");
openAndAddFileToOutputInfo(\%ofile_info_HH, "cmd",  $out_root . ".cmd",  1, "List of executed commands");
openAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".list", 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# open gap info files, if --doalign
if(opt_Get("--doalign", \%opt_HH))  { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "gap_perseq_all",     $out_root . ".gap.perseq-all.txt",     1, "Per sequence gap information (all gaps)");
  openAndAddFileToOutputInfo(\%ofile_info_HH, "gap_perseq_not3",    $out_root . ".gap.perseq-not3.txt",    1, "Per sequence gap information (gaps not a multiple of 3)");
  openAndAddFileToOutputInfo(\%ofile_info_HH, "gap_perseq_special", $out_root . ".gap.perseq-special.txt", 1, "Per sequence gap information (special gaps)");
  openAndAddFileToOutputInfo(\%ofile_info_HH, "gap_pergap_all",     $out_root . ".gap.pergap-all.txt",     1, "Per gap information (all gaps)");
  openAndAddFileToOutputInfo(\%ofile_info_HH, "gap_pergap_not3",    $out_root . ".gap.pergap-not3.txt",    1, "Per gap information (gaps not a multiple of 3)");
  openAndAddFileToOutputInfo(\%ofile_info_HH, "gap_pergap_special", $out_root . ".gap.pergap-special.txt", 1, "Per gap information (special gaps)");
}

# open optional output files
if(opt_Get("--mdlinfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "mdlinfo", $out_root . ".mdlinfo", 1, "Model information (created due to --mdlinfo)");
}
if(opt_Get("--ftrinfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "ftrinfo", $out_root . ".ftrinfo", 1, "Feature information (created due to --ftrinfo)");
}
if(opt_Get("--seqinfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "seqinfo", $out_root . ".seqinfo", 1, "Sequence information (created due to --seqinfo)");
}
if(opt_Get("--errinfo", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "errinfo", $out_root . ".errinfo", 1, "Error information (created due to --errinfo)");
}

# now we have the log file open, output the banner there too
outputBanner($log_FH, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

##############################################
# parse the optional input files, if necessary
##############################################
# --matpept <f>
my @cds2pmatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of primary matpept indices that comprise this CDS
my @cds2amatpept_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of all     matpept indices that comprise this CDS
if($do_matpept) { 
  parseMatPeptSpecFile(opt_Get("--matpept", \%opt_HH), \@cds2pmatpept_AA, \@cds2amatpept_AA, $ofile_info_HH{"FH"});
}
# --specstart <f>
my @specstart_AA = (); # 1st dim: cds index (-1, off-by-one), 2nd dim: value array of allowed start codons for this CDS
if(opt_IsOn("--specstart", \%opt_HH)) { 
  parseSpecStartFile(opt_Get("--specstart", \%opt_HH), \@specstart_AA, $ofile_info_HH{"FH"});
}

###################################################
# make sure the required executables are executable
###################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmscan"}        = $inf_exec_dir   . "cmscan";
$execs_H{"cmalign"}       = $inf_exec_dir   . "cmalign";
$execs_H{"hmmbuild"}      = $hmmer_exec_dir . "hmmbuild";
$execs_H{"hmmalign"}      = $hmmer_exec_dir . "hmmalign";
$execs_H{"esl-reformat"}  = $esl_exec_dir   . "esl-reformat";
$execs_H{"esl_fetch_cds"} = $esl_fetch_cds;
$execs_H{"esl_ssplit"}    = $esl_ssplit;
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################################################################
# Step 0. Read the dnaorg_build.pl consopts file and make sure that it
#         agrees with the options set here.
###########################################################################
my $progress_w = 85; # the width of the left hand column in our progress output, hard-coded
my $start_secs = outputProgressPrior("Verifying options are consistent with options used for dnaorg_build.pl", $progress_w, $log_FH, *STDOUT);
validate_options_are_consistent_with_dnaorg_build($build_root . ".consopts", \%opt_HH, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###########################################################################
# Step 1. Gather and process information on reference genome using Edirect.
###########################################################################
my $progress_str = sprintf("Gathering information on %d sequences using edirect", $nseq);
if(opt_Get("--skipedirect", \%opt_HH)) { 
  $progress_str = sprintf("Processing information on %d sequences fetched earlier using edirect", $nseq);
}
elsif(opt_Get("--infasta", \%opt_HH)) { 
  $progress_str = "Processing input fasta file";
}
$start_secs = outputProgressPrior($progress_str, $progress_w, $log_FH, *STDOUT);

my %cds_tbl_HHA = ();   # CDS data from .cds.tbl file, hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %mp_tbl_HHA = ();    # mat_peptide data from .matpept.tbl file, hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column


# Call the wrapper function that does the following:
#  1) creates the edirect .mat_peptide file, if necessary
#  2) creates the edirect .ftable file
#  3) creates the length file
#  4) parses the edirect .mat_peptide file, if necessary
#  5) parses the edirect .ftable file
#  6) parses the length file

my $orig_infasta_file = $infasta_file;
my $outfasta_file     = (opt_Get("-c", \%opt_HH)) ? $out_root . ".fg.fa" : undef;
if(defined $infasta_file) { 
  # note that we pass in a reference to %ref_seq_info_HA to wrapperGetInfoUsingEdirect()
  # and *not* a reference to %seq_info_HA. We will use %infasta_ref_seq_info_HA to 
  # store information on the reference sequence only.
  wrapperGetInfoUsingEdirect(undef, $ref_accn, $build_root, \%cds_tbl_HHA, \%mp_tbl_HHA, \%infasta_ref_seq_info_HA, \%ofile_info_HH,
                             \%opt_HH, $ofile_info_HH{"FH"}); 
  $nseq = process_input_fasta_file($infasta_file, $outfasta_file, \%seq_info_HA, \%opt_HH, $ofile_info_HH{"FH"}); 
  if(defined $outfasta_file) { # this will be true if -c
    $infasta_file = $outfasta_file; 
  }
}
else { # --infasta not used (default)
  wrapperGetInfoUsingEdirect($listfile, $ref_accn, $out_root, \%cds_tbl_HHA, \%mp_tbl_HHA, \%seq_info_HA, \%ofile_info_HH,
                             \%opt_HH, $ofile_info_HH{"FH"}); 
}

if($do_matpept) {  
  # validate the CDS:mat_peptide relationships that we read from the $matpept input file
  matpeptValidateCdsRelationships(\@cds2pmatpept_AA, \%{$cds_tbl_HHA{$ref_accn}}, \%{$mp_tbl_HHA{$ref_accn}}, opt_Get("-c", \%opt_HH), $seq_info_HA{"accn_len"}[0], $ofile_info_HH{"FH"});
}
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################
# Step 2. Fetch and process the reference genome sequence
##########################################################
my $step_desc = opt_Get("--skipfetch", \%opt_HH) ? "Processing the reference genome from previously fetched sequence file" : "Fetching all sequences and processing the reference genome";
$start_secs = outputProgressPrior($step_desc, $progress_w, $log_FH, *STDOUT);
my %mdl_info_HA = ();  # hash of arrays, values are arrays [0..$nmdl-1];
                       # see dnaorg.pm::validateModelInfoHashIsComplete() for list of all keys
                       # filled in wrapperFetchAndProcessReferenceSequence()
my %ftr_info_HA = ();  # hash of arrays, values are arrays [0..$nftr-1], 
                       # see dnaorg.pm::validateFeatureInfoHashIsComplete() for list of all keys
                       # filled in wrapperFetchAndProcessReferenceSequence()
my $sqfile = undef;    # pointer to the Bio::Easel::SqFile object we'll open in wrapperFetchAllSequencesAndProcessReferenceSequence()

# Call the wrapper function that does the following:
#   1) fetches the sequences listed in @{$seq_info_HAR->{"accn_name"} into a fasta file 
#      and indexes that fasta file, the reference sequence is $seq_info_HAR->{"accn_name"}[0].
#   2) determines information for each feature (strand, length, coordinates, product) in the reference sequence
#   3) determines type of each reference sequence feature ('cds-mp', 'cds-notmp', or 'mp')
#   4) fetches the reference sequence feature and populates information on the models and features
wrapperFetchAllSequencesAndProcessReferenceSequence(\%execs_H, \$sqfile, $out_root, $build_root, 
                                                    ($do_infasta) ? $infasta_ref_seq_info_HA{"accn_name"}[0] : undef,
                                                    ($do_infasta) ? $infasta_ref_seq_info_HA{"accn_len"}[0]  : undef,
                                                    ($do_infasta) ? $infasta_ref_seq_info_HA{"seq_len"}[0]   : undef,
                                                    ($do_infasta) ? $infasta_file                            : undef,
                                                    \%cds_tbl_HHA,
                                                    ($do_matpept) ? \%mp_tbl_HHA      : undef, 
                                                    ($do_matpept) ? \@cds2pmatpept_AA : undef, 
                                                    ($do_matpept) ? \@cds2amatpept_AA : undef, 
                                                    \%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, 
                                                    \%opt_HH, \%ofile_info_HH);

# verify our model, feature, and sequence info hashes are complete, 
# if validateFeatureInfoHashIsComplete() fails then the program will exit with an error message
my $nftr = validateFeatureInfoHashIsComplete  (\%ftr_info_HA, undef, $ofile_info_HH{"FH"}); # nftr: number of features
my $nmdl = validateModelInfoHashIsComplete    (\%mdl_info_HA, undef, $ofile_info_HH{"FH"}); # nmdl: number of homology models
if($nseq != validateSequenceInfoHashIsComplete(\%seq_info_HA, undef, \%opt_HH, $ofile_info_HH{"FH"})) { 
  if(defined $listfile) { 
    DNAORG_FAIL(sprintf("ERROR, number of stored sequences (%d) in seq_info_HA differs from number of accessions read from $listfile (%d)", validateSequenceInfoHashIsComplete(\%seq_info_HA, undef, \%opt_HH, $ofile_info_HH{"FH"}), $nseq), $ofile_info_HH{"FH"});
    # $seq_info_HA won't have any duplicate accessions because it was initially filled by parseListFile()
    # which dies if any duplicates are found. However, if somehow there were duplicates in %seq_info_HA,
    # validateSequenceInfoHashComplete() will die in error.
  }
  else { # --infasta enabled 
    DNAORG_FAIL(sprintf("ERROR, number of stored sequences (%d) in seq_info_HA differs from number of accessions read from $orig_infasta_file (%d)", validateSequenceInfoHashIsComplete(\%seq_info_HA, undef, \%opt_HH, $ofile_info_HH{"FH"}), $nseq), 1, $ofile_info_HH{"FH"});
  }
}    

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

##############################################################################
# Step 3. Verify we have the model files that we need to do homology searches
##############################################################################
my @model_file_A = (); # array of the model files we need
for(my $i = 0; $i < $nmdl; $i++) { 
  my $model_file = $build_root . ".$i.cm";
  if(! -s $model_file) { 
    DNAORG_FAIL("ERROR CM file $model_file should exist but it does not. Did you (successfully) run dnaorg_build.pl?", 1, $ofile_info_HH{"FH"});
  }
  for my $suffix ("i1m", "i1i", "i1f", "i1p") { 
    my $file = $model_file . "." . $suffix;
    if(! -s $file) { 
      DNAORG_FAIL("ERROR CM index file $file should exist but it does not. Did you (successfully) run dnaorg_build.pl?", 1, $ofile_info_HH{"FH"});
    }
  }
  # set mdl_info_HAR->{"cmfile"}[$i]
  $mdl_info_HA{"cmfile"}[$i] = $model_file;
}

# Validate the CMs we are about to use to annotate were actually created
# for the current reference sequence and annotation (same accession
# and features (CDS, etc.)).
# 
# One subtle case occurs when the accession name gets changed because
# a non-RefSeq sequence got promoted to become a RefSeq; in that
# RefSeq promotion case, the accession name will change and the script
# has to be rerun.
if(defined $infasta_file) { 
  $start_secs = outputProgressPrior("Skipping verification that CMs created for current reference $ref_accn (--infasta)", $progress_w, $log_FH, *STDOUT);
}
else { 
  $start_secs = outputProgressPrior("Verifying CMs were created for current reference $ref_accn", $progress_w, $log_FH, *STDOUT);
  validate_cms_built_from_reference(\%mdl_info_HA, \%opt_HH, \%ofile_info_HH);
}


outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###################################################################
# Step 4. (OPTIONAL) Search for origin sequences, if --origin used
###################################################################
# initialize error data structures
my %err_info_HA = (); 
initializeHardCodedErrorInfoHash(\%err_info_HA, $ofile_info_HH{"FH"});

my @err_ftr_instances_AHH = ();
my %err_seq_instances_HH = ();
error_instances_initialize_AHH(\@err_ftr_instances_AHH, \%err_seq_instances_HH, \%err_info_HA, \%ftr_info_HA, $ofile_info_HH{"FH"});

if(opt_IsUsed("--origin", \%opt_HH)) { 
  $start_secs = outputProgressPrior("Identifying origin sequences", $progress_w, $log_FH, *STDOUT);
  find_origin_sequences($sqfile, $origin_seq, \%seq_info_HA, \%err_seq_instances_HH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"}); 
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

#############
# TEMPORARY alternative origin detection method
#############
if(opt_IsUsed("--aorgmodel", \%opt_HH)) { 
  $start_secs = outputProgressPrior("Identifying origin sequences with profile HMM method", $progress_w, $log_FH, *STDOUT);
  aorg_find_origin_sequences($ofile_info_HH{"fullpath"}{"fasta"}, $sqfile, \%execs_H, $out_root, \%seq_info_HA, \%err_seq_instances_HH, \%err_info_HA, \%opt_HH, \%ofile_info_HH); 
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

####################################
# Step 5. Perform homology searches
####################################
my $seq_file = $ofile_info_HH{"fullpath"}{"fasta"};
validateFileExistsAndIsNonEmpty($seq_file, "dnaorg_annotate.pl:main", $ofile_info_HH{"FH"});

# determine how many sequence files we need to split $seq_file into to satisfy <n> sequences
# per job (where n is from --nseq <n>). If more than 1 and --local not used, we split up 
# the sequence file and submit jobs to farm, if only 1 or --local used, we just run it locally
my $nseqfiles = int($nseq / opt_Get("--nseq", \%opt_HH)); 
# int() takes the floor, so there can be a nonzero remainder. We don't add 1 though, 
# because split_fasta_file() will return the actual number of sequence files created
# and we'll use that as the number of jobs subsequently. $nfarmjobs is currently only
# the 'target' number of sequence files that we pass into split_fasta_file().
if($nseqfiles == 0) { $nseqfiles = 1; }

# concatenated tblout file, created by concatenating all of the @tmp_tblout_file_A files
my $tblout_file = $out_root . ".tblout";

# split up the CM file into individual CM files (if necessary)
my @tmp_tblout_file_A = (); # the array of tblout files, we'll remove after we're done, unless --keep
my @tmp_seq_file_A = ();    # the array of sequence files we'll remove after we're done, unless --keep (empty if we run locally)
my @tmp_err_file_A = ();    # the array of error files we'll remove after we're done, unless --keep (empty if we run locally)

if(! opt_Get("--skipscan", \%opt_HH)) { 
  if(($nseqfiles == 1) || (opt_Get("--local", \%opt_HH))) { 
    # run jobs locally
    $start_secs = outputProgressPrior("Running cmscan locally", $progress_w, $log_FH, *STDOUT);
    # run a different cmscan run for each file
    for(my $m = 0; $m < $nmdl; $m++) { 
      my $tmp_tblout_file = $out_root . ".m" . ($m+1) . ".tblout";
      my $tmp_stdout_file = opt_Get("-v", \%opt_HH) ? $out_root . ".m" . ($m+1) . ".cmscan" : "/dev/null";
      my $do_max = ($mdl_info_HA{"length"}[$m] <= (opt_Get("--smallthresh", \%opt_HH))) ? 1 : 0; # do not filter for very short models
      my $do_mid = 0;
      if(! $do_max) { 
        if($mdl_info_HA{"length"}[$m] <= (opt_Get("--midthresh", \%opt_HH))) { 
          $do_max = 1; # set filter thresholds to --mid for middle sized models
        } 
      }
      my $do_big = 0;
      if((! $do_max) && (! $do_mid)) { 
        if($mdl_info_HA{"length"}[$m] >= (opt_Get("--bigthresh",   \%opt_HH))) { 
          $do_big = 1; # use HMM mode for big models
        }
      }
      run_cmscan($execs_H{"cmscan"}, 1, $do_max, $do_mid, $do_big, $mdl_info_HA{"cmfile"}[$m], $seq_file, $tmp_stdout_file, $tmp_tblout_file, \%opt_HH, \%ofile_info_HH); # 1: run locally
      push(@tmp_tblout_file_A, $tmp_tblout_file);
    }
  }
  else { 
    # we need to split up the sequence file, and submit a separate set of cmscan jobs (one per model) for each
    my $nfasta_created = split_fasta_file($execs_H{"esl_ssplit"}, $seq_file, $nseqfiles, \%opt_HH, \%ofile_info_HH);
    # split_fasta_file will return the actual number of fasta files created, 
    # which can differ from the requested amount (which is $nseqfiles) that we pass in. It will put $nseq/$nseqfiles
    # in each file, which often means we have to make $nseqfiles+1 files.

    my $nfarmjobs = $nfasta_created * $nmdl;
    # now submit a job for each
    $start_secs = outputProgressPrior("Submitting $nfarmjobs cmscan jobs to the farm", $progress_w, $log_FH, *STDOUT);
    for(my $s = 1; $s <= $nfasta_created; $s++) { 
      my $tmp_seq_file = $seq_file . "." . $s;
      push(@tmp_seq_file_A,    $tmp_seq_file);
      for(my $m = 0; $m < $nmdl; $m++) { 
        my $tmp_tblout_file = $out_root . ".m" . ($m+1) . ".s" . $s . ".tblout";
        my $tmp_stdout_file = opt_Get("-v", \%opt_HH) ? $out_root . ".m" . ($m+1) . ".s" . $s . ".cmscan" : "/dev/null";
        my $do_max = ($mdl_info_HA{"length"}[$m] <= (opt_Get("--smallthresh", \%opt_HH))) ? 1 : 0; # do not filter for very short models
        my $do_mid = 0;
        if(! $do_max) { 
          if($mdl_info_HA{"length"}[$m] <= (opt_Get("--midthresh", \%opt_HH))) { 
            $do_max = 1; # set filter thresholds to --mid for middle sized models
          } 
        }
        my $do_big = 0;
        if((! $do_max) && (! $do_mid)) { 
          if($mdl_info_HA{"length"}[$m] >= (opt_Get("--bigthresh",   \%opt_HH))) { 
            $do_big = 1; # use HMM mode for big models
          }
        }
        run_cmscan($execs_H{"cmscan"}, 0, $do_max, $do_mid, $do_big, $mdl_info_HA{"cmfile"}[$m],  # 0: do not run locally
                   $tmp_seq_file, $tmp_stdout_file, $tmp_tblout_file, \%opt_HH, \%ofile_info_HH);
        push(@tmp_tblout_file_A, $tmp_tblout_file);
        push(@tmp_err_file_A,    $tmp_tblout_file . ".err"); # this will be the name of the error output file, set in run_cmscan
      }
    }
    outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
    
    # wait for the jobs to finish
    $start_secs = outputProgressPrior(sprintf("Waiting a maximum of %d minutes for all farm jobs to finish", opt_Get("--wait", \%opt_HH)), 
                                      $progress_w, $log_FH, *STDOUT);
    my $njobs_finished = waitForFarmJobsToFinish(\@tmp_tblout_file_A, \@tmp_err_file_A, "[ok]", opt_Get("--wait", \%opt_HH), opt_Get("--errcheck", \%opt_HH), \%{$ofile_info_HH{"FH"}});
    if($njobs_finished != $nfarmjobs) { 
      DNAORG_FAIL(sprintf("ERROR in main() only $njobs_finished of the $nfarmjobs are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", \%opt_HH)), 1, \%{$ofile_info_HH{"FH"}});
    }
    if(! opt_Get("--local", \%opt_HH)) { 
      outputString($log_FH, 1, "# "); # necessary because waitForFarmJobsToFinish() creates lines that summarize wait time and so we need a '#' before 'done' printed by outputProgressComplete()
    }
  } # end of 'else' entered if($nfarmjobs > 1 && ! --local)
    
  # concatenate all the tblout files into one 
  concatenateListOfFiles(\@tmp_tblout_file_A, $tblout_file, "dnaorg_annotate.pl", \%opt_HH, $ofile_info_HH{"FH"});
  
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
} # end of if(! opt_Get(--skipscan", \%opt_HH))

####################################################################
# Step 6. Parse homology search results into usable data structures
####################################################################
$start_secs = outputProgressPrior("Parsing cmscan results", $progress_w, $log_FH, *STDOUT);

my @mdl_results_AAH = ();  # 1st dim: array, 0..$nmdl-1, one per model
                           # 2nd dim: array, 0..$nseq-1, one per sequence
                           # 3rd dim: hash, keys are "p_start", "p_stop", "p_strand", "p_5overhang", "p_3overhang", "p_5seqflush", "p_3seqflush", "p_evalue", "p_fid2ref"
initialize_mdl_results(\@mdl_results_AAH, \%mdl_info_HA, \%seq_info_HA, \%opt_HH, $ofile_info_HH{"FH"});

# parse the cmscan results
parse_cmscan_tblout($tblout_file, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

# calculate the lengths of features
$start_secs = outputProgressPrior("Calculating predicted feature lengths ", $progress_w, $log_FH, *STDOUT);
results_calculate_predicted_lengths(\%mdl_info_HA, \%ftr_info_HA, $nseq, \@mdl_results_AAH, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###########################################################################################
# Step 7. Define names of all per-model and per-feature output files we're about to create
###########################################################################################
define_model_and_feature_output_file_names($out_root, \%mdl_info_HA, \%ftr_info_HA, $ofile_info_HH{"FH"});

################################################
# Step 8. Fetch predicted hits into fasta files
################################################
$start_secs = outputProgressPrior("Fetching cmscan predicted hits into fasta files", $progress_w, $log_FH, *STDOUT);
fetch_hits_given_results($sqfile, "predicted", \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################################
# Step 9. For features modeled by multiple models (e.g. multi-exon CDS)
#         combine all relevant hits into predicted feature sequences
#########################################################################
$start_secs = outputProgressPrior("Combining predicted exons into CDS", $progress_w, $log_FH, *STDOUT);
combine_model_hits("predicted", $seq_info_HA{"seq_name"}, \%mdl_info_HA, \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);
# we need to do this step even if there are no multi-exon CDS, because the function will update
# an important part of %ftr_info_HA indicating which file subsequent steps should use 
# (for models with 1 exon, this will be a file created prior to this (that is, we don't 
# wastefully create a new version of that file for single exon cases.))
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

##################################################################################
# Step 10. For features modeled by multiple other features (e.g. CDS comprised
#          of multiple mature peptides) combine the individual predicted feature
#          sequences into single sequences.
##################################################################################
$start_secs = outputProgressPrior("Combining predicted mature peptides into CDS", $progress_w, $log_FH, *STDOUT);
combine_feature_hits("predicted", $seq_info_HA{"seq_name"}, \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################
# Step 10b. Identify b5e, b5u, b3e, and b3u errors that occur when
#           alignment does not extend to the model boundary. 
#           We need to do this before step 12 which sets trc
#           errors because a sequence with a b5e error cannot
#           also have a trc error (we purposefully avoid checking
#           for early stops). 
#########################################################
$start_secs = outputProgressPrior("Identifying errors associated with incomplete alignment to the model", $progress_w, $log_FH, *STDOUT);
mdl_results_add_b5e_b5u_errors(\%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, 
                               \@err_ftr_instances_AHH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);


#########################################################
# Step 11. Examine each predicted feature to see if:
#          - start codon is valid
#          - where first in-frame stop exists (for all features with valid start codon)
#          Store the relevant information in @err_ftr_instances_AHH either to output later
#          or to examine more later.
#########################################################
$start_secs = outputProgressPrior("Identifying internal starts/stops in coding sequences", $progress_w, $log_FH, *STDOUT);
# Translate predicted CDS/mat_peptide sequences using esl-epn-translate to identify 
# in-frame stop codons.
wrapper_esl_epn_translate_startstop($esl_epn_translate, "predicted", \%ftr_info_HA, (@specstart_AA ? \@specstart_AA : undef), \%err_info_HA, \@err_ftr_instances_AHH, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

####################################################################
# Step 12. For all features that have a valid start but no in-frame
#          stop codon, look for a downstream in-frame stop
####################################################################
# TODO: make a function that performs this step
my %seq_name_index_H = (); # seqname_index_H{$seq_name} = <n>, means that $seq_name is the <n>th sequence name in the @{$seq_name_AR}} array
getIndexHashForArray($seq_info_HA{"seq_name"}, \%seq_name_index_H, $ofile_info_HH{"FH"});

my $ftr_idx;  # counter over features
for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
  foreach my $seq_name (keys %{$err_ftr_instances_AHH[$ftr_idx]{"ext"}}) { 
    my $seq_idx = $seq_name_index_H{$seq_name};
    my $seq_len       = $seq_info_HA{"seq_len"}[$seq_idx];   # differs from accn_len when the sequence is circular and dulicated
    my $accn_len      = $seq_info_HA{"accn_len"}[$seq_idx];  # original length of the sequence without any duplication

    # determine the model nearest to the end of the current feature, for which we have a prediction
    my $mdl_idx = undef;  
    if($ftr_info_HA{"annot_type"}[$ftr_idx] eq "model") { 
      $mdl_idx = $ftr_info_HA{"final_mdl"}[$ftr_idx];
      while(! exists $mdl_results_AAH[$mdl_idx][$seq_idx]{"p_start"}) { 
        $mdl_idx--;
        if($mdl_idx < $ftr_info_HA{"first_mdl"}[$ftr_idx]) { 
          DNAORG_FAIL(sprintf("ERROR ext error exists for $seq_name ftr: %s, but unable to find a model for this feature for which we have a prediction...", 
                              $ftr_info_HA{"out_tiny"}[$ftr_idx]), 1, $ofile_info_HH{"FH"});
        }
      }
    }
    elsif($ftr_info_HA{"annot_type"}[$ftr_idx] eq "multifeature") { 
      my @primary_children_idx_A = ();
      getPrimaryOrAllChildrenFromFeatureInfo(\%ftr_info_HA, $ftr_idx, "primary", \@primary_children_idx_A, $ofile_info_HH{"FH"});
      my $child_idx     = scalar(@primary_children_idx_A)-1;
      my $child_ftr_idx = $primary_children_idx_A[$child_idx];
      $mdl_idx = $ftr_info_HA{"final_mdl"}[$child_ftr_idx];
      while(! exists $mdl_results_AAH[$mdl_idx][$seq_idx]{"p_start"}) { 
        $mdl_idx--;
        if($mdl_idx < $ftr_info_HA{"first_mdl"}[$child_ftr_idx]) { 
          # out of models for this child, decrement child_idx and redefine $child_ftr_idx and $mdl_idx
          $child_idx--;
          if($child_idx < 0) { 
            DNAORG_FAIL(sprintf("ERROR ext error exists for $seq_name ftr: %s, but unable to find a model for this feature for which we have a prediction...", 
                                $ftr_info_HA{"out_tiny"}[$ftr_idx]), 1, $ofile_info_HH{"FH"});
          }
          $child_ftr_idx = $primary_children_idx_A[$child_idx];
          $mdl_idx = $ftr_info_HA{"final_mdl"}[$child_ftr_idx];
        }
      }
    }
    # printf("ext error for $seq_name ftr_idx: $ftr_idx %s mdl_idx: $mdl_idx\n", $ftr_info_HA{"out_tiny"}[$ftr_idx]);

    my $cur_start     = $mdl_results_AAH[$mdl_idx][$seq_idx]{"p_start"};
    my $cur_stop      = $mdl_results_AAH[$mdl_idx][$seq_idx]{"p_stop"};
    my $cur_strand    = $mdl_results_AAH[$mdl_idx][$seq_idx]{"p_strand"};
    my $plen          = (abs($cur_stop - $cur_start)) + 1;
    my $offset        = ($plen % 3); 
    my $posn_to_start = $cur_stop;
    if($cur_strand eq "+") { 
      $posn_to_start -= $offset; # only want to look for downstream stops in-frame with respect to the START codon, not the predicted STOP
      $posn_to_start++; # one past 
      # example: start=10 stop=20 length=11 offset=2 posn_to_start=(20-2+1)=19 because first possible downstream in-frame stop starts at posn 19 
      # (10..12,13..15,16..18,19..21 etc. are in-frame codon positions)
    }
    else { # strand is "-";
      $posn_to_start += $offset; # only want to look for downstream stops in-frame with respect to the START codon, not the predicted STOP
      $posn_to_start--; # one past 
      # example: start=20 stop=10 length=11 offset=2 posn_to_start=(10+2-1)=11 because first possible downstream in-frame stop starts at posn 11
      # (20..18,17..15,14..12,11..9 etc. are in-frame codon positions)
    }

    # check for a special case, where we are not circularized
    my $at_end_of_seq = 0; # set to 1 if we're at the end of the sequence
    if($cur_strand eq "+") { 
      if($posn_to_start > $accn_len) { 
        if(opt_Get("-c", \%opt_HH)) { 
          # genome is circularized, we've duplicated it, it's okay that our possible start is > $accn_len at this point, but we want to subtract $accn_len here
          $posn_to_start -= $accn_len;
        }
        elsif(($posn_to_start-1) == $accn_len) { 
          # a special case
          $at_end_of_seq = 1;
        }
        else {
          DNAORG_FAIL("ERROR when looking for inframe stop for sequence $seq_name, trying to start search at position $posn_to_start but length of sequence is $accn_len", 1, $ofile_info_HH{"FH"}); 
        }
      }
    }
    else { # cur_strand eq "-" 
      if($posn_to_start <= 0) { 
        if(opt_Get("-c", \%opt_HH)) { 
          # genome is circularized, we've duplicated it, it's okay that our possible start is < 0 at this point, but we want to add $accn_len here
          $posn_to_start += $accn_len;
        }
        elsif(($posn_to_start+1) == 1) { 
          # a special case
          $at_end_of_seq = 1;
        }
        else { 
          DNAORG_FAIL("ERROR when looking for inframe stop for sequence $seq_name, trying to start search at position $posn_to_start but length of sequence is $accn_len", 1, $ofile_info_HH{"FH"}); 
        }
      }
    }
    
    my $ext_corr_stop  = undef; # third position of the next in-frame stop beyond where the stop was expected
    my $ext_stop_codon = undef; # which stop codon is used as the next in-frame stop 

    if($at_end_of_seq) { 
      $ext_corr_stop = 0;
    }
    else { 
      ($ext_corr_stop, $ext_stop_codon) = check_for_downstream_stop($sqfile, $seq_name, $posn_to_start, $accn_len, $cur_strand, $ofile_info_HH{"FH"});
    }

    if($ext_corr_stop == 0) { 
      # no stop found: nst error

      my $updated_nst_errmsg = "";
      if($at_end_of_seq) { 
        my (undef, $out_stop) = create_output_start_and_stop($cur_start, $cur_stop, $accn_len, $seq_len, $ofile_info_HH{"FH"});
        $updated_nst_errmsg = "inferred stop codon position (3 nt 3' of $out_stop on $cur_strand strand) is off the end of the sequence";
      }
      error_instances_update(\@err_ftr_instances_AHH, undef, \%err_info_HA, $ftr_idx, "nst", $seq_name, $updated_nst_errmsg, $ofile_info_HH{"FH"});

      # remove ext error "maybe"
      error_instances_remove_maybe(\@err_ftr_instances_AHH, undef, \%err_info_HA, $ftr_idx, "ext", $seq_name, $ofile_info_HH{"FH"});
    }
    else { 
      # stop found: ext error
      $ext_corr_stop -= $offset; # account for offset
      $ext_corr_stop++; # ext_corr_stop is w.r.t to the next posn after the predicted stop (either +1 or -1 depending on strand), we want 
                        # it to be w.r.t the actual predicted stop, so we have to add one.
      error_instances_update(\@err_ftr_instances_AHH, undef, \%err_info_HA, $ftr_idx, "ext", $seq_name, $ext_corr_stop, $ofile_info_HH{"FH"});
      # remove nst error "maybe"
      error_instances_remove_maybe(\@err_ftr_instances_AHH, undef, \%err_info_HA, $ftr_idx, "nst", $seq_name, $ofile_info_HH{"FH"});
    }
  }
}

##############################################################
# Step 13. For all features with either a 'trc' or 'ext'
#          error, add a new stop position to the results_AAH
#          data structure the 'corrected' stop. We utilize
#          the 'cumlen' values in the @mdl_results_AAH to 
#          simplify this for multi-exon features. 
#############################################################

$start_secs = outputProgressPrior("Correcting homology search stop codon predictions to account for observed stop codons", $progress_w, $log_FH, *STDOUT);
results_calculate_corrected_stops(\%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@mdl_results_AAH, \@err_ftr_instances_AHH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#################################################
# Step 14. Identify overlap and adjacency errors
#################################################

$start_secs = outputProgressPrior("Identifying overlap and adjacency errors", $progress_w, $log_FH, *STDOUT);
results_calculate_overlaps_and_adjacencies(\%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@mdl_results_AAH, \@err_ftr_instances_AHH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

####################################################################
# Step 15. Finalize the mdl_results, and fill ftr_results and check
#          for incompatible error combinations.
####################################################################
my @ftr_results_AAH = (); # 1st dim: array, 0..$ftr_idx..$nftr-1, one per model
                          # only defined for $ftr_info_HA{"annot_type"}[$ftr_idx] == "multifeature"
                          # if $ftr_info_HA{"annot_type"}[$ftr_idx] == "model", then $mdl_results_AAH
                          # contains all the results we need.
                          # 2nd dim: array, 0..$nseq-1, one per sequence
                          # 3rd dim: hash, keys are "p_start", "p_stop", "p_strand", "p_5overhang", "p_3overhang", "p_5seqflush", "p_3seqflush", "p_evalue", "fid2ref"

initialize_ftr_results(\@ftr_results_AAH, \%ftr_info_HA, \%seq_info_HA, \%opt_HH, $ofile_info_HH{"FH"});

$start_secs = outputProgressPrior("Finalizing annotations and validating error combinations", $progress_w, $log_FH, *STDOUT);
# report str, nop, b3e, b3u errors, we need to know these before we call ftr_results_calculate()
mdl_results_add_str_nop_ost_b3e_b3u_errors($sqfile, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, 
                                           \@err_ftr_instances_AHH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"});

# calculate out_start, out_stop and out_stop_codon values, we need to know some of these before we call ftr_results_calculate()
mdl_results_calculate_out_starts_and_stops($sqfile, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, $ofile_info_HH{"FH"});

# set most of the multi-feature (e.g. cds-mp) errors
ftr_results_calculate($sqfile, \%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@ftr_results_AAH, \@mdl_results_AAH,
                      \%cds_tbl_HHA, \@err_ftr_instances_AHH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"});

# recalculate out_start and out_stop and out_stop_codon values, they may have changed for final mature peptides in each CDS
mdl_results_calculate_out_starts_and_stops($sqfile, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, $ofile_info_HH{"FH"});

# compare our final annotations to GenBank
mdl_results_compare_to_genbank_annotations(\%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@mdl_results_AAH, 
                                           \%cds_tbl_HHA, ($do_matpept) ? \%mp_tbl_HHA : undef, 
                                           \%opt_HH, $ofile_info_HH{"FH"});

# validate all of our error instances by checking for incompatibilities
# and enforcing required combinations, this function dies if any problems are found
error_instances_validate_all(\@err_ftr_instances_AHH, \%err_seq_instances_HH, \%err_info_HA, \%ftr_info_HA, \%seq_info_HA, \%opt_HH, $ofile_info_HH{"FH"});

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

# dump results, if debugging
# dump_results(*STDOUT, \@mdl_results_AAH, \%mdl_info_HA, \%seq_info_HA, \%opt_HH, \%ofile_info_HH);

# dumpArrayOfHashesOfHashes("Error instances (%err_ftr_instances_AHH)", \@err_ftr_instances_AHH, *STDOUT);

#################################################
# Step 16. Refetch corrected hits into new files
#################################################
$start_secs = outputProgressPrior("Fetching corrected matches into fasta files", $progress_w, $log_FH, *STDOUT);
fetch_hits_given_results($sqfile, "corrected", \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#############################################################################
# Step 17. For features modeled by multiple models (e.g. multi-exon CDS)
#          combine all relevant hits into corrected feature sequences.
#############################################################################
$start_secs = outputProgressPrior("Combining corrected exons into CDS", $progress_w, $log_FH, *STDOUT);
combine_model_hits("corrected", $seq_info_HA{"seq_name"}, \%mdl_info_HA, \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

#########################################################
# Step 18. For features modeled by multiple other features (e.g. CDS
#          comprised of multiple mature peptides), combine the
#          individual corrected feature sequences into single
#          sequences.
#########################################################
$start_secs = outputProgressPrior("Combining corrected mature peptides into CDS", $progress_w, $log_FH, *STDOUT);
combine_feature_hits("corrected", $seq_info_HA{"seq_name"}, \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###############################################################
# Step 19. Translate corrected feature sequences into proteins
###############################################################
if(! opt_Get("--skiptranslate", \%opt_HH)) { 
  $start_secs = outputProgressPrior("Translating corrected nucleotide features into protein sequences", $progress_w, $log_FH, *STDOUT);
  translate_feature_sequences("corrected", "corrected.translated", (@specstart_AA ? \@specstart_AA : undef), \%ftr_info_HA, \%ofile_info_HH);
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

################################################################################
# Step 20. (OPTIONAL) Create multiple alignments of DNA sequences, if --doalign
################################################################################
if(opt_Get("--doalign", \%opt_HH)) { 
  $step_desc = "Aligning and parsing corrected nucleotide hits";
  $start_secs = outputProgressPrior($step_desc, $progress_w, $log_FH, *STDOUT);
  align_hits(\%execs_H, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, \%ofile_info_HH); 
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

####################################################################################
# Step 21. (OPTIONAL) Create multiple alignments of protein sequences, if --doalign
####################################################################################
if(opt_Get("--doalign", \%opt_HH)) { 
  $start_secs = outputProgressPrior("Aligning translated protein sequences", $progress_w, $log_FH, *STDOUT);
  align_protein_sequences(\%execs_H, "corrected.translated", \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

#########################################
# Step 22. Output annotations and errors
#########################################
# open files for writing
openAndAddFileToOutputInfo(\%ofile_info_HH, "tbl",     $out_root . ".tbl",            1, "All annotations in tabular format");
openAndAddFileToOutputInfo(\%ofile_info_HH, "tblsum",  $out_root . ".tbl.summary",    1, "Summary of all annotations");
openAndAddFileToOutputInfo(\%ofile_info_HH, "failtbl", $out_root . ".fail.tbl",       1, "Annotations for all sequences with >= 1 failure in tabular format");
openAndAddFileToOutputInfo(\%ofile_info_HH, "errtbl",  $out_root . ".error.tbl",      1, "Annotations for all sequences with >= 1 error in tabular format");
openAndAddFileToOutputInfo(\%ofile_info_HH, "pererr",  $out_root . ".peraccn.errors", 1, "List of errors, one line per sequence");
openAndAddFileToOutputInfo(\%ofile_info_HH, "allerr",  $out_root . ".all.errors",     1, "List of errors, one line per error");
openAndAddFileToOutputInfo(\%ofile_info_HH, "errsum",  $out_root . ".errors.summary", 1, "Summary of all errors");

my @out_row_header_A  = (); # ref to array of output tokens for column or row headers
my @out_header_exp_A  = (); # same size of 1st dim of @out_col_header_AA and only dim of @out_row_header_A
                            # explanations of each header

###############
# error files #
###############
$start_secs = outputProgressPrior("Generating error code output", $progress_w, $log_FH, *STDOUT);

output_errors_header(\%ftr_info_HA, \%ofile_info_HH);
output_errors_all_sequences(\@err_ftr_instances_AHH, \%err_seq_instances_HH, \%ftr_info_HA, \%seq_info_HA, \%err_info_HA, \%opt_HH, \%ofile_info_HH);

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

############################
# tabular annotation files #
############################
$start_secs = outputProgressPrior("Generating tabular annotation output", $progress_w, $log_FH, *STDOUT);

# fill the data structures with the header information for the tabular annotation file
output_tbl_get_headings(\@out_row_header_A, \@out_header_exp_A, \%mdl_info_HA, \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);

# for each sequence, output the tabular annotation
my $nfail = output_tbl_all_sequences(\%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@mdl_results_AAH, \@ftr_results_AAH, \%opt_HH, \%ofile_info_HH);

# output the explanatory text
output_tbl_explanations(\@out_header_exp_A, \%ofile_info_HH);

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###########################################
# brief summary of annotations and errors #
###########################################
outputString($log_FH, 1, sprintf("#\n# Annotated %d accessions, %d (%6.4f fraction) had at least one annotation 'failure', see %s for all details.\n", 
                                 $nseq, $nfail, ($nfail/$nseq), $ofile_info_HH{"nodirpath"}{"tbl"}));

output_errors_summary($ofile_info_HH{"FH"}{"errsum"}, \@err_ftr_instances_AHH, \%err_seq_instances_HH, \%ftr_info_HA, \%seq_info_HA, \%err_info_HA, 1, \%opt_HH, \%ofile_info_HH); # 1: output to stdout

##################
# gap info files #
##################
if(opt_Get("--doalign", \%opt_HH)) { 
  output_gap_info($ofile_info_HH{"FH"}{"gap_perseq_all"},     $ofile_info_HH{"FH"}{"gap_pergap_all"},     0, 1, 0, 0, \%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, \%ofile_info_HH);
  output_gap_info($ofile_info_HH{"FH"}{"gap_perseq_not3"},    $ofile_info_HH{"FH"}{"gap_pergap_not3"},    0, 0, 1, 0, \%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, \%ofile_info_HH);
  output_gap_info($ofile_info_HH{"FH"}{"gap_perseq_special"}, $ofile_info_HH{"FH"}{"gap_pergap_special"}, 0, 0, 0, 1, \%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, \%ofile_info_HH);
}

################################
# output optional output files #
################################
if(exists $ofile_info_HH{"FH"}{"mdlinfo"}) { 
  dumpInfoHashOfArrays("Model information (%mdl_info_HA)", 0, \%mdl_info_HA, $ofile_info_HH{"FH"}{"mdlinfo"});
}
if(exists $ofile_info_HH{"FH"}{"ftrinfo"}) { 
  dumpInfoHashOfArrays("Feature information (%ftr_info_HA)", 0, \%ftr_info_HA, $ofile_info_HH{"FH"}{"ftrinfo"});
}
if(exists $ofile_info_HH{"FH"}{"seqinfo"}) { 
  dumpInfoHashOfArrays("Sequence information (%seq_info_HA)", 0, \%seq_info_HA, $ofile_info_HH{"FH"}{"seqinfo"});
}
if(exists $ofile_info_HH{"FH"}{"errinfo"}) { 
  dumpInfoHashOfArrays("Error information (%err_info_HA)", 0, \%err_info_HA, $ofile_info_HH{"FH"}{"errinfo"});
}

############
# Conclude #
############

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir_out, \%ofile_info_HH);

###############
# SUBROUTINES #
################################################################
# List of subroutines in this file, divided into categories. 
#
# The names of the subroutines are meant to be descriptive. The
# convention for subroutine names is to use underscores to separate
# words (e.g. 'get_value_from_array()') for subroutines in *this* file as
# opposed to a Perl module where they're named in camel caps
# (e.g. getValueFromArray()).
#
# Subroutines which I consider overly-complicated that are
# future targets for refactoring are labelled with "***":
#
#  Subroutines related to preparing CM files for searches:
#    concatenate_individual_cm_files()
#    validate_cms_built_from_reference()
#
#  Subroutines related to homology searches:
#    run_cmscan()
#    split_fasta_file()
#    split_cm_file()
#    parse_cmscan_tblout()
#
#  Subroutines related to creating fasta files of predicted hits/features: 
#    fetch_hits_given_results()
#    combine_model_hits()
#    combine_feature_hits()
#    combine_sequences()
#
#  Subroutines related to the esl-epn-translate.pl script:
#    parse_esl_epn_translate_startstop_outfile()
#    get_esl_epn_translate_altstart_opt()
#    wrapper_esl_epn_translate_startstop()
#
#  Subroutines related to determining and storing annotations/results:
#    initialize_mdl_results()
#    initialize_ftr_results()
#    results_calculate_predicted_lengths()
#    store_hit()
#    results_calculate_corrected_stops()
#    results_calculate_overlaps_and_adjacencies()
#    mdl_results_add_str_nop_ost_b3e_b3u_errors()
#    mdl_results_calculate_out_starts_and_stops()
#    mdl_results_compare_to_genbank_annotations()
#    ftr_results_calculate() ***
#    dump_results()
#
#  Subroutines related to origin sequences:
#    validate_origin_seq()
#    find_origin_sequences()
#    get_origin_output_for_sequence
#
#  Subroutines related to the error instance data structures:
#    error_instances_initialize_AHH()
#    error_instances_add()
#    error_instances_update()
#    error_instances_remove_maybe()
#    error_instances_remove_not_maybe()
#    error_instances_validate_all()
#
#  Subroutines related to output:
#    output_tbl_explanations()
#    output_tbl_get_headings()
#    output_tbl_get_headings_helper()
#    output_tbl_get_headings_explanation_helper()
#    output_tbl_all_sequences()
#    output_tbl_page_of_sequences()
#    output_errors_header()
#    output_errors_all_sequences()
#    output_errors_summary()
#    output_multifeature_relationships()
#    output_gap_info()
#
#  Miscellaneous subroutines that don't fit in the above categories:
#    find_inframe_stop()
#    combine_ajb_and_aja_strings()
#    compare_to_genbank_annotation
#    count_genbank_annotations
#    translate_feature_sequences
#    align_hits()
#    update_gap_array()
#    find_special_gap()
#    define_model_and_feature_output_file_names()
#    get_mdl_or_ftr_ofile_info_key()
#    align_protein_sequences
#    check_for_downstream_stop()
#    create_output_start_and_stop()
#    seq_name_from_msa_seq_name()
#    accn_name_from_seq_name()
#
#################################################################

#################################################################
#################################################################
#
#  Subroutines related to preparing CM files for searches:
#    validate_cms_built_from_reference()
#
#################################################################
# Subroutine : validate_cms_built_from_reference()
# Incept:      EPN, Mon Feb 29 11:21:11 2016
#
# Purpose:     Validate the CM files in an array were built from
#              the current reference, using information in 
#              $mdl_info_HAR->{"cksum"}.
#
# Arguments: 
#  $mdl_info_HAR:    REF to hash of arrays with information on the models, PRE-FILLED
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:     void
# 
# Dies: If at least one CM was not built from the current reference.
#
################################################################# 
sub validate_cms_built_from_reference { 
  my $sub_name = "validate_cms_built_from_reference()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($mdl_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # validate we have a complete model info hash
  my $nmdl = validateModelInfoHashIsComplete($mdl_info_HAR, undef, $FH_HR);

  my $mismatch_errmsg = ""; # we'll fill this with error messages about any checksum mismatches we find
  my $common_errmsg   = "This may mean that the GenBank annotation for the reference sequence changed since dnaorg_build.pl was run to create the CMs. Rerun dnaorg_build.pl.";

  for(my $i = 0; $i < $nmdl; $i++) { 
    my $model_file = $mdl_info_HAR->{"cmfile"}[$i];
    validateFileExistsAndIsNonEmpty($model_file, $sub_name, $FH_HR); 
    # get the checksum line from the CM file into a file
    my $cksum = `grep ^CKSUM $model_file | awk '{ print \$2 '}`;
    chomp $cksum;
    if($cksum =~ m/\r$/) { chop $cksum; } # remove ^M if it exists
    if($cksum != $mdl_info_HAR->{"checksum"}[$i]) { 
      $mismatch_errmsg .= sprintf("CM #%d checksum %d != alignment checksum: %d\n", $i+1, $cksum, $mdl_info_HAR->{"checksum"}[$i]);
    }
  }
  
  if($mismatch_errmsg ne "") { 
    DNAORG_FAIL("ERROR in $sub_name, checksum mismatch(es):\n$mismatch_errmsg\n$common_errmsg", 1, $FH_HR);
  }

  return;
}

#################################################################
#################################################################
#
#  Subroutines related to preparing CM files for searches:
#    run_cmscan()
#    split_fasta_file()
#    parse_cmscan_tblout()
#
#################################################################
# Subroutine : run_cmscan()
# Incept:      EPN, Mon Feb 29 15:09:22 2016
#
# Purpose:     Run Infernal's cmscan executable using $model_file
#              as the CM file on sequence file $seq_file.
#
# Arguments: 
#  $cmscan:          path to the cmscan executable file
#  $do_local:        '1' to run locally, '0' to submit job to farm
#  $do_max:          '1' to run with --max option, '0' to run with default options
#                    '1' is usually used only when $model_file contains a single 
#                    very short model
#  $do_mid:          '1' to run with --mid option, '0' to not
#  $do_big:          '1' to run with --mxsize 6144 option, '0' to run with default options
#  $model_file:      path to the CM file
#  $seq_file:        path to the sequence file
#  $stdout_file:     path to the stdout file to create, can be "/dev/null", or undef 
#  $tblout_file:     path to the cmscan --tblout file to create
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:     void
# 
# Dies: If at least one CM was not built from the current reference.
#       If $do_max and $do_big are both true.
################################################################# 
sub run_cmscan { 
  my $sub_name = "run_cmscan()";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmscan, $do_local, $do_max, $do_mid, $do_big, $model_file, $seq_file, $stdout_file, $tblout_file, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  if($do_big && $do_max) { 
    DNAORG_FAIL("ERROR in $sub_name, do_big and do_max are both true, only one should be.", 1, $FH_HR);
  }
  if($do_big && $do_mid) { 
    DNAORG_FAIL("ERROR in $sub_name, do_big and do_mid are both true, only one should be.", 1, $FH_HR);
  }
  if($do_max && $do_mid) { 
    DNAORG_FAIL("ERROR in $sub_name, do_max and do_mid are both true, only one should be.", 1, $FH_HR);
  }
  validateFileExistsAndIsNonEmpty($model_file, $sub_name, $FH_HR); 
  validateFileExistsAndIsNonEmpty($seq_file,   $sub_name, $FH_HR);

  my $opts = (opt_Get("-v", $opt_HHR)) ? " " : " --noali ";
  $opts .= " --cpu 0 --tblout $tblout_file --verbose ";
  if($do_max) { # no filtering
    $opts .= "--max -E 0.01 "; # with --max, a lot more FPs get through the filter, so we enforce an E-value cutoff
  }
  elsif($do_mid) { 
    $opts .= " --mid -E 0.1"; # with --mid, more FPs get through the filter, so we enforce an E-value cutoff
  }
  # finally add --nohmmonly if we're not a big model
  if(! $do_big) { # do not use hmm unless model is big
    $opts .= " --nohmmonly ";
  }

  my $cmd = "$cmscan $opts $model_file $seq_file > $stdout_file";

  # remove the tblout file if it exists, this is important because we'll use the existence and
  # final line of this file to determine when the jobs are finished, if it already exists, we'll
  # think the job is finished before it actual is.
  if(-e $tblout_file) { removeFileUsingSystemRm($tblout_file, $sub_name, $opt_HHR, $ofile_info_HHR); }

  # run cmscan, either locally or by submitting jobs to the farm
  if($do_local) { 
    # run locally
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);
  }
  else { 
    # submit job to farm and return
    my $jobname = "s" . removeDirPath($seq_file);
    my $errfile = $tblout_file . ".err";
    if(-e $errfile) { removeFileUsingSystemRm($errfile, $sub_name, $opt_HHR, $ofile_info_HHR); }
    my $farm_cmd = "qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n -l h_rt=288000,h_vmem=8G,mem_free=8G,reserve_mem=8G " . "\"" . $cmd . "\" > /dev/null\n";
    runCommand($farm_cmd, opt_Get("-v", $opt_HHR), $FH_HR);
  }

  return;
}


#################################################################
# Subroutine : split_fasta_file()
# Incept:      EPN, Tue Mar  1 09:30:10 2016
#
# Purpose: Split up a fasta file into <n> smaller files by calling
#          the esl-ssplit perl script.
#
# Arguments: 
#  $esl_ssplit:      path to the esl-ssplit.pl script to use
#  $fasta_file:      fasta file to split up
#  $nfiles:          desired number of files to split $fasta_file into
#  $opt_HHR:         REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
# 
# Returns:    Number of files actually created (can differ from requested
#             amount (which is $nfiles)).
#
# Dies:       if esl-ssplit command fails
#
################################################################# 
sub split_fasta_file { 
  my $sub_name = "split_fasta_file()";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($esl_ssplit, $fasta_file, $nfiles, $opt_HHR, $ofile_info_HHR) = @_;

  # we can only pass $FH_HR to DNAORG_FAIL if that hash already exists
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $outfile = $fasta_file . ".esl-ssplit";
  my $cmd = "$esl_ssplit -v -n $fasta_file $nfiles > $outfile";
  runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);

  # parse output to determine exactly how many files were created:
  # $esl_ssplit will have output exactly 1 line per fasta file it created
  my $nfiles_created = countLinesInFile($outfile, $FH_HR);

  if(! opt_Get("--keep", $opt_HHR)) { 
    runCommand("rm $outfile", opt_Get("-v", $opt_HHR), $FH_HR);
  }

  return $nfiles_created;
}


#################################################################
# Subroutine : parse_cmscan_tblout()
# Incept:      EPN, Tue Mar  1 13:56:46 2016
#
# Purpose:    Parse Infernal 1.1 cmscan --tblout output and store
#             results in $mdl_results_AAH.
#
# Arguments: 
#  $tblout_file:      tblout file to parse
#  $mdl_info_HAR:     REF to hash of arrays with information on the models, PRE-FILLED
#  $seq_info_HAR:     REF to hash of arrays with sequence information, PRE-FILLED
#  $mdl_results_AAHR: REF to results AAH, FILLED HERE
#  $FH_HR:            REF to hash of file handles
#
# Returns:    void
#
# Dies:       if we find a hit to a model or sequence that we don't
#             have stored in $mdl_info_HAR or $seq_name_AR
#
################################################################# 
sub parse_cmscan_tblout { 
  my $sub_name = "parse_cmscan_tblout()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($tblout_file, $mdl_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $FH_HR) = @_;
  
  # make an 'order hash' for the model names and sequence names,
  my %mdlname_index_H = (); # mdlname_index_H{$model_name} = <n>, means that $model_name is the <n>th model in the @{$mdl_info_HAR{*}} arrays
  my %seqname_index_H = (); # seqname_index_H{$seq_name} = <n>, means that $seq_name is the <n>th sequence name in @{$seq_info_HAR{*}} arrays
  getIndexHashForArray($mdl_info_HAR->{"cmname"},   \%mdlname_index_H, $FH_HR);
  getIndexHashForArray($seq_info_HAR->{"seq_name"}, \%seqname_index_H, $FH_HR);

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
      #Maize-streak_r23.NC_001346.ref.mft.4        -         NC_001346:genome-duplicated:NC_001346:1:2689:+:NC_001346:1:2689:+: -          cm        1      819     2527     1709      -    no    1 0.44   0.2  892.0         0 !   -
      my @elA = split(/\s+/, $line);
      my ($mdlname, $seqname, $mod, $mdlfrom, $mdlto, $from, $to, $strand, $score, $evalue) = 
          ($elA[0], $elA[2], $elA[4], $elA[5], $elA[6], $elA[7], $elA[8], $elA[9], $elA[14], $elA[15]);

      if(! exists $mdlname_index_H{$mdlname}) { 
        DNAORG_FAIL("ERROR in $sub_name, do not have information for model $mdlname read in $tblout_file on line $line_ctr", 1, $FH_HR);
      }
      if(! exists $seqname_index_H{$seqname}) { 
        DNAORG_FAIL("ERROR in $sub_name, do not have information for sequence $seqname read in $tblout_file on line $line_ctr", 1, $FH_HR);
      }

      my $mdlidx = $mdlname_index_H{$mdlname}; # model    index for the hit in results_AAH (1st dim of results_AAH)
      my $seqidx = $seqname_index_H{$seqname}; # sequence index for the hit in results_AAH (2nd dim of results_AAH)
      my $mdllen = $mdl_info_HAR->{"length"}[$mdlidx]; # model length, used to determine how far hit is from boundary of the model
      if((! exists $seq_info_HAR->{"accn_len"}[$seqidx]) || (! exists $seq_info_HAR->{"seq_len"}[$seqidx])) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, do not have length information for sequence $seqname, accession %s", $seq_info_HAR->{"accn_name"}[$seqidx]), 1, $FH_HR);
      }
      my $accn_len = $seq_info_HAR->{"accn_len"}[$seqidx]; # accession length, used to exclude storing of hits that start and stop after $seqlen, 
                                                           # which can occur in circular genomes, where we've duplicated the sequence
      my $seq_len  = $seq_info_HAR->{"seq_len"}[$seqidx]; # total length of sequence searched, could be $accn_len or two times $accn_len if -c used

      store_hit($mdl_results_AAHR, $mdlidx, $seqidx, $mdllen, $accn_len, $seq_len, $mdlfrom, $mdlto, $from, $to, $strand, $evalue, $FH_HR);
    }
  }
  close(IN);
  
  return;
}

#################################################################
#################################################################
#
#  Subroutines related to creating fasta files of predicted hits/features: 
#    fetch_hits_given_results()
#    combine_model_hits()
#    combine_feature_hits()
#    combine_sequences()
#
#################################################################
# Subroutine:  fetch_hits_given_results()
# Incept:      EPN, Wed Mar  2 15:25:55 2016
#
# Purpose:    Given the results data structure, fetch all
#             hits to fasta files.
#             As an extra step, fetch any 'appended' sequence
#             for models in which $mdl_info_HAR->{"append_num"}[$mdl_idx]
#             is non-zero. This will go to a separate file,
#             and we'll append it to sequence hits from other 
#             models in combine_feature_hits().
# Arguments: 
#  $sqfile:            REF to Bio::Easel::SqFile object, open sequence file containing sequences,
#                      usually $out_root . ".predicted", or $out_root . ".corrected"
#  $out_key:           key for the output files we'll create here, usually "predicted" or "corrected"
#  $mdl_info_HAR:      REF to hash of arrays with information on the models, PRE-FILLED
#  $seq_info_HAR:      REF to hash of arrays with sequence information, PRE-FILLED
#  $mdl_results_AAHR:  REF to results AAH, PRE-FILLED
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
# 
# Dies: If something is wrong with data in $mdl_info_HAR.
#################################################################
sub fetch_hits_given_results { 
  my $sub_name = "fetch_hits_given_results";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $out_key, $mdl_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $opt_HHR, $ofile_info_HHR) = @_;
    
  my $nmdl = scalar(@{$mdl_results_AAHR});
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $ofile_info_HHR->{"FH"});

  my $mdl_info_file_key        = $out_key . ".hits.fa";
  my $mdl_info_append_file_key = $out_key . ".hits.append.fa";

  for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    my @fetch_AA        = ();
    my @fetch_append_AA = ();
    my $mdl_name = $mdl_info_HAR->{"cmname"}[$mdl_idx];
    my $nseq2fetch        = 0; # number of sequences we'll fetch
    my $nseq2fetch_append = 0; # number of sequences we'll fetch for append file
    my $append_num = ($mdl_info_HAR->{"append_num"}[$mdl_idx]);

    for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      my $seq_name = $seq_info_HAR->{"seq_name"}[$seq_idx];
      if($seq_idx == 0 && (! opt_Get("--infasta", $opt_HHR)) && (! %{$mdl_results_AAHR->[$mdl_idx][$seq_idx]})) { 
        DNAORG_FAIL("ERROR in $sub_name(), no hit from model $mdl_name to the reference sequence $seq_name", 1, $ofile_info_HHR->{"FH"}); 
      }
      if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"}) { 
        # hit exists
        if(! exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_stop"}) { 
          DNAORG_FAIL("ERROR in $sub_name(), no stop value for hit of model $mdl_name to the sequence $seq_name", 1, $ofile_info_HHR->{"FH"}); 
        }
        my $start  = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"};
        # use corrected stop ("c_stop") if it exists
        my $stop    = (exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"c_stop"}) ? $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"c_stop"} : $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_stop"};
        my $strand  = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_strand"};

        # only enter following loop if "prv_trc_flag" is *NOT* set, 
        # if it is then the sequence should not be fetched due to 
        # an early stop (trc error) in a previous model for the same feature
        if(! $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"prv_trc_flag"}) { 
          my $new_name = $seq_name . "/" . $start . "-" . $stop;
          # we always put $start first and $stop second so that we can (usually) tell strand from the name,
          # "+" if start < stop, "-" if start > stop, ambiguous if start==stop
          push(@fetch_AA, [$new_name, $start, $stop, $seq_name]);
          $nseq2fetch++;

          # append sequence, if nec
          if($append_num > 0) { 
            my $append_start;
            my $append_stop;
            if($strand eq "+") { 
              $append_start = $stop + 1;
              $append_stop  = $stop + $append_num;
            }
            else { 
              $append_start = $stop - 1;
              $append_stop  = $stop - $append_num;
            }
            # only do the appending if the full region $append_start..$append_stop 
            # exists
            if(($append_start <= $seq_info_HAR->{"seq_len"}[$seq_idx]) && 
               ($append_stop  <= $seq_info_HAR->{"seq_len"}[$seq_idx])) { 
              my $append_new_name = $seq_name . "/" . $append_start . "-" . $append_stop;
              # we always put $start first and $stop second so that we can (usually) tell strand from the name,
              # "+" if start < stop, "-" if start > stop, ambiguous if start==stop

              # update mdl_results with 'append_start' and 'append_stop'
              $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"append_start"} = $append_start;
              $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"append_stop"}  = $append_stop;
              # printf("added $append_new_name $append_start $append_stop $seq_name to fetch_append_AA\n");
              push(@fetch_append_AA, [$append_new_name, $append_start, $append_stop, $seq_name]);
              
              $nseq2fetch_append++;
            }
            else { 
              if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"append_start"}) { 
                delete $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"append_start"};
              }
              if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"append_stop"}) { 
                delete $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"append_stop"};
              }
            }
          }
        }
      }
    } # end of for loop over $seq_idx

    my $fa_file               = $mdl_info_HAR->{$mdl_info_file_key}[$mdl_idx];
    my $fa_append_file        = $mdl_info_HAR->{$mdl_info_append_file_key}[$mdl_idx];
    my $ofile_info_key        = get_mdl_or_ftr_ofile_info_key("mdl", $mdl_idx, $mdl_info_file_key, $ofile_info_HHR->{"FH"});
    my $ofile_info_append_key = get_mdl_or_ftr_ofile_info_key("mdl", $mdl_idx, $mdl_info_append_file_key, $ofile_info_HHR->{"FH"});

    if($nseq2fetch > 0) { 
      $sqfile->fetch_subseqs(\@fetch_AA, undef, $fa_file);
      # save information on this to the output file info hash
      addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_key, $fa_file, 0, "fasta file with $out_key hits for model " . $mdl_info_HAR->{"out_tiny"}[$mdl_idx]);

      if($nseq2fetch_append > 0) { 
        $sqfile->fetch_subseqs(\@fetch_append_AA, undef, $fa_append_file);
        addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_append_key, $fa_append_file, 0, "fasta file with $out_key appended hits for model " . $mdl_info_HAR->{"out_tiny"}[$mdl_idx]);
      }
    }
    else { 
      # no sequences were fetched update the 
      $mdl_info_HAR->{$mdl_info_file_key}[$mdl_idx]        = "/dev/null"; # indicates to downstream functions that this file does not exist
      $mdl_info_HAR->{$mdl_info_append_file_key}[$mdl_idx] = "/dev/null"; # indicates to downstream functions that this file does not exist
    }
  } # end of for loop over model indices

  return;
}

#################################################################
# Subroutine:  combine_model_hits()
# Incept:      EPN, Thu Mar  3 11:39:17 2016
#
# Purpose:    For all features annotated by models 
#             ($ftr_info_HAR->{"annot_type"}[*] = "model")
#             assign (for single model features) or create
#             (for multiple model features) the hit fasta file
#             for each. Calls 'combine_sequences()' which does
#             much of the work.
#
# Arguments: 
#  $out_key:           key for the output files we'll create here, usually "predicted" or "corrected"
#  $seq_name_AR:       REF to array of sequence names, PRE-FILLED
#  $mdl_info_HAR:      REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:      REF to hash of arrays with information on the features, ADDED TO HERE
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If we have a problem reading the fasta files
#
################################################################# 
sub combine_model_hits { 
  my $sub_name = "combine_model_hits";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_key, $seq_name_AR, $mdl_info_HAR, $ftr_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nftr: number of features
  my $nmdl = validateModelInfoHashIsComplete  ($mdl_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nmdl: number of homology models

  my $mdl_info_file_key        = $out_key . ".hits.fa";
  my $mdl_info_append_file_key = $out_key . ".hits.append.fa";
  my $ftr_info_file_key        = $mdl_info_file_key;
  my $ftr_info_append_file_key = $mdl_info_append_file_key;

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { # we only do this for features annotated by models
      my @tmp_hit_fafile_A = (); # only relevant if this feature has multiple models
      my $mdl_idx        = $ftr_info_HAR->{"first_mdl"}[$ftr_idx];
      my $mdl_hit_fafile = $mdl_info_HAR->{$mdl_info_file_key}[$mdl_idx];
      my $at_least_one_fafile = 0; # set to '1' if at least one fa file exists (is not set to "/dev/null"
      if($mdl_hit_fafile ne "/dev/null") { 
        validateFileExistsAndIsNonEmpty($mdl_hit_fafile, $sub_name, $ofile_info_HHR->{"FH"});
        push(@tmp_hit_fafile_A, $mdl_hit_fafile);
        $at_least_one_fafile = 1;
      }
      #######################################
      # single model (e.g. exon) features
      #######################################
      if($ftr_info_HAR->{"nmodels"}[$ftr_idx] == 1) { 
        # a single model (e.g. exon) gene, we should already have the sequence from fetch_hits
        # this was initialized to something else, redefine it here:
        $ftr_info_HAR->{$ftr_info_file_key}[$ftr_idx] = $mdl_hit_fafile;
      }
      #######################################
      # multi model (e.g. exon) features
      #######################################
      else { 
        # more than one model's hit files need to be combined to make this feature 
        my $ftr_hit_fafile = $ftr_info_HAR->{$ftr_info_file_key}[$ftr_idx];
        for($mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx] + 1; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
          $mdl_hit_fafile = $mdl_info_HAR->{$mdl_info_file_key}[$mdl_idx];
          if($mdl_hit_fafile ne "/dev/null") { 
            validateFileExistsAndIsNonEmpty($mdl_hit_fafile, $sub_name, $ofile_info_HHR->{"FH"});
            push(@tmp_hit_fafile_A, $mdl_hit_fafile);
            $at_least_one_fafile = 1;
          }
        }
        if($at_least_one_fafile) { 
          # combine the sequences into 1 file
          combine_sequences(\@tmp_hit_fafile_A, $seq_name_AR, $ftr_hit_fafile, $opt_HHR, $ofile_info_HHR->{"FH"});
          
          my $ofile_info_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $ftr_info_file_key, $ofile_info_HHR->{"FH"});
          addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_key, $ftr_hit_fafile, 0, "fasta file with $out_key hits for feature " . $ftr_info_HAR->{"out_tiny"}[$ftr_idx] . " from " . $ftr_info_HAR->{"nmodels"}[$ftr_idx] . " combined model predictions");
        }
        else { 
          # no fasta files exist, redefine $ftr_info_HAR->{"$ftr_info_file_key"}[$ftr_idx] to 
          # /dev/null so downstream functions know that it should not exist
          $ftr_info_HAR->{$ftr_info_file_key}[$ftr_idx] = "/dev/null";
        }
      }

      # check if there's a file to append
      my $final_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$ftr_idx];
      my $mdl_ofile_info_append_key = get_mdl_or_ftr_ofile_info_key("mdl", $final_mdl_idx, $mdl_info_append_file_key, $ofile_info_HHR->{"FH"});
      if(exists $ofile_info_HHR->{"fullpath"}{$mdl_ofile_info_append_key}) { 
        # yes, there is
        my $mdl_hit_append_fafile = $ofile_info_HH{"fullpath"}{$mdl_ofile_info_append_key};
        if($mdl_hit_append_fafile ne "/dev/null") { 
          validateFileExistsAndIsNonEmpty($mdl_hit_append_fafile, $sub_name, $ofile_info_HHR->{"FH"});
        }
        # this was initialized to something else, redefine it here:
        $ftr_info_HAR->{$ftr_info_append_file_key}[$ftr_idx] = $mdl_hit_append_fafile;
      }
    } # end of 'if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model")'
  }
  
  return;
}

#################################################################
# Subroutine:  combine_feature_hits()
# Incept:      EPN, Thu Mar  3 12:15:55 2016
#
# Purpose:    For all features that are annotated by combining
#             multiple other features ($ftr_info_HAR->{"annot_type"}[*] 
#             = "multifeature", e.g. CDS made up of mature peptides)
#             create the feature fasta file for each. Calls 
#             'combine_sequences()' which does much of the work.
#
# Arguments: 
#  $out_key:           key for the output files we'll create here, usually "predicted" or "corrected"
#  $seq_name_AR:       REF to array of sequence names, PRE-FILLED
#  $ftr_info_HAR:      REF to hash of arrays with information on the features, ADDED TO HERE
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If we have a problem reading the fasta files
#
################################################################# 
sub combine_feature_hits { 
  my $sub_name = "combine_feature_hits";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }


  my ($out_key, $seq_name_AR, $ftr_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nftr: number of features

  my $ftr_info_file_key        = $out_key . ".hits.fa";
  my $mdl_info_append_file_key = $out_key . ".hits.append.fa";

  # printf("in $sub_name, out_key: $out_key, ftr_info_file_key: $ftr_info_file_key\n");

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") { # we only do this for features annotated by models
      # get the array of primary children feature indices for this feature
      my @primary_children_idx_A = (); # feature indices of the primary children of this feature
      getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "primary", \@primary_children_idx_A, $ofile_info_HHR->{"FH"});
      my @tmp_hit_fafile_A = ();
      my $at_least_one_fafile = 0; # set to 1 once we add a fasta file to @tmp_hit_fafile_A
      my $combined_ftr_hit_fafile = $ftr_info_HAR->{$ftr_info_file_key}[$ftr_idx];
      foreach my $cur_ftr_idx (@primary_children_idx_A) { 
        my $cur_ftr_hit_fafile = $ftr_info_HAR->{$ftr_info_file_key}[$cur_ftr_idx];
        if($cur_ftr_hit_fafile ne "/dev/null") { 
          validateFileExistsAndIsNonEmpty($cur_ftr_hit_fafile, $sub_name, $ofile_info_HHR->{"FH"});
          push(@tmp_hit_fafile_A, $cur_ftr_hit_fafile);
          $at_least_one_fafile = 1;
        }
        # check if this feature has a mandatory file to append
        my $final_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$cur_ftr_idx];
        my $mdl_ofile_info_append_key = get_mdl_or_ftr_ofile_info_key("mdl", $final_mdl_idx, $mdl_info_append_file_key, $ofile_info_HHR->{"FH"});
        if(exists $ofile_info_HHR->{"fullpath"}{$mdl_ofile_info_append_key}) {
          # it does, append it
          my $mdl_hit_append_fafile = $ofile_info_HH{"fullpath"}{$mdl_ofile_info_append_key};
          if($mdl_hit_append_fafile ne "/dev/null") { 
            validateFileExistsAndIsNonEmpty($mdl_hit_append_fafile, $sub_name, $ofile_info_HHR->{"FH"});
            push(@tmp_hit_fafile_A, $mdl_hit_append_fafile);
            $at_least_one_fafile = 1;
          } 
        }
      } # end of 'foreach $cur_ftr_idx'

      if($at_least_one_fafile) { 
        # combine the sequences into 1 file
        combine_sequences(\@tmp_hit_fafile_A, $seq_name_AR, $combined_ftr_hit_fafile, $opt_HHR, $ofile_info_HHR->{"FH"}); 

        my $ofile_info_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $ftr_info_file_key, $ofile_info_HHR->{"FH"});
        addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_key, $combined_ftr_hit_fafile, 0, "fasta file with $out_key hits for feature " . $ftr_info_HAR->{"out_tiny"}[$ftr_idx] . " from " . $ftr_info_HAR->{"nmodels"}[$ftr_idx] . " combined model predictions");
      } # end of 'if($at_least_one_fafile)'
      else { 
        # no fasta files exist, redefine $ftr_info_HAR->{"$ftr_info_file_key"}[$ftr_idx] to 
        # /dev/null so downstream functions know that it should not exist
        $ftr_info_HAR->{$ftr_info_file_key}[$ftr_idx] = "/dev/null";
      }
    }
  } # end of 'for' loop over $ftr_idx
  return;
}


#################################################################
# Subroutine:  combine_sequences()
# Incept:      EPN, Wed Mar  2 16:11:40 2016
#
# Purpose:    Helper function for combine_model_hits() and
#             combine_feature_hits().  Given an array of fasta files,
#             each with a different subsequence from the same parent
#             sequences, create a single new fasta file that has the
#             subsequences concatenated together.  An example is
#             stitching together exons into a CDS.  Uses BioEasel's
#             sqfile module.
#
# Arguments: 
#  $indi_file_AR: REF to array of fasta files to combine
#  $seq_name_AR:  REF to array with order of sequence names
#  $multi_file:   name of multi file to create
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:        REF to hash of file handles
#
# Returns:    void
#
# Dies:       If we have a problem reading the fasta files
#
#################################################################
sub combine_sequences {
  my $sub_name = "combine_sequences";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($indi_file_AR, $seq_name_AR, $multi_file, $opt_HHR, $FH_HR) = @_;
  
  my @sqfile_A                      = (); # array of open Bio::Easel::SqFile objects, one per indi_file_AR element
  my @sqfile_sqname_AA              = (); # 2D array of SqFile sequence names, read from the indi files, [0..$nfiles-1][0..$nseq_in_file_f]
  my $nfiles                        = scalar(@{$indi_file_AR}); # number of input fasta files
  my @seq_name2sqfile_sqname_map_AA = (); # 2D array, 1st dim: 0..$seq_idx..$nseq-1, 2nd dim 0..$file_idx..$nfiles, value: ssi index of $seq_idx subseq in file $f
  my @seq_name_exists_AA            = (); # 2D array, 1st dim: 0..$seq_idx..$nseq-1, 2nd dim: 0..$file_idx..$nfiles-1; value '1' if $seq_idx exists in file $f
  my @seq_name_coords_A             = (); # key: sequence name from @{$seq_name_AR}, value: concatenated set of coordinates of all subseqs in all $nfiles
  my @seq_name_fetch_me_A           = (); # [0..$i..$nseq_name-1]: values differ between first part of function and 2nd part of function:
                                          # in 'for($file_idx..' loop: value is either the most recent (highest) file idx for which seq $i existed 
                                          #                            or -2 if we've determined it is not in a contiguous set of files beginning at 1
                                          # after 'update values in @seq_name_fetch_AA' block: value is "0" if seq $i does not exist in a contiguous
                                          #                                                    block of sequences starting at 1, and "1" if it does
                                          #
                                          # we fetch each sequence that exists in a contiguous set of files in @sqfile_A starting with the first one
                                          # we don't fetch a sequence that doesn't exist in any of the sequences
                                          # we don't fetch a sequence that exists in at least one file in the middle of @sqfile_A

  my @sqname_AA                     = (); # 2D array, 1st dim [0..$file_idx..$nfiles-1], 2nd dim: 0 to number of sequences in file $file_idx

  my $seq_idx;  # counter over seq_name values
  my $file_idx; # counter for files

  # get index hash for @{$seq_name_AR}, this simplifies checking if
  # a given sequence name exists in @{$seq_name_AR}, and getting the idx
  # of that name
  my $nseq_name = scalar(@{$seq_name_AR});
  my %seq_name_idx_H = (); # key: $seq_name, value: idx of $seq_name in @{$seq_info_HAR->{"seq_name"}}
  getIndexHashForArray($seq_name_AR, \%seq_name_idx_H, $FH_HR);
  
  # initialize arrays
  for(my $seq_idx = 0; $seq_idx < $nseq_name; $seq_idx++) { 
    $seq_name_coords_A[$seq_idx]   = "";
    $seq_name_fetch_me_A[$seq_idx] = -1; # set to 0 or 1 in for($file_idx) loop below
    @{$seq_name2sqfile_sqname_map_AA[$seq_idx]} = ();
    for($file_idx = 0; $file_idx < $nfiles; $file_idx++) { 
      $seq_name2sqfile_sqname_map_AA[$seq_idx][$file_idx] = -1; # updated in block below if $seq_idx exists in $file_idx
    }
  }

  for($file_idx = 0; $file_idx < $nfiles; $file_idx++) { 
    if($indi_file_AR->[$file_idx] ne "/dev/null") { 
      # only enter the loop if the file exists (is not "/dev/null")
      # if the file does not exist it has zero sequences in it, so everything else after this
      # loop works as intended, we just have zero sequences added from this file

      validateFileExistsAndIsNonEmpty($seq_file, $indi_file_AR->[$file_idx], $FH_HR);
      
      # create the Bio::Easel object
      # first remove any old .ssi files that may exist
      my $ssi_file = $indi_file_AR->[$file_idx] . ".ssi";
      if(-e $ssi_file) { 
        removeFileUsingSystemRm($ssi_file, $sub_name, $opt_HHR, $FH_HR);
      }
      @{$sqname_AA[$file_idx]} = ();
      my $sqfile_nseq = 0;
      if(-s $indi_file_AR->[$file_idx]) { 
        $sqfile_A[$file_idx] = Bio::Easel::SqFile->new({ fileLocation => $indi_file_AR->[$file_idx] });
        $sqfile_nseq = $sqfile_A[$file_idx]->nseq_ssi;
      }
      else { 
        $sqfile_A[$file_idx] = undef;
        $sqfile_nseq = 0;
      }
      
      # get the names all of sequences in each file
      for(my $sqfile_seq_idx = 0; $sqfile_seq_idx < $sqfile_nseq; $sqfile_seq_idx++) { 
        $sqname_AA[$file_idx][$sqfile_seq_idx] = $sqfile_A[$file_idx]->fetch_seq_name_given_ssi_number($sqfile_seq_idx);
        
        # break down this name into the $seq_name and $coords
        my ($seq_name, $coords) = split("/", $sqname_AA[$file_idx][$sqfile_seq_idx]);
        if(! defined $coords || $coords !~ m/\-/) { 
          DNAORG_FAIL("ERROR in $sub_name, unable to parse sequence name $sqname_AA[$file_idx][$sqfile_seq_idx] into accession and coordinates", 1, $FH_HR);
        }
        if(! exists $seq_name_idx_H{$seq_name}) { 
          DNAORG_FAIL("ERROR in $sub_name, parsed sequence name $seq_name from $sqname_AA[$file_idx][$sqfile_seq_idx] does not exist in our seq_name_A array", 1, $FH_HR);
        }
        
        $seq_idx = $seq_name_idx_H{$seq_name};
        $seq_name2sqfile_sqname_map_AA[$seq_idx][$file_idx] = $sqfile_seq_idx;
        if($seq_name_coords_A[$seq_idx] ne "") { 
          $seq_name_coords_A[$seq_idx] .= ",";
        }
        $seq_name_coords_A[$seq_idx] .= $coords;
        $seq_name_exists_AA[$seq_idx][$file_idx] = 1;
        
        if($seq_name_fetch_me_A[$seq_idx] == ($file_idx-1)) { 
          $seq_name_fetch_me_A[$seq_idx]++; # this can make $seq_name_fetch_me_A[$seq_idx] rise to as high as $nfiles-1,
          # but it really only serves as a flag that this sequence exists in all files
          # starting at the first file ($file_idx == 0) up to the current file, else
          # we would have set this value to -2 in the iteration of the loop corresponding
          # to the file $file_idx in which it doesn't exist
        }
        else { 
          # if we get here, we went through at least one value for $file_idx 
          # in which this sequence did not exist
          $seq_name_fetch_me_A[$seq_idx] = -2; # we'll set this to 0 below
        }
      }
    } # end of 'if($indi_file_AR->[$file_idx] ne "/dev/null")'
  } # end of 'for' loop over file indexes

  # update values in @seq_name_fetch_me_A
  for($seq_idx = 0; $seq_idx < $nseq_name; $seq_idx++) { 
    # 3 possibilities: 
    #   $seq_name_fetch_me_A[$seq_idx] >=  0: if $seq_idx existed in >= 1 files, and did exist in a contiguous subset starting at file_idx 0
    #   $seq_name_fetch_me_A[$seq_idx] == -2: if $seq_idx existed in >= 1 files, but not in a contiguous subset of files starting at file_idx 0
    #   $seq_name_fetch_me_A[$seq_idx] == -1: if $seq_idx existed in    0 files
    $seq_name_fetch_me_A[$seq_idx] = ($seq_name_fetch_me_A[$seq_idx] >= 0) ? 1 : 0;
  } 

  # now for each seq_name that we want to fetch, fetch all subsequences for that 
  # sequence from all the individual files into a new sequence in a new file ($multi_file)
  open(OUT, ">", $multi_file) || die "ERROR unable to open $multi_file for writing";
  for($seq_idx = 0; $seq_idx < $nseq_name; $seq_idx++) { 
    # printf("HEYA in $sub_name creating $multi_file fetch me for sequence %s is %d\n", $seq_name_AR->[$seq_idx], $seq_name_fetch_me_A[$seq_idx]);
    if($seq_name_fetch_me_A[$seq_idx]) { 
      my $seq_name = $seq_name_AR->[$seq_idx];
      print OUT ">" . $seq_name . "/" . $seq_name_coords_A[$seq_idx] . "\n";
      for($file_idx = 0; $file_idx < $nfiles; $file_idx++) { 
        if($seq_name2sqfile_sqname_map_AA[$seq_idx][$file_idx] != -1) { 
          my $sqname = $sqname_AA[$file_idx][$seq_name2sqfile_sqname_map_AA[$seq_idx][$file_idx]];
          my $sqonly = $sqfile_A[$file_idx]->fetch_seq_to_fasta_string($sqname);
          $sqonly =~ s/^\>.+\n//;
          print OUT $sqonly;
        }
      }
    }
  }
  close(OUT);

  # clean up: remove all 'ssi' files we just created
  for($file_idx = 0; $file_idx < $nfiles; $file_idx++) { 
    if(-e $indi_file_AR->[$file_idx] . ".ssi") { 
      removeFileUsingSystemRm($indi_file_AR->[$file_idx] . ".ssi", $sub_name, $opt_HHR, $FH_HR);
    }
  }

  return;
}

#################################################################
#################################################################
#
#  Subroutines related to the esl-epn-translate.pl script:
#    parse_esl_epn_translate_startstop_outfile()
#    get_esl_epn_translate_altstart_opt
#    wrapper_esl_epn_translate_startstop()
#
#################################################################
# Subroutine:  parse_esl_epn_translate_startstop_outfile()
# Incept:      EPN, Fri Mar  4 13:56:56 2016
#
# Purpose:    Parse an output file from esl-epn-translate.pl run
#             with the --startstop option and store the relevant 
#             information we derive from it in @{$err_ftr_instances_AHHR}. 
#
#             Checks for and adds the following error codes for features
#             that with a non-matpept features (type ne "mp", e.g. "cds-mp", "cds-notmp"):
#
#             "str": if predicted start is invalid (detected by esl-epn-translate)
#
#             "stp": "maybe" added for features for which predicted stop is not a 
#                  valid in-frame stop (detected by esl-epn-translate), this is
#                  either removed, or updated later (because an stp error occurs only
#                  when the predicted stop position is not the end of a valid stop
#                  codon -- valid, out-of-frame stop is okay (not an stp error))
#
#             "ext": "maybe" added for features with no in-frame stop, removed or
#                  updated later (not in this function) when we look for an in-frame
#                  stop downstream of the predicted stop (ext and ntr are exclusive)
#
#             "nst": "maybe" added for features with no in-frame stop, removed or
#                  updated later (not in this function) when we look for an in-frame
#                  stop downstream of the predicted stop (ext and ntr are exclusive)
#
#             "trc": added here for features with early in-frame stop (detected by 
#                  esl-epn-translate) with temporary value, the -1 times the number 
#                  of nucleotides the in-frame stop is upstream of the predicted stop.
#                  This value is updated later (not in this function) to a more informative
#                  message for outputting.
#
# Arguments: 
#  $translate_outfile:      path to the file to parse
#  $ftr_idx:                feature index the translate output is for
#  $ftr_info_HAR:           REF to the feature info hash of arrays, PRE-FILLED
#  $err_info_HAR:           REF to the error info HA, PRE-FILLED
#  $err_ftr_instances_AHHR: REF to the error instances AAH, PRE-FILLED
#  $FH_HR:                  REF to hash of file handles
#
#
# Returns:    void
#
# Dies:       If we have trouble parsing $translate_outfile
#
#################################################################
sub parse_esl_epn_translate_startstop_outfile { 
  my $sub_name = "parse_esl_epn_translate_startstop_outfile";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($translate_outfile, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, $FH_HR) = @_;
  
  # is this a mature peptide?
  my $is_matpept = ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp") ? 1 : 0;
  my $append_num = 0; # the number of nucleotides that were appended to this features
                      # sequence past the final predicted nucleotide, we need to 
                      # account for these when we determine the correction
  if((defined $ftr_info_HAR->{"append_num"}) && (defined $ftr_info_HAR->{"append_num"}[$ftr_idx])) { 
    $append_num = $ftr_info_HAR->{"append_num"}[$ftr_idx];
  }

  open(IN, $translate_outfile) || fileOpenFailure($translate_outfile, "main", $!, "reading", $FH_HR);
  
  while(my $line = <IN>) { 
    # example line
    #HQ693465/1-306 1 1 304
    if($line =~ /^(\S+)\/(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/) { 
      my ($seq_name, $coords, $start_is_valid, $stop_is_valid, $first_stop_posn1) = ($1, $2, $3, $4, $5);

      # skip this sequence IFF we have a b5e error already for it, this means 
      # that the alignment does not extend to the 5' boundary of the model but
      # it does extend to the 5' boundary of the sequence (first seq posn of 
      # alignment is 1 (on + strand) or L (on - strand)). In this case we 
      # don't search for a trc error.
      if(! exists $err_ftr_instances_AHHR->[$ftr_idx]{"b5e"}{$seq_name}) { 
        # determine if we have an early stop
        my $cds_len            = dashCoordsStringCommaDelimitedToLength($coords, $sub_name, $FH_HR);
        my $final_codon_posn1  = $cds_len - 2; # final codon position 1 
        my $early_inframe_stop;
        my $corr_len           = 0; # number of nts to correct prediction by, changed if $early_inframe_stop
        if($first_stop_posn1 == $final_codon_posn1) { # first stop is the final codon
          $early_inframe_stop = 0;
        }
        elsif($first_stop_posn1 == 0) { # esl-epn-translate didn't find any in-frame stop
          $early_inframe_stop = 0;
        }
        else { # there is an early stop
          $early_inframe_stop = 1; 
          $corr_len = (-1 * ($final_codon_posn1 - $first_stop_posn1)) + $append_num; # negative because early stop shortens length
        }

        # We now have all of the relevant data on the current
        # CDS/mat_peptide sequence and we need to use to determine
        # what errors each sequence should throw (at least for those
        # that can tell should be thrown at this point in the
        # processing) as well as any corrections to stop predictions
        # that we should make prior to translation (trc errors are
        # thrown when this happens).
        # 
        # There are 4 relevant variables that dictate which errors
        # should be thrown/checked for later and whether a stop
        # correction should be made. The table gives all possible
        # combinations of values of those variables and lists the
        # outcome of each possibility.
        #
        # Variables we know from earlier processing:
        # $is_matpept:     '1' if current feature is a mature peptide, '0' if it is a CDS
        #
        # Variables derived from esl-epn-translate output we're currently parsing
        # $start_is_valid:     '1' if current feature's first 3 nt encode a valid start codon
        # $stop_is_valid:      '1' if current feature's final 3 nt encode a valid stop codon and total feature length is multiple of 3
        # $early_inframe_stop: '1' if an inframe stop exists prior to predicted stop
        #
        # 7 possibilities, each with different outcome (P1-P7):
        #                                                                                       | make correction to |
        # idx | is_matpept | start_is_valid | stop_is_valid | early_inframe_stop ||   errors    | stop coordinate?   |
        # ----|------------|----------------|---------------|--------------------||-------------|--------------------|
        #  P1 |      false |          false |          any  |                any ||         str |                 no |
        #  P2 |      false |           true |        false  |              false ||     stp ext?|     maybe (if ext) |
        #  P3 |      false |           true |        false  |               true ||     stp trc |                yes |
        #  P4 |      false |           true |         true  |              false ||        none |                 no |
        #  P5 |      false |           true |         true  |               true ||         trc |                yes |
        # -----------------------------------------------------------------------||-----------------------------------
        #  P6 |       true |            any |          any  |              false ||        ntr? |                 no |
        #  P7 |       true |            any |          any  |               true ||    trc ntr? |                yes |
        # ------------------------------------------------------------------------------------------------------------
        # 
        # in table above:
        # '?' after error code means that error is possible, we have to check for it later           
        # 'any' means that any value is possible, outcome is unaffected by value
        #
        if(! $is_matpept) { 
          if(! $start_is_valid) { # possibility 1 (P1)
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "str", $seq_name, "", $FH_HR);
            # printf("in $sub_name, feature index $ftr_idx, seq $seq_name $c possibility 1 (str)\n");
          }
          else { 
            # $start_is_valid is 1
            if(! $stop_is_valid) { 
              if(! $early_inframe_stop) { 
                # possibility 2 (P2): stp error, need to check for ext error later

                # add the 3 potential error codes, we'll check again later and possibly remove them
                # the 'stp' error is only a "maybe" because (! $stop_is_valid) implies it's not an
                # *IN-FRAME* valid stop codon, but we only throw 'stp' if the final 3 nt of the prediction
                # are not a valid stop codon, regardless of frame
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "stp", $seq_name, "maybe", $FH_HR);
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "ext", $seq_name, "maybe", $FH_HR);
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "nst", $seq_name, "maybe", $FH_HR);
                #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 2 (stp, maybe ext)\n");
              }
              else { # $early_inframe_stop is 1
                # possibility 3 (P3): stp and trc error

                # add the 2 potential error codes, we'll check again later and possibly remove them
                # the 'stp' error is only a "maybe" because of reason explained above similar case in possibility 2 above
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "stp", $seq_name, "maybe", $FH_HR);
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_name, $corr_len, $FH_HR);
                #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 3 (trc and stp)\n");
              }
            } # end of 'if(! $stop_is_valid)'
            else { # $stop_is_valid is 1
              if(! $early_inframe_stop) { 
                ; 
                # possibility 4 (P4): no errors, do nothing
                #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 4 (no errors)\n");
              }
              else { 
                # possibility 5 (P5): trc error
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_name, $corr_len, $FH_HR);
                #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 5 (trc)\n");
              }
            }              
          }
        } # end of 'if(! $is_matpept)'
        else { # $is_matpept is 1 
          if(! $early_inframe_stop) { 
            ; 
            # possibility 6 (P6): maybe ntr error later, but can't check for it now, do nothing;
            #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 6 (no error)\n");
          }
          else { # $early_inframe_stop is '1'
            # possibility 7 (P7): trc error, maybe ntr error later, but can't check for it now
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_name, $corr_len, $FH_HR);
          }
        }
      } # end of 'if(! exists b5e error)'
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse esl-epn-translate.pl output line:\n$line\n", 1, $FH_HR);
    }
  }
  close(IN);
  
  return;
}

#################################################################
# Subroutine:  get_esl_epn_translate_altstart_opt
# Incept:      EPN, Thu Mar 10 11:23:48 2016
#
# Purpose:    Get the command line option string for esl-epn-translate
#             for the -altstart option given the %{$ftr_info_HAR}
#             the feature index $ftr_idx and the @{$specstart_AAR}.
#             $specstart_AAR can (and often will be) undef, in which
#             case we just return "", indicating the option isn't
#             going to be used. Also if the type of the feature is not
#             a CDS, then we will also return "".
#
# Arguments: 
#  $ftr_info_HAR:      REF to hash of arrays with information on the features, ADDED TO HERE
#  $ftr_idx:           index in %{$ftr_info_HAR} we are interested in
#  $specstart_AAR:     REF to the 2D array with specified start codons
#
# Returns:    Return option string for esl-epn-translate -altstart,
#             e.g. "-altstart GCC". If none, return "".
#
# Dies:       never
# 
################################################################# 
sub get_esl_epn_translate_altstart_opt { 
  my $sub_name = "get_esl_epn_translate_altstart_opt()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $ftr_idx, $specstart_AAR) = @_;

  if(! defined $specstart_AAR) { 
    return ""; 
  }
  elsif(($ftr_info_HA{"type"}[$ftr_idx] ne "cds-mp") && ($ftr_info_HA{"type"}[$ftr_idx] ne "cds-notmp")) { 
    return "";
  }
  else { 
    my $codon_str = "";
    # determine the CDS index of this feature
    my $cds_idx = 0;
    for(my $cur_ftr_idx = 0; $cur_ftr_idx <= $ftr_idx; $cur_ftr_idx++) { 
      if(($ftr_info_HA{"type"}[$cur_ftr_idx] eq "cds-mp") || ($ftr_info_HA{"type"}[$cur_ftr_idx] eq "cds-notmp")) { 
        $cds_idx++; 
      }
    }
    foreach my $codon (@{$specstart_AAR->[($cds_idx-1)]}) { 
      if($codon_str ne "") { 
        $codon_str .= ",";
      }
      $codon_str .= $codon;
    }
    if($codon_str ne "") { 
      return "-altstart $codon_str";
    }
    else { 
      return "";
    }
  }
}

#################################################################
# Subroutine:  wrapper_esl_epn_translate_startstop()
# Incept:      EPN, Thu Mar 10 15:04:53 2016
#
# Purpose:    For each feature's predicted hits in a fasta file,
#             call 'esl-epn-translate' using the -startstop option 
#             to investigate where the in frame stop codons are.
#
# Arguments: 
#  $esl_epn_translate:  path to esl-epn-translate.pl executable
#  $out_key:            key for output files in %{$ftr_info_HAR}
#  $ftr_info_HAR:       REF to hash of arrays with information on the features, ADDED TO HERE
#  $specstart_AAR:      REF to the 2D array with specified start codons, can be undef
#  $err_info_HAR:       REF to the error info hash of arrays, PRE-FILLED
#  $err_ftr_instances_AHHR: REF to error instances AHH, PRE-FILLED with at least trc and ext errors
#  $opt_HHR:            REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:     REF to the 2D hash of output file information
#
# Returns: void
#
################################################################# 
sub wrapper_esl_epn_translate_startstop { 
  my $sub_name = "wrapper_esl_epn_translate_startstop()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($esl_epn_translate, $out_key, $ftr_info_HAR, $specstart_AAR, $err_info_HAR, $err_ftr_instances_AHHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nftr: number of features

  my $ftr_info_fa_file_key  = $out_key . ".hits.fa";
  my $ftr_info_out_file_key = $out_key . ".hits.fa.esl-epn-translate";

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $ftr_hit_fafile            = $ftr_info_HAR->{$ftr_info_fa_file_key}[$ftr_idx];
    my $esl_epn_translate_outfile = $ftr_info_HAR->{$ftr_info_out_file_key}[$ftr_idx];
  
    if($ftr_hit_fafile ne "/dev/null") { # if this is set to /dev/null we know it's not supposed to exist, so we skip this feature

      # deal with alternative starts here
      my $altstart_opt = get_esl_epn_translate_altstart_opt($ftr_info_HAR, $ftr_idx, $specstart_AAR);
      
      if(! opt_Get("--skiptranslate", $opt_HHR)) { 
        # use esl-epn-translate.pl to examine the start and stop codons in each feature sequence
        $cmd = $esl_epn_translate . " $altstart_opt -startstop $ftr_hit_fafile > $esl_epn_translate_outfile";
        runCommand($cmd, opt_Get("-v", $opt_HHR), $ofile_info_HHR->{"FH"});
      }
      else { # --skiptranslate, validate the output file exists
        validateFileExistsAndIsNonEmpty($esl_epn_translate_outfile, $sub_name, $ofile_info_HHR->{"FH"});
      }
      
      # parse the output
      parse_esl_epn_translate_startstop_outfile($esl_epn_translate_outfile, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, $ofile_info_HHR->{"FH"});
      if((! opt_Get("--keep", $opt_HHR)) && (! opt_Get("--skiptranslate", $opt_HHR))) { 
        removeFileUsingSystemRm($esl_epn_translate_outfile, $sub_name, $opt_HHR, $ofile_info_HHR);
      }
      elsif(! opt_Get("--skiptranslate", $opt_HHR)) { 
        my $ofile_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $ftr_info_out_file_key, $ofile_info_HHR->{"FH"});
        addClosedFileToOutputInfo($ofile_info_HHR, $ofile_key, $esl_epn_translate_outfile, 0, sprintf("esl-epn-translate.pl output file for feature %s", $ftr_info_HA{"out_tiny"}[$ftr_idx]));
      }
    }
  }
  return;
}

#################################################################
#################################################################
#  Subroutines related to determining and storing annotations/results:
#    initialize_mdl_results()
#    initialize_ftr_results()
#    results_calculate_predicted_lengths()
#    store_hit()
#    results_calculate_corrected_stops()
#    results_calculate_overlaps_and_adjacencies()
#    mdl_results_add_str_nop_ost_b3e_b3u_errors()
#    mdl_results_calculate_out_starts_and_stops()
#    mdl_results_compare_to_genbank_annotations()
#    ftr_results_calculate
#    dump_results()
#
#################################################################
# Subroutine: initialize_mdl_results()
# Incept:     EPN, Tue Mar 15 05:30:18 2016
#
# Purpose:    Initialize the mdl_results_AAH data structure.
#
# Args:
#  $mdl_results_AAHR: REF to the model results data structure
#                     INITIALIZED HERE
#  $mdl_info_HAR:     REF to hash of arrays with information 
#                     on the models, PRE-FILLED
#  $seq_info_HAR:     REF to hash of arrays with information 
#                     on the sequences, PRE-FILLED
#  $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:            REF to hash of file handles
#
# Returns: void
#
# Dies: If mdl_info_HAR or seq_info_HAR is invalid or incomplete.
#
#################################################################
sub initialize_mdl_results { 
  my $sub_name = "initialize_mdl_results()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_results_AAHR, $mdl_info_HAR, $seq_info_HAR, $opt_HHR, $FH_HR) = @_;

  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR);           # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences

  @{$mdl_results_AAHR} = ();
  for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    @{$mdl_results_AAHR->[$mdl_idx]} = ();
    for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      %{$mdl_results_AAHR->[$mdl_idx][$seq_idx]} = ();
    }
  }

  return;
}

#################################################################
# Subroutine: initialize_ftr_results()
# Incept:     EPN, Tue Mar 15 06:01:32 2016
#
# Purpose:    Initialize the ftr_results_AAH data structure.
#
# Args:
#  $ftr_results_AAHR: REF to the feature results data structure
#                     INITIALIZED HERE
#  $ftr_info_HAR:     REF to hash of arrays with information 
#                     on the features, PRE-FILLED
#  $seq_info_HAR:     REF to hash of arrays with information 
#                     on the sequences, PRE-FILLED
#  $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:            REF to hash of file handles
#
# Returns: void
#
# Dies: If ftr_info_HAR or seq_info_HAR is invalid or incomplete.
#
#################################################################
sub initialize_ftr_results { 
  my $sub_name = "initialize_ftr_results()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_results_AAHR, $ftr_info_HAR, $seq_info_HAR, $opt_HHR, $FH_HR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR);           # nftr: number of features
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences

  @{$ftr_results_AAHR} = ();
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { 
      # for any feature $ftr_idx with 'annot_type' eq 'model', 
      # we will have all the results we need in @mdl_results_AAH,
      # and $ftr_results_AAHR->[$ftr_idx] will remain undef
      $ftr_results_AAHR->[$ftr_idx] = undef;
    }
    else { 
      @{$ftr_results_AAHR->[$ftr_idx]} = ();
      for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
        %{$ftr_results_AAHR->[$ftr_idx][$seq_idx]} = ();
      }
    }
  }
  return;
}

#################################################################
# Subroutine:  results_calculate_predicted_lengths()
# Incept:      EPN, Tue Mar  8 16:12:44 2016
#
# Purpose:    For all results in @{$mdl_results_AAH}, determine the
#             length ("len") and cumulative length ("cumlen") for 
#             each prediction. For model M of feature F, the cumulative
#             length" is the length in nucleotides in the
#             predictions of all other models M' that also model
#             feature F, and are 5' of M, plus the length of M.
#             For example, if M is the third exon of a feature, and has
#             length 201, and the lengths of the first two exons are 99 
#             and 51 nucleotides respectively, the cumulative length of
#             M is 351.
# 
#             This function should be called before results_calculate_corrected_stops(),
#             which will update the length values as necessary for the corrected stops.
#             If, here, any "c_stop" values exist (indicating that results_calculate_corrected_stops()
#             was already called, then we die in error.
#
# Arguments: 
#  $mdl_info_HAR:     REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:     REF to hash of arrays with information on the features, ADDED TO HERE
#  $nseq:             number of sequences we have results for
#  $mdl_results_AAHR: REF to results AAH, ADDED TO HERE
#  $FH_HR:            REF to hash of file handles
#
# Returns:    void
# Dies:       if any "c_stop" results values already exist,
#             indicating that results_calculate_corrected_stops()
#             was already called (which it should *not* have been)
################################################################# 
sub results_calculate_predicted_lengths {
  my $sub_name = "results_calculate_predicted_lengths()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_info_HAR, $ftr_info_HAR, $nseq, $mdl_results_AAHR, $opt_HHR, $FH_HR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nmdl = validateModelInfoHashIsComplete  ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models

  my $start_key   = "p_start";
  my $stop_key    = "p_stop";
  my $strand_key  = "p_strand";
  my $len_key     = "len";
  my $cumlen_key  = "cumlen";
  my $forbid_key  = "c_stop";       # this key should not be set, if it already is (indicating results_calculate_corrected_stops()
                                    # was already called, then we die
  
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { 
      for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
        my $cumlen = 0;
        for(my $mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx]; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
          if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$start_key}) { 
            if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$forbid_key}) { 
              DNAORG_FAIL("ERROR in $sub_name(), results_AAHR->[$mdl_idx][$seq_idx]{$forbid_key} exists, but it shouldn't.", 1, $FH_HR);
            }
            my $len  = abs($mdl_results_AAHR->[$mdl_idx][$seq_idx]{$stop_key} - $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$start_key}) + 1;
            $cumlen += $len;
            $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$len_key}     = $len;
            $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key}  = $cumlen;
          }        
        }
      }
    }
  }
  
  return;
}

#################################################################
# Subroutine : store_hit()
# Incept:      EPN, Tue Mar  1 14:33:38 2016
#
# Purpose:    Helper function for parse_cmscan_tblout().
#             Given info on a hit and a ref to the results AAH, 
#             store info on it. 
#
# Arguments: 
#  $mdl_results_AAHR:  REF to results AAH, FILLED HERE
#  $mdlidx:            model index, 1st dim index in results_AAH to store in
#  $seqidx:            sequence index, 2nd dim index in results_AAH to store in
#  $mdllen:            model length
#  $accn_len:          sequence length (before duplicating, if relevant)
#  $seq_len:           length of sequence actually searched (after duplicating, if relevant)
#  $mdlfrom:           start position of hit
#  $mdlto:             stop position of hit
#  $seqfrom:           start position of hit
#  $seqto:             stop position of hit
#  $strand:            strand of hit
#  $evalue:            E-value of hit
#  $FH_HR:             REF to hash of file handles
#
# Returns:    void
#
# Dies:       never
#
################################################################# 
sub store_hit { 
  my $sub_name = "store_hit()";
  my $nargs_exp = 13;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_results_AAHR, $mdlidx, $seqidx, $mdllen, $accn_len, $seq_len, $mdlfrom, $mdlto, $seqfrom, $seqto, $strand, $evalue, $FH_HR) = @_;

  # only consider hits where either the start or end are less than the total length
  # of the genome. Since we sometimes duplicate all genomes, this gives a simple 
  # rule for deciding which of duplicate hits we'll store 
  if(($seqfrom <= $accn_len) || ($seqto <= $accn_len)) { 

    my $p_5seqflush = undef; # set to '1' if hit is flush with 5' end of sequence (1st nt in seq is 1st nt in hit), else 0
    my $p_3seqflush = undef; # set to '1' if hit is flush with 3' end of sequence (final nt in seq is final nt in hit), else 0

# Code block below is how we used to modify start/stop if hit spanned
# stop..start boundary, now we just store it as it is, and then modify
# the start and/or stop when we output its coordinates. This removes a
# nasty off-by-one with negative indices that kept recurring in
# various situations.
# 
# The commented out code is kept here for reference
#
#    # deal with case where one but not both of from to is > L:
#    if($seqfrom > $L || $seqto > $L) { 
#      $seqfrom -= $L; 
#      $seqto   -= $L; 
#      if($seqfrom < 0)  { $seqfrom--; }
#      if($seqto   < 0)  { $seqto--; }
#    }
    
    # we only store the first hit we see, this is safe because we know 
    # that this will be the lowest E-value
    if(%{$mdl_results_AAHR->[$mdlidx][$seqidx]}) { 
      # a hit for this model:seq pair already exists, make sure it has a lower E-value than the current one
      if($mdl_results_AAHR->[$mdlidx][$seqidx]{"p_evalue"} > $evalue) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, already have hit stored for model index $mdlidx seq index $seqidx with higher evalue (%g > %g), this implies hits were not sorted by E-value...", $mdl_results_AAHR->[$mdlidx][$seqidx]{"evalue"}, $evalue), 1, $FH_HR); 
      }
    }
    else { 
      # no hit yet exists, make one
      # determine if hit extends to very end of sequence on 5' and 3' end
      if($strand eq "+") { 
        $p_5seqflush = (($seqfrom == 1)         || ($seqfrom == ($accn_len+1))) ? 1 : 0; # $accn_len+1 only possible if seq is duplicated
        $p_3seqflush = (($seqto   == $accn_len) || ($seqto   == $seq_len))      ? 1 : 0; # $accn_len == $seq_len if sequence is not duplicated
      }
      else { # strand is -
        $p_5seqflush = (($seqfrom == $accn_len) || ($seqfrom == $seq_len))      ? 1 : 0; # $accn_len == $seq_len if sequence is not duplicated
        $p_3seqflush = (($seqto   == 1)         || ($seqto   == ($accn_len+1))) ? 1 : 0; # $accn_len+1 only possible if seq is duplicated
      }
      %{$mdl_results_AAHR->[$mdlidx][$seqidx]} = ();
      $mdl_results_AAHR->[$mdlidx][$seqidx]{"p_start"}     = $seqfrom;
      $mdl_results_AAHR->[$mdlidx][$seqidx]{"p_stop"}      = $seqto;
      $mdl_results_AAHR->[$mdlidx][$seqidx]{"p_strand"}    = $strand;
      $mdl_results_AAHR->[$mdlidx][$seqidx]{"p_5overhang"} = ($mdlfrom - 1);
      $mdl_results_AAHR->[$mdlidx][$seqidx]{"p_3overhang"} = ($mdllen - $mdlto);
      $mdl_results_AAHR->[$mdlidx][$seqidx]{"p_5seqflush"} = $p_5seqflush;
      $mdl_results_AAHR->[$mdlidx][$seqidx]{"p_3seqflush"} = $p_3seqflush;
      $mdl_results_AAHR->[$mdlidx][$seqidx]{"p_evalue"}    = $evalue;
    }
  }

  return;
}

#################################################################
# Subroutine:  results_calculate_corrected_stops()
# Incept:      EPN, Tue Mar  8 20:23:29 2016
#
# Purpose:    For all results in @{$mdl_results_AAH} for which an 
#             'trc' or 'ext' error exists in \@err_ftr_instances_AHH, 
#             determine what the corrected stop position should be
#             and set it as "c_stop".
#
#             Checks for and updates the following error codes for features
#             for which "annot_type" eq "model":
#
#             "stp": either updates or removes existing error instances of stp
#                  added in parse_esl_epn_translate_startstop_outfile().
#                  Checks to see if predicted stop is a valid stop codon,
#                  if it is, removes the error instance, else it updates it
#                  by making the associated error message more informative.
#
#                  valid in-frame stop (detected by esl-epn-translate), this is
#                  either removed, or updated later (because an stp error is only
#                  when the predicted stop position is not the end of a valid stop
#                  codon -- valid, out-of-frame stop is okay (not an stp error))
#
#             "trc": updates all 'trc' errors originally added in 
#                  parse_esl_epn_translate_startstop_outfile() by making their 
#                  error messages more informative to include information
#                  on the corrected stop position which is calculated in this
#                  function.
#
#             "ext": updates all 'ext' errors originally added in 
#                  parse_esl_epn_translate_startstop_outfile() and then 
#                  later updated in main() (step 11) by making their error
#                  messages more informative to include information on
#                  the corrected stop position which is calculated in this
#                  function.
#
# Arguments: 
#  $mdl_info_HAR:           REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:           REF to hash of arrays with information on the sequences, PRE-FILLED
#  $mdl_results_AAHR:       REF to results AAH, ADDED TO HERE
#  $err_ftr_instances_AHHR: REF to error instances AHH, PRE-FILLED with at least trc and ext errors
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:                REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
################################################################# 
sub results_calculate_corrected_stops {
  my $sub_name = "results_calculate_corrected_stops()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $err_ftr_instances_AHHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;

  # total counts of things
  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences

  # keys in the 3rd dim of @{$mdl_results_AAHR}
  my $start_key     = "p_start";
  my $stop_key      = "p_stop";
  my $strand_key    = "p_strand";
  my $len_key       = "len";
  my $cumlen_key    = "cumlen";
  my $trc_err_key   = "trc_err_flag";
  my $ext_err_key   = "ext_err_flag";
  my $prv_trc_key   = "prv_trc_flag";
  my $corr_stop_key = "c_stop";

  # lengths
  my $newlen;       # the new, corrected length of a prediction
  my $cumlen;       # the stored, cumulative length of a prediction, filled in results_calculate_lengths()
  my $remaininglen; # the current remaining len (<= $newlen) of the new length 
  my $len_corr;     # correction to the length due to trc or ext error, will be negative for trc, positive for ext

  # model indices
  my $first_mdl_idx; # the first model index that models the current feature 
  my $final_mdl_idx; # the final model index that models the current feature 

  # loop through all features and sequences, and correct stop
  # stop predictions for all trc and ext errors
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { 
      # we only deal with features for which annot_type is "model" here, we
      # deal with features with annot_type eq "multifeature" in ftr_results_calculate()
      my $is_matpept = ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp") ? 1 : 0;
      for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
        my $seq_name = $seq_info_HAR->{"seq_name"}[$seq_idx];
        ###########################################
        # block that handles potential stp errors #
        ###########################################
        if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name}) { 
          if($err_ftr_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name} ne "maybe") { 
            DNAORG_FAIL(sprintf("ERROR in $sub_name, stp error with non-maybe value %s for ftr %s seq_name: $seq_name", 
                                $err_ftr_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name}, $ftr_info_HAR->{"out_short"}[$ftr_idx]),
                        1, $FH_HR);
          }
          $final_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$ftr_idx];
          my $stp_err_stop_codon = fetchStopCodon($sqfile, $seq_name, 
                                                  $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"p_stop"}, 
                                                  $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"p_strand"}, $FH_HR);
          if(validateStopCodon($stp_err_stop_codon)) { 
            # it's a valid stop, remove the "maybe"
            error_instances_remove_maybe($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "stp", $seq_name, $FH_HR);
          }
          else { 
            # it is not a valid stop, the error stays and we update the error message
            my $updated_stp_errmsg = sprintf("%s ending at position %d on %s strand", $stp_err_stop_codon, 
                                             $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"p_stop"}, $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"p_strand"}); 
            error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "stp", $seq_name, $updated_stp_errmsg, $FH_HR);
          }
        }        
        # only one of trc or ext can exist, they're 'incompatible' in err_info_HAR
        ###########################################
        # block that handles potential trc errors #
        ###########################################
        if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name}) { 
          $len_corr = $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name};
          if($len_corr >= 0) { 
            DNAORG_FAIL("ERROR in $sub_name, trc error with non-negative correction $len_corr exists for ftr: $ftr_idx seq_name: $seq_name", 1, $FH_HR);
          }
          # determine the new length of the prediction, this is the old length plus the correction (which is negative)
          $first_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx];
          $final_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$ftr_idx];
          # first_mdl_idx should be the first model for which we have a prediction
          while(! exists $mdl_results_AAHR->[$first_mdl_idx][$seq_idx]{"p_start"}) { 
            $first_mdl_idx++; 
            if($first_mdl_idx > $final_mdl_idx) { 
              DNAORG_FAIL(sprintf("ERROR in $sub_name, can't determine first model for feature $ftr_idx (%s) for sequence $seq_name. trc error exists but no models have predictions for this feature.", $ftr_info_HAR->{"out_tiny"}[$ftr_idx]), 1, $FH_HR);
            }
          }
          # final_mdl_idx should be the final model for which we have a prediction
          while(! exists $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"p_start"}) { 
            $final_mdl_idx--; 
            if($final_mdl_idx < $first_mdl_idx) { 
              DNAORG_FAIL(sprintf("ERROR in $sub_name, can't determine final model for feature $ftr_idx (%s) for sequence $seq_name. trc error exists but no models have predictions for this feature.", $ftr_info_HAR->{"out_tiny"}[$ftr_idx]), 1, $FH_HR);
            }
          }

          $newlen = $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{$cumlen_key} + $len_corr;
          # careful: the length correction is with respect to the full feature (potentially multiple models (e.g. exons))
          # so it can be as high as the cumulative length of all model predictions, so we need to adjust the length 
          # by subtracting the correction from the cumulative length, not from the model length
          if($newlen <= 0) { 
            DNAORG_FAIL(sprintf("ERROR in $sub_name, trc error has length correction ($len_corr) longer than predicted length (%d) for ftr: $ftr_idx seq_name: $seq_name", 
                                $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{$cumlen_key}), 1, $FH_HR);
          }
          # determine which model to set corrected stop position for, and set it
          # again, this is only necessary because our correction is with respect to the 
          # feature, and if this feature is multi-model it means we need to figure out
          # which model the corrected stop is going to affect. 
          # 
          # For example: imagine feature 1 is comprised of models 1
          # and 2. Model 1 prediction is 1..300 model 2 prediction is
          # 601..1200.  If the length correction is 660, this means
          # that model 2 is now irrelevant (truncation occurs before
          # it's predicted start, and model 1's prediction should
          # change to 1..240.
          #
          $remaininglen = $newlen;
          for(my $mdl_idx = $first_mdl_idx; $mdl_idx <= $final_mdl_idx; $mdl_idx++) { 
            if(! exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$len_key}) { 
              DNAORG_FAIL("ERROR in $sub_name, results_AAHR->[$mdl_idx][$seq_idx]{$len_key} does not exist, but it should.", 1, $FH_HR); 
            }
            if(! exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key}) { 
              DNAORG_FAIL("ERROR in $sub_name, results_AAHR->[$mdl_idx][$seq_idx]{$len_key} does not exist, but it should.", 1, $FH_HR); 
            }
            $cumlen = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key};
            if($newlen <= $cumlen) { # this is the model to update
              if($mdl_results_AAHR->[$mdl_idx][$seq_idx]{$strand_key} eq "+") { 
                $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$corr_stop_key} = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$start_key} + $remaininglen - 1;
              }
              else { # negative strand
                $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$corr_stop_key} = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$start_key} - $remaininglen + 1;
              }
              # update "len" and "cumlen" and set $trc_err_key
              my $exon_len_corr = ($mdl_results_AAHR->[$mdl_idx][$seq_idx]{$len_key} - $remaininglen);
              $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key} -= $exon_len_corr;
              $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$len_key}     = $remaininglen;
              $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$trc_err_key} = 1;
              my $cur_cumlen = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key};

              # update the errmsg in @{$err_ftr_instances_AHHR} based on 
              # what just figured out about this truncated stop
              my $updated_trc_errmsg = sprintf("homology search predicted %d..%d", 
                                           create_output_start_and_stop($mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"},
                                                                        $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_stop"},
                                                                        $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR));
              if($ftr_info_HAR->{"nmodels"}[$ftr_idx] != 1) { 
                $updated_trc_errmsg .= sprintf(" %s %d of %d", ($is_matpept) ? "segment" : "exon", $mdl_idx - $first_mdl_idx + 1, $ftr_info_HAR->{"nmodels"}[$ftr_idx]);
              }
              $updated_trc_errmsg .= sprintf(" revised to %d..%d (stop shifted %d nt)", 
                                            create_output_start_and_stop($mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"},
                                                                         $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"c_stop"},
                                                                         $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR), 
                                            $exon_len_corr);
              error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_info_HAR->{"seq_name"}[$seq_idx], $updated_trc_errmsg, $FH_HR);
              
              # now for all remaining models for this feature, set $cumlen_key to $cur_cumlen, 
              # and $prv_trc_key to '1' indicating that that model prediction 
              # doesn't really exist due to an early stop in a previous model
              $mdl_idx++;
              while($mdl_idx <= $final_mdl_idx) { 
                if(! exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key}) { 
                  DNAORG_FAIL("ERROR in $sub_name, results_AAHR->[$mdl_idx][$seq_idx]{$len_key} does not exist, but it should.", 1, $FH_HR); 
                }
                $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key}  = $cur_cumlen;
                $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$prv_trc_key} = 1;
                $mdl_idx++;
              }
              # after the above while($mdl_idx) loop, $mdl_idx == $final_mdl_idx+1, so this means 
              # we will break out of the enclosing 'for(my $mdl_idx = $first_mdl_idx; $mdl_idx <= $final_mdl_idx; $mdl_idx++)' loop at 
              # the beginning of the next iteration
            } # end of 'if($newlen <= $cumlen)'
            $remaininglen -= $cumlen; # update $remaininglen
          } # end of 'for(my $mdl_idx' loop
        } # end of 'if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name}'
        #################################
        # block that handles ext errors #
        #################################
        elsif(exists $err_ftr_instances_AHH[$ftr_idx]{"ext"}{$seq_name}) { 
          $len_corr = $err_ftr_instances_AHH[$ftr_idx]{"ext"}{$seq_name};
          if($len_corr <= 0) { 
            DNAORG_FAIL("ERROR in $sub_name, ext error with non-positive correction $len_corr exists for ftr: $ftr_idx seq_name: $seq_name", 1, $FH_HR);
          }
          # ext are easier, always modify stop of final model
          my $mdl_idx = $ftr_info_HAR->{"final_mdl"}[$ftr_idx];
          if(! exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$len_key}) { 
            DNAORG_FAIL("ERROR in $sub_name, results_AAHR->[$mdl_idx][$seq_idx]{$len_key} does not exist, but it should.", 1, $FH_HR); 
          }
          if(! exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key}) { 
            DNAORG_FAIL("ERROR in $sub_name, results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key} does not exist, but it should.", 1, $FH_HR); 
          }
          if($mdl_results_AAHR->[$mdl_idx][$seq_idx]{$strand_key} eq "+") { 
            $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$corr_stop_key} = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$stop_key} + $len_corr;
          }
          else { # negative strand
            $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$corr_stop_key} = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$stop_key} - $len_corr;
          }
          # update "len" and "cumlen"
          $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$len_key}    += $len_corr;
          $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$cumlen_key} += $len_corr;
          $mdl_results_AAHR->[$mdl_idx][$seq_idx]{$ext_err_key} = 1;

          # update the errmsg in @{$err_ftr_instances_AHHR} based on 
          # what we just figured out about this extended stop
          my $updated_ext_errmsg = sprintf("homology search predicted %d..%d", 
                                           create_output_start_and_stop($mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"},
                                                                        $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_stop"},
                                                                        $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR));
          if($ftr_info_HAR->{"nmodels"}[$ftr_idx] != 1) { 
            $updated_ext_errmsg .= sprintf(" %s %d of %d", ($is_matpept) ? "segment" : "exon", $mdl_idx - $ftr_info_HAR->{"first_mdl"}[$ftr_idx] + 1, $ftr_info_HAR->{"nmodels"}[$ftr_idx]);
          }
          $updated_ext_errmsg .= sprintf(" revised to %d..%d (stop shifted %d nt)", 
                                         create_output_start_and_stop($mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"},
                                                                      $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"c_stop"},
                                                                      $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR), 
                                         $len_corr);
          error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "ext", $seq_info_HAR->{"seq_name"}[$seq_idx], $updated_ext_errmsg, $FH_HR);
        }
      } # end of loop over $seq_idx
    } # end of if $ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model"
  } # end of loop over $ftr_idx
  return;
}


#################################################################
# Subroutine:  results_calculate_overlaps_and_adjacencies()
# Incept:      EPN, Sat Mar 12 10:40:29 2016
#
# Purpose:    For all results in @{$mdl_results_AAH}, determine which
#             predictions for the same sequence are adjacent to 
#             and overlap with each other and save those strings
#             in $mdl_results_AAHR. For those that are 
#             inconsistent with the reference, create the 
#             appropriate 'olp', 'aja', or 'ajb' error in 
#             \@err_ftr_instances_AHH, 
#
#             Checks for and adds the following error codes for features
#             with "annot_type" eq "model":
#
#             "olp": adds with error message depicting deviation from reference 
#             "ajb": adds with error message depicting deviation from reference 
#             "aja": adds with error message depicting deviation from reference 
#
# Arguments: 
#  $mdl_info_HAR:       REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:       REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:       REF to hash of arrays with information on the sequences, ADDED TO HERE
#  $mdl_results_AAHR:   REF to results AAH, ADDED TO HERE
#  $err_ftr_instances_AHHR: REF to error instances AHH, PRE-FILLED with at least trc and ext errors
#  $err_info_HAR:       REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:            REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:              REF to hash of file handles
#
# Returns:    void
#
################################################################# 
sub results_calculate_overlaps_and_adjacencies { 
  my $sub_name = "results_calculate_overlaps_and_adjacencies()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $err_ftr_instances_AHHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $ftr_idx; # counter over features
  my $mdl_idx; # counter over models
  my $seq_idx; # counter over sequences
  
  # for each sequence, fill a temporary array with starts, stops and strands
  # then send it to overlapsAndAdjacenciesHelper() to get the adjacency and
  # overlap strings
  for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name = $seq_info_HAR->{"seq_name"}[$seq_idx];
    my @start_A  = ();
    my @stop_A   = ();
    my @strand_A = ();
    for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
      if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"}) { 
        push(@start_A,  $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"});
        push(@stop_A,   (exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"c_stop"} ? 
                         $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"c_stop"} : 
                         $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_stop"}));
        push(@strand_A, $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_strand"});
      }
      else { 
        push(@start_A,  -1);
        push(@stop_A,   -1);
        push(@strand_A, "");
      }
    }
    my @idx_ajb_str_A = (); # [0..$nmdl-1] string of mdl indices describing 'before' adjacencies for each model
    my @idx_aja_str_A = (); # [0..$nmdl-1] string of mdl indices describing 'after'  adjacencies for each model
    my @idx_olp_str_A = (); # [0..$nmdl-1] string of mdl indices describing overlaps for each model
    my @out_ajb_str_A = (); # [0..$nmdl-1] string of mdl descriptions describing 'before' adjacencies for each model
    my @out_aja_str_A = (); # [0..$nmdl-1] string of mdl descriptions describing 'after'  adjacencies for each model
    my @out_olp_str_A = (); # [0..$nmdl-1] string of mdl descriptions describing overlaps for each model
    overlapsAndAdjacenciesHelper($mdl_info_HAR, \@start_A, \@stop_A, \@strand_A, 
                                 $seq_info_HAR->{"seq_len"}[$seq_idx], $seq_info_HAR->{"accn_len"}[$seq_idx], 
                                 \@idx_ajb_str_A, \@idx_aja_str_A, \@idx_olp_str_A, 
                                 \@out_ajb_str_A, \@out_aja_str_A, \@out_olp_str_A, 
                                 $opt_HHR, $FH_HR);
    

    # populate @mdl_results_AAHR, and keep track of per-feature error messages
    my @ftr_olp_err_msg_A = ();
    my @ftr_ajb_err_msg_A = ();
    my @ftr_aja_err_msg_A = ();
    # initialize
    for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      $ftr_olp_err_msg_A[$ftr_idx] = "";
      $ftr_ajb_err_msg_A[$ftr_idx] = "";
      $ftr_aja_err_msg_A[$ftr_idx] = "";
    }
    for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
      my $ftr_idx = $mdl_info_HAR->{"map_ftr"}[$mdl_idx];
      $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"idx_ajb_str"} = $idx_ajb_str_A[$mdl_idx];
      $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"idx_aja_str"} = $idx_aja_str_A[$mdl_idx];
      $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"idx_olp_str"} = $idx_olp_str_A[$mdl_idx];
      $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"out_ajb_str"} = $out_ajb_str_A[$mdl_idx];
      $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"out_aja_str"} = $out_aja_str_A[$mdl_idx];
      $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"out_olp_str"} = $out_olp_str_A[$mdl_idx];

      # construct ajb err message
      if($idx_ajb_str_A[$mdl_idx] ne $mdl_info_HAR->{"idx_ajb_str"}[$mdl_idx]) { 
        my @diff_A = ();
        compareTwoOverlapOrAdjacencyIndexStrings($mdl_info_HAR->{"idx_ajb_str"}[$mdl_idx], 
                                                 $idx_ajb_str_A[$mdl_idx], 
                                                 $nmdl-1,
                                                 \@diff_A, $FH_HR);
        for(my $i = 0; $i < $nmdl; $i++) { 
          if($diff_A[$i] != 0) { 
            $ftr_ajb_err_msg_A[$ftr_idx] .= sprintf("%s%s(%s,%s)", 
                                                    ($ftr_ajb_err_msg_A[$ftr_idx] eq "") ? "" : ",", # need to add a comma only if we're appending
                                                    ($diff_A[$i] eq "-1") ? "-" : "+",              # is it a lost or added adjacency?
                                                    $mdl_info_HAR->{"out_idx"}[$mdl_idx], $mdl_info_HAR->{"out_idx"}[$i]);
          }
        }
      }
      # construct aja err message
      if($idx_aja_str_A[$mdl_idx] ne $mdl_info_HAR->{"idx_aja_str"}[$mdl_idx]) { 
        my @diff_A = ();
        compareTwoOverlapOrAdjacencyIndexStrings($mdl_info_HAR->{"idx_aja_str"}[$mdl_idx], 
                                                 $idx_aja_str_A[$mdl_idx], 
                                                 $nmdl-1,
                                                 \@diff_A, $FH_HR);
        for(my $i = 0; $i < $nmdl; $i++) { 
          if($diff_A[$i] != 0) { 
            $ftr_aja_err_msg_A[$ftr_idx] .= sprintf("%s%s(%s,%s)", 
                                                    ($ftr_aja_err_msg_A[$ftr_idx] eq "") ? "" : ",", # need to add a comma only if we're appending
                                                    ($diff_A[$i] eq "-1") ? "-" : "+",              # is it a lost or added adjacency?
                                                    $mdl_info_HAR->{"out_idx"}[$mdl_idx], $mdl_info_HAR->{"out_idx"}[$i]);
          }
        }
      }

      # construct olp err message
      if($idx_olp_str_A[$mdl_idx] ne $mdl_info_HAR->{"idx_olp_str"}[$mdl_idx]) { 
        my @diff_A = ();
        compareTwoOverlapOrAdjacencyIndexStrings($mdl_info_HAR->{"idx_olp_str"}[$mdl_idx], 
                                                 $idx_olp_str_A[$mdl_idx], 
                                                 $nmdl-1,
                                                 \@diff_A, $FH_HR);
        for(my $i = 0; $i < $nmdl; $i++) { 
          if($diff_A[$i] != 0) { 
            $ftr_olp_err_msg_A[$ftr_idx] .= sprintf("%s%s(%s,%s)", 
                                                    ($ftr_olp_err_msg_A[$ftr_idx] eq "") ? "" : ",", # need to add a comma only if we're appending
                                                    ($diff_A[$i] eq "-1") ? "-" : "+",              # is it a lost or added adjacency?
                                                    $mdl_info_HAR->{"out_idx"}[$mdl_idx], $mdl_info_HAR->{"out_idx"}[$i]);
          }
        }
      }
    } # end of 'for(mdl_idx'
    for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if($ftr_olp_err_msg_A[$ftr_idx] ne "") { 
        error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "olp", $seq_info_HAR->{"seq_name"}[$seq_idx], $ftr_olp_err_msg_A[$ftr_idx], $FH_HR);
      }
      if($ftr_ajb_err_msg_A[$ftr_idx] ne "") { 
        error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "ajb", $seq_info_HAR->{"seq_name"}[$seq_idx], $ftr_ajb_err_msg_A[$ftr_idx], $FH_HR);
      }
      if($ftr_aja_err_msg_A[$ftr_idx] ne "") { 
        error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "aja", $seq_info_HAR->{"seq_name"}[$seq_idx], $ftr_aja_err_msg_A[$ftr_idx], $FH_HR);
      }
    }
  } # end of 'for(my $seq_idx"
  return;
}

#################################################################
# Subroutine:  mdl_results_add_b5e_b5u_errors
# Incept:      EPN, Fri Jan 13 13:26:49 2017
#
# Purpose:    Report 'b5e' and 'b5u'. The reporting of these
#             is decoupled from 'b3e' and 'b3u' which are reported
#             later, because we need to find 'b5e' errors before
#             we look for 'trc' errors because a b5e precludes 
#             a trc.
#
#             Checks for and adds or updates the following error 
#             codes for features with "annot_type" eq "model":
#             
#             "b5e": adds this error, predicted hit not flush with 
#                    model end but flush with sequence end on 5'
#
#             "b5u": adds this error, predicted hit not flush with 
#                    model end and not flush with sequence end on 5'
#
# Arguments: 
#  $sqfile:                 REF to Bio::Easel::SqFile object, open sequence file containing sequences
#  $mdl_info_HAR:           REF to hash of arrays with information on the models, PRE-FILLED
#  $seq_info_HAR:           REF to hash of arrays with information on the sequences, ADDED TO HERE
#  $mdl_results_AAHR:       REF to model results AAH, ADDED TO HERE
#  $err_ftr_instances_AHHR: REF to error instances AHH, PRE-FILLED with at least trc and ext errors
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:                REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies: If we have predicted start and stop coordinates that don't make sense
#       given the lengths of the accession and the sequence we searched.
#
################################################################# 
sub mdl_results_add_b5e_b5u_errors { 
  my $sub_name = "mdl_results_add_b5e_b5u_errors";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $err_ftr_instances_AHHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $mdl_idx; # counter over models
  my $seq_idx; # counter over sequences

  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    my $ftr_idx    = $mdl_info_HAR->{"map_ftr"}[$mdl_idx];
    my $is_first   = $mdl_info_HAR->{"is_first"}[$mdl_idx];
    my $is_final   = $mdl_info_HAR->{"is_final"}[$mdl_idx];
    for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
      my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
      my $accn_len  = $seq_info_HAR->{"accn_len"}[$seq_idx];
      my $seq_len   = $seq_info_HAR->{"seq_len"}[$seq_idx];
      my $mdl_results_HR = \%{$mdl_results_AAHR->[$mdl_idx][$seq_idx]}; # for convenience

      if(exists $mdl_results_HR->{"p_start"}) { 
        if($mdl_results_HR->{"p_5overhang"} != 0) { 
          # and add b5e or b5u error
          if($mdl_results_HR->{"p_5seqflush"} == 1) { # b5e
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "b5e", $seq_name, $mdl_results_HR->{"p_5overhang"} . " nt from 5' end", $FH_HR);
            # special case for which we will not allow a trc error later
          }
          else { # $mdl_results_HR->{"p_5seqflush"} == 0, b5u
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "b5u", $seq_name, $mdl_results_HR->{"p_5overhang"} . " nt from 5' end", $FH_HR);
          }
        }
      } # end of 'else' entered if we have a prediction
    } # end of 'for($seq_idx' loop
  } # end of 'for($mdl_idx' loop
  return;
}      

#################################################################
# Subroutine:  mdl_results_add_str_nop_ost_errors
# Incept:      EPN, Thu Mar 31 13:43:58 2016
#
# Purpose:    Report 'str', 'nop', 'ost', 'b3e', and 'b3u' errors
#             and fill in the following keys in $mdl_results_AAHR:
#             "out_5boundary" and "out_3boundary".
#
#             Checks for and adds or updates the following error 
#             codes for features with "annot_type" eq "model":
#             
#             "nop": adds this error, no model prediction
#             
#             "str": updates this error, originally added in 
#                  parse_esl_epn_translate_startstop_outfile()
#                  by making the error message more informative
#                  to include position and codon of predicted start
#
#             "ost": adds this error, model prediction is on incorrect strand
#
# Arguments: 
#  $sqfile:                 REF to Bio::Easel::SqFile object, open sequence file containing sequences
#  $mdl_info_HAR:           REF to hash of arrays with information on the models, PRE-FILLED
#  $seq_info_HAR:           REF to hash of arrays with information on the sequences, ADDED TO HERE
#  $mdl_results_AAHR:       REF to model results AAH, ADDED TO HERE
#  $err_ftr_instances_AHHR: REF to error instances AHH, PRE-FILLED with at least trc and ext errors
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:                REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies: If we have predicted start and stop coordinates that don't make sense
#       given the lengths of the accession and the sequence we searched.
#
################################################################# 
sub mdl_results_add_str_nop_ost_b3e_b3u_errors { 
  my $sub_name = "mdl_results_add_str_nop_ost_b3e_b3u_errors()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $mdl_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $err_ftr_instances_AHHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $mdl_idx; # counter over models
  my $seq_idx; # counter over sequences

  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    my $ftr_idx    = $mdl_info_HAR->{"map_ftr"}[$mdl_idx];
    my $is_first   = $mdl_info_HAR->{"is_first"}[$mdl_idx];
    my $is_final   = $mdl_info_HAR->{"is_final"}[$mdl_idx];
    for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
      my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
      my $accn_len  = $seq_info_HAR->{"accn_len"}[$seq_idx];
      my $seq_len   = $seq_info_HAR->{"seq_len"}[$seq_idx];
      my $mdl_results_HR = \%{$mdl_results_AAHR->[$mdl_idx][$seq_idx]}; # for convenience

      if(! exists $mdl_results_HR->{"p_start"}) { 
        ############################
        # nop error (no prediction)
        ############################
        error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "nop", $seq_name, "", $FH_HR);
      }        
      else { # $mdl_results_HR->{"p_start"} exists, we have a prediction
        ##########################
        # update str err message #
        ##########################
        if($is_first && (exists $err_ftr_instances_AHHR->[$ftr_idx]{"str"}{$seq_name})) { 
          my $out_start = undef;
          ($out_start, undef) = create_output_start_and_stop($mdl_results_HR->{"p_start"}, $mdl_results_HR->{"p_stop"}, 
                                                             $accn_len, $seq_len, $FH_HR);
          my $updated_str_errmsg = sprintf("%s starting at position %d on strand %s", 
                                           fetchStartCodon($sqfile, $seq_name, $mdl_results_HR->{"p_start"}, $mdl_results_HR->{"p_strand"}, $FH_HR), 
                                           $out_start, $mdl_results_HR->{"p_strand"});
          error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "str", $seq_name, $updated_str_errmsg, $FH_HR);
          $mdl_results_HR->{"str_err_flag"} = 1;
        }
        #######################################
        # set out_5boundary and out_3boundary #
        #######################################
        # update p_5overhang and p_3overhang
        # at this point, $mdl_results_HR->{"p_5overhang"} and $results->{"p_3overhang"} are both integers >= 0
        # 5' block
        if($mdl_results_HR->{"p_5overhang"} == 0) { 
          $mdl_results_HR->{"out_5boundary"} = ".";
        }
        else { 
          if($mdl_results_HR->{"p_5overhang"} > 9) { 
            $mdl_results_HR->{"out_5boundary"} = "+"; 
          }
          else { 
            $mdl_results_HR->{"out_5boundary"} = $mdl_results_HR->{"p_5overhang"};
          }
        }
        # 3' block
        # check for 2 special cases, a trc or ext error, which invalidate the b3e or b3u values
        if(exists $mdl_results_HR->{"trc_err_flag"}) { 
          $mdl_results_HR->{"out_3boundary"} = "t";
        }
        elsif(exists $mdl_results_HR->{"ext_err_flag"}) { 
          $mdl_results_HR->{"out_3boundary"} = "e";
        }
        elsif($mdl_results_HR->{"p_3overhang"} == 0) { 
          $mdl_results_HR->{"out_3boundary"} = ".";
        }
        else { 
          if($mdl_results_HR->{"p_3overhang"} > 9) { 
            $mdl_results_HR->{"out_3boundary"} = "+"; 
          }
          else { 
            $mdl_results_HR->{"out_3boundary"} = $mdl_results_HR->{"p_3overhang"};
          }

          # b3e or b3u error
          if($mdl_results_HR->{"p_3seqflush"} == 1) { # b3e
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "b3e", $seq_name, $mdl_results_HR->{"p_3overhang"} . " nt from 3' end", $FH_HR);
            # special case for which we will not allow a trc error later
          }
          else { # $mdl_results_HR->{"p_3seqflush"} == 0, b3u
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "b3u", $seq_name, $mdl_results_HR->{"p_3overhang"} . " nt from 3' end", $FH_HR);
          }

        }
        if($mdl_results_HR->{"p_strand"} ne $mdl_info_HAR->{"ref_strand"}[$mdl_idx]) { 
          ################################
          # ost error (incorrect strand) #
          ################################
          error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "ost", $seq_name, "", $FH_HR);
        }
      } # end of 'else' entered if we have a prediction
    } # end of 'for($seq_idx' loop
  } # end of 'for($mdl_idx' loop
  return;
}      


#################################################################
# Subroutine:  mdl_results_calculate_out_starts_and_stops
# Incept:      EPN, Thu Mar 31 13:43:51 2016
#
# Purpose:    Compute the 'out_start', 'out_stop', and 'out_stop_codon'
#             values for all model results. These are the nucleotide
#             positions that will be output as annotations.  They may
#             differ from the predicted and corrected (p_start/p_stop
#             or c_start/c_stop) model result values which are in
#             internal coordinate space.  If -c is used, for example,
#             the genome is circular, and the internal coordinate
#             space ranges 1..seqlen*2 but our output coordinate space
#             will be -seqlen..seqlen.
#
# Arguments: 
#  $sqfile:             REF to Bio::Easel::SqFile object, open sequence file containing sequences
#  $mdl_info_HAR:       REF to hash of arrays with information on the models, PRE-FILLED
#  $seq_info_HAR:       REF to hash of arrays with information on the sequences, ADDED TO HERE
#  $mdl_results_AAHR:   REF to model results AAH, ADDED TO HERE
#  $opt_HHR:            REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:              REF to hash of file handles
#
# Returns:    void
#
# Dies: If we have predicted start and stop coordinates that don't make sense
#       given the lengths of the accession and the sequence we searched.
#
################################################################# 
sub mdl_results_calculate_out_starts_and_stops { 
  my $sub_name = "mdl_results_calculate_out_starts_and_stops";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $mdl_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $mdl_idx; # counter over models
  my $seq_idx; # counter over sequences

  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      my $accn_len  = $seq_info_HAR->{"accn_len"}[$seq_idx];
      my $seq_len   = $seq_info_HAR->{"seq_len"}[$seq_idx];
      my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
      my $mdl_results_HR = \%{$mdl_results_AAHR->[$mdl_idx][$seq_idx]}; # for convenience

      if(exists $mdl_results_HR->{"p_start"}) { 
        ($mdl_results_HR->{"out_start"}, $mdl_results_HR->{"out_stop"}) = 
            create_output_start_and_stop($mdl_results_HR->{"p_start"},
                                         exists($mdl_results_HR->{"c_stop"}) ? $mdl_results_HR->{"c_stop"} : $mdl_results_HR->{"p_stop"}, 
                                         $accn_len, $seq_len, $FH_HR);
        $mdl_results_HR->{"out_stop_codon"} = fetchStopCodon($sqfile, $seq_name,
                                                             (defined $mdl_results_HR->{"c_stop"}) ? $mdl_results_HR->{"c_stop"} : $mdl_results_HR->{"p_stop"}, 
                                                             $mdl_results_HR->{"p_strand"}, $FH_HR);

      }
    }
  }

  return;
}

#################################################################
# Subroutine:  mdl_results_compare_to_genbank_annotations
# Incept:      EPN, Thu Mar 31 13:52:15 2016
#
# Purpose:    Compute the 'num_genbank_mdl_annot' values for 
#             the %{$seq_info_HAR} and the 'genbank_mdl_annot_match' 
#             for the %{$mdl_results_AAHR}. These are the number
#             annotated features in Genbank, and whether or not
#             one of our annotated features matches an annotation
#             in GenBank, respectively.
#
# Arguments: 
#  $mdl_info_HAR:       REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:       REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:       REF to hash of arrays with information on the sequences, ADDED TO HERE
#  $mdl_results_AAHR:   REF to model results AAH, ADDED TO HERE
#  $cds_tbl_HHAR:       REF to CDS hash of hash of arrays, PRE-FILLED
#  $mp_tbl_HHAR:        REF to mature peptide hash of hash of arrays, can be undef, else PRE-FILLED
#  $opt_HHR:            REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:              REF to hash of file handles
#
# Returns:    void
#
# Dies: never
#
################################################################# 
sub mdl_results_compare_to_genbank_annotations { 
  my $sub_name = "mdl_results_compare_to_genbank_annotations()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $cds_tbl_HHAR, $mp_tbl_HHAR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $mdl_idx;
  my $seq_idx;

  # first, determine how many genbank annotations we have per sequence
  for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    # determine how many annotations CDS or mature peptide annotations there are in GenBank
    my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
    my $accn_len  = $seq_info_HAR->{"accn_len"}[$seq_idx];
    ($seq_info_HAR->{"num_genbank_mdl_annot"}[$seq_idx], $seq_info_HAR->{"num_genbank_mdl_exon_annot"}[$seq_idx]) = 
        count_genbank_annotations((opt_IsUsed("--matpept", $opt_HHR)) ? $mp_tbl_HHAR->{$accn_name} : $cds_tbl_HHAR->{$accn_name}, 
                                  $accn_len, $opt_HHR, $FH_HR);
  }

  # now check each of our annotations against the GenBank annotations
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    my $ftr_idx    = $mdl_info_HAR->{"map_ftr"}[$mdl_idx];
    my $is_matpept = ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp") ? 1 : 0;
    for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
      my $accn_len  = $seq_info_HAR->{"accn_len"}[$seq_idx];
      my $seq_len   = $seq_info_HAR->{"seq_len"}[$seq_idx];
      my $mdl_results_HR = \%{$mdl_results_AAHR->[$mdl_idx][$seq_idx]}; # for convenience

      if(exists $mdl_results_HR->{"p_start"}) { 
        # check our final annotation against existing annotation in GenBank 
        my $start  = $mdl_results_HR->{"out_start"};
        my $stop   = $mdl_results_HR->{"out_stop"};
        my $strand = $mdl_results_HR->{"p_strand"};
        my $tbl_HAR = ($is_matpept) ? $mp_tbl_HHAR->{$accn_name} : $cds_tbl_HHAR->{$accn_name};
        # TODO: if we get into viruses with many features, will want to do the comparison 
        # in sorted order by start to avoid quadratic time behavior
        if(defined $tbl_HAR) { # if there was annotation for this sequence 
          $mdl_results_HR->{"genbank_mdl_annot_match"} = compare_to_genbank_annotation($start, $stop, $strand, $accn_len, $seq_len, $tbl_HAR, $opt_HHR, $FH_HR);
        }
        else { # annotation doesn't exist, so we don't have a match
          $mdl_results_HR->{"genbank_mdl_annot_match"} = 0;
        }
      } # end of 'else' entered if we have a prediction
    } # end of 'for($seq_idx' loop
  } # end of 'for($mdl_idx' loop

  return;
}      

#################################################################
# Subroutine:  ftr_results_calculate
# Incept:      EPN, Tue Mar 15 06:26:38 2016
#
# Purpose:    Fill the feature results in
#             @{$ftr_results_AAHR}, by filling the following 3rd
#             dim keys for 'cds-mp' feature types: 
#             "out_start", "out_stop", "out_len", "out_start_codon", 
#             "out_stop_codon".
#
#             Checks for and adds the following error codes:
#             "nm3": for all features (annot_type eq "model" or "multifeature")
#             "stp": for features with annot_type eq "multifeature" & type eq "cds-mp"
#             "inp": for features with annot_type eq "multifeature" & type eq "cds-mp"
#             "int": for features with annot_type eq "multifeature" & type eq "cds-mp"
#             "aji": for features with annot_type eq "multifeature" & type eq "cds-mp"
#             "ntr": for features that are children of a multifeature/cds-mp feature
#
# Arguments: 
#  $sqfile:             REF to Bio::Easel::SqFile object, open sequence file containing sequences
#  $mdl_info_HAR:       REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:       REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:       REF to hash of arrays with information on the sequences, ADDED TO HERE
#  $ftr_results_AAHR:   REF to feature results AAH, ADDED TO HERE
#  $mdl_results_AAHR:   REF to finalized model results AAH, PRE-FILLED
#  $cds_tbl_HHAR:       REF to CDS hash of hash of arrays, PRE-FILLED
#  $err_ftr_instances_AHHR: REF to error instances AHH, ADDED TO HERE (potentially)
#  $err_info_HAR:       REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:            REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:              REF to hash of file handles
#
# Returns:    void
#
# Dies: If something unexpected happens: 
#       - an int errmsg exists before it should
#       - an ext error exists with a non-maybe error message before it should
#       - we have trouble dealing with an ext error 
#
################################################################# 
sub ftr_results_calculate { 
  my $sub_name = "ftr_results_calculate()";
  my $nargs_exp = 11;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $ftr_results_AAHR, $mdl_results_AAHR, $cds_tbl_HHAR, $err_ftr_instances_AHHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $ftr_idx;   # counter over features
  my $mdl_idx;   # counter over models
  my $seq_idx;   # counter over sequences
  my $seq_name;  # name of one sequence
  my $seq_len;   # length of the sequence, possibly including doubling
  my $accn_name; # name of accession
  my $accn_len;  # length of the accession, never includes doubling
  my $append_len = 0; # length of appended region
  my $ntr_errmsg = undef; # error message for an ntr error

  # foreach annot_type:multifeature and type:'cds-mp' feature, 
  # determine 'out_start', 'out_stop',
  # 'out_len', 'out_start_codon' and 'out_stop_codon'.
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") &&
       ($ftr_info_HAR->{"type"}[$ftr_idx]       eq "cds-mp")) { 

      # get the primary children array
      my @primary_children_idx_A = (); # feature indices of the primary children of this feature
      getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "primary", \@primary_children_idx_A, $FH_HR);
      my $np_children = scalar(@primary_children_idx_A);

      # get the all children array
      my @all_children_idx_A = (); # feature indices of the primary children of this feature
      getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "all", \@all_children_idx_A, $FH_HR);
      my $na_children = scalar(@all_children_idx_A);

      for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
        $seq_name   = $seq_info_HAR->{"seq_name"}[$seq_idx];
        $seq_len    = $seq_info_HAR->{"seq_len"}[$seq_idx];
        $accn_name  = $seq_info_HAR->{"accn_name"}[$seq_idx];
        $accn_len   = $seq_info_HAR->{"accn_len"}[$seq_idx];

        ($seq_info_HAR->{"num_genbank_ftr_annot"}[$seq_idx], $seq_info_HAR->{"num_genbank_ftr_exon_annot"}[$seq_idx]) = 
            count_genbank_annotations($cds_tbl_HHAR->{$accn_name}, $accn_len, $opt_HHR, $FH_HR);

        # set the str_err_flag, if nec
        if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"str"}{$seq_name}) { 
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"str_err_flag"} = 1;
        }

        # initialize our error-related variables
        my $aji_errmsg          = "";    # filled if we find two mature peptides that should be adjacent but are not
        my $set_start           = 0;     # set to '1' once we've seen the first model with an annotated hit
        my $inp_errmsg          = "";    # a list of model names ("out_tiny") for which we do not have annotations (nop errors), if any
        my $int_errmsg          = "";    # a list of model names ("out_tiny") which are not translated due to trc (ntr errors), if any
        my $cds_out_start       = undef; # start position to output for this CDS 
        my $cds_out_stop        = undef; # stop  position to output for this CDS 
        my $cds_fetch_start     = undef; # start position to *fetch* for this CDS' start codon
        my $cds_fetch_stop      = undef; # stop  position to *fetch* for this CDS' stop codon
        my $start_strand        = undef; # strand start codon is on
        my $stop_strand         = undef; # strand stop codon is on
        my $cds_len             = 0;     # length to output for this CDS
        my $child_had_trc       = 0;     # if we find a child with a trc error, set this to 1

        # step through all primary children of this feature
        for(my $child_idx = 0; $child_idx < $np_children; $child_idx++) { 
          my $child_ftr_idx = $primary_children_idx_A[$child_idx];

          for(my $child_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$child_ftr_idx]; $child_mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]; $child_mdl_idx++) { 
            if(! $child_had_trc) { # if $child_had_trc is true, then we dealt with this feature completely below
              # check to make sure we have a hit annotated for this model
              if(! exists $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_start"}) { 
                if($inp_errmsg ne "") { 
                  $inp_errmsg .= ", ";
                }
                $inp_errmsg .= $mdl_info_HAR->{"out_tiny"}[$child_mdl_idx];
              }
              else { # we do have a hit for $child_mdl_idx in $seq_idx
                $cds_len += $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"len"};
                
                if(! $set_start) { # first model, set cds_out_start
                  $cds_out_start   = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_start"};
                  $cds_fetch_start = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_start"};
                  $start_strand    = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_strand"};
                  $set_start = 1;
                }
                
                # check if we have a trc in this child model, and deal with it if we do
                if(exists $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"trc_err_flag"}) { 
                  $child_had_trc = 1;
                  ####################################################
                  # We have a trc in this child model $child_mdl_idx 
                  # 
                  # Determine the cds stop position to print ($cds_out_stop)
                  # and position in the sequence to fetch the stop ($cds_fetch_stop).
                  # 
                  # If and only if we have all mature peptides up to this point and they 
                  # are all adjacent (that is, we don't have a aji or inp error), 
                  # then this CDS is truncated. We do the following:
                  #
                  # - update the trc errmsg for this CDS in @{$err_ftr_instances_AHHR} 
                  #   based on what just figured out about this truncated stop
                  #   (if we don't have a str error for this CDS)
                  # - set ntr errors for all remaining children 
                  # - set int error for this CDS
                  # - break the loop over all children (we're done with this CDS)
                  ####################################################
                  
                  $cds_out_stop      = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_stop"};
                  $cds_fetch_stop    = (defined $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"c_stop"}) ? 
                      $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"c_stop"} :
                      $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_stop"};
                  $stop_strand    = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_strand"}; 
                  
                  if($inp_errmsg eq "" && $aji_errmsg eq "") { 
                    if(! exists $err_ftr_instances_AHHR->[$ftr_idx]{"str"}{$seq_name}) { 
                      # if we don't have a str error, update the trc error message
                      # use the final childs stop prediction as the predicted stop, if it exists
                      my $final_child_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$primary_children_idx_A[$np_children-1]];
                      my $updated_trc_errmsg  = "";
                      if(exists $mdl_results_AAHR->[$final_child_mdl_idx][$seq_idx]{"p_stop"}) { 
                        my $cds_pred_stop = $mdl_results_AAHR->[$final_child_mdl_idx][$seq_idx]{"p_stop"};
                        $updated_trc_errmsg = sprintf("homology search predicted %d..%d revised to %d..%d (stop shifted %d nt due to early stop in %s)", 
                                                      create_output_start_and_stop($cds_fetch_start, $cds_pred_stop, $accn_len, $seq_len, $FH_HR),
                                                      create_output_start_and_stop($cds_fetch_start, $cds_out_stop,  $accn_len, $seq_len, $FH_HR),
                                                      abs($cds_fetch_stop - $cds_pred_stop), $mdl_info_HAR->{"out_tiny"}[$child_mdl_idx]);
                      }
                      else { 
                        $updated_trc_errmsg = sprintf("homology search predicted %d..? revised to %d..%d (due to early stop in %s)", 
                                                      $cds_out_start, 
                                                      create_output_start_and_stop($cds_out_start,   $cds_out_stop,  $accn_len, $seq_len, $FH_HR),
                                                      $mdl_info_HAR->{"out_tiny"}[$child_mdl_idx]);
                      }
                      
                      if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name}) { 
                        # update it
                        error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_name, $updated_trc_errmsg, $FH_HR);
                      }
                      else { 
                        # it doesn't yet exist, so add the trc error. This
                        # is rare, but we may have no trc error for this
                        # CDS yet, if the predicted child mature peptide
                        # sequences didn't all exist, then we won't have
                        # combined those predictions into a predicted CDS,
                        # and thus we didn't check that predicted CDS for
                        # trc errors in
                        # parse_esl_epn_translate_startstop_outfile().
                        error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_name, $updated_trc_errmsg, $FH_HR);
                      }
                      # set the trc_err_flag for this feature
                      $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"trc_err_flag"} = 1;
                    }
                    # all remaining children get a 'ntr' error,
                    # and the CDS gets an 'int' error, which we need
                    # to build the error message for
                    $ntr_errmsg = sprintf("early stop in mature peptide %d ending at position %d", $child_ftr_idx+1, $cds_out_stop);
                    if($int_errmsg ne "") { 
                      DNAORG_FAIL("ERROR in $sub_name, setting int errmsg for ftr_idx: $ftr_idx due to 'trc', but it is not blank", 1, $FH_HR);;
                    }
                    $child_idx++;
                    my $ntr_err_ct = 0;
                    # for all remaining children: throw 'ntr' and append to 'int' err message
                    while($child_idx < $np_children) {
                      $child_ftr_idx = $primary_children_idx_A[$child_idx];
                      if(! exists $err_ftr_instances_AHHR->[$child_ftr_idx]{"nop"}{$seq_name}) { 
                        error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "ntr", $seq_name, $ntr_errmsg, $FH_HR);
                        $ntr_err_ct++;
                        if($int_errmsg ne "") { 
                          $int_errmsg .= ", ";
                        }
                        $int_errmsg .= $ftr_info_HAR->{"out_tiny"}[$child_ftr_idx];
                      }
                      else { # nop for at least one mdl for this feature
                        if($int_errmsg ne "") { 
                          $int_errmsg .= ", ";
                        }
                        $int_errmsg .= $ftr_info_HAR->{"out_tiny"}[$child_ftr_idx] . "(nop)";
                      }
                      $child_idx++;
                    }
                    if($ntr_err_ct > 0) { 
                      # we set at least one ntr error for mature peptides, set int for this CDS
                      error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "int", $seq_name, $int_errmsg, $FH_HR);
                    }
                    # now $child_idx is $np_children, so this breaks the 'for(child_idx' loop
                  } # end of 'if($inp_errmsg eq "" && $aji_errmsg eq "")
                  else { 
                    # we had a trc in one of the children, but it didn't trigger
                    # a trc in the CDS because we already have an adjacency error
                    # and/or an inp error (one of the mature peptides was not predicted)
                    # which means we cannot confidently say the CDS is truncated.
                    # remove the CDS' trc error IF it exists (it may not if 
                    # there wasn't an early stop due to a frameshift or something
                    # before or after the trc in the child mature peptide).
                    if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name}) { 
                      error_instances_remove_not_maybe($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_name, $FH_HR);
                    }
                  }
                  # deal with a potential stp error, if we have one
                  
                } # end of 'if trc_err_flag'
                else { 
                  ###########################################################
                  # we do not have a trc in this child model $child_mdl_idx #
                  ###########################################################
                  # check if we're adjacent to the next feature, we only need to 
                  # do this if we're the final model for the current feature and
                  # we're not the final child feature
                  if(($child_idx < ($np_children-1)) && 
                     ($child_mdl_idx == $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx])) { 
                    my $nxt_child_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$primary_children_idx_A[($child_idx+1)]];
                    # check if they are adjacent 
                    if(! checkForIndexInOverlapOrAdjacencyIndexString($mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"idx_aja_str"}, $nxt_child_mdl_idx, $FH_HR)) { 
                      if($aji_errmsg ne "") { $aji_errmsg .= ", "; }
                      $aji_errmsg .= sprintf("%s (%s..%s) not adjacent to %s (%s..%s)", 
                                             $mdl_info_HAR->{"out_tiny"}[$child_mdl_idx], 
                                             (defined $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_start"}) ? $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_start"} : "unknown", 
                                             (defined $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_stop"})  ? $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_stop"}  : "unknown", 
                                             $mdl_info_HAR->{"out_tiny"}[$nxt_child_mdl_idx], 
                                             (defined $mdl_results_AAHR->[$nxt_child_mdl_idx][$seq_idx]{"out_start"}) ? $mdl_results_AAHR->[$nxt_child_mdl_idx][$seq_idx]{"out_start"} : "unknown", 
                                             (defined $mdl_results_AAHR->[$nxt_child_mdl_idx][$seq_idx]{"out_stop"})  ? $mdl_results_AAHR->[$nxt_child_mdl_idx][$seq_idx]{"out_stop"}  : "unknown");
                    }        
                  } # end of 'if($child_idx < ($np_children-1))'

                  # if we are the final model of the final child feature, determine the stop position/strand
                  if(($child_idx == ($np_children-1)) && 
                     ($child_mdl_idx == $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx])) { 
                    # the final child, determine the stop position/strand, and deal with ext error if we have one
                    $stop_strand = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_strand"}; 
                    # we should have to append the stop codon
                    # (GenBank annotation of final mature peptides doesn't include the stop codon,
                    #  so it's not covered in our homology model and we have to take special care
                    #  to annotate it.
                    $append_len = 0; 
                    if(exists $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"append_stop"}) { 
                      my ($out_append_start, $out_append_stop) = 
                          create_output_start_and_stop($mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"append_start"}, 
                                                       $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"append_stop"},
                                                       $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR);
                      $cds_out_stop = $out_append_stop;
                      $cds_fetch_stop = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"append_stop"};
                      # and update the cds_len (this is the final child, so this will only happen once)
                      $append_len = abs($mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"append_stop"} - $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"append_start"}) + 1; 
                      $cds_len += $append_len;
                    }
                  }
                }
              } # end of 'else' entered if we don't have a trc error
            } # end of 'if(! $child_had_trc'
          } # end of 'for(my $child_mdl_idx..'
        }  # end of 'for(my $child_idx..'

        #######################################
        # deal with a stp error, if we have one
        #######################################
        if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name}) { 
          # value should be "maybe" 
          if($err_ftr_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name} ne "maybe") { 
            DNAORG_FAIL("ERROR in $sub_name, ext error with non-maybe value for ftr: $ftr_idx seq_name: $seq_name", 1, $FH_HR);
          }
          if((defined $cds_fetch_stop) && (defined $cds_out_stop)) { 
            my $stp_err_stop_codon = fetchStopCodon($sqfile, $seq_name, $cds_fetch_stop, $stop_strand, $FH_HR);
            if(validateStopCodon($stp_err_stop_codon)) { 
              # it's a valid stop, remove the "maybe"
              error_instances_remove_maybe($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "stp", $seq_name, $FH_HR);
            }
            else { 
              # invalid stop, update the error message
              my $updated_stp_errmsg = sprintf("%s ending at position %d on %s strand", $stp_err_stop_codon, $cds_out_stop, $stop_strand); 
              error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "stp", $seq_name, $updated_stp_errmsg, $FH_HR);
            }
          }
          else { 
            # we can't define the stop position, so this is also a stp error with only a "?" error message
            error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "stp", $seq_name, "?", $FH_HR);
          }
        }

        ####################################################################################################
        # deal with an ext error, if we have one and we don't have any adjacency or nop errors in children
        # it's important we do this *after* handling (either validating or removing) a potential 'stp' error
        # in the block again, because we redfine $cds_out_stop, $cds_fetch_stop and $cds_len here.
        ####################################################################################################
        if((exists $err_ftr_instances_AHHR->[$ftr_idx]{"ext"}{$seq_name}) && 
           ($aji_errmsg eq "") && ($inp_errmsg eq "") && (defined $cds_out_stop)) { 
          my $len_corr = $err_ftr_instances_AHH[$ftr_idx]{"ext"}{$seq_name};
          if($len_corr <= 0) { 
            DNAORG_FAIL("ERROR in $sub_name, ext error with non-positive correction $len_corr exists for ftr: $ftr_idx seq_name: $seq_name", 1, $FH_HR);
          }
          my $final_child_ftr_idx       = $primary_children_idx_A[$np_children-1]; # first model for final MP that makes up this CDS
          my $final_first_child_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$final_child_ftr_idx]; # first model for final MP that makes up this CDS
          my $final_final_child_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$final_child_ftr_idx]; # final model for final MP that makes up this CDS
          # make sure everything we assume exists, actually does, if not, there's a coding error somewhere
          if(! exists $mdl_results_AAHR->[$final_first_child_mdl_idx][$seq_idx]{"p_start"}) { 
            DNAORG_FAIL("ERROR in $sub_name, results_AAHR->[$final_first_child_mdl_idx][$seq_idx]{p_start} does not exist, but it should.", 1, $FH_HR); 
          }
          if(! exists $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"len"}) { 
            DNAORG_FAIL("ERROR in $sub_name, results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{len} does not exist, but it should.", 1, $FH_HR); 
          }
          if(! exists $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"cumlen"}) { 
            DNAORG_FAIL("ERROR in $sub_name, results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{cumlen} does not exist, but it should.", 1, $FH_HR); 
          }
          if(! exists $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"append_stop"}) { 
            DNAORG_FAIL("ERROR in $sub_name, results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{append_stop} does not exist, but it should.", 1, $FH_HR); 
          }
          if($mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"p_strand"} eq "+") { 
            $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"c_stop"} = $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"p_stop"} + ($len_corr - $append_len);
            $cds_fetch_stop = $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"c_stop"} + $append_len;
          }
          else { # negative strand
            $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"c_stop"} = $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"p_stop"} - ($len_corr - $append_len);
            $cds_fetch_stop = $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"c_stop"} - $append_len;
          }
          ($mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"out_start"}, $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"out_stop"}) = 
              create_output_start_and_stop($mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"p_start"}, 
                                           $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"c_stop"}, 
                                           $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR);
          # only the final model is affected by an ext error
          $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"len"}    += ($len_corr - $append_len);
          $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"cumlen"} += ($len_corr - $append_len);
          $cds_len += ($len_corr - $append_len);
          
          # get first part of cds error message using current $cds_out_stop, if it exists
          my $updated_cds_ext_errmsg = sprintf("homology search predicted %d..%d", 
                                               create_output_start_and_stop($cds_out_start, $cds_out_stop,
                                                                           $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR));
          # get first part of mp error message 
          my $mp_ext_errmsg = sprintf("homology search predicted %d..%d", 
                                      create_output_start_and_stop($mdl_results_AAHR->[$final_first_child_mdl_idx][$seq_idx]{"p_start"},
                                                                   $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"p_stop"},
                                                                   $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR));
          # now recompute $cds_out_stop

          (undef, $cds_out_stop) = 
              create_output_start_and_stop($cds_fetch_start, # this is irrelevant due to first undef arg
                                           $cds_fetch_stop,
                                           $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR);

          # get second part of CDS error message
          $updated_cds_ext_errmsg .= sprintf(" revised to %d..%d (stop shifted %d nt)", $cds_out_start, $cds_out_stop, $len_corr-$append_len);
          # get second part of MP error message
          $mp_ext_errmsg .= sprintf(" revised to %d..%d (stop shifted %d nt)", 
                                    create_output_start_and_stop($mdl_results_AAHR->[$final_first_child_mdl_idx][$seq_idx]{"p_start"}, 
                                                                 $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"c_stop"}, 
                                                                 $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR), 
                                    $len_corr-$append_len);
          error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "ext", $seq_info_HAR->{"seq_name"}[$seq_idx], $updated_cds_ext_errmsg, $FH_HR);
          error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $final_child_ftr_idx, "ext", $seq_info_HAR->{"seq_name"}[$seq_idx], $mp_ext_errmsg, $FH_HR);
          # and finally, update the results out_3boundary value for this model
          $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"out_3boundary"} = "e";
        }

        # sanity check
        if((defined $start_strand) && (defined $stop_strand) && ($start_strand ne $stop_strand)) { 
          DNAORG_FAIL(sprintf("ERROR in $sub_name, feature $ftr_idx %s for sequence $seq_name has start and stop on different strands", $ftr_info_HAR->{"out_short"}[$ftr_idx]), 1, $FH_HR);
        }

        ###############################################################
        # if we did not find a child with a trc, and we have 
        # a trc error for this CDS, it must be invalid and caused
        # by an aji, nm3, or b5e error or something like it,
        # so we remove it
        ###############################################################
        if(! $child_had_trc) { 
          if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name}) { 
            error_instances_remove_not_maybe($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_name, $FH_HR);
          }
        }

        # add the aji (adjacency) error if necessary
        if($aji_errmsg ne "") { 
          # adjacency error
          error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "aji", $seq_name, $aji_errmsg, $FH_HR);
        }

        # add the inp (INterrupted due to no Prediction) if necessary
        if($inp_errmsg ne "") { 
          error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "inp", $seq_name, $inp_errmsg, $FH_HR);
        }

        # set ftr_results, we can set start if $cds_out_start is defined, 
        if(defined $cds_out_start) { 
          my $start_codon = fetchStartCodon($sqfile, $seq_name, $cds_fetch_start, $start_strand, $FH_HR);
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_start"}       = $cds_out_start;
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_start_codon"} = $start_codon;
          if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"str"}{$seq_name}) { 
            my $updated_str_errmsg = sprintf("%s starting at position %d on strand %s", $start_codon, 
                                             $cds_out_start, $start_strand);
            error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "str", $seq_name, $updated_str_errmsg, $FH_HR);
          }
        }
        else { 
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_start"}       = "?";
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_start_codon"} = "?";
        }
        
        # for the stop and length attributes to be anything but "?" 
        # we require the following: 
        # - no aji error
        # - no inp error
        # - $cds_out_start defined
        # - $cds_out_stop defined
        # 
        if(($aji_errmsg eq "") && ($inp_errmsg eq "") && (defined $cds_out_start) && (defined $cds_out_stop)) { 
          # an aji or inp error means we can't really say where the stop is
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_stop"}       = $cds_out_stop;
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_stop_codon"} = fetchStopCodon($sqfile, $seq_name, $cds_fetch_stop, $stop_strand, $FH_HR);
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_len"}       = $cds_len;
          if(($cds_len % 3) != 0) { 
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "nm3", $seq_name, "$cds_len", $FH_HR);
          }
        }
        else { 
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_stop"}       = "?";
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_stop_codon"} = "?";
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_len"}        = "?";
        }
      
        # check if existing annotation for this CDS exists in %{$cds_tbl_HHAR}
        if((defined $cds_out_start) && (defined $cds_out_stop)) { 
          if(defined ($cds_tbl_HHAR->{$accn_name})) { # if there was annotation for this sequence 
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"genbank_ftr_annot_match"} = 
                compare_to_genbank_annotation($cds_out_start, $cds_out_stop, $start_strand,
                                              $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], 
                                              $cds_tbl_HHAR->{$accn_name}, $opt_HHR, $FH_HR);
          }
          else { # annotation doesn't exist, so we don't have a match
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"genbank_ftr_annot_match"} = 0;
          }
        }
        # one final step: if we have a 'trc' error for this CDS, check the 'all children' array, 
        # and throw 'ntr' errors for any mature peptides encoded by this CDS that are not
        # translated. We did this above for the primary peptides, but here we do it for any
        # non-primary peptides.
        if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name}) { 
          for(my $child_idx = 0; $child_idx < $na_children; $child_idx++) { 
            my $child_ftr_idx = $all_children_idx_A[$child_idx];
            if(! exists $err_ftr_instances_AHHR->[$child_ftr_idx]{"nop"}{$seq_name}) { # there is a prediction for all models for this feature
              my $final_child_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx];
              my $cur_stop = (defined $mdl_results_AAHR->[$final_child_mdl_idx][$seq_idx]{"c_stop"}) ? 
                  $mdl_results_AAHR->[$final_child_mdl_idx][$seq_idx]{"c_stop"} :
                  $mdl_results_AAHR->[$final_child_mdl_idx][$seq_idx]{"p_stop"};
              if(($stop_strand eq "+") && ($cds_fetch_stop < $cur_stop) && (! exists $err_ftr_instances_AHHR->[$child_ftr_idx]{"ntr"}{$seq_name})) { 
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "ntr", $seq_name, $ntr_errmsg, $FH_HR);
              }
              if(($stop_strand eq "-") && ($cds_fetch_stop > $cur_stop) && (! exists $err_ftr_instances_AHHR->[$child_ftr_idx]{"ntr"}{$seq_name})) { 
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "ntr", $seq_name, $ntr_errmsg, $FH_HR);
              }
            }
          }
        }
      } # end of 'for($seq_idx'
    }
  } # end of 'for($ftr_idx'

  # foreach annot_type:'model' feature, look for 
  # 'nm3' errors
  # it's important we do this at the end of this function,
  # after we've potentially redefined the length of the final
  # child features of a multifeature feature, which we do 
  # if we have an ext error.
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { 
      for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
        $seq_name   = $seq_info_HAR->{"seq_name"}[$seq_idx];
        $accn_name  = $seq_info_HAR->{"accn_name"}[$seq_idx];
        $accn_len   = $seq_info_HAR->{"accn_len"}[$seq_idx];
        my $cumlen = undef;
        # go through all models instead of just using the final one in case the final one is not annotated
        for($mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx]; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
          if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"cumlen"}) { 
            $cumlen = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"cumlen"};
          }
        }
        if(defined $cumlen && (($cumlen % 3) != 0)) { 
          # an nm3 error
          error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "nm3", $seq_name, "$cumlen", $FH_HR);
        }
      }
    }
  }

  return;
}       
#################################################################
# Subroutine : dump_results()
# Incept:      EPN, Tue Mar  1 14:54:27 2016
#
# Purpose:    Dump results data structure to $FH. Probably only
#             useful for debugging.
#
# Arguments: 
#  $FH:               file handle to output to
#  $mdl_results_AAHR: REF to results AAH, PRE-FILLED
#  $mdl_info_HAR:     REF to hash of arrays with information on the models, PRE-FILLED
#  $seq_info_HAR:     REF to hash of arrays with sequence information, PRE-FILLED
#  $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:            REF to hash of file handles
#
# Returns:    void
#
# Dies:       never
#
################################################################# 
sub dump_results {
  my $sub_name = "dump_results()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($FH, $mdl_results_AAHR, $mdl_info_HAR, $seq_info_HAR, $opt_HHR, $FH_HR) = @_;

  my $nmdl = scalar(@{$mdl_results_AAHR});
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR);
  
  for(my $m = 0; $m < $nmdl; $m++) { 
    for(my $s = 0; $s < $nseq; $s++) { 
      printf $FH ("model($m): %20s  seq($s): %20s  accn: %10s  len: %10d  ", 
                  $mdl_info_HAR->{"cmname"}[$m],
                  $seq_info_HAR->{"seq_name"}[$s],
                  $seq_info_HAR->{"accn_name"}->[$s],
                  $seq_info_HAR->{"accn_len"}[$s]);
      foreach my $key (sort keys %{$mdl_results_AAHR->[$m][$s]}) { 
        printf $FH ("%s: %s ", $key, $mdl_results_AAHR->[$m][$s]{$key});
      }
      printf $FH ("\n")
    }
  }
  return;
}

#################################################################
#################################################################
#  Subroutines related to origin sequences:
#    validate_origin_seq()
#    find_origin_sequences()
#    get_origin_output_for_sequence
#
#################################################################
# Subroutine: validate_origin_seq
# Incept:     EPN, Tue Mar 15 12:34:21 2016
#
# Synopsis:   Validate an origin sequence passed in
#             as <s> with --origin <s>. It should have 
#             a single '|' in it, which occurs 
#             just before what should be the first nt
#             of the genome. Return the origin offset:
#             the number of nts before the "|".
#
#             For example: "TAATATT|AC"
#             indicates that the final 7 nts of each
#             genome should be "TAATAAT" and the first
#             two should be "AC". In this case the origin
#             offset is 7.
#
# Args:
#   $origin_seq: the origin sequence
#
# Returns:  Origin offset, as explained in synopsis, above.
#
# Dies: if $origin_seq is incorrectly formatted. We die without
#       outputting to any file handles, because this function is 
#       called before the log and cmd files are set-up.
#
####################################
sub validate_origin_seq { 
  my $sub_name = "validate_origin_seq";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($origin_seq) = @_;

  if($origin_seq =~ m/[^ACGT\|\\]/) { 
    die "ERROR with -origin <s>, <s> can only contain characters A, C, G, T or \|";
  }

  my $origin_offset = index($origin_seq, "|");
  if($origin_offset == -1) { 
    die "ERROR with --origin <s>, <s> must contain a single | character immediately before the nucleotide that should be the first nt of the genome";
  }
  elsif($origin_offset < (length($origin_seq)-1)) { # if this isn't the final character of the string, look for another one
    my $second_offset = index($origin_seq, "|", $origin_offset+1);
    if($second_offset != -1) { 
      die "ERROR with --origin <s>, <s> must contain a single | character, $origin_seq has more than one";
    }
  }

  # printf("in $sub_name, $origin_seq returning $origin_offset\n");

  return $origin_offset;
}

#################################################################
# Subroutine: find_origin_sequences
# Incept:     EPN, Tue Mar 15 12:34:30 2016
# 
# Purpose:    Identify all exact occurrences of a sequence in a file
#             of sequences, and a string with the coordinates of the
#             matches in @{$seq_info_HAR->{$key}}. Each sequence
#             name is in @{$seq_info_HAR->{"seq_name"}}.
#
#             Checks for and adds the following error codes: "ori".
#
# Args:   
#  $sqfile:                 Bio::Easel::SqFile object, the sequence file to search in
#  $qseq:                   query sequence we're looking for
#  $seq_info_HAR:           REF to hash of arrays with information 
#                           on the sequences, PRE-FILLED
#  $err_seq_instances_HHR:  REF to the 2D hash of per-sequence errors, initialized here
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:                REF to 2D hash of option values, 
#                           see top of epn-options.pm for description
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies: if --origin option is not used according to %{$opt_HHR}
#       if -c option is not used according to %{$opt_HHR}
#
#################################################################
sub find_origin_sequences { 
  my $sub_name = "find_origin_sequences";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $qseq, $seq_info_HAR, $err_seq_instances_HHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;

  if(! opt_IsUsed("--origin", $opt_HHR)) { 
    DNAORG_FAIL("ERROR in $sub_name, --origin not used...."); 
  }
  if(! opt_IsOn("-c", $opt_HHR)) { 
    DNAORG_FAIL("ERROR in $sub_name, -c not used...."); 
  }

  # fetch each sequence and look for $qseq in it
  # (could (and probably should) make this more efficient...)
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    $seq_info_HAR->{"origin_coords_str"}[$seq_idx] = "";
    my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
    my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
    my $accn_len  = $seq_info_HAR->{"accn_len"}[$seq_idx];
    my $seq_len   = $seq_info_HAR->{"seq_len"}[$seq_idx];
    my $fasta_seq = $sqfile->fetch_seq_to_fasta_string($seq_name, -1);
    my ($header, $seq) = split(/\n/, $fasta_seq);
    my $nfound = 0; # number of occurrences of origin sequence per sequence
    chomp $seq;
    if($seq =~ m/\r$/) { chop $seq; } # remove ^M if it exists

    # now use Perl's index() function to find all occurrences of $qseq
    my $qseq_posn = index($seq, $qseq);
    # if $qseq_posn == -1, then no occurrences were found. In this case we don't store 
    # any entry in coords_HAR for this $accn. The caller needs to know what to do
    # if it sees no entry for an $accn

    while($qseq_posn != -1) { 
      $qseq_posn++;
      if($qseq_posn <= $accn_len) { # note: we've just incremented qseq_posn by 1 in prv line so now it is in 1..length($seq) coords, not 0..length($seq)-1
        my $qseq_start = $qseq_posn;
        my $qseq_stop  = $qseq_posn + length($qseq) - 1;
        # adjust coordinates so they're within 1..L
        ($qseq_start, $qseq_stop) = 
            create_output_start_and_stop($qseq_start, $qseq_stop, $accn_len, $seq_len, $FH_HR);
        if($seq_info_HAR->{"origin_coords_str"}[$seq_idx] ne "") { 
          $seq_info_HAR->{"origin_coords_str"}[$seq_idx] .= ",";
        }
        $seq_info_HAR->{"origin_coords_str"}[$seq_idx] .= $qseq_start . ":" . $qseq_stop;
        $nfound++;
      }
      if($qseq_posn > $accn_len) { 
        $qseq_posn = -1; 
        # this breaks the while loop because we're searching in a duplicated genome, 
        # and we're into the second copy, no need to keep searching the same seq
      }
      else { 
        $qseq_posn = index($seq, $qseq, $qseq_posn);
      }
    }
    # printf("in $sub_name, $seq_name nfound: $nfound\n");
    if($nfound != 1) { 
      error_instances_add(undef, $err_seq_instances_HHR, $err_info_HAR, -1, "ori", $seq_name, $nfound . " occurrences", $FH_HR);
    }
  }
  
  return;
}


#################################################################
# Subroutine: get_origin_output_for_sequence
# Incept:     EPN, Tue Mar 15 13:33:29 2016
#
# Synopsis:   For a given sequence index $seq_idx, determine 
#             the output strings related to the origin sequence
#
# Args: 
#  $seq_info_HAR:  REF to hash of arrays with information 
#                  on the sequences (including origins), PRE-FILLED
#  $seq_idx:       index of sequence we're working on
#  $origin_offset: offset of origin, 1st genome position within origin sequence
#  $FH_HR:         REF to hash of file handles
#
# Returns: 5 values:
#          $oseq_ct:       number of occurrences of origin sequence found
#          $oseq_start:    start position of origin seq if $oseq_ct == 1, else '-'
#          $oseq_stop:     stop  position of origin seq if $oseq_ct == 1, else '-'
#          $oseq_offset:   offset position of origin seq if $oseq_ct == 1, else '-'
#          $oseq_passfail: 'P' if $oseq_ct is 1, else 'F'
#
# Dies: if $seq_info_HAR->{"origin_coords_str"}[$seq_idx] does not exist.
# 
#################################################################
sub get_origin_output_for_sequence { 
  my $sub_name = "get_origin_output_for_sequence";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_info_HAR, $seq_idx, $origin_offset, $FH_HR) = @_;

  if(! exists $seq_info_HAR->{"origin_coords_str"}[$seq_idx]) { 
    DNAORG_FAIL("ERROR in $sub_name, no origin_coords_str (not even empty str) exists in seq_info_HAR for $seq_idx", 1, $FH_HR);
  }
  my $accn_len = $seq_info_HAR->{"accn_len"}[$seq_idx];

  # initializereturn values
  my $oseq_start    = "-"; # changed below if $oseq_ct == 1
  my $oseq_stop     = "-"; # changed below if $oseq_ct == 1
  my $oseq_offset   = "-"; # changed below if $oseq_ct == 1
  my $oseq_passfail = "F"; # changed below if $oseq_ct == 1

  my @coords_A = split(",", $seq_info_HAR->{"origin_coords_str"}[$seq_idx]);
  my $oseq_ct = scalar(@coords_A);

  if($oseq_ct == 1) { 
    ($oseq_start, $oseq_stop) = split(":", $coords_A[0]);
    # printf("HEYA in $sub_name, seq_idx: $seq_idx oseq_start: $oseq_start oseq_stop: $oseq_stop\n");
    $oseq_offset = ($oseq_start < 0) ? ($oseq_start + $origin_offset) : ($oseq_start + $origin_offset - 1);
    # printf("HEYA in $sub_name, seq_idx: $seq_idx oseq_offset: $oseq_offset\n");
    # $oseq_offset is now number of nts to shift origin in counterclockwise direction
    if($oseq_offset > ($accn_len / 2)) { # simpler (shorter distance) to move origin clockwise
      $oseq_offset = $accn_len - $oseq_offset; # note, we don't add 1 here
    }
    else { # simpler to shift origin in counterclockwise direction, we denote this as a negative offset
      $oseq_offset *= -1;
    }
    $oseq_passfail = "P";
  }

  return ($oseq_ct, $oseq_start, $oseq_stop, $oseq_offset, $oseq_passfail);
}



#################################################################
#################################################################
#  Subroutines related to the error instance data structures:
#    error_instances_initialize_AHH()
#    error_instances_add()
#    error_instances_update()
#    error_instances_remove_maybe()
#    error_instances_remove_not_maybe()
#    error_instances_validate_all()
#
#################################################################
# Subroutine:  error_instances_initialize_AHH()
# Incept:      EPN, Fri Mar  4 12:26:42 2016
#
# Purpose:    Initialize the error instances array of arrays of 
#             2 dimensional hashes. The array is [0..$f..$nftr-1] 
#             where $f is a feature index in %{ftr_info_HA}. 
#             The key in the 1st hash dimension is an error code, 
#             the key in the 2nd hash dimension is a sequence 
#             name (from array $seq_info_HAR{}).
#
# Arguments: 
#  $err_ftr_instances_AHHR: REF to the array of 2D hashes of per-feature errors, initialized here
#  $err_seq_instances_HHR:  REF to the 2D hash of per-sequence errors, initialized here
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $ftr_info_HAR:           REF to the feature info hash of arrays, PRE-FILLED
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies:       If err_info_HAR is not complete
#
#################################################################
sub error_instances_initialize_AHH { 
  my $sub_name = "error_instances_initialize_AHH";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_ftr_instances_AHHR, $err_seq_instances_HHR, $err_info_HAR, $ftr_info_HAR, $FH_HR) = @_;
  
  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $FH_HR); 
  my $nerr = validateErrorInfoHashIsComplete($err_info_HAR, undef, $FH_HR); 

  # the per-feature errors
  @{$err_ftr_instances_AHHR} = ();
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    %{$err_ftr_instances_AHHR->[$ftr_idx]} = (); 
    for(my $err_idx = 0; $err_idx < $nerr; $err_idx++) { 
      if($err_info_HAR->{"pertype"}[$err_idx] eq "feature") { 
        %{$err_ftr_instances_AHHR->[$ftr_idx]{$err_info_HAR->{"code"}[$err_idx]}} = ();
      }
    }
  }

  # the per-sequence errors
  %{$err_seq_instances_HHR} = ();
  for(my $err_idx = 0; $err_idx < $nerr; $err_idx++) { 
    if($err_info_HAR->{"pertype"}[$err_idx] eq "sequence") { 
      %{$err_seq_instances_HHR->{$err_info_HAR->{"code"}[$err_idx]}} = ();
    }
  }

  return;
}


#################################################################
# Subroutine:  error_instances_add()
# Incept:      EPN, Tue Mar  8 11:06:18 2016
#
# Purpose:    Add an $err_code error to the @{$err_ftr_instances_AHHR} 
#             for feature index $ftr_idx, sequence name $seq_name.
#
#             If $err_code is an error code for which $err_info_HAR->{"pertype"}[$err_idx]
#             is "sequence", then $err_ftr_instances_AHHR can be undef,
#             because we are updating only $err_seq_instances_HHR.
#
#             If $err_info_HAR->{"pertype"}[$err_idx] is "feature", then 
#             $err_seq_instances_HHR can be undef, because we are updating
#             only $err_ftr_instances_AHHR. 
#
#             $err_idx is derived from $err_code and %{$err_inf_HAR}
#             as follows: $err_info_HAR->{"code"}[$err_idx] ==
#             $err_code.
#             
# Arguments: 
#  $err_ftr_instances_AHHR: REF to per-feature error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "sequence".
#  $err_seq_instances_HHR:  REF to per-sequence error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "feature".
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $ftr_idx:                feature index, -1 if $err_code pertype ($err_info_HAR->{"pertype"}[$err_idx]) is 'sequence'
#  $err_code:               error code we're adding an error for, $err_idx is derived from this 
#  $seq_name:               sequence name
#  $value:                  value to add as $err_ftr_instances_AHHR->[$ftr_idx]{$code}{$seq_name}
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies:       - If we find an error instance incompatibility.
#             - If value is "maybe" but maybes are not allowed for $err_code in %err_info_HAR
#               or a value already exists for that instance
#             - If we already have a value for this $seq_idx+$err_code pair
#               (call error_instances_update() to update a value).
#             - if pertype of $err_code is "feature"  and $err_ftr_instances_AHHR is undef
#             - if pertype of $err_code is "sequence" and $err_seq_instances_HHR is undef
#             - if pertype of $err_code is "sequence" and $ftr_idx is not -1
#
#################################################################
sub error_instances_add { 
  my $sub_name = "error_instances_add()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_ftr_instances_AHHR, $err_seq_instances_HHR, $err_info_HAR, $ftr_idx, $err_code, $seq_name, $value, $FH_HR) = @_;

  # printf("in $sub_name, ftr_idx: $ftr_idx, err_code: $err_code, seq_name: $seq_name, value: $value\n");
  
  my $err_idx = findNonNumericValueInArray($err_info_HAR->{"code"}, $err_code, $FH_HR); 
  if($err_idx == -1) { 
    DNAORG_FAIL("ERROR in $sub_name, unrecognized error code $err_code", 1, $FH_HR);
  }
  # printf("in $sub_name err_code: $err_code err_idx: $err_idx\n");
  
  # check if it's the special "maybe" value, which is only allowed for some error codes 
  if($value eq "maybe") { 
    # make sure it's allowed
    if(! $err_info_HAR->{"maybe_allowed"}[$err_idx]) { 
      DNAORG_FAIL("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, trying to set to maybe but maybes are not allowed for this error code.", 1, $FH_HR);
    }
  }

  # determine if we're a per-feature or per-sequence error
  my $pertype = $err_info_HAR->{"pertype"}[$err_idx];

  if($pertype eq "feature") { 
    if(! defined $err_ftr_instances_AHHR) { 
      DNAORG_FAIL("ERROR in $sub_name error code $err_code is a per-feature error, but err_ftr_instances_AHHR is undefined", 1, $FH_HR);
    }
    # this error shouldn't already exist, unless it's already set to $value, in which case it's okay
    if(exists $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}) { 
      if($err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} ne $value) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, this error already exists as %s (maybe you want to use error_instances_update()?).", $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}), 1, $FH_HR);
      }
    }
    # set the value
    $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} = $value;
  }
  elsif($pertype eq "sequence") { 
    if(! defined $err_seq_instances_HHR) { 
      DNAORG_FAIL("ERROR in $sub_name error code $err_code is a per-sequence error, but err_seq_instances_HHR is undefined", 1, $FH_HR);
    }
    if($ftr_idx != -1) { 
      DNAORG_FAIL("ERROR in $sub_name error code $err_code is a per-sequence error, but passed in ftr_idx is not -1, but $ftr_idx", 1, $FH_HR);
    }
    # this error shouldn't already exist, unless it's already set to $value, in which case it's okay
    if(exists $err_seq_instances_HHR->{$err_code}{$seq_name}) { 
      if($err_seq_instances_HHR->{$err_code}{$seq_name} ne $value) { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, this error already exists as %s (maybe you want to use error_instances_update()?).", $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}), 1, $FH_HR);
      }
    }

    # set the value
    $err_seq_instances_HHR->{$err_code}{$seq_name} = $value; 
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unexpected pertype of $pertype for error $err_code", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine:  error_instances_update()
# Incept:      EPN, Tue Mar 15 14:26:58 2016
#
# Purpose:    Update the value for an already existing $err_code 
#             error in @{$err_ftr_instances_AHHR} 
#             for feature index $ftr_idx, sequence name $seq_name.
#
# Arguments: 
#  $err_ftr_instances_AHHR: REF to per-feature error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "sequence".
#  $err_seq_instances_HHR:  REF to per-sequence error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "feature".
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $ftr_idx:                REF to the feature info hash of arrays, PRE-FILLED
#  $err_code:               error code we're adding an error for
#  $seq_name:               sequence name
#  $value:                  value to add as $err_ftr_instances_AHHR->[$ftr_idx]{$code}{$seq_name}
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies:       - If we find an error instance incompatibility.
#             - If value is "maybe" but maybes are not allowed for $err_code in %err_info_HAR
#               or a value already exists for that instance
#             - If we do not already have a value for this $seq_idx+$err_code pair
#               (call error_instances_add() to add a new value).
#################################################################
sub error_instances_update { 
  my $sub_name = "error_instances_update()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_ftr_instances_AHHR, $err_seq_instances_HHR, $err_info_HAR, $ftr_idx, $err_code, $seq_name, $value, $FH_HR) = @_;
  
  # printf("in $sub_name, ftr_idx: $ftr_idx, err_code: $err_code, seq_name: $seq_name, value: $value\n");

  my $err_idx = findNonNumericValueInArray($err_info_HAR->{"code"}, $err_code, $FH_HR); 
  if($err_idx == -1) { 
    DNAORG_FAIL("ERROR in $sub_name, unrecognized error code $err_code", 1, $FH_HR);
  }
  
  # check if it's the special "maybe" value, which we don't allow 
  # an update to (must be added with error_instances_add())
  if($value eq "maybe") { 
    DNAORG_FAIL("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, trying to set to maybe (you may wnat to use error_instances_add()?).", 1, $FH_HR);
  }

  # determine if we're a per-feature or per-sequence error
  my $pertype = $err_info_HAR->{"pertype"}[$err_idx];

  if($pertype eq "feature") { 
    if(! defined $err_ftr_instances_AHHR) { 
      DNAORG_FAIL("ERROR in $sub_name error code $err_code is a per-feature error, but err_ftr_instances_AHHR is undefined", 1, $FH_HR);
    }
    # this error *should* already exist (we're updating it)
    if(! exists $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, this error does not already exist (maybe you want to use error_instances_add()?).", $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}), 1, $FH_HR);
    }

    # update the value
    $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} = $value;
  }
  elsif($pertype eq "sequence") { 
    if(! defined $err_seq_instances_HHR) { 
      DNAORG_FAIL("ERROR in $sub_name error code $err_code is a per-sequence error, but err_seq_instances_HHR is undefined", 1, $FH_HR);
    }
    if($ftr_idx != -1) { 
      DNAORG_FAIL("ERROR in $sub_name error code $err_code is a per-sequence error, but passed in ftr_idx is not -1, but $ftr_idx", 1, $FH_HR);
    }
    # this error *should* already exist (we're updating it)
    if(! exists $err_seq_instances_HHR->{$err_code}{$seq_name}) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, this error does not already exist (maybe you want to use error_instances_add()?).", $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}), 1, $FH_HR);
    }

    # set the value
    $err_seq_instances_HHR->{$err_code}{$seq_name} = $value; 
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unexpected pertype of $pertype for error $err_code", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine:  error_instances_remove_maybe()
# Incept:      EPN, Tue Mar  8 13:50:22 2016
#
# Purpose:    Remove an error $err_code from @{$err_ftr_instances_AHHR} 
#             for feature index $ftr_idx, sequence name $seq_name,
#             where the current value is "maybe".
#
#             Use error_instances_remove_maybe to remove "maybe" values
#             from error codes that allow maybes.
#
# Arguments: 
#  $err_ftr_instances_AHHR: REF to per-feature error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "sequence".
#  $err_seq_instances_HHR:  REF to per-sequence error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "feature".
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $ftr_idx:                REF to the feature info hash of arrays, PRE-FILLED
#  $err_code:               error code we're adding an error for
#  $seq_name:               sequence name
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies:       - if current value does not exist or is not "maybe"
#             - if $err_code does not allow maybe values
#             - if pertype of $err_code is "feature"  and $err_ftr_instances_AHHR is undef
#             - if pertype of $err_code is "sequence" and $err_seq_instances_HHR is undef
#             - if pertype of $err_code is "sequence" and $ftr_idx is not -1
#
#################################################################
sub error_instances_remove_maybe { 
  my $sub_name = "error_instances_remove_maybe()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_ftr_instances_AHHR, $err_seq_instances_HHR, $err_info_HAR, $ftr_idx, $err_code, $seq_name, $FH_HR) = @_;
  
  # printf("in $sub_name, ftr_idx: $ftr_idx, err_code: $err_code, seq_name: $seq_name\n");

  my $err_idx = findNonNumericValueInArray($err_info_HAR->{"code"}, $err_code, $FH_HR); 
  if($err_idx == -1) { 
    DNAORG_FAIL("ERROR in $sub_name, unrecognized error code $err_code", 1, $FH_HR);
  }
  if(! $err_info_HAR->{"maybe_allowed"}[$err_idx]) { 
    DNAORG_FAIL("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, trying to remove maybe but maybes are not allowed for this error code.", 1, $FH_HR);
  }

  # make sure the current value is "maybe" 

  # determine if we're a per-feature or per-sequence error
  my $pertype = $err_info_HAR->{"pertype"}[$err_idx];
  if($pertype eq "feature") { 
    if(! exists $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}) { 
      DNAORG_FAIL("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, no value exists.", 1, $FH_HR);
    }
    if($err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} ne "maybe") { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, value we are trying to remove is not \"maybe\" but rather %s, you may want to use error_instances_remove_not_maybe()", 
                          $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}), 1, $FH_HR); 
    } 
    # okay, it exists and is "maybe", remove it:
    delete $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name};
  }
  elsif($pertype eq "sequence") { 
    if(! exists $err_seq_instances_HHR->{$err_code}{$seq_name}) { 
      DNAORG_FAIL("ERROR in $sub_name, error code $err_code, seq_name $seq_name, no value exists.", 1, $FH_HR);
    }
    if($ftr_idx != -1) { 
      DNAORG_FAIL("ERROR in $sub_name error code $err_code is a per-sequence error, but passed in ftr_idx is not -1, but $ftr_idx", 1, $FH_HR);
    }
    if($err_seq_instances_HHR->{$err_code}{$seq_name} ne "maybe") { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, error code $err_code, seq_name $seq_name, value we are trying to remove is not \"maybe\" but rather %s, you may want to use error_instances_remove_not_maybe()", 
                          $err_seq_instances_HHR->{$err_code}{$seq_name}), 1, $FH_HR); 
    } 
    # okay, it exists and is "maybe", remove it:
    delete $err_seq_instances_HHR->{$err_code}{$seq_name};
  }

  return
}

#################################################################
# Subroutine:  error_instances_remove_not_maybe()
# Incept:      EPN, Thu Mar 17 09:25:53 2016
#
# Purpose:    Remove an error $err_code from @{$err_ftr_instances_AHHR} 
#             for feature index $ftr_idx, sequence name $seq_name for
#             which the value is not "maybe". 
#
#             Use a different function: error_instances_remove_maybe() 
#             to remove "maybe" values from error codes that allow maybes.
#
# Arguments: 
#  $err_ftr_instances_AHHR: REF to per-feature error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "sequence".
#  $err_seq_instances_HHR:  REF to per-sequence error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "feature".
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $ftr_idx:                REF to the feature info hash of arrays, PRE-FILLED
#  $err_code:               error code we're adding an error for
#  $seq_name:               sequence name
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies:       - if current value does not exist or is "maybe"
#             - if $err_code allows maybe values
#             - if pertype of $err_code is "feature"  and $err_ftr_instances_AHHR is undef
#             - if pertype of $err_code is "sequence" and $err_seq_instances_HHR is undef
#             - if pertype of $err_code is "sequence" and $ftr_idx is not -1
#
#################################################################
sub error_instances_remove_not_maybe { 
  my $sub_name = "error_instances_remove_not_maybe()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_ftr_instances_AHHR, $err_seq_instances_HHR, $err_info_HAR, $ftr_idx, $err_code, $seq_name, $FH_HR) = @_;
  
  # printf("in $sub_name, ftr_idx: $ftr_idx, err_code: $err_code, seq_name: $seq_name\n");

  my $err_idx = findNonNumericValueInArray($err_info_HAR->{"code"}, $err_code, $FH_HR); 
  if($err_idx == -1) { 
    DNAORG_FAIL("ERROR in $sub_name, unrecognized error code $err_code", 1, $FH_HR);
  }
  if($err_info_HAR->{"maybe_allowed"}[$err_idx]) { 
    DNAORG_FAIL("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, trying to remove a value for an error code where maybes are allowed, you may want to use error_intance_remove_maybe.", 1, $FH_HR);
  }

  # determine if we're a per-feature or per-sequence error
  my $pertype = $err_info_HAR->{"pertype"}[$err_idx];
  if($pertype eq "feature") { 
    if(! exists $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}) { 
      DNAORG_FAIL("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, no value exists.", 1, $FH_HR);
    }

    # just to be absolutely sure, check the value is not "maybe" (it shouldn't be since we checked that this err_idx does not
    # allow maybes ($err_info_HAR->{"maybe_allowed"} = 0), but it's possible it is (e.g. a newly written function may have
    # set it to maybe)
    if($err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} eq "maybe") { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name for ftr_idx $ftr_idx, error code $err_code, seq_name $seq_name, value we are trying to remove is \"maybe\", you may want to use error_instances_remove_maybe().", 
                          $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}), 1, $FH_HR); 
    } 
    # okay, it exists and is not "maybe", remove it:
    delete $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name};
  }
  elsif($pertype eq "sequence") { 
    if(! exists $err_seq_instances_HHR->{$err_code}{$seq_name}) { 
      DNAORG_FAIL("ERROR in $sub_name, error code $err_code, seq_name $seq_name, no value exists.", 1, $FH_HR);
    }
    if($ftr_idx != -1) { 
      DNAORG_FAIL("ERROR in $sub_name error code $err_code is a per-sequence error, but passed in ftr_idx is not -1, but $ftr_idx", 1, $FH_HR);
    }
    if($err_seq_instances_HHR->{$err_code}{$seq_name} eq "maybe") { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, error code $err_code, seq_name $seq_name, value we are trying to remove is \"maybe\", you may want use error_instances_remove_maybe().", 
                          $err_seq_instances_HHR->{$err_code}{$seq_name}), 1, $FH_HR); 
    } 
    # okay, it exists and is not "maybe" remove it:
    delete $err_seq_instances_HHR->{$err_code}{$seq_name};
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unexpected pertype of $pertype for error $err_code", 1, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine:  error_instances_validate_all()
# Incept:      EPN, Thu Mar 17 09:59:00 2016
#
# Purpose:    Given all errors in @{$err_ftr_instances_AHHR} and 
#             %{$err_seq_instances_HHR}, check for any combinations
#             of incompatible errors and check for any unfulfilled
#             required combinations of errors and die if any are
#             found. Also check for any "maybe" error messages,
#             and fail if any are found.
#
# Arguments: 
#  $err_ftr_instances_AHHR: REF to per-feature error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "sequence".
#  $err_seq_instances_HHR:  REF to per-sequence error instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $err_info_HAR->{"pertype"}[$err_idx] is "feature".
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $ftr_info_HAR:           REF to the feature info hash of arrays, PRE-FILLED
#  $seq_info_HAR:           REF to the sequence info hash of arrays, PRE-FILLED
#  $opt_HHR:                REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies:       - if we find an incompatible combination of errors for the
#               same sequence/feature pair
#             - if any errors exist without their required other errors for
#               the same sequence/feature pair
#             - if any "maybe" error messages exist
#
#################################################################
sub error_instances_validate_all { 
  my $sub_name = "error_instances_validate_all()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_ftr_instances_AHHR, $err_seq_instances_HHR, $err_info_HAR, $ftr_info_HAR, $seq_info_HAR, $opt_HHR, $FH_HR) = @_;
  
  my $dnaorg_fail_errmsg = ""; # filled as we find incompatibilities

  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $nerr = validateErrorInfoHashIsComplete   ($err_info_HAR, undef, $FH_HR); # nerr: number of errors
  
  my $ftr_idx;
  my $err_idx;
  my $seq_idx;
  my $seq_name;
  my $err_code;

  # validate that no "maybe" values exist
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    for($err_idx = 0; $err_idx < $nerr; $err_idx++) { 
      $err_code = $err_info_HAR->{"code"}[$err_idx];
      foreach $seq_name (keys %{$err_ftr_instances_AHHR->[$ftr_idx]{$err_code}}) { 
        if((exists $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}) && 
           ($err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}) eq "maybe") { 
          if(! $err_info_HAR->{"maybe_allowed"}[$err_idx]) { 
            $dnaorg_fail_errmsg .= sprintf("ERROR in $sub_name, value maybe exists for ftr %s seq %s error $err_code, but maybes not allowed for this error\n", 
                                           $ftr_info_HAR->{"out_tiny"}[$ftr_idx], $seq_name);
          }
          else { # maybes allowed for this error, but we shouldn't have any at this stage
            $dnaorg_fail_errmsg .= sprintf("ERROR in $sub_name, value maybe exists for ftr %s seq %s error $err_code (no maybes should be left at this stage)\n", 
                                           $ftr_info_HAR->{"out_tiny"}[$ftr_idx], $seq_name);
          }
        }
      }
    }
  }

  # validate that there are no incompatibilities
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    for($err_idx = 0; $err_idx < $nerr; $err_idx++) { 
      $err_code = $err_info_HAR->{"code"}[$err_idx];
      my @incompat_A = split(",", $err_info_HAR->{"incompat"}[$err_idx]);
      foreach my $incompat_err_idx (@incompat_A) { 
        my $incompat_err_code = $err_info_HAR->{"code"}[$incompat_err_idx];
        foreach $seq_name (keys %{$err_ftr_instances_AHHR->[$ftr_idx]{$err_code}}) { 
          if(exists $err_ftr_instances_AHHR->[$ftr_idx]{$incompat_err_code}{$seq_name}) { 
            if($incompat_err_idx <= $err_idx) { 
              # this way we only print an error once for each incompatibility b/t 'A and B' (not 'A and B' plus 'B and A')
              $dnaorg_fail_errmsg .= sprintf("ERROR in $sub_name, incompatible error combination $err_code and $incompat_err_code for ftr %s seq %s\n", 
                                             $ftr_info_HAR->{"out_tiny"}[$ftr_idx], $seq_name);
              
            }
          }
        }
      }
    }
  }

  # validate that all required combinations are met
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    for($err_idx = 0; $err_idx < $nerr; $err_idx++) { 
      $err_code = $err_info_HAR->{"code"}[$err_idx];
      my @requires_A = split(",", $err_info_HAR->{"requires"}[$err_idx]);
      foreach my $requires_err_idx (@requires_A) { 
        my $requires_err_code = $err_info_HAR->{"code"}[$requires_err_idx];
        foreach $seq_name (keys %{$err_ftr_instances_AHHR->[$ftr_idx]{$err_code}}) { 
          if(! exists $err_ftr_instances_AHHR->[$ftr_idx]{$requires_err_code}{$seq_name}) { 
            if($requires_err_idx <= $err_idx) { 
              # this way we only print an error once for each faield requirement b/t 'A and B' (not 'A and B' plus 'B and A')
              $dnaorg_fail_errmsg .= sprintf("ERROR in $sub_name, error $err_code exists without the required code $requires_err_code for ftr %s seq %s\n", 
                                             $ftr_info_HAR->{"out_tiny"}[$ftr_idx], $seq_name);
            }
          }
        }
      }
    }
  }

  if($dnaorg_fail_errmsg ne "") { 
    DNAORG_FAIL($dnaorg_fail_errmsg, 1, $FH_HR);
  }

  return;
}

#################################################################
#################################################################
#
#  Subroutines related to output:
#    output_tbl_explanations()
#    output_tbl_get_headings()
#    output_tbl_get_headings_helper()
#    output_tbl_get_headings_explanation_helper()
#    output_tbl_all_sequences()
#    output_tbl_page_of_sequences()
#    output_errors_header()
#    output_errors_all_sequences()
#    output_errors_summary()
#    output_multifeature_relationships()
#    output_gap_info()
#################################################################

#################################################################
# Subroutine:  output_tbl_explanations
# Incept:      EPN, Wed Mar 16 11:19:29 2016
#
# Purpose:    Output the explanatory text for the tabular output
#             files.
#
# Arguments: 
#  $out_header_exp_AR: ref to array of output explanation lines
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
################################################################# 
sub output_tbl_explanations { 
  my $sub_name = "output_tbl_explanations()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_header_exp_AR, $ofile_info_HHR) = @_;

  my $tbl_FH     = $ofile_info_HHR->{"FH"}{"tbl"};
  my $failtbl_FH = $ofile_info_HHR->{"FH"}{"failtbl"};
  my $errtbl_FH  = $ofile_info_HHR->{"FH"}{"errtbl"};

  foreach my $FH ($tbl_FH, $failtbl_FH, $errtbl_FH) { 
    foreach my $line (@{$out_header_exp_AR}) { 
      print $FH $line;
    }
    print $FH "#\n";
    outputDividingLine(undef, $FH); # undef makes outputDividingLine() use its default length for the dividing line
    print $FH "#\n";
  }
  return;
}

#################################################################
# Subroutine : output_tbl_get_headings()
# Incept:      EPN, Thu Mar 10 20:34:07 2016
#
# Purpose:     Fill header data structures with strings for headers
#              in tabular annotation output.
#
#             IMPORTANT: This function must stay in sync with the long
#             block of code in the main script entitled 'Pass through
#             all accessions, and gather and output annotation for
#             each'. Here we define the headers of the output, in the
#             main script we add output for each of those headers, so
#             they must stay in sync.
# Arguments: 
#  $out_row_header_AR:  ref to 1D array of row headers, FILLED HERE 
#                       iff output format is seq-as-cols
#  $out_header_exp_AR:  ref to 1D array of header explanations, each
#                       element is a line to be printed in explanatory
#                       section of the output; FILLED HERE
#  $mdl_info_HAR:       REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:       REF to hash of arrays with information on the features, PRE-FILLED
#  $opt_HHR:            REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:     REF to the 2D hash of output file information
# 
# Returns:     void
# 
# Dies: never
#
################################################################# 
sub output_tbl_get_headings { 
  my $sub_name = "output_tbl_get_headings";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($out_row_header_AR, $out_header_exp_AR, $mdl_info_HAR, $ftr_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;
  
  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience
  my $nmdl = validateModelInfoHashIsComplete  ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $FH_HR); # nftr: number of features

  # determine optional modes
  my $do_fid       = (opt_Get("--doalign", $opt_HHR)) ? 1 : 0; # '1' to do fid output
  my $do_totfid    = (opt_Get("--doalign", $opt_HHR)) ? 1 : 0; # '1' to do fid output
  my $do_ss3       = 1; # '1' to do ss3 output, '0' to skip it
  my $do_stop      = 1; # '1' to do stop output, '0' to skip it
  my $do_mdlb      = 1; # '1' to do model boundary output, '0' to skip it
  my $do_olap      = 1; # '1' to do overlap output, '0' to skip it
  my $do_exist     = (opt_Get("--infasta", $opt_HHR) || opt_Get("--tblnocomp", $opt_HHR)) ? 0 : 1; # '1' to do comparison to existing GenBank annotation, '0' to skip it
  my $do_fullolap  = 0; # '1' to output full overlap strings
  my $do_fulladj   = 0; # '1' to output full adjacency strings
  my $do_matpept   = (numNonNumericValueInArray($ftr_info_HAR->{"type"}, "mp", $FH_HR) > 0) ? 1 : 0;
  my $do_cds_notmp = (numNonNumericValueInArray($ftr_info_HAR->{"type"}, "cds-notmp", $FH_HR) > 0) ? 1 : 0;

  # miscellaneous variables
  my $width_result = 5 + $nmdl + 2; # an important width 
  my $row_div_char = ":"; # divides rows
  my $width;    # width of a field
  my $pad;      # string of all spaces used for pretty formatting
  my $tok1;     # first  level token (line 1 of column headers) 
  my $tok2;     # second level token (line 2 of column headers) 
  my $tok3;     # third  level token (line 3 of column headers) 
  my $tok4;     # fourth level token (line 4 of column headers) 
  my $tok5;     # fifth  level token (line 5 of column headers) 
  my $exp_tok1; # first  level explanation token, only used if has to be different from $tok1
  my $exp_tok4; # fourth level explanation token, only used if has to be different from $tok4
  my @pf_text_A = (); # array of lines to print about pass/fail strings to explanation at the end
  my $pf_idx = 1;     # pass/fail index
  my %need_to_define_H = (); # hash of terms we need to define in explanatory text

##  We store the row headers in a 1D array @{$out_row_header_AR}.
##  We have the same values as in the column headers, but each
##  level is concatenated together per row. Here's the row header
##  information that pertains to the example column header example
##  above.
#
#  idx
#  accession
#  totlen
#  origin sequence:#
#  origin sequence:start
#  origin sequence:stop
#  origin sequence:1stps
#  origin sequence:offst
#  origin sequence:PF
#  CDS #1 [single exon; +]:movement protein:start1
#  CDS #1 [single exon; +]:movement protein:stop1
#  CDS #1 [single exon; +]:movement protein:fid1
#  CDS #1 [single exon; +]:movement protein:md1
#  CDS #1 [single exon; +]:movement protein:length
#  CDS #1 [single exon; +]:movement protein:SS3
#  CDS #1 [single exon; +]:movement protein:stp
#  CDS #1 [single exon; +]:movement protein:PF
#  CDS #2 [2 exons; -]:replication associated protein:start1
#  CDS #2 [2 exons; -]:replication associated protein:stop1
#  CDS #2 [2 exons; -]:replication associated protein:fid1
#  CDS #2 [2 exons; -]:replication associated protein:md1
#  CDS #2 [2 exons; -]:replication associated protein:start2
#  CDS #2 [2 exons; -]:replication associated protein:stop2
#  CDS #2 [2 exons; -]:replication associated protein:fid2
#  CDS #2 [2 exons; -]:replication associated protein:md2
#  CDS #2 [2 exons; -]:replication associated protein:length
#  CDS #2 [2 exons; -]:replication associated protein:SS3
#  CDS #2 [2 exons; -]:replication associated protein:stp
#  CDS #2 [2 exons; -]:replication associated protein:PF
#  GenBank annotation:cds
#  GenBank annotation:exons
#  GenBank annotation:match
#  overlaps?
#  result

  # first, initialize the @{$out_header_exp_AR} with the first two lines:
  push(@{$out_header_exp_AR}, "#\n");
  push(@{$out_header_exp_AR}, "# Explanations of row headings on each page:\n");
  push(@{$out_header_exp_AR}, "#\n");

  # column/row #2: 'idx'
  $tok1 = sprintf("%-4s  ", "");
  $tok2 = $tok1;
  $tok3 = $tok1;
  $tok4 = sprintf("%-4s  ", " idx");
  $tok5 = sprintf("%-4s  ", "----");
  output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok4, undef, undef);
  output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok4, undef, undef, "index of genome in list", $FH_HR);

  # column/row #2: 'accession'
  $tok1 = sprintf("%-19s  ", "");
  $tok2 = $tok1;
  $tok3 = $tok1;
  $tok4 = sprintf("%-19s  ", " accession");
  $tok5 = sprintf("%-19s  ", "-------------------");
  output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); 
  output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok4, undef, undef, "GenBank accession for genomic sequence", $FH_HR);

  # column/row #3: 'totlen'
  $tok1 = sprintf("%-6s", "");
  $tok2 = $tok1;
  $tok3 = $tok1;
  $tok4 = sprintf("%-6s", "totlen");
  $tok5 = sprintf("%-6s", "------");
  output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok4, undef, undef); 
  output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok4, undef, undef, "total length (nt) for accession", $FH_HR);
  output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR); # adds a blank line

  if((opt_IsUsed("--origin", $opt_HHR)) || 
     (opt_IsUsed("--aorgstart", $opt_HHR))) { 
    # column/row #4: 'origin sequence:#'
    $tok1 = sprintf("  %22s", "");
    $tok2 = sprintf("  %22s", "   origin sequence");
    $tok3 = sprintf("  %22s", "----------------------");
    $tok4 = sprintf(" %2s", " #");
    $tok5 = sprintf(" %2s", "--");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "number of occurrences of origin sequence (input with -oseq) in genome", $FH_HR);

    # column/row #5: 'origin sequence:start'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %5s", "start");
    $tok5 = sprintf(" %5s", "-----");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "start position of lone occurrence of origin sequence (if only 1 exists)", $FH_HR);

    # column/row #6: 'origin sequence:stop'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %5s", "stop");
    $tok5 = sprintf(" %5s", "-----");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "stop  position of lone occurrence of origin sequence (if only 1 exists)", $FH_HR);

    # column/row #7: 'origin sequence:1stpos'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %5s", "1stps");
    $tok5 = sprintf(" %5s", "-----");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); 
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "what should be first position of genome, based on origin prediction", $FH_HR);

    # column/row #8: 'origin sequence:offst'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %5s", "offst");
    $tok5 = sprintf(" %5s", "-----");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); 
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "predicted offset of genome, number of nucleotides to shift start (>0: clockwise; <0: counterclockwise)", $FH_HR);

    # column/row #9: 'origin sequence:PF'
    # tok1, tok2, tok3 do not change
    $tok4 = sprintf(" %2s", "PF");
    $tok5 = sprintf(" %2s", "--");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "'P' (for PASS) if there is exactly 1 occurrence of the offset, else 'F' for FAIL", $FH_HR);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR); # adds a blank line
    push(@pf_text_A, "P/F character $pf_idx pertains to the origin sequence test");
    $pf_idx++;
  } # end of 'if(opt_IsUsed("--origin", $opt_HHR))'

  # create columns for 5'UTR, if $do_matpept:
  if($do_matpept) { 
    $width = 6 + 1 + 6 + 1 + 6; #20
    $tok1 = sprintf("  %*s", $width, "");
    $tok2 = sprintf("         %*s", $width, "5' UTR");
    $tok3 = sprintf("  %*s", $width, getMonocharacterString($width, "-", $FH_HR));
    $tok4 = sprintf("  %6s", "start");
    $tok5 = sprintf("  ------");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); 
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "start position of 5' UTR (inferred from other predictions, \"?\" if first mat_peptide is not predicted)", $FH_HR);

    $tok4 = sprintf(" %6s", "stop");
    $tok5 = sprintf(" ------");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "stop  position of 5' UTR (inferred from other predictions, \"?\" if first mat_peptide is not predicted)", $FH_HR);

    $tok4 = sprintf(" %6s", "length");
    $tok5 = sprintf(" ------");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "length of 5' UTR (inferred from other predictions, \"?\" if first mat_peptide is not predicted)", $FH_HR);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR); # adds a blank line
  }

  # create rows for each feature 
  $width = 0;
  my $do_multi_explanation = 1;
  my $do_model_explanation = 1;
  my $nmultifeature = numNonNumericValueInArray($ftr_info_HAR->{"annot_type"}, "multifeature", $FH_HR);
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $ftr_out_short   = $ftr_info_HAR->{"out_short"}[$ftr_idx];
    my $ftr_out_product = $ftr_info_HAR->{"out_product"}[$ftr_idx];
    #####################################################################################
    # block that handles multi-mat_peptide CDS (cds-mp, multifeature) feature annotations
    #####################################################################################
    if(($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") &&
       ($ftr_info_HAR->{"type"}[$ftr_idx]       eq "cds-mp")) { 
      if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") { 
        $width = 6 + 1 + 6 + 1 + 6; #20
        $tok1     = sprintf("  %*s", $width, $ftr_out_short . getMonocharacterString(int($width-length($ftr_out_short)/2), " ", $FH_HR));
        $exp_tok1 = "CDS(MP) #<i>";
        $tok2 = sprintf("  %s", $ftr_out_product); 
        $tok3 = sprintf("  %s", getMonocharacterString($width, "-", $FH_HR));
        $tok4 = sprintf("  %8s", sprintf("%s", "start"));
        $tok5 = sprintf("  %8s", "--------");

        output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
        if($do_multi_explanation) { 
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "start position of CDS(MP) #<i> (inferred from mat_peptides that comprise it,", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef, undef, "\"?\" if first mat_peptide for this CDS is not predicted)", $FH_HR);
        }
        
        $tok4 = sprintf(" %8s", "stop");
        $tok5 = sprintf(" %8s", "------");
        output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
        if($do_multi_explanation) { 
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "stop  position of CDS(MP) #<i> (inferred from mat_peptides that comprise it,", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef, undef, "\"?\" if final mat_peptide for this CDS is not predicted)", $FH_HR);
        }
        
        $tok4 = sprintf(" %6s", "length");
        $tok5 = sprintf(" %6s", "------");
        output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
        if($do_multi_explanation) { 
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "length of CDS(MP) #<i> (\"?\" if any of the mat_peptides that comprise this CDS are not predicted)", $FH_HR);
        }

        $tok4 = sprintf(" %6s", "startc");
        $tok5 = sprintf(" %6s", "------");
        output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
        if($do_multi_explanation) { 
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "start codon of CDS(MP) #<i> (\"?\" if first mat_peptide for this CDS is not predicted)", $FH_HR);
        }
        
        $tok4 = sprintf(" %6s", "stopc");
        $tok5 = sprintf(" %6s", "------");
        output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
        if($do_multi_explanation) { 
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "stop codon of CDS(MP) #<i> (\"?\" if final mat_peptide for this CDS is not predicted)", $FH_HR);
        }
        
        if($do_ss3) { 
          $tok4 = sprintf(" %3s", "ss3");
          $tok5 = sprintf(" %3s", "---");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
          if($do_multi_explanation) { 
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if predicted CDS has a valid start codon, stop codon and is a multiple of 3", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "first  character: '.' if predicted CDS has a valid start codon, '!' if not,", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "                  and '?' if first mat_peptide for this CDS is not predicted", $FH_HR);          
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "second character: '.' if predicted CDS has a valid stop  codon, '!' if not,", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "                      and '?' if final mat_peptide for this CDS is not predicted", $FH_HR);      
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "third  character: '.' if predicted CDS has a length which is a multiple of three, '!' if it is not a", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "                  multiple of 3, and '?' if any of the mat_peptides that comprise it are not predicted.", $FH_HR);
          }
        }
        
        $tok4 = sprintf(" %6s", "PF");
        $tok5 = sprintf(" %6s", "---");
        output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
        if($do_multi_explanation) { 
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if this CDS PASSED ('P') or FAILED ('F')", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  a CDS sequence PASSES ('P') if and only if all of the following", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  conditions are met (else it FAILS):", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (1) it has a valid start codon at beginning of its first mat_peptide", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (2) it has a valid stop  codon immediately after the end of its final mat_peptide", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (3) its length is a multiple of 3", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (4) all of the mat_peptides that comprise it are adjacent", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR);
          if($nmultifeature == 1) { 
            push(@pf_text_A, sprintf("P/F character %d pertains to the lone CDS.", $pf_idx));
          }
          else {
            push(@pf_text_A, sprintf("P/F characters %d to %d pertain to each of the %d CDS, in order.", $pf_idx, $pf_idx + $nmultifeature-1, $nmultifeature));
          }
          $pf_idx += $nmultifeature;
        }
        $do_multi_explanation = 0;
      }
    } # end of 'if' entered if feature is a multifeature cds-mp feature
  
    #############################################
    # block that handles 'annot_type' eq "model"
    # features, these are cds-notmp and mp types
    #############################################
    else { 
      for(my $mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx]; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
        $width += 18;
        my $mdl_exon_idx    = $mdl_info_HAR->{"map_exon"}[$mdl_idx];
        my $is_matpept = ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp") ? 1 : 0;
        if($do_fid)  { $width += 6;  }
        if($do_mdlb) { $width += 4;  }
        if($do_olap) { $width += 11; }
        if($is_matpept)  { $width += 11; }
        if($mdl_info_HAR->{"is_final"}[$mdl_idx]) { 
          $width += 9;
          if($do_ss3)  { $width += 4; }
          if($do_stop) { $width += 4; }
          $tok1     = sprintf("  %*s", $width, $ftr_out_short . getMonocharacterString(int(($width-length($ftr_out_short))/2), " ", $FH_HR));
          $exp_tok1 = "";
          if($do_matpept && $do_cds_notmp) { 
            $exp_tok1 = "{MP,CDS} #<i>";
          }
          elsif($do_matpept && (! $do_cds_notmp)) { 
            $exp_tok1 = "MP #<i>";
          }
          else { # $do_matpept is false
            $exp_tok1 =  "CDS #<i>";
          }
          $tok2 = sprintf("  %s", $ftr_out_product); 
          $tok3 = sprintf("  %s", getMonocharacterString($width, "-", $FH_HR));
          $tok4 = sprintf("  %8s", sprintf("%s%s", "start", $mdl_exon_idx+1));
          $exp_tok4 = "start<j>";
          $tok5 = sprintf("  %8s", "--------");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
          $width = 0; # reset width, this is impt
        }
        else { # not the final exon, still need start coordinate
          $width += 1;
          $tok1 = sprintf("    %s", $ftr_out_short);   # used only for output_tbl_get_headings_helper
          $tok2 = sprintf("    %s", $ftr_out_product); # used only for output_tbl_get_headings_helper
          $tok4 = sprintf("  %8s", sprintf("%s%s", "start", $mdl_exon_idx+1));
          $exp_tok4 = "start<j>";
          $tok5 = sprintf("  %8s", "--------");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
        }
        my $exp_substr = "";
        if($do_matpept && $do_cds_notmp) { 
          $exp_substr = "coding sequence part (or exon) <j> of mat_peptide (or CDS)";
        }
        elsif($do_matpept && (! $do_cds_notmp)) { 
          $exp_substr = "coding sequence part <j> of mat_peptide";
        }
        else { # $do_matpept is false
          $exp_substr = "exon <j> of CDS";
        }
        if($do_model_explanation) { 
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "start position of $exp_substr (\"NP\" if no prediction)", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "enclosed in brackets \"\[e\]\" if start/stop different from all exon start/stops in existing GenBank annotation", $FH_HR);
        }
        
        # stop, fid, and md rows take place for all exons
        # only token 4 changes
        $tok4 = sprintf(" %8s", sprintf("%s%s", "stop", $mdl_exon_idx+1));
        $exp_tok4 = "stop<j>";
        $tok5 = sprintf(" %8s", "--------");
        output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
        if($do_model_explanation) { 
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "stop  position of $exp_substr (\"NP\" if no prediction)", $FH_HR);
          output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "enclosed in brackets \"\[e\]\" if start/stop different from all exon start/stops in existing GenBank annotation", $FH_HR);
        }
        
        if($do_fid) { 
          $tok4 = sprintf(" %5s", sprintf("%s%s", "fid", $mdl_exon_idx+1));
          $exp_tok4 = "fid<j>";
          $tok5 = sprintf(" %5s", "-----");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
          if($do_model_explanation) { 
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "fractional identity between $exp_substr and reference genome (\"NP\" if no prediction)", $FH_HR);
          }
        }

        $exp_substr = $is_matpept ? "mat_peptide coding sequence" : "exon coding sequence";
        if($do_mdlb) { 
          $tok4 = sprintf(" %3s", sprintf("%s%s", "md", $mdl_exon_idx+1));
          $exp_tok4 = "md<j>";
          $tok5 = sprintf(" %3s", "---");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
          if($do_model_explanation) { 
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "annotation indicating if alignment to reference extends to 5' and 3' end of reference $exp_substr.", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "first character pertains to 5' end and second character pertains to 3' end.", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "possible values for each of the two characters:", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "  \".\":   alignment extends to boundary of reference", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "  \"<d>\": alignment truncates <d> nucleotides short of boundary of reference (1 <= <d> <= 9)", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "  \"+\":   alignment truncates >= 10 nucleotides short of boundary of reference", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "  \"t\":   position has been corrected based on predicted, truncated protein sequence", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "           3' position: stop coordinate has been adjusted to first in-frame stop (5' of predicted stop)", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "  \"e\":   position has been corrected based on predicted, extended protein sequence", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "           3' position: stop coordinate has been adjusted to first in-frame stop (3' of predicted stop)", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "  \"-\":   exon/segment is not predicted due to stop codon in earlier exon/segment", $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef,     undef,     undef, "  \"NP\":  (spanning both characters) no prediction", $FH_HR);
          }
        }
        
        if($do_olap) { 
          $tok4 = sprintf(" %10s", sprintf("%s%s", "overlaps", $mdl_exon_idx+1));
          $exp_tok4 = "overlaps<j>";
          $tok5 = sprintf(" %10s", "----------");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
          if($do_model_explanation) { 
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, sprintf("'P' or 'F' followed by list of %s this %s overlaps with", ($is_matpept) ? "mat_peptide" : "exon", ($is_matpept) ? "mat_peptide" : "exon"), $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "first letter is 'P' if agrees exactly with reference, else 'F'", $FH_HR); # adds a second line to explanation
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "\"NP\" if no prediction", $FH_HR);
            $need_to_define_H{"overlap"} = 1;
          }
        }

        if($is_matpept) { 
          $tok4 = sprintf(" %10s", sprintf("%s%s", "adjcnces", $mdl_exon_idx+1));
          $exp_tok4 = "adjcnces<j>";
          $tok5 = sprintf(" %10s", "----------");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
          if($do_model_explanation) { 
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, sprintf("'P' or 'F' followed by list of %s this %s is adjacent with", ($is_matpept) ? "mat_peptide" : "exon", ($is_matpept) ? "mat_peptide" : "exon"), $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "first letter is 'P' if agrees exactly with reference, else 'F'", $FH_HR); # adds a second line to explanation
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "\"NP\" if no prediction", $FH_HR);      
          }
          $need_to_define_H{"adjacent"} = 1;
        }

        $exp_substr = $is_matpept ? "mat_peptide coding sequence" : "CDS";
        if($mdl_info_HAR->{"is_final"}[$mdl_idx]) { 
          $tok4 = sprintf(" %6s", "length");
          $tok5 = sprintf(" %6s", "------");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); 
          if($do_model_explanation) { 
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, sprintf("length of $exp_substr #<i> (all %s summed)", $is_matpept ? "segments" : "exons"), $FH_HR);
          }      

          if((! $is_matpept) && ($do_ss3)) { # skip this in matpept mode, we don't check start/stop of mat_peptides, only CDS, later
            $tok4 = sprintf(" %3s", "ss3");
            $tok5 = sprintf(" %3s", "---");
            output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
            if($do_model_explanation) { 
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if predicted CDS has a valid start codon, stop codon and is a multiple of 3", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "first  character: '.' if predicted CDS has a valid start codon, else '!'", $FH_HR);          
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "second character: '.' if predicted CDS has a valid stop  codon, else '!'", $FH_HR);      
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "third  character: '.' if predicted CDS has a length which is a multiple of three, else '!'", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "\"NP\" if no prediction", $FH_HR);
            }
          }
          if((! $is_matpept) && ($do_stop)) { # skip this in matpept mode, we only check stop of final mat_peptide, later
            $tok4 = sprintf(" %3s", "stp");
            $tok5 = sprintf(" %3s", "---");
            output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
            if($do_model_explanation) { 
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, "the predicted stop codon for this CDS (\"NP\" if no prediction)", $FH_HR);
            }
          }
          
          $tok4 = sprintf(" %2s", "PF");
          $tok5 = sprintf(" %2s", "--");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
          if($do_model_explanation) { 
            if($is_matpept) { 
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if this mat_peptide PASSED ('P') or FAILED ('F')", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  a mat_peptide coding sequence PASSES ('P') if and only if", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  conditions are met (else it FAILS):", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (1) it has a valid start codon or homologous reference mat_peptide", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      does not", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (2) it has a valid stop  codon immediately after its predicted", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      end or reference mat_peptide does not", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (3) its length is a multiple of 3", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (4) has a pairwise alignment to the homologous reference met_peptide that", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      extends to the 5' and 3' boundary of the reference annotation", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (5) overlaps with exact same set of other mat_peptides as the homologous", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      reference mat_peptide", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (6) is adjacent to the exact same set of other mat_peptides as the", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      homologous reference mat_peptide", $FH_HR);
              push(@pf_text_A, sprintf("P/F characters %d to %d pertain to each of the %d mature peptides, in order.", $pf_idx, $pf_idx + $nftr-1, $nftr));
              $pf_idx += $nftr;
              $need_to_define_H{"adjacent"} = 1;
            }
            else { 
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if this CDS PASSED ('P') or FAILED ('F')", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  a CDS sequence PASSES ('P') if and only if all of the following", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  conditions are met (else it FAILS):", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (1) it has a valid start codon at beginning of its first exon", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (2) it has a valid stop  codon at end of its final exon", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (3) its length is a multiple of 3", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (4) all of its exons have a pairwise alignment to the homologous", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      reference exon that extends to the 5' and 3' boundary of the", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      reference annotation.", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (5) all of its exons overlap with exact same set of other exons as the", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      homologous reference CDS", $FH_HR);
              push(@pf_text_A, sprintf("P/F characters %d to %d pertain to each of the %d CDS, in order.", $pf_idx, $pf_idx + $nftr-1, $nftr));
              $pf_idx += $nftr;
              $need_to_define_H{"overlap"} = 1;
            }
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR);
          }
        } # end of 'if($mdl_is_final_AR->[$mdl_idx])'
        $do_model_explanation = 0; # once we see the final exon of the first CDS, we don't need to print CDS explanations anymore
      } # end of 'for(my $mdl_idx)'
    } # end of 'else' entered if we're not a multifeature cds-mp feature
  } # end of 'for(my $ftr_idx)'
  
  # create columns for 3'UTR, if $do_matpept:
  if($do_matpept) { 
    $width = 6 + 1 + 6 + 1 + 6; #20
    $tok1 = sprintf("  %*s", $width, "");
    $tok2 = sprintf("         %*s", $width, "3' UTR");
    $tok3 = sprintf("  %*s", $width, getMonocharacterString($width, "-", $FH_HR));
    $tok4 = sprintf("  %6s", "start");
    $tok5 = sprintf("  ------");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef); 
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "start position of 3' UTR (inferred from other predictions, \"?\" if final mat_peptide is not predicted)", $FH_HR);

    $tok4 = sprintf(" %6s", "stop");
    $tok5 = sprintf(" ------");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "stop  position of 3' UTR (inferred from other predictions, \"?\" if final mat_peptide is not predicted)", $FH_HR);

    $tok4 = sprintf(" %6s", "length");
    $tok5 = sprintf(" ------");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "length of 3' UTR (inferred from other predictions, \"?\" if final mat_peptide is not predicted)", $FH_HR);

    output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR);
  }

  # "totlen"
  $tok1 = sprintf("  %6s", "");
  $tok2 = sprintf("  %6s", "");
  $tok3 = sprintf("  %6s", "");
  $tok4 = sprintf("  %6s", "totlen");
  $tok5 = sprintf("  %6s", "------");
  output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok4, undef, undef);
  output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok4, undef, undef, "total length (nt) for accession (repeated for convenience)", $FH_HR);

  if($do_totfid) { 
    # "totfid"
    $tok1 = sprintf("  %5s", "");
    $tok2 = sprintf("  %5s", "");
    $tok3 = sprintf("  %5s", "");
    $tok4 = sprintf("  %5s", "totfid");
    $tok5 = sprintf("  %5s", "-----");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok4, undef, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok4, undef, undef, "fractional identity of all concatenated pairwise nucleotide alignments for this accession", $FH_HR);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR);
  }

  # existing GenBank annotation
  if($do_exist) { 
    $tok1 = sprintf("  %19s", "");
    $tok2 = sprintf("  %19s", "GenBank annotation");
    $tok3 = sprintf("  %19s", "-------------------");
    $tok4 = sprintf("  %5s", ($do_matpept) ? "mp" : "cds");
    $tok5 = sprintf("  %5s", "-----");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, sprintf("number of %s in the existing GenBank annotation for this accession", ($do_matpept) ? "mat_peptides" : "CDS"), $FH_HR);
    
    $tok4 = sprintf("  %5s", "exons");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "total number of exons in the existing GenBank annotation for this accession", $FH_HR);
    
    $tok4 = sprintf("  %5s", "match");
    output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok2, $tok4, undef);
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok2, $tok4, undef, "number of exons in the existing GenBank annotation for which existing and predicted annotation agree exactly", $FH_HR);
  }
  output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR);

  # result
  $tok1 = sprintf("  %*s",  $width_result, "");
  $tok2 = sprintf("  %*s",  $width_result, "");
  $tok3 = sprintf("  %*s",  $width_result, "");
  $tok4 = sprintf("  %-*s", $width_result, "result");
  $tok5 = sprintf("  %-*s", $width_result, getMonocharacterString($width_result, "-", $FH_HR));
  output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok4, undef, undef);

  output_tbl_get_headings_explanation_helper($out_header_exp_AR, $tok4, undef, undef, "\"PASS\" or \"FAIL\". \"PASS\" if and only if all tests for this accession PASSED ('P')", $FH_HR);
  output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "as indicated in the \"PF\" rows.", $FH_HR);
  output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, sprintf("After that is a string of %d characters, these are the individual P/F results in order.", $pf_idx-1), $FH_HR);
  foreach my $pf_text_str (@pf_text_A) { 
    output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, $pf_text_str, $FH_HR);
  }
  output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "See explanations of the individual P/F values above.", $FH_HR);

  if((defined $need_to_define_H{"overlap"}) || (defined $need_to_define_H{"adjacent"})) {
     output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR);
     push(@{$out_header_exp_AR}, "# Definitions of non-obvious terms above:\n");
     output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR);
     if(defined $need_to_define_H{"overlap"}) { 
       push(@{$out_header_exp_AR}, "# overlap:  two features i and j overlap if they are on both on the same strand and overlap by >= 1 nt.\n");
     }
     if(defined $need_to_define_H{"adjacent"}) { 
       push(@{$out_header_exp_AR}, "# adjacent: two features i and j are adjacent if they are on the same strand and\n");
       push(@{$out_header_exp_AR}, "#           start_i < start_j and stop_i+1 == start_j.\n");
       push(@{$out_header_exp_AR}, "#           (Note that on the positive strand start_i <= stop_i for all i,\n");
       push(@{$out_header_exp_AR}, "#           (and that  on the negative strand start_i >= stop_i for all i)\n");
     }
  }
  return;
}

#################################################################
# Subroutine: output_tbl_get_headings_helper()
# Incept:     EPN, Fri Mar 11 04:50:55 2016
#
# Purpose:   Helper function for output_tbl_get_headings() when
#            used in sequences-as-columns modes. Given up to 3 tokens,
#            add them to the appropriate place in @{$out_col_header_AAR}.
#
# Arguments:
#   $out_col_header_AAR: ref to output column header 2D array
#   $div_char:           divider character to put between tokens
#   $tok1:               token 1, can be undef
#   $tok2:               token 2, can be undef
#   $tok3:               token 3, can be undef
#             
# Returns:  void
# 
# Dies:     Never.
#
#################################################################
sub output_tbl_get_headings_helper { 
  my $sub_name = "output_tbl_get_headings_helper";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_row_header_AR, $div_char, $tok1, $tok2, $tok3) = @_;
  
  # remove whitespace at beginning and end of tokens
  if(defined $tok1) { $tok1 =~ s/^\s+//; $tok1 =~ s/\s+$//; }
  if(defined $tok2) { $tok2 =~ s/^\s+//; $tok2 =~ s/\s+$//; }
  if(defined $tok3) { $tok3 =~ s/^\s+//; $tok3 =~ s/\s+$//; }

  my $toadd = "";
  if(defined $tok1) { 
    $toadd = $tok1; 
  }
  if(defined $tok2) { 
    if($toadd ne "") { $toadd .= $div_char; }
    $toadd .= $tok2; 
  }
  if(defined $tok3) { 
    if($toadd ne "") { $toadd .= $div_char; }
    $toadd .= $tok3; 
  }

  push(@{$out_row_header_AR}, $toadd); 

  return;
}

#################################################################
# Subroutine: output_tbl_get_headings_explanation_helper()
# Incept:     EPN, Thu Mar 10 21:05:36 2016
#
# Purpose:   Helper function for output_tbl_get_headings() for 
#            adding explanatory text to the @{$out_header_exp_AR} array.
#             Given up to 3 tokens that define the header, and one that is the 
#             explanatory text. Can be used in 3 modes:
#
#             Mode 1: at least one of $tok1, $tok2, $tok3 is defined
#                     and $desc is defined 
#                     In this mode, determine header by concatenating all of
#                     $tok1, $tok2, and $tok3 that are defined and use
#                     $desc as the description.
#
#             Mode 2: none of $tok1, $tok2, $tok3 is defined and $desc is
#                     defined.
#                     In this mode, header is blank, and use $desc as 
#                     the description.
#
#             Mode 3: none of $tok1, $tok2, $tok3 is defined and $desc is
#                     not defined either
#                     In this mode, add a blank line to @{$out_row_header_AR}.
#
# Arguments:
#   $out_header_exp_AR:  ref to output column header 2D array
#   $tok1:               token 1, can be undef
#   $tok2:               token 2, can be undef
#   $tok3:               token 3, can be undef
#   $desc:               description text, can be undef
#   $FH_HR:              REF to hash of file handles
#             
# Returns:  void
# 
# Dies:     if desc is not defined but a header token is
#
#################################################################
sub output_tbl_get_headings_explanation_helper { 
  my $sub_name = "output_tbl_get_headings_explanation_helper";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_header_exp_AR, $tok1, $tok2, $tok3, $desc, $FH_HR) = @_;

  my $width = 35;
  # remove whitespace at beginning and end of tokens
  if(defined $tok1) { $tok1 =~ s/^\s+//; $tok1 =~ s/\s+$//; }
  if(defined $tok2) { $tok2 =~ s/^\s+//; $tok2 =~ s/\s+$//; }
  if(defined $tok3) { $tok3 =~ s/^\s+//; $tok3 =~ s/\s+$//; }

  # we don't allow whitespace between headers, the idea is that 
  # each header should be a single white-space delimited token, 
  # so we can more easily manipulate the output file
  my $header = "";
  if(defined $tok1) { 
    $header .= $tok1;
  }
  if(defined $tok2) { 
    if($header ne "") { $header .= ":"; }
    $header .= $tok2;
  }
  if(defined $tok3) { 
    if($header ne "") { $header .= ":"; }
    $header .= $tok3;
  }
  if($header ne "") { 
    $header = "\"" . $header . "\":";
  }

  if(defined $desc) { 
    push(@{$out_header_exp_AR}, sprintf("# %-*s %s\n", $width, $header, $desc)); 
  }
  else { 
    if($header ne "") { 
      DNAORG_FAIL("ERROR in $sub_name, desc is not defined but one of the header tokens is", 1, $FH_HR);
    }
    push(@{$out_header_exp_AR}, "#\n");
  }

  return;
}

#################################################################
# Subroutine: output_tbl_all_sequences()
# Incept:     EPN, Sun Mar 13 21:13:58 2016
#
# Purpose:   Output the tabular annotation for all sequences.
#
# Arguments:
#  $mdl_info_HAR:     REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:     REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:     REF to hash of arrays with information on the sequences, PRE-FILLED
#  $mdl_results_AAHR: REF to model results AAH, PRE-FILLED
#  $ftr_results_AAHR: REF to feature results AAH, PRE-FILLED
#  $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:   REF to the 2D hash of output file information
#             
# Returns:  Number of accessions with >= 1 failure.
# 
# Dies:     never
#
#################################################################
sub output_tbl_all_sequences { 
  my $sub_name = "output_tbl_all_sequences";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $ftr_results_AAHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience
  my $tblsum_FH = $FH_HR->{"tblsum"};
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences

  # data structures necessary for storing output prior to actually printing
  my @ref_out_A   = (); # array of output fields for the reference accession
  my @page_accn_A = (); # [0..$cur_pagesize]2D array, each element is an array of output tokens for one accession
  my @page_out_AA = (); # 2D array, [0..$a..$cur_pagesize-1][0..(ntoks-1)] first dimension is of size $cur_pagesize, 
                        # each element is an array of output tokens for accession $page_accn_A[$a]
  my $cur_pagesize = 0; # current number of accessions we have info for in page_out_AA (size of first dimension in page_out_AA)
                        # when this hits $nseqcol, we dump the output
  my $npages = 0;       # number of pages output
  # analagous data structures for 'failure-only' version of the file
  my @fail_page_accn_A = (); # [0..$cur_fail_pagesize]2D array, each element is an array of output tokens for one accession
  my @fail_page_out_AA = (); # 2D array, [0..$a..$cur_fail_pagesize-1][0..(ntoks-1)] first dimension is of size $cur_pagesize, 
                             # each element is an array of output tokens for accession $fail_page_accn_A[$a]
  my $cur_fail_pagesize = 0; # current number of accessions we have info for in fail_page_out_AA (size of first dimension in fail_page_out_AA)
                             # when this hits $nseqcol, we dump the output
  my $nfail_pages = 0;       # number of fail pages output
  # analagous data structures for 'error-only' version of the file
  my @err_page_accn_A = (); # [0..$cur_err_pagesize]2D array, each element is an array of output tokens for one accession
  my @err_page_out_AA = (); # 2D array, [0..$a..$cur_err_pagesize-1][0..(ntoks-1)] first dimension is of size $cur_pagesize, 
                            # each element is an array of output tokens for accession $err_page_accn_A[$a]
  my $cur_err_pagesize = 0; # current number of accessions we have info for in err_page_out_AA (size of first dimension in err_page_out_AA)
                            # when this hits $nseqcol, we dump the output
  my $nerr_pages = 0;       # number of error pages output

  # variables related to optional output 
  my $do_fid      = (opt_Get("--doalign", $opt_HHR)) ? 1 : 0; # '1' to skip fid output
  my $do_totfid   = (opt_Get("--doalign", $opt_HHR)) ? 1 : 0; # '1' to skip fid output
  my $do_ss3      = 1; # '1' to do ss3 output, '0' to skip it
  my $do_stop     = 1; # '1' to do stop output, '0' to skip it
  my $do_mdlb     = 1; # '1' to do model boundary output, '0' to skip it
  my $do_olap     = 1; # '1' to do overlap output, '0' to skip it 
  my $do_exist    = (opt_Get("--infasta", $opt_HHR) || opt_Get("--tblnocomp", $opt_HHR)) ? 0 : 1; # '1' to do comparison to existing GenBank annotation, '0' to skip it
  my $do_matpept  = (numNonNumericValueInArray($ftr_info_HAR->{"type"}, "mp", $FH_HR) > 0) ? 1 : 0;

  my $nseqcol     = 5; # number of sequences we print per page
  my $do_tblfirst = (opt_Get("--tblfirst", $opt_HHR)) ? 1 : 0; 
  my $act_nseqcol = $do_tblfirst ? $nseqcol-1 : $nseqcol;

  # the possible values for the ss3 output (start/stop/multiple-of-3)
  my $ss3_yes_char    = "."; 
  my $ss3_unsure_char = "?";
  my $ss3_no_char     = "!";

  # miscellaneous variables
  my $at_least_one_fail = undef; # set to 1 if we see a failure, separately for each sequence
  my $pass_fail_char    = undef; # for each possible pass/fail, the 'P' or 'F'
  my $pass_fail_str     = undef; # string of pass/fail characters
  my $tot_nfail         = 0;     # total number of accessions with at least 1 failure

  my $origin_offset = undef;
  if(opt_IsUsed("--origin", $opt_HHR)) { 
    $origin_offset = validate_origin_seq(opt_Get("--origin", $opt_HHR));
  }

  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
    my $seq_len   = $seq_info_HAR->{"seq_len"}[$seq_idx];
    my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
    my $accn_len  = $seq_info_HAR->{"accn_len"}[$seq_idx];
    my @cur_out_A = (); # array of current tokens to print
    my $ngenbank_match = 0; # number of matches with existing annotation
    my $pass_fail_str  = ""; 

    # Create the initial portion of the output line, the accession and length
    push(@cur_out_A, sprintf("%-5d  ", ($seq_idx+1)));
    push(@cur_out_A, sprintf("%-19s  ", $accn_name)); 
    push(@cur_out_A, sprintf("%6d ", $accn_len));

    my ($oseq_ct, $oseq_start, $oseq_stop, $oseq_firstpos, $oseq_offset, $oseq_passfail);
    if((opt_IsUsed("--origin",    $opt_HHR)) || 
       (opt_IsUsed("--aorgmodel", $opt_HHR))) { 

      if(opt_IsUsed("--origin",    $opt_HHR)) { 
        ($oseq_ct, $oseq_start, $oseq_stop, $oseq_offset, $oseq_passfail) = get_origin_output_for_sequence($seq_info_HAR, $seq_idx, $origin_offset, $FH_HR);
        $oseq_firstpos = "?";
      }
      if(opt_IsUsed("--aorgmodel", $opt_HHR)) { 
        ($oseq_ct, $oseq_start, $oseq_stop, $oseq_firstpos, $oseq_offset, $oseq_passfail) = aorg_get_origin_output_for_sequence($seq_info_HAR, $seq_idx, $FH_HR);
      }

      push(@cur_out_A, sprintf("%2d ", $oseq_ct));
      push(@cur_out_A, sprintf("%5s ", $oseq_start));
      push(@cur_out_A, sprintf("%5s ", $oseq_stop));
      push(@cur_out_A, sprintf("%5s ", $oseq_firstpos));
      push(@cur_out_A, sprintf("%5s ", $oseq_offset));
      push(@cur_out_A, sprintf(" %s", $oseq_passfail));
      $pass_fail_str .= $oseq_passfail;
    }

    # we will eventually output the total fractional identity, the fractional identity we would
    # get if we concatenated all features aligned to the reference. 
    my $tot_nid        = 0.; # number of identical positions between this feature and the reference
    my $tot_id_min_len = 0.; # total number of possible identities, this is the sum for all features of the minimum 
                             # of the reference length for the feature and the annotated length for the feature

    # 5' UTR, if nec
    if($do_matpept) { 
      if(! exists $mdl_results_AAHR->[0][$seq_idx]{"p_start"}) { 
        push(@cur_out_A, sprintf("  %6s", "?")); # start
        push(@cur_out_A, sprintf(" %6s", "?"));  # stop
        push(@cur_out_A, sprintf(" %6s", "?"));  # length
      }
      else { # we know that $mdl_results_AAHR->[0][$seq_idx]{"p_start"} exists) { 
        # determine output start and output stop
        my ($cur_start, undef) = create_output_start_and_stop($mdl_results_AAHR->[0][$seq_idx]{"p_start"},
                                                              $mdl_results_AAHR->[0][$seq_idx]{"p_stop"},
                                                              $accn_len, $seq_len, $FH_HR);
        if($mdl_results_AAHR->[0][$seq_idx]{"p_strand"} eq "+") { 
          # positive strand, easy case
          if($cur_start == 1) { 
            push(@cur_out_A, sprintf("  %6d", 0)); # start 
            push(@cur_out_A, sprintf("  %6d", 0)); # stop 
            push(@cur_out_A, sprintf("  %6d", 0)); # length
          }
          else { # 1st matpept does not start at nt 1 (normal case)
            push(@cur_out_A, sprintf("  %6d", 1)); # start 
            push(@cur_out_A, sprintf("  %6d", $cur_start - 1)); # stop
            push(@cur_out_A, sprintf("  %6d", $cur_start - 1)); # length
          }
        }
        elsif($mdl_results_AAHR->[0][$seq_idx]{"p_strand"} eq "-") { 
          # negative strand, more complicated, slightly
          if($cur_start == $accn_len) { 
            push(@cur_out_A, sprintf("  %6d", 0)); # start 
            push(@cur_out_A, sprintf("  %6d", 0)); # stop 
            push(@cur_out_A, sprintf("  %6d", 0)); # length
          }
          else { # 1st feature does not start at nt $accn_len on negative strand
            push(@cur_out_A, sprintf("  %6d", $accn_len)); # start 
            push(@cur_out_A, sprintf("  %6d", $cur_start + 1)); # stop
            push(@cur_out_A, sprintf("  %6d", $accn_len - $cur_start)); # length
          }
        }
        else { # not + or - strand, weird...
          DNAORG_FAIL("ERROR in $sub_name, trying to compute 5' UTR for prediction that exists but is not + or - strand", 1, $FH_HR);
        }
      }
    }

    # go through each feature and collect output tokens and add to @cur_out_A
    my $start_codon_char   = ""; # set below if for models if $is_first
    my $stop_codon_char    = ""; # set below if for models if $is_final
    my $multiple_of_3_char = ""; # set below if for models if $is_final
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      #####################################################################################
      # block that handles multi-mat_peptide CDS (cds-mp, multifeature) feature annotations
      #####################################################################################
      if(($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") &&
         ($ftr_info_HAR->{"type"}[$ftr_idx]       eq "cds-mp")) { 
        # start, stop, length, start_codon, start_codon character
          my $ftr_results_HR = \%{$ftr_results_AAHR->[$ftr_idx][$seq_idx]}; # for convenience
          push(@cur_out_A, sprintf("  %8s", $ftr_results_HR->{"out_start"}));
        push(@cur_out_A, sprintf("  %8s", $ftr_results_HR->{"out_stop"})); 
        push(@cur_out_A, sprintf("  %6s", $ftr_results_HR->{"out_len"}));  
        push(@cur_out_A, sprintf("  %6s", $ftr_results_HR->{"out_start_codon"}));
        push(@cur_out_A, sprintf("  %6s", $ftr_results_HR->{"out_stop_codon"}));
        
        my $start_codon_char   = undef;
        if(($ftr_results_HR->{"out_start"}) eq "?") { 
          $start_codon_char = "?";
        }
        elsif((exists $ftr_results_HR->{"str_err_flag"})) { 
          $start_codon_char = $ss3_no_char;
        }
        else { 
          $start_codon_char = $ss3_yes_char;
        }

        my $stop_codon_char    = undef;
        if(($ftr_results_HR->{"out_stop"}) eq "?") { 
          $stop_codon_char = "?";
        }
        elsif(! validateStopCodon($ftr_results_HR->{"out_stop_codon"})) { 
          $stop_codon_char = $ss3_no_char;
        }
        else { 
          $stop_codon_char = $ss3_yes_char;
        }

        my $multiple_of_3_char = undef;
        if(($ftr_results_HR->{"out_len"}) eq "?") { 
          $multiple_of_3_char = "?";
        }
        elsif(($ftr_results_HR->{"out_len"} % 3) != 0) { 
          $multiple_of_3_char = $ss3_no_char;
        }
        else { 
          $multiple_of_3_char = $ss3_yes_char;
        }
        # determine if this CDS passed or failed
        my $cds_pass_fail = "P";
        if(($start_codon_char ne $ss3_yes_char || $stop_codon_char ne $ss3_yes_char || $multiple_of_3_char ne $ss3_yes_char) ||
           (exists $ftr_results_HR->{"trc_err_flag"})) { 
          $cds_pass_fail = "F";
        }
        push(@cur_out_A, sprintf(" %s%s%s", $start_codon_char, $stop_codon_char, $multiple_of_3_char)); #start,stop,mult_of_3
        push(@cur_out_A, sprintf("  %3s", $cds_pass_fail)); # cds_pass_fail
        $pass_fail_str .= $cds_pass_fail;
      }         

      #############################################
      # block that handles 'annot_type' eq "model"
      # features, these are cds-notmp and mp types
      #############################################
      else { # not a multifeature cds-mp 
        for(my $mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx]; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
          my $is_first = $mdl_info_HAR->{"is_first"}[$mdl_idx]; # is this the first model for feature $ftr_idx?
          my $is_final = $mdl_info_HAR->{"is_final"}[$mdl_idx]; # is this the final model for feature $ftr_idx?
          my $is_matpept = ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp") ? 1 : 0;
          my $mdl_results_HR = \%{$mdl_results_AAHR->[$mdl_idx][$seq_idx]}; # for convenience
          my $ref_olp_str = $mdl_info_HAR->{"out_olp_str"}[$mdl_idx];
          my $ref_adj_str = combine_ajb_and_aja_strings($mdl_info_HAR->{"out_ajb_str"}[$mdl_idx], $mdl_info_HAR->{"out_aja_str"}[$mdl_idx]);
          
          if($is_first) { # reset variables
            $at_least_one_fail  = 0; 
            $start_codon_char   = "";
            $stop_codon_char    = "";
            $multiple_of_3_char = "";
          }

          if(exists $mdl_results_HR->{"p_start"}) { 
            # hit exists
            # check for special case when we had a trc in a previous segment for the same feature
            if((exists $mdl_results_HR->{"prv_trc_flag"}) && ($mdl_results_HR->{"prv_trc_flag"} == 1)) { # flag for a trc error in previous exon
              if($is_first) { 
                DNAORG_FAIL(sprintf("ERROR in $sub_name, found flag for trc error in earlier segment/exon but this is the first segment, $accn_name model: %s", $mdl_info_HAR->{"out_tiny"}[$mdl_idx]), 1, $FH_HR);
              }
              push(@cur_out_A, sprintf("  %8s ", "-"));                      # start 
              push(@cur_out_A, sprintf("%8s", "-"));                         # stop
              if($do_fid)      { push(@cur_out_A, sprintf(" %5s", "-")); }   # fid
              if($do_mdlb)     { push(@cur_out_A, "  " . "--"); }            # mdlb
              if($do_olap)     { push(@cur_out_A, sprintf(" %10s", "-")); }  # olap
              if($is_matpept)  { push(@cur_out_A, sprintf(" %10s", "-")); } # adj
            }
            else { 
              # not special case in which a trc error exists in earlier segment,
              # set annotation we do for all models (regardless of $is_first or $is_final values)
              my $genbank_match = $mdl_results_HR->{"genbank_mdl_annot_match"}; # 1 if existing GenBank annotation matches our annotation
              if($genbank_match) { $ngenbank_match++; }
              push(@cur_out_A, sprintf("  %8s ", (($genbank_match || (! $do_exist)) ? " " . $mdl_results_HR->{"out_start"} . " " : "[" . $mdl_results_HR->{"out_start"} . "]")));
              push(@cur_out_A, sprintf("%8s",    (($genbank_match || (! $do_exist)) ? " " . $mdl_results_HR->{"out_stop"}  . " " : "[" . $mdl_results_HR->{"out_stop"}  . "]")));
              if($do_fid) { push(@cur_out_A, sprintf(" %5.3f", $mdl_results_HR->{"fid2ref"})); } 
              if($do_totfid) { 
                # mdl_results_HR->{"fid2ref"} was computed as the
                # number of identities in the alignment of the
                # annotated feature to the reference feature divided
                # by the minimum of the length of the reference
                # feature and the length of the annotated feature, so
                # to determine total fid (fid we'd have for
                # concatenation of all aligned features to the
                # reference) we need to determine what the minimum
                # length is, as well as infer how many identities we
                # had.
                my $min_len = ($mdl_results_HR->{"len"} < $mdl_info_HAR->{"length"}[$mdl_idx]) ? $mdl_results_HR->{"len"} : $mdl_info_HAR->{"length"}[$mdl_idx];
                $tot_nid += ($mdl_results_HR->{"fid2ref"} * $min_len);
                $tot_id_min_len += $min_len;
              }
              if($do_mdlb) { 
                push(@cur_out_A, "  " . $mdl_results_HR->{"out_5boundary"} . $mdl_results_HR->{"out_3boundary"}); 
                if(($mdl_results_HR->{"out_5boundary"} ne ".") || 
                   (($mdl_results_HR->{"out_3boundary"} ne ".") && 
                    ($mdl_results_HR->{"out_3boundary"} ne "t") && 
                    ($mdl_results_HR->{"out_3boundary"} ne "e"))) { 
                  $at_least_one_fail = 1;
                }
              }
              if($do_olap) { 
                my $out_olp_str = ($mdl_results_HR->{"out_olp_str"} eq "") ? "NONE" : $mdl_results_HR->{"out_olp_str"};
                if($mdl_results_HR->{"out_olp_str"} ne $ref_olp_str) { 
                  $at_least_one_fail = 1; 
                  push(@cur_out_A, sprintf(" %10s", "F:" . $out_olp_str));
                }
                else { 
                  push(@cur_out_A, sprintf(" %10s", "P:" . $out_olp_str));
                }
              }
              
              if($is_matpept) { 
                my $adj_str = combine_ajb_and_aja_strings($mdl_results_HR->{"out_ajb_str"}, $mdl_results_HR->{"out_aja_str"});
                my $out_adj_str = ($adj_str eq "") ? "NONE" : $adj_str;
                if($adj_str ne $ref_adj_str) { 
                  $at_least_one_fail = 1; 
                  push(@cur_out_A, sprintf(" %10s", "F:" . $out_adj_str));
                }
                else { 
                  push(@cur_out_A, sprintf(" %10s", "P:" . $out_adj_str));
                }
              }
              if($is_first) { 
                if($mdl_results_HR->{"str_err_flag"}) { 
                  $at_least_one_fail = 1; 
                }
                $start_codon_char = $mdl_results_HR->{"str_err_flag"} ? $ss3_no_char : $ss3_yes_char;
              }
            } # end of 'else' entered if prv_trc_flag is *not* raised
            
            # now add annotation we only do for the final model of each feature, we do this
            # even if prv_trc_flag is raised
            if($is_final) { 
              if(! $is_matpept) { 
                if(validateStopCodon($mdl_results_HR->{"out_stop_codon"})) { 
                  $stop_codon_char = $ss3_yes_char;
                }
                else { 
                  $stop_codon_char   = $ss3_no_char;
                  $at_least_one_fail = 1;
                }
              } # end of 'if(! $is_matpept)'
              if(($mdl_results_HR->{"cumlen"} % 3) == 0) { 
                $multiple_of_3_char = $ss3_yes_char;
              }
              else { 
                $multiple_of_3_char = $ss3_yes_char;
                $at_least_one_fail = 1;
              }
              push(@cur_out_A, sprintf(" %6d", $mdl_results_HR->{"cumlen"})); 
              
              # add the ss3 (start/stop/multiple of 3 info) if we're not a mature peptide
              if(! $is_matpept) { 
                if($start_codon_char eq "") { die "ERROR $seq_idx ($seq_name) $mdl_idx start_codon_char is blank\n"; }
                if($stop_codon_char  eq "") { die "ERROR $seq_idx ($seq_name) $mdl_idx stop_codon_char is blank\n"; }
                if($multiple_of_3_char  eq "") { die "ERROR $seq_idx $mdl_idx multiple_of_3_char is blank\n"; }
                push(@cur_out_A,  sprintf(" %s%s%s", $start_codon_char, $stop_codon_char, $multiple_of_3_char));
                if($do_stop) { 
                  push(@cur_out_A, sprintf(" %3s", $mdl_results_HR->{"out_stop_codon"}));
                }
              }
              $pass_fail_char = ($at_least_one_fail) ? "F" : "P";
              push(@cur_out_A, sprintf(" %2s", $pass_fail_char));
              $pass_fail_str .= $pass_fail_char;
            } # end of 'if($is_final)'
          } # end of 'if(exists $mdl_results_HR->{"p_start"}'
          else { 
            # no prediction exists
            $at_least_one_fail = 1;
            push(@cur_out_A, sprintf("  %8s ", "NP")); # start position
            push(@cur_out_A, sprintf("%8s",  "NP"));   # stop position
            if($do_fid)  { push(@cur_out_A, sprintf(" %5s", "NP")); } # fid 
            if($do_mdlb) { push(@cur_out_A, "  NP"); } # model boundaries
            if($do_olap) { push(@cur_out_A, "  NP"); } # overlaps
            if($is_matpept)  { push(@cur_out_A, "  NP"); } # adjacencies
            if($is_final) { 
              push(@cur_out_A, sprintf(" %6s", "NP")); # length
              if((! $is_matpept) && ($do_ss3))  { push(@cur_out_A, "  NP"); } # ss3
              if((! $is_matpept) && ($do_stop)) { push(@cur_out_A, sprintf(" %3s", "NP")); } # stop
              $pass_fail_char = "F";
              push(@cur_out_A, sprintf(" %2s", $pass_fail_char));
              $pass_fail_str .= $pass_fail_char;
            }
            if($is_first) { # important to do this so final model has a valid start_codon_char
              $start_codon_char = "NP";
            }
          }
        } # end of 'for(my $mdl_idx'
      } # end of 'else' entered if feature is not a multifeature cds-mp
    } # end of 'for(my $ftr_idx'      

    # 3' UTR, if nec
    if($do_matpept) { 
      if(! exists $mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_start"}) { 
        push(@cur_out_A, sprintf("  %6s", "?")); # start
        push(@cur_out_A, sprintf(" %6s", "?"));  # stop
        push(@cur_out_A, sprintf(" %6s", "?"));  # length
      }
      else { # we know that $mdl_results_AAHR->[0][$seq_idx]{"p_start"} exists) { 
        # determine output start and output stop
        my $cur_stop = undef; # final stop position
        if(exists $mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"append_stop"}) { 
          (undef, $cur_stop) = create_output_start_and_stop($mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"append_start"}, 
                                                            $mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"append_stop"},
                                                            $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR);
        }
        elsif(exists $mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"c_stop"}) { 
          (undef, $cur_stop) = create_output_start_and_stop($mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_start"},  # irrelevant due to the first undef arg
                                                            $mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"c_stop"},
                                                            $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR);
        }
        else { 
          (undef, $cur_stop) = create_output_start_and_stop($mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_start"}, # irrelevant due to the first undef arg
                                                            $mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_stop"}, 
                                                            $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR);
        }
        if($mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_strand"} eq "+") { 
          # positive strand, easy case
          if($cur_stop == $accn_len) { # final model prediction stops at final nt
            push(@cur_out_A, sprintf("  %6d", 0)); # start 
            push(@cur_out_A, sprintf("  %6d", 0)); # stop 
            push(@cur_out_A, sprintf("  %6d", 0)); # length
          }            
          else { 
            push(@cur_out_A, sprintf("  %6d", $cur_stop + 1));           # start
            push(@cur_out_A, sprintf("  %6d", $accn_len));               # stop 
            push(@cur_out_A, sprintf("  %6d", ($accn_len - $cur_stop))); # length
          }
        }
        elsif($mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_strand"} eq "-") { 
          # negative strand, more complicated, slightly
          if($cur_stop == 1) { # final model prediction stops at first nt
            push(@cur_out_A, sprintf("  %6d", 0)); # start 
            push(@cur_out_A, sprintf("  %6d", 0)); # stop 
            push(@cur_out_A, sprintf("  %6d", 0)); # length
          }            
          else { # final feature does not stop at nt 1 on negative strand
            push(@cur_out_A, sprintf("  %6d", 1)); # start 
            push(@cur_out_A, sprintf("  %6d", $cur_stop - 1)); # stop 
            push(@cur_out_A, sprintf("  %6d", $cur_stop - 1)); # length
          }
        }
        else { # not + or - strand, weird...
          DNAORG_FAIL("ERROR in $sub_name, trying to compute 5' UTR for prediction that exists but is not + or - strand", 1, $FH_HR);
        }
      }
    }

    # total length
    push(@cur_out_A, sprintf("  %6d", $accn_len));
    # average fid
    if($do_totfid) { 
      push(@cur_out_A, sprintf("  %5.3f", ($tot_id_min_len > 0) ? ($tot_nid / $tot_id_min_len) : 0.));
    }
    if($do_exist) { 
      # output stats on GenBank annotations and comparison
      push(@cur_out_A, sprintf("  %5d", $seq_info_HAR->{"num_genbank_mdl_annot"}[$seq_idx]));      # number of GenBank annotated features
      push(@cur_out_A, sprintf("  %5d", $seq_info_HAR->{"num_genbank_mdl_exon_annot"}[$seq_idx])); # number of exons in those annotate features
      push(@cur_out_A, sprintf("  %5d", $ngenbank_match)); # number of exons for which our annotation matches start..stop of at least 1 GenBank exon
    }

    my $accn_failed = ($pass_fail_str =~ m/F/) ? 1 : 0;
    my $result_str = ($accn_failed) ? "FAIL" : "PASS";
    $result_str .= " " . $pass_fail_str;
    push(@cur_out_A, sprintf("  %s", $result_str));

    # make a version of $pass_fail_str with spaces for the $tblsum_FH file
    my @pass_fail_A = split("", $pass_fail_str);
    my $pass_fail_str_with_spaces = ($accn_failed) ? "FAIL" : "PASS";
    $pass_fail_str_with_spaces .= " " . $pass_fail_A[0];
    for(my $pf = 1; $pf < scalar(@pass_fail_A); $pf++) { 
      $pass_fail_str_with_spaces .= " " . $pass_fail_A[$pf];
    }    
    print $tblsum_FH ("$accn_name $pass_fail_str_with_spaces\n");

    # actually output the information to the relevant file handles
    if(($seq_idx == 0) && $do_tblfirst) { 
      # copy reference info if this is the reference
      @ref_out_A = @cur_out_A; 
    } 
    else { 
      push(@page_out_AA, [@cur_out_A]);
      $cur_pagesize++;
      if(($accn_failed) || (($seq_idx == 0) && (! opt_IsUsed("--infasta", $opt_HHR)))) { # print ref if --infasta not used
        push(@fail_page_out_AA, [@cur_out_A]);
        $cur_fail_pagesize++;
        $tot_nfail++;
      }        
      if(($seq_info_HAR->{"nerrors"}[$seq_idx] > 0) || (($seq_idx == 0) && (! opt_IsUsed("--infasta", $opt_HHR)))) { # print ref if --infasta not used { 
        push(@err_page_out_AA, [@cur_out_A]);
        $cur_err_pagesize++;
      }        
    }
    if($cur_pagesize == $act_nseqcol) { 
      $npages++;
      output_tbl_page_of_sequences($FH_HR->{"tbl"}, \@out_row_header_A, \@page_out_AA, ($do_tblfirst ? \@ref_out_A : undef), $npages, $FH_HR);
      @page_out_AA = ();
      $cur_pagesize = 0;
    }
    if($accn_failed && 
       ($cur_fail_pagesize == $act_nseqcol)) { 
      $nfail_pages++;
      output_tbl_page_of_sequences($FH_HR->{"failtbl"}, \@out_row_header_A, \@fail_page_out_AA, ($do_tblfirst ? \@ref_out_A : undef), $nfail_pages, $FH_HR);
      @fail_page_out_AA = ();
      $cur_fail_pagesize = 0;
    }
    if(($seq_info_HAR->{"nerrors"}[$seq_idx] > 0) && 
       ($cur_err_pagesize == $act_nseqcol)) { 
      $nerr_pages++;
      output_tbl_page_of_sequences($FH_HR->{"errtbl"}, \@out_row_header_A, \@err_page_out_AA, ($do_tblfirst ? \@ref_out_A : undef), $nerr_pages, $FH_HR);
      @err_page_out_AA = ();
      $cur_err_pagesize = 0;
    }
  } # end of 'for($seq_idx'
  # print final page (if non-empty)
  if($cur_pagesize > 0) { 
    $npages++;
    output_tbl_page_of_sequences($FH_HR->{"tbl"}, \@out_row_header_A, \@page_out_AA, ($do_tblfirst ? \@ref_out_A : undef), $npages, $FH_HR);
  }
  if($cur_fail_pagesize > 0) { 
    $nfail_pages++;
    output_tbl_page_of_sequences($FH_HR->{"failtbl"}, \@out_row_header_A, \@fail_page_out_AA, ($do_tblfirst ? \@ref_out_A : undef), $nfail_pages, $FH_HR);
  }
  if($cur_err_pagesize > 0) { 
    $nerr_pages++;
    output_tbl_page_of_sequences($FH_HR->{"errtbl"}, \@out_row_header_A, \@err_page_out_AA, ($do_tblfirst ? \@ref_out_A : undef), $nerr_pages, $FH_HR);
  }

  return $tot_nfail;
}

#################################################################
# Subroutine: output_tbl_page_of_sequences()
# Incept:     EPN, Mon Mar 14 14:53:30 2016
#
# Purpose:    Output a 'page' of annotation information. 
#
# Args:
#  $FH:            file handle to print to
#  $header_AR:     reference to array of row headers
#  $out_AAR:       reference to the 2D array of output tokens for
#                  current sequences for the page
#  $ref_out_AR:    reference to array of output tokens for 1st column 
#                  (e.g. the reference), undef to skip
#  $page_idx:      page number
#  $FH_HR:         REF to hash of file handles
#
# Returns:    void
#
# Dies: If we don't have the appropriate number of tokens in an output array.
#
#################################################################
sub output_tbl_page_of_sequences { 
  my $sub_name = "output_tbl_page_of_sequences";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($FH, $header_AR, $out_AAR, $ref_out_AR, $page_idx, $FH_HR) = @_;

  my $nseq = scalar(@{$out_AAR});
  if(defined $ref_out_AR) { $nseq++; }
  my $AR; # reference to current array we are printing

  my @cwidth_A = (); # max width of each column (max width of all tokens for each sequence) 
  my ($ntok, $el); # number of tokens, and element of an array
  my $nrow = scalar(@{$header_AR});

  # first, make sure all sequences have same number of output tokens as we have row headers,
  # and determine maximum width of all tokens for each sequence
  for(my $i = 0; $i < $nseq; $i++) { 
    # if $ref_out_AR is defined: first column is reference, then come the other seqs
    my $ip = (defined $ref_out_AR) ? $i-1 : $i;
    $AR = ((defined $ref_out_AR) && ($i == 0)) ? $ref_out_AR : \@{$out_AAR->[$ip]}; 
    $cwidth_A[$i] = 0;
    $ntok = scalar(@{$AR});
    if($ntok != $nrow) { 
      # if you ever get the error in the following line, comment it out and rerun
      # the output will be shifted by some number of tokens and it should be helpful
      # for figuring out what tokens are missing
      DNAORG_FAIL(sprintf("ERROR in $sub_name, we have $nrow headers, but sequence %s has $ntok tokens", $i+1, $ntok), 1, $FH_HR); 
    }
    foreach $el (@{$AR}) { 
      $el =~ s/^\s+//; # remove leading whitespace
      $el =~ s/\s+$//; # remove trailing whitespace
      if(length($el) > $cwidth_A[$i]) { 
        $cwidth_A[$i] = length($el);
      }
    }
  }
  for(my $i = 0; $i < $nseq; $i++) { 
    $cwidth_A[$i] += 2; # add 2 spaces for in-between sequence columns (the two spaces to the right of each column)
  }

  # determine max width of all row headers:
  my $hwidth = 0;
  for(my $r = 0; $r < $nrow; $r++) { 
    if(length($header_AR->[$r]) > $hwidth) { 
      $hwidth = length($header_AR->[$r]);
    }
  }
  $hwidth += 2;

  for(my $r = 0; $r < $nrow; $r++) { 
    printf $FH ("%-*s", $hwidth, $header_AR->[$r]);
    for(my $i = 0; $i < $nseq; $i++) { 
      # if $ref_out_AR is defined: first column is reference, then come the other seqs
      my $ip = (defined $ref_out_AR) ? $i-1 : $i;
      $AR = ((defined $ref_out_AR) && ($i == 0)) ? $ref_out_AR : \@{$out_AAR->[$ip]}; 
      $el = $AR->[$r];
      $el =~ s/\s+//g;
      printf $FH ("%*s", $cwidth_A[$i], $AR->[$r]);
    }
    printf $FH ("\n");
  }
  print $FH "#\n";
  printf $FH ("# end of page %d\n", $page_idx);
  print $FH "#\n";

  return;
}

#################################################################
# Subroutine:  output_errors_header
# Incept:      EPN, Thu Mar 10 18:57:48 2016
#
# Purpose:    Output the header lines for the all errors and 
#             per-accession error files.
#
# Arguments: 
#  $ftr_info_HAR:   REF to hash of arrays with information on the features
#  $ofile_info_HHR: REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
################################################################# 
sub output_errors_header { 
  my $sub_name = "output_errors_header";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $ofile_info_HHR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nftr: number of features

  my $allerr_FH = $ofile_info_HHR->{"FH"}{"allerr"};
  my $pererr_FH = $ofile_info_HHR->{"FH"}{"pererr"};

  # all errors file header
  printf $allerr_FH ("# Each error encountered is printed below, one error per line.\n");
  printf $allerr_FH ("# Each line has four columns with the following labels:\n");
  printf $allerr_FH ("#   \"accn\"         : sequence accession\n");
  printf $allerr_FH ("#   \"idx\"          : feature index, full feature names are listed below for each index\n");
  printf $allerr_FH ("#   \"code\"         : 3 digit error code\n");
  printf $allerr_FH ("#   \"error-message\": error message, possibly with additional information at end enclosed in \"[]\"\n");
  printf $allerr_FH ("#\n");
  printf $allerr_FH ("# List of features:\n");

  # per accession errors file header
  printf $pererr_FH ("# Each accession for which at least one error was found is printed below.\n");
  printf $pererr_FH ("# One line per accession. Each line is formatted as follows:\n");
  printf $pererr_FH ("#   <accession> <idxA>:<errorcodeA1>(,<errorcodeAM>) <idxN>:<errorcodeN1>(,<errorcodeNM>)\n");
  printf $pererr_FH ("# For indices (<idx>) A to N, each with M error codes.\n");
  printf $pererr_FH ("# Each index refers to a 'feature' in the reference accession as defined below.\n");
  printf $pererr_FH ("# If no index exists, then the error pertains to the entire sequence.\n");

  # print feature list
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $str = sprintf("# Feature \#%d: %s %s\n", $ftr_idx+1, $ftr_info_HAR->{"out_short"}[$ftr_idx], $ftr_info_HAR->{"out_product"}[$ftr_idx]);
    print $allerr_FH $str;
    print $pererr_FH $str;
  }

  # a few more all error lines
  printf $allerr_FH ("# \"N/A\" in feature index and desc columns (idx and desc) indicates error pertains to the entire sequence\n");
  printf $allerr_FH ("#       as opposed to a specific feature.\n");
  printf $allerr_FH ("#\n");

  # output multifeature relationships (e.g. CDS made up of mature peptides)
  output_multifeature_relationships($pererr_FH, $ftr_info_HAR, $ofile_info_HH{"FH"});
  output_multifeature_relationships($allerr_FH, $ftr_info_HAR, $ofile_info_HH{"FH"});

  printf $allerr_FH ("%-10s  %3s  %-9s  %4s  error-message\n", "#accn", "idx", "desc", "code");
  printf $pererr_FH ("#\n");

  return;
}

#################################################################
# Subroutine: output_errors_all_sequences()
# Incept:     EPN, Tue Mar 15 14:33:52 2016
#
# Purpose:   Output the errors for all sequences.
#
# Arguments:
#  $err_ftr_instances_AHHR: REF to array of 2D hashes with per-feature errors, PRE-FILLED
#  $err_seq_instances_HHR:  REF to 2D hash with per-sequence errors, PRE-FILLED
#  $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:           REF to hash of arrays with information on the sequences, PRE-FILLED, we add "nerrors" values here
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:                REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:         REF to the 2D hash of output file information
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub output_errors_all_sequences { 
  my $sub_name = "output_errors_all_sequences";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_ftr_instances_AHHR, $err_seq_instances_HHR, $ftr_info_HAR, $seq_info_HAR, $err_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"}; # for convenience
  my $all_FH = $FH_HR->{"allerr"}; # for convenience
  my $per_FH = $FH_HR->{"pererr"}; # for convenience

  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $nerr = getConsistentSizeOfInfoHashOfArrays($err_info_HAR, $FH_HR); 

  my ($seq_idx, $ftr_idx, $err_idx); # counters
  for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_nseqerr = 0; # number of per-sequence errors for this sequence
    my $seq_nftrerr = 0; # number of per-feature errors for this sequence
    my $seq_name    = $seq_info_HAR->{"seq_name"}[$seq_idx];
    my $accn_name   = $seq_info_HAR->{"accn_name"}[$seq_idx];
    my $per_line    = $accn_name . " "; # the line for this sequence to print to $per_FH
    #######################
    # per-sequence errors #
    #######################
    for($err_idx = 0; $err_idx < $nerr; $err_idx++) { 
      if($err_info_HAR->{"pertype"}[$err_idx] eq "sequence") { 
        my $err_code = $err_info_HAR->{"code"}[$err_idx];
        if(exists $err_seq_instances_HHR->{$err_code}{$seq_name}) { 
          # an error exists, output it
          printf $all_FH ("%-10s  %3s  %-9s  %4s  %s%s\n", $accn_name, "N/A", "N/A", $err_code, $err_info_HAR->{"msg"}[$err_idx], 
                          " [" . $err_seq_instances_HHR->{$err_code}{$seq_name} . "]"); 
          $seq_nseqerr++;
          $per_line .= " " . $err_code;
        }
      }  
    }
    ######################
    # per-feature errors #
    ######################
    for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      my $out_tiny = $ftr_info_HAR->{"out_tiny"}[$ftr_idx];
      my $ftr_seq_nftrerr = 0; # number of errors we've seen for this sequence and feature
      for($err_idx = 0; $err_idx < $nerr; $err_idx++) { 
        if($err_info_HAR->{"pertype"}[$err_idx] eq "feature") { 
          my $err_code = $err_info_HAR->{"code"}[$err_idx];
          if(exists $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}) { 
            # an error exists, output it
            printf $all_FH ("%-10s  %3s  %-9s  %4s  %s%s\n", $accn_name, ($ftr_idx+1), $out_tiny, $err_code, $err_info_HAR->{"msg"}[$err_idx], 
                            ($err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} eq "") ? "" : " [" . $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} . "]"); 
            if($ftr_seq_nftrerr == 0) { 
              $per_line .= (" " . ($ftr_idx+1) . ":" . $err_code);
            }
            else { 
              $per_line .= "," . $err_code;
            }
            if($err_code eq "olp") { # special case, add the extra string
              $per_line .= ":" . $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name};
            }
            $ftr_seq_nftrerr++;
            $seq_nftrerr++;
          }
        }
      }
    }
    if($seq_nseqerr + $seq_nftrerr > 0) { 
      print $per_FH $per_line . "\n";
    }
    $seq_info_HAR->{"nerrors"}[$seq_idx] = $seq_nseqerr + $seq_nftrerr;
  } # end of 'for($seq_idx = 0'...
  return;
}

#################################################################
# Subroutine: output_errors_summary()
# Incept:     EPN, Thu Mar 17 12:55:55 2016
#
# Purpose:   Summarize the errors. Create this in a new file
#            and also (optionally) to stdout, if $do_stdout.
#
# Arguments:
#  $FH:                     file handle to print to
#  $err_ftr_instances_AHHR: REF to array of 2D hashes with per-feature errors, PRE-FILLED
#  $err_seq_instances_HHR:  REF to 2D hash with per-sequence errors, PRE-FILLED
#  $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:           REF to hash of arrays with information on the sequences, PRE-FILLED
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $do_stdout:              '1' to output to stdout as well as $FH file handle
#  $opt_HHR:                REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:         REF to the 2D hash of output file information
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub output_errors_summary { 
  my $sub_name = "output_errors_summary()";
  my $nargs_exp = 9;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($FH, $err_ftr_instances_AHHR, $err_seq_instances_HHR, $ftr_info_HAR, $seq_info_HAR, $err_info_HAR, $do_stdout, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"}; # for convenience
  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $nerr = getConsistentSizeOfInfoHashOfArrays($err_info_HAR, $FH_HR); 

  # output header, with explanations
  outputString($FH, $do_stdout, sprintf("# Table below includes counts of error codes\n"));
  outputString($FH, $do_stdout, sprintf("#\n"));
  outputString($FH, $do_stdout, sprintf("# Explanation of columns:\n"));
  outputString($FH, $do_stdout, sprintf("# \"code\"       : the error code\n"));
  outputString($FH, $do_stdout, sprintf("# \"#tot\"       : total number of occurences of code , possibly > 1 for some accessions\n"));
  outputString($FH, $do_stdout, sprintf("# \"#accn\"      : number of accessions with at least 1 occurence of code\n"));
  outputString($FH, $do_stdout, sprintf("# \"fraction...\": fraction of all $nseq accessions with at least 1 occurence of code\n"));
  outputString($FH, $do_stdout, sprintf("#\n"));

  outputString($FH, $do_stdout, sprintf("# Explanation of final two rows beginning with \"total\" and \"any\":\n"));
  outputString($FH, $do_stdout, sprintf("#   \"total\":\"#tot\"   column is total number of error codes reported\n"));
  outputString($FH, $do_stdout, sprintf("#   \"any\":\"#tot\"     column is number of accessions with >= 1 error code\n"));
  outputString($FH, $do_stdout, sprintf("#   \"any\":\"fraction\" column is fraction of accessions with >= 1 error code\n"));
  outputString($FH, $do_stdout, sprintf("#code   #tot  #accn  fraction-of-all-$nseq-accn\n"));
  outputString($FH, $do_stdout, sprintf("#----  -----  -----  ------\n"));

  # for each error code, gather counts
  my $ntot_all  = 0; # total number of errors
  my $naccn_all = 0; # number of accessions with >= 1 error 
  my %seq_counted_all_err_H = (); # for all errors, we add each sequence to this hash if it has >= 1 error
  
  for(my $err_idx = 0; $err_idx < $nerr; $err_idx++) { 
    my $ntot_err  = 0;
    my $naccn_err = 0;
    my $err_code = $err_info_HAR->{"code"}[$err_idx];
    my %seq_counted_this_err_H = (); # for per-feature errors, we add each sequence to this hash 
                                     # as we see it, so we can more efficiently count the errors

    if($err_info_HAR->{"pertype"}[$err_idx] eq "sequence") { 
      foreach my $seq_name (keys %{$err_seq_instances_HHR->{$err_code}}) { 
        $naccn_err++;
        if(! exists $seq_counted_all_err_H{$seq_name}) { 
          $naccn_all++;
          $seq_counted_all_err_H{$seq_name} = 1;
        }
      }
      $ntot_err  = $naccn_err; # only 1 error per sequence for these guys
      $ntot_all += $naccn_err; # only 1 error per sequence for these guys
    }
    elsif($err_info_HAR->{"pertype"}[$err_idx] eq "feature") { 
      for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
        foreach my $seq_name (keys %{$err_ftr_instances_AHHR->[$ftr_idx]{$err_code}}) { 
          $ntot_err++;
          $ntot_all++;
          if(! exists $seq_counted_this_err_H{$seq_name}) { 
            $naccn_err++;
            $seq_counted_this_err_H{$seq_name} = 1;
          }
          if(! exists $seq_counted_all_err_H{$seq_name}) { 
            $naccn_all++;
            $seq_counted_all_err_H{$seq_name} = 1;
          }
        }
      }
    }

    # output stats for this error 
    outputString($FH, $do_stdout, sprintf("%5s  %5d  %5d  %6.4f\n", $err_code, $ntot_err, $naccn_err, ($naccn_err / $nseq)));
  }

  # the 'total' and 'any' lines
  outputString($FH, $do_stdout, sprintf("#----  -----  -----  ------\n"));
  outputString($FH, $do_stdout, sprintf("%5s  %5d  %5s  %6s\n",  "total", $ntot_all,  "-", "-"));
  outputString($FH, $do_stdout, sprintf("%5s  %5d  %5s  %6.4f\n", "any",  $naccn_all, "-", ($naccn_all / $nseq)));

  return;
}

#################################################################
# Subroutine:  output_multifeature_relationships
# Incept:      EPN, Thu Mar 10 19:15:38 2016
#
# Purpose:    Output the multi-feature relationships, e.g. CDS comprised of
#             multiple features.
#
# Arguments: 
#  $FH:             file handle to print to
#  $ftr_info_HAR:   REF to hash of arrays with information on the features
#  $FH_HR:          REF to hash of file handles
#
# Returns:    void
#
################################################################# 
sub output_multifeature_relationships { 
  my $sub_name = "output_multifeature_relationships";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($FH, $ftr_info_HAR, $FH_HR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $FH_HR); # nftr: number of features

  my $nprinted = 0;

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") { # we only do this for features annotated by models
      if($nprinted == 0) { 
        print $FH ("#\n");
        print $FH ("# CDS:MAT_PEPTIDE relationships:\n");
        print $FH ("#\n");
      }
      
      if($ftr_info_HAR->{"type"}[$ftr_idx] ne "cds-mp") { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, unexpected type %s for multifeature, feature #%d %s, expected cds-mp", $ftr_info_HAR->{"type"}[$ftr_idx], $ftr_idx+1, $ftr_info_HAR->{"out_tiny"}[$ftr_idx]), 1, $FH_HR);
      }
      
      # get the array of primary children feature indices for this feature
      my @children_idx_A = (); # feature indices of the primary children of this feature
      getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "primary", \@children_idx_A, $FH_HR);
      printf $FH ("# %s is comprised of the following primary features in order:\n#   ", $ftr_info_HAR->{"out_tiny"}[$ftr_idx]);
      foreach my $child_ftr_idx (@children_idx_A) { 
        printf $FH "%s ", $ftr_info_HAR->{"out_tiny"}[$child_ftr_idx];
      }
      print $FH "\n#\n";
      
      getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "all", \@children_idx_A, $FH_HR);
      printf $FH ("# %s encodes all of the following features in order:\n#   ", $ftr_info_HAR->{"out_tiny"}[$ftr_idx]);
      foreach my $child_ftr_idx (@children_idx_A) { 
        printf $FH "%s ", $ftr_info_HAR->{"out_tiny"}[$child_ftr_idx];
      }
      print $FH "\n#\n";
      $nprinted++;
    }
  }
  outputDividingLine(undef, $FH); # undef makes outputDividingLine() use its default length for the dividing line

  return;
}

#################################################################
# Subroutine: output_gap_info()
# Incept:     EPN, Thu Mar 10 19:15:38 2016
#
# Purpose:    Output the multi-feature relationships, e.g. CDS comprised of
#             multiple features.
#
# Arguments: 
# $perseq_FH:        output file handle to print per-sequence gap info to, undef to not print perseq info
# $pergap_FH:        output file handle to print per-gap info to
# $do_perseq_tbl:    '1' to output per sequence gaps as a table, '0' to print a list
# $do_gap_all:       '1' to output per gap info for all gaps
# $do_gap_not3:      '1' to output per gap info for gaps that are not a multiple of 3, not for all gaps
# $do_gap_special:   '1' to output per gap info for special gaps that are possibly causative of a frameshift
#                    Only 1 of $pergap_all, $pergap_not3, and $pergap_special can be '1'.
# $mdl_info_HAR:     REF to hash of arrays with information on the models, PRE-FILLED
# $ftr_info_HAR:     REF to hash of arrays with information on the features, PRE-FILLED
# $seq_info_HAR:     REF to hash of arrays with information on the sequences, PRE-FILLED
# $mdl_results_AAHR: REF to results AAH, PRE-FILLED
# $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
# $ofile_info_HHR:   REF to the 2D hash of output file information
#
# Returns:    void
#
#################################################################
sub output_gap_info {
  my $sub_name = "output_gap_info";
  my $nargs_exp = 12;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($perseq_FH, $pergap_FH, $do_perseq_tbl, $do_gap_all, $do_gap_not3, $do_gap_special, $mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience

  if($do_gap_all) {
    if($do_gap_not3 || $do_gap_special) { die "ERROR in $sub_name, exactly one of $do_gap_all, $do_gap_not3, and $do_gap_special must be true"; }
  }
  elsif($do_gap_not3) { 
    if($do_gap_all || $do_gap_special) { die "ERROR in $sub_name, exactly one of $do_gap_all, $do_gap_not3, and $do_gap_special must be true"; }
  }
  elsif($do_gap_special) { 
    if($do_gap_all || $do_gap_not3) { die "ERROR in $sub_name, exactly one of $do_gap_all, $do_gap_not3, and $do_gap_special must be true"; }
  }
  else { 
    die "ERROR in $sub_name, exactly one of $do_gap_all, $do_gap_not3, and $do_gap_special must be true"; 
  }

  my @gapstr_AA = ();           # [0..$nseq-1][0..$nftr-1]: string describing all gaps for this sequence and this feature
  my @w_gapstr_A = ();          # [0..$i..$nftr-1] max width of $gapstr_AA[$0..nseq-1][$i] for feature $i over all sequences

  my @tot_gap_length_AA = ();   # [0..$nseq-1][0..$nftr-1]: total number of gap positions for this sequence and this feature
  my @w_tot_gap_length_A = ();  # [0..$i..$nftr-1] max width of $tot_gap_length_AA[$0..nseq-1][$i] for feature $i over all sequences

  my @net_gap_length_AA = ();   # [0..$nseq-1][0..$nftr-1]: net number of gap positions for this sequence and this feature
  my @w_net_gap_length_A = ();  # [0..$i..$nftr-1] max width of $net_gap_length_AA[$0..nseq-1][$i] for feature $i over all sequences

  my @ftr_gapstr_AH      = ();  # [0..$i..$nftr-1]: hash w/key: gapstring for a single position, value: number of times that gapstring occurs in any of @gapstr_AA for feature $i
  
  my $ch_gapstr         = "string";
  my $ch_tot_gap_length = "tot";
  my $ch_net_gap_length = "net";

  my $width_seq = length("#accession");
  my $ftr_idx = 0;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences

  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    $w_gapstr_A[$ftr_idx]         = length($ch_gapstr);
    $w_tot_gap_length_A[$ftr_idx] = length($ch_tot_gap_length);
    $w_net_gap_length_A[$ftr_idx] = length($ch_net_gap_length);
    %{$ftr_gapstr_AH[$ftr_idx]}   = ();
  }

  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    # initialize 

    my $seq_accn = $seq_info_HAR->{"accn_name"}[$seq_idx];
    if(length($seq_accn) > $width_seq) { 
      $width_seq = length($seq_accn);
    }
    @{$gapstr_AA[$seq_idx]} = ();
    my $offset = 0;
    my $ins_idx = 0;
    my $del_idx = 0;
    my $next_ins_str   = undef;
    my $next_del_str   = undef;
    my $next_ins_rfpos = undef;
    my $next_del_rfpos = undef;
    my $next_ins_count = undef;
    my $next_del_count = undef;
    my $gapstr = "";
    my $substr;
    my $tot_gap_length = 0;
    my $net_gap_length = 0;
    $ftr_idx = 0;

    for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
      my @refdel_A = (); 
      my @refins_A = (); 
      if((exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refinsstr"}) && $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refinsstr"} ne "") {
        @refins_A = split(",", $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refinsstr"});
      }
      if((exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refdelstr"}) && $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refdelstr"} ne "") {
        @refdel_A = split(",", $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refdelstr"});
      }
      my $ndel = scalar(@refdel_A);
      my $nins = scalar(@refins_A);
      my $mdl_is_final = $mdl_info_HAR->{"is_final"}[$mdl_idx];
      my $mdl_is_first = $mdl_info_HAR->{"is_first"}[$mdl_idx];
      $ins_idx = 0;
      $del_idx = 0;
      $next_del_str = ($del_idx < $ndel) ? $refdel_A[$del_idx] : undef;
      $next_ins_str = ($ins_idx < $nins) ? $refins_A[$ins_idx] : undef;
      
      while(defined $next_ins_str || defined $next_del_str) { 
        # printf("next_ins_str: %s\n", (defined $next_ins_str) ? $next_ins_str : "undefined");
        # printf("next_del_str: %s\n", (defined $next_del_str) ? $next_del_str : "undefined");
        ($next_del_rfpos, $next_del_count) = (defined $next_del_str) ? split(":", $next_del_str) : (undef, undef);
        ($next_ins_rfpos, $next_ins_count) = (defined $next_ins_str) ? split(":", $next_ins_str) : (undef, undef);
        if(defined $next_del_rfpos) { $next_del_rfpos += $offset; }
        if(defined $next_ins_rfpos) { $next_ins_rfpos += $offset; }
        
        if(defined $next_ins_str && defined $next_del_str) { 
          if($next_del_rfpos <= $next_ins_rfpos) { # delete comes first, print it
            $substr = "D" . $next_del_rfpos . ":" . $next_del_count;
            if($do_gap_all || (($next_del_count % 3) != 0)) { 
              if($gapstr ne "") { $gapstr .= ","; }
              $gapstr .= $substr;
              $tot_gap_length += $next_del_count;
              $net_gap_length -= $next_del_count;
              if(! $do_gap_special) { 
                $ftr_gapstr_AH[$ftr_idx]{$substr}++;
              }
            }
            $del_idx++;
            $next_del_str = ($del_idx < $ndel) ? $refdel_A[$del_idx] : undef;
          }
          elsif($next_ins_rfpos < $next_del_rfpos) { # insert comes first, print it
            $substr = "I" . $next_ins_rfpos . ":" . $next_ins_count;
            if($do_gap_all || (($next_ins_count % 3) != 0)) { 
              if($gapstr ne "") { $gapstr .= ","; }
              $gapstr .= $substr;
              $tot_gap_length += $next_ins_count;
              $net_gap_length += $next_ins_count;
              if(! $do_gap_special) { 
                $ftr_gapstr_AH[$ftr_idx]{$substr}++;
              }
            }
            $ins_idx++;
            $next_ins_str = ($ins_idx < $nins) ? $refins_A[$ins_idx] : undef;
          }
        }
        elsif(defined $next_del_str) { # $next_ins is undefined
          $substr = "D" . $next_del_rfpos . ":" . $next_del_count;
          if($do_gap_all || (($next_del_count % 3) != 0)) { 
            if($gapstr ne "") { $gapstr .= ","; }
            $gapstr .= $substr;
            $tot_gap_length += $next_del_count;
            $net_gap_length -= $next_del_count;
            if(! $do_gap_special) {
              $ftr_gapstr_AH[$ftr_idx]{$substr}++;
            }
          }
          $del_idx++;
          $next_del_str = ($del_idx < $ndel) ? $refdel_A[$del_idx] : undef;
        }
        elsif(defined $next_ins_str) { # $next_del is undefined
          $substr = "I" . $next_ins_rfpos . ":" . $next_ins_count;
          if($do_gap_all | (($next_ins_count % 3) != 0))  { 
            if($gapstr ne "") { $gapstr .= ","; }
            $gapstr .= $substr;
            $tot_gap_length += $next_ins_count;
            $net_gap_length += $next_ins_count;
            if(! $do_gap_special) { 
              $ftr_gapstr_AH[$ftr_idx]{$substr}++;
            }
          }
          $ins_idx++;
          $next_ins_str = ($ins_idx < $nins) ? $refins_A[$ins_idx] : undef;
        }        
      } # end of 'while(defined $next_ins_str || defined $next_del_str) { 
      
      $offset += $mdl_info_HAR->{"length"}[$mdl_idx];
      
      if($mdl_is_final) {
        if($gapstr eq "") { $gapstr = "-"; }

        # important to update $gapstr here, before we store it if $do_gap_special is true
        if($do_gap_special) { 
          # printf("calling find_special_gap with $gapstr\n");
          $gapstr = find_special_gap($gapstr);
          $net_gap_length = $tot_gap_length;
          if($gapstr ne "?" && $gapstr ne "-") { 
            my @el_A = split(",", $gapstr); 
            foreach my $el (@el_A) { 
              $ftr_gapstr_AH[$ftr_idx]{$el}++;
            }
          }
        }

        push(@{$gapstr_AA[$seq_idx]}, $gapstr);
        if(length($gapstr) > $w_gapstr_A[$ftr_idx]) { $w_gapstr_A[$ftr_idx] = length($gapstr); }

        push(@{$tot_gap_length_AA[$seq_idx]}, $tot_gap_length);
        if(length($tot_gap_length) > $w_tot_gap_length_A[$ftr_idx]) { $w_tot_gap_length_A[$ftr_idx] = length($tot_gap_length); }
        $tot_gap_length = 0;

        push(@{$net_gap_length_AA[$seq_idx]}, $net_gap_length);
        if(length($net_gap_length) > $w_net_gap_length_A[$ftr_idx]) { $w_net_gap_length_A[$ftr_idx] = length($net_gap_length); }
        $net_gap_length = 0;

        $gapstr = "";
        $offset = 0;
        $ftr_idx++;

        if(scalar(@w_gapstr_A)         <= $ftr_idx) { $w_gapstr_A[$ftr_idx]         = length($ch_gapstr); }
        if(scalar(@w_tot_gap_length_A) <= $ftr_idx) { $w_tot_gap_length_A[$ftr_idx] = length($ch_tot_gap_length); }
        if(scalar(@w_net_gap_length_A) <= $ftr_idx) { $w_net_gap_length_A[$ftr_idx] = length($ch_net_gap_length); }
        if(scalar(@ftr_gapstr_AH)      <= $ftr_idx) { %{$ftr_gapstr_AH[$ftr_idx]}   = (); }
      }
    }
  }

  #################################################################
  # Output table of all gap info
  # output line 1 of the column headers:
  if($do_gap_all) { 
    printf $perseq_FH ("# List of all gaps in alignment of each sequence:\n#\n");
  }
  elsif($do_gap_not3) { 
    printf $perseq_FH ("# List of all gaps with length that is not a multiple of 3 in alignment of each sequence:\n#\n");
  }
  else { # $do_gap_special 
    printf $perseq_FH ("# List of all gaps that may solely explain a feature's (CDS or mat_peptide) length not being a multiple of 3, for each sequence:\n#\n");
  }

  if(! $do_gap_special) { 
    printf $perseq_FH  ("%-*s  ", $width_seq, "#");
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      my $is_matpept = ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp") ? 1 : 0;
      my $type_idx   = $ftr_info_HAR->{"type_idx"}[$ftr_idx];
      my $w_cur = $w_tot_gap_length_A[$ftr_idx] + 2 + $w_net_gap_length_A[$ftr_idx] + 2 + $w_gapstr_A[$ftr_idx];
      if($ftr_idx > 0) { print $perseq_FH "  "; }
      printf $perseq_FH ("%-*s", $w_cur, ($is_matpept) ? ("mat_peptide#" . $type_idx) : ("CDS#" . $type_idx));
    }
    print $perseq_FH "\n";
    
    # output line 2 (dashes under line 1 of column headers)
    printf $perseq_FH ("%-*s  ", $width_seq, "#");
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      my $w_cur = $w_tot_gap_length_A[$ftr_idx] + 2 + $w_net_gap_length_A[$ftr_idx] + 2 + $w_gapstr_A[$ftr_idx];
      if($ftr_idx > 0) { print $perseq_FH "  "; }
      printf $perseq_FH ("%-*s", $w_cur, getMonocharacterString($w_cur, "=", $FH_HR));
    }    
    print $perseq_FH "\n";
  }

  # output line 3 of the column headers:
  printf $perseq_FH ("%-*s  ", $width_seq, "#accession");
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $is_matpept = ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp") ? 1 : 0;
    my $type_idx   = $ftr_info_HAR->{"type_idx"}[$ftr_idx];
    if($ftr_idx > 0) { print $perseq_FH "  "; }
    if(! $do_gap_special) { 
      printf $perseq_FH ("%-*s  ", $w_tot_gap_length_A[$ftr_idx], $ch_tot_gap_length);
      printf $perseq_FH ("%-*s  ", $w_net_gap_length_A[$ftr_idx], $ch_net_gap_length);
      printf $perseq_FH ("%-*s", $w_gapstr_A[$ftr_idx], $ch_gapstr);
    }
    else { 
      printf $perseq_FH ("%-*s", $w_gapstr_A[$ftr_idx], ($is_matpept) ? ("mat_peptide#" . $type_idx) : ("CDS#" . $type_idx));
    }
  }
  print $perseq_FH "\n";

  # output line 4 (dashes under line 3 of column headers)
  printf $perseq_FH ("%-*s  ", $width_seq, "#" . getMonocharacterString($width_seq-1, "-", $FH_HR));
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_idx > 0) { print $perseq_FH "  "; }
    if(!  $do_gap_special) { 
      printf $perseq_FH ("%-*s  ", $w_tot_gap_length_A[$ftr_idx], getMonocharacterString($w_tot_gap_length_A[$ftr_idx], "-", $FH_HR));
      printf $perseq_FH ("%-*s  ", $w_net_gap_length_A[$ftr_idx], getMonocharacterString($w_net_gap_length_A[$ftr_idx], "-", $FH_HR));
    }
    printf $perseq_FH ("%-*s", $w_gapstr_A[$ftr_idx], getMonocharacterString($w_gapstr_A[$ftr_idx], "-", $FH_HR));
  }
  print $perseq_FH "\n";

  # output actual data, for each sequence
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_accn = $seq_info_HAR->{"accn_name"}[$seq_idx];
    printf $perseq_FH ("%-*s  ", $width_seq, $seq_accn);
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if($ftr_idx > 0) { print $perseq_FH "  "; }
      if($ftr_info_HAR->{"type"}[$ftr_idx] eq "model") { 
        if(! $do_gap_special) { 
          printf $perseq_FH ("%-*s  ", $w_tot_gap_length_A[$ftr_idx], $tot_gap_length_AA[$seq_idx][$ftr_idx]);
          printf $perseq_FH ("%-*s  ", $w_net_gap_length_A[$ftr_idx], $net_gap_length_AA[$seq_idx][$ftr_idx]);
        }
        printf $perseq_FH ("%-*s", $w_gapstr_A[$ftr_idx], $gapstr_AA[$seq_idx][$ftr_idx]);
      }
      else { # feature is not type 'model', we have no gap info
        if(! $do_gap_special) { 
          printf $perseq_FH ("%-*s  ", $w_tot_gap_length_A[$ftr_idx], "0");
          printf $perseq_FH ("%-*s  ", $w_net_gap_length_A[$ftr_idx], "0");
        }
        printf $perseq_FH ("%-*s", $w_gapstr_A[$ftr_idx], "-");
      }
    }
    print $perseq_FH ("\n");
  }

  # output explanatory text
  print $perseq_FH "#\n";
  print $perseq_FH ("# Explanation of the above table:\n");
  if($do_gap_all) { 
    print  $perseq_FH ("# The table includes information on all gaps that exist between all pairwise alignments of\n");
    printf $perseq_FH ("# the reference feature (mature peptide or CDS) and the predicted homologous feature for each sequence.\n");
  }
  elsif($do_gap_not3) { 
    print  $perseq_FH ("# The table includes information on all gaps of lengths that are not multiples of 3 that exist\n");
    printf $perseq_FH ("# between all pairwise alignments of the reference feature (mature peptide or CDS) and the predicted homologous feature for each sequence.\n");
  }
  else { 
    print $perseq_FH ("# The table includes information on some gaps that can solely explain a feature (CDS or mat_peptide) not being a multiple of length 3.\n");
    print $perseq_FH ("# This is (probably) not an exhaustive list of all such gaps.\n");
    print $perseq_FH ("# Specifically it is only gaps X in a feature Y, such that the following criteria are met:\n");
    print $perseq_FH ("#   - length of feature Y is not a multiple of 3\n");
    print $perseq_FH ("#   - if you remove only X from list of all gaps, total feature length of Y becomes a multiple of 3\n");
    print $perseq_FH ("#       with length difference of D with reference feature\n");
    print $perseq_FH ("#    - there are no other gaps Z such that if you remove only Z then length of Y becomes a multiple\n");
    print $perseq_FH ("#       of 3 with length difference D2 from reference where D2 < D.\n");
  }
  print $perseq_FH ("#\n");
  if($do_gap_all || $do_gap_not3) { 
    printf $perseq_FH ("# There are 3 columns under each header \"<%s>#<n> (%s)\" named \"tot\", \"net\",\n", ($do_matpept) ? "mat_peptide" : "CDS", ($do_gap_all) ? "all gaps" : "gaps %3 != 0");
    print $perseq_FH ("# and \"string\".\n");
    print $perseq_FH ("# The \"tot\" columns list the total number of gap positions in either sequence in the pairwise alignment.\n");
    print $perseq_FH ("# The \"net\" columns list the net number of the listed gaps in the pairwise alignment; this is the number\n");
    print $perseq_FH ("#   of gaps in the reference sequence minus the number of gaps in the current sequence (inserts minus deletes)\n");
    print $perseq_FH ("# The \"string\" columns include a list of <n> tokens, each of which describes a gap of length >= 1 nucleotide.\n");
  }
  print  $perseq_FH ("#\n");
  print  $perseq_FH ("# Tokens are in the form: <char><position><length>\n");
  printf $perseq_FH ("#   <char>     is 'I' for an insertion relative to the reference feature (gap in reference sequence)\n");
  printf $perseq_FH ("#              or 'D' for a  deletion  relative to the reference feature (gap in current sequence)\n");
  print  $perseq_FH ("#   <position> is the nucleotide position of the gap in reference coordinates.\n");
  print  $perseq_FH ("#              For insertions this is the reference position after which the insertion occurs.\n");
  print  $perseq_FH ("#              For deletions  this is the first reference position for this deletion.\n");
  print  $perseq_FH ("#   <length>   length of the gap in nucleotides.\n");
  print  $perseq_FH ("#              For insertions this is the number of nucleotides inserted relative to the reference\n");
  print  $perseq_FH ("#              For deletions  this is the number of reference positions deleted.\n");
  print  $perseq_FH ("#\n");
  if($do_gap_special) { 
    print  $perseq_FH ("#\n");
    printf $perseq_FH ("# \"-\" tokens indicate the feature is a multiple of length 3\n");
    printf $perseq_FH ("# \"?\" tokens indicate the feature is not a multiple of length 3, but that no gaps that satisfy our criteria exist.\n");
    print  $perseq_FH ("#\n");
  }

  # Now the per gap information
  if($do_gap_all) { 
    printf $pergap_FH ("# Counts of all gaps:\n#\n");
  }
  elsif($do_gap_not3) { 
    printf $pergap_FH ("# Counts of gaps that have length that is not a multiple of 3:\n#\n");
  }
  else { # $do_gap_special == 1
    printf $pergap_FH ("# Counts of gaps that are special (responsible for net gap length that is not a multiple of 3):\n#\n");
    print $perseq_FH ("# Specifically, these are counts of gaps X in a feature Y (CDS or mat_peptide), such that the following criteria are met:\n");
    print $perseq_FH ("#   - length of feature Y is not a multiple of 3\n");
    print $perseq_FH ("#   - if you remove only X from list of all gaps, total feature length of Y becomes a multiple of 3\n");
    print $perseq_FH ("#       with length difference of D with reference feature\n");
    print $perseq_FH ("#    - there are no other gaps Z such that if you remove only Z then length of Y becomes a multiple\n");
    print $perseq_FH ("#       of 3 with length difference D2 from reference where D2 < D.\n");
    print $perseq_FH ("#\n");
  }
  my $nprinted = 0;
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $is_matpept = ($ftr_info_HAR->{"type"}[$ftr_idx] eq "mp") ? 1 : 0;
    if((scalar(keys %{$ftr_gapstr_AH[$ftr_idx]})) > 0) { 
      foreach my $key (sort keys %{$ftr_gapstr_AH[$ftr_idx]}) { 
        printf $pergap_FH ("%s#" . ($ftr_idx+1) . " " . $key . " " . $ftr_gapstr_AH[$ftr_idx]{$key} . "\n", ($is_matpept) ? "mat_peptide" : "CDS");
        $nprinted++;
      }
      printf $pergap_FH ("#\n");
    }
  }
  if($nprinted == 0) { 
    printf $pergap_FH ("# NONE\n");
  }

  return;
}
#################################################################
#
#  Miscellaneous subroutines that don't fit in the above categories:
#    find_inframe_stop()
#    combine_ajb_and_aja_strings()
#    compare_to_genbank_annotation()
#    count_genbank_annotations()
#    translate_feature_sequences()
#    align_hits()
#    update_gap_array()
#    find_special_gap()
#    define_model_and_feature_output_file_names()
#    get_mdl_or_ftr_ofile_info_key()
#    align_protein_sequences()
#    check_for_downstream_stop()
#    create_output_start_and_stop()
#    seq_name_from_msa_seq_name()
#    accn_name_from_seq_name()
#
#################################################################
# Subroutine:  find_inframe_stop()
# Incept:      EPN, Tue Mar  8 13:39:23 2016
#
# Purpose:    Given a subsequence, find the first inframe stop in it
#             and return the final position of it. Return 0 if none exist
#             in the subsequence.
#
# Arguments: 
#  $sequence:     the actual sequence to search in as a string, no newlines
#
# Returns:    Two values:
#             $ret_posn:   final (3rd) position of first in-frame stop in $sequence.
#                          0 if no in-frame stop exists
#             $stop_codon: in-frame stop codon found, either "TAA", "TAG", "TGA", or "TAR" or undef if
#                          none found
#
# Dies:       never
#
#################################################################
sub find_inframe_stop {
  my $sub_name = "find_inframe_stop()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sequence) = @_;

  # convert to uppercase and convert Us to Ts
  $sequence =~ tr/a-z/A-Z/;
  $sequence =~ tr/U/T/;

  # use substr to get each codon and check if it's a stop
  my $posn = 0;
  my $len  = length($sequence);
  my $ret_posn   = undef;
  my $stop_codon = undef;
  while((! defined $ret_posn) && 
        ($posn+2) < $len) { 
    my $codon = substr($sequence, $posn, 3);
    if($codon eq "TAA") { 
      $stop_codon = "TAA";
      $ret_posn   = $posn+3; # add 3 to get coords into 1..L, (adding 2 would keep us in 0..L-1 coord space)
    }
    elsif($codon eq "TAG") { 
      $stop_codon = "TAG";
      $ret_posn   = $posn+3; # add 3 to get coords into 1..L, (adding 2 would keep us in 0..L-1 coord space)
    }
    elsif($codon eq "TGA") { 
      $stop_codon = "TGA";
      $ret_posn   = $posn+3; # add 3 to get coords into 1..L, (adding 2 would keep us in 0..L-1 coord space)
    }
    elsif($codon eq "TAR") { 
      $stop_codon = "TAR";
      $ret_posn   = $posn+3; # add 3 to get coords into 1..L, (adding 2 would keep us in 0..L-1 coord space)
    }
    $posn += 3;
  }
  
  if(! defined $ret_posn) { 
    $ret_posn = 0; 
  }
  
  return ($ret_posn, $stop_codon);
}

#################################################################
# Subroutine: combine_ajb_and_aja_strings()
# Incept:     EPN, Mon Mar 14 15:10:25 2016
#
# Purpose:    Given a string describing adjacencies-before ($ajb_str)
#             and another describing adjacencies-after ($aja_str),
#             combine them into a single string.
#
# Args:
#  $ajb_str:       adjacencies before string
#  $aja_str:       adjacencies after string
#
# Returns: the combine string
#
#################################################################
sub combine_ajb_and_aja_strings { 
  my $sub_name = "combine_ajb_and_aja_strings";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ajb_str, $aja_str) = @_;

  my $comma_or_not = "";
  if(($ajb_str ne "") && ($aja_str ne "")) { 
    $comma_or_not = ",";
  }
  return $ajb_str . $comma_or_not . $aja_str;
}

#################################################################
# Subroutine:  compare_to_genbank_annotation
# Incept:      EPN, Mon Mar 14 08:42:00 2016
#
# Purpose:    For a given sequence/model pair and start..stop boundaries
#             check if any annotation in %{$tbl_HAR} matches start..stop.
#             Return 1 if it does and 0 if it doesn't.
#
# Arguments: 
#  $pred_start:   predicted start (from our annotation, in -accn_len..accn_len coord space)
#  $pred_stop:    predicted stop (from our annotation, in -accn_len..accn_len coord space)
#  $pred_strand:  predicted strand (from our annotation)
#  $accn_len:     length of the accession's sequence we're interested in
#  $seq_len:      length of the sequence in the file we searched in
#  $tbl_HAR:      REF to hash of arrays for accession we're interested in, PRE-FILLED
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:        REF to hash of file handles
#
# Returns:    '1' if we find a match to $pred_start..$pred_stop to existing (GenBank)
#             annotation stored in %{$tbl_HAR}, else '0'
#
# Dies: If we have a problem getting information we need from $tbl_HAR
#
################################################################# 
sub compare_to_genbank_annotation { 
  my $sub_name = "compare_to_genbank_annotation()";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($pred_start, $pred_stop, $pred_strand, $accn_len, $seq_len, $tbl_HAR, $opt_HHR, $FH_HR) = @_;

  # printf("in $sub_name, pred: $pred_start..$pred_stop $pred_strand\n");

  # get coordinates into the following format:
  # start <= stop
  # coordinate space: -accn_len..accn_len (should already be in this)
  # check we're in -accn_len..accn_len space
  if(($pred_start < ($accn_len * -1)) || ($pred_start > $accn_len)) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, predicted start coordinate: $pred_start not in expected -accn_len..accn_len coordinate space (%d..%d)", (-1 * $accn_len), $accn_len), 1, $FH_HR); 
  }
  if(($pred_stop < ($accn_len * -1)) || ($pred_stop > $accn_len)) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, predicted stop coordinate: $pred_stop not in expected -accn_len..accn_len coordinate space (%d..%d)", (-1 * $accn_len), $accn_len), 1, $FH_HR); 
  }
  if($pred_start > $pred_stop) { 
    my $tmp     = $pred_start;
    $pred_start = $pred_stop;
    $pred_stop  = $tmp;
  }

  my $found_match = 0; # set to '1' if we find a match
  my @coords_A    = (); # coordinates of each annotated feature
  my @len_A       = (); # lengths of each annotated feature (not used)
  # fill @coords_A and len_A
  getLengthsAndCoords($tbl_HAR, \@len_A, \@coords_A, $FH_HR);

  # for each annotated feature, check if we have a match to $pred_start..$pred_stop
  my $nannot = scalar(@len_A);
  for(my $i = 0; $i < $nannot; $i++) { 
    my @starts_A  = ();
    my @stops_A   = ();
    my @strands_A = ();
    my $nsegments = 0;
    startsStopsStrandsFromCoordsLength($coords_A[$i], $accn_len, $opt_HHR, \@starts_A, \@stops_A, \@strands_A, \$nsegments, $FH_HR);
    for(my $j = 0; $j < $nsegments; $j++) { 
      # printf("\tcomparing against GB $starts_A[$j]..$stops_A[$j] $strands_A[$j]\n");
      # get coordinates into the following format:
      # start <= stop
      # coordinate space: -accn_len..accn_len
      my $cur_start = $starts_A[$j];
      my $cur_stop  = $stops_A[$j];
      if($cur_start > $cur_stop) { 
        # swap them (GenBank doesn't follow same order convention that we do, so we always set start < stop)
        my $tmp    = $cur_start;
        $cur_start = $cur_stop;
        $cur_stop  = $tmp;
      }
      ($cur_start, $cur_stop) = create_output_start_and_stop($cur_start, $cur_stop, $accn_len, $seq_len, $FH_HR);
      if(($cur_start == $pred_start) && ($cur_stop == $pred_stop) && ($strands_A[$j] eq $pred_strand)) { 
        $found_match = 1;
        $j = $nsegments+1; # breaks j loop
        $i = $nannot+1;    # breaks i loop
      }
    }
  }

  # printf("returning $found_match\n");
  return $found_match;
}

#################################################################
# Subroutine:  count_genbank_annotations
# Incept:      EPN, Wed Mar 16 08:40:47 2016
#
# Purpose:    Return the number of feature (e.g. CDS) and segment (e.g. exon)
#             annotations in tbl_HAR.
#
# Arguments: 
#  $tbl_HAR:      REF to hash of arrays for accession we're interested in, PRE-FILLED
#  $accn_len:     length of the accession's sequence we're interested in
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:        REF to hash of file handles
#
# Returns:    Two values:
#             1) number of total feature (e.g. CDS) annotations in $tbl_HAR
#             2) number of total segment (e.g. exon) annotations in $tbl_HAR
#
# Dies: If we have a problem getting information we need from $tbl_HAR
#
################################################################# 
sub count_genbank_annotations { 
  my $sub_name = "count_genbank_annotations()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($tbl_HAR, $accn_len, $opt_HHR, $FH_HR) = @_;

  if(! defined $tbl_HAR) { 
    # no annotations
    return (0, 0);
  }
  
  my @coords_A    = (); # coordinates of each annotated feature
  my @len_A       = (); # lengths of each annotated feature (not used)

  # fill @coords_A and len_A
  getLengthsAndCoords($tbl_HAR, \@len_A, \@coords_A, $FH_HR);
  
  my $nannot = scalar(@len_A);
  my $tot_nsegments = 0;
  for(my $i = 0; $i < $nannot; $i++) { 
    my @starts_A  = ();
    my @stops_A   = ();
    my @strands_A = ();
    my $nsegments = 0;
    startsStopsStrandsFromCoordsLength($coords_A[$i], $accn_len, $opt_HHR, \@starts_A, \@stops_A, \@strands_A, \$nsegments, $FH_HR);
    $tot_nsegments += $nsegments;
  }

  return ($nannot, $tot_nsegments);
}

#################################################################
# Subroutine:  translate_feature_sequences
# Incept:      EPN, Thu Mar 10 11:13:25 2016
#
# Purpose:    For each file with corrected feature sequences, 
#             translate them into proteins.
#
# Arguments: 
#  $in_key:         key for the input files we'll translate here, usually "corrected"
#  $out_key:        key for the output files we'll create here, usually "corrected.translated",
#                   this key will be stored in $ftr_info_HAR
#  $specstart_AAR:  REF to 2D array of specified, permissible start codons, can be undef
#  $ftr_info_HAR:   REF to hash of arrays with information on the features, ADDED TO HERE
#  $ofile_info_HHR: REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
################################################################# 
sub translate_feature_sequences { 
  my $sub_name = "translate_feature_sequences()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($in_key, $out_key, $specstart_AAR, $ftr_info_HAR, $ofile_info_HHR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nftr: number of features

  my $in_ftr_info_file_key     = $in_key  . ".hits.fa";
  my $out_ftr_info_file_key    = $out_key . ".hits.fa";

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $nucleotide_fafile = $ftr_info_HAR->{$in_ftr_info_file_key}[$ftr_idx];
    my $protein_fafile    = $ftr_info_HAR->{$out_ftr_info_file_key}[$ftr_idx];

    if($nucleotide_fafile ne "/dev/null") { # if this is set to /dev/null we know it's not supposed to exist, so we skip this feature
      my $opts = "";
      # require a proper start codon and stop codon if this is not a mature peptide
      if($ftr_info_HA{"type"}[$ftr_idx] ne "mp") { 
        $opts = " -reqstart -reqstop ";
      }
      my $altstart_opt = get_esl_epn_translate_altstart_opt($ftr_info_HAR, $ftr_idx, $specstart_AAR);
      
      # use esl-epn-translate.pl to examine the start and stop codons in each feature sequence
      $cmd = $esl_epn_translate . " -endatstop -nostop $opts $altstart_opt $nucleotide_fafile > $protein_fafile";
      runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
      
      # determine the number of >= 1 segments (exons or mature peptides) we put together to make this protein
      my $nsegments = 0;
      if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { 
        $nsegments = $ftr_info_HA{"final_mdl"}[$ftr_idx] - $ftr_info_HA{"first_mdl"}[$ftr_idx] + 1;
      }
      elsif($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") { 
        my @primary_children_idx_A = ();
        getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "primary", \@primary_children_idx_A, $ofile_info_HHR->{"FH"});
        $nsegments = scalar(@primary_children_idx_A);
      }
      else { 
        DNAORG_FAIL("ERROR in $sub_name, feature $ftr_idx with name %s is of unknown annot_type %s\n", $ftr_info_HAR->{"out_tiny"}[$ftr_idx], $ftr_info_HAR->{"annot_type"}[$ftr_idx]);
      } 
      
      my $ofile_info_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $out_ftr_info_file_key, $ofile_info_HHR->{"FH"});
      addClosedFileToOutputInfo(\%ofile_info_HH, $ofile_info_key, $protein_fafile, 0, sprintf("fasta file with translations of corrected hits for feature " . $ftr_info_HA{"out_tiny"}[$ftr_idx] . " composed of %d segments", $nsegments));
    }
  }

  return;
}

#################################################################
# Subroutine:  align_hits()
# Incept:      EPN, Thu Mar 10 13:12:26 2016
#
# Purpose:    For each model, align its corrected hits to the 
#             model. Save information on where gaps are relative
#             to the reference \%{$ref_del_HHA} and \%{$ref_ins_HHA}.
#
# Arguments: 
#  $execs_HR:         REF to a hash with "cmalign" executable paths
#  $mdl_info_HAR:     REF to the hash of arrays with model information
#  $seq_info_HAR:     REF to the hash of arrays with sequence information
#  $mdl_results_AAHR: REF to results AAH, ADDED TO HERE
#  $opt_HHR:          REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:   REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
################################################################# 
sub align_hits {
  my $sub_name = "align_hits";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $mdl_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $ofile_info_HHR->{"FH"}); # nseq: number of sequences

  # get index hash for @{$seq_info_HAR->{"seq_name"}} array
  # this simplifies determining sequence index in @{%seq_info_HAR->{}}
  # arrays for a given sequence name.
  my %seq_name_idx_H = (); # key: $seq_name, value: idx of $seq_name in @{$seq_info_HAR->{"seq_name"}}
  getIndexHashForArray(\@{$seq_info_HAR->{"seq_name"}}, \%seq_name_idx_H, $ofile_info_HHR->{"FH"});

  my $mdl_fa_file_key   = "corrected.hits.fa";
  my $mdl_stk_file_key  = "corrected.hits.stk";

  for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    # validate that we have the fasta file we need
    my $fa_file  = $mdl_info_HAR->{$mdl_fa_file_key}[$mdl_idx];
    validateFileExistsAndIsNonEmpty($fa_file, $sub_name, $ofile_info_HHR->{"FH"});
    my $stk_file = $mdl_info_HAR->{$mdl_stk_file_key}[$mdl_idx];  # Stockholm alignment file model was built from
    my $cmname   = $mdl_info_HAR->{"cmname"}[$mdl_idx];           # name of model 
    my $cm_file  = $mdl_info_HAR->{"cmfile"}[$mdl_idx];           # file name
    validateFileExistsAndIsNonEmpty($cm_file, $sub_name, $ofile_info_HHR->{"FH"});
    
    # create the alignment
    my $mxsize_opt = sprintf("--mxsize %d", opt_Get("--mxsize", $opt_HHR));
    my $cmd = $execs_HR->{"cmalign"} . " --cpu 0 $mxsize_opt $cm_file $fa_file > $stk_file";
    runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HHR->{"FH"});
    # save this file to %{$ofile_info_HHR}
    my $ofile_key = get_mdl_or_ftr_ofile_info_key("mdl", $mdl_idx, $mdl_stk_file_key, $ofile_info_HHR->{"FH"});
    addClosedFileToOutputInfo($ofile_info_HHR, $ofile_key, $stk_file, 0, sprintf("Stockholm alignment of hits for model #%d: %s", $mdl_idx+1, $mdl_info_HAR->{"out_tiny"}[$mdl_idx]));

    # use the alignment to determine:
    # - fractional identity between prediction and reference 
    #   and store in $mdl_results_AAH->[][]{"fid2ref"}
    # - where deletions are relative to reference and 
    #   store in %{$refdel_HHAR}
    # - where insertions are relative to reference and 
    #   store in %{$refins_HHAR}

    # first we need to read in the MSA we just created 
    my $msa = Bio::Easel::MSA->new({
      fileLocation => $stk_file,
      isDna => 1
                                   });  

    # Determine which positions are RF (reference) positions,
    # these are the consensus (nongap) positions of the model.
    # In our case these are positions in the reference sequence
    # the model was built from (as opposed to gaps in the
    # reference sequence, which are non-RF positions).
    my $rfseq = $msa->get_rf();
    my @rfseq_A = split("", $rfseq);
    my $alen  = $msa->alen();    # alignment length
    my @i_am_rfpos_A = (); # 0..$apos..$alen-1: '1' if RF is a nongap at position $apos+1, else '0'
    for(my $apos = 0; $apos < $alen; $apos++) { 
      $i_am_rfpos_A[$apos] = ($rfseq_A[$apos] eq "." || $rfseq_A[$apos] eq "~") ? 0 : 1;
    }          
    
    my $ref_seq_idx = 0; # this will remain '0'
    my $msa_nseq = $msa->nseq;
    for(my $msa_seq_idx = 0; $msa_seq_idx < $msa_nseq; $msa_seq_idx++) { 
      my $msa_seq_name = $msa->get_sqname($msa_seq_idx);
      my $seq_name = seq_name_from_msa_seq_name($msa_seq_name, $ofile_info_HHR->{"FH"}); 
      if(! exists $seq_name_idx_H{$seq_name}) { 
        DNAORG_FAIL("ERROR in $sub_name, sequence $seq_name derived from MSA seq name: $msa_seq_name does not exist in seq_info_HAR", 1, $ofile_info_HHR->{"FH"});
      }
      my $seq_idx = $seq_name_idx_H{$seq_name};
      if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"}) { 
        $mdl_results_AAH[$mdl_idx][$seq_idx]{"fid2ref"} = $msa->pairwise_identity($ref_seq_idx, $msa_seq_idx);
        # printf("storing percent id of $fid2ref_HHR->{$mdl}{$seq} for $mdl $seq\n"); 

        # determine the RF positions that are gaps in this sequence
        # and the positions of the inserted residues in this sequence
        my @tmp_refdel_A = (); # array of deletions, temporary because we'll convert it to a string at end of function
        my @tmp_refins_A = (); # array of deletions, temporary because we'll convert it to a string at end of function
        my $aseqstring  = $msa->get_sqstring_aligned($msa_seq_idx);
        my @aseq_A = split("", $aseqstring);
        my $rfpos = 0;
        for(my $apos = 0; $apos < $alen; $apos++) { 
          if($i_am_rfpos_A[$apos]) { 
            $rfpos++; 
          }
          if($aseq_A[$apos] =~ m/[\.\-]/) { # a gap in the sequence
            if($i_am_rfpos_A[$apos]) { # not a gap in the RF sequence
              # deletion (gap) relative to the reference sequence
              update_gap_array(\@tmp_refdel_A, $rfpos, 1); # 1 informs the subroutine that this is a delete array
            }
          }
          else { # nongap in the sequence
            if(! $i_am_rfpos_A[$apos]) { # gap in the RF sequence
              # insertion in sequence relative to the reference sequence
              update_gap_array(\@tmp_refins_A, $rfpos, 0); # 0 informs the subroutine that this is an insert array
            }
          }
        }
        # printf("printing insert info for $mdl $seq\n");
        # debugPrintGapArray(\@tmp_refins_A);
        # printf("printing delete info for $mdl $seq\n");
        # debugPrintGapArray(\@tmp_refdel_A);

        # convert arrays of refdel and refins to strings, to be stored in results_AAHR
        my $el; # one array element

        # refdel
        $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refdelstr"} = "";
        #printf("HEYA INIT'ED mdl_results_AAHR->[$mdl_idx][$seq_idx]{refdelstr} to empty string\n");
        foreach $el (@tmp_refdel_A) { 
          $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refdelstr"} .= $el . ",";
        }
        $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refdelstr"} =~ s/\,$//; # remove final comma

        # refins
        $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refinsstr"} = "";
        #printf("HEYA INIT'ED mdl_results_AAHR->[$mdl_idx][$seq_idx]{refinsstr} to empty string\n");
        foreach $el (@tmp_refins_A) { 
          $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refinsstr"} .= $el . ",";
        }
        $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"refinsstr"} =~ s/\,$//; # remove final comma
      }
    }
  }

  return;
}

#################################################################
# Subroutine:  update_gap_array()
# Incept:      EPN, Thu Mar 10 13:42:31 2016
#
# Purpose:    Given an rfpos that we have a gap in, update a 
#             array that stores all the gap information for a given
#             sequence. This can be called for an array that holds
#             gaps in the reference (deletions, refdel* data 
#             structures) or in the other sequence (insertions,
#             refins* data structures).
#
#
# Arguments: 
#  $AR:         REF to the array to update
#  $rfpos:      the reference position the gap occurs at or after
#  $is_delete:  '1' if we're updating a delete array, else we're
#               updating an insert array (we do the udpate slightly
#               differently for each type).
#
# Returns:    void
#
################################################################# 
sub update_gap_array {
  my $sub_name = "update_gap_array";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($AR, $rfpos, $is_delete) = @_;

  my $nel = scalar(@{$AR});  # number of elements
  my $same_as_prv = 0;       # insertion or deletion occurs at a position detected previously
  if($nel > 0) { 
    # need to check latest element to see if it pertains to the same rfpos
    my ($prv_rfpos, $prv_cnt) = split(":", $AR->[$nel-1]);
    # printf("\tsplit %s into %s and %d\n", $AR->[$nel-1], $prv_rfpos, $prv_cnt);
    if((  $is_delete) && (($prv_rfpos + $prv_cnt) == $rfpos) ||
       (! $is_delete) && ($prv_rfpos == $rfpos)) { 
      # this is not the first insert/delete at this position, update previous value
      $same_as_prv = 1;
      $AR->[($nel-1)] = $prv_rfpos . ":" . ($prv_cnt+1);
    }
  }
  if(! $same_as_prv) { 
    # either this is the first insert/delete or not the same as the previous one
    push(@{$AR}, $rfpos . ":" . "1");
  }

  return;
}

#################################################################
# Subroutine: find_special_gap()
#
# Purpose: Given a gap string that lists all gaps of a feature in an
#          alignment, determine if the gaps cause the predicted
#          length of the feature to be non-modulo 3. If so,
#          determine if there's exactly one gap that if we remove
#          it, the remaining gaps will make the predicted length of
#          the CDS to be modulo 3.
#
#          I think there can only be 0 or 1 such gap, but I'm not sure.
# 
#          Return two values: 
#          - the length of the speical gap if there is one, or 0 if there is not
#            (this value will be negative if it is a deletion relative to the 
#            reference, and positive if it is an insertion relative to the reference)
#          - the string that describes the special gap, 
#            or '-' if the predicted length of the feature is modulo 3 if we include all gaps
#            or '?' if the predicted length of the feature is not modulo 3 but there is no
#            special gap that can explain the non-modulo-3ness.
#
#          If the net gap length is zero modulo 3, return (0, "-").
#          If the net gap length is nonzero modulo 3, and there's exactly 1 gap that equals the net length modulo 3, return its length and the string..
#          If the net gap length is nonzero modulo 3, and there's not exactly 1 gap that equals the net length modulo 3, return (0, "?");
#            
#          Examples: 
#
#            Input:    $gapstr = "D377:2,I388:2,I1066:1 net gap is 1
#            Returns:  (1, "I1066:1")
#
#            Input:    $gapstr = "D377:2,I388:2", , net gap is 0;
#            Returns:  (0, "-")
#
#            Input:    $gapstr = "D377:5,I388:2", net gap is -3;
#            Returns:  (0, "-")
#
#            Input:    $gapstr = "D377:2,I388:5", net gap is 3;
#            Returns:  (0, "-")
#
#            Input:    $gapstr = "D377:2,I388:3,I1066:1", net gap is 2;
#            Returns:  (0, "?")
#
#            Input:    $gapstr = "D277:2,D377:2,I388:2", $net_gap = 2;
#            Returns:  (0, "?")
#
# Arguments: 
#  $gapstr:  string describing all gaps
#
# Returns:   Two values: as explained above in "Purpose"
# 
################################################################# 
sub find_special_gap {
  my $sub_name = "find_special_gap";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($gapstr) = (@_);
  # printf("in $sub_name gapstr: $gapstr\n");

  if($gapstr eq "-") { # no gaps
    return ("-"); 
  }

  # determine net number of gaps if we include them all
  my @el_A = split(",", $gapstr);
  my $ngaps = scalar(@el_A);
  my @gapsize_A = (); # [0..$ngaps-1], length of each gap
  my $tot_gaps = 0;

  for(my $e = 0; $e < $ngaps; $e++) { 
    my $el = $el_A[$e];
    my ($type_loc, $len) = split(":", $el);

    if(! defined $len) { printf("gapstr: $gapstr, el: $el\n"); }

    if(($len % 3) == 0) { 
      $gapsize_A[$e] = 0;
    }
    else { 
      if($type_loc =~ /^I(\d+)/) { # insert
        $gapsize_A[$e] = $len;
        $tot_gaps += $len;
      }
      elsif($type_loc =~ /^D(\d+)/) { # delete
        $gapsize_A[$e] = -1 * $len;
        $tot_gaps -= $len;
      }
      else { 
        die "ERROR: in $sub_name, unable to parse gap string: $gapstr element $el"; 
      }
    }
    if(($tot_gaps % 3) == 0) { 
      return ("-"); 
    }
  }

  # find gap that gives minimal net gap length
  
  my $min_len = undef;
  for(my $e = 0; $e < $ngaps; $e++) { 
    my $net_gap_len = $tot_gaps - $gapsize_A[$e];
    if(($net_gap_len % 3) == 0) { 
      if(! defined $min_len || ($net_gap_len < $min_len)) { 
        $min_len = $net_gap_len;
      }
    }
  }

  my $nspecial = 0;
  my $special_str = "";
  my $special_len = 0;

  if(! defined $min_len) { # no gaps can be removed and give modulo 3
    return ("?");
  }
  else { # at least 1 gap can be removed and give modulo 3
    for(my $e = 0; $e < $ngaps; $e++) { 
      my $net_gap_len = $tot_gaps - $gapsize_A[$e];
      if(($net_gap_len % 3) == 0) { 
        if($net_gap_len == $min_len) { 
          if($special_str ne "") { $special_str .= ","; }
          $special_str .= $el_A[$e];
        }
      }
    }
    # printf("returning special_str: $special_str\n");
    return $special_str;
  }
}

#################################################################
# Subroutine:  define_model_and_feature_output_file_names()
# Incept:      EPN, Thu Mar 10 13:47:15 2016
#
# Purpose:    Define the output file names for all per-model
#             and per-feature files we will create later in the script.
#
# Arguments: 
#  $out_root:     output root for the file names
#  $mdl_info_HAR: REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR: REF to hash of arrays with information on the features, PRE-FILLED
#  $FH_HR:        REF to hash of file handles
#
# Returns:    void
#
# Dies: if something is wrong with $mdl_info_HAR
#
################################################################# 
sub define_model_and_feature_output_file_names { 
  my $sub_name = "define_model_and_feature_output_file_names";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_root, $mdl_info_HAR, $ftr_info_HAR, $FH_HR) = @_;

  my $nmdl = validateModelInfoHashIsComplete  ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $FH_HR); # nftr: number of features

  my @both_file_types_A = ("predicted.hits.fa", "predicted.hits.append.fa", "corrected.hits.fa", "corrected.hits.append.fa");
  my @mdl_only_file_types_A  = ("corrected.hits.stk");
  my @ftr_only_file_types_A  = ("predicted.hits.fa.esl-epn-translate", "corrected.translated.hits.fa", 
                                "corrected.translated.hits.stk", "corrected.translated.hmm", 
                                "corrected.translated.hmmbuild", "corrected.translated.hmmstk");

  for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    foreach my $file_type (@both_file_types_A, @mdl_only_file_types_A) { 
      $mdl_info_HAR->{$file_type}[$mdl_idx] = $out_root . "." . $mdl_info_HAR->{"filename_root"}[$mdl_idx] . "." . $file_type;
    }
  }

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    foreach my $file_type (@both_file_types_A, @ftr_only_file_types_A) { 
      $ftr_info_HAR->{$file_type}[$ftr_idx] = $out_root . "." . $ftr_info_HAR->{"filename_root"}[$ftr_idx] . "." . $file_type;
    }
  }

  return;
}

#################################################################
# Subroutine:  get_mdl_or_ftr_ofile_info_key()
# Incept:      EPN, Thu Mar 10 13:58:48 2016
#
# Purpose:    Simple function for defining an ofile_info_HH
#             key for a per-model of per-feature file given the
#             index and key in %mdl_info_HA or %ftr_info_HA.
#             Return the ofile_info_HH key.
#
# Arguments: 
#  $mdl_or_ftr: output root for the file names
#  $idx:        index, becomes part of returned string
#  $key:        key in %{$mdl_info_HAR} or %{$ftr_info_HAR}
#               becomes part of the key returned by this function
#  $FH_HR:      REF to hash of file handles, used only if 
#               we need to die (so that open file handles get
#               closed before we exit and error gets printed to the
#               log file)
#
# Returns:    the key for %ofile_info_HH
#
# Dies: if $mdl_or_ftr is not 'mdl' nor 'ftr'
#       if $idx is not a positive integer
################################################################# 
sub get_mdl_or_ftr_ofile_info_key { 
  my $sub_name = "get_mdl_or_ftr_ofile_info_key";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_or_ftr, $idx, $key, $FH_HR) = @_;

  if(($mdl_or_ftr ne "mdl") && ($mdl_or_ftr ne "ftr")) { 
    DNAORG_FAIL("ERROR in $sub_name, mdl_or_ftr should be mdl or ftr, got $mdl_or_ftr", 1, $FH_HR);
  }
  if((! verify_integer($idx)) || ($idx < 0)) { 
    DNAORG_FAIL("ERROR in $sub_name, expected non-negative integer for idx, got $idx", 1, $FH_HR);
  }

  return $mdl_or_ftr . "." . $idx . "." . $key;
}


#################################################################
# Subroutine:  align_protein_sequences
# Incept:      EPN, Thu Mar 10 15:52:44 2016
#
# Purpose:    For each protein sequence feature, fetch the reference
#             sequence, build an HMM from it and align all other
#             proteins to it.
#
# Arguments: 
#  $execs_HR:       REF to a hash with "hmmbuild" and "hmmalign"
#                   executable paths
#  $out_key:        key for the output files we'll create here, usually "corrected.translated",
#                   this key will be stored in $ftr_info_HAR
#  $ftr_info_HAR:   REF to hash of arrays with information on the features, ADDED TO HERE
#  $opt_HHR:        REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR: REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
################################################################# 
sub align_protein_sequences { 
  my $sub_name = "align_protein_sequences()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $out_key, $ftr_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nftr: number of features

  # keys for $ftr_info_HAR, which holds the file names we will use here
  # these were defined in define_model_and_feature_output_file_names()
  my $ftr_info_fa_file_key       = $out_key . ".hits.fa";  # fasta file of all annotated protein sequences, one per feature, already created
  my $ftr_info_hmmstk_file_key   = $out_key . ".hmmstk";   # single sequence alignment of just the reference sequence, one per feature, created here
  my $ftr_info_hmm_file_key      = $out_key . ".hmm";      # HMM built from .hmmstk file, one per feature, created here
  my $ftr_info_hmmbuild_file_key = $out_key . ".hmmbuild"; # hmmbuild output when HMM was built from .hmmstk file, created here
  my $ftr_info_stk_file_key      = $out_key . ".hits.stk"; # all annotated protein sequences, aligned, one per feature, created here by 
                                                           # aligning .hits.fa file with .hmm HMM file

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $fa_file       = $ftr_info_HAR->{$ftr_info_fa_file_key}[$ftr_idx];
    my $stk_file      = $ftr_info_HAR->{$ftr_info_stk_file_key}[$ftr_idx];
    my $hmm_file      = $ftr_info_HAR->{$ftr_info_hmm_file_key}[$ftr_idx];
    my $hmmbuild_file = $ftr_info_HAR->{$ftr_info_hmmbuild_file_key}[$ftr_idx];
    my $hmmstk_file   = $ftr_info_HAR->{$ftr_info_hmmstk_file_key}[$ftr_idx];

    # fetch the reference sequence only 
    # first remove any old .ssi files that may exist
    my $ssi_file = $fa_file . ".ssi";
    if(-e $ssi_file) { 
      removeFileUsingSystemRm($ssi_file, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
    }
    if(! -e $fa_file) { 
      DNAORG_FAIL("ERROR in $sub_name, input protein fasta file to align $fa_file does not exist", 1, $ofile_info_HHR->{"FH"});
    }
    # it's possible that the file exists, but is empty, if none of the predicted sequences are translatable
    if(-s $fa_file) {  
      my $prot_sqfile  = Bio::Easel::SqFile->new({ fileLocation => $fa_file });
      my $ref_prot_str = $prot_sqfile->fetch_consecutive_seqs(1, "", -1, undef); # the reference protein sequence string
      my ($ref_prot_name, $ref_prot_seq) = split(/\n/, $ref_prot_str);
      $ref_prot_name =~ s/^\>//; # remove fasta header line '>'
      
      # write it out to a new stockholm alignment file
      open(OUT, ">", $hmmstk_file) || die "ERROR, unable to open $hmmstk_file for writing.";
      print OUT ("# STOCKHOLM 1.0\n");
      print OUT ("$ref_prot_name $ref_prot_seq\n");
      print OUT ("//\n");
      close OUT;
      
      # build an HMM from this single sequence alignment:
      my $cmd = $execs_HR->{"hmmbuild"} . " $hmm_file $hmmstk_file > $hmmbuild_file";
      runCommand($cmd, opt_Get("-v", $opt_HHR), $ofile_info_HHR->{"FH"}); 
      
      # store the files we just created
      my $ofile_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $ftr_info_hmm_file_key, $ofile_info_HHR->{"FH"});
      addClosedFileToOutputInfo($ofile_info_HHR, $ofile_key, $hmm_file, 0, sprintf("HMM built from the reference protein sequence for feature #%d: %s", $ftr_idx+1, $ftr_info_HAR->{"out_tiny"}[$ftr_idx]));
      
      # remove the hmmbuild and hmmbuild alignment files, unless --keep
      if(! opt_Get("--keep", $opt_HHR)) { 
        removeFileUsingSystemRm($hmmstk_file,   $sub_name, $opt_HHR, $ofile_info_HHR);
        removeFileUsingSystemRm($hmmbuild_file, $sub_name, $opt_HHR, $ofile_info_HHR);
      }
      else { 
        $ofile_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $ftr_info_hmmstk_file_key, $ofile_info_HHR->{"FH"});
        addClosedFileToOutputInfo($ofile_info_HHR, $ofile_key, $hmmstk_file, 0, sprintf("Stockholm alignment of reference protein for feature #%d: %s", $ftr_idx+1, $ftr_info_HAR->{"out_tiny"}[$ftr_idx]));
        $ofile_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $ftr_info_hmmbuild_file_key, $ofile_info_HHR->{"FH"});
        addClosedFileToOutputInfo($ofile_info_HHR, $ofile_key, $hmmbuild_file, 0, sprintf("hmmbuild output file for feature #%d: %s", $ftr_idx+1, $ftr_info_HAR->{"out_tiny"}[$ftr_idx]));
      }
      
      # align all sequences to this HMM
      $cmd = $execs_HR->{"hmmalign"} . " $hmm_file $fa_file > $stk_file";
      runCommand($cmd, opt_Get("-v", $opt_HHR), $ofile_info_HHR->{"FH"});
      
      # store the file in ofile_info_HH
      $ofile_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $ftr_info_stk_file_key, $ofile_info_HHR->{"FH"});
      addClosedFileToOutputInfo($ofile_info_HHR, $ofile_key, $stk_file, 0, sprintf("alignment of all protein sequences to reference for feature #%d: %s", $ftr_idx+1, $ftr_info_HAR->{"out_tiny"}[$ftr_idx]));
      
      # remove the .ssi files, always
      if(-e $fa_file . ".ssi") { 
        removeFileUsingSystemRm($fa_file . ".ssi", $sub_name, $opt_HHR, $ofile_info_HHR);
      }
    }
  }

  return;
}

#################################################################
# Subroutine:  check_for_downstream_stop()
# Incept:      EPN, Tue Mar  8 13:15:47 2016
#
# Purpose:    Check for an in-frame stop codon in the subsequence
#             starting at position $posn on strand $strand in sequence
#             $seq_accn in open sequence file $sqfile.
#
# Arguments: 
#  $sqfile:       open Bio::Easel::SqFile object to fetch sequences from
#  $seq_name:     sequence name to look for stop in
#  $posn:         position to start looking for stop in
#  $maxdist:      total number of nucleotides to look for in-frame stop
#                 $posn+$maxdist is the final nucleotide we will examine
#                 when looking for a stop
#  $strand:       strand on which to look for stop codon
#  $FH_HR:        REF to hash of file handles
#
# Returns:    Two values:
#             $ret_posn:   If we find an in-frame stop before the end of the sequence:
#                            final (3rd) position of first in-frame stop in $sequence
#                            *RELATIVE TO $posn*, e.g. if this value is '3' then the first
#                            3 nt of the sequence starting at $posn is a stop codon.
#                            This value will be <0 if strand is "-"
#                          If we don't find an in-frame stop before the end of the sequence:
#                            this value is 0
#             $stop_codon: in-frame stop codon found, either "TAA", "TAG", "TGA" or "TAR" or undef if
#                          none found
#
# Dies:       if unable to fetch sequence from $sqfile
#
#################################################################
sub check_for_downstream_stop { 
  my $sub_name = "check_for_downstream_stop()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $seq_name, $posn, $maxdist, $strand, $FH_HR) = @_;

  my $seqlen = $sqfile->fetch_seq_length_given_name($seq_name);

  # determine final position we'll look at
  my $final_posn = $posn + $maxdist;
  if($final_posn > $seqlen) { 
    $final_posn = $seqlen; 
  }

  my $begin_posn = $posn;
  my $target_subseq_len = 600; # we'll look at this many nucleotides at a time
  my $inframe_stop_posn = 0;
  my $stop_codon = undef; # type of stop codon found, e.g. "TAA", "TAG", "TGA", or "TAR"
  my $keep_going = 1; # set to 0 if we find a stop codon or we go past $final_posn
  while($keep_going) { 
    my @fetch_AA = (); # the 2D array we use to fetch sequence chunks
    @{$fetch_AA[0]} = ();
    my $end_posn;
    if($strand eq "+") { 
      $end_posn = $begin_posn + ($target_subseq_len - 1);
      if($end_posn > $final_posn) { 
        $end_posn = $final_posn; 
        $keep_going = 0; # breaks loop
      }
    }
    else { 
      $end_posn = $begin_posn - ($target_subseq_len - 1);
      if($end_posn < 1) { 
        $end_posn = 1;
        $keep_going = 0; # breaks loop
      }
    }
    @{$fetch_AA[0]} = ($seq_name . "/" . $begin_posn . "-" . $end_posn, $begin_posn, $end_posn, $seq_name);
    my (undef, $cur_subseq) = split(/\n/, $sqfile->fetch_subseqs(\@fetch_AA, -1, undef)); 
    # look for in-frame stop in $cur_subseq sequence fragment
    ($inframe_stop_posn, $stop_codon) = find_inframe_stop($cur_subseq);

    if($inframe_stop_posn != 0) { 
      $keep_going = 0; # breaks loop
    }
    else { # only update these if we didn't find an in-frame stop
           # if we did find an in-frame stop, then we don't want to update
           # these so we can infer the stop posn below
      $begin_posn = $end_posn;
      if($strand eq "+") { $begin_posn++; }
      else               { $begin_posn--; }
    }
  } # end of while($keep_going)

  if($inframe_stop_posn != 0) {
    # we need to return the stop position relative to $posn that was passed in
    if($strand eq "+") { 
      $inframe_stop_posn += ($begin_posn - 1);
      $inframe_stop_posn -= $posn;
    }
    else { # strand eq "-" 
      $inframe_stop_posn = ($begin_posn - $inframe_stop_posn) + 1;
      $inframe_stop_posn = $posn - $inframe_stop_posn;
      # inframe_stop_posn will be negative
    }
  }

  return($inframe_stop_posn, $stop_codon); # if $inframe_stop_posn is 0, we didn't find an in-frame stop codon, and stop codon will be undef
}

#################################################################
# Subroutine:  create_output_start_and_stop
# Incept:      EPN, Mon Mar 14 09:08:48 2016
#
# Purpose:    Convert predicted start and stop coordinates which are
#             in coordinate space 1..seq_len to 'output' coordinates 
#             in coorindate space 1..accn_len. If seq_len == accn_len,
#             which will be the case if the -c option was not used (that is,
#             if the genome is not circular and thus we did not duplicate it)
#             we do nothing. If seq_len == accn_len * 2, then we have to 
#             check if out_start and out_stop should be different from p_start
#             and p_stop.
#
# Arguments: 
#  $in_start:   predicted start (from our annotation)
#  $in_stop:    predicted (possibly corrected) stop (from our annotation)
#  $accn_len:   length of the GenBank accession's sequence
#  $seq_len:    length of the sequence we searched, either ($accn_len) or (2 * $accn_len)
#  $FH_HR:      REF to hash of file handles
#
# Returns:    Two values:
#             $out_start: the start coordinate to output
#             $out_stop:  the stop coordinate to output
#
# Dies: If (($seq_len != $accn_len) && ($seq_len != (2*$accn_len))
#
################################################################# 
sub create_output_start_and_stop { 
  my $sub_name = "create_output_start_and_stop()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($in_start, $in_stop, $accn_len, $seq_len, $FH_HR) = @_;
  
  # contract check
  if((! verify_integer($in_start)) || ($in_start <= 0)) { 
    DNAORG_FAIL("ERROR in $sub_name, input start coordinate is expected to be a positive integer, got $in_start", 1, $FH_HR);
  }
  if((! verify_integer($in_stop)) || ($in_stop <= 0)) { 
    DNAORG_FAIL("ERROR in $sub_name, input stop coordinate is expected to be a positive integer, got $in_stop", 1, $FH_HR);
  }

  my ($out_start, $out_stop);
  if($seq_len == $accn_len) { # easy case
    $out_start = $in_start;
    $out_stop  = $in_stop;
  }
  elsif($seq_len == (2 * $accn_len)) { # may need to update $in_start and/or $in_stop
    # Three possible cases:
    # case 1: both $in_start and $in_stop are <= $accn_len
    #         action: do nothing in this case
    # case 2: both $in_start and $in_stop are > $accn_len
    #         action: subtract $seq_len from both
    # case 3: one of $in_start is <= $accn_len and the other is > $accn_len
    #         action: subtract $seq_len from both, subtract 1 more from either
    #                 that are now < $seq_len, to account for fact that position '0'
    #                 does not exist.
    if(($in_start <= $accn_len) && ($in_stop <= $accn_len)) { 
      # case 1
      $out_start = $in_start;
      $out_stop  = $in_stop;
    }
    elsif(($in_start > $accn_len) && ($in_stop > $accn_len)) { 
      # case 2
      $out_start = $in_start - $accn_len;
      $out_stop  = $in_stop  - $accn_len;
    }
    elsif(($in_start > $accn_len) || ($in_stop > $accn_len)) { 
      # case 3
      $out_start = $in_start - $accn_len;
      $out_stop  = $in_stop  - $accn_len;
      if($out_start <= 0) { $out_start--; }
      if($out_stop  <= 0) { $out_stop--; }
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unforeseen case in_start: $in_start in_stop: $in_stop seq_len: $seq_len accn_len: $accn_len", 1, $FH_HR);
    }
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, input sequence length ($seq_len) is neither accession length ($accn_len), nor 2 * accession length", 1, $FH_HR);
  }
  return ($out_start, $out_stop);
}

#################################################################
# Subroutine: seq_name_from_msa_seq_name()
# Incept:     EPN, Wed Mar 16 05:15:57 2016
#
# Purpose:   Convert a sequence name from an MSA created in this
#            script to the sequence name it was derived from
#            using the rule: 
#            <msa_seq_name> = <seq_name>/(?:\d+\-\d+)(?:,\d+-\d+)*/
#
# Arguments:
#  $msa_seq_name: the sequence name from the MSA
#  $FH_HR:        REF to hash of file handles
#             
# Returns:  <seq_name> derived from <msa_seq_name>
# 
# Dies:     If <msa_seq_name> is not in the expected format.
#
#################################################################
sub seq_name_from_msa_seq_name { 
  my $sub_name = "seq_name_from_msa_seq_name";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($msa_seq_name, $FH_HR) = @_;

  my $seq_name = undef;
  if($msa_seq_name =~ /^(\S+)\/(?:\d+\-\d+)(?:,\d+-\d+)*$/) { # '?:' just tells Perl not to capture group
    $seq_name = $1;
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unexpected format of msa sequence name: $msa_seq_name", 1, $FH_HR);
  }

  return $seq_name;
}

#################################################################
# Subroutine: accn_name_from_seq_name()
# Incept:     EPN, Wed Mar 16 06:16:25 2016
#
# Purpose:   Convert a sequence name from a sequence file fetched
#            in this script to the accessoin name it was derived
#            from using the rule: 
#            <seq_name> = <accn_name>\:genome.+$/
#
# Arguments:
#  $seq_name: the sequence name
#  $FH_HR:    REF to hash of file handles
#             
# Returns:  <accn_name> derived from <seq_name>
# 
# Dies:     If <seq_name> is not in the expected format.
#
#################################################################
sub accn_name_from_seq_name { 
  my $sub_name = "accn_name_from_seq_name";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $FH_HR) = @_;

  my $accn_name = undef;
  if($seq_name =~ /^(\S+)\:dnaorg.+$/) { 
    $accn_name = $1;
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, unexpected format of sequence file sequence name: $seq_name", 1, $FH_HR);
  }

  return $accn_name;
}

#################################################################
# Subroutine: process_input_fasta_file()
# Incept:     EPN, Wed May 18 10:01:02 2016
#
# Purpose:   Given the name of the input fasta file
#            (specified with --infasta), open the 
#            file, and determine the lengths of all
#            the N sequences in it. Fill 
#            $seq_info_HAR->{"accn_name"}[$i] and
#            $seq_info_HAR->{"accn_len"}[$i]
#            for sequences $i=1..N.
#
#            If -c is enabled, then we also create a 
#            new fasta file, that includes all of the
#            sequences in $infasta_file, duplicated.
#            (That is, if a sequence is:
#             >seq1
#             AGC
#             then it will become
#             >seq1-duplicated
#             AGCAGC
#             in the new file).
#
# Arguments:
#  $infasta_file:  name of the input fasta file
#  $outfasta_file: name of the output fasta file to create with duplicates of seqs from $infasta_file,
#                  undef if we don't want to do this (if -c *not* enabled)
#  $out_root:      root for naming output files
#  $seq_info_HAR:  REF to hash of arrays of sequence information, added to here
#  $opt_HHR:       REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:         REF to hash of file handles
#  
# Returns:  Number of sequences read in $infasta_file.
# 
# Dies: If $seq_info_HAR->{"accn_len"} and $seq_info_HAR->{"accn_name"} 
#       are not both arrays of exactly length 1 (with information on 
#       *only* the reference accession.)
#       
#       If the same sequence exists more than once in the input sequence
#       file.
#
#       If $outfasta_file is defined and -c is not enabled.
#
#       If -c is enabled and $outfasta_file is not defined.
#
#       If the names of any sequences in $infasta_file include the
#       string 'dnaorg', this will cause problems later on because
#       we add this to sequence names sometimes and then parse
#       based on its location, if the original sequence name already
#       has 'dnaorg', we can't guarantee that we'll be able to 
#       parse the modified sequence name correctly.
#
#################################################################
sub process_input_fasta_file { 
  my $sub_name = "process_input_fasta_file";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($infasta_file, $outfasta_file, $seq_info_HAR, $opt_HHR, $FH_HR) = @_;

  my %accn_exists_H = ();  # keeps track of which accessions have been read from the sequence file
  my $err_flag = 0;        # changed to '1' if an accession exists more than once in the input fasta file

  my $do_out = (defined $outfasta_file) ? 1 : 0;
  if($do_out && (! opt_Get("-c", $opt_HHR))) { 
    DNAORG_FAIL("ERROR in $sub_name, outfasta_file is defined, but -c not enabled", 1, $FH_HR);
  }
  if((! $do_out) && (opt_Get("-c", $opt_HHR))) { 
    DNAORG_FAIL("ERROR in $sub_name, outfasta_file is not defined, but -c is enabled", 1, $FH_HR);
  }
    
  if($do_out) { 
    open(OUT, ">", $outfasta_file) || fileOpenFailure($outfasta_file, $sub_name, $!, "writing", $FH_HR);
  }

  my $ssi_file = $infasta_file . ".ssi";
  if(-e $ssi_file) { 
    runCommand("rm $ssi_file", opt_Get("-v", $opt_HHR), $FH_HR);
  }
  my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $infasta_file }); # the sequence file object
  my $nseq = $sqfile->nseq_ssi;

  for(my $i = 0; $i < $nseq; $i++) { 
    my $next_fasta_str = $sqfile->fetch_consecutive_seqs(1, "", -1);
    # get name and length of the sequence
    if($next_fasta_str =~ /^\>(\S+).*\n(\S+)\n+$/) { 
      my ($accn, $seq) = ($1, $2);
      my $len = length($seq);
      if($accn =~ m/dnaorg/) { 
        DNAORG_FAIL("ERROR in $sub_name, sequence $accn read from $infasta_file includes the string \"dnaorg\", this is not allowed", 1, $FH_HR);
      }
      push(@{$seq_info_HAR->{"accn_name"}}, $accn);
      push(@{$seq_info_HAR->{"accn_len"}}, $len);

      if($do_out) { 
        # output the sequence if necessary
        my $outaccn = $accn . ":dnaorg-duplicated:$accn:1:$len:+:$accn:1:$len:+"; 
        # this mimics how dnaorg.pm:fetchSequencesUsingEslFetchCds names sequences, if that 
        # changes, this should be changed to match
        print OUT ">" . $outaccn . "\n";
        # print out sequence in 100-nt-per-line chunks
        print OUT $seq . "\n";
        print OUT $seq . "\n";
      }

      if(exists $accn_exists_H{$accn}) {
        $accn_exists_H{$accn}++;
        $err_flag = 1;
      }
      else { 
        $accn_exists_H{$accn} = 1;
      }
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse fasta sequence string $next_fasta_str", 1, $FH_HR);
    }
  }

  if($do_out) { 
    close OUT;
  }

  # fail if we found an accession more than once in $infasta_file
  my $errmsg = "";
  if($err_flag) { 
    $errmsg = "ERROR in $sub_name, the following accessions exist more than once in $infasta_file\nEach accession should exist only once.\n";
    for(my $i = 0; $i < scalar(@{$seq_info_HAR->{"accn_name"}}); $i++) { 
      my $accn = $seq_info_HAR->{"accn_name"}[$i];
      if((exists $accn_exists_H{$accn}) && ($accn_exists_H{$accn} > 1)) { 
        $errmsg .= "\t$accn ($accn_exists_H{$accn} occurrences)\n";
        delete $accn_exists_H{$accn};
      }
    }
    DNAORG_FAIL($errmsg, 1, $FH_HR);
  }

  return $nseq;
}

#################################################################
# Subroutine: validate_options_are_consistent_with_dnaorg_build()
# Incept:     EPN, Fri May 27 12:59:19 2016
#
# Purpose:   Given the name of the log file and checksum
#            file created by dnaorg_build.pl, read it and check 
#            that all options that need to be consistent 
#            between dnaorg_build.pl and dnaorg_annotate.pl
#            are consistent, and any files that need to have
#            the same checksums actually do.
#
# Arguments:
#  $consopts_file:     name of the dnaorg_build consopts file
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:             REF to hash of file handles
# 
# Returns:  void
# 
# Dies: If $consopts_file doesn't exist, or we can't parse it.
#
#       If an option enabled in dnaorg_build.pl that needs to 
#       be consistently used in dnaorg_annotate.pl is not, or
#       vice versa.
#
#       If an option enabled in dnaorg_build.pl that takes a file
#       as input has a different checksum for that file than 
#       
#################################################################
sub validate_options_are_consistent_with_dnaorg_build { 
  my $sub_name = "validate_options_are_consistent_with_dnaorg_build";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($consopts_file, $opt_HHR, $FH_HR) = @_;

  # read the consopts file
  if(! -e $consopts_file) { 
    DNAORG_FAIL("ERROR in $sub_name, consopts file $consopts_file does not exist.\nThis file should have been created by dnaorg_build.pl.\nYou probably need to rerun dnaorg_build.pl if it was run before May 31, 2016.", 1, $FH_HR);
  }
  open(IN, $consopts_file) || fileOpenFailure($consopts_file, $sub_name, $!, "reading", $FH_HR);
  my $line_ct = 0;
  my $no_build_c_opt         = 1; # changed to 0 below if -c was used by dnaorg_build.pl
  my $no_build_nomatpept_opt = 1; # changed to 0 below if --nomatpept was used by dnaorg_build.pl
  my $no_build_matpept_opt   = 1; # changed to 0 below if --matpept was used by dnaorg_build.pl

  while(my $line = <IN>) { 
    chomp $line;
    $line_ct++;
    if(($line eq "none") && ($line_ct == 1)) { 
      ; # this is fine, none of the options that need to be consistent were set by dnaorg_build.pl
    }
    elsif($line =~ /^\-c$/) { 
      $no_build_c_opt = 0;
      if((! defined (opt_Get("-c", $opt_HHR))) || (opt_Get("-c", $opt_HHR) != 1)) { 
        DNAORG_FAIL("ERROR, the -c option was used when dnaorg_build.pl was run (according to file $consopts_file).\nYou must also use it with dnaorg_annotate.pl.", 1, $FH_HR);
      }
    }
    elsif($line =~ /^\-\-nomatpept$/) { # first string is file name, second is md5 checksum (obtained with 'md5sum' executable)
      $no_build_nomatpept_opt = 0;
      if((! defined (opt_Get("--nomatpept", $opt_HHR))) || (opt_Get("--nomatpept", $opt_HHR) != 1)) { 
        DNAORG_FAIL("ERROR, the --nomatpept option was used when dnaorg_build.pl was run (according to file $consopts_file).\nYou must also use it with dnaorg_annotate.pl.", 1, $FH_HR);
      }
    }
    elsif($line =~ /^\-\-matpept\s+\S+\s+(\S+)$/) { # first string is file name, second is md5 checksum (obtained with 'md5sum' executable)
      my $build_matpept_cksum = $1;
      $no_build_matpept_opt = 0;
      if(! opt_IsUsed("--matpept", $opt_HHR)) { 
        DNAORG_FAIL("ERROR, the --matpept option was used when dnaorg_build.pl was run (according to file $consopts_file).\nYou must also use it with dnaorg_annotate.pl.", 1, $FH_HR);
      }
      else { # make sure checksum matches
        my $annotate_matpept_file = opt_Get("--matpept", $opt_HHR);
        if(! -s $annotate_matpept_file) { 
          DNAORG_FAIL("ERROR, the file $annotate_matpept_file specified with the --matpept option does not exist.", 1, $FH_HR);
        }          
        my $annotate_matpept_cksum = md5ChecksumOfFile($annotate_matpept_file, $sub_name, $opt_HHR, $FH_HR);
        if($build_matpept_cksum ne $annotate_matpept_cksum) { 
          DNAORG_FAIL("ERROR, the file $annotate_matpept_file specified with the --matpept file does not appear to be identical to the file used\nwith dnaorg_build.pl. The md5 checksums of the two files differ: dnaorg_build.pl: $build_matpept_cksum, dnaorg_annotate.pl: $annotate_matpept_cksum", 1, $FH_HR);
        }
      }
    }        
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse line from consopts file $consopts_file:\n$line\n", 1, $FH_HR);
    }
  }
  close(IN);

  # now for any options that were not read from $consopts_file, make sure they are also
  # not enabled here for dnaorg_annotate.pl
  if($no_build_c_opt && (opt_Get("-c", $opt_HHR))) { 
    DNAORG_FAIL("ERROR, the -c option was not used when dnaorg_build.pl was run (according to file $consopts_file).\nYou must also not use it with dnaorg_annotate.pl, or you need to rerun dnaorg_build.pl with -c.", 1, $FH_HR);
  }    
  if($no_build_nomatpept_opt && (opt_Get("--nomatpept", $opt_HHR))) { 
    DNAORG_FAIL("ERROR, the --nomatpept option was not used when dnaorg_build.pl was run (according to file $consopts_file).\nYou must also not use it with dnaorg_annotate.pl, or you need to rerun dnaorg_build.pl with --nomatpept.", 1, $FH_HR);
  }    
  if($no_build_matpept_opt && (opt_IsUsed("--matpept", $opt_HHR))) {  
    DNAORG_FAIL("ERROR, the --matpept option was not used when dnaorg_build.pl was run (according to file $consopts_file).\nYou must also not use it with dnaorg_annotate.pl, or you need to rerun dnaorg_build.pl with --matpept.", 1, $FH_HR);
  }    

  # if we get here, all options are consistent
  return;
}

################################################
## TEMPORARY SUBROUTINES for the --aorg* options
################################################

#################################################################
# Subroutine: aorg_find_origin_sequences
# Incept:     EPN, Thu Jul 28 14:46:31 2016
# 
# Purpose:    TEMPORARY function for identifying origin sequences 
#             using 'alternative' method -- a profile HMM.
#
#             Checks for and adds the following error codes: "ori".
#
# Args:   
#  $fasta_file:             the fasta file to look for origins in
#                           should contain duplicated genome sequences
#                           (so if --infasta <f> used, probably not <f>)
#  $sqfile:                 Bio:Easel object for $fasta_file
#  $execs_HR:               ref to hash with paths to executables
#  $out_root:               root for naming output files
#  $seq_info_HAR:           REF to hash of arrays with information 
#                           on the sequences, PRE-FILLED
#  $err_seq_instances_HHR:  REF to the 2D hash of per-sequence errors, initialized here
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:                REF to 2D hash of option values, 
#                           see top of epn-options.pm for description
#  $ofile_info_HHR:         REF to the 2D hash of output file information
#
# Returns:    void
#
#################################################################
sub aorg_find_origin_sequences { 
  my $sub_name = "aorg_find_origin_sequences";
  my $nargs_exp = 9;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($fasta_file, $sqfile, $execs_HR, $out_root, $seq_info_HAR, $err_seq_instances_HHR, $err_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $nseq = scalar(@{$seq_info_HAR->{"accn_name"}});

  my $aorg_model  = opt_Get("--aorgmodel",  $opt_HHR);
  my $aorg_start  = opt_Get("--aorgstart",  $opt_HHR);
  my $aorg_len    = opt_Get("--aorglen",    $opt_HHR);
  my $aorg_offset = opt_Get("--aorgoffset", $opt_HHR);
  
  # run cmscan locally
  my $tblout_file = $out_root . ".aorg.tblout";
  my $cmscan_file = $out_root . ".aorg.cmscan";
  my $opts = " --cpu 0 --tblout $tblout_file --verbose ";
  $opts .= " --nohmmonly --F1 0.02 --F2 0.001 --F2b 0.001 --F3 0.00001 --F3b 0.00001 --F4 0.0002 --F4b 0.0002 --F5 0.0002 --F6 0.0001 "; 

  my $cmd = $execs_HR->{"cmscan"} . " $opts $aorg_model $fasta_file > $cmscan_file";

  runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);

  addClosedFileToOutputInfo(\%ofile_info_HH, "aorgtblout",  "$tblout_file",     1, "tblout file from cmscan for origin identification");
  addClosedFileToOutputInfo(\%ofile_info_HH, "aorgcmscan",  "$cmscan_file",     1, "standard output file from cmscan for origin identification");

  my %hit_HH  = (); # 2D hash of top hits, 1st dim key is sequence name, 2nd is attribute, e.g. "start"    
  aorg_parse_cmscan_tblout_s2($tblout_file, $seq_info_HAR, \%hit_HH, $opt_HHR, $FH_HR);

  # Fetch all hits into a fasta file that we can align with cmalign
  my @fetch_AA = ();
  for(my $i = 0; $i < $nseq; $i++) { 
    my $seqname = $seq_info_HAR->{"seq_name"}[$i];
    if(exists $hit_HH{$seqname}) { 
      my $newname = $seqname . "/" . $hit_HH{$seqname}{"start"} . "-" . $hit_HH{$seqname}{"stop"};
      push(@fetch_AA, [ $newname, $hit_HH{$seqname}{"start"}, $hit_HH{$seqname}{"stop"}, $seqname ]);
      #printf("pushed $newname\n");
    }
  }
  my $out_fasta_file = $out_root . ".cmscan.fa";
  $sqfile->fetch_subseqs(\@fetch_AA, undef, $out_fasta_file);
  addClosedFileToOutputInfo(\%ofile_info_HH, "aorgoutfasta", "$out_fasta_file", 1, "cmscan hits in fasta format");

  # Align all hits with cmalign
  my $out_stk_file     = $out_root . ".cmalign.stk";
  my $out_cmalign_file = $out_root . ".cmalign";
  $cmd = $execs_HR->{"cmalign"} . " -o $out_stk_file $aorg_model $out_fasta_file > $out_cmalign_file";

  runCommand($cmd, 0, $ofile_info_HH{"FH"});
  addClosedFileToOutputInfo(\%ofile_info_HH, "outstk",     "$out_stk_file",     1, "alignment of cmscan hits");
  addClosedFileToOutputInfo(\%ofile_info_HH, "outcmalign", "$out_cmalign_file", 1, "cmalign output");

  # determine which positions the origin sequence is in each sequence
  my $ori_msa = Bio::Easel::MSA->new   ({ 
    fileLocation => $out_stk_file,
    isDna => 1
                                        });
  my $ori_start_rfpos    = $aorg_start;
  my $ori_stop_rfpos     = $aorg_start + $aorg_len - 1;
  my $ori_offset_rfpos   = $aorg_start + $aorg_offset - 1;
  my $ori_start_apos     = $ori_msa->rfpos_to_aligned_pos($ori_start_rfpos);
  my $ori_stop_apos      = $ori_msa->rfpos_to_aligned_pos($ori_stop_rfpos);
  my $ori_offset_apos    = $ori_msa->rfpos_to_aligned_pos($ori_offset_rfpos);

  my $ori_ppthresh = opt_Get("--aorgppthresh", $opt_HHR);
  #printf("ori_start_apos:  $ori_start_apos\n");
  #printf("ori_stop_apos:   $ori_stop_apos\n");
  #printf("ori_offset_apos: $ori_offset_apos\n");
  
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name = $seq_info_HAR->{"seq_name"}[$seq_idx];
    my $has_hit = (defined $hit_HH{$seq_name}) ? 1 : 0;
    $seq_info_HAR->{"origin_coords_str"}[$seq_idx] = ""; # initialize
    $seq_info_HAR->{"origin_offset"}[$seq_idx]     = ""; # initialize
    my $found_origin = 0;
    
    if($has_hit) { 
      my $start = $hit_HH{$seq_name}{"start"};
      my $stop  = $hit_HH{$seq_name}{"stop"};
      my $msa_seqname = $seq_name . "/" . $start . "-" . $stop;
      my $msa_sqidx = $ori_msa->get_sqidx($msa_seqname);
      if($msa_sqidx == -1) { 
        DNAORG_FAIL("ERROR, unable to find $msa_seqname in $out_stk_file", 1, $ofile_info_HH{"FH"});
      }
      # determine average pp of aligned residues between $ori_start_apos and $ori_stop_apos (inclusive)
      # and enforce our posterior probability threshold
      my ($pp_avg, $pp_cnt) = $ori_msa->get_pp_avg($msa_sqidx, $ori_start_apos, $ori_stop_apos);
      # if($pp_cnt == 0) { printf("no aligned residues for origin prediction for $seq_name\n"); }
      # if($pp_avg < $ori_ppthresh) { printf("avg PP of aligned origin below threshold for $seq_name ($pp_avg < $ori_ppthresh)\n"); }
      if(($pp_cnt > 0) && ($pp_avg >= $ori_ppthresh)) { 

        # determine the unaligned positions the origin spans in the alignment file $out_stk_file
        my ($cur_ori_start_uapos, $cur_ori_start_apos) = $ori_msa->aligned_to_unaligned_pos($msa_sqidx, $ori_start_apos, 1); # '1': if ori_start_apos is a gap, return position of 1st non-gap nucleotide after it 
        my ($cur_ori_stop_uapos,  $cur_ori_stop_apos)  = $ori_msa->aligned_to_unaligned_pos($msa_sqidx, $ori_stop_apos,  0); # '0': if ori_start_apos is a gap, return position of 1st non-gap nucleotide before it 
        
        # determine the unaligned position that maps to the special 'sequence start' position of the origin
        # (this is --aorgoffset position of the origin), which will become position number 1 when we renumber
        # sequences
        my ($cur_ori_offset_uapos, $cur_ori_offset_apos) = $ori_msa->aligned_to_unaligned_pos($msa_sqidx, $ori_offset_apos, 1); # '1': if ori_offset_apos is a gap, return position of 1st non-gap nucleotide after it 
        
        # do we have at least 1 nucleotide predicted in the origin positions?
        if(($cur_ori_start_apos < $ori_stop_apos) && ($cur_ori_stop_apos > $ori_start_apos))  { 
          # yes, we do:
          # determine strand 
          if($start <= $stop) { 
            # positive strand
            $cur_ori_start_uapos  += $start - 1;
            $cur_ori_stop_uapos   += $start - 1;
            $cur_ori_offset_uapos += $start - 1;
          }
          else { 
            # negative strand
            $cur_ori_start_uapos   = $start - $cur_ori_start_uapos  + 1;
            $cur_ori_stop_uapos    = $start - $cur_ori_stop_uapos   + 1;
            $cur_ori_offset_uapos  = $start - $cur_ori_offset_uapos + 1;
          }
          # adjust coordinates so they're within 1..L
          my ($out_ori_start_uapos, $out_ori_stop_uapos) = 
              create_output_start_and_stop($cur_ori_start_uapos, $cur_ori_stop_uapos, $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR);
          
          my $out_ori_offset_uapos;
          ($out_ori_offset_uapos, undef) = 
              create_output_start_and_stop($cur_ori_offset_uapos, $cur_ori_stop_uapos, $seq_info_HAR->{"accn_len"}[$seq_idx], $seq_info_HAR->{"seq_len"}[$seq_idx], $FH_HR);
          #printf("cur_ori_start_uapos:  $cur_ori_start_uapos\n");
          #printf("cur_ori_offset_uapos: $cur_ori_offset_uapos\n");
          #printf("cur_ori_stop_uapos:   $cur_ori_stop_uapos\n\n");
          
          #printf("out_ori_start_uapos:  $out_ori_start_uapos\n");
          #printf("out_ori_offset_uapos: $out_ori_offset_uapos\n");
          #printf("out_ori_stop_uapos:   $out_ori_stop_uapos\n\n");
          
          $seq_info_HAR->{"origin_coords_str"}[$seq_idx] .= $out_ori_start_uapos . ":" . $out_ori_stop_uapos;
          $seq_info_HAR->{"origin_offset"}[$seq_idx]      = $out_ori_offset_uapos;
          $found_origin = 1;
        } 
      } # end of 'if(($pp_cnt > 0) && ($pp_avg >= $ori_ppthresh))'
    } # end of 'if($has_hit)'
    if(! $found_origin) { 
      error_instances_add(undef, $err_seq_instances_HHR, $err_info_HAR, -1, "ori", $seq_name, "0 occurrences", $FH_HR);
    }
  } # end of 'for' loop over $seq_idx


  return;
}

#################################################################
# Subroutine : aorg_parse_cmscan_tblout_s2()
# Incept:      EPN, Tue Jul 12 08:54:07 2016
#
# Arguments: 
#  $tblout_file:   tblout file to parse
#  $seq_info_HAR:  REF to hash of arrays with information 
#                  on the sequences, PRE-FILLED
#  $hit_HHR:       REF to 2D hash of top hits, 1st dim key is sequence name, 2nd is attribute, e.g. "start"    
#  $opt_HHR:       ref to 2D options hash
#  $FH_HR:         REF to hash of file handles
#
# Returns:    void
#
#################################################################
sub aorg_parse_cmscan_tblout_s2 { 
  my $sub_name = "aorg_parse_cmscan_tblout_s2()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($tblout_file, $seq_info_HAR, $hit_HHR, $opt_HHR, $FH_HR) = @_;
  
  my %seqname_index_H = (); # seqname_index_H{$seq_name} = <n>, means that $seq_name is the <n>th sequence name in @{$seq_info_HAR{*}} arrays
  getIndexHashForArray($seq_info_HAR->{"seq_name"}, \%seqname_index_H, $FH_HR);

  open(IN, $tblout_file) || fileOpenFailure($tblout_file, $sub_name, $!, "reading", $FH_HR);

  my $ethresh = opt_Get("--aorgethresh", $opt_HHR);
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

      my $seqidx = $seqname_index_H{$seqname}; # sequence index for the hit in results_AAH (2nd dim of results_AAH)
      my $seqlen = $seq_info_HAR->{"accn_len"}[$seqidx];

      # only consider hits where either the start or end are less than the total length
      # of the genome. Since we sometimes duplicate all genomes, this gives a simple 
      # rule for deciding which of duplicate hits we'll store 
      if(($seqfrom <= $seqlen) || ($seqto <= $seqlen)) { 
        if($evalue <= $ethresh) { # enforce E-value threshold
          if(! exists $hit_HHR->{$seqname}) { # takes only the top hit
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
}


#################################################################
# Subroutine: aorg_get_origin_output_for_sequence
# Incept:     EPN, Fri Jul 29 14:45:43 2016
#
# Synopsis:   For a given sequence index $seq_idx, determine 
#             the output strings related to the origin sequence
#             using the new profile HMM based method.
#
# Args: 
#  $seq_info_HAR:  REF to hash of arrays with information 
#                  on the sequences (including origins), PRE-FILLED
#  $seq_idx:       index of sequence we're working on
#  $FH_HR:         REF to hash of file handles
#
# Returns: 5 values:
#          $oseq_ct:       number of occurrences of origin sequence found
#          $oseq_start:    start position of origin seq if $oseq_ct == 1, else '-'
#          $oseq_stop:     stop  position of origin seq if $oseq_ct == 1, else '-'
#          $oseq_firstpos: what should be the first position of the seq if $oseq_ct == 1, else '-'
#          $oseq_offset:   offset position of origin seq if $oseq_ct == 1, else '-'
#          $oseq_passfail: 'P' if $oseq_ct is 1, else 'F'
#
# Dies: if $seq_info_HAR->{"origin_coords_str"}[$seq_idx] does not exist.
# 
#################################################################
sub aorg_get_origin_output_for_sequence { 
  my $sub_name = "aorg_get_origin_output_for_sequence";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_info_HAR, $seq_idx, $FH_HR) = @_;

  if(! exists $seq_info_HAR->{"origin_coords_str"}[$seq_idx]) { 
    DNAORG_FAIL("ERROR in $sub_name, no origin_coords_str (not even empty str) exists in seq_info_HAR for $seq_idx", 1, $FH_HR);
  }
  my $accn_len = $seq_info_HAR->{"accn_len"}[$seq_idx];

  # initializereturn values
  my $oseq_start    = "-"; # changed below if $oseq_ct == 1
  my $oseq_stop     = "-"; # changed below if $oseq_ct == 1
  my $oseq_firstpos = "-"; # changed below if $oseq_ct == 1
  my $oseq_offset   = "-"; # changed below if $oseq_ct == 1
  my $oseq_passfail = "F"; # changed below if $oseq_ct == 1

  my @coords_A = split(",", $seq_info_HAR->{"origin_coords_str"}[$seq_idx]);
  my $oseq_ct = scalar(@coords_A);

  if($oseq_ct == 1) { 
    ($oseq_start, $oseq_stop) = split(":", $coords_A[0]);
    $oseq_firstpos = $seq_info_HAR->{"origin_offset"}[$seq_idx];
    if($oseq_firstpos < 0) { 
      $oseq_firstpos = $accn_len + $oseq_firstpos + 1;
    }
    $oseq_offset = $oseq_firstpos; 

    # $oseq_offset is now what should be the first position of the sequence, 
    # if this is > 0, we subtract 1 because that's how many we need to shift counterclockwise 
    # direction to put $oseq_offset as position 1 (imagine if $oseq_offset is 1, then we need to shift '0', not '1')
    # if it's negative then it's already the correct shift (due to the off-by-one
    # introduced because there is no position 0)
    if($oseq_offset > 0) { $oseq_offset -= 1; }
    # printf("HEYA in $sub_name, seq_idx: $seq_idx oseq_offset: $oseq_offset\n");
    # $oseq_offset is now number of nts to shift origin in counterclockwise direction
    if($oseq_offset > ($accn_len / 2)) { # simpler (shorter distance) to move origin clockwise
      $oseq_offset = $accn_len - $oseq_offset; # note, we don't add 1 here
    }
    else { # simpler to shift origin in counterclockwise direction, we denote this as a negative offset
      $oseq_offset *= -1;
    }

    $oseq_passfail = "P";
  }

  return ($oseq_ct, $oseq_start, $oseq_stop, $oseq_firstpos, $oseq_offset, $oseq_passfail);
}
