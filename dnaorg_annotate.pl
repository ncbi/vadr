#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Thu Feb 18 12:48:16 2016 [dnaorg_annotate.pl split off from dnaorg_annotate_genomes.pl]
#
use strict;
use warnings;
#use diagnostics;
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
# Step 19. Run BLASTX: all full length sequences and all corrected nucleotide 
#          features versus all proteins
# Step 20. Add zft errors for sequences with zero annotated features
# Step 21. Translate corrected feature sequences into proteins
# Step 22. (OPTIONAL) Create multiple alignments of DNA sequences, if --doalign.
# Step 23. (OPTIONAL) Create multiple alignments of protein sequences, if --doalign.
# Step 24. Output annotations and errors.
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
# 4. mdl_results_add_str_nop_ost_lsc_dup_b3e_b3u_errors()
# 5. ftr_results_calculate()
# 6. find_origin_sequences()
# 7. MAIN (not a function but rather the main body of the script):
# 8. mdl_results_add_b5e_b5u_errors()
# 9. ftr_results_add_b5e_errors()
#
#              annot_type
#          -------------------
# errcode  model  multifeature sequence
# -------  -----  ------------ --------
# nop      4      N/A          N/A
# nm3      5      5            N/A
# b5e      8,9    N/A          N/A
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
# mtr      5      N/A          N/A
# nst      1,7    1,7          N/A
# ost      4      N/A          N/A
# aji      N/A    5            N/A
# int      N/A    5            N/A
# inp      N/A    5            N/A
# ori      N/A    N/A          6
# 
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
my $esl_epn_translate = $dnaorgdir . "/esl-epn-translate/esl-epn-translate.pl";
my $esl_ssplit        = $dnaorgdir . "/Bio-Easel/scripts/esl-ssplit.pl";
my $blast_exec_dir    = "/usr/bin/"; # HARD-CODED FOR NOW

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
my $g = 0; # option group

# Add all options to %opt_HH and @opt_order_A.
# This section needs to be kept in sync (manually) with the &GetOptions call below
#     option            type       default               group   requires incompat    preamble-output                                   help-output    
opt_Add("-h",           "boolean", 0,                        0,    undef, undef,      undef,                                            "display this help",                                  \%opt_HH, \@opt_order_A);
$opt_group_desc_H{++$g} = "REQUIRED options";
opt_Add("--infasta",     "string",  undef,                  $g,    "*",   "*",       "fasta file with sequences to anntoate is <s>",    "fasta file with sequences to annotate is <s>",       \%opt_HH, \@opt_order_A);
opt_Add("--refaccn",     "string",  undef,                  $g,    "*",   "*",       "reference accession to annotate based on is <s>", "reference accession to annotate based on is <s>",    \%opt_HH, \@opt_order_A);
opt_Add("--dirout",      "string",  undef,                  $g,    "*",   "*",       "output directory to create is <s>",               "output directory to create is <s>",                  \%opt_HH, \@opt_order_A);
opt_Add("--dirbuild",    "string",  undef,                  $g,    "*",   "*",       "output directory created by dnaorg_build.pl",     "output directory created by dnaorg_build.pl is <s>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "basic options";
opt_Add("-c",           "boolean", 0,                       $g,    undef, undef,      "genome is closed (circular)",                                 "genome is closed (circular)",                       \%opt_HH, \@opt_order_A);
opt_Add("-f",           "boolean", 0,                       $g,    undef,undef,       "forcing directory overwrite",                                 "force; if dir from --dirout exists, overwrite it",   \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                       $g,    undef, undef,      "be verbose",                                                  "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--origin",     "string",  undef,                   $g,     "-c", undef,      "identify origin seq <s> in genomes",                          "identify origin seq <s> in genomes, put \"|\" at site of origin (\"|\" must be escaped, i.e. \"\\|\"", \%opt_HH, \@opt_order_A);
opt_Add("--matpept",    "string",  undef,                   $g,    undef, undef,      "using pre-specified mat_peptide info",                        "read mat_peptide info in addition to CDS info, file <s> explains CDS:mat_peptide relationships", \%opt_HH, \@opt_order_A);
opt_Add("--nomatpept",  "boolean", 0,                       $g,    undef,"--matpept", "ignore mat_peptide annotation",                               "ignore mat_peptide information in reference annotation", \%opt_HH, \@opt_order_A);
opt_Add("--xfeat",      "string",  undef,                   $g,    undef, undef,      "use models of additional qualifiers",                         "use models of additional qualifiers in string <s>", \%opt_HH, \@opt_order_A);  
opt_Add("--dfeat",      "string",  undef,                   $g,    undef, undef,      "annotate additional qualifiers as duplicates", "annotate qualifiers in <s> from duplicates (e.g. gene from CDS)",  \%opt_HH, \@opt_order_A);  
opt_Add("--specstart",  "string",  undef,                   $g,    undef, undef,      "using pre-specified alternate start codons",                  "read specified alternate start codons per CDS from file <s>", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                       $g,    undef, undef,      "leaving intermediate files on disk",                          "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for alternative modes";
#        option               type   default                group  requires incompat                        preamble-output                                                      help-output    

$opt_group_desc_H{++$g} = "options for tuning protein validation with blastx";
#        option               type   default                group  requires incompat   preamble-output                                                                                 help-output    
opt_Add("--xalntol",     "integer",  5,                      $g,     undef, undef,     "max allowed difference in nucleotides b/t nucleotide and blastx start/end predictions is <n>", "max allowed difference in nucleotides b/t nucleotide and blastx start/end postions is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--xindeltol",   "integer",  27,                     $g,     undef, undef,     "max allowed nucleotide insertion and deletion length in blastx validation is <n>",             "max allowed nucleotide insertion and deletion length in blastx validation is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--xlonescore",  "integer",  80,                     $g,     undef, undef,     "minimum score for a lone blastx hit (not supported by a CM hit) to cause an error ",           "minimum score for a lone blastx (not supported by a CM hit) to cause an error is <n>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for modifying which errors are reported";
#       option               type   default                group  requires incompat preamble-output                                     help-output    
opt_Add("--classerrors","string",  0,                      $g,    undef,   undef,   "read per-sequence classification errors from <s>", "read per-sequence classification errors from <s>", \%opt_HH, \@opt_order_A);
opt_Add("--allolp",     "boolean", 0,                      $g,    undef,   undef,   "report all olp errors, do not skip due to nop",    "report all olp errors, even when other feature is not predicted (nop error)", \%opt_HH, \@opt_order_A);
opt_Add("--alladj",     "boolean", 0,                      $g,    undef,   undef,   "report all adj errors, do not skip due to nop",    "report all aja/ajb errors, even when other feature is not predicted (nop error)", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for modifying cmalign runs";
#        option               type   default                group  requires incompat   preamble-output                                                                help-output    
opt_Add("--mxsize",     "integer", 8000,                    $g,    undef, undef,      "set max allowed dp matrix size --mxsize value for cmalign calls to <n> Mb",    "set max allowed dp matrix size --mxsize value for cmalign calls to <n> Mb", \%opt_HH, \@opt_order_A);
opt_Add("--tau",        "real",    1E-7,                    $g,    undef, undef,      "set the initial tau value for cmalign to <x>",                                 "set the initial tau value for cmalign to <x>", \%opt_HH, \@opt_order_A);
opt_Add("--noglocal",   "boolean", 0,                       $g,    undef, "-c",       "do not run cmalign in glocal mode (run in local mode)",                        "do not run cmalign in glocal mode (run in local mode)", \%opt_HH, \@opt_order_A);
opt_Add("--nofixedtau", "boolean", 0,                       $g,    undef, "-c",       "do not fix the tau value when running cmalign, allow it to decrease if nec",   "do not fix the tau value when running cmalign, allow it to decrease if nec", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options related to parallelization on compute farm";
#        option               type   default                group  requires incompat   preamble-output                                                                                 help-output    
opt_Add("--local",      "boolean", 0,                       $g,    undef, undef,      "run cmscan locally instead of on farm",                       "run cmscan locally instead of on farm", \%opt_HH, \@opt_order_A);
opt_Add("--errcheck",   "boolean", 0,                       $g,    undef,"--local",   "consider any farm stderr output as indicating a job failure", "consider any farm stderr output as indicating a job failure", \%opt_HH, \@opt_order_A);
opt_Add("--nkb",        "integer", 50,                      $g,    undef,"--local",   "number of KB of sequence for each cmscan farm job",           "set target number of KB of sequences for each cmscan farm job to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--maxnjobs",   "integer", 2500,                    $g,    undef,"--local",   "maximum allowed number of jobs for compute farm",             "set max number of jobs to submit to compute farm to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--wait",       "integer", 500,                     $g,    undef,"--local",   "allow <n> minutes for cmscan jobs on farm",                   "allow <n> wall-clock minutes for cmscan jobs on farm to finish, including queueing time", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for skipping/adding optional stages";
#       option               type   default                group  requires incompat preamble-output                             help-output    
opt_Add("--doalign",    "boolean", 0,                      $g,    undef,   undef,   "create nucleotide and protein alignments", "create nucleotide and protein alignments", \%opt_HH, \@opt_order_A);
opt_Add("--amxsize",    "integer", 2048,                   $g, "--doalign",undef,   "with --doalign, set --mxsize <n> to <n>",  "with --doalign, set --mxsize <n> for cmalign to <n>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options that modify the tabular output file";
#       option               type   default                group  requires incompat preamble-output                                               help-output    
opt_Add("--tblfirst",    "boolean", 0,                     $g,    undef,   undef,   "put first accession first on each .tbl page",               "include annotation for first accession on each page of .tbl output file", \%opt_HH, \@opt_order_A);
opt_Add("--tblnocomp",   "boolean", 0,                     $g,    undef,   undef,   "do not compare annotations to existing GenBank annotation", "do not include information comparing predicted annotations to existing GenBank annotations", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "optional output files";
#       option       type       default                  group  requires incompat  preamble-output                          help-output    
opt_Add("--mdlinfo",    "boolean", 0,                       $g,    undef, undef, "output internal model information",     "create file with internal model information",   \%opt_HH, \@opt_order_A);
opt_Add("--ftrinfo",    "boolean", 0,                       $g,    undef, undef, "output internal feature information",   "create file with internal feature information", \%opt_HH, \@opt_order_A);
opt_Add("--seqinfo",    "boolean", 0,                       $g,    undef, undef, "output internal sequence information",  "create file with internal sequence information", \%opt_HH, \@opt_order_A);
opt_Add("--errinfo",    "boolean", 0,                       $g,    undef, undef, "output internal error information",     "create file with internal error information", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for skipping stages and using files from earlier, identical runs, primarily useful for debugging";
#     option               type       default            group   requires    incompat                    preamble-output                                            help-output    
opt_Add("--skipedirect",   "boolean", 0,                    $g,   undef,      "-f,--nkb,--maxnjobs,--local,--wait", "skip the edirect steps, use existing results",           "skip the edirect steps, use data from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipfetch",     "boolean", 0,                    $g,   undef,      "-f,--nkb,--maxnjobs,--local,--wait", "skip the sequence fetching steps, use existing results", "skip the sequence fetching steps, use files from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skipalign",     "boolean", 0,                    $g,   undef,      "-f,--nkb,--maxnjobs,--local,--wait", "skip the cmalign step, use existing results",             "skip the cmscan step, use results from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--skiptranslate", "boolean", 0,                    $g,"--skipalign",  undef,                      "skip the translation steps, use existing results",       "skip the translation steps, use results from an earlier run of the script", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "TEMPORARY options for the alternative method of identifying origin sequences";
#     option               type       default            group   requires                                   incompat     preamble-output                                                         help-output    
opt_Add("--aorgmodel",     "string",  undef,                $g,   "-c,--aorgstart,--aorgoffset,--aorglen",   "--origin",  "use alternative origin method with model <s>",                         "use alternative origin method with origin model in <s>", \%opt_HH, \@opt_order_A);
opt_Add("--aorgstart",     "integer", 0,                    $g,   "-c,--aorgmodel,--aorgoffset,--aorglen",   "--origin",  "origin begins at position <n> in --aorgmodel model",                   "origin begins at position <n> in --aorgmodel model",     \%opt_HH, \@opt_order_A);
opt_Add("--aorgoffset",    "integer", 0,                    $g,   "-c,--aorgmodel,--aorgstart,--aorglen",    "--origin",  "first position of genome sequence is position <n> in origin sequence", "first position of genome sequence is position <n> in origin sequence", \%opt_HH, \@opt_order_A);
opt_Add("--aorglen",       "integer", 0,                    $g,   "-c,--aorgmodel,--aorgstart,--aorgoffset", "--origin",  "length of origin sequence is <n>",                                     "length of origin sequence is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--aorgethresh",   "real",    1.0,                  $g,   "-c,--aorgmodel,--aorgstart,--aorgoffset", "--origin",  "E-value threshold for origin detection is <x>",                        "E-value threshold for origin detection is <x>", \%opt_HH, \@opt_order_A);
opt_Add("--aorgppthresh",  "real",    0.6,                  $g,   "-c,--aorgmodel,--aorgstart,--aorgoffset", "--origin",  "average PP threshold for origin detection is <x>",                     "average PP threshold for origin detection is <x>", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_annotate.pl\n";
$usage      .= "\t--infasta  <sequence fasta file to annotate>              (REQUIRED option)\n";
$usage      .= "\t--refaccn  <reference accession to annotate with>         (REQUIRED option)\n";
$usage      .= "\t--dirbuild <path to directory created by dnaorg_build.pl> (REQUIRED option)\n";
$usage      .= "\t--dirout   <path to output directory to create>           (REQUIRED option)\n";
$usage      .= "\t[additional options]\n";
$usage      .= "\n";
my $synopsis = "dnaorg_annotate.pl :: annotate sequences based on a reference annotation";
my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# REQUIRED options
                'infasta=s'    => \$GetOptions_H{"--infasta"},
                'refaccn=s'    => \$GetOptions_H{"--refaccn"},
                'dirbuild=s'   => \$GetOptions_H{"--dirbuild"},
                'dirout=s'     => \$GetOptions_H{"--dirout"},
# basic options
                'c'            => \$GetOptions_H{"-c"},
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'origin=s'     => \$GetOptions_H{"--origin"},
                'matpept=s'    => \$GetOptions_H{"--matpept"},
                'nomatpept'    => \$GetOptions_H{"--nomatpept"},
                'xfeat=s'      => \$GetOptions_H{"--xfeat"},
                'dfeat=s'      => \$GetOptions_H{"--dfeat"},
                'specstart=s'  => \$GetOptions_H{"--specstart"},
                'keep'         => \$GetOptions_H{"--keep"},
# options for tuning protein validation with blastx
                'xalntol=s'    => \$GetOptions_H{"--xalntol"},
                'xindeltol=s'  => \$GetOptions_H{"--xindeltol"},
                'xlonescore=s' => \$GetOptions_H{"--xlonescore"},
# options for modifying which errors are reported
                'classerrors=s' => \$GetOptions_H{"--classerrors"},
                'allolp'        => \$GetOptions_H{"--allolp"},
                'alladj'        => \$GetOptions_H{"--alladj"},
# options for changing search sensitivity modes
                'mxsize=s'     => \$GetOptions_H{"--mxsize"},
                'tau=s'        => \$GetOptions_H{"--tau"},
                'noglocal'     => \$GetOptions_H{"--noglocal"},
                'nofixedtau'   => \$GetOptions_H{"--nofixedtau"},
# options related to parallelization
                'local'        => \$GetOptions_H{"--local"}, 
                'errcheck'     => \$GetOptions_H{"--errcheck"},  
                'nkb=s'        => \$GetOptions_H{"--nkb"}, 
                'maxnjobs=s'   => \$GetOptions_H{"--maxnjobs"}, 
                'wait=s'       => \$GetOptions_H{"--wait"},
# options for skipping/adding optional stages
                'doalign'      => \$GetOptions_H{"--doalign"},
                'amxsize=s'    => \$GetOptions_H{"--amxsize"},
# options that affect tabular output file
                'tblfirst'     => \$GetOptions_H{"--tblfirst"},
                'tblnocomp'    => \$GetOptions_H{"--tblnocomp"},
# optional output files
                'mdlinfo'      => \$GetOptions_H{"--mdlinfo"},
                'ftrinfo'      => \$GetOptions_H{"--ftrinfo"}, 
                'seqinfo'      => \$GetOptions_H{"--seqinfo"}, 
                'errinfo'      => \$GetOptions_H{"--errinfo"},
# options for skipping stages, using earlier results
                'skipedirect'   => \$GetOptions_H{"--skipedirect"},
                'skipfetch'     => \$GetOptions_H{"--skipfetch"},
                'skipalign'     => \$GetOptions_H{"--skipalign"},
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
my $version       = "0.45";
my $releasedate   = "Feb 2019";

# make *STDOUT file handle 'hot' so it automatically flushes whenever we print to it
# it is printed to
select *STDOUT;
$| = 1;

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date, $dnaorgdir);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 0) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do dnaorg_annotate.pl -h\n\n";
  exit(1);
}

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

my $orig_infasta_file = opt_Get("--infasta", \%opt_HH);
my $dir_build  = opt_Get("--dirbuild", \%opt_HH);  # this will be undefined unless --dirbuild set on cmdline
my $dir_out    = opt_Get("--dirout",   \%opt_HH);  # this will be undefined unless --dirout set on cmdline
my $do_matpept = opt_IsOn("--matpept", \%opt_HH);  # this will be '0' unless --matpept set on cmdline 
my $do_infasta = 1;

$dir_out =~ s/\/$//; # remove final '/' if there is one
$dir_build =~ s/\/$//; # remove final '/' if there is one

if($dir_out eq $dir_build) { 
  DNAORG_FAIL("ERROR, with --dirout <s1> and --dirbuild <s2>, <s1> and <s2> must be different directories", 1, undef);
}
# make sure build directory exists
if(! -d $dir_build) {
  DNAORG_FAIL("ERROR, directory $dir_build (specified with --dirbuild) does not exist.\nDid you run \"dnaorg_build.pl --dirout $dir_build\" yet? If not, you need to do that first.", 1, undef);
}

###############
# Preliminaries
###############
my $cmd;               # a command to run with runCommand()
my @early_cmd_A = ();  # array of commands we run before our log file is opened
my %seq_info_HA = ();  # hash of arrays, values are arrays with index range [0..$nseq-1];
                       # 1st dim keys are "seq_name", "accn_name", "len"
                       # $seq_info_HA{"accn_name"}[0] is our reference accession
@{$seq_info_HA{"accn_name"}} = ();

my %infasta_ref_seq_info_HA = ();  # hash of arrays, for reference sequence information. 
                                   # only used if --infasta used. Actually only stores information
                                   # on 1 sequence, so could be just a hash, but it is a hash of 
                                   # single element arrays so that it is the same type of data
                                   # structure as %seq_info_HA so we can pass it into 
                                   # functions (namely wrapperGetInfoUsingEdirect) in place
                                   # of %seq_info_HA.
                                   # 1st dim keys are "seq_name", "accn_name", "len"
                                   # $infasta_ref_seq_info_HA{"accn_name"}[0] is our reference accession

my $nseq = 0;
my $ref_accn = opt_Get("--refaccn", \%opt_HH);
stripVersion(\$ref_accn);
if($dir_out eq $ref_accn) { 
  DNAORG_FAIL("ERROR, with --dirout <s1> and --refaccn <s2>, <s1> and <s2> must be different", 1, undef);
}

# remove dirout if -f used
# check if one of the skip options was used (begin with --skip) 
# if so, try to use it. Else tell user to either rerun with -f
# or delete it.
if(-d $dir_out) { 
  $cmd = "rm -rf $dir_out";
  if(opt_Get("-f", \%opt_HH)) { # -f used, always remove it
    runCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); 
  }
  else { # dirout exists but -f not used
    if(! ((opt_IsUsed("--skipedirect",   \%opt_HH)) || 
          (opt_IsUsed("--skipfetch",     \%opt_HH)) || 
          (opt_IsUsed("--skipalign",      \%opt_HH)) || 
          (opt_IsUsed("--skiptranslate", \%opt_HH)))) { 
      die "ERROR directory named $dir_out (specified with --dirout) already exists. Remove it, or use -f to overwrite it."; 
    }
    # if a --skip option is used, we just press on
  }
}
elsif(-e $dir_out) { 
  $cmd = "rm $dir_out";
  if(opt_Get("-f", \%opt_HH)) { runCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR a file named $dir_out (specified with --dirout) already exists. Remove it, or use -f to overwrite it."; }
}

# if $dir_out does not exist, create it
if(! -d $dir_out) {
  $cmd = "mkdir $dir_out";
  runCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd);
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
outputBanner($log_FH, $version, $releasedate, $synopsis, $date, $dnaorgdir);
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

# --classerrors <f>
my $do_class_errors = (opt_IsUsed("--classerrors", \%opt_HH)) ? 1 : 0;
my $class_errors_file      = undef;
my %class_errors_per_seq_H = ();
if($do_class_errors) { 
  $class_errors_file = opt_Get("--classerrors", \%opt_HH);
  if(! -s $class_errors_file) { 
    die "ERROR file $class_errors_file specified with --classerrors does not exist or is empty"; 
  }
  parse_class_errors_list_file($class_errors_file, \%class_errors_per_seq_H, $ofile_info_HH{"FH"});
}

###################################################
# make sure the required executables are executable
###################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmscan"}            = $inf_exec_dir   . "cmscan";
$execs_H{"cmalign"}           = $inf_exec_dir   . "cmalign";
$execs_H{"hmmbuild"}          = $hmmer_exec_dir . "hmmbuild";
$execs_H{"hmmalign"}          = $hmmer_exec_dir . "hmmalign";
$execs_H{"esl-reformat"}      = $esl_exec_dir   . "esl-reformat";
$execs_H{"esl_fetch_cds"}     = $esl_fetch_cds;
$execs_H{"esl-epn-translate"} = $esl_epn_translate;
$execs_H{"esl-ssplit"}        = $esl_ssplit;
$execs_H{"blastx"}            = $blast_exec_dir . "blastx";
$execs_H{"parse_blastx"}      = $dnaorgdir . "/dnaorg_scripts/parse_blastx.pl";
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

###########################################################################
# Step 0. Preliminaries:
#         - Read the dnaorg_build.pl consopts file and make sure that it
#           agrees with the options set here.
#         - Initialize error-related data structures.
#            
###########################################################################
my $progress_w = 85; # the width of the left hand column in our progress output, hard-coded
my $start_secs = outputProgressPrior("Verifying options are consistent with options used for dnaorg_build.pl", $progress_w, $log_FH, *STDOUT);
validate_options_are_consistent_with_dnaorg_build($build_root . ".consopts", \%opt_HH, $ofile_info_HH{"FH"});
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

# initialize error related data structures
my %err_info_HA = (); 
initializeHardCodedErrorInfoHash(\%err_info_HA, $ofile_info_HH{"FH"});

###########################################################################
# Step 1. Gather and process information on reference genome using Edirect.
###########################################################################
my $progress_str = undef;
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
my %xfeat_tbl_HHHA = (); # xfeat (eXtra feature) data from feature table file, hash of hash of hashes of arrays
                         # 1D: qualifier name, e.g. 'gene'
                         # 2D: key: accession
                         # 3D: key: column name in gene ftable file
                         # 4D: per-row values for each column
my %dfeat_tbl_HHHA = (); # dfeat (Duplicate feature) data from feature table file, hash of hash of hashes of arrays
                         # 1D: qualifier name, e.g. 'gene'
                         # 2D: key: accession
                         # 3D: key: column name in gene ftable file
                         # 4D: per-row values for each column

# parse --xfeat option if necessary and initiate hash of hash of arrays for each comma separated value
my $do_xfeat = 0;
if(opt_IsUsed("--xfeat", \%opt_HH)) { 
  $do_xfeat = 1;
  my $xfeat_str = opt_Get("--xfeat", \%opt_HH);
  foreach my $xfeat (split(",", $xfeat_str)) { 
    %{$xfeat_tbl_HHHA{$xfeat}} = ();
  }
}
# parse --dfeat option if necessary
my $do_dfeat = 0;
if(opt_IsUsed("--dfeat", \%opt_HH)) { 
  $do_dfeat = 1;
  my $dfeat_str = opt_Get("--dfeat", \%opt_HH);
  foreach my $dfeat (split(",", $dfeat_str)) { 
    %{$dfeat_tbl_HHHA{$dfeat}} = ();
    if(exists $xfeat_tbl_HHHA{$dfeat}) {
      DNAORG_FAIL("ERROR, with --xfeat <s1> and --dfeat <s2>, no qualifier names can be in common between <s1> and <s2>, found $dfeat", 1, $ofile_info_HH{"FH"});
    }
  }
}

# Call the wrapper function that does the following:
#  1) creates the edirect .mat_peptide file, if necessary
#  2) creates the edirect .ftable file
#  3) creates the length file
#  4) parses the edirect .mat_peptide file, if necessary
#  5) parses the edirect .ftable file
#  6) parses the length file
# make a copy of the fasta file in the current directory
my $infasta_file = $out_root . ".in.fa";
runCommand("cp $orig_infasta_file $infasta_file", opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});;
# note that we pass in a reference to %ref_seq_info_HA to wrapperGetInfoUsingEdirect()
# and *not* a reference to %seq_info_HA. We will use %infasta_ref_seq_info_HA to 
# store information on the reference sequence only.
wrapperGetInfoUsingEdirect(undef, $ref_accn, $build_root, \%cds_tbl_HHA, \%mp_tbl_HHA, \%xfeat_tbl_HHHA, \%dfeat_tbl_HHHA, \%infasta_ref_seq_info_HA, \%ofile_info_HH,
                           \%opt_HH, $ofile_info_HH{"FH"}); 
$nseq = process_input_fasta_file($infasta_file, \%seq_info_HA, \%opt_HH, $ofile_info_HH{"FH"}); 

if($do_matpept) {  
  # validate the CDS:mat_peptide relationships that we read from the $matpept input file
  matpeptValidateCdsRelationships(\@cds2pmatpept_AA, \%{$cds_tbl_HHA{$ref_accn}}, \%{$mp_tbl_HHA{$ref_accn}}, opt_Get("-c", \%opt_HH), $seq_info_HA{"len"}[0], $ofile_info_HH{"FH"});
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
                                                    ($do_infasta) ? $infasta_ref_seq_info_HA{"len"}[0]       : undef,
                                                    ($do_infasta) ? $infasta_file                            : undef,
                                                    \%cds_tbl_HHA,
                                                    ($do_matpept) ? \%mp_tbl_HHA      : undef, 
                                                    ($do_xfeat)   ? \%xfeat_tbl_HHHA  : undef,
                                                    ($do_dfeat)   ? \%dfeat_tbl_HHHA  : undef,
                                                    ($do_matpept) ? \@cds2pmatpept_AA : undef, 
                                                    ($do_matpept) ? \@cds2amatpept_AA : undef, 
                                                    \%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, 
                                                    \%opt_HH, \%ofile_info_HH);

# verify our model, feature, and sequence info hashes are complete, 
# if validateFeatureInfoHashIsComplete() fails then the program will exit with an error message
my $nftr = validateFeatureInfoHashIsComplete  (\%ftr_info_HA, undef, $ofile_info_HH{"FH"}); # nftr: number of features
my $nmdl = validateModelInfoHashIsComplete    (\%mdl_info_HA, undef, $ofile_info_HH{"FH"}); # nmdl: number of homology models
if($nseq != validateSequenceInfoHashIsComplete(\%seq_info_HA, undef, \%opt_HH, $ofile_info_HH{"FH"})) { 
  DNAORG_FAIL(sprintf("ERROR, number of stored sequences (%d) in seq_info_HA differs from number of accessions read from $orig_infasta_file (%d)", validateSequenceInfoHashIsComplete(\%seq_info_HA, undef, \%opt_HH, $ofile_info_HH{"FH"}), $nseq), 1, $ofile_info_HH{"FH"});
}    
# also verify that we have all the blastx db files that we need
validateBlastDbExists($build_root . ".prot.fa", undef);
for(my $tmp_f = 0; $tmp_f < $nftr; $tmp_f++) { 
  if(($ftr_info_HA{"type"}[$tmp_f] eq "cds-mp") || 
     ($ftr_info_HA{"type"}[$tmp_f] eq "cds-notmp")) { 
    validateBlastDbExists(($build_root . ".f" . $tmp_f . ".prot.fa"), $ofile_info_HH{"FH"});
  }
}

# now that we have the ftr_info_HA filled, we can initialize the error data structures
my @err_ftr_instances_AHH = ();
my %err_seq_instances_HH = ();
error_instances_initialize_AHH(\@err_ftr_instances_AHH, \%err_seq_instances_HH, \%err_info_HA, \%ftr_info_HA, $ofile_info_HH{"FH"});

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###################################################################
# Step 3. Verify we have the model file that we need to run cmalign
###################################################################
my $model_file = $build_root . ".cm";
if(! -s $model_file) { 
  DNAORG_FAIL("ERROR CM file $model_file should exist but it does not. Did you (successfully) run dnaorg_build.pl?", 1, $ofile_info_HH{"FH"});
}
for(my $i = 0; $i < $nmdl; $i++) { 
  # set mdl_info_HAR->{"cmfile"}[$i]
  $mdl_info_HA{"cmfile"}[$i] = $model_file;
}

###################################################################
# Step 4. (OPTIONAL) Search for origin sequences, if --origin used
###################################################################
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

#########################
# Step 5. Align sequences
#########################
my $seq_file = $ofile_info_HH{"fullpath"}{"fasta"};
validateFileExistsAndIsNonEmpty($seq_file, "dnaorg_annotate.pl:main", $ofile_info_HH{"FH"});
my $tot_len_nt = sumArray(\@{$seq_info_HA{"len"}});
my $mdl_file = $mdl_info_HA{"cmfile"}[0];
$mdl_file =~ s/\.\d+\.cm$/.cm/; 

my $ifile_file = $out_root . ".ifile"; # concatenated --ifile output files, created by concatenating all of the individual 
                                       # ifile files in cmalignOrNhmmscanWrapper()
my @stk_file_A        = (); # array of all stk output files created in cmalignOrNhmmscanWrapper()
my @overflow_seq_A    = (); # array of sequences that fail cmalign b/c required matrix was too big
my @overflow_mxsize_A = (); # array of required matrix sizes for each sequence in @overflow_seq_A
my $cmalign_stdout_file = $out_root . ".cmalign.stdout";
my $cmalign_ifile_file  = $out_root . ".cmalign.ifile";
cmalignOrNhmmscanWrapper(\%execs_H, 1, $out_root, $seq_file, $tot_len_nt, $progress_w, 
                         $mdl_file, \@stk_file_A, \@overflow_seq_A, \@overflow_mxsize_A, \%opt_HH, \%ofile_info_HH);

my $ndmo_errors = scalar(@overflow_seq_A);
if($ndmo_errors > 0) { 
  add_dmo_errors(\@overflow_seq_A, \@overflow_mxsize_A, \%err_seq_instances_HH, \%err_info_HA, \%opt_HH, \%ofile_info_HH);
}

####################################################################
# Step 6. Parse homology search results into usable data structures
####################################################################
$start_secs = outputProgressPrior("Parsing cmalign results", $progress_w, $log_FH, *STDOUT);

my @mdl_results_AAH = ();  # 1st dim: array, 0..$nmdl-1, one per model
                           # 2nd dim: array, 0..$nseq-1, one per sequence
                           # 3rd dim: hash, keys are "p_start", "p_stop", "p_strand", "p_5overhang", "p_3overhang", "p_5seqflush", "p_3seqflush", "p_evalue", "p_fid2ref"
initialize_mdl_results(\@mdl_results_AAH, \%mdl_info_HA, \%seq_info_HA, \%opt_HH, $ofile_info_HH{"FH"});

# make an 'order hash' for the sequence names,
my %seq_name_index_H = (); # seq_name_index_H{$seq_name} = <n>, means that $seq_name is the <n>th sequence name in the @{$seq_name_AR}} array
getIndexHashForArray($seq_info_HA{"seq_name"}, \%seq_name_index_H, $ofile_info_HH{"FH"});

if($nseq > $ndmo_errors) { # at least 1 sequence was aligned
  parse_cmalign_ifile($cmalign_ifile_file, \%seq_name_index_H, \%seq_info_HA, $ofile_info_HH{"FH"});

  # parse the cmalign alignments
  for(my $a = 0; $a < scalar(@stk_file_A); $a++) { 
    if(-s $stk_file_A[$a]) { # skip empty alignments, which will exist for any r1 run that fails
      parse_cmalign_stk($stk_file_A[$a], \%seq_name_index_H, \%seq_info_HA, \%mdl_info_HA, \@mdl_results_AAH, $ofile_info_HH{"FH"});
    }
  }
}
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
my $combined_model_seqname_maxlen = 
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
ftr_results_add_b5e_errors(\%ftr_info_HA, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, 
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
wrapper_esl_epn_translate_startstop($execs_H{"esl-epn-translate"}, "predicted", \%ftr_info_HA, (@specstart_AA ? \@specstart_AA : undef), \%err_info_HA, \@err_ftr_instances_AHH, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

####################################################################
# Step 12. For all features that have a valid start but no in-frame
#          stop codon, look for a downstream in-frame stop
####################################################################
# TODO: make a function that performs this step

my $ftr_idx;  # counter over features
for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
  foreach my $seq_name (keys %{$err_ftr_instances_AHH[$ftr_idx]{"ext"}}) { 
    my $seq_idx = $seq_name_index_H{$seq_name};
    my $seq_len = $seq_info_HA{"len"}[$seq_idx]; 

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
    #printf("ext error for $seq_name ftr_idx: $ftr_idx %s mdl_idx: $mdl_idx\n", $ftr_info_HA{"out_tiny"}[$ftr_idx]);

    my $cur_start     = $mdl_results_AAH[$mdl_idx][$seq_idx]{"p_start"};
    my $cur_stop      = $mdl_results_AAH[$mdl_idx][$seq_idx]{"p_stop"};
    my $cur_strand    = $mdl_results_AAH[$mdl_idx][$seq_idx]{"p_strand"};
    my $cur_cumlen    = $mdl_results_AAH[$mdl_idx][$seq_idx]{"cumlen"};
    my $offset        = ($cur_cumlen % 3); 
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

    # check for a special case
    my $at_end_of_seq = 0; # set to 1 if we're at the end of the sequence
    if($cur_strand eq "+") { 
      if($posn_to_start > $seq_len) { 
        if(($posn_to_start-1) == $seq_len) { 
          # a special case
          $at_end_of_seq = 1;
        }
        else {
          DNAORG_FAIL("ERROR when looking for inframe stop for sequence $seq_name, trying to start search at position $posn_to_start but length of sequence is $seq_len", 1, $ofile_info_HH{"FH"}); 
        }
      }
    }
    else { # cur_strand eq "-" 
      if($posn_to_start <= 0) { 
        if(($posn_to_start+1) == 1) { 
          # a special case
          $at_end_of_seq = 1;
        }
        else { 
          DNAORG_FAIL("ERROR when looking for inframe stop for sequence $seq_name, trying to start search at position $posn_to_start but length of sequence is $seq_len", 1, $ofile_info_HH{"FH"}); 
        }
      }
    }
    
    my $ext_corr_stop  = undef; # third position of the next in-frame stop beyond where the stop was expected
    my $ext_stop_codon = undef; # which stop codon is used as the next in-frame stop 

    if($at_end_of_seq) { 
      $ext_corr_stop = 0;
    }
    else { 
      ($ext_corr_stop, $ext_stop_codon) = check_for_downstream_stop($sqfile, $seq_name, $posn_to_start, $seq_len, $cur_strand, $ofile_info_HH{"FH"});
    }

    if($ext_corr_stop == 0) { 
      # no stop found: nst error

      my $updated_nst_errmsg = "";
      if($at_end_of_seq) { 
        my (undef, $out_stop) = create_output_start_and_stop($cur_start, $cur_stop, $seq_len, $ofile_info_HH{"FH"});
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
mdl_results_add_str_nop_ost_lsc_dup_b3e_b3u_errors($sqfile, \%ftr_info_HA, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, 
                                                   \@err_ftr_instances_AHH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"});

# report b3e errors for cds-mp multifeature features
ftr_results_add_b3e_errors(\%ftr_info_HA, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, 
                           \@err_ftr_instances_AHH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"});

# calculate out_start, out_stop and out_stop_codon values, we need to know some of these before we call ftr_results_calculate()
mdl_results_calculate_out_starts_and_stops($sqfile, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, $ofile_info_HH{"FH"});

# set most of the multi-feature (e.g. cds-mp) errors
ftr_results_calculate($sqfile, \%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@ftr_results_AAH, \@mdl_results_AAH,
                      \%cds_tbl_HHA, \@err_ftr_instances_AHH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"});


# recalculate out_start and out_stop and out_stop_codon values, they may have changed for final mature peptides in each CDS
mdl_results_calculate_out_starts_and_stops($sqfile, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, $ofile_info_HH{"FH"});

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


##########################################################################################################
# Step 19. Run BLASTX: all full length sequences and all corrected nucleotide features versus all proteins
##########################################################################################################
$start_secs = outputProgressPrior("Running and parsing BLASTX", $progress_w, $log_FH, *STDOUT);
my $seq_nodesc_file = $ofile_info_HH{"fullpath"}{"fastanodesc"};
run_blastx_and_summarize_output(\%execs_H, $out_root, $seq_nodesc_file, $build_root, \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);
outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

# calculate the blastx related information
ftr_results_calculate_blastx($ofile_info_HH{"fullpath"}{"blastx-summary"}, \%ftr_info_HA, \%seq_info_HA, \%seq_name_index_H, \@ftr_results_AAH, \%opt_HH, \%ofile_info_HH);

# add xi* and mxi errors
openAndAddFileToOutputInfo(\%ofile_info_HH, "blasttbl", $out_root . ".blastx.tbl", 1, "information on blast and CM hits for CDS features in tabular format");
my $tmp_len = maxLengthScalarValueInArray($seq_info_HA{"seq_name"});
if($tmp_len > $combined_model_seqname_maxlen) { 
  $combined_model_seqname_maxlen = $tmp_len; 
}
ftr_results_add_blastx_errors($ofile_info_HH{"FH"}{"blasttbl"}, $combined_model_seqname_maxlen, 
                              \%ftr_info_HA, \%seq_info_HA, \@ftr_results_AAH, \@mdl_results_AAH, 
                              \@err_ftr_instances_AHH, \%err_info_HA, \%opt_HH, $ofile_info_HH{"FH"});

####################################################################
# Step 20. Add zft errors for sequences with zero annotated features
####################################################################
# add per-sequence 'zft' errors (zero annotated features)
add_zft_errors(\@err_ftr_instances_AHH, \%err_seq_instances_HH, \%ftr_info_HA, \%seq_info_HA, \%err_info_HA, \@mdl_results_AAH, \%opt_HH, \%ofile_info_HH);

###############################################################
# Step 21. Translate corrected feature sequences into proteins
###############################################################
if(! opt_Get("--skiptranslate", \%opt_HH)) { 
  $start_secs = outputProgressPrior("Translating corrected nucleotide features into protein sequences", $progress_w, $log_FH, *STDOUT);
  translate_feature_sequences($execs_H{"esl-epn-translate"}, "corrected", "corrected.translated", (@specstart_AA ? \@specstart_AA : undef), \%ftr_info_HA, \%ofile_info_HH);
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

################################################################################
# Step 22. (OPTIONAL) Create multiple alignments of DNA sequences, if --doalign
################################################################################
if(opt_Get("--doalign", \%opt_HH)) { 
  $step_desc = "Aligning and parsing corrected nucleotide hits";
  $start_secs = outputProgressPrior($step_desc, $progress_w, $log_FH, *STDOUT);
  align_hits(\%execs_H, \%mdl_info_HA, \%seq_info_HA, \@mdl_results_AAH, \%opt_HH, \%ofile_info_HH); 
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

####################################################################################
# Step 23. (OPTIONAL) Create multiple alignments of protein sequences, if --doalign
####################################################################################
if(opt_Get("--doalign", \%opt_HH)) { 
  $start_secs = outputProgressPrior("Aligning translated protein sequences", $progress_w, $log_FH, *STDOUT);
  align_protein_sequences(\%execs_H, "corrected.translated", \%ftr_info_HA, \%opt_HH, \%ofile_info_HH);
  outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

#########################################
# Step 24. Output annotations and errors
#########################################
# open files for writing
openAndAddFileToOutputInfo(\%ofile_info_HH, "tbl",            $out_root . ".tbl",               1, "All annotations in tabular format");
openAndAddFileToOutputInfo(\%ofile_info_HH, "tblsum",         $out_root . ".tbl.summary",       1, "Summary of all annotations");
openAndAddFileToOutputInfo(\%ofile_info_HH, "failtbl",        $out_root . ".fail.tbl",          1, "Annotations for all sequences with >= 1 failure in tabular format");
openAndAddFileToOutputInfo(\%ofile_info_HH, "errtbl",         $out_root . ".error.tbl",         1, "Annotations for all sequences with >= 1 error in tabular format");
openAndAddFileToOutputInfo(\%ofile_info_HH, "pererr",         $out_root . ".peraccn.errors",    1, "List of errors, one line per sequence");
openAndAddFileToOutputInfo(\%ofile_info_HH, "allerr",         $out_root . ".all.errors",        1, "List of errors, one line per error");
openAndAddFileToOutputInfo(\%ofile_info_HH, "errsum",         $out_root . ".errors.summary",    1, "Summary of all errors");
openAndAddFileToOutputInfo(\%ofile_info_HH, "pass_ftbl",      $out_root . ".ap.sqtable",        1, "Sequin feature table output for passing sequences");
openAndAddFileToOutputInfo(\%ofile_info_HH, "fail_ftbl",      $out_root . ".af.sqtable",        1, "Sequin feature table output for failing sequences (minimal)");
openAndAddFileToOutputInfo(\%ofile_info_HH, "long_ftbl",      $out_root . ".long.sqtable",      1, "Sequin feature table output for failing sequences (verbose)");
openAndAddFileToOutputInfo(\%ofile_info_HH, "pass_list",      $out_root . ".ap.seqlist",        1, "list of passing sequences");
openAndAddFileToOutputInfo(\%ofile_info_HH, "fail_list",      $out_root . ".af.seqlist",        1, "list of failing sequences");
openAndAddFileToOutputInfo(\%ofile_info_HH, "errors_list",    $out_root . ".errlist",           1, "list of errors in the sequence tables");
if(opt_IsUsed("--classerrors", \%opt_HH)) { 
  openAndAddFileToOutputInfo(\%ofile_info_HH, "fail_co_list",   $out_root . ".af-co.seqlist",   1, "list of failing sequences that would have passed if not for a classification error");
}

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
output_tbl_all_sequences(\%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \@mdl_results_AAH, \@ftr_results_AAH, \%opt_HH, \%ofile_info_HH);

# output the explanatory text
output_tbl_explanations(\@out_header_exp_A, \%ofile_info_HH);

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

######################
# feature table file #
######################
$start_secs = outputProgressPrior("Generating feature table output", $progress_w, $log_FH, *STDOUT);

my $npass = output_feature_tbl_all_sequences(\@err_ftr_instances_AHH, \%err_seq_instances_HH, \%mdl_info_HA, \%ftr_info_HA, \%seq_info_HA, \%err_info_HA, \@mdl_results_AAH, \@ftr_results_AAH, (($do_class_errors) ? \%class_errors_per_seq_H : undef), \%opt_HH, \%ofile_info_HH);

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

###########################################
# brief summary of annotations and errors #
###########################################
# OLD WAY, where we cared about how many sequences had errors:
#outputString($log_FH, 1, sprintf("#\n# Annotated %d accessions, %d (%6.4f fraction) had at least one annotation 'failure', see %s for all details.\n", 
#                                 $nseq, $nfail, ($nfail/$nseq), $ofile_info_HH{"nodirpath"}{"tbl"}));

outputString($log_FH, 1, sprintf("#\n# Annotated %d accessions:\n# %6d PASS (%5.3f) listed in $out_root.ap.seqlist\n# %6d FAIL (%5.3f) listed in $out_root.af.seqlist\n", $nseq, $npass, $npass/$nseq, ($nseq-$npass), ($nseq-$npass)/$nseq));
output_errors_summary($ofile_info_HH{"FH"}{"errsum"}, \@err_ftr_instances_AHH, \%err_seq_instances_HH, \%ftr_info_HA, \%seq_info_HA, \%err_info_HA, 0, \%opt_HH, \%ofile_info_HH); # 1: output to stdout

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
#    mdl_results_add_str_nop_ost_lsc_dup_b3e_b3u_errors()
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
#  Subroutines related to parsing cmalign output
#    parse_cmalign_ifile()
#
#################################################################
# Subroutine : parse_cmalign_ifile()
# Incept:      EPN, Thu Jan 31 13:06:54 2019
#
# Purpose:    Parse Infernal 1.1 cmalign --ifile output and store
#             results in %{$seq_info_HR}.
#
# Arguments: 
#  $ifile_file:        ifile file to parse
#  $seq_name_index_HR: REF to hash of arrays with sequence index information in seq_info_HAR, PRE-FILLED
#                      seq_name_index_H{$seq_name} = <n>, means that $seq_name is the <n>th sequence name 
#                      in @{$seq_info_HAR{*}} arrays
#  $seq_info_HAR:      REF to hash of arrays with sequence information, PRE-FILLED
#  $FH_HR:             REF to hash of file handles
#
# Returns:    void
#
# Dies:       if we find a hit to a model or sequence that we don't
#             have stored in $mdl_info_HAR or $seq_name_AR
#
################################################################# 
sub parse_cmalign_ifile { 
  my $sub_name = "parse_cmalign_ifile()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($ifile_file, $seq_name_index_HR, $seq_info_HAR, $FH_HR) = @_;
  
  # initialize the arrays
  my $nseq = scalar(keys %{$seq_name_index_HR});
  @{$seq_info_HAR->{"ifile_spos"}} = ();
  @{$seq_info_HAR->{"ifile_epos"}} = ();
  @{$seq_info_HAR->{"ifile_ins"}}  = ();
  for(my $i = 0; $i < $nseq; $i++) { 
    $seq_info_HAR->{"ifile_spos"}[$i] = -1;
    $seq_info_HAR->{"ifile_epos"}[$i] = -1;
    $seq_info_HAR->{"ifile_ins"}[$i]  = "";
  }

  open(IN, $ifile_file) || fileOpenFailure($ifile_file, $sub_name, $!, "reading", $FH_HR);

  my $line_ctr = 0;  # counts lines in ifile_file
  while(my $line = <IN>) { 
    $line_ctr++;
    if($line !~ m/^\#/) { 
      chomp $line;
      if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
      # 2 types of lines, those with 2 tokens, and those with 4 or more tokens
      #norovirus.NC_039477 7567
      #gi|669176088|gb|KM198574.1| 7431 17 7447  2560 2539 3  2583 2565 3
      my @el_A = split(/\s+/, $line);
      if   (scalar(@el_A) == 2) { ; } # ignore these lines
      elsif(scalar(@el_A) >= 4) { 
        my $nel = scalar(@el_A); 
        if((($nel - 4) % 3) != 0) { # check number of elements makes sense
          DNAORG_FAIL("ERROR in $sub_name, unexpected number of elements ($nel) in ifile line in $ifile_file on line $line_ctr:\n$line\n", 1, $FH_HR);
        }          
        my ($seqname, $seq_len, $spos, $epos) = ($el_A[0], $el_A[1], $el_A[2], $el_A[3]);
        if(! exists $seq_name_index_HR->{$seqname}) { 
          DNAORG_FAIL("ERROR in $sub_name, do not have information for sequence $seqname read in $ifile_file on line $line_ctr", 1, $FH_HR);
        }
        my $seqidx = $seq_name_index_HR->{$seqname}; # sequence index for the hit in results_AAH (2nd dim of results_AAH)
        if(! exists $seq_info_HAR->{"len"}[$seqidx]) { 
          DNAORG_FAIL(sprintf("ERROR in $sub_name, do not have length information for sequence $seqname, accession %s", $seq_info_HAR->{"accn_name"}[$seqidx]), 1, $FH_HR);
        }
        if($seq_len != $seq_info_HAR->{"len"}[$seqidx]) { 
          DNAORG_FAIL(sprintf("ERROR in $sub_name, conflicting length information for sequence $seqname in ifile, accession %s", $seq_info_HAR->{"accn_name"}[$seqidx]), 1, $FH_HR);
        }          

        # create the insert string
        my $insert_str = "";
        for(my $el_idx = 4; $el_idx < scalar(@el_A); $el_idx += 3) { 
          $insert_str .= $el_A[$el_idx] . ":" . $el_A[$el_idx+1] . ":" . $el_A[$el_idx+2] . ";"; 
        }
        
        $seq_info_HAR->{"ifile_spos"}[$seqidx] = $spos;
        $seq_info_HAR->{"ifile_epos"}[$seqidx] = $epos;
        $seq_info_HAR->{"ifile_ins"}[$seqidx]  = $insert_str;
      }
    }
  }
  close(IN);
  
  return;
}

#################################################################
# Subroutine : parse_cmalign_stk()
# Incept:      EPN, Thu Jan 31 13:06:54 2019
#
# Purpose:    Parse Infernal 1.1 cmalign stockholm alignment file
#             and store results in @{$mdl_results_AAHR}. 
#             and %{$seq_file_inserts_HAR}.
#
# Arguments: 
#  $stk_file:         stockholm alignment file to parse
#  $seq_name_index_HR:   REF to hash of arrays with sequence index information in seq_info_HAR, PRE-FILLED
#                     seq_name_index_H{$seq_name} = <n>, means that $seq_name is the <n>th sequence name 
#                     in @{$seq_info_HAR{*}} arrays
#  $seq_info_HAR:     REF to hash of arrays with sequence information, PRE-FILLED
#  $mdl_info_HAR:     REF to hash of arrays with information on the models, PRE-FILLED
#  $mdl_results_AAHR: REF to results AAH, FILLED HERE
#  $FH_HR:            REF to hash of file handles
#
# Returns:    void
#
# Dies:
#
################################################################# 
sub parse_cmalign_stk { 
  my $sub_name = "parse_cmalign_stk()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($stk_file, $seq_name_index_HR, $seq_info_HAR, $mdl_info_HAR, $mdl_results_AAHR, $FH_HR) = @_;

  # create an ESL_MSA object from the alignment
  # open and validate file
  my $msa = Bio::Easel::MSA->new({
    fileLocation => $stk_file,
    isDna => 1});  

  # build a map of aligned positions to model RF positions and vice versa, only need to do this once per alignment
  my $alen = $msa->alen;
  my @rf2a_A = (); # [1..$rfpos..$rflen] = $apos;  rf position $rfpos maps to alignment position $apos [1..$alen]  ($rf2a_A[0] = -1  (dummy value))
  $rf2a_A[0] = -1; 
  my $rf_str = $msa->get_rf;
  my @rf_A = split("", $rf_str);
  if(scalar(@rf_A) != $alen) { 
    DNAORG_FAIL(sprintf("ERROR in $sub_name, unexpected alignment length mismatch $alen != %d\n", scalar(@rf_A)), 1, $FH_HR);
  }
  my $rfpos = 0; # nongap RF (model) position [1..$rflen]
  my $apos  = 0; # alignment position [1..$alen]
  for($apos = 1; $apos <= $alen; $apos++) { 
    if($rf_A[($apos-1)] ne ".") { 
      # nongap RF (model position)
      $rfpos++;
      $rf2a_A[$rfpos] = $apos;
    }
  }
  my $rflen = $rfpos;

  # move through each sequence in the alignment and determine its boundaries for each model region
  my $nseq = $msa->nseq; 
  my $nmdl = scalar(@{$mdl_results_AAHR});
  # for each sequence, go through all models and fill in the start and stop (unaligned seq) positions
  for(my $i = 0; $i < $nseq; $i++) { 
    my $seqname = $msa->get_sqname($i);
    if(! exists $seq_name_index_HR->{$seqname}) { 
      DNAORG_FAIL("ERROR in $sub_name, do not have information for sequence $seqname from alignment in $stk_file", 1, $FH_HR);
    }
    my $seqidx = $seq_name_index_HR->{$seqname}; # sequence index for the hit in results_AAH (2nd dim of results_AAH)
    if(! exists $seq_info_HAR->{"len"}[$seqidx]) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, do not have length information for sequence $seqname, accession %s", $seq_info_HAR->{"accn_name"}[$seqidx]), 1, $FH_HR);
    }
    my $seq_len = $seq_info_HAR->{"len"}[$seqidx]; # sequence length

    # fill sequence-specific arrays
    # insert info from seq_info_HAR (read previously from cmalign --ifile output)
    my @rf2ipos_A = (); # [0..$rfpos..$rflen] = $uapos; rf position $rfpos has an insert immediately after it start at unaligned position $uapos [1..$seq_len], 0 for none
    my @rf2ilen_A = (); # [0..$rfpos..$rflen] = $ilen;  rf position $rfpos has an insert immediately after it of length $ilen, 0 for none
    # fill insert arrays, need to do this once per sequence
    # initialize insert arrays
    for($rfpos = 0; $rfpos <= $rflen; $rfpos++) { 
      $rf2ipos_A[$rfpos] = -1;
      $rf2ilen_A[$rfpos] = -1;
    }
    if($seq_info_HAR->{"ifile_ins"}[$seqidx] ne "") { 
      my @ins_A = split(";", $seq_info_HAR->{"ifile_ins"}[$seqidx]); 
      foreach my $ins_tok (@ins_A) { 
        #printf("ins_tok: $ins_tok\n");
        if($ins_tok =~ /^(\d+)\:(\d+)\:(\d+)$/) { 
          my ($i_rfpos, $i_uapos, $i_len) = ($1, $2, $3);
          $rf2ipos_A[$i_rfpos] = $i_uapos;
          $rf2ilen_A[$i_rfpos] = $i_len;
          #printf("rf2ipos_A[%5d]: %5d rf2ilen_A[%5d]: %5d\n", $i_rfpos, $i_uapos, $i_rfpos, $i_len);
        }
        else { 
          DNAORG_FAIL("ERROR in $sub_name, failed to parse insert information read from ifile for $seqname:\n" . $seq_info_HAR->{"ifile_ins"}[$seqidx], 1, $FH_HR);
        }
      }
    }    

    # alignment info, filled once per sequence then used for all model spans
    my @min_rfpos_after_A  = (); # [0..$rfpos..rflen+1]: $rfpos2 is minimum rfpos >= $rfpos that is not a gap 
                                 #                       OR is a gap but has an insert *after* it (before $rfpos2+1)
                                 #                       possible values are [-1,1..rflen] (-1 if none) with two *exceptions*:
                                 #                       element [0]       is special: can be 0 if inserts exist after position 0 (before position 1)
                                 #                       element [rflen+1] is special: always -1
    my @max_rfpos_before_A = (); # [0..$rfpos..rflen+1]: $rfpos2 is maximum rfpos >= $rfpos that is not a gap 
                                 #                       OR is a gap but has an insert *before* it (after $rfpos2-1)
                                 #                       possible values are [-1,1..rflen] (-1 if none) with two *exceptions*
                                 #                       element [0]       is special: always -1
                                 #                       element [rflen+1] is special: can be rfpos+1 if inserts exist after rflen (before position rflen+1)

    my @min_uapos_after_A  = (); # [0..$rfpos..rflen+1]: minimum unaligned position for current sequence 
                                 #                       that aligns at or inserts *after*  $min_rfpos_after_A[$rfpos]
                                 #                       -1 if $min_rfpos_after_A[$rfpos]  == -1
    my @max_uapos_before_A = (); # [0..$rfpos..rflen+1]: maximum unaligned position for current sequence 
                                 #                       that aligns at or inserts *before* $max_rfpos_before_A[$rfpos]
                                 #                       -1 if $max_rfpos_before_A[$rfpos] == -1

    my @rfpos_pp_A = ();         # [0..$rfpos..rflen+1]: posterior probability character for current sequence at RF position $rfpos
                                 #                       '.' if sequence is a gap at that RF position $rfpos
                                 #                       special values: $rfpos_pp_A[0] = -1, $rfpos_pp_A[$rflen+1] = -1

    $rfpos = 0;    # model positions (nongap RF position)
    my $uapos = 0; # unaligned sequence position (position in actual sequence)
    # initialize
    for($rfpos = 0; $rfpos <= ($rflen+1); $rfpos++) { 
      $min_rfpos_after_A[$rfpos]  = -1;
      $max_rfpos_before_A[$rfpos] = -1;
      $min_uapos_after_A[$rfpos]  = -1;
      $max_uapos_before_A[$rfpos] = -1;
      $rfpos_pp_A[$rfpos]         = ".";
    }
    # get aligned sequence, length will be alen
    my $sqstring_aligned = $msa->get_sqstring_aligned($i);
    my $ppstring_aligned = $msa->get_ppstring_aligned($i);
    if(length($sqstring_aligned) != $alen) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, fetched aligned seqstring of unexpected length (%d, not %d)\n$sqstring_aligned\n", length($sqstring_aligned), $alen), 1, $FH_HR);
    }
    if(length($ppstring_aligned) != $alen) { 
      DNAORG_FAIL(sprintf("ERROR in $sub_name, fetched aligned posterior probability string of unexpected length (%d, not %d)\n$sqstring_aligned\n", length($ppstring_aligned), $alen), 1, $FH_HR);
    }
    my @sq_A = split("", $sqstring_aligned);
    my @pp_A = split("", $ppstring_aligned);
    #printf("sq_A size: %d\n", scalar(@sq_A));

    # first pass, from right to left to fill $min_**pos_after arrays, and rf
    my $min_rfpos = -1;
    my $min_uapos = $seq_len+1;
    for($rfpos = $rflen; $rfpos >= 0; $rfpos--) { 
      # is this position a gap for this sequence in the alignment? 
      $apos = $rf2a_A[$rfpos];
      #printf("rfpos: %5d  apos: %5d  min_rfpos: %5d  min_uapos: %5d\n", $rfpos, $apos, $min_rfpos, $min_uapos);
      if(($rfpos > 0) && ($sq_A[($apos-1)] ne ".") && ($sq_A[($apos-1)] ne "-")) { 
        # sequence is NOT a gap at $rfpos
        $min_rfpos = $rfpos;
        $min_uapos--;
        $min_rfpos_after_A[$rfpos] = $min_rfpos;
        $min_uapos_after_A[$rfpos] = $min_uapos;
        if($rf2ipos_A[$rfpos] != -1) { 
          $min_uapos -= $rf2ilen_A[$rfpos]; # subtract inserts
        }
        $rfpos_pp_A[$rfpos] = $pp_A[($apos-1)];
        #printf("rfpos: $rfpos apos: $apos set rfpos_pp_A[$rfpos] to " . $pp_A[($apos-1)] . "\n");
      }
      elsif($rf2ipos_A[$rfpos] != -1) { 
        # sequence is a gap at $rfpos but has an insert AFTER it
        $min_rfpos = $rfpos;
        $min_uapos -= $rf2ilen_A[$rfpos]; # subtract inserts
        $min_rfpos_after_A[$rfpos] = $min_rfpos;
        $min_uapos_after_A[$rfpos] = $min_uapos;
      }
      elsif($min_rfpos != -1) { 
        # sequence is a gap at $rfpos and has no insert AFTER it
        # but we have seen at least 1 nucleotide thus far (at rfpos $min_rfpos)
        $min_rfpos_after_A[$rfpos] = $min_rfpos;
        $min_uapos_after_A[$rfpos] = $min_uapos;
      }
      # else we haven't seen any nucleotides yet, leave values as -1 
    }
    if($min_uapos != 1) { 
      DNAORG_FAIL("ERROR in $sub_name, failed to account for all nucleotides when parsing alignment for $seqname, pass 1 (min_uapos should be 1 but it is $min_uapos)", 1, $FH_HR);
    }      

    # second pass, from left to right to fill $max_**pos_before arrays:
    my $max_rfpos = -1;
    my $max_uapos = 0;
    for($rfpos = 1; $rfpos <= ($rflen+1); $rfpos++) { 
      # is this position a gap for this sequence in the alignment? 
      $apos = $rf2a_A[$rfpos];
      if(($rfpos <= $rflen) && ($sq_A[($apos-1)] ne ".") && ($sq_A[($apos-1)] ne "-")) { 
        # sequence is NOT a gap at $rfpos
        $max_rfpos = $rfpos;
        $max_uapos++;
        $max_rfpos_before_A[$rfpos] = $max_rfpos;
        $max_uapos_before_A[$rfpos] = $max_uapos;
        if($rf2ipos_A[($rfpos-1)] != -1) { 
          $max_uapos += $rf2ilen_A[($rfpos-1)]; # add inserts
        }          
      }
      elsif($rf2ipos_A[($rfpos-1)] != -1) { 
        # sequence is a gap at $rfpos but has an insert BEFORE it
        $max_rfpos  = $rfpos;
        $max_uapos += $rf2ilen_A[($rfpos-1)]; # add inserts
        $max_rfpos_before_A[$rfpos] = $max_rfpos;
        $max_uapos_before_A[$rfpos] = $max_uapos;
      }
      elsif($max_rfpos != -1) { 
        # sequence is a gap at $rfpos and has no insert BEFORE it
        # but we have seen at least 1 nucleotide thus far (at rfpos $max_rfpos)
        $max_rfpos_before_A[$rfpos] = $max_rfpos;
        $max_uapos_before_A[$rfpos] = $max_uapos;
      }
      # else we haven't seen any nucleotides yet, leave values as -1 
    }
    if($max_uapos != $seq_len) { 
      DNAORG_FAIL("ERROR in $sub_name, failed to account for all nucleotides when parsing alignment for $seqname, pass 2 (max_uapos should be $seq_len but it is $max_uapos)", 1, $FH_HR);
    }      
    
    # DEBUG PRINT
#    printf("***************************************************\n");
#    printf("DEBUG print $seqname\n");
#    for($rfpos = 0; $rfpos <= ($rflen+1); $rfpos++) { 
#      printf("rfpos[%5d] min_rf_after_A: %5d  min_ua_after_A: %5d  max_rf_before_A: %5d  max_rf_before_A: %5d\n", 
#             $rfpos, 
#             $min_rfpos_after_A[$rfpos],
#             $min_uapos_after_A[$rfpos],
#             $max_rfpos_before_A[$rfpos],
#             $max_uapos_before_A[$rfpos]);
#    }
#    printf("***************************************************\n");

    # given model span s..e
    # if strand eq "+"
    #   if C[rfpos] > D[rfpos] then no hit (A[rfpos] should be > B[rfpos])
    #   else (C[rfpos] <= D[rfpos]) 
    #        hit uaseq span is from C[rfpos] to D[rfpos]
    #        hit rf span is from A[rfpos] to B[rfpos]

    # now we have all the info we need for this sequence to determine sequence boundaries for each model region
    for(my $m = 0; $m < $nmdl; $m++) { 
      my $mdl_start_rfpos = $mdl_info_HAR->{"ref_start"}[$m];
      my $mdl_stop_rfpos  = $mdl_info_HAR->{"ref_stop"}[$m];
      my $mdl_strand      = $mdl_info_HAR->{"ref_strand"}[$m];

#      printf("model $m $mdl_start_rfpos..$mdl_stop_rfpos\n");
#      $rfpos = $mdl_start_rfpos;
#      printf("\trfpos[%5d] min_rf_after_A: %5d  min_ua_after_A: %5d  max_rf_before_A: %5d  max_ua_before_A: %5d\n", 
#             $rfpos, 
#             $min_rfpos_after_A[$rfpos],
#             $min_uapos_after_A[$rfpos],
#             $max_rfpos_before_A[$rfpos],
#             $max_uapos_before_A[$rfpos]);
#      $rfpos = $mdl_stop_rfpos;
#      printf("\trfpos[%5d] min_rf_after_A: %5d  min_ua_after_A: %5d  max_rf_before_A: %5d  max_ua_before_A: %5d\n", 
#             $rfpos, 
#             $min_rfpos_after_A[$rfpos],
#             $min_uapos_after_A[$rfpos],
#             $max_rfpos_before_A[$rfpos],
#             $max_uapos_before_A[$rfpos]);

      my $start_rfpos = -1; # model position of start of this model region for this aligned sequence, stays at -1 if none
      my $stop_rfpos  = -1; # model position of stop  of this model region for this aligned sequence, stays at -1 if none
      my $start_uapos = -1; # unaligned position of start of this model region for this aligned sequence, stays at -1 if none
      my $stop_uapos  = -1; # unaligned position of stop  of this model region for this aligned sequence, stays at -1 if none
      my $p_5seqflush = undef;
      my $p_3seqflush = undef;

      # this should work regardless of strand
      if(($min_rfpos_after_A[$mdl_start_rfpos] != -1) && 
         ($max_rfpos_before_A[$mdl_stop_rfpos] != -1)) { 

        $start_uapos = $min_uapos_after_A[$mdl_start_rfpos];
        $stop_uapos  = $max_uapos_before_A[$mdl_stop_rfpos];

        $start_rfpos = $min_rfpos_after_A[$mdl_start_rfpos];
        $stop_rfpos  = $max_rfpos_before_A[$mdl_stop_rfpos];

        if($mdl_strand eq "+") { 
          $p_5seqflush = ($start_uapos == 1)        ? 1 : 0;
          $p_3seqflush = ($stop_uapos  == $seq_len) ? 1 : 0;
        }
        else { 
          $p_5seqflush = ($start_uapos == $seq_len) ? 1 : 0;
          $p_3seqflush = ($stop_uapos  == 1)        ? 1 : 0;
        }

        %{$mdl_results_AAHR->[$m][$seqidx]} = ();
        $mdl_results_AAHR->[$m][$seqidx]{"p_start"}     = $start_uapos;
        $mdl_results_AAHR->[$m][$seqidx]{"p_stop"}      = $stop_uapos;
        $mdl_results_AAHR->[$m][$seqidx]{"p_strand"}    = $mdl_strand;
        $mdl_results_AAHR->[$m][$seqidx]{"p_5overhang"} = $start_rfpos - $mdl_start_rfpos;
        $mdl_results_AAHR->[$m][$seqidx]{"p_3overhang"} = $mdl_stop_rfpos - $stop_rfpos;
        $mdl_results_AAHR->[$m][$seqidx]{"p_5seqflush"} = $p_5seqflush;
        $mdl_results_AAHR->[$m][$seqidx]{"p_3seqflush"} = $p_3seqflush;
        $mdl_results_AAHR->[$m][$seqidx]{"p_nhits"}     = 1;
        $mdl_results_AAHR->[$m][$seqidx]{"p_startgap"}  = ($rfpos_pp_A[$mdl_start_rfpos] eq ".") ? 1  : 0;
        $mdl_results_AAHR->[$m][$seqidx]{"p_stopgap"}   = ($rfpos_pp_A[$mdl_stop_rfpos]  eq ".") ? 1  : 0;
        $mdl_results_AAHR->[$m][$seqidx]{"p_startpp"}   = ($rfpos_pp_A[$mdl_start_rfpos] eq ".") ? -1 : convert_pp_char_to_pp_avg($rfpos_pp_A[$mdl_start_rfpos], $FH_HR);
        $mdl_results_AAHR->[$m][$seqidx]{"p_stoppp"}    = ($rfpos_pp_A[$mdl_stop_rfpos]  eq ".") ? -1 : convert_pp_char_to_pp_avg($rfpos_pp_A[$mdl_stop_rfpos], $FH_HR);

        printf("model: $mdl_start_rfpos to $mdl_stop_rfpos\n");
        foreach my $key ("p_start", "p_stop", "p_strand", "p_5overhang", "p_3overhang", "p_5seqflush", "p_3seqflush", "p_startgap", "p_stopgap", "p_startpp", "p_stoppp") { 
          printf("stored $m $seqidx $key $mdl_results_AAHR->[$m][$seqidx]{$key}\n");
        }
      }
    } # end of 'for(my $m = 0; $m < $nmdl; $m++)'
  } # end of 'for(my $i = 0; $m < $nseq; $m++)'
  undef $msa;

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
            if(($append_start <= $seq_info_HAR->{"len"}[$seq_idx]) && 
               ($append_stop  <= $seq_info_HAR->{"len"}[$seq_idx]) && 
               ($append_start >= 1) && 
               ($append_stop  >= 1)) { 
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
# Returns:    Length of longest sequence name in output file
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

  my $ret_seqname_maxlen = 0;

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

          my $cur_seqname_maxlen = combine_sequences(\@tmp_hit_fafile_A, $seq_name_AR, $ftr_hit_fafile, 
                            0, # do not require sequences be in a contiguous subset of files beginning with the first one to be combined, 
                               # allow any contiguous subset of files
                            $opt_HHR, $ofile_info_HHR->{"FH"});
          
          my $ofile_info_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $ftr_info_file_key, $ofile_info_HHR->{"FH"});
          addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_key, $ftr_hit_fafile, 0, "fasta file with $out_key hits for feature " . $ftr_info_HAR->{"out_tiny"}[$ftr_idx] . " from " . $ftr_info_HAR->{"nmodels"}[$ftr_idx] . " combined model predictions");
          if($cur_seqname_maxlen > $ret_seqname_maxlen) { $ret_seqname_maxlen = $cur_seqname_maxlen; }
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

  return $ret_seqname_maxlen;
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
# Returns:    Length of longest sequence name in output file
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

  my $ret_seqname_maxlen = 0;

  my $ftr_info_file_key        = $out_key . ".hits.fa";
  my $mdl_info_append_file_key = $out_key . ".hits.append.fa";

  # printf("in $sub_name, out_key: $out_key, ftr_info_file_key: $ftr_info_file_key\n");

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") { # we only do this for features created by combining other features
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
        my $cur_seqname_maxlen = combine_sequences(\@tmp_hit_fafile_A, $seq_name_AR, $combined_ftr_hit_fafile, 
                                                   0, # do not require sequences be in a contiguous subset of files beginning with the first one to be combined, 
                                                   # allow any contiguous subset of files
                                                   $opt_HHR, $ofile_info_HHR->{"FH"}); 
        my $ofile_info_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, $ftr_info_file_key, $ofile_info_HHR->{"FH"});
        addClosedFileToOutputInfo($ofile_info_HHR, $ofile_info_key, $combined_ftr_hit_fafile, 0, "fasta file with $out_key hits for feature " . $ftr_info_HAR->{"out_tiny"}[$ftr_idx] . " from " . $ftr_info_HAR->{"nmodels"}[$ftr_idx] . " combined model predictions");
        if($cur_seqname_maxlen > $ret_seqname_maxlen) { $ret_seqname_maxlen = $cur_seqname_maxlen; }
      } # end of 'if($at_least_one_fafile)'
      else { 
        # no fasta files exist, redefine $ftr_info_HAR->{"$ftr_info_file_key"}[$ftr_idx] to 
        # /dev/null so downstream functions know that it should not exist
        $ftr_info_HAR->{$ftr_info_file_key}[$ftr_idx] = "/dev/null";
      }
    }
  } # end of 'for' loop over $ftr_idx

  return $ret_seqname_maxlen;
}


#################################################################
# Subroutine:  combine_sequences()
# Incept:      EPN, Wed Mar  2 16:11:40 2016
#
# Purpose:    Helper function for combine_model_hits() and
#             combine_feature_hits(). Given an array of fasta files,
#             each with a different subsequence from the same parent
#             sequences, create a single new fasta file that has the
#             subsequences concatenated together.  An example is
#             stitching together exons into a CDS.  Uses BioEasel's
#             sqfile module.
#
# Arguments: 
#  $indi_file_AR:  REF to array of fasta files to combine
#  $seq_name_AR:   REF to array with order of sequence names
#  $multi_file:    name of multi file to create
#  $require_first: '1' to require each sequence to keep be in a contiguous set starting with file 1
#                  '0' to require each sequence to keep be in a contiguous set starting with any file
#  $opt_HHR:      REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:        REF to hash of file handles
#
# Returns:    Length of longest sequence name created
#
# Dies:       If we have a problem reading the fasta files
#
#################################################################
sub combine_sequences {
  my $sub_name = "combine_sequences";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($indi_file_AR, $seq_name_AR, $multi_file, $require_first, $opt_HHR, $FH_HR) = @_;
  
  my @sqfile_A                      = (); # array of open Bio::Easel::SqFile objects, one per indi_file_AR element
  my @sqfile_sqname_AA              = (); # 2D array of SqFile sequence names, read from the indi files, [0..$nfiles-1][0..$nseq_in_file_f]
  my $nfiles                        = scalar(@{$indi_file_AR}); # number of input fasta files
  my @seq_name2sqfile_sqname_map_AA = (); # 2D array, 1st dim: 0..$seq_idx..$nseq-1, 2nd dim 0..$file_idx..$nfiles, value: ssi index of $seq_idx subseq in file $f
  my @seq_name_exists_AA            = (); # 2D array, 1st dim: 0..$seq_idx..$nseq-1, 2nd dim: 0..$file_idx..$nfiles-1; value '1' if $seq_idx exists in file $f
  my @seq_name_coords_A             = (); # key: sequence name from @{$seq_name_AR}, value: concatenated set of coordinates of all subseqs in all $nfiles
  my @seq_name_fetch_me_A           = (); # [0..$i..$nseq_name-1]: value is "0" if seq $i should be included in the $multi_file, '0' if not
                                          # How we decide depends on value of $require_all input variable:
                                          #
                                          # IF $require_first is '1':
                                          # we only include sequences that exist in a contiguous set of files in @sqfile_A starting with the first one
                                          #
                                          # IF $require_first is '0':
                                          # we fetch each sequence that exists only in a contiguous set of files in @sqfile_A starting with any file
                                          #
                                          # To determine this we use @seq_name_exists_AA to create a string of 0s and 1s indicating which
                                          # files it is in.
                                          # example A: 01110: sequence exists in files 2, 3, and 4 but not 1 and 5
                                          # example B: 11110: sequence exists in files 1, 2, 3, and 4 but not 5
                                          # example C: 11010: sequence exists in files 1, 2, and 4 but not 3 and 5
                                          # If $require_first is '1', only example B would be included
                                          # If $require_first is '0', examples A and B would be included
            
  my @sqname_AA                     = (); # 2D array, 1st dim [0..$file_idx..$nfiles-1], 2nd dim: 0 to number of sequences in file $file_idx

  my $ret_seqname_maxlen = 0;

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
    $seq_name_fetch_me_A[$seq_idx] = 0;
    @{$seq_name2sqfile_sqname_map_AA[$seq_idx]} = ();
    for($file_idx = 0; $file_idx < $nfiles; $file_idx++) { 
      $seq_name2sqfile_sqname_map_AA[$seq_idx][$file_idx] = -1; # updated in block below if $seq_idx exists in $file_idx
      $seq_name_exists_AA[$seq_idx][$file_idx] = 0;
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
        my ($is_nse, $seq_name, $start, $end, $str) = nseBreakdown($sqname_AA[$file_idx][$sqfile_seq_idx]);
        my $coords = $start . "-" . $end;
        if((! $is_nse) || $coords !~ m/\-/) { 
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
      }
    } # end of 'if($indi_file_AR->[$file_idx] ne "/dev/null")'
  } # end of 'for' loop over file indexes

  # update values in @seq_name_fetch_me_A based on values in @seq_name_exists_AA
  for($seq_idx = 0; $seq_idx < $nseq_name; $seq_idx++) { 
    my $seq_name_exists_str = join("", @{$seq_name_exists_AA[$seq_idx]});
    # remove trailing 0s
    $seq_name_exists_str =~ s/0+$//;
    # remove leading 0s (only if (! $require_first))
    if(! $require_first) { $seq_name_exists_str =~ s/^0+//; }
    # determine if we have a string of >=1 '1' left, with zero zeroes
    if((length($seq_name_exists_str) > 0) && ($seq_name_exists_str !~ m/0/)) { 
      # sequence exists in a contiguous subset of the sequence files
      $seq_name_fetch_me_A[$seq_idx] = 1;
    }
    else { 
      $seq_name_fetch_me_A[$seq_idx] = 0;
    }
  } 

  # now for each seq_name that we want to fetch, fetch all subsequences for that 
  # sequence from all the individual files into a new sequence in a new file ($multi_file)
  open(OUT, ">", $multi_file) || die "ERROR unable to open $multi_file for writing";
  for($seq_idx = 0; $seq_idx < $nseq_name; $seq_idx++) { 
    # printf("HEYA in $sub_name creating $multi_file fetch me for sequence %s is %d\n", $seq_name_AR->[$seq_idx], $seq_name_fetch_me_A[$seq_idx]);
    if($seq_name_fetch_me_A[$seq_idx]) { 
      my $seq_name = $seq_name_AR->[$seq_idx];
      print OUT ">" . $seq_name . "/" . $seq_name_coords_A[$seq_idx] . "\n";
      my $cur_seqname_len = length($seq_name . "/" . $seq_name_coords_A[$seq_idx]);
      if($cur_seqname_len > $ret_seqname_maxlen) { $ret_seqname_maxlen = $cur_seqname_len; }
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

  return $ret_seqname_maxlen;
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
#             Checks for and adds the following error codes for CDS features
#             (type is "cds-mp" or "cds-notmp"):
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
#                  stop downstream of the predicted stop (ext and mtr are exclusive)
#
#             "nst": "maybe" added for features with no in-frame stop, removed or
#                  updated later (not in this function) when we look for an in-frame
#                  stop downstream of the predicted stop (ext and mtr are exclusive)
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
  my $is_cds     = featureTypeIsCds($ftr_info_HAR->{"type"}[$ftr_idx]);
  my $is_matpept = featureTypeIsMaturePeptide($ftr_info_HAR->{"type"}[$ftr_idx]);
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
        # $is_cds:     '1' if current feature is a CDS, '0' if it is not (e.g. mature peptide, $xfeat, or $dfeat)
        # $is_matpept: '1' if current feature is a CDS, '0' if it is not (e.g. mature peptide, $xfeat, or $dfeat)
        #
        # Variables derived from esl-epn-translate output we're currently parsing
        # $start_is_valid:     '1' if current feature's first 3 nt encode a valid start codon
        # $stop_is_valid:      '1' if current feature's final 3 nt encode a valid stop codon and total feature length is multiple of 3
        # $early_inframe_stop: '1' if an inframe stop exists prior to predicted stop
        #
        # 7 possibilities, each with different outcome (P1-P7):
        #                                                                                         | make correction to |
        # idx | is_cds|is_mp | start_is_valid | stop_is_valid | early_inframe_stop ||   errors    | stop coordinate?   |
        # ----|--------------|----------------|---------------|--------------------||-------------|--------------------|
        #  P1 | true  |false |          false |          any  |                any ||         str |                 no |
        #  P2 | true  |false |           true |        false  |              false ||     stp ext?|     maybe (if ext) |
        #  P3 | true  |false |           true |        false  |               true ||     stp trc |                yes |
        #  P4 | true  |false |           true |         true  |              false ||        none |                 no |
        #  P5 | true  |false |           true |         true  |               true ||         trc |                yes |
        # -----------------------------------------------------------------------||-------------------------------------
        #  P6 | false | true |            any |          any  |              false ||        mtr? |                 no |
        #  P7 | false | true |            any |          any  |               true ||    trc mtr? |                yes |
        # --------------------------------------------------------------------------------------------------------------
        # 
        # in table above:
        # '?' after error code means that error is possible, we have to check for it later           
        # 'any' means that any value is possible, outcome is unaffected by value
        #
        if($is_cds) { 
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
                # printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 2 (stp, maybe ext)\n");
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
                # printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 5 (trc)\n");
              }
            }              
          }
        } # end of 'if($is_cds)'
        elsif($is_matpept) { 
          if(! $early_inframe_stop) { 
            ; 
            # possibility 6 (P6): maybe mtr error later, but can't check for it now, do nothing;
            #printf("in $sub_name, feature index $ftr_idx, seq $seq_name, possibility 6 (no error)\n");
          }
          else { # $early_inframe_stop is '1'
            # possibility 7 (P7): trc error, maybe mtr error later, but can't check for it now
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
        runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
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
#    mdl_results_add_str_nop_ost_lsc_dup_b3e_b3u_errors()
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
      my $is_matpept = featureTypeIsMaturePeptide($ftr_info_HAR->{"type"}[$ftr_idx]);
      my $is_cds     = featureTypeIsCds($ftr_info_HAR->{"type"}[$ftr_idx]);
      for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
        my $seq_name = $seq_info_HAR->{"seq_name"}[$seq_idx];

        # only proceed past here if we have a stp, trc or ext error
        if((exists $err_ftr_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name}) ||
           (exists $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name}) ||
           (exists $err_ftr_instances_AHH[$ftr_idx]{"ext"}{$seq_name})) { 

          # determine first and final model indices, we'll use these differently depending 
          # on if we're a stp, trc or ext, but we need to know these for all three cases
          $first_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx];
          $final_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$ftr_idx];
          # first_mdl_idx should be the first model for which we have a prediction
          while(! exists $mdl_results_AAHR->[$first_mdl_idx][$seq_idx]{"p_start"}) { 
            $first_mdl_idx++; 
            if($first_mdl_idx > $final_mdl_idx) { 
              DNAORG_FAIL(sprintf("ERROR in $sub_name, can't determine first model for feature $ftr_idx (%s) for sequence $seq_name. stp or trc or ext error exists but no models have predictions for this feature.", $ftr_info_HAR->{"out_tiny"}[$ftr_idx]), 1, $FH_HR);
            }
          }
          # final_mdl_idx should be the final model for which we have a prediction
          while(! exists $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"p_start"}) { 
            $final_mdl_idx--; 
            if($final_mdl_idx < $first_mdl_idx) { 
              DNAORG_FAIL(sprintf("ERROR in $sub_name, can't determine final model for feature $ftr_idx (%s) for sequence $seq_name. stp or trc or ext error exists but no models have predictions for this feature.", $ftr_info_HAR->{"out_tiny"}[$ftr_idx]), 1, $FH_HR);
            }
          }
          ###########################################
          # block that handles potential stp errors #
          ###########################################
          if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name}) { 
            if($err_ftr_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name} ne "maybe") { 
              DNAORG_FAIL(sprintf("ERROR in $sub_name, stp error with non-maybe value %s for ftr %s seq_name: $seq_name", 
                                  $err_ftr_instances_AHHR->[$ftr_idx]{"stp"}{$seq_name}, $ftr_info_HAR->{"out_short"}[$ftr_idx]),
                          1, $FH_HR);
            }
            my $stp_err_stop_codon = fetchStopCodon($sqfile, $seq_name, 
                                                    $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"p_stop"}, 
                                                    $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"p_strand"}, 1, $FH_HR); 
                                                    # '1' above: it's okay if codon we fetch is less than 3 nt due to being off end of the sequence
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
                                                                              $seq_info_HAR->{"len"}[$seq_idx], $FH_HR));
                if($ftr_info_HAR->{"nmodels"}[$ftr_idx] != 1) { 
                  $updated_trc_errmsg .= sprintf(" %s %d of %d", ($is_cds) ? "exon" : "segment", $mdl_idx - $first_mdl_idx + 1, $ftr_info_HAR->{"nmodels"}[$ftr_idx]);
                }
                $updated_trc_errmsg .= sprintf(" revised to %d..%d (stop shifted %d nt)", 
                                               create_output_start_and_stop($mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"},
                                                                            $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"c_stop"},
                                                                            $seq_info_HAR->{"len"}[$seq_idx], $FH_HR), 
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
            my $mdl_idx = $final_mdl_idx;
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
                                                                          $seq_info_HAR->{"len"}[$seq_idx], $FH_HR));
            if($ftr_info_HAR->{"nmodels"}[$ftr_idx] != 1) { 
              $updated_ext_errmsg .= sprintf(" %s %d of %d", ($is_cds) ? "exon" : "segment", $mdl_idx - $ftr_info_HAR->{"first_mdl"}[$ftr_idx] + 1, $ftr_info_HAR->{"nmodels"}[$ftr_idx]);
            }
            $updated_ext_errmsg .= sprintf(" revised to %d..%d (stop shifted %d nt)", 
                                           create_output_start_and_stop($mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"},
                                                                        $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"c_stop"},
                                                                        $seq_info_HAR->{"len"}[$seq_idx], $FH_HR), 
                                           $len_corr);
            error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "ext", $seq_info_HAR->{"seq_name"}[$seq_idx], $updated_ext_errmsg, $FH_HR);
          } # end of ext block
        } # end of if statement entered if stp, trc or ext
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

  my $do_all_olp = opt_Get("--allolp", $opt_HHR);
  my $do_all_adj = opt_Get("--allolp", $opt_HHR);
  
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
                                 $seq_info_HAR->{"len"}[$seq_idx], 
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
      if(($start_A[$mdl_idx] != -1) || # this model has a prediction (does not have a 'nop' error)
         ($do_all_adj)) {              # --alladj option used, report all aja/ajb errors for models without predictions (with 'nop')
        if($idx_ajb_str_A[$mdl_idx] ne $mdl_info_HAR->{"idx_ajb_str"}[$mdl_idx]) { 
          my @diff_A = ();
          compareTwoOverlapOrAdjacencyIndexStrings($mdl_info_HAR->{"idx_ajb_str"}[$mdl_idx], 
                                                   $idx_ajb_str_A[$mdl_idx], 
                                                   $nmdl-1,
                                                   \@diff_A, $FH_HR);
          if(! $do_all_adj) { # --alladj option not used, so we ignore any aja/ajb errors with models without predictions (with 'nop' errors)
            for(my $a = 0; $a < $nmdl; $a++) { 
              if($start_A[$a] == -1) { $diff_A[$a] = 0; } # now, in next for($i) loop diff values for model $i where $i has no prediction won't be printed
            }
          }
          for(my $i = 0; $i < $nmdl; $i++) { 
            if($diff_A[$i] != 0) { 
              $ftr_ajb_err_msg_A[$ftr_idx] .= sprintf("%s%s(%s,%s)", 
                                                      ($ftr_ajb_err_msg_A[$ftr_idx] eq "") ? "" : ",", # need to add a comma only if we're appending
                                                      ($diff_A[$i] eq "-1") ? "-" : "+",              # is it a lost or added adjacency?
                                                      $mdl_info_HAR->{"out_idx"}[$mdl_idx], $mdl_info_HAR->{"out_idx"}[$i]);
            }
          }
        }
      }
      # construct aja err message
      if(($start_A[$mdl_idx] != -1) || # this model has a prediction (does not have a 'nop' error)
         ($do_all_adj)) {              # --alladj option used, report all aja/ajb errors for models without predictions (with 'nop')
        if($idx_aja_str_A[$mdl_idx] ne $mdl_info_HAR->{"idx_aja_str"}[$mdl_idx]) { 
          my @diff_A = ();
          compareTwoOverlapOrAdjacencyIndexStrings($mdl_info_HAR->{"idx_aja_str"}[$mdl_idx], 
                                                 $idx_aja_str_A[$mdl_idx], 
                                                 $nmdl-1,
                                                 \@diff_A, $FH_HR);
          if(! $do_all_adj) { # --alladj option not used, so we ignore any aja/ajb errors with models without predictions (with 'nop' errors)
            for(my $a = 0; $a < $nmdl; $a++) { 
              if($start_A[$a] == -1) { $diff_A[$a] = 0; } # now, in next for($i) loop diff values for model $i where $i has no prediction won't be printed
            }
          }
          for(my $i = 0; $i < $nmdl; $i++) { 
            if($diff_A[$i] != 0) { 
              $ftr_aja_err_msg_A[$ftr_idx] .= sprintf("%s%s(%s,%s)", 
                                                      ($ftr_aja_err_msg_A[$ftr_idx] eq "") ? "" : ",", # need to add a comma only if we're appending
                                                      ($diff_A[$i] eq "-1") ? "-" : "+",              # is it a lost or added adjacency?
                                                      $mdl_info_HAR->{"out_idx"}[$mdl_idx], $mdl_info_HAR->{"out_idx"}[$i]);
            }
          }
        }
      }

      # construct olp err message
      if(($start_A[$mdl_idx] != -1) || # this model has a prediction (does not have a 'nop' error)
         ($do_all_olp)) {              # --allolp option used, report all olp errors for models without predictions (with 'nop')
        if($idx_olp_str_A[$mdl_idx] ne $mdl_info_HAR->{"idx_olp_str"}[$mdl_idx]) { 
          my @diff_A = ();
          compareTwoOverlapOrAdjacencyIndexStrings($mdl_info_HAR->{"idx_olp_str"}[$mdl_idx], 
                                                   $idx_olp_str_A[$mdl_idx], 
                                                   $nmdl-1,
                                                   \@diff_A, $FH_HR);
          if(! $do_all_olp) { # --allolp option not used, so we ignore any olp errors with models without predictions (with 'nop' errors)
            for(my $a = 0; $a < $nmdl; $a++) { 
              if($start_A[$a] == -1) { $diff_A[$a] = 0; } # now, in next for($i) loop diff values for model $i where $i has no prediction won't be printed
            }
          }
          for(my $i = 0; $i < $nmdl; $i++) { 
            if($diff_A[$i] != 0) { 
            $ftr_olp_err_msg_A[$ftr_idx] .= sprintf("%s%s(%s,%s)", 
                                                    ($ftr_olp_err_msg_A[$ftr_idx] eq "") ? "" : ",", # need to add a comma only if we're appending
                                                    ($diff_A[$i] eq "-1") ? "-" : "+",              # is it a lost or added adjacency?
                                                    $mdl_info_HAR->{"out_idx"}[$mdl_idx], $mdl_info_HAR->{"out_idx"}[$i]);
            }
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
# Subroutine:  ftr_results_add_b5e_errors
# Incept:      EPN, Wed Feb 21 12:19:14 2018
#
# Purpose:    Report 'b5e' for features with 'annot_type' of 'multifeature'
#             and type of 'cds-mp.
#             Uses mdl_results to do this. If the first (5'-most) primary
#             child of a multifeature feature should have a 'b5e' error
#             (by looking at p_5overhang and p_5seqflush values in mdl_results)
#             then the parent multifeature should too.
#
#             Checks for and adds or updates the following error 
#             codes for features with "annot_type" eq "multifeature' and
#             type 'cds-mp':
#             
#             "b5e": adds this error, predicted hit of 5'-most model 
#                    not flush with model end but flush with sequence end on 5'
#                    OR predicted hit of 5'-most model *is* flush with model end
#                    but also flush with sequence end on 5' and *is not* the first
#                    model of the first child of this multifeature parent feature.
#
#
# Arguments: 
#  $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
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
################################################################# 
sub ftr_results_add_b5e_errors { 
  my $sub_name = "ftr_results_add_b5e_errors";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $mdl_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $err_ftr_instances_AHHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $ftr_idx;   # counter over features
  my $mdl_idx;   # counter over models
  my $seq_idx;   # counter over sequences
  my $seq_name;  # name of one sequence

  # foreach annot_type:multifeature and type:'cds-mp' feature, 
  # determine if the first primary child model with a prediction has a 'b5e' error, if so add a 'b5e' error to its parent cds-mp feature
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
        my $seen_hit = 0;      # set to '1' once we've seen the first model with an annotated hit
        $seq_name = $seq_info_HAR->{"seq_name"}[$seq_idx];
        my $b5e_flag = 0;

        # step through all primary children of this feature
        for(my $child_idx = 0; $child_idx < $np_children; $child_idx++) { 
          my $child_ftr_idx = $primary_children_idx_A[$child_idx];
          for(my $child_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$child_ftr_idx]; $child_mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]; $child_mdl_idx++) { 
            my $mdl_results_HR = \%{$mdl_results_AAHR->[$child_mdl_idx][$seq_idx]}; # for convenience
            # two scenarios in which we can get a b5e error for the cds-mp parent
            # 1. first model with a prediction has p_5seqflush == 1 and p_5overhang != 0
            # 2. first model with a prediction has p_5seqflush == 1 and p_5overhang == 0 
            #    and is not the first model of the first child
            if((! $seen_hit) && # haven't seen a hit yet
               (exists $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_start"}) && # has a prediction (not 'nop')
               ($mdl_results_HR->{"p_5seqflush"} == 1)) {  # prediction extends to 5' boundary of sequence
              if($mdl_results_HR->{"p_5overhang"} != 0) {
                # 1. first model with a prediction has p_5seqflush == 1 and p_5overhang != 0
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "b5e", $seq_name, $mdl_results_HR->{"p_5overhang"} . " nt from 5' end of mature peptide \#" . ($child_idx+1) . " of $np_children", $FH_HR);
                $b5e_flag = 1;
              }
              elsif(($mdl_results_HR->{"p_5overhang"} == 0) && (($child_idx > 0) || ($child_mdl_idx > $ftr_info_HAR->{"first_mdl"}[$child_ftr_idx]))) { 
                # 2. first model with a prediction has p_5seqflush == 1 and p_5overhang == 0 
                #    and is not the first model of the first child
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "b5e", $seq_name, "$child_idx expected mature peptides not observed on 5' end", $FH_HR);
                $b5e_flag = 1;
              }
              $seen_hit  = 1;
              $child_idx = $np_children; # breaks 'for(my $child_idx' loop;
            }
          }
        }
        # if we added a b5e, step through all (not just primary) children of this feature and add m3e
        if($b5e_flag) { 
          for(my $child_idx = 0; $child_idx < $na_children; $child_idx++) { 
            my $child_ftr_idx = $all_children_idx_A[$child_idx];
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "m5e", $seq_name, "", $FH_HR);
          }
        }
      } # end of 'for($seq_idx' loop
    }
  } # end of 'for($ftr_idx' loop

  return;
}      

#################################################################
# Subroutine:  ftr_results_add_b3e_errors
# Incept:      EPN, Wed Feb 21 13:30:07 2018
#
# Purpose:    Report 'b3e' for features with 'annot_type' of 'multifeature'
#             and type of 'cds-mp.
#             Uses mdl_results to do this. If the final (3'-most) primary
#             child of a multifeature feature should have a 'b3e' error
#             (by looking at p_3overhang and p_3seqflush values in mdl_results)
#             then the parent multifeature should too.
#
#             Checks for and adds or updates the following error 
#             codes for features with "annot_type" eq "multifeature' and
#             type 'cds-mp':
#             
#             "b3e": adds this error, predicted hit of 3'-most model 
#                    not flush with model end but flush with sequence end on 3'
#                    OR predicted hit of 3'-most model *is* flush with model end
#                    but also flush with sequence end on 3' and *is not* the first
#                    model of the first child of this multifeature parent feature.
#
#
# Arguments: 
#  $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
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
################################################################# 
sub ftr_results_add_b3e_errors { 
  my $sub_name = "ftr_results_add_b3e_errors";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_HAR, $mdl_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $err_ftr_instances_AHHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $ftr_idx;   # counter over features
  my $mdl_idx;   # counter over models
  my $seq_idx;   # counter over sequences
  my $seq_name;  # name of one sequence

  # foreach annot_type:multifeature and type:'cds-mp' feature, 
  # determine if the first primary child model with a prediction has a 'b5e' error, if so add a 'b5e' error to its parent cds-mp feature
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
        my $seen_hit = 0;      # set to '1' once we've seen the first model with an annotated hit
        $seq_name = $seq_info_HAR->{"seq_name"}[$seq_idx];
        my $b3e_flag = 0;
        # step through all primary children of this feature
        for(my $child_idx = $np_children-1; $child_idx >= 0; $child_idx--) { 
          my $child_ftr_idx = $primary_children_idx_A[$child_idx];
          for(my $child_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]; $child_mdl_idx >= $ftr_info_HAR->{"first_mdl"}[$child_ftr_idx]; $child_mdl_idx--) { 
            my $mdl_results_HR = \%{$mdl_results_AAHR->[$child_mdl_idx][$seq_idx]}; # for convenience
            # two scenarios in which we can get a b3e error for the cds-mp parent
            # 1. 3'-most model with a prediction has p_3seqflush == 1 and p_3overhang != 0 (and no trc or ext error)
            # 2. 3'-most model with a prediction has p_3seqflush == 1 and p_3overhang == 0 (and no trc or ext error)
            #    and is not the final model of the final child
            if((! $seen_hit) && # haven't seen a hit yet
               (exists $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_start"}) && # has a prediction (not 'nop')
               (! exists $mdl_results_HR->{"trc_err_flag"}) && # no trc error for this model
               (! exists $mdl_results_HR->{"ext_err_flag"}) && # no ext error for this model
               ($mdl_results_HR->{"p_3seqflush"} == 1)) {  # prediction extends to 3' boundary of sequence
              if($mdl_results_HR->{"p_3overhang"} != 0) {
                # 1. first model with a prediction has p_3seqflush == 1 and p_3overhang != 0
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "b3e", $seq_name, $mdl_results_HR->{"p_3overhang"} . " nt from 3' end of mature peptide \#" . ($child_idx+1) . " of $np_children", $FH_HR);
                $b3e_flag = 1;
              }
              elsif(($mdl_results_HR->{"p_3overhang"} == 0) && (($child_idx < ($np_children-1)) || ($child_mdl_idx < $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]))) { 
                # 2. first model with a prediction has p_3seqflush == 1 and p_3overhang == 0 
                #    and is not the final model of the final child
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "b3e", $seq_name, ($np_children-1)-$child_idx . " expected mature peptides not observed on 3' end", $FH_HR);
                $b3e_flag = 1;
              }
              $seen_hit = 1;
              $child_idx = -1; # breaks 'for(my $child_idx' loop;
            }
          }
        }
        # if we added a b3e, step through all (not just primary) children of this feature and add m3e
        if($b3e_flag) { 
          for(my $child_idx = 0; $child_idx < $na_children; $child_idx++) { 
            my $child_ftr_idx = $all_children_idx_A[$child_idx];
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "m3e", $seq_name, "", $FH_HR);
          }
        }
      } # end of 'for($seq_idx' loop
    }
  } # end of 'for($ftr_idx' loop

  return;
}      

#################################################################
# Subroutine:  ftr_results_add_blastx_errors
# Incept:      EPN, Tue Oct 23 15:54:50 2018
#
# Purpose:    Report 'x**' and 'mxi' errors for features of type 'cds-notmp' and 'cds-mp'
#             Uses ftr_results to do this. Possible x** errors are:
#
#             "xnh": adds this error if blastx validation of a CDS prediction fails due to
#                    no blastx hits
#             "xos": adds this error if blastx validation of a CDS prediction fails due to
#                    strand mismatch between CM and blastx prediction
#             "x5l": adds this error if blastx validation of a CDS prediction fails due to
#                    BLASTX alignment being too long on 5' end (extending past CM alignment by > 0 nt)
#             "x5s": adds this error if blastx validation of a CDS prediction fails due to
#                    BLASTX alignment being too short on 5' end (more than $xalntol shorter than CM)
#             "x3l": adds this error if blastx validation of a CDS prediction fails due to
#                    BLASTX alignment being too long on 3' end (extending past CM alignment by > 0 nt)
#             "x3s": adds this error if blastx validation of a CDS prediction fails due to
#                    BLASTX alignment being too short on 3' end (more than $xalntol shorter than CM)
#             "xin": adds this error if blastx validation of a CDS prediction fails due to
#                    too long of an insert
#             "xde": adds this error if blastx validation of a CDS prediction fails due to
#                    too long of a delete
#             "xtr": adds this error if blastx validation of a CDS prediction fails due to
#                    an in-frame stop codon in the blastx alignment
#             "xnn": adds this error if blastx has a prediction of sufficient score 
#                    for a feature for which there is no CM/nucleotide based prediction
#
# Arguments: 
#  $out_FH:                 file handle to output blast table to 
#  $query_width:            max length of any query name
#  $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:           REF to hash of arrays with information on the sequences, PRE-FILLED
#  $ftr_results_AAHR:       REF to feature results AAH, PRE-FILLED
#  $mdl_results_AAHR:       REF to model results AAH, PRE-FILLED
#  $err_ftr_instances_AHHR: REF to error instances AHH, ADDED TO HERE
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:                REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies: if we have a multi-segment cds-notmp feature and are unable to find a p
################################################################# 
sub ftr_results_add_blastx_errors { 
  my $sub_name = "ftr_results_add_blastx_errors";
  my $nargs_exp = 10;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_FH, $query_width, $ftr_info_HAR, $seq_info_HAR, $ftr_results_AAHR, $mdl_results_AAHR, $err_ftr_instances_AHHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nftr = validateFeatureInfoHashIsComplete ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $ftr_idx;   # counter over features
  my $seq_idx;   # counter over sequences
  my $seq_name;  # name of one sequence

  my $aln_tol     = opt_Get("--xalntol",    $opt_HHR); # maximum allowed difference between start/end point prediction between CM and blastx
  my $indel_tol   = opt_Get("--xindeltol",  $opt_HHR); # maximum allowed insertion and deletion length in blastx output
  my $min_x_score = opt_Get("--xlonescore", $opt_HHR); # minimum score for a lone hit (no corresponding CM prediction) to be considered

  my $seq_name_width    = maxLengthScalarValueInArray($seq_info_HAR->{"seq_name"});
  my $ftr_product_width = maxLengthScalarValueInArray($ftr_info_HAR->{"out_product"});
  if($seq_name_width    < length("#sequence")) { $seq_name_width    = length("#sequence"); }
  if($ftr_product_width < length("product"))   { $ftr_product_width = length("product"); }
  if($query_width       < length("bquery"))    { $query_width   = length("bquery"); }

  my @out_per_seq_AA = (); # [0..$nseq-1]: per-cds feature output lines for each sequence, we output at end

  # foreach type:'cds-mp' or 'cds-notmp' feature, 
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-mp") ||
       ($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-notmp")) { 

      my @all_children_idx_A = ();
      my $na_children = 0;
      if($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-mp") { 
        # get the all children array
        @all_children_idx_A = (); # feature indices of the primary children of this feature
        getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "all", \@all_children_idx_A, $FH_HR);
        $na_children = scalar(@all_children_idx_A);
      }

      for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
        $seq_name = $seq_info_HAR->{"seq_name"}[$seq_idx];
        my $ftr_results_HR = \%{$ftr_results_AAHR->[$ftr_idx][$seq_idx]}; # for convenience

        # determine if we have any CM predictions for any models related to this feature
        my $xnn_err_possible = 1; 
        if(check_for_defined_pstart_in_mdl_results($seq_idx, $ftr_idx, $ftr_info_HAR, $mdl_results_AAHR, $FH_HR)) { 
          $xnn_err_possible = 0; # we have at least one prediction for this feature, we can't have a xnn error
        }

        my %err_str_H = ();   # added to as we find errors below, possible keys are:
                              # "xnn", "xnh", "xws", "xin", "xde", "xst", "xnn", "x5u", "x3u"

        # initialize 
        my $p_start        = undef; # predicted start  from CM 
        my $p_stop         = undef; # predicted stop   from CM 
        my $p_strand       = undef; # predicted strand from CM 
        my $x_start        = undef; # predicted start  from blastx
        my $x_stop         = undef; # predicted stop   from blastx
        my $x_start2print  = undef; # predicted start from blastx, to output
        my $x_stop2print   = undef; # predicted stop  from blastx, to output
        my $x_strand       = undef; # predicted strand from blastx
        my $x_maxins       = undef; # maximum insert from blastx
        my $x_maxdel       = undef; # maximum delete from blastx
        my $x_trcstop      = undef; # premature stop from blastx
        my $x_score        = undef; # raw score from blastx
        my $x_query        = undef; # query name from blastx hit
        my $x_qlen         = undef; # length of query sequence, if $x_feature_flag == 1
        my $x_feature_flag = 0; # set to '1' if $x_query is a fetched feature sequence, not a full length input sequence

        my $start_diff = undef; # difference in start values between CM and blastx
        my $stop_diff  = undef; # difference in start values between CM and blastx
        my $p_has_stop = undef; # '1' if predicted CM stop ends with a stop codon, else '0'
        
        # first, determine predicted start/stop/strand from CM and blastx for this feature
        # and get the p_start and p_stop values from "out_start", "out_stop"
        # type == cds-mp
        if($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-mp") { 
          if((defined $ftr_results_HR->{"out_start"})  && 
             ($ftr_results_HR->{"out_start"} ne "?") && 
             (defined $ftr_results_HR->{"out_stop"})  && 
             ($ftr_results_HR->{"out_stop"} ne "?")) { 
            $p_start    = $ftr_results_HR->{"out_start"};
            $p_stop     = $ftr_results_HR->{"out_stop"};
            $p_strand   = ($p_start <= $p_stop) ? "+" : "-";
            $p_has_stop = ((defined $ftr_results_HR->{"out_stop_codon"}) && 
                           (validateStopCodon($ftr_results_HR->{"out_stop_codon"}))) ? 1 : 0;
          }
        }
        # type == cds-notmp
        elsif($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-notmp") { 
          my $first_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx];
          my $final_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$ftr_idx];
          if($first_mdl_idx == $final_mdl_idx) { # single exon 
            if((defined $mdl_results_AAHR->[$first_mdl_idx][$seq_idx]{"out_start"}) && 
               ($mdl_results_AAHR->[$first_mdl_idx][$seq_idx]{"out_start"} ne "?")  && 
               (defined $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"out_stop"}) && 
               ($mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"out_stop"} ne "?")) { 
              $p_start = $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"out_start"};
              $p_stop  = $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"out_stop"};
              $p_strand = ($p_start <= $p_stop) ? "+" : "-";
              $p_has_stop = ((defined $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"out_stop_codon"}) &&
                             (validateStopCodon($mdl_results_AAHR->[$final_mdl_idx][$seq_idx]{"out_stop_codon"}))) ? 1 : 0;
            }
          }
          else { 
            for(my $cur_mdl_idx = $first_mdl_idx; $cur_mdl_idx <= $final_mdl_idx; $cur_mdl_idx++) { 
              # looking for first valid p_start, and final valid p_stop
              # these can (and probably will) be different segments, but if so we make sure they are on the same strand
              if((! defined $p_start) && 
                 (defined $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_start"}) && 
                 ($mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_start"} ne "?") && 
                 (defined $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop"}) && 
                 ($mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop"} ne "?")) { 
                $p_start = $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_start"};
                $p_stop  = $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop"};
                $p_strand = ($p_start <= $p_stop) ? "+" : "-";
                $p_has_stop = ((defined $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop_codon"}) &&
                               (validateStopCodon($mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop_codon"}))) ? 1 : 0;
              }
              elsif((defined $p_start) && 
                    (defined $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_start"}) && 
                    ($mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_start"} ne "?") && 
                    (defined $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop"}) && 
                    ($mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop"} ne "?")) { 
                if((($p_strand eq "+") && ($mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_start"} <= $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop"})) ||
                   (($p_strand eq "-") && ($mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_start"}  > $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop"}))) { 
                  # same strand as current p_start, p_stop, p_strand
                  $p_stop  = $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop"};
                  $p_has_stop = ((defined $mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop_codon"}) &&
                                 (validateStopCodon($mdl_results_AAHR->[$cur_mdl_idx][$seq_idx]{"out_stop_codon"}))) ? 1 : 0;
                }
              }
            } # end of 'for(my $cur_mdl_idx = $first_mdl_idx; $cur_mdl_idx <= $final_mdl_idx; $cur_mdl_idx++)'
            if((! defined $p_start) && (! $xnn_err_possible)) { 
              DNAORG_FAIL("ERROR, in $sub_name, feature with multiple segments, unable to find a valid start and stop to compare to blastx but xnn_err_possible was set to true", 1, $FH_HR);
            }
          }
#         if((defined $p_start) &&
#             (defined $p_stop)  && 
#             (defined $p_has_stop)) { 
#            printf("HEYA TEMP p_start: $p_start p_stop: $p_stop p_has_stop: $p_has_stop\n");
#          }
        }

        if((defined $ftr_results_HR->{"x_start"}) && 
           (defined $ftr_results_HR->{"x_stop"})) { 
          $x_start   = $ftr_results_HR->{"x_start"};
          $x_stop    = $ftr_results_HR->{"x_stop"};
          $x_strand  = $ftr_results_HR->{"x_strand"};
          $x_query   = $ftr_results_HR->{"x_query"};
          if(defined $ftr_results_HR->{"x_maxins"}) { 
            $x_maxins  = $ftr_results_HR->{"x_maxins"};
          }
          if(defined $ftr_results_HR->{"x_maxdel"}) { 
            $x_maxdel  = $ftr_results_HR->{"x_maxdel"};
          }
          if(defined $ftr_results_HR->{"x_trcstop"}) { 
            $x_trcstop = $ftr_results_HR->{"x_trcstop"};
          }
          if(defined $ftr_results_HR->{"x_score"}) { 
            $x_score = $ftr_results_HR->{"x_score"};
          }
          # determine if the query is a full length sequence, or a fetched sequence feature:
          (undef, undef, $x_qlen) = helper_blastx_breakdown_query($x_query, $seq_name, undef, $FH_HR); 
          # helper_blastx_breakdown_query() will exit if $x_query is unparseable
          # first two undefs: seqname after coords_str is removed, and coords_str
          # $x_qlen will be undefined if $x_query is a full sequence name name
          $x_feature_flag = (defined $x_qlen) ? 1 : 0; 
          #printf("HEYA seq_name: $seq_name ftr: $ftr_idx x_query: $x_query x_feature_flag: $x_feature_flag x_start: $x_start x_stop: $x_stop x_score: $x_score\n");
        }

        if(($xnn_err_possible) && (! defined $p_start) && (defined $x_start))  { # no CM prediction but there is a blastx prediction
          if((defined $x_score) && ($x_score >= $min_x_score)) { 
            $err_str_H{"xnn"} = "blastx hit from $x_start to $x_stop with score $x_score, but no CM hit";
          }
        }

        #if(defined $p_start) { 
        #  printf("HEYAA seq $seq_idx ftr_idx $ftr_idx " . $ftr_info_HAR->{"type"}[$ftr_idx] . " p_start: $p_start p_stop: $p_stop p_strand: $p_strand\n");
        #}
        #else { 
        #  printf("HEYAA seq $seq_idx ftr_idx $ftr_idx no p_start\n");
        #}

        # if we have a prediction from the CM, so we should check for xip errors
        if(defined $p_start) { 
          # check for xnh: lack of prediction failure
          if(! defined $x_start) { 
            $err_str_H{"xnh"} = "no blastx hit";
          }
          else { # $x_start is defined, we can compare CM and blastx predictions
            # check for xos: strand mismatch failure, differently depending on $x_feature_flag
            if(((  $x_feature_flag) && ($x_strand eq "-")) || 
               ((! $x_feature_flag) && ($p_strand ne $x_strand))) { 
              $err_str_H{"xos"} = "strand mismatch between nucleotide-based and blastx-based predictions";
            }
            else { 
              # determine $start_diff and $stop_diff, differently depending on if hit
              # was to the full sequence or a fetched features (true if $x_feature_flag == 1)
              if($x_feature_flag) { 
                $start_diff = $x_start - 1; 
                $stop_diff  = $x_qlen - $x_stop;
                $x_start2print = sprintf("$p_start %s $start_diff", ($p_strand eq "+") ? "+" : "-");
                $x_stop2print  = sprintf("$p_stop %s $stop_diff",  ($p_strand eq "+") ? "-" : "+");
              }
              else { 
                $start_diff = abs($p_start - $x_start);
                $stop_diff  = abs($p_stop  - $x_stop);
                $x_start2print = $x_start;
                $x_stop2print  = $x_stop;
              }
              # check for 'x5l': only for non-feature seqs blastx alignment extends outside of nucleotide/CM alignment on 5' end
              if((! $x_feature_flag) && 
                 ((($p_strand eq "+") && ($x_start < $p_start)) || 
                  (($p_strand eq "-") && ($x_start > $p_start)))) { 
                $err_str_H{"x5l"} = "blastx alignment extends outside CM alignment on 5' end (strand:$p_strand CM:$p_start blastx:$x_start2print)";
              }

              # check for 'x5s': blastx 5' end too short, not within $aln_tol nucleotides
              if(! exists $err_str_H{"x5l"}) { # only add x5s if x5l does not exist
                if($start_diff > $aln_tol) { 
                  $err_str_H{"x5s"} = "start positions differ by $start_diff > $aln_tol (strand:$p_strand CM:$p_start blastx:$x_start2print)";
                }                
              }

              # check for 'x3l': blastx alignment extends outside of nucleotide/CM alignment on 3' end
              if((! $x_feature_flag) && 
                 ((($p_strand eq "+") && ($x_stop  > $p_stop)) || 
                  (($p_strand eq "-") && ($x_stop  < $p_stop)))) { 
                $err_str_H{"x3l"} = "blastx alignment extends outside CM alignment on 3' end (strand:$p_strand CM:$p_stop blastx:$x_stop2print)";
              }

              # check for 'x3s': blastx 3' end too short, not within $aln_tol nucleotides
              # for the stop coordinates, we do this differently if the CM prediction 
              # includes the stop codon or not, if it does, we allow 3 more positions different
              my $cur_aln_tol = undef;
              my $cur_stop_str = undef;
              if((defined $p_has_stop) && ($p_has_stop == 1)) { 
                $cur_aln_tol  = $aln_tol + 3;
                $cur_stop_str = "valid stop codon";
              }
              else { 
                $cur_aln_tol  = $aln_tol;
                $cur_stop_str = "no valid stop codon";
              }
              if(! exists $err_str_H{"x3l"}) { # only add x3s if x3l does not exist
                if($stop_diff > $cur_aln_tol) { 
                  $err_str_H{"x3s"} = "stop positions differ by $stop_diff > $cur_aln_tol (strand:$p_strand CM:$p_stop blastx:$x_stop2print, $cur_stop_str in CM prediction)";
                }
              }

              # check for 'xin': too big of an insert
              if((defined $x_maxins) && ($x_maxins > $indel_tol)) { 
                $err_str_H{"xin"} = "longest blastx predicted insert of length $x_maxins > $indel_tol";
              }

              # check for 'xde': too big of a deletion
              if((defined $x_maxdel) && ($x_maxdel > $indel_tol)) { 
                $err_str_H{"xde"} = "longest blastx predicted delete of length $x_maxdel > $indel_tol";
              }

              # check for 'xtr': blast predicted truncation
              if(defined $x_trcstop) { 
                $err_str_H{"xtr"} = "blastx alignment includes stop codon ($x_trcstop)";
              }
            }
          }
        } # end of 'if(defined $p_start)'
        my $err_flag = 0;
        my $output_err_str = "";
        foreach my $err_code (sort keys %err_str_H) { 
          error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, $err_code, $seq_name, sprintf("%s: %s", $ftr_info_HAR->{"out_product"}[$ftr_idx], $err_str_H{$err_code}), $FH_HR);
          $err_flag = 1;
          if($output_err_str ne "") { $output_err_str .= ","; }
          $output_err_str .= $err_code;
        }
        if($output_err_str eq "") { $output_err_str = "-"; }
        # if we added an error, step through all (not just primary) children of this feature and add mxi
        if((defined $p_start) && ($err_flag) && ($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-mp")) { 
          for(my $child_idx = 0; $child_idx < $na_children; $child_idx++) { 
            my $child_ftr_idx = $all_children_idx_A[$child_idx];
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "mxi", $seq_name, 
                                sprintf("MP: %s, CDS %s", 
                                        $ftr_info_HAR->{"out_product"}[$child_ftr_idx],
                                        $ftr_info_HAR->{"out_product"}[$ftr_idx]), 
                                $FH_HR);
          }
        }

        # seqname out_product CM blast CM-start CM-stop blastx-start blastx-stop blastx-score startdiff stopdiff blastx-maxins blastx-maxdel blastx-trcstop xip xnn
        if($ftr_idx == 0) { 
          @{$out_per_seq_AA[$seq_idx]} = ();
        }
        # determine if we should output the best blast hit
        # we only do this if there's also a CM hit OR
        # we're above score of at least $min_x_score
        my $x_hit_printable = 0;
        if((defined $x_start) && 
           ((defined $p_start) || ($x_score >= $min_x_score))) { 
          $x_hit_printable = 1;
        }
        push(@{$out_per_seq_AA[$seq_idx]}, sprintf("%-*s  %-*s  %6s  %6s  %-*s  %8s  %7s  %7s  %7s  %7s  %7s  %7s  %7s  %7s  %7s  %7s  %-s\n", 
                                                   $seq_name_width,    $seq_name, 
                                                   $ftr_product_width, $ftr_info_HAR->{"out_product"}[$ftr_idx],
                                                   (defined $p_start)                       ? "yes"       : "no",   # CM prediction?
                                                   (defined $x_start  && $x_hit_printable)  ? "yes"       : "no",   # blastx prediction? 
                                                   $query_width,                                   
                                                   (defined $x_query)                       ? $x_query    : "-",    # query name 
                                                   (defined $x_query && $x_feature_flag)    ? "feature"   : "full", # hit to feature sequence or full sequence?
                                                   (defined $p_start)                       ? $p_start    : "-",    # CM-start
                                                   (defined $p_stop)                        ? $p_stop     : "-",    # CM-stop
                                                   (defined $x_start  && $x_hit_printable)  ? $x_start    : "-",    # blastx-start
                                                   (defined $x_stop   && $x_hit_printable)  ? $x_stop     : "-",    # blastx-stop
                                                   (defined $x_score)                       ? $x_score    : "-",    # blastx-score
                                                   (defined $start_diff)                    ? $start_diff : "-",    # start-diff
                                                   (defined $stop_diff)                     ? $stop_diff  : "-",    # stop-diff
                                                   (defined $x_maxins  && $x_hit_printable) ? $x_maxins   : "-",    # blastx-maxins
                                                   (defined $x_maxdel  && $x_hit_printable) ? $x_maxdel   : "-",    # blastx-maxdel
                                                   (defined $x_trcstop && $x_hit_printable) ? $x_trcstop  : "-",    # blastx-maxdel
                                                   $output_err_str));
            


      } # end of 'for($seq_idx' loop
    }
  } # end of 'for($ftr_idx' loop

  printf $out_FH ("#sequence: sequence name\n");
  printf $out_FH ("#product:  CDS product name\n");
  printf $out_FH ("#cm?:      is there a CM (nucleotide-based) prediction/hit? above threshold\n");
  printf $out_FH ("#blast?:   is there a blastx (protein-based) prediction/hit? above threshold\n");
  printf $out_FH ("#bquery:   name of blastx query name\n");
  printf $out_FH ("#feature?: 'feature' if blastx query was a fetched feature, from CM prediction\n");
  printf $out_FH ("#          'full'    if blastx query was a full input sequence\n");
  printf $out_FH ("#cmstart:  start position of CM (nucleotide-based) prediction\n");
  printf $out_FH ("#cmstop:   stop  position of CM (nucleotide-based) prediction\n");
  printf $out_FH ("#bxstart:  start position of blastx top HSP\n");
  printf $out_FH ("#bxstop:   stop  position of blastx top HSP\n");
  printf $out_FH ("#bxscore:  raw score of top blastx HSP (if one exists, even if it is below threshold)\n");
  printf $out_FH ("#startdf:  difference between cmstart and bxstart\n");
  printf $out_FH ("#stopdf:   difference between cmstop and bxstop\n");
  printf $out_FH ("#bxmaxin:  maximum insert length in top blastx HSP\n");
  printf $out_FH ("#bxmaxde:  maximum delete length in top blastx HSP\n");
  printf $out_FH ("#bxtrc:    position of stop codon in top blastx HSP, if there is one\n");
  printf $out_FH ("#errors:   list of errors for this sequence, - if none\n");
  printf $out_FH ("%-*s  %-*s  %6s  %6s  %-*s  %8s  %7s  %7s  %7s  %7s  %7s  %7s  %7s  %7s  %7s  %7s  %-s\n",
                  $seq_name_width,    "#sequence", 
                  $ftr_product_width, "product",
                  "cm?", "blast?", 
                  $query_width, "bquery", 
                  "feature?", "cmstart", "cmstop", "bxstart", "bxstop", "bxscore", "startdf", "stopdf", "bxmaxin", "bxmaxde", "bxtrc", "errors");

         

  # now go back and output per seq
  for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    foreach my $out_line (@{$out_per_seq_AA[$seq_idx]}) { 
      print $out_FH $out_line;
    }
  }

  return;
}    

#################################################################
# Subroutine:  mdl_results_add_str_nop_ost_lsc_dup_b3e_b3u_errors
# Incept:      EPN, Thu Mar 31 13:43:58 2016
#
# Purpose:    Report 'str', 'nop', 'ost', 'lsc', 'dup', 'b3e', and 'b3u' errors
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
#  $ftr_info_HAR:           REF to hash of arrays with information on the features, PRE-FILLED
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
sub mdl_results_add_str_nop_ost_lsc_dup_b3e_b3u_errors { 
  my $sub_name = "mdl_results_add_str_nop_ost_lsc_dup_b3e_b3u_errors()";
  my $nargs_exp = 9;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile, $ftr_info_HAR, $mdl_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $err_ftr_instances_AHHR, $err_info_HAR, $opt_HHR, $FH_HR) = @_;
  
  # total counts of things
  my $nmdl = validateModelInfoHashIsComplete   ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nseq = validateSequenceInfoHashIsComplete($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $mdl_idx; # counter over models
  my $seq_idx; # counter over sequences

  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    #printf("\tmdl_idx: $mdl_idx\n");
    my $ftr_idx    = $mdl_info_HAR->{"map_ftr"}[$mdl_idx];
    my $is_first   = $mdl_info_HAR->{"is_first"}[$mdl_idx];
    my $is_final   = $mdl_info_HAR->{"is_final"}[$mdl_idx];
    for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
      my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
      my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
      my $seq_len   = $seq_info_HAR->{"len"}[$seq_idx];
      my $mdl_results_HR = \%{$mdl_results_AAHR->[$mdl_idx][$seq_idx]}; # for convenience

      if(exists $mdl_results_HR->{"p_start"}) { 
        ##########################
        # update str err message #
        ##########################
        if($is_first && (exists $err_ftr_instances_AHHR->[$ftr_idx]{"str"}{$seq_name})) { 
          my $out_start = undef;
          ($out_start, undef) = create_output_start_and_stop($mdl_results_HR->{"p_start"}, $mdl_results_HR->{"p_stop"}, 
                                                             $seq_len, $FH_HR);
          my $updated_str_errmsg = sprintf("%s starting at position %d on strand %s", 
                                           fetchStartCodon($sqfile, $seq_name, $mdl_results_HR->{"p_start"}, $mdl_results_HR->{"p_strand"}, 1, $FH_HR), 
                                           # '1' above: it's okay if codon we fetch is less than 3 nt due to being off end of the sequence
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
        if((exists $mdl_results_HR->{"p_score"}) && ($mdl_results_HR->{"p_score"} < 0.)) { 
          #########################
          # lsc error (low score) #
          #########################
          if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"lsc"}{$seq_name}) { # lsc error already exists, update it
            error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "lsc", $seq_name, $err_ftr_instances_AHHR->[$ftr_idx]{"lsc"}{$seq_name} . ", " . $mdl_results_HR->{"p_score"} . " bits", $FH_HR);
          }
          else { # first lsc error
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "lsc", $seq_name, sprintf("%s: %s bits", $ftr_info_HAR->{"out_product"}[$ftr_idx], $mdl_results_HR->{"p_score"}), $FH_HR);
          }
        }
        if($mdl_results_HR->{"p_nhits"} > 1) { 
          ##############################
          # dup error (duplicate hits) #
          ##############################
          if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"dup"}{$seq_name}) { # dup error already exists, update it
            error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "dup", $seq_name, $err_ftr_instances_AHHR->[$ftr_idx]{"dup"}{$seq_name} . ", " . $mdl_results_HR->{"p_nhits"} . " hits", $FH_HR);
          }
          else { # first dup error
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "dup", $seq_name, sprintf("%s: %d regions", $ftr_info_HAR->{"out_product"}[$ftr_idx], $mdl_results_HR->{"p_nhits"}), $FH_HR);
          }
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
      my $seq_len   = $seq_info_HAR->{"len"}[$seq_idx];
      my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
      my $mdl_results_HR = \%{$mdl_results_AAHR->[$mdl_idx][$seq_idx]}; # for convenience

      if(exists $mdl_results_HR->{"p_start"}) { 
        ($mdl_results_HR->{"out_start"}, $mdl_results_HR->{"out_stop"}) = 
            create_output_start_and_stop($mdl_results_HR->{"p_start"},
                                         exists($mdl_results_HR->{"c_stop"}) ? $mdl_results_HR->{"c_stop"} : $mdl_results_HR->{"p_stop"}, 
                                         $seq_len, $FH_HR);
        $mdl_results_HR->{"out_stop_codon"} = fetchStopCodon($sqfile, $seq_name,
                                                             (defined $mdl_results_HR->{"c_stop"}) ? $mdl_results_HR->{"c_stop"} : $mdl_results_HR->{"p_stop"}, 
                                                             $mdl_results_HR->{"p_strand"}, 1, $FH_HR);
                                              # '1' above: it's okay if codon we fetch is less than 3 nt due to being off end of the sequence

      }
    }
  }

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
#             "nm3": for features with annot_type eq ("model" or "multifeature") & type eq "cds-notmp" OR "cds-mp" OR "mp"
#             "stp": for features with annot_type eq "multifeature" & type eq "cds-mp"
#             "inp": for features with annot_type eq "multifeature" & type eq "cds-mp"
#             "int": for features with annot_type eq "multifeature" & type eq "cds-mp"
#             "aji": for features with annot_type eq "multifeature" & type eq "cds-mp"
#             "mtr": for features that are children of a multifeature/cds-mp feature
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
  my $append_len = 0; # length of appended region
  my $mtr_errmsg = undef; # error message for an mtr error

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
        $seq_len    = $seq_info_HAR->{"len"}[$seq_idx];
        $accn_name  = $seq_info_HAR->{"accn_name"}[$seq_idx];

        # set the str_err_flag, if nec
        if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"str"}{$seq_name}) { 
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"str_err_flag"} = 1;
        }

        # initialize our error-related variables
        my $aji_errmsg          = "";    # filled if we find two mature peptides that should be adjacent but are not
        my $set_start           = 0;     # set to '1' once we've seen the first model with an annotated hit
        my $inp_errmsg          = "";    # a list of model names ("out_tiny") for which we do not have annotations (nop errors), if any
        my $int_errmsg          = "";    # a list of model names ("out_tiny") which are not translated due to trc (mtr errors), if any
        my $cds_out_start       = undef; # start position to output for this CDS 
        my $cds_out_stop        = undef; # stop  position to output for this CDS 
        my $cds_fetch_start     = undef; # start position to *fetch* for this CDS' start codon
        my $cds_fetch_stop      = undef; # stop  position to *fetch* for this CDS' stop codon
        my $start_strand        = undef; # strand start codon is on
        my $stop_strand         = undef; # strand stop codon is on
        my $cds_len             = 0;     # length to output for this CDS
        my $child_had_trc       = 0;     # if we find a child with a trc error, set this to 1
        my $mn3_flag            = 0;     # set to '1' if we set a bad nm3 (no b5e or b3e) error for mother CDS
        my $maj_flag            = 0;     # set to '1' if we set a aji error for mother CDS
        my $mit_flag            = 0;     # set to '1' if we set a int error for mother CDS
        my $mip_flag            = 0;     # set to '1' if we set a inp error for mother CDS

        my $seen_prv_b3e        = 0;     # set to '1' if we see a b3e error for a MP for the mother CDS

        # first make sure that we have at least 1 prediction in a primary child, if not we won't have prediction for this CDS so skip
        my $cds_any_pred_flag = 0;
        for(my $child_idx = 0; $child_idx < $np_children; $child_idx++) { 
          my $child_ftr_idx = $primary_children_idx_A[$child_idx];
          for(my $tmp_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$child_ftr_idx]; $tmp_mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]; $tmp_mdl_idx++) { 
            if(exists $mdl_results_AAHR->[$tmp_mdl_idx][$seq_idx]{"p_start"}) { 
              $cds_any_pred_flag = 1;
            }
          }
        }

        if($cds_any_pred_flag) { 
          # step through all primary children of this feature
          for(my $child_idx = 0; $child_idx < $np_children; $child_idx++) { 
            my $child_ftr_idx = $primary_children_idx_A[$child_idx];
            if((exists $err_ftr_instances_AHHR->[$child_ftr_idx]{"b3e"}{$seq_name}) || # we have a b3e for this mat_peptide
               ($mdl_results_AAHR->[$ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]][$seq_idx]{"p_3seqflush"})) { # rare case: we don't have a b3e but final nt of prediction is final nt of sequence, so for our purposes here, we do have a b3e
              $seen_prv_b3e = 1;
            }

            for(my $child_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$child_ftr_idx]; $child_mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]; $child_mdl_idx++) { 
              if(! $child_had_trc) { # if $child_had_trc is true, then we dealt with this feature completely below
                # check to make sure we have a hit annotated for this model
                # if not, we trigger an inp UNLESS we don't expect a hit due to a relevant b5e or relevant b3e
                # we have a relevant b5e if: the CDS has a b5e AND $set_start has NOT yet been set
                # we have a relevant b3e if: either this MP or a previous one has a b3e ($seen_prv_b3e == 1)
                # we do not have a relevant b5e if: $set_start has been set
                # we do not have a relevant b3e if: this MP does not have a b3e and no previous MP had a b3e
                if(! exists $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_start"}) { 
                  # we might need to trigger an inp, check if we have a relevant b5e or b3e
                  my $have_relevant_b5e = ((! $set_start) && (exists $err_ftr_instances_AHHR->[$ftr_idx]{"b5e"}{$seq_name})) ? 1 : 0;
                  my $have_relevant_b3e = ($seen_prv_b3e) ? 1 : 0;
                  if((! $have_relevant_b5e) && (! $have_relevant_b3e)) { 
                    if($inp_errmsg ne "") { 
                      $inp_errmsg .= ", ";
                    }
                    $inp_errmsg .= $mdl_info_HAR->{"out_tiny"}[$child_mdl_idx];
                  }
                }
                else { # we do have a hit for $child_mdl_idx in $seq_idx
                  $cds_len += $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"len"};
                  
                  if(! $set_start) { # first model, set cds_out_start
                    $cds_out_start   = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_start"};
                    $cds_fetch_start = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_start"};
                    $start_strand    = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_strand"};
                    $set_start = 1;
                  }
                  # (as of v0.28: always update the stop, remove this to revert to pre-v0.28 behavior
                  #  in regards to stop output in tbl) 
                  $cds_out_stop = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_stop"};
                  $cds_fetch_stop    = (defined $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"c_stop"}) ? 
                      $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"c_stop"} :
                      $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_stop"};
                  $stop_strand    = $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"p_strand"}; 
                  
                  # check if we have a trc in this child model, and deal with it if we do
                  if(exists $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"trc_err_flag"}) { 
                    $child_had_trc = 1;
                    error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "ctr", $seq_name, 
                                        sprintf("mat_peptide %s includes trc error", $ftr_info_HAR->{"out_product"}[$mdl_info_HAR->{"map_ftr"}[$child_mdl_idx]]), $FH_HR);
                    ####################################################
                    # We have a trc in this child model $child_mdl_idx 
                    # 
                    # Add the 'ctr' error for this CDS (Child has TRc)
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
                    # - set mtr errors for all remaining children 
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
                                                        create_output_start_and_stop($cds_fetch_start, $cds_pred_stop, $seq_len, $FH_HR),
                                                        create_output_start_and_stop($cds_fetch_start, $cds_out_stop,  $seq_len, $FH_HR),
                                                        abs($cds_fetch_stop - $cds_pred_stop), $mdl_info_HAR->{"out_tiny"}[$child_mdl_idx]);
                        }
                        else { 
                          $updated_trc_errmsg = sprintf("homology search predicted %d..? revised to %d..%d (due to early stop in %s)", 
                                                        $cds_out_start, 
                                                        create_output_start_and_stop($cds_out_start, $cds_out_stop, $seq_len, $FH_HR),
                                                        $mdl_info_HAR->{"out_tiny"}[$child_mdl_idx]);
                        }
                        
                        if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name}) { 
                          # update it
                          error_instances_update($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_name, $updated_trc_errmsg, $FH_HR);
                          # set the trc_err_flag for this feature
                          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"trc_err_flag"} = 1;
                        }
                        else { 
                          # it doesn't yet exist, so add the trc error IF 
                          # we don't have a b5e for this guy. This
                          # is rare, but we may have no trc error for this
                          # CDS yet, if the predicted child mature peptide
                          # sequences didn't all exist, then we won't have
                          # combined those predictions into a predicted CDS,
                          # and thus we didn't check that predicted CDS for
                          # trc errors in
                          # parse_esl_epn_translate_startstop_outfile().
                          if(! exists $err_ftr_instances_AHHR->[$ftr_idx]{"b5e"}{$seq_name}) { 
                            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "trc", $seq_name, $updated_trc_errmsg, $FH_HR);
                            # set the trc_err_flag for this feature
                            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"trc_err_flag"} = 1;
                          }
                        }
                      }
                      # all remaining children get a 'mtr' error,
                      # and the CDS gets an 'int' error, which we need
                      # to build the error message for
                      $mtr_errmsg = sprintf("early stop in mat_peptide %s ending at position %d", $ftr_info_HAR->{"out_product"}[$child_ftr_idx], $cds_out_stop);
                      if($int_errmsg ne "") { 
                        DNAORG_FAIL("ERROR in $sub_name, setting int errmsg for ftr_idx: $ftr_idx due to 'trc', but it is not blank", 1, $FH_HR);
                      }
                      $child_idx++;
                      my $mtr_err_ct = 0;
                      # for all remaining children: throw 'mtr' and append to 'int' err message
                      while($child_idx < $np_children) {
                        $child_ftr_idx = $primary_children_idx_A[$child_idx];
                        # check if we have predictions for any (probably just 1) of the models for this feature
                        my $any_pred_flag = 0;
                        for(my $tmp_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$child_ftr_idx]; $tmp_mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]; $tmp_mdl_idx++) { 
                          if(exists $mdl_results_AAHR->[$tmp_mdl_idx][$seq_idx]{"p_start"}) { 
                            $any_pred_flag = 1;
                          }
                        }
                        # if we do have predictions for any of the models for this feature, add mtr error
                        if($any_pred_flag) { 
                          error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "mtr", $seq_name, $mtr_errmsg, $FH_HR);
                          $mtr_err_ct++;
                          if($int_errmsg ne "") { 
                            $int_errmsg .= ", ";
                          }
                          $int_errmsg .= $ftr_info_HAR->{"out_tiny"}[$child_ftr_idx];
                        }
                        $child_idx++;
                      }
                      if($mtr_err_ct > 0) { 
                        # we set at least one mtr error for mature peptides, set int for this CDS
                        error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "int", $seq_name, $int_errmsg, $FH_HR);
                        $mit_flag = 1; # this causes 'mit' in children MPs
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
                    # we're not the final child feature, and we haven't seen a b3e 
                    # for this CDS yet
                    if(($child_idx < ($np_children-1)) && 
                       ($child_mdl_idx == $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]) && 
                       (! $seen_prv_b3e)) { 
                      my $nxt_child_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$primary_children_idx_A[($child_idx+1)]];
                      # check if they are adjacent 
                      if(! checkForIndexInOverlapOrAdjacencyIndexString($mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"idx_aja_str"}, $nxt_child_mdl_idx, $FH_HR)) { 
                        if($aji_errmsg ne "") { $aji_errmsg .= ", "; }
                        $aji_errmsg .= sprintf("%s (%s..%s) not adjacent to %s (%s..%s)", 
                                               $ftr_info_HAR->{"out_product"}[$mdl_info_HAR->{"map_ftr"}[$child_mdl_idx]], 
                                               (defined $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_start"}) ? $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_start"} : "unknown", 
                                               (defined $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_stop"})  ? $mdl_results_AAHR->[$child_mdl_idx][$seq_idx]{"out_stop"}  : "unknown", 
                                               $ftr_info_HAR->{"out_product"}[$mdl_info_HAR->{"map_ftr"}[$nxt_child_mdl_idx]], 
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
                                                         $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
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
              my $stp_err_stop_codon = fetchStopCodon($sqfile, $seq_name, $cds_fetch_stop, $stop_strand, 1, $FH_HR);
              # '1' above: it's okay if codon we fetch is less than 3 nt due to being off end of the sequence
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
                                             $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
            # only the final model is affected by an ext error
            $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"len"}    += ($len_corr - $append_len);
            $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"cumlen"} += ($len_corr - $append_len);
            $cds_len += ($len_corr - $append_len);
            
            # get first part of cds error message using current $cds_out_stop, if it exists
            my $updated_cds_ext_errmsg = sprintf("homology search predicted %d..%d", 
                                                 create_output_start_and_stop($cds_out_start, $cds_out_stop,
                                                                              $seq_info_HAR->{"len"}[$seq_idx], $FH_HR));
            # get first part of mp error message 
            my $mp_ext_errmsg = sprintf("homology search predicted %d..%d", 
                                        create_output_start_and_stop($mdl_results_AAHR->[$final_first_child_mdl_idx][$seq_idx]{"p_start"},
                                                                     $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"p_stop"},
                                                                     $seq_info_HAR->{"len"}[$seq_idx], $FH_HR));
            # now recompute $cds_out_stop
            (undef, $cds_out_stop) = 
                create_output_start_and_stop($cds_fetch_start, # this is irrelevant due to first undef arg
                                             $cds_fetch_stop,
                                             $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);

            # get second part of CDS error message
            $updated_cds_ext_errmsg .= sprintf(" revised to %d..%d (stop shifted %d nt)", $cds_out_start, $cds_out_stop, $len_corr-$append_len);
            # get second part of MP error message
            $mp_ext_errmsg .= sprintf(" revised to %d..%d (stop shifted %d nt)", 
                                      create_output_start_and_stop($mdl_results_AAHR->[$final_first_child_mdl_idx][$seq_idx]{"p_start"}, 
                                                                   $mdl_results_AAHR->[$final_final_child_mdl_idx][$seq_idx]{"c_stop"}, 
                                                                   $seq_info_HAR->{"len"}[$seq_idx], $FH_HR), 
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
            $maj_flag = 1; # this causes 'maj' in children MPs
          }

          # add the inp (INterrupted due to no Prediction) if necessary
          if($inp_errmsg ne "") { 
            error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "inp", $seq_name, $inp_errmsg, $FH_HR);
            $mip_flag = 1; # this causes 'mip' in children MPs
          }

          # set ftr_results, we can set start if $cds_out_start is defined, 
          if(defined $cds_out_start) { 
            my $start_codon = fetchStartCodon($sqfile, $seq_name, $cds_fetch_start, $start_strand, 1, $FH_HR);
            # '1' above: it's okay if codon we fetch is less than 3 nt due to being off end of the sequence
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
          
          if((defined $cds_out_start) && (defined $cds_out_stop)) { 
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"ftbl_out_start"} = $cds_out_start;
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"ftbl_out_stop"}  = $cds_out_stop;
          }
          else { 
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"ftbl_out_start"}       = "?";
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"ftbl_out_start_codon"} = "?";
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
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_stop_codon"} = fetchStopCodon($sqfile, $seq_name, $cds_fetch_stop, $stop_strand, 1, $FH_HR);
            # '1' above: it's okay if codon we fetch is less than 3 nt due to being off end of the sequence
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"out_len"}       = $cds_len;
            if(($cds_len % 3) != 0) { 
              error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $ftr_idx, "nm3", $seq_name, "$cds_len", $FH_HR);
              # if this is a 'bad' nm3 (no b5e or b3e), then update $mn3_flag
              if((! exists $err_ftr_instances_AHHR->[$ftr_idx]{"b5e"}{$seq_name}) && 
                 (! exists $err_ftr_instances_AHHR->[$ftr_idx]{"b3e"}{$seq_name})) { 
                $mn3_flag = 1;
              }
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
                                                $seq_info_HAR->{"len"}[$seq_idx], 
                                                $cds_tbl_HHAR->{$accn_name}, $opt_HHR, $FH_HR);
            }
            else { # annotation doesn't exist, so we don't have a match
              $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"genbank_ftr_annot_match"} = 0;
            }
          }
          # one final step: if we have a 'trc' error for this CDS, check the 'all children' array, 
          # and throw 'mtr' errors for any mature peptides encoded by this CDS that are not
          # translated. We did this above for the primary peptides, but here we do it for any
          # non-primary peptides.
          if(exists $err_ftr_instances_AHHR->[$ftr_idx]{"trc"}{$seq_name}) { 
            for(my $child_idx = 0; $child_idx < $na_children; $child_idx++) { 
              my $child_ftr_idx = $all_children_idx_A[$child_idx];
              # check if we have predictions for all of the models for this feature
              my $all_pred_flag = 1;
              for(my $tmp_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$child_ftr_idx]; $tmp_mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx]; $tmp_mdl_idx++) { 
                if(! exists $mdl_results_AAHR->[$tmp_mdl_idx][$seq_idx]{"p_start"}) { 
                  $all_pred_flag = 0;
                }
              }
              if($all_pred_flag) { # there is a prediction for all models for this feature
                my $final_child_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$child_ftr_idx];
                my $cur_stop = (defined $mdl_results_AAHR->[$final_child_mdl_idx][$seq_idx]{"c_stop"}) ? 
                    $mdl_results_AAHR->[$final_child_mdl_idx][$seq_idx]{"c_stop"} :
                    $mdl_results_AAHR->[$final_child_mdl_idx][$seq_idx]{"p_stop"};
                if(($stop_strand eq "+") && ($cds_fetch_stop < $cur_stop) && (! exists $err_ftr_instances_AHHR->[$child_ftr_idx]{"mtr"}{$seq_name})) { 
                  error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "mtr", $seq_name, $mtr_errmsg, $FH_HR);
                }
                if(($stop_strand eq "-") && ($cds_fetch_stop > $cur_stop) && (! exists $err_ftr_instances_AHHR->[$child_ftr_idx]{"mtr"}{$seq_name})) { 
                  error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "mtr", $seq_name, $mtr_errmsg, $FH_HR);
                }
              }
            }
          }
          #######################################################################################
          # add mn3, maj, mit, mip errors for all children if mother CDS had corresponding error
          #######################################################################################
          if($mn3_flag || $maj_flag || $mit_flag || $mip_flag) { 
            for(my $child_idx = 0; $child_idx < $na_children; $child_idx++) { 
              my $child_ftr_idx = $all_children_idx_A[$child_idx];
              if($mn3_flag) { 
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "mn3", $seq_name, 
                                    sprintf("MP: %s, CDS %s", 
                                            $ftr_info_HAR->{"out_product"}[$child_ftr_idx],
                                            $ftr_info_HAR->{"out_product"}[$ftr_idx]), 
                                    $FH_HR);
              }
              if($maj_flag) { 
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "maj", $seq_name, 
                                    sprintf("MP: %s, CDS %s", 
                                            $ftr_info_HAR->{"out_product"}[$child_ftr_idx],
                                            $ftr_info_HAR->{"out_product"}[$ftr_idx]), 
                                    $FH_HR);
              }
              if($mit_flag) { 
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "mit", $seq_name, 
                                    sprintf("MP: %s, CDS %s", 
                                            $ftr_info_HAR->{"out_product"}[$child_ftr_idx],
                                            $ftr_info_HAR->{"out_product"}[$ftr_idx]), 
                                    $FH_HR);
              }
              if($mip_flag) { 
                error_instances_add($err_ftr_instances_AHHR, undef, $err_info_HAR, $child_ftr_idx, "mip", $seq_name, 
                                    sprintf("MP: %s, CDS %s", 
                                            $ftr_info_HAR->{"out_product"}[$child_ftr_idx],
                                            $ftr_info_HAR->{"out_product"}[$ftr_idx]), 
                                    $FH_HR);
              }
            }
          }
        } # end of 'if($cds_any_pred_flag)'
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
      my $is_cds     = featureTypeIsCds($ftr_info_HAR->{"type"}[$ftr_idx]); 
      my $is_matpept = featureTypeIsMaturePeptide($ftr_info_HAR->{"type"}[$ftr_idx]); 
      if($is_cds || $is_matpept) { 
        for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
          $seq_name   = $seq_info_HAR->{"seq_name"}[$seq_idx];
          $accn_name  = $seq_info_HAR->{"accn_name"}[$seq_idx];
          $seq_len    = $seq_info_HAR->{"len"}[$seq_idx];
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
                  $seq_info_HAR->{"len"}[$s]);
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
    my $seq_len   = $seq_info_HAR->{"len"}[$seq_idx];
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
      if($qseq_posn <= $seq_len) { # note: we've just incremented qseq_posn by 1 in prv line so now it is in 1..length($seq) coords, not 0..length($seq)-1
        my $qseq_start = $qseq_posn;
        my $qseq_stop  = $qseq_posn + length($qseq) - 1;
        # adjust coordinates so they're within 1..L
        ($qseq_start, $qseq_stop) = 
            create_output_start_and_stop($qseq_start, $qseq_stop, $seq_len, $FH_HR);
        if($seq_info_HAR->{"origin_coords_str"}[$seq_idx] ne "") { 
          $seq_info_HAR->{"origin_coords_str"}[$seq_idx] .= ",";
        }
        $seq_info_HAR->{"origin_coords_str"}[$seq_idx] .= $qseq_start . ":" . $qseq_stop;
        $nfound++;
      }
      $qseq_posn = index($seq, $qseq, $qseq_posn);
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
  my $seq_len = $seq_info_HAR->{"len"}[$seq_idx];

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
    if($oseq_offset > ($seq_len / 2)) { # simpler (shorter distance) to move origin clockwise
      $oseq_offset = $seq_len - $oseq_offset; # note, we don't add 1 here
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

  my $tmp_errmsg = "ftr_idx: $ftr_idx, err_code: $err_code, seq_name: $seq_name, value: $value\n";
#  if($err_code eq "olp" || 
#     $err_code eq "ajb" || 
#     $err_code eq "aja" || 
#     $err_code eq "inp" || 
#     $err_code eq "lsc" || 
#     $err_code eq "dup" || 
#     $err_code eq "aji" || 
#     $err_code eq "maj" || 
#     $err_code eq "mip" || 
#     $err_code eq "ost" || 
#     $err_code eq "ori") { 
#    DNAORG_FAIL("ERROR in $sub_name, $tmp_errmsg", 1, $FH_HR);
#  }
  
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
  my $do_xfeat     = (numNonNumericValueInArray($ftr_info_HAR->{"type"}, "xfeat", $FH_HR) > 0) ? 1 : 0;
  my $do_dfeat     = (numNonNumericValueInArray($ftr_info_HAR->{"type"}, "dfeat", $FH_HR) > 0) ? 1 : 0;

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
    } # end of 'if' entered if feature is a multifeature cds-mp feature
  
    #############################################
    # block that handles 'annot_type' eq "model"
    # features, these are cds-notmp and mp types
    #############################################
    elsif($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { 
      for(my $mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx]; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
        $width += 18;
        my $mdl_exon_idx = $mdl_info_HAR->{"map_exon"}[$mdl_idx];
        my $is_matpept   = featureTypeIsMaturePeptide($ftr_info_HAR->{"type"}[$ftr_idx]);
        my $is_cds       = featureTypeIsCds($ftr_info_HAR->{"type"}[$ftr_idx]);
        my $is_xfeat     = featureTypeIsExtraFeature($ftr_info_HAR->{"type"}[$ftr_idx]);
        my $is_dfeat     = featureTypeIsDuplicateFeature($ftr_info_HAR->{"type"}[$ftr_idx]);
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
          if($do_matpept && $do_cds_notmp && ($do_xfeat || $do_dfeat)) { 
            $exp_tok1 = "{MP,CDS,other} #<i>";
          }
          elsif($do_matpept && $do_cds_notmp) { 
            $exp_tok1 = "{MP,CDS} #<i>";
          }
          elsif($do_matpept && (! $do_cds_notmp) && ($do_xfeat || $do_dfeat)) { 
            $exp_tok1 = "{MP,other} #<i>";
          }
          elsif($do_matpept && (! $do_cds_notmp)) { 
            $exp_tok1 = "MP #<i>";
          }
          elsif($do_cds_notmp && ($do_xfeat || $do_dfeat)) { 
            $exp_tok1 = "{CDS,other} #<i>";
          }
          else { 
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
        if((! $do_xfeat) && (! $do_dfeat)) { 
          if($do_matpept && $do_cds_notmp) { 
            $exp_substr = "coding sequence part (or exon) <j> of mat_peptide (or CDS)";
          }
          elsif($do_matpept && (! $do_cds_notmp)) { 
            $exp_substr = "coding sequence part <j> of mat_peptide";
          }
          else { # $do_matpept is false
            $exp_substr = "exon <j> of CDS";
          }
        }
        else { # $do_xfeat or $do_dfeat is TRUE
          if($do_matpept && $do_cds_notmp) { 
            $exp_substr = "coding sequence part (or exon or segment) <j> of mat_peptide (or CDS or other feature)";
          }
          elsif($do_matpept && (! $do_cds_notmp)) { 
            $exp_substr = "coding sequence part (or segment) <j> of mat_peptide (or other feature)";
          }
          else { # $do_matpept is false
            $exp_substr = "exon (or segment) <j> of CDS (or other feature)";
          }
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

        $exp_substr = "sequence";
        if($is_matpept) { 
          $exp_substr = "mat_peptide coding sequence";
        }
        if($is_cds) { 
          $exp_substr = "exon coding sequence";
        }
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
        
        my $expl_str = "sequence";
        if($is_matpept) { $expl_str = "mat_peptide"; }
        if($is_cds)     { $expl_str = "exon"; }
        if($do_olap) { 
          $tok4 = sprintf(" %10s", sprintf("%s%s", "overlaps", $mdl_exon_idx+1));
          $exp_tok4 = "overlaps<j>";
          $tok5 = sprintf(" %10s", "----------");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
          if($do_model_explanation) { 
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, sprintf("'P' or 'F' followed by list of features this feature overlaps with $expl_str"), $FH_HR);
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
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $exp_tok4, undef, sprintf("'P' or 'F' followed by list of mat_peptides this mat_peptide is adjacent with $expl_str"), $FH_HR);
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "first letter is 'P' if agrees exactly with reference, else 'F'", $FH_HR); # adds a second line to explanation
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "\"NP\" if no prediction", $FH_HR);      
          }
          $need_to_define_H{"adjacent"} = 1;
        }

        $exp_substr = "sequence";
        if($is_matpept) { $exp_substr = "mat_peptide coding sequence"; }
        if($is_cds)     { $exp_substr = "CDS"; }
        if($mdl_info_HAR->{"is_final"}[$mdl_idx]) { 
          $tok4 = sprintf(" %6s", "length");
          $tok5 = sprintf(" %6s", "------");
          output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4); 
          if($do_model_explanation) { 
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, sprintf("length of $exp_substr #<i> (all %s summed)", $is_cds ? "exons" : "segments"), $FH_HR);
          }      

          if($is_cds && $do_ss3) { # skip this in matpept mode, we don't check start/stop of mat_peptides, only CDS, later
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
          if($is_cds && $do_stop) { # skip this in matpept mode, we only check stop of final mat_peptide, later
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
            elsif($is_cds) { 
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
            else { 
              my $cur_type = $ftr_info_HAR->{"type_ftable"}[$ftr_idx];
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, $exp_tok1, $tok4, undef, "annotation indicating if this $cur_type PASSED ('P') or FAILED ('F')", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  a $cur_type sequence PASSES ('P') if and only if all of the following", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  conditions are met (else it FAILS):", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (1) all of its segments have a pairwise alignment to the homologous", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      reference exon that extends to the 5' and 3' boundary of the", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      reference annotation.", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "  (5) all of its segments overlap with exact same set of other segments as the", $FH_HR);
              output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, "      homologous reference $cur_type", $FH_HR);
              push(@pf_text_A, sprintf("P/F characters %d to %d pertain to each of the %d $cur_type, in order.", $pf_idx, $pf_idx + $nftr-1, $nftr));
              $pf_idx += $nftr;
              $need_to_define_H{"overlap"} = 1;
            }
            output_tbl_get_headings_explanation_helper($out_header_exp_AR, undef, undef, undef, undef, $FH_HR);
          }
        } # end of 'if($mdl_is_final_AR->[$mdl_idx])'
        $do_model_explanation = 0; # once we see the final exon of the first CDS, we don't need to print CDS explanations anymore
      } # end of 'for(my $mdl_idx)'
    } # end of 'else' entered if we're not a multifeature cds-mp feature but do have annot_type eq 'model'
    ################################################
    # block that handles 'annot_type' eq "duplicate"
    ################################################
    elsif($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "duplicate") { 
      $width = 6 + 1 + 6 + 1 + 6; #20
      $tok1     = sprintf("  %*s", $width, $ftr_out_short . getMonocharacterString(int($width-length($ftr_out_short)/2), " ", $FH_HR));
      $exp_tok1 = $ftr_out_short;
      $tok2 = sprintf("  %s", $ftr_out_product); 
      $tok3 = sprintf("  %s", getMonocharacterString($width, "-", $FH_HR));
      $tok4 = sprintf("  %8s", sprintf("%s", "start"));
      $tok5 = sprintf("  %8s", "--------");
      
      output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
      
      $tok4 = sprintf(" %8s", "stop");
      $tok5 = sprintf(" %8s", "------");
      output_tbl_get_headings_helper($out_row_header_AR,  $row_div_char, $tok1, $tok2, $tok4);
    } # end of block that handles 'annot_type' eq 'duplicate'
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
    my $seq_len   = $seq_info_HAR->{"len"}[$seq_idx];
    my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
    my @cur_out_A = (); # array of current tokens to print
    my $ngenbank_match = 0; # number of matches with existing annotation
    my $pass_fail_str  = ""; 

    # Create the initial portion of the output line, the accession and length
    push(@cur_out_A, sprintf("%-5d  ", ($seq_idx+1)));
    push(@cur_out_A, sprintf("%-19s  ", $accn_name)); 
    push(@cur_out_A, sprintf("%6d ", $seq_len));

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
                                                              $seq_len, $FH_HR);
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
          if($cur_start == $seq_len) { 
            push(@cur_out_A, sprintf("  %6d", 0)); # start 
            push(@cur_out_A, sprintf("  %6d", 0)); # stop 
            push(@cur_out_A, sprintf("  %6d", 0)); # length
          }
          else { # 1st feature does not start at nt $seq_len on negative strand
            push(@cur_out_A, sprintf("  %6d", $seq_len)); # start 
            push(@cur_out_A, sprintf("  %6d", $cur_start + 1)); # stop
            push(@cur_out_A, sprintf("  %6d", $seq_len - $cur_start)); # length
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
      elsif($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { 
        for(my $mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx]; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
          my $is_first = $mdl_info_HAR->{"is_first"}[$mdl_idx]; # is this the first model for feature $ftr_idx?
          my $is_final = $mdl_info_HAR->{"is_final"}[$mdl_idx]; # is this the final model for feature $ftr_idx?
          my $is_matpept = featureTypeIsMaturePeptide($ftr_info_HAR->{"type"}[$ftr_idx]);
          my $is_cds     = featureTypeIsCds($ftr_info_HAR->{"type"}[$ftr_idx]);
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
              if($is_matpept)  { push(@cur_out_A, sprintf(" %10s", "-")); }  # adj
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
              if($is_cds) { 
                if(validateStopCodon($mdl_results_HR->{"out_stop_codon"})) { 
                  $stop_codon_char = $ss3_yes_char;
                }
                else { 
                  $stop_codon_char   = $ss3_no_char;
                  $at_least_one_fail = 1;
                }
              } # end of 'if($is_cds)'
              if(($mdl_results_HR->{"cumlen"} % 3) == 0) { 
                $multiple_of_3_char = $ss3_yes_char;
              }
              else { 
                $multiple_of_3_char = $ss3_yes_char;
                $at_least_one_fail = 1;
              }
              push(@cur_out_A, sprintf(" %6d", $mdl_results_HR->{"cumlen"})); 
              
              # add the ss3 (start/stop/multiple of 3 info) if we're not a mature peptide
              if($is_cds) { 
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
              if($is_cds && $do_ss3)  { push(@cur_out_A, "  NP"); } # ss3
              if($is_cds && $do_stop) { push(@cur_out_A, sprintf(" %3s", "NP")); } # stop
              $pass_fail_char = "F";
              push(@cur_out_A, sprintf(" %2s", $pass_fail_char));
              $pass_fail_str .= $pass_fail_char;
            }
            if($is_first) { # important to do this so final model has a valid start_codon_char
              $start_codon_char = "NP";
            }
          }
        } # end of 'for(my $mdl_idx'
      } # end of 'elsif($ftr_info_HAR->{"annot_type"} eq "model") {' entered if feature is not a multifeature cds-mp but has annot_type == "model"
      elsif($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "duplicate") {
        my $src_idx = $ftr_info_HAR->{"source_idx"}[$ftr_idx];
        if($src_idx == -1) {
          DNAORG_FAIL("ERROR in $sub_name, feature index $ftr_idx has annot_type of duplicate, but has source_idx of -1", 1, $ofile_info_HHR->{"FH"});
        }
        my $start_pos = "-";
        my $stop_pos  = "-";
        if(($ftr_info_HAR->{"annot_type"}[$src_idx] eq "multifeature") && 
           ($ftr_info_HAR->{"type"}[$src_idx]       eq "cds-mp")) {
          $start_pos = $ftr_results_AAHR->[$src_idx][$seq_idx]->{"out_start"};
          $stop_pos  = $ftr_results_AAHR->[$src_idx][$seq_idx]->{"out_stop"};
        }
        elsif($ftr_info_HAR->{"annot_type"}[$src_idx] eq "model") {
          my $first_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$src_idx];
          my $final_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$src_idx];
          if(exists $mdl_results_AAHR->[$first_mdl_idx][$seq_idx]->{"out_start"}) { 
            $start_pos = $mdl_results_AAHR->[$first_mdl_idx][$seq_idx]->{"out_start"};
          }
          if(exists $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]->{"out_stop"}) { 
            $stop_pos = $mdl_results_AAHR->[$final_mdl_idx][$seq_idx]->{"out_stop"};
          }
        }
        push(@cur_out_A, sprintf("  %8s ", $start_pos)); # start position
        push(@cur_out_A, sprintf("%8s",    $stop_pos));  # stop position
      }        
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
                                                            $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
        }
        elsif(exists $mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"c_stop"}) { 
          (undef, $cur_stop) = create_output_start_and_stop($mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_start"},  # irrelevant due to the first undef arg
                                                            $mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"c_stop"},
                                                            $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
        }
        else { 
          (undef, $cur_stop) = create_output_start_and_stop($mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_start"}, # irrelevant due to the first undef arg
                                                            $mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_stop"}, 
                                                            $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
        }
        if($mdl_results_AAHR->[($nmdl-1)][$seq_idx]{"p_strand"} eq "+") { 
          # positive strand, easy case
          if($cur_stop == $seq_len) { # final model prediction stops at final nt
            push(@cur_out_A, sprintf("  %6d", 0)); # start 
            push(@cur_out_A, sprintf("  %6d", 0)); # stop 
            push(@cur_out_A, sprintf("  %6d", 0)); # length
          }            
          else { 
            push(@cur_out_A, sprintf("  %6d", $cur_stop + 1));          # start
            push(@cur_out_A, sprintf("  %6d", $seq_len));               # stop 
            push(@cur_out_A, sprintf("  %6d", ($seq_len - $cur_stop))); # length
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
    push(@cur_out_A, sprintf("  %6d", $seq_len));
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
# Subroutine: output_feature_tbl_all_sequences()
# Incept:     EPN, Tue Dec  5 13:49:17 2017 [rewritten Tue Oct 30 05:59:04 2018]
# Purpose:    Output the feature table for all sequences.
#
# Arguments:
#  $err_ftr_instances_AHHR:  REF to array of 2D hashes with per-feature errors, PRE-FILLED
#  $err_seq_instances_HHR:   REF to 2D hash with per-sequence errors, PRE-FILLED
#  $mdl_info_HAR:            REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:            REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:            REF to hash of arrays with information on the sequences, PRE-FILLED
#  $err_info_HAR:            REF to the error info hash of arrays, PRE-FILLED
#  $mdl_results_AAHR:        REF to model results AAH, PRE-FILLED
#  $ftr_results_AAHR:        REF to feature results AAH, PRE-FILLED
#  $opt_HHR:                 REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:          REF to the 2D hash of output file information
#             
# Returns:  Number of sequences that 'pass'.
# 
# Dies:     never
#
#################################################################
sub output_feature_tbl_all_sequences { 
  my $sub_name = "output_feature_tbl_all_sequences";
  my $nargs_exp = 11;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_ftr_instances_AHHR, $err_seq_instances_HHR, $mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, 
      $err_info_HAR, $mdl_results_AAHR, $ftr_results_AAHR,
      $class_errors_per_seq_HR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience
  my $pass_ftbl_FH    = $FH_HR->{"pass_ftbl"};   # feature table for PASSing sequences
  my $fail_ftbl_FH    = $FH_HR->{"fail_ftbl"};   # feature table for FAILing sequences
  my $long_ftbl_FH    = $FH_HR->{"long_ftbl"};   # long feature table for all sequences
  my $pass_list_FH    = $FH_HR->{"pass_list"};   # list of PASSing seqs
  my $fail_list_FH    = $FH_HR->{"fail_list"};    # list of FAILing seqs
  my $fail_co_list_FH = $FH_HR->{"fail_co_list"}; # list of FAILing seqs, that would have passed if there were no class errors
  my $errors_FH       = $FH_HR->{"errors_list"}; # list of errors 
  print $errors_FH "#sequence\terror\tfeature\terror-description\n";

  my $ret_npass = 0;  # number of sequences that pass, returned from this subroutine

  # validate input and get counts of things
  my $nmdl = validateModelInfoHashIsComplete    ($mdl_info_HAR, undef, $FH_HR); # nmdl: number of homology models
  my $nftr = validateFeatureInfoHashIsComplete  ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nseq = validateSequenceInfoHashIsComplete ($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences
  my $nerr = getConsistentSizeOfInfoHashOfArrays($err_info_HAR, $FH_HR); # nerr: number of different error codes

  # define the hard-coded type priority hash, which defines the order of types in the feature table output, lower is higher priority
  my %type_priority_H = ();
  $type_priority_H{"gene"}         = 0;
  $type_priority_H{"CDS"}          = 1;
  $type_priority_H{"misc_feature"} = 2;
  $type_priority_H{"mat_peptide"}  = 3;
  my $npriority = scalar(keys %type_priority_H);
  # reference for the way we sort the information we collect for the feature table
  #https://stackoverflow.com/questions/10395383/sorting-an-array-of-hash-by-multiple-keys-perl      

  my $qval_sep = ";;"; # value separating multiple qualifier values in a single element of $ftr_info_HAR->{$key}[$ftr_idx]
  # NOTE: $qval_sep == ';;' is hard-coded value for separating multiple qualifier values for the same 
  # qualifier (see dnaorg.pm::edirectFtableOrMatPept2SingleFeatureTableInfo

  # main loop: for each sequence
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
    my $seq_len   = $seq_info_HAR->{"len"}[$seq_idx];
    my $accn_name = $seq_info_HAR->{"accn_name"}[$seq_idx];
    
    # class for this sequence, possible values:

    my @ftout_AH      = (); # array of hashes with output for feature table, kept in a hash so we can sort before outputting
    my $ftidx         = 0;  # index in @ftout_AH
    my $min_coord     = -1; # minimum coord in this feature
    my $cur_min_coord = -1; # minimum coord in this segment
    my %fidx2idx_H    = (); # key is feature index $fidx, value is $ftidx index in @ftout_AH that $fidx corresponds to
    my $i;

    my @seq_error_A = (); # all errors for this sequence
    my @seq_note_A  = (); # all notes for this sequence

    my $missing_codon_start_flag = 0; # set to 1 if a feature for this sequence should have a codon_start value added but doesn't

    # first check for per-sequence errors
    my $seq_err_str = helper_ftable_get_seq_error_code_strings($seq_name, $err_seq_instances_HHR, $err_info_HAR, $FH_HR);
    processSequenceErrorsForFTable($seq_err_str, $seq_name, $err_info_HAR, $err_seq_instances_HHR, \@seq_error_A, $FH_HR);

    # loop over features
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      my $feature_type = $ftr_info_HAR->{"type_ftable"}[$ftr_idx];
      my $orig_feature_type = $feature_type; 

      # determine if this feature can be ignored or cannot be ignored (should be annotated), 
      # it can be ignored if:
      # - we do not have a defined "p_start" for any model associated with this feature
      #   (this is equivalent to having a 'nop' error for all models associated with this feature)
      # - we do not have a 'xnn' error for this feature
      my $is_duplicate  = ($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "duplicate") ? 1 : 0;
      my $dup_src_ftidx = undef; # for duplicate features, set to index in $ftout_AH that corresponds to source for this feature

      my $defined_pstart = 0;
      my $xnn_flag       = 0;
      my $do_ignore      = 1; 
      if($is_duplicate) { 
        if($ftr_info_HAR->{"source_idx"}[$ftr_idx] == -1) { 
          DNAORG_FAIL("ERROR in $sub_name, feature index $ftr_idx has annot_type of duplicate, but has source_idx of -1", 1, $FH_HR);
        }
        my $src_idx = $ftr_info_HAR->{"source_idx"}[$ftr_idx];
        if(exists $fidx2idx_H{$src_idx}) { 
          $dup_src_ftidx = $fidx2idx_H{$src_idx}; 
          $do_ignore = 0;
        }
      }
      else { # not a duplicate feature
        $defined_pstart = check_for_defined_pstart_in_mdl_results($seq_idx, $ftr_idx, $ftr_info_HAR, $mdl_results_AAHR, $FH_HR);
        $xnn_flag       = (exists $err_ftr_instances_AHHR->[$ftr_idx]{"xnn"}{$seq_name}) ? 1 : 0;
        $do_ignore      = ($defined_pstart || $xnn_flag) ? 0 : 1;
      }

      if(! $do_ignore) { 
        # per-feature values that are modified if we have an exception that covers all errors for this feature
        my $do_misc_feature = 0;
        my $do_start_carrot = 0;
        my $do_stop_carrot  = 0;
        my $do_pred_stop    = 0;

        my $ftr_coords_str = ""; # string of coordinates for this feature
        my $ftr_out_str    = ""; # output string for this feature
        my $long_out_str   = ""; # long feature table output string for this feature
        my @ftr_long_output_A = (); # array of strings with complete error messages, for long feature table
        my @ftr_note_A    = ();  # notes for this feature/sequence combination
        my $ftr_tiny = $ftr_info_HAR->{"out_tiny"}[$ftr_idx];

        if(! $is_duplicate) { # only look up errors if we're not a duplicate feature
          # fill an array and strings with all errors for this sequence/feature combo
          my $ftr_err_str = helper_ftable_get_ftr_error_code_strings($seq_name, $ftr_idx, $err_ftr_instances_AHHR, $err_info_HAR, \@ftr_long_output_A, $FH_HR);
          $do_start_carrot = ($ftr_err_str =~ m/b5e/) ? 1 : 0;
          $do_stop_carrot  = ($ftr_err_str =~ m/b3e/) ? 1 : 0;
          $do_pred_stop = processFeatureErrorsForFTable($ftr_err_str, $seq_name, $ftr_idx, $ftr_info_HAR, $err_info_HAR, $err_ftr_instances_AHHR, 
                                                        \@ftr_note_A, \@seq_error_A, $FH_HR);
          if(scalar(@ftr_note_A) > 0) { 
            $do_misc_feature = 1;
            $feature_type = "misc_feature";
            push(@seq_note_A, @ftr_note_A);
          }
          else { 
            $do_misc_feature = 0; 
          }
        } # end of 'if(! $is_duplicate)'

        # determine coordinates for the feature differently depending on:
        # if we are: 
        # - a duplicate ($is_duplicate)
        # - we have at least one prediction ($defined_pstart) 
        # - not a duplicate and don't have at least one prediction (in this case, $xnn_flag will be up, or else we would have ignored this feature)
        if($is_duplicate) { 
          $ftr_coords_str  = $ftout_AH[$dup_src_ftidx]{"coords"};
          $min_coord       = $ftout_AH[$dup_src_ftidx]{"mincoord"};
          $do_start_carrot = $ftout_AH[$dup_src_ftidx]{"carrot"};
        }
        elsif(! $defined_pstart) { 
          # $xnn_flag must be TRUE
          if(! $xnn_flag) { # sanity check
            DNAORG_FAIL("ERROR in $sub_name, sequence $accn_name feature $ftr_idx $ftr_tiny is not being ignored but defined_pstart and xnn_flag are both FALSE - shouldn't happen.", 1, $ofile_info_HHR->{"FH"});
          }
          $ftr_coords_str = helper_ftable_get_coords_xnn_flag($seq_idx, $ftr_idx, $do_start_carrot, $do_stop_carrot, \$min_coord, $seq_info_HAR, $ftr_results_AAHR, $FH_HR);
        }
        else { # $is_duplicate is '0' and $defined_pstart is '1'
          # this function will create coordinates differently depending on feature type, do_pred_stop, etc.
          $ftr_coords_str = helper_ftable_get_coords_standard($seq_idx, $ftr_idx, $do_start_carrot, $do_stop_carrot, $do_pred_stop, \$min_coord, 
                                                              $mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $ftr_results_AAHR, $FH_HR);
        }
        # convert coordinate string to output string
        $ftr_out_str = helper_ftable_coords_to_out_str($ftr_coords_str, $feature_type, $FH_HR);
        
        # add qualifiers: product, gene, exception and codon_start (if !duplicate)
        if(! $do_misc_feature) { 
          $ftr_out_str .= helper_ftable_add_qualifier_from_ftr_info($ftr_idx, "out_product",   $qval_sep, $ftr_info_HAR, $FH_HR);
          $ftr_out_str .= helper_ftable_add_qualifier_from_ftr_info($ftr_idx, "out_gene",      $qval_sep, $ftr_info_HAR, $FH_HR);
          $ftr_out_str .= helper_ftable_add_qualifier_from_ftr_info($ftr_idx, "out_exception", $qval_sep, $ftr_info_HAR, $FH_HR);
          # check for existence of "x_frame" value for all CDS, but only actually output them if we have a start_carrot
          if((! $is_duplicate) && 
             (($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-mp") || ($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-notmp"))) { 
            my $tmp_str = helper_ftable_add_qualifier_from_ftr_results($seq_idx, $ftr_idx, "x_frame", "codon_start", $ftr_results_AAHR, $FH_HR);
            if($tmp_str eq "") { 
              # we didn't have an x_frame value for this CDS, so raise a flag
              # we check later that if the sequence passes that this flag 
              # is *NOT* raised, if it is, something went wrong and we die
              $missing_codon_start_flag = 1; 
            } 
            if($do_start_carrot) { # only add the codon_start if we have a start_carrot
              $ftr_out_str .= $tmp_str;
            }
          }
        }
        # add notes and full error messages (if !duplicate)
        if(! $is_duplicate) { 
          foreach my $note_value (@ftr_note_A) { 
            $ftr_out_str .= sprintf("\t\t\t%s\t%s\n", "note", $note_value);
          }
          # add the full error messages as 'notes' to the long output string, which will be output to the long feature table
          foreach my $long_line (@ftr_long_output_A) { 
            $long_out_str .= sprintf("\t\t\t%s\t%s\n", "note", $long_line);
          }
        }

        # push to the output hash
        %{$ftout_AH[$ftidx]} = ();
        $ftout_AH[$ftidx]{"carrot"}          = $do_start_carrot;
        $ftout_AH[$ftidx]{"mincoord"}        = $min_coord;
        $ftout_AH[$ftidx]{"type_priority"}   = (exists $type_priority_H{$orig_feature_type}) ? $type_priority_H{$orig_feature_type} : $npriority;
        $ftout_AH[$ftidx]{"coords"}          = $ftr_coords_str;
        $ftout_AH[$ftidx]{"output"}          = $ftr_out_str;
        $ftout_AH[$ftidx]{"long_output"}     = $ftr_out_str . $long_out_str;
        $fidx2idx_H{$ftr_idx} = $ftidx;
        $ftidx++;
      } # end of 'if(! $do_ignore)'
    } # end of 'for(my $ftr_idx...'

    #######################################
    # OUTPUT section 
    #######################################
    # done with this sequence, determine what type of output we will have 
    my $has_class_errors = 0;
    if((defined $class_errors_per_seq_HR) &&
       (exists $class_errors_per_seq_HR->{$seq_name})) { 
      $has_class_errors = 1; 
    }

    my $cur_noutftr  = scalar(@ftout_AH);
    my $cur_nerror   = scalar(@seq_error_A);
    my $cur_nnote    = scalar(@seq_note_A);

    # sort output
    if($cur_noutftr > 0) { 
      @ftout_AH = sort { $a->{"mincoord"}      <=> $b->{"mincoord"} or 
                             $b->{"carrot"}        <=> $a->{"carrot"}   or
                             $a->{"type_priority"} <=> $b->{"type_priority"} 
      } @ftout_AH;
    }              

    # sequences only pass if:
    # - at least one feature is annotated ($cur_noutftr > 0)
    # - zero notes and errors
    my $do_pass = (($cur_noutftr > 0) && ($cur_nnote == 0) && ($cur_nerror == 0) && ($has_class_errors == 0)) ? 1 : 0;

    # sanity check, if we have no notes, errors and didn't skip anything, we should also have set codon_start for all features
    if($cur_noutftr > 0 && $cur_nnote == 0 && $cur_nerror == 0 && ($missing_codon_start_flag)) { 
      DNAORG_FAIL("ERROR in $sub_name, sequence $accn_name set to PASS, but at least one CDS had no codon_start set - shouldn't happen.", 1, $ofile_info_HHR->{"FH"});
    }
              
    if($do_pass) { 
      # print to both the passing feature table file and the long feature table file
      $ret_npass++;
      print $pass_list_FH $accn_name . "\n";
      print $pass_ftbl_FH ">Feature $accn_name\n";
      print $long_ftbl_FH ">Feature $accn_name\n";
      for($i = 0; $i < scalar(@ftout_AH); $i++) { 
        # print 
        print $pass_ftbl_FH $ftout_AH[$i]{"output"};
        print $long_ftbl_FH $ftout_AH[$i]{"long_output"};
      }
    }
    else { # $do_pass == 0
      print $fail_list_FH $accn_name . "\n";
      print $fail_ftbl_FH ">Feature $accn_name\n";
      print $long_ftbl_FH ">Feature $accn_name\n";
      for($i = 0; $i < scalar(@ftout_AH); $i++) { 
        # print 
        print $fail_ftbl_FH $ftout_AH[$i]{"output"};
        print $long_ftbl_FH $ftout_AH[$i]{"long_output"};
      }
      if(($cur_nerror > 0) || ($has_class_errors)) { 
        print $fail_ftbl_FH "\nAdditional note(s) to submitter:\n"; 
        print $long_ftbl_FH "\nAdditional note(s) to submitter:\n"; 
        for(my $e = 0; $e < scalar(@seq_error_A); $e++) { 
          my $error_line = $seq_error_A[$e];
          print $fail_ftbl_FH "ERROR: " . $error_line . "\n"; 
          print $long_ftbl_FH "ERROR: " . $error_line . "\n"; 
          if($error_line =~ /([^\:]+)\:\s\(([^\)]+)\)\s*(.+)$/) {
            print $errors_FH ($accn_name . "\t" . $1 . "\t" . $2 . "\t" . $3 . "\n");
          }
          else {
            DNAORG_FAIL("ERROR in $sub_name, unable to split error_line for output: $error_line", 1, $ofile_info_HHR->{"FH"});
          }
        }
        if($has_class_errors) { 
          print $errors_FH ($class_errors_per_seq_HR->{$seq_name});
          if(($cur_nerror == 0) && (defined $fail_co_list_FH)) { 
            print $fail_co_list_FH $accn_name . "\n";
          }
          my $ftable_error_str = error_list_output_to_ftable_errors($seq_name, $class_errors_per_seq_HR->{$seq_name}, $ofile_info_HHR->{"FH"});
          print $fail_ftbl_FH $ftable_error_str;
          print $long_ftbl_FH $ftable_error_str;
        }
      } # end of 'if($cur_nerror > 0) || $has_class_errors'
    }
  } # end of loop over sequences

  return $ret_npass;
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
          printf $all_FH ("%-10s  %3s  %-9s  %4s  %s%s\n", $accn_name, "N/A", "N/A", $err_code, $err_info_HAR->{"desc"}[$err_idx], 
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
            printf $all_FH ("%-10s  %3s  %-9s  %4s  %s%s\n", $accn_name, ($ftr_idx+1), $out_tiny, $err_code, $err_info_HAR->{"desc"}[$err_idx], 
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
      my $is_cds     = featureTypeIsCds($ftr_info_HAR->{"type"}[$ftr_idx]);
      my $is_matpept = featureTypeIsMaturePeptide($ftr_info_HAR->{"type"}[$ftr_idx]);
      my $type_idx   = $ftr_info_HAR->{"type_idx"}[$ftr_idx];
      my $w_cur = $w_tot_gap_length_A[$ftr_idx] + 2 + $w_net_gap_length_A[$ftr_idx] + 2 + $w_gapstr_A[$ftr_idx];
      if($ftr_idx > 0) { print $perseq_FH "  "; }
      printf $perseq_FH ("%-*s", $w_cur, $ftr_info_HAR->{"type_ftable"} . "#");
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
    my $is_cds     = featureTypeIsCds($ftr_info_HAR->{"type"}[$ftr_idx]);
    my $is_matpept = featureTypeIsMaturePeptide($ftr_info_HAR->{"type"}[$ftr_idx]);
    my $type_idx   = $ftr_info_HAR->{"type_idx"}[$ftr_idx];
    if($ftr_idx > 0) { print $perseq_FH "  "; }
    if(! $do_gap_special) { 
      printf $perseq_FH ("%-*s  ", $w_tot_gap_length_A[$ftr_idx], $ch_tot_gap_length);
      printf $perseq_FH ("%-*s  ", $w_net_gap_length_A[$ftr_idx], $ch_net_gap_length);
      printf $perseq_FH ("%-*s", $w_gapstr_A[$ftr_idx], $ch_gapstr);
    }
    else { 
      printf $perseq_FH ("%-*s", $w_gapstr_A[$ftr_idx], $ftr_info_HAR->{"type_ftable"} . "#");
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
    my $is_cds     = featureTypeIsCds($ftr_info_HAR->{"type"}[$ftr_idx]);
    my $is_matpept = featureTypeIsMaturePeptide($ftr_info_HAR->{"type"}[$ftr_idx]);
    if((scalar(keys %{$ftr_gapstr_AH[$ftr_idx]})) > 0) { 
      foreach my $key (sort keys %{$ftr_gapstr_AH[$ftr_idx]}) { 
        printf $pergap_FH ("%s#" . ($ftr_idx+1) . " " . $key . " " . $ftr_gapstr_AH[$ftr_idx]{$key} . "\n", $ftr_info_HAR->{"type_ftable"});
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
# Subroutine:  translate_feature_sequences
# Incept:      EPN, Thu Mar 10 11:13:25 2016
#
# Purpose:    For each file with corrected feature sequences, 
#             translate them into proteins.
#
# Arguments: 
#  $esl_epn_translate: executable for esl-epn-translate.pl
#  $in_key:            key for the input files we'll translate here, usually "corrected"
#  $out_key:           key for the output files we'll create here, usually "corrected.translated",
#                      this key will be stored in $ftr_info_HAR
#  $specstart_AAR:     REF to 2D array of specified, permissible start codons, can be undef
#  $ftr_info_HAR:      REF to hash of arrays with information on the features, ADDED TO HERE
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
################################################################# 
sub translate_feature_sequences { 
  my $sub_name = "translate_feature_sequences()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($esl_epn_translate, $in_key, $out_key, $specstart_AAR, $ftr_info_HAR, $ofile_info_HHR) = @_;

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
      runCommand($cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HH{"FH"});
      
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
    my $mxsize_opt = sprintf("--mxsize %d", opt_Get("--amxsize", $opt_HHR));
    my $cmd = $execs_HR->{"cmalign"} . " --cpu 0 $mxsize_opt $cm_file $fa_file > $stk_file";
    runCommand($cmd, opt_Get("-v", \%opt_HH), 0, $ofile_info_HHR->{"FH"});
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
    if($ftr_info_HAR->{"annot_type"}[$ftr_idx] ne "duplicate") { 
      foreach my $file_type (@both_file_types_A, @ftr_only_file_types_A) { 
        $ftr_info_HAR->{$file_type}[$ftr_idx] = $out_root . "." . $ftr_info_HAR->{"filename_root"}[$ftr_idx] . "." . $file_type;
      }
    }
    else { 
      foreach my $file_type (@both_file_types_A, @ftr_only_file_types_A) { 
        $ftr_info_HAR->{$file_type}[$ftr_idx] = "/dev/null";
      }
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
      runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"}); 
      
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
      runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
      
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
#             in coordinate space 1..accn_len. If seq_len == accn_len,
#             which will be the case if the -c option was not used (that is,
#             if the genome is not circular and thus we did not duplicate it)
#             we do nothing. If seq_len == accn_len * 2, then we have to 
#             check if out_start and out_stop should be different from p_start
#             and p_stop.
#
# Arguments: 
#  $in_start:   predicted start (from our annotation)
#  $in_stop:    predicted (possibly corrected) stop (from our annotation)
#  $seq_len:    length of the sequence we searched, either ($accn_len) or (2 * $accn_len)
#  $FH_HR:      REF to hash of file handles
#
# Returns:    Two values:
#             $out_start: the start coordinate to output
#             $out_stop:  the stop coordinate to output
#
# Dies: never
#
################################################################# 
sub create_output_start_and_stop { 
  my $sub_name = "create_output_start_and_stop()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($in_start, $in_stop, $seq_len, $FH_HR) = @_;
  
  # contract check
  if((! verify_integer($in_start)) || ($in_start <= 0)) { 
    DNAORG_FAIL("ERROR in $sub_name, input start coordinate is expected to be a positive integer, got $in_start", 1, $FH_HR);
  }
  if((! verify_integer($in_stop)) || ($in_stop <= 0)) { 
    DNAORG_FAIL("ERROR in $sub_name, input stop coordinate is expected to be a positive integer, got $in_stop", 1, $FH_HR);
  }

  my ($out_start, $out_stop);
  $out_start = $in_start;
  $out_stop  = $in_stop;

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
#            $seq_info_HAR->{"len"}[$i]
#            for sequences $i=1..N.
#
# Arguments:
#  $infasta_file:  name of the input fasta file
#  $out_root:      root for naming output files
#  $seq_info_HAR:  REF to hash of arrays of sequence information, added to here
#  $opt_HHR:       REF to 2D hash of option values, see top of epn-options.pm for description
#  $FH_HR:         REF to hash of file handles
#  
# Returns:  Number of sequences read in $infasta_file.
# 
# Dies: If $seq_info_HAR->{"len"} and $seq_info_HAR->{"accn_name"} 
#       are not both arrays of exactly length 1 (with information on 
#       *only* the reference accession.)
#       
#       If the same sequence exists more than once in the input sequence
#       file.
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
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($infasta_file, $seq_info_HAR, $opt_HHR, $FH_HR) = @_;

  my %accn_exists_H = ();  # keeps track of which accessions have been read from the sequence file
  my $err_flag = 0;        # changed to '1' if an accession exists more than once in the input fasta file

  my $ssi_file = $infasta_file . ".ssi";
  if(-e $ssi_file) { 
    runCommand("rm $ssi_file", opt_Get("-v", $opt_HHR), 0, $FH_HR);
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
      push(@{$seq_info_HAR->{"len"}}, $len);

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


  # parse the consopts file
  my %consopts_used_H    = (); # key option used in dnaorg_buid.pl from consopts file, value argument used in dnaorg_build.pl
  my %consopts_notused_H = (); # key option in consopts file, value argument used in dnaorg_build.pl
  my %consmd5_H  = ();         # key option used in dnaorg_build.pl in consopts file, md5 checksum value of the file name argument used in dnaorg_build.pl
  parseConsOptsFile($consopts_file, \%consopts_used_H, \%consopts_notused_H, \%consmd5_H, $FH_HR);

  # make sure options are consistent with what we read in the consopts file
  my $opt;
  my $optfile;
  my $optfile_md5;
  my $optarg;
  foreach $opt (sort keys (%consopts_used_H)) { 
    if(! opt_IsUsed($opt, $opt_HHR)) { 
      DNAORG_FAIL("ERROR, the $opt option was used when dnaorg_build.pl was run (according to file $consopts_file).\nYou must also use it with dnaorg_annotate.pl.", 1, $FH_HR);
    }
    # option was used in both dnaorg_build.pl and dnaorg_annotate.pl, 
    # if it has a consmd5 value, check those are the same (in those 
    # cases we don't require argument is identical (files can have different names
    # as long as their md5s are identical)
    if($consmd5_H{$opt} ne "") { 
      my $optfile = opt_Get($opt, $opt_HHR);
      if(! -s $optfile) { 
        DNAORG_FAIL("ERROR, the file $optfile specified with the $opt option does not exist.", 1, $FH_HR);
      }          
      $optfile_md5 = md5ChecksumOfFile($optfile, $sub_name, $opt_HHR, $FH_HR);
      if($consmd5_H{$opt} ne $optfile_md5) { 
        DNAORG_FAIL("ERROR, the file $optfile specified with the $opt option does not appear to be identical to the file used\nwith dnaorg_build.pl. The md5 checksums of the two files differ: dnaorg_build.pl: " . $consmd5_H{$opt} . " dnaorg_annotate.pl: " . $optfile_md5, 1, $FH_HR);
      }
    }
    else { 
      # no md5 value, so we verify that option arguments are identical, if there is an argument
      if($consopts_used_H{$opt} ne "") { 
        $optarg = opt_Get($opt, $opt_HHR);
        if($consopts_used_H{$opt} ne $optarg) { 
          DNAORG_FAIL("ERROR, the option argument string $optarg specified with the $opt option does not appear to be identical to the argument string used\nwith the $opt option when dnaorg_build.pl was run, which was " . $consopts_used_H{$opt}, 1, $FH_HR);
        }
      }
    }
  }
  # all options that were used by dnaorg_build were also used by dnaorg_annotate,
  # now check that all options NOT used by dnaorg_build were also not used
  # by dnaorg_annotate

  foreach $opt (sort keys (%consopts_notused_H)) { 
    if(opt_IsUsed($opt, $opt_HHR)) { 
      DNAORG_FAIL("ERROR, the $opt option was not used when dnaorg_build.pl was run (according to file $consopts_file).\nYou must also not use it with dnaorg_annotate.pl, or you need to rerun dnaorg_build.pl with -c.", 1, $FH_HR);
    }
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
#  $opts .= " --nohmmonly --F1 0.02 --F2 0.001 --F2b 0.001 --F3 0.00001 --F3b 0.00001 --F4 0.0002 --F4b 0.0002 --F5 0.0002 --F6 0.0001 "; 
  $opts .= " --nohmmonly --FZ 30 ";

  my $cmd = $execs_HR->{"cmscan"} . " $opts $aorg_model $fasta_file > $cmscan_file";

  runCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

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

  runCommand($cmd, 0, 0, $ofile_info_HH{"FH"});
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
              create_output_start_and_stop($cur_ori_start_uapos, $cur_ori_stop_uapos, $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
          
          my $out_ori_offset_uapos;
          ($out_ori_offset_uapos, undef) = 
              create_output_start_and_stop($cur_ori_offset_uapos, $cur_ori_stop_uapos, $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
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
  
  my %seq_name_index_H = (); # seq_name_index_H{$seq_name} = <n>, means that $seq_name is the <n>th sequence name in @{$seq_info_HAR{*}} arrays
  getIndexHashForArray($seq_info_HAR->{"seq_name"}, \%seq_name_index_H, $FH_HR);

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

      my $seqidx = $seq_name_index_H{$seqname}; # sequence index for the hit in results_AAH (2nd dim of results_AAH)
      my $seqlen = $seq_info_HAR->{"len"}[$seqidx];

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
  my $seq_len = $seq_info_HAR->{"len"}[$seq_idx];

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
      $oseq_firstpos = $seq_len + $oseq_firstpos + 1;
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
    if($oseq_offset > ($seq_len / 2)) { # simpler (shorter distance) to move origin clockwise
      $oseq_offset = $seq_len - $oseq_offset; # note, we don't add 1 here
    }
    else { # simpler to shift origin in counterclockwise direction, we denote this as a negative offset
      $oseq_offset *= -1;
    }

    $oseq_passfail = "P";
  }

  return ($oseq_ct, $oseq_start, $oseq_stop, $oseq_firstpos, $oseq_offset, $oseq_passfail);
}

#################################################################
# Subroutine:  run_blastx_and_summarize_output()
# Incept:      EPN, Thu Oct  4 15:25:00 2018
#
# Purpose:    For each fasta file of predicted hits, run them as
#             blastx queries against the appropriate target blast DBs.
#
# Arguments: 
#  $execs_HR:          REF to a hash with "hmmbuild" and "hmmalign"
#                      executable paths
#  $out_root:          output root for the file names
#  $query_file:        query fasta file (input fasta file to script, same fasta file cmscan searches against)
#  $build_root:        path to build dir
#  $ftr_info_HAR:      REF to hash of arrays with information on the features
#  $compatibility_AHR: REF to array of hashes, 0..$f..$nftr-1, key: query sequence name in $blastx_query_file that 
#                      is 'compatible' with feature $f; that query sequence is an allowed hit for feature $f, FILLED HERE
#  $opt_HHR:           REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If blastx fails.
#
################################################################# 
sub run_blastx_and_summarize_output {
  my $sub_name = "run_blastx_and_summarize_output";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $out_root, $query_file, $build_root, $ftr_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  my $nftr = validateFeatureInfoHashIsComplete($ftr_info_HAR, undef, $ofile_info_HHR->{"FH"}); # nftr: number of features
  my $ncds = 0; 

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-mp") || 
       ($ftr_info_HAR->{"type"}[$ftr_idx] eq "cds-notmp")) { 
      $ncds++;
    }
  }

  # run blastx once on the full sequence file
  my $cur_db_file = $build_root . ".prot.fa";
  my $blastx_out_file = $out_root . ".blastx.out";
#  my $blastx_cmd = $execs_HR->{"blastx"} . " -query $query_file -db $cur_db_file -seg no -num_descriptions $ncds -num_alignments $ncds -out $blastx_out_file";
  my $blastx_cmd = $execs_HR->{"blastx"} . " -query $query_file -db $cur_db_file -seg no -out $blastx_out_file";
  runCommand($blastx_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $ofile_info_key = get_mdl_or_ftr_ofile_info_key("ftr", $ftr_idx, "predicted.hits.fa", $ofile_info_HHR->{"FH"});
    if(exists $ofile_info_HH{"fullpath"}{$ofile_info_key}) { 
      # printf("ftr_idx: $ftr_idx, out_tiny: %s ofile_info_key: $ofile_info_key %s\n", $ftr_info_HAR->{"out_tiny"}[$ftr_idx], $ofile_info_HHR->{"fullpath"}{$ofile_info_key});
      my $cur_query_file      = $ofile_info_HH{"fullpath"}{$ofile_info_key};
      my $cur_db_file         = $build_root . ".f" . $ftr_idx . ".prot.fa";
      my $cur_blastx_out_file = $out_root . ".f" . $ftr_idx . ".blastx.out";

      # run blast for this feature
      #$blastx_cmd = $execs_HR->{"blastx"} . " -query $cur_query_file -db $cur_db_file -seg no -num_descriptions 1 -num_alignments 1 -out $cur_blastx_out_file";
      $blastx_cmd = $execs_HR->{"blastx"} . " -query $cur_query_file -db $cur_db_file -seg no -out $cur_blastx_out_file";
      runCommand($blastx_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});

      # concatenate the blastx output for this feature to the growing blastx output for all blastx runs
      my $concat_cmd = "cat $cur_blastx_out_file >> $blastx_out_file";
      runCommand($concat_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
      #if(! opt_Get("--keep", $opt_HHR)) { 
      if(0) { 
        removeFileUsingSystemRm($cur_blastx_out_file, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"}); 
      }
    }
  }
  addClosedFileToOutputInfo($ofile_info_HHR, "blastx-out", $blastx_out_file, 0, "blastx output");

  # now summarize its output
  my $blastx_summary_file = $out_root . ".blastx.summary.txt";
  my $parse_cmd = $execs_HR->{"parse_blastx"} . " --input $blastx_out_file > $blastx_summary_file";
  runCommand($parse_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
  addClosedFileToOutputInfo($ofile_info_HHR, "blastx-summary", $blastx_summary_file, 0, "parsed (summarized) blastx output");

  return;
}

#################################################################
# Subroutine:  ftr_results_calculate_blastx()
# Incept:      EPN, Thu Oct  4 15:25:00 2018
#              [modified from subroutine of same name by Alejandro Schaffer, compare_predictions.pl]
#
# Purpose:    Parse blastx summary file and store results
#             in ftr_results.
#
# Arguments: 
#  $blastx_summary_file: path to blastx summary file to parse
#  $ftr_info_HAR:        REF to hash of arrays with information on the features
#  $seq_info_HAR:        REF to hash of sequence information
#  $seq_name_idx_HR;     REF to hash, key: $seq_name, value: idx of $seq_name in @{$seq_info_HAR->{"seq_name"}} and $ftr_info_HAR
#  $ftr_results_AAHR:    REF to feature results AAH, ADDED TO HERE
#  $opt_HHR:             REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:      REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If blastx fails.
#
################################################################# 
sub ftr_results_calculate_blastx { 
  my $sub_name = "ftr_results_calculate_blastx";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastx_summary_file, $ftr_info_HAR, $seq_info_HAR, $seq_name_idx_HR, $ftr_results_AAHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  validateFileExistsAndIsNonEmpty($blastx_summary_file, $sub_name, $FH_HR);
  my $nftr = validateAndGetSizeOfInfoHashOfArrays($ftr_info_HAR, undef, $FH_HR);

  open(IN, $blastx_summary_file) || fileOpenFailure($blastx_summary_file, $sub_name, $!, "reading", $FH_HR);
  # first step, determine which sequence each hit corresponds to
  my $query      = undef;
  my $target     = undef;
  my $hsp_idx    = undef;
  my $seq_idx    = undef; 
  my $ftr_idx    = undef;
  my $no_coords_query = undef; # name of query without coords, if query sequence is a predicted feature, e.g. "NC_002549.1/6039-8068" for query "NC_002549.1/6039-8068/10-885,885-2030"
  my $coords_str = undef; # string of coordinates if this is a predicted feature, e.g. "10-885,885-2030" for "NC_002549.1/6039-8068/10-885,885-2030"
  my $top_score_flag = 0;   # set to '1' if current hit is top scoring one for this sequence/feature pair

  while(my $line = <IN>) { 
    chomp $line;
    if($line ne "END_MATCH") { 
      my @el_A = split(/\t/, $line);
      if(scalar(@el_A) != 2) { 
        DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, did not read exactly 2 tab-delimited tokens in line $line", 1, $FH_HR);
      }
      my ($key, $value) = (@el_A);
      if($key eq "QACC") { 
        $query = $value;
        # determine what sequence it is
        my $matching_query = undef;
        ($matching_query, undef, undef) = helper_blastx_breakdown_query($query, undef, $seq_name_idx_HR, $FH_HR); 
        # helper_blastx_breakdown_query() will exit if $query is unparseable
        # two undefs: coords_str and query length, both irrelevant here
        $seq_idx = $seq_name_idx_HR->{$matching_query};
      }
      elsif($key eq "HACC") { 
        if(! defined $query) { 
          DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read HACC line before QACC line\n", 1, $FH_HR);
        }
        $target = $value;
        # determine what feature it is
        if($target =~ /(\S+)\/(\S+)/) { 
          my ($accn, $coords) = ($1, $2);
          # find it in @{$ftr_info_HAR->{"ref_coords"}}
          $ftr_idx = blastxDbSeqNameToFtrIdx($target, $ftr_info_HAR, $FH_HR); # will die if problem parsing $target, or can't find $ftr_idx
        }
      }
      elsif($key eq "HSP") { 
        if((! defined $query) || (! defined $ftr_idx)) { 
          DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read HSP line before one or both of QACC and HACC lines\n", 1, $FH_HR);
        }
        # printf("HEYA BLASTX HSP $key $value\n");
        if($value =~ /^(\d+)$/) { 
          $hsp_idx = $value;
          #printf("HEYA BLASTX set hsp_idx to $hsp_idx\n");
        }
        else { 
          DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary HSP line $line", 1, $FH_HR);
        }
      }
      elsif($key eq "QRANGE") { 
        if($value eq "..") { # special case, no hits, silently move on
          ;
        }
        elsif((! defined $query) || (! defined $ftr_idx) || (! defined $hsp_idx)) { 
          DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read QRANGE line before one or more of QACC, HACC, or HSP lines\n", 1, $FH_HR);
        }
        elsif($top_score_flag) { 
          if($value =~ /^(\d+)..(\d+)$/) { 
            my ($blast_start, $blast_stop) = ($1, $2);
            my $blast_strand = ($blast_start <= $blast_stop) ? "+" : "-";
            my ($out_start, $out_stop) = create_output_start_and_stop($blast_start, $blast_stop,
                                                                      $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
            # unless -c was used: $xstart will equal $blast_start and $xstop will equal $blast_stop
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_start"}  = $out_start;
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_stop"}   = $out_stop;
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_strand"} = $blast_strand;
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_query"}  = $query;
            #printf("HEYA BLASTX set ftr_results_AAHR->[$ftr_idx][$seq_idx]{x_start}  to " . $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_start"} . "\n");
            #printf("HEYA BLASTX set ftr_results_AAHR->[$ftr_idx][$seq_idx]{x_stop}   to " . $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_stop"} . "\n");
            #printf("HEYA BLASTX set ftr_results_AAHR->[$ftr_idx][$seq_idx]{x_strand} to " . $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_strand"} . "\n");
          }
          else { 
            DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary QRANGE line $line", 1, $FH_HR);
          }
        }
      }
      elsif($key eq "MAXIN") { 
        if((! defined $query) || (! defined $ftr_idx) || (! defined $hsp_idx)) { 
          DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read MAXIN line before one or more of QACC, HACC, or HSP lines\n", 1, $FH_HR);
        }
        if($top_score_flag) { 
          if($value =~ /^(\d+)$/) { 
            my $maxins = $1;
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_maxins"} = $maxins;
            #printf("HEYA BLASTX set ftr_results_AAHR->[$ftr_idx][$seq_idx]{x_maxins} to " . $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_maxins"} . "\n");
          }
          else { 
            DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary MAXIN line $line", 1, $FH_HR);
          }
        }
      }
      elsif($key eq "MAXDE") { 
        if((! defined $query) || (! defined $ftr_idx) || (! defined $hsp_idx)) { 
          DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read MAXDE line before one or more of QACC, HACC, or HSP lines\n", 1, $FH_HR);
        }
        if($top_score_flag) {
          if($value =~ /^(\d+)$/) { 
            my $maxdel = $1;
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_maxdel"} = $maxdel;
            #printf("HEYA BLASTX set ftr_results_AAHR->[$ftr_idx][$seq_idx]{x_maxdel} to " . $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_maxdel"} . "\n");
          }
          else { 
            DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary MAXDE line $line", 1, $FH_HR);
          }
        }
      }
      elsif($key eq "FRAME") { 
        if((! defined $query) || (! defined $ftr_idx) || (! defined $hsp_idx)) { 
          DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read FRAME line before one or more of QACC, HACC, or HSP lines\n", 1, $FH_HR);
        }
        if($top_score_flag) { 
          if($value =~ /^[\+\-]([123])$/) { 
            my $frame = $1;
            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_frame"} = $frame;
            #printf("HEYA BLASTX set ftr_results_AAHR->[$ftr_idx][$seq_idx]{x_frame} to " . $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_frame"} . "\n");
          }
          else { 
            DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary FRAME line $line ($key $value)", 1, $FH_HR);
          }
        }
      }
      elsif($key eq "STOP") { 
        if((! defined $query) || (! defined $ftr_idx) || (! defined $hsp_idx)) { 
          DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read STOP line before one or more of QACC, HACC, or HSP lines\n", 1, $FH_HR);
        }
        if($top_score_flag) { 
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_trcstop"} = $value;
          #printf("HEYA BLASTX set ftr_results_AAHR->[$ftr_idx][$seq_idx]{x_trcstop} to " . $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_trcstop"} . "\n");
        }
      }
      elsif($key eq "SCORE") { 
        if((! defined $query) || (! defined $ftr_idx) || (! defined $hsp_idx)) { 
          DNAORG_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read SCORE line before one or more of QACC, HACC, or HSP lines\n", 1, $FH_HR);
        }
        # is this sequence compatible and the highest scoring hit for this feature for this sequence? 
        if((! exists $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_score"}) || # first hit, so must be highest
           ($value > $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_score"})) { # highest scoring hit
          $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_score"} = $value;
          $top_score_flag = 1;
          #printf("HEYA BLASTX set ftr_results_AAHR->[$ftr_idx][$seq_idx]{x_score} to " . $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_score"} . "\n");
        }
        else { 
          $top_score_flag = 0;
        }
      }
      # Add elsif($key eq "") blocks here to store more values from the blastx.summary file
    }
  }
  close(IN);

  return 0;
}

#################################################################
# Subroutine:  check_for_defined_pstart_in_mdl_results()
# Incept:      EPN, Mon Oct 29 14:44:38 2018
#
# Purpose:    Given a feature and sequence index
#             in ftr_results.
#
# Arguments: 
#  $seq_idx:             sequence index
#  $ftr_idx:             feature index
#  $ftr_info_HAR:        REF to hash of arrays with information on the features
#  $mdl_results_AAHR:    REF to model results AAH
#  $FH_HR:               REF to hash of file handles
#
# Returns:    '1' if at least one model for this feature has a p_start value defined
#
# Dies:       If blastx fails.
#
################################################################# 
sub check_for_defined_pstart_in_mdl_results { 
  my $sub_name = "check_for_defined_pstart_in_mdl_results";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_idx, $ftr_idx, $ftr_info_HAR, $mdl_results_AAHR, $FH_HR) = @_;

  if($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { 
    for(my $mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx]; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
      if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]->{"p_start"}) { 
        return 1; 
      }
    }
  }
  elsif(($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") &&
        ($ftr_info_HAR->{"type"}[$ftr_idx]       eq "cds-mp")) { 
    my @all_children_idx_A = (); # feature indices of the all children of this feature
    getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "all", \@all_children_idx_A, $FH_HR);
    my $na_children = scalar(@all_children_idx_A);
    for(my $child_idx = 0; $child_idx < $na_children; $child_idx++) { 
      # call this subroutine recursively, but do a sanity check firs
      if($ftr_info_HAR->{"annot_type"}[$child_idx] ne "model") { 
        DNAORG_FAIL(sprintf("ERROR in $sub_name, about to make recursive call for feature idx $child_idx of annot_type %s != model", $ftr_info_HAR->{"annot_type"}[$child_idx]), 1, $FH_HR);
      }
      if(check_for_defined_pstart_in_mdl_results($seq_idx, $child_idx, $ftr_info_HAR, $mdl_results_AAHR, $FH_HR)) { 
        return 1; 
      }
    }
  }
  return 0; # if we get here, we have 0 p_start values in any models for this feature
}

#################################################################
# Subroutine:  helper_ftable_get_ftr_error_code_strings()
# Incept:      EPN, Tue Oct 30 11:56:39 2018
#
# Purpose:    Given a sequence name and feature index, construct
#             a string of all errors for that sequence/feature pair
#             and return it. Also add newline-terminated strings,
#             one per error, to @{$ftr_long_output_AR} for eventual
#             output to the long feature table.
#
# Arguments: 
#  $seq_name:               name of sequence
#  $ftr_idx:                feature index
#  $err_ftr_instances_AHHR: REF to array of 2D hashes with per-feature errors, PRE-FILLED
#  $err_info_HAR:           REF to the error info hash of arrays, PRE-FILLED
#  $ftr_long_output_AR:     REF to array of long output strings, FILLED HERE
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    A string with all err codes for this sequence/feature combo concatenated and 
#             separated by commas.
#
#
################################################################# 
sub helper_ftable_get_ftr_error_code_strings { 
  my $sub_name = "helper_ftable_get_ftr_error_code_strings";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $ftr_idx, $err_ftr_instances_AHHR, $err_info_HAR, $ftr_long_output_AR, $FH_HR) = @_;

  my $nerr = scalar(@{$err_info_HAR->{"code"}});

  my $ret_err_str = "";
  for(my $err_idx = 0; $err_idx < $nerr; $err_idx++) { 
    if($err_info_HAR->{"pertype"}[$err_idx] eq "feature") { 
      my $err_code = $err_info_HAR->{"code"}[$err_idx];
      if(exists $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name}) { 
        if($ret_err_str ne "") { $ret_err_str .= ","; }
        $ret_err_str .= $err_code;
        push(@{$ftr_long_output_AR}, sprintf("%3s error code: %s%s", 
                                             $err_code, 
                                             $err_info_HAR->{"desc"}[$err_idx], 
                                             ($err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} eq "") ? "" : " [" . $err_ftr_instances_AHHR->[$ftr_idx]{$err_code}{$seq_name} . "]")); 
      }
    }
  }
  
  return $ret_err_str;
}

#################################################################
# Subroutine:  helper_ftable_get_seq_error_code_strings()
# Incept:      EPN, Thu Jan 24 12:03:50 2019
#

# Purpose:    Given a sequence name, construct a string of all 
#             per-sequence errors and return it.
#
# Arguments: 
#  $seq_name:              name of sequence
#  $err_seq_instances_HHR: REF to 2D hashes with per-sequence errors, PRE-FILLED
#  $err_info_HAR:          REF to the error info hash of arrays, PRE-FILLED
#  $FH_HR:                 REF to hash of file handles
#
# Returns:    A string with all per-sequence err codes for this sequence 
#             concatenated and separated by commas.
#
################################################################# 
sub helper_ftable_get_seq_error_code_strings { 
  my $sub_name = "helper_ftable_get_seq_error_code_strings";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $err_seq_instances_HHR, $err_info_HAR, $FH_HR) = @_;

  my $nerr = scalar(@{$err_info_HAR->{"code"}});

  my $ret_err_str = "";
  for(my $err_idx = 0; $err_idx < $nerr; $err_idx++) { 
    if($err_info_HAR->{"pertype"}[$err_idx] eq "sequence") { 
      my $err_code = $err_info_HAR->{"code"}[$err_idx];
      if(exists $err_seq_instances_HHR->{$err_code}{$seq_name}) { 
        if($ret_err_str ne "") { $ret_err_str .= ","; }
        $ret_err_str .= $err_code;
      }
    }
  }
  
  return $ret_err_str;
}

#################################################################
# Subroutine:  helper_ftable_get_coords_standard
# Incept:      EPN, Tue Oct 30 12:59:13 2018
#
# Purpose:    Given a sequence name and feature index, construct
#             a feature table coordinate string, possibly of 
#             multiple lines, one per segment.
#
# Arguments: 
#  $seq_idx:           sequence index
#  $ftr_idx:           feature index
#  $do_start_carrot:   '1' to do carrot for first start
#  $do_stop_carrot:    '1' to do carrot for final stop
#  $do_pred_stop:      '1' to use predicted stop (p_stop) instead of corrected one
#  $ret_min_coord:     REF to minimum coordinate, to fill
#  $mdl_info_HAR:      REF to hash of arrays with information on the models, PRE-FILLED
#  $ftr_info_HAR:      REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:      REF to hash of arrays with information on the sequences, PRE-FILLED
#  $mdl_results_AAHR:  REF to model results AAH, PRE-FILLED
#  $ftr_results_AAHR:  REF to feature results AAH, PRE-FILLED
#  $FH_HR:             REF to hash of file handles
#
# Returns:    A string that gives the coordinates for the seq_idx/ftr_idx
#             pair in feature table format.
#
# Dies:       if x_start or x_stop does not exist in the ftr_results_AAHR->[$ftr_idx][$seq_idx] hash
################################################################# 
sub helper_ftable_get_coords_standard { 
  my $sub_name = "helper_ftable_get_coords_standard";
  my $nargs_exp = 12;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_idx, $ftr_idx, $do_start_carrot, $do_stop_carrot, $do_pred_stop, $ret_min_coord, $mdl_info_HAR, $ftr_info_HAR, $seq_info_HAR, $mdl_results_AAHR, $ftr_results_AAHR, $FH_HR) = @_;

  my $ftr_results_HR = \%{$ftr_results_AAHR->[$ftr_idx][$seq_idx]}; # for convenience

  my @start_A = ();
  my @stop_A  = ();

  if(($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "multifeature") &&
     ($ftr_info_HAR->{"type"}[$ftr_idx]       eq "cds-mp")) { 

    my $ftbl_out_start = $ftr_results_HR->{"ftbl_out_start"};
    my $ftbl_out_stop  = $ftr_results_HR->{"ftbl_out_stop"};

    if($do_pred_stop) { # need to overwrite $ftbl_out_stop
      my @cur_primary_children_idx_A = (); # feature indices of the primary children of this feature
      getPrimaryOrAllChildrenFromFeatureInfo($ftr_info_HAR, $ftr_idx, "primary", \@cur_primary_children_idx_A, $FH_HR);
      my $first_ftr_idx = $cur_primary_children_idx_A[0];
      my $final_ftr_idx = $cur_primary_children_idx_A[(scalar(@cur_primary_children_idx_A)-1)];
      my $first_child_mdl_idx = $ftr_info_HAR->{"first_mdl"}[$first_ftr_idx];
      my $final_child_mdl_idx = $ftr_info_HAR->{"final_mdl"}[$final_ftr_idx];
      # go through and determine final model between first_child_mdl_idx and final_child_mdl_idx that has a valid p_stop value, 
      # there must be at least 1, then set ftbl_out_stop to that (p_start is actually irrelevant here)
      my $pred_start = undef;
      my $pred_stop  = undef;
      for(my $m = $first_child_mdl_idx; $m <= $final_child_mdl_idx; $m++){ 
        if((defined $mdl_results_AAHR->[$m][$seq_idx]{"p_start"}) && 
           (defined $mdl_results_AAHR->[$m][$seq_idx]{"p_stop"})  && 
           ($mdl_results_AAHR->[$m][$seq_idx]{"p_start"} ne "?")  && 
           ($mdl_results_AAHR->[$m][$seq_idx]{"p_stop"} ne "?")) { 
          $pred_start = $mdl_results_AAHR->[$m][$seq_idx]{"p_start"};
          $pred_stop  = $mdl_results_AAHR->[$m][$seq_idx]{"p_stop"};
          if(($mdl_results_AAHR->[$m][$seq_idx]{"p_strand"} eq "+") && ($pred_stop > $seq_info_HAR->{"len"})) { $pred_stop = $seq_info_HAR->{"len"}; }
          if(($mdl_results_AAHR->[$m][$seq_idx]{"p_strand"} eq "-") && ($pred_stop < 1))                      { $pred_stop = 1; }
        }
      }
      if((! defined $pred_start) || (! defined $pred_stop)) { 
        DNAORG_FAIL("ERROR in $sub_name, feature type is multifeature, do_pred_stop is 1 but unable to find start/stop - shouldn't happen", 1, $FH_HR);
      }
      (undef, $ftbl_out_stop) = create_output_start_and_stop($pred_start, # irrelevant due to the first undef arg
                                                             $pred_stop, 
                                                             $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
    }

    if(($ftbl_out_start eq "?") || ($ftbl_out_stop eq "?")) { 
      return ""; # indicating we don't have a full CDS prediction for this seq/feature pair
    }
    push(@start_A, $ftbl_out_start);
    push(@stop_A,  $ftbl_out_stop);

    return helper_ftable_start_stop_arrays_to_coords(\@start_A, \@stop_A, $do_start_carrot, $do_stop_carrot, $ret_min_coord, $FH_HR);
  }
  elsif($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "model") { # not a multifeature cds-mp but annotated by models
    if(! check_for_defined_pstart_in_mdl_results($seq_idx, $ftr_idx, $ftr_info_HAR, $mdl_results_AAHR, $FH_HR)) { 
      DNAORG_FAIL("ERROR in $sub_name, feature annot_type is model but no models have any predicitions - shouldn't happen", 1, $FH_HR);
    }
    for(my $mdl_idx = $ftr_info_HAR->{"first_mdl"}[$ftr_idx]; $mdl_idx <= $ftr_info_HAR->{"final_mdl"}[$ftr_idx]; $mdl_idx++) { 
      if(exists $mdl_results_AAHR->[$mdl_idx][$seq_idx]->{"p_start"}) { 
        my $out_start = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"out_start"};
        my $out_stop  = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"out_stop"};
        if((! defined $out_start) || ($out_start eq "?")) { 
          DNAORG_FAIL("ERROR in $sub_name, feature annot_type is model but out_start is not defined or ? - shouldn't happen", 1, $FH_HR);
        }
        if((! defined $out_stop) || ($out_stop eq "?")) { 
          DNAORG_FAIL("ERROR in $sub_name, feature annot_type is model but out_stop is not defined or ? - shouldn't happen", 1, $FH_HR);
        }
        if($do_pred_stop) { 
          my $pred_start = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"};
          my $pred_stop  = $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_stop"} + $mdl_info_HAR->{"append_num"}[$mdl_idx];
          (undef, $out_stop) = create_output_start_and_stop($mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_start"}, # irrelevant due to the first undef arg
                                                            $mdl_results_AAHR->[$mdl_idx][$seq_idx]{"p_stop"}, 
                                                            $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);
        }
        push(@start_A, $out_start);
        push(@stop_A,  $out_stop);
      }
    }
    return helper_ftable_start_stop_arrays_to_coords(\@start_A, \@stop_A, $do_start_carrot, $do_stop_carrot, $ret_min_coord, $FH_HR);
  }
  else { 
    DNAORG_FAIL("ERROR in $sub_name, seq idx $seq_idx, ftr_idx $ftr_idx feature is not cds-mp or annot_type model - shouldn't happen", 1, $FH_HR);
  }

  return ""; # never reached
}
#################################################################
# Subroutine:  helper_ftable_get_coords_xnn_flag
# Incept:      EPN, Tue Oct 30 12:26:05 2018
#
# Purpose:    Given a sequence name and feature index, construct
#             a feature table coordinate string, possibly of 
#             multiple lines, one per segment, for the special
#             case that this feature has a 'xnn' error: blastx 
#             prediction but no CM prediction.
#
# Arguments: 
#  $seq_idx:           sequence index
#  $ftr_idx:           feature index
#  $do_start_carrot:   '1' to do carrot for first start
#  $do_stop_carrot:    '1' to do carrot for final stop
#  $ret_min_coord:     REF to minimum coordinate, to fill
#  $seq_info_HAR:      REF to hash of arrays with information on the sequences, PRE-FILLED
#  $ftr_results_AAHR:  REF to feature results AAH, PRE-FILLED
#  $FH_HR:             REF to hash of file handles
#
# Returns:    A string that gives the coordinates for the seq_idx/ftr_idx
#             pair in feature table format.
#
# Dies:       if x_start or x_stop does not exist in the ftr_results_AAHR->[$ftr_idx][$seq_idx] hash
################################################################# 
sub helper_ftable_get_coords_xnn_flag { 
  my $sub_name = "helper_ftable_get_coords_xnn_flag";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_idx, $ftr_idx, $do_start_carrot, $do_stop_carrot, $ret_min_coord, $seq_info_HAR, $ftr_results_AAHR, $FH_HR) = @_;

  # NOTE: for 'xnn' errors, the x_start and x_stop are always set at the feature level
  if((! exists $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_start"}) ||
     (! exists $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_stop"})) { 
    DNAORG_FAIL("ERROR in $sub_name, ftr_results_AAHR->[$ftr_idx][$seq_idx]{x_start|x_stop} does not exists", 1, $FH_HR);
  }

  my ($out_start, $out_stop) = create_output_start_and_stop($ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_start"},
                                                            $ftr_results_AAHR->[$ftr_idx][$seq_idx]{"x_stop"},
                                                            $seq_info_HAR->{"len"}[$seq_idx], $FH_HR);


  my @start_A = ($out_start);
  my @stop_A  = ($out_stop);

  return helper_ftable_start_stop_arrays_to_coords(\@start_A, \@stop_A, $do_start_carrot, $do_stop_carrot, $ret_min_coord, $FH_HR);
}

#################################################################
# Subroutine:  helper_ftable_start_stop_arrays_to_coords()
# Incept:      EPN, Tue Oct 30 12:39:59 2018
#
# Purpose:    Given refs to two arrays of start and stop coordinates,
#             construct coordinate strings in feature table format.
#
# Arguments: 
#  $start_AR:        REF to array of start coordinates
#  $stop_AR:         REF to array of stop coordinates
#  $do_start_carrot: '1' to do carrot for first start
#  $do_stop_carrot:  '1' to do carrot for final stop
#  $ret_min_coord:   REF to minimum coordinate, to fill
#  $FH_HR:           REF to hash of file handles
#
# Returns:    A string that gives the coordinates in feature table format.
#             Or "" if $start_AR->[0] and/or $stop_AR->[0] is "?" and size of those arrays is 1
#
# Dies: if either @{$start_AR} or @{$stop_AR} are empty
#       if @{$start_AR} and @{$stop_AR} are different sizes
#       if $start_AR->[$i] and/or $stop_AR->[$i] eq "?"
#
################################################################# 
sub helper_ftable_start_stop_arrays_to_coords { 
  my $sub_name = "helper_ftable_start_stop_arrays_to_coords";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($start_AR, $stop_AR, $do_start_carrot, $do_stop_carrot, $ret_min_coord, $FH_HR) = @_;

  my $ret_coords_str = "";
  my $ncoord = scalar(@{$start_AR});
  my $min_coord = undef;
  if($ncoord == 0) { 
    DNAORG_FAIL("ERROR in $sub_name, start_A array is empty", 1, $FH_HR);
  }
  if($ncoord != scalar(@{$stop_AR})) { # sanity check
    DNAORG_FAIL("ERROR in $sub_name, start_A array and stop_A arrays are different sizes", 1, $FH_HR);
  }
  for(my $c = 0; $c < $ncoord; $c++) { 
    my $is_first = ($c == 0)           ? 1 : 0;
    my $is_final = ($c == ($ncoord-1)) ? 1 : 0;
    my $start = $start_AR->[$c];
    my $stop  = $stop_AR->[$c];
    if($start eq "?" || $stop eq "?") { # should only happen for cds-mp features
      DNAORG_FAIL("ERROR in $sub_name, start ($start) and/or stop ($stop) are equal to ? - shouldn't happen", 1, $FH_HR);
    }
    else {
      if((! defined $min_coord) || ($start < $min_coord)) { $min_coord = $start; }
      if((! defined $min_coord) || ($stop  < $min_coord)) { $min_coord = $stop;  }
      $ret_coords_str .= sprintf("%s%d\t%s%d\n", 
                                 ($do_start_carrot && $is_first) ? "<" : "", $start, 
                                 ($do_stop_carrot  && $is_final) ? ">" : "", $stop);
    }
  }

  $$ret_min_coord = $min_coord;
  return $ret_coords_str;
}

#################################################################
# Subroutine:  helper_ftable_coords_to_out_str()
# Incept:      EPN, Tue Oct 30 13:37:51 2018
#
# Purpose:    Convert a string with >= 1 feature table lines of
#             coordinates to an output string for a feature.
#
# Arguments: 
#  $coords_str:      REF to array of start coordinates
#  $feature_type:    name of feature type, e.g. "CDS"
#  $FH_HR:           REF to hash of file handles
#
# Returns:    The output string, which is the coords string with the
#             feature type added at the appropriate place.
#
# Dies: if $coords_str is empty
#
################################################################# 
sub helper_ftable_coords_to_out_str { 
  my $sub_name = "helper_ftable_coords_to_out_str";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($coords_str, $feature_type, $FH_HR) = @_;

  my $ret_out_str = "";

  if($coords_str eq "") { 
    DNAORG_FAIL("ERROR in $sub_name, coords_str is empty - shouldn't happen", 1, $FH_HR);
  }
  my @coords_A = split("\n", $coords_str);

  for(my $cidx = 0; $cidx < scalar(@coords_A); $cidx++) { 
    $ret_out_str .= $coords_A[$cidx];
    if($cidx == 0) { $ret_out_str .= "\t" . $feature_type; }
    $ret_out_str .= "\n";
  }
 
  return $ret_out_str;
}

#################################################################
# Subroutine:  helper_ftable_add_qualifier_from_ftr_info()
# Incept:      EPN, Tue Oct 30 13:41:58 2018
#
# Purpose:    Add a qualifier line to a string that will be 
#             part of a feature table output.
#
# Arguments: 
#  $ftr_idx:      feature index
#  $key:          key in ftr_info_HAR
#  $qval_sep:     characters that separate multiple qualifiers in $ftr_info_HAR->{$key}[$ftr_idx]
#  $ftr_info_HAR: REF to hash of arrays with information on the features, PRE-FILLED
#  $FH_HR:        REF to hash of file handles
#
# Returns:    Strings to append to the feature table.
#
# Dies: never
#
################################################################# 
sub helper_ftable_add_qualifier_from_ftr_info {
  my $sub_name = "helper_ftable_add_qualifier_from_ftr_info";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_idx, $key, $qval_sep, $ftr_info_HAR, $FH_HR) = @_;
  
  my $ret_str = "";
  if((exists $ftr_info_HAR->{$key}[$ftr_idx]) && ($ftr_info_HAR->{$key}[$ftr_idx] ne "")) { 
    my $qualifier_name = featureInfoKeyToFeatureTableQualifierName($key, $FH_HR);
    my @qval_A = split($qval_sep, $ftr_info_HAR->{$key}[$ftr_idx]); 
    foreach my $qval (@qval_A) { 
      $ret_str .= sprintf("\t\t\t%s\t%s\n", $qualifier_name, $qval);
    }
  }
  return $ret_str;
}

#################################################################
# Subroutine:  helper_ftable_add_qualifier_from_ftr_results()
# Incept:      EPN, Tue Oct 30 13:52:19 2018
#
# Purpose:    Add a qualifier line to a string that will be 
#             part of a feature table output.
#
# Arguments: 
#  $seq_idx:           sequence index
#  $ftr_idx:           feature index
#  $results_key:       key in ftr_results_AAHR
#  $qualifier:         name for the qualifier
#  $ftr_results_AAHR:  REF to feature results AAH, PRE-FILLED
#  $FH_HR:             REF to hash of file handles
#
# Returns:    "" if $ftr_results_AAHR->[$ftr_idx][$seq_idx]{$results_key} does not exist
#             else a string for the feature table
#
# Dies: never
#
################################################################# 
sub helper_ftable_add_qualifier_from_ftr_results {
  my $sub_name = "helper_ftable_add_qualifier_from_ftr_results";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_idx, $ftr_idx, $results_key, $qualifier, $ftr_results_AAHR, $FH_HR) = @_;

  my $ret_str = "";

  if(defined $ftr_results_AAHR->[$ftr_idx][$seq_idx]{$results_key}) { 
    $ret_str = sprintf("\t\t\t%s\t%s\n", $qualifier, $ftr_results_AAHR->[$ftr_idx][$seq_idx]{$results_key});
  }
  return $ret_str;
}

#################################################################
# Subroutine:  parse_class_errors_list_file
# Incept:      EPN, Wed Dec 12 13:44:21 2018
#
# Purpose:    Parse the --classerrors input file
#
# Arguments: 
#  $in_file:       file to parse with per-sequence classification errors
#  $errors_seq_HR: ref to hash of classification errors per sequence, filled here
#                  with unmodified lines from $in_file
#  $FH_HR:         ref to hash of file handles
#
# Returns:    "" if $ftr_results_AAHR->[$ftr_idx][$seq_idx]{$results_key} does not exist
#             else a string for the feature table
#
# Dies: never
#
################################################################# 
sub parse_class_errors_list_file { 
  my $sub_name = "parse_class_errors_list_file";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($in_file, $errors_seq_HR, $FH_HR) = @_;

  open(IN, $in_file) || fileOpenFailure($in_file, $sub_name, $!, "reading", $FH_HR);
  my $line;
  while($line = <IN>) { 
    if($line !~ m/^\#/) { 
      chomp $line;
      my @el_A = split(/\t/, $line);
      ##sequence	error	feature	error-description
      #MG763368.1	Unexpected Classification	*sequence*	NC 001959,NC 029647 was specified, but NC 039476 is predicted
      if(scalar(@el_A) != 4) { 
        foreach my $tok (@el_A) { printf("tok: $tok\n"); }
        DNAORG_FAIL("ERROR in $sub_name, did not find exactly 2 tokens in line $line", 1, $FH_HR);
      }
      my $seq = $el_A[0];
      $errors_seq_HR->{$seq} .= $line . "\n";
    }
  }
  close(IN);

  return;
}


#################################################################
# Subroutine: error_list_output_to_ftable_errors()
# Incept:     EPN, Wed Dec 12 14:02:24 2018
#
# Purpose:    Given output from a error list file, with 4 tab-delimited
#             tokens per new-line delimited string, return a string
#             that can be output to a feature table with the 
#             same information.
#
#             The input will be new-line delimited strings, each of 
#             which will have 4 tab-delimited tokens:
#             <sequence-name>
#             <error-name>
#             <feature-name>
#             <error-description>
#
# Arguments:
#   $in_seqname: expected name of sequence
#   $errliststr: error string to convert
#   $FH_HR:      ref to hash of file handles, including 'log'
#             
# Returns:    $err_ftbl_str: feature table strings in the format described above.
#
# Dies: If $errliststr is not in required format, or if it includes
#       data for a sequence other than <$seqname>
#
#################################################################
sub error_list_output_to_ftable_errors { 
  my $sub_name  = "error_list_output_to_ftable_errors";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_seqname, $errliststr, $FH_HR) = (@_);

  my @errliststr_A = split(/\n/, $errliststr);

  my $retstr = "";
  foreach my $errlistline (@errliststr_A) { 
    my @el_A = split(/\t/, $errlistline); 
    if(scalar(@el_A) != 4) { 
      DNAORG_FAIL("ERROR in $sub_name, did not read exactly 4 tab-delimited tokens in line: $errlistline", 1, $FH_HR);
    }
    my ($seqname, $error_name, $feature_name, $error_desc) = (@el_A);
    if($in_seqname ne $seqname) { 
      DNAORG_FAIL("ERROR in $sub_name, read unexpected sequence name $seqname != $in_seqname in line: $errlistline", 1, $FH_HR);
    }
    $retstr .= sprintf("ERROR: $error_name: %s$error_desc\n", ($feature_name eq "*sequence*") ? "" : "($feature_name) ");
  }

  return $retstr;
}

#################################################################
# Subroutine: helper_blastx_breakdown_query()
# Incept:     EPN, Wed Dec 19 12:05:02 2018
#
# Purpose:    Given a query name from blastx output, determine if
#             it is for an entire input sequence, or a feature
#             sequence fetched from an input sequence. 
#             The way we tell is by looking up $in_query in $seq_info_HAR.
#             If it exists as a key, then the query is an entire input
#             sequence. If it does not exist as a key, then it should
#             be <seqname>/<coords_str>, and <seqname> should be 
#             a key in $seq_info_HAR.
#
# Arguments:
#   $in_query:     query name from blastx output 
#   $exp_seqname:  expected sequence name, of undef if unknown
#                  can only be undef if $seq_HR is undefined
#   $seq_HR:       ref to hash, keys are possible sequence names, values are irrelevant
#   $FH_HR:        ref to hash of file handles, including 'log'
#             
# Returns:  Three values:
#           <seqname>:    name of input sequence this query corresponds to
#           <coords_str>: undef if $in_query == <seqname>, else the coords
#                         string for the fetched feature from <seqname>
#           <len>:        undef if $in_query == <seqname>, total length of 
#                         coordinate ranges listed in <coords_str>
#
# Dies: If $in_query is not a valid key in $seq_info_HAR, and we can't
#       break it down into a valid <seqname> and <coords_str>.
#
#################################################################
sub helper_blastx_breakdown_query {
  my $sub_name  = "helper_blastx_breakdown_query";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_query, $exp_seqname, $seq_HR, $FH_HR) = (@_);

  my $ret_len = 0;
  my $ret_seqname = undef;
  my $ret_coords_str = undef;

  # contract check
  if((! defined $exp_seqname) && (! defined $seq_HR)) { 
    DNAORG_FAIL("ERROR in $sub_name, contract failure: both <exp_seqname> and <seq_HR> are undefined", 1, $FH_HR);
  }
  if((defined $exp_seqname) && (defined $seq_HR)) { 
    DNAORG_FAIL("ERROR in $sub_name, contract failure: both <exp_seqname> and <seq_HR> are defined", 1, $FH_HR);
  }

  if((defined $exp_seqname) && ($in_query eq $exp_seqname)) {   
    return($in_query, undef, undef);
  }
  elsif((! defined $exp_seqname) && (exists $seq_HR->{$in_query})) { 
    return($in_query, undef, undef);
  }
  else { 
    # sequence name does not exist, as is
    # example: NC_002549.1/6039-8068/1-885,885-2030
    # sequence is NC_002549.1/6039-8068 but we fetched the feature 1-885,885-2030 from it
    if($in_query =~ /^(\S+)\/([\d\,\-]+$)/) { 
      ($ret_seqname, $ret_coords_str) = ($1, $2);
      if((defined $exp_seqname) && ($ret_seqname ne $exp_seqname)) { 
        DNAORG_FAIL("ERROR in $sub_name, unexpected sequence name $ret_seqname != $exp_seqname, derived from $in_query", 1, $FH_HR); 
      }
      if((! defined $exp_seqname) && (! exists $seq_HR->{$ret_seqname})) { 
        DNAORG_FAIL("ERROR in $sub_name, unable to find input sequence with name $in_query or $ret_seqname", 1, $FH_HR); 
      }

      # if we get here $ret_seqname is a valid sequence
      my @range_A = split(",", $ret_coords_str);
      foreach my $range (@range_A) { 
        if($range =~ /^(\d+)\-(\d+)$/) { 
          $ret_len += abs($1 - $2) + 1;
        }
        else { 
          DNAORG_FAIL("ERROR in $sub_name, unable to parse coords string $ret_coords_str (element $range) in blastx query sequence $in_query", 1, $FH_HR); 
        }
      }
    }
    else { 
      DNAORG_FAIL("ERROR in $sub_name, unable to parse blastx query sequence name $in_query", 1, $FH_HR); 
    }
  }
  return ($ret_seqname, $ret_coords_str, $ret_len);
}

#################################################################
# Subroutine: add_zft_errors()
# Incept:     EPN, Thu Jan 24 12:31:16 2019
# Purpose:    Adds zft errors for sequences with 0 predicted features
#
# Arguments:
#  $err_ftr_instances_AHHR:  REF to array of 2D hashes with per-feature errors, PRE-FILLED
#  $err_seq_instances_HHR:   REF to 2D hash with per-sequence errors, PRE-FILLED
#  $ftr_info_HAR:            REF to hash of arrays with information on the features, PRE-FILLED
#  $seq_info_HAR:            REF to hash of arrays with information on the sequences, PRE-FILLED
#  $err_info_HAR:            REF to the error info hash of arrays, PRE-FILLED
#  $mdl_results_AAHR:        REF to model results AAH, PRE-FILLED
#  $opt_HHR:                 REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:          REF to the 2D hash of output file information
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub add_zft_errors { 
  my $sub_name = "add_zft_errors";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($err_ftr_instances_AHHR, $err_seq_instances_HHR, $ftr_info_HAR, $seq_info_HAR, 
      $err_info_HAR, $mdl_results_AAHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience
  my $nftr = validateFeatureInfoHashIsComplete  ($ftr_info_HAR, undef, $FH_HR); # nftr: number of features
  my $nseq = validateSequenceInfoHashIsComplete ($seq_info_HAR, undef, $opt_HHR, $FH_HR); # nseq: number of sequences

  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name  = $seq_info_HAR->{"seq_name"}[$seq_idx];
    my $seq_nftr = 0; # number of annotated features for this sequence

    # loop over features
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      # determine if this feature can be ignored or cannot be ignored (should be annotated), 
      # it can be ignored if:
      # - we do not have a defined "p_start" for any model associated with this feature
      #   (this is equivalent to having a 'nop' error for all models associated with this feature)
      # - we do not have a 'xnn' error for this feature
      my $is_duplicate  = ($ftr_info_HAR->{"annot_type"}[$ftr_idx] eq "duplicate") ? 1 : 0;
      my $defined_pstart = 0;
      my $xnn_flag       = 0;
      my $do_ignore      = 1; 
      if(! $is_duplicate) { 
        $defined_pstart = check_for_defined_pstart_in_mdl_results($seq_idx, $ftr_idx, $ftr_info_HAR, $mdl_results_AAHR, $FH_HR);
        $xnn_flag       = (exists $err_ftr_instances_AHHR->[$ftr_idx]{"xnn"}{$seq_name}) ? 1 : 0;
        $do_ignore      = ($defined_pstart || $xnn_flag) ? 0 : 1;
      }
      if(! $do_ignore) { 
        $seq_nftr++;
        $ftr_idx = $nftr; # breaks for $ftr_idx loop
      } 
    }
    if($seq_nftr == 0) { 
      error_instances_add(undef, $err_seq_instances_HHR, $err_info_HAR, -1, "zft", $seq_name, "-", $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine: add_dmo_errors()
# Incept:     EPN, Thu Feb  7 11:54:56 2019
# Purpose:    Adds dmo errors for sequences listed in the array @overflow_seq_A, if any.
#
# Arguments:
#  $overflow_seq_AR:         REF to array of sequences that failed due to matrix overflows, pre-filled
#  $overflow_mxsize_AR:      REF to array of required matrix sizes for each sequence that failed due to matrix overflows, pre-filled
#  $err_seq_instances_HHR:   REF to 2D hash with per-sequence errors, PRE-FILLED
#  $err_info_HAR:            REF to the error info hash of arrays, PRE-FILLED
#  $opt_HHR:                 REF to 2D hash of option values, see top of epn-options.pm for description
#  $ofile_info_HHR:          REF to the 2D hash of output file information
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub add_dmo_errors { 
  my $sub_name = "add_dmo_errors";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($overflow_seq_AR, $overflow_mxsize_AR, $err_seq_instances_HHR, $err_info_HAR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience

  my $noverflow = scalar(@{$overflow_seq_AR});
  for(my $s = 0; $s < $noverflow; $s++) { 
    error_instances_add(undef, $err_seq_instances_HHR, $err_info_HAR, -1, "dmo", $overflow_seq_AR->[$s], "required matrix size: $overflow_mxsize_AR->[$s] Mb", $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: convert_pp_char_to_pp_avg()
# Incept:     EPN, Thu Feb  7 11:54:56 2019
# Purpose:    Convert a cmalign alignment PP value to the average posterior 
#             probability it represents.
#
# Arguments:
#  $ppchar:  the PP character from the alignment, cannot be a gap
#  $FH_HR:   ref to hash of file handles, including 'log'
#             
# Returns:  average value for $ppchar
# 
# Dies:     if $ppchar is not ['0'-'9'] or '*'
#
#################################################################
sub convert_pp_char_to_pp_avg { 
  my $sub_name = "convert_pp_char_to_pp_avg";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ppchar, $FH_HR) = @_;

  if($ppchar eq "*") { return 0.975; }
  if($ppchar eq "9") { return 0.90; }
  if($ppchar eq "8") { return 0.80; }
  if($ppchar eq "7") { return 0.70; }
  if($ppchar eq "6") { return 0.60; }
  if($ppchar eq "5") { return 0.50; }
  if($ppchar eq "4") { return 0.40; }
  if($ppchar eq "3") { return 0.30; }
  if($ppchar eq "2") { return 0.20; }
  if($ppchar eq "1") { return 0.10; }
  if($ppchar eq "0") { return 0.25; }

  DNAORG_FAIL("ERROR in $sub_name, invalid PP char: $ppchar", 1, $FH_HR); 

  return 0.; # NEVER REACHED
}
