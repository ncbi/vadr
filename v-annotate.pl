#!/usr/bin/env perl
# EPN, Wed May  1 10:18:55 2019 [renamed to v-annotate.pl]
# EPN, Thu Feb 18 12:48:16 2016 [dnaorg_annotate.pl split off from dnaorg_annotate_genomes.pl]
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
#
use strict;
use warnings;
use Getopt::Long qw(:config no_auto_abbrev);
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

require "vadr.pm"; 
require "vadr_seed.pm"; 
require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seq.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

#######################################################################################
# What this script does: 
#
# Given an input sequence file, each sequence is compared against a
# library of homology models (CMs) and classified and annotated and
# determined to PASS or FAIL.  Certain types of unexpected results are
# detected and reported in the output as 'alerts'. There are two
# classes of alerts: those that cause a sequence to FAIL and those
# that do not. 
#
# The script proceeds through four main stages:
#
# (1) classification: each sequence is compared against the CM library
#     using a relatively fast HMM scoring algorithm, and the highest
#     scoring model in the library is defined as the winning CM for
#     that sequence.
#
# (2) coverage determination: each sequence is compared against its
#     winning model (only) for a second time using a more expensive HMM
#     scoring algorithm that is local with respect to the model and
#     sequence. This stage allows statistics related to the coverage of
#     the sequence and model to be determined, and some alerts can be
#     reported based on those statisics.
#
# (3) alignment/annotation: each sequence is aligned to its winning
#     model using a still more expensive CM algorithm that takes into
#     account secondary structure in the model (if any). This algorithm
#     is aligns the full sequence either locally or globally with
#     respect to the model. Features are then annotated based on the
#     alignment coordinates and the known feature coordinates in the 
#     model (supplied via the modelinfo file). 
#   
# (4) protein validation: CDS features are then validated via
#     blastx or hmmer by comparing predicted feature spans from (3) to
#     pre-computed BLAST or HMMER databases for the model. Alerts can
#     be reported based on the blast/hmmer results. 
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Important options that change this behavior:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 
# Seeded alignment for faster processing (originally developed for
# SARS-CoV-2 sequences), enabled with the -s option:
#
# Stage 1 and 2 are performed by a single blastn search instead of two
# rounds of cmsearch.  The top scoring HSP is identified and used
# (after potentially trimming) to 'seed' the alignment of that
# sequence by cmalign. This blastn alignmetn seed region is considered
# fixed and only the sequence before and after it (plus 100nt of
# overlap on each side, controllable with --s_overhang) is aligned
# separately by cmalign. The up to three alignments are then joined to
# get the final alignment.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 
# Replacement of stretches of Ns with the -r option:
# 
# Stage 1 and 2 are run as pre-processing steps to identify stretches
# of each sequence that are N-rich and replace the Ns in them with the
# expected nucleotide from the best-matching model, where possible.
# After Ns are replaced, the new sequence with Ns replaced is analyzed
# by all four stages a new (so stages 1 and 2 are performed twice, once
# on the input sequences and once on those same sequences but with Ns
# replaced. 
#
# By default, stages 1 and 2 are performed with blastn which in
# anecdotal (but not systematic) testing seems less likely then cmsearch
# to extend alignments through stretches of Ns (although this is
# probably controllable to an extent with command line options).
# However, with --r_prof, preprocessing stages 1 and 2 are performed
# with cmsearch.
#
#######################################################################################
#
# Alert identification:
#
# This script identifies and outputs a list of all alerts in each of
# the sequences it annotates. Each type of alert has an associated
# five letter alert 'code'. The list of alert codes is in the
# vadr.pm perl module file in the subroutine: vdr_AlertInfoInitialize().
#
# A table of alerts is output by v-annotate.pl when the --alt_list
# option is used.
# 
# List of subroutines in which alerts are detected and added:
#  1. alert_add_ambgnt5s_ambgnt3s()
#     ambgnt5s, ambgnt3s (2)
#
#  2. add_classification_alerts()
#     noannotn, lowscore, indfclas, qstsbgrp, qstgroup, incsbgrp, incgroup, revcompl, lowcovrg, biasdseq (10)
#
#  3. alert_add_unexdivg()
#     unexdivg (1)
#
#  4. parse_stk_and_add_alignment_cds_and_mp_alerts()
#     indf5gap, indf5lcc, indf5lcn, indf3gap, indf3lcc, indf3lcn, deletinf, deletins, deletina (9)
#
#  5. fetch_features_and_add_cds_and_mp_alerts_for_one_sequence()
#     mutstart, unexleng, mutendcd, mutendex, mutendns, cdsstopn, ambgnt5c, ambgnt3c, ambgnt5f, ambgnt3f (10)
#
#  6. add_protein_validation_alerts()
#     indfantn, indfstrp, indf5plg, indf5pst, indf3plg, indf3pst, insertnp, deletinp, cdsstopp, indfantp (10)
#
#  7. alert_add_parent_based()
#     peptrans (1)
# 
#  8. add_low_similarity_alerts_for_one_sequence()
#     lowsim5c, lowsim3c, lowsimic, lowsim5n, lowsim3n, lowsimin, lowsim5s, lowsim3s, lowsimis (9)
# 
#  9. add_frameshift_alerts_for_one_sequence()
#     fsthicft, fsthicfi, fstloft, fstlocfi, fstukcft, fstukcfi (6)
#
# 10. join_alignments_and_add_unjoinbl_alerts()
#     unjoinbl (1)
#
# 12. output_feature_table()
#     noftrann, noftrant, ftskipfl (1)
#
#######################################################################################
# make sure required environment variables are set
my $env_vadr_scripts_dir  = utl_DirEnvVarValid("VADRSCRIPTSDIR");
my $env_vadr_model_dir    = utl_DirEnvVarValid("VADRMODELDIR");
my $env_vadr_blast_dir    = utl_DirEnvVarValid("VADRBLASTDIR");
my $env_vadr_infernal_dir = utl_DirEnvVarValid("VADRINFERNALDIR");
my $env_vadr_hmmer_dir    = utl_DirEnvVarValid("VADRHMMERDIR");
my $env_vadr_easel_dir    = utl_DirEnvVarValid("VADREASELDIR");
my $env_vadr_bioeasel_dir = utl_DirEnvVarValid("VADRBIOEASELDIR");
my $env_vadr_fasta_dir    = utl_DirEnvVarValid("VADRFASTADIR");

my %execs_H = (); # hash with paths to all required executables
$execs_H{"cmalign"}       = $env_vadr_infernal_dir . "/cmalign";
$execs_H{"cmemit"}        = $env_vadr_infernal_dir . "/cmemit";
$execs_H{"cmfetch"}       = $env_vadr_infernal_dir . "/cmfetch";
$execs_H{"cmsearch"}      = $env_vadr_infernal_dir . "/cmsearch";
$execs_H{"hmmfetch"}      = $env_vadr_hmmer_dir    . "/hmmfetch";
$execs_H{"hmmscan"}       = $env_vadr_hmmer_dir    . "/hmmscan";
$execs_H{"hmmsearch"}     = $env_vadr_hmmer_dir    . "/hmmsearch";
$execs_H{"esl-alimerge"}  = $env_vadr_easel_dir    . "/esl-alimerge";
$execs_H{"esl-alimanip"}  = $env_vadr_easel_dir    . "/esl-alimanip";
$execs_H{"esl-reformat"}  = $env_vadr_easel_dir    . "/esl-reformat";
$execs_H{"esl-seqstat"}   = $env_vadr_easel_dir    . "/esl-seqstat";
$execs_H{"esl-translate"} = $env_vadr_easel_dir    . "/esl-translate";
$execs_H{"esl-ssplit"}    = $env_vadr_bioeasel_dir . "/scripts/esl-ssplit.pl";
$execs_H{"blastx"}        = $env_vadr_blast_dir    . "/blastx";
$execs_H{"blastn"}        = $env_vadr_blast_dir    . "/blastn";
$execs_H{"makeblastdb"}   = $env_vadr_blast_dir    . "/makeblastdb";
$execs_H{"parse_blast"}   = $env_vadr_scripts_dir  . "/parse_blast.pl";
$execs_H{"glsearch"}      = $env_vadr_fasta_dir    . "/glsearch36";
utl_ExecHValidate(\%execs_H, undef);

#########################################################
# Command line and option processing using sqp_opts.pm
#
# opt_HH: 2D hash:
#         1D key: option name (e.g. "-h")
#         2D key: string denoting type of information 
#                 (one of "type", "default", "group", "requires", "incompatible", "preamble", "help")
#         value:  string explaining 2D key:
#                 "type":         "boolean", "string", "integer" or "real"
#                 "default":      default value for option
#                 "group":        integer denoting group number this option belongs to
#                 "requires":     string of 0 or more other options this option requires to work, each separated by a ','
#                 "incompatible": string of 0 or more other options this option is incompatible with, each separated by a ','
#                 "preamble":     string describing option for preamble section (beginning of output from script)
#                 "help":         string describing option for help section (printed if -h used)
#                 "setby":        '1' if option set by user, else 'undef'
#                 "value":        value for option, can be undef if default is undef
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
#     option            type       default group   requires incompat    preamble-output                                   help-output    
opt_Add("-h",           "boolean", 0,          0,    undef, undef,      undef,                                            "display this help",                                  \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "basic options";
#     option            type       default group   requires incompat    preamble-output                                                                            help-output    
opt_Add("-f",           "boolean", 0,         $g,    undef, undef,      "force directory overwrite",                                                               "force; if output dir exists, overwrite it",   \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,         $g,    undef, undef,      "be verbose",                                                                              "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--atgonly",    "boolean", 0,         $g,    undef, undef,      "only consider ATG a valid start codon",                                                   "only consider ATG a valid start codon", \%opt_HH, \@opt_order_A);
opt_Add("--minpvlen",   "integer", 30,        $g,    undef, undef,      "min CDS/mat_peptide/gene length for feature table output and protein validation is <n>",  "min CDS/mat_peptide/gene length for feature table output and protein validation is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--nkb",        "integer", 300,       $g,    undef,  undef,     "number of KB of sequence for each alignment job and/or chunk is <n>",                     "number of KB of sequence for each alignment job and/or chunk is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,         $g,    undef, undef,      "leaving intermediate files on disk",                                                      "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for specifying classification";
#        option               type   default  group  requires incompat    preamble-output                                                     help-output    
opt_Add("--group",         "string",  undef,     $g,     undef, undef,     "set expected classification of all seqs to group <s>",             "set expected classification of all seqs to group <s>",            \%opt_HH, \@opt_order_A);
opt_Add("--subgroup",      "string",  undef,     $g, "--group", undef,     "set expected classification of all seqs to subgroup <s>",          "set expected classification of all seqs to subgroup <s>",         \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling severity of alerts";
#        option               type   default  group  requires incompat         preamble-output                                                                               help-output    
opt_Add("--alt_list",     "boolean",  0,         $g,     undef, undef,         "output summary of all alerts and exit",                                                      "output summary of all alerts and exit",                                \%opt_HH, \@opt_order_A);
opt_Add("--alt_pass",      "string",  undef,     $g,     undef, undef,         "specify that alert codes in <s> do not cause FAILure",                                       "specify that alert codes in comma-separated <s> do not cause FAILure", \%opt_HH, \@opt_order_A);
opt_Add("--alt_fail",      "string",  undef,     $g,     undef, undef,         "specify that alert codes in <s> cause FAILure",                                              "specify that alert codes in comma-separated <s> do cause FAILure",     \%opt_HH, \@opt_order_A);
opt_Add("--alt_mnf_yes",   "string",  undef,     $g,     undef,"--ignore_mnf", "alert codes in <s> for 'misc_not_failure' features cause misc_feature-ization, not failure", "alert codes in <s> for 'misc_not_failure' features cause misc_feature-ization, not failure", \%opt_HH, \@opt_order_A);
opt_Add("--alt_mnf_no",    "string",  undef,     $g,     undef,"--ignore_mnf", "alert codes in <s> for 'misc_not_failure' features cause failure, not misc_feature-ization", "alert codes in <s> for 'misc_not_failure' features cause failure, not misc-feature-ization", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for ignoring specific keys in the input model info (.minfo) file";
#        option               type        default  group requires incompat  preamble-output                                                                               help-output    
opt_Add("--ignore_mnf",       "boolean",  0,       $g,     undef, undef,    "ignore non-zero 'misc_not_failure' values in .minfo file, set to 0 for all features/models", "ignore non-zero 'misc_not_feature' values in .minfo file, set to 0 for all features/models", \%opt_HH, \@opt_order_A);
opt_Add("--ignore_isdel",     "boolean",  0,       $g,     undef, undef,    "ignore non-zero 'is_deletable' values in .minfo file, set to 0 for all features/models",     "ignore non-zero 'is_deletable' values in .minfo file, set to 0 for all features/models", \%opt_HH, \@opt_order_A);
opt_Add("--ignore_afset",     "boolean",  0,       $g,     undef, undef,    "ignore 'alternative_ftr_set' and 'alternative_ftr_set_subn' values in .minfo file",          "ignore 'alternative_ftr_set' and 'alternative_ftr_set_subn' values in .minfo file", \%opt_HH, \@opt_order_A);
opt_Add("--ignore_afsetsubn", "boolean",  0,       $g,     undef, undef,    "ignore 'alternative_ftr_set_subn' values in .minfo file",                                    "ignore 'alternative_ftr_set_subn' values in .minfo file", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options related to model files";
#        option               type default  group  requires incompat   preamble-output                                                                   help-output    
opt_Add("-m",           "string",  undef,      $g,    undef, undef,       "use CM file <s> instead of default",                                             "use CM file <s> instead of default", \%opt_HH, \@opt_order_A);
opt_Add("-a",           "string",  undef,      $g, "--pv_hmmer",undef,    "use protein HMM file <s> instead of default",                                    "use protein HMM file <s> instead of default", \%opt_HH, \@opt_order_A);
opt_Add("-i",           "string",  undef,      $g,    undef, undef,       "use model info file <s> instead of default",                                     "use model info file <s> instead of default", \%opt_HH, \@opt_order_A);
opt_Add("-n",           "string",  undef,      $g,    undef, undef,       "use blastn db file <s> instead of default",                                      "use blastn db file <s> instead of default",  \%opt_HH, \@opt_order_A);
opt_Add("-x",           "string",  undef,      $g,    undef, undef,       "blastx dbs are in dir <s>, instead of default",                                  "blastx dbs are in dir <s>, instead of default", \%opt_HH, \@opt_order_A);
opt_Add("--mkey",       "string","calici",     $g,    undef,"-m,-i,-a",   ".cm, .minfo, blastn .fa files in \$VADRMODELDIR start with key <s>, not 'vadr'", ".cm, .minfo, blastn .fa files in \$VADRMODELDIR start with key <s>, not 'vadr'",  \%opt_HH, \@opt_order_A);
opt_Add("--mdir",       "string",  undef,      $g,    undef, undef,       "model files are in directory <s>, not in \$VADRMODELDIR",                        "model files are in directory <s>, not in \$VADRMODELDIR",  \%opt_HH, \@opt_order_A);
opt_Add("--mlist",      "string",  undef,      $g,    undef, "-s",        "only use models listed in file <s>",                                             "only use models listed in file <s>",  \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling output feature table";
#        option               type   default group  requires incompat    preamble-output                                                               help-output    
opt_Add("--nomisc",       "boolean",  0,        $g,    undef,   undef,      "in feature table for failed seqs, never change feature type to misc_feature",                 "in feature table for failed seqs, never change feature type to misc_feature", \%opt_HH, \@opt_order_A);
opt_Add("--notrim",       "boolean",  0,        $g,    undef,   undef,      "in feature table, don't trim coords due to ambiguities (for any feature types)",              "in feature table, don't trim coords due to ambiguities (for any feature types)",   \%opt_HH, \@opt_order_A);
opt_Add("--noftrtrim",    "string",   undef,    $g,    undef,"--notrim",    "in feature table, don't trim coords due to ambiguities for ftr types in comma-delimited <s>", "in feature table, don't trim coords due to ambiguities for feature types in comma-delmited <s>",  \%opt_HH, \@opt_order_A);
opt_Add("--noprotid",     "boolean",  0,        $g,    undef,   undef,      "in feature table, don't add protein_id for CDS and mat_peptides",                             "in feature table, don't add protein_id for CDS and mat_peptides",         \%opt_HH, \@opt_order_A);
opt_Add("--forceprotid",  "boolean",  0,        $g,    undef,"--noprotid",  "in feature table, force protein_id value to be sequence name, then idx",                      "in feature table, force protein_id value to be sequence name, then idx",  \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling thresholds related to alerts";
#       option          type         default  group   requires incompat           preamble-output                                                                    help-output    
opt_Add("--lowsc",      "real",      0.3,       $g,   undef,   undef,            "lowscore/LOW_SCORE bits per nucleotide threshold is <x>",                         "lowscore/LOW_SCORE bits per nucleotide threshold is <x>",                         \%opt_HH, \@opt_order_A);
opt_Add("--indefclass", "real",      0.03,      $g,   undef,   undef,            "indfclas/INDEFINITE_CLASSIFICATION bits per nucleotide diff threshold is <x>",    "indfcls/INDEFINITE_CLASSIFICATION bits per nucleotide diff threshold is <x>",     \%opt_HH, \@opt_order_A);
opt_Add("--incspec",    "real",      0.2,       $g,   undef,   undef,            "inc{group,subgrp}/INCORRECT_{GROUP,SUBGROUP}' bits/nt threshold is <x>",          "inc{group,subgrp}/INCORRECT_{GROUP,SUBGROUP} bits/nt threshold is <x>",           \%opt_HH, \@opt_order_A);
opt_Add("--lowcov",     "real",      0.9,       $g,   undef,   undef,            "lowcovrg/LOW_COVERAGE fractional coverage threshold is <x>",                      "lowcovrg/LOW_COVERAGE fractional coverage threshold is <x>",                      \%opt_HH, \@opt_order_A);
opt_Add("--dupregolp",  "integer",   20,        $g,   undef,   undef,            "dupregin/DUPLICATE_REGIONS minimum model overlap is <n>",                         "dupregin/DUPLICATE_REGIONS minimum model overlap is <n>",                         \%opt_HH, \@opt_order_A);
opt_Add("--dupregsc",   "real",      10,        $g,   undef,   undef,            "dupregin/DUPLICATE_REGIONS minimum bit score is <x>",                             "dupregin/DUPLICATE_REGIONS minimum bit score is <x>",                             \%opt_HH, \@opt_order_A);
opt_Add("--indefstr",   "real",      25,        $g,   undef,   undef,            "indfstrn/INDEFINITE_STRAND minimum weaker strand bit score is <x>",               "indfstrn/INDEFINITE_STRAND minimum weaker strand bit score is <x>",               \%opt_HH, \@opt_order_A);
opt_Add("--lowsim5seq", "integer",  15,         $g,   undef,   undef,            "lowsim5s/LOW_SIMILARITY_START minimum length is <n>",                             "lowsim5s/LOW_SIMILARITY_START minimum length is <n>",                             \%opt_HH, \@opt_order_A);
opt_Add("--lowsim3seq", "integer",  15,         $g,   undef,   undef,            "lowsim3s/LOW_SIMILARITY_END minimum length is <n>",                               "lowsim3s/LOW_SIMILARITY_END minimum length is <n>",                               \%opt_HH, \@opt_order_A);
opt_Add("--lowsimiseq", "integer",   1,         $g,   undef,   undef,            "lowsimis/LOW_SIMILARITY (internal) minimum length is <n>",                        "lowsimi/LOW_SIMILARITY (internal) minimum length is <n>",                         \%opt_HH, \@opt_order_A);
opt_Add("--lowsim5ftr", "integer",   5,         $g,   undef,   undef,            "lowsim5{c,n}/LOW_FEATURE_SIMILARITY_START minimum length is <n>",                 "lowsim5{c,n}/LOW_FEATURE_SIMILARITY_START minimum length is <n>",                 \%opt_HH, \@opt_order_A);
opt_Add("--lowsim3ftr", "integer",   5,         $g,   undef,   undef,            "lowsim3{c,n}/LOW_FEATURE_SIMILARITY_END minimum length is <n>",                   "lowsim3{c,n}/LOW_FEATURE_SIMILARITY_END minimum length is <n>",                   \%opt_HH, \@opt_order_A);
opt_Add("--lowsimiftr", "integer",   1,         $g,   undef,   undef,            "lowsimi{c,n}/LOW_FEATURE_SIMILARITY (internal) minimum length is <n>",            "lowsimi{c,n}/LOW_FEATURE_SIMILARITY (internal) minimum length is <n>",            \%opt_HH, \@opt_order_A);
opt_Add("--biasfract",  "real",      0.25,      $g,   undef,   undef,            "biasdseq/BIASED_SEQUENCE fractional threshold is <x>",                            "biasdseq/BIASED_SEQUENCE fractional threshold is <x>",                            \%opt_HH, \@opt_order_A);
opt_Add("--nmiscftrthr","integer",   3,         $g,   undef,   undef,            "nmiscftr/TOO_MANY_MISC_FEATURES reported if <n> or more misc_features",           "nmiscftr/TOO_MANY_MISC_FEATURES reported if <n> or more misc_features",           \%opt_HH, \@opt_order_A);
opt_Add("--indefann",   "real",      0.8,       $g,   undef,   undef,            "indf{5,3}lc{c,n}/INDEFINITE_ANNOTATION_{START,END} non-mat_peptide min allowed post probability is <x>",         "indf{5,3}lc{c,n}/'INDEFINITE_ANNOTATION_{START,END} non-mat_peptide min allowed post probability is <x>", \%opt_HH, \@opt_order_A);
opt_Add("--indefann_mp","real",      0.6,       $g,   undef,   undef,            "indf{5,3}lc{c,n}/INDEFINITE_ANNOTATION_{START,END} mat_peptide min allowed post probability is <x>",             "indf{5,3}lc{c,n}/'INDEFINITE_ANNOTATION_{START,END} mat_peptide min allowed post probability is <x>", \%opt_HH, \@opt_order_A);
opt_Add("--fstminntt",  "integer",    4,        $g,   undef,   undef,            "fst{hi,lo,uk}cft/POSSIBLE_FRAMESHIFT{_{HIGH,LOW}_CONF,} max allowed terminal nt length w/o alert is <n>",   "fst{hi,lo,uk}cft/POSSIBLE_FRAMESHIFT{_{HIGH,LOW}_CONF,} max allowed terminal nt length w/o alert is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--fstminnti",  "integer",    6,        $g,   undef,   undef,            "fst{hi,lo,uk}cfi/POSSIBLE_FRAMESHIFT{_{HIGH,LOW}_CONF,} max allowed internal nt length w/o alert is <n>",   "fst{hi,lo,uk}cfi/POSSIBLE_FRAMESHIFT{_{HIGH,LOW}_CONF,} max allowed internal nt length w/o alert is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--fsthighthr", "real",      0.8,       $g,   undef,"--glsearch",        "fsthicf{t,i}/POSSIBLE_FRAMESHIFT_HIGH_CONF minimum average probability for alert is <x>",                  "fsthicf{t,i}/POSSIBLE_FRAMESHIFT_HIGH_CONF minimum average probability for alert is <x>", \%opt_HH, \@opt_order_A);
opt_Add("--fstlowthr",  "real",      0.0,       $g,   undef,"--glsearch",        "fstlocf{t,i}/POSSIBLE_FRAMESHIFT_LOW_CONF minimum average probability for alert is <x>",                   "fstlocf{t,i}/POSSIBLE_FRAMESHIFT_LOW_CONF minimum average probability for alert is <x>", \%opt_HH, \@opt_order_A);
opt_Add("--xalntol",    "integer",   5,         $g,   undef,   undef,            "indf{5,3}{st,lg}/INDEFINITE_ANNOTATION_{START,END} max allowed nt diff blastx start/end is <n>",   "indf{5,3}{st,lg}/INDEFINITE_ANNOTATION_{START,END} max allowed nt diff blastx start/end is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--xmaxins",    "integer",   27,        $g,   undef,"--pv_skip,--pv_hmmer", "insertnp/INSERTION_OF_NT max allowed nucleotide insertion length in blastx validation is <n>",     "insertnp/INSERTION_OF_NT max allowed nucleotide insertion length in blastx validation is <n>",   \%opt_HH, \@opt_order_A);
opt_Add("--xmaxdel",    "integer",   27,        $g,   undef,"--pv_skip,--pv_hmmer", "deletinp/DELETION_OF_NT max allowed nucleotide deletion length in blastx validation is <n>",       "deletinp/DELETION_OF_NT max allowed nucleotide deletion length in blastx validation is <n>",     \%opt_HH, \@opt_order_A);
opt_Add("--nmaxins",    "integer",   27,        $g,   undef,   undef,            "insertnn/INSERTION_OF_NT max allowed nucleotide (nt) insertion length in CDS nt alignment is <n>", "insertnn/INSERTION_OF_NT max allowed nucleotide (nt) insertion length in CDS nt alignment is <n>",   \%opt_HH, \@opt_order_A);
opt_Add("--nmaxdel",    "integer",   27,        $g,   undef,   undef,            "deletinn/DELETION_OF_NT max allowed nucleotide (nt) deletion length in CDS nt alignment is <n>",   "deletinn/DELETION_OF_NT max allowed nucleotide (nt) deletion length in CDS nt alignment is <n>",     \%opt_HH, \@opt_order_A);
opt_Add("--xlonescore",  "integer",  80,        $g,   undef,"--pv_skip,--pv_hmmer", "indfantp/INDEFINITE_ANNOTATION min score for a blastx hit not supported by CM analysis is <n>",    "indfantp/INDEFINITE_ANNOTATION min score for a blastx hit not supported by CM analysis is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--hlonescore",  "integer",  10,        $g,"--pv_hmmer","--pv_skip",        "indfantp/INDEFINITE_ANNOTATION min score for a hmmer hit not supported by CM analysis is <n>",     "indfantp/INDEFINITE_ANNOTATION min score for a hmmer hit not supported by CM analysis is <n>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling cmalign alignment stage";
#        option               type default group  requires incompat   preamble-output                                                                help-output    
opt_Add("--mxsize",     "integer", 16000,     $g,    undef,"--glsearch", "set max allowed memory for cmalign to <n> Mb",                                 "set max allowed memory for cmalign to <n> Mb", \%opt_HH, \@opt_order_A);
opt_Add("--tau",        "real",    1E-3,      $g,    undef,"--glsearch", "set the initial tau value for cmalign to <x>",                                 "set the initial tau value for cmalign to <x>", \%opt_HH, \@opt_order_A);
opt_Add("--nofixedtau", "boolean", 0,         $g,    undef,"--glsearch", "do not fix the tau value when running cmalign, allow it to increase if nec",   "do not fix the tau value when running cmalign, allow it to decrease if nec", \%opt_HH, \@opt_order_A);
opt_Add("--nosub",      "boolean", 0,         $g,    undef,"--glsearch", "use alternative alignment strategy for truncated sequences",                   "use alternative alignment strategy for truncated sequences", \%opt_HH, \@opt_order_A);
opt_Add("--noglocal",   "boolean", 0,         $g,"--nosub","--glsearch", "do not run cmalign in glocal mode (run in local mode)",                        "do not run cmalign in glocal mode (run in local mode)", \%opt_HH, \@opt_order_A);
opt_Add("--cmindi",     "boolean", 0,         $g,    undef, "--nkb,--glsearch", "force cmalign to align one seq at a time",                              "force cmalign to align on seq at a time", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling glsearch alignment stage as alternative to cmalign";
#        option               type default group  requires incompat   preamble-output                                                                help-output    
opt_Add("--glsearch",     "boolean", 0,         $g,"--glsearch", undef,      "align with glsearch from the FASTA package, not to a cm with cmalign",         "align with glsearch from the FASTA package, not to a cm with cmalign", \%opt_HH, \@opt_order_A);
opt_Add("--gls_match",    "integer", 5,         $g,"--glsearch", undef,      "set glsearch match score to <n> > 0 with glsearch -r option",                  "set glsearch match score to <n> > 0 with glsearch -r option", \%opt_HH, \@opt_order_A);
opt_Add("--gls_mismatch", "integer", -3,        $g,"--glsearch", undef,      "set glsearch mismatch score to <n> < 0 with glsearch -r option",               "set glsearch mismatch score to <n> < 0 with glsearch -r option", \%opt_HH, \@opt_order_A);
opt_Add("--gls_gapopen",  "integer", -17,       $g,"--glsearch", undef,      "set glsearch gap open score to <n> < 0 with glsearch -f option",               "set glsearch gap open score to <n> < 0 with glsearch -f option", \%opt_HH, \@opt_order_A);
opt_Add("--gls_gapextend","integer", -4,        $g,"--glsearch", undef,      "set glsearch gap extend score to <n> < 0 with glsearch -g option",             "set glsearch gap extend score to <n> < 0 with glsearch -g option", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for controlling blastx protein validation stage";
#        option               type   default  group  requires incompat            preamble-output                                                                                 help-output    
opt_Add("--xmatrix",     "string",   undef,      $g,     undef,"--pv_skip,--pv_hmmer", "use the matrix <s> with blastx (e.g. BLOSUM45)",                                                "use the matrix <s> with blastx (e.g. BLOSUM45)", \%opt_HH, \@opt_order_A);
opt_Add("--xdrop",       "integer",  25,         $g,     undef,"--pv_skip,--pv_hmmer", "set the xdrop value for blastx to <n>",                                                         "set the xdrop value for blastx to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--xnumali",     "integer",  20,         $g,     undef,"--pv_skip,--pv_hmmer", "number of alignments to keep in blastx output and consider if --xlongest is <n>",               "number of alignments to keep in blastx output and consider if --xlongest is <n>", \%opt_HH, \@opt_order_A);
opt_Add("--xlongest",    "boolean",  0,          $g,     undef,"--pv_skip,--pv_hmmer", "keep the longest blastx hit, not the highest scoring one",                                      "keep the longest blastx hit, not the highest scoring one", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for using hmmer instead of blastx for protein validation";
#     option          type       default group   requires    incompat   preamble-output                                     help-output    
opt_Add("--pv_hmmer",    "boolean", 0,        $g,     undef,  "--pv_skip", "use hmmer for protein validation, not blastx",     "use hmmer for protein validation, not blastx", \%opt_HH, \@opt_order_A);
opt_Add("--h_max",    "boolean", 0,        $g, "--pv_hmmer",  "--pv_skip", "use --max option with hmmsearch",                  "use --max option with hmmsearch", \%opt_HH, \@opt_order_A);
opt_Add("--h_minbit", "real",    -10,      $g, "--pv_hmmer",  "--pv_skip", "set minimum hmmsearch bit score threshold to <x>", "set minimum hmmsearch bit score threshold to <x>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options related to blastn-derived seeded alignment acceleration";
#        option               type   default group   requires  incompat        preamble-output                                                        help-output    
opt_Add("-s",             "boolean",      0,   $g,      undef, undef,          "use top-scoring HSP from blastn to seed the alignment",               "use top-scoring HSP from blastn to seed the alignment", \%opt_HH, \@opt_order_A);
opt_Add("--s_blastnws",   "integer",      7,   $g,       "-s", undef,          "for -s, set blastn -word_size <n> to <n>",                            "for -s, set blastn -word_size <n> to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--s_blastnrw",   "integer",      1,   $g,       "-s", undef,          "for -s, set blastn -reward <n> to <n>",                               "for -s, set blastn -reward <n> to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--s_blastnpn",   "integer",     -2,   $g,       "-s", undef,          "for -s, set blastn -penalty <n> to <n>",                              "for -s, set blastn -penalty <n> to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--s_blastngo",   "integer",      2,   $g,       "-s","--s_blastngdf", "for -s, set blastn -gapopen <n> to <n>",                              "for -s, set blastn -gapopen <n> to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--s_blastnge",   "real",         1,   $g,       "-s","--s_blastngdf", "for -s, set blastn -gapextend <x> to <x>",                            "for -s, set blastn -gapextend <x> to <x>", \%opt_HH, \@opt_order_A);
opt_Add("--s_blastngdf",  "boolean",      0,   $g,       "-s", undef,          "for -s, don't use -gapopen/-gapextend w/blastn (use default values)", "for -s, don't use -gapopen/-gapextend w/blastn (use default values)", \%opt_HH, \@opt_order_A);
opt_Add("--s_blastnsc",   "real",      50.0,   $g,       "-s", undef,          "for -s, set blastn minimum HSP score to consider to <x>",             "for -s, set blastn minimum HSP score to consider to <x>", \%opt_HH, \@opt_order_A);
opt_Add("--s_blastntk",   "boolean",      0,   $g,       "-s", undef,          "for -s, set blastn option -task blastn",                              "for -s, set blastn option -task blastn", \%opt_HH, \@opt_order_A);
opt_Add("--s_blastnxd",   "integer",    110,   $g,       "-s", undef,          "for -s, set blastn option -xdrop_gap_final <n> to <n>",               "for -s, set blastn -xdrop_gap_final <n> to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--s_minsgmlen",  "integer",     10,   $g,       "-s", undef,          "for -s, set minimum length of ungapped region in HSP seed to <n>",    "for -s, set minimum length of ungapped region in HSP seed to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--s_allsgm",     "boolean",      0,   $g,       "-s", "--s_minsgmlen", "for -s, keep full HSP seed, do not enforce minimum segment length",  "for -s, keep full HSP seed, do not enforce minimum segment length", \%opt_HH, \@opt_order_A);
opt_Add("--s_ungapsgm",   "boolean",      0,   $g,       "-s", "--s_minsgmlen,--s_allsgm", "for -s, only keep max length ungapped segment of HSP",    "for -s, only keep max length ungapped segment of HSP", \%opt_HH, \@opt_order_A);
opt_Add("--s_startstop",  "boolean",      0,   $g,       "-s", "--s_ungapsgm", "for -s, allow seed to include gaps in start/stop codons",             "for -s, allow seed to include gaps in start/stop codons", \%opt_HH, \@opt_order_A);
opt_Add("--s_overhang",   "integer",    100,   $g,       "-s", undef,          "for -s, set length of nt overhang for subseqs to align to <n>",       "for -s, set length of nt overhang for subseqs to align to <n>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options related to replacing Ns with expected nucleotides";
#        option               type   default group requires incompat  preamble-output                                                              help-output    
opt_Add("-r",             "boolean",      0,   $g,   undef, undef,    "replace stretches of Ns with expected nts, where possible",                 "replace stretches of Ns with expected nts, where possible",                \%opt_HH, \@opt_order_A);
opt_Add("--r_minlen",     "integer",      5,   $g,    "-r", undef,    "minimum length subsequence to replace Ns in is <n>",                        "minimum length subsequence to replace Ns in is <n>",                       \%opt_HH, \@opt_order_A);
opt_Add("--r_minfract5",     "real",   0.25,   $g,    "-r", undef,    "minimum fraction of Ns in subseq at 5' end to trigger replacement is <x>",  "minimum fraction of Ns in subseq at 5' end to trigger replacement is <x>", \%opt_HH, \@opt_order_A);
opt_Add("--r_minfract3",     "real",   0.25,   $g,    "-r", undef,    "minimum fraction of Ns in subseq at 3' end to trigger replacement is <x>",  "minimum fraction of Ns in subseq at 3' end to trigger replacement is <x>", \%opt_HH, \@opt_order_A);
opt_Add("--r_minfracti",     "real",    0.5,   $g,    "-r", undef,    "minimum fraction of Ns in internal subseq to trigger replacement is <x>",   "minimum fraction of Ns in internal subseq to trigger replacement is <x>",  \%opt_HH, \@opt_order_A);
opt_Add("--r_fetchr",     "boolean",      0,   $g,    "-r", undef,    "fetch features for output fastas from seqs w/Ns replaced, not originals",   "fetch features for output fastas from seqs w/Ns replaced, not originals", \%opt_HH, \@opt_order_A);
opt_Add("--r_cdsmpr",     "boolean",      0,   $g,    "-r", undef,    "detect CDS and MP alerts in sequences w/Ns replaced, not originals",        "detect CDS and MP alerts in sequences w/Ns replaced, not originals",      \%opt_HH, \@opt_order_A);
opt_Add("--r_pvorig",     "boolean",      0,   $g,    "-r", undef,    "use original sequences for protein validation step, not replaced seqs",     "use original sequences for protein validation, not replaced seqs",        \%opt_HH, \@opt_order_A);
opt_Add("--r_prof",       "boolean",      0,   $g,    "-r", undef,    "use slower profile methods, not blastn, to identify Ns to replace",         "use slower profile methods, not blastn, to identify Ns to replace",       \%opt_HH, \@opt_order_A);
opt_Add("--r_list",       "string",   undef,   $g,    "-r", undef,    "with -r, only use models listed in file <s> for N replacement stage",       "with -r, only use models listed in file <s> for N replacement stage",     \%opt_HH, \@opt_order_A);
opt_Add("--r_only",       "string",   undef,   $g,    "-r","--r_list","with -r, only use model named <s> for N replacement stage",                 "with -r, only use model named <s> for N replacement stage",               \%opt_HH, \@opt_order_A);
opt_Add("--r_blastnws",   "integer",      7,   $g,    "-r", undef,          "for -r, set blastn -word_size <n> to <n>",                            "for -r, set blastn -word_size <n> to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--r_blastnrw",   "integer",      1,   $g,    "-r", undef,          "for -r, set blastn -reward <n> to <n>",                               "for -r, set blastn -reward <n> to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--r_blastnpn",   "integer",     -2,   $g,    "-r", undef,          "for -r, set blastn -penalty <n> to <n>",                              "for -r, set blastn -penalty <n> to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--r_blastngo",   "integer",      2,   $g,    "-r","--r_blastngdf", "for -r, set blastn -gapopen <n> to <n>",                              "for -r, set blastn -gapopen <n> to <n>", \%opt_HH, \@opt_order_A);
opt_Add("--r_blastnge",   "real",         1,   $g,    "-r","--r_blastngdf", "for -r, set blastn -gapextend <x> to <x>",                            "for -r, set blastn -gapextend <x> to <x>", \%opt_HH, \@opt_order_A);
opt_Add("--r_blastngdf",  "boolean",      0,   $g,    "-r", undef,          "for -r, don't use -gapopen/-gapextend w/blastn (use default values)", "for -r, don't use -gapopen/-gapextend w/blastn (use default values)", \%opt_HH, \@opt_order_A);
opt_Add("--r_blastnsc",   "real",      50.0,   $g,    "-r", undef,          "for -r, set blastn minimum HSP score to consider to <x>",             "for -r, set blastn minimum HSP score to consider to <x>", \%opt_HH, \@opt_order_A);
opt_Add("--r_blastntk",   "boolean",      0,   $g,    "-r", undef,          "for -r, set blastn option -task blastn",                              "for -r, set blastn option -task blastn", \%opt_HH, \@opt_order_A);
opt_Add("--r_blastnxd",   "integer",    110,   $g,    "-r", undef,          "for -r, set blastn option -xdrop_gap_final <n> to <n>",               "for -r, set blastn -xdrop_gap_final <n> to <n>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options related to splitting input file into chunks and processing each chunk separately";
#     option            type       default  group   requires incompat    preamble-output                                                          help-output    
opt_Add("--split",      "boolean", 0,          $g,    undef,  "-p",       "split input file into chunks, run each chunk separately",              "split input file into chunks, run each chunk separately", \%opt_HH, \@opt_order_A);
opt_Add("--cpu",        "integer", 1,          $g,    undef, undef,       "parallelize across <n> CPU workers (requires --split or --glsearch)",  "parallelize across <n> CPU workers (requires --split or --glsearch)", \%opt_HH, \@opt_order_A);
opt_Add("--sidx",       "integer", 1,          $g,    undef,"--split",    "start sequence indexing at <n> in tabular output files",               "start sequence indexing at <n> in tabular output files", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options related to parallelization on compute farm";
#     option            type       default  group   requires incompat    preamble-output                                                help-output    
opt_Add("-p",           "boolean", 0,          $g,    undef,  undef,      "parallelize cmsearch/cmalign on a compute farm",              "parallelize cmsearch/cmalign on a compute farm", \%opt_HH, \@opt_order_A);
opt_Add("-q",           "string",  undef,      $g,     "-p",  undef,      "use qsub info file <s> instead of default",                   "use qsub info file <s> instead of default", \%opt_HH, \@opt_order_A);
opt_Add("--errcheck",   "boolean", 0,          $g,     "-p",  undef,      "consider any farm stderr output as indicating a job failure", "consider any farm stderr output as indicating a job failure", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options related to splitting input and parallelization on compute farm";
opt_Add("--wait",       "integer", 500,        $g,    undef,  undef,      "allow <n> minutes for jobs on farm",                          "allow <n> wall-clock minutes for jobs on farm to finish, including queueing time", \%opt_HH, \@opt_order_A);
opt_Add("--maxnjobs",   "integer", 2500,       $g,    undef,  undef,      "maximum allowed number of jobs for compute farm",             "set max number of jobs to submit to compute farm to <n>", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "options for skipping stages";
#     option               type       default group   requires    incompat                        preamble-output                                            help-output    
opt_Add("--pv_skip",       "boolean", 0,         $g,   undef,      undef,                         "do not perform blastx-based protein validation",          "do not perform blastx-based protein validation", \%opt_HH, \@opt_order_A);
opt_Add("--align_skip",    "boolean", 0,         $g,   undef,      "-f",                          "skip the alignment step, use existing results",           "skip the alignment step, use results from an earlier run of the script", \%opt_HH, \@opt_order_A);
opt_Add("--val_only",      "boolean", 0,         $g,   undef,      undef,                         "validate CM and other input files and exit",              "validate CM and other input files and exit", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "optional output files";
#       option       type       default   group  requires incompat  preamble-output                                                      help-output    
opt_Add("--out_stk",        "boolean", 0,    $g,    undef,"--keep", "output per-model full length stockholm alignments (.stk)",          "output per-model full length stockholm alignments (.stk)",          \%opt_HH, \@opt_order_A);
opt_Add("--out_afa",        "boolean", 0,    $g,    undef,"--keep", "output per-model full length fasta alignments (.afa)",              "output per-model full length fasta alignments (.afa)",              \%opt_HH, \@opt_order_A);
opt_Add("--out_rpstk",      "boolean", 0,    $g,     "-r","--keep", "with -r, output stockholm alignments of seqs with Ns replaced",     "with -r, output stockholm alignments of seqs with Ns replaced",     \%opt_HH, \@opt_order_A);
opt_Add("--out_rpafa",      "boolean", 0,    $g,     "-r","--keep", "with -r, output fasta alignments of seqs with Ns replaced",         "with -r, output fasta alignments of seqs with Ns replaced",         \%opt_HH, \@opt_order_A);
opt_Add("--out_fsstk",      "boolean", 0,    $g,    undef,"--keep", "additionally output frameshift stockholm alignment files",          "additionally output frameshift stockholm alignment files",          \%opt_HH, \@opt_order_A);
opt_Add("--out_allfasta",   "boolean", 0,    $g,    undef,"--keep", "additionally output fasta files of features",                       "additionally output fasta files of features",                       \%opt_HH, \@opt_order_A);
opt_Add("--out_nofasta",    "boolean", 0,    $g,    undef,"--keep,--out_allfasta", "do not output fasta files of passing/failing seqs",  "do not output fasta files of passing/failing seqs",                 \%opt_HH, \@opt_order_A);
opt_Add("--out_debug",      "boolean", 0,    $g,    undef,"--split","dump voluminous info from various data structures to output files", "dump voluminous info from various data structures to output files", \%opt_HH, \@opt_order_A);

$opt_group_desc_H{++$g} = "other expert options";
#       option            type          default     group  requires incompat  preamble-output                                                          help-output    
opt_Add("--execname",     "string",  undef,         $g,    undef,   undef,    "define executable name of this script as <s>",                           "define executable name of this script as <s>", \%opt_HH, \@opt_order_A);        
opt_Add("--alicheck",     "boolean", 0,             $g,    undef,   undef,    "for debugging, check aligned sequence vs input sequence for identity",   "for debugging, check aligned sequence vs input sequence for identity", \%opt_HH, \@opt_order_A);
opt_Add("--noseqnamemax", "boolean", 0,             $g,    undef,   undef,    "do not enforce a maximum length of 50 for sequence names (GenBank max)", "do not enforce a maximum length of 50 for sequence names (GenBank max)", \%opt_HH, \@opt_order_A);
opt_Add("--minbit",       "real",    -10,           $g,    undef,   undef,    "set minimum cmsearch bit score threshold to <x>",                        "set minimum cmsearch bit score threshold to <x>", \%opt_HH, \@opt_order_A);
opt_Add("--origfa",       "boolean", 0,             $g,    undef,   undef,    "do not copy fasta file prior to analysis, use original",                 "do not copy fasta file prior to analysis, use original", \%opt_HH, \@opt_order_A);
opt_Add("--msub",         "string",  undef,         $g,    undef,   undef,    "read model substitution file from <s>",                                  "read model substitution file from <s>", \%opt_HH, \@opt_order_A);        
opt_Add("--xsub",         "string",  undef,         $g,    undef,   undef,    "read blastx db substitution file from <s>",                              "read blastx db substitution file from <s>", \%opt_HH, \@opt_order_A);
opt_Add("--nodcr",        "boolean", 0,             $g,    undef,   undef,    "do not doctor alignments to shift gaps in start/stop codons",            "do not doctor alignments to shift gaps in start/stop codons", \%opt_HH, \@opt_order_A);
opt_Add("--forcedcrins",  "boolean", 0,             $g,"--cmindi",  undef,    "force insert type alignment doctoring, requires --cmindi",               "force insert type alignment doctoring, requires --cmindi", \%opt_HH, \@opt_order_A);
opt_Add("--xnoid",        "boolean", 0,             $g,    undef,"--pv_hmmer,--pv_skip", "ignore blastx hits that are full length and 100% identical",  "ignore blastx hits that are full length and 100% identical", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay = 
    &GetOptions('h'             => \$GetOptions_H{"-h"}, 
# basic options
                'f'             => \$GetOptions_H{"-f"},
                'v'             => \$GetOptions_H{"-v"},
                'atgonly'       => \$GetOptions_H{"--atgonly"}, 
                'minpvlen=s'    => \$GetOptions_H{"--minpvlen"},
                'nkb=s'         => \$GetOptions_H{"--nkb"}, 
                'keep'          => \$GetOptions_H{"--keep"},
# options for specifiying classification
                'group=s'       => \$GetOptions_H{"--group"},
                'subgroup=s'    => \$GetOptions_H{"--subgroup"},
# options for controlling which alerts cause failure
                "alt_list"         => \$GetOptions_H{"--alt_list"},
                "alt_pass=s"       => \$GetOptions_H{"--alt_pass"},
                "alt_fail=s"       => \$GetOptions_H{"--alt_fail"},
                "alt_mnf_yes=s"    => \$GetOptions_H{"--alt_mnf_yes"},
                "alt_mnf_no=s"     => \$GetOptions_H{"--alt_mnf_no"},
                "ignore_mnf"       => \$GetOptions_H{"--ignore_mnf"},
                "ignore_isdel"     => \$GetOptions_H{"--ignore_isdel"},
                "ignore_afset"     => \$GetOptions_H{"--ignore_afset"},
                "ignore_afsetsubn" => \$GetOptions_H{"--ignore_afsetsubn"},
# options related to model files
                'm=s'           => \$GetOptions_H{"-m"}, 
                'a=s'           => \$GetOptions_H{"-a"}, 
                'i=s'           => \$GetOptions_H{"-i"}, 
                'n=s'           => \$GetOptions_H{"-n"}, 
                'x=s'           => \$GetOptions_H{"-x"}, 
                'mkey=s'        => \$GetOptions_H{"--mkey"}, 
                'mdir=s'        => \$GetOptions_H{"--mdir"}, 
                'mlist=s'       => \$GetOptions_H{"--mlist"}, 
# options for controlling output feature tables
                "nomisc"        => \$GetOptions_H{"--nomisc"},
                "notrim"        => \$GetOptions_H{"--notrim"},
                "noftrtrim=s"   => \$GetOptions_H{"--noftrtrim"},
                "noprotid"      => \$GetOptions_H{"--noprotid"},
                "forceprotid"   => \$GetOptions_H{"--forceprotid"},
# options for controlling alert thresholds
                "lowsc=s"       => \$GetOptions_H{"--lowsc"},
                'indefclass=s'  => \$GetOptions_H{"--indefclass"},
                'incspec=s'     => \$GetOptions_H{"--incspec"},  
                "lowcov=s"      => \$GetOptions_H{"--lowcov"},
                'dupregolp=s'   => \$GetOptions_H{"--dupregolp"},  
                'dupregsc=s'    => \$GetOptions_H{"--dupregsc"},  
                'indefstr=s'    => \$GetOptions_H{"--indefstr"},  
                'lowsim5seq=s'  => \$GetOptions_H{"--lowsim5seq"},
                'lowsim3seq=s'  => \$GetOptions_H{"--lowsim3seq"},
                'lowsimiseq=s'  => \$GetOptions_H{"--lowsimiseq"},
                'lowsim5ftr=s'  => \$GetOptions_H{"--lowsim5ftr"},
                'lowsim3ftr=s'  => \$GetOptions_H{"--lowsim3ftr"},
                'lowsimiftr=s'  => \$GetOptions_H{"--lowsimiftr"},
                'biasfract=s'   => \$GetOptions_H{"--biasfract"},  
                'nmiscftrthr=s' => \$GetOptions_H{"--nmiscftrthr"},  
                'indefann=s'    => \$GetOptions_H{"--indefann"},  
                'indefann_mp=s' => \$GetOptions_H{"--indefann_mp"},  
                'fstminntt=s'   => \$GetOptions_H{"--fstminntt"},
                'fstminnti=s'   => \$GetOptions_H{"--fstminnti"},
                'fsthighthr=s'  => \$GetOptions_H{"--fsthighthr"},
                'fstlowthr=s'   => \$GetOptions_H{"--fstlowthr"},
                'xalntol=s'     => \$GetOptions_H{"--xalntol"},
                'xmaxins=s'     => \$GetOptions_H{"--xmaxins"},
                'xmaxdel=s'     => \$GetOptions_H{"--xmaxdel"},
                'nmaxins=s'     => \$GetOptions_H{"--nmaxins"},
                'nmaxdel=s'     => \$GetOptions_H{"--nmaxdel"},
                'xlonescore=s'  => \$GetOptions_H{"--xlonescore"},
                'hlonescore=s'  => \$GetOptions_H{"--hlonescore"},
# options for controlling cmalign alignment stage 
                'mxsize=s'      => \$GetOptions_H{"--mxsize"},
                'tau=s'         => \$GetOptions_H{"--tau"},
                'nofixedtau'    => \$GetOptions_H{"--nofixedtau"},
                'nosub'         => \$GetOptions_H{"--nosub"},
                'noglocal'      => \$GetOptions_H{"--noglocal"},
                'cmindi'        => \$GetOptions_H{"--cmindi"},
# options for controlling glsearch alignment stage 
                'glsearch'       => \$GetOptions_H{"--glsearch"},
                'gls_match=s'    => \$GetOptions_H{"--gls_match"},
                'gls_mismatch=s' => \$GetOptions_H{"--gls_mismatch"},
                'gls_gapopen=s'  => \$GetOptions_H{"--gls_gapopen"},
                'gls_gapextend=s'=> \$GetOptions_H{"--gls_gapextend"},
# options for controlling protein blastx protein validation stage
                'xmatrix=s'     => \$GetOptions_H{"--xmatrix"},
                'xdrop=s'       => \$GetOptions_H{"--xdrop"},
                'xnumali=s'     => \$GetOptions_H{"--xnumali"},
                'xlongest'      => \$GetOptions_H{"--xlongest"},
# options for using hmmer instead of blastx for protein validation
                'pv_hmmer'      => \$GetOptions_H{"--pv_hmmer"},
                'h_max'         => \$GetOptions_H{"--h_max"},
                'h_minbit=s'    => \$GetOptions_H{"--h_minbit"},
# options related to blastn-based acceleration
                's'             => \$GetOptions_H{"-s"},
                's_blastnws=s'  => \$GetOptions_H{"--s_blastnws"},
                's_blastnrw=s'  => \$GetOptions_H{"--s_blastnrw"},
                's_blastnpn=s'  => \$GetOptions_H{"--s_blastnpn"},
                's_blastngo=s'  => \$GetOptions_H{"--s_blastngo"},
                's_blastnge=s'  => \$GetOptions_H{"--s_blastnge"},
                's_blastngdf'   => \$GetOptions_H{"--s_blastngdf"},
                's_blastnsc=s'  => \$GetOptions_H{"--s_blastnsc"},
                's_blastntk'    => \$GetOptions_H{"--s_blastntk"},
                's_blastnxd=s'  => \$GetOptions_H{"--s_blastnxd"},
                's_minsgmlen=s' => \$GetOptions_H{"--s_minsgmlen"},
                's_allsgm'      => \$GetOptions_H{"--s_allsgm"},
                's_ungapsgm'    => \$GetOptions_H{"--s_ungapsgm"},
                's_startstop'   => \$GetOptions_H{"--s_startstop"},
                's_overhang=s'  => \$GetOptions_H{"--s_overhang"},
# options related to replacing Ns with expected nucleotides
                'r'             => \$GetOptions_H{"-r"},
                'r_minlen=s'    => \$GetOptions_H{"--r_minlen"},
                'r_minfract5=s' => \$GetOptions_H{"--r_minfract5"},
                'r_minfract3=s' => \$GetOptions_H{"--r_minfract3"},
                'r_minfracti=s' => \$GetOptions_H{"--r_minfracti"},
                'r_fetchr'      => \$GetOptions_H{"--r_fetchr"},
                'r_cdsmpr'      => \$GetOptions_H{"--r_cdsmpr"},
                'r_pvorig'      => \$GetOptions_H{"--r_pvorig"},
                'r_prof'        => \$GetOptions_H{"--r_prof"},
                'r_list=s'      => \$GetOptions_H{"--r_list"},
                'r_only=s'      => \$GetOptions_H{"--r_only"},
                'r_blastnws=s'  => \$GetOptions_H{"--r_blastnws"},
                'r_blastnrw=s'  => \$GetOptions_H{"--r_blastnrw"},
                'r_blastnpn=s'  => \$GetOptions_H{"--r_blastnpn"},
                'r_blastngo=s'  => \$GetOptions_H{"--r_blastngo"},
                'r_blastnge=s'  => \$GetOptions_H{"--r_blastnge"},
                'r_blastngdf'   => \$GetOptions_H{"--r_blastngdf"},
                'r_blastnsc=s'  => \$GetOptions_H{"--r_blastnsc"},
                'r_blastntk'    => \$GetOptions_H{"--r_blastntk"},
                'r_blastnxd=s'  => \$GetOptions_H{"--r_blastnxd"},
# options related to splitting
                'split'         => \$GetOptions_H{"--split"},
                'cpu=s'         => \$GetOptions_H{"--cpu"}, 
                'sidx=s'        => \$GetOptions_H{"--sidx"}, 
# options related to parallelization
                'p'             => \$GetOptions_H{"-p"},
                'q=s'           => \$GetOptions_H{"-q"},
                'errcheck'      => \$GetOptions_H{"--errcheck"},
# options related to -p or --split
                'wait=s'        => \$GetOptions_H{"--wait"},
                'maxnjobs=s'    => \$GetOptions_H{"--maxnjobs"},
# options for skipping stages
                'pv_skip'       => \$GetOptions_H{"--pv_skip"},
                'align_skip'    => \$GetOptions_H{"--align_skip"},
                'val_only'      => \$GetOptions_H{"--val_only"},
# optional output files
                'out_stk'       => \$GetOptions_H{"--out_stk"}, 
                'out_afa'       => \$GetOptions_H{"--out_afa"}, 
                'out_rpstk'     => \$GetOptions_H{"--out_rpstk"}, 
                'out_rpafa'     => \$GetOptions_H{"--out_rpafa"}, 
                'out_fsstk'     => \$GetOptions_H{"--out_fsstk"}, 
                'out_allfasta'  => \$GetOptions_H{"--out_allfasta"}, 
                'out_nofasta'   => \$GetOptions_H{"--out_nofasta"}, 
                'out_debug'     => \$GetOptions_H{"--out_debug"},
# other expert options
                'execname=s'    => \$GetOptions_H{"--execname"},
                'alicheck'      => \$GetOptions_H{"--alicheck"},
                'noseqnamemax'  => \$GetOptions_H{"--noseqnamemax"},
                'minbit=s'      => \$GetOptions_H{"--minbit"},
                'origfa'        => \$GetOptions_H{"--origfa"},
                'msub=s'        => \$GetOptions_H{"--msub"},
                'xsub=s'        => \$GetOptions_H{"--xsub"},
                'nodcr'         => \$GetOptions_H{"--nodcr"},
                'forcedcrins'   => \$GetOptions_H{"--forcedcrins"},
                'xnoid'         => \$GetOptions_H{"--xnoid"});

my $total_seconds = -1 * ofile_SecondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $execname_opt  = $GetOptions_H{"--execname"};
my $executable    = (defined $execname_opt) ? $execname_opt : "v-annotate.pl";
my $usage         = "Usage: $executable [-options] <fasta file to annotate> <output directory to create>\n";
my $synopsis      = "$executable :: classify and annotate sequences using a model library";
my $date          = scalar localtime();
my $version       = "1.4";
my $releasedate   = "Dec 2021";
my $pkgname       = "VADR";

# make *STDOUT file handle 'hot' so it automatically flushes whenever we print to it
# it is printed to
select *STDOUT;
$| = 1;

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, undef);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

my $do_keep       = opt_Get("--keep", \%opt_HH);
my $do_replace_ns = opt_Get("-r", \%opt_HH);
my $do_nofasta    = opt_Get("--out_nofasta", \%opt_HH);

#######################################
# deal with --alt_list option, if used
#######################################
# initialize error related data structures, we have to do this early, so we can deal with --alt_list opition
my %alt_info_HH = (); 
vdr_AlertInfoInitialize(\%alt_info_HH, undef);
if(opt_IsUsed("--alt_list",\%opt_HH)) { 
  alert_list_option(\%alt_info_HH, $pkgname, $version, $releasedate);
  exit 0;
}

# check that number of command line args is correct
if(scalar(@ARGV) != 2) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do $executable -h\n\n";
  exit(1);
}

my ($orig_in_fa_file, $dir) = (@ARGV);

# enforce that --alt_pass and --alt_fail options are valid
if((opt_IsUsed("--alt_pass", \%opt_HH)) || (opt_IsUsed("--alt_fail", \%opt_HH))) { 
  alert_pass_fail_options(\%alt_info_HH, \%opt_HH);
}

# enforce that --alt_mnf_yes and --alt_mnf_no options are valid
if((opt_IsUsed("--alt_mnf_yes", \%opt_HH)) || (opt_IsUsed("--alt_mnf_no", \%opt_HH))) { 
  alert_misc_not_failure_options(\%alt_info_HH, \%opt_HH);
}

# enforce that --fsthighthr and --fstlowthr values make sense
if(opt_Get("--fsthighthr", \%opt_HH) < opt_Get("--fstlowthr", \%opt_HH)) { 
  if((opt_IsUsed("--fsthighthr", \%opt_HH)) && (opt_IsUsed("--fstlowthr", \%opt_HH))) { 
    die "ERROR if using --fsthighthr <x1> and --fstlowthr <x2>, <x1> must be > <x2>";
  }
  elsif(opt_IsUsed("--fsthighthr", \%opt_HH)) { 
    die "ERROR if using --fsthighthr <x>, <x> must be > " . opt_Get("--fstlowthr", \%opt_HH);
  }
  elsif(opt_IsUsed("--fstlowthr", \%opt_HH)) { 
    die "ERROR if using --fstlowthr <x>, <x> must be < " . opt_Get("--fsthighthr", \%opt_HH);
  }
  else {
    die "ERROR, default value for --fsthighthr (" . opt_Get("--fsthighthr", \%opt_HH) . ") is less than default value for --fstlowthr (" . opt_Get("--fstlowthr", \%opt_HH) . ")";
  }
}

# check for option requirements that sqp_opts is not sophisticated enough
# to check for:
if(opt_IsUsed("--wait", \%opt_HH)) {
  if((! opt_IsUsed("-p", \%opt_HH)) && (! opt_IsUsed("--split", \%opt_HH))) {
    die "ERROR, --wait only makes sense in combination with -p or --split";
  }
}
if(opt_IsUsed("--maxnjobs", \%opt_HH)) {
  if((! opt_IsUsed("-p", \%opt_HH)) && (! opt_IsUsed("--split", \%opt_HH))) {
    die "ERROR, --maxnjobs only makes sense in combination with -p or --split";
  }
}
if(opt_IsUsed("--cpu", \%opt_HH)) {
  if((! opt_IsUsed("--glsearch", \%opt_HH)) && (! opt_IsUsed("--split", \%opt_HH))) {
    die "ERROR, --cpu only makes sense in combination with --glsearch or --split";
  }
}

# if --split and --out_afa,   we require --out_stk 
# if --split and --out_rpafa, we require --out_fpstk 
# this is because we can't merge afa alignments, only stockholm alignments
# so if we want to merge afa alignments we need the stockholm equivalents
if(opt_IsUsed("--split", \%opt_HH)) { 
  if((opt_IsUsed("--out_afa", \%opt_HH)) && (! opt_IsUsed("--out_stk", \%opt_HH))) { 
    die "ERROR, with --split and --out_afa, --out_stk is also required";
  }
  if((opt_IsUsed("--out_rpafa", \%opt_HH)) && (! opt_IsUsed("--out_rpstk", \%opt_HH))) { 
    die "ERROR, with --split and --out_rpafa, --out_rpstk is also required";
  }
}

# if--nmiscftrthr <n> is used, <n> must be >= 2
if((opt_IsUsed("--nmiscftrthr", \%opt_HH)) && (opt_Get("--nmiscftrthr", \%opt_HH) < 2)) { 
  die "ERROR, with --nmiscftrthr <n>, <n> must be >= 2";
}

#######################################################
# determine if we are running blastx, hmmer, and blastn
#######################################################
# set defaults, and change if nec
my $do_pv_blastx = 1; 
my $do_pv_hmmer  = 0;
if(opt_Get("--pv_skip", \%opt_HH)) { 
  $do_pv_blastx = 0;
  $do_pv_hmmer  = 0;
}
elsif(opt_Get("--pv_hmmer", \%opt_HH)) { 
  $do_pv_blastx = 0;
  $do_pv_hmmer  = 1;
}

my $do_blastn_rpn = (opt_Get("-r", \%opt_HH) && (! opt_Get("--r_prof", \%opt_HH))) ? 1 : 0;
my $do_blastn_cls = opt_Get("-s", \%opt_HH) ? 1 : 0;
my $do_blastn_cdt = opt_Get("-s", \%opt_HH) ? 1 : 0;
my $do_blastn_ali = opt_Get("-s", \%opt_HH) ? 1 : 0;
my $do_blastn_any = ($do_blastn_rpn || $do_blastn_cls || $do_blastn_cdt || $do_blastn_ali) ? 1 : 0;
# we have separate flags for each blastn stage even though
# they are all turned on/off with -s in case future changes
# only need some but not all

my $do_glsearch = opt_Get("--glsearch",  \%opt_HH) ? 1 : 0;

#############################
# create the output directory
#############################
my $cmd;               # a command to run with utl_RunCommand()
my @early_cmd_A = ();  # array of commands we run before our log file is opened

if($dir !~ m/\/$/) { $dir =~ s/\/$//; } # remove final '/' if it exists
if(-d $dir) { 
  $cmd = "rm -rf $dir";
  if(opt_Get("-f", \%opt_HH)) { utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR directory named $dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $dir) { 
  $cmd = "rm $dir";
 if(opt_Get("-f", \%opt_HH)) { utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef); push(@early_cmd_A, $cmd); }
  else                       { die "ERROR a file named $dir already exists. Remove it, or use -f to overwrite it."; }
}

# create the dir
$cmd = "mkdir $dir";
utl_RunCommand($cmd, opt_Get("-v", \%opt_HH), 0, undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root = $dir . "/" . $dir_tail . ".vadr";

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("sequence file", "output directory");
my @arg_A      = ($orig_in_fa_file, $dir);
my %extra_H    = ();
$extra_H{"\$VADRSCRIPTSDIR"}  = $env_vadr_scripts_dir;
$extra_H{"\$VADRMODELDIR"}    = $env_vadr_model_dir;
$extra_H{"\$VADRINFERNALDIR"} = $env_vadr_infernal_dir;
$extra_H{"\$VADREASELDIR"}    = $env_vadr_easel_dir;
$extra_H{"\$VADRBIOEASELDIR"} = $env_vadr_bioeasel_dir;
if($do_pv_blastx || $do_blastn_any) { 
  $extra_H{"\$VADRBLASTDIR"} = $env_vadr_blast_dir;
}
if($do_glsearch) { 
  $extra_H{"\$VADRFASTADIR"} = $env_vadr_fasta_dir;
}
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
                         # 2D keys (at least initially)
                         #  "log":  log file of what's output to stdout
                         #  "cmd":  command file with list of all commands executed
                         #  "list": file with list of all output files created

# open the log and command files 
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "log",  $out_root . ".log",      1, 1, "Output printed to screen");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "cmd",  $out_root . ".cmd",      1, 1, "List of executed commands");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".filelist", 1, 1, "List and description of all output files");
my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
my $FH_HR  = $ofile_info_HH{"FH"};
# output files are all open, if we exit after this point, we'll need
# to close these first.

# open optional output files
if(opt_Get("--out_debug", \%opt_HH)) { 
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "ftrinfo",         $out_root . ".ftrinfo",         1, 1, "per-model feature ftr_info_HAH data (created due to --out_debug)");
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "sgminfo",         $out_root . ".sgminfo",         1, 1, "per-model segment sgm_info_HAH data (created due to --out_debug)");
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "altinfo",         $out_root . ".altinfo",         1, 1, "per-alert-code alt_info_HH data (created due to --out_debug)");
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "stgresults",      $out_root . ".stgresults",      1, 1, "per-sequence stg_results_HHH data (created due to --out_debug)");
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "ftrresults",      $out_root . ".ftrresults",      1, 1, "per-sequence, per-feature ftr_results_HHAH data (created due to --out_debug)");
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "sgmresults",      $out_root . ".sgmresults",      1, 1, "per-sequence, per-segment sgm_results_HHAH data (created due to --out_debug)");
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "altseqinstances", $out_root . ".altseqinstances", 1, 1, "per-sequence-alert alt_seq_instances_HH data (created due to --out_debug)");
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "altftrinstances", $out_root . ".altftrinstances", 1, 1, "per-feature-alert alt_ftr_instances_HHH data (created due to --out_debug)");
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "clsoutput",       $out_root . ".clsoutput",       1, 1, "per-sequence cls_output_HH data (created due to --out_debug)");
  if(opt_Get("-s", \%opt_HH)) { 
    ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "sdaoutput",     $out_root . ".sdaoutput",       1, 1, "per-sequence sda_output_HH data (created due to --out_debug and -s)");
  }
  if(opt_Get("-r", \%opt_HH)) { 
    ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "rpnoutput",     $out_root . ".rpnoutput",       1, 1, "per-sequence rpn_output_HH data (created due to --out_debug and -r)");
  }
}

# now we have the log file open, output the banner there too
ofile_OutputBanner($log_FH, $pkgname, $version, $releasedate, $synopsis, $date, \%extra_H);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

my $progress_w = opt_Get("--split", \%opt_HH) ? 105 : 87; # the width of the left hand column in our progress output, hard-coded
my $start_secs = ofile_OutputProgressPrior("Validating input", $progress_w, $log_FH, *STDOUT);

my @to_remove_A   = (); # list of files to remove at end of subroutine, if --keep not used

###########################################
# Validate that we have all the files we need:
# fasta file
utl_FileValidateExistsAndNonEmpty($orig_in_fa_file, "input fasta sequence file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty

my $opt_mdir_used  = opt_IsUsed("--mdir", \%opt_HH);
my $opt_mkey_used  = opt_IsUsed("--mkey", \%opt_HH);
my $opt_mlist_used = opt_IsUsed("--mlist", \%opt_HH);
my $opt_rlist_used = opt_IsUsed("--r_list", \%opt_HH);
my $opt_m_used     = opt_IsUsed("-m", \%opt_HH);
my $opt_a_used     = opt_IsUsed("-a", \%opt_HH);
my $opt_i_used     = opt_IsUsed("-i", \%opt_HH);
my $opt_n_used     = opt_IsUsed("-n", \%opt_HH);
my $opt_x_used     = opt_IsUsed("-x", \%opt_HH);
my $opt_q_used     = opt_IsUsed("-q", \%opt_HH);
my $opt_msub_used  = opt_IsUsed("--msub", \%opt_HH);
my $opt_xsub_used  = opt_IsUsed("--xsub", \%opt_HH);

my $model_key      = opt_Get("--mkey", \%opt_HH); # special case, default value is set in option definition

my $model_dir      = ($opt_mdir_used)  ? opt_Get("--mdir",     \%opt_HH) : $env_vadr_model_dir;
my $model_list     = ($opt_mlist_used) ? opt_Get("--mlist",    \%opt_HH) : undef;
my $replace_list   = ($opt_rlist_used) ? opt_Get("--r_list",   \%opt_HH) : undef;
my $cm_file        = ($opt_m_used)     ? opt_Get("-m",         \%opt_HH) : $model_dir . "/" . $model_key . ".cm";
my $hmm_pt_file    = ($opt_a_used)     ? opt_Get("-a",         \%opt_HH) : $model_dir . "/" . $model_key . ".hmm";
my $minfo_file     = ($opt_i_used)     ? opt_Get("-i",         \%opt_HH) : $model_dir . "/" . $model_key . ".minfo";
my $blastn_db_file = ($opt_n_used)     ? opt_Get("-n",         \%opt_HH) : $model_dir . "/" . $model_key . ".fa";
my $blastx_db_dir  = ($opt_x_used)     ? opt_Get("-x",         \%opt_HH) : $model_dir;
my $qsubinfo_file  = ($opt_q_used)     ? opt_Get("-q",         \%opt_HH) : $env_vadr_scripts_dir . "/vadr.qsubinfo";
my $msub_file      = ($opt_msub_used)  ? opt_Get("--msub",     \%opt_HH) : undef;
my $xsub_file      = ($opt_xsub_used)  ? opt_Get("--xsub",     \%opt_HH) : undef;
my $cm_extra_string       = "";
my $pthmm_extra_string    = "";
my $minfo_extra_string    = "";
my $blastn_extra_string   = "";
my $blastx_extra_string   = "";
my $qsubinfo_extra_string = "";

if($opt_mdir_used)  { $cm_extra_string       .= " --mdir"; $pthmm_extra_string .= " --mdir"; $minfo_extra_string .= " --mdir"; $blastn_extra_string .= " --mdir"; } 
if($opt_mkey_used)  { $cm_extra_string       .= " --mkey"; $pthmm_extra_string .= " --mkey"; $minfo_extra_string .= " --mkey"; $blastn_extra_string .= " --mkey"; } 
if($opt_m_used)     { $cm_extra_string       .= " -m"; }
if($opt_a_used)     { $pthmm_extra_string    .= " -a"; }
if($opt_i_used)     { $minfo_extra_string    .= " -i"; }
if($opt_n_used)     { $blastn_extra_string   .= " -n"; }
if($opt_x_used)     { $blastx_extra_string   .= " -x"; }
if($opt_q_used)     { $qsubinfo_extra_string .= " -q"; }

# check for minfo file which we always need
utl_FileValidateExistsAndNonEmpty($minfo_file,  sprintf("model info file%s",  ($minfo_extra_string  eq "") ? "" : ", due to $minfo_extra_string"), undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty

# only check for cm file if we need it
if((! $do_glsearch) || (opt_Get("--r_prof", \%opt_HH)) || (opt_Get("--val_only", \%opt_HH))) { 
  utl_FileValidateExistsAndNonEmpty($cm_file,  sprintf("CM file%s",  ($cm_extra_string  eq "") ? "" : ", due to $cm_extra_string"), undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
  for my $sfx (".i1f", ".i1i", ".i1m", ".i1p") { 
    utl_FileValidateExistsAndNonEmpty($cm_file . $sfx, "cmpress created $sfx file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
  }
  # cm file must end in .cm, it's how cmalign_or_glsearch*() subroutines
  # determine if they should run cmalign or glsearch.
  if($cm_file !~ m/\.cm$/) { 
    ofile_FAIL("ERROR, CM file name must end in '.cm', but $cm_file does not", $cm_file, 1, $FH_HR);
  }
}

# only check for blastn db file if we need it
if(($do_blastn_any) || ($do_replace_ns) || ($do_glsearch)) { # we always need this file if $do_replace_ns (-r) because we fetch the consensus model sequence from it
  utl_FileValidateExistsAndNonEmpty($blastn_db_file, sprintf("blastn db file%s", ($blastn_extra_string eq "") ? "" : ", due to $blastn_extra_string"), undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
  foreach my $sfx (".nhr", ".nin", ".nsq", ".ndb", ".not", ".nto", ".ntf") { 
    utl_FileValidateExistsAndNonEmpty($blastn_db_file . $sfx, "blastn $sfx file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
  }
  if(($do_glsearch) || (defined $replace_list) || (opt_IsUsed("--r_only", \%opt_HH))) { 
    foreach my $sfx (".ssi") { # for fetching seqs from
      utl_FileValidateExistsAndNonEmpty($blastn_db_file . $sfx, "easel $sfx file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
    }
  }
}

# only check for blastx db if we need it
if($do_pv_blastx) { 
  $blastx_db_dir =~ s/\/$//; # remove trailing '/'
  if(! -d $blastx_db_dir) { 
    ofile_FAIL(sprintf("ERROR, blast db directory $blastx_db_dir%s does not exist", $blastx_extra_string), 1, $FH_HR);
  }
}

# only check for protein hmm file if we need it
if($do_pv_hmmer) { 
  utl_FileValidateExistsAndNonEmpty($hmm_pt_file, sprintf("HMM file%s", ($pthmm_extra_string eq "") ? "" : ", due to $cm_extra_string"), undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
  for my $sfx (".h3f", ".h3i", ".h3m", ".h3p") { 
    utl_FileValidateExistsAndNonEmpty($hmm_pt_file . $sfx, "hmmpress created $sfx file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
  }
}

# only check for qsubinfo file if we need it
my ($qsub_prefix, $qsub_suffix) = (undef, undef);
if(opt_IsUsed("-p", \%opt_HH)) { 
  utl_FileValidateExistsAndNonEmpty($qsubinfo_file,  sprintf("qsub info file%s",  ($qsubinfo_extra_string  eq "") ? "" : ", specified with $qsubinfo_extra_string"), undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
  # parse the qsubinfo file
  ($qsub_prefix, $qsub_suffix) = vdr_ParseQsubFile($qsubinfo_file, $ofile_info_HH{"FH"});
}

# only check for model list file if --mlist used
if(defined $model_list) { 
  utl_FileValidateExistsAndNonEmpty($model_list, "model list file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
}

# only check for model substitution file if --msub used
if(defined $msub_file) { 
  utl_FileValidateExistsAndNonEmpty($msub_file, "model substitution file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
}

# only check for blastx db substitution file if --xsub used
if(defined $xsub_file) { 
  utl_FileValidateExistsAndNonEmpty($xsub_file, "blastx db substitution file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
}

# only check for -r model list file if --r_list used
if(defined $replace_list) { 
  utl_FileValidateExistsAndNonEmpty($replace_list, "replacement model list file", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty
}

###########################
# Parse the model info file
###########################
my @mdl_info_AH  = (); # array of hashes with model info
my %ftr_info_HAH = (); # hash of array of hashes with feature info 
my %sgm_info_HAH = (); # hash of array of hashes with segment info 

my @reqd_mdl_keys_A = ("name", "length");
my @reqd_ftr_keys_A = ("type", "coords");
utl_FileValidateExistsAndNonEmpty($minfo_file, "model info file", undef, 1, $FH_HR);
vdr_ModelInfoFileParse($minfo_file, \@reqd_mdl_keys_A, \@reqd_ftr_keys_A, \@mdl_info_AH, \%ftr_info_HAH, $FH_HR);

# validate %mdl_info_AH
my $nmdl = utl_AHValidate(\@mdl_info_AH, \@reqd_mdl_keys_A, "ERROR reading model info from $minfo_file", $FH_HR);
my $mdl_idx;
# verify feature coords make sense and parent_idx_str is valid
for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
  my $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
  vdr_FeatureInfoValidateCoords(\@{$ftr_info_HAH{$mdl_name}}, $mdl_info_AH[$mdl_idx]{"length"}, $FH_HR); 
}

# if --group or --subgroup used, make sure at least one model has that group/subgroup
my $exp_group    = opt_Get("--group", \%opt_HH);
my $exp_subgroup = opt_Get("--subgroup", \%opt_HH);
if(opt_IsUsed("--group", \%opt_HH)) { 
  if(utl_AHCountKeyValue(\@mdl_info_AH, "group", $exp_group) == 0) { 
    ofile_FAIL("ERROR with --group $exp_group, did not read any models with group defined as $exp_group in model info file:\n$minfo_file", 1, $FH_HR);
  }
}
if(opt_IsUsed("--subgroup", \%opt_HH)) { 
  if(! defined $exp_group) {
    # opt_ValidateSet() will have enforced --subgroup requires --group, but we check again 
    ofile_FAIL("ERROR with --subgroup, the --group option must also be used", 1, $FH_HR);
  }
  if(utl_AHCountKeyValue(\@mdl_info_AH, "subgroup", $exp_subgroup) == 0) { 
    ofile_FAIL("ERROR with --group $exp_group and --subgroup $exp_subgroup,\ndid not read any models with group defined as $exp_group and subgroup defined as $exp_subgroup in model info file:\n$minfo_file", 1, $FH_HR);
  }
}

# if --mlist used ($model_list will be defined) validate all models listed 
# in $model_list are in model info file
if(defined $model_list) { 
  my $err_msg = "";
  my %mdl_name_H      = (); # key is a model name, value is 1
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name_H{$mdl_info_AH[$mdl_idx]{"name"}} = 1;
  }
  my @mdl_name_list_A = (); # all models read from $model_list file
  utl_FileLinesToArray($model_list, 1, \@mdl_name_list_A, $FH_HR);
  foreach my $mdl (@mdl_name_list_A) { 
    if(! defined $mdl_name_H{$mdl}) { 
      $err_msg .= "$mdl\n";
    }
  }
  if($err_msg ne "") { 
    ofile_FAIL("ERROR, the following models listed in $model_list do not exist in model info file $minfo_file:\n$err_msg\n", 1, $FH_HR);
  }
}

# if --msub used ($msub_file) validate all models listed 
# in $msub_file are in model info file
my %mdl_sub_H = ();
if(defined $msub_file) { 
  my $err_msg = validate_and_parse_sub_file($msub_file, \@mdl_info_AH, \%mdl_sub_H, $FH_HR);
  if($err_msg ne "") { 
    ofile_FAIL("ERROR, problem parsing file $msub_file (--msub):\n$err_msg\n", 1, $FH_HR);
  }
}

# if --xsub used ($xsub_file) validate all models listed 
# in $xsub_file are in model info file
my %blastx_sub_H = ();
if(defined $xsub_file) { 
  my $err_msg = validate_and_parse_sub_file($xsub_file, \@mdl_info_AH, \%blastx_sub_H, $FH_HR);
  if($err_msg ne "") { 
    ofile_FAIL("ERROR, problem parsing file $xsub_file (--xsub):\n$err_msg\n", 1, $FH_HR);
  }
}

my @ftr_reqd_keys_A = ("type", "coords");
for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
  my $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
  utl_AHValidate(\@{$ftr_info_HAH{$mdl_name}}, \@ftr_reqd_keys_A, "ERROR reading feature info for model $mdl_name from $minfo_file", $FH_HR);
  vdr_FeatureInfoImputeLength(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  vdr_FeatureInfoInitializeParentIndexStrings(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  vdr_FeatureInfoValidateParentIndexStrings(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  vdr_FeatureInfoImpute3paFtrIdx(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  vdr_FeatureInfoImputeOutname(\@{$ftr_info_HAH{$mdl_name}});
  vdr_FeatureInfoInitializeMiscNotFailure(\@{$ftr_info_HAH{$mdl_name}}, opt_Get("--ignore_mnf", \%opt_HH), $FH_HR);
  vdr_FeatureInfoInitializeIsDeletable(\@{$ftr_info_HAH{$mdl_name}}, opt_Get("--ignore_isdel", \%opt_HH), $FH_HR);
  vdr_FeatureInfoInitializeAlternativeFeatureSet(\@{$ftr_info_HAH{$mdl_name}}, opt_Get("--ignore_afset", \%opt_HH), $FH_HR);
  vdr_FeatureInfoInitializeAlternativeFeatureSetSubstitution(\@{$ftr_info_HAH{$mdl_name}}, (opt_Get("--ignore_afset", \%opt_HH) || opt_Get("--ignore_afsetsubn", \%opt_HH)), $FH_HR);
  vdr_FeatureInfoValidateMiscNotFailure(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  vdr_FeatureInfoValidateIsDeletable(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  vdr_FeatureInfoValidateAlternativeFeatureSet(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  vdr_FeatureInfoValidateAlternativeFeatureSetSubstitution(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
  vdr_SegmentInfoPopulate(\@{$sgm_info_HAH{$mdl_name}}, \@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
}

# if there are any CDS features, validate that the BLAST db files we need exist, if nec
if($do_pv_blastx) { 
  for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    my $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
    my $ncds = vdr_FeatureInfoCountType(\@{$ftr_info_HAH{$mdl_name}}, "CDS"); 
    if($ncds > 0) { 
      if(! defined $mdl_info_AH[$mdl_idx]{"blastdb"}) { 
        ofile_FAIL("ERROR, model $mdl_name has $ncds CDS features, but \"blastdb\" is not defined in model info file:\n$minfo_file\n", 1, $FH_HR);
      }
      my $blastx_db = $blastx_db_dir . "/" . $mdl_info_AH[$mdl_idx]{"blastdb"};
      foreach my $sfx ("", ".phr", ".pin", ".psq", ".pto", ".ptf", ".pot", ".pdb") { 
        if(! -s ($blastx_db . $sfx)) { 
          ofile_FAIL("ERROR, required blastx_db file $blastx_db" . $sfx . " for model $mdl_name does not exist in directory $blastx_db_dir.\nUse -x to specify a different directory.\n", 1, $FH_HR);
        }
      }
      $mdl_info_AH[$mdl_idx]{"blastdbpath"} = $blastx_db;
    }
  }
}
# for any features with if there are any CDS features, validate that the BLAST db files we need exist

# if --val_only used, validate CM file has all models from $minfo_file, then exit
# note this is the only way in which we validate the CM file
if(opt_Get("--val_only", \%opt_HH)) { 
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  my $start_secs = ofile_OutputProgressPrior("Validating CM file contains all models from model info file", $progress_w, $log_FH, *STDOUT);
  # make sure $cm_file includes CMs for all models we just read in $minfo_file, unless --cmval_skip used
  my $cm_name_file = $out_root . ".cm.namelist";
  my $grep_cmd = "grep ^NAME $cm_file | sed 's/^NAME *//' > $cm_name_file";
  utl_RunCommand($grep_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
  my %cm_name_H = ();
  utl_FileLinesToHash($cm_name_file, 1, \%cm_name_H, $FH_HR);
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    my $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
    if(! exists $cm_name_H{$mdl_name}) { 
      ofile_FAIL("ERROR, read model named $mdl_name in model info file ($minfo_file)\nbut a model with that name does not exist in the CM file ($cm_file)", 1, $FH_HR);
    }
  }
  push(@to_remove_A, $cm_name_file);

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  $total_seconds += ofile_SecondsSinceEpoch();
  ofile_OutputConclusionAndCloseFilesOk($total_seconds, $dir, \%ofile_info_HH);
  exit(0); 
}

###########################################
# Copy and validate the input sequence file
###########################################
my $in_fa_file        = undef;
my $blastn_in_fa_file = undef;
if(opt_Get("--origfa", \%opt_HH)) { 
  # --origfa: analyze original fasta file, do not copy it
  $in_fa_file = $orig_in_fa_file;
  if(-e $in_fa_file . ".ssi") { unlink $in_fa_file . ".ssi"}; # remove SSI file if it exists, it may be out of date
}
else { 
  # default: copy original fasta file and analyze that, but don't just copy it, 
  # use 'esl-reformat fasta', this was introduced to sidestep some mysterious 
  # SSI related issues
  $in_fa_file = $out_root . ".in.fa";
  utl_RunCommand($execs_H{"esl-reformat"} . " fasta $orig_in_fa_file > $in_fa_file", opt_Get("-v", \%opt_HH), 0, $FH_HR);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "cp.in.fasta", $in_fa_file, $do_keep, $do_keep, "copy of input fasta file");
  push(@to_remove_A, $in_fa_file);
  push(@to_remove_A, $in_fa_file . ".ssi");
}
if(($do_replace_ns) || ($do_blastn_any) || ($do_glsearch)) { 
  # need a copy of the input fasta file that does not have 
  # descriptions because blast{n,x} AND glsearch do not output sequences 
  # and descriptions in a parseable way (see github issue #4)
  # (actually we don't really need this if --r_prof but we 
  # make it anyway)
  $blastn_in_fa_file = $out_root . ".blastn.fa";
  sqf_FastaFileRemoveDescriptions($in_fa_file, $blastn_in_fa_file, \%ofile_info_HH);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastn.in.fasta", $blastn_in_fa_file, $do_keep, $do_keep, "copy of input fasta file with descriptions removed for blastn");
  push(@to_remove_A, $blastn_in_fa_file);
}  

my $seqstat_file = $out_root . ".seqstat";
my @seq_name_A = (); # [0..$i..$nseq-1]: name of sequence $i in input file
my %seq_len_H = ();  # key: sequence name (guaranteed to be unique), value: seq length
utl_RunCommand($execs_H{"esl-seqstat"} . " --dna -a $in_fa_file > $seqstat_file", opt_Get("-v", \%opt_HH), 0, $FH_HR);
ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "seqstat", $seqstat_file, 1, 1, "esl-seqstat -a output for input fasta file");
sqf_EslSeqstatOptAParse($seqstat_file, \@seq_name_A, \%seq_len_H, $FH_HR);

# make sure that no sequence names exceed our max_length, unless --noseqnamemax used
my $max_seqname_length = 50; # hard-coded
my $lcl_max_seqname_length = $max_seqname_length + length("lcl|"); # NCBI allows length 54 if it starts with lcl|
if(! opt_Get("--noseqnamemax", \%opt_HH)) { 
  foreach my $seq_name (@seq_name_A) { 
    if($seq_name =~ /^lcl\|/) { 
      if(length($seq_name) > $lcl_max_seqname_length) { 
        ofile_FAIL("ERROR, at least one sequence name that begins with lcl| exceeds the maximum GenBank allowed length of $lcl_max_seqname_length\nfor seq names that start 'lcl|' (otherwise max length is $max_seqname_length):\n$seq_name\nTo bypass this restriction, rerun with the --noseqnamemax option enabled.\n", 1, $FH_HR);
      }
    }
    else { 
      if(length($seq_name) > $max_seqname_length) { 
        ofile_FAIL("ERROR, at least one sequence name exceeds the maximum GenBank allowed length of $max_seqname_length:\n$seq_name\nTo bypass this restriction, rerun with the --noseqnamemax option enabled.\n", 1, $FH_HR);
      }
    }
  }
}

# pre-processing complete
###############################

###############################
# if --split, split up the sequence file into chunks and run each chunk separately
my $do_split = opt_Get("--split", \%opt_HH);
if($do_split) {
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  my $ncpu = opt_Get("--cpu", \%opt_HH);

  $start_secs = ofile_OutputProgressPrior(sprintf("Splitting sequence file into chunks to run independently%s",
                                                  ($ncpu > 1) ? " in parallel on $ncpu processors" : ""),
                                          $progress_w, $FH_HR->{"log"}, *STDOUT);

  my $tot_len_nt  = utl_HSumValues(\%seq_len_H);
  my $nchunk_estimate = vdr_SplitNumSeqFiles($tot_len_nt, \%opt_HH);
  my $nchunk = 1; # rewritten if $nchunk_estimate > 1
  my $nseq = scalar(@seq_name_A);
  my @nseqs_per_chunk_A = (); # [0..$nchunk-1] number of sequences in each chunked fasta file

  # update nchunk_estimate to make parallelization as efficient as possible
  if($nseq == 1) { 
    ; # do nothing, we won't split into chunks
  }
  elsif($nseq <= $ncpu) { 
    $nchunk_estimate = -1; # this will put one sequence per chunk
  }
  elsif($nchunk_estimate < $ncpu) { 
    $nchunk_estimate = $ncpu; # sets number of chunks to number of cpus
  }
  else { 
    # ensure number of chunks is a multiple of number of cpus
    while(($nchunk_estimate % $ncpu) != 0) { 
      $nchunk_estimate++;
    }
  } 

  if($nchunk_estimate != 1) { 
    $nchunk = vdr_SplitFastaFile($execs_H{"esl-ssplit"}, $in_fa_file, $nchunk_estimate, \@nseqs_per_chunk_A, \%opt_HH, \%ofile_info_HH);
    # vdr_SplitFastaFile will return the actual number of fasta files created, 
    # which can differ from the requested amount (which is $nchunk_estimate) that we pass in. 
  }
  else { 
    # write_v_annotate_scripts_for_split_mode() knows about the expected 
    # fasta file name in this case (there is no .1 suffix)
    $nseqs_per_chunk_A[0] = scalar(@seq_name_A); # all seqs will be in only seq file
  }


  # write $ncpu scripts that will execute the $nchunk v-annotate.pl jobs
  my @chunk_outdir_A   = (); # output directory names for $nchunk v-annotate.pl jobs
  my @cpu_out_file_AH  = (); # holds name of output files that vdr_WaitForFarmJobsToFinish() will check
                             # to see when all jobs are complete, will be filled in write_v_annotate_scripts_for_split_mode()
  my $out_root_no_vadr = $dir . "/" . $dir_tail;
  my $script_cmd = write_v_annotate_scripts_for_split_mode($nchunk, $ncpu, $in_fa_file, $out_root_no_vadr, \@nseqs_per_chunk_A, 
                                                           \@chunk_outdir_A, \@cpu_out_file_AH, \@to_remove_A, \%opt_HH, \%ofile_info_HH);
  my $nscript = scalar(@cpu_out_file_AH);
  
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  my $ncpu2print = $ncpu;
  if($ncpu2print > $nscript) { $ncpu2print = $nscript; }
  $start_secs = ofile_OutputProgressPrior(sprintf("Executing $nscript script%s to process $nchunk partition(s) of all %d sequence(s)",
                                                  ($ncpu2print > 1) ? "s in parallel on $ncpu processors" : "",
                                                  scalar(@seq_name_A)), 
                                          $progress_w, $FH_HR->{"log"}, *STDOUT);

  # execute the $ncpu scripts
  utl_RunCommand($script_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);


  my $nscripts_finished = 1; # the final script has finished
  if($nscript > 1) { # we may need to wait for the rest of the jobs
    $nscripts_finished = vdr_WaitForFarmJobsToFinish(0, # we're not running cmalign
                                                     0, # do not exit if any err files are written to (blastx outputs warnings to error files sometimes)
                                                     "out", 1, 5, \@cpu_out_file_AH, undef, undef, "[ok]", \%opt_HH, 
                                                     $ofile_info_HH{"FH"});
    if($nscripts_finished != $nscript) { 
      ofile_FAIL(sprintf("ERROR only $nscripts_finished of the $nscript --split scripts are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", \%opt_HH)), 1, $ofile_info_HH{"FH"});
    }
    ofile_OutputString($log_FH, 1, "# "); # necessary because waitForFarmJobsToFinish() creates lines that summarize wait time and so we need a '#' before 'done' printed by ofile_OutputProgressComplete()
  }
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  
  # merge all per-chunk files together
  $start_secs = ofile_OutputProgressPrior("Merging and finalizing output", $progress_w, $FH_HR->{"log"}, *STDOUT);

  # deal with .cmd file first, this one of the more complicated cases
  my $cmd_file = $ofile_info_HH{"fullpath"}{"cmd"};
  my @cmd_filelist_A = (); # array of chunk cmd files to concatenate
  vdr_MergeOutputGetFileList($out_root_no_vadr, ".cmd", 1, \@cmd_filelist_A, \@chunk_outdir_A, $FH_HR);
  my $chunk_cmd_file = $out_root_no_vadr . ".vadr.cmd.chunk";
  utl_ConcatenateListOfFiles(\@cmd_filelist_A, $chunk_cmd_file, "v-annotate.pl:main", \%opt_HH, $FH_HR);
  # close cmd file so we can append to it
  close $ofile_info_HH{"FH"}{"cmd"};
  my $cat_cmd = "cat $chunk_cmd_file >> $cmd_file";
  utl_RunCommand($cat_cmd, opt_Get("-v", \%opt_HH), 0, undef); # pass undef instead of $FH_HR because $cmd file is closed
  utl_RunCommand("echo $cat_cmd >> $cmd_file", opt_Get("-v", \%opt_HH), 0, undef); # record the command we just ran to the $cmd file
  # reopen, so we can add to it again before exiting
  if(! open($ofile_info_HH{"FH"}{"cmd"}, ">>", $ofile_info_HH{"fullpath"}{"cmd"})) { 
    ofile_FAIL("ERROR in v-annotate.pl:main(), unable to re-open " . $ofile_info_HH{"fullpath"}{"cmd"} . " for appending", 1, $FH_HR);
  }
  push(@to_remove_A, $chunk_cmd_file); 

  if(($do_keep) || (opt_Get("--out_fsstk", \%opt_HH))) { 
    vdr_MergeFrameshiftStockholmFiles($out_root_no_vadr, \@mdl_info_AH, \%ftr_info_HAH, \%sgm_info_HAH, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  }
  if(($do_keep) || (opt_Get("--out_stk", \%opt_HH)) || (opt_Get("--out_afa", \%opt_HH)) || (opt_Get("--out_rpstk", \%opt_HH)) || (opt_Get("--out_rpafa", \%opt_HH))) { 
    vdr_MergeAlignments($out_root_no_vadr, \%execs_H, \@mdl_info_AH, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  }
  if(($do_keep) || (opt_Get("--out_allfasta", \%opt_HH))) { 
    vdr_MergePerFeatureFastaFiles($out_root_no_vadr, \@mdl_info_AH, \%ftr_info_HAH, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  }
  
  my $do_check_exists = 1; # require that all files we are expecting to concatenate below exist, if not exit with error message
  vdr_MergeOutputConcatenateOnly($out_root_no_vadr, ".pass.tbl",  "pass_tbl",    "5 column feature table output for passing sequences",  $do_check_exists, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  vdr_MergeOutputConcatenateOnly($out_root_no_vadr, ".fail.tbl",  "fail_tbl",    "5 column feature table output for failing sequences",  $do_check_exists, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  vdr_MergeOutputConcatenateOnly($out_root_no_vadr, ".pass.list", "pass_list",   "list of passing sequences",                            $do_check_exists, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  vdr_MergeOutputConcatenateOnly($out_root_no_vadr, ".fail.list", "fail_list",   "list of failing sequences",                            $do_check_exists, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  vdr_MergeOutputConcatenateOnly($out_root_no_vadr, ".alt.list",  "alerts_list", "list of alerts in the feature tables",                 $do_check_exists, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  if(! opt_Get("--out_nofasta", \%opt_HH)) { 
    vdr_MergeOutputConcatenateOnly($out_root_no_vadr, ".pass.fa",  "pass_fa", "fasta file with passing sequences", $do_check_exists, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
    vdr_MergeOutputConcatenateOnly($out_root_no_vadr, ".fail.fa",  "fail_fa", "fasta file with failing sequences", $do_check_exists, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  }

  # merge files for which we take special care to preserve spacing
  my $nlines_preserve_spacing = 100; # we preserve spacing for up to 100 lines
  my @head_AA = ();  # 2D array with header strings
  my @cljust_A = (); # '1'/'0' array for whether each column is left-justified or not

  helper_tabular_fill_header_and_justification_arrays("ant", \@head_AA, \@cljust_A, $FH_HR);
  vdr_MergeOutputConcatenatePreserveSpacing($out_root_no_vadr, ".sqa", "ant", "per-sequence tabular annotation summary file", $do_check_exists, $nlines_preserve_spacing, "  ", 1, \@head_AA, \@cljust_A, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  
  helper_tabular_fill_header_and_justification_arrays("cls", \@head_AA, \@cljust_A, $FH_HR);
  vdr_MergeOutputConcatenatePreserveSpacing($out_root_no_vadr, ".sqc", "cls", "per-sequence tabular classification summary file", $do_check_exists, $nlines_preserve_spacing, "  ", 1, \@head_AA, \@cljust_A, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);

  helper_tabular_fill_header_and_justification_arrays("ftr", \@head_AA, \@cljust_A, $FH_HR);
  vdr_MergeOutputConcatenatePreserveSpacing($out_root_no_vadr, ".ftr", "ftr", "per-feature tabular summary file", $do_check_exists, $nlines_preserve_spacing, "  ", 1, \@head_AA, \@cljust_A, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);

  helper_tabular_fill_header_and_justification_arrays("sgm", \@head_AA, \@cljust_A, $FH_HR);
  vdr_MergeOutputConcatenatePreserveSpacing($out_root_no_vadr, ".sgm", "sgm", "per-model-segment tabular summary file", $do_check_exists, $nlines_preserve_spacing, "  ", 1, \@head_AA, \@cljust_A, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);

  helper_tabular_fill_header_and_justification_arrays("alt", \@head_AA, \@cljust_A, $FH_HR);
  vdr_MergeOutputConcatenatePreserveSpacing($out_root_no_vadr, ".alt", "alt", "per-alert tabular summary file", $do_check_exists, $nlines_preserve_spacing, "  ", 1, \@head_AA, \@cljust_A, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);

  my $zero_alt = vdr_MergeOutputAlcTabularFile ($out_root_no_vadr, \%alt_info_HH,"alert count tabular summary file", $do_check_exists, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);

  vdr_MergeOutputMdlTabularFile ($out_root_no_vadr, "per-model tabular summary file", $do_check_exists, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);

  helper_tabular_fill_header_and_justification_arrays("dcr", \@head_AA, \@cljust_A, $FH_HR);
  vdr_MergeOutputConcatenatePreserveSpacing($out_root_no_vadr, ".dcr", "dcr", "alignment doctoring tabular summary file", $do_check_exists, $nlines_preserve_spacing, "  ", 0, \@head_AA, \@cljust_A, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);

  if($do_blastn_ali) {
    helper_tabular_fill_header_and_justification_arrays("sda", \@head_AA, \@cljust_A, $FH_HR);
    vdr_MergeOutputConcatenatePreserveSpacing($out_root_no_vadr, ".sda", "sda", "seed alignment summary file (-s)", $do_check_exists, $nlines_preserve_spacing, "  ", 1, \@head_AA, \@cljust_A, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  }
  if($do_replace_ns) { 
    helper_tabular_fill_header_and_justification_arrays("rpn", \@head_AA, \@cljust_A, $FH_HR);
    vdr_MergeOutputConcatenatePreserveSpacing($out_root_no_vadr, ".rpn", "rpn", "replaced stretches of Ns summary file (-r)", $do_check_exists, $nlines_preserve_spacing, "  ", 1, \@head_AA, \@cljust_A, \@chunk_outdir_A, \%opt_HH, \%ofile_info_HH);
  }

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  output_mdl_and_alc_files_and_remove_temp_files($zero_alt, \@to_remove_A, \%opt_HH, \%ofile_info_HH);

  # remove per-chunk directories, unless --keep
  if(! opt_Get("--keep", \%opt_HH)) { 
    foreach my $chunk_outdir (@chunk_outdir_A) { 
      utl_RunCommand("rm $chunk_outdir/*", opt_Get("-v", \%opt_HH), 0, $FH_HR); 
      rmdir $chunk_outdir; # we call this using Perl's 'rmdir' instead of utl_RunCommand() because sometimes a subdir 
                           # won't get removed (b/c it's not empty) for reasons I don't understand (it seems like it 
                           # actually is empty), using Perl's rmdir doesn't print an error if dir is not removed.
      print $cmd_FH "rmdir $chunk_outdir\n";
    }
  }

  $total_seconds += ofile_SecondsSinceEpoch();
  ofile_OutputConclusionAndCloseFilesOk($total_seconds, $dir, \%ofile_info_HH);

  exit 0;
} # end of 'if($do_split)'
###############################

# open the sequence file into a Bio::Easel::SqFile object
my $in_sqfile  = Bio::Easel::SqFile->new({ fileLocation => $in_fa_file }); # the sequence file object
my $rpn_sqfile = undef;
# open the blastn_db sequence file too, if we need it
my $blastn_db_sqfile = undef;
my $r_subset_blastn_db_file = undef;
if(($do_blastn_any) || ($do_replace_ns) || ($do_glsearch)) { 
  $blastn_db_sqfile = Bio::Easel::SqFile->new({ fileLocation => $blastn_db_file });

  # deal with the --r_list or --r_only option (they are incompatible so we can only have one or the other)
  if((defined $replace_list) || (opt_IsUsed("--r_only", \%opt_HH))) { 
    # create the blastn db file that we will use with -r, which must
    # be a subset of seqs in $blastn_db_file (if not, script will exit in error in the fetch_seqs_given_names
    # call below
    my @r_subset_seq_name_A = (); # array of sequence names read from --r_list file
    if(defined $replace_list) { 
      utl_FileLinesToArray($replace_list, 1, \@r_subset_seq_name_A, $FH_HR);
    }
    else { # opt_IsUsed("--r_only", \%opt_HH)
      @r_subset_seq_name_A = (opt_Get("--r_only", \%opt_HH));
    }
    $r_subset_blastn_db_file = $out_root . ".r_subset.blastn.fa";
    $blastn_db_sqfile->fetch_seqs_given_names(\@r_subset_seq_name_A, 60, $r_subset_blastn_db_file);
    sqf_BlastDbCreate($execs_H{"makeblastdb"}, "nucl", $r_subset_blastn_db_file, \%opt_HH, $FH_HR);
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "r_subset.blastn.fa",  $r_subset_blastn_db_file,          0, $do_keep, "--r_list blastn db fasta file");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "r_subset.blastn.nhr", $r_subset_blastn_db_file . ".nhr", 0, $do_keep, "--r_list blastn db .nhr file");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "r_subset.blastn.nin", $r_subset_blastn_db_file . ".nin", 0, $do_keep, "--r_list blastn db .nin file");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "r_subset.blastn.nsq", $r_subset_blastn_db_file . ".nsq", 0, $do_keep, "--r_list blastn db .nsq file");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "r_subset.blastn.ndb", $r_subset_blastn_db_file . ".ndb", 0, $do_keep, "--r_list blastn db .ndb file");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "r_subset.blastn.not", $r_subset_blastn_db_file . ".not", 0, $do_keep, "--r_list blastn db .not file");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "r_subset.blastn.ntf", $r_subset_blastn_db_file . ".ntf", 0, $do_keep, "--r_list blastn db .ntf file");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "r_subset.blastn.nto", $r_subset_blastn_db_file . ".nto", 0, $do_keep, "--r_list blastn db .nto file");
    push(@to_remove_A, $r_subset_blastn_db_file);
    push(@to_remove_A, $r_subset_blastn_db_file . ".nhr");
    push(@to_remove_A, $r_subset_blastn_db_file . ".nin");
    push(@to_remove_A, $r_subset_blastn_db_file . ".nsq");
    push(@to_remove_A, $r_subset_blastn_db_file . ".ndb");
    push(@to_remove_A, $r_subset_blastn_db_file . ".not");
    push(@to_remove_A, $r_subset_blastn_db_file . ".ntf");
    push(@to_remove_A, $r_subset_blastn_db_file . ".nto");
  }
} 

# Initialize the classification results
my %alt_seq_instances_HH = (); # 2D key with info on all instances of per-sequence alerts 
                               # key1: sequence name, key2 alert code, value: alert message

# Add ambgnt5s and ambgnt3s alerts, if any
alert_add_ambgnt5s_ambgnt3s(\$in_sqfile, \@seq_name_A, \%seq_len_H, \%alt_seq_instances_HH, \%alt_info_HH, \%opt_HH, \%ofile_info_HH);
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

my %stg_results_HHH = (); # key 1: sequence name, 
                          # key 2: ("{pre,std}.cls.1","{pre,std}.cls.2","{pre,std}.cls.eg","{pre,std}.cdt.bs", "{pre,std}.cdt.os")
                          #        where 'pre' implies pre-processing replace-ambiguities stage, and 'std' implies standard stage
                          # key 3: ("model", "coords", "bstrand", "score", "bias")

###############################################################
# If -r, pre-processing to identify and replace stretches of Ns
###############################################################
# -r related output for .rpn file
my %rpn_output_HH = (); # 2D key with info to output related to the  option
                        # key1: sequence name, key2 various stats (see output_tabular())
my $rpn_fa_file = undef;
if($do_replace_ns) { 
  my %seq_replaced_H = ();
  my %mdl_seq_name_HA = ();
  classification_stage(\%execs_H, "rpn.cls", $cm_file,
                       ((defined $r_subset_blastn_db_file) ? $r_subset_blastn_db_file : $blastn_db_file),
                       $blastn_in_fa_file, \%seq_len_H,
                       $qsub_prefix, $qsub_suffix, \@mdl_info_AH, \%stg_results_HHH, 
                       $out_root, $progress_w, \@to_remove_A, \%opt_HH, \%ofile_info_HH);
  coverage_determination_stage(\%execs_H, "rpn.cdt", $cm_file, \$in_sqfile, \@seq_name_A, \%seq_len_H,
                               $qsub_prefix, $qsub_suffix, \@mdl_info_AH, \%mdl_seq_name_HA, undef, \%stg_results_HHH, 
                               ((opt_IsUsed("--msub", \%opt_HH)) ? \%mdl_sub_H : undef),
                               $out_root, $progress_w, \@to_remove_A, \%opt_HH, \%ofile_info_HH);

  my $start_secs = ofile_OutputProgressPrior(sprintf("Replacing Ns based on results of %s-based pre-processing", "blastn"), $progress_w, $FH_HR->{"log"}, *STDOUT);
  my $rpn_subset_fa_file = $out_root . ".rpn.sub.fa";
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "rpn.sub.fa", $rpn_subset_fa_file, 0, $do_keep, "fasta file with sequences for which Ns were replaced");
  push(@to_remove_A, $rpn_subset_fa_file);
  push(@to_remove_A, $rpn_subset_fa_file.".ssi");
  
  # for each model with seqs to align to, create the sequence file and run cmalign/glsearch
  my $mdl_name;
  my $nseq_replaced = 0;
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
    if(defined $mdl_seq_name_HA{$mdl_name}) { 
      my $tblout_file = $ofile_info_HH{"fullpath"}{"rpn.cdt.$mdl_name.tblout"};
      $nseq_replaced += parse_cdt_tblout_file_and_replace_ns($tblout_file, \$in_sqfile, \$blastn_db_sqfile, \@mdl_info_AH, $mdl_name, $mdl_idx,
                                                             \@seq_name_A, \%seq_len_H, \%seq_replaced_H, \%rpn_output_HH, $out_root, \%opt_HH, \%ofile_info_HH);
    }
  }
  close($ofile_info_HH{"FH"}{"rpn.sub.fa"}); 

  if($nseq_replaced > 0) { 
    my $rpn_subset_sqfile = Bio::Easel::SqFile->new({ fileLocation => $rpn_subset_fa_file }); # the sequence file object

    # create a new file with original sequences for those that did not have any Ns replaced
    # and doctored sequences with Ns replaced for those that did
    # sequences are named identically and in the same order as in the input fasta file
    $rpn_fa_file = $out_root . ".rpn.fa";
    ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "rpn.fa", $rpn_fa_file, 0, $do_keep, sprintf("fasta file with all sequences, %d with Ns replaced", $nseq_replaced));
    push(@to_remove_A, $rpn_fa_file);
    push(@to_remove_A, $rpn_fa_file.".ssi");
    my $fa_FH = $ofile_info_HH{"FH"}{"rpn.fa"};
    foreach my $seq_name (@seq_name_A) { 
      if(defined $seq_replaced_H{$seq_name}) { 
        print $fa_FH $rpn_subset_sqfile->fetch_seq_to_fasta_string($seq_name);
      }
      else { 
        print $fa_FH $in_sqfile->fetch_seq_to_fasta_string($seq_name);
      }
    }
    close($ofile_info_HH{"FH"}{"rpn.fa"}); 

    # replace main $fa_file we will work with subsequently
    $rpn_sqfile = Bio::Easel::SqFile->new({ fileLocation => $rpn_fa_file });
  }
  ofile_OutputProgressComplete($start_secs, undef, $FH_HR->{"log"}, *STDOUT);
} # end of 'if($do_replace_ns)'

# determine the fasta file we'll use for remaining analysis:
# --
# if((! $do_blastn_cls) && (! $do_glsearch)): 
#   if(defined $rpn_fa_file): use $rpn_fa_file (will only be true if $do_replace_ns AND at least one seq was replaced)
#   else:                     use $in_fa_file
# NOTE: even if $do_replace_ns, if no seqs were replaced then $rpn_fa_file will be undef
# --
# if($do_blastn_cls || $do_glsearch): 
#   if(defined $rpn_fa_file) { (will only be true if $do_replace_ns AND at least one seq was replaced)
#       need to make copy of $rpn_fa_file of $in_fa_file
#       with descriptions removed to avoid the issue with 
#       ambiguity between seq names and seq descriptions 
#       in blastn parsing (github issue #4).
#   else { # ! defined $rpn_fa_file
#       use $blastn_in_fa_file, a copy 
#       of seq file with descriptions removed that we already
#       made above
# --
my $fa_file_for_analysis = undef; 
my $blastn_rpn_fa_file   = undef;
if((! $do_blastn_cls) && (! $do_glsearch)) { 
  $fa_file_for_analysis = (defined $rpn_fa_file) ? $rpn_fa_file : $in_fa_file;
}
else { 
  if(defined $rpn_fa_file) { 
    # make copy with descriptions removed
    $blastn_rpn_fa_file = $out_root . ".blastn.rpn.fa";
    sqf_FastaFileRemoveDescriptions($rpn_fa_file, $blastn_rpn_fa_file, \%ofile_info_HH);
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "blastn.rpn.fa", $blastn_rpn_fa_file, 0, $do_keep, "copy of input fasta file with Ns replaced with descriptions removed for blastn");
    push(@to_remove_A, $blastn_rpn_fa_file);
    push(@to_remove_A, $blastn_rpn_fa_file.".ssi");
    $fa_file_for_analysis = $blastn_rpn_fa_file;
  }
  else { 
    $fa_file_for_analysis = $blastn_in_fa_file;
  }
}  
my $glsearch_sqfile = undef;
if($do_glsearch) { 
  $glsearch_sqfile = Bio::Easel::SqFile->new({ fileLocation => $fa_file_for_analysis }); # the sequence file object
}

# set up sqfile values for analysis, feature fetching and cds and mp alerts
# this is independent of whether we are doing blastn or not because these 
# are the seqfiles used for fetching not blastn analysis
my $sqfile_for_cds_mp_alerts_R = ((defined $rpn_sqfile) && (opt_Get("--r_cdsmpr", \%opt_HH)))   ? \$rpn_sqfile : \$in_sqfile; # sqfile we'll fetch from to analyze CDS and mature peptide features
my $sqfile_for_pv_R            = ((defined $rpn_sqfile) && (! opt_Get("--r_pvorig", \%opt_HH))) ? \$rpn_sqfile : \$in_sqfile; # sqfile we'll fetch from for protein validation
my $sqfile_for_output_fastas_R = ((defined $rpn_sqfile) && (opt_Get("--r_fetchr", \%opt_HH)))   ? \$rpn_sqfile : \$in_sqfile; # sqfile we'll fetch from to make per-feature output fastas 

# if --glsearch we need to do alignment with description-less sequence file (github issue #33)
# $glsearch_sqfile will be either $blastn_rpn_fa_file or $blastn_in_fa_file (see block that
# defined $fa_file_for_analysis above
my $sqfile_for_analysis_R = undef;
if(defined $glsearch_sqfile) { 
  $sqfile_for_analysis_R = \$glsearch_sqfile;
}
elsif(defined $rpn_sqfile) { 
  $sqfile_for_analysis_R = \$rpn_sqfile;
}
else { 
  $sqfile_for_analysis_R = \$in_sqfile;
}

# determine if we need to create separate files with cds seqs for the protein validation stage
# if -r and we replaced at least one sequence, we do (actually, for some combinations of 
# --r_cdsmpr and --r_fetchr we don't need to, but we do anyway because it's more complicated to
# check than it is worth)
my $do_separate_cds_fa_files_for_protein_validation = (($do_replace_ns) && (defined $rpn_sqfile)) ? 1 : 0; 
    
####################################
# Classification: cmsearch round 1
####################################
classification_stage(\%execs_H, "std.cls", $cm_file, $blastn_db_file, $fa_file_for_analysis, \%seq_len_H,
                     $qsub_prefix, $qsub_suffix, \@mdl_info_AH, \%stg_results_HHH, 
                     $out_root, $progress_w, \@to_remove_A, \%opt_HH, \%ofile_info_HH);

###########################################
# Coverage determination: cmsearch round 2
###########################################
my %mdl_cls_ct_H = ();  # key is model name $mdl_name, value is number of sequences classified to this model
coverage_determination_stage(\%execs_H, "std.cdt", $cm_file, $sqfile_for_analysis_R, \@seq_name_A, \%seq_len_H,
                             $qsub_prefix, $qsub_suffix, \@mdl_info_AH, undef, \%mdl_cls_ct_H, \%stg_results_HHH, 
                             ((opt_IsUsed("--msub", \%opt_HH)) ? \%mdl_sub_H : undef),
                             $out_root, $progress_w, \@to_remove_A, \%opt_HH, \%ofile_info_HH);


############################
# Add classification alerts
############################
# add classification errors based on stg_results_HHH
# keep track of seqs to annotate per model
my %cls_output_HH = (); # 2D key with info to output derived from the classification stage
                        # key1: sequence name, key2 one of: "score", "scpnt", "scdiff", "bstrand", "scov", "mcov", "model1", "model2"
add_classification_alerts(\%alt_seq_instances_HH, \%seq_len_H, \@mdl_info_AH, \%alt_info_HH, \%stg_results_HHH, \%cls_output_HH, \%opt_HH, \%ofile_info_HH);

##################
# Align sequences
##################
# create per-model files again, because some classification alerts can cause sequences not to be annotated 
# (e.g. revcompl (reverse complement))
# clear out %mdl_seq_name_HA, %mdl_seq_len_H, %seq2mdl_H, we'll refill them 
my %mdl_seq_name_HA = ();
my %mdl_seq_len_H   = ();
my %seq2mdl_H       = ();
my %mdl_ant_ct_H = ();  # key is model name, value is number of sequences annotated with that model

populate_per_model_data_structures_given_classification_results(\@seq_name_A, \%seq_len_H, "std.aln", \%stg_results_HHH, 
                                                                \%cls_output_HH, \%alt_info_HH, \%alt_seq_instances_HH,
                                                                ((opt_IsUsed("--msub", \%opt_HH)) ? \%mdl_sub_H : undef),
                                                                \%mdl_seq_name_HA, \%mdl_seq_len_H, \%mdl_ant_ct_H, 
                                                                \%seq2mdl_H, $FH_HR);

# file names and data structures necessary for the alignment stage
my %stk_file_HA       = ();     # hash of arrays of all stk output files created in cmalignOrNhmmscanWrapper()
                                # 1D key is $mdl_name
my @overflow_seq_A    = ();     # array of sequences that fail cmalign b/c required matrix was too big
my @overflow_mxsize_A = ();     # array of required matrix sizes for each sequence in @overflow_seq_A
my %ftr_results_HHAH  = ();     # 1st dim: hash, keys are model names
                                # 2nd dim: hash, keys are sequence names
                                # 3rd dim: array, 0..$nsgm-1, one per segment
                                # 4th dim: hash of feature results, keys are:
                                # keys include "n_start", "n_stop", "n_stop_c", "n_strand", "n_5trunc", "n_3trunc"
                                # "p_qstart", "p_qstop", "p_hstart", "p_hstop", "p_strand", "p_query", "p_ins", p_del", "p_trcstop", "p_score", and more
my %sgm_results_HHAH  = ();     # 1st dim: hash, keys are model names
                                # 2nd dim: hash, keys are sequence names
                                # 3rd dim: array, 0..$nsgm-1, one per segment
                                # 4th dim: hash, keys are "sstart", "sstop", "mstart", "mstop", "strand", "5seqflush", "3seqflush", "5trunc", "3trunc"
my %alt_ftr_instances_HHH = (); # hash of arrays of hashes
                                # key1: sequence name, key2: feature index, key3: alert code, value: alert message
my %mdl_unexdivg_H = ();        # key is model name, value is number of unexdivg alerts thrown for that model in alignment stage
my %mdl_unexdivg_HA = ();       # key is model name, value is array of seqs that had unexdivg alerts thrown for that model in alignment stage

my $cur_mdl_fa_file;         # fasta file with sequences to align to current model
my $cur_mdl_align_fa_file;   # fasta file with sequences to align to current model
my $cur_mdl_nseq;            # number of sequences assigned to model
my $cur_mdl_nalign;          # number of sequences we are aligning for current model will be $cur_mdl_nseq unless -s
my $cur_mdl_tot_seq_len;     # sum of total number of nucleotides we are aligning
my %dcr_output_HAH = ();     # hash of array of hashes with info to output related to rare sequences for which the alignment was doctored
                             # more info on this in output_tabular() comments/code

# -s related output for .sda file
my %sda_output_HH = (); # 2D key with info to output related to the -s option
# per-model variables only used if -s used
my %sda_mdl_H     = ();  # key is sequence name, value is mdl coords of seed from blastn alignment
my %sda_seq_H     = ();  # key is sequence name, value is seq coords of seed from blastn alignment
my %seq2subseq_HA = ();  # hash of arrays, key 1: sequence name, array is list of subsequences fetched for this sequence
my %subseq2seq_H  = ();  # hash, key: subsequence name, value is sequence it derives from
my %subseq_len_H  = ();  # key is name of subsequence, value is length of that subsequence

# for each model with seqs to align to, create the sequence file and run cmalign/glsearch
my $mdl_name;

for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
  $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};

  if(defined $mdl_seq_name_HA{$mdl_name}) { 
    %sda_mdl_H     = ();
    %sda_seq_H     = ();
    %seq2subseq_HA = ();
    %subseq2seq_H  = ();
    %subseq_len_H  = ();
    
    $cur_mdl_fa_file = $out_root . "." . $mdl_name . ".a.fa";
    $cur_mdl_nseq = scalar(@{$mdl_seq_name_HA{$mdl_name}});

    my $glsearch_db_file = undef;
    if($do_glsearch) { # create the glsearch db file with a single sequence
      $glsearch_db_file = $out_root . "." . $mdl_name . ".glsearch.fa";
      my @glsearch_seqname_A = ($mdl_name);
      $blastn_db_sqfile->fetch_seqs_given_names(\@glsearch_seqname_A, 60, $glsearch_db_file);
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . ".glsearch.library", $glsearch_db_file, 0, $do_keep, sprintf("glsearch library file for model $mdl_name"));
      push(@to_remove_A, $glsearch_db_file);
    }

    # fetch seqs (we need to do this even if we are not going to send the full seqs to cmalign/glsearch (e.g if $do_blastn_ali))
    $$sqfile_for_analysis_R->fetch_seqs_given_names(\@{$mdl_seq_name_HA{$mdl_name}}, 60, $cur_mdl_fa_file);
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . ".a.fa", $cur_mdl_fa_file, 0, $do_keep, sprintf("%sinput seqs that match best to model $mdl_name", ($do_replace_ns) ? "replaced " : ""));
    push(@to_remove_A, $cur_mdl_fa_file); 

    # set info on seqs we will align, we do this different if $do_blastn_ali or not
    if(! $do_blastn_ali) { 
      $cur_mdl_align_fa_file = $cur_mdl_fa_file;
      $cur_mdl_nalign = $cur_mdl_nseq;
      $cur_mdl_tot_seq_len = utl_HSumValuesSubset(\%seq_len_H, \@{$mdl_seq_name_HA{$mdl_name}});
    }
    else {
      # $do_blastn_ali == 1
      # create the fasta file with sets of subsequences that omit well-defined regions from blastn alignment
      my $indel_file = $ofile_info_HH{"fullpath"}{"std.cdt.$mdl_name.indel"};
      my @subseq_AA = ();
      $cur_mdl_align_fa_file = $out_root . "." . $mdl_name . ".a.subseq.fa";
      my ($start_codon_coords, $stop_codon_coords) = vdr_FeatureInfoCdsStartStopCodonCoords(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR);
      parse_blastn_indel_file_to_get_subseq_info($indel_file, \@{$mdl_seq_name_HA{$mdl_name}}, \%seq_len_H, 
                                                 $mdl_name, $start_codon_coords, $stop_codon_coords, 
                                                 \@subseq_AA, \%sda_mdl_H, \%sda_seq_H, 
                                                 \%seq2subseq_HA, \%subseq2seq_H, \%subseq_len_H, 
                                                 \%opt_HH, \%ofile_info_HH);
      $cur_mdl_nalign = scalar(@subseq_AA);
      if($cur_mdl_nalign > 0) { 
        $$sqfile_for_analysis_R->fetch_subseqs(\@subseq_AA, 60, $cur_mdl_align_fa_file);
        my $subseq_key = $mdl_name . ".a.subseq.fa";
        ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $subseq_key, $cur_mdl_align_fa_file, 0, $do_keep, sprintf("subsequences to align with %s for model $mdl_name (created due to -s)", ($do_glsearch) ? "glsearch" : "cmalign"));
        push(@to_remove_A, $ofile_info_HH{"fullpath"}{$subseq_key});
        $cur_mdl_tot_seq_len = utl_HSumValues(\%subseq_len_H);
      }
    }

    # run cmalign/glsearch
    @{$stk_file_HA{$mdl_name}} = ();
    if($cur_mdl_nalign > 0) { 
      cmalign_or_glsearch_wrapper(\%execs_H, $qsub_prefix, $qsub_suffix, \$blastn_db_sqfile,
                                  ($do_glsearch ? $glsearch_db_file : $cm_file), 
                                  $mdl_name, $cur_mdl_align_fa_file, $out_root, "", $cur_mdl_nalign,
                                  $cur_mdl_tot_seq_len, $progress_w, \@{$stk_file_HA{$mdl_name}}, 
                                  \@overflow_seq_A, \@overflow_mxsize_A, \%opt_HH, \%ofile_info_HH);
    }

    if($do_blastn_ali) {
      # join alignments of subsequences and update all variables
      my $start_secs = ofile_OutputProgressPrior(sprintf("Joining alignments from %s and blastn for model $mdl_name ($cur_mdl_nseq seq%s)",
                                                         ($do_glsearch ? "glsearch" : "cmalign"),
                                                         ($cur_mdl_nseq > 1) ? "s" : ""), $progress_w, $FH_HR->{"log"}, *STDOUT);
      
      # first, replace any overflow info we have on subseqs to be for full seqs and remove them from the list of seqs to align
      my @join_seq_name_A = (); # array of full seqs we'll try to join alignments for, this is all seqs except those with overflows
      if(scalar(@overflow_seq_A) == 0) { 
        @join_seq_name_A = @{$mdl_seq_name_HA{$mdl_name}};
      }
      else { # at least one overflow
        my @full_overflow_seq_A    = ();  
        my @full_overflow_mxsize_A = ();
        update_overflow_info_for_joined_alignments(\@overflow_seq_A, \@overflow_mxsize_A, \%subseq2seq_H, \@full_overflow_seq_A, \@full_overflow_mxsize_A);
        @overflow_seq_A    = @full_overflow_seq_A;
        @overflow_mxsize_A = @full_overflow_mxsize_A;

        # fill @join_seq_name_A with only seqs that do not have an overflow
        my %full_seq_overflow_H = ();
        my $full_seq;
        foreach $full_seq (@full_overflow_seq_A) { 
          $full_seq_overflow_H{$full_seq} = 1; 
        }
        @join_seq_name_A = ();
        foreach $full_seq (@{$mdl_seq_name_HA{$mdl_name}}) { 
          if(! defined $full_seq_overflow_H{$full_seq}) { 
            push(@join_seq_name_A, $full_seq);
          }
        }
      }

      my @joined_stk_file_A = ();   # array of joined stk files created by join_alignments_and_add_unjoinbl_alerts()
      my @unjoinbl_seq_name_A = (); # array of seqs with unjoinbl alerts
      if(scalar(@join_seq_name_A > 0)) { 
        join_alignments_and_add_unjoinbl_alerts($sqfile_for_analysis_R, \$blastn_db_sqfile, \%execs_H, 
                                                $do_glsearch, $cm_file,
                                                \@join_seq_name_A, \%seq_len_H, 
                                                \@mdl_info_AH, $mdl_idx, \%sda_mdl_H, \%sda_seq_H, 
                                                \%seq2subseq_HA, \%subseq_len_H, \@{$stk_file_HA{$mdl_name}}, 
                                                \@joined_stk_file_A, \%sda_output_HH,
                                                \%alt_seq_instances_HH, \%alt_info_HH,
                                                \@unjoinbl_seq_name_A, $out_root, \%opt_HH, \%ofile_info_HH);
      }
      push(@to_remove_A, (@{$stk_file_HA{$mdl_name}}));
      if(defined $ofile_info_HH{"fullpath"}{($mdl_name . ".cseq.fa")}) { push(@to_remove_A, $ofile_info_HH{"fullpath"}{($mdl_name . ".cseq.fa")}); }
      ofile_OutputProgressComplete($start_secs, undef, $FH_HR->{"log"}, *STDOUT);
      @{$stk_file_HA{$mdl_name}} = @joined_stk_file_A;

      # check for unjoinbl alerts, if we have any re-align the full seqs
      my $cur_unjoinbl_nseq = scalar(@unjoinbl_seq_name_A);
      if($cur_unjoinbl_nseq > 0) {
        # at least one sequence had 'unjoinbl' alert, align the full seqs
        # create fasta file
        my $unjoinbl_mdl_fa_file = $out_root . "." . $mdl_name . ".uj.a.fa";
        $$sqfile_for_analysis_R->fetch_seqs_given_names(\@unjoinbl_seq_name_A, 60, $unjoinbl_mdl_fa_file);
        ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . ".uj.a.fa", $unjoinbl_mdl_fa_file, 0, $do_keep, sprintf("%sinput seqs that match best to model $mdl_name with unjoinbl alerts", ($do_replace_ns) ? "replaced " : ""));
        $cur_mdl_tot_seq_len = utl_HSumValuesSubset(\%seq_len_H, \@unjoinbl_seq_name_A);
        cmalign_or_glsearch_wrapper(\%execs_H, $qsub_prefix, $qsub_suffix, \$blastn_db_sqfile,
                                    ($do_glsearch ? $glsearch_db_file : $cm_file), 
                                    $mdl_name, $unjoinbl_mdl_fa_file, $out_root, "uj.", $cur_unjoinbl_nseq,
                                    $cur_mdl_tot_seq_len, $progress_w, \@{$stk_file_HA{$mdl_name}}, 
                                    \@overflow_seq_A, \@overflow_mxsize_A, \%opt_HH, \%ofile_info_HH);
        # append insert file we just made to larger join insert file (if we created it (we may not have if all seqs had overflow error))
        my $cur_ifile = sprintf("%s.%s.uj.align.ifile", $out_root, $mdl_name);
        if(-s $cur_ifile) { 
          my $concat_cmd = sprintf("cat $cur_ifile >> %s.%s.jalign.ifile", $out_root, $mdl_name);
          utl_RunCommand($concat_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
        }
      }
    }

    # add unexdivg errors: sequences that were too divergent to align (cmalign was unable to align with a DP matrix of allowable size)
    $mdl_unexdivg_H{$mdl_name} = scalar(@overflow_seq_A);
    if($mdl_unexdivg_H{$mdl_name} > 0) { 
      alert_add_unexdivg(\@overflow_seq_A, \@overflow_mxsize_A, \%alt_seq_instances_HH, \%alt_info_HH, \%opt_HH, \%ofile_info_HH);
      @{$mdl_unexdivg_HA{$mdl_name}} = @overflow_seq_A; # make copy to hash of arrays for processing later
    }
  }
}

######################################
# Parse cmalign or glsearch alignments
######################################
$start_secs = ofile_OutputProgressPrior("Determining annotation", $progress_w, $log_FH, *STDOUT);

for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
  $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
  my $mdl_tt   = (defined $mdl_info_AH[$mdl_idx]{"transl_table"}) ? $mdl_info_AH[$mdl_idx]{"transl_table"} : 1; # default to standard genetic code
  if(defined $mdl_seq_name_HA{$mdl_name}) { 
    my $mdl_nseq = scalar(@{$mdl_seq_name_HA{$mdl_name}});
    initialize_ftr_or_sgm_results_for_model(\@{$mdl_seq_name_HA{$mdl_name}}, \@{$ftr_info_HAH{$mdl_name}}, \%{$ftr_results_HHAH{$mdl_name}}, $FH_HR);
    initialize_ftr_or_sgm_results_for_model(\@{$mdl_seq_name_HA{$mdl_name}}, \@{$sgm_info_HAH{$mdl_name}}, \%{$sgm_results_HHAH{$mdl_name}}, $FH_HR);
    my %seq_inserts_HH = ();
    my $align_stdout_file = $out_root . "." . $mdl_name . ".align.stdout";
    my $align_ifile_file  = sprintf("%s.%s.align.ifile", $out_root, $mdl_name);
    if($do_blastn_ali) { # use a different ifile
      push(@to_remove_A, $align_ifile_file);
      $align_ifile_file  = sprintf("%s.%s.jalign.ifile", $out_root, $mdl_name);
    }

    # pre-calculate some info per-feature so we don't need to do it per-sequence
    my $nftr = scalar(@{$ftr_info_HAH{$mdl_name}});
    my @ftr_fileroot_A = (); # for naming output files for each feature
    my @ftr_outroot_A  = (); # for describing output files for each feature
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      $ftr_fileroot_A[$ftr_idx] = vdr_FeatureTypeAndTypeIndexString(\@{$ftr_info_HAH{$mdl_name}}, $ftr_idx, ".");
      $ftr_outroot_A[$ftr_idx]  = vdr_FeatureTypeAndTypeIndexString(\@{$ftr_info_HAH{$mdl_name}}, $ftr_idx, "#");
    }

    # parse the cmalign --ifile file
    if($mdl_nseq > $mdl_unexdivg_H{$mdl_name}) { # at least 1 sequence was aligned
      vdr_CmalignParseInsertFile($align_ifile_file, \%seq_inserts_HH, undef, undef, undef, undef, \%{$ofile_info_HH{"FH"}});
      push(@to_remove_A, ($align_stdout_file, $align_ifile_file));
    }

    # parse the stk alignments
    for(my $a = 0; $a < scalar(@{$stk_file_HA{$mdl_name}}); $a++) { 
      if(-s $stk_file_HA{$mdl_name}[$a]) { # skip empty alignments, which may exist if all seqs were not alignable
        parse_stk_and_add_alignment_cds_and_mp_alerts($stk_file_HA{$mdl_name}[$a], \$in_sqfile, $mdl_tt,
                                                      \%seq_len_H, \%seq_inserts_HH, \@{$sgm_info_HAH{$mdl_name}},
                                                      \@{$ftr_info_HAH{$mdl_name}}, \%alt_info_HH, \%stg_results_HHH,
                                                      \%{$sgm_results_HHAH{$mdl_name}}, \%{$ftr_results_HHAH{$mdl_name}}, 
                                                      \%alt_seq_instances_HH, \%alt_ftr_instances_HHH, \%dcr_output_HAH,
                                                      $mdl_name, \@ftr_fileroot_A, \@ftr_outroot_A, 
                                                      $$sqfile_for_cds_mp_alerts_R, $$sqfile_for_output_fastas_R, $$sqfile_for_pv_R,
                                                      $do_separate_cds_fa_files_for_protein_validation, \@to_remove_A,
                                                      $out_root, \%opt_HH, \%ofile_info_HH);
      }
      push(@to_remove_A, ($stk_file_HA{$mdl_name}[$a]));
    }

    # Create option-defined output alignments, if any. 
    if(opt_Get("--keep", \%opt_HH) || opt_Get("--out_stk", \%opt_HH) || opt_Get("--out_afa", \%opt_HH) || opt_Get("--out_rpstk", \%opt_HH) || opt_Get("--out_rpafa", \%opt_HH)) { 
      if(scalar(@{$stk_file_HA{$mdl_name}}) > 0) { 
        output_alignments(\%execs_H, \$in_sqfile, \@{$stk_file_HA{$mdl_name}}, $mdl_name, \%rpn_output_HH, $out_root, \@to_remove_A, \%opt_HH, \%ofile_info_HH);
      }
    }
  }
}

# For any seqs with unexdivg alerts, add low similarity alerts, these
# were not aligned so we don't have alignment information (@ua2rf_A)
for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
  $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
  if(defined $mdl_unexdivg_HA{$mdl_name}) { 
    foreach my $seq_name (@{$mdl_unexdivg_HA{$mdl_name}}) { 
      add_low_similarity_alerts_for_one_sequence($seq_name, \%seq_len_H, undef,
                                                 \@{$ftr_info_HAH{$mdl_name}}, \@{$sgm_info_HAH{$mdl_name}}, \%alt_info_HH, 
                                                 \%stg_results_HHH, \%{$sgm_results_HHAH{$mdl_name}}, \%{$ftr_results_HHAH{$mdl_name}}, 
                                                 \%alt_seq_instances_HH, \%alt_ftr_instances_HHH, \%opt_HH, \%ofile_info_HH);
    }
  }
}

ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);


#########################################################################################
# Run BLASTX: all full length sequences and all fetched CDS features versus all proteins
#########################################################################################
if($do_pv_blastx) { 
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
    if(defined $mdl_seq_name_HA{$mdl_name}) { 
      my $ncds = vdr_FeatureInfoCountType(\@{$ftr_info_HAH{$mdl_name}}, "CDS"); 
      if($ncds > 0) { # only run blast for models with >= 1 CDS
        my $nseq = scalar(@{$mdl_seq_name_HA{$mdl_name}});
        # determine blastx db file to use
        # by default this will be the blastx db for current model but if --xsub is used
        # we may substitute another blastx db for this model
        my $blastx_db_mdl_idx = $mdl_idx;
        if((opt_IsUsed("--xsub", \%opt_HH)) && (defined $blastx_sub_H{$mdl_name})) { 
          # now find the index of $blastx_db_mdl_name in @mdl_info_AH
          $blastx_db_mdl_idx = 0; 
          while(($blastx_db_mdl_idx < $nmdl) && ($mdl_info_AH[$blastx_db_mdl_idx]{"name"} ne $blastx_sub_H{$mdl_name})) { 
            $blastx_db_mdl_idx++;
          }
          if($blastx_db_mdl_idx > $nmdl) { 
            ofile_FAIL("ERROR unable to find blastx model subsitution for $mdl_name (was looking for $blastx_sub_H{$mdl_name} in model info on second pass", 1, $FH_HR);
          }
        }
        my $blastx_db_file = $mdl_info_AH[$blastx_db_mdl_idx]{"blastdbpath"};
        if(! defined $blastx_db_file) { 
          ofile_FAIL("ERROR, path to BLAST DB is unknown for model $mdl_name", 1, $FH_HR);
        }

        $start_secs = ofile_OutputProgressPrior(sprintf("Validating proteins with blastx ($mdl_name: $nseq seq%s)", ($nseq > 1) ? "s" : ""), $progress_w, $log_FH, *STDOUT);
        run_blastx_and_summarize_output(\%execs_H, $out_root, \%{$mdl_info_AH[$mdl_idx]}, \@{$ftr_info_HAH{$mdl_name}}, $blastx_db_file, 
                                        $do_separate_cds_fa_files_for_protein_validation, 
                                        \%opt_HH, \%ofile_info_HH);
        push(@to_remove_A, 
             ($ofile_info_HH{"fullpath"}{$mdl_name . ".pv-blastx-fasta"},
              $ofile_info_HH{"fullpath"}{$mdl_name . ".blastx-out"},
              $ofile_info_HH{"fullpath"}{$mdl_name . ".blastx-summary"}));

        # if --xsub used and we have a substitute db for blastx, use that
        my $ftr_info_blastx_HR = ((opt_IsUsed("--xsub", \%opt_HH)) && (defined $blastx_sub_H{$mdl_name})) ? 
            \@{$ftr_info_HAH{$blastx_sub_H{$mdl_name}}} : \@{$ftr_info_HAH{$mdl_name}};
        parse_blastx_results($ofile_info_HH{"fullpath"}{($mdl_name . ".blastx-summary")}, \@{$mdl_seq_name_HA{$mdl_name}}, \%seq_len_H, 
                             $ftr_info_blastx_HR, \%{$ftr_results_HHAH{$mdl_name}}, \%opt_HH, \%ofile_info_HH);
        
        add_protein_validation_alerts($mdl_name, \@{$mdl_seq_name_HA{$mdl_name}}, \%seq_len_H, \@{$ftr_info_HAH{$mdl_name}}, \%alt_info_HH, 
                                      \%{$ftr_results_HHAH{$mdl_name}}, \%alt_ftr_instances_HHH, \%opt_HH, \%{$ofile_info_HH{"FH"}});
        ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
      }
    }
  }
} # end of 'if($do_pv_blastx)'
elsif(! $do_pv_hmmer) { 
  $start_secs = ofile_OutputProgressPrior("Skipping protein validation step (--pv_skip)", $progress_w, $log_FH, *STDOUT);
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
}

############################################################################################################
# Run hmmsearch: all full length sequences and all fetched CDS features versus best-matching protein profile
############################################################################################################
if($do_pv_hmmer) { 
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
    if(defined $mdl_seq_name_HA{$mdl_name}) { 
      my $ncds = vdr_FeatureInfoCountType(\@{$ftr_info_HAH{$mdl_name}}, "CDS"); 
      if($ncds > 0) { # only run blast for models with >= 1 CDS
        my $nseq = scalar(@{$mdl_seq_name_HA{$mdl_name}});
        $start_secs = ofile_OutputProgressPrior(sprintf("Validating proteins with hmmsearch ($mdl_name: $nseq seq%s)", ($nseq > 1) ? "s" : ""), $progress_w, $log_FH, *STDOUT);
        run_esl_translate_and_hmmsearch(\%execs_H, $out_root, \%{$mdl_info_AH[$mdl_idx]}, \@{$ftr_info_HAH{$mdl_name}}, 
                                        $do_separate_cds_fa_files_for_protein_validation, \%opt_HH, \%ofile_info_HH);
        parse_hmmer_domtblout($ofile_info_HH{"fullpath"}{($mdl_name . ".domtblout")}, 0, \@{$mdl_seq_name_HA{$mdl_name}}, \%seq_len_H, 
                                  \@{$ftr_info_HAH{$mdl_name}}, \%{$ftr_results_HHAH{$mdl_name}}, \%opt_HH, \%ofile_info_HH);
        add_protein_validation_alerts($mdl_name, \@{$mdl_seq_name_HA{$mdl_name}}, \%seq_len_H, \@{$ftr_info_HAH{$mdl_name}}, \%alt_info_HH, 
                                      \%{$ftr_results_HHAH{$mdl_name}}, \%alt_ftr_instances_HHH, \%opt_HH, \%{$ofile_info_HH{"FH"}});
        ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
        push(@to_remove_A, 
             ($ofile_info_HH{"fullpath"}{$mdl_name . ".pv.hmmer.esl_translate.aa.fa"}, 
              $ofile_info_HH{"fullpath"}{$mdl_name . ".pv.hmmer.fa"}, 
              $ofile_info_HH{"fullpath"}{$mdl_name . ".hmmsearch"},
              $ofile_info_HH{"fullpath"}{$mdl_name . ".hmmlist"},
              $ofile_info_HH{"fullpath"}{$mdl_name . ".domtblout"}));
      }                
    }
  }
} # end of 'if($do_pv_hmmer)'

##############################################################
# Add alerts to children features that have parents with 
# fatal alerts for specific feature combinations:
# (currently only one such parent/child type relationship
#  but could be expanded):
#      parent_ftr_type  child_ftr_type child_alert
#      ---------------  -------------- -----------
#      CDS              mat_peptide    peptrans
##############################################################
for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
  $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
  if(defined $mdl_seq_name_HA{$mdl_name}) { 
    alert_add_parent_based(\@{$mdl_seq_name_HA{$mdl_name}}, \@{$ftr_info_HAH{$mdl_name}}, \%alt_info_HH, \%{$ftr_results_HHAH{$mdl_name}}, 
                           \%alt_ftr_instances_HHH, "CDS", "mat_peptide", "peptrans", "VADRNULL", \%opt_HH, \%{$ofile_info_HH{"FH"}});
  }
}

##############################################################
# Potentially remove annotation and alerts for some features
# if we have 'alternative_ftr_set' information in the mdl info
##############################################################
for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
  $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
#  if((0) && (defined $mdl_seq_name_HA{$mdl_name})) { 
 if(defined $mdl_seq_name_HA{$mdl_name}) { 
    if(vdr_FeatureInfoValidateAlternativeFeatureSet(\@{$ftr_info_HAH{$mdl_name}}, $FH_HR)) { 
      # we'll enter this 'if' if any features have non-empty alternative_ftr_set
      
      # get children info for all features, we'll need this when picking features
      my @i_am_child_A = ();
      my @children_AA  = ();
      my $nchildren = vdr_FeatureInfoChildrenArrayOfArrays(\@{$ftr_info_HAH{$mdl_name}}, undef, \@i_am_child_A, \@children_AA, $FH_HR);
                                                           
      # first pick features from sets that are not composed of any children
      # this will remove features in alternative_ftr_sets that are not picked *and* their children
      pick_features_from_all_alternatives(\@{$mdl_seq_name_HA{$mdl_name}}, \@{$ftr_info_HAH{$mdl_name}}, \%alt_info_HH, 
                                          \%{$ftr_results_HHAH{$mdl_name}}, \%alt_ftr_instances_HHH, 
                                          0, \@i_am_child_A, \@children_AA, 
                                          \%opt_HH, \%{$ofile_info_HH{"FH"}});
      # now pick features from sets that are not composed of any children (that we have left)
      if($nchildren > 0) { 
        pick_features_from_all_alternatives(\@{$mdl_seq_name_HA{$mdl_name}}, \@{$ftr_info_HAH{$mdl_name}}, \%alt_info_HH, 
                                            \%{$ftr_results_HHAH{$mdl_name}}, \%alt_ftr_instances_HHH, 
                                            1, \@i_am_child_A, \@children_AA,
                                            \%opt_HH, \%{$ofile_info_HH{"FH"}});
      }
    }
  }
}

################################
# Output annotations and alerts
################################

# Output feature table first, because we may add noftranc alerts for
# sequences for which zero of the annotated features are output to 
# the feature table file (because they are too short). 

######################
# feature table file #
######################

# open files for writing
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "pass_tbl",       $out_root . ".pass.tbl",       1, 1, "5 column feature table output for passing sequences");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "fail_tbl",       $out_root . ".fail.tbl",       1, 1, "5 column feature table output for failing sequences");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "pass_list",      $out_root . ".pass.list",      1, 1, "list of passing sequences");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "fail_list",      $out_root . ".fail.list",      1, 1, "list of failing sequences");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "alerts_list",    $out_root . ".alt.list",       1, 1, "list of alerts in the feature tables");

$start_secs = ofile_OutputProgressPrior("Generating feature table output", $progress_w, $log_FH, *STDOUT);
my $npass = output_feature_table(\%mdl_cls_ct_H, \@seq_name_A, \%ftr_info_HAH, \%sgm_info_HAH, \%alt_info_HH, 
                                 \%stg_results_HHH, \%ftr_results_HHAH, \%sgm_results_HHAH, \%alt_seq_instances_HH,
                                 \%alt_ftr_instances_HHH, 
                                 ((opt_IsUsed("--msub", \%opt_HH)) ? \%mdl_sub_H : undef),
                                 \$in_sqfile, $out_root, \%opt_HH, \%ofile_info_HH);
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

########################
# tabular output files #
########################

# open files for writing
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "ant",      $out_root . ".sqa", 1, 1, "per-sequence tabular annotation summary file");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "cls",      $out_root . ".sqc", 1, 1, "per-sequence tabular classification summary file");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "ftr",      $out_root . ".ftr", 1, 1, "per-feature tabular summary file");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "sgm",      $out_root . ".sgm", 1, 1, "per-model-segment tabular summary file");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "mdl",      $out_root . ".mdl", 1, 1, "per-model tabular summary file");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "alt",      $out_root . ".alt", 1, 1, "per-alert tabular summary file");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "alc",      $out_root . ".alc", 1, 1, "alert count tabular summary file");
ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "dcr",      $out_root . ".dcr", 1, 1, "alignment doctoring tabular summary file");
if($do_blastn_ali) {
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "sda",    $out_root . ".sda", 1, 1, "seed alignment summary file (-s)");
}
if($do_replace_ns) { 
  ofile_OpenAndAddFileToOutputInfo(\%ofile_info_HH, "rpn",    $out_root . ".rpn", 1, 1, "replaced stretches of Ns summary file (-r)");
}

$start_secs = ofile_OutputProgressPrior("Generating tabular output", $progress_w, $log_FH, *STDOUT);
my ($zero_cls, $zero_alt) = output_tabular(\@mdl_info_AH, \%mdl_cls_ct_H, \%mdl_ant_ct_H, \@seq_name_A, \%seq_len_H, 
                                           \%ftr_info_HAH, \%sgm_info_HAH, \%alt_info_HH, \%cls_output_HH, \%ftr_results_HHAH, \%sgm_results_HHAH, 
                                           \%alt_seq_instances_HH, \%alt_ftr_instances_HHH, \%dcr_output_HAH,
                                           ($do_blastn_ali) ? \%sda_output_HH : undef,
                                           ($do_replace_ns) ? \%rpn_output_HH : undef,
                                           ((opt_IsUsed("--msub", \%opt_HH)) ? \%mdl_sub_H : undef),
                                           \%opt_HH, \%ofile_info_HH);
ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);


################################
# output optional output files #
################################
if(exists $ofile_info_HH{"FH"}{"ftrinfo"}) { 
  utl_HAHDump($ofile_info_HH{"desc"}{"ftrinfo"}, \%ftr_info_HAH, $ofile_info_HH{"FH"}{"ftrinfo"});
}
if(exists $ofile_info_HH{"FH"}{"sgminfo"}) { 
  utl_HAHDump($ofile_info_HH{"desc"}{"sgminfo"}, \%sgm_info_HAH, $ofile_info_HH{"FH"}{"sgminfo"});
}
if(exists $ofile_info_HH{"FH"}{"altinfo"}) { 
  vdr_AlertInfoDump(\%alt_info_HH, $ofile_info_HH{"FH"}{"altinfo"});
}
if(exists $ofile_info_HH{"FH"}{"stgresults"}) { 
  utl_HHHDump($ofile_info_HH{"desc"}{"stgresults"}, \%stg_results_HHH, $ofile_info_HH{"FH"}{"stgresults"});
}
if(exists $ofile_info_HH{"FH"}{"ftrresults"}) { 
  utl_HHAHDump($ofile_info_HH{"desc"}{"ftrresults"}, \%ftr_results_HHAH, $ofile_info_HH{"FH"}{"ftrresults"});
}
if(exists $ofile_info_HH{"FH"}{"sgmresults"}) { 
  utl_HHAHDump($ofile_info_HH{"desc"}{"sgmresults"}, \%sgm_results_HHAH, $ofile_info_HH{"FH"}{"sgmresults"});
}
if(exists $ofile_info_HH{"FH"}{"altseqinstances"}) { 
  utl_HHDump($ofile_info_HH{"desc"}{"altseqinstances"}, \%alt_seq_instances_HH, $ofile_info_HH{"FH"}{"altseqinstances"});
}
if(exists $ofile_info_HH{"FH"}{"altftrinstances"}) { 
  utl_HHHDump($ofile_info_HH{"desc"}{"altftrinstances"}, \%alt_ftr_instances_HHH, $ofile_info_HH{"FH"}{"altftrinstances"});
}
if(exists $ofile_info_HH{"FH"}{"clsoutput"}) { 
  utl_HHDump($ofile_info_HH{"desc"}{"clsoutput"}, \%cls_output_HH, $ofile_info_HH{"FH"}{"clsoutput"});
}
if(exists $ofile_info_HH{"FH"}{"rpnoutput"}) { 
  utl_HHDump($ofile_info_HH{"desc"}{"rpnoutput"}, \%rpn_output_HH, $ofile_info_HH{"FH"}{"rpnoutput"});
}
if(exists $ofile_info_HH{"FH"}{"sdaoutput"}) { 
  utl_HHDump($ofile_info_HH{"desc"}{"sdaoutput"}, \%sda_output_HH, $ofile_info_HH{"FH"}{"sdaoutput"});
}

############
# Conclude #
############
output_mdl_and_alc_files_and_remove_temp_files($zero_alt, \@to_remove_A, \%opt_HH, \%ofile_info_HH);

$total_seconds += ofile_SecondsSinceEpoch();
ofile_OutputConclusionAndCloseFilesOk($total_seconds, $dir, \%ofile_info_HH);

exit 0;

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
#################################################################
#
#################################################################
#
# Subroutines related to classification and coverage determination:
# classification_stage()
# coverage_determination_stage()
# cmsearch_wrapper
# cmsearch_run
# cmsearch_parse_sorted_tblout
# cmsearch_store_hit
# add_classification_errors
# populate_per_model_data_structures_given_classification_results
#
# Subroutines related to cmalign and alignment:
# cmalign_or_glsearch_wrapper
# cmalign_or_glsearch_wrapper_helper
# cmalign_or_glsearch_run
# parse_stk_and_add_alignment_cds_and_mp_alerts 
# add_frameshift_alerts_for_one_sequence
# cmalign_store_overflow
# fetch_features_and_add_cds_and_mp_alerts_for_one_sequence
# sqstring_check_start
# sqstring_find_stops 
# add_low_similarity_alerts_for_one_sequence
# frameshift_determine_span
#
# Subroutines related to protein validation (blastx or hmmer):
# make_protein_validation_fasta_file()
# add_protein_validation_alerts
# run_blastx_and_summarize_output
# parse_blastx_results 
# run_esl_translate_and_hmmsearch()
# parse_hmmer_domtblout()
# convert_esl_translate_to_blastx_frame()
# get_hmm_list_for_model
# helper_protein_validation_breakdown_source
# helper_blastx_breakdown_max_indel_str
# helper_blastx_max_indel_token_to_alt_coords
# helper_protein_validation_db_seqname_to_ftr_idx 
# helper_protein_validation_check_overlap()
#
# Other subroutines related to alerts: 
# alert_list_option
# alert_pass_fail_options
# alert_sequence_instance_add 
# alert_feature_instance_add 
# alert_sequence_instance_fetch
# alert_feature_instance_fetch
# alert_add_unexdivg 
# alert_instances_check_prevents_annot
#
# Subroutines for creating tabular output:
# output_tabular
# helper_tabular_ftr_results_strand
# helper_tabular_ftr_results_trunc_string
# helper_tabular_sgm_results_trunc_string
# helper_tabular_replace_spaces
#
# Subroutines for creating feature table output:
# output_feature_table
# output_parent_child_relationships 
# helper_ftable_coords_from_nt_prediction 
# helper_ftable_coords_prot_only_prediction 
# helper_ftable_start_stop_strand_arrays_to_coords 
# helper_ftable_coords_to_out_str 
# helper_ftable_add_qualifier_from_ftr_info
# helper_ftable_add_qualifier_from_ftr_results
# helper_ftable_class_model_for_sequence
# helper_ftable_process_sequence_alerts
# helper_ftable_process_feature_alerts
#
# Other output-related subroutines:
# helper_output_sequence_alert_strings()
# helper_output_feature_alert_strings()
# output_alignments()
# msa_replace_sequences()
# 
# Subroutines related to -r:
# parse_cdt_tblout_file_and_replace_ns()
# 
# Miscellaneous subroutines:
# initialize_ftr_or_sgm_results()
# convert_pp_char_to_pp_avg ()
# group_subgroup_string_from_classification_results()
# get_accession_from_ncbi_seq_name() 
# check_for_tabular_ftr_feature_prediction()
# check_for_valid_ftbl_feature_prediction()
# check_if_sequence_passes()
# check_if_sequence_was_annotated()
# check_for_feature_alert_codes()
# helper_sort_hit_array()
# get_5p_most_sgm_idx_with_results()
# get_3p_most_sgm_idx_with_results()
# 
#################################################################
# 
# Subroutines related to classification and coverage determination:
# classification_stage()
# coverage_determination_stage()
# cmsearch_wrapper
# cmsearch_run
# cmsearch_parse_sorted_tblout
# cmsearch_store_hit
# add_classification_errors
# populate_per_model_data_structures_given_classification_results
#
#################################################################
# Subroutine:  classification_stage()
# Incept:      EPN, Mon Apr 13 17:07:28 2020
#
# Purpose:    Wrapper that does all the steps of the classification
#             stage.
#
# Arguments: 
#  $execs_HR:          REF to a hash with "blastx" and "parse_blastx.pl""
#  $stg_key:           stage key, "rpn.cls" or "std.cls"
#  $cm_file:           path to the CM file 
#  $blastn_db_file:    path to blastn db
#  $fa_file:           fasta file
#  $seq_len_HR:        REF to hash of sequence lengths, pre-filled
#  $qsub_prefix:       prefix for qsub calls
#  $qsub_suffix:       suffix for qsub calls
#  $mdl_info_AHR:      REF to array of hashes with model info     
#  $stg_results_HHHR:  REF to 3D hash of classification results, PRE-FILLED
#  $out_root:          root name for output file names
#  $progress_w:        width for outputProgressPrior output
#  $to_remove_AR:      REF to array of files to remove before exiting (unless --keep)
#  $opt_HHR:           REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If esl-alimerge fails.
#
################################################################# 
sub classification_stage { 
  my $sub_name = "classification_stage";
  my $nargs_exp = 15;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $stg_key, $cm_file, $blastn_db_file, $fa_file, $seq_len_HR, 
      $qsub_prefix, $qsub_suffix, $mdl_info_AHR, $stg_results_HHHR, 
      $out_root, $progress_w, $to_remove_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $ncpu = opt_Get("--cpu", $opt_HHR);
  if($ncpu == 0) { $ncpu = 1; }

  if(($stg_key ne "rpn.cls") && ($stg_key ne "std.cls")) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized stage key: $stg_key, should be rpn.cls or std.cls", 1, $FH_HR);
  }

  my $nseq = scalar(keys %{$seq_len_HR});
  # will we use blastn or cmsearch?
  my $do_blastn  = 0; 
  # default is to use blastn in rpn.cls, but don't if --r_prof used
  if((! opt_Get("--r_prof", $opt_HHR)) && ($stg_key eq "rpn.cls")) { $do_blastn = 1; }
  # default is to not use blastn in std.cls, but do if -s used
  if((  opt_Get("-s",       $opt_HHR)) && ($stg_key eq "std.cls")) { $do_blastn = 1; }

  if($do_blastn) { # -s: use blastn for classification
    run_blastn_and_summarize_output($execs_HR, $blastn_db_file, $fa_file, $out_root, $stg_key,
                                    $nseq, $progress_w, $opt_HHR, $ofile_info_HHR);
    parse_blastn_results($ofile_info_HHR->{"fullpath"}{"$stg_key.blastn.summary"}, $seq_len_HR, 
                         undef, undef, $out_root, $stg_key, $opt_HHR, $ofile_info_HHR);
    push(@{$to_remove_AR}, 
         $ofile_info_HHR->{"fullpath"}{"$stg_key.blastn.out"},
         $ofile_info_HHR->{"fullpath"}{"$stg_key.blastn.summary"},
         $ofile_info_HHR->{"fullpath"}{"$stg_key.blastn.pretblout"});
  }
  else { # default: use cmsearch for classification
    my $cmsearch_opts = " -T " . opt_Get("--minbit", $opt_HHR) . " --cpu $ncpu --trmF3 --noali --hmmonly"; 
    my $tot_len_nt  = utl_HSumValues($seq_len_HR);
    cmsearch_wrapper($execs_HR, $qsub_prefix, $qsub_suffix,
                               $cm_file, undef, $fa_file, $cmsearch_opts, 
                               $out_root, $stg_key, $nseq, $tot_len_nt, 
                               $progress_w, $opt_HHR, $ofile_info_HHR);
  } 
  # sort into a new file by score
  my $tblout_key  = $stg_key . ".tblout"; # set in cmsearch_wrapper() or above if $do_blastn_cls
  my $stdout_key  = $stg_key . ".stdout"; # set in cmsearch_wrapper()
  my $err_key     = $stg_key . ".err"; # set in cmsearch_wrapper()
  my $tblout_file = $ofile_info_HH{"fullpath"}{$tblout_key};
  my $sort_tblout_file = $tblout_file . ".sort";
  my $sort_tblout_key  = $tblout_key  . ".sort";
  utl_FileValidateExistsAndNonEmpty($tblout_file, "$stg_key stage tblout output", undef, 1, \%{$ofile_info_HH{"FH"}}); # '1' says: die if it doesn't exist or is empty

  my $sort_cmd = "grep -v ^\# $tblout_file | sed 's/  */ /g' | sort -k 1,1 -k 3,3rn > $sort_tblout_file"; 
  # the 'sed' call replaces multiple spaces with a single one, because sort is weird about multiple spaces sometimes
  utl_RunCommand($sort_cmd, opt_Get("-v", \%opt_HH), 0, $FH_HR);
  ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $sort_tblout_key, $sort_tblout_file, 0, $do_keep, "stage $stg_key sorted tblout file");
  push(@{$to_remove_AR}, 
       ($tblout_file, 
        $ofile_info_HH{"fullpath"}{$stdout_key},
        $ofile_info_HH{"fullpath"}{$err_key}, 
        $sort_tblout_file));
  
  # parse the sorted tblout file
  cmsearch_parse_sorted_tblout($sort_tblout_file, $stg_key,
                               $mdl_info_AHR, $stg_results_HHHR, $opt_HHR, $FH_HR);

  return;
}

#################################################################
# Subroutine:  coverage_determination_stage()
# Incept:      EPN, Tue Apr 14 06:53:31 2020
#
# Purpose:    Wrapper that does all the steps of the coverage
#             determination stage.
#
# Arguments: 
#  $execs_HR:          REF to a hash with "blastx" and "parse_blastx.pl""
#  $stg_key:           stage key, "rpn.cdt" or "std.cdt"
#  $cm_file:           path to the CM file 
#  $sqfile_R:          REF to Bio::Easel::SqFile object from main fasta file
#  $seq_name_AR:       REF to array of sequence names, pre-filled
#  $seq_len_HR:        REF to hash of sequence lengths, pre-filled
#  $qsub_prefix:       prefix for qsub calls
#  $qsub_suffix:       suffix for qsub calls
#  $mdl_info_AHR:      REF to array of hashes with model info     
#  $mdl_seq_name_HAR:  REF to hash of arrays of sequences per model, can be undef if stage_key is "std.cdt"
#  $mdl_cls_ct_HR:     REF to hash of counts of seqs assigned to each model, can be undef if stage_key is "rpn.cdt"
#  $stg_results_HHHR:  REF to 3D hash of classification results, PRE-FILLED
#  $mdl_sub_HR:        REF to hash with model substitution info, should be undef unless --msub used
#  $out_root:          root name for output file names
#  $progress_w:        width for outputProgressPrior output
#  $to_remove_AR:      REF to array of files to remove before exiting (unless --keep)
#  $opt_HHR:           REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If $stg_key is unrecognized
#
################################################################# 
sub coverage_determination_stage { 
  my $sub_name = "coverage_determination_stage";
  my $nargs_exp = 18;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $stg_key, $cm_file, $sqfile_R, $seq_name_AR, $seq_len_HR, 
      $qsub_prefix, $qsub_suffix, $mdl_info_AHR, $mdl_seq_name_HAR, $mdl_cls_ct_HR, $stg_results_HHHR, 
      $mdl_sub_HR, $out_root, $progress_w, $to_remove_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $ncpu = opt_Get("--cpu", $opt_HHR);
  if($ncpu == 0) { $ncpu = 1; }

  if(($stg_key ne "rpn.cdt") && ($stg_key ne "std.cdt")) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized stage key: $stg_key, should be rpn.cdt or std.cdt", 1, $FH_HR);
  }
  if(($stg_key eq "std.cdt") && (! defined $mdl_cls_ct_HR)) { 
    ofile_FAIL("ERROR in $sub_name, stage key: $stg_key, but mdl_cls_ct_HR is undef", 1, $FH_HR);
  }
  if(($stg_key eq "rpn.cdt") && (! defined $mdl_seq_name_HAR)) { 
    ofile_FAIL("ERROR in $sub_name, stage key: $stg_key, but mdl_seq_name_HAR is undef", 1, $FH_HR);
  }

  my $nmdl = scalar(@{$mdl_info_AHR});
  my $mdl_name;            # a model name
  my $seq_name;            # a sequence name
  my @tblout_key_A  = ();  # array of round 2 search tblout keys in %ofile_info_HH
  my @tblout_file_A = ();  # array of round 2 search tblout files 
  my $nseq          = scalar(@{$seq_name_AR});

  # will we use blastn or cmsearch?
  my $do_blastn  = 0; 
  # default is to use blastn in rpn.cdt, but don't if --r_prof used
  if((! opt_Get("--r_prof", $opt_HHR)) && ($stg_key eq "rpn.cdt")) { $do_blastn = 1; }
  # default is to not use blastn in std.cdt, but do if -s used
  if((  opt_Get("-s",       $opt_HHR)) && ($stg_key eq "std.cdt")) { $do_blastn = 1; }

  # fill per-model data structures based on classification reesults
  my %mdl_seq_name_HA  = ();  # key is model name $mdl_name, array is of seq names classified to model $mdl_name 
  my %mdl_seq_len_H    = ();  # key is model name $mdl_name, value is summed length of all seqs in @{$mdl_seq_HA{$mdl_name}
  my %seq2mdl_H        = ();  # key is sequence name $seq_name, value is $mdl_name of model this sequence is classified to
  my $local_mdl_seq_name_HAR = (defined $mdl_seq_name_HAR) ? $mdl_seq_name_HAR : \%mdl_seq_name_HA;  # to deal with fact that mdl_seq_name_HAR may be undef
  populate_per_model_data_structures_given_classification_results($seq_name_AR, $seq_len_HR, $stg_key, $stg_results_HHHR, undef, undef, undef, 
                                                                  $mdl_sub_HR, $local_mdl_seq_name_HAR, \%mdl_seq_len_H, $mdl_cls_ct_HR, \%seq2mdl_H, $FH_HR);
  
  # for each model, fetch the sequences classified to it 
  my $nmdl_cdt = 0; # number of models coverage determination stage is called for
  my @cls_mdl_name_A = (); # array of model names that have >=1 sequences classified to them
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name = $mdl_info_AH[$mdl_idx]{"name"};
    if(defined $local_mdl_seq_name_HAR->{$mdl_name}) { 
      my $mdl_fa_file = $out_root . "." . $mdl_name . ".fa";
      push(@{$to_remove_AR}, $mdl_fa_file);
      $$sqfile_R->fetch_seqs_given_names(\@{$local_mdl_seq_name_HAR->{$mdl_name}}, 60, $mdl_fa_file);
      push(@cls_mdl_name_A, $mdl_name);
      $nmdl_cdt++;
    }
  }
  
  # get the coverage determination results:
  #   if default mode: perform the per model round 2 searches to get sequence coverage
  #   if blastn mode:  parse our blastn results a second time to get model
  #                    specific tblout files to use instead of cmsearch tblout
  #                    files
  if($do_blastn) { 
    my $stg_desc = "";
    if($stg_key eq "rpn.cdt") { 
      $stg_desc = sprintf("Preprocessing for N replacement: coverage determination from blastn results ($nseq seq%s)", ($nseq > 1) ? "s" : "");
    }
    else { # stg_key eq "std.cdt"
      $stg_desc = sprintf("Determining sequence coverage from blastn results ($nseq seq%s)", ($nseq > 1) ? "s" : "");
    }
    $start_secs = ofile_OutputProgressPrior($stg_desc, $progress_w, $log_FH, *STDOUT);
    my $blastn_summary_key = ($stg_key eq "rpn.cdt") ? "rpn.cls.blastn.summary" : "std.cls.blastn.summary";
    parse_blastn_results($ofile_info_HHR->{"fullpath"}{$blastn_summary_key}, $seq_len_HR, 
                         \%seq2mdl_H, \@cls_mdl_name_A, $out_root, $stg_key, $opt_HHR, $ofile_info_HHR);
    # keep track of the tblout output files:
    foreach $mdl_name (@cls_mdl_name_A) { 
      my $tblout_key = "$stg_key.$mdl_name.tblout";
      push(@tblout_key_A,  $tblout_key);
      push(@tblout_file_A, $ofile_info_HH{"fullpath"}{$tblout_key});
      push(@to_remove_A, 
           ($ofile_info_HH{"fullpath"}{$tblout_key}, 
            $ofile_info_HH{"fullpath"}{"$stg_key.$mdl_name.indel"}));
    }
    ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
  else { # default, not (! $do_blastn) 
    my $cmsearch_opts = " -T " . opt_Get("--minbit", $opt_HHR) . " --cpu $ncpu --hmmonly "; # cmsearch options for round 2 searches to determine coverage

    if(! opt_Get("-v", \%opt_HH)) { $cmsearch_opts .= " --noali "; }
    foreach $mdl_name (@cls_mdl_name_A) { 
      my $mdl_fa_file = $out_root . "." . $mdl_name . ".fa";
      cmsearch_wrapper(\%execs_H, $qsub_prefix, $qsub_suffix,
                                 $cm_file, $mdl_name, $mdl_fa_file, $cmsearch_opts, 
                                 $out_root, $stg_key, scalar(@{$local_mdl_seq_name_HAR->{$mdl_name}}), 
                                 $mdl_seq_len_H{$mdl_name}, $progress_w, \%opt_HH, \%ofile_info_HH);
      my $tblout_key = "$stg_key.$mdl_name.tblout"; # set in cmsearch_wrapper()
      my $stdout_key = "$stg_key.$mdl_name.stdout"; # set in cmsearch_wrapper()
      my $err_key    = "$stg_key.$mdl_name.err";    # set in cmsearch_wrapper()
      push(@tblout_key_A,  $tblout_key);
      push(@tblout_file_A, $ofile_info_HH{"fullpath"}{$tblout_key});
      push(@to_remove_A, 
           ($ofile_info_HH{"fullpath"}{$tblout_key}, 
            $ofile_info_HH{"fullpath"}{$stdout_key}, 
            $ofile_info_HH{"fullpath"}{$err_key}));
    }
  }

  # sort the coverage determination search results, we concatenate all model's tblout files and sort them
  my $sort_tblout_key  = "$stg_key.tblout.sort";
  my $sort_tblout_file = $out_root . "." . $sort_tblout_key;
  if($nmdl_cdt > 0) { # only sort output if we ran coverage determination stage for at least one model
    my $sort_cmd = "cat " . join(" ", @tblout_file_A) . " | grep -v ^\# | sed 's/  */ /g' | sort -k 1,1 -k 15,15rn -k 16,16g > $sort_tblout_file"; 
    # the 'sed' call replaces multiple spaces with a single one, because sort is weird about multiple spaces sometimes
    utl_RunCommand($sort_cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
    ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $sort_tblout_key, $sort_tblout_file, 0, $do_keep, "stage $stg_key sorted tblout file");
    push(@{$to_remove_AR}, $sort_tblout_file);

    # parse cmsearch round 2 tblout data
    cmsearch_parse_sorted_tblout($sort_tblout_file, $stg_key,
                                 $mdl_info_AHR, $stg_results_HHHR, $opt_HHR, $FH_HR);
  }
  return;
}

#################################################################
# Subroutine:  cmsearch_wrapper()
# Incept:      EPN, Mon Mar 18 14:44:46 2019
#
# Purpose:     Run one or more cmsearch jobs on the farm
#              or locally, after possibly splitting up the input
#              sequence file with vdr_SplitFastaFile and 
#              then calling vdr_CmalignOrCmsearchWrapperHelper(). 
#
# Arguments: 
#  $execs_HR:        ref to executables with "esl-ssplit" and "cmsearch"
#                    defined as keys
#  $qsub_prefix:     qsub command prefix to use when submitting to farm, undef if running locally
#  $qsub_suffix:     qsub command suffix to use when submitting to farm, undef if running locally
#  $mdl_file:        name of model file to use
#  $mdl_name:        name of model to fetch from $mdl_file (undef to not fetch)
#  $seq_file:        name of sequence file with all sequences to run against
#  $opt_str:         option string for cmsearch run
#  $out_root:        string for naming output files
#  $stg_key:         stage key, "rpn.cls" or "std.cls" for classification (cmsearch),
#                    or "rpn.cdt" or "std.cdt" for coverage determination (cmsearch)
#  $nseq:            number of sequences in $seq_file
#  $tot_len_nt:      total length of all nucleotides in $seq_file
#  $progress_w:      width for outputProgressPrior output
#  $opt_HHR:         REF to 2D hash of option values, see top of sqp-opts.pm for description
#  $ofile_info_HHR:  REF to 2D hash of output file information
#
# Returns:     void
# 
# Dies: If an executable doesn't exist, or cmalign or nhmmscan or esl-ssplit
#       command fails if we're running locally
################################################################# 
sub cmsearch_wrapper { 
  my $sub_name = "cmsearch_wrapper";
  my $nargs_expected = 14;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $qsub_prefix, $qsub_suffix, 
      $mdl_file, $mdl_name, $seq_file, $opt_str, 
      $out_root, $stg_key, $nseq, $tot_len_nt, 
      $progress_w, $opt_HHR, $ofile_info_HHR) = @_;

  my $log_FH = $ofile_info_HHR->{"FH"}{"log"}; # for convenience
  my $do_parallel = opt_Get("-p",     $opt_HHR);
  my $do_keep     = opt_Get("--keep", $opt_HHR);
  my $ncpu = opt_Get("--cpu", $opt_HHR);
  if($ncpu == 0) { $ncpu = 1; }

  if(($stg_key ne "rpn.cls") && ($stg_key ne "rpn.cdt") && 
     ($stg_key ne "std.cls") && ($stg_key ne "std.cdt")) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized stage key: $stg_key, should be rpn.cls, rpn.cdt, std.cls, or std.cdt", 1, $FH_HR);
  }

  # set up output file names
  my @seq_file_A  = (); # [0..$nr-1]: name of sequence file for this run
  my @out_file_AH = (); # array of hashes ([0..$nr-1]) of output files for cmsearch runs
  my $nseq_files = 1; # number of sequence files/runs 

  # split up sequence file, if -p and we're going to do more than one job
  my $do_split = 0; # set to '1' if we run vdr_SplitFastaFile
  if($do_parallel) { 
    # -p used: we need to split up the sequence file, and submit a separate 
    # cmsearch job for each
    my $targ_nseqfiles = vdr_SplitNumSeqFiles($tot_len_nt, $opt_HHR);
    if($targ_nseqfiles > 1) { # we are going to split up the fasta file 
      $do_split = 1;
      $nseq_files = vdr_SplitFastaFile($execs_HR->{"esl-ssplit"}, $seq_file, $targ_nseqfiles, undef, $opt_HHR, $ofile_info_HHR);
      # vdr_SplitFastaFile will return the actual number of fasta files created, 
      # which can differ from the requested amount (which is $targ_nseqfiles) that we pass in. 
      for(my $i = 0; $i < $nseq_files; $i++) { 
        $seq_file_A[$i] = $seq_file . "." . ($i+1);
      }
    }
    else { # targ_nseqfiles is 1, no need to split
      $seq_file_A[0] = $seq_file;
    }
  }
  else { # -p not used
    $seq_file_A[0] = $seq_file;
  }
    
  # determine description of the stage we are about to do
  my $stg_desc = "";
  if($stg_key eq "rpn.cls") { 
    $stg_desc = ($do_parallel) ? 
        sprintf("Preprocessing for N replacement: cmsearch classification job farm submission ($nseq seq%s, $nseq_files job%s)", (($nseq > 1) ? "s" : ""), (($nseq_files > 1) ? "s" : "")) :
        sprintf("Preprocessing for N replacement: cmsearch classification ($nseq seq%s)", ($nseq > 1) ? "s" : "");
  }
  elsif($stg_key eq "std.cls") { 
    $stg_desc = ($do_parallel) ? 
        "Submitting $nseq_files cmsearch classification job(s) to the farm" : 
        sprintf("Classifying sequences ($nseq seq%s)", ($nseq > 1) ? "s" : "");
  }
  elsif($stg_key eq "rpn.cdt") {
    $stg_desc = ($do_parallel) ? 
        sprintf("Preprocessing for N replacement: cmsearch coverage determination job farm submission ($mdl_name: $nseq seq%s, $nseq_files job%s)", (($nseq > 1) ? "s" : ""), (($nseq_files > 1) ? "s" : "")) :
        sprintf("Preprocessing for N replacement: cmsearch coverage determination ($mdl_name: $nseq seq%s)", ($nseq > 1) ? "s" : "");
  }
  else { 
    $stg_desc = ($do_parallel) ? 
        sprintf("Submitting $nseq_files cmsearch coverage determination job(s) ($mdl_name: $nseq seq%s) to the farm", ($nseq > 1) ? "s" : "") :
        sprintf("Determining sequence coverage ($mdl_name: $nseq seq%s)", ($nseq > 1) ? "s" : "");
  }
  my $start_secs = ofile_OutputProgressPrior($stg_desc, $progress_w, $log_FH, *STDOUT);
  # run cmsearch
  my $out_key;
  my @out_keys_A = ("stdout", "err", "tblout", "sh");
  for(my $s = 0; $s < $nseq_files; $s++) { 
    %{$out_file_AH[$s]} = (); 
    foreach my $out_key (@out_keys_A) { 
      $out_file_AH[$s]{$out_key} = $out_root . "." . $stg_key . ".s" . $s . "." . $out_key;
    }
    cmsearch_run($execs_HR, $qsub_prefix, $qsub_suffix, $mdl_file, $mdl_name, 
                           $seq_file_A[$s], $opt_str, \%{$out_file_AH[$s]}, $opt_HHR, $ofile_info_HHR);   
  }

  if($do_parallel) { 
    ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
    # wait for the jobs to finish
    $start_secs = ofile_OutputProgressPrior(sprintf("Waiting a maximum of %d minutes for all farm jobs to finish", opt_Get("--wait", $opt_HHR)), 
                                            $progress_w, $log_FH, *STDOUT);
    my $njobs_finished = vdr_WaitForFarmJobsToFinish(0, # we're not running cmalign
                                                     (opt_Get("--errcheck", $opt_HHR)), 
                                                     "tblout", 15, 15, \@out_file_AH, undef, undef, "[ok]", $opt_HHR, $ofile_info_HHR->{"FH"});
    if($njobs_finished != $nseq_files) { 
      ofile_FAIL(sprintf("ERROR in $sub_name only $njobs_finished of the $nseq_files are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", $opt_HHR)), 1, $ofile_info_HHR->{"FH"});
    }
    ofile_OutputString($log_FH, 1, "# "); # necessary because waitForFarmJobsToFinish() creates lines that summarize wait time and so we need a '#' before 'done' printed by ofile_OutputProgressComplete()
  }

  # concatenate files into one
  foreach $out_key (@out_keys_A) { 
    if(($do_parallel) || (($out_key ne "err") && ($out_key ne "sh"))) { # .err and .sh files don't exist if (! $do_parallel)
      my $concat_key  = sprintf("%s.%s%s", $stg_key, (defined $mdl_name) ? $mdl_name . "." : "", $out_key);                                
      my $concat_file = $out_root . "." . $concat_key;
      my @concat_A = ();
      utl_ArrayOfHashesToArray(\@out_file_AH, \@concat_A, $out_key);
      if(scalar(@concat_A) > 0) { 
        utl_ConcatenateListOfFiles(\@concat_A, $concat_file, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
        # utl_ConcatenateListOfFiles() removes individual files unless --keep enabled
        ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $concat_key, $concat_file, 0, $do_keep, sprintf("stage $stg_key $out_key file%s", (defined $mdl_name) ? "for model $mdl_name" : ""));
      }
    }
  }

  # remove sequence files if we created any
  if(($do_split) && (! opt_Get("--keep", $opt_HHR))) { 
    utl_FileRemoveList(\@seq_file_A, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
  }

  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  return;
}

#################################################################
# Subroutine:  cmsearch_run()
# Incept:      EPN, Wed Feb  6 12:38:11 2019
#
# Purpose:     Run Infernal's cmsearch executable using $mdl_file
#              as the model file on sequence file $seq_file, either
#              locally or on the farm.
#
# Arguments: 
#  $execs_HR:         hash with paths to cmsearch and cmfetch
#  $qsub_prefix:      qsub command prefix to use when submitting to farm, undef if running locally
#  $qsub_suffix:      qsub command suffix to use when submitting to farm, undef if running locally
#  $mdl_file:         path to the CM file
#  $mdl_name:         name of model to fetch from $mdl_file (undef to not fetch)
#  $seq_file:         path to the sequence file
#  $opt_str:          option string for cmsearch runs
#  $out_file_HR:      ref to hash of output files to create
#                     required keys: "stdout", "tblout", "err", "sh"
#  $opt_HHR:          REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:   REF to 2D hash of output file information
# 
# Returns:     void
# 
################################################################# 
sub cmsearch_run {
  my $sub_name = "cmsearch_run()";
  my $nargs_expected = 10;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($execs_HR, $qsub_prefix, $qsub_suffix, 
      $mdl_file, $mdl_name, $seq_file, $opt_str, 
      $out_file_HR, $opt_HHR, $ofile_info_HHR) = @_;
  
  my $FH_HR       = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $do_parallel = opt_Get("-p", $opt_HHR) ? 1 : 0;

  my $stdout_file = $out_file_HR->{"stdout"};
  my $tblout_file = $out_file_HR->{"tblout"};
  my $err_file    = $out_file_HR->{"err"};
  my $sh_file     = $out_file_HR->{"sh"};
  if(! defined $stdout_file) { ofile_FAIL("ERROR in $sub_name, stdout output file name is undefined", 1, $FH_HR); }
  if(! defined $tblout_file) { ofile_FAIL("ERROR in $sub_name, tblout output file name is undefined", 1, $FH_HR); }
  if(! defined $err_file)    { ofile_FAIL("ERROR in $sub_name, err    output file name is undefined", 1, $FH_HR); }
  if(-e $stdout_file) { unlink $stdout_file; }
  if(-e $tblout_file) { unlink $tblout_file; }
  if(-e $err_file)    { unlink $err_file; }
  
  utl_FileValidateExistsAndNonEmpty($mdl_file, "CM file", $sub_name, 1, $FH_HR); 
  utl_FileValidateExistsAndNonEmpty($seq_file, "sequence file", $sub_name, 1, $FH_HR);

  if(! (defined $opt_str)) { $opt_str = ""; }
  $opt_str .= " --tblout $tblout_file"; 

  my $cmd = undef;
  if(defined $mdl_name) { 
    $cmd = $execs_HR->{"cmfetch"} . " $mdl_file $mdl_name | " . $execs_HR->{"cmsearch"} . " $opt_str - $seq_file > $stdout_file";
  }
  else { 
    if(opt_IsUsed("--mlist", $opt_HHR)) { 
      # opt_Get("--mlist", $opt_HHR) includes the subset of models we want to use
      $cmd = $execs_HR->{"cmfetch"} . " -f $mdl_file " . opt_Get("--mlist", $opt_HHR) . " | " . $execs_HR->{"cmsearch"} . " $opt_str - $seq_file > $stdout_file";
    }
    else { 
      $cmd = $execs_HR->{"cmsearch"} . " $opt_str $mdl_file $seq_file > $stdout_file";
    }
  }
  if($do_parallel) { 
    my $job_name = "J" . utl_RemoveDirPath($seq_file);
    my $nsecs  = opt_Get("--wait", $opt_HHR) * 60.;
    my $mem_gb = (opt_Get("--mxsize", $opt_HHR) / 1000.); # use --mxsize * 1000 (8 Gb by default)
    if($mem_gb < 16.) { $mem_gb = 16.; } # set minimum of 16 Gb
    #vdr_SubmitJob($cmd, $qsub_prefix, $qsub_suffix, $job_name, $err_file, $mem_gb, $nsecs, $opt_HHR, $ofile_info_HHR);
    vdr_SubmitJobAsScript($cmd, $qsub_prefix, $qsub_suffix, $job_name, $sh_file, $err_file, $mem_gb, $nsecs, $opt_HHR, $ofile_info_HHR);
  }
  else { 
    utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);
  }
  return; 
}

#################################################################
# Subroutine:  cmsearch_parse_sorted_tblout()
# Incept:      EPN, Wed Mar 20 13:30:16 2019
#
# Purpose:     Parse a sorted cmsearch tblout output file and 
#              store results in %{$results_HHHR}. This is 
#              done differently depending on if $round is 1 or 2.
#
# Arguments: 
#  $tblout_file:   name of sorted tblout file to parse
#  $stg_key:       stage key, "rpn.cls" or "std.cls" for classification (cmsearch),
#                  or "rpn.cdt" or "std.cdt" for coverage determination (cmsearch)
#  $mdl_info_AHR:  ref to model info array of hashes
#  $results_HHHR:  ref to results 3D hash
#  $opt_HHR:       ref to options 2D hash
#  $FH_HR:         ref to file handle hash
# 
# Returns: void
# 
# Dies: if there's a problem parsing the tblout file
#      
################################################################# 
sub cmsearch_parse_sorted_tblout {
  my $sub_name = "cmsearch_parse_sorted_tblout()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($tblout_file, $stg_key, $mdl_info_AHR, $results_HHHR, $opt_HHR, $FH_HR) = @_;

  my $opt_mlist_used = opt_IsUsed("--mlist", $opt_HHR);

  if(($stg_key ne "rpn.cls") && ($stg_key ne "rpn.cdt") && 
     ($stg_key ne "std.cls") && ($stg_key ne "std.cdt")) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized stage key: $stg_key, should be rpn.cls, rpn.cdt, std.cls, or std.cdt", 1, $FH_HR);
  }

  # determine if we have an expected group and subgroup
  # (--subgroup requires --group)
  # and if so, fill hashes of model groups and subgroups, 
  # (this data is in @{$mdl_info_HAR} already, we fill 
  # these hashes only for convenience)
  my $exp_group    = opt_Get("--group", \%opt_HH);
  my $exp_subgroup = opt_Get("--subgroup", \%opt_HH);
  # fill hashes of model groups and subgroups, for convenience
  my %mdl_group_H    = ();
  my %mdl_subgroup_H = ();
  my $nmdl = scalar(@{$mdl_info_AHR});
  if(($stg_key eq "rpn.cls") || ($stg_key eq "std.cls")) { 
    for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
      my $mdl_name = $mdl_info_AHR->[$mdl_idx]{"name"};
      if(defined $mdl_info_AHR->[$mdl_idx]{"group"}) { 
        $mdl_group_H{$mdl_name} = $mdl_info_AHR->[$mdl_idx]{"group"}; 
      }
      if(defined $mdl_info_AHR->[$mdl_idx]{"subgroup"}) { 
        $mdl_subgroup_H{$mdl_name} = $mdl_info_AHR->[$mdl_idx]{"subgroup"}; 
      }
    }
  }

  my $seq;      # sequence name
  my $model;    # model name
  my $score;    # bit score
  my $m_from;   # model from position
  my $m_to;     # model to position
  my $s_from;   # sequence from position
  my $s_to;     # sequence to position
  my $strand;   # strand of hit
  my $bias;     # bias score
  my $key2;     # key in 2nd dim of %results
  my $group;    # group for current sequence, can be undef
  my $subgroup; # subgroup for current sequence, can be undef
  open(IN, $tblout_file) || ofile_FileOpenFailure($tblout_file, $sub_name, 1, "reading", $FH_HR);

  while(my $line = <IN>) { 
    my $HR = undef; # ref to 3rd dimension hash in %results_HHH we will update for this hit
    if($line !~ m/^\#/) { 
      chomp $line;
      $seq    = undef;
      $model  = undef;
      $score  = undef;
      $m_from = undef;
      $m_to   = undef;
      $s_from = undef;
      $s_to   = undef;
      $strand = undef;
      my @el_A = split(/\s+/, $line);
      if(($stg_key eq "rpn.cls") || ($stg_key eq "std.cls")) { # cmsearch --trmF3 tblout 
        # cmsearch output
        #sequence                      modelname  score  start    end strand bounds ovp      seqlen
        ##---------------------------  ---------  ----------------- ------ ------ ------ ------ ------ --- -----------
        #gi|1215708385|gb|KY594653.1|   NC_039477 275.8      1    301      +     []  *          301
        if(scalar(@el_A) != 9) { 
          ofile_FAIL("ERROR parsing $tblout_file for stage $stg_key, unexpected number of space-delimited tokens on line $line", 1, $FH_HR); 
        }
        ($seq, $model, $score, $s_from, $s_to, $strand) = ($el_A[0], $el_A[1], $el_A[2], $el_A[3], $el_A[4], $el_A[5]);
        # if we were parsing cmscan --trmF3 it would look like this:
        #modelname sequence                      score  start    end strand bounds ovp      seqlen
        ##--------- ---------------------------- ------ ------ ------ ------ ------ --- -----------
        #NC_039477  gi|1215708385|gb|KY594653.1|  275.8      1    301      +     []  *          301
        # and we would parse like this:
        #($model, $seq, $score, $s_from, $s_to, $strand) = ($el_A[0], $el_A[1], $el_A[2], $el_A[3], $el_A[4], $el_A[5]);
        if(! defined $results_HHHR->{$seq}) { 
          %{$results_HHHR->{$seq}} = ();
        }
        $group    = (defined $mdl_group_H{$model})    ? $mdl_group_H{$model}    : undef;
        $subgroup = (defined $mdl_subgroup_H{$model}) ? $mdl_subgroup_H{$model} : undef;

        # determine if we are going to store this hit, and to what 2D keys 
        # (the following code only works because we know we are sorted by score)
        my $is_1   = 0; # should we store this in $results_HHHR->{$seq}{"{rpn,std}.cls.1"} ? 
        my $is_2   = 0; # should we store this in $results_HHHR->{$seq}{"{rpn,std}.cls.2"} ? 
        my $is_eg  = 0; # should we store this in $results_HHHR->{$seq}{"{rpn,std}.cls.eg"} ? 
        my $is_esg = 0; # should we store this in $results_HHHR->{$seq}{"{rpn,std}.cls.esg"} ? 

        # store hit as cls.1 only if it's the first hit seen to any model, or it's an additional hit to the r1.1 model
        # SHOULD WE ONLY BE STORING TOP HIT IN CLASSIFICATION?
        if((! defined $results_HHHR->{$seq}{"$stg_key.1"}) || # first (best) hit for this sequence 
           (($results_HHHR->{$seq}{"$stg_key.1"}{"model"} eq $model) && 
            ($results_HHHR->{$seq}{"$stg_key.1"}{"bstrand"} eq $strand))) { # additional hit for r1.1 sequence/model/strand trio
          $is_1 = 1; 
        }
        # store hit as r1.2 only if it's an additional hit to the r1.2 model
        # or there is no r1.2 hit yet, and r1.1 model and r1.2 model are not in the same subgroup
        # (if either r1.1 or r1.2 models do not have a subgroup, they are considered to NOT be in the same subgroup)
        elsif((defined $results_HHHR->{$seq}{"$stg_key.2"}) && 
              (($results_HHHR->{$seq}{"$stg_key.2"}{"model"} eq $model) && 
               ($results_HHHR->{$seq}{"$stg_key.2"}{"bstrand"} eq $strand))) { # additional hit for r1.2 sequence/model/strand trio
          $is_2 = 1;
        }
        elsif((! defined $results_HHHR->{$seq}{"$stg_key.2"}) && # no r1.2 hit yet exists AND
              ((! defined $results_HHHR->{$seq}{"$stg_key.1"}{"subgroup"}) ||       # (r1.1 hit has no subgroup OR
               (! defined $subgroup) ||                                       #  this hit has no subgroup OR
               ($results_HHHR->{$seq}{"$stg_key.1"}{"subgroup"} ne $subgroup))) {   #  both have subgroups but they differ)
          $is_2 = 1;
        }               

        # determine if we are going to store this hit as best in 'group' and/or 'subgroup'
        # to the expected group and/or expected subgroup
        if((defined $exp_group) && (defined $group) && ($exp_group eq $group)) { 
          if((! defined $results_HHHR->{$seq}{"$stg_key.eg"}) || # first (top) hit for this sequence to this group
             (($results_HHHR->{$seq}{"$stg_key.eg"}{"model"} eq $model) && 
              ($results_HHHR->{$seq}{"$stg_key.eg"}{"bstrand"} eq $strand))) { # additional hit for r1.eg sequence/model/strand trio
            $is_eg = 1; 
          }
          if((defined $exp_subgroup) && (defined $subgroup) && ($exp_subgroup eq $subgroup)) { 
            if((! defined $results_HHHR->{$seq}{"$stg_key.esg"}) || # first (top) hit for this sequence to this subgroup
               (($results_HHHR->{$seq}{"$stg_key.esg"}{"model"} eq $model) && 
                ($results_HHHR->{$seq}{"$stg_key.esg"}{"bstrand"} eq $strand))) { # additional hit for r1.esg sequence/model/strand trio
              $is_esg = 1; 
            }
          }
        }
        if($is_1)   { cmsearch_store_hit(\%{$results_HHHR->{$seq}{"$stg_key.1"}},   $model, $score, $strand, $bias, $s_from, $s_to, $m_from, $m_to, $group, $subgroup, $FH_HR); }
        if($is_2)   { cmsearch_store_hit(\%{$results_HHHR->{$seq}{"$stg_key.2"}},   $model, $score, $strand, $bias, $s_from, $s_to, $m_from, $m_to, $group, $subgroup, $FH_HR); }
        if($is_eg)  { cmsearch_store_hit(\%{$results_HHHR->{$seq}{"$stg_key.eg"}},  $model, $score, $strand, $bias, $s_from, $s_to, $m_from, $m_to, $group, $subgroup, $FH_HR); }
        if($is_esg) { cmsearch_store_hit(\%{$results_HHHR->{$seq}{"$stg_key.esg"}}, $model, $score, $strand, $bias, $s_from, $s_to, $m_from, $m_to, $group, $subgroup, $FH_HR); }
      }
      else { # $stg_key eq "rpn.cdt" or "std.cdt": coverage determination, cmsearch --tblout output
        ##target name                 accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
        ##--------------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
        #gi|1215708385|gb|KY594653.1| -         NC_039477            -         hmm     5089     5389        1      301      +     -    6 0.52   0.0  268.1   3.4e-85 !   Norovirus GII.4 isolate Hu/GII/CR7410/CHN/2014 VP1 gene, partial cds
        if(scalar(@el_A) < 18) { 
          ofile_FAIL("ERROR parsing $tblout_file for stage $stg_key, unexpected number of space-delimited tokens on line $line", 1, $FH_HR); 
        }
        ($seq, $model, $m_from, $m_to, $s_from, $s_to, $strand, $bias, $score) = ($el_A[0], $el_A[2], $el_A[5], $el_A[6], $el_A[7], $el_A[8], $el_A[9], $el_A[13], $el_A[14]);
        # determine if we are going to store this hit
        if((! defined $results_HHHR->{$seq}{"$stg_key.bs"}) || # first (top) hit for this sequence, 
           (($results_HHHR->{$seq}{"$stg_key.bs"}{"model"}   eq $model) && 
            ($results_HHHR->{$seq}{"$stg_key.bs"}{"bstrand"} eq $strand))) { # additional hit for this sequence/model/strand trio
          cmsearch_store_hit(\%{$results_HHHR->{$seq}{"$stg_key.bs"}},  $model, $score, $strand, $bias, $s_from, $s_to, $m_from, $m_to, undef, undef, $FH_HR); # undefs are for group and subgroup which are irrelevant in coverage determination stage
        }
        elsif((! defined $results_HHHR->{$seq}{"$stg_key.os"}) || # first (top) hit on OTHER strand for this sequence, 
              (($results_HHHR->{$seq}{"$stg_key.os"}{"model"}   eq $model) && 
               ($results_HHHR->{$seq}{"$stg_key.os"}{"bstrand"} eq $strand))) { # additional hit for this sequence/model/strand trio
          cmsearch_store_hit(\%{$results_HHHR->{$seq}{"$stg_key.os"}},  $model, $score, $strand, $bias, $s_from, $s_to, $m_from, $m_to, undef, undef, $FH_HR); # undefs are for group and subgroup which are irrelevant in coverage determination stage
        }
      }
    }
  }

  return; 
}

#################################################################
# Subroutine:  cmsearch_store_hit()
# Incept:      EPN, Thu Mar 21 14:55:25 2019
#
# Purpose:     Store information on a hit in %{$HR}. 
#              Initialize %{$HR} if this is the first hit.
#
# Arguments: 
#  $HR:            hash to store hit info in
#  $model:         model name
#  $score:         score of hit
#  $strand:        strand of hit
#  $bias:          bias score
#  $s_from:        start position on sequence
#  $s_to:          end position on sequence
#  $m_from:        start position on model
#  $m_to:          end position on model
#  $group:         group for $model, can be undef
#  $subgroup:      subgroup for $model, can be undef
#  $FH_HR:         ref to file handle hash
# 
# Returns: void
# 
# Dies: if not the first hit, but $HR->{"model"} ne $model
#      
################################################################# 
sub cmsearch_store_hit { 
  my $sub_name = "cmsearch_store_hit()";
  my $nargs_expected = 12;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($HR, $model, $score, $strand, $bias, $s_from, $s_to, $m_from, $m_to, $group, $subgroup, $FH_HR) = @_;

  if(scalar(keys(%{$HR})) == 0) { # initialize
    $HR->{"model"}    = $model; # all hits stored in this hash will be to $model
    $HR->{"bstrand"}  = $strand; # strand of 'best hit' (first hit)
    $HR->{"s_coords"} = "";
    $HR->{"score"}    = "";
    if(defined $m_from) { $HR->{"m_coords"} = ""; }
    if(defined $bias)   { $HR->{"bias"}     = ""; }
  }
  else { 
    if($HR->{"model"} ne $model) { 
      ofile_FAIL("ERROR in $sub_name, trying to add additional hit to a different model", 1, $FH_HR);
    }
    $HR->{"s_coords"} .= ",";
    $HR->{"score"}    .= ",";
    if(defined $m_from) { $HR->{"m_coords"} .= ","; }
    if(defined $bias)   { $HR->{"bias"}     .= ","; }
  }
  $HR->{"s_coords"} .= $s_from . ".." . $s_to . ":" . $strand;
  $HR->{"score"}    .= sprintf("%.1f", $score);
  if(defined $m_from)   { $HR->{"m_coords"} .= $m_from . ".." . $m_to . ":+"; }
  if(defined $bias)     { $HR->{"bias"}     .= sprintf("%.1f", $bias); }
  if(defined $group)    { $HR->{"group"}     = $group; }
  if(defined $subgroup) { $HR->{"subgroup"}  = $subgroup; }

  return;
}

#################################################################
# Subroutine: add_classification_alerts()
# Incept:     EPN, Thu Mar 21 05:59:00 2019
# Purpose:    Adds c_* alerts for sequences based on classification
#             results in %{$stg_results_HHHR}.
#
# Types of 
# Arguments:
#  $alt_seq_instances_HHR:   REF to 2D hash with per-sequence alerts, added to here
#  $seq_len_HR:              REF to hash of sequence lengths
#  $mdl_info_AHR:            REF to array of hashes with information on the sequences, PRE-FILLED
#  $alt_info_HHR:            REF to the alert info hash of arrays, PRE-FILLED
#  $stg_results_HHHR:        REF to 3D hash with classification search results, PRE-FILLED
#  $cls_output_HHR:          REF to 2D hash of classification output info, FILLED HERE 
#  $opt_HHR:                 REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:          REF to the 2D hash of output file information
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub add_classification_alerts { 
  my $sub_name = "add_classification_alerts";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($alt_seq_instances_HHR, $seq_len_HR, $mdl_info_AHR, $alt_info_HHR, $stg_results_HHHR, $cls_output_HHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience
  my $nseq = scalar(keys (%{$seq_len_HR}));

  # create the model index hash which gives index in $mdl_info_AHR[] 
  # for a given model name, this allows us to find model length given model name
  my %mdl_idx_H = ();
  utl_IdxHFromAH(\%mdl_idx_H, $mdl_info_AHR, "name", $sub_name, $FH_HR);
   
  # determine if we have an expected group and subgroup
  # (--subgroup requires --group)
  # and if so, fill hashes of model groups and subgroups, 
  # (this data is in @{$mdl_info_AHR} already, we fill 
  # these hashes onlyfor convenience)
  my $exp_group    = opt_Get("--group",    \%opt_HH);
  my $exp_subgroup = opt_Get("--subgroup", \%opt_HH);

  # get thresholds
  my $small_value    = 0.00000001; # for handling precision issues
  my $lowcov_opt     = opt_Get("--lowcov",     $opt_HHR) - $small_value;
  my $lowsc_opt      = opt_Get("--lowsc",      $opt_HHR) - $small_value;
  my $indefstr_opt   = opt_Get("--indefstr",   $opt_HHR) + $small_value;
  my $indefclass_opt = opt_Get("--indefclass", $opt_HHR) - $small_value;
  my $biasfract_opt  = opt_Get("--biasfract",  $opt_HHR) + $small_value;
  my $incspec_opt    = opt_Get("--incspec",    $opt_HHR) + $small_value;
  my $dupregolp_opt  = opt_Get("--dupregolp",  $opt_HHR);
  my $dupregsc_opt   = opt_Get("--dupregsc",   $opt_HHR) - $small_value;

  my $lowcov_opt2print     = sprintf("%.3f", opt_Get("--lowcov",     $opt_HHR));
  my $lowsc_opt2print      = sprintf("%.3f", opt_Get("--lowsc",      $opt_HHR));
  my $indefstr_opt2print   = sprintf("%.1f", opt_Get("--indefstr",   $opt_HHR));
  my $indefclass_opt2print = sprintf("%.3f", opt_Get("--indefclass", $opt_HHR));
  my $biasfract_opt2print  = sprintf("%.3f", opt_Get("--biasfract",  $opt_HHR));
  my $incspec_opt2print    = sprintf("%.3f", opt_Get("--incspec",    $opt_HHR));
  my $dupregsc_opt2print   = sprintf("%.1f", opt_Get("--dupregsc",   $opt_HHR));

  # get info on position-specific indfstrn exceptions, if any
  my @dupregin_exc_AA = ();
  my @indfstrn_exc_AA = ();
  my $nmdl = scalar(@{$mdl_info_AHR});
  for($mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    @{$dupregin_exc_AA[$mdl_idx]} = ();
    @{$indfstrn_exc_AA[$mdl_idx]} = ();
    vdr_ModelInfoCoordsListValueBreakdown($mdl_info_AHR, $mdl_idx, "dupregin_exc", \@{$dupregin_exc_AA[$mdl_idx]}, $FH_HR);
    vdr_ModelInfoCoordsListValueBreakdown($mdl_info_AHR, $mdl_idx, "indfstrn_exc", \@{$indfstrn_exc_AA[$mdl_idx]}, $FH_HR);
  }

  my $alt_scoords; # sequence coordinates related to an alert
  my $alt_mcoords; # model    coordinates related to an alert

  # if we used blastn for the cdt stage, we may have overlapping hits in sequence coords, 
  # this is relevant if/when we call helper_sort_hit_array for the dupregin alert below
  my $do_blastn_cdt = opt_Get("-s", \%opt_HH) ? 1 : 0;

  %{$cls_output_HHR} = ();
  foreach my $seq_name (sort keys(%{$seq_len_HR})) { 
    my $seq_len  = $seq_len_HR->{$seq_name};
    my $mdl_name = undef;
    my $mdl_idx  = undef;
    my $mdl_len  = undef;
    my %score_H  = (); # key is $stg_results_HHHR 2D key (search category), value is summed score
    my %scpnt_H  = (); # key is $stg_results_HHHR 2D key (search category), value is summed length
    my $alt_str = "";
    %{$cls_output_HHR->{$seq_name}} = ();

    # check for noannotn alert: 3 possibilities
    # 1) no hits in round 1 search (most common cause of noannotn)
    # 2) >= 1 hits in -r       classification stage (rpn.cls.1)  but 0 hits in standard classification stage (std.cls.1) (rare)
    # 3) >= 1 hits in standard classification stage (std.cdt.bs) but 0 hits in coverage determination stage (std.cdt.bs) (rare)
    if((! defined $stg_results_HHHR->{$seq_name}) || # case 1
       ((defined $stg_results_HHHR->{$seq_name}) &&
        (defined $stg_results_HHHR->{$seq_name}{"rpn.cls.1"}) &&
        (! defined $stg_results_HHHR->{$seq_name}{"std.cls.1"})) || # case 2
       ((defined $stg_results_HHHR->{$seq_name}) &&
        (defined $stg_results_HHHR->{$seq_name}{"std.cls.1"}) &&
        (! defined $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}))) { # case 3
      alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "noannotn", $seq_name, "VADRNULL", $FH_HR);
    }
    else { 
      if(! defined $stg_results_HHHR->{$seq_name}{"std.cls.1"}) { 
        ofile_FAIL("ERROR in $sub_name, seq $seq_name should have but does not have any cls.1 hits", 1, $FH_HR);
      }
      # determine model name and length
      $mdl_name = $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"model"};
      $mdl_idx  = $mdl_idx_H{$mdl_name};
      $mdl_len  = $mdl_info_AHR->[$mdl_idx]{"length"};
      foreach my $rkey (keys (%{$stg_results_HHHR->{$seq_name}})) { 
        my @score_A = split(",", $stg_results_HHHR->{$seq_name}{$rkey}{"score"});
        $score_H{$rkey} = utl_ASum(\@score_A);
        $scpnt_H{$rkey} = $score_H{$rkey} / vdr_CoordsLength($stg_results_HHHR->{$seq_name}{$rkey}{"s_coords"}, $FH_HR);
      }
      my $have_cdt_bs = (defined $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}) ? 1 : 0;

      my $scpnt2print = sprintf("%.3f", $scpnt_H{"std.cls.1"});
      $cls_output_HHR->{$seq_name}{"scpnt"} = $scpnt2print;
      $cls_output_HHR->{$seq_name}{"score"} = sprintf("%.1f", $score_H{"std.cls.1"});

      # low score (lowscore)
      if($scpnt_H{"std.cls.1"} < $lowsc_opt) { 
        alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR,  "lowscore", $seq_name, $scpnt2print . "<" . $lowsc_opt2print . " bits/nt", $FH_HR);
      }

      # indefinite classification (indfclas))
      if(defined $scpnt_H{"std.cls.2"}) { 
        my $diffpnt = $scpnt_H{"std.cls.1"} - $scpnt_H{"std.cls.2"};
        my $diffpnt2print = sprintf("%.3f", $diffpnt);
        $cls_output_HHR->{$seq_name}{"scdiff"}  = sprintf("%.1f", $score_H{"std.cls.1"} - $score_H{"std.cls.2"});
        $cls_output_HHR->{$seq_name}{"diffpnt"} = $diffpnt2print;
        my $group_str = "best group/subgroup: " . 
            group_subgroup_string_from_classification_results($stg_results_HHHR->{$seq_name}{"std.cls.1"}) . 
            ", second group/subgroup: " . 
            group_subgroup_string_from_classification_results($stg_results_HHHR->{$seq_name}{"std.cls.2"});

        if($diffpnt < $indefclass_opt) { 
          alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "indfclas", $seq_name, $diffpnt2print . "<" . $indefclass_opt2print . " bits/nt, " . $group_str, $FH_HR);
        }
      }

      # incorrect group (incgroup) 
      # - $exp_group must be defined and group of cls.1 must be undef or != $exp_group
      # - no hits in cls.eg (no hits to group) (incgroup)
      # OR 
      # - hit(s) in cls.eg but scpernt diff between
      #   cls.eg and cls.1 exceeds incspec_opt (incgroup)
      #
      # questionable group (qstgroup)
      # - $exp_group must be defined 
      # - hit(s) in cls.eg but scpernt diff between
      #   cls.eg and cls.1 does not exceed incspec_opt (qstgroup)
      #
      my $igr_flag = 0;
      my $qgr_flag = 0;
      if((defined $exp_group) && # $exp_group defined
         ((! defined $stg_results_HHHR->{$seq_name}{"std.cls.1"}{"group"}) || # cls.1 group undefined
          ($stg_results_HHHR->{$seq_name}{"std.cls.1"}{"group"} ne $exp_group))) { # cls.1 group != $exp_group
        # $exp_group is defined AND 
        # (cls.1 group undefined OR cls.1 group != $exp_group)
        if(! defined $scpnt_H{"std.cls.eg"}) { 
          # no hit to $exp_group
          alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "incgroup", $seq_name, 
                                      "no hits to expected group $exp_group, best model group/subgroup: " . 
                                      group_subgroup_string_from_classification_results($stg_results_HHHR->{$seq_name}{"std.cls.1"}), $FH_HR);
          $igr_flag = 1;
        }
        else { 
          # at least one hit to $exp_group exists and cls.1's group undef or != $exp_group 
          # so we either have a incgroup or qstgroup alert
          my $diff = $scpnt_H{"std.cls.1"} - $scpnt_H{"std.cls.eg"};
          my $diff2print = sprintf("%.3f", $diff);
          if($diff > $incspec_opt) { 
            $alt_str = "$diff2print > $incspec_opt2print bits/nt diff, best model group/subgroup: " . group_subgroup_string_from_classification_results($stg_results_HHHR->{$seq_name}{"std.cls.1"});
            alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "incgroup", $seq_name, $alt_str, $FH_HR);
            $igr_flag = 1;
          }
          else { 
            $alt_str = "$diff2print bits/nt diff, best model group/subgroup: " . group_subgroup_string_from_classification_results($stg_results_HHHR->{$seq_name}{"std.cls.1"});
            alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "qstgroup", $seq_name, $alt_str, $FH_HR);
            $qgr_flag = 1;
          }
        }
      }

      # incorrect subgroup (c_sgr) 
      # - $exp_subgroup must be defined and subgroup of cls.1 must be undef or != $exp_subgroup
      # - incgroup not already reported
      # - no hits in cls.esg (no hits to group) (incsbgrp)
      # OR 
      # - hit(s) in cls.esg but scpernt diff between
      #   cls.esg and cls.1 exceeds incspec_opt (incsbgrp)
      #
      # questionable subgroup (qstsbgrp)
      # - $exp_subgroup must be defined 
      # - incgroup not already reported
      # - qstgroup not already reported
      # - hit(s) in cls.esg but scpernt diff between
      #   cls.esg and cls.1 does not exceed incspec_opt (qstsbgrp)
      #
      if((! $igr_flag) # incgroup alert not reported
         && (defined $exp_subgroup) && # exp_sugroup defined  
         ((! defined $stg_results_HHHR->{$seq_name}{"std.cls.1"}{"subgroup"}) || # cls.1 subgroup undefined
          ($stg_results_HHHR->{$seq_name}{"std.cls.1"}{"subgroup"} ne $exp_subgroup))) { # cls.1 subgroup != $exp_subgroup
        # incgroup alert not reported AND
        # $exp_subgroup is defined AND 
        # (cls.1 subgroup undefined OR cls.1 subgroup != $exp_subgroup)
        if(! defined $scpnt_H{"std.cls.esg"}) { 
          alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "incsbgrp", $seq_name, "no hits to expected subgroup $exp_subgroup", $FH_HR);
        }
        else { 
          my $diff = $scpnt_H{"std.cls.1"} - $scpnt_H{"std.cls.esg"};
          my $diff2print = sprintf("%.3f", $diff);
          if($diff > $incspec_opt) { 
            $alt_str = $diff2print . ">" . $incspec_opt2print . " bits/nt diff, best model group/subgroup: " . group_subgroup_string_from_classification_results($stg_results_HHHR->{$seq_name}{"std.cls.1"});
            alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "incsbgrp", $seq_name, $alt_str, $FH_HR);
          }
          elsif(! $qgr_flag) {
            $alt_str = "$diff2print bits/nt diff, best model group/subgroup: " . group_subgroup_string_from_classification_results($stg_results_HHHR->{$seq_name}{"std.cls.1"});
            alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "qstsbgrp", $seq_name, $alt_str, $FH_HR);
          }
        }
      }

      # classification alerts that depend on round 2 results
      if($have_cdt_bs) { 
        my @bias_A   = split(",", $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"bias"});
        my $bias_sum = utl_ASum(\@bias_A);
        my $bias_fract = undef;
        if(($bias_sum <= $small_value) || ($score_H{"std.cdt.bs"} < $small_value)) { 
         # (a) bias is 0 or negative OR (b) score_H{"std.cdt.bs"} is 0 
         # (note: the sum ($score_H{"std.cdt.bs"} + $bias_sum) can't be 0 or negative if a or b isn't true) 
          $bias_fract = 0.; # bias fraction doesn't really make sense if score is negative
                            # this also ensures we don't try to divide by 0 (i.e. if ($score_H{"std.cdt.bs"} + $bias_sum) == 0)
        }
        else {
          $bias_fract = $bias_sum / ($score_H{"std.cdt.bs"} + $bias_sum);
        }
        my $nhits = scalar(@bias_A);
        $cls_output_HHR->{$seq_name}{"nhits"}   = $nhits;
        $cls_output_HHR->{$seq_name}{"bias"}    = $bias_sum;
        $cls_output_HHR->{$seq_name}{"bstrand"} = $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"bstrand"};
        my $s_len = vdr_CoordsLength($stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"}, $FH_HR);
        my $m_len = vdr_CoordsLength($stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"m_coords"}, $FH_HR);
        my $scov = $s_len / $seq_len;
        my $scov2print = sprintf("%.3f", $scov);
        my $mcov2print = sprintf("%.3f", $m_len / $mdl_len);
        $cls_output_HHR->{$seq_name}{"scov"} = $scov2print;
        $cls_output_HHR->{$seq_name}{"mcov"} = $mcov2print;
      
        # reverse complement (revcompl)
        if($stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"bstrand"} eq "-") { 
          $alt_scoords = "seq:" . $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"} . ";";
          $alt_mcoords = "mdl:" . $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"m_coords"} . ";";
          my @bstrand_score_A = split(",", $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"score"});
          $alt_str = sprintf("best_hit_score:%.1f", $bstrand_score_A[0]);
          alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "revcompl", $seq_name, $alt_scoords . $alt_mcoords . $alt_str, $FH_HR);
        }

        # low coverage (lowcovrg)
        if($scov < $lowcov_opt) { 
          $alt_scoords = "seq:" . vdr_CoordsMissing($stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"}, $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"bstrand"}, $seq_len, $FH_HR) . ";";
          $alt_mcoords = "mdl:VADRNULL;";
          $alt_str = sprintf("%s<%s", $scov2print, $lowcov_opt2print);
          alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "lowcovrg", $seq_name, $alt_scoords . $alt_mcoords . $alt_str, $FH_HR);
        }

        # high bias (biasdseq) 
        if($bias_fract > $biasfract_opt) { 
          my $bias_fract2print = sprintf("%.3f", $bias_fract);
          alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "biasdseq", $seq_name, $bias_fract2print . ">" . $biasfract_opt2print, $FH_HR);
        }

        # inconsistent hits: multiple strands (indfstrn) 
        if(defined $stg_results_HHHR->{$seq_name}{"std.cdt.os"}) { 
          my @ostrand_score_A = split(",", $stg_results_HHHR->{$seq_name}{"std.cdt.os"}{"score"});
          my $top_ostrand_score = $ostrand_score_A[0];
          if($top_ostrand_score > $indefstr_opt) { 
            my @ostrand_sstart_A  = ();
            my @ostrand_sstop_A   = ();
            my @ostrand_sstrand_A = ();
            my @ostrand_mstart_A  = ();
            my @ostrand_mstop_A   = ();
            my @ostrand_mstrand_A = ();
            vdr_FeatureStartStopStrandArrays($stg_results_HHHR->{$seq_name}{"std.cdt.os"}{"s_coords"}, \@ostrand_sstart_A, \@ostrand_sstop_A, \@ostrand_sstrand_A, $FH_HR);
            vdr_FeatureStartStopStrandArrays($stg_results_HHHR->{$seq_name}{"std.cdt.os"}{"m_coords"}, \@ostrand_mstart_A, \@ostrand_mstop_A, \@ostrand_mstrand_A, $FH_HR);
            # check if this is an exempted region
            my $exempted_region = 0;
            foreach my $exc_coords (@{$indfstrn_exc_AA[$mdl_idx]}) { 
              if(vdr_CoordsCheckIfSpans($exc_coords, vdr_CoordsSegmentCreate($ostrand_mstart_A[0], $ostrand_mstop_A[0], $ostrand_mstrand_A[0], $FH_HR), $FH_HR)) { 
                $exempted_region = 1;
              }
            }
            if(! $exempted_region) { 
              $alt_scoords = "seq:" . vdr_CoordsSegmentCreate($ostrand_sstart_A[0], $ostrand_sstop_A[0], $ostrand_sstrand_A[0], $FH_HR) . ";";
              $alt_mcoords = "mdl:" . vdr_CoordsSegmentCreate($ostrand_mstart_A[0], $ostrand_mstop_A[0], $ostrand_mstrand_A[0], $FH_HR) . ";";
              $alt_str = sprintf("score:%.1f>%s", $top_ostrand_score, $indefstr_opt2print);
              alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "indfstrn", $seq_name, $alt_scoords . $alt_mcoords . $alt_str, $FH_HR);
            }
          }
        }

        # inconsistent hits: duplicate regions (dupregin) 
        $alt_str = "";
        if($nhits > 1) { 
          my @m_start_A  = ();
          my @m_stop_A   = ();
          my @m_strand_A = ();
          my @s_start_A  = ();
          my @s_stop_A   = ();
          my @s_strand_A = ();
          $alt_scoords = "";
          $alt_mcoords = "";
          vdr_FeatureStartStopStrandArrays($stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"m_coords"}, \@m_start_A, \@m_stop_A, \@m_strand_A, $FH_HR);
          vdr_FeatureStartStopStrandArrays($stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"}, \@s_start_A, \@s_stop_A, \@s_strand_A, $FH_HR);
          my @dupreg_score_A = split(",", $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"score"});
          for(my $i = 0; $i < $nhits; $i++) { 
            if($dupreg_score_A[$i] > $dupregsc_opt) { 
              for(my $j = $i+1; $j < $nhits; $j++) { 
                if($dupreg_score_A[$j] > $dupregsc_opt) { 
                  my ($mdl_noverlap, $mdl_overlap_str) = vdr_CoordsSegmentOverlap(
                    vdr_CoordsSegmentCreate($m_start_A[$i], $m_stop_A[$i], $m_strand_A[$i], $FH_HR), 
                    vdr_CoordsSegmentCreate($m_start_A[$j], $m_stop_A[$j], $m_strand_A[$j], $FH_HR), $FH_HR);
                  my ($seq_noverlap, $seq_overlap_str) = vdr_CoordsSegmentOverlap(
                    vdr_CoordsSegmentCreate($s_start_A[$i], $s_stop_A[$i], $s_strand_A[$i], $FH_HR), 
                    vdr_CoordsSegmentCreate($s_start_A[$j], $s_stop_A[$j], $s_strand_A[$j], $FH_HR), $FH_HR);
                  my $noverlap = $mdl_noverlap - $seq_noverlap;
                  if($noverlap >= $dupregolp_opt) { 
                    $alt_scoords .= sprintf("%s%s,%s", 
                                            ($alt_scoords eq "") ? "seq:" : ",",
                                            vdr_CoordsSegmentCreate($s_start_A[$i], $s_stop_A[$i], $s_strand_A[$i], $FH_HR), 
                                            vdr_CoordsSegmentCreate($s_start_A[$j], $s_stop_A[$j], $s_strand_A[$j], $FH_HR));
                    $alt_mcoords .= sprintf("%s%s,%s", 
                                            ($alt_mcoords eq "") ? "mdl:" : ",",
                                            vdr_CoordsSegmentCreate($m_start_A[$i], $m_stop_A[$i], $m_strand_A[$i], $FH_HR), 
                                            vdr_CoordsSegmentCreate($m_start_A[$j], $m_stop_A[$j], $m_strand_A[$j], $FH_HR));
                    $mdl_overlap_str =~ s/\-/\.\./; # replace '-' with '..', e.g. '10-15' to '10..15'
                    $mdl_overlap_str .= ":" . $m_strand_A[$i]; # both $i and $j are same strand
                    # check if this is an exempted region
                    my $exempted_region = 0;
                    foreach my $exc_coords (@{$dupregin_exc_AA[$mdl_idx]}) { 
                      if(vdr_CoordsCheckIfSpans($exc_coords, $mdl_overlap_str, $FH_HR)) { 
                        $exempted_region = 1;
                      }
                    }
                    if(! $exempted_region) { # only report if not exempted
                      $alt_str .= sprintf("%s%s (len %d>=%d) hits %d (%.1f bits) and %d (%.1f bits)",
                                          ($alt_str eq "") ? "" : ", ",
                                          $mdl_overlap_str, $noverlap, $dupregolp_opt, 
                                          ($i+1), $dupreg_score_A[$i], 
                                          ($j+1), $dupreg_score_A[$j]);
                    }
                  }
                }
              }
            }
          }
          if($alt_str ne "") { 
            alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "dupregin", $seq_name, $alt_scoords . ";" . $alt_mcoords . ";" . $alt_str, $FH_HR);
          }
        }
      
        # inconsistent hits: wrong hit order (discontn)
        if($nhits > 1) { 
          my $i;
          my @seq_hit_order_A = (); # array of sequence boundary hit indices in sorted order [0..nhits-1] values are in range 1..nhits
          my @mdl_hit_order_A = (); # array of model    boundary hit indices in sorted order [0..nhits-1] values are in range 1..nhits
          my @seq_hit_coords_A = split(",", $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"});
          my @mdl_hit_coords_A = split(",", $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"m_coords"});
          my $seq_hit_order_str = undef;
          # if blastn was used, we allow overlaps in the seq hits because blastn can report these but cmsearch cannot
          if($do_blastn_cdt) { 
            $seq_hit_order_str = helper_sort_hit_array(\@seq_hit_coords_A, \@seq_hit_order_A, 1, $FH_HR); # 1 means duplicate values in best array are not allowed
          }
          else { 
            $seq_hit_order_str = helper_sort_hit_array(\@seq_hit_coords_A, \@seq_hit_order_A, 0, $FH_HR); # 0 means duplicate values in best array are not allowed
          }
          my $mdl_hit_order_str = helper_sort_hit_array(\@mdl_hit_coords_A, \@mdl_hit_order_A, 1, $FH_HR); # 1 means duplicate values in best array are allowed
          # check if the hits are out of order we don't just check for equality of the
          # two strings because it's possible (but rare) that there could be duplicates in the model
          # order array and sequence order array, so we need to allow for that.
          my $out_of_order_flag = 0;
          for($i = 0; $i < $nhits; $i++) { 
            my $x = $mdl_hit_order_A[$i];
            my $y = $seq_hit_order_A[$i];
            # check to see if hit $i is same order in both mdl and seq coords
            # or if it is not, it's okay if it is identical to the one that is
            # example: 
            # hit 1 seq 1..10,+   model  90..99,+
            # hit 2 seq 11..20,+  model 100..110,+
            # hit 3 seq 21..30,+  model 100..110,+
            # seq order: 1,2,3
            # mdl order: 1,3,2 (or 1,2,3) we want both to be ok (not FAIL)
            if($x ne $y) { # hits are not the same order
              my $mdl_identical_flag = ($mdl_hit_coords_A[($x-1)] eq $mdl_hit_coords_A[($y-1)]) ? 1 : 0;
              my $seq_identical_flag = ($seq_hit_coords_A[($x-1)] eq $seq_hit_coords_A[($y-1)]) ? 1 : 0;
              if($mdl_identical_flag && $seq_identical_flag) { 
                ofile_FAIL("ERROR in $sub_name, found two hits identical in both seq and mdl coords for seq $seq_name seq_coords: " . $seq_hit_coords_A[($x-1)] . ", mdl_coords: " . $mdl_hit_coords_A[($x-1)], 1, $FH_HR);
              }
              if((! $mdl_identical_flag) && (! $seq_identical_flag)) { 
                # hit is not identical in either mdl or seq coords to hit in correct order
                $out_of_order_flag = 1;
                $i = $nhits; # breaks 'for i' loop, slight optimization
              }
            }
          }
          if($out_of_order_flag) { 
#            $alt_str = "seq order: " . $seq_hit_order_str . "(" . $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"} . ")";
#            $alt_str .= ", model order: " . $mdl_hit_order_str . "(" . $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"m_coords"} . ")";
            $alt_str = sprintf("%s%sseq order: %s, mdl order: %s",
                               "seq:" . $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"} . ";", 
                               "mdl:" . $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"m_coords"} . ";", 
                               $seq_hit_order_str, 
                               $mdl_hit_order_str);
            alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "discontn", $seq_name, $alt_str, $FH_HR);
          }
        }
      }
      
      # finally fill $cls_output_HHR->{$seq_name} info related to models and groups
      if(defined $stg_results_HHHR->{$seq_name}{"std.cls.1"}) { 
        $cls_output_HHR->{$seq_name}{"model1"}    = $stg_results_HHHR->{$seq_name}{"std.cls.1"}{"model"};
        $cls_output_HHR->{$seq_name}{"group1"}    = $stg_results_HHHR->{$seq_name}{"std.cls.1"}{"group"};    # can be undef
        $cls_output_HHR->{$seq_name}{"subgroup1"} = $stg_results_HHHR->{$seq_name}{"std.cls.1"}{"subgroup"}; # can be undef
      }
      if(defined $stg_results_HHHR->{$seq_name}{"std.cls.2"}) { 
        $cls_output_HHR->{$seq_name}{"model2"}    = $stg_results_HHHR->{$seq_name}{"std.cls.2"}{"model"};
        $cls_output_HHR->{$seq_name}{"group2"}    = $stg_results_HHHR->{$seq_name}{"std.cls.2"}{"group"};    # can be undef
        $cls_output_HHR->{$seq_name}{"subgroup2"} = $stg_results_HHHR->{$seq_name}{"std.cls.2"}{"subgroup"}; # can be undef
      }
      # save -r data (which may differ from std.cls data, especially if --r_list or --r_only)
      if(defined $stg_results_HHHR->{$seq_name}{"rpn.cls.1"}) { 
        $cls_output_HHR->{$seq_name}{"rpn.model1"} = $stg_results_HHHR->{$seq_name}{"rpn.cls.1"}{"model"};
      }
    } # else entered if we didn't report a noannotn alert
  } # end of foreach seq loop

  return;
}

#################################################################
# Subroutine:  populate_per_model_data_structures_given_classification_results
# Incept:      EPN, Tue Apr  2 11:33:29 2019
#
# Purpose:    Given classification results, fill several 'per-model'
#             data structures with sequences that match best to that
#             model. Also fill the 'per-sequence' hash %{$seq2mdl_HR}
#             that maps each sequence to the model it is classified to.
# 
#             This subroutine is called twice, once following round
#             1 classification and once following round 2 classification
#             and does extra work following round 2. We determine
#             which round of classification was just performed
#             based on if $alt_seq_instances_HHR is defined or not.
#             If it is defined we just performed round 2, else round 1.
#
# Arguments: 
#  $seq_name_AR:           ref to sequence names
#  $seq_len_HR:            ref to hash of sequence lengths
#  $stg_key:               stage key, "rpn.cdt" or "std.cdt" if $alt_seq_instances_HHR is undef
#                          or "std.aln" if $alt_seq_instances_HHR is defined
#  $stg_results_HHHR:      ref to 3D hash of classification results, PRE-FILLED
#  $cls_output_HHR:        ref to 2D hash of classification results to output, PRE-FILLED
#                          can be undef if $alt_seq_instances_HHR is undef
#  $alt_info_HHR:          ref to 2d hash of alert info, PRE-FILLED
#                          can be undef if $alt_seq_instances_HHR is undef
#  $alt_seq_instances_HHR: ref to hash of per-seq alert instances, PRE-FILLED
#                          undef to use stg_results_HHHR->{}{"cls.1"} # round 1 search
#  $mdl_sub_HR:            ref to hash of of model substitutions, PRE-FILLED, should be undef unless --msub used
#  $mdl_seq_name_HAR:      ref to hash of arrays of sequences per model, FILLED HERE
#  $mdl_seq_len_HR:        ref to hash of summed sequence length per model, FILLED HERE
#  $mdl_ct_HR:             ref to hash of number of sequences per model, FILLED HERE, can be undef if $stg_key eq "rpn.cdt"
#  $seq2mdl_HR:            ref to hash of classified model per sequence, FILLED HERE
#  $FH_HR:                 ref to hash of file handles
#
# Returns:   void
#
#
################################################################# 
sub populate_per_model_data_structures_given_classification_results {
  my $sub_name = "populate_per_model_data_structures_given_classification_results";
  my $nargs_exp = 13;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name_AR, $seq_len_HR, $stg_key, $stg_results_HHHR, 
      $cls_output_HHR, $alt_info_HHR, $alt_seq_instances_HHR, $mdl_sub_HR, 
      $mdl_seq_name_HAR, $mdl_seq_len_HR, $mdl_ct_HR, $seq2mdl_HR, $FH_HR) = @_;

  # determine what stage results we are using and check that combo of 
  # stage key and defined/undefined of %{$alt_seq_instances_HHR} is valid
  my $cls_2d_key = undef;
  if(($stg_key eq "rpn.cdt") || ($stg_key eq "std.cdt")) { 
    if(defined $alt_seq_instances_HHR) { 
      ofile_FAIL("ERROR in $sub_name, stage key is rpn.cdt or std.cdt but alt_seq_instances_HHR is not defined", 1, $FH_HR);
    }
    $cls_2d_key = ($stg_key eq "rpn.cdt") ? "rpn.cls.1" : "std.cls.1";
  }
  elsif($stg_key eq "std.aln") { 
    if(! defined $alt_seq_instances_HHR) { 
      ofile_FAIL("ERROR in $sub_name, stage key is std.aln but alt_seq_instances_HHR is not defined", 1, $FH_HR);
    }
    $cls_2d_key = "std.cdt.bs";
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unrecognized stage key: $stg_key, should be rpn.cdt, std.cdt or std.aln", 1, $FH_HR);
  }
  if(($stg_key ne "rpn.cdt") && (! defined $mdl_ct_HR)) { 
    ofile_FAIL("ERROR in $sub_name, stage key is not rpn.cdt but mdl_ct_HR is not defined", 1, $FH_HR);
  }

  my $nseq = scalar(@{$seq_name_AR});
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name = $seq_name_AR->[$seq_idx];
    my $mdl_name = ((defined $stg_results_HHHR->{$seq_name}) && 
                    (defined $stg_results_HHHR->{$seq_name}{$cls_2d_key}) && 
                    (defined $stg_results_HHHR->{$seq_name}{$cls_2d_key}{"model"})) ? 
                    $stg_results_HHHR->{$seq_name}{$cls_2d_key}{"model"} : undef;
    if(defined $mdl_name) { 
      # check if we need to substitute another model for this one (--msub option)
      if((defined $mdl_sub_HR) && (defined $mdl_sub_HR->{$mdl_name})) { 
        $mdl_name = $mdl_sub_HR->{$mdl_name};
      }
      # determine if we are going to add this sequence to our per-model hashes, depending on what round
      my $add_seq = 0;
      if(($cls_2d_key eq "std.cdt.1") || ($cls_2d_key eq "rpn.cdt.1")) { 
        $add_seq = 1; # always add seq after classification round
      }
      else { 
        # if $cls_2d_key is "std.aln", check if this sequence has any alerts that prevent annotation
        $add_seq = (alert_instances_check_prevents_annot($seq_name, $alt_info_HHR, $alt_seq_instances_HHR, $FH_HR)) ? 0 : 1;
        # update "annot" key
        $cls_output_HHR->{$seq_name}{"annot"} = ($add_seq) ? 1 : 0;
      }
      if($add_seq) { 
        if(! defined $mdl_seq_name_HAR->{$mdl_name}) { 
          @{$mdl_seq_name_HAR->{$mdl_name}} = ();
          $mdl_seq_len_HR->{$mdl_name} = 0;
          if(defined $mdl_ct_HR) { $mdl_ct_HR->{$mdl_name} = 0; }
        }
        push(@{$mdl_seq_name_HAR->{$mdl_name}}, $seq_name);
        $mdl_seq_len_HR->{$mdl_name} += $seq_len_HR->{$seq_name};
        if(defined $mdl_ct_HR) { $mdl_ct_HR->{$mdl_name}++; }
        $seq2mdl_HR->{$seq_name} = $mdl_name;
      }
    }
  }

  return;
}

#################################################################
#
# Subroutines related to cmalign and alignment:
# cmalign_or_glsearch_wrapper
# cmalign_or_glsearch_wrapper_helper
# cmalign_or_glsearch_run
# parse_stk_and_add_alignment_cds_and_mp_alerts 
# cmalign_store_overflow
# fetch_features_and_add_cds_and_mp_alerts_for_one_sequence
# sqstring_check_start
# sqstring_find_stops 
#
#################################################################
# Subroutine:  cmalign_or_glsearch_wrapper()
# Incept:      EPN, Mon Mar 18 14:20:56 2019
#
# Purpose:     Run one or more cmalign or glsearch jobs on the farm
#              or locally, after possibly splitting up the input
#              sequence file with vdr_SplitFastaFile and 
#              then calling cmalign_or_glsearch_wrapper_helper().
#              We run glsearch if $mdl_file ends in '.hmm' else we run
#              cmalign.
# 
#              We may have to do two rounds of sequence file splitting
#              and job running/submission because there is an error
#              case in cmalign that we want to be able to detect. That
#              error case is when the sequence requires too much
#              memory to align. In order to catch those error cases we
#              need to run each offending sequence individually, so
#              our strategy for cmalign is:
#
#              Split full fasta file up using default method and run
#              >= 1 cmalign jobs. If any of those R runs fail, then 
#              split up run R's sequence file into >= 1 files with
#              exactly 1 sequence in them. One or more of those should
#              fail and that reveals which specific sequences are
#              causing the memory overflow.
# 
# Arguments: 
#  $execs_HR:              ref to executables with "esl-ssplit" and "cmalign"
#                          defined as keys
#  $qsub_prefix:           qsub command prefix to use when submitting to farm, undef if running locally
#  $qsub_suffix:           qsub command suffix to use when submitting to farm, undef if running locally
#  $blastn_db_sqfile_R:    ref to Bio::Easel::SqFile object with glsearch target seqs (model seqs)
#  $mdl_file:              name of model file to use (if ends with .fa, run glsearch, else run cmalign)
#  $mdl_name:              name of model to fetch from $mdl_file (undef to not fetch)
#  $seq_file:              name of sequence file with all sequences to run against
#  $out_root:              string for naming output files
#  $extra_key:             extra key for output file names, "" or "uj." (latter for seqs with unjoinbl alerts)
#  $nseq:                  total number of all seqs in $seq_file
#  $tot_len_nt:            total length of all nucleotides in $seq_file
#  $progress_w:            width for outputProgressPrior output
#  $stk_file_AR:           ref to array of stockholm files created here, FILLED HERE
#  $overflow_seq_AR:       ref to array of sequences that failed due to matrix overflows, FILLED HERE
#  $overflow_mxsize_AR:    ref to array of required matrix sizes for each sequence that failed due to matrix overflows, FILLED HERE
#  $opt_HHR:               ref to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:        ref to 2D hash of output file information
#
# Returns:     void
# 
# Dies: If an executable doesn't exist, or cmalign or nhmmscan or esl-ssplit
#       command fails if we're running locally
################################################################# 
sub cmalign_or_glsearch_wrapper { 
  my $sub_name = "cmalign_or_glsearch_wrapper";
  my $nargs_expected = 17;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $qsub_prefix, $qsub_suffix, $blastn_db_sqfile_R,
      $mdl_file, $mdl_name, $seq_file, $out_root, $extra_key, 
      $nseq, $tot_len_nt, $progress_w, $stk_file_AR, $overflow_seq_AR, 
      $overflow_mxsize_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $do_glsearch = ($mdl_file =~ m/\.fa$/) ? 1 : 0;
  my $nfasta_created = 0; # number of fasta files created by esl-ssplit
  my $log_FH      = $ofile_info_HHR->{"FH"}{"log"}; # for convenience
  my $start_secs; # timing start
  my $do_parallel = opt_Get("-p", $opt_HHR);
  my $do_keep     = opt_Get("--keep", $opt_HHR);
  my $do_cmindi   = opt_Get("--cmindi", $opt_HHR);
  @{$overflow_seq_AR} = (); # we will fill this with names of sequences that fail cmalign because
                            # the matrix required to align them is too big

  # set up output file names for concatenation
  my @concat_keys_A = (); # %r{1,2}_out_file_HAR keys we are going to concatenate files for
  my %concat_HA = ();     # hash of arrays of all files to concatenate together
  my $out_key;            # key for an output file: e.g. "stdout", "ifile", "tfile", "tblout", "err", "sh"
  
  push(@concat_keys_A, "stdout"); 
  push(@concat_keys_A, "ifile"); 
  if($do_parallel || $do_glsearch) { 
    push(@concat_keys_A, "err"); 
  }
  if($do_parallel) { 
    push(@concat_keys_A, "sh"); 
  }
  foreach $out_key (@concat_keys_A) { 
    @{$concat_HA{$out_key}} = ();
  }    

  my $nr1 = 0; # number of runs in round 1 (one per sequence file we create)
  my @r1_out_file_AH = (); # array of hashes ([0..$nr1-1]) of output files for cmalign round 1 runs
  my @r1_success_A   = (); # [0..$nr1-1]: '1' if this run finishes successfully, '0' if not
  my @r1_mxsize_A    = (); # [0..$nr1-1]: if $r1_success_A[$r1_i] is '0', required size for dp mx, else '0'
  my @r1_seq_file_A  = (); # [0..$nr1-1]: name of sequence file for this run
  my $r1_i;                # counter over round 1 runs
  my $r1_do_split = 0;     # set to '1' if we split up fasta file
  # we need to split up the sequence file, and submit a separate set of cmalign jobs for each
  my $targ_nseqfiles = vdr_SplitNumSeqFiles($tot_len_nt, $opt_HHR);
  if(($targ_nseqfiles > 1) || ($do_cmindi)) { # we are going to split up the fasta file 
    if($do_cmindi) { $targ_nseqfiles = -1; } # makes vdr_SplitFastaFile create 1 file per seq
    $r1_do_split = 1;
    $nr1 = vdr_SplitFastaFile($execs_HR->{"esl-ssplit"}, $seq_file, $targ_nseqfiles, undef, $opt_HHR, $ofile_info_HHR);
    # vdr_SplitFastaFile will return the actual number of fasta files created, 
    # which can differ from the requested amount (which is $targ_nseqfiles) that we pass in. 
    for($r1_i = 0; $r1_i < $nr1; $r1_i++) { # update sequence file names
      $r1_seq_file_A[$r1_i] = $seq_file . "." . ($r1_i+1);
    }
  }
  else { # not going to split up the sequence file
    $nr1 = 1;
    $r1_seq_file_A[0] = $seq_file;
  }
  
  cmalign_or_glsearch_wrapper_helper($execs_HR, $mdl_file, $mdl_name, $out_root, 1, $extra_key, $nseq, $progress_w, 
                                     \@r1_seq_file_A, \@r1_out_file_AH, \@r1_success_A, \@r1_mxsize_A, 
                                     $opt_HHR, $ofile_info_HHR);

  my $nr2            = 0;  # number of round 2 runs (sequence files)
  my @r2_out_file_AH = (); # array of hashes ([0..$nr2-1]) of output files for cmalign round 2 runs
  my @r2_success_A   = (); # [0..$nr2-1]: '1' if this run finishes successfully, '0' if not
  my @r2_mxsize_A    = (); # [0..$nr2-1]: if $r2_success_A[$r2_i] is '0', required size for dp mx, else '0'
  my @r2_seq_file_A  = (); # [0..$nr2-1]: name of sequence file for this run
  my $r2_i;                # counter over round 2 runs

  # go through each run:
  # if it finished successfully record its output files to concatenate later
  # if it did not finish successfully rerun all of its sequences (if $do_cmalign)
  for($r1_i = 0; $r1_i < $nr1; $r1_i++) { 
    if($r1_success_A[$r1_i]) { 
      # run finished successfully
      # if $do_glsearch, create the stockholm output file and insert file
      if($do_glsearch) { 
        my $glsearch_nstk = vdr_GlsearchFormat3And9CToStockholmAndInsertFile($r1_out_file_AH[$r1_i]{"stdout"}, 
                                                                             $r1_out_file_AH[$r1_i]{"stk"},
                                                                             $r1_out_file_AH[$r1_i]{"ifile"},
                                                                             $blastn_db_sqfile_R, 
                                                                             $mdl_name, $opt_HHR, $ofile_info_HHR);
        # add each individual sequence alignment file to our list of stk files
        # we use one alignment file per sequence with --glsearch because only 1
        # seq can exist in an alignment to deal with 'insert' doctor cases
        for(my $z = 1; $z <= $glsearch_nstk; $z++) { 
          push(@{$stk_file_AR}, $r1_out_file_AH[$r1_i]{"stk"} . "." . $z);
        }
      }
      else { # $do_glsearch is 0, using cmalign, so we have 1 alignment file for all seqs
        push(@{$stk_file_AR}, $r1_out_file_AH[$r1_i]{"stk"});
      }
      foreach $out_key (@concat_keys_A) { 
        push(@{$concat_HA{$out_key}}, $r1_out_file_AH[$r1_i]{$out_key});
      }
    }
    else { 
      # run did not finish successfully
      # split this sequence file up into multiple files with only 1 sequence each, 
      my $cur_nr2 = vdr_SplitFastaFile($execs_HR->{"esl-ssplit"}, $r1_seq_file_A[$r1_i], -1, undef, $opt_HHR, $ofile_info_HHR);
      if($cur_nr2 == 1) { 
        # special case, r1 sequence file had only 1 sequence, so we know the culprit
        # and don't need to rerun cmalign
        cmalign_store_overflow($r1_seq_file_A[$r1_i], $r1_mxsize_A[$r1_i], $overflow_seq_AR, $overflow_mxsize_AR, $ofile_info_HHR->{"FH"}); 
        unlink $r1_seq_file_A[$r1_i] . ".1"; # remove the file we just created 
        $nr2 = 0;
      }
      else { 
        # r1 sequence file had > 1 sequence, we need to run each sequence independently through cmalign
        for($r2_i = 0; $r2_i < $cur_nr2; $r2_i++) { 
          push(@r2_seq_file_A, $r1_seq_file_A[$r1_i] . "." . ($r2_i+1));
        }
        $nr2 += $cur_nr2;
      }
    }
  }

  # do all round 2 runs
  if($nr2 > 0) { 
    if($do_glsearch) { # shouldn't happen if $do_glsearch
      ofile_FAIL("ERROR in $sub_name, running glsearch but trying to run stage 2 for at least 1 seq which should only happen if running cmalign", 1, $FH_HR); 
    }
    cmalign_or_glsearch_wrapper_helper($execs_HR, $mdl_file, $mdl_name, $out_root, 2, $extra_key, $nr2, $progress_w, 
                                       \@r2_seq_file_A, \@r2_out_file_AH, \@r2_success_A, \@r2_mxsize_A, 
                                       $opt_HHR, $ofile_info_HHR);
    # go through all round 2 runs: 
    # if it finished successfully record its output files to concatenate later
    # if it did not finish successfully, record the name of the sequence and mxsize required
    for($r2_i = 0; $r2_i < $nr2; $r2_i++) { 
      if($r2_success_A[$r2_i]) { 
        # run finished successfully
        foreach my $out_key (@concat_keys_A) { 
          push(@{$concat_HA{$out_key}}, $r2_out_file_AH[$r2_i]{$out_key});
        }
        push(@{$stk_file_AR}, $r2_out_file_AH[$r2_i]{"stk"});
      }
      else { 
        # run did not finish successfully
        cmalign_store_overflow($r2_seq_file_A[$r2_i], $r2_mxsize_A[$r2_i], $overflow_seq_AR, $overflow_mxsize_AR, $ofile_info_HHR->{"FH"}); 
      }
    }
    # remove sequence files if --keep not used
    if(! opt_Get("--keep", $opt_HHR)) { 
      utl_FileRemoveList(\@r2_seq_file_A, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
    }
  }
    
  # concatenate files into one 
  foreach $out_key (@concat_keys_A) { 
    if(scalar(@{$concat_HA{$out_key}}) > 0) { 
      my $concat_file = sprintf($out_root . ".%s%salign.$out_key", ((defined $mdl_name) ? $mdl_name . "." : ""), $extra_key);                                
      utl_ConcatenateListOfFiles($concat_HA{$out_key}, $concat_file, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
      # utl_ConcatenateListOfFiles() removes individual files unless --keep enabled
      my $out_root_key = sprintf(".concat.%s%salign.$out_key", ((defined $mdl_name) ? $mdl_name . "." : ""), $extra_key);
      ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $out_root_key, $concat_file, 0, $do_keep, sprintf("align $out_key file%s", (defined $mdl_name) ? "for model $mdl_name" : ""));
    }
  }
  # remove sequence files 
  if(($r1_do_split) && (! opt_Get("--keep", $opt_HHR))) { 
    utl_FileRemoveList(\@r1_seq_file_A, $sub_name, $opt_HHR, $ofile_info_HHR->{"FH"});
  }

  return;
}

#################################################################
# Subroutine:  cmalign_or_glsearch_wrapper_helper()
# Incept:      EPN, Wed Mar 20 06:20:51 2019
#
# Purpose:     Run one or more cmalign or glsearch jobs on the farm 
#              or locally. We run glsearch if $mdl_file ends in '.hmm' 
#              else we run cmalign.
#
#              Helper subroutine for cmalign_or_glsearch_wrapper()
#              see that sub's "Purpose" for more details.
#
# Arguments: 
#  $execs_HR:              ref to hash with paths to cmalign, cmsearch and cmfetch
#  $mdl_file:              name of model file to use (if ends with .fa use glsearch else run cmalign)
#  $mdl_name:              name of model to fetch from $mdl_file (undef to not fetch)
#  $out_root:              string for naming output files
#  $round:                 round we are on, "1" or "2"
#  $extra_key:             extra key for output file names, "" or "uj." (latter for seqs with unjoinbl alerts)
#  $nseq:                  total number of sequences in all seq files in @{$seq_file_AR}
#  $progress_w:            width for ofile_OutputProgress* subroutines
#  $seq_file_AR:           ref to array of sequence file names for each cmalign/nhmmscan call, PRE-FILLED
#  $out_file_AHR:          ref to array of hashes of output file names, FILLED HERE 
#  $success_AR:            ref to array of success values
#                          $success_AR->[$j] set to '1' if job finishes successfully
#                                            set to '0' if job fails due to mx overflow (must be cmalign)
#  $mxsize_AR:             ref to array of required matrix sizes
#                          $mxsize_AR->[$j] set to value readh from cmalign output, if $success_AR->[$j] == 0
#                                           else set to '0'
#  $opt_HHR:               REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:        REF to 2D hash of output file information
#
# Returns:     void
# 
# Dies: If an executable doesn't exist, or command fails (and its not a cmalign allowed failure)
#
################################################################# 
sub cmalign_or_glsearch_wrapper_helper { 
  my $sub_name = "cmalign_or_glsearch_wrapper_helper";
  my $nargs_expected = 14;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $mdl_file, $mdl_name, $out_root, $round, $extra_key, $nseq, $progress_w, 
      $seq_file_AR, $out_file_AHR, $success_AR, $mxsize_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $log_FH      = $ofile_info_HHR->{"FH"}{"log"}; # for convenience
  my $do_parallel = opt_Get("-p", $opt_HHR) ? 1 : 0;
  my $do_glsearch = ($mdl_file =~ m/\.fa$/) ? 1 : 0;
  my $nseq_files  = scalar(@{$seq_file_AR});

  # determine description of the runs we are about to do, 
  # depends on $do_parallel, $round, and ($progress_w < 0), and 
  my $stg_desc = "";
  if($do_parallel) { 
    $stg_desc = sprintf("Submitting $nseq_files %s job(s) ($mdl_name: $nseq %sseq%s) to the farm%s", 
                        ($do_glsearch) ? "glsearch" : "cmalign", 
                        ($extra_key eq "uj.") ? "unjoinbl " : "", 
                        ($nseq > 1) ? "s" : "",
                        ($round == 1) ? "" : " to find seqs too divergent to annotate");
  }
  else { 
    $stg_desc = sprintf("Aligning %ssequences ($mdl_name: $nseq seq%s)%s", 
                        ($extra_key eq "uj.") ? "unjoinbl " : "", 
                        ($nseq > 1) ? "s" : "",
                        ($round == 1) ? "" : " to find seqs too divergent to annotate");
  }
  my $start_secs = ofile_OutputProgressPrior($stg_desc, $progress_w, $log_FH, *STDOUT);

  my $key; # a file key
  my $s;   # counter over sequence files
  my @out_keys_A = ("stdout", "err", "ifile", "stk", "sh");
  @{$out_file_AHR} = ();
  for(my $s = 0; $s < $nseq_files; $s++) { 
    %{$out_file_AHR->[$s]} = (); 
    foreach $key (@out_keys_A) { 
      $out_file_AHR->[$s]{$key} = $out_root . "." . $mdl_name . "." . $extra_key . "align.r" . $round . ".s" . $s . "." . $key;
    }
    $success_AR->[$s] = cmalign_or_glsearch_run($execs_HR, $qsub_prefix, $qsub_suffix, 
                                                $mdl_file, $mdl_name, $seq_file_AR->[$s], \%{$out_file_AHR->[$s]},
                                                (defined $mxsize_AR) ? \$mxsize_AR->[$s] : undef, 
                                                $opt_HHR, $ofile_info_HHR);   
    # if we are running parallel, ignore the success return values from the cmalign_or_glsearch_run subroutine
    # vdr_WaitForFarmJobsToFinish() will fill these later
    if($do_parallel) { $success_AR->[$s] = 0; }
  }
  ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

  if($do_parallel) { 
    if((opt_Exists("--align_skip", $opt_HHR)) && (opt_Get("--align_skip", $opt_HHR))) { 
      for($s = 0; $s < $nseq_files; $s++) { 
        $success_AR->[$s] = 1; 
      }
    }
    else { 
      # --align_skip not enabled
      # wait for the jobs to finish
      $start_secs = ofile_OutputProgressPrior(sprintf("Waiting a maximum of %d minutes for all farm jobs to finish", opt_Get("--wait", $opt_HHR)), 
                                              $progress_w, $log_FH, *STDOUT);
      my $njobs_finished = vdr_WaitForFarmJobsToFinish(($do_glsearch ? 0 : 1), # are we are doing cmalign?
                                                       (opt_Get("--errcheck", $opt_HHR)), 
                                                       "stdout", 15, 15, 
                                                       $out_file_AHR,
                                                       $success_AR, 
                                                       $mxsize_AR,  
                                                       ($do_glsearch ? "GLSEARCH" : ""), # value is irrelevant for cmalign
                                                       $opt_HHR, $ofile_info_HHR->{"FH"});
      if($njobs_finished != $nseq_files) { 
        ofile_FAIL(sprintf("ERROR in $sub_name only $njobs_finished of the $nseq_files are finished after %d minutes. Increase wait time limit with --wait", opt_Get("--wait", $opt_HHR)), 1, $ofile_info_HHR->{"FH"});
      }
      ofile_OutputString($log_FH, 1, "# "); # necessary because waitForFarmJobsToFinish() creates lines that summarize wait time and so we need a '#' before 'done' printed by ofile_OutputProgressComplete()
    }

    ofile_OutputProgressComplete($start_secs, undef, $log_FH, *STDOUT);
  }
  
  return;
}

#################################################################
# Subroutine:  cmalign_or_glsearch_run()
# Incept:      EPN, Wed Feb  6 12:30:08 2019
#
# Purpose:     Run Infernal's cmalign or FASTA's glsearch executable 
#              using $mdl_file as the model file on sequence file 
#              $seq_file, either locally or on the farm. We run
#              glsearch if $mdl_file ends in '.hmm' else we run
#              cmalign.
#              
#              If job does not finish successfully, we need to 
#              parse the stderr output (which we redirect to stdout)
#              and see if it failed because of a specific type of
#              error, because the required DP matrix exceeded the
#              size limit. That error looks like this:
#        
#              Error: HMM banded truncated alignment mxes need 60.75 Mb > 2.00 Mb limit.
#              Use --mxsize, --maxtau or --tau.
#              
# Arguments: 
#  $execs_HR:         ref to hash with paths to cmalign, cmfetch and glsearch
#  $qsub_prefix:      qsub command prefix to use when submitting to farm, undef if running locally
#  $qsub_suffix:      qsub command suffix to use when submitting to farm, undef if running locally
#  $mdl_file:         path to the CM or HMM file
#  $mdl_name:         name of model to fetch from $mdl_file (undef to not fetch)
#  $seq_file:         path to the sequence file
#  $out_file_HR:      ref to hash of output files to create
#                     required keys: "stdout", "ifile", "stk", "err", "sh"
#  $ret_mxsize_R:     REF to required matrix size, only filled if return value is '0'
#  $opt_HHR:          REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:   REF to 2D hash of output file information
# 
# Returns:     '1' upon success, which occurs if
#                  job run on farm and submission went ok
#                  job run locally and finished without error
#              '0' upon allowed failure, which occurs if
#                  job run locally and fails because of too big a required matrix
#
# Dies: upon unallowed failure, which occurs if
#                  job run on farm and submission failed
#                  job run locally and finished with unallowed failure
# 
################################################################# 
sub cmalign_or_glsearch_run { 
  my $sub_name = "cmalign_or_glsearch_run()";
  my $nargs_expected = 10;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($execs_HR, $qsub_prefix, $qsub_suffix, 
      $mdl_file, $mdl_name, $seq_file, 
      $out_file_HR, $ret_mxsize_R, $opt_HHR, $ofile_info_HHR) = @_;

  if(defined $ret_mxsize_R) { 
    $$ret_mxsize_R = 0; # overwritten below if nec
  }

  my $FH_HR       = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;
  my $do_parallel = opt_Get("-p", $opt_HHR) ? 1 : 0;
  my $do_glsearch = ($mdl_file =~ m/\.fa$/) ? 1 : 0;
  my $ncpu = opt_Get("--cpu", $opt_HHR);
  if($ncpu == 0) { $ncpu = 1; }

  my $stdout_file = $out_file_HR->{"stdout"};
  my $ifile_file  = $out_file_HR->{"ifile"};
  #my $tfile_file  = $out_file_HR->{"tfile"};
  my $stk_file    = $out_file_HR->{"stk"};
  my $err_file    = $out_file_HR->{"err"};
  my $sh_file     = $out_file_HR->{"sh"};
  if(! defined $stdout_file) { ofile_FAIL("ERROR in $sub_name, stdout output file name is undefined", 1, $FH_HR); }
  if(! defined $err_file)    { ofile_FAIL("ERROR in $sub_name, err    output file name is undefined", 1, $FH_HR); }
  if(! defined $sh_file)     { ofile_FAIL("ERROR in $sub_name, sh     output file name is undefined", 1, $FH_HR); }
  if(! $do_glsearch) { # only need stdout file and ifile file if running glsearch
    if(! defined $stk_file)    { ofile_FAIL("ERROR in $sub_name, stk    output file name is undefined", 1, $FH_HR); }
    if(! defined $ifile_file)  { ofile_FAIL("ERROR in $sub_name, ifile  output file name is undefined", 1, $FH_HR); }
  }
  if((! opt_Exists("--align_skip", $opt_HHR)) || (! opt_Get("--align_skip", $opt_HHR))) { 
    if(-e $stdout_file) { unlink $stdout_file; }
    if(-e $err_file)    { unlink $err_file; }
    if(-e $sh_file)     { unlink $sh_file; }
    if(! $do_glsearch) { 
      if(-e $stk_file)    { unlink $stk_file; }
      if(-e $ifile_file)  { unlink $ifile_file; }
    }
  }
  utl_FileValidateExistsAndNonEmpty($mdl_file, sprintf("%s file", $do_glsearch ? "nucleotide model fasta" : "CM"), $sub_name, 1, $FH_HR); 
  utl_FileValidateExistsAndNonEmpty($seq_file, "sequence file", $sub_name, 1, $FH_HR);

  my $cmd = undef;

  # determine cmalign options based on command line options
  if($do_glsearch) { 
    my $gls_opts = sprintf("-r +%s/%s -f %s -g %s", opt_Get("--gls_match", $opt_HHR), opt_Get("--gls_mismatch", $opt_HHR), opt_Get("--gls_gapopen", $opt_HHR), opt_Get("--gls_gapextend", $opt_HHR));
    $cmd = "cat $seq_file | " . $execs_HR->{"glsearch"} . " $gls_opts -T $ncpu -m 3,9C -z -1 -n -3 -d 1 - $mdl_file > $stdout_file 2>$err_file";
  }
  else { # running cmalign
    my $cmalign_mxsize = sprintf("%.2f", (opt_Get("--mxsize", $opt_HHR) / 4.)); # empirically cmalign can require as much as 4X the amount of memory it thinks it does, this is a problem to fix in infernal
    my $opts = sprintf(" --dnaout --verbose --cpu $ncpu --ifile $ifile_file -o $stk_file --tau %s --mxsize $cmalign_mxsize", opt_Get("--tau", $opt_HHR));
    # add --sub and --notrunc unless --nosub used
    if(! opt_Get("--nosub", $opt_HHR)) { 
      $opts .= " --sub --notrunc"; 
    }
    # add -g unless --noglocal used
    if(! opt_Get("--noglocal", $opt_HHR)) { 
      $opts .= " -g"; 
    }
    if(! opt_Get("--nofixedtau", $opt_HHR)) { 
      $opts .= " --fixedtau"; 
    }
    if(defined $mdl_name) { 
      $cmd = $execs_HR->{"cmfetch"} . " $mdl_file $mdl_name | " . $execs_HR->{"cmalign"} . " $opts - $seq_file > $stdout_file 2>&1";
    }
    else { 
      $cmd = $execs_HR->{"cmalign"} . " $opts $mdl_file $seq_file > $stdout_file 2>&1";
    }
  }

  my $success = 1;
  if($do_parallel) { 
    my $job_name = "J" . utl_RemoveDirPath($seq_file);
    my $nsecs  = opt_Get("--wait", $opt_HHR) * 60.;
    my $mem_gb = opt_Get("--mxsize", $opt_HHR) / 1000.;
    if($mem_gb < 16.) { $mem_gb = 16.; } # set minimum of 16 Gb
    if((! opt_Exists("--align_skip", $opt_HHR)) || (! opt_Get("--align_skip", $opt_HHR))) { 
      vdr_SubmitJobAsScript($cmd, $qsub_prefix, $qsub_suffix, $job_name, $sh_file, $err_file, $mem_gb, $nsecs, $opt_HHR, $ofile_info_HHR);
    }
  }
  else { 
    if((! opt_Exists("--align_skip", $opt_HHR)) || (! opt_Get("--align_skip", $opt_HHR))) { 
      utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 1, $FH_HR); # 1 says: it's okay if job fails
    }
    # command has completed
    if($do_glsearch) { # glsearch: final line of stdout file should have 'GLSEARCH' in it
      my $final_line = `tail -n 1 $stdout_file`;
      chomp $final_line;
      if($final_line =~ m/\r$/) { chop $final_line; } # remove ^M if it exists
      $success = ($final_line =~ m/GLSEARCH/) ? 1 : 0;
    }
    else { # cmalign: check for the error in the stdout, or a final line of 'CPU' indicating that it worked.
      $success = vdr_CmalignCheckStdOutput($stdout_file, $ret_mxsize_R, $FH_HR);
      if($success == -1) { # indicates job did not finish properly, this shouldn't happen because utl_RunCommand() didn't die
        ofile_FAIL("ERROR in $sub_name, cmalign failed in a bad way, see $stdout_file for error output", 1, $ofile_info_HHR->{"FH"});
      }
    }
  }
  
  return $success; 
}

#################################################################
# Subroutine : parse_stk_and_add_alignment_cds_and_mp_alerts()
# Incept:      EPN, Thu Jan 31 13:06:54 2019
#
# Purpose:    Parse Infernal 1.1 cmalign stockholm alignment file
#             and store results in @{$mdl_results_AAHR}. 
# 
#             Added post v1.1.3 with --glsearch was added to deal with
#             case where some SARS-CoV-2 seqs were failing due to a
#             indel near start/stop codon: Potentially doctors (modifies) the
#             stockholm alignment for each sequence if (and this
#             should be rare) a gap exists at first position of start
#             or final position of stop codon, or single insertion occurs
#             after first stop position or before final stop position.
#             Only doctors if doctoring will create a valid start/stop.
#             See toy examples in comments at beginning of 
#             doctoring_check_new_codon_validity() subroutine.
#
#             Also, potentially re-doctors already doctored aligned
#             sequence exactly once if the initial doctoring disrupted
#             a different valid start/stop.
#
#             Detects and adds the following alerts to 
#             @{$alt_ftr_instances_HHHR}:
#             indf5gap: gap at 5' boundary of model span for a feature segment
#             indf3gap: gap at 5' boundary of model span for a feature segment
#             indf5lcc: low posterior prob at 5' boundary of model span for a coding feature segment
#             indf5lcn: low posterior prob at 5' boundary of model span for a noncoding feature segment
#             indf3lcc: low posterior prob at 3' boundary of model span for a coding feature segment
#             indf3lcn: low posterior prob at 3' boundary of model span for a noncoding feature segment
#
# Arguments: 
#  $stk_file:                  stockholm alignment file to parse
#  $in_sqfile_R:               REF to Bio::Easel::SqFile object from input fasta file, can be undef unless --alicheck used
#  $mdl_tt:                    the translation table ('1' for standard)
#  $seq_len_HR:                REF to hash of sequence lengths, PRE-FILLED
#  $seq_inserts_HHR:           REF to hash of hashes with sequence insert information, PRE-FILLED
#  $sgm_info_AHR:              REF to array of hashes with information on the model segments, PRE-FILLED
#  $ftr_info_AHR:              REF to array of hashes with information on the features, PRE-FILLED
#  $alt_info_HHR:              REF to hash of hashes with information on the errors, PRE-FILLED
#  $stg_results_HHHR:          REF to 3D hash of classification results, PRE-FILLED
#  $sgm_results_HAHR:          REF to results HAH, FILLED HERE
#  $ftr_results_HAHR:          REF to feature results HAH, possibly ADDED TO HERE
#  $alt_seq_instances_HHR:     REF to array of hash with per-sequence alerts, ADDED TO HERE
#  $alt_ftr_instances_HHHR:    REF to error instances HHH, ADDED TO HERE
#  $dcr_output_HAHR:           REF to hash of array of hashes with info on doctored seqs to output, ADDED TO HERE
#  $mdl_name:                  model name this alignment pertains to
#  $ftr_fileroot_AR:           REF to array of per-feature file root values, pre-calc'ed and passed in so we don't need to do it per-seq
#  $ftr_outroot_AR:            REF to array of per-feature output root values, pre-calc'ed and passed in so we don't need to do it per-seq
#  $sqfile_for_cds_mp_alerts:  REF to Bio::Easel::SqFile object, open sequence file with sequences
#                              to fetch CDS and mat_peptides from to analyze for possible alerts 
#  $sqfile_for_output_fastas:  REF to Bio::Easel::SqFile object, open sequence file with sequences
#                              that we'll fetch feature sequences from to output to per-feature fasta files
#  $sqfile_for_pv:             REF to Bio::Easel::SqFile object, open sequence file with sequences
#                              that we'll fetch feature sequences from for protein validation
#  $do_separate_cds_fa_files:  '1' to output two sets of cds files, one with fetched features from $sqfile_for_output_fastas
#                              and one for the protein validation stage fetched from $sqfile_for_cds_mp_alerts
#  $to_remove_AR:              REF to array of files to remove before exiting, possibly added to here if $do_separate_cds_fa_files
#  $out_root:                  string for naming output files
#  $opt_HHR:                   REF to 2D hash of option values
#  $ofile_info_HHR:            REF to 2D hash of output file information
#
# Returns:    void
#
# Dies:
#
################################################################# 
sub parse_stk_and_add_alignment_cds_and_mp_alerts { 
  my $sub_name = "parse_stk_and_add_alignment_cds_and_mp_alerts()";
  my $nargs_exp = 25;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($stk_file, $in_sqfile_R, $mdl_tt, $seq_len_HR, $seq_inserts_HHR, $sgm_info_AHR, 
      $ftr_info_AHR, $alt_info_HHR, $stg_results_HHHR, $sgm_results_HAHR, $ftr_results_HAHR, 
      $alt_seq_instances_HHR, $alt_ftr_instances_HHHR, $dcr_output_HAHR, 
      $mdl_name, $ftr_fileroot_AR, $ftr_outroot_AR, 
      $sqfile_for_cds_mp_alerts, $sqfile_for_output_fastas, $sqfile_for_pv,
      $do_separate_cds_fa_files, $to_remove_AR, $out_root, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = \%{$ofile_info_HHR->{"FH"}};
  my $pp_thresh_non_mp = opt_Get("--indefann",    $opt_HHR); # threshold for non-mat_peptide features
  my $pp_thresh_mp     = opt_Get("--indefann_mp", $opt_HHR); # threshold for mat_peptide features
  my $do_alicheck      = opt_Get("--alicheck",    $opt_HHR); # check aligned sequences are identical to those fetched from $sqfile (except maybe Ns if -r) 
  my $do_replace_ns    = opt_Get("-r",            $opt_HHR); # only relevant if $do_alicheck
  my $do_glsearch      = opt_Get("--glsearch",    $opt_HHR); # we won't have PP values if this is 1
  my $do_nodcr         = opt_Get("--nodcr",       $opt_HHR); # do not doctor alignment to correct start/stop codons
  my $do_forcedcrins   = opt_Get("--forcedcrins", $opt_HHR); # force doctoring of insert type, must be 1 seq per alignment (--forcedcrins requires --cmindi)
  my $small_value = 0.000001; # for checking if PPs are below threshold
  my $nftr = scalar(@{$ftr_info_AHR});
  my $nsgm = scalar(@{$sgm_info_AHR});

  # create an ESL_MSA object from the alignment
  # open and validate file
  my $msa = Bio::Easel::MSA->new({
    fileLocation => $stk_file,
    isDna => 1});  

  if(($do_alicheck) && (! defined $in_sqfile_R)) { 
    ofile_FAIL("ERROR in $sub_name, --alicheck used but no sqfile provided", 1, $FH_HR);
  }

  # build a map of aligned positions to model RF positions and vice versa, only need to do this once per alignment
  my @rf2a_A = (); # [1..$rfpos..$rflen] = $apos;  rf position $rfpos maps to alignment position $apos [1..$alen]  ($rf2a_A[0] = -1  (dummy value))
  my $rflen = msa_create_rfpos_to_apos_map($msa, \@rf2a_A, $FH_HR);
  my $rfpos = 0; # nongap RF (model) position [1..$rflen]
  my $apos  = 0; # alignment position [1..$alen]
  my $alen  = $msa->alen;

  # 'doctor' flags, which keep track of when we need to doctor the alignment in an attempt
  # to fix indels in first position of CDS starts and final positions of CDS stops
  # We can doctor each sequence up to twice (second doctoring will actually revert previous one)
  # but no more, else we'll enter an infinite loop (relevant comments below include string 'infinite')
  my $seq_doctor_ctr   = 0; # incremented when if we need to doctor the sequence because it had 
                            # an indel in first/final position of a CDS, if this is going to 
                            # exceed 2 we stop doctoring
  my $seq_doctor_flag  = 0; # set to 1 if we should doctor the current sequence
  my $msa_doctor_flag  = 0; # set to 1 if we end up doctoring any sequence, if 1 at end
                            # we have to rewrite the stockholm MSA file to save doctored changes

  # move through each sequence in the alignment and determine its boundaries for each model region
  my $nseq = $msa->nseq; 
  if($do_glsearch && ($nseq != 1)) { 
    ofile_FAIL("ERROR in $sub_name, --glsearch enabled but $nseq > 1 seqs in alignment for parsing", 1, $FH_HR);
  }
  if($do_forcedcrins && ($nseq != 1)) { 
    ofile_FAIL("ERROR in $sub_name, --forcedcrins enabled but $nseq > 1 seqs in alignment for parsing", 1, $FH_HR);
  }

  # for each sequence, go through all segments and fill in the start and stop (unaligned seq) positions
  for(my $i = 0; $i < $nseq; $i++) { 
    my $seq_name = $msa->get_sqname($i);
    if(! exists $seq_len_HR->{$seq_name}) { 
      ofile_FAIL("ERROR in $sub_name, do not have length information for sequence $seq_name from alignment in $stk_file", 1, $FH_HR);
    }
    my $seq_len = $seq_len_HR->{$seq_name};
    if(! defined $seq_inserts_HHR->{$seq_name}{"ins"}) { 
      ofile_FAIL("ERROR in $sub_name, do not have insert information for sequence $seq_name from alignment in $stk_file", 1, $FH_HR);
    }
    my $seq_ins = $seq_inserts_HHR->{$seq_name}{"ins"}; # string of inserts
    my @do_dcr_idx_A = (); # array of indices in $dcr_output_HAHR->{$seq_name} that we will actually doctor alignment for
    $seq_doctor_flag = 0;

    @{$sgm_results_HAHR->{$seq_name}} = ();
    
    # arrays that hold per-alert info that we defer until after the 'for(sgm_idx...' block 
    # just in case we have to doctor the alignment and reevaluate the sequence
    # (we don't want to have reported any alerts for a seq we are going to reevaluate after doctoring)
    my @alt_code_A    = (); # array of alert codes to add for this sequence after for(sgm... block
    my @alt_str_A     = (); # array of alert strings to add for this sequence after for(sgm... block
    my @alt_scoords_A = (); # array of alert sequence coords strings to add for this sequence after for(sgm... block
    my @alt_mcoords_A = (); # array of alert model coords strings to add for this sequence after for(sgm... block
    my @alt_ftr_A     = (); # array of alert ftr_idx to add for this sequence after for(sgm... block

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
    
    if($seq_ins ne "") { 
      my @ins_A = split(";", $seq_inserts_HHR->{$seq_name}{"ins"});
      foreach my $ins_tok (@ins_A) { 
        if($ins_tok =~ /^(\d+)\:(\d+)\:(\d+)$/) { 
          my ($i_rfpos, $i_uapos, $i_len) = ($1, $2, $3);
          $rf2ipos_A[$i_rfpos] = $i_uapos;
          $rf2ilen_A[$i_rfpos] = $i_len;
          #printf("rf2ipos_A[%5d]: %5d rf2ilen_A[%5d]: %5d\n", $i_rfpos, $i_uapos, $i_rfpos, $i_len);
        }
        else { 
          ofile_FAIL("ERROR in $sub_name, failed to parse insert information read from ifile for $seq_name:\n" . $seq_inserts_HHR->{$seq_name}{"ifile_ins"}, 1, $FH_HR);
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
                                 #                       or '?' if $do_glsearch (--glsearch)
                                 #                       '.' if sequence is a gap at that RF position $rfpos (even if $do_glsearch)
                                 #                       special values: $rfpos_pp_A[0] = -1, $rfpos_pp_A[$rflen+1] = -1
                                 # 

    my @ua2rf_A = ();            # [1..$uapos..$ualen]: reference position that unaligned position $uapos aligns to 
                                 #                      if $ua2rf_A[$uapos] <  0, $uapos inserts *after* ref posn (-1 * $ua2rf_A[$uapos])
                                 #                      if $ua2rf_A[$uapos] == 0, $uapos inserts *before* ref posn 1
                                 #                      $ua2rf_A[0] is invalid (set to 0)

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
    $ua2rf_A[0] = 0; # invalid
    # get aligned sequence, length will be alen
    my $sqstring_aligned = $msa->get_sqstring_aligned($i);
    my $ppstring_aligned = ($do_glsearch) ? undef : $msa->get_ppstring_aligned($i);
    if(length($sqstring_aligned) != $alen) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, fetched aligned seqstring of unexpected length (%d, not %d)\n$sqstring_aligned\n", length($sqstring_aligned), $alen), 1, $FH_HR);
    }
    if((! $do_glsearch) && (length($ppstring_aligned) != $alen)) { 
      ofile_FAIL(sprintf("ERROR in $sub_name, fetched aligned posterior probability string of unexpected length (%d, not %d)\n$sqstring_aligned\n", length($ppstring_aligned), $alen), 1, $FH_HR);
    }

    # check de-aligned sequence is identical to input sequence fetched from $in_sqfile IF --alicheck
    if($do_alicheck) { 
      my $in_sqstring = $$in_sqfile_R->fetch_seq_to_sqstring($seq_name);
      my $ua_sqstring = $sqstring_aligned;
      $ua_sqstring =~ s/\W//g;
      my $ua_len = length($ua_sqstring);
      if($ua_len != length($in_sqstring)) { 
        ofile_FAIL(sprintf("ERROR in $sub_name, checking aligned vs input sequences due to --alicheck, lengths for $seq_name do not match: input: %d, aligned: %d\n", length($in_sqstring), $ua_len), 1, $FH_HR);
      }
      $in_sqstring =~ tr/a-z/A-Z/; # uppercase-ize
      $ua_sqstring =~ tr/a-z/A-Z/; # uppercase-ize
      my @ua_sqstring_A = split("", $ua_sqstring);
      my @in_sqstring_A = split("", $in_sqstring);
      for(my $z = 0; $z < $ua_len; $z++) { 
        if(($in_sqstring_A[$z] ne $ua_sqstring_A[$z]) && 
           ((! $do_replace_ns) || ($in_sqstring_A[$z] ne "N"))) { 
          ofile_FAIL(sprintf("ERROR in $sub_name, checking aligned vs input sequences due to --alicheck, seq $seq_name position %d different%s, input: %s, aligned: %s\n", 
                             ($z+1), ($do_replace_ns ? " and not N" : ""), $in_sqstring_A[$z], $ua_sqstring_A[$z]), 
                     1, $FH_HR);
        }
      }
    }

    my @sq_A = split("", $sqstring_aligned);
    my @pp_A = ();
    if(! $do_glsearch) { 
      @pp_A = split("", $ppstring_aligned);
    }

    # first pass, from right to left to fill $min_**pos_after arrays, and rf
    my $min_rfpos = -1;
    my $min_uapos = $seq_len+1;
    for($rfpos = $rflen; $rfpos >= 0; $rfpos--) { 
      $apos = $rf2a_A[$rfpos];
      my $nongap_rf    = (($rfpos > 0) && ($sq_A[($apos-1)] ne ".") && ($sq_A[($apos-1)] ne "-")) ? 1 : 0;
      my $insert_after = ($rf2ipos_A[$rfpos] != -1) ? 1 : 0;
      if($nongap_rf || $insert_after) { 
        $min_rfpos = $rfpos;
      }
      if($insert_after) { 
        $min_uapos -= $rf2ilen_A[$rfpos]; # subtract inserts between $rfpos and ($rfpos+1)
        for(my $tmp_uapos = $min_uapos; $tmp_uapos < ($min_uapos + $rf2ilen_A[$rfpos]); $tmp_uapos++) { 
          $ua2rf_A[$tmp_uapos] = -1 * $rfpos; # note if $rfpos is 0, this will be 0
        }
      }
      if($nongap_rf) { 
        $min_uapos--;
        $ua2rf_A[$min_uapos] = $rfpos;
        $rfpos_pp_A[$rfpos] = ($do_glsearch) ? "?" : $pp_A[($apos-1)];
      }
      $min_rfpos_after_A[$rfpos] = $min_rfpos;
      $min_uapos_after_A[$rfpos] = $min_uapos;
      #printf("rfpos: %5d  apos: %5d  min_rfpos: %5d  min_uapos: %5d\n", $rfpos, $apos, $min_rfpos, $min_uapos);
    }
    if($min_uapos != 1) { 
      ofile_FAIL("ERROR in $sub_name, failed to account for all nucleotides when parsing alignment for $seq_name, pass 1 (min_uapos should be 1 but it is $min_uapos)", 1, $FH_HR);
    }      

    # second pass, from left to right to fill $max_**pos_before arrays:
    my $max_rfpos = -1;
    my $max_uapos = 0;
    for($rfpos = 1; $rfpos <= ($rflen+1); $rfpos++) { 
      $apos = $rf2a_A[$rfpos];
      my $nongap_rf     = (($rfpos <= $rflen) && ($sq_A[($apos-1)] ne ".") && ($sq_A[($apos-1)] ne "-")) ? 1 : 0;
      my $insert_before = ($rf2ipos_A[($rfpos-1)] != -1) ? 1 : 0;
      if($nongap_rf || $insert_before) { 
        $max_rfpos = $rfpos;
      }
      if($insert_before) { 
        $max_uapos += $rf2ilen_A[($rfpos-1)]; # subtract inserts between ($rfpos-1) and $rfpos
      }
      if($nongap_rf) { 
        $max_uapos++;
        $rfpos_pp_A[$rfpos] = ($do_glsearch) ? "?" : $pp_A[($apos-1)];
      }
      $max_rfpos_before_A[$rfpos] = $max_rfpos;
      $max_uapos_before_A[$rfpos] = $max_uapos;

      #if($rfpos <= $rflen) { 
      #  printf("rfpos: %5d  apos: %5d  max_rfpos: %5d  max_uapos: %5d\n", $rfpos, $apos, $max_rfpos, $max_uapos);
      #}
    }
    if($max_uapos != $seq_len) { 
      ofile_FAIL("ERROR in $sub_name, failed to account for all nucleotides when parsing alignment for $seq_name, pass 2 (max_uapos should be $seq_len but it is $max_uapos)", 1, $FH_HR);
    }      

    # given model span s..e
    # if strand eq "+"
    #   if C[rfpos] > D[rfpos] then no hit (A[rfpos] should be > B[rfpos])
    #   else (C[rfpos] <= D[rfpos]) 
    #        hit uaseq span is from C[rfpos] to D[rfpos]
    #        hit rf span is from A[rfpos] to B[rfpos]

    # now we have all the info we need for this sequence to determine sequence boundaries for each model segment
    my $sgm_idx; 
    my $ftr_idx;
    my %ftr_deletinf_alt_msg_HA = (); # key is $ftr_idx, value is an array of alert messages for deletinf alerts, one per 
                                      # segment for $ftr_idx that is completely deleted. This is rare and *not identifying*
                                      # these is github issue 21
    for($sgm_idx = 0; $sgm_idx < $nsgm; $sgm_idx++) { 
      my $sgm_start_rfpos = $sgm_info_AHR->[$sgm_idx]{"start"};
      my $sgm_stop_rfpos  = $sgm_info_AHR->[$sgm_idx]{"stop"};
      my $sgm_strand      = $sgm_info_AHR->[$sgm_idx]{"strand"};
      $ftr_idx = $sgm_info_AHR->[$sgm_idx]{"map_ftr"};
      my $ftr_pp_thresh = (vdr_FeatureTypeIsMatPeptide($ftr_info_AHR, $ftr_idx)) ? $pp_thresh_mp : $pp_thresh_non_mp;
      my $ftr_pp_msg    = (vdr_FeatureTypeIsMatPeptide($ftr_info_AHR, $ftr_idx)) ? " (mat_peptide feature)" : "";

#####################################
# Debugging print block
#      printf("segment $sgm_idx $sgm_start_rfpos..$sgm_stop_rfpos\n");
#      $rfpos = $sgm_start_rfpos;
#      printf("\trfpos[%5d] min_rf_after_A: %5d  min_ua_after_A: %5d  max_rf_before_A: %5d  max_ua_before_A: %5d\n", 
#             $rfpos, 
#             $min_rfpos_after_A[$rfpos],
#             $min_uapos_after_A[$rfpos],
#             $max_rfpos_before_A[$rfpos],
#             $max_uapos_before_A[$rfpos]);
#      $rfpos = $sgm_stop_rfpos;
#      printf("\trfpos[%5d] min_rf_after_A: %5d  min_ua_after_A: %5d  max_rf_before_A: %5d  max_ua_before_A: %5d\n", 
#             $rfpos, 
#             $min_rfpos_after_A[$rfpos],
#             $min_uapos_after_A[$rfpos],
#             $max_rfpos_before_A[$rfpos],
#             $max_uapos_before_A[$rfpos]);
#####################################

      my $start_rfpos = -1; # model position of start of this model region for this aligned sequence, stays at -1 if none
      my $stop_rfpos  = -1; # model position of stop  of this model region for this aligned sequence, stays at -1 if none
      my $start_uapos = -1; # unaligned position of start of this model region for this aligned sequence, stays at -1 if none
      my $stop_uapos  = -1; # unaligned position of stop  of this model region for this aligned sequence, stays at -1 if none
      my $p_5seqflush = undef;
      my $p_3seqflush = undef;

      # determine start and stop position differently based on strand
      $start_rfpos = ($sgm_strand eq "+") ? $min_rfpos_after_A[$sgm_start_rfpos] : $max_rfpos_before_A[$sgm_start_rfpos];
      $stop_rfpos  = ($sgm_strand eq "+") ? $max_rfpos_before_A[$sgm_stop_rfpos] : $min_rfpos_after_A[$sgm_stop_rfpos];

      $start_uapos = ($sgm_strand eq "+") ? $min_uapos_after_A[$sgm_start_rfpos] : $max_uapos_before_A[$sgm_start_rfpos];
      $stop_uapos  = ($sgm_strand eq "+") ? $max_uapos_before_A[$sgm_stop_rfpos] : $min_uapos_after_A[$sgm_stop_rfpos];

      # determine if we have a valid annotation for this segment
      my $is_valid = 1; # assume we do, and check for 3 cases in which we don't below

      if(($start_rfpos == -1) || ($stop_rfpos == -1)) { 
        # alignment doesn't span segment RF positions
        $is_valid = 0; 
      }
      elsif((($sgm_strand eq "+") && ($start_rfpos > $stop_rfpos)) || # complete segment deleted on + strand
            (($sgm_strand eq "-") && ($start_rfpos < $stop_rfpos))) { # complete segment deleted on - strand
        $is_valid = 0; 
        # keep track that this segment is completely deleted, by constructing its 
        # alert message for a possible deletinf alert. However we can't report it yet
        # because if all segments for this feature are deleted we will report a 
        # deletins or deletina (per-sequence) alert instead. So we just store the possible
        # deletinf alert here in %ftr_deletinf_alt_msg_HA and then deal with it after
        # the 'for($sgm_idx=0..$nsgm-1)' block below
        my $ftr_nsgm = ($ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"} - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}) + 1;
        my $alt_scoords = "seq:VADRNULL;"; # feature is deleted, no sequence info available
        my $alt_mcoords = "mdl:" . vdr_CoordsSegmentCreate($sgm_start_rfpos, $sgm_stop_rfpos, $sgm_strand, $FH_HR) . ";"; 
        my $alt_msg = $alt_scoords . $alt_mcoords;
        $alt_msg .= ($ftr_nsgm > 1) ? 
            sprintf("segment %d of %d deleted", ($sgm_idx - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}) + 1, $ftr_nsgm) : 
            "complete single segment feature deleted";
        if(! defined $ftr_deletinf_alt_msg_HA{$ftr_idx}) { 
          @{$ftr_deletinf_alt_msg_HA{$ftr_idx}} = (); 
        }
        push(@{$ftr_deletinf_alt_msg_HA{$ftr_idx}}, $alt_msg);
      }  

      if($is_valid) { 
        if($sgm_strand eq "+") { 
          $p_5seqflush = ($start_uapos == 1)        ? 1 : 0;
          $p_3seqflush = ($stop_uapos  == $seq_len) ? 1 : 0;
        }
        else { 
          $p_5seqflush = ($start_uapos == $seq_len) ? 1 : 0;
          $p_3seqflush = ($stop_uapos  == 1)        ? 1 : 0;
        }

        %{$sgm_results_HAHR->{$seq_name}[$sgm_idx]} = ();
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"}    = $start_uapos;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstop"}     = $stop_uapos;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"mstart"}    = $start_rfpos;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"mstop"}     = $stop_rfpos;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"}    = $sgm_strand;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"5seqflush"} = $p_5seqflush;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"3seqflush"} = $p_3seqflush;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"5trunc"}    = ($p_5seqflush && ($start_rfpos != $sgm_start_rfpos)) ? 1 : 0;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"3trunc"}    = ($p_3seqflush && ($stop_rfpos  != $sgm_stop_rfpos))  ? 1 : 0;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"startgap"}  = ($rfpos_pp_A[$sgm_start_rfpos] eq ".") ? 1  : 0;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"stopgap"}   = ($rfpos_pp_A[$sgm_stop_rfpos]  eq ".") ? 1  : 0;
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"startpp"}   = ($rfpos_pp_A[$sgm_start_rfpos] eq ".") ? -1 : ($do_glsearch ? "?" : convert_pp_char_to_pp_avg($rfpos_pp_A[$sgm_start_rfpos], $FH_HR));
        $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"stoppp"}    = ($rfpos_pp_A[$sgm_stop_rfpos]  eq ".") ? -1 : ($do_glsearch ? "?" : convert_pp_char_to_pp_avg($rfpos_pp_A[$sgm_stop_rfpos], $FH_HR));

        # Check for special case where CDS starts/stops with a gap
        # or has a single insert after first position of a start,
        # or before final position of a stop.
        # If so, we may actually doctor the alignment and then 
        # rerun the main loop over segments by setting $seq_doctor_flag to 1
        #
        # Start/stop codons that start/end with a gap are referred to as
        # 'delete' type doctorings in the code.
        # Start/stop codons that have a single insertion in them
        # are referred to as 'insert' type doctorings in the code.
        #
        # We deal with each type differently. 
        #
        # Both types of doctorings only actually take place if doing them
        # will lead to a valid start/stop codon.
        # 
        # Delete type notes:
        # To fix a start or stop codon that starts/ends with a gap, 
        # we will doctor the alignment by swapping gap with closest nt
        # in proper 5'/3' direction (depending on start/stop and strand)
        # if: 
        # - first/final RF position of segment aligns to a gap
        # - this is first/final segment of a CDS
        # - there exists another nt to swap with
        # - the nt we will swap with is not an insert (this would invalidate insert info)
        # - the swap will result in a valid start/stop codon, taking
        #   strand into account
        # 
        # Insert type notes:
        # To fix a start or stop codon that has a single insert, 
        # we will doctor the alignment by swapping the adjacent non-gap 
        # RF position with gap RF position in the proper 5'/3' direction 
        # (depending on start/stop and strand)
        # if: 
        # - single nt insertion after first start RF position or
        #   single nt insertion before final stop RF position
        # - this is first/final segment of a CDS
        # - there exists another RF position to swap with
        # - the swap will result in a valid start/stop codon, taking
        #   strand into account
        # 
        # See 8 toy examples of doctorings for the 8 possible combinations of
        # delete/insert type, start/stop, +/- strand 
        # in the comments of the doctoring_check_new_codon_validity() subroutine.
        #
        # doctoring_check_new_codon_validity() does most of the work. It 
        # stores the information on what positions to doctor in %{$dcr_output_HAHR}
        # and then we actually do the doctoring after the 'for(sgm..' loop when all 
        # such doctorings have been collected.
        # 
        # NOTE: there are some situations involving multiple gaps where this will likely
        # not fix the problem and you'll still get an error when a valid start/stop exists,
        # possibly with non-standard translation tables.
        #
        my $dcr_delete_or_insert = ""; # set to 'delete' or 'insert' if we find doctoring possibility of type 'delete' or 'insert'
        my $do_dcr_idx; 
        if(! $do_nodcr) { # if --nodcr we never doctor
          # printf("sgm_start_rfpos: $sgm_start_rfpos, rf2ilen_A[$sgm_start_rfpos] $rf2ilen_A[$sgm_start_rfpos]\n");
          # check for gap at start of start codon that we can try to fix (delete type doctoring)
          # or insert near start position that we can try to fix (insert type doctoring)
          if((vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx) && ($sgm_info_AHR->[$sgm_idx]{"is_5p"})) &&  # this is first segment of a CDS
             ((($sgm_strand eq "+") && ($start_uapos > 1)) || (($sgm_strand eq "-") && ($start_uapos < $seq_len)))) { # we have an nt to swap with
            $dcr_delete_or_insert = ""; # set to 'insert' or 'delete' if we identify an insertion or deletion type doctoring below
            ###############################################################
            # check for two cases where we would doctor the alignment to fix a start: 
            # 'delete' type: gap in first RF position of start (deletion)
            # 'insert' type: insert after first RF position of start (insertion)
            if(($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"startgap"}) && # first RF position of segment aligns to a gap
               ((($sgm_strand eq "+") && ($rf2ilen_A[($sgm_start_rfpos-1)] == -1)) ||  # + strand: we won't be swapping with an insert
                (($sgm_strand eq "-") && ($rf2ilen_A[($sgm_start_rfpos)]   == -1)))) { # - strand: we won't be swapping with an insert
              $dcr_delete_or_insert = "delete";
            }
            elsif((($sgm_strand eq "+") && ($rf2ilen_A[($sgm_start_rfpos)]   == 1)) || # + strand: insert after first RF start position
                  (($sgm_strand eq "-") && ($rf2ilen_A[($sgm_start_rfpos-1)] == 1))) { # - strand: insert before first RF start position
              $dcr_delete_or_insert = "insert";
            }
            if(($dcr_delete_or_insert eq "delete") || 
               (($dcr_delete_or_insert eq "insert") && ($do_glsearch || $do_forcedcrins))) { 
              $do_dcr_idx = doctoring_check_new_codon_validity($dcr_delete_or_insert, "start", $sgm_strand, 
                                                               $seq_name, $mdl_name, $mdl_tt, $ftr_idx, 
                                                               $start_uapos, $sgm_start_rfpos, $sqstring_aligned, 
                                                               $seq_doctor_ctr, \@rf2a_A, $dcr_output_HAHR, 
                                                               $opt_HHR, $FH_HR);
              if($do_dcr_idx != -1) { 
                push(@do_dcr_idx_A, $do_dcr_idx);
                $seq_doctor_flag = 1;
                $msa_doctor_flag = 1;
              }
            }
          }
          # check for gap at end of stop codon that we can try to fix (delete type doctoring)
          # or insert near stop position that we can try to fix (insert type doctoring)
          if((vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx) && ($sgm_info_AHR->[$sgm_idx]{"is_3p"})) && # this is final segment of a CDS
             ((($sgm_strand eq "+") && ($stop_uapos < $seq_len)) || (($sgm_strand eq "-") && ($stop_uapos > 1)))) {  # we have an nt to swap with
            $dcr_delete_or_insert = ""; # set to 'insert' or 'delete' if we identify an insertion or deletion type doctoring below
            ###############################################################
            # check for two cases where we would doctor the alignment to fix a stop:
            # 'delete' type: gap in final RF position of stop (deletion)
            # 'insert' type: insert before final RF position of stop (insertion)
            if(($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"stopgap"}) && # final RF position of segment aligns to a gap
                ((($sgm_strand eq "+") && ($rf2ilen_A[$sgm_stop_rfpos]     == -1)) ||  # + strand: we won't be swapping with an insert
                 (($sgm_strand eq "-") && ($rf2ilen_A[($sgm_stop_rfpos-1)] == -1)))) { # - strand: we won't be swapping with an insert
              $dcr_delete_or_insert = "delete";
            }
            elsif((($sgm_strand eq "+") && ($rf2ilen_A[($sgm_stop_rfpos-1)] == 1))  || # + strand: insert before final RF stop position
                  (($sgm_strand eq "-") && ($rf2ilen_A[($sgm_stop_rfpos)]   == 1))) {  # - strand: insert after final RF stop position
              $dcr_delete_or_insert = "insert";
            }
            if(($dcr_delete_or_insert eq "delete") || 
               (($dcr_delete_or_insert eq "insert") && ($do_glsearch || $do_forcedcrins))) { 
              $do_dcr_idx = doctoring_check_new_codon_validity($dcr_delete_or_insert, "stop", $sgm_strand, 
                                                               $seq_name, $mdl_name, $mdl_tt, $ftr_idx, 
                                                               $stop_uapos, $sgm_stop_rfpos, $sqstring_aligned, 
                                                               $seq_doctor_ctr, \@rf2a_A, $dcr_output_HAHR, 
                                                               $opt_HHR, $FH_HR);
              if($do_dcr_idx != -1) { 
                push(@do_dcr_idx_A, $do_dcr_idx);
                $seq_doctor_flag = 1;
                $msa_doctor_flag = 1;
              }
            }
          }
        } # end of 'if(! $do_nodcr)'

        # store info on alerts we will report later, if nec
        if(! $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"5trunc"}) { 
          if($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"startgap"}) { 
            push(@alt_code_A,    "indf5gap");
            my $tmp_alt_str = vdr_FeatureSummarizeSegment($ftr_info_AHR, $sgm_info_AHR, $sgm_idx);
            if($tmp_alt_str eq "") { 
              $tmp_alt_str = "VADRNULL";
            }
            else { 
              $tmp_alt_str =~ s/^\,\s+//; # ", segment 1 of 2" -> "segment 1 of 2"
            }
            push(@alt_str_A,     $tmp_alt_str);
            push(@alt_scoords_A, sprintf("seq:%s;", vdr_CoordsSinglePositionSegmentCreate($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"}, 
                                                                                          $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"}, $FH_HR)));
            push(@alt_mcoords_A, sprintf("mdl:%s;", vdr_CoordsSinglePositionSegmentCreate($sgm_start_rfpos, $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"}, $FH_HR)));
            push(@alt_ftr_A,     $ftr_idx);
          } 
          elsif((! $do_glsearch) && (($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"startpp"} - $ftr_pp_thresh) < (-1 * $small_value))) { # only check PP if it's not a gap
            # report indf5lcc or indf5lcn
            # indf5lcc: if this segment is 5'-most segment of a CDS or 5'-most segment of a feature that has same start position as a CDS
            # indf5lcn: if not indf5lcc
            if(($sgm_info_AHR->[$sgm_idx]{"is_5p"}) && # segment is 5'-most segment in feature
               ((vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) || # feature is a CDS
                (vdr_SegmentStartIdenticalToCds($ftr_info_AHR, $sgm_info_AHR, $sgm_idx, $FH_HR)))) { # segment has start same as a CDS start
              push(@alt_code_A, "indf5lcc");
            }
            else { 
              push(@alt_code_A, "indf5lcn");
            }
            push(@alt_str_A,     sprintf("%.2f<%.2f%s%s", $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"startpp"}, $ftr_pp_thresh, $ftr_pp_msg, vdr_FeatureSummarizeSegment($ftr_info_AHR, $sgm_info_AHR, $sgm_idx)));
            push(@alt_scoords_A, sprintf("seq:%s;", vdr_CoordsSinglePositionSegmentCreate($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"}, 
                                                                                         $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"}, $FH_HR)));
            push(@alt_mcoords_A, sprintf("mdl:%s;", vdr_CoordsSinglePositionSegmentCreate($sgm_start_rfpos, $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"}, $FH_HR)));
            push(@alt_ftr_A, $ftr_idx);
          }
        }
        if(! $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"3trunc"}) { 
          if($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"stopgap"}) { 
            push(@alt_code_A,    "indf3gap");
            my $tmp_alt_str = vdr_FeatureSummarizeSegment($ftr_info_AHR, $sgm_info_AHR, $sgm_idx);
            if($tmp_alt_str eq "") { 
              $tmp_alt_str = "VADRNULL";
            }
            else { 
              $tmp_alt_str =~ s/^\,\s+//; # ", segment 1 of 2" -> "segment 1 of 2"
            }
            push(@alt_str_A,     $tmp_alt_str);
            push(@alt_scoords_A, sprintf("seq:%s;", vdr_CoordsSinglePositionSegmentCreate($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstop"}, 
                                                                                         $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"}, $FH_HR)));
            push(@alt_mcoords_A, sprintf("mdl:%s;", vdr_CoordsSinglePositionSegmentCreate($sgm_stop_rfpos, $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"}, $FH_HR)));
            push(@alt_ftr_A,     $ftr_idx);
          }
          elsif((! $do_glsearch) && (($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"stoppp"} - $ftr_pp_thresh) < (-1 * $small_value))) { # only check PP if it's not a gap
            # report indf3lcc or indf3lcn
            # indf5lcc: if this segment is 3'-most segment of a CDS or 3'-most segment of a feature that has same stop position as a CDS
            # indf5lcn: if not indf5lcc
            if(($sgm_info_AHR->[$sgm_idx]{"is_3p"}) && # segment is 3'-most segment in feature
               ((vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) || # feature is a CDS
                (vdr_SegmentStopIdenticalToCds($ftr_info_AHR, $sgm_info_AHR, $sgm_idx, $FH_HR)))) { # segment has stop same as a CDS stop
              push(@alt_code_A, "indf3lcc");
            }
            else { 
              push(@alt_code_A, "indf3lcn");
            }
            push(@alt_str_A,     sprintf("%.2f<%.2f%s%s", $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"stoppp"}, $ftr_pp_thresh, $ftr_pp_msg, vdr_FeatureSummarizeSegment($ftr_info_AHR, $sgm_info_AHR, $sgm_idx)));
            push(@alt_scoords_A, sprintf("seq:%s;", vdr_CoordsSinglePositionSegmentCreate($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstop"}, 
                                                                                         $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"}, $FH_HR)));
            push(@alt_mcoords_A, sprintf("mdl:%s;", vdr_CoordsSinglePositionSegmentCreate($sgm_stop_rfpos, $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"}, $FH_HR)));
            push(@alt_ftr_A,     $ftr_idx);
          }
        }

        # Debugging print block
        #printf("segment $sgm_idx: $sgm_start_rfpos to $sgm_stop_rfpos\n");
        #foreach my $key ("sstart", "sstop", "mstart", "mstop", "strand", "5seqflush", "3seqflush", "5trunc", "3trunc", "startgap", "stopgap", "startpp", "stoppp") { 
        #printf("\tstored $key $sgm_results_HAHR->{$seq_name}[$sgm_idx]{$key}\n");
        #}
      }
    } # end of 'for(my $sgm_idx = 0; $sgm_idx < $nsgm; $sgm_idx++)'

    # Check if we should doctor, we do this if:
    # 1) we have at least one CDS start/stop that should be doctored 
    # 2) we have either not doctored this seq yet, or only doctored it once
    #    if we've doctored it once, this second round will revert to the original
    #    if we kept doctoring, we'd get into an infinite loop,
    #    this helps us in the following case:
    #
    #    seq GCGTAA-TG   (initial)
    #    RF  GCGTAAATG
    #       stop^  ^start
    #
    #    seq GCGTA-ATG   (after first doctoring)
    #    RF  GCGTAAATG
    #       stop^  ^start
    #
    #    seq GCGTAA-TG   (after second doctoring - which reverts it)
    #    RF  GCGTAAATG
    #       stop^  ^start
    # 
    if(($seq_doctor_flag) && ($seq_doctor_ctr <= 1)) { 
      $seq_doctor_ctr++; 
      for(my $dcr_idx = 0; $dcr_idx < scalar(@do_dcr_idx_A); $dcr_idx++) { 
        my $dcr_output_idx = $do_dcr_idx_A[$dcr_idx];
        if($dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"dcr_type"} eq "delete") { 
          $msa->swap_gap_and_closest_residue($i, 
                                             $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"indel_apos"}, 
                                             $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"before"});
        }
        else { # $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"dcr_type"} eq "insert"
          # rewrite RF by swapping insert and start position, remember this can only happen if we have 1 seq in the alignment
          # because we checked that $do_glsearch is true above, and that if $do_glsearch is true, then we only have 1 alignment
          # but we do another sanity check here
          if($nseq != 1) { 
            ofile_FAIL("ERROR in $sub_name, trying to perform doctoring of insert type, but have more than 1 seq in alignment", 1, $FH_HR);
          }
          my $msa_ss_cons = undef;
          if($msa->has_ss_cons) { $msa_ss_cons = $msa->get_ss_cons; }
          my ($new_rf, $new_sscons, $rf_errmsg) = swap_gap_and_adjacent_nongap_in_rf($msa->get_rf, $msa_ss_cons, 
                                                                                     $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"indel_apos"}, 
                                                                                     $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"before"});
          if($rf_errmsg ne "") { 
            ofile_FAIL("ERROR in $sub_name, trying to rewrite RF for doctored alignment (insert type):\n$rf_errmsg\n", 1, $FH_HR);
          }
          $msa->set_rf($new_rf);
          if(defined $new_sscons) { $msa->set_ss_cons($new_sscons); }

          # update the rf2a_A map
          msa_create_rfpos_to_apos_map($msa, \@rf2a_A, $FH_HR);

          # update the insert information
          my $orig_ins_tok = ($dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"before"}) ? 
              sprintf("%d:%d:1", $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"rfpos"},     $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"new_seq_uapos"}) : 
              sprintf("%d:%d:1", $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"rfpos"} - 1, $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"new_seq_uapos"});
          my $new_ins_tok  = ($dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"before"}) ? 
              sprintf("%d:%d:1", $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"rfpos"} - 1, $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"new_seq_uapos"} - 1) : 
              sprintf("%d:%d:1", $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"rfpos"},     $dcr_output_HAHR->{$seq_name}[$dcr_output_idx]{"new_seq_uapos"} + 1);
          $seq_inserts_HHR->{$seq_name}{"ins"} = vdr_ReplaceInsertTokenInInsertString($seq_inserts_HHR->{$seq_name}{"ins"}, $orig_ins_tok, $new_ins_tok, $FH_HR)
        }
      }
      $i--; # makes it so we'll reevaluate this sequence in next iteration of the loop
    }
    else { 
      # usual case: we did not doctor the alignment (or we already doctored it twice)
      $seq_doctor_ctr = 0; # reset this (we don't want to reset this in main loop above)

      # report any indf{5,3}{gap,loc} alerts for this sequence that we stored in loop above
      for(my $alt_idx = 0; $alt_idx < scalar(@alt_code_A); $alt_idx++) { 
        alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, $alt_code_A[$alt_idx], $seq_name, $alt_ftr_A[$alt_idx], 
                                   sprintf("%s%s%s", $alt_scoords_A[$alt_idx], $alt_mcoords_A[$alt_idx], $alt_str_A[$alt_idx]),  
                                   $FH_HR);
      }
      
      # report any deletinf/deletins/deletina alerts
      for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
        if(defined $ftr_deletinf_alt_msg_HA{$ftr_idx}) { 
          my $nsgm_alt = scalar(@{$ftr_deletinf_alt_msg_HA{$ftr_idx}});
          my $nsgm_tot = vdr_FeatureNumSegments($ftr_info_AHR, $ftr_idx);
          if($nsgm_alt == $nsgm_tot) { 
            # all segments are deleted, report either deletins or
            # deletina (per-sequence) alert, we do NOT report any
            # deletinf alerts, one reason is there is no feature
            # annotation for $ftr_idx in this case
            my $alt_scoords = "seq:VADRNULL;"; # feature is deleted, no sequence info available
            my $alt_mcoords = "mdl:" . $ftr_info_AHR->[$ftr_idx]{"coords"} . ";";
            # determine which code to use depending on is_deletable value from ftr_info
            my $alt_code    = $ftr_info_AHR->[$ftr_idx]{"is_deletable"} ? "deletina" : "deletins";
            alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, $alt_code, $seq_name, 
                                        sprintf("%s%s%s feature number %s: %s",
                                                $alt_scoords, $alt_mcoords,
                                                $ftr_info_AHR->[$ftr_idx]{"type"}, 
                                                vdr_FeatureTypeIndex($ftr_info_AHR, $ftr_idx), 
                                                $ftr_info_AHR->[$ftr_idx]{"outname"}), 
                                        $FH_HR);
          }
          else { 
            # at least one but not all segments are deleted, report 1 or more deletinf (per-feature)
            # alerts, we do NOT report a deletins alert because this feature is 
            # annotated, just not all segments are.
            # NOTE: this won't happen if a segment is not annotated because it is truncated
            # away due to a sequence terminus (i.e. should exist before/after start/end of sequence)
            foreach my $alt_msg (@{$ftr_deletinf_alt_msg_HA{$ftr_idx}}) { 
              alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, "deletinf", $seq_name, $ftr_idx, $alt_msg, $FH_HR);
            }
          }
        }
      }
      # detect and report any frameshifts for this sequence
      add_frameshift_alerts_for_one_sequence($msa, $seq_name, $i, \@rf2a_A, \@rfpos_pp_A, \@rf2ilen_A, 
                                             \@max_uapos_before_A, \@{$sgm_info_HAH{$mdl_name}},
                                             \@{$ftr_info_HAH{$mdl_name}}, \%alt_info_HH, 
                                             \%{$sgm_results_HHAH{$mdl_name}}, \%{$ftr_results_HHAH{$mdl_name}}, 
                                             \%alt_ftr_instances_HHH, $mdl_name, $out_root, \%opt_HH, \%ofile_info_HH);

      # fetch features, detect and add cds and mp alerts for this sequence
      # we have to do this here becuase we need @ua2rf_A map of unaligned positions 
      # to RF positions to report model positions for alerts
      fetch_features_and_add_cds_and_mp_alerts_for_one_sequence($sqfile_for_cds_mp_alerts, $sqfile_for_output_fastas, $sqfile_for_pv,
                                                                $do_separate_cds_fa_files, $mdl_name, $mdl_tt, 
                                                                $seq_name, $seq_len_HR, $ftr_info_AHR, $sgm_info_AHR, $alt_info_HHR, 
                                                                $sgm_results_HAHR, $ftr_results_HAHR, $alt_ftr_instances_HHHR, 
                                                                \@ua2rf_A, $ftr_fileroot_AR, $ftr_outroot_AR, 
                                                                $to_remove_AR, $opt_HHR, $ofile_info_HHR);
      
      # add low similarity alerts for this sequence
      # we have to do this here becuase we need @ua2rf_A map of unaligned positions 
      # to RF positions to report model positions for alerts
      add_low_similarity_alerts_for_one_sequence($seq_name, \%seq_len_H, \@ua2rf_A, 
                                                 $ftr_info_AHR, $sgm_info_AHR, $alt_info_HHR, 
                                                 $stg_results_HHHR, $sgm_results_HAHR, $ftr_results_HAHR, 
                                                 $alt_seq_instances_HHR, $alt_ftr_instances_HHHR, $opt_HHR, $ofile_info_HHR);
      
      # update the model coords for ambgnt5s/ambgnt3s seqs, now that we have the alignment
      alert_sequence_instance_update_mdl_coords($alt_seq_instances_HHR, $alt_info_HHR, "ambgnt5s", $seq_name, \@ua2rf_A, $FH_HR);
      alert_sequence_instance_update_mdl_coords($alt_seq_instances_HHR, $alt_info_HHR, "ambgnt3s", $seq_name, \@ua2rf_A, $FH_HR);

    } # end of 'else' entered if ! $doctor_flag
  } # end of 'for(my $i = 0; $i < $nseq; $i++)'

  if($msa_doctor_flag) { 
    $msa->write_msa($stk_file, "pfam", 0);
  }

  undef $msa;
  return;
}

#################################################################
# Subroutine : add_frameshift_alerts_for_one_sequence()
# Incept:      EPN, Mon Mar  2 19:54:25 2020
#
# Purpose:    Given information about a parsed alignment of a single
#             sequence, detect frameshift alerts and report them. Also
#             create frameshift-annotated stockholm files if --keep or
#             --out_fsstk.
#             
#             Detects and adds the following alerts to 
#             @{$alt_ftr_instances_AAHR}:
#
#             If --glsearch NOT used:
#             fsthicft: CDS has a possible terminal frameshift, high confidence
#             fsthicfi: CDS has a possible internal frameshift, high confidence
#             fstlocft: CDS has a possible terminal frameshift, low confidence
#             fstlocfi: CDS has a possible internal frameshift, low confidence
#
#             If --glsearch IS used:
#             fstukcft: CDS has a possible terminal frameshift, unknown confidence
#             fstukcfi: CDS has a possible internal frameshift, unknown confidence
#
# Arguments: 
#  $msa:                    the ESL_MSA alignment object
#  $seq_name:               the sequence we are analyzing
#  $seq_idx:                index of sequence in $msa
#  $rf2a_AR:                REF to array: [1..$rfpos..$rflen] = $apos;  
#                           rf position $rfpos maps to alignment position $apos [1..$alen]  
#                           ($rf2a_A[0] = -1  (dummy value))
#  $rfpos_pp_AR             REF to array: [0..$rfpos..rflen+1]: posterior probability 
#                           character for current sequence at RF position $rfpos
#                           '.' if sequence is a gap at that RF position $rfpos
#                           special values: $rfpos_pp_A[0] = -1, $rfpos_pp_A[$rflen+1] = -1
#                           undef if --glsearch used ($do_glsearch == 1 in caller)
#  $rf2ilen_AR:             REF to array: [1..$rfpos..$rflen] = $apos; rf position $rfpos maps to 
#                           alignment position $apos [1..$alen]  ($rf2a_A[0] = -1 (dummy value))
#  $max_uapos_before_AR:    REF to array; [0..$rfpos..rflen+1]: maximum unaligned position for 
#                           current sequence that aligns at or inserts *before* 
#                           $max_rfpos_before_A[$rfpos], -1 if $max_rfpos_before_A[$rfpos] == -1
#  $sgm_info_AHR:           REF to hash of arrays with information on the model segments, PRE-FILLED
#  $ftr_info_AHR:           REF to hash of arrays with information on the features, PRE-FILLED
#  $alt_info_HHR:           REF to hash of hashes with information on the errors, PRE-FILLED
#  $sgm_results_HAHR:       REF to results HAH, FILLED HERE
#  $ftr_results_HAHR:       REF to feature results HAH, possibly ADDED TO HERE
#  $alt_ftr_instances_HHHR: REF to error instances HAH, ADDED TO HERE
#  $mdl_name:               model name this alignment pertains to
#  $out_root:               string for naming output files
#  $opt_HHR:                REF to 2D hash of option values
#  $ofile_info_HHR:         REF to 2D hash of output file information
#
# Returns:    void
#
# Dies:
#
################################################################# 
sub add_frameshift_alerts_for_one_sequence { 
  my $sub_name = "add_frameshift_alerts_for_one_sequence";
  my $nargs_exp = 17;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($msa, $seq_name, $seq_idx, $rf2a_AR, $rfpos_pp_AR, $rf2ilen_AR, $max_uapos_before_AR, 
      $sgm_info_AHR, $ftr_info_AHR, $alt_info_HHR, $sgm_results_HAHR, $ftr_results_HAHR, 
      $alt_ftr_instances_HHHR, $mdl_name, $out_root, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = \%{$ofile_info_HHR->{"FH"}};

  my $do_output_frameshift_stk = ((opt_Get("--keep", $opt_HHR)) || (opt_Get("--out_fsstk", $opt_HHR))) ? 1 : 0;
  my $fst_min_ntt    = opt_Get("--fstminntt",   $opt_HHR); # maximum allowed nt length of unexpected frame at 3' end without a frameshift alert 
  my $fst_min_nti    = opt_Get("--fstminnti",   $opt_HHR); # maximum allowed nt length of internal region unexpected frame without a frameshift alert 
  my $fst_high_ppthr = opt_Get("--fsthighthr",  $opt_HHR); # minimum average probability for fsthicnf frameshift alert  
  my $fst_low_ppthr  = opt_Get("--fstlowthr",   $opt_HHR); # minimum average probability for fslowcnf frameshift alert 
  my $nmaxins        = opt_Get("--nmaxins",     $opt_HHR); # maximum allowed insertion length in nucleotide alignment
  my $nmaxdel        = opt_Get("--nmaxdel",     $opt_HHR); # maximum allowed deletion length in nucleotide alignment
  my $do_glsearch    = opt_Get("--glsearch",    $opt_HHR); # we won't have PP values if this is 1
  my $small_value = 0.000001; # for checking if PPs are below threshold
  my $nftr = scalar(@{$ftr_info_AHR});
  my $ftr_idx;
  my $alert_scoords = undef; # sequence coords string for an alert
  my $alert_mcoords = undef; # model coords string for an alert

  # get info on position-specific insert and delete maximum exceptions, and frameshift regions, if there are any
  my @nmaxins_exc_AH  = ();
  my @nmaxdel_exc_AH  = ();
  my @fs_exc_AA       = ();
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    %{$nmaxins_exc_AH[$ftr_idx]} = ();
    %{$nmaxdel_exc_AH[$ftr_idx]} = ();
    @{$fs_exc_AA[$ftr_idx]} = ();
    vdr_FeaturePositionSpecificValueBreakdown($ftr_info_AHR, $ftr_idx, "nmaxins_exc", \%{$nmaxins_exc_AH[$ftr_idx]}, $FH_HR);
    vdr_FeaturePositionSpecificValueBreakdown($ftr_info_AHR, $ftr_idx, "nmaxdel_exc", \%{$nmaxdel_exc_AH[$ftr_idx]}, $FH_HR);
    vdr_FeatureCoordsListValueBreakdown($ftr_info_AHR, $ftr_idx, "frameshift_exc", \@{$fs_exc_AA[$ftr_idx]}, $FH_HR);
  }

  # for each CDS: determine frame, and report frameshift alerts
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) { 
      my $frame_stok_str = ""; # string of ';' delimited tokens that describe subsequence stretches that imply the same frame
      my $frame_mtok_str = ""; # string of ';' delimited tokens that describe model stretches that correspond to $frame_stok_str
      my @frame_ct_A = (0, 0, 0, 0); # [0..3], number of RF positions that 'vote' for each candidate frame (frame_ct_A[0] is invalid and will stay as 0)
      my $ftr_strand = undef; # strand for this feature
      my $ftr_sstart = undef; # starting sequence position of this CDS feature
      my $ftr_sstop  = undef; # ending   sequence position of this CDS feature
      my $ftr_mstart = undef; # starting model position of this CDS feature that $ftr_sstart pertains to
      my $ftr_mstop  = undef; # ending   model position of this CDS feature that $ftr_sstop pertains to
      my $ftr_start_rfpos = undef; # start model position of this CDS (regardless of where sequence alignment to the CDS starts)
      my $ftr_stop_rfpos  = undef; # stop  model position of this CDS (regardless of where sequence alignment to the CDS stops)
      my $nsgm = 0; # number of segments for this CDS
      my @gr_frame_str_A = (); # [0..$nsgm-1] GR annotation of frame per-position per CDS segment, only relevant if a frameshift alert occurs for this CDS
      my @sgm_idx_A = (); # array of segment indices that are covered by this seq/CDS
      my $rf_diff = 0;  # number of rf positions seen since first rf position aligned to a nt for current CDS
      my $ua_diff = 0;  # number of nt seen since first nt for current CDS
      my $F_0 = undef;  # frame of initial nongap RF position for current CDS
      my $full_ppstr = undef; # unaligned posterior probability string for this sequence, only defined if nec (if frameshift alert is reported)
      my @cds_alt_str_A = ();
      my $first_sgm_idx = get_5p_most_sgm_idx_with_results($ftr_info_AHR, $sgm_results_HAHR, $ftr_idx, $seq_name);
      my $final_sgm_idx = get_3p_most_sgm_idx_with_results($ftr_info_AHR, $sgm_results_HAHR, $ftr_idx, $seq_name);
      my $tmp_ftr_len = 0;
      if($first_sgm_idx != -1) { 
        for(my $sgm_idx = $first_sgm_idx; $sgm_idx <= $final_sgm_idx; $sgm_idx++) { 
          #check if sgm is valid, it's possible this segment was completely deleted and thus has no results
          # *even if other segments in this feature* do have results (this is related to the github issue #21 bug)
          if((defined $sgm_results_HAHR->{$seq_name}) && 
             (defined $sgm_results_HAHR->{$seq_name}[$sgm_idx]) && 
             (defined $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"})) { 
            push(@sgm_idx_A, $sgm_idx); # store this segment index
            my $is_first_sgm = ($sgm_idx == $first_sgm_idx) ? 1 : 0;
            my $is_final_sgm = ($sgm_idx == $final_sgm_idx) ? 1 : 0;
            my $gr_frame_str = ""; # GR annotation of frame per-position for this CDS segment, only relevant if a frameshift alert occurs for this CDS
            my $sgm_results_HR = $sgm_results_HAHR->{$seq_name}[$sgm_idx]; # for convenience
            my $sgm_start_rfpos = $sgm_info_AHR->[$sgm_idx]{"start"};
            my $sgm_stop_rfpos  = $sgm_info_AHR->[$sgm_idx]{"stop"};
            if(! defined $ftr_start_rfpos) { $ftr_start_rfpos = $sgm_start_rfpos; }
            $ftr_stop_rfpos  = $sgm_stop_rfpos;
            my $sgm_strand   = $sgm_info_AHR->[$sgm_idx]{"strand"};
            my $sstart = $sgm_results_HR->{"sstart"}; # sequence position this segment starts at
            my $sstop  = $sgm_results_HR->{"sstop"};  # sequence position this segment stops at
            $tmp_ftr_len += abs($sstop - $sstart) + 1;
            my $mstart = ($sgm_idx == $first_sgm_idx) ? $sgm_results_HR->{"mstart"} : $sgm_start_rfpos; 
            my $mstop  = ($sgm_idx == $final_sgm_idx) ? $sgm_results_HR->{"mstop"}  : $sgm_stop_rfpos; 
            my $strand = $sgm_results_HR->{"strand"};
            my $cur_delete_len = 0; # current length of deletion
            if(! defined $ftr_sstart) { $ftr_sstart = $sstart; }
            if(! defined $ftr_mstart) { $ftr_mstart = $mstart; }
            $ftr_sstop = $sstop;
            $ftr_mstop = $mstop;
            if(! defined $F_0) { 
              $F_0 = vdr_FrameAdjust(1, abs($mstart - $sgm_start_rfpos), $FH_HR);
              # $F_0 is frame of initial nongap RF position for this CDS 
            } 

            # sanity checks about strand
            if((defined $ftr_strand) && ($ftr_strand ne $strand)) { 
              ofile_FAIL("ERROR, in $sub_name, different segments of same CDS feature have different strands ... can't deal", 1, $FH_HR);
            }
            $ftr_strand = $strand;
            if($strand ne $sgm_strand) { 
              ofile_FAIL("ERROR, in $sub_name, predicted strand for segment inconsistent with strand from segment info", 1, $FH_HR);
            }
            my $F_prv = undef;     # frame of previous RF position 
            my $uapos_prv = undef; # unaligned sequence position that aligns to previous RF position
            my $rfpos_prv = undef; # previous RF position
            if(($strand ne "+") && ($strand ne "-")) { 
              ofile_FAIL("ERROR, in $sub_name, strand is neither + or -", 1, $FH_HR);
            }
            # for each RF position covered by the predicted segment
            # we want to deal with both + and - strands with same code block, 
            # so can't use a simple for loop 
            my $rfpos = $mstart;
            my $uapos = undef;
            while(($strand eq "+" && $rfpos <= $mstop) || 
                  ($strand eq "-" && $rfpos >= $mstop)) { 
              $rf_diff++; # number of RF positions seen since first nt in this CDS
              if($rfpos_pp_AR->[$rfpos] ne ".") { 
                # this rfpos is not aligned to a gap in the sequence
                # determine uapos, the unaligned sequence position that aligns to RF pos $rfpos
                # $max_uapos_before_AR->[$rfpos] actually gives you the maximum unaligned seq position that 
                # aligns at or inserts before $rfpos, but we know it aligns at $rfpos because we just 
                # checked that it's not a gap (rfpos_pp_A[$rfpos] is not a gap)
                $uapos = $max_uapos_before_AR->[$rfpos]; 
                $ua_diff++; # increment number of nucleotides seen since first nt in this CDS
                my $F_cur = vdr_FrameAdjust($F_0, ($rf_diff - $ua_diff), $FH_HR); # frame implied by current nt aligned to current rfpos
                if($strand eq "+") { $gr_frame_str .= $F_cur; }
                else               { $gr_frame_str  = $F_cur . $gr_frame_str; } # prepend for negative string
                $frame_ct_A[$F_cur]++;
                if((! defined $F_prv) || ($F_cur != $F_prv)) { 
                  # frame changed, 
                  # first complete the previous frame 'token' that described the contiguous subsequence that was in the previous frame
                  if(defined $F_prv) { 
                    $frame_stok_str .= $uapos_prv . "," . "D" . (abs($rfpos - $rfpos_prv) - 1) . ",0;"; # 0: not end of segment
                    $frame_mtok_str .= $rfpos_prv . ",0;" # 0: not end of segment
                    # (($rfpos-$rfpos_prv)-1) part is number of deleted reference positions we just covered
                  } 
                  # and begin the next frame 'token' that will describe the contiguous subsequence that is in the previous frame
                  my $nins = (defined $F_prv) ? (abs($uapos - $uapos_prv) - 1) : 0;
                  my $ins_str = "I" . $nins;
                  $frame_stok_str .= $F_cur . "," . $ins_str . "," . $uapos . "..";
                  $frame_mtok_str .= $F_cur . "," . $rfpos . "..";
                }
                $uapos_prv = $uapos;
                $rfpos_prv = $rfpos;
                $F_prv     = $F_cur;
                my $local_rfpos   = ($strand eq "+") ? ($rfpos - $cur_delete_len) : ($rfpos + $cur_delete_len);
                my $local_nmaxdel = defined ($nmaxdel_exc_AH[$ftr_idx]{$local_rfpos}) ? $nmaxdel_exc_AH[$ftr_idx]{$local_rfpos} : $nmaxdel;
                if($cur_delete_len > $local_nmaxdel) { 
                  $alert_scoords = sprintf("seq:%s;", ($strand eq "+") ? 
                                           vdr_CoordsSegmentCreate($uapos-1, $uapos-1, $strand, $FH_HR) : 
                                           vdr_CoordsSegmentCreate($uapos+1, $uapos+1, $strand, $FH_HR));
                  $alert_mcoords = sprintf("mdl:%s;", ($strand eq "+") ? 
                                           vdr_CoordsSegmentCreate($local_rfpos, $rfpos - 1, $strand, $FH_HR) : 
                                           vdr_CoordsSegmentCreate($local_rfpos, $rfpos + 1, $strand, $FH_HR));
                  alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, "deletinn", $seq_name, $ftr_idx, 
                                             sprintf("%s%s%d>%d", 
                                                     $alert_scoords, $alert_mcoords, $cur_delete_len, $local_nmaxdel), $FH_HR);
                }
                $cur_delete_len = 0;
              }
              else { # rf position is a gap, add 'd' GR frame annotation
                if($strand eq "+") { $gr_frame_str .= "d"; }
                else               { $gr_frame_str =  "d" . $gr_frame_str; } # prepend for negative strand
                $cur_delete_len++;
              }
              # add 'i' GR frame annotation for inserts that occur after (or before if neg strand) this rfpos, if any
              if($strand eq "+") { 
                if(($rfpos < $mstop) && ($rf2ilen_AR->[$rfpos] > 0)) { 
                  for(my $ipos = 0; $ipos < $rf2ilen_AR->[$rfpos]; $ipos++) { 
                    $gr_frame_str .= "i"; 
                    $ua_diff++; # increment number of seq positions seen
                  }
                }
              }
              else { # negative strand, look for inserts that occur before this position
                if(($rfpos > $mstop) && ($rf2ilen_AR->[($rfpos-1)] > 0)) { 
                  for(my $ipos = 0; $ipos < $rf2ilen_AR->[($rfpos-1)]; $ipos++) { 
                    $gr_frame_str =  "i" . $gr_frame_str; # prepend for negative strand
                    $ua_diff++; # increment number of seq positions seen
                  }
                }
              }
              # add insertnn alert, if nec
              my $local_nmaxins = defined ($nmaxins_exc_AH[$ftr_idx]{$rfpos}) ? $nmaxins_exc_AH[$ftr_idx]{$rfpos} : $nmaxins;
              if($rf2ilen_AR->[$rfpos] > $local_nmaxins) { 
                $alert_scoords = sprintf("seq:%s;", ($strand eq "+") ? 
                                         vdr_CoordsSegmentCreate($uapos+1, $uapos+1 + $rf2ilen_AR->[$rfpos]-1, $strand, $FH_HR) : 
                                         vdr_CoordsSegmentCreate($uapos-1, $uapos-1 - $rf2ilen_AR->[$rfpos]+1, $strand, $FH_HR));
                $alert_mcoords = sprintf("mdl:%s;", vdr_CoordsSegmentCreate($rfpos, $rfpos, $strand, $FH_HR));
                alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, "insertnn", $seq_name, $ftr_idx, 
                                           sprintf("%s%s%d>%d", 
                                                   $alert_scoords, $alert_mcoords, $rf2ilen_AR->[$rfpos], $local_nmaxins), $FH_HR);
              }

              # increment or decrement rfpos
              if($strand eq "+") { $rfpos++; } 
              else               { $rfpos--; }
            }
            # complete final frame token
            my $nterm_del = (abs($rfpos - $rfpos_prv) - 1);
            $frame_stok_str .= $uapos . "," . "D" . $nterm_del . ",1;"; # the '1' indicates the end of a segment
            $frame_mtok_str .= sprintf("%d,1;", (($strand eq "+") ? ($rfpos - $nterm_del - 1) : ($rfpos + $nterm_del + 1))); # we already incremented/decremented rfpos for next segment

            $nsgm++;
            push(@gr_frame_str_A, $gr_frame_str);
            my $local_rfpos   = ($strand eq "+") ? ($rfpos - $cur_delete_len) : ($rfpos + $cur_delete_len);
            my $local_nmaxdel = defined ($nmaxdel_exc_AH[$ftr_idx]{$local_rfpos}) ? $nmaxdel_exc_AH[$ftr_idx]{$local_rfpos} : $nmaxdel;
            if($cur_delete_len > $local_nmaxdel) { 
              $alert_scoords = sprintf("seq:%s;", ($strand eq "+") ? 
                                       vdr_CoordsSegmentCreate($uapos-1, $uapos-1, $strand, $FH_HR) : 
                                       vdr_CoordsSegmentCreate($uapos+1, $uapos+1, $strand, $FH_HR));
              $alert_mcoords = sprintf("mdl:%s;", ($strand eq "+") ? 
                                       vdr_CoordsSegmentCreate($local_rfpos, $rfpos - 1, $strand, $FH_HR) : 
                                       vdr_CoordsSegmentCreate($local_rfpos, $rfpos + 1, $strand, $FH_HR));
              alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, "deletinn", $seq_name, $ftr_idx, 
                                         sprintf("%s%s%d>%d", 
                                                 $alert_scoords, $alert_mcoords, $cur_delete_len, $local_nmaxdel), $FH_HR);
            }
          } # end of 'if' entered if segment has valid results 
        } # end of 'for(my $sgm_idx = $first_sgm_idx; $sgm_idx <= $final_sgm_idx; $sgm_idx++) {' 
      } # end of 'if($first_sgm_idx != -1)'

      # debug print block
      #printf("frame_ct_A[1]: $frame_ct_A[1]\n");
      #printf("frame_ct_A[2]: $frame_ct_A[2]\n");
      #printf("frame_ct_A[3]: $frame_ct_A[3]\n");
      #printf("frame_stok_str: $frame_stok_str\n");
      #printf("frame_mtok_str: $frame_mtok_str\n");

      # deconstruct $frame_stok_str, looking for potential frameshifts, 
      # we combine any subseqs not in the expected frame together and
      # then check if any (possibly joined) unexpected frame subseqs
      # are long enough to trigger an alert
      my @frame_stok_A = split(";", $frame_stok_str);
      my @frame_mtok_A = split(";", $frame_mtok_str);
      my $nframe_stok = scalar(@frame_stok_A); # number of sequence frame tokens
      my $nframe_mtok = scalar(@frame_mtok_A); # number of model frame tokens
      if($nframe_stok != $nframe_mtok) { 
        ofile_FAIL("ERROR, in $sub_name, different numbers of sequence and model frame tokens, internal coding error: frame_stok_str: $frame_stok_str, frame_mtok_str: $frame_mtok_str", 1, $FH_HR);
      }

      if($nframe_stok >= 1) { 
        # determine and store dominant frame, the frame with maximum count in @frame_ct_A, frame_ct_A[0] will be 0
        # determine and store expected frame, the expected frame of the CDS
        # expected_frame is the predicted frame of the first position unless
        # the CDS is 5' truncated AND the length of the first region is less than $fst_min_nti length,
        # in which case we set expected_frame to dominant frame. This avoids frameshift calls
        # when first region is very short and we're 5' truncated (when we are 5' truncated we are not 
        # as confident about what the expected frame is b/c we don't have a start codon)
        my $first_span_slen = undef;
        if($frame_stok_A[0] =~ /^[123],I\d+,(\d+)\.\.(\d+),D\d+\,[01]$/) { 
          my ($first_sstart, $first_sstop) = ($1, $2); 
          $first_span_slen = abs($first_sstop - $first_sstart) + 1;
        }
        else { 
          ofile_FAIL("ERROR, in $sub_name, unable to parse frame_stok, internal coding error: $frame_stok_A[0]", 1, $FH_HR);
        }
        my $is_5p_trunc = $sgm_results_HAHR->{$seq_name}[$first_sgm_idx]{"5trunc"};
        my $is_3p_trunc = $sgm_results_HAHR->{$seq_name}[$final_sgm_idx]{"3trunc"};
        my $dominant_frame = utl_AArgMax(\@frame_ct_A);
        my $expected_frame = (($is_5p_trunc) && ($first_span_slen < $fst_min_nti)) ? $dominant_frame : $F_0;
        $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"n_codon_start_expected"} = $expected_frame;
        $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"n_codon_start_dominant"} = $dominant_frame;

        if($nframe_stok > 1) { # if there's only one frame_stok, we can't have a frameshift
          my $prv_sstop   = undef; # last sequence position in the previous frame token
          my $prv_mstop   = undef; # last model    position in the previous frame token
          my $prv_frame  = undef; # frame of previous frame token 
          my $prv_exp_span_sstart  = undef; # last sequence position in the previous expected frame token
          my $prv_exp_span_mstart  = undef; # last model    position in the previous expected frame token
          my $prv_exp_span_sstop   = undef; # last sequence position in the previous expected frame token
          my $prv_exp_span_mstop   = undef; # last model    position in the previous expected frame token
          my $shifted_span_sstart  = undef; # first sequence position of a unexpected frame subseq 
          my $shifted_span_sstop   = undef; # final sequence position of a unexpected frame subseq 
          my $shifted_span_slen    = undef; # length in sequence of a unexpected frame subseq
          my $shifted_span_mstart  = undef; # first model    position of a unexpected frame subseq
          my $shifted_span_mstop   = undef; # final model    position of a unexpected frame subseq
          my $shifted_span_mlen    = undef; # length in model of a unexpected frame subseq
          my $shifted_frame        = undef; # frame of shifted frame token 
          my $exp_span_slen        = undef; # length in sequence of an expected frame subseq
          my @cur_indel_str_A      = ();    # array of indel mutations related to a frameshift
          my $prv_tok_sgm_end_flag = 0;     # flag for previous token being special token indicating end of a segment
          my $is_terminal          = 0;     # set to 1 if the frameshifted region includes 3'-most (final) nt of CDS feature, else 0
          my $f;                            # counter over frame tokens

          # now step through each subseq that has a different frame and report frameshift alerts when necessary
          for($f = 0; $f < $nframe_stok; $f++) { 
            #printf("f: $f frame_stok: %s\n", $frame_stok_A[$f]);
            #printf("f: $f frame_mtok: %s\n", $frame_mtok_A[$f]);
            if($frame_stok_A[$f] =~ /^([123]),I(\d+),(\d+)\.\.(\d+),D(\d+)\,([01])$/) { 
              my ($cur_frame,  $cur_ninsert, $cur_sstart, $cur_sstop, $cur_ndelete, $cur_sgmend) = ($1, $2, $3, $4, $5, $6); 
              my ($cur_mframe, $cur_mstart, $cur_mstop, $cur_msgmend);
              # add to growing list of inserts, if nec
              # we do this before we report an alert because insert info 
              # in the current frame token is relevant to the alert we may be about to report
              # we add to delete info *after* we report an alert because delete
              # info in this frame token is relevant to the next alert we may report
              
              if($frame_mtok_A[$f] =~ /^([123])\,(\d+)\.\.(\d+)\,([01])$/) { 
                ($cur_mframe, $cur_mstart, $cur_mstop, $cur_msgmend) = ($1, $2, $3, $4);
                if($cur_frame != $cur_mframe) { 
                  ofile_FAIL("ERROR, in $sub_name, different frame in sequence and model frame tokens, internal coding error:\nstok:$frame_stok_A[$f]\nmtok:$frame_mtok_A[$f]", 1, $FH_HR);
                }
              }              
              else { 
                ofile_FAIL("ERROR, in $sub_name, unable to parse frame_mtok, internal coding error: $frame_mtok_A[$f]", 1, $FH_HR);
              }

              # deal with inserted positions between previous frame token and this one
              my $cur_insert_str = undef;
              if($f > 0) { 
                # add any inserted positions between previous frame token and this one to insert_str
                if($ftr_strand eq "+") { 
                  if((($prv_sstop + 1) < ($cur_sstart)) && (! $prv_tok_sgm_end_flag)) { # at least one inserted nt and previous token was not a segment end
                    if(($prv_sstop + 1) == ($cur_sstart - 1)) { # exactly one inserted nt
                      $cur_insert_str = sprintf("insert,S:%d(%d),M:%d;", ($prv_sstop + 1), 1, $prv_mstop); 
                    }
                    else { # more than one inserted nt, specify the range
                      $cur_insert_str = sprintf("insert,S:%d..%d(%d),M:%d;", $prv_sstop+1, $cur_sstart-1, (abs(($prv_sstop+1) - ($cur_sstart-1))+1), $prv_mstop); 
                    }
                  }
                }
                else { # negative strand
                  if((($prv_sstop - 1) > ($cur_sstart)) && (! $prv_tok_sgm_end_flag)) { # at least one inserted nt and previous token was not a segment end
                    if(($prv_sstop - 1) == ($cur_sstart + 1)) { # exactly one inserted nt
                      $cur_insert_str = sprintf("insert,S:%d(%d),M:%d;", ($prv_sstop - 1), 1, $prv_mstop); 
                    }
                    else { # more than one inserted nt, specify the range
                      $cur_insert_str = sprintf("insert,S:%d..%d(%d),M:%d;", $prv_sstop-1, $cur_sstart+1, (abs(($prv_sstop-1) - ($cur_sstart+1))+1), $prv_mstop); 
                    }
                  }
                }
              }
              if(defined $cur_insert_str) { 
                push(@cur_indel_str_A, $cur_insert_str);
                # printf("\npushing insert str $cur_insert_str to cur_indel_str_A (new size: %d)\n", scalar(@cur_indel_str_A));
              }
                
              # Determine if we may have a frameshift alert
              # Two possible cases:
              # Case 1: this subseq is in expected frame, but previous was not (that is, it's not the first frame_stok ($f != 0))
              # Case 2: this subseq is not in expected frame and it's the final one ($f == ($nframe_stok - 1))
              # We purposefully do not detect frameshifts when the current subseq is in the expected frame and 
              # the previous frame was ! expected because we need to be able to detect the case when there is
              # two subseqs in a row that are both unexpected frame (e.g. 1111133333222221111111)
              # The only way to detect this and correctly report its span is when you are back in the expected
              # frame at the end of the CDS. 
              if((($cur_frame == $expected_frame) && ($f > 0) && ($prv_frame != $expected_frame)) ||  # Case 1
                 (($cur_frame != $expected_frame) && ($f == ($nframe_stok-1)))) {  # Case 2
                $is_terminal   = 0; # set to '1' below if frameshift region includes 3'-most nt of CDS feature
                $shifted_frame = undef; # will save shifted frame for alert output
                
                # determine $shifted_span_sstart: the first position of the unexpected frame subseq
                if(defined $prv_exp_span_sstop) { 
                  # we've seen at least one expected frame segment,
                  # start of the unexpected stretch is 1 nt 3' of that
                  $shifted_span_sstart = ($ftr_strand eq "+") ? $prv_exp_span_sstop + 1 : $prv_exp_span_sstop - 1;
                  $shifted_span_mstart = ($ftr_strand eq "+") ? $prv_exp_span_mstop + 1 : $prv_exp_span_mstop - 1;
                }
                else { 
                  # we haven't seen a expected frame segment yet, 
                  # span start is first nt of CDS ($ftr_sstart)
                  $shifted_span_sstart = $ftr_sstart; 
                  $shifted_span_mstart = $ftr_mstart; 
                }
                # determine $shifted_span_sstop: the final position of the unexpected frame subseq
                if(($cur_frame != $expected_frame) && ($f == ($nframe_stok-1))) { 
                  # (case 2) this subseq is not in expected frame and it's the final one ($f == ($nframe_stok - 1))
                  # so final nt of the unexpected stretch is the final nt of the CDS ($ftr_sstop) 
                  $shifted_span_sstop = $ftr_sstop;
                  $shifted_span_mstop = $ftr_mstop;
                  $is_terminal = 1; 
                  $shifted_frame = $cur_frame;
                }
                else { 
                  # (case 1) previous frame token was a unexpected frame, so final nt of that unexpected stretch
                  # is 1 nt 5' of start of current frame token
                  $shifted_span_sstop = ($ftr_strand eq "+") ? $cur_sstart - 1 : $cur_sstart + 1;
                  $shifted_span_mstop = ($ftr_strand eq "+") ? $cur_mstart - 1 : $cur_mstart + 1;
                  $shifted_frame = $prv_frame;
                }
                $shifted_span_slen = abs($shifted_span_sstop - $shifted_span_sstart) + 1;
                $shifted_span_mlen = abs($shifted_span_mstop - $shifted_span_mstart) + 1;
                $exp_span_slen     = abs($prv_exp_span_sstop - $prv_exp_span_sstart) + 1;
                
                # check if this is an exempted region
                my $exempted_region = 0;
                foreach my $exc_coords (@{$fs_exc_AA[$ftr_idx]}) { 
                  if(vdr_CoordsCheckIfSpans($exc_coords, vdr_CoordsSegmentCreate($shifted_span_mstart, $shifted_span_mstop, $ftr_strand, $FH_HR), $FH_HR)) { 
                    $exempted_region = 1;
                  }
                }
                if(! $exempted_region) { 
                  if(((  $is_terminal) && ($shifted_span_slen >= $fst_min_ntt)) || 
                     ((! $is_terminal) && ($shifted_span_slen >= $fst_min_nti))) { 
                    # above our length threshold, if $do_glsearch, we always report this, if not it depends on the avg PP value

                    # determine causative, restorative and intermediate strings
                    my ($causative_indel_str, $intermediate_indel_str, $restorative_indel_str) = (undef, undef, undef);
                    my $nindel_str = scalar(@cur_indel_str_A);
                    if($nindel_str > 0) { 
                      $causative_indel_str = $cur_indel_str_A[0]; # first indel is the causative one
                      my $nindel_final_intermediate = undef;
                      if($cur_frame == $expected_frame) { # case 1: we just restored frame, so we know we have a restorative mutation
                        if($nindel_str > 1) { 
                          $restorative_indel_str = $cur_indel_str_A[($nindel_str-1)];
                        }
                        $nindel_final_intermediate = $nindel_str-2; # if this is < 1 we have zero intermediate mutations
                      }
                      else { # case 2: we got to end of CDS and did not restore frame
                        $nindel_final_intermediate = $nindel_str-1; # if this is < 1 we have zero intermediate mutations
                      }
                      for(my $z = 1; $z <= $nindel_final_intermediate; $z++) { 
                        $intermediate_indel_str .= $cur_indel_str_A[$z];
                      }
                    }

                    # determine frame summary string 
#                    my ($frame_sum_str, $length_sum_str) = 
                    my ($frame_sum_str, $length_sum_str, $tmp_length_sum) = 
                        determine_frame_and_length_summary_strings(\@frame_stok_A, 
                                                                   (($cur_frame != $expected_frame) && ($f == ($nframe_stok-1))) ? $f : ($f-1),
                                                                   $expected_frame, $is_5p_trunc, $is_3p_trunc, $FH_HR);

                    if($do_glsearch) { # we don't have PP values, so all frameshifts are treated equally
                      my $alt_code = ($is_terminal) ? "fstukcft" : "fstukcfi";
                      my $alt_scoords_tok =      vdr_CoordsSegmentCreate(    $shifted_span_sstart,      $shifted_span_sstop,  $ftr_strand, $FH_HR);
                      my $alt_mcoords = "mdl:" . vdr_CoordsSegmentCreate(abs($shifted_span_mstart), abs($shifted_span_mstop), $ftr_strand, $FH_HR) . ";";
                      my $alt_scoords = "seq:" . $alt_scoords_tok . ";";
                      my $alt_str  = sprintf("%s%s", $alt_scoords, $alt_mcoords);
                      $alt_str .= sprintf("cause:%s", $causative_indel_str);
                      if(defined $restorative_indel_str)  { $alt_str .= sprintf(" restore:%s", $restorative_indel_str); }
                      $alt_str .= sprintf(" frame:%s;", $frame_sum_str);
                      $alt_str .= sprintf(" length:%s;", $length_sum_str);
                      if(defined $intermediate_indel_str) { $alt_str .= sprintf(" intermediate:%s",   $intermediate_indel_str); }
                      alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, $alt_code, $seq_name, $ftr_idx, $alt_str, $FH_HR);
                      push(@cds_alt_str_A, $alt_str);
                    }
                    else { # $do_glsearch is 0 so we have PP values and we examine them to determine type of frameshift
                      # this *may* be a fstlocf{t,i} or fsthicf{t,i} alert, depending on the average PP of the shifted region
                      # determine average posterior probability of unexpected frame subseq
                      if(! defined $full_ppstr) { 
                        $full_ppstr = $msa->get_ppstring_aligned($seq_idx); 
                        $full_ppstr =~ s/[^0123456789\*]//g; # remove gaps, so we have 1 character in $full_ppstr per nt in the sequence
                      }
                      my $shifted_span_ppstr = ($ftr_strand eq "+") ? 
                          substr($full_ppstr, $shifted_span_sstart - 1, ($shifted_span_slen)) : 
                          substr($full_ppstr, $shifted_span_sstop  - 1, ($shifted_span_slen));
                      my $exp_span_ppstr = ($ftr_strand eq "+") ? 
                          substr($full_ppstr, $prv_exp_span_sstart - 1, ($exp_span_slen)) : 
                          substr($full_ppstr, $prv_exp_span_sstop  - 1, ($exp_span_slen));
                      my $shifted_span_avgpp;
                      my $exp_span_avgpp;
                      ($shifted_span_avgpp, undef) = Bio::Easel::MSA->get_ppstr_avg($shifted_span_ppstr);
                      ($exp_span_avgpp,     undef) = Bio::Easel::MSA->get_ppstr_avg($exp_span_ppstr);
                      #printf("shifted_span: $shifted_span_sstart..$shifted_span_sstop $shifted_span_avgpp, $shifted_span_ppstr\n");
                      #printf("exp_span: $prv_exp_span_sstart..$prv_exp_span_sstop $exp_span_avgpp, $exp_span_ppstr\n");
                      my $min_span_avgpp = utl_Min($shifted_span_avgpp, $exp_span_avgpp);
                      if($min_span_avgpp > ($fst_low_ppthr - $small_value)) { # we have a fstlocnf or fsthicnf alert
                        my $hi_alt_code = ($is_terminal) ? "fsthicft" : "fsthicfi";
                        my $lo_alt_code = ($is_terminal) ? "fstlocft" : "fstlocfi";
                        my $alt_scoords_tok =      vdr_CoordsSegmentCreate(    $shifted_span_sstart,      $shifted_span_sstop,  $ftr_strand, $FH_HR);
                        my $alt_mcoords = "mdl:" . vdr_CoordsSegmentCreate(abs($shifted_span_mstart), abs($shifted_span_mstop), $ftr_strand, $FH_HR) . ";";
                        my $alt_scoords = "seq:" . $alt_scoords_tok . ";";
                        my $alt_str  = sprintf("%s%s", $alt_scoords, $alt_mcoords);
                        $alt_str .= sprintf("cause:%s", $causative_indel_str);
                        if(defined $restorative_indel_str)  { $alt_str .= sprintf(" restore:%s", $restorative_indel_str); }
                        $alt_str .= sprintf(" frame:%s;", $frame_sum_str);
                        $alt_str .= sprintf(" length:%s;", $length_sum_str);
                        $alt_str .= sprintf(" shifted_avgpp:%.3f;", $shifted_span_avgpp);
                        $alt_str .= sprintf(" exp_avgpp:%.3f;", $exp_span_avgpp);
                        if(defined $intermediate_indel_str) { $alt_str .= sprintf(" intermediate:%s",   $intermediate_indel_str); }
                        my $is_hicnf = ($min_span_avgpp > ($fst_high_ppthr - $small_value)) ? 1 : 0;
                        alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, 
                                                   ($is_hicnf) ? $hi_alt_code : $lo_alt_code,
                                                   $seq_name, $ftr_idx, $alt_str, $FH_HR);
                        push(@cds_alt_str_A, $alt_str);
                      }
                    }
                  }
                }
                @cur_indel_str_A = (); # reset this so we can start storing indels relevant to next fs (if we have one)
              } # end of 2 case if entered if we have a frameshift alert

              # deal with deleted positions, using an array, when we report a FS we determine which is 
              # causative, intermediate and restorative.
              my $cur_delete_str = undef;
              if($f != ($nframe_stok-1)) { 
                if($cur_ndelete > 0) { 
                  if($ftr_strand eq "+") { 
                    if($cur_ndelete == 1) { 
                      $cur_delete_str = sprintf("delete,S:%d,M:%d(%d);", $cur_sstop, ($cur_mstop+1), $cur_ndelete);
                    }
                    else { 
                      $cur_delete_str = sprintf("delete,S:%d,M:%d..%d(%d);", $cur_sstop, ($cur_mstop+1), ($cur_mstop+$cur_ndelete), $cur_ndelete);
                    }
                  }
                  else { # negative strand
                    if($cur_ndelete == 1) { 
                      $cur_delete_str = sprintf("delete,S:%d,M:%d(%d);", $cur_sstop, ($cur_mstop-1), $cur_ndelete);
                    }
                    else { 
                      $cur_delete_str = sprintf("delete,S:%d,M:%d..%d(%d);", $cur_sstop, ($cur_mstop-1), ($cur_mstop-$cur_ndelete), $cur_ndelete);
                    }
                  }
                }
              }
              if(defined $cur_delete_str) { 
                push(@cur_indel_str_A, $cur_delete_str);
                #printf("\npushing delete str $cur_delete_str to cur_indel_str_A (new size: %d)\n", scalar(@cur_indel_str_A));
              }

              # keep track of previous values we may need in next loop iteration
              if($cur_frame == $expected_frame) { 
                $prv_exp_span_sstart = $cur_sstart; 
                $prv_exp_span_sstop  = $cur_sstop; 
                $prv_exp_span_mstart = $cur_mstart; 
                $prv_exp_span_mstop  = $cur_mstop; 
              }
              $prv_sstop = $cur_sstop;
              $prv_mstop = $cur_mstop;
              $prv_frame = $cur_frame;
              $prv_tok_sgm_end_flag = $cur_sgmend;
            } # end if statement that parses $frame_stok_A[$f]
            else { 
              ofile_FAIL("ERROR, in $sub_name, unable to parse frame_stok, internal coding error: $frame_stok_A[$f]", 1, $FH_HR);
            }
          } # end of 'for(my $f = 0; $f < $nframe_stok; $f++) {'
        } # end of 'if($nframe_stok > 1)'
      } # end of 'if($nframe_stok >= 1)'

      if(scalar(@cds_alt_str_A) > 0) { 
        # create and output a stockholm file for each segment of this seq/CDS 
        # remove all sequences other than the one we want
        my $msa_nseq = $msa->nseq;
        for(my $s = 0; $s < $nsgm; $s++) { 
          my $sgm_idx      = $sgm_idx_A[$s];
          my $gr_frame_str = $gr_frame_str_A[$s];
          my $sgm_results_HR = $sgm_results_HAHR->{$seq_name}[$sgm_idx]; # for convenience
          my $sgm_start_rfpos = $sgm_info_AHR->[$sgm_idx]{"start"};
          my $sgm_stop_rfpos  = $sgm_info_AHR->[$sgm_idx]{"stop"};
          my $mstart = ($sgm_idx == $first_sgm_idx) ? $sgm_results_HR->{"mstart"} : $sgm_start_rfpos; 
          my $mstop  = ($sgm_idx == $final_sgm_idx) ? $sgm_results_HR->{"mstop"}  : $sgm_stop_rfpos; 
          my $sgm_strand = $sgm_info_AHR->[$sgm_idx]{"strand"};

          my @cds_sgm_seq_A = ();
          for(my $i2 = 0; $i2 < $msa_nseq; $i2++) { $cds_sgm_seq_A[$i2] = 0; }
          $cds_sgm_seq_A[$seq_idx] = 1; # keep this one seq
          my $cds_sgm_msa = $msa->sequence_subset(\@cds_sgm_seq_A);
          my $alen = $cds_sgm_msa->alen;
          $cds_sgm_msa->addGC_rf_column_numbers(); # number RF columns

          # remove columns outside the CDS segment
          my @cds_sgm_col_A = (); # [0..$alen-1], 0 to remove this column, 1 to keep (if within CDS) 
          my $cds_sgm_apos_start = ($sgm_strand eq "+" ? $rf2a_AR->[$mstart] : $rf2a_AR->[$mstop])  - 1; # -1 puts it into 0..alen-1 coords
          my $cds_sgm_apos_stop  = ($sgm_strand eq "+" ? $rf2a_AR->[$mstop]  : $rf2a_AR->[$mstart]) - 1; # -1 puts it into 0..alen-1 coords
          for(my $a = 0;                      $a <  $cds_sgm_apos_start; $a++) { $cds_sgm_col_A[$a] = 0; } # before CDS
          for(my $a = $cds_sgm_apos_start;    $a <= $cds_sgm_apos_stop;  $a++) { $cds_sgm_col_A[$a] = 1; } # CDS
          for(my $a = $cds_sgm_apos_stop + 1; $a <  $alen;               $a++) { $cds_sgm_col_A[$a] = 0; } # after CDS
          $cds_sgm_msa->column_subset(\@cds_sgm_col_A);

          # remove all gap columns
          $cds_sgm_msa->remove_all_gap_columns(1); # 1: don't delete any nongap RF columns
          
          # add GR annotation
          # printf("in $sub_name, seq_name: $seq_name\n");
          # printf("gr_frame_str len: " . length($gr_frame_str) . "\n");
          # print("$gr_frame_str\n");
          # printf("alen: %d\n", $cds_sgm_msa->alen());
          # printf("adding GR (CS) for $seq_name\n");
          $cds_sgm_msa->addGR("CS", 0, $gr_frame_str);
          
          # output alignment, if nec
          if($do_output_frameshift_stk) { 
            my $cds_and_sgm_idx = vdr_FeatureTypeAndTypeIndexString($ftr_info_AHR, $ftr_idx, ".") . "." . ($sgm_idx - $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"} + 1);
            my $stk_file_name   = $out_root . "." . $mdl_name . "." . $cds_and_sgm_idx . ".frameshift.stk";
            # add comment to stockholm file
            my $comment = "Alignment of CDS " . vdr_FeatureTypeIndex($ftr_info_AHR, $ftr_idx);
            $comment .= " segment " . vdr_FeatureRelativeSegmentIndex($ftr_info_AHR, $ftr_idx, $sgm_idx);
            $comment .= " of " . vdr_FeatureNumSegments($ftr_info_AHR, $ftr_idx);
            $comment .= " for sequence " . $cds_sgm_msa->get_sqname(0); 
            $comment .= " to model $mdl_name with at least one frameshift alert (possibly in a different segment for multi-segment CDS).";
            $cds_sgm_msa->addGF("CC", $comment);
            $comment  = "GR CS annotation indicates the codon_start value each nongap RF position implies.";
            $cds_sgm_msa->addGF("CC", $comment);
            $comment  = "Changes from the expected codon_start value indicate possibly frameshifted regions.";
            $cds_sgm_msa->addGF("CC", $comment);
            for(my $c = 0; $c < scalar(@cds_alt_str_A); $c++) { 
              my ($instance_scoords, $instance_mcoords, $instance_detail) = alert_instance_parse($cds_alt_str_A[$c]);
              $cds_sgm_msa->addGS("FS." . ($c+1), "seq_coords:$instance_scoords; mdl_coords: $instance_mcoords; $instance_detail", 0); # 0: seq idx
            }
            # change RF to DNA
            my $cds_sgm_msa_rf = $cds_sgm_msa->get_rf();
            seq_SqstringDnaize(\$cds_sgm_msa_rf);
            $cds_sgm_msa->set_rf($cds_sgm_msa_rf);
            # output to potentially already existent alignment file
            $cds_sgm_msa->write_msa($stk_file_name, "stockholm", 1); # 1: append to file if it exists
            my $stk_file_key = $mdl_name . "." . $cds_and_sgm_idx . ".frameshift.stk";
            if(! defined $ofile_info_HHR->{"fullpath"}{$stk_file_key}) { 
              ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $stk_file_key, $stk_file_name, 1, 1, "Stockholm file for >= 1 possible frameshifts for $cds_and_sgm_idx for model $mdl_name");
            }
            undef $cds_sgm_msa;
          }
        }
      }
    } # end of 'if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx))'
  } # end of 'for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++)'


  undef $msa;

  return;
}

#################################################################
# Subroutine:  determine_frame_and_length_summary_strings()
# Incept:      EPN, Mon Nov 15 12:54:06 2021
#
# Purpose:     Determine a string (<$frame_str>) that summarizes a
#              frameshift within the context of frames of all subseqs
#              in a CDS, and another string (<$length_str>) that
#              summarizes the lengths of each region in a different
#              frame
#
#              Examples:                                                frame_str  length_str
#              non-truncated CDS w/internal frameshift in frame 2:      "1(2)1"    "21,(13),102"
#              5' truncated CDS that starts in frame 2 and shifts to 3: "<2(3)"    "<102,(43)"
#              same as above but also 3' truncated:                     "<2(3)>"   "<102,(43)>"
#
# Arguments: 
#  $frame_stok_AR:      reference to array of tokens of sequence subregions in different frames
#  $idx:                index of frameshifted subregion in $frame_stok_AR
#  $expected_frame:     expected frame, shifted regions are not in this frame
#  $is_5p_trunc:        '1' if the relevant CDS is truncated on 5' end, 0 if not
#  $is_3p_trunc:        '1' if the relevant CDS is truncated on 3' end, 0 if not
#  $FH_HR:              ref to file handle hash
# 
# Returns:     Two values:
#              $frame_str:  string describing frameshifted region in context of full CDS
#              $length_str: string describing length of frameshifted region in context of full CDS
#
# Dies: if unable to parse a token of @{$frame_stok_AR}
#
################################################################# 
sub determine_frame_and_length_summary_strings {
  my $sub_name = "determine_frame_and_length_summary_strings";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($frame_stok_AR, $idx, $expected_frame, $is_5p_trunc, $is_3p_trunc, $FH_HR) = @_;

  my $nframe_stok = scalar(@{$frame_stok_AR});

  my $prv_frame             = undef; # frame of previous token looked at
  my $prv_sum_length        = 0;     # summed length of previous token(s) all in same frame
  my $cur_length            = 0;     # length of current token
  my $parentheses_open_flag = 0;     # set to '1' when we should add left parentheses before next frame/length substr token
  my $parentheses_end_flag  = 0;     # set to '1' when we should add right parentheses after next frame/length substr token
  my $save_ninsert          = 0;     # number of inserts seen since last unexpected frame region
  my $ntok_added            = 0;     # number of substring tokens added to $ret_frame_str and $ret_length_str
  my $length_tok            = undef; # substring token to add to $ret_length_str
  my $frame_tok             = undef; # substring token to add to $ret_frame_str

  # variables parsed from $frame_stok_AR->[$f]
  my $cur_frame   = undef;
  my $cur_ninsert = undef;
  my $cur_sstart  = undef;
  my $cur_sstop   = undef; 
  my $cur_ndelete = undef; 
  my $cur_sgmend  = undef;

  # start with '>' if 5' truncated
  my $ret_frame_str  = ($is_3p_trunc) ? ">" : "";
  my $ret_length_str = ($is_3p_trunc) ? ">" : "";

  my $tmp_ret_len_sum = 0;
  # we go backwards so we can more easily handle cases where there are two shifted non-expected regions adjacent to each other
  for(my $f = ($nframe_stok-1); $f >= 0; $f--) { 
    if($frame_stok_AR->[$f] =~ /^([123]),I(\d+),(\d+)\.\.(\d+),D(\d+)\,([01])$/) { 
      ($cur_frame, $cur_ninsert, $cur_sstart, $cur_sstop, $cur_ndelete, $cur_sgmend) = ($1, $2, $3, $4, $5, $6); 

      if((defined $prv_frame) && ($cur_frame != $prv_frame)) { 
        # add to length and frame str

        # potentially add right parentheses
        $length_tok = $prv_sum_length;
        $frame_tok  = $prv_frame;
        if($parentheses_end_flag) { 
          $length_tok .= ")";
          $frame_tok  .= ")";
          $parentheses_end_flag  = 0;
          $parentheses_open_flag = 1;
        }
        # potentially prepend left parentheses 
        if(($parentheses_open_flag) && ($cur_frame == $expected_frame)) { 
          $length_tok = "(" . $length_tok; 
          $frame_tok  = "(" . $frame_tok;
          $parentheses_open_flag = 0;
        }
        if($ntok_added > 0) { 
          $length_tok .= ":";
        }
        $ret_length_str = $length_tok . $ret_length_str;
        $ret_frame_str  = $frame_tok  . $ret_frame_str;
        $ntok_added++;
        $tmp_ret_len_sum += $prv_sum_length;
        $prv_sum_length = 0;
      }
      if($f == $idx) { 
        $parentheses_end_flag = 1;
      }

      # determine $cur_length
      $cur_length = abs($cur_sstop - $cur_sstart) + 1;
      if($cur_frame != $expected_frame) { 
        $cur_length += $cur_ninsert;
        $cur_length += $save_ninsert;
        $save_ninsert = 0;
      }
      else { # $cur_frame == $expected_frame
        $save_ninsert += $cur_ninsert;
      }

      $prv_frame       = $cur_frame;
      $prv_sum_length += $cur_length;
    }
    else { 
      ofile_FAIL("ERROR, in $sub_name, unable to parse frame_stok, internal coding error: " . $frame_stok_AR->[$f], 1, $FH_HR);
    }
  }

  # add final frame and length

  # potentially add right parenthesis
  $length_tok = $prv_sum_length;
  $frame_tok  = $prv_frame;
  if($parentheses_end_flag) { 
    $length_tok .= ")";
    $frame_tok  .= ")";
    $parentheses_end_flag  = 0;
    $parentheses_open_flag = 1;
  }
  # potentially prepend left parentheses 
  if(($parentheses_open_flag) && ($cur_frame == $expected_frame)) { 
    $length_tok = "(" . $length_tok; 
    $frame_tok  = "(" . $frame_tok;
    $parentheses_open_flag = 0;
  }
  if($ntok_added > 0) { 
    $length_tok .= ":";
  }
  $ret_length_str = $length_tok . $ret_length_str;
  $ret_frame_str  = $frame_tok  . $ret_frame_str;

  $tmp_ret_len_sum += $prv_sum_length;

  # prepend '<' if 5' truncated
  if($is_5p_trunc) { 
    $ret_frame_str  = "<" . $ret_frame_str;
    $ret_length_str = "<" . $ret_length_str;
  }

  return ($ret_frame_str, $ret_length_str, $tmp_ret_len_sum);
}

#################################################################
# Subroutine:  cmalign_store_overflow()
# Incept:      EPN, Wed Feb 13 16:04:53 2019
#
# Purpose:     Store information on a sequence that has caused
#              a DP matrix memory overflow. 
#              
# Arguments: 

#  $seq_file:           the sequence file with the single sequence that failed in it
#  $mxsize:             matrix size to add to @{$overflow_mxsize_AR}
#  $overflow_seq_AR:    ref to array of sequences that failed due to matrix overflows, to add to
#  $overflow_mxsize_AR: ref to array of required matrix sizes for each sequence that failed due to matrix overflows, to add to
#  $FH_HR:              ref to file handle hash
# 
# Returns:     void
#
# Dies: if there's some problem opening the sequence file
#
################################################################# 
sub cmalign_store_overflow { 
  my $sub_name = "cmalign_store_overflow";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seq_file, $mxsize, $overflow_seq_AR, $overflow_mxsize_AR, $FH_HR) = @_;

  my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $seq_file }); # the sequence file object

  my ($overflow_seq, undef) = $sqfile->fetch_next_seq_name_length();
  push(@{$overflow_seq_AR},    $overflow_seq);
  push(@{$overflow_mxsize_AR}, $mxsize);

  $sqfile = undef;

  return;
}

#################################################################
# Subroutine: fetch_features_and_add_cds_and_mp_alerts_for_one_sequence()
# Incept:     EPN, Thu May 20 09:33:04 2021
#             EPN, Fri Feb 22 14:25:49 2019 (prv version for all sequences)
#
# Purpose:   For each sequence, fetch each feature sequence, and 
#            detect mutstart, cdsstopn, mutendcd, mutendns, mutendex, and unexleng alerts 
#            where appropriate. For cdsstopn alerts, correct the predictions
#            and fetch the corrected feature.
#
# Arguments:
#  $sqfile_for_cds_mp_alerts:  REF to Bio::Easel::SqFile object, open sequence file with sequences
#                              to fetch CDS and mat_peptides from to analyze for possible alerts 
#  $sqfile_for_output_fastas:  REF to Bio::Easel::SqFile object, open sequence file with sequences
#                              that we'll fetch feature sequences from to output to per-feature fasta files
#  $sqfile_for_pv:             REF to Bio::Easel::SqFile object, open sequence file with sequences
#                              that we'll fetch feature sequences from for protein validation
#  $do_separate_cds_fa_files:  '1' to output two sets of cds files, one with fetched features from $sqfile_for_output_fastas
#                              and one for the protein validation stage fetched from $sqfile_for_cds_mp_alerts
#  $mdl_tt:                    the translation table ('1' for standard)
#  $seq_name:                  name of sequence we are processing
#  $seq_len_HR:                REF to hash of sequence lengths, PRE-FILLED
#  $ftr_info_AHR:              REF to hash of arrays with information on the features, PRE-FILLED
#  $sgm_info_AHR:              REF to hash of arrays with information on the model segments, PRE-FILLED
#  $alt_info_HHR:              REF to the alert info hash of arrays, PRE-FILLED
#  $sgm_results_HAHR:          REF to model segment results HAH, pre-filled
#  $ftr_results_HAHR:          REF to feature results HAH, added to here
#  $alt_ftr_instances_HHHR:    REF to array of 2D hashes with per-feature alerts, PRE-FILLED
#  $ua2rf_AR:                  REF to array that maps unaligned positions to reference positions
#                              [1..$uapos..$ualen]: reference position that unaligned position $uapos aligns to 
#                              if $ua2rf_A[$uapos] <  0, $uapos inserts *after* ref posn (-1 * $ua2rf_A[$uapos])
#                              if $ua2rf_A[$uapos] == 0, $uapos inserts *before* ref posn 1
#                              $ua2rf_A[0] is invalid (set to 0)
#  $ftr_fileroot_AR:           REF to array of per-feature file root values, pre-calc'ed and passed in so we don't need to do it per-seq
#  $ftr_outroot_AR:            REF to array of per-feature output root values, pre-calc'ed and passed in so we don't need to do it per-seq
#  $to_remove_AR:              REF to array of files to remove before exiting, possibly added to here if $do_separate_cds_fa_files
#  $opt_HHR:                   REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:            REF to the 2D hash of output file information
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub fetch_features_and_add_cds_and_mp_alerts_for_one_sequence { 
  my $sub_name = "fetch_features_and_add_cds_and_mp_alerts_for_one_sequence";
  my $nargs_exp = 20;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqfile_for_cds_mp_alerts, $sqfile_for_output_fastas, $sqfile_for_pv,
      $do_separate_cds_fa_files, $mdl_name, $mdl_tt, 
      $seq_name, $seq_len_HR, $ftr_info_AHR, $sgm_info_AHR, $alt_info_HHR, 
      $sgm_results_HAHR, $ftr_results_HAHR, $alt_ftr_instances_HHHR, 
      $ua2rf_AR, $ftr_fileroot_AR, $ftr_outroot_AR, $to_remove_AR, 
      $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"}; # for convenience
  my $nftr = scalar(@{$ftr_info_AHR});
  my $nsgm = scalar(@{$sgm_info_AHR});

  my $atg_only    = opt_Get("--atgonly", $opt_HHR);
  my $do_keep     = opt_Get("--keep", $opt_HHR);
  my $do_allfasta = ($do_keep || opt_Get("--out_allfasta", $opt_HHR)) ? 1 : 0;

  my $ftr_idx;

  my $sqfile_for_cds_mp_alerts_path = $sqfile_for_cds_mp_alerts->path;
  my $sqfile_for_output_fastas_path = $sqfile_for_output_fastas->path;
  my $sqfile_for_pv_path            = $sqfile_for_pv->path;

  my $seq_len  = $seq_len_HR->{$seq_name};

  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    my $ftr_is_cds_or_mp = vdr_FeatureTypeIsCdsOrMatPeptide($ftr_info_AHR, $ftr_idx);
    my $ftr_is_cds       = vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx);
    my $ftr_matches_cds  = vdr_FeatureTypeIsCdsOrIdStartAndStop($ftr_info_AHR, $ftr_idx);
    my $ftr_is_mp        = vdr_FeatureTypeIsMatPeptide($ftr_info_AHR, $ftr_idx);
    my $ftr_type_idx     = $ftr_fileroot_AR->[$ftr_idx];
    my $ftr_sqstring_alt = ""; # fetched from sqfile_for_cds_mp_alerts
    my $ftr_sqstring_out = ""; # fetched from sqfile_for_output_fastas
    my $ftr_sqstring_pv  = ""; # fetched from sqfile_for_pv
    my $ftr_seq_name = undef;
    my @ftr2org_pos_A = (); # [1..$ftr_pos..$ftr_len] original sequence position that corresponds to this position in the feature
    $ftr2org_pos_A[0] = -1; # invalid
    my $ftr_len = 0;
    my $ftr_strand = undef;
    my $ftr_start  = undef; # predicted start for the feature
    my $ftr_stop   = undef; # predicted stop  for the feature
    my $ftr_stop_c = undef; # corrected stop  for the feature, stays undef if no correction needed (no 'trc' or 'ext')
    my $ftr_ofile_key    = $mdl_name . ".pfa." . $ftr_idx;
    my $pv_ftr_ofile_key = $mdl_name . ".pfa." . $ftr_idx . ".pv";
    my $ftr_results_HR = \%{$ftr_results_HAHR->{$seq_name}[$ftr_idx]}; # for convenience
    # printf("in $sub_name, set ftr_results_HR to ftr_results_HAHR->{$seq_name}[$ftr_idx]\n");
    my $ftr_5ablen     = 0; # number of consecutive nt starting at ftr_start (on 5' end) that are ambiguous (commonly 0)
    my $ftr_3ablen     = 0; # number of consecutive nt ending   at ftr_stop  (on 3' end) that are ambiguous (commonly 0)
    my $ftr_5ablen_pv  = 0; # number of consecutive nt starting at ftr_start (on 5' end) that are ambiguous (commonly 0) in protein validation sqstring
    my $ftr_3ablen_pv  = 0; # number of consecutive nt ending   at ftr_stop  (on 3' end) that are ambiguous (commonly 0) in protein validation sqstring
    my $ftr_start_non_ab    = undef; # sequence position of first non-N on 5' end, commonly $ftr_start, -1 if complete feature is ambiguities
    my $ftr_stop_non_ab     = undef; # sequence position of first non-N on 3' end, commonly $ftr_stop, -1 if complete feature is ambiguities
    my $ftr_start_non_ab_pv = undef; # sequence position of first non-N on 5' end in protein validation sqstring, commonly $ftr_start, -1 if complete feature is ambiguities
    my $ftr_stop_non_ab_pv  = undef; # sequence position of first non-N on 3' end in protein validation sqstring, commonly $ftr_stop, -1 if complete feature is ambiguities
    my $ftr_scoords = undef; # coords string with sequence coordinates of all segments of the feature
    my $ftr_mcoords = undef; # coords string with model    coordinates of all segments of the feature

    my %alt_str_H = (); # added to as we find alerts below
    # ambgnt5c, ambgnt3c, ambgnt5f, ambgnt3f, mutstart, unexleng, mutendcd, mutendex, mutendns, cdsstopn
    my $alt_scoords; # sequence coordinates related to an alert
    my $alt_mcoords; # model    coordinates related to an alert
    my $alt_codon;   # codon string for an alert message

    # determine if this feature is 5' and/or 3' truncated
    # we do this outside the main loop since the logic is a bit complex:
    # - a feature is 5' truncated if:
    #   (A) its 5'-most segment with results is not the 5'-most segment of the feature
    #      (regardless of whether its 5'-most segment is truncated or not) 
    #   OR
    #   (B) its 5'-most feature is truncated
    # - and vice versa for 3' truncation
    #         
    my $ftr_is_5trunc = undef;
    my $ftr_is_3trunc = undef;
    my $first_sgm_idx = get_5p_most_sgm_idx_with_results($ftr_info_AHR, $sgm_results_HAHR, $ftr_idx, $seq_name);
    my $final_sgm_idx = get_3p_most_sgm_idx_with_results($ftr_info_AHR, $sgm_results_HAHR, $ftr_idx, $seq_name);
    if($first_sgm_idx != -1) { 
      $ftr_is_5trunc = 
          (($first_sgm_idx != $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}) || # (A) above
           ($sgm_results_HAHR->{$seq_name}[$first_sgm_idx]{"5trunc"}))
          ? 1 : 0;
      $ftr_is_3trunc = (($final_sgm_idx != $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"}) || # (A) above
                        ($sgm_results_HAHR->{$seq_name}[$final_sgm_idx]{"3trunc"})) 
          ? 1 : 0;
    }

    # main loop over segments
    for(my $sgm_idx = $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}; $sgm_idx <= $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"}; $sgm_idx++) { 
      if((defined $sgm_results_HAHR->{$seq_name}) && 
         (defined $sgm_results_HAHR->{$seq_name}[$sgm_idx]) && 
         (defined $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"})) { 
        my $sgm_results_HR = $sgm_results_HAHR->{$seq_name}[$sgm_idx]; # for convenience
        my ($sstart, $sstop, $sstrand) = ($sgm_results_HR->{"sstart"}, $sgm_results_HR->{"sstop"}, $sgm_results_HR->{"strand"});
        my ($mstart, $mstop)           = ($sgm_results_HR->{"mstart"}, $sgm_results_HR->{"mstop"});
        
        # only update start and 5trunc value if this is the first segment annotated
        if(! defined $ftr_start) { $ftr_start = $sstart; }
        $ftr_stop = $sstop; # always update $ftr_stop
        
        # set feature strand if this is the first segment annotated
        # else for cds/mp validate it hasn't changed and fail if it has
        # or update strand to "!" if not cds/mp and it has changed
        if(! defined $ftr_strand) { 
          $ftr_strand = $sstrand; 
        }
        elsif($ftr_strand ne $sstrand) { 
          # mixture of strands on different segments, this shouldn't happen if we're a CDS or mat_peptide
          if($ftr_is_cds_or_mp) { 
            # this 'shouldn't happen' for a CDS or mature peptide, all segments should be the sames strand
            ofile_FAIL("ERROR, in $sub_name, different model segments have different strands for a CDS or MP feature $ftr_idx", 1, undef);
          }
          # mixture of strands, set to "!" 
          $ftr_strand = "!";
        }
        
        # update $ftr_sqstring_alt, $ftr_sqstring_out, $ftr_sqstring_pv, $ftr_seq_name, $ftr_len, @ftr2org_pos_A, and @ftr2sgm_idx_A
        my $sgm_len = abs($sstop - $sstart) + 1;

        ###############
        # fetch the sequence string, only do this when nec, some files may be identical and we take advantage of that
        # alt
        my $cur_ftr_sqstring_alt = $sqfile_for_cds_mp_alerts->fetch_subseq_to_sqstring($seq_name, $sstart, $sstop, ($sstrand eq "-"));
        $ftr_sqstring_alt .= $cur_ftr_sqstring_alt;

        # out
        my $cur_ftr_sqstring_out = undef;
        if($sqfile_for_output_fastas_path eq $sqfile_for_cds_mp_alerts_path) { # no need to fetch again
          $ftr_sqstring_out .= $cur_ftr_sqstring_alt;
        }
        else { # fetch
          $cur_ftr_sqstring_out = $sqfile_for_output_fastas->fetch_subseq_to_sqstring($seq_name, $sstart, $sstop, ($sstrand eq "-"));
          $ftr_sqstring_out .= $cur_ftr_sqstring_out;
        }

        # pv
        if($sqfile_for_pv_path eq $sqfile_for_cds_mp_alerts_path) { # no need to fetch again
          $ftr_sqstring_pv .= $cur_ftr_sqstring_alt;
        }
        elsif($sqfile_for_pv_path eq $sqfile_for_output_fastas_path) { # no need to fetch again
          $ftr_sqstring_pv .= $cur_ftr_sqstring_out;
        }
        else { 
          $ftr_sqstring_pv .= $sqfile_for_pv->fetch_subseq_to_sqstring($seq_name, $sstart, $sstop, ($sstrand eq "-"));
        }
        ###############

        if(! defined $ftr_seq_name) { 
          $ftr_seq_name = $seq_name . "/" . $ftr_type_idx . "/"; 
        }
        else { 
          $ftr_seq_name .= ",";
        }
        $ftr_seq_name .= $sstart . ".." . $sstop . ":" . $sstrand;
        
        if($ftr_is_cds_or_mp) { 
          # update ftr2org_pos_A, if nec
          my $sgm_offset = 0;
          for(my $sgm_offset = 0; $sgm_offset < $sgm_len; $sgm_offset++) { 
            $ftr2org_pos_A[$ftr_len + $sgm_offset + 1] = ($sstrand eq "-") ? $sstart - $sgm_offset : $sstart + $sgm_offset;
            # slightly wasteful in certain cases, if $ftr_is_5trunc && $ftr_is_3trunc then we won't use this
          }
        }
        $ftr_scoords = vdr_CoordsAppendSegment($ftr_scoords, vdr_CoordsSegmentCreate($sstart, $sstop, $sstrand, $FH_HR));
        $ftr_mcoords = vdr_CoordsAppendSegment($ftr_mcoords, vdr_CoordsSegmentCreate($mstart, $mstop, $sstrand, $FH_HR));
        $ftr_len += $sgm_len;
      } # end of 'if(defined $sgm_results_HAHR->{$seq_name}...'
    } # end of 'for(my $sgm_idx = $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}...

    # printf("in $sub_name seq_name: $seq_name ftr_idx: $ftr_idx ftr_len: $ftr_len ftr_start: $ftr_start ftr_stop: $ftr_stop\n");
    if($ftr_len > 0) { 
      # we had a prediction for at least one of the segments for this feature

      # determine the position of the first and final N or n in ftr_sqstring_alt and ftr_sqstring_pv
      # we use ftr_sqstring_alt values for alerts
      # we use ftr_sqstring_pv  values later during protein validation to adjust protein/nucleotide difference tolerance at ends
      $ftr_5ablen = helper_feature_terminal_ambiguities($ftr_sqstring_alt, 0, $ftr_is_5trunc, # 0: not reversed
                                                        $ftr_start, $ftr_stop, $ftr_strand, $ftr_scoords, $ftr_len, 
                                                        $ftr_is_cds, $ftr_matches_cds, $mdl_tt, $atg_only, \%alt_str_H, $ua2rf_AR, $FH_HR);
      $ftr_start_non_ab = ($ftr_5ablen != $ftr_len) ? vdr_CoordsRelativeSingleCoordToAbsolute($ftr_scoords, ($ftr_5ablen + 1), $FH_HR) : -1;

      # same drill for ftr_sqstring_pv
      $ftr_5ablen_pv = helper_feature_terminal_ambiguities($ftr_sqstring_pv, 0, $ftr_is_5trunc, # 0: not reversed
                                                           $ftr_start, $ftr_stop, $ftr_strand, $ftr_scoords, $ftr_len, 
                                                           $ftr_is_cds, $ftr_matches_cds, $mdl_tt, $atg_only, undef, $ua2rf_AR, $FH_HR); # undef \%alt_str_H, don't report alerts for pv 

      my $rev_ftr_sqstring_alt = reverse($ftr_sqstring_alt);
      $ftr_3ablen = helper_feature_terminal_ambiguities($rev_ftr_sqstring_alt, 1, $ftr_is_3trunc, # 1: reversed
                                                        $ftr_start, $ftr_stop, $ftr_strand, $ftr_scoords, $ftr_len, 
                                                        $ftr_is_cds, $ftr_matches_cds, $mdl_tt, $atg_only, \%alt_str_H, $ua2rf_AR, $FH_HR);
      $ftr_stop_non_ab = ($ftr_3ablen != $ftr_len) ? vdr_CoordsRelativeSingleCoordToAbsolute($ftr_scoords, ($ftr_len - $ftr_3ablen), $FH_HR) : -1;

      # same drill for ftr_sqstring_pv
      my $rev_ftr_sqstring_pv = reverse($ftr_sqstring_pv);
      $ftr_3ablen_pv = helper_feature_terminal_ambiguities($rev_ftr_sqstring_pv, 1, $ftr_is_3trunc, # 1: reversed
                                                           $ftr_start, $ftr_stop, $ftr_strand, $ftr_scoords, $ftr_len, 
                                                           $ftr_is_cds, $ftr_matches_cds, $mdl_tt, $atg_only, undef, $ua2rf_AR, $FH_HR); # undef \%alt_str_H, don't report alerts for pv 

      # output the sequence
      if(! exists $ofile_info_HHR->{"FH"}{$ftr_ofile_key}) { 
        ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, $ftr_ofile_key,  $out_root . "." . $mdl_name . "." . $ftr_fileroot_AR->[$ftr_idx] . ".fa", ($do_allfasta ? 1 : 0), ($do_allfasta ? 1 : 0), "model $mdl_name feature " . $ftr_outroot_AR->[$ftr_idx] . " predicted seqs");
        if(! $do_allfasta) { 
          push(@{$to_remove_AR}, $ofile_info_HHR->{"fullpath"}{$ftr_ofile_key});
        }
      }
      print { $ofile_info_HHR->{"FH"}{$ftr_ofile_key} } (">" . $ftr_seq_name . "\n" . 
                                                         seq_SqstringAddNewlines($ftr_sqstring_out, 60)); 
      if(($do_separate_cds_fa_files) && ($ftr_is_cds)) { 
        if(! exists $ofile_info_HHR->{"FH"}{$pv_ftr_ofile_key}) { 
          my $separate_cds_fa_file = $out_root . "." . $mdl_name . "." . $ftr_fileroot_AR->[$ftr_idx] . ".pv.fa"; 
          ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, $pv_ftr_ofile_key, $separate_cds_fa_file, 0, $do_keep, "model $mdl_name feature " . $ftr_outroot_AR->[$ftr_idx] . " predicted seqs for protein validation");
          push(@{$to_remove_AR}, $separate_cds_fa_file);
        }
        print { $ofile_info_HHR->{"FH"}{$pv_ftr_ofile_key} } (">" . $ftr_seq_name . "\n" . 
                                                              seq_SqstringAddNewlines($ftr_sqstring_pv, 60)); 
      }
      # deal with mutstart for all CDS that are not 5' truncated
      if(! $ftr_is_5trunc) {
        # feature is not 5' truncated, look for a start codon if it's a CDS
        # and no ambgnt5c or ambgcd5c alert already reported
        if(($ftr_is_cds) && (! defined $alt_str_H{"ambgnt5c"}) && (! defined $alt_str_H{"ambgcd5c"})) { 
          if(($ftr_len >= 3) && (! sqstring_check_start($ftr_sqstring_alt, $mdl_tt, $atg_only, $FH_HR))) { 
            $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($ftr2org_pos_A[1], $ftr2org_pos_A[3], $ftr_strand, $FH_HR) . ";";
            $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[($ftr2org_pos_A[1])]), abs($ua2rf_AR->[($ftr2org_pos_A[3])]), $ftr_strand, $FH_HR) . ";";
            $alt_codon = substr($ftr_sqstring_alt, 0, 3);
            $alt_codon =~ tr/a-z/A-Z/;
            $alt_str_H{"mutstart"} = sprintf("%s%s%s", $alt_scoords, $alt_mcoords, $alt_codon);
          }
        }
      }
      # deal with mutendcd for all CDS that are not 3' truncated BUT are 5' truncated
      if((! $ftr_is_3trunc) && ($ftr_is_5trunc)) { 
        # feature is not 3' truncated, but it is 5' truncated, look for a stop codon if it's a CDS
        # and no ambgnt3c or ambgcd3c alert already reported
        if(($ftr_is_cds) && (! defined $alt_str_H{"ambgnt3c"}) && (! defined $alt_str_H{"ambgcd3c"})) {
          if(($ftr_len >= 3) && (! sqstring_check_stop($ftr_sqstring_alt, $mdl_tt, $FH_HR))) { 
            $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($ftr2org_pos_A[($ftr_len-2)], $ftr2org_pos_A[$ftr_len], $ftr_strand, $FH_HR) . ";";
            $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[($ftr2org_pos_A[($ftr_len-2)])]), abs($ua2rf_AR->[($ftr2org_pos_A[$ftr_len])]), $ftr_strand, $FH_HR) . ";";
            $alt_codon = substr($ftr_sqstring_alt, -3, 3);
            $alt_codon =~ tr/a-z/A-Z/;
            $alt_str_H{"mutendcd"} = sprintf("%s%s%s", $alt_scoords, $alt_mcoords, $alt_codon);
          }
        }
      }      

      # deal with all CDS that are not 5' truncated and not 3' truncated
      # check for unexleng alert for non-truncated CDS/mat_peptides
      if(($ftr_is_cds_or_mp) && (! $ftr_is_5trunc) && (! $ftr_is_3trunc)) { 
        # feature is not truncated on either end, look for stop codons
        if(($ftr_len % 3) != 0) { 
          # not a multiple of 3, unexleng alert 
          $alt_scoords  = "seq:" . $ftr_scoords . ";";
          $alt_mcoords  = "mdl:" . $ftr_mcoords . ";"; 
          $alt_str_H{"unexleng"} = $alt_scoords . $alt_mcoords . $ftr_len;
        }
      }

      # if CDS: look for all valid in-frame stops 
      if($ftr_is_cds) { 
        if(! defined $ftr_results_HR->{"n_codon_start_expected"}) { 
          ofile_FAIL("ERROR in $sub_name, sequence $seq_name CDS feature (ftr_idx: $ftr_idx) has no codon_start info", 1, $FH_HR);
        }
        my $cds_codon_start = (defined $ftr_results_HR->{"n_codon_start_expected"}) ? $ftr_results_HR->{"n_codon_start_expected"} : 1;
        # make a new sqstring $ftr_sqstring_alt_stops from $ftr_sqstring_alt to search for stops in that:
        # - is identical to $ftr_sqstring_alt                   if codon_start == 1
        # - has the first     nt removed from $ftr_sqstring_alt if codon_start == 2
        # - has the first two nt removed from $ftr_sqstring_alt if codon_start == 3
        #my $n_nt_skipped_at_5p_end = ($ftr_is_5trunc) ? ($ftr_results_HR->{"n_codon_start_expected"} - 1) : 0;
        my $n_nt_skipped_at_5p_end = $cds_codon_start - 1;
        my $ftr_sqstring_alt_stops = substr($ftr_sqstring_alt, $n_nt_skipped_at_5p_end);
        my $ftr_len_stops = length($ftr_sqstring_alt_stops);
        if($ftr_len_stops >= 3) { 
          # check for mutendcd alert (final 3 nt are a valid stop) if ! 3' truncated
          if((! $ftr_is_3trunc) && (! sqstring_check_stop($ftr_sqstring_alt_stops, $mdl_tt, $FH_HR)) && (! defined $alt_str_H{"ambgnt3c"}) && (! defined $alt_str_H{"ambgcd3c"})) { 
            $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($ftr2org_pos_A[($ftr_len_stops + $n_nt_skipped_at_5p_end - 2)], $ftr2org_pos_A[($ftr_len_stops + $n_nt_skipped_at_5p_end)], $ftr_strand, $FH_HR) . ";";
            $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[($ftr2org_pos_A[($ftr_len_stops + $n_nt_skipped_at_5p_end - 2)])]), abs($ua2rf_AR->[($ftr2org_pos_A[$ftr_len_stops + $n_nt_skipped_at_5p_end])]), $ftr_strand, $FH_HR) . ";";
            $alt_codon = substr($ftr_sqstring_alt_stops, -3, 3);
            $alt_codon =~ tr/a-z/A-Z/;
            $alt_str_H{"mutendcd"} = sprintf("%s%s%s", $alt_scoords, $alt_mcoords, $alt_codon);
          }
          # look for cdsstopn, mutendex and mutendns ONLY IF 
          # CDS is NOT 5' truncated OR 
          # n_codon_start_expected == n_codon_start_dominant
          if((! $ftr_is_5trunc) || 
             ((defined $ftr_results_HR->{"n_codon_start_expected"}) && 
              (defined $ftr_results_HR->{"n_codon_start_dominant"}) && 
              ($ftr_results_HR->{"n_codon_start_expected"} == $ftr_results_HR->{"n_codon_start_dominant"}))) { 
            my @ftr_nxt_stp_A = ();
            sqstring_find_stops($ftr_sqstring_alt_stops, $mdl_tt, \@ftr_nxt_stp_A, $FH_HR);
            if($ftr_nxt_stp_A[1] != $ftr_len_stops) { 
              # first in-frame stop codon 3' of $ftr_start is not $ftr_stop
              # We will need to add an alert, (exactly) one of:
              # 'mutendex': no stop exists in $ftr_sqstring_alt_stops, but one does 3' of end of $ftr_sqstring_alt_stops
              # 'mutendns': no stop exists in $ftr_sqstring_alt_stops, and none exist 3' of end of $ftr_sqstring_alt_stops either
              # 'cdsstopn': an early stop exists in $ftr_sqstring_alt_stops
              if((! $ftr_is_3trunc) && ($ftr_nxt_stp_A[1] == 0)) { 
                # there are no valid in-frame stops in $ftr_sqstring_alt_stops
                # if we are not 3' truncated then we have a 'mutendns' or 'mutendex' alert, to find out which 
                # we need to fetch the sequence ending at $fstop to the end of the sequence 
                if((($ftr_strand eq "+") && ($ftr_stop < $seq_len)) ||
                   (($ftr_strand eq "-") && ($ftr_stop > 1))) { 
                  # we have some sequence left 3' of ftr_stop
                  # *careful* we don't always want to fetch starting at next nt after predicted 
                  # stop, because we have to make sure 1st position of $ext_sqstring is frame 1
                  # because sqstring_find_stops() expects this.
                  my $ext_sqstring = undef;
                  my $ext_sqstring_start = undef;
                  if($ftr_strand eq "+") { 
                    if   (($ftr_len_stops % 3) == 0) { $ext_sqstring_start = $ftr_stop+1; } # first in-frame stop can start at next posn
                    elsif(($ftr_len_stops % 3) == 1) { $ext_sqstring_start = $ftr_stop;   } # first in-frame stop can start at final posn
                    elsif(($ftr_len_stops % 3) == 2) { $ext_sqstring_start = $ftr_stop-1; } # first in-frame stop can start at prev posn
                    $ext_sqstring = $sqfile_for_cds_mp_alerts->fetch_subseq_to_sqstring($seq_name, $ext_sqstring_start, $seq_len, 0); 
                  }
                  else { # negative strand
                    if   (($ftr_len_stops % 3) == 0) { $ext_sqstring_start = $ftr_stop-1; } # first in-frame stop can start at next posn
                    elsif(($ftr_len_stops % 3) == 1) { $ext_sqstring_start = $ftr_stop;   } # first in-frame stop can start at final posn
                    elsif(($ftr_len_stops % 3) == 2) { $ext_sqstring_start = $ftr_stop+1; } # first in-frame stop can start at prev posn
                    $ext_sqstring = $sqfile_for_cds_mp_alerts->fetch_subseq_to_sqstring($seq_name, $ext_sqstring_start, 1, 1);
                  }
                  my @ext_nxt_stp_A = ();
                  sqstring_find_stops($ext_sqstring, $mdl_tt, \@ext_nxt_stp_A, $FH_HR);
                  if($ext_nxt_stp_A[1] != 0) { 
                    # there is an in-frame stop codon, mutendex alert
                    if((! defined $alt_str_H{"ambgnt3c"}) && (! defined $alt_str_H{"ambgcd3c"})) { # report it only if !ambgnt3c and !ambgcd3c
                      # determine what position it is
                      $ftr_stop_c = ($ftr_strand eq "+") ? ($ext_sqstring_start + ($ext_nxt_stp_A[1] - 1)) : ($ext_sqstring_start - ($ext_nxt_stp_A[1] - 1));
                      if($ftr_strand eq "+") { 
                        $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($ftr_stop_c-2, $ftr_stop_c, $ftr_strand, $FH_HR) . ";";
                        $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[($ftr_stop_c-2)]), abs($ua2rf_AR->[$ftr_stop_c]), $ftr_strand, $FH_HR) . ";";
                      }
                      else {
                        $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($ftr_stop_c+2, $ftr_stop_c, $ftr_strand, $FH_HR) . ";";
                        $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[($ftr_stop_c+2)]), abs($ua2rf_AR->[$ftr_stop_c]), $ftr_strand, $FH_HR) . ";";
                      }
                      $alt_codon = substr($ext_sqstring, ($ext_nxt_stp_A[1]-3), 3);
                      $alt_codon =~ tr/a-z/A-Z/;
                      $alt_str_H{"mutendex"} = sprintf("%s%s%s", $alt_scoords, $alt_mcoords, $alt_codon);
                    }
                  }
                } # end of 'if((($ftr_strand eq "+") && ($ftr_stop < $seq_len)) || ($ftr_strand eq "-") && ($ftr_stop > 1)))'
                if(! defined $ftr_stop_c) { 
                  # if we get here, either $ftr_stop == $seq_len (and there was no more seq to check for a stop codon)
                  # or we checked the sequence but didn't find any
                  # either way, we have a mutendns alert:
                  if((! defined $alt_str_H{"ambgnt3c"}) && (! defined $alt_str_H{"ambgcd3c"})) { # report it only if !ambgnt3c and !ambgcd3c
                    $ftr_stop_c = "?"; # special case, we don't know where the stop is, but we know it's not $ftr_stop;
                    # we don't provide scoords or mcoords for this alert
                    $alt_str_H{"mutendns"} = "VADRNULL";
                  }
                }
              } # end of 'if((! $ftr_is_3trunc) && ($ftr_nxt_stp_A[1] == 0) {' 
              elsif($ftr_nxt_stp_A[1] != 0) { 
                # there is an early stop (cdsstopn) in $ftr_sqstring_alt_stops
                if($ftr_nxt_stp_A[1] > $ftr_len_stops) { 
                  # this shouldn't happen, it means there's a bug in sqstring_find_stops()
                  ofile_FAIL("ERROR, in $sub_name, problem identifying stops in feature sqstring for ftr_idx $ftr_idx, found a stop at position that exceeds feature length", 1, undef);
                }
                $ftr_stop_c = $ftr2org_pos_A[($ftr_nxt_stp_A[1] + $n_nt_skipped_at_5p_end)];
                if($ftr_strand eq "+") { 
                  $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($ftr_stop_c-2, $ftr_stop_c, $ftr_strand, $FH_HR) . ";";
                  $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[($ftr_stop_c-2)]), abs($ua2rf_AR->[$ftr_stop_c]), $ftr_strand, $FH_HR) . ";";
                }
                else {
                  $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($ftr_stop_c+2, $ftr_stop_c, $ftr_strand, $FH_HR) . ";";
                  $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[($ftr_stop_c+2)]), abs($ua2rf_AR->[$ftr_stop_c]), $ftr_strand, $FH_HR) . ";";
                }
                $alt_codon = substr($ftr_sqstring_alt_stops, $ftr_nxt_stp_A[1]-3, 3);
                $alt_codon =~ tr/a-z/A-Z/;
                $alt_str_H{"cdsstopn"} = sprintf("%s%s%s, shifted S:%d,M:%d", $alt_scoords, $alt_mcoords, $alt_codon, abs($ftr_stop-$ftr_stop_c), abs(abs($ua2rf_AR->[$ftr_stop]) - abs($ua2rf_AR->[$ftr_stop_c])));
              } # end of 'elsif($ftr_nxt_stp_A[1] != 0)'
            } # end of 'if($ftr_nxt_stp_A[1] != $ftr_len_stops) {' 
          } # end of big if entered if CDS is not truncated or dominant frame and expected frame are identical
        } # end of 'if($ftr_len_stops >= 3)'
      } # end of 'if($ftr_is_cds) {' 
    
      # actually add the alerts
      foreach my $alt_code (sort keys %alt_str_H) { 
        alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, $alt_code, $seq_name, $ftr_idx, $alt_str_H{$alt_code}, $FH_HR);
      }
      
      # if we are a mature peptide, make sure we are adjacent to the next one, if there is one
      if($ftr_is_mp && ($ftr_info_AHR->[$ftr_idx]{"3pa_ftr_idx"} != -1)) { 
        my $ftr_3pa_idx = $ftr_info_AHR->[$ftr_idx]{"3pa_ftr_idx"};
        my $sgm_5p_idx  = $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"};     
        # yes, the 3'-most segment of $ftr_idx is the 5'-most mature peptide we are interested in
        my $sgm_3p_idx  = $ftr_info_AHR->[$ftr_3pa_idx]{"5p_sgm_idx"}; 
        # and, yes, the 5'-most segment of the 3' adjacent $ftr_idx is the 3'-most mature peptide we're interested in
        
        my $sgm_5p_valid = ((defined $sgm_results_HAHR->{$seq_name}) && 
                            (defined $sgm_results_HAHR->{$seq_name}[$sgm_5p_idx]) && 
                            (defined $sgm_results_HAHR->{$seq_name}[$sgm_5p_idx]{"sstart"})) ? 1 : 0;
        my $sgm_3p_valid = ((defined $sgm_results_HAHR->{$seq_name}) && 
                            (defined $sgm_results_HAHR->{$seq_name}[$sgm_3p_idx]) && 
                            (defined $sgm_results_HAHR->{$seq_name}[$sgm_3p_idx]{"sstart"})) ? 1 : 0;
        my $sgm_5p_3flush = ($sgm_5p_valid && $sgm_results_HAHR->{$seq_name}[$sgm_5p_idx]{"3seqflush"}) ? 1 : 0;
        my $sgm_3p_5flush = ($sgm_3p_valid && $sgm_results_HAHR->{$seq_name}[$sgm_3p_idx]{"5seqflush"}) ? 1 : 0;
        
        my $stop_5p  = ($sgm_5p_valid) ? $sgm_results_HAHR->{$seq_name}[$sgm_5p_idx]{"sstop"}  : undef;
        my $start_3p = ($sgm_3p_valid) ? $sgm_results_HAHR->{$seq_name}[$sgm_3p_idx]{"sstart"} : undef;
        
        # Three ways we can get a 'pepadjcy' alert: 
        if($sgm_5p_valid && $sgm_3p_valid) { # both are valid 
          if((abs($stop_5p - $start_3p)) != 1) { # they're not adjacent
            # 1) both mature peptides are annotated but not adjacent, alert on $ftr_idx
            $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate((($ftr_strand eq "+") ? $stop_5p + 1 : $stop_5p - 1), (($ftr_strand eq "+") ? $start_3p - 1 : $start_3p + 1), $ftr_strand, $FH_HR) . ";";
            $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$stop_5p]), abs($ua2rf_AR->[$start_3p]), $ftr_strand, $FH_HR) . ";";
            alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, "pepadjcy", $seq_name, $ftr_idx, 
                                       sprintf("%s%sVADRNULL", $alt_scoords, $alt_mcoords), 
                                       $FH_HR);
          }
        }
        elsif(($sgm_5p_valid) && (! $sgm_3p_valid) && (! $sgm_5p_3flush)) { 
          # 2) 5' mature peptide is annotated and ends before end of sequence, but 3' mature peptide is not annotated, alert for $ftr_idx
          $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($stop_5p, $stop_5p, $ftr_strand, $FH_HR) . ";";
          $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$stop_5p]), abs($ua2rf_AR->[$stop_5p]), $ftr_strand, $FH_HR) . ";";
          alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, "pepadjcy", $seq_name, $ftr_idx, 
                                     sprintf("%s%sfeature stops at seq position $stop_5p on %s strand which is not terminal but expected 3'-adjacent feature is not annotated", $alt_scoords, $alt_mcoords, $sgm_results_HAHR->{$seq_name}[$sgm_5p_idx]{"strand"}),
                                     $FH_HR);
        }
        elsif(($sgm_3p_valid) && (! $sgm_5p_valid) && (! $sgm_3p_5flush)) { 
          # 3) 3' mature peptide is annotated and starts after start of sequence, but 5' mature peptide is not annotated, alert for $ftr_3pa_idx (NOT $ftr_idx)
          $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($start_3p, $start_3p, $ftr_strand, $FH_HR) . ";";
          $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$start_3p]), abs($ua2rf_AR->[$start_3p]), $ftr_strand, $FH_HR) . ";";
          alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, "pepadjcy", $seq_name, $ftr_3pa_idx, 
                                     sprintf("%s%sfeature starts at seq position $start_3p on %s strand which is not terminal but expected 5'-adjacent feature is not annotated", $alt_scoords, $alt_mcoords, $sgm_results_HAHR->{$seq_name}[$sgm_3p_idx]{"strand"}),
                                     $FH_HR);
        }
      } # end of 'if($ftr_is_mp && ($ftr_info_AHR->[$ftr_idx]{"3pa_ftr_idx"} != -1))'

      # update %ftr_results_HR
      $ftr_results_HR->{"n_strand"}          = $ftr_strand;
      $ftr_results_HR->{"n_start"}           = $ftr_start;
      $ftr_results_HR->{"n_stop"}            = $ftr_stop;
      $ftr_results_HR->{"n_stop_c"}          = (defined $ftr_stop_c) ? $ftr_stop_c : $ftr_stop;
      $ftr_results_HR->{"n_5trunc"}          = $ftr_is_5trunc;
      $ftr_results_HR->{"n_3trunc"}          = $ftr_is_3trunc;
      $ftr_results_HR->{"n_5ablen"}          = $ftr_5ablen;
      $ftr_results_HR->{"n_3ablen"}          = $ftr_3ablen;
      $ftr_results_HR->{"n_5ablen_pv"}       = $ftr_5ablen_pv;
      $ftr_results_HR->{"n_3ablen_pv"}       = $ftr_3ablen_pv;
      $ftr_results_HR->{"n_start_non_ab"}    = $ftr_start_non_ab;
      $ftr_results_HR->{"n_stop_non_ab"}     = $ftr_stop_non_ab;
      $ftr_results_HR->{"n_start_non_ab_pv"} = $ftr_start_non_ab_pv;
      $ftr_results_HR->{"n_stop_non_ab_pv"}  = $ftr_stop_non_ab_pv;
      $ftr_results_HR->{"n_len"}             = $ftr_len;
      $ftr_results_HR->{"n_scoords"}         = $ftr_scoords;
      $ftr_results_HR->{"n_mcoords"}         = $ftr_mcoords;
    } # end of 'if($ftr_len > 0)'
  } # end of 'for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { '
  
  return;
}


#################################################################
# Subroutine: sqstring_check_start()
# Incept:     EPN, Sat Feb 23 10:18:22 2019
#
# Arguments:
#  $sqstring: the sequence string
#  $tt:       the translation table ('1' for standard)
#  $atg_only: only allow 'ATG' for start, regardless of translation table 
#  $FH_HR:    REF to hash of file handles
#  
# Returns: '1' if $sqstring starts with a valid
#           start codon on the positive strand
#           '0' if not
# 
#################################################################
sub sqstring_check_start {
  my $sub_name = "sqstring_check_start";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqstring, $tt, $atg_only, $FH_HR) = @_;

  my $sqlen = length($sqstring);
  if($sqlen < 3) { return 0; } 

  my $start_codon = substr($sqstring, 0, 3);
  $start_codon =~ tr/a-z/A-Z/; # convert to uppercase
  $start_codon =~ tr/U/T/;     # convert to DNA

  return seq_CodonValidateStartCapDna($start_codon, $tt, $atg_only);
}

#################################################################
# Subroutine: sqstring_check_stop()
# Incept:     EPN, Wed Jan 22 14:47:37 2020
#
# Arguments:
#  $sqstring: the sequence string
#  $tt:       the translation table ('1' for standard)
#  $FH_HR:    REF to hash of file handles
#  
# Returns: '1' if $sqstring ends with a valid
#           stop codon on the positive strand
#           '0' if not
# 
#################################################################
sub sqstring_check_stop {
  my $sub_name = "sqstring_check_stop";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqstring, $tt, $FH_HR) = @_;

  my $sqlen = length($sqstring);
  if($sqlen < 3) { return 0; } 

  my $stop_codon = substr($sqstring, -3, 3);
  $stop_codon =~ tr/a-z/A-Z/; # convert to uppercase
  $stop_codon =~ tr/U/T/;     # convert to DNA

  return seq_CodonValidateStopCapDna($stop_codon, $tt);

}

#################################################################
# Subroutine: sqstring_find_stops()
# Incept:     EPN, Fri Feb 22 14:52:53 2019
#
# Purpose:   Find all occurences of stop codons in an 
#            input sqstring (on the positive strand only),
#            and update the input array.
#           
#            This subroutine could be easily modified to 
#            find the nearest stop codon in each frame, 
#            but currently it is only concerned with 
#            frame 1.           
#
#            We determine what frame each position is by asserting
#            that the first position is frame 1 (the frame for the
#            first position of a start codon). Previously, up until
#            version 1.1.3 we assumed the final position was frame 
#            3 but this led to uninformative error messages related to 
#            stop codons.
#
# Arguments:
#  $sqstring:       the sequence string
#  $tt:             the translation table ('1' for standard)
#  $nxt_stp_AR:     [1..$i..$sqlen] = $x; closest stop codon at or 3' of position
#                   $i in frame 1 on positive strand *ends* at position $x; 
#                   '0' if there are none.
#                   special values: 
#                   $nxt_stp_AR->[0] = -1
#  $FH_HR:          REF to hash of file handles
#             
# Returns:  void, updates arrays that are not undef
# 
# Dies:     never
#
#################################################################
sub sqstring_find_stops { 
  my $sub_name = "sqstring_find_stops";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqstring, $tt, $nxt_stp_AR, $FH_HR) = @_;
  
  @{$nxt_stp_AR} = ();
  $nxt_stp_AR->[0] = -1;

  # create the @sqstring_A we'll use to find the stops
  $sqstring =~ tr/a-z/A-Z/; # convert to uppercase
  $sqstring =~ tr/U/T/;     # convert to DNA
  my $sqlen = length($sqstring);
  # add a " " as first character of $sqstart so we get nts in elements 1..$sqlen 
  my @sqstring_A = split("", " " . $sqstring);

  my $i       = undef; # [1..$i..$sqlen]: sequence position
  my $frame   = undef; # 1, 2, or 3
  my $cstart  = undef; # [1..$sqlen] start position of current codon
  my $codon   = undef; # the codon
  my $cur_stp = 0;     # closest stop codon at or 3' of current position
                       # in frame 1 *ends* at position $cur_stp;

  # pass over sequence from right to left, filling @{$nxt_stp_AR}
  $cur_stp = 0;
  for($i = ($sqlen-2); $i >= 1; $i--) { 
    if(($i % 3) == 1) { # starting position of a codon, frame == 1
      $codon = $sqstring_A[$i] . $sqstring_A[($i+1)] . $sqstring_A[($i+2)];
      if(seq_CodonValidateStopCapDna($codon, $tt)) { 
        $cur_stp = $i+2;
      }
    }
    $nxt_stp_AR->[$i] = $cur_stp;
  }
  $nxt_stp_AR->[($sqlen-1)] = 0;
  $nxt_stp_AR->[$sqlen]     = 0;

#  for($i = 1; $i <= $sqlen; $i++) { 
#    printf("position %5d: nxt_stp: %5d\n", $i, $nxt_stp_AR->[$i]);
#  }

  return;
}

#################################################################
# Subroutine: add_low_similarity_alerts_for_one_sequence()
# Incept:    EPN, Thu Jun  3 14:58:10 2021 
#            EPN, Mon Apr 29 13:29:37 2019 [version for an array of seqs]
#
# Purpose:   For a given sequence, check for and report any of the following alerts:
#            low similarity per-sequence alerts          (lowsim5s, lowsim3s, lowsimis) and
#            low similarity per-coding-feature alerts    (lowsim5c, lowsim3c, lowsimic) and 
#            low similarity per-noncoding-feature alerts (lowsim5n, lowsim3n, lowsimin).
#
# Arguments:
#  $seq_name:               name of the sequence
#  $seq_len_HR:             REF to hash of sequence lengths, PRE-FILLED
#  $ua2rf_AR:               REF to array that maps unaligned positions to reference positions
#                           [1..$uapos..$ualen]: reference position that unaligned position $uapos aligns to 
#                           if $ua2rf_A[$uapos] <  0, $uapos inserts *after* ref posn (-1 * $ua2rf_A[$uapos])
#                           if $ua2rf_A[$uapos] == 0, $uapos inserts *before* ref posn 1
#                           $ua2rf_A[0] is invalid (set to 0)
#                           can be undef, if so, we just don't report model coords in alert instances
#  $ftr_info_AHR:           REF to hash of arrays with information on the features, PRE-FILLED
#  $sgm_info_AHR:           REF to hash of arrays with information on the model segments, PRE-FILLED
#  $alt_info_HHR:           REF to the alert info hash of arrays, PRE-FILLED
#  $stg_results_HHHR:       REF to 3D hash of classification results, PRE-FILLED
#  $sgm_results_HAHR:       REF to model segment results HAH, pre-filled
#  $ftr_results_HAHR:       REF to feature results HAH, added to here
#  $alt_seq_instances_HHR:  REF to array of hash with per-sequence alerts, PRE-FILLED
#  $alt_ftr_instances_HHHR: REF to array of 2D hashes with per-feature alerts, PRE-FILLED
#  $opt_HHR:                REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:         REF to the 2D hash of output file information
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub add_low_similarity_alerts_for_one_sequence { 
  my $sub_name = "add_low_similarity_alerts_for_one_sequence";
  my $nargs_exp = 13;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $seq_len_HR, $ua2rf_AR, $ftr_info_AHR, $sgm_info_AHR, $alt_info_HHR, 
      $stg_results_HHHR, $sgm_results_HAHR, $ftr_results_HAHR, $alt_seq_instances_HHR, $alt_ftr_instances_HHHR, 
      $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"}; # for convenience

  my $nftr = scalar(@{$ftr_info_AHR});

  my $terminal_seq_5_min_length = opt_Get("--lowsim5seq", $opt_HHR); # minimum length of terminal missing region that triggers a lowsim5s alert
  my $terminal_seq_3_min_length = opt_Get("--lowsim3seq", $opt_HHR); # minimum length of terminal missing region that triggers a lowsim3s alert
  my $internal_seq_min_length   = opt_Get("--lowsimiseq", $opt_HHR); # minimum length of internal missing region that triggers a lowsimis alert
  my $terminal_ftr_5_min_length = opt_Get("--lowsim5ftr", $opt_HHR); # minimum length of terminal missing region in a feature that triggers a lowsim5f alert
  my $terminal_ftr_3_min_length = opt_Get("--lowsim3ftr", $opt_HHR); # minimum length of terminal missing region in a feature that triggers a lowsim3f alert
  my $internal_ftr_min_length   = opt_Get("--lowsimiftr", $opt_HHR); # minimum length of internal missing region in a feature that triggers a lowsimif alert

  my $alt_msg        = undef; # message for reporting an alert
  my $alt_scoords    = undef; # sequence coords string for alert message
  my $alt_mcoords    = undef; # model    coords string for alert message
  my $soverlap_start = undef; # sequence position of start of overlap region
  my $soverlap_stop  = undef; # sequence position of end   of overlap region

  # set $min_length as minimum of the 5 length thresholds
  my $min_length = $internal_ftr_min_length;
  if($min_length > $internal_seq_min_length)   { $min_length = $internal_seq_min_length; }
  if($min_length > $terminal_ftr_5_min_length) { $min_length = $terminal_ftr_5_min_length; }
  if($min_length > $terminal_ftr_3_min_length) { $min_length = $terminal_ftr_3_min_length; }
  if($min_length > $terminal_seq_5_min_length) { $min_length = $terminal_seq_5_min_length; }
  if($min_length > $terminal_seq_3_min_length) { $min_length = $terminal_seq_3_min_length; }

  my $seq_len  = $seq_len_HR->{$seq_name};
  # determine number of nucleotides not covered by r2.bs search stage 
  # at 5' and 3' ends
  if((defined $stg_results_HHHR->{$seq_name}) && 
     (defined $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}) && 
     (defined $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"score"})) { 
    my @tmp_A = split(",", $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"score"}); # only do this to get nhits
    my $nhits = scalar(@tmp_A); 
    my $missing_coords = "";
    my $bstrand = $stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"bstrand"};
    if($nhits == 1) { 
      # only 1 hit
      my $min_coord = vdr_CoordsMin($stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"}, $FH_HR);
      my $max_coord = vdr_CoordsMax($stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"}, $FH_HR);
      if($min_coord != 1) { 
        if($bstrand eq "+") { $missing_coords = vdr_CoordsSegmentCreate(1, $min_coord-1, "+", $FH_HR); }
        else                { $missing_coords = vdr_CoordsSegmentCreate($min_coord-1, 1, "-", $FH_HR); }
      }
      if($max_coord != $seq_len) { 
        if($missing_coords ne "") { $missing_coords .= ","; }
        if($bstrand eq "+") { $missing_coords .= vdr_CoordsSegmentCreate($max_coord+1, $seq_len, "+", $FH_HR); }
        else                { $missing_coords .= vdr_CoordsSegmentCreate($seq_len, $max_coord+1, "-", $FH_HR); }
      }
    }
    else { 
      # multiple hits
      $missing_coords .= vdr_CoordsMissing($stg_results_HHHR->{$seq_name}{"std.cdt.bs"}{"s_coords"}, $bstrand, $seq_len, $FH_HR);
    }
    if($missing_coords ne "") { 
      my @missing_coords_A = split(",", $missing_coords);
      foreach my $missing_coords_tok (@missing_coords_A) { 
        my ($start, $stop, undef) = vdr_CoordsSegmentParse($missing_coords_tok, $FH_HR);
        my $length = abs($start - $stop) + 1;
        if($length >= $min_length) { 
          # length is greater than the minimum of all alert length reporting thresholds
          if($bstrand eq "+") { 
            my $is_start   = ($start == 1)        ? 1 : 0;
            my $is_end     = ($stop  == $seq_len) ? 1 : 0;
            # does this overlap with a feature by at least minimum overlap length threshold? 
            my $ftr_overlap_flag = 0;
            for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
              # determine if this feature qualifies as a 'coding' feature for purposes of the alert
              # it does if it is a CDS, mat_peptide or has identical coords to a CDS or mat_peptide
              # - feature is a CDS or mat_peptide OR has identical coordinates to a CDS or mat_peptide
              my $ftr_matches_coding = vdr_FeatureTypeIsCdsOrMatPeptideOrIdStartAndStop($ftr_info_AHR, $ftr_idx);
              my $ftr_results_HR = $ftr_results_HAHR->{$seq_name}[$ftr_idx]; # for convenience
              if((defined $ftr_results_HR->{"n_start"}) || (defined $ftr_results_HR->{"p_qstart"})) { 
                my $f_start  = (defined $ftr_results_HR->{"n_start"}) ? $ftr_results_HR->{"n_start"}  : $ftr_results_HR->{"p_qstart"};
                my $f_stop   = (defined $ftr_results_HR->{"n_start"}) ? $ftr_results_HR->{"n_stop"}   : $ftr_results_HR->{"p_qstop"};
                my $f_strand = (defined $ftr_results_HR->{"n_start"}) ? $ftr_results_HR->{"n_strand"} : $ftr_results_HR->{"p_strand"};
                if($f_strand eq $bstrand) { 
                  my $noverlap = undef;
                  my $overlap_reg = "";
                  my $start1 = utl_Min($start,   $stop);
                  my $stop1  = utl_Max($start,   $stop);
                  my $start2 = utl_Min($f_start, $f_stop);
                  my $stop2  = utl_Max($f_start, $f_stop);
                  ($noverlap, $overlap_reg) = seq_Overlap($start1, $stop1, $start2, $stop2, $FH_HR);
                  if($noverlap > 0) { 
                    #printf("is_start: $is_start, is_end: $is_end, length: $length\n");
                    # for 5'/3'/internal cases: only actually report an alert for non-CDS and non-MP features
                    # because CDS and MP are independently validated by blastx (unless --pv_skip)
                    my ($soverlap_start, $soverlap_stop, $alt_scoords, $alt_mcoords);
                    if($overlap_reg =~ /^(\d+)\-(\d+)$/) { 
                      $soverlap_start = ($f_strand eq "+") ? $1 : $2;
                      $soverlap_stop  = ($f_strand eq "+") ? $2 : $1;
                    }
                    else { 
                      ofile_FAIL("ERROR, in $sub_name, unable to parse overlap region: $overlap_reg", 1, $FH_HR);
                    }
                    $alt_scoords = "seq:" . vdr_CoordsSegmentCreate($soverlap_start, $soverlap_stop, $f_strand, $FH_HR) . ";"; 
                    if(defined $ua2rf_AR) { 
                      $alt_mcoords = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$soverlap_start]), abs($ua2rf_AR->[$soverlap_stop]), $f_strand, $FH_HR) . ";"; 
                    }
                    else { 
                      $alt_mcoords = "mdl:VADRNULL;";
                    }
                    $alt_msg = sprintf("%s%s%d nt overlap b/t low similarity region of length %d (%d..%d) and annotated feature (%d..%d)",
                                       $alt_scoords, $alt_mcoords, $noverlap, $length, $start, $stop, $f_start, $f_stop);
                    if(($is_start) && ($noverlap >= $terminal_ftr_5_min_length)) { 
                      $ftr_overlap_flag = 1;
                      alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, ($ftr_matches_coding ? "lowsim5c" : "lowsim5n"), $seq_name, $ftr_idx, $alt_msg, $FH_HR);
                    }
                    if(($is_end) && ($noverlap >= $terminal_ftr_3_min_length)) { 
                      $ftr_overlap_flag = 1;
                      alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, ($ftr_matches_coding ? "lowsim3c" : "lowsim3n"), $seq_name, $ftr_idx, $alt_msg, $FH_HR);
                    }
                    if((! $is_start) && (! $is_end) && ($noverlap >= $internal_ftr_min_length)) { 
                      $ftr_overlap_flag = 1;
                      alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, ($ftr_matches_coding ? "lowsimic" : "lowsimin"), $seq_name, $ftr_idx, $alt_msg, $FH_HR);
                    }
                  }
                }
              }
            } # end of 'for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++)'
            if(! $ftr_overlap_flag) { # no features overlapped above length threshold, potentially report lowsim5s, lowsim3s, or lowsimis
              $alt_scoords = "seq:" . vdr_CoordsSegmentCreate($start, $stop, "+", $FH_HR) . ";"; 
              if(defined $ua2rf_AR) { 
                $alt_mcoords = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$start]), abs($ua2rf_AR->[$stop]), "+", $FH_HR) . ";"; 
              }
              else { 
                $alt_mcoords = "mdl:VADRNULL;";
              }
              $alt_msg = sprintf("%s%slow similarity region of length %d", 
                                 $alt_scoords, $alt_mcoords, $length);
              if(($is_start) && ($length >= $terminal_seq_5_min_length)) { 
                alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "lowsim5s", $seq_name, $alt_msg, $FH_HR);
              }
              if(($is_end) && ($length >= $terminal_seq_3_min_length)) { 
                alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "lowsim3s", $seq_name, $alt_msg, $FH_HR);
              }
              if((! $is_start) && (! $is_end) && ($length >= $internal_seq_min_length)) { 
                alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "lowsimis", $seq_name, $alt_msg, $FH_HR);
              }
            }
          }
        } # end of 'if($length >= $min_length)'
      } # end of 'foreach my $missing_coords_tok (@missing_coords_A)'
    }
  }
  return;
}

#################################################################
#
# Subroutines related to blastx:
# add_protein_validation_alerts
# run_blastx_and_summarize_output
# parse_blastx_results 
# helper_protein_validation_breakdown_source
# helper_blastx_breakdown_max_indel_str
# helper_protein_validation_db_seqname_to_ftr_idx 
#
#################################################################

#################################################################
# Subroutine: make_protein_validation_fasta_file()
# Incept:     EPN, Wed Mar 18 12:01:36 2020
#
# Purpose:    Create a fasta file to use as input for blastx or 
#             hmmsearch in the protein validation stage.
#             This file consists of the full length sequences 
#             in the file $mdl_fa_file (these are typically the
#             set of sequences that match best to a specific model)
#             and the predicted CDS sequences for each sequence
#             in $mdl_fa_file.
#      
#             If $do_pv_blastx is '1' We remove the descriptions from the
#             sequences in $mdl_fa_file because blastx parsing is more
#             difficult if descriptions are included.
#
# Arguments: 
#  $out_fa_file:              name of output fasta file to create
#  $mdl_name:                 name of model
#  $do_pv_blastx:             '1' if we are going to run blastx, else '0'
#  $do_separate_cds_fa_files: '1' if we output a separate file for the protein validation stage
#  $ftr_info_AHR:             REF to array of hashes with feature info 
#  $opt_HHR:                  REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:           REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       
#
################################################################# 
sub make_protein_validation_fasta_file {
  my $sub_name = "make_protein_validation_fasta_file";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($out_fa_file, $mdl_name, $do_pv_blastx, $do_separate_cds_fa_files, $ftr_info_AHR, $opt_HHR, $ofile_info_HHR) = (@_);

  my $ofile_info_key = $mdl_name . ".a.fa";
  my $mdl_fa_file = $ofile_info_HH{"fullpath"}{$ofile_info_key};
  # printf("in $sub_name, ofile_info_key: $ofile_info_key, mdl_fa_file: $mdl_fa_file\n");
  my $nftr = scalar(@{$ftr_info_AHR});

  if($do_pv_blastx) { 
    sqf_FastaFileRemoveDescriptions($mdl_fa_file, $out_fa_file, $ofile_info_HHR);
  }
  else { 
    utl_RunCommand("cp $mdl_fa_file $out_fa_file", opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
  }
  # now add the predicted CDS sequences
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) { 
      my $ofile_info_key = $mdl_name . ".pfa." . $ftr_idx;
      if($do_separate_cds_fa_files) { 
        $ofile_info_key .= ".pv";
      }
      if(exists $ofile_info_HH{"fullpath"}{$ofile_info_key}) { 
        utl_RunCommand("cat " . $ofile_info_HH{"fullpath"}{$ofile_info_key} . " >> $out_fa_file", opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
      }
    }
  }

  return;
}

#################################################################
# Subroutine:  add_protein_validation_alerts
# Incept:      EPN, Wed Mar 18 06:14:21 2020 [generalized to work for hmmer too]
#              EPN, Tue Oct 23 15:54:50 2018
#
# Purpose:    Report protein validation related errors for features of type 'cds'
#             using data stored in earlier parsing of blast or hmmer results 
#             in @{$ftr_results_AAH} (filled in parse_blastx_results(), or
#             parse_hmmer_domtblout).
#
#             Types of alerts added are:
#             "indfantp": adds this alert if blastx/hmmer has a prediction 
#                         for a feature for which there is no CM/nucleotide based prediction
#             "indfantn": adds this alert if protein validation of a CDS prediction fails due to
#                         no blastx hits
#             "indfstrp": adds this alert if protein validation of a CDS prediction fails due to
#                         strand mismatch between CM and blastx prediction
#             "indf5plg": adds this alert if protein validation of a CDS prediction fails due to
#                         protein alignment being too long on 5' end (extending past CM alignment by > 0 nt)
#             "indf5pst": adds this alert if protein validation of a CDS prediction fails due to
#                         protein alignment being too short on 5' end (more than $xalntol shorter than CM)
#             "indf3plg": adds this alert if protein validation of a CDS prediction fails due to
#                         protein alignment being too long on 3' end (extending past CM alignment by > 0 nt)
#             "indf3pst": adds this alert if protein validation of a CDS prediction fails due to
#                         protein alignment being too short on 3' end (more than $xalntol shorter than CM)
#            
#             Alerts only possible added if blastx used (not possible if hmmer used):
#             "insertnp": adds this alert if protein validation of a CDS prediction fails due to
#                         too long of an insert (only possibly reported if BLASTX used)
#             "deletinp": adds this alert if protein validation of a CDS prediction fails due to
#                         too long of a delete (only possibly reported if BLASTX used)
#             "cdsstopp": adds this alert if protein validation of a CDS prediction fails due to
#                         an in-frame stop codon in the protein alignment (only possibly reported
#                         if blastx used
#
# Arguments: 
#  $mdl_name:               name of model we are adding alerts for
#  $seq_name_AR:            REF to array of sequence names, PRE-FILLED
#  $seq_len_HR:             REF to hash of of sequence lengths, PRE-FILLED
#  $ftr_info_AHR:           REF to array of hashes with information on the features, PRE-FILLED
#  $alt_info_HHR:           REF to array of hashes with information on the alerts, PRE-FILLED
#  $ftr_results_HAHR:       REF to feature results HAH, PRE-FILLED
#  $alt_ftr_instances_HHHR: REF to alert instances HAH, ADDED TO HERE
#  $opt_HHR:                REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
################################################################# 
sub add_protein_validation_alerts { 
  my $sub_name = "add_protein_validation_alerts";
  my $nargs_expected = 9;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($mdl_name, $seq_name_AR, $seq_len_HR, $ftr_info_AHR, $alt_info_HHR, $ftr_results_HAHR, $alt_ftr_instances_HHHR, $opt_HHR, $FH_HR) = @_;
  
  my $do_pv_hmmer = opt_Get("--pv_hmmer", $opt_HHR) ? 1 : 0;

  my $nseq = scalar(@{$seq_name_AR});
  my $nftr = scalar(@{$ftr_info_AHR});
  my $seq_idx;   # counter over sequences
  my $seq_name;  # name of one sequence
  my $ftr_idx;   # counter over features
  
  my $aln_tol  = opt_Get("--xalntol",   $opt_HHR); # maximum allowed difference between start/end point prediction between CM and blastx
  my $xmaxins  = opt_Get("--xmaxins",   $opt_HHR); # maximum allowed insertion length in blastx output
  my $xmaxdel  = opt_Get("--xmaxdel",   $opt_HHR); # maximum allowed deletion length in blastx output
  my $minpvlen = opt_Get("--minpvlen", $opt_HHR);

  # determine if we are trimming CDS eventually in the ftr table output, if
  # so we'll change the tolerance on endpoints between nucleotide and protein
  # if the CDS has ambiguities at beginning/end
  my $do_notrim   = opt_Get("--notrim",   $opt_HHR); # 1 to never trim any features
  my %noftrtrim_H = (); # key is feature type read from --noftrtrim <s> option, value is 1 to not trim start/end due to ambiguities
  if(opt_IsUsed("--noftrtrim", $opt_HHR)) { 
    my @noftrtrim_A  = split(",", opt_Get("--noftrtrim", $opt_HHR));
    foreach my $ftr_type (@noftrtrim_A) { $noftrtrim_H{$ftr_type} = 1; }
  }
  my $do_cds_trim = (($do_notrim) || (defined $noftrtrim_H{"CDS"})) ? 0 : 1; 
  
  # get info on position-specific insert and delete maximum exceptions if there are any
  # skip this if we are using hmmer instead of blastx b/c we don't check for inserts/deletes
  # with hmmer
  my @xmaxins_exc_AH = ();
  my @xmaxdel_exc_AH = ();
  if(! $do_pv_hmmer) { 
    for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      %{$xmaxins_exc_AH[$ftr_idx]} = ();
      %{$xmaxdel_exc_AH[$ftr_idx]} = ();

      vdr_FeaturePositionSpecificValueBreakdown($ftr_info_AHR, $ftr_idx, "xmaxins_exc", \%{$xmaxins_exc_AH[$ftr_idx]}, $FH_HR);
      vdr_FeaturePositionSpecificValueBreakdown($ftr_info_AHR, $ftr_idx, "xmaxdel_exc", \%{$xmaxdel_exc_AH[$ftr_idx]}, $FH_HR);
    }
  }

  # for each sequence, for each feature, detect and report alerts
  for($seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    # for each feature
    $seq_name = $seq_name_AR->[$seq_idx];
    if($seq_len_HR->{$seq_name} >= $minpvlen) { 
      for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
        if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) { 
          my $ftr_strand = vdr_FeatureSummaryStrand($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR);
          my $ftr_results_HR = \%{$ftr_results_HAHR->{$seq_name}[$ftr_idx]}; # for convenience
          # printf("in $sub_name, set ftr_results_HR to ftr_results_HAHR->{$seq_name}[$ftr_idx] ");
          my %alt_str_H = ();   # added to as we find alerts below, possible keys are:
          # "indfantp", "indfantn", "indfstrp", "indf5plg", "indf5pst", "indf3plg", "indf3pst", "insertnp", "deletinp", "cdsstopp"
          
          # initialize 
          my $n_start        = undef; # predicted start  from CM 
          my $n_stop         = undef; # predicted stop   from CM 
          my $n_strand       = undef; # predicted strand from CM 
          my $n_len          = undef; # predicted length from CM (summed over all segments)
          my $n_5ablen_pv    = undef; # number of ambiguities at 5' end of CDS
          my $n_3ablen_pv    = undef; # number of ambiguities at 3' end of CDS
          my $n_scoords      = undef; # coordinates of nucleotide prediction in sequence positions
          my $n_mcoords      = undef; # coordinates of nucleotide prediction in model positions
          my $p_qstart       = undef; # predicted query start from blastx
          my $p_qstop        = undef; # predicted query stop  from blastx
          my $p_hstart       = undef; # predicted subject start from blastx
          my $p_hstop        = undef; # predicted subject stop  from blastx
          my $p_sstart       = undef; # blastx predicted sequence start, in 1..seqlen sequence coords
          my $p_sstop        = undef; # blastx predicted sequence stop, in 1..seqlen sequence coords
          my $p_mstart       = undef; # blastx predicted model start, in 1..mdllen model coords
          my $p_mstop        = undef; # blastx predicted model stop,  in 1..mdllen model coords
          my $p_strand       = undef; # predicted strand from blastx
          my $p_ins          = undef; # insert string from blastx
          my $p_del          = undef; # delete string from blastx
          my $p_trcstop      = undef; # premature stop from blastx
          my $p_score        = undef; # raw score from blastx
          my $p_query        = undef; # query name from blastx hit
          my $p_qlen         = undef; # length of query sequence, if $p_feature_flag == 1
          my $p_hlen         = undef; # length of blastx hit
          my $p_qseq_name    = undef; # query seq name parsed out of blast query $p_query
          my $p_qftr_idx     = undef; # feature idx a blast query pertains to, parsed out of blast query $p_query
          my $p_blastx_feature_flag = 0; # set to '1' if $do_pv_hmmer is 0 and $p_query is a fetched feature sequence, not a full length input sequence
          my $p_ftr_scoords  = undef; # set to sequence coords for current features if query is a feature (not full seq)
          
          my $start_diff = undef; # difference in start values between CM and blastx
          my $stop_diff  = undef; # difference in start values between CM and blastx

          my $alt_scoords = undef; # sequence coords for an alert
          my $alt_mcoords = undef; # model coords for an alert
          
          if(defined $ftr_results_HR->{"n_start"}) { 
            $n_start     = $ftr_results_HR->{"n_start"};
            $n_stop      = $ftr_results_HR->{"n_stop"};
            $n_strand    = $ftr_results_HR->{"n_strand"};
            $n_len       = $ftr_results_HR->{"n_len"};
            $n_5ablen_pv = $ftr_results_HR->{"n_5ablen_pv"};
            $n_3ablen_pv = $ftr_results_HR->{"n_3ablen_pv"};
            if(defined $ftr_results_HR->{"n_scoords"}) { $n_scoords = $ftr_results_HR->{"n_scoords"}; }
            if(defined $ftr_results_HR->{"n_mcoords"}) { $n_mcoords = $ftr_results_HR->{"n_mcoords"}; }
          }

          # only proceed if we have a nucleotide prediction >= min length OR
          # we have no nucleotide prediction
          if(((defined $n_start) && ($n_len >= $minpvlen)) || (! defined $n_start)) { 
            if((defined $ftr_results_HR->{"p_qstart"}) && 
               (defined $ftr_results_HR->{"p_qstop"})) { 
              $p_qstart  = $ftr_results_HR->{"p_qstart"};
              $p_qstop   = $ftr_results_HR->{"p_qstop"};
              $p_strand  = $ftr_results_HR->{"p_strand"};
              $p_query   = $ftr_results_HR->{"p_query"};
              $p_hlen    = $ftr_results_HR->{"p_len"};
              if(defined $ftr_results_HR->{"p_ins"})     { $p_ins     = $ftr_results_HR->{"p_ins"};  }
              if(defined $ftr_results_HR->{"p_del"})     { $p_del     = $ftr_results_HR->{"p_del"};  }
              if(defined $ftr_results_HR->{"p_trcstop"}) { $p_trcstop = $ftr_results_HR->{"p_trcstop"}; }
              if(defined $ftr_results_HR->{"p_score"})   { $p_score   = $ftr_results_HR->{"p_score"};   }
              if(defined $ftr_results_HR->{"p_hstart"})  { $p_hstart  = $ftr_results_HR->{"p_hstart"}; }
              if(defined $ftr_results_HR->{"p_hstop"})   { $p_hstop   = $ftr_results_HR->{"p_hstop"};   }

              # determine if the query is a full length sequence, or a fetched sequence feature:
              ($p_qseq_name, $p_qftr_idx, $p_qlen, $p_ftr_scoords) = helper_protein_validation_breakdown_source($p_query, $seq_len_HR, $FH_HR); 
              if($p_qseq_name ne $seq_name) { 
                ofile_FAIL("ERROR, in $sub_name, unexpected query name parsed from $p_query (parsed $p_qseq_name, expected $seq_name)", 1, $FH_HR);
              }
              $p_blastx_feature_flag = ((! $do_pv_hmmer) && ($p_qftr_idx ne "")) ? 1 : 0;
              # printf("seq_name: $seq_name ftr: $ftr_idx p_query: $p_query p_qlen: $p_qlen p_blastx_feature_flag: $p_blastx_feature_flag p_qstart: $p_qstart p_qstop: $p_qstop p_score: $p_score\n");
            }

            # add alerts as needed:
            # check for indfantp
            if((! defined $n_start) && (defined $p_qstart) && (defined $p_score))  { 
              # no nucleotide-based prediction but there is a protein-based blastx prediction
              # only add this if length meets our minimum
              if($p_hlen >= $minpvlen) { 
                if(! $p_blastx_feature_flag) { # this should always be true because if n_start is not defined then there was no 
                                               # nucleotide feature to blastx against, but this is a rare alert so to be safe we require it here
                  $alt_scoords = "seq:" . vdr_CoordsSegmentCreate($p_qstart, $p_qstop, $p_strand, $FH_HR) . ";";
                  $alt_mcoords = "mdl:";
                  if((defined $p_hstart) && (defined $p_hstop)) { 
                    # get subject nucleotide coords
                    # always create in + strand first, vdr_CoordsProteinRelativeToAbsolute requires it
                    my $tmp_alt_mcoords = vdr_CoordsProteinRelativeToAbsolute($ftr_info_AHR->[$ftr_idx]{"coords"},
                                                                              vdr_CoordsSegmentCreate($p_hstart, $p_hstop, "+", $FH_HR), $FH_HR);
                    if($p_strand eq "+") { # just append
                      $alt_mcoords .= $tmp_alt_mcoords . ";";
                    }
                    else { # append rev comp
                      $alt_mcoords .= vdr_CoordsReverseComplement($tmp_alt_mcoords, 0, $FH_HR) . ";"; # 0: don't do carrots
                    }
                  }
                  else { 
                    $alt_mcoords .= "VADRNULL;";
                  }
                }
                else { # $p_blastx_feature_flag is true
                  $alt_scoords = "seq:VADRNULL;";
                  $alt_mcoords = "mdl:VADRNULL;";
                }
                # so we can just use $p_hstart/$p_hstop for model positions
                $alt_str_H{"indfantp"} = $alt_scoords . $alt_mcoords . "raw_score:$p_score";
              }
            }
            if(defined $n_start) { 
              my $cur_5aln_tol = $aln_tol;
              my $cur_3aln_tol = $aln_tol;

              if($do_cds_trim) { 
                # adjust the tolerance to allow the ambiguities at the ends to be missed by the protein validation step
                $cur_5aln_tol += $n_5ablen_pv;
                $cur_3aln_tol += $n_3ablen_pv;
                # if the tolerance is within 3 of the full length, reset it to the default
                if($cur_5aln_tol >= ($n_len - 3)) { 
                  $cur_5aln_tol = $aln_tol;
                }
                if($cur_3aln_tol >= ($n_len - 3)) { 
                  $cur_3aln_tol = $aln_tol;
                }
              }

              # check for indfantn
              if(! defined $p_qstart) { 
                $alt_scoords = (defined $n_scoords) ? 
                    "seq:" . $n_scoords . ";" : "seq:VADRNULL;";
                $alt_mcoords = (defined $n_mcoords) ? 
                    "mdl:" . $n_mcoords . ";" : "mdl:VADRNULL;";
                $alt_str_H{"indfantn"} = $alt_scoords . $alt_mcoords . "VADRNULL";
              }
              else { 
                # we have both $n_start and $p_qstart, we can compare CM and blastx predictions

                # check for indfstrp: strand mismatch failure
                if($n_strand ne $p_strand) { 
                  # first calculate model coords, this is calc'ed same way regardless of value of $p_blastx_feature_flag
                  $alt_mcoords = "mdl:";
                  if((defined $p_hstart) && (defined $p_hstop)) { 
                    # always create in + strand first, vdr_CoordsProteinRelativeToAbsolute requires it
                    my $tmp_alt_mcoords = vdr_CoordsProteinRelativeToAbsolute($ftr_info_AHR->[$ftr_idx]{"coords"}, vdr_CoordsSegmentCreate($p_hstart, $p_hstop, "+", $FH_HR), $FH_HR);
                    if($p_strand eq "+") { # just append
                      $alt_mcoords .= $tmp_alt_mcoords . ";";
                    }
                    else { # append rev comp
                      $alt_mcoords .= vdr_CoordsReverseComplement($tmp_alt_mcoords, 0, $FH_HR) . ";"; # 0: don't do carrots
                    }
                  }
                  else { 
                    $alt_mcoords .= "VADRNULL;";
                  }
                  # calculate seq coords, this is calc'ed differently depending on value of $p_blastx_feature_flag
                  $alt_scoords = "seq:";
                  if($p_blastx_feature_flag) { 
                    if(defined $n_scoords) { 
                      $alt_scoords .= vdr_CoordsRelativeToAbsolute($n_scoords, vdr_CoordsSegmentCreate($p_qstart, $p_qstop, $p_strand, $FH_HR), $FH_HR) . ";";
                    }
                    else { 
                      $alt_scoords .= "VADRNULL;";
                    }
                  }
                  else { # $p_blastx_feature_flag is 0
                    $alt_scoords .= vdr_CoordsSegmentCreate($p_qstart, $p_qstop, $p_strand, $FH_HR) . ";";
                  }
                  $alt_str_H{"indfstrp"} = $alt_scoords . $alt_mcoords . "VADRNULL";
                }
                else { 
                  # we have both $n_start and $p_qstart and predictions on the same strand
                  # determine if predictions are 'close enough' in terms of sequence positions
                  # calcuate $start_diff and $stop_diff, differently depending on if hit
                  # was to the full sequence or a fetched feature (true if $p_blastx_feature_flag == 1)
                  if($p_blastx_feature_flag) { 
                    $start_diff = $p_qstart - 1; 
                    $stop_diff  = $p_qlen - $p_qstop;
                    $p_sstart = ($n_strand eq "+") ? $n_start + $start_diff : $n_start - $start_diff;
                    $p_sstop  = ($n_strand eq "+") ? $n_stop  - $stop_diff  : $n_stop  + $stop_diff;
                    # printf("p_blastx_feature_flag: $p_blastx_feature_flag, p_qstart: $p_qstart, n_start: $n_start p_qstop: $p_qstop, n_stop: $n_stop, start_diff: $start_diff, stop_diff: $stop_diff\n");
                 }
                  else { 
                    $start_diff = abs($n_start - $p_qstart);
                    $stop_diff  = abs($n_stop  - $p_qstop);
                    $p_sstart = $p_qstart;
                    $p_sstop  = $p_qstop;
                  }
                  # check for 'indf5plg': only for non-feature seqs, blastx alignment extends outside of nucleotide/CM alignment on 5' end
                  if((! $p_blastx_feature_flag) && 
                     ((($n_strand eq "+") && ($p_sstart < $n_start)) || 
                      (($n_strand eq "-") && ($p_sstart > $n_start)))) { 
                    $alt_str_H{"indf5plg"} = sprintf("%s%sVADRNULL", 
                                                     "seq:" . vdr_CoordsSegmentCreate($p_sstart, (($n_strand eq "+") ? $n_start-1 : $n_start+1), $n_strand, $FH_HR) . ";", 
                                                     "mdl:" . vdr_CoordsSinglePositionSegmentCreate(vdr_Feature5pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR), $n_strand, $FH_HR) . ";");
                  }
                  # check for 'indf5pst': blastx 5' end too short, not within $cur_5aln_tol nucleotides
                  if(! exists $alt_str_H{"indf5plg"}) { # only add indf5pst if indf5plg does not exist
                    if($start_diff > $cur_5aln_tol) { 
                      $alt_str_H{"indf5pst"} = sprintf("%s%s%d>%d", 
                                                       "seq:" . vdr_CoordsSegmentCreate($n_start, (($n_strand eq "+") ? $p_sstart-1 : $p_sstart+1), $n_strand, $FH_HR) . ";", 
                                                       "mdl:" . vdr_CoordsSinglePositionSegmentCreate(vdr_Feature5pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR), $n_strand, $FH_HR) . ";",
                                                       $start_diff, $cur_5aln_tol);
                    }                
                  }
                  # check for 'indf3plg': blastx alignment extends outside of nucleotide/CM alignment on 3' end
                  if((! $p_blastx_feature_flag) && 
                     ((($n_strand eq "+") && ($p_qstop  > $n_stop)) || 
                      (($n_strand eq "-") && ($p_qstop  < $n_stop)))) { 
                    $alt_str_H{"indf3plg"} = sprintf("%s%sVADRNULL", 
                                                     "seq:" . vdr_CoordsSegmentCreate((($n_strand eq "+") ? $n_stop+1 : $n_stop-1), $p_sstop, $n_strand, $FH_HR) . ";", 
                                                     "mdl:" . vdr_CoordsSinglePositionSegmentCreate(vdr_Feature3pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR), $n_strand, $FH_HR) . ";");
                  }
                  # check for 'indf3pst': blastx 3' end too short, not within $cur_3aln_tol nucleotides
                  # for the stop coordinates, we do this differently if the nucleotide prediction 
                  # includes the stop codon or not, if it does, we allow 3 more positions different
                  my $cur_aln_tol = undef;
                  my $cur_stop_str = undef;
                  my $n_has_stop_codon = 1;
                  if(($ftr_results_HR->{"n_3trunc"}) || 
                     (defined (alert_feature_instance_fetch($alt_ftr_instances_HHHR, $seq_name, $ftr_idx, "mutendcd"))) ||
                     (defined (alert_feature_instance_fetch($alt_ftr_instances_HHHR, $seq_name, $ftr_idx, "ambgnt3c"))) || 
                     (defined (alert_feature_instance_fetch($alt_ftr_instances_HHHR, $seq_name, $ftr_idx, "ambgcd3c")))) { 
                    $n_has_stop_codon = 0; 
                  }                    
                  if($n_has_stop_codon) { 
                    $cur_aln_tol  = $cur_3aln_tol + 3;
                    $cur_stop_str = "valid stop codon";
                  }
                  else { 
                    $cur_aln_tol  = $cur_3aln_tol;
                    $cur_stop_str = "no valid stop codon";
                  }
                  if(! exists $alt_str_H{"indf3plg"}) { # only add indf3pst if indf3plg does not exist
                    if($stop_diff > $cur_aln_tol) { 
                      $alt_str_H{"indf3pst"} = sprintf("%s%s%d>%d", 
                                                       "seq:" . vdr_CoordsSegmentCreate((($n_strand eq "+") ? $p_sstop+1 : $p_sstop-1), $n_stop, $n_strand, $FH_HR) . ";", 
                                                       "mdl:" . vdr_CoordsSinglePositionSegmentCreate(vdr_Feature3pMostPosition($ftr_info_AHR->[$ftr_idx]{"coords"}, $FH_HR), $n_strand, $FH_HR) . ";",
                                                       $stop_diff, $cur_aln_tol);
                      if(! defined (alert_feature_instance_fetch($alt_ftr_instances_HHHR, $seq_name, $ftr_idx, "unexleng"))) { 
                        $alt_str_H{"indf3pst"} .= ", $cur_stop_str in nucleotide-based prediction";
                      }
                    }
                  }
                  # check for 'insertnp': too long of an insert
                  if(defined $p_ins) { 
                    my @p_ins_qpos_A = ();
                    my @p_ins_spos_A = ();
                    my @p_ins_len_A  = ();
                    my $nins = helper_blastx_breakdown_max_indel_str($p_ins, \@p_ins_qpos_A, \@p_ins_spos_A, \@p_ins_len_A, $FH_HR);
                    for(my $ins_idx = 0; $ins_idx < $nins; $ins_idx++) { 
                      my $local_xmaxins = defined ($xmaxins_exc_AH[$ftr_idx]{$p_ins_spos_A[$ins_idx]}) ? $xmaxins_exc_AH[$ftr_idx]{$p_ins_spos_A[$ins_idx]} : $xmaxins;
                      if($p_ins_len_A[$ins_idx] > $local_xmaxins) { 
                        if(defined $alt_str_H{"insertnp"}) { $alt_str_H{"insertnp"} .= ":VADRSEP:"; } # we are adding another instance
                        else                               { $alt_str_H{"insertnp"}  = ""; } # initialize
                        ($alt_scoords, $alt_mcoords) = helper_blastx_max_indel_token_to_alt_coords(1, # $is_insert
                                                                                                   $p_ins_qpos_A[$ins_idx], $p_ins_spos_A[$ins_idx], $p_ins_len_A[$ins_idx], 
                                                                                                   $p_blastx_feature_flag, $p_ftr_scoords, $ftr_info_AHR->[$ftr_idx]{"coords"}, $ftr_strand, $p_strand, 
                                                                                                   $seq_len_HR->{$seq_name}, $FH_HR);
                        $alt_str_H{"insertnp"} .= sprintf("%s%s%d>%d", 
                                                          $alt_scoords, $alt_mcoords, $p_ins_len_A[$ins_idx], $local_xmaxins);
                        
                      }
                    }
                  }
                  # check for 'deletinp': too long of a deletion
                  if(defined $p_del) { 
                    my @p_del_qpos_A = ();
                    my @p_del_spos_A = ();
                    my @p_del_len_A  = ();
                    my $ndel = helper_blastx_breakdown_max_indel_str($p_del, \@p_del_qpos_A, \@p_del_spos_A, \@p_del_len_A, $FH_HR);
                    for(my $del_idx = 0; $del_idx < $ndel; $del_idx++) { 
                      my $local_xmaxdel = defined ($xmaxdel_exc_AH[$ftr_idx]{$p_del_spos_A[$del_idx]}) ? $xmaxdel_exc_AH[$ftr_idx]{$p_del_spos_A[$del_idx]} : $xmaxdel;
                      if($p_del_len_A[$del_idx] > $local_xmaxdel) { 
                        if(defined $alt_str_H{"deletinp"}) { $alt_str_H{"deletinp"} .= ":VADRSEP:"; } # we are adding another instance
                        else                               { $alt_str_H{"deletinp"} = ""; }           # initialize
                        ($alt_scoords, $alt_mcoords) = helper_blastx_max_indel_token_to_alt_coords(0, # $is_insert
                                                                                                   $p_del_qpos_A[$del_idx], $p_del_spos_A[$del_idx], $p_del_len_A[$del_idx], 
                                                                                                   $p_blastx_feature_flag, $p_ftr_scoords, $ftr_info_AHR->[$ftr_idx]{"coords"}, $ftr_strand, $p_strand, 
                                                                                                   $seq_len_HR->{$seq_name}, $FH_HR);
                        $alt_str_H{"deletinp"} .= sprintf("%s%s%d>%d", 
                                                          $alt_scoords, $alt_mcoords, $p_del_len_A[$del_idx], $local_xmaxdel);
                      }
                    }
                  }
                  # check for 'cdsstopp': blast predicted truncation
                  if(defined $p_trcstop) { 
                    my $cdsstopp_alt_str = $p_trcstop; # $p_trcstop is sequence, model coords and detail ("seq:<coords>;mdl:<coords>;<blastx_target>")
                    # laboriously change detail to '-' if it matches model name, idea being including model name in these cases
                    # is redundant and could needlessly confuse user, but for models built from multiple sequences (e.g. cox1)
                    # we want to keep it in there because model coords pertain to the specific blastx query
                    my ($cdsstopp_scoords, $cdsstopp_mcoords, $cdsstopp_detail) = alert_instance_parse($cdsstopp_alt_str);
                    my $beginning_of_detail = ($cdsstopp_detail eq "VADRNULL") ? "VADRNULL" : substr($cdsstopp_detail, 0, length($mdl_name));
                    if($beginning_of_detail eq $mdl_name) { 
                      $cdsstopp_detail = "VADRNULL"; 
                    }
                    elsif($beginning_of_detail ne "VADRNULL") { 
                      $cdsstopp_detail = "mdl-coords_wrt:" . $cdsstopp_detail;
                    }
                    $alt_str_H{"cdsstopp"} = "seq:" . $cdsstopp_scoords . ";mdl:" . $cdsstopp_mcoords . ";" . $cdsstopp_detail;
                  }
                }
              }
            } # end of 'if(defined $n_start)'

            # actually add the alerts
            foreach my $alt_code (sort keys %alt_str_H) { 
              my @alt_str_A = split(":VADRSEP:", $alt_str_H{$alt_code});
              foreach my $alt_str (@alt_str_A) { 
                alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, $alt_code, $seq_name, $ftr_idx, $alt_str, $FH_HR);
              }
            }
          } # end of 'if(((defined $n_start) && ($n_len >= $minpvlen)) || (! defined $n_start))'
        } # end of 'if(featureTypeIsCds($ftr_info_AHR, $ftr_idx))'
      } # end of 'for($ftr_idx' loop
    } # end of 'if($seq_len_HR->{$seq_name} >= $minpvlen)'
  } # end of 'for($seq_idx' loop
  return;
}
#################################################################
# Subroutine:  run_blastx_and_summarize_output()
# Incept:      EPN, Thu Oct  4 15:25:00 2018
#
# Purpose:    For each fasta file of predicted hits, run them as
#             blastx queries against the appropriate target blast DBs
#             and summarize the output with parse_blast.pl.
#
# Arguments: 
#  $execs_HR:                  REF to a hash with "blastx" and "parse_blastx.pl""
#                              executable paths
#  $out_root:                  output root for the file names
#  $mdl_info_HR:               REF to hash of arrays of model info
#  $ftr_info_AHR:              REF to hash of arrays with information on the features, PRE-FILLED
#  $blastx_db_file             blastx DB file to use (usually $mdl_info_HR->{"blastdbpath"}, but can be diff if --xsub used)
#  $do_separate_cds_fa_files:  '1' if we output a separate file for the protein validation stage
#  $opt_HHR:                   REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:            REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If blastx fails.
#
################################################################# 
sub run_blastx_and_summarize_output {
  my $sub_name = "run_blastx_and_summarize_output";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $out_root, $mdl_info_HR, $ftr_info_AHR, $blastx_db_file,$do_separate_cds_fa_files, $opt_HHR, $ofile_info_HHR) = @_;

  my $do_keep = opt_Get("--keep", $opt_HHR);
  my $nftr = scalar(@{$ftr_info_AHR});
  my $mdl_name = $mdl_info_HR->{"name"};
  my $ncpu = opt_Get("--cpu", $opt_HHR);
  if($ncpu == 0) { $ncpu = 1; }

  # make a query fasta file for blastx, consisting of full length
  # sequences (with sequence descriptions removed because they can
  # affect the output and mess up our parsing if they are too long)
  # AND all the predicted CDS sequences
  my $blastx_query_fa_file = $out_root . "." . $mdl_name . ".pv.blastx.fa";
  make_protein_validation_fasta_file($blastx_query_fa_file, $mdl_name,  1, $do_separate_cds_fa_files, $ftr_info_AHR, $opt_HHR, $ofile_info_HHR);
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".pv-blastx-fasta", $blastx_query_fa_file, 0, opt_Get("--keep", \%opt_HH), "sequences for protein validation for model $mdl_name");
  
  # run blastx 
  my $blastx_options = "";
  if(defined $mdl_info_HR->{"transl_table"}) { 
    $blastx_options .= " -query_gencode " . $mdl_info_HR->{"transl_table"};
  }
  if(opt_IsUsed("--xmatrix", $opt_HHR)) { 
    $blastx_options .= " -matrix " . opt_Get("--xmatrix", $opt_HHR); 
  }
  if(opt_IsUsed("--xdrop", $opt_HHR)) { 
    my $xdrop_opt = opt_Get("--xdrop", $opt_HHR);
    $blastx_options .= " -xdrop_ungap $xdrop_opt -xdrop_gap $xdrop_opt -xdrop_gap_final $xdrop_opt";
  }
  my $xnumali = opt_Get("--xnumali", $opt_HHR);

  my $blastx_out_file = $out_root . "." . $mdl_name . ".blastx.out";
  my $blastx_cmd = $execs_HR->{"blastx"} . " -num_threads $ncpu -num_alignments $xnumali -query $blastx_query_fa_file -db $blastx_db_file -seg no -out $blastx_out_file" . $blastx_options;
  utl_RunCommand($blastx_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".blastx-out", $blastx_out_file, 0, $do_keep, "blastx output for model $mdl_name");

  # now summarize its output
  my $blastx_summary_file = $out_root . "." . $mdl_name . ".blastx.summary.txt";
  my $parse_cmd = $execs_HR->{"parse_blast"} . " --program x --input $blastx_out_file > $blastx_summary_file";
  utl_RunCommand($parse_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".blastx-summary", $blastx_summary_file, 0, $do_keep, "parsed (summarized) blastx output");

  return;
}

#################################################################
# Subroutine:  parse_blastx_results()
# Incept:      EPN, Thu Oct  4 15:25:00 2018
#              [modified from subroutine of same name by Alejandro Schaffer, compare_predictions.pl]
#
# Purpose:    Parse blastx summary file and store results
#             in ftr_results.
#
# Arguments: 
#  $blastx_summary_file: path to blastx summary file to parse
#  $seq_name_AR:         REF to array of sequence names
#  $seq_len_HR:          REF to hash of sequence lengths
#  $ftr_info_AHR:        REF to array of hashes with feature info 
#  $ftr_results_HAHR:    REF to feature results HAH, ADDED TO HERE
#  $opt_HHR:             REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:      REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If blastx fails.
#
################################################################# 
sub parse_blastx_results { 
  my $sub_name = "parse_blastx_results";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastx_summary_file, $seq_name_AR, $seq_len_HR, $ftr_info_AHR, $ftr_results_HAHR, $opt_HHR, $ofile_info_HHR) = @_;

  # printf("in $sub_name, blastx_summary_file: $blastx_summary_file\n");
  
  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  utl_FileValidateExistsAndNonEmpty($blastx_summary_file, "blastx summary file", $sub_name, 1, $FH_HR);
  my $nftr = scalar(@{$ftr_info_AHR});
  my $do_xlongest = opt_Get("--xlongest", $opt_HHR) ? 1 : 0;

  # create a hash mapping ftr_type_idx strings to ftr_idx:
  my %ftr_type_idx2ftr_idx_H = ();
  vdr_FeatureInfoMapFtrTypeIndicesToFtrIndices($ftr_info_AHR, \%ftr_type_idx2ftr_idx_H, $FH_HR);

  # check for special case of there only being 1 CDS feature
  # if so, we don't need to follow the idiom that the blast 
  # target name is <protein-accession>/<coords-str>
  # (that is, when the blast db was created it did not
  # need to have its sequences named this way)
  my $ncds = vdr_FeatureInfoCountType($ftr_info_AHR, "CDS");
  if($ncds == 0) { 
    ofile_FAIL("ERROR in $sub_name, unable to find a CDS in model info, no need for blastx steps...", 1, $FH_HR);
  }
  my $ftr_idx_lone_cds    = ($ncds == 1) ? $ftr_type_idx2ftr_idx_H{"CDS.1"} : -1;
  my $ftr_strand_lone_cds = ($ncds == 1) ? vdr_FeatureSummaryStrand($ftr_info_AHR->[$ftr_idx_lone_cds]{"coords"}, $FH_HR) : undef;

  open(IN, $blastx_summary_file) || ofile_FileOpenFailure($blastx_summary_file, $sub_name, $!, "reading", $FH_HR);
  
  my $line_idx   = 0;
  my $minpvlen   = opt_Get("--minpvlen",  $opt_HHR);
  my $xlonescore = opt_Get("--xlonescore", $opt_HHR);
  my $do_xnoid   = opt_Get("--xnoid", $opt_HHR);
  my $seq_name   = undef; # sequence name this hit corresponds to 
  my $q_len      = undef; # length of query sequence
  my $q_coords   = undef; # coordinates of of query sequence, parsed from QACC line
  my $q_ftr_idx  = undef; # feature index query pertains to [0..$nftr-1] OR -1: a special case meaning query is full sequence (not a fetched CDS feature)
  my $t_ftr_idx  = undef; # feature index target (subject) pertains to [0..$nftr-1]
  my $t_strand   = undef; # summary strand of feature $t_ftr_idx
  my %cur_H = (); # values for current hit (HSP)
  
  # 
  # Order of lines in <IN>:
  # -----per-query/target-block---
  # QACC
  # QDEF   ignored
  # MATCH  ignored
  # HACC  
  # HDEF   ignored
  # SLEN   ignored
  # ------per-HSP-block------
  # HSP   
  # RAWSCORE 
  # EVALUE ignored
  # HLEN   ignored
  # IDENT  ignored
  # GAPS   ignored
  # FRAME
  # QSTOP  not always present
  # HSTOP  not always present
  # DEL    not always present
  # MAXDE  not always present
  # INS    not always present
  # MAXINS not always present
  # QRANGE
  # SRANGE
  # ------per-HSP-block------
  # END_MATCH

  while(my $line = <IN>) { 
    chomp $line;
    $line_idx++;
    if($line ne "END_MATCH") { 
      my @el_A = split(/\t/, $line);
      if(scalar(@el_A) != 2) { 
        ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, did not read exactly 2 tab-delimited tokens in line $line", 1, $FH_HR);
      }
      my ($key, $value) = (@el_A);
      if($key eq "QACC") { 
        $cur_H{$key} = $value;
        # determine what sequence it is
        my $q_ftr_type_idx; # feature type and index string, from $seq_name if not a full sequence (e.g. "CDS.4")
        ($seq_name, $q_ftr_type_idx, $q_len, $q_coords) = helper_protein_validation_breakdown_source($value, $seq_len_HR, $FH_HR); 
        # helper_protein_validation_breakdown_source() will die if $query is unparseable
        # determine what feature this query corresponds to
        $q_ftr_idx = ($q_ftr_type_idx eq "") ? -1 : $ftr_type_idx2ftr_idx_H{$q_ftr_type_idx};
        if(! defined $q_ftr_idx) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, problem parsing QACC line unable to determine feature index from query $value", 1, $FH_HR);
        }
      }
      elsif($key eq "HACC") { 
        if(! defined $cur_H{"QACC"}) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read HACC line before QACC line (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
        # determine what feature it is, unless we only have 1 CDS in which case we assume it's that 1 CDS
        if($ftr_idx_lone_cds != -1) { 
          $t_ftr_idx = $ftr_idx_lone_cds;
          $t_strand  = $ftr_strand_lone_cds;
        }
        elsif($value =~ /(\S+)\/(\S+)/) { 
          my ($accn, $coords) = ($1, $2);
          # find it in @{$ftr_info_AHR} (or set to lone CDS if there is only 1
          ($t_ftr_idx, $t_strand) = helper_protein_validation_db_seqname_to_ftr_idx($value, $ftr_info_AHR, $FH_HR); # will die if problem parsing $target, or can't find $t_ftr_idx
        }
        else {
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse HACC line $line", 1, $FH_HR);
        }
      }
      elsif($key eq "SLEN") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read SLEN line before QACC line (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        if($value !~ /^(\d+)$/) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary SLEN line $line", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
      }
      elsif($key eq "HSP") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read HSP line before one or both of QACC and HACC lines (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
        if($value !~ /^(\d+)$/) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary HSP line $line", 1, $FH_HR);
        }
      }
      elsif($key eq "RAWSCORE") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"}) || (! defined $cur_H{"HSP"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read RAWSCORE line before one or more of QACC, HACC, or HSP lines (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
        if($value !~ /^(\d+)$/) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary RAWSCORE line $line", 1, $FH_HR);
        }
      }
      elsif($key eq "HLEN") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"}) || (! defined $cur_H{"HSP"}) || (! defined $cur_H{"RAWSCORE"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read HLEN line before one or more of QACC, HACC, HSP, or RAWSCORE lines (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        if($value !~ /^(\d+)$/) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary HLEN line $line", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
      }
      elsif($key eq "IDENT") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"}) || (! defined $cur_H{"HSP"}) || (! defined $cur_H{"RAWSCORE"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read IDENT line before one or more of QACC, HACC, HSP, or RAWSCORE lines (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        if($value !~ /^(\d+)$/) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary IDENT line $line", 1, $FH_HR);
        }
        $cur_H{$key} = $value;
      }
      elsif($key eq "FRAME") { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"}) || (! defined $cur_H{"HSP"}) || (! defined $cur_H{"RAWSCORE"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read FRAME line before one or more of QACC, HACC, HSP, or RAWSCORE lines (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        if($value =~ /^[\+\-]([123])$/) { 
          $cur_H{$key} = $1;
        }
        else { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary FRAME line $line ($key $value)", 1, $FH_HR);
        }
      }
      elsif(($key eq "QSTOP") || ($key eq "HSTOP") || ($key eq "DEL") || ($key eq "INS")) { 
        if((! defined $cur_H{"QACC"}) || (! defined $cur_H{"HACC"}) || (! defined $cur_H{"HSP"}) || (! defined $cur_H{"RAWSCORE"}) || (! defined $cur_H{"FRAME"})) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read $key line before one or more of QACC, HACC, HSP, RAWSCORE or FRAME lines (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        if(($value ne "") && ($value ne "BLASTNULL")) { 
          $cur_H{$key} = $value;
        } 
      }
      elsif($key eq "SRANGE") { 
        # we don't require all of QACC, HACC, HSP, RAWSCORE and FRAME even though we should have them
        # sometimes we don't (may be a bug in parse-blastx.pl), we only require QACC
        if(! defined $cur_H{"QACC"}) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read $key line before QACC line (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        if($value eq "..") { # special case, no hits, silently move on
          ;
        }
        elsif(! defined $cur_H{"RAWSCORE"}) { # special case, no RAWSCORE lines yet seen (may be a bug in parse-blastx.pl?), silently move on
          ;
        }
        else { 
          # determine if we should store this hit
          if($value =~ /^\d+..\d+$/) { 
            $cur_H{$key} = $value;
          }
          else { 
            ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary QRANGE line $line", 1, $FH_HR);
          }
        }
      }
      elsif($key eq "QRANGE") { 
        # we don't require all of QACC, HACC, HSP, RAWSCORE and FRAME even though we should have them
        # sometimes we don't (may be a bug in parse-blastx.pl), we only require QACC
        if(! defined $cur_H{"QACC"}) { 
          ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, read $key line before QACC line (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
        }
        if($value eq "..") { # special case, no hits, silently move on
          ;
        }
        elsif(! defined $cur_H{"RAWSCORE"}) { # special case, no RAWSCORE lines yet seen (may be a bug in parse-blastx.pl?), silently move on
          ;
        }
        else { 
          # determine if we should store this hit
          if($value =~ /^(\d+)..(\d+)$/) { 
            my ($blast_qstart, $blast_qstop) = ($1, $2);
            my $blast_strand = ($blast_qstart <= $blast_qstop) ? "+" : "-";
            
            # should we store this query/target/hit trio?
            # we do if A, B, C, D are all TRUE and at least one of X or Y or Z is TRUE
            #  A. this query/target pair is compatible (query is full sequence or correct CDS feature) 
            #  B. if --xlongest not used: this is the highest scoring hit for this feature for this sequence (query/target pair)? 
            #     if --xlongest is  used: this is the longest hit (query coords) for this feature for this sequence (query/target pair)? 
            #  C. query length (full length seq or predicted CDS) is at least <x> nt from --minpvlen
            #  D. --xnoid NOT used OR alignment is not full length and 100% identical (--xnoid is an option used for
            #     testing training sets where we want to ignore self hits)
            # 
            #  X. query is a single predicted CDS feature (not the full sequence), in this case checking for overlap
            #     won't work because blast coords are relative to CDS not full sequence, and we already know it
            #     must overlap because query/target are compatible and query is not full seq (A)
            #  Y. hit score is above minimum (--xlonescore)
            #  Z. query is the full sequence ($q_ftr_idx == -1) and blast hit overlaps by at least 1 nt 
            #     with a nucleotide prediction for the current target CDS

            my $blast_hit_qlen = abs($blast_qstart - $blast_qstop) + 1;
            my $a_true = (($q_ftr_idx == -1) || ($q_ftr_idx == $t_ftr_idx)) ? 1 : 0; # query is full sequence OR query is fetched CDS that pertains to target
            my $b_true = undef;
            if(! $do_xlongest) { 
              $b_true = ((! defined $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_score"}) ||  # first hit, so must be highest score
                         ($cur_H{"RAWSCORE"} > $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_score"})) ? 1 : 0; # highest scoring hit
            }
            else { 
              $b_true = ((! defined $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_score"}) ||  # first hit, so must be longest
                         ($blast_hit_qlen > $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_len"})) ? 1 : 0; # longest hit
            }

            my $c_true = ($q_len >= $minpvlen) ? 1 : 0; # length >= --minpvlen
            my $d_true = $do_xnoid ? 0 : 1;
            if(! $d_true) { 
              if(((defined $cur_H{"SLEN"}) && (defined $cur_H{"HLEN"}) && (defined $cur_H{"IDENT"})) && 
                 ($cur_H{"SLEN"} == $cur_H{"HLEN"}) && ($cur_H{"SLEN"} == $cur_H{"IDENT"})) { 
                $d_true = 0; # 100% identical full length blast hit, --xnoid is enabled so we want to skip this
              }
              else { 
                $d_true = 1; 
              }
            }
            if($a_true && $b_true && $c_true && $d_true) { 
              my $x_true = ($q_ftr_idx != -1) ? 1 : 0; # D is true if query is a single CDS as opposed to full sequence
              my $y_true = ($cur_H{"RAWSCORE"} >= $xlonescore) ? 1 : 0; # E is true if raw score exceeds --xlonescore
              my $z_true = 0; 
              # only bother determining $z_true if $x_true and $y_true are both false (0)
              if((! $x_true) && (! $y_true)) { 
                if((defined $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"n_strand"}) &&
                   ($ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"n_strand"} eq $blast_strand)) { 
                  my $noverlap = helper_protein_validation_check_overlap($ftr_results_HAHR->{$seq_name}[$t_ftr_idx], $blast_qstart, $blast_qstop, $blast_strand, $FH_HR);
                  if($noverlap > 0) { $z_true = 1; }
                }
              }
              if($x_true || $y_true || $z_true) { 
                # store the hit
                $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_qstart"} = $blast_qstart;
                $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_qstop"}  = $blast_qstop;
                $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_strand"} = $blast_strand;
                $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_len"}    = $blast_hit_qlen;
                $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_query"}  = $cur_H{"QACC"};
                $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_score"}  = $cur_H{"RAWSCORE"};
                $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_frame"}  = $cur_H{"FRAME"};
                # p_strand is actually complicated to set:
                # if q_ftr_idx == -1, query is the full sequence in which case we use blast_strand from output
                # if q_ftr_idx != -1, query is a single CDS so blast_strand will be + if same strand as target, - if opposite 
                if($q_ftr_idx == -1) { 
                  $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_strand"} = $blast_strand;
                }
                else { 
                  # query is a single CDS, if blast_strand is "+" this means that it agrees with 
                  # strand of the CDS, if blast_strand is "-" then it disagrees
                  if(! defined $t_strand) { 
                    ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to set strand of hit because t_strand is undef", 1, $FH_HR);
                  }
                  if($blast_strand eq "+") { 
                    $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_strand"} = $t_strand;
                  } 
                  elsif($t_strand eq "+") { # blast_strand is "-"
                    $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_strand"} = "-";
                  }
                  elsif($t_strand eq "-") { # blast_strand is "-"
                    $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_strand"} = "+";
                  }
                  else { # blast_strand is "-", t_strand is not "+" or "-", should be "?"
                    $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_strand"} = "?";
                  }
                }
                if(defined $cur_H{"INS"}) { 
                  $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_ins"} = $cur_H{"INS"};
                }
                else { 
                  $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_ins"} = undef;
                }
                if(defined $cur_H{"DEL"}) { 
                  $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_del"} = $cur_H{"DEL"};
                }
                else { 
                  $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_del"} = undef;
                }
                if((defined $cur_H{"SRANGE"}) && ($cur_H{"SRANGE"} =~ /^(\d+)..(\d+)$/)) { 
                  ($ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_hstart"}, $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_hstop"}) = 
                      ($1, $2);
                }
                else { 
                  ofile_FAIL("ERROR in $sub_name, trying to store hit but don't have defined and parseable SRANGE line (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
                }
                if(defined $cur_H{"QSTOP"}) { 
                  if(! defined $cur_H{"HSTOP"}) { 
                    ofile_FAIL("ERROR in $sub_name, read QSTOP but not HSTOP line (seq: $seq_name, line: $line_idx)\n", 1, $FH_HR);
                  }
                  # we want to store only the first (5'-most) stop codon
                  # query is the sequence, subject is the model (protein sequence)
                  my $first_qstop = $cur_H{"QSTOP"};
                  $first_qstop =~ s/\;.*$//; # remove first ';' and anything after it
                  if($q_ftr_idx != -1) { # query is a feature, not the full sequence, so we need to determine coords in full seq
                    if(defined $q_coords) { # $q_coords should be defined but this is just a sanity check
                      $first_qstop = vdr_CoordsRelativeToAbsolute($q_coords, $first_qstop, $FH_HR);
                    }
                    else { 
                      $first_qstop = "-"
                    }
                  }
                  my $first_hstop = $cur_H{"HSTOP"};
                  $first_hstop =~ s/\;.*$//; # remove first ';' and anything after it
                  $first_hstop = vdr_CoordsProteinRelativeToAbsolute($ftr_info_AHR->[$t_ftr_idx]{"coords"}, $first_hstop, $FH_HR);
                  my $hacc_accn = $cur_H{"HACC"};
                  if($hacc_accn =~ /(\S+)\/\S+/) { 
                    $hacc_accn = $1;
                  }
                  $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_trcstop"} = "seq:" . $first_qstop . ";mdl:" . $first_hstop . ";" . $hacc_accn;
                }
                else { 
                  $ftr_results_HAHR->{$seq_name}[$t_ftr_idx]{"p_trcstop"} = undef;
                }
              }
            }
          }
          else { 
            ofile_FAIL("ERROR in $sub_name, reading $blastx_summary_file, unable to parse blastx summary QRANGE line $line", 1, $FH_HR);
          }
        } # end of 'else' entered if QRANGE is NOT ".."
        # reset variables
        my $save_qacc = $cur_H{"QACC"}; 
        my $save_hacc = $cur_H{"HACC"};
        %cur_H = (); 
        $cur_H{"QACC"} = $save_qacc; # save QACC
        $cur_H{"HACC"} = $save_hacc; # save HACC
      } # end of 'elsif $key eq "QRANGE"
    } # end of 'if($line ne "END_MATCH")'
    else { # $line eq "END_MATCH"
      # reset variables
      my $save_qacc = $cur_H{"QACC"};
      %cur_H = (); 
      $cur_H{"QACC"} = $save_qacc; # save QACC
      $t_ftr_idx = undef;
    }
  } # end of 'while($my $line = <IN>)'
  close(IN);

  return 0;
}

#################################################################
# Subroutine:  run_esl_translate_and_hmmsearch()
# Incept:      EPN, Mon Mar  9 14:51:56 2020
#
# Purpose:    Translate each CDS nt fasta file, and run hmmsearch
#             against the resulting protein sequences.
#
# Arguments: 
#  $execs_HR:                 REF to a hash with "blastx" and "parse_blastx.pl""
#                             executable paths
#  $out_root:                 output root for the file names
#  $mdl_info_HR:              REF to hash of model info
#  $ftr_info_AHR:             REF to hash of arrays with information on the features, PRE-FILLED
#  $do_separate_cds_fa_files: '1' if we output a separate file for the protein validation stage
#  $opt_HHR:                  REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:           REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If esl-translate fails.
#
################################################################# 
sub run_esl_translate_and_hmmsearch { 
  my $sub_name = "run_esl_translate_and_hmmsearch";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $out_root, $mdl_info_HR, $ftr_info_AHR, $do_separate_cds_fa_files, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $do_keep = opt_Get("--keep", $opt_HHR);
  my $nftr = scalar(@{$ftr_info_AHR});
  my $mdl_name = $mdl_info_HR->{"name"};

  my $model_domtblout_file = $out_root . "." . $mdl_name . ".hmmscan.domtblout";
  # make a query fasta file for blastx, consisting of full length
  # sequences (with sequence descriptions removed because they can
  # affect the output and mess up our parsing if they are too long)
  # AND all the predicted CDS sequences
  my $pv_fa_file = $out_root . "." . $mdl_name . ".pv.hmmer.fa";
  make_protein_validation_fasta_file($pv_fa_file, $mdl_name,  0, $do_separate_cds_fa_files, $ftr_info_AHR, $opt_HHR, $ofile_info_HHR); # 0: not doing blastx
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".pv.hmmer.fa", $pv_fa_file, 0, opt_Get("--keep", \%opt_HH), "sequences for protein validation for model $mdl_name");

  # now esl-translate it
  my $esl_translate_opts = "-l 1 -W "; # -W avoids presumed bug in esl-translate for input seqs of length < 3 (see 20_0420_vadr_1p1_models/00LOG.txt)
  if(defined $mdl_info_HR->{"transl_table"}) { 
    $esl_translate_opts .= " -c " . $mdl_info_HR->{"transl_table"};
  }
  my $esl_translate_prot_fa_file = $out_root . "." . $mdl_name . ".pv.hmmer.esl_translate.aa.fa";
  my $esl_translate_cmd = $execs_HR->{"esl-translate"} . " $esl_translate_opts $pv_fa_file > $esl_translate_prot_fa_file";
  utl_RunCommand($esl_translate_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".pv.hmmer.esl_translate.aa.fa", $esl_translate_prot_fa_file, 0, $do_keep, "esl-translate output for protein validation sequences for model $mdl_name");

  # run hmmsearch against it using only those HMMs that pertain to this model
  my @hmm_name_A = (); 
  get_hmm_list_for_model($ftr_info_AHR, $mdl_name, \@hmm_name_A, $FH_HR);
  my $hmm_list_file = $out_root . "." . $mdl_name . ".pv.hmmlist";
  utl_AToFile(\@hmm_name_A, $hmm_list_file, 1, $FH_HR);

  my $hmmsearch_opts = "";
  my $hmmsearch_out_file       = $out_root . "." . $mdl_name . ".hmmsearch.out";
  my $hmmsearch_domtblout_file = $out_root . "." . $mdl_name . ".hmmsearch.domtblout";
  #my $hmmsearch_stk_file       = $out_root . "." . $mdl_name . ".hmmsearch.stk";
  #$hmmsearch_opts .= " -A $hmmsearch_stk_file";
  if(opt_Get("--h_max", $opt_HHR)) { $hmmsearch_opts .= " --max"; }
  my $hmmfetch_cmd  = $execs_HR->{"hmmfetch"}  . " -f $hmm_pt_file $hmm_list_file | ";
  my $hmmsearch_minbit = opt_Get("--h_minbit", $opt_HHR);
  my $hmmsearch_cmd = $hmmfetch_cmd . " " . $execs_HR->{"hmmsearch"} . " --domT $hmmsearch_minbit -T $hmmsearch_minbit --domtblout $hmmsearch_domtblout_file $hmmsearch_opts - $esl_translate_prot_fa_file > $hmmsearch_out_file";
  utl_RunCommand($hmmsearch_cmd, opt_Get("-v", $opt_HHR), 0, $ofile_info_HHR->{"FH"});
        
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".hmmlist",        $hmm_list_file,            0, $do_keep, "list of hmm files fetched for model $mdl_name");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".hmmsearch",      $hmmsearch_out_file,       0, $do_keep, "hmmsearch standard output for model $mdl_name");
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".domtblout",      $hmmsearch_domtblout_file, 0, $do_keep, "hmmsearch --domtblout output for model $mdl_name");

  return;
}

#################################################################
# Subroutine:  parse_hmmer_domtblout()
# Incept:      EPN, Tue Mar 10 06:23:54 2020
#
# Purpose:    Parse hmmsearch or hmmscan --domtblout file and store
#             the results in ftr_results.
#
# Arguments: 
#  $domtblout_file:      hmmsearch or hmmscan  --domtblout file to parse
#  $do_hmmscan:          "1" if output is from hmmscan, "0" if from hmmsearch
#  $seq_name_AR:         REF to array of sequence names
#  $seq_len_HR:          REF to hash of sequence lengths
#  $ftr_info_AHR:        REF to array of hashes with feature info 
#  $ftr_results_HAHR:    REF to feature results HAH, ADDED TO HERE
#  $opt_HHR:             REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:      REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If there's a problem parsing the domtblout out file 
#             because the format is unexpected.
#
################################################################# 
sub parse_hmmer_domtblout {
  my $sub_name = "parse_hmmer_domtblout";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($domtblout_file, $do_hmmscan, $seq_name_AR, $seq_len_HR, $ftr_info_AHR, $ftr_results_HAHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  utl_FileValidateExistsAndNonEmpty($domtblout_file, "hmmer --domtblout file", $sub_name, 1, $FH_HR);
  my $nftr = scalar(@{$ftr_info_AHR});

  open(IN, $domtblout_file) || ofile_FileOpenFailure($domtblout_file, $sub_name, $!, "reading", $FH_HR);
  
  my $line_idx   = 0;
  my $hminntlen  = opt_Get("--minpvlen",  $opt_HHR); # yes, it should be hminntlen = --minpvlen
  my $hlonescore = opt_Get("--hlonescore", $opt_HHR);
  my $seq_name   = undef; # sequence name this hit corresponds to 
  my $q_len      = undef; # length of query sequence
  my $q_ftr_idx  = undef; # feature index query pertains to, [0..$nftr-1]

  # create a hash mapping ftr_type_idx strings to ftr_idx:
  my %ftr_type_idx2ftr_idx_H = ();
  vdr_FeatureInfoMapFtrTypeIndicesToFtrIndices($ftr_info_AHR, \%ftr_type_idx2ftr_idx_H, $FH_HR);

##                                                                                --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
## target name        accession   tlen query name               accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
##------------------- ---------- -----     -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
#orf31                -            208 cox1.tt5.arthropoda.cds1 -            511  5.4e-130  426.4  16.0   1   1  9.6e-132  6.6e-130  426.1  16.0    35   239     1   205     1   207 1.00 source=lcl|Seq1674/CDS.1/1..627:+ coords=3..626 length=208 frame=3 desc=

# white-space separated fields:
#  0: target name : orf<d>         : irrelevant
#  1: accession   : -              : irrelevant
#  2: tlen        : <d>            : length of translated protein
#  3: query name  : <model>.cds<n> : <n> tells us which CDS this pertains to
#  4: accession   : -              : irrelevant
#  5: qlen        : <d>            : aa length of predicted CDS
#  6: E-value     : <g>            : irrelevant, E-value for full sequence
#  7: score       : <f>            : irrelevant, bit score for full sequence
#  8: bias        : <f>            : irrelevant, bias score for full sequence
#  9: #           : <n>            : irrelevant, index of domain/hit in sequence
# 10: of          : <n>            : irrelevant, total number of domains/hits in sequence
# 11: c-Evalue    : <g>            : conditional E-value for this domain/hit
# 12: i-Evalue    : <g>            : independent E-value for this domain/hit
# 13: score       : <f>            : bit score for this domain/hit
# 14: bias        : <f>            : bias score for full sequence
# 15: hmm from    : <n>            : model start position of HMM alignment (not envelope)
# 16: hmm to      : <n>            : model end position of HMM alignment (not envelope)
# 17: ali from    : <n>            : protein sequence start position of HMM alignment (not envelope)
# 18: ali to      : <n>            : protein sequence end position of HMM alignment (not envelope)
# 19: env from    : <n>            : protein sequence start position of HMM envelope
# 20: env to      : <n>            : protein sequence end position of HMM envelope
# 21: acc         : <f>            : irrelvant, average pp of aligned protein
# description can be broken down into parseable tokens ONLY because we
# know how esl-translate creates this description, at least up to the desc= field
# 22: source=     : source=<s>     : <s> is <s1>\/CDS.<d>\/<s2> where <s1> is seq name, 
#                                  : <d> is CDS idx and <s2> is coords string of predicted CDS coords
# 23: coords=     : coords=<s>     : <d>..<d> coordinates of translated protein in nucleotide space
#                                  : from esl-translate
# 24: length=     : length=<d>     : <d> is length of translated protein, redundant with field 2, tlen
# 25: frame=      : frame=<d>      : frame of translated protein from esl-translate
# 26: desc=       : desc=<s>       : irrelevant, <s> is first field of description from nt seq input to esl-translate
# 27+: may or may not exist, irrelevant anyway
#      will exist only if a description existed for the nt seq
#      input to esl-translate, that is - if we added descriptions 
#      the seq in the esl-translate input file
#  
  my $seq_ftr_idx = undef; # feature index sequence being evaluated pertains to, [0..$nftr-1] 
                           # OR -1: a special case meaning sequence is full sequence (not a fetched CDS feature)
  my $mdl_ftr_idx = undef; # feature index the HMM pertains to [0..$nftr-1]
  while(my $line = <IN>) { 
    chomp $line;
    $line_idx++;
    if($line !~ m/^#/) { 
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) < 27) { 
        ofile_FAIL("ERROR in $sub_name, reading $domtblout_file, did not read at least 25 space-delimited tokens in data line\n$line", 1, $FH_HR);
      }
      my ($orig_seq_name, $orig_seq_len, $mdl_name, $mdl_len, $hit_ieval, $hit_score, $hit_bias, $hmm_from, $hmm_to, $ali_from, $ali_to, $env_from, $env_to, 
          $source_str, $coords_str, $qlen_str, $frame_str) = 
              ($el_A[0], $el_A[2], $el_A[3], $el_A[5], $el_A[12], $el_A[13], $el_A[14], $el_A[15], $el_A[16], $el_A[17], $el_A[18], $el_A[19], $el_A[20],
               $el_A[22], $el_A[23],  $el_A[24], $el_A[25]);
      if($do_hmmscan) { # swap query/target names and lengths
        utl_Swap(\$orig_seq_name, \$mdl_name);
        utl_Swap(\$orig_seq_len,  \$mdl_len);
      }
      ($mdl_ftr_idx, undef) = helper_protein_validation_db_seqname_to_ftr_idx($mdl_name, $ftr_info_AHR, $FH_HR); # will die if problem parsing $target, or can't find $t_ftr_idx

      # further parse some of the tokens
      my ($seq_ftr_type_idx, $seq_len, $source_coords, $source_val);
      if($source_str =~ /^source=(\S+)$/) { 
        $source_val = $1;
        ($seq_name, $seq_ftr_type_idx, $seq_len, $source_coords) = helper_protein_validation_breakdown_source($source_val, $seq_len_HR, $FH_HR); 
        if($seq_ftr_type_idx eq "") {
          $seq_ftr_idx = -1; # hit is to full sequence
        }
        else { 
          $seq_ftr_idx = $ftr_type_idx2ftr_idx_H{$seq_ftr_type_idx};
        }
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, reading $domtblout_file, unable to parse source= token ($source_str) on line\n$line", 1, $FH_HR);
      }      

      my ($orf_start, $orf_end) = (undef, undef);
      if($coords_str =~ /^coords=(\d+)\.\.(\d+)$/) { 
         ($orf_start, $orf_end) = ($1, $2);
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, reading $domtblout_file, unable to parse coords= token $coords_str on line\n$line", 1, $FH_HR);
      }      

      my $frame = undef;
      if($frame_str =~ /^frame=([123456])$/) { 
         $frame = $1;
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, reading $domtblout_file, unable to parse frame= token $frame_str on line\n$line", 1, $FH_HR);
      }      

      my $seq_strand = undef;
      if($seq_ftr_idx == -1) { 
        $seq_strand = "+";
        $source_coords = "1.." . $seq_len . ":+";
      }
      else { 
        $seq_strand = vdr_FeatureSummaryStrand($ftr_info_AHR->[$seq_ftr_idx]{"coords"}, $FH_HR);
        # $source_coords was defined by helper_protein_validation_breakdown_source() call above
      }
      # debugging print statements:
      #print("line:$line\n");
      #print("\torig_seq_len: $orig_seq_len\n");
      #print("\tmdl_name:     $mdl_name\n");
      #print("\tmdl_len:      $mdl_len\n");
      #print("\thit_ieval:  $hit_ieval\n");
      #print("\thit_score:  $hit_score\n");
      #print("\thit_bias:   $hit_bias\n");
      #print("\thmm_from:   $hmm_from\n");
      #print("\thmm_to:     $hmm_to\n");
      #print("\tali_from:   $ali_from\n");
      #print("\tali_to:     $ali_to\n");
      #print("\tenv_from:   $env_from\n");
      #print("\tenv_to:     $env_to\n");
      #print("\tsource_str: $source_str\n");
      #print("\tsource_val: $source_val\n");
      #print("\tcoords_str: $coords_str\n");
      #print("\tqlen_str:   $qlen_str\n");
      #print("\tframe_str:  $frame_str\n");
      #print("\t\tseq_ftr_idx:  $seq_ftr_idx\n");
      #print("\t\tsource_coords: $source_coords\n");
      #print("\t\torf_start:  $orf_start\n");
      #print("\t\torf_end:    $orf_end\n");
      #print("\t\tframe:      $frame\n");
      #print("\t\tseq_len:    $seq_len\n");
      #print("\t\tseq_name:   $seq_name\n");
      #print("\t\tseq_ftr_type_idx:   $seq_ftr_type_idx\n");

      my $seq_source_len_nt = vdr_CoordsLength($source_coords, $FH_HR); # length of predicted CDS or full sequence
      my $orf_coords = sprintf("%d..%d:%s", $orf_start, $orf_end, ($orf_start < $orf_end) ? "+" : "-");
      #print("\t\t\tseq_source_len_nt:   $seq_source_len_nt\n");
      #print("\t\t\torf_coords:          $orf_coords\n");

      # convert orf coordinates from relative nt coords within $source_coords to absolute coords (1..seqlen)
      my $hmmer_orf_nt_coords = vdr_CoordsRelativeToAbsolute($source_coords, $orf_coords, $FH_HR);
      #print("\t\t\thmmer_orf_nt_coords: $hmmer_orf_nt_coords (nt)\n");

      # convert hmmer env amino acid coordinates from relative aa coords within $hmmer_orf_nt_coords to absolute coords (1..seqlen)
      my $env_aa_strand = ($env_from < $env_to) ? "+" : "-";
      # always set protein coords to +, vdr_CoordsProteinRelativeToAbsolute() requires it, we revcomp afterwards if necessary
      my $env_aa_coords = ($env_aa_strand eq "+") ? sprintf("%d..%d:+", $env_from, $env_to) : sprintf("%d..%d:+", $env_to, $env_from);
      my $hmmer_env_nt_coords  = vdr_CoordsProteinRelativeToAbsolute($hmmer_orf_nt_coords, $env_aa_coords, $FH_HR);
      if($env_aa_strand eq "-") { # rev comp
        $hmmer_env_nt_coords = vdr_CoordsReverseComplement($hmmer_env_nt_coords, 0, $FH_HR); # 0: don't do carrots
      }

      #print("\t\t\tseq_source_len_nt:   $seq_source_len_nt\n");
      #print("\t\t\torf_coords:          $orf_coords\n");
      #print("\t\t\thmmer_orf_nt_coords: $hmmer_orf_nt_coords (nt)\n");
      #print("\t\t\tenv_aa_coords:       $env_aa_coords\n");
      #print("\t\t\thmmer_env_nt_coords: $hmmer_env_nt_coords (nt)\n");

      my $hmmer_nsgm = 0;
      my @hmmer_start_A  = (); # fill these below, only if nec
      my @hmmer_stop_A   = (); # fill these below, only if nec
      my @hmmer_strand_A = (); # fill these below, only if nec
      my $hmmer_summary_strand = undef; 

      # should we store this model/sequence/hit trio?
      # we do if A, B, and C are all TRUE and one or both of D or E is TRUE
      #  A. this model/sequence pair is compatible (sequence is full sequence or correct CDS feature for this model) 
      #  B. this is the highest scoring hit for this model for this sequence
      #  C. query length (full length seq or predicted CDS) is at least <x> nt from --minpvlen
      # 
      #  D. hit score is above minimum (--hlonescore)
      #  E. hit overlaps by at least 1 nt with a nucleotide prediction
      my $a_true = (($seq_ftr_idx == -1) || ($seq_ftr_idx == $mdl_ftr_idx)) ? 1 : 0; # sequence is full sequence OR model is CDS that pertains to sequence
      my $b_true = ((! defined $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"p_score"}) ||  # first hit, so must be highest score
                    ($hit_score > $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"p_score"})) ? 1 : 0; # highest scoring hit
      my $c_true = ($seq_source_len_nt >= $hminntlen) ? 1 : 0; # length >= --minpvlen

      if($a_true && $b_true && $c_true) { 
        my $x_true = ($hit_score >= $hlonescore) ? 1 : 0;
        my $y_true = 0; 
        # only bother determining $y_true if $x_true is 0
        if(! $x_true) { 
          $hmmer_summary_strand = vdr_FeatureSummaryStrand($hmmer_env_nt_coords, $FH_HR);
          if((defined $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"n_strand"}) &&
             ($ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"n_strand"} eq $hmmer_summary_strand)) { 
            $hmmer_nsgm = vdr_FeatureStartStopStrandArrays($hmmer_env_nt_coords, \@hmmer_start_A, \@hmmer_stop_A, \@hmmer_strand_A, $FH_HR);
            my $noverlap = helper_protein_validation_check_overlap($ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx], 
                                                                   $hmmer_start_A[0], $hmmer_stop_A[($hmmer_nsgm-1)], $hmmer_summary_strand, $FH_HR);
            if($noverlap > 0) { $y_true = 1; }
          }
        }
        if($x_true || $y_true) { 
          if($hmmer_nsgm == 0) { # if != 0, we already called FeatureStartStopStrandArrays() above
            vdr_FeatureStartStopStrandArrays($hmmer_env_nt_coords, \@hmmer_start_A, \@hmmer_stop_A, \@hmmer_strand_A, $FH_HR);
          }
          if(! defined $hmmer_summary_strand) { # if defined, we already calculated it above
            $hmmer_summary_strand = vdr_FeatureSummaryStrand($hmmer_env_nt_coords, $FH_HR);
          }
          $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"p_qstart"}  = $hmmer_start_A[0];
          $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"p_qstop"}   = $hmmer_stop_A[($hmmer_nsgm-1)];
          $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"p_strand"} = $hmmer_summary_strand;
          $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"p_len"}    = vdr_CoordsLength($hmmer_env_nt_coords, $FH_HR);
          $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"p_query"}  = $source_val;
          $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"p_score"}  = $hit_score;
          $ftr_results_HAHR->{$seq_name}[$mdl_ftr_idx]{"p_frame"}  = convert_esl_translate_to_blastx_frame($frame, $FH_HR);
        } # end of 'if($x_true || $y_true)'
      } # end of 'if($a_true && $b_true && $c_true)'
    } # end of 'if($line !~ m/^\#/)'
  } # end of 'while($my $line = <IN>)'
  close(IN);

  return 0;
}

#################################################################
# Subroutine: convert_esl_translate_to_blastx_frame()
# Incept:     EPN, Wed Mar 18 09:31:35 2020
#
# Purpose:    Convert a esl-translate frame value to blastx
#             input  return-value
#             1      +1
#             2      +2
#             3      +3
#             4      -1
#             5      -2
#             6      -3
#
# Arguments: 
#  $esl_translate_frame: value 1-6
#  $FH_HR:               ref to hash of file handles
#
# Returns:    blastx frame value: +1, +2, +3, -1, -2, -3
#
# Dies:       If esl_translate_frame is not 1, 2, 3, 4, 5, or 6
#
################################################################# 
sub convert_esl_translate_to_blastx_frame {
  my $sub_name = "convert_esl_translate_to_blastx_frame";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($esl_translate_frame, $FH_HR) = (@_);

  if($esl_translate_frame !~ m/^[123456]$/) { 
    ofile_FAIL("ERROR in $sub_name, unexpected esl-translate frame value of $esl_translate_frame", 1, $FH_HR);
  }
  if($esl_translate_frame <= 3) { 
    # return "+" . $esl_translate_frame;
    return $esl_translate_frame;
  }
  # return "-" . ($esl_translate_frame - 3);
  return ($esl_translate_frame - 3);
}

#################################################################
# Subroutine: get_hmm_list_for_model()
# Incept:     EPN, Wed Mar 18 14:21:54 2020
#
# Purpose:    Fill an array with names of all protein HMMs in the
#             HMM library that pertain to a model. There is one
#             per CDS for that model.
#      
# Arguments: 
#  $ftr_info_AHR:   REF to array of hashes with information on the features, PRE-FILLED
#  $mdl_name:       name of model
#  $hmm_name_AR:    REF to array of HMM names, FILLED HERE 
#  $FH_HR:             REF to hash of file handles
#
# Returns:    void, adds to @{$hmm_name_AR}
#
# Dies:       If HMM has zero CDS features
#
################################################################# 
sub get_hmm_list_for_model { 
  my $sub_name = "get_hmm_list_for_model()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $mdl_name, $hmm_name_AR, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});
  my $nhmm = 0;
  for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS")) { 
      push(@{$hmm_name_AR}, $mdl_name . "/" . $ftr_info_AHR->[$ftr_idx]{"coords"});
      $nhmm++;
    }
  }
  if($nhmm == 0) { 
    ofile_FAIL("ERROR in $sub_name, no CDS features exist for model $mdl_name", 1, $FH_HR);
  }
  return;
}

#################################################################
# Subroutine: helper_protein_validation_breakdown_source()
# Incept:     EPN, Wed Dec 19 12:05:02 2018
#
# Purpose: Given a query name from blastx output or a source name
#          (<str> from source=<str>) from hmmer output, determine if it is
#          for an entire input sequence, or a feature sequence
#          fetched from an input sequence.  The way we tell is by
#          looking up $in_source in $seq_len_HR.  If it exists as a
#          key, then the source is an entire input sequence. If it
#          does not exist as a key, then it should be
#          "<seq_name>/<ftr_type_idx>/<coords>", and <seq_name>
#          should be a key in $seq_info_HAR.
#
# Arguments:
#   $in_source:    source string, query name from blastx output or <str>
#                  from source=<str> in description of esl-translate output
#                  of a nucleotide sequence
#   $seq_len_HR:   ref to hash, keys are possible sequence names, values are sequence lengths
#   $FH_HR:        ref to hash of file handles, including 'log'
#             
# Returns:  3 values:
#           <seq_name>:     name of input sequence $in_source corresponds to
#           <ftr_type_idx>: "" if $in_source eq <seq_name>, else a string 
#                           of feature type and index separated by a '.'
#                           e.g. 'CDS.4', derived from $in_source.
#           <len>:          length of sequence ($seq_len_HR->{<seq_name>} if 
#                           $in_source eq <seq_name> else length derived from
#                           coordinate ranges listed in <coords_str>.
#           <coords_str>:   <coords>str> parsed from $in_source
#
# Dies: If $in_source is either not a valid <seq_name> (key in %{$seq_len_HR})
#       and $in_source cannot be broken down into a valid <seq_name><ftr_type_idx><coords_str>
#       
#
#################################################################
sub helper_protein_validation_breakdown_source {
  my $sub_name  = "helper_protein_validation_breakdown_source";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_source, $seq_len_HR, $FH_HR) = (@_);

  # printf("\nin $sub_name, in_source: $in_source\n");
  
  # check for simple case, $in_source is a sequence name key in %{$seq_len_HR}
  if(defined $seq_len_HR->{$in_source}) { 
    return($in_source, "", $seq_len_HR->{$in_source}, "");
  }

  # else, break down $in_source into:
  # <seq_name>/<ftr_type_idx>/<coords>
  # example: NC_002549.1/CDS.4/6039..6923:+,6923..8068:+
  # sequence is NC_002549.1 but we fetched feature CDS#4 as coordinates 6039..6293:+,6293..8068:+ from it
  my $ret_seq_name     = undef;
  my $ret_ftr_type_idx = undef;
  my $ret_len          = undef;
  my $coords_str       = undef;
  if($in_source =~ /^(\S+)\/([^\/]+\.\d+)\/([^\/]+$)/) { 
    ($ret_seq_name, $ret_ftr_type_idx, $coords_str) = ($1, $2, $3);
    if(! defined $seq_len_HR->{$ret_seq_name}) { 
      ofile_FAIL("ERROR in $sub_name, problem parsing query $in_source, unexpected sequence name $ret_seq_name (does not exist in input sequence length hash)", 1, $FH_HR); 
    }
    $ret_len = vdr_CoordsLength($coords_str, $FH_HR);
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse query $in_source", 1, $FH_HR); 
  }

  return ($ret_seq_name, $ret_ftr_type_idx, $ret_len, $coords_str);
}

#################################################################
# Subroutine: helper_blastx_breakdown_max_indel_str()
# Incept:     EPN, Tue Apr  2 06:35:13 2019
#
# Purpose: Given a string of one or more indel strings in the format:
#          "Q<d1>:S<d2>[+-]<d3>", separated by ";" if more than one.
#          fill arrays with <d1>, <d2>, and <d3>.
#
# Arguments:
#   $in_str:       input max indel string returned from parse_blastx.pl
#   $qpos_AR:      ref to array of query positions to fill
#   $spos_AR:      ref to array of subject positions to fill
#   $len_AR:       ref to array of lengths to fill
#   $FH_HR:        ref to hash of file handles, including 'log'
#             
# Returns:  number of indel strings parsed (number of elements added to @{$qpos_AR}, 
#           @{$spos_AR} and @{$len_AR}.
#
# Dies: If $in_str is not parseable
#
#################################################################
sub helper_blastx_breakdown_max_indel_str {
  my $sub_name  = "helper_blastx_breakdown_max_indel_str";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($in_str, $qpos_AR, $spos_AR, $len_AR, $FH_HR) = (@_);

  # printf("in $sub_name, in_str: $in_str\n");

  my @str_A = split(";", $in_str); 
  foreach my $str (@str_A) { 
    if($str =~ /^Q(\d+)\:S(\d+)[\+\-](\d+)$/) { 
      if(defined $qpos_AR) { push(@{$qpos_AR}, $1); }
      if(defined $spos_AR) { push(@{$spos_AR}, $2); }
      if(defined $len_AR)  { push(@{$len_AR},  $3); }
    }
    else { 
      ofile_FAIL("ERROR, in $sub_name, unable to parse indel string $str parsed out of $in_str", 1, $FH_HR);
    }
  }

  return scalar(@str_A);
}

#################################################################
# Subroutine: helper_blastx_max_indel_token_to_alt_coords
# Incept:     EPN, Tue May 25 18:42:57 2021
#
# Purpose: Given Given a string of one or more indel strings in the format:
#          "Q<d1>:S<d2>[+-]<d3>", separated by ";" if more than one.
#          fill arrays with <d1>, <d2>, and <d3>.
#
# Arguments:
#   $is_insert:     '1' if this is an insert, '0' if delete
#   $spos:          sequence position of indel
#   $aa_mpos:       amino acid position of indel
#   $len:           length of indel in nucleotides
#   $is_feature:    '1' if sequence is a feature so $ins_spos is relative to the feature coords $ftr_coords
#                   '0' if sequence is a the full sequence (not a feature) so $ins_spos is absolute coords
#   $ftr_scoords:   feature sequence coordinates 
#   $ftr_mcoords:   feature model coordinates 
#   $ftr_strand:    feature summary strand
#   $blastx_strand: strand of blastx query 
#   $seq_len:       total length of sequence
#   $FH_HR:         ref to hash of file handles, including 'log'
#             
# Returns:  2 values:
#           $alt_scoords: alert sequence coords string
#           $alt_mcoords: alert model    coords string
#
#################################################################
sub helper_blastx_max_indel_token_to_alt_coords {
  my $sub_name  = "helper_blastx_max_indel_token_to_alt_coords";
  my $nargs_expected = 11;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($is_insert, $spos, $aa_mpos, $len, $is_feature, $ftr_scoords, $ftr_mcoords, $ftr_strand, $blastx_strand, $seq_len, $FH_HR) = (@_);

  #printf("in $sub_name\n\tspos: $spos\n\taa_mpos: $aa_mpos\n\tlen: $len\n\tis_feature: $is_feature\n\tftr_scoords: $ftr_scoords\n\tftr_mcoords: $ftr_mcoords\n\tblastx_strand: $blastx_strand\n\tseq_len: $seq_len\n\n");

  my $absolute_scoords = undef; # absolute sequence coords
  my $relative_scoords = undef; # relative sequence coords
  my $absolute_mpos    = undef; # absolute model (nt) position
  my $absolute_mstrand = undef; # model strand
  my $absolute_mcoords = undef; # absolute query (sequence) coords
  my $alt_scoords      = undef; # to return: alert sequence coords string
  my $alt_mcoords      = undef; # to return: alert model    coords string

  if($is_insert) { 
    # determine sequence coordinates
    # determine coordinates in full sequence, this is complicated by fact that blastx query can be full seq or a single feature
    if(($is_feature && ($ftr_strand eq "+")) || ((! $is_feature) && ($blastx_strand eq "+"))) { 
      # blast alignment is to positive strand
      $relative_scoords = vdr_CoordsSegmentCreate($spos, $spos + $len - 1, $blastx_strand, $FH_HR);
    }
    else { 
      # blast alignment is to negative strand
      $relative_scoords = vdr_CoordsSegmentCreate($spos, $spos - $len + 1, $blastx_strand, $FH_HR);
    }
    $absolute_scoords = ($is_feature) ? 
        $ftr_scoords : vdr_CoordsSegmentCreate(1, $seq_len, "+", $FH_HR); # the full sequence
    $alt_scoords = sprintf("seq:%s;", vdr_CoordsRelativeToAbsolute($absolute_scoords, $relative_scoords, $FH_HR));
    
    # determine model coordinates
    $absolute_mcoords = vdr_CoordsProteinRelativeToAbsolute($ftr_mcoords,
                                                            vdr_CoordsSegmentCreate($aa_mpos, $aa_mpos, "+", $FH_HR), $FH_HR);
    # $absolute_mcoords will now be a full codon, but we only want the final position
    (undef, $absolute_mpos, $absolute_mstrand) = vdr_CoordsSegmentParse($absolute_mcoords, $FH_HR);
    $absolute_mcoords = vdr_CoordsSegmentCreate($absolute_mpos, $absolute_mpos, $absolute_mstrand, $FH_HR); # yes, we want same start/end
    $alt_mcoords = sprintf("mdl:%s;", $absolute_mcoords);
  }
  else { # delete
    # determine sequence coordinates
    $relative_scoords = vdr_CoordsSegmentCreate($spos, $spos, $blastx_strand, $FH_HR);
    $absolute_scoords = ($is_feature) ? 
        $ftr_scoords : vdr_CoordsSegmentCreate(1, $seq_len, "+", $FH_HR); # the full sequence
    $alt_scoords = sprintf("seq:%s;", vdr_CoordsRelativeToAbsolute($absolute_scoords, $relative_scoords, $FH_HR));
    
    # determine model coordinates
    my $aa_len = $len / 3; 
    $absolute_mcoords = vdr_CoordsProteinRelativeToAbsolute($ftr_mcoords,
                                                            vdr_CoordsSegmentCreate($aa_mpos+1, $aa_mpos+$aa_len, "+", $FH_HR), $FH_HR);
    $alt_mcoords = sprintf("mdl:%s;", $absolute_mcoords);
  }
  
  return($alt_scoords, $alt_mcoords);
}

#################################################################
# Subroutine: helper_protein_validation_db_seqname_to_ftr_idx()
# Incept:     EPN, Tue Dec 18 13:27:50 2018
#
# Purpose:    Find the feature $ftr_idx that corresponds to the blastx
#             db sequence that was named with the convention:
#
#             <protein-accession>/<coords-str>
#
#             Where <coords-str> is identical to $ftr_info_AHR->{"ref_coords"}[$ftr_idx].
#
# Arguments: 
#  $blastx_seqname: sequence name
#  $ftr_info_AHR:   ref to the feature info array of hashes 
#  $FH_HR:          ref to hash of file handles
#
# Returns:    Two values:
#             <$ftr_idx>:    index of matching feature
#             <$ftr_strand>: summary strand of feature <$ftr_idx>
#
# Dies:       If we find zero features that match to this sequence
#             If we find more than 1 features that match to this sequence
#
################################################################# 
sub helper_protein_validation_db_seqname_to_ftr_idx { 
  my $sub_name = "helper_protein_validation_db_seqname_to_ftr_idx";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($blastx_seqname, $ftr_info_AHR, $FH_HR) = @_;

  my $nftr = scalar(@{$ftr_info_AHR});

  my $ret_ftr_idx = undef;
  my $ret_strand  = undef;

  if($blastx_seqname =~ /(\S+)\/(\S+)/) { 
    my ($accn, $coords) = ($1, $2);
    # find it in @{$ftr_info_AHR->{"coords"}}
    $ret_strand = vdr_FeatureSummaryStrand($coords, $FH_HR);
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if(($ftr_info_AHR->[$ftr_idx]{"type"} eq "CDS")) { 
        if($ftr_info_AHR->[$ftr_idx]{"coords"} eq $coords) { 
          if(defined $ret_ftr_idx) { # found more than 1 features that match
            ofile_FAIL("ERROR in $sub_name, found blastx db sequence with coords that match two features, ftr_idx: $ftr_idx and $ret_ftr_idx", 1, $FH_HR);
          }                  
          $ret_ftr_idx = $ftr_idx;
        }
      }
    }
    if(! defined $ret_ftr_idx) { # did not find match
      ofile_FAIL("ERROR in $sub_name, did not find matching feature for blastx db sequence $blastx_seqname", 1, $FH_HR);
    }
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, unable to parse blastx db sequence name $blastx_seqname", 1, $FH_HR); 
  }

  return ($ret_ftr_idx, $ret_strand);
}

#################################################################
# Subroutine: helper_protein_validation_check_overlap()
# Incept:     EPN, Sat Mar 21 09:14:34 2020
#
# Purpose:    Check if a protein validation hit overlaps with 
#             a nucleotide prediction.
#      
# Arguments: 
#  $ftr_results_HR: REF to hash of $ftr_results for a specific sequence 
#                   and feature index
#  $pv_start:       protein validation predicted start position
#  $pv_stop:        protein validation predicted stop position
#  $pv_strand:      protein validation predicted strand
#  $FH_HR:          REF to hash of file handles
#
# Returns:    number of nucleotide overlap, 0 if none
#
# Dies:       If seq_Overlap() has a problem or if $pv_strand is not "+" or "-"
#
################################################################# 
sub helper_protein_validation_check_overlap { 
  my $sub_name = "helper_protein_validation_check_overlap";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_results_HR, $pv_start, $pv_stop, $pv_strand, $FH_HR) = @_;

  my $noverlap = 0;
  if($pv_strand eq "+") { 
    ($noverlap, undef) = seq_Overlap($ftr_results_HR->{"n_start"},
                                     $ftr_results_HR->{"n_stop"},
                                     $pv_start, $pv_stop, $FH_HR);
  }
  elsif($pv_strand eq "-") { 
    ($noverlap, undef) = seq_Overlap($ftr_results_HR->{"n_stop"},
                                     $ftr_results_HR->{"n_start"},
                                     $pv_stop, $pv_start, $FH_HR);
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, pv_strand is $pv_strand (expected to be + or -)\n", 1, $FH_HR);
  }

  return $noverlap;
}

#################################################################
#
# Other subroutines related to alerts: 
# alert_list_option
# alert_feature_instance_add 
# alert_sequence_instance_add 
# alert_add_unexdivg 
# alert_instances_check_prevents_annot
#
#################################################################
# Subroutine:  alert_list_option()
# Incept:      EPN, Wed Apr  3 11:11:26 2019
#
# Purpose:    Handle the --alt_list option by outputting a list 
#             of all alerts and whether they cause failure or not
#             by default to stdout.
#
# Arguments: 
#  $alt_info_HHR:  REF to the alert info hash of arrays, PRE-FILLED
#  $pkgname:       package name ("VADR"
#  $version:       version
#  $releasedate:   release date
#
# Returns:    void
#
# Dies:       never
#
#################################################################
sub alert_list_option { 
  my $sub_name = "alert_list_option()"; 
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($alt_info_HHR, $pkgname, $version, $releasedate) = @_;
  
  my $div_line = utl_StringMonoChar(60, "#", undef) . "\n";
  
  # determine order of codes to print
  my @code_A = ();
  alert_order_arrays($alt_info_HHR, \@code_A, undef, undef);
  
  my $idx = 0;
  my $code = undef;

  my @data_AA  = ();
  my @head_AA  = ();
  my @bcom_A   = ();

  @{$head_AA[0]} = ("",    "",       "",    "misc, not",      "",            "");
  @{$head_AA[1]} = ("",    "alert",  "",    "failure",        "short",       "long");
  @{$head_AA[2]} = ("idx", "code",   "S/F", "(if in .minfo)", "description", "description");


  push(@bcom_A, $div_line);
  push(@bcom_A, "#\n");
  push(@bcom_A, "# $pkgname $version ($releasedate)\n");
  push(@bcom_A, "#\n");
  push(@bcom_A, "# Alert codes that ALWAYS cause a sequence to FAIL, and cannot be\n");
  push(@bcom_A, "# listed in --alt_pass, --alt_fail, --alt_mnf_yes, or --alt_mnf_no\n");
  push(@bcom_A, "# option strings:\n#\n");
  $idx = 0;
  foreach $code (@code_A) { 
    if($alt_info_HHR->{$code}{"always_fails"}) { 
      $idx++;
      push(@data_AA, [$idx, $code, 
                      ($alt_info_HHR->{$code}{"pertype"} eq "sequence" ? "S" : "F"), 
                      "never",
                      helper_tabular_replace_spaces($alt_info_HHR->{$code}{"sdesc"}), 
                      $alt_info_HHR->{$code}{"ldesc"}]);
    }
  }
  ofile_TableHumanOutput(\@data_AA, \@head_AA, undef, \@bcom_A, undef, "  ", "-", "#", "", "", "-", *STDOUT, undef, undef);

  @bcom_A  = ();
  @data_AA = ();
  push(@bcom_A, "#\n");
  push(@bcom_A, $div_line);
  push(@bcom_A, "#\n");
  push(@bcom_A, "# Alert codes that cause a sequence to FAIL by default, but can be set\n");
  push(@bcom_A, "# to not FAIL a sequence by listing the code as part of a comma separated\n");
  push(@bcom_A, "# string of codes in <s> with the --alt_pass <s> option:\n#\n");
  $idx = 0;
  foreach $code (@code_A) { 
    if(($alt_info_HHR->{$code}{"causes_failure"}) && 
       (! $alt_info_HHR->{$code}{"always_fails"})) { 
      $idx++;
      my $misc_not_fail_str = "never";
      if($alt_info_HHR->{$code}{"pertype"} eq "feature") { 
        $misc_not_fail_str = $alt_info_HHR->{$code}{"misc_not_failure"} ? "yes" : "no";
      }
      push(@data_AA, [$idx, $code, 
                      ($alt_info_HHR->{$code}{"pertype"} eq "sequence" ? "S" : "F"), 
                      $misc_not_fail_str,
                      helper_tabular_replace_spaces($alt_info_HHR->{$code}{"sdesc"}), 
                      $alt_info_HHR->{$code}{"ldesc"}]);
    }
  }
  ofile_TableHumanOutput(\@data_AA, \@head_AA, undef, \@bcom_A, undef, "  ", "-", "#", "", "", "-", *STDOUT, undef, undef);

  @bcom_A  = ();
  @data_AA = ();
  push(@bcom_A, "#\n");
  push(@bcom_A, $div_line);
  push(@bcom_A, "#\n");
  push(@bcom_A, "# Alert codes that do not cause a sequence to FAIL by default, but can be set\n");
  push(@bcom_A, "# to FAIL a sequence by listing the code as part of a comma separated\n");
  push(@bcom_A, "# string of codes in <s> with the --alt_fail <s> option:\n#\n");
  $idx = 0;
  foreach $code (@code_A) { 
    if((! $alt_info_HHR->{$code}{"causes_failure"}) && 
       (! $alt_info_HHR->{$code}{"always_fails"})) { 
      $idx++;
      my $misc_not_fail_str = "never";
      if($alt_info_HHR->{$code}{"pertype"} eq "feature") { 
        $misc_not_fail_str = $alt_info_HHR->{$code}{"misc_not_failure"} ? "yes" : "no";
      }
      push(@data_AA, [$idx, $code, 
                      ($alt_info_HHR->{$code}{"pertype"} eq "sequence" ? "S" : "F"), 
                      $misc_not_fail_str,
                      helper_tabular_replace_spaces($alt_info_HHR->{$code}{"sdesc"}), 
                      $alt_info_HHR->{$code}{"ldesc"}]);
    }
  }
  ofile_TableHumanOutput(\@data_AA, \@head_AA, undef, \@bcom_A, undef, "  ", "-", "#", "", "", "-", *STDOUT, undef, undef);

  return;
}

#################################################################
# Subroutine:  alert_pass_fail_options()
# Incept:      EPN, Wed Apr  3 12:51:25 2019
#
# Purpose:    Handle the --alt_pass and --alt_fail options by 
#             parsing their strings, determining if they are valid
#             and updating the "causes_failure" values in 
#             %{$alt_info_HHR}.
#
# Arguments: 
#  $alt_info_HHR:   REF to the alert info hash of arrays, PRE-FILLED
#  $opt_HHR:        REF to 2D hash of option values
#
# Returns:    void
#
# Dies:       if --alt_pass or --alt_fail option strings are invalid
#
#################################################################
sub alert_pass_fail_options { 
  my $sub_name = "alert_pass_fail_options()"; 
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($alt_info_HHR, $opt_HHR) = @_;
  
  my @pass_A = ();
  my @fail_A = ();
  if(opt_IsUsed("--alt_pass", $opt_HHR)) { 
    @pass_A = split(",", opt_Get("--alt_pass", $opt_HHR));
  }
  if(opt_IsUsed("--alt_fail", $opt_HHR)) { 
    @fail_A = split(",", opt_Get("--alt_fail", $opt_HHR));
  }

  my $die_str = "";
  my $alt_code = undef;

  # --alt_pass codes
  foreach my $alt_code (@pass_A) { 
    if(! defined $alt_info_HHR->{$alt_code}) { 
      $die_str .= "alert code $alt_code is invalid (does not exist)\n";
    }
    elsif($alt_info_HHR->{$alt_code}{"always_fails"}) { 
      $die_str .= "alert code $alt_code always causes failure, it cannot be listed in --alt_pass string\n";
    }
    else { 
      vdr_AlertInfoSetCausesFailure($alt_info_HHR, $alt_code, 0, undef);
    }
  }

  # --alt_fail codes
  foreach my $alt_code (@fail_A) { 
    if(! defined $alt_info_HHR->{$alt_code}) { 
      $die_str .= "alert code $alt_code is invalid (does not exist)\n";
    }
    else { 
      vdr_AlertInfoSetCausesFailure($alt_info_HHR, $alt_code, 1, undef);
    }
  }

  if($die_str ne "") { 
    $die_str .= "Use the --alt_list to see a list possible alert codes\nto use with --alt_pass and --alt_fail.\n";
    ofile_FAIL("ERROR processing --alt_fail and/or --alt_pass options:\n$die_str", 1, undef);
  }
  
  return;
}

#################################################################
# Subroutine:  alert_misc_not_failure_options()
# Incept:      EPN, Tue Feb  9 13:43:36 2021
#
# Purpose:    Handle the --alt_mnf_yes and --alt_mnf_no options by 
#             parsing their strings, determining if they are valid
#             and updating the "misc_not_fail" values in 
#             %{$alt_info_HHR}. Also make sure they make sense with
#             the "causes_failure" values in %{$alt_info_HHR}: 
#             "misc_not_failure" can only be '1' for an alert code 
#             if "cause_failure" is '1'.
#
#             This subroutine should be called *after* 
#             alert_pass_fail_options().
#
# Arguments: 
#  $alt_info_HHR:   REF to the alert info hash of arrays, PRE-FILLED
#  $opt_HHR:        REF to 2D hash of option values
#
# Returns:    void
#
# Dies:       if --alt_mnf_yes or --alt_mnf_no option strings are invalid
#             if trying to set alert "misc_not_failure" value conflicts
#             with "causes_failure" value as described in "Purpose".
#
#################################################################
sub alert_misc_not_failure_options { 
  my $sub_name = "alert_misc_not_failure_options()"; 
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($alt_info_HHR, $opt_HHR) = @_;
  
  my @yes_A = ();
  my @no_A  = ();
  if(opt_IsUsed("--alt_mnf_yes", $opt_HHR)) { 
    @yes_A = split(",", opt_Get("--alt_mnf_yes", $opt_HHR));
  }
  if(opt_IsUsed("--alt_mnf_no", $opt_HHR)) { 
    @no_A = split(",", opt_Get("--alt_mnf_no", $opt_HHR));
  }

  my $die_str = "";
  my $alt_code = undef;

  # --alt_mnf_yes codes
  foreach my $alt_code (@yes_A) { 
    if(! defined $alt_info_HHR->{$alt_code}) { 
      $die_str .= "alert code $alt_code is invalid (does not exist)\n";
    }
    elsif($alt_info_HHR->{$alt_code}{"always_fails"}) { 
      $die_str .= "alert code $alt_code always causes failure, it cannot be listed in --alt_mnf_yes string\n";
    }
    elsif(! $alt_info_HHR->{$alt_code}{"causes_failure"}) { 
      $die_str .= "alert code $alt_code does *not* cause failure (either by default or due to --alt_fail option), it cannot be listed in --alt_mnf_yes string\n";
    }
    else { 
      vdr_AlertInfoSetMiscNotFailure($alt_info_HHR, $alt_code, 1, undef);
    }
  }

  # --alt_mnf_no codes
  foreach my $alt_code (@no_A) { 
    if(! defined $alt_info_HHR->{$alt_code}) { 
      $die_str .= "alert code $alt_code is invalid (does not exist)\n";
    }
    else { 
      vdr_AlertInfoSetMiscNotFailure($alt_info_HHR, $alt_code, 0, undef);
    }
  }

  if($die_str ne "") { 
    $die_str .= "Use the --alt_list to see a list possible alert codes\nto use with --alt_pass and --alt_fail.\n";
    ofile_FAIL("ERROR processing --alt_mnf_yes and/or --alt_mnf_no options:\n$die_str", 1, undef);
  }
  
  return;
}

#################################################################
# Subroutine:  alert_order_arrays()
# Incept:      EPN, Thu Apr  4 12:59:41 2019
#
# Purpose:    Fill arrays with the order of alerts for outputting.
#             Alerts that DO NOT cause failure are listed first, 
#             then alerts that DO cause failure. 
#
# Arguments: 
#  $alt_info_HHR:  ref to the alert info hash of arrays, PRE-FILLED
#  $all_order_AR:  ref to array of order of all alerts, FILLED HERE if defined
#  $seq_order_AR:  ref to array of order of per-sequence alerts, FILLED HERE if defined
#  $ftr_order_AR:  ref to array of order of per-feature alerts, FILLED HERE if defined
#
# Returns:    void
#
# Dies:       never
#
#################################################################
sub alert_order_arrays { 
  my $sub_name = "alert_order_arrays()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($alt_info_HHR, $all_order_AR, $seq_order_AR, $ftr_order_AR) = @_;
  
  my @all_A = ();
  my @all_pass_A = ();
  my @all_fail_A = ();
  my @seq_pass_A = ();
  my @seq_fail_A = ();
  my @ftr_pass_A = ();
  my @ftr_fail_A = ();
  my $alt_code = undef;
  # first get array using "order" key
  foreach $alt_code (sort keys (%{$alt_info_HHR})) { 
    $all_A[($alt_info_HHR->{$alt_code}{"order"})] = $alt_code;
  }

  foreach $alt_code (@all_A) { 
    if($alt_info_HHR->{$alt_code}{"causes_failure"}) { 
      push(@all_fail_A, $alt_code); 
      if($alt_info_HHR->{$alt_code}{"pertype"} eq "sequence") { 
        push(@seq_fail_A, $alt_code); 
      }
      if($alt_info_HHR->{$alt_code}{"pertype"} eq "feature") { 
        push(@ftr_fail_A, $alt_code); 
      }
    }
    else { # does not cause failure
      push(@all_pass_A, $alt_code); 
      if($alt_info_HHR->{$alt_code}{"pertype"} eq "sequence") { 
        push(@seq_pass_A, $alt_code); 
      }
      if($alt_info_HHR->{$alt_code}{"pertype"} eq "feature") { 
        push(@ftr_pass_A, $alt_code); 
      }
    }
  }

  if(defined $all_order_AR) { @{$all_order_AR} = (@all_pass_A, @all_fail_A); }
  if(defined $seq_order_AR) { @{$seq_order_AR} = (@seq_pass_A, @seq_fail_A); }
  if(defined $ftr_order_AR) { @{$ftr_order_AR} = (@ftr_pass_A, @ftr_fail_A); }

  return;
}

#################################################################
# Subroutine:  alert_sequence_instance_add()
# Incept:      EPN, Thu Apr  4 11:45:05 2019
#
# Purpose:    Add an $alt_code alert to the @{$alt_seq_instances_HHR} 
#             for sequence name $seq_name.
#
# Arguments: 
#  $alt_seq_instances_HHR:  REF to per-sequence alert instances to add to, ADDED TO HERE (maybe),
#                           can be undef if $alt_info_HHR->{$alt_code}{"pertype"} is "feature".
#  $alt_info_HHR:           REF to the alert info hash of arrays, PRE-FILLED
#  $alt_code:               alert code we're adding an alert for
#  $seq_name:               sequence name
#  $value:                  value to add as $alt_seq_instances_HHR->{$seq_name}{$code}
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies:       - If we find an alert instance incompatibility.
#             - if pertype of $alt_code is not "sequence"
#
#################################################################
sub alert_sequence_instance_add { 
  my $sub_name = "alert_sequence_instance_add()";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($alt_seq_instances_HHR, $alt_info_HHR, $alt_code, $seq_name, $value, $FH_HR) = @_;

  # my $tmp_altmsg = "alt_code: $alt_code, seq_name: $seq_name, value: $value\n";
  # printf("in $sub_name $tmp_altmsg\n");
  
  if(! defined $alt_seq_instances_HHR) { 
    ofile_FAIL("ERROR in $sub_name, alt_seq_instances_HHR is undefined", 1, $FH_HR);
  }
  if(! defined $alt_info_HHR->{$alt_code}) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized alert code $alt_code", 1, $FH_HR);
  }
  if($alt_info_HHR->{$alt_code}{"pertype"} ne "sequence") { 
    ofile_FAIL("ERROR in $sub_name alert code $alt_code is not a per-sequence alert", 1, $FH_HR);
  }
  if(! defined $value) { 
    ofile_FAIL("ERROR in $sub_name, value is undefined", 1, $FH_HR);
  }
  if($value eq "") { 
    ofile_FAIL("ERROR in $sub_name, value is empty string", 1, $FH_HR);
  }

  if(! defined $alt_seq_instances_HHR->{$seq_name}) { 
    %{$alt_seq_instances_HHR->{$seq_name}} = (); 
  }

  # if this alert already exists (rare case), add to it
  if(defined $alt_seq_instances_HHR->{$seq_name}{$alt_code}) { 
    $alt_seq_instances_HHR->{$seq_name}{$alt_code} .= ":VADRSEP:" . $value; 
  }
  else { # if it doesn't already exist (normal case), create it
    $alt_seq_instances_HHR->{$seq_name}{$alt_code} = $value; 
  }

  return;
}

#################################################################
# Subroutine: alert_sequence_instance_update_mdl_coords()
# Incept:     EPN, Tue Jun  8 14:21:01 2021
# Purpose:    Updates model coordinates for a per-sequence alert
#             given the @{$ua2rf_AR} map of unaligned positions to 
#             model coordinates.
#
# Arguments:
#  $alt_seq_instances_HHR: REF to 2D hash with per-sequence alerts, PRE-FILLED
#  $alt_info_HHR           REF to the alert info hash of arrays, PRE-FILLED
#  $alt_code:              alert code we're adding an alert for
#  $seq_name:              sequence name
#  $ua2rf_AR:              REF to array that maps unaligned positions to reference positions
#                          [1..$uapos..$ualen]: reference position that unaligned position $uapos aligns to 
#                          if $ua2rf_A[$uapos] <  0, $uapos inserts *after* ref posn (-1 * $ua2rf_A[$uapos])
#                          if $ua2rf_A[$uapos] == 0, $uapos inserts *before* ref posn 1
#                          $ua2rf_A[0] is invalid (set to 0)
#  $FH_HR:                 REF to hash of file handles
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub alert_sequence_instance_update_mdl_coords { 
  my $sub_name = "alert_sequence_instance_update_mdl_coords";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($alt_seq_instances_HHR, $alt_info_HHR, $alt_code, $seq_name, $ua2rf_AR, $FH_HR) = @_;

  if(! defined $alt_seq_instances_HHR) { 
    ofile_FAIL("ERROR in $sub_name, alt_seq_instances_HHR is undefined", 1, $FH_HR);
  }
  if(! defined $alt_info_HHR->{$alt_code}) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized alert code $alt_code", 1, $FH_HR);
  }
  if($alt_info_HHR->{$alt_code}{"pertype"} ne "sequence") { 
    ofile_FAIL("ERROR in $sub_name alert code $alt_code is not a per-sequence alert", 1, $FH_HR);
  }
  if(! defined $ua2rf_AR) { 
    ofile_FAIL("ERROR in $sub_name, ua2rf_AR is undefined", 1, $FH_HR);
  }

  if(defined $alt_seq_instances_HHR->{$seq_name}) { 
    my $alt_instance = undef;
    $alt_instance = alert_sequence_instance_fetch($alt_seq_instances_HHR, $seq_name, $alt_code);
    if(defined $alt_instance) { 
      my $new_alt_instance = "";
      my @instance_str_A = split(":VADRSEP:", $alt_instance);
      my $instance_ctr = 0;
      foreach my $instance_str (@instance_str_A) { 
        my ($instance_scoords, $instance_mcoords, $instance_detail) = alert_instance_parse($instance_str);
        # alert_instance_parse removes 'seq:' and ';', 'mdl:', and ';'
        my $new_mcoords = "mdl:";
        if($instance_scoords ne "VADRNULL") { 
          my @start_A  = (); 
          my @stop_A   = (); 
          my @strand_A = ();
          vdr_FeatureStartStopStrandArrays($instance_scoords, \@start_A, \@stop_A, \@strand_A, $FH_HR);
          for(my $i = 0; $i < scalar(@start_A); $i++) { 
            if($i > 0) { $new_mcoords .= ","; }
            $new_mcoords .= vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$start_A[$i]]), 
                                                    abs($ua2rf_AR->[$stop_A[$i]]), 
                                                    $strand_A[$i], $FH_HR);
          }
          $new_mcoords .= ";";
        }
        else { 
          $new_mcoords .= "VADRNULL;";
        }
        if($instance_ctr > 0) { $new_alt_instance .= ":VADRSEP:"; }
        $new_alt_instance .= "seq:" . $instance_scoords . ";" . $new_mcoords . $instance_detail;
      }
      # update the alert instance
      $alt_seq_instances_HHR->{$seq_name}{$alt_code} = $new_alt_instance;
    }
  }

  return;
}

#################################################################
# Subroutine:  alert_feature_instance_add()
# Incept:      EPN, Thu Apr  4 11:40:09 2019
#
# Purpose:    Add an $alt_code alert to the @{$alt_ftr_instances_HHHR} 
#             for sequence name $seq_name and feature index $ftr_idx.
#
# Arguments: 
#  $alt_ftr_instances_HHHR: REF to per-feature alert instances to add to, ADDED TO HERE
#  $alt_info_HHR:           REF to the alert info hash of arrays, PRE-FILLED
#  $alt_code:               alert code we're adding an alert for
#  $seq_name:               sequence name
#  $ftr_idx:                feature index
#  $value:                  value to add as $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}{$code}
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    void
#
# Dies:       - If we find an alert instance incompatibility.
#             - if pertype of $alt_code is not "feature"
#
#################################################################
sub alert_feature_instance_add { 
  my $sub_name = "alert_feature_instance_add()";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($alt_ftr_instances_HHHR, $alt_info_HHR, $alt_code, $seq_name, $ftr_idx, $value, $FH_HR) = @_;

  # my $tmp_altmsg = "ftr_idx: $ftr_idx, alt_code: $alt_code, seq_name: $seq_name, value: $value\n";
  # printf("in $sub_name $tmp_altmsg\n");
  
  if(! defined $alt_ftr_instances_HHHR) { 
    ofile_FAIL("ERROR in $sub_name, but alt_ftr_instances_HHHR is undefined", 1, $FH_HR);
  }
  if(! defined $alt_info_HHR->{$alt_code}) { 
    ofile_FAIL("ERROR in $sub_name, unrecognized alert code $alt_code", 1, $FH_HR);
  }
  if($alt_info_HHR->{$alt_code}{"pertype"} ne "feature") { 
    ofile_FAIL("ERROR in $sub_name alert code $alt_code is not a per-feature alert", 1, $FH_HR);
  }
  if(! defined $value) { 
    ofile_FAIL("ERROR in $sub_name, value is undefined", 1, $FH_HR);
  }
  if($value eq "") { 
    ofile_FAIL("ERROR in $sub_name, value is empty string", 1, $FH_HR);
  }

  if(! defined $alt_ftr_instances_HHHR->{$seq_name}) { 
    %{$alt_ftr_instances_HHHR->{$seq_name}} = (); 
  }
  if(! defined $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}) { 
    %{$alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}} = (); 
  }

  # if this alert already exists (rare case), add to it
  if(defined $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}{$alt_code}) { 
    $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}{$alt_code} .= ":VADRSEP:" . $value; 
  }
  else { # if it doesn't already exist (normal case), create it
    $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}{$alt_code} = $value; 
  }

  return;
}

#################################################################
# Subroutine: alert_sequence_instance_fetch()
# Incept:     EPN, Fri Apr  5 12:02:29 2019
# Purpose:    Return $alt_seq_instance_HHR->{$seq_name}{$alt_code} 
#             if it is defined, else returns undef.
#
# Arguments:
#  $alt_seq_instances_HHHR: ref to 2D hash of sequence alert instances
#  $seq_name:               the sequence name
#  $alt_code:               the alert code.
#
# Returns:  $alt_seq_instance_HHHR->{$seq_name}{$alt_code} 
#           undef if that is not defined
#
# Dies:     never
#
#################################################################
sub alert_sequence_instance_fetch { 
  my $sub_name = "alert_sequence_instance_fetch";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($alt_seq_instances_HHHR, $seq_name, $alt_code) = (@_);

  if((defined $alt_seq_instances_HHHR) && 
     (defined $alt_seq_instances_HHHR->{$seq_name}) && 
     (defined $alt_seq_instances_HHHR->{$seq_name}{$alt_code})) { 
    return $alt_seq_instances_HHHR->{$seq_name}{$alt_code}; # defined
  }

  return undef; # not defined
}

#################################################################
# Subroutine: alert_feature_instance_fetch()
# Incept:     EPN, Fri Apr  5 11:59:00 2019
# Purpose:    Return $alt_ftr_instance_HHHR->{$seq_name}{$ftr_idx}{$alt_code} 
#             if it is defined, else returns undef.
#
# Arguments:
#  $alt_ftr_instances_HHHR: ref to 3D hash of feature alert instances
#  $seq_name:               the sequence name
#  $ftr_idx:                the feature index
#  $alt_code:               the alert code.
#
# Returns:  $alt_ftr_instance_HHHR->{$seq_name}{$ftr_idx}{$alt_code} 
#           undef if that is not defined
#
# Dies:     never
#
#################################################################
sub alert_feature_instance_fetch { 
  my $sub_name = "alert_feature_instance_fetch";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($alt_ftr_instances_HHHR, $seq_name, $ftr_idx, $alt_code) = (@_);

  if((defined $alt_ftr_instances_HHHR) && 
     (defined $alt_ftr_instances_HHHR->{$seq_name}) && 
     (defined $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}) && 
     (defined $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}{$alt_code})) { 
    return $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}{$alt_code}; # defined
  }

  return undef; # not defined
}

#################################################################
# Subroutine: alert_instance_parse()
# Incept:     EPN, Tue May 18 09:28:13 2021
# Purpose:    Parse a (feature or sequence) alert instance 
#             string and return:
#             <scoords> (seq coords string)
#             <mcoords> (model coords string)
#             <detail>  (detail on the coords string)
#
# Arguments:
#  $alt_instance_str: the alert instance string
#
# Returns:  3 values: 
#           <scoords>: sequence coords, or 'VADRNULL' if none
#           <mcoords>: model coords,    or 'VADRNULL' if none
#           <detail>:  detail string,   or 'VADRNULL' if none
#
# Dies:     never
#
#################################################################
sub alert_instance_parse { 
  my $sub_name = "alert_instance_parse";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($alt_instance_str) = (@_);
  
  my $scoords = "VADRNULL";
  my $mcoords = "VADRNULL";
  my $detail  = "VADRNULL";

  # printf("in $sub_name, alt_instance_str: $alt_instance_str\n");

  if($alt_instance_str =~ /^seq\:([^\;]+);mdl\:([^\;]+);(.*)$/) { 
    ($scoords, $mcoords, $detail) = ($1, $2, $3);
  }
  elsif($alt_instance_str eq "VADRNULL") { 
    $detail = "VADRNULL";
  }
  else { 
    $detail = $alt_instance_str;
  }

  # printf("returning: scoords: $scoords mcoords: $mcoords detail: $detail\n");

  return ($scoords, $mcoords, $detail);
}

#################################################################
# Subroutine: alert_add_parent_based()
# Incept:     EPN, Fri Mar 27 06:37:18 2020
# Purpose:    Adds alerts to children features that have 
#             parents features with fatal alerts.
#
# Arguments:
#  $seq_name_AR:             REF to array of sequence names, PRE-FILLED
#  $ftr_info_AHR:            REF to array of hashes with information on the features, PRE-FILLED
#  $alt_info_HHR:            REF to array of hashes with information on the alerts, PRE-FILLED
#  $ftr_results_HAHR:        REF to feature results HAH, PRE-FILLED
#  $alt_ftr_instances_HHHR:  REF to array of 2D hashes with per-feature alerts, PRE-FILLED
#  $parent_type:             feature type of parent (e.g. "CDS")
#  $child_type:              feature type of child  (e.g. "mat_peptide")
#  $alt_code:                alert code (e.g. "peptrans")
#  $alt_msg:                 message for alert code instance (often "VADRNULL")
#  $opt_HHR:                 REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $FH_HR:                   REF to hash of file handles, including 'log'
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub alert_add_parent_based { 
  my $sub_name = "alert_add_parent_based";
  my $nargs_exp = 11;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name_AR, $ftr_info_AHR, $alt_info_HHR, $ftr_results_HAHR, 
      $alt_ftr_instances_HHHR, $parent_type, $child_type, $alt_code, $alt_msg, $opt_HHR, $FH_HR) = @_;

  my $nseq = scalar(@{$seq_name_AR});
  my $nftr = scalar(@{$ftr_info_AHR});

  # get children info for all features, we'll use this in the loop below
  my @children_AA = ();
  vdr_FeatureInfoChildrenArrayOfArrays($ftr_info_AHR, $child_type, undef, \@children_AA, $FH_HR);

  # get array of feature indices that are of type $parent_type
  my @parent_ftr_idx_A = ();
  my $ftr_idx; 
  for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
    if(($ftr_info_AHR->[$ftr_idx]{"type"} eq $parent_type) && # feature $ftr_idx is of type $parent_type
       (scalar(@{$children_AA[$ftr_idx]}) > 0)) {             # feature $ftr_idx has at least one child of type $child_type
      push(@parent_ftr_idx_A, $ftr_idx);
    }
  }

  my $nparent_ftr = scalar(@parent_ftr_idx_A); # number of features of type $parent_type with >= 1 child of type $child_type
  if($nparent_ftr == 0) { 
    # no features of type $parent_idx with children of type $child_type, 
    # so no children to add alerts for, return
    return;
  }

  # get array of all fatal feature alert types 
  # note: we don't care about which feature/alert pairs have "misc_not_failure" here,
  # because if a parent becomes a misc_feature (instead of failing) we still want
  # children to get the corresponding alerts (e.g. peptrans)
  # they may (or may not) escape failure themselves if they have "misc_not_failure"
  # but we want the alerts reported regardless
  my @fatal_alt_codes_A = (); # array of all alert codes with "pertype" eq "feature" and "causes_failure" == 1
  foreach my $alt_code (sort keys (%{$alt_info_HHR})) { 
    if(($alt_info_HHR->{$alt_code}{"pertype"} eq "feature") && 
       ($alt_info_HHR->{$alt_code}{"causes_failure"} == 1)) { 
      push(@fatal_alt_codes_A, $alt_code);
    }
  }

  # for each sequence:
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name = $seq_name_AR->[$seq_idx];
    if(defined $alt_ftr_instances_HHHR->{$seq_name}) { 
      # at least one feature alert exists for this sequence
      # for each feature that is a parent of type $parent_type with
      # at least one child of type $child_type:
      foreach my $parent_ftr_idx (@parent_ftr_idx_A) { 
        if(defined $alt_ftr_instances_HHHR->{$seq_name}{$parent_ftr_idx}) {
          # at least one feature alert exists for this parent feature in this sequence
          # check if any of the alerts for this parent feature are fatal
          my $have_fatal = check_for_feature_alert_codes($alt_info_HHR, \@fatal_alt_codes_A, $alt_ftr_instances_HHHR->{$seq_name}{$parent_ftr_idx});
          if($have_fatal) { 
            # at least one fatal alert for this parent feature, add child alert if it doesn't already exist
            my $nchildren = scalar(@{$children_AA[$parent_ftr_idx]});
            for(my $child_idx = 0; $child_idx < $nchildren; $child_idx++) { 
              my $child_ftr_idx = $children_AA[$parent_ftr_idx][$child_idx];
              if((! defined $alt_ftr_instances_HHHR->{$seq_name}) ||
                 (! defined $alt_ftr_instances_HHHR->{$seq_name}{$child_ftr_idx}) ||
                 (! defined $alt_ftr_instances_HHHR->{$seq_name}{$child_ftr_idx}{$alt_code})) { 
                alert_feature_instance_add($alt_ftr_instances_HHHR, $alt_info_HHR, $alt_code, $seq_name, $child_ftr_idx, $alt_msg, $FH_HR);
              }
            }
          }
        }
      }
    }
  }
  return;
}

#################################################################
# Subroutine: alert_add_unexdivg()
# Incept:     EPN, Thu Feb  7 11:54:56 2019
# Purpose:    Adds unexdivg alerts for sequences listed in the array @overflow_seq_A, if any.
#
# Arguments:
#  $overflow_seq_AR:         REF to array of sequences that failed due to matrix overflows, pre-filled
#  $overflow_mxsize_AR:      REF to array of required matrix sizes for each sequence that failed due to matrix overflows, pre-filled
#  $alt_seq_instances_HHR:   REF to 2D hash with per-sequence alerts, PRE-FILLED
#  $alt_info_HHR:            REF to the alert info hash of arrays, PRE-FILLED
#  $opt_HHR:                 REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:          REF to the 2D hash of output file information
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub alert_add_unexdivg { 
  my $sub_name = "alert_add_unexdivg";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($overflow_seq_AR, $overflow_mxsize_AR, $alt_seq_instances_HHR, $alt_info_HHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience
  
  my $noverflow = scalar(@{$overflow_seq_AR});
  for(my $s = 0; $s < $noverflow; $s++) { 
    alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "unexdivg", $overflow_seq_AR->[$s], "required matrix size: $overflow_mxsize_AR->[$s] Mb", $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: alert_add_ambgnt5s_ambgnt3s()
# Incept:     EPN, Fri Apr 17 10:31:22 2020
# Purpose:    Adds ambgnt5s and ambgnt3s alerts for seqs with 
#             an N as the first/final nucleotide
#
# Arguments:
#  $in_sqfile_R:           REF to Bio::Easel::SqFile object from input fasta file
#  $seq_name_AR:           REF to array of sequence names, PRE-FILLED
#  $seq_len_HR:            REF to hash of sequence lengths, PRE-FILLED
#  $alt_seq_instances_HHR: REF to 2D hash with per-sequence alerts, PRE-FILLED
#  $alt_info_HHR:          REF to the alert info hash of arrays, PRE-FILLED
#  $opt_HHR:               REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR: REF to the 2D hash of output file information
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub alert_add_ambgnt5s_ambgnt3s {
  my $sub_name = "alert_add_ambgnt5s_ambgnt3s";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($in_sqfile_R, $seq_name_AR, $seq_len_HR, $alt_seq_instances_HHR, $alt_info_HHR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience
  my $nseq = scalar(@seq_name_A);

  foreach my $seq_name (@{$seq_name_AR}) { 
    if(! defined $seq_len_HR->{$seq_name}) { 
      ofile_FAIL("ERROR in $sub_name, do not have length info for sequence $seq_name", 1, $FH_HR);
    }
    my $seq_len = $seq_len_HR->{$seq_name};
    my $first_nt = $$in_sqfile_R->fetch_subseq_to_sqstring($seq_name,        1,        1, 0); # 0: do not reverse complement
    my $final_nt = $$in_sqfile_R->fetch_subseq_to_sqstring($seq_name, $seq_len, $seq_len, 0); # 0: do not reverse complement
    my $sqstring = undef;
    my $pos_retval = undef;
    if(($first_nt eq "N") || ($first_nt eq "n")) { 
      # determine first non-N
      $sqstring = $$in_sqfile_R->fetch_seq_to_sqstring($seq_name);  
      my $nlen = count_terminal_ambiguities_in_sqstring($sqstring);
      my $alt_scoords = "seq:" . vdr_CoordsSegmentCreate(1, $nlen, "+", $FH_HR) . ";";
      my $alt_mcoords = "mdl:VADRNULL;"; # this will be updated later in parse_stk_and_add_alignment_cds_and_mp_alerts()
      alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "ambgnt5s", $seq_name, 
                                  sprintf("%s%sVADRNULL", $alt_scoords, $alt_mcoords), $FH_HR); 
    }
    if(($final_nt eq "N") || ($final_nt eq "n")) { 
      if(! defined $sqstring) { $sqstring = $$in_sqfile_R->fetch_seq_to_sqstring($seq_name); }
      my $rev_sqstring = reverse($sqstring);
      my $nlen = count_terminal_ambiguities_in_sqstring($rev_sqstring);
      my $alt_scoords = "seq:" . vdr_CoordsSegmentCreate(($seq_len - $nlen + 1), $seq_len, "+", $FH_HR) . ";";
      my $alt_mcoords = "mdl:VADRNULL;"; # this will be updated later in parse_stk_and_add_alignment_cds_and_mp_alerts()
      alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "ambgnt3s", $seq_name, 
                                  sprintf("%s%sVADRNULL", $alt_scoords, $alt_mcoords), $FH_HR);
    }
  }

  return;
}

#################################################################
# Subroutine:  alert_instances_check_prevents_annot()
# Incept:      EPN, Fri Mar 22 11:57:32 2019
#

# Purpose:    Given a sequence name, determine if any alerts 
#             stored in %{$alt_seq_instances_HHR} exist for it
#             that prevent its annotation.
#
# Arguments: 
#  $seq_name:              name of sequence
#  $alt_info_HHR:          REF to the alert info hash of arrays, PRE-FILLED
#  $alt_seq_instances_HHR: REF to 2D hashes with per-sequence alerts, PRE-FILLED
#  $FH_HR:                 REF to hash of file handles
#
# Returns:    '1' if >=1 alerts prevent annotation
#             '0' if 0   alerts prevent annotation
#
################################################################# 
sub alert_instances_check_prevents_annot { 
  my $sub_name = "alert_instances_check_prevents_annot";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $alt_info_HHR, $alt_seq_instances_HHR, $FH_HR) = @_;

  if(defined $alt_seq_instances_HHR->{$seq_name}) { 
    foreach my $alt_code (keys (%{$alt_seq_instances_HHR->{$seq_name}})) { 
      if(($alt_info_HHR->{$alt_code}{"pertype"} eq "sequence") && 
         ($alt_info_HHR->{$alt_code}{"prevents_annot"} == 1)) { 
        return 1;
      }
    }
  }
  
  return 0;
}

#################################################################
# Subroutine:  alert_feature_instances_count_fatal()
# Incept:      EPN, Wed Oct 13 12:13:15 2021
#

# Purpose:    Given a sequence name and feature index, count the
#             number of fatal alerts stored in %{$alt_ftr_instances_HHHR}
#             and return that number.
#
# Arguments: 
#  $seq_name:               name of sequence
#  $ftr_idx:                feature index
#  $alt_info_HHR:           REF to the alert info hash of arrays, PRE-FILLED
#  $alt_ftr_instances_HHHR: REF to 3D hashes with per-feature alerts, PRE-FILLED
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    Number of fatal alerts for this sequence and feature
#
################################################################# 
sub alert_feature_instances_count_fatal {
  my $sub_name = "alert_feature_instances_count_fatal";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $ftr_idx, $alt_info_HHR, $alt_ftr_instances_HHHR, $FH_HR) = @_;

  my $nfatal = 0;
  if(defined $alt_ftr_instances_HHHR->{$seq_name}) { 
    if(defined $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}) { 
      foreach my $alt_code (keys (%{$alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}})) { 
        if($alt_info_HHR->{$alt_code}{"causes_failure"}) { 
          my @instance_str_A = split(":VADRSEP:", alert_feature_instance_fetch($alt_ftr_instances_HHHR, $seq_name, $ftr_idx, $alt_code));
          $nfatal += scalar(@instance_str_A);
        }
      }
    }
  }
  
  return $nfatal;
}

#################################################################
#
# Subroutines for creating tabular output:
# output_tabular
# helper_tabular_ftr_results_strand
# helper_tabular_ftr_results_trunc_string
# helper_tabular_ftr_results_5N_string
# helper_tabular_sgm_results_trunc_string
# helper_tabular_get_ftr_alert_strings
# helper_tabular_get_seq_alert_strings
# helper_tabular_replace_spaces
# 
#################################################################
# Subroutine: output_tabular()
# Incept:     EPN, Mon Mar  4 21:02:12 2019
# Purpose:    Output tabular files.
#
# Arguments:
#  $mdl_info_AHR:            REF to array of hashes with model info
#  $mdl_cls_ct_HR:           REF to hash with counts of seqs classified per model
#  $mdl_ant_ct_HR:           REF to hash with counts of seqs annotated per model
#  $seq_name_AR:             REF to array of sequence names
#  $seq_len_HR:              REF to hash of sequence lengths
#  $ftr_info_HAHR:           REF to hash of arrays with information on the features, PRE-FILLED
#  $sgm_info_HAHR:           REF to hash of arrays with information on the segments, PRE-FILLED
#  $alt_info_HHR:            REF to the alert info hash of arrays, PRE-FILLED
#  $cls_output_HHR:          REF to 2D hash of classification results to output, PRE-FILLED
#  $ftr_results_HAHR:        REF to feature results AAH, PRE-FILLED
#  $sgm_results_HAHR:        REF to model results AAH, PRE-FILLED
#  $alt_seq_instances_HHR:   REF to 2D hash with per-sequence alerts, PRE-FILLED
#  $alt_ftr_instances_HHHR:  REF to array of 2D hashes with per-feature alerts, PRE-FILLED
#  $dcr_output_HAHR:         REF to hash of array of hashes with info on doctored seqs to output, PRE-FILLED, most seqs will be undef
#  $sda_output_HHR:          REF to 2D hash of -s related results to output, PRE-FILLED, undef unless -s
#  $rpn_output_HHR:          REF to 2D hash of -r related results to output, PRE-FILLED, undef unless -r
#  $mdl_sub_HR:              REF to hash of of model substitutions, PRE-FILLED, undef unless --msub used
#  $opt_HHR:                 REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:          REF to the 2D hash of output file information
#             
# Returns:  Two values:
#           $zero_classifications: '1' if no sequences were classified, else '0'
#           $zero_alerts:          '1' if no alerts were reported, else '0'
# 
# Dies:     never
#
#################################################################
sub output_tabular { 
  my $sub_name = "output_tabular";
  my $nargs_exp = 19;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_info_AHR, $mdl_cls_ct_HR, $mdl_ant_ct_HR, 
      $seq_name_AR, $seq_len_HR, 
      $ftr_info_HAHR, $sgm_info_HAHR, $alt_info_HHR, 
      $cls_output_HHR, $ftr_results_HHAHR, $sgm_results_HHAHR, $alt_seq_instances_HHR, 
      $alt_ftr_instances_HHHR, $dcr_output_HAHR, $sda_output_HHR, $rpn_output_HHR,
      $mdl_sub_HR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience

  # if --glsearch we won't have PP values
  my $do_glsearch = opt_Get("--glsearch", $opt_HHR) ? 1 : 0;

  # deal with --sidx offset
  my $sidx_offset = opt_Get("--sidx", $opt_HHR) - 1;
  my $do_headers  = ($sidx_offset == 0) ? 1 : 0; # if --sidx value is not 1 do not print comment header lines

  # validate input and determine maximum counts of things
  my $nseq = scalar(@{$seq_name_AR});
  my $nalt = scalar(keys %{$alt_info_HHR});
  my $nmdl = scalar(@{$mdl_info_AHR});
  my ($mdl_idx, $ftr_idx, $sgm_idx);
  my $mdl_name;

  # determine order of alert codes to print
  my $alt_code;
  my @alt_code_A     = (); # all alerts in output order
  my @seq_alt_code_A = (); # per-sequence alerts in output order
  my @ftr_alt_code_A = (); # per-feature  alerts in output order
  alert_order_arrays($alt_info_HHR, \@alt_code_A, \@seq_alt_code_A, \@ftr_alt_code_A);

  my $nalt_seq_possible = scalar(@seq_alt_code_A);
  my %alt_ct_H     = (); # total number of alerts per alert code
  my %alt_seq_ct_H = (); # number of seqs with at least one instance of each alert code
  foreach $alt_code (@alt_code_A) { 
    $alt_ct_H{$alt_code} = 0;
    $alt_seq_ct_H{$alt_code} = 0;
  }

  # make hash of model names to indices in mdl_info_AHR, we'll 
  # use this when creating the per-model table
  my %mdl_idx_H = ();
  utl_IdxHFromAH(\%mdl_idx_H, $mdl_info_AHR, "name", $sub_name, $FH_HR);
  my %mdl_pass_ct_H = (); # number of seqs per model that pass
  my %mdl_fail_ct_H = (); # number of seqs per model that fail

  # define the header columns
  my @head_ant_AA = ();
  my @data_ant_AA = ();
  my @clj_ant_A   = ();
  helper_tabular_fill_header_and_justification_arrays("ant", \@head_ant_AA, \@clj_ant_A, $FH_HR);

  my @head_cls_AA = ();
  my @data_cls_AA = ();
  my @clj_cls_A   = ();
  helper_tabular_fill_header_and_justification_arrays("cls", \@head_cls_AA, \@clj_cls_A, $FH_HR);

  my @head_ftr_AA = ();
  my @data_ftr_AA = ();
  my @clj_ftr_A   = ();
  helper_tabular_fill_header_and_justification_arrays("ftr", \@head_ftr_AA, \@clj_ftr_A, $FH_HR);

  my @head_sgm_AA = ();
  my @data_sgm_AA = ();
  my @clj_sgm_A   = ();
  helper_tabular_fill_header_and_justification_arrays("sgm", \@head_sgm_AA, \@clj_sgm_A, $FH_HR);

  my @head_alt_AA = ();
  my @data_alt_AA = ();
  my @clj_alt_A   = ();
  helper_tabular_fill_header_and_justification_arrays("alt", \@head_alt_AA, \@clj_alt_A, $FH_HR);

  my @head_alc_AA = ();
  my @data_alc_AA = ();
  my @clj_alc_A   = ();
  helper_tabular_fill_header_and_justification_arrays("alc", \@head_alc_AA, \@clj_alc_A, $FH_HR);

  my @head_mdl_AA = ();
  my @data_mdl_AA = ();
  my @clj_mdl_A   = ();
  helper_tabular_fill_header_and_justification_arrays("mdl", \@head_mdl_AA, \@clj_mdl_A, $FH_HR);

  my @head_dcr_AA = ();
  my @data_dcr_AA = ();
  my @clj_dcr_A   = ();
  helper_tabular_fill_header_and_justification_arrays("dcr", \@head_dcr_AA, \@clj_dcr_A, $FH_HR);

  # optional .sda file
  my $do_sda = opt_Get("-s", $opt_HHR) ? 1 : 0;
  my @head_sda_AA = ();
  my @data_sda_AA = ();
  my @clj_sda_A   = ();
  helper_tabular_fill_header_and_justification_arrays("sda", \@head_sda_AA, \@clj_sda_A, $FH_HR);

  # optional .rpn file
  my $do_rpn = opt_Get("-r", $opt_HHR) ? 1 : 0;
  my @head_rpn_AA = ();
  my @data_rpn_AA = ();
  my @clj_rpn_A   = ();
  helper_tabular_fill_header_and_justification_arrays("rpn", \@head_rpn_AA, \@clj_rpn_A, $FH_HR);

  my $zero_classifications = 1; # set to '0' below if we have >= 1 seqs that are classified ($seq_mdl1 ne "-")

  # main loop: for each sequence
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name  = $seq_name_AR->[$seq_idx];
    my $seq_len   = $seq_len_HR->{$seq_name};
    my $seq_nftr_alt    = 0;
    my $seq_nseq_alt    = 0;
    my $seq_nftr_annot  = 0;
    my $seq_nftr_5trunc = 0;
    my $seq_nftr_3trunc = 0;
    my $seq_idx2print   = $seq_idx + $sidx_offset + 1;
    my $nftr = 0;
 
   # get per-sequence info from %{$cls_output_HHR->{$seq_name}}
    my $cls_output_HR = (defined $cls_output_HHR->{$seq_name}) ? \%{$cls_output_HHR->{$seq_name}} : undef;
    my $seq_score   = ((defined $cls_output_HR) && (defined $cls_output_HR->{"score"}))     ? $cls_output_HR->{"score"}     : "-";
    my $seq_scpnt   = ((defined $cls_output_HR) && (defined $cls_output_HR->{"scpnt"}))     ? $cls_output_HR->{"scpnt"}     : "-";
    my $seq_scov    = ((defined $cls_output_HR) && (defined $cls_output_HR->{"scov"}))      ? $cls_output_HR->{"scov"}      : "-";
    my $seq_mcov    = ((defined $cls_output_HR) && (defined $cls_output_HR->{"mcov"}))      ? $cls_output_HR->{"mcov"}      : "-";
    my $seq_bias    = ((defined $cls_output_HR) && (defined $cls_output_HR->{"bias"}))      ? $cls_output_HR->{"bias"}      : "-";
    my $seq_nhits   = ((defined $cls_output_HR) && (defined $cls_output_HR->{"nhits"}))     ? $cls_output_HR->{"nhits"}     : "-";
    my $seq_strand  = ((defined $cls_output_HR) && (defined $cls_output_HR->{"bstrand"}))   ? $cls_output_HR->{"bstrand"}   : "_";
    my $seq_scdiff  = ((defined $cls_output_HR) && (defined $cls_output_HR->{"scdiff"}))    ? $cls_output_HR->{"scdiff"}    : "-";
    my $seq_diffpnt = ((defined $cls_output_HR) && (defined $cls_output_HR->{"diffpnt"}))   ? $cls_output_HR->{"diffpnt"}   : "-";
    my $seq_mdl1    = ((defined $cls_output_HR) && (defined $cls_output_HR->{"model1"}))    ? $cls_output_HR->{"model1"}    : "-";
    my $seq_mdl2    = ((defined $cls_output_HR) && (defined $cls_output_HR->{"model2"}))    ? $cls_output_HR->{"model2"}    : "-";
    my $seq_grp1    = ((defined $cls_output_HR) && (defined $cls_output_HR->{"group1"}))    ? $cls_output_HR->{"group1"}    : "-";
    my $seq_grp2    = ((defined $cls_output_HR) && (defined $cls_output_HR->{"group2"}))    ? $cls_output_HR->{"group2"}    : "-";
    my $seq_subgrp1 = ((defined $cls_output_HR) && (defined $cls_output_HR->{"subgroup1"})) ? $cls_output_HR->{"subgroup1"} : "-";
    my $seq_subgrp2 = ((defined $cls_output_HR) && (defined $cls_output_HR->{"subgroup2"})) ? $cls_output_HR->{"subgroup2"} : "-";

    my $seq_mdl_rpn = ((defined $cls_output_HR) && (defined $cls_output_HR->{"rpn.model1"})) ? $cls_output_HR->{"rpn.model1"} : "-";

    my $sda_output_HR = (($do_sda) && (defined $sda_output_HHR->{$seq_name})) ? \%{$sda_output_HHR->{$seq_name}} : undef;
    my $sda_seq       = (($do_sda) && (defined $sda_output_HR->{"sda_seq"}))  ? $sda_output_HR->{"sda_seq"} : "-";
    my $sda_mdl       = (($do_sda) && (defined $sda_output_HR->{"sda_mdl"}))  ? $sda_output_HR->{"sda_mdl"} : "-";
    my $sda_fract     = (($do_sda) && (defined $sda_output_HR->{"sda_seq"}))  ? vdr_CoordsLength($sda_output_HR->{"sda_seq"}, $FH_HR) / $seq_len : "-";
    my $sda_5p_seq    = (($do_sda) && (defined $sda_output_HR->{"5p_seq"}))  ? $sda_output_HR->{"5p_seq"} : "-";
    my $sda_5p_mdl    = (($do_sda) && (defined $sda_output_HR->{"5p_mdl"}))  ? $sda_output_HR->{"5p_mdl"} : "-";
    my $sda_5p_fract  = (($do_sda) && (defined $sda_output_HR->{"5p_seq"}))  ? vdr_CoordsLength($sda_output_HR->{"5p_seq"}, $FH_HR) / $seq_len : "-";
    my $sda_3p_seq    = (($do_sda) && (defined $sda_output_HR->{"3p_seq"}))  ? $sda_output_HR->{"3p_seq"} : "-";
    my $sda_3p_mdl    = (($do_sda) && (defined $sda_output_HR->{"3p_mdl"}))  ? $sda_output_HR->{"3p_mdl"} : "-";
    my $sda_3p_fract  = (($do_sda) && (defined $sda_output_HR->{"3p_seq"}))  ? vdr_CoordsLength($sda_output_HR->{"3p_seq"}, $FH_HR) / $seq_len : "-";
     
    my $rpn_output_HR      = (($do_rpn) && (defined $rpn_output_HHR->{$seq_name}))       ? \%{$rpn_output_HHR->{$seq_name}} : undef;
    my $rpn_nnt_n_tot      = (($do_rpn) && (defined $rpn_output_HR->{"nnt_n_tot"}))      ? $rpn_output_HR->{"nnt_n_tot"}      : "-";
    my $rpn_nnt_n_rp_tot   = (($do_rpn) && (defined $rpn_output_HR->{"nnt_n_rp_tot"}))   ? $rpn_output_HR->{"nnt_n_rp_tot"}   : "-";
    my $rpn_nnt_n_rp_fract = (($do_rpn) && (defined $rpn_output_HR->{"nnt_n_rp_fract"})) ? $rpn_output_HR->{"nnt_n_rp_fract"} : "-";
    my $rpn_ngaps_tot      = (($do_rpn) && (defined $rpn_output_HR->{"ngaps_tot"}))      ? $rpn_output_HR->{"ngaps_tot"}      : "-";
    my $rpn_ngaps_int      = (($do_rpn) && (defined $rpn_output_HR->{"ngaps_int"}))      ? $rpn_output_HR->{"ngaps_int"}      : "-";
    my $rpn_ngaps_rp       = (($do_rpn) && (defined $rpn_output_HR->{"ngaps_rp"}))       ? $rpn_output_HR->{"ngaps_rp"}       : "-";
    my $rpn_ngaps_rp_full  = (($do_rpn) && (defined $rpn_output_HR->{"ngaps_rp_full"}))  ? $rpn_output_HR->{"ngaps_rp_full"}  : "-";
    my $rpn_ngaps_rp_part  = (($do_rpn) && (defined $rpn_output_HR->{"ngaps_rp_part"}))  ? $rpn_output_HR->{"ngaps_rp_part"}  : "-";
    my $rpn_nnt_rp_full    = (($do_rpn) && (defined $rpn_output_HR->{"nnt_rp_full"}))    ? $rpn_output_HR->{"nnt_rp_full"}    : "-";
    my $rpn_nnt_rp_part    = (($do_rpn) && (defined $rpn_output_HR->{"nnt_rp_part"}))    ? $rpn_output_HR->{"nnt_rp_part"}    : "-";
    my $rpn_coords         = (($do_rpn) && (defined $rpn_output_HR->{"coords"}) && ($rpn_output_HR->{"coords"} ne "")) ? $rpn_output_HR->{"coords"} : "-";

    # deal with --msub and model subs, we need this ONLY for mdl_{pass,fail}_ct_H and for ftr_info_HAHR in check_if_sequence_passes
    my $tmp_mdl = $seq_mdl1;
    if(($seq_mdl1 ne "-") && (defined $mdl_sub_HR) && (defined $mdl_sub_HR->{$seq_mdl1})) { 
      $tmp_mdl = $mdl_sub_HR->{$seq_mdl1};
    }

    my $seq_pass_fail = (check_if_sequence_passes($seq_name, (($tmp_mdl ne "-") ? \@{$ftr_info_HAHR->{$tmp_mdl}} : undef), $alt_info_HHR, $alt_seq_instances_HHR, $alt_ftr_instances_HHHR, $FH_HR)) ? "PASS" : "FAIL";
    my $seq_annot     = (check_if_sequence_was_annotated($seq_name, $cls_output_HHR)) ? "yes" : "no";

    if($seq_mdl1 ne "-") { 
      if(! defined $mdl_pass_ct_H{$tmp_mdl}) { $mdl_pass_ct_H{$tmp_mdl} = 0; }
      if(! defined $mdl_fail_ct_H{$tmp_mdl}) { $mdl_fail_ct_H{$tmp_mdl} = 0; }
      if($seq_pass_fail eq "PASS") { $mdl_pass_ct_H{$tmp_mdl}++; }
      if($seq_pass_fail eq "FAIL") { $mdl_fail_ct_H{$tmp_mdl}++; }
      $zero_classifications = 0;
    }

    my $ftr_nprinted   = 0; # number of total features printed for this sequence 
    my $sgm_nprinted   = 0; # number of total segments printed for this sequence
    my $alt_nprinted   = 0; # number of total alerts printed for this sequence
    my $alt_nftr       = 0; # number of features we've printed at least 1 alert for
    my $alt_nseqftr    = 0; # number of features we've printed for this sequence/feature combo
    my %alt_seqcode_H  = (); # key is alert code: '1' if we've printed >= 1 alerts for this sequence/code combo

    # print per-sequence alerts, if any
    if(defined $alt_seq_instances_HHR->{$seq_name}) { 
      foreach $alt_code (@seq_alt_code_A) { 
        my $alt_instance = alert_sequence_instance_fetch($alt_seq_instances_HHR, $seq_name, $alt_code);
        if(defined $alt_instance) { 
          if(($alt_nprinted == 0) && ((scalar(@data_alt_AA) > 0) || (! $do_headers))) { 
            push(@data_alt_AA, []);  # push empty array --> blank line 
            # if (!$do_headers) for --split, we add blank line before first data line to mimic non-split output
          }
          if(! defined $alt_seqcode_H{$alt_code}) { 
            $alt_seq_ct_H{$alt_code}++; 
            $alt_seqcode_H{$alt_code} = 1;
          }
          if($alt_nseqftr == 0) { 
            $alt_nftr++;
          }
          my @instance_str_A = split(":VADRSEP:", $alt_instance);
          foreach my $instance_str (@instance_str_A) { 
            my ($instance_scoords, $instance_mcoords, $instance_detail) = alert_instance_parse($instance_str);
            $alt_nseqftr++;
            $alt_ct_H{$alt_code}++;
            my $alt_idx2print = $seq_idx2print . "." . $alt_nftr . "." . $alt_nseqftr;
            push(@data_alt_AA, [$alt_idx2print, $seq_name, $seq_mdl1, "-", "-", "-", $alt_code, 
                                $alt_info_HHR->{$alt_code}{"causes_failure"} ? "yes" : "no", 
                                helper_tabular_replace_spaces($alt_info_HHR->{$alt_code}{"sdesc"}), 
                                ($instance_scoords eq "VADRNULL") ? "-" : $instance_scoords, 
                                ($instance_scoords eq "VADRNULL") ? "-" : vdr_CoordsLength($instance_scoords, $FH_HR),
                                ($instance_mcoords eq "VADRNULL") ? "-" : $instance_mcoords, 
                                ($instance_mcoords eq "VADRNULL") ? "-" : vdr_CoordsLength($instance_mcoords, $FH_HR),
                                $alt_info_HHR->{$alt_code}{"ldesc"} . (($instance_detail eq "VADRNULL") ? " [-]" : " [" . $instance_detail . "]")]);
            $alt_nprinted++;
          }
        }
      }
    }

    if($seq_mdl1 ne "-") { 
      $nftr = scalar(@{$ftr_info_HAHR->{$seq_mdl1}});
      my $ftr_info_AHR = \@{$ftr_info_HAHR->{$seq_mdl1}}; # for convenience
      for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
        if((defined $ftr_results_HHAHR->{$seq_mdl1}) && 
           (defined $ftr_results_HHAHR->{$seq_mdl1}{$seq_name}) && 
           (defined $ftr_results_HHAHR->{$seq_mdl1}{$seq_name}[$ftr_idx])) { 
          my $ftr_results_HR = $ftr_results_HHAHR->{$seq_mdl1}{$seq_name}[$ftr_idx]; # for convenience
          my $ftr_idx2print = $seq_idx2print . "." . ($seq_nftr_annot + 1);
          if((defined $ftr_results_HR->{"n_start"}) || (defined $ftr_results_HR->{"p_qstart"})) { 
            $seq_nftr_annot++;
            my $ftr_name = $ftr_info_AHR->[$ftr_idx]{"outname"};
            my $ftr_name2print = helper_tabular_replace_spaces($ftr_name);
            my $ftr_parent_idx = ((defined $ftr_info_AHR->[$ftr_idx]{"parent_idx_str"}) && ($ftr_info_AHR->[$ftr_idx]{"parent_idx_str"} ne "GBNULL")) ? ($ftr_info_AHR->[$ftr_idx]{"parent_idx_str"}+1) : -1;
            my $ftr_type = $ftr_info_AHR->[$ftr_idx]{"type"};
            my $ftr_strand   = helper_tabular_ftr_results_strand($ftr_info_AHR, $ftr_results_HR, $ftr_idx);
            my $ftr_trunc    = helper_tabular_ftr_results_trunc_string($ftr_results_HR);
            my $ftr_5ablen   = (defined $ftr_results_HR->{"n_5ablen"})   ? $ftr_results_HR->{"n_5ablen"}   : "-";
            my $ftr_3ablen   = (defined $ftr_results_HR->{"n_3ablen"})   ? $ftr_results_HR->{"n_3ablen"}   : "-";
            my $ftr_n_start  = (defined $ftr_results_HR->{"n_start"})   ? $ftr_results_HR->{"n_start"}   : "-";
            my $ftr_n_stop   = (defined $ftr_results_HR->{"n_stop"})    ? $ftr_results_HR->{"n_stop"}    : "-";
            my $ftr_n_stop_c = (defined $ftr_results_HR->{"n_stop_c"})  ? $ftr_results_HR->{"n_stop_c"}  : "-";
            if(($ftr_n_stop_c ne "-") && ($ftr_n_stop_c ne "?") && ($ftr_n_stop_c == $ftr_n_stop)) { $ftr_n_stop_c = "-"; }
            my $ftr_p_qstart  = (defined $ftr_results_HR->{"p_qstart"})   ? $ftr_results_HR->{"p_qstart"}   : "-";
            my $ftr_p_qstop   = (defined $ftr_results_HR->{"p_qstop"})    ? $ftr_results_HR->{"p_qstop"}    : "-";
            my $ftr_p_qstop_c = (defined $ftr_results_HR->{"p_trcstop"}) ? $ftr_results_HR->{"p_trcstop"} : "-";
            if($ftr_p_qstop_c ne "-") { 
              $ftr_p_qstop_c =~ s/;.*$//; # keep only first early stop position
            }
            my $ftr_p_score = (defined $ftr_results_HR->{"p_score"})  ? $ftr_results_HR->{"p_score"} : "-";
            if(((defined $ftr_results_HR->{"n_5trunc"})  && ($ftr_results_HR->{"n_5trunc"})) ||      # feature is 5' truncated due to sequence end
               ((defined $ftr_results_HR->{"n_5ablen"}) && ($ftr_results_HR->{"n_5ablen"} > 0))) { # feature starts with >= 1 N
              $seq_nftr_5trunc++; 
            }
            if(((defined $ftr_results_HR->{"n_3trunc"})  && ($ftr_results_HR->{"n_3trunc"})) ||      # feature is 3' truncated due to sequence end
               ((defined $ftr_results_HR->{"n_3ablen"}) && ($ftr_results_HR->{"n_3ablen"} > 0))) { # feature ends with >= 1 N
               $seq_nftr_3trunc++; 
            }
            
            my $ftr_alt_str = helper_output_feature_alert_strings($seq_name, $ftr_idx, 1, $alt_info_HHR, \@ftr_alt_code_A, $alt_ftr_instances_HHHR, $FH_HR);
            if($ftr_alt_str ne "") { 
              $seq_nftr_alt++; 
              $seq_nftr_alt += $ftr_alt_str =~ tr/,//; # plus 1 for each comma
            }
            
            my $s_coords_str   = ""; # sequence coordinate string for feature
            my $m_coords_str   = ""; # model    coordinate string for feature
            my $ftr_nsgm_annot = 0;
            my $ftr_len_by_sgm = 0;
            my $ftr_first_sgm  = $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"};
            my $ftr_final_sgm  = $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"};
            my $ftr_nsgm       = $ftr_final_sgm - $ftr_first_sgm + 1;
            for(my $sgm_idx = $ftr_first_sgm; $sgm_idx <= $ftr_final_sgm; $sgm_idx++) { 
              if((defined $sgm_results_HHAHR) && 
                 (defined $sgm_results_HHAHR->{$seq_mdl1}) && 
                 (defined $sgm_results_HHAHR->{$seq_mdl1}{$seq_name}) && 
                 (defined $sgm_results_HHAHR->{$seq_mdl1}{$seq_name}[$sgm_idx]) && 
                 (defined $sgm_results_HHAHR->{$seq_mdl1}{$seq_name}[$sgm_idx]{"sstart"})) { 
                $ftr_nsgm_annot++;
                my $sgm_idx2print = $seq_idx2print . "." . $seq_nftr_annot . "." . $ftr_nsgm_annot;
                my $sgm_results_HR = $sgm_results_HHAHR->{$seq_mdl1}{$seq_name}[$sgm_idx]; # for convenience
                my $sgm_sstart = $sgm_results_HR->{"sstart"};
                my $sgm_sstop  = $sgm_results_HR->{"sstop"};
                my $sgm_mstart = $sgm_results_HR->{"mstart"};
                my $sgm_mstop  = $sgm_results_HR->{"mstop"};
                my $sgm_slen   = abs($sgm_sstart - $sgm_sstop) + 1;
                my $sgm_mlen   = abs($sgm_mstart - $sgm_mstop) + 1;
                my $sgm_strand = $sgm_results_HR->{"strand"};
                my $sgm_trunc  = helper_tabular_sgm_results_trunc_string($sgm_results_HR);
                my ($sgm_pp5, $sgm_pp3);
                if($do_glsearch) { 
                  $sgm_pp5    = ($sgm_results_HR->{"startpp"} ne "?") ? "-" : "?";
                  $sgm_pp3    = ($sgm_results_HR->{"stoppp"}  ne "?") ? "-" : "?";
                }
                else { 
                  $sgm_pp5    = ($sgm_results_HR->{"startpp"} == -1) ? "-" : $sgm_results_HR->{"startpp"};
                  $sgm_pp3    = ($sgm_results_HR->{"stoppp"}  == -1) ? "-" : $sgm_results_HR->{"stoppp"};
                }
                my $sgm_gap5   = ($sgm_results_HR->{"startgap"}) ? "yes" : "no";
                my $sgm_gap3   = ($sgm_results_HR->{"stopgap"})  ? "yes" : "no";
                
                if($s_coords_str ne "") { $s_coords_str .= ","; }
                if($m_coords_str ne "") { $m_coords_str .= ","; }
                $s_coords_str .= $sgm_sstart . ".." . $sgm_sstop . ":" . $sgm_strand;
                $m_coords_str .= $sgm_mstart . ".." . $sgm_mstop . ":" . $sgm_strand;
                $ftr_len_by_sgm += abs($sgm_sstart - $sgm_sstop) + 1;
                
                if(($sgm_nprinted == 0) && ((scalar(@data_sgm_AA) > 0) || (! $do_headers))) { 
                  push(@data_sgm_AA, []); # empty array -> blank line
                  # if (!$do_headers) for --split, we add blank line before first data line to mimic non-split output
                }
                push(@data_sgm_AA, [$sgm_idx2print, $seq_name, $seq_len, $seq_pass_fail, $seq_mdl1, $ftr_type, $ftr_name2print, ($ftr_idx+1), 
                                    $ftr_nsgm, ($sgm_idx-$ftr_first_sgm+1), $sgm_sstart, $sgm_sstop, $sgm_mstart, $sgm_mstop, $sgm_slen, $sgm_strand, 
                                    $sgm_trunc, $sgm_pp5, $sgm_pp3, $sgm_gap5, $sgm_gap3]);
                $sgm_nprinted++;
              }
            } # end of 'for' loop over sgms

            my $ftr_nsgm_noannot = $ftr_nsgm - $ftr_nsgm_annot;
            if($ftr_len_by_sgm == 0) { $ftr_len_by_sgm = "-"; }
            if($ftr_alt_str eq "")   { $ftr_alt_str = "-"; }

            if(($ftr_nprinted == 0) && ((scalar(@data_ftr_AA) > 0) || (! $do_headers))) { 
              push(@data_ftr_AA, []); # empty array -> blank line
              # if (!$do_headers) for --split, we add blank line before first data line to mimic non-split output
            } 
            if($s_coords_str eq "") { $s_coords_str = "-"; } # will happen only for protein-validation only predictions
            if($m_coords_str eq "") { $m_coords_str = "-"; } # will happen only for protein-validation only predictions
            push(@data_ftr_AA, [$ftr_idx2print, $seq_name, $seq_len, $seq_pass_fail, $seq_mdl1, $ftr_type, $ftr_name2print, $ftr_len_by_sgm, 
                                ($ftr_idx+1), $ftr_parent_idx, $ftr_strand, $ftr_n_start, $ftr_n_stop, $ftr_n_stop_c, $ftr_trunc, $ftr_5ablen, $ftr_3ablen, 
                                $ftr_p_qstart, $ftr_p_qstop, $ftr_p_qstop_c, $ftr_p_score, $ftr_nsgm_annot, $ftr_nsgm_noannot, 
                                $s_coords_str, $m_coords_str, $ftr_alt_str]);
            $ftr_nprinted++;

            # print per-feature alerts, if any
            $alt_nseqftr = 0;
            if((defined $alt_ftr_instances_HHHR->{$seq_name}) && 
               (defined $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx})) { 
              foreach my $alt_code (@ftr_alt_code_A) { 
                my $alt_instance = alert_feature_instance_fetch($alt_ftr_instances_HHHR, $seq_name, $ftr_idx, $alt_code);
                if(defined $alt_instance) { 
                  if(($alt_nprinted == 0) && ((scalar(@data_alt_AA) > 0) || (! $do_headers))) { 
                    push(@data_alt_AA, []); # empty array -> blank line
                    # if (!$do_headers) for --split, we add blank line before first data line to mimic non-split output
                  }
                  if(! defined $alt_seqcode_H{$alt_code}) { 
                    $alt_seq_ct_H{$alt_code}++; 
                    $alt_seqcode_H{$alt_code} = 1;
                  }
                  if($alt_nseqftr == 0) { 
                    $alt_nftr++;
                  }
                  my @instance_str_A = split(":VADRSEP:", $alt_instance);
                  foreach my $instance_str (@instance_str_A) { 
                    my ($instance_scoords, $instance_mcoords, $instance_detail) = alert_instance_parse($instance_str);
                    $alt_nseqftr++;
                    $alt_ct_H{$alt_code}++;
                    my $alt_idx2print = $seq_idx2print . "." . $alt_nftr . "." . $alt_nseqftr;
                    push(@data_alt_AA, [$alt_idx2print, $seq_name, $seq_mdl1, $ftr_type, $ftr_name2print, ($ftr_idx+1), $alt_code, 
                                        vdr_FeatureAlertCausesFailure($ftr_info_AHR, $alt_info_HHR, $ftr_idx, $alt_code) ? "yes" : "no", 
                                        helper_tabular_replace_spaces($alt_info_HHR->{$alt_code}{"sdesc"}), 
                                        ($instance_scoords eq "VADRNULL") ? "-" : $instance_scoords, 
                                        ($instance_scoords eq "VADRNULL") ? "-" : vdr_CoordsLength($instance_scoords, $FH_HR),
                                        ($instance_mcoords eq "VADRNULL") ? "-" : $instance_mcoords, 
                                        ($instance_mcoords eq "VADRNULL") ? "-" : vdr_CoordsLength($instance_mcoords, $FH_HR),
                                        $alt_info_HHR->{$alt_code}{"ldesc"} . (($instance_detail eq "VADRNULL") ? " [-]" : " [" . $instance_detail . "]")]);
                    $alt_nprinted++;
                  }
                }
              }
            }
          }
        }
      }
    }

    my $seq_alt_str = helper_output_sequence_alert_strings($seq_name, 1, $alt_info_HHR, \@seq_alt_code_A, $alt_seq_instances_HHR, $FH_HR);
    if($seq_alt_str ne "") { 
      $seq_nseq_alt = 1;
      $seq_nseq_alt += $seq_alt_str =~ tr/,//; # plus 1 for each comma
    }
    my $seq_nftr_notannot = $nftr - $seq_nftr_annot;
    if($seq_alt_str eq "")   { $seq_alt_str  = "-"; }
    if($seq_annot   eq "no") { $seq_nftr_annot = $seq_nftr_notannot = $seq_nftr_5trunc = $seq_nftr_3trunc = $seq_nftr_alt = "-"; }

    push(@data_ant_AA, [$seq_idx2print, $seq_name, $seq_len, $seq_pass_fail, $seq_annot, $seq_mdl1, $seq_grp1, $seq_subgrp1, 
                        $seq_nftr_annot, $seq_nftr_notannot, $seq_nftr_5trunc, $seq_nftr_3trunc, $seq_nftr_alt, $seq_alt_str]);
    
    push(@data_cls_AA, [$seq_idx2print, $seq_name, $seq_len, $seq_pass_fail, $seq_annot, $seq_mdl1, 
                            helper_tabular_replace_spaces($seq_grp1), 
                            helper_tabular_replace_spaces($seq_subgrp1), 
                            $seq_score, $seq_scpnt, $seq_scov, $seq_mcov, $seq_bias, $seq_nhits, $seq_strand, $seq_mdl2, 
                            helper_tabular_replace_spaces($seq_grp2), 
                            helper_tabular_replace_spaces($seq_subgrp2), 
                            $seq_scdiff, $seq_diffpnt, $seq_alt_str]);

    if(defined $dcr_output_HAHR->{$seq_name}) { 
      my $ndcr = scalar(@{$dcr_output_HAHR->{$seq_name}});
      for(my $dcr_idx = 0; $dcr_idx < $ndcr; $dcr_idx++) { 
        my $dcr_idx2print = sprintf("%d.%d", $seq_idx2print, ($dcr_idx+1));
        my $dcr_mdl_name = $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"mdl_name"};
        my $dcr_ftr_idx  = $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"ftr_idx"};
        my $dcr_ftr_name = $ftr_info_HAHR->{$dcr_mdl_name}[$dcr_ftr_idx]{"outname"};
        my $dcr_ftr_name2print = helper_tabular_replace_spaces($dcr_ftr_name);
        push(@data_dcr_AA, [($dcr_idx2print, $seq_name, $dcr_mdl_name,
                             $ftr_info_HAHR->{$dcr_mdl_name}[$dcr_ftr_idx]{"type"}, 
                             $dcr_ftr_name2print, $dcr_ftr_idx, 
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"dcr_type"},
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"rfpos"},
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"indel_apos"}, 
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"orig_seq_uapos"}, 
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"new_seq_uapos"}, 
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"codon_type"}, 
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"codon_coords"}, 
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"orig_codon"}, 
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"new_codon"}, 
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"dcr_iter"}, 
                             $dcr_output_HAHR->{$seq_name}[$dcr_idx]{"did_swap"})]);
      }
      push(@data_dcr_AA, []); # empty array -> blank line
    }

    if($do_sda) {
      my $sda_fract2print = ($sda_fract ne "-") ? sprintf("%.3f", $sda_fract) : "-";
      my $sda_5p_fract2print  = ($sda_5p_fract  ne "-") ? sprintf("%.3f", $sda_5p_fract)  : "-";
      my $sda_3p_fract2print  = ($sda_3p_fract  ne "-") ? sprintf("%.3f", $sda_3p_fract)  : "-";
      push(@data_sda_AA, [$seq_idx2print, $seq_name, $seq_len, $seq_mdl1, $seq_pass_fail,
                          $sda_seq, $sda_mdl, $sda_fract2print, 
                          $sda_5p_seq, $sda_5p_mdl, $sda_5p_fract2print, 
                          $sda_3p_seq, $sda_3p_mdl, $sda_3p_fract2print]);
    }
    if($do_rpn) {
      my $rpn_nnt_n_rp_fract2print = (($rpn_nnt_n_rp_fract ne "-") && ($rpn_nnt_n_tot ne "-") && ($rpn_nnt_n_tot > 0)) ? 
          sprintf("%.3f", $rpn_nnt_n_rp_fract) : "-";
      push(@data_rpn_AA, [$seq_idx2print, $seq_name, $seq_len, $seq_mdl_rpn, $seq_pass_fail,
                          $rpn_nnt_n_tot, $rpn_nnt_n_rp_tot, $rpn_nnt_n_rp_fract2print,
                          $rpn_ngaps_tot, $rpn_ngaps_int, $rpn_ngaps_rp, 
                          $rpn_ngaps_rp_full, $rpn_ngaps_rp_part,
                          $rpn_nnt_rp_full, $rpn_nnt_rp_part,
                          $rpn_coords]);
    }
  }

  # add data to the alert count table
  my $alt_idx = 0;
  my $zero_alerts = 1; # set to '0' below if we have >= 1 alerts
  my $sum_alt_ct     = 0;
  my $sum_alt_seq_ct = 0;
  my $alc_sep_flag = 0;
  foreach my $alt_code (@alt_code_A) { 
    if($alt_ct_H{$alt_code} > 0) { 
      if(! $alt_info_HHR->{$alt_code}{"causes_failure"}) { 
        $alc_sep_flag = 1; 
      }
      if(($alt_info_HHR->{$alt_code}{"causes_failure"}) && $alc_sep_flag) { 
        # print separation line between alerts that cause and do not cause failure
        push(@data_alc_AA, []); # separator line
        $alc_sep_flag = 0; 
      }
      $alt_idx++;
      # determine if this alert has "misc_not_failure" for at least 1 feature in 1 model,
      # if so, we'll add a '*' 
      my $misc_not_failure_flag = 0;
      my $cur_nmdl = scalar(@{$mdl_info_AHR});
      for(my $cur_mdl_idx = 0; $cur_mdl_idx < $cur_nmdl; $cur_mdl_idx++) { 
        my $cur_mdl_name = $mdl_info_AHR->[$cur_mdl_idx]{"name"};
        my $cur_nftr = scalar(@{$ftr_info_HAHR->{$cur_mdl_name}});
        for(my $cur_ftr_idx = 0; $cur_ftr_idx < $cur_nftr; $cur_ftr_idx++) { 
          if(vdr_FeatureAlertIsMiscNotFailure($ftr_info_HAHR->{$cur_mdl_name}, $alt_info_HHR, $cur_ftr_idx, $alt_code)) { 
            $misc_not_failure_flag = 1;
            $cur_ftr_idx = $cur_nftr; # breaks ftr loop
            $cur_mdl_idx = $cur_nmdl; # breaks mdl loop
          }
        }
      }
      my $causes_failure_str = $alt_info_HH{$alt_code}{"causes_failure"} ? "yes" : "no";
      if(($causes_failure_str eq "yes") && ($misc_not_failure_flag)) { 
        $causes_failure_str .= "*";
      }
      push(@data_alc_AA, [$alt_idx, $alt_code, 
                          $causes_failure_str,
                          helper_tabular_replace_spaces($alt_info_HH{$alt_code}{"sdesc"}), 
                          $alt_info_HH{$alt_code}{"pertype"}, 
                          $alt_ct_H{$alt_code}, $alt_seq_ct_H{$alt_code},
                          $alt_info_HHR->{$alt_code}{"ldesc"}]);
      $sum_alt_ct     += $alt_ct_H{$alt_code};
      $sum_alt_seq_ct += $alt_seq_ct_H{$alt_code};
      $zero_alerts = 0;
    }
  }
  if(! $zero_alerts) { 
    push(@data_alc_AA, []); # separator line
  }

  # add data to the model table
  my @mdl_tbl_order_A = (sort { $mdl_cls_ct_HR->{$b} <=> $mdl_cls_ct_HR->{$a} or 
                                    $a cmp $b 
                         } keys (%{$mdl_cls_ct_HR}));
  my $mdl_tbl_idx = 0;
  my $sum_mdl_cls_ct     = 0;
  my $sum_mdl_pass_ct    = 0;
  my $sum_mdl_fail_ct    = 0;
  my $sum_mdl_noannot_ct = 0;
  foreach $mdl_name (@mdl_tbl_order_A) { 
    if($mdl_cls_ct_HR->{$mdl_name} > 0) { 
      $mdl_tbl_idx++;
      $mdl_idx = $mdl_idx_H{$mdl_name};
      my $mdl_ant_ct = (defined $mdl_ant_ct_HR->{$mdl_name}) ? $mdl_ant_ct_HR->{$mdl_name} : 0;
      push(@data_mdl_AA, [$mdl_tbl_idx, $mdl_name, 
                          (defined $mdl_info_AHR->[$mdl_idx]{"group"})    ? $mdl_info_AHR->[$mdl_idx]{"group"}    : "-", 
                          (defined $mdl_info_AHR->[$mdl_idx]{"subgroup"}) ? $mdl_info_AHR->[$mdl_idx]{"subgroup"} : "-", 
                          $mdl_cls_ct_HR->{$mdl_name},
                          (defined $mdl_pass_ct_H{$mdl_name}) ? $mdl_pass_ct_H{$mdl_name} : 0, 
                          (defined $mdl_fail_ct_H{$mdl_name}) ? $mdl_fail_ct_H{$mdl_name} : 0]); 
      $sum_mdl_cls_ct     += $mdl_cls_ct_HR->{$mdl_name};
      $sum_mdl_pass_ct    += (defined $mdl_pass_ct_H{$mdl_name}) ? $mdl_pass_ct_H{$mdl_name} : 0;
      $sum_mdl_fail_ct    += (defined $mdl_fail_ct_H{$mdl_name}) ? $mdl_fail_ct_H{$mdl_name} : 0;
      $sum_mdl_noannot_ct += $mdl_cls_ct_HR->{$mdl_name} - $mdl_ant_ct;
    }
  }
  # add mdl summary line
  push(@data_mdl_AA, []); # separator line
  push(@data_mdl_AA, ["-", "*all*", "-", "-", 
                      $sum_mdl_cls_ct, 
                      $sum_mdl_pass_ct, 
                      $sum_mdl_fail_ct]);
  # add no-model line (sequences that were not classified)
  push(@data_mdl_AA, ["-", "*none*", "-", "-", 
                      $nseq - $sum_mdl_cls_ct, 
                      0, 
                      $nseq - $sum_mdl_cls_ct]);
  push(@data_mdl_AA, []); # separator line

  # output the tables:
  # if we are doing headers, we always call ofile_TableHumanOutput()
  if($do_headers) { 
    ofile_TableHumanOutput(\@data_ant_AA, \@head_ant_AA, \@clj_ant_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"ant"}, undef, $FH_HR);
    ofile_TableHumanOutput(\@data_cls_AA, \@head_cls_AA, \@clj_cls_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"cls"}, undef, $FH_HR);
    ofile_TableHumanOutput(\@data_ftr_AA, \@head_ftr_AA, \@clj_ftr_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"ftr"}, undef, $FH_HR);
    ofile_TableHumanOutput(\@data_sgm_AA, \@head_sgm_AA, \@clj_sgm_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"sgm"}, undef, $FH_HR);
    ofile_TableHumanOutput(\@data_alt_AA, \@head_alt_AA, \@clj_alt_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"alt"}, undef, $FH_HR);
    ofile_TableHumanOutput(\@data_alc_AA, \@head_alc_AA, \@clj_alc_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"alc"}, undef, $FH_HR);
    ofile_TableHumanOutput(\@data_mdl_AA, \@head_mdl_AA, \@clj_mdl_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"mdl"}, undef, $FH_HR);
    ofile_TableHumanOutput(\@data_dcr_AA, \@head_dcr_AA, \@clj_dcr_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"dcr"}, undef, $FH_HR);
    if($do_sda) {
      ofile_TableHumanOutput(\@data_sda_AA, \@head_sda_AA, \@clj_sda_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"sda"}, undef, $FH_HR);
    }
    if($do_rpn) {
      ofile_TableHumanOutput(\@data_rpn_AA, \@head_rpn_AA, \@clj_rpn_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"rpn"}, undef, $FH_HR);
    }
  }
  else { 
    # if we are not doing headers, we only call ofile_TableHumanOutput() if we have data, because ofile_TableHumanOutput() requires at least one of header and data arrays be non-empty
    if(scalar(@data_ant_AA) > 0) { ofile_TableHumanOutput(\@data_ant_AA, undef, \@clj_ant_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"ant"}, undef, $FH_HR); }
    if(scalar(@data_cls_AA) > 0) { ofile_TableHumanOutput(\@data_cls_AA, undef, \@clj_cls_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"cls"}, undef, $FH_HR); }
    if(scalar(@data_ftr_AA) > 0) { ofile_TableHumanOutput(\@data_ftr_AA, undef, \@clj_ftr_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"ftr"}, undef, $FH_HR); }
    if(scalar(@data_sgm_AA) > 0) { ofile_TableHumanOutput(\@data_sgm_AA, undef, \@clj_sgm_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"sgm"}, undef, $FH_HR); }
    if(scalar(@data_alt_AA) > 0) { ofile_TableHumanOutput(\@data_alt_AA, undef, \@clj_alt_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"alt"}, undef, $FH_HR); }
    if(scalar(@data_alc_AA) > 0) { ofile_TableHumanOutput(\@data_alc_AA, undef, \@clj_alc_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"alc"}, undef, $FH_HR); }
    if(scalar(@data_mdl_AA) > 0) { ofile_TableHumanOutput(\@data_mdl_AA, undef, \@clj_mdl_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"mdl"}, undef, $FH_HR); }
    if(scalar(@data_dcr_AA) > 0) { ofile_TableHumanOutput(\@data_dcr_AA, undef, \@clj_dcr_A, undef, undef, "  ", "-", "#", "#", "", 0, $FH_HR->{"dcr"}, undef, $FH_HR); }
    if(($do_sda) && (scalar(@data_sda_AA) > 0)) {
      ofile_TableHumanOutput(\@data_sda_AA, undef, \@clj_sda_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"sda"}, undef, $FH_HR);
    }
    if(($do_rpn) && (scalar(@data_rpn_AA) > 0)) {
      ofile_TableHumanOutput(\@data_rpn_AA, undef, \@clj_rpn_A, undef, undef, "  ", "-", "#", "#", "", 1, $FH_HR->{"rpn"}, undef, $FH_HR);
    }
  }
  return ($zero_classifications, $zero_alerts);
}
    
#################################################################
# Subroutine : helper_tabular_ftr_results_strand()
# Incept:      EPN, Tue Mar  5 09:10:36 2019
#
# Purpose:    Return string describing feature strand based on 
#             $ftr_results_HR->{"n_strand"} && $ftr_results_HR->{"p_strand"}
#
# Arguments: 
#  $ftr_info_AHR:    REF to array of hashes with information on the features, PRE-FILLED
#  $ftr_results_HR:  REF to hash of feature results for one sequence/feature pair
#  $ftr_idx:         feature index
# 
# Returns:    "+", "-" or "?" (if n_strand && p_strand disagree or neither exists)
#
# Dies:       never
#
################################################################# 
sub helper_tabular_ftr_results_strand { 
  my $sub_name = "helper_tabular_ftr_results_strand()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $ftr_results_HR, $ftr_idx) = (@_);

  if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx)) { 
    if((! defined $ftr_results_HR->{"n_strand"}) && 
       (! defined $ftr_results_HR->{"p_strand"})) { 
      # neither defined
      return "?";
    }
    if((defined $ftr_results_HR->{"n_strand"}) && 
       (defined $ftr_results_HR->{"p_strand"})) { 
      if($ftr_results_HR->{"n_strand"} ne $ftr_results_HR->{"p_strand"}) { 
        # both defined, but they disagree
        return "?";
      }
      # both defined and both agree
      return $ftr_results_HR->{"n_strand"}; 
    }
    if(defined $ftr_results_HR->{"n_strand"}) { 
      # only n_strand defined
      return $ftr_results_HR->{"n_strand"}; 
    }
    # only p_strand defined
    return $ftr_results_HR->{"p_strand"}; 
  } # end of 'if(vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx))'
  else { 
    return (defined $ftr_results_HR->{"n_strand"}) ? $ftr_results_HR->{"n_strand"} : "?";
  }
}

#################################################################
# Subroutine : helper_tabular_ftr_results_trunc_string()
# Incept:      EPN, Tue Mar  5 09:10:36 2019
#
# Purpose:    Return string describing truncation mode based on 
#             $ftr_results_HR->{"n_5trunc"} && $ftr_results_HR->{"n_3trunc"}
#
# Arguments: 
#  $ftr_results_HR:   REF to hash of feature results for one sequence/feature pair
#
# Returns:    "5'&3'", "5'", "3'", "no", or "?"
#
# Dies:       never
#
################################################################# 
sub helper_tabular_ftr_results_trunc_string { 
  my $sub_name = "helper_tabular_ftr_results_trunc_string()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_results_HR) = (@_);

  if((! defined $ftr_results_HR->{"n_5trunc"}) || 
     (! defined $ftr_results_HR->{"n_3trunc"})) { 
    return "?";
  }
  if(($ftr_results_HR->{"n_5trunc"}) && 
     ($ftr_results_HR->{"n_3trunc"})) { 
    return "5'&3'";
  }
  elsif($ftr_results_HR->{"n_5trunc"}) { 
    return "5'";
  }
  elsif($ftr_results_HR->{"n_3trunc"}) { 
    return "3'";
  }
  return "no";
}

#################################################################
# Subroutine : helper_tabular_sgm_results_trunc_string()
# Incept:      EPN, Tue Mar  5 15:02:34 2019
#
# Purpose:    Return string describing truncation mode based on 
#             $sgm_results_HR->{"5trunc"} && $sgm_results_HR->{"3trunc"}
#
# Arguments: 
#  $sgm_results_HR:   REF to hash of model results for one sequence/model pair
#
# Returns:    "5'&3'", "5'", "3'", "no", or "?"
#
# Dies:       never
#
################################################################# 
sub helper_tabular_sgm_results_trunc_string { 
  my $sub_name = "helper_tabular_sgm_results_trunc_string()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sgm_results_HR) = (@_);

  if((! defined $sgm_results_HR->{"5trunc"}) || 
     (! defined $sgm_results_HR->{"3trunc"})) { 
    return "?";
  }
  if(($sgm_results_HR->{"5trunc"}) && 
     ($sgm_results_HR->{"3trunc"})) { 
    return "5'&3'";
  }
  elsif($sgm_results_HR->{"5trunc"}) { 
    return "5'";
  }
  elsif($sgm_results_HR->{"3trunc"}) { 
    return "3'";
  }
  return "no";
}

#################################################################
# Subroutine:  helper_tabular_replace_spaces
# Incept:      EPN, Fri Mar 29 13:57:35 2019
#
# Purpose:    Given a string, replace any spaces in it with "_".
#
# Arguments: 
#  $str:       string to remove spaces from
#
# Returns:    $str with spaces converted to "_"
#
################################################################# 
sub helper_tabular_replace_spaces { 
  my $sub_name = "helper_tabular_replace_spaces";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($str) = @_;

  $str =~ s/\s/\_/g;  # replace each space with "_"

  return $str;
}

#################################################################
#
# Subroutines for creating feature table output:
# output_feature_table
# output_parent_child_relationships 
# helper_ftable_coords_from_nt_prediction 
# helper_ftable_coords_prot_only_prediction 
# helper_ftable_start_stop_strand_arrays_to_coords 
# helper_ftable_coords_to_out_str 
# helper_ftable_add_qualifier_from_ftr_info
# helper_ftable_add_qualifier_from_ftr_results
# helper_ftable_class_model_for_sequence
# helper_ftable_process_feature_alerts
# helper_ftable_process_sequence_alerts
#
#################################################################

#################################################################
# Subroutine: output_feature_table()
# Incept:     EPN, Tue Dec  5 13:49:17 2017 [rewritten Tue Oct 30 05:59:04 2018]
# Purpose:    Output the feature table for all sequences.
#
# Arguments:
#  $mdl_cls_ct_HR:           REF to hash with counts of seqs classified per model
#  $seq_name_AR:             REF to hash of arrays with information on the sequences, PRE-FILLED
#  $ftr_info_HAHR:           REF to hash of arrays with information on the features, PRE-FILLED
#  $sgm_info_HAHR:           REF to hash of arrays with information on the segments, PRE-FILLED
#  $alt_info_HHR:            REF to the alert info hash of arrays, PRE-FILLED
#  $stg_results_HHHR:        REF to 3D hash of classification results, PRE-FILLED
#  $ftr_results_HAHR:        REF to feature results AAH, PRE-FILLED
#  $sgm_results_HAHR:        REF to model results AAH, PRE-FILLED
#  $alt_seq_instances_HHR:   REF to 2D hash with per-sequence alerts, PRE-FILLED
#  $alt_ftr_instances_HHHR:  REF to array of 2D hashes with per-feature alerts, PRE-FILLED
#  $mdl_sub_HR:              REF to hash of of model substitutions, PRE-FILLED, should be undef unless --msub used
#  $in_sqfile_R:             REF to Bio::Easel::SqFile object of input fasta, to create .pass.fa and .fail.fa files with
#  $out_root:                output root for the output fasta file names
#  $opt_HHR:                 REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:          REF to the 2D hash of output file information
#             
# Returns:  Number of sequences that 'pass'.
# 
# Dies:     never
#
#################################################################
sub output_feature_table { 
  my $sub_name = "output_feature_table";
  my $nargs_exp = 15;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($mdl_cls_ct_HR, $seq_name_AR, $ftr_info_HAHR, $sgm_info_HAHR, $alt_info_HHR, 
      $stg_results_HHHR, $ftr_results_HHAHR, $sgm_results_HHAHR, $alt_seq_instances_HHR, 
      $alt_ftr_instances_HHHR, $mdl_sub_HR, $in_sqfile_R, $out_root, $opt_HHR, $ofile_info_HHR) = @_;

  my $do_pv_blastx  = (opt_Get("--pv_skip", $opt_HHR) || opt_Get("--pv_hmmer", $opt_HHR)) ? 0 : 1;
  my $do_nofasta    = opt_Get("--out_nofasta", $opt_HHR) ? 1 : 0;
  my $sidx_offset   = opt_Get("--sidx", $opt_HHR) - 1;
  my $do_headers    = ($sidx_offset == 0) ? 1 : 0; # if --sidx value is not 1 do not print comment header lines

  my $FH_HR = $ofile_info_HHR->{"FH"}; # for convenience
  my $pass_ftbl_FH = $FH_HR->{"pass_tbl"};     # feature table for PASSing sequences
  my $fail_ftbl_FH = $FH_HR->{"fail_tbl"};     # feature table for FAILing sequences
  my $pass_list_FH = $FH_HR->{"pass_list"};    # list of PASSing seqs
  my $fail_list_FH = $FH_HR->{"fail_list"};    # list of FAILing seqs
  my $alerts_FH    = $FH_HR->{"alerts_list"};  # list of alerts
  if($do_headers) { 
    print $alerts_FH "#sequence\tmodel\tfeature-type\tfeature-name\terror\tseq-coords\tmdl-coords\terror-description\n";
  }

  my $ret_npass = 0;  # number of sequences that pass, returned from this subroutine
  my $mdl_name = undef;
  my $ftr_idx = undef;

  my $nseq = scalar(@{$seq_name_AR}); # nseq: number of sequences
  my $nalt = scalar(keys %{$alt_info_HHR});

  my $do_nomisc       = opt_Get("--nomisc",       $opt_HHR); # 1 to never output misc_features
  my $do_noprotid     = opt_Get("--noprotid",     $opt_HHR); # 1 to never output protein_id qualifiers
  my $do_forceprotid  = opt_Get("--forceprotid",  $opt_HHR); # 1 to never modify sequence name for protein_id qualifiers
  my $do_noseqnamemax = opt_Get("--noseqnamemax", $opt_HHR); # 1 to allow protein_id values of any length
  my $nmiscftr_thr    = opt_Get("--nmiscftrthr",  $opt_HHR); # max num misc_features allowed via misc_not_failure w/o nmiscftr alert
  my $max_protein_id_length = 50; # hard-coded

  my $do_notrim   = opt_Get("--notrim",   $opt_HHR); # 1 to never trim any features
  my %noftrtrim_H = (); # key is feature type read from --noftrtrim <s> option, value is 1 to not trim start/end due to ambiguities
  if(opt_IsUsed("--noftrtrim", $opt_HHR)) { 
    my @noftrtrim_A  = split(",", opt_Get("--noftrtrim", $opt_HHR));
    foreach my $ftr_type (@noftrtrim_A) { $noftrtrim_H{$ftr_type} = 1; }
  }
  # may want to add 5'UTR and 3'UTR to %noftrtrim_H by default in future
  # don't add mat_peptide to this, even though we don't trim mat_peptide
  # coords, that happens because we trim the parent ftr for child features
  # like mat_peptides, search for 'my $trim_idx' below for details

  # determine order of alert codes to print
  my $alt_code;
  my @seq_alt_code_A = (); # per-sequence alerts in output order
  my @ftr_alt_code_A = (); # per-feature  alerts in output order
  alert_order_arrays($alt_info_HHR, undef, \@seq_alt_code_A, \@ftr_alt_code_A);

  # define the hard-coded type priority hash, which defines the order of types in the feature table output, lower is higher priority
  my %type_priority_H = ();
  $type_priority_H{"gene"}         = 0;
  $type_priority_H{"CDS"}          = 1;
  $type_priority_H{"misc_feature"} = 2;
  $type_priority_H{"mat_peptide"}  = 3;
  my $npriority = scalar(keys %type_priority_H);
  # reference for the way we sort the information we collect for the feature table
  #https://stackoverflow.com/questions/10395383/sorting-an-array-of-hash-by-multiple-keys-perl      

  my $qval_sep = ":GPSEP:"; # value separating multiple qualifier values in a single element of $ftr_info_HAHR->{$mdl_name}[$ftr_idx]{}
  # NOTE: $qval_sep == ':GPSEP:' is hard-coded value for separating multiple qualifier values for the same 
  # qualifier (see vadr.pm::vdr_GenBankStoreQualifierValue)

  # two hash of arrays 1D keys: model names, values are arrays
  # we only fill these for each model as we need it, so as not 
  # to wastefully fill these for models for which no seqs have been assigned
  my %ftr_min_len_HA     = (); # hash of arrays with minimum valid length per model/feature, 1D keys are model names, 2D elements are feature indices

  my @pass_fa_A = (); # array of seq names for pass fasta file
  my @fail_fa_A = (); # array of seq names for fail fasta file

  # main loop: for each sequence
  for(my $seq_idx = 0; $seq_idx < $nseq; $seq_idx++) { 
    my $seq_name = $seq_name_AR->[$seq_idx];
    my $seq_ntabftr = 0; # number of features for this sequence annotated in tabular .ftr file (may have shorter features than are permitted in .ftbl file)
    my $seq_idx2print = $seq_idx + $sidx_offset + 1; # will be $seq_idx unless --sidx used and set to > 1
    
    my @ftout_AH      = (); # array of hashes with output for feature table, kept in a hash so we can sort before outputting
    my $ftidx         = 0;  # index in @ftout_AH
    my %ftr_idx2ftout_idx_H = (); # key is feature index $fidx, value is $ftidx index in @ftout_AH that $fidx corresponds to
    my $i;

    my @seq_alert_A = (); # all alerts for this sequence
    my @seq_note_A  = (); # all notes for this sequence
    my $cur_nmiscftr = 0; # number of features turned into misc_features due to misc_not_failure attributes

    # first check for per-sequence alerts
    my $seq_alt_str = helper_output_sequence_alert_strings($seq_name, 0, $alt_info_HHR, \@seq_alt_code_A, $alt_seq_instances_HHR, $FH_HR);
    my $prevents_annot_flag = helper_ftable_process_sequence_alerts($seq_alt_str, $seq_name, $alt_info_HHR, $alt_seq_instances_HHR, \@seq_alert_A, $FH_HR);

    $mdl_name = helper_ftable_class_model_for_sequence($stg_results_HHHR, $seq_name);
    if(defined $mdl_name) { 
      if((defined $mdl_sub_HR) && (defined $mdl_sub_HR->{$mdl_name})) { 
        $mdl_name = $mdl_sub_HR->{$mdl_name};
      }
      my $ftr_info_AHR     = \@{$ftr_info_HAHR->{$mdl_name}};     # for convenience
      my $ftr_results_HAHR = \%{$ftr_results_HHAHR->{$mdl_name}}; # for convenience
      my $sgm_results_HAHR = \%{$sgm_results_HHAHR->{$mdl_name}}; # for convenience
      my $nftr = scalar(@{$ftr_info_AHR});
      # variables related to protein_id qualifiers for CDS and mat_peptides
      my $nprotein_id = 0; # index of protein_id qualifier, incremented as they are added
      my %ftr_idx2protein_id_idx_H = (); # key is a feature index that is a CDS, value is protein_id index for that feature

      # fill @{$ftr_min_len_HA{$mdl_name}} for this model, if it's not already filled
      if(! defined $ftr_min_len_HA{$mdl_name}) { 
        @{$ftr_min_len_HA{$mdl_name}} = ();
        for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
          my $parent_ftr_idx = vdr_FeatureParentIndex($ftr_info_AHR, $ftr_idx);
          $ftr_min_len_HA{$mdl_name}[$ftr_idx] = (vdr_FeatureTypeIsCdsOrMatPeptideOrGene($ftr_info_AHR, $ftr_idx)) ?
              opt_Get("--minpvlen", $opt_HHR) : 1;
        }
      }

      for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
        if(check_for_tabular_ftr_feature_prediction(\%{$ftr_results_HAHR->{$seq_name}[$ftr_idx]})) { 
          $seq_ntabftr++; 
        }
        if(check_for_valid_ftbl_feature_prediction(\%{$ftr_results_HAHR->{$seq_name}[$ftr_idx]}, $ftr_min_len_HA{$mdl_name}[$ftr_idx])) { 
          # initialize
          my $feature_type            = $ftr_info_AHR->[$ftr_idx]{"type"}; # type of feature, e.g. 'CDS' or 'mat_peptide' or 'gene'
          my $orig_feature_type       = $feature_type;                     # original feature type ($feature_type could be changed to misc_feature)
          my $is_5trunc_term_or_n     = 0;  # '1' if first segment of this feature is truncated at the 5' end due to sequence terminus or ambiguities
          my $is_3trunc_term_or_n     = 0;  # '1' if final segment of this feature is truncated at the 3' end due to sequence terminus or ambiguities
          my $is_misc_feature         = 0;  # '1' if this feature turns into a misc_feature due to alert(s)
          my $ftr_ftbl_coords_str     = "";    # string of coordinates for this feature
          my $ftr_ftbl_coords_len     = undef; # length of feature, in feature table coords 
                                               # (possibly shorter than actual feature length due to truncations due to ambiguities)
          my $ftr_out_str             = ""; # output string for this feature
          my $is_cds_or_mp            = vdr_FeatureTypeIsCdsOrMatPeptide($ftr_info_AHR, $ftr_idx);
          my $is_cds                  = vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx);
          my $parent_ftr_idx          = vdr_FeatureParentIndex($ftr_info_AHR, $ftr_idx); # will be -1 if has no parents
          my $parent_is_cds           = ($parent_ftr_idx == -1) ? 0 : vdr_FeatureTypeIsCds($ftr_info_AHR, $parent_ftr_idx);
          my $is_cds_or_parent_is_cds = ($is_cds || $parent_is_cds) ? 1 : 0;
          my $min_coord               = undef; # minimum coord in this feature
          my $cds_codon_start         = undef; # codon start value, only set for CDS

          my $defined_n_start   = (defined $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"n_start"}) ? 1: 0;
          my $defined_p_qstart   = (defined $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"p_qstart"}) ? 1: 0;
          my $ftr_is_trimmable  = (($do_notrim) || (defined $noftrtrim_H{$feature_type})) ? 0 : 1; # should we possible trim this feature due to ambiguities

          # sanity check
          if($is_cds && $parent_is_cds) { 
            ofile_FAIL("ERROR in $sub_name, feature $ftr_idx is a CDS and its parent is a CDS", 1, $FH_HR);
          }          
          
          # determine coordinates for the feature
          if(! $defined_n_start) { 
            # $defined_p_qstart must be TRUE
            ($ftr_ftbl_coords_str, $ftr_ftbl_coords_len, $min_coord, 
             $is_5trunc_term_or_n, $is_3trunc_term_or_n) =
                 helper_ftable_coords_prot_only_prediction($seq_name, $ftr_idx, $ftr_results_HAHR, $FH_HR);
            # note: $is_5trunc_{term,n} will always be 0, we don't try to do truncations for protein only predictions
            # note: $is_3trunc_{term,n} will alwyas be 0, we don't try to do truncations for protein only predictions
          }
          else { # $defined_n_start is '1'
            # fill $ftr_start_non_ab and $ftr_stop_non_ab if this feature may be trimmed due to N 
            my $ftr_start_non_ab = undef;
            my $ftr_stop_non_ab  = undef;
            if($ftr_is_trimmable) { 
              my $trim_idx = ($parent_is_cds) ? $parent_ftr_idx : $ftr_idx; # use parent if parent is a cds (e.g. mat_peptides)
              $ftr_start_non_ab = $ftr_results_HAHR->{$seq_name}[$trim_idx]{"n_start_non_ab"};
              $ftr_stop_non_ab  = $ftr_results_HAHR->{$seq_name}[$trim_idx]{"n_stop_non_ab"};
              #printf("set ftr_start_non_ab for ftr: $ftr_idx based on trim_idx: $trim_idx to %s\n", (defined $ftr_start_non_ab) ? $ftr_start_non_ab : "undef");
              #printf("set ftr_stop_non_ab  for ftr: $ftr_idx based on trim_idx: $trim_idx to %s\n", (defined $ftr_stop_non_ab)  ? $ftr_stop_non_ab : "undef");
            }
            ($ftr_ftbl_coords_str, $ftr_ftbl_coords_len, $min_coord, 
             $is_5trunc_term_or_n, $is_3trunc_term_or_n) = 
                 helper_ftable_coords_from_nt_prediction($seq_name, $ftr_idx, $ftr_start_non_ab, $ftr_stop_non_ab, 
                                                         $ftr_info_AHR, \%{$sgm_results_HHAHR->{$mdl_name}}, $FH_HR);
          }
          if($ftr_ftbl_coords_str ne "") { # if $ftr_ftbl_coords_str is "", we won't output the feature because it was entirely ambiguities
            # fill an array and strings with all alerts for this sequence/feature combo
            my $ftr_alt_str = helper_output_feature_alert_strings($seq_name, $ftr_idx, 0, $alt_info_HHR, \@ftr_alt_code_A, $alt_ftr_instances_HHHR, $FH_HR);
            my ($have_fatal_alt, $have_misc_alt) = helper_ftable_process_feature_alerts($ftr_alt_str, $seq_name, $ftr_idx, $ftr_info_AHR, $alt_info_HHR, $alt_ftr_instances_HHHR, \@seq_alert_A, $FH_HR);
            # should we make this a misc_feature?
            # yes if:
            # - --nomisc not enabled OR we have >=1 'misc_not_failure' feature/alert but zero fatal alerts
            # AND 
            # - feature type is not one of our hard-coded list of feature types that never get misc_feature-ized
            if($have_fatal_alt || $have_misc_alt) { 
              if((! $do_nomisc) || ((! $have_fatal_alt) && ($have_misc_alt))) { 
                if(vdr_FeatureTypeCanBecomeMiscFeature($ftr_info_AHR, $ftr_idx)) { 
                  $is_misc_feature = 1;
                  $feature_type = "misc_feature";
                  if((! $have_fatal_alt) && ($have_misc_alt)) { 
                    $cur_nmiscftr++; 
                  }
                }
              }
            }
            # determine codon_start if CDS
            if($is_cds) { 
              if(! $defined_n_start) { 
                $cds_codon_start = 1; # protein only prediction, codon start must be 1
              }
              else { 
                # n_start is defined, we have a nt prediction, we should have n_codon_start_expected
                # sanity check
                if(! defined $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"n_codon_start_expected"}) { 
                  ofile_FAIL("ERROR in $sub_name, sequence $seq_name CDS feature (ftr_idx: $ftr_idx) has no codon_start info", 1, $FH_HR);
                }
                $cds_codon_start = $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"n_codon_start_expected"};
                # if we trimmed the CDS start due to ambiguities update frame for that

                if(($ftr_is_trimmable) &&
                   (defined $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"n_5ablen"}) && 
                   ($ftr_results_HAHR->{$seq_name}[$ftr_idx]{"n_5ablen"} > 0)) { 
                  $cds_codon_start = vdr_FrameAdjust($cds_codon_start, $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"n_5ablen"}, $FH_HR);
                }
              } # end of else entered if n_start defined (codon_start block)
            } # end of 'if($is_cds)' entered to determine codon_start
            
            # convert coordinate string to output string
            $ftr_out_str = helper_ftable_coords_to_out_str($ftr_ftbl_coords_str, $feature_type, $FH_HR);
            
            # add qualifiers: product, gene, exception and codon_start
            if(! $is_misc_feature) { 
              $ftr_out_str .= helper_ftable_add_qualifier_from_ftr_info($ftr_idx, "product", $qval_sep, $ftr_info_AHR, $FH_HR);
              if(! $is_cds_or_mp) { 
                $ftr_out_str .= helper_ftable_add_qualifier_from_ftr_info($ftr_idx, "gene", $qval_sep, $ftr_info_AHR, $FH_HR);
              }

              my $ftr_nsgm = $ftr_ftbl_coords_str =~ tr/\n//; # counts number of lines of ftr_ftbl_coords_str (this is number of segments)
              if($ftr_nsgm > 1) { # only annotate ribsomal_slippage if more than one segments exist
                $ftr_out_str .= helper_ftable_add_qualifier_from_ftr_info($ftr_idx, "ribosomal_slippage", $qval_sep, $ftr_info_AHR, $FH_HR);
              }
              # have to be a little careful with 'exception' because there's a special case: 
              # "exception":"ribosomal slippage" should only be added if we have > 1 segment
              my $exception_str = helper_ftable_add_qualifier_from_ftr_info($ftr_idx, "exception", $qval_sep, $ftr_info_AHR, $FH_HR);
              if(($exception_str =~ /\t\t\texception\tribosomal slippage\n/) && ($ftr_nsgm == 1)) { 
                # remove ribosomal slippage if it exists
                $exception_str =~ s/\t\t\texception\tribosomal slippage\n//;
              }
              $ftr_out_str .= $exception_str;
              
              # add ncRNA_class qualifiers, if any
              $ftr_out_str .= helper_ftable_add_qualifier_from_ftr_info($ftr_idx, "ncRNA_class", $qval_sep, $ftr_info_AHR, $FH_HR);
              
              # add note qualifiers, if any
              $ftr_out_str .= helper_ftable_add_qualifier_from_ftr_info($ftr_idx, "note", $qval_sep, $ftr_info_AHR, $FH_HR);

              # if CDS, append the codon start only if we are truncated
              if($is_cds && $is_5trunc_term_or_n) { 
                $ftr_out_str .= helper_ftable_add_qualifier_specified($ftr_idx, "codon_start", $cds_codon_start, $FH_HR);
              }

              if((! $do_noprotid) && ($is_cds_or_parent_is_cds)) { 
                # add protein_id if we are a cds or parent is a cds
                # determine index for th protein_id qualifier
                my $protein_id_ftr_idx = ($is_cds) ? $ftr_idx : $parent_ftr_idx; # if !$is_cds, parent must be cds
                my $protein_id_idx = undef;
                # determine index for this protein
                if(defined $ftr_idx2protein_id_idx_H{$protein_id_ftr_idx}) { 
                  # the CDS itself or at least one mat_peptide with this
                  # CDS as its parent was already output, so use the same
                  # index that feature used
                  $protein_id_idx = $ftr_idx2protein_id_idx_H{$protein_id_ftr_idx};
                }
                else { 
                  # no index for this CDS yet exists, create it
                  $nprotein_id++;
                  $protein_id_idx = $nprotein_id;
                  $ftr_idx2protein_id_idx_H{$protein_id_ftr_idx} = $protein_id_idx;
                }
                
                # determine the protein_id value
                # - this cannot exceed $max_protein_id_length (50) characters as per GenBank rules (see github issue #12)
                #   so we shorten it to 50 characters if necessary UNLESS --forceprotid OR --noseqnamemax are used in which 
                #   case we assume user doesn't care about GenBank maximum
                # - first we try <seqname>_<index_of_protein_id_for_this_seq>, if this is <= $max_protein_id_length then we use that,
                #   if not, then we add a new suffix "_seq<seqidx>_<index_of_protein_id_for_this_seq>" at prepend the 
                #   first $max_protein_id_length - length(suffix) characters of the sequence name to it
                my $protein_id_value = sprintf("%s" . "_" . "%d", (($do_forceprotid) ? $seq_name : get_accession_from_ncbi_seq_name($seq_name)), $protein_id_idx);
                if((! $do_forceprotid) && (! $do_noseqnamemax)) { # neither --forceprotid and --noseqnamemax used
                  # make sure length of protein_id_value doesn't exceed the maximum, if so, shorten it.
                  if((length($protein_id_value)) > $max_protein_id_length) { 
                    my $new_sfx = sprintf("...seq%d_%d", $seq_idx2print, $protein_id_idx);
                    my $len_new_sfx = length($new_sfx);
                    if($len_new_sfx > $max_protein_id_length) { 
                      ofile_FAIL("ERROR in $sub_name, suffix being used to prevent protein id from exceeding $max_protein_id_length characters is itself more than $max_protein_id_length characters:\n$new_sfx\n", 1, $FH_HR);
                    }
                    my $alt_seq_name = get_accession_from_ncbi_seq_name($seq_name);
                    if((length($alt_seq_name) + $len_new_sfx) <= $max_protein_id_length) { 
                      $protein_id_value = $alt_seq_name . $new_sfx;
                    }
                    else { 
                      $protein_id_value = substr($alt_seq_name, 0, ($max_protein_id_length - $len_new_sfx)) . $new_sfx;
                    }
                  }
                }
                $ftr_out_str .= helper_ftable_add_qualifier_specified($ftr_idx, "protein_id", $protein_id_value, $FH_HR);
              }
            }
            else { # we are a misc_feature, add the 'similar to X' note
              $ftr_out_str .= sprintf("\t\t\t%s\t%s\n", "note", "similar to " . $ftr_info_AHR->[$ftr_idx]{"outname"});
            }
            
            # push to the output hash
            %{$ftout_AH[$ftidx]} = ();
            $ftout_AH[$ftidx]{"5trunc_term_or_n"} = ($is_5trunc_term_or_n) ? 1 : 0;
            $ftout_AH[$ftidx]{"3trunc_term_or_n"} = ($is_3trunc_term_or_n) ? 1 : 0;
            $ftout_AH[$ftidx]{"mincoord"}         = $min_coord;
            $ftout_AH[$ftidx]{"type_priority"}    = (exists $type_priority_H{$orig_feature_type}) ? $type_priority_H{$orig_feature_type} : $npriority;
            $ftout_AH[$ftidx]{"coords"}           = $ftr_ftbl_coords_str;
            $ftout_AH[$ftidx]{"output"}           = $ftr_out_str;
            $ftout_AH[$ftidx]{"codon_start"}      = (defined $cds_codon_start) ? $cds_codon_start : -1;
            $ftout_AH[$ftidx]{"ftbl_len"}         = $ftr_ftbl_coords_len;
            $ftr_idx2ftout_idx_H{$ftr_idx} = $ftidx;
            $ftidx++;
          } # end of 'if($ftr_ftbl_coords_str ne "")'
        } # end of 'if(check_for_valid_ftbl_feature_prediction('
      } # end of 'for(my $ftr_idx...'

      # Final step before outputting: 
      # Possibly remove some output from @ftout_AH before outputting
      # We remove output for the following features:
      # 1. CDS that are too short to encode a single AA or
      #    are only 1 AA which is the stop codon.
      # 2. mat_peptides that are too short to encode a single AA
      # 3. any feature that has a parent that does not have its own
      #    feature output
      # 
      # This is mainly necessary because feature table feature lengths
      # can differ from actual feature lengths due to ambiguities at
      # the beginning and end of features, which can make some of them
      # too short to encode even one AA.
      #
      # We can't do this pruning earlier because mat_peptides need 
      # to have info (specifically length and codon_start info) 
      # from their parent CDS in order to determine if they should
      # be removed. And we do not enforce that children have to 
      # come after their parents in ftr_info_AHR so we are not 
      # guaranteed we have this info until we've generated it for
      # all features.
      
      # initialize
      my @remove_me_A = (); # [0..$i..(scalar(@ftout_AH)-1)]: 1 to remove hash of output with index $i
      my $pre_remove_noutftr = scalar(@ftout_AH);
      for($ftidx = 0; $ftidx < $pre_remove_noutftr; $ftidx++) { 
        $remove_me_A[$ftidx] = 0;
      }
      
      # remove output for CDS and MPs that are too short
      for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
        if(defined $ftr_idx2ftout_idx_H{$ftr_idx}) {
          $ftidx = $ftr_idx2ftout_idx_H{$ftr_idx};
          my $is_cds  = vdr_FeatureTypeIsCds($ftr_info_AHR, $ftr_idx) ? 1 : 0;
          my $is_mp   = vdr_FeatureTypeIsMatPeptide($ftr_info_AHR, $ftr_idx) ? 1 : 0;
          my $parent_ftr_idx = vdr_FeatureParentIndex($ftr_info_AHR, $ftr_idx); # will be -1 if no parent
          if($is_cds) {
            my $ftr_ftidx           = $ftr_idx2ftout_idx_H{$ftr_idx};
            my $ftbl_len            = $ftout_AH[$ftr_ftidx]{"ftbl_len"};
            my $is_3trunc_term_or_n = $ftout_AH[$ftr_ftidx]{"3trunc_term_or_n"}; # 3' truncated due to sequence terminus and/or ambigs
            my $codon_start         = $ftout_AH[$ftr_ftidx]{"codon_start"};
            # is it too short? 4 cases:
            if(($ftbl_len < 3)                            || # less than 1 AA, regardless of frame
               (($ftbl_len == 3) && ($codon_start != 1))  || # less than 1 AA, frame 2 or 3
               (($ftbl_len == 4) && ($codon_start == 3))  || # less than 1 AA, frame 3
               (($ftbl_len <= 5) && (! $is_3trunc_term_or_n))) { # only a stop codon
              $remove_me_A[$ftidx] = 1;
            }
          }
          elsif(($is_mp) && # mat_peptide
                ($parent_ftr_idx != -1) && # with parent
                (defined $ftr_idx2ftout_idx_H{$parent_ftr_idx})) { # parent has output 
            # (if parent does not have output we'll remove this MP in loop below that checks for parentless output)
            my $ftr_ftidx         = $ftr_idx2ftout_idx_H{$ftr_idx};
            my $parent_ftr_ftidx  = $ftr_idx2ftout_idx_H{$parent_ftr_idx};
            my $ftbl_len          = $ftout_AH[$ftr_ftidx]{"ftbl_len"};
            my $codon_start       = $ftout_AH[$parent_ftr_ftidx]{"codon_start"};
            # is it too short? 3 cases:
            if(($ftbl_len < 3)                           || # less than 1 AA, regardless of frame
               (($ftbl_len == 3) && ($codon_start != 1)) || # less than 1 AA, frame 2 or 3
               (($ftbl_len == 4) && ($codon_start == 3))) {  # less than 1 AA, frame 3 (don't need to check for only stop codon case for MPs
              $remove_me_A[$ftidx] = 1;
            }
          }
        }
      }

      # go back through and remove output for any feature which has a parent that does not have output itself
      for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
        if(defined $ftr_idx2ftout_idx_H{$ftr_idx}) { 
          $ftidx = $ftr_idx2ftout_idx_H{$ftr_idx};
          if(! $remove_me_A[$ftidx]) { # we're not removing this output yet
            my $parent_ftr_idx = vdr_FeatureParentIndex($ftr_info_AHR, $ftr_idx); # will be -1 if no parent
            if($parent_ftr_idx != -1) { # we have a parent
              if((! defined $ftr_idx2ftout_idx_H{$parent_ftr_idx}) ||     # parent has no output
                 ($remove_me_A[$ftr_idx2ftout_idx_H{$parent_ftr_idx}])) { # parent has output but we are removing it
                $remove_me_A[$ftidx] = 1;
              }
            }
          }
        }
      }
      
      # actually remove the output
      for($ftidx = ($pre_remove_noutftr-1); $ftidx >= 0; $ftidx--) { 
        # descending so we can use splice without messing up indices
        if($remove_me_A[$ftidx]) { 
          splice(@ftout_AH, $ftidx, 1);
        }
      }
    } # end of 'if(defined $mdl_name)'

    #######################################
    # OUTPUT section 
    #######################################
    # done with this sequence, determine what type of output we will have 
    my $cur_noutftr = scalar(@ftout_AH);

    # possibly add nmiscftr, noftrann, and noftrant alerts
    if(($cur_nmiscftr > 0) && ($cur_nmiscftr >= $nmiscftr_thr)) { 
      # more than the maximum allowed number of misc_features have been created due to the misc_not_failure attributes
      alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "nmiscftr", $seq_name, sprintf("seq:VADRNULL;mdl:VADRNULL;%d>=%d", $cur_nmiscftr, $nmiscftr_thr), $FH_HR);
      # only add to @seq_alert_A if nmiscftr is fatal
      if($alt_info_HHR->{"nmiscftr"}{"causes_failure"}) { 
        push(@seq_alert_A, sprintf("%s: (*sequence*) S:-; M:-; %s %s", $alt_info_HHR->{"nmiscftr"}{"sdesc"}, $alt_info_HHR->{"nmiscftr"}{"ldesc"}, sprintf("[%d>=%d]", $cur_nmiscftr, $nmiscftr_thr)));
      }
    }
    if($seq_ntabftr == 0) { 
      # no features annotated, even in eventual .ftr tabular file
      # first check to see if any per-sequence alerts that prevent annotation have already been reported,
      # if so, we don't report noftrann because there is no annotation at all for this sequence
      if(! $prevents_annot_flag) { 
        alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "noftrann", $seq_name, "VADRNULL", $FH_HR);
        # set @seq_alert_A to this lone alert
        @seq_alert_A = ();
        push(@seq_alert_A, sprintf("%s: (*sequence*) S:-; M:-; %s %s", $alt_info_HHR->{"noftrann"}{"sdesc"}, $alt_info_HHR->{"noftrann"}{"ldesc"}, "[-]"));
      }
    }
    elsif($cur_noutftr == 0) { 
      # >= 1 features annotated in eventual .ftr tbl file, but zero in this .tbl file because 
      # they were all too short
      alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "noftrant", $seq_name, "VADRNULL", $FH_HR);
      # set @seq_alert_A to this lone alert
      @seq_alert_A = ();
      push(@seq_alert_A, sprintf("%s: (*sequence*) S:-; M:-; %s %s", $alt_info_HHR->{"noftrant"}{"sdesc"}, $alt_info_HHR->{"noftrant"}{"ldesc"}, "[-]"));
    }

    # determine if the sequence will pass
    # the sequences only pass if:
    # 1) at least one feature is annotated ($cur_noutftr > 0)
    # 2) zero fatal sequence alerts and feature alerts for features annotated in the feature table
    # 3) zero fatal feature alerts for features NOT annotated in the feature table (because they are too short)
    # we need criteria 3 so that all seqs that FAIL in tabular output files also FAIL in feature table files 
    # (e.g. are listed in fail.tbl instead of .pass.tbl)
    my $do_pass = (($cur_noutftr > 0) && ((scalar(@seq_alert_A)) == 0)) ? 1 : 0; # checks only that criteria 1 and 2 are met

    # next line checks for criteria 3 (only if first 2 criteria have been met)
    if($do_pass && (! check_if_sequence_passes($seq_name, ((defined $mdl_name) ? \@{$ftr_info_HAHR->{$mdl_name}} : undef), $alt_info_HHR, $alt_seq_instances_HHR, $alt_ftr_instances_HHHR, $FH_HR))) { 
      # sequence has zero fatal sequence alerts *and* zero fatal feature alerts for all features output to feature table
      # BUT at least one fatal alert for a feature NOT output to feature table (e.g. a feature that is too short to meet
      # minimum length requirements for the feature table). We throw a special alert here (ftskipfl) for this
      # so the sequence fails both in tabular and feature table output.
      alert_sequence_instance_add($alt_seq_instances_HHR, $alt_info_HHR, "ftskipfl", $seq_name, "seq:VADRNULL;mdl:VADRNULL;see .ftr and .alt output files for details", $FH_HR);
      # set @seq_alert_A to this lone alert
      push(@seq_alert_A, sprintf("%s: (*sequence*) S:-; M:-; %s %s", $alt_info_HHR->{"ftskipfl"}{"sdesc"}, $alt_info_HHR->{"ftskipfl"}{"ldesc"}, "[see .ftr and .alt output files for details]"));
      $do_pass = 0; # this seq fails
    }

    # sort output
    if($cur_noutftr > 0) { 
      @ftout_AH = sort { $a->{"mincoord"}             <=> $b->{"mincoord"} or 
                             $b->{"5trunc_term_or_n"} <=> $a->{"5trunc_term_or_n"} or
                             $a->{"3trunc_term_or_n"} <=> $b->{"3trunc_term_or_n"} or
                             $a->{"type_priority"}    <=> $b->{"type_priority"} 
      } @ftout_AH;
    }              
      
    if($do_pass) { 
      # print to the passing feature table file
      $ret_npass++;
      print $pass_list_FH $seq_name . "\n";
      push(@pass_fa_A, $seq_name);
      print $pass_ftbl_FH ">Feature $seq_name\n";
      for($i = 0; $i < scalar(@ftout_AH); $i++) { 
        print $pass_ftbl_FH $ftout_AH[$i]{"output"};
      }
    }
    else { # $do_pass == 0
      print $fail_list_FH $seq_name . "\n";
      push(@fail_fa_A, $seq_name);
      print $fail_ftbl_FH ">Feature $seq_name\n";
      for($i = 0; $i < scalar(@ftout_AH); $i++) { 
        print $fail_ftbl_FH $ftout_AH[$i]{"output"};
      }
      if((scalar(@seq_alert_A)) > 0) { 
        my $mdl_name2print = (defined $mdl_name) ? $mdl_name : "-";
        print $fail_ftbl_FH "\nAdditional note(s) to submitter:\n"; 
        for(my $e = 0; $e < scalar(@seq_alert_A); $e++) { 
          my $error_line = $seq_alert_A[$e];
          my $error_line2print = $error_line; # modified below
          if($error_line =~ /(^[^\:]+)\:\s+\((.*)\)\s+S\:([^\;]+)\;\s+M\:([^\;]+)\;\s+(.+)$/) {
            my ($tmp_sdesc, $tmp_ftr_type_and_name, $tmp_scoords, $tmp_mcoords, $tmp_edesc) = ($1, $2, $3, $4, $5);
            my ($tmp_ftr_type, $tmp_ftr_name) = ("-", "*sequence*"); # 
            if ($tmp_ftr_type_and_name ne "*sequence*") { 
              if($tmp_ftr_type_and_name =~ /([^\:]+)\:(.*)$/) { 
                ($tmp_ftr_type, $tmp_ftr_name) = ($1, $2);
              }
              else { 
                ofile_FAIL("ERROR in $sub_name, unable to split ftr_type_and_name:$tmp_ftr_type_and_name for output from line $error_line", 1, $ofile_info_HHR->{"FH"});
              }
            }
            print $alerts_FH    $seq_name . "\t" . $mdl_name2print . "\t" . $tmp_ftr_type . "\t" . $tmp_ftr_name . "\t" . $tmp_sdesc . "\t" . $tmp_scoords . "\t" . $tmp_mcoords . "\t" . $tmp_edesc . "\n";
            print $fail_ftbl_FH "ERROR: " . $tmp_sdesc . ": (" . $tmp_ftr_type_and_name . ") " . $tmp_edesc . "; seq-coords:" . $tmp_scoords . "; mdl-coords:" . $tmp_mcoords . "; mdl:" . $mdl_name2print . ";\n";
          }
          else {
            ofile_FAIL("ERROR in $sub_name, unable to split alert_line for output: $error_line", 1, $ofile_info_HHR->{"FH"});
          }
        }
      } # end of 'if($cur_nalert > 0)'
    }
  } # end of loop over sequences

  if(! $do_nofasta) { 
    my $pass_fa_file = $out_root . ".pass.fa";
    my $fail_fa_file = $out_root . ".fail.fa";
    $$in_sqfile_R->fetch_seqs_given_names(\@pass_fa_A, 60, $pass_fa_file);
    $$in_sqfile_R->fetch_seqs_given_names(\@fail_fa_A, 60, $fail_fa_file);
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "pass_fa", $pass_fa_file, 1, 1, "fasta file with passing sequences");
    ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, "fail_fa", $fail_fa_file, 1, 1, "fasta file with failing sequences");
  }

  return $ret_npass;
}

#################################################################
# Subroutine:  helper_ftable_coords_from_nt_prediction
# Incept:      EPN, Tue Oct 30 12:59:13 2018
#
# Purpose:    Given a sequence name and feature index, construct
#             a feature table coordinate string, possibly of 
#             multiple lines, one per segment.
#
# Arguments: 
#  $seq_name:          sequence name
#  $ftr_idx:           feature index
#  $start_non_ab:       first position of feature that is not an N
#  $stop_non_ab:        final position of feature that is not an N
#  $ret_min_coord:     REF to minimum coordinate, to fill
#  $ftr_info_AHR:      REF to array of hashes with information on the features, PRE-FILLED
#  $sgm_results_HAHR:  REF to segment results HAH, PRE-FILLED
#  $FH_HR:             REF to hash of file handles
#
# Returns:    Five values:
#             $ftr_ftbl_coords_str: string that gives the coordinates for this feature in feature table format
#                                   "" if entire feature is ambigs, in this case we won't output it to feature table
#             $ftr_ftbl_coords_len: length in nt of feature in output coords for feature table
#             $min_coord:           minimum coordinate for feature
#             $is_5trunc:           '1' if first segment is truncated on 5' end due to sequence terminus or ambigs
#             $is_3trunc:           '1' if final segment is truncated on 3' end due to sequence terminus or ambigs
#
# Dies: if either @{$start_AR} or @{$stop_AR} are empty
#       if $start_non_ab is -1 but stop_non_ab is not
#
################################################################# 
sub helper_ftable_coords_from_nt_prediction { 
  my $sub_name = "helper_ftable_coords_from_nt_prediction";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $ftr_idx, $start_non_ab, $stop_non_ab, $ftr_info_AHR, $sgm_results_HAHR, $FH_HR) = @_;

  # arrays with per-sgm info
  my @start_A     = ();
  my @stop_A      = ();
  my @strand_A    = ();
  my @is_5trunc_A = ();
  my @is_3trunc_A = ();
  
  for(my $sgm_idx = $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}; $sgm_idx <= $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"}; $sgm_idx++) { 
    if(defined $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"}) { 
      push(@start_A,     $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"});
      push(@stop_A,      $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstop"});
      push(@strand_A,    $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"strand"});
      push(@is_5trunc_A, $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"5trunc"});
      push(@is_3trunc_A, $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"3trunc"});
    }
  }
  return helper_ftable_start_stop_strand_arrays_to_coords(\@start_A, \@stop_A, \@strand_A, \@is_5trunc_A, \@is_3trunc_A, 
                                                          $start_non_ab, $stop_non_ab, $FH_HR);
}

#################################################################
# Subroutine:  helper_ftable_coords_prot_only_prediction
# Incept:      EPN, Tue Oct 30 12:26:05 2018
#
# Purpose:    Given a sequence name and feature index, construct
#             a feature table coordinate string, possibly of 
#             multiple lines, one per segment, for the special
#             case that this feature has a 'indfantp' alert: blastx 
#             prediction but no CM prediction.
#
# Arguments: 
#  $seq_name:         sequence name
#  $ftr_idx:          feature index
#  $ftr_results_HAHR: REF to feature results AAH, PRE-FILLED
#  $FH_HR:            REF to hash of file handles
#
# Returns:    Five values:
#             $ftr_ftbl_coords_str: string that gives the coordinates for this feature in feature table format
#             $ftr_ftbl_coords_len: length in nt of feature in output coords for feature table
#             $min_coord:           minimum coordinate for feature
#             $is_5trunc:           always 0, we don't try to truncate for protein only predictions
#             $is_3trunc:           always 0, we don't try to truncate for protein only predictions
# 
# Dies:       if p_qstart or p_qstop does not exist in the ftr_results_HAHR->{$seq_name}[$ftr_idx] hash
#             If $N_{5,3}trunc is > 0 and corresponding strand is not "+" or "-"
################################################################# 
sub helper_ftable_coords_prot_only_prediction { 
  my $sub_name = "helper_ftable_coords_prot_only_prediction";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $ftr_idx, $ftr_results_HAHR, $FH_HR) = @_;

  # NOTE: for 'indfantp' alerts, the p_qstart and p_qstop are always set at the feature level
  if((! exists $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"p_qstart"}) ||
     (! exists $ftr_results_HAHR->{$seq_name}[$ftr_idx]{"p_qstop"})) { 
    ofile_FAIL("ERROR in $sub_name, ftr_results_HAHR->{$seq_name}[$ftr_idx]{p_qstart|p_qstop} does not exists", 1, $FH_HR);
  }

  my @start_A     = ($ftr_results_HAHR->{$seq_name}[$ftr_idx]{"p_qstart"});
  my @stop_A      = ($ftr_results_HAHR->{$seq_name}[$ftr_idx]{"p_qstop"});
  my @strand_A    = ($ftr_results_HAHR->{$seq_name}[$ftr_idx]{"p_strand"});
  my @is_5trunc_A = (0); # can't detect truncation for protein predictions, currently
  my @is_3trunc_A = (0); # can't detect truncation for protein predictions, currently

  return helper_ftable_start_stop_strand_arrays_to_coords(\@start_A, \@stop_A, \@strand_A, \@is_5trunc_A, \@is_3trunc_A, 
                                                          undef, undef, $FH_HR);
}

#################################################################
# Subroutine:  helper_ftable_start_stop_strand_arrays_to_coords()
# Incept:      EPN, Tue Oct 30 12:39:59 2018
#
# Purpose:    Given refs to two arrays of start and stop coordinates,
#             construct coordinate strings in feature table format.
#
# Arguments: 
#  $start_AR:      REF to array of start coordinates, one per sgm
#  $stop_AR:       REF to array of stop coordinates, one per sgm
#  $strand_AR:     REF to array of strands, one per sgm
#  $is_5trunc_AR:  REF to array of is_5trunc values, one per sgm
#  $is_3trunc_AR:  REF to array of is_3trunc values, one per sgm
#  $start_non_ab:   first position of feature that is not an N (may be > $stop_non_ab)
#  $stop_non_ab:    final position of feature that is not an N (may be < $start_non_ab)
#  $FH_HR:         REF to hash of file handles
#
# Returns:    Five values:
#             $ftr_ftbl_coords_str: string that gives the coordinates for this feature in feature table format
#                                   "" if entire feature is ambigs, in this case we won't output it to feature table
#             $ftr_ftbl_coords_len: length in nt of feature in output coords for feature table
#             $min_coord:           minimum coordinate for feature
#             $is_5trunc:           '1' if first segment is truncated on 5' end due to sequence terminus or ambigs
#             $is_3trunc:           '1' if final segment is truncated on 3' end due to sequence terminus or ambigs
#
# Dies: if either @{$start_AR} or @{$stop_AR} are empty
#       if $start_non_ab is -1 but stop_non_ab is not
#
################################################################# 
sub helper_ftable_start_stop_strand_arrays_to_coords { 
  my $sub_name = "helper_ftable_start_stop_strand_arrays_to_coords";
  my $nargs_exp = 8;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($start_AR, $stop_AR, $strand_AR, $is_5trunc_AR, $is_3trunc_AR, $start_non_ab, $stop_non_ab, $FH_HR) = @_;

  # return values
  my $ret_ftr_ftbl_coords_str = "";
  my $ret_ftr_ftbl_coords_len = 0;
  my $ret_min_coord = undef; # minimum coordinate output to table
  my $ret_is_5trunc_term_or_n_first_sgm = undef; # set to '1' if first segment is 5' truncated due to sequence terminus or ambigs, '0' if not
  my $ret_is_3trunc_term_or_n_final_sgm = undef; # set to '1' if final segment is 3' truncated due to sequence terminus or ambigs, '0' if not

  my $ncoord = scalar(@{$start_AR});
  if($ncoord == 0) { 
    ofile_FAIL("ERROR in $sub_name, start_A array is empty", 1, $FH_HR);
  }

  my ($min_non_ab, $max_non_ab) = (undef, undef);
  if((defined $start_non_ab) && (defined $stop_non_ab)) { 
    ($min_non_ab, $max_non_ab) = ($start_non_ab, $stop_non_ab);
    if($min_non_ab > $max_non_ab) { utl_Swap(\$min_non_ab, \$max_non_ab); }
  }

  for(my $c = 0; $c < $ncoord; $c++) { 
    my $is_first = ($c == 0)           ? 1 : 0;
    my $is_final = ($c == ($ncoord-1)) ? 1 : 0;
    my $start     = $start_AR->[$c];
    my $stop      = $stop_AR->[$c];
    my $strand    = $strand_AR->[$c];
    my $is_5trunc_term = $is_5trunc_AR->[$c]; # segment is 5' truncated due to sequence terminus
    my $is_3trunc_term = $is_3trunc_AR->[$c]; # segment is 3' truncated due to sequence terminus

    # potentially modify start/stop based on $start_non_ab and $stop_non_ab
    my $is_5trunc_n = 0; # set to 1 below if start position is truncated due to ambigs
    my $is_3trunc_n = 0; # set to 1 below if stop  position is truncated due to ambigs
    my $add_this_sgm = 1; # set to 0 below if full sgm is ambigs, in which case we don't add it

    if((defined $start_non_ab) && (defined $stop_non_ab)) { 
      if($start_non_ab == -1) { # this means entire segment is ambigs
        if($stop_non_ab != -1) { # sanity check 
          ofile_FAIL("ERROR, in $sub_name, start_non_ab is -1 but stop_non_ab is not ($stop_non_ab)", 1, $FH_HR);
        }
        $add_this_sgm = 0;
      }
      else { 
        # get min/max between start/stop and start_non_ab/stop_non_ab 
        # to make the following complicated checks of cases a little less complicated
        my ($min, $max) = ($start, $stop);
        if($min > $max) { utl_Swap(\$min, \$max); }

        if(($min > $max_non_ab) ||  # $min_non_ab <= $max_non_ab < $min       <= $max
           ($max < $min_non_ab)) {  # $min       <= $max       < $min_non_ab <= $max_non_ab
          # full sgm is starts/ends before $min_non_ab or after $max_non_ab, don't output it
          $add_this_sgm = 0; # don't add it
        }
        else { 
          if($min < $min_non_ab) { # minimum starts before min_non_ab
            if($start == $min) { $start = $min_non_ab; $is_5trunc_n = 1; }
            if($stop  == $min) { $stop  = $min_non_ab; $is_3trunc_n = 1; }
          }
          if($max > $max_non_ab) { # maximum ends after max_non_ab
            if($start == $max) { $start = $max_non_ab; $is_5trunc_n = 1; }
            if($stop  == $max) { $stop  = $max_non_ab; $is_3trunc_n = 1; }
          }
        }
      }
    }
    if($add_this_sgm) { 
      if((! defined $ret_min_coord) || ($start < $ret_min_coord)) { $ret_min_coord = $start; }
      if($stop < $ret_min_coord) { $ret_min_coord = $stop; }

      if (! defined $ret_is_5trunc_term_or_n_first_sgm) { 
        $ret_is_5trunc_term_or_n_first_sgm = ($is_5trunc_term || $is_5trunc_n) ? 1 : 0;
      }
      # always update final_sgm value
      $ret_is_3trunc_term_or_n_final_sgm = ($is_3trunc_term || $is_3trunc_n) ? 1 : 0;

      $ret_ftr_ftbl_coords_str .= sprintf("%s%d\t%s%d\n", 
                                 ($is_5trunc_term || $is_5trunc_n) ? "<" : "", $start, 
                                 ($is_3trunc_term || $is_3trunc_n) ? ">" : "", $stop);
      $ret_ftr_ftbl_coords_len += abs($stop - $start) + 1;
    }
  }
  if(! defined $ret_min_coord)                     { $ret_min_coord = -1; } # irrelevant, caller's responsibility to handle this
  if(! defined $ret_is_5trunc_term_or_n_first_sgm) { $ret_is_5trunc_term_or_n_first_sgm =  0; } # irrelevant, caller's responsibility to handle this
  if(! defined $ret_is_3trunc_term_or_n_final_sgm) { $ret_is_3trunc_term_or_n_final_sgm =  0; } # irrelevant, caller's responsibility to handle this
  if(! defined $ret_min_coord) { $ret_min_coord = -1; } # irrelevant, caller's responsibility to handle this

  return ($ret_ftr_ftbl_coords_str, $ret_ftr_ftbl_coords_len, $ret_min_coord, 
          $ret_is_5trunc_term_or_n_first_sgm, $ret_is_3trunc_term_or_n_final_sgm);
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
    ofile_FAIL("ERROR in $sub_name, coords_str is empty - shouldn't happen", 1, $FH_HR);
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
#             part of a feature table output, where the 
#             qualifier value is obtained from
#             @{$ftr_info_AHR}.
#
# Arguments: 
#  $ftr_idx:      feature index
#  $key:          key in ftr_info_AHR
#  $qval_sep:     characters that separate multiple qualifiers in $ftr_info_AHR->[$ftr_idx]{$key}
#  $ftr_info_AHR: REF to array of hashes with information on the features, PRE-FILLED
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

  my ($ftr_idx, $key, $qval_sep, $ftr_info_AHR, $FH_HR) = @_;
  
  my $ret_str = "";
  if((defined $ftr_info_AHR) && 
     (defined $ftr_info_AHR->[$ftr_idx]) && 
     (defined $ftr_info_AHR->[$ftr_idx]{$key})) { 
    my @qval_A = split($qval_sep, $ftr_info_AHR->[$ftr_idx]{$key});
    foreach my $qval (@qval_A) { 
      my @qval_A = split(":GBSEP:", $qval);
      foreach my $qval_line (@qval_A) { 
        if($qval_line eq "GBNULL") { 
          $ret_str .= sprintf("\t\t\t%s\n", $key);
        }
        else { 
          $ret_str .= sprintf("\t\t\t%s\t%s\n", $key, $qval_line);
        }
      }
    }
  }

  return $ret_str;
}

#################################################################
# Subroutine:  helper_ftable_add_qualifier_from_ftr_results()
# Incept:      EPN, Tue Oct 30 13:52:19 2018
#
# Purpose:    Add a qualifier line to a string that will be 
#             part of a feature table output, where the
#             qualifier value is obtained from 
#             %{$ftr_results_HAHR}.
#
# Arguments: 
#  $seq_name:          sequence name
#  $ftr_idx:           feature index
#  $results_key:       key in ftr_results_HAHR
#  $qualifier:         name for the qualifier
#  $ftr_results_HAHR:  REF to feature results HAH, PRE-FILLED
#  $FH_HR:             REF to hash of file handles
#
# Returns:    "" if $ftr_results_HAHR->{$seq_name}[$ftr_idx]{$results_key} does not exist
#             else a string for the feature table
#
# Dies: never
#
################################################################# 
sub helper_ftable_add_qualifier_from_ftr_results {
  my $sub_name = "helper_ftable_add_qualifier_from_ftr_results";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $ftr_idx, $results_key, $qualifier, $ftr_results_HAHR, $FH_HR) = @_;

  my $ret_str = "";

  if((defined $ftr_results_HAHR) && 
     (defined $ftr_results_HAHR->{$seq_name}) && 
     (defined $ftr_results_HAHR->{$seq_name}[$ftr_idx]) && 
     (defined $ftr_results_HAHR->{$seq_name}[$ftr_idx]{$results_key})) { 
    if($ftr_results_HAHR->{$seq_name}[$ftr_idx]{$results_key} eq "") { 
      $ret_str = sprintf("\t\t\t%s\n", $qualifier);
    }
    else { 
      $ret_str = sprintf("\t\t\t%s\t%s\n", $qualifier, $ftr_results_HAHR->{$seq_name}[$ftr_idx]{$results_key});
    }
  }
  return $ret_str;
}

#################################################################
# Subroutine:  helper_ftable_add_qualifier_specified()
# Incept:      EPN, Thu Mar 26 14:19:27 2020
#
# Purpose:    Add a qualifier line to a string that will be 
#             part of a feature table output, where that
#             qualifier and its value are specified by the
#             caller.
#
# Arguments: 
#  $ftr_idx:      feature index
#  $qualifier:    name for the qualifier
#  $value:        value for the qualifier (undef for no value)
#  $FH_HR:        REF to hash of file handles
#
# Returns:    String to append to the feature table.
#
# Dies: never
#
################################################################# 
sub helper_ftable_add_qualifier_specified {
  my $sub_name = "helper_ftable_add_qualifier_specified";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_idx, $qualifier, $value, $FH_HR) = @_;

  my $ret_str = "";

  if((defined $value) && ($value ne "")) { 
    $ret_str = sprintf("\t\t\t%s\t%s\n", $qualifier, $value);
  }
  else { 
    $ret_str = sprintf("\t\t\t%s\n", $qualifier);
  }

  return $ret_str;
}

#################################################################
# Subroutine:  helper_ftable_class_model_for_sequence()
# Incept:      EPN, Tue Mar 26 14:47:53 2019
#

# Purpose:    Given a sequence name and a reference of the classification
#             results 3D hash, determine what model the sequence was
#             classified to.
#
# Arguments: 
#  $stg_results_HHHR:      ref to classification results 3D hash
#  $seq_name:              name of sequence
#
# Returns:    model name $seq_name is classified to, or undef 
#             if none
#
################################################################# 
sub helper_ftable_class_model_for_sequence {
  my $sub_name = "helper_ftable_class_model_for_sequence";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($stg_results_HHHR, $seq_name) = @_;

  return ((defined $stg_results_HHHR->{$seq_name}) && 
          (defined $stg_results_HHHR->{$seq_name}{"std.cls.1"}) && 
          (defined $stg_results_HHHR->{$seq_name}{"std.cls.1"}{"model"})) ? 
          $stg_results_HHH{$seq_name}{"std.cls.1"}{"model"} : undef;
}

#################################################################
# Subroutine: helper_ftable_process_sequence_alerts()
# Incept:     EPN, Thu Jan 24 12:09:24 2019
#
# Purpose:    Given a string of per-sequence alerts that correspond
#             to a specific sequence, use the %{$alt_info_HHR} and
#             process that string to determine what (if any) 
#             alerts should be added to the feature table
#             for this sequence. Note that we do not add any 'notes'
#             as we possibly could in helper_ftable_process_feature_alerts()
#             because we are dealing with the full sequence and not
#             a feature for a sequence.
#
# Arguments:
#   $alt_code_str:           string of alerts, comma separated, can be ""
#   $seq_name:               name of sequence
#   $alt_info_HHR:           REF to hash of hashes with information on the alerts, PRE-FILLED
#   $alt_seq_instances_HHR:  REF to 2D hashes with per-sequence alerts, PRE-FILLED
#   $ret_alert_AR:           REF to array of alerts, possibly added to here (not created)
#   $FH_HR:                  REF to hash of file handles, including "log" and "cmd"
# 
# Returns: '1' if any alert added to $ret_alert_AR has the 'prevents_annot'
#          flag as '1', else '0'. Returns '0' if no alerts are added to 
#          $ret_alert_AR.
#
# Dies: Never
#################################################################
sub helper_ftable_process_sequence_alerts { 
  my $sub_name = "helper_ftable_process_sequence_alerts";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_code_str, $seq_name, $alt_info_HHR, $alt_seq_instances_HHR, $ret_alert_AR, $FH_HR) = (@_);

  if($alt_code_str eq "") { 
    return 0; 
  }

  my $ret_prevents_annot_flag = 0;

  # create a hash of all alerts in the input $alt_str, and also verify they are all valid errors
  my %input_alt_code_H = (); # $input_err_code_H{$alt_code} = 1 if $alt_code is in $alt_code_str
  my $alt_code; 
  foreach $alt_code (split(",", $alt_code_str)) { 
    if(! defined $alt_info_HHR->{$alt_code}) { 
      ofile_FAIL("ERROR in $sub_name, input error of $alt_code in string $alt_code_str is invalid", 1, $FH_HR);
    }
    $input_alt_code_H{$alt_code} = 1; 
  }

  my $do_report = 0; # '1' if we should report this alert in the feature table, '0' if not
  foreach $alt_code (sort keys (%input_alt_code_H)) { 
    $do_report = $alt_info_HHR->{$alt_code}{"causes_failure"}; # only report alerts that cause failure in the feature table
    # check if this alert is invalidated by another we will also report
    if(($do_report) && ($alt_info_HHR->{$alt_code}{"ftbl_invalid_by"} ne "")) { 
      my @invalid_by_alt_code_A = split(",", $alt_info_HHR->{$alt_code}{"ftbl_invalid_by"});
      foreach my $alt_code2 (@invalid_by_alt_code_A) {
        if(($alt_info_HHR->{$alt_code2}{"causes_failure"}) && 
           (exists $input_alt_code_H{$alt_code2})) { 
          $do_report = 0; # $alt_code is invalidated by $alt_code2, $alt_code2 causes failure and is also present in $alt_code_str
        }
      }
    }
    if($do_report) { 
      # determine if prevents_annot flag is '1'
      if($alt_info_HHR->{$alt_code}{"prevents_annot"}) { 
        $ret_prevents_annot_flag = 1;
      }
      # we could have more than one instance of this sequence/alert pair
      my @instance_str_A = split(":VADRSEP:", $alt_seq_instances_HHR->{$seq_name}{$alt_code});
      foreach my $instance_str (@instance_str_A) { 
        my ($instance_scoords, $instance_mcoords, $instance_detail) = alert_instance_parse($instance_str);
        # alert_instance_parse removes 'seq:' and ';', 'mdl:', and ';'
        my $alert_str = sprintf("%s: (*sequence*) S:%s; M:%s; %s %s", 
                                $alt_info_HHR->{$alt_code}{"sdesc"}, 
                                ($instance_scoords eq "VADRNULL") ? "-" : $instance_scoords, 
                                ($instance_mcoords eq "VADRNULL") ? "-" : $instance_mcoords, 
                                $alt_info_HHR->{$alt_code}{"ldesc"}, 
                                ($instance_detail ne "VADRNULL") ? "[" . $instance_detail . "]" : "[-]");
        # only add the alert, if an identical alert does not already exist in @{$ret_alert_AR}
        my $idx = utl_AFindNonNumericValue($ret_alert_AR, $alert_str, $FH_HR);
        if($idx == -1) { 
          push(@{$ret_alert_AR}, $alert_str); 
        }
      }
    }
  }

  return $ret_prevents_annot_flag;
}

#################################################################
# Subroutine: helper_ftable_process_feature_alerts()
# Incept:     EPN, Thu Nov  1 12:10:34 2018
#
# Purpose:    Given a string of alerts that correspond to a specific
#             sequence and feature, use the %{$alt_info_HHR} and
#             process that string to determine what (if any) 
#             alert/error messages should be added to the feature table
#             for this seq/feature pair by adding a string to 
#             @{$ret_alert_AR}. 
# 
#             As of v1.1.3 this subroutine also takes into account
#             if a feature/alert combination has the 'misc_not_failure'
#             attribute. If so, we do not add it to the @{$ret_alert_AR}
#             array, but the $ret_misc_flag (second return value)
#             will be 1. See 'Returns' for details.
#
#
# Arguments:
#   $alt_code_str:           string of errors, comma separated, can be ""
#   $seq_name:               name of sequence
#   $ftr_idx:                feature index
#   $ftr_info_AHR:           REF to array of hashes with information on the features, PRE-FILLED
#   $alt_info_HHR:           REF to hash of hashes with information on the errors, PRE-FILLED
#   $alt_ftr_instances_HHHR: REF to hash of array of hashes with per-feature errors, PRE-FILLED
#   $ret_alert_AR:           REF to array of errors, possibly added to here (not created)
#   $FH_HR:                  REF to hash of file handles, including "log" and "cmd"
# 
# Returns: Two values:
#   $ret_fatal_flag: '1' if any of the alerts in $alt_code_str are fatal (and do not have 'misc_not_failure' attribute)
#   $ret_misc_flag:  '1' if any of the alerts in $alt_code_str have 'misc_not_failure' attribute
#
# Dies: Never
#################################################################
sub helper_ftable_process_feature_alerts { 
  my $sub_name = "helper_ftable_process_feature_alerts";
  my $nargs_expected = 8;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
 
  my ($alt_code_str, $seq_name, $ftr_idx, $ftr_info_AHR, $alt_info_HHR, $alt_ftr_instances_HHHR, $ret_alert_AR, $FH_HR) = (@_);

  my $ret_fatal_flag = 0; # will be '1' if any of the alerts in $alt_code_str are fatal and don't have 'misc_not_failure' attribute
  my $ret_misc_flag  = 0; # will be '1' if any of the alerts in $alt_code_str have 'misc_not_failure' attribute

  if($alt_code_str eq "") { 
    return ($ret_fatal_flag, $ret_misc_flag); 
  }

  # printf("in $sub_name $seq_name $ftr_idx, $alt_code_str\n");

  # create a hash of all alerts in the input $alt_str
  my %input_alt_code_H = (); # $input_slt_code_H{$alt_code} = 1 if $alt_code is in $alt_code_str
  my $alt_code; 
  
  foreach $alt_code (split(",", $alt_code_str)) { 
    if(! defined $alt_info_HHR->{$alt_code}) { 
      ofile_FAIL("ERROR in $sub_name, input error of $alt_code in string $alt_code_str is invalid", 1, $FH_HR);
    }
    $input_alt_code_H{$alt_code} = 1; 
  }

  my $is_fatal  = 0; # '1' if this alert is fatal
  my $is_misc   = 0; # '1' if this feature/alert has 'misc_not_failure' attribute
  my $do_report = 0; # '1' if we should report this alert in the feature table, '0' if not, we don't report all fatal alerts, some we skip to avoid duplicates
  foreach $alt_code (sort keys (%input_alt_code_H)) { 
    $is_fatal = vdr_FeatureAlertCausesFailure($ftr_info_AHR, $alt_info_HHR, $ftr_idx, $alt_code) ? 1 : 0;
    $is_misc  = vdr_FeatureAlertIsMiscNotFailure($ftr_info_AHR, $alt_info_HHR, $ftr_idx, $alt_code) ? 1 : 0;
    if($is_fatal) { $ret_fatal_flag = 1; }
    if($is_misc)  { $ret_misc_flag  = 1; }
    my $do_report = $is_fatal; 
    # check if this alert is invalidated by another we will also report
    if(($do_report) && ($alt_info_HHR->{$alt_code}{"ftbl_invalid_by"} ne "")) { 
      my @invalid_by_alt_code_A = split(",", $alt_info_HHR->{$alt_code}{"ftbl_invalid_by"});
      foreach my $alt_code2 (@invalid_by_alt_code_A) {
        if(($alt_info_HHR->{$alt_code2}{"causes_failure"}) && 
           (exists $input_alt_code_H{$alt_code2})) { 
          $do_report = 0; # $alt_code is invalidated by $alt_code2, $alt_code2 causes failure and is also present in $alt_code_str
          # printf("\t\t\tinvalidated by $alt_code2\n");
        }
      }
    }
    if($do_report) { 
      # we could have more than one instance of this sequence/feature/alert trio
      my @instance_str_A = split(":VADRSEP:", $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}{$alt_code});
      foreach my $instance_str (@instance_str_A) { 
        my ($instance_scoords, $instance_mcoords, $instance_detail) = alert_instance_parse($instance_str);
        my $alert_str = sprintf("%s: (%s:%s) S:%s; M:%s; %s %s", 
                                $alt_info_HHR->{$alt_code}{"sdesc"}, 
                                $ftr_info_AHR->[$ftr_idx]{"type"}, 
                                $ftr_info_AHR->[$ftr_idx]{"outname"}, 
                                ($instance_scoords eq "VADRNULL") ? "-" : $instance_scoords, 
                                ($instance_mcoords eq "VADRNULL") ? "-" : $instance_mcoords, 
                                $alt_info_HHR->{$alt_code}{"ldesc"}, 
                                ($instance_detail ne "VADRNULL") ? "[" . $instance_detail . "]" : "[-]");
        # only add the alert, if an identical alert does not already exist in @{$ret_alert_AR}
        my $idx = utl_AFindNonNumericValue($ret_alert_AR, $alert_str, $FH_HR);
        if($idx == -1) { 
          push(@{$ret_alert_AR}, $alert_str); 
        }
      }
    }
  }

  return ($ret_fatal_flag, $ret_misc_flag);
}

#################################################################
#
# Other output-related subroutines:
# helper_output_sequence_alert_strings
# helper_output_feature_alert_strings
# output_alignments()
# msa_replace_sequences()
# 
#################################################################

#################################################################
# Subroutine:  helper_output_sequence_alert_strings()
# Incept:      EPN, Thu Apr  4 12:37:59 2019
#
# Purpose:    Given a sequence name, construct a string of all
#             per-sequence alerts for that sequence and return it. 
#
# Arguments: 
#  $seq_name:              name of sequence
#  $include_sdesc:         '1' to return a string with "sdesc", '0' not to
#  $alt_info_HHR:          REF to the alert info hash of arrays, PRE-FILLED
#  $alt_code_AR:           ref to alert codes in order to check
#  $alt_seq_instances_HHR: REF to array of 2D hashes with per-feature alerts, PRE-FILLED
#  $FH_HR:                 REF to hash of file handles
#
# Returns:    A string with all per-sequence alt codes for this sequence 
#             concatenated and separated by commas.
#
################################################################# 
sub helper_output_sequence_alert_strings { 
  my $sub_name = "helper_output_sequence_alert_strings";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $include_sdesc, $alt_info_HHR, $alt_code_AR, $alt_seq_instances_HHR, $FH_HR) = @_;

  my $ret_alt_str = "";
  if(defined $alt_seq_instances_HHR->{$seq_name}) { 
    foreach my $alt_code (@{$alt_code_AR}) { 
      if(defined $alt_seq_instances_HHR->{$seq_name}{$alt_code}) { 
        if($ret_alt_str ne "") { $ret_alt_str .= ","; }
        $ret_alt_str .= ($include_sdesc) ? 
            sprintf("%s(%s)", helper_tabular_replace_spaces($alt_info_HHR->{$alt_code}{"sdesc"}), $alt_code) : 
            $alt_code;
      }
    }
  }
  
  return $ret_alt_str;
}

#################################################################
# Subroutine:  helper_output_feature_alert_strings()
# Incept:      EPN, Thu Apr  4 12:31:51 2019
#
# Purpose:    Given a sequence name and feature index, construct
#             a string of all per-feature alerts for that 
#             sequence/feature pair and return it. 
#
# Arguments: 
#  $seq_name:               name of sequence
#  $ftr_idx:                feature index
#  $include_sdesc:          '1' to return a string with "sdesc", '0' not to
#  $alt_info_HHR:           REF to the alert info hash of arrays, PRE-FILLED
#  $alt_code_AR:            ref to alert codes in order to check
#  $alt_ftr_instances_HHHR: REF to array of 2D hashes with per-feature alerts, PRE-FILLED
#  $FH_HR:                  REF to hash of file handles
#
# Returns:    A string with all per-feature alt codes for this 
#             sequence/feature combo concatenated and 
#             separated by commas.
#
################################################################# 
sub helper_output_feature_alert_strings { 
  my $sub_name = "helper_output_feature_alert_strings";
  my $nargs_exp = 7;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $ftr_idx, $include_sdesc, $alt_info_HHR, $alt_code_AR, $alt_ftr_instances_HHHR, $FH_HR) = @_;

  my $ret_alt_str = "";
  if((defined $alt_ftr_instances_HHHR->{$seq_name}) && 
     (defined $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx})) { 
    foreach my $alt_code (@{$alt_code_AR}) { 
      if(defined $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}{$alt_code}) { 
        if($ret_alt_str ne "") { $ret_alt_str .= ","; }
        $ret_alt_str .= ($include_sdesc) ? 
            sprintf("%s(%s)", helper_tabular_replace_spaces($alt_info_HHR->{$alt_code}{"sdesc"}), $alt_code) : 
            $alt_code;
      }
    }
  }
  
  return $ret_alt_str;
}

#################################################################
# Subroutine:  output_alignments()
# Incept:      EPN, Sun Apr 19 08:16:36 2020
#
# Purpose:    Merge and output alignments in stockholm, aligned fasta
#             or both for a single model. The logic differs
#             significantly depending on whether -r is used or not.
#             If -r is not enabled: 
#             
#             If -r is enabled: --out_stk and --out_afa need to have
#             original sequences in them, not the replaced sequences 
#             which cmalign aligned, so we replace back to the originals
#             by fetching from {$$in_sqfile_R}. And --out_rpstk and 
#             --out_rpafa need to have replaced sequences in them.
# 
#             Regarding efficiency: if both --out_stk and --out_afa 
#             are used, we could get away with only calling 
#             msa_replace_sequences() once but we call it twice 
#             because it would complicate this already complicated code logic.
#             (And that situation should be rare, a savvy user could 
#             use esl-reformat on one or the other.)
#
# Arguments: 
#  $execs_HR:          REF to a hash with "blastx" and "parse_blastx.pl""
#  $in_sqfile_R:       REF to Bio::Easel::SqFile object from input fasta file
#  $stk_file_AR:       REF to array of stockholm files to merge to get full alignment
#  $mdl_name:          name of model this alignment is to
#  $rpn_output_HHR:    REF to 2D hash of -r related results to output, used to 
#                      determine which sequences had Ns replaced in them, will be undef unless -r
#  $out_root:          root name for output file names
#  $to_remove_AR:      REF to array of files to eventually remove, possibly ADDED TO HERE 
#  $opt_HHR:           REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:    REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If esl-alimerge fails.
#
################################################################# 
sub output_alignments { 
  my $sub_name = "output_alignments";
  my $nargs_exp = 9;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $in_sqfile_R, $stk_file_AR, $mdl_name, $rpn_output_HHR, $out_root, $to_remove_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  my $do_out_stk   = opt_Get("--out_stk",   $opt_HHR);
  my $do_out_afa   = opt_Get("--out_afa",   $opt_HHR);
  my $do_out_rpstk = opt_Get("--out_rpstk", $opt_HHR);
  my $do_out_rpafa = opt_Get("--out_rpafa", $opt_HHR);
  my $do_keep      = opt_Get("--keep", $opt_HHR);
  if($do_keep) { 
    $do_out_stk = 1;
    $do_out_afa = 1;
    if(opt_Get("-r", $opt_HHR)) { 
      $do_out_rpstk = 1;
      $do_out_rpafa = 1;
    }      
  }

  my $stk_list_file = $out_root . "." . $mdl_name . ".align.stk.list";
  utl_AToFile($stk_file_AR, $stk_list_file, 1, $FH_HR);
      
  if(! opt_Get("-r", $opt_HHR)) { 
    # default, no replacements happened
    if($do_out_stk) { 
      my $out_rfrna_stk_file = $out_root . "." . $mdl_name . ".rfrna.align.stk";
      sqf_EslAlimergeListRun($execs_H{"esl-alimerge"}, $stk_list_file, "", $out_rfrna_stk_file, "stockholm", $opt_HHR, $FH_HR);
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . ".rfrna.align.stk", $out_rfrna_stk_file, 0, $do_keep, sprintf("model $mdl_name full sequence alignment with original RF line (stockholm)"));
      if(! $do_keep) { push(@{$to_remove_AR}, $out_rfrna_stk_file); }
      # for stockholm we need to replace RNA RF with DNA
      my $msa = Bio::Easel::MSA->new({
        fileLocation => $out_rfrna_stk_file,
        forceText => 1});  
      my $rna_rf = $msa->get_rf();
      seq_SqstringDnaize(\$rna_rf);
      $msa->set_rf($rna_rf);
      $msa->addGC_rf_column_numbers();
      $msa->capitalize_based_on_rf();

      my $out_stk_file = $out_root . "." . $mdl_name . ".align.stk";
      $msa->write_msa($out_stk_file, "stockholm", 0); # 0: do not append to file if it exists
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . ".align.stk", $out_stk_file, 1, 1, sprintf("model $mdl_name full sequence alignment (stockholm)"));
    }
    if($do_out_afa) { 
      my $out_afa_file = $out_root . "." . $mdl_name . ".align.afa";
      sqf_EslAlimergeListRun($execs_H{"esl-alimerge"}, $stk_list_file, "", $out_afa_file, "afa", $opt_HHR, $FH_HR);
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . ".align.afa", $out_afa_file, 1, 1, sprintf("model $mdl_name full sequence alignment (afa)"));
      # for afa, no RF line so don't need to replace with DNA
    }
  }
  else { 
    # -r enabled
    if(($do_out_stk) || ($do_out_rpstk)) { 
      my $out_rfrna_rpstk_file = $out_root . "." . $mdl_name . ".rfrna.align.rpstk";
      my $out_stk_file         = $out_root . "." . $mdl_name . ".align.stk";
      sqf_EslAlimergeListRun($execs_H{"esl-alimerge"}, $stk_list_file, "--dna", $out_rfrna_rpstk_file, "stockholm", $opt_HHR, $FH_HR);
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . "rfrna.align.rpstk", $out_rfrna_rpstk_file, 0, $do_keep, sprintf("model $mdl_name full replaced sequence alignment with original RF line (stockholm)"));
      if(! $do_keep) { push(@{$to_remove_AR}, $out_rfrna_rpstk_file); }
      # for stockholm we need to replace RNA RF with DNA
      my $msa = Bio::Easel::MSA->new({
        fileLocation => $out_rfrna_rpstk_file,
        forceText => 1});  
      my $rna_rf = $msa->get_rf();
      seq_SqstringDnaize(\$rna_rf);
      $msa->set_rf($rna_rf);
      $msa->addGC_rf_column_numbers();
      $msa->capitalize_based_on_rf();

      my $out_rpstk_file = $out_root . "." . $mdl_name . ".align.rpstk";
      $msa->write_msa($out_rpstk_file, "stockholm", 0); # 0: do not append to file if it exists
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . ".align.rpstk", $out_rpstk_file, $do_out_rpstk, $do_out_rpstk, sprintf("model $mdl_name replaced sequence alignment (stockholm)"));

      if($do_out_stk) { 
        # swap replaced sequences back with original sequences in the alignment
        msa_replace_sequences($execs_HR, $out_rpstk_file, $out_stk_file, $in_sqfile_R, $rpn_output_HHR, $mdl_name, 
                              "stockholm", "stockholm", $to_remove_AR, $opt_HHR, $ofile_info_HHR);
        ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . "align.stk", $out_stk_file, 1, 1, sprintf("model $mdl_name full original sequence alignment (stockholm)"));
      }
      if(! $do_out_rpstk) { push(@{$to_remove_AR}, $out_rpstk_file); }
    }
    if(($do_out_afa) || ($do_out_rpafa)) { 
      my $out_rpafa_file = $out_root . "." . $mdl_name . ".align.rpafa";
      my $out_afa_file   = $out_root . "." . $mdl_name . ".align.afa";
      sqf_EslAlimergeListRun($execs_H{"esl-alimerge"}, $stk_list_file, "--dna", $out_rpafa_file, "afa", $opt_HHR, $FH_HR);
      ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . "align.rpafa", $out_rpafa_file, $do_out_rpafa, $do_out_rpafa, sprintf("model $mdl_name full replaced sequence alignment (afa)"));
      # for afa, no RF line so don't need to replace with DNA
      if($do_out_afa) { 
        # swap replaced sequences back with original sequences in the alignment
        msa_replace_sequences($execs_HR, $out_rpafa_file, $out_afa_file, $in_sqfile_R, $rpn_output_HHR, $mdl_name,
                              "afa", "afa", $to_remove_AR, $opt_HHR, $ofile_info_HHR);
        ofile_AddClosedFileToOutputInfo(\%ofile_info_HH, $mdl_name . "align.afa", $out_afa_file, 1, 1, sprintf("model $mdl_name full original sequence alignment (afa)"));
      }
      if(! $do_out_rpafa) { push(@{$to_remove_AR}, $out_rpafa_file); }
    }
  } # end of block for -r

  push(@{$to_remove_AR}, $stk_list_file);

  return;
}

#################################################################
# Subroutine:  msa_replace_sequences()
# Incept:      EPN, Sun Apr 19 07:50:59 2020
#
# Purpose:    Given an alignment file with a set of sequences 
#             and a sequence file including sequences of the same
#             names and lengths as those in the alignment file, 
#             replace the aligned sequences in the alignment file
#             with the sequences from the sequence file and output
#             the new alignment.
#
# Arguments: 
#  $execs_HR:       REF to a hash with "blastx" and "parse_blastx.pl""
#  $aln_file:       name of alignment file to replace seqs in
#  $out_aln_file:   name of alignment file to create with replaced seqs
#  $in_sqfile_R:    REF to Bio::Easel::SqFile object from input fasta file
#  $rpn_output_HHR: REF to 2D hash of -r related results to output, used to 
#                   determine which sequences had Ns replaced in them
#  $mdl_name:       name of model these seqs were aligned to
#  $informat:       input format, must be "stockholm" or "afa"
#  $outformat:      output format, must be "stockholm" or "afa"
#  $to_remove_AR:   REF to array of files to eventually remove, possibly ADDED TO HERE 
#  $opt_HHR:        REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR: REF to 2D hash of output file information, ADDED TO HERE
#
# Returns:    void
#
# Dies:       If problem parsing alignment or replacing sequences
#
################################################################# 
sub msa_replace_sequences { 
  my $sub_name = "msa_replace_sequences";
  my $nargs_exp = 11;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($execs_HR, $aln_file, $out_aln_file, $in_sqfile_R, $rpn_output_HHR, $mdl_name, $informat, $outformat, $to_remove_AR, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR = (defined $ofile_info_HHR->{"FH"}) ? $ofile_info_HHR->{"FH"} : undef;

  # we need $informat in the file names and keys because we may make these files twice, one from afa and once from stk
  my $tmp_pfam_aln_file     = $aln_file . "." . $informat . ".pfam";
  my $tmp_pfam_new_aln_file = $aln_file . "." . $informat . ".new.pfam";

  my $do_keep = opt_Get("--keep", $opt_HHR);

  sqf_EslReformatRun($execs_HR->{"esl-reformat"}, "-d", $aln_file, $tmp_pfam_aln_file, $informat, "pfam", $opt_HHR, $FH_HR);
  ofile_AddClosedFileToOutputInfo($ofile_info_HHR, $mdl_name . ".tmp.$informat.pfam", $tmp_pfam_aln_file, 0, $do_keep, "pfam formatted alignment of replaced (non-original) seqs for model $mdl_name");

  ofile_OpenAndAddFileToOutputInfo($ofile_info_HHR, $mdl_name . ".tmp.$informat.new.pfam", $tmp_pfam_new_aln_file, 0, $do_keep, "pfam formatted alignment of re-replaced (original) seqs for model $mdl_name");
  my $out_FH = $ofile_info_HHR->{"FH"}{"$mdl_name.tmp.$informat.new.pfam"};

  if(! $do_keep) { 
    push(@{$to_remove_AR}, $tmp_pfam_aln_file);
    push(@{$to_remove_AR}, $tmp_pfam_new_aln_file);
  }

  my $line_ctr = 0;
  my $alen = undef;
  my $uc_alnchar = undef; # an aligned uppercased nt
  my $uc_uachar  = undef; # an unaligned uppercased nt
  open(IN, $tmp_pfam_aln_file) || ofile_FileOpenFailure($tmp_pfam_aln_file, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    chomp $line;
    $line_ctr++;
    if(($line !~ m/^\#/) && ($line =~ m/^\S+\s+\S+$/)) { 
      $line =~ /^(\S+)\s+(\S+)$/;
      my ($seq_name, $aln_sqstring) = ($1, $2);
      my $rp_aln_sqstring = "";
      if((defined $rpn_output_HHR->{$seq_name}) && 
         (defined $rpn_output_HHR->{$seq_name}{"nnt_n_rp_tot"}) && 
         ($rpn_output_HHR->{$seq_name}{"nnt_n_rp_tot"} ne "-") && 
         ($rpn_output_HHR->{$seq_name}{"nnt_n_rp_tot"} > 0)) { 
        # at least 1 N replaced 
        my $in_ua_sqstring = $$in_sqfile_R->fetch_seq_to_sqstring($seq_name);
        my @in_ua_sqstring_A = split("", $in_ua_sqstring);
        my @aln_sqstring_A = split("", $aln_sqstring);
        my $uapos = 0;
        my $uachar = undef;
        if(! defined $alen) { 
          $alen = scalar(@aln_sqstring_A); 
        }
        elsif($alen != scalar(@aln_sqstring_A)) { 
          ofile_FAIL("ERROR in $sub_name, not all aligned sequences are the same length, failed on line $line_ctr:\n$line\n", 1, $FH_HR);
        }
        for(my $apos = 0; $apos < $alen; $apos++) { 
          if($aln_sqstring_A[$apos] =~ m/[A-Za-z]/) { 
            $rp_aln_sqstring .= $in_ua_sqstring_A[$uapos];
            
            # extra sanity check that would be removed if we weren't only replacing Ns
            $uc_alnchar = $aln_sqstring_A[$apos];
            $uc_uachar  = $in_ua_sqstring_A[$uapos];
            $uc_alnchar =~ tr/a-z/A-Z/;
            $uc_uachar  =~ tr/a-z/A-Z/;
            if(($uc_uachar ne "N") && ($uc_uachar ne $uc_alnchar)) { 
              ofile_FAIL(sprintf("ERROR in $sub_name, for $seq_name, replacing alignment position %d with unaligned position %d, but unaligned char is %s (not N or n) and aligned char is %s, they are expected to match", $apos+1, $uapos+1, $in_ua_sqstring_A[$uapos], $aln_sqstring_A[$apos]), 1, $FH_HR);
            }
            
            $uapos++;
          }
          else { 
            $rp_aln_sqstring .= $aln_sqstring_A[$apos];
          }
        }
        print $out_FH $seq_name . " " . $rp_aln_sqstring . "\n";
      } # end of 'if' entered if >= N was replaced for $seq_name
      else { 
        print $out_FH $line . "\n";
      }
    } # end of 'if(($line =~ m/^\#/) && ($line =~ m/^\S+\s+\S+$/))' { 
    else { 
      print $out_FH $line . "\n";
    }
  }
  close(IN);

  # convert newly created pfam file to desired output format and add numbering of RF columns if outformat is stockholm/pfam
  my $alimanip_opts = " --informat stockholm --dna --outformat $outformat";
  if(($outformat eq "stockholm") || ($outformat eq "pfam")) { 
    $alimanip_opts .= " --num-rf";
  }
  my $cmd = $execs_HR->{"esl-alimanip"} . $alimanip_opts . " " . $tmp_pfam_new_aln_file . " > " . $out_aln_file;
  utl_RunCommand($cmd, opt_Get("-v", $opt_HHR), 0, $FH_HR);

  return;
}


#################################################################
#
# Subroutines related to -r:
# parse_cdt_tblout_file_and_replace_ns()
#
#################################################################
#################################################################
# Subroutine:  parse_cdt_tblout_file_and_replace_ns()
# Incept:      EPN, Tue Apr 14 16:34:50 2020
#
# Purpose:     Parse a tblout file from the coverage determination
#              stage and greedily determine (based on higher score
#              first) determine the set of non-overlapping hits,
#              and 'missing' regions of sequence not covered by that
#              set of hits. For each missing region determine if it
#              satisfies the minimum criteria for being replaced
#              (length >= --r_minlen, fraction_ns >= --r_minfract{5,3,i}, 
#              missing length of sequence region == missing length of 
#              model region) and if so replace all Ns in that region 
#              with the expected nt at each corresponding position.
#              Then output that new sequence to a fasta file.
#
#              The three --r_minfract{5,3,i} options control fraction
#              of Ns at 5' end, 3' end and internal regions 
#              independently.
#
# Arguments: 
#  $tblout_file:           tblout file from a 'cdt' stage for a single model
#  $sqfile_R:              REF to Bio::Easel::SqFile object from main fasta file
#  $blastn_db_sqfile_R:    REF to Bio::Easel::SqFile object for blastn db 
#  $mdl_info_AHR:          REF to model info array of hashes, possibly added to here 
#  $exp_mdl_name:          name of model we expect on all lines of $indel_file
#  $mdl_idx:               index of $exp_mdl_name in $mdl_info_AHR
#  $seq_name_AR:           REF to array of sequences we want to parse indel info for
#  $seq_len_HR:            REF to hash of sequence lengths
#  $seq_replaced_HR:       REF to hash, key is sequence name, value is 1 if this seq was replaced
#  $rpn_output_HHR:        REF to 2D hash with information to output to .rpn tabular file, ADDED TO HERE
#  $out_root:              string for naming output files
#  $opt_HHR:               REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:        REF to 2D hash of output file information, ADDED TO HERE
#                         
# Returns:    Number of sequences that had Ns replaced and were output to fasta file
#
# Dies:       if unable to parse $indel_file
#
################################################################# 
sub parse_cdt_tblout_file_and_replace_ns { 
  my $sub_name = "parse_cdt_tblout_file_and_replace_ns";
  my $nargs_exp = 13;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($tblout_file, $sqfile_R, $blastn_db_sqfile_R, $mdl_info_AHR, $exp_mdl_name, $mdl_idx, 
      $seq_name_AR, $seq_len_HR, $seq_replaced_HR, $rpn_output_HHR, $out_root, $opt_HHR, $ofile_info_HHR) = @_;

  my $FH_HR  = $ofile_info_HHR->{"FH"};

  my $r_minlen_opt    = opt_Get("--r_minlen", $opt_HHR);
  my $small_value     = 0.00000001;
  my $r_minfract5_opt = opt_Get("--r_minfract5", $opt_HHR) - $small_value;
  my $r_minfract3_opt = opt_Get("--r_minfract3", $opt_HHR) - $small_value;
  my $r_minfracti_opt = opt_Get("--r_minfracti", $opt_HHR) - $small_value;
  my $do_keep         = opt_Get("--keep", $opt_HHR);
  my %tblout_coords_HAH = (); # hash of arrays of hashes 
                              # key is seq name
                              # value is array of hashes with hash keys: "seq_coords", "mdl_coords", "seq_start"
  my @processed_seq_name_A = (); # array of sequences read from the file, in order

  my $fa_FH = $FH_HR->{"rpn.sub.fa"};
  if(! defined $fa_FH) { 
    ofile_FAIL("ERROR in $sub_name, file handle for outputting fasta file with replaced sequences is undefined", 1, $FH_HR);
  }
  my $nseq_output = 0; # number of seqs written to the fasta file
  my $mdl_len = $mdl_info_AHR->[$mdl_idx]{"length"}; 

  # variables related to the model consensus sequence, 
  # these are only filled if nec (if we do a N-stretch-replacment for >= 1 seq)
  my $mdl_consensus_sqstring   = (defined $mdl_info_AHR->[$mdl_idx]{"cseq"}) ? $mdl_info_AHR->[$mdl_idx]{"cseq"} : undef;
  my @mdl_consensus_sqstring_A = (); 

  open(IN, $tblout_file) || ofile_FileOpenFailure($tblout_file, $sub_name, $!, "reading", $FH_HR);
  while(my $line = <IN>) { 
    if($line !~ m/^#/) { 
      chomp $line; 
      # example from cmsearch tblout
      #target name  accession query name   accession mdl     mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
      #------------ --------- ------------ --------- ---     -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
      #MT281530.1   -         NC_045512            - hmm         8276   27807      8232    27769      +     -    6 0.37 583.9 20047.7         0 !   -
      #
      # example from blastn-based tblout (converted) 
      #MT281530.1   -         NC_045512            -  blastn      8277   27806      8233    27768      +     -    -    -   0.0 36027.0       0.0 ?   -
      chomp $line;
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) < 18) {
        ofile_FAIL("ERROR in $sub_name, unable to parse tblout line, unexpected number (fewer than 18) of tokens:\n$line\n", 1, $FH_HR);
      }
      my ($seq_name, $mdl_name, $mdl_start, $mdl_stop, $seq_start, $seq_stop, $seq_strand) = ($el_A[0], $el_A[2], $el_A[5], $el_A[6], $el_A[7], $el_A[8], $el_A[9]);
      if(! defined $seq_len_HR->{$seq_name}) {
        ofile_FAIL("ERROR in $sub_name, unrecognized sequence $seq_name on line:\n$line\n", 1, $FH_HR);
      }          
      if($mdl_name ne $exp_mdl_name) { 
        ofile_FAIL("ERROR in $sub_name, unexpected model $mdl_name (expected $exp_mdl_name) on line:\n$line\n", 1, $FH_HR);
      }          
      # ignore hits on negative strand, this is okay because any seqs with best hit on negative strand 
      # will get revcompl alerts and *not* be annotated (aligned) anyway
      if($seq_strand eq "+") { 
        # add this hit to the growing model and seq coords strings if
        # it does not *completely* overlap with any of the segments so far added
        # Note: until version 1.2.1 we disallowed any overlap, but this prevented
        # the replacement of Ns in some regions in seqs with chance
        # similarity at sequence ends surrounding a model deletion
        # from being found and replaced correctly. See github issue #37.
        # Now we allow overlaps as long as they're not 100%, and then downstream
        # code that analyzes each missing region can handle overlaps as long 
        # as none are complete. If none are complete this guarantees that
        # if we sort by start position we will also be sorted by stop position.
        my $found_overlap = 0;
        my $ncoords = (defined $tblout_coords_HAH{$seq_name}) ? scalar(@{$tblout_coords_HAH{$seq_name}}) : 0;
        if(defined $tblout_coords_HAH{$seq_name}) { 
          for(my $i = 0; $i < $ncoords; $i++) { 
            my ($noverlap, undef) = seq_Overlap($seq_start, $seq_stop, $tblout_coords_HAH{$seq_name}[$i]{"seq_start"}, $tblout_coords_HAH{$seq_name}[$i]{"seq_stop"}, $FH_HR);  
            my $min_length = utl_Min((abs($seq_stop-$seq_start)+1), (abs($tblout_coords_HAH{$seq_name}[$i]{"seq_stop"} - $tblout_coords_HAH{$seq_name}[$i]{"seq_start"})+1));
            if($noverlap >= $min_length) { # actually max this can be should be $min_length
              $found_overlap = 1;
              $i = $ncoords; # breaks loop
            }
          }
        }
        else { 
          @{$tblout_coords_HAH{$seq_name}} = ();
          push(@processed_seq_name_A, $seq_name);
        }
        if(! $found_overlap) {
          %{$tblout_coords_HAH{$seq_name}[$ncoords]} = ();
          $tblout_coords_HAH{$seq_name}[$ncoords]{"seq_start"} = $seq_start;
          $tblout_coords_HAH{$seq_name}[$ncoords]{"seq_stop"}  = $seq_stop;
          $tblout_coords_HAH{$seq_name}[$ncoords]{"mdl_start"} = $mdl_start;
          $tblout_coords_HAH{$seq_name}[$ncoords]{"mdl_stop"}  = $mdl_stop;
          #printf("added S:$seq_start..$seq_stop M:$mdl_start..$mdl_stop\n");
        }
      } # end of 'if($seq_strand eq "+")'
    } # end of 'if($line !~ m/^#/)'
  } # end of 'while(my $line = <IN>)'
  close(IN);

  # for each sequence, determine the regions that are not covered by the set of nonoverlapping hits
  foreach my $seq_name (@processed_seq_name_A) { 
    my $seq_len = $seq_len_HR->{$seq_name};
    my @cur_seq_tblout_coords_AH = @{$tblout_coords_HAH{$seq_name}};
    @cur_seq_tblout_coords_AH = sort { 
      $a->{"seq_start"} <=> $b->{"seq_start"} 
    } @cur_seq_tblout_coords_AH;

    # initialize rpn_output_HHR for this sequence, to output later to .rpn file in output_tabular
    $rpn_output_HHR->{$seq_name}{"nnt_n_tot"}      = 0;  # total number of Ns
    $rpn_output_HHR->{$seq_name}{"nnt_n_rp_tot"}   = 0;  # total number of Ns replaced
    $rpn_output_HHR->{$seq_name}{"nnt_n_rp_fract"} = 0;  # fraction of Ns that are replaced
    $rpn_output_HHR->{$seq_name}{"ngaps_tot"}      = 0;  # total number of missing regions
    $rpn_output_HHR->{$seq_name}{"ngaps_int"}      = 0;  # number of internal missing regions
    $rpn_output_HHR->{$seq_name}{"ngaps_rp"}       = 0;  # number of regions in which at least 1 N is replaced
    $rpn_output_HHR->{$seq_name}{"ngaps_rp_full"}  = 0;  # number of regions in which all Ns are replaced
    $rpn_output_HHR->{$seq_name}{"ngaps_rp_part"}  = 0;  # number of regions in which not all Ns are replaced
    $rpn_output_HHR->{$seq_name}{"nnt_rp_full"}    = 0;  # number of N nts replaced in regions in which all Ns are replaced
    $rpn_output_HHR->{$seq_name}{"nnt_rp_part"}    = 0;  # number of N nts replaced in regions in which all Ns are replaced
    $rpn_output_HHR->{$seq_name}{"coords"}         = ""; # pseudo-coordinate string describing number of Ns replaced per region

    # get start and stop arrays for all seq and mdl coords (remember all strands are +)
    my $ncoords = scalar(@cur_seq_tblout_coords_AH);
    my @seq_start_A = ();
    my @mdl_start_A = ();
    my @seq_stop_A = ();
    my @mdl_stop_A = ();
    my $i;
    for($i = 0; $i < $ncoords; $i++) { 
      $seq_start_A[$i] = $cur_seq_tblout_coords_AH[$i]{"seq_start"};
      $seq_stop_A[$i]  = $cur_seq_tblout_coords_AH[$i]{"seq_stop"};
      $mdl_start_A[$i] = $cur_seq_tblout_coords_AH[$i]{"mdl_start"};
      $mdl_stop_A[$i]  = $cur_seq_tblout_coords_AH[$i]{"mdl_stop"};
      #printf("set seq_stop_A[$i] to $seq_stop_A[$i]\n");
    }

    # determine missing regions
    my @missing_seq_start_A = ();
    my @missing_seq_stop_A  = ();
    my @missing_mdl_start_A = ();
    my @missing_mdl_stop_A  = ();
    # flags used only to making sure $rpn_output_HHR->{$seq_name}{ngaps_tot} is accurate
    my $too_many_nt_5p_flag = 0; # set to '1' if missing region on 5' end extends past end of model (too many nts on 5' end)
    my $too_many_nt_3p_flag = 0; # set to '1' if missing region on 3' end extends past end of model (too many nts on 5' end)
    # check for missing sequence before first aligned region, infer first model position
    if($seq_start_A[0] != 1) { 
      # printf("$seq_name %10d..%10d is not covered\n", 1, $seq_start_A[0]-1);
      my $missing_seq_len = ($seq_start_A[0]-1) - 1 + 1;
      my $cur_missing_mdl_start = (($mdl_start_A[0]-1) - $missing_seq_len) + 1;
      # only add this missing region if it doesn't extend past end of model
      if($cur_missing_mdl_start >= 1) { 
        push(@missing_seq_start_A, 1);
        push(@missing_seq_stop_A,  $seq_start_A[0]-1);
        push(@missing_mdl_start_A, $cur_missing_mdl_start);
        push(@missing_mdl_stop_A, $mdl_start_A[0]-1);
      }
      else {
        $too_many_nt_5p_flag = 1;
      }
    }
    # check for missing sequence in between each aligned region
    for($i = 0; $i < ($ncoords-1); $i++) { 
      #printf("$seq_name %10d..%10d is not covered (mdl: %10d..%10d)\n", $seq_stop_A[$i]+1, $seq_start_A[($i+1)]-1, $mdl_stop_A[$i]+1, $mdl_start_A[($i+1)]-1);
      push(@missing_seq_start_A, $seq_stop_A[$i]+1);
      push(@missing_seq_stop_A,  $seq_start_A[($i+1)]-1);
      push(@missing_mdl_start_A, $mdl_stop_A[$i]+1);
      push(@missing_mdl_stop_A,  $mdl_start_A[($i+1)]-1);
      $rpn_output_HHR->{$seq_name}{"ngaps_int"}++;
    }
    # check for missing sequence after final aligned region, 
    # infer final model position, if it's longer than our model then 
    # the region is not the correct length so we don't attempt to 
    # replace this region. An alternative would be to replace to 
    # the end of the model, but I think that's too aggressive.
    if($seq_stop_A[($ncoords-1)] != $seq_len) { 
      #printf("$seq_name %10d..%10d is not covered\n", $seq_stop_A[($ncoords-1)], $seq_len);
      my $missing_seq_len = $seq_len - ($seq_stop_A[($ncoords-1)]+1) + 1;
      my $cur_missing_mdl_stop = ($mdl_stop_A[$i]+1) + ($missing_seq_len - 1);
      #printf("seq_stop_A[(ncoords-1)] +1 : " . ($seq_stop_A[($ncoords-1)]+1) . "\n");
      #printf("missing_seq_len:      $missing_seq_len\n");
      #printf("mdl_stop_A[$i]:      " . $mdl_stop_A[$i] . "\n");
      #printf("cur_missing_mdl_stop: $cur_missing_mdl_stop\n");
      #printf("mdl_len:              $mdl_len\n");
      if($cur_missing_mdl_stop <= $mdl_len) { 
        # only add this missing region if it doesn't extend past end of model
        push(@missing_seq_start_A, $seq_stop_A[($ncoords-1)]+1);
        push(@missing_seq_stop_A,  $seq_len);
        push(@missing_mdl_start_A, $mdl_stop_A[$i]+1);
        push(@missing_mdl_stop_A,  $cur_missing_mdl_stop);
      }
      else {
        $too_many_nt_3p_flag = 1;
      }
    }
    my $nmissing = scalar(@missing_seq_start_A);
    $rpn_output_HHR->{$seq_name}{"ngaps_tot"} = $nmissing;
    if($too_many_nt_5p_flag) { $rpn_output_HHR->{$seq_name}{"ngaps_tot"}++; }
    if($too_many_nt_3p_flag) { $rpn_output_HHR->{$seq_name}{"ngaps_tot"}++; }
    
    # first pass through all missing regions to determine if any should be replaced
    # because they meet minimum replacement thresholds:
    # - length of sequence region and model region must be identical
    #   (otherwise we wouldn't know what nt to replace Ns with)
    # - length of sequence region is at or above minimum from --r_minlen
    # - fraction of Ns in sequence region is at or above minimum from --r_minfract
    my $replaced_sqstring = "";
    my $original_seq_start = 1; # updated to position after last replaced string when we do a replacement
    my $nreplaced_regions = 0;
    my $nreplaced_nts = 0; # number of N nts replaced
    my $n_tot         = 0; # total number of Ns in the sequence
    my $seq_desc      = "";    # fetched sequence description, if any
    if($nmissing > 0) { # at least one missing region
      my $fasta_seq = $$sqfile_R->fetch_seq_to_fasta_string($seq_name, -1); # -1 puts entire sequence into second line of $fasta_sqstring
      my $fetched_seq_name = undef; # name of fetched sequence, should eq $seq_name
      my $sqstring         = "";    # fetched sqstring
      if($fasta_seq =~ /^>(\S+)(\s*[^\n]*)\n(\S+)\n$/) { 
        ($fetched_seq_name, $seq_desc, $sqstring) = ($1, $2, $3);
        # sanity check
        if($fetched_seq_name ne $seq_name) { 
          ofile_FAIL("ERROR in $sub_name, tried to fetch sequence $seq_name but fetched $fetched_seq_name", 1, $FH_HR); 
        }
      }
      else { 
        ofile_FAIL("ERROR in $sub_name, unable to parse fetched sequence fasta:\n$fasta_seq\n", 1, $FH_HR);
      }
      for($i = 0; $i < $nmissing; $i++) {
        my $missing_seq_len = $missing_seq_stop_A[$i] - $missing_seq_start_A[$i] + 1;
        my $missing_mdl_len = $missing_mdl_stop_A[$i] - $missing_mdl_start_A[$i] + 1;
        my $cur_r_minfract_opt = $r_minfracti_opt; # set to 5' or 3' threshold below if nec
        if($missing_seq_start_A[$i] == 1) { 
          $cur_r_minfract_opt = $r_minfract5_opt;
        }
        if($missing_seq_stop_A[$i] == $seq_len) { 
          $cur_r_minfract_opt = $r_minfract3_opt;
        }
        if(($missing_seq_len == $missing_mdl_len) && ($missing_seq_len >= $r_minlen_opt)) { 
          my $missing_sqstring = substr($sqstring, ($missing_seq_start_A[$i]-1), $missing_seq_len);
          $missing_sqstring =~ tr/[a-z]/[A-Z]/; # uppercaseize
          my $count_n = $missing_sqstring =~ tr/N//;
          my $fract_n = $count_n / $missing_seq_len;
          if($fract_n >= $cur_r_minfract_opt) { 
            # replace Ns in this region with expected nt
            # 
            # get the model consensus sequence if we don't have it already
            $rpn_output_HHR->{$seq_name}{"ngaps_rp"}++;
            $rpn_output_HHR->{$seq_name}{"coords"} .= "S:" . $missing_seq_start_A[$i] . ".." . $missing_seq_stop_A[$i] . ",";
            $rpn_output_HHR->{$seq_name}{"coords"} .= "M:" . $missing_mdl_start_A[$i] . ".." . $missing_mdl_stop_A[$i] . ",";
            $rpn_output_HHR->{$seq_name}{"coords"} .= "N:" . $count_n . "/" . $missing_seq_len . ";";
            if(! defined $mdl_consensus_sqstring) { 
              $mdl_info_AHR->[$mdl_idx]{"cseq"} = $$blastn_db_sqfile_R->fetch_seq_to_sqstring($exp_mdl_name);
              $mdl_consensus_sqstring = $mdl_info_AHR->[$mdl_idx]{"cseq"};
            }
            # fill in non-replaced region since previous replacement 
            # (or 5' chunk up to replacement start if this is the first replacement, 
            #  in this case $original_seq_start will be its initialized value of 1)
            if($missing_seq_start_A[$i] != 1) { # if $missing_seq_start_A[$i] is 1, there's no chunk 5' of the missing region to fetch
              $replaced_sqstring .= $$sqfile_R->fetch_subseq_to_sqstring($seq_name, $original_seq_start, $missing_seq_start_A[$i] - 1, 0); # 0: do not reverse complement
            }
            if($count_n eq $missing_seq_len) { 
              # region to replace is entirely Ns, easy case
              # replace with substr of model cseq
              $replaced_sqstring .= substr($mdl_consensus_sqstring, $missing_mdl_start_A[$i] - 1, $missing_mdl_len);
              $nreplaced_nts += $missing_seq_len;
              $rpn_output_HHR->{$seq_name}{"ngaps_rp_full"}++;
              $rpn_output_HHR->{$seq_name}{"nnt_rp_full"} += $missing_seq_len;
            }
            else { 
              # region to replace is not entirely Ns, more laborious case
              # replace only Ns with model positions
              $rpn_output_HHR->{$seq_name}{"ngaps_rp_part"}++;
              if(scalar(@mdl_consensus_sqstring_A) == 0) { # if != 0 we already have this
                @mdl_consensus_sqstring_A = split("", $mdl_consensus_sqstring); 
              }
              my @missing_sqstring_A = split("", $missing_sqstring);
              for(my $spos = 0; $spos < $missing_seq_len; $spos++) { 
                if($missing_sqstring_A[$spos] eq "N") { 
                  #printf("replacing missing_sqstring_A[$spos] with mdl_consensus_sqstring_A[%d + %d - 1 = %d] which is %s\n", $missing_mdl_start_A[$i], $spos, $missing_mdl_start_A[$i] + $spos - 1, $mdl_consensus_sqstring_A[($missing_mdl_start_A[$i] + $spos - 1)]);
                  $replaced_sqstring .= $mdl_consensus_sqstring_A[($missing_mdl_start_A[$i] + $spos - 1)];
                  $nreplaced_nts++;
                  $rpn_output_HHR->{$seq_name}{"nnt_rp_part"}++;
                }
                else { 
                  $replaced_sqstring .= $missing_sqstring_A[$spos];
                }
              }
            }
            $original_seq_start = $missing_seq_stop_A[$i] + 1;
            $nreplaced_regions++;
          } # end of 'if($fract_n >= $cur_r_minfract_opt)
        }
      } # end of 'for($i = 0; $i < nmissing; $i++);'
    } # end of 'if($nmissing > 0)'
    # if we have generated a replacement sqstring, we need to finish it off if necessary
    # with final region of the sequence after the final replaced region
    if($replaced_sqstring ne "") { 
      if($original_seq_start <= $seq_len) { 
        $replaced_sqstring .= $$sqfile_R->fetch_subseq_to_sqstring($seq_name, $original_seq_start, $seq_len, 0); # 0: do not reverse complement
      }
      if(length($replaced_sqstring) != $seq_len) { 
        ofile_FAIL(sprintf("ERROR in $sub_name, trying to replace at least one region in $seq_name, but failed, unexpected length %d should be $seq_len", length($replaced_sqstring)), 1, $FH_HR);
      }
      $n_tot  = ($replaced_sqstring =~ tr/N//);
      $n_tot += ($replaced_sqstring =~ tr/n//);
      $n_tot += $nreplaced_nts;
      $rpn_output_HHR->{$seq_name}{"nnt_n_tot"}      = $n_tot;
      $rpn_output_HHR->{$seq_name}{"nnt_n_rp_tot"}   = $nreplaced_nts;
      $rpn_output_HHR->{$seq_name}{"nnt_n_rp_fract"} = $nreplaced_nts / $n_tot;
      printf $fa_FH (">%s%s\n%s\n", $seq_name, $seq_desc, $replaced_sqstring);
      $seq_replaced_HR->{$seq_name} = 1;
      $nseq_output++;
    } # end of 'if($replaced_sqstring)'
    else { # no Ns replaced
      my $full_sqstring = $$sqfile_R->fetch_seq_to_sqstring($seq_name);
      my $n_tot  = ($full_sqstring =~ tr/N//);
      $n_tot += ($full_sqstring =~ tr/n//);
      $rpn_output_HHR->{$seq_name}{"nnt_n_tot"}      = $n_tot;
      $rpn_output_HHR->{$seq_name}{"nnt_n_rp_tot"}   = 0;
      $rpn_output_HHR->{$seq_name}{"nnt_n_rp_fract"} = 0.;
    }
  } # end of 'foreach my $seq_name'

  return $nseq_output;
}

#################################################################
#
# Miscellaneous subroutines:
# initialize_ftr_or_sgm_results()
# convert_pp_char_to_pp_avg()
# group_subgroup_string_from_classification_results()
# check_for_tabular_ftr_feature_prediction()
# check_for_valid_ftbl_feature_prediction()
# check_if_sequence_passes()
# check_if_sequence_was_annotated()
# helper_sort_hit_array()
# get_5p_most_sgm_idx_with_results()
# get_3p_most_sgm_idx_with_results()
# check_for_feature_alert_codes()
# get_accession_from_ncbi_seq_name() 
#
#################################################################

#################################################################
# Subroutine: initialize_ftr_or_sgm_results_for_model()
# Incept:     EPN, Tue Mar 15 06:01:32 2016
#
# Purpose:    Initialize the {ftr,sgm}_results_HAH data structure.
#             (for one model).
#
# Args:
#  $seq_name_AR:      REF to hash of arrays with information 
#                     on the sequences, PRE-FILLED
#  $info_AHR:         REF to array of hashes with information 
#                     on the features or segments, PRE-FILLED
#  $results_HAHR:     REF to the feature or segments results data 
#                     structure, INITIALIZED HERE
#  $FH_HR:            REF to hash of file handles
#
# Returns: void
#
#
#################################################################
sub initialize_ftr_or_sgm_results_for_model { 
  my $sub_name = "initialize_ftr_or_sgm_results_for_model()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name_AR, $info_AHR, $results_HAHR, $FH_HR) = @_;

  my $nftr_or_sgm = scalar(@{$info_AHR});

  %{$results_HAHR} = ();
  foreach my $seq_name (@{$seq_name_AR}) { 
    @{$results_HAHR->{$seq_name}} = ();
    for(my $ftr_or_sgm_idx = 0; $ftr_or_sgm_idx < $nftr_or_sgm; $ftr_or_sgm_idx++) { 
      %{$results_HAHR->{$seq_name}[$ftr_or_sgm_idx]} = ();
    }
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
  if($ppchar eq "0") { return 0.05; }

  ofile_FAIL("ERROR in $sub_name, invalid PP char: $ppchar", 1, $FH_HR); 

  return 0.; # NEVER REACHED
}

#################################################################
# Subroutine: group_subgroup_string_from_classification_results()
# Incept:     EPN, Wed Apr  3 10:25:32 2019
# Purpose:    Return a string describing a model's group
#             and subgroup 
#             probability it represents.
#
# Arguments:
#  $results_HR: hash potentially with keys "group" and "subgroup"
#               (typically 3rd dim of stg_results_HHH,
#                e.g. stg_results_HHHR->{$seq_name}{"std.cls.1"})
#             
# Returns:  "*NONE*" if "group" key not defined in %{$results_HR}
#           $results_HR->{"group"}/*NONE* if "group" key defined but "subgroup" key undef
#           $results_HR->{"group"}/$results_HR->{"subgroup"} if both "group" and "subgroup" keys defined
#
# Dies:     never
#
#################################################################
sub group_subgroup_string_from_classification_results { 
  my $sub_name = "group_subgroup_string_from_classification_results";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($results_HR) = (@_);

  if(! defined $results_HR->{"group"}) { 
    return "*NONE*";
  }
  elsif(! defined $results_HR->{"subgroup"}) {
    return $results_HR->{"group"} . "/*NONE*";
  }
  else { 
    return $results_HR->{"group"} . "/" . $results_HR->{"subgroup"};
  }

  return; # NEVER REACHED
}

#################################################################
# Subroutine: check_for_tabular_ftr_feature_prediction()
# Incept:     EPN, Wed Apr  3 13:40:42 2019
# Purpose:    Return '1' if we have a prediction for >= 1 
#             features we will output to .ftr tabular file, 
#             else return '0'.
#
# Arguments:
#  $results_HR:     hash potentially with keys "n_start", "p_qstart", "n_len";
#  $min_len:        minimum length for the feature, can be 0
#             
# Returns:  1 if a valid feature prediction exists, else 0
#
# Dies:     never
#
#################################################################
sub check_for_tabular_ftr_feature_prediction { 
  my $sub_name = "check_for_tabular_ftr_feature_prediction";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($results_HR) = (@_);

  if((defined $results_HR->{"n_start"} || 
      defined $results_HR->{"p_qstart"})) { 
    return 1;
  }

  return 0;
}

#################################################################
# Subroutine: check_for_valid_ftbl_feature_prediction()
# Incept:     EPN, Wed Apr  3 13:40:42 2019
# Purpose:    Return '1' if we have a valid prediction we could
#             output to a feature table, else return '0'.
#
# Arguments:
#  $results_HR:     hash potentially with keys "n_start", "p_qstart", "n_len";
#  $min_len:        minimum length for the feature, can be 0
#             
# Returns:  1 if a valid feature prediction exists, else 0
#
# Dies:     never
#
#################################################################
sub check_for_valid_ftbl_feature_prediction { 
  my $sub_name = "check_for_valid_ftbl_feature_prediction";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($results_HR, $min_len) = (@_);

  if((defined $results_HR->{"n_start"} || 
      defined $results_HR->{"p_qstart"}) && 
     ((! defined $results_HR->{"n_len"}) || 
      ($results_HR->{"n_len"} >= $min_len))) { 
    return 1;
  }

  return 0;
}

#################################################################
# Subroutine: check_if_sequence_passes()
# Incept:     EPN, Wed Apr  3 14:09:05 2019
# Purpose:    Check if a sequence should PASS or not based
#             on alert instances. 
#
# Arguments:
#  $seq_name:               sequence name
#  $ftr_info_AHR:           ref to array of hashes with information on the features, 
#                           can be undef if no model assigned to this seq
#  $alt_info_HHR:           ref to 2D hash of alert info
#  $alt_seq_instances_HHR:  ref to 2D hash of sequence alert instances
#  $alt_ftr_instances_HHHR: ref to 3D hash of feature alert instances
#  $FH_HR:                  ref to hash of file handles, including 'log'
#             
# Returns:  '1' if the sequence should pass, else '0'
#
# Dies:     If alt_ftr_instances_HHHR->{$seq_name} has alerts but 
#           ftr_info_AHR is undefined
#
#################################################################
sub check_if_sequence_passes { 
  my $sub_name = "check_if_sequence_passes";
  my $nargs_exp = 6;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $ftr_info_AHR, $alt_info_HHR, $alt_seq_instances_HHR, $alt_ftr_instances_HHHR, $FH_HR) = (@_);

  if((! defined $alt_seq_instances_HHR->{$seq_name}) && 
     (! defined $alt_ftr_instances_HHHR->{$seq_name})) { 
    return 1; # no alerts exist, sequence PASSes
  }
  
  my $alt_code;
  if(defined $alt_seq_instances_HHR->{$seq_name}) { 
    foreach $alt_code (keys (%{$alt_seq_instances_HHR->{$seq_name}})) { 
      if($alt_info_HHR->{$alt_code}{"causes_failure"}) { 
        return 0;  # a sequence alert that causes failure
      }
    }
  }

  my $ftr_idx;
  if(defined $alt_ftr_instances_HHHR->{$seq_name}) { 
    foreach $ftr_idx (keys (%{$alt_ftr_instances_HHHR->{$seq_name}})) { 
      foreach $alt_code (keys (%{$alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}})) { 
        if(! defined $ftr_info_AHR) { 
          # shouldn't happen caller should pass valid $ftr_info_AHR if any feature alerts exist          
          # because that means >=1 feature was annotated so a model should have been selected
          ofile_FAIL("ERROR in $sub_name, trying to check feature alert but ftr_info_AHR is undefined", 1, $FH_HR); 
        }
        if(vdr_FeatureAlertCausesFailure($ftr_info_AHR, $alt_info_HHR, $ftr_idx, $alt_code)) { 
          return 0;  # a feature alert that causes failure
        }
      }
    }
  }

  return 1; # no alerts that cause failures, sequence passes

}

#################################################################
# Subroutine: check_if_sequence_was_annotated()
# Incept:     EPN, Wed Apr  3 14:41:58 2019
# Purpose:    Check if a sequence was annotated or not.
#
# Arguments:
#  $seq_name:           sequence name
#  $cls_output_HHR:     ref to classification output
#             
# Returns:  '1' if the sequence was annotated, else '0'
#
# Dies:     never
#
#################################################################
sub check_if_sequence_was_annotated { 
  my $sub_name = "check_if_sequence_was_annotated";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name, $cls_output_HHR) = (@_);

  if((defined $cls_output_HHR->{$seq_name}) && 
     (defined $cls_output_HHR->{$seq_name}{"annot"}) && 
     ($cls_output_HHR->{$seq_name}{"annot"})) { 
    return 1; # yes, it was annotated
  }

  return 0; # no, it wasn't annotated
}

#################################################################
# Subroutine: helper_sort_hit_array()
# Incept:     EPN, Tue Apr 25 06:23:42 2017 [ribovore]
#
# Purpose:    Sort an array of regions of hits.
#
# Args:
#  $tosort_AR:   ref of array to sort, PRE-FILLED
#  $order_AR:    ref to array of original indices corresponding to @{$tosort_AR}, [1..$nhit] (not 0..$nhit-1) FILLED HERE
#  $allow_dups:  '1' to allow duplicates in $tosort_AR, '0' to not and die if
#                they're found
#  $FH_HR:       ref to hash of file handles, including "cmd"
#
# Returns:  string indicating the order of the elements in $tosort_AR in the sorted
#           array.
#
# Dies:     - if some of the regions in @{$tosort_AR} are on different strands
#             or are in the wrong format
#           - if there are duplicate values in $tosort_AR and $allow_dups is 0
#
#################################################################
sub helper_sort_hit_array { 
  my $sub_name = "helper_sort_hit_array";
  my $nargs_expected = 4;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($tosort_AR, $order_AR, $allow_dups, $FH_HR) = @_;

  my ($i, $j); # counters

  my $nel = scalar(@{$tosort_AR});

  if($nel == 1) { ofile_FAIL("ERROR in $sub_name, nel is 1 (should be > 1)", 1, $FH_HR); }

  # make a temporary hash and sort it by value, and 
  # die if we see a different strand along the way
  my %hash = ();
  my $bstrand;
  for($i = 0; $i < $nel; $i++) { 
    my($start, $stop, $strand) = vdr_CoordsSegmentParse($tosort_AR->[$i], $FH_HR);
    $hash{($i+1)} = $start . "." . $stop;
    if($i == 0) { 
      $bstrand = $strand;
    }
    if(($i > 0) && ($strand ne $bstrand)) { 
      ofile_FAIL("ERROR in $sub_name, not all regions are on same strand, region 1: $tosort_AR->[0] $bstrand, region " . $i+1 . ": $tosort_AR->[$i] $strand", 1, $FH_HR);
    }
  }
  # the <=> comparison function means sort numerically ascending
  @{$order_AR} = (sort {$hash{$a} <=> $hash{$b}} (keys %hash));

  # now that we have the sorted order, we can easily check for dups
  if(! $allow_dups) { 
    for($i = 1; $i < $nel; $i++) { 
      if($hash{$order_AR->[($i-1)]} eq $hash{$order_AR->[$i]}) { 
        ofile_FAIL("ERROR in $sub_name, duplicate values exist in the array: " . $hash{$order_AR->[$i]} . " appears twice", 1, $FH_HR); 
      }
    }
  }

  # reverse array if strand is "-"
  if($bstrand eq "-") { 
    @{$order_AR} = reverse @{$order_AR};
  }

  # construct return string
  my $ret_str = $order_AR->[0];
  for($i = 1; $i < $nel; $i++) { 
    $ret_str .= "," . $order_AR->[$i];
  }

  return $ret_str;
}

#################################################################
# Subroutine: get_5p_most_sgm_idx_with_results()
# Incept:     EPN, Mon Feb 24 15:11:47 2020
# Purpose:    Return segment index $sgm_idx of 5'-most segment for 
#             feature $ftr_idx that has results for $seq_name 
#             defined ($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"}
#
# Arguments:
#  $ftr_info_AHR:       REF to hash of arrays with information on the features, PRE-FILLED
#  $sgm_results_HAHR:   REF to results HAH, PRE-FILLED
#  $ftr_idx:            feature index
#  $seq_name:           sequence name
#             
# Returns:  index of 5'-most segment that has results for this ftr/seq
#           -1 if none
# Dies:     never
#
#################################################################
sub get_5p_most_sgm_idx_with_results { 
  my $sub_name = "get_5p_most_sgm_idx_with_results";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $sgm_results_HAHR, $ftr_idx, $seq_name) = (@_);

  for(my $sgm_idx = $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}; $sgm_idx <= $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"}; $sgm_idx++) { 
    if((defined $sgm_results_HAHR->{$seq_name}) && 
       (defined $sgm_results_HAHR->{$seq_name}[$sgm_idx]) && 
       (defined $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"})) { 
      return $sgm_idx;
    }
  }

  return -1; # none found
}

#################################################################
# Subroutine: get_3p_most_sgm_idx_with_results()
# Incept:     EPN, Mon Feb 24 15:15:38 2020
# Purpose:    Return segment index $sgm_idx of 3'-most segment for 
#             feature $ftr_idx that has results for $seq_name 
#             defined ($sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"}
#
# Arguments:
#  $ftr_info_AHR:       REF to hash of arrays with information on the features, PRE-FILLED
#  $sgm_results_HAHR:   REF to results HAH, PRE-FILLED
#  $ftr_idx:            feature index
#  $seq_name:           sequence name
#             
# Returns:  index of 5'-most segment that has results for this ftr/seq
#           -1 if none
# Dies:     never
#
#################################################################
sub get_3p_most_sgm_idx_with_results { 
  my $sub_name = "get_3p_most_sgm_idx_with_results";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_info_AHR, $sgm_results_HAHR, $ftr_idx, $seq_name) = (@_);

  # loop 3' -> 5'
  for(my $sgm_idx = $ftr_info_AHR->[$ftr_idx]{"3p_sgm_idx"}; $sgm_idx >= $ftr_info_AHR->[$ftr_idx]{"5p_sgm_idx"}; $sgm_idx--) { 
    if((defined $sgm_results_HAHR->{$seq_name}) && 
       (defined $sgm_results_HAHR->{$seq_name}[$sgm_idx]) && 
       (defined $sgm_results_HAHR->{$seq_name}[$sgm_idx]{"sstart"})) { 
      return $sgm_idx;
    }
  }

  return -1; # none found
}


#################################################################
# Subroutine: check_for_feature_alert_codes()
# Incept:     EPN, Fri Mar 27 07:01:41 2020
# Purpose:    Return '1' if at least one alert in an array
#             of alert codes exists in $alt_ftr_instances_HR,
#             else return 0.
#
# Arguments:
#  $alt_info_HHR:         REF to the alert info hash of arrays, PRE-FILLED
#  $alt_code_AR:          array of feature codes we are interested in
#  $alt_ftr_instances_HR: REF to hash of feature alert instances, key is 
#                         alert code, value is message for that alert code
#                         can be undef if not alerts exist
#             
# Returns:  '1' if any of the alert codes in @{$alt_code_AR} exist
#           as keys in %{$alt_ftr_instances_HR}, else '0'
#
# Dies:     never
#
#################################################################
sub check_for_feature_alert_codes { 
  my $sub_name = "check_for_feature_alert_codes";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($alt_info_HHR, $alt_code_AR, $alt_ftr_instances_HR) = (@_);

  if((! defined $alt_ftr_instances_HR) || 
     (! defined $alt_code_AR)) { 
    return 0;
  }

  foreach my $alt_code (@{$alt_code_AR}) { 
    if(defined $alt_ftr_instances_HR->{$alt_code}) { 
      return 1; 
    }
  }

  return 0;
}

#################################################################
# Subroutine: get_accession_from_ncbi_seq_name
# Incept:     EPN, Fri Apr  3 14:39:06 2020
# Purpose:    Given a sequence name that may or may not be 
#             in NCBI format with /^[^\|]*|(<accession>)\|.* return the
#             <accession> (where accession can have any character except |
#             If the sequence name is not in this NCBI format, replace
#             any underscores in it with _ and return that.
#             This subroutine is the main change between v1.0.5 and v1.0.6.
#
# Arguments:
#  $seq_name: sequence name to extract accession from
#             
# Returns:  <accession> extracted from $seq_name or $seq_name
#
# Dies:     never
#
#################################################################
sub get_accession_from_ncbi_seq_name { 
  my $sub_name = "accession_from_ncbi_seq_name";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name) = (@_);
  if($seq_name =~ /^[^\|]*\|([^\|]+)/) { 
    return $1; 
  }
  $seq_name =~ s/\|/\_/g;

  return $seq_name;
}

#################################################################
# Subroutine: validate_and_parse_sub_file
# Incept:     EPN, Fri Sep 25 15:06:04 2020
# Purpose:    Validate a file passed with option --msub or --xsub.
#             Each line should have 2 white-space separated tokens.
#             Both tokens are model names. Token 2 models will substitute
#             for token 1 models.
#
# Arguments:
#  $sub_file:     name of file (we already checked that it exists)
#  $mdl_info_AHR: REF to model info array of hashes, possibly added to here 
#  $sub_HR:       REF to hash of model substitutions, filled here
#  $FH_HR:        REF to hash of file handles, including "cmd"
#             
# Returns:  $err_msg, "" if no errors, if ne "", caller should exit in error
#
# Dies:     if $sub_file does not exist or is empty
#
#################################################################
sub validate_and_parse_sub_file { 
  my $sub_name = "validate_and_parse_sub_file";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sub_file, $mdl_info_AHR, $sub_HR, $FH_HR) = (@_);

  my $err_msg = "";
  my %mdl_name_H = (); # key is a model name, value is 1
  my $nmdl = scalar(@{$mdl_info_AHR});

  # create a hash of all model names, to make it easy to look them up
  for(my $mdl_idx = 0; $mdl_idx < $nmdl; $mdl_idx++) { 
    $mdl_name_H{$mdl_info_AH[$mdl_idx]{"name"}} = 1;
  }
  my @sub_line_A = (); # all lines read from $sub_file
  utl_FileLinesToArray($sub_file, 1, \@sub_line_A, $FH_HR);
  foreach my $line (@sub_line_A) { 
    if($line =~ /^(\S+)\s+(\S+)$/) { 
      my ($mdl1, $mdl2) = ($1, $2);
      if(! defined $mdl_name_H{$mdl1}) { 
        $err_msg .= "unexpected model name: $mdl1\n";
      }
      if(! defined $mdl_name_H{$mdl2}) { 
        $err_msg .= "unexpected model name: $mdl2\n";
      }
      $sub_HR->{$mdl1} = $mdl2;
    }
    else { 
      $err_msg .= "unable to parse line: $line\n";
    }
  }

  return $err_msg;
}

#################################################################
# Subroutine: write_v_annotate_scripts_for_split_mode
# Incept:     EPN, Mon Mar  8 06:42:56 2021
# Purpose:    Write one or more shell scripts that when executed
#             will run v-annotate.pl one or more times on chunks
#             of the original input fasta file.
#
# Arguments:
#  $nchunk:             number of fasta files we have created from original
#  $ncpu:               number of scripts to write
#  $in_fa_file:         main fasta file that was split up
#  $out_root_no_vadr:   string for naming output files, without 'vadr' suffix
#  $nseqs_per_chunk_AR: number of sequences per chunked fasta file, PRE-FILLED
#  $chunk_outdir_AR:    [0..$nchunk-1] REF to array of output directories created for
#                       each chunk, FILLED HERE
#  $cpu_out_file_AHR:   [0..$ncpu-1], REF to array of hashes of output files names, FILLED HERE 
#  $to_remove_AR:       REF to array of files to remove eventually, unless --keep
#  $opt_HHR:            REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR:     REF to 2D hash of output file information, ADDED TO HERE
#             
# Returns:  String that is a command to run all scripts created 
#           in this subroutine.
#
# Dies:     if unable to write scripts
#
#################################################################
sub write_v_annotate_scripts_for_split_mode { 
  my $sub_name = "write_v_annotate_scripts_for_split_mode";
  my $nargs_exp = 10;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($nchunk, $ncpu, $in_fa_file, $out_root_no_vadr, $nseqs_per_chunk_AR, $chunk_outdir_AR, $cpu_out_file_AHR, $to_remove_AR, $opt_HHR, $ofile_info_HHR) = (@_);

  my $FH_HR = $ofile_info_HHR->{"FH"};

  # determine original v-annotate.pl command and modify it as necessary for chunks
  my $v_annotate_plus_opts = get_command_and_opts($opt_HHR, $ofile_info_HHR);
  # remove --split option (must exist)
  if($v_annotate_plus_opts =~ /\s+\-\-split\s*/) { 
    $v_annotate_plus_opts =~ s/\s+\-\-split\s*/ /;
  }
  else { 
    ofile_FAIL("ERROR in $sub_name, did not find --split in original v-annotate.pl command", 1, $FH_HR);
  }
  # remove --cpu and --maxnjob option (may or may not exist)
  if($v_annotate_plus_opts =~ /\s+\-\-cpu\s+\d+\s*/) { 
    $v_annotate_plus_opts =~ s/\s+\-\-cpu\s+\d+\s*/ /;
  }
  # remove --maxnjob option (may or may not exist)
  if($v_annotate_plus_opts =~ /\s+\-\-maxnjobs\s+\d+\s*/) { 
    $v_annotate_plus_opts =~ s/\s+\-\-maxnjobs\s+\d+\s*/ /;
  }
  $v_annotate_plus_opts =~ s/\s+$//; # remove trailing whitespace if we created it

  # printf("in $sub_name, root command with opts:\n$v_annotate_plus_opts\n");

  # open all $ncpu script output files at the beginning
  my @out_FH_A         = (); # [0..$fidx..$ncpu-1] file handle for file $fidx
  my @out_scriptname_A = (); # [0..$fidx..$ncpu-1] file name for file $fidx
  my @out_ncmd_A       = (); # [0..$fidx..$ncpu-1] num commands output to file $fidx
  my $fidx;
  for($fidx = 0; $fidx < $ncpu; $fidx++) { 
    $out_scriptname_A[$fidx] = $out_root_no_vadr . ".annotate." . ($fidx+1) . ".sh";
    $out_ncmd_A[$fidx] = 0;
    open($out_FH_A[$fidx], ">", $out_scriptname_A[$fidx]) || ofile_FileOpenFailure($out_scriptname_A[$fidx], $sub_name, $!, "writing", $FH_HR);
  }

  # go through each chunk and put it in a script, alternating between all scripts
  my $sidx = 1;
  my $sidx_opt = "";
  @{$cpu_out_file_AHR} = (); # tricky: need to fill "out" in i loop over chunks and "err" in fidx loop over cpus
  for(my $i = 1; $i <= $nchunk; $i++) { 
    my $fasta_file = ($nchunk == 1) ? $in_fa_file : $in_fa_file . "." . $i;
    my $out_dir    = $out_root_no_vadr . "." . $i;
    my $out_file   = $out_root_no_vadr . "." . $i . ".out";
    $fidx = ($i-1) % $ncpu;
    $out_ncmd_A[$fidx]++;
    # determine --sidx option (note --sidx is incompatible with --split so we can assume --sidx for 
    # *this* execution of v-annotate.pl with --split has not used --sidx)
    $sidx_opt = "--sidx $sidx";
    if(! defined $cpu_out_file_AHR->[$fidx]) { 
      %{$cpu_out_file_AHR->[$fidx]} = ();
    }
    $cpu_out_file_AHR->[$fidx]{"out"} = $out_file; # may overwrite previous one
    my $FH = $out_FH_A[$fidx];
    print $FH "$v_annotate_plus_opts $sidx_opt $fasta_file $out_dir > $out_file\n";
    push(@{$chunk_outdir_AR}, $out_dir); # save chunk directory
    $sidx += $nseqs_per_chunk_AR->[($i-1)]; # update number of sequences for next command

    push(@{$to_remove_AR}, $fasta_file);
    push(@{$to_remove_AR}, $out_file);
  }

  my $script_cmd = "";
  # output first script last, so we can run it in foreground (not background)
  # we don't run final script last in foreground because it is most likely to be the 
  # smallest sequence subset

  my $err_file = undef;
  for($fidx = 1; $fidx < $ncpu; $fidx++) { 
    if($out_ncmd_A[$fidx] > 0) { 
      $err_file = $out_scriptname_A[$fidx] . ".err";
      $cpu_out_file_AHR->[$fidx]{"err"} = $err_file;
      $script_cmd .= "sh " . $out_scriptname_A[$fidx] . " > /dev/null 2> $err_file &\n"; 
    }
    push(@{$to_remove_AR}, $out_scriptname_A[$fidx]);
    push(@{$to_remove_AR}, $err_file);
    close $out_FH_A[$fidx];
  }
  # add first script at the end to run in foreground, see comment above
  if($out_ncmd_A[0] == 0) { 
    ofile_FAIL("ERROR in $sub_name, first script file $out_scriptname_A[0] has 0 commands", 1, $FH_HR);
  }
  $fidx = 0;
  $err_file = $out_scriptname_A[$fidx] . ".err";
  $cpu_out_file_AHR->[$fidx]{"err"} = $err_file;
  $script_cmd .= "sh " . $out_scriptname_A[$fidx] . " > /dev/null 2> $err_file\n"; 
  push(@{$to_remove_AR}, $out_scriptname_A[$fidx]);
  push(@{$to_remove_AR}, $err_file);
  close $out_FH_A[$fidx];

  return $script_cmd;
}

#################################################################
# Subroutine: get_command_and_opts
# Incept:     EPN, Mon Mar  8 07:02:31 2021
# Purpose:    Return a string that is the command used to execute
#             v-annotate.pl along with all command line options
#             but not input args (fasta file and output dir).
#
# Arguments:
#  $opt_HHR:        REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR: REF to 2D hash of output file information, ADDED TO HERE
#             
# Returns:  command string
#
# Dies:     never
#
#################################################################
sub get_command_and_opts { 
  my $sub_name = "get_command_and_opts";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($opt_HHR, $ofile_info_HHR) = (@_);

  my $FH_HR = $ofile_info_HHR->{"FH"};
  
  my $cmd = $0; 

  foreach my $optname ( keys %{$opt_HHR}) { 
    if(opt_IsUsed($optname, $opt_HHR)) { 
      if($cmd ne "") { $cmd .= " "; }
      $cmd .= $optname;
      if($opt_HHR->{$optname}{"type"} ne "boolean") { 
        $cmd .= " " . opt_Get($optname, $opt_HHR);
      }
    }
  }

  return $cmd;
}

#################################################################
# Subroutine: swap_gap_and_adjacent_nongap_in_rf
# Incept:     EPN, Fri Mar 12 08:06:40 2021
# Purpose:    Given an aligned msa->RF string from an MSA, swap
#             two adjacent positions in that string and return 
#             the new string. 
#             
#             Purposefully not put into Bio-Easel because manipulating
#             RF has other consequences this subroutine doesn't
#             deal with, and is only called when we know certain
#             things about the alignment, such as there is only
#             1 sequence in it, which makes manipulating the 
#             RF less fraught.
#
# Arguments:
#  $orig_rf:     original RF string
#  $orig_sscons: original SS_cons string, can be undef
#  $gap_apos:    position that is a gap in $orig_rf ('.' character)
#                that we will swap with adjacent nongap
#  $do_before:   '1' to swap with first RF position before gap
#                '0' to swap with first RF position after gap  
#             
# Returns:  Three values:
#           1. new RF string, or "" if error encountered
#           2. new SS_cons string, or "" if error encountered, or undef if $orig_sscons is undef
#           3. Error message, or "" if no error encountered
#              Errors occur if 
#                - $gap_apos is not a gap in $orig_rf
#                -    $do_before  and position before $gap_apos in $orig_rf is also a gap or doesn't exist
#                - (! $do_before) and position after  $gap_apos in $orig_rf is also a gap or doesn't exist
#
# Dies: Never
#
#################################################################
sub swap_gap_and_adjacent_nongap_in_rf { 
  my $sub_name = "swap_gap_and_adjacent_nongap_in_rf";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($orig_rf, $orig_sscons, $gap_apos, $do_before) = (@_);
  my $err_str = "";

  my $rflen = length($orig_rf);

  # printf("in $sub_name, gap_apos: $gap_apos, do_before: $do_before, rflen: $rflen\n");
  # printf("orig_rf  $orig_rf\n");

  # contract checks
  if($gap_apos < 1) { 
    return("", "", "ERROR in $sub_name, gap_apos $gap_apos is negative.\n");
  }
  if($gap_apos > $rflen) { 
    return("", "", "ERROR in $sub_name, gap_apos $gap_apos exceeds RF length $rflen.\n");
  }
  if(($gap_apos == 1) && ($do_before)) { 
    return("", "", "ERROR in $sub_name, gap_apos is 1 and do_before is 1, can't do swap.\n");
  }
  if(($gap_apos == $rflen) && (! $do_before)) { 
    return("", "", "ERROR in $sub_name, gap_apos is final RF position $rflen and do_before is 0, can't do swap.\n");
  }

  my $orig_gap_rf    = substr($orig_rf, ($gap_apos-1), 1);
  my $orig_nongap_rf = ($do_before) ? 
      substr($orig_rf, ($gap_apos-2), 1) : 
      substr($orig_rf, $gap_apos,     1);
  if($orig_gap_rf ne ".") { 
    return("", "", "ERROR in $sub_name, original RF position $gap_apos is $orig_gap_rf, expected a '.' (gap) character.\n");
  }
  if($orig_nongap_rf !~ m/\w/) { 
    if($do_before) { 
      return("", "", "ERROR in $sub_name, original RF position before $gap_apos is $orig_nongap_rf, expected a nongap character.\n");
    }
    else { 
      return("", "", "ERROR in $sub_name, original RF position after $gap_apos is $orig_nongap_rf, expected a nongap character.\n");
    }
  }

  # we are not as strict with checks to the SS_cons
  my $orig_gap_ss    = undef;
  my $orig_nongap_ss = undef;
  if(defined $orig_sscons) { 
    $orig_gap_ss    = substr($orig_sscons, ($gap_apos-1), 1);
    $orig_nongap_ss = ($do_before) ? 
      substr($orig_sscons, ($gap_apos-2), 1) : 
      substr($orig_sscons, $gap_apos,     1);
  }

  # if we get here we can do the swap
  my $ret_rf = "";
  my $ret_sscons = (defined $orig_sscons) ? "" : undef;
  if($do_before) { 
    $ret_rf .= substr($orig_rf, 0, ($gap_apos-2));
    if(defined $ret_sscons) { $ret_sscons .= substr($orig_sscons, 0, ($gap_apos-2)); } 
    # printf("ret_rf = substr(orig_rf, 0, %d)\n", ($gap_apos-2));
    # print("ret_rf 0 $ret_rf\n");

    $ret_rf .= $orig_gap_rf . $orig_nongap_rf;
    if(defined $ret_sscons) { $ret_sscons .= $orig_gap_ss . $orig_nongap_ss; } 
    # print("ret_rf .= $orig_gap_rf + $orig_nongap_rf\n");
    # print("ret_rf 1 $ret_rf\n");

    $ret_rf .= substr($orig_rf, $gap_apos);
    if(defined $ret_sscons) { $ret_sscons .= substr($orig_sscons, $gap_apos); }
    # print("ret_rf .= substr(orig_rf, $gap_apos)\n");
    # print("ret_rf 2 $ret_rf\n");
  }
  else {  # ! $do_before
    $ret_rf .= substr($orig_rf, 0, ($gap_apos-1));
    if(defined $ret_sscons) { $ret_sscons .= substr($orig_sscons, 0, ($gap_apos-1)); } 
    # printf("ret_rf = substr(orig_rf, 0, %d)\n", ($gap_apos-2));
    # print("ret_rf 0 $ret_rf\n");

    $ret_rf .= $orig_nongap_rf . $orig_gap_rf;
    if(defined $ret_sscons) { $ret_sscons .= $orig_nongap_ss . $orig_gap_ss; } 
    # print("ret_rf .= $orig_nongap_rf + $orig_gap_rf\n");
    # print("ret_rf 1 $ret_rf\n");

    $ret_rf .= substr($orig_rf, $gap_apos+1);
    if(defined $ret_sscons) { $ret_sscons .= substr($orig_sscons, $gap_apos+1); }
    # print("ret_rf .= substr(orig_rf, $gap_apos)\n");
    # print("ret_rf 2 $ret_rf\n");
  }

  if(length($ret_rf) != $rflen) { 
    return("", "", "ERROR in $sub_name, did not create RF string of correct length, should be $rflen, created " . length($ret_rf) . "\n");
  }

  return ($ret_rf, $ret_sscons, "");
}

#################################################################
# Subroutine: msa_create_rfpos_to_apos_map
# Incept:     EPN, Wed Mar 17 14:43:19 2021
# Purpose:    Fill @{$rf2a_AR} array indicating what alignment 
#             position each RF position maps to.
#
# Arguments:
#  $msa:        the alignment
#  $rf2a_AR:    [1..$rfpos..$rflen] = $apos;  rf position $rfpos maps to alignment position $apos [1..$alen]  
#               ($rf2a_A[0] = -1  (dummy value))
#               FILLED HERE
#  $FH_HR:      ref to hash of file handles
#             
# Returns:  nongap length of RF annotation
#
# Dies: If RF string does not match $msa->alen
#
#################################################################
sub msa_create_rfpos_to_apos_map {
  my $sub_name = "msa_create_rfpos_to_apos_map";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($msa, $rf2a_AR, $FH_HR) = (@_);

  my $alen = $msa->alen;

  @{$rf2a_AR} = ();
  $rf2a_AR->[0] = -1; 
  my $rf_str = $msa->get_rf;
  my @rf_A = split("", $rf_str);
  if(scalar(@rf_A) != $alen) { 
    ofile_FAIL(sprintf("ERROR in $sub_name, alignment length and RF length mismatch $alen != %d\n", scalar(@rf_A)), 1, $FH_HR);
  }
  my $rfpos = 0; # nongap RF (model) position [1..$rflen]
  my $apos  = 0; # alignment position [1..$alen]
  for($apos = 1; $apos <= $alen; $apos++) { 
    if($rf_A[($apos-1)] ne ".") { 
      # nongap RF (model position)
      $rfpos++;
      $rf2a_AR->[$rfpos] = $apos;
    }
  }
  return $rfpos;
}

#################################################################
# Subroutine: doctoring_check_new_codon_validity
# Incept:     EPN, Wed Mar 17 16:07:26 2021
# Purpose:    Update output information in %{$dcr_output_HAHR} for
#             a possible doctoring of a start codon. Check if 
#             doctoring the start codon will make it valid
#             and if so, fill information on how to do that doctoring
#             for caller in %{$seq_dcr_info_HAR}.
#
# Arguments:
#  $delete_or_insert: 'insert' or 'delete' doctoring type
#  $start_or_stop:    'start' or 'stop' for codon type being (possibly) doctored
#  $strand:            strand of segment
#  $seq_name:          name of sequence
#  $mdl_name:          name of model 
#  $mdl_tt:            translation table for model
#  $ftr_idx:           feature index the doctoring pertains to
#  $uapos:             unaligned position of first start or final stop position (pre-doctoring)
#  $rfpos:             reference position of first start or final stop position (doesn't change with doctoring)
#  $sqstring_aligned:  aligned sequence string
#  $seq_doctor_ctr:    how many times this sequence has already had >= 1 segment doctored
#  $rf2a_AR:           reference to map of reference positions to alignment positions, PRE-FILLED
#  $dcr_output_HAHR:   reference to output information for .dcr file, ADDED TO HERE
#  $opt_HHR:           reference to 2D hash of option values, see top of sqp_opts.pm for description
#  $FH_HR:             REF to hash of file handles, including "cmd"
#   
# Returns:  '1' if start or stop codon will be valid after doctoring and 
#           doctoring limit hasn't been reached (we can only doctor a 
#           sequence twice). If '1' $seq_dcr_info_HAR will have been added to
#           and caller will actually do the doctoring of the alignment.
#
# Alignment doctoring examples: 
#
#################################################
# Before doctoring:    |      After doctoring:  |
# =================    |      ----------------  |
#------------------------------------------------
# 1. delete, + strand, start                    |
#                      |                        |
#     12 34            |        1 234           |
# seq AA-TG            |    seq A-ATG           |
# RF  AAATG            |    RF  AAATG           |
#     12345            |        12345           |
# (start) uapos: 3     | (start) uapos: 2       |
# -----------------------------------------------
# indel_apos:    3 (does not change)            |
# start_rfpos:   3 (does not change)            |
# swap_before:   1 (true, does not change)      |
# -----------------------------------------------
# 2. delete, + strand, stop                     |
#                      |                        |
#     1234 56          |        12345 6         |
# seq GGTA-AC          |    seq GGTAA-C         |
# RF  GGTAAAC          |    RF  GGTAAAC         |
#     1234567          |        1234567         |
# (stop) uapos: 4      | (stop) uapos: 5        |
# -----------------------------------------------
# indel_apos:    5 (does not change)            |
# stop_rfpos:    5 (does not change)            |
# swap_before:   0 (false, does not change)     |
# -----------------------------------------------
# 3. delete, - strand, start                    |
#                      |                        |
#     1234 567         |        12345 67        |
# seq GGCA-TTT         |    seq GGCAT-TT        |
# RF  GGCATTTT         |    RF  GGCATTTT        |
#     12345678         |        12345678        |
#       gta            |          gta           |
# (start) uapos: 4     | (start) uapos: 5       |
# -----------------------------------------------
# indel_apos:    5 (does not change)            |
# start_rfpos:   5 (does not change)            |
# swap_before:   0 (false, does not change)     |
# -----------------------------------------------
# 4. delete, - strand, stop                     |
#                      |          gat           |
#     1234 567         |        1 23456         |
# seq GC-TACC          |    seq G-CTACC         |
# RF  GCTTACC          |    RF  GCTTACC         |
#     1234567          |        1234567         |
#      gaat            |          aat           |
# (stop) uapos: 3      | (stop) uapos: 2        |
# -----------------------------------------------
# indel_apos:    3 (does not change)            |
# stop_rfpos:    3 (does not change)            |
# swap_before:   1 (true, does not change)      |
# -----------------------------------------------
# ===============================================
# ===============================================
# Before doctoring:    |      After doctoring:  |
# =================    |      ----------------  |
#------------------------------------------------
# 5. insert, + strand, start                    |
#                      |                        |
#     123456           |        123456          |
# seq AAAaTG           |    seq AAaATG          |
# RF  AAA.TG           |    RF  AA.ATG          |
#     123 45           |        12 345          |
# (start) uapos: 3     | (start) uapos: 4       |
# ins_str: 3:4:1       | ins_str: 2:3:1         |
# -----------------------------------------------
# indel_apos:  4 (does not change)              |
# start_rfpos: 3 (insert after, does not change)|
# swap_before: 1 (true, does not change)        |
# -----------------------------------------------
# 6. insert, + strand, stop                     |
#                      |                        |
#     1234567          |        1234567         |
# seq GTAaACG          |    seq GTAAaCG         |
# RF  GTA.ACG          |    RF  GTAA.CG         |
#     123 456          |        1234.56         |
# (stop) uapos: 5      | (stop) uapos: 4        |
# ins_str: 3:4:1       | ins_str: 4:5:1         |
# -----------------------------------------------
# indel_apos: 5 (does not change)               |
# stop_rfpos: 5 (insert before, does not change)|
# swap_before:0 (false, does not change)        |
# -----------------------------------------------
# 7. insert, - strand, start                    |
#                      |                        |
#     12345678         |        12345678        |
# seq GGCAtTTT         |    seq GGCATtTT        |
# RF  GGCA.TTT         |    RF  GGCAT.TT        |
#     1234 567         |        12345.67        |
#       gt a           |          gta           |
# (start) uapos: 6     | (start) uapos: 5       |
# ins_str: 4:5:1       | ins_str: 5:6:1         |
# -----------------------------------------------
# indel_apos:  5 (does not change)              |
# start_rfpos: 5 (insert before, doesn't change)|
# swap_before: 0 (false, does not change)       |
# -----------------------------------------------
# 8. insert, - strand, stop                     |
#                      |           aat          |
#     12345678         |        12345678        |
# seq GGCtTACC         |    seq GGcTTACC        |
# RF  GGC.TACC         |    RF  GG.CTACC        |
#     123.4567         |        12 34567        |
#       g at           |           gat          |
# (stop) uapos: 3      | (stop) uapos: 4        |
# ins_str: 3:4:1       | ins_str: 2:3:1         |
# -----------------------------------------------
# indel_apos:  4 (does not change)              |
# stop_rfpos:  3 (insert after, doesn't change) |
# swap_before: 1 (true, does not change)        |
# -----------------------------------------------
# 
# It is possible that doctoring one start/stop will
# break another, in which case we re-doctor to fix the 
# break. This is done by the caller
# parse_stk_and_add_alignment_alerts_cds_and_mp_alerts() 
# which takes care to only allow 2 rounds of doctoring else we
# could get into an infinite loop.
# 
# Example of situation where doctoring breaks a different
# start/stop and so a second doctoring has to take place
# to undo the first:
# 
# before 1st doctoring (and after 2nd):
# seq         ACTAAA-TGTCTGA
# RF          ACTAAAATGTCTGA
#                   ^
#                   start codon
#
# after 1st doctoring, before 2nd doctoring:
# seq         ACTA-AGTGTCTGA
# RF          ACTAAAGTGTCTGA
#               ^
#               stop codon
#
#################################################################
sub doctoring_check_new_codon_validity { 
  my $sub_name = "doctoring_check_new_codon_validity";
  my $nargs_exp = 15;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my($delete_or_insert, $start_or_stop, $strand, $seq_name, $mdl_name, $mdl_tt, $ftr_idx, 
     $uapos, $rfpos, $sqstring_aligned, $seq_doctor_ctr, $rf2a_AR, $dcr_output_HAHR, $opt_HHR, $FH_HR) = (@_);

  # to check we need to get full unaligned sqstring, this is expensive, but should be rare
  my $ua_sqstring = $sqstring_aligned;
  $ua_sqstring =~ s/\W//g;

  # store information on this to dcr_output for eventual output in output_tabular()
  if(! defined $dcr_output_HAHR->{$seq_name}) { 
    @{$dcr_output_HAHR->{$seq_name}} = ();
  }
  my $ndcr = scalar(@{$dcr_output_HAHR->{$seq_name}});
  %{$dcr_output_HAHR->{$seq_name}[$ndcr]} = ();
  # fill information that is independent of start/stop, delete/insert, +/- strand
  $dcr_output_HAHR->{$seq_name}[$ndcr]{"mdl_name"}       = $mdl_name;
  $dcr_output_HAHR->{$seq_name}[$ndcr]{"ftr_idx"}        = $ftr_idx;
  $dcr_output_HAHR->{$seq_name}[$ndcr]{"dcr_type"}       = $delete_or_insert;
  $dcr_output_HAHR->{$seq_name}[$ndcr]{"rfpos"}          = $rfpos;
  $dcr_output_HAHR->{$seq_name}[$ndcr]{"orig_seq_uapos"} = $uapos;
  $dcr_output_HAHR->{$seq_name}[$ndcr]{"codon_type"}     = $start_or_stop;
  $dcr_output_HAHR->{$seq_name}[$ndcr]{"dcr_iter"}       = $seq_doctor_ctr+1;
  $dcr_output_HAHR->{$seq_name}[$ndcr]{"did_swap"}       = "no"; # possibly changed to "yes" below

  my $new_codon_is_valid = 0; # possibly set to '1' below
  my $orig_codon = undef; # needed so we can reverse complement if nec
  my $new_codon = undef;  # needed so we can reverse complement if nec
  ###################################################
  # start block
  if($start_or_stop eq "start") { 
    # fill in start-specific info that is insert/delete dependent and possibly strand-dependent
    if($delete_or_insert eq "delete") { 
      # start, delete
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"indel_apos"} = $rf2a_AR->[$rfpos];
      if($strand eq "+") { 
        # start, delete, + strand
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}     = substr($ua_sqstring, $uapos-2, 3);
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"} = $uapos-1;
      } 
      else { 
        # start, delete, - strand
        $new_codon = substr($ua_sqstring, $uapos-2, 3);
        seq_SqstringReverseComplement(\$new_codon);
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}     = $new_codon;
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"} = $uapos+1;
      }
    }
    else { 
      # start, insert
      if($strand eq "+") { 
        # start, insert, + strand
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}     = substr($ua_sqstring, $uapos, 3);
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"} = $uapos+1;
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"indel_apos"}    = $rf2a_AR->[$rfpos]+1;
      }
      else { 
        # start, insert, - strand
        $new_codon = substr($ua_sqstring, $uapos-4, 3);
        seq_SqstringReverseComplement(\$new_codon);
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}     = $new_codon;
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"} = $uapos-1;
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"indel_apos"}    = $rf2a_AR->[$rfpos]-1;
      }
    }

    # fill in start-specific info that is strand-dependent but insert/delete independent
    if($strand eq "+") { 
      # start, insert or delete, + strand
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"orig_codon"}   = substr($ua_sqstring, $uapos-1, 3);
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"codon_coords"} = vdr_CoordsSegmentCreate($dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"}, 
                                                                                     $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"}+2, 
                                                                                     "+", $FH_HR);
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"before"} = 1; 
    }
    else { # strand is -
      # start, insert or delete, - strand
      $orig_codon = substr($ua_sqstring, $uapos-3, 3);
      seq_SqstringReverseComplement(\$orig_codon);
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"orig_codon"}   = $orig_codon;
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"codon_coords"} = vdr_CoordsSegmentCreate($dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"}, 
                                                                                     $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"}-2, 
                                                                                     "-", $FH_HR);
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"before"} = 0; 
    }
    $new_codon_is_valid = sqstring_check_start($dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}, $mdl_tt, (opt_Get("--atgonly", $opt_HHR)), $FH_HR);
  }
  ###################################################
  # stop block
  else { # $start_or_stop eq "stop"
    # fill in start-specific info that is insert/delete dependent and possibly strand-dependent
    if($delete_or_insert eq "delete") {
      # stop, delete
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"indel_apos"} = $rf2a_AR->[$rfpos];
      if($strand eq "+") { 
        # stop, delete, + strand
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}     = substr($ua_sqstring, $uapos-2, 3);
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"} = $uapos+1;
      }
      else { 
        # stop, delete, - strand
        $new_codon = substr($ua_sqstring, $uapos-2, 3);
        seq_SqstringReverseComplement(\$new_codon);
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}     = $new_codon;
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"} = $uapos-1;
      }
    }
    else { 
      # stop, insert
      if($strand eq "+") { 
        # stop, insert, + strand
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}     = substr($ua_sqstring, $uapos-4, 3);
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"indel_apos"}    = $rf2a_AR->[$rfpos]-1;
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"} = $uapos-1;
      }
      else { 
        # stop, insert, - strand
        $new_codon = substr($ua_sqstring, $uapos, 3);
        seq_SqstringReverseComplement(\$new_codon);
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}     = $new_codon;
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"indel_apos"}    = $rf2a_AR->[$rfpos]+1;
        $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"} = $uapos+1;
      }
    }

    # fill in stop-specific info that is strand-dependent but insert/delete independent
    if($strand eq "+") { 
      # stop, insert or delete, + strand
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"orig_codon"}   = substr($ua_sqstring, $uapos-3, 3);
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"codon_coords"} = vdr_CoordsSegmentCreate($dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"}-2, 
                                                                                     $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"}, 
                                                                                     "+", $FH_HR);
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"before"} = 0; 
    }
    else { 
      # stop, insert or delete, - strand
      $orig_codon = substr($ua_sqstring, $uapos-1, 3);
      seq_SqstringReverseComplement(\$orig_codon);
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"orig_codon"}   = $orig_codon;
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"codon_coords"} = vdr_CoordsSegmentCreate($dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"}+2, 
                                                                                     $dcr_output_HAHR->{$seq_name}[$ndcr]{"new_seq_uapos"}, "-", $FH_HR);
      $dcr_output_HAHR->{$seq_name}[$ndcr]{"before"} = 1; 
    }
    $new_codon_is_valid = sqstring_check_stop($dcr_output_HAHR->{$seq_name}[$ndcr]{"new_codon"}, $mdl_tt, $FH_HR);
  }
  ###################################################

  if(($new_codon_is_valid) && ($seq_doctor_ctr <= 1)) { 
    $dcr_output_HAHR->{$seq_name}[$ndcr]{"did_swap"}  = "yes";
    return $ndcr;
  }

  return -1; # codon was not valid or $seq_doctor_ctr >= 2
}

#################################################################
# Subroutine: output_mdl_and_alc_files_and_remove_temp_files()
# Incept:     EPN, Mon Mar 22 16:32:35 2021
# Purpose:    Output the mdl and alc files and remove all files
#             if (@{$to_remove_A}) unless --keep. 
#
# Arguments:
#  $zero_alt:       '1' if zero alerts were output (in which case we don't output alc file)
#  $to_remove_AR:   ref to array of files to remove
#  $opt_HHR:        ref to 2D hash of option values, see top of sqp_opts.pm for description
#  $ofile_info_HHR: ref to 2D hash of output file information, added to here
#             
# Returns:  void
#
#################################################################
sub output_mdl_and_alc_files_and_remove_temp_files { 
  my $sub_name = "output_mdl_and_alc_files_and_remove_temp_files";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($zero_alt, $to_remove_AR, $opt_HHR, $ofile_info_HHR) = (@_);

  # close the two files we may output to stdout and the log
  close($ofile_info_HHR->{"FH"}{"mdl"}); 
  close($ofile_info_HHR->{"FH"}{"alc"}); 
  
  my $FH_HR  = $ofile_info_HH{"FH"};
  
  my @conclude_A = ();
  push(@conclude_A, "#");
  push(@conclude_A, "# Summary of classified sequences:");
  push(@conclude_A, "#");
  my @file_A = ();
  utl_FileLinesToArray($ofile_info_HHR->{"fullpath"}{"mdl"}, 1, \@file_A, $FH_HR);
  push(@conclude_A, @file_A);
  push(@conclude_A, "#");
  if($zero_alt) { 
    push(@conclude_A, "# Zero alerts were reported.");
  }
  else { 
    push(@conclude_A, "# Summary of reported alerts:");
    push(@conclude_A, "#");
    my @file_A = ();
    utl_FileLinesToArray($ofile_info_HHR->{"fullpath"}{"alc"}, 1, \@file_A, $FH_HR);
    push(@conclude_A, @file_A);
  }
  
  foreach my $line (@conclude_A) { 
    ofile_OutputString($FH_HR->{"log"}, 1, $line . "\n");
  }
  
# remove unwanted files, unless --keep
  if(! opt_Get("--keep", $opt_HHR)) { 
    my @to_actually_remove_A = (); # sanity check: make sure the files we're about to remove actually exist
    my %to_actually_remove_H = (); # sanity check: to make sure we don't try to delete 
    foreach my $to_remove_file (@to_remove_A) { 
      if((defined $to_remove_file) && (-e $to_remove_file) && (! defined $to_actually_remove_H{$to_remove_file})) { 
        push(@to_actually_remove_A, $to_remove_file); 
        $to_actually_remove_H{$to_remove_file} = 1; 
      }
    }
    utl_FileRemoveList(\@to_actually_remove_A, $sub_name, $opt_HHR, $FH_HR);
  }

  return;
}

#################################################################
# Subroutine: helper_tabular_fill_header_and_justification_arrays
# Incept:     EPN, Mon May 24 15:28:52 2021
# Purpose:    For a given tabular output file, fill the header and
#             left justification arrays.
#
# Arguments:
#  $ofile_key: key for tabular file in ofile_info_HHR
#  $head_AAR:  ref to 2D array of header values
#  $clj_AR:    ref to '1'/'0' array of indicating if a column is left justified or not
#  $FH_HR:     ref to hash of file handles, including "cmd"
#             
# Returns:  void
#
#################################################################
sub helper_tabular_fill_header_and_justification_arrays { 
  my $sub_name = "helper_tabular_fill_header_and_justification_arrays";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ofile_key, $head_AAR, $clj_AR, $FH_HR) = (@_);

  undef @{$head_AAR};
  undef @{$clj_AR};
  @{$head_AAR} = ();
  @{$clj_AR}   = ();
  
  if($ofile_key eq "ant") {
    @{$head_AAR->[0]} = ("seq", "seq",  "seq", "",    "",    "best",  "",    "sub", "",    "",    "",    "",    "",      "seq");
    @{$head_AAR->[1]} = ("idx", "name", "len", "p/f", "ant", "model", "grp", "grp", "nfa", "nfn", "nf5", "nf3", "nfalt", "alerts");
    @{$clj_AR}        = (1,     1,      0,     1,     1,     1,       1,     1,     0,     0,     0,     0,     0,       1);
  }
  elsif($ofile_key eq "cls") {
    @{$head_AAR->[0]} = ("seq", "seq",  "seq", "",    "",    "",       "",     "sub",  "",      "",      "seq", "mdl", "",     "num",  "",    "",       "",     "sub",  "score", "diff/", "seq");
    @{$head_AAR->[1]} = ("idx", "name", "len", "p/f", "ant", "model1", "grp1", "grp1", "score", "sc/nt", "cov", "cov", "bias", "hits", "str", "model2", "grp2", "grp2", "diff",  "nt",    "alerts");
    @{$clj_AR}        = (1,     1,      0,     1,     1,     1,        1,      1,      0,       0,       0,     0,     0,      0,      0,     1,        1,      1,      0,       0,       1);
  }
  elsif($ofile_key eq "ftr") {
    @{$head_AAR->[0]} = ("",    "seq",  "seq", "",    "",      "ftr",  "ftr",  "ftr", "ftr", "par", "",    "",       "",     "",        "",    "",     "",     "",       "",     "",        "",     "",    "",    "seq",    "model",  "ftr");
    @{$head_AAR->[1]} = ("idx", "name", "len", "p/f", "model", "type", "name", "len", "idx", "idx", "str", "n_from", "n_to", "n_instp", "trc", "5'N",  "3'N",  "p_from", "p_to", "p_instp", "p_sc", "nsa", "nsn", "coords", "coords", "alerts");
    @{$clj_AR}        = (1,     1,      0,     1,     1,       1,      1,      0,     0,     0,     0,     0,        0,      0,         1,     0,      0,      0,        0,      0,         0,      0,     0,     0,        0,        1);
  }
  elsif($ofile_key eq "sgm") {
    @{$head_AAR->[0]} = ("",    "seq",  "seq", "",    "",      "ftr",  "ftr",  "ftr", "num", "sgm", "seq",  "seq", "mdl",  "mdl", "sgm", "",    "",    "5'", "3'", "5'",  "3'");
    @{$head_AAR->[1]} = ("idx", "name", "len", "p/f", "model", "type", "name", "idx", "sgm", "idx", "from", "to",  "from", "to",  "len", "str", "trc", "pp", "pp", "gap", "gap");
    @{$clj_AR}     = (1,     1,      0,     1,     1,       1,      1,      0,     0,     0,     0,      0,     0,      0,     0,     0,     1,     0,    0,    1,     1);
  }
  elsif($ofile_key eq "alt") {
    @{$head_AAR->[0]} = ("",    "seq",  "",      "ftr",  "ftr",  "ftr", "alert", "",     "alert",       "seq",    "seq", "mdl",    "mdl", "alert");
    @{$head_AAR->[1]} = ("idx", "name", "model", "type", "name", "idx", "code",  "fail", "description", "coords", "len", "coords", "len", "detail");
    @{$clj_AR  }      = (1,     1,      1,       1,      1,      0,     1,       1,      1,          0,        0,     0,        0,     1);
  }
  elsif($ofile_key eq "alc") {
    @{$head_AAR->[0]} = ("",    "alert",  "causes",  "short",       "per",  "num",   "num",  "long");
    @{$head_AAR->[1]} = ("idx", "code",   "failure", "description", "type", "cases", "seqs", "description");
    @{$clj_AR}        = (1,     1,        1,          1,            0,      0,      0,        1);
  }
  elsif($ofile_key eq "mdl") {
    @{$head_AAR->[0]} = ("",    "",      "",      "",         "num",  "num",  "num");
    @{$head_AAR->[1]} = ("idx", "model", "group", "subgroup", "seqs", "pass", "fail");
    @{$clj_AR}        = (1,     1,       1,       1,          0,      0,      0);
  }
  elsif($ofile_key eq "dcr") {
    @{$head_AAR->[0]} = ("",    "seq",   "mdl",  "ftr",  "ftr",  "ftr",  "dcr",  "model", "indel",      "orig",       "new", "codon",  "codon",  "orig",   "new",   "dcr",   "did");
    @{$head_AAR->[1]} = ("idx", "name", "name", "type", "name",  "idx", "type",    "pos",  "apos", "seq-uapos", "seq-uapos",  "type", "coords", "codon", "codon",  "iter", "swap?");
    @{$clj_AR}        = (1,     1,      1,      1,      1,       0,     1,       0,       0,       0,           0,           1,       0,        0,       0,        0,      1);
  }
  elsif($ofile_key eq "sda") {
    @{$head_AAR->[0]} = ("seq", "seq",    "seq", "",      "",      "seed",      "seed",     "seed",     "5'unaln", "5'unaln", "5'unaln",  "3'unaln", "3'unaln", "3'unaln");
    @{$head_AAR->[1]} = ("idx", "name",   "len", "model", "p/f",   "seq",       "mdl",      "fraction", "seq",     "mdl",     "fraction", "seq",     "mdl",     "fraction");
    @{$clj_AR}        = (1,     1,        0,     1,       1,       0,           0,          0,          0,         0,         0,          0,         0,         0);
  }  
  elsif($ofile_key eq "rpn") {
    @{$head_AAR->[0]} = ("seq", "seq",    "seq", "",      "",      "num_Ns",  "num_Ns", "fract_Ns", "ngaps", "ngaps",  "ngaps",   "ngaps",   "ngaps",   "nnt",     "nnt",     "replaced_coords");
    @{$head_AAR->[1]} = ("idx", "name",   "len", "model", "p/f",   "tot",     "rp",     "rp",       "tot",   "int",    "rp",      "rp-full", "rp-part", "rp-full", "rp-part", "seq(S),mdl(M),#rp(N);");
    @{$clj_AR}        = (1,     1,        0,     1,       1,       0,         0,        0,          0,       0,        0,         0,         0,         0,         0,         1);
  }  
  else {
    ofile_FAIL("ERROR in $sub_name, unrecognized ofile_key $ofile_key", 1, $FH_HR);
  }
  return;
}

#################################################################
# Subroutine: pick_features_from_all_alternatives()
# Incept:     EPN, Wed Oct 13 12:00:26 2021
# Purpose:    For features that have a non-empty 'alternative_ftr_set'
#             value, choose one representative from all alternatives
#             and delete annotation and alerts from all other alternatives.
#             This subroutine should be called twice, first with the
#             $only_children_flag set as '0' to only pick features from 
#             sets that are not children, then again with $only_children_flag 
#             set as '1' to only pick features from sets that *are* children.
#             This is important because when we remove results/alerts for
#             parents we also remove results/alerts for their children.
#             This only works because we validate (in vdr_FeatureInfoValidateAlternativeFeatureSet)
#             - that for any alternative_ftr_set set that includes >= 1 feature
#               child, all members of that set are children, with the same parent.
#             And in vdr_FeatureInfoParentIndexStrings we validate that
#             all features that are parents are not themselves children
#             of any feature.
#
# Arguments:
#  $seq_name_AR:             REF to array of sequence names, PRE-FILLED
#  $ftr_info_AHR:            REF to array of hashes with information on the features, PRE-FILLED
#  $alt_info_HHR:            REF to array of hashes with information on the alerts, PRE-FILLED
#  $ftr_results_HAHR:        REF to feature results HAH, PRE-FILLED
#  $alt_ftr_instances_HHHR:  REF to array of 2D hashes with per-feature alerts, PRE-FILLED
#  $only_children_flag:      '1' to only remove sets that are all children, '0' to remove sets that
#                            are not all children, see 'Purpose'
#  $i_am_child_AR:           REF [0..$nftr-1] to array of 1/0 values on whether each feature is a child, PRE-FILLED
#  $children_AAR:            REF to array of arrays of children feature indices, FILLED HERE, can be undef, PRE-FILLED
#  $opt_HHR:                 REF to 2D hash of option values, see top of sqp_opts.pm for description
#  $FH_HR:                   REF to hash of file handles, including 'log'
#             
# Returns:  void
# 
# Dies:     never
#
#################################################################
sub pick_features_from_all_alternatives { 
  my $sub_name = "pick_features_from_all_alternatives";
  my $nargs_exp = 10;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($seq_name_AR, $ftr_info_AHR, $alt_info_HHR, $ftr_results_HAHR, $alt_ftr_instances_HHHR, 
      $only_children_flag, $i_am_child_AR, $children_AAR, $opt_HHR, $FH_HR) = @_;

  my $nseq = scalar(@{$seq_name_AR});
  my $nftr = scalar(@{$ftr_info_AHR});

  my $ftr_idx; 
  my $ftr_idx2; 
  foreach my $seq_name (@{$seq_name_AR}) { 
    my %sets_completed_H = (); # key is name of a set, value is '1' if we've already completed that set
    if((defined $ftr_results_HAHR->{$seq_name}) || (defined $alt_ftr_instances_HHHR->{$seq_name})) { 
      for($ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
        if(((defined $ftr_results_HAHR->{$seq_name})       && (defined $ftr_results_HAHR->{$seq_name}[$ftr_idx])) || 
           ((defined $alt_ftr_instances_HHHR->{$seq_name}) && (defined $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx}))) { 
          # if(  only_children_flag), we only pick from a set for children
          # if(! only_children_flag), we only pick from a set for non-children
          if(((  $only_children_flag) && (  $i_am_child_AR->[$ftr_idx])) || 
             ((! $only_children_flag) && (! $i_am_child_AR->[$ftr_idx]))) { 
            my $set = $ftr_info_AHR->[$ftr_idx]{"alternative_ftr_set"};
            if(($set ne "") && (! defined $sets_completed_H{$set})) { 
              my @ftr_set_A = ($ftr_idx);
              for($ftr_idx2 = $ftr_idx+1; $ftr_idx2 < $nftr; $ftr_idx2++) { # can start at $ftr_idx+1 b/c we've already covered earlier ftr_idx values  
                if($ftr_info_AHR->[$ftr_idx2]{"alternative_ftr_set"} eq $set) { 
                  push(@ftr_set_A, $ftr_idx2);
                }
              }
              # find 'winner' out of all possible features in @ftr_set_A
              my $nset = scalar(@ftr_set_A);
              if($nset == 1) { 
                ofile_FAIL("ERROR in $sub_name, alt feature set with value $set only has 1 member post-validation", 1, $FH_HR);
              }
              my $nfatal         = undef;
              my $winner_nfatal  = undef;
              my $winner_set_idx = undef;
              my $ftr_set_idx    = undef;
              for($ftr_set_idx = 0; $ftr_set_idx < $nset; $ftr_set_idx++) { 
                $ftr_idx2 = $ftr_set_A[$ftr_set_idx];
                # count fatal alerts for substitute if we have one
                # printf("\tftr_set_idx: $ftr_set_idx, ftr_idx2: $ftr_idx2\n");
                if((defined $ftr_info_AHR->[$ftr_idx2]{"alternative_ftr_set_subn"}) && 
                   $ftr_info_AHR->[$ftr_idx2]{"alternative_ftr_set_subn"} ne "") {
                  $ftr_idx2 = $ftr_info_AHR->[$ftr_idx2]{"alternative_ftr_set_subn"};
                  # printf("\tsubstituted $ftr_idx2 due to alternative_ftr_set_subn values\n");
                }
                my $nfatal = alert_feature_instances_count_fatal($seq_name, $ftr_idx2, $alt_info_HHR, $alt_ftr_instances_HHHR, $FH_HR);
                if($nfatal == 0) { 
                  # no fatal alerts, this is our winner, break the loop          
                  $winner_nfatal  = $nfatal;
                  $winner_set_idx = $ftr_set_idx; 
                  $ftr_set_idx    = $nset; # breaks loop
                }
                elsif((! defined $winner_nfatal) || ($nfatal < $winner_nfatal)) { 
                  # >= 1 fatal alerts, but minimum yet seen, this is our current winner but keep looking
                  $winner_nfatal  = $nfatal;
                  $winner_set_idx = $ftr_set_idx; 
                }
              }
              # go through and remove results and alerts for all non-winners in this set and their children
              # printf("winner_set_idx is $winner_set_idx, ftr_idx: $ftr_set_A[$winner_set_idx]\n");
              for($ftr_set_idx = 0; $ftr_set_idx < $nset; $ftr_set_idx++) { 
                if($ftr_set_idx != $winner_set_idx) { 
                  $ftr_idx2 = $ftr_set_A[$ftr_set_idx];
                  %{$ftr_results_HAHR->{$seq_name}[$ftr_idx2]} = ();
                  %{$alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx2}} = ();
                  undef $ftr_results_HAHR->{$seq_name}[$ftr_idx2];
                  undef $alt_ftr_instances_HHHR->{$seq_name}{$ftr_idx2};
                  my $nchildren = scalar(@{$children_AAR->[$ftr_idx2]}); 
                  # nchildren will always be '0' if $only_children_flag is '1' because 
                  # children can't have children, enforced in vdr_FeatureInfoValidateParentIndexStrings()
                  for(my $child_idx = 0; $child_idx < $nchildren; $child_idx++) { 
                    my $child_ftr_idx = $children_AAR->[$ftr_idx2][$child_idx];
                    %{$ftr_results_HAHR->{$seq_name}[$child_ftr_idx]} = ();
                    %{$alt_ftr_instances_HHHR->{$seq_name}{$child_ftr_idx}} = ();
                    undef $ftr_results_HAHR->{$seq_name}[$child_ftr_idx];
                    undef $alt_ftr_instances_HHHR->{$seq_name}{$child_ftr_idx};
                  }
                }
              }
              $sets_completed_H{$set} = 1;
            }
          }
        }
      }
    }
  }
  return;
}

#################################################################
# Subroutine: count_terminal_ambiguities_in_sqstring
# Incept:     EPN, Tue Nov 23 11:21:38 2021
# Purpose:    Count the number of terminal non-ACGTUs at the beginning 
#             of a sqstring and return it.
#
# Arguments:
#  $sqstring:       sqstring to count ambiguities in
#             
# Returns:  number of ambiguities at beginning of $sqstring
#
#################################################################
sub count_terminal_ambiguities_in_sqstring {
  my $sub_name = "count_terminal_ambiguities_in_sqstring";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqstring) = (@_);

  #printf("in $sub_name, sqstring: $sqstring\n");

  my $pos_retval = undef;
  $sqstring =~ m/[ACGTUacgtu]/g; 
  $pos_retval = pos($sqstring); # returns position of first non-ambiguous nt

  # if $pos_retval is undef entire sqstring is ambiguities
  my $ret_val = undef;
  if(defined $pos_retval) { 
    $ret_val = ($pos_retval - 1);
  }
  else { 
    $ret_val = length($sqstring); 
  }
  return $ret_val;
}

#################################################################
# Subroutine: check_for_ambiguous_nts_in_sqstring
# Incept:     EPN, Tue Nov 23 13:00:16 2021
# Purpose:    Return 1 if there are any characters except canonical
#             upper or lowercase ACGTU in the input string, else
#             return 0;
#
# Arguments:
#  $sqstring:       sqstring to check
#             
# Returns:  '1' if any nts are ambiguous, '0' if not
#
#################################################################
sub check_for_ambiguous_nts_in_sqstring {
  my $sub_name = "check_for_ambiguous_nts_in_sqstring";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($sqstring) = (@_);

  return ($sqstring =~ m/[^ACGTUacgtu]/) ? 1 : 0
}

#################################################################
# Subroutine: helper_feature_terminal_ambiguities
# Incept:     EPN, Wed Dec  8 14:14:35 2021
#
# Purpose:    Count the number of terminal ambiguities at beginning of a
#             feature sqstring and report alerts to %{$alt_str_HR} if
#             defined and if necessary.
#
#             Returns $ablen: the number of consecutive ambiguities at
#             beginning of a sqstring (with caveat about start/stop
#             codon ambiguities explained below ***).
#
#             This subroutine may be called with a reversed feature
#             sqstring. If so $is_reversed will be 1 and alerts
#             relevant to the 3' end will be reported.
#
#             ***If (! $is_trunc): if an ambiguous nt exists in the
#             first 3 positions (start or stop codon) but not all of
#             the first 3 positions are ambiguities consider all 3 of
#             those positions as ambiguities for the purposes of
#             calculating $ablen. This is done to deal with a quirk of
#             NCBI GenBank annotation. A protein translation cannot
#             start with a start or stop codon that has an ambiguity
#             in a non-truncated CDS (i.e. start/stop coordinate
#             prefixed with </>). By considering the first/final 3nt
#             ambiguities, we trim and make it truncated, avoiding
#             this GenBank specific issue.
#
# Arguments:
#  $ftr_sqstring:    sqstring to count ambiguities in
#  $is_reversed:     '1' if sqstring was reversed (so we are dealing with 3' end), '0' if not (5' end)                
#  $is_trunc:        '1' if feature is truncated 
#  $ftr_start:       start position of feature, sequence coords
#  $ftr_stop:        stop position of feature, sequence coords
#  $ftr_strand:      strand of feature
#  $ftr_scoords:     full coords strings (potentially multiple segments of feature)
#  $ftr_len:         total length of feature in nucleotides
#  $ftr_is_cds:      '1' if feature is a CDS, '0' if not
#  $ftr_matches_cds: '1' if feature has same start/stop as a CDS (or is a CDS), '0' if not
#  $tt:              translation table
#  $atg_only:        '1' if only acceptable start is ATG
#  $alt_str_HR:      ref to alert string to add to, undef means don't add to it
#  $ua2rf_AR:        ref to array mapping sequence positions to model positions
#  $FH_HR:           ref to hash of file handles
#
# Returns:  number of ambiguities at beginning of $sqstring, see ***caveat above.
#
#################################################################
sub helper_feature_terminal_ambiguities {
  my $sub_name = "helper_feature_terminal_ambiguities";
  my $nargs_exp = 15;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($ftr_sqstring, $is_reversed, $is_trunc, $ftr_start, $ftr_stop, $ftr_strand, $ftr_scoords, $ftr_len, 
      $ftr_is_cds, $ftr_matches_cds, $tt, $atg_only, $alt_str_HR, $ua2rf_AR, $FH_HR) = (@_);

  my $ret_ablen = count_terminal_ambiguities_in_sqstring($ftr_sqstring);

  # report ambgnt5c/ambgnt5f/ambgnt3c/ambgnt3f
  if(($ret_ablen != 0) && (defined $alt_str_HR)) { 
    if($is_reversed) { 
      # dealing with 3' end
      my $ambg_alt     = ($ftr_is_cds) ? "ambgnt3c" : "ambgnt3f";
      my $ftr_first_ab = vdr_CoordsRelativeSingleCoordToAbsolute($ftr_scoords, ($ftr_len - $ret_ablen + 1), $FH_HR);
      my $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($ftr_first_ab, $ftr_stop, $ftr_strand, $FH_HR) . ";";
      my $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$ftr_first_ab]), abs($ua2rf_AR->[$ftr_stop]), $ftr_strand, $FH_HR) . ";";
      $alt_str_HR->{$ambg_alt} = sprintf("%s%sVADRNULL", $alt_scoords, $alt_mcoords);
    }
    else { 
      # dealing with 5' end 
      my $ambg_alt     = ($ftr_is_cds) ? "ambgnt5c" : "ambgnt5f";
      my $ftr_final_ab = vdr_CoordsRelativeSingleCoordToAbsolute($ftr_scoords, $ret_ablen, $FH_HR);
      my $alt_scoords  = "seq:" . vdr_CoordsSegmentCreate($ftr_start, $ftr_final_ab, $ftr_strand, $FH_HR) . ";";
      my $alt_mcoords  = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$ftr_start]), abs($ua2rf_AR->[$ftr_final_ab]), $ftr_strand, $FH_HR) . ";";
      $alt_str_HR->{$ambg_alt} = sprintf("%s%sVADRNULL", $alt_scoords, $alt_mcoords);
    }
  }

  # Caveat: if we have a CDS or gene with identical start/stop to a CDS that
  #         - is not truncated
  #         - has ambiguities in first 3nt (this the stop codon if $is_reversed and start codon)
  #         - is not a valid start/stop despite ambiguities (e.g. TRA is valid stop in trans table 1)
  #         then we count the start/stop codon as 3 ambiguities and find the first non-N after the first 3
  if((($ftr_is_cds || $ftr_matches_cds)) && (! $is_trunc) && ($ret_ablen == 0) && ($ret_ablen != $ftr_len)) { 
    my $start_or_stop_codon = "";
    if($is_reversed) { 
      $start_or_stop_codon = reverse(substr($ftr_sqstring, 0, 3));
    }
    else { 
      $start_or_stop_codon = substr($ftr_sqstring, 0, 3);
    }
    # check if it has any ambiguities
    if(check_for_ambiguous_nts_in_sqstring($start_or_stop_codon)) { 
      my $codon_is_valid = 0;
      # check if it is a valid start/stop
      if($ftr_len >= 3) { 
        $start_or_stop_codon =~ tr/a-z/A-Z/; # convert to uppercase
        $start_or_stop_codon =~ tr/U/T/;     # convert to DNA
        if($is_reversed) { 
          $codon_is_valid = seq_CodonValidateStopCapDna($start_or_stop_codon, $tt);
        }
        else { 
          $codon_is_valid = seq_CodonValidateStartCapDna($start_or_stop_codon, $tt, $atg_only);
        }
      }
      if(! $codon_is_valid) { 
        my $codon_len = utl_Min(3, $ftr_len);
        $ret_ablen = $codon_len + count_terminal_ambiguities_in_sqstring(substr($ftr_sqstring, 3));

        if((defined $alt_str_HR) && ($ftr_is_cds)) { 
          if($is_reversed) { 
            # dealing with 3' end
            my $ftr_codon_start = vdr_CoordsRelativeSingleCoordToAbsolute($ftr_scoords, ($ftr_len - $codon_len + 1), $FH_HR);
            my $alt_scoords     = "seq:" . vdr_CoordsSegmentCreate($ftr_codon_start, $ftr_stop, $ftr_strand, $FH_HR) . ";";
            my $alt_mcoords     = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$ftr_codon_start]), abs($ua2rf_AR->[$ftr_stop]), $ftr_strand, $FH_HR) . ";";
            $alt_str_HR->{"ambgcd3c"} = sprintf("%s%sVADRNULL", $alt_scoords, $alt_mcoords);
          }
          else { 
            # dealing with 5' end
            my $ftr_codon_end = vdr_CoordsRelativeSingleCoordToAbsolute($ftr_scoords, $codon_len, $FH_HR);
            my $alt_scoords   = "seq:" . vdr_CoordsSegmentCreate($ftr_start, $ftr_codon_end, $ftr_strand, $FH_HR) . ";";
            my $alt_mcoords   = "mdl:" . vdr_CoordsSegmentCreate(abs($ua2rf_AR->[$ftr_start]), abs($ua2rf_AR->[$ftr_codon_end]), $ftr_strand, $FH_HR) . ";";
            $alt_str_HR->{"ambgcd5c"} = sprintf("%s%sVADRNULL", $alt_scoords, $alt_mcoords);
          }
        }        
      }
    }
  }

  return $ret_ablen;
}
