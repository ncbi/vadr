#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Mon Feb  1 15:07:43 2016 [dnaorg_build.pl split off from dnaorg_annotate_genomes.pl]
# LES, Mon Jul 25 09:25    2016 [dnaorg_refseq_assign.pl split off from dnaorg_build.pl]
#
# Some code in here is also adopted from dnaorg.pm

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

use Data::Dumper;
use File::Slurp;

require "dnaorg.pm"; 
require "epn-options.pm";

#######################################################################################
# What this script does: 
#
# Preliminaries: 
#   - process options
#   - create the output directory
#   - output program banner and open output files
#   - parse the input files
#
# Step 1. Create an HMM library of RefSeqs
#
# Step 2. Run nhmmscan for all accessions
#
# Step 3. Parse the nhmmscan output to determine the proper RefSeq for each accession
#
# Step 4. Generate ntlists for each RefSeq
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
my $hmmer_exec_dir    = $dnaorgdir . "/hmmer-3.1b2/src/";
my $esl_exec_dir      = $dnaorgdir . "/infernal-1.1.2/easel/miniapps/";


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
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",           "force; if dir <reference accession> exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                            "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);
opt_Add("--dirout",     "string",  undef,                    1,    undef, undef,      "output directory specified as <s>",     "specify output directory as <s>, not <ref accession>", \%opt_HH, \@opt_order_A);
opt_Add("--keep",       "boolean", 0,                        1,    undef, undef,      "leaving intermediate files on disk",    "do not remove intermediate files, keep them all on disk", \%opt_HH, \@opt_order_A);


# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_refseq_assign.pl [-options] <RefSeq list> <sequence list>\n";
my $synopsis = "dnaorg_refseq_assign.pl :: Given sequences, decides which RefSeq's ntlist to add them to";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'dirout=s'     => \$GetOptions_H{"--dirout"},
                'keep'         => \$GetOptions_H{"--keep"}
                );

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.20";
my $releasedate   = "Nov 2017";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 2) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do dnaorg_build.pl -h\n\n";
  exit(1);
}
my ($ref_list, $seq_list) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);


my $dir = opt_Get("--dirout", \%opt_HH);          # this will be undefined unless -d set on cmdline

#############################
# create the output directory
#############################
my $cmd;              # a command to run with runCommand()
my @early_cmd_A = (); # array of commands we run before our log file is opened
# check if the $dir exists, and that it contains the files we need
# check if our output dir $symbol exists
if(! defined $dir) { 
  $dir = $seq_list . "-ntlists";
}
else { 
  if($dir !~ m/\/$/) { $dir =~ s/\/$//; } # remove final '/' if it exists
}
if(-d $dir) { 
  $cmd = "rm -rf $dir";
  if(opt_Get("-f", \%opt_HH)) { runCommand($cmd, opt_Get("-v", \%opt_HH), undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR directory named $dir already exists. Remove it, or use -f to overwrite it."; }
}
if(-e $dir) { 
  $cmd = "rm $dir";
  if(opt_Get("-f", \%opt_HH)) { runCommand($cmd, opt_Get("-v", \%opt_HH), undef); push(@early_cmd_A, $cmd); }
  else                        { die "ERROR a file named $dir already exists. Remove it, or use -f to overwrite it."; }
}

# create the dir
$cmd = "mkdir $dir";
runCommand($cmd, opt_Get("-v", \%opt_HH), undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $out_root = $dir . "/" . $dir_tail;

#############################################
# output program banner and open output files
#############################################
# output preamble
my @arg_desc_A = ("list of refseqs", "list of sequences");
my @arg_A      = ($ref_list, $seq_list);
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
                         # 2D keys:
                         #  "log": log file of what's output to stdout
                         #  "cmd": command file with list of all commands executed
                         #  "list": list and description of output files
                         #  Per accession '$seq':
                         #                       $seq.tbl
                         #                       $seq.results

# open the log and command files 
openAndAddFileToOutputInfo(\%ofile_info_HH, "log", $out_root . ".log", 1, "Output printed to screen");
openAndAddFileToOutputInfo(\%ofile_info_HH, "cmd", $out_root . ".cmd", 1, "List of executed commands");
openAndAddFileToOutputInfo(\%ofile_info_HH, "list", $out_root . ".list", 1, "List and description of all output files");

my $log_FH = $ofile_info_HH{"FH"}{"log"};
my $cmd_FH = $ofile_info_HH{"FH"}{"cmd"};
# output files are all open, if we exit after this point, we'll need
# to close these first.


# now we have the log file open, output the banner there too
outputBanner($log_FH, $version, $releasedate, $synopsis, $date);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $cmd_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}

my $do_keep = opt_Get("--keep", \%opt_HH); # should we leave intermediates files on disk, instead of removing them?
my @files2rm_A = ();                       # will be filled with files to remove, --keep was not enabled

###################################################
# make sure the required executables are executable
###################################################
my %execs_H = (); # hash with paths to all required executables
$execs_H{"nhmmscan"}      = $hmmer_exec_dir . "nhmmscan";
$execs_H{"hmmbuild"}      = $hmmer_exec_dir . "hmmbuild";
$execs_H{"hmmpress"}      = $hmmer_exec_dir . "hmmpress";
$execs_H{"esl-reformat"}  = $esl_exec_dir   . "esl-reformat";
validateExecutableHash(\%execs_H, $ofile_info_HH{"FH"});

#################################################################################
#################################################################################
#
# LES Jul 26 2016
#
# ASSUMPTIONS: $ref_list    is the file name of a list of RefSeq accns
#              $seq_list    is a file which contains the accession numbers of all the
#                           sequences which are to be assigned to ntlists (one accn #
#                           per line)
#
#################################################################################

my $progress_w = 80; # the width of the left hand column in our progress output, hard-coded                                                     
my $start_secs = outputProgressPrior("Creating RefSeq HMM Library", $progress_w, $log_FH, *STDOUT);

# create fasta files for refseqs
createFastas($ref_list, $out_root, \%opt_HH, $ofile_info_HH{"FH"});

# eliminate blank lines in $ref_list
my $blank_match = '^\s*$';
$cmd = "sed -i '/$blank_match/d' $ref_list";
runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});

# open and parse list of RefSeqs
open(REFSEQ_LIST, $ref_list) || fileOpenFailure($ref_list, $0, $!, "reading", $ofile_info_HH{"FH"});
my @refseqs_A = (); #array of refseqs
while (<REFSEQ_LIST>) {
    my $ref_tmp = $_;
    chomp $ref_tmp;
    $ref_tmp =~ s/\s*$//;
    push(@refseqs_A, $ref_tmp);
}

## DEBUG
#print Dumper(\@refseqs_A);

# create HMM library
my $ref_library = $out_root . ".hmm";

# make a .stk and .hmm file for each RefSeq
foreach (@refseqs_A) {
    my $fa_file = $out_root . ".$_.fasta";
    my $stk_file = $out_root . ".$_.stk";
    $cmd = $execs_H{"esl-reformat"} . " --informat afa stockholm $fa_file > $stk_file";
    runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});

    my $hmm_file = $out_root . ".$_.hmm";
    $cmd = $execs_H{"hmmbuild"} . " $hmm_file $stk_file > /dev/null";  # dumps the output of hmmbuild
    runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
    
    $cmd = "cat $hmm_file >> $ref_library";   # adds each individual hmm to the hmm library
    runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});

    # clean up intermediate files
    if(! $do_keep) {
	push(@files2rm_A, $fa_file);
	push(@files2rm_A, $stk_file);
	push(@files2rm_A, $hmm_file);
    }
}

addClosedFileToOutputInfo(\%ofile_info_HH, "HMMLib", $ref_library, 1, "Library of HMMs of RefSeqs. Not press'd");

$cmd = $execs_H{"hmmpress"} . " $ref_library > /dev/null";
runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

########### RUN nhmmscan and generate output files ####################################################################################################

$start_secs = outputProgressPrior("Running nhmmscan on sequences", $progress_w, $log_FH, *STDOUT);

# create fasta files for all accessions
createFastas($seq_list, $out_root, \%opt_HH, $ofile_info_HH{"FH"});

# open, get rid of empty lines, and parse list of sequences
$cmd = "sed -i '/$blank_match/d' $seq_list";
runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
open(SEQ_LIS, $seq_list) || fileOpenFailure($seq_list, $0, $!, "reading", $ofile_info_HH{"FH"});     
my @seqs_to_assign_A = (); #array that will contain all the accessions in $seq_list (all the accns to be assigned to a ntlist)
while (<SEQ_LIS>) {
    chomp;
    push(@seqs_to_assign_A, $_);
}

# creates a tblout summary of the nhmmscan results and a file with nhmmscan results for each accession
nhmmscanSeqs($execs_H{"nhmmscan"}, $out_root, \@seqs_to_assign_A, $ref_library, \@files2rm_A, \%ofile_info_HH); 


# clean up intermediate files (hmmpress binaries)
if(! $do_keep) {
    push(@files2rm_A, "$ref_library.h3m");
    push(@files2rm_A, "$ref_library.h3f");
    push(@files2rm_A, "$ref_library.h3i");
    push(@files2rm_A, "$ref_library.h3p");
}


outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);


########################################################################################################################################################
#
# PARSE nhmmscan files
# loop through sequences
# loop through lines
#
# Chooses a RefSeq for each sequence
# Creates: .matches.info file giving information on each match
# Creates: ntlists for each RefSeq
#
########################################################################################################################################################

$start_secs = outputProgressPrior("Creating match info file", $progress_w, $log_FH, *STDOUT);

# Generate match information file header                                                                               
my $match_file = $out_root . ".matches.info";
open(MATCH_INFO, "> $match_file") || fileOpenFailure($match_file, $0, $!, "writing", $ofile_info_HH{"FH"});
print MATCH_INFO "########################################################################################################################################\n";
print MATCH_INFO "#\n";
print MATCH_INFO "# Query:       Accession number of the sequence \n";
print MATCH_INFO "# RefSeq:      The RefSeq that the sequence was assigned to\n";
print MATCH_INFO "# Bit score:   Bit score of hit to 'RefSeq'\n";
print MATCH_INFO "# E-val:       E-val of hit to 'RefSeq'\n";
print MATCH_INFO "# Coverage:    The percentage of the query that the hit to 'RefSeq' covers (Hit length)/(Query length)\n";
print MATCH_INFO "# Bias:        TODO\n";
print MATCH_INFO "# # Hits:      The number of individual hits to this RefSeq (the program combines stats such as bit score and covg. from separate hits)\n";
print MATCH_INFO "# H2: RefSeq:  The RefSeq that had the second strongest hit\n";
print MATCH_INFO "# Bit Diff:    The amount by which the bit score of the 'RefSeq' hit is greater than that of the 'H2: RefSeq' hit\n";
print MATCH_INFO "# Covg. Diff:  The amount by which the coverage of the 'RefSeq' hit is greater than that of the 'H2: RefSeq' hit\n";
print MATCH_INFO "# Num. Correct Hits: The amount of times 'Exp. RefSeq' produced a hit\n";
print MATCH_INFO "#\n";
print MATCH_INFO "########################################################################################################################################\n";
print MATCH_INFO "\n";
print MATCH_INFO "\n";
print MATCH_INFO "#H ";
my @header = ("Query   ","RefSeq   ","Bit score","E-val","Coverage","Bias","# Hits","H2: RefSeq","Bit Diff","Covg. Diff \n");
print MATCH_INFO join("\t\t", @header);
print MATCH_INFO "# ";
print MATCH_INFO "--------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
print MATCH_INFO "\n";


# Generate data structures to build ntlists from
my %ntlist_HA = (); # Hash of arrays containing each RefSeq's ntlist
                    # Key:    RefSeq accession #
                    # Value:  Array of accessions for which have been assigned to this RefSeq
foreach (@refseqs_A) {
    $ntlist_HA{$_} = (); # initialize each RefSeq's value to an empty array
}

# initialize the non-assigned list to an empty array
$ntlist_HA{"non-assigned"} = ();


# Generate data structures to build a sequence profile for further evaluation                       
my %hit_info_HAA = (); # Hash of hash of arrays containing nhmmscan info for each hit for each sequence
                       # 1D Key:    Sequence accession #                                                                             
                       # 1D Array:  Order of hit (1st nhmmscan hit = index 0, 2nd = index 1, etc)                                                                                                      
                       # Value:     Array w/ relevant nhmmscan info for the hit (see comments below)


foreach my $seq (@seqs_to_assign_A) {
    my $seq_tbl_file = $out_root . ".$seq.tblout";
    my $eligible_hits = 0; # Boolean set to false
#    # debug
#    if(-e $seq_tbl_file) {
#	print "$seq_tbl_file exists \n";
#    } else {
#	print "doesn't exist \n";
#    }
#    if(-s $seq_tbl_file) {
#	print "$seq_tbl_file is non-empty \n";
#    } else {
#	print "is empty \n";
#    }

    open(NHMMSCANTBL, "< $seq_tbl_file") || fileOpenFailure($seq_tbl_file, $0, $!, "reading", $ofile_info_HH{"FH"});

    # initialize $seq in the hash
    if( ! exists($hit_info_HAA{$seq}) ) {
	@{$hit_info_HAA{$seq}} = ();

    }

    my @hit_specs_A;
    my $counter = 0;
    my @hit_output_A;
    my $refseq;

    while(my $line = <NHMMSCANTBL>) {

	if($line !~ m/^\#/) {               #if this line is not commented
	    $eligible_hits = 1; # Boolean set to true - there are eligible hits for this sequence

	    ##########################################################################
	    # In @hit_specs_A - array containing all the output info for the hit on this line
	    #
	    # index         feature
	    #-----------------------
	    # 0             target name
	    # 1             accession (not sure what this refers to, shows up as a - for every hit in my output)
	    # 2             query name (name from > line of fasta file)
	    # 3             accession (again, ?)
	    # 4             hmmfrom
	    # 5             hmm to
	    # 6             alifrom
	    # 7             ali to
	    # 8             env from
	    # 9             env to
	    # 10            modlen (entire length of model seq)
	    # 11            strand (+,-)
	    # 12            E-value
	    # 13            bit score
	    # 14            bias
	    # 15            description of target (showing up as - for every hit in my output)
	    ##########################################################################

	    @hit_specs_A = split(/\s+/, $line);
            $refseq = $hit_specs_A[0];
            $refseq =~ s/^.+\.//;

            # ########################################################################
            # In @hit_output_A - array containing all the output for this hit
            #
            # index         feature
            #--------------------------
            # 0             Query accn
            # 1             RefSeq accn
            # 2             Bit-score
            # 3             E-value
            # 4             Coverage (Hit length)/(Query length)
            # 5             Bias
	    # 6             Number of hits to this RefSeq
            ###########################################################################
	    # check if this RefSeq has appeared in a previous hit for this sequence                                              
            my $first_index = undef;  # the first occuring index of a previous hit, if one exists                                                             
            # if there are previous hits, search through them to see if this RefSeq has been hit before
	    if( @{$hit_info_HAA{$seq}} != 0 ) {
                for(my $i=0;  $i < scalar(@{$hit_info_HAA{$seq}});  $i++) {
                    if( @{@{$hit_info_HAA{$seq}}[$i]}[1] eq $refseq ) {       # if index $i's RefSeq is the same as this one    
                        $first_index = $i;
                    }
                }
            }

	    ## Find length of $seq - will be used to find Coverage                                                                                                                                                                                       
	    my $full_nhmmscan_file = $out_root . "." . $seq . ".nhmmscan.out";
	    my $seq_len = `grep 'Query:' $full_nhmmscan_file`;
	    # $seq_len is now in the format "Query:   gi|nums|db|accn#| [L=length]                                                                                      
	    $seq_len =~ s/Query:\s+\S+\s+//;
	    # $seq_len is now in the format "[L=length]"                                                                                                                          
	    $seq_len =~ s/\[L=//;
	    $seq_len =~ s/\]$//;
		
	    # Deciding hit length from alito alifrom - will be used to find Coverage
	    my $hit_len = $hit_specs_A[7] - $hit_specs_A[6];
		

	    # if this is the first hit to this RefSeq
	    if(! defined($first_index)) {
		#debug
		#print "first occurence of $seq \n";

		@hit_output_A = (); # Array to be outputed to match info file                                                                                                                          
		push(@hit_output_A, $seq);
		push(@hit_output_A, $refseq);
           
		push(@hit_output_A, $hit_specs_A[13]); # add bit score
		my $e_val = sprintf("%.7f", $hit_specs_A[12]);
		push(@hit_output_A, $e_val); # add E-val                                                                                                            

		# add coverage				
		my $coverage = $hit_len/$seq_len;
		$coverage = sprintf("%.7f", $coverage);
		push(@hit_output_A, $coverage);
		
		push(@hit_output_A, $hit_specs_A[14]);  # add bias                                                                                                                            
		push(@hit_output_A, 1);                 # initialize 'Number of hits' to 1

                                                                                                                               
		@{$hit_info_HAA{$seq}}[$counter] = ();
		@{$hit_info_HAA{$seq}[$counter]} = @hit_output_A;
		
		$counter++;

	    }else {    # if this is a hit to a sequence that has already been hit
		#@hit_output_A = \@{$hit_info_HAA{$seq}[$first_index]};
  
		#debug
		#print "second, third, etc hit of $seq \n";

		# TODO - check with Eric!
		@{$hit_info_HAA{$seq}[$first_index]}[2] += $hit_specs_A[13]; # Add bit score
		@{$hit_info_HAA{$seq}[$first_index]}[3] += $hit_specs_A[12]; # Add E-val
		@{$hit_info_HAA{$seq}[$first_index]}[5] += $hit_specs_A[14]; # Add bias
		@{$hit_info_HAA{$seq}[$first_index]}[6] ++;                  # increase 'Number of hits' by 1

		##############################################
		# Calculate and add coverage
		my $prev_covg = @{$hit_info_HAA{$seq}[$first_index]}[4];
		my $prev_hit_len = $prev_covg*$seq_len;
		my $coverage = ($prev_hit_len + $hit_len) / $seq_len;
		$coverage = sprintf("%.7f", $coverage);
		@{$hit_info_HAA{$seq}[$first_index]}[4] = $coverage;
		#############################################
	    }  

	}
    }
    # debug - prob not necessary, but won't hurt
    close NHMMSCANTBL || die "couldn't close NHMMSCANTBL. oops! \n";

    if( $eligible_hits ) { # if there was at least one valid hit for this sequence
	                   # then parse and output the best hit's info
	
	# TODO
	#      use %hit_info_HAA to decide what the correct RefSeq is
	#
	
	
	# sort order of hits by bit score
	@{$hit_info_HAA{$seq}} = sort { $a->[2] <=> $b->[2] } @{$hit_info_HAA{$seq}};
	
	# Assigns the sequence to the refseq hit with the the highest bit score
	#my $max_index = 0; # contains the index of the hit with the highest bit score - initialized at index 0
	#my $second_index = 0;
	#for( my $hit=0; $hit< scalar(@{$hit_info_HAA{$seq}}) ; $hit++  ) {
	#	if( @{$hit_info_HAA{$seq}[$hit]}[2] > @{$hit_info_HAA{$seq}[$max_index]}[2] ) {
	#	    $max_index = $hit;
	#	}
	#}
	
	
	
	# DEBUG
	#print "This is right before the error\n";
	#print Dumper(\%hit_info_HAA);
	
	# assign hit_output_A to the assigned refseq hit's info, assign $refseq to the assigned refseq
	@hit_output_A = @{@{$hit_info_HAA{$seq}}[-1]}; 
	$refseq = $hit_output_A[1];
	
	##############################################################################
	# Add comparison indicies to @hit_output_A
	#
	# index        feature
	# ----------------------
	# 7            RefSeq accn # for hit 2
	# 8            difference between bit scores of hit 1 and hit 2
	# 9            difference between coverage of hit 1 and hit 2
	##############################################################################
	
	# TODO - hardcode-y 
	if( scalar(@{$hit_info_HAA{$seq}}) >= 2) {     # if there is more than one hit  
	    push(@hit_output_A, @{@{$hit_info_HAA{$seq}}[-2]}[1]);
	
	    my $bit_diff = $hit_info_HAA{$seq}[-1][2] - $hit_info_HAA{$seq}[-2][2];
	    $bit_diff = sprintf("%9.1f", $bit_diff);
	    my $cov_diff = $hit_info_HAA{$seq}[-1][4] - $hit_info_HAA{$seq}[-2][4];
	    $cov_diff = sprintf("%.7f", $cov_diff);
	    push(@hit_output_A, $bit_diff);
	    push(@hit_output_A, $cov_diff);
	
	} else { # if there was only one hit, there's no comparison
	    push(@hit_output_A, "-----");
	    push(@hit_output_A, "-----");
	    push(@hit_output_A, "-----");
	}

	# add hit info to .matches.info file
	print MATCH_INFO join("\t\t", @hit_output_A);
	print MATCH_INFO "\n";
	
	# Add this seq to the hash key that corresponds to it's RefSeq
	push(@{$ntlist_HA{$refseq}}, $seq);


    } else { # if there were no eligible hits
	@hit_output_A = ($seq,"--------","--------","--------","--------","--------","--------","--------","--------","--------");
	print MATCH_INFO join("\t\t", @hit_output_A);
	print MATCH_INFO "\n";

	# add this sequence to the file that lists non-assigned sequences
	push(@{$ntlist_HA{"non-assigned"}}, $seq);

    }

}

addClosedFileToOutputInfo(\%ofile_info_HH, "match-info", $match_file, 1, "Table with statistics for each match");

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);


$start_secs = outputProgressPrior("Creating ntlists and other output files", $progress_w, $log_FH, *STDOUT);

# Generate ntlist files
foreach my $refseq (@refseqs_A) {
    my $ntlist_file = $out_root . ".$refseq.ntlist";

    open(NTLIST, "> $ntlist_file")  || fileOpenFailure($ntlist_file, $0, $!, "writing", $ofile_info_HH{"FH"});
    print NTLIST "$refseq\n";
    foreach (@{$ntlist_HA{$refseq}}) {
	print NTLIST "$_\n";
    }
    
    addClosedFileToOutputInfo(\%ofile_info_HH, "$refseq.ntlist", $ntlist_file, 1, "ntlist for $refseq");

}

# Generate file that lists non-assigned sequences
my $non_assigned_file = $out_root . ".non-assigned";

open(NALIST, "> $non_assigned_file")  || fileOpenFailure($non_assigned_file, $0, $!, "writing", $ofile_info_HH{"FH"});
foreach (@{$ntlist_HA{"non-assigned"}}) {
    print NALIST "$_\n";
}

addClosedFileToOutputInfo(\%ofile_info_HH, "non-assigned", $non_assigned_file, 1, "List of sequences not assigned to a RefSeq");


# add $ref_list and #seq_list to output direcotry
my $out_ref_list = $out_root . ".all.refseqs";
my $out_seq_list = $out_root . ".all.seqs";
$cmd = "cat $ref_list > $out_ref_list";
runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});
$cmd = "cat $seq_list > $out_seq_list";
runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});


# eliminate blank lines in $out_seq_list
$cmd = "sed -i '/$blank_match/d' $out_seq_list";
runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});

# eliminate blank lines in $out_ref_list
$cmd = "sed -i '/$blank_match/d' $out_ref_list";
runCommand($cmd, opt_Get("-v", \%opt_HH), $ofile_info_HH{"FH"});


addClosedFileToOutputInfo(\%ofile_info_HH, "RefSeqs", $out_ref_list, 1, "List of RefSeqs in the HMM library");
addClosedFileToOutputInfo(\%ofile_info_HH, "Seqs", $out_seq_list, 1, "List of sequences that were sorted into ntlists");

outputProgressComplete($start_secs, undef, $log_FH, *STDOUT);

##########################################################################################################################################################

##########
# Conclude
##########

# if not --keep, then remove unnessecary files (while checking to make sure that you're not removing the same file twice
my %seen;

foreach my $file2rm (@files2rm_A) {
    next if $seen{$file2rm}++;

    runCommand("rm $file2rm", 0, $ofile_info_HH{"FH"});
}

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir, \%ofile_info_HH);
exit 0;



#################################################################################
#
# Sub name:  createFastas()
#
# Author:    Lara Shonkwiler
# Date:      2016 Aug 01
#
# Purpose:   Given a list of accessions, creates a fasta file for each one
# 
# Arguments: $accn_list    list of all accessions for fastas to be made for
#            $out_root     directory path
#            $opt_HHR      reference to hash of options
#            $FH_HR        reference to hash of file handles
#
#
# Returns:   void
#
#################################################################################
sub createFastas {
    my $sub_name  = "createFastas()";
    my $nargs_expected = 4;
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

    my ($accn_list, $out_root, $opt_HHR, $FH_HR) = (@_);

    # create a file containing the concatenated fasta files for all refseqs
    my $all_refseq_fastas = $out_root . ".all.fa";
    $cmd = "cat $accn_list | epost -db nuccore -format acc | efetch -format fasta > $all_refseq_fastas";
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);

    # separate that file into individual fasta files
    # This section contains code adopted from break_fasta.pl (created by Alejandro Schaffer)
    # This section also contains code adopted from entrez_names.pl (created by Lara Shonkwiler)
    my $infile; #input FASTA file with many sequences                                                                                        
    my $outfile; #output FASTA file with one sequence                                                                                                                                                                 
    my $nofile = 1;
    my $nextline;
    my $same_sequence = 2;
    my $state = $nofile;
    my $id;
    my $new_file_name;

    open(SEQS, $all_refseq_fastas) or die "Cannot open $all_refseq_fastas\n";
    while(defined($nextline = <SEQS>)) {
	chomp($nextline);
	if ($nofile == $state) {
	    #$id = $nextline;
	    ($id) = ($nextline =~ m/^>(\S+)/);               # get rid of > and words after space                                                                                                                              
	    $id =~ s/^gi\|?\d+\|\w+\|//; # get rid of everything before the accession number                                                                                                                  
            # special case to deal with accessions that being with pdb\|, e.g. 'pdb|5TSN|T' which is the 
            # name of the sequence that gets fetched for the accession 5TSN|T
            if($id =~ /^pdb\|(\S+)\|(\S+)$/) { 
              $id = $1 . "_" . $2;
            }
	    else { 
              $id =~ s/\|$//;               # get rid of end | or                                                                                                                                  
              $id =~ s/\|\w+//;             # get rid of end | and everything after |
            }        
	    # gets rid of version number
            $id =~ s/\.\d+$//;

	    $new_file_name = $out_root . "." . $id . "." . "fasta";
	    open(ONEFILE, ">$new_file_name") or die "Cannot open $new_file_name\n";
	    $state = $same_sequence;
	    print ONEFILE "$nextline";
	    print ONEFILE "\n";
	}
	else {
	    if ($same_sequence == $state) {
		if ($nextline =~ m/>/ ) {
		    close(ONEFILE);
		    #$id = $nextline;
		    ($id) = ($nextline =~ m/^>(\S+)/);               # get rid of > and words after space                                                                                                                              
		    $id =~ s/^gi\|?\d+\|\w+\|//; # get rid of everything before the accession number                                                                                                     
                    # special case to deal with accessions that being with pdb\|, e.g. 'pdb|5TSN|T' which is the 
                    # name of the sequence that gets fetched for the accession 5TSN|T
                    if($id =~ /^pdb\|(\S+)\|(\S+)$/) { 
                      $id = $1 . "_" . $2;
                    }
                    else { 
                      $id =~ s/\|$//;               # get rid of end | or                                                                                                                                  
                      $id =~ s/\|\w+//;             # get rid of end | and everything after |
                    }        
		    
		    # gets rid of version number                                                                                                                                                                                                         
		    if($id =~ m/\.\d+$/) {
		     $id =~ s/\.\d+$//;
		    }
 
		    $new_file_name = $out_root . "." . $id . "." . "fasta";
		    open(ONEFILE, ">$new_file_name") or die "Cannot open $new_file_name\n";
		    $state = $same_sequence;
		}
		print ONEFILE "$nextline";
		print ONEFILE "\n";
	    }
	}
    }
    close(ONEFILE);
    close(SEQS);

    # get rid of concatenated fasta file
    $cmd = "rm $all_refseq_fastas";
    runCommand($cmd, opt_Get("-v", $opt_HHR), $FH_HR);
}


#######################################################################################
#
# Sub name:  nhmmscanSeqs()
#
# Author:    Lara Shonkwiler (adopted from code by Alejandro Schaffer)
# Date:      2016 Aug 01
#
# Purpose:   Calls nhmmscan for each seq against the given HMM library
# 
# Arguments: $nhmmscan            - path to hmmer-3.1b2 nhmmscan executable
#            $out_root            - path to the output directory
#            $seqs_to_assign_AR   - unique identifier that distinguishes this command from others
#            $ref_library         - HMM library of RefSeqs
#            $files2rm_AR         - REF to a list of files to remove at the end of the
#                                   program
#            $ofile_info_HHR      - REF to the output file info hash
#            
#            
#
# Returns:   void
#
#######################################################################################
sub nhmmscanSeqs {
    my $sub_name  = "nhmmscanSeqs()";
    my $nargs_expected = 6;
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

    my ($nhmmscan, $out_root, $seqs_to_assign_AR, $ref_library, $files2rm_AR, $ofile_info_HHR) = (@_);

    my $qsub_script_name = $out_root . ".qsub.script";                                                                                                                                         
    open(QSUB, ">$qsub_script_name") || fileOpenFailure($qsub_script_name, $sub_name, $!, "writing", $ofile_info_HHR->{"FH"});
    foreach my $seq (@{$seqs_to_assign_AR}) {
	my $seq_tbl_file = $out_root . ".$seq.tblout";
	my $seq_results_file = $out_root . ".$seq.nhmmscan.out";
	my $fa_file = $out_root . ".$seq.fasta";

	# --noali , --cpu 0
	my $cmd = "$nhmmscan --noali --cpu 0 --tblout $seq_tbl_file $ref_library $fa_file > $seq_results_file";

      	my $output_script_name = "$out_root" . "." . "$seq" . "\.qsub\.csh";                                                                                                                               
	open(SCRIPT, ">$output_script_name") or die "Cannot open 4 $output_script_name\n";
	print QSUB "chmod +x $output_script_name\n";
	print QSUB "qsub $output_script_name\n";
	print SCRIPT "source /etc/profile\n";
	print SCRIPT "\#!/bin/tcsh\n";
	print SCRIPT "\#\$ -P unified\n";
	print SCRIPT "\n";
	print SCRIPT "\# list resource request options\n";
	print SCRIPT "\#\$  -l h_rt=288000,h_vmem=32G,mem_free=32G,reserve_mem=32G  \n";
	# old version of above line: -l h_vmem=32G,reserve_mem=32G,mem_free=32G
	# new version of above line: -l h_rt=288000,h_vmem=32G,mem_free=32G,reserve_mem=32G
	print SCRIPT "\n";
	print SCRIPT "\# split stdout and stderr files (default is they are joined into one file)\n";
	print SCRIPT "\#\$ -j n\n";
	print SCRIPT "\n";
	print SCRIPT "\#define stderr file\n";
	my $error_file_name = "$out_root." . "$seq" . "\.qsub\.err";                                                                                                                               
	print SCRIPT "\#\$ -e $error_file_name\n";
	print SCRIPT "\# define stdout file\n";
	my $diagnostic_file_name = "$out_root" . "." . "$seq" . "\.qsub\.out";                                                                                                                               
	print SCRIPT "\#\$ -o $diagnostic_file_name\n";
	print SCRIPT "\n";
	print SCRIPT "\# job is re-runnable if SGE fails while it's running (e.g. the host reboots)\n";
	print SCRIPT "\#\$ -r y\n";
	print SCRIPT "\# stop email from being sent at the end of the job\n";
	print SCRIPT "\#\$ -m n\n";
	print SCRIPT "\n";
	print SCRIPT "\# trigger NCBI facilities so runtime enviroment is similar to login environment\n";
	print SCRIPT "\#\$ -v SGE_FACILITIES\n";
	print SCRIPT "\n";
	print SCRIPT "echo \"Running qsub\"\n";                                                                                                                               
	print SCRIPT "\n";
	print SCRIPT "$cmd";                                                                                                                               
	print SCRIPT "\n";
	close(SCRIPT);                                                                                                                               
	system("chmod +x $output_script_name");                                                                                                          

	addClosedFileToOutputInfo($ofile_info_HHR, "$seq.tbl", $seq_tbl_file, 1, "Table summarizing nhmmscan results for $seq against $ref_library");
	#addClosedFileToOutputInfo(\%ofile_info_HH, "$seq.results", $seq_results_file, 1, "nhmmscan results for searching $ref_library with $seq");                                         

	# clean up leftover fasta files and qsub files
	if(! $do_keep) {
	    push(@{$files2rm_AR}, $fa_file);
	    push(@{$files2rm_AR}, $seq_results_file);
	    push(@{$files2rm_AR}, $output_script_name);
	    push(@{$files2rm_AR}, $error_file_name);
	    push(@{$files2rm_AR}, $diagnostic_file_name);
	}
    
    }

    close(QSUB);                                                                                                                               
    system("chmod +x $qsub_script_name");                                                                                                                               
    system("$qsub_script_name");

    # wait until all nhmmscan jobs are done before proceeding

#    # only checks last seq - causes null array ref error
#    my $last_seq = @{$seqs_to_assign_AR}[-1];
#    my $last_seq_file = $out_root . "." . $last_seq . ".tblout";
#    my $done = 0; # Boolean - set to false
#    my $last_line;
#
#   until( $done ) {
#	sleep(10);
#
#	if(-s $last_seq_file){
#	    $last_line = `tail -n1 $last_seq_file`;
#	    if( $last_line =~ m/\[ok\]/ ) {
#		$done = 1;
#	    }
#	}
#    }


#   slow, but accurate version
#    my $last_line;
#    
#  for my $seq (@{$seqs_to_assign_AR}) {
#	my $done = 0; # Boolean - set to false
#	
#	until( $done ) {
#	    sleep(10);
#	    #debug
#	    print "waiting on $seq";
#
#	    my $seq_file = $out_root . "." . $seq . ".tblout";
#	    if(-s $seq_file){
#		$last_line = `tail -n1 $seq_file`;
#		if( $last_line =~ m/\[ok\]/ ) {
#		    $done = 1;
#		}
#	    }
#	}
#    }
    
    

    # This one seems to work well
    my $last_line;
    
    my $counter = 0;
    my $end = scalar(@{$seqs_to_assign_AR});
 
    until($counter == $end){
	$counter = 0;

	for my $seq (@{$seqs_to_assign_AR}) {
	    my $seq_file = $out_root . "." . $seq . ".tblout";
	    my $full_seq_file = $out_root . "." . $seq . ".nhmmscan.out";
	    if( (-s $seq_file) && (-s $full_seq_file)) {
		$last_line = `tail -n1 $seq_file`;
		if( $last_line =~ m/\[ok\]/ ) {
		    $counter++;
		}
	    }
	}

	sleep(10);
    }


    
}
