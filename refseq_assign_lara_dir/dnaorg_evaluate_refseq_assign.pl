########################################################################################################
#
# Program:               dnaorg_evaluate_refseq_assign.pl
# Author:                Lara Shonkwiler
# Date:                  Aug 04, 2016
#
# Description:           Outputs a file which shows the nhmmscans for the correct hit, and the difference 
#                        between that hit and the next best hit
#
# Inputs:        $results      - Folder which contains the results for a run of dnaorg_refseq_assign.pl                                                             
#                $standards    - Folder containing 'gold standard' ntlists
#                                ntlists must be named as <accn_#>.ntlist
#                                                                                     
# Outputs:       $hit_info     - File which shows detailed nhmmscan and comparison  info for each hit
#
#
#
#
########################################################################################################               

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;

use Data::Dumper;
use File::Slurp;

require "dnaorg.pm";
require "epn-options.pm";


########################################################                                                                                                                                                                                                     
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
my $usage    = "Usage: dnaorg_evaluate_refseq_assign.pl [-options] <dnaorg_refseq_assign results folder> <folder containing 'gold standard' ntlists\n";
my $synopsis = "dnaorg_evaluate_refseq_assign.pl :: Outputs info which helps the user evaluate dnaorg_refseq_assign.pl's performance";

my $options_okay =
    &GetOptions('h'            => \$GetOptions_H{"-h"},
# basic options                                                                                                                                                                                                 
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
                'dirout=s'     => \$GetOptions_H{"--dirout"},
                'keep'         => \$GetOptions_H{"--keep"},
    );

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time                                                 
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.42";
my $releasedate   = "Dec 2018";

# print help and exit if necessary                                                                                                                                                                                                                           
if((! $options_okay) || ($GetOptions_H{"-h"})) {
    outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date, undef);
    opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
    if(! $options_okay) { die "ERROR, unrecognized option;"; }
    else                { exit 0; } # -h, exit with 0 status                                                                                                                                        
}

# check that number of command line args is correct                                                                       
if(scalar(@ARGV) != 2) {
    print "Incorrect number of command line arguments.\n";
    print $usage;
    print "\nTo see more help on available options, do dnaorg_evaluate_refseq_assign.pl -h\n\n";
    exit(1);
}
my ($results, $standards) = (@ARGV);

# set options in opt_HH                                                                                                                                                      
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)                                                                                                                                   
opt_ValidateSet(\%opt_HH, \@opt_order_A);


my $dir        = opt_Get("--dirout", \%opt_HH);          # this will be undefined unless -d set on cmdline                                                  

#############################                                                                                                                              
# create the output directory                                                                                                                                                 
#############################                                                                                                                                             
my $cmd;              # a command to run with runCommand()                                                                              
my @early_cmd_A = (); # array of commands we run before our log file is opened                                                           
# check if the $dir exists, and that it contains the files we need                          
# check if our output dir $symbol exists                                                 
if(! defined $dir) {
    $dir = $results . "-evaluation";
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
my @arg_desc_A = ("results of dnaorg_refseq_assign", "'gold standard' ntlists");
my @arg_A      = ($results, $standards);
outputBanner(*STDOUT, $version, $releasedate, $synopsis, $date, undef);
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
outputBanner($log_FH, $version, $releasedate, $synopsis, $date, undef);
opt_OutputPreamble($log_FH, \@arg_desc_A, \@arg_A, \%opt_HH, \@opt_order_A);

# output any commands we already executed to $log_FH                                                                                                                                                
foreach $cmd (@early_cmd_A) {
    print $cmd_FH $cmd . "\n";
}

my $do_keep = opt_Get("--keep", \%opt_HH); # should we leave intermediates files on disk, instead of removing them?                                                                  
my @files2rm_A = ();                       # will be filled with files to remove, --keep was not enabled                                                                                                                                         

#####################################################################################################
#####################################################################################################
#####################################################################################################

# Regex that matches a blank line
my $blank_match = '^\s*$';

# Generate data structures to build a sequence profile for further evaluation                            
my %hit_info_HAA = ();  # Hash of hash of arrays containing nhmmscan info for each hit for each sequence       
                        # 1D Key:    Sequence accession #                                                    
                        # 1D Array:  Order of hit (1st nhmmscan hit = index 0, 2nd = index 1, etc)                                                                                                                  
                        # Value:     Array w/ relevant nhmmscan info for the hit (see comments below)

# Find $results' $out_root
$dir_tail = $results;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $results_out_root = $results . "/" . $dir_tail;

# Find $standards' $out_root
$dir_tail = $standards;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $standards_out_root = $standards . "/" . $dir_tail;

# Variable which can be used w/ regexes to trim whitespace
my $r_trim = '\s*$';
my $l_trim = '^\s*';




# Find expected RefSeq for each sequence
my %standard_refseq_H = (); # Hash
                            # Key:   Sequence accession number
                            # Value: Proper RefSeq

my $refs_file = $results_out_root . ".all.refseqs";
open(REFS, "< $refs_file") || fileOpenFailure($refs_file, $0, $!, "reading", $ofile_info_HH{"FH"});

while(my $refseq = <REFS>) {
    chomp $refseq;
    $refseq =~ s/$l_trim//;
    $refseq =~ s/$r_trim//;

    my $standards_file = $standards_out_root . ".$refseq.ntlist";
    open(ST_NTLIST, "< $standards_file")  || fileOpenFailure($standards_file, $0, $!, "reading", $ofile_info_HH{"FH"});
    my @ntlist = <ST_NTLIST>;
    #shift(@ntlist); # get rid of RefSeq at the top of the file                                                                                                             
    foreach my $seq (@ntlist) {
        chomp $seq;
        $seq =~ s/$l_trim//;
        $seq =~ s/$r_trim//;

        $standard_refseq_H{$seq} = $refseq;
    }
}

#print Dumper(\%standard_refseq_H);

#############################################################################################
# Fill %hit_info_HAA for each sequence, and add specific info to the output file
#############################################################################################

# Obtain a list of all the sequnces
#my $seqs_file = $results_out_root . ".all.seqs";
#open(SEQS, $seqs_file) || fileOpenFailure($seqs_file, $0, $!, "reading", $ofile_info_HH{"FH"});


my @header = ("Query   ","RefSeq   ","Exp. RefSeq","Pass?","Bit score","E-val","Coverage","Bias","H2: RefSeq","Bit Diff","Covg. Diff, Num. Correct Hits \n");


my $match_info_file = $results_out_root . ".matches.info"; # match info file from results
my $eval_file = $out_root . ".evaluation"; # file to output evaluation to

# copy match info to output file
#$cmd = "cat $match_info_file > $eval_file";
#runCommand($cmd, opt_Get("-v", \%opt_HH), undef);

# open match info for reading, open evaluation file for writing
open(EVALUATION, "> $eval_file") || fileOpenFailure($eval_file, $0, $!, "writing", $ofile_info_HH{"FH"});
open(MATCH_TBL, $match_info_file) || fileOpenFailure($match_info_file, $0, $!, "reading", $ofile_info_HH{"FH"});
    
# debug
#$cmd = "cat $seq_tbl_file";
#runCommand($cmd, opt_Get("-v", \%opt_HH), undef);

my @hit_specs_A;
my @hit_output_A;

while(my $line = <MATCH_TBL>){
	chomp $line;
	if($line !~ m/^\#/ && $line !~ m/$blank_match/) {               #if this line is not commented

	    @hit_specs_A = split(/\s+/, $line);
	    my $query = $hit_specs_A[0];
	    my $refseq = $hit_specs_A[1];
	    
	    # ########################################################################
            # In @hit_output_A - array containing all the output for this hit
	    #
	    # index         feature
	    #--------------------------
	    # 0             Query accn #
	    # 1             RefSeq accn #
	    # 2             Expected RefSeq accn #
	    # 3             Is the assignment correct? P(pass) or F(fail)
	    # 4             Bit-score
	    # 5             E-value
	    # 6             Coverage (Hit length)/(Query length)
	    # 7             Bias
	    # 8             RefSeq accn # of second hit
	    # 9             Difference in bit scores for first and second hit
	    # 10            Difference in coverage for first and second hit
	    # 11            Number of Correct hits to the intended RefSeq
	    ###########################################################################

	    @hit_output_A = @hit_specs_A; # Array to be outputed to the evaluation file - contains the match info, and expected RefSeq info
	    
	    # think it would be best to splice stuff in here

	    # check if the correct refseq was assigned
      	    my $correct;
	    if(! exists $standard_refseq_H{$query}) {  # if this sequence is not included in the 'standard' ntlists, sets the 'standard' RefSeq to -------
		$standard_refseq_H{$query} = "-------";
	    }

	    if($refseq eq $standard_refseq_H{$query}){
		$correct = "P";
	    }else {
		$correct = "F";
	    }



	    splice(@hit_output_A, 2, 0, $standard_refseq_H{$query}, $correct);

	    print EVALUATION join("\t\t", @hit_output_A);
	    print EVALUATION "\n";

	    #debug
	    #print Dumper(\%hit_info_HAA);

	} elsif ( $line =~ m/^#H/ ) { # if this is the header line, inserts headers for the two added fields
	    my @header_A = split(/\t\t/, $line);
	    
	    splice(@header_A, 2, 0, "Exp RefSeq", "Pass?");
	    
	    print EVALUATION join("\t\t", @header_A);
	    print EVALUATION "\n";

	} else {      # if this was a commented or blank line, simply copy it into the output file
	    print EVALUATION $line . "\n";

	    # inserts information into the explanation table about the two added fields
	    if($line =~ m/^# RefSeq/) {
		print EVALUATION "# Exp RefSeq: The RefSeq that this sequence should have been assigned to\n";
		print EVALUATION "# Pass?: Was the sequence assigned to the right RefSeq? P (pass) or F (fail)\n";
	    }
	}
}

    
#print Dumper(\%hit_info_HAA);

addClosedFileToOutputInfo(\%ofile_info_HH, "evaluation", $eval_file, 1, "Table summarizing the success of dnaorg_refseq_assign.pl");

$total_seconds += secondsSinceEpoch();
outputConclusionAndCloseFiles($total_seconds, $dir, \%ofile_info_HH);
exit 0;



