#!/usr/bin/env perl
# EPN, Mon Aug 10 10:39:33 2015 [development began on dnaorg_annotate_genomes.pl]
# EPN, Mon Feb  1 15:07:43 2016 [dnaorg_build.pl split off from dnaorg_annotate_genomes.pl]
# LES, Wed Jul  6 15:10    2016 [dnaorg_get_matpepts.pl split off from dnaorg_build.pl]
 
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;
use Data::Dumper;

require "dnaorg.pm"; 
require "epn-options.pm";

#######################################################################################
# What this script does: 
#
# Preliminaries: 
#   - process options
#   - create the output directory
#   - output program banner and open output files
#   - parse optional input files, if necessary
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
my $inf_exec_dir      = $dnaorgdir . "/infernal-1.1.2/src/";
my $hmmer_exec_dir    = $dnaorgdir . "/hmmer-3.1b2/src/";
my $esl_exec_dir      = $dnaorgdir . "/infernal-1.1.2/easel/miniapps/";
my $esl_fetch_cds     = $dnaorgdir . "/esl-fetch-cds/esl-fetch-cds.pl";

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
opt_Add("-f",           "boolean", 0,                        1,    undef, undef,      "forcing directory overwrite",           "force; if dir <reference accession>.get_mat_pept_file exists, overwrite it", \%opt_HH, \@opt_order_A);
opt_Add("-v",           "boolean", 0,                        1,    undef, undef,      "be verbose",                            "be verbose; output commands to stdout as they're run", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $usage    = "Usage: dnaorg_build.pl [-options] <reference accession>\n";
my $synopsis = "dnaorg_build.pl :: build homology models for features of a reference sequence";

my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                'v'            => \$GetOptions_H{"-v"},
    );

my $total_seconds = -1 * secondsSinceEpoch(); # by multiplying by -1, we can just add another secondsSinceEpoch call at end to get total time
my $executable    = $0;
my $date          = scalar localtime();
my $version       = "0.23";
my $releasedate   = "Jan 2018";

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
    exit(1);
}
my ($ref_accn) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);


#############################
# create the output directory
#############################
my $cmd;              # a command to run with runCommand()
my @early_cmd_A = (); # array of commands we run before our log file is opened
# check if the $dir exists, and that it contains the files we need
# check if our output dir $symbol exists
my $dir; #originally $dir could be defined through a command line option; I've modified it so that it always has to be defined by $ref_accn
if(! defined $dir) { 
  $dir = $ref_accn . ".get_mat_pept_file"; #Makes output directory accn#.get_mat_peptide_file
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
# changed opt_Get to 1 (be verbose)
runCommand($cmd, 1, undef);
push(@early_cmd_A, $cmd);

my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir - this is the file folder that all output files will be in
my $out_root = $dir . "/" . $dir_tail . ".dnaorg_get_matpepts";

#############################################
# output program banner and open output files
#############################################

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


# output any commands we already executed to $log_FH
foreach $cmd (@early_cmd_A) { 
  print $cmd_FH $cmd . "\n";
}



###########################################################################
# Step 1. Gather and process information on reference genome using Edirect.
###########################################################################

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
mat_pept_wrapperGetInfoUsingEdirect(undef, $ref_accn, $out_root, \%cds_tbl_HHA, \%mp_tbl_HHA, \%ofile_info_HH, $ofile_info_HH{"FH"}); # 1st argument is undef because we are only getting info for $ref_accn


########### - DEBUGGING - #############
# Print out contents of mp_tbl_HHA and cds_tbl_HHA
#print "CDS tbl_HHA:";
#print Dumper(\%cds_tbl_HHA);
#print "Mat Pept tbl_HHA:";
#print Dumper(\%mp_tbl_HHA);

#print $cds_tbl_HHA{$ref_accn}{'coords'}[0];



########################################################################
# Create matpept.in file
########################################################################

# first array contains the starting coordinate of each CDS/matpept, in order of retrieval - second array contains corresponding stop coordinate, third contains number of exons for that CDS
my @cds_start_coords = ();
my @cds_end_coords = ();
my @cds_num_exons = ();

my @mp_start_coords = (); 
my @mp_end_coords = ();
my @mp_num_exons = ();

# Array which contains the skipped over/repeated part from a join() if there is one in the CDS/matpept. Otherwise that index contains a 0
my @cds_interruptions = ();
my @mp_interruptions = ();

# fills the 4 cds arrays, then fills the 4 mp arrays
foreach my $coord_set (@{$cds_tbl_HHA{$ref_accn}{'coords'}}) {
    mat_pept_startsStopsFromCoords($coord_set,\@cds_start_coords,\@cds_end_coords,\@cds_num_exons, \@cds_interruptions, $ofile_info_HH{"FH"}); 
}
foreach my $coord_set (@{$mp_tbl_HHA{$ref_accn}{'coords'}}) {
    mat_pept_startsStopsFromCoords($coord_set,\@mp_start_coords,\@mp_end_coords,\@mp_num_exons, \@mp_interruptions, $ofile_info_HH{"FH"});
}


###DEBUGGING###
#print "CDS coords: start, end, num_exons, interruptions: ";
#print Dumper(\@cds_start_coords);
#print Dumper(\@cds_end_coords);
#print Dumper(\@cds_num_exons);
#print Dumper(\@cds_interruptions);

#print "MP coords: start, end, num_exons, interruptions: ";
#print Dumper(\@mp_start_coords);
#print Dumper(\@mp_end_coords);
#print Dumper(\@mp_num_exons);
#print Dumper(\@mp_interruptions);
###############

#print scalar(@cds_start_coords);

# Arrays of arrays which store @primary/@all arrays for each CDS
my @cds_primary_AA = (); #1D - CDS#
                         #2D - @primary (array of primary matpepts in that cds)
my @cds_all_AA = ();




# iterate through each cds (i = the cds's # in the matpept file - 1)
for (my $i=0; $i<scalar(@cds_start_coords); $i++) {
    my $start_cds = $cds_start_coords[$i];
    my $end_cds = $cds_end_coords[$i];
    my @primary = (); # array of the indexes of the primary matpepts in this cds (values are 1 less than the matpept file indices)
    my @all = (); # indexes of all matpepts in this cds
    my $curr_primary = undef; # current primary matpept

    # check for interruptions
    my $interrupt_mp = undef; # index of matpept which is interrupted
    if($cds_interruptions[$i] ne "0") {
	# finds the index of the matpept which is interrupted
	for(my $k=0; $k<scalar(@mp_interruptions); $k++) {           
	    if($mp_interruptions[$k] eq $cds_interruptions[$i]) { #if this matpept has the same interruption as the cds
		$interrupt_mp = $k;
		last;
	    }
	}
    }

    # check to make sure that interrupted matpept is the last matpept in the cds                                                                                                   
    if(defined($interrupt_mp)  &&  $mp_end_coords[$interrupt_mp] != ($end_cds - 3) ) {   
        die "Error in $0 main: There is a frame-shifted matpept which is not the last matpept of its CDS. There is not yet code for this case"; }   


    # fill up @primary and @all for this cds
    for (my $j=0; $j<scalar(@mp_start_coords); $j++) {
        # if first mp; assuming here that the primary matpept is the first matpept listed
	# also checks to make sure that if the CDS is interrupted and this matpept starts at the same point as the interrupted matpept, doesn't start at this matpept (waits so that it can start at the interrupted matpept)
        if($mp_start_coords[$j] == $start_cds && !defined($primary[0])  &&  (!defined($interrupt_mp) || $mp_start_coords[$j] != $mp_start_coords[$interrupt_mp] || $j == $interrupt_mp) ) { 
	    $primary[0] = $j;
	    $all[0] = $j;
	    $curr_primary = $j;

	}elsif(defined($primary[0])) { # if not the first mp	   
	    # add matpept to all, and primary if appropriate
	    push(@all, $j);

	    if($mp_start_coords[$j] > $mp_end_coords[$curr_primary]) { # if primary (if this mp starts after the last primary matpept ended)
		# check to make sure primary matpepts are properly aligned
		unless($mp_start_coords[$j] == ($mp_end_coords[$curr_primary] + 1)) {
		    die "Error in $0 main at matpept @{[$j+1]}: Primary matpepts not properly aligned."; }
		
		push(@primary, $j);
		$curr_primary = $j;

       	     }else { # if not primary
		 # check that matpept ends are aligned
		 unless($curr_primary == ($j-1)  ||  $mp_start_coords[$j] == ($mp_end_coords[$j-1] + 1) ) {
		     die "Error in $0 main at matpept @{[$j+1]}: Matpepts not properly aligned."; }
		 
	     }

	}
	
	 # check for last matpept in cds
	 # if the end of this matpept is the end of the cds, or if this is the matpept before the interrupted matpept or if this is the interrupted matpept (happens when entire CDS is the interrupted matpept)
	 if($mp_end_coords[$j] == ($end_cds-3)   ||   (defined($interrupt_mp) &&  $mp_end_coords[$j] == $mp_start_coords[$interrupt_mp]-1 )){  
	     #print "ENTERED for $i at mp $j \n";
	     # if this is the last matpept, or if the next matpept starts after the end of cds (checks to make sure this is not just the last primary matpept, but the last all matpept)
	    if( ($j+1) == scalar(@mp_end_coords)  ||  $mp_start_coords[$j+1] > $end_cds  ||  $mp_start_coords[$j+1]  ) {
		#print "ENTERED for $i at mp $j \n";
		last;
	    }
	 }
    }

    # if there was an interruption, adds last, interrupted mp
    if(defined($interrupt_mp) && $primary[-1] != $interrupt_mp) {
	push(@primary, $interrupt_mp);
	push(@all, $interrupt_mp);
    }

    # add @primary and @all to the appropriate CDS index in @cds_primary_AA
    push(@cds_primary_AA, \@primary);
    push(@cds_all_AA, \@all);
}


##########- DEBUGGING - ####################
#print "Primary mp list: ";
#print Dumper(\@cds_primary_AA);
#print "All mp list: ";
#print Dumper(\@cds_all_AA);

#####################################
#Format and output the matpept.in file

open MPIN, "> $ref_accn.matpept.in"
    or die "Cannot open a file to put matpept.in info in";

print MPIN "# This file explains how CDS and mat_peptide annotation for $ref_accn are related\n";  #TODO: add $ref_accn scientific name
print MPIN "#\n";
print MPIN "# Format of lines in this file:\n";
print MPIN "#<CDS-idx> <'primary' OR 'all'> <mat_peptide-1-idx>:<mat_peptide-2-idx>:<mat_peptide-n-idx>\n";
print MPIN "#\n";
print MPIN "# 'primary' lines: these define the 'primary' peptides in order, the CDS <CDS-idx> is comprised\n";
print MPIN "#                  of the peptides listed in final token, which are contiguous, start of first\n";
print MPIN "#                  mat_peptide to stop of final mat_peptide is one contiguous subsequence.\n";
print MPIN "# 'all' lines:     these define all the peptides that are ultimately derived from CDS <CDS-idx>.\n";
print MPIN "#                  It will be a superset of the primary line for this index but will additionally\n";
print MPIN "#                  include mat_peptides that are secondarily cleaved from the primary mat_peptides.\n";
print MPIN "#\n";
print MPIN "#\n";

# for each cds
for(my $cds=0; $cds < scalar(@cds_primary_AA); $cds++) {
    # print out list of primary matpepts
    print MPIN "" . $cds+1 . " primary\t";
    for (my $i=0; $i< scalar(@{$cds_primary_AA[$cds]}); $i++){
	print MPIN @{$cds_primary_AA[$cds]}[$i] + 1;

	if($i != scalar(@{$cds_primary_AA[$cds]})-1)  {print MPIN ":";}
    }
    print MPIN "\n";

    #print out list of all matpepts
    print MPIN "" . $cds+1 . " all\t\t";
    for (my $i=0; $i < scalar(@{$cds_all_AA[$cds]}); $i++){
	print MPIN @{$cds_all_AA[$cds]}[$i] + 1;

	if($i != scalar(@{$cds_all_AA[$cds]})-1)  {print MPIN ":";};
    }
    print MPIN "\n";


}
close(MPIN);


##########
# Conclude
##########

# close any open file handles                                                                                                                               
my $key2d;
foreach $key2d (keys ($ofile_info_HH{"FH"})) {
  if(defined $ofile_info_HH{"FH"}{$key2d}) {
    close $ofile_info_HH{"FH"}{$key2d};
  }
}

exit 0;

################################################################################
# Wrapper sub which collects info from Edirect
#
#
# LES 7/08/16
# Adapted from wrapperGetInfoUsingEdirect in dnaorg.pm
#
#
################################################################################
sub mat_pept_wrapperGetInfoUsingEdirect {
    my $sub_name = "mat_pept_wrapperGetInfoUsingEdirect";
    my $nargs_expected = 7;
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name, entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

    my ($listfile, $ref_accn, $out_root, $cds_tbl_HHAR, $mp_tbl_HHAR, $ofile_info_HHR, $FH_HR) = @_;

    my $cmd; # a command to be run by runCommand()                                                                                                                                                                                                  

###
    # listfile comes in undef, so $have_listfile is always 0
    my $have_listfile = (defined $listfile) ? 1 : 0;
###

  # We create the .mat_peptide file first because we will die with an                                                                                                                                                                             
  # error if mature peptide info exists and neither --matpept nor                                                                                                                                                                                 
  # --nomatpept was used (and we want to die as early as possible in the                                                                                                                                                                          
  # script to save the user's time)                                                                                                                                                                                                               
  #                                                                                                                                                                                                                                               
  # 1) create the edirect .mat_peptide file                                                                                                                                                                                      
    my $mp_file = $out_root . ".mat_peptide";	
	
        $cmd = "esearch -db nuccore -query $ref_accn";
       	$cmd .= " | efetch -format gpc | xtract -insd mat_peptide INSDFeature_location product > $mp_file";
        # right now I have second arg ($be_verbose) set to 1, meaning that the name of the cmd will be printed to stdout (change to 0 to not print)
	runCommand($cmd, 1, $FH_HR);

	# Checks to make sure edirect does have matpept info
	if(! -s  $mp_file) {
	    DNAORG_FAIL("ERROR, in $sub_name, --matpept enabled but no mature peptide information exists.", 1, $FH_HR);
	}

        addClosedFileToOutputInfo($ofile_info_HHR, "mp", $mp_file, 1, "Mature peptide information obtained via edirect");
		 

  # 2) create the edirect .ftable file
  # create the edirect ftable file                                                                                                                                                                                                              
    my $ft_file  = $out_root . ".ftable";
    
        $cmd = "esearch -db nuccore -query $ref_accn";
	$cmd .= " | efetch -format ft > $ft_file";
        # hard-coded one makes the program print out the name of the command as it is run 
	runCommand($cmd, 1, $FH_HR);
	addClosedFileToOutputInfo($ofile_info_HHR, "ft", $ft_file, 1, "Feature table obtained via edirect");
    

  # 4) parse the edirect .mat_peptide file
    	edirectFtableOrMatPept2SingleFeatureTableInfo($mp_file, 1, "mat_peptide", $mp_tbl_HHAR, $FH_HR); # 1: it is a mat_peptide file                                                                                                                    
    
	if (! exists ($mp_tbl_HHAR->{$ref_accn})) {
	    DNAORG_FAIL("ERROR in $sub_name, --matpept enabled, but no mature peptide information stored for reference accession", 1, $FH_HR);
	}
    
  # 5) parse the edirect .ftable file                                                                                                                                                                                                                     
    edirectFtableOrMatPept2SingleFeatureTableInfo($ft_file, 0, "CDS", $cds_tbl_HHAR, $FH_HR); # 0: it's not a mat_peptide file                                                                                                                            
   if(! exists ($cds_tbl_HHAR->{$ref_accn})) {
	DNAORG_FAIL("ERROR in $sub_name, no CDS information stored for reference accession", 1, $FH_HR);
    }

}

##################################################################
#
# Subroutine: mat_pept_startsStopsFromCoords()
# Incept:     EPN, Thu Feb 11 14:22:54 2016                                                                                                                     
#             Copied by LES into dnaorg_get_matpepts from dnaorg.pm, and modified, Tue, Jul 12
#                                                                                                                                                      
# Purpose:    Extract the starts and stops from a coords string.                                                                               
#                                                                                                                                                 
# Args:                                                                                                                                                                
#   $coords:    the coords string                                                                                                                                                                     
#   $starts_AR: REF to array to fill with start positions, FILLED HERE                                                                                                                                
#   $stops_AR:  REF to array to fill with stop positions, FILLED HERE                                                                                                                  
#   $nexons_R:  REF to scalar that fill with the number of exons, FILLED HERE                                                                                                                         
#   $FH_HR:     REF to hash of file handles, including "log" and "cmd"                                                                                                                                    
#                                                                                                                                                                                    
# Returns:      void; but fills @{$starts_AR}, @{$stops_AR}, and $$nexons_R.                                                                                                                                 
#                                                                                                                                                                                         
# Dies:         if we can't parse $coords                                                                                                                                                                   
#################################################################                                                                                                                                   
sub mat_pept_startsStopsFromCoords {
    my $sub_name = "mat_pept_startsStopsFromCoords()";
    my $nargs_expected = 6;
    if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); }

    my ($coords, $starts_AR, $stops_AR, $nexons_AR, $interruption_AR, $FH_HR) = @_;

    my $orig_coords = $coords;
    # Examples:                                                                                                                                                                                                       
    # complement(2173412..2176090)                                                                                                                                                                                                                       
    # complement(join(226623..226774, 226854..229725))    
    # remove 'complement('  ')'                                                                                                                                                                                                                   
    $coords =~ s/^complement\(//;
    $coords =~ s/\)$//;

    # remove 'join('  ')'                                                                                                                                                                                                                                   
    $coords =~ s/^join\(//;
    $coords =~ s/\)$//;

    my @el_A = split(/\s*\,\s*/, $coords);
    push(@{$nexons_AR}, scalar(@el_A));

    # if there are multiple frameshifts for one feature, die                                                                                                                             
    if(@{$nexons_AR}[-1] > 2) {
        die "ERROR in $sub_name, more than two exons in a feature."; }


    # if there are more than one exons, marks an interruption, and then changes $coords so that it only shows the start and end of the entire feature (no interruptions)
    if(scalar(@el_A)>1){
	# extracts the end coord of the first exon for $interruption, and the beg coord of the first exon for $coords
        $el_A[0] =~ m/(\d*)..(\d*)/;
       	# $interruption will be a string of the format 'num..num'
	my $interruption = $2;
	$coords = $1;

	# extracts the begginning coord of the second exon for $interruption, and the end coord of the second exon for $coords
	$el_A[1] =~ m/(\d*)..(\d*)/;
	$interruption .= ".." . $1;
	$coords .= ".." . $2;

	# adds the $interruption to that array
	push(@{$interruption_AR}, $interruption); 

	
    }elsif(scalar(@el_A)==1){
	push(@{$interruption_AR}, 0);

    }else {
	die "ERROR in $sub_name, no exons found";
    }

    # adds the stop and start coords to the appropriate arrays
    if($coords =~ m/^(\d+)\.\.(\d+)$/) {
	    push(@{$starts_AR}, $1);
	    push(@{$stops_AR},  $2);
    }
    elsif($coords =~ m/^(\d+)$/) { # a single nucleotide                                                                                                                                                                                              
	    push(@{$starts_AR}, $1);
	    push(@{$stops_AR},  $1);
    }
    else {
	    DNAORG_FAIL("ERROR in $sub_name, unable to parse coordinates $orig_coords", 1, $FH_HR);
    }
    

    return;
}
