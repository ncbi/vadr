#!/usr/bin/env perl
# EPN, Wed Aug  7 15:12:40 2024
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::MSA;
use Bio::Easel::SqFile;
use LWP::Simple; 
use Mozilla::CA;

require "vadr.pm"; 
require "sqp_opts.pm";
require "sqp_ofile.pm";
require "sqp_seqfile.pm";
require "sqp_utils.pm";

#######################################################################################
# What this script does: 
#
# - Reads in a vadr .tbl (or .ftr) file
# - Outputs a GFF file with the same information in the .tbl (or .ftr) file
# 
#######################################################################################

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
my $executable    = "annotate-tbl2gff.pl";
my $usage         = "Usage: $executable [-options]\n\t<path to v-annotate.pl output .tbl file>\n";
my $synopsis      = "$executable :: add a single protein to a VADR blastx protein database";
my $date          = scalar localtime();
my $version       = "1.6.4";
my $releasedate   = "Jun 2024";
my $pkgname       = "VADR";

# print help and exit if necessary
if((! $options_okay) || ($GetOptions_H{"-h"})) { 
  ofile_OutputBanner(*STDOUT, $pkgname, $version, $releasedate, $synopsis, $date, undef);
  opt_OutputHelp(*STDOUT, $usage, \%opt_HH, \@opt_order_A, \%opt_group_desc_H);
  if(! $options_okay) { die "ERROR, unrecognized option;"; }
  else                { exit 0; } # -h, exit with 0 status
}

# check that number of command line args is correct
if(scalar(@ARGV) != 1) {   
  print "Incorrect number of command line arguments.\n";
  print $usage;
  print "\nTo see more help on available options, do $executable -h\n\n";
  exit(1);
}
my ($in_tbl_file) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

#####################
# Input the tbl file 
#####################
my %ftr_info_HAH = (); # hash of array of hashes with feature info 

my $progress_w = 50;
my $start_secs = ofile_OutputProgressPrior("Parsing input feature table file", $progress_w, undef, *STDOUT);

my @seq_order_A = (); # array of sequence names read in order from feature table file
local_sqf_FeatureTableParse($in_tbl_file, \%ftr_info_HAH, \@seq_order_A, 0, undef);

#utl_HAHDump("ftr_info_HAH", \%ftr_info_HAH, *STDOUT);

FeatureInfoSetIdAndParentForGff(\%ftr_info_HAH, \@seq_order_A, undef);
GffOutput(\%ftr_info_HAH, \@seq_order_A, undef, \%opt_HH, undef);

##########
# Conclude
##########

exit 0;

#################################################################
# Subroutine: local_sqf_FeatureTableParse()
# Incept:     EPN, Thu Aug  8 10:41:39 2024
#             EPN, Mon May 20 09:36:09 2019
#
# Synopsis: Parse a INSDC feature table format file.
#
# Arguments:
#  $infile:         feature table file to parse
#  $ftr_info_HAHR:  feature information, filled here
#                   1D key: accession
#                   2D:     feature index
#                   3D key: qualifer, value: qualifier value
#  $seq_order_AR:   REF to array of sequence names, added to here if defined, can be undef
#  $do_short_accn:  '1' to store information using short accessions if possible, '0' not to
#                   e.g. short accession for 'gi|126364580|dbj|AB271840.1|' is 'AB271840.1'
#  $FH_HR:          REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       if we have trouble parsing the file
#             if $allow_incomplete is '1' and we read an incomplete feature
#
# Reference: https://www.ncbi.nlm.nih.gov/Sequin/table.html
#################################################################
sub local_sqf_FeatureTableParse { 
  my $sub_name = "local_sqf_FeatureTableParse";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $ftr_info_HAHR, $seq_order_AR, $do_short_accn, $FH_HR) = @_;

  open(IN, $infile) || ofile_FileOpenFailure($infile, $sub_name, $!, "reading", $FH_HR);

  my $long_accver = undef;   # full accession, e.g. ref|NC_000883.2|
  my $acc         = undef;   # accession, e.g. 'NC_000883.2'
  my $ver         = undef;   # version, e.g. '2'
  my $qname       = undef;   # a qualifier name,  e.g. 'organism'
  my $qval        = undef;   # a qualifier value, e.g. 'Paramecia bursaria Chlorella virus 1'
  my $feature     = undef;   # a feature name, e.g. "CDS", or "gene"
  my $ftr_idx     = -1;      # number of features read for current sequence
  my $coords      = undef;   # coordinates
  my $strand      = undef;   # strand of last segment read 
  my $trunc5      = undef;   # '1' if current feature is 5' truncated (start carrot, e.g. NC_031327:"<3281..4207")
  my $trunc3      = undef;   # '1' if current feature is 5' truncated (start carrot, e.g. "3281..>4207")
  my $line_idx    = 0;       # count of number of lines read in ftable
  my $prv_was_accn           = 0; # set to '1' if previous line was an accession line
  my $prv_was_coords_feature = 0; # set to '1' if previous line was a coordinates line with a feature name
  my $prv_was_coords_only    = 0; # set to '1' if previous line was a coordinates line without a feature name
  my $prv_was_quals          = 0; # set to '1' if previous line was a qualifier_name qualifier value line

  while(my $line = <IN>) { 
    $line_idx++;
    chomp $line;
    if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
    if($line =~ m/\w/) { 
      # parse each of the 6 line types differently
      # -------------------------------------------------------
      if($line =~ /^ERROR\:/)                                { ; } # special v-annotate.pl output line, ignore 
      elsif($line =~ /^Additional note\(s\) to submitter\:/) { ; } # special v-annotate.pl output line, ignore 
      elsif($line =~ /^\>Feature\s+(\S+)$/) { 
        # ACCESSION LINE
        # example:
        #>Feature ref|NC_001359.1|    
        # or
        #>Feature anyseqname
        $long_accver = $1;

        # store name
        if(defined $seq_order_AR) { push(@{$seq_order_AR}, $long_accver); }
        # accession line can occur after any other line type, so we don't have to check if line order makes sense for this case

        # if our previous line was coords_feature or coords_only, we need to store the feature from that previous line
        if(($prv_was_coords_feature) || ($prv_was_coords_only)) { 
          $ftr_idx++;
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",   $feature, $FH_HR);
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "coords", $coords,  $FH_HR);
          if($trunc5) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc5", 1, $FH_HR); }
          if($trunc3) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc3", 1, $FH_HR); }
        }

        # determine accession and version, e.g. NC_001359.1 in above example
        if($do_short_accn && ($long_accver =~ /[^\|]*\|([^\|]+)\.(\d+)\|/)) { 
          $acc = $1;
          $ver = $2;
        }
        else { 
          $acc = $long_accver;
          #ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unable to parse header line:\n$line\n", 1, $FH_HR);
        }
        @{$ftr_info_HAHR->{$acc}} = (); # initialize array
        $feature = undef; 
        $coords  = undef;
        $trunc5  = undef;
        $trunc3  = undef;
        $strand  = undef;
        $ftr_idx = -1;

        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 1;
        $prv_was_coords_feature = 0;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 0;
        #printf("set prv_was_accn\n");
      }
      # -------------------------------------------------------
      elsif($line =~ /^(\<?)(\d+\^?)\t(\>?)(\d+)\t(\S+)$/) { 
        # COORDINATES LINE WITH A FEATURE NAME (coords_feature)
        # example:
        # 230   985     gene
        my ($start_carrot, $start_coord, $stop_carrot, $stop_coord, $tmp_feature) = ($1, $2, $3, $4, $5);
        # coords_feature line can occur after any other line type, so we don't have to check if line order makes sense for this case

        # if our previous line was coords_feature or coords_only, we need to store the feature from that previous line
        if(($prv_was_coords_feature) || ($prv_was_coords_only)) { 
          $ftr_idx++;
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",   $feature, $FH_HR);
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "coords", $coords,  $FH_HR);
          if($trunc5) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc5", 1, $FH_HR); }
          if($trunc3) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc3", 1, $FH_HR); }
        }
        $feature = $tmp_feature;
        $trunc5 = ($start_carrot eq "<") ? 1 : 0;
        $trunc3 = ($stop_carrot  eq ">") ? 1 : 0;

        $strand = ($start_coord <= $stop_coord) ? "+" : "-"; # sets + strand for single nt
        $coords = $start_coord . ".." . $stop_coord . ":" . $strand;

        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 0;
        $prv_was_coords_feature = 1;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 0;
        #printf("set prv_was_coords_feature\n");
      }
      # -------------------------------------------------------
      elsif($line =~ /^(\<?)(\d+)\t(\>?)(\d+)$/) {  
        # COORDINATES LINE WITHOUT A FEATURE NAME (coords_only) 
        # example:
        # 154   183
        my ($start_carrot, $start_coord, $stop_carrot, $stop_coord) = ($1, $2, $3, $4);

        # a coords_only line can only occur after a coords_feature line or coords_only line, 
        # check to make sure that's the case
        if($prv_was_coords_feature || # previous line was a coords line with a feature (common)
           $prv_was_coords_only) {    # previous line was a coords line without a feature (common)
          # line order makes sense, keep going...

          if($start_carrot ne "") { 
            ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read start_carrot indicating 5' truncated feature from coordinate line that wasn't the first coordinate line for the feature, line:\n$line\n", 1, $FH_HR);
          }
          $trunc3 = ($stop_carrot eq ">") ? 1 : 0;

          if($start_coord == $stop_coord) { 
            ; # single nt, use strand of previous segment
          }
          else { 
            $strand = ($start_coord < $stop_coord) ? "+" : "-";
          }
          $coords .= "," . $start_coord . ".." . $stop_coord . ":" . $strand;

          # update '$prv_*' values that we use to make sure line order makes sense
          $prv_was_accn           = 0;
          $prv_was_coords_feature = 0;
          $prv_was_coords_only    = 1;
          $prv_was_quals          = 0;
          #printf("set prv_was_coords_only\n");
        }
        else { # line order is unexpected
          ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unexpected line order (coords_only line), line:\n$line\n", 1, $FH_HR);
        }
      }
      # -------------------------------------------------------
      elsif(($line =~ /^\t\t\t[^\t]+\t[^\t]+$/) || 
            ($line =~ /^\t\t\t[^\t]+$/)) { 
        # QUALIFIER LINE
        # examples:
        #gene       AR1
        #locus_tag  PhyvvsAgp1

        # before parsing it, do two sanity checks
        if(! defined $coords)  { 
          ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read qualifier line but coords is not yet defined, line:\n$line\n", 1, $FH_HR);
        }
        if(! defined $feature)  { 
          ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, read qualifier line but feature is not yet defined, line:\n$line\n", 1, $FH_HR);
        }
        # does line order make sense?
        if($prv_was_coords_feature || 
           $prv_was_coords_only    ||
           $prv_was_quals) { 
          # line order makes sense, keep going...
          if(! $prv_was_quals) { 
            # first quals line for this feature, store the feature
            $ftr_idx++;
            sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",   $feature, $FH_HR);
            sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "coords", $coords,  $FH_HR);
            if($trunc5) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc5", 1, $FH_HR); }
            if($trunc3) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc3", 1, $FH_HR); }
          }
          # parse the line
          if($line =~ /^\t\t\t([^\t]+)\t([^\t]+)$/) { 
            ($qname, $qval) = ($1, $2);
          }
          elsif($line =~ /^\t\t\t([^\t]+)$/) { 
            ($qname, $qval) = ($1, "");
          }
          else { 
            ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unable to parse qualifier line, line:\n$line\n", 1, $FH_HR);
          }
          sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, $qname, $qval, $FH_HR);
        } # end of 'if() that checks line order makes sense
        else { # unexpected line order
          ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unexpected line order (quals line), line:\n$line\n", 1, $FH_HR);
        }          

        # update '$prv_*' values that we use to make sure line order makes sense
        $prv_was_accn           = 0;
        $prv_was_coords_feature = 0;
        $prv_was_coords_only    = 0;
        $prv_was_quals          = 1;
        #printf("set prv_was_quals\n");
      }
      # -------------------------------------------------------
      else { 
        ofile_FAIL("ERROR in $sub_name, problem parsing $infile at line $line_idx, unable to parse line, line:\n$line\n", 1, $FH_HR);
      }
      # -------------------------------------------------------
    }
  }
  if(! defined $acc) { 
    ofile_FAIL("ERROR in $sub_name, problem parsing $infile, did not read any accession lines\n", 1, $FH_HR);
  }

  # add the final feature, if we haven't already, we can tell based on previous line type
  if(($prv_was_coords_feature) || ($prv_was_coords_only)) { 
    $ftr_idx++;
    sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "type",   $feature, $FH_HR);
    sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "coords", $coords,  $FH_HR);
    if($trunc5) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc5", 1, $FH_HR); }
    if($trunc3) { sqf_StoreQualifierValue(\@{$ftr_info_HAHR->{$acc}}, $ftr_idx, "trunc3", 1, $FH_HR); }
  }

  return;
}

#################################################################
# Subroutine: FeatureInfoSetIdAndParentForGff()
# Incept:     EPN, Thu Aug  8 14:07:34 2024
#
# Synopsis: Given feature info, further populate it with ID and
#           parent information in preparation for GFF output. Parent 
#           information is determined using the protein_id key. 
#
# Arguments:
#  $ftr_info_HAHR: feature info
#  $seq_order_AR:  ref to array with order of sequences to output, if undef, output in sorted order
#  $FH_HR:         ref to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       never
#################################################################
sub FeatureInfoSetIdAndParentForGff { 
  my $sub_name = "FeatureInfoSetIdAndParentForGff";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_HAHR, $seq_order_AR, $FH_HR) = @_;
  my @seq_A = ();
  if(! defined $seq_order_AR) {
    @seq_A = sort keys %{$ftr_info_HAHR};
  }
  else {
    @seq_A = @{$seq_order_AR};
  }
  
  # first pass through to determine ID value and
  # populate %protein_id2cds_id_H
  my %protein_id2cds_id_H = (); # key is protein_id of a CDS, value is that CDS' ID value
  my %ftr_id_idx_H = (); # key is feature name, value is index for this feature for current accn
  foreach my $seq (@seq_order_A) {
    if(! defined $ftr_info_HAHR->{$seq}) {
      ofile_FAIL("ERROR in $sub_name, no feature information for sequence $seq", 1, $FH_HR);
    }
    my $nftr = scalar(@{$ftr_info_HAHR->{$seq}});
    
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      my $ftr_type = $ftr_info_HAHR->{$seq}[$ftr_idx]->{"type"};   
      if(! defined $ftr_id_idx_H{$ftr_type}) {
        $ftr_id_idx_H{$ftr_type} = 1;
      }
      else {
        $ftr_id_idx_H{$ftr_type}++;
      }
      $ftr_info_HAHR->{$seq}[$ftr_idx]{"GFF_ID"} = $ftr_type . $ftr_id_idx_H{$ftr_type};
      # keep track of the id for this protein id, so we can set parent indices for 
      if((vdr_FeatureTypeIsCds($ftr_info_HAHR->{$seq}, $ftr_idx)) &&
         (defined $ftr_info_HAHR->{$seq}[$ftr_idx]{"protein_id"})) {
        $protein_id2cds_id_H{$ftr_info_HAHR->{$seq}[$ftr_idx]{"protein_id"}} = $ftr_info_HAHR->{$seq}[$ftr_idx]{"GFF_ID"};
        printf("HEYA set protein_id2cds_id_H{ftr_info_HAHR->{$seq}[$ftr_idx]{protein_id}} to " . $ftr_info_HAHR->{$seq}[$ftr_idx]{"GFF_ID"} .  "\n");
      }
    }
  }
  
  # second pass through to set parent values for mat_peptide and sig_peptide features
  foreach my $seq (@seq_order_A) {
    my $nftr = scalar(@{$ftr_info_HAHR->{$seq}});
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if(($ftr_info_HAHR->{$seq}[$ftr_idx]{"type"} eq "mat_peptide") ||
         ($ftr_info_HAHR->{$seq}[$ftr_idx]{"type"} eq "sig_peptide") && 
         (defined $ftr_info_HAHR->{$seq}[$ftr_idx]{"protein_id"})) {
        my $protein_id = $ftr_info_HAHR->{$seq}[$ftr_idx]{"protein_id"};
        if(defined $protein_id2cds_id_H{$protein_id}) {
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"GFF_Parent"} = $protein_id2cds_id_H{$protein_id};
          printf("HEYA set ftr_info_HAHR->{$seq}[$ftr_idx]{GFF_Parent} to $protein_id2cds_id_H{$protein_id}\n");
        }
      }
    }
  }    
  return;
}
  
#################################################################
# Subroutine: GffOutput()
# Incept:     EPN, Wed Aug  7 15:35:24 2024
#
# Synopsis: Output a GFF file given a ftr_info_HAH.
#           Format ref: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#
# Arguments:
#  $ftr_info_HAHR: feature info
#  $seq_order_AR:  ref to array with order of sequences to output, if undef, output in sorted order
#  $out_FH:        file handle to output to 
#  $opt_HHR:       REF to 2D hash of option values, see top of epn-options.pm for description, PRE-FILLED
#  $FH_HR:         REF to hash of file handles, including "log" and "cmd"
#
# Returns:    void
#
# Dies:       never
#################################################################
sub GffOutput { 
  my $sub_name = "GffOutput";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_HAHR, $seq_order_AR, $out_FH, $opt_HHR, $FH_HR) = @_;
  my @seq_A = ();
  if(! defined $seq_order_AR) {
    @seq_A = sort keys %{$ftr_info_HAHR};
  }
  else {
    @seq_A = @{$seq_order_AR};
  }
  
  my $source = "VADR:v" . $version . ":v-annotate.pl";
  foreach my $seq (@seq_order_A) {
    if(! defined $ftr_info_HAHR->{$seq}) {
      ofile_FAIL("ERROR in $sub_name, no feature information for sequence $seq", 1, $FH_HR);
    }
    my $nftr = scalar(@{$ftr_info_HAHR->{$seq}});
    
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      my $ftr_type = $ftr_info_HAHR->{$seq}[$ftr_idx]->{"type"};   
      my @start_A  = ();
      my @stop_A   = ();
      my @strand_A = ();
      vdr_FeatureStartStopStrandArrays($ftr_info_HAHR->{$seq}[$ftr_idx]{"coords"}, \@start_A, \@stop_A, \@strand_A, $FH_HR);
      my $nsgm = scalar(@start_A);
      my @ftr_lines_A = (); # we'll store the feature lines here, and output them in forward order or reverse order depending on strand
      my $summary_strand = vdr_FeatureSummaryStrand($ftr_info_HAHR->{$seq}[$ftr_idx]{"coords"}, $FH_HR);
      my $id = $ftr_info_HAHR->{$seq}[$ftr_idx]{"GFF_ID"};
      if(! defined $id) {
        ofile_FAIL("ERROR in $sub_name, no ID set for feature for sequence $seq", 1, $FH_HR);
      }
      my $key_values = "";
      foreach my $key (sort keys(%{$ftr_info_HAHR->{$seq}[$ftr_idx]})) { 
        if(($key ne "GFF_ID") && ($key ne "coords") && ($key ne "type") && ($key ne "protein_id")) { 
          my $key2print = $key;
          $key2print =~ s/^GFF\_//;
          $key_values .= $key2print . "=" . $ftr_info_HAHR->{$seq}[$ftr_idx]{$key} . ";";
        }
      }
      for(my $sgm_idx = 0; $sgm_idx < $nsgm; $sgm_idx++) {
        my ($start, $stop, $strand) = (undef, undef, undef);
        if($strand_A[$sgm_idx] eq "+") {
          ($start, $stop, $strand) = ($start_A[$sgm_idx], $stop_A[$sgm_idx], $strand_A[$sgm_idx]);
        }
        else {
          ($stop, $start, $strand) = ($start_A[$sgm_idx], $stop_A[$sgm_idx], $strand_A[$sgm_idx]);
        }
        my $attributes = "ID:" . $id . ";" . $key_values;
        push(@ftr_lines_A,
             sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
                     $seq,                                 # token 1: 'sequence' (sequence name)
                     $source,                              # token 2: 'source'
                     $ftr_type,                            # token 3: 'feature' (model name) you may want to change this to 'ncRNA'
                     $start,                               # token 4: 'start' in coordinate space [1..seqlen], must be <= 'end'
                     $stop,                                # token 5: 'end' in coordinate space [1..seqlen], must be >= 'start'
                     ".",                                  # token 6: 'score' 
                     $strand,                              # token 7: 'strand' ('+' or '-')
                     ".",                                  # token 8: 'phase' irrelevant for noncoding RNAs
                     $attributes));                         # token 9: attributes, currently only E-value, unless --all, --none or --desc
      } # end of 'for(my $sgm_idx = 0'...

      # output lines for this feature
      my $line_idx; 
      if($summary_strand ne "-") { # summary_strand is -
        for($line_idx = 0; $line_idx < scalar(@ftr_lines_A); $line_idx++) {
          print $ftr_lines_A[$line_idx];
        }
      }
      else { # summary_strand is -
        for($line_idx = (scalar(@ftr_lines_A)-1); $line_idx >= 0; $line_idx--) {
          print $ftr_lines_A[$line_idx];
        }
      }
    }
  }
  
  return;
}

