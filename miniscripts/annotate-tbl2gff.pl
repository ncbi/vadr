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
my $g = 1;
$opt_group_desc_H{$g} = "basic options";
#     option            type       default  group   requires incompat    preamble-output                                                          help-output    
opt_Add("-h",           "boolean", 0,          $g,    undef, undef,      undef,                                                                   "display this help",                                   \%opt_HH, \@opt_order_A);
opt_Add("-s",           "string",  undef,      $g,    undef, undef,      "define source field as <s>",                                            "define source field as <s>",                          \%opt_HH, \@opt_order_A);
opt_Add("--ftr",        "boolean",  0,         $g,    undef, undef,      "input file is a .ftr file not a .tbl file",                             "input file is a .ftr file not a .tbl file",           \%opt_HH, \@opt_order_A);
$opt_group_desc_H{++$g} = "options for controlling which features or qualifiers are output in attributes field";
opt_Add("--qall",       "boolean",  0,         $g,    undef,  undef,     "output info for all qualifiers (except those in --qskip)",              "output info for all qualifiers (except those in --qskip)", \%opt_HH, \@opt_order_A);
opt_Add("--fskip",      "string",   undef,     $g,    undef,  undef,     "do not output features in comma separated string <s>",                  "do not output features in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--qskip",      "string",   undef,     $g,    undef,  undef,     "do not output qualifiers in comma separated string <s>",                "do not output qualifiers in comma separated string <s>", \%opt_HH, \@opt_order_A);
opt_Add("--noaddgene",  "boolean",  0,         $g,    undef,  undef,     "do not add gene qualifiers from gene features to overlapping features", "do not add gene qualifiers from gene features to overlapping features", \%opt_HH, \@opt_order_A);

# This section needs to be kept in sync (manually) with the opt_Add() section above
my %GetOptions_H = ();
my $options_okay = 
    &GetOptions('h'            => \$GetOptions_H{"-h"}, 
# basic options
                'f'            => \$GetOptions_H{"-f"},
                's=s'          => \$GetOptions_H{"-s"},
                'ftr'          => \$GetOptions_H{"--ftr"},
# options for controlling which features or qualifiers are output in attributes field
                'qall'         => \$GetOptions_H{"--qall"},
                'fskip=s'      => \$GetOptions_H{"--fskip"},
                'qskip=s'      => \$GetOptions_H{"--qskip"},
                'noaddgene'    => \$GetOptions_H{"--noaddgene"});

my $total_seconds = -1 * ofile_SecondsSinceEpoch(); # by multiplying by -1, we can just add another ofile_SecondsSinceEpoch call at end to get total time
my $execname_opt  = $GetOptions_H{"--execname"};
my $executable    = "annotate-tbl2gff.pl";
my $usage         = "Usage: $executable [-options]\n\t<path to v-annotate.pl output .tbl file>\n";
my $synopsis      = "$executable :: convert a v-annotate.pl .tbl or .ftr output file to GFF\n";
my $date          = scalar localtime();
my $version       = "1.7";
my $releasedate   = "Sep 2025";
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
my ($in_file) = (@ARGV);

# set options in opt_HH
opt_SetFromUserHash(\%GetOptions_H, \%opt_HH);

# validate options (check for conflicts)
opt_ValidateSet(\%opt_HH, \@opt_order_A);

# process --fskip and --qskip options
my %fskip_H    = (); # feature types to not output
my %qskip_H    = (); # qualifier values to not output in attributes field
if(opt_IsUsed("--fskip", \%opt_HH)) { utl_ExistsHFromCommaSepString(\%fskip_H, opt_Get("--fskip", \%opt_HH)); }
if(opt_IsUsed("--qskip", \%opt_HH)) { utl_ExistsHFromCommaSepString(\%qskip_H, opt_Get("--qskip", \%opt_HH)); }

#############################
# input the .tbl or .ftr file 
#############################
my %ftr_info_HAH = (); # hash of array of hashes with feature info 

my @seq_order_A = (); # array of sequence names read in order from feature table file
if(opt_Get("--ftr", \%opt_HH)) {
  ftr_parse($in_file, \%ftr_info_HAH, \@seq_order_A);
  feature_info_set_idx_and_parent_from_parent_idx_for_gff(\%ftr_info_HAH, \@seq_order_A);
}
else {
  local_feature_table_parse($in_file, \%ftr_info_HAH, \@seq_order_A, 0, undef);
  feature_info_set_id_and_parent_from_protein_id_for_gff(\%ftr_info_HAH, \@seq_order_A);
}
#utl_HAHDump("ftr_info_HAH", \%ftr_info_HAH, *STDOUT);

################
# output the GFF
################
gff_output(\%ftr_info_HAH, \@seq_order_A, *STDOUT, \%opt_HH, undef);

exit 0;

#################################################################
# Subroutine: local_feature_table_parse()
# Incept:     EPN, Thu Aug  8 10:41:39 2024
#             EPN, Mon May 20 09:36:09 2019
#
# Synopsis: Parse a INSDC feature table format file.
#           Very similar to sequips sqf_FeatureTableParse() but
#           with some modifications required to handle vadr output
#           .tbl files (e.g. 'Additional notes to submitter' lines)
#           and to keep track of sequence order.
#           We may eventually want to put this version into sequip.
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
sub local_feature_table_parse { 
  my $sub_name = "local_feature_table_parse";
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
# Subroutine: ftr_parse()
# Incept:     EPN, Wed Aug 14 11:38:06 2024
#             
#
# Synopsis: Parse a vadr ftr files into a %ftr_info_HAHR
#
# Arguments:
#  $infile:         feature table file to parse
#  $ftr_info_HAHR:  feature information, filled here
#                   1D key: accession
#                   2D:     feature index
#                   3D key: qualifer, value: qualifier value
#  $seq_order_AR:   REF to array of sequence names, added to here if defined, can be undef
#
# Returns:    void
#
# Dies:       if we have trouble parsing the file
#
# Reference: https://github.com/ncbi/vadr/blob/master/documentation/formats.md#explanation-of-ftr-suffixed-output-files
#################################################################
sub ftr_parse { 
  my $sub_name = "ftr_parse";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($infile, $ftr_info_HAHR, $seq_order_AR) = @_;

  open(IN, $infile) || ofile_FileOpenFailure($infile, $sub_name, $!, "reading", undef);

  my $line_ctr = 0;
  while(my $line = <IN>) { 
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    $line_ctr++;
    chomp $line;
    if($line =~ m/\r$/) { chop $line; } # remove ^M if it exists
    if(($line !~ m/^\#/) && ($line =~ m/\w/)) { 
      #      seq                               seq                  ftr          ftr                       ftr  ftr  par                                                                                                              seq                  model  ftr   
      #idx   name                              len  p/f   model     type         name                      len  idx  idx  str  n_from  n_to  n_instp  trc    5'N  3'N  p_from  p_to         p_instp  p_sc  nsa  nsn                coords                 coords  alerts
      #----  -------------------------------  ----  ----  --------  -----------  -----------------------  ----  ---  ---  ---  ------  ----  -------  -----  ---  ---  ------  ----  --------------  ----  ---  ---  --------------------  ---------------------  ------
      my @el_A = split(/\s+/, $line);
      my $ftr_idx;
      if(scalar(@el_A) == 26) { 
        my ($seqidx, $seq, $seqlen, $pf, $model, $ftr_type, $ftr_name, $ftr_len, $mdl_ftr_idx, $par_idx, $strand, $n_from, $n_to, $n_instp, $trc, $n5, $n3, $p_from, $p_to, $p_instp, $p_sc, $nsa, $nsn, $seq_coords, $mdl_coords, $alerts) = (@el_A);
        if($seq_coords ne "-") { # silently skip features with no coords, this can happen if a CDS has an indfantp alert
          if(! defined $ftr_info_HAHR->{$seq}) {
            @{$ftr_info_HAHR->{$seq}} = ();
            $ftr_idx = 0;
            push(@{$seq_order_AR}, $seq);
          }
          else {
            $ftr_idx = scalar(@{$ftr_info_HAHR->{$seq}});
          }
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"seq_len"}   = $seqlen;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"passfail"}  = $pf;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"model"}     = $model;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"type"}      = $ftr_type;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"name"}      = $ftr_name;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"ftr_len"}   = $ftr_len;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"mdl_ftr_idx"} = $mdl_ftr_idx;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"par_idx"}   = $par_idx;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"sumstrand"} = $strand;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"n_from"}    = $n_from;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"n_to"}      = $n_to;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"n_instp"}   = $n_instp;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"trc"}       = $trc;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"n5"}        = $n5;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"n3"}        = $n3;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"p_from"}    = $p_from;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"p_to"}      = $p_to;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"p_instp"}   = $p_instp;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"p_sc"}      = $p_sc;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"nsa"}       = $nsn;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"coords"}    = $seq_coords;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"mdlcoords"} = $mdl_coords;
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"alerts"}    = $alerts;
        }
      }
      else {
        die "ERROR, did not read expected number of tokens (26) on line $line_ctr of $infile\nline: $line\n"; 
      }
    }
  } 
  close(IN);
  
  return;
}

#################################################################
# Subroutine: feature_info_set_idx_and_parent_from_parent_idx_for_gff()
# Incept:     EPN, Wed Aug 14 14:42:58 2024
#
# Synopsis: Given feature info read from a .ftr file that includes
#           parent_idx qualifiers, further populate it with ID and
#           parent information in preparation for GFF output.
#
# Arguments:
#  $ftr_info_HAHR: feature info
#  $seq_order_AR:  ref to array with order of sequences to output, if undef, output in sorted order
#
# Returns:    void
#
# Dies:       never
#################################################################
sub feature_info_set_idx_and_parent_from_parent_idx_for_gff { 
  my $sub_name = "feature_info_set_idx_and_parent_from_parent_idx_for_gff";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_HAHR, $seq_order_AR) = @_;
  my @seq_A = ();
  if(! defined $seq_order_AR) {
    @seq_A = sort keys %{$ftr_info_HAHR};
  }
  else {
    @seq_A = @{$seq_order_AR};
  }
  
  # first pass through to determine ID value and
  # populate %protein_id2cds_id_H
  my %mdl_ftr_idx2cds_id_HH = (); # key1 is seq name, key 2 is model ftr_idx of a CDS (read from .ftr input file), value is that CDS' ID value
  my %ftr_id_idx_H = (); # key is feature name, value is index for this feature for current accn
  foreach my $seq (@seq_order_A) {
    if(! defined $ftr_info_HAHR->{$seq}) {
      ofile_FAIL("ERROR in $sub_name, no feature information for sequence $seq", 1, undef);
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
         (defined $ftr_info_HAHR->{$seq}[$ftr_idx]{"mdl_ftr_idx"})) {
        $mdl_ftr_idx2cds_id_HH{$seq}{$ftr_info_HAHR->{$seq}[$ftr_idx]{"mdl_ftr_idx"}} = $ftr_info_HAHR->{$seq}[$ftr_idx]{"GFF_ID"};
      }
    }
  }
  
  # second pass through to set parent values for mat_peptide and sig_peptide features
  foreach my $seq (@seq_order_A) {
    my $nftr = scalar(@{$ftr_info_HAHR->{$seq}});
    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      if(($ftr_info_HAHR->{$seq}[$ftr_idx]{"type"} eq "mat_peptide") ||
         ($ftr_info_HAHR->{$seq}[$ftr_idx]{"type"} eq "sig_peptide") && 
         (defined $ftr_info_HAHR->{$seq}[$ftr_idx]{"par_idx"})) {
        my $parent_idx = $ftr_info_HAHR->{$seq}[$ftr_idx]{"par_idx"};
        if(defined $mdl_ftr_idx2cds_id_HH{$seq}{$parent_idx}) {
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"GFF_Parent"} = $mdl_ftr_idx2cds_id_HH{$seq}{$parent_idx};
        }
      }
    }
  }    
  return;
}

#################################################################
# Subroutine: feature_info_set_id_and_parent_from_protein_id_for_gff()
# Incept:     EPN, Thu Aug  8 14:07:34 2024
#
# Synopsis: Given feature info read from a .tbl file that includes
#           protein_id qualifiers, further populate it with ID and
#           parent information in preparation for GFF output.
#
# Arguments:
#  $ftr_info_HAHR: feature info
#  $seq_order_AR:  ref to array with order of sequences to output, if undef, output in sorted order
#
# Returns:    void
#
# Dies:       never
#################################################################

sub feature_info_set_id_and_parent_from_protein_id_for_gff { 
  my $sub_name = "feature_info_set_id_and_parent_from_protein_id_for_gff";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_HAHR, $seq_order_AR) = @_;
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
      ofile_FAIL("ERROR in $sub_name, no feature information for sequence $seq", 1, undef);
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
        if((defined $protein_id) && (defined $protein_id2cds_id_H{$protein_id})) {
          $ftr_info_HAHR->{$seq}[$ftr_idx]{"GFF_Parent"} = $protein_id2cds_id_H{$protein_id};
        }
      }
    }
  }    
  return;
}

#################################################################
# Subroutine: gff_output()
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
#
# Returns:    void
#
# Dies:       never
#################################################################
sub gff_output { 
  my $sub_name = "gff_output";
  my $nargs_expected = 5;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  
  my ($ftr_info_HAHR, $seq_order_AR, $out_FH, $opt_HHR) = @_;

  my %fskip_H    = (); # feature types to not output
  my %qskip_H    = (); # qualifier values to not output in attributes field
  if(opt_IsUsed("--fskip", $opt_HHR)) { utl_ExistsHFromCommaSepString(\%fskip_H, opt_Get("--fskip", $opt_HHR)); }
  if(opt_IsUsed("--qskip", $opt_HHR)) { utl_ExistsHFromCommaSepString(\%qskip_H, opt_Get("--qskip", $opt_HHR)); }
  my $do_ftr  = (opt_IsUsed("--ftr",   $opt_HHR)) ? 1 : 0;
  my $do_qall = (opt_IsUsed("--qall",  $opt_HHR)) ? 1 : 0;

  my @seq_A = ();
  if(! defined $seq_order_AR) {
    @seq_A = sort keys %{$ftr_info_HAHR};
  }
  else {
    @seq_A = @{$seq_order_AR};
  }
  
  my $source = (opt_IsUsed("-s", $opt_HHR)) ? opt_Get("-s", $opt_HHR) : "VADR:v" . $version;
  foreach my $seq (@seq_order_A) {
    if(! defined $ftr_info_HAHR->{$seq}) {
      ofile_FAIL("ERROR in $sub_name, no feature information for sequence $seq", 1, undef);
    }
    my $nftr = scalar(@{$ftr_info_HAHR->{$seq}});

    # Note: we do not impute gene qualifiers, based on overlap, because this can fail due to ambiguity
    # (example is flu seq CY089961.1 which has a CDS that overlaps with two genes (M1 and M2, model:CY002009).
    # v-annotate.pl will put gene qualifiers in the output tbl if --forcegene is used.

    for(my $ftr_idx = 0; $ftr_idx < $nftr; $ftr_idx++) { 
      my $ftr_type = $ftr_info_HAHR->{$seq}[$ftr_idx]->{"type"};   
      if(! defined $fskip_H{$ftr_type}) { 
        my @start_A  = ();
        my @stop_A   = ();
        my @strand_A = ();
        if(! defined $ftr_info_HAHR->{$seq}[$ftr_idx]{"coords"}) {
          die "ERROR coords undefined for seq $seq ftr_idx: $ftr_idx\n";
        }
        vdr_FeatureStartStopStrandArrays($ftr_info_HAHR->{$seq}[$ftr_idx]{"coords"}, \@start_A, \@stop_A, \@strand_A, undef);
        my $nsgm = scalar(@start_A);
        my @ftr_lines_A = (); # we'll store the feature lines here, and output them in forward order or reverse order depending on strand
        my $summary_strand = vdr_FeatureSummaryStrand($ftr_info_HAHR->{$seq}[$ftr_idx]{"coords"}, undef);
        my $id = $ftr_info_HAHR->{$seq}[$ftr_idx]{"GFF_ID"};
        if(! defined $id) {
          ofile_FAIL("ERROR in $sub_name, no ID set for feature for sequence $seq", 1, undef);
        }
        $id = "ID:" . $id;
        if(defined $qskip_H{"ID"}) {
          $id = "";
        }
        my $key_values = "";
        foreach my $key (sort keys(%{$ftr_info_HAHR->{$seq}[$ftr_idx]})) { 
          if(($key ne "GFF_ID") &&
             (($key ne "coords") || $do_qall) &&
             (($key ne "type")   || $do_qall) &&
             (($key ne "protein_id") || $do_qall)) { 
            my $key2print = $key;
            $key2print =~ s/^GFF\_//;
            if(! defined $qskip_H{$key2print}) { 
              $key_values .= $key2print . "=" . $ftr_info_HAHR->{$seq}[$ftr_idx]{$key} . ";";
            }
          }
        }
        my $codon_start = undef;
        if(($ftr_type eq "CDS") && (! $do_ftr)) {
          # --ftr NOT used, so we should have codon_start info for the first segment
          # if it is ! defined then we can assume it is 1. We need to determine
          # the 'phase' value for each segment
          # if codon_start == 1, then phase = 0
          # if codon_start == 2, then phase = 1
          # if codon_start == 3, then phase = 2
          $codon_start = (defined $ftr_info_HAHR->{$seq}[$ftr_idx]{"codon_start"}) ?
              $ftr_info_HAHR->{$seq}[$ftr_idx]{"codon_start"} : 1;
        }
          
        my $cum_length = 0;
        for(my $sgm_idx = 0; $sgm_idx < $nsgm; $sgm_idx++) {
          my ($start, $stop, $strand) = (undef, undef, undef);
          if($strand_A[$sgm_idx] eq "+") {
            ($start, $stop, $strand) = ($start_A[$sgm_idx], $stop_A[$sgm_idx], $strand_A[$sgm_idx]);
          }
          else {
            ($stop, $start, $strand) = ($start_A[$sgm_idx], $stop_A[$sgm_idx], $strand_A[$sgm_idx]);
          }
          my $attributes = $id . ";" . $key_values;
          my $phase = ".";
          if(defined $codon_start) {
            $phase = vdr_FrameAdjust($codon_start, $cum_length, undef) - 1;
          }
          push(@ftr_lines_A,
               sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
                       $seq,                                 # token 1: 'sequence' (sequence name)
                       $source,                              # token 2: 'source' (VADR:v1.6.4)
                       $ftr_type,                            # token 3: 'feature'
                       $start,                               # token 4: 'start' in coordinate space [1..seqlen], must be <= 'end'
                       $stop,                                # token 5: 'end' in coordinate space [1..seqlen], must be >= 'start'
                       ".",                                  # token 6: 'score' 
                       $strand,                              # token 7: 'strand' ('+' or '-')
                       $phase,                               # token 8: 'phase'
                       $attributes));                        # token 9: attributes, remaining qualifiers read from .tbl/.ftr input file
          $cum_length += vdr_CoordsLength(vdr_CoordsSegmentCreate($start_A[$sgm_idx], $stop_A[$sgm_idx], $strand_A[$sgm_idx], undef), undef);
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
      } # end of 'if(! defined $fskip_H{$ftr_type})'
    }
  }
  
  return;
}


