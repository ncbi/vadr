#!/usr/bin/env perl

# Note: This script was generated with the assistance of GitHub Copilot using the Gemini 3.1 Pro (Preview) model.

use strict;
use warnings;
use Getopt::Long qw(:config no_auto_abbrev);
use LWP::Simple;
use HTML::Entities qw(decode_entities);

my $taxid     = "";
my $outprefix = "";
my $api_key   = "";
my $maxdate   = "";
my $mindate   = "";
my $quiet     = 0;

GetOptions(
    "taxid=s"   => \$taxid,
    "out=s"     => \$outprefix,
    "api_key=s" => \$api_key,
    "maxdate=s" => \$maxdate,
    "mindate=s" => \$mindate,
    "q|quiet"   => \$quiet,
);

if (!$taxid || !$outprefix) {
    print STDERR "Usage: $0 --taxid <taxid> --out <output_prefix> [--api_key <key>] [--maxdate <YYYY/MM/DD>] [--mindate <YYYY/MM/DD>]\n";
    print STDERR "Example: $0 --taxid 138951 --out evD_metadata --api_key ABC123DEF456\n";
    exit(1);
}

if ($maxdate && $maxdate !~ /^\d{4}\/\d{2}\/\d{2}$/) {
    die "ERROR: --maxdate format invalid; expected YYYY/MM/DD, got $maxdate\n";
}
if ($mindate && $mindate !~ /^\d{4}\/\d{2}\/\d{2}$/) {
    die "ERROR: --mindate format invalid; expected YYYY/MM/DD, got $mindate\n";
}

my $base_esearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
my $base_esummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi";

# Avoid patent sequences, synthetic constructs, environmental samples, and
# UNVERIFIED GenBank submissions (whose DEFINITION lines are prefixed with
# "UNVERIFIED:"). UNVERIFIED entries contain inferred annotation rather than
# curator-verified data and should not be used as training/test data.
my $query = "txid$taxid\[Organism\] NOT (gbdiv_syn\[prop\] OR gbdiv_pat\[prop\] OR env_sample\[prop\] OR UNVERIFIED\[Title\])";

my $api_params = "";
if ($api_key) { $api_params .= "&api_key=$api_key"; }

if ($maxdate || $mindate) {
    my $min_part = $mindate || "0001/01/01";
    my $max_part = $maxdate || "3000/12/31";
    $query .= " AND $min_part:$max_part\[pdat\]";
}

# Step 1: esearch with usehistory to get WebEnv
my $esearch_url = "$base_esearch?db=nuccore&term=" . _url_encode($query) . "&usehistory=y&retmax=0" . $api_params;
print "- Querying NCBI esearch for taxonomy ID $taxid...\n" unless $quiet;
if (!$quiet && ($mindate || $maxdate)) {
    my @bounds = ();
    push @bounds, "mindate=$mindate" if $mindate;
    push @bounds, "maxdate=$maxdate" if $maxdate;
    print "- Date bounds: " . join(" ", @bounds) . "\n";
}

my $esearch_res = get($esearch_url);
if (!defined($esearch_res)) {
    die "ERROR: Failed to fetch from NCBI esearch. It may be temporarily down or blocking connections.\n";
}

my $count = 0;
my $query_key = "";
my $webenv = "";

if ($esearch_res =~ /<Count>(\d+)<\/Count>/) {
    $count = $1;
}
if ($esearch_res =~ /<QueryKey>([^<]+)<\/QueryKey>/) {
    $query_key = $1;
}
if ($esearch_res =~ /<WebEnv>([^<]+)<\/WebEnv>/) {
    $webenv = $1;
}

if ($count == 0 || !$query_key || !$webenv) {
    print STDERR "No sequences found for taxid $taxid or failed to parse esearch XML.\n";
    exit(0);
}

print "  Found $count candidate sequences. Preparing metadata fetch...\n" unless $quiet;

# Step 2: Fetch summaries in batches of 500
my $retstart = 0;
my $retmax   = 500;
my $out_tsv  = "$outprefix.tsv";

open(my $out_fh, ">", $out_tsv) or die "ERROR: Cannot open $out_tsv for writing: $!";
print $out_fh "Accession\tLength\tCreateDate\tSerotype\tGenotype\tIsolate\tTitle\n";

# If user provided an API key, we can afford a higher rate limit. Without one, sleep more.
my $sleep_time = ($api_key) ? 0.35 : 1.0; 

while ($retstart < $count) {
    my $batch_url = "$base_esummary?db=nuccore&query_key=$query_key&WebEnv=$webenv&retstart=$retstart&retmax=$retmax&version=2.0" . $api_params;
    print "  Fetching records $retstart - " . ($retstart + $retmax - 1 > $count - 1 ? $count - 1 : $retstart + $retmax - 1) . "\n" unless $quiet;
    
    my $esummary_res = get($batch_url);
    if (!defined($esummary_res)) {
         print STDERR "WARNING: Failed to fetch esummary batch at retstart=$retstart. Will retry once...\n";
         sleep(3);
         $esummary_res = get($batch_url);
         if(!defined($esummary_res)) {
             die "ERROR: Failed twice to fetch batch at retstart=$retstart. Aborting.\n";
         }
    }

    # Split into DocumentSummary blocks
    my @docs = split(/<DocumentSummary\s+[^>]+>/, $esummary_res);
    shift @docs; # First element is the header before the first document

    foreach my $doc (@docs) {
        my $acc      = ($doc =~ /<AccessionVersion>([^<]+)<\/AccessionVersion>/) ? $1 : "";
        my $slen     = ($doc =~ /<Slen>(\d+)<\/Slen>/) ? $1 : "";
        my $cdate    = ($doc =~ /<CreateDate>([^<]+)<\/CreateDate>/) ? $1 : "";
        my $subtype  = ($doc =~ /<SubType>([^<]+)<\/SubType>/) ? $1 : "";
        my $subname  = ($doc =~ /<SubName>([^<]+)<\/SubName>/) ? $1 : "";
        my $title    = ($doc =~ /<Title>([^<]*)<\/Title>/) ? $1 : "";

        # Clean HTML entities if any
        $subtype = defined(&decode_entities) ? decode_entities($subtype) : $subtype;
        $subname = defined(&decode_entities) ? decode_entities($subname) : $subname;
        $title   = defined(&decode_entities) ? decode_entities($title)   : $title;
        # Strip TSV-breaking whitespace from title
        $title =~ s/[\t\r\n]+/ /g;

        my $serotype = "";
        my $genotype = "";
        my $isolate  = "";

        if ($subtype && $subname) {
            my @types = split(/\|/, $subtype);
            my @names = split(/\|/, $subname);
            
            for (my $i = 0; $i < scalar(@types); $i++) {
                if (defined($types[$i]) && defined($names[$i])) {
                    if (lc($types[$i]) eq 'serotype') {
                        $serotype = $names[$i];
                    } elsif (lc($types[$i]) eq 'genotype') {
                        $genotype = $names[$i];
                    } elsif (lc($types[$i]) eq 'isolate') {
                        $isolate = $names[$i];
                    }
                }
            }
        }

        if ($acc) {
            print $out_fh "$acc\t$slen\t$cdate\t$serotype\t$genotype\t$isolate\t$title\n";
        }
    }

    $retstart += $retmax;
    sleep($sleep_time); # polite to NCBI rate limits
}

close($out_fh);
print "- Finished successfully! Metadata saved to: $out_tsv\n" unless $quiet;

exit(0);

sub _url_encode {
    my $str = shift;
    $str =~ s/([^A-Za-z0-9])/sprintf("%%%02X", ord($1))/seg;
    return $str;
}
