#!/usr/bin/env perl

use strict;

use Data::Dumper;
use Getopt::Long;


# This code macthes the outputs from the blast and the exact matches into a single text file. 
my ($pepMatchResults, $fullSeqMatchResults, $peptideTab, $distinctPeptidesFile, $outputFile);
&GetOptions("pepMatchResults=s"=> \$pepMatchResults,
            "fullSeqMatchResults=s" => \$fullSeqMatchResults,
            "peptideTab=s" => \$peptideTab,
            "outputFile=s"=> \$outputFile,
            "distinctPeptides=s" => \$distinctPeptidesFile,
    ) ;


my %distinctPeptides;
open(PEPTIDES, $distinctPeptidesFile) or die "Cannot open file $distinctPeptidesFile for reading: $!";
while(<PEPTIDES>) {
    chomp;
    $distinctPeptides{$_} = 1;
}
close PEPTIDES;

my %peptideToIedb;
open(PEPTAB, $peptideTab) or die "Cannot open file $peptideTab for reading: $!";
while(<PEPTAB>) {
    chomp;
    my @a = split(/\t/, $_);
    my $peptide = $a[3];
    next unless($distinctPeptides{$peptide});

    my $iedbId = $a[1];
    $peptideToIedb{$peptide} = $iedbId;
}
close PEPTAB;


my %fullSeqMatchResults;
open(FULLRES, $fullSeqMatchResults) or die "Cannot open file $fullSeqMatchResults for reading: $!";
while(<FULLRES>) {
    chomp;
    my ($peptide, $protein, $iedb, $taxonMatch, $proteinMatch) = split(",", $_);

    my $key = "${peptide}_${protein}";

    $fullSeqMatchResults{$key}->{taxon_match} = $taxonMatch;
    $fullSeqMatchResults{$key}->{full_protein_match} = $proteinMatch;

}
close FULLRES;

open(OUT, ">$outputFile") or die "Cannot open output file $outputFile for writing: $!";
open(PEPMATCH, $pepMatchResults) or die "Cannot open file $pepMatchResults for reading: $!";


while(<PEPMATCH>) {
    chomp;
    my ($peptide, $protein, $mismatches, $start, $end) = split(",", $_);

    # pepMatch preprocess db adds a ".1 suffix to proteins; remove it
    $protein =~ s/\.\d$//;

    my $iedbId = $peptideToIedb{$peptide};
    die "Could not map peptide to iedb id " unless $iedbId;

    my $key = "${peptide}_${protein}";

    my $taxonMatch = 0;
    my $fullProteinMatch = 0;

    if($fullSeqMatchResults{$key}) {
        $taxonMatch = 1 if $fullSeqMatchResults{$key}->{taxon_match};
        $fullProteinMatch = 1 if $fullSeqMatchResults{$key}->{full_protein_match};
    }

    print OUT join("\t", $peptide, $protein, $iedbId, $taxonMatch, $fullProteinMatch, $mismatches, $start, $end, "\n");
}


close PEPMATCH;
close OUT;


1;
