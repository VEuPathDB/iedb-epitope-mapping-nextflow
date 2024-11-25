#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;


# This code macthes the outputs from the blast and the exact matches into a single text file. 
my ($peptideMatchFile, $blastFile, $outFile);
&GetOptions("exactMatchFile=s"=> \$peptideMatchFile,
            "blastFile=s"=> \$blastFile,
            "outputFile=s"=> \$outFile,
    ) ;

unless (-e $peptideMatchFile && -e $blastFile && $outFile) {
    &usage("Both input files must exist;  output file must be declared")
}

sub usage {
    my ($e) = @_;

    print STDERR "mergeBlastAndExactMatch.pl --exactMatchFile <FILE> --blastFile <FILE> --outputFile OUT\n";
    die $e if($e);
}

sub readExactMatchesAsHash {
    my ($peptideMatchFile) = @_;

    open(my $pepHandle, $peptideMatchFile) or die "Could not open file '$peptideMatchFile' $!";

    my %peptideHash;
    while (my $row = <$pepHandle>) {
        chomp $row;
        my @counts_list = split("\t", $row);

        my $proteinID = $counts_list[0];
        my $peptideID = $counts_list[6];
        my $pepMatch = $counts_list[1];
        my $proteinMatch = $counts_list[2];
        my $taxonMatch = $counts_list[3];
        my $pepLen = $counts_list[4];
        my $alignment = $counts_list[5];
        my $matchStart = $counts_list[8];
        my $matchEnd = $counts_list[9];

        my $key = $proteinID . "_" . $peptideID;

        $peptideHash{$key}{protein} = $proteinID;
        $peptideHash{$key}{peptide} = $peptideID;
        $peptideHash{$key}{pepMatch} = $pepMatch;
        $peptideHash{$key}{proteinMatch} = $proteinMatch;
        $peptideHash{$key}{taxonMatch} = $taxonMatch;
        $peptideHash{$key}{matchStart} = $matchStart;
        $peptideHash{$key}{matchEnd} = $matchEnd;
        $peptideHash{$key}{pepLen} = $pepLen;
        $peptideHash{$key}{alignment} = $alignment;
    }

    close $pepHandle;

    my @colNames = ("pepMatch",
                    "proteinMatch",
                    "taxonMatch",
                    "pepLen",
                    "matchStart",
                    "matchEnd",
                    "alignment");

    return \%peptideHash, \@colNames;
}

my ($peptideHashRef, $colNames) = &readExactMatchesAsHash($peptideMatchFile);

open(my $blastHandle, $blastFile) or die "Could not open file '$blastFile' $!";
open(OUT, '>', $outFile) or die $!;

while (my $row = <$blastHandle>) {
    chomp $row;
    my @a = split("\t", $row);
    my $proteinID = shift @a;
    my $peptideID = shift @a;
    my $key = $proteinID . "_" . $peptideID;

    my @peptideExactMatchValues = map { $peptideHashRef->{$key}->{$_} } @$colNames;

    if ($peptideID) {
        print OUT $proteinID . "\t" . $peptideID . "\t" . join("\t", @peptideExactMatchValues[0..2]) . "\t". join("\t", @a) . "\n";
    }

    $peptideHashRef->{$key}->{foundBlast} = 1;
}

foreach my $pepProteinKey (keys %{$peptideHashRef}) {
    next if($peptideHashRef->{$pepProteinKey}->{foundBlast});

    my $peptideID = $peptideHashRef->{$pepProteinKey}->{peptide};
    my $proteinId = $peptideHashRef->{$pepProteinKey}->{protein};

    my @peptideExactMatchValues = map { $peptideHashRef->{$pepProteinKey}->{$_} } @$colNames;

    my $start =  @peptideExactMatchValues[4];
    my $end =  @peptideExactMatchValues[5];
    my $seq = @peptideExactMatchValues[6];
    print OUT $proteinId . "\t" . $peptideID . "\t" . join("\t", @peptideExactMatchValues[0..3]) . "\t". "\t". "\t". $start . "\t" . $end . "\t". "\t". "\t". $seq . "\t". $seq . "\t" .$seq . "\n";
}
close $blastHandle;

1;
