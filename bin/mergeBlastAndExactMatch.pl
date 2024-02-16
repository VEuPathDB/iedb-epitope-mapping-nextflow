#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

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


        #TODO:  what is in col 6?
        my $proteinID = $counts_list[0];
        my $peptideID = $counts_list[5];
        my $pepMatch = $counts_list[1];
        my $proteinMatch = $counts_list[2];
        my $taxonMatch = $counts_list[3];
        my $alignment = $counts_list[4];
        my $matchStart = $counts_list[7];
        my $matchEnd = $counts_list[8];

        my $key = $proteinID . "_" . $peptideID;

        $peptideHash{$key}{protein} = $proteinID;
        $peptideHash{$key}{peptide} = $peptideID;
        $peptideHash{$key}{pepMatch} = $pepMatch;
        $peptideHash{$key}{proteinMatch} = $proteinMatch;
        $peptideHash{$key}{taxonMatch} = $taxonMatch;
        $peptideHash{$key}{matchStart} = $matchStart;
        $peptideHash{$key}{matchEnd} = $matchEnd;
        $peptideHash{$key}{alignment} = $alignment;
    }

    close $pepHandle;

    my @colNames = ("pepMatch",
                    "proteinMatch",
                    "taxonMatch",
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


    print OUT $proteinID . "\t" . $peptideID . "\t" . join("\t", @peptideExactMatchValues) . "\t" . join("\t", @a) . "\n";

    $peptideHashRef->{$key}->{foundBlast} = 1;
}
close $blastHandle;

foreach my $pepProteinKey (keys %{$peptideHashRef}) {
    next if($peptideHashRef->{$pepProteinKey}->{foundBlast});

    my $peptideId = $peptideHashRef->{$pepProteinKey}->{peptide};
    my $proteinId = $peptideHashRef->{$pepProteinKey}->{protein};

    my @peptideExactMatchValues = map { $peptideHashRef->{$pepProteinKey}->{$_} } @$colNames;

    # we are not writing out empty columns; will this cause problems downstream?
    print OUT $proteinId . "\t" . $peptideId . "\t" . join("\t", @peptideExactMatchValues) . "\n";
}


1;
