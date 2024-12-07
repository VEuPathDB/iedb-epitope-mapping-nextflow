#!/usr/bin/env perl

use strict;

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


# NOTE: when joining these 2 files, by proteinSourceId + epitopeId we can first
# read the blast results into memory.  this is an "outer join" operation.
# take everything from the exact match file and append the blast data where available

my ($blastHashRef) = &readBlastsAsHash($blastFile);

&readExactMatchesAndPrint($peptideMatchFile, $blastHashRef, $outFile);



sub usage {
    my ($e) = @_;

    print STDERR "mergeBlastAndExactMatch.pl --exactMatchFile <FILE> --blastFile <FILE> --outputFile OUT\n";
    die $e if($e);
}

sub readBlastsAsHash {
    my ($blastFile) = @_;

    my @columns = ("proteinDefline", "epitopeId", "peptideLength", "bitScore", "evalue", "hitStart", "hitEnd", "hitIdentity", "alignmentLength", "querySequence", "hitSequence", "alignment");

    open(my $blastHandle, $blastFile) or die "Could not open file '$blastFile' $!";

    my %blastHash;

    while (my $row = <$blastHandle>) {
        chomp $row;
        next unless $row;

        my @v = split("\t", $row);

        my %hash;
        @hash{@columns} = @v;

        my $proteinDefline = $hash{proteinDefline};

        if(!$proteinDefline) {
                    print Dumper \%hash;
        exit;

        }
        my ($proteinId) = $proteinDefline =~ /([^\s]+)/;

        my $epitopeId = $hash{epitopeId};

        my $key = $proteinId . "_" . $epitopeId;

        $blastHash{$key} = \%hash;
    }
    close $blastHandle;

    return \%blastHash;
}


sub readExactMatchesAndPrint {
    my ($peptideMatchFile, $blastHashRef, $outFile) = @_;

    open(OUT, '>', $outFile) or die "Cannot open file $outFile for writing:  $!";

    open(my $pepHandle, $peptideMatchFile) or die "Could not open file '$peptideMatchFile' $!";

    my @columns = ('proteinId', 'peptideMatch', 'proteinMatch', 'taxonMatch', 'peptideLength', 'peptide', 'epitopeId', 'epitopeAccession', 'matchStart', 'matchEnd', 'epitopeTaxon');

    while (my $row = <$pepHandle>) {
        chomp $row;
        next unless $row;

        my @v = split("\t", $row);

        my %hash;
        @hash{@columns} = @v;

        my $proteinId = $hash{proteinId};
        my $epitopeId = $hash{epitopeId};

        my $key = $proteinId . "_" . $epitopeId;

        my ($hitStart, $hitEnd, $hitIdentity, $alignmentLength, $bitScore, $evalue, $alignment);
        if($blastHashRef->{$key}) {
            $hitStart = $blastHashRef->{$key}->{hitStart};
            $hitEnd = $blastHashRef->{$key}->{hitEnd};
            $hitIdentity = $blastHashRef->{$key}->{hitIdentity};
            $alignmentLength = $blastHashRef->{$key}->{alignmentLength};
            $bitScore = $blastHashRef->{$key}->{bitScore};
            $evalue = $blastHashRef->{$key}->{evalue};
            $alignment = $blastHashRef->{$key}->{alignment};
        }

        print OUT $hash{proteinId}
            , "\t", $hash{epitopeId}
            , "\t", $hash{peptideMatch}
            , "\t", $hash{proteinMatch}
            , "\t", $hash{taxonMatch}
            , "\t", $hash{matchStart}
            , "\t", $hash{matchEnd}
            , "\t", $hitStart
            , "\t", $hitEnd
            , "\t", $hitIdentity
            , "\t", $alignmentLength
            , "\t", $bitScore
            , "\t", $evalue
            , "\t", $alignment
            , "\n";
    }

    close $pepHandle;
    close OUT;

}


1;
