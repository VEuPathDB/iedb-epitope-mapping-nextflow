#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my ($exactMatch, $blast);
&GetOptions("exactMatchFiles=s"=> \$exactMatch,
            "blastOutput=s"=> \$blast,
           ) ;
die("Please provide both input and out put") unless ($exactMatch & $blast);

my $peptideMatchFile = $exactMatch; # $ARGV[0]; #"/Users/saikouybah/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/IedbAnalysis/Result/PeptideGene.txt";

my $blastFile = $blast; # $ARGV[1]; #"/Users/saikouybah/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/IedbAnalysis/Result/BlastOut/AnnotatedProteins.txt";

my $outFile = "./EpitopesSearchResults.txt";

sub loadEpitopeps{
    open(my $pepHandle, $peptideMatchFile) or die "Could not open file '$peptideMatchFile' $!";
    open(FH ,'>' ,$outFile) or die $!;

    my %peptideHash = ();
    while (my $row = <$pepHandle>) {
        chomp $row;
        my @counts_list = split("\t", $row); #split /\s+/,$row;
        my $proteinID = $counts_list[0];
        my $peptideID = $counts_list[3];
        my $matchType = $counts_list[1];
        my $matchStart = $counts_list[5];
        my $matchEnd = $counts_list[6];
        my $key = $proteinID . "_" . $peptideID;
        $peptideHash{$key}{protein} = $proteinID;
        $peptideHash{$key}{peptide} = $peptideID;
        $peptideHash{$key}{MatchType} = $matchType;
        $peptideHash{$key}{matchStart} = $matchStart;
        $peptideHash{$key}{matchEnd} = $matchEnd;

    }

    open(my $blastHandle, $blastFile) or die "Could not open file '$blastFile' $!";

    while (my $row = <$blastHandle>) {
        chomp $row;
        my @counts_list = split("\t", $row);
        my $proteinID = $counts_list[0];
        my $peptideID = $counts_list[1];
        my $PercentageIdentity = $counts_list[2];
        my $alignmentLength = $counts_list[3];
        my $matchStart = $counts_list[6];
        my $matchEnd = $counts_list[7];
        
        my $key = $proteinID . "_" . $peptideID;

        if (exists($peptideHash{$key})){
            $peptideHash{$key}{blastIdentity} = $PercentageIdentity;
            $peptideHash{$key}{alignmentLenght} = $alignmentLength;
        } else {
            $peptideHash{$key}{protein} = $proteinID;
            $peptideHash{$key}{peptide} = $peptideID;
            $peptideHash{$key}{blastIdentity} = $PercentageIdentity;
            $peptideHash{$key}{alignmentLenght} = $alignmentLength;
            $peptideHash{$key}{matchStart} = $matchStart;
            $peptideHash{$key}{matchEnd} = $matchEnd;

        }

    }
    # print Dumper {%peptideHash}

    foreach my $key (keys %peptideHash) {
        my @currentList = ();

        if (exists($peptideHash{$key}{protein})){
            push @currentList, $peptideHash{$key}{protein};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{peptide})){
            push @currentList, $peptideHash{$key}{peptide};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{MatchType})){
            push @currentList, $peptideHash{$key}{MatchType};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{matchStart})){
            push @currentList, $peptideHash{$key}{matchStart};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{matchEnd})){
            push @currentList, $peptideHash{$key}{matchEnd};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{blastIdentity})){
            push @currentList, $peptideHash{$key}{blastIdentity};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{alignmentLenght})){
            push @currentList, $peptideHash{$key}{alignmentLenght};
        } else {
            push @currentList, " ";
        }

        print FH ($currentList[0], "\t", $currentList[1],  "\t", $currentList[3], "\t", $currentList[4], "\t", $currentList[5], "\t", $currentList[6], "\t", $currentList[2],"\n")
    }
}

loadEpitopeps()