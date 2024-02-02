#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my ($exactMatch, $blast, $ResulstOut);
&GetOptions("exactMatchFiles=s"=> \$exactMatch,
            "blastOutput=s"=> \$blast,
            "ResulstOut=s"=> \$ResulstOut,
           ) ;
die("Please provide both input and out put") unless ($exactMatch & $blast);

my $peptideMatchFile = $exactMatch;

my $blastFile = $blast;

my $outFile = $ResulstOut;

sub processEpitopeps{
    open(my $pepHandle, $peptideMatchFile) or die "Could not open file '$peptideMatchFile' $!";
    open(FH ,'>' ,$outFile) or die $!;

    my %peptideHash = ();
    while (my $row = <$pepHandle>) {
        chomp $row;
        my @counts_list = split("\t", $row);
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

    open(my $blastHandle, $blastFile) or die "Could not open file '$blastFile' $!";

    while (my $row = <$blastHandle>) {
        chomp $row;
        my @counts_list = split("\t", $row);
        my $proteinID = $counts_list[0];
        my $peptideID = $counts_list[1];
        my $pepLen = $counts_list[2];
        my $bitScore = $counts_list[3];
        my $matchStart = $counts_list[4];
        my $matchEnd = $counts_list[5];
        my $indetity = $counts_list[6];
        my $alignmentLength = $counts_list[7];
        my $refSeq = $counts_list[8];
        my $hitSeq = $counts_list[9];
        my $alignment = $counts_list[10];
       
        
        my $key = $proteinID . "_" . $peptideID;

        if (exists($peptideHash{$key})){
            $peptideHash{$key}{pepLen} = $pepLen;
            $peptideHash{$key}{bitScore} = $bitScore;
            $peptideHash{$key}{indetity} = $indetity;
            $peptideHash{$key}{alignmentLength} = $alignmentLength;
            $peptideHash{$key}{refSeq} = $refSeq;
            $peptideHash{$key}{hitSeq} = $hitSeq;
        } else {
            $peptideHash{$key}{protein} = $proteinID;
            $peptideHash{$key}{peptide} = $peptideID;
            $peptideHash{$key}{pepLen} = $pepLen;
            $peptideHash{$key}{bitScore} = $bitScore;
            $peptideHash{$key}{indetity} = $indetity;
            $peptideHash{$key}{matchStart} = $matchStart;
            $peptideHash{$key}{matchEnd} = $matchEnd;
            $peptideHash{$key}{alignmentLength} = $alignmentLength;
            $peptideHash{$key}{refSeq} = $refSeq;
            $peptideHash{$key}{hitSeq} = $hitSeq;
            $peptideHash{$key}{alignment} = $alignment;

        }

    }

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
        if (exists($peptideHash{$key}{pepLen})){
            push @currentList, $peptideHash{$key}{pepLen};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{bitScore})){
            push @currentList, $peptideHash{$key}{bitScore};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{indetity})){
            push @currentList, $peptideHash{$key}{indetity};
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
        if (exists($peptideHash{$key}{alignmentLength})){
            push @currentList, $peptideHash{$key}{alignmentLength};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{refSeq})){
            push @currentList, $peptideHash{$key}{refSeq};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{hitSeq})){
            push @currentList, $peptideHash{$key}{hitSeq};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{alignment})){
            push @currentList, $peptideHash{$key}{alignment};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{pepMatch})){
            push @currentList, $peptideHash{$key}{pepMatch};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{proteinMatch})){
            push @currentList, $peptideHash{$key}{proteinMatch};
        } else {
            push @currentList, " ";
        }
        if (exists($peptideHash{$key}{taxonMatch})){
            push @currentList, $peptideHash{$key}{taxonMatch};
        } else {
            push @currentList, " ";
        }
        

        print FH ($currentList[0], "\t", $currentList[1],  "\t", $currentList[2], "\t", $currentList[3],  "\t", $currentList[4], "\t", $currentList[5], "\t", $currentList[6], "\t", 
        , $currentList[7], "\t", $currentList[8],"\t" , $currentList[9], "\t", $currentList[10], "\t", $currentList[11],  "\t", $currentList[12], "\t", $currentList[13], "\n")
    }
}

processEpitopeps()