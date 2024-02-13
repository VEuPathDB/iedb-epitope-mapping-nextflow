#!/usr/bin/env perl

use strict;
use XML::LibXML;
use Data::Dumper;
use open qw(:std :utf8);
use Getopt::Long;

my ($xmlFile);
my ($outputFile);
&GetOptions("inputXml=s"=> \$xmlFile,
            "outputFile=s"=> \$outputFile,
           ) ;
die("Please provide both input and out put") unless ($xmlFile & $outputFile);


my $xml = XML::LibXML->load_xml(location => $xmlFile);

my @queries = $xml->findnodes('/BlastOutput/BlastOutput_iterations/Iteration');


my $outfile = $outputFile;

#FIXME:  should use ">" NOT ">>"
open(FH, '>>', $outfile) or die $!;

foreach my $query (@queries) {

   my $refName = ($query->findvalue('./Iteration_query-def')); 
   my @refNameSplit = split / /, $refName;
   my $pepID = $refNameSplit[0];

   my @peptideID = ($query->findnodes('./Iteration_hits/Hit')); 
   foreach my $hit (@peptideID) {
        my $name = $hit->findvalue('./Hit_def');
        my $pepLen = $hit->findvalue('./Hit_len');
        my @Hit_hsps = ($hit->findnodes('./Hit_hsps/Hsp'));
        
        foreach my $hsp (@Hit_hsps) {
            my $bitScore = $hsp->findvalue('./Hsp_bit-score');
            my $evalue = $hsp->findvalue('./Hsp_evalue');
            my $hitStart = $hsp->findvalue('./Hsp_hit-from');
            my $hitEnd = $hsp->findvalue('./Hsp_hit-to');
            my $hitIdentity = $hsp->findvalue('./Hsp_identity');
            my $aligmentLength = $hsp->findvalue('./Hsp_align-len');
            my $querySequence = $hsp->findvalue('./Hsp_qseq');
            my $hitSequence = $hsp->findvalue('./Hsp_hseq');
            my $alignment = $hsp->findvalue('./Hsp_midline');
            print FH ($name, "\t", $pepID, "\t", $pepLen, "\t", $bitScore, "\t", $hitStart, "\t", $hitEnd, "\t", $hitIdentity, "\t", $aligmentLength, "\t", $querySequence,
            "\t", $hitSequence, "\t", $alignment, "\n");
           
        }

   }

}
