#!/usr/bin/env perl

use strict;
use XML::LibXML;
use Data::Dumper;
use open qw(:std :utf8);
use Getopt::Long;

my ($xmlFile, $outFasta);
&GetOptions("inputXml=s"=> \$xmlFile,
            "output=s"=> \$outFasta,
           ) ;
die("Please provide both input and out put") unless ($xmlFile & $outFasta);

my $xml = XML::LibXML->load_xml(location => $xmlFile);


my @epitopes = $xml->findnodes('/References/Reference/Epitopes/Epitope') ;

open(FH, '>>', $outFasta) or die $!;

foreach my $epitope (@epitopes) {
  my ($EpitopeId, $epitopeName, $GenBankId, $SourceOrganismId, $LinearSequence);
  $epitopeName = ($epitope->findvalue('./EpitopeName'));
  $EpitopeId = ($epitope->findvalue('./EpitopeId'));
  my ($EpitopeStructure) = $epitope->findnodes('EpitopeStructure/FragmentOfANaturalSequenceMolecule') ;
  next unless $EpitopeStructure;
  $GenBankId = ($epitope->findvalue('./EpitopeStructure/FragmentOfANaturalSequenceMolecule/SourceMolecule/GenBankId'));
  $SourceOrganismId = ($epitope->findvalue('./EpitopeStructure/FragmentOfANaturalSequenceMolecule/SourceOrganismId'));
  $LinearSequence = ($epitope->findvalue('./EpitopeStructure/FragmentOfANaturalSequenceMolecule/LinearSequence'));

  print FH (">",$epitopeName, "|" ,$EpitopeId , "|", $GenBankId, ,"|", "$SourceOrganismId","\n");
  print FH ($LinearSequence, "\n")

}

close(FH)