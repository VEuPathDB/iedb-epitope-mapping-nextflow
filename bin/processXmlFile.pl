#!/usr/bin/perl

use strict;
use XML::LibXML;
use Data::Dumper;
use open qw(:std :utf8);
use Getopt::Long;


my $filename = $ARGV[0];

my $xml = XML::LibXML->load_xml(location => $filename);


my @epitopes = $xml->findnodes('/References/Reference/Epitopes/Epitope') ;

my $outFile = $ARGV[1];
open(FH, '>>', $outFile) or die $!;

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