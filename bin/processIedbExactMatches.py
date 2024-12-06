#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys, getopt, re


# This piece of code check if peptide are present in a given protein, weather it is an exact match, does the peptide source protein (sequence) match a refence protein 
# and does the peptide source protein Taxon (NCBI) matches a given NCBI taxon.


# Function below get the peptide start and end.
def getPeptideMatches(peptide, refSeq):
    match=re.search(peptide, refSeq)
    matchStart = match.start() + 1
    matchEnd = match.end()
    return(matchStart, matchEnd)

class Epitope:
        """
        This class define an epitope instance to contain information about an epitope. It has the iedb epitope ID, protein ID (protein source), 
        peptide taxon (NCBI taxon of the peptide source protein) and the peptide sequence. 
        """
        def __init__(self, peptideID, proteinID, peptideTaxon, peptide, sequence):
            self._peptideID = peptideID
            self._proteinID = proteinID
            self._peptideTaxon = peptideTaxon
            self._peptide = peptide
            self._sequence = sequence
        
        @property
        def peptideID(self):
            return self._peptideID   

        @property
        def proteinID(self):
            return self._proteinID
        
        @property
        def peptideTaxon(self):
            return self._peptideTaxon  

        @property
        def peptide(self):
            return self._peptide
        
        @property
        def sequence(self):
            return self._sequence

def main(argv):
    refProteome = ''
    epitopeProtein = ''
    epitopetab = ''
    peptideMatchOutput = ''
    filteredPeptideFasta = ''
    refTaxa = ''

    try:
        opts, args = getopt.getopt(argv,"hr:e:l:t:p:o:",["refProteome=","epitopeProtein=","epitopetab=", "refTaxa=", "peptideMatchOutput=", "filteredPeptideFasta="])
    except getopt.GetoptError:
        print ('processIedbExactMatches.py -r <refProteome> -e <epitopeProtein> -l <epitopetab> -t <refTaxa> -p <peptideMatchOutput> -o <filteredPeptideFasta>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('processIedbExactMatches.py -r <refProteome> -e <epitopeProtein> -l <epitopetab> -t <refTaxa> -p <peptideMatchOutput> -o <filteredPeptideFasta>')
            sys.exit()
        elif opt in ("-r", "--refProteome"):
            refProteome = arg
        elif opt in ("-e", "--epitopeProtein"):
            epitopeProtein = arg
        elif opt in ("-l", "--epitopetab"):
            epitopetab = arg
        elif opt in ("-t", "--refTaxa"):
            refTaxa = arg
        elif opt in ("-p", "--peptideMatchOutput"):
            peptideMatchOutput = arg

    outPut = open(peptideMatchOutput, 'w')

    peptideTab = open(epitopetab)
    
    referenceTaxa = []
    with open(refTaxa) as taxaFile:
        for line in taxaFile:
            currentLine = line.strip()
            referenceTaxa.append(currentLine)


    pepProtDict = {}
    for pepProtSeq in SeqIO.parse(epitopeProtein, "fasta"):
        pepProtDict[pepProtSeq.id] = pepProtSeq.seq
    
    epitopeDict = {}
    try:
        with open(epitopetab) as pepFile:
            for lines in pepFile:
                line = lines.strip()
                peptidesProperties = line.split("\t")

                peptideId = peptidesProperties[1]
                iedbTaxon = peptidesProperties[2]
                peptide = peptidesProperties[3]
                protID = peptidesProperties[0]
                sequence = pepProtDict.get(protID)

                epitopeDict[peptideId] = Epitope(peptideId, protID, iedbTaxon, peptide, sequence)
    except FileNotFoundError:
        print(f"File {epitopetab} not found!", file=sys.stderr)

    for refSeq in SeqIO.parse(refProteome, "fasta"):    

        for epitopeKey in epitopeDict:

            epitope = epitopeDict[epitopeKey]

            pepLen = len(epitope.peptide)

            peptideMatch = 0
            proteinMatch = 0
            TaxonMatch = 0
            matchStart = ''
            matchEnd = ''

            if epitope.peptideTaxon in referenceTaxa:
                TaxonMatch = 0

            if epitope.peptide in refSeq.seq:
                peptideMatch = 1
                matchStart, matchEnd = getPeptideMatches(epitope.peptide, str(refSeq.seq))
                if refSeq.seq == epitope.sequence:
                    proteinMatch = 1

            if peptideMatch == 0 and proteinMatch == 0 and TaxonMatch == 0:
                pass
            else:
                print(f'{refSeq.id}\t{peptideMatch}\t{proteinMatch}\t{TaxonMatch}\t{pepLen}\t{epitope.peptide}\t{epitope.peptideID}\t{epitope.proteinID}\t{matchStart}\t{matchEnd}\t{epitope.peptideTaxon}', file=outPut)
    
    outPut.close()

    peptideTab.close()

if __name__ == "__main__":
    main(sys.argv[1:])
