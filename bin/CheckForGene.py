#!/usr/bin/env python

from Bio import SeqIO
import sys, getopt, re

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
        print ('CheckForGene.py -r <refProteome> -e <epitopeProtein> -l <epitopetab> -t <refTaxa> -p <peptideMatchOutput> -o <filteredPeptideFasta>')   
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('CheckForGene.py -r <refProteome> -e <epitopeProtein> -l <epitopetab> -t <refTaxa> -p <peptideMatchOutput> -o <filteredPeptideFasta>')
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
        elif opt in ("-o", "--filteredPeptideFasta"):
            filteredPeptideFasta = arg
   
    outPut = open(peptideMatchOutput, 'w')
    fastaOut = open(filteredPeptideFasta, 'w')
    peptideTab = open(epitopetab)
    
    referenceTaxa = [ ]

    try:
        with open(refTaxa) as taxaFile:
            for line in taxaFile:
        
                currentLine = line.strip()
                referenceTaxa.append(int(currentLine))
    except FileNotFoundError:
        print(print(f"File {refTaxa} not found!", file=sys.stderr))

    referenceTaxa = set(referenceTaxa)

    proteinDict = {}
    for refSeq in SeqIO.parse(epitopeProtein, "fasta"):
        proteinDict[refSeq.id] = refSeq.seq

    class Epitope:
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
    
    epitopeDict = {}
    try:
        with open(epitopetab) as pepFile:
            for lines in pepFile:
                line = lines.strip()
                peptidesProperties = line.split("\t")

                peptideId = peptidesProperties[1]
                iedbTaxon = int(peptidesProperties[2])
                peptide = peptidesProperties[3]
                protID = peptidesProperties[0]
                sequence = proteinDict.get(protID)
                if iedbTaxon in referenceTaxa: 
                    print(">", peptideId, sep="",file=fastaOut)
                    print(peptide, file=fastaOut)

                c1 = Epitope(peptideId, protID, iedbTaxon, peptide, sequence)
                epitopeDict[peptideId] = c1
    except FileNotFoundError:
         print(print(f"File {epitopetab} not found!", file=sys.stderr))

    same = 1
    different = 0

    for refSeq in SeqIO.parse(refProteome, "fasta"):    
        
        for epitopeInstance in epitopeDict:
            peptideID = epitopeDict[epitopeInstance].peptideID
            peptideProteinSeq = epitopeDict[epitopeInstance].sequence
            peptideTaxon = epitopeDict[epitopeInstance].peptideTaxon
            proteinID = epitopeDict[epitopeInstance].proteinID
            peptideSeq = epitopeDict[epitopeInstance].peptide
        
            if peptideTaxon in referenceTaxa:
                if refSeq.seq == peptideProteinSeq and str(peptideSeq) in refSeq.seq:
                    match=(re.search(str(peptideSeq), str(refSeq.seq)))
                    matchStart = match.start() + 1
                    matchEnd = match.end()
                    print(refSeq.id,"\t" ,same, "\t", same,"\t", same ,"\t", peptideSeq, "\t", peptideID,  "\t", proteinID, "\t", matchStart, "\t", matchEnd, "\t", peptideTaxon, "\t", file=outPut, sep="")

                elif refSeq.seq != peptideProteinSeq and peptideSeq in refSeq.seq:
                    match=(re.search(str(peptideSeq), str(refSeq.seq)))
                    matchStart = match.start() + 1
                    matchEnd = match.end()
                    print(refSeq.id, "\t", same, "\t", different,"\t", same, "\t", peptideSeq, "\t", peptideID, "\t", proteinID, "\t", matchStart, "\t", matchEnd,  "\t", peptideTaxon, "\t", file=outPut, sep="")
                    
            else:
                if refSeq.seq == peptideProteinSeq and peptideSeq in refSeq.seq:
                    match=(re.search(str(peptideSeq), str(refSeq.seq)))
                    matchStart = match.start() + 1
                    matchEnd = match.end()
                    print(refSeq.id,"\t" , same, "\t", same,"\t", different ,"\t", peptideSeq, "\t", peptideID,  "\t", proteinID, "\t", matchStart, "\t", matchEnd, "\t", peptideTaxon, file=outPut, sep="")
                elif refSeq.seq != peptideProteinSeq and peptideSeq in refSeq.seq:
                    match=(re.search(str(peptideSeq), str(refSeq.seq)))
                    matchStart = match.start() + 1
                    matchEnd = match.end()
                    print(refSeq.id, "\t", same, "\t", different,"\t", different , "\t", peptideSeq, "\t", peptideID, "\t", proteinID, "\t", matchStart, "\t", matchEnd, "\t", peptideTaxon, file=outPut, sep="")
                    
    
    outPut.close()
    fastaOut.close()
    peptideTab.close()

if __name__ == "__main__":
   main(sys.argv[1:])
