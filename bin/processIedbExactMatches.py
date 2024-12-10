#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys, getopt, re


def main(argv):
    refProteome = ''
    epitopeProtein = ''
    epitopetab = ''
    peptideMatchOutput = ''
    refTaxa = ''

    try:
        opts, args = getopt.getopt(argv,"hr:e:l:t:p:s:f",["refProteome=","epitopeProtein=","epitopetab=", "refTaxa=", "peptideMatchOutput="])
    except getopt.GetoptError:
        print ('processIedbExactMatches.py -r <refProteome> -e <epitopeProtein> -l <epitopetab> -t <refTaxa> -p <peptideMatchOutput>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('processIedbExactMatches.py -r <refProteome> -e <epitopeProtein> -l <epitopetab> -t <refTaxa> -p <peptideMatchOutput>')
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

#    peptideTab = open(epitopetab)
    
    referenceTaxa = {}
    with open(refTaxa) as taxaFile:
        for line in taxaFile:
            currentLine = line.strip()
            referenceTaxa[currentLine] = 1


    pepProtDict = {}
    for pepProtSeq in SeqIO.parse(epitopeProtein, "fasta"):
        pepProtDict[pepProtSeq.id] = pepProtSeq.seq

    refSeqDict = {}
    for refSeq in SeqIO.parse(refProteome, "fasta"):
        refSeqDict[refSeq.id] = refSeq.seq

    try:
        with open(epitopetab) as pepFile:
            for lines in pepFile:
                line = lines.strip()
                peptidesProperties = line.split("\t")

                peptideId = peptidesProperties[1]
                peptide = peptidesProperties[3]

                accessionNumber = peptidesProperties[0]

                iedbTaxon = peptidesProperties[2]
                sequence = pepProtDict.get(accessionNumber)

                TaxonMatch = 0
                if iedbTaxon in referenceTaxa:
                    TaxonMatch = 1

                    seenPerfectMatch = False
                    for refSeqId, refSeqSeq in refSeqDict.items():
                        if refSeqSeq == sequence:
                            seePerfectMatch = True
                            print(f'{peptide},{refSeqId},{peptideId},{TaxonMatch},1', file=outPut)

                    if seenPerfectMatch is False:
                        print(f'{peptide},,{peptideId},{TaxonMatch},0', file=outPut)

    except FileNotFoundError:
        print(f"File {epitopetab} not found!", file=sys.stderr)

    outPut.close()

if __name__ == "__main__":
    main(sys.argv[1:])
