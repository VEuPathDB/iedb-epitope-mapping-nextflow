#!/usr/bin/env python

from Bio import SeqIO
import sys, getopt

def main(argv):
    refProteome = ''
    epitopeProtein = ''
    epitopetab = ''
    outfile = 'PeptideGene.txt'
    outfasta = 'Peptides.fasta'

    try:
        opts, args = getopt.getopt(argv,"hr:e:l:",["refProteome=","epitopeProtein=","epitopetab="])
    except getopt.GetoptError:
        print ('CheckForGene.py -r <refProeteom> -e <refProeteom> -l <epitopetab>')   
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('CheckForGene.py -r <refProeteom> -e <refProeteom> -l <epitopetab>')
            sys.exit()
        elif opt in ("-r", "--refProteome"):
         refProteome = arg
        elif opt in ("-e", "--epitopeProtein"):
            epitopeProtein = arg
        elif opt in ("-l", "--epitopetab"):
            epitopetab = arg
   
    outPut = open(outfile, 'w')
    fastaOut = open(outfasta, 'w')
    peptideTab = open(epitopetab)

    peptideDic = {}
    peptideNames = {}
    for lines in peptideTab:
        line = lines.strip()
        peptidesProperties = line.split("\t")
        protein = peptidesProperties[0]
        peptideId = peptidesProperties[1]
        peptide = peptidesProperties[3]
        print(">", peptideId, sep="",file=fastaOut)
        print(peptide, file=fastaOut)

        if protein in peptideDic:
            peptideDic[protein].append(peptide)
        else:
            peptideDic[protein] = [peptide]

        if peptide in peptideNames:
            peptideNames[peptide] = peptideId
        else:
            peptideNames[peptide] = peptideId

    for refSeq in SeqIO.parse(refProteome, "fasta"):
        for pepSeq in SeqIO.parse(epitopeProtein, "fasta"):
            peptideList = peptideDic.get(pepSeq.id)
            for pep in peptideList:
                pepID = peptideNames.get(pep)
                if refSeq.seq == pepSeq.seq and pep in refSeq.seq:
                    print(refSeq.id,"\t" ,"Peptide and protein match","\t", pep, "\t", pepID,  "\t", pepSeq.id, file=outPut, sep="")
                elif pep in refSeq.seq:
                   print(refSeq.id, "\t", "Only peptide match", "\t", pep, "\t", pepID, "\t", pepSeq.id, file=outPut, sep="")

    outPut.close()
    fastaOut.close()

if __name__ == "__main__":
   main(sys.argv[1:])
