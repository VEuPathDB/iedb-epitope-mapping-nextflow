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
    for line in peptideTab:
        peptidesProperties = line.split()
        protein = peptidesProperties[0]
        peptideId = peptidesProperties[1]
        peptide = peptidesProperties[3]
        print(">", protein, "_", peptideId, sep="",file=fastaOut)
        print(peptide, file=fastaOut)

        if protein in peptideDic:
            peptideDic[protein].append(peptide)
        else:
            peptideDic[protein] = [peptide]    

    for refSeq in SeqIO.parse(refProteome, "fasta"):
        for pepSeq in SeqIO.parse(epitopeProtein, "fasta"):
            if refSeq.seq == pepSeq.seq: 
                peptideList = peptideDic.get(pepSeq.id)
                for pep in peptideList:
                    if pep in refSeq.seq:
                        print(refSeq.id,"\t" ,"Exact match","\t", pep, file=outPut)
                    else:
                        print(refSeq.id, "\t", "Not exact match", "\t", pep, file=outPut)

        else:
            for key in peptideDic:
                pepList = peptideDic[key]
                for pep in pepList:
                    if pep in refSeq.seq:
                        print(refSeq.id, "\t",  "Has a peptide in gene", "\t", pep, file=outPut)

    outPut.close()
    fastaOut.close()

if __name__ == "__main__":
   main(sys.argv[1:])
