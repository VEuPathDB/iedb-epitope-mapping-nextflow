#!/usr/bin/env python

from Bio import SeqIO
import sys, getopt, re

def main(argv):
    refProteome = ''
    epitopeProtein = ''
    epitopetab = ''
    outfile = 'PeptideGene.txt'
    outfasta = 'Peptides.fasta'
    refTaxa = 0

    try:
        opts, args = getopt.getopt(argv,"hr:e:l:t:",["refProteome=","epitopeProtein=","epitopetab=", "refTaxa="])
    except getopt.GetoptError:
        print ('CheckForGene.py -r <refProeteom> -e <epitopeProtein> -l <epitopetab> -t <refTaxa>')   
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('CheckForGene.py -r <refProeteom> -e <epitopeProtein> -l <epitopetab> -t <refTaxa>')
            sys.exit()
        elif opt in ("-r", "--refProteome"):
            refProteome = arg
        elif opt in ("-e", "--epitopeProtein"):
            epitopeProtein = arg
        elif opt in ("-l", "--epitopetab"):
            epitopetab = arg
        elif opt in ("-t", "--refTaxa"):
            refTaxa = arg
   
    outPut = open(outfile, 'w')
    fastaOut = open(outfasta, 'w')
    peptideTab = open(epitopetab)
    referenceTaxa = refTaxa
    peptideDic = {}
    peptideNames = {}
    peptideTaxa = {}
    for lines in peptideTab:
        line = lines.strip()
        peptidesProperties = line.split("\t")
        protein = peptidesProperties[0]
        peptideId = peptidesProperties[1]
        taxa = peptidesProperties[2]
        peptide = peptidesProperties[3]
        if int(referenceTaxa) == int(taxa):
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
        
        if peptide in peptideTaxa:
            peptideTaxa[peptide] = taxa
        else:
            peptideTaxa[peptide] = taxa

    for refSeq in SeqIO.parse(refProteome, "fasta"):
        for pepSeq in SeqIO.parse(epitopeProtein, "fasta"):
            peptideList = peptideDic.get(pepSeq.id)
           
            for pep in peptideList:
                pepID = peptideNames.get(pep)
                taxa = peptideTaxa.get(pep)
                if int(taxa) == int(referenceTaxa):
                    if refSeq.seq == pepSeq.seq and pep in refSeq.seq:
                        match=(re.search(str(pep), str(refSeq.seq)))
                        matchStart = match.start() + 1
                        matchEnd = match.end()
                        print(refSeq.id,"\t" ,"Yes", "\t", "Yes","\t", "Yes" ,"\t", pep, "\t", pepID,  "\t", pepSeq.id, "\t", matchStart, "\t", matchEnd, file=outPut, sep="")
                        #"Peptide and protein match; same strain"
                    elif refSeq.seq != pepSeq.seq and pep in refSeq.seq:
                        match=(re.search(str(pep), str(refSeq.seq)))
                        matchStart = match.start() + 1
                        matchEnd = match.end()
                        print(refSeq.id, "\t", "Yes", "\t", "No","\t", "Yes", "\t", pep, "\t", pepID, "\t", pepSeq.id, "\t", matchStart, "\t", matchEnd, file=outPut, sep="")
                        #"Only peptide match; same strain"
                else:
                    if refSeq.seq == pepSeq.seq and pep in refSeq.seq:
                        match=(re.search(str(pep), str(refSeq.seq)))
                        matchStart = match.start() + 1
                        matchEnd = match.end()
                        print(refSeq.id,"\t" , "Yes", "\t", "Yes","\t", "No" ,"\t", pep, "\t", pepID,  "\t", pepSeq.id, "\t", matchStart, "\t", matchEnd, file=outPut, sep="")
                        #"Peptide and protein match; different strain"
                    elif refSeq.seq != pepSeq.seq and pep in refSeq.seq:
                        match=(re.search(str(pep), str(refSeq.seq)))
                        matchStart = match.start() + 1
                        matchEnd = match.end()
                        print(refSeq.id, "\t", "Yes", "\t", "No","\t", "No" , "\t", pep, "\t", pepID, "\t", pepSeq.id, "\t", matchStart, "\t", matchEnd, file=outPut, sep="")
                        #"Only peptide match; different strain"
        
    outPut.close()
    fastaOut.close()
    peptideTab.close()

if __name__ == "__main__":
   main(sys.argv[1:])
