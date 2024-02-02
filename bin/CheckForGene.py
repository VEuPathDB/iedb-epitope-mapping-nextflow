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
    

    # FIXME:  make this a set instead of list
    referenceTaxa = [ ]

    try:
        with open(refTaxa) as taxaFile:
            for line in taxaFile:
        
                currentLine = line.strip()
                referenceTaxa.append(int(currentLine))
    except FileNotFoundError:
        print(print(f"File {refTaxa} not found!", file=sys.stderr))

    referenceTaxa = set(referenceTaxa)
    
    peptideDic = {}
    peptideNames = {}
    peptideTaxa = {}
    for lines in peptideTab:
        line = lines.strip()
        peptidesProperties = line.split("\t")
        # reanme this iedbProtein
        protein = peptidesProperties[0]
        peptideId = peptidesProperties[1]
        # iedbTaxon not taxa
        taxa = int(peptidesProperties[2])
        peptide = peptidesProperties[3]
        if taxa in referenceTaxa: 
            print(">", peptideId, sep="",file=fastaOut)
            print(peptide, file=fastaOut)



        # TODO: make a Peptide class which has instance variables for taxon,iedbid and proteinAccession
        if protein in peptideDic:
            peptideDic[protein].append(peptide)
        else:
            peptideDic[protein] = [peptide]

        # FIXME:  check but error if peptide is found more than once
        if peptide in peptideNames:
            peptideNames[peptide] = peptideId
        else:
            peptideNames[peptide] = peptideId

        # FIXME:  always set. no need for if
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
                # FIXME:  set all variables in if/else and print one time at bottom
                if taxa in referenceTaxa: 
                    if refSeq.seq == pepSeq.seq and pep in refSeq.seq:
                        match=(re.search(str(pep), str(refSeq.seq)))
                        matchStart = match.start() + 1
                        matchEnd = match.end()
                        print(refSeq.id,"\t" ,1, "\t", 1,"\t", 1 ,"\t", pep, "\t", pepID,  "\t", pepSeq.id, "\t", matchStart, "\t", matchEnd, "\t", taxa, "\t", file=outPut, sep="")

                    elif refSeq.seq != pepSeq.seq and pep in refSeq.seq:
                        match=(re.search(str(pep), str(refSeq.seq)))
                        matchStart = match.start() + 1
                        matchEnd = match.end()
                        print(refSeq.id, "\t", 1, "\t", 0,"\t", 1, "\t", pep, "\t", pepID, "\t", pepSeq.id, "\t", matchStart, "\t", matchEnd,  "\t", taxa, "\t", file=outPut, sep="")
                        
                else:
                    if refSeq.seq == pepSeq.seq and pep in refSeq.seq:
                        match=(re.search(str(pep), str(refSeq.seq)))
                        matchStart = match.start() + 1
                        matchEnd = match.end()
                        print(refSeq.id,"\t" , 1, "\t", 1,"\t", 0 ,"\t", pep, "\t", pepID,  "\t", pepSeq.id, "\t", matchStart, "\t", matchEnd, "\t", taxa, file=outPut, sep="")
                p
                    elif refSeq.seq != pepSeq.seq and pep in refSeq.seq:
                        match=(re.search(str(pep), str(refSeq.seq)))
                        matchStart = match.start() + 1
                        matchEnd = match.end()
                        print(refSeq.id, "\t", 1, "\t", 0,"\t", 0 , "\t", pep, "\t", pepID, "\t", pepSeq.id, "\t", matchStart, "\t", matchEnd, "\t", taxa, file=outPut, sep="")
                        
        
    outPut.close()
    fastaOut.close()
    peptideTab.close()

if __name__ == "__main__":
   main(sys.argv[1:])
