#!/usr/bin/env python

import os
import sys

# short peptides need different pepMatch preprocessing (kmer=2)
shortPeptideCutoff = 5

file = sys.argv[1]
taxonFile = sys.argv[2]
outputFile = sys.argv[3]

peptideMatchTaxaOutputFastaFile = sys.argv[4]
smallPeptideOutputFastaFile = sys.argv[5]
otherPeptideOutputFastaFile = sys.argv[6]

peptideTabFh = open(file)
taxonFH = open(taxonFile)

keepTaxa = []

for line in taxonFH:
  keepTaxa.append(int(line.strip()))

outPut = open(outputFile, 'w')
pepMatchTaxaOutputFasta = open(peptideMatchTaxaOutputFastaFile, 'w')
smallFastaOut = open(smallPeptideOutputFastaFile, 'w')
otherFastaOut = open(otherPeptideOutputFastaFile, 'w')

print("##Accessions File", file=outPut)

for line in peptideTabFh:
    if not line.strip():
        continue

    splitLine = line.split("\t")
    iedbTaxId = int(splitLine[2])

    accessionNumber = splitLine[0]

    epitopeId = splitLine[1]
    peptide = splitLine[3]

    if len(peptide) < shortPeptideCutoff:
      print(f'>{epitopeId}\n{peptide}', file=smallFastaOut)
    elif iedbTaxId in keepTaxa:
      print(accessionNumber, file=outPut)
      print(f">{epitopeId}\n{peptide}", file=pepMatchTaxaOutputFasta)
    else:
      print(f'>{epitopeId}\n{peptide}', file=otherFastaOut)

peptideTabFh.close()
outPut.close()

