#!/usr/bin/env python

import os
import sys

file = sys.argv[1]
taxonFile = sys.argv[2]
outputFile = sys.argv[3]

peptideMatchTaxaOutputFastaFile = sys.argv[4]
smallPeptideOutputFastaFile = sys.argv[5]
peptideOutputFastaFile = sys.argv[6]

peptideTabFh = open(file)
taxonFH = open(taxonFile)

keepTaxa = []

for line in taxonFH:
  keepTaxa.append(int(line.strip()))

outPut = open(outputFile, 'w')
pepMatchTaxaOutputFasta = open(peptideMatchTaxaOutputFastaFile, 'w')
smallFastaOut = open(smallPeptideOutputFastaFile, 'w')
fastaOut = open(peptideOutputFastaFile, 'w')

print("##Accessions File", file=outPut)

for line in peptideTabFh:
    if not line.strip():
        continue

    splitLine = line.split("\t")
    iedbTaxId = int(splitLine[2])

    accessionNumber = splitLine[0]

    epitopeId = splitLine[1]
    peptide = splitLine[3]

    if len(peptide) < 5:
      print(f'>{epitopeId}\n{peptide}', file=smallFastaOut)
    else:
      print(f'>{epitopeId}\n{peptide}', file=fastaOut)

    if iedbTaxId in keepTaxa:
      print(accessionNumber, file=outPut)
      print(f">{splitLine[1]}\n{splitLine[3]}", file=pepMatchTaxaOutputFasta)

peptideTabFh.close()
outPut.close()

