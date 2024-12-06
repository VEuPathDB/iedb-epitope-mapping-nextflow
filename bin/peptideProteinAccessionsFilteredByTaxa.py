#!/usr/bin/env python

import os
import sys

file = sys.argv[1]
taxonFile = sys.argv[2]
outputFile = sys.argv[3]
peptideOutputFastaFile = sys.argv[4]

peptideTabFh = open(file)
taxonFH = open(taxonFile)

keepTaxa = []

for line in taxonFH:
  keepTaxa.append(int(line.strip()))

outPut = open(outputFile, 'w')
pepOutputFasta = open(peptideOutputFastaFile, 'w')

print("##Accessions File", file=outPut)

for line in peptideTabFh:
    if not line.strip():
        continue

    splitLine = line.split("\t")
    iedbTaxId = int(splitLine[2])

    accessionNumber = splitLine[0]
    if iedbTaxId in keepTaxa:
      print(accessionNumber, file=outPut)
      print(f">{splitLine[1]}\n{splitLine[3]}", file=pepOutputFasta)

peptideTabFh.close()
outPut.close()

