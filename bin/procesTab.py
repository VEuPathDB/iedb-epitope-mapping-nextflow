#!/usr/bin/env python

import os
import sys

file = sys.argv[1]
peptideTab = open(file)
taxonFile = sys.argv[2]
taxonFH = open(taxonFile)

taxonList = []

for line in taxonFH:
  taxonList.append(int(line.strip()))

ncbiDict = {}
for line in peptideTab:
  tax = int(line.split("\t")[2])
  ncbiID = line.split("\t")[0].split('.')[0]
  if tax in taxonList:
    ncbiDict[ncbiID] = ncbiID

outPut = open("peptideProteins.txt", 'w') 
for key in ncbiDict:
  print(key, file=outPut)

peptideTab.close()
outPut.close()


com = """
  fileItemString=$(cat  peptideProteins.txt | tr "\n" ",")
  efetch -db protein -id $fileItemString -format fasta >> pepProtein.fasta
  sleep 5
"""
os.system(com)