#!/usr/bin/env python

similarityFile = "/Users/saikouybah/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/IedbAnalysis/Result/PeptideGene.txt"
blastFile = "/Users/saikouybah/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/IedbAnalysis/Result/BlastOut/AnnotatedProteins.txt"

blastData = open(blastFile)
exactMatchData = open(similarityFile)


exactMatchDic = {}
for line in exactMatchData:
    blastValue = line.strip().split("\t")
    geneName = blastValue[0] 
    peptide = blastValue[3]
    dic_value = geneName + "_" + peptide
    exactMatchDic[dic_value] = blastValue



blastResultDic = {}
for line in blastData:
    blastValue = line.strip().split("\t")
    geneName = blastValue[0] 
    peptide = blastValue[2]
    dic_value = geneName + "_" + peptide
    print(dic_value, line.strip())