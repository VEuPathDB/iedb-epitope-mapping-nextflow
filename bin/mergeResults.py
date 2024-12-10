#!/usr/bin/env python3

import argparse,re

def parseArguments():
    parser = argparse.ArgumentParser(description="Matches outputs from pepMatch and exact full Protein matches into a single gff file.")
    parser.add_argument("--pepMatchResults", required=True, help="Path to peptide match results")
    parser.add_argument("--fullSeqMatchResults", required=True, help="Path to full sequence match results")
    parser.add_argument("--peptideTab", required=True, help="Path to peptide tab file")
    parser.add_argument("--outputFile", required=True, help="Path to output file")
    parser.add_argument("--distinctPeptides", required=True, help="Path to distinct peptides file")
    parser.add_argument("--nonTaxaShortPeptideCutoff", required=True, help="Peptides less than this length from other taxa will not be output")
    return parser.parse_args()

def parseDistinctPeptides(file):
    rv = {}
    with open(file) as fh:
        for line in fh:
            peptide = line.strip()
            rv[peptide] = 1
    fh.close()
    return(rv)

def parsePeptideToIedb(file, keep):
    rv = {}
    with open(file) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            peptide = fields[3]

            if peptide in keep:
                iedbId = fields[1]
                rv[peptide] = iedbId
    fh.close()
    return(rv)

def parseFullSeqMatchResults(file):
    rv = {}
    with open(file) as fh:
        for line in fh:
            peptide, protein, iedb, taxonMatch, proteinMatch = line.strip().split(",")

            if protein == "":
                protein = "NA"

            key = peptide + "_" + protein;

            vals = [taxonMatch, proteinMatch]
            rv[key] = vals
    fh.close()
    return(rv)

def main(args):

    gffSource = "veupathdb"
    gffType = "peptide"

    nonTaxaShortPeptideCutoff = int(args.nonTaxaShortPeptideCutoff)

    distinctPeptides = parseDistinctPeptides(args.distinctPeptides)
    peptideToIedb = parsePeptideToIedb(args.peptideTab, distinctPeptides)
    fullSeqMatchResults = parseFullSeqMatchResults(args.fullSeqMatchResults)

    removeSuffixPattern = r"\.\d+$"

    out = open(args.outputFile, 'w')
    with open(args.pepMatchResults) as fh:
        for line in fh:
            peptide, protein, mismatches, start, end = line.strip().split(",")

            # pepMatch preprocess db adds a ".1 suffix to proteins; remove it
            trimmedProtein = re.sub(removeSuffixPattern, "", protein)

            # this will throw KeyError exception if not found
            iedbId = peptideToIedb[peptide]

            key = peptide + "_" + trimmedProtein;
            altKey = peptide + "_NA"

            taxonMatch = 0
            fullProteinMatch = 0
            if key in fullSeqMatchResults:
                tm, pm = fullSeqMatchResults[key]
                taxonMatch = tm
                fullProteinMatch = pm
            elif altKey in fullSeqMatchResults:
                tm, pm = fullSeqMatchResults[altKey]
                taxonMatch = tm
                fullProteinMatch = pm
            else:
                taxonMatch = 0
                fullProteinMatch = 0

            gffId = trimmedProtein + "_" + start + "_" + end
            attributes = f"ID={gffId};iedb={iedbId};matchesTaxon={taxonMatch};matchesFullLengthProtein={fullProteinMatch};mismatches={mismatches};peptide={peptide}"

            gffRow = f"{trimmedProtein}\t{gffSource}\t{gffType}\t{start}\t{end}\t.\t.\t.\t{attributes}"

            if not (taxonMatch != 1 and len(peptide) < nonTaxaShortPeptideCutoff):
                print(gffRow, file=out)

    fh.close()

    out.close()

if __name__ == "__main__":
    main(parseArguments())
