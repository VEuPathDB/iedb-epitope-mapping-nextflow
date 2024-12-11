# Immune Epitope analyses
This repository contains scripts and processes for the identification of genes whose protein products contains epitope sequences identified by  Immune Epitope Database and Analysis Resource (IEDB).


The analysis begin with processing of the epitopes taking a tab file containing the epitopes and reference proteome.  We divide the peptides into 3 categories:  very small peptides which require separate preprocessing of the reference genome, peptides annotated as matching the input taxon which we will allow up to one mismatch in the alignment, and the rest which we will do exact matching.  We also retrieve the original protein sequences from genbank for those peptides which match the input taxon.  Full peptide sequence are also matched against the reference proteome.

iedb pepMatch (https://github.com/IEDB/PEPMatch) is used for peptide matching.

