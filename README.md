# Immune Epitope analyses
This repository contains scripts and processes for the identification of genes whose protein products contains epitope sequences identified by  Immune Epitope Database and Analysis Resource (IEDB).


The analysis begin with processing of the epitopes taking a tab file containing the epitopes and reference proteome.  We divide the peptides into 3 categories:  very small peptides which require separate preprocessing of the reference genome, peptides annotated as matching the input taxon which we will allow up to one mismatch in the alignment, and the rest which we will do exact matching.  We also retrieve the original protein sequences from genbank for those peptides which match the input taxon.  Full peptide sequence are also matched against the reference proteome.

iedb PEPMatch (https://github.com/IEDB/PEPMatch) is used for peptide matching.


**<p align=left>Get Started</p>**
To run the work the following dependencies need to be downloaded and installed. 

* Docker 
> `https://docs.docker.com/engine/install/`
* Nextflow
> `curl https://get.nextflow.io | bash`

* Pull the GitHub with the command below
> `git pull https://github.com/VEuPathDB/antismash-nextflow.git`
<br> Then do > `nextflow run main.nf`  

* Alternatively the workflow can be run directly from the online which pull the repo and run it.
> `nextflow run VEuPathDB/antismash-nextflow -with-trace -c <config_file> -r main`
The `-c <config_file>` is the nextflow.config file, an example can found in the repo directory


***<p align=center>Nextflow workflow diagram</p>*** 
``` mermaid
flowchart TB
    subgraph " "
    v0["taxonID"]
    v2["peptideTabfile"]
    v8["refFasta"]
    v9["peptideTabfile"]
    v11["fasta"]
    v16["fasta"]
    v21["fasta"]
    end
    subgraph epitopeMapping
    v1([fetchTaxon])
    v3([peptideProteinAccessionsFilteredByTaxa])
    v5([fetchProtein])
    v10([iedbExactMatches])
    subgraph smallExactPepMatch
    v12([preprocess])
    v13([match])
    v14([filterResults])
    end
    subgraph exactPepMatch
    v17([preprocess])
    v18([match])
    v19([filterResults])
    end
    subgraph inexactForTaxaPeptidesPepMatch
    v22([preprocess])
    v23([match])
    v24([filterResults])
    end
    v28([mergeResultsFiles])
    v30([indexResults])
    v4(( ))
    v6(( ))
    v15(( ))
    v20(( ))
    v25(( ))
    
    end
    v0 --> v1
    v1 --> v3
    v1 --> v10
    v2 --> v3
    v3 --> v13
    v3 --> v4
    v3 --> v15
    v3 --> v20
    v4 --> v5
    v5 --> v6
    v8 --> v10
    v9 --> v10
    v6 --> v10
    v10 --> v28
    v11 --> v12
    v12 --> v13
    v13 --> v14
    v14 --> v25
    v16 --> v17
    v17 --> v18
    v15 --> v18
    v18 --> v19
    v19 --> v25
    v21 --> v22
    v22 --> v23
    v20 --> v23
    v23 --> v24
    v24 --> v25
    v25 --> v28
    v28 --> v30
```