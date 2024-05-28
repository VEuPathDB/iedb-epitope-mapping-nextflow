# Immune Epitope analyses
This repository contains scripts and processes for the identification of genes whose protein products contains epitope sequences identified by  Immune Epitope Database and Analysis Resource (IEDB).


The analysis begin with processing of the epitopes taking a tab file containing the epitopes and reference proteome. The first step identify epitopes that have exact match in the reference proteome and weather the gene source of the epitope is identical to the reference protein in the proteome. During this step a fasta file of the epitopes is generated to blast against the reference proteome. 

The Second step uses Blast to identify epitopes hits in the reference proteome. Default blast option is NCBI plastp and user can set blast to diamond in the nextflow config file. The blast is perform only in cases where the epitope source taxon agree with the reference taxa. 

The steps in the nextflow graph is shown below. 


**<p align=left>Get Started</p>**
To run the work the following dependencies need to be install
* Docker
> `https://docs.docker.com/engine/install/`
* Nextflow
> `curl https://get.nextflow.io | bash`

* The pull the git hub repo using the following command
> `git pull https://github.com/VEuPathDB/iedb-epitope-mapping-nextflow.git`

* Alternatively the workflow can be run directly using nextflow which pull down the repo. 
> `nextflow run VEuPathDB/iedb-epitope-mapping-nextflow -with-trace -c  <config_file> -r main`



***<p align=center>Nextflow workflow</p>*** 

```mermaid

flowchart TD
    p0((Channel.fromPath))
    p1((Channel.fromPath))
    p2((Channel.fromPath))
    p3[epitopesBlast:peptideSimilarity]
    p4[epitopesBlast:makeBlastDatabase]
    p5[epitopesBlast:blastSeq]
    p6[epitopesBlast:processXml]
    p7[epitopesBlast:mergeeResultsFiles]
    p8(( ))
    p0 -->|refFasta| p3
    p1 -->|peptidesTab| p3
    p2 -->|peptidesGeneFasta| p3
    p3 --> p4
    p3 --> p7
    p4 --> p5
    p0 -->|refFasta| p5
    p5 --> p6
    p6 --> p7
    p7 --> p8


```
