# Immune Epitope analyses
This repository consists of scripts and procedure for the identification of genes whose protein products show similarity to an epitope identified by  Immune Epitope Database and Analysis Resource (IEDB).


The analysis begin with processing of the Epitopes talking a tab file containing the Epitopes and reference proteome. The first step identify epitopes that have exact match in the reference proteome and weather the gene source of the epitope is identical to the reference protein in the proteome. 

The Second step use Blast to identify hit of the Epitopes in the reference proteome. 

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
