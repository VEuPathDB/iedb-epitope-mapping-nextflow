#!/usr/bin/env bash

set -euo pipefail

mkdir BlastDB

cat ${fasta} > BlastDB/${sample_base}.fa

makeblastdb -in BlastDB/${sample_base}.fa   -title "Cookbook demo" -dbtype prot  