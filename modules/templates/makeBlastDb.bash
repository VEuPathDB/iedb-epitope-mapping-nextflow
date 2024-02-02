#!/usr/bin/env bash

set -euo pipefail

#mkdir BlastDB
cat ${fasta} | sed -e 's/ le.*//g' > ${sample_base}.fasta

makeblastdb -in ${sample_base}.fasta   -dbtype prot  