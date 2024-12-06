#!/usr/bin/env bash

set -euo pipefail

mkdir db
cp ${fasta} db/${sample_base}.fasta

# TODO:  Not sure what this sed was doing here!
# cat ${fasta} | sed -e 's/ le.*//g' > db/${sample_base}.fasta

makeblastdb -in db/${sample_base}.fasta   -dbtype prot  
