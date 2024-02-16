#!/usr/bin/env bash

set -euo pipefail

mkdir db
cat ${fasta} | sed -e 's/ le.*//g' > db/${sample_base}.fasta

makeblastdb -in db/${sample_base}.fasta   -dbtype prot  