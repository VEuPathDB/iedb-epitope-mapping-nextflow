#!/usr/bin/env bash

set -euo pipefail

blastp -query ${query} -task blastp-short -word_size 2  -evalue 100 -max_target_seqs 10000 -outfmt 5 -db ${db}/*fasta -out ${sample_base}.xml

