#!/usr/bin/env bash

set -euo pipefail

blastp -query ${query} -num_descriptions 3 -num_alignments 3 -outfmt 6 -db ${db}/*fa -out ${sample_base}.txt
#-outfmt 6 
# -v == -num_descriptions
