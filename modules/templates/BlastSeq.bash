#!/usr/bin/env bash

set -euo pipefail

blastp -query ${query} -task blastp-short -word_size 4  -num_descriptions 3 -outfmt 6 -db ${db}/*fa -out ${sample_base}.txt

