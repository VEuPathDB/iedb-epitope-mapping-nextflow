#!/usr/bin/env bash

set -euo pipefail

blastp -query ${query} -outfmt 6 -db ${db}/*fa -out ${sample_base}.txt