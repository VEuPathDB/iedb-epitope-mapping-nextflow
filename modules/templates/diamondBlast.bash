#!/usr/bin/env bash

set -euo pipefail

diamond blastp -q ${query} -d ${db} -f 5 -o ${sample_base}.xml