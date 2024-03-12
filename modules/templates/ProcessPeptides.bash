#!/usr/bin/env bash

set -euo pipefail


CheckForGene.py -r ${refFasta} \
    -l ${pepTab}  \
    -e ${pepProtfasta} \
    -t ${taxon} \
    -p $params.peptideMatchResults \
    -o $params.peptidesFilteredBySpeciesFasta
