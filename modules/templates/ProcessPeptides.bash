#!/usr/bin/env bash

set -euo pipefail

CheckForGene.py -r ${refFasta} \
    -l ${pepTab}  \
    -e ${pepProtfasta} \
    -t ${taxon} \
    -p ${peptideMatchResultsOutput} \
    -o ${peptidesFilteredBySpeciesFastaOutput}
