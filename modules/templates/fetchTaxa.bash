#!/usr/bin/env bash

set -euo pipefail

esearch -db taxonomy \
        -query "txid${taxonID}[Subtree]" \
        | efetch -format xml | xtract \
        -pattern TaxaSet -block "*/Taxon" -tab "\n" \
        -element TaxId > TaxaList.txt
