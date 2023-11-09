#!/usr/bin/env bash

set -euo pipefail

file1=\$(cat ${url})

mkdir iedb_export

curl -JL  -o iedb_export.zip  "\${file1}" 

unzip iedb_export.zip -d iedb_export