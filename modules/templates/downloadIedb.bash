#!/usr/bin/env bash

set -euo pipefail

file1=\$(cat ${url})

curl -OJL   "\${file1}" -o iedb_export1

unzip 