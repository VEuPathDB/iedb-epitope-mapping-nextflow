#!/usr/bin/env bash

set -euo pipefail

diamond makedb --in ${fasta} -d ${sample_base}