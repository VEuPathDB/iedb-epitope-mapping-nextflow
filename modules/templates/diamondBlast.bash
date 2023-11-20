#!/usr/bin/env bash

set -euo pipefail

diamond blastp -q ${query} -d ${db} -o ${sample_base}.txt