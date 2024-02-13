#!/usr/bin/env bash

set -euo pipefail

parseXML.pl --inputXml ${xml} \
    --outputFile ${sample_base}.txt
