#!/usr/bin/env bash

set -euo pipefail

mergeBlastAndExactMatch.pl  --exactMatchFile ${exactMatch} \
                            --blastFile ${balstOutput} \
                            --outputFile ${peptideMatchBlastCombiedResults}
