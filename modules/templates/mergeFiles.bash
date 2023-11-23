#!/usr/bin/env bash

set -euo pipefail

mergeBlastAndExactMatch.pl --exactMatchFiles ${exactMatch} --blastOutput ${balst}
