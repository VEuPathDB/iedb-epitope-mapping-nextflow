#!/usr/bin/env bash

set -euo pipefail

sed -e s/"<References.*"/"<References>"/ $xmlFile > temp.xml

processXmlFile.pl ./temp.xml ${sample_base}.fa