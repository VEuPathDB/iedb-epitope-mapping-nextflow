#!/usr/bin/env bash

set -euo pipefail

sed -e s/"<References.*"/"<References>"/ $xmlFile > temp.xml

processXmlFile.pl --inputXml ./temp.xml --output ${sample_base}.fa