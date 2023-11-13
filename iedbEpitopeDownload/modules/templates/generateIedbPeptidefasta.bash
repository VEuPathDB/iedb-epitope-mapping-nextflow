#!/usr/bin/env bash

set -euo pipefail

sed -e s/"<References.*"/"<References>"/ $xmlFile > temp.xml

processXmlFile.pl --inputXml ./temp.xml --output temp.fa 
sed -e 's/ //g' temp.fa >  ${sample_base}.fa

rm temp.fa