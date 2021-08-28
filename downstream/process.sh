#!/bin/bash
set -euo pipefail

rm -rf outputs
mkdir outputs
mkdir outputs/tables
mkdir outputs/pictures

python3 extract.py
echo "Setup done"

for f in outputs/bracken_output/*kreport; do
    python kreport_to_krona.py < $f > ${f%.kreport}.krona
done
ktImportText -o outputs/taxonomy.html outputs/bracken_output/*.krona
echo "Taxonomy kronagraph generated"

python3 process.py
echo "Analysis done"

python3 statistics.py
echo "Statistics done"

zip -r outputs.zip outputs
