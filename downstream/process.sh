#!/bin/bash
set -euo pipefail

python3 ./extract.py
echo "Setup done"

python3 ./generate_abundance_table.py > abundances.csv
echo "Abundance table generated"

python3 ./calc_richness.py > richness.csv
echo "Richness calculated"

python3 ./rarefaction.py
echo "Rarefaction generated"

python3 ./graphs.py
echo "Richness summary visualization complete"

python3 ./pairwise.py
echo "Pairwise analysis done"

python3 ./enterotypes.py
echo "Basic enterotype analysis done"
