#!/bin/bash
set -euo pipefail

conda activate R
STAG_MWC_OUTPUT_PATH=/home/groups/lu_kpmi/outputs/
../../kraken-biom/kraken_biom.py $STAG_MWC_OUTPUT_PATH/*/kraken2/*_bracken_species.kreport -o abundance.biom --fmt json
Rscript alpha_beta_diversity.R
