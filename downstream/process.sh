#!/bin/bash
set -euo pipefail

python3 extract.py
echo "Setup done"

python3 analysis.py
echo "Analysis done"
