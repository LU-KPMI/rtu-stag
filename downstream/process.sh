#!/bin/bash
set -euo pipefail

python3 extract.py
echo "Setup done"

python3 process.py
echo "Analysis done"

python3 statistics.py
echo "Statistics done"
