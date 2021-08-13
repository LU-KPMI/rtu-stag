#!/bin/bash
set -euo pipefail

rm -rf outputs
mkdir outputs
mkdir outputs/tables
mkdir outputs/pictures

python3 extract.py
echo "Setup done"

python3 process.py
echo "Analysis done"

python3 statistics.py
echo "Statistics done"

zip -r outputs.zip outputs
