#!/bin/bash
set -exo pipefail -o nounset

date=2016_06_15
mkdir -p results/$date/
python exac-regions.py | bgzip -c > results/$date/regs.bed.gz
python resid-plot.py results/$date/regs.bed.gz > results/$date/resids.txt
