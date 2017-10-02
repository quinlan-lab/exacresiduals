#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ./%j-%N.out
#SBATCH -e ./%j-%N.err
#SBATCH --time=1:00:00
#set -exo pipefail -o nounset

bash regions.sh -c -w -v gnomAD10x.5 -d 10 -d 0.5 -x data/segmental.bed.gz -x data/self-chains.gt90.bed.gz
