#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ../%j-%N.out
#SBATCH -e ../%j-%N.err
#SBATCH --time=6:00:00
#set -exo pipefail -o nounset

bash regions.sh -c -w -d 2017_04_09 #noasjfin #nontcga
