#!/bin/bash
#SBATCH --account=quinlan-kp
#SBATCH --partition=quinlan-kp
#SBATCH -o ../%j-%N.out
#SBATCH -e ../%j-%N.err
#SBATCH --time=1:00:00
#set -exo pipefail -o nounset

bash regions.sh -c -w -d exacv1newweight #newweight30x.5 #exacv1newweight #parallel #exacv1final
