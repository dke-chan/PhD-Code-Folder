#!/bin/bash -e

#SBATCH --job-name    whaleFit
#SBATCH --time        1:00:00
#SBATCH --mem         2000MB
#SBATCH --partition   large
#SBATCH --mail-user   dcha704@aucklanduni.ac.nz
#SBATCH --mail-type   END

module load R/4.1.0-gimkl-2020a
module load admb

srun R --slave < 04B-core-scr-fit.R
