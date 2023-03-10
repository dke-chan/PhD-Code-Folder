#!/bin/bash -e

#SBATCH --job-name    fascrSim
#SBATCH --time        0:30:00
#SBATCH --mem         1500MB
#SBATCH --partition   large
#SBATCH --mail-user   dcha704@aucklanduni.ac.nz
#SBATCH --mail-type   END

module load R/4.1.0-gimkl-2020a
module load admb

srun R --slave < 02B-fascrNew-Runner.R
