#!/bin/bash -e

#SBATCH --job-name    ascrDisperse
#SBATCH --time        1:00:00
#SBATCH --mem         1000MB
#SBATCH --partition   large
#SBATCH --output      consoleOutputs/ascrDisperse_%j.out
#SBATCH --mail-user   dcha704@aucklanduni.ac.nz
#SBATCH --mail-type   END

module load R/4.1.0-gimkl-2020a
module load admb

srun R --slave < scripts/2x2_Multi_A.r
