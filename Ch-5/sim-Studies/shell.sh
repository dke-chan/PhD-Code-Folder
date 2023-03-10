module load R/4.1.0-gimkl-2020a

export nSim=1000
export nCore=250

mkdir -p ./results ./consoleOutputs

R --slave --no-save < ./seedGen.r --args ${nSim} ${nCore}
R --slave --no-save < ./libraryBuilder.r

# sbatch --array=1-${nCore} slurm.sl
