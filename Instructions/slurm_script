#!/bin/bash
## this is sample SLURM script for Stellar (stellar.princeton.edu)
## Stellar has 96 physical cores per node
#SBATCH --nodes=1
## Ask for a total of 96 MPI tasks
#SBATCH --ntasks=96
## Set a run time limit of 4 hours
#SBATCH --time=4:00:00
## sends mail when process begins, and when it ends. Make sure you define your email
#SBATCH --mail-type=begin 
#SBATCH --mail-type=end 
#SBATCH --mail-user=your@email.aaa 

# Load compiler and MPI environments
# On the Princeton University systems, the modules need to be spelled out
# We use the Intel-MPI library on Stellar
module load intel/2021.1.2
module load intel-mpi/intel/2021.3.1

# Set PETSC_DIR environment variable
# This location should work but you can also build your own version of PETSc
export PETSC_DIR=/scratch/gpfs/ethier/STELLAR/Software/INTEL_MPI_2021

# Add $PETSC_DIR/lib to shared libraries search path
export LD_LIBRARY_PATH=${PETSC_DIR}/lib:${LD_LIBRARY_PATH}

# Run edipic2d
srun -n 96 ./edipic2d >& output.log

