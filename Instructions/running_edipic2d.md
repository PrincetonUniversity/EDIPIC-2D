# Running EDIPIC-2D

The easiest way to get started with EDIPIC-2D is by copying one of the examples found in
the top directory of the repository and modifying it to match the parameters of your system.
The example directories are:

  - [input_data_sample](../input_data_sample)
  - [input_data_sample_periodic_SEE](../input_data_sample_periodic_SEE)
  - [input_data_sample_periodic_ebeam](../input_data_sample_periodic_ebeam)
  - [iinput_data_files_waveform](../input_data_files_waveform)

A description of the input parameters can be found in the **Doc** directory:

  - [EDIPIC2D_input_data_description](../Doc/EDIPIC2D_input_data_description_1.pdf)

Most high performance computing (HPC) clusters are equipped with a fast parallel filesystem
where simulations are to be carried out. It is recommended that you copy the example directory
to that location (often called **scratch**) as well as the executable **edipic2d**.

```
cd /gpfs/scratch/my_project   # for example
cp -r ~/EDIPIC-2D/input_data_sample_periodic_ebeam EXAMPLE_EBEAM
cd EXAMPLE_EBEAM
cp ~/EDIPIC-2D/src/edipic2d .     # Copy executable
```

If you are running on the PPPL or Princeton University clusters, you will need to use the
[SLURM](https://slurm.schedmd.com/documentation.html) workload managing system (or "batch"
system) to submit a job. Chances are that your cluster is also running SLURM, in which
case you would use the same type of job script. EDIPIC-2D is an MPI code so it typically
runs on 32 CPU cores (32 MPI tasks). Here is an example job script to launch a 32-task
simulation:

```
#!/bin/bash
## this is sample SLURM script
## Ask for 2 nodes on the cluster
#SBATCH --nodes=2
## Ask for a total of 32 MPI tasks
#SBATCH --ntasks=32
## Ste a run time limit of 4 hours
#SBATCH --time=4:00:00
# sends mail when process begins, and when it ends. Make sure you define your email 
#SBATCH --mail-type=begin 
#SBATCH --mail-type=end 
#SBATCH --mail-user=your@email.aaa 

# Load compiler and openmpi environments
module load intel/2021.1
module load openmpi/intel-2021.1/4.1.0

# Set PETSC_DIR environment variable
export PETSC_DIR=location_of_PETSc_installation

# Add $PETSC_DIR/lib to shared libraries search path
export LD_LIBRARY_PATH=${PETSC_DIR}/lib:${LD_LIBRARY_PATH}

# Run edipic2d
srun -n 32 ./edipic2d >& output.log
```

Here is a job script to run with 96 MPI taks on stellar.princeton.edu: [stellar slurm_script](./slurm_script)

To submit this job to the queuing system, copy the slurm script to your run directory and do:

```
sbatch slurm_script
```

\*\*NOTE\*\* You always need to set **PETSC_DIR** and **LD_LIBRARY_PATH** in order to run the code successfully. Also make sure that you have the file **petsc.rc**, which contains run parameters
for PETSc. There should be a copy in all the example directories.

The examples in the repository take less than 15 minutes to run so you can reduce the time
limit to `--time=15:00` or run interactively.

To run interactively on a system that is not using SLURM, try:

```
mpirun -np 32 ./edipic2d >& output.log &
tail -f output.log    # to see the content of output.log as data is being generated
```

