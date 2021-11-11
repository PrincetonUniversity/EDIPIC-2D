# Installing EDIPIC-2D

EDIPIC-2D is a Fortran code so you need access to a Fortran compiler. The current makefile
works with the Intel Fortran compiler as well as with GNU gfortran. The makefile determines
which compiler is being used by looking at the output of the "mpifort -show" command
(see MPI section below).

EDIPIC-2D requires a few libraries in order to run. These are :
 - MPI ([OpenMPI](https://www.open-mpi.org/), [MPICH](https://www.mpich.org/), [Intel-MPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html), etc.)
 - [PETSc](https://www.mcs.anl.gov/petsc/)
 - [HYPRE](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods)
 - BLAS/LAPACK

## MPI

MPI stands for "Message Passing Interface" and it is the most widely used distributed parallel communication library in high performance computing. If you are working on a Linux cluster there is probably an MPI library already installed. Make sure that you have access to the commands **mpif90** and **mpicc**. These are wrapper scripts that come with all MPI distributions and they take care of finding the MPI include files and linking with the proper MPI libraries. The **Makefile** included with EDIPIC-2D uses the `mpifort -show` command to determine the Fortran compiler being used. For example ::

```
$ mpif90 -show
gfortran -I/usr/pppl/gcc/9.3-pkgs/openmpi-4.0.3/include -pthread -I/usr/pppl/gcc/9.3-pkgs/openmpi-4.0.3/lib -Wl,-rpath -Wl,/usr/pppl/gcc/9.3-pkgs/openmpi-4.0.3/lib -Wl,--enable-new-dtags -L/usr/pppl/gcc/9.3-pkgs/openmpi-4.0.3/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi
```

This shows that the **gfortran** compiler is being used as well as the open source MPI library **openmpi**. 

If MPI is not installed on your system, you should ask your system administrator to install it since it requires "superuser" privileges. 

## PETSc + HYPRE

The distributed field equations in EDIPIC-2D are solved in parallel using a
multi-grid algorithm from the [HYPRE](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) library as preconditioner,
and GMRES for the final solve. The calls and communication between distributed
tasks are handled by the [PETSc](https://www.mcs.anl.gov/petsc/) framework.
While HPYRE is a separate library, it can be downloaded and built as part of
the PETSc installation. For a quick guide on how to install PETSc+HYPRE, see
the following document:

  [Installing PETSc+HYPRE](./installing_PETSc.md)


## Compiling EDIPIC-2D

After setting the required environment variables, EDIPIC-2D can be compiled by using `make` in the `src` directory ::

```
export PETSC_DIR=path_to_petsc_installation
export LD_LIBRARY_PATH=${PETSC_DIR}/lib:${LD_LIBRARY_PATH}
cd src
make
```

If everything proceeds without error you will get an executable file named `edipic2d`.


## Building EDIPIC-2D on the **STELLAR** cluster at Princeton University

General information about Stellar can be found [HERE](https://researchcomputing.princeton.edu/systems/stellar).

The preferred compiler and MPI library for building EDIPIC-2D on Stellar are the Intel OneAPI versions (as of 11/10/2021, there seems to be a problem with the OpenMPI library but Intel-MPI works fine). They can be accessed by loading the following modules:

```
module purge
module load intel/2021.1.2
module load intel-mpi/intel/2021.3.1
```
Compiled versions of the PETSc and HYPRE libraries can be found in **/scratch/gpfs/ethier/STELLAR/Software/INTEL_MPI_2021**. You just need to set the following environment variables for compiling and running EDIPIC-2D :

```
export PETSC_DIR=/scratch/gpfs/ethier/STELLAR/Software/INTEL_MPI_2021
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/${PETSC_DIR}/lib
```
Alternatively, you can build PETSc and Hypre from source by following these instructions: [Installing PETSc+HYPRE](./installing_PETSc.md)


