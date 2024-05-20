# Installing PETSc and HYPRE

The HYPRE library can be installed as part of the PETSc installation, which is the recommended approach. It is also the case for the BLAS/LAPACK libraries although HYPRE installs its own versions.

The requirements for installing PETSc+HYPRE are:
 - C/C++ and Fortran compilers
 - A working version of MPI built with the compilers above (use `mpicc -show` to verify)


## Downloading PETSc

EDIPIC-2D currently uses version 3.14 of PETSc so it is recommended to download and installed this version. There are two options for downloading:

### Option 1: Clone GIT repository

```
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc
git checkout v3.14.6
```

### Option 2: Gzipped tar file

The easiest way to get the tar file is by using `wget` if it is available on your system (`curl -O` on MacOS)

```
wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.14.6.tar.gz
tar -zxvf petsc-3.14.6.tar.gz
cd petsc-3.14.6
```

## Configuring PETSc

The first thing to do before running `configure` is to decide where the final libraries and include files will be installed. Ideally you want this location to be accessible by the whole group of users who will be running EDIPIC-2D. This location/directory is the value of `--prefix` below. 
Run the following configure command within the petsc distribution top directory after making sure to replace `<path_to_final_installation>` with a valid path location and `<compiler_name_and_version>` with a judicious name (for example `gcc12_mpich`). It is useful to put it in a shell script:

```
   ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            COPTFLAGS='-O2' CXXOPTFLAGS='-O2' FOPTFLAGS='-O2' \
            --with-blaslapack-dir=$MKLROOT \
            --download-hypre \
            --prefix=<path_to_final_installation>  \
            PETSC_ARCH=<compiler_name_and_version>  \
            PETSC_DIR=$PWD
```


## Compiling

Once the **configure** step is successful, run the following **make** commands:

```
   make PETSC_DIR=$PWD PETSC_ARCH=compiler_name_and_version all
   make PETSC_DIR=$PWD PETSC_ARCH=compiler_name_and_version install
```


## Example: Building PETSc+HYPRE on the PPPL cluster

Users of the PPPL cluster manage their environment with the `module` tool.
To load the default Intel compiler and OpenMPI library, one simply uses:

```
   module purge            # clear environment
   module load intel       # load default Intel development tools and libraries
   module load openmpi     # load OpenMPI library built with Intel compilers

   cd software
   git clone -b release https://gitlab.com/petsc/petsc.git petsc
   cd petsc
   git checkout v3.14.6    # get version 3.14.6

   ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            COPTFLAGS='-O2' CXXOPTFLAGS='-O2' FOPTFLAGS='-O2' \
            --with-blaslapack-dir=$MKLROOT \
            --download-hypre \
            --prefix=${HOME}/software  \
            PETSC_ARCH=INTEL_2019  \
            PETSC_DIR=${PWD}

   make PETSC_DIR=${PWD} PETSC_ARCH=INTEL_2019 all
   make PETSC_DIR=${PWD} PETSC_ARCH=INTEL_2019 install
```

Notice that we added the line `--with-blaslapack-dir=$MKLROOT` to the configure
step. The Intel development software comes with the highly optimized
**Math Kernel Library (MKL)**, which includes highly optimized versions of
the **BLAS** and **LAPACK** libraries. These library functions are used by
PETSc and HYPRE and the lowest level so it is important to use the fastest
version possible. The **MKLROOT** environment variable contains the location
of the MKL library on the system. It is set when loading the `intel` module.


## Another example: Building PETSc+HYPRE on STELLAR at Princeton University

Users of the Princeton University clusters also manage their environment with the `module` tool.
On STELLAR we use the Intel compiler and Intel-MPI library:

```
   module purge                            # clear environment
   module load intel/2021.1.2              # load default Intel development tools and libraries
   module load intel-mpi/intel/2021.3.1    # load Intel-MPI library

   cd software
   git clone -b release https://gitlab.com/petsc/petsc.git petsc
   cd petsc
   git checkout v3.14.6    # get version 3.14.6

   ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            COPTFLAGS='-O2' CXXOPTFLAGS='-O2' FOPTFLAGS='-O2' \
            --with-blaslapack-dir=$MKLROOT \
            --download-hypre \
            --prefix=${HOME}/software  \
            PETSC_ARCH=INTEL_MPI_2021  \
            PETSC_DIR=${PWD}

   make PETSC_DIR=${PWD} PETSC_ARCH=INTEL_MPI_2021 all
   make PETSC_DIR=${PWD} PETSC_ARCH=INTEL_MPI_2021 install
```

Notice that we added the line `--with-blaslapack-dir=$MKLROOT` to the configure
step. The Intel development software comes with the highly optimized
**Math Kernel Library (MKL)**, which includes highly optimized versions of
the **BLAS** and **LAPACK** libraries. These library functions are used by
PETSc and HYPRE and the lowest level so it is important to use the fastest
version possible. The **MKLROOT** environment variable contains the location
of the MKL library on the system. It is set when loading the `intel` module.

