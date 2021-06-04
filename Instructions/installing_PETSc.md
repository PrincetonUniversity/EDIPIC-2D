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
wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.14.tar.gz
tar -zxvf petsc-3.14.tar.gz
cd petsc-3.14
```

## Configuring PETsc

The first thing to do before running `configure` is to decide where the final libraries and include files will be installed. Ideally you want this location to be accessible by the whole group of people who will be running EDIPIC-2D. This location/directory is the value of `--prefix` below. 
Run the following configure command within the petsc distribution top directory. It is useful to put it in a shell script:

```
   ./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
            COPTFLAGS='-O2' CXXOPTFLAGS='-O2' FOPTFLAGS='-O2' \
            --with-blaslapack-dir=$MKLROOT \
            --download-hypre \
            --prefix=*path_to_final_installation*  \
            PETSC_ARCH=*compiler_name_and_version*  \
            PETSC_DIR=$PWD
```


## Compiling

Once the **configure** step is successful, run the following **make** commands:

```
   make PETSC_DIR=$PWD PETSC_ARCH=compiler_name_and_version all
   make PETSC_DIR=$PWD PETSC_ARCH=compiler_name_and_version install
```

