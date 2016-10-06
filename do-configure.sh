#!/bin/bash

# Youu can invoke this shell script with additional command-line
# arguments.  They will be passed directly to CMake.
#
EXTRA_ARGS=$@

#
# Each invocation of CMake caches the values of build options in a
# CMakeCache.txt file.  If you run CMake again without deleting the
# CMakeCache.txt file, CMake won't notice any build options that have
# changed, because it found their original values in the cache file.
# Deleting the CMakeCache.txt file before invoking CMake will insure
# that CMake learns about any build options you may have changed.
#
rm -rf CMakeCache.txt CMakeFiles

#
# A sampling of CMake build options.
#
# CMAKE_INSTALL_PREFIX: Where to install Trilinos.
# CMAKE_BUILD_TYPE: DEBUG or RELEASE.
# CMAKE_CXX_COMPILER: The C++ compiler to use when building Trilinos.
# CMAKE_C_COMPILER: The C compiler to use when building Trilinos.
#   Some parts of Trilinos are implemented in C.
# CMAKE_Fortran_COMPILER: The Fortran compiler to use when building
#   Trilinos.  Some parts of Trilinos are implemented in Fortran.
# HAVE_GCC_ABI_DEMANGLE: Setting this option to ON improves 
#   debugging messages.
# CMAKE_VERBOSE_MAKEFILE: Set to OFF (or FALSE) if you prefer a quiet
#   build.
# Trilinos_ENABLE_ALL_PACKAGES: If you like, you can build _all_ of
#   Trilinos, but you don't have to.
# Trilinos_ENABLE_Epetra: If ON, build Epetra.
# Trilinos_ENABLE_Triutils: If ON, bulid Triutils.
# Trilinos_ENABLE_TESTS: If ON, build the tests for all packages that
#   are to be built.
# Trilinos_ENABLE_EXAMPLES: If ON, build the examples for all
#   packages that are to be built.
#
# Added -g flag to enable debugging/profiling
#
# Optional MKL BLAS / LAPACK setting   
#    -D TPL_ENABLE_MKL:BOOL=ON \
#    -D MKL_LIBRARY_DIRS:STRING="/opt/intel/Compiler/16.0/current/mkl/lib/intel64" \
#    -D BLAS_LIBRARY_DIRS:STRING="/opt/intel/Compiler/16.0/current/mkl/lib/intel64;/opt/intel/Compiler/16.0/current/compiler/lib/intel64;/usr/lib64" \
#    -D BLAS_LIBRARY_NAMES:STRING="mkl_intel_lp64; mkl_intel_thread; mkl_core; iomp5; pthread" \
#    -D LAPACK_LIBRARY_DIRS:STRING="/opt/intel/Compiler/16.0/current/mkl/lib/intel64;/opt/intel/Compiler/16.0/current/compiler/lib/intel64;/usr/lib64" \
#    -D LAPACK_LIBRARY_NAMES:STRING="mkl_intel_lp64; mkl_intel_thread; mkl_core; iomp5; pthread" \
#   
# export PATH=/opt/gcc/gcc-5.2.0/bin/:$PATH
# export LD_LIBRARY_PATH=/opt/gcc/gcc-5.2.0/lib64:$LD_LIBRARY_PATH


PREFIX=~/Documents/TeaLeaf_Trilinos/libs/trilinos

cmake \
    -D CMAKE_INSTALL_PREFIX:PATH=$PREFIX \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D CMAKE_CXX_COMPILER:FILEPATH=mpicxx \
    -D CMAKE_C_COMPILER:FILEPATH=mpicc \
    -D CMAKE_Fortran_COMPILER:FILEPATH=mpif90  \
    -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
    -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=FALSE \
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
    -D Trilinos_ENABLE_Tpetra:BOOL=ON \
    -D Trilinos_ENABLE_Belos:BOOL=ON \
    -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
    -D Trilinos_ENABLE_MueLu:BOOL=ON \
    -D Trilinos_ENABLE_ML:BOOL=ON \
    -D Trilinos_ENABLE_Epetra:BOOL=ON \
    -D Trilinos_ENABLE_Xpetra:BOOL=ON \
    -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
    -D Trilinos_ENABLE_TESTS:BOOL=OFF \
    -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D Trilinos_ENABLE_OpenMP:BOOL=ON \
    -D Trilinos_ENABLE_TEUCHOS_TIME_MONITOR:BOOL=ON \
    -D BLAS_LIBRARY_DIRS:FILEPATH=/usr/lib/lapack/ \
    -D LAPACK_LIBRARY_DIRS:FILEPATH=/usr/lib/libblas/ \
    -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D CMAKE_CXX_FLAGS:STRING="-DMPICH_IGNORE_CXX_SEEK -g" \
    $EXTRA_ARGS \
    ../

make
make install
