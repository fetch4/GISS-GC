#!/usr/bin/bash
# ============================================================================ #
# Activate the Python and spack environments used by GISS-GC.                  #
# ============================================================================ #

# Envronment variables for GISS modelE
# NOTE: Path may need to be edited for your system
export GISS_HOME=${SOFTWARE}/GISS-GC
export ModelE_Support=${HOME}/run/giss-gc
mkdir -p ${ModelE_Support}
# Environment variables for compiler
export CC=gcc
export CXX=g++
export FC=gfortran
export F90=${FC}
export F77=${FC}
# Misc. enviroment variables
export F_UFMTENDIAN=big
export KMP_STACKSIZE=100000000
export OMP_NUM_THREADS=1

# Spack setup
spack env activate -p giss-gc
export MPI_ROOT=${HOME}/software/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.4.0/openmpi-4.1.6-s3fu5gvaasgjy4jecnb6rvemx7oofexx

# Environment variables for passing NetCDF-C paths to GEOS-Chem
export NETCDF_HOME=$(nc-config --prefix)
export GC_BIN=${NETCDF_HOME}/bin
export GC_INCLUDE=${NETCDF_HOME}/include
export GC_LIB=${NETCDF_HOME}/lib

# Environment variables for passing NetCDF-Fortran paths to GEOS-Chem
export NETCDF_F_HOME=$(nf-config --prefix)
export GC_F_BIN=${NETCDF_F_HOME}/bin
export GC_F_INCLUDE=${NETCDF_F_HOME}/include
export GC_F_LIB=${NETCDF_F_HOME}/lib

# GEOS-Chem input data
export ROOT=${DATA}/gcclassic/ExtData/HEMCO/

# Put tools in the path
export PATH=${SOFTWARE}/tools/mk_diags:${PATH}
