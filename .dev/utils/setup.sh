#!/usr/bin/bash
# ============================================================================ #
# Activate the Python and spack environments used by GISS-GC.                  #
# ============================================================================ #

# Environment variables for GISS modelE
# NOTE: Path may need to be edited for your system
export GISS_HOME=${HOME}/software/GISS-GC
# NOTE: Path may need to be edited for your system
export ModelE_Support=${HOME}/run/giss-gc
mkdir -p ${ModelE_Support}
# Environment variables for compiler
export CC=gcc      # NOTE: C compiler may need to be modified for your system
export CXX=g++     # NOTE: C++ compiler may need to be modified for your system
export FC=gfortran # NOTE: Fortran compiler may need to be modified for your system
export F90=${FC}
export F77=${FC}
# Misc. enviroment variables
export F_UFMTENDIAN=big
export KMP_STACKSIZE=100000000
export OMP_NUM_THREADS=1

# Spack setup
# NOTE: This section may need to be edited for your setup
spack env activate -p giss-gc
MPIF90=$(find ${SPACK_ENV} -name mpif90 | head -n 1)
export MPI_ROOT=${MPIF90%/bin/mpif90}

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
# NOTE: Path may need to be edited for your system
export GC_INPUTS=${DATA}/gcclassic
export ROOT=${GC_INPUTS}/ExtData/HEMCO/

# Put tools in the path
export PATH=${SOFTWARE}/tools/mk_diags:${PATH}
