## Building GISS-GC on a laptop with a Linux operating system

#### Step 1: Basic setup

> [!INFO] Do this the first time only.

Choose the location where GISS-GC will exist (e.g., `~/software/GISS-GC`) and use it to define the `${GISS_HOME}` environment variable:
```sh
export GISS_HOME=${HOME}/software/GISS-GC
```

Checkout the GISS-GC repository into the `${GISS_HOME}` location and set up its submodules:
```sh
git clone git@github.com:fetch4/GISS-GC.git ${GISS_HOME}
cd ${GISS_HOME}
git checkout develop
git submodule init
git submodule update
```

In the directory, create a script `setup.sh` for defining modules to be loaded and various paths required for the GISS-GC build. Copy and paste the following code block into the script, noting that a few entries will require modification:
```sh
#!/bin/bash

# ============================================================================ #
# Activate the Python and Spack environments used by GISS-GC.                  #
# ============================================================================ #

# Environment variables for GISS modelE
export SOFTWARE=${HOME}/software          # NOTE: Edit if you're using a different location for software
export GISS_HOME=${SOFTWARE}/GISS-GC      # NOTE: This should be consistent with the location you cloned GISS-GC to
export ModelE_Support=${HOME}/run/giss-gc # NOTE: Edit if you're using a different run directory location
mkdir -p ${ModelE_Support}/exec
mkdir -p ${ModelE_Support}/huge_space
mkdir -p ${ModelE_Support}/prod_decks
mkdir -p ${ModelE_Support}/prod_input_files
mkdir -p ${ModelE_Support}/prod_runs
# Environment variables for compiler
export CC=gcc       # NOTE: Edit if you're using a different C compiler
export CXX=g++      # NOTE: Edit if you're using a different C++ compiler
export FC=gfortran  # NOTE: Edit if you're using a different Fortran compiler
export F90=gfortran # NOTE: Edit if you're using a different Fortran compiler
export F77=gfortran # NOTE: Edit if you're using a different Fortran compiler
# Misc. environment variables
export F_UFMTENDIAN=big
export KMP_STACKSIZE=100000000
export OMP_NUM_THREADS=36

# Spack environment # NOTE: Edit MPI_ROOT to match the path to your Spack-installed MPI distribution
spack env activate -p giss-gc
export MPI_ROOT=${SOFTWARE}/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.4.0/openmpi-4.1.6-s3fu5gvaasgjy4jecnb6rvemx7oofexx

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

# GEOS-Chem input data # NOTE: Edit if you are storing this somewhere else
export ROOT=${ModelE_Support}/prod_input_files/ExtData/HEMCO/

# Add modelE tools to path
export PATH=${HOME}/software/tools/mk_diags:${PATH}
```

> [!WARNING] Don't `source` the script yet because we need to set up Spack first.


#### Step 2: Setup Spack environment for GISS-GC

> [!INFO] Do this the first time only.

Install the Spack package manager using the [instructions](https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html) online, if you don't already have it. I put it in `~/software/spack`.

Create a Spack configuration in `${GISS_HOME}/spack.yaml` and copy and paste the
following code block into it. Note that you may need to adjust the compiler,
compiler version, and operating system.
```yaml
spack:
  packages:
    all:
      compiler: [gcc@11.4.0]
  container:
    images:
      os: ubuntu:22.04
      spack: '0.21'
  specs:
  - cmake
  - gmake
  - hdf5
  - netcdf-c+mpi+parallel-netcdf
  - netcdf-fortran
  - openmpi
  concretizer:
    unify: true
  config:
    install_missing_compilers: true
  compilers:
  - compiler:
      spec: gcc@=11.4.0
      paths:
        cc: /bin/gcc
        cxx: /bin/g++
        f77: /bin/gfortran
        fc: /bin/gfortran
      flags: {}
      operating_system: ubuntu22.04
      target: x86_64
      modules: []
      environment: {}
      extra_rpaths: []
```
Note that this will create a file `~/.spack/linux/compilers.yaml` containing the
compiler specification.

Create a Spack environment off the YAML file with
```sh
cd ${GISS_HOME}
spack env create giss-gc spack.yaml
```

With the Spack environment active, install the packages required by the Spack
environment with
```sh
spack install
```
This may take an hour or so.

Run
```sh
source ${GISS_HOME}/setup.sh
```
and check that `${MPI_ROOT}`, `${NETCDF_HOME}`, and `${NETCDF_F_HOME}` (for example) are valid paths. The activation of the Spack environment should result in your command line prompt being prefaced with `[giss-gc]`, possibly in a different colour.


#### Step 3: Build GISS tools and put them in the path

> [!INFO] Do this the first time only.

GISS comes with some tools (e.g., for post-processing diagnostics) which need to be compiled.
We already added this location to the path in `setup.sh`.
```sh
# source ${GISS_HOME}/setup.sh
cp -rp ${GISS_HOME}/model/mk_diags ${HOME}/software/tools
cd ${HOME}/software/tools/mk_diags
/bin/bash compscr
```

#### Step 4: Create `~/.modelErc`

> [!INFO] Do this the first time only.

GISS is configured using a global `~/.modelErc` file. This is created as follows:
```sh
# source ${GISS_HOME}/setup.sh
cd ${GISS_HOME}/decks
make config COMPILER=${FC} ModelE_Support=${ModelE_Support}
```

Tweak `~/.modelErc` as follows:
```diff
# This file contains global options for modelE. 
# By default it assumes that the directory structure for modelE runs
# is set under /scratch/jwallwo2/run .

## Directory structure ##

# DECKS_REPOSITORY - a directory for permanenet storage of run info.
# All rundecks that you create will be copied to this directory. 
-DECKS_REPOSITORY=/scratch/jwallwo2/run/prod_decks
+DECKS_REPOSITORY=${ModelE_Support}/prod_decks

# CMRUNDIR - directory to which all run directories will be linked.
# This directory will be searched by most scripts for locations of 
# specific runs.
-CMRUNDIR=/scratch/jwallwo2/run/prod_runs
+CMRUNDIR=${ModelE_Support}/prod_runs

# GCMSEARCHPATH - directory to search for gcm input files.
# All necessary input files should be copied or linked to this directory.
-GCMSEARCHPATH=/scratch/jwallwo2/run/prod_input_files
+CMSEARCHPATH=${ModelE_Support}/prod_input_files

# EXECDIR - path to directory with modelE scripts and with some
# executables. This directory should contain the scripts from modelE/exec.
-EXECDIR=/scratch/jwallwo2/run/exec
+EXECDIR=${ModelE_Support}/exec

# SAVEDISK - a directory where all run directories (which will contain
# all output files such as rsf, acc etc.) will be created. This should
# be big enough to accomodate all model output.
-SAVEDISK=/scratch/jwallwo2/run/huge_space
+SAVEDISK=${ModelE_Support}/huge_space

## External libraries ##

# Some of these options can be provided by environment modules (if you 
# use them). Specify here only what is necessary. Options specified 
# here will overwrite options proviided by environment modules.

# NETCDFHOME - path to location of netcdf installation directory.
-# NETCDFHOME=/opt/netcdf/3.6.3
+NETCDFHOME=${NETCDF_F_HOME}

# MPI - set to YES if you want to compile the model for parallel 
# execution on multiple CPU cores. Keep in mind, that functional 
# MPI library should be installed on your computer and its type 
# and location should be specified below.
# This option can be overwritten from the compile line.
-MPI=NO
+MPI=YES

# MPIDISTR - the MPI distribution you are using. Currently supported 
# distributions are: 'intel, 'openmpi', 'mpich2', 'mvapich2', 'SCALI',
# 'mpt' 
-# MPIDISTR=openmpi
+MPIDISTR=openmpi

# MPIDIR - path to the MPI installation directory. (Needs to be set
# only if compiler can't find correct MPI library and include files by
# default)
-# MPIDIR=/opt/openmpi
+MPIDIR=${MPI_ROOT}

# MPILIBDIR - path to the location of MPI library. Set it only if 
# it is different from the default $MPIDIR/lib
# MPILIBDIR=/opt/openmpi/lib

# MPIINCLUDEDIR - path to location of MPI include files. Set it only
# if it is different from the default $MPIDIR/include
# MPIINCLUDEDIR=/opt/openmpi/include

# ESMF5_DIR - path to the installation directory of ESMF (version 5)
# library. (Required only for Cubed Sphere simulations)
# ESMF5_DIR=

# ESMF_BOPT - optimization level of ESMF library. (Should only be used
# togeteher with ESMF5_DIR)
# ESMF_BOPT=O
+ESMF=NO

## Architecture and compiler

# ABI - Application Binary Interfaces. This variable specifies the
# architecture you are using. The valid values are '64' and '32'. 
# On most modern systems you should use '64'. Use '32' if your
# hardware or compiler support only 32-bit binaries.
ABI=64

# COMPILER - specifies the Fortran compiler you are using. Currently
# only 'intel' and 'gfortran' are supported. ('nag' has partial
# support on development branch.) If you are using Modules for
# Environment Management, then this variable may already be set in the
# environment. In this case you don't need to set it here.
-# COMPILER=gfortran
+COMPILER=gfortran

## General User Preferences ##

# MAILTO - email address of the user. When the program ends/crashes
# all notifications will be sent to this address. Leave empty 
# or unset if you don't want to receive these emails
MAILTO=

# UMASK - the value of 'umask' you want to use for model runs. The files
# inside the run directory will have permissions set according to this
# mask.
UMASK=002

# OVERWRITE - can "gmake rundeck" overwrite files already in repository?
# (i.e. in the directory DECKS_REPOSITORY)
OVERWRITE=NO

# OUTPUT_TO_FILES - if set to YES all errors and warnings will be sent
# to files with the names <source_name>.ERR
OUTPUT_TO_FILES=NO

# VERBOSE_OUTPUT - if set to YES gmake will show compilation commands
# and some other information. Otherwise most of the output will be
# suppressed
VERBOSE_OUTPUT=YES
```

#### Step 5: Get the rundeck

> [!INFO] Do this the first time only.

GISS uses 'rundecks' for configuring the model, including setting up preprocessors. These should live within the `${GISS_HOME}/decks` subdirectory. Currently, we are using the `GISS_GC_14.R` rundeck.
This can be found at `.github/rundecks/GISS_GC_14.R`.


#### Step 6: Download the data

> [!INFO] Do this the first time only.

We aren't actually going to be able to run the model on our laptops, so it's sufficient to just create empty files with the expected names:

```sh
export DATA=${ModelE_Support}/prod_input_files
touch ${DATA}/CD144X90.ext.nc
touch ${DATA}/cloud.epsilon4.72x46
touch ${DATA}/CO2profile.Jul16-2017.txt
touch ${DATA}/CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
touch ${DATA}/dH2O_by_CH4_monthly
touch ${DATA}/GHG.CMIP6.1-2014.txt
touch ${DATA}/GIC.144X90.DEC01.1.ext_1.nc
touch ${DATA}/GLMELT_144X90_gas.OCN.nc
touch ${DATA}/H2Ocont_MT_CKD
touch ${DATA}/Irrig144x90_1848to2100_FixedFuture_v3.nc
touch ${DATA}/ISCCP.tautables
touch ${DATA}/LWCorrTables33k
touch ${DATA}/LWTables33k_lowH2O_CO2_O3_planck_1-800
touch ${DATA}/miescatpar.abcdv2
touch ${DATA}/MSU_SSU_RSS_weights.txt
touch ${DATA}/NCARIC.144x90.D7712010_ext.nc
touch ${DATA}/o3_2010_shindell_144x90x49_April1850.nc
touch ${DATA}/oct2003.relhum.nr.Q633G633.table
touch ${DATA}/OST_144x90.1876-1885avg.CMIP6.nc
touch ${DATA}/RD_Fd.nc
touch ${DATA}/RD_Fd.names.txt
touch ${DATA}/REG2X2.5
touch ${DATA}/S144X900098M.ext.nc
touch ${DATA}/sgpgxg.table8
touch ${DATA}/SICE_144x90.1876-1885avg.CMIP6.nc
touch ${DATA}/soil_textures_top30cm_2x2.5
touch ${DATA}/soilcarb_top30cm_2x2.5.nc
touch ${DATA}/solar.CMIP6official.ann1850-2299_with_E3_fastJ.nc
touch ${DATA}/STRATAER.VOL.1850-2014_CMIP6_hdr
touch ${DATA}/top_index_144x90_a.ij.ext.nc
touch ${DATA}/topcld.trscat8
touch ${DATA}/V144x90_EntMM16_height_trimmed_scaled_ext.nc
touch ${DATA}/V144x90_EntMM16_lai_max_trimmed_scaled_ext.nc
touch ${DATA}/V144x90_EntMM16_lai_trimmed_scaled_ext.nc
touch ${DATA}/V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.ext.nc
touch ${DATA}/Z2HX2fromZ1QX1N.BS1.nc
touch ${DATA}/ZSIfac_144x90.1876-1885avg.CMIP6.nc
touch ${DATA}/ZVAR2X25A.nc
mkdir -p ${DATA}/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCA
mkdir -p ${DATA}/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCB
mkdir -p ${DATA}/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCdalbsn
mkdir -p ${DATA}/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/DUST
mkdir -p ${DATA}/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/NIT
mkdir -p ${DATA}/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/O3
mkdir -p ${DATA}/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/OCA
mkdir -p ${DATA}/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/SSA
mkdir -p ${DATA}/cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/SUL
mkdir -p ${DATA}/nudging/merra2
touch ${DATA}/nudging/merra2/uwnd.2014.MERRA2onGISSE2.nc4
touch ${DATA}/nudging/merra2/vwnd.2014.MERRA2onGISSE2.nc4
touch ${DATA}/nudging/merra2/uwnd.2015.MERRA2onGISSE2.nc4
touch ${DATA}/nudging/merra2/vwnd.2015.MERRA2onGISSE2.nc4
touch ${DATA}/nudging/merra2/uwnd.2016.MERRA2onGISSE2.nc4
touch ${DATA}/nudging/merra2/vwnd.2016.MERRA2onGISSE2.nc4
```
#### Step 7: Compile

> [!CHECK] Do this every time.

Create an executable `build.sh` script in `${GISS_HOME}` with the following contents:
```sh
#!/bin/bash

# Default value for MP
MP=NO

# Check if an argument is provided
if [ $# -gt 0 ]; then
  MP=$1
fi

# Check the RUNID environment variable is set
if [ -z "${RUNID}" ]; then
  echo "Error: RUNID environment variable is not set."
  exit 1
fi

cd decks/
# make clean
make -j setup RUN=${RUNID} F90=mpif90 GC=YES MP=${MP} MPI=YES
cd -
```
and run it with `./build.sh`.

Notes:
* Run with `./build.sh MP=YES` to compile with OpenMP shared memory parallelism
* Uncomment the `make clean` line to build from scratch.

#### Step 8: Run the model

> [!DANGER] Do not attempt to run the model on your laptop - it will require far too much RAM.
