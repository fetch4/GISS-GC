#!/bin/env bash
# ============================================================================ #
# Build GISS Model E with (or without) GEOS-Chem support.                      #
# ============================================================================ #

set -e

# Default values
OPENMP=false
GISS_ONLY=false
MECH=fullchem
CLASSIC=false
DEBUG=false
FRESH=false

# Function to display help text
show_help() {
  echo "Usage: $0 [MECH=fullchem|carbon|Hg|custom] [--openmp] [--giss-only] [--classic]"
  echo "          [-f] [--debug]"
  echo
  echo "Options:"
  echo "  MECH=<mechanism>  Set the chemical mechanism (defaults to fullchem)."
  echo "  --openmp          Compile with OpenMP enabled."
  echo "  --giss-only       Build without GEOS-Chem coupling."
  echo "  --classic   Build without GCClassic as the driver, rather than Model E."
  echo "  --debug           Run with debugging turned on."
  echo "  -f                Fresh rebuild of the model."
  echo "  --help            Show this help message and exit."
}

# Check for --help option
if [ "$1" = "--help" ]; then
  show_help
  exit 0
fi

# Parse arguments
for arg in "$@"; do
  case $arg in
  MECH=*)
    MECH="${arg#*=}"
    ;;
  --openmp)
    OPENMP=true
    shift
    ;;
  --giss-only)
    GISS_ONLY=true
    shift
    ;;
  --classic)
    CLASSIC=true
    shift
    ;;
  --debug)
    DEBUG=true
    shift
    ;;
  -f)
    FRESH=true
    shift
    ;;
  *)
    echo "Unknown argument: $arg"
    show_help
    exit 1
    ;;
  esac
done

if [ "${GISS_ONLY}" = true ]; then
  if [ "${CLASSIC}" = true ]; then
    echo "--giss-only and --classic are mutually exclusive"
    exit 1
  fi
  GC=NO
  RUNID=GISS_ONLY
else
  GC=YES
  RUNID=GISS_GC_14
fi

# Print the values for verification
echo "MECH=${MECH}"
echo "OPENMP=${OPENMP}"
echo "GC=${GC}"
echo "GISS_ONLY=${GISS_ONLY}"
echo "CLASSIC=${CLASSIC}"
echo "FRESH=${FRESH}"
echo "DEBUG=${DEBUG}"

# Conditionally fresh rebuild of the model
cd ${GISS_HOME}/decks/
if [ "${FRESH}" = true ]; then
  make clean OVERWRITE=YES
  make clean_all OVERWRITE=YES
fi

if [ "${DEBUG}" = true ]; then
  TYPE=Debug
else
  TYPE=Release
fi

# Compile
if [ "${CLASSIC}" = true ]; then
  # Set up GCClassic rundir
  cd "${GISS_HOME}/model/geos-chem/src/GEOS-Chem/run/GCClassic"
  ./createRunDir.sh <<<"1
  1
  1
  2
  2
  ${GCCLASSIC_RUNDIR}

  n"
  # Build GCClassic
  cd ${GCCLASSIC_RUNDIR}
  BUILD_DIR=build
  if [ "${DEBUG}" = true ]; then
    BUILD_DIR=${BUILD_DIR}_debug
  fi
  if [ "${FRESH}" = true ]; then
    rm -rf ${BUILD_DIR}
  fi
  mkdir -p ${BUILD_DIR}
  cd ${BUILD_DIR}
  cmake ../CodeDir -DRUNDIR=.. -DCMAKE_BUILD_TYPE=${TYPE} -DMECH=${MECH}
  make -j10
  RUNDIR=${GCCLASSIC_RUNDIR}
else
  if [ "${DEBUG}" = true ]; then
    ln -s -f ${GISS_HOME}/.github/rundecks/${RUNID}.R $(pwd)/${RUNID}_DEBUG.R
    RUNID="${RUNID}_DEBUG"
    make -j setup RUN=${RUNID} F90=mpif90 GC=${GC} MP=${OPENMP} MPI=YES MECH=${MECH} \
      TYPE=${TYPE} DEBUG=YES COMPILE_WITH_TRAPS=YES TRACEBACK=YES OVERWRITE=YES
  else
    ln -s -f ${GISS_HOME}/.github/rundecks/${RUNID}.R $(pwd)/${RUNID}.R
    make -j setup RUN=${RUNID} F90=mpif90 GC=${GC} MP=${OPENMP} MPI=YES MECH=${MECH} \
      TYPE=${TYPE} OVERWRITE=YES
  fi
  RUNDIR=${ModelE_Support}/huge_space/${RUNID}
fi

# Configuration
if [ "${GISS_ONLY}" = false ]; then
  # Copy over configuration files
  CONFIG=${GISS_HOME}/.dev/config
  for DIR in ${GISS_HOME} ${RUNDIR}; do
    ln -s -f ${CONFIG}/geoschem_config.yml ${DIR}/geoschem_config.yml
    ln -s -f ${CONFIG}/HEMCO_Config.rc ${DIR}/HEMCO_Config.rc
    ln -s -f ${CONFIG}/HEMCO_Diagn.rc ${DIR}/HEMCO_Diagn.rc
    ln -s -f ${CONFIG}/HISTORY.rc ${DIR}/HISTORY.rc
    ln -s -f ${CONFIG}/species_database.yml ${DIR}/species_database.yml
  done
  # Create output directories
  mkdir -p ${RUNDIR}/OutputDir
  # Setup restarts
  mkdir -p ${HUGE_SPACE}/Restarts
  # NOTE: The restart file will need to have been saved in the following location
  ln -s -f ${GC_INPUTS}/ExtData/GEOSCHEM_RESTARTS/GC_14.3.0/GEOSChem.Restart.20160701_0000z.LATEST.nc4 \
    ${HUGE_SPACE}/Restarts/GEOSChem.Restart.20160701_0000z.nc4
  # Edit HEMCO_Config to say whether we are running with or without meteorology
  if [ "${CLASSIC}" = true ]; then
    sed -i "s/METEOROLOGY            :       false/METEOROLOGY            :       true /" ${CONFIG}/HEMCO_Config.rc
  else
    sed -i "s/METEOROLOGY            :       true /METEOROLOGY            :       false/" ${CONFIG}/HEMCO_Config.rc
  fi
fi
