#!/bin/env bash
# ============================================================================ #
# Build GISS Model E with (or without) GEOS-Chem support.                      #
# ============================================================================ #

set -e

# Default values
FRESH=false
OPENMP=false
GISS_ONLY=false
MECH=fullchem
DEBUG=false

# Function to display help text
show_help() {
  echo "Usage: $0 [MECH=fullchem|carbon|Hg|custom] [--openmp] [--giss-only] [-f] [--debug]"
  echo
  echo "Options:"
  echo "  MECH=<mechanism>  Set the chemical mechanism (defaults to fullchem)."
  echo "  --openmp          Compile with OpenMP enabled."
  echo "  --giss-only       Build without GEOS-Chem coupling."
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
  -f)
    FRESH=true
    shift
    ;;
  --giss-only)
    GISS_ONLY=true
    shift
    ;;
  --debug)
    DEBUG=true
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
echo "FRESH=${FRESH}"
echo "DEBUG=${DEBUG}"

# Conditionally fresh rebuild of the model
cd ${GISS_HOME}/decks/
if [ "${FRESH}" = true ]; then
  make clean OVERWRITE=YES
  make clean_all OVERWRITE=YES
fi

# Compile
if [ "${DEBUG}" = true ]; then
  ln -s -f ${GISS_HOME}/.github/rundecks/${RUNID}.R $(pwd)/${RUNID}_DEBUG.R
  RUNID="${RUNID}_DEBUG"
  make -j setup RUN=${RUNID} F90=mpif90 GC=${GC} MP=${OPENMP} MPI=YES MECH=${MECH} \
    TYPE=Debug DEBUG=YES COMPILE_WITH_TRAPS=YES TRACEBACK=YES OVERWRITE=YES
else
  ln -s -f ${GISS_HOME}/.github/rundecks/${RUNID}.R $(pwd)/${RUNID}.R
  make -j setup RUN=${RUNID} F90=mpif90 GC=${GC} MP=${OPENMP} MPI=YES MECH=${MECH} \
    TYPE=Release OVERWRITE=YES
fi

if [ "${GISS_ONLY}" = false ]; then
  # Copy over configuration files
  HUGE_SPACE=${ModelE_Support}/huge_space/${RUNID}
  CONFIG=${GISS_HOME}/.dev/config
  for DIR in ${GISS_HOME} ${HUGE_SPACE}; do
    ln -s -f ${CONFIG}/geoschem_config.yml ${DIR}/geoschem_config.yml
    ln -s -f ${CONFIG}/HEMCO_Config.rc ${DIR}/HEMCO_Config.rc
    ln -s -f ${CONFIG}/HEMCO_Diagn.rc ${DIR}/HEMCO_Diagn.rc
    ln -s -f ${CONFIG}/HISTORY.rc ${DIR}/HISTORY.rc
    ln -s -f ${CONFIG}/species_database.yml ${DIR}/species_database.yml
  done
  # Create output directories
  mkdir -p ${HUGE_SPACE}/OutputDir
  # Setup restarts
  mkdir -p ${HUGE_SPACE}/Restarts
  # NOTE: The restart file will need to have been saved in the following location
  ln -s -f ${GC_INPUTS}/ExtData/GEOSCHEM_RESTARTS/GC_14.3.0/GEOSChem.Restart.20160701_0000z.LATEST.nc4 \
    ${HUGE_SPACE}/Restarts/GEOSChem.Restart.20160701_0000z.nc4
fi
