#!/usr/bin/bash
# ============================================================================ #
# Build GISS Model E with (or without) GEOS-Chem support.                      #
# ============================================================================ #

set -e

# Default values
FRESH=false
OPENMP=false
GISS_ONLY=false
DEBUG=false

# Function to display help text
show_help() {
  echo "Usage: $0 [MP=YES|NO] [--openmp] [--giss-only] [-f] [--debug]"
  echo
  echo "Options:"
  echo "  --help      Show this help message and exit."
  echo "  --openmp    Compile with OpenMP enabled."
  echo "  --giss-only Build without GEOS-Chem coupling."
  echo "  -f          Fresh rebuild of the model."
  echo "  --debug     Run with debugging turned on."
}

# Check for --help option
if [ "$1" = "--help" ]; then
  show_help
  exit 0
fi

# Parse arguments
for arg in "$@"; do
  case $arg in
  NP=*)
    NP="${arg#*=}"
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
echo "NP=${NP}"
echo "OPENMP=${OPENMP}"
echo "GC=${GC}"
echo "FRESH=${FRESH}"
echo "DEBUG=${DEBUG}"

# Conditionally fresh rebuild of the model
cd ${GISS_HOME}/decks/
ln -s -f ${GISS_HOME}/.github/rundecks/${RUNID}.R $(pwd)/${RUNID}.R
if [ "${FRESH}" = true ]; then
  make clean
  make clean_all
fi

# Compile
if [ "${DEBUG}" = true ]; then
  ln -s -f $(pwd)/${RUNID}.R $(pwd)/${RUNID}_DEBUG.R
  RUNID="${RUNID}_DEBUG"
  git apply ${GISS_HOME}/.dev/utils/DEBUGGING_FLAGS.patch
  make -j setup RUN=${RUNID} F90=mpif90 GC=${GC} MP=${OPENMP} MPI=YES MECH=carbon \
    TYPE=Debug DEBUG=YES COMPILE_WITH_TRAPS=YES TRACEBACK=YES OVERWRITE=YES
  git apply -R ${GISS_HOME}/.dev/utils/DEBUGGING_FLAGS.patch
else
  make -j setup RUN=${RUNID} F90=mpif90 GC=${GC} MP=${OPENMP} MPI=YES MECH=carbon \
    TYPE=Release
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
  ln -s -f ${GC_INPUTS}/ExtData/GEOSCHEM_RESTARTS/GC_14.3.0/GEOSChem.Restart.fullchem.20190701_0000z.nc4 \
    ${HUGE_SPACE}/Restarts/GEOSChem.Restart.20190701_0000z.nc4
fi
