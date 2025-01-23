#!/bin/env bash
# ============================================================================ #
# Run GISS Model E with (or without) GEOS-Chem support.                        #
# ============================================================================ #

set -e

# Default values
NP=1
GISS_ONLY=false
DEBUG=false

# Function to display help text
show_help() {
  echo "Usage: $0 [NP=<integer>] [--giss-only]"
  echo
  echo "Arguments:"
  echo "  NP          Set number of MPI processes (default: 1). Must be an integer."
  echo
  echo "Options:"
  echo "  --help      Show this help message and exit."
  echo "  --giss-only Build without GEOS-Chem coupling."
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

# Set RUNID appropriately
if [ "${GISS_ONLY}" = true ]; then
  RUNID=GISS_ONLY
else
  RUNID=GISS_GC_14
fi
if [ "${DEBUG}" = true ]; then
  ln -s -f $(pwd)/${RUNID}.R $(pwd)/${RUNID}_DEBUG.R
  RUNID="${RUNID}_DEBUG"
fi

# Navigate to the run directory and run the model for one hour
cd ${ModelE_Support}/prod_runs/${RUNID}
./${RUNID}ln
MP_SET_NUM_THREADS="${NP}" ./${RUNID} -i I -cold-restart &
tail -f ${RUNID}.PRT
