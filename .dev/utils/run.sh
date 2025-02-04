#!/bin/env bash
# ============================================================================ #
# Run GISS Model E with (or without) GEOS-Chem support.                        #
# ============================================================================ #

set -e

# Default values
NP=1
GISS_ONLY=false
CLASSIC=false
COLD_RESTART=false
DEBUG=false

# Function to display help text
show_help() {
  echo "Usage: $0 [NP=<integer>] [--giss-only] [--classic] [--cold-restart] [--debug] [--help]"
  echo
  echo "Arguments:"
  echo "  NP          Set number of MPI processes (default: 1). Must be an integer."
  echo
  echo "Options:"
  echo "  --giss-only     Build without GEOS-Chem coupling."
  echo "  --classic       Build without GCClassic as the driver, rather than Model E."
  echo "  --cold-restart  Run for a single hour from the restart files to generate a checkpoint."
  echo "  --debug         Run with debugging turned on."
  echo "  --help          Show this help message and exit."
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
  --classic)
    CLASSIC=true
    shift
    ;;
  --cold-restart)
    COLD_RESTART=true
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

# Print the values for verification
echo "GISS_ONLY=${GISS_ONLY}"
echo "CLASSIC=${CLASSIC}"
echo "COLD_RESTART=${COLD_RESTART}"
echo "DEBUG=${DEBUG}"
echo "NP=${NP}"

# Check for unset environment variables
if [ -z ${ModelE_Support+x} ]; then
  echo "ModelE_Support is unset. Exiting."
  exit 0
fi
if [ "${CLASSIC}" = true ]; then
  if [ -z ${GCCLASSIC_RUNDIR+x} ]; then
    echo "GCCLASSIC_RUNDIR is unset. Exiting."
    exit 0
  fi
fi

if [ "${CLASSIC}" = true ]; then
  if [ "${NP}" != "1" ]; then
    echo "GCClassic only runs in serial"
    exit 1
  fi
  if [ "${COLD_RESTART}" = true ]; then
    echo "GCClassic does not support cold restart"
    exit 1
  fi
  cd "${GCCLASSIC_RUNDIR}"
  if [ "${DEBUG}" = true ]; then
    ./build_debug/bin/gcclassic
  else
    ./build/bin/gcclassic
  fi
else
  # Set RUNID appropriately
  if [ "${GISS_ONLY}" = true ]; then
    RUNID=GISS_ONLY
  else
    RUNID=GISS_GC_14
  fi
  if [ "${DEBUG}" = true ]; then
    ln -s -f "$(pwd)/${RUNID}.R" "$(pwd)/${RUNID}_DEBUG.R"
    RUNID="${RUNID}_DEBUG"
  fi
  echo "RUNID=${RUNID}"

  # Navigate to the run directory
  cd "${ModelE_Support}/prod_runs/${RUNID}"
  ./${RUNID}ln
  if [ "${COLD_RESTART}" = true ]; then
    # Run the model for one hour
    ./${RUNID} -np "${NP}" -i I -cold-restart -l cold-restart.log
  else
    # Pick up from a checkpoint and run the model for the full duration
    ./${RUNID} -np "${NP}" -i I
  fi
fi
