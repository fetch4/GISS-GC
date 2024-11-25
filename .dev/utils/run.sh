#!/usr/bin/sh

# Default values
NP=1
GISS_ONLY=false

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

# Navigate to the run directory and run the model for one hour
cd ${ModelE_Support}/prod_runs/${RUNID}
./${RUNID}ln
mpiexec -np ${NP} ./${RUNID}.exe -i I -cold-restart
