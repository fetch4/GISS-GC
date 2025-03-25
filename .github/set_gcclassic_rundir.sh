#!/bin/bash
# ============================================================================ #
# Script for setting the GCClassic rundir in the CI                            #
# ============================================================================ #

set -eu

mkdir -p "${GCCLASSIC_RUNDIR}"
cd "${GISS_HOME}/model/geos-chem/src/GEOS-Chem/run/GCClassic"
./createRunDir.sh <<<"1
1
1
2
2
${GCCLASSIC_RUNDIR}

n"
