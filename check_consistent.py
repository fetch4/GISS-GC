#!/usr/bin/env python3
# """
# Python script for checking the consistency of HISTORY.rc files for GEOS-Chem and GISS
# Model E.
#
# Simply execute as a script:
#   ./check_consistent.py
# """

filenames = {"gc": "HISTORY.rc", "giss": "HISTORY_ModelE.rc"}
lines = {}

for model, filename in filenames.items():
    with open(filename, "r") as f:
        lines[model] = f.readlines()

    print(f"Number of lines in {filename:20s}: {len(lines[model])}")

# TODO: Actually check the files are consistent
