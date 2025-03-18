#!/usr/bin/env python3
# """
# Python script for checking the consistency of HISTORY.rc files for GEOS-Chem and GISS
# Model E.
#
# Simply execute as a script:
#   ./check_consistent.py
# """

# Define the two models 'gc' and 'giss' and their corresponding configuration files
filenames = {"gc": "HISTORY.rc", "giss": "HISTORY_ModelE.rc"}
models = tuple(filenames.keys())


def other_model(model):
    """
    Given model 'gc' or 'giss', return the other one.

    :arg model: model name
    :type model: str
    """
    return "gc" if model == "giss" else "giss"


# Read the two configuration files and report the number of lines in each
lines = {}
for model, filename in filenames.items():
    with open(filename, "r") as f:
        lines[model] = f.readlines()
    print(f"Number of lines in {filename:20s}: {len(lines[model])}")

# Check that the COLLECTIONS used in each file are consistent
collections = {}
for model in models:
    found = False
    collections[model] = set()
    for line in lines[model]:
        if line.startswith("COLLECTIONS:"):
            found = True
            line = line.replace("COLLECTIONS:", "")
        if not found:
            continue
        line = line.replace(" ", "").replace("\n", "").replace("'", "").replace(",", "")
        if line.startswith("::"):
            break
        if line.startswith("#"):
            continue
        collections[model].add(line)
for model in models:
    other = other_model(model)
    if not collections[model].issubset(collections[other]):
        raise ValueError(
            f"{filenames[model]} contains collections not in {filenames[other]}:"
            f" {collections[model].difference(collections[other])}"
        )

# TODO: Further checks
