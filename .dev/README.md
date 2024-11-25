# Developer utils

**CODE IN THIS SUBDIRECTORY IS NOT INTENDED FOR GENERAL USERS**

First, ensure that any existing `setup.sh`, `build.sh`, or `run.sh` scripts are backed up / moved to another location.
Alternatively, skip any step that isn't needed.

Starting from the location where you cloned GISS-GC (i.e., `${GISS_HOME}`):
```
cp .dev/utils/setup.sh .
source setup.sh

cp .dev/utils/build.sh .
./build.sh --help
# This should show some options. Build with the options you want.

cp .dev/utils/run.sh .
./run.sh --help
# This should show some options. Run with the options you want.
```
