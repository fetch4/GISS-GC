# Developer Zone

**NOTE: CODE IN THIS SUBDIRECTORY IS NOT INTENDED FOR GENERAL USERS**

## Developer utilities

The `.dev/utils` subdirectory contains three shell scripts that streamline the
development process, as described in the following.

The first one is the `setup.sh` script, which sets up your environment for
building and running GISS-GC. This includes setting environment variables
and activating a Spack environment. It's likely that this script will need
some modification to be usable on your system. (Open the script and check the
parts marked by `NOTE`.) To use the `setup.sh` script, navigate to the directory
where you cloned GISS-GC (which is set as `${GISS_HOME}` in the script) and run
```
cp .dev/utils/setup.sh .
source setup.sh
```

The second utility is the `build.sh` script, which builds the model with given
configuration options. From the same location, run
```
cp .dev/utils/build.sh .
./build.sh --help
```
This should print some help text to the screen that describes what options may
be passed to the script and what they mean. The options can be combined. For
example, running
```
./build.sh --giss-only -f
```
will build GISS Model E *without* GEOS-Chem support and it will remove any
existing builds before doing so. Note that if you run the build script once with
the `--giss-only` option and once without then it will create two separate
builds.

The third utility is the `run.sh` script, which runs a model configuration that
has been built.

Again, you may run the following to get help text on how to use this script:
```
cp .dev/utils/run.sh .
./run.sh --help
```
Note that in order to run the model with `--giss-only` or not, you will need to
have built the model in the same way.

## Configuration scripts

The `.dev/config` subdirectory contains five configuration files for GEOS-Chem
and the HEMCO emissions module that provide good starting points for GISS-GC
development. (Again, these are not intended for general users as they may not
produce scientifically interesting output.) The configuration scripts are
automatically picked up by the developer utility scripts mentioned in the
previous section.

The configuration scripts are:
* `geoschem_config.yml` - top level configuration file for GEOS-Chem.
* `HEMCO_Config.rc` - top level configuration file for HEMCO.
* `HISTORY.rc` - output configuration for GEOS-Chem.
* `HEMCO_Diagn.rc` - output configuration for HEMCO.
* `species_database.yml` - configuration file for species used by GEOS-Chem.
