# Tools Directory Readme 

The tools directory contains: 

* bash/ - a folder with scripts that convert inp to nml input files 
* ftn2src.sh - a script that converts ftn/ folder with switch preprocesssing to src/ folder with CPP if defs
* switch2cpp.awk - the script for converting ftn -> CPP 

# To convert ftn to src 

From the model/tools directory run the following script: 
```bash
sh ftn2src.sh
```

Your ftn directory should now be converted to a src directory.

# Cleaning Build Artifacts with CMake

To clean up and remove all object files and executable binaries (including those associated with unit tests) generated during a CMake build:

### 1. Target Clean
To remove object files, intermediate build artifacts, and compiled executables (including model binaries and unit test executables built by CMake targets) while preserving the CMake configuration:

From the top-level repository directory (assuming the build directory is `build`):
```bash
cmake --build build --target clean
```

Alternatively, from within the `build/` directory (when using a Makefile generator):
```bash
make clean
```

### 2. Complete Build Clean (Recommended)
Because CMake performs out-of-source builds, all generated object files, executables (model binaries and unit test executables), CMake cache files, and intermediate build artifacts reside inside the build directory. To perform a complete cleanup and remove all object files, executables, unit test artifacts, and CMake configurations:

```bash
rm -rf build
```

This completely removes all generated object and executable files and restores the workspace to a pristine state.
