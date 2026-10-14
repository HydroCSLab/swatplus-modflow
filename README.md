# SWAT+MODFLOW

This repository provides scripts for building SWAT+MODFLOW with GNU Fortran using the [original source code](https://zenodo.org/records/14674981/files/SourceCode_2024_11_07.zip). The original source files are intentionally preserved unchanged. The build script records and applies the compatibility patches during the build, making the changes explicit and reproducible.

The source code comes from [SWAT+MODFLOW model: code and example models](https://zenodo.org/records/14674981), published by Ryan Bailey on Zenodo (version v1, January 16, 2025). The repository contains an unmodified copy of [SourceCode_2024_11_07.zip](https://zenodo.org/records/14674981/files/SourceCode_2024_11_07.zip).

## Build the included source code

From the repository directory, run:

```bash
./build.sh
```

## Download and build the original source code

To download the original source archive from Zenodo and build it, run:
```bash
./get-build.sh
```
