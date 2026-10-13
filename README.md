# SWAT+MODFLOW

This repository starts from the original SWAT+MODFLOW source archive and records GNU Fortran compatibility changes as separate commits.

## Build

Requires GNU Fortran, Meson, Ninja, and GNU Make. From the repository root:

```sh
make build
```

The executable is `build/swatplus-modflow`. A later `make build` recompiles only sources affected by changes. `make clean` removes the build directory. `make` lists available targets.

Meson is the compiler build backend because Fortran `USE` statements need module dependency ordering across the SWAT+, MODFLOW, Linkage, and RT3D directories. The Makefile provides the lab's simple workflow interface. You can also use Meson directly:

```sh
meson setup build
meson compile -C build
```

This build currently targets GNU Fortran. The flags for fixed and free source line lengths and legacy procedure arguments are defined in `meson.build`.

## Build environment

All 727 Fortran source files compiled and linked successfully on October 12,
2026, at commit `f6e1037`. Integration commit `5cc0917` has the same file
contents. The toolchain and environment versions were recorded on October
13, 2026:

| Component | Version |
| --- | --- |
| GNU Fortran (GCC) | 15.2.0 |
| Meson | 1.8.3 |
| Ninja | 1.13.1 |
| GNU Make | 4.4.1 |
| Python (used by Meson) | 3.12.11 |
| GNU binutils linker (`ld.bfd`) | 2.45-slack151 |
| Linux kernel | 6.12.41 |
| Architecture | x86_64 |

The successful build used the repository defaults:

```text
buildtype=debug
fortran_std=legacy
-ffixed-line-length-none
-ffree-line-length-none
-fallow-argument-mismatch
```

The source changes address compiler-specific extensions and type/interface
errors reported by GNU Fortran 15.2.0. The table and options above provide a
reference environment for reproducing the build.

To record the tools on another machine:

```sh
gfortran --version
meson --version
ninja --version
make --version
python3 --version
ld --version
uname -r
uname -m
```

Meson records the selected compiler in
`build/meson-info/intro-compilers.json` and toolchain configuration details in
`build/meson-logs/meson-log.txt`.
