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
