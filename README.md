# Sparse2D
[![build](https://github.com/CosmoStat/Sparse2D/workflows/CI/badge.svg)](https://github.com/CosmoStat/Sparse2D/actions?query=workflow%3ACI)
[![cpp](https://img.shields.io/badge/language-C%2B%2B-red)](https://isocpp.org/std/the-standard)[![bigmac](https://img.shields.io/github/v/release/sfarrens/bigmac?label=BigMac&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAYAAAAfSC3RAAAABGdBTUEAALGPC%2FxhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAeGVYSWZNTQAqAAAACAAFARIAAwAAAAEAAQAAARoABQAAAAEAAABKARsABQAAAAEAAABSASgAAwAAAAEAAgAAh2kABAAAAAEAAABaAAAAAAAAAEgAAAABAAAASAAAAAEAAqACAAQAAAABAAAADqADAAQAAAABAAAADgAAAAABAtxLAAAACXBIWXMAAAsTAAALEwEAmpwYAAACZmlUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyIKICAgICAgICAgICAgeG1sbnM6ZXhpZj0iaHR0cDovL25zLmFkb2JlLmNvbS9leGlmLzEuMC8iPgogICAgICAgICA8dGlmZjpPcmllbnRhdGlvbj4xPC90aWZmOk9yaWVudGF0aW9uPgogICAgICAgICA8dGlmZjpSZXNvbHV0aW9uVW5pdD4yPC90aWZmOlJlc29sdXRpb25Vbml0PgogICAgICAgICA8ZXhpZjpDb2xvclNwYWNlPjE8L2V4aWY6Q29sb3JTcGFjZT4KICAgICAgICAgPGV4aWY6UGl4ZWxYRGltZW5zaW9uPjUwPC9leGlmOlBpeGVsWERpbWVuc2lvbj4KICAgICAgICAgPGV4aWY6UGl4ZWxZRGltZW5zaW9uPjUwPC9leGlmOlBpeGVsWURpbWVuc2lvbj4KICAgICAgPC9yZGY6RGVzY3JpcHRpb24%2BCiAgIDwvcmRmOlJERj4KPC94OnhtcG1ldGE%2BCrh6Oa8AAAI%2BSURBVCgVhVLPaxNBFH7zY7tpkmKp2UgrwSrShkSNxYpFRXJqNVVQJBWpF1E8Fjx5HfwTPFXoQemtlYKKetCSFAsiQpHSVBQ1LWKlJK3Nj26y2d0ZZ7bkUhQHHvvt%2B7733rezD8GuMzWVJmkjjqCQEwBpAGPZw2hk2t0l3XkVApBgDP%2BVlEkmOaVp8h5QCYRATgBYenrrDOdiiGhaBDgHV6AfRPDX8csTbxXf1KImENLikt7xWNPwKKVUClQ3aUC44EiTtu284D7zWmJocltROJtlRHVaJIF77UEyum3adrlSU%2BGUq3WnVKlbrmtLBa%2F43aKjtOpQgKwHWnTS3XBsIMgFDkiTo6RWcL%2Bvhdp2feHYlUfXlVB9K0KMU5jz6mD%2BTW7tUCwCWtCHdIJBI0g24BxhFzcweq9ULx9c0FNjzFKYFnKGdykHOtu6a7MfYD1sgRMMQbXhhxbNxbptQ22xdFqJU2OvLAaAZXBsxAvezZZa%2FT87zw9Ab%2FQonOvX4XbqCwyeNPFhw8f7Lx5JzI%2BPTMpaQxXJSYisJlfwShb49KnPN2qrjT5zK2CVKyG0WQnzfGEfLzptIles1rO%2Fv58whxtdhVlrZo4BpUZs5%2F%2Fd7I1O9NT2DuOqHb7z6xkMlqPgJw48QV%2FhftdZGgl2rCda3fFlKIOq8WwyJn0zaeHbwz35%2FMKlmY1PxwcC%2B0PSGrwz1zautvd8PGjEnqO%2Bu1tNreK8w8S%2F1%2B2%2FGrktKJNJ0mQGqFp2FQqrnOKaDdTzD77QAxggF14iAAAAAElFTkSuQmCC)](https://github.com/sfarrens/bigmac)

Sparse2D provides an array of sparsity-based signal processing tools and a convenient C++ library for performing various wavelet transforms.

You can find out more about the applications of Sparse2D on the [CosmoStat website](http://www.cosmostat.org).

The core Sparse2D libraries are used as a backend for [PySAP](https://github.com/CEA-COSMIC/pysap).

## Docker installation

### Pull the Docker image

If you have [Docker](https://www.docker.com/) installed, you can pull the latest build of the [Sparse2D image](https://github.com/cosmostat/Sparse2D/pkgs/container/sparse2d) as follows:

```bash
docker pull ghcr.io/cosmostat/sparse2d:master
```

No further installation is required.

### Run a Docker container

To run a container on data in your current working directory, simply run:

```bash
docker run -v ${PWD}:/workdir --rm ghcr.io/cosmostat/sparse2d:master <EXECUTABLE> <ARGUMENTS>
```

where `<EXECUTABLE>` is one of the Sparse2D binaries and `<ARGUMENTS>` are the corresponding command line arguments for this executable. The reference to `${PWD}` can be replaced by the path to any directory on your system.

For example, to run a bspline wavelet transform on a FITS image called `myfile.fits` you would run:

```bash
docker run -v ${PWD}:/workdir --rm ghcr.io/cosmostat/sparse2d:master mr_transform -t 2 myfile.fits myoutput.mr
```

> Tip: If you don't want to constantly write the full Docker run command you can create an alias *e.g.*:
> ```bash
> alias sparse2d="docker run -v ${PWD}:/workdir --rm ghcr.io/cosmostat/sparse2d:master"
> ```
> then you can simply run *e.g.*:
> ```bash
> sparse2d mr_transform -h
> ```

### Launch a Jupyter notebook

You can also run a Jupyter notebook with a Docker container as the backend in order to use the `pysparse` Python bindings to some Sparse2D tools.

```bash
docker run -p 8888:8888 -v ${PWD}:/workdir --rm ghcr.io/cosmostat/sparse2d:master notebook
```

## Homebrew installation

### Standard installation

Sparse2D can be installed on macOS using [Homebrew](https://brew.sh/).

```bash
brew tap cosmostat/science
brew install sparse2d
```

The [Homebrew formula](https://github.com/CosmoStat/homebrew-science/blob/master/sparse2d.rb) handles all of the required dependencies.

### Build options

By default, the Homebrew formula builds the full Sparse2D package including the `pysparse` Python bindings. Some options are available to limit of the scope of the build. These options can be listed as follows (after tapping `cosmostat/science`):

```bash
brew info sparse2d
```

### Python bindings

The `pysparse` bindings will be saved to `/usr/local/opt/sparse2d/python` by default. You will need to add this to your `PYTHONPATH` in order to access the bindings. For example, for Bash you would run the following:

```bash
PYTHONPATH="/usr/local/opt/sparse2d/python:$PYTHONPATH"
```

Note that `pysparse` will be built using the Python executable installed by Homebrew. These bindings will only work with the same version of Python.

## Building from Source

### Dependencies

In order to build Sparse2D from source, you will need to ensure you have installed all of the following dependencies. If possible, please use a package management tool to properly install them (*e.g.* `apt` on Ubuntu or `brew` on macOS).

#### Required

   - C/C++ compiler (*e.g.* `gcc` or `clang`)
   - [Catch2](https://github.com/catchorg/Catch2) (>= v3)
   - [CMake](http://www.cmake.org) (>= v3.12)
   - [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/)
   - [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)

#### Optional

   - [Armadillo](https://arma.sourceforge.net/)
     (not required if `ONLY_SPARSE=ON` or `ONLY_INPAINT=ON`)
   - [BigMac](https://github.com/sfarrens/bigmac)
     (only required if using macOS `clang`)(>= v0.0.6)
   - [FFTW](https://www.fftw.org/)
     (not required if `ONLY_SPARSE=ON` and `USE_FFTW=OFF`)
   - [GSL](https://www.gnu.org/software/gsl/) (not required if `ONLY_SPARSE=ON`)
   - [HEALPix](https://healpix.sourceforge.io/)
     (not required if `ONLY_SPARSE=ON` or `ONLY_INPAINT=ON`)
   - [libomp](https://openmp.llvm.org/) (only required if using macOS `clang`)
   - [Pybind11](https://pybind11.readthedocs.io/)
     (not required if `BUILD_PYBIND=OFF` or `ONLY_INPAINT=ON`)
   - [Python](https://www.python.org/)
     (not required if `BUILD_PYBIND=OFF` or `ONLY_INPAINT=ON`)

### Full Sparse2D build

Download the [latest release of Sparse2D](https://github.com/CosmoStat/Sparse2D/releases).

Create a build directory inside the root directory of the Sparse2D package:

```bash
cd Sparse2D
mkdir build
cd build
```

Build Sparse2D:

```bash
cmake ..
make
make install
```

> Tip: You can significantly increase the speed of compilation by using the
> `--jobs` (or `-j`) option for `make`, which builds the targets in parallel.
> For example, to use 8 cores you would run `make -j 8`.

### Build options

Sparse2D supports the following CMake build options:

- `BUILD_MSVST` (default `ON`): Build the MSVST package
- `BUILD_MISC` (default `ON`): Build the MISC package
- `BUILD_MR` (default `ON`): Build the MR package
- `BUILD_MC` (default `ON`): Build the MC package
- `BUILD_MGA` (default `ON`): Build the MGA package
- `BUILD_MWIR` (default `ON`): Build the MWIR package
- `BUILD_DICLEARN` (default `ON`): Build the DICLEARN package
- `BUILD_MRS` (default `ON`): Build the MRS package
- `BUILD_ASTRO_WL` (default `ON`): Build the ASTRO_WL package
- `BUILD_ASTRO_GAL` (default `ON`): Build the ASTRO_GAL package
- `BUILD_PYBIND` (default `ON`): Build Python bindings
- `ONLY_SPARSE` (default `OFF`): Only build the SPARSE package
- `ONLY_INPAINT` (default `OFF`): Only build the packages required for inpainting
- `USE_FFTW` (default `ON`): Use FFTW libraries (optional for sparse package only)
- `BUILD_CFITSIO` (default `OFF`): Build CFITSIO from source
- `BUILD_FFTW3` (default `OFF`): BUILD FFTW3 from source
- `BUILD_HEALPIX` (default `OFF`): BUILD HEALPix from source
- `BUILD_GSL` (default `OFF`): BUILD GSL from source (not implemented)
- `BUILD_DEPS` (default `OFF`): BUILD all dependencies from source

To use these options prepend `-D` when running the `cmake` command, *e.g.*:

```bash
cmake .. -DONLY_SPARSE=ON
```

### Non-default compiler

Finally, if you wish to build using a compiler other than the default on your
system (*e.g.* `gcc` on macOS) you can do so as follows:

```bash
CC=gcc CXX=g++ cmake ..
```

## Python bindings

By default Sparse2D will build `pysparse`, which includes Python bindings to some Sparse2D tools. This file will need to be in your `PYTHONPATH` in order for these bindings to available outside of the build directory.

`pysparse` can be imported in a given Python session as follows:

```python
import pysparse
```

To see what tools are available in this module run the `help` function.

```python
help(pysparse)
```

For example, to run a bspline wavelet transform on an image (stored in memory as a Numpy array) called `myimage` you would run:

```python
import pysparse
mrt = pysparse.MRTransform(2)
result = mrt.transform(myimage)
```
