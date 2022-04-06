# Sparse2D
[![build](https://github.com/CosmoStat/Sparse2D/workflows/CI/badge.svg)](https://github.com/CosmoStat/Sparse2D/actions?query=workflow%3ACI)
[![cpp](https://img.shields.io/badge/language-C%2B%2B-red)](https://isocpp.org/std/the-standard)[![bigmac](https://img.shields.io/github/v/release/sfarrens/bigmac?label=BigMac&logo=data%3Aimage%2Fpng%3Bbase64%2CiVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAYAAAAfSC3RAAAABGdBTUEAALGPC%2FxhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAeGVYSWZNTQAqAAAACAAFARIAAwAAAAEAAQAAARoABQAAAAEAAABKARsABQAAAAEAAABSASgAAwAAAAEAAgAAh2kABAAAAAEAAABaAAAAAAAAAEgAAAABAAAASAAAAAEAAqACAAQAAAABAAAADqADAAQAAAABAAAADgAAAAABAtxLAAAACXBIWXMAAAsTAAALEwEAmpwYAAACZmlUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyIKICAgICAgICAgICAgeG1sbnM6ZXhpZj0iaHR0cDovL25zLmFkb2JlLmNvbS9leGlmLzEuMC8iPgogICAgICAgICA8dGlmZjpPcmllbnRhdGlvbj4xPC90aWZmOk9yaWVudGF0aW9uPgogICAgICAgICA8dGlmZjpSZXNvbHV0aW9uVW5pdD4yPC90aWZmOlJlc29sdXRpb25Vbml0PgogICAgICAgICA8ZXhpZjpDb2xvclNwYWNlPjE8L2V4aWY6Q29sb3JTcGFjZT4KICAgICAgICAgPGV4aWY6UGl4ZWxYRGltZW5zaW9uPjUwPC9leGlmOlBpeGVsWERpbWVuc2lvbj4KICAgICAgICAgPGV4aWY6UGl4ZWxZRGltZW5zaW9uPjUwPC9leGlmOlBpeGVsWURpbWVuc2lvbj4KICAgICAgPC9yZGY6RGVzY3JpcHRpb24%2BCiAgIDwvcmRmOlJERj4KPC94OnhtcG1ldGE%2BCrh6Oa8AAAI%2BSURBVCgVhVLPaxNBFH7zY7tpkmKp2UgrwSrShkSNxYpFRXJqNVVQJBWpF1E8Fjx5HfwTPFXoQemtlYKKetCSFAsiQpHSVBQ1LWKlJK3Nj26y2d0ZZ7bkUhQHHvvt%2B7733rezD8GuMzWVJmkjjqCQEwBpAGPZw2hk2t0l3XkVApBgDP%2BVlEkmOaVp8h5QCYRATgBYenrrDOdiiGhaBDgHV6AfRPDX8csTbxXf1KImENLikt7xWNPwKKVUClQ3aUC44EiTtu284D7zWmJocltROJtlRHVaJIF77UEyum3adrlSU%2BGUq3WnVKlbrmtLBa%2F43aKjtOpQgKwHWnTS3XBsIMgFDkiTo6RWcL%2Bvhdp2feHYlUfXlVB9K0KMU5jz6mD%2BTW7tUCwCWtCHdIJBI0g24BxhFzcweq9ULx9c0FNjzFKYFnKGdykHOtu6a7MfYD1sgRMMQbXhhxbNxbptQ22xdFqJU2OvLAaAZXBsxAvezZZa%2FT87zw9Ab%2FQonOvX4XbqCwyeNPFhw8f7Lx5JzI%2BPTMpaQxXJSYisJlfwShb49KnPN2qrjT5zK2CVKyG0WQnzfGEfLzptIles1rO%2Fv58whxtdhVlrZo4BpUZs5%2F%2Fd7I1O9NT2DuOqHb7z6xkMlqPgJw48QV%2FhftdZGgl2rCda3fFlKIOq8WwyJn0zaeHbwz35%2FMKlmY1PxwcC%2B0PSGrwz1zautvd8PGjEnqO%2Bu1tNreK8w8S%2F1%2B2%2FGrktKJNJ0mQGqFp2FQqrnOKaDdTzD77QAxggF14iAAAAAElFTkSuQmCC)](https://github.com/sfarrens/bigmac)


Sparse2D provides an array of sparsity based tools and a convenient C++ library for performing various wavelet transforms.

This package is part of the iSAP suite, available on the [CosmoStat website](http://www.cosmostat.org/software/isap).

Python bindings to Sparse2D are provided in [PySAP](https://github.com/CEA-COSMIC/pysap). PySAP handles Sparse2D installation internally.

## Installation instructions

### Prerequisites

   - C/C++ compiler
   - [CMake](http://www.cmake.org)
   - [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) ( >V3.31)
   - [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
   - [BigMac](https://github.com/sfarrens/bigmac) (For macOS AppleClang compiler)

Please use a package management tool to properly install cfistio and the other dependencies on linux (*e.g.* `apt` on Ubuntu or `brew` on macOS).

### Docker

If you have [Docker](https://www.docker.com/) installed, you can pull the latest build of the Sparse2D image from [Docker Hub](TBD).

```bash
$ docker pull cosmostat/sparse2d
```

No further installation is required.

To run this image on data in your current working directory, simply run:

```bash
$ docker run -v ${PWD}:/workdir -it cosmostat/sparse2d /bin/bash -c "cd workdir && <EXECUTABLE>"
```

where `<EXECUTABLE>` is one of the Sparse2D binaries. The reference to `${PWD}` can be replaced by the path to any directory on your system and options can be passed to `<EXECUTABLE>` inside the double quotes.

### Homebrew

Sparse2D can be built on macOS using [Homebrew](https://brew.sh/).

```bash
$ brew tap cosmostat/science
$ brew install sparse2d
```

The Homebrew formula handles all of the required dependencies.

### Building from Source

#### Basic Sparse2D

Clone the repository:

```
$ git clone https://github.com/CosmoStat/Sparse2D.git
$ cd Sparse2D
```

Create a build directory for Sparse2D:

```bash
$ mkdir build
$ cd build
```

Build Sparse2D:

```bash
$ cmake ..
$ make
$ make install
```

#### MSVST

Additional MSVST binaries can be build with the following option:

```bash
$ cmake .. -DBUILD_MSVST=ON
```

#### FFTW

Sparse2D can also be built using pre-installed [FFTW](http://www.fftw.org/) libraries with the following option:

```bash
$ cmake .. -DUSE_FFTW=ON
```

Alternatively, Sparse2D can build FFTW libraries from source as follows:

```bash
$ cmake .. -DUSE_FFTW=ON -DBUILD_FFTW=ON
```

#### nFFT

Sparse2D can also build [nFFT](https://github.com/NFFT/nfft) libraries with the following options:

```bash
$ cmake .. -DBUILD_FFTW=ON -DBUILD_NFFT=ON
```

#### Non-default Compiler

Finally, if you wish to build using a compiler other than the default on your
system (*e.g.* `gcc` on macOS) you can do so as follows:

```bash
$ CC=gcc CXX=g++ cmake ..
```

## Usage

The two main executables of the package are `mr_transform` and `mr_filter`, see the instructions in the [README](./README) folder for usage examples.
