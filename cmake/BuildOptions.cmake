# ------------- #
# Build options #
# ------------- #

# This module sets the default build options for Sparse2D. The following
# options can be passed to the CMake command, e.g. cmake -DBUILD_MSVST=OFF.
# Note that some subpackages depend on others so disabling the build of a given
# package may raise an error. The ONLY_SPARSE option will disable the build of
# all of the subpackages except for the core of Sparse2D. The BUILD_DEPS option
# will ignore pre-installed dependencies on the system and build them from
# source. The values of the build options can be seen by running CMake with
# --loglevel=VERBOSE or by passing the -L option.

# Set default build options
option(BUILD_MSVST "BUILD_MSVST" ON)          # Build the MSVST package
option(BUILD_MISC "BUILD_MISC" ON)            # Build the MISC package
option(BUILD_MR "BUILD_MR" ON)                # Build the MR package
option(BUILD_MC "BUILD_MC" ON)                # Build the MC package
option(BUILD_MGA "BUILD_MGA" ON)              # Build the MGA package
option(BUILD_MWIR "BUILD_MWIR" ON)            # Build the MWIR package
option(BUILD_DICLEARN "BUILD_DICLEARN" ON)    # Build the DICLEARN package
option(BUILD_MRS "BUILD_MRS" ON)              # Build the MRS package
option(BUILD_ASTRO_WL "BUILD_ASTRO_WL" ON)    # Build the ASTRO_WL package
option(BUILD_ASTRO_GAL "BUILD_ASTRO_GAL" ON)  # Build the ASTRO_GAL package
option(BUILD_PYBIND "BUILD_PYBIND" ON)        # Build Python bindings
option(ONLY_SPARSE "ONLY_SPARSE" OFF)         # Only build the SPARSE package
option(USE_FFTW "USE_FFTW" ON)                # Use FFTW libraries
option(BUILD_CFITSIO "BUILD_CFITSIO" OFF)     # Build CFITSIO from source
option(BUILD_FFTW3 "BUILD_FFTW3" OFF)         # BUILD FFTW3 from source
option(BUILD_HEALPIX "BUILD_HEALPIX" OFF)     # BUILD HEALPix from source
option(BUILD_GSL "BUILD_GSL" OFF)             # BUILD GSL from source (not implemented)
option(BUILD_DEPS "BUILD_DEPS" OFF)           # BUILD all dependencies from source

# If ONLY_SPARSE is set disable all other package builds
if(ONLY_SPARSE)
  set(BUILD_MGA OFF)
  set(BUILD_MSVST OFF)
  set(BUILD_MISC OFF)
  set(BUILD_MR OFF)
  set(BUILD_MC OFF)
  set(BUILD_MWIR OFF)
  set(BUILD_DICLEARN OFF)
  set(BUILD_MRS OFF)
  set(BUILD_ASTRO_WL OFF)
  set(BUILD_ASTRO_GAL OFF)
endif()

# If BUILD_DEPS is set build all dependencies from source
if(BUILD_DEPS)
  set(BUILD_CFITSIO ON)
  set(BUILD_FFTW3 ON)
  set(BUILD_HEALPIX ON)
endif()

# ------------------- #
# Check build options #
# ------------------- #

set(report_msg "Please check your build options and try running CMake again.")

if(BUILD_MR AND (NOT BUILD_MISC OR NOT USE_FFTW))
  message(FATAL_ERROR "The MR package depends on MISC and FFTW. ${report_msg}")
endif()

if(BUILD_MC AND NOT BUILD_MR)
  message(FATAL_ERROR "The MC package depends on MISC and MR. ${report_msg}")
endif()

if(BUILD_MGA AND NOT BUILD_MC)
  message(FATAL_ERROR "The MGA package depends on MISC, MR and MC. ${report_msg}")
endif()

if(BUILD_MWIR AND NOT BUILD_MISC)
  message(FATAL_ERROR "The MWIR package depends on MISC. ${report_msg}")
endif()

if(BUILD_DICLEARN AND NOT BUILD_MR)
  message(FATAL_ERROR "The DICLEARN package depends on MISC and MR. ${report_msg}")
endif()

if(BUILD_ASTRO_WL AND (NOT BUILD_MRS OR NOT USE_FFTW))
  message(FATAL_ERROR "The ASTRO_WL package depends on MRS and FFTW. ${report_msg}")
endif()

if(BUILD_ASTRO_GAL AND (NOT BUILD_MISC OR NOT USE_FFTW))
  message(FATAL_ERROR "The ASTRO_GAL package depends on MISC and FFTW. ${report_msg}")
endif()

if(BUILD_GSL)
  message(FATAL_ERROR "Sorry, build commands for GSL have not been implemented.")
endif()

# -------------------- #
# Report build options #
# -------------------- #

message(VERBOSE "Sparse2D build summary:")
message(VERBOSE "  BUILD_MSVST: ${BUILD_MSVST}")
message(VERBOSE "  BUILD_MISC: ${BUILD_MISC}")
message(VERBOSE "  BUILD_MR: ${BUILD_MR}")
message(VERBOSE "  BUILD_MC: ${BUILD_MC}")
message(VERBOSE "  BUILD_MGA: ${BUILD_MGA}")
message(VERBOSE "  BUILD_MWIR: ${BUILD_MWIR}")
message(VERBOSE "  BUILD_DICLEARN: ${BUILD_DICLEARN}")
message(VERBOSE "  BUILD_MRS: ${BUILD_MRS}")
message(VERBOSE "  BUILD_ASTRO_WL: ${BUILD_ASTRO_WL}")
message(VERBOSE "  BUILD_ASTRO_GAL: ${BUILD_ASTRO_GAL}")
message(VERBOSE "  BUILD_PYBIND: ${BUILD_PYBIND}")
message(VERBOSE "  ONLY_SPARSE: ${ONLY_SPARSE}")
message(VERBOSE "  USE_FFTW: ${USE_FFTW}")
message(VERBOSE "  BUILD_CFITSIO: ${BUILD_CFITSIO}")
message(VERBOSE "  BUILD_FFTW3: ${BUILD_FFTW3}")
message(VERBOSE "  BUILD_HEALPIX: ${BUILD_HEALPIX}")
message(VERBOSE "  BUILD_GSL: ${BUILD_GSL}")
message(VERBOSE "  BUILD_DEPS: ${BUILD_DEPS}")
