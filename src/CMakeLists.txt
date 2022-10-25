# -------------- #
# Build packages #
# -------------- #

# The module contains the build commands for each of the Sparse2D subpackages.
# Any new or modified files should be automatically detected and handled. If
# a new subpackage (i.e. a new directory in src) is added then a new section
# should be added to this file with the appropriate build command. All of the
# subpackage builds are optional except for sparse, which contains the core
# Sparse2D code. New or modified Python binding will also be handled
# automatically as long as they are located in the src/python_binding
# directory. Note that the build functions (i.e. build_package and
# build_pybind) are defined in cmake/BuildPkg.cmake and
# cmake/BuildPybind.cmake.

# Define basic libraries required for building all packages
set(basic_libs ${CFITSIO_LIBRARIES} ${FFTW3_LIBRARIES})

# ------ #
# SPARSE #
# ------ #

# List package libraries (in static dependency order)
set(sparse_libs mga2d sparse3d sparse2d sparse1d tools)
# List package dependencies
# - CFITSIO
# - FFTW
set(sparse_deps ${basic_libs})
# Build package and return package libs
build_package(sparse "${sparse_libs}" "${sparse_deps}")

# Define the core libraries required for building packages that use sparse
set(core_libs ${sparse_libs} ${basic_libs})

# ----- #
# MSVST #
# ----- #

if(BUILD_MSVST)

  # List package libraries (in static dependency order)
  set(msvst_libs msvst)
  # List package dependencies
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(msvst_deps ${core_libs})
  # Build package and return package libs
  build_package(msvst "${msvst_libs}" "${msvst_deps}" )

endif()

# ---- #
# MISC #
# ---- #

if(BUILD_MISC)

  # List package libraries (in static dependency order)
  set(misc_libs misc3d misc2d misc1d)
  # List package dependencies
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(misc_deps ${core_libs})
  # Build package
  build_package(misc "${misc_libs}" "${misc_deps}")

endif()

# -- #
# MR #
# -- #

if(BUILD_MR)

  # List package libraries (in static dependency order)
  set(mr_libs mr2d mr1d)
  # List package dependencies
  # - MISC
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(mr_deps ${misc_libs} ${core_libs})
  # Build package
  build_package(mr "${mr_libs}" "${mr_deps}")

endif()

# -- #
# MC #
# -- #

if(BUILD_MC)

  # List package libraries (in static dependency order)
  set(mc_libs mc1d mc2d)
  # List package dependencies
  # - MR
  # - MISC
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(mc_deps ${mr_libs} ${misc_libs} ${core_libs})
  # Build package
  build_package(mc "${mc_libs}" "${mc_deps}")

endif()

# --- #
# MGA #
#---- #

if(BUILD_MGA)

  # List package libraries (in static dependency order)
  set(mga_libs mga1d mga3d)
  # List package dependencies
  # - MC
  # - MR
  # - MISC
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(mga_deps ${mc_libs} ${mr_libs} ${misc_libs} ${core_libs})
  # Build package and return package libs
  build_package(mga "${mga_libs}" "${mga_deps}")

endif()

# ---- #
# MWIR #
# ---- #

if(BUILD_MWIR)

  # List package libraries (in static dependency order)
  set(mwir_libs mwir)
  # List package dependencies
  # - MISC
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(mwir_deps ${misc_libs} ${core_libs})
  # Build package
  build_package(mwir "${mwir_libs}" "${mwir_deps}")

endif()

# ------------------- #
# DICTIONARY LEARNING #
# ------------------- #

# NOTE: Need to to include GSL and Armadillo libraries

if(BUILD_DICLEARN)

  # List package libraries (in static dependency order)
  set(diclearn_libs diclearn)
  # List package dependencies
  # - MR
  # - MISC
  # - GSL
  # - Armadillo
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(diclearn_deps
    ${mr_libs}
    ${misc_libs}
    ${GSL_LIBRARIES}
    ${ARMADILLO_LIBRARIES}
    ${core_libs}
  )
  # Build package
  build_package(diclearn "${diclearn_libs}" "${diclearn_deps}")

endif()

# --- #
# MRS #
# --- #

if(BUILD_MRS)

  # List package libraries (in static dependency order)
  set(mrs_libs mrs mrsp)
  # List package dependencies
  # - HEALPix
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(mrs_deps ${HEALPIX_LIBRARIES} ${core_libs})
  # Build package
  build_package(mrs "${mrs_libs}" "${mrs_deps}" )

endif()

# -------- #
# Astro WL #
# -------- #

if(BUILD_ASTRO_WL)

  # List package libraries (in static dependency order)
  set(astro_wl_libs divcurlwt wl wls)
  # List package dependencies
  # - MRS
  # - HEALPix
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(astro_wl_deps ${mrs_libs} ${HEALPIX_LIBRARIES} ${core_libs})
  # Build package
  build_package(astro_wl "${astro_wl_libs}" "${astro_wl_deps}")

endif()

# --------- #
# Astro Gal #
# --------- #

if(BUILD_ASTRO_GAL)

  # List package libraries (in static dependency order)
  set(astro_gal_libs gal)
  # List package dependencies
  # - MISC
  # - SPARSE
  # - CFITSIO
  # - FFTW
  set(astro_gal_deps ${misc_libs} ${core_libs})
  # Build package
  build_package(astro_gal "${astro_gal_libs}" "${astro_gal_deps}")

endif()

# --------------- #
# Python Bindings #
# --------------- #

if(BUILD_PYBIND)

  # List binding dependencies
  # - MRS
  # - HEALPix
  # - SPARSE
  # - CFITSIO
  # - FFTW
  # Build bindings
  build_pybind("${core_libs}" FALSE)
  # Note: Uncomment the following line when pymrs issues have been resolved
  # if(ONLY_SPARSE)
  #   build_pybind("${core_libs}" FALSE)
  # else()
  #   set(pybind_deps ${mrs_libs} ${HEALPIX_LIBRARIES} ${core_libs})
  #   build_pybind("${pybind_deps}" TRUE)
  # endif()

endif()

# ----------#
# Test libs #
# ----------#

set(test_libs ${core_libs} PARENT_SCOPE)
