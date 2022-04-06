#################
# Sparse2D 2017 #
#################

cmake_minimum_required(VERSION 3.12.0)
project(sparse2d)

# Include modules in CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules/")

# Load CMake tools
include(ExternalProject)
include(FindPkgConfig)

# Load custom CMake functions
include(functions)

# Set build type
if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "CMake Build Type: ${CMAKE_BUILD_TYPE}")

# ----------------- #
# Find Dependencies #
# ----------------- #

# Locate OpenMP
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  find_package(BigMac 0.0.1 REQUIRED)
else()
  find_package(OpenMP REQUIRED)
endif()

# Locate CFITSIO
find_pkg(CFITSIO cfitsio)

# Set build location for external modules
set(MODULE_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/module_build/")

# Set compilation flags
if("${CMAKE_CXX_FLAGS}" STREQUAL "")
  set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -w -std=c++17 -DNO_DISP_IO -fPIC -Wno-write-strings\
    ${OpenMP_CXX_FLAGS} ${BigMac_NOWARN}"
  )
endif()
set(CMAKE_CXX_FLAGS_RELEASE
  "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -g0 -fomit-frame-pointer"
)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -O0")

# ---- #
# FFTW #
# -----#

# Optional use of FFTW
option(USE_FFTW "Use FFTW for Fourier transformation" OFF)

# Optional build of FFTW
option(BUILD_FFTW "Build FFTW libraries" OFF)
include(BuildFFTW)

# ---- #
# NFFT #
# -----#

# Optional build of NFFT
option(BUILD_NFFT "Download and build NFFT" OFF)
include(BuildNFFT)

# ------------ #
# Sparse2D STD #
# -------------#

# Build standard Sparse2D libs
set(sparse2d_libs mga2d sparse3d sparse2d sparse1d tools)
foreach(library ${sparse2d_libs})
  build_lib(${library})
endforeach()

# Find all Sparse2D targets
find_targets(sparse2d_targets src cc)

# Build binaries
foreach(program ${sparse2d_targets})
  build_bin(${program} "${sparse2d_libs}" src cc)
endforeach(program)

# Install libraries
install(TARGETS ${sparse2d_libs} DESTINATION lib EXPORT sparse2d_libs)
install(EXPORT sparse2d_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${sparse2d_targets} DESTINATION bin)

# ----- #
# MSVST #
# ------#

option(BUILD_MSVST "Build MSVST library" OFF)
message(STATUS "MSVST Build: ${BUILD_MSVST}")

if(BUILD_MSVST)

  # Build MSVST lib
  build_lib(msvst)
  set(msvst_libs msvst tools)

  # Find all MSVST targets
  find_targets(msvst_targets src/msvst cc)

  # Build binaries
  foreach(program ${msvst_targets})
    build_bin(${program} "${msvst_libs}" src/msvst cc)
  endforeach(program)

  # Install binaraies
  install(TARGETS ${msvst_targets} DESTINATION bin)

endif()

# ---------- #
# UNIT TESTS #
# -----------#

option(RUN_TESTS "Build and run unit tests" ON)
message(STATUS "Run Tests: ${RUN_TESTS}")

if(RUN_TESTS)

  # Add unit tests
  enable_testing()

  # Find all unit test targets
  find_targets(unit_tests tests cpp)

  # Build binaries
  foreach(program ${unit_tests})
    build_bin(${program} "${sparse2d_libs}" tests cpp)
    add_test(${program} ${program})
  endforeach(program)

endif()
