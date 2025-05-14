############
# Sparse2D #
############

# -------------- #
# Project Set Up #
# -------------- #

cmake_minimum_required(VERSION 3.12.0)
project(
	Sparse2D
	VERSION 3.0
	DESCRIPTION "Sparsity-based signal processing tools"
	LANGUAGES CXX
)

# Set policy for later versions of CMake
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
		cmake_policy(SET CMP0135 NEW)
endif()

# Prohibit in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
	message(FATAL_ERROR "In-source build is not permitted!")
endif()

# Include custom modules in CMake module path
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Load CMake tools
include(ExternalProject)
include(FetchContent)
include(FindPkgConfig)

# Load custom CMake functions: find_pkg, build_package, build_pybind
include(FindPkg)
include(BuildPackage)
include(BuildPybind)
include(BuildTests)

# Load Sparse2D build options
include(BuildOptions)

# Set default build type as Release
if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "CMake Build Type: ${CMAKE_BUILD_TYPE}")

# Set RPATH options
set(CMAKE_BUILD_RPATH "@loader_path/../lib")
set(CMAKE_INSTALL_RPATH "@loader_path/../lib")

# ----------------- #
# Find Dependencies #
# ----------------- #

if(BUILD_PYBIND)

	# Locate Python
	set(Python_FIND_STRATEGY LOCATION)
	find_package(Python COMPONENTS Interpreter Development)
	include_directories(${Python_INCLUDE_DIRS})
	link_directories(${Python_LIBRARY_DIRS})

	# Source Pybind11
	include(FetchPybind11)

endif()

# Use BigMac for local clang build
if(
	"${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"
	AND NOT "${CMAKE_CXX_COMPILER}" MATCHES "Homebrew"
)
  find_package(BigMac 0.0.6 REQUIRED)
endif()

# Locate OpenMP
find_package(OpenMP REQUIRED)

# Locate Catch2
find_package(Catch2 3)

# Locate CFITSIO
find_pkg(CFITSIO cfitsio)

if(USE_FFTW)

	# Locate FFTW
	find_pkg(FFTW3 "fftw3 fftw3f")

endif()

if(BUILD_MR)

	# Locate GSL
	find_pkg(GSL gsl)

endif()

if(BUILD_MRS)

	# Locate HEALPix
	find_pkg(HEALPIX healpix_cxx)

endif()

if(BUILD_DICLEARN)

	# Locate GSL
	find_pkg(GSL gsl)

	# Locate Armadillo
	find_package(Armadillo)
	include_directories(${ARMADILLO_INCLUDE_DIRS})

endif()

# ------------------ #
# Compilation Set Up #
# ------------------ #

# Set CXX options
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set CXX compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BigMac_NOWARN} -DNO_DISP_IO -w -Wno-everything -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION")
if(USE_FFTW)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FFTW")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")

# -------------- #
# Build Sources  #
# -------------- #

# Build packages in the source directory
add_subdirectory(src)

# ---------- #
# Unit Tests #
# -----------#

if(Catch2_FOUND)

	# Enable unit Tests
	enable_testing()

	# Build unit tests
	add_subdirectory(tests)

endif()
