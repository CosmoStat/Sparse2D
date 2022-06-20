#################
# Sparse2D 2017 #
#################
# to run in debug mode, call cmake with command:
#       cmake -DCMAKE_BUILD_TYPE=Debug


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
    set(OpenMP_INCLUDE_PATH "/usr/local/include")
    set(OpenMP_LIB_PATH "/usr/local/lib")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -lomp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_CXX_LIBRARIES "${OpenMP_LIB_PATH}/libomp.dylib")
    include_directories(${OpenMP_INCLUDE_PATH})
    link_directories(${OpenMP_LIB_PATH})
else()
  find_package(OpenMP REQUIRED)
endif()

# Locate Python
find_package(Python COMPONENTS Interpreter Development)
include_directories(${Python_INCLUDE_DIRS})
link_directories(${Python_LIBRARY_DIRS})

# Locate CFITSIO
find_pkg(CFITSIO cfitsio)

# Locate GSL using pkg-config
pkg_check_modules(GSL gsl)
if (GSL_FOUND)
    message("GSL package: == YES == ")
    include_directories(${GSL_INCLUDE_DIRS})
    link_directories(${GSL_LIBRARY_DIRS})
else()
    message("GSL not found ... ")
endif(GSL_FOUND)

# Set build location for external modules
set(MODULE_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/module_build/")

# Set compilation flags
if("${CMAKE_CXX_FLAGS}" STREQUAL "")
#   set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -w -std=c++17 -DNO_DISP_IO -fPIC -Wno-write-strings ${OpenMP_CXX_FLAGS} ${BigMac_NOWARN}")
#   set(CMAKE_CXX_FLAGS "-std=c++11 -DNO_DISP_IO -ggdb3 -fPIC -O2 -ffast-math -fomit-frame-pointer ${APPLE_FLAGS} ${OpenMP_CXX_FLAGS} -Wno-write-strings -DNDEBUG")
    set(CMAKE_CXX_FLAGS "-std=c++17 -DNO_DISP_IO -ggdb3 -fPIC -O2 -ffast-math -fomit-frame-pointer ${BigMac_NOWARN} ${OpenMP_CXX_FLAGS} -Wno-write-strings -DNDEBUG")
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS} -g0 -fomit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} -O0")

# Flags for MAC OSX
if(APPLE)
    set(CMAKE_MACOSX_RPATH 1)
    set(APPLE_FLAGS "-DMACOS")
  else(APPLE)
    set(APPLE_FLAGS "")
endif(APPLE)
  
# ---- #
# FFTW #
# -----#

# Optional use of FFTW
option(USE_FFTW "Use FFTW for Fourier transformation" OFF)

# Locate FFTW using pkg-config
pkg_check_modules(FFTW fftw3)
if (FFTW_FOUND)
    include_directories(${FFTW_INCLUDE_DIRS})
    link_directories(${FFTW_LIBRARY_DIRS})
    message(STATUS "FFTW: ${FFTW_INCLUDE_DIRS}")
endif(FFTW_FOUND)

 # Optional build of FFTW
 # option(BUILD_FFTW "Build FFTW libraries" ON)
 # include(BuildFFTW)
 # message(STATUS "BUILD_FFTW Build: ${BUILD_FFTW}")

# ---- #
# NFFT #
# -----#

# Optional build of NFFT
option(BUILD_NFFT "Download and build NFFT" OFF)
include(BuildNFFT)

# ------------  COMPILATION OPTIONS ------------  #

set(DO_NOT_COMPILE_FORTRAN_BIN 1)  # we have some fortran code to be included
set(DO_NOT_COMPILE_HEALPIX_BIN 0)
set(DO_NOT_COMPILE_COSMOSTAT_LIBS 0)  # if set to 1, most COSMOSTAT C++ code is not compiled, only what is required
                                      # by pysap is compiled
set(DO_NOT_COMPILE_PYTHON_LIBS 0)     # to have python libraries compiled, we need both COSMOSTAT and PYTHON to be set to 0

if(DO_NOT_COMPILE_COSMOSTAT_LIBS)
    option(BUILD_MSVST "Build MSVST library" ON)
    option(BUILD_MISC "Build BUILD_MISC" ON)
    option(BUILD_MR "Build MR" ON)
    option(BUILD_MC "Build MC" ON)
    option(BUILD_MGA "Build MGA" ON)
    option(BUILD_MWIR "Build MWIR library" ON)
    option(BUILD_ASTRO_GAL "Build ASTRO GAL library" ON)
    if (GSL_FOUND)
        option(BUILD_DICLEARN "Build Dictionary Learning library" ON)
    else()
        option(BUILD_DICLEARN "Build Dictionary Learning library" OFF)
    endif()
    if(DO_NOT_COMPILE_HEALPIX_BIN)
        option(BUILD_MRS "Build BUILD_MRS" OFF)
        option(BUILD_ASTRO_WL "Build ASTRO WL library" OFF)
    else()
        # ------- #
        # HEALPIX #
        # --------#
        find_pkg(HEALPIX healpix_cxx>=3.80)

        if (HEALPIX_FOUND)
            include_directories(${HEALPIX_INCLUDE_DIRS})
            link_directories(${HEALPIX_LIBRARY_DIRS})
            message("HEALPIX package: == YES == ")
            message("DIRLIB = ${HEALPIX_LINK_LIBRARIES}")
            message("INCLIB = ${HEALPIX_INCLUDE_DIRS}")
            message("${HEALPIX_LIBRARIES}")
        else(HEALPIX_FOUND)
            message("======================================================")
            message("NO HEALPIX package. MRS binaries will not be compiled.")
            message("======================================================")
        endif(HEALPIX_FOUND)
    endif(DO_NOT_COMPILE_HEALPIX_BIN)
    
    if(DO_NOT_COMPILE_PYTHON_LIBS)
        option(BUILD_BINDING "Build Python Bindings" OFF)
    else()
        option(BUILD_BINDING "Build Python Bindings" ON)
    endif(DO_NOT_COMPILE_PYTHON_LIBS)

endif(DO_NOT_COMPILE_COSMOSTAT_LIBS)

message("======================================================")
message(STATUS "MSVST: ${BUILD_MSVST}")
message(STATUS "MISC Build: ${BUILD_MISC}")
message(STATUS "MRS Build: ${BUILD_MRS}")
message(STATUS "MR Build: ${BUILD_MR}")
message(STATUS "MC Build: ${BUILD_MC}")
message(STATUS "MGA Build: ${BUILD_MGA}")
message(STATUS "MWIR Build: ${BUILD_MWIR}")
message(STATUS "DICO_LEARN Build: ${BUILD_DICLEARN}")
message(STATUS "ASTRO_WL Build: ${BUILD_ASTRO_WL}")
message(STATUS "ASTRO_GAL Build: ${BUILD_ASTRO_GAL}")
message(STATUS "BINDING Build: ${BUILD_BINDING}")

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

# ------------ #
# MISC STD #
# -------------#

if(BUILD_MISC)
# Build standard Sparse2D libs
set(misc_libs misc1d misc2d misc3d)
foreach(library ${misc_libs})
  build_lib1(${library} src/misc cc h)
endforeach()

set(misc_alllibs ${misc_libs} ${sparse2d_libs} ${FFTW_LIBRARIES})

# Find all Sparse2D targets
find_targets(misc_targets1 src/misc/main2d cc)
find_targets(misc_targets2 src/misc/main3d cc)

# Build binaries
foreach(program ${misc_targets1})
  build_bin(${program} "${misc_alllibs}" src/misc/main2d cc)
endforeach(program)
foreach(program ${misc_targets2})
  build_bin(${program} "${misc_alllibs}" src/misc/main3d cc)
endforeach(program)

# Install libraries
install(TARGETS ${misc_libs} DESTINATION lib EXPORT sparse2d_libs)
install(EXPORT sparse2d_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${misc_targets1} DESTINATION bin)
install(TARGETS ${misc_targets2} DESTINATION bin)
install(TARGETS ${misc_targets3} DESTINATION bin)

endif()

# -- #
# MR #
# -- #

if(BUILD_MR)

# Build standard Sparse2D libs
set(mr_libs mr1d mr2d)
foreach(library ${mr_libs})
  build_lib1(${library} src/mr cc h)
endforeach()

set(mr_alllibs ${mr_libs} ${misc_alllibs})

# Find all Sparse2D targets
find_targets(mr_targets1 src/mr/mrmain1d cc)
find_targets(mr_targets2 src/mr/mrmain2d cc)
find_targets(mr_targets3 src/mr/mrmain3d cc)

# Build binaries
foreach(program ${mr_targets1})
  build_bin(${program} "${mr_alllibs}" src/mr/mrmain1d cc)
endforeach(program)
foreach(program ${mr_targets2})
  build_bin(${program} "${mr_alllibs}" src/mr/mrmain2d cc)
endforeach(program)
foreach(program ${mr_targets3})
  build_bin(${program} "${mr_alllibs}" src/mr/mrmain3d cc)
endforeach(program)

# Install libraries
install(TARGETS ${mr_libs} DESTINATION lib EXPORT mr_libs)
install(EXPORT mr_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${mr_targets1} DESTINATION bin)
install(TARGETS ${mr_targets2} DESTINATION bin)
install(TARGETS ${mr_targets3} DESTINATION bin)

endif()

# -- #
# MC #
# -- #

if(BUILD_MC)
 # Build standard Sparse2D libs
set(mc_libs mc1d mc2d)
foreach(library ${mc_libs})
  build_lib1(${library} src/mc cc h)
endforeach()

set(mc_alllibs ${mc_libs} ${mr_alllibs})

# Find all Sparse2D targets
find_targets(mc_targets1 src/mc/mcmain1d cc)
find_targets(mc_targets2 src/mc/mcmain2d cc)

# Build binaries
foreach(program ${mc_targets1})
  build_bin(${program} "${mc_alllibs}" src/mc/mcmain1d cc)
endforeach(program)
foreach(program ${mc_targets2})
  build_bin(${program} "${mc_alllibs}" src/mc/mcmain2d cc)
endforeach(program)
 
# Install libraries
install(TARGETS ${mc_libs} DESTINATION lib EXPORT mc_libs)
install(EXPORT mc_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${mc_targets1} DESTINATION bin)
install(TARGETS ${mc_targets2} DESTINATION bin)
endif()

# --- #
# MGA #
# --- #

if(BUILD_MGA)

 # Build standard Sparse2D libs
set(mga_libs mga1d mga3d)
foreach(library ${mga_libs})
  build_lib1(${library} src/mga cc h)
endforeach()

set(mga_alllibs ${mga_libs} ${mc_alllibs})

# Find all Sparse2D targets
find_targets(mga_targets1 src/mga/mgamain1d cc)
find_targets(mga_targets2 src/mga/mgamain2d cc)
find_targets(mga_targets3 src/mga/mgamain3d cc)

# Build binaries
foreach(program ${mga_targets1})
  build_bin(${program} "${mga_alllibs}" src/mga/mgamain1d cc)
endforeach(program)
foreach(program ${mga_targets2})
  build_bin(${program} "${mga_alllibs}" src/mga/mgamain2d cc)
endforeach(program)
foreach(program ${mga_targets3})
  build_bin(${program} "${mga_alllibs}" src/mga/mgamain3d cc)
endforeach(program)

# Install libraries
install(TARGETS ${mga_libs} DESTINATION lib EXPORT mga_libs)
install(EXPORT mga_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${mga_targets1} DESTINATION bin)
install(TARGETS ${mga_targets2} DESTINATION bin)
install(TARGETS ${mga_targets3} DESTINATION bin)



endif()

# ---- #
# MWIR #
# ---- #

if(BUILD_MWIR)

# Build standard Sparse2D libs
set(mwir_libs mwir)
foreach(library ${mwir_libs})
  build_lib1(${library} src/mwir cc h)
endforeach()

set(mwir_alllibs ${mwir_libs} ${mc_alllibs})

# Find all Sparse2D targets
find_targets(mwir_targets1 src/mwir/mwirmain cc)
 
# Build binaries
foreach(program ${mwir_targets1})
  build_bin(${program} "${mwir_alllibs}" src/mwir/mwirmain cc)
endforeach(program)
 

# Install libraries
install(TARGETS ${mwir_libs} DESTINATION lib EXPORT mwir_libs)
install(EXPORT mwir_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${mwir_targets1} DESTINATION bin)

endif()

# ------------------- #
# DICTIONARY LEARNING #
# ------------------- #

# NOTE: Need to to include GSL and Armadillo libraries
#       MSVST must be compiled as well

if(BUILD_DICLEARN)

  # Locate GSL using pkg-config
  # pkg_check_modules(GSL gsl)
  # find_pkg(GSL gsl)

  # Locate the Armadillo library
  find_package(Armadillo)
  include_directories(${ARMADILLO_INCLUDE_DIRS})

# Build standard Sparse2D libs
set(diclearn_libs diclearn)
foreach(library ${diclearn_libs})
  build_lib1(${library} src/diclearn cc h)
endforeach()

set(diclearn_alllibs ${diclearn_libs} ${mc_alllibs} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES})

# Find all Sparse2D targets
find_targets(diclearn_targets1 src/diclearn/diclearnmain cc)
 
# Build binaries
foreach(program ${diclearn_targets1})
  build_bin(${program} "${diclearn_alllibs}" src/diclearn/diclearnmain cc)
endforeach(program)
 

# Install libraries
install(TARGETS ${diclearn_libs} DESTINATION lib EXPORT diclearn_libs)
install(EXPORT diclearn_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${diclearn_targets1} DESTINATION bin)

endif()

# --- #
# MRS #
# --- #

if(BUILD_MRS)
  set(mrs_libs mrs mrsp)
  foreach(library ${mrs_libs})
  build_lib1(${library} src/mrs cc h)
endforeach()

# set(mrs_alllibs ${mrs_libs} ${mc_alllibs} ${HEALPIX_LIBRARIES}  ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES})
set(mrs_alllibs ${mrs_libs} ${mc_alllibs} ${HEALPIX_LIBRARIES})

# Find all Sparse2D targets
find_targets(mrs_targets1 src/mrs/mrsmain cc)
find_targets(mrs_targets2 src/mrs/mrspmain cc)

# Build binaries
foreach(program ${mrs_targets1})
  build_bin(${program} "${mrs_alllibs}" src/mrs/mrsmain cc)
endforeach(program)
 foreach(program ${mrs_targets2})
  build_bin(${program} "${mrs_alllibs}" src/mrs/mrspmain cc)
endforeach(program)

# Install libraries
install(TARGETS ${mrs_libs} DESTINATION lib EXPORT mrs_libs)
install(EXPORT mrs_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${mrs_targets1} DESTINATION bin)
install(TARGETS ${mrs_targets2} DESTINATION bin)

endif()

# -------- #
# Astro WL #
# -------- #

if(BUILD_ASTRO_WL)
 set(wl_libs wl divcurlwt)
  foreach(library ${wl_libs})
  build_lib1(${library} src/astro/weaklensing/ cc h)
endforeach()

set(wl_alllibs ${wl_libs} ${mc_alllibs})

# Find all Sparse2D targets
find_targets(wl_targets1 src/astro/weaklensing/mainwl cc)

# Build binaries
foreach(program ${wl_targets1})
  build_bin(${program} "${wl_alllibs}" src/astro/weaklensing/mainwl cc)
endforeach(program)

# Install libraries
install(TARGETS ${wl_libs} DESTINATION lib EXPORT wl_libs)
install(EXPORT wl_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${wl_targets1} DESTINATION bin)
 
#### second library using healpix
if(BUILD_MRS)
  set(wl_libs2 wls)
  foreach(library ${wl_libs2})
  build_lib1(${library} src/astro/weaklensing/ cc h)
endforeach()

set(wl_alllibs2 ${wl_libs2} ${wl_libs} ${mrs_alllibs} ${HEALPIX_LIBRARIES})

# Find all Sparse2D targets
find_targets(wl_targets2 src/astro/weaklensing/mainwls cc)

# Build binaries
 foreach(program ${wl_targets2})
  build_bin(${program} "${wl_alllibs2}" src/astro/weaklensing/mainwls cc)
endforeach(program)

# Install libraries
install(TARGETS ${wl_libs2} DESTINATION lib EXPORT wl_libs2)
install(EXPORT wl_libs2 DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${wl_targets2} DESTINATION bin)
endif()
endif()


# --------- #
# Astro Gal #
# --------- #

if(BUILD_ASTRO_GAL)

 set(gal_libs gal)
foreach(library ${gal_libs})
  build_lib1(${library} src/astro/galaxies cc h)
endforeach()

set(gal_alllibs ${gal_libs} ${mc_alllibs})

# Find all Sparse2D targets
find_targets(gal_targets1 src/astro/galaxies/galmain cc)
 
# Build binaries
foreach(program ${gal_targets1})
  build_bin(${program} "${gal_alllibs}" src/astro/galaxies/galmain cc)
endforeach(program)

# Install libraries
install(TARGETS ${gal_libs} DESTINATION lib EXPORT gal_libs)
install(EXPORT gal_libs DESTINATION lib/cmake)

# Install binaraies
install(TARGETS ${gal_targets1} DESTINATION bin)

endif()


# --------------- #
# Python Bindings #
# --------------- #

if(BUILD_BINDING)
  find_package(PythonInterp REQUIRED)
  if(NOT DEFINED PYTHON_LIBRARIES OR NOT DEFINED PYTHON_INCLUDE_DIRS)
    find_package(PythonLibs REQUIRED)
  else()
    message(STATUS "Using Python Lib: ${PYTHON_LIBRARIES}")
    message(STATUS "Using Python Inc: ${PYTHON_INCLUDE_DIRS}")
  endif()
  include_directories(${PYTHON_INCLUDE_DIRS})
  link_directories(${PYTHON_LIBRARY_DIRS})
  
  # Define binding targets
  set(binding_targets pymrs pysparse)
  
  # Build and link the module
  foreach(program ${binding_targets})

    add_library(${program} SHARED src/python_binding/${program}.cpp)

    if(APPLE)
      target_link_libraries(${program}  ${OpenMP_CXX_LIBRARIES} ${mrs_alllibs})
    else(APPLE)
      target_link_libraries(${program} ${OpenMP_CXX_LIBRARIES} ${mrs_alllibs} ${cfitsio_LIBRARIES} ${PYTHON_LIBRARIES})
    endif(APPLE)

    # add_dependencies(${program} cfitsio)
    set_target_properties(${program} PROPERTIES SUFFIX .so)
    set_target_properties(${program} PROPERTIES PREFIX "")

    if(APPLE)
      set_target_properties(${program} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
    endif(APPLE)

  endforeach(program)
  
  # set(pybind_deps ${mrs_alllibs} ${HEALPIX_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
  # Build bindings
  # build_pybind_list("${binding_targets}"  #List of bindings
  #  "${pybind_deps}"      #List of binding depdendencies
  #  src/python_binding        #Path to bindings
  #  cpp)                   #Source extension
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
