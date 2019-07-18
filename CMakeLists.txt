#################
# Sparse2D 2017 #
#################

cmake_minimum_required(VERSION 3.0.0)
include(ExternalProject)
include(FindPkgConfig)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/modules/")

project(sparse2d)

  # Set build type
  if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release")
  endif()
  message(STATUS "CMake Build Type: ${CMAKE_BUILD_TYPE}")

  # Set the default installation path
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/sparse2d"
        CACHE PATH "Setting new default path" FORCE)
  endif()
  message(STATUS "Installing to: ${CMAKE_INSTALL_PREFIX}")

  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
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

  # Locate CFITSIO using pkg-config or use command line arguments to configure
  # CFITSIO
  if(CFITSIO_LIBRARIES STREQUAL "" OR NOT DEFINED CFITSIO_LIBRARIES OR
     CFITSIO_LIBRARY_DIRS STREQUAL "" OR NOT DEFINED CFITSIO_LIBRARY_DIRS OR
     CFITSIO_INCLUDE_DIRS STREQUAL "" OR NOT DEFINED CFITSIO_INCLUDE_DIRS)
    pkg_check_modules(CFITSIO REQUIRED cfitsio)
  else()
    message(STATUS "Use manually configured cfitsio")
    message(STATUS "  includes: ${CFITSIO_INCLUDE_DIRS}")
    message(STATUS "  libs: ${CFITSIO_LIBRARY_DIRS}")
    message(STATUS "  flags: ${CFITSIO_LIBRARIES}")
  endif()
  include_directories(${CFITSIO_INCLUDE_DIRS})
  link_directories(${CFITSIO_LIBRARY_DIRS})

  # Build location for external modules
  set(MODULE_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/module_build/")

  # Optional build of FFTW3
  option(USE_FFTW "Use FFTW for Fourier transformation" OFF)
  include(BuildFFTW3)

  # Optional build of NFFT
  option(COMPILE_NFFT "Get and compile NFFT" OFF)
  include(BuildNFFT)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_DISP_IO -fPIC ${OpenMP_CXX_FLAGS} -Wno-write-strings ${FFTW_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g0 -O2 -fomit-frame-pointer -DNDEBUG")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
  message(STATUS "Using CXX Flags for ${CMAKE_SYSTEM}")

  # Build tools library
  FILE(GLOB src_lib1  "${PROJECT_SOURCE_DIR}/src/libtools/*.cc")
  include_directories("${PROJECT_SOURCE_DIR}/src/libtools")
  add_library(tools STATIC ${src_lib1})
  target_link_libraries(tools ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
  if(USE_FFTW)
    add_dependencies(tools fftw3)
  endif(USE_FFTW)

  # Build sparse1d library
  FILE(GLOB src_lib2 "${PROJECT_SOURCE_DIR}/src/libsparse1d/*.cc")
  include_directories("${PROJECT_SOURCE_DIR}/src/libsparse1d")
  add_library(sparse1d STATIC ${src_lib2})
  target_link_libraries(sparse1d ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
  if(USE_FFTW)
    add_dependencies(sparse1d fftw3)
  endif(USE_FFTW)

  # Build sparse2d library
  FILE(GLOB src_lib2 "${PROJECT_SOURCE_DIR}/src/libsparse2d/*.cc")
  include_directories("${PROJECT_SOURCE_DIR}/src/libsparse2d")
  add_library(sparse2d STATIC ${src_lib2})
  target_link_libraries(sparse2d ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
  if(USE_FFTW)
    add_dependencies(sparse2d fftw3)
  endif(USE_FFTW)

  option(SPARSE3D "Build Sparse3D library" ON)
  if(SPARSE3D)
    # Build sparse3d library
    FILE(GLOB src_lib2 "${PROJECT_SOURCE_DIR}/src/libsparse3d/*.cc")
    include_directories("${PROJECT_SOURCE_DIR}/src/libsparse3d")
    add_library(sparse3d STATIC ${src_lib2})
    target_link_libraries(sparse3d ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
    if(USE_FFTW)
      add_dependencies(sparse3d fftw3)
    endif(USE_FFTW)
  endif(SPARSE3D)
  message(STATUS "Sparse3D Build: ${SPARSE3D}")

  # Build mga2d library
  FILE(GLOB src_mgalib2 "${PROJECT_SOURCE_DIR}/src/libmga2d/*.cc")
  include_directories("${PROJECT_SOURCE_DIR}/src/libmga2d")
  add_library(mga2d STATIC ${src_mgalib2})
  target_link_libraries(mga2d ${CFITSIO_LIBRARIES} ${FFTW_LD_FLAGS})
  if(USE_FFTW)
    add_dependencies(mga2d fftw3)
  endif(USE_FFTW)

  option(BUILD_MSVST "Build MSVST library" OFF)
  if(BUILD_MSVST)
    FILE(GLOB msvst_src_lib "${PROJECT_SOURCE_DIR}/src/libmsvst/*.cc")
    include_directories("${PROJECT_SOURCE_DIR}/src/libmsvst")
    add_library(msvst STATIC ${msvst_src_lib} )
  endif(BUILD_MSVST)
  message(STATUS "MSVST Build: ${BUILD_MSVST}")

  # Compile and link executables
  set(BINMR2D mr_transform mr_recons mr_filter mr_deconv mr_info cur_contrast
      cur_deconv cur_filter cur_stat cur_trans im1d_convert im1d_info im1d_stf
      im1d_tfreq mr1d_detect mr1d_trans im1d_der_loggamma im1d_period im1d_tend
      mr1d_filter)
  if(SPARSE3D)
    set(BINMR2D ${BINMR2D} mr3d_trans mr3d_filter mr3d_recons mr3d_stat
        mr2d1d_trans)
  endif(SPARSE3D)

  foreach(program ${BINMR2D})
    add_executable(${program} ${PROJECT_SOURCE_DIR}/src/${program}.cc)
    target_link_libraries(${program} mga2d sparse2d sparse1d tools)
    if(SPARSE3D)
      target_link_libraries(${program} mga2d sparse3d sparse2d sparse1d tools)
    endif(SPARSE3D)
  endforeach(program)

  if(BUILD_MSVST)
    set(MSVST msvst_uwt2d msvst_iwt2d msvst_2d1d msvst_iwt2d_coupled)
    foreach(program ${MSVST})
      add_executable(${program} ${PROJECT_SOURCE_DIR}/src/${program}.cc)
      target_link_libraries(${program} msvst tools)
    endforeach(program)
    FILE(GLOB inc_lib "${PROJECT_SOURCE_DIR}/src/libmsvst/*.h")
    INSTALL(FILES ${inc_lib} DESTINATION include/sparse2d)
    INSTALL(TARGETS ${MSVST} DESTINATION bin)
  endif(BUILD_MSVST)

  # Install headers
  FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/libsparse1d/*.h")
  INSTALL(FILES ${inc_lib} DESTINATION include/sparse2d)
  FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/libsparse2d/*.h")
  INSTALL(FILES ${inc_lib} DESTINATION include/sparse2d)
  FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/libmga2d/*.h")
  INSTALL(FILES ${inc_lib} DESTINATION include/sparse2d)
  FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/libtools/*.h")
  INSTALL(FILES ${inc_lib} DESTINATION include/sparse2d)
  if(SPARSE3D)
    FILE(GLOB inc_lib  "${PROJECT_SOURCE_DIR}/src/libsparse3d/*.h")
    INSTALL(FILES ${inc_lib} DESTINATION include/sparse2d)
  endif(SPARSE3D)

  # Install library
  INSTALL(TARGETS sparse1d sparse2d mga2d tools DESTINATION lib)
  if(SPARSE3D)
    INSTALL(TARGETS sparse3d DESTINATION lib)
  endif(SPARSE3D)

  # install  executables
  INSTALL(TARGETS ${BINMR2D} DESTINATION bin)

  # ---------- #
  # UNIT TESTS #
  # -----------#

  option(RUN_TESTS "Build and run unit tests" ON)
  if(RUN_TESTS)
    # Add unit tests
    enable_testing()
    set(UNIT_TESTS mr_transform_test)

    foreach(program ${UNIT_TESTS})
      add_executable(${program} ${PROJECT_SOURCE_DIR}/tests/${program}.cpp)
      target_link_libraries(${program} mga2d sparse2d sparse1d tools)
      if(SPARSE3D)
        target_link_libraries(${program} mga2d sparse3d sparse2d sparse1d tools)
      endif(SPARSE3D)
      add_test(${program} ${program})
    endforeach(program)
  endif(RUN_TESTS)
