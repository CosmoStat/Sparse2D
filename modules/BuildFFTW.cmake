#=============#
# Build FFTW  #
#=============#

# Set FFTW Version
set(FFTWVersion 3.3.8)
set(FFTWHash 8aac833c943d8e90d51b697b27d4384d)

# Set FFTW libraries to use
set(fftw_lib_list "-lfftw3f_omp -lfftw3_omp -lfftw3f -lfftw3 -lm")

# Find FFTW
if(USE_FFTW AND NOT BUILD_FFTW)
  find_package(FFTW3 REQUIRED)
  message(STATUS "Checking for module 'fftw3'")
  message(STATUS "  Found fftw3, version ${FFTW3_VERSION}")
  set(FFTW_LD_FLAGS "-L ${FFTW3_LIBRARY_DIRS} ${fftw_lib_list}")
elseif(USE_FFTW AND BUILD_FFTW)
  set(FFTW_LD_FLAGS "-L ${MODULE_BUILD_DIR}/lib ${fftw_lib_list}")
endif()

# Set FFTW flags
if(USE_FFTW)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FFTW")
endif()

# Build FFTW
if(BUILD_FFTW)

  # Set C/C++ compiler and flags
  set(FFTW_COMPILE
      CC=${CMAKE_C_COMPILER}
      CXX=${CMAKE_CXX_COMPILER}
      ${BigMac_FFTW})

  # Set FFTW configuration flags
  set(FFTW_CONFIG_FLAGS
      --prefix=${MODULE_BUILD_DIR}
      --enable-maintainer-mode
      --enable-shared
      --enable-threads
      --enable-sse2
      --enable-openmp)

  # Download and build FFTW
  ExternalProject_Add(fftw
      URL               http://www.fftw.org/fftw-${FFTWVersion}.tar.gz
      URL_HASH          MD5=${FFTWHash}
      SOURCE_DIR        ${MODULE_BUILD_DIR}fftw
      BUILD_IN_SOURCE   1
      CONFIGURE_COMMAND ./configure ${FFTW_COMPILE} ${FFTW_CONFIG_FLAGS}
      COMMAND           make && make install
      COMMAND           ./configure ${FFTW_COMPILE} ${FFTW_CONFIG_FLAGS} --enable-single
      COMMAND           make && make install
      COMMAND           make clean
  )

endif()

message(STATUS "Use FFTW: ${USE_FFTW}")
message(STATUS "FFTW Build: ${BUILD_FFTW}")
