# ------------ #
# Build FFTW3  #
# ------------ #

# Set FFTW Version
set(FFTWVersion 3.3.10)
set(FFTWHash 8ccbf6a5ea78a16dbc3e1306e234cc5c)

# Set C/C++ compiler and flags
set(FFTW_COMPILE
  CC=${CMAKE_C_COMPILER}
  CXX=${CMAKE_CXX_COMPILER}
  ${BigMac_FFTW}
)

# Set FFTW3 configuration flags
set(FFTW_CONFIG_FLAGS
  --prefix=${MODULE_BUILD_DIR}
  --enable-maintainer-mode
  --enable-shared
  --enable-threads
  --enable-sse2
  --enable-openmp
)

# Download and build FFTW3
ExternalProject_Add(fftw3_build
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

# Set FFTW3 variables
set(FFTW3_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/lib/)
set(FFTW3_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/include/)
set(FFTW3_LIBRARIES -lfftw3 -lfftw3f)
