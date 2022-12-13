# ------------- #
# Build CfitsIO #
# ------------- #

# Set CFITSIO Version
set(CFITSIOVersion 4.1.0)
set(cfitsioSHA256 b367c695d2831958e7166921c3b356d5dfa51b1ecee505b97416ba39d1b6c17a)

# Set C/C++ compiler and flags
set(CFITSIO_COMPILE
  CC=${CMAKE_C_COMPILER}
  CXX=${CMAKE_CXX_COMPILER}
)

# Set CFITSIO configuration flags
set(CFITSIO_CONFIG_FLAGS
  --prefix=${CMAKE_BINARY_DIR}
  --enable-reentrant
)

# Download and build CFITSIO
ExternalProject_Add(cfitsio_build
    PREFIX cfitsio
    URL  http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-${CFITSIOVersion}.tar.gz
    URL_HASH  SHA256=${cfitsioSHA256}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ${CFITSIO_COMPILE} ./configure ${CFITSIO_CONFIG_FLAGS}
    BUILD_COMMAND make shared
    INSTALL_COMMAND make install
)

# Set CFITSIO variables
set(CFITSIO_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/lib/ CACHE INTERNAL "")
set(CFITSIO_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/include/ CACHE INTERNAL "")
set(CFITSIO_LIBRARIES -lcfitsio CACHE INTERNAL "")
