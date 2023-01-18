# -------------- #
# Build HEALPix  #
# -------------- #

# Set HEALPix Version
set(HEALPixVersion 3.82)
set(HEALPixFileDate 2022Jul28)
set(HEALPixHash f71ffbccd3009f5e602fb48b4fa134d5)

# Set C/C++ compiler and flags
set(HEALPIX_COMPILE
  CC=${CMAKE_C_COMPILER}
  CXX=${CMAKE_CXX_COMPILER}
  FITSDIR=${CFITSIO_LIBRARY_DIRS}
  FITSINC=${CFITSIO_INCLUDE_DIRS}
)

# Set HEALPix configuration flags
set(HEALPIX_CONFIG_FLAGS
  -L
  --auto=profile,sharp,c,cxx
)

# Download and build HEALPix
ExternalProject_Add(healpix_build
    URL               https://sourceforge.net/projects/healpix/files/Healpix_${HEALPixVersion}/Healpix_${HEALPixVersion}_${HEALPixFileDate}.tar.gz
    URL_HASH          MD5=${HEALPixHash}
    SOURCE_DIR        ${CMAKE_BINARY_DIR}/healpix
    BUILD_IN_SOURCE   1
    CONFIGURE_COMMAND ${HEALPIX_COMPILE} ./configure ${HEALPIX_CONFIG_FLAGS}
    BUILD_COMMAND     make
    COMMAND           make clean
    INSTALL_COMMAND   cmake -E echo "Skipping install step."
    DEPENDS           cfitsio_build
)

# Set HEALPix variables
set(HEALPIX_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/healpix/lib/ CACHE INTERNAL "")
set(HEALPIX_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/healpix/include/healpix_cxx/ CACHE INTERNAL "")
set(HEALPIX_LIBRARIES -lhealpix_cxx CACHE INTERNAL "")
