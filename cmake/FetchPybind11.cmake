# -------------- #
# Fetch Pybind11 #
# -------------- #

# Set Pybind11 Version
set(Pybind11Version 2.10.3)
set(Pybind11SHA256 5d8c4c5dda428d3a944ba3d2a5212cb988c2fae4670d58075a5a49075a6ca315)

# Set Pybdin11 fetch location
FetchContent_Declare(
    pybind11
    URL  https://github.com/pybind/pybind11/archive/refs/tags/v${Pybind11Version}.tar.gz
    URL_HASH  SHA256=${Pybind11SHA256}
)

# Fetch and source Pybind11
FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()