# -------------- #
# Fetch Pybind11 #
# -------------- #

# Set Pybind11 Version
set(Pybind11Version 2.13.6)
set(Pybind11SHA256 e08cb87f4773da97fa7b5f035de8763abc656d87d5773e62f6da0587d1f0ec20)

# Set Pybdin11 fetch location
FetchContent_Declare(
    pybind11
    URL  https://github.com/pybind/pybind11/archive/refs/tags/v${Pybind11Version}.tar.gz
    URL_HASH  SHA256=${Pybind11SHA256}
)

# Fetch and source Pybind11
FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_MakeAvailable(pybind11)
endif()