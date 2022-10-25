# ----------- #
# Build tests #
# ----------- #

# This module includes several functions for building Sparse2D unit tests

# Build a Sparse2D package with a list of depencencies
function(build_tests lib_deps_list)

  # Source the relevant headers
  foreach(libname ${lib_deps_list})
    # Include headers
    include_headers(${libname})
  endforeach()

  # Set the full path to the package
  set(test_path "${PROJECT_SOURCE_DIR}/tests")

  # Find test targets
  file(GLOB test_targets CONFIGURE_DEPENDS "${test_path}/*.c*")

  # Build the tests
  foreach(target ${test_targets})
    # Extract the target name
    get_filename_component(tarname ${target} NAME_WLE)
    # Add target executable
    add_executable(${tarname} "${target}")
    # Link libraries to target
    target_link_libraries(${tarname} ${lib_deps_list} OpenMP::OpenMP_CXX Catch2::Catch2WithMain)
    # Add unit test
    add_test(${tarname} ${tarname})
  endforeach()

endfunction()

# Incliude headers
function(include_headers libname)

  # Set the full path to the library directories
  set(lib_path "${PROJECT_SOURCE_DIR}/src/*/lib")

  # Find the library subdirectories
  file(GLOB lib_dir CONFIGURE_DEPENDS "${lib_path}${libname}")

  # Include the path to the header files
  include_directories(${lib_dir})

endfunction()
