# --------------------- #
# Build python bindings #
# --------------------- #

# This module includes several functions for building Python bindings to
# Sparse2D C++ routines. Note that the build_pybind_target function is only
# meant to be run internally. The build_pybind simply requires a list of
# libraries that the bindings depends on (e.g. cfitsio, sparse, etc.) as an
# input.

# Build Python bindings
function(build_pybind lib_deps_list build_all)

  # Set the full path to the bindings
  set(pybind_path "${PROJECT_SOURCE_DIR}/src/python_binding")

  # Find binding targets
  if(build_all)
    file(GLOB pybind_targets CONFIGURE_DEPENDS "${pybind_path}/*.c*")
  else()
    set(pybind_targets "${pybind_path}/pysparse.cpp")
  endif()

  # Loop through the list of targets
  foreach(pybind_target ${pybind_targets})
    # Build python binding
    build_pybind_target(${pybind_target} "${lib_deps_list}")
  endforeach()

endfunction()

# Build Python binding target
function(build_pybind_target target libs)

  # Extract the target name
  get_filename_component(tarname ${target} NAME_WLE)

  # Add the target library
  add_library(${tarname} SHARED ${target})

  # Link libraries to target
  target_link_libraries(${tarname} "${libs}" OpenMP::OpenMP_CXX pybind11::headers)

  # Set system dependend properties
  if(APPLE)
    set_target_properties(${tarname} PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
  else(APPLE)
    target_link_libraries(${tarname} ${Python_LIBRARIES})
  endif(APPLE)

  # Set target prefix and suffix properties
  pybind11_extension(${tarname})
  pybind11_strip(${tarname})

  # Set the install path for bindings
  install(TARGETS ${tarname} DESTINATION ${PYBIND_INSTALL_PATH})

endfunction()
