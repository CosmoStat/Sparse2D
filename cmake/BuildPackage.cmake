# ------------- #
# Build package #
# ------------- #

# This module includes several functions for building Sparse2D subpackage
# libraries and main executable binaries. Note that most of the functions are
# are only meant to be used internally by other functions. The only function
# intended for external use is build_package, which expects the package name
# (e.g. msvst), a list of libraries that the packages depends on (e.g.
# cfitsio, sparse, etc.) and a variable name for the libraries created (e.g.
# msvst_libs) as inputs.

# Build a Sparse2D package with a list of depencencies
function(build_package package_name package_libs lib_deps_list)

  # Set the full path to the package
  set(package_path "${PROJECT_SOURCE_DIR}/src/${package_name}")

  # Find the library and main subdirectories
  file(GLOB lib_dir_list CONFIGURE_DEPENDS "${package_path}/lib*")
  file(GLOB main_dir_list CONFIGURE_DEPENDS "${package_path}/*main*")

  # Build the libraries
  build_lib_list("${lib_dir_list}")

  # Define all the libraries needed for the package build
  set(package_all_libs ${package_libs} ${lib_deps_list})

  # Build the mains
  build_main_list("${main_dir_list}" "${package_all_libs}")

endfunction()

# Build library
function(build_lib library)

  # Extract the library name
  get_filename_component(lib_dir_name ${library} NAME)
  string(REPLACE "lib" "" libname ${lib_dir_name})

  # Find the source and header files in the library
  file(GLOB src_${libname} CONFIGURE_DEPENDS "${library}/*.c*")
  file(GLOB inc_${libname} CONFIGURE_DEPENDS "${library}/*.h*")

  # Include the path to the header files
  include_directories("${library}")

  # Add the library
  add_library(${libname} STATIC ${src_${libname}})
  set_property(TARGET ${libname} PROPERTY POSITION_INDEPENDENT_CODE ON)
  if(BUILD_CFITSIO)
    add_dependencies(${libname} cfitsio_build)
  endif()

  # Set the install path for header files and libraries
  install(FILES ${inc_${libname}} DESTINATION include/sparse2d)
  install(TARGETS ${libname} DESTINATION lib)

endfunction()

# Build list of libraries
function(build_lib_list library_list)

  # Loop through the list of libraries
  foreach(library ${library_list})
    # Build library
    build_lib(${library})
  endforeach()

endfunction()

# Build binary
function(build_bin target libs)

  # Extract the target name
  get_filename_component(tarname ${target} NAME_WLE)

  # Add target executable
  add_executable(${tarname} "${target}")

  # Link libraries to target
  target_link_libraries(${tarname} ${libs} OpenMP::OpenMP_CXX)

  # Set the install target binaries
  install(TARGETS ${tarname} DESTINATION bin)

endfunction()

# Build list of binaries
function(build_bin_list target_list libs)

  # Loop through the list of targets
  foreach(target ${target_list})
    # Build binary
    build_bin(${target} "${libs}")
  endforeach()

endfunction()

# Build list of main directories
function(build_main_list main_list libs)

  foreach(main ${main_list})
    # Find main targets
    file(GLOB main_targets CONFIGURE_DEPENDS "${main}/*.c*")
    # Build the binaries
    build_bin_list("${main_targets}" "${libs}")
  endforeach()

endfunction()
