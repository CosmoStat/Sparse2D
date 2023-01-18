# -------- #
# Find PKG #
# -------- #

# This module defines the find_pkg function that expects the package name and
# corresponding module name(s) as inputs. This function will automatically
# link the header and library directories for the package. If the
# BUILD_<package name> option is set the package will be built from source. A
# report of the package variables set for building can be seen by running CMake
# with --loglevel=VERBOSE.

# Locate pkgconfig package
function(find_pkg package module)

  # Build from source or find package
  if(BUILD_${package})

    # Include build package module
    message(STATUS "${package} will be built from source.")
    set(${package}_FOUND TRUE)
    include("Build${package}")

  else()

    # Seach for the package module(s)
    pkg_check_modules(${package} QUIET "${module}")

  endif()

  # Make sure the package was found
  if(NOT ${package}_FOUND)

    message(
      FATAL_ERROR "${package} not found! Install this package or set \
      BUILD_${package}=ON. Note that some depenencies will take a long time \
      to build."
    )

  endif()

  # Report package contents
  message(VERBOSE "  Package: ${package}")
  message(VERBOSE "  ${package}_INCLUDE_DIRS: ${${package}_INCLUDE_DIRS}")
  message(VERBOSE "  ${package}_LIBRARY_DIRS: ${${package}_LIBRARY_DIRS}")
  message(VERBOSE "  ${package}_LIBRARIES: ${${package}_LIBRARIES}")

  # Include paths to package header and library files
  include_directories(${${package}_INCLUDE_DIRS})
  link_directories(${${package}_LIBRARY_DIRS})

endfunction()
