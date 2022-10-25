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

  # Seach for the package module(s)
  pkg_check_modules(${package} QUIET "${module}")

  # Build from source or report package variables found
  if(BUILD_${package})

    message(STATUS "${package} will be built from source.")
    include("Build${package}")

  elseif(${package}_FOUND)

    message(STATUS "Found ${package}: version ${${package}_VERSION}")
    message(VERBOSE "  ${package}_INCLUDE_DIRS: ${${package}_INCLUDE_DIRS}")
    message(VERBOSE "  ${package}_LIBRARY_DIRS: ${${package}_LIBRARY_DIRS}")
    message(VERBOSE "  ${package}_LIBRARIES: ${${package}_LIBRARIES}")

  else()

    message(
      FATAL_ERROR "${package} not found! Install this package or set \
      BUILD_DEPS=ON. Note that some depenencies will take a long time to \
      build."
    )

  endif()

  # Include paths to package header and library files
  include_directories(${${package}_INCLUDE_DIRS})
  link_directories(${${package}_LIBRARY_DIRS})

endfunction()
