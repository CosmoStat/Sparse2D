# Locate pkgconfig package
function(find_pkg package module)
  if(${package}_LIBRARIES STREQUAL "" OR NOT DEFINED ${package}_LIBRARIES OR
     ${package}_LIBRARY_DIRS STREQUAL "" OR NOT DEFINED ${package}_LIBRARY_DIRS OR
     ${package}_INCLUDE_DIRS STREQUAL "" OR NOT DEFINED ${package}_INCLUDE_DIRS)
    pkg_check_modules(${package} REQUIRED ${module})
  endif()
  include_directories(${${package}_INCLUDE_DIRS})
  link_directories(${${package}_LIBRARY_DIRS})
  if(${package}_FOUND)
    message(STATUS "  include dirs: ${${package}_INCLUDE_DIRS}")
    message(STATUS "  lib dirs: ${${package}_LIBRARY_DIRS}")
    message(STATUS "  libs: ${${package}_LIBRARIES}")
  endif()
endfunction()

# Extract target names from source files
function(find_targets targets target_path ext)
  file(GLOB src_targets "${PROJECT_SOURCE_DIR}/${target_path}/*.${ext}")
  list(TRANSFORM src_targets REPLACE "${PROJECT_SOURCE_DIR}/${target_path}/" "")
  list(TRANSFORM src_targets REPLACE ".${ext}" "")
  set(${targets} "${src_targets}" PARENT_SCOPE)
endfunction()

# Build library
function(build_lib library)
  file(GLOB src_${library} "${PROJECT_SOURCE_DIR}/src/lib${library}/*.cc")
  file(GLOB inc_${library} "${PROJECT_SOURCE_DIR}/src/lib${library}/*.h")
  include_directories("${PROJECT_SOURCE_DIR}/src/lib${library}")
  add_library(${library} STATIC ${src_${library}})
  target_link_libraries(${library} ${CFITSIO_LIBRARIES} ${fftw_lib_list})
  if(BUILD_FFTW)
    add_dependencies(${library} fftw)
  endif()
  INSTALL(FILES ${inc_${library}} DESTINATION include/sparse2d)
endfunction()

# Build library
function(build_lib1 library lib_path src_ext head_ext)
  file(GLOB src_${library} "${PROJECT_SOURCE_DIR}/${lib_path}/lib${library}/*.${src_ext}")
  file(GLOB inc_${library} "${PROJECT_SOURCE_DIR}/${lib_path}/lib${library}/*.${head_ext}")
  include_directories("${PROJECT_SOURCE_DIR}/${lib_path}/lib${library}")
  add_library(${library} STATIC ${src_${library}})
  target_link_libraries(${library} ${CFITSIO_LIBRARIES} ${fftw_lib_list})
  if(BUILD_FFTW)
    add_dependencies(${library} fftw)
  endif()
  INSTALL(FILES ${inc_${library}} DESTINATION include/sparse2d)
endfunction()


# Build binary
function(build_bin program libs target_path ext)
  add_executable(${program} "${PROJECT_SOURCE_DIR}/${target_path}/${program}.${ext}")
  target_link_libraries(${program} ${CFITSIO_LIBRARIES} ${fftw_lib_list} ${libs})
endfunction()

