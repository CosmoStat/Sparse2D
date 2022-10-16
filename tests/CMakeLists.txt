# ---------- #
# Unit tests #
# ---------- #

# The module contains build commands for the unit tests. Any new or modified
# files should be automatically detected and handled. Note that the build
# functions (i.e. build_tests) are defined in cmake/BuildTests.cmake.

# Build unit tests
build_tests("${test_libs}")
