cmake_minimum_required (VERSION 3.5.1)

project(demo CXX)
set(CMAKE_CXX_STANDARD 17)
option( BUILD_STATIC "Set to ON to include static versions of the library" OFF)

find_package(OpenFHE CONFIG REQUIRED)
if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif ()

# Remove any existing -Werror flag from OpenFHE_CXX_FLAGS
string(REPLACE "-Werror" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Add -Wno-error to ensure warnings are not treated as errors
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")

set( CMAKE_CXX_FLAGS ${OpenFHE_CXX_FLAGS} )

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${OpenFHE_INCLUDE} )
include_directories( ${OpenFHE_INCLUDE}/third-party/include )
include_directories( ${OpenFHE_INCLUDE}/core )
include_directories( ${OpenFHE_INCLUDE}/pke )
include_directories( ${OpenFHE_INCLUDE}/binfhe )
include_directories(degrees PUBLIC /usr/local/include /usr/local/include/NTL)

link_directories( ${OpenFHE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} )
    link_libraries(ntl gmp pthread m)
endif()

# Set the MPFR_DIR to the directory containing MPFRConfig.cmake or mpfr-config.cmake
# set(MPFR_DIR "/usr/include/boost/math/bindings/mpfr.hpp")

# Append the cmake/ directory to the CMAKE_MODULE_PATH
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set custom NTL path (for local user space installation)
set(NTL_LIBRARY_DIR /home/s2086034/.local/lib)
set(NTL_INCLUDE_DIR /home/s2086034/.local/include)

# Find the NTL library and headers in local user space
find_library(NTL_LIBRARY NAMES ntl PATHS ${NTL_LIBRARY_DIR} REQUIRED)
find_path(NTL_INCLUDE_DIR NAMES NTL/ZZ.h PATHS ${NTL_INCLUDE_DIR} REQUIRED)

# Include directories for NTL
include_directories(${NTL_INCLUDE_DIR})

# link_directories(${MPFR_LIBRARY_DIRS} ${GMP_LIBRARY_DIRS})

### ADD YOUR EXECUTABLE(s) HERE
### add_executable( EXECUTABLE-NAME SOURCES )cd
###
### EXAMPLE:
### add_executable( test demo-simple-example.cpp )

set(SOURCES
    src/utils-basics.cpp
    src/utils-eval.cpp
    src/utils-matrices.cpp
    src/utils-ptxt.cpp
    src/ranking.cpp
    src/ranking_compos.cpp
    src/minimum.cpp
    src/sorting.cpp
    src/demo.cpp
    src/func.cpp
    src/optimized_degrees.cpp
    src/utils-remez.cpp
    # src/MinimaxComp_degrees/run_compare.cpp
    # src/MinimaxComp_degrees/func.cpp
    # src/MinimaxComp_degrees/optimized_degrees.cpp
    # src/remez_cpp.cpp
    )

add_executable( program ${SOURCES} )

# Link NTL and optionally GMP to your target
target_link_libraries(program ${NTL_LIBRARY})

# Optionally, add target-specific flags
target_compile_options(program PRIVATE -Wno-error)
