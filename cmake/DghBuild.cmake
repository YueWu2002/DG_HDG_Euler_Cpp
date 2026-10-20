include_guard(GLOBAL)
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    find_path(DGH_EIGEN_INCLUDE_DIR Eigen/Core
        PATHS "${CMAKE_CURRENT_LIST_DIR}/../.deps/eigen-3.4.0" REQUIRED)
    add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${DGH_EIGEN_INCLUDE_DIR}")
endif()
set(DGH_BLAS_BACKEND "AUTO" CACHE STRING "Eigen backend: AUTO, NONE, MKL, BLAS, or OPENBLAS")
set_property(CACHE DGH_BLAS_BACKEND PROPERTY STRINGS AUTO NONE MKL BLAS OPENBLAS)
string(TOUPPER "${DGH_BLAS_BACKEND}" DGH_BLAS_BACKEND)
if(NOT DGH_BLAS_BACKEND MATCHES "^(AUTO|NONE|MKL|BLAS|OPENBLAS)$")
    message(FATAL_ERROR "DGH_BLAS_BACKEND must be AUTO, NONE, MKL, BLAS, or OPENBLAS")
endif()
add_library(dgh_build_options INTERFACE)
target_link_libraries(dgh_build_options INTERFACE Eigen3::Eigen)
target_compile_features(dgh_build_options INTERFACE cxx_std_17)
# Debug keeps Eigen assertions and allocation guards; Release disables Eigen assertions.
target_compile_definitions(dgh_build_options INTERFACE
    EIGEN_RUNTIME_NO_MALLOC EIGEN_DONT_PARALLELIZE EIGEN_STACK_ALLOCATION_LIMIT=524288
    $<$<CONFIG:Release>:EIGEN_NO_DEBUG>)
# Opt in to tuning instructions for the machine doing the build. Keep this
# disabled by default so binaries remain portable across CPU generations.
option(DGH_NATIVE_ARCH "Tune generated instructions for the build host CPU" OFF)
if(DGH_NATIVE_ARCH)
    if(CMAKE_CXX_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
        target_compile_options(dgh_build_options INTERFACE -march=native)
        message(STATUS "Native CPU optimization: ON (-march=native)")
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
        target_compile_options(dgh_build_options INTERFACE -march=native)
        message(STATUS "Native CPU optimization: ON (-march=native)")
    else()
        message(FATAL_ERROR
            "DGH_NATIVE_ARCH=ON is unsupported for C++ compiler ${CMAKE_CXX_COMPILER_ID}")
    endif()
else()
    message(STATUS "Native CPU optimization: OFF (portable compiler target)")
endif()
# IntelLLVM (icpx) and classic Intel (icpc) are supported CMake C++ IDs.
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND NOT MSVC)
    target_compile_options(dgh_build_options INTERFACE -Wall $<$<CONFIG:Debug>:-UNDEBUG>)
elseif(MSVC)
    target_compile_options(dgh_build_options INTERFACE /W4 $<$<CONFIG:Debug>:/UNDEBUG>)
else()
    target_compile_options(dgh_build_options INTERFACE -Wall -Wextra -Wpedantic $<$<CONFIG:Debug>:-UNDEBUG>)
endif()
# Backend modules use the options target defined above. Their order matters:
# LAPACK reuses the selected BLAS/MKL libraries, as in the original project.
include("${CMAKE_CURRENT_LIST_DIR}/DghBlas.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/DghLapack.cmake")
# Eigen handles block LU and the direct reference. No external sparse solver.
add_library(dgh_fastgl STATIC "${CMAKE_CURRENT_LIST_DIR}/../vendor/fastgl/fastgl.cpp")
target_compile_features(dgh_fastgl PRIVATE cxx_std_17)
target_include_directories(dgh_fastgl PUBLIC "${CMAKE_CURRENT_LIST_DIR}/../vendor/fastgl")

function(dgh_add_solver target)
    add_executable(${target} ${ARGN})
    target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(${target} PRIVATE dgh_build_options dgh_fastgl)
endfunction()

function(dgh_add_common target)
    target_include_directories(${target} PUBLIC "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../common/include")
    target_sources(${target} PRIVATE
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../common/src/linalg/block_graph.cpp"
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../common/src/linalg/block_ilu.cpp"
        "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../common/src/linalg/gmres.cpp")
endfunction()

# SIMD directives without enabling threaded OpenMP. Match the reference build.
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(dgh_build_options INTERFACE -qopenmp-simd)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(dgh_build_options INTERFACE -fopenmp-simd)
endif()
include(CheckIPOSupported)
check_ipo_supported(RESULT DGH_IPO_SUPPORTED OUTPUT DGH_IPO_ERROR)
if(DGH_IPO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()
