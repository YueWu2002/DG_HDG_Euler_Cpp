# Optional Eigen dense-decomposition backend. Eigen's own headers implement
# the LAPACKE adapter; solver .cpp/.hpp files keep using ordinary Eigen LU.
# Linking LAPACK alone is insufficient: LAPACKE supplies its C entry points.
set(DGH_LAPACK_BACKEND "AUTO" CACHE STRING
    "Eigen dense decomposition backend: AUTO, NONE, MKL, or LAPACKE")
set_property(CACHE DGH_LAPACK_BACKEND PROPERTY STRINGS AUTO NONE MKL LAPACKE)
string(TOUPPER "${DGH_LAPACK_BACKEND}" DGH_LAPACK_BACKEND)
if(NOT DGH_LAPACK_BACKEND MATCHES "^(AUTO|NONE|MKL|LAPACKE)$")
    message(FATAL_ERROR "DGH_LAPACK_BACKEND must be AUTO, NONE, MKL, or LAPACKE")
endif()

include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

# Exercise Eigen's row/column-major LU and LLT APIs used by the solver and
# self-tests, rather than testing for one bare LAPACK symbol. This catches
# incomplete LAPACKE libraries and incompatible links before compiling.
# The CMake check compiles and links; it does not execute the numerical solve.
function(dgh_probe_eigen_lapack result)
    cmake_push_check_state(RESET)
    set(CMAKE_REQUIRED_DEFINITIONS -DEIGEN_USE_LAPACKE)
    set(CMAKE_REQUIRED_LIBRARIES Eigen3::Eigen ${ARGN})
    unset(${result} CACHE)
    check_cxx_source_compiles("#include <Eigen/Dense>
        #ifndef EIGEN_USE_LAPACKE
        #error Eigen disabled the requested LAPACKE adapter
        #endif
        int main() {
          Eigen::MatrixXd a(2,2);
          a << 3.0, 1.0, 1.0, 2.0;
          Eigen::PartialPivLU<Eigen::MatrixXd> lu(a);
          Eigen::VectorXd b(2); b << 1.0, 2.0;
          Eigen::VectorXd x = lu.solve(b);
          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> row(a);
          Eigen::PartialPivLU<decltype(row)> row_lu(row);
          Eigen::VectorXd z = row_lu.solve(b);
          Eigen::LLT<Eigen::MatrixXd> chol(a);
          Eigen::VectorXd y = chol.solve(b);
          return x.size() != 2 || y.size() != 2 || z.size() != 2;
        }" ${result})
    set(${result} "${${result}}" PARENT_SCOPE)
    cmake_pop_check_state()
endfunction()

set(DGH_RESOLVED_LAPACK_BACKEND "NONE")
if(DGH_LAPACK_BACKEND STREQUAL "MKL" AND NOT DGH_RESOLVED_BACKEND STREQUAL "MKL")
    message(FATAL_ERROR
        "DGH_LAPACK_BACKEND=MKL requires -DDGH_BLAS_BACKEND=MKL (or AUTO selecting MKL)")
endif()
if(DGH_LAPACK_BACKEND STREQUAL "LAPACKE" AND DGH_RESOLVED_BACKEND STREQUAL "MKL")
    message(FATAL_ERROR
        "MKL BLAS is selected; use DGH_LAPACK_BACKEND=MKL to avoid mixing LAPACK vendors")
endif()

if((DGH_LAPACK_BACKEND STREQUAL "AUTO" OR DGH_LAPACK_BACKEND STREQUAL "MKL") AND
   DGH_RESOLVED_BACKEND STREQUAL "MKL")
    # MKL::MKL supplies both the F77 LAPACK and LAPACKE LP64 entry points.
    dgh_probe_eigen_lapack(DGH_MKL_LAPACKE_LINKS ${DGH_BLAS_TARGET})
    if(DGH_MKL_LAPACKE_LINKS)
        set(DGH_RESOLVED_LAPACK_BACKEND "MKL")
        set(DGH_LAPACK_TARGETS ${DGH_BLAS_TARGET})
    elseif(DGH_LAPACK_BACKEND STREQUAL "MKL")
        message(FATAL_ERROR "MKL was found, but its Eigen/LAPACKE compile-link check failed")
    else()
        message(STATUS "MKL LAPACKE unavailable; using Eigen native dense decompositions")
    endif()
elseif(DGH_LAPACK_BACKEND STREQUAL "LAPACKE" OR
       (DGH_LAPACK_BACKEND STREQUAL "AUTO" AND NOT DGH_RESOLVED_BACKEND STREQUAL "MKL"))
    # Eigen supplies lapacke.h itself. The installed liblapacke and LAPACK
    # (with its BLAS dependency) must still provide the actual functions.
    find_library(DGH_LAPACKE_LIBRARY NAMES lapacke)
    if(DGH_LAPACKE_LIBRARY)
        find_package(LAPACK QUIET)
    endif()
    if(DGH_LAPACKE_LIBRARY AND LAPACK_FOUND)
        set(_dgh_lapack_libraries "${DGH_LAPACKE_LIBRARY}" LAPACK::LAPACK)
        if(DGH_BLAS_TARGET)
            list(APPEND _dgh_lapack_libraries ${DGH_BLAS_TARGET})
        endif()
        dgh_probe_eigen_lapack(DGH_LAPACKE_LINKS ${_dgh_lapack_libraries})
        if(DGH_LAPACKE_LINKS)
            set(DGH_RESOLVED_LAPACK_BACKEND "LAPACKE")
            set(DGH_LAPACK_TARGETS ${_dgh_lapack_libraries})
        endif()
        unset(_dgh_lapack_libraries)
    endif()
    if(DGH_LAPACK_BACKEND STREQUAL "LAPACKE" AND
       NOT DGH_RESOLVED_LAPACK_BACKEND STREQUAL "LAPACKE")
        message(FATAL_ERROR
            "LAPACKE requested, but a linkable liblapacke + LAPACK Eigen adapter was not found")
    endif()
endif()

if(NOT DGH_RESOLVED_LAPACK_BACKEND STREQUAL "NONE")
    target_compile_definitions(dgh_build_options INTERFACE EIGEN_USE_LAPACKE)
    target_link_libraries(dgh_build_options INTERFACE ${DGH_LAPACK_TARGETS})
    message(STATUS "Eigen LAPACK linked: YES (${DGH_RESOLVED_LAPACK_BACKEND})")
else()
    message(STATUS "Eigen LAPACK linked: NO (using Eigen native dense decompositions)")
endif()
target_compile_definitions(dgh_build_options INTERFACE
    DGH_DENSE_LAPACK_BACKEND="${DGH_RESOLVED_LAPACK_BACKEND}")
message(STATUS "Eigen LAPACK backend: ${DGH_RESOLVED_LAPACK_BACKEND} (requested ${DGH_LAPACK_BACKEND})")
