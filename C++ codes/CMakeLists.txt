cmake_minimum_required(VERSION 3.20)
project(BayesProject)

# Set C++ Standard to 17 (Recommended for modern Rcpp/Armadillo)
set(CMAKE_CXX_STANDARD 17)


# --- 1. Setup for R ---
if(APPLE)
    set(R_EXECUTABLE "/Library/Frameworks/R.framework/Resources/bin/R")
else()
    find_program(R_EXECUTABLE R)
endif()

if(NOT EXISTS "${R_EXECUTABLE}")
    message(FATAL_ERROR "Could not find R executable at ${R_EXECUTABLE}. Please check your installation.")
else()
    message(STATUS "R executable set to: ${R_EXECUTABLE}")
endif()

# Query R for its Home Directory (RHOME)
execute_process(
        COMMAND ${R_EXECUTABLE} RHOME
        OUTPUT_VARIABLE R_HOME
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "R Home directory: ${R_HOME}")

set(R_INCLUDE_DIRS "${R_HOME}/include")
set(R_LIBRARY_DIR "${R_HOME}/lib")


# --- 2. Find required R packages(Rcpp, Armadillo, BH) ---
# We define a macro to ask R where specific packages are installed.
# This ensures portability across different computers (Mac, Linux, Windows).

macro(find_r_package_include PACKAGE_NAME VAR_NAME)
    execute_process(
            COMMAND ${R_EXECUTABLE} --vanilla --slave -e "cat(system.file('include', package='${PACKAGE_NAME}'))"
            OUTPUT_VARIABLE ${VAR_NAME}
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if("${${VAR_NAME}}" STREQUAL "")
        message(FATAL_ERROR "Error: R package '${PACKAGE_NAME}' not found. Please run install.packages('${PACKAGE_NAME}') in RStudio.")
    else()
        message(STATUS "Found ${PACKAGE_NAME} include path: ${${VAR_NAME}}")
    endif()
endmacro()

# Find the specific paths for your dependencies
find_r_package_include("Rcpp" RCPP_INC)
find_r_package_include("RcppArmadillo" ARMA_INC)
find_r_package_include("BH" BH_INC)


# --- 2.1. Find BLAS & LAPACK (essential for Armadillo Standalone) ---
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


# --- 3. Project configuration and targets ---
# Include all discovered directories so the compiler can find the headers
include_directories(
        ${R_INCLUDE_DIRS}   # Core R headers
        ${RCPP_INC}         # Rcpp headers
        ${ARMA_INC}         # Armadillo headers
        ${BH_INC}           # Boost headers
        "."                 # Current directory (for BayesCore.h)
)

# Add the R library directory to the linker search path
link_directories(${R_LIBRARY_DIR})

# --- 3.1. TARGET 1: debug executable (for CLion testing) ---
# This target allows you to run and debug your C++ logic directly in CLion.
# It requires 'main.cpp' to exist.
# usage: Select 'BayesDebug' in CLion and press Run.
add_executable(BayesDebug
        BayesCore.h
        BayesCore.cpp
        main.cpp
)

# Link against the R library (libR).
# This is crucial for BLAS/LAPACK functions used by Armadillo.
target_link_libraries(BayesDebug
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
)

# --- 3.2. TARGET 2: shared library (for RcppInterface check) ---
# This target builds the library that R would load.
# In CLion, this is useful to check 'RcppInterface.cpp' for syntax errors.
# Note: You cannot "run" this target directly, but you can build it.
add_library(BayesRcppLib SHARED
        RcppInterface.cpp
        BayesCore.h
        BayesCore.cpp
        main.cpp
)

target_link_libraries(BayesRcppLib R)