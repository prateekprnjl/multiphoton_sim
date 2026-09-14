cmake_minimum_required(VERSION 3.20)

project(main LANGUAGES CXX)

option(ENABLE_MPI "Enable MPI parallel support" ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE PkgConfig::PETSC)

if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(main PRIVATE MPI::MPI_CXX)
    target_compile_definitions(main PRIVATE USE_MPI)
endif()
