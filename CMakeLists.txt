cmake_minimum_required(VERSION 3.20)

project(main LANGUAGES CXX)

option(ENABLE_MPI "Enable MPI parallel support" ON)

find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC-COMPLEX REQUIRED IMPORTED_TARGET PETSc)

add_executable(dipole src/dipole.cpp)

target_link_libraries(dipole PRIVATE PkgConfig::PETSC-COMPLEX)

if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(dipole PRIVATE MPI::MPI_CXX)
    target_compile_definitions(dipole PRIVATE USE_MPI)
endif()
