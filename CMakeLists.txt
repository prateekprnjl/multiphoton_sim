cmake_minimum_required(VERSION 3.20)

project(main LANGUAGES CXX)

enable_testing()
option(ENABLE_MPI "Enable MPI parallel support" ON)

find_package(PkgConfig REQUIRED)
find_package(Catch2 3 REQUIRED)

pkg_check_modules(PETSC-COMPLEX REQUIRED IMPORTED_TARGET PETSc)

add_executable(main 
    main.cpp 
    src/bound_continuum.cpp
    )

target_include_directories(main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(main PRIVATE PkgConfig::PETSC-COMPLEX)

if(ENABLE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(main PRIVATE MPI::MPI_CXX)
    target_compile_definitions(main PRIVATE USE_MPI)
endif()

# --- TESTING --- #
add_executable(test_bound_continuum 
    tests/test_main.cpp 
    tests/test_bound_continuum.cpp 
    src/bound_continuum.cpp
    )

target_include_directories(test_bound_continuum PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(test_bound_continuum PRIVATE PkgConfig::PETSC-COMPLEX)
target_link_libraries(test_bound_continuum PRIVATE Catch2::Catch2)

if(ENABLE_MPI)
    target_link_libraries(test_bound_continuum PRIVATE MPI::MPI_CXX)
    target_compile_definitions(test_bound_continuum PRIVATE USE_MPI)
endif()

add_test(
    NAME test_bound_continuum 
    COMMAND test_bound_continuum)
