cmake_minimum_required(VERSION 3.20)

project(main LANGUAGES CXX)

enable_testing()
option(ENABLE_MPI "Enable MPI parallel support" ON)

find_package(PkgConfig REQUIRED)
find_package(Catch2 3 REQUIRED)

pkg_check_modules(PETSC-COMPLEX REQUIRED IMPORTED_TARGET PETSc)

add_executable(main 
    main.cpp 
    src/basis.cpp
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
add_executable(test_multi_ph
    tests/test_main.cpp 
    tests/test_basis.cpp
    tests/test_bound_continuum.cpp 
    src/basis.cpp
    src/bound_continuum.cpp
    )

target_include_directories(test_multi_ph PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(test_multi_ph PRIVATE PkgConfig::PETSC-COMPLEX)
target_link_libraries(test_multi_ph PRIVATE Catch2::Catch2)

if(ENABLE_MPI)
    target_link_libraries(test_multi_ph PRIVATE MPI::MPI_CXX)
    target_compile_definitions(test_multi_ph PRIVATE USE_MPI)
endif()

add_test(
    NAME test_multi_ph 
    COMMAND test_multi_ph)
