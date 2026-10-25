#include <basis.hpp>
#include <bound_continuum.hpp>

#include <catch2/catch_test_macros.hpp>
#include <petscsys.h>

TEST_CASE("Bound-continuum file loads correctly"){
    BoundContinuum bc;
    Basis bas;

    REQUIRE(bas.load("../input/scat1.info") == PETSC_SUCCESS);
    REQUIRE(bc.load("../input/DipoleTransAmpVel_z_FromBoxState_1_1", bas.number_of_channels()) == PETSC_SUCCESS);

    PetscInt len_energies;
    PetscInt len_dipoles_x, len_dipoles_y;

    REQUIRE(VecGetSize(bc.energies(), &len_energies) == PETSC_SUCCESS);
    REQUIRE(MatGetSize(bc.dipoles(), &len_dipoles_x, &len_dipoles_y) == PETSC_SUCCESS);

    CHECK(bc.grid_size() == static_cast<PetscInt>(len_energies));
    CHECK(bc.number_of_channels() * bc.grid_size() == static_cast<PetscInt>(len_dipoles_x * len_dipoles_y));
}
