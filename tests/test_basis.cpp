#include <basis.hpp>
#include <bound_continuum.hpp>

#include <catch2/catch_test_macros.hpp>
#include <petscsys.h>

TEST_CASE("Scat-info file loads correctly"){
    Basis bas;
    BoundContinuum bc;

    REQUIRE(bas.load("../input/scat1.info") == PETSC_SUCCESS);
    REQUIRE(bc.load("../input/DipoleTransAmpVel_z_FromBoxState_1_1") == PETSC_SUCCESS);

    CHECK(bc.number_of_channels() == bas.number_of_channels());
    CHECK(bas.cs().size() == bas.l_m().size());

    /* Total number of (l, m) entries = the number of channels */
    std::size_t n_entries = 0;
    for (const auto& row : bas.l_m()){
        n_entries += row.size(); 
    }

    CHECK(bas.number_of_channels() == n_entries);
}