#include "basis.hpp"
#include "bound_continuum.hpp"

#include "cross_sections.hpp"

#include <petscsys.h>

int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, nullptr, "TDSE calculation.");
    {
        Basis basis;
        BoundContinuum bc;
        CrossSections sigma;

        PetscCall(basis.load("../input/scat1.info"));
        PetscCall(bc.load("../input/DipoleTransAmpVel_z_FromBoxState_1_1"));

        /* Observables - Cross sections (optional)*/
        PetscCall(sigma.calculate(bc, "../input/H_Eigenvalues"));
        PetscCall(sigma.cation_cross_sections(basis));
        PetscCall(sigma.total_cross_sections());
        PetscCall(sigma.write("../output/cross_sections.dat", basis, bc));

        PetscPrintf(PETSC_COMM_WORLD, "Grid size: %d\n", static_cast<PetscInt>(bc.grid_size()));
        PetscPrintf(PETSC_COMM_WORLD, "Number of cations: %d\n", static_cast<PetscInt>(basis.cs().size()));
        PetscPrintf(PETSC_COMM_WORLD, "Maximum channels: %d\n", static_cast<PetscInt>(basis.number_of_channels()));
    }

    PetscFinalize();

    return 0;
}