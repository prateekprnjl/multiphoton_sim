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
        PetscCall(sigma.calculate_all("../input/H_Eigenvalues", basis, bc, "../output/cross_sections.dat"));

        PetscPrintf(PETSC_COMM_WORLD, "Number of cations: %d\n", static_cast<PetscInt>(basis.cs().size()));
        PetscPrintf(PETSC_COMM_WORLD, "Maximum channels: %d\n", static_cast<PetscInt>(basis.number_of_channels()));
        PetscPrintf(PETSC_COMM_WORLD, "Grid size: %d\n", static_cast<PetscInt>(bc.grid_size()));
    }

    PetscFinalize();

    return 0;
}