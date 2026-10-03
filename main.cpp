#include "basis.hpp"
#include "bound_continuum.hpp"

#include <petscsys.h>

int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, nullptr, "TDSE calculation.");
    {
        Basis bas;
        BoundContinuum bc;

        PetscCall(bas.load("../input/scat1.info"));
        PetscCall(bc.load("../input/DipoleTransAmpVel_z_FromBoxState_1_1"));

        PetscPrintf(PETSC_COMM_WORLD, "Grid size: %d\n", static_cast<PetscInt>(bc.grid_size()));
        PetscPrintf(PETSC_COMM_WORLD, "Number of cations: %d\n", static_cast<PetscInt>(bas.cs().size()));
        PetscPrintf(PETSC_COMM_WORLD, "Maximum channels: %d\n", static_cast<PetscInt>(bas.number_of_channels()));
    }

    PetscFinalize();

    return 0;
}