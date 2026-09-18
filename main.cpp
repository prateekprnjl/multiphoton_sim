#include "bound_continuum.hpp"

#include <petscsys.h>

int main(int argc, char **argv)
{
    PetscInitialize(&argc, &argv, nullptr, "TDSE calculation.");
    {
        BoundContinuum bc;

        PetscCall(bc.load("../input/DipoleTransAmpVel_z_FromBoxState_1_1"));

        PetscPrintf(PETSC_COMM_WORLD, "Grid size: %d\n", static_cast<PetscInt>(bc.grid_size()));
        PetscPrintf(PETSC_COMM_WORLD, "Maximum channels: %d\n", static_cast<PetscInt>(bc.number_of_channels()));
    }

    PetscFinalize();

    return 0;
}