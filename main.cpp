#include <iostream>
#include <petscsys.h>

int main(int argc, char **argv){
    PetscCall(PetscInitialize(&argc, &argv, nullptr, nullptr));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "Hello, world!\n"));
    PetscCall(PetscFinalize());

    return 0;
}
