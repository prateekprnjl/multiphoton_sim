#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include <cstddef>
#include <petscvec.h>

struct Hamiltonian {
    std::vector<std::vector<PetscComplex>> H;
};

#endif