#ifndef BOUND_CONTINUUM_HPP
#define BOUND_CONTINUUM_HPP

#include <cstddef>
#include <petscmat.h>
#include <petscvec.h>

class BoundContinuum {

public:

    PetscErrorCode load(const char* filename, std::size_t number_of_channels);

    std::size_t grid_size() const;
    std::size_t number_of_channels() const;

    Vec energies() const;
    Mat dipoles() const;

    ~BoundContinuum();

private:

    PetscInt n_grid = 0;
    PetscInt max_channels_in_dipoles = 0;

    Vec energies_ = nullptr;
    Mat dipoles_ = nullptr;
};

#endif