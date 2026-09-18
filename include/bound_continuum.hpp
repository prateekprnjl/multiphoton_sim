#ifndef BOUND_CONTINUUM_HPP
#define BOUND_CONTINUUM_HPP

#include <cstddef>
#include <petscvec.h>

class BoundContinuum {

public:

    PetscErrorCode load(const char* filename);

    std::size_t grid_size() const;
    std::size_t number_of_channels() const;

    Vec energies() const;
    Vec dipoles() const;

    ~BoundContinuum();

private:

    std::size_t n_grid = 0;
    std::size_t n_channels = 0;

    Vec energies_ = nullptr;
    Vec dipoles_ = nullptr;
};

#endif