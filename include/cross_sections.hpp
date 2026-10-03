#ifndef CROSS_SECTIONS_HPP
#define CROSS_SECTIONS_HPP

#include <cstddef>
#include <petscvec.h>

class Basis;
class BoundContinuum;

class CrossSections {

public:

    PetscErrorCode calculate(const BoundContinuum& bc, const char* filename);
    PetscErrorCode cation_cross_sections(const Basis& basis);
    PetscErrorCode total_cross_sections();
    PetscErrorCode write(const char* filename, const Basis& basis, const BoundContinuum& bc) const;

    Vec total() const;

    ~CrossSections();

private:

    double E_0 = 0.0;
    
    std::size_t n_grid_ = 0;
    std::size_t n_channels_ = 0;

    Vec sigma_ = nullptr;
    Vec total_sigma_ = nullptr;

    std::vector<Vec> cation_sigma_;
};

#endif