#ifndef BASIS_HPP
#define BASIS_HPP

#include <cstddef>
#include <vector>
#include <petscvec.h>

class Basis{

public:

    PetscErrorCode load(const char* filename);

    std::size_t number_of_channels() const;
    std::size_t cation_for_channel(std::size_t channel) const;
    
    const std::vector<std::size_t>& cs() const;
    const std::vector<std::vector<std::pair<int, int>>>& l_m() const;

    ~Basis();

private:

    std::size_t n_channels = 0;
    std::vector<std::size_t> cs_;
    std::vector<std::vector<std::pair<int, int>>> l_m_;
};

#endif