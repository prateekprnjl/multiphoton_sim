#include "basis.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

#include <petscsys.h>

PetscErrorCode Basis::load(const char* filename){

    /* Open scattering states label file */
    std::ifstream scat(filename);

    if (!scat){
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Scattering label file missing.\n");
    }
    
    /* Reading values from dipole file w/ header and variable number of channels */
    std::string line;
    std::size_t i=0;
    while(std::getline(scat, line)){

        std::stringstream ss(line);
        std::string label;

        /* Header line "number_of_channels energy"*/
        if (i==0){
            if (!(ss >> n_channels)){
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, "Invalid Scat info header.\n");
            }
            ++i;
            continue;
        }

        /* Reading channel labels from 2A.3X4.-1 where cation_state = 3, (l, m) = (4, -1) */
        ss >> label;
        
        std::size_t a = label.find('A');
        std::size_t x = label.find('X');
        std::size_t dot = label.find('.', x);

        std::size_t cation = std::stoi(label.substr(a + 2, x - a - 1));

        int l = std::stoi(label.substr(x + 1, dot - x - 1));
        int m = std::stoi(label.substr(dot + 1));

        /* Check if cation state exists, else create a new row */
        if (cs_.empty() || cs_.back() != cation) {
            cs_.push_back(cation);
            l_m_.push_back({});
        }
        l_m_.back().push_back({l, m}); 

        ++i;
    }

    return PETSC_SUCCESS;
}

std::size_t Basis::number_of_channels() const{
    return n_channels;
}

const std::vector<std::size_t>& Basis::cs() const{
    return cs_;
}

const std::vector<std::vector<std::pair<int, int>>>& Basis::l_m() const{
    return l_m_;
}
