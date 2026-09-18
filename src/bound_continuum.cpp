#include "bound_continuum.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>

#include <petscsys.h>

PetscErrorCode BoundContinuum::load(const char* filename){

    /* Open the bound-continuum dipoles */
    std::ifstream dipoles_b_c(filename);

    if (!dipoles_b_c){
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Bound-Continuum file missing.\n");
    }

    /* Temporary C++ vectors for reading from file */
    std::vector<double> energy_values;
    std::vector<std::vector<std::complex<double>>> dipole_values;
    
    /* Reading values from dipole file w/ header and variable number of channels */
    std::string line;
    std::size_t i=0;
    while(std::getline(dipoles_b_c, line)){

        std::stringstream ss(line);

        /* Header line "# number_of_lines"*/
        if (i==0){
            char hash_c;
            if (!(ss >> hash_c >> n_grid) || hash_c != '#'){
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, "Invalid Bound-Continuum Dipole header.\n");
            }
            ++i;
            continue;
        }

        /* Dipole lines format: 
        Energy number_of_channels real_part_1 complex_part_1 real_part_2 complex_part_2 ... complex_part_{num_of_channels} */
        double energy;
        std::size_t channels;

        if (!(ss >> energy >> channels)){
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, "Broken energy-channel in Bound-Continuum Dipole file @ Line: %d\n", static_cast<PetscInt>(i));
        }

        /* Store energy */
        energy_values.push_back(energy);

        /* Track maximum number of channels */
        if (channels > n_channels){
            n_channels = channels;
        }

        /* Reading complex values*/
        std::vector<std::complex<double>> values(channels);

        for (std::size_t j=0; j < channels; j++){
            double real_part, imaginary_part;

            if (!(ss >> real_part >> imaginary_part)){
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, "Broken channel value line in Bound-Continuum Dipole file @ Line: %d\n", static_cast<PetscInt>(i));
            }

            values[j] = std::complex<double> (real_part, imaginary_part);
        }
        dipole_values.push_back(values);
        ++i;
    }

    /* Pad with 0 + 0j for missing channels */
    for (auto& row : dipole_values){
        row.resize(n_channels, std::complex<double>(0.0, 0.0));
    }

    /* Convert to PETSc Complex Vectors */
    VecCreate(PETSC_COMM_WORLD, &energies_);
    VecSetSizes(energies_, PETSC_DECIDE, n_grid);
    VecSetFromOptions(energies_);

    PetscInt total = static_cast<PetscInt>(n_grid * n_channels);

    VecCreate(PETSC_COMM_WORLD, &dipoles_);
    VecSetSizes(dipoles_, PETSC_DECIDE, total);
    VecSetFromOptions(dipoles_);

    for (std::size_t i = 0; i < n_grid; ++i){
        
        VecSetValue(energies_, static_cast<PetscInt>(i), energy_values[i], INSERT_VALUES);
        
        for(std::size_t j = 0; j < n_channels; ++j){
            PetscInt index = static_cast<PetscInt>(i * n_channels + j);
            PetscScalar value = dipole_values[i][j];
            VecSetValue(dipoles_, index, value, INSERT_VALUES);
        }
    }
    VecAssemblyBegin(energies_);
    VecAssemblyEnd(energies_);
    VecAssemblyBegin(dipoles_);
    VecAssemblyEnd(dipoles_);

    return PETSC_SUCCESS;
}

std::size_t BoundContinuum::grid_size() const{
    return n_grid;
}

std::size_t BoundContinuum::number_of_channels() const{
    return n_channels;
}

Vec BoundContinuum::energies() const{
    return energies_;
}

Vec BoundContinuum::dipoles() const{
    return dipoles_;
}

BoundContinuum::~BoundContinuum(){
    if (energies_ != nullptr){
        VecDestroy(&energies_);
    }
    if (dipoles_ != nullptr){
        VecDestroy(&dipoles_);
    }
}
