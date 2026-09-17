#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>
#include <petscsys.h>
#include <petscvec.h>

std::string line;
std::size_t n_grid;
std::size_t n_channels;
Vec bound_continuum;
Vec energies;

int main(){
    PetscInitialize(nullptr, nullptr, nullptr, nullptr);

    /* Open the bound-continuum dipoles */
    std::ifstream dipoles_b_c("../input/DipoleTransAmpVel_z_FromBoxState_1_1");

    if (!dipoles_b_c){
        std::cerr << "Bound-Continuum file missing.\n";

        PetscFinalize();
        return 1;
    }

    /* Temporary C++ vectors for reading from file */
    std::vector<double> energy_values;
    std::vector<std::vector<std::complex<double>>> dipole_values;
    
    /* Reading values from dipole file w/ header and variable number of channels */
    std::size_t i=0;
    while(std::getline(dipoles_b_c, line)){

        std::stringstream ss(line);

        /* Header line "# number_of_lines"*/
        if (i==0){
            char hash_c;
            if (!(ss >> hash_c >> n_grid || hash_c != '#')){
                std::cerr << "Invalid Bound-Continuum Dipole header.\n";
                
                PetscFinalize();
                return 1;
            }
            ++i;
            continue;
        }

        /* Dipole lines format: 
        Energy number_of_channels real_part_1 complex_part_1 real_part_2 complex_part_2 ... complex_part_{num_of_channels} */

        double energy;
        std::size_t channels;

        if (!(ss >> energy >> channels)){
            std::cerr << "Broken energy-channel in Bound-Continuum Dipole file @ Line" << i << ".\n";

            PetscFinalize();
            return 1;
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
                std::cerr << "Broken channel value line in Bound-Continuum Dipole file @ Line" << i << ".\n";

                PetscFinalize();
                return 1;
            }

            values[j] = std::complex<double> (real_part, imaginary_part);
        }
        dipole_values.push_back(values);
        ++i;
    }

    /* Check if number_of_energies = number_of_lines_in_file */
    if (energy_values.size() != n_grid){
        std::cerr << "Header mentions:" << n_grid << ", actual number of lines in file:" << energy_values.size() << ".\n";

        PetscFinalize();
        return 1;
    }

    /* Pad with 0 + 0j for missing channels */
    for (auto& row : dipole_values){
        row.resize(n_channels, std::complex<double>(0.0, 0.0));
    }

    /* Convert to PETSc Complex Vectors */
    VecCreate(PETSC_COMM_WORLD, &energies);
    VecSetSizes(energies, PETSC_DECIDE, n_grid);
    VecSetFromOptions(energies);

    PetscInt total = static_cast<PetscInt>(n_grid * n_channels);

    VecCreate(PETSC_COMM_WORLD, &bound_continuum);
    VecSetSizes(bound_continuum, PETSC_DECIDE, total);
    VecSetFromOptions(bound_continuum);

    for (std::size_t i = 0; i < n_grid; ++i){
        
        VecSetValue(energies, static_cast<PetscInt>(i), energy_values[i], INSERT_VALUES);
        
        for(std::size_t j = 0; j < n_channels; ++j){
            PetscInt index = static_cast<PetscInt>(i * n_channels + j);
            PetscScalar value = dipole_values[i][j];
            VecSetValue(bound_continuum, index, value, INSERT_VALUES);
        }
    }
    VecAssemblyBegin(energies);
    VecAssemblyEnd(energies);
    VecAssemblyBegin(bound_continuum);
    VecAssemblyEnd(bound_continuum);

    /* Sanity checks */
    PetscInt energy_size, dipole_size;

    VecGetSize(energies, &energy_size);
    VecGetSize(bound_continuum, &dipole_size);

    PetscPrintf(PETSC_COMM_WORLD, "Number of energies read: %d.\n", energy_size);
    PetscPrintf(PETSC_COMM_WORLD, "Number of maximum channels: %d.\n", static_cast<PetscInt>(n_channels));
    PetscPrintf(PETSC_COMM_WORLD, "Number of dipoles read (n_can x n_grid): %d.\n", dipole_size);

    /* Cleanup */
    VecDestroy(&energies);
    VecDestroy(&bound_continuum);

    PetscFinalize();

    return 0;
}
