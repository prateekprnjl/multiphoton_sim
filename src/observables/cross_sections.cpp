#include "basis.hpp"
#include "bound_continuum.hpp"
#include "cross_sections.hpp"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <cmath>

#include <petscsys.h>

PetscErrorCode CrossSections::calculate(const char* file_in, const BoundContinuum& bc){

    const double MBn = 28.002852053;    // squared Bohr radius in SI units
    const double sol = 1.0/137.0;           // Speed of light in Hartree units

    /* Read the initial state energy E_0 */
    std::ifstream energy_file(file_in);
    if (!energy_file){
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "H_Eigenvalues file missing.\n");
    }
    std::size_t i=0;
    std::string line;
    std::getline(energy_file, line); std::getline(energy_file, line);
    std::stringstream ss(line);
    if (!(ss >> i >> E_0)){
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, "Invalid format of H_Eigenvalues file.\n");
    }

    /* Collect data for cross-section calculation */
    n_grid_ = bc.grid_size();
    n_channels_ = bc.number_of_channels();
    Vec energies = bc.energies();
    Vec dipoles = bc.dipoles();

    /* Create vector to store cross-sections */
    VecCreate(PETSC_COMM_WORLD, &sigma_);
    VecSetSizes(sigma_, PETSC_DECIDE, n_grid_ * n_channels_);
    VecSetFromOptions(sigma_);

    const PetscScalar* energy_array;
    const PetscScalar* dipole_array;
    PetscScalar* sigma_array;

    VecGetArrayRead(energies, &energy_array);
    VecGetArrayRead(dipoles, &dipole_array);
    VecGetArray(sigma_, &sigma_array);

    /* Calculate Cross-sections */
    for (std::size_t i = 0; i < n_grid_; ++i) {

        const double omega = PetscRealPart(energy_array[i]) - E_0;

        for (std::size_t j = 0; j < n_channels_; ++j) {

            const std::size_t index = i * n_channels_ + j;
            const std::complex<double> dipole = dipole_array[index];

            const double sigma = MBn * sol * 4.0 * M_PI * M_PI * std::norm(dipole) / (3.0 * omega);
            sigma_array[index] = sigma;
        }
    }

    VecRestoreArrayRead(energies, &energy_array);
    VecRestoreArrayRead(dipoles, &dipole_array);
    VecRestoreArray(sigma_, &sigma_array);

    return PETSC_SUCCESS;
}

PetscErrorCode CrossSections::cation_cross_sections(const Basis& basis){
    const std::size_t n_cations = basis.cs().size();

    /* Total number of points */
    PetscInt n_total;
    PetscCall(VecGetSize(sigma_, &n_total));

    /* Make a n-cation column vector for individual cation cross-sections */
    cation_sigma_.resize(n_cations, nullptr);
    for (std::size_t c = 0; c < n_cations; ++c) {
        PetscCall(VecCreate(PETSC_COMM_WORLD, &cation_sigma_[c]));
        PetscCall(VecSetSizes(cation_sigma_[c], PETSC_DECIDE, n_grid_));
        PetscCall(VecSetFromOptions(cation_sigma_[c]));
        PetscCall(VecSet(cation_sigma_[c], 0.0));
    }

    /* Add in values for each cation state for different energy points */
    const PetscScalar* sigma_array;
    PetscCall(VecGetArrayRead(sigma_, &sigma_array));
    for (std::size_t i = 0; i < n_grid_; ++i) {
        for (std::size_t j = 0; j < n_channels_; ++j) {
            const std::size_t index = i * n_channels_ + j;

            /* Use the helper function to find the cation state for each channel number */
            const std::size_t cation = basis.cation_for_channel(j);

            /* ADD to the existing value at that energy point */
            PetscCall(VecSetValue(cation_sigma_[cation], static_cast<PetscInt>(i), sigma_array[index],
                ADD_VALUES));
        }
    }

    PetscCall(VecRestoreArrayRead(sigma_, &sigma_array));

    for (Vec& vec : cation_sigma_) {
        PetscCall(VecAssemblyBegin(vec));
        PetscCall(VecAssemblyEnd(vec));
    }
    return PETSC_SUCCESS;
}

PetscErrorCode CrossSections::total_cross_sections(){
    if (total_sigma_ != nullptr) {
        PetscCall(VecDestroy(&total_sigma_));
    }

    PetscCall(VecCreate(PETSC_COMM_WORLD, &total_sigma_));
    PetscCall(VecSetSizes(total_sigma_, PETSC_DECIDE, n_grid_));
    PetscCall(VecSetFromOptions(total_sigma_));
    PetscCall(VecSet(total_sigma_, 0.0));

    for (const Vec& cation : cation_sigma_) {
        PetscCall(VecAXPY(total_sigma_, 1.0, cation));
    }

    return PETSC_SUCCESS;
}

PetscErrorCode CrossSections::write(const char* file_out, const Basis& basis, const BoundContinuum& bc) const{
    const double eV = 27.2114079527;    // Hartree to eV conversion
    std::ofstream file(file_out);

    if (!file) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Could not open cross-section output file.\n");
    }

    /* Header */
    file << "# Energy(eV)";

    for (std::size_t c = 0; c < basis.cs().size(); ++c) {
        file << "\tCation" << basis.cs()[c];
    }

    file << "\tTotal\n";

    /* Access energies */
    const PetscScalar* energy_array;
    PetscCall(VecGetArrayRead(bc.energies(), &energy_array));

    /* Access cation cross sections */
    std::vector<const PetscScalar*> cation_arrays(cation_sigma_.size());
    for (std::size_t c = 0; c < cation_sigma_.size(); ++c) {
        PetscCall(VecGetArrayRead(cation_sigma_[c], &cation_arrays[c]));
    }
    const PetscScalar* total_array;
    PetscCall(VecGetArrayRead(total_sigma_, &total_array));

    for (std::size_t i = 0; i < n_grid_; ++i) {
        const double energy_eV = (PetscRealPart(energy_array[i]) - E_0) * eV;
        file << energy_eV;
        for (std::size_t c = 0; c < cation_sigma_.size(); ++c) {
            file << '\t' << PetscRealPart(cation_arrays[c][i]);
        }
        file << '\t' << PetscRealPart(total_array[i])<< '\n';
    }

    /* Restore arrays */
    PetscCall(VecRestoreArrayRead(bc.energies(), &energy_array));
    for (std::size_t c = 0; c < cation_sigma_.size(); ++c) {
        PetscCall(VecRestoreArrayRead(cation_sigma_[c], &cation_arrays[c]));
    }
    PetscCall(VecRestoreArrayRead(total_sigma_, &total_array));

    return PETSC_SUCCESS;
}

PetscErrorCode CrossSections::calculate_all(const char* file_in, const Basis& basis, 
        const BoundContinuum& bc, const char* file_out){
    CrossSections::calculate(file_in, bc);
    CrossSections::cation_cross_sections(basis);
    CrossSections::total_cross_sections();
    CrossSections::write(file_out, basis, bc);
    
    return PETSC_SUCCESS;
}

CrossSections::~CrossSections(){
    if (sigma_ != nullptr){
        VecDestroy(&sigma_);
    }
}
