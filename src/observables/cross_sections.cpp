#include "basis.hpp"
#include "bound_continuum.hpp"
#include "cross_sections.hpp"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>

#include <petscsys.h>
#include <petscmat.h>


PetscErrorCode CrossSections::calculate(const char* file_in, const BoundContinuum& bc){
    PetscFunctionBeginUser;

    const double MBn = 28.002852053;
    const double sol = 1.0 / 137.0;

    // Read initial-state energy on rank 0 and broadcast it.
    PetscMPIInt rank;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

    PetscInt read_error = 0;

    if (rank == 0) {
        std::ifstream energy_file(file_in);
        std::string line;
        std::size_t index;

        if (!energy_file) {
            read_error = 1;
        } 
        else if (!std::getline(energy_file, line) || !std::getline(energy_file, line)) {
            read_error = 2;
        } 
        else {
            std::stringstream ss(line);

            if (!(ss >> index >> E_0)) {
                read_error = 2;
            }
        }
    }

    PetscCallMPI(MPI_Bcast(&read_error, 1, MPIU_INT, 0, PETSC_COMM_WORLD));

    PetscCheck(read_error == 0, PETSC_COMM_WORLD, read_error == 1 ? PETSC_ERR_FILE_OPEN : PETSC_ERR_FILE_READ, "Could not read the initial-state energy.");

    PetscCallMPI(MPI_Bcast(&E_0, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD));

    n_grid_ = bc.grid_size();
    n_channels_ = bc.number_of_channels();

    Vec energies = bc.energies();
    Mat dipoles = bc.dipoles();

    // Make the energy values available on every rank.
    Vec energies_all;
    VecScatter scatter;

    PetscCall(VecScatterCreateToAll(energies, &scatter, &energies_all));
    PetscCall(VecScatterBegin(scatter, energies, energies_all, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(scatter, energies, energies_all, INSERT_VALUES, SCATTER_FORWARD));

    const PetscScalar* energy_array;
    PetscCall(VecGetArrayRead(energies_all, &energy_array));

    // Create the cross-section matrix.
    if (sigma_ != nullptr) {
        PetscCall(MatDestroy(&sigma_));
    }

    PetscCall(MatCreate(PETSC_COMM_WORLD, &sigma_));
    PetscCall(MatSetSizes(sigma_, PETSC_DECIDE, PETSC_DECIDE, static_cast<PetscInt>(n_grid_), static_cast<PetscInt>(n_channels_)));
    PetscCall(MatSetFromOptions(sigma_));
    PetscCall(MatSetUp(sigma_));

    // Process only locally owned dipole rows.
    PetscInt row_start, row_end;
    PetscCall(MatGetOwnershipRange(dipoles, &row_start, &row_end));

    std::vector<PetscInt> columns(n_channels_);
    std::vector<PetscScalar> dipole_row(n_channels_);

    for (std::size_t j = 0; j < n_channels_; ++j) {
        columns[j] = static_cast<PetscInt>(j);
    }

    for (PetscInt i = row_start; i < row_end; ++i) {
        PetscCall(MatGetValues(dipoles, 1, &i, static_cast<PetscInt>(n_channels_), columns.data(), dipole_row.data()));

        const double omega = PetscRealPart(energy_array[i]) - E_0;

        for (std::size_t j = 0; j < n_channels_; ++j) {

            const double cross_section = MBn * sol * 4.0 * M_PI * M_PI * PetscAbsScalar(dipole_row[j]) *
                PetscAbsScalar(dipole_row[j]) / (3.0 * omega);

            PetscCall(MatSetValue(sigma_, i, static_cast<PetscInt>(j), cross_section, INSERT_VALUES));
        }
    }

    PetscCall(MatAssemblyBegin(sigma_, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(sigma_, MAT_FINAL_ASSEMBLY));

    PetscCall(VecRestoreArrayRead(energies_all, &energy_array));
    PetscCall(VecScatterDestroy(&scatter));
    PetscCall(VecDestroy(&energies_all));

    PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode CrossSections::cation_cross_sections(const Basis& basis){
    PetscFunctionBeginUser;

    const std::size_t n_cations = basis.cs().size();

    cation_sigma_.resize(n_cations, nullptr);

    for (std::size_t c = 0; c < n_cations; ++c) {
        PetscCall(VecCreate(PETSC_COMM_WORLD, &cation_sigma_[c]));
        PetscCall(VecSetSizes(cation_sigma_[c], PETSC_DECIDE, static_cast<PetscInt>(n_grid_)));
        PetscCall(VecSetFromOptions(cation_sigma_[c]));
        PetscCall(VecSet(cation_sigma_[c], 0.0));
    }

    PetscInt row_start, row_end;
    PetscCall(MatGetOwnershipRange(sigma_, &row_start, &row_end));

    std::vector<PetscInt> columns(n_channels_);
    std::vector<PetscScalar> row(n_channels_);

    for (std::size_t j = 0; j < n_channels_; ++j) {
        columns[j] = static_cast<PetscInt>(j);
    }

    for (PetscInt i = row_start; i < row_end; ++i) {
        PetscCall(MatGetValues(sigma_, 1, &i, static_cast<PetscInt>(n_channels_), columns.data(), row.data()));

        for (std::size_t j = 0; j < n_channels_; ++j) {
            const std::size_t cation = basis.cation_for_channel(j);

            PetscCall(VecSetValue(cation_sigma_[cation], i, row[j], ADD_VALUES));
        }
    }

    for (Vec& vec : cation_sigma_) {
        PetscCall(VecAssemblyBegin(vec));
        PetscCall(VecAssemblyEnd(vec));
    }

    PetscFunctionReturn(PETSC_SUCCESS);
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
    PetscFunctionBeginUser;

    const double eV = 27.2114079527;

    PetscMPIInt rank;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

    // Gather energies to rank 0.
    VecScatter energy_scatter;
    Vec energy_root = nullptr;

    PetscCall(VecScatterCreateToZero(bc.energies(), &energy_scatter, &energy_root));

    PetscCall(VecScatterBegin(energy_scatter, bc.energies(), energy_root, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(energy_scatter, bc.energies(), energy_root, INSERT_VALUES, SCATTER_FORWARD));

    // Gather each cation's cross-section vector to rank 0.
    std::vector<VecScatter> cation_scatters(cation_sigma_.size());
    std::vector<Vec> cation_root(cation_sigma_.size(), nullptr);

    for (std::size_t c = 0; c < cation_sigma_.size(); ++c) {
        PetscCall(VecScatterCreateToZero(cation_sigma_[c], &cation_scatters[c], &cation_root[c]));

        PetscCall(VecScatterBegin(cation_scatters[c], cation_sigma_[c], cation_root[c], INSERT_VALUES, SCATTER_FORWARD));
        PetscCall(VecScatterEnd(cation_scatters[c], cation_sigma_[c], cation_root[c], INSERT_VALUES, SCATTER_FORWARD));
    }

    // Gather the total cross-section vector to rank 0.
    VecScatter total_scatter;
    Vec total_root = nullptr;

    PetscCall(VecScatterCreateToZero(total_sigma_, &total_scatter, &total_root));

    PetscCall(VecScatterBegin(total_scatter, total_sigma_, total_root, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(total_scatter, total_sigma_, total_root, INSERT_VALUES, SCATTER_FORWARD));

    // Only rank 0 writes the file.
    PetscInt write_error = 0;

    if (rank == 0) {
        std::ofstream file(file_out);

        if (!file) {
            write_error = 1;
        } else {
            // Header
            file << "# Energy(eV)";

            for (std::size_t c = 0; c < basis.cs().size(); ++c) {
                file << "\tCation" << basis.cs()[c];
            }

            file << "\tTotal\n";

            // Access the gathered arrays.
            const PetscScalar* energy_array;
            PetscCall(VecGetArrayRead(energy_root, &energy_array));

            std::vector<const PetscScalar*> cation_arrays(cation_root.size());

            for (std::size_t c = 0; c < cation_root.size(); ++c) {
                PetscCall(VecGetArrayRead(cation_root[c], &cation_arrays[c]));
            }

            const PetscScalar* total_array;
            PetscCall(VecGetArrayRead(total_root, &total_array));

            for (std::size_t i = 0; i < n_grid_; ++i) {
                const double energy_eV = (PetscRealPart(energy_array[i]) - E_0) * eV;

                file << energy_eV;

                for (std::size_t c = 0;
                     c < cation_arrays.size(); ++c) {
                    file << '\t' << PetscRealPart(cation_arrays[c][i]);
                }

                file << '\t' << PetscRealPart(total_array[i]) << '\n';
            }

            PetscCall(VecRestoreArrayRead(energy_root, &energy_array));

            for (std::size_t c = 0; c < cation_root.size(); ++c) {
                PetscCall(VecRestoreArrayRead(cation_root[c], &cation_arrays[c]));
            }

            PetscCall(VecRestoreArrayRead(total_root, &total_array));

            if (!file.good()) {
                write_error = 2;
            }
        }
    }

    // Communicate file-writing errors to every rank.
    PetscCallMPI(MPI_Bcast(&write_error, 1, MPIU_INT, 0, PETSC_COMM_WORLD));

    // Destroy the scatter contexts and gathered vectors.
    PetscCall(VecScatterDestroy(&energy_scatter));
    PetscCall(VecDestroy(&energy_root));

    for (std::size_t c = 0; c < cation_sigma_.size(); ++c) {
        PetscCall(VecScatterDestroy(&cation_scatters[c]));

        if (cation_root[c] != nullptr) {
            PetscCall(VecDestroy(&cation_root[c]));
        }
    }

    PetscCall(VecScatterDestroy(&total_scatter));
    PetscCall(VecDestroy(&total_root));

    PetscCheck(write_error == 0, PETSC_COMM_WORLD, write_error == 1 ? PETSC_ERR_FILE_OPEN : PETSC_ERR_FILE_WRITE, "Failed to write cross-section output file.");

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode CrossSections::calculate_all(const char* file_in, const Basis& basis, 
        const BoundContinuum& bc, const char* file_out){
    CrossSections::calculate(file_in, bc);
    CrossSections::cation_cross_sections(basis);
    CrossSections::total_cross_sections();
    CrossSections::write(file_out, basis, bc);
    
    return PETSC_SUCCESS;
}

CrossSections::~CrossSections()
{
    if (sigma_ != nullptr) {
        MatDestroy(&sigma_);
    }

    for (Vec& vec : cation_sigma_) {
        if (vec != nullptr) {
            VecDestroy(&vec);
        }
    }

    if (total_sigma_ != nullptr) {
        VecDestroy(&total_sigma_);
    }
}