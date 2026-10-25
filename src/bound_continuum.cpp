#include "bound_continuum.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <complex>
#include <vector>

#include <petscsys.h>
#include <petscvec.h>
#include <petscsys.h>

PetscErrorCode BoundContinuum::load(const char* filename, std::size_t number_of_channels)
{
    PetscFunctionBeginUser;

    MPI_Comm comm = PETSC_COMM_WORLD;
    PetscMPIInt rank;
    PetscCallMPI(MPI_Comm_rank(comm, &rank));

    const PetscInt nch = static_cast<PetscInt>(number_of_channels);

    /* Rank 0 reads the file. All ranks receive the parsed data. */
    PetscInt parse_error = 0;
    PetscInt grid_size = 0;

    std::vector<double> energy_values;
    std::vector<double> dipole_parts;

    if (rank == 0) {
        std::ifstream file(filename);

        if (!file) {
            parse_error = 1;
        } 
        else {
            std::string line;
            std::size_t line_number = 0;
            std::size_t declared_grid_size = 0;
            std::size_t last_channels = 0;

            bool header_read = false;

            while (std::getline(file, line)) {
                ++line_number;

                // Header: # number_of_lines
                if (!header_read) {
                    std::stringstream ss(line);
                    char hash;
                    std::size_t count;

                    if (!(ss >> hash >> count) || hash != '#' || count == 0 || count > static_cast<std::size_t>(PETSC_MAX_INT)) {
                        parse_error = 2;
                        break;
                    }

                    declared_grid_size = count;
                    header_read = true;

                    energy_values.reserve(count);
                    dipole_parts.reserve(2 * count * number_of_channels);

                    continue;
                }

                // Ignore blank lines.
                if (line.find_first_not_of(" \t\r") == std::string::npos) {
                    continue;
                }

                std::stringstream ss(line);
                double energy;
                std::size_t channels;

                if (!(ss >> energy >> channels)) {
                    parse_error = 3;
                    break;
                }

                std::vector<double> row(2 * number_of_channels, 0.0);

                for (std::size_t j = 0; j < channels; ++j) {
                    double re, im;
                    if (!(ss >> re >> im)) {
                        parse_error = 3;
                        break;
                    }
                    row[2 * j]     = re;
                    row[2 * j + 1] = im;
                }

                if (parse_error != 0) {
                    break;
                }

                energy_values.push_back(energy);
                dipole_parts.insert(dipole_parts.end(), row.begin(), row.end());
                last_channels = channels;
            }

            // Validate the complete file.
            if (parse_error == 0) {
                if (!header_read || energy_values.size() != declared_grid_size) {
                    parse_error = 5;
                } 
                else if (last_channels != number_of_channels) {
                    parse_error = 6;
                }
            }

            if (parse_error == 0) {
                grid_size = static_cast<PetscInt>(declared_grid_size);
            }
        }
    }

    {
    /* Check if file loaded correctly */
    PetscCallMPI(MPI_Bcast(&parse_error, 1, MPIU_INT, 0, comm));
    PetscCheck(parse_error == 0, comm, PETSC_ERR_FILE_READ,
               "Failed to read the Bound-Continuum dipole file. "
               "Error code: %" PetscInt_FMT
               ". Check the file, its header, row data, "
               "and final channel count.",
               parse_error);
    }

    PetscCallMPI(MPI_Bcast(&grid_size, 1, MPIU_INT, 0, comm));

    // Broadcast energies and dipole components.
    const PetscMPIInt energy_count = static_cast<PetscMPIInt>(grid_size);
    const PetscMPIInt component_count = static_cast<PetscMPIInt>(2 * static_cast<std::size_t>(grid_size) * number_of_channels);

    if (rank != 0) {
        energy_values.resize(static_cast<std::size_t>(grid_size));
        dipole_parts.resize(static_cast<std::size_t>(component_count));
    }

    PetscCallMPI(MPI_Bcast(energy_values.data(), energy_count, MPI_DOUBLE, 0, comm));
    PetscCallMPI(MPI_Bcast(dipole_parts.data(), component_count, MPI_DOUBLE, 0, comm));

    // Store metadata.
    n_grid = static_cast<std::size_t>(grid_size);
    max_channels_in_dipoles = number_of_channels;

    // Create the distributed energy vector.
    PetscCall(VecCreate(comm, &energies_));
    PetscCall(VecSetSizes(energies_, PETSC_DECIDE, grid_size));
    PetscCall(VecSetFromOptions(energies_));

    PetscInt energy_start, energy_end;
    PetscCall(VecGetOwnershipRange(energies_, &energy_start, &energy_end));

    for (PetscInt i = energy_start; i < energy_end; ++i) {
        PetscCall(VecSetValue(energies_, i, static_cast<PetscScalar>(energy_values[static_cast<std::size_t>(i)]), INSERT_VALUES));
    }

    PetscCall(VecAssemblyBegin(energies_));
    PetscCall(VecAssemblyEnd(energies_));
    PetscCall(MatCreate(comm, &dipoles_));
    PetscCall(MatSetSizes(dipoles_, PETSC_DECIDE, PETSC_DECIDE, grid_size, nch));
    PetscCall(MatSetFromOptions(dipoles_));
    PetscCall(MatSetUp(dipoles_));

    PetscInt row_start, row_end;
    PetscCall(MatGetOwnershipRange(dipoles_, &row_start, &row_end));

    // Each rank inserts only its locally owned rows.
    for (PetscInt i = row_start; i < row_end; ++i) {
        for (PetscInt j = 0; j < nch; ++j) {
            const std::size_t k = 2 * (static_cast<std::size_t>(i) * number_of_channels + static_cast<std::size_t>(j));

            const PetscScalar value = PetscCMPLX(dipole_parts[k], dipole_parts[k + 1]);

            PetscCall(MatSetValue(dipoles_, i, j, value, INSERT_VALUES));
        }
    }

    PetscCall(MatAssemblyBegin(dipoles_, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(dipoles_, MAT_FINAL_ASSEMBLY));

    PetscFunctionReturn(PETSC_SUCCESS);
}

std::size_t BoundContinuum::grid_size() const{
    return n_grid;
}

std::size_t BoundContinuum::number_of_channels() const{
    return max_channels_in_dipoles;
}

Vec BoundContinuum::energies() const{
    return energies_;
}

Mat BoundContinuum::dipoles() const{
    return dipoles_;
}

BoundContinuum::~BoundContinuum(){
    if (energies_ != nullptr){
        VecDestroy(&energies_);
    }
    if (dipoles_ != nullptr){
        MatDestroy(&dipoles_);
    }
}
