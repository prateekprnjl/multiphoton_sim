/* TDSE Parallel Solver 
   PetSc + MPI implementation 
*/

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <petsc.h>
#include <petscts.h>
#include <petscvec.h>
#include <petscmat.h>

#define ROOT 0
#define PI 3.1415926535

Vec E;			/* Energy vector */
Mat A; 			/* Interaction matrix */
Vec psi0;		/* Initial wavefunction	*/
Vec psi;		/* Working wavefunction	*/
int rank, size;		/* MPI */

/* Laser parameters */
typedef struct{
	double omega;
	double amplitude;

	double t_ini;
	double t_end;
} LaserPulse;

/* Pump-Probe field setup */
typedef struct{
	LaserPulse pump;
	LaserPulse probe;

	double delay;
	double t_final;
} SimulationParameters;

void initPETSc(int argc,char **argv);
void loadMatrix(void);
void loadInitialState(void);

void solveTDSE(MPI_Comm subcomm, SimulationParameters job);
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec u,Vec f,void *ctx);
double externalPotential(double t, SimulationParameters *job);

void writeOutput(MPI_Comm subcomm,int delayIndex);
void cleanup(void);

int main(){
	int cpusPerJob;
	initPETSc(argc,argv);
	parseInputs(argc,argv,&cpusPerJob);
	loadMatrix();
	loadInitialState();

	cleanup();
	PetscFinalize();

	return 0;
}

