// === TDSE Parallel Solver === //
// === PetSc + MPI implementation === //

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <petsc.h>

#define PI 3.1415926535
#define ROOT 0
#define MAX_ITER 1000

Vec psi;
Vec E;
Mat A;
double omega;
int rank;

typedef struct{
	double amplitude;
	double omega;
	double t_ini;
	double t_end;
} LaserPulse;

void initPETSc();
void initialize();
void loadInitialState(int vec);
void buildHamiltonian();
double externalPotential(double t);
void solveTDSE();
void writeOutput();

int main(){
	initPETSc();
	loadHamiltonian();
	loadInitialState(vec);
	buildHamiltonian();
	solveTDSE();
	writeOutput();
	cleanup();
}

