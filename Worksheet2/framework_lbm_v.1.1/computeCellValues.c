#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "helper.h"

/*
 * Computes the density within the current cell (according to eq.(9)) 
 * and stores the result into the density.
 */
void computeDensity(const double *const currentCell, double *density){
	int i;

	*density=.0;
	for(i=0; i<Q; i++)
		*density = *density + currentCell[i];
}

/*
 * Computes the velocity within the current cell (according to eq.(9)) 
 * and stores the result into the velocity.
 */
void computeVelocity(const double * const currentCell, const double * const density, double *velocity){
	int i, d;

	for(d=0; d<DIM; d++)
		velocity[d] = .0;

	for(i=0; i<Q; i++)
		for(d=0; d<DIM; d++)
			velocity[d] += currentCell[i] * (double)LATTICEVELOCITIES[i][d];
			
	for(d=0; d<DIM; d++)
		velocity[d] /= (*density);
}

/*
 * Computes the equilibrium distribution (according to eq.(10)) from the density and the velocity
 * and stores the result into the feq.
 */
void computeFeq(const double * const density, const double * const velocity, double *feq){
	unsigned long i;

	for(i=0; i<Q; i++)
		feq[i] = LATTICEWEIGHTS[i] * (*density) * (1.0 + DOTP(LATTICEVELOCITIES[i], velocity) / SQ(C_S)
			+ SQ( DOTP(LATTICEVELOCITIES[i], velocity) ) / ( 2.0 * SQ(SQ(C_S)) )
			- DOTP(velocity, velocity) / ( 2.0 * SQ(C_S) ));
}



