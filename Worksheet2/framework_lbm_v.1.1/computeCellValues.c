#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "helper.h"

void computeDensity(const double *const currentCell, double *density){
	int i;

	*density=.0;
	for(i=0; i<Q; i++)
		*density = *density + currentCell[i];
}

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

void computeFeq(const double * const density, const double * const velocity, double *feq){
	unsigned long i;

	for(i=0; i<Q; i++)
		feq[i] = LATTICEWEIGHTS[i] * (*density) * (1.0 + DOTP(LATTICEVELOCITIES[i], velocity) / SQ(C_S)
			+ SQ( DOTP(LATTICEVELOCITIES[i], velocity) ) / ( 2.0 * SQ(SQ(C_S)) )
			- DOTP(velocity, velocity) / ( 2.0 * SQ(C_S) ));
}



