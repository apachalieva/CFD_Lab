#include "collision.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include "helper.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
	int i;

	for(i=0ul; i<Q; i++)
		currentCell[i] = currentCell[i] + (feq[i] - currentCell[i]) / (*tau);
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){
	double density, velocity[3], feq[Q], *currentCell;
	int x, y, z;

	for(z=1; z<xlength+1; z++)
		for(y=1; y<xlength+1; y++)
			for(x=1; x<xlength+1; x++){
				/* get the current cell */
				currentCell = collideField + Q * ( x + (xlength+2) * y + SQ(xlength+2) * z );
			
				/* compute collision */
				computeDensity( currentCell, &density );
				computeVelocity( currentCell, &density, velocity );
				computeFeq( &density, velocity,feq );
				computePostCollisionDistributions( currentCell, tau, feq );
			    }
}
