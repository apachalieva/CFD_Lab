#include "collision.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include "helper.h"

void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
	unsigned long i;

	for(i=0ul; i<Q; i++)
		currentCell[i] = currentCell[i] + (feq[i] - currentCell[i]) / *tau;
}

void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){
	double density, velocity[3], feq[Q], *currentCell;
	unsigned long x, y, z, i;

	for(z=1ul; z<xlength+1; z++){
		for(y=1ul; y<xlength+1; y++){
			for(x=1ul; x<xlength+1; x++){
				/* get the current cell */
			    for( i = 0; i < Q; ++i ){ /*you sure about this?*/
				currentCell = &collideField[ Q * ( x + (xlength+2) * y + SQ(xlength+2) * z ) + i ];
			
				/* compute collision */
				computeDensity( currentCell, &density );
				computeVelocity( currentCell, &density, velocity );
				computeFeq( &density, velocity,feq );
				computePostCollisionDistributions( currentCell, tau, feq );
			    }
				
			}
		}
	}
}

