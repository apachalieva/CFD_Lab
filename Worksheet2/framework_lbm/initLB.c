#include "initLB.h"
#include "LBDefinitions.h"
#include "helper.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
    char *szFileName="cavity100.dat";
    if(argc==2) szFileName=argv[1];
    else ERROR("No parameter file specified");

    /* reading from the specified file */
    read_int(szFileName,"xlength",xlength);
    read_double(szFileName,"tau",tau);
    read_double(szFileName,"velocityWall",velocityWall);
    read_double(szFileName,"velocityWall",velocityWall+1);
    read_double(szFileName,"velocityWall",velocityWall+2);
    read_int(szFileName,"timesteps",timesteps);
    read_int(szFileName,"timestepsPerPlotting",timestepsPerPlotting);

  return 1;
}




void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
	unsigned long x, y, z, i, ix;

	for(z=0ul; z < xlength+2; z++)
		for(y=0ul; y < xlength+2; y++)
			for(x=0ul; x < xlength+2; x++){
				ix = x + xlength * y + SQ(xlength) * z;

				/* pdf intial condition */
				for(i=0ul; i < Q; i++){
					collideField[Q*ix+i] = LATTICEWEIGHTS[i];
					streamField[Q*ix+i] = LATTICEWEIGHTS[i];
				}

				/* set boundaries properties */
				if( y == xlength+1 )
					flagField[ix] = MOVING_WALL;
				else if( x==0 || x==xlength+1 || y==0 || z==0 || z==xlength+1 )
					flagField[ix] = NO_SLIP;
				else
					flagField[ix] = FLUID;
			}
}

