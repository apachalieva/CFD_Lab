#include "initLB.h"
#include "LBDefinitions.h"
#include "helper.h"

/*
 * Read the parameters - domainsize xlength, relaxation time tau, velocityWall, 
 * number of time steps and interval between subsequent VTK outputs. 
 */
int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){

    double *velocityWall_0;
    double *velocityWall_1;
    double *velocityWall_2;
    char *szFileName;
   
    if( argc == 2 )
    {
      szFileName = argv[ 1 ];
    }
    else
    {
      printf( "WARRNING: No data-file given. Use the cavity.dat as default.\n" );
      szFileName = "cavity.dat";
    }
    read_int( szFileName, "xlength", xlength );
    read_int( szFileName, "timesteps", timesteps );
    read_int( szFileName, "timestepsPerPlotting", timestepsPerPlotting );
    
    read_double( szFileName, "tau", tau );
    velocityWall_0 = &velocityWall[ 0 ];
    velocityWall_1 = &velocityWall[ 1 ];
    velocityWall_2 = &velocityWall[ 2 ];
    read_double( szFileName, "velocityWall_0", velocityWall_0 );
    read_double( szFileName, "velocityWall_1", velocityWall_1 );
    read_double( szFileName, "velocityWall_2", velocityWall_2 );
    
    return 1;
}

/*
 * Initialise the fields with velocity = 0 and density = 1
 * The flagField stores information about the geometry: 
 * FLUID = 0, NO_SLIP = 1, MOVING_WALL = 2  
 */
void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
	int x, y, z, i, ix;

	for(z=0; z < xlength+2; z++)
		for(y=0; y < xlength+2; y++)
			for(x=0; x < xlength+2; x++){
				ix = x + (xlength+2) * y + SQ(xlength+2) * z;

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

