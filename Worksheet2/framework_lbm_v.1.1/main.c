#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"
#include <time.h>


int main (int argc, char *argv[]){
	double *collideField=NULL, *streamField=NULL, tau, velocityWall[3];
	int t, xlength, timesteps, timestepsPerPlotting, *flagField=NULL;
	clock_t start = clock();

	readParameters( &xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv);
	
	/* 
	 * Allocate memory blocks: collideField, streamField and flagField
	 */
	collideField = malloc(Q * CUBE(xlength+2) * sizeof(double));
	streamField  = malloc(Q * CUBE(xlength+2) * sizeof(double));
	flagField    = malloc(    CUBE(xlength+2) * sizeof(double));
	
	/*
	 * Initialise the fields with velocity = 0 and density = 1
	 * The flagField stores information about the geometry: 
	 * FLUID = 0, NO_SLIP = 1, MOVING_WALL = 2  
	 */
	initialiseFields(collideField, streamField, flagField, xlength);

	for(t=0; t < timesteps; t++){
	  	
		double *swap=NULL;
		doStreaming(collideField, streamField, flagField, xlength);

		swap = collideField;
		collideField = streamField;
		streamField = swap;
		doCollision(collideField,flagField,&tau,xlength);
		
		treatBoundary(collideField,flagField,velocityWall,xlength);

		if (t%timestepsPerPlotting==0){
			writeVtkOutput(collideField, flagField, argv[1], t, xlength);
		}


	}
	
	
	float seconds = (float)(clock() - start) / CLOCKS_PER_SEC;
	printf("MLUPS %f\n", (float)(CUBE(xlength+2) / 1000000) / seconds);	
	
	
	/* 
	 * Free the memory allocated with malloc() 
	 */
	free(collideField);
	free(streamField);
	free(flagField);
	
	return 0;
}





#endif

