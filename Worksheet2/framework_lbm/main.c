#ifndef _MAIN_C_
#define _MAIN_C_

#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "LBDefinitions.h"



int main (int argc, char *argv[]){
	double *collideField=NULL, *streamField=NULL, tau, velocityWall[3];
	int t, xlength, timesteps, timestepsPerPlotting, *flagField=NULL;

	readParameters( &xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc, argv);

	collideField = malloc(Q * CUBE(xlength+2) * sizeof(double));
	streamField  = malloc(Q * CUBE(xlength+2) * sizeof(double));
	flagField    = malloc(    CUBE(xlength+2) * sizeof(double));


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

	return 0;
}





#endif

