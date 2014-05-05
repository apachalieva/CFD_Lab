#include "boundary.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "initLB.h"


inline int inv_velocity(int i){
	return Q-i-1;
}

inline int getCell(int x, int y, int z, int xlength){
	 return Q * (x + (xlength+2) * y + SQ(xlength+2) * z);
}

inline void treatCase(double *collideField, int* flagField, int x0, int x1, int y0, int y1, int z0, int z1, int const * const vel, int n, int xlength, const double * const wallVelocity){
	int x, y, z, i;
	double density;

	for(z=z0; z<=z1; z++)
		for(y=y0; y<=y1; y++)
			for(x=x0; x<=x1; x++)
					for(i=0; i<n; i++){
						collideField[getCell(x, y, z, xlength) + vel[i]] =
								collideField[ getCell(	LATTICEVELOCITIES[vel[i]][0] + x,
														LATTICEVELOCITIES[vel[i]][1] + y,
														LATTICEVELOCITIES[vel[i]][2] + z,
														xlength)
														+ inv_velocity(vel[i])
											];


						if(flagField[x + (xlength+2) * y + SQ(xlength+2) * z] == MOVING_WALL){
							computeDensity(collideField+ getCell(	LATTICEVELOCITIES[vel[i]][0] + x,
													LATTICEVELOCITIES[vel[i]][1] + y,
													LATTICEVELOCITIES[vel[i]][2] + z,
													xlength),
											&density);

							collideField[getCell(x, y, z, xlength) + vel[i]] += 2.0 * LATTICEWEIGHTS[vel[i]] * density * DOTP(LATTICEVELOCITIES[vel[i]], wallVelocity) / SQ(C_S);

						}


					}
}


void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){

	int vel[19];


	/* z=0 boundary */

	/* 			y=0 boundary */
	vel[0] = 18;
	treatCase(collideField, flagField, 1, xlength, 0, 0, 0, 0, vel, 1, xlength, wallVelocity );

	/* 			y=xlength+1 boundary */
	vel[0] = 14;
	treatCase(collideField, flagField, 1, xlength, xlength+1, xlength+1, 0, 0, vel, 1, xlength, wallVelocity );

	/* 			x=0 boundary */
	vel[0] = 17;
	treatCase(collideField, flagField, 0, 0, 1, xlength, 0, 0, vel, 1, xlength, wallVelocity );

	/* 			x=xlength+1 boundary */
	vel[0] = 15;
	treatCase(collideField, flagField, xlength+1, xlength+1, 1, xlength, 0, 0, vel, 1, xlength, wallVelocity );

	/* 			inner face boundary */
	vel[0] = 14;
	vel[1] = 15;
	vel[2] = 16;
	vel[3] = 17;
	vel[4] = 18;
	treatCase(collideField, flagField, 1, xlength, 1, xlength, 0, 0, vel, 5, xlength, wallVelocity );

	/* z=xlength+1 boundary */

	/* 			y=0 boundary */
	vel[0] = 4;
	treatCase(collideField, flagField, 1, xlength, 0, 0, xlength+1, xlength+1, vel, 1, xlength, wallVelocity );

	/* 			y=xlength+1 boundary */
	vel[0] = 0;
	treatCase(collideField, flagField, 1, xlength, xlength+1, xlength+1, xlength+1, xlength+1, vel, 1, xlength, wallVelocity );

	/* 			x=0 boundary */
	vel[0] = 3;
	treatCase(collideField, flagField, 0, 0, 1, xlength, xlength+1, xlength+1, vel, 1, xlength, wallVelocity );

	/* 			x=xlength+1 boundary */
	vel[0] = 1;
	treatCase(collideField, flagField, xlength+1, xlength+1, 1, xlength, xlength+1, xlength+1, vel, 1, xlength, wallVelocity );

	/* 			inner face boundary */
	vel[0] = 0;
	vel[1] = 1;
	vel[2] = 2;
	vel[3] = 3;
	vel[4] = 4;
	treatCase(collideField, flagField, 1, xlength, 1, xlength, xlength+1, xlength+1, vel, 5, xlength, wallVelocity );


	/* x=0 boundary */

	/* 			y=0 boundary */
	vel[0] = 13;
	treatCase(collideField, flagField, 0, 0, 0, 0, 1, xlength, vel, 1, xlength, wallVelocity );

	/* 			y=xlength+1 boundary */
	vel[0] = 7;
	treatCase(collideField, flagField, 0, 0, xlength+1, xlength+1, 1, xlength, vel, 1, xlength, wallVelocity );

	/* 			z=0 boundary */
	vel[0] = 17;
	treatCase(collideField, flagField, 0, 0, 1, xlength, 0, 0, vel, 1, xlength, wallVelocity );

	/* 			z=xlength+1 boundary */
	vel[0] = 3;
	treatCase(collideField, flagField, 0, 0, 1, xlength, xlength+1, xlength+1, vel, 1, xlength, wallVelocity );

	/* 			inner face boundary */
	vel[0] = 13;
	vel[1] = 7;
	vel[2] = 17;
	vel[3] = 3;
	vel[4] = 10;
	treatCase(collideField, flagField, 0, 0, 1, xlength, 1, xlength, vel, 5, xlength, wallVelocity );



	/* x=xlength+1 boundary */

	/* 			y=0 boundary */
	vel[0] = 11;
	treatCase(collideField, flagField, xlength+1, xlength+1, 0, 0, 1, xlength, vel, 1, xlength, wallVelocity );

	/* 			y=xlength+1 boundary */
	vel[0] = 5;
	treatCase(collideField, flagField, xlength+1, xlength+1, xlength+1, xlength+1, 1, xlength, vel, 1, xlength, wallVelocity );

	/* 			z=0 boundary */
	vel[0] = 15;
	treatCase(collideField, flagField, xlength+1, xlength+1, 1, xlength, 0, 0, vel, 1, xlength, wallVelocity );

	/* 			z=xlength+1 boundary */
	vel[0] = 1;
	treatCase(collideField, flagField, xlength+1, xlength+1, 1, xlength, xlength+1, xlength+1, vel, 1, xlength, wallVelocity );

	/* 			inner face boundary */
	vel[0] = 11;
	vel[1] = 5;
	vel[2] = 15;
	vel[3] = 1;
	vel[4] = 8;
	treatCase(collideField, flagField, 0, 0, 1, xlength, 1, xlength, vel, 5, xlength, wallVelocity );

	/* y=0 boundary */

	/* 			x=0 boundary */
	vel[0] = 13;
	treatCase(collideField, flagField, 0, 0, 0, 0, 1, xlength, vel, 1, xlength, wallVelocity );

	/* 			x=xlength+1 boundary */
	vel[0] = 11;
	treatCase(collideField, flagField, xlength+1, xlength+1, 0, 0, 1, xlength, vel, 1, xlength, wallVelocity );

	/* 			z=0 boundary */
	vel[0] = 18;
	treatCase(collideField, flagField, 1, xlength, 0, 0, 0, 0, vel, 1, xlength, wallVelocity );

	/* 			z=xlength+1 boundary */
	vel[0] = 4;
	treatCase(collideField, flagField, 1, xlength, 0, 0, xlength+1, xlength+1, vel, 1, xlength, wallVelocity );

	/* 			inner face boundary */
	vel[0] = 13;
	vel[1] = 11;
	vel[2] = 18;
	vel[3] = 4;
	vel[4] = 12;
	treatCase(collideField, flagField, 1, xlength, 0, 0, 1, xlength, vel, 5, xlength, wallVelocity );

	/* y=xlength+1 boundary */

	/* 			x=0 boundary */
	vel[0] = 7;
	treatCase(collideField, flagField, 0, 0, xlength+1, xlength+1, 1, xlength, vel, 1, xlength, wallVelocity );

	/* 			x=xlength+1 boundary */
	vel[0] = 5;
	treatCase(collideField, flagField, xlength+1, xlength+1, xlength+1, xlength+1, 1, xlength, vel, 1, xlength, wallVelocity );

	/* 			z=0 boundary */
	vel[0] = 14;
	treatCase(collideField, flagField, 1, xlength, xlength+1, xlength+1, 0, 0, vel, 1, xlength, wallVelocity );

	/* 			z=xlength+1 boundary */
	vel[0] = 0;
	treatCase(collideField, flagField, 1, xlength, xlength+1, xlength+1, xlength+1, xlength+1, vel, 1, xlength, wallVelocity );

	/* 			inner face boundary */
	vel[0] = 7;
	vel[1] = 5;
	vel[2] = 14;
	vel[3] = 0;
	vel[4] = 6;
	treatCase(collideField, flagField, 1, xlength, xlength+1, xlength+1, 1, xlength, vel, 5, xlength, wallVelocity );
}
