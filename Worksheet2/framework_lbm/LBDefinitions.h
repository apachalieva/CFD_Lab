#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#define Q 19 		/* Number of velocities */
#define DIM 3		/* Number of dimensions */

  static const int LATTICEVELOCITIES[19][3] = {   {0,-1,-1}, {-1,0,-1}, {0,0,-1}, {1,0,-1}, {0,1,-1}, {-1,-1,0}, {0,-1,0}, {1,-1,0}, {-1,0,0}, {0,0,0}, {1,0,0}, {-1,1,0}, {0,1,0}, {1,1,0}, {0,-1,1}, {-1,0,1}, {0,0,1}, {1,0,1}, {0,1,1} };
  static const double LATTICEWEIGHTS[19] = { 1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 2.0/36.0, 1.0/36.0, 2.0/36.0, 12.0/36.0,  2.0/36.0,  1.0/36.0,  2.0/36.0,  1.0/36.0, 1.0/36.0,  1.0/36.0,  2.0/36.0,  1.0/36.0,  1.0/36.0 };
  static const double C_S = 0.577350269; /* 1/sqrt(3) */

#endif

