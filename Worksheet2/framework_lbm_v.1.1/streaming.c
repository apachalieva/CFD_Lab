#include "streaming.h"
#include "LBDefinitions.h"
#include "helper.h"
/*#include<stdio.h>*/

void doStreaming( double *collideField, double *streamField, int *flagField, int xlength ){
 
  /* 
   *This worksheet is restricted to the D3Q19 model - three-dimentional model with 19 lattice velocities ci 
   */		
    unsigned long x, y, z, i;			/* Helper variables for moving through the velocities in our 3D model */
    
    for( z = 1lu; z < xlength + 1; ++z ){
      for( y = 1lu; y < xlength + 1; ++y ){
	for( x = 1lu; x < xlength + 1; ++x ){
	  /* Proceed, only if the state of the cells is FLUID = 0 */
	  
	  if( flagField[ x + (xlength+2) * y + SQ(xlength+2) * z ] == 0 ){
	   for( i = 0lu; i < Q; ++i ){
	     streamField[ Q * ( x + (xlength+2) * y + SQ(xlength+2) * z ) + i ] = collideField[ Q * ( ( x - LATTICEVELOCITIES[i][0] ) + 
												  (xlength+2) * ( y - LATTICEVELOCITIES[i][1] ) + 
												  SQ(xlength+2) * ( z - LATTICEVELOCITIES[i][2] ) ) + 
												  i ];
	    }
	  }
	}
      }
    }
}