#include "boundary.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"
#include "helper.h"

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
  /* TODO Me */
  
  int x, y, z, i, j;
  
  double density, *currentCell;
  int INV_LATTICEVELOCITIES[19][3];
  
  for( i = 0; i < Q; ++i ){
    for( j = 0; j < DIM; ++j ){
      if( LATTICEVELOCITIES[i][j] == 1 || LATTICEVELOCITIES[i][j] == -1 ){
	INV_LATTICEVELOCITIES[i][j] = -LATTICEVELOCITIES[i][j];
      }
      else{
	INV_LATTICEVELOCITIES[i][j] = 0;
      }
    }
  }
/* 
 * Printing the inverse velocities c_i
 */
/*
  for( i = 0lu; i < Q; ++i ){
    printf( "inverse_c[%lu] = ", i );
    for( j = 0lu; j < DIM; ++j ){
      printf( "%d, ", inverse_c[i][j]  );
    }
   printf( "\n " );
  }
*/

  for( z = 1; z < xlength + 1; ++z ){
    for( y = 1; y < xlength + 1; ++y ){
      for( x = 1; x < xlength + 1; ++x ){
	 for( i = 0; i < Q; ++i ){
	    
	  currentCell = &collideField[ Q * ( x + (xlength+2) * y + SQ(xlength+2) * z ) + i ];
	
	  /* NO-SLIP = 1 boundary conditions */
	  if( flagField[ x + (xlength+2) * y + SQ(xlength+2) * z ] == 1 ){
	     *currentCell = collideField[ Q * ( ( x + INV_LATTICEVELOCITIES[i][0] ) + (xlength+2) * ( y + INV_LATTICEVELOCITIES[i][1] ) + SQ(xlength+2) * ( z + INV_LATTICEVELOCITIES[i][2] ) ) + i ];
	     
	  }
	  else if( flagField[  x + (xlength+2) * y + SQ(xlength+2) * z  ] == 2 ){
	     computeDensity( currentCell, &density );
	     
	     *currentCell = collideField[ Q * ( ( x + INV_LATTICEVELOCITIES[i][0] ) + (xlength+2) * ( y + INV_LATTICEVELOCITIES[i][1] ) + SQ(xlength+2) * ( z + INV_LATTICEVELOCITIES[i][2] ) ) + i ] 
					+ 2 * LATTICEWEIGHTS[i] * density * ( LATTICEVELOCITIES[i][0]*wallVelocity[0]+LATTICEVELOCITIES[i][1]*wallVelocity[1]+LATTICEVELOCITIES[i][2]*wallVelocity[2])*3;
	  } 
	}
      }
    }
  } 
}

