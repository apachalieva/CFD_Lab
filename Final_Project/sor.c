#include "sor.h"
#include "helper.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  const int fluid_cells,
  double **P,
  double **RS,
  int    **Flag,
  double *res,
  char* problem,
  double dp
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));


  
  
  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
      for(j = 1; j<=jmax; j++) {
		  /* SOR computation limited to the fluid cells ( C_F - fluid cell )*/

			if ( IS_FLUID(Flag[ i ][ j ]) )
				P[i][j] =  (1.0-omg)*P[i][j] + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);

      }
  }


  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++)
      for(j = 1; j <= jmax; j++)
	  /* Residual computation limited to the fluid cells ( C_F - fluid cell ) */
	      rloc += (IS_FLUID(Flag[ i ][ j ])) *
	      	  	  ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
		          ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);


  /* Residual devided only by the number of fluid cells instead of imax*jmax */
  rloc = rloc/((double)fluid_cells);
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;




  /* set homogenous neumann boundary conditions for pressure in the horizontal walls */
   for(i = 1; i <= imax; i++) {
 	P[i][0] = P[i][1];
 	P[i][jmax+1] = P[i][jmax];
   }

 if (strcmp(problem,"shear")==0 || strcmp(problem,"karman")==0|| strcmp(problem,"step")==0){
 	/* pressure differece driven flow */
 	for (j=1; j<=jmax; j++){
 		P[0][j]=2.0*dp-P[1][j]; 					/* set left pressure dirichlet condition to p_w = dp */
 		P[imax+1][j]=-P[imax][j]; 					/* set right pressure dirichlet condition to p_w = 0 */
 	}
 } else {
 	  /* set homogenous neumann boundary conditions for pressure in the vertical walls */
 	  for(j = 1; j <= jmax; j++) {
 	    P[0][j] = P[1][j];
 	    P[imax+1][j] = P[imax][j];
 	  }
 }

   /* Boundary conditions for the obstacle cells */
   for( i = 1; i <= imax; i++ ){
 	for( j = 1; j <= jmax; j++ ){
 	    if( IS_BOUNDARY(Flag[i][j])){

 	 		/* Boundary conditions for obstacles with North-Eastern fluid cell */
 	 		if( ( Flag[ i ][ j ] & B_NE ) == B_NE ){
 	 		    P[ i ][ j ] = 0.5*( P[ i+1 ][ j ] + P[ i ][ j+1 ] );
 	 		} else

 	 		/* Boundary conditions for obstacles with North-Western fluid cell */
 	 		if( ( Flag[ i ][ j ] & B_NW ) == B_NW ){
 	 		    P[ i ][ j ] = 0.5*( P[ i-1 ][ j ] + P[ i ][ j+1 ] );
 	 		} else

 	 		/* Boundary conditions for obstacles with South-Eastern fluid cell */
 	 		if( ( Flag[ i ][ j ] & B_SE ) == B_SE ){
 	 		    P[ i ][ j ] = 0.5*( P[ i+1 ][ j ] + P[ i ][ j-1 ] );
 	 		} else

 	 		/* Boundary conditions for obstacles with South-Western fluid cell */
 	 		if( ( Flag[ i ][ j ] & B_SW ) == B_SW ){
 	 		    P[ i ][ j ] = 0.5*( P[ i-1 ][ j ] + P[ i ][ j-1 ] );
 	 		} else
 		/* Boundary conditions for obstacles with Northern fluid cell */
 		if( ( Flag[ i ][ j ] & B_N ) == B_N ){
 		    P[ i ][ j ] = P[ i ][ j+1 ];
 		} else
 		/* Boundary conditions for obstacles with Southern fluid cell */
 		if( ( Flag[ i ][ j ] & B_S ) == B_S ){
 		    P[ i ][ j ] = P[ i ][ j-1 ];
 		} else
 		/* Boundary conditions for obstacles with Western fluid cell */
 		if( ( Flag[ i ][ j ] & B_W ) == B_W ){
 		    P[ i ][ j ] = P[ i-1 ][ j ];
 		} else
 		/* Boundary conditions for obstacles with Eastern fluid cell */
 		if( ( Flag[ i ][ j ] & B_E ) == B_E ){
 		    P[ i ][ j ] = P[ i+1 ][ j ];
 		}

 	    }
 	}
     }
}

