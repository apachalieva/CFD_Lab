/*
 * boundary_val.c
 *
 */
#include <stdio.h>
#include <string.h>
#include "boundary_val.h"
#include "helper.h"
#include <math.h>

/**
 * No-slip boundary conditions
 */
void no_slip(
  int 	imax, 
  int 	jmax, 
  double **U, 
  double **V, 
  double **K, 
  double **E, 
  int    c
){
	int i,j;

	switch (c){
	case 0:{	/* left boundary */
		for (j=1; j<=jmax; j++){
			U[0][j] = .0;		/* set zero to u-values exactly at the boundary */
			V[0][j] = -V[1][j];	/* set the v-values not at the boundary so that the mean in the boundary is zero */
			K[0][j] = 0.000001;	/* set 0.000001 to k-values at the boundary */
			E[0][j] = E[1][j];	/* homogeneous Neumann */
			}
		}
	break;
	case 1:{	/* right boundary */
		for (j=1; j<=jmax; j++){
			U[imax][j]   =.0;		/* set zero to u-values exactly at the boundary */
			V[imax+1][j] = -V[imax][j];	/* set the v-values not at the boundary so that the mean in the boundary is zero */
			K[imax+1][j] = 0.000001;	/* set 0.000001 to k-values at the boundary */
			E[imax+1][j] = E[imax][j];	/* homogeneous Neumann */
			}
		}
	break;
	case 2:{	/* bottom boundary */
		for (i=1; i<=imax; i++){
			U[i][0] = -U[i][1];	/* set the u-values not at the boundary so that the mean in the boundary is zero */
			V[i][0] = .0;		/* set zero to v-values exactly at the boundary */
			K[i][0] = 0.000001;	/* set 0.000001 to k-values at the boundary */
			E[i][0] = E[i][1];	/* homogeneous Neumann */
			}
		}
	break;
	case 3:{	/* top boundary */
		for (i=1; i<=imax; i++){
			U[i][jmax+1] = -U[i][jmax];	/* set the u-values not at the boundary so that the mean in the boundary is zero */
			V[i][jmax]   = .0;		/* set zero to v-values exactly at the boundary */
			K[i][jmax+1] = 0.000001;	/* set 0.000001 to k-values at the boundary */
			E[i][jmax+1] = E[i][jmax];	/* homogeneous Neumann */
			}
		}
	break;
	}

}

/**
 * Free-slip boundary conditions
 */
void free_slip(
  int    imax, 
  int    jmax, 
  double **U, 
  double **V, 
  double **K, 
  double **E, 
  int    c
){
	int i,j;

	switch (c){
	case 0:{	/* left boundary */
		for (j=1; j<=jmax; j++){
			U[0][j] = .0;		/* set zero to u-values exactly at the boundary */
			V[0][j] = V[1][j];	/* set the v-values not at the boundary so that the derivative normal to the boundary is zero */
			K[0][j] = K[1][j];	/* homogeneous Neumann */
			E[0][j] = E[1][j];	/* homogeneous Neumann */
			}
		}
	break;
	case 1:{	/* right boundary */
		for (j=1; j<=jmax; j++){
			U[imax][j]   = 0;		/* set zero to u-values exactly at the boundary */
			V[imax+1][j] = V[imax][j];	/* set the v-values not at the boundary so that the derivative normal to the boundary is zero */
			K[imax+1][j] = K[imax][j];	/* homogeneous Neumann */
			E[imax+1][j] = E[imax][j];	/* homogeneous Neumann */
			}
		}
	break;
	case 2:{	/* bottom boundary */
		for (i=1; i<=imax; i++){
			U[i][0] = U[i][1];	/* set the u-values not at the boundary so that the derivative normal to the boundary is zero */
			V[i][0] = 0;		/* set zero to v-values exactly at the boundary */
			K[i][0] = K[i][1];	/* homogeneous Neumann */
			E[i][0] = E[i][1];	/* homogeneous Neumann */
			}
		}
	break;
	case 3:{	/* top boundary */
		for (i=1; i<=imax; i++){
			U[i][jmax+1] = U[i][jmax];	/* set the u-values not at the boundary so that the derivative normal to the boundary is zero */
			V[i][jmax]   = 0;		/* set zero to v-values exactly at the boundary */
			K[i][jmax+1] = K[i][jmax];	/* homogeneous Neumann */
			E[i][jmax+1] = E[i][jmax];	/* homogeneous Neumann */
			}

		}
	break;
	}

}

/**
 * Outflow boundary conditions
 */
void outflow(
  int    imax, 
  int    jmax, 
  double **U, 
  double **V, 
  double **K, 
  double **E, 
  int    c
){
	int i,j;

	switch (c){
	case 0:{	/* left boundary */
		for (j=1; j<=jmax; j++){
			U[0][j] = U[1][j];	/* set the normal derivative of u to zero */
			V[0][j] = V[1][j];	/* set the normal derivative of v to zero */
			K[0][j] = K[1][j];	/* set the normal derivative of k to zero */
			E[0][j] = E[1][j];	/* set the normal derivative of eps to zero */
			}
		}
	break;
	case 1:{	/* right boundary */
		for (j=1; j<=jmax; j++){
			U[imax][j]   = U[imax-1][j];	/* set the normal derivative of u to zero */
			V[imax+1][j] = V[imax][j];	/* set the normal derivative of v to zero */
			K[imax+1][j] = K[imax][j];	/* set the normal derivative of k to zero */
			E[imax+1][j] = E[imax][j];	/* set the normal derivative of eps to zero */
			}
		}
	break;
	case 2:{	/* bottom boundary */
		for (i=1; i<=imax; i++){
			U[i][0] = U[i][1];	/* set the normal derivative of u to zero */
			V[i][0] = V[i][1];	/* set the normal derivative of v to zero */
			K[i][0] = K[i][1];	/* set the normal derivative of k to zero */
			E[i][0] = E[i][1];	/* set the normal derivative of eps to zero */
			}
		}
	break;
	case 3:{	/* top boundary */
		for (i=1; i<=imax; i++){
			U[i][jmax+1] = U[i][jmax];	/* set the normal derivative of u to zero */
			V[i][jmax]   = V[i][jmax-1];	/* set the normal derivative of v to zero */
			K[i][jmax+1] = K[i][jmax];	/* set the normal derivative of k to zero */
			E[i][jmax+1] = E[i][jmax];	/* set the normal derivative of eps to zero */
			}
		}
	break;
	}

}

/**
 * Set the boundary conditions depending on the chosen model 
 */
void boundaryvalues(
  int    imax,
  int    jmax,
  double **U,
  double **V,
  double **K,
  double **E,
  int    *b, 
  int    **Flag )
{
	int i, j;
	int c, bound_now;

	for( c = 0; c < 4; c++ ){
		bound_now = b[c];
		/* treating different cases of boundaries
		 * inflow treated separately
		 * */
		switch(bound_now){
		case 1: no_slip( imax, jmax, U, V, K, E, c);
		break;
		case 3: outflow( imax, jmax, U, V, K, E, c);
		break;
		default: free_slip( imax, jmax, U, V, K, E, c);
		break;
		}
	}
	
	/* Boundary conditions for the obstacle cells */
	for( i = 1; i <= imax; i++ )
	    for( j = 1; j <= jmax; j++ )
			if( IS_BOUNDARY(Flag[i][j]) ){
				  /* Boundary conditions for obstacles with North-Eastern fluid cell */
				  if( ( Flag[ i ][ j ] & B_NE ) == B_NE ){
					  U[ i ][ j ] 	= .0;
					  V[ i ][ j ] 	= .0;
					  U[ i-1 ][ j ] = -U[ i-1 ][ j+1 ];
					  V[ i ][ j-1 ] = -V[ i+1 ][ j-1 ];
					  K[ i ][ j ] 	= 0.00001;
			 		  E[ i ][ j ] 	= 0.5*( E[ i+1 ][ j ] + E[ i ][ j+1 ] );
				  } else

				  /* Boundary conditions for obstacles with North-Western fluid cell */
				  if( ( Flag[ i ][ j ] & B_NW ) == B_NW ){
					  U[ i-1 ][ j ] = .0;
					  V[ i ][ j ] 	= .0;
					  U[ i ][ j ] 	= -U[ i ][ j+1 ];
					  V[ i ][ j-1 ] = -V[ i-1 ][ j-1 ];
					  K[ i ][ j ] 	= 0.00001;
					  E[ i ][ j ] 	= 0.5*( E[ i-1 ][ j ] + E[ i ][ j+1 ] );
				  } else

				  /* Boundary conditions for obstacles with South-Eastern fluid cell */
				  if( ( Flag[ i ][ j ] & B_SE ) == B_SE ){
					  U[ i ][ j ] 	= .0;
					  V[ i ][ j-1 ] = .0;
					  U[ i-1 ][ j ] = -U[ i-1 ][ j-1 ];
					  V[ i ][ j ]	= -V[ i+1 ][ j ];
					  K[ i ][ j ] 	= 0.00001;
					  E[ i ][ j ] 	= 0.5*( E[ i+1 ][ j ] + E[ i ][ j-1 ] );
				  } else

				  /* Boundary conditions for obstacles with South-Western fluid cell */
				  if( ( Flag[ i ][ j ] & B_SW ) == B_SW ){
					  U[ i-1 ][ j ] = .0;
					  V[ i ][ j-1 ] = .0;
					  U[ i ][ j ] 	= -U[ i ][ j-1 ];
					  V[ i ][ j ] 	= -V[ i-1 ][ j ];
					  K[ i ][ j ] 	= 0.00001;
					  E[ i ][ j ] 	= 0.5*( E[ i-1 ][ j ] + E[ i ][ j-1 ] );
				  } else


				  /* Boundary conditions for obstacles with Northern fluid cell */
				  if( ( Flag[ i ][ j ] & B_N ) == B_N ){
					  V[ i ][ j ] 	= .0;
					  U[ i ][ j ] 	= -U[ i ][ j+1 ];
					  U[ i-1 ][ j ] = -U[ i-1 ][ j+1 ];
					  K[ i ][ j ] 	= 0.00001;
					  E[ i ][ j ] 	= E[ i ][ j+1 ];
				  } else
				  /* Boundary conditions for obstacles with Southern fluid cell */
			  if(( Flag[ i ][ j ] & B_S ) == B_S ){
					  V[ i ][ j-1 ] = .0;
					  U[ i ][ j ] 	= -U[ i ][ j-1 ];
					  U[ i-1 ][ j ] = -U[ i-1 ][ j-1 ];
					  K[ i ][ j ]	= 0.00001;
			 		  E[ i ][ j ] 	= E[ i ][ j-1 ];
				  }else
				  /* Boundary conditions for obstacles with Western fluid cell */
				  if( ( Flag[ i ][ j ] & B_W ) == B_W ){
					  U[ i-1 ][ j ] = .0;
					  V[ i ][ j ] 	= -V[ i-1 ][ j ];
					  V[ i ][ j-1 ] = -V[ i-1 ][ j-1 ];
					  K[ i ][ j ] 	= 0.00001;
			 		  E[ i ][ j ] 	= E[ i-1 ][ j ];
				  } else
				  /* Boundary conditions for obstacles with Eastern fluid cell */
				  if( ( Flag[ i ][ j ] & B_E ) == B_E ){
					  U[ i ][ j ]	= .0;
					  V[ i ][ j ] 	= -V[ i+1 ][ j ];
					  V[ i ][ j-1 ] = -V[ i+1 ][ j-1 ];
					  K[ i ][ j ] 	= 0.00001;
					  E[ i ][ j ] 	= E[ i+1 ][ j ];
				  }


			}
}


/* fuction for
 * INFLOW boundary condition
 */
void spec_boundary_val( 
  char   *problem, 
  int    imax, 
  int    jmax, 
  double **U, 
  double **V, 
  double **K, 
  double **E, 
  double Re, 
  double dp, 
  double cn, 
  double ylength
){
	int j;
	for (j=1; j<=jmax; j++){

		if(dp==0){
			/* inflow velocity */
			U[0][j]=1.0;			/* normalized inflow velocity */
			V[0][j]=-V[1][j];		/* setting the average equal to 0 */
		}else{
			/* pressure driven flow */			
			U[0][j]= U[1][j];		/* homogeneous neumann conditions */
			V[0][j]= V[1][j];		/* homogeneous neumann conditions */
		}

		K[0][j] = K[1][j];			/* homogeneous neumann conditions */
		E[0][j] = E[1][j];			/* homogeneous neumann conditions */
	}
}

