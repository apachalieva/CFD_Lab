/*
 * boundary_val.c
 *
 *  Created on: May 11, 2014
 *      Author: mauro
 */
#include <stdio.h>
#include <string.h>
#include "boundary_val.h"
#include "helper.h"
#include <math.h>

void no_slip(int imax, int jmax, double** U, double** V, double** K, double **E, int c){

	int i,j;

	switch (c){
	case 0:{	/* left boundary */
		for (j=1; j<=jmax; j++){
			U[0][j]=.0;						/* set zero to u-values exactly at the boundary */
			V[0][j] = -V[1][j];				/* set the v-values not at the boundary so that the mean in the boundary is zero */
			/*K[0][j]=2*0.0001-K[1][j];*/
			K[0][j]=0.000001;
			E[0][j]=E[1][j];
			}
		}
	break;
	case 1:{	/* right boundary */
		for (j=1; j<=jmax; j++){
			U[imax][j]=.0;					/* set zero to u-values exactly at the boundary */
			V[imax+1][j] = -V[imax][j];		/* set the v-values not at the boundary so that the mean in the boundary is zero */
			/*K[imax+1][j]=2*0.0001-K[imax][j];*/
			K[imax+1][j]=0.000001;
			E[imax+1][j]=E[imax][j];
			}
		}
	break;
	case 2:{	/* bottom boundary */
		for (i=1; i<=imax; i++){
			U[i][0]=-U[i][1];				/* set the u-values not at the boundary so that the mean in the boundary is zero */
			V[i][0]=.0;						/* set zero to v-values exactly at the boundary */
			/*K[i][0]=2*0.0001-K[i][1];*/
			K[i][0]=0.000001;
			E[i][0]=E[i][1];
			}

		}
	break;
	case 3:{	/* top boundary */
		for (i=1; i<=imax; i++){
			U[i][jmax+1]=-U[i][jmax];		/* set the u-values not at the boundary so that the mean in the boundary is zero */
			V[i][jmax]=.0;						/* set zero to v-values exactly at the boundary */
			K[i][jmax+1]=0.000001;
			/*K[i][jmax+1]=2*0.0001 -K[i][jmax];*/
			E[i][jmax+1]=E[i][jmax];
			}

		}
	break;
	}

}

void free_slip(int imax, int jmax, double** U, double** V, double** K, double **E, int c){

	int i,j;

	switch (c){
	case 0:{	/* left boundary */
		for (j=1; j<=jmax; j++){
			U[0][j]=.0;						/* set zero to u-values exactly at the boundary */
			V[0][j] = V[1][j];				/* set the v-values not at the boundary so that the derivative normal to the boundary is zero */
			K[0][j]=K[1][j];				/* check this */
			E[0][j]=E[1][j];
			}
		}
	break;
	case 1:{	/* right boundary */
		for (j=1; j<=jmax; j++){
			U[imax][j]=0;					/* set zero to u-values exactly at the boundary */
			V[imax+1][j] = V[imax][j];		/* set the v-values not at the boundary so that the derivative normal to the boundary is zero */
			K[imax+1][j]=K[imax][j];		/* check this */
			E[imax+1][j]=E[imax][j];
			}
		}
	break;
	case 2:{	/* bottom boundary */
		for (i=1; i<=imax; i++){
			U[i][0]=U[i][1];				/* set the u-values not at the boundary so that the derivative normal to the boundary is zero */
			V[i][0]=0;						/* set zero to v-values exactly at the boundary */
			K[i][0]=K[i][1];				/* check this */
			E[i][0]=E[i][1];
			}

		}
	break;
	case 3:{	/* top boundary */
		for (i=1; i<=imax; i++){
			U[i][jmax+1]=U[i][jmax];		/* set the u-values not at the boundary so that the derivative normal to the boundary is zero */
			V[i][jmax]=0;					/* set zero to v-values exactly at the boundary */
			K[i][jmax+1]=K[i][jmax];		/* check this */
			E[i][jmax+1]=E[i][jmax];
			}

		}
	break;
	}

}

void outflow(int imax, int jmax, double** U, double** V, double** K, double **E, int c){

	int i,j;

	switch (c){
	case 0:{	/* left boundary */
		for (j=1; j<=jmax; j++){
			U[0][j] = U[1][j];				/* set the normal derivative of u to zero */
			V[0][j] = V[1][j];				/* set the normal derivative of v to zero */
			K[0][j] = K[1][j];
			E[0][j] = E[1][j];
			}
		}
	break;
	case 1:{	/* right boundary */
		for (j=1; j<=jmax; j++){
			U[imax][j]=U[imax-1][j];				/* set the normal derivative of u to zero */
			V[imax+1][j] = V[imax][j];				/* set the normal derivative of v to zero */
			K[imax+1][j]=K[imax][j];
			E[imax+1][j]=E[imax][j];
			}
		}
	break;
	case 2:{	/* bottom boundary */
		for (i=1; i<=imax; i++){
			U[i][0]=U[i][1];				/* set the normal derivative of u to zero */
			V[i][0] = V[i][1];				/* set the normal derivative of v to zero */
			K[i][0]=K[i][1];
			E[i][0]=E[i][1];
			}

		}
	break;
	case 3:{	/* top boundary */
		for (i=1; i<=imax; i++){
			U[i][jmax+1]=U[i][jmax];				/* set the normal derivative of u to zero */
			V[i][jmax] = V[i][jmax-1];				/* set the normal derivative of v to zero */
			K[i][jmax+1]=K[i][jmax];
			E[i][jmax+1]=E[i][jmax];
			}
		}
	break;
	}

}


void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  double **K,
  double **E,
  int *b, 
  int **Flag )
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
			if( Flag[i][j] < C_F ){
				  /* Boundary conditions for obstacles with North-Eastern fluid cell */
				  if( ( Flag[ i ][ j ] & B_NE ) == B_NE ){
					  U[ i ][ j ] = .0;
					  V[ i ][ j ] = .0;
					  U[ i-1 ][ j ] = -U[ i-1 ][ j+1 ];
					  V[ i ][ j-1 ] = -V[ i+1 ][ j-1 ];
					  /*K[ i ][ j ] = -0.5*( K[ i+1 ][ j ] + K[ i ][ j+1 ] );*/
					  K[ i ][ j ] = -0.5*( K[ i+1 ][ j ] + K[ i ][ j+1 ] ) + 0.0002;
					  /*K[ i ][ j ] = 0.0001;*/
			 		  E[ i ][ j ] = 0.5*( E[ i+1 ][ j ] + E[ i ][ j+1 ] );
				  } else

				  /* Boundary conditions for obstacles with North-Western fluid cell */
				  if( ( Flag[ i ][ j ] & B_NW ) == B_NW ){
					  U[ i-1 ][ j ] = .0;
					  V[ i ][ j ] = .0;
					  U[ i ][ j ] = -U[ i ][ j+1 ];
					  V[ i ][ j-1 ] = -V[ i-1 ][ j-1 ];
					  K[ i ][ j ] = -0.5*( K[ i-1 ][ j ] + K[ i ][ j+1 ] )+0.0002;
					  /*K[ i ][ j ] = 0.0001;*/
					  E[ i ][ j ] = 0.5*( E[ i-1 ][ j ] + E[ i ][ j+1 ] );
				  } else

				  /* Boundary conditions for obstacles with South-Eastern fluid cell */
				  if( ( Flag[ i ][ j ] & B_SE ) == B_SE ){
					  U[ i ][ j ] = .0;
					  V[ i ][ j-1 ] = .0;
					  U[ i-1 ][ j ] = -U[ i-1 ][ j-1 ];
					  V[ i ][ j ] = -V[ i+1 ][ j ];
					  K[ i ][ j ] = -0.5*( K[ i+1 ][ j ] + K[ i ][ j-1 ] )+0.0002;
					  /*K[ i ][ j ] = 0.0001;*/
					  E[ i ][ j ] = 0.5*( E[ i+1 ][ j ] + E[ i ][ j-1 ] );
				  } else

				  /* Boundary conditions for obstacles with South-Western fluid cell */
				  if( ( Flag[ i ][ j ] & B_SW ) == B_SW ){
					  U[ i-1 ][ j ] = .0;
					  V[ i ][ j-1 ] = .0;
					  U[ i ][ j ] = -U[ i ][ j-1 ];
					  V[ i ][ j ] = -V[ i-1 ][ j ];
					  K[ i ][ j ] = -0.5*( K[ i-1 ][ j ] + K[ i ][ j-1 ] )+0.0002;
					  /*K[ i ][ j ] = 0.0001;*/
					  E[ i ][ j ] = 0.5*( E[ i-1 ][ j ] + E[ i ][ j-1 ] );
				  } else


				  /* Boundary conditions for obstacles with Northern fluid cell */
				  if( ( Flag[ i ][ j ] & B_N ) == B_N ){
					  V[ i ][ j ] = .0;
					  U[ i ][ j ] = -U[ i ][ j+1 ];
					  U[ i-1 ][ j ] = -U[ i-1 ][ j+1 ];
					  /*K[ i ][ j ] = -K[ i ][ j+1 ];*/
					  K[ i ][ j ] = 2*0.0001 -K[ i ][ j+1 ];
					  E[ i ][ j ] = E[ i ][ j+1 ];
				  } else
				  /* Boundary conditions for obstacles with Southern fluid cell */
			  if(( Flag[ i ][ j ] & B_S ) == B_S ){
					  V[ i ][ j-1 ] = .0;
					  U[ i ][ j ] = -U[ i ][ j-1 ];
					  U[ i-1 ][ j ] = -U[ i-1 ][ j-1 ];
					  /*K[ i ][ j ] = -K[ i ][ j-1 ];*/
					  K[ i ][ j ] = 0.0002 -K[ i ][ j-1 ];
			 		  E[ i ][ j ] = E[ i ][ j-1 ];
				  }else
				  /* Boundary conditions for obstacles with Western fluid cell */
				  if( ( Flag[ i ][ j ] & B_W ) == B_W ){
					  U[ i-1 ][ j ] = .0;
					  V[ i ][ j ] = -V[ i-1 ][ j ];
					  V[ i ][ j-1 ] = -V[ i-1 ][ j-1 ];
					  /*K[ i ][ j ] = -K[ i-1 ][ j ];*/
					  K[ i ][ j ] = 0.0002 -K[ i-1 ][ j ];
			 		  E[ i ][ j ] = E[ i-1 ][ j ];
				  } else
				  /* Boundary conditions for obstacles with Eastern fluid cell */
				  if( ( Flag[ i ][ j ] & B_E ) == B_E ){
					  U[ i ][ j ] = .0;
					  V[ i ][ j ] = -V[ i+1 ][ j ];
					  V[ i ][ j-1 ] = -V[ i+1 ][ j-1 ];
					  /*K[ i ][ j ] = -K[ i+1 ][ j ];*/
					  K[ i ][ j ] = 0.0002 -K[ i+1 ][ j ];
					  E[ i ][ j ] = E[ i+1 ][ j ];
				  }


			}
}


/* fuction for
 * INFLOW boundary condition
 */
void spec_boundary_val( char* problem, int imax, int jmax, double **U, double **V, double **k, double **eps, double Re, double dp, double cn, double ylength){
/* supposing the three different problems:
 * karman = Karman vortex street
 * shear = plane shear flow
 * steo = flow over a step
 * we deal with these problems
 */
	int j;
	double kin;
	if (strcmp(problem,"karman")==0)
		/* printf("setting the left boundary to velocity : u=1, v=0;\n"); */
		for (j=1; j<=jmax; j++){
			/*U[0][j]=1.0;
			V[0][j]=-V[1][j]; */		/* setting the average equal to 0 */
			U[0][j]= U[1][j];
			V[0][j]= V[1][j];
			kin = 0.003*SQ(U[0][j]);
			k[0][j]=2.*kin-k[1][j];
			eps[0][j]=2.*cn*sqrt(fabs(kin)*SQ(kin))/0.03/ylength-eps[1][j];

		}
	else if(strcmp(problem,"shear")==0)
			for (j=1; j<=jmax; j++){
				U[0][j]= U[1][j];
				V[0][j]= V[1][j]; 		/* setting the average equal to 0 */
				kin = 0.003*SQ(U[0][j]);
				k[0][j]=2.*kin-k[1][j];
				eps[0][j]=2.*cn*sqrt(fabs(kin)*SQ(kin))/0.03/ylength-eps[1][j];
				/*k[0][j]=.5;*/
			}
	else if(strcmp(problem,"step")==0){
				 /* printf("setting the left boundary: lower half = step, upper half: u=1, v=0;\n"); */
				if (jmax%2!=0)
					printf("odd number of cells on the vertical boundary: asymmetric problem!");

				for (j = 1; j<=jmax/2; j++){
						 	U[0][j]= 0.0;
							V[0][j]= -V[1][j]; 		/* setting the average equal to 0 */
							kin = 0.003*SQ(U[0][j]);
							k[0][j]=2.*kin-k[1][j];
							eps[0][j]=2.*cn*sqrt(fabs(kin)*SQ(kin))/0.03/ylength-eps[1][j];
					}

				for (j = (jmax/2+1); j<=jmax; j++){
						 	U[0][j]= 1.0;
							V[0][j]= -V[1][j]; 		/* setting the average equal to 0 */
							kin = 0.003*SQ(U[0][j]);
							k[0][j]=2.*kin-k[1][j];
							eps[0][j]=2.*cn*sqrt(kin*SQ(kin))/0.03/ylength-eps[1][j];
						}
			}
}

