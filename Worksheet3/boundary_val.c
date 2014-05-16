/*
 * boundary_val.c
 *
 *  Created on: May 11, 2014
 *      Author: mauro
 */
#include <stdio.h>
#include <string.h>
#include "boundary_val.h"

void no_slip(int imax, int jmax, double** U, double** V,int c){

	int i,j;

	switch (c){
	case 0:{	/* left boundary */
		for (j=1; j<=jmax; j++){
			U[0][j]=0;						/* set zero to u-values exactly at the boundary */
			V[0][j] = -V[1][j];				/* set the v-values not at the boundary so that the mean in the boundary is zero */
			}
		}
	break;
	case 1:{	/* right boundary */
		for (j=1; j<=jmax; j++){
			U[imax][j]=0;					/* set zero to u-values exactly at the boundary */
			V[imax+1][j] = -V[imax][j];		/* set the v-values not at the boundary so that the mean in the boundary is zero */

			}
		}
	break;
	case 2:{	/* bottom boundary */
		for (i=1; i<=imax; i++){
			U[i][0]=-U[i][1];				/* set the u-values not at the boundary so that the mean in the boundary is zero */
			V[i][0]=0;						/* set zero to v-values exactly at the boundary */
			}

		}
	break;
	case 3:{	/* top boundary */
		for (i=1; i<=imax; i++){
			U[i][jmax+1]=-U[i][jmax];		/* set the u-values not at the boundary so that the mean in the boundary is zero */
			V[i][jmax]=0;						/* set zero to v-values exactly at the boundary */
			}

		}
	break;
	}

}

void free_slip(int imax, int jmax, double** U, double** V,int c){

	int i,j;

	switch (c){
	case 0:{	/* left boundary */
		for (j=1; j<=jmax; j++){
			U[0][j]=0;						/* set zero to u-values exactly at the boundary */
			V[0][j] = V[1][j];				/* set the v-values not at the boundary so that the derivative normal to the boundary is zero */
			}
		}
	break;
	case 1:{	/* right boundary */
		for (j=1; j<=jmax; j++){
			U[imax][j]=0;					/* set zero to u-values exactly at the boundary */
			V[imax+1][j] = V[imax][j];		/* set the v-values not at the boundary so that the derivative normal to the boundary is zero */
			}
		}
	break;
	case 2:{	/* bottom boundary */
		for (i=1; i<=imax; i++){
			U[i][0]=U[i][1];				/* set the u-values not at the boundary so that the derivative normal to the boundary is zero */
			V[i][0]=0;						/* set zero to v-values exactly at the boundary */
			}

		}
	break;
	case 3:{	/* top boundary */
		for (i=1; i<=imax; i++){
			U[i][jmax+1]=U[i][jmax];		/* set the u-values not at the boundary so that the derivative normal to the boundary is zero */
			V[i][jmax]=0;					/* set zero to v-values exactly at the boundary */
			}

		}
	break;
	}

}

void outflow(int imax, int jmax, double** U, double** V,int c){

	int i,j;

	switch (c){
	case 0:{	/* left boundary */
		for (j=1; j<=jmax; j++){
			U[0][j]=U[1][j];				/* set the normal derivative of u to zero */
			V[0][j] = V[1][j];				/* set the normal derivative of v to zero */
			}
		}
	break;
	case 1:{	/* right boundary */
		for (j=1; j<=jmax; j++){
			U[imax][j]=U[imax-1][j];				/* set the normal derivative of u to zero */
			V[imax+1][j] = V[imax][j];				/* set the normal derivative of v to zero */
			}
		}
	break;
	case 2:{	/* bottom boundary */
		for (i=1; i<=imax; i++){
			U[i][0]=U[i][1];				/* set the normal derivative of u to zero */
			V[i][0] = V[i][1];				/* set the normal derivative of v to zero */
			}

		}
	break;
	case 3:{	/* top boundary */
		for (i=1; i<=imax; i++){
			U[i][jmax+1]=U[i][jmax];				/* set the normal derivative of u to zero */
			V[i][jmax] = V[i][jmax-1];				/* set the normal derivative of v to zero */
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
  int *b)
{
	int c,bound_now;

	for(c=0;c<4;c++)
		bound_now=b[c];
		/* treating different cases of boundaries
		 * inflow treated separately
		 * */
		switch(bound_now){
		case 1: no_slip( imax, jmax, U, V, c);
		break;
		case 3: outflow( imax, jmax, U, V, c);
		break;
		default: free_slip( imax, jmax, U, V, c);
		break;
		}

}


/* fuction for
 * INFLOW boundary condition
 */
void spec_boundary_val( char* problem, int imax, int jmax, double **U, double **V, double **P, double u_in, double v_in, double dp){
/* supposing the three different problems:
 * karman = Karman vortex street
 * shear = plane shear flow
 * step = flow over a step
 * we deal with these problems
 */

	int j;

	if (strcmp(problem,"shear")==0){
		/* pressure differece driven flow */

		for (j=1; j<=jmax; j++){
			P[0][j]=2.0*dp-P[1][j]; 					/* set left pressure dirichlet condition to p_w = dp */
			P[imax+1][j]=-P[imax][j]; 					/* set right pressure dirichlet condition to p_w = 0 */
		}
	} else if (strcmp(problem,"karman")==0 || strcmp(problem,"step")==0){
		/* printf("setting the left boundary to inflow velocity : u=u_in, v=v_in;\n"); */
		for (j=1; j<=jmax; j++){
			U[0][j]=u_in;
			V[0][j]=-V[1][j] + 2.0 * v_in; 				/* setting the average equal to v_in */
		}
	}

}


