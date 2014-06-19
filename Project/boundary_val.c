/*
 * boundary_val.c
 *
 *  Created on: Apr 16, 2014
 *      Author: davide
 */


void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
){
	int i,j;

	for(i=1; i<=imax; i++){
		V[i][0] = 0;		/* set zero to v-values exactly at the boundary */
		V[i][jmax] = 0;		/* set zero to v-values exactly at the boundary */

		U[i][0] = -U[i][1];				/* set the u-values not at the boundary so that the mean in the boundary is zero */
		U[i][jmax+1] = 2.0-U[i][jmax];

	}

	for(j=1; j<=jmax; j++){
		U[0][j] = 0;		/* set zero to u-values exactly at the boundary */
		U[imax][j] = 0;		/* set zero to u-values exactly at the boundary */

		V[0][j] = -V[1][j];				/* set the v-values not at the boundary so that the mean in the boundary is zero */
		V[imax+1][j] = -V[imax][j];		/* set the v-values not at the boundary so that the mean in the boundary is zero */
	}

}

