/*
 * boundary_val.c
 *
 *  Created on: Apr 16, 2014
 *      Author: davide
 */


void boundaryvalues(
  int il,
  int ir,
  int jt,
  int jb,
  int imax,
  int jmax,
  double **U,
  double **V
){
	int i,j;

	if(il==1)
		for(j=jb-1; j<=jt+1; j++){
			U[0][j] = 0;		/* set zero to u-values exactly at the boundary */
			V[0][j] = -V[1][j];				/* set the v-values not at the boundary so that the mean in the boundary is zero */
		}

	if(ir==imax)
		for(j=jb-1; j<=jt+1; j++){
			U[imax][j] = 0;		/* set zero to u-values exactly at the boundary */
			V[imax+1][j] = -V[imax][j];		/* set the v-values not at the boundary so that the mean in the boundary is zero */
		}

	if(jb==1)
		for(i=il-1; i<=ir+1; i++){
			V[i][0] = 0;		/* set zero to v-values exactly at the boundary */
			U[i][0] = -U[i][1];				/* set the u-values not at the boundary so that the mean in the boundary is zero */
		}

	if(jt==jmax)
		for(i=il-1; i<=ir+1; i++){
			V[i][jmax] = 0;		/* set zero to v-values exactly at the boundary */
			U[i][jmax+1] = 2.0-U[i][jmax];
		}
}

