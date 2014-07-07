#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double dx,
  double dy,
  double **U,
  double **V,
  double **K,
  double **W,
  double nu,
  int *b,
  int **Flag );


/* fuction for
 * INFLOW boundary condition
 */
void spec_boundary_val( char* problem, int imax, int jmax,double **U, double **V, double **K, double **E, double Re, double dp, double cn, double ylength);

#endif

