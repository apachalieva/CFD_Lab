#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  int *boundrs,		/* vector representing different types of boundaries */
  int **Flag		/* Flag field, describing the geometry of the problem */
);

/* setting the boundary values for k or eps */
void boundaryvalues_k_eps(int imax, int jmax, double** k,int* boundaries,int** Flag);

/* fuction for
 * INFLOW boundary condition
 */
void spec_boundary_val( char* problem, int imax, int jmax,double **U, double **V, double **k, double **eps, double Re, double dp, double h);

#endif
