#include "sor.h"
#include <math.h>

void sor(
  double omg,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res,
  double dp,
  char *problem			/* type of problem. It is needed to specify the pressure boundary conditions */
) {
  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  rloc = rloc/(imax*jmax);
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;





  	  /* set homogenous neumann boundary conditions for pressure in the horizontal walls */
	  for(i = 1; i <= imax; i++) {
		P[i][0] = P[i][1];
		P[i][jmax+1] = P[i][jmax];
	  }

	if (strcmp(problem,"shear")==0){
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



}

