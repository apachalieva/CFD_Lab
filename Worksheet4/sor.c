#include "sor.h"
#include <math.h>
#include "parallel.h"
#include <mpi.h>

void print_matr2(char* name, int n, double** M, int il,int  ir,int  jt,int  jb, int ioffs, int joffs){
	int i,j;
	printf("\n\n Matrix %s [%d] \n", name, n);

	for(j = jt+1; j>=jb-1-joffs; j--){
	  for(i = il-1-ioffs; i <= ir+1; i++)
	    	printf("%f\t", M[i][j]);
	    printf("\n");
	    }
}

void sor(
  double omg,
  double dx,
  double dy,
  int il, int ir, int jt, int jb, int rank_l, int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv, MPI_Status *status, int chunk,
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res
) {
  int i,j;
  double rloc, residual;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));


  /* set boundary values */
  if(il==1)
  		for(j=jb-1; j<=jt+1; j++){
  		    P[0][j] = P[1][j];
  		}

  	if(ir==imax)
  		for(j=jb-1; j<=jt+1; j++){
  			P[imax+1][j] = P[imax][j];
  		}

  	if(jb==1)
  		for(i=il-1; i<=ir+1; i++){
  			P[i][0] = P[i][1];
  		}

  	if(jt==jmax)
  		for(i=il-1; i<=ir+1; i++){
  			P[i][jmax+1] = P[i][jmax];
  		}

  /* SOR iteration */
  for(i = il; i <= ir; i++)
    for(j = jb; j<=jt; j++)
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);


	/*if(rank_r==3){
		print_matr2("---------- P", -1, P, il, ir, jt, jb, 0, 0);
	}*/
  /* call the function pressure_comm(..) here, with the right parameters */
  pressure_comm(P, il, ir, jt, jb, rank_l, rank_r, rank_b, rank_t,  bufSend, bufRecv, status, chunk);
	/*if(rank_r==3){
		print_matr2("P", -1, P, il, ir, jt, jb, 0, 0);
	}*/

  /* compute the residual */
  rloc = 0;
  for(i = il; i <= ir; i++) {
    for(j = jb; j <= jt; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  rloc = rloc/(imax*jmax);


  /* compute and spread the sum of local residuals to all processors */
  MPI_Allreduce(&rloc, &residual, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  rloc = sqrt(residual);

  /* set residual */
  *res = rloc;
}

