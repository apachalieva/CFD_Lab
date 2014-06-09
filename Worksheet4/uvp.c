#include <stdlib.h>
#include "helper.h"
#include "uvp.h"
#include "parallel.h"


#define SQ(a) ((a)*(a))

/* central difference approximation of the second derivative in x */
double d2dx(double **m, int i, int j, double dx){
	return (m[i+1][j] - 2*m[i][j] + m[i-1][j]) / (SQ(dx));
}

/* central difference approximation of the second derivative in y */
double d2dy(double **m, int i, int j, double dy){
	return (m[i][j+1] - 2*m[i][j] + m[i][j-1]) / (SQ(dy));
}

/* forward difference approximation of the first derivative in x */
double ddx(double **m, int i, int j, double dx){
	return (m[i+1][j] - m[i][j]) / dx;
}

/* approximation of the first derivative of the square of u in x */
double du2dx(double **m, int i, int j, double dx, double alpha){
	return (
			SQ(m[i][j]+m[i+1][j]) - SQ(m[i-1][j]+m[i][j])
			+ alpha * ( abs(m[i][j]+m[i+1][j]) * (m[i][j]-m[i+1][j]) -  abs(m[i-1][j]+m[i][j]) * (m[i-1][j]-m[i][j]) )
	                       )/dx/4.0;
}

/* approximation of the first derivative of the square of v in y */
double dv2dy(double **m, int i, int j, double dy, double alpha){
	return (
			SQ(m[i][j]+m[i][j+1]) - SQ(m[i][j-1]+m[i][j])
			+ alpha * ( abs(m[i][j]+m[i][j+1]) * (m[i][j]-m[i][j+1]) -  abs(m[i][j-1]+m[i][j]) * (m[i][j-1]-m[i][j]))
	                       )/dy/4.0;
}


double duvdx(double ** u, double **v, int i, int j, double dx, double alpha){
	return (
			(u[i][j]+u[i][j+1]) * (v[i][j]+v[i+1][j]) - (u[i-1][j]+u[i-1][j+1]) * (v[i-1][j]+v[i][j])
			+ alpha * ( abs(u[i][j]+u[i][j+1]) * (v[i][j]-v[i+1][j]) -  abs(u[i-1][j]+u[i-1][j+1]) * (v[i-1][j]-v[i][j]) )
	                       )/dx/4.0;
}

double duvdy(double **u, double **v, int i, int j, double dy, double alpha){
	return (
			(u[i][j]+u[i][j+1]) * (v[i][j]+v[i+1][j]) - (u[i][j-1]+u[i][j]) * (v[i][j-1]+v[i+1][j-1])
			+ alpha * ( abs(v[i][j]+v[i+1][j]) * (u[i][j]-u[i][j+1]) -  abs(v[i][j-1]+v[i+1][j-1]) * (u[i][j-1]-u[i][j]) )
	                       )/dy/4.0;
}


void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int il,
  int ir,
  int jt,
  int jb,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
){

	int i,j;


	/*boundary*/
	if(il==1)
		for(j=jb-1; j<=jt+1; j++){
			F[0][j] = U[0][j];
		}

	if(ir==imax)
		for(j=jb-1; j<=jt+1; j++){
			F[imax][j] = U[imax][j];
		}

	if(jb==1)
		for(i=il-1; i<=ir+1; i++){
			G[i][0] = V[i][0];
		}

	if(jt==jmax)
		for(i=il-1; i<=ir+1; i++){
			G[i][jmax] = V[i][jmax];
		}


	/* boundary of computation: compute values on the local domain boundary just if it is not the global domain boundary */
	int ilb = il==1? il : il-1;
	int irb = ir==imax? ir-1 : ir;

	/* inner values */
	for(i=ilb; i<=irb; i++)
		for(j=jb; j<=jt; j++)
			F[i][j] = U[i][j] + dt * (
					(d2dx(U,i,j,dx) + d2dy(U,i,j,dy))/Re - du2dx(U, i, j, dx, alpha) - duvdy(U,V,i,j,dy, alpha) + GX
					);

	/* boundary of computation: compute values on the local domain boundary just if it is not the global domain boundary */
	int jbb = jb==1? jb : jb-1;
	int jtb = jt==jmax? jt-1 : jt;

	for(i=il; i<=ir; i++)
		for(j=jbb; j<=jtb; j++)
			G[i][j] = V[i][j] + dt * (
					(d2dx(V,i,j,dx) + d2dy(V,i,j,dy))/Re - duvdx(U, V, i, j, dx, alpha) - dv2dy(V,i,j,dy, alpha) + GY
					);

}


void calculate_rs(
  double dt,
  double dx,
  double dy,
  int il,
  int ir,
  int jt,
  int jb,
  double **F,
  double **G,
  double **RS
){
	int i,j;

	for(i=il; i<=ir; i++)
		for(j=jb; j<=jt; j++)
			RS[i][j] = 1/dt*((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy) ;

}




void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int il, int ir, int jt, int jb, double *bufSend, double *bufRecv,
  int imax,
  int jmax,
  double **U,
  double **V
){

	/* local U max */
	bufSend[0] = fmatrix_max(U, il-2, ir+1, jb-1, jt+1);
	/* local V max */
	bufSend[1] = fmatrix_max(V, il-1, ir+1, jb-2, jt+1);

	/* compute and spread the maximum value in U and V */
	MPI_Allreduce(bufSend, bufRecv, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	/* all processors compute the dt */
	*dt = tau * fmin( fmin(Re/2/(1/(dx*dx) + 1/(dy*dy)) , dx/bufRecv[0] 	), dy/bufRecv[1] );
}


void calculate_uv(
  double dt,
  double dx,
  double dy,
  int il, int ir, int jt, int jb, int rank_l, int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv, MPI_Status *status, int chunk,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
){
	int i,j;



	int ilb = il==1? il : il-1;
	int irb = ir==imax? ir-1 : ir;

	/* inner values */
	for(i=ilb; i<=irb; i++)
		for(j=jb; j<=jt; j++)
			U[i][j] = F[i][j] - dt/dx*(P[i+1][j]-P[i][j]);

	int jbb = jb==1? jb : jb-1;
	int jtb = jt==jmax? jt-1 : jt;

	/* inner values */
	for(i=il; i<=ir; i++)
		for(j=jbb; j<=jtb; j++)
			V[i][j] = G[i][j] - dt/dy*(P[i][j+1]-P[i][j]);

	/* comunicate u and v */
	uv_comm(U, V, il, ir, jt, jb, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, status, chunk);

}

