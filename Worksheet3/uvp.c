/*
 * uvp.c
 *
 *  Created on: Apr 16, 2014
 *      Author: davide
 */

#include <stdlib.h>
#include "helper.h"


#define SQ(a) ((a)*(a))

/* central difference approximation of the second derivative in x */
inline double d2dx(double **m, int i, int j, double dx){
	return (m[i+1][j] - 2*m[i][j] + m[i-1][j]) / (SQ(dx));
}

/* central difference approximation of the second derivative in y */
inline double d2dy(double **m, int i, int j, double dy){
	return (m[i][j+1] - 2*m[i][j] + m[i][j-1]) / (SQ(dy));
}

/* forward difference approximation of the first derivative in x */
inline double ddx(double **m, int i, int j, double dx){
	return (m[i+1][j] - m[i][j]) / dx;
}

/* approximation of the first derivative of the square of u in x */
inline double du2dx(double **m, int i, int j, double dx, double alpha){
	return (
			SQ(m[i][j]+m[i+1][j]) - SQ(m[i-1][j]+m[i][j])
			+ alpha * ( abs(m[i][j]+m[i+1][j]) * (m[i][j]-m[i+1][j]) -  abs(m[i-1][j]+m[i][j]) * (m[i-1][j]-m[i][j]) )
	                       )/dx/4.0;
}

/* approximation of the first derivative of the square of v in y */
inline double dv2dy(double **m, int i, int j, double dy, double alpha){
	return (
			SQ(m[i][j]+m[i][j+1]) - SQ(m[i][j-1]+m[i][j])
			+ alpha * ( abs(m[i][j]+m[i][j+1]) * (m[i][j]-m[i][j+1]) -  abs(m[i][j-1]+m[i][j]) * (m[i][j-1]-m[i][j]))
	                       )/dy/4.0;
}


inline double duvdx(double ** u, double **v, int i, int j, double dx, double alpha){
	return (
			(u[i][j]+u[i][j+1]) * (v[i][j]+v[i+1][j]) - (u[i-1][j]+u[i-1][j+1]) * (v[i-1][j]+v[i][j])
			+ alpha * ( abs(u[i][j]+u[i][j+1]) * (v[i][j]-v[i+1][j]) -  abs(u[i-1][j]+u[i-1][j+1]) * (v[i-1][j]-v[i][j]) )
	                       )/dx/4.0;
}

inline double duvdy(double **u, double **v, int i, int j, double dy, double alpha){
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
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
){

	int i,j;

	/*boundary*/
	for(j=1; j<=jmax; j++){
		F[0][j] = U[0][j];
		F[imax][j] = U[imax][j];
	}

	for(i=1; i<=imax; i++){
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}

	/* inner values */
	for(i=1; i<=imax-1; i++)
		for(j=1; j<=jmax; j++)
			F[i][j] = U[i][j] + dt * (
					(d2dx(U,i,j,dx) + d2dy(U,i,j,dy))/Re - du2dx(U, i, j, dx, alpha) - duvdy(U,V,i,j,dy, alpha) + GX
					);

	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax-1; j++)
			G[i][j] = V[i][j] + dt * (
					(d2dx(V,i,j,dx) + d2dy(V,i,j,dy))/Re - duvdx(U, V, i, j, dx, alpha) - dv2dy(V,i,j,dy, alpha) + GY
					);

}


void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
){
	int i,j;

	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
			RS[i][j] = 1/dt*((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy) ;

}




void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
){


	*dt = tau * fmin(
					fmin(Re/2/(1/(dx*dx) + 1/(dy*dy)) ,
							dx/fmatrix_max(U,0,imax+1,0,jmax+1) 	),
								dy/fmatrix_max(V,0,imax+1,0,jmax+1)
							);

}


void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
){
	int i,j;

	for(i=1; i<=imax-1; i++)
		for(j=1; j<=jmax; j++)
			U[i][j] = F[i][j] - dt/dx*(P[i+1][j]-P[i][j]);

	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax-1; j++)
			V[i][j] = G[i][j] - dt/dy*(P[i][j+1]-P[i][j]);


}

