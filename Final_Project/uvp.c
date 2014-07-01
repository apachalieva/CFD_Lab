#include "helper.h"
#include <math.h>

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

/* forward difference approximation of the first derivative in y */
inline double ddy(double **m, int i, int j, double dy){
	return (m[i][j+1] - m[i][j]) / dy;
}

/* approximation of the first derivative of the square of u in x */
inline double du2dx(double **m, int i, int j, double dx, double alpha){
	return (
			SQ(m[i][j]+m[i+1][j]) - SQ(m[i-1][j]+m[i][j])
			+  alpha * ( fabs(m[i][j]+m[i+1][j]) * (m[i][j]-m[i+1][j]) -  fabs(m[i-1][j]+m[i][j]) * (m[i-1][j]-m[i][j]) )
            )/dx/4.0;
}

/* approximation of the first derivative of the square of v in y */
inline double dv2dy(double **m, int i, int j, double dy, double alpha){
	return (
			SQ(m[i][j]+m[i][j+1]) - SQ(m[i][j-1]+m[i][j])
			+ alpha * ( fabs(m[i][j]+m[i][j+1]) * (m[i][j]-m[i][j+1]) -  fabs(m[i][j-1]+m[i][j]) * (m[i][j-1]-m[i][j]))
	                       )/dy/4.0;
}


inline double duvdx(double ** u, double **v, int i, int j, double dx, double alpha){
	return (
			(u[i][j]+u[i][j+1]) * (v[i][j]+v[i+1][j]) - (u[i-1][j]+u[i-1][j+1]) * (v[i-1][j]+v[i][j])
			+ alpha * ( fabs(u[i][j]+u[i][j+1]) * (v[i][j]-v[i+1][j]) -  fabs(u[i-1][j]+u[i-1][j+1]) * (v[i-1][j]-v[i][j]) )
	                       )/dx/4.0;
}

inline double duvdy(double **u, double **v, int i, int j, double dy, double alpha){
	return (
			(u[i][j]+u[i][j+1]) * (v[i][j]+v[i+1][j]) - (u[i][j-1]+u[i][j]) * (v[i][j-1]+v[i+1][j-1])
			+ alpha * ( fabs(v[i][j]+v[i+1][j]) * (u[i][j]-u[i][j+1]) -  fabs(v[i][j-1]+v[i+1][j-1]) * (u[i][j-1]-u[i][j]) )
	                       )/dy/4.0;
}


inline double Rt(double **k, double **e, double nu, int i, int j){
	/*printf("KA-------- %d %d: %f %f %f % f\n", i,j,SQ(k[i][j]), nu , e[i][j], SQ(k[i][j]) / nu / e[i][j]);*/
	return SQ(k[i][j]) / nu / e[i][j];
}

inline double get_delta(int i, int j, int imax, int jmax, double dx, double dy){
	return MIN(
					MIN(  (i-0.5)*dx, (imax-(i-0.5))*dx  ),
					MIN(  (j-0.5)*dy,  (jmax-(j-0.5))*dy )
		);
}

double Rd(double **k, double nu, double delta, int i, int j){
	/*printf("** %f, %f, %f, %f \n", k[i][j],nu,delta,sqrt(fabs(k[i][j])) * delta / nu);*/
	return sqrt(fabs(k[i][j])) * delta / nu;
}

double fnu(double **k, double **e, double nu, double delta, int i, int j){
	/*return SQ( 1-exp( -0.0165 * Rd(k,nu,delta,i,j) ) )  * (1 + 20.5 / Rt(k,e,nu,i,j) );*/
	/*printf("	%f: %f\n", -0.0165 * Rd(k,nu,delta,i,j), exp( -0.0165 * Rd(k,nu,delta,i,j) ));*/

	return 1.0;
	/*printf("KA-----* %f %f %f % f\n", -0.0165 * Rd(k,nu,delta,i,j), exp( -0.0165 * Rd(k,nu,delta,i,j) ), Rt(k,e,nu,i,j), ( 1- exp( -0.0165 * Rd(k,nu,delta,i,j) ) )  *( 1-exp( -0.0165 * Rd(k,nu,delta,i,j) ) )  * (1 + 20.5 / Rt(k,e,nu,i,j) ));*/
	/*return ( 1- exp( -0.0165 * Rd(k,nu,delta,i,j) ) )  *( 1-exp( -0.0165 * Rd(k,nu,delta,i,j) ) )  * (1 + 20.5 / Rt(k,e,nu,i,j) );*/
}

inline double f1(double **k, double **e, double nu, double delta, int i, int j){
	return 1.0;
	/*return 1 + (0.05 / fnu(k,e,nu,delta,i,j) ) * (0.05 / fnu(k,e,nu,delta,i,j) ) * (0.05 / fnu(k,e,nu,delta,i,j) );*/
}

inline double f2(double **k, double **e, double nu, double delta, int i, int j){
	return 1.0;
	/*return 1 - exp( - SQ( Rt(k,e,nu,i,j) ) );*/
}

double visc_t(double **k, double **e, double nu, double cn, double delta, int i, int j){
	/*printf("KA---- %f %f %f %f \n", cn, fnu(k,e,nu,delta,i,j), SQ(k[i][j]) , e[i][j]);*/
	return cn * fnu(k,e,nu,delta,i,j) * SQ(k[i][j]) / e[i][j];
}

inline double visc(double **k, double **e, double nu, double cn, double delta, int i, int j){
	return nu + visc_t(k,e,nu,cn,delta,i,j);
}

inline double visc2(double **k, double **e, double nu, double cn, double delta, int i, int j){
	return nu + visc_t(k,e,nu,cn,delta,i,j)/1.3;
}

inline double visc_corner(double **k, double **e, double nu, double cn, double delta, int i, int j){
	return ( visc(k,e,nu,cn,delta, i+1, j+1) + visc(k,e,nu,cn,delta, i, j+1) + visc(k,e,nu,cn,delta, i+1, j) + visc(k,e,nu,cn,delta, i, j) ) / 4.0;
}

inline double dndudxdx(double **k, double **e, double nu, double cn, double delta, double **u, int i, int j, double dx){
	return (visc(k,e,nu,cn,delta, i+1, j) * (u[i+1][j] - u[i][j]) - visc(k,e,nu,cn,delta, i, j ) * (u[i][j] - u[i-1][j])) / SQ(dx);
}

inline double dndvdydy(double **k, double **e, double nu, double cn, double delta, double **v, int i, int j, double dy){
	return (visc(k,e,nu,cn,delta, i, j+1) * (v[i][j+1] - v[i][j]) - visc(k,e,nu,cn,delta, i, j ) * (v[i][j] - v[i][j-1])) / SQ(dy);
}


inline double dndudypdvdxdy(double **k, double **e, double nu, double cn, double delta, double **u, double **v, int i, int j, double dx, double dy){
	return (
			  visc_corner(k,e,nu,cn,delta, i,j  )  * ( (u[i][j+1] - u[i][j  ])/dy + (v[i+1][j  ] - v[i][j  ])/dx )
			- visc_corner(k,e,nu,cn,delta, i,j-1)  * ( (u[i][j  ] - u[i][j-1])/dy + (v[i+1][j-1] - v[i][j-1])/dx )
			) /  dy;
}

inline double dndudypdvdxdx(double **k, double **e, double nu, double cn, double delta, double **u, double **v, int i, int j, double dx, double dy){
	return (
			  visc_corner(k,e,nu,cn,delta, i   ,j)  * ( (u[i  ][j+1] - u[i  ][j  ])/dy + (v[i+1][j  ] - v[i  ][j  ])/dx )
			- visc_corner(k,e,nu,cn,delta, i-1,j )  * ( (u[i-1][j+1] - u[i-1][j  ])/dy + (v[i  ][j  ] - v[i-1][j  ])/dx )
			) /  dx;
}


inline double dvisct_dx( double **k, double **e, double nu, double cn, double delta, double dx, int i, int j ){
	double c1 = (visc( k, e, nu, cn, delta, i, j )+visc( k, e, nu, cn, delta, i+1, j ))/2;
	double c2 = (visc( k, e, nu, cn, delta, i-1, j )+visc( k, e, nu, cn, delta, i, j ))/2;
	double c3 = c1 * (k[i+1][j]-k[i][j]) - c2 *(k[i][j]-k[i-1][j]);
	c3=c3/SQ(dx);

	/*printf("KA---- %d %d: %f %f %f %f %f %f %f \n", i, j, c1, c2, c3, k[i+1][j], k[i][j], k[i][j], k[i-1][j] );*/

	return c3;
}

inline double dvisct_dy( double **k, double **e, double nu, double cn, double delta, double dy, int i, int j ){
	double c1 = (visc( k, e, nu, cn, delta, i, j )+visc( k, e, nu, cn, delta, i, j+1 ))/2;
	double c2 = (visc( k, e, nu, cn, delta, i, j-1 )+visc( k, e, nu, cn, delta, i, j ))/2;
	double c3 = c1 * (k[i][j+1]-k[i][j]) - c2 * (k[i][j]-k[i][j-1]);
	c3=c3/SQ(dy);

	/*printf("KA---- %f %f %f \n", c1, c2, c3 );*/




	return c3;
}

inline double dvisctEP_dx( double **k, double **e, double nu, double cn, double delta, double dx, int i, int j ){
	return ( ((visc_t( k, e, nu, cn, delta, i, j )+visc_t( k, e, nu, cn, delta, i+1, j ))/2)*(e[i+1][j]-e[i][j]) - 
	         ((visc_t( k, e, nu, cn, delta, i-1, j )+visc_t( k, e, nu, cn, delta, i, j ))/2)*(e[i][j]-e[i-1][j]) )/SQ(dx);
}

inline double dvisctEP_dy( double **k, double **e, double nu, double cn, double delta, double dy, int i, int j ){
	return ( ((visc_t( k, e, nu, cn, delta, i, j )+visc_t( k, e, nu, cn, delta, i, j+1 ))/2)*(e[i][j+1]-e[i][j]) -  ((visc_t( k, e, nu, cn, delta, i, j-1 )+visc_t( k, e, nu, cn, delta, i, j ))/2)*(e[i][j]-e[i][j-1]) )/SQ(dy);
}

double d_fnu_visctEP_dx( double **k, double **e, double nu, double cn, double delta, double dx, int i, int j ){

	double mean_f_1 = (fnu( k, e, nu, delta, i, j )+fnu( k, e, nu, delta, i+1, j ))/2;
	double mean_nu_1 = (visc2( k, e, nu, cn, delta, i, j )+visc2( k, e, nu, cn, delta, i+1, j ))/2;
	/*double mean_nu_1 = visc_t( k, e, nu, cn, delta, i, j);*/

	double c1 = mean_f_1 * mean_nu_1 * (e[i+1][j]-e[i][j]);

	double mean_f_2 = (fnu( k, e, nu, delta, i-1, j )+fnu( k, e, nu, delta, i, j ))/2;
	double mean_nu_2 = (visc2( k, e, nu, cn, delta, i-1, j )+visc2( k, e, nu, cn, delta, i, j ))/2;
	/*double mean_nu_2 = visc_t( k, e, nu, cn, delta, i, j );*/

	double c2 = mean_f_2*mean_nu_2*(e[i][j]-e[i-1][j]);

	double res = (c1-c2)/SQ(dx);


	return res;
}

double d_fnu_visctEP_dy( double **k, double **e, double nu, double cn, double delta, double dy, int i, int j ){
	/*return (  ((fnu( k, e, nu, delta, i, j )+fnu( k, e, nu, delta, i, j+1 ))/2)*((visc_t( k, e, nu, cn, delta, i, j )+visc_t( k, e, nu, cn, delta, i, j+1 ))/2)*(e[i][j+1]-e[i][j]) - ((fnu( k, e, nu, delta, i, j-1 )+fnu( k, e, nu, delta, i, j ))/2)*((visc_t( k, e, nu, cn, delta, i, j-1 )+visc_t( k, e, nu, cn, delta, i, j ))/2)*(e[i][j]-e[i][j-1]) )/SQ(dy); */

	double c1 = ((fnu( k, e, nu, delta, i, j )+fnu( k, e, nu, delta, i, j+1 ))/2);
	double c2 = ((visc2( k, e, nu, cn, delta, i, j )+visc2( k, e, nu, cn, delta, i, j+1 ))/2);

	double c3 = ((fnu( k, e, nu, delta, i, j-1 )+fnu( k, e, nu, delta, i, j ))/2);
	double c4 = ((visc2( k, e, nu, cn, delta, i, j-1 )+visc2( k, e, nu, cn, delta, i, j ))/2);

	double res = (  c1*c2*(e[i][j+1]-e[i][j]) - c3 *c4 *(e[i][j]-e[i][j-1]) )/(dy*dy);

	/*printf("d_fnu_visctEP_dy ---- %d %d: %f %f %f %f %f \n", i, j, c1, c2, c3, c4, res);*/

	return res;
}

/* backward discretization */
inline double dUdx( double **u, double dx, int i, int j ){
	return ( u[i][j] - u[i-1][j] ) / dx;
}

inline double dUdy( double **u, double dy, int i, int j ){
	/*printf("KA dudy ------ %d %d: %f %f %f %f %f\n", i, j, u[i][j], u[i][j+1], u[i-1][j+1], u[i][j-1], u[i-1][j-1]);*/
	return ( (u[i][j+1]+u[i-1][j+1]) - (u[i][j-1]+u[i-1][j-1]) ) / (4*dy);
}

inline double dVdx( double **v, double dx, int i, int j ){
	return ( (v[i+1][j]+v[i+1][j-1]) - (v[i-1][j]+v[i-1][j-1]) ) / (4*dx);
}

inline double dVdy( double **v, double dy, int i, int j ){
	return (v[i][j]-v[i][j-1]) / dy;
}

inline double gradU( double **u, double **v, double dx, double dy, int i, int j ){
	double c1 = 4*( dUdx( u, dx, i, j) )*( dUdx( u, dx, i, j) );
	double c2 = 2*( dUdy( u, dy, i, j) + dVdx( v, dx, i, j ) )*( dUdy( u, dy, i, j) +dVdx( v, dx, i, j ) );
	double c3 = 4*( dVdy( v, dy, i, j ) )*( dVdy( v, dy, i, j ) );

	/*printf("KA grad ---- %d %d: %f %f %f %f %f %f %f \n", i, j, c1, c2, c3, dUdx( u, dx, i, j), dUdy( u, dy, i, j), dVdx( v, dx, i, j ) , dVdy( v, dy, i, j ) );*/

	return c1 + c2 + c3;
}
/* alpha [0, 1] determines a weighted average of discretizing with central differences and the donor-cell discretization */
inline double dUkedx( double **u, double **ke, double dx, double alpha, int i, int j ){

	double c1 = u[i][j]*(ke[i][j]+ke[i+1][j])/2;
	double c2 = u[i-1][j]*(ke[i-1][j]+ke[i][j])/2;

	double c3 = fabs(u[i][j])*(ke[i][j]-ke[i+1][j])/2;
	double c4 = fabs(u[i-1][j])*(ke[i-1][j]-ke[i][j])/2;

	double res = (c1-c2)/dx +alpha/dx*(c3-c4);

	return res;
}

/* alpha [0, 1] determines a weighted average of discretizing with central differences and the donor-cell discretization */
inline double dVkedy( double **v, double **ke, double dy, double alpha, int i, int j ){
	double c1 = v[i][j]*(ke[i][j]+ke[i][j+1])/2;
	double c2 = v[i][j-1]*(ke[i][j-1]+ke[i][j])/2;
	double c3 = fabs(v[i][j])*(ke[i][j]-ke[i][j+1])/2;
	double c4 = fabs(v[i][j-1])*(ke[i][j-1]-ke[i][j])/2;
	double res = (c1-c2)/dy+alpha/dy*(c3-c4);


	return res;
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
  double **G,
  double **K,
  double **E,
  double nu,
  double cn,
  int **Flag
){
	int i,j;
	double delta;
	
	/* boundary conditions */
	for(j=1; j<=jmax; j++){
		F[ 0 ][ j ] = U[ 0 ][ j ];
		F[ imax ][ j ] = U[ imax ][ j ];
	}

	for(i=1; i<=imax; i++){
		G[i][0] = V[i][0];
		G[i][jmax] = V[i][jmax];
	}
	
	/* Boundary conditions for the obstacle cells */
	for( i = 1; i <= imax; i++ )
	    for( j = 1; j <= jmax; j++ )
		if( IS_BOUNDARY(Flag[i][j])  ){
		      /* Boundary conditions for obstacles with North-Eastern fluid cell */
		      if( ( Flag[ i ][ j ] & B_NE ) == B_NE ){
			  F[ i ][ j ] = U[ i ][ j ]; 
			  G[ i ][ j ] = V[ i ][ j ];
		      } else
		      
		      /* Boundary conditions for obstacles with North-Western fluid cell */
		      if( ( Flag[ i ][ j ] & B_NW ) == B_NW ){
			  F[ i-1 ][ j ] = U[ i-1 ][ j ]; 
			  G[ i ][ j ] = V[ i ][ j ];
		      } else
		      
		      /* Boundary conditions for obstacles with South-Eastern fluid cell */
		      if( ( Flag[ i ][ j ] & B_SE ) == B_SE ){
			  F[ i ][ j ] = U[ i ][ j ]; 
			  G[ i ][ j-1 ] = V[ i ][ j-1 ];
		      } else
		      
		      /* Boundary conditions for obstacles with South-Western fluid cell */
		      if( ( Flag[ i ][ j ] & B_SW ) == B_SW ){
			  F[ i-1 ][ j ] = U[ i-1 ][ j ]; 
			  G[ i ][ j-1 ] = V[ i ][ j-1 ];
		      }else
		      /* Boundary conditions for obstacles with Northern fluid cell */
		      if( ( Flag[ i ][ j ] & B_N ) == B_N ){
			  G[ i ][ j ] = V[ i ][ j ];
		      } else
		      /* Boundary conditions for obstacles with Southern fluid cell */
		      if( ( Flag[ i ][ j ] & B_S ) == B_S ) {
			  G[ i ][ j-1 ] = V[ i ][ j-1 ];
		      } else
		      /* Boundary conditions for obstacles with Western fluid cell */
		      if( ( Flag[ i ][ j ] & B_W ) == B_W ){
			  F[ i-1 ][ j ] = U[ i-1 ][ j ];
		      } else
		      /* Boundary conditions for obstacles with Eastern fluid cell */
		      if( ( Flag[ i ][ j ] & B_E ) == B_E ){
			  F[ i ][ j ] = U[ i ][ j ];
		      }
		}

	/* inner values */
	for(i=1; i<=imax-1; i++)
		for(j=1; j<=jmax; j++)
			if( IS_FLUID(Flag[i][j]) && IS_FLUID(Flag[i+1][j]) ){
				delta = get_delta(i, j, imax, jmax, dx, dy);

				F[i][j] = U[i][j] + dt * (
						2.0 * dndudxdx(K, E, nu, cn, delta, U, i, j, dx)
						+ dndudypdvdxdy(K, E, nu, cn, delta, U, V, i, j, dx, dy)
						/*- 2.0 * ddx(K,i,j,dx) / 3.0*/
						- du2dx(U, i, j, dx, alpha)
						- duvdy(U,V,i,j,dy, alpha)
						+ GX
						);
			}

	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax-1; j++)
			if( IS_FLUID(Flag[i][j]) && IS_FLUID(Flag[i][j+1]) ){
				delta = get_delta(i, j, imax, jmax, dx, dy);

				G[i][j] = V[i][j] + dt * (
						dndudypdvdxdx(K, E, nu, cn, delta, U, V, i, j, dx, dy)
						+ 2.0 * dndvdydy(K, E, nu, cn, delta, V, i, j, dy)
						/*- 2.0 * ddy(K,i,j,dy) / 3.0*/
						- duvdx(U,V,i,j,dx, alpha)
						- dv2dy(V, i, j, dy, alpha)
						+ GY
						);
			}

}


void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS,
  int **Flag
){
	int i,j;

	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
				RS[i][j] = IS_FLUID(Flag[i][j]) * 1/dt*((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy) ;

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
  double **P,
  int **Flag
){
	int i,j;

	for(i=1; i<=imax-1; i++)
		for(j=1; j<=jmax; j++)
			if(  IS_FLUID(Flag[i][j]) && IS_FLUID(Flag[i+1][j]) )
				U[i][j] = F[i][j] - dt/dx*(P[i+1][j]-P[i][j]);

	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax-1; j++)
			if( IS_FLUID(Flag[i][j]) && IS_FLUID(Flag[i][j+1]) )
				V[i][j] = G[i][j] - dt/dy*(P[i][j+1]-P[i][j]);


}


void comp_KAEP( 
  double Re, 
  double nu,
  double cn, 
  double ce, 
  double c1, 
  double c2,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U, 
  double **V, 
  double **KA, 
  double **EP, 
  double GX, 
  double GY, 
  int **Flag 
){
	int i, j;
	double delta;
	
	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
			if(IS_FLUID(Flag[i][j])	) {
				delta = get_delta(i, j, imax, jmax, dx, dy);

				/*printf("KA-- %d %d: %f %f %f %f %f %f %f %f %f %f \n", i,j,KA[i][j], dt,	dvisct_dx( KA, EP, nu, cn, delta, dx, i, j )
									, dvisct_dy( KA, EP, nu, cn, delta, dy, i, j )
								    ,dUkedx( U, KA, dx, alpha, i, j )
								    , dVkedy( V, KA, dy, alpha, i, j )
							        , visc_t( KA, EP, nu, cn, delta, i, j )
							         ,gradU( U, V, dx, dy, i, j )
							         , EP[i][j]
							         , KA[i][j] + dt*( dvisct_dx( KA, EP, nu, cn, delta, dx, i, j ) + dvisct_dy( KA, EP, nu, cn, delta, dy, i, j )- dUkedx( U, KA, dx, alpha, i, j )- dVkedy( V, KA, dy, alpha, i, j ) + 0.5 * visc_t( KA, EP, nu, cn, delta, i, j )* gradU( U, V, dx, dy, i, j )- EP[i][j]));
							         */


				KA[i][j] = KA[i][j] + dt*(
						dvisct_dx( KA, EP, nu, cn, delta, dx, i, j )
						+ dvisct_dy( KA, EP, nu, cn, delta, dy, i, j )
					    - dUkedx( U, KA, dx, alpha, i, j )
					    - dVkedy( V, KA, dy, alpha, i, j )
				        + 0.5 * visc_t( KA, EP, nu, cn, delta, i, j )
				              * gradU( U, V, dx, dy, i, j )
				        - EP[i][j]);

				/*printf("EP-- %d %d --> %f %f %f %f %f %f %f %f %f %f \n",i,j,EP[i][j] , dt
						  ,d_fnu_visctEP_dx( KA, EP, nu, cn, delta, dx, i, j )
		                , d_fnu_visctEP_dy( KA, EP, nu, cn, delta, dy, i, j )
		                , dUkedx( U, EP, dx, alpha, i, j )
		                , dVkedy( V, EP, dy, alpha, i, j )
		                , ((c1*f1(KA, EP, nu, delta, i, j))/2),
		                KA[i][j],
		                gradU(U, V, dx, dy, i, j)
		                , c2*f2( KA, EP, nu, delta, i, j )
				);*/
			}


	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
			if(IS_FLUID(Flag[i][j])	) {
				delta = get_delta(i, j, imax, jmax, dx, dy);
				EP[i][j] = EP[i][j] + dt*(
													  /*(ce/cn)**/ d_fnu_visctEP_dx( KA, EP, nu, cn, delta, dx, i, j )
									                + /*(ce/cn)**/ d_fnu_visctEP_dy( KA, EP, nu, cn, delta, dy, i, j )
									                - dUkedx( U, EP, dx, alpha, i, j )
									                - dVkedy( V, EP, dy, alpha, i, j )
									                + c1*f1(KA, EP, nu, delta, i, j)/2
									                	*KA[i][j]
									                	*gradU(U, V, dx, dy, i, j)
									                - c2*f2( KA, EP, nu, delta, i, j )
														*SQ(EP[i][j])
														/KA[i][j]
											);
			}
}

