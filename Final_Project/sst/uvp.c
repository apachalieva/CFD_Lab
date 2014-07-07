#include "uvp.h"
#include "helper.h"
#include <math.h>


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
	double c1 = 4.0*( dUdx( u, dx, i, j) )*( dUdx( u, dx, i, j) );
	double c2 = 2.0*( dUdy( u, dy, i, j) + dVdx( v, dx, i, j ) )*( dUdy( u, dy, i, j) +dVdx( v, dx, i, j ) );
	double c3 = 4.0*( dVdy( v, dy, i, j ) )*( dVdy( v, dy, i, j ) );

	/*printf("KA grad ---- %d %d: %f %f %f %f %f %f %f \n", i, j, c1, c2, c3, dUdx( u, dx, i, j), dUdy( u, dy, i, j), dVdx( v, dx, i, j ) , dVdy( v, dy, i, j ) );*/

	return c1 + c2 + c3;
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

inline double ddx_fw(double** M, double dx, int i, int j){
	return (M[i+1][j]-M[i][j])/dx;
}

inline double ddx_bw(double** M, double dx, int i, int j){
	return (M[i][j]-M[i-1][j])/dx;
}

inline double ddy_fw(double** M, double dy, int i, int j){
	return (M[i][j+1]-M[i][j])/dy;
}

inline double ddy_bw(double** M, double dy, int i, int j){
	return (M[i][j]-M[i][j-1])/dy;
}

inline double dUdx_hedge( double **u, double dx, int i, int j ){
	return ( (u[i][j+1]+u[i][j]) - (u[i-1][j+1]+u[i-1][j]) ) / (2*dx);
}
inline double dUdy_hedge( double **u, double dy, int i, int j ){
	return ( (u[i][j+1]+u[i-1][j+1]) - (u[i][j]+u[i-1][j]) ) / (2*dy);
}
inline double dUdy_vedge( double **u, double dy, int i, int j ){
	return ( u[i][j+1] - u[i][j-1] ) / (2*dy);
}
inline double dVdx_hedge( double **v, double dx, int i, int j ){
	return  (v[i+1][j]-v[i-1][j] ) / (2*dx);
}
inline double dVdx_vedge( double **v, double dx, int i, int j ){
	return ( (v[i][j]+v[i][j-1]) - (v[i-1][j]+v[i-1][j-1])  ) / (2*dx);
}
inline double dVdy_vedge( double **v, double dy, int i, int j ){
	return ( (v[i-1][j]+v[i][j]) - (v[i-1][j-1]+v[i][j-1]) ) / (2*dy);
}

inline double dkdw(double** K, double** W, int i, int j, double dx, double dy){
	double dkdx, dkdy, dwdx, dwdy, app1;
	if(i!=0){
		dkdx=ddx_bw(K,dx,i,j);
		dwdx=ddx_bw(W,dx,i,j);
	}
	else{
		dkdx=ddx_fw(K,dx,i,j);
		dwdx=ddx_fw(W,dx,i,j);
	}
	if(j!=0){
		dkdy=ddy_bw(K,dy,i,j);
		dwdy=ddy_bw(W,dy,i,j);
	}
	else{
		dkdy=ddy_fw(K,dy,i,j);
		dwdy=ddy_fw(W,dy,i,j);
	}
	app1 = dkdx*dwdx+dkdy*dwdy;
	return app1;
}

inline double CDkw(double** K, double** W, int i, int j,double dx, double dy){
	double app= 2*sigma_w2/W[i][j]*dkdw(K,W,i,j, dx, dy);
	double app2 =fmax(app, pow(10,-10));
	return app2;
}

inline double f1(double** K, double** W, double** delta, double nu,int i, int j,double dx, double dy){ 	/* putting rho=1 */
	double app1 = CDkw(K,W,i,j,dx,dy);
	double app2 = sqrt(K[i][j])/(beta_star*W[i][j]*delta[i][j]);
	double app3 = 500.0*nu/(SQ(delta[i][j])*W[i][j]);
	double app4 = 4*sigma_w2*K[i][j]/(app1*SQ(delta[i][j]) );
	double app5 = fmin( fmax( app2,  app3) , app4 );
	double app6 = tanhf( pow(app5, 4)  );

	return app6;
}

inline double f2(double** K, double** W, double** delta, double nu, int i, int j,double dx, double dy){
	double app1 =  2.0*sqrt(K[i][j])/(beta_star*W[i][j]*delta[i][j]);
	double app2 = 500.0*nu/(SQ(delta[i][j])*W[i][j] ) ;
	double app3 = fmax(app1, app2);
	double app4 = tanhf( app3*app3 );

	return app4;
}



inline double S_funct(int i, int j, double **U, double **V, double dx, double dy){
	double dudy, dvdx;
	if (i!=0){
		dvdx=(V[i][j]-V[i-1][j])/dx;
	}
	else{
		dvdx=(V[i+1][j]-V[i][j])/dx;
	}
	if (j!=0){
		dudy=(U[i][j]-U[i][j-1])/dy;
	}
	else{
		dudy=(U[i][j+1]-U[i][j])/dy;
	}

	return fabs(0.5*(dudy-dvdx));
}

/*
 * for S:
 * 0.5*fabs(dudy(lv, parameters) - dvdx(lv, parameters))
 */

inline double visc_t(double **K, double **W, double** U, double** V, double nu, double** delta, int i, int j, double dx, double dy){
	/*printf("KA---- %f %f %f %f \n", cn, fnu(k,e,nu,delta,i,j), SQ(k[i][j]) , e[i][j]);*/
	/*
	if(i!=0 && j!=0 && i!=i)
		s=sqrt(0.5*gradU( U, V, dx, dy, i, j ));
	*/
	double app1 = a1 * K[i][j] / fmax( a1*W[i][j], S_funct( i, j, U, V, dx, dy) *f2( K, W, delta, nu, i, j, dx, dy));
	/*
	printf("visc_t --> [%d, %d] %f %f %f %f %f",i,j,a1*K[i][j],a1*W[i][j],S_funct( i, j, U, V, dx, dy),f2( K, W, delta, nu, i, j, dx, dy), fmax( a1*W[i][j], S_funct( i, j, U, V, dx, dy) *f2( K, W, delta, nu, i, j, dx, dy)) );
	*/
	return app1;
}

inline double sigma_k(double** K, double** W, double** delta, double nu,int i, int j,double dx, double dy){
	double app1 = f1( K, W, delta, nu,i,j, dx, dy);
	double app2 = sigma_k1*app1 +sigma_k2*(1-app1);
	return app2;
}

inline double sigma_w(double** K, double** W, double** delta, double nu,int i, int j,double dx, double dy){
	double app1 = f1( K, W, delta, nu,i,j, dx, dy);
	double app2 = sigma_w1*app1+sigma_w2*(1-app1);
	return app2;
}

inline double visc(double **k, double **w, double** U, double **V, double nu, double** delta, int i, int j, double dx, double dy){
	double app1 = nu + sigma_k(k,w,delta,nu,i,j,dx,dy)*visc_t(k,w,U,V,nu,delta, i, j, dx, dy);
	return app1;
}

inline double visc2(double **k, double **w, double** U, double** V, double nu, double** delta, int i, int j, double dx, double dy){
	double app1 =nu + sigma_w(k,w,delta,nu,i,j,dx,dy)*visc_t(k,w,U,V,nu,delta, i, j, dx, dy);
	return app1;
}

inline double visc3(double **k, double **w, double** U, double** V, double nu, double** delta, int i, int j, double dx, double dy){
	double app1 =nu + visc_t(k,w,U,V,nu,delta, i, j, dx, dy);
	return app1;
}

inline double visc_corner(double **k, double **e,double** U, double** V, double nu, double** delta, int i, int j, double dx, double dy){
	double app1=( visc3(k,e,U,V,nu,delta, i+1, j+1, dx, dy) + visc3(k,e,U,V,nu,delta, i, j+1, dx, dy) + visc3(k,e,U,V,nu,delta, i+1, j, dx, dy) + visc3(k,e,U,V,nu,delta, i, j, dx, dy) ) / 4.0;

	return app1;
}

inline double dndudxdx(double **k, double **e, double nu, double** delta, double **u, double** V, int i, int j, double dx, double dy){
	double app1 = (visc3(k,e,u, V, nu,delta, i+1, j, dx, dy) * (u[i+1][j] - u[i][j]) - visc3(k,e,u, V, nu,delta, i, j, dx, dy) * (u[i][j] - u[i-1][j])) / SQ(dx);
	return app1;
}

inline double dndvdydy(double **k, double **e, double nu, double** delta, double** U, double **v, int i, int j, double dx, double dy){
	double app1 = (visc3( k, e, U, v, nu, delta, i, j+1, dx, dy) * (v[i][j+1] - v[i][j]) - visc3(k,e,U, v, nu,delta, i, j, dx, dy) * (v[i][j] - v[i][j-1])) / SQ(dy);
	return app1;
}


inline double dndudypdvdxdy(double **k, double **e, double nu, double** delta, double **u, double **v, int i, int j, double dx, double dy){
	return (
			  visc_corner(k,e,u,v,nu,delta, i,j , dx , dy )  * ( (u[i][j+1] - u[i][j  ])/dy + (v[i+1][j  ] - v[i][j  ])/dx )
			- visc_corner(k,e,u,v,nu,delta, i,j-1, dx , dy)  * ( (u[i][j  ] - u[i][j-1])/dy + (v[i+1][j-1] - v[i][j-1])/dx )
			) /  dy;
}

inline double dndudypdvdxdx(double **k, double **e, double nu, double** delta, double **u, double **v, int i, int j, double dx, double dy){
	return (
			  visc_corner(k,e,u,v,nu,delta, i   ,j, dx , dy)  * ( (u[i  ][j+1] - u[i  ][j  ])/dy + (v[i+1][j  ] - v[i  ][j  ])/dx )
			- visc_corner(k,e,u,v,nu,delta, i-1,j, dx , dy)  * ( (u[i-1][j+1] - u[i-1][j  ])/dy + (v[i  ][j  ] - v[i-1][j  ])/dx )
			) /  dx;
}


inline double dvisct_dx( double **k, double **e, double** U, double** V, double nu, double** delta, double dx, double dy, int i, int j ,int choice){
	double app1,app2,app3;
	if (choice==1){
		app1 = (visc( k, e, U, V, nu, delta, i, j , dx, dy)+visc( k, e, U, V, nu, delta, i+1, j, dx, dy))/2;
		app2 = (visc( k, e, U, V, nu, delta, i-1, j, dx, dy)+visc( k, e, U, V, nu, delta, i, j, dx, dy))/2;
		app3 = app1 * (k[i+1][j]-k[i][j]) - app2 *(k[i][j]-k[i-1][j]);
	}
	else{
		app1 = (visc2( k, e, U, V, nu, delta, i, j , dx, dy)+visc2( k, e, U, V, nu, delta, i+1, j, dx, dy))/2;
		app2 = (visc2( k, e, U, V, nu, delta, i-1, j, dx, dy)+visc2( k, e, U, V, nu, delta, i, j, dx, dy))/2;
		app3 = app1 * (e[i+1][j]-e[i][j]) - app2 *(e[i][j]-e[i-1][j]);
	}
	app3=app3/SQ(dx);

	return app3;
}

inline double dvisct_dy( double **k, double **e, double** U, double** V, double nu, double** delta, double dx, double dy, int i, int j, int choice){
	double app1, app2, app3;
	if (choice==1){
		app1 = (visc( k, e, U, V, nu, delta, i, j, dx, dy)+visc( k, e, U, V, nu, delta, i, j+1, dx, dy))/2;
		app2 = (visc( k, e, U, V, nu, delta, i, j-1, dx, dy)+visc( k, e, U, V, nu, delta, i, j, dx, dy))/2;
		app3 = app1 * (k[i][j+1]-k[i][j]) - app2 * (k[i][j]-k[i][j-1]);
	}
	else {
		app1 = (visc2( k, e, U, V, nu, delta, i, j, dx, dy)+visc2( k, e, U, V, nu, delta, i, j+1, dx, dy))/2;
		app2 = (visc2( k, e, U, V, nu, delta, i, j-1, dx, dy)+visc2( k, e, U, V, nu, delta, i, j, dx, dy))/2;
		app3 = app1 * (e[i][j+1]-e[i][j]) - app2 * (e[i][j]-e[i][j-1]);
	}
	app3=app3/SQ(dy);

	return app3;
}




inline double Pk(double** K, double** W, double** U, double** V, double** delta, double nu, double dx, double dy, int i, int j){
	double app1 = 0.5*visc_t( K, W,U,V, nu, delta, i, j, dx, dy) * gradU( U, V, dx, dy, i, j);
	return app1;
}

inline double P_tilda_k(double** K, double** W, double** U, double** V, double** delta, double nu, double dx, double dy, int i, int j){
	double app1 = Pk(K,W,U,V,delta,nu,dx,dy,i,j);
	double app2 = 10.0*beta_star*K[i][j]*W[i][j];
	double app3 = fmin( app1 , app2 );
	return app3;
}

/* alpha [0, 1] determines a weighted average of discretizing with central differences and the donor-cell discretization */
inline double dUkedx( double **u, double **ke, double dx, double alpha, int i, int j ){

	double app1 = u[i][j]*(ke[i][j]+ke[i+1][j])/2;
	double app2 = u[i-1][j]*(ke[i-1][j]+ke[i][j])/2;

	double app3 = fabs(u[i][j])*(ke[i][j]-ke[i+1][j])/2;
	double app4 = fabs(u[i-1][j])*(ke[i-1][j]-ke[i][j])/2;

	double app5 = (app1 - app2)/dx +alpha/dx*(app3 - app4);

	return app5;
}

inline double dVkedy( double **v, double **ke, double dy, double alpha, int i, int j ){
	double app1 = v[i][j]*(ke[i][j]+ke[i][j+1])/2;
	double app2 = v[i][j-1]*(ke[i][j-1]+ke[i][j])/2;
	double app3 = fabs(v[i][j])*(ke[i][j]-ke[i][j+1])/2;
	double app4 = fabs(v[i][j-1])*(ke[i][j-1]-ke[i][j])/2;
	double app5 = (app1-app2)/dy+alpha/dy*(app3-app4);

	return app5;
}

/* alpha [0, 1] determines a weighted average of discretizing with central differences and the donor-cell discretization */


inline double Alpha(double** K, double** W, double** delta, double nu,int i, int j,double dx, double dy){
	double app1=f1( K, W, delta, nu, i, j, dx, dy);
	double app2 = alpha_1*app1+(1.0-app1)*alpha_2;
	return app2;
}

inline double Beta(double** K, double** W, double** delta, double nu,int i, int j,double dx, double dy){
	double app1=f1( K, W, delta, nu, i, j, dx, dy);
	double app2=beta_1*app1+(1.0-app1)*beta_2;
	return app2;
}

/**
 * Calculation of F and G
 */
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
  double** delta,
  int **Flag
){
	int i,j;
	
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
			if( IS_FLUID(Flag[i][j]) && IS_FLUID(Flag[i+1][j]) )
				F[i][j] = U[i][j] + dt * (
						2.0 * dndudxdx(K, E, nu, delta, U, V, i, j, dx, dy)
						+ dndudypdvdxdy(K, E, nu, delta, U, V, i, j, dx, dy)
						- du2dx(U, i, j, dx, alpha)
						- duvdy(U,V,i,j,dy, alpha)
						+ GX
						);


	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax-1; j++)
			if( IS_FLUID(Flag[i][j]) && IS_FLUID(Flag[i][j+1]) )
				G[i][j] = V[i][j] + dt * (
						dndudypdvdxdx(K, E, nu, delta, U, V, i, j, dx, dy)
						+ 2.0 * dndvdydy(K, E, nu, delta, U, V, i, j, dx, dy)
						- duvdx(U,V,i,j,dx, alpha)
						- dv2dy(V, i, j, dy, alpha)
						+ GY
						);

}

/**
 * Calculation of the right-hand side of the pressure equation
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int    imax,
  int    jmax,
  double **F,
  double **G,
  double **RS,
  int    **Flag
){
	int i,j;

	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
				RS[i][j] = IS_FLUID(Flag[i][j]) * 1/dt*((F[i][j]-F[i-1][j])/dx + (G[i][j]-G[i][j-1])/dy) ;

}

/**
 * Calculation of the u and v values
 */
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


void comp_KAW(
  double Re, 
  double nu,
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
  double **W,
  double GX, 
  double GY, 
  double** delta,
  int **Flag 
){
	int i, j;
	
	/* Calculation of the k - turbulent kinetic energy */
	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
			if(IS_FLUID(Flag[i][j])	)
				KA[i][j] = KA[i][j] + dt*(
						- dUkedx( U, KA, dx, alpha, i, j )
						- dVkedy( V, KA, dy, alpha, i, j )
						+ dvisct_dx( KA, W, U, V, nu, delta, dx, dy, i, j, 1)
						+ dvisct_dy( KA, W, U, V, nu, delta, dx, dy, i, j, 1)
						- beta_star*KA[i][j]*W[i][j]
					    + P_tilda_k( KA, W, U, V, delta, nu, dx, dy, i, j)
				       );

	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
			if(IS_FLUID(Flag[i][j])	) 
				W[i][j] = W[i][j] + dt*(
													- dUkedx( U, W, dx, alpha, i, j)
									                - dVkedy( V, W, dy, alpha, i, j)
									                + Alpha( KA, W, delta, nu, i, j, dx, dy)
									                * (P_tilda_k( KA, W, U, V, delta, nu, dx, dy, i, j)/ visc_t( KA, W, U, V, nu, delta, i, j, dx, dy))
									                - Beta( KA, W, delta, nu, i, j, dx, dy) * SQ(W[i][j])
									                + dvisct_dx( KA, W, U, V, nu, delta, dx, dy, i, j, 0)
													+ dvisct_dy( KA, W, U, V, nu, delta, dx, dy, i, j, 0)
													+ 2.0*(1.0-f1( KA, W, delta, nu, i, j, dx, dy) ) * sigma_w2 / W[i][j]
													* dkdw( KA, W, i, j, dx, dy)
											);
}





void comp_surface_force(
  double Re,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **P,
  int **Flag,
  double *Fu,
  double *Fv
){

	int i,j;
	double pmin;

	*Fu=.0;
	*Fv=.0;

	pmin = P[1][1];
	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
			pmin = MIN(pmin, P[i][j]);

	for(i=1; i<=imax; i++)
		for(j=1; j<=jmax; j++)
			if(IS_BOUNDARY(Flag[i][j])	) {

				  
				  if( ( Flag[ i ][ j ] & B_N ) == B_N ){
					  *Fu += dx*(dUdy_hedge(U, dy, i, j) + dVdx_hedge(V,dx,i,j)) / Re;
					  *Fv += dx*( (2.0* ddy_bw(V,dy,i,j+1)) /Re - (P[i][j+1]-pmin) );
				  }
				
				  if( ( Flag[ i ][ j ] & B_S ) == B_S ) {
					  *Fu += -dx*(dUdy_hedge(U, dy, i, j-1) + dVdx_hedge(V,dx,i,j-1)) / Re;
					  *Fv += -dx*( (2.0* ddy_bw(V,dy,i,j-1)) /Re  - (P[i][j-1]-pmin) );
				  } 	
				 
				  if( ( Flag[ i ][ j ] & B_W ) == B_W ){
					  *Fu += -dy*( (2.0 * ddx_bw(U, dx, i-1, j)) / Re - (P[i-1][j]-pmin) );
					  *Fv += -dy*(dVdx_vedge(V, dx, i, j) + dUdy_vedge(U, dy, i-1, j)) /Re;
				  }
				  
				  if( ( Flag[ i ][ j ] & B_E ) == B_E ){
					  *Fu += dy*( (2.0 * ddx_bw(U, dx, i+1, j)) / Re - (P[i+1][j]-pmin) );
					  *Fv += dy*(dVdx_vedge(V, dx, i+1, j) + dUdy_vedge(U, dy, i, j)) /Re;
				  }

			}

}

