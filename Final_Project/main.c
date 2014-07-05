#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PARAMF "cavity.dat"
#define VISUAF "visual/sim"

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */




void has_nan(double **m, int imax, int jmax){
	int i,j;

	for(i=0; i<=imax+1; i++)
			for(j=0; j<=jmax+1; j++){
				if (m[i][j] != m[i][j])
					printf("[%d, %d] NaN\n", i,j);
				if (isinf(m[i][j]))
					printf("[%d, %d] Inf\n", i,j);
			}
}


int main(int argc, char** args){
	char pgm[50];
	double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, t, res, dp, nu;
	double **U, **V, **P, **F, **G, **RS;
	double **K, **E;					/* turbulent kinetic energy k, dissipation rate epsilon*/
	double Fu, Fv;						/* force integration variables */
	double KI, EI, cn, ce, c1, c2; 			/* K and E: Initial values for k and epsilon */
	int n, it, imax, jmax, itermax, pb;
	int fluid_cells;		/* Number of fluid cells in our geometry */
	char problem[10];		/* Problem name, file name */
	int boundaries[4];
	char *fname;

	int **Flag;			/* Flagflield matrix */

	if(argc>=2) 	fname=args[1];
	else 			fname = PARAMF;

	read_parameters(fname, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, boundaries, &dp, &pb, &KI, &EI, &cn, &ce, &c1, &c2, pgm);

	nu=1.0 / Re;
	KI = 0.003 * SQ(UI);
	EI = sqrt(KI*SQ(KI)) * cn / 0.03 / ylength;

	/* setting of the problem */
	switch (pb){
		case 0:	strcpy(problem,"karman");
		break;
		case 1:	strcpy(problem,"shear");
		break;
		case 2:	strcpy(problem,"step");
		break;
		default: strcpy(problem,"none");
		}


	/* Allocate Flag matrix */
	Flag = imatrix( 0, imax+1, 0, jmax+1 );

	/* should we change the dimension of the matrices in order to save space? */
	U = matrix ( 0 , imax+1 , 0 , jmax+1 );
	V = matrix ( 0 , imax+1 , 0 , jmax+1 );
	P = matrix ( 0 , imax+1 , 0 , jmax+1 );

	F = matrix ( 0 , imax , 0 , jmax );
	G = matrix ( 0 , imax , 0 , jmax );
	RS = matrix ( 0 , imax , 0 , jmax );

	K=matrix ( 0 , imax+1 , 0 , jmax+1 );
	E=matrix ( 0 , imax+1 , 0 , jmax+1 );
	
	init_flag( pgm, imax, jmax, &fluid_cells, Flag );
	init_uvp(UI, VI, PI, KI, EI, imax, jmax, U, V, P, K, E, Flag, problem);

	printf("Problem: %s\n", problem );
	printf( "xlength = %f, ylength = %f\n", xlength, ylength );
	printf( "imax = %d, jmax = %d\n", imax, jmax );
	printf( "dt = %f, dx = %f, dy = %f %f\n", dt, dx, dy, dt/dx/dy);
	printf( "Number of fluid cells = %d\n", fluid_cells );
	printf( "Reynolds number: %f\n", Re);

	t=.0;
	n=0;

	while( t <= t_end ){
		/*if( tau > 0 ) calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);*/

		boundaryvalues( imax, jmax, U, V, K, E, boundaries, Flag );

		/*has_nan(U,imax,jmax);
		has_nan(V,imax,jmax);
		has_nan(K,imax,jmax);
		has_nan(E,imax,jmax);*/

		/* special inflow boundaries, including k and eps */
		spec_boundary_val( problem, imax, jmax, U, V, K, E, Re, dp, cn, ylength);

		/*has_nan(U,imax,jmax);
		has_nan(V,imax,jmax);
		has_nan(K,imax,jmax);
		has_nan(E,imax,jmax);*/

		comp_KAEP(Re, nu, cn, ce, c1, c2, alpha, dt, dx, dy, imax, jmax, U, V, K, E, GX, GY, Flag);

		/*has_nan(U,imax,jmax);
		has_nan(V,imax,jmax);
		has_nan(K,imax,jmax);
		has_nan(E,imax,jmax);*/

		/* calculate new values for F and G */
		calculate_fg( Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, K, E, nu, cn, Flag );

		/*has_nan(F,imax-1,jmax-1);
		has_nan(G,imax-1,jmax-1);*/

		/* calculate right hand side */
		calculate_rs( dt, dx, dy, imax, jmax, F, G, RS, Flag );
		/*has_nan(RS,imax-1,jmax-1);*/

		it = 0;
		res = 10000.0;
		while( it < itermax && fabs(res) > eps ){
			sor( omg, dx, dy, imax, jmax, fluid_cells, P, RS, Flag, &res, problem, dp );
			/*has_nan(P,imax,jmax);*/
			it++;
		}

		printf("[%4d: %f] dt: %f, sor iterations: %4d \n", n, t, dt, it);

		calculate_uv( dt, dx, dy, imax, jmax, U, V, F, G, P, Flag );

		/*has_nan(U,imax,jmax);
		has_nan(V,imax,jmax);*/

		t += dt;
		n++;
	}

	write_vtkFile( VISUAF, 1, xlength, ylength, imax, jmax, dx, dy, U, V, P, K, E, Flag);

	comp_surface_force( Re, dx, dy, imax, jmax, U, V, P, Flag, &Fu, &Fv);

	printf("Problem: %s\n", problem );
	printf( "xlength = %f, ylength = %f\n", xlength, ylength );
	printf( "imax = %d, jmax = %d\n", imax, jmax );
	printf( "dt = %f, dx = %f, dy = %f %f\n", dt, dx, dy, dt/dx/dy);
	printf( "Number of fluid cells = %d\n", fluid_cells );
	printf( "Reynolds number: %f\n", Re);

	printf( "Drag force = %f Lift force = %f\n", Fu, Fv);

	/*int i,j;
	for(j = 0; j < jmax+1; j++) {
		    for(i = 0; i < imax+1; i++) {
		      fprintf(stderr, "%f ", (U[i][j] + U[i][j+1]) * 0.5);
		      if(i!=imax) fprintf(stderr, ", ");
		    }
		    fprintf(stderr, "\n");
		  }

	 printf("\n\n");
	 for(j = 0; j < jmax+1; j++) {
		    for(i = 0; i < imax+1; i++) {
		      fprintf(stderr, "%f ", (V[i][j] + V[i+1][j]) * 0.5 );
		      if(i!=imax) fprintf(stderr, ", ");
		    }
		    fprintf(stderr, "\n");
		  }*/


	/* free memory */
	free_matrix(U,0,imax+1,0,jmax+1);
	free_matrix(V,0,imax+1,0,jmax+1);
	free_matrix(P,0,imax+1,0,jmax+1);
	free_matrix(K,0,imax+1,0,jmax+1);
	free_matrix(E,0,imax+1,0,jmax+1);


	free_matrix(F,0,imax,0,jmax);
	free_matrix(G,0,imax,0,jmax);
	free_matrix(RS,0,imax,0,jmax);

	free_imatrix( Flag, 0, imax+1, 0, jmax+1 );

	return 0;
}
