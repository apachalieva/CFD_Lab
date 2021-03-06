#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define VISUA_FOLDER "./visual/"
#define CONFIGS_FOLDER "./../configs/"

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


int main(int argc, char** args){
	
	double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, t, res, dp, nu;
	double **U, **V, **P, **F, **G, **RS;
	double **K;					/* turbulent kinetic energy k */
	double **W; 				/* omega: dissipation frequency eps/K */
	double **delta;				/* matrix of the distances */
	double Fu, Fv;				/* force integration variables */
	double KI, EI, WI, cn, ce, c1, c2; 			/* K and E: Initial values for k and epsilon */

	int n, it, imax, jmax, itermax, pb, boundaries[4];
	int fluid_cells;		/* Number of fluid cells in our geometry */
	int **Flag;			/* Flagflield matrix */
	
	char vtkname[200];
	char pgm[200];
	char problem[10];		/* Problem name */
	char fname[200];

	if(argc<2){
		printf("No parameter file specified. Terminating...\n");
		exit(1);
	}

	sprintf(fname, "%s%s", CONFIGS_FOLDER, args[1]);
	printf("%s\n",fname);

	read_parameters(fname, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &itermax, &eps, boundaries, &dp, &pb, &KI, &WI, &cn, &ce, &c1, &c2, pgm, &nu, problem);

	/* Allocate Flag matrix */
	Flag = imatrix( 0, imax+1, 0, jmax+1 );
	/* Allocate delta matrix */
	delta = matrix( 0, imax+1, 0, jmax+1 );

	U = matrix ( 0 , imax+1 , 0 , jmax+1 );
	V = matrix ( 0 , imax+1 , 0 , jmax+1 );
	P = matrix ( 0 , imax+1 , 0 , jmax+1 );

	F = matrix ( 0 , imax , 0 , jmax );
	G = matrix ( 0 , imax , 0 , jmax );
	RS = matrix ( 0 , imax , 0 , jmax );

	K=matrix ( 0 , imax+1 , 0 , jmax+1 );
	W=matrix ( 0 , imax+1 , 0 , jmax+1 );
	
	/* Initialize values to the Flag, u, v and p */
	init_flag(CONFIGS_FOLDER,pgm, imax, jmax, &fluid_cells, Flag );
	init_delta(Flag,delta,imax,jmax,dx,dy);
	init_uvp(UI, VI, PI, KI, WI, imax, jmax, U, V, P, K, W, Flag, problem);

	printf("Problem: %s\n", problem );
	printf( "xlength = %f, ylength = %f\n", xlength, ylength );
	printf( "imax = %d, jmax = %d\n", imax, jmax );
	printf( "dt = %f, dx = %f, dy = %f\n", dt, dx, dy);
	printf( "Number of fluid cells = %d\n", fluid_cells );
	printf( "Reynolds number: %f\n\n", Re);

	t=.0;
	n=0;

	while( t <= t_end ){
		boundaryvalues( imax, jmax, dx, dy, U, V, K, W, nu, boundaries, Flag );

		/* special inflow boundaries, including k and eps */
		spec_boundary_val( problem, imax, jmax, U, V, K, W, Re, dp, cn, ylength);

		/* calculate new values for k and w */
		comp_KAW(Re, nu, ce, c1, c2, alpha, dt, dx, dy, imax, jmax, U, V, K, W, GX, GY, delta, Flag);

		/* calculate new values for F and G */
		calculate_fg( Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G, K, W, nu, delta,Flag );

		/* calculate right hand side */
		calculate_rs( dt, dx, dy, imax, jmax, F, G, RS, Flag );

		it = 0;
		res = 10000.0;
		while( it < itermax && fabs(res) > eps ){
			sor( omg, dx, dy, imax, jmax, fluid_cells, P, RS, Flag, &res, problem, dp );
			it++;
		}

		printf("[%5d: %f] dt: %f, sor iterations: %4d \n", n, t, dt, it);

		/* calculate new values for u and v */
		calculate_uv( dt, dx, dy, imax, jmax, U, V, F, G, P, Flag );

		t += dt;
		n++;
	}

	sprintf(vtkname, "%s%s", VISUA_FOLDER,  args[1]);
	write_vtkFile( vtkname, 1, xlength, ylength, imax, jmax, dx, dy, U, V, P, K, W, Flag);
	
	comp_surface_force( Re, dx, dy, imax, jmax, U, V, P, Flag, &Fu, &Fv);

	printf( "\nProblem: %s\n", problem );
	printf( "xlength = %f, ylength = %f\n", xlength, ylength );
	printf( "imax = %d, jmax = %d\n", imax, jmax );
	printf( "dt = %f, dx = %f, dy = %f\n", dt, dx, dy);
	printf( "Number of fluid cells = %d\n", fluid_cells );
	printf( "Reynolds number: %f\n", Re);
	printf( "Drag force = %f Lift force = %f\n", Fu, Fv);

	/* free memory */
	free_matrix(U,0,imax+1,0,jmax+1);
	free_matrix(V,0,imax+1,0,jmax+1);
	free_matrix(P,0,imax+1,0,jmax+1);
	free_matrix(K,0,imax+1,0,jmax+1);
	free_matrix(W,0,imax+1,0,jmax+1);
	free_matrix(delta,0,imax+1,0,jmax+1);


	free_matrix(F,0,imax,0,jmax);
	free_matrix(G,0,imax,0,jmax);
	free_matrix(RS,0,imax,0,jmax);

	free_imatrix(Flag,0,imax+1,0,jmax+1);

	return 0;
}
