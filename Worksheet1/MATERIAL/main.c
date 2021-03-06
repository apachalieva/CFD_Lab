#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>


#define PARAMF "cavity100.dat"
#define VISUAF "visual/cavity"


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
int main(int argn, char** args){
	double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, t, res;
	double **U, **V, **P, **F, **G, **RS;
	int n, it, imax, jmax, itermax,step=0;

	read_parameters(PARAMF, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value );

	U=matrix ( 0 , imax+1 , 0 , jmax+1 );
	V=matrix ( 0 , imax+1 , 0 , jmax+1 );
	P=matrix ( 0 , imax+1 , 0 , jmax+1 );

	F=matrix ( 0 , imax , 0 , jmax );
	G=matrix ( 0 , imax , 0 , jmax );
	RS=matrix ( 0 , imax , 0 , jmax );


	init_uvp(UI, VI, PI, imax, jmax, U, V, P);

	t=.0;
	n=0;

	while(t<t_end){
		if(tau>0) calculate_dt(Re, tau, &dt, dx, dy, imax, jmax, U, V);

		boundaryvalues(imax, jmax, U, V);
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G);
		calculate_rs(dt, dx, dy, imax, jmax, F, G, RS);

		it = 0;
		res=10000.0;
		while(it < itermax && fabs(res) > eps){
			sor(omg, dx, dy, imax, jmax, P, RS, &res);
			it++;
		}

		printf("[%d: %f] dt: %f, sor iterations: %d \n", n, t, dt, it);

		calculate_uv(dt, dx, dy, imax, jmax, U, V, F, G, P);
		
		t += dt;
		n++;
		
		if(step*dt_value <= t){
			/* output vtk file for visualization */
			write_vtkFile( VISUAF, n, xlength, ylength, imax, jmax, dx, dy, U, V, P );
			step++;
		}		
		
	}

	/*printf("U[imax/2][7*jmax/8]: %f", U[imax/2][7*jmax/8]);
	*/


	free_matrix(U,0,imax+1,0,jmax+1);
	free_matrix(V,0,imax+1,0,jmax+1);
	free_matrix(P,0,imax+1,0,jmax+1);

	free_matrix(F,0,imax,0,jmax);
	free_matrix(G,0,imax,0,jmax);
	free_matrix(RS,0,imax,0,jmax);

	return 0;
}
