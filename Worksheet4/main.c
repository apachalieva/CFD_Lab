#include "helper.h"
#include "visual.h"
#include "init.h"
#include "uvp.h"
#include "boundary_val.h"
#include "sor.h"
#include <stdio.h>
#include "parallel.h"



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

int main(int argc, char** argv){

	double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, t, res;
	double **U, **V, **P, **F, **G, **RS, *bufSend, *bufRecv;
	int n, it, step, imax, jmax, itermax;
	int il, ir, jb, jt, rank_l, rank_r, rank_t, rank_b, omg_i, omg_j, myrank, nproc, iproc, jproc;			/* new int */
	MPI_Status status;
	char strbuf[650];

	   /* initialisation */
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &nproc );     /* asking for the number of processes  */
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );    /* asking for the local process id   */

	/* read parameters from file */
	read_parameters(PARAMF, &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &iproc, &jproc );

	/* compute and spread computation parameters */
	init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jt, &jb, &rank_l, &rank_r, &rank_b, &rank_t, &omg_i, &omg_j, nproc);

	if(myrank==0){
		sprintf(strbuf, "Processors: %d", nproc);
		Program_Message(strbuf);
	}

	sprintf(strbuf, "myrank %d, il %d, ir %d, jb %d, jt %d, rank_l %d, rank_r %d, rank_b %d, rank_t %d, omg_i %d, omg_j %d", myrank, il, ir, jb, jt, rank_l,rank_r, rank_b, rank_t, omg_i, omg_j);
	Program_Message(strbuf);

	/* changes in the dimensions of the matrices */
	U=matrix ( il-2 , ir+1 , jb-1 , jt+1 );
	V=matrix ( il-1 , ir+1 , jb-2 , jt+1 );
	P=matrix ( il-1 , ir+1 , jb-1 , jt+1 );

	F=matrix ( il-2 , ir+1 , jb-1 , jt+1 );
	G=matrix ( il-1 , ir+1 , jb-2 , jt+1 );
	RS=matrix ( il , ir , jb , jt );

	init_uvp(UI, VI, PI, il, ir, jt, jb, U, V, P);

	/* allocate communication buffers once with the maximum needed dimension */
	bufSend=malloc(sizeof(*bufSend)*  max(ir-il+2 + ir-il+1, jt-jb+1 + jt-jb+2 ) );
	bufRecv=malloc(sizeof(*bufRecv)*  max(ir-il+2 + ir-il+1, jt-jb+1 + jt-jb+2 ) );

	t=.0;
	n=0;
	step=0;

	while(t<t_end){
		if(tau>0) calculate_dt(Re, tau, &dt, dx, dy, il, ir, jt, jb, bufSend, bufRecv, imax, jmax, U, V);

		boundaryvalues( il, ir, jt, jb, imax, jmax, U, V);
		calculate_fg(Re, GX, GY, alpha, dt, dx, dy, il, ir, jt, jb, imax, jmax, U, V, F, G);
		calculate_rs(dt, dx, dy, il, ir, jt, jb, F, G, RS);

		it=0;
		res=10000.0;
		while(it < itermax && fabs(res) > eps) {
			sor(omg, dx, dy, il, ir, jt, jb, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, 0, imax, jmax, P, RS, &res);
			it++;
		}


		if (myrank==0){
			sprintf(strbuf, "[t %f, dt %f, sor_iterations %d] sor terminated with res = %f ", t, dt, it, res);
			Program_Message(strbuf);
		}

		/* compute and calculate uv */
		calculate_uv(dt, dx, dy, il, ir, jt, jb, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, &status, 0, imax, jmax, U, V, F, G, P);

		t += dt;
		n++;

		if(step*dt_value <= t){
			/* write vtk file using the interval specified in the parameters file */
			write_vtkFile(VISUAF, n, xlength, ylength, il, ir, jt, jb, myrank, dx, dy, U, V, P );
			step++;
		}

	}


	/* modification in the dimensions */
	free_matrix (U, il-2 , ir+1 , jb-1 , jt+1 );
	free_matrix (V, il-1 , ir+1 , jb-2 , jt+1 );
	free_matrix(P, il-1 , ir+1 , jb-1 , jt+1 );

	free_matrix(F, il-2 , ir+1 , jb-1 , jt+1 );
	free_matrix(G, il-1 , ir+1 , jb-2 , jt+1 );
	free_matrix(RS, il , ir , jb , jt );

	/* deallocate comunication buffers */
	free(bufSend);
	free(bufRecv);

	Programm_Stop("finished");

	return 0;
}
