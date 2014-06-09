#include "parallel.h"
#include <mpi.h>



void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}


void init_parallel (int iproc, int jproc, int imax, int jmax, int *myrank, int *il, int *ir, int *jt, int *jb, int *rank_l,
		int *rank_r, int *rank_b, int *rank_t, int *omg_i, int *omg_j, int num_proc){

	int vars[10];

	if(*myrank==0){

		int i_block, j_block, i, j, now, jb_end;

		i_block=imax/iproc;
		j_block=jmax/jproc;

		int il_base = 1, jb_base=1;


		for( j=1; j<=jproc; j++){

			/* calculate end of partition in the y direction */
			jb_end = jb_base + j_block;
			if(j <= jmax%jproc)
				jb_end++;

			il_base = 1;

			for (i=1; i<=iproc; i++){

				/* processor we are considering */
				now= i-1 +(j-1)*iproc;

				vars[0] = il_base; /*il*/

				il_base += i_block;
				if(i <= imax%iproc)
					il_base++;

				vars[1] = il_base-1; /*ir*/

				vars[2] = jb_base;	 /* jb */
				vars[3] = jb_end-1; /* jt*/

				vars[4] = 	i>1 ? now-1 : MPI_PROC_NULL;	/*rank_l*/
				vars[5] = 	i<iproc ? now+1 : MPI_PROC_NULL;	/*rank_r*/
				vars[6] = 	j>1 ? now-iproc : MPI_PROC_NULL;	/*rank_b*/
				vars[7] = 	j<jproc ? now+iproc : MPI_PROC_NULL;	/*rank_t*/
				vars[8] = 	i;	/*omg_i*/
				vars[9] = 	j;	/*omg_j*/

				if(i!=1 || j!=1){
					MPI_Send( vars, 10, MPI_INT, now, 0, MPI_COMM_WORLD );
				}else{

					/* the parameters are of the master process */
					*il = vars[0];
					*ir = vars[1];
					*jb = vars[2];
					*jt = vars[3];
					*rank_l = vars[4];
					*rank_r = vars[5];
					*rank_b = vars[6];
					*rank_t = vars[7];
					*omg_i = vars[8];
					*omg_j  = vars[9];

				}

			}

			jb_base = jb_end;
		}

	} else {

		MPI_Status status;
		MPI_Recv( vars, 10, MPI_INT, 0, 0,MPI_COMM_WORLD, &status );

		/* assign variable in order */
		*il = vars[0];
		*ir = vars[1];
		*jb = vars[2];
		*jt = vars[3];
		*rank_l = vars[4];
		*rank_r = vars[5];
		*rank_b = vars[6];
		*rank_t = vars[7];
		*omg_i = vars[8];
		*omg_j  = vars[9];

	}

}


void pressure_comm(double **P, int il,int ir, int jt, int jb, int rank_l, int rank_r, int rank_b, int rank_t, double *bufSend, double* bufRecv,MPI_Status *status){
	int i,j;

	/* send left, receive right */
	if(rank_l!=MPI_PROC_NULL)
		for (j=jb; j<=jt; j++)
			bufSend[j-jb]=P[il][j];

	MPI_Sendrecv(bufSend, jt-jb+1 , MPI_DOUBLE, rank_l, 1,bufRecv, jt-jb+1, MPI_DOUBLE,rank_r, 1, MPI_COMM_WORLD , status);

	/* setting right boundary with the new values */
	if(rank_r!=MPI_PROC_NULL)
		for (j=jb; j<=jt; j++)
			P[ir+1][j]=bufRecv[j-jb];

	/* sending right, receive left */
	if(rank_r!=MPI_PROC_NULL)
		for (j=jb; j<=jt; j++)
			bufSend[j-jb]=P[ir][j];

	MPI_Sendrecv(bufSend, jt-jb+1 , MPI_DOUBLE, rank_r, 2, bufRecv, jt-jb+1, MPI_DOUBLE,rank_l, 2, MPI_COMM_WORLD , status);

	/* setting left boundary with the new values */
	if(rank_l!=MPI_PROC_NULL)
		for (j=jb; j<=jt; j++)
			P[il-1][j]=bufRecv[j-jb];

	/* send top, receive bottom */
	if(rank_t!=MPI_PROC_NULL)
		for (i=il; i<=ir; i++)
			bufSend[i-il]=P[i][jt];

	MPI_Sendrecv(bufSend, ir-il+1 , MPI_DOUBLE, rank_t, 3, bufRecv, ir-il+1, MPI_DOUBLE,rank_b, 3, MPI_COMM_WORLD , status);

	/* setting bottom boundary with the new values */
	if(rank_b!=MPI_PROC_NULL)
		for (i=il; i<=ir; i++)
			P[i][jb-1]=bufRecv[i-il];

	/* send bottom, receive top */
	if(rank_b!=MPI_PROC_NULL)
		for (i=il; i<=ir; i++)
			bufSend[i-il]=P[i][jb];

	MPI_Sendrecv(bufSend, ir-il+1 , MPI_DOUBLE, rank_b, 4, bufRecv, ir-il+1, MPI_DOUBLE,rank_t, 4, MPI_COMM_WORLD , status);

	/* setting top boundary with the new values */
	if(rank_t!=MPI_PROC_NULL)
		for (i=il; i<=ir; i++)
			P[i][jt+1]=bufRecv[i-il];

}

void uv_comm(double **U, double **V, int il, int ir, int jt, int jb, int rank_l, int rank_r, int rank_b, int rank_t, double *bufSend, double *bufRecv, MPI_Status *status){

	int i,j, u_bsize, v_bsize, bsize;

	u_bsize =  jt-jb+1;
	v_bsize =  jt-jb+2;
	bsize = u_bsize + v_bsize;

	/* send left, receive right */
	if(rank_l!=MPI_PROC_NULL){
		for (j=0; j<v_bsize; j++)
			bufSend[j] = V[il][jb+j-1];
		for (j=0; j<u_bsize; j++)
			bufSend[j+v_bsize] = U[il][jb+j];
	}

	MPI_Sendrecv(bufSend, bsize , MPI_DOUBLE, rank_l, 5,bufRecv, bsize, MPI_DOUBLE,rank_r, 5, MPI_COMM_WORLD , status);

	/* setting right boundary with the new values */
	if(rank_r!=MPI_PROC_NULL){
		for (j=0; j<v_bsize; j++)
			V[ir+1][jb+j-1] = bufRecv[j];
		for (j=0; j<u_bsize; j++)
			U[ir+1][jb+j] = bufRecv[j+v_bsize];
	}


	/* sending right, receive left */
	if(rank_r!=MPI_PROC_NULL){
		for (j=0; j<v_bsize; j++)
			bufSend[j] = V[ir][jb+j-1];
		for (j=0; j<u_bsize; j++)
			bufSend[j+v_bsize] = U[ir-1][jb+j];
	}

	MPI_Sendrecv(bufSend, bsize, MPI_DOUBLE, rank_r, 6, bufRecv, bsize, MPI_DOUBLE,rank_l, 6, MPI_COMM_WORLD , status);

	/* setting left boundary with the new values */
	if(rank_l!=MPI_PROC_NULL){
		for (j=0; j<v_bsize; j++)
			V[il-1][jb+j-1] = bufRecv[j];
		for (j=0; j<u_bsize; j++)
			U[il-2][jb+j] = bufRecv[j+v_bsize];
	}


	u_bsize =  ir-il+2;
	v_bsize =  ir-il+1;
	bsize = u_bsize + v_bsize;

	/* send top, receive bottom */
	if(rank_t!=MPI_PROC_NULL){
		for (i=0; i<v_bsize; i++)
			bufSend[i] = V[il+i][jt-1];
		for (i=0; i<u_bsize; i++)
			bufSend[i+v_bsize] = U[il-1+i][jt];
	}

	MPI_Sendrecv(bufSend, bsize , MPI_DOUBLE, rank_t, 7, bufRecv, bsize, MPI_DOUBLE,rank_b, 7, MPI_COMM_WORLD , status);

	/* setting bottom boundary with the new values */
	if(rank_b!=MPI_PROC_NULL){
		for (i=0; i<v_bsize; i++)
			V[il+i][jb-2] = bufRecv[i];
		for (i=0; i<u_bsize; i++)
			U[il-1+i][jb-1] = bufRecv[i+v_bsize];
	}


	/* send bottom, receive top */
	if(rank_b!=MPI_PROC_NULL){
		for (i=0; i<v_bsize; i++)
			bufSend[i] = V[il+i][jb];
		for (i=0; i<u_bsize; i++)
			bufSend[i+v_bsize] = U[il-1+i][jb];
	}

	MPI_Sendrecv(bufSend, bsize , MPI_DOUBLE, rank_b, 8, bufRecv, bsize, MPI_DOUBLE,rank_t, 8, MPI_COMM_WORLD , status);

	/* setting top boundary with the new values */
	if(rank_t!=MPI_PROC_NULL){
		for (i=0; i<v_bsize; i++)
			V[il+i][jt+1] = bufRecv[i];
		for (i=0; i<u_bsize; i++)
			U[il-1+i][jt+1] = bufRecv[i+v_bsize];
	}
}

