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


void init_parallel (int iproc, int jproc, int imax, int jmax, int *myrank, int *il, int *ir, int *jb, int *jt, int *rank_l,
		int *rank_r, int *rank_b, int *rank_t, int *omg_i, int *omg_j, int num_proc){

	MPI_Status status;
	int i_block,j_block,i,j,app_il,app_ir,app_jb,app_jt,o_i,o_j, r_b, r_t, r_l, r_r, now;

	i_block=imax/iproc;		/* attention to the division  */
	j_block=jmax/jproc;		/* attention to the division  */


	if(*myrank==0){
		/*implementation for master thread*/
		/*
		*omg_i=1;
		*omg_j=1;
		*il=1;
		*ir=i_block;
		*jb=1;
		*jt=j_block;
		*/


		for( j=1; i<=jproc; j++){
			for (i=1; j<=iproc; i++){

				/* processor we are considering */
				now= i-1 +(j-1)*iproc;

				o_i=i;
				o_j=j;
				MPI_Send( &o_i, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( omg_i, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );
				MPI_Send( &o_j, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( omg_j, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );


				app_il=( num_proc%iproc ) * i_block + 1;
				app_ir=app_il + i_block -1;
				app_jb=( num_proc%jproc ) * j_block + 1;
				app_jt=app_jb + j_block -1;

				/* NB problem if not exact division !!! */

				MPI_Send( &app_il, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( il, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );

				MPI_Send( &app_ir, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( ir, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );

				MPI_Send( &app_jb, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( jb, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );

				MPI_Send( &app_jt, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( jt, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );

				/* ranks */
				r_l=now-1;
				r_r=now+1;
				r_t=now+iproc;
				r_b=now-iproc;

				if (i==1)
					r_l=MPI_PROC_NULL;		/* no left neighbor */

				if (i==iproc)
					r_r=MPI_PROC_NULL;		/* no right neighbor */

				if (j==1)
					r_b=MPI_PROC_NULL;		/* no bottom neighbor */

				if (j==jproc)
					r_t=MPI_PROC_NULL;		/* no top neighbor */

				/*
				MPI_Send( &r_l, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( rank_l, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );

				MPI_Send( &r_r, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( rank_r, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );

				MPI_Send( &r_t, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( rank_t, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );

				MPI_Send( &r_b, 1, MPI_INT, now, 1, MPI_COMM_WORLD );
				MPI_Recv( rank_b, 1, MPI_INT, 0, 1,MPI_COMM_WORLD, &status );

				*/

				MPI_Sendrecv(&r_l, 1 , MPI_INT, now, 1,MPI_COMM_WORLD, rank_l, 1, MPI_INT, 0, 1, MPI_COMM_WORLD , status);
				MPI_Sendrecv(&r_r, 1 , MPI_INT, now, 1,MPI_COMM_WORLD, rank_r, 1, MPI_INT, 0, 1, MPI_COMM_WORLD , status);
				MPI_Sendrecv(&r_b, 1 , MPI_INT, now, 1,MPI_COMM_WORLD, rank_b, 1, MPI_INT, 0, 1, MPI_COMM_WORLD , status);
				MPI_Sendrecv(&r_t, 1 , MPI_INT, now, 1,MPI_COMM_WORLD, rank_t, 1, MPI_INT, 0, 1, MPI_COMM_WORLD , status);

				/* CHECK sendrecive */
			}
		}

	}

}

void pressure_comm(double **P, int il,int ir, int jb, int jt, int rank_l, int rank_r, int rank_b, int rank_t, double *bufSend, double* bufRecv,MPI_Status *status, int chunck){

	/* 1- Check function send- receive
	 * 2- What the hell is chunk? Why do we need it??
	 *  */
	int i,j;

	bufSend=malloc(sizeof(*bufSend)* ( jt-jb+1 ) );
	bufRecv=malloc(sizeof(*bufSend)* ( jt-jb+1 ) );

	/* left */
	if(rank_l!=MPI_PROC_NULL){
		for (j=jb; j<=jt; j++){
			bufSend[j-jb]=P[il][j];
		}
		MPI_Sendrecv(bufSend, jt-jb+1 , MPI_DOUBLE, rank_l, 1,MPI_COMM_WORLD,bufRecv, jt-jb+1, MPI_DOUBLE,rank_l+1, 1, MPI_COMM_WORLD , status);
		/* setting left boundary with the new values */
		for (j=jb; j<=jt; j++){
					P[il-1][j]=bufRecv[j-jb];
				}

	}

	/* right */
		if(rank_r!=MPI_PROC_NULL){
			for (j=jb; j<=jt; j++){
				bufSend[j-jb]=P[ir][j];
			}
			MPI_Sendrecv(bufSend, jt-jb+1 , MPI_DOUBLE, rank_r, 1,MPI_COMM_WORLD, bufRecv, jt-jb+1, MPI_DOUBLE,rank_r-1, 1, MPI_COMM_WORLD , status);
			/* setting right boundary with the new values */
			for (j=jb; j<=jt; j++){
						P[ir+1][j]=bufRecv[j-jb];
					}

		}

	free(bufSend);
	free(bufRecv);

	bufSend=malloc(sizeof(*bufSend)* ( ir-il+1 ) );
	bufRecv=malloc(sizeof(*bufSend)* ( ir-il+1 ) );

		/* top */
		if(rank_t!=MPI_PROC_NULL){
			for (i=il; i<=ir; i++){
				bufSend[i-il]=P[i][jt];
			}
			MPI_Sendrecv(bufSend, ir-il+1 , MPI_DOUBLE, rank_t, 1,MPI_COMM_WORLD, bufRecv, ir-il+1, MPI_DOUBLE,rank_t-chuck, 1, MPI_COMM_WORLD , status);
			/* needed here chunk?? it would be equal to iproc...*/
			/* setting top boundary with the new values */
			for (i=il; i<=ir; i++){
							P[i][jt+1]=bufRecv[i-il];
						}
		}

		/* bottom */
		if(rank_b!=MPI_PROC_NULL){
			for (i=il; i<=ir; i++){
				bufSend[i-il]=P[i][jb];
			}
			MPI_Sendrecv(bufSend, ir-il+1 , MPI_DOUBLE, rank_b, 1,MPI_COMM_WORLD, bufRecv, ir-il+1, MPI_DOUBLE,rank_b+chunck, 1, MPI_COMM_WORLD , status);
			/* setting bottom boundary with the new values */
			for (i=il; i<=ir; i++){
							P[i][jb-1]=bufRecv[i-il];
						}
		}
		free(bufSend);
		free(bufRecv);

}
