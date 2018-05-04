/* mm.c 
 * Sample matrix multiplication in C 
 * Author: Philip Papadopoulos
 * Email: ppapadopoulos@ucsd.edu
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define DATA 10
#define DONE 0

double cpu_time ( void );
void mult(double **, double **, double **, int, int, int);
void multM(double ****, double ****, double ****, int, int, int, int, int, int*, double*);
void add(double **, double *, int);
void zero(double *, int);
void fill_m(double**, double**, double*, int, int, int);
void fill_w(double**, double*, int, int, int, double);
void unload_w(double*, double**, double**, int);
void sendToWorker(int, double*, MPI_Request*, int, int);
int waitForOne(double****, MPI_Request*, int, double*, int, double*);
void waitForAll(double**** , MPI_Request*, int, int, int, double*, double*);
void usage( void );


/* ===========================================================================
 * Usage:  mm <N> <trials> 
 * New Usage: mm <N> <M> <P> <trials>
 */
int main(int argc, char * argv[]){
  int blkSize=atoi(argv[1]);
	int N=atoi(argv[2]); 
  int M=atoi(argv[3]);
  int P=atoi(argv[4]);
  FILE *fp=NULL;
	char* outfile=argv[5];
  int comm_sz;
  int my_rank;
	int i;
  int j;
  double startTime;
  double endTime;

	double ****A, ****B, ****C;
	double ***Alinear, ***Blinear, ***Clinear;
  int* blksPerRank;
  double * blksTime;

	if (argc != 6){
    if(my_rank == 0){
      usage();
    }
    MPI_Finalize();
    return -1;
  }


  // set up MPI components
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // master code
  if(my_rank == 0){
    fp = fopen(outfile, "w");
    // A is NxM, B is MxP, C is NxP
	  A = (double ****) calloc(N, sizeof (double ***));
	  B = (double ****) calloc(M, sizeof (double ***));
	  C = (double ****) calloc(N, sizeof (double ***));
	  Alinear = calloc(N * M, sizeof (double**));
	  Blinear = calloc(M * P, sizeof (double**));
	  Clinear = calloc(N * P, sizeof (double**));
    blksPerRank = calloc(comm_sz, sizeof(int));
    blksTime = calloc(comm_sz, sizeof(double));

	  /* populate A, B ,C so that we can use  A[i][j] addressing */
	  for (i=0; i<N ; i++){ 
		  A[i] = Alinear + i * M ;
		  C[i] = Clinear + i * P;
	  }
    for (i=0; i<M; i++){
		  B[i] = Blinear + i * P;
    }
  
    double * blkLiner;
    // zero fill matrix C
    for (i=0; i<N*P; i++){
      Clinear[i] = calloc(blkSize, sizeof(double*));
      blkLiner = calloc(blkSize * blkSize, sizeof(double));
      for (j=0; j<blkSize; j++){
        Clinear[i][j] = blkLiner + j * blkSize;
      }
    }
  

    // fill random values into A and B
	  for (i=0; i<N*M; i++){
		  Alinear[i] = calloc(blkSize, sizeof(double*));
		  blkLiner  = calloc(blkSize * blkSize, sizeof(double));
      for (j=0; j<blkSize; j++){
		    Alinear[i][j] = blkLiner + j * blkSize;
      }
      for(j=0; j<blkSize*blkSize; j++){
        blkLiner[j] = drand48();
      }
	  }
    for (i=0; i<M*P; i++){
		  Blinear[i] = calloc(blkSize, sizeof(double*));
		  blkLiner  = calloc(blkSize * blkSize, sizeof(double));
      for(j=0; j<blkSize; j++){
        Blinear[i][j] = blkLiner + j * blkSize;
      }
      for(j=0; j<blkSize*blkSize; j++){
        blkLiner[j] = drand48();
      }
    }
  
  // debug print, remember to comment out
  int m,n;
  double** pblock;
  fprintf(fp, "%d,%d,%d", N, M, P);
  fprintf(fp, "Below is A:\n");
  for (i=0; i<N; i++){
    for (j=0; j<M; j++){
      fprintf(fp, "A[%d][%d] is block:\n", i, j);
      pblock = A[i][j];
      for(m=0; m<blkSize; m++){
        for(n=0; n<blkSize; n++){
          fprintf(fp, "%f ", pblock[m][n]);
        }
        fprintf(fp, "\n");
      }
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\nBelow is B:\n");
  for (i=0; i<M; i++){
    for (j=0; j<P; j++){
      fprintf(fp, "B[%d][%d] is block:\n", i, j);
      pblock = B[i][j];
      for(m=0; m<blkSize; m++){
        for(n=0; n<blkSize; n++){
          fprintf(fp, "%f ", pblock[m][n]);
        }
        fprintf(fp, "\n");
      }
    }
    fprintf(fp,"\n");
  }

    startTime = cpu_time();
	  multM(A,B,C,N,M,P,blkSize, comm_sz, blksPerRank, blksTime);
    endTime = cpu_time();

  // debug print, remeber to comment out 
  fprintf(fp, "Below is C:\n");
  for (i=0; i<N; i++){
    for (j=0; j<P; j++){
      pblock = C[i][j];
      fprintf(fp, "C[%d][%d] is:\n", i, j);
      for(m=0; m<blkSize; m++){
        for(n=0; n<blkSize; n++){
          fprintf(fp, "%f ", pblock[m][n]);
        }
        fprintf(fp, "\n");
      }
    }
    fprintf(fp, "\n");
  }

    printf("Parameters: %d, %d, %d, %d, %d\n", blkSize, N, M, P, comm_sz-1);
    for(i=1; i<comm_sz; i++){
      printf("Worker %d:  %d, %15.8f\n", i, blksPerRank[i], blksTime[i]);
    }
    printf("Total Time: %15.8f\n", endTime-startTime);


    fclose(fp);
	  free(A);
	  free(Alinear);
	  free(B);
	  free(Blinear);
	  free(C);
	  free(Clinear);
  }

  // Worker set
  else{
    MPI_Status status;
    //MPI_Request * request = alloca((blkSize*blkSize*blkSize) * sizeof(MPI_Request));
	  double** X = (double **) calloc(blkSize, sizeof (double *));
	  double** Y = (double **) calloc(blkSize, sizeof (double *));
	  double** Z = (double **) calloc(blkSize, sizeof (double *));
    for(i=0; i<blkSize; i++){
      X[i] = (double*) calloc(blkSize, sizeof(double));
      Y[i] = (double*) calloc(blkSize, sizeof(double));
      Z[i] = (double*) calloc(blkSize, sizeof(double));
    }
    double* recvLiner = calloc(blkSize*blkSize*2+2, sizeof (double));
    double* sendLiner = calloc(blkSize*blkSize+3, sizeof (double));
    int done = 0;
    int tag;
    int row, col;
    int index=0;
    while(!done){

      //printf("WORKER[%d]: task_%d before receive\n", my_rank, index);

      MPI_Recv(&recvLiner[0], blkSize*blkSize*2+2, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      //printf("WORKER[%d]: task_%d after receive\n", my_rank, index);
      tag = status.MPI_TAG;
      // receive the done message
      if(tag == DONE){
        // free
        //printf("WORKER[%d] receive done\n", my_rank);
        done = 1;
      }
      // receive a task
      else{
        row = recvLiner[blkSize*blkSize*2];
        col = recvLiner[blkSize*blkSize*2+1];

  
        /*printf("\nWORKER[%d] send to row:%d, col:%d\n", my_rank, row, col);
        for(i=0; i<blkSize*blkSize; i++){
          printf("%f ", recvLiner[i]);
        }
        printf("\n the row, col is:\n");
        for(i=blkSize*blkSize*2; i<blkSize*blkSize*2+2; i++){
          printf("%f ", recvLiner[i]);
        }
        printf("\n");*/



        unload_w(recvLiner, X, Y, blkSize);
        startTime = cpu_time();
        mult(X, Y, Z, blkSize, blkSize, blkSize);
        endTime = cpu_time();
        double duration = endTime - startTime;


        fill_w(Z, sendLiner, blkSize, row, col, duration);


        /* debug, remember to comment out
        printf("WORKER[%d] calculate matrix\n", my_rank);
        for(i=0; i<blkSize*blkSize; i++){
          printf("%f ", recvLiner[i]);
        }
        printf("\n");*/

        //MPI_Isend(&recvLiner[0], blkSize*blkSize+2, MPI_DOUBLE, 0, DATA, MPI_COMM_WORLD, &request[index]);
        MPI_Send(&sendLiner[0], blkSize*blkSize+3, MPI_DOUBLE, 0, DATA, MPI_COMM_WORLD);
        index++;
      }
      //printf("WORKER[%d]: task_%d after computation\n", my_rank, index);
    }
    //MPI_Waitall(index, request, &status);
    //printf("WORKER [%d] after wait\n", my_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    free(X);
    free(Y);
    free(Z);
    //for(i=0; i<blkSize; i++){
    //  free(X[i]);
    //}
    free(recvLiner);
    free(sendLiner);

  }

  MPI_Finalize();
	return 0;
}


void fill_w(double**temp, double*liner, int size, int row, int col, double duration){
  int i,j;
  for(i=0; i<size; i++){
    for(j=0; j<size; j++){
      liner[i*size+j] = temp[i][j];
    }
  }
  liner[size*size] = row;
  liner[size*size+1] = col;
  liner[size*size+2] = duration;
}


void unload_w(double*liner, double**A, double**B, int size){
  int i,j;
  for(i=0; i<size; i++){
    for(j=0; j<size; j++){
      A[i][j] = liner[i*size+j];
    }
  }
  for(i=0; i<size; i++){
    for(j=0; j<size; j++){
      B[i][j] = liner[size*size + i*size+j];
    }
  }
}


/** mult
 * Multiply matrices, C = A x B.
 * Assumes Matrices are already allocated.
 * values are not checked for sanity
*/
void mult(double **A, double **B, double **C, int N, int M, int P)
{
	int i, row, col;
	double sum;

  /* printf("in mult, A is:\n");
  for(row=0; row<N; row++){
    for(col=0; col<M; col++){
      printf("%f ", A[row][col]);
    }
    printf("\n");
  }

  printf("in mult, B is:\n");
  for(row=0; row<M; row++){
    for(col=0; col<P; col++){
      printf("%f ", B[row][col]);
    }
    printf("\n");
  }*/

	for (row = 0; row < N; row ++){
		for (col = 0; col < P; col ++){
			sum = 0.0;
			for (i = 0; i < M ; i++){
				sum += A[row][i] * B[i][col];
			}
			C[row][col] = sum;
		}
	}

  /* printf("in mult, C is:\n");
  for(row=0; row<N; row++){
    for(col=0; col<P; col++){
      printf("%f ", C[row][col]);
    }
    printf("\n");
  }*/

}

void multM(double ****A, double ****B, double ****C, int N, int M, int P, int blkSize, int np, int* blksPerRank, double*blksTime){
  int i, row, col;
  double** temp = (double **) calloc(blkSize, sizeof (double *));
	double* tempLiner = calloc(blkSize*blkSize, sizeof (double));
  // fill in the send buffer
  double* sendLiner = calloc(blkSize*blkSize*2+2, sizeof (double));
  double* recvLiner = calloc(blkSize*blkSize+3, sizeof (double));

  for (i=0; i<blkSize; i++){
		temp[i] = tempLiner + i * blkSize;
  }

  int startPhase = np-1;
  int nextWorker;
  int activeRequests = 0;
  MPI_Request * recvRqsts = alloca((np+1) * sizeof(MPI_Request));
  for(i=0; i<np; i++){
    recvRqsts[i] = MPI_REQUEST_NULL;
  }
  
  // iterative distribute task
  for (row = 0; row < N; row ++){
    for (col = 0; col < P; col ++){
      zero(tempLiner, blkSize);
      for (i = 0; i < M; i++){

        //printf("MASTER[%d]: i=%d\n", 0, i);
        //printf("M is %d\n", M);
        // start phase, handout mult task
        if (startPhase >0){
          nextWorker = startPhase --;
        }
        else{
          //printf("MASTER[%d]: before receive, %d\n", 0, i);
          nextWorker = waitForOne(C, recvRqsts, np, recvLiner, blkSize, blksTime);
          activeRequests --;
          //printf("MASTER[%d]: after receive, %d\n", 0, i);
        }
        
        // send out task
        fill_m(A[row][i], B[i][col], sendLiner, blkSize, row, col);

        /*printf("\nMASTER send to row:%d, col:%d\n", row, col);
        for(j=0; j<blkSize*blkSize; j++){
          printf("%f ", recvLiner[j]);
        }
        printf("\n the lower part is:\n");
        for(j=blkSize*blkSize; j<blkSize*blkSize*2+2; j++){
          printf("%f ", sendLiner[j]);
        }
        printf("\n");*/

        //printf("MASTER[%d]: before send, %d\n", 0, i);
        sendToWorker(nextWorker, sendLiner, recvRqsts, blkSize, i);
        //printf("MASTER[%d]: after send, %d\n", 0, i);
        activeRequests ++;

        // record how many blks worker processed
        blksPerRank[nextWorker]++;
      }
    }
  }
  waitForAll(C, recvRqsts, activeRequests, np, blkSize, recvLiner, blksTime);



  // bcast the done message
  for(i=0; i<np; i++){
    recvRqsts[i] = MPI_REQUEST_NULL;
  }
  for(i=0; i<np; i++){
    //printf("send done message to worker:%d\n", i);
    MPI_Isend(&sendLiner[0], blkSize*blkSize*2+2, MPI_DOUBLE, i, DONE, MPI_COMM_WORLD, &recvRqsts[i]);
  }
  //printf("after send done message\n");
  //MPI_Waitall(np, recvRqsts, &status);
  MPI_Barrier(MPI_COMM_WORLD);

  // free
  free(temp);
  free(tempLiner);
  //free(recvRqsts);
  free(sendLiner);
  free(recvLiner);
}


void sendToWorker(int nextWorker, double* sendLiner, MPI_Request * recvRqsts, int blkSize, int i){
  MPI_Send(&sendLiner[0], blkSize*blkSize*2+2, MPI_DOUBLE, nextWorker, DATA, MPI_COMM_WORLD);
}

void waitForAll(double**** C, MPI_Request* recvRqsts, int activeRequests, int np, int blkSize, double*recvLiner, double* blksTime){
  //printf("MASTER[%d]: enter waitForAll\n", 0);
  //printf("activeRequests is %d\n", activeRequests);
  int row, col;
  int source;
  int i;
  double duration;
  MPI_Status status;
  //MPI_Waitall(np, recvRqsts, &status);
  
  for(i=0; i<activeRequests; i++){
    //printf("in iteration %d\n", i);
    MPI_Recv(&recvLiner[0], blkSize*blkSize+3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    row = (int) recvLiner[blkSize*blkSize];
    col = (int) recvLiner[blkSize*blkSize+1];
    duration = recvLiner[blkSize*blkSize+2];
    add(C[row][col], recvLiner, blkSize);

    source = status.MPI_SOURCE;
    blksTime[source] += duration;
    
  /* debug print
  int i;
  printf("MASTER Recv matrix\n");
  for(i=0; i<blkSize*blkSize; i++){
    printf("%f ", recvLiner[i]);
  }
  printf("\n");*/

  }
  //printf("MASTER[%d]: exit waitForAll\n", 0);
}


int waitForOne(double**** C, MPI_Request* recvRqsts, int np, double* recvLiner, int blkSize, double* blksTime){
  //printf("MASTER[%d]: enter waitForOne\n", 0);
  int row, col;
  int source;
  double duration;
  MPI_Status status;
  //MPI_Waitall(np, recvRqsts, &status);

  MPI_Recv(&recvLiner[0], blkSize*blkSize+3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

  /* debug print
  int i;
  printf("MASTER Recv matrix\n");
  for(i=0; i<blkSize*blkSize; i++){
    printf("%f ", recvLiner[i]);
  }
  printf("\n");*/

  row = (int) recvLiner[blkSize*blkSize];
  col = (int) recvLiner[blkSize*blkSize+1];
  duration = recvLiner[blkSize*blkSize+2];
  

  add(C[row][col], recvLiner, blkSize);
  source = status.MPI_SOURCE;
  blksTime[source] += duration;

  //printf("MASTER[%d]: exit waitForOne, source:%d\n", 0, source);
  return source; 
}



void add(double** C, double *liner, int size){
  int i,j;
  for(i=0; i<size; i++){
    for(j=0; j<size; j++){
      C[i][j] += liner[i*size+j];
    }
  }
}


void fill_m(double**A, double**B, double*liner, int size, int row, int col){
  int i,j;
  for(i=0; i<size; i++){
    for(j=0; j<size; j++){
      liner[i*size+j] = A[i][j];
    }
  }
  for(i=0; i<size; i++){
    for(j=0; j<size; j++){
      liner[size*size + i*size+j] = B[i][j];
    }
  }
  liner[size*size*2] = (double)row;
  liner[size*size*2+1] = (double)col;
}




void zero(double * matrixLiner, int size){
  int i;
  for(i=0; i<size*size; i++){
    matrixLiner[i] = 0;
  }
}

/* usage
 * ***/
void usage(void)
{
	fprintf (stderr, "usage: mm <N> <trials>\n");
	exit(-1);
}

/******************************************************************************/
double cpu_time ( void )
/******************************************************************************/
/*
Purpose:
	CPU_TIME returns the current reading on the CPU clock.
Licensing:
	This code is distributed under the GNU LGPL license. 
Modified:
	06 June 2005
Author:
	John Burkardt
Parameters:
	Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
*/
{
	double value;
	value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;
	return value;
}
/* vim: set ts=4
*/
