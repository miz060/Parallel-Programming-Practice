/**
 * File name: misc.c
 * File description: define helper method for mmmw to reemove redundancy
 * Name:  Mingcheng Zhu
 * Email: zhumc11@gmail.com
 * Date:  Feb 10, 2018
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include "misc.h"

/* fill_w
 * purpose: called by worker to compress input matrix temp into a matrix liner
 * param: input matrix temp, output liner, matrix size, and its row&col in master,
 *        worker computation time duration
 */
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


/* unload_w
 * purpose: called by worker to unload input liner into two matrix 
 * param: input liner, output matrix A and B, matrix size
 */
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
void mult(double **A, double **B, double **C, int N, int M, int P){
  int i, row, col;
  double sum;

  for (row = 0; row < N; row ++){
    for (col = 0; col < P; col ++){
      sum = 0.0;
      for (i = 0; i < M ; i++){
        sum += A[row][i] * B[i][col];
      }
      C[row][col] = sum;
    }
  }
}


/* multM
 * purpose: do block matrix multiplication of AxB, save result to C
 * param: input matrix A,B, output matrix C, matrix size N, M, P, blkSize, 
 *        number of processes np, arrary to record worker's task#,
 *        array to record worker's computation time
 */
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

        // start phase, handout mult task
        if (startPhase >0){
          nextWorker = startPhase --;
        }
        else{
          nextWorker = waitForOne(C, recvRqsts, np, recvLiner, blkSize, blksTime);
          activeRequests --;
        }
        
        // send out task
        fill_m(A[row][i], B[i][col], sendLiner, blkSize, row, col);
        sendToWorker(nextWorker, sendLiner, recvRqsts, blkSize, i);
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
    MPI_Isend(&sendLiner[0], blkSize*blkSize*2+2, MPI_DOUBLE, i, DONE, MPI_COMM_WORLD, &recvRqsts[i]);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // free
  free(temp);
  free(tempLiner);
  free(sendLiner);
  free(recvLiner);
}


/* sendToWorker
 * purpose: called by master to send mult liner to workers
 */
void sendToWorker(int nextWorker, double* sendLiner, MPI_Request * recvRqsts, int blkSize, int i){
  MPI_Send(&sendLiner[0], blkSize*blkSize*2+2, MPI_DOUBLE, nextWorker, DATA, MPI_COMM_WORLD);
}


/* waitForAll
 * purpose: called by master to wait for all workers to finish its work 
 *          and add their results to the master matrix
 */
void waitForAll(double**** C, MPI_Request* recvRqsts, int activeRequests, int np, int blkSize, double*recvLiner, double* blksTime){
  int row, col;
  int source;
  int i;
  double duration;
  MPI_Status status;
  
  for(i=0; i<activeRequests; i++){
    MPI_Recv(&recvLiner[0], blkSize*blkSize+3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    row = (int) recvLiner[blkSize*blkSize];
    col = (int) recvLiner[blkSize*blkSize+1];
    duration = recvLiner[blkSize*blkSize+2];
    // add to the master matrix
    add(C[row][col], recvLiner, blkSize);

    source = status.MPI_SOURCE;
    blksTime[source] += duration;
  }
}


/* waitForOne
 * purpose: called by master to wait for a worker to finish work 
 *          and add its result to the master matrix
 */
int waitForOne(double**** C, MPI_Request* recvRqsts, int np, double* recvLiner, int blkSize, double* blksTime){
  int row, col;
  int source;
  double duration;
  MPI_Status status;

  MPI_Recv(&recvLiner[0], blkSize*blkSize+3, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  row = (int) recvLiner[blkSize*blkSize];
  col = (int) recvLiner[blkSize*blkSize+1];
  duration = recvLiner[blkSize*blkSize+2];
  
  // add to the master matrix
  add(C[row][col], recvLiner, blkSize);
  source = status.MPI_SOURCE;
  blksTime[source] += duration;

  return source; 
}


/* add
 * purpose: do matrix addition of matrix liner B on A
 * param: added matrix A, input matrix liner B, matrix size
 */
void add(double** C, double *liner, int size){
  int i,j;
  for(i=0; i<size; i++){
    for(j=0; j<size; j++){
      C[i][j] += liner[i*size+j];
    }
  }
}


/* fill_m
 * purpose: called by master to compress input matrix A and B into a matrix liner
 * param: input matrix A and B, output liner, matrix size, row and col in master matrix
 */
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


/* zero
 * purpose: zero fill input matrix
 * param: matrix to zero fill, matrix size
 */
void zero(double * matrixLiner, int size){
  int i;
  for(i=0; i<size*size; i++){
    matrixLiner[i] = 0;
  }
}


/* usage
 * purpose: print out the usage message to stdout
 */
void usage(void){
  fprintf (stderr, "usage: mmmw <blocksize> <N> <M> <P> <output file>\n");
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
