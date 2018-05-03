/**
 * CSE 160: Programming Assignment 3
 * File name: mmmw.c
 * File description: block matrix multiplication with master-worker algo
 * Name:  Mingcheng Zhu
 * PID:   A92047564
 * Email: miz060@ucsd.edu
 * Date:  Feb 10, 2018
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include "misc.h"

#define DATA 10
#define DONE 0

/* 
 * Usage:  mmmw <blocksize> <N> <M> <P> <output file> 
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

  // set up MPI components
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // command line validation
  if (argc != 6){
    if(my_rank == 0){
      usage();
    }
    MPI_Finalize();
    return -1;
  }

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
  
    // zero fill matrix C
    double * blkLiner;
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
  
    // print A and B matrix to outfile
    int n,m;
    int ni,mj;
    double** pblock;
    fprintf(fp, "%d,%d,%d\n", N*blkSize, M*blkSize, P*blkSize);
    for (i=0; i<N*blkSize; i++){
      for (j=0; j<M*blkSize-1; j++){
        n = i/blkSize;
        ni = i%blkSize;
        m = j/blkSize;
        mj = j%blkSize;
        pblock = A[n][m];
        fprintf(fp, "%15.8f,", pblock[ni][mj]);
      }
      n = i/blkSize;
      ni = i%blkSize;
      m = j/blkSize;
      mj= blkSize-1;
      pblock = A[n][m];
      fprintf(fp, "%15.8f", pblock[ni][mj]);
      fprintf(fp, "\n");
    }
    for (i=0; i<M*blkSize; i++){
      for (j=0; j<P*blkSize-1; j++){
        n = i/blkSize;
        ni = i%blkSize;
        m = j/blkSize;
        mj = j%blkSize;
        pblock = B[n][m];
        fprintf(fp, "%15.8f,", pblock[ni][mj]);
      }
      n = i/blkSize;
      ni = i%blkSize;
      m = j/blkSize;
      mj= blkSize-1;
      pblock = B[n][m];
      fprintf(fp, "%15.8f", pblock[ni][mj]);
      fprintf(fp, "\n");
    }

    // time the computation
    startTime = cpu_time();
    multM(A,B,C,N,M,P,blkSize, comm_sz, blksPerRank, blksTime);
    endTime = cpu_time();

    // print C matrix to outfile
    for (i=0; i<N*blkSize; i++){
      for (j=0; j<P*blkSize-1; j++){
        n = i/blkSize;
        ni = i%blkSize;
        m = j/blkSize;
        mj = j%blkSize;
        pblock = C[n][m];
        fprintf(fp, "%15.8f,", pblock[ni][mj]);
      }
      n = i/blkSize;
      ni = i%blkSize;
      m = j/blkSize;
      mj= blkSize-1;
      pblock = C[n][m];
      fprintf(fp, "%15.8f", pblock[ni][mj]);
      fprintf(fp, "\n");
    }

    // print stdout message
    printf("Parameters: %d, %d, %d, %d, %d\n", blkSize, N, M, P, comm_sz-1);
    for(i=1; i<comm_sz; i++){
      printf("Worker %d:  %d, %15.8f\n", i, blksPerRank[i], blksTime[i]);
    }
    printf("Total Time: %15.8f\n", endTime-startTime);

    // clean up
    fclose(fp);
    for (i=0; i<N*P; i++){
      free(Clinear[i][0]);
      free(Clinear[i]);
    }
    for (i=0; i<N*M; i++){
      free(Alinear[i][0]);
      free(Alinear[i]);
    }
    for (i=0; i<M*P; i++){
      free(Blinear[i][0]);
      free(Blinear[i]);
    }
    free(A);
    free(Alinear);
    free(B);
    free(Blinear);
    free(C);
    free(Clinear);
    free(blksPerRank);
    free(blksTime);
  }

  // Worker set
  else{
    MPI_Status status;
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
      MPI_Recv(&recvLiner[0], blkSize*blkSize*2+2, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      tag = status.MPI_TAG;
      // receive the done message
      if(tag == DONE){
        done = 1;
      }
      // receive a task
      else{
        row = recvLiner[blkSize*blkSize*2];
        col = recvLiner[blkSize*blkSize*2+1];

        unload_w(recvLiner, X, Y, blkSize);
        startTime = cpu_time();
        mult(X, Y, Z, blkSize, blkSize, blkSize);
        endTime = cpu_time();
        double duration = endTime - startTime;

        fill_w(Z, sendLiner, blkSize, row, col, duration);
        MPI_Send(&sendLiner[0], blkSize*blkSize+3, MPI_DOUBLE, 0, DATA, MPI_COMM_WORLD);
        index++;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(i=0; i<blkSize; i++){
      free(X[i]);
      free(Y[i]);
      free(Z[i]);
    }
    free(X);
    free(Y);
    free(Z);
    free(recvLiner);
    free(sendLiner);
  }
  MPI_Finalize();
  return 0;
}




