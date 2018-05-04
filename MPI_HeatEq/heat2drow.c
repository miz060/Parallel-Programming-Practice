/** 
 * File name: heat2drow.c
 * File description: row decomposition version of the heat equation program
 * Name:  Mingcheng Zhu
 * Email: zhumc11@gmail.com
 * Date:  Jan 28, 2018
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "svalidate.h"
#include "misc.h"
#include <mpi.h>

#define DATA 0
#define PRINT 1

int main ( int argc, char *argv[] );

/* main
 * purpose: main function for the heat equation program
 * param:   argc - cmd argument counts, argv - cmd argument array
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
int main(int argc, char* argv[]){

  double ctime;
  double ctime1;
  double ctime2;
  double global_diff;
  double eps;
  FILE *fp=NULL;
  int M;
  int N;
  int i;
  int iterations;
  int iterations_print;
  int j;
  double mean;
  char *output_file;
  double **u;
  double **w;
  double Tl,Tr,Tt,Tb;

  int comm_sz;
  int my_rank;
  // set up MPI components
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // command line validation
  if (argc != 9){
    if(my_rank == 0){
      usage();
    }
    MPI_Finalize();
    return -1;
  }
    
  if(!cLineValid(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7])){
    if(my_rank == 0){
      usage();
    }
    MPI_Finalize();
    return -1;
  }

  // variable initialization
  M = atoi(argv[1]); 
  N = atoi(argv[2]); 
  Tl = atof(argv[3]);
  Tr = atof(argv[4]);
  Tt = atof(argv[5]);
  Tb = atof(argv[6]);
  eps = atof(argv[7]);
  output_file = argv[8];


  // set up variables for process transmission
  int pre_rank;
  int next_rank;
  pre_rank = (comm_sz+my_rank-1)%comm_sz;
  next_rank = (my_rank+1)%comm_sz;

  // divide the work set
  int length = N/comm_sz;
  int extra = N%comm_sz;
  int other_length = length;
  if(my_rank ==0){
    if(extra != 0){
      length = length+extra;
    }
  }
  // check if easy division is possible
  if (comm_sz > N){
    if(my_rank == 0){
      usage();
    }
    MPI_Finalize();
    return -1;
  }

  if(my_rank == 0){
    printf ( "HEAT2DROW\n" );
    printf ( "  C version\n" );
    printf ( "  A program to solve for the steady state temperature distribution\n" );
    printf ( "  over a rectangular plate.\n" );
    printf ( "  Spatial grid of %d by %d points.\n", N, M );
    printf ( "\n" );
  }

  // allocate memory for this process
  u = (double **) malloc(length*sizeof(double *));
  w = malloc(length* sizeof(double *));
  for (i = 0; i < length; i ++) {
	u[i] = (double *) malloc(M  * sizeof(double));
	w[i] = (double *) malloc(M  * sizeof(double));
  }

  // print user info as process 0
  if(my_rank == 0){
    printf ( "\n" );
    printf ( "  The iteration will be repeated until the change is <= %G\n", eps );
    printf ( "  Boundary Temperatures top: %G  bottom: %G  left: %G  right: %G\n", Tt, Tb, Tl, Tr );
    printf ( "  The steady state solution will be written to '%s'.\n", output_file );
  }
  global_diff = eps;
  // set interior solution to be mean of the boundary value
  mean = (Tt*N + Tb*N + Tl*(M-2) + Tr*(M-2)) / (double)(2*N + 2*M-4);

  // set the top and down boundry values as boundry process
  // hanlde process 0 top row
  if(my_rank == 0){
    for(j=0; j<M; j++){
      w[0][j] = Tt;
    }
    if(length!=1){
      for(j=1; j<M-1; j++){
        w[length-1][j] = mean;
      }
      w[length-1][0] = Tl;
      w[length-1][M-1] = Tr;
    }
  }
  // handle process comm_sz-1 bottom row
  else if(my_rank == comm_sz-1){
    for(j=0; j<M; j++){
      w[length-1][j] = Tb;
    }
    if(length!=1){
      for(j=1; j<M-1; j++){
        w[0][j] = mean;
      }
      w[0][0] = Tl;
      w[0][M-1] = Tr;
    }
  }
  else{
    for(j=0; j<M; j++){
      w[0][j] = mean;
      w[length-1][j] = mean;
    }
    w[0][0] = Tl;
    w[0][M-1] = Tr;
    w[length-1][0] = Tl;
    w[length-1][M-1] = Tr;
  }

  // initialize the interior solution to be the mean and also set up boundry value for Tl and Tr
  for(i=1; i<length-1; i++){
    w[i][0] = Tl;
    w[i][M-1] = Tr;
  } 
  for(i=1; i<length-1; i++){
    for(j=1; j<M-1; j++){
      w[i][j] = mean;
    }
  }

  MPI_Status status[4];
  MPI_Request request[4];
  MPI_Request creq[2];
  MPI_Status cstats[2];

  // iterate until the new solution differs from the old solution U no more tha EPSILON.
  iterations = 0;
  iterations_print = 1;

  // time the calculation as process 0
  if(my_rank == 0){
    printf("\n");
    printf(" Iteration Change\n");
    printf("\n");
    ctime1 = cpu_time();
  }

  // single processor case: just behave like a serial program
  if(comm_sz == 1){
    serial(N, M, u, w, mean, eps, Tl, Tr, Tt, Tb);

    // write the solution to the output file.
    serial_write(output_file, N, M, w, fp);
    printf ( "HEAT2DROW:\n" );
    printf ( "  Normal end of execution.\n" );
    MPI_Finalize();
    return 0;
  }

  // actual parallel case 
  while(eps <= global_diff){
    // save the old solution in u.
    for(i=0; i<length; i++){
      for(j=0; j<M; j++){
        u[i][j] = w[i][j];
      }
    }
    // send lowest row to next process
    double* bottom = malloc(M*sizeof(double));
    double* top = malloc(M*sizeof(double));
    if(my_rank!=comm_sz-1 && my_rank!=0){
      MPI_Isend(&u[length-1][0], M, MPI_DOUBLE, next_rank, iterations, MPI_COMM_WORLD, &request[0]);
      MPI_Irecv(&bottom[0], M, MPI_DOUBLE, next_rank, iterations, MPI_COMM_WORLD, &request[1]);

      MPI_Isend(&u[0][0], M, MPI_DOUBLE, pre_rank, iterations, MPI_COMM_WORLD, &request[2]);
      MPI_Irecv(&top[0], M, MPI_DOUBLE, pre_rank, iterations, MPI_COMM_WORLD, &request[3]);
      // wait for transmission to finish
      MPI_Waitall(4, request, status);
    }
    else if(my_rank==0){
      MPI_Isend(&u[length-1][0], M, MPI_DOUBLE, next_rank, iterations, MPI_COMM_WORLD, &creq[0]);
      MPI_Irecv(&bottom[0], M, MPI_DOUBLE, next_rank, iterations, MPI_COMM_WORLD, &creq[1]);
      // wait for transmission to finish
      MPI_Waitall(2, creq, cstats);
    }
    else if(my_rank==comm_sz-1){
      MPI_Isend(&u[0][0], M, MPI_DOUBLE, pre_rank, iterations, MPI_COMM_WORLD, &creq[0]);
      MPI_Irecv(&top[0], M, MPI_DOUBLE, pre_rank, iterations, MPI_COMM_WORLD, &creq[1]);
      // wait for transmission to finish
      MPI_Waitall(2, creq, cstats);
    } 
  
    // do the calculation
    global_diff = calculation(u, w, M, my_rank, comm_sz, length, top, bottom);
    MPI_Allreduce(&global_diff, &global_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    iterations++;
    if( iterations == iterations_print){
      if(my_rank == 0){
        printf ( "  %8d  %f\n", iterations, global_diff);
        iterations_print = 2 * iterations_print;
      }
    }
    free(bottom);
    free(top);
  }
  // print execution info as process 0
  if(my_rank == 0){
    ctime2 = cpu_time();
    ctime = ctime2 - ctime1;
    printf ( "\n" );
    printf ( "  %8d  %f\n", iterations, global_diff );
    printf ( "\n" );
    printf ( "  Error tolerance achieved.\n" );
    printf ( "  CPU time = %f\n", ctime );
  }
      

  MPI_Status p_status; 
  // send the partial matrix to process zero
  if(my_rank!=0){
    char* str = (char*) malloc(8*M*length*sizeof(char));
    strcpy(str, "");
    for(i=0; i<length; i++){
      for(j=0; j<M; j++){
        char temp[8];
        sprintf(temp, "%6.2f ", w[i][j]);
        strcat(str, temp);
      }
      strcat(str, "\n");
    }  
    MPI_Send(str, 8*M*length, MPI_CHAR, 0, PRINT, MPI_COMM_WORLD);
    free(str);
  }
  // collect and print the whole matrix as process zero
  else{
    fp = fopen(output_file, "w");
    fprintf(fp, "%d\n", N);
    fprintf(fp, "%d\n", M);

    // process 0 self printing
    for ( i = 0; i < length; i++ ){
      for ( j = 0; j < M; j++){
        fprintf ( fp, "%6.2f ", w[i][j] );
      }
      fputc ( '\n', fp);
    }
    // process 0 receive partial matrix from others
    int id;
    for(id=1; id<comm_sz; id++){
      char* part = (char*) malloc(8*M*length*sizeof(char));
      MPI_Recv(part, 8*M*other_length, MPI_CHAR, id, PRINT, MPI_COMM_WORLD, &p_status);

      // print the partial table to the file
      fprintf(fp, part);
      free(part);
    }
    fclose(fp); 
    printf ( "\n" );
    printf ("  Solution written to the output file '%s'\n", output_file );

    printf ( "\n" );
    printf ( "HEAT2DROW:\n" );
    printf ( "  Normal end of execution.\n" );
  }

  // clean up
  MPI_Barrier(MPI_COMM_WORLD);
  for (i = 0; i < length; i ++) {
	  free(u[i]);
	  free(w[i]);
  }
  free(u);
  free(w);
  MPI_Finalize();
  return 0;
}
    
    













