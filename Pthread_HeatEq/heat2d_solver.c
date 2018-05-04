#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "heat2d_solver.h"

/** 
 * CSE 160: Programming Assignment 4
 * File name: heat2d_solver.c
 * File description: defines the solver function for heat2d.c.
 *                   use row decomposition to do the calculation.
 *                   parallelism achieved through pthread.
 * Name:  Mingcheng Zhu
 * PID:   A92047564
 * Email: miz060@ucsd.edu
 * Date:  Feb 25, 2018
 */

// global variables
int M;
int N;
double eps;
int print;
volatile double **u;
double *tol;
int thread_count;
int iterations;
int iterations_print;
int i,j;
double diff;


// mutex to guard diff
pthread_mutex_t* mutex_diff;
// barrier to synchronize threads within each iteration
pthread_barrier_t* barrier;


/* calculation
 *  rank - the rank of the current thread
 *  purpose: the thread function for head2dSolve. 
 *    distribute the matrix u to different threads and
 *    synchronize them to collaborately solve the heat 
 *    distribution problem.
 */
// thread function
void *calculation (void* rank){
  long my_rank = (long) rank;
  int height = M/thread_count;
  int extra = M%thread_count;

  int first_row = my_rank*height;
  int last_row = (my_rank+1)*height - 1;
  if(my_rank == thread_count - 1){
    last_row = last_row + extra;
  }
  double local_diff;
  int local_iterations = 0;

  double *rowPrev = calloc(N, sizeof(double));
  double *rowCurr = calloc(N, sizeof(double));
  double *rowTmp;

  double *up_neignbor = calloc(N, sizeof(double));
  double *down_neighbor = calloc(N, sizeof(double));

  // threaded iteration
  while ( eps <= diff ){
    int i,j;

    // fetching data from blocks of other threads
    if(my_rank == 0 && thread_count != 1){
      memcpy(down_neighbor, (double*)u[last_row+1], N*sizeof(double));
    }
    else if(my_rank == thread_count-1 && thread_count != 1){
      memcpy(up_neignbor, (double*)u[first_row-1], N*sizeof(double));
    }
    else if(thread_count != 1){
      memcpy(up_neignbor, (double*)u[first_row-1], N*sizeof(double));
      memcpy(down_neighbor, (double*)u[last_row+1], N*sizeof(double));
    }

    local_diff = 0.0;
    // wait till all threads finish
    pthread_barrier_wait(barrier);

    //Initialize the first row based on local data and reset global diff
    memcpy(rowCurr, (double*)u[first_row], N*sizeof(double));
    pthread_mutex_lock(mutex_diff);
    diff = 0.0;
    pthread_mutex_unlock(mutex_diff);
    pthread_barrier_wait(barrier);
    // wait till all threads finish

    /*
       Determine the new estimate of the solution at the interior points.
       The new solution W is the average of north, south, east and west 
       neighbors.  */

    // mutilple rows per thread case
    if(height != 1 || ((my_rank == thread_count-1) && extra>0)){
      // update the first row
      if(my_rank !=0){
        for ( j = 1; j < N - 1; j++ ){
          u[first_row][j] = (up_neignbor[j] + u[first_row+1][j] + 
              rowCurr[j-1] + rowCurr[j+1] ) / 4.0;

          double delta = fabs(rowCurr[j] - u[first_row][j]);
          if ( local_diff < delta ) {
            local_diff = delta; 
          }
        }
      }

      // update and calculate local diff of the inner rows
      for ( i = first_row+1; i < last_row; i++ ){
        /* swap rowPrev and rowCurr pointers. Save the current row */
        rowTmp = rowPrev; rowPrev=rowCurr; rowCurr=rowTmp;
        memcpy(rowCurr, (double*)u[i], N*sizeof(double));

        for ( j = 1; j < N - 1; j++ ){
          u[i][j] = (rowPrev[j] + u[i+1][j] + 
              rowCurr[j-1] + rowCurr[j+1] ) / 4.0;

          double delta = fabs(rowCurr[j] - u[i][j]);
          if ( local_diff < delta ){
            local_diff = delta; 
          }
        }
      }

      // update the bottom row
      if(my_rank != thread_count-1){
        rowTmp = rowPrev; rowPrev=rowCurr; rowCurr=rowTmp;
        memcpy(rowCurr, (double*)u[last_row], N*sizeof(double));
        for ( j = 1; j < N - 1; j++ ){
          u[last_row][j] = (rowPrev[j] + down_neighbor[j] + 
              rowCurr[j-1] + rowCurr[j+1] ) / 4.0;

          double delta = fabs(rowCurr[j] - u[last_row][j]);
          if ( local_diff < delta ){
            local_diff = delta; 
          }
        }
      }
    }

    // single row per thread case
    else{
      if(my_rank != 0 && my_rank != thread_count-1){
        for(j=1; j<N-1; j++){
          u[first_row][j] = (up_neignbor[j] + down_neighbor[j] + rowCurr[j-1] + rowCurr[j+1] ) / 4.0;

          double delta = fabs(rowCurr[j] - u[last_row][j]);
          if ( local_diff < delta ){
            local_diff = delta; 
          }
        }
      }
    }

    // calculate the global diff
    pthread_mutex_lock(mutex_diff);
    if(diff < local_diff){
      diff = local_diff;
    }
    pthread_mutex_unlock(mutex_diff);

    // wait till all thread finishes
    pthread_barrier_wait(barrier);

    local_iterations++;

    // thread 0 do the printing
    if(my_rank == 0){
      if ( print && local_iterations == iterations_print ){
        printf ( "  %8d  %f\n", local_iterations, diff );
        iterations_print *= 2;
      }
    }
  }

  if(my_rank == 0){
    iterations = local_iterations;
  }

  // clean up
  free(rowCurr);
  free(rowPrev);
  free(up_neignbor);
  free(down_neighbor);
}


/* heat2dSolve 
 * 	M - number of rows (input)
 * 	N - number of cols (output)
 * 	u - temperature distribution(input/output)
 *	eps - tolerance
 *	print - print iteration information (boolean)
 *
 * 	returns
 * 	    - number of iterations
 * 	    - u contains the final temperature distribution 
 */
int heat2dSolve(int M_, int N_, double eps_, int print_, volatile double **u_, double *tol_, pthread_t* thread_handles, int thread_count_){

  // set up the global variale based on input parameters
  M = M_;
  N = N_;
  eps = eps_;
  print = print_;
  u = u_;
  tol = tol_;
  thread_count = thread_count_;

  iterations_print = 1;
  diff = 2.0 * eps;

  if (print) 
    printf( "\n Iteration  Change\n" );

  // initialize the lock and barrier
  mutex_diff = malloc(sizeof(pthread_mutex_t));
  barrier = malloc(sizeof(pthread_barrier_t));
  pthread_mutex_init(mutex_diff, NULL);
  pthread_barrier_init(barrier, NULL, thread_count);


  // start threading
  long thread;
  for(thread = 0; thread < thread_count; thread++){

    pthread_create(&thread_handles[thread], NULL, calculation, (void*) thread);
  }

  // join all threads
  for(thread = 0; thread < thread_count; thread++){
    pthread_join(thread_handles[thread], NULL);
  }

  /* memory cleanup */
  pthread_mutex_destroy(mutex_diff);
  pthread_barrier_destroy(barrier);
  free(mutex_diff);
  free(barrier);
  *tol = diff;
  return iterations;
}
