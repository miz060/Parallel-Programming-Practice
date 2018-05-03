/** 
 * CSE 160: Programming Assignment 2
 * File name: misc.c
 * File description: define helper method for heat2d programs to remove redundancy.
 * Name:  Mingcheng Zhu
 * PID:   A92047564
 * Email: miz060@ucsd.edu
 * Date:  Jan 28, 2018
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include "misc.h"
#include "svalidate.h"


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

/* usage
 * purpose: print usage info of heat2d to stderr
 */
int usage(){
	fprintf(stderr, "usage: heat2d M N Tl Tr Tt Tb eps file\n");
}


/* cLineValid
 * purpose: validate command line arguments for the calling program
 * param:   argv1,2,3,4,5,6,7 - command line arguments of the calling program
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
int cLineValid(char* argv1, char* argv2, char* argv3, 
               char* argv4, char* argv5, char* argv6,
               char* argv7){
  
  // command line validation
  if(!isInteger(trim(argv1)) || !isInteger(trim(argv2))){
    return FALSE;
  }
  int M = atoi(argv1);
  int N = atoi(argv2);
  if(M<10 || N<10) return FALSE;

  if(!isDouble(trim(argv3)) || !isDouble(trim(argv4)) ||
     !isDouble(trim(argv5)) || !isDouble(trim(argv6))){
    return FALSE;
  }

  if(!isFloat(trim(argv7))) return FALSE;
  double eps = atof(argv7);
  if(eps>=1.0 || eps<=0) return FALSE;

  return TRUE;
}


/* serial
 * purpose: do serial computation of the convergence iterations
 * param:   M, N, u, w, mean, eps, Tl, Tr, Tt, Tb - the relative variables
 *          in the calling program
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
int serial(int M, int N, double** u, double** w, double mean, double eps, 
           double Tl, double Tr, double Tt, double Tb){

  // iterate until the  new solution W differs from the old solution U
  // by no more than EPSILON.
  double global_diff = eps;
  int iterations = 0 ; 
  int iterations_print = 1;
  int i;
  int j;
  double ctime;
  double ctime1;
  double ctime2;
  // set up the inital matrix
  for ( i = 1; i < M - 1; i++ ){
    w[i][0] = Tl;
    w[i][N-1] = Tr;
  }
  for ( j = 0; j < N; j++ ){
    w[0][j] = Tt;
    w[M-1][j] = Tb;
  }
  for ( i = 1; i < M - 1; i++ ){
    for ( j = 1; j < N - 1; j++ ){
      w[i][j] = mean;
    }
  }

  ctime1 = cpu_time ( );
  while ( eps <= global_diff){
    // save the old solution in U.
    for ( i = 0; i < M; i++){
      for ( j = 0; j < N; j++ ){
        u[i][j] = w[i][j];
      }
    }

    // Determine the new estimate of the solution at the interior points.
    // The new solution W is the average of north, south, east and west neighbors.
    global_diff = 0.0;
    for ( i = 1; i < M - 1; i++ ){
      for ( j = 1; j < N - 1; j++ ){
        w[i][j] = ( u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] ) / 4.0;

        if ( global_diff < fabs ( w[i][j] - u[i][j] ) ){
          global_diff = fabs ( w[i][j] - u[i][j] );
        }
      }
    }
    // print the iteration info
    iterations++;
    if ( iterations == iterations_print ){
      printf ( "  %8d  %f\n", iterations, global_diff);
      iterations_print = 2 * iterations_print;
    }
  }
  // get the calculation time
  ctime2 = cpu_time ( );
  ctime = ctime2 - ctime1;

  printf ( "\n" );
  printf ( "  %8d  %f\n", iterations, global_diff );
  printf ( "\n" );
  printf ( "  Error tolerance achieved.\n" );
  printf ( "  CPU time = %f\n", ctime );
  return 0;
}


/* serial_write
 * purpose: do serial write of a complete matrix to the outfile
 * param:   outfile, M, N, w, fp - the matching variables in the calling program
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
int serial_write(char* output_file, int M, int N, double**w, FILE* fp){

  int i;
  int j;
  // print M and N value
  fp = fopen ( output_file, "w" );
  fprintf ( fp, "%d\n", M );
  fprintf ( fp, "%d\n", N );
  // print the data matirx
  for ( i = 0; i < M; i++ ){
    for ( j = 0; j < N; j++){
      fprintf ( fp, "%6.2f ", w[i][j] );
    }
    fputc ( '\n', fp);
  }
  fclose ( fp );
  printf ( "\n" );
  printf ("  Solution written to the output file '%s'\n", output_file );
  printf ( "\n" );

  return 0;
}


/* serial_write_s
 * purpose: do serial write of a complete matrix to the outfile
 *          the same as serial_write, except w would be a char* matrix
 *          (this one is special for heat2dcol.c
 * param:   outfile, M, N, w, fp - the matching variables in the calling program
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
int serial_write_s(char* output_file, int M, int N, char***w, FILE* fp){

  int i;
  int j;
  fp = fopen ( output_file, "w" );
  // print M and N value
  fprintf ( fp, "%d\n", M );
  fprintf ( fp, "%d\n", N );
  // print the data matirx
  for ( i = 0; i < M; i++ ){
    for ( j = 0; j < N; j++){
      fprintf ( fp, "%6s ", w[i][j] );
    }
    fputc ( '\n', fp);
  }
  fclose ( fp );
  printf ( "\n" );
  printf ("  Solution written to the output file '%s'\n", output_file );
  printf ( "\n" );

  return 0;
}


/* calculation
 * purpose: do serial calculation of the matrix convergence
 * param:   u, w, my_rank, comm_sz, length, top, bottom - the matching variables in the 
 *          calling program
 *          X - the array length of the partial matrix. N for heat2drow, M for heat2dcol
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
double calculation(double** u, double**w, int X, int my_rank,
                   int comm_sz, int length, double* top, double* bottom){
  double global_diff = 0.0;
  int i;
  int j;
  // inner unit calculation
  for(i=1; i<length-1; i++){
    for(j=1; j<X-1; j++){
      w[i][j] = ( u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] ) / 4.0;

      if (global_diff < fabs(w[i][j] - u[i][j])){
        global_diff = fabs ( w[i][j] - u[i][j] );
      }
    }
  }
  // each process has single row case
  if(my_rank!=0 && my_rank!= comm_sz-1 && length == 1){
    for(j=1; j<X-1; j++){
      w[0][j] = (top[j] + bottom[j] + u[0][j-1] + u[0][j+1]) / 4.0;
      if (global_diff < fabs(w[0][j] - u[0][j])){
        global_diff = fabs ( w[0][j] - u[0][j] );
      }
    }
  }
  // mutilple rows per process, non-process 0 case
  // avoid w[-1][j]
  if(my_rank!=0 && length!=1){
    for(j=1; j<X-1; j++){
      w[0][j] = (top[j] + u[1][j] + u[0][j-1] + u[0][j+1]) / 4.0;
      if (global_diff < fabs(w[0][j] - u[0][j])){
        global_diff = fabs ( w[0][j] - u[0][j] );
      }
    }
  }
  // mutilple rows per process, non-process comm_sz-1 case
  // avoid w[comm_sz][j]
  if(my_rank!=(comm_sz-1)&&length != 1){
    for(j=1; j<X-1; j++){
      w[length-1][j] = (u[length-2][j] + bottom[j] + u[length-1][j-1] + u[length-1][j+1]) /4.0;
      if (global_diff < fabs(w[length-1][j] - u[length-1][j])){
        global_diff = fabs ( w[length-1][j] - u[length-1][j] );
      }
    }
  }
  return global_diff;
}





