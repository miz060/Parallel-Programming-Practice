/** 
 * CSE 160: Programming Assignment 2
 * File name: heat2d.c
 * File description: serial version of the heat equation program
 * Name:  Mingcheng Zhu
 * PID:   A92047564
 * Email: miz060@ucsd.edu
 * Date:  Jan 28, 2018
 */
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include "svalidate.h"
# include "misc.h"

int main ( int argc, char *argv[] );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HEATED_PLATE.

  Discussion:

    This code solves the steady state heat equation on a rectangular region.

    The physical region, and the boundary conditions, are suggested
    by this diagram;

                   W = 0
             +------------------+
             |                  |
    W = Tl   |                  | W = Tr 
             |                  |
             +------------------+
                   W = 0

    The region is covered with a grid of M by N nodes, and an N by N
    array W is used to record the temperature.  The correspondence between
    array indices and locations in the region is suggested by giving the
    indices of the four corners:

                  I = 0
          [0][0]-------------[0][N-1]
             |                  |
      J = 0  |                  |  J = N-1
             |                  |
        [M-1][0]-----------[M-1][N-1]
                  I = M-1

    The steady state solution to the discrete heat equation satisfies the
    following condition at an interior grid point:

      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )

    where "Central" is the index of the grid point, "North" is the index
    of its immediate neighbor to the "north", and so on.
   
    Given an approximate solution of the steady state heat equation, a
    "better" solution is given by replacing each interior point by the
    average of its 4 neighbors - in other words, by using the condition
    as an ASSIGNMENT statement:

      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )

    If this process is repeated often enough, the difference between successive 
    estimates of the solution will go to zero.

    This program carries out such an iteration, using a tolerance specified by
    the user, and writes the final estimate of the solution to a file that can
    be used for graphic processing.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    22 July 2008

  Author:

    Original C version by Michael Quinn.
    Modifications by John Burkardt.


    More modifications by Philip Papadopoulos

  Reference:

    Michael Quinn,
    Parallel Programming in C with MPI and OpenMP,
    McGraw-Hill, 2004,
    ISBN13: 978-0071232654,
    LC: QA76.73.C15.Q55.

  Parameters:

    Commandline argument 1,  M  number of rows
    Commandline argument 2,  N  number of columns 
    Commandline argument 3, double Tleft, Temperature along the left boundary 
    Commandline argument 4, double Tright, Temperature along the right boundary.  
    Commandline argument 5, double Ttop, Temperature along the top boundary.  
    Commandline argument 6, double Tbottom, Temperature along the bottom boundary.  
    Commandline argument 7, double EPSILON, the error tolerance.  

    Commandline argument 8, char *OUTPUT_FILE, the name of the file into which
    the steady state solution is written when the program has completed.

  Local parameters:

    Local, double DIFF, the norm of the change in the solution from one iteration
    to the next.

    Local, double MEAN, the average of the boundary values, used to initialize
    the values of the solution in the interior.

    Local, double U[M][N], the solution at the previous iteration.

    Local, double W[M][N], the solution computed at the latest iteration.
*/
{
  double eps;
  FILE *fp=NULL;
  int M;
  int N;
  int i;
  int j;
  double mean;
  char *output_file;
  double **u;
  double **w;
  double Tl,Tr,Tt,Tb;

  // command line validation
  if (argc != 9){ 
    usage();
    return -1;
  }
  if(!cLineValid(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7])){
    usage();
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

  printf ( "HEAT2D\n" );
  printf ( "  C version\n" );
  printf ( "  A program to solve for the steady state temperature distribution\n" );
  printf ( "  over a rectangular plate.\n" );
  printf ( "  Spatial grid of %d by %d points.\n", M, N );
  printf ( "\n" );

  // allocate memory for the partial matrix
  u = (double **) malloc(M*sizeof(double *));
  w = malloc(M * sizeof(double *));
  for (i = 0; i < M; i ++) {
	u[i] = (double *) malloc(N  * sizeof(double));
	w[i] = (double *) malloc(N  * sizeof(double));
  }
  /** Note: u[i][j] = *( *(arr +i) + j) **/
  printf ( "\n" );
  printf ( "  The iteration will be repeated until the change is <= %G\n", eps );
  printf ( "  Boundary Temperatures top: %G  bottom: %G  left: %G  right: %G\n", Tt, Tb, Tl, Tr );
  printf ( "  The steady state solution will be written to '%s'.\n", output_file );

  //Set the boundary values, which don't change. 
  for ( i = 1; i < M - 1; i++ )
  {
    w[i][0] = Tl;
    w[i][N-1] = Tr;
  }
  for ( j = 0; j < N; j++ )
  {
    w[0][j] = Tt;
    w[M-1][j] = Tb;
  }
  // Average the boundary values, to come up with a reasonable
  // initial value for the interior.
  mean = 0.0;
  for ( i = 1; i < M - 1; i++ )
  {
    mean += w[i][0];
    mean += w[i][N-1];
  }
  for ( j = 0; j < N; j++ )
  {
    mean += w[0][j];
    mean += w[M-1][j];
  }
  mean = mean / ( double ) ( 2 * M + 2 * N - 4 );

  // Initialize the interior solution to the mean value.
  for ( i = 1; i < M - 1; i++ )
  {
    for ( j = 1; j < N - 1; j++ )
    {
      w[i][j] = mean;
    }
  }

  // iterate until the  new solution W differs from the old solution U
  // by no more than EPSILON.
  printf("\n");
  printf(" Iteration Change\n");
  printf("\n");
  serial(M, N, u, w, mean, eps, Tl, Tr, Tt, Tb);

  // Write the solution to the output file.
  serial_write(output_file, M, N, w, fp);
  printf ( "HEAT2D:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;

}

