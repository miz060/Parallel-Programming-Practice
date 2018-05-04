/* mm.c 
 * Sample matrix multiplication in C 
 * Author: Philip Papadopoulos
 * Email: ppapadopoulos@ucsd.edu
 */

/* Modified: Mingcheng Zhu
 * Email: zhumc11@gmail.com
 * Dtae:  Feb 10, 2016
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
double cpu_time ( void );
void mult(double **, double **, double **, int, int, int);
void usage_2( void );
/* ===========================================================================
 * Usage: mm <N> <M> <P> <trials>
 */
int main(int argc, char * argv[]){
	if (argc != 5) usage_2();
	int N=atoi(argv[1]); 
	int M=atoi(argv[2]);
	int P=atoi(argv[3]);
	int rounds=atoi(argv[4]);

	double **A, **B, **C;
	double *Alinear, *Blinear, *Clinear;
  	// A is NxM, B is MxP, C is NxP
	A = (double **) calloc(N, sizeof (double **));
	B = (double **) calloc(M, sizeof (double **));
	C = (double **) calloc(N, sizeof (double **));
	Alinear = calloc(N * M, sizeof (double));
	Blinear = calloc(M * P, sizeof (double));
	Clinear = calloc(N * P, sizeof (double));
	int i;

	/* populate A, B ,C so that we can use  A[i][j] addressing */
	for (i=0; i<N ; i++){ 
		A[i] = Alinear + i * M ;
		C[i] = Clinear + i * P;
	}
	for (i=0; i<M; i++){
		B[i] = Blinear + i * P;
	}

	// fill random values into A and B
	for (i=0; i<N*M; i++){
		Alinear[i] = drand48();
	}
	for (i=0; i<M*P; i++){
		Blinear[i] = drand48();
	}

	// time the computation
	double startTime, endTime;
	startTime = cpu_time();
	for (i = 0; i < rounds; i ++)
	{
		mult(A,B,C,N,M,P);
	}
	endTime = cpu_time(); 
	double avgTime = (endTime - startTime)/(double) rounds;
	printf("%dx%d mult %dx%d, in %f seconds (%d trials) \n", 
		N, M, M, P, avgTime, rounds);

	// clean up
	free(A);
	free(Alinear);
	free(B);
	free(Blinear);
	free(C);
	free(Clinear);
	return 0;
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
			for (i =0 ; i < M ; i++){
				sum += A[row][i] * B[i][col];
			}
			C[row][col] = sum;
		}
	}
}

/* usage
 * ***/
void usage_2(void){
	fprintf (stderr, "usage: mm <N> <M> <P> <trials>\n");
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
