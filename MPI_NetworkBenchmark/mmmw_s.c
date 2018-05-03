/* mm.c 
 * Sample matrix multiplication in C 
 * Author: Philip Papadopoulos
 * Email: ppapadopoulos@ucsd.edu
 * UCSD Course: cse160 
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
double cpu_time ( void );
void mult(double **, double **, double **, int, int, int);
void multM(double ****, double ****, double ****, int, int, int, int);
void add(double **, double **, int);
void zero(double *, int);
void usage( void );


/* ===========================================================================
 * Usage:  mm <N> <trials> 
 * New Usage: mm <N> <M> <P> <trials>
 */
int main(int argc, char * argv[]){
	if (argc != 6) usage();
  int blkSize=atoi(argv[1]);
	int N=atoi(argv[2]); 
  int M=atoi(argv[3]);
  int P=atoi(argv[4]);
	char* outfile=argv[5];

	double ****A, ****B, ****C;
	double ***Alinear, ***Blinear, ***Clinear;
  // A is NxM, B is MxP, C is NxP
	A = (double ****) calloc(N, sizeof (double ***));
	B = (double ****) calloc(M, sizeof (double ***));
	C = (double ****) calloc(N, sizeof (double ***));
	Alinear = calloc(N * M, sizeof (double**));
	Blinear = calloc(M * P, sizeof (double**));
	Clinear = calloc(N * P, sizeof (double**));
	int i;
  int j;
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
  printf("Below is A:\n");
  for (i=0; i<N; i++){
    for (j=0; j<M; j++){
      printf("A[%d][%d] is block:\n", i, j);
      pblock = A[i][j];
      for(m=0; m<blkSize; m++){
        for(n=0; n<blkSize; n++){
          printf("%f ", pblock[m][n]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }
  printf("\nBelow is B:\n");
  for (i=0; i<M; i++){
    for (j=0; j<P; j++){
      printf("B[%d][%d] is block:\n", i, j);
      pblock = B[i][j];
      for(m=0; m<blkSize; m++){
        for(n=0; n<blkSize; n++){
          printf("%f ", pblock[m][n]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

	multM(A,B,C,N,M,P,blkSize);
  printf("mutiplied %dx%d matrix block with %dx%d matrix block with blockSize %d\n", N,M,M,P, blkSize);

  // debug print, remeber to comment out 
  printf("Below is C:\n");
  for (i=0; i<N; i++){
    for (j=0; j<P; j++){
      pblock = C[i][j];
      printf("C[%d][%d] is:\n", i, j);
      for(m=0; m<blkSize; m++){
        for(n=0; n<blkSize; n++){
          printf("%f ", pblock[m][n]);
        }
        printf("\n");
      }
    }
    printf("\n");
  }

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

void multM(double ****A, double ****B, double ****C, int N, int M, int P, int blkSize){
  int i, row, col;
  double** temp = (double **) calloc(blkSize, sizeof (double *));
	double* tempLiner = calloc(blkSize*blkSize, sizeof (double));
  for (i=0; i<blkSize; i++){
		temp[i] = tempLiner + i * blkSize;
  }
  
  for (row = 0; row < N; row ++){
    for (col = 0; col < P; col ++){
      zero(tempLiner, blkSize);
      for (i = 0; i < M; i++){
        mult(A[row][i], B[i][col], temp, blkSize, blkSize, blkSize);
        add(C[row][col], temp, blkSize);
      }
    }
  }
}

void add(double** C, double **temp, int size){
  int i,j;
  for(i=0; i<size; i++){
    for(j=0; j<size; j++){
      C[i][j] += temp[i][j];
    }
  }
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
