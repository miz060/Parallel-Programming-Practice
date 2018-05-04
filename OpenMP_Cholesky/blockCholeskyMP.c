#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <omp.h>
/** 
 * File name: blockCholeskyMP.c
 * File description: parallel block cholesky LU factorization
 * Name:  Mingcheng Zhu
 * Email: zhumc11@gmail.com
 * Date:  Mar 18, 2018
 */
#define ARGS 4
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define uSECtoSEC 1.0E-6
#define THRESH 1e-13
#define SCALE 100.0
#define TRUE 1
#define FALSE 0

// Macro to define index into a linear array for 2D indexing. Stored 
// row by row.
#define IDX(i,j,n) ((i*n)+j)

int thread_count;
/* return a clock value with usec precision */
double get_clock() {
  struct timeval tv;
  int status;
  status = gettimeofday(&tv, (void *) 0);
  if(status<0)
    printf("gettimeofday error");
  return (tv.tv_sec * 1.0 + tv.tv_usec * 1.0E-6);
}

  void
printMatrix(double *A, int N){
  int i,j;
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++)
      printf("%lf ", A[IDX(i,j,N)]);
    printf("\n");
  }
}


/* Multiply A*A^T.  
 * double *result - result matrix
 * double * A - source matrix
 * int N - size of matrix (N x N);
 * int lowerT - A is lower triangular
 */
int multT(double *result, double *A, int N, int lowerT){
  //printf("inside multT\n");
  int i,j,k;
  //bzero(result, N*N*sizeof(double));
  # pragma omp parallel for num_threads(thread_count) schedule(dynamic)
  for(i = 0; i < N; i++){
    /* Result is symmetric, just compute upper triangle */
    # pragma omp parallel for num_threads(thread_count) schedule(dynamic)
    for(j = i; j < N; j++){
      double sum = 0.0;
      /* if A is lower Triangular don't multiply zeroes */
      # pragma omp parallel for num_threads(thread_count) \
        reduction(+: sum) schedule(static)
      for(k = 0; k < (!lowerT ? N : j+1) ; k++){
        //printf("product of A21(%d,%d):%f, A21(%d,%d):%f\n", i,k,A[IDX(i,k,N)], j,k,A[IDX(j,k,N)]);
        sum += A[IDX(i,k,N)] * A[IDX(j,k,N)];
      }
      result[IDX(i,j,N)] = sum;
      result[IDX(j,i,N)] = sum; /*enforce symmetry */
    }
  }
  return 0;
}


/** mult1
 * Multiply matrices, C = A x B.
 * Assumes Matrices are already allocated.
 * values are not checked for sanity
 */
/* Multiply L21*L21^T.  
 * double *result - result matrix
 * double * A - source matrix
 * int N - size of matrix (N x N);
 * int blk - the block size
 */
int multT1(double *result, double *A, int N, int blk){
  //printf("inside multT\n");
  int i,j,k;
  //bzero(result, N*N*sizeof(double));
  # pragma omp parallel for num_threads(thread_count) schedule(dynamic)
  for(i = 0; i < N; i++){
    /* Result is symmetric, just compute upper triangle */
    for(j = i; j < N; j++) {
      double sum = 0.0;
      /* if A is lower Triangular don't multiply zeroes */
      for(k = 0; k < blk; k++){
        //printf("product of L21(%d,%d):%f, L21(%d,%d):%f\n", i,k,A[IDX(i,k,blk)], j,k,A[IDX(j,k,blk)]);
        sum += A[IDX(i,k,blk)] * A[IDX(j,k,blk)];
      }
      result[IDX(i,j,N)] = sum;
      result[IDX(j,i,N)] = sum; /*enforce symmetry */
    }
  }
  return 0;
}


/* Sub A22 with L21xL21^T
 * double * A- A22
 * double * B- L21xL21^T
 * int N - size of matrix (N x N);
 * int itr - the current iteration in block cholesky
 * int blk - the block size
 */
int sub(double *A, double *B, int N, int itr, int blk){
  int i,j,k;
  # pragma omp parallel for num_threads(thread_count) schedule(static) collapse(2)
  for(i=0; i<N-blk*(itr+1); i++){
    for(j=0; j<N-blk*(itr+1); j++){
      //printf("A(%d,%d) -= B(%d,%d)\n", (i+blk*(itr+1)), (j+blk*(itr+1)), i,j);
      A[IDX((i+blk*(itr+1)),(j+blk*(itr+1)),N)] -= B[IDX(i,j,(N-blk*(itr+1)))];
    }
  }
}


/* Validate that A ~= L*L^T 
 * double * A -  NxN symmetric positive definite matrix
 * double * L -  NxN lower triangular matrix such that L*L^T ~= A
 * int N      -  size of matrices
 * thresh     -  threshold considered to be zero (1e-14 good number)
 *
 * Returns # of elements where the residual abs(A(i,j) - LL^T(i,j)) > thresh
 */
int validate(double *A, double * L, int N, double thresh){
  double *R = malloc(N*N * sizeof(double));
  multT(R,L,N,TRUE);
  int badcount = 0;
  int i,j;
  double rdiff; /* relative difference */
  # pragma omp parallel for num_threads(thread_count) schedule(static) collapse(2)
  for (i = 0 ; i < N; i++){
    for(j = 0; j < N; j ++){
      rdiff = fabs((R[IDX(i,j,N)] - A[IDX(i,j,N)])/A[IDX(i,j,N)]);
      if (rdiff > thresh){
        printf("(%d,%d):R(i,j):%f,A(i,j):%f (delta: %1.15f)\n",
            i,j,R[IDX(i,j,N)],A[IDX(i,j,N)],
            rdiff);

        badcount++; 
      }
    }
  }
  //fprintf(stderr, "badcount=%d\n", badcount);
  free(R);
  return badcount;
}


/* Initialize the N X N  array with Random numbers
 * In such a way that the resulting matrix in Symmetric
 * Positive definite (hence Cholesky factorization is valid)
 * args:
 * 	int N - size of the array
 * 	int trueRandom  - nonzero, seed with time of day, 0 don't seed
 *	double **A - N x N double array, allocated by caller
 */
void init_array(int N, int trueRandom, double *A) {
  int i,j,k;
  struct drand48_data rbuf;
  if (trueRandom)
    srand48_r((long int) time(NULL),&rbuf);
  else
    srand48_r(1L,&rbuf);

  double *B = calloc(N * N, sizeof(double));

  //printf("Random number generation\n");
# pragma omp parallel for num_threads(thread_count) schedule(static) collapse(2)
  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      drand48_r(&rbuf,&B[IDX(i,j,N)]);
      B[IDX(i,j,N)] *= SCALE;
    }
  }
  //printf("done random number generation\n");
  // printMatrix(B,N);
  /* Compute B*B^T to get symmetric, positive definite*/
  multT(A,B,N,0);
  free (B);
}


/* Compute the Cholesky Decomposition of A 
 * L - NxN result matrix, Lower Triangular L*L^T = A
 * A - NxN symmetric, positive definite matrix A
 * N - size of matrices;
 */
// i is j and j is i in the formula
void cholesky(double *L, double *A, int N, int N_g){
  int i,j,k;
  //bzero(L,N*N*sizeof(double));
  double temp;
  for (i = 0; i < N; i++){
    for (j = 0; j < (i+1); j++) {
      temp = 0;
      /* Inner product of ith row of L, jth row of L */
      for (k = 0; k < j; k++){
        temp += L[IDX(i,k,N)] * L[IDX(j,k,N)];
      }
      if (i == j){
        //printf("modify (%d,%d)\n",i,j);
        L[IDX(i,j,N)] = sqrt(A[IDX(i,i,N_g)] - temp);
      }
      else {
        //printf("modify (%d,%d)\n",i,j);
        L[IDX(i,j,N)] = (A[IDX(i,j,N_g)] - temp)/ L[IDX(j,j,N)];
      }
    }
  }
}


/* Compute the Cholesky Decomposition of A, block version
 * L - NxN result matrix, Lower Triangular L*L^T = A
 * A - NxN symmetric, positive definite matrix A
 * N - size of matrices;
 * blk - the block size
 * itr - the starting iteration
 */
void blockCholesky(double *L, double*A, int N, int blk, int itr){
  int i,j,k;
  double temp;
  double* L21L21T;
  double* transpose;
  double* L11;
  double* L21;

  // local work space
  L11 = (double *)malloc(blk*blk* sizeof(double));
  L21L21T = (double *)malloc(N*N* sizeof(double));
  L21 = (double *)malloc(N*blk* sizeof(double));

  for(itr=0; itr<N/blk; itr++){

    cholesky(L11, A+N*blk*itr+blk*itr, blk, N);

    // calculate bottom rectangular matrix L21
    for(j=0; j<blk; j++){
      // col by col
      //printf("calculating %d col of L\n", j);
      for(i=0; i<N-blk*(itr+1); i++){
        temp = A[IDX((i+(itr+1)*blk), (j+itr*blk), N)];
        for(k=0; k<=j-1; k++){
          temp -= L11[IDX(j,k,blk)] * L21[IDX(i,k,blk)];
        }
        //printf("Calculating L21(%d,%d) itr is %d\n", i, j, itr);
        //printf("i+(itr+1)*blk is %d, j+itr*blk is %d\n", (i+(itr+1)*blk), (j+itr*blk));
        //printf("A(3,0) is %f\n", A[IDX(3,0,N)]);
        //printf("A21(%d,%d) is %f\n",i,j, A[IDX((i+(itr+1)*blk), (j+itr*blk), N)]);
        //printf("L11(%d,%d) is %f\n",j,j,L11[IDX(j,j,N)]);
        // write L21 to L
        L21[IDX(i,j,blk)] =  temp / L11[IDX(j,j,blk)];
      }
    }

    // update A22
    multT1(L21L21T, L21, (N-blk*(itr+1)), blk);
    //printf("AFTER  transp[3,0] is %f\n", transpose[IDX(3,0,(N-blk*(itr+1)))]);

    sub(A, L21L21T, N, itr, blk); 

    double* L11_ = L+N*blk*itr+blk*itr;
    double* L21_ = L11_+N*blk;
    # pragma omp parallel for num_threads(thread_count) schedule(dynamic)
    for(i=0; i<(N-blk*(itr+1)); i++){
      # pragma omp parallel for num_threads(thread_count) schedule(static)
      for(j=0; j<blk; j++){
        L21_[IDX(i,j,N)] = L21[IDX(i,j,blk)];
      }
    }

    # pragma omp parallel for num_threads(thread_count) schedule(static)
    for(i=0; i<blk; i++){
      for(j=0; j<=i; j++){
        L11_[IDX(i,j,N)] = L11[IDX(i,j,blk)];
      }
    }
  }
  /*printf("L final is \n");
    for(i=0; i<N; i++){
    for(j=0; j<N; j++){
    printf("%15.8f ", L[IDX(i,j,N)]);
    }
    printf("\n");
    }*/
  free(L21L21T);
  free(L11);
  free(L21);
}


int main(int argc, char* argv[]) {
  int n;
  int blk;
  int i,j,k;
  double ts, te; /* Starting time and ending time */
  double *A, *L, *A_copy;
  double temp;
  if(argc < ARGS || argc > ARGS){
    fprintf(stderr,"Wrong # of arguments.\nUsage: %s <K> <blksize> <nthreads>\n",argv[0]);
    return -1;
  }
  n = atoi(argv[1]);
  blk = atoi(argv[2]);
  thread_count = atoi(argv[3]);
  int N = n*blk;
  A = (double *)malloc(N*N*sizeof(double));
  A_copy = (double *)calloc(N*N,sizeof(double));
  L = (double *)calloc(N*N,sizeof(double));
  printf(" %dX%d block matrix OpenMP with blk=%d, %d threads\n", n,n, blk, thread_count);

  ts=get_clock();
  init_array(N,0,A);
  te=get_clock();
  printf("intialize:      %f seconds\n", te - ts);
  // printf("Initial matrix:\n");
  // printMatrix(A,n);

  /*printf("A is \n");
    for(i=0; i<N; i++){
    for(j=0; j<N; j++){
    printf("%15.8f ", A[IDX(i,j,N)]);
    }
    printf("\n");
    }*/

  // step 1: copy the element of A into L
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      A_copy[IDX(i,j,N)] = A[IDX(i,j,N)];
    }
  }
  // step 2: set the upper triangular to be zero
  /*for(i=0; i<N; i++){
    for(j=i+1; j<N; j++){
    L[IDX(i,j,N)] = 0;
    }
    }*/

  ts=get_clock();
  /*Serial decomposition*/
  blockCholesky(L,A,N, blk, 0);
  te=get_clock();
  printf("block cholesky: %f seconds\n", te - ts);
  //printf("Decomposed matrix:\n");
  //printMatrix(L,n);
  int badcount = 0;
  ts=get_clock();
  badcount = validate(A_copy,L,N,THRESH);
  te=get_clock();
  printf("validation:     %f seconds\n", te - ts);
  if (badcount == 0)
    printf("solution validates\n");
  else
    printf("solution is invalid, %d elements above threshold\n",badcount);
  free(A);
  free(A_copy);
  free(L);
  return badcount;
}
