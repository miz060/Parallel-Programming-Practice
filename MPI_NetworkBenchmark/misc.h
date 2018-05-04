/**
 * File name: misc.h
 * File description: define helper method for mmmw to reemove redundancy
 * Name:  Mingcheng Zhu
 * Email: zhumc11@gmail.com
 * Date:  Feb 10, 2018
 */
#define DATA 10
#define DONE 0


/* cpu_time
 * function by John Burkardt
 */
extern double cpu_time ( void );

/* mult
 * purpose: do matrix multiplication of AxB, save result to C
 * param: input matrix A,B, output matrix C, matrix size N, M, P
 */
extern void mult(double **, double **, double **, int, int, int);


/* multM
 * purpose: do block matrix multiplication of AxB, save result to C
 * param: input matrix A,B, output matrix C, matrix size N, M, P, blkSize, 
 *        number of processes np, arrary to record worker's task#,
 *        array to record worker's computation time
 */
extern void multM(double ****, double ****, double ****, int, int, int, int, int, int*, double*);


/* add
 * purpose: do matrix addition of matrix liner B on A
 * param: added matrix A, input matrix liner B, matrix size
 */
extern void add(double **, double *, int);


/* zero
 * purpose: zero fill input matrix
 * param: matrix to zero fill, matrix size
 */
extern void zero(double *, int);


/* fill_m
 * purpose: called by master to compress input matrix A and B into a matrix liner
 * param: input matrix A and B, output liner, matrix size, row and col in master matrix
 */
extern void fill_m(double**, double**, double*, int, int, int);


/* fill_w
 * purpose: called by worker to compress input matrix temp into a matrix liner
 * param: input matrix temp, output liner, matrix size, and its row&col in master,
 *        worker computation time duration
 */
extern void fill_w(double**, double*, int, int, int, double);


/* unload_w
 * purpose: called by worker to unload input liner into two matrix 
 * param: input liner, output matrix A and B, matrix size
 */
extern void unload_w(double*, double**, double**, int);


/* sendToWorker
 * purpose: called by master to send mult liner to workers
 */
extern void sendToWorker(int, double*, MPI_Request*, int, int);


/* waitForOne
 * purpose: called by master to wait for a worker to finish work 
 *          and add its result to the master matrix
 */
extern int waitForOne(double****, MPI_Request*, int, double*, int, double*);


/* waitForAll
 * purpose: called by master to wait for all workers to finish its work 
 *          and add their results to the master matrix
 */
extern void waitForAll(double**** , MPI_Request*, int, int, int, double*, double*);


/* usage
 * purpose: print out the usage message to stdout
 */
extern void usage( void );

