#include <pthread.h>
/** 
 * File name: heat2d_solver.h
 * File description: declares the solver function for heat2d.c.
 *                   use row decomposition to do the calculation.
 *                   parallelism achieved through pthread.
 * Name:  Mingcheng Zhu
 * Email: zhumc11@gmail.com
 * Date:  Feb 25, 2018
 */


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
int heat2dSolve(int M, int N, double eps, int print, volatile double **u, double *tol, pthread_t* thread_handles, int thread_count);
