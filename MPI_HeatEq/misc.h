/** 
 * File name: misc.h
 * File description: define helper method for heat2d programs to remove redundancy.
 *                   header file for misc.c
 * Name:  Mingcheng Zhu
 * Email: zhumc11@gmail.com
 * Date:  Jan 28, 2018
 */
#define FALSE 0
#define TRUE 1
/* Function prototypes  */

/* usage
 * purpose: print usage info of heat2d to stderr
 */
extern int usage();


/* cLineValid
 * purpose: validate command line arguments for the calling program
 * param:   argv1,2,3,4,5,6,7 - command line arguments of the calling program
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
extern int cLineValid(char* argv1, char* argv2, char* argv3, 
               char* argv4, char* argv5, char* argv6,
               char* argv7);


/* serial
 * purpose: do serial computation of the convergence iterations
 * param:   M, N, u, w, mean, eps, Tl, Tr, Tt, Tb - the relative variables
 *          in the calling program
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
extern int serial(int M, int N, double** u, double** w, double mean, double eps, 
           double Tl, double Tr, double Tt, double Tb);


/* function by John Burkardt */
extern double cpu_time ( void );


/* serial_write
 * purpose: do serial write of a complete matrix to the outfile
 * param:   outfile, M, N, w, fp - the matching variables in the calling program
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
extern int serial_write(char* output_file, int M, int N, double**w, FILE* fp);


/* serial_write_s
 * purpose: do serial write of a complete matrix to the outfile
 *          the same as serial_write, except w would be a char* matrix
 *          (this one is special for heat2dcol.c
 * param:   outfile, M, N, w, fp - the matching variables in the calling program
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
extern int serial_write_s(char* output_file, int M, int N, char***w, FILE* fp);


/* calculation
 * purpose: do serial calculation of the matrix convergence
 * param:   u, w, my_rank, comm_sz, length, top, bottom - the matching variables in the 
 *          calling program
 *          X - the array length of the partial matrix. N for heat2drow, M for heat2dcol
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
extern double calculation(double** u, double**w, int X, int my_rank, int comm_sz, 
                          int length, double* top, double* bottom);
