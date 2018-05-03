/**
 * CSE 160: Programming Assignment 3
 * File name: msgbench.c
 * File description: benchmark message passing
 * Name:  Mingcheng Zhu
 * PID:   A92047564
 * Email: miz060@ucsd.edu
 * Date:  Feb 10, 2018
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define DATA_0 10
#define DATA_1 11

double cpu_time ( void );

/* usage
 * purpose: print out the usage message
 */
int usage(){
  fprintf(stderr, "usage: msgbench <msgsize> <trials> [bidirectional]\n");
  return 1;
}

/* main
 * purpose: main function for the message passing benchmarking
 * param:   argc - cmd argument counts, argv - cmd argument array
 * return exit status, 0 if normal, non-zero if error occurs
 */
int main ( int argc, char *argv[] );
int main(int argc, char* argv[]){

  int size;
  int round;
  int bidirectional=0;
  int i;
  int j;
  char *send0;
  char *recv0;
  char *send1;
  char *recv1;
  int comm_sz;
  int my_rank;
  double startTime, endTime, avgTime;
  double speed;
  // set up MPI components
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  MPI_Request request[4];
  MPI_Status status[4];
  
  // command line validation
  if (argc != 3 && argc != 4){
    if(my_rank == 0){
      usage();
    }
    MPI_Finalize();
    return -1;
  }

  // check if bidirectional
  if(argc == 4){
    bidirectional = 1;
  }
    
  size = atoi(argv[1]); 
  round = atoi(argv[2]); 

  // set up the message buffers
  send0 = (char*) calloc(size, sizeof (char));
  recv0 = (char*) calloc(size, sizeof (char));
  send1 = (char*) calloc(size, sizeof (char));
  recv1 = (char*) calloc(size, sizeof (char));

  // populate the send buffer to be non-zero as process 0
  for (i=0; i<size; i++){
    send0[i] = random()%100;
    send1[i] = random()%101;
  }

  // process 0 case
  if(my_rank == 0){
    // process 0 start the timer
    startTime = cpu_time();

    for(j=0; j<round; j++){
      for(i=0; i<4; i++){
        request[i] = MPI_REQUEST_NULL;
      }

      if(!bidirectional){
        MPI_Isend(&send0[0], size, MPI_CHAR, 1, DATA_0+2*j, MPI_COMM_WORLD, &request[0]);
        MPI_Irecv(&recv0[0], size, MPI_CHAR, 1, DATA_0+2*j, MPI_COMM_WORLD, &request[1]);

        MPI_Waitall(2, request, status);
      }
      else{
        MPI_Isend(&send0[0], size, MPI_CHAR, 1, DATA_0+2*j, MPI_COMM_WORLD, &request[0]);
        MPI_Irecv(&recv1[0], size, MPI_CHAR, 1, DATA_1+2*j, MPI_COMM_WORLD, &request[1]);
        MPI_Waitall(2, request, status);

        MPI_Isend(&recv1[0], size, MPI_CHAR, 1, DATA_1+2*j, MPI_COMM_WORLD, &request[2]);
        MPI_Irecv(&recv0[0], size, MPI_CHAR, 1, DATA_0+2*j, MPI_COMM_WORLD, &request[3]);
        MPI_Waitall(2, &request[2], status);
      }
    }

    // process 0 stop the timer
    endTime = cpu_time();
    avgTime = (endTime - startTime)/(double) round;
    speed = (double)(2*size)/avgTime;
    if(avgTime != 0){
      if(!bidirectional){
        printf("%d bytes (%d trials) %f time, speed: %f bytes/s\n", size, round, avgTime, speed);
      }
      else{
        speed = (double)(4*size)/avgTime;
        printf("%d bytes (%d trials) %f time, speed: %f bytes/s (bidirectional)\n", size, round, avgTime, speed);
      }
    }
    else{
      if(!bidirectional){
        printf("%d bytes (%d trials) %f time, speed: %f bytes/s, msgSize is too small, please try a larger size\n", 
        size, round, avgTime, speed);
      }
      else{
        printf("%d bytes (%d trials) %f time, speed: %f bytes/s, msgSize is too small, please try a larger size (bidirectional)\n", 
        size, round, avgTime, speed);
      }
    }
  }

  // process 1 case
  else if(my_rank == 1){
    // process 1 receive and send
    for(j=0; j<round; j++){

      for(i=0; i<4; i++){
        request[i] = MPI_REQUEST_NULL;
      }
        
      if(!bidirectional){
        MPI_Irecv(&recv0[0], size, MPI_CHAR, 0, DATA_0+2*j, MPI_COMM_WORLD, &request[0]);
        MPI_Isend(&recv0[0], size, MPI_CHAR, 0, DATA_0+2*j, MPI_COMM_WORLD, &request[1]);
        MPI_Waitall(2, request, status);
      }
      else{
        MPI_Isend(&send1[0], size, MPI_CHAR, 0, DATA_1+2*j, MPI_COMM_WORLD, &request[0]);
        MPI_Irecv(&recv0[0], size, MPI_CHAR, 0, DATA_0+2*j, MPI_COMM_WORLD, &request[1]);
        MPI_Waitall(2, request, status);

        MPI_Isend(&recv0[0], size, MPI_CHAR, 0, DATA_0+2*j, MPI_COMM_WORLD, &request[2]);
        MPI_Irecv(&recv1[0], size, MPI_CHAR, 0, DATA_1+2*j, MPI_COMM_WORLD, &request[3]);
        MPI_Waitall(2, &request[2], status);
      }
    }
  
  }

  free(send0);
  free(send1);
  free(recv0);
  free(recv1);
  MPI_Finalize();
  return 0;

}
    
    
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









