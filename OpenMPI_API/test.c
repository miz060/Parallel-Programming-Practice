/*
 * =====================================================================================
 *
 *       Filename:  test.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/01/18 16:58:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingcheng Zhu (A92047564), miz060@ucsd.edu
 *        Company:  ucsd
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <assert.h>
#include <string.h>
#include "mp.h"

#define TAG 9

void *func(void * arg)
{
  int i,j;
  int id = MP_Rank(pthread_self());
  id = MP_Rank(pthread_self());
  MSGRQST requests[2];
  char* rbuf = malloc(3*sizeof(char));
  char* sbuf = malloc(3*sizeof(char));
  
  for(j=0; j<2; j++){
  printf("ITERATION %d starts\n", j);
  if(id == 0){
    sbuf[0] = 'a'+(char)j*3;
    sbuf[1] = 'b'+(char)j*3;
    sbuf[2] = 'c'+(char)j*3;
    
    printf("Thread[%d] send message %d,%d,%d, iteration %d\n", id, sbuf[0], sbuf[1], sbuf [2], j);
    //printf("Thread[%d] before iSend\n", id);
    iSend(sbuf, 3, 0, 1, TAG+j, &requests[0]);
    //printf("Thread[%d] before msgWait\n", id);

    //sleep(1);
    //msgWait(&requests[0]);
    //msgWait(&MSGRQST_NULL);
    //printf("Thread[%d] finish msgWait, iteration %d\n", id, j);
    //free(sbuf);
  }
  else if(id == 1){

    /*printf("Thread[%d] receive buf before iRecv (iteration %d):", id, j);
    for(i=0; i<3; i++){
      printf("%d ", rbuf[i]);
    }
    printf("\n");*/

    //printf("Thread[%d] before iRecv\n", id);
    iRecv(rbuf, 3, ANY_SRC, 1, TAG+j, &requests[0]);    
    //printf("Thread[%d] before msgWait\n", id);
    //sleep(1);
    //msgWait(&requests[0]);
    //msgWait(&MSGRQST_NULL);
  }

  if(id ==0){
    char* ackrbuf = malloc(3*sizeof(char));
    //printf("Thread[%d] before iRecv ACK\n", id);
    iRecv(ackrbuf, 3, ANY_SRC, 0, 100+j, &requests[1]);
    //printf("Thread[%d] before msgWait ACK\n", id);
    msgWaitAll(&requests[0], 2);
    //sleep(2);
    printf("Thread[%d] got ack message of len %d from thread %d with tag %d: (iteration %d)", id, getLen(&requests[1]), getRank(&requests[1]), getTag(&requests[1]),j);
    for(i=0; i<3; i++){
      printf("%d ", ackrbuf[i]);
    }
    printf("\n");
    //free(rbuf);
    //free(ackrbuf);
  }
  else if(id==1){

    //printf("Thread[%d] before iSend ACK\n", id);
    iSend(rbuf, 3, 1, 0, 100+j, &requests[1]);
    //printf("Thread[%d] before msgWait ACK\n", id);
    msgWaitAll(&requests[0], 2);
    //free(rbuf);
    printf("Thread[%d] receive message of len %d from thread %d with tag %d: (iteration %d)", id, getLen(&requests[0]), getRank(&requests[0]), getTag(&requests[0]), j);
    for(i=0; i<3; i++){
      printf("%d ", rbuf[i]);
    }
    printf("\n");
  }
  }
}


main(int argc, char *argv[])
{
  int i;
  int nthread = atoi(argv[1]);
	MP_Init(nthread);

  pthread_t* thread_handles = malloc(nthread*sizeof(pthread_t));

	/* create the benchmark threads. wait for both to complete */
	for (i=0; i < nthread; i++)
	{
		pthread_create(thread_handles+i,NULL, func, NULL);
	}
	for (i=0; i < nthread; i++) 
	{
		pthread_join(thread_handles[i],NULL);
	}
  printf("Hello from main thread\n");
  free(thread_handles);
	MP_Finalize();
}

