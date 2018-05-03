#include <pthread.h>
#include <sched.h>
#include <semaphore.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "cs160mp.h"
/** 
 * CSE 160: Programming Assignment 5
 * File name: cs160mp.c
 * File description: Implement the message passing API with pthread
 * Name:  Mingcheng Zhu
 * PID:   A92047564
 * Email: miz060@ucsd.edu
 * Date:  Mar 4, 2018
 */


// global variables
MSGRQST ** sendQueue;
MSGRQST ** recvQueue;

// condition variables for broadcast and wait
pthread_cond_t* conditions;
pthread_mutex_t* mutexes;

// mutex to guard the global queue;
pthread_mutex_t* queue_mutex;
pthread_mutex_t* id_mutex;

// global varible with no need to be synchronized
int thread_count;
int sqLength;
int rqLength;
int id_generator = 0;
pthread_t* ids;
int search_match(int src, int dest, int tag, int type, int id, int finish);


/* MP_Init
 * purpose: Initialize any internal data structures needed to implement message passing.
 *          Called once per program.
 * param:   nthreads - number of threads created
 */
int MP_Init(int nthreads){
  int i;
  // this is an arbitary queue size. Will be realloced if not large enough
  sqLength = 1000;
  rqLength = 1000;
  sendQueue = calloc(sqLength, sizeof(MSGRQST*));
  recvQueue = calloc(sqLength, sizeof(MSGRQST*));
  ids = calloc(nthreads, sizeof(pthread_t));
    
  // set up mutexes
  queue_mutex = malloc(sizeof(pthread_mutex_t));
  id_mutex = malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(queue_mutex, NULL);
  pthread_mutex_init(id_mutex, NULL);

  // set up the conditional variable for msgWait
  conditions = malloc(nthreads * sizeof(pthread_cond_t));
  mutexes = malloc(nthreads * sizeof(pthread_mutex_t));
  for(i=0; i<nthreads; i++){
    pthread_cond_init(&conditions[i], NULL);
    pthread_mutex_init(&mutexes[i], NULL);
  }
  thread_count = nthreads;
  return 0;
}


/* MP_Size
 * purpose: Return the number of threads passed to MP_Init
 */
int MP_Size(){
  return thread_count;
}


/* MP_Rank
 * purpose: Give a rank (0 .. MP_Size()) to the thread with a particular id
 * return:  the given thread rank
 */
int MP_Rank(pthread_t thread){
  int i;
  pthread_mutex_lock(id_mutex); 
  // check if the thead has an id assigned
  for(i=0; i<thread_count; i++){
    if(ids[i] == thread){
      pthread_mutex_unlock(id_mutex);
      return i;
    }
  }
  // assign an id to the thread
  ids[id_generator] = thread;
  i = id_generator;
  id_generator++;
  pthread_mutex_unlock(id_mutex);
  return i;
}


/* MP_Finalize
 * purpose: Clean up all allocated memory. 
 *          Indicates no more message passing will be performed
 *          Called once after all message passing is complete
 */
int MP_Finalize(){
  pthread_mutex_destroy(queue_mutex);
  pthread_mutex_destroy(id_mutex);
  // free MSGRQST's memory
  int i;
  for(i=0; i<sqLength; i++){
    free(sendQueue[i]);
  }
  for(i=0; i<rqLength; i++){
    free(recvQueue[i]);
  }
  
  for(i=0; i<thread_count; i++){
    pthread_cond_destroy(&conditions[i]);
    pthread_mutex_destroy(&mutexes[i]);
  }
  free(mutexes);
  free(conditions);

  free(queue_mutex);
  free(id_mutex);
  free(sendQueue);
  free(recvQueue);
  free(ids);
  return 0;
}


#define STATUS_NORMAL 0
#define STATUS_ERR -1
#define SEND 0
#define RECV 1


/* iSend
 * purpose: Send a memory buffer (no type) of paticular length from src to dest using specific tag.
 *          Asynchronous
 * param:   buf - the buffer with the message to be sent.
 *          len - the length of the message to be sent.
 *          src - the source thread's rank.
 *          dest- the destination thread's rank.
 *          tag - the message passing tag
 *          request - pointer to the MSGRQST of this send operation.
 */
int iSend(void *buf, int len, int src, int dest, int tag, MSGRQST * request){
  pthread_mutex_lock(queue_mutex); 
  MSGRQST* lrequest = malloc(sizeof(MSGRQST));
  request -> src = src;
  request -> dest = dest;
  request -> tag = tag;
  request -> len = len;
  request -> type = SEND;
  request -> message = buf;
  lrequest -> src = src;
  lrequest -> dest = dest;
  lrequest -> tag = tag;
  lrequest -> len = len;
  lrequest -> type = SEND;
  lrequest -> message = buf;
  lrequest -> finish = 0;

  int i,j;
  while(1){
    for(i=0; i<sqLength; i++){
      // find an empty spot on the queue
      if(sendQueue[i] == NULL){
        sendQueue[i] = lrequest;
  
        // wake up the early recv thread waiting on the current (src) thread
        pthread_mutex_lock(&mutexes[src]);
        pthread_cond_broadcast(&conditions[src]);

        pthread_mutex_unlock(queue_mutex);
        pthread_mutex_unlock(&mutexes[src]);
        return STATUS_NORMAL;
      }
    }
    // alloc more memory if needed
    for(j=0; j<sqLength; j++){
      free(sendQueue[j]);
    }
    sqLength += 1000;
    sendQueue = realloc(sendQueue, sqLength*sizeof(MSGRQST*));
    memset(sendQueue, 0, sqLength*sizeof(MSGRQST*));
  }
  // should never be reached
  pthread_mutex_unlock(queue_mutex);
  return STATUS_ERR;  
}


/* iRecv
 * purpose: Recv into memory buffer (no type) of size len 
 *          from src (may be wildcard) to dest using specific tag (may be wildcard).
 *          Asynchronous.
 * param:   buf - the buffer to store the received message.
 *          len - the length of receive buffer
 *          src - the source thread's rank.
 *          dest- the destination thread's rank.
 *          tag - the message passing tag
 *          request - pointer to the MSGRQST of this receive operation.
 */
int iRecv(void *buf, int len, int src, int dest, int tag, MSGRQST * request){
  pthread_mutex_lock(queue_mutex);
  MSGRQST* lrequest = malloc(sizeof(MSGRQST));
  request -> src = src;
  request -> dest = dest;
  request -> tag = tag;
  request -> len = len;
  request -> type = RECV;
  request -> message = buf;
  lrequest -> src = src;
  lrequest -> dest = dest;
  lrequest -> tag = tag;
  lrequest -> len = len;
  lrequest -> type = RECV;
  lrequest -> message = buf;
  lrequest -> finish = 0;

  int i,j;
  while(1){
    for(i=0; i<rqLength; i++){
      // find an empty spot on the queue
      if(recvQueue[i] == NULL){
        recvQueue[i] = lrequest;

        // wake up the early send thread waiting on the current (dest) thread
        pthread_mutex_lock(&mutexes[dest]);
        pthread_cond_broadcast(&conditions[dest]);

        pthread_mutex_unlock(queue_mutex);
        pthread_mutex_unlock(&mutexes[dest]);
        return STATUS_NORMAL;
      }
    }
    // alloc more memory if needed
    for(j=0; j<rqLength; j++){
      free(recvQueue[j]);
    }
    rqLength += 1000;
    recvQueue = realloc(recvQueue, rqLength*sizeof(MSGRQST*));
    memset(recvQueue, 0, rqLength*sizeof(MSGRQST*));
  }
  // should never be reached
  pthread_mutex_unlock(queue_mutex);
  return STATUS_ERR;  
}


/* msgWait
 * purpose: Wait for a specific message to complete
 * param:   request - the request to be waited.
 */
int msgWait( MSGRQST * request){
  pthread_mutex_lock(queue_mutex);
  int src = request -> src;
  int dest = request -> dest;
  int len = request -> len;
  int tag = request -> tag;
  int type = request -> type;
  void* message = request -> message;
  int i;
  int match = -1;
  int self;

  // return on MSGRQST_NULL
  if(src == -1){
    pthread_mutex_unlock(queue_mutex);
    return STATUS_NORMAL;
  }
  // return if the message passing is finished
  if(request -> finish == 1){
    pthread_mutex_unlock(queue_mutex);
    return STATUS_NORMAL;
  }

  if(type == SEND){
    //printf("Thread[%d] start looking for a recv request\n", src);
    // search for the recv request
    self = search_match(src, dest, tag, SEND, src, 1);
    // the other thread has finished the request for this one
    if(self == -1){
      match = search_match(src, dest, tag, SEND, dest, 0);
      request -> len = recvQueue[match] -> len;
      pthread_mutex_unlock(queue_mutex);
      return STATUS_NORMAL;
    }
    match = search_match(src, dest, tag, RECV, src, 1);

    // do the finishing job as the late one
    if(match != -1){
      //printf("Thread[%d] is the late send, match is %d\n", src, match);
      len = recvQueue[match] -> len;
      request -> len = len;
      // do the message transmission
      memcpy(recvQueue[match]->message, message, len*sizeof(char));
      sendQueue[self] -> finish = 1;
      recvQueue[match] -> finish = 1;
      pthread_mutex_unlock(queue_mutex);
      return STATUS_NORMAL;
    }   

    // keep wait for the signal as the early one
    pthread_mutex_lock(&mutexes[dest]);
    pthread_mutex_unlock(queue_mutex);
    while(match == -1){
      pthread_cond_wait(&conditions[dest], &mutexes[dest]);
      pthread_mutex_unlock(&mutexes[dest]);
      pthread_mutex_lock(queue_mutex);
      pthread_mutex_lock(&mutexes[dest]);
      
      // check again when waked up
      match = search_match(src, dest, tag, RECV, src, 0);
      // release the lock if still not matched
      if(match == -1){
        pthread_mutex_unlock(queue_mutex);
      }
    }
    if(recvQueue[match] -> finish != 1){
      //printf("Thread[%d] do the writing as waked thread\n", src);
      len = recvQueue[match] -> len;
      request -> len = len;
      memcpy(recvQueue[match]->message, message, len*sizeof(char));
      
      sendQueue[self] -> finish = 1;
      recvQueue[match] -> finish = 1;
    }
    pthread_mutex_unlock(queue_mutex);
    pthread_mutex_unlock(&mutexes[dest]);
    return STATUS_NORMAL;
  }

  else if(type == RECV){

    //printf("Thread[%d] start looking for a send request\n", dest);
    // search for the send request
    self = search_match(src, dest, tag, RECV, dest, 1);
    // the other thread has finished the request for this one
    if(self == -1){
      if(src == ANY_SRC){
        match = search_match(src, dest, tag, SEND, dest, 0);
        request -> src = sendQueue[match] -> src;
      }
      if(tag == ANY_TAG){
        match = search_match(src, dest, tag, SEND, dest, 0);
        request -> tag = sendQueue[match] -> tag;
      }
      pthread_mutex_unlock(queue_mutex);
      return STATUS_NORMAL;
    }
    match = search_match(src, dest, tag, SEND, dest, 1);

    // do the finishing job as the late one
    if(match != -1){
      if(src == ANY_SRC){
        request -> src = sendQueue[match] -> src;
      }
      if(tag == ANY_TAG){
        request -> tag = sendQueue[match] -> tag;
      }
      //printf("Thread[%d] is the late recv, match is %d\n", dest, match);
      // do the message transmission
      memcpy(message, sendQueue[match]->message, len*sizeof(char));
      recvQueue[self] -> finish = 1;
      sendQueue[match] -> finish = 1;
      pthread_mutex_unlock(queue_mutex);
      return STATUS_NORMAL;
    }   
    // keep wait for the signal as the early one
    if(src != ANY_SRC){
      pthread_mutex_lock(&mutexes[src]);
    }
    pthread_mutex_unlock(queue_mutex);
    while(match == -1){
      if(src != ANY_SRC){
        pthread_cond_wait(&conditions[src], &mutexes[src]);
        pthread_mutex_unlock(&mutexes[src]);
        pthread_mutex_lock(queue_mutex);
        pthread_mutex_lock(&mutexes[src]);
        // check again if waked up
        match = search_match(src, dest, tag, SEND, dest, 0);
      }
      else{
        // busy wait if ANY_SRC is accepted
        usleep(5);
        sched_yield();
        pthread_mutex_lock(queue_mutex);
        match = search_match(src, dest, tag, SEND, dest, 0);
      }
      // release the lock if still not matched
      if(match == -1){
        pthread_mutex_unlock(queue_mutex);
      }
    }
    if(src == ANY_SRC){
      // update the received src
      request -> src = sendQueue[match] -> src;
    }
    if(tag == ANY_TAG){
      // update the received tag
      request -> tag = sendQueue[match] -> tag;
    }
    if(sendQueue[match] -> finish != 1){
      //printf("Thread[%d] do the writing as waked thread\n", dest);
      memcpy(message, sendQueue[match]->message, len*sizeof(char));
      recvQueue[self] -> finish = 1;
      sendQueue[match] -> finish = 1;
    }
    pthread_mutex_unlock(queue_mutex);
    if(src != ANY_SRC){
      pthread_mutex_unlock(&mutexes[src]);
    }
    return STATUS_NORMAL;
  }
  return STATUS_NORMAL;
}


/* msgWaitAll
 * purpose: Wait for an array of message request to complete
 * param:   requests - the array of message request to be completed
 */
int msgWaitAll( MSGRQST *requests, int nqrsts){
  int i;
  for(i=0; i<nqrsts; i++){
    msgWait(&requests[i]);
  }
  return STATUS_NORMAL;
}


/* getRank
 * purpose: Read the rank of a (completed) message request
 * param:   request - the request to be read.
 */
int getRank(MSGRQST *request){
  return request -> src;
}


/* getTag
 * purpose: Read the tag of a (completed) message request
 * param:   request - the request to be read.
 */
int getTag(MSGRQST *request){
  return request -> tag;
}


/* getLen
 * purpose: Read the length of a (completed) message request
 * param:   request - the request to be read.
 */
int getLen(MSGRQST *request){
  return request -> len;
}


/* search_match
 * purpose: Find the MSGRQST matching the given info on the message request queue
 * param:   src - the source of the message request
 *          dest- the destination of the message request
 *          tag - the tag of the message request
 *          type- the type of the message request, SEND or RECV
 *          id -  the rank of the calling thread
 *          finish - whether to avoid finished message, 0 as don't avoid, 1 as avoid
 * return:  the index of the matching message request on the relative queue
 *          -1 if not found
 */
int search_match(int src, int dest, int tag, int type, int id, int finish){
  int i;
  if(type == RECV){
    //printf("Thread[%d] enter search_match looking for recv\n", id);
    for(i=0; i<rqLength; i++){
      if(recvQueue[i]){
        if(recvQueue[i] -> dest == dest &&
         (recvQueue[i] -> src == src || recvQueue[i] -> src == ANY_SRC) &&
         (recvQueue[i] -> tag == tag || recvQueue[i] -> tag == ANY_TAG)){
          if(finish == 0){
            //printf("Thread[%d] search_match found recv (might be finished) \n", id);
            return i;
          }
          if(recvQueue[i] -> finish != 1){
            //printf("Thread[%d] search_match found recv\n", id);
            return i;
          }      
        }
      }
    }
    //printf("Thread[%d] looking for recv request no result\n", id);
    return -1;
  }
  else if(type == SEND){
    //printf("Thread[%d] enter search_match looking for send\n", id);
    for(i=0; i<sqLength; i++){
      //printf("send search i is %d\n", i);
      if(sendQueue[i]){
        if(sendQueue[i] -> dest == dest &&
         (sendQueue[i] -> src == src || src == ANY_SRC) &&
         (sendQueue[i] -> tag == tag || tag == ANY_TAG)){
          if(finish == 0){
            //printf("Thread[%d] search_match found send (might be finished)\n", id);
            return i;
          }
          if(sendQueue[i] -> finish != 1){
            //printf("Thread[%d] search_match found send\n", id);
            return i;      
          }
        }
      }
    }
    //printf("Thread[%d] looking for send request no result\n", id);
    return -1;
  }
  return -1;
}
