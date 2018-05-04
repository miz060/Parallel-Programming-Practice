#include <pthread.h>
/** 
 * File name: mp.h
 * File description: Definitions for POSIX thread Message Passing
 * Name:  Mingcheng Zhu
 * Email: zhumc11@gmail.com
 * Date:  Mar 4, 2018
 */

#ifndef __CS160MP
#define __CS160MP
#define ANY_TAG -2
#define ANY_SRC -2
#define MPTRUE 1
#define MPFALSE 0

/* definition of the message request struct */
typedef struct
{
	int src;
  int dest;
  int tag;
  int len;
  int finish;
  int write;
  int type;
  void* message;
} MSGRQST; 
static MSGRQST MSGRQST_NULL = {-1};

/* you MUST NOT change the definitions of API below */
/* MP_Init
 * purpose: Initialize any internal data structures needed to implement message passing.
 *          Called once per program.
 * param:   nthreads - number of threads created
 */
int MP_Init(int nthreads);


/* MP_Size
 * purpose: Return the number of threads passed to MP_Init
 */
int MP_Size();


/* MP_Rank
 * purpose: Give a rank (0 .. MP_Size()) to the thread with a particular id
 * return:  the given thread rank
 */
int MP_Rank(pthread_t thread);


/* MP_Finalize
 * purpose: Clean up all allocated memory. 
 *          Indicates no more message passing will be performed
 *          Called once after all message passing is complete
 */
int MP_Finalize();


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
int iSend(void *buf, int len, int src, int dest, int tag, MSGRQST * request);


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
int iRecv(void *buf, int len, int src, int dest, int tag, MSGRQST * request);


/* msgWait
 * purpose: Wait for a specific message to complete
 * param:   request - the request to be waited.
 */
int msgWait( MSGRQST * request);


/* msgWaitAll
 * purpose: Wait for an array of message request to complete
 * param:   requests - the array of message request to be completed
 */
int msgWaitAll( MSGRQST *requests, int nqrsts);


/* getRank
 * purpose: Read the rank of a (completed) message request
 * param:   request - the request to be read.
 */
int getRank(MSGRQST *request);


/* getTag
 * purpose: Read the tag of a (completed) message request
 * param:   request - the request to be read.
 */
int getTag(MSGRQST *request);


/* getLen
 * purpose: Read the length of a (completed) message request
 * param:   request - the request to be read.
 */
int getLen(MSGRQST *request);

#endif
