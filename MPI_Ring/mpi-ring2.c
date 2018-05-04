/** 
 * File name: mpi-ring2.c
 * File description: practice for mpi package
 * Name:  Mingcheng Zhu
 * Email: zhumc11@gmail.com
 * Date:  Jan 19, 2018
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include "svalidate.h"

#define PROCID_P 1
#define SEEDVAL_P 2
#define S_TAG 0
#define P_TAG 1


/* main
 * purpose: main function for executable mpi-ring2
 * param:   argc - cmd argument counts, argv - cmd argument array
 * return:  exit status, 0 if normal, non-zero if error occurs
 */
int main(int argc, char* argv[]){

	int debug = 0;
	float power = 1.2;

	MPI_Status status;
	int comm_sz;
	int my_rank;
	int pre_rank;
	int next_rank;
	
	// set up MPI components
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// special case: numproc = 1
	if(comm_sz==1){
		MPI_Finalize();
		return 0;
	}

	// validate command line arguments
	// check arugument counts
	if(argc != 3){
		if(my_rank == 0){
			fprintf(stderr, "usage: mpi-ring2 <procid> <seedval>\n");
		}
		MPI_Finalize();
		return -1;
	}
	char* pid_arg = trim(argv[PROCID_P]);
	char* seed_arg = trim(argv[SEEDVAL_P]);

	// check if procid can be interpreted as an int
	// check if seedval can be interpreted as a float
	if(!isInteger(pid_arg) || !isFloat(seed_arg)){
		if(my_rank == 0){
			fprintf(stderr, "usage: mpi-ring2 <procid> <seedval>\n");
		}
		MPI_Finalize();
		return -1;
	}
	float seedval = atof(seed_arg);

	// check if seedval>=0
	if(seedval<0){
		if(my_rank == 0){
			fprintf(stderr, "usage: mpi-ring2 <procid> <seedval>\n");
		}
		MPI_Finalize();
		return -1;
	}
	// command line arguments validation ends

	int start_id = atoi(pid_arg);
	// procid>=0 case: flow through the ring in clockwise fashion
	if(start_id>=0){
		// normalize start_id
		while(start_id >= comm_sz){
			start_id -= comm_sz;
		}
		pre_rank = (comm_sz+my_rank-1)%comm_sz;
		next_rank = (my_rank+1)%comm_sz;
	}
	// procid<0 case: flow through the ring in counter clockwise fashion
	else{
		// normalize start_id
		while(start_id<0){
			start_id += comm_sz;
		}
		pre_rank = (my_rank+1)%comm_sz;
		next_rank = (comm_sz+my_rank-1)%comm_sz;
	}
	if(debug){
		printf("the start id is %d\n", start_id);
	}

	float message;
	float value;
	// start_id process
	if(my_rank == start_id){
		
		// send seed value to next process
		MPI_Send(&seedval, 1, MPI_FLOAT, next_rank, S_TAG, MPI_COMM_WORLD);

		// receive calculated value from previous process
		MPI_Recv(&message, 1, MPI_FLOAT, pre_rank, S_TAG, MPI_COMM_WORLD, &status);
		if(debug){
			printf("process [%d] receive float: %f\n", my_rank, message);
		}
		value = powf(message+(float)my_rank, power);
		if(debug){
			printf("process [%d] calculate value: %f\n", my_rank, value);
		}

		// send the final result to process 0 for printing
		MPI_Send(&value, 1, MPI_FLOAT, 0, P_TAG, MPI_COMM_WORLD);
	}

	// inner-ring process
	else{
		// receive calculated value from previous process
		MPI_Recv(&message, 1, MPI_FLOAT, pre_rank, S_TAG, MPI_COMM_WORLD, &status);
		if(debug){
			printf("process [%d] receive float: %f\n", my_rank, message);
		}
		value = powf(message+(float)my_rank, power);
		if(debug){
			printf("process [%d] calculate value: %f\n", my_rank, value);
		}

		// send the calculated value to next process
		MPI_Send(&value, 1, MPI_FLOAT, next_rank, S_TAG, MPI_COMM_WORLD);
	}

	// process 0 printing
	if(my_rank == 0){
		printf("start %d seed %f nproc %d\n", start_id, seedval, comm_sz);
		float result;
		// receive final result from start_id process
		MPI_Recv(&result, 1, MPI_FLOAT, start_id, P_TAG, MPI_COMM_WORLD, &status);
		printf("%f\n", result);
	}

	MPI_Finalize();
	return 0;
}	
