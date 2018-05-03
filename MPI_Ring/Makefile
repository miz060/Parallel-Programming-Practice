# Sample Makefile for CS160
CC = cc

default: svtest 

svalidate.o : svalidate.c svalidate.h
	$(CC) -c svalidate.c

mpi-ring2: mpi-ring2.c svalidate.o
	mpicc -g -Wall -o mpi-ring2 mpi-ring2.c svalidate.o -lm

svtest: svtest.c svalidate.o mpi-ring2
	$(CC) -o svtest svtest.c svalidate.o

clean:
	- /bin/rm svalidate.o
	- /bin/rm svtest 
	- /bin/rm mpi-ring2 
