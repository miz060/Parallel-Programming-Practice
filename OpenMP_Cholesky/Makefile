# Makefile for CSE160 P6
# Mingcheng Zhu, A92047564, miz060@ucsd.edu

CC = gcc
CFLAGS = -g

default: blockCholeskyMP

cholesky: cholesky.c cholesky.o
	$(CC) $(CFLAGS) -fopenmp -o cholesky -O3 cholesky.c cs160validate.c -lm 

blockCholesky: blockCholesky.c blockCholesky.o cs160validate
	$(CC) $(CFLAGS) -fopenmp -o blockCholesky -O3 blockCholesky.c cs160validate.c -lm 

blockCholeskyMP: blockCholeskyMP.c blockCholeskyMP.o cs160validate
	$(CC) $(CFLAGS) -fopenmp -o blockCholeskyMP -O3 blockCholeskyMP.c cs160validate.c -lm 

cs160validate: cs160validate.c cs160validate.h
	$(CC) $(CFLAGS) -c -O3 cs160validate.c -lm 

#choleskyS: choleskyS.c choleskyS.o
#	$(CC) $(CFLAGS) -o choleskyS -O3 choleskyS.c cs160validate.c -lm 

clean:
	-/bin/rm *o
