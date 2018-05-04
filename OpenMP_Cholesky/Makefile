# Mingcheng Zhu, zhumc11@gmail.com

CC = gcc
CFLAGS = -g

default: blockCholeskyMP

cholesky: cholesky.c cholesky.o
	$(CC) $(CFLAGS) -fopenmp -o cholesky -O3 cholesky.c -lm 

blockCholesky: blockCholesky.c blockCholesky.o
	$(CC) $(CFLAGS) -fopenmp -o blockCholesky -O3 blockCholesky.c -lm 

blockCholeskyMP: blockCholeskyMP.c blockCholeskyMP.o
	$(CC) $(CFLAGS) -fopenmp -o blockCholeskyMP -O3 blockCholeskyMP.c -lm 

clean:
	-/bin/rm *o
