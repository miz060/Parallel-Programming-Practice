# Makefile for CSE160 PA5
# Mingcheng Zhu, A92047564, miz060@ucsd.edu

CC = gcc
CFLAGS = -g

default: libcs160mp.a

cs160mp.o: cs160mp.c cs160mp.h
	$(CC)  $(CFLAGS) -c cs160mp.c -pthread

libcs160mp.a: cs160mp.o
	ar rcs libcs160mp.a cs160mp.o

msgbench: msgbench.c libcs160mp.a
	$(CC)  $(CFLAGS) -o msgbench msgbench.c libcs160mp.a -lpthread

test: test.c libcs160mp.a
	$(CC)  $(CFLAGS) -o test test.c libcs160mp.a -lpthread

clean:
	-/bin/rm *o
	-/bin/rm libcs160mp.a
	
