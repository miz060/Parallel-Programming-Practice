# Mingcheng Zhu, zhumc11@gmail.com

CC = gcc
CFLAGS = -g

default: libmp.a

mp.o: mp.c mp.h
	$(CC)  $(CFLAGS) -c mp.c -pthread

libmp.a: mp.o
	ar rcs libmp.a mp.o

msgbench: msgbench.c libmp.a
	$(CC)  $(CFLAGS) -o msgbench msgbench.c libmp.a -lpthread

test: test.c libmp.a
	$(CC)  $(CFLAGS) -o test test.c libmp.a -lpthread

clean:
	-/bin/rm *o
	-/bin/rm libmp.a
	
