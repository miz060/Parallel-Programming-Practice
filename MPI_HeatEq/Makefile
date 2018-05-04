# Mingcheng Zhu, zhumc11@gmail.com
CC = cc

default: heat2d

svalidate.o : svalidate.c svalidate.h
	$(CC) -c svalidate.c

misc.o : misc.c misc.h svalidate.o
	$(CC) -c misc.c svalidate.o -lm

heat2drow: heat2drow.c svalidate.o misc.o
	mpicc -g -Wall -o heat2drow heat2drow.c svalidate.o misc.o -lm

heat2dcol: heat2dcol.c svalidate.o misc.o heat2drow
	mpicc -g -Wall -o heat2dcol heat2dcol.c svalidate.o misc.o -lm

heat2d: heat2d.c svalidate.o misc.o heat2drow heat2dcol
	$(CC) -g -Wall -o heat2d heat2d.c svalidate.o misc.o -lm

clean:
	- /bin/rm svalidate.o
	- /bin/rm misc.o
	- /bin/rm heat2d
	- /bin/rm heat2drow
	- /bin/rm heat2dcol
