# Mingcheng Zhu, zhumc11@gmail.com
CC = cc

default: mmmw

misc.o : misc.c misc.h
	mpicc -c misc.c

mm: mm.c
	gcc -O3 -g -Wall -o mm mm.c

mmmw: mmmw.c misc.o
	mpicc -g -Wall -o mmmw mmmw.c misc.o

msgbench: msgbench.c
	mpicc -g -Wall -o msgbench msgbench.c

clean:
	- /bin/rm mm
	- /bin/rm mmmw
	- /bin/rm msgbench
	- /bin/rm misc.o
