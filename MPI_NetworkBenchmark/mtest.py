#!/usr/bin/env python
# Reads a matrix file of the format specified in PR3 for CSE160 
# Computes the maximum difference of results from code and numpy
# anything less than 1e-7 should be considered "equal"
# invoke as ./mtest.py <outputfile>
import  numpy
import sys
def convF(x):
	return float(x.strip())
def convD(x):
	return int(x.strip())

fname = sys.argv[1]
f = open(fname,"r")
N,M,P = map(convD, f.readline().strip().split(','))
print N,M,P
A=[]
B=[]
C=[]
for i in range(0,N):
	row = map(convF,f.readline().strip().split(','))
	A.append(row)
for i in range(0,M):
	row = map(convF,f.readline().strip().split(','))
	B.append(row)
for i in range(0,N):
	row = map(convF,f.readline().strip().split(','))
	C.append(row)
An=numpy.matrix(A)
Bn=numpy.matrix(B)
Cn=numpy.matrix(C)
diff = Cn-An*Bn
print "maximum difference", numpy.matrix.max(abs(diff))
