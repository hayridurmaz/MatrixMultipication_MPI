make: rowmv.c
	mpicc -fopenmp -g -ggdb -Wall -lm -o rowmv rowmv.c 