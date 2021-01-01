make: rowmv.c
  ifeq ($(DEBUG),1)
	gcc -fopenmp -g -ggdb -lm -DDEBUG=1 -o rowmv rowmv.c 
  else
	gcc -fopenmp -g -ggdb -lm -o rowmv rowmv.c 
  endif