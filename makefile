make: rowmv.c
  ifeq ($(DEBUG),1)
	gcc -fopenmp -g -Wall -Wextra -ggdb -lm -DDEBUG=1 -o rowmv rowmv.c 
  else
	gcc -fopenmp -g -Wall -Wextra -ggdb -lm -o rowmv rowmv.c 
  endif