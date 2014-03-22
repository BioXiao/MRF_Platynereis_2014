all:
	gcc EM.h EM.c -o EM -pedantic -Wall -lm
	
windows:
	gcc EM.h EM.c -o EM -pedantic -Wall -lm -ansi	
