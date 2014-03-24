all:
	gcc EM.c -o EM -pedantic -Wall -lm
	
windows:
	gcc EM.c -o EM -pedantic -Wall -lm -ansi
