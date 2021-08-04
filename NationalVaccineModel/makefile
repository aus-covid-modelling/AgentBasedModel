CC = g++ -std=c++11 -Wpedantic -O3
IBM = ../vaccine_abm/*.cpp


run_build:
	$(CC) -g -c $(IBM)
	$(CC) -g -o Run main_national.cpp *.o

clean: 
	rm  *.o Run