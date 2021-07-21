CC = g++ -std=c++11 -Wpedantic
IBM = ../vaccine_abm/*.cpp


run_build:
	$(CC) -c $(IBM)
	$(CC) -o Run main_national.cpp *.o

clean: 
	rm  *.o Run