CC = g++ -std=c++11 -Wpedantic -O3 -I ../nlohmann
IBM = ../vaccine_abm/*.cpp


run_build:
	$(CC) -g -c $(IBM)
	$(CC) -g -o Run main_national_interim_report_ALTERED.cpp *.o
	# $(CC) -g -o Run main_national_initcond.cpp *.o
	# $(CC) -g -o Run main_national_interim_report.cpp *.o

clean: 
	rm  *.o Run