flush.exe: flush_main.o basic_funcs.o
	icpc -std=c++11 -openmp flush_main.o basic_funcs.o -o flush.exe

flush_main.o: ../flush_main.cpp
	icpc -c -std=c++11 -openmp ../flush_main.cpp

basic_funcs.o: ../basic_funcs.cpp ../basic_funcs.h
	icpc -c -std=c++11 -openmp ../basic_funcs.cpp
