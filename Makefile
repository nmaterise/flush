exe: flush_main.o basic_funcs.o
	icpc flush_main.o basic_funcs.o

flush_main.o: flush_main.cpp basic_funcs.h 
	icpc -std=c++11 -openmp -I ~/Eigen/Eigen -I ~/Eigen/unsupported/Eigen flush_main.cpp basic_funcs.h

basic_funcs.o: basic_funcs.cpp basic_funcs.h
	icpc -std=c++11 -openmp -I ~/Eigen/Eigen -I ~/Eigen/unsupported/Eigen -C basic_funcs.cpp
