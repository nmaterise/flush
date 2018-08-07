INC = -I C:\Users\Bebotron\Eigen\Eigen -I C:\Users\Bebotron\Eigen\unsupported\Eigen

flush.exe: flush_main.o basic_funcs.o
	g++ -std=c++11 -openmp flush_main.o basic_funcs.o -o flush.exe

flush_main.o: flush_main.cpp
	g++ -c -std=c++11 -openmp $(INC) flush_main.cpp -o flush_main.o

basic_funcs.o: basic_funcs.cpp basic_funcs.h
	g++ -c -std=c++11 -openmp $(INC) basic_funcs.cpp -o basic_funcs.o

# flush.exe: flush_main.o basic_funcs.o
# 	icpc -std=c++11 -openmp flush_main.o basic_funcs.o -o flush.exe

# flush_main.o: ../flush_main.cpp
# 	icpc -c -std=c++11 -openmp ../flush_main.cpp

# basic_funcs.o: ../basic_funcs.cpp ../basic_funcs.h
# 	icpc -c -std=c++11 -openmp ../basic_funcs.cpp
