CC = g++
CPPFLAGS = -g -Wall -O2
LDFLAGS =

CPX_INCDIR  = /opt/CPLEX_Studio/cplex/include
CPX_LIBDIR  = /opt/CPLEX_Studio/cplex/lib/x86-64_linux/static_pic
CPX_LDFLAGS = -lcplex -lm -pthread

OBJ = main.o 

%.o: %.cpp
		$(CC) $(CPPFLAGS) -I$(CPX_INCDIR) -c $^ -o $@

main: $(OBJ)
		$(CC) $(CPPFLAGS) $(OBJ) -o main -L$(CPX_LIBDIR) $(CPX_LDFLAGS)

copy:
		cp driller-extended.cpp main.cpp

clean:
		rm -rf $(OBJ) main

.PHONY: clean
