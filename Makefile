# Linear-algebra library
EIGEN = /usr/include/eigen3/

CXX          = g++
CLINKER     = g++

CXXFLAGS      = -I $(EIGEN) -I. -fopenmp -O3 -fdiagnostics-color=always -std=c++17
#LIBS        = -lm
DEPEND= makedepend

SRC1        =main.cpp
EXECS1      =main

main:$(OBJS1)
	$(CLINKER) -g $(SRC1) $(OPTFLAGS) -o $(EXECS1) $(CXXFLAGS)

clean:
	/bin/rm -f *.o *~ $(EXECS1)

.c.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp