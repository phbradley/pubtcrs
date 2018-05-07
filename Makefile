# the compiler
CC = g++

## compiler flags:

# debugging
#CCFLAGS  = -g -std=c++11 -Wall

# production
CCFLAGS  = -O3 -std=c++11 -Wall

# EDIT THIS to point to the location where you installed boost
# this should be the main directory containing boost/ doc/ libs/ etc as subdirectories
BOOSTDIR = /home/pbradley/include/boost_1_67_0

# 
INCLUDES = -I ./include/  -I $(BOOSTDIR)

# recompile if any .hh files changed
HHS = src/*.hh

all: bin/pgen bin/neighbors bin/tcrdists bin/correlations

bin/pgen:  src/pgen.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o bin/pgen src/pgen.cc

bin/neighbors: src/neighbors.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o bin/neighbors src/neighbors.cc

bin/tcrdists:  src/tcrdists.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o bin/tcrdists src/tcrdists.cc

bin/correlations:  src/correlations.cc  $(HHS)
	$(CC) $(CCFLAGS) $(INCLUDES) -o bin/correlations src/correlations.cc

clean:
	-rm bin/*
