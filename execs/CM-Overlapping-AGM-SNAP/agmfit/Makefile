#
# Makefile for non-Microsoft compilers
#

## Linux  (uncomment the 2 lines below for compilation on Linux)
CXXFLAGS += -std=c++98 -Wall -fopenmp
LDFLAGS += -lrt -fopenmp
SNAPPATH = ../lib/

## CygWin (uncomment the 2 lines below for compilation on CygWin)
#CXXFLAGS += -Wall
#LDFLAGS += 

## Main application file
MAIN = agmfitmain

all: $(MAIN)

opt: CXXFLAGS += -O4 
opt: LDFLAGS += -O4
opt: $(MAIN)

# COMPILE
$(MAIN): $(MAIN).cpp Snap.o 
	g++ $(LDFLAGS) -o $(MAIN) $(MAIN).cpp $(SNAPPATH)agm/agm.cpp $(SNAPPATH)agm/agmfit.cpp Snap.o -I$(SNAPPATH)GLib -I$(SNAPPATH)snap -I$(SNAPPATH)agm 

Snap.o: 
	g++ -c $(CXXFLAGS) $(SNAPPATH)snap/Snap.cpp -I$(SNAPPATH)GLib -I$(SNAPPATH)snap

clean:
	rm -f *.o  $(MAIN)  $(MAIN).exe
	rm -rf Debug Release
