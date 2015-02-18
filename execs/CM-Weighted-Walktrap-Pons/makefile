CC=g++
CFLAGS= -Wall -O3

walktrap : walktrap.o communities.o graph.o heap.o
	$(CC) -o $@ $^ $(CFLAGS)

all : walktrap

%.o : %.cpp
	$(CC) -c $< $(CFLAGS) 
clean:
	rm *.o
	

