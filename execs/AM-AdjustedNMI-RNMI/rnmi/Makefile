CC=g++
#CFLAGS= -O3 -pg -Wall
UNAME= $(shell uname)
ifeq ($(UNAME), Linux)
	CFLAGS=-O3
#	CFLAGS= 
endif
ifeq ($(UNAME), Darwin)
#	CFLAGS=-fast
	CFLAGS=-O3
endif
MYBINDIR=.
OFILE := zrg.o
OUTPUT := nmi rnmi
all: ${OUTPUT} ${OFILE}

${OUTPUT}: ${OFILE}
	${CC} -o ${MYBINDIR}/$@ ${CFLAGS} $@.cpp ${OFILE}
%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@
clean:
	rm -f *.o
