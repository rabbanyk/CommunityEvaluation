###### LOCAL CUSTOMIZATION ####################

## c compiler
CC = g++

## c compiler flags
CFLAGS = -Wall -O3 -DNDEBUG
#CFLAGS = -Wall -g

# arpack, lapack, blas libraries and where to find them 
ARPACK_LIB = -larpack
LAPACK_LIB = -llapack
BLAS_LIB = -lblas
LDFLAGS = -L/usr/lib/sse2/atlas

###### END OF CUSTOMIZATION ###################

LDLIBS = ${ARPACK_LIB} ${LAPACK_LIB} ${BLAS_LIB}

OBJS = commfind.o eig.o
SRC = commfind.c eig.c
HDRS = commfind.h
DISTFILES = README COPYING INSTALL CHANGELOG example.net

commfind: ${OBJS}

${OBJS}: ${HDRS}

.PHONY: clean tar dist

clean:
	rm -f *.o *~ core
tar:
	tar zcvf commfind.tgz Makefile ${SRC} ${HDRS}

dist :
	tar zcvf commfind.tgz Makefile ${SRC} ${HDRS} ${DISTFILES}
