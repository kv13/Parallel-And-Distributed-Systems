SHELL := /bin/bash
CC = gcc-7
MPICC  = mpicc

.PHONY: lib

.PHONY: clean

lib:
	mkdir -p lib
	cd src; $(CC) -c knnring_sequential.c -lopenblas -lm; cd ..
	cd src; ar rcs ../lib/knnring_sequential.a knnring_sequential.o; cd ..

v2_rand:
	cd src; $(MPICC) knnring_mpi.c vptree_ring.c v2_mpi_interface.c -lm -lopenblas -O2 -o v2; cd ..

v2:
	cd src; $(MPICC) knnring_mpi.c vptree_ring.c readlib.c v2_mpi2_2.c -lm -lopenblas -O2 -o v2_2; cd ..

v1_rand:
	cd src; $(MPICC) knnring_mpi.c v1_mpi2.c -lm -lopenblas -O2 -o v1; cd ..

v1:
	cd src; $(MPICC) knnring_mpi.c readlib.c v1_mpi2_2.c -lm -lopenblas -O2 -o v1_2; cd ..


all: v1 v1_rand v2 v2_rand

clean:
	cd src; rm -f v1 v2 v1_2 v2_2



