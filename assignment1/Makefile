CC=gcc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS = -O2


default:	all

V3:		triangles_library.o mmio.o src/V3_seq.c
	$(CC) $(CFLAGS) triangles_library.o mmio.o -o V3_seq.out src/V3_seq.c -no-pie

HADAMARD:	triangles_library.o mmio.o src/HADAMARD_seq.c
	$(CC) $(CFLAGS) triangles_library.o mmio.o -o HADAMARD_seq.out src/HADAMARD_seq.c -no-pie

V3_cilk:	triangles_library.o mmio.o src/V3_cilk.c
	$(CILKCC) $(CFLAGS) triangles_library.o mmio.o -o V3_cilk.out src/V3_cilk.c -fcilkplus

HADAMARD_cilk:		triangles_library.o mmio.o src/HADAMARD_cilk.c
	$(CILKCC) $(CFLAGS) triangles_library.o mmio.o -o HADAMARD_cilk.out src/HADAMARD_cilk.c -fcilkplus

V3_openmp:	triangles_library.o mmio.o src/V3_openmp.c
	$(CC) $(CFLAGS) triangles_library.o mmio.o -o V3_openmp.out src/V3_openmp.c -no-pie -fopenmp

HADAMARD_openmp:	triangles_library.o mmio.o src/HADAMARD_openmp.c
	$(CC) $(CFLAGS) triangles_library.o mmio.o -o HADAMARD_openmp.out src/HADAMARD_openmp.c -no-pie -fopenmp

V3_pthreads:	triangles_library.o mmio.o src/V3_pthread.c
	$(CC) $(CFLAGS) triangles_library.o mmio.o -o V3_pthreads.out src/V3_pthread.c -no-pie -lpthread -lm

HADAMARD_pthreads:	triangles_library.o mmio.o src/HADAMARD_pthread.c
	$(CC) $(CFLAGS) triangles_library.o mmio.o -o HADAMARD_pthreads.out src/HADAMARD_pthread.c -no-pie -lpthread -lm

mmio.o: lib/mmio.c lib/mmio.h
	$(CC) $(CFLAGS) -c lib/mmio.c

triangles_library.o:	lib/triangles_library.h lib/triangles_library.c
	$(CC) $(CFLGS) -c lib/triangles_library.c

.PHONY:	clean


all: V3 HADAMARD V3_cilk HADAMARD_cilk V3_openmp HADAMARD_openmp V3_pthreads HADAMARD_pthreads

clean:
	rm -f V3_seq.out HADAMARD_seq.out V3_cilk.out HADAMARD_cilk.out V3_openmp.out HADAMARD_openmp.out V3_pthreads.out HADAMARD_pthreads.out
