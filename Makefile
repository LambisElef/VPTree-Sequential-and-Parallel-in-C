SHELL := /bin/bash
CC = gcc-7

vptree:
	cp inc/vptree.h src/vptree_sequential.c ./
	$(CC) -O2 -c vptree_sequential.c vptree.h
	ar -rcs vptree_sequential.a vptree_sequential.o
	mv vptree_sequential.a lib/
	rm vptree_sequential.c vptree_sequential.o

	cp inc/vptree.h src/vptree_pthreads.c ./
	$(CC) -O2 -c -pthread vptree_pthreads.c vptree.h
	ar -rcs vptree_pthreads.a vptree_pthreads.o
	mv vptree_pthreads.a lib/
	rm vptree_pthreads.c vptree_pthreads.o

	cp inc/vptree.h src/vptree_openmp.c ./
	$(CC) -O2 -c -fopenmp vptree_openmp.c vptree.h
	ar -rcs vptree_openmp.a vptree_openmp.o
	mv vptree_openmp.a lib/
	rm vptree_openmp.c vptree_openmp.o

	cp inc/vptree.h src/vptree_cilk.c ./
	$(CC) -O2 -c -fcilkplus vptree_cilk.c vptree.h
	ar -rcs vptree_cilk.a vptree_cilk.o
	mv vptree_cilk.a lib/
	rm vptree_cilk.c vptree_cilk.o
