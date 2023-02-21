all: kd-tree.c kd-tree_MPI.c 
	gcc -g -O1 -c util_functions.c -o util_functions.o
	gcc -g -O1 -fopenmp -c  -std=c99 kd-tree.c -o kd-tree.o
	gcc -g -O1 -fno-omit-frame-pointer -fopenmp -std=c99 util_functions.o kd-tree.o -o kd-tree_OMP.x -lm
	gcc -g -O1 -c -std=c99 util_functions.c -o util_functions.o
	gcc -g -O1 -c -std=c99 kd-tree.c -o kd-tree.o
	gcc -g -O1 -fno-omit-frame-pointer -std=c99 util_functions.o kd-tree.o -o kd-tree_serial.x -lm
	gcc -g -O1 -c -std=c99 util_functions.c -o util_functions.o
	mpicc -g -O1 -c -std=c99 kd-tree_MPI.c -o kd-tree.o
	mpicc -g -O1 -fno-omit-frame-pointer -std=c99 util_functions.o kd-tree.o -o kd-tree_MPI.x -lm
	gcc -g -O1 -std=c99 -fopenmp data_sample.c -o data_sample.x

kdtree: kd-tree.c util_functions.c util_functions.h
	gcc -g -Wall -O1 -c util_functions.c -o util_functions.o
	gcc -g -Wall -O1 -fopenmp -c  -std=c99 kd-tree.c -o kd-tree.o
	gcc -g -Wall -O1 -fno-omit-frame-pointer -fopenmp -std=c99 util_functions.o kd-tree.o -o kd-tree_OMP.x -lm

kdtrees: kd-tree.c util_functions.c util_functions.h
	gcc -g -Wall -c -O1 -std=c99 util_functions.c -o util_functions.o
	gcc -g -Wall -c -O1-std=c99 kd-tree.c -o kd-tree.o
	gcc -g -Wall -O1 -fno-omit-frame-pointer -std=c99 util_functions.o kd-tree.o -o kd-tree_serial.x -lm

kdtreempi: kd-tree_MPI.c util_functions.c util_functions.h
	gcc -g -c -O1 -std=c99 util_functions.c -o util_functions.o
	mpicc -g -O1 -c -std=c99 kd-tree_MPI.c -o kd-tree.o
	mpicc -g -O1 -fno-omit-frame-pointer -std=c99 util_functions.o kd-tree.o -o kd-tree_MPI.x -lm

clean: 
	rm *.o
	rm *.x

datasample: data_sample.c
	gcc -g -Wall -std=c99 -fopenmp data_sample.c -o data_sample.x

benchmark: 
	gcc -g -O1 -c util_functions.c -o util_functions.o
	gcc -g -O1 -fopenmp -c -std=c99 kd-tree_test.c -o kd-tree.o
	gcc -g -O1 -fno-omit-frame-pointer -fopenmp -std=c99 util_functions.o kd-tree.o -o kd-tree_OMP_test.x -lm
	gcc -g -O1 -c -std=c99 util_functions.c -o util_functions.o
	gcc -g -O1 -c -std=c99 kd-tree_test.c -o kd-tree.o
	gcc -g -O1 -fno-omit-frame-pointer -std=c99 util_functions.o kd-tree.o -o kd-tree_serial_test.x -lm
	gcc -g -c -O1 -std=c99 util_functions.c -o util_functions.o
	mpicc -g -O1 -c -std=c99 kd-tree_MPI_test.c -o kd-tree.o
	mpicc -g -O1 -fno-omit-frame-pointer -std=c99 util_functions.o kd-tree.o -o kd-tree_MPI_test.x -lm

cleancsv:
	rm benchmarks/*.csv
	
	
