# A parallel implementation of k-D tree

This repository the code implemented for Assignment 2 of the course of Foundations of High Performance Computing.

## A general overview

A k-D tree is a particular type of binary tree proposed by Friedman et al. in 1976 ([original paper](https://www.researchgate.net/publication/220493118_An_Algorithm_for_Finding_Best_Matches_in_Logarithmic_Expected_Time)) as a solution for the nearest neighbor problem, where each node of the tree represents a partition of a given input dataset. 
The aim of the project is to implement such algorithm using a parallel approach to the computation and exploiting both shared and distributed memory.

## Folder content

The file present in the folder are the following:
- `data_sample.c`: source code for simple synthetic data generation according to the constraints included in `Report`
- `kd-tree.c`: source code of the serial and OpenMP implementations
- `kd-tree_MPI.c`: source code for the MPI implementation
- `Makefile` (see below)
- `Report`: report of the project including analysis of the algorithms used, details of the implementations and benchmark results.
- `test.sh`: bash script to test the three implementations
- `util_functions.c` and `util_functions.h`: library used by all the implementations

## Makefile

The Makefile is used to properly compile the sources, the rules implemented are:
- `make all`: compile all the files present in the folder
- `make <serial/omp/mpi>`: compile the selected version of the program
- `make datasample`: compile the data generator code
- `make clean`: remove all the objects and executable files

## Basic usage

First, to generate a dataset using the provided code the syntax is
```bash
./data_sample.x [number of points to generate] [number of dimensions fo each point] > dataset.csv
```
To run one of the implementations once the dataset is generated one can use
- Serial
``` bash
./kd-tree_serial.x [number of points] [number of threads]
```
For the serial implementation the number of threads will not have any effect.

- OpenMP
``` bash
./kd-tree_OMP.x [number of points] [number of threads]
```

- MPI
``` bash
mpirun -np [number of processes] kd-tree_MPI.x [number of points]
```
As explained in the report, MPI will only allow the use of number of processes in powers of two (2,4,8, ..., $2^N$).

Furthermore, if all the three implementations are compiled one can use the `test.sh` script to test them, using the syntax

```bash
./test.sh [number of points] [number of threads / processes]
```





