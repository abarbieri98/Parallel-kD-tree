#!/bin/bash

# 1 = n_points
# 2 = n_processes / n_threads

echo "-- Serial Time --"
./kd-tree_serial.x $1 $2
echo "-- Parallel Time (OMP)--"
./kd-tree_OMP.x $1 $2 
echo "-- Parallel Time (MPI)"
mpirun -np $2 kd-tree_MPI.x $1 2>/dev/null
