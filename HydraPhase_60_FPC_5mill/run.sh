#!/bin/bash
# Number of MPI ranks
NP=2

# Launch with MPI, each rank gets a GPU based on its local rank
mpirun -np $NP bash -c 'export ACC_DEVICE_NUM=$OMPI_COMM_WORLD_LOCAL_RANK; ./a.out'

