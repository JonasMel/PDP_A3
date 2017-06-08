#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 10:00
#SBATCH -p node -n 64

module load gcc openmpi
make loopy


mpirun -np 1 ./loopy 1 100000000 10 &&

mpirun -np 2 ./loopy 2 100000000 10 &&

mpirun -np 4 ./loopy 4 100000000 10 &&

mpirun -np 8 ./loopy 8 100000000 10 &&

mpirun -np 16 ./loopy 16 100000000 10 &&

mpirun -np 32 ./loopy 32 100000000 10 &&

mpirun -np 64 ./loopy 64 100000000 10
