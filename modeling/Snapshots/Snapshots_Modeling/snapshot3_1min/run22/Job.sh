#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
#$ -N test1
#$ -pe smp 16
#$ -l h_rt=48:00:00

module load Sali
module load mpi/openmpi-x86_64
module load imp
module load python3/scikit/0.21.3

mpirun -np $NSLOTS python3 static_snapshot.py 3 1min
