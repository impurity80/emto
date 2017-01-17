#!/bin/bash
#
#$ -q all.q
#$ -pe mpi_32 32 
#$ -N step-1.585
#$ -M impurity@oz
#$ -m bea
#$ -o run.out
#$ -e run.err
#$ -V
#$ -cwd
#

mpirun -np 32 python step0.py
