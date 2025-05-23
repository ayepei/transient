#!/bin/bash
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -J transient-28core-test 
#SBATCH --exclusive


# cd /work1/wangjh/yumiao/xuerong/OpenMOC-transient_update/transient_update/build
# make && ./run/OPENMOC-transient
cd /home/jrkkkkkk/hpc/transient_update/build
make && ./run/OPENMOC-transient