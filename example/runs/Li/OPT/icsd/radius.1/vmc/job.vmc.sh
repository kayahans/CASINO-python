#!/bin/bash
#PBS -N vmc
#PBS -o vmc.o
#PBS -e vmc.e
#PBS -l walltime=12:00:00
#PBS -l nodes=8:ppn=12

module purge
module load mvapich2-2.0/intel

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/Compiler/11.1/072/lib/intel64:/opt/intel/mkl/10.2.5.035/lib/em64t
source /opt/intel/Compiler/11.1/072/bin/ifortvars.sh intel64
    
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
mpirun_rsh -np 96 -hostfile $PBS_NODEFILE casino