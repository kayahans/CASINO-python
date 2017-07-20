#!/bin/bash
#PBS -N dft
#PBS -o dft.o
#PBS -e dft.e
#PBS -l walltime=12:00:00
#PBS -l nodes=2:ppn=12

module purge
module load mvapich2-2.0/intel

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/Compiler/11.1/072/lib/intel64:/opt/intel/mkl/10.2.5.035/lib/em64t
source /opt/intel/Compiler/11.1/072/bin/ifortvars.sh intel64
if [ -f /usr/share/Modules/init/bash ]; then
        source /usr/share/Modules/init/bash
else
        echo "Could not source environment modules!"
        exit 1
fi
    
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
mpirun_rsh -np 24 -hostfile $PBS_NODEFILE /home/apps/espresso-5.0.3/bin/pw.x -pw2casino < input.qe > output.qe