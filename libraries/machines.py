from casino_python import settings
import copy
import re

class Job:

    def __init__(self, name = 'job', cores = None, nodes = None, app = None, time=24):
        self.name = str(name)
        self.machine = Machines(settings.machinename)
        self.outfile= self.name+'.o'
        self.errfile= self.name+'.e'
        self.time=walltime(time)
        self.app = app
        self.run_command=self.machine.run_command
        if cores is None and nodes is None:
            raise ValueError('Number of nodes or cores must be provided in the Job declaration')
        if cores is None:
            self.nodes=nodes
            self.cores=self.machine.cpu_per_node*self.nodes
        if nodes is None:
            self.cores=cores
            self.nodes=cores/self.machine.cpu_per_node
        if app is None:
            raise ValueError('app must be provided in Job declaration')
        self.script=write_job(self, app)


def walltime(time):

    if re.search(r"([0-9][0-9]:[0-5][0-9]:[0-5][0-9])", str(time)):
        return str(time)
    elif re.search("([0-9][0-9])", str(time)):
        return str(time)+':00:00'
    else:
        raise ValueError('Time format is incorrect')

def write_job(job,app):

    assert isinstance(job, Job)

    if job.machine.scheduler == 'pbs':
        c = '#!/bin/bash\n'
        c += '#PBS -N {0}\n'.format(job.name)
        c += '#PBS -o {0}\n'.format(job.outfile)
        c += '#PBS -e {0}\n'.format(job.errfile)
        c += '#PBS -l walltime={0}\n'.format(job.time)
        c += '#PBS -l nodes={0}:ppn={1}\n'.format(job.nodes, job.machine.cpu_per_node)
        c += '{0}'.format(job.machine.declarations)
        c += '''
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
'''
        c += '{0} -np {1} -hostfile $PBS_NODEFILE {2}'.format(job.machine.mpirun, job.cores, app)
        return c



#################
# Define new computers here
# Script is created
#
#

cruller=dict(
    name='cruller',
    scheduler='pbs',
    mpirun='mpirun_rsh',
    declarations='''
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
    ''',
    cpu_per_node=12,
    max_node=38,
    max_time=96,
    run_command='qsub'
)

known_machines = [cruller]

class Machines:
    """"""
    def __init__(self,name=None):
        for machine in known_machines:
            if name == machine['name']:
                self.scheduler = machine['scheduler']
                self.mpirun = machine['mpirun']
                self.declarations = machine['declarations']
                self.cpu_per_node = machine['cpu_per_node']
                self.max_node = machine['max_node']
                self.max_time = machine['max_time']
                self.run_command = machine['run_command']

class Cruller:
    name = 'Cruller'





