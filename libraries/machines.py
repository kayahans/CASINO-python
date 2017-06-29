
class Job:

    machine = None

    intel_compilers = set(['intel', 'icc', 'icpc', 'ifort'])

    def __init__(self, name = 'job', machine = None, cores = None, nodes = None, app = None, compiler = None):
        self.name = name
        self.machine = machine
        self.cores = cores
        self.nodes = nodes

        if compiler in self.intel_compilers:
            self.compiler = 'intel'


class Machines:
    """"""
    known_machines = ["cruller"],

    def __init__(self,scheduler, mpirun, declarations, cpu_per_node, max_node, max_time):


    @staticmethod
    def machine_init():

class Cruller:
    name = 'Cruller'

    def write_job_header(self, job):
        if job.queue is None:
            job.queue = 'batch'
        # end if
        c = '#!/bin/bash\n'
        c += '#PBS -N {0}\n'.format(job.name)
        c += '#PBS -o {0}\n'.format(job.outfile)
        c += '#PBS -e {0}\n'.format(job.errfile)
        c += '#PBS -l walltime={0}\n'.format(job.pbs_walltime())
        c += '#PBS -l nodes={0}\n'.format(job.nodes)
        c += '''
    echo $PBS_O_WORKDIR
    cd $PBS_O_WORKDIR
    '''
        return c



