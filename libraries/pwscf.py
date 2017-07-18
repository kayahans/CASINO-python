import sys,os,copy
from collections import OrderedDict
from units import atomic_weight
import fnmatch
from machines import Job
from casino_python import settings,System
from error_handler import error, warning

default_input_control = OrderedDict(
    header='&CONTROL',
    title='generic.qe',
    calculation='scf',
    pseudo_dir='./',
    wf_collect='True'
)

default_input_system = OrderedDict(
    header='&SYSTEM',
    ecutwfc=None,
    input_dft=None,
    nspin=2,
    occupations='smearing',
    degauss=0.02,
    starting_magnetization=0.7,
    nosym='True',
    noinv='True',
    nosym_evc='True',
    ibrav=0,
    ntyp=None,
    nat=None
)

default_input_electrons = OrderedDict(
    header='&ELECTRONS',
    startingwfc='random',
    electron_maxstep=500
)

default_cell_params = OrderedDict(
    header='CELL_PARAMETERS angstrom',
    lattice_params = None
)

default_atomic_species = OrderedDict(
    header='ATOMIC_SPECIES',
    atoms_weights_psps = None
)

default_atomic_positions = OrderedDict(
    header='ATOMIC_POSITIONS angstrom',
    cart_coords = None
)

default_k_points = OrderedDict(
    header = 'K_POINTS crystal',
    num_kpoints = 0,
    kgrid = None
)

class Pwscf:

    def __init__(self, system=None, job=None, **kwargs):

        assert isinstance(system, System)
        assert isinstance(job, Job)
        #Inputs
        self.system = system
        self.job=Job(name=job.name,cores=job.cores,time=job.time,app=job.app+' < input.qe > output.qe')
        self.kwargs = kwargs

        #Class derived parameters
        self.rundir = system.rundir + '/qe'
        self.inputfile = self.rundir + '/input.qe'
        self.outputfile = self.rundir + '/output.qe'
        self.complete = False

        #Prepare input parameters - new attr
        self.process_system_inputs()
        self.process_pwscf_inputs()

        #Check if calculation already finished
        self.complete=False
        self.running=False
        self.control_pwscf()

        # If complete, do nothing; else write input.qe
        if not self.running and not self.complete:
            self.write_pwscf_inputs()
            self.execute()

    def process_pwscf_inputs(self):

        self.input_control = default_input_control.copy()
        self.input_system = default_input_system.copy()
        self.input_electrons = default_input_electrons.copy()
        self.input_cell_params = default_cell_params.copy()
        self.input_atomic_species = default_atomic_species.copy()
        self.input_atomic_positions = default_atomic_positions.copy()
        self.input_k_points = default_k_points.copy()

        for key, value in self.kwargs.iteritems():
            if key in self.input_control:
                self.input_control[key] = value
            if key in self.input_system:
                self.input_system[key] = value
            if key in self.input_electrons:
                self.input_electrons[key] = value
            if key in self.input_cell_params:
                self.input_cell_params[key] = value
            if key in self.input_atomic_species:
                self.input_atomic_species[key] = value
            if key in self.input_atomic_positions:
                self.input_atomic_positions[key] = value
            if key in self.input_k_points:
                self.input_k_points[key] = value

        del self.kwargs

    def process_system_inputs(self):

        system = self.system
        structure=self.system.structure
        unique_atoms = set(structure.species)
        natoms = len(structure.species)

        awp = []
        max_ecut = -1

        for specie in unique_atoms:  # Check here later
            if not os.listdir(system.runpspdir) == []:
                for file in os.listdir(system.runpspdir):

                    if fnmatch.fnmatch(file, specie + '*'):
                        awp.append([str(specie), atomic_weight(specie), file])
                    else:
                        awp.append([str(specie), atomic_weight(specie), specie + "." + settings.pspname + ".upf"])
            else:
                awp.append([str(specie), atomic_weight(specie), specie + "." + settings.pspname + ".upf"])

            if not os.path.exists(system.runpspdir + '/' + specie + '.' + settings.pspname + '.upf'):
                os.symlink(
                    settings.pspdir + '/' + settings.pspname + '/' + specie + '/' + specie + '.' + settings.pspname + '.upf',
                    system.runpspdir + '/' + specie + '.' + settings.pspname + '.upf')
            if settings.psp_dict[specie] > max_ecut:
                max_ecut = settings.psp_dict[specie]

        atom_coords = []
        for i, atom in enumerate(structure.species):
            atom_coords.append([atom, structure.coords[i]])

        self.kwargs = dict(self.kwargs, ntyp=len(unique_atoms), nat=natoms,
                                    lattice_params=structure.lattice, atoms_weights_psps=awp,
                                    cart_coords=atom_coords,
                                    pseudo_dir=system.runpspdir,
                                    num_kpoints=len(structure.kgrid),
                                    kgrid=copy.deepcopy(structure.kgrid),
                                    ecutwfc=copy.copy(max_ecut))

    def transition_metal(self):
        t_metals = set(['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', \
                        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', \
                        'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', \
                        'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn'])
        atoms = self.system.structure.species

        unique_atoms = set(self.system.structure.species)
        if list(t_metals & unique_atoms):
            t_metal_atom = list(t_metals & unique_atoms)[0]
            t_metal_atom_index = unique_atoms.index(t_metal_atom)
            return t_metal_atom_index
        else:
            return 1

    def write_pwscf_inputs(self):

        if not os.path.exists(self.rundir):
            os.mkdir(self.rundir)

        with open(self.rundir + '/input.qe', mode='w') as f:

            # Write Control Block
            f.write(self.input_control["header"] + "\n")
            for key, value in self.input_control.iteritems():

                if key is 'header':
                    pass
                else:
                    if isinstance(value, int) or isinstance(value, float):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    elif isinstance(value, bool):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    else:
                        f.write("{0:30} = '{1}' \n".format(str(key), str(value)))
            f.write("/" + "\n")

            # Write System Block
            f.write(self.input_system["header"] + "\n")
            for key, value in self.input_system.iteritems():

                if key is 'header':
                    pass
                elif key is 'starting_magnetization':
                    key='starting_magnetization({0})'.format(self.transition_metal())
                    f.write("{0:30} = {1} \n".format(key, str(value)))
                else:
                    if isinstance(value, int) or isinstance(value, float):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    elif isinstance(value, bool):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    else:
                        f.write("{0:30} = '{1}' \n".format(str(key), str(value)))
            f.write("/" + "\n")

            # Write Electrons Block
            f.write(self.input_electrons["header"] + "\n")
            for key, value in self.input_electrons.iteritems():

                if key is 'header':
                    pass
                else:
                    if isinstance(value, int) or isinstance(value, float):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    else:
                        f.write("{0:30} = '{1}' \n".format(str(key), str(value)))
            f.write("/" + "\n")

            # Write CELL_PARAMETERS Block
            f.write(self.input_cell_params["header"] + "\n")
            for key, value in self.input_cell_params.iteritems():

                if key is 'header':
                    pass
                elif key is 'lattice_params':
                    for row in value:
                        f.write(
                            '{0:15}  {1:15}  {2:15}'.format(str(format(row[0], '.10f')), str(format(row[1], '.10f')),
                                                            str(format(row[2], '.10f'))) + '\n')

                else:
                    if isinstance(value, int) or isinstance(value, float):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    elif isinstance(value, bool):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    else:
                        f.write("{0:30} = '{1}' \n".format(str(key), str(value)))
            f.write("\n")

            # Write ATOMIC SPECIES BLOCK
            f.write(self.input_atomic_species["header"] + "\n")
            for key, value in self.input_atomic_species.iteritems():

                if key is 'header':
                    pass
                elif key is 'atoms_weights_psps':
                    for row in value:
                        f.write('{0:5}  {1:15}  {2:15}'.format(row[0], str(format(row[1], '.10f')), row[2]) + '\n')
                else:
                    if isinstance(value, int) or isinstance(value, float):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    elif isinstance(value, bool):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    else:
                        f.write("{0:30} = '{1}' \n".format(str(key), str(value)))
            f.write("\n")

            # Write Atomic Positiostarting_magnetizationns Block
            f.write(self.input_atomic_positions["header"] + "\n")
            for key, value in self.input_atomic_positions.iteritems():

                if key is 'header':
                    pass
                elif key is 'cart_coords':
                    for row in value:
                        f.write('{0:5}'.format(row[0]))
                        f.write('{0:15}  {1:15}  {2:15}'.format(str(format(row[1][0], '.10f')),
                                                                str(format(row[1][1], '.10f')),
                                                                str(format(row[1][2], '.10f'))) + '\n')
                else:
                    if isinstance(value, int) or isinstance(value, float):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    elif isinstance(value, bool):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    else:
                        f.write("{0:30} = '{1}' \n".format(str(key), str(value)))
            f.write("\n")

            # Write K POINTS Block
            f.write(self.input_k_points["header"] + "\n")
            f.write(str(self.input_k_points["num_kpoints"]) + "\n")

            for key, value in self.input_k_points.iteritems():

                if key is 'header':
                    pass
                elif key is 'num_kpoints':
                    pass
                elif key is 'kgrid':
                    a = len(value)

                    f.write('\n'.join(' '.join(format(cell, '.10f') for cell in row) for row in value))

                    # for row in value:
                    # f.write('{0:15} {1:15} {2:15} {3:15}'.format(str(format(row[0], '.10f')), str(format(row[1], '.10f')), str(format(row[2], '.10f')), str(format(kgrid_weights, '.10f'))) + '\n')
                else:
                    if isinstance(value, int) or isinstance(value, float):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    elif isinstance(value, bool):
                        f.write("{0:30} = {1} \n".format(str(key), str(value)))
                    else:
                        f.write("{0:30} = '{1}' \n".format(str(key), str(value)))

            f.write("\n")

    def control_pwscf(self):

        if os.path.exists(self.outputfile):
            self.running=True
            if os.path.exists(self.rundir + '/pwscf.save/data-file.xml'):
                with open(self.outputfile) as f:
                    for line in f:
                        if line.strip() == 'iteration #  {0}'.format(self.input_electrons['electron_maxstep']):
                            self.complete = False
                            break
                        if line.strip() == 'JOB DONE.':
                            self.complete = True
                            self.running = False

    def execute(self):
        cur_dir=os.getcwd()
        os.chdir(self.rundir)
        script_name='job.'+self.job.name+'.sh'
        with open(script_name, 'w') as f:
            f.write(self.job.script)


        os.system(self.job.run_command + " " + script_name)
        os.chdir(cur_dir)

def generate_pwscf(**kwargs):
    pwscf = Pwscf(**kwargs)
    return pwscf

