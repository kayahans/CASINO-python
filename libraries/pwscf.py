import sys,os,copy
from collections import OrderedDict
from units import atomic_weight
import fnmatch
from casino_python import settings
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
    header='&CELL cubic',
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

    def __init__(self, system=None, **kwargs):

        if system is None:
            print "system keyword must be used in generate_pwscf block"
            sys.exit()

        self.system = system
        self.kwargs = dict()

        for scell_num, str in system.scells.iteritems():
            self.kwargs.update({scell_num :kwargs})

        self.process_system_inputs()
        self.process_pwscf_inputs()
        del self.kwargs
        self.write_pwscf_inputs()



    #self.k_point_lattice = gen_k_points_by_density(k_point_density)

    def process_pwscf_inputs(self):

        self.input_control = dict()
        self.input_system= dict()
        self.input_electrons= dict()
        self.input_cell_params= dict()
        self.input_atomic_species= dict()
        self.input_atomic_positions= dict()
        self.input_k_points= dict()


        for scell_num, scell_kwargs in self.kwargs.iteritems():
            self.input_control[scell_num] = default_input_control.copy()
            self.input_system[scell_num] = default_input_system.copy()
            self.input_electrons[scell_num] = default_input_electrons.copy()
            self.input_cell_params[scell_num] = default_cell_params.copy()
            self.input_atomic_species[scell_num] = default_atomic_species.copy()
            self.input_atomic_positions[scell_num] = default_atomic_positions.copy()
            self.input_k_points[scell_num] = default_k_points.copy()

            for key, value in scell_kwargs.iteritems():
                if key in self.input_control[scell_num]:
                    self.input_control[scell_num][key] = value
                if key in self.input_system[scell_num]:
                    self.input_system[scell_num][key] = value
                if key in self.input_electrons[scell_num]:
                    self.input_electrons[scell_num][key] = value
                if key in self.input_cell_params[scell_num]:
                    self.input_cell_params[scell_num][key] = value
                if key in self.input_atomic_species[scell_num]:
                    self.input_atomic_species[scell_num][key] = value
                if key in self.input_atomic_positions[scell_num]:
                    self.input_atomic_positions[scell_num][key] = value
                if key in self.input_k_points[scell_num]:
                    self.input_k_points[scell_num][key] = value

    def process_system_inputs(self):

        system = self.system
        for str_num, structure in system.structures.iteritems():

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
                    max_ecut=settings.psp_dict[specie]

            atom_coords = []
            for i, atom in enumerate(structure.species):
                atom_coords.append([atom, structure.coords[i]])

            self.kwargs[str_num] = dict(self.kwargs[str(str_num)], ntyp = len(unique_atoms), nat = natoms, lattice_params = structure.lattice, atoms_weights_psps = awp,
                                        cart_coords = atom_coords,
                                        pseudo_dir = system.runpspdir,
                                        num_kpoints = len(structure.kgrid),
                                        kgrid = copy.deepcopy(structure.kgrid),
                                        ecutwfc=copy.copy(max_ecut))

    def transition_metal(self):
        t_metals = set(['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', \
                        'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', \
                        'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', \
                        'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn'])
        atoms = self.system.structures.species

        unique_atoms = set(self.system.structures.species)
        if list(t_metals & unique_atoms):
            t_metal_atom = list(t_metals & unique_atoms)[0]
            t_metal_atom_index = atoms.index(t_metal_atom)
            self.kwargs["starting_magnetization"]



    def write_pwscf_inputs(self):

        rundirs = self.system.rundirs
        assert isinstance(rundirs, dict), "run directories must be of type dict"
        for runnum, rundir in rundirs.iteritems():
            rundir +=  '/qe'
            if not os.path.exists(rundir):
                os.makedirs(rundir)

            with open(rundir+'/input.qe', mode='w') as f:

                #Write Control Block
                f.write( self.input_control[runnum]["header"] + "\n" )
                for key, value in self.input_control[runnum].iteritems():

                    if key is 'header':
                        pass
                    else:
                        if isinstance(value, int) or isinstance(value,float):
                            print isinstance(value)
                            f.write( "{0:30} = {1} \n".format(str(key), str(value)) )
                        else:
                            f.write( "{0:30} = '{1}' \n".format(str(key), str(value)) )
                f.write( "/" + "\n" )

                #Write System Block
                f.write( self.input_system[runnum]["header"] + "\n" )
                for key, value in self.input_system[runnum].iteritems():

                    if key is 'header':
                        pass
                    else:
                        if isinstance(value, int) or isinstance(value,float):
                            f.write( "{0:30} = {1} \n".format(str(key), str(value)) )
                        else:
                            f.write( "{0:30} = '{1}' \n".format(str(key), str(value)) )
                f.write( "/" + "\n" )

                #Write Electrons Block
                f.write( self.input_electrons[runnum]["header"] + "\n" )
                for key, value in self.input_electrons[runnum].iteritems():

                    if key is 'header':
                        pass
                    else:
                        if isinstance(value, int) or isinstance(value, float):
                            f.write( "{0:30} = {1} \n".format(str(key), str(value)) )
                        else:
                            f.write( "{0:30} = '{1}' \n".format(str(key), str(value)) )
                f.write( "/" + "\n" )

                #Write Cell Block
                f.write( self.input_cell_params[runnum]["header"] + "\n" )
                for key, value in self.input_cell_params[runnum].iteritems():

                    if key is 'header':
                        pass
                    elif key is 'lattice_params':
                        for row in value:
                            f.write('{0:15}  {1:15}  {2:15}'.format(str(format(row[0], '.10f')), str(format(row[1], '.10f')), str(format(row[2], '.10f'))) + '\n')

                    else:
                        if isinstance(value, int) or isinstance(value, float):
                            f.write( "{0:30} = {1} \n".format(str(key), str(value)) )
                        else:
                            f.write( "{0:30} = '{1}' \n".format(str(key), str(value)) )
                f.write( "/" + "\n" )

                #Write ATOMIC SPECIES BLOCK
                f.write( self.input_atomic_species[runnum]["header"] + "\n" )
                for key, value in self.input_atomic_species[runnum].iteritems():

                    if key is 'header':
                        pass
                    elif key is 'atoms_weights_psps':
                        for row in value:
                            f.write('{0:5}  {1:15}  {2:15}'.format(row[0], str(format(row[1], '.10f')), row[2]) + '\n')
                    else:
                        if isinstance(value, int) or isinstance(value, float):
                            f.write( "{0:30} = {1} \n".format(str(key), str(value)) )
                        else:
                            f.write( "{0:30} = '{1}' \n".format(str(key), str(value)) )
                f.write( "/" + "\n" )

                #Write Atomic Positions Block
                f.write( self.input_atomic_positions[runnum]["header"] + "\n" )
                for key, value in self.input_atomic_positions[runnum].iteritems():

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
                            f.write( "{0:30} = {1} \n".format(str(key), str(value)) )
                        else:
                            f.write( "{0:30} = '{1}' \n".format(str(key), str(value)) )
                f.write( "/" + "\n" )

                #Write K POINTS Block
                f.write( self.input_k_points[runnum]["header"] + "\n" )
                f.write( str(self.input_k_points[runnum]["num_kpoints"]) + "\n" )
                for key, value in self.input_k_points[runnum].iteritems():

                    if key is 'header':
                        pass
                    elif key is 'num_kpoints':
                        pass
                    elif key is 'kgrid':
                        a=len(value)
                        kgrid_weights = 1.0 / len(value)

                        for row in value:
                            f.write('{0:15} {1:15} {2:15} {3:15}'.format(str(format(row[0], '.10f')), str(format(row[1], '.10f')), str(format(row[2], '.10f')), str(format(kgrid_weights, '.10f'))) + '\n')
                    else:
                        if isinstance(value, int) or isinstance(value, float):
                            f.write( "{0:30} = {1} \n".format(str(key), str(value)) )
                        else:
                            f.write( "{0:30} = '{1}' \n".format(str(key), str(value)) )



def generate_pwscf(**kwargs):
    pwscf = Pwscf(**kwargs)
    return pwscf