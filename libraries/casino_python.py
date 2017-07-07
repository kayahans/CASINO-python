# from machines import Machine
import os
import socket
import glob
import fnmatch
import copy
from vasp import Vasp
from error_handler import error,warning
import structure

# Default settings
settings_list = dict(rootdir="./",
                    pspdir="./",
                    pspname="default",
                    pspcutoffs=280,
                    machinename="workstation",
                    results_only = 1)

class Settings(object):
    """Settings here"""

    def __call__(self, **kwargs):

        for key,value in kwargs.iteritems():
            if key in settings_list:
                settings_list[key]=value

        for key,value in settings_list.iteritems():
            setattr(self,key,value)

        self.rootdir = os.path.abspath(self.rootdir)
        if not os.path.exists(self.rootdir):
            error("rootdir does not exist")

        self.pspdir = os.path.abspath(self.pspdir)
        if not os.path.exists(self.pspdir):
            error("pspdir:(" +self.pspdir+ ") does not exist")

        self.pspname = str(self.pspname)
        self.machinename = self.machinename
        if not isinstance(self.pspcutoffs, int):
            if not os.path.exists(self.pspcutoffs):
                self.pspcutoffs = self.pspdir + '/' + self.pspname + '/' + self.pspcutoffs
                if not os.path.exists(self.pspcutoffs):
                    error("pspcutoffs file does not exist")

                psp_dict = dict()

                with open(self.pspcutoffs, 'r') as inf:
                    for line in inf:
                        (key, val) = line.split()
                        psp_dict[str(key)] = int(val)
                self.psp_dict = psp_dict
                inf.close()



        if not socket.gethostname() == self.machinename:
            message = "System hostname (" + str(
                socket.gethostname()) + ") is different from input hostname (" + self.machinename + "), but running anyway!"
            warning(message)



settings = Settings()

class System():
    """ Input file locations """
    def __init__(self, name=None, source=None, supercell_size=None, structdir=None, folded = True, real_or_complex='Real', mindistance=16):


        assert isinstance(name, str), error("System settings are wrong")
        assert isinstance(source, str), error("System settings are wrong")

        self.structdir = os.path.abspath(structdir + "/" + name + "/" + source + "/")
        self.scell_size=supercell_size
        self.inputfile=None
        self.real_or_complex = real_or_complex
        self.mindistance = mindistance

        vaspInput=False

        for file in os.listdir(self.structdir):
            if fnmatch.fnmatch(file, '*.vasp'):
                self.inputfile = os.path.abspath(self.structdir + "/" + file)
                vaspInput = True

        if self.inputfile is not None:
            # Reads only vasp input for now
            if vaspInput:

                unit_cell = Vasp.read_poscar(self.inputfile)
                unit_cell = unit_cell.structure

                #Transfer variables from System to structure

                assert isinstance(unit_cell, structure.Structure)
                self.scell = unit_cell.get_supercell(supercell_size)

                if real_or_complex == 'Real':
                    unit_cell.kgrid=unit_cell.gen_real_kpts_lattice(self.scell)
                if real_or_complex == 'Complex':
                    unit_cell.kgrid=unit_cell.gen_complex_kpts_lattice(self.scell,mindistance)

                self.unitcell=unit_cell
                self.structure = []

                if not folded:
                    self.structure=unit_cell.make_scell(self.scell)
                else:
                    self.structure=copy.deepcopy(unit_cell)

        for file in glob.glob(self.inputfile):
            self.inputfile = file

        self.runpspdir = os.path.abspath(
            settings.rootdir + '/runs/' + name + '/' + settings.pspname + '/' + "psps")

        if not os.path.exists(self.runpspdir):
            os.makedirs(self.runpspdir)

        self.rundir = settings.rootdir + '/runs/' + name + '/' + settings.pspname + '/' + source + '/radius.' + str(supercell_size)

        if not os.path.exists(self.rundir):
            os.makedirs(self.rundir)


