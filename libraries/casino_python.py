# from machines import Machine
import os
import socket
import glob
import fnmatch
import copy
from vasp import Vasp
from error_handler import error,warning

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
            setattr(self,key,value)

        self.rootdir = os.path.abspath(self.rootdir)
        self.pspdir = os.path.abspath(self.pspdir)
        if not os.path.exists(self.rootdir):
            error("rootdir does not exist")

        if not os.path.exists(self.pspdir):
            error("pspdir:(" +self.pspdir+ ") does not exist")

        if not isinstance(self.pspcutoffs, int):
            if not os.path.exists(self.pspcutoffs):
                self.pspcutoffs = self.pspdir + '/' + self.pspname + '/' + self.pspcutoffs
                if not os.path.exists(self.pspcutoffs):
                    error("pspcutoffs file does not exist")
        else:
            psp_dict = []
            with open(self.pspcutoffs, 'r') as inf:
                for line in inf:
                    psp_dict.append(eval(line))
            self.psp_dict = psp_dict


        if not socket.gethostname() == self.machinename:
            message = "System hostname (" + str(
                socket.gethostname()) + ") is different from input hostname (" + self.machinename + "), but running anyway!"
            warning(message)


settings = Settings()

class System():
    """ Input file locations """
    def __init__(self, name=None, source=None, supercells=None, structdir=None, folded = True):

        assert isinstance(name, str), Settings.error("System settings are wrong")
        assert isinstance(source, str), Settings.error("System settings are wrong")

        self.structdir = os.path.abspath(structdir + "/" + name + "/" + source + "/")

        vaspInput = False
        if self.structdir is not None:
            vaspInput = True
            # Reads only vasp input for now
            for file in os.listdir(self.structdir):
                if fnmatch.fnmatch(file, '*.vasp'):
                    self.inputfile = os.path.abspath(self.structdir + "/" + file)

            if vaspInput is True:
                unit_cell = Vasp.read_poscar(self.inputfile)
                unit_cell = unit_cell.structure
                self.scells = unit_cell.get_supercells(supercells)

                self.structures = dict()

                if not folded:
                    for key,value in self.scells.iteritems():
                        self.structures.update({key : unit_cell.make_scell(value)})
                else:
                    for key, value in self.scells.iteritems():
                        unit_cell.kgrid = copy.deepcopy(unit_cell.gen_kpts_lattice(value))
                        self.structures.update({key : copy.deepcopy(unit_cell)})


        for file in glob.glob(self.inputfile):
            self.inputfile = file

        self.runpspdir = os.path.abspath(
            settings.rootdir + '/runs/' + name + '/' + settings.pspname + '/' + "psps")

        if not os.path.exists(self.runpspdir):
            os.makedirs(self.runpspdir)

        #for atoms in self.unique_atoms:
        #    assert isinstance(atoms, str), Settings.error("Atoms in POSCAR are not strings")
        #    os.symlink((settings.pspdir + '/' + atoms.title() + '.upf'),
        #               (self.matpspdir + '/' + atoms.title() + settings.pspname + '.upf'))
        #    os.symlink((settings.pspdir + '/' + atoms.lower() + '_pp.data'),
        #               (self.matpspdir + '/' + atoms.title() + settings.pspname + '_data'))

        rundirs = dict()
        for key, value in self.scells.iteritems():
            temp_rundir = settings.rootdir + '/runs/' + name + '/' + settings.pspname + '/' + source + '/' + str(key)
            rundirs.update({key : temp_rundir})
            if not os.path.exists(temp_rundir):
                os.makedirs(temp_rundir)

        self.rundirs = rundirs
