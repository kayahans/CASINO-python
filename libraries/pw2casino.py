import pwscf
import os
from xml_functions import Pwxml

class kpoints:
    def __init__(self, num=None, nbnds_up=None, nbnds_down=None, neu=None, ned=None, coords=None,start=None, end=None,blips=None):
        self.num = num
        self.nbnds_up=nbnds_up
        self.nbnds_down=nbnds_down
        self.neu=neu,
        self.ned=ned,
        self.coords=coords
        self.start=start
        self.end=end
        self.blips=blips

class pw2casino:

    def __init__(self, dft=None):
        assert isinstance(dft, pwscf.Pwscf)

        #Inputs
        self.dft=dft

        #Class derived parameters
        self.rundir=self.dft.system.rundir + '/qe_wfns/'
        self.bwfn = str(dft.rundir) + '/pwscf.bwfn.data'
        self.twists=[]
        self.neu=dict()
        self.ned=dict()

        self.complete=False
        self.running=False
        self.dependencies=True

        self.control_pw2casino()

        if self.dependencies:
            print self.rundir + ' not ready!'
        else:
            if not self.complete and not self.running:
                print self.rundir + ' is running now!'
                self.xml = Pwxml(self.dft.rundir)
                print 'data-file.xml is read'
                self.divide_wfn()
                self.complete=True
            else:
                with open(self.rundir+'/summary.txt', 'r') as f:
                    for line in f:
                        if not line.startswith("#"):
                            self.neu.update({self.dft.system.rundir + '/qe_wfns/{0}'.format(line.split()[0]): line.split()[1]})
                            self.ned.update({self.dft.system.rundir + '/qe_wfns/{0}'.format(line.split()[0]): line.split()[2]})
                            self.twists.append(self.dft.system.rundir + '/qe_wfns/{0}'.format(line.split()[0]))

    def divide_wfn(self):
        sys_dir = self.dft.system.rundir
        structure = self.dft.system.structure
        k_grid = self.dft.input_k_points.values()

        kpoint_s = 'k-point # ; # of bands (up spin/down spin);            k-point coords (au)'
        k_list = []
        header = []
        new = True
        first = True

        print 'Reading pwscf.bwfn.data'

        first=True
        with open(self.bwfn) as f:
            for line in f:
                if first:
                    if line.split() == kpoint_s:
                        first=False
                    else:
                        header += line
                else:
                    if line.split() == kpoint_s.split():
                        new=True
                        k_info = line.split()
                        k_list.append(
                            kpoints(num=k_info[0], nbnds_up=k_info[1], nbnds_down=k_info[2], coords=k_info[3:]))

                    else:
                        new=False
        print header

    def control_pw2casino(self):

        if self.dft.complete:
            self.dependencies=False
            if os.path.exists(self.rundir):
                if os.path.exists(self.dft.rundir+'/pwscf.bwfn.data'):
                    if len(os.listdir(self.rundir)) > 8:
                        self.complete = True

                if os.path.exists(self.rundir+'/summary.txt'):
                    self.complete = True








