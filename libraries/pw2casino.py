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

        self.control_pw2casino()
        self.xml = Pwxml(self.dft.rundir)

        if not self.complete:
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
        first = True
        with open(self.bwfn) as f:
            old = f.readlines()
            for i in range(0, len(old)):
                line = old[i]
                if old[i].split() == kpoint_s.split():
                    if first:
                        first = False
                        header = old[0:i]

                    k_info = old[i + 1].split()
                    k_list.append(kpoints(num=k_info[0], nbnds_up=k_info[1], nbnds_down=k_info[2], coords=k_info[3:],
                                          start=i + 2))
                    if len(k_list) > 1:
                        k_list[-2].end = i
            k_list[-1].end = len(old) - 1

            header[-1] = '\t' + str(self.dft.system.scell_size) + '\n'
            if not os.path.exists(self.rundir):
                os.mkdir(self.rundir)
            with open(self.rundir+'/summary.txt', 'w') as f:
                f.write('# twist_wavefunction_name number_of_spin_up_electrons number_of_spin_down_electrons')
                for index, kpt in enumerate(k_list):
                    assert isinstance(kpt, kpoints)
                    self.twists.append(sys_dir + '/qe_wfns/bwfn.{0:0>3}.data'.format(index))
                    header[0] = 'bwfn.{0:0>3}.data'.format(index) + '\t' + self.dft.input_control["title"]
                    self.neu.update({sys_dir + '/qe_wfns/bwfn.{0:0>3}.data'.format(index): self.xml.up_nelect[index]})
                    self.ned.update({sys_dir + '/qe_wfns/bwfn.{0:0>3}.data'.format(index): self.xml.down_nelect[index]})
                    f.write('bwfn.{0:0>3}.data'.format(index) + ' ' + str(self.xml.up_nelect[index]) + ' ' + str(self.xml.down_nelect[index])+'\n')
                    with open(sys_dir + '/qe_wfns/bwfn.{0:0>3}.data'.format(index), 'w') as output:
                        for item in header:
                            output.write(item)
                        output.write(kpoint_s + '\n')
                        output.write('1' + '\t' + k_list[index].nbnds_up + '\t' + k_list[index].nbnds_down + '\t')
                        for item in kpt.coords:
                            output.write(item + '\t')
                        output.write('\n')
                        for item in old[kpt.start:kpt.end]:
                            output.write(item)

    def control_pw2casino(self):

        if self.dft.complete:
            if os.path.exists(self.rundir):
                if os.path.exists(self.dft.rundir+'/pwscf.bwfn.data'):
                    if len(os.listdir(self.rundir)) > 8:
                        self.complete = True

                if os.path.exists(self.rundir+'/summary.txt'):
                    self.complete = True








