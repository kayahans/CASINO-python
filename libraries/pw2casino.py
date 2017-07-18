import pwscf
import os
from xml_functions import Pwxml
import sys

class kpoints:
    def __init__(self, num=None, nbnds_up=None, nbnds_down=None, neu=None, ned=None, coords=None,blips=None):
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
        kpoint_s=kpoint_s.split()

        k_list = []
        header = []
        after_kpts_str = False

        count=False
        numkpts=0
        index=0
        files=[]
        print 'Reading pwscf.bwfn.data'

        first=True

        if not os.path.exists(self.rundir):
            os.mkdir(self.rundir)

        multi_twist=False
        nscell=1
        after_kpts_str = False

        if self.dft.system.scell_size > 1:
            multi_twist=True
            nscell=self.dft.system.scell_size

        with open(self.bwfn) as f:
            with open(self.rundir + '/summary.txt', 'w') as g:
                g.write('# twist_wavefunction_name number_of_spin_up_electrons number_of_spin_down_electrons')
                for line in f:
                    line = line.split()

                    if first:
                        if line == 'Number of k-points'.split():
                            count = True
                            header += ' '.join(line) + '\n'
                        elif count == True:
                            numkpts = int(str(line[0])) / int(self.dft.system.scell_size)
                            files = [open(sys_dir + '/qe_wfns/bwfn.{0:0>3}.data'.format(x), 'w') for x in range(1, numkpts)]
                            count = False
                            header += ' '.join('\t' + str(self.dft.system.scell_size) + '\n')
                            first=False

                        else:
                            header += ' '.join(line)+ '\n'
                    else:
                        if multi_twist:
                            pass
                        else:
                            if line == kpoint_s:

                                prt = (index * 100) / numkpts
                                sys.stdout.write(str(prt) + ' percent complete' + '\r')
                                sys.stdout.flush()

                                index += 1
                                print index
                                for item in header:
                                    files[index - 1].write(item)

                                files[index - 1].write(' '.join(kpoint_s) + '\n')
                                # files[index - 1].write(str(self.dft.system.scell_size) + '\t' + k_list[index].nbnds_up + '\t' + k_list[index].nbnds_down + '\t' + '\n')

                                self.twists.append(sys_dir + '/qe_wfns/bwfn.{0:0>3}.data'.format(index))
                                header[0] = 'bwfn.{0:0>3}.data'.format(index) + '\t' + self.dft.input_control["title"]
                                self.neu.update(
                                    {sys_dir + '/qe_wfns/bwfn.{0:0>3}.data'.format(index): self.xml.up_nelect[index]})
                                self.ned.update(
                                    {sys_dir + '/qe_wfns/bwfn.{0:0>3}.data'.format(index): self.xml.down_nelect[index]})
                                g.write(
                                    'bwfn.{0:0>3}.data'.format(index) + ' ' + str(
                                        self.xml.up_nelect[index]) + ' ' + str(
                                        self.xml.down_nelect[index]) + '\n')
                                after_kpts_str = True

                            elif after_kpts_str:
                                print line
                                line[0]=str(1)
                                files[index - 1].write(' '.join(line) + '\n')
                                after_kpts_str = False

                            else:

                                files[index - 1].write(' '.join(line) + '\n')
        print ""

    def control_pw2casino(self):

        if self.dft.complete:
            self.dependencies=False
            if os.path.exists(self.rundir):
                if os.path.exists(self.dft.rundir+'/pwscf.bwfn.data'):
                    if len(os.listdir(self.rundir)) > 8:
                        self.complete = True

                if os.path.exists(self.rundir+'/summary.txt'):
                    self.complete = True








