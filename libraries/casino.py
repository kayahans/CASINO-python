import sys,os,copy
from error_handler import error, warning
from pwscf import Pwscf
from pw2casino import pw2casino
import fnmatch
import numpy as np
import subprocess
import re

# Input files are grouped into five
# default_qmc_system
# default_qmc_run
# default_vmc_opt
# default_dmc
# default_qmc_expectation_values

default_qmc_system = dict(
    header="#SYSTEM",
    neu=None,                           #*! Number of up electrons (Integer)
    ned=None,                           #*! Number of down electrons (Integer)
    periodic=True,                      #*! Periodic boundary conditions (Boolean)
    atom_basis_type="blip",             #*! Basis set type (Text)
    write_binary_blips=True,            #*! Write binary blips
    psi_s="slater",                     #*! Type of [anti]symmetrizing wfn (Text)
    complex_wf=False,                   #*! Wave function real or complex (Boolean)
    interaction="ewald",
)

default_qmc_run =dict(
    header="#RUN",
    newrun=True,                        #*! New run or continue old (Boolean)
    testrun=False,                      #*! Test run flag (Boolean)
    vm_reweight=False,                  #*! Whether to reweight during variance minimization.  F is more stable
    use_jastrow=True,                   #*! Use a Jastrow function (Boolean)
    backflow=False,                     #*! Use backflow corrections (Boolean)
    expot=False,                        #*! Use external potential (Boolean)
    timing_info=False,                  #*! Activate subroutine timers (Boolean)
    esupercell=False,                   #*! Energy/supercell in output (Boolean)
    neighprint=0,                       #*! Neighbour analysis (Integer)
    mpc_cutoff=30,                      #*! G vector cutoff for MPC (Physical)
    finite_size_corr=False,             #*! Eval. finite size correction (Boolean)
    forces=False,                       #*! Evaluate forces on atoms (Boolean)
)
default_vmc_opt = dict(
    header="#VMC",
    runtype="vmc_opt",                  #*! Type of calculation (Text)
    vmc_equil_step=15000,               #*! Number of equilibration steps (Integer)
    vmc_nstep=15000,                    #*! Number of steps (Integer)
    vmc_nblock=1,                       #*! Number of checkpoints (Integer)
    vmc_nconfig_write=15000,            #*! Number of configs to write (Integer)
    vmc_decorr_period=10,               #*! Set to a high value (e.g. 10) for configuration generation or 0 for automatic decorrelation
    opt_method="madmin",                #*! Mean absolute deviation of energies
    # opt_method       : varmin         #*! Variance of energies
    # opt_method       : varmin_linjas  #*! Variance of energies linear Jastrow parameters only
    opt_cycles=5,                       #*! Number of optimization cycles (Integer)
    opt_jastrow=True,                   #*! Optimize Jastrow factor (Boolean)
    opt_det_coeff=False,                #*! Optimize determinant coeffs (Boolean)
    opt_backflow=False,                 #*! Optimize backflow parameters (Boolean)
    opt_orbitals=False,                 #*! Optimize orbital parameters (Boolean)
)
default_dmc=dict(
    header="#DMC",
    runtype="vmc_dmc",                  #*! Type of calculation (Text)
    vmc_equil_step=15000,               #*! Number of equilibration steps (Integer)
    vmc_nstep=15000,                    #*! Number of steps (Integer)
    vmc_nblock=1,                       #*! Number of checkpoints (Integer)
    vmc_nconfig_write=15000,            #*! Number of configs to write (Integer)
    vmc_decorr_period=10,               #*! Set to a high value (e.g. 10) for configuration generation or 0 for automatic decorrelation
    dmc_equil_nstep= 500,               #*! Number of steps (Integer)
    dmc_equil_nblock= 1,                #*! Number of checkpoints (Integer)
    dmc_stats_nstep= 2500,              #*! Number of steps (Integer)
    dmc_stats_nblock= 5,                #*! Number of checkpoints (Integer)
    dmc_target_weight= 1152,            #*! Number of configs (DMC) (Integer)
    dmc_trip_weight= 3456,              #*! When the block is restarted.  Should be about 3x target weight
    dtdmc=0.01,                         #*! DMC time step (Double Precision)
    use_tmove=True,                     #*! Casula nl pp for DMC (Boolean)
)

default_qmc_expectation_values=dict(
    header="#Expectation values",
    density=False,                      #*! Accumulate density (Boolean)
    spin_density=False,                 #*! Accumulate spin densities (Boolean)
    pair_corr=False,                    #*! Accumulate rec. space PCF (Boolean)
    pair_corr_sph=False,                #*! Accumulate sph. real space PCF (Boolean)
    loc_tensor=False,                   #*! Accumulate localization tensor (Boolean)
    structure_factor=False,             #*! Accumulate structure factor (Boolean)
    struc_factor_sph=False,             #*! Accumulate sph. struc. factor (Boolean)
    onep_density_mat=False,             #*! Accumulate 1p density matrix (Boolean)
    twop_density_mat=False,             #*! Accumulate 2p density matrix (Boolean)
    cond_fraction=False,                #*! Accumulate cond fraction (Boolean)
    dipole_moment=False,                #*! Accumulate elec. dipole moment (Boolean)
    expval_cutoff=30.0,                 #*! G vector cutoff for expval (Physical)
    permit_den_symm = False,            #*! Symmetrize QMC charge data (Boolean)
    qmc_density_mpc = False,            #*! Use QMC density in MPC int (Boolean)
    int_sf = False                      #*! Calc ee int from strucfac (Boolean)
)

default_jastrow_parameters=dict(
    order_u=8,
    order_chi=8,
    order_f=3,
    three_body=False
)

class Casino:

    def __init__(self, dft=None, qmc_prev=None, qmc=None, job=None, psi=None, **kwargs):

        if dft is None:
            print "Previous DFT calculation must be provided"
            sys.exit()

        assert isinstance(dft, Pwscf), error("In generate_casino block, dft must be a Pwscf object")
        assert isinstance(psi, pw2casino)

        # Inputs
        self.dft = dft
        self.psi = psi
        self.qmc_prev = qmc_prev
        self.job = job
        self.rundir = []

        if qmc == 'vmc_opt' or qmc == 'vmc_dmc' or qmc == None:
            self.qmc = qmc
            if qmc == 'vmc_opt' or qmc == None:
                self.rundir=self.dft.system.rundir + '/vmc'
            else:
                self.rundir=self.dft.system.rundir + '/dmc'
        else:
            raise ValueError, 'qmc can be either "vmc_opt" or "vmc_dmc" '

        self.kwargs = kwargs

        # Class derived parameters
        self.system = dft.system
        self.kwargs.update({'scell_num': self.system.scell_size})
        self.qmc_inputs = []



        self.dependencies = False
        self.running = False
        self.complete = False

        self.control_casino()

        if not self.running and not self.complete and self.dependencies_complete:
            if self.dft.system.real_or_complex == 'Complex':
                self.kwargs.update({'complex_wf': True})

            if qmc_prev is None:
                # Then, no DMC calculation is performed until now, so perform VMC
                if not os.path.exists(self.rundir):
                    os.mkdir(self.rundir)

                self.kwargs.update({'opt_jastrow': True})
                self.correlation = self.write_correlation()

                # self.correlation=generate_default_correlation()
                self.qmc_defaults = [default_qmc_system, default_qmc_run, default_vmc_opt,
                                     default_qmc_expectation_values]
            else:
                assert isinstance(qmc_prev, Casino), error(
                    "In generate_casino block, qmc_prev must be a Casino object or None object")
                if qmc_prev.qmc == "vmc_opt":
                    if self.qmc == "vmc_opt":
                        # Copying the correlation data file from a smaller cell vmc_opt file

                        if not os.path.exists(self.rundir):
                            os.mkdir(self.rundir)

                        self.kwargs.update({'opt_jastrow': True})
                        self.correlation = self.transfer_correlation()

                        self.qmc_defaults = [default_qmc_system, default_qmc_run, default_vmc_opt,
                                             default_qmc_expectation_values]

                    elif self.qmc == "vmc_dmc":
                        # Previous calculation is VMC, then, this is a DMC calculation, jastrow already optimized

                        if not os.path.exists(self.rundir):
                            os.mkdir(self.rundir)

                        self.kwargs.update({'opt_jastrow': False})
                        self.type = "vmc_dmc"
                        self.correlation = self.best_correlation(qmc_prev)
                        self.qmc_defaults = [default_qmc_system, default_qmc_run, default_dmc,
                                             default_qmc_expectation_values]

            if self.qmc == 'vmc_opt':

                    for file in os.listdir(self.system.runpspdir):
                        if fnmatch.fnmatch(file, '*_pp.data'):
                            os.symlink(self.system.runpspdir + '/' + file, self.rundir)

                    ran_num = 2
                    filename = self.dft.system.rundir + '/qe_wfns/bwfn.{0:0>3}.data'.format(ran_num)
                    self.kwargs.update({'neu': self.psi.neu[filename]})
                    self.kwargs.update({'ned': self.psi.ned[filename]})
                    self.process_casino_inputs()
                    self.outputfile = self.rundir + '/out'

                    self.running = False
                    self.complete = False

                    if not self.complete and not self.running:
                        if not os.path.exists(self.rundir + 'bwfn.data'):
                            if not os.path.realpath(self.psi.twists[ran_num]) == os.path.realpath(
                                            self.rundir + '/bwfn.data'):
                                os.symlink(self.psi.twists[ran_num], self.rundir + '/bwfn.data')

                        self.write_casino_inputs(self.rundir)
                        self.execute()

            elif self.qmc == 'vmc_dmc':

                self.rundir_list = []

                for num, psi in enumerate(self.psi.twists):
                    if not os.path.exists(self.rundir + '/{0:0>3}'.format(num)):
                        os.mkdir(self.rundir + '/{0:0>3}'.format(num))
                        self.rundir_list.append(self.rundir + '/{0:0>3}'.format(num))
                        for file in os.listdir(self.system.runpspdir):
                            if fnmatch.fnmatch(file, '*_pp.data'):
                                os.symlink(self.system.runpspdir + '/' + file, self.rundir_list[-1])

                        if not os.path.realpath(self.psi.twists[num]) == os.path.realpath(
                                                self.rundir + '/{0:0>3}'.format(num) + '/bwfn.data'):
                            os.symlink(self.psi.twists[num], self.rundir + '/{0:0>3}'.format(num) + '/bwfn.data')

                        filename = self.dft.system.rundir + '/qe_wfns/bwfn.{0:0>3}.data'.format(num)
                        self.kwargs.update({'neu': self.psi.neu[filename]})
                        self.kwargs.update({'ned': self.psi.ned[filename]})
                        self.outputfile = self.rundir + '/{0:0>3}'.format(num) + '/out'

                        self.process_casino_inputs()
                        self.write_casino_inputs(self.rundir + '/{0:0>3}'.format(num))


        elif not self.dependencies_complete:
            print self.rundir + ' dependencies are not complete'
        elif self.running:
            print self.rundir + ' is still running'
        elif self.complete:
            print self.rundir + ' is running'



    def process_casino_inputs(self):

        qmc_inputs=copy.copy(self.qmc_defaults)

        for item in qmc_inputs:
            assert isinstance(item, dict)
            for key,value in self.kwargs.iteritems():
                if key in item.keys():
                    item[key]=value

        self.qmc_inputs=qmc_inputs

    def write_casino_inputs(self,rundir):


        '''Write casino input file'''
        with open(rundir+'/input', mode='w') as f:
            f.write("# -------------------#\n")
            f.write("# CASINO input file #\n")
            f.write("# -------------------#\n")
            f.write('\n')
            first=True
            for item in self.qmc_inputs:
                keylist=item.keys()
                keylist.sort()
                f.write(str(item['header'])+'\n')
                f.write('\n')
                for key in keylist:
                    if not key == 'header':
                        f.write(str(key)+' = '+str(item[key])+'\n')
                if first:
                    f.write('\n')
                    f.write("%block scell_matrix\n")
                    f.write(' '.join(map(str, self.system.scell))+'\n')
                    f.write("%block scell_matrix\n")
                    first=False
                f.write('\n')

    def write_correlation(self):

        params = default_jastrow_parameters

        structure = self.dft.system.structure
        unique_atoms = list(set(structure.species))

        correlation=[]

        # Correlation.data header

        correlation +='''START JASTROW
 Title
   Simple Jastrow
 Truncation order
   3

'''
        # U parameters in Jastrow

        correlation +='''START U TERM
 Number of sets
   1
 START SET 1
 Spherical harmonic l,m
   0 0
 Expansion order
   {ORDER}
 Spin-dep params (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud)
   1
 Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)
   0.d0                            1
 Parameter values  ;  Optimizable (0=NO; 1=YES)
 END SET 1
 END U TERM
 '''.format(ORDER=params["order_u"])

        # Chi parameters in Jastrow

        correlation_chi_header = '''
START CHI TERM
Number of sets ; labelling (1->atom in s. cell; 2->atom in p. cell; 3->species)
    {NSETS} 1
'''.format(NSETS=len(unique_atoms))

        correlation += correlation_chi_header

        for index, atom in enumerate(unique_atoms):
            atom_labels = ' '.join([str(i + 1) for i, val in enumerate(structure.species) if val == atom])
            natoms = (len(atom_labels) + 1) / 2

            correlation += """START SET {SET_NUM}
 Spherical harmonic l,m
    0 0
 Number of atoms in set
    {NUM_ATOMS}
 Labels of the atoms in this set
    {ATOM_LABELS}
 Impose electron-nucleus cusp (0=NO; 1=YES)
    0
 Expansion order
    {ORDER}
 Spin-dep params (0->u=d; 1->u/=d)
    0
 Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)
    0.d0                            1
 Parameter values  ;  Optimizable (0=NO; 1=YES)
END SET {SET_NUM}
""".format(SET_NUM=index + 1, NUM_ATOMS=natoms, ATOM_LABELS=atom_labels, ORDER=8)


        correlation += "END CHI TERM\n\n"

        # F parameters (3-body)

        if params["three_body"]:

            correlation_f_header = '''
START F TERM
Number of sets ; labelling (1->atom in s. cell; 2->atom in p. cell; 3->species)
    {NSETS} 1
'''.format(NSETS=len(unique_atoms))

            correlation +=correlation_f_header

            for index, atom in enumerate(unique_atoms):

                correlation +="""START SET {SET_NUM}
 Number of atoms in set
   {NUM_ATOMS}
 Labels of the atoms in this set
   {ATOM_LABELS}
 Prevent duplication of u term (0=NO; 1=YES)
   0
 Prevent duplication of chi term (0=NO; 1=YES)
   0
 Electron-nucleus expansion order
   {ORDER}
 Electron-electron expansion order
   {ORDER}
 Spin-dep params (0->uu=dd=ud; 1->uu=dd/=ud; 2->uu/=dd/=ud)
   1
 Cutoff (a.u.)     ;  Optimizable (0=NO; 1=YES)
   0.d0                            1
 Parameter values  ;  Optimizable (0=NO; 1=YES)
 END SET {SET_NUM}
""".format(SET_NUM=index+1, NUM_ATOMS=natoms, ATOM_LABELS=atom_labels, ORDER=3)


            correlation += " END F TERM \n\n"

        path = os.path.dirname(__file__)
        make_p_stars_path=os.path.realpath(path+"/../executables/make_p_stars")
        scell = self.system.scell

        if os.path.exists(make_p_stars_path):

            if len(scell) == 3:
                scell = np.array([[scell[0], 0, 0], [0, scell[1], 0], [0, 0, scell[2]]])
            elif len(scell) == 9:
                scell = np.array(scell)
            else:
                raise ValueError, 'Supercell matrix has wrong dimensions in make_p_stars utility'


            with open(self.rundir + '/star_input', 'w') as f:
                f.write(str(3)+'\n')
                f.write('O\n')
                np.savetxt(f,structure.lattice,fmt='%-10.6f')
                f.write(str(4)+'\n') #number of stars
                np.savetxt(f,scell,fmt='%-10.6f')

            command_str='{MAKE_P_STARS} < {STAR_INPUT} | grep -v Program | grep -v "*" | grep -v Please | grep -v "MAKE" | grep -v "==" |  egrep -v "^[[:space:]]*$|^#"'.format(MAKE_P_STARS=make_p_stars_path, STAR_INPUT=self.rundir+'/star_input')

            os.remove(self.rundir+'/star_input')

            process = subprocess.Popen(command_str, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            correlation_p, err = process.communicate()

            correlation +=correlation_p
        else:
            print "make_p_stars utility not found, correlation.data is created without P terms"


        correlation +="\nEND JASTROW \n"
        correlation=''.join(correlation)
        with open(self.rundir+'/correlation.data', 'w') as f:
            f.write(correlation)


    def transfer_correlation(self, prev_correlation):

        pass
        with open(self.rundir+'/correlation.data', 'w') as f:
            with open(prev_correlation, 'r') as g:
                modify_u_f=False
                modify_p=False
                for line in f.split(os.linesep):

                    # U parameters are transfered as is

                    if re.match('START CHI TERM', line) or re.match('START F TERM', line):
                        modify_u_f = True
                        modify_p = False
                    elif re.match('START P TERM'):
                        modify_u_f = False
                        modify_p = True
                    else:
                        if modify_u_f:
                            pass
                        elif modify_p:
                            pass
                        else:
                            f.write(''.join(line))

    def control_casino(self):

        if not self.psi.complete or not self.qmc_prev.complete:
            self.dependencies_complete=False
            self.running=False
            self.complete=False

        else:
            self.dependencies_complete=True
            if os.path.exists(self.outputfile):
                self.running = True
                with open(self.outputfile) as f:
                    for line in f:
                        if line.strip() == 'FINAL RESULT:':
                            self.complete = True
                            self.running = False
                            break
                        else:
                            self.complete = False


    def execute(self):
        cur_dir=os.getcwd()
        os.chdir(self.rundir)
        script_name='job.'+self.job.name+'.sh'
        with open(script_name, 'w') as f:
            f.write(self.job.script)

        #os.system(self.job.run_command + " " + script_name)
        os.chdir(cur_dir)

def generate_casino(**kwargs):
    qmc = Casino(**kwargs)
    return qmc