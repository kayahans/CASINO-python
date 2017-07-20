#! /usr/bin/env python
from casino import generate_casino
from casino_python import settings,System
from machines import Job
from pw2casino import  pw2casino
from pwscf import generate_pwscf

settings(rootdir="./",
         pspdir="./Psps",
         pspname="OPT",
         pspcutoffs="./psp.txt",
         machinename="cruller"
         )

vmc_old = None
for scell in [1,2]:
    generic = System(
        name="Li",
        structdir="./Structs",
        source="icsd",
        supercell_size=scell,
        folded=True,
        real_or_complex='Complex',
        mindistance=16
        #real_or_complex='Real'
    )
    scf = generate_pwscf(
        system = generic,
        input_dft = 'lda',
        job=Job(nodes=2,time=12,name='dft',app='/home/apps/espresso-5.0.3/bin/pw.x -pw2casino')

    )
    psi = pw2casino(
        dft=scf
    )
    vmc = generate_casino(
        dft=scf,
        psi=psi,
        qmc='vmc_opt',
        qmc_prev=None,
        job=Job(nodes=8, time=12, name='vmc', app='casino')
    )
    vmc_old=vmc
    dmc = generate_casino(
        dft=scf,
        psi=psi,
        qmc='vmc_dmc',
        qmc_prev=vmc,
        job=Job(nodes=80, time=12, name='vmc', app='casino')
    )

