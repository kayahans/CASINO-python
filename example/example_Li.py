#! /usr/bin/env python
from casino_python import settings,System
from pwscf import generate_pwscf
from casino import generate_casino
from analyze_results import results
from pw2casino import  pw2casino
from machines import Job

settings(rootdir="./",
         pspdir="./Psps",
         pspname="OPT",
         pspcutoffs="./psp.txt",
         machinename="cruller"
         )


for scell in [1,2,4,8]:
    generic = System(
        name="Li",
        structdir="./Structs",
        source="icsd",
        supercell_size=scell,
        folded=True,
        real_or_complex='Real',
        mindistance=16
        #real_or_complex='Real'
    )
    scf = generate_pwscf(
        system = generic,
        input_dft = 'lda',
        job=Job(nodes=2,time=12,name='dft',app='pwscf.x -pw2casino')

    )
    psi = pw2casino(
        dft=scf
    )
    vmc = generate_casino(
        dft=scf,
        psi=psi,
        qmc_prev=None,
        job=Job(nodes=8, time=12, name='vmc', app='casino')
    )
    dmc = generate_casino(
        dft=scf,
        psi=psi,
        qmc_prev=vmc,
        job=Job(nodes=80, time=12, name='vmc', app='casino')
    )

