#! /usr/bin/env python
from casino_python import settings,System
from pwscf import generate_pwscf
from casino import generate_casino
from analyze_results import results
from pw2casino import  pw2casino

settings(rootdir="./",
         pspdir="./Psps",
         pspname="OPT",
         pspcutoffs="./psp.txt",
         machinename="cruller",
         )


sims = []
for scell in [1,2,3,4]:
    generic = System(
        name="NiO2",
        structdir="./Structs",
        source="icsd",
        supercell_size=[scell],
        folded=True,
        real_or_complex='Complex',
        mindistance=16
        #real_or_complex='Real'
    )

    scf = generate_pwscf(
        system=generic,
        input_dft='lda'
    )

    psi = pw2casino(
    dft=scf
    )

#vmc = generate_casino(
#    psi=psi,
#    qmc_prev=None,
#)

#dmc = generate_casino(
#    psi=psi,
#    qmc_prev=vmc
#)