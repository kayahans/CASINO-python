#! /usr/bin/env python
from casino_python import settings,System
import pwscf
from casino import generate_casino
from analyze_results import results
from pw2casino import  pw2casino
from pprint import pprint

settings(rootdir="./",
         pspdir="./Psps",
         pspname="OPT",
         pspcutoffs="./psp.txt",
         machinename="cruller",
         )


#print pprint(vars(settings))

sims = []
for scell in [1]:
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
    scf = pwscf.generate_pwscf(
    system=generic,
    input_dft='lda'
    )
    sims.append(scf)
    psi = pw2casino(
        dft=scf
    )
    sims.append(psi)
    vmc = generate_casino(
        dft=scf,
        psi=psi,
        qmc_prev=None,
    )
    sims.append(vmc)
    dmc = generate_casino(
        dft=scf,
        psi=psi,
        qmc_prev=vmc
    )
#for sim in sims:
#    print pprint(vars(sim))






#dmc = generate_casino(
#    psi=psi,
#    qmc_prev=vmc
#)