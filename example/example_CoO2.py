#! /usr/bin/env python
from casino_python import settings,System
from pwscf import generate_pwscf
from analyze_results import results

settings(rootdir="./",
         pspdir="./Psps",
         pspname="OPT",
         pspcutoffs="psp.txt",
         #pspcutoffs=280
         machinename="cruller",
             )

generic = System(
        name="CoO2",
        structdir="./Structs",
        source="icsd",
        supercells=[1,2,3,4],
        folded=True,
        real_or_complex='Complex',
        mindistance=16
)

sims = []
scf = generate_pwscf(
        system=generic,
        input_dft='lda'
)

#print generic.structures["radius.1"].kgrid