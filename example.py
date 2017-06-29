#! /usr/bin/env python
from casino_python import settings,System
from pwscf import generate_pwscf
from analyze_results import results

settings(rootdir="/Users/kayahan/Desktop",
         pspdir="/Users/kayahan/Desktop",
         pspname="OPT",
         pspcutoffs="/Users/kayahan/Desktop/psp.txt",
         machinename="cruller",
         )

generic = System(
        name="CoO",
        structdir="/Users/kayahan/Desktop/Structs",
        source="icsd",
        supercells=[1,2,3,4],
        folded=True,
    )

sims = []
scf = generate_pwscf(
        system=generic,
        input_dft='lda'

    )

