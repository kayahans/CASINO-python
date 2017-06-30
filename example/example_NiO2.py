#! /usr/bin/env python
from casino_python import settings,System
from pwscf import generate_pwscf
from analyze_results import results
from vasp import  Vasp
import os
import subprocess
import re
import numpy as np
import structure
k=Vasp()
print k
settings(rootdir="./",
         pspdir="./Psps",
         pspname="OPT",
         pspcutoffs="./psp.txt",
         machinename="cruller",
         )

generic = System(
        name="NiO2",
        structdir="./Structs",
        source="icsd",
        supercells=[1],
        folded=True,
        real_or_complex='Complex'
        #real_or_complex='Real'
)



sims = []
scf = generate_pwscf(
        system=generic,
        input_dft='lda'
)