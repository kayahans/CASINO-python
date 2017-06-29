import xml.etree.cElementTree as ET
import numpy as np
from pwscf import Pwscf

def results(sim_list):


    pwscf_sims = object

    for sim in sim_list:
        if isinstance(sim, Pwscf):
            pwscf_sims.append(sim)


            print dir(sim)
            print dir(sim.system)
            print sim.system.rundirs

