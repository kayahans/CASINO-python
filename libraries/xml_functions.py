#! /usr/bin/env python

import xml.etree.cElementTree as ET
import numpy as np
import fnmatch
import os,sys
import json
from structure import Structure





def get_dft_correction(qedir="."):
    print np.version.version
    qedir = qedir + "/pwscf.save/"

    os.chdir(qedir)

    folder_list = []
    datafile = "data-file.xml"
    datafile_tree = ET.ElementTree(file=datafile)
    eig = datafile_tree.find("EIGENVALUES")

    weights = dict()
    spin_up_energies = dict()
    spin_down_energies = dict()

    for kpoint in eig.getchildren():
        up = kpoint.find("DATAFILE.1")
        up_file = up.attrib["iotk_link"]
        print up_file
        up_energy = get_energies(up_file)

        spin_up_energies.update({int(up_file[4:8]): up_energy})

        down = kpoint.find("DATAFILE.2")

        down_file = down.attrib["iotk_link"]
        down_energy = get_energies(down_file)

        spin_down_energies.update({int(down_file[4:8]): down_energy})

        weight = kpoint.findtext("WEIGHT")

        if int(up_file[4:8]) == int(down_file[4:8]):
            weights.update({int(up_file[4:8]):float(weight)})


    # Normalize weights

    sum_weight = 0
    for kpoint, weight in weights.iteritems():
        sum_weight += weight

    for kpoint, weight in weights.iteritems():
        weights[kpoint] = weight / sum_weight

    total_energy = 0

    for kpoint, weight in weights.iteritems():
        total_energy += weights[kpoint]*(spin_up_energies[kpoint]+spin_down_energies[kpoint])


    #print json.dumps(spin_up_energies, indent=2)
    #print json.dumps(spin_down_energies, indent=2)
    #print json.dumps(weights, indent=2)
    print total_energy



def get_energies(file=None):
    tree = ET.ElementTree(file=file)
    eig = tree.findtext("EIGENVALUES")
    list = []
    for line in eig.splitlines():
        assert isinstance(line, str)
        line=line.strip()
        if line != '':
            list.append(float(line))
    eig_arr = np.array(list)

    occ = tree.findtext("OCCUPATIONS")
    list = []
    for line in occ.splitlines():
        assert isinstance(line, str)
        line = line.strip()
        if line != '':
            list.append(float(line))
    occ_arr = np.array(list)
    k_point_energy = np.dot(eig_arr, occ_arr)
    return k_point_energy

if len(sys.argv) > 1:
    get_dft_correction(qedir=sys.argv[1])
else:
    get_dft_correction(qedir="/Users/kayahan/Desktop/runs/Ni/OPT/icsd/radius.1/qe")