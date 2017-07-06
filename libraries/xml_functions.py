#! /usr/bin/env python

import xml.etree.cElementTree as ET
import numpy as np
import fnmatch
import os,sys
import json
from structure import Structure


class Pwxml:

    def __init__(self,qedir="."):
        self.qedir=qedir
        self.get_dft_correction(qedir=qedir)

    def get_dft_correction(self, qedir="."):

        qedir += "/pwscf.save/"

        datafile = qedir + "/data-file.xml"
        datafile_tree = ET.ElementTree(file=datafile)
        eig = datafile_tree.find("EIGENVALUES")

        weights = dict()
        self.spin_up_energies = dict()
        self.spin_down_energies = dict()
        self.up_nelect = []
        self.down_nelect = []
        self.total_energy=0

        for knum, kpoint in enumerate(eig.getchildren()):
            up = kpoint.find("DATAFILE.1")
            up_file = qedir+str(up.attrib["iotk_link"])
            (up_energy,neu) = get_energies(up_file)


            self.spin_up_energies.update({int(up_file[-19:-15]): up_energy})
            self.up_nelect.append(neu)

            down = kpoint.find("DATAFILE.2")
            down_file = qedir+str(down.attrib["iotk_link"])
            (down_energy,ned) = get_energies(down_file)
            print (neu,ned)
            self.spin_down_energies.update({int(down_file[-19:-15]): down_energy})
            self.down_nelect.append(ned)

            weight = kpoint.findtext("WEIGHT")

            if int(up_file[-19:-15]) == int(down_file[-19:-15]):
                weights.update({int(up_file[-19:-15]): float(weight)})

        # Normalize weights

        sum_weight = 0
        for kpoint, weight in weights.iteritems():
            sum_weight += weight

        for kpoint, weight in weights.iteritems():
            weights[kpoint] = weight / sum_weight


        for kpoint, weight in weights.iteritems():
            self.total_energy += weights[kpoint] * (self.spin_up_energies[kpoint] + self.spin_down_energies[kpoint])

        # print json.dumps(spin_up_energies, indent=2)
        # print json.dumps(spin_down_energies, indent=2)
        # print json.dumps(weights, indent=2)
        return self.total_energy


def get_energies(filename=None):
    tree = ET.ElementTree(file=filename)
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

    occ_arr = np.asarray(list)

    nelect=(occ_arr > 0.5).sum()
    print nelect
    k_point_energy = np.dot(eig_arr, occ_arr)
    return k_point_energy,nelect

#if len(sys.argv) > 1:
#    get_dft_correction(qedir=sys.argv[1])
#else:
#    get_dft_correction(qedir="/Users/kayahan/Desktop/runs/Ni/OPT/icsd/radius.1/qe")