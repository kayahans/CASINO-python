import structure
import numpy as np
import os,sys
import copy
from error_handler import error, warning

class Vasp:

    def __init__(self, comment=None, scale=None, struct=None, iscartesian=True):

        if comment is None:
            self.comment = "Default VASP inp_str"
        else:
            self.comment = comment

        if scale is None:
            self.scale = 1.0
        else:
            self.scale = scale

        assert isinstance(struct, structure.Structure), error("structure type is incorrect")

        self.structure = struct
        self.iscartesian = iscartesian




    @staticmethod
    def read_poscar(filename):
        """
Reads VASP POSCAR files, and returns VASP object, currently only cartesian coordinates are supported
Input VASP file must be in the following format (Cartesian POSCAR output using VESTA, http://jp-minerals.org/vesta/en/, will suffice).
Whitespaces are not important:
Ni1 O2
1.0
        4.8754000664         0.0000000000         0.0000000000
        0.0000000000         2.8141000271         0.0000000000
       -3.2680774749         0.0000000000         4.5253056414
   Ni    O
    2    4
Cartesian
     0.000000000         0.000000000         0.000000000
     2.437700033         1.407050014         0.000000000
     1.648646319         0.000000000         1.045345631
    -0.041323873         0.000000000         3.479960010
     4.086346498         1.407050014         1.045345631
    -2.479023761         1.407050014         3.479960010
"""
        filename=os.path.abspath(filename)
        with open(filename) as f:
            lines = f.read().splitlines()

        comment = lines[0]
        scale   = lines[1]
        lattice = [[float(x) for x in lines[2].split()], [float(x) for x in lines[3].split()], [float(x) for x in lines[4].split()]]

        if not isinstance(lines[5], str):
            """If the fifth line is string, it means that there are no species information"""
            error("No species information in the POSCAR")

        species_strs = lines[5].split()
        species_counts=[int(x) for x in lines[6].split()]

        if len(species_strs) != len(species_counts):
            error("Number of species and number of species strings do not match")

        species = []

        for i in range(0, len(species_counts)):
            species += [species_strs[i]] * species_counts[i]

        iscartesian = []
        if lines[7].startswith("C") or lines[7].startswith("c"):
            iscartesian = True
        elif lines[7].startswith("R") or lines[7].startswith("D") or lines[7].startswith("r") or lines[7].startswith("d"):
            iscartesian = False
        else:
            error("Coordinates are not Cartesian or Reciprocal!")

        coords = []
        if iscartesian:
            for i in range(8,len(lines)):

                coords += [[float(x) for x in lines[i].split()[0:3]]]
        else:
            for i in range(8,len(lines)):
                coords += [[float(x) for x in lines[i].split()[0:3]]]

            #Dot product with the lattice, always print Cartesian
            crd=np.dot(np.array(coords),np.array(lattice))
            coords=crd.tolist()

        return Vasp(comment=comment, scale=scale, struct=structure.Structure(lattice=lattice, species=species, coords=coords), iscartesian=iscartesian)

    @staticmethod
    def write_poscar(inp_str, outfile='POSCAR'):
        """Writes POSCAR file from Structure or VASP object"""

        out_str = []
        assert isinstance(inp_str, structure.Structure) or isinstance(inp_str,Vasp), "'write_poscar' function input must be of type 'Structure' or 'Vasp'"

        # Look for VASP object first, since it is superclass
        if isinstance(inp_str, Vasp):
            out_str = copy.deepcopy(inp_str)
        elif isinstance(inp_str,structure.Structure):
            out_str = Vasp(comment='Default VASP POSCAR', scale='1.0', struct=structure.Structure(lattice=inp_str.lattice, coords=inp_str.coords, species=inp_str.species), iscartesian=True)


        dir = os.getcwd()
        filename = dir + '/' + outfile
        with open(filename, 'w') as f:
            f.write(out_str.comment + '\n')
            f.write(str(out_str.scale)+ '\n')
            for row in out_str.structure.lattice:
                f.write('{0:15}  {1:15}  {2:15}'.format(str(format(row[0], '.10f')), str(format(row[1], '.10f')),
                                                        str(format(row[2], '.10f'))) + '\n')

            unique_elements = set(out_str.structure.species)

            f.write('\t'+ '\t'.join(unique_elements)+ '\n')

            el_num = [0] * len(unique_elements)

            for i, element in enumerate(unique_elements):
                el_num[i]=out_str.structure.species.count(element)

            f.write('\t'+ '\t'.join(str(x) for x in el_num) + '\n')
            if out_str.iscartesian:
                f.write("Cartesian" + '\n')
            for row in out_str.structure.coords:
                 f.write('{0:>12}  {1:>12}  {2:>12}'.format(str(row[0]), str(row[1]), str(row[2])) + '\n')
            f.close()

    def get_structure(self):
        """Return structure object from Structure Class"""
        return self.structure