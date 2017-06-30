import numpy as np
from error_handler import error, warning
import vasp
import os,re
import subprocess
import copy

class Structure:

    def __init__(self, lattice=None, species=None, coords=None, real_or_complex='Real', mindistance=16):


        if lattice is None:
            self.lattice = []
            print "No lattice is set up"
        else:
            self.lattice = lattice

        if species is None:
            self.species = []
            print "No species in Structure"
        else:
            self.species = species

        if coords is None:
            self.coords = []
            print "No coords in Structure"
        else:
            self.coords = coords


        if real_or_complex == 'Complex':
            self.kgrid = self.gen_complex_kpts_lattice([1, 1, 1])
        elif real_or_complex == 'Real':
            self.kgrid = self.gen_real_kpts_lattice([1,1,1])
        else:
            print "No kgrid"

            self.mindistance=mindistance

        if len(species) != len(coords):
            print "Number of species is not equal to number of coordinates"
        elif len(lattice) != 3:
            print "Structure is not a 3D lattice"
        else:
            assert isinstance(lattice, list), error("")
            assert isinstance(species, list),error("")
            assert isinstance(coords, list),error("")


    @staticmethod
    def rec_lattice_pts_in_scell(scell_matrix=None):
        """Significant portion is taken from the lattice_points_in_supercell function of pymatgen"""

        if scell_matrix is None:
            return
        elif isinstance(scell_matrix, list):
            scell_matrix = np.array(scell_matrix)
        else:
            assert isinstance(scell_matrix, np.ndarray)

        if scell_matrix.size == 3:
            scell_matrix = np.array([[scell_matrix[0], 0, 0], [0, scell_matrix[1], 0], [0, 0, scell_matrix[2]]])
        elif scell_matrix.size == 9:
            scell_matrix = np.array(scell_matrix)
        else:
            print "Invalid scell matrix array size"

        diagonals = np.array(
            [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1],
             [1, 1, 0], [1, 1, 1]])

        d_points = np.dot(diagonals, scell_matrix)

        mins = np.min(d_points, axis=0)
        maxes = np.max(d_points, axis=0) + 1

        ar = np.arange(mins[0], maxes[0])[:, None] * \
             np.array([1, 0, 0])[None, :]
        br = np.arange(mins[1], maxes[1])[:, None] * \
             np.array([0, 1, 0])[None, :]
        cr = np.arange(mins[2], maxes[2])[:, None] * \
             np.array([0, 0, 1])[None, :]

        all_points = ar[:, None, None] + br[None, :, None] + cr[None, None, :]
        all_points = all_points.reshape((-1, 3))

        frac_points = np.dot(all_points, np.linalg.inv(scell_matrix))

        tvects = frac_points[np.all(frac_points < 1 - 1e-10, axis=1)
                             & np.all(frac_points >= -1e-10, axis=1)]
        assert len(tvects) == round(abs(np.linalg.det(scell_matrix)))

        # Tolerance
        tol = 1.0e-6
        tvects[abs(tvects) < tol] = 0.

        return tvects

    def get_lattice(self):
        return self.lattice

    def get_species(self):
        return self.species

    def get_coords(self):
        return self.coords

    def make_scell(self,scell_matrix):
        """Make supercell using 3x1 or 3x3 scell_matrix"""

        if isinstance(scell_matrix,list):
            scell_matrix=np.array(scell_matrix)
        else:
            assert isinstance(scell_matrix, np.ndarray)
            print "Wrong scell format in make_scell"

        if scell_matrix.size == 3:
            scell_matrix = np.array(((scell_matrix[0], 0, 0), (0, scell_matrix[1], 0), (0, 0, scell_matrix[2])))
        elif scell_matrix.size == 9:
            scell_matrix = np.array(scell_matrix)
        else:
            print "Invalid scell matrix array size"

        new_lattice = np.matmul(scell_matrix, np.array(self.get_lattice()))

        # Get reciprocal lattice points
        rec_lattice_pts = np.array(self.rec_lattice_pts_in_scell(scell_matrix))
        real_lattice_pts = []

        for i, lat_pt in enumerate(rec_lattice_pts):
            lat_pt = np.array(lat_pt[0]*new_lattice[0] + lat_pt[1]*new_lattice[1] + lat_pt[2]*new_lattice[2])
            lat_pt = lat_pt.reshape(1,3)

            if i is 0:
                real_lattice_pts=lat_pt
            else:
                real_lattice_pts = np.concatenate((real_lattice_pts, lat_pt) ,axis=0)

        coords = np.array(self.get_coords())
        species = np.array(self.get_species())

        num_coords = coords.shape[0]

        new_coords = []
        new_species = []

        for i, real_lattice_pt in enumerate(real_lattice_pts):

            temp_coords = coords + np.tile(real_lattice_pt, (num_coords,1))
            temp_species = species.copy()

            lattice_origin = np.array([0,0,0],dtype=float)

            if np.array_equal(real_lattice_pt, lattice_origin):
                new_coords = coords.copy()
                new_species = species.copy()
            else:
                new_coords = np.append(new_coords, temp_coords, axis=0)
                new_species = np.append(new_species, temp_species,axis=0)

        #Sort species and coordinates
        species_sorted_index = new_species.argsort()
        new_species = new_species[species_sorted_index]
        new_coords = new_coords[species_sorted_index]

        #Convert to list for Structure class
        new_lattice=new_lattice.tolist()
        new_species=new_species.tolist()
        new_coords=new_coords.tolist()

        return Structure(lattice=new_lattice, species=new_species,coords=new_coords)

    def get_supercells(self, scell_list):
        # Find the most spherical diagonal supercell
        result = dict()
        axes = self.lattice
        a = np.matrix(axes)
        b_keep = []

        if isinstance(scell_list, int):
            scell_list = [scell_list]

        for scellsize in scell_list:
            max_lattice = 100000
            sa_range = (range(-scellsize, scellsize + 1, 1))
            sa_range = [x for x in sa_range if x != 0]
            for sa in sa_range:
                sb_range = (range(-abs(scellsize / sa), abs(scellsize / sa) + 1))
                sb_range = [x for x in sb_range if x != 0]
                if (isinstance(item, int) for item in sb_range):
                    for sb in sb_range:
                        sc_range = (range(-abs(scellsize / (sa * sb)), abs(scellsize / (sa * sb)) + 1))
                        sc_range = [x for x in sc_range if x != 0]
                        if (isinstance(item, int) for item in sc_range):
                            for sc in sc_range:
                                if abs(sa * sb * sc) == scellsize:
                                    b = np.array([sa, sb, sc])

                                    b_tmp = np.multiply(a, b)
                                    b_tmp = np.multiply(b_tmp, b_tmp)
                                    b_tmp = np.sum(b_tmp, axis=1)
                                    b_tmp = np.sqrt(b_tmp)
                                    max_b = np.max(b_tmp)

                                    if max_b <= max_lattice:
                                        max_lattice = max_b
                                        b_keep = b

            result.update({"radius." + str(abs(sa * sb * sc)): b_keep.tolist()})
        return result

    def gen_real_kpts_lattice(self,scell_matrix):
        """Generates gamma centered reciprocal lattice on real lattice for DMC calculations with scell_matrix
        Currently works successfully with diagonal lattices, but non-diagonal lattices must be tested to make sure"""

        if isinstance(scell_matrix,list):
            scell_matrix=np.array(scell_matrix)
        else:
            assert isinstance(scell_matrix, np.ndarray)
            print "Wrong scell format in make_scell"

        if scell_matrix.size == 3:
            scell_matrix = np.array(((scell_matrix[0], 0, 0), (0, scell_matrix[1], 0), (0, 0, scell_matrix[2])))
        elif scell_matrix.size == 9:
            scell_matrix = np.array(scell_matrix)
        else:
            print "Invalid scell matrix array size"

        #Reciprocal lattice points inside the Brillouin zone for supercell lattice vectors
        rec_lattice_pts = np.array(self.rec_lattice_pts_in_scell(scell_matrix))

        #Real grid for [1,1,1] unit cell
        unitcell_k_grid = np.array([[0.,0.,0.],[0.5,0.,0.],[0.,0.5,0.],[0.,0.,0.5],[0.5,0.5,0.],[0.5,0.,0.5],[0.,0.5,0.5],[0.5,0.5,0.5]])
        unitcell_k_weights = np.array([[0.125],[0.125],[0.125],[0.125],[0.125],[0.125],[0.125],[0.125]])

        #Scale the grid for supercells
        scaled_unitcell_k_grid = (np.dot(unitcell_k_grid, np.linalg.inv(scell_matrix)))
        scell_size = np.linalg.det(scell_matrix)
        scaled_unitcell_k_weights = unitcell_k_weights/scell_size

        rec_grid = []
        rec_weights =[]

        for rec_lat_pt in rec_lattice_pts:

            lattice_origin = np.array([0, 0, 0], dtype=float)

            if np.array_equal(rec_lat_pt, lattice_origin):
                rec_grid = scaled_unitcell_k_grid
                rec_weights = scaled_unitcell_k_weights
            else:
                rec_grid = np.vstack((rec_grid, rec_lat_pt + scaled_unitcell_k_grid))
                rec_weights = np.vstack((rec_weights, scaled_unitcell_k_weights))

        return np.hstack((rec_grid,rec_weights))

    def gen_complex_kpts_lattice(self,scell_matrix):
        """Lattice based on monkhorst-pack kgrid"""


        if isinstance(scell_matrix,list):
            scell_matrix=np.array(scell_matrix)
        else:
            assert isinstance(scell_matrix, np.ndarray)
            print "Wrong scell format in make_scell"

        if scell_matrix.size == 3:
            scell_matrix = np.array(((scell_matrix[0], 0, 0), (0, scell_matrix[1], 0), (0, 0, scell_matrix[2])))
        elif scell_matrix.size == 9:
            scell_matrix = np.array(scell_matrix)
        else:
            print "Invalid scell matrix array size"

        #Reciprocal lattice points inside the Brillouin zone for supercell lattice vectors
        rec_lattice_pts = np.array(self.rec_lattice_pts_in_scell(scell_matrix))
        #Real grid for [1,1,1] unit cell

        a=vasp.Vasp(struct=self)

        vasp.Vasp.write_poscar(a)
        f = open('PRECALC', 'w')
        f.write("MINDISTANCE={0}\n".format(self.mindistance))
        f.write("INCLUDEGAMMA=AUTO\n")
        f.close()



        process = subprocess.Popen(
            "curl -s 'http://muellergroup.jhu.edu:8080/PreCalcServer/PreCalcServlet?format=vasp&messagelist=TRUE&clientversion=C2016.06.06' --form fileupload=@PRECALC --form fileupload=@POSCAR",
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = process.communicate()


        os.remove("POSCAR")
        os.remove("PRECALC")
        match = False
        unitcell_k_grid = np.array([-1,-1,-1], dtype=float)
        unitcell_k_weights = np.array([-1], dtype=float)

        for line in out.split(os.linesep):
            if re.match('Fractional', line):
                match = True
            elif re.match('Tetrahedra', line):
                match = False
            elif match:
                line = line.split()


                if np.array_equal(unitcell_k_grid, np.array([-1,-1,-1])):

                    unitcell_k_grid = np.array([float(line[0]), float(line[1]), float(line[2])])
                else:
                    unitcell_k_grid=np.vstack([unitcell_k_grid, np.array([float(line[0]), float(line[1]), float(line[2])])])


                if np.array_equal(unitcell_k_weights, np.array([-1])):
                    unitcell_k_weights=np.array([float(line[3])])
                else:
                    unitcell_k_weights = np.vstack([unitcell_k_weights, np.array([float(line[3])])])


        #Scale the grid for supercells
        unitcell_k_weights = unitcell_k_weights / unitcell_k_weights.sum()
        scaled_unitcell_k_grid = (np.dot(unitcell_k_grid, np.linalg.inv(scell_matrix)))
        scell_size = np.linalg.det(scell_matrix)
        scaled_unitcell_k_weights = unitcell_k_weights/scell_size

        rec_grid = []
        rec_weights =[]

        for rec_lat_pt in rec_lattice_pts:

            lattice_origin = np.array([0, 0, 0], dtype=float)

            if np.array_equal(rec_lat_pt, lattice_origin):
                rec_grid = scaled_unitcell_k_grid
                rec_weights = scaled_unitcell_k_weights
            else:
                rec_grid = np.vstack((rec_grid, rec_lat_pt + scaled_unitcell_k_grid))
                rec_weights = np.vstack((rec_weights, scaled_unitcell_k_weights))

        return np.hstack((rec_grid,rec_weights))