# Simple-minded reader for NAMD files

"""Reader for NAMD binary coordinate (and velocities) format."""


import numpy as np
import simtk.unit as u
from simtk.unit.quantity import Quantity

class NAMDBin:

    def __init__(self, filename):
        with open(filename, 'rb') as namdbin:
            self.n_atoms = np.fromfile(namdbin, dtype=np.int32, count=1)[0]
            coord_double = np.fromfile(namdbin,
                                       dtype=np.float64,
                                       count=self.n_atoms * 3)
        self.pos = coord_double.reshape(self.n_atoms,3)

    def getPositions(self):
        return Quantity(self.pos, u.angstrom)
