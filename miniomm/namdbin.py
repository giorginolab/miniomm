# Simple-minded reader for NAMD files

"""Reader for NAMD binary coordinate (and velocities) format."""


import numpy as np
import simtk.unit as u
from simtk.unit.quantity import Quantity

# self.data is a Numpy array in Angstroms

class NAMDBin:

    def __init__(self, init_data):
        if type(init_data) == str:
            self.read_file(init_data)
        else:
            self.read_data(init_data)

    def read_file(self, filename):
        with open(filename, 'rb') as namdbin:
            self.n_atoms = np.fromfile(namdbin, dtype=np.int32, count=1)[0]
            coord_double = np.fromfile(namdbin,
                                       dtype=np.float64,
                                       count=self.n_atoms * 3)
        self.data = coord_double.reshape(self.n_atoms, 3)

    def read_data(self, pos):
        X = pos.value_in_unit(u.angstrom)
        self.n_atoms = X.shape[0]
        self.data = X

    def write_file(self,filename):
        with open(filename, "wb") as namdbin:
            np.array(self.n_atoms).astype('int32').tofile(namdbin)
            self.data.astype('float64').tofile(namdbin)

    def getPositions(self):
        return Quantity(self.data, u.angstrom)

    def getVelocities(self):
        # https://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2006-2007/0884.html
        PDBVELFACTOR = 20.45482706 
        return Quantity(self.data * PDBVELFACTOR, u.angstrom/u.picosecond)


