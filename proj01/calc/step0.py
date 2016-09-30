
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *

bulk = FaceCenteredCubic('Fe', latticeconstant = 3.6)

calc = EMTO()

bulk.set_calculator(calc)
bulk.get_potential_energy()

# view(bulk)
