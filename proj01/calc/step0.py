
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *

bulk = FaceCenteredCubic('Fe', latticeconstant = 3.6)

calc = EMTO()

os.system('export OMP_NUM_THREADS=1')
os.system('export OMP_STACKSIZE=400m')
os.system('mkdir work')
os.chdir('work')

bulk.set_calculator(calc)
bulk.get_potential_energy()

os.chdir('..')

# view(bulk)
