
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *
from ase.lattice import bulk

atoms = bulk('Cu', 'fcc', a=3.65)
atoms.set_initial_magnetic_moments([1.0])

calc = EMTO()
calc.set(dir='work-0',
         lat=2, # fcc lattice calculation
         amix=0.05,
         afm='F',
         kpts=[1,13,1])

atoms.set_calculator(calc)
p = atoms.get_potential_energy()
v = atoms.get_volume()

print v, p , 'eV'

