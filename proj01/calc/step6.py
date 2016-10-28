
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *
from ase.lattice import bulk

atoms = bulk('Fe', 'bcc', a=2.8)
atoms.set_tags([1])
atoms.set_initial_magnetic_moments([1])

alloys = []
alloys.append(Alloy(1, 'Cr', 0.15, 0.0))
alloys.append(Alloy(1, 'Ni', 0.15, 0.0))
alloys.append(Alloy(1, 'Fe', 0.35, 1.0))
alloys.append(Alloy(1, 'Fe', 0.35, -1.0))

calc = EMTO()
calc.set_alloys(alloys)

calc.set(dir='work-5',
         lat=3, # bcc lattice calculation
         amix=0.05,
         afm='F',
         kpts=[5,5,5])

atoms.set_calculator(calc)
p = atoms.get_potential_energy()
v = atoms.get_volume()

print v, p , 'eV'

