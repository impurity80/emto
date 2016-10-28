
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *
from ase.lattice import bulk

atoms = bulk('Cu', 'fcc', a=3.65)
atoms.set_tags([1])

alloys = []
alloys.append(Alloy(1, 'Cr', 0.15, 0.0))
alloys.append(Alloy(1, 'Ni', 0.15, 0.0))
alloys.append(Alloy(1, 'Fe', 0.35, 1.0))
alloys.append(Alloy(1, 'Fe', 0.35, -1.0))

calc = EMTO()
calc.set(dir='work-1',
         amix=0.05,
         afm='F',
         lat=2,
         kpts=[1,13,1])

calc.set_alloys(alloys)

atoms.set_calculator(calc)
p = atoms.get_potential_energy()
v = atoms.get_volume()

print v, p , 'eV'

