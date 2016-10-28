
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *
from ase.lattice import bulk

a0 = 3.60/np.sqrt(2)
c0 = np.sqrt(8/3.0)*a0
atoms = bulk('Cu', 'hcp', a=a0, c=c0)
atoms.set_tags([1,1])

alloys = []
alloys.append(Alloy(1, 'Cr', 0.15, 0.0))
alloys.append(Alloy(1, 'Ni', 0.15, 0.0))
alloys.append(Alloy(1, 'Fe', 0.35, 1.0))
alloys.append(Alloy(1, 'Fe', 0.35, -1.0))

calc = EMTO()
calc.set(dir='work-3',
         lat=4,
         amix=0.05,
         afm='F',
         kpts=[1, 13, 1])

calc.set_alloys(alloys)

atoms.set_calculator(calc)
p = atoms.get_potential_energy()
v = atoms.get_volume()

print v, p , 'eV'

