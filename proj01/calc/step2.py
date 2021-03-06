
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *
from ase.lattice import bulk

a0 = 3.60/np.sqrt(2)
c0 = np.sqrt(8/3.0)*a0
atoms = bulk('Cu', 'hcp', a=a0, c=c0)

calc = EMTO()
calc.set(dir='work-2',
         amix=0.05,
         afm='F',
         lat=4,
         kpts=[1, 13, 1])

atoms.set_calculator(calc)
p = atoms.get_potential_energy()
v = atoms.get_volume()

print v, p , 'eV'

