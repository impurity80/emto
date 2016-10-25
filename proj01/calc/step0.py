
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *

atoms = FaceCenteredCubic('Cu', latticeconstant = 3.65, directions=[[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]])
# bulk.set_initial_magnetic_moments([5.0])

# alloy = alloy + Atom('Fe', [0.2, 0.0, 0.0], magmom=-5.0, index=0)
# view(alloy)

alloys = []
alloys.append(Alloy(0, 'Fe', 0.35, 1.0))
alloys.append(Alloy(0, 'Fe', 0.35, -1.0))
alloys.append(Alloy(0, 'Cr', 0.15, 0.0))
alloys.append(Alloy(0, 'Ni', 0.15, 0.0))

calc = EMTO()
calc.set(dir='work-1',
         ncpa=20,
         mnta=4,
         amix=0.05,
         afm='F')

calc.set_alloys(alloys)

atoms.set_calculator(calc)
p = atoms.get_potential_energy()
v = atoms.get_volume()

print v, p , 'eV'

# bulk.read_energy_PBE()


# view(bulk)
