
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *
from ase.lattice import bulk

# atoms = FaceCenteredCubic('Cu', latticeconstant = 3.65, directions=[[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]])

atoms = bulk('Cu', 'fcc', a=3.65)
#atoms.set_tags([0])
#a0 = 3.65/np.sqrt(2)
#c0 = np.sqrt(8/3.0)*a0
#atoms = bulk('Cu', 'hcp', a=a0, c=c0)
#atoms.set_tags([0,1])

# print atoms.get_cell()


#cell = [[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]]

#atoms.set_cell(cell)

#view(atoms)

#a0 = 3.6/np.sqrt(2)
#c0 = np.sqrt(8/3.0)*a0
#atoms = bulk('Cu', 'hcp', a=a0, c=c0)

# bulk.set_initial_magnetic_moments([5.0])

# alloy = alloy + Atom('Fe', [0.2, 0.0, 0.0], magmom=-5.0, index=0)
# view(alloy)

#alloys = []
#alloys.append(Alloy(0, 'Cr', 0.15, 0.0))
#alloys.append(Alloy(0, 'Ni', 0.15, 0.0))
#alloys.append(Alloy(0, 'Fe', 0.35, 1.0))
#alloys.append(Alloy(0, 'Fe', 0.35, -1.0))

calc = EMTO()
calc.set(dir='work-0',
         ncpa=20,
         mnta=4,
         amix=0.05,
         afm='F',
         lat=2,
         kpts=[1,13,1])

#calc.set_alloys(alloys)

atoms.set_calculator(calc)
p = atoms.get_potential_energy()
v = atoms.get_volume()

print v, p , 'eV'

# bulk.read_energy_PBE()


# view(bulk)
