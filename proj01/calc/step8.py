
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *
from ase.lattice import bulk

fcc = bulk('Fe', 'fcc', a=3.60)
fcc.set_tags([1])

a0 = 3.60/np.sqrt(2)
c0 = np.sqrt(8/3.0)*a0
hcp = bulk('Fe', 'hcp', a=a0, c=c0)
hcp.set_tags([1,1])

#bcc = bulk('Fe', 'bcc', a=2.80)
#bcc.set_tags([1])
# bcc.set_initial_magnetic_moments([1])

alloys = []
alloys.append(Alloy(1, 'Cr', 0.15, 0.0))
alloys.append(Alloy(1, 'Ni', 0.15, 0.0))
alloys.append(Alloy(1, 'Fe', 0.35, 1.0))
alloys.append(Alloy(1, 'Fe', 0.35, -1.0))

#alloys.append(Alloy(2, 'Cr', 0.01, 0.0))
#alloys.append(Alloy(2, 'Ni', 0.01, 0.0))
#alloys.append(Alloy(2, 'Fe', 0.98, 1.0))

calc = EMTO()
calc.set_alloys(alloys)

#calc.set(dir='work8/bcc',
#         ncpa=20,
#         amix=0.05,
#         afm='F',
#         lat=3,
#         kpts=[5,5,5])
#bcc.set_calculator(calc)
#bcc_energy = bcc.get_potential_energy()
#bcc_volume = bcc.get_volume()

calc.set(dir='work8/hcp',
         lat=4,
         ncpa=20,
         amix=0.05,
         afm='F',
         kpts=[5, 5, 5])

hcp.set_calculator(calc)
hcp_energy = hcp.get_potential_energy()/2.0
hcp_volume = hcp.get_volume()

print 'hcp_energy : ', hcp_energy

calc.set(dir='work8/fcc',
         lat=2,
         ncpa=20,
         amix=0.05,
         afm='F',
         kpts=[5, 5, 5])
fcc.set_calculator(calc)
fcc_energy = fcc.get_potential_energy()
fcc_volume = fcc.get_volume()

print 'fcc_energy : ', fcc_energy

print 'delta_e: ', hcp_energy-fcc_energy


