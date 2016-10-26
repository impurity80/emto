
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

# atoms = FaceCenteredCubic('Fe', latticeconstant = 3.6, directions=[[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]])

a = 3.6
atoms = bulk('Cu', 'fcc', a=3.6)
cell = atoms.get_cell().copy()

OPTIONS = np.linspace(0.98, 1.02, 101)

energies = []
volumes = []

alloys = []
alloys.append(Alloy(0, 'Fe', 0.35, 1.0))
alloys.append(Alloy(0, 'Fe', 0.35, -1.0))
alloys.append(Alloy(0, 'Cr', 0.15, 0.0))
alloys.append(Alloy(0, 'Ni', 0.15, 0.0))

for opt in OPTIONS:

    atoms.set_cell([opt*cell[0],opt*cell[1],opt*cell[2]],scale_atoms=True)

    calc = EMTO()
    calc.set(dir='work-{0:0.4f}'.format(opt),
             ncpa=20,
             mnta=4,
             amix=0.05,
             afm='F', # ferromagnetic calculation
             iprim = 0,
             lat = 2
             )
    calc.set_alloys(alloys)

    atoms.set_calculator(calc)
    p = atoms.get_potential_energy()
    v = atoms.get_volume()
    if p < 0 :
        energies.append(p)
        volumes.append(v)

    #    volumes.append(4.0/3.0*pi*(opt*0.529177)**3.0)

eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
eos.plot('eos.png')

plt.plot(OPTIONS, energies)
plt.xlabel('OPTIONS')
plt.ylabel('Energy (eV/atom)')

plt.savefig('step1.png')
plt.show('step1.png')

# view(bulk)
