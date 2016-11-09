
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

id = 0
curr_dir = os.getcwd()
os.system('mkdir result')
result = '{0}/result/result-{1}.txt'.format(curr_dir,id)
os.system('rm {0}'.format(result))

save(result, 'bct-afm calculation')

l = 3.6
a = l/sqrt(2)
c = l

atoms = Atoms('Fe2',
             scaled_positions=[
                 (0.0, 0.0, 0),
                 (0.5, 0.5, 0.5)],
             cell=[a, a, c],
             pbc=(1, 1, 1))

atoms.set_tags([1,2])

view(atoms)

calc = EMTO()

alloys = []
alloys.append(Alloy(1, 'Fe', 1.0, 1.0))
alloys.append(Alloy(2, 'Fe', 1.0, -1.0))

calc.set(dir='work-{0}'.format(id),
         lat = 6, # face centered cubic
         ncpa=20,
         amix=0.05,
         afm='F', # ferromagnetic calculation
         kpts=[13,13,13],
         )
calc.set_alloys(alloys)

atoms.set_calculator(calc)
p = atoms.get_potential_energy()
v = atoms.get_volume()

save(result, '{0} {1} '.format(v,p))
save(result, '------------------------')
