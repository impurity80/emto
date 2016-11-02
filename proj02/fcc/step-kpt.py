
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

id = 'kpt'
curr_dir = os.getcwd()
os.system('mkdir result')
result = '{0}/result/result-{1}.txt'.format(curr_dir,id)
os.system('rm {0}'.format(result))

save(result, 'fcc calculation')

atoms = bulk('Fe', 'fcc', a=3.6*0.98)
atoms.set_tags([1])

# OPTIONS = np.linspace(0.97, 1.02, 11)
OPTIONS = range(5,20,2)

energies = []
volumes = []

alloys = []
alloys.append(Alloy(1, 'Fe', 0.4, 1.0))
alloys.append(Alloy(1, 'Fe', 0.4, -1.0))
alloys.append(Alloy(1, 'Mn', 0.2, 0.0))

for opt in OPTIONS:

#    atoms.set_cell([opt*cell[0],opt*cell[1],opt*cell[2]],scale_atoms=True)

    calc = EMTO()
    calc.set(dir='work-{1}/opt-{0:0.4f}'.format(opt,id),
             lat = 2, # face centered cubic
             ncpa=20,
             amix=0.05,
             afm='F', # ferromagnetic calculation
             kpts=[opt,opt,opt],
             )
    calc.set_alloys(alloys)

    atoms.set_calculator(calc)
    p = atoms.get_potential_energy()
    v = atoms.get_volume()

    save(result, '{0} {1} {2}'.format(opt,v,p))
    save(result, '------------------------')

    if p < 0 :
        energies.append(p)
        volumes.append(v)

    #    volumes.append(4.0/3.0*pi*(opt*0.529177)**3.0)

save(result, OPTIONS)
save(result, volumes)
save(result, energies)
plt.plot(OPTIONS, energies)
plt.xlabel('OPTIONS')
plt.ylabel('Energy (eV/atom)')

plt.savefig('step-{0}.png'.format(id))
os.system('mv step-{0}.png result'.format(id))
plt.show('step-{0}.png'.format(id))


# view(bulk)

