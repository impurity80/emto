
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

id = 5
curr_dir = os.getcwd()
os.system('mkdir result')
result = '{0}/result/result-{1}.txt'.format(curr_dir,id)
os.system('rm {0}'.format(result))

save(result, 'hcp calculation {0}'.format(id))

a0 = 3.60/np.sqrt(2)
c0 = np.sqrt(8/3.0)*a0*0.98
atoms = bulk('Fe', 'hcp', a=a0, c=c0)
atoms.set_tags([1,1])

cell = atoms.get_cell().copy()

OPTIONS = np.linspace(0.96, 1.02, 7)

energies = []
volumes = []

alloys = []
alloys.append(Alloy(1, 'Fe', 0.5, 1.0))
alloys.append(Alloy(1, 'Fe', 0.5, -1.0))

for opt in OPTIONS:

    atoms.set_cell([opt*cell[0],opt*cell[1],opt*cell[2]],scale_atoms=True)

    calc = EMTO()
    calc.set(dir='work-{1}/opt-{0:0.4f}'.format(opt,id),
             lat = 4, # hcp
             ncpa=20,
             amix=0.05,
             afm='F', # ferromagnetic calculation
             kpts=[13,13,13],
             )
    calc.set_alloys(alloys)

    atoms.set_calculator(calc)
    p = atoms.get_potential_energy()/2
    v = atoms.get_volume()/2

    save(result, '{0} {1} {2}'.format(opt,v,p))
    save(result, '------------------------')

    if p < 0 :
        energies.append(p)
        volumes.append(v)

    #    volumes.append(4.0/3.0*pi*(opt*0.529177)**3.0)


eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
eos.plot('eos-{0}.png'.format(id))
os.system('mv eos-{0}.png result'.format(id))
# plt.show('eos-{0}.png'.format(id))

save(result, '{0} {1} {2} {3}'.format(v0, e0, B, (4.0*v0)**(1.0/3.0)))
save(result, OPTIONS)
save(result, volumes)
save(result, energies)
plt.plot(OPTIONS, energies)
plt.xlabel('OPTIONS')
plt.ylabel('Energy (eV/atom)')

plt.savefig('step-{0}.png'.format(id))
os.system('mv step-{0}.png result'.format(id))
# plt.show('step-{0}.png'.format(id))


# view(bulk)
