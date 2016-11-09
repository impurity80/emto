
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

id = 1
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

cell = atoms.get_cell()

OPTIONS = np.linspace(0.96, 1.02, 7)

energies = []
volumes = []

alloys = []
alloys.append(Alloy(1, 'Fe', 1.0, 1.0))
alloys.append(Alloy(2, 'Fe', 1.0, -1.0))

for opt in OPTIONS:

    atoms.set_cell([opt*cell[0],opt*cell[1],opt*cell[2]],scale_atoms=True)

    calc = EMTO()
    calc.set(dir='work-{1}/opt-{0:0.4f}'.format(opt,id),
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
plt.show('eos-{0}.png'.format(id))

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
