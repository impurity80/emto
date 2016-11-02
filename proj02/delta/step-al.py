
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

id = 'al'
curr_dir = os.getcwd()
os.system('mkdir result')
result = '{0}/result/result-{1}.txt'.format(curr_dir,id)
os.system('rm {0}'.format(result))

save(result, 'delta calculation {0}'.format(id))

a = 3.6
fcc = bulk('Fe', 'fcc', a=3.6)
fcc.set_tags([1])

a0 = 3.60*0.98/np.sqrt(2)
c0 = np.sqrt(8/3.0)*a0
hcp = bulk('Fe', 'hcp', a=a0, c=c0)
hcp.set_tags([1,1])

OPTIONS = np.linspace(0, 0.1, 11)

save(result, OPTIONS)

energies = []
volumes = []

for opt in OPTIONS:

    alloys = []

    fe = 1.0-0.2-opt

    alloys.append(Alloy(1, 'Mn', 0.2, 0.0))
    alloys.append(Alloy(1, 'Al', opt, 0.0))
    alloys.append(Alloy(1, 'Fe', fe/2.0, 1.0))
    alloys.append(Alloy(1, 'Fe', fe/2.0, -1.0))

    calc = EMTO()
    calc.set(dir='work-{1}/opt-{0:0.4f}/hcp'.format(opt,id),
             lat = 4, # hcp
             ncpa=20,
             amix=0.05,
             afm='F', # ferromagnetic calculation
             kpts=[5,5,5],
             )
    calc.set_alloys(alloys)

    hcp.set_calculator(calc)
    hcp_p = hcp.get_potential_energy()/2
    hcp_v = hcp.get_volume()/2

    save(result, 'hcp result : {0} {1} {2}'.format(opt,hcp_v,hcp_p))

    calc = EMTO()
    calc.set(dir='work-{1}/opt-{0:0.4f}/fcc'.format(opt,id),
             lat = 2, # fcc
             ncpa=20,
             amix=0.05,
             afm='F', # ferromagnetic calculation
             kpts=[5,5,5],
             )
    calc.set_alloys(alloys)

    fcc.set_calculator(calc)
    fcc_p = fcc.get_potential_energy()
    fcc_v = fcc.get_volume()

    save(result, 'fcc result : {0} {1} {2}'.format(opt,fcc_v,fcc_p))

    save(result, 'energy differece : {0}'.format(hcp_p - fcc_p))
    energies.append(hcp_p-fcc_p)

    save(result, '------------------------')

    #    volumes.append(4.0/3.0*pi*(opt*0.529177)**3.0)

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
