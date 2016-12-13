
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

id = '2'

curr_dir = os.getcwd()
os.system('mkdir result')
result = '{0}/result/result-{1}.txt'.format(curr_dir,id)
os.system('rm {0}'.format(result))

save(result, 'dhcp calculation')

OPTIONS = range(1,20) #np.linspace(0.90, 1.00, 11)

print OPTIONS

dhcp_volumes = []
dhcp_energies = []

for opt in OPTIONS:
    a0 = 3.60 * 0.9 / np.sqrt(2)
    c0 = np.sqrt(8 / 3.0) * a0
    hcp = bulk('Fe', 'hcp', a=a0, c=c0)

    dhcp = hcp*(1,1,2)

    dhcp[3].position[0] = dhcp.get_cell()[0][0]/2
    dhcp[3].position[1] = dhcp[1].position[1]/2

    dhcp.set_tags([1, 1, 1, 1])

    alloys = []
    alloys.append(Alloy(1, 'Fe', 0.4, 1.0))
    alloys.append(Alloy(1, 'Mn', 0.1, 1.0))
    alloys.append(Alloy(1, 'Fe', 0.4, -1.0))
    alloys.append(Alloy(1, 'Mn', 0.1, -1.0))

    calc = EMTO()
    calc.set(dir='work-{0}/opt-{1}'.format(id, opt),
             lat=opt,
             ncpa=20,
             amix=0.05,
             afm='F',  # ferromagnetic calculation
             kpts=[13, 13, 7]
             )
    calc.set_alloys(alloys)

    dhcp.set_calculator(calc)
    dhcp_p = dhcp.get_potential_energy() / dhcp.get_number_of_atoms()
    dhcp_v = dhcp.get_volume() / dhcp.get_number_of_atoms()

    print dhcp_v, dhcp.get_number_of_atoms()

    dhcp_volumes.append(dhcp_v)
    dhcp_energies.append(dhcp_p)

    save(result, 'dhcp result : {0} {1} {2}'.format(opt, dhcp_v, dhcp_p))

    save(result, '------------------------')

eos = EquationOfState(dhcp_volumes, dhcp_energies)
bct_v0, bct_e0, bct_B = eos.fit()
eos.plot('eos-dhcp-{0}.png'.format(id))
os.system('mv eos-dhcp-{0}.png'.format(id))
