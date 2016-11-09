
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

id = 3

curr_dir = os.getcwd()
os.system('mkdir result')
result = '{0}/result/result-{1}.txt'.format(curr_dir,id)
os.system('rm {0}'.format(result))

result_all = '{0}/result/result_all-{1}.txt'.format(curr_dir,id)
os.system('rm {0}'.format(result_all))

save(result, 'delta calculation {0}'.format(id))

BCT_OPTIONS = np.linspace(0.97, 1.03, 3)
HCP_OPTIONS = np.linspace(0.97, 1.03, 3)

bct_volumes = []
bct_energies = []
hcp_volumes = []
hcp_energies = []

mn = 0.2
fe = 1-mn

save(result, 'bct calculation {0}'.format(id))

for opt in BCT_OPTIONS:

    l = 3.6*opt
    a = l/sqrt(2)
    c = l

    bct = Atoms('Fe2',
                 scaled_positions=[
                     (0.0, 0.0, 0),
                     (0.5, 0.5, 0.5)],
                 cell=[a, a, c],
                 pbc=(1, 1, 1))

    bct.set_tags([1,2])

    alloys = []
    alloys.append(Alloy(1, 'Fe', fe, 1.0))
    alloys.append(Alloy(1, 'Mn', mn, 1.0))
    alloys.append(Alloy(2, 'Fe', fe, -1.0))
    alloys.append(Alloy(2, 'Mn', mn, -1.0))

    calc = EMTO()
    calc.set(dir='work-{1}/opt-{0:0.4f}/fcc'.format(opt,id),
             lat = 6, # fcc
             ncpa=20,
             amix=0.05,
             afm='F', # ferromagnetic calculation
             kpts=[13,13,13]
             )
    calc.set_alloys(alloys)

    bct.set_calculator(calc)
    bct_p = bct.get_potential_energy()/2
    bct_v = bct.get_volume()/2

    bct_volumes.append(bct_v)
    bct_energies.append(bct_p)

    save(result, 'fcc result : {0} {1} {2}'.format(opt,bct_v,bct_p))

    save(result, '------------------------')

eos = EquationOfState(bct_volumes, bct_energies)
bct_v0, bct_e0, bct_B = eos.fit()
eos.plot('eos-{0}.png'.format(id))
os.system('mv eos-{0}.png result'.format(id))
# plt.show('eos-{0}.png'.format(id))

save(result, '{0} {1} {2} {3}'.format(bct_v0, bct_e0, bct_B, (4.0*bct_v0)**(1.0/3.0)))
save(result, BCT_OPTIONS)
save(result, bct_volumes)
save(result, bct_energies)
plt.plot(BCT_OPTIONS, bct_energies)
plt.xlabel('OPTIONS')
plt.ylabel('Energy (eV/atom)')

plt.savefig('step-{0}.png'.format(id))
os.system('mv step-{0}.png result'.format(id))

save(result, 'hcp calculation {0}'.format(id))

for opt in HCP_OPTIONS:
    a0 = 3.60 * opt / np.sqrt(2)
 #   c0 = np.sqrt(8 / 3.0) * a0
    c0 = a0*1.585
    hcp = bulk('Fe', 'hcp', a=a0, c=c0)
    hcp.set_tags([1, 1])

    alloys.append(Alloy(1, 'Mn', opt / 2, 1.0))
    alloys.append(Alloy(1, 'Mn', opt / 2, -1.0))
    alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
    alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))

    calc = EMTO()
    calc.set(dir='work-{1}/opt-{0:0.4f}/fcc'.format(opt, id),
             lat=4,  # fcc
             ncpa=20,
             amix=0.05,
             afm='F',  # ferromagnetic calculation
             kpts=[13, 13, 13]
             )
    calc.set_alloys(alloys)

    hcp.set_calculator(calc)
    hcp_p = hcp.get_potential_energy() / 2
    hcp_v = hcp.get_volume() / 2

    hcp_volumes.append(bct_v)
    hcp_energies.append(bct_p)

    save(result, 'hcp result : {0} {1} {2}'.format(opt, hcp_v, hcp_p))

    save(result, '------------------------')

eos = EquationOfState(hcp_volumes, hcp_energies)
hcp_v0, hcp_e0, hcp_B = eos.fit()
eos.plot('eos-{0}-hcp.png'.format(id))
os.system('mv eos-{0}-hcp.png result'.format(id))
# plt.show('eos-{0}.png'.format(id))

save(result, '{0} {1} {2} {3}'.format(hcp_v0, hcp_e0, hcp_B, (4.0 * hcp_v0) ** (1.0 / 3.0)))
save(result, HCP_OPTIONS)
save(result, hcp_volumes)
save(result, hcp_energies)
plt.plot(HCP_OPTIONS, hcp_energies)
plt.xlabel('OPTIONS')
plt.ylabel('Energy (eV/atom)')

plt.savefig('step-{0}-hcp.png'.format(id))
os.system('mv step-{0}-hcp.png result'.format(id))

save(result, bct_energies)
save(result, hcp_energies)
save(result, '{0} {1} {2}'.format(bct_v0, bct_e0, bct_B))
save(result, '{0} {1} {2}'.format(hcp_v0, hcp_e0, hcp_B))

save(result_all, '{0} | {1} | {2} {3} {4} {5} {6} {7}'.format(bct_energies, hcp_energies, bct_v0, bct_e0, bct_B, hcp_v0, hcp_e0, hcp_B, hcp_e0-bct_e0))

