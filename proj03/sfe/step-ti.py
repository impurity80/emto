
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

id = 'Ti'

curr_dir = os.getcwd()
os.system('mkdir result')
result = '{0}/result/result-{1}.txt'.format(curr_dir,id)
os.system('rm {0}'.format(result))

result_summary = '{0}/result/result_summary-{1}.txt'.format(curr_dir,id)
os.system('rm {0}'.format(result_summary))

save(result, 'SFE calculation {0}'.format(id))
save(result_summary, 'SFE calculation {0} summary'.format(id))

OPTIONS = np.linspace(0.0, 0.1, 11)

for conc in OPTIONS:

    save(result, 'concentration : {0}'.format(conc))

    mn = 0.2
    fe = 1-mn-conc

    l = 3.6
    a = l / sqrt(2)
    c = l

    bct = Atoms('Fe2',
                scaled_positions=[
                    (0.0, 0.0, 0),
                    (0.5, 0.5, 0.5)],
                cell=[a, a, c],
                pbc=(1, 1, 1))

    bct.set_tags([1, 2])

    alloys = []
    alloys.append(Alloy(1, 'Fe', fe, 1.0))
    alloys.append(Alloy(1, 'Mn', mn, 1.0))
    alloys.append(Alloy(1, id, conc, 1.0))

    alloys.append(Alloy(2, 'Fe', fe, -1.0))
    alloys.append(Alloy(2, 'Mn', mn, -1.0))
    alloys.append(Alloy(2, id, conc, -1.0))

    calc = EMTO()
    calc.set(dir='work-{0}/opt-{1:0.4f}/bct'.format(id, conc),
             lat=6,
             ncpa=20,
             amix=0.05,
             afm='F',  # ferromagnetic calculation
             kpts=[13, 13, 13]
             )
    calc.set_alloys(alloys)

    bct.set_calculator(calc)
    bct_p = bct.get_potential_energy() / 2
    bct_v = bct.get_volume() / 2

    save(result, 'fcc result : {0} {1} '.format(bct_v, bct_p))

    a0 = 3.60 / np.sqrt(2)
    c0 = np.sqrt(8 / 3.0) * a0
    # c0 = a0 * 1.585
    hcp = bulk('Fe', 'hcp', a=a0, c=c0)
    hcp.set_tags([1, 1])

    alloys = []
    alloys.append(Alloy(1, 'Mn', mn / 2, 1.0))
    alloys.append(Alloy(1, 'Mn', mn / 2, -1.0))
    alloys.append(Alloy(1, id, conc / 2, 1.0))
    alloys.append(Alloy(1, id, conc / 2, -1.0))
    alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
    alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))

    calc = EMTO()
    calc.set(dir='work-{0}/opt-{1:0.4f}/hcp'.format(id, conc),
             lat=4,
             ncpa=20,
             amix=0.05,
             afm='F',  # ferromagnetic calculation
             kpts=[13, 13, 13]
             )
    calc.set_alloys(alloys)

    hcp.set_calculator(calc)
    hcp_p = hcp.get_potential_energy() / 2
    hcp_v = hcp.get_volume() / 2

    save(result, 'hcp result : {0} {1} '.format(hcp_v, hcp_p))

    dhcp = hcp*(1,1,2)

    dhcp[3].position[0] = dhcp.get_cell()[0][0]/2
    dhcp[3].position[1] = dhcp[1].position[1]/2

    dhcp.set_tags([1, 1, 1, 1])

    alloys = []
    alloys.append(Alloy(1, 'Mn', mn / 2, 1.0))
    alloys.append(Alloy(1, 'Mn', mn / 2, -1.0))
    alloys.append(Alloy(1, id, conc / 2, 1.0))
    alloys.append(Alloy(1, id, conc / 2, -1.0))
    alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
    alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))

    calc = EMTO()
    calc.set(dir='work-{0}/opt-{1:0.4f}/dhcp'.format(id, conc),
             lat=4,
             ncpa=20,
             amix=0.05,
             afm='F',  # ferromagnetic calculation
             kpts=[13, 13, 7]
             )
    calc.set_alloys(alloys)

    dhcp.set_calculator(calc)
    dhcp_p = dhcp.get_potential_energy() / 4
    dhcp_v = dhcp.get_volume() / 4

    save(result, 'dhcp result : {0} {1} '.format(dhcp_v, dhcp_p))

    save(result, '------------------------')

    save(result_summary, '{0} {1} {2} {3} {4} {5}'.format(conc, bct_p, hcp_p, dhcp_p, hcp_p-bct_p, hcp_p+2*dhcp_p-3*hcp_p))

