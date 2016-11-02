
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from emto import *
from ase.lattice import bulk

id = 1
curr_dir = os.getcwd()
os.system('mkdir result')
result1 = '{0}/result/result1.txt'.format(curr_dir)
result2 = '{0}/result/result2.txt'.format(curr_dir)

save(result1, 'Energy calculation start')
save(result2, 'cr, ni, fcc, bcc, hcp, bcc-fcc, hcp-fcc')

XOPTIONS = np.linspace(0.0, 0.15, 16)
YOPTIONS = np.linspace(0.0, 0.15, 16)

save(result1, 'XOPTIONS: {0}'.format(XOPTIONS))
save(result1, 'YOPTIONS: {0}'.format(YOPTIONS))

for xopt in XOPTIONS:
    for yopt in YOPTIONS:

        save(result1, 'option: {0} {1}'.format(xopt, yopt))
        cr = xopt
        ni = yopt
        work_dir = 'work-{0}/xopt-{1}/yopt-{2}/'.format(id,xopt,yopt)

        fe = (1.0-cr-ni)/2.0

        a0 = 3.60/np.sqrt(2)
        c0 = np.sqrt(8/3.0)*a0
        hcp = bulk('Fe', 'hcp', a=a0, c=c0)
        hcp.set_tags([1,1])

        fcc = bulk('Fe', 'fcc', a=3.60)
        fcc.set_tags([1])

        bcc = bulk('Fe', 'bcc', a=2.80)
        bcc.set_tags([1])
        # bcc.set_initial_magnetic_moments([1])

        alloys = []
        alloys.append(Alloy(1, 'Cr', cr, 0.0))
        alloys.append(Alloy(1, 'Ni', ni, 0.0))
        alloys.append(Alloy(1, 'Fe', fe, 1.0))
        alloys.append(Alloy(1, 'Fe', fe, -1.0))

        calc = EMTO()
        calc.set_alloys(alloys)

        calc.set(dir= work_dir+'bcc',
                 ncpa=20,
                 amix=0.05,
                 afm='F',
                 lat=3,
                 kpts=[5,5,5])
        bcc.set_calculator(calc)
        bcc_energy = bcc.get_potential_energy()
        bcc_volume = bcc.get_volume()

        save(result1, 'BCC_energy : {0} eV/atom'.format(bcc_energy) )

        calc.set(dir=work_dir+'hcp',
                 ncpa=20,
                 amix=0.05,
                 afm='F',
                 lat=4,
                 kpts=[5, 5, 5])

        hcp.set_calculator(calc)
        hcp_energy = hcp.get_potential_energy()/2.0
        hcp_volume = hcp.get_volume()

        save(result1, 'HCP_energy : {0} eV/atom'.format(hcp_energy) )

        calc.set(dir=work_dir+'fcc',
                 lat=2)
        fcc.set_calculator(calc)
        fcc_energy = fcc.get_potential_energy()
        fcc_volume = fcc.get_volume()

        save(result1, 'FCC_energy : {0} eV/atom'.format(fcc_energy) )

        save(result2, '{0},{1},{2},{3},{4},{5},{6}'.format(cr, ni, fcc_energy, hcp_energy, bcc_energy, hcp_energy-fcc_energy, bcc_energy-fcc_energy))

#print fcc_energy, hcp_energy, bcc_energy
#print hcp_energy-fcc_energy, bcc_energy-fcc_energy


