
import csv
import os
from ase import Atom, Atoms
from ase.lattice.cubic import *
from ase.visualize import view
from numpy import *
from emto import *
from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice import bulk

name = 'dlm'

curr_dir = os.getcwd()
os.system('mkdir eos')
os.system('mkdir result')
result = '{0}/result/result-{1}.txt'.format(curr_dir,name)
os.system('rm {0}'.format(result))

result_all = '{0}/result/result_summary-{1}.csv'.format(curr_dir,name)
os.system('rm {0}'.format(result_all))

save(result, 'delta calculation {0}'.format(name))
save(result_all, 'delta calculation {0}'.format(name))

csvfile = open('mole.csv', 'rb')
buffer = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)


for row in buffer:
    id = int(row[0])

    c = row[1]
    mn = round(row[2]/2.0, 3 )*2.0
    ni = round(row[3]/2.0, 3 )*2.0
    cr = round(row[4]/2.0, 3 )*2.0
    al = round(row[5]/2.0, 3 )*2.0
    si = round(row[6]/2.0, 3 )*2.0
    mo = round(row[7]/2.0, 3 )*2.0
    co = round(row[8]/2.0, 3 )*2.0
    cu = round(row[9]/2.0, 3 )*2.0
    nb = round(row[10]/2.0, 3 )*2.0
    ti = round(row[11]/2.0, 3 )*2.0
    v = round(row[12]/2.0, 3 )*2.0
    w = round(row[13]/2.0, 3 )*2.0

    print row
    print mn, ni, cr

    fe = 1-mn-ni-cr-al-si-mo-co-cu-nb-ti-v-w

    save(result, 'alloy id {0}'.format(id))

    OPTIONS = np.linspace(0.98, 1.02, 9)

    volumes = []
    energies = []

    save(result, 'nonmagnetic calculate {0}'.format(id))

    for opt in OPTIONS:
        l = 3.59 * opt
        fcc = bulk('Fe', 'fcc', a=l)

        fcc.set_tags([1])

        alloys = []
        alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
        alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))

        if mn > 1e-7:
            alloys.append(Alloy(1, 'Mn', mn / 2, 1.0))
            alloys.append(Alloy(1, 'Mn', mn / 2, -1.0))
        if ni > 1e-7:
            alloys.append(Alloy(1, 'Ni', ni / 2, 1.0))
            alloys.append(Alloy(1, 'Ni', ni / 2, -1.0))
        if cr > 1e-7:
            alloys.append(Alloy(1, 'Cr', cr / 2, 1.0))
            alloys.append(Alloy(1, 'Cr', cr / 2, -1.0))
        if al > 1e-7:
            alloys.append(Alloy(1, 'Al', al / 2, 1.0))
            alloys.append(Alloy(1, 'Al', al / 2, -1.0))
        if si > 1e-7:
            alloys.append(Alloy(1, 'Si', si / 2, 1.0))
            alloys.append(Alloy(1, 'Si', si / 2, -1.0))
        if mo > 1e-7:
            alloys.append(Alloy(1, 'Mo', mo / 2, 1.0))
            alloys.append(Alloy(1, 'Mo', mo / 2, -1.0))
        if co > 1e-7:
            alloys.append(Alloy(1, 'Co', co / 2, 1.0))
            alloys.append(Alloy(1, 'Co', co / 2, -1.0))
        if cu > 1e-7:
            alloys.append(Alloy(1, 'Cu', cu / 2, 1.0))
            alloys.append(Alloy(1, 'Cu', cu / 2, -1.0))
        if nb > 1e-7:
            alloys.append(Alloy(1, 'Nb', nb / 2, 1.0))
            alloys.append(Alloy(1, 'Nb', nb / 2, -1.0))
        if ti > 1e-7:
            alloys.append(Alloy(1, 'Ti', ti / 2, 1.0))
            alloys.append(Alloy(1, 'Ti', ti / 2, -1.0))
        if v > 1e-7:
            alloys.append(Alloy(1, 'V', v / 2, 1.0))
            alloys.append(Alloy(1, 'V', v / 2, -1.0))
        if w > 1e-7:
            alloys.append(Alloy(1, 'W', w / 2, 1.0))
            alloys.append(Alloy(1, 'W', w / 2, -1.0))

        calc = EMTO()
        calc.set(dir='work/{1}/alloy-{2}/opt-{0:0.4f}'.format(opt, name, id),
                 lat=2,
                 ncpa=20,
                 amix=0.05,
                 afm='F',
                 kpts=[13, 13, 13]
                 )
        calc.set_alloys(alloys)

        fcc.set_calculator(calc)
        nm_e = fcc.get_potential_energy()
        nm_v = fcc.get_volume()

        if nm_e < -0.001:
            volumes.append(nm_v)
            energies.append(nm_e)

        save(result, '{3} result : {0} {1} {2}'.format(opt, nm_v, nm_e, name))


    print volumes, energies


    temp_volumes = []
    temp_energies = []
    pivot = energies[0]
    for v, e in zip(volumes, energies):
        if e-pivot > -0.04 and e-pivot < 0.00:
            temp_volumes.append(v)
            temp_energies.append(e)

    eos = EquationOfState(temp_volumes, temp_energies)

    v0, e0, B = eos.fit()
    eos.plot('eos/{1}-{0}.png'.format(id,name))

    save(result, '{0} {1} {2} {3}'.format(v0, e0, B, (4.0 * v0) ** (1.0 / 3.0)))
    save(result, OPTIONS)
    save(result, volumes)
    save(result, energies)

    save(result, '------------------------')

    save(result_all, '{0}, {1}, {2}, {3}, {4}, {5}'.format(id, e0, v0, B, volumes, energies ))
#    save(result_all, '{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}, {13}, {14}, {15}, {16}, {17}, {18}, {19}, {20}, {21}, {22} '.format(id, hcp_e0-bct_e0, hcp_e0-fcc_e0, hcp_e0-fccf_e0, fcc_e0-bct_e0, fccf_e0-bct_e0, row, fcc_v0, fcc_e0, fcc_B, fccf_v0, fccf_e0, fccf_B, bct_v0, bct_e0, bct_B, hcp_v0, hcp_e0, hcp_B, fcc_energies, fccf_energies, bct_energies, hcp_energies))
