
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
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

print rank, size

name = 'afm'

curr_dir = os.getcwd()
os.system('mkdir eos')
os.system('mkdir result')
result = '{0}/result/result-{1}-{2}.txt'.format(curr_dir,name, rank)
os.system('rm {0}'.format(result))

result_all = '{0}/result/result_summary-{1}.csv'.format(curr_dir,name)
os.system('rm {0}'.format(result_all))

save(result, 'delta calculation {0}'.format(name))
save(result_all, 'delta calculation {0}'.format(name))

csvfile = open('mole.csv', 'rb')
buffer = csv.reader(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)


for row in buffer:
    id = int(row[0])

    if id%size==rank:
        c = row [1]
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

        save(result, 'fm calculate {0}'.format(id))

        for opt in OPTIONS:
            l = 3.59 * opt
            a = l / sqrt(2)
            c = l

            atoms = Atoms('Fe2',
                          scaled_positions=[
                              (0.0, 0.0, 0),
                              (0.5, 0.5, 0.5)],
                          cell=[a, a, c],
                          pbc=(1, 1, 1))

            atoms.set_tags([1, 2])

            alloys = []
            alloys.append(Alloy(1, 'Fe', fe, 1.0))
            alloys.append(Alloy(2, 'Fe', fe, -1.0))

            if mn > 1e-7:
                alloys.append(Alloy(1, 'Mn', mn, 1.0))
                alloys.append(Alloy(2, 'Mn', mn, -1.0))
            if ni > 1e-7:
                alloys.append(Alloy(1, 'Ni', ni, 1.0))
                alloys.append(Alloy(2, 'Ni', ni, -1.0))
            if cr > 1e-7:
                alloys.append(Alloy(1, 'Cr', cr, 1.0))
                alloys.append(Alloy(2, 'Cr', cr, -1.0))
            if al > 1e-7:
                alloys.append(Alloy(1, 'Al', al, 1.0))
                alloys.append(Alloy(2, 'Al', al, -1.0))
            if si > 1e-7:
                alloys.append(Alloy(1, 'Si', si, 1.0))
                alloys.append(Alloy(2, 'Si', si, -1.0))
            if mo > 1e-7:
                alloys.append(Alloy(1, 'Mo', mo, 1.0))
                alloys.append(Alloy(2, 'Mo', mo, -1.0))
            if co > 1e-7:
                alloys.append(Alloy(1, 'Co', co, 1.0))
                alloys.append(Alloy(2, 'Co', co, -1.0))
            if cu > 1e-7:
                alloys.append(Alloy(1, 'Cu', cu, 1.0))
                alloys.append(Alloy(2, 'Cu', cu, -1.0))
            if nb > 1e-7:
                alloys.append(Alloy(1, 'Nb', nb, 1.0))
                alloys.append(Alloy(2, 'Nb', nb, -1.0))
            if ti > 1e-7:
                alloys.append(Alloy(1, 'Ti', ti, 1.0))
                alloys.append(Alloy(2, 'Ti', ti, -1.0))
            if v > 1e-7:
                alloys.append(Alloy(1, 'V', v, 1.0))
                alloys.append(Alloy(2, 'V', v, -1.0))
            if w > 1e-7:
                alloys.append(Alloy(1, 'W', w, 1.0))
                alloys.append(Alloy(2, 'W', w, -1.0))

            calc = EMTO()
            calc.set(dir='work/{1}/alloy-{2}/opt-{0:0.4f}'.format(opt, name, id),
                     lat=6,
                     ncpa=20,
                     amix=0.05,
                     afm='F',
                     kpts=[13, 13, 13]
                     )
            calc.set_alloys(alloys)

            atoms.set_calculator(calc)
            nm_e = atoms.get_potential_energy()/atoms.get_number_of_atoms()
            nm_v = atoms.get_volume()/atoms.get_number_of_atoms()

            if nm_e < -0.001:
                volumes.append(nm_v)
                energies.append(nm_e)

            save(result, '{3} result : {0} {1} {2}'.format(opt, nm_v, nm_e, name))

        print volumes, energies

        temp_volumes = []
        temp_energies = []
        pivot = energies[0]
        for v, e in zip(volumes, energies):
            if e-pivot > -0.06 and e-pivot < 0.00:
                temp_volumes.append(v)
                temp_energies.append(e)

        if len(temp_volumes) > 1:
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
