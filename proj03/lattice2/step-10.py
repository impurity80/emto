
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

name = 'files10'

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

    if id%16==10:
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

        save(result, 'alloy id  {0}'.format(id))

        BCT_OPTIONS = np.linspace(0.98, 1.02, 9)
        HCP_OPTIONS = np.linspace(0.96, 1.01, 11)

        bct_volumes = []
        bct_energies = []
        bct_volumes2 = []
        bct_energies2 = []

        hcp_volumes = []
        hcp_energies = []
        hcp_volumes2 = []
        hcp_energies2 = []

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
            calc.set(dir='work-{1}/alloy-{2}/bct/opt-{0:0.4f}'.format(opt,name,id),
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

            if bct_p < -0.001:
                bct_volumes.append(bct_v)
                bct_energies.append(bct_p)

            save(result, 'fcc result : {0} {1} {2}'.format(opt,bct_v,bct_p))

            save(result, '------------------------')

        pivot = bct_energies[0]
        for b, e in zip(bct_volumes, bct_energies) :
            if e-pivot > -0.04 and e-pivot < 0.02:
                bct_volumes2.append(b)
                bct_energies2.append(e)

        eos = EquationOfState(bct_volumes2, bct_energies2)
        bct_v0, bct_e0, bct_B = eos.fit()
        eos.plot('eos/bct-{0}.png'.format(id))
        os.system('mv eos/bct-{0}.png result'.format(id))
        # plt.show('eos-{0}.png'.format(id))

        save(result, '{0} {1} {2} {3}'.format(bct_v0, bct_e0, bct_B, (4.0*bct_v0)**(1.0/3.0)))
        save(result, BCT_OPTIONS)
        save(result, bct_volumes)
        save(result, bct_energies)
    #    plt.plot(BCT_OPTIONS, bct_energies)
    #    plt.xlabel('OPTIONS')
    #    plt.ylabel('Energy (eV/atom)')

        plt.savefig('step-bct-{0}.png'.format(id))
        os.system('mv step-bct-{0}.png result'.format(id))

        save(result, 'hcp calculation {0}'.format(id))

        for opt in HCP_OPTIONS:
            a0 = 3.60 * opt / np.sqrt(2)
         #   c0 = np.sqrt(8 / 3.0) * a0
            c0 = a0*1.585
            hcp = bulk('Fe', 'hcp', a=a0, c=c0)
            hcp.set_tags([1, 1])

            alloys = []
            alloys.append(Alloy(1, 'Fe', fe/2, 1.0))
            alloys.append(Alloy(1, 'Fe', fe/2, -1.0))

            if mn > 1e-7:
                alloys.append(Alloy(1, 'Mn', mn/2, 1.0))
                alloys.append(Alloy(1, 'Mn', mn/2, -1.0))
            if ni > 1e-7:
                alloys.append(Alloy(1, 'Ni', ni/2, 1.0))
                alloys.append(Alloy(1, 'Ni', ni/2, -1.0))
            if cr > 1e-7:
                alloys.append(Alloy(1, 'Cr', cr/2, 1.0))
                alloys.append(Alloy(1, 'Cr', cr/2, -1.0))
            if al > 1e-7:
                alloys.append(Alloy(1, 'Al', al/2, 1.0))
                alloys.append(Alloy(1, 'Al', al/2, -1.0))
            if si > 1e-7:
                alloys.append(Alloy(1, 'Si', si/2, 1.0))
                alloys.append(Alloy(1, 'Si', si/2, -1.0))
            if mo > 1e-7:
                alloys.append(Alloy(1, 'Mo', mo/2, 1.0))
                alloys.append(Alloy(1, 'Mo', mo/2, -1.0))
            if co > 1e-7:
                alloys.append(Alloy(1, 'Co', co/2, 1.0))
                alloys.append(Alloy(1, 'Co', co/2, -1.0))
            if cu > 1e-7:
                alloys.append(Alloy(1, 'Cu', cu/2, 1.0))
                alloys.append(Alloy(1, 'Cu', cu/2, -1.0))
            if nb > 1e-7:
                alloys.append(Alloy(1, 'Nb', nb/2, 1.0))
                alloys.append(Alloy(1, 'Nb', nb/2, -1.0))
            if ti > 1e-7:
                alloys.append(Alloy(1, 'Ti', ti/2, 1.0))
                alloys.append(Alloy(1, 'Ti', ti/2, -1.0))
            if v > 1e-7:
                alloys.append(Alloy(1, 'V', v/2, 1.0))
                alloys.append(Alloy(1, 'V', v/2, -1.0))
            if w > 1e-7:
                alloys.append(Alloy(1, 'W', w/2, 1.0))
                alloys.append(Alloy(1, 'W', w/2, -1.0))
        #    alloys = []
        #    alloys.append(Alloy(1, 'Mn', mn / 2, 1.0))
        #    alloys.append(Alloy(1, 'Mn', mn / 2, -1.0))
        #    alloys.append(Alloy(1, id , conc / 2, 1.0))
        #    alloys.append(Alloy(1, id , conc / 2, -1.0))
        #    alloys.append(Alloy(1, 'Fe', fe / 2, 1.0))
        #    alloys.append(Alloy(1, 'Fe', fe / 2, -1.0))

            calc = EMTO()
            calc.set(dir='work-{1}/alloy-{2}/hcp/opt-{0:0.4f}'.format(opt, name, id),
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

            if hcp_p < -0.001 :
                hcp_volumes.append(hcp_v)
                hcp_energies.append(hcp_p)

            save(result, 'hcp result : {0} {1} {2}'.format(opt, hcp_v, hcp_p))

            save(result, '------------------------')

        eos = EquationOfState(hcp_volumes, hcp_energies)
        hcp_v0, hcp_e0, hcp_B = eos.fit()
        eos.plot('eos/hcp-{0}.png'.format(id))
        os.system('mv eos/hcp-{0}.png result'.format(id))
        # plt.show('eos-{0}.png'.format(id))

        save(result, '{0} {1} {2} {3}'.format(hcp_v0, hcp_e0, hcp_B, (4.0 * hcp_v0) ** (1.0 / 3.0)))
        save(result, HCP_OPTIONS)
        save(result, hcp_volumes)
        save(result, hcp_energies)
       # plt.plot(HCP_OPTIONS, hcp_energies)
       # plt.xlabel('OPTIONS')
       # plt.ylabel('Energy (eV/atom)')

       # plt.savefig('step-hcp-{0}.png'.format(id))
       # os.system('mv step-hcp-{0}.png result'.format(id))

        save(result, bct_energies)
        save(result, hcp_energies)
        save(result, '{0} {1} {2}'.format(bct_v0, bct_e0, bct_B))
        save(result, '{0} {1} {2}'.format(hcp_v0, hcp_e0, hcp_B))

        save(result_all, '{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11} '.format(id, row, bct_energies, hcp_energies, bct_v0, bct_e0, bct_B, hcp_v0, hcp_e0, hcp_B, hcp_v0-bct_v0, hcp_e0-bct_e0))
