
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

name = 'files2'

curr_dir = os.getcwd()
# os.system('mkdir eos')
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

    if id%12==2:
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

        FCC_OPTIONS = np.linspace(0.98, 1.02, 9)
        FCCF_OPTIONS = np.linspace(0.995, 1.035, 9)
        BCT_OPTIONS = np.linspace(0.98, 1.02, 9)
        HCP_OPTIONS = np.linspace(0.97, 1.01, 9)

        fcc_volumes = []
        fcc_energies = []
        fcc_volumes2 = []
        fcc_energies2 = []

        fccf_volumes = []
        fccf_energies = []
        fccf_volumes2 = []
        fccf_energies2 = []

        bct_volumes = []
        bct_energies = []
        bct_volumes2 = []
        bct_energies2 = []

        hcp_volumes = []
        hcp_energies = []
        hcp_volumes2 = []
        hcp_energies2 = []

        save(result, 'fccf calculate {0}'.format(id))

        for opt in FCCF_OPTIONS:
            l = 3.59 * opt
            fccf = bulk('Fe', 'fcc', a=l)
            fccf.set_tags([1])

            alloys = []
            alloys.append(Alloy(1, 'Fe', fe, 1.0))

            if mn > 1e-7:
                alloys.append(Alloy(1, 'Mn', mn , 1.0))
            if ni > 1e-7:
                alloys.append(Alloy(1, 'Ni', ni , 1.0))
            if cr > 1e-7:
                alloys.append(Alloy(1, 'Cr', cr , 1.0))
            if al > 1e-7:
                alloys.append(Alloy(1, 'Al', al , 1.0))
            if si > 1e-7:
                alloys.append(Alloy(1, 'Si', si , 1.0))
            if mo > 1e-7:
                alloys.append(Alloy(1, 'Mo', mo , 1.0))
            if co > 1e-7:
                alloys.append(Alloy(1, 'Co', co , 1.0))
            if cu > 1e-7:
                alloys.append(Alloy(1, 'Cu', cu , 1.0))
            if nb > 1e-7:
                alloys.append(Alloy(1, 'Nb', nb , 1.0))
            if ti > 1e-7:
                alloys.append(Alloy(1, 'Ti', ti , 1.0))
            if v > 1e-7:
                alloys.append(Alloy(1, 'V', v , 1.0))
            if w > 1e-7:
                alloys.append(Alloy(1, 'W', w , 1.0))

            calc = EMTO()
            calc.set(dir='work/alloy-{2}/fccf/opt-{0:0.4f}'.format(opt, name, id),
                     lat=2,
                     ncpa=20,
                     amix=0.05,
                     afm='F',  # ferromagnetic calculation
                     kpts=[13, 13, 13]
                     )
            calc.set_alloys(alloys)

            fccf.set_calculator(calc)
            fccf_p = fccf.get_potential_energy()
            fccf_v = fccf.get_volume()

            if fccf_p < -0.001:
                fccf_volumes.append(fccf_v)
                fccf_energies.append(fccf_p)

            save(result, 'fccf result : {0} {1} {2}'.format(opt, fccf_v, fccf_p))

        pivot = fccf_energies[0]
        for b, e in zip(fccf_volumes, fccf_energies) :
            if e-pivot > -0.038 and e-pivot < 0.01:
                fccf_volumes2.append(b)
                fccf_energies2.append(e)

        eos = EquationOfState(fccf_volumes2, fccf_energies2)
        fccf_v0, fccf_e0, fccf_B = eos.fit()
        eos.plot('eos/fccf-{0}.png'.format(id))
        os.system('mv eos/fccf-{0}.png result'.format(id))

        save(result, '{0} {1} {2} {3}'.format(fccf_v0, fccf_e0, fccf_B, (4.0 * fccf_v0) ** (1.0 / 3.0)))
        save(result, FCCF_OPTIONS)
        save(result, fccf_volumes)
        save(result, fccf_energies)

        #---------------------------------------------------------------------------------

        save(result, 'fcc calculate {0}'.format(id))

        for opt in FCC_OPTIONS:
            l = 3.59*opt
            fcc = bulk('Fe', 'fcc', a=l)
            fcc.set_tags([1])

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

            calc = EMTO()
            calc.set(dir='work/alloy-{2}/fcc/opt-{0:0.4f}'.format(opt, name, id),
                     lat=2,
                     ncpa=20,
                     amix=0.05,
                     afm='F',  # ferromagnetic calculation
                     kpts=[13, 13, 13]
                     )
            calc.set_alloys(alloys)

            fcc.set_calculator(calc)
            fcc_p = fcc.get_potential_energy()
            fcc_v = fcc.get_volume()

            if fcc_p < -0.001:
                fcc_volumes.append(fcc_v)
                fcc_energies.append(fcc_p)

            save(result, 'fcc result : {0} {1} {2}'.format(opt,fcc_v,fcc_p))

        pivot = fcc_energies[0]
        for b, e in zip(fcc_volumes, fcc_energies) :
            if e-pivot > -0.038 and e-pivot < 0.01:
                fcc_volumes2.append(b)
                fcc_energies2.append(e)

        eos = EquationOfState(fcc_volumes2, fcc_energies2)
        fcc_v0, fcc_e0, fcc_B = eos.fit()
        eos.plot('eos/fcc-{0}.png'.format(id))
        os.system('mv eos/fcc-{0}.png result'.format(id))

        save(result, '{0} {1} {2} {3}'.format(fcc_v0, fcc_e0, fcc_B, (4.0 * fcc_v0) ** (1.0 / 3.0)))
        save(result, FCC_OPTIONS)
        save(result, fcc_volumes)
        save(result, fcc_energies)


        #---------------------------------------------------------------------------------
        save(result, 'bct calculation {0}'.format(id))

        for opt in BCT_OPTIONS:

            l = 3.59*opt
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
            calc.set(dir='work/alloy-{2}/bct/opt-{0:0.4f}'.format(opt,name,id),
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

            save(result, 'bct result : {0} {1} {2}'.format(opt,bct_v,bct_p))

        pivot = bct_energies[0]
        for b, e in zip(bct_volumes, bct_energies) :
            if e-pivot > -0.038 and e-pivot < 0.01:
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

    #    plt.savefig('step-bct-{0}.png'.format(id))
    #    os.system('mv step-bct-{0}.png result'.format(id))

    # ---------------------------------------------------------------------------------------------
        save(result, 'hcp calculation {0}'.format(id))

        for opt in HCP_OPTIONS:
            a0 = 3.59 * opt / np.sqrt(2)
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
            calc.set(dir='work/alloy-{2}/hcp/opt-{0:0.4f}'.format(opt, name, id),
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

     #   save(result, bct_energies)
     #   save(result, hcp_energies)
     #   save(result, '{0} {1} {2}'.format(bct_v0, bct_e0, bct_B))
     #   save(result, '{0} {1} {2}'.format(hcp_v0, hcp_e0, hcp_B))

        save(result, '------------------------')

        save(result_all, '{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}, {13}, {14}, {15}, {16}, {17}, {18}, {19}, {20}, {21}, {22} '.format(id, hcp_e0-bct_e0, hcp_e0-fcc_e0, hcp_e0-fccf_e0, fcc_e0-bct_e0, fccf_e0-bct_e0, row, fcc_v0, fcc_e0, fcc_B, fccf_v0, fccf_e0, fccf_B, bct_v0, bct_e0, bct_B, hcp_v0, hcp_e0, hcp_B, fcc_energies, fccf_energies, bct_energies, hcp_energies))
