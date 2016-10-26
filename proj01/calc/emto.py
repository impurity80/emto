from __future__ import print_function

import os
import sys
import re
from ase.calculators.general import Calculator
from os.path import join, isfile, islink
from numpy import *

import numpy as np

import ase
import ase.io
from ase.utils import devnull

from ase.calculators.singlepoint import SinglePointCalculator

class Alloy():
    def __init__(self, index, symbol, conc, split):
        self.index = index
        self.symbol = symbol
        self.conc = conc
        self.split = split

element_keys = [
    'Al', # 2d metal
    'Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn', # 3d transition metal
    'Y', # 4d transition metal
]

elements = {}
for key in element_keys:
    elements[key] = None

elements['Al'] =  'Iz=  13 Norb=  7 Ion=  0 Config= 3s2_3p1\n' \
                  'n      1  2  2  2  3  3  3\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2\n' \
                  'Occup  2  2  2  4  2  0  1\n' \
                  'Valen  0  0  0  0  1  1  1\n'
elements['Ti'] =  'Iz=  22 Norb= 10 Ion=  0 Config= 3d2_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  2  0  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1  1\n'
elements['V'] =   'Iz=  23 Norb= 10 Ion=  0 Config= 3d3_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  3  0  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1  1\n'
elements['Cr'] =  'Iz=  24 Norb= 10 Ion=  0 Config= 3d5_4s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  1  1\n' \
                  'Valen  0  0  0  0  0  0  0  1  1  1\n'
elements['Mn'] =  'Iz=  25 Norb= 10 Ion=  0 Config= 3d5_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  1  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1  1\n'
elements['Fe'] =  'Iz=  26 Norb= 10 Ion=  0 Config= 3d6_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  2  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1  1\n'
elements['Co'] =  'Iz=  27 Norb= 10 Ion=  0 Config= 3d7_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  3  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1  1\n'
elements['Ni'] =  'Iz=  28 Norb= 10 Ion=  0 Config= 3d8_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  4  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1  1\n'
elements['Cu'] =  'Iz=  29 Norb= 10 Ion=  0 Config= 3d10_4s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  1\n' \
                  'Valen  0  0  0  0  0  0  0  1  1  1\n'
elements['Zn'] =  'Iz=  30 Norb= 10 Ion=  0 Config= 3d10_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1  1\n'
elements['Y'] =  'Iz=  39 Norb= 15 Ion=  0 Config= 4d1_5s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  1  0  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'

common_keys = [
    'dir',
    'lat',
    'nq',
    'iprim',
    'nqr2',
    'nl',
    'lamda',
    'amax',
    'bmax',
]

bmdl_keys = [
]

kgrn_keys = [
    'niter',
    'nlin',
    'nprn',
    'ncpa',
    'nt',
    'mnta',
    'mode',
    'frc',
    'dos',
    'ops',
    'afm',
    'crt',
    'lmaxh',
    'lmaxt',
    'nfi',
    'fixg',
    'shf',
    'sofc',
    'kmsh',
    'ibz',
    'nkx',
    'nky',
    'nkz',
    'fbz',
    'kmsh2',
    'ibz2',
    'nkx2',
    'nky2',
    'nkz2',
    'zmsh',
    'nz1',
    'nz2',
    'nz3',
    'nres',
    'nzd',
    'depth',
    'imagz',
    'eps',
    'elim',
    'amix',
    'efmiz',
    'vmtz',
    'mmom',
    'tole',
    'tolef',
    'tolcpa',
    'tfermi',
   # 'sws',
    'nsws',
    'dsws',
    'alpcpa'
]

class EMTO(Calculator):
    name = 'EMTO'

    def __init__(self, restart=None, **kwargs):

        self.alloys = []

        self.common_params = {}
        self.bmdl_params = {}
        self.kstr_params = {}
        self.shape_params = {}
        self.kgrn_params = {}
        self.kfcd_params = {}

        for key in common_keys:
            self.common_params[key] = None
        for key in kgrn_keys:
            self.kgrn_params[key] = None

        self.energy_pbe = 0

        self.common_params['dir'] = 'work'
    #    self.common_params['nq'] = 1
        self.common_params['lat'] = 1
        self.common_params['iprim'] = 0
        self.common_params['nqr2'] = 0
        self.common_params['nl'] = 7
        self.common_params['lamda'] = 2.50
        self.common_params['amax'] = 4.50
        self.common_params['bmax'] = 4.50
        self.common_params['nghbp'] = 13

        self.kgrn_params['nlin'] = 31
        self.kgrn_params['nprn'] = 0
        self.kgrn_params['ncpa'] = 20
        self.kgrn_params['nt'] = 1
        self.kgrn_params['mnta'] = 4
        self.kgrn_params['mode'] = '3D'
        self.kgrn_params['frc'] = 'N'
        self.kgrn_params['dos'] = 'N'
        self.kgrn_params['ops'] = 'N'
        self.kgrn_params['afm'] = 'F'
        self.kgrn_params['crt'] = 'M'
        self.kgrn_params['lmaxh'] = 8
        self.kgrn_params['lmaxt'] = 4
        self.kgrn_params['nfi'] = 31
        self.kgrn_params['fixg'] = 2
        self.kgrn_params['shf'] = 0
        self.kgrn_params['sofc'] = 'N'
        self.kgrn_params['kmsh'] = 'G'
        self.kgrn_params['ibz'] = 2
        self.kgrn_params['nkx'] = 0
        self.kgrn_params['nky'] = 13
        self.kgrn_params['nkz'] = 0
        self.kgrn_params['fbz'] = 'N'
        self.kgrn_params['kmsh2'] = 'G'
        self.kgrn_params['ibz2'] = 1
        self.kgrn_params['nkx2'] = 4
        self.kgrn_params['nky2'] = 0
        self.kgrn_params['nkz2'] = 51
        self.kgrn_params['zmsh'] = 'C'
        self.kgrn_params['nz1'] = 16
        self.kgrn_params['nz2'] = 16
        self.kgrn_params['nz3'] = 8
        self.kgrn_params['nres'] = 4
        self.kgrn_params['nzd'] = 200
        self.kgrn_params['depth'] = 1.0
        self.kgrn_params['imagz'] = 0.02
        self.kgrn_params['eps'] = 0.200
        self.kgrn_params['elim'] = -1.0
        self.kgrn_params['amix'] = 0.05
        self.kgrn_params['efmix'] = 1.0
        self.kgrn_params['vmtz'] = 0.0
        self.kgrn_params['mmom'] = 0.0
        self.kgrn_params['tole'] = '1.d-08' # 1e-08
        self.kgrn_params['tolef'] = '1.d-07' # 1e-07
        self.kgrn_params['tolcpa'] = '1.d-06' # 1e-06
        self.kgrn_params['tfermi'] = 5000.0
    #    self.kgrn_params['sws'] = 2.69
        self.kgrn_params['nsws'] = 1
        self.kgrn_params['dsws'] = 0.05
        self.kgrn_params['alpcpa'] = 0.6020

    def set(self, **kwargs):
        for key in kwargs:
            if key in self.common_params:
                self.common_params[key] = kwargs[key]
            if key in self.bmdl_params:
                self.bmdl_params[key] = kwargs[key]
            if key in self.kgrn_params:
                self.kgrn_params[key] = kwargs[key]

    def set_kgrn(self, **kwargs):
        for key in kwargs:
            self.kgrn_params[key] = kwargs[key]

    def initialize(self, atoms):
        self.natoms = len(atoms)

        if len(self.alloys) > 0:
            a = 0
    #    self.alloys = []
    #    for atom in atoms:
    #        self.alloys.append(Alloy(atom.index, atom.symbol, 1.0, 1.0))

    def calculation_required(selfself, atoms, quantities):
        return True

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        return self.energy_pbe

    def update(self, atoms):
        self.calculate(atoms)
#        self.calculate(atoms)

    def calculate2(self, atoms):
        os.chdir(self.common_params['dir'])
        self.set_results(atoms)

    def calculate(self, atoms):
        os.system('rm -r {0}'.format(self.common_params['dir']))
        os.system('mkdir {0}'.format(self.common_params['dir']))
        os.chdir(self.common_params['dir'])

        self.initialize(atoms)
        self.write_bmdl(atoms)
        os.system(os.environ['BMDL_CMD'])
#        os.system('time -f " Elapsed time = %E" bmdl < bmdl.dat')
        self.write_kstr(atoms)
        os.system(os.environ['KSTR_CMD'])
#        os.system('time -f " Elapsed time = %E" kstr < kstr.dat')
        self.write_shape(atoms)
        os.system(os.environ['SHAPE_CMD'])
#        os.system('time -f " Elapsed time = %E" shape < shape.dat')
        self.write_kgrn(atoms)
        os.system(os.environ['KGRN_CMD'])
#        os.system('time -f " Elapsed time = %E" kgrn_cpa < kgrn.dat')

        os.system('cp bmdl.mdl common.mdl')
        os.system('cp shape.shp common.shp')
        os.system('cp kstr.tfh common.tfh')
        os.system('cp kgrn.chd kfcd.chd')

        self.write_kfcd(atoms)
        os.system(os.environ['KFCD_CMD'])
#        os.system('time -f " Elapsed time = %E" kfcd_cpa < kfcd.dat')

        self.set_results(atoms)

        os.chdir('..')

    def set_results(self, atoms):
        self.read(atoms)

    def read(self, atoms):
        self.read_energy()

    def write_bmdl(self, atoms):
        bmdl = open('bmdl.dat', 'w')
        bmdl.write('BMDL      HP......=N                               22 Jan 08\n')
        bmdl.write('JOBNAM...=bmdl       MSGL.=  1 NPRN.=  0\n')
        bmdl.write('DIR001=./\n')
        bmdl.write('DIR006=\n')
        bmdl.write('Madelung potential\n')

        bmdl.write('NL.....={:2d}\n'.format(self.common_params['nl']))

        bmdl.write('LAMDA....={:10.2f} '.format(self.common_params['lamda']))
        bmdl.write('AMAX....={:10.2f} '.format(self.common_params['amax']))
        bmdl.write('BMAX....={:10.2f}\n'.format(self.common_params['bmax']))

        bmdl.write('NQ....={:3d} '.format(atoms.get_number_of_atoms()))
        bmdl.write('LAT...={:2d} '.format(self.common_params['lat']))
        bmdl.write('IPRIM.={:2d} '.format(self.common_params['iprim']))
        bmdl.write('NQR2..={:2d}\n'.format(self.common_params['nqr2']))

        if self.common_params['iprim']==0:
            cell = atoms.get_cell()

            a = np.linalg.norm(cell[0])
            b = np.linalg.norm(cell[1])
            c = np.linalg.norm(cell[2])

            bmdl.write('A........={:10.3f} '.format(a/a))
            bmdl.write('B.......={:10.3f} '.format(b/a))
            bmdl.write('C.......={:10.3f}\n'.format(c/a))

            if self.common_params['lat']==2: # fcc
                l = cell[1][0]*2
            elif self.common_params['lat']==3: # bcc
                l = cell[0][0]*2
            else:
                l = cell[0][0]

            cell = atoms.get_cell().copy()/l
            bmdl.write('BSX......={:10.7f} '.format(cell[0][0]))
            bmdl.write('BSY.....={:10.7f} '.format(cell[0][1]))
            bmdl.write('BSZ.....={:10.7f}\n'.format(cell[0][2]))
            bmdl.write('BSX......={:10.7f} '.format(cell[1][0]))
            bmdl.write('BSY.....={:10.7f} '.format(cell[1][1]))
            bmdl.write('BSZ.....={:10.7f}\n'.format(cell[1][2]))
            bmdl.write('BSX......={:10.7f} '.format(cell[2][0]))
            bmdl.write('BSY.....={:10.7f} '.format(cell[2][1]))
            bmdl.write('BSZ.....={:10.7f}\n'.format(cell[2][2]))

            for atom in atoms:
                bmdl.write('QX.......={:10.7f} '.format(atom.position[0]))
                bmdl.write('QY......={:10.7f} '.format(atom.position[1]))
                bmdl.write('QZ......={:10.7f}\n'.format(atom.position[2]))
        else:
            bmdl.write('A........=     1.000 B.......=     1.000 C.......=     1.000\n')
            bmdl.write('ALFA.....=      90.0 BETA....=      90.0 GAMMA...=      90.0\n')
            bmdl.write('QX(1)....=       0.0 QY(1)...=       0.0 QZ(1)...=       0.0\n')

    def write_kstr(self, atoms):
        kstr = open('kstr.dat', 'w')
        kstr.write('KSTR      HP......=N                               22 Jan 08\n')
        kstr.write('JOBNAM...=kstr       MSGL.=  0 MODE...=B STORE..=Y HIGH...=Y\n')
        kstr.write('DIR001=./\n')
        kstr.write('DIR006=./\n')
        kstr.write('Slope matrices, fcc (spdf), (kappa*w)^2= 0.0\n')
        kstr.write('NL.....= 4 NLH...=11 NLW...= 9 NDER..= 6 ITRANS= 3 NPRN..= 0\n')
        kstr.write('(K*W)^2..=  0.000000 DMAX....=    1.7000 RWATS...=      0.10\n')

        kstr.write('NQ....={:3d} '.format(atoms.get_number_of_atoms()))
        kstr.write('LAT...={:2d} '.format(self.common_params['lat']))
        kstr.write('IPRIM.={:2d} '.format(self.common_params['iprim']))
        kstr.write('NGHBP.={:2d} '.format(self.common_params['nghbp']))
        kstr.write('NQR2..={:2d}\n'.format(self.common_params['nqr2']))

        if self.common_params['iprim']==0:
            cell = atoms.get_cell()

            a = np.linalg.norm(cell[0])
            b = np.linalg.norm(cell[1])
            c = np.linalg.norm(cell[2])

            kstr.write('A........={:10.3f} '.format(a / a))
            kstr.write('B.......={:10.3f} '.format(b / a))
            kstr.write('C.......={:10.3f}\n'.format(c / a))

            if self.common_params['lat'] == 2:  # fcc
                l = cell[1][0] * 2
            elif self.common_params['lat'] == 3:  # bcc
                l = cell[0][0] * 2
            else:
                l = cell[0][0]

            cell = atoms.get_cell().copy() / l
            kstr.write('BSX......={:10.7f} '.format(cell[0][0]))
            kstr.write('BSY.....={:10.7f} '.format(cell[0][1]))
            kstr.write('BSZ.....={:10.7f}\n'.format(cell[0][2]))
            kstr.write('BSX......={:10.7f} '.format(cell[1][0]))
            kstr.write('BSY.....={:10.7f} '.format(cell[1][1]))
            kstr.write('BSZ.....={:10.7f}\n'.format(cell[1][2]))
            kstr.write('BSX......={:10.7f} '.format(cell[2][0]))
            kstr.write('BSY.....={:10.7f} '.format(cell[2][1]))
            kstr.write('BSZ.....={:10.7f}\n'.format(cell[2][2]))

            for atom in atoms:
                kstr.write('QX.......={:10.7f} '.format(atom.position[0]))
                kstr.write('QY......={:10.7f} '.format(atom.position[1]))
                kstr.write('QZ......={:10.7f}\n'.format(atom.position[2]))
        else:
            kstr.write('A........=     1.000 B.......=     1.000 C.......=     1.000\n')
            kstr.write('ALFA.....=      90.0 BETA....=      90.0 GAMMA...=      90.0\n')
            kstr.write('QX(1)....=       0.0 QY(1)...=       0.0 QZ(1)...=       0.0\n')

        kstr.write('a/w(IQ)..= 0.70 0.70 0.70 0.70\n')

        kstr.write('LAMDA....={:10.4f} '.format(self.common_params['lamda']))
        kstr.write('AMAX....={:10.4f} '.format(self.common_params['amax']))
        kstr.write('BMAX....={:10.4f}\n'.format(self.common_params['bmax']))

    def write_shape(self, atoms):
        shape = open('shape.dat', 'w')
        shape.write('SHAPE     HP......=N                             22 Jan 08\n')
        shape.write('JOBNAM...=shape      MSGL.=  1\n')
        shape.write('FOR001=kstr.tfh\n')
        shape.write('DIR002=./\n')
        shape.write('DIR006=./\n')
        shape.write('Lmax..= 30 NSR..=129 NFI..= 11\n')
        shape.write('NPRN..=  0 IVEF.=  3\n')

    def write_kgrn(self, atoms):
        kgrn = open('kgrn.dat', 'w')
        kgrn.write('KGRN                                               13 Oct 12\n')
        kgrn.write('JOBNAM=kgrn\n')
        kgrn.write('STRT..=  A MSGL.=  0 EXPAN.= S FCD..=  Y FUNC..= SCA\n')
        kgrn.write('FOR001=kstr.tfh\n')
        kgrn.write('FOR001=kstr.tfh\n')
        kgrn.write('DIR002=./\n')
        kgrn.write('DIR003=./\n')
        kgrn.write('FOR004=bmdl.mdl\n')
        kgrn.write('DIR006=\n')
        kgrn.write('DIR009=./\n')
        kgrn.write('DIR010=./\n')
        kgrn.write('DIR011=/tmp/\n')
        kgrn.write('Self-consistent KKR calculation\n')

        kgrn.write('Band: 10 lines\n')

        kgrn.write('NITER.={:3d} '.format(self.kgrn_params['niter']))
        kgrn.write('NLIN.={:3d} '.format(self.kgrn_params['nlin']))
        kgrn.write('NPRN.={:3d} '.format(self.kgrn_params['nprn']))
        kgrn.write('NCPA.={:3d} '.format(self.kgrn_params['ncpa']))
        kgrn.write('NT...={:3d} '.format(self.kgrn_params['nt']))
        kgrn.write('MNTA.={:3d}\n'.format(self.kgrn_params['mnta']))

        kgrn.write('MODE..={:>3} '.format(self.kgrn_params['mode']))
        kgrn.write('FRC..={:>3} '.format(self.kgrn_params['frc']))
        kgrn.write('DOS..={:>3} '.format(self.kgrn_params['dos']))
        kgrn.write('OPS..={:>3} '.format(self.kgrn_params['ops']))
        kgrn.write('AFM..={:>3} '.format(self.kgrn_params['afm']))
        kgrn.write('CRT..={:>3}\n'.format(self.kgrn_params['crt']))

        kgrn.write('Lmaxh.={:3d} '.format(self.kgrn_params['lmaxh']))
        kgrn.write('Lmaxt={:3d} '.format(self.kgrn_params['lmaxt']))
        kgrn.write('NFI..={:3d} '.format(self.kgrn_params['nfi']))
        kgrn.write('FIXG.={:3d} '.format(self.kgrn_params['fixg']))
        kgrn.write('SHF..={:3d} '.format(self.kgrn_params['shf']))
        kgrn.write('SOFC.={:>3}\n'.format(self.kgrn_params['sofc']))

        kgrn.write('KMSH...={:>2} '.format(self.kgrn_params['kmsh']))
        kgrn.write('IBZ..={:3d} '.format(self.kgrn_params['ibz']))
        kgrn.write('NKX..={:3d} '.format(self.kgrn_params['nkx']))
        kgrn.write('NKY..={:3d} '.format(self.kgrn_params['nky']))
        kgrn.write('NKZ..={:3d} '.format(self.kgrn_params['nkz']))
        kgrn.write('FBZ..={:>3}\n'.format(self.kgrn_params['fbz']))

        kgrn.write('KMSH2..={:>2} '.format(self.kgrn_params['kmsh2']))
        kgrn.write('IBZ2.={:3d} '.format(self.kgrn_params['ibz2']))
        kgrn.write('NKX2.={:3d} '.format(self.kgrn_params['nkx2']))
        kgrn.write('NKY2.={:3d} '.format(self.kgrn_params['nky2']))
        kgrn.write('NKZ2.={:3d}\n'.format(self.kgrn_params['nkz2']))

        kgrn.write('ZMSH...={:>2} '.format(self.kgrn_params['zmsh']))
        kgrn.write('NZ1..={:3d} '.format(self.kgrn_params['nz1']))
        kgrn.write('NZ2..={:3d} '.format(self.kgrn_params['nz2']))
        kgrn.write('NZ3..={:3d} '.format(self.kgrn_params['nz3']))
        kgrn.write('NRES.={:3d} '.format(self.kgrn_params['nres']))
        kgrn.write('NZD.={:4d}\n'.format(self.kgrn_params['nzd']))

        kgrn.write('DEPTH..={:7.3f} '.format(self.kgrn_params['depth']))
        kgrn.write('IMAGZ.={:7.3f} '.format(self.kgrn_params['imagz']))
        kgrn.write('EPS...={:7.3f} '.format(self.kgrn_params['eps']))
        kgrn.write('ELIM..={:7.3f}\n'.format(self.kgrn_params['elim']))

        kgrn.write('AMIX...={:7.3f} '.format(self.kgrn_params['amix']))
        kgrn.write('EFMIX.={:7.3f} '.format(self.kgrn_params['efmix']))
        kgrn.write('VMTZ..={:7.3f} '.format(self.kgrn_params['vmtz']))
        kgrn.write('MMOM..={:7.3f}\n'.format(self.kgrn_params['mmom']))

        kgrn.write('TOLE...={:>7} '.format(self.kgrn_params['tole']))
        kgrn.write('TOLEF.={:>7} '.format(self.kgrn_params['tolef']))
        kgrn.write('TOLCPA={:>7} '.format(self.kgrn_params['tolcpa']))
        kgrn.write('TFERMI={:7.1f} (K)\n'.format(self.kgrn_params['tfermi']))

        sws = (atoms.get_volume()/len(atoms)*3.0/4.0/pi)**(1.0/3.0)/0.529177
        kgrn.write('SWS......={0:8.6f}   '.format(sws))
        kgrn.write('NSWS.={:3d} '.format(self.kgrn_params['nsws']))
        kgrn.write('DSWS..={:7.2f} '.format(self.kgrn_params['dsws']))
        kgrn.write('ALPCPA={:7.4f}\n'.format(self.kgrn_params['alpcpa']))

        kgrn.write('Setup: 2 + NQ*NS lines\n')
        kgrn.write('EFGS...=  0.000 HX....=  0.100 NX...=  5 NZ0..=  6 STMP..= N\n')
        kgrn.write('Symb   IQ IT ITA NZ  CONC   Sm(s)  S(ws) WS(wst) QTR SPLT Fix\n')

        i = 1
        for alloy in self.alloys:
            kgrn.write(alloy.symbol + '      1  1' + '{0:3d}'.format(i))
            nz = elements[alloy.symbol].split(' ')[2]
            kgrn.write('{0:4d}  '.format(int(nz)))
            kgrn.write('{:4.3f}'.format(alloy.conc))
            kgrn.write('  1.000  1.000  1.000  0.0 ')
            kgrn.write('{:4.1f}  N\n'.format(alloy.split))
            i = i+1

    #    for atom in atoms:
    #        kgrn.write(atom.symbol + '      1  1  1  ')
    #        nz = elements[atom.symbol].split(' ')[2]
    #        kgrn.write(nz)
    #        kgrn.write('  1.000  1.000  1.000  1.000  0.0  0.0  N\n')

    #    kgrn.write('Cu      1  1  1  26  1.000  1.000  1.000  1.000  0.0  0.0  N\n')
        kgrn.write('Atom:  4 lines + NT*NTA*6 lines\n')
        kgrn.write('IEX...=  4 NP..= 251 NES..= 15 NITER=100 IWAT.=  0 NPRNA=  0\n')
        kgrn.write('VMIX.....=  0.300000 RWAT....=  3.500000 RMAX....= 20.000000\n')
        kgrn.write('DX.......=  0.030000 DR1.....=  0.002000 TEST....=  1.00E-12\n')
        kgrn.write('TESTE....=  1.00E-12 TESTY...=  1.00E-12 TESTV...=  1.00E-12\n')

        for alloy in self.alloys:
            kgrn.write(alloy.symbol + '\n')
            kgrn.write(elements[alloy.symbol])
    #    for atom in atoms:
    #        kgrn.write(atom.symbol + '\n')
    #        kgrn.write(elements[atom.symbol])

    def write_kfcd(self, atoms):
        kfcd = open('kfcd.dat', 'w')
        kfcd.write('KFCD      MSGL..=  1                               22 Jan 08\n')
        kfcd.write('JOBNAM...=kfcd\n')
        kfcd.write('STRNAM...=common\n')
        kfcd.write('DIR001=./\n')
        kfcd.write('DIR002=./\n')
        kfcd.write('DIR003=./\n')
        kfcd.write('DIR004=./\n')
        kfcd.write('DIR006=./\n')
        kfcd.write('Lmaxs.= 30 NTH..= 41 NFI..= 81 FPOT..= N\n')
        kfcd.write('OVCOR.=  Y UBG..=  N NPRN.=  0\n')

    def set_alloys(self, alloys):
        self.alloys = alloys

    def read_energy(self, all=None):
        file = open('kfcd.prn','r')
        lines = file.readlines()
        file.close()

        for line in lines:
            if line.rfind('TOT-PBE') > -1:
                self.energy_pbe = float(line.split('(Ry)')[0].split('TOT-PBE')[1].strip())*13.6058
            #    self.energy_pbe = float(line.split(' ')[8].strip())*13.6058
            #    self.energy_pbe = self.energy_pbe*13.6058 # Ry -> eV conversion




