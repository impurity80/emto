from __future__ import print_function

import os
import sys
import re
from ase.calculators.general import Calculator
from os.path import join, isfile, islink

import numpy as np

import ase
import ase.io
from ase.utils import devnull

from ase.calculators.singlepoint import SinglePointCalculator

class EMTO(Calculator):
    name = 'EMTO'

    def __init__(self, restart=None, **kwargs):
        self.bmdl_params = {}
        self.kstr_params = {}
        self.shape_params = {}
        self.kgrn_params = {}
        self.kfcd_params = {}

    def set(self, **kargs):
        a = 0

    def update(self, atoms):
        b = 0

    def initialize(self, atoms):
        self.natoms = len(atoms)

    def calculation_required(selfself, atoms, quantities):
        return True

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)

    def update(self, atoms):
        self.calculate(atoms)

    def calculate(self, atoms):
        self.initialize(atoms)
        self.write_bmdl(atoms)
        os.system('time -f " Elapsed time = %E" bmdl < bmdl.dat')
        self.write_kstr(atoms)
        os.system('time -f " Elapsed time = %E" kstr < kstr.dat')
        self.write_shape(atoms)
        os.system('time -f " Elapsed time = %E" shape < shape.dat')
        self.write_kgrn(atoms)
        os.system('time -f " Elapsed time = %E" kgrn_cpa < kgrn.dat')

        os.system('mv bmdl.mdl common.mdl')
        os.system('mv shape.shp common.shp')
        os.system('mv kstr.tfh common.tfh')
        os.system('mv kgrn.chd kfcd.chd')

        self.write_kfcd(atoms)
        os.system('time -f " Elapsed time = %E" kfcd_cpa < kfcd.dat')

    def write_bmdl(self, atoms):
        bmdl = open('bmdl.dat', 'w')
        bmdl.write('BMDL      HP......=N                               22 Jan 08\n')
        bmdl.write('JOBNAM...=bmdl       MSGL.=  1 NPRN.=  0\n')
        bmdl.write('DIR001=./\n')
        bmdl.write('DIR006=\n')
        bmdl.write('Madelung potential for fcc bulk\n')
        bmdl.write('NL.....= 7\n')
        bmdl.write('LAMDA....=      2.50 AMAX....=      4.50 BMAX....=      4.50\n')
        bmdl.write('NQ....=  1 LAT...= 2 IPRIM.= 1 NQR2..= 0\n')
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
        kstr.write('NQ3...=  1 LAT...= 2 IPRIM.= 0 NGHBP.=13 NQR2..= 0\n')
        kstr.write('A........= 1.0000000 B.......= 1.0000000 C.......= 1.0000000\n')
        kstr.write('BSX......= 0.5000000 BSY.....= 0.5000000 BSZ.....= 0.0000000\n')
        kstr.write('BSX......= 0.0000000 BSY.....= 0.5000000 BSZ.....= 0.5000000\n')
        kstr.write('BSX......= 0.5000000 BSY.....= 0.0000000 BSZ.....= 0.5000000\n')
        kstr.write('QX(IQ)...= 0.0000000 QY......= 0.0000000 QZ......= 0.0000000\n')
        kstr.write('a/w(IQ)..= 0.70 0.70 0.70 0.70\n')
        kstr.write('LAMDA....=    2.5000 AMAX....=    4.5000 BMAX....=    4.5000\n')

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
        kgrn.write('Self-consistent KKR calculation for fcc Cu\n')
        kgrn.write('Band: 10 lines\n')
        kgrn.write('NITER.= 50 NLIN.= 31 NPRN.=  0 NCPA.=  7 NT...=  1 MNTA.=  1\n')
        kgrn.write('MODE..= 3D FRC..=  N DOS..=  N OPS..=  N AFM..=  P CRT..=  M\n')
        kgrn.write('Lmaxh.=  8 Lmaxt=  4 NFI..= 31 FIXG.=  2 SHF..=  0 SOFC.=  N\n')
        kgrn.write('KMSH...= G IBZ..=  2 NKX..=  0 NKY..= 13 NKZ..=  0 FBZ..=  N\n')
        kgrn.write('KMSH2..= G IBZ2.=  1 NKX2.=  4 NKY2.=  0 NKZ2.= 51\n')
        kgrn.write('ZMSH...= C NZ1..= 16 NZ2..=  8 NZ3..=  8 NRES.=  4 NZD.=1500\n')
        kgrn.write('DEPTH..=  1.000 IMAGZ.=  0.020 EPS...=  0.200 ELIM..= -1.000\n')
        kgrn.write('AMIX...=  0.100 EFMIX.=  1.000 VMTZ..=  0.000 MMOM..=  0.000\n')
        kgrn.write('OLE...= 1.d-07 TOLEF.= 1.d-07 TOLCPA= 1.d-06 TFERMI=  500.0 (K)\n')
        kgrn.write('SWS......=2.686842   NSWS.=  1 DSWS..=   0.05 ALPCPA= 0.6020\n')
        kgrn.write('Setup: 2 + NQ*NS lines\n')
        kgrn.write('EFGS...=  0.000 HX....=  0.100 NX...=  5 NZ0..=  6 STMP..= N\n')
        kgrn.write('Symb   IQ IT ITA NZ  CONC   Sm(s)  S(ws) WS(wst) QTR SPLT Fix\n')
        kgrn.write('Cu      1  1  1  29  1.000  1.000  1.000  1.000  0.0  0.0  N\n')
        kgrn.write('Atom:  4 lines + NT*NTA*6 lines\n')
        kgrn.write('IEX...=  4 NP..= 251 NES..= 15 NITER=100 IWAT.=  0 NPRNA=  0\n')
        kgrn.write('VMIX.....=  0.300000 RWAT....=  3.500000 RMAX....= 20.000000\n')
        kgrn.write('DX.......=  0.030000 DR1.....=  0.002000 TEST....=  1.00E-12\n')
        kgrn.write('TESTE....=  1.00E-12 TESTY...=  1.00E-12 TESTV...=  1.00E-12\n')
        kgrn.write('Cu\n')
        kgrn.write('Iz=  29 Norb= 10 Ion=  0 Config= 3d10_4s1\n')
        kgrn.write('n      1  2  2  2  3  3  3  3  3  4\n')
        kgrn.write('Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n')
        kgrn.write('Occup  2  2  2  4  2  2  4  4  6  1\n')
        kgrn.write('Valen  0  0  0  0  0  0  0  1  1  1\n')

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

