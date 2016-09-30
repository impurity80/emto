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

    def write_bmdl(self, atoms):
        bmdl = open('bmdl.dat', 'w')
        bmdl.write('BMDL      HP......=N                               22 Jan 08\n')
        bmdl.write('JOBNAM...=fcc        MSGL.=  1 NPRN.=  0\n')
        bmdl.write('DIR001=./\n')
        bmdl.write('DIR006=\n')
        bmdl.write('Madelung potential for fcc bulk\n')
        bmdl.write('NL.....= 7\n')
        bmdl.write('LAMDA....=      2.50 AMAX....=      4.50 BMAX....=      4.50\n')
        bmdl.write('NQ....=  1 LAT...= 2 IPRIM.= 1 NQR2..= 0\n')
        bmdl.write('A........=     1.000 B.......=     1.000 C.......=     1.000\n')
        bmdl.write('ALFA.....=      90.0 BETA....=      90.0 GAMMA...=      90.0\n')
        bmdl.write('QX(1)....=       0.0 QY(1)...=       0.0 QZ(1)...=       0.0\n')


