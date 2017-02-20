
# *******************************************************
#  Copyright (C) Korea Institute of Materials Science - All Rights Reserved
#  Unauthorized copying of this file, via any medium is strictly prohibited
#  Proprietary and confidential
#  Written by Jae Hoon Jang <jhjang@kims.re.kr>, January 2017
#  Version : 0.1
# *******************************************************/

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

def save( filename, arg ):
    f = open(filename, 'a+t')
    f.write('{0} \n'.format(arg))
    f.close()

class Alloy():
    def __init__(self, id, symbol, conc, split):
        self.id = id
        self.symbol = symbol
        self.conc = conc
        self.split = split

#element_keys = [
#    'Al', 'Si', 'P', 'S', # 2d metal
#    'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn', # 3d transition metal
#    'Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd', # 4d transition metal
#
#]

elements = {}
#for key in element_keys:
#    elements[key] = None

elements['B'] =  'Iz=   5 Norb=  3 Ion=  0 Config= 2s2_2p1\n' \
                  'n      1  2  2\n' \
                  'Kappa -1 -1  1\n' \
                  'Occup  2  2  1\n' \
                  'Valen  0  1  1\n'
elements['C'] =  'Iz=   6 Norb=  3 Ion=  0 Config= 2s2_2p2\n' \
                  'n      1  2  2\n' \
                  'Kappa -1 -1  1\n' \
                  'Occup  2  2  2\n' \
                  'Valen  0  1  1\n'
elements['N'] =  'Iz=   7 Norb=  4 Ion=  0 Config= 2s2_2p3\n' \
                  'n      1  2  2  2\n' \
                  'Kappa -1 -1  1 -2\n' \
                  'Occup  2  2  2  1\n' \
                  'Valen  0  1  1  1\n'
elements['O'] =  'Iz=   8 Norb=  4 Ion=  0 Config= 2s2_2p4\n' \
                  'n      1  2  2  2\n' \
                  'Kappa -1 -1  1 -2\n' \
                  'Occup  2  2  2  2\n' \
                  'Valen  0  1  1  1\n'
elements['Mg'] =  'Iz=  12 Norb=  5 Ion=  0 Config= 3s2\n' \
                  'n      1  2  2  2  3\n' \
                  'Kappa -1 -1  1 -2 -1\n' \
                  'Occup  2  2  2  4  2\n' \
                  'Valen  0  0  0  0  1\n'
elements['Al'] =  'Iz=  13 Norb=  6 Ion=  0 Config= 3s2_3p1\n' \
                  'n      1  2  2  2  3  3\n' \
                  'Kappa -1 -1  1 -2 -1  1\n' \
                  'Occup  2  2  2  4  2  1\n' \
                  'Valen  0  0  0  0  1  1\n'
elements['Si'] =  'Iz=  14 Norb=  6 Ion=  0 Config= 3s2_3p2\n' \
                  'n      1  2  2  2  3  3\n' \
                  'Kappa -1 -1  1 -2 -1  1\n' \
                  'Occup  2  2  2  4  2  2\n' \
                  'Valen  0  0  0  0  1  1\n'
elements['P'] =  'Iz=  15 Norb=  6 Ion=  0 Config= 3s2_3p3\n' \
                  'n      1  2  2  2  3  3\n' \
                  'Kappa -1 -1  1 -2 -1 -2\n' \
                  'Occup  2  2  2  4  2  3\n' \
                  'Valen  0  0  0  0  1  1\n'
elements['S'] =  'Iz=  16 Norb=  6 Ion=  0 Config= 3s2_3p4\n' \
                  'n      1  2  2  2  3  3\n' \
                  'Kappa -1 -1  1 -2 -1 -2\n' \
                  'Occup  2  2  2  4  2  4\n' \
                  'Valen  0  0  0  0  1  1\n'
elements['Sc'] =  'Iz=  21 Norb=  9 Ion=  0 Config= 3d1_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  1  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1\n'
elements['Ti'] =  'Iz=  22 Norb=  9 Ion=  0 Config= 3d2_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  2  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1\n'
elements['V'] =   'Iz=  23 Norb=  9 Ion=  0 Config= 3d3_4s2\n' \
                  'n      1  2  2  2  3  3  3  3  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  3  2\n' \
                  'Valen  0  0  0  0  0  0  0  1  1\n'
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
elements['Y'] =  'Iz=  39 Norb= 14 Ion=  0 Config= 4d1_5s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  1  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n'
elements['Zr'] =  'Iz=  40 Norb= 14 Ion=  0 Config= 4d2_5s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  2  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n'
elements['Nb'] =  'Iz=  41 Norb= 14 Ion=  0 Config= 4d4_5s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  1\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n'
elements['Mo'] =  'Iz=  42 Norb= 15 Ion=  0 Config= 4d5_5s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  1  1\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['Tc'] =  'Iz=  43 Norb= 15 Ion=  0 Config= 4d6_5s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  2  1\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['Ru'] =  'Iz=  44 Norb= 15 Ion=  0 Config= 4d7_5s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  3  1\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['Rh'] =  'Iz=  45 Norb= 15 Ion=  0 Config= 4d8_5s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  4  1\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['Pd'] =  'Iz=  46 Norb= 14 Ion=  0 Config= 4d10_5s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n'
elements['Ag'] =  'Iz=  47 Norb= 15 Ion=  0 Config= 4d10_5s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  1\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['Cd'] =  'Iz=  48 Norb= 15 Ion=  0 Config= 4d10_5s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['In'] =  'Iz=  49 Norb= 16 Ion=  0 Config= 4d10_5s25p1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  1\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n'
elements['Sn'] =  'Iz=  50 Norb= 16 Ion=  0 Config= 4d10_5s25p2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n'
elements['La'] =  'Iz=  57 Norb= 19 Ion=  0 Config= 5d1_6s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  4  1  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n'
elements['Ce'] =  'Iz=  58 Norb= 19 Ion=  0 Config= 4f2_6s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  4 -1  1 -2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  2  4  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1\n'
elements['Pr'] =  'Iz=  59 Norb= 19 Ion=  0 Config= 4f3_6s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -1  1 -2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  3  2  2  4  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1\n'
elements['Hf'] =  'Iz=  72 Norb= 21 Ion=  0 Config= 5d2_6s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  2  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n'
elements['Ta'] =  'Iz=  73 Norb= 21 Ion=  0 Config= 5d3_6s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  3  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n'
elements['W'] =  'Iz=  74 Norb= 21 Ion=  0 Config= 5d4_6s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n'
elements['Re'] =  'Iz=  75 Norb= 22 Ion=  0 Config= 5d5_6s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  1  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['Os'] =  'Iz=  76 Norb= 22 Ion=  0 Config= 5d6_6s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  2  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['Ir'] =  'Iz=  77 Norb= 21 Ion=  0 Config= 5d9\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  5\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n'
elements['Pt'] =  'Iz=  78 Norb= 22 Ion=  0 Config= 5d9_6s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  5  1\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['Au'] =  'Iz=  79 Norb= 22 Ion=  0 Config= 5d10_6s1\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  1\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'
elements['Hg'] =  'Iz=  80 Norb= 22 Ion=  0 Config= 5d10_6s2\n' \
                  'n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n' \
                  'Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n' \
                  'Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2\n' \
                  'Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n'


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



kstr_keys = [
    'nl',
    'nlh',
    'nlw',
    'nder',
    'nprn',
    'dmax',
    'rwats'
    'kw2'
]

kgrn_keys = [
    'strt',
    'msgl',
    'expan',
    'fcd',
    'func',
    'niter',
    'nlin',
    'nprn',
    'ncpa',
#    'nt',
#    'mnta',
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
    'kpts'
#    'nkx',
#    'nky',
#    'nkz',
    'fbz',
    'kmsh2',
    'ibz2',
    'kpts2'
#    'nkx2',
#    'nky2',
#    'nkz2',
    'zmsh',
    'zpts',
#    'nz1',
#    'nz2',
#    'nz3',
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

        self.energy_lda = 0
        self.energy_pbe = 0
        self.energy_p07 = 0
        self.energy_am5 = 0
        self.energy_lag = 0

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

        self.kstr_params['nl'] = 4
        self.kstr_params['nlh'] = 11
        self.kstr_params['nlw'] = 9
        self.kstr_params['nder'] = 6
        self.kstr_params['itrans'] = 3
        self.kstr_params['nprn'] = 0
        self.kstr_params['dmax'] = 1.70
        self.kstr_params['rwats'] = 0.10
        self.kstr_params['kw2'] = 0.0

        self.kgrn_params['strt'] = 'A'
        self.kgrn_params['msgl'] = 0
        self.kgrn_params['expan'] = 'S'
        self.kgrn_params['fcd'] = 'Y'
        self.kgrn_params['func'] = 'SCA'

        self.kgrn_params['niter'] = 100
        self.kgrn_params['nlin'] = 31
        self.kgrn_params['nprn'] = 0
        self.kgrn_params['ncpa'] = 20
    #    self.kgrn_params['nt'] = 1
    #    self.kgrn_params['mnta'] = 4
        self.kgrn_params['mode'] = '3D'
        self.kgrn_params['frc'] = 'Y'
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
    #    self.kgrn_params['ibz'] = 2
    #    self.kgrn_params['nkx'] = 0
    #    self.kgrn_params['nky'] = 13
    #    self.kgrn_params['nkz'] = 0
        self.kgrn_params['kpts'] = [1,1,1]
        self.kgrn_params['fbz'] = 'N'
        self.kgrn_params['kmsh2'] = 'G'
        self.kgrn_params['ibz2'] = 1
    #    self.kgrn_params['nkx2'] = 4
    #    self.kgrn_params['nky2'] = 0
    #    self.kgrn_params['nkz2'] = 51
        self.kgrn_params['kpts2'] = [1,1,1]
        self.kgrn_params['zmsh'] = 'C'

        self.kgrn_params['zpts'] = [16,16,16]
    #    self.kgrn_params['nz1'] = 16
    #    self.kgrn_params['nz2'] = 16
    #    self.kgrn_params['nz3'] = 8
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
        self.kgrn_params['tfermi'] = 500.0
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
            if key in self.kstr_params:
                self.kstr_params[key] = kwargs[key]
            if key in self.kgrn_params:
                self.kgrn_params[key] = kwargs[key]

    def set_kgrn(self, **kwargs):
        for key in kwargs:
            self.kgrn_params[key] = kwargs[key]

    def initialize(self, atoms):
        self.natoms = len(atoms)

        for atom in atoms:
            if atom.tag==0:
                self.alloys.append(Alloy(atom.tag, atom.symbol, 1.0, atom.magmom))

    #    self.alloys = []
    #    for atom in atoms:
    #        self.alloys.append(Alloy(atom.index, atom.symbol, 1.0, 1.0))

    def calculation_required(self, atoms, quantities):
        return True

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        return self.energy_pbe

    def get_potential_energies():
        return [self.energy_lda, self.energy_pbe, self.energy_p07, self.energy_am5, self.energy_lag]

    def update(self, atoms):
        self.calculate(atoms)
#        self.calculate(atoms)

    def calculate2(self, atoms):
        os.chdir(self.common_params['dir'])
        self.set_results(atoms)

    def calculate(self, atoms):
        cur_dir = os.getcwd()
        os.system('rm -r {0}'.format(self.common_params['dir']))
        os.system('mkdir -p {0}'.format(self.common_params['dir']))
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

        os.chdir(cur_dir)

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

        cell = atoms.get_cell()

        a = np.linalg.norm(cell[0])
        b = np.linalg.norm(cell[1])
        c = np.linalg.norm(cell[2])

        bmdl.write('A........={:10.7f} '.format(a / a))
        bmdl.write('B.......={:10.7f} '.format(b / a))
        bmdl.write('C.......={:10.7f}\n'.format(c / a))

        if self.common_params['iprim']==0:

            if len(atoms) == 1 : #primitive
                l = cell[0][1]*2
            else:
                l = cell[0][0]

       #     if self.common_params['lat']==2: # fcc
       #         l = cell[1][0]*2
       #     elif self.common_params['lat']==3: # bcc
       #         l = cell[0][0]*2
       #     else:
       #         l = cell[0][0]

            cell = atoms.get_cell().copy()/l
            bmdl.write('BSX......={:10.7f} '.format(cell[0][0]))
            bmdl.write('BSY.....={:10.7f} '.format(cell[1][0]))
            bmdl.write('BSZ.....={:10.7f}\n'.format(cell[2][0]))
            bmdl.write('BSX......={:10.7f} '.format(cell[0][1]))
            bmdl.write('BSY.....={:10.7f} '.format(cell[1][1]))
            bmdl.write('BSZ.....={:10.7f}\n'.format(cell[2][1]))
            bmdl.write('BSX......={:10.7f} '.format(cell[0][2]))
            bmdl.write('BSY.....={:10.7f} '.format(cell[1][2]))
            bmdl.write('BSZ.....={:10.7f}\n'.format(cell[2][2]))

            for atom in atoms:
                bmdl.write('QX.......={:10.7f} '.format(atom.position[0]/a))
                bmdl.write('QY......={:10.7f} '.format(atom.position[1]/a))
                bmdl.write('QZ......={:10.7f}\n'.format(atom.position[2]/a))
        else:
            bmdl.write('ALPHA....=      90.0 BETA....=      90.0 GAMMA...=      90.0\n')
            bmdl.write('QX(1)....=       0.0 QY(1)...=       0.0 QZ(1)...=       0.0\n')

    def write_kstr(self, atoms):
        kstr = open('kstr.dat', 'w')
        kstr.write('KSTR      HP......=N                               22 Jan 08\n')
        kstr.write('JOBNAM...=kstr       MSGL.=  0 MODE...=B STORE..=Y HIGH...=Y\n')
        kstr.write('DIR001=./\n')
        kstr.write('DIR006=./\n')
        kstr.write('Slope matrices, fcc (spdf), (kappa*w)^2= 0.0\n')

        kstr.write('NL.....={:2} '.format(self.kstr_params['nl']))
        kstr.write('NLH...={:2} '.format(self.kstr_params['nlh']))
        kstr.write('NLW...={:2} '.format(self.kstr_params['nlw']))
        kstr.write('NDER..={:2} '.format(self.kstr_params['nder']))
        kstr.write('ITRANS={:2} '.format(self.kstr_params['itrans']))
        kstr.write('NPRN..={:2}\n'.format(self.kstr_params['nprn']))

        kstr.write('(K*W)^2..={:10.7f} '.format(self.kstr_params['kw2']))
        kstr.write('DMAX....={:10.7f} '.format(self.kstr_params['dmax']))
        kstr.write('RWATS...={:10.7f}\n'.format(self.kstr_params['rwats']))

        kstr.write('NQ....={:3d} '.format(atoms.get_number_of_atoms()))
        kstr.write('LAT...={:2d} '.format(self.common_params['lat']))
        kstr.write('IPRIM.={:2d} '.format(self.common_params['iprim']))
        kstr.write('NGHBP.={:2d} '.format(self.common_params['nghbp']))
        kstr.write('NQR2..={:2d}\n'.format(self.common_params['nqr2']))

        cell = atoms.get_cell()

        a = np.linalg.norm(cell[0])
        b = np.linalg.norm(cell[1])
        c = np.linalg.norm(cell[2])

        kstr.write('A........={:10.7f} '.format(a / a))
        kstr.write('B.......={:10.7f} '.format(b / a))
        kstr.write('C.......={:10.7f}\n'.format(c / a))

        if self.common_params['iprim']==0:

            if len(atoms) == 1:  # primitive
                l = cell[0][1] * 2
            else:
                l = cell[0][0]

                #     if self.common_params['lat']==2: # fcc
                #         l = cell[1][0]*2
                #     elif self.common_params['lat']==3: # bcc
                #         l = cell[0][0]*2
                #     else:
                #         l = cell[0][0]

            cell = atoms.get_cell().copy() / l
            kstr.write('BSX......={:10.7f} '.format(cell[0][0]))
            kstr.write('BSY.....={:10.7f} '.format(cell[1][0]))
            kstr.write('BSZ.....={:10.7f}\n'.format(cell[2][0]))
            kstr.write('BSX......={:10.7f} '.format(cell[0][1]))
            kstr.write('BSY.....={:10.7f} '.format(cell[1][1]))
            kstr.write('BSZ.....={:10.7f}\n'.format(cell[2][1]))
            kstr.write('BSX......={:10.7f} '.format(cell[0][2]))
            kstr.write('BSY.....={:10.7f} '.format(cell[1][2]))
            kstr.write('BSZ.....={:10.7f}\n'.format(cell[2][2]))

            for atom in atoms:
                kstr.write('QX.......={:10.7f} '.format(atom.position[0]/a))
                kstr.write('QY......={:10.7f} '.format(atom.position[1]/a))
                kstr.write('QZ......={:10.7f}\n'.format(atom.position[2]/a))
        else:
        #    kstr.write('A........=     1.000 B.......=     1.000 C.......=     1.000\n')
            kstr.write('ALPHA....=      90.0 BETA....=      90.0 GAMMA...=      90.0\n')
            kstr.write('QX(1)....=       0.0 QY(1)...=       0.0 QZ(1)...=       0.0\n')

        for atom in atoms:
            if atom.symbol == 'C':
                kstr.write('a/w(IQ)..= 0.60 0.60 0.60 0.60\n')
            else :
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

        kgrn.write('STRT..={:>3} '.format(self.kgrn_params['strt']))
        kgrn.write('MSGL.={:3d} '.format(self.kgrn_params['msgl']))
        kgrn.write('EXPAN.={:>2} '.format(self.kgrn_params['expan']))
        kgrn.write('FCD..={:>3} '.format(self.kgrn_params['fcd']))
        kgrn.write('FUNC..={:>4}\n'.format(self.kgrn_params['func']))

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
        kgrn.write('NT...={:3d} '.format(atoms.get_number_of_atoms()))

        kgrn.write('MNTA.={:3d}\n'.format(len([alloy.conc for alloy in self.alloys if alloy.id == 1])))

#         kgrn.write('MNTA.={:3d}\n'.format(len(self.alloys)))

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
        kgrn.write('IBZ..={:3d} '.format(self.common_params['lat']))
        kgrn.write('NKX..={:3d} '.format(self.kgrn_params['kpts'][0]))
        kgrn.write('NKY..={:3d} '.format(self.kgrn_params['kpts'][1]))
        kgrn.write('NKZ..={:3d} '.format(self.kgrn_params['kpts'][2]))
        kgrn.write('FBZ..={:>3}\n'.format(self.kgrn_params['fbz']))

        kgrn.write('KMSH2..={:>2} '.format(self.kgrn_params['kmsh2']))
        kgrn.write('IBZ2.={:3d} '.format(self.kgrn_params['ibz2']))
        kgrn.write('NKX2.={:3d} '.format(self.kgrn_params['kpts2'][0]))
        kgrn.write('NKY2.={:3d} '.format(self.kgrn_params['kpts2'][1]))
        kgrn.write('NKZ2.={:3d}\n'.format(self.kgrn_params['kpts2'][2]))

        kgrn.write('ZMSH...={:>2} '.format(self.kgrn_params['zmsh']))
        kgrn.write('NZ1..={:3d} '.format(self.kgrn_params['zpts'][0]))
        kgrn.write('NZ2..={:3d} '.format(self.kgrn_params['zpts'][1]))
        kgrn.write('NZ3..={:3d} '.format(self.kgrn_params['zpts'][2]))
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

        for atom in atoms:
            i = 1
            for alloy in self.alloys:
                if alloy.id == atom.tag:
                    kgrn.write('{0:<2}'.format(alloy.symbol))
                    kgrn.write('{0:7d}{0:3d}'.format(atom.index+1))
                    kgrn.write('{0:3d}'.format(i))
                    nz = elements[alloy.symbol].split(' ')[2]
                    kgrn.write('{0:4d}  '.format(int(nz)))
                    kgrn.write('{:4.3f}'.format(alloy.conc))
                    kgrn.write('  1.000  1.000  1.000  0.0 ')
                    kgrn.write('{:4.1f}  N\n'.format(alloy.split))
                    i = i+1

             #   print(alloy.id, alloy.symbol, alloy.conc)

#        i = 1
#        for atom in atoms:

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

        for atom in atoms:
            for alloy in self.alloys:
                if alloy.id == atom.tag:
                    kgrn.write(alloy.symbol + '\n')
                    kgrn.write(elements[alloy.symbol])

    #        for alloy in self.alloys:
    #            kgrn.write(alloy.symbol + '\n')
    #            kgrn.write(elements[alloy.symbol])
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
            if line.rfind('TOT-LDA') > -1:
                self.energy_lda = float(line.split('(Ry)')[0].split('TOT-LDA')[1].strip())*13.6058
            if line.rfind('TOT-P07') > -1:
                self.energy_p07 = float(line.split('(Ry)')[0].split('TOT-P07')[1].strip())*13.6058
            if line.rfind('TOT-AM5') > -1:
                self.energy_am5 = float(line.split('(Ry)')[0].split('TOT-AM5')[1].strip())*13.6058
            if line.rfind('TOT-LAG') > -1:
                self.energy_lag = float(line.split('(Ry)')[0].split('TOT-LAG')[1].strip()) * 13.6058
            #    self.energy_pbe = float(line.split(' ')[8].strip())*13.6058
            #    self.energy_pbe = self.energy_pbe*13.6058 # Ry -> eV conversion




