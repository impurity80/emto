#!/bin/bash

cd bmdl
time bmdl < fcc.dat

cd ..
cd kstr
time kstr < fcc.dat

cd ..
cd shape
time shape < fcc.dat

cd ..
cd kgrn
time kgrn_cpa < cu1.dat

cd ..
cd kfcd
time kfcd_cpa < cu1.dat

cd ..
