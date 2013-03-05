#!/usr/bin/env python
'''
Values for L1pT vs recoPt from fitTau.py
'''
def getVals(thresh,FixPlat,VarBins):
 if thresh ==0.85:
  if not FixPlat and VarBins:
   mUIso=0.752941
   bUIso=28.8971
   mU=1.07279
   bU=21.1397
   mC=1.64265
   bC=6.63235

 if thresh == 0.9:
  if FixPlat and not VarBins:
   mUIso=0.691
   bUIso=39.0
   mU=1.25
   bU=20.5
   mC=1.48
   bC=13.4
  if not FixPlat and not VarBins:
   mUIso=0.846
   bUIso=34.2
   mU=1.14
   bU=22.9
   mC=1.81
   bC=6.98
  if not FixPlat and VarBins:
#   mUIso=0.688971
#   bUIso=35.8235
#   mU=1.20221
#   bU=20.3603
#   mC=1.74706
#   bC=8.47794
   mUIso=1.29412
   bUIso=18.7059
   mU=1.14559
   bU=22.4044
   mC=1.74706
   bC=8.47794

 if thresh==0.95:
  if FixPlat and not VarBins:
   mUIso=0.692
   bUIso=42.6
   mU=1.51
   bU=19.1
   mC=1.38
   bC=22.7
  if not FixPlat and not VarBins:
   mUIso=0.792
   bUIso=37.7
   mU=1.25
   bU=24.4
   mC=2.00
   bC=2.72
  if not FixPlat and VarBins:
   mUIso=0.721324
   bUIso=39.9412
   mU=1.35368
   bU=21.1838
   mC=2.00147
   bC=5.07353

 return mUIso,bUIso,mU,bU,mC,bC
 
# From Jim
#mUIso = 0.8
#bUIso = 46
#mU = 1.2
#bU = 4.5
#mC = 1.07
#bC = 11
