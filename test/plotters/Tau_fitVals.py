#!/usr/bin/env python
'''
Values for L1pT vs recoPt from fitTau.py
'''
def getVals(thresh,FixPlat,VarBins,AbsEff):
 if thresh ==0.85:
  if not FixPlat and VarBins:
   mUIso=0.752941
   bUIso=28.8971
   mU=1.07279
   bU=21.1397
   mC=0.692308 
   bC=2.92308
   pUIso=0.702131383694
   pU=0.977643518384
   pC=0.691563180805

 if thresh == 0.9:
  if not FixPlat and VarBins:
   mUIso=1.25368
   bUIso=19.8088
   mU=1.14559
   bU=22.4044
   mC=0.736364 
   bC=2.85455
   pUIso=0.73991148688
   pU=0.974768978201
   pC=0.626672832606

 if thresh==0.95:
  if not FixPlat and VarBins:
   mUIso=0.695588
   bUIso=39.7794
   mU=1.28382
   bU=23.4412
   mC=0.745455 
   bC=4.25455
   pUIso=0.704664365692
   pU=0.971301127603
   pC=0.691563180805

 if thresh==0.4 and AbsEff:
  if not FixPlat:
   mUIso=0.738235
   bUIso=18.6618
   mU=0.697059
   bU=15.8529
   mC=0.739011 
   bC=-1.6978
   pUIso=0.71131612275
   pU=0.969063841137
   pC=0.691563180805
  
 if thresh==0.5 and AbsEff:
  if not FixPlat:
   mUIso=0.752941
   bUIso=23.0221
   mU=0.773529
   bU=16.6765
   mC=0.81044
   bC=-2.1044
   pUIso=0.71131612275
   pU=0.969063841137
   pC=0.691563180805
  
 if thresh==0.6 and AbsEff:
  if not FixPlat:
   mUIso=0.738971
   bUIso=29.6985
   mU=0.85
   bU=17.5
   mC=1.36667 
   bC=-15.6889
   pUIso=0.71131612275
   pU=0.969063841137
   pC=0.691563180805
  
 return mUIso,bUIso,pUIso,mU,bU,pU,mC,bC,pC
