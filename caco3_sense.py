# -*- coding: utf-8 -*-
import sys
import os
crntdir = os.getcwd()
if not crntdir in sys.path:sys.path.append(crntdir)
from caco3 import caco3_main
import numpy as np

ccflxi      = 12e-6         # CaCO3 rain flux in mol cm-2 yr-1  
om2cc       = 0.7           # OM/CaCO3 rain ratio 
dep         = 3.5           # water depth in km  
dt          = 1e8           # time step in yr  
fl          = 'sense'       # simulation name 
biot        = 'fickian'     # bioturbation style  
oxonly      = False         # Oxic only for OM degradation?  
runmode     = 'sense'       # simulation mode
co2chem     = 'co2'         # co2 chemistry
sparse      = False         # use sparse matrix solver?
showiter    = False         # show every iteration step?
rrlist      = np.array([0.0,0.5,0.6666,1.0,1.5])
nz          = 25
for i in range(nz):
    for j in range(10):
        for k in range(5):
            for l in range(2):
                ccflxi = int((j+1)*6)*1e-6
                om2cc = rrlist[k]
                dep =(i+1.)*6.0/(1.0*nz)
                if l==0:oxonly=True  
                if l==1:oxonly=False
                caco3_main(ccflxi,om2cc,dep,dt,fl,biot
                           ,oxonly,runmode,co2chem,sparse,showiter)
