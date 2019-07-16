# -*- coding: utf-8 -*-
"""
signal tracking diagenesis model
"""
import os
if os.path.exists('./caco3mod.dll'):import caco3mod

ccflxi      = 12e-6                 # caco3 rain flux in mol cm-2 yr-1
om2cc       = 0.7                   # om/caco3 rain ratio 
dep         = 3.5                   # water depth in km 
dtinput     = 1e2                   # time step in yr 
runname     = 'test'                # file name 
oxonly      = False                 # oxic-only OM degradation model or not 
biotmode    = 'fickian'             # biomixing styles, choose from 'fickian','nobio','turbo2','labs'

caco3mod.caco3_python(ccflxi,om2cc,dep,dtinput,runname,oxonly,biotmode) 
