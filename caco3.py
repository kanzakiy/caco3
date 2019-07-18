# -*- coding: utf-8 -*-
"""
signal tracking diagenesis model
"""
import numpy as np
import os
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
#import scipy.linalg as la
import scipy.linalg.lapack as lapack
if os.path.exists('./mocsy.dll'):import mocsy

np.set_printoptions(formatter={'float': '{:.2e}'.format})

def calceq1(tmp,sal,dep):
    tmp_k=np.float64(tmp+273.15)
    pres=np.float64(dep*100.)
    coeff=np.zeros(6,dtype=np.float64)
    coeff[0]=13.4191e0
    coeff[1]=0.0331e0
    coeff[2]=-5.33e-5
    coeff[3]=-530.1228e0
    coeff[4]=-6.103e0
    coeff[5]=-2.06950e0
    calceq1=-126.34048e0+6320.813e0/tmp_k+19.568224*np.log(tmp_k)
    calceq1+=coeff[0]*sal**0.5+coeff[1]*sal+coeff[2]*sal**2. \
              +(coeff[3]*sal**0.5+coeff[4]*sal)/tmp_k  \
              +coeff[5]*sal**0.5*np.log(tmp_k)
    calceq1=10.**(-calceq1)
    # calceq1*=np.exp((-(-15.82e0-0.0219e0*tmp)*pres
                     # +0.5e0*(1.13e-3-0.1475e-3*tmp)*pres*pres)/83.131e0/tmp_k)
    calceq1*=np.exp((-(-25.50e0+0.1271e0*tmp)*pres
                     +0.5e0*(-3.08e-3+0.0877e-3*tmp)*pres*pres)/83.131e0/tmp_k)
    return calceq1

def calceq2(tmp,sal,dep):
    tmp_k=np.float64(tmp+273.15)
    pres=np.float64(dep*100.)
    coeff=np.zeros(6,dtype=np.float64)
    coeff[0]=21.0894e0
    coeff[1]=0.1248e0
    coeff[2]=-0.0003687e0
    coeff[3]=-772.483e0
    coeff[4]=-20.051e0
    coeff[5]=-3.32254e0
    calceq2=-90.18333e0+5143.692e0/tmp_k+14.613358*np.log(tmp_k)
    calceq2+=coeff[0]*sal**0.5+coeff[1]*sal+coeff[2]*sal**2. \
              +(coeff[3]*sal**0.5+coeff[4]*sal)/tmp_k  \
              +coeff[5]*sal**0.5*np.log(tmp_k)
    calceq2=10.**(-calceq2)
    # calceq2*=np.exp((-(-25.50e0+0.1271e0*tmp)*pres
                     # +0.5e0*(-3.08e-3+0.0877e-3*tmp)*pres*pres)/83.131e0/tmp_k)
    calceq2*=np.exp((-(-15.82e0-0.0219e0*tmp)*pres
                     +0.5e0*(1.13e-3-0.1475e-3*tmp)*pres*pres)/83.131e0/tmp_k)
    return calceq2

def calceqcc(tmp,sal,dep):
    tmp_k=np.float64(tmp+273.15e0)
    pres = np.float64(dep*100e0)
    calceqcc = -171.9065e0 - 0.077993e0*tmp_k + 2839.319e0/tmp_k \
               +71.595e0*np.log10(tmp_k) \
               +(-0.77712e0+0.0028426e0*tmp_k+178.34e0/tmp_k)*sal**0.5e0 \
               -0.07711e0*sal+0.0041249e0*sal**1.5e0
    calceqcc = 10e0**calceqcc
    calceqcc*=np.exp((-(-48.76e0+0.5304e0*tmp)*pres
                   +0.5e0*(-11.76e-3+0.3692e-3*tmp)*pres*pres)/83.131e0/tmp_k)
    return calceqcc

def calceqw(tmp,sal,dep):
    tmp_k=np.float64(tmp+273.15e0)
    pres = np.float64(dep*100e0)
    calceqw= 148.96502e0 -13847.26e0/tmp_k -23.6521e0*np.log(tmp_k) \
        + (118.67e0/tmp_k -5.977e0 + 1.0495e0*np.log(tmp_k))*sal**0.5e0 - 0.01615e0*sal
    calceqw=np.exp(calceqw)
    calceqw*=np.exp((-(-25.60e0+0.2324e0*tmp-3.6246e-3*tmp*tmp)*pres+0.5e0*(-5.13e-3+0.0794e-3*tmp)*pres*pres)/83.131e0/tmp_k)
    return calceqw

def chk_caco3_therm():
    temp = 2e0
    sal = 35e0
    n = 100
    for i in range(n):
        dep=6e0*(i+1)/np.float64(n)
        keqcc=calceqcc(temp,sal,dep)
        print dep,keqcc
    return 
    
def chk_caco3_therm_sbrtns():
    nz = 100
    temp = 2e0 # temperature in Celsius 
    sal = 35e0 # salinity wt o/oo 
    dep = 4.0e0 # depth in km
    dic = np.ones(nz,dtype=np.float64)*2211e0*1e-6/1e3  # 2211 uM converted to mol/cm3; all the same value at different grids
    alk = np.ones(nz,dtype=np.float64)*2285e0*1e-6/1e3
    # calling subroutine to calculate all aqueous co2 species and pH
    co2,hco3,co3,ph,info = calcspecies(temp,sal,dep,dic,alk,nz)
    if info!=0:  # if ph cannot be calculated in the subroutine, info=1 is returned 
        print  'error in calcspecies'
    # calling subroutine to calculate derivatives of co3 conc. wrt dic and alk (defined as dco3_ddic and dco3_dalk, respectively)
    dco3_ddic,dco3_dalk,info = calcdevs(temp,sal,dep,dic,alk,nz)
    if info!=0: # if ph cannot be calculated in the subroutine, info=1 is returned 
        print  'error in calcdevs'
    # printing results on screen; if you want check with MATLAB version by copy and paste
    for iz in range(nz):
        print dic[iz],alk[iz],co2[iz],hco3[iz],co3[iz],ph[iz],dco3_ddic[iz],dco3_dalk[iz]
        # note that the above concentrations (alk,dic,co2,hco3,co3) are all in mol/cm3
        # ph is actually H+ concentration in mol/L 
        # dco3_dalk and dco3_ddic should be dimensionless
    return
    
def calcspecies(tmp,sal,dep,dic,alk,nz):
    info = 0
    k1=calceq1(tmp,sal,dep)
    k2=calceq2(tmp,sal,dep)
    a=np.zeros(nz,dtype=np.float64)
    b=np.zeros(nz,dtype=np.float64)
    c=np.zeros(nz,dtype=np.float64)
    pro=np.zeros(nz,dtype=np.float64)
    co2=np.zeros(nz,dtype=np.float64)
    hco3=np.zeros(nz,dtype=np.float64)
    co3=np.zeros(nz,dtype=np.float64)
    a[:]=1e0
    b[:]=(1.-dic[:]/alk[:])*k1
    c[:]=(1.-2.*dic[:]/alk[:])*k1*k2
    pro[:]=(-b[:]+(b[:]*b[:]-4.*a*c[:])**0.5)*0.5
    if any(pro<0.):
           print '... unable to calculate ph'
#           print pro
#           print dic
#           print alk
           info = 1
##    for iz in range(nz):
##        p=np.poly1d([a[iz],b[iz],c[iz]])
##        for ip in range(p.roots[np.isreal(p.roots)].shape[0]):
##            if p.roots[np.isreal(p.roots)][ip]>0:
##                pro[iz]=p.roots[np.isreal(p.roots)][ip]
##                break
    co2[:]=alk[:]/(k1/pro[:]+2.*k1*k2/pro[:]/pro[:])
    hco3[:]=alk[:]/(1.+2.*k2/pro[:])
    co3[:]=alk[:]/(pro[:]/k2+2.)
    db_dalk=np.zeros(nz,dtype=np.float64)
    dc_dalk=np.zeros(nz,dtype=np.float64)
    db_ddic=np.zeros(nz,dtype=np.float64)
    dc_ddic=np.zeros(nz,dtype=np.float64)
    dph_dalk=np.zeros(nz,dtype=np.float64)
    dph_ddic=np.zeros(nz,dtype=np.float64)
    dco3_dalk=np.zeros(nz,dtype=np.float64)
    dco3_ddic=np.zeros(nz,dtype=np.float64)
    db_dalk[:] = k1*(-1.)*dic[:]*(-1.)/alk[:]/alk[:]
    dc_dalk[:] = k1*k2*(-2.)*dic[:]*(-1.)/alk[:]/alk[:]
    db_ddic[:] = k1*(-1./alk[:])
    dc_ddic[:] = k1*k2*(-2./alk[:])
    dph_dalk[:] = -0.5*db_dalk[:] + 0.5*0.5*(b[:]*b[:]-4.*c[:])**(-0.5)\
                  *(2.*b[:]*db_dalk[:] - 4.*dc_dalk[:])
    dph_ddic[:] = -0.5*db_ddic[:] + 0.5*0.5*(b[:]*b[:]-4.*c[:])**(-0.5)\
                  *(2.*b[:]*db_ddic[:] - 4.*dc_ddic[:])
    dco3_dalk[:] = 1./(pro[:]/k2+2.) + alk[:]*(-1.)/((pro[:]/k2+2.)**2.)\
                   *(1./k2)*dph_dalk[:]
    dco3_ddic[:] = 0./(pro[:]/k2+2.) + alk[:]*(-1.)/((pro[:]/k2+2.)**2.)\
                   *(1./k2)*dph_ddic[:]
    return co2,hco3,co3,pro,dco3_ddic,dco3_dalk,info
    
def calcco2h2o(tmp,sal,dep,dic,alk,nz):
    info=0
    k1=calceq1(tmp,sal,dep)
    k2=calceq2(tmp,sal,dep)
    kw=calceqw(tmp,sal,dep)
    a=np.zeros(nz,dtype=np.float64)
    b=np.zeros(nz,dtype=np.float64)
    c=np.zeros(nz,dtype=np.float64)
    d=np.zeros(nz,dtype=np.float64)
    e=np.zeros(nz,dtype=np.float64)
    ph=np.zeros(nz,dtype=np.float64)
    a[:]=1e0
    b[:]=alk[:]+k1
    c[:]=-(dic[:]-alk[:])*k1-kw+k1*k2
    d[:]=-(2e0*dic[:]-alk[:])*k1*k2-k1*kw
    e[:]=-k1*k2*kw
    for iz in range(nz):
        p=np.poly1d([a[iz],b[iz],c[iz],d[iz],e[iz]])
        for ip in range(p.roots[np.isreal(p.roots)].shape[0]):
            if p.roots[np.isreal(p.roots)][ip]>0:
                ph[iz]=p.roots[np.isreal(p.roots)][ip]
                break
    co2=np.zeros(nz,dtype=np.float64)
    hco3=np.zeros(nz,dtype=np.float64)
    co3=np.zeros(nz,dtype=np.float64)
    co2[:] = dic[:]/(1e0+k1/ph[:]+k1*k2/ph[:]/ph[:])
    hco3[:] = dic[:]/(ph[:]/k1+1e0+k2/ph[:])
    co3[:] = dic[:]/(ph[:]*ph[:]/k1/k2+ph[:]/k2+1e0)
    da_dalk=np.zeros(nz,dtype=np.float64)
    db_dalk=np.zeros(nz,dtype=np.float64)
    dc_dalk=np.zeros(nz,dtype=np.float64)
    dd_dalk=np.zeros(nz,dtype=np.float64)
    de_dalk=np.zeros(nz,dtype=np.float64)
    da_ddic=np.zeros(nz,dtype=np.float64)
    db_ddic=np.zeros(nz,dtype=np.float64)
    dc_ddic=np.zeros(nz,dtype=np.float64)
    dd_ddic=np.zeros(nz,dtype=np.float64)
    de_ddic=np.zeros(nz,dtype=np.float64)
    da_dalk[:]=0e0
    db_dalk[:]=1e0
    dc_dalk[:]=k1
    dd_dalk[:]=k1*k2
    de_dalk[:]=0e0
    da_ddic[:]=0e0
    db_ddic[:]=0e0
    dc_ddic[:]=-k1
    dd_ddic[:]=-2e0*k1*k2
    de_ddic[:]=0e0
    dph_dalk=np.zeros(nz,dtype=np.float64)
    dph_ddic=np.zeros(nz,dtype=np.float64)
    dco3_dalk=np.zeros(nz,dtype=np.float64)
    dco3_ddic=np.zeros(nz,dtype=np.float64)
    dph_dalk[:] = -(da_dalk[:]*ph[:]**4e0+db_dalk[:]*ph[:]**3e0+dc_dalk[:]*ph[:]**2e0+dd_dalk[:]*ph[:]+de_dalk[:])  \
        /(a[:]*4e0*ph[:]**3e0+b[:]*3e0*ph[:]**2e0+c[:]*2e0*ph[:]+d[:])
    dph_ddic[:] = -(da_ddic[:]*ph[:]**4e0+db_ddic[:]*ph[:]**3e0+dc_ddic[:]*ph[:]**2e0+dd_ddic[:]*ph[:]+de_ddic[:])  \
        /(a[:]*4e0*ph[:]**3e0+b[:]*3e0*ph[:]**2e0+c[:]*2e0*ph[:]+d[:])        
    dco3_dalk[:] = 0e0/(ph[:]*ph[:]/k1/k2+ph[:]/k2+1e0) \
        + dic[:]*(-1e0)/((ph[:]*ph[:]/k1/k2+ph[:]/k2+1e0)**2e0)*(2e0*ph[:]/k1/k2+1e0/k2)*dph_dalk[:] 
    dco3_ddic[:] = 1e0/(ph[:]*ph[:]/k1/k2+ph[:]/k2+1e0) \
        + dic[:]*(-1e0)/((ph[:]*ph[:]/k1/k2+ph[:]/k2+1e0)**2e0)*(2e0*ph[:]/k1/k2+1e0/k2)*dph_ddic[:] 
    return co2,hco3,co3,ph,dco3_ddic,dco3_dalk,info

def test_co2h2o():
    nz = 100
    dic=np.zeros(nz,dtype=np.float64)
    alk=np.zeros(nz,dtype=np.float64)
    dep = 4.
    sal = 35. 
    temp = 2.
    dic[:] = 2285e0#*1e-6
    alk[:] = 2211e0#*1e-6
    co2,hco3,co3,ph,dco3_ddic,dco3_dalk,info = calcco2h2o(temp,sal,dep,dic,alk,nz)
    print co2
    print ''
    print hco3
    print ''
    print co3
    print ''
    print ph
    print ''
    print dco3_ddic
    print ''
    print dco3_dalk
    print '' 
    
def calcdevs(tmp,sal,dep,dic,alk,nz):
    info = 0
    k1=np.float64(calceq1(tmp,sal,dep))
    k2=np.float64(calceq2(tmp,sal,dep))
    a=np.float64(1.)
    b=np.zeros(nz,dtype=np.float64)
    c=np.zeros(nz,dtype=np.float64)
    pro=np.zeros(nz,dtype=np.float64)
    co2=np.zeros(nz,dtype=np.float64)
    hco3=np.zeros(nz,dtype=np.float64)
    co3=np.zeros(nz,dtype=np.float64)
    b[:]=(1.-dic[:]/alk[:])*k1
    c[:]=(1.-2.*dic[:]/alk[:])*k1*k2
    pro[:]=(-b[:]+(b[:]*b[:]-4.*a*c[:])**0.5)*0.5
    if any(pro<0.):
           print '... unable to calculate ph'
           info = 1
    co2 = np.zeros(nz,dtype=np.float64)
    hco3 = np.zeros(nz,dtype=np.float64)
    co3 = np.zeros(nz,dtype=np.float64)
    co2[:]=alk[:]/(k1/pro[:]+2.*k1*k2/pro[:]/pro[:])
    hco3[:]=alk[:]/(1.+2.*k2/pro[:])
    co3[:]=alk[:]/(pro[:]/k2+2.)
    db_dalk=np.zeros(nz,dtype=np.float64)
    dc_dalk=np.zeros(nz,dtype=np.float64)
    db_ddic=np.zeros(nz,dtype=np.float64)
    dc_ddic=np.zeros(nz,dtype=np.float64)
    dph_dalk=np.zeros(nz,dtype=np.float64)
    dph_ddic=np.zeros(nz,dtype=np.float64)
    dco3_dalk=np.zeros(nz,dtype=np.float64)
    dco3_ddic=np.zeros(nz,dtype=np.float64)
    db_dalk[:] = k1*(-1.)*dic[:]*(-1.)/alk[:]/alk[:]
    dc_dalk[:] = k1*k2*(-2.)*dic[:]*(-1.)/alk[:]/alk[:]
    db_ddic[:] = k1*(-1./alk[:])
    dc_ddic[:] = k1*k2*(-2./alk[:])
    dph_dalk[:] = -0.5*db_dalk[:] + 0.5*0.5*(b[:]*b[:]-4.*c[:])**(-0.5)\
                  *(2.*b[:]*db_dalk[:] - 4.*dc_dalk[:])
    dph_ddic[:] = -0.5*db_ddic[:] + 0.5*0.5*(b[:]*b[:]-4.*c[:])**(-0.5)\
                  *(2.*b[:]*db_ddic[:] - 4.*dc_ddic[:])
    dco3_dalk[:] = 1./(pro[:]/k2+2.) + alk[:]*(-1.)/((pro[:]/k2+2.)**2.)\
                   *(1./k2)*dph_dalk[:]
    dco3_ddic[:] = 0./(pro[:]/k2+2.) + alk[:]*(-1.)/((pro[:]/k2+2.)**2.)\
                   *(1./k2)*dph_ddic[:]
    return dco3_ddic,dco3_dalk,info

def co2sys_mocsy(nz,alk,dic,tempi,depi,sali):
    # following is copied and pasted from test_mocsy.py and then modified
    # Define input data (typical values at depth from 0 to 5000 meters)
    temp = np.repeat(tempi, nz).astype('float32')
    depth = np.repeat (depi, nz).astype('float32')
    sal = np.repeat(sali, nz).astype('float32')
    sil = phos = np.repeat(0.0, nz).astype('float32')
    Patm = np.repeat(1.0, nz).astype('float32')
    # optK1K2 = 'l'
    optcon  = 'mol/m3'  # input concentrations are in MOL/m3
    optt    = 'Tinsitu' # input temperature, variable 'temp' is actually IN SITU temp [°C]
    optp    = 'm'      # input variable 'depth' is in 'DECIBARS'
    optb    = 'l10'
    optk1k2 = 'l'
    optkf   = 'dg'
    # Create output arrays
    # --------------------
    # computed variables at 6 input points
    lat    = np.zeros((nz,)).astype('float32')
    ph     = np.zeros((nz,)).astype('float32')
    pco2   = np.zeros((nz,)).astype('float32')
    fco2   = np.zeros((nz,)).astype('float32')
    co2    = np.zeros((nz,)).astype('float32')
    hco3   = np.zeros((nz,)).astype('float32')
    co3    = np.zeros((nz,)).astype('float32')
    OmegaA = np.zeros((nz,)).astype('float32')
    OmegaC = np.zeros((nz,)).astype('float32')
    BetaD  = np.zeros((nz,)).astype('float32')
    rhoSW  = np.zeros((nz,)).astype('float32')
    p      = np.zeros((nz,)).astype('float32')
    tempis = np.zeros((nz,)).astype('float32')
    # values of derivatives w/ respect to 6 input variables and at 6 input points
    ph_deriv     = np.zeros((6*nz,)).astype('float32')
    pco2_deriv   = np.zeros((6*nz,)).astype('float32')
    OmegaA_deriv = np.zeros((6*nz,)).astype('float32')
    ph_deriv     = ph_deriv.reshape ((6,nz), order='F')
    pco2_deriv   = pco2_deriv.reshape ((6,nz), order='F')
    OmegaA_deriv = OmegaA_deriv.reshape ((6,nz), order='F')
    # Run mocsy.vars()
    # Notice that option names are all lowercase
    ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis = \
        mocsy.mvars (temp, sal, alk, dic, sil, phos, Patm, depth, lat
            ,optcon=optcon, optt=optt, optp=optp, optb=optb, optk1k2=optk1k2, optkf=optkf)
    # print mocsy results
    # -------------------
    # print "pH     pCO2   fCO2     CO2*       HCO3-       CO32-      OmegaA OmegaC  R    Density Press  Temperature"
    # for i in range (0, 3):
        # print ph[i], pco2[i], fco2[i], co2[i], hco3[i], co3[i], OmegaA[i], OmegaC[i], BetaD[i], rhoSW[i], p[i], tempis[i]
    # Compute automatic derivatives (using automatic differentiation)
    ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, hco3_deriv, co3_deriv, OmegaA_deriv, OmegaC_deriv = \
        mocsy.mderivauto(temp, sal, alk, dic, sil, phos, Patm, depth, lat,                                     # INPUT
                        optcon=optcon, optt=optt, optp=optp, optb=optb, optk1k2=optk1k2, optkf=optkf)  # INPUT OPTIONS
    # Compute buffer factors from Egleston
    # pco2_deriv[2,:] are derivatives of pCO2 w/ respect to DIC
    # pco2_deriv[1,:] are        ...          w/ respect to Alk
    gamma_DIC = pco2 / pco2_deriv[1,:]
    gamma_Alk = pco2 / pco2_deriv[0,:]
    beta_DIC  = -1. / (np.log(10.) * ph_deriv[1,:])
    beta_Alk  = -1. / (np.log(10.) * ph_deriv[0,:])
    # Here, we use Omega of Aragonite (use of Calcite would have been equaly valid)
    omega_DIC = OmegaA / OmegaA_deriv[1,:]
    omega_Alk = OmegaA / OmegaA_deriv[0,:]
    # print ""
    # print "gamma_DIC     gamma_Alk     beta_DIC     beta_Alk    omega_DIC    omega_Alk"
    # for i in range (0, 3):
        # print gamma_DIC[i], gamma_Alk[i], beta_DIC[i], beta_Alk[i], omega_DIC[i], omega_Alk[i]
    # Print derivatives of pH with respect to phosphate, silicate, temperature and salinity
    # print ""
    # print "dpH/dPhos  dpH/dSil  dpH/dT   dpH/dS"
    # print "pH/µMol     pH/µMol  pH/°C   pH/psu"
    # for i in range (0, 3):
        # print ph_deriv[2,i], ph_deriv[3,i], ph_deriv[4,i], ph_deriv[5,i]
    return co2,hco3,co3,10.**-ph,OmegaC,OmegaC_deriv[1,:],OmegaC_deriv[0,:]

def calcupwindscheme(w,nz):
    # ------------ determine calculation scheme for advection 
    up=np.zeros(nz,dtype=np.float64)
    dwn=np.zeros(nz,dtype=np.float64)
    cnr=np.zeros(nz,dtype=np.float64)
    adf=np.ones(nz,dtype=np.float64)
    up[:] = 0
    dwn[:]=0
    cnr[:] =0
    adf[:]=1e0
    for iz in range(nz ):
        if iz==0: 
            if w[iz]>=0e0 and w[iz+1]>=0e0:up[iz] = 1
            elif w[iz]<=0e0 and w[iz+1]<=0e0:dwn[iz] = 1
            else:   #  where burial sign changes  
                if not w[iz]*w[iz+1] <=0e0:
                    print 'error'
                    input()
                cnr[iz] = 1
        elif iz==nz-1:
            if w[iz]>=0e0 and w[iz-1]>=0e0:up[iz] = 1
            elif w[iz]<=0e0 and w[iz-1]<=0e0:dwn[iz] = 1
            else: 
                if not w[iz]*w[iz-1] <=0e0:
                    print 'error'
                    input()
                cnr[iz] = 1
        else :
            if w[iz] >=0e0:
                if w[iz+1]>=0e0 and w[iz-1]>=0e0:up[iz] = 1
                else:cnr[iz] = 1
            else:  
                if w[iz+1]<=0e0 and w[iz-1]<=0e0:dwn[iz] = 1
                else:cnr[iz] = 1
    if np.sum(up[:]+dwn[:]+cnr[:])!=nz:
        print 'error',np.sum(up),np.sum(dwn),np.sum(cnr)
        input()
    for iz in range(nz-1):
        if cnr[iz]==1 and cnr[iz+1]==1:
            if w[iz]>=0e0 and w[iz+1] < 0e0:
                corrf = np.float64(5e0)  #  This assignment of central advection term helps conversion especially when assuming turbo2 mixing 
                cnr[iz+1]=abs(w[iz]**corrf)/(abs(w[iz+1]**corrf)+abs(w[iz]**corrf))
                cnr[iz]=abs(w[iz+1]**corrf)/(abs(w[iz+1]**corrf)+abs(w[iz]**corrf))
                dwn[iz+1]=1e0-cnr[iz+1]
                up[iz]=1e0-cnr[iz]
        if cnr[iz]==1 and cnr[iz+1]==1: 
            if w[iz]< 0e0 and w[iz+1] >= 0e0:
                cnr[iz+1]=0
                cnr[iz]=0
                up[iz+1]=1
                dwn[iz]=1
                adf[iz]=abs(w[iz+1])/(abs(w[iz+1])+abs(w[iz]))
                adf[iz+1]=abs(w[iz])/(abs(w[iz+1])+abs(w[iz]))
    return up,dwn,cnr,adf
    
def makegrid(beta,ztot,nz):
    z = np.zeros(nz,dtype=np.float64)
    dz = np.zeros(nz,dtype=np.float64)
    eta = np.linspace(ztot/nz,ztot,nz)
    for iz in range(nz):
        if iz==0:dz[iz]=ztot*np.log((beta+(eta[iz]/ztot)**2) 
                                 /(beta-(eta[iz]/ztot)**2)) \
                                 /np.log((beta+1.)/(beta-1.))
        else:dz[iz] = ztot*np.log((beta+(eta[iz]/ztot)**2) 
                                  /(beta-(eta[iz]/ztot)**2)) \
                                  /np.log((beta+1.)/(beta-1.)) - np.sum(dz[:iz])
    for iz in range(nz):
        if iz==0: z[iz]=dz[iz]*0.5
        else: z[iz] = z[iz-1]+dz[iz-1]*0.5 + 0.5*dz[iz]
    return dz,z
    
def getporosity(z,nz):
    poro = np.zeros(nz,dtype=np.float64)
    calgg = 0.0
    pore_max = 1.- ( 0.483 + 0.45 * calgg) / 2.5  # porosity at the bottom 
    exp_pore = 0.25*calgg + 3.0 *(1.-calgg)# scale depth of e-fold porosity decrease 
    poro[:] = np.exp(-z[:]/exp_pore) * (1.0-pore_max) + pore_max 
    return poro

def dep2age(z,dz,nz,w):
    age = np.zeros(nz,dtype=np.float64)
    dage = np.zeros(nz,dtype=np.float64)
    dage[:] = dz[:]/w[:]
    for iz in range(nz):
        if iz==0: age[iz]=dage[iz]*0.5
        else: age[iz]=age[iz-1]+dage[iz-1]*0.5+0.5*dage[iz]
    return age
    
def coefs(cai,temp,nz,nspcc,poro,komi,kcci,size,sal,dep):
    dif_dic = np.zeros(nz,dtype=np.float64)
    dif_alk = np.zeros(nz,dtype=np.float64)
    dif_o2 = np.zeros(nz,dtype=np.float64)
    kom = np.zeros(nz,dtype=np.float64)
    kcc = np.zeros((nz,nspcc),dtype=np.float64)
    co3sat = 0.
    ff = np.zeros(nz,dtype=np.float64)
    ff[:] = poro[:]*poro[:]
    dif_dic0 = (151.69 + 7.93*temp) # cm2/yr at 2 oC (Huelse et al. 2018)
    dif_alk0 = (151.69 + 7.93*temp) # cm2/yr  (Huelse et al. 2018)
    dif_o20 =  (348.62 + 14.09*temp) 
    dif_dic[:] = dif_dic0*ff[:]  # reflecting tortuosity factor 
    dif_alk[:] = dif_alk0*ff[:]
    dif_o2[:] = dif_o20*ff[:]
    kom[:] = komi  # assume reference values for all reaction terms 
    kcc[:,:] = kcci
    if size:
        kcc[:,0:4]=kcci*10.
    # keq1 = calceq1(temp,sal,dep)
    # keq2 = calceq2(temp,sal,dep)
    keqcc = calceqcc(temp,sal,dep)
    co3sat = keqcc/cai
    return dif_dic,dif_alk,dif_o2,kom,kcc,co3sat
    
def calc_zox(
    oxic,anoxic,nz,o2x,o2th,komi,ztot,z,o2i,dz  # input
    ):   
    """calculation of zox""" 
    tol = 1e-6
    zox = 0e0
    o2_penetrated = True
    for iz in range(nz):
        if o2x[iz]<=0e0: 
            o2_penetrated = False
            break
    if o2_penetrated: # oxygen never gets less than 0 
        zox = ztot # zox is the bottom depth 
    else: 
        if iz==0: # calculating zox interpolating at z=0 with SWI conc. and at z=z(iz) with conc. o2x(iz)
            zox = (z[iz]*o2i*1e-6/1e3 + 0e0*np.abs(o2x[iz]))/(o2i*1e-6/1e3+np.abs(o2x[iz]))
        elif iz==1:  
            zox = z[iz-1] - o2x[iz-1]/((o2i*1e-6/1e3 - o2x[iz-1])/(0e0-z[iz-1]))
        else:     # calculating zox interpolating at z=z(iz-1) with o2x(iz-1) and at z=z(iz) with conc. o2x(iz)
            zox = z[iz-1] - o2x[iz-1]/((o2x[iz-2] - o2x[iz-1])/(z[iz-2]-z[iz-1]))
    # calculation of kom 
    kom = np.zeros(nz,dtype=np.float64)
    kom_ox = np.zeros(nz,dtype=np.float64)
    kom_an = np.zeros(nz,dtype=np.float64)
    izox = -100 
    if anoxic: 
        kom[:] = komi
        for iz in range(nz):
            if z[iz]+0.5e0*dz[iz]<=zox: 
                kom_ox[iz]=komi
                if iz> izox:  izox = iz
            elif z[iz]+0.5e0*dz[iz]>zox and z[iz]-0.5e0*dz[iz]< zox: 
                kom_ox[iz]=komi* (1e0- ( (z[iz]+0.5e0*dz[iz]) - zox)/dz[iz])
                kom_an[iz]=komi* (( (z[iz]+0.5e0*dz[iz]) - zox)/dz[iz])
                if iz> izox: izox = iz
            elif z[iz]-0.5e0*dz[iz]>=zox: 
                kom_an[iz]=komi
        if not (abs(kom_ox[:]+kom_an[:]-kom[:])/komi<tol).all(): 
            print 'error: calc kom',kom_ox,kom_an,kom
            input()
    else:
        for iz in range(nz):
            if z[iz]+0.5e0*dz[iz]<=zox: 
                kom_ox[iz]=komi
                if iz> izox: izox = iz
            elif z[iz]+0.5e0*dz[iz]>zox and z[iz]-0.5e0*dz[iz]< zox: 
                kom_ox[iz]=komi* (1e0- ( (z[iz]+0.5e0*dz[iz]) - zox)/dz[iz])
                if iz> izox :izox = iz
            elif z[iz]-0.5e0*dz[iz]>=zox :
                continue
        kom[:] = kom_ox[:]
    return izox,kom,zox,kom_ox,kom_an
    
def omcalc( 
    kom,omx     # in&output
    ,oxic,anoxic,o2x,om,nz,sporo,sporoi,sporof,o2th,komi  # input 
    ,w,wi,dt,up,dwn,cnr,adf,trans,nspcc,labs,turbo2,nonlocal,omflx,poro,dz  # input 
    ):
    nsp = 1
    nmx=nz*nsp
    amx=np.zeros((nmx,nmx),dtype=np.float64)
    ymx=np.zeros((nmx),dtype=np.float64)
    for iz in range(nz):
        row = (iz)*nsp
        if iz==0:
            ymx[row] = \
                + sporo[iz]*(-om[iz])/dt \
                - omflx/dz[iz]
            amx[row,row] = (
                # time change term 
                + sporo[iz]*(1e0)/dt 
                # advection terms 
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*1e0-sporoi*wi*0e0)/dz[iz]   
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*0e0
                                   -sporo[iz]*w[iz]*1e0)/dz[iz]   
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*0e0
                                   -sporoi*wi*0e0)/dz[iz]   
                #  rxn term 
                + sporo[iz]*kom[iz]   
                )
            # matrix filling at grid iz but for unknwon at grid iz + 1
            # (here only advection terms) 
            amx[row,row+nsp] =  (
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*1e0-sporo[iz]*w[iz]*0e0)/dz[iz]   
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*1e0-sporoi*wi*0e0)/dz[iz]   
                )
        elif iz == nz-1:# : ! need to reflect lower boundary; none; but must care that iz + 1 does not exist 
            ymx[row] = 0e0   \
                + sporo[iz]*(-om[iz])/dt 
            amx[row,row] = (
                # time change term 
                + sporo[iz]*(1e0)/dt 
                # advection terms 
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*1e0-sporo[iz-1]*w[iz-1]*0e0)/dz[iz]  
                + adf[iz]*dwn[iz]*(sporof*w[iz]*1e0-sporo[iz]*w[iz]*1e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(sporof*w[iz]*1e0-sporo[iz-1]*w[iz-1]*0e0)/dz[iz]  
                # rxn term 
                + sporo[iz]*kom[iz]   
                )
            # filling matrix at grid iz but for unknown at grid iz-1 (only advection terms) 
            amx[row,row-nsp] = ( 
                + adf[iz]*up[iz]*(sporof*w[iz]*0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(sporof*w[iz]*0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                )
        else:# # do not have to care about boundaries 
            ymx[row] = 0e0  \
                + sporo[iz]*(0e0-om[iz])/dt 
            amx[row,row] = (
                # time change term 
                + sporo[iz]*(1e0)/dt 
                # advection terms 
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*1e0-sporo[iz-1]*w[iz-1]*0e0)/dz[iz]  
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*0e0-sporo[iz]*w[iz]*1e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*0e0-sporo[iz-1]*w[iz-1]*0e0)/dz[iz]  
                # rxn term 
                + sporo[iz]*kom[iz]   
                )
            # filling matrix at grid iz but for unknown at grid iz+1 (only advection terms) 
            amx[row,row+nsp] =  (
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*1e0-sporo[iz]*w[iz]*0e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*1e0-sporo[iz-1]*w[iz-1]*0e0)/dz[iz]  
                )
            # filling matrix at grid iz but for unknown at grid iz-1 (only advection terms) 
            amx[row,row-nsp] =  (
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                )
        # diffusion terms are reflected with transition matrices 
        if turbo2[0] or labs[0]: 
            for iiz in range( nz):
                col =  (iiz)*nsp 
                if trans[iiz,iz,0]==0e0: continue
                amx[row,col] += -trans[iiz,iz,0]/dz[iz]*dz[iiz]*(1e0-poro[iiz])
        else :
            for iiz in range( nz):
                col =  (iiz)*nsp 
                if trans[iiz,iz,0]==0e0: continue
                amx[row,col] += -trans[iiz,iz,0]/dz[iz]
    ymx = - ymx 
    # kai = np.linalg.solve(amx, ymx)
    # kai = la.solve(amx, ymx)
    lu, piv, kai, info = lapack.dgesv(amx, ymx)
    ymx[:]=kai[:].copy()
    omx[:] = ymx[:].copy() # now passing the solution to unknowns omx 
    return omx
    
def calcflxom(  
    sporo,om,omx,dt,w,dz,z,nz,turbo2,labs,nonlocal,poro
    ,up,dwn,cnr,adf,rho,mom,trans,kom,sporof,sporoi,wi,nspcc,omflx,workdir  # input 
    ):
    """ calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations) """
    omadv = np.float64(0e0)
    omdec = np.float64(0e0)
    omdif = np.float64(0e0)
    omrain = np.float64(0e0)
    omtflx = np.float64(0e0)
    for iz in range(nz ):
        if iz == 0: 
            omtflx +=  sporo[iz]*(omx[iz]-om[iz])/dt*dz[iz] 
            omadv +=  adf[iz]*up[iz]*(sporo[iz]*w[iz]*omx[iz]-0e0)/dz[iz]*dz[iz]  \
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*omx[iz+1]-sporo[iz]*w[iz]*omx[iz])/dz[iz]*dz[iz]  \
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*omx[iz+1]-0e0)/dz[iz]*dz[iz]  
            omdec +=  sporo[iz]*kom[iz]*omx[iz]*dz[iz]
            omrain += - omflx/dz[iz]*dz[iz]
        elif iz == nz-1: 
            omtflx = omtflx + sporo[iz]*(omx[iz]-om[iz])/dt*dz[iz]
            omadv = omadv + adf[iz]*up[iz]*(sporo[iz]*w[iz]*omx[iz]-sporo[iz-1]*w[iz-1]*omx[iz-1])/dz[iz]*dz[iz] \
                + adf[iz]*dwn[iz]*(sporof*w[iz]*omx[iz]-sporo[iz]*w[iz]*omx[iz])/dz[iz]*dz[iz] \
                + adf[iz]*cnr[iz]*(sporof*w[iz]*omx[iz]-sporo[iz-1]*w[iz-1]*omx[iz-1])/dz[iz]*dz[iz] 
            omdec = omdec + sporo[iz]*kom[iz]*omx[iz]*dz[iz]
        else :
            omtflx = omtflx + sporo[iz]*(omx[iz]-om[iz])/dt*dz[iz]
            omadv = omadv + adf[iz]*up[iz]*(sporo[iz]*w[iz]*omx[iz]-sporo[iz-1]*w[iz-1]*omx[iz-1])/dz[iz]*dz[iz]  \
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*omx[iz+1]-sporo[iz]*w[iz]*omx[iz])/dz[iz]*dz[iz]  \
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*omx[iz+1]-sporo[iz-1]*w[iz-1]*omx[iz-1])/dz[iz]*dz[iz] 
            omdec = omdec + sporo[iz]*kom[iz]*omx[iz]*dz[iz]
        if turbo2[0] or labs[0]:  
            for  iiz in range( nz):
                if trans[iiz,iz,0]==0e0: continue
                omdif = omdif -trans[iiz,iz,0]/dz[iz]*dz[iiz]*(1e0-poro[iiz])*dz[iz]*omx[iiz]
        else:
            for iiz in range(nz):
                if trans[iiz,iz,0]==0e0: continue
                omdif = omdif -trans[iiz,iz,0]/dz[iz]*dz[iz]*omx[iiz]  # check previous versions 
    omres = omadv + omdec + omdif + omrain + omtflx # this is residual flux should be zero equations are exactly satisfied 
    flg_restart = False
    if any(omx<0e0):   # if negative om conc. is detected, need to stop  
        print 'negative om, stop'
        file_tmp=open(workdir+'NEGATIVE_OM.txt','w')
        for  iz in range(nz):
            print >> file_tmp, z[iz],omx[iz]*mom/rho[iz]*100e0,w[iz],up[iz],dwn[iz],cnr[iz],adf[iz]
        file_tmp.close()
        flg_restart = True
        # input()
    if any(np.isnan(omx)) :  # if NAN, ... the same ... stop
        print 'nan om, stop'
        print omx
        flg_restart = True
        # input()
    return omadv,omdec,omdif,omrain,omres,omtflx,flg_restart

def o2calc_ox(  
    izox,nz,poro,o2,kom,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
    ):
    nsp = 1
    nmx=nsp*nz
    amx=np.zeros((nmx,nmx),dtype=np.float64)
    ymx=np.zeros((nmx),dtype=np.float64)
    amx[:] = 0e0
    ymx[:] = 0e0 
    for iz in range(nz): 
        row =  (iz)*nsp 
        if iz == 0 : # be careful about upper boundary 
            ymx[row] = ( 
                + poro[iz]*(0e0-o2[iz])/dt 
                - ((poro[iz]*dif_o2[iz]+poro[iz+1]*dif_o2[iz+1])*0.5e0*(0e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                - poro[iz]*dif_o2[iz]*(0e0-o2i*1e-6/1e3)/dz[iz])/dz[iz]  
                # rxn term 
                + sporo[iz]*ox2om*kom[iz]*omx[iz]  
                )
            amx[row,row] = (
                + poro[iz]*(1e0)/dt 
                # diffusion term 
                - ((poro[iz]*dif_o2[iz]+poro[iz+1]*dif_o2[iz+1])*0.5e0*(-1e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                - poro[iz]*dif_o2[iz]*(1e0)/dz[iz])/dz[iz]
                )
            # filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term) 
            amx[row,row+nsp] = (
                - ((poro[iz]*dif_o2[iz]+poro[iz+1]*dif_o2[iz+1])*0.5e0*(1e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                - 0e0)/dz[iz]
                )
        elif iz == nz-1 : # be careful about lower boundary 
            ymx[row] = (0e0 
                # time change term 
                + poro[iz]*(0e0-o2[iz])/dt 
                # diffusion term 
                - (0e0 - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(0e0)/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] 
                # rxn term 
                + sporo[iz]*ox2om*kom[iz]*omx[iz]  
                )
            amx[row,row] = ( 
                # time change term 
                + poro[iz]*(1e0)/dt 
                # diffusion term 
                - (0e0 - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(1e0)/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] 
                )
            # filling matrix at grid iz but for unknown at grid iz-1 (only diffusion term) 
            amx[row,row-nsp] = ( 
                - (0e0 - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(-1e0)/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] 
                ) 
        else :
            ymx[row] = ( 0e0
                # time change term 
                + poro[iz]*(0e0-o2[iz])/dt 
                # diffusion term 
                - (0.5e0*(poro[iz+1]*dif_o2[iz+1]+poro[iz]*dif_o2[iz])*(0e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(0e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]  
                # rxn term 
                + sporo[iz]*ox2om*kom[iz]*omx[iz]  
                )
            amx[row,row] = (
                # time change term 
                + poro[iz]*(1e0)/dt 
                # diffusion term 
                - (0.5e0*(poro[iz+1]*dif_o2[iz+1]+poro[iz]*dif_o2[iz])*(-1e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(1e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]  
                )
            # filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term) 
            amx[row,row+nsp] = (
                - (0.5e0*(poro[iz+1]*dif_o2[iz+1]+poro[iz]*dif_o2[iz])*(1e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                - 0e0)/dz[iz]  
                )
            # filling matrix at grid iz but for unknown at grid iz-1 (only diffusion term) 
            amx[row,row-nsp]= (
                - (0e0 
                - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(-1e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz] 
                )
    ymx = - ymx  # sign change; see above for the case of om 
    # ymx = np.linalg.solve(amx, ymx)
    # kai = la.solve(amx, ymx)
    lu, piv, kai, info = lapack.dgesv(amx, ymx)
    ymx[:]=kai[:].copy()
    o2x[:] = ymx[:].copy() # passing solutions to unknowns 
    return o2x

def calcflxo2_ox( 
    nz,sporo,kom,omx,dz,poro,dif_o2,dt,o2,o2x,ox2om,o2i  # input
    ):
    o2dec = np.float64(0e0)
    o2dif = np.float64(0e0)
    o2tflx = np.float64(0e0)
    for iz in range(nz):
        if iz == 0 : 
            o2dec = o2dec + sporo[iz]*ox2om*kom[iz]*omx[iz]*dz[iz]
            o2tflx = o2tflx + (o2x[iz]-o2[iz])/dt*dz[iz]*poro[iz]
            o2dif = o2dif - ((poro[iz]*dif_o2[iz]+poro[iz+1]*dif_o2[iz+1])*0.5e0*(o2x[iz+1]-o2x[iz])/(0.5e0*(dz[iz]+dz[iz+1])) \
                - poro[iz]*dif_o2[iz]*(o2x[iz]-o2i*1e-6/1e3)/dz[iz])/dz[iz] *dz[iz]
        elif iz == nz-1 : 
            o2dec = o2dec + (1e0-poro[iz])*ox2om*kom[iz]*omx[iz]/poro[iz]*dz[iz]*poro[iz]
            o2tflx = o2tflx + (o2x[iz]-o2[iz])/dt*dz[iz]*poro[iz]
            o2dif = o2dif \
                - (0e0 - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(o2x[iz]-o2x[iz-1])/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] \
                *dz[iz]
        else :
            o2dec = o2dec + (1e0-poro[iz])*ox2om*kom[iz]*omx[iz]/poro[iz]*dz[iz]*poro[iz]
            o2tflx = o2tflx + (o2x[iz]-o2[iz])/dt*dz[iz]*poro[iz]
            o2dif = o2dif \
                - (0.5e0*(poro[iz+1]*dif_o2[iz+1]+poro[iz]*dif_o2[iz])*(o2x[iz+1]-o2x[iz])/(0.5e0*(dz[iz+1]+dz[iz])) \
                - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(o2x[iz]-o2x[iz-1])/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]  \
                *dz[iz]
    o2res = o2dec + o2dif + o2tflx  # residual flux 
    return o2dec,o2dif,o2tflx,o2res
    
def o2calc_sbox(  
    izox,nz,poro,o2,kom,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
    ):
#    print izox
    # if izox!=nz-1:iz_zero = izox  # modification made only for Python 
    # else:iz_zero = izox-1  # modification made only for Python 
    iz_zero = izox
    nsp = 1
    nmx=nsp*nz
    amx=np.zeros((nmx,nmx),dtype=np.float64)
    ymx=np.zeros((nmx),dtype=np.float64)
    amx[:] = 0e0
    ymx[:] = 0e0 
    for iz in range(nz ):
        row = (iz)*nsp 
        if iz == 0 : 
            ymx[row] = ( 
                # time change 
                + poro[iz]*(0e0-o2[iz])/dt 
                # diffusion 
                - ((poro[iz]*dif_o2[iz]+poro[iz+1]*dif_o2[iz+1])*0.5e0*(0e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                - poro[iz]*dif_o2[iz]*(0e0-o2i*1e-6/1e3)/dz[iz])/dz[iz]  
                # rxn 
                + sporo[iz]*ox2om*kom[iz]*omx[iz]  
                )
            amx[row,row] = (
                # time change 
                + poro[iz]*(1e0)/dt 
                # diffusion 
                - ((poro[iz]*dif_o2[iz]+poro[iz+1]*dif_o2[iz+1])*0.5e0*(-1e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                - poro[iz]*dif_o2[iz]*(1e0)/dz[iz])/dz[iz]
                )
            # filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term) 
            amx[row,row+nsp] = (
                - ((poro[iz]*dif_o2[iz]+poro[iz+1]*dif_o2[iz+1])*0.5e0*(1e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                - 0e0)/dz[iz]
                )
        elif iz>0 and iz<= iz_zero : 
            ymx[row] = ( 0e0
                # time change 
                + poro[iz]*(0e0-o2[iz])/dt 
                # diffusion 
                - (0.5e0*(poro[iz+1]*dif_o2[iz+1]+poro[iz]*dif_o2[iz])*(0e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(0e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]  
                # rxn 
                + sporo[iz]*ox2om*kom[iz]*omx[iz]  
                )
            amx[row,row] = (
                # time change 
                + poro[iz]*(1e0)/dt 
                # diffusion
                - (0.5e0*(poro[iz+1]*dif_o2[iz+1]+poro[iz]*dif_o2[iz])*(-1e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(1e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]  
                )
            # filling matrix at grid iz but for unknown at grid iz+1 (only diffusion term) 
            amx[row,row+nsp] = (
                - (0.5e0*(poro[iz+1]*dif_o2[iz+1]+poro[iz]*dif_o2[iz])*(1e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                - 0e0)/dz[iz]  
                )
            # filling matrix at grid iz but for unknown at grid iz-1 (only diffusion term) 
            amx[row,row-nsp] = (
                - (0e0 
                - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(-1e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz] 
                )
        elif iz> iz_zero :  # at lower than zox; zero conc. is forced 
            ymx[row] = ( 0e0
                )
            amx[row,row] = (
                + 1e0 
                )
    ymx = - ymx  # change signs 
    # kai = np.linalg.solve(amx, ymx) # solving
    # kai = la.solve(amx, ymx) # solving
    lu, piv, kai, info = lapack.dgesv(amx, ymx)
    ymx[:]=kai[:].copy()
    o2x[:] = ymx[:].copy() # passing solution to variable 
    return o2x 
    
def calcflxo2_sbox( 
    nz,sporo,kom,omx,dz,poro,dif_o2,dt,o2,o2x,izox,ox2om,o2i  # input
    ):
    # if izox!=nz-1:iz_zero = izox  # modification made only for Python 
    # else:iz_zero = izox-1  # modification made only for Python 
    iz_zero = izox
    # fluxes relevant to oxygen 
    o2dec = np.float64(0e0)
    o2dif = np.float64(0e0)
    o2tflx = np.float64(0e0)
    for iz in range(nz ):
        if iz == 0: 
            o2dec = o2dec + sporo[iz]*ox2om*kom[iz]*omx[iz]*dz[iz]
            o2tflx = o2tflx + (o2x[iz]-o2[iz])/dt*dz[iz]*poro[iz]
            o2dif = o2dif - ((poro[iz]*dif_o2[iz]+poro[iz+1]*dif_o2[iz+1])*0.5e0*(o2x[iz+1]-o2x[iz])/(0.5e0*(dz[iz]+dz[iz+1])) 
                - poro[iz]*dif_o2[iz]*(o2x[iz]-o2i*1e-6/1e3)/dz[iz])/dz[iz] *dz[iz]
        elif iz>0 and iz<=iz_zero : 
            o2dec = o2dec + (1e0-poro[iz])*ox2om*kom[iz]*omx[iz]/poro[iz]*dz[iz]*poro[iz]
            o2tflx = o2tflx + (o2x[iz]-o2[iz])/dt*dz[iz]*poro[iz]
            o2dif = o2dif \
                - (0.5e0*(poro[iz+1]*dif_o2[iz+1]+poro[iz]*dif_o2[iz])*(o2x[iz+1]-o2x[iz])/(0.5e0*(dz[iz+1]+dz[iz])) 
                - 0.5e0*(poro[iz]*dif_o2[iz]+poro[iz-1]*dif_o2[iz-1])*(o2x[iz]-o2x[iz-1])/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]  \
                *dz[iz]
    o2res = o2dec + o2dif + o2tflx  # residual flux 
    return o2dec,o2dif,o2tflx,o2res

def calccaco3sys(  #
    ccx,dicx,alkx,rcc,dt  # in&output
    ,nspcc,dic,alk,dep,sal,temp,labs,turbo2,nonlocal,sporo,sporoi,sporof,poro,dif_alk,dif_dic # input
    ,w,up,dwn,cnr,adf,dz,trans,cc,oxco2,anco2,co3sat,kcc,ccflx,ncc,omega,nz,tol,sparse,fact  # input
    ,dici,alki,ccx_th,showiter,w_pre,co2chem,mcc,rho,workdir
    ):
    drcc_dco3=np.zeros((nz,nspcc),dtype=np.float64)
    drcc_dcc=np.zeros((nz,nspcc),dtype=np.float64)
    drcc_ddic=np.zeros((nz,nspcc),dtype=np.float64)
    drcc_dalk=np.zeros((nz,nspcc),dtype=np.float64)
    dco3_dalk=np.zeros((nz),dtype=np.float64)
    dco3_ddic=np.zeros((nz),dtype=np.float64)
    error = 1e4
    itr = 0
    nsp = 2 + nspcc  # now considered species are dic, alk and nspcc of caco3 
    nmx = nz*nsp  # col (and row) of matrix; the same number of unknowns 
    flg_restart = False
    while error > tol:
        if flg_restart: break # exit from loop for caco3 system iteration 
        if sparse:amx=lil_matrix((nmx,nmx),dtype=np.float64)
        else:amx=np.zeros((nmx,nmx),dtype=np.float64)
        amx[:,:]=0e0
        ymx=np.zeros(nmx,dtype=np.float64)
        # calling subroutine from caco3_therm.f90 to calculate aqueous co2 species 
        if co2chem != 'mocsy':
            if co2chem == 'co2':
                co2x,hco3x,co3x,prox,dco3_ddic,dco3_dalk,info = calcspecies(temp,sal,dep,dicx,alkx,nz)
                # print info
                # calcspecies(temp,sal,dep,dicx,alkx)
                if info==1: # which means error in calculation 
                    dt=dt/10e0
                    dicx[:]=dic[:].copy()
                    alkx[:]=alk[:].copy()
                    ccx[:,:]=cc[:,:].copy()
                    w[:]=w_pre[:].copy()
                    # upwind(w)
                    up,dwn,cnr,adf = calcupwindscheme(w,nz)
                    print 'location 1'
                    flg_restart = True 
                    continue
                    print 'stop'
                    input()
            elif co2chem=='co2h2o':
                co2x,hco3x,co3x,prox,dco3_ddic,dco3_dalk,info = calcco2h2o(temp,sal,dep,dicx,alkx,nz)
                if info==1: # which means error in calculation 
                    dt=dt/10e0
                    dicx[:]=dic[:].copy()
                    alkx[:]=alk[:].copy()
                    ccx[:,:]=cc[:,:].copy()
                    w[:]=w_pre[:].copy()
                    # upwind(w)
                    up,dwn,cnr,adf = calcupwindscheme(w,nz)
                    print 'location 2'
                    flg_restart = True 
                    continue
                    print 'stop'
                    input()
            for isp in range(nspcc):
                # calculation of dissolution rate for individual species 
                rcc[:,isp] = kcc[:,isp]*ccx[:,isp]*abs(1e0-co3x[:]*1e3/co3sat)**ncc*((1e0-co3x[:]*1e3/co3sat)>0e0).astype(float)
                # calculation of derivatives of dissolution rate wrt conc. of caco3 species, dic and alk 
                drcc_dcc[:,isp] = kcc[:,isp]*abs(1e0-co3x[:]*1e3/co3sat)**ncc*((1e0-co3x[:]*1e3/co3sat)>0e0).astype(float)
                drcc_dco3[:,isp] = kcc[:,isp]*ccx[:,isp]*ncc*abs(1e0-co3x[:]*1e3/co3sat)**(ncc-1e0)  \
                    *((1e0-co3x[:]*1e3/co3sat)>0e0).astype(float)*(-1e3/co3sat)
                drcc_ddic[:,isp] = drcc_dco3[:,isp]*dco3_ddic[:]
                drcc_dalk[:,isp] = drcc_dco3[:,isp]*dco3_dalk[:]
        elif co2chem == 'mocsy': 
            co2x,hco3x,co3x,prox,ohmega,dohmega_ddic,dohmega_dalk = co2sys_mocsy(nz,alkx*1e6,dicx*1e6,temp,dep*1e3,sal)
            co2x = co2x/1e6
            hco3x = hco3x/1e6
            co3x = co3x/1e6
            dohmega_ddic = dohmega_ddic*1e6
            dohmega_dalk = dohmega_dalk*1e6
            drcc_dohmega = np.zeros((nz,nspcc),dtype=np.float64)
            for isp in range(nspcc):
                # calculation of dissolution rate for individual species 
                rcc[:,isp] = kcc[:,isp]*ccx[:,isp]*abs(1e0-ohmega[:])**ncc*((1e0-ohmega[:])>0e0).astype(float) 
                # calculation of derivatives of dissolution rate wrt conc. of caco3 species, dic and alk 
                drcc_dcc[:,isp] = kcc[:,isp]*abs(1e0-ohmega[:])**ncc*((1e0-ohmega[:])>0e0).astype(float)
                drcc_dohmega[:,isp] = kcc[:,isp]*ccx[:,isp]*ncc*abs(1e0-ohmega[:])**(ncc-1e0)  \
                    *((1e0-ohmega[:])>0e0).astype(float)*(-1e0)
                drcc_ddic[:,isp] = drcc_dohmega[:,isp]*dohmega_ddic[:]
                drcc_dalk[:,isp] = drcc_dohmega[:,isp]*dohmega_dalk[:]
        for iz in range(nz):
            row =(iz)*nsp 
            if iz == 0: # when upper condition must be taken account; *** comments for matrix filling are given only in this case 
                for isp in range(nspcc):  # multiple caco3 species 
                    # put f(x) for isp caco3 species 
                    ymx[row+isp] = \
                        + sporo[iz]*(ccx[iz,isp]-cc[iz,isp])/dt \
                        - ccflx[isp]/dz[iz] \
                        + adf[iz]*up[iz]*(sporo[iz]*w[iz]*ccx[iz,isp]-0e0)/dz[iz]  \
                        + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*ccx[iz+1,isp]-sporo[iz]*w[iz]*ccx[iz,isp])/dz[iz]  \
                        + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*ccx[iz+1,isp]-0e0)/dz[iz]  \
                        + sporo[iz]*rcc[iz,isp]
                    # derivative of f(x) wrt isp caco3 conc. at grid iz in ln 
                    amx[row+isp,row+isp] = (
                        + sporo[iz]*(1e0)/dt 
                        + adf[iz]*up[iz]*(sporo[iz]*w[iz]*1e0-0e0)/dz[iz]   
                        + adf[iz]*dwn[iz]*(0e0-sporo[iz]*w[iz]*1e0)/dz[iz]  
                        + sporo[iz]* drcc_dcc[iz,isp]  
                        )* ccx[iz,isp] 
                    # derivative of f(x) wrt isp caco3 conc. at grid iz+1 in ln 
                    amx[row+isp,row+isp+nsp] =  (
                        + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*1e0-0e0)/dz[iz]  
                        + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*1e0-0e0)/dz[iz]  
                        )*ccx[iz+1,isp]
                    # derivative of f(x) wrt dic conc. at grid iz in ln
                    amx[row+isp,row+nspcc] = (
                        + sporo[iz]*drcc_ddic[iz,isp]  
                        )*dicx[iz]
                    # derivative of f(x) wrt alk conc. at grid iz in ln
                    amx[row+isp,row+nspcc+1] = (
                        + sporo[iz]*drcc_dalk[iz,isp]  
                        )*alkx[iz]
                    # DIC
                    # derivative of f(x) for dic at iz wrt isp caco3 conc. at grid iz in ln
                    amx[row+nspcc,row+isp] = (
                        - (1e0-poro[iz])*drcc_dcc[iz,isp]  
                        )*ccx[iz,isp]*fact
                    # ALK 
                    # derivative of f(x) for alk at iz wrt isp caco3 conc. at grid iz in ln
                    amx[row+nspcc+1,row+isp] = (
                        - 2e0* (1e0-poro[iz])*drcc_dcc[iz,isp]  
                        )*ccx[iz,isp]*fact
                #  DIC 
                # put f(x) for dic at iz  
                ymx[row+nspcc] = ( 
                    + poro[iz]*(dicx[iz]-dic[iz])/dt 
                    - ((poro[iz]*dif_dic[iz]+poro[iz+1]*dif_dic[iz+1])*0.5e0*(dicx[iz+1]-dicx[iz])/(0.5e0*(dz[iz]+dz[iz+1])) 
                    - poro[iz]*dif_dic[iz]*(dicx[iz]-dici*1e-6/1e3)/dz[iz])/dz[iz]  
                    - oxco2[iz] 
                    - anco2[iz] 
                    - (1e0-poro[iz])*np.sum(rcc[iz,:])  
                    )*fact
                # put derivative of f(x) for dic at iz wrt dic at iz in ln 
                amx[row+nspcc,row+nspcc] = (
                    + poro[iz]*(1e0)/dt 
                    - ((poro[iz]*dif_dic[iz]+poro[iz+1]*dif_dic[iz+1])*0.5e0*(-1e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                    - poro[iz]*dif_dic[iz]*(1e0)/dz[iz])/dz[iz]
                    - (1e0-poro[iz])*np.sum(drcc_ddic[iz,:])  
                    )*dicx[iz]*fact
                # put derivative of f(x) for dic at iz wrt dic at iz+1 in ln 
                amx[row+nspcc,row+nspcc+nsp] = (
                    - ((poro[iz]*dif_dic[iz]+poro[iz+1]*dif_dic[iz+1])*0.5e0*(1e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                    - 0e0)/dz[iz]
                    )*dicx[iz+1]*fact
                # put derivative of f(x) for dic at iz wrt alk at iz in ln 
                amx[row+nspcc,row+nspcc+1] = ( 
                    - (1e0-poro[iz])*np.sum(drcc_dalk[iz,:])  
                    )*alkx[iz]*fact
                # ALK
                # put f(x) for alk at iz  
                ymx[row+nspcc+1] = (
                    + poro[iz]*(alkx[iz]-alk[iz])/dt 
                    - ((poro[iz]*dif_alk[iz]+poro[iz+1]*dif_alk[iz+1])*0.5e0*(alkx[iz+1]-alkx[iz])/(0.5e0*(dz[iz]+dz[iz+1])) 
                    - poro[iz]*dif_alk[iz]*(alkx[iz]-alki*1e-6/1e3)/dz[iz])/dz[iz] 
                    - anco2[iz] 
                    - 2e0* (1e0-poro[iz])*np.sum(rcc[iz,:])  
                    )*fact
                # put derivative of f(x) for alk at iz wrt alk at iz in ln 
                amx[row+nspcc+1,row+nspcc+1] = (
                    + poro[iz]*(1e0)/dt 
                    - ((poro[iz]*dif_alk[iz]+poro[iz+1]*dif_alk[iz+1])*0.5e0*(-1e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                    - poro[iz]*dif_alk[iz]*(1e0)/dz[iz])/dz[iz]  
                    - 2e0* (1e0-poro[iz])*np.sum(drcc_dalk[iz,:])  
                    )*alkx[iz]*fact
                # put derivative of f(x) for alk at iz wrt alk at iz+1 in ln 
                amx[row+nspcc+1,row+nspcc+1+nsp] = (
                    - ((poro[iz]*dif_alk[iz]+poro[iz+1]*dif_alk[iz+1])*0.5e0*(1e0)/(0.5e0*(dz[iz]+dz[iz+1])) 
                    - 0e0)/dz[iz]
                    )*alkx[iz+1]*fact
                # put derivative of f(x) for alk at iz wrt dic at iz in ln 
                amx[row+nspcc+1,row+nspcc] = (
                    - 2e0* (1e0-poro[iz])*np.sum(drcc_ddic[iz,:])  
                    )*dicx[iz]*fact
            elif iz == nz-1: # need be careful about lower boundary condition; no diffusive flux from the bottom  
                for isp in range(nspcc):
                    ymx[row+isp] = \
                        + sporo[iz]*(ccx[iz,isp]-cc[iz,isp])/dt \
                        + adf[iz]*up[iz]*(sporo[iz]*w[iz]*ccx[iz,isp]-sporo[iz-1]*w[iz-1]*ccx[iz-1,isp])/dz[iz]  \
                        + adf[iz]*cnr[iz]*(sporof*w[iz]*ccx[iz,isp]-sporo[iz-1]*w[iz-1]*ccx[iz-1,isp])/dz[iz]  \
                        + adf[iz]*dwn[iz]*(sporof*w[iz]*ccx[iz,isp]-sporo[iz]*w[iz]*ccx[iz,isp])/dz[iz]  \
                        + sporo[iz]*rcc[iz,isp]
                    amx[row+isp,row+isp] = (
                        + sporo[iz]*(1e0)/dt 
                        + adf[iz]*up[iz]*(sporo[iz]*w[iz]*1e0-0e0)/dz[iz]  
                        + adf[iz]*cnr[iz]*(sporof*w[iz]*1e0-0e0)/dz[iz]  
                        + adf[iz]*dwn[iz]*(sporof*w[iz]*1e0-sporo[iz]*w[iz]*1e0)/dz[iz]  
                        + sporo[iz]*drcc_dcc[iz,isp]   
                        )*ccx[iz,isp]
                    amx[row+isp,row+isp-nsp] = ( 
                        + adf[iz]*up[iz]*(0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                        + adf[iz]*cnr[iz]*(0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                        )*ccx[iz-1,isp]
                    amx[row+isp,row+nspcc] = (
                        + sporo[iz]*drcc_ddic[iz,isp] 
                        )*dicx[iz]
                    amx[row+isp,row+nspcc+1] = (
                        + sporo[iz]*drcc_dalk[iz,isp] 
                        )*alkx[iz]
                    #DIC 
                    amx[row+nspcc,row+isp] = (
                        - sporo[iz]*drcc_dcc[iz,isp]  
                        )*ccx[iz,isp]*fact
                    #ALK 
                    amx[row+nspcc+1,row+isp] = (
                        - 2e0*sporo[iz]*drcc_dcc[iz,isp]  
                        )*ccx[iz,isp]*fact
                # DIC
                ymx[row+nspcc] = (
                    + poro[iz]*(dicx[iz]-dic[iz])/dt 
                    - (0e0 - 0.5e0*(poro[iz]*dif_dic[iz]+poro[iz-1]*dif_dic[iz-1])*(dicx[iz]-dicx[iz-1])  
                        /(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] 
                    - oxco2[iz] 
                    - anco2[iz] 
                    - sporo[iz]*np.sum(rcc[iz,:])  
                    )*fact
                amx[row+nspcc,row+nspcc] = ( 
                    + poro[iz]*(1e0)/dt 
                    - (0e0 - 0.5e0*(poro[iz]*dif_dic[iz]+poro[iz-1]*dif_dic[iz-1])*(1e0)/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] 
                    - sporo[iz]*np.sum(drcc_ddic[iz,:])  
                    )*dicx[iz]*fact
                amx[row+nspcc,row+nspcc-nsp] = ( 
                    - (0e0 - 0.5e0*(poro[iz]*dif_dic[iz]+poro[iz-1]*dif_dic[iz-1])*(-1e0)/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] 
                    ) * dicx[iz-1]*fact
                amx[row+nspcc,row+nspcc+1] = (
                    - sporo[iz]*np.sum(drcc_dalk[iz,:])  
                    )*alkx[iz]*fact
                # ALK 
                ymx[row+nspcc+1] = ( 
                    + poro[iz]*(alkx[iz]-alk[iz])/dt 
                    - (0e0 - 0.5e0*(poro[iz]*dif_alk[iz]+poro[iz-1]*dif_alk[iz-1])*(alkx[iz]-alkx[iz-1])/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] 
                    - anco2[iz] 
                    - 2e0*sporo[iz]*np.sum(rcc[iz,:])  
                    )*fact
                amx[row+nspcc+1,row+nspcc+1] = ( 
                    + poro[iz]*(1e0)/dt 
                    - (0e0 - 0.5e0*(poro[iz]*dif_alk[iz]+poro[iz-1]*dif_alk[iz-1])*(1e0)/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] 
                    - 2e0*sporo[iz]*np.sum(drcc_dalk[iz,:])  
                    )*alkx[iz]*fact
                amx[row+nspcc+1,row+nspcc+1-nsp] = ( 
                    - (0e0 - 0.5e0*(poro[iz]*dif_alk[iz]+poro[iz-1]*dif_alk[iz-1])*(-1e0)/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz] 
                    ) * alkx[iz-1]*fact
                amx[row+nspcc+1,row+nspcc] = (
                    - 2e0*sporo[iz]*np.sum(drcc_ddic[iz,:])  
                    )*dicx[iz]*fact
            else: #  do not have to be careful abount boundary conditions 
                for isp in range(nspcc):
                    ymx[row+isp] = \
                        + sporo[iz]*(ccx[iz,isp]-cc[iz,isp])/dt \
                        + adf[iz]*up[iz]*(sporo[iz]*w[iz]*ccx[iz,isp]-sporo[iz-1]*w[iz-1]*ccx[iz-1,isp])/dz[iz]  \
                        + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*ccx[iz+1,isp]-sporo[iz]*w[iz]*ccx[iz,isp])/dz[iz]  \
                        + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*ccx[iz+1,isp]-sporo[iz-1]*w[iz-1]*ccx[iz-1,isp])/dz[iz]  \
                        + sporo[iz]*rcc[iz,isp]
                    amx[row+isp,row+isp] = (
                        + sporo[iz]*(1e0)/dt 
                        + adf[iz]*up[iz]*(sporo[iz]*w[iz]*1e0-0e0)/dz[iz]  
                        + adf[iz]*dwn[iz]*(0e0-sporo[iz]*w[iz]*1e0)/dz[iz]  
                        + sporo[iz]*drcc_dcc[iz,isp]  
                        )*ccx[iz,isp]
                    amx[row+isp,row+isp+nsp] =  (
                        + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*1e0-0e0)/dz[iz]  
                        + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*1e0-0e0)/dz[iz]  
                        )*ccx[iz+1,isp]
                    amx[row+isp,row+isp-nsp] =  (
                        + adf[iz]*up[iz]*(0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                        + adf[iz]*cnr[iz]*(0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                        )*ccx[iz-1,isp]
                    amx[row+isp,row+nspcc] = (
                        + sporo[iz]*drcc_ddic[iz,isp]  
                        )*dicx[iz]
                    amx[row+isp,row+nspcc+1] = (
                        + sporo[iz]*drcc_dalk[iz,isp] 
                        )*alkx[iz]
                    # DIC 
                    amx[row+nspcc,row+isp] = (
                        - sporo[iz]*drcc_dcc[iz,isp]  
                        )*ccx[iz,isp]*fact
                    # ALK 
                    amx[row+nspcc+1,row+isp] = (
                        - 2e0*sporo[iz]*drcc_dcc[iz,isp]  
                        )*ccx[iz,isp]*fact 
                # DIC 
                ymx[row+nspcc] = ( 
                    + poro[iz]*(dicx[iz]-dic[iz])/dt 
                    - (0.5e0*(poro[iz+1]*dif_dic[iz+1]+poro[iz]*dif_dic[iz])*(dicx[iz+1]-dicx[iz])/(0.5e0*(dz[iz+1]+dz[iz])) 
                    - 0.5e0*(poro[iz]*dif_dic[iz]+poro[iz-1]*dif_dic[iz-1])*(dicx[iz]-dicx[iz-1])/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]  
                    - oxco2[iz] 
                    - anco2[iz] 
                    - sporo[iz]*np.sum(rcc[iz,:])  
                    )*fact
                amx[row+nspcc,row+nspcc] = (
                    + poro[iz]*(1e0)/dt 
                    - (0.5e0*(poro[iz+1]*dif_dic[iz+1]+poro[iz]*dif_dic[iz])*(-1e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                    - 0.5e0*(poro[iz]*dif_dic[iz]+poro[iz-1]*dif_dic[iz-1])*(1e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]  
                    - sporo[iz]*np.sum(drcc_ddic[iz,:])  
                    )*dicx[iz]*fact
                amx[row+nspcc,row+nspcc+nsp] = (
                    - (0.5e0*(poro[iz+1]*dif_dic[iz+1]+poro[iz]*dif_dic[iz])*(1e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                    - 0e0)/dz[iz]  
                    )*dicx[iz+1]*fact
                amx[row+nspcc,row+nspcc-nsp] = (
                    - (0e0 
                    - 0.5e0*(poro[iz]*dif_dic[iz]+poro[iz-1]*dif_dic[iz-1])*(-1e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz] 
                    )*dicx[iz-1]*fact
                amx[row+nspcc,row+nspcc+1] = (
                    - sporo[iz]*np.sum(drcc_dalk[iz,:])  
                    )*alkx[iz]*fact
                # ALK 
                ymx[row+nspcc+1] = (
                    + poro[iz]*(alkx[iz]-alk[iz])/dt 
                    - (0.5e0*(poro[iz+1]*dif_alk[iz+1]+poro[iz]*dif_alk[iz])*(alkx[iz+1]-alkx[iz])/(0.5e0*(dz[iz+1]+dz[iz])) 
                    - 0.5e0*(poro[iz]*dif_alk[iz]+poro[iz-1]*dif_alk[iz-1])*(alkx[iz]-alkx[iz-1])/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz] 
                    - anco2[iz] 
                    - 2e0*sporo[iz]*np.sum(rcc[iz,:])  
                    ) *fact
                amx[row+nspcc+1,row+nspcc+1] = (
                    + poro[iz]*(1e0)/dt 
                    - (0.5e0*(poro[iz+1]*dif_alk[iz+1]+poro[iz]*dif_alk[iz])*(-1e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                    - 0.5e0*(poro[iz]*dif_alk[iz]+poro[iz-1]*dif_alk[iz-1])*(1e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]  
                    - 2e0*sporo[iz]*np.sum(drcc_dalk[iz,:])  
                    )*alkx[iz]*fact
                amx[row+nspcc+1,row+nspcc+1+nsp] = ( 
                    - (0.5e0*(poro[iz+1]*dif_alk[iz+1]+poro[iz]*dif_alk[iz])*(1e0)/(0.5e0*(dz[iz+1]+dz[iz])) 
                    - 0e0)/dz[iz]  
                    )*alkx[iz+1]*fact
                amx[row+nspcc+1,row+nspcc+1-nsp] = (
                    - (0e0 
                    - 0.5e0*(poro[iz]*dif_alk[iz]+poro[iz-1]*dif_alk[iz-1])*(-1e0)/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz] 
                    )*alkx[iz-1]*fact
                amx[row+nspcc+1,row+nspcc] = (
                    - 2e0*sporo[iz]*np.sum(drcc_ddic[iz,:])  
                    )*dicx[iz]*fact 
            # diffusion terms are filled with transition matrices 
            for isp in range(nspcc):
                if turbo2[isp+2] or labs[isp+2]: 
                    for iiz in range(nz):
                        col = (iiz)*nsp
                        if trans[iiz,iz,isp+2]==0e0: continue
                        amx[row+isp,col+isp] = amx[row+isp,col+isp] \
                            - trans[iiz,iz,isp+2]/dz[iz]*dz[iiz]*(1e0-poro[iiz])*ccx[iiz,isp]
                        ymx[row+isp] = ymx[row+isp] \
                            - trans[iiz,iz,isp+2]/dz[iz]*dz[iiz]*(1e0-poro[iiz])*ccx[iiz,isp]
                else:
                    for iiz in range( nz):
                        col = (iiz)*nsp
                        if trans[iiz,iz,isp+2]==0e0: continue
                        amx[row+isp,col+isp] = amx[row+isp,col+isp] -trans[iiz,iz,isp+2]/dz[iz]*ccx[iiz,isp]
                        ymx[row+isp] = ymx[row+isp] - trans[iiz,iz,isp+2]/dz[iz]*ccx[iiz,isp]
        ymx = - ymx  # because I put f(x) into ymx (=B), minus sign need be added 
        if sparse:
            amx=amx.tocsr()
            kai=spsolve(amx,ymx)
        else:
            # kai = np.linalg.solve(amx, ymx)
            # kai = la.solve(amx, ymx) # solving
            lu, piv, kai, info = lapack.dgesv(amx, ymx)
        # if itr ==1:time.sleep(5)
        ymx[:] = kai[:].copy()  
        if any(np.isnan(ymx)):
            print 'NAN in ymx'
            dt=dt/10e0
            dicx[:]=dic[:].copy()
            alkx[:]=alk[:].copy()
            ccx[:,:]=cc[:,:].copy()
            w[:]=w_pre[:].copy()
            # upwind(w)
            up,dwn,cnr,adf = calcupwindscheme(w,nz)
            print 'location 3'
            flg_restart = True 
            break 
            file_tmp = open(workdir+'chk_ymx_pre.txt','w')
            for iz in range( nmx):
                print>>file_tmp, ymx[iz]
            file_tmp.close()
            input()
        for iz in range( nz ):
            row = (iz)*nsp
            # print iz,row, ymx[row:row+nsp]
            for isp in range(nspcc):
                if ymx[row+isp]>10e0: # this help conversion 
                    ccx[iz,isp] *= 1.5e0
                    if np.isinf(ccx[iz,isp]): 
                        print iz,'too large ccx'
                        input()
                elif ymx[row+isp]<-10e0: # this help conversion  
                    ccx[iz,isp] *= 0.5e0
                    if np.isinf(ccx[iz,isp]): 
                        print iz,'too large ccx'
                        input()
                else:
                    # print iz,isp,ymx[row+isp],np.exp(ymx[row+isp])
                    ccx[iz,isp] *= np.exp(ymx[row+isp])
                    if np.isinf(ccx[iz,isp]): 
                        print iz,'too large ccx'
                        input()
                if ccx[iz,isp]<ccx_th: # too small trancate value and not be accounted for error 
                    ccx[iz,isp]=ccx_th
                    ymx[row+isp] = 0e0
            if ymx[row+nspcc]>10e0: 
                dicx[iz]*=1.5e0
            elif ymx[row+nspcc]<-10e0: 
                dicx[iz]*=0.5e0
            else :
                # print ymx[row+nspcc],np.exp(ymx[row+nspcc])
                dicx[iz] *=np.exp(ymx[row+nspcc])
            if ymx[row+nspcc+1]>10e0:
                alkx[iz] *=1.5e0
            elif ymx[row+nspcc+1]<-10e0: 
                alkx[iz] *=0.5e0
            else :
                # print ymx[row+nspcc+1],np.exp(ymx[row+nspcc+1])
                alkx[iz] *=np.exp(ymx[row+nspcc+1])
            if dicx[iz]<1e-100: ymx[row+nspcc] = 0e0
            if alkx[iz]<1e-100: ymx[row+nspcc+1] = 0e0
        error = np.max(np.exp(np.abs(ymx[:]))) - 1e0
        itr = itr + 1
        if itr > 5000:
            dt=dt/10e0
            dicx[:]=dic[:].copy()
            alkx[:]=alk[:].copy()
            ccx[:,:]=cc[:,:].copy()
            w[:]=w_pre[:].copy()
            # upwind(w)
            up,dwn,cnr,adf = calcupwindscheme(w,nz)
            print '... unable to converge'
            flg_restart = True 
            continue
            print 'stop'
        if showiter:
            print 'co2 iteration',itr,error
            print  'cc :',np.sum(ccx[0:nz:nz/5,:]*mcc[:],axis=1)/rho[0:nz:nz/5]*100e0
            print  'dic:',dicx[0:nz:nz/5]#*1e3
            print  'alk:',alkx[0:nz:nz/5]#*1e3
            print  '   ..... multiple cc species ..... '
            for isp in range(nspcc ):
                print '{:03}'.format(isp+1),ccx[0:nz:nz/5,isp]*mcc[isp]/rho[0:nz:nz/5]*100e0
        if (ccx<0e0).any():
            print 'negative ccx, stop'
            print ccx
            input()
        if (np.isnan(ccx)).any():
            print 'nan om, stop'
            print ccx
            input()
        if (dicx<0e0).any():
            print 'negative dicx, stop'
            print dicx
            input()
        if (np.isnan(dicx)).any():
            print 'nan dic, stop'
            print dicx
            input()
        if (alkx<0e0).any():
            print 'negative alk, stop'
            print alkx
        if (np.isnan(alkx)).any():
            print 'nan alk, stop'
            print alkx
            input()
    return ccx,dicx,alkx,rcc,dt,flg_restart,w
    
def calcflxcaco3sys(  
    dw # inoutput
     ,nspcc,ccx,cc,dt,dz,rcc,adf,up,dwn,cnr,w,dif_alk,dif_dic,dic,dicx,alk,alkx,oxco2,anco2,trans    # input
     ,turbo2,labs,nonlocal,sporof,it,nz,poro,sporo,ccflx,dici,alki,mvcc,tol,workdir        # input
     ): 
    cctflx =np.zeros((nspcc),dtype=np.float64)
    ccdis = np.zeros((nspcc),dtype=np.float64) 
    ccdif = np.zeros((nspcc),dtype=np.float64)
    ccadv = np.zeros((nspcc),dtype=np.float64)
    ccrain = np.zeros((nspcc),dtype=np.float64)
    ccres = np.zeros((nspcc),dtype=np.float64)
    dictflx = np.float64(0e0) 
    dicdis = np.float64(0e0) 
    dicdif = np.float64(0e0) 
    dicdec = np.float64(0e0) 
    dicres = np.float64(0e0)
    alktflx = np.float64(0e0) 
    alkdis = np.float64(0e0) 
    alkdif = np.float64(0e0) 
    alkdec = np.float64(0e0) 
    alkres = np.float64(0e0)
    for iz in range(nz):
        if iz == 0: 
            for isp in range(nspcc):
                cctflx[isp] = cctflx[isp] + (1e0-poro[iz])*(ccx[iz,isp]-cc[iz,isp])/dt *dz[iz]
                ccdis[isp] = ccdis[isp]  + (1e0-poro[iz])*rcc[iz,isp] *dz[iz]
                ccrain[isp] = ccrain[isp] - ccflx[isp]/dz[iz]*dz[iz]
                ccadv[isp] = ccadv[isp] + adf[iz]*up[iz]*(sporo[iz]*w[iz]*ccx[iz,isp]-0e0)/dz[iz] * dz[iz] \
                    + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*ccx[iz+1,isp]-sporo[iz]*w[iz]*ccx[iz,isp])/dz[iz] * dz[iz]  \
                    + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*ccx[iz+1,isp]-0e0)/dz[iz] * dz[iz]
            #  DIC 
            dictflx = dictflx +(dicx[iz]-dic[iz])/dt*dz[iz]*poro[iz] 
            dicdif = dicdif - ((poro[iz]*dif_dic[iz]+poro[iz+1]*dif_dic[iz+1])*0.5e0*(dicx[iz+1]-dicx[iz])/(0.5e0*(dz[iz]+dz[iz+1])) \
                - poro[iz]*dif_dic[iz]*(dicx[iz]-dici*1e-6/1e3)/dz[iz])/dz[iz]*dz[iz]
            dicdec = dicdec - oxco2[iz]*dz[iz] - anco2[iz]*dz[iz] 
            dicdis = dicdis - np.sum(rcc[iz,:])*sporo[iz]*dz[iz] 
            # ALK
            alktflx = alktflx + (alkx[iz]-alk[iz])/dt*dz[iz]*poro[iz]
            alkdif = alkdif - ((poro[iz]*dif_alk[iz]+poro[iz+1]*dif_alk[iz+1])*0.5e0*(alkx[iz+1]-alkx[iz])/(0.5e0*(dz[iz]+dz[iz+1])) \
                - poro[iz]*dif_alk[iz]*(alkx[iz]-alki*1e-6/1e3)/dz[iz])/dz[iz]*dz[iz]
            alkdec = alkdec - anco2[iz]*dz[iz] 
            alkdis = alkdis - 2e0* sporo[iz]*np.sum(rcc[iz,:])*dz[iz] 
        elif iz == nz-1: 
            for isp in range(nspcc):
                cctflx[isp] = cctflx[isp] + sporo[iz]*(ccx[iz,isp]-cc[iz,isp])/dt *dz[iz]
                ccdis[isp] = ccdis[isp]  + sporo[iz]*rcc[iz,isp] *dz[iz]
                ccadv[isp] = ccadv[isp] \
                    + adf[iz]*up[iz]*(sporo[iz]*w[iz]*ccx[iz,isp]-sporo[iz-1]*w[iz-1]*ccx[iz-1,isp])/dz[iz] * dz[iz]  \
                    + adf[iz]*cnr[iz]*(sporof*w[iz]*ccx[iz,isp]-sporo[iz-1]*w[iz-1]*ccx[iz-1,isp])/dz[iz] * dz[iz]  \
                    + adf[iz]*dwn[iz]*(sporof*w[iz]*ccx[iz,isp]-sporo[iz]*w[iz]*ccx[iz,isp])/dz[iz] * dz[iz]  
            # DIC
            dictflx = dictflx +(dicx[iz]-dic[iz])/dt*dz[iz]*poro[iz] 
            dicdif = dicdif - (0e0 
                - 0.5e0*(poro[iz]*dif_dic[iz]+poro[iz-1]*dif_dic[iz-1])*(dicx[iz]-dicx[iz-1])/(0.5e0*(dz[iz-1]+dz[iz])) 
                )/dz[iz]*dz[iz]
            dicdec = dicdec - oxco2[iz]*dz[iz] - anco2[iz]*dz[iz] 
            dicdis = dicdis - sporo[iz]*np.sum(rcc[iz,:])*dz[iz] 
            # ALK 
            alktflx = alktflx + (alkx[iz]-alk[iz])/dt*dz[iz]*poro[iz]
            alkdif = alkdif - (0e0 
                - 0.5e0*(poro[iz]*dif_alk[iz]+poro[iz-1]*dif_alk[iz-1])*(alkx[iz]-alkx[iz-1])/(0.5e0*(dz[iz-1]+dz[iz])))/dz[iz]*dz[iz]
            alkdec = alkdec - anco2[iz]*dz[iz]
            alkdis = alkdis - 2e0* sporo[iz]*np.sum(rcc[iz,:])*dz[iz]
        else :
            for isp in range(nspcc):
                cctflx[isp] = cctflx[isp] + sporo[iz]*(ccx[iz,isp]-cc[iz,isp])/dt *dz[iz]
                ccdis[isp] = ccdis[isp]  + sporo[iz]*rcc[iz,isp] *dz[iz]
                ccadv[isp] = ccadv[isp] \
                    + adf[iz]*up[iz]*(sporo[iz]*w[iz]*ccx[iz,isp]-sporo[iz-1]*w[iz-1]*ccx[iz-1,isp])/dz[iz] * dz[iz]  \
                    + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*ccx[iz+1,isp]-sporo[iz]*w[iz]*ccx[iz,isp])/dz[iz] * dz[iz]  \
                    + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*ccx[iz+1,isp]-sporo[iz-1]*w[iz-1]*ccx[iz-1,isp])/dz[iz] * dz[iz]  
            # DIC 
            dictflx = dictflx +(dicx[iz]-dic[iz])/dt*dz[iz]*poro[iz] 
            dicdif = dicdif - (0.5e0*(poro[iz+1]*dif_dic[iz+1]+poro[iz]*dif_dic[iz])*(dicx[iz+1]-dicx[iz])/(0.5e0*(dz[iz+1]+dz[iz])) 
                - 0.5e0*(poro[iz]*dif_dic[iz]+poro[iz-1]*dif_dic[iz-1])*(dicx[iz]-dicx[iz-1])/(0.5e0*(dz[iz]+dz[iz-1])) 
                )/dz[iz]*dz[iz]
            dicdec = dicdec - oxco2[iz]*dz[iz] - anco2[iz]*dz[iz] 
            dicdis = dicdis - sporo[iz]*np.sum(rcc[iz,:])*dz[iz] 
            # ALK 
            alktflx = alktflx + (alkx[iz]-alk[iz])/dt*dz[iz]*poro[iz]
            alkdif = alkdif - (0.5e0*(poro[iz+1]*dif_alk[iz+1]+poro[iz]*dif_alk[iz])*(alkx[iz+1]-alkx[iz])/(0.5e0*(dz[iz+1]+dz[iz])) 
                - 0.5e0*(poro[iz]*dif_alk[iz]+poro[iz-1]*dif_alk[iz-1])*(alkx[iz]-alkx[iz-1])/(0.5e0*(dz[iz]+dz[iz-1])))/dz[iz]*dz[iz]
            alkdec = alkdec - anco2[iz]*dz[iz]
            alkdis = alkdis - 2e0* sporo[iz]*np.sum(rcc[iz,:])*dz[iz]
        for isp in range(nspcc):
            if labs[isp+2] or  turbo2[isp+2]: 
                for iiz in range( nz):
                    if trans[iiz,iz,isp+2]==0e0:continue
                    ccdif[isp] = ccdif[isp] -trans[iiz,iz,isp+2]/dz[iz]*dz[iiz]*(1e0-poro[iiz])*dz[iz]*ccx[iiz,isp]
            else :
                for iiz in range(nz):
                    if trans[iiz,iz,isp+2]==0e0: continue
                    ccdif[isp] = ccdif[isp] -trans[iiz,iz,isp+2]/dz[iz]*dz[iz]*ccx[iiz,isp]
            if labs[isp+2] or  turbo2[isp+2]: 
                for iiz in range(nz):
                    if trans[iiz,iz,isp+2]==0e0: continue
                    dw[iz] = dw[iz] - mvcc[isp]*(-trans[iiz,iz,isp+2]/dz[iz]*dz[iiz]*(1e0-poro[iiz])*ccx[iiz,isp])
            else :
                if nonlocal[isp+2]:
                    for iiz in range( nz):
                        if trans[iiz,iz,isp+2]==0e0: continue
                        dw[iz] = dw[iz] -mvcc[isp]*(-trans[iiz,iz,isp+2]/dz[iz]*ccx[iiz,isp])
        dw[iz] = dw[iz] -(1e0-poro[iz])*mvcc[isp]*np.sum(rcc[iz,:])
    # residual fluxes 
    ccres[:] = cctflx[:] +  ccdis[:] +  ccdif[:] + ccadv[:] + ccrain[:]
    dicres = dictflx + dicdis + dicdif + dicdec 
    alkres = alktflx + alkdis + alkdif + alkdec 
    if abs(alkres)/np.max([abs(alktflx),abs(alkdis) ,abs(alkdif) , abs(alkdec)]) > tol*10e0:   
    # if residula fluxes are relatively large, record just in case  
        print 'not enough accuracy in co2 calc:stop',abs(alkres)/np.max([abs(alktflx),abs(alkdis) ,abs(alkdif) , abs(alkdec)])
        file_err=open(workdir+'errlog.txt','a')
        print >> file_err, 'not enough accuracy in co2 calc:stop',abs(alkres)/np.max([abs(alktflx),abs(alkdis) ,abs(alkdif) , abs(alkdec)])
    return cctflx,ccflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres \
        ,dictflx,dicdis,dicdif,dicres,dicdec  \
        ,dw
        
def claycalc( 
    ptx
    ,nz,sporo,pt,dt,w,dz,detflx,adf,up,dwn,cnr,trans  # input
    ,nspcc,labs,turbo2,nonlocal,poro,sporof,msed,workdir     # intput
    ):
    nsp = 1 #  only consider clay
    nmx = nz*nsp  # matrix is linear and solved like om and o2, so see comments there for calculation procedures 
    amx = np.zeros((nmx,nmx),dtype=np.float64)
    ymx = np.zeros((nmx),dtype=np.float64)
    for iz in range(nz ):
        row = (iz)*nsp 
        if iz == 0: 
            ymx[row] = \
                + sporo[iz]*(-pt[iz])/dt \
                - detflx/msed/dz[iz]
            amx[row,row] = (
                + sporo[iz]*(1e0)/dt 
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*1e0-0e0)/dz[iz]   
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*0e0-sporo[iz]*w[iz]*1e0)/dz[iz]   
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*0e0-0e0)/dz[iz]   
                )            
            amx[row,row+nsp] =  (
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*1e0-sporo[iz]*w[iz]*0e0)/dz[iz]   
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*1e0-0e0)/dz[iz]   
                )
        elif iz == nz-1:
            ymx[row] = \
                + sporo[iz]*(-pt[iz])/dt 
            amx[row,row] = (
                + sporo[iz]*(1e0)/dt 
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*1e0-0e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(sporof*w[iz]*1e0-0e0)/dz[iz]  
                + adf[iz]*dwn[iz]*(sporof*w[iz]*1e0-sporo[iz]*w[iz]*1e0)/dz[iz]  
                )
            amx[row,row-nsp] = ( 
                + adf[iz]*up[iz]*(0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                )
        else :
            ymx[row] = \
                + sporo[iz]*(-pt[iz])/dt 
            amx[row,row] = (
                + sporo[iz]*(1e0)/dt 
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*1e0-0e0)/dz[iz]  
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*0e0-sporo[iz]*w[iz]*1e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*0e0-0e0)/dz[iz]  
                )
            amx[row,row+nsp] =  (
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*1e0-sporo[iz]*w[iz]*0e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*1e0-0e0)/dz[iz]  
                )
            amx[row,row-nsp] =  (
                + adf[iz]*up[iz]*(0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                + adf[iz]*cnr[iz]*(0e0-sporo[iz-1]*w[iz-1]*1e0)/dz[iz]  
                )
        if labs[1] or turbo2[1]:
            for iiz in range(nz):
                col = (iiz)*nsp
                if trans[iiz,iz,1]==0e0: continue
                amx[row,col] = amx[row,col] -trans[iiz,iz,1]/dz[iz]*dz[iiz]*(1e0-poro[iiz])
        else :
            for iiz in range(nz):
                col = (iiz)*nsp
                if trans[iiz,iz,1]==0e0: continue
                amx[row,col] = amx[row,col] -trans[iiz,iz,1]/dz[iz]
    ymx = - ymx
    if any(np.isnan(ymx)):
        print 'NAN in ymx:pt'
        file_tmp=open(workdir+ 'chk_ymx_pre_pt.txt', 'w')
        for iz in range(nmx):
            print>> file_tmp,ymx[iz]
        file_tmp.close()
        input()
    # kai = np.linalg.solve(amx, ymx)
    # kai = la.solve(amx, ymx) # solving
    lu, piv, kai, info = lapack.dgesv(amx, ymx)
    ymx[:] = kai[:].copy()
    if np.isnan(amx).any():
        print 'NAN in amx:pt'
        file_tmp = open(workdir+'chk_amx_pt.txt', 'w')
        for iz in range( nmx):
            print>> file_tmp,amx[iz,:]
        file_tmp.close()
        input()
    if np.isnan(ymx).any(): 
        print 'NAN in ymx:pt'
        file_tmp = open(workdir+'chk_ymx_pt.txt', 'w')
        for iz in range( nmx):
            print>>file_tmp,ymx[iz]
        file_tmp.close()
    ptx[:] = ymx[:].copy()
    return ptx

def calcflxclay( 
    dw         # in&output
    ,nz,sporo,ptx,pt,dt,dz,detflx,w,adf,up,dwn,cnr,sporof,trans,nspcc,turbo2,labs,nonlocal,poro          #  input
    ,msed,mvsed # input
    ):
    pttflx = 0e0 
    ptdif = 0e0 
    ptadv = 0e0 
    ptres = 0e0
    ptrain = 0e0
    for  iz in range(nz ):
        if iz == 0:
            pttflx = pttflx + sporo[iz]*(ptx[iz]-pt[iz])/dt*dz[iz]
            ptrain = ptrain - detflx/msed
            ptadv = ptadv \
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*ptx[iz]-0e0)/dz[iz]*dz[iz]  \
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*ptx[iz+1]-sporo[iz]*w[iz]*ptx[iz])/dz[iz]*dz[iz]  \
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*ptx[iz+1]-0e0)/dz[iz]*dz[iz]  
        elif iz == nz-1:
            pttflx = pttflx + (1e0-poro[iz])*(ptx[iz]-pt[iz])/dt*dz[iz]
            ptadv = ptadv \
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*ptx[iz]-sporo[iz-1]*w[iz-1]*ptx[iz-1])/dz[iz]*dz[iz]  \
                + adf[iz]*cnr[iz]*(sporof*w[iz]*ptx[iz]-sporo[iz-1]*w[iz-1]*ptx[iz-1])/dz[iz]*dz[iz]  \
                + adf[iz]*dwn[iz]*(sporof*w[iz]*ptx[iz]-sporo[iz]*w[iz]*ptx[iz])/dz[iz]*dz[iz]  
        else :
            pttflx = pttflx + (1e0-poro[iz])*(ptx[iz]-pt[iz])/dt*dz[iz]
            ptadv = ptadv \
                + adf[iz]*up[iz]*(sporo[iz]*w[iz]*ptx[iz]-sporo[iz-1]*w[iz-1]*ptx[iz-1])/dz[iz]*dz[iz]  \
                + adf[iz]*dwn[iz]*(sporo[iz+1]*w[iz+1]*ptx[iz+1]-sporo[iz]*w[iz]*ptx[iz])/dz[iz]*dz[iz]  \
                + adf[iz]*cnr[iz]*(sporo[iz+1]*w[iz+1]*ptx[iz+1]-sporo[iz-1]*w[iz-1]*ptx[iz-1])/dz[iz]*dz[iz]
        if turbo2[1] or labs[1]:
            for iiz in range( nz):
                if trans[iiz,iz,1]==0e0: continue
                ptdif = ptdif -trans[iiz,iz,1]*ptx[iiz]/dz[iz]*dz[iiz]*dz[iz]
        else :
            for iiz in range( nz):
                if trans[iiz,iz,1]==0e0: continue
                ptdif = ptdif -trans[iiz,iz,1]*ptx[iiz]/dz[iz]    \
                    *dz[iz]
        if turbo2[1] or labs[1]:
            for iiz in range( nz):
                if trans[iiz,iz,1]==0e0: continue
                dw[iz] = dw[iz] - mvsed*(-trans[iiz,iz,1]*ptx[iiz]/dz[iz]*dz[iiz]*(1e0-poro[iiz]))
        else :
            if nonlocal[1]: 
                for iiz in range( nz):
                    if trans[iiz,iz,1]==0e0: continue
                    dw[iz] = dw[iz] - mvsed*(-trans[iiz,iz,1]*ptx[iiz]/dz[iz])
    ptres = pttflx + ptdif + ptadv + ptrain
    return pttflx,ptdif,ptadv,ptres,ptrain

def make_transmx(  
    labs,nspcc,turbo2,nobio,dz,sporo,nz,z,zmlref,size  # input
    ):
    nonlocal = np.zeros(nspcc+2,dtype=bool)
    translabs=np.zeros((nz,nz),dtype=np.float64)
    if any(labs):
        translabs = np.loadtxt('./labs-mtx.txt')
    if True:
        translabs *= 365.25/10.*1./3.
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    zml = np.zeros(nspcc+2)
    zml[:]=zmlref
    zrec = 1.1*np.max(zml)
    zrec2 = 2.0*np.max(zml)
    if size:
        zml[2+0]=20.
        zml[2+1]=20.
        zml[2+2]=20.
        zml[2+3]=20.
        zrec = 1.1*np.min(zml)
        zrec2= 1.1*np.max(zml)
    for iz in range(nz):
        if z[iz]<=zrec:izrec = iz
        if z[iz]<=zrec2:izrec2 = iz

    nonlocal[:]=False
    trans = np.zeros((nz,nz,nspcc+2),dtype=np.float64)
    for isp in range(nspcc+2):
        if turbo2[isp] or labs[isp]:nonlocal[isp]=True
        dbio=np.zeros(nz,dtype=np.float64)
        for iz in range(nz):
            if z[iz]<=zml[isp]:
                dbio[iz]=0.15
                izml = iz
            else:
                dbio[iz]=0.
        transdbio = np.zeros((nz,nz),dtype=np.float64)
        for iz in range(0,izml+1):
            if iz==0:
                transdbio[iz,iz]=0.5*(sporo[iz]*dbio[iz]+sporo[iz+1]*dbio[iz+1])\
                                  *(-1.)/(0.5*(dz[iz]+dz[iz+1]))
                transdbio[iz+1,iz]=0.5*(sporo[iz]*dbio[iz]+sporo[iz+1]*dbio[iz+1])\
                                    *(1.)/(0.5*(dz[iz]+dz[iz+1]))
            elif iz==izml:
                transdbio[iz,iz]=0.5*(sporo[iz]*dbio[iz]+sporo[iz-1]*dbio[iz-1])\
                                  *(-1.)/(0.5*(dz[iz]+dz[iz-1]))
                transdbio[iz-1,iz]=0.5*(sporo[iz]*dbio[iz]+sporo[iz-1]*dbio[iz-1])\
                                    *(1.)/(0.5*(dz[iz]+dz[iz-1]))
            else:
                transdbio[iz,iz]=0.5*(sporo[iz]*dbio[iz]+sporo[iz-1]*dbio[iz-1])\
                                  *(-1.)/(0.5*(dz[iz]+dz[iz-1]))\
                                  + 0.5*(sporo[iz]*dbio[iz]+sporo[iz+1]*dbio[iz+1])\
                                  *(-1.)/(0.5*(dz[iz]+dz[iz+1]))
                transdbio[iz-1,iz]=0.5*(sporo[iz]*dbio[iz]+sporo[iz-1]*dbio[iz-1])\
                                    *(1.)/(0.5*(dz[iz]+dz[iz-1]))
                transdbio[iz+1,iz]=0.5*(sporo[iz]*dbio[iz]+sporo[iz+1]*dbio[iz+1])\
                                    *(1.)/(0.5*(dz[iz]+dz[iz+1]))
        transturbo2 = np.zeros((nz,nz),dtype=np.float64)
        probh = 0.001
        transturbo2[:izml+1,:izml+1] = probh
        for iz in range(0,izml+1):
            transturbo2[iz,iz]=-probh*izml

        if turbo2[isp]: translabs = transturbo2

        trans[:,:,isp] = transdbio[:,:].copy()

        if nonlocal[isp]: trans[:,:,isp] = translabs[:,:].copy()
        if nobio[isp]: trans[:,:,isp] = 0.

    if not any(nonlocal):
        for isp in range(nspcc+2-1):
            if not (trans[:,:,isp+1]==trans[:,:,isp]).all():nonlocal[:] = True
    return trans,izrec,izrec2,izml,nonlocal

def recordprofile(itrec 
    ,nz,z,age,pt,ptx,msed,w,wi,rho,frt 
    ,cc,mcc,dic,alk,co3,co3sat,rcc,pro,dicx,alkx,co3x,prox,nspcc
    ,o2x,oxco2,anco2
    ,om,omx,mom
    ,ccx
    ,d13c_ocni,d18o_ocni,d13c_blk,d18o_blk
    ,up,dwn,cnr,adf
    ,workdir
    ): 
    if itrec == 0: 
        file_tmp=open(workdir+'ptx-'+'{:03}'.format(0)+'.txt','w')
        for iz in range(nz):
            print >> file_tmp,z[iz],age[iz],pt[iz]*msed/2.5*100.,0.,1.,wi
        file_tmp.close()
        file_tmp=open(workdir+'ccx-'+'{:03}'.format(0)+'.txt','w')
        for iz in range(nz):
            print >> file_tmp,z[iz],age[iz],np.sum(cc[iz,:]*mcc[:])/2.5*100.,dic[iz]*1e3\
                  ,alk[iz]*1e3,co3[iz]*1e3-co3sat,np.sum(rcc[iz,:]),-np.log10(pro[iz])
        file_tmp.close()
        file_tmp=open(workdir+'o2x-'+'{:03}'.format(0)+'.txt','w')
        for iz in range(nz):
            print >> file_tmp,z[iz],age[iz],o2x[iz]*1e3,oxco2[iz],anco2[iz]
        file_tmp.close()
        file_tmp=open(workdir+'omx-'+'{:03}'.format(0)+'.txt','w')
        for iz in range(nz):
            print >> file_tmp,z[iz],age[iz],om[iz]*mom/2.5*100.
        file_tmp.close()
        file_tmp=open(workdir+'ccx_sp-'+'{:03}'.format(0)+'.txt','w')
        for iz in range(nz):
            # print >> file_tmp,z[iz],age[iz],ccx[iz,:]*mcc/2.5*100.
            print >> file_tmp,z[iz],age[iz],str(ccx[iz,:]*mcc[:]/2.5*100.)[1:-1]
        file_tmp.close()
        file_tmp=open(workdir+'sig-'+'{:03}'.format(0)+'.txt','w')
        for iz in range(nz):
            print >> file_tmp,z[iz],d13c_ocni,d18o_ocni
        file_tmp.close()
        file_tmp=open(workdir+'bur-'+'{:03}'.format(0)+'.txt','w')
        for iz in range(nz):
            print >> file_tmp,z[iz],age[iz],w[iz],up[iz],dwn[iz],cnr[iz],adf[iz]
        file_tmp.close()
    else:
        file_tmp = open(workdir+'ptx-'+ '{:03}'.format(itrec)+'.txt','w') 
        for iz in range(nz):
            print>>file_tmp,z[iz],age[iz],ptx[iz]*msed/rho[iz]*100e0,rho[iz],frt[iz]  ,w[iz]
        file_tmp.close()
        file_tmp = open(workdir+'ccx-'+ '{:03}'.format(itrec)+'.txt','w') 
        for iz in range(nz):
            print>>file_tmp, z[iz],age[iz],np.sum(ccx[iz,:]*mcc[:])/rho[iz]*100e0, dicx[iz]*1e3, alkx[iz]*1e3  \
                , co3x[iz]*1e3-co3sat, np.sum(rcc[iz,:]),-np.log10(prox[iz]) 
        file_tmp.close()
        file_tmp = open(workdir+'omx-'+ '{:03}'.format(itrec)+'.txt','w') 
        for iz in range(nz):
            print>>file_tmp,z[iz],age[iz],omx[iz]*mom/rho[iz]*100e0
        file_tmp.close()
        file_tmp = open(workdir+'o2x-'+ '{:03}'.format(itrec)+'.txt','w') 
        for iz in range(nz):
            print>>file_tmp,z[iz],age[iz],o2x[iz]*1e3, oxco2[iz], anco2[iz]
        file_tmp.close()
        file_tmp = open(workdir+'ccx_sp-'+ '{:03}'.format(itrec)+'.txt','w') 
        for iz in range(nz):
            # print>>file_tmp,z[iz],age[iz],ccx[iz,0:nspcc]*mcc/rho[iz]*100e0
            print>>file_tmp,z[iz],age[iz], str(ccx[iz,0:nspcc]*mcc[0:nspcc]/rho[iz]*100.)[1:-1]
        file_tmp.close()
        file_tmp = open(workdir+'sig-'+ '{:03}'.format(itrec)+'.txt','w') 
        for iz in range(nz):
            print>>file_tmp,z[iz],age[iz],d13c_blk[iz],d18o_blk[iz]
        file_tmp.close()
        file_tmp = open(workdir+'bur-'+ '{:03}'.format(itrec)+'.txt','w') 
        for iz in range(nz):
            print>>file_tmp,z[iz],age[iz],w[iz],up[iz],dwn[iz],cnr[iz],adf[iz]
        file_tmp.close()
    return 
    
def timestep(time,time_spn,time_trs,time_aft
    ,nt_spn,nt_trs,nt_aft
    ,flg_restart,dt
    ):
    if time <= time_spn:
        # dt=time_spn/np.float64(nt_spn)
        if not flg_restart:
            if time+dt>time_spn:
                dt=time_trs/np.float64(nt_trs)
            else:
                dt=time_spn/np.float64(nt_spn)
    elif time>time_spn and time<=time_spn+time_trs:
        if not flg_restart: dt=time_trs/np.float64(nt_trs)
    elif time>time_spn+time_trs:
        if not flg_restart:dt=time_trs/np.float64(nt_aft) #! not too large time step
    return dt

def signal_flx(time,time_spn,time_trs,time_aft
    ,d13c_ocni,d18o_ocni,d13c_ocnf,d18o_ocnf,ccflxi,ccflx,track2,size,biotest
    ,d13c_sp,d18o_sp,flxfini,flxfinf,nspcc,it
    ):
    if time <= time_spn:
        d13c_ocn = d13c_ocni
        d18o_ocn = d18o_ocni
        ccflx[:] = 0.
        ccflx[0] = ccflxi
        if track2:
            cntsp=0
            d18o_sp[cntsp]=d18o_ocn
            d13c_sp[cntsp]=d13c_ocn
            ccflx[:]=0.
            ccflx[cntsp]=ccflxi
        if size:
            ccflx[:]=0.
            flxfin = flxfini
            ccflx[0]=(1.-flxfin)*ccflxi
            ccflx[4]=flxfin*ccflxi
    elif time>time_spn and time<=time_spn+time_trs:
        d13c_ocn = d13c_ocni + (time-time_spn)*(d13c_ocnf-d13c_ocni)/time_trs
        d18o_ocn = d18o_ocni + (time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs
        if time-time_spn<=time_trs/2.:
            d18o_ocn = d18o_ocni + (time-time_spn)*(d18o_ocnf-d18o_ocni)\
                       /time_trs*2.   
            flxfin = flxfini + (time-time_spn)*(flxfinf-flxfini)/time_trs*2.           
        else: 
            d18o_ocn = 2.*d18o_ocnf - d18o_ocni \
                       -(time-time_spn)*(d18o_ocnf-d18o_ocni)/time_trs*2.
            flxfin = 2.*flxfinf - flxfini \
                     - (time-time_spn)*(flxfinf-flxfini)/time_trs*2.
        if not biotest:
            if time-time_spn<=time_trs/10.:
                d13c_ocn = d13c_ocni \
                           + (time-time_spn)*(d13c_ocnf-d13c_ocni)\
                           /time_trs*10e0
            elif time-time_spn>time_trs/10e0 \
                 and time-time_spn<=time_trs/10e0*9e0:
                d13c_ocn = d13c_ocnf
            elif  time-time_spn>time_trs/10e0*9e0: 
                d13c_ocn = 10e0*d13c_ocnf - 9e0*d13c_ocni \
                           - (time-time_spn)*(d13c_ocnf-d13c_ocni)\
                           /time_trs*10e0
        if not(d13c_ocn>=d13c_ocnf and d13c_ocn<=d13c_ocni):
            #check if calculated d13c and d18o are within the assumed ranges  
            print 'error in d13c',d13c_ocn
            input()
        flxfrc = np.zeros(nspcc)
        flxfrc2 = np.zeros(nspcc)
        flxfrc[0] = abs(d13c_ocnf-d13c_ocn)\
                    /(abs(d13c_ocnf-d13c_ocn)+abs(d13c_ocni-d13c_ocn))
        flxfrc[1] = abs(d13c_ocni-d13c_ocn)\
                    /(abs(d13c_ocnf-d13c_ocn)+abs(d13c_ocni-d13c_ocn))
        flxfrc[2] = abs(d18o_ocnf-d18o_ocn)\
                    /(abs(d18o_ocnf-d18o_ocn)+abs(d18o_ocni-d18o_ocn))
        flxfrc[3] = abs(d18o_ocni-d18o_ocn)\
                    /(abs(d18o_ocnf-d18o_ocn)+abs(d18o_ocni-d18o_ocn))
        while True:
            # case flxfrc2[0]=0
            flxfrc2[:]=0e0
            flxfrc2[1] = flxfrc[0]
            flxfrc2[2] = flxfrc[2]
            flxfrc2[3] = flxfrc[1]-flxfrc2[2]
            if all(flxfrc2>=0e0):break 
            flxfrc2[:]=0e0
            # case flxfrc2(2)=0
            flxfrc2[0] = flxfrc[0]
            flxfrc2[2] = flxfrc[2]-flxfrc2[0]
            flxfrc2[3] = flxfrc[3]
            if all(flxfrc2>=0e0):break 
            flxfrc2[:]=0e0
            # case flxfrc2(3)=0
            flxfrc2[0] = flxfrc[2]
            flxfrc2[1] = flxfrc[0]-flxfrc2[0]
            flxfrc2[3] = flxfrc[1]
            if all(flxfrc2>=0e0):break 
            flxfrc2[:]=0e0
            # case flxfrc2(4)=0
            flxfrc2[1] = flxfrc[3]
            flxfrc2[0] = flxfrc[0]-flxfrc2[1]
            flxfrc2[2] = flxfrc[1]
            if all(flxfrc2>=0e0):break 
            print 'error' #! should not come here 
            input()
        if track2:
            if time-time_spn<=time_trs/10e0:
                if it%int(5000/(nspcc-2))==0:    
                    cntsp=cntsp+1 # new species assinged 
                    d18o_sp[cntsp]=d18o_ocn
                    d13c_sp[cntsp]=d13c_ocn
                    ccflx[:] = 0e0
                    ccflx[cntsp] = ccflxi
            elif time-time_spn>time_trs/10e0 \
                 and time-time_spn<=time_trs/10e0*9e0:
                if it%int(5000/(nspcc-2))==0:   
                    cntsp=cntsp+1
                    d18o_sp[cntsp]=d18o_ocn
                    d13c_sp[cntsp]=d13c_ocn
                    ccflx[:] = 0e0
                    ccflx[cntsp] = ccflxi
            elif  time-time_spn>time_trs/10e0*9e0: 
                if it%int(5000/(nspcc-2))==0:   
                    cntsp=cntsp+1
                    d18o_sp[cntsp]=d18o_ocn
                    d13c_sp[cntsp]=d13c_ocn
                    ccflx[:] = 0e0
                    ccflx[cntsp] = ccflxi
        for isp in range(nspcc):
            ccflx[isp]=flxfrc2[isp]*ccflxi
        if size:
            for isp in range(nspcc):
                if isp<=3:ccflx[isp]=flxfrc2[isp]*ccflxi*(1.-flxfin)
                else:ccflx[isp]=flxfrc2[isp-4]*ccflxi*flxfin
    elif time>time_spn+time_trs:
        d13c_ocn = d13c_ocni # now again initial values 
        d18o_ocn = d18o_ocni
        ccflx[:] = 0e0
        ccflx[0] = ccflxi
        if track2:
            if cntsp+1!=nspcc-1: # checking used caco3 species number
                # is enough and necessary 
                print 'fatal error in counting',cntsp
                input()
            d18o_sp[cntsp+1]=d18o_ocn
            d13c_sp[cntsp+1]=d13c_ocn
            ccflx[:] = 0e0
            ccflx[cntsp+1] = ccflxi
        if size:
            ccflx[:]=0.
            flxfin = flxfini
            ccflx[0] = (1e0-flxfin)*ccflxi  # fine species 
            ccflx[4] = flxfin*ccflxi # coarse species
        if biotest: 
            d13c_ocn = d13c_ocnf   # finish with final value 
            d18o_ocn = d18o_ocni
            ccflx[:] = 0e0
            ccflx[2] = ccflxi
    return d13c_ocn,d18o_ocn,d13c_sp,d18o_sp,ccflx
    
def bdcnd(time,time_spn,time_trs,time_aft,depi,biotest,depf):
    if time <= time_spn:
        dep = depi
    elif time>time_spn and time<=time_spn+time_trs:
        if not biotest:
            if time-time_spn<=time_trs/10.:
                dep = depi + (depf-depi)*(time-time_spn)/time_trs*10e0
            elif time-time_spn>time_trs/10e0 \
                 and time-time_spn<=time_trs/10e0*9e0:
                dep = depf
            elif  time-time_spn>time_trs/10e0*9e0: 
                dep = 10e0*depf-9e0*depi \
                      - (depf-depi)*(time-time_spn)/time_trs*10e0
        else:
            dep = depi
    elif time>time_spn+time_trs:
        dep = depi
    return dep 

def caco3_main(ccflxi,om2cc,dep,dt,fl,biot,oxonly,runmode,co2chem,sparse,showiter):    
    # ------------------- parameter inputs
    nz = np.int32(100)              # total grid number 
    nspcc = np.int32(4)             # caco3 species 
    ccflxi = np.float64(ccflxi)     # total caco3 rain flux [umol cm2 yr-1]
    om2cc = np.float64(om2cc)       # om/caco3 rain ratio 
    rhocc = np.float64(2.71)        # density of caco3 [g cm-3]
    rhosed = np.float64(2.6)        # density of clay [g cm-3]
    rhoom = np.float64(1.2)         # density of om [g cm-3]
    mcc = np.float64(100.)          # molar mass of caco3 [g mol-1]
    msed = np.float64(258.16)       # molar mass of clay [g mol-1]
    mom = np.float64(30.)           # molar mass of om [g mol-1]
    ox2om = np.float64(1.3)         # o2/om mole ratio consumed upon oxic degradation of om 
    d2yr = np.float64(365.25)       # days/year 
    kcci = np.float64(1.*d2yr)      # ref. caco3 dissoltion rate const. [yr-1]
    komi = np.float64(0.06)         # ref. om degradation rate const. [yr-1]
    ncc = np.float64(4.5)           # reaction order wrt caco3 dissolution 
    fact = np.float64(1e-3)         # arbitrary factor to facilitate calculation 
    temp = np.float64(2.)           # temperature [C]
    poroi = np.float64(0.8)         # ref. porosity
    sal = np.float64(35.)           # salinity [o/oo]
    dep = np.float64(dep)           # water depth [km]
    cai = np.float64(10.3e-3)       # seawater ca conc. [mol kg-1]
    o2i = np.float64(165.)          # seawater o2 conc. [umol kg-1]
    alki = np.float64(2285.)        # seawater alk conc. [umol kg-1]
    dici = np.float64(2211.)        # seaawter dic conc. [umol kg-1]
    o2th = np.float64(0.)           # threshold o2 conc. [mol cm-3]
    zmlref = np.float64(12.)        # ref. mixled layer depth [cm]
    ccx_th = np.float64(1e-300)     # threshold caco3 conc. [mol cm-3]
    omx_th = np.float64(1e-300)     # threshold om conc. [mol cm-3]
    ztot = np.float64(500.)         # total depth of sediment [cm]
    nrec = np.int32(15)             # total number of recording 
    tol = np.float64(1e-6)          # tolerance value of reltaive error
    dt = np.float64(dt)             # time step [yr]
    itr_w_max = 20                  # maximum iteration number for w 
    # swithces 
    oxic = True                     # allowing oxic degradation of om 
    anoxic = True                   # anoxic degradation of om 
    sense = False                   # not tracking signals
    biotest = False                 # 5kyr signal tracking
    size = False                    # two sizes of caco3 species
    track2 = False                  # signal tracking by method 2
    isotrack = False                # direct isotopologue tracking
    showiter = showiter             # show each iteration 
    sparse = sparse                 # use sparse matrix solver for co2 system
    # switches for mixing 
    allturbo2 = False               # all turbo2-like mixing
    alllabs = False                 # all labs mixing
    allnobio = False                # no bioturbation for all solid species 
    # switches for co2 chemistry 
    co2chem = co2chem               # dic = co2+hco3+co3; alk = hco3+2*co3
    # working directory 
    workdir = 'C:/Users/YK/Desktop/Sed_res'
    # file nane 
    if len(fl)==0:filename = '-no_name_specified' 
    else: filename = '-'+fl
    # choose bioturbation style (Fickian mixing as default)
    if biot == 'nobio':allnobio = True
    elif biot == 'turbo2':allturbo2 = True
    elif biot == 'labs': alllabs = True
    # consider only oxic degrdation of OM if selected
    if oxonly: anoxic = False
    # option of model run 
    if runmode == 'sense':sense = True
    elif runmode== 'biotest': biotest = True
    elif runmode == 'size': size = True
    elif runmode == 'track2': track2 = True       
    elif runmode == 'isotrack': isotrack = True
    # ===========  checking something 
    # chk_caco3_therm()
    # chk_caco3_therm_sbrtns()
    # ------------------- dependent parameter calculations
    nspcc = np.int32(4) # default species number for tracking 2 isotope signals 
    if sense:nspcc =np.int32(1)
    if size:nspcc = np.int32(8)
    if track2:nspcc = np.int32(42)
    if isotrack:nspcc = np.int32(5)
    # use a shallow sediment for sensitivity analysis
    if sense:ztot = 50.
    # assign IDs to isotopologues; reference isotope ratios  
    if isotrack: 
        i12c16o=1;i12c18o=2;i13c16o=3;i13c18o=4;i14c=5
        r18o_pdb = 0.0020672 ; r17o_pdb = 0.0003859 ; r13c_pdb = 0.011180
    # mixing properties 
    nobio = np.zeros(nspcc+2,dtype=bool)
    turbo2 = np.zeros(nspcc+2,dtype=bool)
    labs = np.zeros(nspcc+2,dtype=bool)
    if allturbo2: turbo2[:] = True
    if alllabs: labs[:] = True
    if allnobio: nobio[:] = True
    #  preparing directory to store results 
    workdir += '/test-translabs/profiles/python/multi/'
    if not anoxic: workdir += 'ox'
    else: workdir += 'oxanox'
    if any(labs): workdir += '_labs'
    if any(turbo2): workdir += '_turbo2'
    if any(nobio): workdir += '_nobio'
    workdir += '/'
    workdir += 'cc-'+'{:.1e}'.format(ccflxi)+'_rr-'+'{:.1e}'.format(om2cc)
    if sense: workdir += '_dep-'+'{:.1e}'.format(dep)
    else: workdir += filename
    if not os.path.exists(workdir):os.makedirs(workdir)
    workdir += '/'
    file_ptflx=open(workdir+'ptflx.txt','w')
    file_ccflx=open(workdir+'ccflx.txt','w')
    file_omflx=open(workdir+'omflx.txt','w')
    file_o2flx=open(workdir+'o2flx.txt','w')
    file_dicflx=open(workdir+'dicflx.txt','w')
    file_alkflx=open(workdir+'alkflx.txt','w')
    file_err=open(workdir+'errlog.txt','w')
    file_bound=open(workdir+'bound.txt','w')
    file_totfrac=open(workdir+'frac.txt','w')
    file_sigmly=open(workdir+'sigmly.txt','w')
    file_sigmlyd=open(workdir+'sigmlyd.txt','w')
    file_sigbtm=open(workdir+'sigbtm.txt','w')
    for isp in range(nspcc):
        f=open(workdir+'ccflx-sp_'+'{:03}'.format(isp+1)+'.txt','w')
        f.close()
    # making grid
    beta = 1. + 5e-11
    dz,z = makegrid(beta,ztot,nz)
    ### FUNDAMENTAL PARAMETERS 
    omflx = np.float64(om2cc*ccflxi)
    detflx = np.float64(1./9.*ccflxi*mcc)
    ccflx = np.zeros((nspcc),dtype=np.float64)
    ccflx[:] = ccflxi/float(nspcc)
    # ccflx[0]=ccflxi
    mcc = np.ones(nspcc,dtype=np.float64)*100.
    if isotrack: 
        mcc[i12c16o] = 40.078+12.+16.*3.
        mcc[i12c18o] = 40.078+12.+16.*2.+18.
        mcc[i13c16o] = 40.078+13.+16.*3.
        mcc[i13c18o] = 40.078+13.+16.*2.+18.
        mcc[i14c] = 40.078+14.+16.*3.
    mvom = np.float64(mom/rhoom)
    mvsed = np.float64(msed/rhosed)
    mvcc = np.ones(nspcc,dtype=np.float64)
    mvcc[:] = mcc[:]/rhocc
    sporo = np.zeros(nz,dtype=np.float64)
    poro = getporosity(z,nz)
    porof = poro[nz-1]  # this assumes zero-porosity gradient at the depth; these choices do not affect the calculation 
    sporof = 1.-porof  #  volume fraction of solids at bottom depth
    sporoi = 1.-poroi # volume fraction of solids at the seawater-sediment interface (SWI)
    sporo = 1. - poro  #  volume fraction of solids
    wi = (detflx/msed*mvsed + np.sum(ccflx[:]*mvcc[:]))/(1.-poroi)
    # here reference point of iteration
    w = np.zeros(nz,dtype=np.float64)
    wx = np.zeros(nz,dtype=np.float64)
    wxx = np.zeros(nz,dtype=np.float64)
    w_pre = np.zeros(nz,dtype=np.float64)
    w[:] = wi
    age = dep2age(z,dz,nz,w)
    # upwind(w)
    up,dwn,cnr,adf = calcupwindscheme(w,nz)
    #---------------
    zox = np.float64(10.)
    #------ recording time 
    # ++++ tracking experiment
    time_spn = np.float64(ztot/wi*50.)  # spinup
    time_trs = np.float64(50e3)   # transient event 
    time_aft = np.float64(time_trs *3.)
    if sense:
        time_trs = np.float64(0.)
        time_aft = np.float64(0.)
    if biotest:
        time_trs = np.float64(5e3)
        time_aft = np.float64(time_trs*10.)
    rectime = np.zeros(nrec)
    for itrec in range(0,nrec/3):
        rectime[itrec]=(itrec+1)*time_spn/float(nrec/3.)
    for itrec in range(nrec/3,nrec*2/3):
        rectime[itrec]=rectime[nrec/3-1]+(itrec-nrec/3+1)*time_trs/float(nrec/3.)
    for itrec in range(nrec*2/3,nrec):
        rectime[itrec]=rectime[nrec/3*2-1]+(itrec-nrec/3*2+1)*time_aft/float(nrec/3.)
    cntrec=0
    np.savetxt(workdir+'rectime.txt',rectime)
    depi = 3.5
    depf = dep
    flxfini = 0.5
    flxfinf = 0.9
    #//////// isotopes //////////
    d13c_blk = np.zeros(nz,dtype=np.float64)
    d13c_blkc = np.zeros(nz,dtype=np.float64)
    d13c_blkf = np.zeros(nz,dtype=np.float64)
    d18o_blk = np.zeros(nz,dtype=np.float64)
    d18o_blkc = np.zeros(nz,dtype=np.float64)
    d18o_blkf = np.zeros(nz,dtype=np.float64)
    d13c_ocni = 2.
    d13c_ocnf = -1.
    d18o_ocni = 1.
    d18o_ocnf = -1.
    d13c_sp = np.zeros(nspcc)
    d18o_sp = np.zeros(nspcc)
    if not sense:
        d13c_sp[0]=d13c_ocni
        d18o_sp[0]=d18o_ocni
        d13c_sp[1]=d13c_ocni
        d18o_sp[1]=d18o_ocnf
        d13c_sp[2]=d13c_ocnf
        d18o_sp[2]=d18o_ocni
        d13c_sp[3]=d13c_ocnf
        d18o_sp[3]=d18o_ocnf
    else:
        d18o_sp[:]=0.
        d13c_sp[:]=0.
    if size:
        d13c_sp[0+4]=d13c_ocni
        d18o_sp[0+4]=d18o_ocni
        d13c_sp[1+4]=d13c_ocni
        d18o_sp[1+4]=d18o_ocnf
        d13c_sp[2+4]=d13c_ocnf
        d18o_sp[2+4]=d18o_ocni
        d13c_sp[3+4]=d13c_ocnf
        d18o_sp[3+4]=d18o_ocnf
    #!!!!!!!!! TRANSITION MATRIX !!!!!!!!!!!!!!!!!!!
    trans,izrec,izrec2,izml,nonlocal = make_transmx(  
         labs,nspcc,turbo2,nobio,dz,sporo,nz,z,zmlref,size  # input
        )
    # ~~~~~~~~~~~~~~ diffusion & reaction~~~~~~~~~~~~~~
    dif_dic,dif_alk,dif_o2,kom,kcc,co3sat = coefs(cai,temp,nz,nspcc,poro,komi,kcci,size,sal,dep)
    # ~~~~~~~~~~~~~~~~ initial conditions 
    cc = np.zeros((nz,nspcc),dtype=np.float64)
    dic = np.zeros(nz,dtype=np.float64)
    alk=np.zeros(nz,dtype=np.float64)
    co2=np.zeros(nz,dtype=np.float64)
    hco3=np.zeros(nz,dtype=np.float64)
    co3=np.zeros(nz,dtype=np.float64)
    pro=np.zeros(nz,dtype=np.float64)
    info = 0
    cc[:,:] = 1e-8
    dic[:]=dici*1e-6/1e3
    alk[:]=alki*1e-6/1e3
    if co2chem=='co2': co2,hco3,co3,pro,dco3_ddic,dco3_dalk,info = calcspecies(temp,sal,dep,dic,alk,nz)
    elif co2chem=='co2h2o':co2,hco3,co3,pro,dco3_ddic,dco3_dalk,info = calcco2h2o(temp,sal,dep,dic,alk,nz)
    elif co2chem=='mocsy':
        co2,hco3,co3,pro,ohmega,dohmega_ddic,dohmega_dalk = co2sys_mocsy(nz,alk*1e6,dic*1e6,temp,dep*1e3,sal)
        co2 = co2/1e6
        hco3 = hco3/1e6
        co3 = co3/1e6
    pt=np.zeros(nz,dtype=np.float64)
    om=np.zeros(nz,dtype=np.float64)
    o2=np.zeros(nz,dtype=np.float64)
    pt[:] = 1e-8
    om[:] = 1e-8
    o2[:] = o2i*1e-6/1e3
    ccx=np.zeros((nz,nspcc),dtype=np.float64)
    dicx=np.zeros(nz,dtype=np.float64)
    alkx=np.zeros(nz,dtype=np.float64)
    co2x=np.zeros(nz,dtype=np.float64)
    hco3x=np.zeros(nz,dtype=np.float64)
    co3x=np.zeros(nz,dtype=np.float64)
    prox=np.zeros(nz,dtype=np.float64)
    ptx=np.zeros(nz,dtype=np.float64)
    omx=np.zeros(nz,dtype=np.float64)
    o2x=np.zeros(nz,dtype=np.float64)
    ccx[:,:]=cc[:,:].copy()
    dicx[:]=dic[:].copy()
    alkx[:]=alk[:].copy()
    hco3x[:]=hco3[:].copy()
    co3x[:]=co3[:].copy()
    co3i = co3[0]
    if co2chem=='mocsy':        
        cai = (0.02128e0/40.078e0) * sal/1.80655e0
        co3sat = co3i*1e3/ohmega[0]
    ptx[:]=pt[:].copy()
    omx[:]=om[:].copy()
    o2x[:]=o2[:].copy()
    frt=np.zeros(nz,dtype=np.float64)
    rho=np.zeros(nz,dtype=np.float64)
    rcc=np.zeros((nz,nspcc),dtype=np.float64)
    omega=np.zeros((nz,nspcc),dtype=np.float64)
    oxco2 = np.zeros(nz,dtype=np.float64)
    anco2 = np.zeros(nz,dtype=np.float64)
    for isp in range(nspcc):
        rcc[:,isp]=kcc[:,isp]*ccx[:,isp]*abs(1.-co3[:]*1e3/co3sat)**ncc\
                    *((1.-co3[:]*1e3/co3sat)>0.).astype(float)
    # recording
    recordprofile(0 
        ,nz,z,age,pt,ptx,msed,w,wi,rho,frt 
        ,cc,mcc,dic,alk,co3,co3sat,rcc,pro,dicx,alkx,co3x,prox,nspcc
        ,o2x,oxco2,anco2
        ,om,omx,mom
        ,ccx
        ,d13c_ocni,d18o_ocni,d13c_blk,d18o_blk
        ,up,dwn,cnr,adf
        ,workdir
        )
    #START OF TIME INTEGRAL 
    time=0.
    it=1
    flg_restart = False
    nt_spn=400;nt_trs=5000;nt_aft=1000
    dt_save = dt
    while True: # time loop 
        if not sense:
            dt = timestep(time,time_spn,time_trs,time_aft
                ,nt_spn,nt_trs,nt_aft
                ,flg_restart,dt
                ) 
            d13c_ocn,d18o_ocn,d13c_sp,d18o_sp,ccflx = signal_flx(time,time_spn,time_trs,time_aft
                ,d13c_ocni,d18o_ocni,d13c_ocnf,d18o_ocnf,ccflxi,ccflx,track2,size,biotest
                ,d13c_sp,d18o_sp,flxfini,flxfinf,nspcc,it
                )
            dep = bdcnd(time,time_spn,time_trs,time_aft,depi,biotest,depf)
        else:
            dt = dt_save
            d13c_ocn = 0. 
            d18o_ocn = 0.
        d13c_flx = np.sum(d13c_sp[:]*ccflx[:])/ccflxi
        d18o_flx = np.sum(d18o_sp[:]*ccflx[:])/ccflxi
        if track2:
            if abs(d13c_flx - d13c_ocn)>tol or abs(d18o_flx - d18o_ocn)>tol:
                print 'error in assignment of proxy'
                print >> file_err, 'error in assignment of proxy'\
                      ,d18o_ocn,d13c_ocn,d18o_flx,d13c_flx
                input()
        ## === temperature & pressure and associated boundary changes ====
        # if temperature is changed during signal change event this affect
        # diffusion coeff etc. 
        dif_dic,dif_alk,dif_o2,kom,kcc,co3sat = coefs(cai,temp,nz,nspcc,poro,komi,kcci,size,sal,dep)
        
        if it==1:
            print >> file_bound, '#time  d13c_ocn  d18o_ocn, fluxes of cc:'\
                  ,np.linspace(1,nspcc,nspcc),'temp  dep  sal  dici  alki  o2i'
        if not size:
            print>> file_bound, time, d13c_ocn, d18o_ocn, str(ccflx)[1:-1]\
                    ,temp, dep, sal,dici,alki, o2i
        else:
            print>> file_bound, time, d13c_ocn, d18o_ocn\
                    , np.sum(ccflx[0:4]),np.sum(ccflx[4:8]),ccflx[:]\
                    ,temp, dep, sal,dici,alki, o2i
        itr_w = 0
        err_w_min = 1e4
#        itr_f = 0
#        err_f_min = 1e4
#        dfrt_df = 0.
#        d2frt_df2 = 0.
        err_f = 0.
#        err_fx = 0.
        w_pre[:]=w[:].copy()
        #  here should be a reference point for iteration of w
        while True:  # burial loop 
            print 'it :',it,dt
            dw = np.zeros(nz,dtype=np.float64)
##            oxco2[:]=0.
##            anco2[:]=0.
            itr =0
            error = 1e4
            minerr = 1e4
            izox = nz
            # om & o2 iteration wrt zox
#            while error>tol:  # om & o2 loop 
            while True:  # om & o2 loop 
                izox,kom,zox,kom_ox,kom_an = calc_zox( 
                    oxic,anoxic,nz,o2x,o2th,komi,ztot,z,o2i,dz  # input
                    )
                omx = omcalc( 
                    kom,omx     # in&output
                    ,oxic,anoxic,o2x,om,nz,sporo,sporoi,sporof,o2th,komi  # input 
                    ,w,wi,dt,up,dwn,cnr,adf,trans,nspcc,labs,turbo2,nonlocal,omflx,poro,dz  # input 
                    ) 
                omadv,omdec,omdif,omrain,omres,omtflx,flg_restart = calcflxom(  
                    sporo,om,omx,dt,w,dz,z,nz,turbo2,labs,nonlocal,poro,up,dwn,cnr
                    ,adf,rho,mom,trans,kom,sporof,sporoi,wi,nspcc,omflx,workdir  # input 
                    )
                if flg_restart:
                    dt = dt/10e0
                    w[:]=w_pre[:]
                    up,dwn,cnr,adf = calcupwindscheme(w,nz)
                    print 'must restart after om calc'
                    break
                    print 'stop'
                #~~~~~~~~~~~~~~~~~ O2 calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                # reset matrices 
                # if izox == nz -1 : # fully oxic; lower boundary condition ---> no diffusive out flow 
                o2x = o2calc_ox(  
                    izox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
                    )
                if showiter: print  'o2 :',itr, o2x[0:nz:nz/10]
                #  fluxes relevant to o2 (at the same time checking the satisfaction of difference equations) 
                o2dec,o2dif,o2tflx,o2res = calcflxo2_ox( 
                    nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt,o2,o2x,ox2om,o2i  # input
                    )
                if (o2x[:]>=0e0).all() and izox==nz-1: 
                    iizox_errmin = nz-1
                    all_oxic = True
                elif any(o2x[:]<0e0): 
                    all_oxic = False
                    error_o2min = 1e4
                    iizox_errmin = izox
                    for iizox in range(nz):   ## if oxygen is depleted within calculation domain, lower boundary changes to zero concs.
#                        print iizox                        
                        if iizox<nz-1: 
                            o2x = o2calc_sbox(  
                                iizox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
                                )
                        elif iizox==nz-1: 
                            o2x = o2calc_ox(  
                                izox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
                                )
                        if showiter: print  'o2 :',iizox, o2x[0:nz:nz/10]
                        if (o2x[:]>=0e0).all():
                            if np.abs(o2x[max(iizox-1,0)])<error_o2min: 
                                error_o2min = np.abs(o2x[max(iizox-1,0)])
                                iizox_errmin = iizox
                    if iizox_errmin<nz-1:
                        o2x = o2calc_sbox(  
                            iizox_errmin,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
                            )
                        if showiter: print  'o2 :',iizox_errmin, o2x[0:nz:nz/10]
                        o2dec,o2dif,o2tflx,o2res = calcflxo2_sbox( 
                            nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt,o2,o2x,iizox_errmin,ox2om,o2i  # input
                            )
                    elif iizox_errmin==nz-1:
                        o2x = o2calc_ox(  
                            izox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
                            )
                        if showiter: print  'o2 :',iizox_errmin, o2x[0:nz:nz/10]
                        o2dec,o2dif,o2tflx,o2res = calcflxo2_ox( 
                            nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt,o2,o2x,ox2om,o2i  # input
                            )
                    iizox_errmin,kom_dum,zox,kom_ox_dum,kom_an_dum = calc_zox( 
                        oxic,anoxic,nz,o2x,o2th,komi,ztot,z,o2i,dz  # input
                        )
                    if showiter: print iizox_errmin
#                    if iizox_errmin2!=iizox_errmin: 
#                        print 'iizox_errmin2!=iizox_errmin',iizox_errmin2,iizox_errmin
#                        input()
                error = abs(izox-iizox_errmin)  #  difference 
                if showiter: print  'zox',itr,izox, iizox_errmin
                if izox==iizox_errmin: 
                    if not all_oxic:
                        if izox<nz-1:
                            o2x = o2calc_sbox(  
                                izox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
                                )
                            if showiter: print  'o2 :',izox, o2x[0:nz:nz/10]
                            o2dec,o2dif,o2tflx,o2res = calcflxo2_sbox( 
                                nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt,o2,o2x,izox,ox2om,o2i  # input
                                )
                        elif izox == nz-1: 
                            o2x = o2calc_ox(  
                                izox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
                                )
                            if showiter: print  'o2 :',izox, o2x[0:nz:nz/10]
                            o2dec,o2dif,o2tflx,o2res = calcflxo2_ox( 
                                nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt,o2,o2x,ox2om,o2i  # input
                                )
                    break 
                if error < minerr :  
                    minerr = error 
                else: 
                    if izox < nz-1 and iizox_errmin == nz-1:  
                        o2x = o2calc_sbox(  
                            izox,nz,poro,o2,kom_ox,omx,sporo,dif_o2,dz,dt,o2i,ox2om,o2x # input
                            )
                        o2dec,o2dif,o2tflx,o2res = calcflxo2_sbox( 
                            nz,sporo,kom_ox,omx,dz,poro,dif_o2,dt,o2,o2x,izox,ox2om,o2i  # input
                            )
                        break
                itr = itr + 1
            #~~  O2 calculation END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if flg_restart:
                flg_timereset = True
                continue 
#            for iz in range(nz):
#                if o2x[iz] > o2th:
#                    oxco2[iz] = (1e0-poro[iz])*kom[iz]*omx[iz]  # aerobic respiration 
#                else :
#                    # o2x[iz] = o2th
#                    if anoxic: 
#                        anco2[iz] = (1e0-poro[iz])*kom[iz]*omx[iz]  # anaerobic respiration
            oxco2[:] = (1e0-poro[:])*kom_ox[:]*omx[:]
            anco2[:] = (1e0-poro[:])*kom_an[:]*omx[:]
            for iz in range(nz):
                dw[iz] = dw[iz] -(1e0-poro[iz])*mvom*kom[iz]*omx[iz]  ## burial rate change need reflect volume change caused by chemical reactions 
                # as well as non-local mixing 
                if turbo2[0] or labs[0]: 
                    for iiz in range( nz):
                        if trans[iiz,iz,0]==0e0: continue
                        dw[iz] = dw[iz] - mvom*(-trans[iiz,iz,0]/dz[iz]*dz[iiz]*(1e0-poro[iiz])*omx[iiz])
                else :
                    if nonlocal[0]: 
                        for iiz in range( nz):
                            if trans[iiz,iz,0]==0e0: continue
                            dw[iz] = dw[iz] - mvom*(-trans[iiz,iz,0]/dz[iz]*omx[iiz])
            for iz in range(nz):
                if omx[iz]<omx_th: omx[iz]=omx_th  ## truncated at minimum value 
            ##  ~~~~~~~~~~~~~~~~~~~~~~ CaCO3 solid, ALK and DIC  calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ccx,dicx,alkx,rcc,dt,flg_restart,w = calccaco3sys(  #
                ccx,dicx,alkx,rcc,dt  # in&output
                ,nspcc,dic,alk,dep,sal,temp,labs,turbo2,nonlocal,sporo,sporoi,sporof,poro,dif_alk,dif_dic # input
                ,w,up,dwn,cnr,adf,dz,trans,cc,oxco2,anco2,co3sat,kcc,ccflx,ncc,omega,nz,tol,sparse,fact 
                ,dici,alki,ccx_th,showiter,w_pre,co2chem,mcc,rho,workdir  # input
                )
            # ~~~~  End of calculation iteration for CO2 species ~~~~~~~~~~~~~~~~~~~~
            # update aqueous co2 species 
            # calcspecies(temp,sal,dep,dicx,alkx,co2x,hco3x,co3x,prox,info)
            # if flg_restart: break # exit loop for w iteration  
            if flg_restart: 
                flg_timereset = True
                continue # exit loop for w iteration  
            # calcspecies(temp,sal,dep,dicx,alkx)
            if co2chem=='co2':co2x,hco3x,co3x,prox,dco3_ddic,dco3_dalk,info = calcspecies(temp,sal,dep,dicx,alkx,nz)
            elif co2chem=='co2h2o':co2x,hco3x,co3x,prox,dco3_ddic,dco3_dalk,info = calcco2h2o(temp,sal,dep,dicx,alkx,nz)
            elif co2chem=='mocsy':
                co2x,hco3x,co3x,prox,ohmega,dohmega_ddic,dohmega_dalk = co2sys_mocsy(nz,alkx*1e6,dicx*1e6,temp,dep*1e3,sal)
                co2x = co2x/1e6
                hco3x = hco3x/1e6
                co3x = co3x/1e6
            if info==1:  
                dt=dt/10e0
                print 'stop'
                input()
            # calculation of fluxes relevant to caco3 and co2 system
            cctflx,ccflx,ccdis,ccdif,ccadv,ccrain,ccres,alktflx,alkdis,alkdif,alkdec,alkres \
                ,dictflx,dicdis,dicdif,dicres,dicdec  \
                ,dw = \
                calcflxcaco3sys(  
                    dw # inoutput
                     ,nspcc,ccx,cc,dt,dz,rcc,adf,up,dwn,cnr,w,dif_alk,dif_dic,dic,dicx,alk,alkx,oxco2,anco2,trans    # input
                     ,turbo2,labs,nonlocal,sporof,it,nz,poro,sporo,ccflx,dici,alki,mvcc,tol,workdir        # input
                     )
            # ~~~~ calculation clay  ~~~~~~~~~~~~~~~~~~
            ptx = claycalc( 
                ptx
                ,nz,sporo,pt,dt,w,dz,detflx,adf,up,dwn,cnr,trans  # input
                ,nspcc,labs,turbo2,nonlocal,poro,sporof,msed,workdir     # intput
                )
            pttflx,ptdif,ptadv,ptres,ptrain = calcflxclay( 
                dw         # in&output
                ,nz,sporo,ptx,pt,dt,dz,detflx,w,adf,up,dwn,cnr,sporof,trans,nspcc,turbo2,labs,nonlocal,poro           #  input
                ,msed,mvsed
                )
            if it==0: 
                print >>file_o2flx,'time, o2dec, o2dif, o2tflx, o2res'
                print >>file_omflx, 'time, omtflx, omadv, omdec, omdif, omrain, omres'
                print >>file_ccflx, 'time, cctflx, ccdis, ccdif, ccadv, ccrain, ccres' 
                print >>file_dicflx, 'time, dictflx, dicdis, dicdif, dicdec,  dicres' 
                print >>file_alkflx,  'time, alktflx, alkdis, alkdif, alkdec, alkres' 
                print>>file_ptflx, 'time, pttflx, ptdif, ptadv, ptrain, ptres'
                for isp in range(nspcc):
                    f=open(workdir+'ccflx-sp_'+'{:03}'.format(isp+1)+'.txt','a')
                    print>>f,'time, cctflx, ccdis, ccdif, ccadv, ccrain, ccres' 
                    f.close()
            print>>file_o2flx, time,o2dec, o2dif,o2tflx,o2res
            print>>file_omflx, time, omtflx, omadv, omdec, omdif, omrain, omres
            print>>file_ccflx,time,np.sum(cctflx[:]), np.sum(ccdis[:]), np.sum(ccdif[:]), np.sum(ccadv[:]), np.sum(ccrain[:]), np.sum(ccres[:])
            print>>file_dicflx,time,dictflx, dicdis, dicdif, dicdec,  dicres 
            print>>file_alkflx, time,alktflx, alkdis, alkdif, alkdec, alkres 
            print>>file_ptflx, time, pttflx, ptdif, ptadv, ptrain, ptres
            for isp in range(nspcc):
                f=open(workdir+'ccflx-sp_'+'{:03}'.format(isp+1)+'.txt','a')
                print>>f,time,cctflx[isp], ccdis[isp], ccdif[isp], ccadv[isp], ccrain[isp], ccres[isp]
                f.close()
            ## ~~~~~~~~~End of clay calculation 
#            err_fx = np.max(np.abs(frt[:] - 1e0))  # recording previous error in total vol. fraction of solids 
            for iz in range(nz ):
                rho[iz] = omx[iz]*mom + ptx[iz]*msed +  np.sum(ccx[iz,:]*mcc[:])  # calculating bulk density 
                frt[iz] = omx[iz]*mvom + ptx[iz]*mvsed + np.sum(ccx[iz,:]*mvcc[:])  # calculation of total vol. fraction of solids 
            err_f = np.max(np.abs(frt[:] - 1e0))  # new error in total vol. fraction (must be 1 in theory) 
#            if err_f < err_fx: err_f_min = err_f  # recording minimum error 
            wx[:]=w[:]
            wi = (detflx/msed*mvsed + np.sum(ccflx[:]*mvcc[:]) +omflx*mvom)/(1e0-poroi)  # upper value; (1e0-poroi) is almost meaningless, see below 
            for iz in range(nz):
                if iz==0: 
                    w[iz]=((1e0-poroi)*wi + dw[iz]*dz[iz])/(1e0-poro[iz])
                else :
                    w[iz]=((1e0-poro[iz-1])*w[iz-1] + dw[iz]*dz[iz])/(1e0-poro[iz])
            error = 1e4
            # ------------ determine calculation scheme for advection 
            # upwind(w)
            up,dwn,cnr,adf = calcupwindscheme(w,nz)
            itr_w = itr_w + 1  # couting iteration for w 
            err_w = np.max(np.abs((w[:]-wx[:])/wx[:]))  # relative difference of w 
            if err_w<err_w_min: 
                err_w_min= err_w  # recording minimum relative difference of  w 
                wxx[:] = wx[:].copy()  # recording w which minimizes deviation of total sld fraction from 1 
            if itr_w>itr_w_max:   # if iteration gets too many, force to end with optimum w where error is minimum
                if itr_w==itr_w_max+1: 
                    w[:] = wxx[:].copy()   
                    continue
                elif itr_w==itr_w_max+2: 
                    w[:] = wxx[:].copy()
                    print>>file_err, 'not converging w',time, err_w, err_w_min
                    break
            if err_w > tol: continue
            else: break
        # if flg_restart:continue # going back
        if sense:
            if err_f < tol: break  # if total vol. fraction is near enough to 1, steady-state solution is obtained 
        ## depth -age conversion 
        age = dep2age(z,dz,nz,w)
        if any(rho<0e0):  # if negative density stop ....
            print 'negative density'
            file_tmp = open(workdir+'NEGATIVE_RHO.txt', 'w')
            for iz in range( nz):
                print>> file_tmp,z[iz],rho[iz],w[iz],up[iz],dwn[iz],cnr[iz],adf[iz]
            file_tmp.close()
            input()
        #/////// ISOTOPES /////
        for iz in range(nz ):
            d18o_blk[iz] = np.sum(d18o_sp[:]*ccx[iz,:])/np.sum(ccx[iz,:])
            d13c_blk[iz] = np.sum(d13c_sp[:]*ccx[iz,:])/np.sum(ccx[iz,:])
            if size:
                d18o_blkf[iz] = np.sum(d18o_sp[0:4]*ccx[iz,0:4])/np.sum(ccx[iz,0:4])
                d13c_blkf[iz] = np.sum(d13c_sp[0:4]*ccx[iz,0:4])/np.sum(ccx[iz,0:4])
                d18o_blkc[iz] = np.sum(d18o_sp[4:8]*ccx[iz,4:8])/np.sum(ccx[iz,4:8])
                d13c_blkc[iz] = np.sum(d13c_sp[4:8]*ccx[iz,4:8])/np.sum(ccx[iz,4:8])
        ##### PRINTING RESULTS ##################################
        if time>=rectime[cntrec]: 
            recordprofile(cntrec+1 
                ,nz,z,age,pt,ptx,msed,w,wi,rho,frt 
                ,cc,mcc,dic,alk,co3,co3sat,rcc,pro,dicx,alkx,co3x,prox,nspcc
                ,o2x,oxco2,anco2
                ,om,omx,mom
                ,ccx
                ,d13c_ocni,d18o_ocni,d13c_blk,d18o_blk
                ,up,dwn,cnr,adf
                ,workdir
                )
            cntrec = cntrec + 1
            if cntrec == nrec: break
        print  'time   :',time, np.max(np.abs(frt[:] - 1e0))
        print '~~~~ conc ~~~~'
        print  'z  :',z[0:nz:nz/5]
        print  'om :',omx[0:nz:nz/5]*mom/rho[0:nz:nz/5]*100e0
        print  'o2 :',o2x[0:nz:nz/5]#*1e3
        print  'cc :',np.sum(ccx[0:nz:nz/5,:]*mcc[:],axis=1)/rho[0:nz:nz/5]*100e0
        print  'dic:',dicx[0:nz:nz/5]#*1e3
        print  'alk:',alkx[0:nz:nz/5]#*1e3
        print  'ph :',-np.log10(prox[0:nz:nz/5])
        print  'sed:',ptx[0:nz:nz/5]*msed/rho[0:nz:nz/5]*100e0
        print  '   ..... multiple cc species ..... '
        for isp in range(nspcc ):
            print '{:03}'.format(isp+1),ccx[0:nz:nz/5,isp]*mcc[isp]/rho[0:nz:nz/5]*100e0
        print  '++++ flx ++++'
        print  'flx:',np.array(['tflx   ','adv    ','dif    ','omrxn  ','ccrxn  ','rain   ','res    '])
        print  'om :',np.array([ omtflx, omadv,  omdif, omdec,0e0,omrain, omres])
        print  'o2 :',np.array([o2tflx,0e0, o2dif,o2dec, 0e0,0e0,o2res])
        print  'cc :',np.array([np.sum(cctflx[:]),  np.sum(ccadv[:]), np.sum(ccdif[:]),0e0,np.sum(ccdis[:]), np.sum(ccrain[:]), np.sum(ccres[:]) ])
        print  'dic:',np.array([dictflx, 0e0,dicdif, dicdec,  dicdis, 0e0,dicres ])
        print  'alk:',np.array([alktflx, 0e0, alkdif, alkdec, alkdis, 0e0, alkres ])
        print  'sed:',np.array([pttflx, ptadv,ptdif,  0e0, 0e0, ptrain, ptres])
        print  '   ..... multiple cc species ..... '
        for isp in range(nspcc ):
            print  '{:03}'.format(isp+1), np.array([cctflx[isp], ccadv[isp], ccdif[isp],0e0,ccdis[isp], ccrain[isp], ccres[isp]])
        print '==== burial etc ===='
        print  'z  :',z[0:nz:nz/5]
        print  'w  :',w[0:nz:nz/5]
        print  'rho:',rho[0:nz:nz/5]
        print  'frc:',frt[0:nz:nz/5]
        print ''
        ## in theory, o2dec/ox2om + alkdec = dicdec = omdec (in absolute value)
        if om2cc != 0e0:
            if  np.abs((o2dec/ox2om - alkdec + dicdec)/dicdec) > tol:
                print  abs((o2dec/ox2om + alkdec - dicdec)/dicdec) ,o2dec/ox2om,alkdec,dicdec
                print >> file_err, time, dt \
                    , np.abs((o2dec/ox2om + alkdec - dicdec)/dicdec),o2dec/ox2om,alkdec,dicdec
        print>>file_totfrac, time,np.max(np.abs(frt[:] - 1e0))
        if size :
            if all(w>=0e0):  # not recording when burial is negative 
                print>>file_sigmly,time-age[izrec],d13c_blk[izrec],d18o_blk[izrec] \
                    ,np.sum(ccx[izrec,:]*mcc[:])/rho[izrec]*100e0,ptx[izrec]*msed/rho[izrec]*100e0
                print>>file_sigmlyd, time-age[izrec2],d13c_blk[izrec2],d18o_blk[izrec2] \
                    ,np.sum(ccx[izrec2,:]*mcc[:])/rho[izrec2]*100e0,ptx[izrec2]*msed/rho[izrec2]*100e0
                print>>file_sigbtm,  time-age[nz-1],d13c_blk[nz-1],d18o_blk[nz-1] \
                    ,np.sum(ccx[nz-1,:]*mcc[:])/rho[nz-1]*100e0,ptx[nz-1]*msed/rho[nz-1]*100e0
        else :
            if all(w>=0e0): # not recording when burial is negative 
                print>>file_sigmly, time-age[izrec],d13c_blk[izrec],d18o_blk[izrec] \
                    ,np.sum(ccx[izrec,:]*mcc[:])/rho[izrec]*100e0,ptx[izrec]*msed/rho[izrec]*100e0  \
                    ,d13c_blkf[izrec],d18o_blkf[izrec],np.sum(ccx[izrec,0:4]*mcc[0:4])/rho[izrec]*100e0  \
                    ,d13c_blkc[izrec],d18o_blkc[izrec],np.sum(ccx[izrec,4:8]*mcc[4:8])/rho[izrec]*100e0  
                print>>file_sigmlyd, time-age[izrec2],d13c_blk[izrec2],d18o_blk[izrec2] \
                    ,np.sum(ccx[izrec2,:]*mcc[:])/rho[izrec2]*100e0,ptx[izrec2]*msed/rho[izrec2]*100e0  \
                    ,d13c_blkf[izrec2],d18o_blkf[izrec2],np.sum(ccx[izrec2,0:4]*mcc[0:4])/rho[izrec2]*100e0  \
                    ,d13c_blkc[izrec2],d18o_blkc[izrec2],np.sum(ccx[izrec2,4:8]*mcc[4:8])/rho[izrec2]*100e0  
                print>>file_sigbtm, time-age[nz-1],d13c_blk[nz-1],d18o_blk[nz-1] \
                    ,np.sum(ccx[nz-1,:]*mcc[:])/rho[nz-1]*100e0,ptx[nz-1]*msed/rho[nz-1]*100e0 \
                    ,d13c_blkf[nz-1],d18o_blkf[nz-1],np.sum(ccx[nz-1,0:4]*mcc[0:4])/rho[nz-1]*100e0  \
                    ,d13c_blkc[nz-1],d18o_blkc[nz-1],np.sum(ccx[nz-1,4:8]*mcc[4:8])/rho[nz-1]*100e0  
        # before going to next time step, update variables 
        time = time + dt
        it = it + 1
        o2[:] = o2x[:].copy()
        om[:] = omx[:].copy()
        cc[:] = ccx[:].copy()
        dic[:] = dicx[:].copy()
        alk[:] = alkx[:].copy()
        pt[:] = ptx[:].copy()
        if flg_timereset:
            if not sense and time<time_spn: 
                time = 0.
                it = 1
                print '... because of flg_restart, time is forced back to the beginning ...'
            flg_timereset = False
        
    file_tmp=open(workdir+'sp-trace.txt','w') 
    for isp  in range(nspcc):
        print>>file_tmp,isp+1,d13c_sp[isp],d18o_sp[isp]
    file_tmp.close()
    file_ptflx.close()
    file_ccflx.close()
    file_omflx.close()
    file_o2flx.close()
    file_dicflx.close()
    file_alkflx.close()
    file_err.close()
    file_bound.close()
    file_totfrac.close()
    file_sigmly.close()
    file_sigmlyd.close()
    file_sigbtm.close()
    workdir = 'C:/Users/YK/Desktop/Sed_res/'
    workdir += 'test-translabs/res/'
    workdir += 'python/'
    workdir += 'multi/'
    if not anoxic:
        workdir += 'ox'
    else :
        workdir += 'oxanox'
    if any(labs): workdir += '-labs'
    if any(turbo2): workdir += '-turbo2'
    if any(nobio): workdir += '-nobio'
    if not os.path.exists(workdir):os.makedirs(workdir)
    workdir += '/'
    file_tmp=open(workdir+'lys_sense_'
        +'cc-'+'{:.1e}'.format(ccflxi)+'_rr-'+'{:.1e}'.format(om2cc)  
        +'.txt','a') 
    print>>file_tmp,1e6*(co3i*1e3-co3sat), np.sum(ccx[0,:]*mcc[:])/rho[0]*100e0, frt[0]  \
        ,np.sum(ccx[nz-1,:]*mcc[:])/rho[nz-1]*100e0, frt[nz-1],np.sum(ccx[izml,:]*mcc[:])/rho[izml]*100e0, frt[izml]
    file_tmp.close()
    file_tmp=open(workdir+'ccbur_sense_' 
        +'cc-'+'{:.1e}'.format(ccflxi)+'_rr-'+'{:.1e}'.format(om2cc)  
        +'.txt','a') 
    print>>file_tmp,1e6*(co3i*1e3-co3sat), 1e6*np.sum(ccadv[:])
    file_tmp.close()

def getinput():
    co2chem     = raw_input('Enter how to calculate CO2 chemistry (co2 or mocsy): ')  
    ccflxi      = raw_input('Enter CaCO3 rain flux in mol cm-2 yr-1: ')  
    om2cc       = raw_input('Enter OM/CaCO3 rain ratio: ') 
    dep         = raw_input('Enter water depth in km: ')  
    dt          = raw_input('Enter time step in yr: ')  
    fl          = raw_input('Enter simulation name: ') 
    biot        = raw_input('Enter bioturbation style (nobio, fickian, labs or turbo2): ') 
    oxonly      = raw_input('Oxic only for OM degradation? (True or False): ') 
    runmode     = raw_input('Enter simulation mode (sense, diss. exp., size, biotest or track2): ') 
    if len(co2chem)==0:                # no input
        co2chem = 'co2'                # default
    if len(ccflxi)==0:                # no input
        ccflxi = 12e-6                # default
    else:
        ccflxi = eval(ccflxi)
    if len(om2cc)==0:                 # no input
        om2cc = 0.7                   # default
    else:
        om2cc = eval(om2cc)
    if len(dep)==0:                   # no input
        dep = 3.5                     # default 
    else:
        dep = eval(dep)
    if len(dt)==0:                    # no input
        dt = 1000.                    # default
    else:
        dt = eval(dt)
    if len(fl)==0:                    # no input
        fl = '50kyr-dis4_5'           # default or you can specify the name here
    if len(oxonly)==0:                # no input
        oxonly = False                # default 
    else:
        oxonly = eval(oxonly)
    # co2chem = 'co2'
    # co2chem = 'co2h2o'              # dic = co2+hco3+co3; alk = hco3+2*co3+oh-h
    sparse = True
    showiter = False
    return ccflxi,om2cc,dep,dt,fl,biot,oxonly,runmode,co2chem,sparse,showiter

if __name__ == '__main__':
    ccflxi,om2cc,dep,dt,fl,biot,oxonly,runmode,co2chem,sparse,showiter = getinput()
    caco3_main(ccflxi,om2cc,dep,dt,fl,biot,oxonly,runmode,co2chem,sparse,showiter)
