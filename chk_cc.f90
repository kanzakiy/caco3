program chk_cc
! just check the calculation of conc. and flx of caco3, dic and alk using subroutines in caco3_test_mod_v5_6.f90

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran -cpp -I/path/to/working/directory chk_cc.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! you can check conc. and flx of caco3, alk, dic 
! please copy and paste results to whatever file to be compared with results with MATLAB version 
! NOTE: here decomposition of om is ignored. The test including decomposition of om is made in another file 
! because burial is not adjusted, cc conc. can goes > 100 wt% in this test 

use globalvariables
implicit none 
integer(kind=4) interval  ! choose value between 1 to nz 

interval =10 ! choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
! e.g., if interval = 5, om conc. at 5 depths are shown   

workdir = './' ! working directory 

#include <defines.h>

#ifdef allnobio 
nobio = .true.
#elif defined allturbo2 
turbo2 = .true.
#elif defined alllabs 
labs = .true.
#endif 

#ifdef oxonly
anoxic = .false. 
#endif


dep = 4.0d0 ! km water depth; note that temperature and salinitiy has initially assumed values in globalvariables.mod  

beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z)

call getporosity() ! assume porosity profile 

!!!!!!!!!!!!! flx assignement and initial guess for burial rate !!!!!!!!!!!!!!!!!!!!!!
call flxstat()  ! assume fluxes of om, cc and clay, required to calculate burial velocity
! print*,om2cc,ccflxi,detflx,omflx,sum(ccflx)
! molar volume (cm3 mol-1) needed for burial rate calculation 
mvom = mom/rhoom  ! om
mvsed = msed/rhosed ! clay 
mvcc = mcc/rhocc ! caco3
call burial_pre() ! initial guess of burial profile, requiring porosity profile  
call dep2age() ! depth -age conversion 
call calcupwindscheme() ! determine factors for upwind scheme to represent burial advection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call make_transmx()

call coefs(temp,sal,dep)  ! need to specify diffusion coefficient as well as om decomposition rate const. etc.

!!   INITIAL CONDITIONS !!!!!!!!!!!!!!!!!!! 
cc = 1d-8   ! assume an arbitrary low conc. 
dic = dici*1d-6/1d3 ! mol/cm3; factor is added to change uM to mol/cm3 
alk = alki*1d-6/1d3 ! mol/cm3

call calcspecies(dic,alk,temp,sal,dep,pro,co2,hco3,co3,nz,infosbr)   ! calculation of individual co2 species at initial conditions

! passing to transient variables 
ccx = cc
dicx = dic
alkx = alk 

! this may not be necessary as these individual species assume equilibrium 
co2x = co2
hco3x = hco3
co3x = co3

time = 0d0 ! model time [yr]
it = 1 ! integration count 
nt = 10 ! total integration 
dt = 1d6 ! time step [yr]

rho = 2.5d0 ! assume here density (this is going to be calculated based on solid phase composition )

oxco2 = 0d0  ! oxic degradation of om; here assumed 0 at all time and depth  
anco2 = 0d0  ! anoxic degradation of om; here assumed 0 at all time and depth  

do it=1,nt
    print'(A,i0,A,E11.3,A,E11.3,A)','(it,dt,time)  (',it,',',dt,',',time,')'
    
    call calccaco3sys()
    if (flg_500) then
        print*,'error after calccaco3sys'
        stop
    endif 
    
    call calcspecies(dicx,alkx,temp,sal,dep,prox,co2x,hco3x,co3x,nz,infosbr)
    if (infosbr==1) then 
        print*,'error after calcspecies after calccaco3sys'
        stop
    endif 
    
    ! calculation of fluxes relevant to caco3 and co2 system
    call calcflxcaco3sys()
    
    write(dumchr(2),'(i0)') interval
    dumchr(1)="(A,"//trim(adjustl(dumchr(2)))//"E11.3"//")"
    ! showing results on screen
    print*,'~~~~ conc ~~~~'
    print dumchr(1), 'z  :',(z(iz),iz=1,nz,nz/interval)
    print dumchr(1), 'cc :',(sum(ccx(iz,:)*mcc)/rho(iz)*100d0,iz=1,nz,nz/interval)
    print dumchr(1), 'dic:',(dicx(iz)*1d3,iz=1,nz,nz/interval)
    print dumchr(1), 'alk:',(alkx(iz)*1d3,iz=1,nz,nz/interval)
    print*, '   ..... multiple cc species ..... '
    write(dumchr(2),'(i0)') interval
    dumchr(1)="(i0.3,':',"//trim(adjustl(dumchr(2)))//"E11.3"//")"
    do isp=1,nspcc 
        print dumchr(1),isp,(ccx(iz,isp)*mcc/rho(iz)*100d0,iz=1,nz,nz/interval)
    enddo
    print*,'++++ flx ++++'
    print'(7A11)', 'tflx','adv','dif','omrxn','ccrxn','rain','res'
    print'(A,7E11.3)', 'cc :',sum(cctflx),  sum(ccadv), sum(ccdif),0d0,sum(ccdis), sum(ccrain), sum(ccres) 
    print'(A,7E11.3)', 'dic:',dictflx, 0d0,dicdif, dicdec,  dicdis, 0d0,dicres 
    print'(A,7E11.3)', 'alk:',alktflx, 0d0, alkdif, alkdec, alkdis, 0d0, alkres 
    print*, '   ..... multiple cc species ..... '
    do isp=1,nspcc 
        print'(i0.3,":",7E11.3)',isp,cctflx(isp), ccadv(isp), ccdif(isp),0d0,ccdis(isp), ccrain(isp), ccres(isp) 
    enddo
    
    print*,''
    print*,''
    print*,''
    
    cc = ccx
    dic = dicx
    alk = alkx
    time = time +dt
    
enddo

endprogram 