program chk_om
! just check the calculation of om conc. and flx using subroutines in caco3_test_mod_v5_6.f90

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran -cpp -I/path/to/working/directory chk_om.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! you can check om calculation including om conc. and flx 
! please copy and paste results to whatever file to be compared with results with MATLAB version 
! NOTE: here oxygen conc. has to be assume. The model calculate om and o2 iteratively, but not here. 

use globalvariables
implicit none 
integer(kind=4) interval  ! choose value between 1 to nz 

interval =10 ! choose a value between 1 to nz; om depth profile is shown with this interval; e.g., if inteval = nz, om conc. at all depths are shown
! e.g., if interval = 5, om conc. at 5 depths are shown   

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

beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z)

! call getporosity() ! assume porosity profile 
call getporosity(  &
     poro,porof,sporo,sporof,sporoi & ! output
     ,z,nz  & ! input
     )

!!!!!!!!!!!!! flx assignement and initial guess for burial rate !!!!!!!!!!!!!!!!!!!!!!
! call flxstat()  ! assume fluxes of om, cc and clay, required to calculate burial velocity
call flxstat(  &
    omflx,detflx,ccflx  & ! output
    ,om2cc,ccflxi,mcc,nspcc  & ! input 
    )
! print*,om2cc,ccflxi,detflx,omflx,sum(ccflx)
! molar volume (cm3 mol-1) needed for burial rate calculation 
mvom = mom/rhoom  ! om
mvsed = msed/rhosed ! clay 
mvcc = mcc/rhocc ! caco3
! call burial_pre() ! initial guess of burial profile, requiring porosity profile  
call burial_pre(  &
    w,wi  & ! output
    ,detflx,ccflx,nspcc,nz  & ! input 
    )
! call dep2age() ! depth -age conversion 
call dep2age(  &
    age &  ! output 
    ,dz,w,nz  &  ! input
   )
! call calcupwindscheme() ! determine factors for upwind scheme to represent burial advection
call calcupwindscheme(  &
    up,dwn,cnr,adf & ! output 
    ,w,nz   & ! input &
    )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call make_transmx()
call make_transmx(  &
    trans,izrec,izrec2,izml,nonlocal  & ! output 
    ,labs,nspcc,turbo2,nobio,dz,sporo,nz,z  & ! input
    )

om = 1d-8  ! assume an arbitrary low conc. 
o2 = o2i*1d-6/1d3 ! o2 conc. in uM converted to mol/cm3

omx = om
o2x = o2

time = 0d0 ! model time [yr]
it = 1 ! integration count 
nt = 10 ! total integration 
dt = 100d0 ! time step [yr]

rho = 2.5d0 ! assume here 

do it=1,nt
    print'(A,i0,A,E11.3,A,E11.3,A)','(it,dt,time)  (',it,',',dt,',',time,')'
    
    ! call omcalc() ! om conc. calculation 
    call omcalc( &
        omx,izox  & ! output 
        ,kom   &  ! in&output
        ,oxic,anoxic,o2x,om,nz,sporo,sporoi,sporof &! input 
        ,w,wi,dt,up,dwn,cnr,adf,trans,nspcc,labs,turbo2,nonlocal,omflx,poro,dz &! input 
        ) 
    ! calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
    ! call calcflxom()
    call calcflxom(  &
        omadv,omdec,omdif,omrain,omflx,omres,omtflx  & ! output 
        ,sporo,om,omx,dt,w,dz,z,nz,turbo2,labs,nonlocal,poro,up,dwn,cnr,adf,rho,mom,trans,kom,sporof,sporoi,wi,nspcc  & ! input 
        )
    
    write(dumchr(2),'(i0)') interval
    dumchr(1)="(A,"//trim(adjustl(dumchr(2)))//"E11.3"//")"
    ! showing results on screen
    print*,'~~~~ conc ~~~~'
    print dumchr(1), 'z  :',(z(iz),iz=1,nz,nz/interval)
    print dumchr(1), 'om :',(omx(iz)*mom/rho(iz)*100d0,iz=1,nz,nz/interval)
    print*,'++++ flx ++++'
    print'(7A11)', 'tflx','adv','dif','omrxn','ccrxn','rain','res'
    print'(A,7E11.3)', 'om :', omtflx, omadv,  omdif, omdec,0d0,omrain, omres
    
    print*,''
    print*,''
    print*,''
    
    om = omx
    time = time +dt
    
enddo

endprogram 