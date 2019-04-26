program chk_bur
! just check the calculation of coefficients for sediment burial using subroutines in caco3_test_mod_v5_6.f90

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran -cpp -I/path/to/working/directory chk_bur.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! <<<< Maybe first only consider Fickian mixing so switch off all macros in defines.h (you can switch on 'test') >>>>>

! you can check coefficients for sediment burial 
! please copy and paste results to whatever file to be compared with results with MATLAB version 

use globalvariables
implicit none 
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
! stop
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

print*,'z   ,w  ,age    ,up    ,dwn ,cnr    ,adf    '
! showing parameters relevant to burial on screen 
do iz=1,nz 
    print*,z(iz),w(iz),age(iz),up(iz),dwn(iz),cnr(iz) ,adf(iz)
    ! depth [cm], burial rate [cm yr-1], grid age [yr], and factors for burial advection; 
    ! up, dwn, cnr uses (i-1,i),(i,i+1) and (i-1,i+1), respectively for burial at i. 
    ! adf is the factor to account for mass balance when cnr is not zero.
enddo

endprogram 