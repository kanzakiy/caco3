program chk_om
! just check the calculation of om conc. and flx using subroutines in caco3_test_mod_v5_6.f90

! (1) gfortran -c caco3_therm.f90
! (2) gfortran -c -cpp -I/path/to/working/directory caco3_test_mod_v5_6.f90 
! (3) gfortran -cpp -I/path/to/working/directory chk_om.f90 caco3_test_mod_v5_6.o caco3_therm.o -lopenblas -g -fcheck=all
! (4) ./a.out

! you can om calculation including om conc. and flx 
! please copy and paste results to whatever file to be compared with results with MATLAB version 
! NOTE: here oxygen conc. has to be assume. The model calculate om and o2 iteratively, but not here. 

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

call getporosity() ! assume porosity profile 

call make_transmx()

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
    
    call omcalc() ! om conc. calculation 
    ! calculating the fluxes relevant to om diagenesis (and checking the calculation satisfies the difference equations )
    call calcflxom()
    
    ! showing results on screen
    print*,'~~~~ conc ~~~~'
    print'(A,5E11.3)', 'z  :',(z(iz),iz=1,nz,nz/5)
    print'(A,5E11.3)', 'om :',(omx(iz)*mom/rho(iz)*100d0,iz=1,nz,nz/5)
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